// shared/database.nf - Reference database bundle resolution and validation.
//
// A "bundle" is a directory holding every reference the pipeline needs for one
// reference set, in a fixed layout:
//
//   <db>/
//   |-- db.yml                 optional metadata; every key optional
//   |-- blastn_db.fa           nt reference, headers accession|subtype|segment
//   |-- nextclade/             one directory per segment, uppercase (HA/, NA/, ...)
//   `-- snpeff/
//       |-- snpEff.config      its data.dir is overridden at runtime by -dataDir
//       `-- data/<genome>/     snpEffectPredictor.bin etc.
//
// The filenames are fixed on purpose. A bundle is identified by its layout, so this
// file is the only place that knows where anything inside one lives; --db is the only
// path the caller supplies. --flumut_db points at a second, independent bundle whose
// nextclade/ datasets are the source coordinate system for the optional flumut merge.
//
// Bundles may be partial. Each workflow declares the components it consumes and only
// those are validated: SHORT_READ requires 'blastn' (and picks up snpeff annotation
// when the bundle happens to carry it), PAFTOOLS requires 'blastn' and 'snpeff',
// NEXTCLADE requires 'nextclade'.

import org.yaml.snakeyaml.Yaml

// The canonical segment set. Nextclade datasets are resolved using these labels directly,
// and coordinate_conversion shards its map by exactly these labels.
@groovy.transform.Field
List DB_SEGMENTS = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']


def db_layout_help() {
    return """
Expected bundle layout:

  <db>/
    db.yml                                       (optional metadata)
    blastn_db.fa                                 (nt reference; headers accession|subtype|segment)
    nextclade/{PB2,PB1,PA,HA,NP,NA,M,NS}/        (reference.fasta + genome_annotation.gff3 each)
    snpeff/snpEff.config
    snpeff/data/<genome>/snpEffectPredictor.bin

Not every component is needed by every workflow; see docs/database_bundles.md.
""".stripIndent()
}


// snakeyaml turns an unquoted `version: 2025-07-07` into a java.util.Date. Render it back
// in UTC -- formatting in local time would shift the day and silently misreport the DB
// version in the provenance record.
def db_version_string(value) {
    if (value == null) {
        return 'unknown'
    }
    if (value instanceof Date) {
        def formatter = new java.text.SimpleDateFormat('yyyy-MM-dd')
        formatter.setTimeZone(java.util.TimeZone.getTimeZone('UTC'))
        return formatter.format(value)
    }
    return value.toString()
}


def read_db_metadata(root) {
    def manifest = root.resolve('db.yml')
    if (!manifest.exists()) {
        return [:]
    }
    def parsed = null
    try {
        parsed = new Yaml().load(manifest.text)
    }
    catch (Exception e) {
        error "ERROR: could not parse ${manifest}: ${e.message}"
    }
    if (parsed == null) {
        return [:]
    }
    if (!(parsed instanceof Map)) {
        error "ERROR: ${manifest} must contain a YAML mapping, found ${parsed.getClass().simpleName}."
    }
    return parsed
}


// Segments whose Nextclade dataset directory is missing or incomplete. Checked up front
// because the alternative is a confusing failure inside nextclade_run_batch, several
// minutes into a run, for a dataset that was never going to work.
def incomplete_nextclade_segments(nextclade_dir) {
    return DB_SEGMENTS.findAll { segment ->
        def dataset = nextclade_dir.resolve(segment)
        !dataset.isDirectory() ||
            !dataset.resolve('reference.fasta').exists() ||
            !dataset.resolve('genome_annotation.gff3').exists()
    }
}


// Resolve a bundle root into the paths and metadata the workflows consume, validating
// the components named in `required` (a subset of ['blastn', 'nextclade', 'snpeff']).
def resolve_db(String path, List required, String param_name = '--db') {

    if (!path || path == 'NO_FILE') {
        error "ERROR: ${param_name} must be set to a reference database bundle.\n" + db_layout_help()
    }

    def root = file(path)
    if (!root.exists()) {
        error "ERROR: ${param_name} ${path} does not exist."
    }
    if (!root.isDirectory()) {
        error "ERROR: ${param_name} ${path} is not a directory; it must be a bundle root.\n" + db_layout_help()
    }

    def blastn        = root.resolve('blastn_db.fa')
    def nextclade     = root.resolve('nextclade')
    def snpeff_dir    = root.resolve('snpeff')
    def snpeff_config = snpeff_dir.resolve('snpEff.config')
    def snpeff_data   = snpeff_dir.resolve('data')

    def present = [
        blastn:    blastn.exists(),
        nextclade: nextclade.isDirectory(),
        snpeff:    snpeff_config.exists() && snpeff_data.isDirectory(),
    ]

    def problems = []
    required.each { component ->
        if (present[component]) {
            return
        }
        switch (component) {
            case 'blastn':
                problems << "  missing nucleotide BLAST reference: ${blastn}"
                break
            case 'nextclade':
                // The single most likely mistake is pointing at the dataset directory
                // itself rather than at the bundle that contains it.
                if (root.resolve('HA').isDirectory()) {
                    problems << "  missing ${nextclade} -- ${path} looks like a Nextclade dataset\n" +
                                "  directory rather than a bundle root (it contains HA/). Pass its parent instead:\n" +
                                "    ${param_name} ${root.parent}"
                }
                else {
                    problems << "  missing Nextclade dataset directory: ${nextclade}"
                }
                break
            case 'snpeff':
                if (!snpeff_config.exists()) {
                    problems << "  missing SnpEff config: ${snpeff_config}"
                }
                if (!snpeff_data.isDirectory()) {
                    problems << "  missing SnpEff data directory: ${snpeff_data}"
                }
                break
            default:
                problems << "  unknown database component '${component}'"
        }
    }

    // Reported alongside any structural problems above rather than only after they are
    // fixed, so one failed run surfaces everything wrong with the bundle.
    if (present.nextclade && required.contains('nextclade')) {
        def incomplete = incomplete_nextclade_segments(nextclade)
        if (incomplete) {
            problems << "  Nextclade datasets missing or incomplete for: ${incomplete.join(', ')}\n" +
                        "  (each needs ${nextclade}/<segment>/{reference.fasta,genome_annotation.gff3})"
        }
    }

    if (problems) {
        error "ERROR: incomplete reference database bundle at ${root}:\n" +
              problems.join('\n') + '\n' + db_layout_help()
    }

    def meta = read_db_metadata(root)
    def name = meta.name ?: root.name

    return [
        root:               root,
        name:               name,
        version:            db_version_string(meta.version),
        // Named separately from `name` so a bundle can be renamed on disk without
        // changing the reference labels already written into published coordinate maps.
        nextclade_ref_name: meta.nextclade_ref_name ?: name,
        blastn:             blastn,
        nextclade:          nextclade,
        snpeff_dir:         snpeff_dir,
        snpeff_config:      snpeff_config,
        snpeff_data:        snpeff_data,
        has_blastn:         present.blastn,
        has_nextclade:      present.nextclade,
        has_snpeff:         present.snpeff,
        meta:               meta,
    ]
}


// The flumut bundle is resolved independently of --db: it changes far less often, and
// only the optional merge step reads it. Only its nextclade/ datasets are consumed.
def resolve_flumut_db(String path) {
    if (!path || path == 'NO_FILE') {
        error "ERROR: --flumut_db must be set to the FluMut reference bundle when " +
              "--merge_flumut_mutations is enabled.\n" + db_layout_help()
    }
    return resolve_db(path, ['nextclade'], '--flumut_db')
}


def db_summary(Map db) {
    def components = [
        db.has_blastn    ? 'blastn'    : null,
        db.has_nextclade ? 'nextclade' : null,
        db.has_snpeff    ? 'snpeff'    : null,
    ].findAll().join(', ')

    return "${db.name} (version: ${db.version}; components: ${components ?: 'none'})\n" +
           "                      ${db.root}"
}
