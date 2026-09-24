// shared/sample_id.nf - the one place a consensus or sample sheet sample_id is derived.
//
// The value this produces is final. It keys every channel and names every file
// from prepare_consensus onward, so nothing downstream may re-derive or edit it:
// processes take it as val(sample_id) rather than parsing it back out of a name.
//
// The rule, applied to a sample sheet ID or to a filename with all extensions dropped:
//   1. strip a trailing _consensus (the .consensus form is already gone with the
//      extensions, e.g. SAMPLE-A.consensus.fasta -> SAMPLE-A)
//   2. underscores to dashes, since underscores separate the segment label in
//      the canonical {sample_id}_{SEGMENT} header. Two consumers read the ID back
//      as everything before the FIRST underscore and would truncate a name that
//      kept one: the saveAs in modules/nextclade/nextclade.nf and the sample_name
//      column in bin/filter_alignments.py
//   3. strip the submitter prefix some partner organisations prepend,
//      e.g. 24-1234-IP7-SAMPLE-A -> SAMPLE-A

def clean_sample_id(String name) {
    name.replaceFirst(/(?i)_consensus$/, '')
        .replaceAll(/_/, '-')
        .replaceFirst(/^\d{2}-\d+(?:-IP\d+)?-/, '')
}

// Consensus FASTA discovery for the PAFTOOLS and NEXTCLADE branches. Each file
// holds the segments of one sample.
def consensus_fasta_channel(String input_dir) {
    def ch = Channel.fromPath(input_dir + "/*.{fasta,fa,fna}", checkIfExists: true)
                    .map { fasta -> [clean_sample_id(fasta.simpleName), fasta] }

    // Two filenames collapsing to one sample_id would share an output directory and
    // cross-join downstream, silently. Discovery completes at startup, so this aborts
    // the run before any task can finish.
    ch.map { id, fasta -> id }.toList().subscribe { ids ->
        def duplicates = ids.countBy{ it }.findAll{ it.value > 1 }.keySet()
        if (duplicates) {
            error "ERROR: input filenames collapse to duplicate sample IDs: ${duplicates.join(', ')}"
        }
    }

    return ch
}

// Sample sheet discovery for every workflow: an ID column plus the given file columns,
// ['R1', 'R2'] or ['FASTA']. Parsed eagerly, so a bad sheet fails before any task runs.
// Relative paths resolve against the launch directory.
def samplesheet_channel(String path, List<String> columns) {
    def rows = file(path, checkIfExists: true).splitCsv(header: true, strip: true)
    if (!rows) {
        error "ERROR: sample sheet ${path} has no rows."
    }
    def missing = (['ID'] + columns) - rows[0].keySet()
    if (missing) {
        error "ERROR: sample sheet ${path} is missing column(s): ${missing.join(', ')}"
    }
    // file('') throws without naming the row, so blank cells are caught here first.
    rows.eachWithIndex { row, i ->
        def blank = (['ID'] + columns).findAll { !row[it] }
        if (blank) {
            error "ERROR: sample sheet ${path} line ${i + 2} has no value for: ${blank.join(', ')}"
        }
    }

    // Cleaning cannot fix '|' (the SHORT_READ segment_regex delimiter), '/' or whitespace.
    def invalid = rows*.ID.findAll { !(clean_sample_id(it) ==~ /[A-Za-z0-9.-]+/) }
    if (invalid) {
        error "ERROR: sample IDs in ${path} may contain only letters, digits, '.', '-' and '_': ${invalid.join(', ')}"
    }
    def ids = rows.collect { clean_sample_id(it.ID) }
    def duplicates = ids.countBy{ it }.findAll{ it.value > 1 }.keySet()
    if (duplicates) {
        error "ERROR: IDs in ${path} collapse to duplicate sample IDs: ${duplicates.join(', ')}"
    }

    // [sample_id, file, ...] in the order of columns, e.g. [sample_id, R1, R2]
    return Channel.fromList([ids, rows].transpose().collect { id, row ->
        [id] + columns.collect { file(row[it], checkIfExists: true) }
    })
}
