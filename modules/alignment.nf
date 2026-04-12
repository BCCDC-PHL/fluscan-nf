// alignment.nf - Read alignment and consensus generation processes

process map_reads {

    tag {sample_id}

    errorStrategy 'ignore'

	label 'heavy'

    input: 
    tuple val(sample_id), path(reads_1), path(reads_2), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit:provenance

    script:
	"""
    printf -- "- process_name: map_reads\\n" > ${sample_id}-bwa-provenance.yml
    printf -- "  tools: \\n  - tool_name: bwa\\n    tool_version: \$(bwa 2>&1 | head -n3 | tail -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-bwa-provenance.yml

    bwa index ${reference}
	bwa mem -t ${task.cpus} ${reference} ${reads_1} ${reads_2} > ${sample_id}.sam
	"""
}


process map_reads_long {

    tag { "${name}" }

    errorStrategy 'ignore'
    label 'heavy'

    input:
    tuple val(sample_id), val(segment), path(consensus_fasta), path(reference_fasta)

    output:
    tuple val(sample_id), val(segment), path("${name}.paf"), path(reference_fasta), emit: paf
    tuple val(sample_id), path("${provenance_file}"), emit: provenance

    script:
    provenance_file = "${sample_id}_map_reads_long_provenance.yml"
    name = "${sample_id}_${segment}"
    """
    printf -- "- process_name: map_reads_long\\n" > ${provenance_file}
    printf -- "  tools: \\n" >> ${provenance_file}
    printf -- "  - tool_name: minimap2\\n    tool_version: \$(minimap2 --version 2>&1 | head -n1)\\n" >> ${provenance_file}

    minimap2 -cx asm20 --cs ${reference_fasta} ${consensus_fasta} > ${name}.paf
    """
}


process sort_filter_sam {

	tag {sample_id}

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*{depth.tsv}", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*sort.filter.bam", mode: 'symlink'


    input: 
    tuple val(sample_id), path(samfile)

    output:
    tuple val(sample_id), path("${sample_id}.sort.filter.bam"), path("${sample_id}*bai"), emit: bam
    tuple val(sample_id), path("${sample_id}_depth.tsv"), emit: depth
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
	"""
    printf -- "- process_name: sort_filter_sam\\n" > ${sample_id}-samtools-provenance.yml
    printf -- "  tools: \\n  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-samtools-provenance.yml

	samtools view -f 3 -F 2828 -q 30 -h ${samfile} | samtools sort -o ${sample_id}.sort.filter.bam 
    samtools index ${sample_id}.sort.filter.bam 

    samtools depth -a ${sample_id}.sort.filter.bam  > ${sample_id}_depth.tsv
	"""
}


process make_pileup {

	tag {sample_id}

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*tsv", mode: 'copy'

    input: 
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    tuple val(sample_id), path("${sample_id}*pileup.tsv"), emit: pileup

    script:
	"""
    make_pileup.py --input ${bam} --output ${sample_id}_pileup.tsv
	"""
}


process mask_low_coverage {
	tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*.bed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bamfile), path(bam_index)
	
    output:
    tuple val(sample_id), path("${sample_id}_low_coverage.bed"), emit: masked_bed
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: mask_low_coverage\\n" > ${sample_id}-bedtools-provenance.yml
    printf -- "  tools: \\n  - tool_name: bedtools\\n    tool_version: \$(bedtools --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}-bedtools-provenance.yml

    bedtools genomecov -bga -ibam ${bamfile} |  awk '\$4 < ${params.min_depth} {{print}}' | awk 'BEGIN{FS=OFS="\\t"} {print \$1,\$2+1,\$3+1,\$4}' > ${sample_id}_low_coverage.bed
    """
}


process make_consensus {

    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*consensus.fasta", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_index), path(reference), path(mask_file)
	
    output:
    tuple val(sample_id), path("${sample_id}*consensus.fasta"), emit: consensus
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: make_consensus\\n" > ${sample_id}-bcftools-provenance.yml
    printf -- "  tools: \\n  - tool_name: bcftools\\n    tool_version: \$(bcftools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-bcftools-provenance.yml

    bcftools consensus -m ${mask_file} -f ${reference} ${vcf} |  sed 's/^>/>${sample_id}|/' > ${sample_id}_consensus.fasta

    """
}
