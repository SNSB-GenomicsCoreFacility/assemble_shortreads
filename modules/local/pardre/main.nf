process PARDRE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pardre:2.2.5--h2aad775_4':
        'biocontainers/pardre:2.2.5--h2aad775_4' }"
    publishDir("${params.outdir}/${index_out}", mode:'copy')

    input:// TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    tuple val(meta), path(reads)
    val (index_out)

    output:
    tuple val(meta), path("*.gz"), emit: reads
    tuple val(meta), path("*_log.txt"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_dedup_padre"
    def input = " -i ${reads[0]} -p ${reads[1]} "
    def output = " -o ${prefix}_1.fq.gz -r ${prefix}_2.fq.gz "
    """
    ParDRe \\
        $input \\
        $output \\
        -t $task.cpus \\
        $args \\
        &>${prefix}_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pardre: 2.1.5-PC170509
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pardre: 2.1.5-PC170509
    END_VERSIONS
    """
}
