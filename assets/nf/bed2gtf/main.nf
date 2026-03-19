process BED2GTF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' : 
        'ghcr.io/alejandrogzi/bed2gtf:latest' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta1), path(isoforms)

    output:
    tuple val(meta), path("*.gtf"),    optional: true, emit: gtf
    tuple val(meta), path("*.gtf.gz"), optional: true, emit: gtf_gz
    tuple val(meta), path("*.gff"),    optional: true, emit: gff
    tuple val(meta), path("*.gff.gz"), optional: true, emit: gff_gz
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def format = task.ext.format ?: 'gtf'
    def prefix = task.ext.prefix ?: bed.baseName + '.' + format
    def isoforms = isoforms ? "--isoforms ${isoforms}" : ''
    """
    bed2gtf \\
        $args \\
        $isoforms \\
        -T ${task.cpus} \\
        -i ${bed} \\
        -o $prefix
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bed2gtf: \$(bed2gtf --version | sed -e "s/bed2gtf v//g")
    END_VERSIONS
    """

    stub:
    """
    touch *.gtf
    touch *.gtf.gz
    touch *.gff
    touch *.gff.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bed2gtf: \$(bed2gtf --version | sed -e "s/bed2gtf v//g")
    END_VERSIONS
    """
}
