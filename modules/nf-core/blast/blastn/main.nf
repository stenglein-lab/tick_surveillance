process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    // if using conda                                                             
    conda "${moduleDir}/environment.yml"                                          
                                                                                  
    // if using singularity                                                       
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/blast:2.16.0--hc155240_2"
    } else {                                                                      
        container "quay.io/biocontainers/blast:2.16.0--hc155240_2"                
    }   

    input:
    tuple val(meta), path(fasta)
    path  db

    output:
    tuple val(meta), path('*.blastn.txt'), emit: txt
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    DB=`find -L ./ -name "*.ndb" | sed 's/\\.ndb\$//'`
    blastn \\
        -num_threads $task.cpus \\
        -db \$DB \\
        -query $fasta \\
        $args \\
        -out ${prefix}.blastn.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
