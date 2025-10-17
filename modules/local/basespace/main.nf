process BASESPACE {
    tag "$run_id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/csgenetics/basespace-cli:0.1' :
        'quay.io/csgenetics/basespace-cli@sha256:514422008a5e32ba3795c2d61992dd3c9aaf589156440cd378e5a3157d82eb7d' }"

    input:
    val run_id
    path config

    output:
    path "*/*fastq.gz", emit: fastqs
    path "*/*json"    , emit: jsons

    script:
    """
    cp /bin/bs .
    RUNID=\$(./bs list run | grep $run_id | awk '{print \$4}')

    #this creates the list of all the lanes for downloading
    ./bs list dataset --input-run \$RUNID | awk '{print \$2;}' | grep 'OG' > ${run_id}.prefix.txt
    sed -i '1,3d' ${run_id}.prefix.txt


    for PREFIX in \$(cat ${run_id}.prefix.txt); do
        ID=\$(./bs list dataset --input-run \$RUNID | grep \$PREFIX | awk '{print \$4;}')
        echo \$PREFIX \$ID ">>" $run_id
        ./bs download dataset ---input-run \$RUNID -i \$ID -o \$PREFIX
    done
    
    """
}