process DOWNLOAD_TAXONKIT_DB {
    tag "${db_name}"
    label 'process_low'
    storeDir params.taxonkit_db_dir // Cache the database

    input:
    val db_name

    output:
    path("taxonkit_dbs"), emit: db_files

    script:
    db_dir = params.taxonkit_db_dir
    """
    set -euo pipefail

    mkdir -p taxonkit_dbs

    # Reuse the cached dump only when the two files we actually read are present
    # AND non-empty. A -f test is not enough: a half-finished extract leaves
    # zero-byte .dmp files behind, and an empty names.dmp resolves nothing while
    # looking exactly like a cache hit.
    if [ -s "${db_dir}/taxonkit_dbs/names.dmp" ] \\
    && [ -s "${db_dir}/taxonkit_dbs/nodes.dmp" ]; then
        echo "Files already exist — linking to work dir"
        for f in "${db_dir}"/taxonkit_dbs/*; do
            ln -s "\$f" taxonkit_dbs/
        done
    else
        echo "Downloading fresh taxonomy database..."
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -xzf taxdump.tar.gz -C taxonkit_dbs
        rm taxdump.tar.gz
    fi

    # Never publish a cache that would silently resolve nothing.
    for f in names.dmp nodes.dmp; do
        if [ ! -s "taxonkit_dbs/\$f" ]; then
            echo "ERROR: taxonkit_dbs/\$f is missing or empty after setup" >&2
            exit 1
        fi
    done
    """

    stub:
    """
    mkdir -p taxonkit_dbs
    touch taxonkit_dbs/citations.dmp \\
        taxonkit_dbs/delnodes.dmp \\
        taxonkit_dbs/division.dmp \\
        taxonkit_dbs/gencode.dmp \\
        taxonkit_dbs/images.dmp \\
        taxonkit_dbs/merged.dmp \\
        taxonkit_dbs/names.dmp \\
        taxonkit_dbs/nodes.dmp \\
        taxonkit_dbs/gc.prt
    """
}
