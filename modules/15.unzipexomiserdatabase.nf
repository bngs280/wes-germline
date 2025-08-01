process unzip_database {
    input:
    val genome
    
    output:
    path "unzipped_db", emit: unzipped_db
    
    script:
    """
    mkdir -p /usr/src/app/exomiser/db/unzipped_db
    7z x -aoa /usr/src/app/exomiser/db/2406_phenotype.zip -o/usr/src/app/exomiser/db/unzipped_db
    
    if [[ "$genome" == "hg19" ]]; then
        echo "Unzipping hg19 database..."
        unzip -o /usr/src/app/exomiser/db/2406_hg19.zip -d /usr/src/app/exomiser/db/unzipped_db
    elif [[ "$genome" == "hg38" ]]; then
        echo "Unzipping hg38 database..."
        unzip -o /usr/src/app/exomiser/db/2406_hg38.zip -d /usr/src/app/exomiser/db/unzipped_db
    else
        echo "Invalid genome specified: $genome"
        exit 1
    fi
    ln -s /usr/src/app/exomiser/db/unzipped_db unzipped_db
    """
}
