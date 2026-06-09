#!/bin/bash
echo "Start database building script"
makeblastdb
    -in unite2024ITS.fasta
    -dbtype nucl
    -out unite_its_database

blastn
    -query query.fasta
    -db ~/databases/unite/unite_its_database
    -out results.xml #TODO
    -outfmt 5
    -evalue e_value_threshold #TODO
echo "End database building script"
