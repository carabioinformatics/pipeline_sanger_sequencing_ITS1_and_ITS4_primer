#!/bin/bash

mkdir -p Output
#window_size = $1
#quality_threshold = $2
#taxonomy_mode = $3
#e_value_threshold = $4

for forward in Input/*ITS1___1.ab1
do
    # Extract sample name
    sample_name=$(basename "$forward" | cut -d'_' -f2)

    # Find reverse filename
    reverse=$(find . -name "*_${sample_name}_ITS4___1.ab1")

    # Get consensus file name
    consensus_name=$"INT26_${sample_name}_CONSENSUS.fas"

    # Start printing
    echo "================================================="
    echo "Processing sample: $sample_name"
    echo "Forward input file: $forward"
    echo "Reverse input file: $reverse"

    # Check if reverse file exists
    if [[ -z "$reverse" ]]; then
        echo "Error: Reverse file not found for $sample_name."
        continue
    fi

    # Make relevant directories
    mkdir -p ./Output/$sample_name/BLAST
    mkdir -p ./Output/$sample_name/Cleanup
    mkdir -p ./Output/$sample_name/Final
    blast_attachment="./Output/$sample_name/BLAST"
    cleanup_attachment="./Output/$sample_name/Cleanup"
    final_attachment="./Output/$sample_name/Final"

    # Run scrips now
    python3 ./src/TRIMMING_AND_CONSENSUS.py $forward $reverse $1 $2 $cleanup_attachment $final_attachment $blast_attachment
    python3 ./src/BLAST_AND_TAXONOMY.py ./Output/$sample_name/BLAST/$consensus_name $3 $4 $cleanup_attachment $final_attachment $blast_attachment

    # Ending
    echo "$sample_name complete!"
    echo "================================================="
done

echo "Pipeline finished"