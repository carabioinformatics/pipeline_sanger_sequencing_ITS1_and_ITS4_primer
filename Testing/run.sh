#!/bin/bash

echo "Start main run script"
#./build_database.sh

mkdir -p Output
#window_size = $1
#quality_threshold = $2
#taxonomy_mode = $3
#e_value_threshold = $4
timing_report="./Output/INT26_TIME_REPORT.txt"
echo - "Sample name\t\tCleanup time (s)\t\tBLAST time (s)\t\tTotal time (s)" > "$timing_report"

for forward in ./Input/*ITS1___1.ab1
do
    total_start_time=$(date +%s.%N)
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
    mkdir -p ./Output/Comparison
    blast_attachment="./Output/$sample_name/BLAST"
    cleanup_attachment="./Output/$sample_name/Cleanup"
    final_attachment="./Output/$sample_name/Final"
    comparison_attachment="./Output/$sample_name/Comparison"

    # Run scrips now
    cleanup_start_time=$(date +%s.%N)
    python3 ./src/TRIMMING_AND_CONSENSUS.py $forward $reverse $1 $2 $cleanup_attachment $final_attachment $blast_attachment
    cleanup_end_time=$(date +%s.%N)
    #blastn -query $1 -db ./database/unite/unite_its_database -out $2 -outfmt 5
    BLAST_start_time=$(date +%s.%N)
    python3 ./src/BLAST_AND_TAXONOMY.py ./Output/$sample_name/BLAST/$consensus_name $3 $4 $cleanup_attachment $final_attachment $blast_attachment $comparison_attachment
    BLAST_end_time=$(date +%s.%N)

    # Timing report printing
    cleanup_time=$(awk "BEGIN {print $cleanup_end_time - $cleanup_start_time}")
    BLAST_time=$(awk "BEGIN {print $BLAST_end_time - $BLAST_start_time}")
    total_end_time=$(date +%s.%N)
    total_time=$(awk "BEGIN {print $total_end_time - $total_start_time}")
    echo - "$sample_name\t\t$cleanup_time\t\t$BLAST_time\t\t$total_time" >> "$timing_report"

    # Ending
    echo "$sample_name complete!"
    echo "================================================="
done

echo "Pipeline finished"
pytho3 TIME_REPORT_PLOTTING.py $timing_report
echo "End main run script"