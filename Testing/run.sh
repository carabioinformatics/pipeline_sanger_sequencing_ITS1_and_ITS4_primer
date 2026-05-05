#!/bin/bash

#window_size = $1
#quality_threshold = $2
#taxonomy_mode = $3
#e_value_threshold = $4

python3 ./src/TRIMMING_AND_CONSENSUS.py "D05_CLS3_ITS1___1.ab1" "G05_CLS3_ITS4___1.ab1" $1 $2
python3 ./src/BLAST_AND_TAXONOMY.py "INT26_CLS3_CONSENSUS.fas" $3 $4