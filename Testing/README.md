########### Author and project details ##########
Author: Cara Louw
Supervisor: Dr Annelise Botes at Stellenbosch University

Project: 
Process Sanger sequence data of the ITS region from fungal samples to identify which species was in the sample. The forward and reverse reading of the regions are first trimmed to remove low quality regions using a sliding window algorithm. It is then aligned using ClustalO. A consensus sequence of the aligned region is used in a BLASTn search. This program performs the BLASTn search and then makes a summary of the most relevant results. This can then be used to identify what the species of the sequenced sample is. This is used to make the process more consistent, accurate and efficient.

Sample preparation and processing: 
Samples were collected from different Ostrich farms and egg incubation rooms and were then cultured so that fungal species could grow. The samples underwent DNA extraction and then a PCR reaction using ITS1 (forward) and ITS4 (reverse) primers to confirm that fungal DNA can be found in the sample. If there was positive amplification, the PCR reaction is repeated. The PCR product of about 600 bp is purified and then sent to CAF Stellenbosch Campus for Sanger sequencing. 

############# Directory structure ############
~/pipeline_sanger_sequencing_ITS1_and_ITS4_primer/
|
Testing/
|
|-- run.sh
|-- build_database.sh
|-- README.md
|-- .gitignore
|-- src
|   |-- BLAST_AND_TAXONOMY.py
|   |-- REVERSE_COMPLEMENT_AND_ALIGNMENT.py
|   |-- TRIMMING_AND_CONSENSUS.py
|-- Input (needs to be created by user)
|   |-- *.ab1 files (forward and reverse files)
|-- Output (gets created by run.sh file)
|   |-- Sample1 (created by run.sh for each sample)
|       |-- BLAST
|       |-- Cleanup
|       |-- Final
|   |-- Sample2 (created by run.sh for each sample)
|       |-- BLAST
|       |-- Cleanup
|       |-- Final

################ How to run ##################
- In root directory
- ./run.sh window_size quality_threshold_for_trimming taxonomy_mode e_value_threshold
- What I usually use: ./run.sh 5 20 0 0.001
- What each parameter means
    Window size:
        This is used when trimming the raw sequences. The raw sequences is split into windows with this size. The average phred score of this window is taken. This this average score is not good, then this region is considered for trimming. This is more accurate than looking at individual nucleotide phred scores.
    Quality threshold for trimming:
        Above what treshold the phred score should be to not be trimmed off.
    Taxonomy mode:
        Mode 0: Does complete BLASTn search and makes a summary of the results.
        Mode 1: Does a complete BLASTn search and makes a summary of the results. But also then counts how many hits were found for each species.
    E-value threshold:
        Instead of downloading the entire results list for a BLASTn search, only searches with a lower E-value than the threshold is saved in this report.

########### Download requirements #################
- Python 3.10 worked best for me
- Biopython (pip install biopython)
- Miniconda (https://www.anaconda.com/docs/getting-started/miniconda/install/windows-gui-install)
- ClustalO ()
- UNITE2024ITS.fasta file fingerprint: md5:a2deecb84d0f322a7cde137a1c8c67f3
For updating UNITE database: 
(Replaces all old database files with updated ones. Code can run normally again)
    - Download the new FASTA
    - Run this in the bash file:
    makeblastdb \
        -in new_unite_release.fasta \
        -dbtype nucl \
        -out unite_its

############### Database structure #################
- In terminal, outside conda, execute
    mkdir -p ~/databases/unite
    cd ~/databases/unite
    Download the UNITE FASTA file into this directory
    Run ./build_databash.sh script only at the start to build database
- Structure after building
~/databases/unite/
|
|-- unite2024ITS.fasta
|-- unite_its.nhr
|-- unite_its.nin
|-- unite_its.nsq
|-- ...