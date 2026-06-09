echo "Hello"
echo "Start database building script"
makeblastdb -in /home/caral/c/c/databases/unite/unite2024ITS.fasta -dbtype nucl -out /home/caral/c/c/databases/unite/unite_its_database
#"C:\C\databases\unite\unite2024ITS.fasta"
#blastn -query query.fasta -db ~/c/c/databases/unite/unite_its_database -out results.xml -outfmt 5 -evalue e_value_threshold
echo "End database building script"