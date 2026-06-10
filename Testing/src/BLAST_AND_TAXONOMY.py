import sys
import time
from datetime import datetime
import socket
import subprocess
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Entrez
from urllib.error import URLError
Entrez.email = "26869993@sun.ac.za"

############################ TODO ################################
#TODO Direct stderr to seperate file, not just terminal
#TODO Make .sp more general

def main():
    """
        Takes 3 arguments, 
        sys.argv[1] = Name of input consensus file in format INT26_{}_CONSENSUS.fas
        sys.argv[2] = Mode: 0 - Just complete BLAST results and BLAST summary
                            1 - BLAST complete results, BLAST summary and taxonomy counting
        sys.argv[3] = E-value threshold
    """
    print("CL: BLAST search and taxonomy starting")
    start_time = time.perf_counter()
    ################### Arguments initialisation ##############################
    s_input_file = sys.argv[1]
    i_mode = int(sys.argv[2])
    f_e_value_threshold = float(sys.argv[3])
    if (f_e_value_threshold == None):
        f_e_value_threshold = 0.001
    cleanup_directory_extension = sys.argv[4]
    final_directory_extension = sys.argv[5]
    blast_directory_extension = sys.argv[6]
    ######################### Starting values #################################
    counter = 0
    i_name_start = 0
    i_name_end = 0
    s_file_identifier = ""
    s_blast_complete = ""
    s_blast_summary = ""
    s_taxonomy_counter = ""
    s_performance_report = ""
    table_tax_count = {}
    #################### Establish how files will be named #####################
    i_name_start = s_input_file.find("_") + 1
    i_name_end = s_input_file.find("_", i_name_start)
    if ((i_name_start != -1) and (i_name_end != -1)):
        s_file_identifier = s_input_file[i_name_start:i_name_end]
        s_blast_complete = blast_directory_extension + "/INT26_" + s_file_identifier + "_BLAST_COMPLETE.xml"
        s_blast_summary = final_directory_extension + "/INT26_" + s_file_identifier + "_BLAST_SUMMARY.txt"
        s_taxonomy_counter = final_directory_extension + "/INT26_" + s_file_identifier + "_TAXONOMY_COUNTER.xml"
        s_performance_report = final_directory_extension + "/INT26_" + s_file_identifier + "_PERFORMANCE_REPORT.txt"
    else: 
        sys.stderr.write("Error: Invalid starting file argument. Does not follow structure of INT26_{}_CONSENSUS.fas\n")
        sys.exit()
        # with open("Error_file.txt", "w") as std_error:
        #     std_error.write("Invalid starting file argument. Does not follow structure of INT26_{}_CONSENSUS.fas")
        # std_error.close()
    ################### Start a BLASTn search using input file ##################
    blastn_unite(s_input_file, s_blast_complete)
    print("CL: BLASTn search complete.")
    ############ Make summary of BLASTn report using recent output file #########
    with open(s_blast_summary, "w") as summary_out:
        with open(s_blast_complete) as result_handle:
            #out_handle.write(result_handle.read())
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < f_e_value_threshold: #filter by E-value
                            accession = alignment.accession

                            tax_record = get_taxonomy(accession)
                            taxonomy = parse_taxonomy(tax_record)
                            if ((taxonomy != None) and (taxonomy['species'] != "Penicillium sp.")):
                                summary_out.write(f"Accession: {alignment.accession}\t\t")
                                #summary_out.write(f"Genus: {taxonomy['genus']}")
                                summary_out.write(f"Species: {taxonomy['species']}\t\t\t\t")
                                summary_out.write(f"E-value: {hsp.expect}\n")
                            #print(f"Match: {alignment.title}")
                            #print(f"E-value: {hsp.expect}")
                            if (i_mode == 1):
                                counter += 1
                                table_tax_count = get_taxonomy_count(table_tax_count, taxonomy)

    #print("Counter: " + str(counter))
    #summary_out.write("Counter: " + str(counter))
    summary_out.close()
    print("CL: Done with BLASTn summary document.")
    ################### Print taxonomy counter report ##########################
    if (i_mode == 1):
        with open(s_taxonomy_counter, "w") as taxonomy_counter_output:
            taxonomy_counter_output.write("Species name\t\t\t\t\tCount of species hit\t\t'%'of hits for species\n")
            for key, value in table_tax_count.items():
                taxonomy_counter_output.write(key + "\t\t\t\t" + str(value) + "\t\t\t\t\t" + str(value/counter) + "\n")
        print("CL: Done with taxonomy counter.")
    #################### Make performance report ###############################
    end_time = time.perf_counter()
    with open(s_performance_report, "a") as performance_output:
        performance_output.write("------------------------------ BLASTn report --------------------------------\n")
        performance_output.write(f"Date and time of BLASTn search end:\t{datetime.now():%Y-%m-%d %H:%M}\n")
        performance_output.write(f"Elapsed time for BLASTn search: \t{end_time-start_time:.2f} seconds\n")
        performance_output.write("Input file received: " + s_input_file + "\n")
        performance_output.write("Mode chosen: " + str(i_mode) + "\n")
        performance_output.write(f"E-value threshold: E-value < {f_e_value_threshold:.4f} \n")
        performance_output.write("------------------------------ BLASTn report -------------------------------\n")
        #performance_output.write("Computer where search was run: " + socket.gethostbyname())
    performance_output.close()
    print("CL: BLAST search and taxonomy ending.")

def blastn_ncbi(s_input_file, s_blast_complete):
    record = SeqIO.read(s_input_file, format="fasta")
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
    except URLError as e:
        sys.stderr.write("Error: Connection error occured\n")
        sys.exit()
    except Exception as e:
        sys.stderr.write(f"Error: {e}")
        sys.exit()
    ############## Write complete BLASTn report to output file ##################
    with open(s_blast_complete, "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()
    print("CL: Writing of complete BLASTn report done.")
    
def blastn_unite(s_input_file, s_blast_complete):
    UNITE_database = "./database/unite/unite_its_database"
    command = ["blastn", "-query", s_input_file, "-db", UNITE_database, "-out", s_blast_complete, "-outfmt", "5"]
    #"-evalue", str(f_e_value_threshold), "-max_target_seqs", "10"]
    
    blast_result = subprocess.run(command, capture_output=True, text=True)
    print("return code: " + str(blast_result.returncode))
    print("STDERR: " + blast_result.stderr)

def get_taxonomy(accession):
    """
        Does an Entrez search to find the taxonomy of the received accession number.
        Returns array tax_records that contains the taxonomy of specificied accession number.
    """
    tax_handle = Entrez.efetch(db="nucleotide", id=accession, retmode="xml")
    records = Entrez.read(tax_handle)
    tax_handle.close()

    #extracting ID
    #tax_id = records[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_value"]
    features = records[0]["GBSeq_feature-table"]
    
    tax_id = None
    for feature in features:
        if "GBFeature_quals" in feature:
            for qual in feature["GBFeature_quals"]:
                if qual["GBQualifier_name"] == "db_xref" and "taxon:" in qual["GBQualifier_value"]:
                    tax_id = qual["GBQualifier_value"].split(":")[1]

    if tax_id is None:
        raise ValueError("Tax_id not found")

    #fetch taxonomy
    tax_handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
    tax_records = Entrez.read(tax_handle)
    tax_handle.close()

    return tax_records[0]

def parse_taxonomy(tax_records):
    """
        Organise data from tax_records to be more readible.
        Returns taxonomy array that you can search using keywords like 'species'
    """
    linaeage = tax_records["LineageEx"]

    taxonomy = {
        "species": tax_records.get("ScientificName", "NA"),
        "genus": "NA",
        "family": "NA",
        "order": "NA",
        "class": "NA",
        "phylum": "NA",
        "kingdom": "NA"
    }

    for rank in linaeage:
        if rank["Rank"] in taxonomy:
            taxonomy[rank["Rank"]] = rank["ScientificName"]

    return taxonomy

def get_taxonomy_count(table_tax_count, taxonomy):
    """
        Gets counter of how much each species match occurs.
        Returns arr_tax_count where you can search for species name and get count.
    """
    cur_value = 0
    new_hash = 0

    if (taxonomy['species']) in table_tax_count:
        cur_value = table_tax_count.get(taxonomy['species'])
        cur_value += 1
        table_tax_count[taxonomy['species']] = cur_value
    else:
        table_tax_count.update({taxonomy['species']: 1})

    return table_tax_count
if __name__ == "__main__" : main()