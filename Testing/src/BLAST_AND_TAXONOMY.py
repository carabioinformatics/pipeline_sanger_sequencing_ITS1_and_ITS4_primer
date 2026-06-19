import sys
import time
from datetime import datetime
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
#TODO Make comparison between different databases

class Hit:
    def __init__(self, species, database_type, num_of_hits):
        self.species = species
        self.database_type = database_type #0 - NCBI, 1 - UNITE
        self.num_hits = num_of_hits
        self.arr_bit_scores = []
        self.max_bit_score = 0
        self.arr_coverage = []
        self.arr_identity = []
        self.arr_scores = []
        self.avg_scores = 0
    
    def get_species(self):
        return self.species
    def add_num_hits(self):
        self.num_hits += 1
    def get_num_hits(self):
        return self.num_hits
    def update_average_score(self):
        sum = 0.00
        for i in range(0, len(self.arr_scores)):
            sum += self.arr_scores[i]
        self.avg_scores = sum/self.num_hits
    def get_average_score(self):
        return self.avg_scores
    def add_score(self, temp_score):
        self.arr_scores.append(temp_score)
        self.update_average_score()
    def update_score(self, index, score):
        self.arr_scores[index] = score
        self.update_average_score()
    def add_bit_score(self, temp_bit_score):
        self.arr_bit_scores.append(temp_bit_score)
        max_bit_score = self.get_max_bit_score()
        if (temp_bit_score > max_bit_score):
            self.max_bit_score = temp_bit_score
    def get_max_bit_score(self):
        max = 0
        for i in range(0, len(self.arr_bit_scores)):
            if (self.arr_bit_scores[i] > max):
                max = self.arr_bit_scores[i]
        return max
    def get_bit_score(self, index):
        return self.arr_bit_scores[index]
    def add_coverage(self, temp_coverage):
        self.arr_coverage.append(temp_coverage)
    def add_identity(self, temp_identity):
        self.arr_identity.append(temp_identity)


def main():
    """
        Takes 3 arguments, 
        sys.argv[1] = Name of input consensus file in format INT26_{}_CONSENSUS.fas
        sys.argv[2] = Mode: 0 - Just complete BLAST results and BLAST summary
                            1 - BLAST complete results, BLAST summary and taxonomy counting
        sys.argv[3] = E-value threshold
        sys.argv[4] = String extension added to files that need to be saved in the cleanup directory
        sys.argv[5] = String extension added to files that need to be saved in the final directory
        sys.argv[6] = String extension added to files that need to be saved in the blast directory
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
    comparison_directory_extension = sys.argv[7]
    ######################### Starting values #################################
    i_name_start = 0
    i_name_end = 0
    s_file_identifier = ""
    s_blast_complete_ncbi = ""
    s_blast_summary_ncbi = ""
    s_blast_complete_unite = ""
    s_blast_summary_unite = ""
    s_taxonomy_counter_ncbi = ""
    s_taxonomy_counter_unite = ""
    s_performance_report = ""
    s_comparison_report = ""
    ncbi_start_time = 0
    ncbi_end_time = 0
    unite_start_time = 0
    unite_end_time = 0
    comparison_start_time = 0
    comparison_end_time = 0
    #################### Establish how files will be named #####################
    i_name_start = s_input_file.find("_") + 1
    i_name_end = s_input_file.find("_", i_name_start)
    if ((i_name_start != -1) and (i_name_end != -1)):
        s_file_identifier = s_input_file[i_name_start:i_name_end]
        s_blast_complete_ncbi = blast_directory_extension + "/INT26_" + s_file_identifier + "_BLAST_COMPLETE_NCBI.xml"
        s_blast_summary_ncbi = final_directory_extension + "/INT26_" + s_file_identifier + "_BLAST_SUMMARY_NCBI.txt"
        s_blast_complete_unite = blast_directory_extension + "/INT26_" + s_file_identifier + "_BLAST_COMPLETE_UNITE.xml"
        s_blast_summary_unite = final_directory_extension + "/INT26_" + s_file_identifier + "_BLAST_SUMMARY_UNITE.txt"
        s_taxonomy_counter_ncbi = final_directory_extension + "/INT26_" + s_file_identifier + "_TAXONOMY_COUNTER_NCBI.xml"
        s_taxonomy_counter_unite = final_directory_extension + "/INT26_" + s_file_identifier + "_TAXONOMY_COUNTER_UNITE.xml"
        s_performance_report = final_directory_extension + "/INT26_" + s_file_identifier + "_PERFORMANCE_REPORT.txt"
        s_comparison_report = comparison_directory_extension + "/INT26_" + s_file_identifier = "_COMPARISON_REPORT.txt"
    else: 
        sys.stderr.write("Error: Invalid starting file argument. Does not follow structure of INT26_{}_CONSENSUS.fas\n")
        sys.exit()
        # with open("Error_file.txt", "w") as std_error:
        #     std_error.write("Invalid starting file argument. Does not follow structure of INT26_{}_CONSENSUS.fas")
        # std_error.close()
    ################### Perform BLAST on different databases ##################
    ncbi_start_time = time.perf_counter()
    arr_hits_ncbi = []
    arr_hits_ncbi = ncbi(s_input_file, s_blast_summary_ncbi, s_blast_complete_ncbi, s_taxonomy_counter_ncbi, f_e_value_threshold, i_mode)
    ncbi_end_time = time.perf_counter()
    unite_start_time = time.perf_counter()
    arr_hits_unite = unite(s_input_file, s_blast_summary_unite, s_blast_complete_unite, s_taxonomy_counter_unite, f_e_value_threshold, i_mode)
    unite_end_time = time.perf_counter()
    #################### Make performance report ###############################
    end_time = time.perf_counter()
    with open(s_performance_report, "a") as performance_output:
        performance_output.write("------------------------------ BLASTn report --------------------------------\n")
        performance_output.write(f"Date and time of BLASTn search end:\t{datetime.now():%Y-%m-%d %H:%M}\n")
        performance_output.write(f"Elapsed time for total search: \t{end_time-start_time:.2f} seconds\n")
        performance_output.write(f"Elapsed time for NCBI search: \t{ncbi_end_time-ncbi_start_time:.2f} seconds\n")
        performance_output.write(f"Elapsed time for UNITE search: \t{unite_end_time-unite_start_time:.2f} seconds\n")
        performance_output.write("Input file received: " + s_input_file + "\n")
        performance_output.write("Mode chosen: " + str(i_mode) + "\n")
        performance_output.write(f"E-value threshold: E-value < {f_e_value_threshold:.4f} \n")
        performance_output.write("------------------------------ BLASTn report -------------------------------\n")
        #performance_output.write("Computer where search was run: " + socket.gethostbyname())
    performance_output.close()
    print("CL: BLAST search and taxonomy ending.")
    print("CL: Begin comparison between databases.")
    comparison_start_time = time.perf_counter()
    comparison(comparison_start_time, s_comparison_report, arr_hits_ncbi, arr_hits_unite, f_e_value_threshold)
    print("CL: End comparison between databases.")

def ncbi(s_input_file, s_blast_summary_ncbi, s_blast_complete_ncbi, s_taxonomy_counter_ncbi, f_e_value_threshold, i_mode):
    """
        This method handles all computation related to BLASTn searching on the NCBI database. First does a 
        blastn search and then builds up a summary report that can report the taxonomy found in the search.
    """
    print("CL: Start with BLASTn summary for NCBI.")
    table_tax_count = {}
    counter = 0
    arr_hits_ncbi = []
    obj_hit = None
    # blastn_ncbi(s_input_file, s_blast_complete_ncbi)
    ############ Make summary of BLASTn report using NCBI DATABASE #########
    with open(s_blast_summary_ncbi, "w") as summary_out:
        with open(s_blast_complete_ncbi) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if ((((hsp.identities/hsp.align_length)*100) > 70) and (((hsp.align_length/record.query_length) * 100)> 70)):
                            print("bit score: " + str(hsp.bits))
                            # print("coverage: " + str((hsp.align_length/record.query_length) * 100))
                            # print("identity: " + str((hsp.identities/hsp.align_length) * 100))
                            #(hsp.expect < f_e_value_threshold) and 
                            accession = alignment.accession
                            tax_record = get_taxonomy_ncbi(accession)
                            taxonomy = parse_taxonomy_ncbi(tax_record)
                            print("CL: NCBI taxonomy done")
                            if ((taxonomy != None)):
                                summary_out.write(f"Accession: {alignment.accession}\t\t")
                                summary_out.write(f"Species: {taxonomy['species']}\t\t\t\t")
                                summary_out.write(f"E-value: {hsp.expect}\n")
                            else: 
                                print("taxonomy == NONE")
                            obj_hit = find_hit_obj_in_arr(arr_hits_ncbi, taxonomy['species'])
                            print("finding done")
                            if (obj_hit == None):
                                print("New")
                                obj_hit = add_new_obj(taxonomy['species'], 0, hsp.bits, (hsp.align_length/record.query_length) * 100, (hsp.identities/hsp.align_length)*100)
                                arr_hits_ncbi.append(obj_hit)
                            else:
                                print("Additional")
                                #TODO, don't think this will update the array, only the obj
                                add_additional_entry(obj_hit, hsp.bits, (hsp.align_length/record.query_length) * 100, (hsp.identities/hsp.align_length)*100)
                            if (i_mode == 1):
                                counter += 1
                                table_tax_count = get_taxonomy_count_ncbi(table_tax_count, taxonomy)
    summary_out.close()
    print("CL: Done with BLASTn summary for NCBI.")
    ####################### Sort species with counts #################################
    sorted_table_tax_count = dict(sorted(table_tax_count.items(), key = lambda x: x[1], reverse=True))
    ################### Print NCBI taxonomy counter report ##########################
    if (i_mode == 1):
        with open(s_taxonomy_counter_ncbi, "w") as taxonomy_counter_output:
            taxonomy_counter_output.write("Species name\t\t\t\t\tCount of species hit\t\t'%'of " + str(counter) +" hits for species\n")
            for key, value in sorted_table_tax_count.items():
                taxonomy_counter_output.write(key + "\t\t\t\t" + str(value) + "\t\t\t\t\t" + str(value/counter*100) + "\n")
        print("CL: Done with NCBI taxonomy counter.")
    print(len(arr_hits_ncbi))
    # for a in range(0, len(arr_hits_ncbi)):
    #     for y in range(0, arr_hits_ncbi[a].get_num_hits()):
    #         print(arr_hits_ncbi[a].get_species())
    #         print(arr_hits_ncbi[a].get_bit_score(y))
    return arr_hits_ncbi
        
def unite(s_input_file, s_blast_summary_unite, s_blast_complete_unite, s_taxonomy_counter_unite, f_e_value_threshold, i_mode):
    """
        This method handles all computation related to BLASTn searching on the locally downloaded UNITE database. 
        First does a blastn search on the local database and then builds up a summary report that can report the 
        taxonomy found in the search.
    """
    print("CL: Start with BLASTn summary for UNITE.")
    #blastn_unite(s_input_file, s_blast_complete_unite)
    table_tax_count = {}
    counter = 0
    arr_hits_unite = []
    ############ Make summary of BLASTn report using UNITE DATABASE #########
    with open(s_blast_summary_unite, "w") as summary_out:
        with open(s_blast_complete_unite) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if (hsp.expect < f_e_value_threshold): # and ((hsp.identities/hsp.align_length)*100 > 90) and ((hsp.align_length/record.query_length) *100 > 85)):
                            # e_value and %identity and %coverage
                            species = get_taxonomy_unite(alignment.title)
                            if (species != None):
                                summary_out.write(f"Accession: {alignment.accession}\t\t")
                                summary_out.write(f"Species: {species}\t\t")
                                summary_out.write(f"E-value {hsp.expect}\n")
                            # print("bit score: " + str(hsp.bits))
                            # print("identity: " + str((hsp.align_length/record.query_length) * 100))
                            # print("coverage: " + str((hsp.identities/hsp.align_length) * 100))
                            obj_hit = find_hit_obj_in_arr(arr_hits_unite, species)
                            if (obj_hit == None):
                                obj_hit = add_new_obj(species, 1, hsp.bits, (hsp.align_length/record.query_length) * 100, (hsp.identities/hsp.align_length)*100)
                                arr_hits_unite.append(obj_hit)
                            else:
                                #TODO, don't think this will update the array, only the obj
                                add_additional_entry(obj_hit, hsp.bits, (hsp.align_length/record.query_length) * 100, (hsp.identities/hsp.align_length) * 100)
                            if (i_mode == 1):
                                counter += 1
                                table_tax_count = get_taxonomy_count_unite(table_tax_count, species)
    summary_out.close()
    for t in range(0, len(arr_hits_unite)):
        for a in range(0, arr_hits_unite[t].get_num_hits()):
            score = calculate_score(arr_hits_unite[t].get_bit_score(a), arr_hits_unite[t].get_max_bit_score())
            arr_hits_unite[t].update_score(a, score)
    print("CL: Done with BLASTn summary for UNITE.")
    ####################### Sort species with counts #################################
    sorted_table_tax_count = dict(sorted(table_tax_count.items(), key = lambda x: x[1], reverse=True))
    ################### Print UNITE taxonomy counter report ##########################
    if (i_mode == 1):
        with open(s_taxonomy_counter_unite, "w") as taxonomy_counter_output:
            taxonomy_counter_output.write("Species name\t\t\t\t\tCount of species hit\t\t'%'of " + str(counter) + " hits for species\n")
            for key, value in sorted_table_tax_count.items():
                taxonomy_counter_output.write(key + "\t\t\t\t" + str(value) + "\t\t\t\t\t" + str(value/counter*100) + "\n")
        print("CL: Done with UNITE taxonomy counter.")
    print(len(arr_hits_unite))
    # for a in range(0, len(arr_hits_unite)):
    #     for y in range(0, arr_hits_unite[a].get_num_hits()):
    #         print(arr_hits_unite[a].get_species())
    #         print(arr_hits_unite[a].get_bit_score(y))
    return arr_hits_unite

def blastn_ncbi(s_input_file, s_blast_complete):
    """
        Does the BLASTn search online using the NCBI database.
    """
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
    print("CL: Writing of complete NCBI BLASTn report done.")
    
def blastn_unite(s_input_file, s_blast_complete):
    """
        Does a BLASTn search locally using the UNITE database.
    """
    UNITE_database = "./database/unite/unite_its_database"
    command = ["blastn", "-query", s_input_file, "-db", UNITE_database, "-out", s_blast_complete, "-outfmt", "5"]
    #"-evalue", str(f_e_value_threshold), "-max_target_seqs", "10"]
    
    blast_result = subprocess.run(command, capture_output=True, text=True)
    print("CL: Writing of complete UNITE BLASTn report done.")

def get_taxonomy_ncbi(accession):
    """
        Does an Entrez search to find the taxonomy of the received accession number.
        Returns array tax_records that contains the taxonomy of specificied accession number.
    """
    tax_handle = Entrez.efetch(db="nucleotide", id=accession, retmode="xml")        
    print("1")
    records = Entrez.read(tax_handle)
    print("2")
    tax_handle.close()
    print("first entrez")
    #extracting ID
    #tax_id = records[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_value"]
    features = records[0]["GBSeq_feature-table"]
    if (features == None):
        print("None")
    
    tax_id = None
    for feature in features:
        if "GBFeature_quals" in feature:
            for qual in feature["GBFeature_quals"]:
                if qual["GBQualifier_name"] == "db_xref" and "taxon:" in qual["GBQualifier_value"]:
                    tax_id = qual["GBQualifier_value"].split(":")[1]
    print("Here")
    if tax_id is None:
        raise ValueError("Tax_id not found")

    #fetch taxonomy
    tax_handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
    tax_records = Entrez.read(tax_handle)
    tax_handle.close()
    print("second entrez")
    return tax_records[0]

def parse_taxonomy_ncbi(tax_records):
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

def get_taxonomy_count_ncbi(table_tax_count, taxonomy):
    """
        Gets counter of how much each species match occurs.
        Returns arr_tax_count where you can search for species name and get count.
    """
    cur_value = 0

    if (taxonomy['species']) in table_tax_count:
        cur_value = table_tax_count.get(taxonomy['species'])
        cur_value += 1
        table_tax_count[taxonomy['species']] = cur_value
    else:
        table_tax_count.update({taxonomy['species']: 1})

    return table_tax_count

def get_taxonomy_unite(title):
    """
        Takes the title of the fasta file entry and finds the species name in the string. Returns only the species name.
    """
    title_parts = title.split("|")
    id = None
    taxonomy = None

    for p in title_parts:
        if (p.startswith("UDB")):
            id = p
        if ("s__") in p:
            count = p.find("s__")
            taxonomy = p[count+3:]
    
    return taxonomy

def get_taxonomy_count_unite(table_tax_count, taxonomy):
    """
        Gets counter of how much each species match occurs.
        Returns arr_tax_count where you can search for species name and get count.
    """
    cur_value = 0

    if (taxonomy) in table_tax_count:
        cur_value = table_tax_count.get(taxonomy)
        cur_value += 1
        table_tax_count[taxonomy] = cur_value
    else:
        table_tax_count.update({taxonomy: 1})

    return table_tax_count

def comparison(comparison_start_time, s_comparison_report, arr_hits_ncbi, arr_hits_unite, f_e_value_threshold):
    arr_hits = []
    with open(s_comparison_report, "a") as comparison_output:
        comparison_output.write("------------------------------ Comparison report --------------------------------\n")
        # with open(s_blast_complete_ncbi) as ncbi_input:
        #     with open(s_blast_complete_unite) as unite_input:
        #         ncbi_blast_records = NCBIXML.parse(ncbi_input)
        #         unite_blast_records = NCBIXML.parse(unite_input)
        #         for ncbi_record in ncbi_blast_records:
        #             for ncbi_alignment in ncbi_record.alignments:
        #                 for ncbi_hsp in ncbi_alignment.hsps:
        #                     if ncbi_hsp.expect < f_e_value_threshold: #filter by E-value
        #                         obj_hit = Hit(species, 0, 0)
        comparison_output.write(f"Date and time of BLASTn search end:\t{datetime.now():%Y-%m-%d %H:%M}\n")
        comparison_end_time = time.perf_counter()
        comparison_output.write(f"Elapsed time for total comparison: \t{comparison_end_time-comparison_start_time:.2f} seconds\n")
        comparison_output.write("------------------------------ Comparison report --------------------------------\n")
    comparison_output.close()

def calculate_score(bit_score, max_bit_score):
    score = 0
    score = bit_score/max_bit_score
    return score

def add_new_obj(species, database_type, bit_score, coverage, identity):
    obj_hit = Hit(species, database_type, 1)
    obj_hit.add_bit_score(bit_score)
    obj_hit.add_coverage(coverage)
    obj_hit.add_identity(identity)
    obj_hit.add_score(0)
    return obj_hit

def add_additional_entry(obj_hit, bit_score, coverage, identity):
    obj_hit.add_num_hits()
    obj_hit.add_bit_score(bit_score)
    obj_hit.add_coverage(coverage)
    obj_hit.add_identity(identity)
    obj_hit.add_score(0)
    return obj_hit

def find_hit_obj_in_arr(arr_hits, species):
    obj_hit = Hit("",0,0)
    for b in range(0, len(arr_hits)):
        obj_hit = arr_hits[b]
        if (obj_hit.get_species() == species):
            return obj_hit
    return None

if __name__ == "__main__" : main()