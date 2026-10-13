# Import libraries 
import os
import re
import csv
from Bio import SeqIO
from Bio.Seq import Seq


# --------Definition of functions
#------------------------ Function to extract information about integrons from a gbk from IntegronFinder2
def extract_integron_info_gbk(genome_acc, gbk_files):
    """
    Extracts integron information from .gbk files.

    Args:
    ------
        genome_acc (str): The accession number of the genome.
        gbk_files (list): A list containing the file paths to the .gbk files.

    Return:
    ------
        integron_info (list of dict): A list of dictionaries containing the extracted information. 
            For CALIN integrons, the dictionary contains:
                - "acc" (str): Genome accession number.
                - "calin" (bool): True.
                - "integron_id" (str): The ID of the integron.
                - "integron_type" (str): The type of the integron.
                - "locus_id" (str): The replicon/locus ID.
            For non-CALIN integrons, it also includes:
                - "calin" (bool): False.
                - "start" (int): Start position of the integrase feature.
                - "end" (int): End position of the integrase feature.
                - "strand" (int): Strand of the integrase feature (1 for +, -1 for -).
                - "sequence" (Bio.Seq.Seq): The extracted sequence of the integrase.
    """
    # List to store integron information
    integron_info = []
    # Iterate through the .gbk files
    for gbk_file in gbk_files: # (more than one file can be generated from 1 genome)
        record = SeqIO.read(gbk_file, "genbank") # Read .gbk file
        locus_id = record.id # Store replicon id

        integrons = [] # list to store integron information for the current record
        # Two iterations:
        # 1) To iterate through integrons
        for feature1 in record.features:
                    if feature1.type == "integron":
                        # Get integron id (in case no value, empty list)
                        # [0] to get the first value of the list
                        integron_id = feature1.qualifiers.get("integron_id", [""])[0]
                        # Get integron type
                        type_integron = feature1.qualifiers.get("integron_type", [""])[0]
                        start = int(feature1.location.start) # Start of Int
                        end = int(feature1.location.end) # End of Int
                        # Store the integron information into the dictionary
                        integrons.append({"id": integron_id,
                                          "type": type_integron,
                                          "start": start,
                                          "end": end})
                        # If the last integron is a CALIN (no integrase), store its information
                        if integrons[-1]["type"] == "CALIN":
                            # Store info about calin-integron
                            integron_info.append({"acc": genome_acc, 
                                                  "integron_id": integrons[-1]["id"], # Last id
                                                  "integron_type": integrons[-1]["type"], # Last type
                                                  "locus_id": locus_id,
                                                  "calin": True})
         # 1) To iterate through integrases
        for feature2 in record.features:
            if feature2.type != "integrase": # If the feature is not a integrase, continue
                continue

            start = int(feature2.location.start) # Start of Integrase
            end = int(feature2.location.end) # End of Integrase

            # If the integron corresponds to the current integrase, takes the integron information
            curr_integron = None
            for each_integron in integrons:
                if each_integron["start"] < end and start < each_integron["end"]:
                    curr_integron = each_integron # Updates the variable with the current integron information
                    break  

            if not curr_integron: # Control to avoid errors in case no integron is found for the integrase
                curr_integron = {"id": "", "type":""}
            
            strand = feature2.location.strand # Strand of Integrase
            seq = feature2.extract(record.seq) # Sequence of Integrase
            # Store info about integron with integrase
            integron_info.append({"acc": genome_acc, 
                                  "integron_id": curr_integron["id"],
                                  "integron_type": curr_integron["type"],
                                  "locus_id": locus_id,
                                  "calin": False,
                                  "start": start, 
                                  "end": end,         
                                  "strand": strand,
                                  "sequence": seq})   
        
    return integron_info

# Regular Expressions of the defined pc variants
pc_variants_regex = {
    "PcS": ("TTGACA[ACGT]{17}TAAACT"),
    "PcW": ("TGGACA[ACGT]{17}TAAGCT"),
    "PcWTGN-10": ("TGGACA[ACGT]{14}TG[ACGT]TAAGCT"),
    "PcH1": ("TGGACA[ACGT]{17}TAAACT"),
    "PcH1TGN-10": ("TGGACA[ACGT]{14}TG[ACGT]TAAACT"),
    "PcH2": ("TTGACA[ACGT]{17}TAAGCT"),
    "PcH2TGN-10": ("TTGACA[ACGT]{14}TG[ACGT]TAAGCT"),
    "TGGGCA-N14-TGn-TAAGCT": ("TGGGCA[ACGT]{14}TG[ACGT]TAAGCT"),
    "PcSS": ("TTGATA[ACGT]{17}TAAACT"),
    "Pcln42": ("TTGGCA[ACGT]{17}TAAACT"),
    "Pcln116": ("TTGACA[ACGT]{17}TGAACT"),
    "PcPUO": ("TCGACA[ACGT]{17}TAAACT"),
    "P2": ("TTGTTA[ACGT]{17}TACAGT"),
    "P2m1":("TTGTTA[ACGT]{17}GACAGT"),
    "P2m2": ("TTGTTA[ACGT]{17}TACACA")
}

#------------------------ Function to search and classify Pc in Integrase sequences
def classify_pc_regex(sequence, variants = pc_variants_regex,
                      hit_start = None, hit_strand = None):
    """
    Searches and classifies Pc promoter variants within an integrase sequence 
    using regular expressions. 
    
    Args:
    ------
    sequence (str): The nucleotide sequence of the integrase to be analyzed.
    variants (dict, optional): A dictionary containing:
                                -key: names of the Pc variant.
                                -value: their corresponding regular expression patterns. 
                               Defaults to pc_variants_regex.
    hit_start (int, optional): Start of the BLAST hit in the contig (1-based).
    hit_strand (str, optional): Strand of the BLAST hit ("+" or "-").
                                If both are given, genomic coordinates are added.
    
    Return:
    ------
    results (dict): A dictionary containing the classification results of the matched Pc variant. 
                    Returns None if no match is found. 
    """
    # The Pc is in the opposite direction of Integrase, therefore the reverse complement is analyzed
    rc_seq = Seq(sequence).reverse_complement()
    seq = str(rc_seq) # Sequence to string

    # Initialize a variable to store results
    results = None 
    # Iterates over the name of the Pc and the RE of the sequence
    for name_pc_variant, pc_sequence in variants.items():
        match = re.search(pc_sequence, seq) # First match of this variant
        if match: # Later variants overwrite earlier ones 
            #### (keep TGN-10 after its generic variant)
            results = {"Pc_variant": name_pc_variant,
                       "start": match.start(),
                       "end": match.end(),
                       "sequence": match.group()}
        # If no variant matched, NA
        if not results:
            results = {"Pc_variant": "NA", 
                       "start": "NA", 
                       "end": "NA", 
                       "sequence": seq}
    if hit_start is not None:
        if hit_strand is not None:
            # Update coordinates to genomic
            results = update_genomic_coords(results, 
                                            hit_start, hit_strand, len(sequence))
    # Information to keep
    cols_to_keep = ["Pc_variant", "genome_start", "genome_end", "sequence"]
    pc_info = {} # Dictionary to store final results
    for key, value in results.items():
        if key in cols_to_keep:
            pc_info[key] = value

    return pc_info

#------------------------ Function to update the genomic coordinates of the Pc
def update_genomic_coords(pc_info, hit_start, hit_strand, seq_len):
    """
    Adds the genomic coordinates of the Pc (1-based, inclusive) to pc_info.

    Args:
    ------
    pc_info (dict): Output of classify_pc_regex, with "start" and "end" relative
                    to the reverse complement of the analyzed sequence.
    hit_start (int): Start of the BLAST hit in the contig (1-based, lowest coordinate).
    hit_strand (str): Strand of the BLAST hit ("+" or "-").
    seq_len (int): Length of the analyzed sequence.

    Return:
    ------
    pc_info (dict): The same dictionary with "genome_start" and "genome_end" added
                    ("NA" if no Pc was found).
    """
    # If no Pc found, no coordinates to update
    if pc_info["start"] == "NA": 
        pc_info["genome_start"] = "NA"
        pc_info["genome_end"] = "NA"
        return pc_info

    # If the sequence is in the + strand, it was rc previously
    if hit_strand == "+":
        # start = start position on the genome  + length - end of the match
        pc_info["genome_start"] = hit_start + seq_len - pc_info["end"]
        # end = start position on the genome  + length - start of the match
        pc_info["genome_end"] = hit_start + seq_len - pc_info["start"] - 1
        # -1 to fix 0-based
    
    else: # If the sequence is in the - strand, is in the original orientation
        # start = start position on the genome + start
        pc_info["genome_start"] = hit_start + pc_info["start"]
        # end = end position on the genome + end -1 (0-based)
        pc_info["genome_end"] = hit_start + pc_info["end"] - 1
    return pc_info

#------------------------ Function to store information from dictionaries into a csv file
def info_to_csv(dict_info, path_out_file):
    """
    Appends the information from a dictionary into a CSV file. 
    
    Args:
    ------
    dict_info (dict): A dictionary containing the data to be stored. 
    path_out_file (str): The file path where the CSV file will be saved.
    """
    # Open the file in appending mode
    with open(path_out_file, "a", newline="") as out_csv:
        csvwriter = csv.writer(out_csv) # Writer
        # If the file is new, writes as header the keys of the dictionary
        if os.path.getsize(path_out_file) == 0:
            csvwriter.writerow(dict_info.keys())
        # Writes the values of the dictionary
        csvwriter.writerow(dict_info.values())
    return

#------------------------ Function to merge general and Pc information into a single dictionary
def adjust_Pc_info(general_info, # From extract_integron_info_gbk() function
                   Pc_info, # From classify_pc_regex() function
                   path_out_file, flag_Pc_seq = True): # Flag to store the sequence of the Pc 
    """
    Merges general integron information and Pc variant classification details into a single 
    dictionary and appends the resulting data to a CSV file.

    Args:
    ------
    general_info (dict): A dictionary containing general integron information.
                         It must include the keys "acc", "integron_id", and "locus_id".
    Pc_info (dict): A dictionary containing the classification results of the Pc variant. 
                    It must include the keys "start", "end", and "Pc_variant", and 
                    "sequence" if the sequence is to be stored.
    path_out_file (str): The file path where the merged information will be saved as a CSV.
    flag_Pc_seq (bool, optional): A flag indicating whether to include the exact Pc sequence. 
                                  Defaults to True.
    """
    # Extract values from general_info dictionary
    acc_genome = general_info["acc"]
    integron_id = general_info["integron_id"]
    locus_id = general_info["locus_id"]

    # Extract values from Pc_info dictionary
    start = Pc_info["start"]
    end = Pc_info["end"]
    Pc_variant = Pc_info["Pc_variant"]

    # Final dictionary with all the information
    final_results = {"acc_genome" : acc_genome,
                     "integron_id" : integron_id,
                     "locus_id" : locus_id,
                     "start_Pc": start,
                     "end_Pc": end,
                     "Pc_variant": Pc_variant}
    
    # If flag is True, the sequence of the Pc is added
    if flag_Pc_seq:
        final_results["sequence"] = Pc_info["sequence"]

    # Stores the final dictionary into the csv file
    info_to_csv(final_results, path_out_file)
    return 

#------------------------ Function to sort general information into a dictionary
def adjust_calin_info(general_info, # From extract_integron_info_gbk() function
                      path_out_file):
    """
    Formats and extracts specific general information for CALIN integrons 
    into a new dictionary, and appends the resulting data to a CSV file.

    Args:
    ------
    general_info (dict): A dictionary containing general integron information. 
                         It must include the keys "acc_genome", "integron_id", and "locus_id".
    path_out_file (str): The file path where the extracted information will be saved as a CSV.
    """
    
    # Extract values from general_info dictionary
    acc_genome = general_info["acc"]
    integron_id = general_info["integron_id"]
    locus_id = general_info["locus_id"]

    # Final dictionary with all the information
    final_results = {"acc_genome" : acc_genome,
                     "integron_id" : integron_id,
                     "locus_id" : locus_id}

    # Stores the final dictionary into the csv file
    info_to_csv(final_results, path_out_file)

    return 

