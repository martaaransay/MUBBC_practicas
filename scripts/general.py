# Import functions from the pipeline module 
from pipeline.downloadgenomes import datasets_flags, download_genomes
from pipeline.generate_db import download_species_genomes, generate_multifasta, generate_db
from pipeline.useful_funct import contig_to_genome
from pipeline.parallel_download import process_all
from pipeline.blastn import run_blastn, process_hits
from pipeline.analyze_Pc import classify_pc_regex, info_to_csv


if __name__ == "__main__":
    n_genomes = 100000
    
    db_created = False
    if not db_created:
        flags = datasets_flags(assembly_level="complete,chromosome,contig",
                               exclude_atypical=True,
                               mag=False)
    
        # species to download
        # names = {"Ecoli" : "Escherichia coli", 
        #         "Klebsiella" : "Klebsiella pneumoniae",
        #         "Pseudomonas" : "Pseudomonas aeruginosa",
        #         "Acinetobacter" : "Acinetobacter baumannii",
        #         "Enterobacter" : "Enterobacter"}
        
        ecoli = {"Ecoli" : "Escherichia coli"}
        ## For downloading normal genomes --> download_genomes() directly
        genomes_dir = f"../data_practicas/Ecoli_complete"
        dirs = download_species_genomes(species = ecoli, 
                                        outdir = genomes_dir, 
                                        flags = flags,
                                        target_num = n_genomes)
        
        failed = process_all(dirs, total_workers=30, n_split=500)
        if failed: 
            print(f"Failed to rehydrate: {failed}")
            exit(0)

    multifasta_path = f"../data_practicas/Ecoli_complete/multifasta_genomes.fasta"
       
    generate_multifasta(path_to_genomes = f"../data_practicas/Ecoli_complete/*/*/ncbi_dataset/data/*/*.fna",
                        outfile = multifasta_path)
    generate_db(multifasta = multifasta_path,
                        out_db = f"../data_practicas/Ecoli_complete/db/db")

    
    query_file = "data_practicas/Int1/Int1_Pc_pattern.fa"


    extra_columns = {"bioproject": "assemblyInfo.bioprojectLineage.bioprojects.accession",
                     "biosample": "assemblyInfo.biosample.accession",
                     "organism_name": "organism.organismName",
                     "assembly_level": "assemblyInfo.assemblyLevel"}
    
    result, t = run_blastn(query = query_file, 
                           output = f"../results/final/Ecoli_complete/blast/blast.tsv",
                           db = f"../data_practicas/Ecoli_complete/db/db",
                           max_target = len(contig_to_genome(genomes_dir))) # max num of seqs to keep: - mapping length

    hits_df = process_hits(result,
                           dir_genomes = f"../data_practicas/Ecoli_complete",
                           multifasta = multifasta_path,
                           min_length = 1000,
                           extra_columns = extra_columns)

    for hit in hits_df.to_dict("records"):
        # Classify each sequence to a defined Pc variant
        pc_info = classify_pc_regex(hit["sequence"],
                                    hit_start = hit["start"],
                                    hit_strand = hit["strand"])
        
        # Add info for: genome, seq_id and species_folder
        row_info = {"genome": hit["genome"],   
                    "contig": hit["seq_id"],  
                    "species": hit["species_folder"]}

        row_info.update(pc_info) # Add the info about the classification of Pc

        
        for col in extra_columns.keys(): # Add extra cols
            row_info[col] = hit[col]

        info_to_csv(row_info, f"../results/final/Ecoli_complete/final.csv")


