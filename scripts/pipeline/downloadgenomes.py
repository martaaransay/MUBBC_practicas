# Import libraries 
import subprocess 
import shlex 
import json 
import os 
import random

# --------Definition of functions
#------------------------ Function to define flags for genome download (datasets download)
def datasets_flags(assembly_source = None, # RefSeq | GenBank | None =all
                   assembly_level = None, # comma-separated (string): chromosome | complete | contig | scaffold | None =all
                   annotated = None, # To exclude annotated genomes  (≠ None)
                   exclude_atypical = None, # Excludes atypical assemblies (≠ None) 
                   exclude_multi_isolate = None, # Excludes assemblies from multi-isolate projects  (≠ None)
                   mag = None, # If True or "only", only uses MAGS; if false or "exclude", excludes them
                   ):
    """
    Builds and validates a list of flags for the NCBI datasets download command.

    Args
    --------
    assembly_source (str): Assembly source. 
                           Allowed values: 'RefSeq' or 'GenBank' (eg. assembly_source='RefSeq').
    assembly_level (str): Assembly level - comma-separated without blank.
                          Allowed values: 'chromosome', 'complete', 'contig', 'scaffold' (eg. assembly_source='chromosome,complete').
    annotated (bool): If evaluated to ≠ False/None, includes annotated genomes ('--annotated').
    exclude_atypical (bool): If evaluated to ≠ False/None, excludes atypical assemblies ('--exclude-atypical').
    exclude_multi_isolate (bool): If evaluated to ≠ False/None, excludes multi-isolate projects ('--exclude-multi-isolate').
    mag (bool | str): Filter for metagenomic genomes. 
                      Allowed values: True/'only', False/'exclude' or 'all' (eg. mag='only').

    Return
    --------
    flags (list): List of strings with the arguments to be passed to subprocess.
    """
    
    flags = [] # List of flags to fill

    if assembly_source: # If an argument was provided, adds it to the flags list
        # datasets has a checkpoint control in case the argument is not valid
        flags += ["--assembly-source", assembly_source] # Add the argument to the list
        
    if assembly_level: 
        # datasets has a checkpoint control in case the argument is not valid
        flags += ["--assembly-level", assembly_level]

    if annotated:
        flags.append("--annotated")

    if exclude_atypical:
        flags.append("--exclude-atypical")

    if exclude_multi_isolate:
        flags.append("--exclude-multi-isolate")

    if mag is not None: 
        if mag is True:
            flags += ["--mag", "only"]
        elif mag is False:
            flags += ["--mag", "exclude"]
        else:
            flags += ["--mag", mag]

    return flags


#------------------------ Function to define taxon/accesion 
def taxon_acc_flag(source, # taxon | accession | inputfile
                   name # List of values (source="taxon" | source="accession")\
                        # or file name (source="inputfile")
                   ):
    """
    Builds the argument list (taxon or accession info) for the NCBI datasets download command.

    Args
    --------
    source (str): Input type for the search. 
                  Allowed values: 'taxon', 'accession' or 'inputfile'.
    name (str): Values to search (comma-separated string) or file path if the source is 'inputfile'.

    Return
    --------
    source_flags (list): List of strings with the unpacked source arguments for datasets.
    """
    # Control of valid arguments (inputfile to internal control is added)
    if source not in ["taxon", "accession", "inputfile"]: 
        raise ValueError(f"One of the valid options was not provided: taxon | accession | inputfile \
                        \n (provided: {source})")
    
    if source == "inputfile": # File
        if not os.path.isfile(name): # Check in case the file does not exist
            raise FileNotFoundError(f"The file {name} was not found.")
            
        flag = ["accession", "--inputfile", name] # Flags to add with the name of the file
    
    else:
        names_flag = [] 
        for each_name in name.split(","): # Single element or several comma-separated
            names_flag.append(each_name)

        flag = [source, *names_flag] # * unpacks, so it does not generate a nested list
    return flag

#------------------------ Function to filter genomes
def quality_filter_genomes(source, name, # Arguments for the taxon_acc_flag(source, name) function
                           filename, # Filename to keep the same directory
                           log_dir,
                           min_n50 = 2000, 
                           max_contamination = 5, 
                           min_completeness = 90,
                           target_num = None,
                           seed = 1,
                           flags = None): 
    """
    Queries genome metadata in NCBI datasets, filters by quality metrics, and saves valid accessions in a text file.
    
    Args
    ------------
    source (str): Source type for the query.
    name (str): Value or identifier associated with the source.
    filename (str): Name of the resulting .zip file. 
    log_dir (str): Directory where logs files will be saved.
    min_n50 (int, optional): Minimum acceptable value for the scaffold N50. 
                                   Default is 2000.
    max_contamination (float | int, optional): Maximum tolerated contamination percentage according to CheckM. 
                                               Default is 5%.
    min_completeness (float | int, optional): Minimum required completeness percentage according to CheckM. 
                                               Default is 90%.
    target_num (int, optional): Maximum number of genomes to download.
                                Default is None (no limit).
    seed (int, optional): Seed for random sampling of genomes if target_num is specified.
    flags (list, optional): List of additional flags to add to the 'datasets summary' command. 
                            Default is None.

    Return
    --------------
    accessions_tuple (tuple): Tuple with the generated file name if there are genomes that pass the filters
                              None if there are no valid results.
    """
    print("Filtering by quality metrics genomes.")

    # Path to a log file (to store errors)
    log_file = os.path.join(log_dir, "datasets_summary_errors.log")

    if not flags: # If no flags were added, empty list
        flags = []

    # To get the JSON with the genome data 
    summary_cmd = ["datasets", "summary", "genome", 
                    *taxon_acc_flag(source, name), *flags, "--as-json-lines"]
    
    print(f"Executing: \n{shlex.join(summary_cmd)}")
    valid_accs = [] # List to store valid accessions 
    control_stop = False # Flag to control if the execution is early stopped

    try:
        with open(log_file, "a") as file:
            # Flush the buffer so the header is written to the file before
            file.write("-------------- datasets summary log --------------\n")
            file.flush()

            # Start the process 
            process = subprocess.Popen(summary_cmd, 
                                    stdout = subprocess.PIPE, # Capture the output
                                    stderr = file,
                                    text = True) # Returns the output as text instead of bytes
            
            for line in process.stdout: # Reads the output line by line
                if not line.strip(): continue # Skips empty lines
                datos = json.loads(line) # Deserialization
                acc = datos.get("accession") # Get the accession
                    
                # Quality metrics:
                #------------------- N50
                assembly_stats = datos.get("assembly_stats", {}) # Get stats (if not present, empty dict)
                n50 = assembly_stats.get("scaffold_n50", 0)
                
                #------------------- CheckM info
                checkm_info = datos.get("checkm_info", {})
                # If there is no checkm info, the default values are set
                # Contamination 
                contamination = checkm_info.get("contamination", 5) 
                # Completeness
                completeness = checkm_info.get("completeness", 90) 

                # Check if it meets the selected filters 
                if n50 >= min_n50 and contamination <= max_contamination and completeness >= min_completeness:
                    valid_accs.append(acc)
                    # If the number of accessions exceeds twice the target number, terminate the process to avoid unnecessary processing
                    if target_num and len(valid_accs) >= (target_num * 2): # twice to allow random sampling
                        control_stop = True 
                        process.terminate() # Stop the process
                        break

            process.wait() # Wait for the process to finish
        if process.returncode != 0 and not control_stop:
            print(f"Error executing datasets summary.\nDetails in {log_file}.")
            return
        
    except Exception as e: 
        print(f"Error executing datasets summary because python failed.\nDetails:\n{e}.")
        return

    if not valid_accs:
        print("No genome meets the quality criteria.")
        return
    
    if target_num: 
        if len(valid_accs) > target_num: # Limit the number of genomes
            random.seed(seed) # Set seed
            valid_accs = random.sample(valid_accs, target_num) # Random sampling
        elif len(valid_accs) < target_num: # Not enough genomes
            print(f"Only {len(valid_accs)} genomes met the quality criteria ")

    # Store the valid accessions into a text file
    accs_file_name = f"accessions_{os.path.splitext(os.path.basename(filename))[0]}.txt" # Name of the file
    accs_file_path = os.path.join(os.path.dirname(filename), accs_file_name) # Path of the file

    with open(accs_file_path, "w") as f: # Create a file with the valid accs (one per line)
        for each_acc in valid_accs:
            f.write(f"{each_acc}\n")
        print(f"The accessions to be analyzed have been saved in the file '{accs_file_path}'.")

    accessions_tuple = ("inputfile", accs_file_path)

    return accessions_tuple

##------------------------ Function to unzip genomes
def unzip(log_dir, filename = "ncbi_dataset.zip"): 
    """
    Unzips an NCBI ZIP file.

    Args
    -------
    log_dir (str): Directory where log files will be saved.
    filename (str, optional): Path to the ZIP file to decompress. 
                              Default is "ncbi_dataset.zip".
    
    Return
    -------
    target_dir (str): Path to the directory where the ZIP file was decompressed.
                      If the file does not exist or an error occurs, returns None.
    """
    if not os.path.isfile(filename): # Check that the file exists
        raise FileNotFoundError(f"The file {filename} was not found")

    ############# ------------------- unzip
    # Path to a log file (to store errors)
    log_file = os.path.join(log_dir, "unzip_errors.log")

    out_dir = os.path.dirname(filename) or "."  # If there is no path in the filename, returns "."
    new_file = os.path.basename(filename)
    new_folder = os.path.splitext(new_file)[0]
    # Obtain the target directory
    target_dir = os.path.join(out_dir, new_folder)

    unzip_cmd = ["unzip",
                 "-q", # Quiet
                 "-o", # Overwrite files to avoid errors
                 filename,
                 "-d", target_dir] # Directory where it is unzipped
    print("Unzipping the downloaded genomes.\n")
    print(f"Executing: \n{shlex.join(unzip_cmd)}\n")
    try:
        subprocess.run(unzip_cmd, 
                       check = True, # Check for errors
                       capture_output = True, # Capture the output (stdout + stderr)
                       text = True)
    except subprocess.CalledProcessError as e: 
        print(f"Error executing unzip.\nDetails in {log_file}.")
        # Error stored into the log file
        with open(log_file, "a") as file:
            file.write("-------------- unzip log --------------\n")
            file.write(f"{e.stderr}\n")
        return None
    return target_dir


##------------------------ Function to download genomes
def download_genomes(source, name, # Arguments for the taxon_acc_flag(source, name) function
                     log_dir,
                     quality_filter_flag = False, # Flag to filter by quality metrics is needed
                     filename = "ncbi_dataset.zip", # Name of the file that will contain the genomes 
                                                  # Can have a path to specify a target directory
                     flags = None, # Flags for the datasets command
                     target_num = None,
                     seed = 1,
                     n_tries = 5
                     ): 
    """
    Executes the NCBI datasets command to download genomes (--dehydrated), allowing prior quality filtering.

    Args
    --------
    source (str): Source type for the query.
    name (str): Value or identifier associated with the source.
    log_dir (str): Directory where log files will be saved.
    quality_filter_flag (bool, optional): If evaluates to True, executes quality_filter_genomes() function.
                                          Default is False.
    filename (str, optional): Name of the resulting .zip file.  
                              Default is "ncbi_dataset.zip".
    flags (list, optional): List of additional flags to add to the 'datasets download' command. 
                            Default is None.
    target_num (int, optional): Maximum number of genomes to download.
                                Default is None (no limit).
    seed (int, optional): Seed for random sampling of genomes if target_num is specified.
                          Default is 1.
    n_tries (int, optional): Number of tries to try the download in case of error.
                             Default is 5.
    
    Return
    --------
    target_dir (str): Path to the directory where the ZIP file was decompressed.
                      If the file does not exist or an error occurs, returns None.
    """
    # Path to a log file (to store errors)
    log_file = os.path.join(log_dir, "datasets_download_errors.log")
    # Check that the directories exist to prevent errors
    out_dir = os.path.dirname(filename)
    if out_dir: # Prevents errors if the path has no folders 
        os.makedirs(out_dir, exist_ok=True) 
    if target_num:
        quality_filter_flag = True
    if quality_filter_flag:
        filter_result = quality_filter_genomes(source, name, filename, # Executes the filtering function
                                               target_num = target_num, 
                                               seed = seed, flags = flags) 
        if not filter_result: # If there is no accession, the datasets command cannot be used
            # The created directories are removed:
            if out_dir and os.path.exists(out_dir):
                if not os.listdir(out_dir):
                    # Removes empty directories (parent directories too) 
                    os.removedirs(out_dir)
            return
        source, name = filter_result # Unpack tuple

    # If no flags were added, empty list
    if flags is None: flags = []

    ########### -------------------- datasets download
    download_cmd = ["datasets", "download", "genome", 
                    *taxon_acc_flag(source, name), 
                    "--dehydrated", "--filename", 
                    filename, *flags]
    
    print(f"Executing: \n{shlex.join(download_cmd)}")
    for n in range(1, n_tries + 1):
        try:
            subprocess.run(download_cmd, 
                           check = True, capture_output = True, text = True) 
            break
        # Capture the exception
        except subprocess.CalledProcessError as e:
            print(f"Error executing datasets download on try {n}/{n_tries}.\nDetails in {log_file}.")
            # Error stored into the log file
            with open(log_file, "a") as file:
                file.write("-------------- datasets download log --------------")
                file.write(f"try {n}/{n_tries}")
                file.write("----------------------------------------------------")
                file.write(f"{e.stderr}\n")
                
            # The created directories are removed (provided they are completely empty)
            if n == n_tries:
                if out_dir and os.path.exists(out_dir):
                    if not os.listdir(out_dir): 
                        # Removes empty directories (parent directories too)
                        os.removedirs(out_dir) 
                return None
        
    target_dir = unzip(log_dir, filename) # Executes unzip function to decompress the ZIP file
    if not target_dir: # If the unzip function fails, return None
        return None
    
    return target_dir