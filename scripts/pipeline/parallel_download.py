# Import libraries 
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed 
from tqdm import tqdm


# --------Definition of functions

#------------------------ Function to create a recovery file to allow the modification of the fetch file
def create_recovery_file(fetch_directory):
    """
    Creates or identifies the recovery file to track download progress.
    
    Args
    ------
    fetch_directory (str): The path to the directory containing the fetch file.
        
    Return
    ------
    recovery_file (str): The path to 'recovery.txt' if successful
                         None if it failed. 
    """
    # Path to the fetch file
    fetch_file = os.path.join(fetch_directory, "fetch.txt")
    # Path to the recovery file
    recovery_file = os.path.join(fetch_directory, "recovery.txt")
    
    # If the recovery file already exists, it does nothing (control)
    if os.path.isfile(recovery_file): pass
    # If the fetch file already exists, renames it as the recovery file
    elif os.path.isfile(fetch_file):
        os.replace(fetch_file, recovery_file)  # Overwrites the file
    else: # Control to avoid errors if the download failed
        print("Neither fetch file nor recovery file")
        return None
    
    return recovery_file

#------------------------ Function to split the fetch file into parts
def splitting_fetch_file(recovery_file, n_split = 100):
    """
    Reads the recovery file and splits it.
    
    Args
    ------
    recovery_file (str): The path to the recovery text file.
    n_split (int): The number of lines to extract for the current batch.
                   Default to 100.
    
    Return
    ------
    splits (list): The first n_split lines to be processed now.
    recovery (list): The remaining lines to be processed later.
    """
    splits = [] # List to store the first n_split lines to be processed
    recovery = [] # List to store the remaining lines
    # Open the recovery file in read mode
    with open(recovery_file, "r") as f:
        for line in f:
            splits.append(line) # Stores the first n_split lines
            if len(splits) == n_split: 
                recovery = f.readlines() # Store the rest of the lines
                break

    return splits, recovery

#------------------------ Function to execute datasets rehydrate of a splitted fetch file
def processing(directory, log_dir, n_split = 100, 
               max_workers = 18, n_tries = 5,
               position = 0, lab = None):
    """
    Processes the genome downloads in chunks.

    Args
    ------
    directory (str): The target root directory for the genomes.
    log_dir (str): Directory where logs will be saved.
    n_split (int): The number of genomes to download per batch.
    max_workers (int): Number of concurrent connections for the datasets tool.
    n_tries (int): Maximum number of execution attempts per chunk.
    position (int, opt): Row of the progress bar (needed when several run in parallel).
                         Default to 0.
    lab (str, opt): Label of the progress bar.
                    Default to None.
        
    Return
    ------
    bool: True if the entire directory was processed successfully.
          False if it failed.
    """
    print("-------------------- Rehydrating genomes --------------------")

    # Directory where is the fetch file
    fetch_directory = os.path.join(directory, "ncbi_dataset")
    # Creates the recovery file
    recovery_file = create_recovery_file(fetch_directory)
    if recovery_file is None: # Control to avoid errors
        return False
    
    # Path to the new fetch file (with the current lines to rehydrate)
    current_fetch_file = os.path.join(fetch_directory, "fetch.txt")
    # Path to a failed file (with the lines that have falied)
    failed_file = os.path.join(fetch_directory, "failed.txt")
    # Path to a log file (to store errors)
    log_file = os.path.join(log_dir, "rehydrate_errors.log")
    
    
    # Command to execute datasets rehydrate
    rehydrate_cmd = ["datasets", "rehydrate", 
                     "--directory", directory, # Directory to store genomes
                     "--max-workers", str(max_workers)] # Number of workers

    # Total number of lines to process (for the progress bar)
    with open(recovery_file, "r") as file:
        total = 0
        for _ in file:
            total += 1

    n_failed = 0
    if not lab: # If no lab provided, name of the directory
        lab = os.path.basename(directory)
    
    pbar = tqdm(total = total, # Total number of files to process
                desc = lab, # Label of the bar
                position = position, # Row where the bar is drawn to avoid overlapping between threads
                unit = "files",  # Unit shown in the speed (files/s)
                leave = True) # Keeps the bar on screen when it finishes

    while True:
        # Split the fetch file
        splitted, recovery = splitting_fetch_file(recovery_file, n_split)
        if not splitted: # If no more lines in the fetch file, stops rehydrating
            break
        # Writes the current lines into the current fetch file
        with open(current_fetch_file, "w") as f:
            f.writelines(splitted)
        
        success = False # Flag to check succession
        # Control to avoid rehydrate errors: In case of error, tries again the rehydratation
        for n in range(1, n_tries + 1):

            try:
                subprocess.run(rehydrate_cmd, 
                               check=True, capture_output=True, text=True)
                success = True
                break
            except subprocess.CalledProcessError as e:
                tqdm.write(f"Error executing datasets rehydrate on a concrete split on try {n}/{n_tries}.\
                             Details in {log_file}.")
                # Error stored into the log file
                with open(log_file, "a") as file:
                    file.write("-------------- datasets rehydrate log --------------\n")
                    file.write(f"try {n}/{n_tries}\n")
                    file.write("----------------------------------------------------\n")
                    file.write(f"{e.stderr}\n")
                
        if not success: 
            tqdm.write(f"Chunk failed after {n_tries} tries.")
            n_failed += len(splitted)
            pbar.set_postfix(failed = n_failed) # Shows failed files in the bar
            with open(failed_file, "a") as f:
                f.writelines(splitted)
        
        # Writes the rest of the lines into the recovery file
        with open(recovery_file, "w") as f:
            f.writelines(recovery)
        pbar.update(len(splitted)) # Advance the bar

    pbar.close()
    # If there were failed splits, writes them into the failed file to continue the execution
    if os.path.isfile(failed_file):
        print(f">>> WARNING: failed splits in '{directory}'.")
        os.replace(failed_file, recovery_file)
        return False
    return True 


#------------------------ Function to process several datasets rehydrate 
def process_all(dirs, log_dir, total_workers = 18, 
                n_split = 100, n_tries = 5):
    """
    Processes multiple species directories concurrently using a thread pool.
    
    Args
    ------
    dirs (dict): A dictionary where keys are species names and values are directory paths.
    log_dir (str): Directory where logs will be saved.
    total_workers (int): Total number of workers available to distribute across all species.
    n_split (int): The chunk size or number of splits to pass to the processing function.
    n_tries (int): The maximum number of retry attempts for processing.
    
    Return
    ------
    failed (list): A list containing the names of the species that failed processing.
    """
    # Distribute workers among directories to process    
    n_species = len(dirs) 
    per_dir = max(1, 
                  total_workers // n_species) 

    failed = []  # Failed species

    # Create a thread pool (one thread per species)
    with ThreadPoolExecutor(max_workers = n_species) as executor:

        # Queue one task per specie
        tasks = {}  # Dictionary to store the tasks for each specie
        for idx, (specie, directory) in enumerate(dirs.items()):

            # Submit executes the function using a free thread
            fut_task = executor.submit(processing, directory, log_dir, n_split, per_dir, n_tries,
                                       idx, specie)
            # fut_task is a "Future"
            # Stores the task and its specie
            tasks[fut_task] = specie

        # Store results when tasks are completed 
        for task in as_completed(tasks):
            name_species = tasks[task]      
            state_processing = task.result() # Result of the execution

            if not state_processing: # state_processing is bool
                failed.append(name_species) # If False, exeuction has failed

    # When the with block finishes, ensures that all threads have finished
    return failed
