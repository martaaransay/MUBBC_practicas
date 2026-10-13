# Import libraries 
import json
import os


# --------Definition of functions

#------------------------ Function to convert a jsonl file into a dictionary
def load_assembly_reports(report_path):
    """
    Loads an NCBI assembly_data_report.jsonl file into a dictionary.

    Args
    ------
    report_path (str): The path to the assembly_data_report.jsonl file.

    Return
    ------
    reports (dict): Dictionary {genome accession : record (dict)}.
    """
    reports = {} # Dictionary to store the records
    with open(report_path, "r") as f:
        for line in f: # Iterates over each line in the file
            line = line.strip() # Removes leading whitespace
            if line:
                record = json.loads(line) # Converts the line from JSON to a dictionary
                # Stores the record in the reports dictionary with the genome accession as the key
                reports[record["accession"]] = record
    return reports

#------------------------ Function to extract values from a JSON record
def get_nested(data, path):
    """
    Gets a value from a nested JSON record using a dotted path.

    Args
    ------
    data (dict): JSON record.
    path (str): Dotted path to the value.

    Return
    ------
    value (str): The value found. 
                 If several unique values are found, they are joined with ";". 
                 "NA" if the path does not exist.
    """
    values = [data] # List to store the current level of data
    for key in path.split("."): # Iterates over each key in the dotted path
        next_values = [] # List to store the next level of data
        for value in values: # Iterates over each value in the current level
            if isinstance(value, dict): # If the value is a dictionary
                if key in value: # If the key exists in the dictionary, add its value to next_values
                    result = value[key]
                    if isinstance(result, list): # If the results is a list, extend next_values with the list
                        next_values.extend(result)
                    else: # If the result is not a list, append it to next_values
                        next_values.append(result)
        
        values = next_values # Update values to the next level of data
    # If no values are found, return "NA"
    if not values: 
        return "NA"
    
    unique_values = []  # List to store the unique values
    for v in values:
        v_str = str(v)  # Value to string
        if v_str not in unique_values:  # Skip duplicates
            unique_values.append(v_str)
    final_results = ";".join(unique_values)  # Join the unique values with ";"
    
    return final_results
    

#------------------------ Function to add species and extra columns to the hits dataframe
def add_genome_metadata(hits_df, extra_columns = {}):
    """
    Adds the species and, optionally, extra columns taken from the
    assembly_data_report.jsonl of the genome each hit belongs to.

    Args
    ------
    hits_df (pd.DataFrame): Hits table with the columns genome and genome_dir.
    extra_columns (dict, optional): Extra columns to add, as
                                    {name: dotted path in the JSON}.
                                    Defaults to {}.

    Return
    ------
    hits_df (pd.DataFrame): A copy of the input with the species column
                            and one column per entry of extra_columns 
                            ("NA" if not found).
    """
    # Copy of the dataframe to not modify the original
    hits_df = hits_df.copy()
    
    # Store the assembly_data_report.jsonl of each genome_dir in a dictionary to avoid reading it multiple times
    reports_jsonl = {} 

    species = [] # List to store the species of each hit
    extras = {} # Dictionary to store the extra columns of each hit
    for col in extra_columns: 
        extras[col] = [] # Initialize the extra columns with empty lists
         
    # Iterate over each hit (gets the path of the genome and its accession)
    ##### General structure --> {species name}/{species name}/ncbi_dataset/data/{genome accession}
    for genome_dir, genome in zip(hits_df["genome_dir"], hits_df["genome"]): # zip to iterate over two columns simultaneously
        data_dir = os.path.dirname(genome_dir) # {species name}/{species name}/ncbi_dataset/data
        ncbi_dataset_dir = os.path.dirname(data_dir) # {species name}/{species name}/ncbi_dataset
        species_dir = os.path.dirname(ncbi_dataset_dir)
        # Append to the species list the name of the current species
        species.append(os.path.basename(species_dir)) # {species name}

        # If extra columns are requested, load the assembly_data_report.jsonl of the current species 
        if extra_columns: 
            if data_dir not in reports_jsonl: # Loaded only once per species
                jsonl_file = os.path.join(data_dir, "assembly_data_report.jsonl")
                # Stores the report in the dictionary with the data_dir as the key
                reports_jsonl[data_dir] = load_assembly_reports(jsonl_file)
            # For each hit, get the record of the genome and extract the values of the extra columns
            record = reports_jsonl[data_dir].get(genome, {}) # If not found, returns an empty dictionary
            
            for col, path in extra_columns.items(): # Iterates over each extra column and its dotted path
                # Addds to the dictionary the value of the extra column for the current hit
                extras[col].append(get_nested(record, path))
    
    # Adds the species and extra columns to the hits dataframe
    hits_df["species_folder"] = species
    for col, values in extras.items():
        hits_df[col] = values
    return hits_df