# Import libraries 
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn.layout.venn3 import DefaultLayoutAlgorithm
from matplotlib.patches import Patch

# Plotting

#------------------------ Function to plot e-values distribution
def plot_evalue(logdata, labels, 
                title = "E-value Distribution", outdir = "results/plots", 
                filename = "evalue_distribution.png"):
    """
    Saves a boxplot of the -log10 of e-values for each query.
    
    Args
    ------
    logata (list): A list of lists containing the -log10 of e-values for each query.
    labels (list): A list of the names of the queries.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(7, 5)) # Creates axes for the plot
    # Creates the boxplot with the data and labels defined
    ax.boxplot(logdata, tick_labels=labels) 
    # Adjusting labels / title:
    ax.set_ylabel("-log10(E-value)") # Y label
    ax.set_title(title) # Title of the plot
    plt.xticks(rotation=20) # Rotates the x-axis labels
    plt.tight_layout() # Adjusts the layout to avoid overlapping
    # Saves the plot
    plt.savefig(os.path.join(outdir, filename)) 
    plt.close(fig)
    return

#------------------------ Function to plot bitscore distribution
def plot_bitscore(bitscore_data, labels, 
                  title = "BitScore Distribution", outdir = "results/plots", 
                  filename = "bitscore_distribution.png"):
    """
    Saves a boxplot of the bitscores for each query.
        
    Args
    ------
    bitscore_data (list): A list of lists containing the bitscores for each query.
    labels (list): A list of the names of the queries.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    
    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5)) # Creates axes for the plot
    # Creates the boxplot with the data and labels defined
    ax.boxplot(bitscore_data, tick_labels=labels)
    # Adjusting labels / title:
    ax.set_ylabel("BitScore") # Y label
    ax.set_title(title) # Title of the plot
    plt.xticks(rotation=20) # Rotates the x-axis labels
    plt.tight_layout() # Adjusts the layout to avoid overlapping
    # Saves the plot
    plt.savefig(os.path.join(outdir, filename))
    plt.close(fig)
    return

#------------------------ Function to plot e-value thresholds
def plot_evalue_thresholds(dfs, times,
                           evalues = np.logspace(1, -300, 200),
                           colors = ["#963FB0", "#41B883", "#1A6D70"],
                           title = "E-Values thresholds", outdir = "results/plots",
                           filename = "evalue_thresholds.png"):
    """
    Saves a plot of the number of hits and unique genomes for each query
    at different e-value thresholds.

    Args
    -----
    dfs (dict): query name -> blastn dataframe.
    times (dict): query name -> execution time.
    evalues (list, optional): List of e-value thresholds to consider.
                              Defaults to np.logspace(1, -300, 200).
    colors (list, optional): List of colors for each query.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    """
    os.makedirs(outdir, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # Iterate over each query and its corresponding dataframe and color
    for (query, df), color in zip(dfs.items(), colors): 
        # Lists to store the number of hits and unique genomes for each e-value threshold
        n_hits = []
        n_unique = []

        for evalue in evalues: # Iterates over thresholds
            evalue_mask = df["Evalue"] <= evalue # Hits that pass the threshold
            n_hits.append(evalue_mask.sum()) # Counts the number of hits 
            n_unique.append(df.loc[evalue_mask, "seq_id"].nunique()) # Counts the number of seq_id 

        label = f"{query} - {times[query]:.2f} secs"
        axes[0].plot(evalues, n_hits, label = label, color = color)
        axes[1].plot(evalues, n_unique, label = label, color = color)

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("E-value threshold")
        ax.legend()
    axes[0].set_ylabel("Total hits")
    axes[1].set_ylabel("Unique hits (seq_id)")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, filename))
    plt.close(fig)

#------------------------ Function to plot the comparison between methods
def plot_hits_comparison(summary, total_times,
                         colors = ["#35C40E", "#F791F4", "#F1BE3B"],
                         title = "Method comparison", outdir = "results/plots",
                         filename = "method_comparison.png"):
    """
    Saves a plot comparing the total hits, unique hits, and execution times across three methods.

    Args
    -----
    summary (dict): A dictionary where keys are the names of the methods and
                    values are tuples with the total hits and unique hits.
    total_times (dict): A dictionary where keys are the names of the methods and
                        values are the execution times.
    colors (list, optional): List of colors for each method.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok=True)

    methods = list(summary.keys()) # List of methods
    n = len(methods) # Number of methods
    # Lists to store results 
    totals = []
    uniques = []
    times = []
    # Colors to plot
    dict_colors = {}
    
    for i, meth in enumerate(methods): # Store data from each method
        dict_colors[meth] = colors[i]# Set a color for each method
        totals.append(summary[meth][0])
        uniques.append(summary[meth][1])
        times.append(total_times[meth])
        

    fig, axes = plt.subplots(1, 2, figsize=(13,5)) # Creates axes for the plot
    x = np.arange(2) # Two comparisons: total and unique
    width = 0.8 / n # Width of each bar

    for i, m in enumerate(methods):
        values = [totals[i], uniques[i]] # Pairs to plot for each method
        center = (n - 1) / 2 # Middle bar 
        offset = width * (i - center) # Horizontal shift
        # Plot the bars
        axes[0].bar(x + offset, # Shift
                    values, width = width, 
                    label = m, color = dict_colors[m])
    # Labels, titles and axes
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(["Total hits", "Unique hits"])
    axes[0].set_ylabel("Number of hits")
    axes[0].set_title("Total hits vs. Unique hits")
    axes[0].legend(title = "Method")

    # Subplot 2: "Time of execution"
    # Background plot to plot times
    axes[1].plot(range(n), times, "-", color = "gray")
    for i, m in enumerate(methods):
        # Plot time of each method
        axes[1].plot(i, times[i], "o", color = dict_colors[m])
        # Annotate the time of each method
        axes[1].annotate(f"{times[i]:.2f}s", (i, times[i]), 
                         textcoords = "offset points", xytext = (0, 10),
                          ha="center")
    
    # Labels, titles and axes
    axes[1].set_xticks(range(n))
    axes[1].set_xticklabels(methods)
    axes[1].set_ylabel("Time (secs)")
    axes[1].set_title("Time of execution")
    axes[1].margins(y=0.15)  # deja hueco arriba para las anotaciones

    plt.suptitle(title) # Set the principal title

    plt.tight_layout()# Adjusts the layout to avoid overlapping
    # Saves the plot
    plt.savefig(os.path.join(outdir, filename))
    plt.close(fig)

    return
    

#------------------------ Function to analyze the overlapping
def plot_hits_overlap(sets_dict, 
                      colors = ["#35C40E", "#F791F4", "#F1BE3B"], 
                      title = "Hits overlapping", outdir = "results/plots",
                      filename = "hits_overlap.png"):
    """
    Saves a plot showing the overlapping of hits across three methods using a Venn diagram.

    Args
    -----
    sets_dict (dict): A sets dictionary where keys are the names of the methods 
                      and values are sets of hits.
    colors (list, optional): List of colors for each method.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    """
    os.makedirs(outdir, exist_ok=True)
    
    names = list(sets_dict.keys()) # Sets of names of methods
    sets = list(sets_dict.values()) # Sets of hits for each method


    fig, ax = plt.subplots(figsize=(8, 8)) # Creates axes for the plot
    # Creates the Venn diagram with the sets and their corresponding labels and colors
    venn3(sets, set_labels=("", "", ""), # Labels for each set
          set_colors = colors, # Colors for each set
          alpha = 0.5, ax = ax, # Transparency and axes
          # Plot with the regions with the same size
          layout_algorithm = DefaultLayoutAlgorithm(fixed_subset_sizes = (1,) *7)) # 7 regions
    
    handles = [] # List to store the legend handles
    for color, meth, set in zip(colors, names, sets): # Iterate over each method and its corresponding color
        # Create a patch for the legend with the corresponding color and label
        handle = Patch(facecolor = color, alpha = 0.5, 
                       label = f"{meth} (n={len(set)})")
        handles.append(handle)

    # Add the legend to the plot with the created handles
    ax.legend(handles = handles, loc = "upper right")

    ax.set_title(title) # Set the title of the plot
    plt.tight_layout() # Adjusts the layout to avoid overlapping
    plt.savefig(os.path.join(outdir, filename))
    plt.close(fig)
    return

#------------------------ Function to plot the distribution of Pc variants by species
def plot_pc_by_species(csv_path,
                       title = "Distribution of Pc variants by Species",
                       outdir = None, filename = "pc_variants_by_species.png"):
    """
    Plots a stacked bar chart with the number of hits of each Pc variant,
    split by species.

    Args
    ------
    csv_path (str): The path to the final.csv file (needs the columns Pc_variant and species).
    title (str): Title of the plot.
    outdir (str, optional): Directory where the plot is saved. 
                            Default to None (the plot is shown instead of saved).
    filename (str): Name of the output file. 
                    Default to "pc_variants_by_species.png".

    Return
    ------
    count_table (pd.DataFrame): Number of hits per Pc variant (rows) and species (columns).
    """
    
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip() # Removes spaces from the column names

    # Cross table: number of occurrences of each Pc variant in each species
    count_table = pd.crosstab(df["Pc_variant"], df["species"])

    fig, ax = plt.subplots(figsize = (12, 7))
    # Colors to plot
    colors = ["#984FB3", "#AAAA00", "#FF5733", "#FF8D1A", 
              "#FFE400", "#FF96E8", "#96FFE5", "#4F80B3"]
    # Stacked bars (one color per species)
    count_table.plot(kind = "bar", stacked = True, ax = ax,
                     color = colors[:len(count_table.columns)], 
                     edgecolor = "black")

    ax.set_xlabel("Pc variant")
    ax.set_ylabel("Count")
    ax.set_title(title)

    plt.setp(ax.get_xticklabels(), rotation = 45, ha = "right")
    # Legend outside the plot to not cover the bars
    ax.legend(title = "Species", bbox_to_anchor = (1.05, 1), loc = "upper left")
    fig.tight_layout()

    if outdir: # Saves the plot
        os.makedirs(outdir, exist_ok = True)
        fig.savefig(os.path.join(outdir, filename), bbox_inches = "tight")
        plt.close(fig)
    else: # Shows the plot
        plt.show()

    return count_table