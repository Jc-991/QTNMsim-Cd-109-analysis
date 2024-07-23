#!/usr/bin/env python
# coding: utf-8

# # Raw Data Analysis
# This script is used to find the histogram for the selected partcles generated with their corresonding values. The particles are classified according to their Monte Carlo Particle number (PDG). The PDG for electron, electron neutrino and photon (gamma-rays) are 11, 12 and 22 respectively. 

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import argparse


# In[10]:


def rawhist(filename, column, eval=None, gval=None):
    '''
    This function reads the csv file and selects the targeted column for electrons, neutrinos, and gamma-rays, 
    and makes it into a histogram.
    
    filename - the name of the csv file
    column - name of the selected column
    eval - a list of x-coordinates where vertical lines should be added for electrons (default is None)
    gval - a list of x-coordinates where vertical lines should be added for gamma-rays (default is None)
    '''
    df = pd.read_csv(filename)

    # Ensure column names are stripped of leading/trailing whitespace
    df.columns = df.columns.str.strip()
    
    # Define a dictionary to handle different PDG codes and their respective titles and labels
    pdg_info = {
        11: {"name": "electrons", "color": "red", "lines": eval, "line_label": "Expected peaks of electrons"},
        12: {"name": "electon neutrinos", "color": "green", "lines": eval, "line_label": "Expected peaks of neutrinos"},
        22: {"name": "gamma-rays", "color": "blue", "lines": gval, "line_label": "Expected peaks of gamma rays"}
    }

    for pdg_code, info in pdg_info.items():
        # Filter the DataFrame for rows where 'PDG' column is equal to the current pdg_code
        filtered_df = df[df['PDG'] == pdg_code]
        
        # Check if the filtered DataFrame is not empty
        if not filtered_df.empty:
            # Select the targeted column from the filtered DataFrame
            Tcolumn = filtered_df[column.strip()]

            plt.figure(figsize=(8, 6))
            plt.hist(Tcolumn, bins=1000, edgecolor='black')
            plt.yscale('log')  # make the y scale logarithmic
            plt.title(f'Distribution of {column.strip()} for {info["name"]}')
            plt.xlabel(column.strip())
            plt.ylabel('Frequency')
            plt.grid(True)

            # Flag to ensure the label is added only once
            label_added = False

            # Add vertical lines if specified
            if info["lines"]:
                for x_value in info["lines"]:
                    if not label_added:
                        plt.axvline(x=x_value, color=info["color"], linestyle='--', linewidth=1.5, label=info["line_label"])
                        label_added = True
                    else:
                        plt.axvline(x=x_value, color=info["color"], linestyle='--', linewidth=1.5)

            plt.legend()

            # Save the plot with an appropriate filename
            plt.savefig(f'{column.strip()}_plot_{info["name"]}.png')
            plt.show()
        else:
            print(f"No data available for PDG code {pdg_code} in column {column.strip()}")



# In[ ]:




