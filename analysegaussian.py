#!/usr/bin/env python
# coding: utf-8

# # Gaussian Distribution and Data Analysis
# This script should mainly be used for the find gaussian distribution for the raw and data analysis for the targeted peaks hence analyse the precision of the calibration.

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks, peak_widths
import scipy.stats as stats
from scipy.optimize import fsolve
import argparse


# ## 1. Kernel Density Estimtion
# The code in this section is mainly to use the kernel density estimation method to find the smooth curve for the non-uniform data. Please note that the fit parameter can change the fit level for PDG so please find appropriate parameter for the case.
# 

# In[1]:


def ehist_kde(filename, column, m):
    '''
    This function reads the csv file and plot the corresponding Kernal density estimation distributions
    filename - name of the file
    column - the targeted column in csv file 
    m - fit parameter
    '''
    df = pd.read_csv(filename)
    df.columns = df.columns.str.strip()
    filtered_df = df[df['PDG'] == 11]

    if not filtered_df.empty:
        Tcolumn = filtered_df[column.strip()]

        plt.figure(figsize=(8, 6))
        
        # Plot histogram
        plt.hist(Tcolumn, bins=300, density=True, edgecolor='black', alpha=0.6, label='Data histogram')

        # KDE estimation
        kde = gaussian_kde(Tcolumn)
        x = np.linspace(Tcolumn.min(), Tcolumn.max(), 1000)
        kde.set_bandwidth(bw_method=kde.factor / 3)  # change the fitness level of the kde estimation
        kde_values = kde(x)
        
        # Plot KDE
        plt.plot(x, kde_values, 'r-', linewidth=2, label='KDE fit')

        plt.title(f'Distribution of {column.strip()} for electron (t = 1000 um)')
        plt.xlabel(column.strip())
        plt.ylabel('Density')
        plt.legend()
        
        plt.grid(True)
        #plt.savefig(f'{column.strip()}shell_plot_r=50.png')
        plt.savefig(f'pointr_1.png')
        plt.show()
    else:
        print(f"No data available for PDG code 11 in column {column.strip()}")


# ## 2. Peaks and Width of the KDE fit
# 
# In this section, it is mainly to find the each peak and analyse the width for each one. You might have to input manually the range of your tergeted peaks. Here the relative height ratio r is the ratio of the height used to find the width for each peak.

# In[50]:


def ehist_kde_peaks(filename, column, m, r, min, max):
    '''
    This function plots the histogram, KDE, and finds the width of peaks at 0.75 height.
    m - the fitness level of the KDE curve and data points.
    r - ratio of relative height
    min - left boundary of the targeted range
    max - right boundary of the targeted range
    '''
    # Read the CSV file
    df = pd.read_csv(filename)
    df.columns = df.columns.str.strip()
    
    # Filter the dataframe for PDG code 11
    filtered_df = df[df['PDG'] == 11]
    print(f'the ratio of relative height is {r}')
    print(f'the fit parameter is set to {m}')
    print()
    
    if not filtered_df.empty:
        Tcolumn = filtered_df[column.strip()]

        plt.figure(figsize=(8, 6))
        
        # Plot histogram
        plt.hist(Tcolumn, bins=300, density=True, edgecolor='black', alpha=0.6, label='Data histogram')

        # KDE estimation
        kde = gaussian_kde(Tcolumn)
        n = 2000
        x = np.linspace(Tcolumn.min(), Tcolumn.max(), n)
        kde.set_bandwidth(bw_method=kde.factor / m)  # change the fitness level of the kde estimation
        kde_values = kde(x)

        # Plot KDE
        plt.plot(x, kde_values, 'r-', linewidth=2, label='KDE')

        # Find peaks in KDE
        peaks, _ = find_peaks(kde_values)
        peak_x = x[peaks]
        peak_y = kde_values[peaks]

        # Find intersections and calculate widths
        for i, px in enumerate(peak_x):
            if min <= px <= max:
                rel_h = r * peak_y[i]
                # Define function for fsolve
                def func(x1):
                    return kde(x1) - rel_h
            
                # Find roots using fsolve with initial guesses
                try:
                    root_left = fsolve(func, px - 1)[0]
                    root_right = fsolve(func, px + 1)[0]
                
                    # Calculate width
                    width = root_right - root_left

                    # Print width
                    print(f'Peak {i+1}: the width is {width:.6f}')

                    # Plot intersection points and relative height line
                    plt.plot(root_left, rel_h, "yo")
                    plt.plot(root_right, rel_h, "go")
                    plt.hlines(rel_h, root_left, root_right, colors='green', linestyles='dashed')
                    plt.axvline(x=min, color='green', linestyle='--', linewidth=1.5)
                    plt.axvline(x=max, color='green', linestyle='--', linewidth=1.5)

                except RuntimeError:
                    print(f"Failed to find intersection for peak {i+1}. Skipping...")
                    

        # Plot peaks
        plt.plot(peak_x, peak_y, "bo", label='Peaks')

        plt.title(f'Distribution of {column.strip()} ')
        plt.xlabel(column.strip())
        plt.ylabel('Density')
        plt.legend()
        plt.grid(True)
        #plt.xlim(15,21)
        #plt.ylim(0, 0.01)
        #plt.savefig('example_kde_width')
        plt.show()
    else:
        print(f"No data available for PDG code 11 in column {column.strip()}")


# ## 3. Analysis for the width and thickness
# This function is used for analysis for the result from previous section. The outcome of this function includes: yields per 100k events, positions of the peaks, the standard error for the peak position and standard error. The method to calculate the standard error is to calculate the standard deviation using equation: 
# $$\sigma = \sqrt{\frac{1}{N} \sum_{i=1}^{N} (x_i - \mu)^2}  $$
# Where N is number of yields and $\mu$ is the mean for the data. The standard errors for the peak positions and width are calculated according to following functions:
# $$ \text{SEM}_{\text{peak}} = \frac{\sigma}{\sqrt{N}} $$
# $$ \text{SEM}_{\text{width}} = \frac{\sigma}{\sqrt{2N}} $$
# Please not that in this function, you have to input manually your target range for your peaks.

# In[46]:


def ehist_kde_error(filename, column, m, min, max):
    '''
    This function plots the histogram, KDE, and finds the width of peaks at 0.75 height.
    m - the fitness level of the KDE curve and data points.
    min - left boundary of the targeted range
    max - right boundary of the targeted range
    '''
    # Read the CSV file
    
    df = pd.read_csv(filename)
    df.columns = df.columns.str.strip()
    
    # Filter the dataframe for PDG code 11
    filtered_df = df[df['PDG'] == 11]
    print(f'The fit parameter is set to {m}')
    print()
    
    if not filtered_df.empty:
        Tcolumn = filtered_df[column.strip()]
        frequency, T_edge = np.histogram(Tcolumn, bins = 2000)
        
        # KDE estimation
        kde = gaussian_kde(Tcolumn)
        x = np.linspace(Tcolumn.min(), Tcolumn.max(), 2000)
        kde.set_bandwidth(bw_method=kde.factor / m)  # change the fitness level of the kde estimation
        kde_values = kde(x)

       ##########################
        #calculate the uncertainty of the fit
        
        #Find the mean in the range of (15-19)
        kdetot = []
        for xi in range(15,19):
            kdei = kde(xi)
            kdetot.append(kdei)
        

        #total number of events in the targeted range 
        n = [ke for ke in Tcolumn if min <= ke < max]
        mean = np.mean(np.array(n))
        print(f'the total number of events is {len(n)}')
        

        #standard deviation 
        if n:
            std = np.sqrt(np.mean((np.array(n)-mean)**2))
        print(f'std {std:.6f}')
        #standard error of the mean (SEM)
        sempeak = std / np.sqrt(len(n))
        semwid = std / np.sqrt(2*len(n))

        #find the peaks
        peaks, _ = find_peaks(kde_values)
        peak_x = x[peaks]
        peak_y = kde_values[peaks]
        #plt.plot(x, kde_values, 'r-', linewidth=2, label='KDE')
        #plt.plot(peak_x, peak_y, "bo", label='Peaks')
        
        print()
        for i, px in enumerate(peak_x):
            if min <= px <= max:
                print(f'the x position of the peak is {px:.6f}')
                print(f'the y position of the peak is {peak_y[i]:.6f}')
        print()
        print(f'the standard error for the peak position is {sempeak:.6f}')
        print(f'the standard error for the width is {semwid:.6f}')


# In[ ]:




