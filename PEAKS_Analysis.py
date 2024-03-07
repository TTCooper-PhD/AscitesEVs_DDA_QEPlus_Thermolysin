import csv
import json
import math
import matplotlib.pyplot as plt 
from matplotlib_venn import venn2,venn3
import os
import pandas as pd
import seaborn as sns
import numpy as np
import statsmodels,statistics

# Thermolysin Analysis (PEAKS):

# Load the CSV file
thermolysin_df = pd.read_csv('PEAKS\\Thermolysin\\peaks_ptm.protein-peptides.csv')
tryplysc_df = pd.read_csv('PEAKS\\TrypLysC\\peaks_ptm.protein-peptides.csv')


#

def process_dataframe(df, enzyme):
    print(f"Processing dataframe for {enzyme}")
    print(f"Columns: {df.columns}")
    # Split the data by Donor and Fraction
    donors = []
    fractions = []
    for col in df.columns:
        if 'Area' in col:
            parts = col.split('_')
            if len(parts) >= 2:
                donors.append(parts[0].split(' ')[1])  # Splitting on space and underscore to get the donor number
                fractions.append(parts[1])

    # Initialize results dictionary
    results = {}
    # Process each donor and fraction
    for donor in set(donors):
        for fraction in set(fractions):
            column_name = f'Area {donor}_{fraction}'
            if column_name in df.columns:
                # Filter data for current donor and fraction
                filtered_df = df[df[column_name].notna() & (df[column_name] > 0)]
                # Calculate metrics
                peptide_ids = filtered_df['Peptide'].nunique()
                unique_peptides = filtered_df[filtered_df['Unique'] == True]['Peptide'].nunique()
                avg_charge_state = filtered_df['z'].mean()
                avg_peptide_length = filtered_df['Length'].mean()
                avg_mass = filtered_df['Mass'].mean()
                avg_m_z = filtered_df['m/z'].mean()
                # Store results
                results[(donor, fraction)] = {
                    'Number of Peptide IDs': peptide_ids,
                    'Number of Unique Peptides': unique_peptides,
                    'Average Charge State': avg_charge_state,
                    'Average Peptide Length': avg_peptide_length,
                    'Average Mass': avg_mass,
                    'Average m/z': avg_m_z,
                }
                print(f"Processed {enzyme}, Donor: {donor}, Fraction: {fraction}")
    return results


thermolysin_results = process_dataframe(thermolysin_df, "Thermolysin")
tryplysc_results = process_dataframe(tryplysc_df, "TrypLysC")

def plot_results(results1, enzyme1, color1, results2, enzyme2, color2):
    fractions = ['Frac1', 'Frac2', 'Frac3', 'Frac4', 'Frac5', 'Frac6']
    donors = ['Donor1', 'Donor2', 'Donor3']

    # For Peptide Count
    fig, ax = plt.subplots()
    for results, enzyme, color in [(results1, enzyme1, color1), (results2, enzyme2, color2)]:
        total_peptide_counts = []
        unique_peptide_counts = []
        for frac in fractions:
            peptide_counts_per_donor = [results.get((donor, frac), {'Number of Peptide IDs': 0})['Number of Peptide IDs'] for donor in donors]
            unique_peptide_counts_per_donor = [results.get((donor, frac), {'Number of Unique Peptides': 0})['Number of Unique Peptides'] for donor in donors]

            total_peptide_counts.append(np.mean(peptide_counts_per_donor))
            unique_peptide_counts.append(np.mean(unique_peptide_counts_per_donor))

        ax.plot(fractions, total_peptide_counts, color=color, linestyle='-', marker='o', label=f'Total Peptides {enzyme}')
        ax.plot(fractions, unique_peptide_counts, color=color, linestyle='--', marker='o', label=f'Unique Peptides {enzyme}')
    ax.set_xlabel('Fraction')
    ax.set_ylabel('Peptide Count')
    ax.legend()
    plt.title('Peptide Analysis')
    plt.savefig('Peptide_Analysis.png', dpi=300)
    plt.show()

    # For Average Charge State
    fig, ax = plt.subplots()
    for results, enzyme, color in [(results1, enzyme1, color1), (results2, enzyme2, color2)]:
        avg_charge_states = []
        std_dev_charge_states = []
        for frac in fractions:
            charge_states_per_donor = [results.get((donor, frac), {'Average Charge State': 0})['Average Charge State'] for donor in donors]

            avg_charge_states.append(np.mean(charge_states_per_donor))
            std_dev_charge_states.append(np.std(charge_states_per_donor))

        ax.errorbar(fractions, avg_charge_states, yerr=std_dev_charge_states, color=color, linestyle='-', marker='o', label=f'Average Charge State {enzyme}')
    ax.set_xlabel('Fraction')
    ax.set_ylabel('Average Charge State')
    ax.legend()
    plt.title('Average Charge State Analysis')
    plt.savefig('Average_Charge_State_Analysis.png', dpi=300)
    plt.show()


plot_results(tryplysc_results, "TrypLysC", 'g', thermolysin_results, "Thermolysin", 'purple')


thermolysin_pro = pd.read_csv('PEAKS\\Thermolysin\\peaks_ptm.proteins.csv')
tryplysc_pro = pd.read_csv('PEAKS\\TrypLysC\\peaks_ptm.proteins.csv')



def calculate_total_proteins(df, enzyme):
    donors = ["Donor1", "Donor2", "Donor3"]
    total_proteins = []
    
    for donor in donors:
        # Filter the columns related to the current donor
        area_columns = [col for col in df.columns if donor in col and "Area" in col]
        # Filter the rows where at least one of the Area columns for this donor is not null
        donor_df = df[df[area_columns].notna().any(axis=1)]
        # Count the unique proteins for this donor
        total_proteins.append(donor_df["Accession"].nunique())
    
    return total_proteins

# Calculate total proteins for each enzyme
thermolysin_proteins = calculate_total_proteins(thermolysin_pro, "Thermolysin")
tryplysc_proteins = calculate_total_proteins(tryplysc_pro, "TrypLysC")

# Calculate means and standard deviations
means = [np.mean(thermolysin_proteins), np.mean(tryplysc_proteins)]
std_devs = [np.std(thermolysin_proteins), np.std(tryplysc_proteins)]

# Plot
plt.bar(["Thermolysin", "TrypLysC"], means, yerr=std_devs, color=["purple", "green"])
plt.ylabel('Total Proteins')
plt.title('Total Proteins Identified per Donor Sample')
plt.savefig('Proteins Identified.png', dpi=300)
plt.show()



def calculate_average_coverage(df, enzyme):
    donors = ["Donor1", "Donor2", "Donor3"]
    average_coverages = []
    
    for donor in donors:
        # Filter the columns related to the current donor
        coverage_columns = [col for col in df.columns if donor in col and "Coverage(%)" in col]
        # Filter the rows where at least one of the Coverage columns for this donor is not null and not zero
        donor_df = df[df[coverage_columns].notna().any(axis=1) & df[coverage_columns].ne(0).any(axis=1)]
        # Calculate the average coverage for this donor
        average_coverages.append(donor_df[coverage_columns].mean().mean())
    
    return average_coverages

# Calculate average coverage for each enzyme
thermolysin_coverages = calculate_average_coverage(thermolysin_pro, "Thermolysin")
tryplysc_coverages = calculate_average_coverage(tryplysc_pro, "TrypLysC")

# Calculate means and standard deviations
means = [np.mean(thermolysin_coverages), np.mean(tryplysc_coverages)]
std_devs = [np.std(thermolysin_coverages), np.std(tryplysc_coverages)]

# Plot
plt.bar(["Thermolysin", "TrypLysC"], means, yerr=std_devs, color=["purple", "green"])
plt.ylabel('Average Coverage (%)')
plt.title('Average Coverage per Donor Sample')
plt.savefig('Coverage.png', dpi=300)
plt.show()



# Filter for proteins identified in at least one donor sample (non-null Area)
thermolysin_pro_filtered = thermolysin_pro[thermolysin_pro.filter(regex='Area').notna().any(axis=1)]
tryplysc_pro_filtered = tryplysc_pro[tryplysc_pro.filter(regex='Area').notna().any(axis=1)]

# Get unique proteins from both datasets
thermolysin_proteins = set(thermolysin_pro_filtered['Accession'])
tryplysc_proteins = set(tryplysc_pro_filtered['Accession'])

# Create the Venn diagram
venn2([thermolysin_proteins, tryplysc_proteins], set_labels = ('Thermolysin', 'TrypLysC'))

plt.show()


# Filter for proteins identified in each donor (non-null Area)
donor_columns = ['Area Donor1_Frac1', 'Area Donor1_Frac2', 'Area Donor1_Frac3', 'Area Donor1_Frac4', 'Area Donor1_Frac5', 'Area Donor1_Frac6',
                 'Area Donor2_Frac1', 'Area Donor2_Frac2', 'Area Donor2_Frac3', 'Area Donor2_Frac4', 'Area Donor2_Frac5', 'Area Donor2_Frac6',
                 'Area Donor3_Frac1', 'Area Donor3_Frac2', 'Area Donor3_Frac3', 'Area Donor3_Frac4', 'Area Donor3_Frac5', 'Area Donor3_Frac6']

for donor in ['Donor1', 'Donor2', 'Donor3']:
    donor_cols = [col for col in donor_columns if donor in col]
    thermolysin_donors = thermolysin_pro[thermolysin_pro[donor_cols].notna().any(axis=1)]['Accession'].nunique()
    tryplysc_donors = tryplysc_pro[tryplysc_pro[donor_cols].notna().any(axis=1)]['Accession'].nunique()
    print(f"For {donor}, {thermolysin_donors} proteins were identified with Thermolysin and {tryplysc_donors} proteins were identified with TrypLysC.")


# List of donors
donors = ['Donor1', 'Donor2', 'Donor3']

# List to store protein counts
protein_counts_thermolysin = []
protein_counts_tryplysc = []

# Filter for proteins identified in each donor (non-null Area)
donor_columns = ['Area Donor1_Frac1', 'Area Donor1_Frac2', 'Area Donor1_Frac3', 'Area Donor1_Frac4', 'Area Donor1_Frac5', 'Area Donor1_Frac6',
                 'Area Donor2_Frac1', 'Area Donor2_Frac2', 'Area Donor2_Frac3', 'Area Donor2_Frac4', 'Area Donor2_Frac5', 'Area Donor2_Frac6',
                 'Area Donor3_Frac1', 'Area Donor3_Frac2', 'Area Donor3_Frac3', 'Area Donor3_Frac4', 'Area Donor3_Frac5', 'Area Donor3_Frac6']

for donor in donors:
    donor_cols = [col for col in donor_columns if donor in col]
    protein_counts_thermolysin.append(thermolysin_pro[thermolysin_pro[donor_cols].notna().any(axis=1)]['Accession'].nunique())
    protein_counts_tryplysc.append(tryplysc_pro[tryplysc_pro[donor_cols].notna().any(axis=1)]['Accession'].nunique())

# Set the width of the bars
barWidth = 0.3

# Set position of bar on X axis
r1 = np.arange(len(protein_counts_thermolysin))
r2 = [x + barWidth for x in r1]

# Create bar plot
plt.figure(figsize=(8,6))
plt.bar(r1, protein_counts_thermolysin, color='purple', width=barWidth, edgecolor='grey', label='Thermolysin')
plt.bar(r2, protein_counts_tryplysc, color='green', width=barWidth, edgecolor='grey', label='TrypLysC')

# Adding xticks
plt.xlabel('Donors', fontweight='bold')
plt.ylabel('Number of Proteins', fontweight='bold')
plt.xticks([r + barWidth/2 for r in range(len(protein_counts_thermolysin))], donors)

plt.legend()
plt.savefig("ProNumber.png",dpi=300)
plt.show()





# Create a DataFrame from the dictionary
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Create a DataFrame from the dictionary
df_dict = {
    "Thermolysin_Shotgun":[165,179,186], 
    "TrypLysC_Shotgun":[201,205,217], 
    "Thermolysin_SCX":[257,268,254],
    "TrypLysC_SCX":[325,330,327],
    "Thermolysin_DIA+SCX":[239,250,249],
    "Thermolysin_DIA+LibFree":[139,138,139],
    "TrypLysC_DIA+SCX":[389,400,397],
    "TrypLysC_DIA+LibFree":[234,235,236]
}



# Convert dictionary to dataframe
df = pd.DataFrame(df_dict)

# Calculate the mean and standard deviation
means = df.mean()
errors = df.std()

# Create a new dataframe for means and errors
df_mean_error = pd.DataFrame({'mean': means, 'error': errors})

# Add enzyme information for color mapping
df_mean_error['Enzyme'] = ['Thermolysin', 'TrypLysC', 'Thermolysin', 'TrypLysC', 'Thermolysin', 'Thermolysin', 'TrypLysC', 'TrypLysC']

# Define color mapping
colors = df_mean_error['Enzyme'].map({'TrypLysC': 'green', 'Thermolysin': 'purple'})

# Create bar plot
plt.figure(figsize=(10,6))
df_mean_error['mean'].plot(kind='bar', yerr=df_mean_error['error'], color=colors, capsize=4)

plt.title('Average Values with Standard Error')
plt.ylabel('# of Proteins')
plt.xticks(rotation=90)
plt.savefig("BarPlot.png",dpi=300)
plt.show()


# Convert dictionary to dataframe
df = pd.DataFrame(df_dict)

# Apply the specified transformations
for key in df_dict.keys():
    if key in ['Thermolysin_SCX', 'TrypLysC_SCX']:
        df[key] = df[key] / (6*75)
    else:
        df[key] = df[key] / 75

# Calculate the mean and standard deviation
means = df.mean()
errors = df.std()

# Create a new dataframe for means and errors
df_mean_error = pd.DataFrame({'mean': means, 'error': errors})

# Add enzyme information for color mapping
df_mean_error['enzyme'] = ['Thermolysin', 'TrypLysC', 'Thermolysin', 'TrypLysC', 'Thermolysin', 'Thermolysin', 'TrypLysC', 'TrypLysC']

# Define color mapping
colors = df_mean_error['enzyme'].map({'TrypLysC': 'green', 'Thermolysin': 'purple'})

# Create bar plot
plt.figure(figsize=(10,6))
df_mean_error['mean'].plot(kind='bar', yerr=df_mean_error['error'], color=colors, capsize=4)

plt.title('Average Values with Standard Error (Normalized)')
plt.ylabel('Protein IDs per Minute')
plt.xticks(rotation=90)
plt.savefig("MachineTime.png",dpi=300)
plt.show()