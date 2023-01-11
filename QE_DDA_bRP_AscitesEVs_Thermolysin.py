import os 
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plot


colours={"TrypLysC":"green","Thermolysin":"purple"}
colour_scheme="PrGn"
dpi=300


### MAX QUANT ANALYSES
def MQ_get_samples(MQ_df):
    x=[x.split("LFQ intensity ")[1] for x in MQ_df.columns if "LFQ" in x]
    print(x)
    return x
#MQ Files
os.listdir("MQ")
os.listdir("MQ/Grouped")
os.listdir("MQ/Seperated")
p2_MQ_proteins=os.path.join("MQ","Seperated","proteinGroups.txt")
p2_MQ_peptides=os.path.join("MQ","Seperated","peptides.txt")
MQ_proteins = pd.read_csv(p2_MQ_proteins, delimiter = "\t")
MQ_peptides = pd.read_csv(p2_MQ_peptides, delimiter="\t")

old_names = MQ_proteins.columns
new_names = [name.replace('FT', 'X') for name in old_names]
MQ_proteins.rename(columns={old_name: new_name for old_name, new_name in zip(old_names, new_names)}, inplace=True)
MQ_samples=MQ_get_samples(MQ_proteins)

MQ_Pro={} #MaxQuant Dictionary of Protein IDs, Sample is Parent Key
MQ_Pep= {} #MaxQuant Dictionary of Peptide IDs, Sample is Parent Key
pro_grab=["Gene names","Sequence"]
pep_grab=[]
i=1
for sample in MQ_samples:
    x=[x for x in MQ_proteins if sample in x]
    print(i, sample)
    i+=1
    print(len(x))
    
    
\



## PEAK ANALYSES
#PEAKS Files
os.listdir("PEAKS")


## FragPipe Analyses
#FragPipe Files
os.listdir("FragPipe")

#Updated





