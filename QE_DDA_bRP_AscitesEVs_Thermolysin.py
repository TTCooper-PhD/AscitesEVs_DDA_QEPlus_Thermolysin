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
#MQ Files
os.listdir("MQ")
os.listdir("MQ/Grouped")
os.listdir("MQ/Seperated")
p2_MQ_proteins=os.path.join(os.getcwd, "MQ","Seperated","proteinGroups.txt")
df = pd.read_csv(p2_MQ_proteins, delimiter = "\t")

class MQ():
    
def MQ_get_samples(MQ_df):
    x=[x.split("LFQ intensity ")[1] for x in MQ_df.columns if "LFQ" in x]
    print(x)

def MQ_

get_MQ_samples(df)



## PEAK ANALYSES
#PEAKS Files
os.listdir("PEAKS")


## FragPipe Analyses
#FragPipe Files
os.listdir("FragPipe")

#Updated





