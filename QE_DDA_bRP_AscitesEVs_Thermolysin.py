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
p2_MQ_proteins
df = pd.read_csv("C:\Users\tcoop\Desktop\AscitesEVs_DDA_QEPlus_Thermolysin\MQ_\Seperated\proteinGroups.txt", delimiter = "\t")
print(df)


## PEAK ANALYSES
#PEAKS Files
os.listdir("PEAKS")


## FragPipe Analyses
#FragPipe Files
os.listdir("FragPipe")

#Updated





