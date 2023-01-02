#pd.options.display.max_rows = None
from itertools import combinations
import re
import Peptide_Chef as Chef
from pyteomics import parser, electrochem
import math
import pandas as pd
import numpy as np

#Statistics
import statsmodels.api as sm
from statsmodels.formula.api import ols

#Figure Generation
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib import offsetbox
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d import axes3d
from IPython.display import Image, display
import seaborn as sns
from adjustText import adjust_text
import glob

#Venn Diagrams
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles


cmap = 'PRGn'
fmt='eps'
dpi=600


    
    

#generate library of human proteome
url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
menu=Chef.CookBook(homebrew=False,takeout=True,url=url, measure=True, target="Peptide")
                #### Top EV Proteins (cite Kugerastski et al, 2021 Nature Cell Biology)

ascites_vesiclepedia=pd.read_excel(io="Vesiclepedia_Ascites.xlsx")
EV_List=list(set(ascites_vesiclepedia["Content"]))
                
# EV_List= ["SDCB1","AAAT","4F2","CLH1","GBB1","CD47","GBB2","ITB1","BASI",
#           "AT1A1","RAP1B","GNAI3","LG3BP","GNAI2","RAB10",
#           "PDC6I","RS27A","G3PT","TS101","RRAS","EF2","RAN","CD9","CD63","CD81","HS90A","KPYM"]
ev_menu=menu.loc[menu['Gene'].isin(EV_List)]
ev_dict = dict(zip(ev_menu['Peptide'], ev_menu['UniprotID']))


## Trypsin (In Silico Digestion)
missed_sites=[0,1,2,3,4,5]
enzyme="trypsin"
Tryp_Peptides=[]
for site in missed_sites:
    x=f"Tryp_{site}_dig"
    z=Chef.ButcherShop(ev_menu,target="Peptide",identifier="UniprotID", rule=enzyme,min_charge=0.01,missed=site,exception=None,min_length=6,max_length=65)
    globals()[str(x)]=Chef.Deli(z,meat_package=True)
    Tryp_Peptides.append(globals()[str(x)])
i=0
columns=["0","1","2","3","4","5"]
Trypsin_Coverage=[]
Tryp_Silico_Coverage=pd.DataFrame()
for df in Tryp_Peptides:
    name=f"{columns[i]}" 
    globals()[str(name)]={}
    for peptide, gene in ev_dict.items():
        temp=tuple(df.loc[df.gene == gene]["peptide"])
        x=Chef.Pep2Pro(peptide,temp)
        globals()[str(name)][gene]=x
    temp=list(globals()[str(name)].values())
    Tryp_Silico_Coverage[name]= np.array(temp)
    i+=1
                                        ## Coverage Figure 
                                        
sns.set(rc={'figure.figsize':(8,4)})
sns.set_style("ticks")
meanlineprops = dict(linestyle='--', linewidth=2.5, color='black')
boxprops = dict( linestyle='-',linewidth=2, color='#CDCDCD')
plt.tick_params(axis='x', colors='black', length=10, width=2)
plt.tick_params(axis='y', colors='black', length=10, width=2)
sns.stripplot(data=Tryp_Silico_Coverage,color='green',s=5)
sns.boxplot(data=Tryp_Silico_Coverage,showmeans=True, meanline=True,
            boxprops=boxprops,meanprops=meanlineprops,showfliers=True,color="black",linewidth=1.5)
plt.ylabel("Sequence Coverage",fontsize=22)
plt.xlabel("Missed Cleavage Sites",fontsize=22)
plt.xticks(fontsize=14,rotation=0)
plt.yticks(fontsize=14)
plt.savefig('TrypLysC_EVs_SC.eps',dpi=dpi,bbox_inches="tight")
plt.show() 
    
## Thermolysin (in silico digestion)

missed_sites=[0,1,2,3,4,5]
enzyme="thermolysin"
Thermo_Peptides=[]
for site in missed_sites:
    x=f"Thermo_{site}_dig"
    z=Chef.ButcherShop(ev_menu,target="Peptide",identifier="UniprotID", rule=enzyme,min_charge=0.01,missed=site,exception=None,min_length=6,max_length=65)
    globals()[str(x)]=Chef.Deli(z,meat_package=True)
    Thermo_Peptides.append(globals()[str(x)])
i=0
columns=["0","1","2","3","4","5"]
Thermolysin_Coverage=[]
Thermo_Silico_Coverage=pd.DataFrame()
for df in Thermo_Peptides:
    name=f"{columns[i]}" 
    globals()[str(name)]={}
    for peptide, gene in ev_dict.items():
        temp=tuple(df.loc[df.gene == gene]["peptide"])
        x=Chef.Pep2Pro(peptide,temp)
        globals()[str(name)][gene]=x
    temp=list(globals()[str(name)].values())
    Thermo_Silico_Coverage[name]= np.array(temp)
    i+=1 
                ## Coverage Figure 
sns.set(rc={'figure.figsize':(8,4)})
sns.set_style("ticks")
meanlineprops = dict(linestyle='--', linewidth=2.5, color='black')
boxprops = dict( linestyle='-',linewidth=2, color='#CDCDCD')
plt.tick_params(axis='x', colors='black', length=10, width=2)
plt.tick_params(axis='y', colors='black', length=10, width=2)
sns.stripplot(data=Thermo_Silico_Coverage,color='purple',s=5)
sns.boxplot(data=Thermo_Silico_Coverage,showmeans=True, meanline=True,
            boxprops=boxprops,meanprops=meanlineprops,showfliers=True,color="black",linewidth=1.5)
plt.ylabel("Sequence Coverage",fontsize=22)
plt.xlabel("Missed Cleavage Sites",fontsize=22)
plt.xticks(fontsize=14,rotation=0)
plt.savefig('Thermolysin_EVs_SC.eps',dpi=dpi,bbox_inches="tight")
plt.show()

### Plot # of Peptides per Missed Cleavage


Tryp_Pep_Counts=[]
for df in Tryp_Peptides:
    log = math.log2(len(df))
    Tryp_Pep_Counts.append(log)
    
Thermo_Pep_Counts=[]

for df in Thermo_Peptides:
    log = math.log2(len(df))
    Thermo_Pep_Counts.append(log)
    
plt.figure(figsize=(8,4))
plt.plot(missed_sites,Tryp_Pep_Counts,color="green")
plt.plot(missed_sites,Thermo_Pep_Counts,color="purple")
plt.xlabel('Missed Cleavages per Proteins', fontsize=18)
plt.ylabel('Log2(Number of Peptides)', fontsize=18)
plt.xticks(missed_sites,fontsize=14)
plt.yticks(np.arange(12,20,1),fontsize=14)
plt.savefig("Silico_EV_Peptides.eps",dpi=600,bbox_inches='tight')
plt.show()

#### Peptide Characteristics

enzymes=["TrypLysC","Thermolysin"]
numb=[0,2,5]
ev_pep_df_list=[]
dfs=[Tryp_Peptides,Thermo_Peptides]
x=0

for df in dfs:
    enzyme=enzymes[x]
    x+=1
    for i in numb:
        temp=df[i].loc[:,("peptide","Length","z","Mass","m/z")]
        temp['ID'] = np.resize(f"{enzyme}_{i}",len(temp))
        ev_pep_df_list.append(temp)
ev_pep_df=pd.concat(ev_pep_df_list,ignore_index=True)
Chef.Marinate(ev_pep_df,"peptide",length=0,IPC=True,Hydro=True,GRAVY=True,NeutralZ=True,Peptide_Inspector=False)
ev_colors={"TrypLysC_0":"green",
           "TrypLysC_2":"green",
           "TrypLysC_5":"green",
           "Thermolysin_0":"purple",
           "Thermolysin_2":"purple",
           "Thermolysin_5":"purple",
          }
sns.set(rc={'figure.figsize':(8,4)})
sns.set_style("ticks")
meanlineprops = dict(linestyle='--', linewidth=2.5, color='yellow')
boxprops = dict( linestyle='-',linewidth=2, color='#CDCDCD')
plt.tick_params(axis='x', colors='black', length=10, width=2)
plt.tick_params(axis='y', colors='black', length=10, width=2)
sns.stripplot(x="ID", y="Hydro_Sum",palette=ev_colors,data=ev_pep_df,s=1.5)
sns.boxplot(x="ID", y="Hydro_Sum",data=ev_pep_df,palette=ev_colors,showmeans=True, meanline=True,
            boxprops=boxprops,meanprops=meanlineprops,showfliers=False,color="black")
# plt.title('Sequence Coverage of All Proteins',fontname='Times New Roman',fontweight='bold',fontsize=20,pad=30,backgroundcolor='#cbe7e3',color='black',style='italic');
plt.ylabel("Peptide GRAVY Score",fontsize=24)
plt.xlabel("Missed Cleavages",fontsize=24)
plt.xticks(fontsize=14,rotation=90)
plt.yticks(np.arange(-5,5,1),fontsize=14)
plt.axhline(y=0.3, linestyle="--",color='black',linewidth=2)
plt.axhline(y=-0.4, linestyle="--",color='black',linewidth=2)
plt.savefig("ev_pep_GRAVY.eps",dpi=dpi,bbox_inches="tight")
plt.show()

sns.set(rc={'figure.figsize':(8,4)})
sns.set_style("ticks")
meanlineprops = dict(linestyle='--', linewidth=2.5, color='yellow')
boxprops = dict( linestyle='-',linewidth=2, color='#CDCDCD')
plt.tick_params(axis='x', colors='black', length=10, width=2)
plt.tick_params(axis='y', colors='black', length=10, width=2)
sns.stripplot(x="ID", y="Length",palette=ev_colors,data=ev_pep_df,s=1.5)
sns.boxplot(x="ID", y="Length",data=ev_pep_df,palette=ev_colors,showmeans=True, meanline=True,
            boxprops=boxprops,meanprops=meanlineprops,showfliers=False,color="black")
# plt.title('Sequence Coverage of All Proteins',fontname='Times New Roman',fontweight='bold',fontsize=20,pad=30,backgroundcolor='#cbe7e3',color='black',style='italic');
plt.ylabel("Peptide Length",fontsize=24)
plt.xlabel("Missed Cleavages",fontsize=24)
plt.xticks(fontsize=14,rotation=90)
plt.yticks(fontsize=14)
plt.axhline(y=7, linestyle="--",color='black',linewidth=2)
plt.savefig("ev_pep_Length.eps",dpi=dpi,bbox_inches="tight")
plt.show()

sns.set(rc={'figure.figsize':(8,4)})
sns.set_style("ticks")
meanlineprops = dict(linestyle='--', linewidth=2.5, color='yellow')
boxprops = dict( linestyle='-',linewidth=2, color='#CDCDCD')
plt.tick_params(axis='x', colors='black', length=10, width=2)
plt.tick_params(axis='y', colors='black', length=10, width=2)
sns.stripplot(x="ID", y="z",palette=ev_colors,data=ev_pep_df,s=1.5)
sns.boxplot(x="ID", y="z",data=ev_pep_df,showmeans=True,palette=ev_colors, meanline=True,
            boxprops=boxprops,meanprops=meanlineprops,showfliers=False,color="black")
# plt.title('Sequence Coverage of All Proteins',fontname='Times New Roman',fontweight='bold',fontsize=20,pad=30,backgroundcolor='#cbe7e3',color='black',style='italic');
plt.ylabel("Peptide Charge",fontsize=24)
plt.xlabel("Missed Cleavages",fontsize=24)
plt.xticks(fontsize=14,rotation=90)
plt.yticks(np.arange(0,21,2),fontsize=14)
plt.axhline(y=2, linestyle="--",color='black',linewidth=2)
plt.savefig("ev_pep_Charge.eps",dpi=dpi,bbox_inches="tight")
plt.show()

                ## Peptides per Protein - EV in silico figures
                
numb=[0,2,5]
ev_pep_count_list=[]
dfs=[Tryp_Peptides,Thermo_Peptides]
enzymes=["TrypLysC","Thermolysin"]
x=0

for df in dfs:
    enzyme=enzymes[x]
    x+=1
    for i in numb:
        temp=df[i]
        temp = temp.groupby('gene').first().reset_index()
        temp = temp.loc[:,("gene","counts")]
        temp['ID'] = np.resize(f"{enzyme}_{i}",len(temp))
        ev_pep_count_list.append(temp)
ev_pep_count=pd.concat(ev_pep_count_list,ignore_index=True)
sns.set(rc={'figure.figsize':(8,4)})
sns.set_style("ticks")
meanlineprops = dict(linestyle='--', linewidth=2.5, color='yellow')
boxprops = dict( linestyle='-',linewidth=2, color='#CDCDCD')
plt.tick_params(axis='x', colors='black', length=10, width=2)
plt.tick_params(axis='y', colors='black', length=10, width=2)
sns.stripplot(x="ID", y="counts",palette=ev_colors,data=ev_pep_count,s=3)
sns.boxplot(x="ID", y="counts",data=ev_pep_count,palette=ev_colors,showmeans=True, meanline=True,
            boxprops=boxprops,meanprops=meanlineprops,showfliers=False,color="black")
# plt.title('Sequence Coverage of All Proteins',fontname='Times New Roman',fontweight='bold',fontsize=20,pad=30,backgroundcolor='#cbe7e3',color='black',style='italic');
plt.ylabel("Peptides per Protein",fontsize=24)
plt.xlabel("Missed Cleavages",fontsize=24)
plt.xticks(fontsize=14,rotation=90)
# plt.yticks(np.arange(0,21,2),fontsize=14)
# plt.axhline(y=2, linestyle="--",color='black',linewidth=2)
plt.savefig("ev_pep_count.eps",dpi=dpi,bbox_inches="tight")
plt.show()

