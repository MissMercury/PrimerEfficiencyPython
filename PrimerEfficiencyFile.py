# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Import excelfile and select columns
qpcr = pd.read_excel('Primer_Eff_Test.xls', sheet_name = 'Results', skiprows = 36)
qpcr = qpcr.iloc[:, np.r_[3, 4, 8, 22, 24, 26]]

#Rename columns
column_names = ["Sample","Gene","CT","Tm1","Tm2","Tm3"]
qpcr.columns = column_names

#CT to numeric and remove undetermined
qpcr = qpcr[qpcr.Sample != 'Blanc']
qpcr.Sample = pd.to_numeric(qpcr.Sample)
qpcr.CT = pd.to_numeric(qpcr.CT)
qpcr['Gene'] = qpcr['Gene'].str.replace(' ', '_')
qpcr = qpcr[qpcr.Tm2 != "NaN"]
qpcr = qpcr[qpcr.Tm3 != "NaN"]


#plot
unique_gene = pd.DataFrame(qpcr["Gene"].unique())
unique = dict(tuple(qpcr.groupby('Gene')))

for key in unique:
    unique[key].plot(kind = 'scatter', x = 'Sample', y = 'CT')
    plt.title([key], loc= 'center')
    plt.show()
    

#Linear model
LinearRegression = pd.DataFrame()

for key in unique:
    result = list(stats.linregress(unique[key].Sample, unique[key].CT))
    
    LinearRegression = pd.DataFrame(LinearRegression, columns = ['slope', 'intercept', 'r_value', 'p_value', 'std_err'])
    LinearRegression.loc[key] = result 


#Efficiency

LinearRegression['Efficiency'] = LinearRegression.slope / (10/3)
LinearRegression['r2'] = LinearRegression.r_value ** 2
LinearRegression['Gene'] = qpcr.Gene.unique()
LinearRegression['PassFail'] = '' 
    
for i in LinearRegression.Gene.unique():
    PassFail = ((LinearRegression.Efficiency[i] >= 0.9) & (LinearRegression.Efficiency[[i]] <= 1.1) & (LinearRegression.r2[[i]] >= 0.975)).bool
    LinearRegression.PassFail[i] = bool(PassFail)


        
        
        
       

