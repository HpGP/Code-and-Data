import pandas as pd
output = open('HpGP_6g_map.txt','w') #open file
data = pd.read_excel("HpGP_coancestry_populations.xlsx") #open dataset
freq_grp=data[["COUNTRY","DAPC groups"]] #subset
popul=data["COUNTRY"] #only populations
pop=[] #pop array
ngroups = 6 #number of groups to count
for el in range(len(popul)): 
    if popul[el] not in pop:
        pop.append(popul[el])
for r in pop: 
    filteredPop = freq_grp.loc[freq_grp["COUNTRY"]==r] 
    riga = r + '\t' #pop name + tab
    for i in range(1, ngroups + 1): #groups
        filterGroup = filteredPop.loc[freq_grp["DAPC groups"]==i] 
        riga+= str(len(filterGroup)) + '\t' #add lenght of filtered group + tab
    output.write(riga+'\n')
    output.flush()
output.close()