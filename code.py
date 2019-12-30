import pandas as pd
from scipy import stats
import math
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------------------------------------------#

# solution 1 is done as asked in the question.
# solution 2 histogram is attached
# solution 3 n0 is taken same as n because of the histogram observed and the discussion done in the class.
# solution 4 is done based on the p values.
# solution 5 and 6 are done as per asked
# solution 7 is :Women Smokers vs non-Smokers increase Genes  ['SULT1A1', 'AOC2', 'CYP2S1', 'HNF4A', 'PNKP', 'PTPN6', 'HLA-C', 'HLA-E', 'HLA-G']
#                Women Smokers vs non-Smokers decrease Genes  ['AADAC', 'AS3MT', 'IFNG', 'KLRC2', 'PRF1']
#                Men Smokers vs non-Smokers increase Genes  ['AADAC', 'HNF4A', 'AS3MT', 'IFNG', 'KLRC2', 'PRF1', 'HLA-E', 'HLA-G']
#                Men Smokers vs non-Smokers descrease Genes  ['SULT1A1', 'AOC2', 'CYP2S1', 'HNF4A', 'PNKP', 'PTPN6', 'HLA-C', 'HLA-E', 'HLA-G']

# --------------------------------------------------------------------------------------------------------------#

data1 = pd.read_csv('../data/Raw Data_GeneSpring.txt',sep='\t')
data2 = pd.read_csv('../data/NKCellCytotoxicity.txt', sep='\t')
#data3 = pd.read_csv('/content/gdrive/My Drive/Colab Notebooks/data/FreeRadicalResponse.txt', header = None)
#data4 = pd.read_csv('/content/gdrive/My Drive/Colab Notebooks/data/XenobioticMetabolism1.txt', header = None)
#data5 = pd.read_csv('/content/gdrive/My Drive/Colab Notebooks/data/DNARepair1.txt', header = None)

data1_list= data1.iloc[:,1:49].values.tolist()
gene=data1.iloc[:,49:50].values.tolist()



A=np.zeros((48,4))

for i in range (4):
    for j in range(12):
        A[j+(i*12)][i]=1

B=np.zeros((48,4))

for i in range (24):
    B[i][0]=1
for i in range (24):
    B[24+i][1]=1
for i in range (12):
    B[i][2]=1
for i in range (12):
    B[i+24][2]=1
for i in range (12):
    B[i+12][3]=1
for i in range (12):
    B[i+36][3]=1

A=np.matrix(A)
B=np.matrix(B)

# for i in data1.iterrows():
#     print(i)

F=[]
for j in data1_list:
    i= np.matrix(j)
    temp= np.dot(np.dot(A,np.linalg.pinv(np.dot(A.T,A))),A.T)
    temp1=np.dot(np.dot(B,np.linalg.pinv(np.dot(B.T,B))),B.T)
    temp2=np.subtract(temp , temp1)
    #print(i)
    #print(temp2.shape,i.shape)
    numerator = np.dot(np.dot(i , temp2) , i.T)
    I=np.identity(48)
    denominator=np.dot(np.dot(i , np.subtract(I,temp)),i.T)
    val=float(numerator/denominator)*((48-np.linalg.matrix_rank(A))/(np.linalg.matrix_rank(A)-np.linalg.matrix_rank(B)))
    F.append(val)
#print(F)

one=48-np.linalg.matrix_rank(A)
two=np.linalg.matrix_rank(A)-np.linalg.matrix_rank(B)
p=1-stats.f.cdf(F,two,one)
# print(p)
update_p=[]
for i in p:
    if (i<1):
        update_p.append(i)

plt.hist(update_p,bins=100)
plt.xlabel("updated_p_value")
plt.ylabel("Frequency")
plt.savefig('hist.png')
plt.show()


gene_dist=[]
for i in range(len(p)):
    if p[i] < 0.05:
        gene_dist.append(gene[i][0])
#print(len(gene_dist))

with open('../data/NKCellCytotoxicity.txt','r') as file:
    with open('../data/temp.txt','w') as output:
        next(file)
        next(file)
        for row in file:
            output.write(row)
data2 = pd.read_csv('../data/temp.txt',header=None)
data2=data2.loc[1:]
q2=set(data2[0])
#print(len(q2))
qq=set(gene_dist)
#print(len(qq))
qqq2=q2.intersection(qq)
#print(len(qqq2))

with open('../data/FreeRadicalResponse.txt','r') as file:
    with open('../data/temp.txt','w') as output:
        next(file)
        next(file)
        for row in file:
            output.write(row)
data3 = pd.read_csv('../data/temp.txt',header=None)
q3=set(data3[0])
#print(len(q3))
qq=set(gene_dist)
#print(len(qq))
qqq3=q3.intersection(qq)
#print(len(qqq3))

with open('../data/XenobioticMetabolism1.txt','r') as file:
    with open('../data/temp.txt','w') as output:
        next(file)
        next(file)
        for row in file:
            output.write(row)
data4 = pd.read_csv('../data/temp.txt',header=None)
q4=set(data4[0])
#print(len(q4))
qqq4=q4.intersection(qq)
#print(len(qqq4))

with open('../data/DNARepair1.txt','r') as file:
    with open('../data/temp.txt','w') as output:
        next(file)
        next(file)
        for row in file:
            output.write(row)
data5 = pd.read_csv('../data/temp.txt',header=None)
q5=set(data5[0])
#print(len(q5))
qqq5=q5.intersection(qq)
#print(len(qqq5))

total_int=list(qqq2) + list(qqq3) + list(qqq4) + list(qqq5)
#print(len(set(total_int)))


a=[]
for i in range(len(total_int)):
    for j in range(len(gene)):
        #print(1)
        if (total_int[i] == gene[j][0]):
            a.append(j)

list1=[]
list2=[]
list3=[]
list4=[]
for i in a:
    l1=np.mean(data1_list[i][:12])
    l2=np.mean(data1_list[i][12:24])
    l3=np.mean(data1_list[i][24:36])
    l4=np.mean(data1_list[i][36:48])
    if(l1>l2):
        list2.extend(gene[i])
    else:
        list1.extend(gene[i])
    if(l3>l4):
        list4.extend(gene[i])
    else:
        list3.extend(gene[i])

print('Women Smokers vs non-Smokers increase Genes' )
print(set(list3))
print('Women Smokers vs non-Smokers decrease Genes' )
print(set(list4))
print('Men Smokers vs non-Smokers increase Genes')
print(set(list1))
print('Men Smokers vs non-Smokers decrease Genes')
print(set(list2))

