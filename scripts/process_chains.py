#!/usr/bin/python
# nano process_chains.py

import re
import os
import glob
import sys
import pandas as pd
import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from scipy.stats import gaussian_kde
import matplotlib.patheffects as path_effects
from itertools import product
import arviz


def hpd(data, level):
    """
    Return highest posterior density interval from a list,
    given the percent posterior density interval required.
    """
    d = list(data)
    d.sort()
    nData = len(data)
    nIn = int(round(level * nData))
    if nIn < 2 :
        return None
    #raise RuntimeError("Not enough data. N data: %s"%(len(data)))
    i = 0
    r = d[i+nIn-1] - d[i]
    for k in range(len(d) - (nIn - 1)) :
        rk = d[k+nIn-1] - d[k]
        if rk < r :
            r = rk
            i = k
    assert 0 <= i <= i+nIn-1 < len(d)
    return (d[i], d[i+nIn-1])


def effectiveSampleSize(data, stepSize = 1) :
    """ Effective sample size, as computed by BEAST Tracer."""
    samples = len(data)
    assert len(data) > 1,"no stats for short sequences"
    maxLag = min(samples//3, 1000)
    gammaStat = [0,]*maxLag
    varStat = 0.0;
    if type(data) != np.ndarray :
        data = np.array(data)
    normalizedData = data - data.mean()
    for lag in range(maxLag) :
        v1 = normalizedData[:samples-lag]
        v2 = normalizedData[lag:]
        v = v1 * v2
        gammaStat[lag] = sum(v) / len(v)
        if lag == 0 :
            varStat = gammaStat[0]
        elif lag % 2 == 0 :
            s = gammaStat[lag-1] + gammaStat[lag]
            if s > 0 :
                varStat += 2.0*s
            else :
                break
    # auto correlation time
    act = stepSize * varStat / gammaStat[0]
    # effective sample size
    ess = (stepSize * samples) / act
    return ess

def DIC(dev):
    """ takes the deviance at each sample and returns the deviance information criterion value"""
    DIC=0.5*np.var(dev)+np.mean(dev)
    return DIC

def colourDict(data,cmap=mpl.cm.viridis,sort=False):
    """
    returns a dictionary of unique data_items:colour hex code, normalised
    takes a list of items and a cmap as mpl.cm.name. If sort = False, the function takes an ordered list of unique elements
    """
    cmap=cmap # default viridis
    data_unique=data if sort == False else list(set(data)) # includes nan values as a key. That is desirable sometimes
    norm = mpl.colors.Normalize(vmin=0, vmax=len(data_unique))
    colors = [cmap(norm(value)) for value in range(len(data_unique))]
    return dict(zip(data_unique,colors))

### print file argument information
input_path=os.getcwd()
output_path=os.getcwd()

print('Working directory and output path:\n','input: ',input_path,'\noutput:',output_path)

# lazy and running out of time
del sys.argv[0]
filearray=sys.argv

print('string model list')
print([str(x).split('-')[1] for x in filearray])
print('set')
print(set([str(x).split('-')[1] for x in filearray]))
print('keys')

models=[y for y in set([str(x).split('-')[1] for x in filearray])]
print(models)

print('Number of files: ', len(filearray))
print('Files: ',[str(x) for x in filearray])

print('Number of models: ', len(models))
print('Models: ',[str(x) for x in models])


# get file headers. All files must be the same model and thus shuld have the same header
with open(os.path.join(input_path,'%s'%(filearray[0])), 'r') as f:
    header=f.readline().rstrip().split('\t')

# print(header)

# create empty dataframe
alldata=pd.DataFrame(columns=header+['chain','model']) # new columns in other models will be automatically added anyways

# for file in filearray:

for index in range(len(filearray)):
    filename=glob.glob(os.path.join(input_path,'%s'%(filearray[index])))
    model=filearray[index].split('-')[1]
    print('Loading: ',filename)
    temp=pd.read_csv(filename[0],sep='\t')
    temp.loc[temp.index,'chain']=int(filearray[index].split('-')[2].replace('.txt.out',''))
    temp.loc[temp.index,'model']=str(model)
    alldata=pd.concat([alldata,temp],sort=False) # all columns should be in the same order
    print('Rows by now: %s  '%(len(alldata)),'Rows added: %s  '%(len(temp)))

print('All data loaded with %s rows'%(len(alldata)))

### chain plots all chains separatedly for better assessment
spinescol='#999999'

chain_coldict=colourDict([int(x) for x in alldata['chain'].unique()],cmap=mpl.cm.viridis,sort=True)

print('plotting variables with separate chains')

for mod in models:
    temp1=alldata[alldata['model']==mod]
    colnames=list(temp1.columns)
    for x in ['index','Intercept','chain','model']:
        colnames.remove(x)
    cols=len(colnames)
    rows=len(temp1['chain'].unique())
    plt.figure(figsize=(30*cols,3*rows),facecolor='w')
    G = gridspec.GridSpec(rows,cols,hspace=0.05,wspace=0.05) # rows, columns

    for colindex,col in enumerate(colnames):
        for chainindex,chain in enumerate(alldata['chain'].unique()):
            temp2=temp1[temp1['chain']==chain]
            # print('indices col,row',colindex,chainindex)
            # print('data',col,chain)
            ax=plt.subplot(G[chainindex,colindex],facecolor='none') # columns,rows colindex,chainindex
            ax.plot([x for x in range(len(temp2))],temp2[col],c=chain_coldict[chain],lw=0.8)
            plt.ylabel('Chain %s'%(int(chainindex)), fontsize=14)
            plt.xlabel(str(col), fontsize=14)
            plt.yticks(fontsize = 12, color=spinescol); plt.xticks(fontsize = 12, color=spinescol)
            plt.title(mod, fontdict={'fontsize':15}, loc='center')
            plt.grid(axis='y')

    plt.savefig('./%s_indiv_chains.jpg'%(mod),dpi=200)
    plt.close()

# lal chains overlaped
print('plotting variables with overlapping chains')

for mod in models:
    temp1=alldata[alldata['model']==mod]
    colnames=list(temp1.columns)
    for x in ['index','Intercept','chain','model']:
        colnames.remove(x)
    cols=len(colnames)

    plt.figure(figsize=(30,3*cols),facecolor='w')
    G = gridspec.GridSpec(cols,1,hspace=0.05,wspace=0.05)

    for colindex,col in enumerate(colnames):
        ax=plt.subplot(G[colindex,0],facecolor='none')
        for chain in alldata['chain'].unique():
            temp2=temp1[temp1['chain']==chain]
            ax.plot([x for x in range(len(temp2))],temp2[col],c=chain_coldict[chain],lw=0.8)
            plt.ylabel(str(col), fontsize=14)
            plt.yticks(fontsize = 12, color=spinescol); plt.xticks(fontsize = 12, color=spinescol)
            plt.title(mod, fontdict={'fontsize':15}, loc='center')
            plt.grid(axis='y')
        colindex+=1

    plt.savefig('./%s_all_chains.jpg'%(mod),dpi=200)
    plt.close()

print('plots saved')

print('calculating HPDs, ESS, and data proportions')

num_cols=[x for x in alldata.columns]
for x in ['index','Intercept','chain','model']:
    num_cols.remove(x)

alldescdict={}
for mod in models:
    temp1=alldata[alldata['model']==mod][num_cols].dropna(axis=1)
    desc=temp1.describe(percentiles=[.025, .5, .975]).T.copy(deep=True)
    desc['median']=temp1.median()
    for col in list(temp1.columns):
        desc.loc[col,'ESS']=effectiveSampleSize(list(temp1[col].dropna().values))
        desc.loc[col,'HPD_025']=hpd(list(temp1[col].dropna().values),0.95)[0]
        desc.loc[col,'HPD_095']=hpd(list(temp1[col].dropna().values),0.95)[1]
        desc.loc[col,'samples_above0']=len(temp1[temp1[col]>0])
        desc.loc[col,'samples_total']=len(temp1[col])
        desc.loc[col,'prop_above0']=len(temp1[temp1[col]>0])/len(temp1[col])
    alldescdict[mod]=desc

alldesclist=pd.concat(alldescdict.values(),keys=alldescdict.keys(), sort=False)
alldesclist.reset_index(inplace=True)
alldesclist.to_csv('./Allmodels_chains_summary.txt',sep='\t') # save it to file
print('done')
