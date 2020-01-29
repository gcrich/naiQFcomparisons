#
# globalFit_NaI-Tl.py
#
# gcrich
#
# made to perform global fits to lit and new QF data for NaI[Tl]
#
from __future__ import print_function
import numpy as np
import csv as csvlib
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import matplotlib.path as path
from scipy.stats import (multivariate_normal, norm)
from scipy.optimize import curve_fit
import random
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.font_manager as font_manager
from numpy.polynomial.polynomial import polyval



fontProps=font_manager.FontProperties()

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


rcParams.update({'figure.autolayout':True})
rcParams['font.family']='sans-serif'
rcParams['font.sans-serif']='Helvetica'
rcParams['legend.fontsize']=16
rcParams['axes.labelsize']=16
rcParams['axes.titlesize']=16
rcParams['xtick.labelsize']=14
rcParams['ytick.labelsize']=14




dataSetNames = ['chagani', 'collar', 'simon', 'spooner', 'gerbier', 'xu',
                'tovey', 'tunl']
dataFilenames= ['./data/qfData_NaI-Tl_' + name + '.dat' for name in dataSetNames]
dataLabels=['Chagani et al.', 'Collar', 'Simon et al.', 'Spooner et al.',
            'Gerbier et al.', 'Xu et al.', 'Tovey et al.', 'TUNL'];

colors=["#377eb8","#e41a1c","#984ea3","#4daf4a","#ff7f00","#a65628","#666666",
        "#f781bf"]

alphas=[0.4, 0.4, 1.0, 1.0, 1.0, 0.8,1.0, 1.0, 1.0]
#markers=['o','s','v','^','o','d','s','x']
markers=['D','p','H','d','o','^','s','x', '^']


iColors = ["#ff7f00","#4daf4a","#666666"]
iMarkers = ['o','d','s']
iMarksizes = [5.5,5.5,7]



datasets = []
under30datasets = []
fullDatasets_noUC = []
datasets_xuOnly = []
under30datasets_xuOnly = []
qfList = []
qfErrorHiList = []
qfErrorLiList = []
eRecoilList = []
eRecoilHiList = []
eRecoilLoList = []
dataFieldNames = ['eRecoil', 'eRecError', 'eRecErrorHi', 'eRecErrorLo', 'qf', 'qfError', 'qfErrorHi', 'qfErrorLo']
for filename,datasetname in zip(dataFilenames,dataSetNames):
    with open(filename, 'r') as qfDatafile:
        csvreader = csvlib.DictReader((row for row in qfDatafile if not row.startswith('#')), delimiter=',',
                                      fieldnames=dataFieldNames,
                                      skipinitialspace=True)
        datalist = []
        datalist_under30 = []
        for row in csvreader:
            newEntry = [float(row[dataname]) for dataname in dataFieldNames]
            if datasetname == 'tunl':
                continue
            datalist.append(newEntry)
            if newEntry[0] < 100 and newEntry[0] > 20:
                datalist_under30.append(newEntry)
        under30datasets.append(np.array(datalist_under30))
        newArray = np.array(datalist)
        if datasetname =='xu':
            newArray = newArray[1:]
            datasets_xuOnly.append(newArray)
            under30datasets_xuOnly.append(np.array(datalist_under30))
        datasets.append(newArray)
        if datasetname != 'collar':
            fullDatasets_noUC.append(newArray)
            
            
# read iodine
iodineNames = ['gerbier', 'spooner', 'tovey']
iodineLabels = ['Gerbier et al.', 'Spooner et al.', 'Tovey et al.']

iodineFilenames = ['./data/iodineRecoils/qfData-' + name + '.dat' for name in iodineNames]
iodineDatasets = []
for filename, datasetname in zip(iodineFilenames, iodineNames):
    with open(filename,'r') as qfDatafile:
        csvreader = csvlib.DictReader((row for row in qfDatafile if not row.startswith('#')), delimiter=',',
                                      fieldnames=dataFieldNames,
                                      skipinitialspace=True)
        datalist=[]
        for row in csvreader:
            newEntry = [float(row[dataname]) for dataname in dataFieldNames]
            datalist.append(newEntry)
        iodineDatasets.append(np.array(datalist))
        


pltArgs={'fontname':'Helvetica','size':14}
fig, ax = plt.subplots(figsize=(8.5,5.25))
zords=[10,10,8,8, 10, 10, 10, 10]
marksizes=[5.5,8,7,7, 7, 7, 7, 7]
for data, color, dataname,mark,datasetname, marksize in zip(datasets, colors, dataLabels,markers, dataSetNames, marksizes):
    #if datasetname == 'collar':
    #    continue
    if datasetname=='tunl':
        continue
#    if datasetname != 'xu' and datasetname != 'tunl':
#        continue
    if all(data[:,7] == 0.) and all(data[:,6] == 0.):
        datayerr = data[:,5]
    else:
        datayerr = [data[:,7], data[:,6]]

    if all(data[:,2] == 0.) and all(data[:,3] == 0.):
        dataxerr = data[:,1]
    else:
        dataxerr = [data[:,3], data[:,2]]
    ax.errorbar(data[:,0], data[:,4], marker=mark, color=color, xerr=dataxerr,
                 yerr=datayerr,
                 label=dataname, capsize=0, linestyle='None',
                 markeredgewidth=0, markersize=marksize )


ax.set_ylim(0.0, 60)
ax.set_xscale('log',nonposx='clip')
ax.legend(loc='upper left', 
          frameon=False,
          ncol=2,
          prop={'family':'Helvetica'},
          numpoints=1,
          labelspacing=0.25)

ax.set_ylabel('Quenching factor (%)')
ax.set_xlabel('Nuclear recoil energy (keV)')


plt.draw()
plt.savefig('globalQF_NaI-Tl.pdf')
#plt.show()




pltArgs={'fontname':'Helvetica','size':14}
fig, ax = plt.subplots(figsize=(8.5,5.25))
zords=[10,10,8]

for data, color, dataname,mark,datasetname, marksize in zip(iodineDatasets, iColors, iodineLabels, iMarkers, iodineNames, iMarksizes):
    #if datasetname == 'collar':
    #    continue
    if datasetname=='tunl':
        continue
#    if datasetname != 'xu' and datasetname != 'tunl':
#        continue
#    if all(data[:,7] == 0.) and all(data[:,6] == 0.):
#        datayerr = data[:,5]
#    else:
#        datayerr = [data[:,7], data[:,6]]
#
#    if all(data[:,2] == 0.) and all(data[:,3] == 0.):
#        dataxerr = data[:,1]
#    else:
#        dataxerr = [data[:,3], data[:,2]]
    datayerr = data[:,5] * 100
    dataxerr = data[:,1]
    ax.errorbar(data[:,0], 100*data[:,4], marker=mark, color=color, xerr=dataxerr,
                 yerr=datayerr,
                 label=dataname, capsize=0, linestyle='None',
                 markeredgewidth=0, markersize=10 )


ax.set_ylim(0.0, 14)
ax.set_xlim(0, 140)
#ax.set_xscale('log',nonposx='clip')
ax.legend(loc='upper right', 
          frameon=False,
          ncol=1,
          prop={'family':'Helvetica'},
          numpoints=1,
          labelspacing=0.25)

ax.set_ylabel('Quenching factor (%)')
ax.set_xlabel('Nuclear recoil energy (keV)')

plt.draw()
plt.savefig('iodineQF.pdf')