#!/usr/bin/env python3

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas
import sys
import os.path

homozygous = False
for x in sys.argv:
  if x == '--homozygous':
    homozygous = True
  elif x == "-h":
    print("variant_heatmap.py [--homozygous]")
    exit(0)

if (not os.path.exists("output_sorted")):
  print("Cant find output_sorted directory")
  exit(0)

covs = [30, 50, 100, 200, 500, 1000, 10000]
covLabels = [str(x) for x in covs]
covLabels[len(covLabels) - 1] = 'All'
countFileList = ['output_sorted/counts_'+str(x)+'.txt' for x in covs]

if homozygous:
  countFileList = ['output_sorted/counts_homo_'+str(x)+'.txt' for x in covs]
numSamples = len(countFileList)

frequencies = []

for i in range(0, numSamples):
  dict = {}
  with open(countFileList[i]) as f:
    for line in f:
      vals = line.split()
      dict[int(vals[0])] = float(vals[1])
  frequencies.append(dict)
  print(dict)
  
allKeys = {}
keyList = []
for dict in frequencies:
  for a in dict:
    if not a in allKeys:
      keyList.append(a)
    allKeys[a] = True
keyList.sort()
print(keyList)

heatmap = []

for i in range(0, numSamples):
  dict = frequencies[i]
  heatmaprow = []
  for val in keyList:
    if val in dict:
      heatmaprow.append(dict[val])
    else:
      heatmaprow.append(0.0)
  heatmap.append(heatmaprow)
print(heatmap)

heatmapPandas = pandas.DataFrame(heatmap)
sns.set()
fig, ax = plt.subplots(figsize=(20,6))         # Sample figsize in inches
colors = sns.cubehelix_palette(11)

allFull = True
for h in heatmap:
  if min(h) < 1.0:
    allFull = False
if allFull:
  colors = colors[10:]
sns.heatmap(heatmapPandas, xticklabels=True, yticklabels=True, cmap=colors, ax = ax)
#PiYG is also decent 
ax.set_xticklabels(keyList, fontdict={'fontsize': 8})
ax.set_yticklabels(covLabels)
ax.hlines(np.arange(0, numSamples+1), *ax.get_xlim())
ax.vlines(np.arange(0, len(keyList)+1), *ax.get_ylim())
#ax.xaxis.tick_top()
plt.xticks(rotation=90)
plt.xlabel('SNP Position', labelpad = 10)
plt.ylabel('Coverage Filter')
#plt.savefig('snps.png')
plt.title('SNP Support in Downsampled Reads (JHU004)')
if homozygous:
  plt.title('Homozygous SNP Support in Downsampled Reads (JHU004)')
plt.gcf().subplots_adjust(bottom=0.15)
fn = ''
if len(sys.argv) > 1:
  fn = sys.argv[1]
if len(fn) > 0:
  plt.savefig(fn)
else:
  plt.show()
    
  
