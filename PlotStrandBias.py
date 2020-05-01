#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys

if (len(sys.argv) != 5):
  print("PlotStrandBias.py sbfull.txt sb_old.txt sb_new.txt out.png");
  exit(0);

allFn = sys.argv[1]
oldFn = sys.argv[2]
newFn = sys.argv[3]
outFn = sys.argv[4]

allSb = []
oldSb = []
newSb = []

with open(allFn) as f:
  for line in f:
    allSb.append(float(line))

with open(oldFn) as f:
  for line in f:
    oldSb.append(float(line))

with open(newFn) as f:
  for line in f:
    newSb.append(float(line))
    
sns.set()

fig, axs = plt.subplots(ncols=3)
sns.distplot(allSb, kde = False, ax=axs[0])
sns.distplot(oldSb, kde = False, ax=axs[1])
sns.distplot(newSb, kde = False, ax=axs[2])
fig.suptitle('Strand Bias')
axs[0].title.set_text('Full Dataset')
axs[1].title.set_text('Sample (old)')
axs[2].title.set_text('Sample (even_strand)')
axs[1].set_xlabel('+ strand read proportion')
axs[0].get_yaxis().set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[2].get_yaxis().set_visible(False)
plt.subplots_adjust(wspace=0.35)
plt.savefig(outFn)
