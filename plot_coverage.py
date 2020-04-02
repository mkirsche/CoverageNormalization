#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np

if (len(sys.argv) != 2):
  print("plot_coverage.py oldcov.newcov.txt")
  exit(0)

fn = sys.argv[1]
with open(fn) as f:
  oldcovsample = []
  newcovsample = []
  oldcov = []
  newcov = []
  lineNum = 0
  for line in f:
    tokens = line.split()
    if lineNum%100 == 0:
      oldcovsample.append(int(tokens[0]))
      newcovsample.append(int(tokens[1]))
    oldcov.append(int(tokens[0]))
    newcov.append(int(tokens[1]))
    lineNum += 1
    
sns.set()
xs = np.arange(len(oldcovsample)) * 100
plot_ = sns.barplot(x=xs, y = oldcovsample, linewidth = 0)
plt.ylabel('Coverage')
plt.xlabel('Position (bp)')
plt.title('Old coverage')
for ind, label in enumerate(plot_.get_xticklabels()):
    if ind % 50 == 0:  # every 5000th label is kept
        label.set_visible(True)
    else:
        label.set_visible(False)
plt.savefig('oldcovscatter.png')

plt.clf()
plt.cla()

plot_ = sns.barplot(x=xs, y = newcovsample, linewidth = 0)
plt.ylabel('Coverage')
plt.xlabel('Position (bp)')
plt.title('New coverage')
for ind, label in enumerate(plot_.get_xticklabels()):
    if ind % 50 == 0:  # every 5000th label is kept
        label.set_visible(True)
    else:
        label.set_visible(False)
plt.savefig('newcovscatter.png')

plt.clf()
plt.cla()

sns.distplot(oldcov, kde = False)
plt.ylabel('Number of bases')
plt.xlabel('Coverage')
plt.title('Old coverage')
plt.savefig('oldcovhist.png')

plt.clf()
plt.cla()

sns.distplot(newcov, kde = False)
plt.ylabel('Number of bases')
plt.xlabel('Coverage')
plt.title('New coverage')
plt.savefig('newcovhist.png')
