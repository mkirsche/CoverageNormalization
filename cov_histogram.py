#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys

if (len(sys.argv) != 2):
  print("cov_histogram.py oldcov.newcov.txt")
  exit(0)


fn = sys.argv[1]
with open(fn) as f:
  oldcov = []
  newcov = []
  for line in f:
    tokens = line.split()
    oldcov.append(int(tokens[0]))
    newcov.append(int(tokens[1]))
    
sns.distplot(oldcov)
plt.title('Old coverage')
plt.savefig('oldcov.png')

plt.clf()
plt.cla()

sns.distplot(newcov)
plt.title('New coverage')
plt.savefig('newcov.png')
