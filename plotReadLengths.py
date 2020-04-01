import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys

allFn = sys.argv[1]
sampleFn = sys.argv[2]

allLengths = []
sampleLengths = []

with open(allFn) as f:
  for line in f:
    allLengths.append(int(line))
    
with open(sampleFn) as f:
  for line in f:
    sampleLengths.append(int(line))
    
sns.distplot(allLengths)
plt.title('All Read Lengths')
plt.savefig('oldreadlengths.png')

plt.clf()
plt.cla()

sns.distplot(sampleLengths)
plt.title('Sample Read Lengths')
plt.savefig('samplereadlengths.png')
