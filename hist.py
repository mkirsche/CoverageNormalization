import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys


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
