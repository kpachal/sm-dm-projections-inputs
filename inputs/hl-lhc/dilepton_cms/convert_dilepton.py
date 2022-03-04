import numpy as np
import matplotlib.pyplot as plt
import statistics as stats

points_full = {}
with open('ZPrime_limitCard_HLLHC_Exp.txt') as f :
    for line in f :
        tokens = line.strip().split()
        mass = tokens[0].strip(',')
        if not mass in points_full.keys() :
            points_full[mass] = []
        points_full[mass].append(float(tokens[1].strip(',')))
xvals = []
yvals = []
for xval in points_full.keys() :
    xvals.append(float(xval))
    yval = stats.median(points_full[xval])
    yvals.append(yval*1928*1.076923)

# Plot it
print(xvals)
print(yvals)
fig,ax=plt.subplots(1,1)
plt.xlim(1250, 8000)
#plt.ylim(1e-7, 10)
plt.yscale('log')
plt.plot(xvals,yvals)
plt.savefig('dilepton_cms_hl-lhc_validation.pdf',bbox_inches='tight')