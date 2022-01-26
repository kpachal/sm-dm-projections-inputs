import numpy as np
import matplotlib.pyplot as plt

# Get from initial file
xvals = []
yvals = []
zvals = []
with open('limits_3000ifb.txt') as f:
    found_table = False
    for line in f:
        if line.startswith("med_mass") :
            found_table = True
        if not found_table or line.startswith("med_mass") :
            continue

        # Now extract numbers. Columns 0, 1, 3 are x, y z
        tokens = line.strip().split()
        xvals.append(float(tokens[0].strip(',')))
        yvals.append(float(tokens[1].strip(',')))
        zvals.append(float(tokens[3].strip(',')))
xvals_np = np.array(xvals)
yvals_np = np.array(yvals)
zvals_np = np.array(zvals)

# Check that these are what we expect by plotting them
levels = range(26)
fig,ax=plt.subplots(1,1)
plt.xlim(0, 3000)
plt.ylim(0, 1000)
ax.set_aspect(abs(3000./1000.))
cp = ax.tricontourf(xvals_np, yvals_np, zvals_np, levels=levels, cmap='Blues_r')
fig.colorbar(cp)
xexcl,yexcl = [],[]
xnon,ynon = [],[]
for x,y,z in zip(xvals_np,yvals_np,zvals_np) :
    if z < 1. : 
        xexcl.append(x)
        yexcl.append(y)
    else :
        xnon.append(x)
        ynon.append(y)
ax.scatter(xnon,ynon,color='red', marker='o',facecolors='none',linewidths=2)
ax.scatter(xexcl,yexcl,color='white', marker='o',facecolors='none',linewidths=2)
ax.set_xlabel("m$_{ZA}$ [GeV]")
ax.set_ylabel("m$_{\chi}$ [GeV]")
# Add my own contours
ax.tricontour(xvals, yvals, zvals,levels=[1],colors=['w'],linewidths=[2])
plt.savefig('monojet_hl-lhc_validation.pdf',bbox_inches='tight')

# And finally save a new output file
np.savez('hl-lhc-monojet.npz',xvals_np,yvals_np,zvals_np)
