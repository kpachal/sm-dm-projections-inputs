import numpy as np
import matplotlib.pyplot as plt

xvals = []
yvals = []
with open('Universal_Zprime_FCC_hh.txt') as f :
    found_table = False
    for line in f:
        if line.strip().startswith("(GeV)") :
            found_table = True
            continue
        elif not found_table :
            continue

        # Now extract numbers. Columns 0, 1 are mass, xsec * BR * ACC
        tokens = line.strip().split()
        xvals.append(float(tokens[0].strip(',')))
        yvals.append(float(tokens[1].strip(',')))
xvals_np = np.array(xvals)
yvals_np = np.array(yvals)

xvals_obs = []
yvals_obs = []
with open('discovery_limit_xsec_FCC_hh.txt') as f :
    found_table = False
    for line in f:
        if line.strip().startswith("(GeV)") :
            found_table = True
            continue
        elif not found_table :
            continue

        # Now extract numbers. Column 0 is mass and 2 is exclusion limit.
        tokens = line.strip().split()
        xvals_obs.append(float(tokens[0].strip(',')))
        yvals_obs.append(float(tokens[2].strip(',')))
xvals_obs_np = np.array(xvals_obs)
yvals_obs_np = np.array(yvals_obs)

# Check that these are what we expect by plotting them
fig,ax=plt.subplots(1,1)
plt.xlim(4500, 75000)
plt.ylim(1e-10, 10)
plt.yscale('log')
plt.plot(xvals_np,yvals_np)
plt.plot(xvals_obs_np, yvals_obs_np)
plt.savefig('dijet_fcc-hh_validation.pdf',bbox_inches='tight')

# And finally save a new output file
np.savez('fcc-hh-dijet.npz',xvals_theory=xvals_np,yvals_theory=yvals_np,xvals_obs=xvals_obs_np,yvals_obs=yvals_obs_np)
