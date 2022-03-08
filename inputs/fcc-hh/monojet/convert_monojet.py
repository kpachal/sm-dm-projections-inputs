import ROOT
import numpy as np
import matplotlib.pyplot as plt


# Things to potentially change when I get more info from Phil
#ntuple = "phil_inputs/results_100_v13.root"
ntuple = "phil_inputs/results_100_100_VA_W_v2_update.root"

# Luminosity I want the equivalent of
target_luminosity = 3000 # ifb

# What is the equivalent in the files
luminosity = 100

# Things I know:
# obs is definitely nothing; I want exp
# xsec: all options look odd to me
# But I can reproduce phil's plots if I do

############
# Get tree and leaves
infile = ROOT.TFile.Open(ntuple,"READ")
tree = infile.Get("limit")
leaves = tree.GetListOfLeaves()

# define dynamically a python class containing root Leaves objects
class PyListOfLeaves(dict) :
    pass

# create an istance
pyl = PyListOfLeaves()

# Set up leaves
for i in range(0,leaves.GetEntries() ) :
    leaf = leaves.At(i)
    name = leaf.GetName()
    print("Getting leaf",name)
    # add dynamically attribute to my class 
    pyl.__setattr__(name,leaf)

# Now cycle through tree and interpret.
# This dict for storing outputs
limits = {"vector" : {}, "axial" : {}}
for entry in range(tree.GetEntries()):
    tree.GetEntry(entry)
    # get values from the tree using Python class pyl which contains leaves

    proc = pyl.fProc.GetValue()
    mdm = pyl.fMass.GetValue()
    mmed = pyl.fMed.GetValue()
    obs = pyl.fObs.GetValue()
    exp = pyl.fExp.GetValue()
    idn = pyl.fId.GetValue()
    gq = pyl.fGQ.GetValue()
    gdm = pyl.fGDM.GetValue()
    # think Phil said fXSY
    xs = pyl.fXS.GetValue() 

    # Want id==1 for merged 1 and 2j
    if not idn==1 : continue

    if not (proc == 800 or proc == 801) :
        continue
    elif proc==800 :
        processname = "vector"
    elif proc==801 :
        processname = "axial"

    if not (gq, gdm) in limits[processname].keys() :
        limits[processname][(gq, gdm)] = {"mmed" : [], "mdm" : [], "obs" : [], "exp" : [], "xs" : []}
    limits[processname][(gq, gdm)]["mmed"].append(mmed)
    limits[processname][(gq, gdm)]["mdm"].append(mdm)
    limits[processname][(gq, gdm)]["obs"].append(obs)
    limits[processname][(gq, gdm)]["exp"].append(exp)
    limits[processname][(gq, gdm)]["xs"].append(xs)
infile.Close()

# For cleaning out nan values
def clean_grid(xvals, yvals, zvals) :
  xclean = xvals[np.logical_not(np.logical_or(np.isnan(zvals),np.isinf(zvals)))]
  yclean = yvals[np.logical_not(np.logical_or(np.isnan(zvals),np.isinf(zvals)))]
  zclean = zvals[np.logical_not(np.logical_or(np.isnan(zvals),np.isinf(zvals)))]
  return xclean, yclean, zclean

# Now loop over our two processes.
for processname in ['vector','axial'] :

    zvals_obs = []
    zvals_exp = []

    # Check whether we ever have anything other than simple coupling limits in obs/exp dict
    if len(limits[processname].keys()) > 1 :
        print("Non simple limits - abort and do something else!")
        exit(1)
    original_couplings = list(limits[processname].keys())[0]

    datadict = limits[processname][original_couplings]

    # Good to go. Get observed and expected limits.
    xvals = datadict["mmed"]
    yvals = datadict["mdm"]
    obs_list = datadict["obs"]
    exp_list = datadict["exp"]
    xs_list = datadict["xs"]
    for obs, exp, xsec in zip(obs_list,exp_list,xs_list) :
        scalefactor = np.sqrt(float(luminosity)/float(target_luminosity))
        zvals_obs.append(obs*scalefactor/xsec)
        zvals_exp.append(exp*scalefactor/xsec)

    xvals_np = np.array(xvals)
    yvals_np = np.array(yvals)
    zvals_obs_np = np.array(zvals_obs)
    zvals_exp_np = np.array(zvals_exp)

    # Occasional nan: clean up
    xvals_np, yvals_np, zvals_exp_np = clean_grid(xvals_np, yvals_np, zvals_exp_np)

    # Make some validation plots
    for zname, zlist in zip(["exp"],[zvals_exp_np]) : # just do exp, obs is nothing
        levels = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5]#range(26)
        fig,ax=plt.subplots(1,1)
        xhigh = 12000.
        yhigh = 5000.
        plt.xlim(0, xhigh)
        plt.ylim(0, yhigh)
        ax.set_aspect(abs(xhigh/yhigh))
        cp = ax.tricontourf(xvals_np, yvals_np, zlist, levels=levels, cmap='Blues_r')
        fig.colorbar(cp)
        xexcl,yexcl = [],[]
        xnon,ynon = [],[]
        for x,y,z in zip(xvals_np,yvals_np,zlist) :
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
        ax.tricontour(xvals_np, yvals_np, zlist,levels=[1],colors=['y'],linewidths=[2])
        plt.savefig('monojet_fcc-hh_validation_{0}_{1}.pdf'.format(zname,processname),bbox_inches='tight')

        # And finally save a new output file
        np.savez('fcc-hh-monojet_{0}.npz'.format(processname),xvals=xvals_np,yvals=yvals_np,zvals=zlist)