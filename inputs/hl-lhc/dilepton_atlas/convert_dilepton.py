## ATLAS

import ROOT
import numpy as np
import matplotlib.pyplot as plt

# Get experimental limit
infile = ROOT.TFile.Open("LimitInterpolator_CL95_14TeV.root","READ")
exp_limit = infile.Get("ExpLimit_ll")
infile.Close()
xvals_obs = []
yvals_obs = []
for i in range(exp_limit.GetN()) :
    xvals_obs.append(1000*exp_limit.GetPointX(i)) # convert to GeV
    yvals_obs.append(exp_limit.GetPointY(i)*1e-3) # convert to pb
xvals_obs = np.array(xvals_obs)
yvals_obs = np.array(yvals_obs)

# Get theory curves
infile_ee = ROOT.TFile.Open("DMCrossSectionGraphs_axial_ee.root","READ")
th_ee_dm1TeV = infile_ee.Get("graphFidXsec_1p0000")
infile_ee.Close()
infile_mumu = ROOT.TFile.Open("DMCrossSectionGraphs_axial_mumu.root","READ")
th_mumu_dm1TeV = infile_mumu.Get("graphFidXsec_1p0000")
infile_mumu.Close()
xvals_th = []
yvals_th = []
for i in range(th_ee_dm1TeV.GetN()) :
    x_ee = th_ee_dm1TeV.GetPointX(i)
    x_mumu = th_mumu_dm1TeV.GetPointX(i)
    if x_ee != x_mumu :
        print("Problem!!")
        print(x_ee, x_mumu)
        exit(1)
    y_tot = th_ee_dm1TeV.GetPointY(i) + th_mumu_dm1TeV.GetPointY(i)
    xvals_th.append(1000*x_ee) # convert to GeV
    yvals_th.append(y_tot*1e-3) # convert to pb
xvals_th = np.array(xvals_th)
yvals_th = np.array(yvals_th)

# Plot to validate
fig,ax=plt.subplots(1,1)
plt.xlim(2000, 12500)
plt.ylim(1e-7, 10)
plt.yscale('log')
plt.plot(xvals_obs,yvals_obs)
plt.plot(xvals_th, yvals_th)
plt.savefig('dilepton_hl-lhc_validation.pdf',bbox_inches='tight')

# And save to file
np.savez('hl-lhc-dilepton.npz',xvals_theory=xvals_th,yvals_theory=yvals_th,xvals_obs=xvals_obs,yvals_obs=yvals_obs)

# And a combined root file that matches
mygraph_obs = ROOT.TGraph()
for ix, iy in zip(xvals_obs, yvals_obs) :
    mygraph_obs.SetPoint(mygraph_obs.GetN(), ix, iy)
mygraph_th = ROOT.TGraph()
for ix, iy in zip(xvals_th,yvals_th) :
    mygraph_th.SetPoint(mygraph_th.GetN(),ix,iy)
outfile = ROOT.TFile.Open("hl-lhc-dilepton.root","RECREATE")
outfile.cd()
mygraph_obs.Write("expected_limit")
mygraph_th.Write("theory")
outfile.Close()

