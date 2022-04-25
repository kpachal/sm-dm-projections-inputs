import numpy as np
import matplotlib.pyplot as plt
import statistics as stats
import ROOT

# Input: units = picobarns I think.

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

# Theory curve: take one with no acceptance factor from Jared and convert to pb.
# Will need to extract and convert. Use the 1 so interpretation matches ATLAS.
infile_ee = ROOT.TFile.Open("DMCrossSectionGraphs_axial_ee.root","READ")
th_ee_dm1TeV = infile_ee.Get("graphFidXsec_1p0000")
infile_ee.Close()
infile_mumu = ROOT.TFile.Open("DMCrossSectionGraphs_axial_mumu.root","READ")
th_mumu_dm1TeV = infile_mumu.Get("graphFidXsec_1p0000")
th_mumu_dm1TeV_fullxsec = infile_mumu.Get("graphXsec_1p0000")
infile_mumu.Close()

acceptance = ROOT.TGraph()
th_ee_dm1TeV_fullxsec = ROOT.TGraph()
th_tot_fullxsec = ROOT.TGraph()
xvals_th = []
yvals_th = []
for point in range(th_ee_dm1TeV.GetN()) :
    xval = th_ee_dm1TeV.GetPointX(point)
    yval_ee = th_ee_dm1TeV.GetPointY(point)*1e-3 # convert to pb
    th_mumu_full = th_mumu_dm1TeV_fullxsec.Eval(xval)*1e-3 # convert to pb
    th_mumu_fiducial = th_mumu_dm1TeV.Eval(xval)*1e-3 # convert to pb
    acc = th_mumu_fiducial/th_mumu_full
    xval_gev = xval*1000. # convert to GeV
    acceptance.SetPoint(point,xval_gev,acc)
    xsec_full = th_mumu_full+(yval_ee/acc)
    th_ee_dm1TeV_fullxsec.SetPoint(point,xval_gev,yval_ee/acc)
    th_tot_fullxsec.SetPoint(point,xval_gev,xsec_full)
    xvals_th.append(xval_gev)
    yvals_th.append(xsec_full)
xvals_th = np.array(xvals_th)
yvals_th = np.array(yvals_th)

# Plot to validate
print(xvals)
print(yvals)
fig,ax=plt.subplots(1,1)
plt.xlim(1250, 8000)
#plt.ylim(1e-7, 10)
plt.yscale('log')
plt.plot(xvals,yvals)
plt.plot(xvals_th,yvals_th)
plt.savefig('dilepton_cms_hl-lhc_validation.pdf',bbox_inches='tight')

# And save to file
np.savez('hl-lhc-dilepton-cms.npz',xvals_theory=xvals_th,yvals_theory=yvals_th,xvals_obs=xvals,yvals_obs=yvals)

# Output root file
mygraph_obs = ROOT.TGraph()
for ix, iy in zip(xvals, yvals) :
    mygraph_obs.SetPoint(mygraph_obs.GetN(),ix,iy) # To same units
outfile = ROOT.TFile.Open("hl-lhc-dilepton-cms.root","RECREATE")
outfile.cd()
acceptance.Write("acceptance")
th_ee_dm1TeV_fullxsec.Write("graphXsec_1p0000")
th_tot_fullxsec.Write("total_xsec_nb")
mygraph_obs.Write("expected_limit_nb")
outfile.Close()
