import sys,os
import numpy as np
import pickle
from scipy import interpolate
import ROOT

from couplingscan.scan import *
from couplingscan.rescaler import *
from couplingscan.limitparsers import *
from basic_plotter import *

# Adjust to wherever you put the inputs
input_path = "../inputs/"

# Set up plot location for these ones
plot_path = "plots/validation"

# Couplings to test: 
# turns out "go big or go home" does actually hit a limit here. So we're going to do a few versions of this.
# based on explicitly what we are trying to target
test_coupling_scenarios = {
  "gq_lim" : {
    "test_gq" : np.logspace(np.log10(0.0001),0,101),
    "test_gdm" : [0.0, 0.05, 0.1, 0.2, 1.0],
    "test_gl" : [0.0],
  },
  "gdm_lim" : {
    "test_gq" : [0.01, 0.02, 0.05, 0.1, 0.25],
    "test_gdm" : np.logspace(np.log10(0.0001),0,101),
    "test_gl" : [0.0]
  },
  "gl_lim" : {
    "test_gq" : [0.01, 0.1, 0.25],
    "test_gdm" : [0.0, 1.0],
    "test_gl" : np.logspace(np.log10(0.0001),0,101),
  }
}

# Begin main code
####################

# For cleaning out nan values
def clean_grid(xvals, yvals, zvals) :
  if zvals.size == 0 :
    return np.array([]), np.array([]), np.array([])
  xclean = xvals[np.logical_and(np.logical_and(np.isfinite(xvals),np.isfinite(yvals)),np.isfinite(zvals))]
  yclean = yvals[np.logical_and(np.logical_and(np.isfinite(xvals),np.isfinite(yvals)),np.isfinite(zvals))]
  zclean = zvals[np.logical_and(np.logical_and(np.isfinite(xvals),np.isfinite(yvals)),np.isfinite(zvals))]
  return xclean, yclean, zclean

plotlims = {"mDM" : {'hl-lhc' : 2200, 'fcc-hh' : 2200},
            "MMed" : {'hl-lhc' : 7500, 'fcc-hh' : 15000}}

# Now process everything
for xaxis in ["mDM", "MMed"] :

  # Axis ranges and scanned ratios vary with which
  # item is being scanned.
  xmax = 5000 if "mDM" in xaxis else 15000
  m_x_test = np.linspace(0,xmax,1001)
  m_x_test[0] = 1.0

  test_mass_scenarios = { 
    # x axis is mDM: MMed should be a multiple of it (or fixed)
    "mDM" :
    {
      "ratio2p5_fixedMMed" : [i*2.5 for i in m_x_test],
      "ratio3_fixedMMed" : [i*3.0 for i in m_x_test],
      "ratio5_fixedMMed" : [i*5.0 for i in m_x_test],
      "ratio10_fixedMMed" : [i*10.0 for i in m_x_test],
      "ratio100_fixedMMed" : [i*100.0 for i in m_x_test],
      "ratio1000_fixedMMed" : [i*1000.0 for i in m_x_test],
    #  "mMed50GeV_fixedMMed" : [50.0 for i in mDM_test],
    #  "mMed1000GeV_fixedMMed" : [1000.0 for i in mDM_test]      
    },
    # x axis is MMed: mDM should be a fraction of it (or fixed, incl. decoupled)
    "MMed" :
    {
      "dmDecoupled" : [100000 for i in m_x_test],
      "dm1GeV" : [1.0 for i in m_x_test],
      "DPLike" : [i/3.0 for i in m_x_test],
      "dm10MeV" : [0.01 for i in m_x_test],
      "dm100MeV" : [0.1 for i in m_x_test],
      "dm10GeV" : [10.0 for i in m_x_test],
      "dm100GeV" : [100.0 for i in m_x_test]
    }    
  }

  for collider in ["hl-lhc","fcc-hh"] : # hl-lhc done already.

    print(collider)

    # Monojet
    if "hl-lhc" in collider :
      monojet_data = np.load(input_path+'/{0}/monojet/{0}-monojet.npz'.format(collider))
    else :
      monojet_data = np.load(input_path+'/{0}/monojet/{0}-monojet_axial.npz'.format(collider))
    monojet_mmed_ax = monojet_data['xvals']
    monojet_mdm_ax = monojet_data['yvals']
    monojet_exdepth_A1 = monojet_data['zvals']
    # Already a grid, so can go straight to equivalent scan with this one.
    # Make a rescaler so we can get the V in one go.
    monojet_scan_A1 = DMAxialModelScan(mmed=monojet_mmed_ax, mdm=monojet_mdm_ax,
      gq=0.25, gdm=1.0, gl=0.0)
    rescaler_fromA1_monojetgrid = Rescaler(monojet_scan_A1, monojet_exdepth_A1, 0.5)

    # Dijet
    dijet_data = np.load(input_path+'/{0}/dijet/{0}-dijet.npz'.format(collider))
    dijet_mmed_th = dijet_data['xvals_theory']
    dijet_xsec_th = dijet_data['yvals_theory']
    dijet_mmed_obs = dijet_data['xvals_obs']
    dijet_xsec_obs = dijet_data['yvals_obs']
    # Need to extract initial exclusion depth in an appropriate baseline scan
    dijet_xsec_limit = CrossSectionLimit_Dijet(
      mmed_limit=dijet_mmed_obs,
      xsec_limit=dijet_xsec_obs,
      mmed_theory=dijet_mmed_th,
      xsec_theory=dijet_xsec_th,
      mdm=10000,
      gq=0.25,
      gdm=0.0,
      gl=0.0,
      coupling='vector',
      max_intrinsic_width=0.12 # As specified by Robert
    )

    # Dilepton, if doing
    if "hl-lhc" in collider :
        dilepton_data = np.load(input_path+'/{0}/dilepton_atlas/{0}-dilepton.npz'.format(collider))
        dilepton_mmed_th = dilepton_data['xvals_theory']
        dilepton_xsec_th = dilepton_data['yvals_theory']
        dilepton_mmed_obs = dilepton_data['xvals_obs']
        dilepton_xsec_obs = dilepton_data['yvals_obs']
        # Need to make a dictionary for the y limits from data,
        # because of width dependence. Let's just set it to 0.08
        dilepton_obs_lists = {0.08 : dilepton_xsec_obs}

        # Extract initial exclusion depth in an appropriate baseline scan
        # Parameters obtained from Jared: this is a 1 TeV DM particle
        # with couplings 0.1, 1, 0.1 and an axial model.
        dilepton_xsec_limit = CrossSectionLimit_Dilepton(
          mmed_limit=dilepton_mmed_obs,
          xsec_limit=dilepton_obs_lists,
          mmed_theory=dilepton_mmed_th,
          xsec_theory=dilepton_xsec_th,
          mdm=1000,
          gq=0.1,
          gdm=1.0,
          gl=0.1,
          coupling='axial'
        )

    # With mono-x, two paths: already done for FCC, but needs computing for HL-LHC.
    # The computed version
    if "hl-lhc" in collider :  

      # Recall we only want to convert mono-x limits between models once, since it's slow.
      # So we'll go to V1 and then get other vector models from there.
      # Scans and rescaler we'll use in monojet
      monojet_mmed_vec = monojet_mmed_ax
      monojet_mdm_vec = monojet_mdm_ax
      V1_scan_monojetgrid = DMVectorModelScan(mmed=monojet_mmed_vec, mdm=monojet_mdm_vec,gq=0.25, gdm=1.0, gl=0.0)
      monojet_exdepth_V1 = rescaler_fromA1_monojetgrid.rescale_by_hadronic_xsec_monox(0.25, 1.0, 0.0,'vector')[(0.25,1.0,0.0)]
      # TEST ONLY
      #monojet_exdepth_V1 = rescaler_fromA1_monojetgrid.rescale_by_parton_level_xsec_monox(0.25, 1.0, 0.0,'vector')[(0.25,1.0,0.0)]
    
    else :
      monojet_data_vec = np.load(input_path+'/{0}/monojet/{0}-monojet_vector.npz'.format(collider))
      monojet_mmed_vec = monojet_data_vec['xvals']
      monojet_mdm_vec = monojet_data_vec['yvals']
      monojet_exdepth_V1 = monojet_data_vec['zvals']
    
    # Mirror monojet here, since it doesn't actually go all the way to the axis.
    # New
    use_mmed_mono = np.concatenate((monojet_mmed_vec,monojet_mmed_vec))
    use_mdm_mono = np.concatenate((monojet_mdm_vec,-1.0*monojet_mdm_vec))
    use_exdepth_mono = np.concatenate((monojet_exdepth_V1,monojet_exdepth_V1))
    monojet_mmed_vec = use_mmed_mono
    monojet_mdm_vec = use_mdm_mono
    monojet_exdepth_V1 = use_exdepth_mono

    # Storage
    monojet_contours_axial = {}
    dijet_contours_axial = {}
    monojet_contours_vector = {}
    dijet_contours_vector = {}
    # Add dilepton matching all of the above for hl-lhc only
    dilepton_contours_axial = {}
    dilepton_contours_vector = {}

    # Now we're going to loop over our scenarios, going straight to contours and plots.
    for hypothesis, test_masses in test_mass_scenarios[xaxis].items() :

      # For dijet and dilepton.
      if "mDM" in xaxis :
        line_scan_A2 = DMAxialModelScan(
          mmed=test_masses,
          mdm=m_x_test,
          gq=0.1, gdm=1.0, gl=0.01)
      else :
        line_scan_A2 = DMAxialModelScan(
          mmed=m_x_test,
          mdm=test_masses,
          gq=0.1, gdm=1.0, gl=0.01)
      initial_depths_dijet = dijet_xsec_limit.extract_exclusion_depths(line_scan_A2)
      rescaler_fromA2_dijet = Rescaler(line_scan_A2,initial_depths_dijet,0.12)
      if "hl-lhc" in collider :
        initial_depths_dilepton = dilepton_xsec_limit.extract_exclusion_depths(line_scan_A2)
      else :
        if "mDM" in xaxis : initial_depths_dilepton = np.array([np.nan for i in test_masses])
        else : initial_depths_dilepton = np.array([np.nan for i in m_x_test])
      rescaler_fromA2_dilepton = Rescaler(line_scan_A2,initial_depths_dilepton,0.1)

      # For monojet: interpolate from existing limit grid (separate A and V)
      (use_monojetx, use_monojety) = (test_masses, m_x_test) if "mDM" in xaxis else (m_x_test, test_masses)
      line_scan_A1 = DMAxialModelScan(mmed=use_monojetx,mdm=use_monojety,gq=0.25,gdm=1.0,gl=0.0)
      initial_depths_axial_monojet = interpolate.griddata((monojet_mmed_ax,monojet_mdm_ax), monojet_exdepth_A1, (use_monojetx,use_monojety))
      rescaler_fromA1_monojet = Rescaler(line_scan_A1, initial_depths_axial_monojet,0.4)
      line_scan_V1 = DMVectorModelScan(mmed=use_monojetx,mdm=use_monojety,gq=0.25,gdm=1.0,gl=0.0)
      initial_depths_vector_monojet = interpolate.griddata((monojet_mmed_vec,monojet_mdm_vec), monojet_exdepth_V1, (use_monojetx,use_monojety))
      rescaler_fromV1_monojet = Rescaler(line_scan_V1, initial_depths_vector_monojet,0.4)

      for test_coupling, coupling_dict in test_coupling_scenarios.items() :
        test_gq = coupling_dict["test_gq"]
        test_gdm = coupling_dict["test_gdm"]
        test_gl = coupling_dict["test_gl"]

        # Loop through non-primary couplings and pull out plots
        if "gq_lim" in test_coupling :
          test_couplings = test_gq
          coupling_index = 0
          other_first, other_second = np.meshgrid(test_gdm, test_gl)
          others_tag = "gdm{0}_gl{1}"
        elif "gdm_lim" in test_coupling :
          test_couplings = test_gdm
          coupling_index = 1
          other_first, other_second = np.meshgrid(test_gq, test_gl)
          others_tag = "gq{0}_gl{1}"
        else :
          test_couplings = test_gl
          coupling_index = 2
          other_first, other_second = np.meshgrid(test_gq, test_gdm)
          others_tag = "gq{0}_gdm{1}"

        limitplot_x, limitplot_y = np.meshgrid(m_x_test,test_couplings)
        full_limits_dijet_axial = rescaler_fromA2_dijet.rescale_by_br_quarks(test_gq, test_gdm, test_gl,'axial')
        full_limits_dijet_vector = rescaler_fromA2_dijet.rescale_by_br_quarks(test_gq, test_gdm, test_gl,'vector')
        full_limits_dilepton_axial = rescaler_fromA2_dilepton.rescale_by_br_leptons(test_gq, test_gdm, test_gl,'axial')
        full_limits_dilepton_vector = rescaler_fromA2_dilepton.rescale_by_br_leptons(test_gq, test_gdm, test_gl,'vector')
        full_limits_monojet_axial = rescaler_fromA1_monojet.rescale_by_propagator(test_gq, test_gdm, test_gl,'axial')
        full_limits_monojet_vector = rescaler_fromV1_monojet.rescale_by_propagator(test_gq, test_gdm, test_gl,'vector')

        # Now we're looping over the other two couplings
        for other_one, other_two in zip(other_first.flatten(), other_second.flatten()) :
          these_depths_dij_ax = []
          these_depths_dij_vec = []
          these_depths_mono_ax = []
          these_depths_mono_vec = []
          these_depths_dil_ax = []
          these_depths_dil_vec = []

          for ci in test_couplings :
            full_couplings = [other_one, other_two]
            full_couplings.insert(coupling_index, ci)

            these_depths_dij_ax.append(full_limits_dijet_axial[tuple(full_couplings)])
            these_depths_dij_vec.append(full_limits_dijet_vector[tuple(full_couplings)])
            these_depths_mono_ax.append(full_limits_monojet_axial[tuple(full_couplings)])
            these_depths_mono_vec.append(full_limits_monojet_vector[tuple(full_couplings)])
            these_depths_dil_ax.append(full_limits_dilepton_axial[tuple(full_couplings)])
            these_depths_dil_vec.append(full_limits_dilepton_vector[tuple(full_couplings)])          

          # Plots with contours, and extract contours for storage
          thiskey = collider+"_axial_{0}_{1}_".format(test_coupling,hypothesis)+others_tag.format(other_one, other_two)

          xdij_ax, ydij_ax, zdij_ax = clean_grid(limitplot_x.flatten(), limitplot_y.flatten(), np.array(these_depths_dij_ax).flatten())
          drawContourPlotRough([[xdij_ax, ydij_ax, zdij_ax]], addPoints = False, this_tag = thiskey+"_dijet",plot_path = plot_path, xhigh=plotlims[xaxis][collider],yhigh=0.5,vsCoupling=True)
          dijet_contours_axial[thiskey] = get_contours(xdij_ax, ydij_ax, zdij_ax)[0]

          xmono_ax, ymono_ax, zmono_ax = clean_grid(limitplot_x.flatten(), limitplot_y.flatten(), np.array(these_depths_mono_ax).flatten())
          drawContourPlotRough([[xmono_ax, ymono_ax, zmono_ax]], addPoints = False, this_tag = thiskey+"_monojet",plot_path = plot_path, xhigh=plotlims[xaxis][collider],yhigh=0.5,vsCoupling=True)
          monojet_contours_axial[thiskey] = get_contours(xmono_ax, ymono_ax, zmono_ax)[0]

          xdil_ax, ydil_ax, zdil_ax = clean_grid(limitplot_x.flatten(), limitplot_y.flatten(), np.array(these_depths_dil_ax).flatten())
          drawContourPlotRough([[xdil_ax, ydil_ax, zdil_ax]], addPoints = False, this_tag = thiskey+"_dilepton",plot_path = plot_path, xhigh=plotlims[xaxis][collider],yhigh=0.5,vsCoupling=True)
          dilepton_contours_axial[thiskey] = get_contours(xdil_ax, ydil_ax, zdil_ax)[0]

          thiskey = collider+"_vector_{0}_{1}_".format(test_coupling,hypothesis)+others_tag.format(other_one, other_two)

          xdij_vec, ydij_vec, zdij_vec = clean_grid(limitplot_x.flatten(), limitplot_y.flatten(), np.array(these_depths_dij_vec).flatten())
          drawContourPlotRough([[xdij_vec, ydij_vec, zdij_vec]], addPoints = False, this_tag = thiskey+"_dijet",plot_path = plot_path, xhigh=plotlims[xaxis][collider],yhigh=0.5,vsCoupling=True)
          dijet_contours_vector[thiskey] = get_contours(xdij_vec, ydij_vec, zdij_vec)[0]

          xmono_vec, ymono_vec, zmono_vec = clean_grid(limitplot_x.flatten(), limitplot_y.flatten(), np.array(these_depths_mono_vec).flatten())
          drawContourPlotRough([[xmono_vec, ymono_vec, zmono_vec]], addPoints = False, this_tag = thiskey+"_monojet",plot_path = plot_path, xhigh=plotlims[xaxis][collider],yhigh=0.5,vsCoupling=True)
          monojet_contours_vector[thiskey] = get_contours(xmono_vec, ymono_vec, zmono_vec)[0]

          xdil_vec, ydil_vec, zdil_vec = clean_grid(limitplot_x.flatten(), limitplot_y.flatten(), np.array(these_depths_dil_vec).flatten())
          drawContourPlotRough([[xdil_vec, ydil_vec, zdil_vec]], addPoints = False, this_tag = thiskey+"_dilepton",plot_path = plot_path, xhigh=plotlims[xaxis][collider],yhigh=0.5,vsCoupling=True)
          dilepton_contours_vector[thiskey] = get_contours(xdil_vec, ydil_vec, zdil_vec)[0]

    # Save output in a clean way so that paper plot making script can be separate without re-running
    with open("vector_exclusion_contours_coupling{0}_{1}.pkl".format(xaxis,collider), "wb") as poly_file:
      out_dict = {"dijet" : dijet_contours_vector,
                  "monojet" : monojet_contours_vector,
                  "dilepton" : dilepton_contours_vector}
      pickle.dump(out_dict, poly_file, pickle.HIGHEST_PROTOCOL)    
    with open("axial_exclusion_contours_coupling{0}_{1}.pkl".format(xaxis,collider), "wb") as poly_file:
      out_dict = {"dijet" : dijet_contours_axial,
                  "monojet" : monojet_contours_axial,
                  "dilepton" : dilepton_contours_axial}
      pickle.dump(out_dict, poly_file, pickle.HIGHEST_PROTOCOL)

    # And some TGraphs for ease
    big_ol_dict = {
      "vector" : {"dijet" : dijet_contours_vector,
                  "monojet" : monojet_contours_vector,
                  "dilepton" : dilepton_contours_vector},
      "axial" : {"dijet" : dijet_contours_axial,
                  "monojet" : monojet_contours_axial,
                  "dilepton" : dilepton_contours_axial}
    }
    for model, middict in big_ol_dict.items() :
      outfilename = "{0}_exclusion_contours_coupling{1}_{2}.root".format(model,xaxis,collider)
      outfile = ROOT.TFile.Open(outfilename, "RECREATE")
      outfile.cd()    
      for analysis, smalldict in middict.items() :
        for key, contour in smalldict.items() :
          igraph = ROOT.TGraph()
          if not contour : continue
          for icontour in contour :
            for x, y in list(icontour.exterior.coords) :
              igraph.AddPoint(x,y)
          outname = "{0}_{1}".format(key,analysis)
          igraph.Write(outname)
    outfile.Close()