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
    "test_gq" : np.logspace(np.log10(0.001),0,101),
    "test_gdm" : [0.0, 1.0],
    "test_gl" : [0.0],
  },
  "gdm_lim" : {
    "test_gq" : [0.01, 0.1, 0.25],
    "test_gdm" : np.logspace(np.log10(0.001),0,101),
    "test_gl" : [0.0]
  },
  "gl_lim" : {
    "test_gq" : [0.01, 0.1, 0.25],
    "test_gdm" : [0.0, 1.0],
    "test_gl" : np.logspace(np.log10(0.001),0,101),
  }
}

# Should be able to use this for everything
mMed_test = np.linspace(0,5000,201)
print("mMed_test:")
print(mMed_test)
test_mass_scenarios = {
  "dmDecoupled" : [100000 for i in mMed_test],
  "dmLight" : [1.0 for i in mMed_test],
  "DPLike" : [i/3.0 for i in mMed_test]
}

plotlims = {'hl-lhc' : 7500, 'fcc-hh' : 15000}

# Begin main code
####################

# For cleaning out nan values
def clean_grid(xvals, yvals, zvals) :
  xclean = xvals[np.logical_not(np.logical_or(np.isnan(zvals),np.isinf(zvals)))]
  yclean = yvals[np.logical_not(np.logical_or(np.isnan(zvals),np.isinf(zvals)))]
  zclean = zvals[np.logical_not(np.logical_or(np.isnan(zvals),np.isinf(zvals)))]
  return xclean, yclean, zclean

# Now process everything
for collider in ["hl-lhc","fcc-hh"] : # hl-lhc done already.

  print(collider)

  # Grids to use for limits that start 1D.
  # Big enough to stay valid for FCC
  target_mmed_vals = np.linspace(100,15000,150)
  target_mdm_vals = np.linspace(0,5000,101)
  target_mdm_vals = np.append(target_mdm_vals,100000) # Decoupled
  target_gvals = np.logspace(np.log10(0.01),0,101)
  target_xgrid, target_ygrid = np.meshgrid(target_mmed_vals,target_mdm_vals)

  # Start by creating an A2 scan and a rescaler. We'll use this
  # for all our limits that are 1D inputs.
  target_scan_A2 = DMAxialModelScan(mmed=target_xgrid.flatten(),mdm=target_ygrid.flatten(),
    gq=0.1, gdm=1.0, gl=0.1)  
  rescaler_fromA2 = Rescaler(target_scan_A2)

  # Monojet
  if "hl-lhc" in collider :
    monojet_data = np.load(input_path+'/{0}/monojet/{0}-monojet.npz'.format(collider))
  else :
    monojet_data = np.load(input_path+'/{0}/monojet/{0}-monojet_axial.npz'.format(collider))
  monojet_mmed = monojet_data['xvals']
  monojet_mdm = monojet_data['yvals']
  monojet_exdepth_A1 = monojet_data['zvals']
  # Already a grid, so can go straight to equivalent scan with this one.
  # Model and couplings taken from documentation in ATL-PHYS-PUB-2018-043:
  # it starts as A1, so we'll make our scan and rescaler from that.
  monojet_scan_A1 = DMAxialModelScan(mmed=monojet_mmed, mdm=monojet_mdm,
    gq=0.25, gdm=1.0, gl=0.0)
  rescaler_fromA1_monojetgrid = Rescaler(monojet_scan_A1)

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
  dijet_exdepth_A2 = dijet_xsec_limit.extract_exclusion_depths(target_scan_A2)

  # Dilepton, if doing
  if "hl-lhc" in collider :
      dilepton_atlas_data = np.load(input_path+'/{0}/dilepton_atlas/{0}-dilepton.npz'.format(collider))
      dilepton_atlas_mmed_th = dilepton_atlas_data['xvals_theory']
      dilepton_atlas_xsec_th = dilepton_atlas_data['yvals_theory']
      dilepton_atlas_mmed_obs = dilepton_atlas_data['xvals_obs']
      dilepton_atlas_xsec_obs = dilepton_atlas_data['yvals_obs']
      # Need to make a dictionary for the y limits from data,
      # because of width dependence. Let's just set it to 0.08
      dilepton_atlas_obs_lists = {0.08 : dilepton_atlas_xsec_obs}

      # Extract initial exclusion depth in an appropriate baseline scan
      # Parameters obtained from Jared: this is a 1 TeV DM particle
      # with couplings 0.1, 1, 0.1 and an axial model.
      dilepton_xsec_limit = CrossSectionLimit_Dilepton(
        mmed_limit=dilepton_atlas_mmed_obs,
        xsec_limit=dilepton_atlas_obs_lists,
        mmed_theory=dilepton_atlas_mmed_th,
        xsec_theory=dilepton_atlas_xsec_th,
        mdm=1000,
        gq=0.1,
        gdm=1.0,
        gl=0.1,
        coupling='axial'
      )
      dilepton_exdepth_A2 = dilepton_xsec_limit.extract_exclusion_depths(target_scan_A2)

  # With mono-x, two paths: already done for FCC, but needs computing for HL-LHC.
  # The computed version
  if "hl-lhc" in collider :  

    # Recall we only want to convert mono-x limits between models once, since it's slow.
    # So we'll go to V1 and then get other vector models from there.
    # Scans and rescaler we'll use in monojet
    V1_scan_monojetgrid = DMVectorModelScan(mmed=monojet_mmed, mdm=monojet_mdm,gq=0.25, gdm=1.0, gl=0.0)
    rescaler_fromV1_monojetgrid = Rescaler(V1_scan_monojetgrid)
    # And the actual scale factors: this is the slow bit
    monojet_sfs_A1toV1 = rescaler_fromA1_monojetgrid.rescale_by_hadronic_xsec_monox(0.25, 1.0, 0.0,'vector')[(0.25,1.0,0.0)]
    # TEST ONLY
    #monojet_sfs_A1toV1 = rescaler_fromA1_monojetgrid.rescale_by_parton_level_xsec_monox(0.25, 1.0, 0.0,'vector')[(0.25,1.0,0.0)]
  
  else :
    monojet_data_vec = np.load(input_path+'/{0}/monojet/{0}-monojet_vector.npz'.format(collider))
    monojet_mmed_vec = monojet_data['xvals']
    monojet_mdm_vec = monojet_data['yvals']
    monojet_exdepth_V1 = monojet_data['zvals']
    # Straight to rescaled
    V1_scan_monojetgrid = DMVectorModelScan(mmed=monojet_mmed_vec, mdm=monojet_mdm_vec,
      gq=0.25, gdm=1.0, gl=0.0)
    rescaler_fromV1_monojetgrid = Rescaler(V1_scan_monojetgrid)

  # Storage
  monojet_contours_axial = {}
  dijet_contours_axial = {}
  monojet_contours_vector = {}
  dijet_contours_vector = {}
  # Add dilepton matching all of the above for hl-lhc only
  dilepton_contours_axial = {}
  dilepton_contours_vector = {}
  dilepton_exclusiondepths_axial = {}
  dilepton_exclusiondepths_vector = {}

  # Now we're going to loop over our scenarios, going straight to contours and plots.
  for test_scenario in test_coupling_scenarios.keys() :

    test_gq = test_coupling_scenarios[test_scenario]["test_gq"]
    test_gdm = test_coupling_scenarios[test_scenario]["test_gdm"]
    test_gl = test_coupling_scenarios[test_scenario]["test_gl"]

    # Collect a range of interesting axial model scale factors for both of these, including A1 and A2
    # Naming keeps track of what we started from, since SFs will be different depending on the 
    # input scenario.
    dijet_sfs_allaxial_fromA2 = rescaler_fromA2.rescale_by_br_quarks(test_gq,test_gdm, test_gl,'axial')
    monojet_sfs_allaxial_fromA1 = rescaler_fromA1_monojetgrid.rescale_by_propagator(test_gq,test_gdm,test_gl,'axial')

    # Compute actual exclusion depths for axial
    monojet_exclusiondepths_axial = {k : monojet_exdepth_A1/v for k, v in monojet_sfs_allaxial_fromA1.items()}
    dijet_exclusiondepths_axial = {k : dijet_exdepth_A2/v for k, v in dijet_sfs_allaxial_fromA2.items()}

    # And some vector scans.
    dijet_sfs_allvector_fromA2 = rescaler_fromA2.rescale_by_br_quarks(test_gq,test_gdm,test_gl,'vector')
    monojet_sfs_allvector_fromV1 = rescaler_fromV1_monojetgrid.rescale_by_propagator(test_gq,test_gdm,test_gl,'vector')

    # Compute actual exclusion depths
    dijet_exclusiondepths_vector = {k : dijet_exdepth_A2/v for k, v in dijet_sfs_allvector_fromA2.items()}
    if "hl-lhc" in collider : 
      monojet_exclusiondepths_vector = {k : monojet_exdepth_A1/(monojet_sfs_A1toV1*v) for k, v in monojet_sfs_allvector_fromV1.items()}
    else :
      monojet_exclusiondepths_vector = {k : monojet_exdepth_V1/v for k, v in monojet_sfs_allvector_fromV1.items()}

    # Dilepton if doing
    if "hl-lhc" in collider :
      dilepton_sfs_allaxial_fromA2 = rescaler_fromA2.rescale_by_br_leptons(test_gq,test_gdm,test_gl,'axial')
      dilepton_sfs_allvector_fromA2 = rescaler_fromA2.rescale_by_br_leptons(test_gq,test_gdm,test_gl,'vector')
      dilepton_exclusiondepths_axial = {k : dilepton_exdepth_A2/v for k, v in dilepton_sfs_allaxial_fromA2.items()}
      dilepton_exclusiondepths_vector = {k : dilepton_exdepth_A2/v for k, v in dilepton_sfs_allvector_fromA2.items()}

    # For each test scenario, we actually have several sub-tests based on the grid of the non-scanned couplings.
    if "gq_lim" in test_scenario :
      test_couplings = test_gq
      coupling_index = 0
      other_first, other_second = np.meshgrid(test_gdm, test_gl)
      others_tag = "gdm{0}_gl{1}"
    elif "gdm_lim" in test_scenario :
      test_couplings = test_gdm
      coupling_index = 1
      other_first, other_second = np.meshgrid(test_gq, test_gl)
      others_tag = "gq{0}_gl{1}"
    else :
      test_couplings = test_gl
      coupling_index = 2
      other_first, other_second = np.meshgrid(test_gq, test_gdm)
      others_tag = "gq{0}_gdm{1}"
    for other_one, other_two in zip(other_first.flatten(), other_second.flatten()) :

      # Now need to extract contours. This will take interpolation: in addition to all the coupling scenarios
      # we have a few different slices we want to make through mMed, mDM space.
      for hypothesis in test_mass_scenarios.keys() :
        print("Beginning hypothesis",hypothesis)
        mDM_test = test_mass_scenarios[hypothesis]
        xvals = []
        yvals = []
        zvals_mono_axial = []
        zvals_dijet_axial = []
        zvals_dilepton_axial = []
        zvals_mono_vector = []
        zvals_dijet_vector = []
        zvals_dilepton_vector = []
        # Monojet needs mirroring because it doesn't actually go to zero. Double everything.
        use_xvals_mono = np.concatenate((V1_scan_monojetgrid.mmed,V1_scan_monojetgrid.mmed))
        use_yvals_mono = np.concatenate((V1_scan_monojetgrid.mdm,-1.0*V1_scan_monojetgrid.mdm))
        for coupling in test_couplings :
          full_couplings = [other_one, other_two]
          full_couplings.insert(coupling_index, coupling)
          
          # Universal across analyses because we are doing the same interpolations
          xvals += list(mMed_test)
          yvals += [coupling for i in mMed_test]      

          # Monojet extraction
          depths_mono_axial = monojet_exclusiondepths_axial[tuple(full_couplings)]
          zvals_mono_axial_raw = interpolate.griddata((use_xvals_mono, use_yvals_mono), np.concatenate((depths_mono_axial,depths_mono_axial)),(mMed_test,mDM_test),method='linear')
          zvals_mono_axial += list(zvals_mono_axial_raw)
          depths_mono_vector = monojet_exclusiondepths_vector[tuple(full_couplings)]
          zvals_mono_vector_raw = interpolate.griddata((use_xvals_mono, use_yvals_mono), np.concatenate((depths_mono_vector,depths_mono_vector)),(mMed_test,mDM_test),method='linear')
          zvals_mono_vector += list(zvals_mono_vector_raw)

          # Dijet extraction
          depths_dijet_axial = dijet_exclusiondepths_axial[tuple(full_couplings)]
          zvals_dijet_axial_raw = interpolate.griddata((target_scan_A2.mmed, target_scan_A2.mdm), depths_dijet_axial, (mMed_test, mDM_test),method='linear')
          zvals_dijet_axial += list(zvals_dijet_axial_raw)
          depths_dijet_vector = dijet_exclusiondepths_vector[tuple(full_couplings)]
          zvals_dijet_vector_raw = interpolate.griddata((target_scan_A2.mmed, target_scan_A2.mdm), depths_dijet_vector, (mMed_test, mDM_test),method='linear')
          zvals_dijet_vector += list(zvals_dijet_vector_raw)    

          # Dilepton extraction and storage
          if "hl-lhc" in collider :
            depths_dilepton_axial = dilepton_exclusiondepths_axial[tuple(full_couplings)]
            zvals_dilepton_axial_raw = interpolate.griddata((target_scan_A2.mmed, target_scan_A2.mdm), depths_dilepton_axial, (mMed_test, mDM_test),method='linear')
            zvals_dilepton_axial += list(zvals_dilepton_axial_raw)
            depths_dilepton_vector = dilepton_exclusiondepths_vector[tuple(full_couplings)]
            zvals_dilepton_vector_raw = interpolate.griddata((target_scan_A2.mmed, target_scan_A2.mdm), depths_dilepton_vector, (mMed_test, mDM_test),method='linear')
            zvals_dilepton_vector += list(zvals_dilepton_vector_raw)

        thiskey = "axial_{0}_{1}_".format(test_scenario,hypothesis)+others_tag.format(other_one, other_two)

        cleanx_mono, cleany_mono, cleanz_mono = clean_grid(np.array(xvals), np.array(yvals), np.array(zvals_mono_axial))
        if cleanz_mono.size > 0 :
          # Make a quick plot to check this looks sane
          drawContourPlotRough([[cleanx_mono, cleany_mono, cleanz_mono]], addPoints = False, this_tag = thiskey+"_monojet",plot_path = plot_path, xhigh=3000.,yhigh=0.5,vsCoupling=True)
          monojet_contours_axial[thiskey] = get_contours(cleanx_mono, cleany_mono, cleanz_mono)[0]

        cleanx_dijet, cleany_dijet, cleanz_dijet = clean_grid(np.array(xvals), np.array(yvals), np.array(zvals_dijet_axial))
        if cleanz_dijet.size > 0 :
          # Make a quick plot to check this looks sane
          drawContourPlotRough([[cleanx_dijet, cleany_dijet, cleanz_dijet]], addPoints = False, this_tag = thiskey+"_dijet",plot_path = plot_path, xhigh=3000.,yhigh=0.5,vsCoupling=True)
          dijet_contours_axial[thiskey] = get_contours(cleanx_dijet, cleany_dijet, cleanz_dijet)[0]

        cleanx_dilepton, cleany_dilepton, cleanz_dilepton = clean_grid(np.array(xvals), np.array(yvals), np.array(zvals_dilepton_axial))
        if cleanz_dilepton.size > 0 :
          # Make a quick plot to check this looks sane
          drawContourPlotRough([[cleanx_dilepton, cleany_dilepton, cleanz_dilepton]], addPoints = False, this_tag = thiskey+"_dilepton",plot_path = plot_path, xhigh=3000.,yhigh=0.5,vsCoupling=True)
          dilepton_contours_axial[thiskey] = get_contours(cleanx_dilepton, cleany_dilepton, cleanz_dilepton)[0]

        thiskey = "vector_{0}_{1}_".format(test_scenario,hypothesis)+others_tag.format(other_one, other_two)

        cleanx_mono, cleany_mono, cleanz_mono = clean_grid(np.array(xvals), np.array(yvals), np.array(zvals_mono_vector))
        if cleanz_mono.size > 0 :
          # Make a quick plot to check this looks sane
          drawContourPlotRough([[cleanx_mono, cleany_mono, cleanz_mono]], addPoints = False, this_tag = thiskey+"_monojet",plot_path = plot_path, xhigh=3000.,yhigh=0.5,vsCoupling=True)
          monojet_contours_vector[thiskey] = get_contours(cleanx_mono, cleany_mono, cleanz_mono)[0]

        cleanx_dijet, cleany_dijet, cleanz_dijet = clean_grid(np.array(xvals), np.array(yvals), np.array(zvals_dijet_vector))
        if cleanz_dijet.size > 0 :
          # Make a quick plot to check this looks sane
          drawContourPlotRough([[cleanx_dijet, cleany_dijet, cleanz_dijet]], addPoints = False, this_tag = thiskey+"_dijet",plot_path = plot_path, xhigh=3000.,yhigh=0.5,vsCoupling=True)
          dijet_contours_vector[thiskey] = get_contours(cleanx_dijet, cleany_dijet, cleanz_dijet)[0]

        cleanx_dilepton, cleany_dilepton, cleanz_dilepton = clean_grid(np.array(xvals), np.array(yvals), np.array(zvals_dilepton_vector))
        if cleanz_dilepton.size > 0 :
          # Make a quick plot to check this looks sane
          drawContourPlotRough([[cleanx_dilepton, cleany_dilepton, cleanz_dilepton]], addPoints = False, this_tag = thiskey+"_dilepton",plot_path = plot_path, xhigh=3000.,yhigh=0.5,vsCoupling=True)
          dilepton_contours_vector[thiskey] = get_contours(cleanx_dilepton, cleany_dilepton, cleanz_dilepton)[0]          

  # Save output in a clean way so that paper plot making script can be separate without re-running
  with open("vector_exclusion_depths_couplingmass_{0}.pkl".format(collider), "wb") as outfile_vec_depths :
    out_dict = {"dijet" : dijet_exclusiondepths_vector,
                "monojet" : monojet_exclusiondepths_vector,
                "dilepton" : dilepton_exclusiondepths_vector}
    pickle.dump(out_dict, outfile_vec_depths)
  with open("axial_exclusion_depths_couplingmass_{0}.pkl".format(collider), "wb") as outfile_axial_depths :
    out_dict = {"dijet" : dijet_exclusiondepths_axial,
                "monojet" : monojet_exclusiondepths_axial,
                "dilepton" : dilepton_exclusiondepths_axial}
    pickle.dump(out_dict, outfile_axial_depths)    
  with open("vector_exclusion_contours_couplingmass_{0}.pkl".format(collider), "wb") as poly_file:
    out_dict = {"dijet" : dijet_contours_vector,
                "monojet" : monojet_contours_vector,
                "dilepton" : dilepton_contours_vector}
    pickle.dump(out_dict, poly_file, pickle.HIGHEST_PROTOCOL)    
  with open("axial_exclusion_contours_couplingmass_{0}.pkl".format(collider), "wb") as poly_file:
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
    outfile = ROOT.TFile.Open("{0}_exclusion_contours_couplingmass_{1}.root".format(model,collider), "RECREATE")
    outfile.cd()    
    for analysis, smalldict in middict.items() :
      for key, contour in smalldict.items() :
        igraph = ROOT.TGraph()
        for icontour in contour :
          for x, y in list(icontour.exterior.coords) :
            igraph.AddPoint(x,y)
        outname = "{0}_{1}".format(key,analysis)
        igraph.Write(outname)
  outfile.Close()