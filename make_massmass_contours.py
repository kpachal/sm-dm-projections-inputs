import sys,os
import numpy as np
import pickle

from couplingscan.scan import *
from couplingscan.rescaler import *
from couplingscan.limitparsers import *
from basic_plotter import *

# Adjust to wherever you put the inputs
input_path = "/Users/katherinepachal/Code/sm-dm-projections-inputs/"

# Set up plot location for these ones
plot_path = "plots/validation"

# Couplings to test: you'll get a grid of all combinations
test_gq = [0.25,0.15,0.1,0.01]
test_gdm = 1.0
test_gl = [0.1,0.05,0.01,0.0]

# Begin main code
####################

# For cleaning out nan values
def clean_grid(xvals, yvals, zvals) :
  xclean = xvals[np.logical_not(np.logical_or(np.isnan(zvals),np.isinf(zvals)))]
  yclean = yvals[np.logical_not(np.logical_or(np.isnan(zvals),np.isinf(zvals)))]
  zclean = zvals[np.logical_not(np.logical_or(np.isnan(zvals),np.isinf(zvals)))]
  return xclean, yclean, zclean

# Now process everything
for collider in ["hl-lhc","fcc-hh"] :

  print(collider)

  # Grids to use for limits that start 1D
  target_xvals = np.linspace(100,10000,100)
  target_yvals = np.linspace(0,1200,25)
  target_xgrid, target_ygrid = np.meshgrid(target_xvals,target_yvals)

  # Start by creating an A2 scan and a rescaler. We'll use this
  # for all our limits that are 1D inputs.
  target_scan_A2 = DMAxialModelScan(mmed=target_xgrid.flatten(),mdm=target_ygrid.flatten(),
    gq=0.1, gdm=1.0, gl=0.1)  
  rescaler_fromA2 = Rescaler(target_scan_A2)

  # Monojet
  monojet_data = np.load(input_path+'/{0}/monojet/{0}-monojet.npz'.format(collider))
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

  # Collect a range of interesting axial model scale factors for both of these, including A1 and A2
  # Naming keeps track of what we started from, since SFs will be different depending on the 
  # input scenario.
  dijet_sfs_allaxial_fromA2 = rescaler_fromA2.rescale_by_br_quarks(test_gq,test_gdm, test_gl,'axial')
  monojet_sfs_allaxial_fromA1 = rescaler_fromA1_monojetgrid.rescale_by_propagator(test_gq,test_gdm,test_gl,'axial')

  # Compute actual exclusion depths for axial
  monojet_exclusiondepths_axial = {k : monojet_exdepth_A1/v for k, v in monojet_sfs_allaxial_fromA1.items()}
  dijet_exclusiondepths_axial = {k : dijet_exdepth_A2/v for k, v in dijet_sfs_allaxial_fromA2.items()}

  # Extract contours
  monojet_contours_axial = {}
  dijet_contours_axial = {}
  for coupling in monojet_exclusiondepths_axial.keys() :
    monojet_depth = monojet_exclusiondepths_axial[coupling]
    monojet_contours_axial[coupling] = get_contours(monojet_mmed, monojet_mdm, monojet_depth)[0]
    dijet_depth = dijet_exclusiondepths_axial[coupling]
    dijet_contours_axial[coupling] = get_contours(target_scan_A2.mmed, target_scan_A2.mdm, dijet_depth)[0]

  # Now let's do some vector scans.
  # Recall we only want to convert mono-x limits between models once, since it's slow.
  # So we'll go to V1 and then get other vector models from there.
  # Scans and rescaler we'll use in monojet
  V1_scan_monojetgrid = DMVectorModelScan(mmed=monojet_mmed, mdm=monojet_mdm,gq=0.25, gdm=1.0, gl=0.0)
  rescaler_fromV1_monojetgrid = Rescaler(V1_scan_monojetgrid)
  # And the actual scale factors: this is the slow bit
  monojet_sfs_A1toV1 = rescaler_fromV1_monojetgrid.rescale_by_hadronic_xsec_monox(0.25, 1.0, 0.0,'vector')[(0.25,1.0,0.0)]
  # For fast tests - NOT REAL, DON'T KEEP
  #monojet_sfs_A1toV1 = rescaler_fromV1_monojetgrid.rescale_by_parton_level_xsec_monox(0.25, 1.0, 0.0,'vector')[(0.25,1.0,0.0)]

  # Debug: this should match one of the plots that comes later.
  debug_monojetv1 = monojet_exdepth_A1/monojet_sfs_A1toV1
  drawContourPlotRough([[monojet_mmed, monojet_mdm, debug_monojetv1]], addPoints = False, this_tag = "monojet_V1",plot_path = plot_path)

  # And now we can collect a range of interesting vector model scale factors, including V1 and V2.
  dijet_sfs_allvector_fromA2 = rescaler_fromA2.rescale_by_br_quarks(test_gq,test_gdm,test_gl,'vector')
  monojet_sfs_allvector_fromV1 = rescaler_fromV1_monojetgrid.rescale_by_propagator(test_gq,test_gdm,test_gl,'vector')

  # Compute actual exclusion depths
  monojet_exclusiondepths_vector = {k : monojet_exdepth_A1/(monojet_sfs_A1toV1*v) for k, v in monojet_sfs_allvector_fromV1.items()}
  dijet_exclusiondepths_vector = {k : dijet_exdepth_A2/v for k, v in dijet_sfs_allvector_fromA2.items()}

  # Extract contours
  monojet_contours_vector = {}
  dijet_contours_vector = {}  
  for coupling in monojet_exclusiondepths_vector.keys() :
    monojet_depth = monojet_exclusiondepths_vector[coupling]
    monojet_contours_vector[coupling] = get_contours(monojet_mmed, monojet_mdm, monojet_depth)[0]
    dijet_depth = dijet_exclusiondepths_vector[coupling]
    dijet_contours_vector[coupling] = get_contours(target_scan_A2.mmed, target_scan_A2.mdm, dijet_depth)[0]

  # Add dilepton matching all of the above for hl-lhc only
  dilepton_contours_axial = {}
  dilepton_contours_vector = {}
  dilepton_exclusiondepths_axial = {}
  dilepton_exclusiondepths_vector = {}
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

      # Get exclusion depths
      dilepton_sfs_allaxial_fromA2 = rescaler_fromA2.rescale_by_br_leptons(test_gq,test_gdm,test_gl,'axial')
      dilepton_sfs_allvector_fromA2 = rescaler_fromA2.rescale_by_br_leptons(test_gq,test_gdm,test_gl,'vector')
      dilepton_exclusiondepths_axial = {k : dilepton_exdepth_A2/v for k, v in dilepton_sfs_allaxial_fromA2.items()}
      dilepton_exclusiondepths_vector = {k : dilepton_exdepth_A2/v for k, v in dilepton_sfs_allvector_fromA2.items()}

      # Get contours
      for coupling in dilepton_exclusiondepths_axial.keys() :
        dilepton_depth = dilepton_exclusiondepths_axial[coupling]
        # Check for nans: if all nan, just continue
        cleanx, cleany, cleanz = clean_grid(target_scan_A2.mmed, target_scan_A2.mdm,dilepton_depth)
        if (cleanz.size>0) : dilepton_contours_axial[coupling] = get_contours(cleanx, cleany, cleanz)[0]
      for coupling in dilepton_exclusiondepths_vector.keys() :
        dilepton_depth = dilepton_exclusiondepths_vector[coupling]
        cleanx, cleany, cleanz = clean_grid(target_scan_A2.mmed, target_scan_A2.mdm,dilepton_depth)
        if (cleanz.size>0) : dilepton_contours_vector[coupling] = get_contours(cleanx, cleany, cleanz)[0]

      # TODO add CMS when available ...

  # Make some rough plots to validate everything 
  for coupling in dijet_exclusiondepths_axial.keys() :
    grid_list_axial = [
      [monojet_mmed, monojet_mdm, monojet_exclusiondepths_axial[coupling]],
      [target_scan_A2.mmed, target_scan_A2.mdm, dijet_exclusiondepths_axial[coupling]]
    ]
    grid_list_vector = [
      [monojet_mmed, monojet_mdm, monojet_exclusiondepths_vector[coupling]],
      [target_scan_A2.mmed, target_scan_A2.mdm, dijet_exclusiondepths_vector[coupling]]
    ]
    if "hl-lhc" in collider and coupling[2] > 0.0 :
      clean_x_axial, clean_y_axial, clean_z_axial = clean_grid(target_scan_A2.mmed, target_scan_A2.mdm, dilepton_exclusiondepths_axial[coupling])
      grid_list_axial.append([clean_x_axial, clean_y_axial, clean_z_axial])
      clean_x_vector, clean_y_vector, clean_z_vector = clean_grid(target_scan_A2.mmed, target_scan_A2.mdm, dilepton_exclusiondepths_vector[coupling])
      grid_list_vector.append([clean_x_vector, clean_y_vector, clean_z_vector])

    drawContourPlotRough(grid_list_axial, addPoints = False, this_tag = "axial_gq{0}_gdm{1}_gl{2}".format(coupling[0],coupling[1],coupling[2]),plot_path = plot_path)
    drawContourPlotRough(grid_list_vector, addPoints = False, this_tag = "vector_gq{0}_gdm{1}_gl{2}".format(coupling[0],coupling[1],coupling[2]),plot_path = plot_path)

  # Save output in a clean way so that paper plot making script can be separate without re-running
  with open("vector_exclusion_depths.pkl", "wb") as outfile_vec_depths :
    out_dict = {"dijet" : dijet_exclusiondepths_vector,
                "monojet" : monojet_exclusiondepths_vector,
                "dilepton" : dilepton_exclusiondepths_vector}
    pickle.dump(out_dict, outfile_vec_depths)
  with open("axial_exclusion_depths.pkl", "wb") as outfile_axial_depths :
    out_dict = {"dijet" : dijet_exclusiondepths_axial,
                "monojet" : monojet_exclusiondepths_axial,
                "dilepton" : dilepton_exclusiondepths_axial}
    pickle.dump(out_dict, outfile_axial_depths)    
  with open("vector_exclusion_contours.pkl", "wb") as poly_file:
    out_dict = {"dijet" : dijet_contours_vector,
                "monojet" : monojet_contours_vector,
                "dilepton" : dilepton_contours_vector}
    pickle.dump(out_dict, poly_file, pickle.HIGHEST_PROTOCOL)    
  with open("axial_exclusion_contours.pkl", "wb") as poly_file:
    out_dict = {"dijet" : dijet_contours_axial,
                "monojet" : monojet_contours_axial,
                "dilepton" : dilepton_contours_axial}
    pickle.dump(out_dict, poly_file, pickle.HIGHEST_PROTOCOL)