import pickle
from shapely.geometry import Polygon as shapely_pol

from basic_plotter import *

#####################
# Fig 1: individual analysis limits.
# 1a) as in initial report but with dilepton beginning at 200 GeV.
fig1_contours = {}
for collider in ['hl-lhc', 'fcc-hh'] :
    for couplings in [(0.1, 1.0, 0.1),(0.1, 1.0, 0.0)] :
        if not collider in fig1_contours.keys() : fig1_contours[collider] = {}
        if not couplings in fig1_contours[collider].keys() : fig1_contours[collider][couplings] = {}        
        with open('vector_exclusion_contours_{0}.pkl'.format(collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)
            for signature in ['dilepton','dijet','monojet'] :
                # e.g. no coupling to leptons, skip:
                if couplings not in loaded_polygons[signature].keys() :
                    continue            
                exclusions = loaded_polygons[signature][couplings]
                fig1_contours[collider][couplings][signature] = exclusions

#####################
# Make fig 1a
contours_1a_dict = fig1_contours["hl-lhc"][(0.1, 1.0, 0.1)]
# extend contour for dilepton to the left
low_box = shapely_pol([[-15,0],[3000,0],[3000,1205],[-15,1205]])
new_dilepton = merge_exclusions([[low_box],contours_1a_dict["dilepton"]])
contours_1a_dict["dilepton"] = new_dilepton
# And go
analyses_1a = contours_1a_dict.keys()
contours_1a = [contours_1a_dict[i] for i in analyses_1a]
label_line_1a = "Vector\nHL-LHC, g$_{1}$={4}\ng$_{0}$={3}, g$_{2}$={5}".format("q","\chi","l",0.1,1.0,0.1)
drawMassMassPlot(contours_1a,[i.capitalize() for i in analyses_1a],this_tag="vector_massmass_hl-lhc",plot_path = "plots/efreport/", addText=label_line_1a, xhigh=7500, yhigh=1200, use_colourscheme=True,gradient_fill=[True,False,False])

#####################
# Make fig 1b
contours_1b_solid = fig1_contours["fcc-hh"][(0.1, 1.0, 0.0)]
# extend contour for dijet to the left
low_box = shapely_pol([[-15,0],[5050,0],[5050,4000],[-15,4000]])
new_dijet = merge_exclusions([[low_box],contours_1b_solid["dijet"]])
contours_1b_solid["dijet"] = new_dijet
# And go
contours_1b_dashed = fig1_contours["hl-lhc"][(0.1, 1.0, 0.0)]
legend_lines = ["{0}, FCC-hh".format(i.capitalize()) for i in contours_1b_solid.keys()]
legend_lines += ["{0}, HL-LHC".format(i.capitalize()) for i in contours_1b_dashed.keys()]
label_line_1b = "Vector, g$_{1}$={4}\ng$_{0}$={3}, g$_{2}$={5}".format("q","\chi","l",0.1,1.0,0.1)
drawMassMassPlot(list(contours_1b_solid.values())+list(contours_1b_dashed.values()),legend_lines,this_tag="vector_massmass_comparecolliders",plot_path = "plots/efreport/", addText=label_line_1b, xhigh=15000, yhigh=3500, use_colourscheme=True, dash_contour = [False, False, True, True], gradient_fill=[True,False,False,False],text_spot = [0.55,0.55])