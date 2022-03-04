import pickle

from basic_plotter import *

test_gq = [0.25,0.15,0.1,0.01]
test_gl = [0.1,0.05,0.01,0.0]

# Load pickle files with polygons
for model in ['vector','axial'] :
    with open('{0}_exclusion_contours.pkl'.format(model), "rb") as poly_file:
        loaded_polygons = pickle.load(poly_file)
        
        # Grid of plots: 
        for gq in test_gq :
            contours_list_couplingscan = []
            for gl in test_gl :
                contours_list = []
                legend_lines = []
                for signature in ['dijet','monojet','dilepton'] :
                    # e.g. no coupling to leptons, skip:
                    if (gq, 1.0, gl) not in loaded_polygons[signature].keys() :
                        continue
                    exclusions = loaded_polygons[signature][(gq, 1.0, gl)]
                    contours_list.append(exclusions)
                    legend_lines.append(signature)
                # First set of plots: 3 contours, one plot for every coupling combo
                label_line =  "{0}, g$_{5}$={2}\ng$_{4}$={1}, g$_{6}$={3}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,1.0,gl,"q","\chi","l")
                drawMassMassPlot(contours_list, legend_lines, this_tag = "gq{0}_gdm1.0_gl{1}".format(gq, gl), plot_path = "plots/massmass", addText=label_line)
            # Second set of plots: merge all contours; fix one coupling and vary the other
        # Need second set of plots with gl fixed instead:



            
