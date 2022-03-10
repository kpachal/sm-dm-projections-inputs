import pickle

from basic_plotter import *

test_gq = [0.25,0.15,0.1,0.05,0.01]
test_gl = [0.1,0.05,0.01,0.0]

plotlims = {'hl-lhc' : (7500,1200), 'fcc-hh' : (15000,5000)}

# Load pickle files with polygons
for collider in ['hl-lhc', 'fcc-hh'] :
    for model in ['vector','axial'] :
        with open('{0}_exclusion_contours_{1}.pkl'.format(model,collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)
            
            # Grid of plots: 
            for gq in test_gq :
                contours_list_couplingscan = []
                legend_lines_couplingscan = []
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
                    label_line =  "{0}\n{7}, g$_{5}$={2}\ng$_{4}$={1}, g$_{6}$={3}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,1.0,gl,"q","\chi","l",collider.upper())
                    drawMassMassPlot(contours_list, legend_lines, this_tag = model+"_gq{0}_gdm1.0_gl{1}".format(gq, gl), plot_path = "plots/massmass/"+collider, addText=label_line, xhigh=plotlims[collider][0], yhigh=plotlims[collider][1])
                    full_polygons = merge_exclusions(contours_list)
                    contours_list_couplingscan.append(full_polygons)
                    legend_lines_couplingscan.append("g$_{0}$={1}".format("l",gl))
                # Second set of plots: merge all contours; fix gq and vary gl.
                # Note this is not meaningful where we don't have dilepton projections - skip then.
                label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,1.0,"q","\chi",collider.upper())
                drawMassMassPlot(contours_list_couplingscan,legend_lines_couplingscan, this_tag = model+"_gq{0}_gdm1.0".format(gq), plot_path = "plots/massmass/"+collider, addText = label_line,is_scaling=True, xhigh=plotlims[collider][0], yhigh=plotlims[collider][1])
            # Need second set of plots with gl fixed instead:
            for gl in test_gl :
                contours_list_couplingscan = []
                legend_lines_couplingscan = []
                for gq in test_gq :
                    contours_list = []
                    for signature in ['dijet','monojet','dilepton'] :
                        # e.g. no coupling to leptons, skip:
                        if (gq, 1.0, gl) not in loaded_polygons[signature].keys() :
                            continue
                        exclusions = loaded_polygons[signature][(gq, 1.0, gl)]
                        contours_list.append(exclusions)
                    full_polygons = merge_exclusions(contours_list)
                    contours_list_couplingscan.append(full_polygons)
                    legend_lines_couplingscan.append("g$_{0}$={1}".format("q",gq))
                label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}".format(("Axial-vector" if 'axial' in model else "Vector"),1.0,gl,"\chi","l",collider.upper())
                drawMassMassPlot(contours_list_couplingscan,legend_lines_couplingscan, this_tag = model+"_gl{0}_gdm1.0".format(gl), plot_path = "plots/massmass/"+collider, addText = label_line,is_scaling=True, xhigh=plotlims[collider][0], yhigh=plotlims[collider][1])            
