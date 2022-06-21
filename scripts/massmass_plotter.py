import pickle

from basic_plotter import *

test_gq = [0.25,0.1,0.05,0.02,0.01]
test_gdm = [1.0, 0.2, 0.1, 0.05, 0.0]
test_gl = [0.1,0.05,0.01,0.0]

plotlims = {'hl-lhc' : (7500,1200), 'fcc-hh' : (15000,5000)}

# Load pickle files with polygons
for collider in ['hl-lhc', 'fcc-hh'] :
    for model in ['vector','axial'] :
        with open('{0}_exclusion_contours_{1}.pkl'.format(model,collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)
            
            # Grid of plots: 
            for gdm in test_gdm :
                for gq in test_gq :
                    contours_list_couplingscan = []
                    legend_lines_couplingscan = []
                    for gl in test_gl :
                        contours_list = []
                        legend_lines = []
                        for signature in ['dilepton','dijet','monojet'] :
                            # e.g. no coupling to leptons, skip:
                            if (gq, gdm, gl) not in loaded_polygons[signature].keys() :
                                continue
                            exclusions = loaded_polygons[signature][(gq, gdm, gl)]
                            if exclusions : 
                                contours_list.append(exclusions)
                                legend_lines.append(signature)
                        # First set of plots: 3 contours, one plot for every coupling combo
                        label_line =  "{0}\n{7}, g$_{5}$={2}\ng$_{4}$={1}, g$_{6}$={3}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,gdm,gl,"q","\chi","l",collider.upper())
                        drawMassMassPlot(contours_list, legend_lines, this_tag = model+"_gq{0}_gdm{1}_gl{2}".format(gq, gdm, gl), plot_path = "plots/massmass/"+collider, addText=label_line, xhigh=plotlims[collider][0], yhigh=plotlims[collider][1], use_colourscheme=True)
                        full_polygons = merge_exclusions(contours_list)
                        if full_polygons :
                            contours_list_couplingscan.append(full_polygons)
                            legend_lines_couplingscan.append("g$_{0}$={1}".format("l",gl))
                    # Second set of plots: merge all contours; fix gq and vary gl.
                    # Note this is not meaningful where we don't have dilepton projections - skip then.
                    label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,gdm,"q","\chi",collider.upper())
                    drawMassMassPlot(contours_list_couplingscan,legend_lines_couplingscan, this_tag = model+"_gq{0}_gdm{1}".format(gq,gdm), plot_path = "plots/massmass/"+collider, addText = label_line,is_scaling=True, transluscent=True,xhigh=plotlims[collider][0], yhigh=plotlims[collider][1])
                    # And for analyses separately
                    for signature in ['dijet','monojet','dilepton'] :
                        label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}\n{6}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,gdm,"q","\chi",collider.upper(),signature.capitalize())
                        contours_list = []
                        legend_lines = []
                        for gl in test_gl :
                            if (gq, gdm, gl) not in loaded_polygons[signature].keys() : continue
                            exclusions = loaded_polygons[signature][(gq, gdm, gl)]
                            if exclusions :
                                contours_list.append(exclusions)
                                legend_lines.append("g$_{0}$={1}".format("l",gl))
                        drawMassMassPlot(contours_list,legend_lines,this_tag=model+"_gq{0}_gdm{1}_{2}".format(gq,gdm,signature),plot_path = "plots/massmass/"+collider, addText = label_line,is_scaling=True, transluscent=True,xhigh=plotlims[collider][0], yhigh=plotlims[collider][1])
                # Need second set of plots with gl fixed instead:
                for gl in test_gl :
                    contours_list_couplingscan = []
                    legend_lines_couplingscan = []
                    for gq in test_gq :
                        contours_list = []
                        for signature in ['dijet','monojet','dilepton'] :
                            # e.g. no coupling to leptons, skip:
                            if (gq, gdm, gl) not in loaded_polygons[signature].keys() : continue
                            exclusions = loaded_polygons[signature][(gq, gdm, gl)]
                            if exclusions :
                                contours_list.append(exclusions)
                        full_polygons = merge_exclusions(contours_list)
                        if full_polygons :
                            contours_list_couplingscan.append(full_polygons)
                            legend_lines_couplingscan.append("g$_{0}$={1}".format("q",gq))
                    label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}".format(("Axial-vector" if 'axial' in model else "Vector"),gdm,gl,"\chi","l",collider.upper())
                    drawMassMassPlot(contours_list_couplingscan,legend_lines_couplingscan, this_tag = model+"_gl{0}_gdm{1}".format(gl,gdm), plot_path = "plots/massmass/"+collider, addText = label_line,is_scaling=True,transluscent=True, xhigh=plotlims[collider][0], yhigh=plotlims[collider][1]) 
                    for signature in ['dijet','monojet','dilepton'] :
                        label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}\n{6}".format(("Axial-vector" if 'axial' in model else "Vector"),gdm,gl,"\chi","l",collider.upper(),signature.capitalize())
                        contours_list = []
                        legend_lines = []
                        for gq in test_gq :
                            if (gq, gdm, gl) not in loaded_polygons[signature].keys() : continue
                            exclusions = loaded_polygons[signature][(gq, gdm, gl)]
                            if exclusions :
                                contours_list.append(exclusions)
                                legend_lines.append("g$_{0}$={1}".format("q",gq))
                        drawMassMassPlot(contours_list,legend_lines,this_tag = model+"_gl{0}_gdm{1}_{2}".format(gl,gdm,signature), plot_path = "plots/massmass/"+collider, addText = label_line,is_scaling=True,transluscent=True, xhigh=plotlims[collider][0], yhigh=plotlims[collider][1])
            # And now third set of plots with gq and gl fixed:
            for gq in test_gq :
                for gl in test_gl :
                    contours_list_couplingscan = []
                    legend_lines_couplingscan = []
                    for gdm in test_gdm :
                        contours_list = []
                        legend_lines = []                    
                        for signature in ['dijet','monojet','dilepton'] :
                            # e.g. no coupling to leptons, skip:
                            if (gq, gdm, gl) not in loaded_polygons[signature].keys() :
                                continue                            
                            exclusions = loaded_polygons[signature][(gq, gdm, gl)]
                            if exclusions :
                                contours_list.append(exclusions)
                        if all(not i for i in contours_list) : continue
                        full_polygons = merge_exclusions(contours_list)
                        if full_polygons :
                            contours_list_couplingscan.append(full_polygons)
                            legend_lines_couplingscan.append("g$_{0}$={1}".format("\chi",gdm))
                    label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,gl,"q","l",collider.upper())
                    drawMassMassPlot(contours_list_couplingscan,legend_lines_couplingscan, this_tag = model+"_gq{0}_gl{1}".format(gq,gl), plot_path = "plots/massmass/"+collider, addText = label_line,is_scaling=True, transluscent=True,xhigh=plotlims[collider][0], yhigh=plotlims[collider][1])
                    for signature in ['dijet','monojet','dilepton'] :
                        label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}\n{5}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,gl,"q","l",collider.upper(),signature.capitalize())
                        contours_list = []
                        legend_lines = []
                        for gdm in test_gdm :
                            if (gq, gdm, gl) not in loaded_polygons[signature].keys() : continue
                            exclusions = loaded_polygons[signature][(gq, gdm, gl)]
                            if exclusions :
                                contours_list.append(exclusions)
                                legend_lines.append("g$_{0}$={1}".format("\chi",gdm))
                        drawMassMassPlot(contours_list,legend_lines, this_tag = model+"_gq{0}_gl{1}_{2}".format(gq,gl,signature), plot_path = "plots/massmass/"+collider, addText = label_line,is_scaling=True, transluscent=True,xhigh=plotlims[collider][0], yhigh=plotlims[collider][1])