import pickle

from basic_plotter import *

test_gq = [0.01, 0.1, 0.25]
test_gdm = [1.0,0.8,0.6,0.4,0.2]
test_gl = [0.1,0.05,0.01,0.0]

plotlims = {'hl-lhc' : (7500,1e-2,0.5), 'fcc-hh' : (15000,1e-3,0.5)}

masslines = {
  "dmDecoupled" : "m$_{\chi}$ = 100 TeV",
  "DPLike" : r"m$_{\chi}$ = m$_{\rm med}/3$",
  "dmLight" : "m$_{\chi}$ = 1 GeV",
}
ylabels = {
    "gl_lim" : r"g$_{l}$",
    "gq_lim" : r"g$_{q}$",
    "gdm_lim" : r"g$_{\chi}$",
}

def transform_coupling(couplingstring) : 
    outstring = couplingstring.replace("gq",r"$g_{q} = $").replace("gl",r"$g_{l} = $")
    outstring = outstring.replace("gdm","g$_{\chi}$ = ")
    return outstring

# Load pickle files with polygons.
for collider in ['hl-lhc'] : #, 'fcc-hh'] :
    for model in ['vector','axial'] :
        print("Starting model",model)
        with open('{0}_exclusion_contours_couplingmass_{1}.pkl'.format(model,collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)

            # Make new dict re-sorted by what coupling the limit is with respect to
            dict_bycoupling = {}

            for signature in ['dijet','monojet','dilepton'] :
                signature_contours = loaded_polygons[signature]
                for key,contour in signature_contours.items() :
                    tokens = key.split("_")
                    print(tokens)
                    couplingscan = "_".join([tokens[2],tokens[3]])
                    if not couplingscan in dict_bycoupling.keys() :
                        dict_bycoupling[couplingscan] = {}
                    dmhypothesis = tokens[4]
                    if not dmhypothesis in dict_bycoupling[couplingscan].keys() :
                        dict_bycoupling[couplingscan][dmhypothesis] = {}
                    couplings = "_".join([tokens[5],tokens[6]])
                    if not couplings in dict_bycoupling[couplingscan][dmhypothesis].keys() :
                        dict_bycoupling[couplingscan][dmhypothesis][couplings] = {}
                    if contour: dict_bycoupling[couplingscan][dmhypothesis][couplings][signature] = contour
            
            # And now we loop.
            for couplingscan in dict_bycoupling.keys() :
                all_couplingsets = []
                for dmhypothesis in dict_bycoupling[couplingscan].keys() :
                    for coupling_set, analysis_list in dict_bycoupling[couplingscan][dmhypothesis].items() :
                        if not coupling_set in all_couplingsets : all_couplingsets.append(coupling_set)
                        # First set of plots: all analyses contributing to a particular limit scenario
                        couplings_separate = coupling_set.split("_")
                        tag = "{0}_{1}_{2}_{3}".format(model,couplingscan,dmhypothesis,coupling_set)
                        label_line =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),masslines[dmhypothesis],transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                        legend_lines = []
                        contours_list = []
                        for name, contours in analysis_list.items() :
                            if contours : 
                                legend_lines.append(name)
                                contours_list.append(contours)
                        print(contours_list)
                        print(legend_lines)                                
                        drawCouplingMassPlot(contours_list, legend_lines, this_tag = tag, plot_path = "plots/couplingmass/"+collider, ylabel=ylabels[couplingscan], addText=label_line, xhigh=plotlims[collider][0], ylow=plotlims[collider][1], yhigh=plotlims[collider][2],use_colourscheme=True)
            
                # Another version: all hypotheses with same couplings on one plot
                for coupling_set in all_couplingsets :
                    couplings_separate = coupling_set.split("_")
                    contours_list_dijet = {}
                    contours_list_monojet = {}
                    contours_list_dilepton = {}
                    hypotheses_simple = ["dmDecoupled","DPLike","dmLight"]
                    for dmhypothesis in hypotheses_simple :
                        analysis_list = dict_bycoupling[couplingscan][dmhypothesis][coupling_set]
                        if "dijet" in analysis_list.keys() : contours_list_dijet[dmhypothesis] = analysis_list["dijet"]
                        if "monojet" in analysis_list.keys() : contours_list_monojet[dmhypothesis] = analysis_list["monojet"]
                        if "dilepton" in analysis_list.keys() : contours_list_dilepton[dmhypothesis] = analysis_list["dilepton"]
                    
                    tag = "{0}_{1}_{2}".format(model,couplingscan,coupling_set)
                    label_line_dijet =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),"Dijet",transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                    label_line_monojet =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),"Monojet",transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                    label_line_dilepton =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),"Dilepton",transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                    drawCouplingMassPlot([contours_list_dijet[i] for i in hypotheses_simple if i in contours_list_dijet.keys()], [masslines[i] for i in hypotheses_simple if i in contours_list_dijet.keys()], this_tag = tag+"_dijet", plot_path = "plots/couplingmass/"+collider, ylabel=ylabels[couplingscan], addText=label_line_dijet, xhigh=plotlims[collider][0], ylow=plotlims[collider][1], yhigh=plotlims[collider][2])
                    drawCouplingMassPlot([contours_list_monojet[i] for i in hypotheses_simple if i in contours_list_monojet.keys()], [masslines[i] for i in hypotheses_simple if i in contours_list_monojet.keys()], this_tag = tag+"_monojet", plot_path = "plots/couplingmass/"+collider, ylabel=ylabels[couplingscan], addText=label_line_monojet, xhigh=plotlims[collider][0], ylow=plotlims[collider][1], yhigh=plotlims[collider][2])
                    drawCouplingMassPlot([contours_list_dilepton[i] for i in hypotheses_simple if i in contours_list_dilepton.keys()], [masslines[i] for i in hypotheses_simple if i in contours_list_dilepton.keys()], this_tag = tag+"_dilepton", plot_path = "plots/couplingmass/"+collider, ylabel=ylabels[couplingscan], addText=label_line_dilepton, xhigh=plotlims[collider][0], ylow=plotlims[collider][1], yhigh=plotlims[collider][2])
