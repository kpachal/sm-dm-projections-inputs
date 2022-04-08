import pickle

from basic_plotter import *

test_gq = [0.25,0.15,0.1,0.05,0.01]
test_gdm = [1.0,0.8,0.6,0.4,0.2]
test_gl = [0.1,0.05,0.01,0.0]

plotlims = {'hl-lhc' : 
    {"dijet" : (7500,1e-2,0.5), 
    "monojet" : (3000,1e-2,0.5),
    "dilepton" : (7500,1e-2,0.5),
    "merged" : (3500,1e-2,0.5)
    },
    'fcc-hh' : 
    {"dijet" : (15000,1e-3,0.5),
    "monojet" : (6000,1e-3,0.5),
    "dilepton" : (15000,1e-3,0.5),
    "merged" : (15000,1e-3,0.5)
    }
}

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

xlabel = r"m$_{\rm med}$ [GeV]"

def get_legend_line(couplingstring) :
    if "gdm" in couplingstring :
        number = couplingstring.replace("gdm","")
        line = "g$_{\chi}$="+number
    elif "gl" in couplingstring :
        number = couplingstring.replace("gl","")
        line = r"g$_{l}$="+number
    elif "gq" in couplingstring :
        number = couplingstring.replace("gq","")
        line = r"g$_{q}$="+number
    return line

def transform_coupling(couplingstring) : 
    outstring = couplingstring.replace("gq",r"$g_{q} = $").replace("gl",r"$g_{l} = $")
    outstring = outstring.replace("gdm","g$_{\chi}$ = ")
    return outstring

# Load pickle files with polygons.
for collider in ['hl-lhc', 'fcc-hh'] :
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
                        tag = "mMed_{0}_{1}_{2}_{3}".format(model,couplingscan,dmhypothesis,coupling_set)
                        label_line =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),masslines[dmhypothesis],transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                        legend_lines = []
                        contours_list = []
                        for name, contours in analysis_list.items() :
                            if contours : 
                                legend_lines.append(name)
                                contours_list.append(contours)                              
                        drawCouplingMassPlot(contours_list, legend_lines, this_tag = tag, plot_path = "plots/couplingmass/"+collider, ylabel=ylabels[couplingscan], addText=label_line, xhigh=plotlims[collider]["merged"][0], ylow=plotlims[collider]["merged"][1], yhigh=plotlims[collider]["merged"][2],use_colourscheme=True)

                    # Now let's overlay varying couplings couplings.
                    contours_list_dijet = {}
                    contours_list_monojet = {}
                    contours_list_dilepton = {}
                    contours_list_merged = {}
                    for coupling_set in all_couplingsets :
                        couplings_separate = coupling_set.split("_")
                        analysis_list = dict_bycoupling[couplingscan][dmhypothesis][coupling_set]
                        merge_list = []

                        # Everything will enter the dictionary twice, cross-listed, so it's easier to make full plots.
                        for (analysis, dictionary) in [("dijet",contours_list_dijet),("monojet",contours_list_monojet),("dijet",contours_list_dijet)] :
                            if couplings_separate[0] not in dictionary.keys() :
                                dictionary[couplings_separate[0]] = {}
                            if couplings_separate[1] not in dictionary.keys() :
                                dictionary[couplings_separate[1]] = {}
                            if analysis in analysis_list.keys() : 
                                dictionary[couplings_separate[0]][couplings_separate[1]] = analysis_list[analysis]
                                dictionary[couplings_separate[1]][couplings_separate[0]] = analysis_list[analysis]
                                merge_list.append(analysis_list[analysis])
                        merged_contours = merge_exclusions(merge_list)
                        if not merged_contours : continue
                        if couplings_separate[0] not in contours_list_merged.keys() :
                            contours_list_merged[couplings_separate[0]] = {}
                        if couplings_separate[1] not in contours_list_merged.keys() :
                            contours_list_merged[couplings_separate[1]] = {}
                        contours_list_merged[couplings_separate[0]][couplings_separate[1]] = merged_contours
                        contours_list_merged[couplings_separate[1]][couplings_separate[0]] = merged_contours

                    # Hold one of the two constant and display others overlapping
                    for (analysis, dictionary) in [("dijet",contours_list_dijet),("monojet",contours_list_monojet),("dijet",contours_list_dijet), ("merged", contours_list_merged)] :
                        for fixedcoupling in dictionary.keys():
                            legend_lines = []
                            overlays = []
                            tag = "mMed_{0}_{1}_{2}_{3}_{4}".format(model,couplingscan,dmhypothesis, fixedcoupling, analysis)
                            ordered_couplings = list(dictionary[fixedcoupling].keys())
                            ordered_couplings.sort(reverse=True)
                            for variedcoupling in ordered_couplings :
                                contour = dictionary[fixedcoupling][variedcoupling]
                                overlays.append(contour)
                                legend_lines.append(get_legend_line(variedcoupling))
                            label_line = "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),analysis.capitalize(),masslines[dmhypothesis], transform_coupling(fixedcoupling))
                            drawCouplingMassPlot(overlays, legend_lines, this_tag = tag, plot_path = "plots/couplingmass/"+collider, xlabel=xlabel, ylabel=ylabels[couplingscan], addText=label_line, xhigh=plotlims[collider][analysis][0], ylow=plotlims[collider][analysis][1], yhigh=plotlims[collider][analysis][2], is_scaling=True, transluscent=True)                        
            
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
                    
                    tag = "mMed_{0}_{1}_{2}".format(model,couplingscan,coupling_set)
                    label_line_dijet =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),"Dijet",transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                    label_line_monojet =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),"Monojet",transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                    label_line_dilepton =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),"Dilepton",transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                    drawCouplingMassPlot([contours_list_dijet[i] for i in hypotheses_simple if i in contours_list_dijet.keys()], [masslines[i] for i in hypotheses_simple if i in contours_list_dijet.keys()], this_tag = tag+"_dijet", plot_path = "plots/couplingmass/"+collider, ylabel=ylabels[couplingscan], addText=label_line_dijet, xhigh=plotlims[collider]["dijet"][0], ylow=plotlims[collider]["dijet"][1], yhigh=plotlims[collider]["dijet"][2])
                    drawCouplingMassPlot([contours_list_monojet[i] for i in hypotheses_simple if i in contours_list_monojet.keys()], [masslines[i] for i in hypotheses_simple if i in contours_list_monojet.keys()], this_tag = tag+"_monojet", plot_path = "plots/couplingmass/"+collider, ylabel=ylabels[couplingscan], addText=label_line_monojet, xhigh=plotlims[collider]["monojet"][0], ylow=plotlims[collider]["monojet"][1], yhigh=plotlims[collider]["monojet"][2])
                    drawCouplingMassPlot([contours_list_dilepton[i] for i in hypotheses_simple if i in contours_list_dilepton.keys()], [masslines[i] for i in hypotheses_simple if i in contours_list_dilepton.keys()], this_tag = tag+"_dilepton", plot_path = "plots/couplingmass/"+collider, ylabel=ylabels[couplingscan], addText=label_line_dilepton, xhigh=plotlims[collider]["dilepton"][0], ylow=plotlims[collider]["dilepton"][1], yhigh=plotlims[collider]["dilepton"][2])
