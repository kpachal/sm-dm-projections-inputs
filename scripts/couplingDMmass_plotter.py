import pickle

from basic_plotter import *

plotlims = {'hl-lhc' : 
    {"dijet" : (2500,1e-2,0.5), 
    "monojet" : (1000,1e-2,0.5),
    "dilepton" : (2500,1e-2,0.5),
    "merged" : (2500,1e-2,0.5)
    },
    'fcc-hh' : 
    {"dijet" : (5000,1e-3,0.5),
    "monojet" : (2000,1e-3,0.5),
    "dilepton" : (5000,1e-3,0.5),
    "merged" : (5000,1e-3,0.5)
    }
}

masslines = {
  "ratio2p5_fixedMMed" : r"m$_{\rm med}$ = 2.5 m$_{\chi}$",
  "ratio3_fixedMMed" : r"m$_{\rm med}$ = 3 m$_{\chi}$",
  "ratio5_fixedMMed" : r"m$_{\rm med}$ = 5 m$_{\chi}$",
  "ratio10_fixedMMed" : r"m$_{\rm med}$ = 10 m$_{\chi}$",
  "ratio100_fixedMMed" : r"m$_{\rm med}$ = 100 m$_{\chi}$",
  "ratio1000_fixedMMed" : r"m$_{\rm med}$ = 1000 m$_{\chi}$",    
  "mMed50GeV_fixedMMed" : r"m$_{\rm med}$ = 50 GeV",
  "mMed1000GeV_fixedMMed" : r"m$_{\rm med}$ = 1000 GeV"  
}
ylabels = {
    "gl_lim" : r"g$_{l}$",
    "gq_lim" : r"g$_{q}$",
    "gdm_lim" : r"g$_{\chi}$",
}

xlabel = "m$_\chi$ [GeV]"

def transform_coupling(couplingstring) : 
    outstring = couplingstring.replace("gq",r"$g_{q} = $").replace("gl",r"$g_{l} = $")
    outstring = outstring.replace("gdm","g$_{\chi}$ = ")
    return outstring

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

# Load pickle files with polygons.
for collider in ['hl-lhc', 'fcc-hh'] :
    for model in ['vector','axial'] :
        print("Starting model",model)
        #with open('{0}_exclusion_contours_couplingDMmass_{1}.pkl'.format(model,collider), "rb") as poly_file:
        with open('{0}_exclusion_contours_couplingmDM_{1}.pkl'.format(model,collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)

            # Make new dict re-sorted by what coupling the limit is with respect to
            dict_bycoupling = {}

            for signature in ['dijet','monojet','dilepton'] :
                signature_contours = loaded_polygons[signature]
                for key,contour in signature_contours.items() :
                    tokens = key.split("_")
                    couplingscan = "_".join([tokens[2],tokens[3]])
                    if not couplingscan in dict_bycoupling.keys() :
                        dict_bycoupling[couplingscan] = {}
                    dmhypothesis = "_".join([tokens[4],tokens[5]])
                    if not dmhypothesis in dict_bycoupling[couplingscan].keys() :
                        dict_bycoupling[couplingscan][dmhypothesis] = {}
                    couplings = "_".join([tokens[6],tokens[7]])
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
                        tag = "mDM_{0}_{1}_{2}_{3}".format(model,couplingscan,dmhypothesis.split("_")[0],coupling_set)
                        label_line =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),masslines[dmhypothesis],transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                        legend_lines = []
                        contours_list = []
                        for name, contours in analysis_list.items() :
                            if contours : 
                                legend_lines.append(name)
                                contours_list.append(contours)
                        drawCouplingMassPlot(contours_list, legend_lines, this_tag = tag, plot_path = "plots/couplingmass/"+collider, xlabel=xlabel, ylabel=ylabels[couplingscan], addText=label_line, xhigh=plotlims[collider]["merged"][0], ylow=plotlims[collider]["merged"][1], yhigh=plotlims[collider]["merged"][2],use_colourscheme=True)
            
                    # No need for multiple hypotheses on one plot here: we have only one, so let's just do it.
                    # Instead let's overlay other couplings.
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
                            tag = "mDM_{0}_{1}_{2}_{3}_{4}".format(model,couplingscan,dmhypothesis, fixedcoupling, analysis)
                            ordered_couplings = list(dictionary[fixedcoupling].keys())
                            ordered_couplings.sort(reverse=True)
                            for variedcoupling in ordered_couplings :
                                contour = dictionary[fixedcoupling][variedcoupling]
                                overlays.append(contour)
                                legend_lines.append(get_legend_line(variedcoupling))
                            label_line = "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),analysis.capitalize(),masslines[dmhypothesis], transform_coupling(fixedcoupling))
                            drawCouplingMassPlot(overlays, legend_lines, this_tag = tag, plot_path = "plots/couplingmass/"+collider, xlabel=xlabel, ylabel=ylabels[couplingscan], addText=label_line, xhigh=plotlims[collider][analysis][0], ylow=plotlims[collider][analysis][1], yhigh=plotlims[collider][analysis][2], is_scaling=True, transluscent=True)