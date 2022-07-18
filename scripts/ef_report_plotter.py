import pickle
from shapely.geometry import Polygon as shapely_pol

from basic_plotter import *

#####################
# Plot helpers
masslines = {
  "dmDecoupled" : "m$_{\chi}$ = 100 TeV",
  "DPLike" : r"m$_{\chi}$ = m$_{\rm med}/3$",
  "dm1GeV" : "m$_{\chi}$ = 1 GeV",
  "ratio2p5" : r"m$_{\rm med}$ = 2.5 m$_{\chi}$",
  "ratio3" : r"m$_{\rm med}$ = 3 m$_{\chi}$",
  "ratio5" : r"m$_{\rm med}$ = 5 m$_{\chi}$",
  "ratio10" : r"m$_{\rm med}$ = 10 m$_{\chi}$",
  "ratio100" : r"m$_{\rm med}$ = 100 m$_{\chi}$",
  "ratio1000" : r"m$_{\rm med}$ = 1000 m$_{\chi}$"  
}

ylabels = {
    "gl" : r"g$_{l}$",
    "gq" : r"g$_{q}$",
    "gdm" : r"g$_{\chi}$",
    "gl_lim" : r"g$_{l}$",
    "gq_lim" : r"g$_{q}$",
    "gdm_lim" : r"g$_{\chi}$",    
}

def get_couplingmass_legend_line(couplingstring) :
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

def format_coupling(coupling,value) : 
    outstring = coupling.replace("gq",r"$g_{q} = $").replace("gl",r"$g_{l} = $")
    outstring = outstring.replace("gdm","g$_{\chi}$ = ")+"{0}".format(value)
    return outstring

def transform_coupling(couplingstring) : 
    outstring = couplingstring.replace("gq",r"$g_{q} = $").replace("gl",r"$g_{l} = $")
    outstring = outstring.replace("gdm","g$_{\chi}$ = ")
    return outstring

#####################
# Fig 1: individual analysis limits.
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

#####################
# Fig 2: individual analysis limits, coupling-mass.
dmhypothesis = "dm1GeV"
fig2_contours = {}
fig2_couplings = {"gq" : {"gdm" : 1.0, "gl" : 0.0},
                "gdm" : {"gq" : 0.1, "gl" : 0.0},
                "gl" : {"gq" : 0.25, "gdm" : 1.0}}
for collider in ['hl-lhc', 'fcc-hh'] :
    for versus in ["gq","gdm","gl"] :
        if not collider in fig2_contours.keys() : fig2_contours[collider] = {}
        if not "versus_{0}".format(versus) in fig2_contours[collider].keys() : fig2_contours[collider]["versus_{0}".format(versus)] = {}
        key_couplings = fig2_couplings[versus]
        gs,vals = list(key_couplings.keys()),list(key_couplings.values())
        key = "{0}_vector_{1}_lim_{6}_{2}{3}_{4}{5}".format(collider,versus,gs[0],vals[0],gs[1],vals[1],dmhypothesis)
        with open('vector_exclusion_contours_couplingmass_{0}.pkl'.format(collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)
            for signature in ['dilepton','dijet','monojet'] :
                if key in loaded_polygons[signature].keys() :
                    fig2_contours[collider]["versus_{0}".format(versus)][signature] = loaded_polygons[signature][key]

#####################
# Make fig 2a
# TODO not yet happy with an approach to gradient fill here
contours_2a_solid = fig2_contours["fcc-hh"]["versus_gq"]
contours_2a_dashed = fig2_contours["hl-lhc"]["versus_gq"]
legend_lines = []
plot_contours_2a = []
for name,contours in contours_2a_solid.items() :
    if contours :
        legend_lines.append("{0}, FCC-hh".format(name.capitalize()))
        plot_contours_2a.append(contours)
for name,contours in contours_2a_dashed.items() :
    if contours :
        legend_lines.append("{0}, HL-LHC".format(name.capitalize()))
        plot_contours_2a.append(contours)
coupling1_2a = format_coupling(list(fig2_couplings["gq"].keys())[0],list(fig2_couplings["gq"].values())[0])
coupling2_2a = format_coupling(list(fig2_couplings["gq"].keys())[1],list(fig2_couplings["gq"].values())[1])
label_line_2a =  "Vector, {0}\n{1}, {2}".format(masslines[dmhypothesis],coupling1_2a,coupling2_2a)
drawCouplingMassPlot(plot_contours_2a, legend_lines, this_tag = "vector_couplingmass_comparecolliders_gqlimit", plot_path = "plots/efreport/", ylabel=ylabels["gq"], addText=label_line_2a, xhigh=12000, ylow=1e-3, yhigh=0.5,use_colourscheme=True, dash_contour = [False, False, True, True],legend_spot='lower right',text_spot=[0.48,0.47], gradient_fill=[True,False,False,False])

#####################
# Make fig 2b
contours_2b_solid = fig2_contours["fcc-hh"]["versus_gdm"]
# extend contour for dijet to the left
low_box = shapely_pol([[-15,1e-4],[5050,1e-4],[5050,0.5],[-15,0.5]])
new_dijet = merge_exclusions([[low_box],contours_2b_solid["dijet"]])
contours_2b_solid["dijet"] = new_dijet
# Continue
contours_2b_dashed = fig2_contours["hl-lhc"]["versus_gdm"]
legend_lines = []
plot_contours_2b = []
for name,contours in contours_2b_solid.items() :
    if contours :
        legend_lines.append("{0}, FCC-hh".format(name.capitalize()))
        plot_contours_2b.append(contours)
for name,contours in contours_2b_dashed.items() :
    if contours :
        legend_lines.append("{0}, HL-LHC".format(name.capitalize()))
        plot_contours_2b.append(contours)
coupling1_2b = format_coupling(list(fig2_couplings["gdm"].keys())[0],list(fig2_couplings["gdm"].values())[0])
coupling2_2b = format_coupling(list(fig2_couplings["gdm"].keys())[1],list(fig2_couplings["gdm"].values())[1])
label_line_2b =  "Vector, {0}\n{1}, {2}".format(masslines[dmhypothesis],coupling1_2b,coupling2_2b)
drawCouplingMassPlot(plot_contours_2b, legend_lines, this_tag = "vector_couplingmass_comparecolliders_gdmlimit", plot_path = "plots/efreport/", ylabel=ylabels["gdm"], addText=label_line_2b, xhigh=12000, ylow=1e-3, yhigh=0.5,use_colourscheme=True, dash_contour = [False, False, True, True],legend_spot='lower right',text_spot=[0.48,0.47], gradient_fill=[True,False,False,False])

#####################
# Make fig 2c
# TODO not yet happy with an approach to gradient fill here
contours_2c = fig2_contours["hl-lhc"]["versus_gl"]
# extend contour for dilepton to the left
#low_box = shapely_pol([[-15,0],[3000,0],[3000,1205],[-15,1205]])
#new_dilepton = merge_exclusions([[low_box],contours_1a_dict["dilepton"]])
#contours_1a_dict["dilepton"] = new_dilepton
legend_lines = []
plot_contours_2c = []
for name,contours in contours_2c.items() :
    if contours :
        legend_lines.append(name.capitalize())
        plot_contours_2c.append(contours)
label_line_2c =  "Vector\nHL-LHC, {0}\n{1}, {2}".format(masslines[dmhypothesis],coupling1_2b,coupling2_2b)
drawCouplingMassPlot(plot_contours_2c, legend_lines, this_tag = "vector_couplingmass_gllimit", plot_path = "plots/efreport/", ylabel=ylabels["gl"], addText=label_line_2c, xhigh=9000, ylow=5e-3, yhigh=0.5,use_colourscheme=True, legend_spot='lower right',gradient_fill=[True,False,False])

#####################
# Next figs: Like fig 2 but demonstrating the coupling scan.
dict_couplingmass_bycoupling = {}
all_couplingsets = []
for collider in ['hl-lhc', 'fcc-hh'] :
    with open('vector_exclusion_contours_couplingmass_{0}.pkl'.format(collider), "rb") as poly_file:
        loaded_polygons = pickle.load(poly_file)

        for signature in ['dijet','monojet','dilepton'] :
            signature_contours = loaded_polygons[signature]
            for key,contour in signature_contours.items() :
                tokens = key.split("_")
                couplingscan = "_".join([tokens[2],tokens[3]])
                if not couplingscan in dict_couplingmass_bycoupling.keys() :
                    dict_couplingmass_bycoupling[couplingscan] = {}
                dmhypothesis = tokens[4]
                if not dmhypothesis in dict_couplingmass_bycoupling[couplingscan].keys() :
                    dict_couplingmass_bycoupling[couplingscan][dmhypothesis] = {}
                if not collider in dict_couplingmass_bycoupling[couplingscan][dmhypothesis].keys() :
                    dict_couplingmass_bycoupling[couplingscan][dmhypothesis][collider] = {}
                couplings = "_".join([tokens[5],tokens[6]])
                if not couplings in dict_couplingmass_bycoupling[couplingscan][dmhypothesis][collider].keys() :
                    dict_couplingmass_bycoupling[couplingscan][dmhypothesis][collider][couplings] = {}
                if not couplings in all_couplingsets : all_couplingsets.append(couplings)
                if contour: dict_couplingmass_bycoupling[couplingscan][dmhypothesis][collider][couplings][signature] = contour

# Get straight to the heart of the six plots we want.
xlabel = r"m$_{\rm med}$ [GeV]"
for couplingscan in ["gq_lim","gdm_lim"] :
    for dmhypothesis in ["dm1GeV","DPLike"] :

        # Group and split by couplings
        contours_dijet = {}
        contours_monojet = {}
        for coupling_set in all_couplingsets :
            couplings_separate = coupling_set.split("_")
            analysis_list_hllhc = None
            analysis_list_fcchh = None
            if coupling_set in dict_couplingmass_bycoupling[couplingscan][dmhypothesis]["hl-lhc"].keys() :
                analysis_list_hllhc = dict_couplingmass_bycoupling[couplingscan][dmhypothesis]["hl-lhc"][coupling_set]
            if coupling_set in dict_couplingmass_bycoupling[couplingscan][dmhypothesis]["fcc-hh"].keys() :
                analysis_list_fcchh = dict_couplingmass_bycoupling[couplingscan][dmhypothesis]["fcc-hh"][coupling_set]
            # Everything will enter the dictionary twice, cross-listed, so it's easier to make full plots.
            for (analysis, dictionary) in [("dijet",contours_dijet),("monojet",contours_monojet)] :
                inner_dict = {}
                if analysis_list_hllhc and analysis in analysis_list_hllhc.keys() : 
                    inner_dict['hl-lhc'] = analysis_list_hllhc[analysis]
                if analysis_list_fcchh and analysis in analysis_list_fcchh.keys() :
                    inner_dict['fcc-hh'] = analysis_list_fcchh[analysis]
                if couplings_separate[0] not in dictionary.keys() :
                    dictionary[couplings_separate[0]] = {}
                dictionary[couplings_separate[0]][couplings_separate[1]] = inner_dict
                if couplings_separate[1] not in dictionary.keys() :
                    dictionary[couplings_separate[1]] = {}
                dictionary[couplings_separate[1]][couplings_separate[0]] = inner_dict

        # Hold one of the two constant and display others overlapping
        for (analysis, dictionary) in [("dijet",contours_dijet),("monojet",contours_monojet)] :
            # Only fixed coupling we care about is gl = 0.
            fixedcoupling = "gl0.0"
            legend_lines_solid = []
            contours_solid = []
            contours_dashed = []
            legend_lines_dashed = []
            tag = "vector_{0}_{1}_{2}".format(couplingscan, dmhypothesis, analysis)
            ordered_couplings = list(dictionary[fixedcoupling].keys())
            ordered_couplings.sort(reverse=True)
            for variedcoupling in ordered_couplings :
                these_contours = dictionary[fixedcoupling][variedcoupling]
                if 'hl-lhc' in these_contours.keys() :
                    contours_dashed.append(these_contours['hl-lhc'])
                    legend_lines_dashed.append(get_couplingmass_legend_line(variedcoupling)+", HL-LHC")
                if 'fcc-hh' in these_contours.keys() :
                    contours_solid.append(these_contours['fcc-hh'])
                    legend_lines_solid.append(get_couplingmass_legend_line(variedcoupling)+", FCC-hh")
            label_line = "Vector, {0}\n{1}, {2}".format(analysis.capitalize(),masslines[dmhypothesis], transform_coupling(fixedcoupling))
            drawCouplingMassPlot(contours_solid, legend_lines_solid, dashed_lines=contours_dashed, dashed_legends=legend_lines_dashed, this_tag = tag, plot_path = "plots/efreport/", xlabel=xlabel, ylabel=ylabels[couplingscan], addText=label_line, xhigh=(15000 if 'dijet' in analysis else 6000), ylow=1e-3, yhigh=0.5, is_scaling=True, legend_spot='outside', transluscent=True)  
        

#####################
## Figures: DD fixed mass ratios with FCC on top.
## Ignoring axial-vector entirely here.
dict_dd_bycoupling = {}
all_couplingsets = []
for collider in ['hl-lhc', 'fcc-hh'] :
    with open('vector_DD_contours_fixedMassRatio_{0}.pkl'.format(collider), "rb") as poly_file:
        loaded_polygons = pickle.load(poly_file)
        # Only care right now about gq integrated out of x axis and gl = 0.
        relevant_limits = loaded_polygons['gq_lim']
        for dmhypothesis_long in relevant_limits.keys() :
            dmhypothesis = dmhypothesis_long.split("_")[0]
            for couplingset, lims in relevant_limits[dmhypothesis_long].items() :
                if not "gl0.0" in couplingset : continue
                gdm = couplingset.split("_")[0]
                for signature in ['dijet','monojet','dilepton'] :
                    if not signature in lims.keys() :
                        contour = []
                    else : contour = lims[signature]
                    if not signature in dict_dd_bycoupling.keys() :
                        dict_dd_bycoupling[signature] = {}
                    if not dmhypothesis in dict_dd_bycoupling[signature].keys() :
                        dict_dd_bycoupling[signature][dmhypothesis] = {}
                    if not collider in dict_dd_bycoupling[signature][dmhypothesis].keys() :
                        dict_dd_bycoupling[signature][dmhypothesis][collider] = {}
                    dict_dd_bycoupling[signature][dmhypothesis][collider][gdm] = contour

# We want 1 plot per angle.
# No overlapping signatures here,
# and only care about monojet right now.
xlabel = r"m$_{\rm med}$ [GeV]"
set_of_couplings = []
for dmhypothesis, contour_set in dict_dd_bycoupling["monojet"].items() :
    contours_solid = []
    contours_dashed = []
    legend_lines_solid = []
    legend_lines_dashed = []
    for gdm,contour in contour_set["hl-lhc"].items() :
        contours_dashed.append(contour)
        legend_lines_dashed.append(get_couplingmass_legend_line(gdm)+", HL-LHC")
        if not gdm in set_of_couplings : set_of_couplings.append(gdm)
    for gdm,contour in contour_set["fcc-hh"].items() :
        contours_solid.append(contour)
        legend_lines_solid.append(get_couplingmass_legend_line(gdm)+", FCC-hh")
        if not gdm in set_of_couplings : set_of_couplings.append(gdm)
    use_ylabel = "$\sigma_{SI}$ ($\chi$-nucleon) [cm$^2$]"
    tag_line = "{0}_{1}_gl0.0_compareColliders".format("vector",dmhypothesis)
    label_line_ddcompare="Vector, Monojet\n{0}, {1}".format(masslines[dmhypothesis], transform_coupling("gl0.0"))    
    drawDDPlot(contours_solid,legend_lines_solid, this_tag = tag_line, plot_path = "plots/efreport/", addText = label_line_ddcompare, ylabel=use_ylabel, is_scaling=True, transluscent=True, xhigh=4000, ylow=1e-50, yhigh=1e-37,dashed_lines=contours_dashed, dashed_legends=legend_lines_dashed)

# Second version: all angles on one plot for a given coupling.
for gdm in set_of_couplings :
    if gdm=="gdm0.0" : continue
    contours_solid = []
    contours_dashed = []
    legend_lines_solid = []
    legend_lines_dashed = []
    for dmhypothesis, contour_set in dict_dd_bycoupling["monojet"].items() :
        if gdm in contour_set["hl-lhc"].keys() :
            contours_dashed.append(contour_set["hl-lhc"][gdm])
            legend_lines_dashed.append(masslines[dmhypothesis]+", HL-LHC")
        if gdm in contour_set["fcc-hh"].keys() :
            contours_solid.append(contour_set["fcc-hh"][gdm])
            legend_lines_solid.append(masslines[dmhypothesis]+", FCC-hh")
    use_ylabel = "$\sigma_{SI}$ ($\chi$-nucleon) [cm$^2$]"
    tag_line = "vector_{0}_gl0.0_compareColliders".format(gdm)
    label_line_ddcompare="Vector, Monojet\n{0}, {1}".format(transform_coupling(gdm), transform_coupling("gl0.0"))
    drawDDPlot(contours_solid,legend_lines_solid, this_tag = tag_line, plot_path = "plots/efreport/", addText = label_line_ddcompare, ylabel=use_ylabel, is_scaling=True, transluscent=True, xhigh=4000, ylow=1e-50, yhigh=1e-37,dashed_lines=contours_dashed, dashed_legends=legend_lines_dashed)

    

#####################
# Figures: coupling limit versus mA/mchi
# Curves should already be in dict_couplingmass_bycoupling.

# Rearrange so that we have this grouped by analysis, then other couplings, then collider, then hypothesis.
dict_couplingmass_byhypothesis = {}
for hypothesis,dicti in dict_couplingmass_bycoupling["gdm_lim"].items() :
    for collider, dictj in dicti.items() :
        for couplingset,dictk in dictj.items() :
            for analysis,contours in dictk.items() :
                if not analysis in dict_couplingmass_byhypothesis.keys() :
                    dict_couplingmass_byhypothesis[analysis] = {}
                if not couplingset in dict_couplingmass_byhypothesis[analysis].keys() :
                    dict_couplingmass_byhypothesis[analysis][couplingset] = {}
                if not collider in dict_couplingmass_byhypothesis[analysis][couplingset].keys() :
                    dict_couplingmass_byhypothesis[analysis][couplingset][collider] = {}
                dict_couplingmass_byhypothesis[analysis][couplingset][collider][hypothesis] = contours

# Version 1: just do g_DM^2 y axis for validation.
print(dict_couplingmass_byhypothesis)
test_mdm_list = ["dm1GeV","dm10GeV","dm100GeV"]
xlabel = r"m$_{\rm med}$/m$_{\chi}$"
root_tests = []
for analysis, dicta in dict_couplingmass_byhypothesis.items() :
    if "monojet" not in analysis : continue
    for couplingset, dictb in dicta.items() :
        for collider, dictc in dictb.items() :
            legendlines = []
            contour_collection = []

            for hypothesis in test_mdm_list :
                if hypothesis not in dictc.keys() :
                    continue
                original_contours = dictc[hypothesis]
                converted_contours = []
                mdm = eval(hypothesis.replace("dm","").replace("GeV",""))
                for contour in original_contours :
                    original_vertices = list(contour.exterior.coords)
                    new_vertices = []
                    for (x, y) in original_vertices :
                        gdm_2 = y**2
                        xratio = x/mdm
                        #print(xratio,gdm_2)
                        new_vertices.append((xratio,gdm_2))
                    new_contour = shapely_pol(new_vertices)
                    converted_contours.append(new_contour)
                contour_collection.append(converted_contours)
                legend_lines.append(hypothesis)
            tag_line = "newplot_{0}_{1}_{2}".format(analysis,couplingset,collider)
            label_line_new="Vector, {0}\n{0}, {1}".format(analysis.capitalize(),couplingset, collider)   
            drawDDPlot(contour_collection,legendlines, this_tag = tag_line, plot_path = "plots/efreport/", addText = label_line_new, ylabel=r"g$_{DM}^{2}$", is_scaling=True, transluscent=True, xhigh=1000, ylow=0.0001, yhigh=1)

# Version 2: convert to "y" variable using Phil's formula.