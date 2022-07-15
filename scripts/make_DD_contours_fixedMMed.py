import numpy as np
from shapely.geometry import Polygon as shapely_pol
import pickle
from basic_plotter import *
import ROOT
import itertools as it
import sys
sys.path.insert(1, '../inputs/directdetection')
import collect_dd

# Define an option for each here
def test_mass_scenarios(scenario,mymdm) : 
  if scenario == "ratio2p5_fixedMMed" : return 2.5*mymdm
  elif scenario == "ratio3_fixedMMed" : return 3.*mymdm
  elif scenario == "ratio5_fixedMMed" : return 5.*mymdm
  elif scenario == "ratio10_fixedMMed" : return 10.*mymdm
  elif scenario == "ratio100_fixedMMed" : return 100.*mymdm
  elif scenario == "ratio1000_fixedMMed" : return 1000.*mymdm

masslines = {
  "ratio2p5" : r"m$_{\rm med}$ = 2.5 m$_{\chi}$",
  "ratio3" : r"m$_{\rm med}$ = 3 m$_{\chi}$",
  "ratio5" : r"m$_{\rm med}$ = 5 m$_{\chi}$",
  "ratio10" : r"m$_{\rm med}$ = 10 m$_{\chi}$",
  "ratio100" : r"m$_{\rm med}$ = 100 m$_{\chi}$",
  "ratio1000" : r"m$_{\rm med}$ = 1000 m$_{\chi}$"
}

sd_proton = collect_dd.get_sd_proton()
sd_neutron = collect_dd.get_sd_neutron()
spin_independent = collect_dd.get_spin_independent()

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

# eqn 4.10 https://arxiv.org/pdf/1603.04156.pdf#page12
def calculate_sd(gq, gdm, gl, mMed, mdm) :

    mn = 0.939 # GeV
    val = 2.4e-42 * (gq*gdm/0.25)**2 * (1000./mMed)**4 * (mn*mdm/(mn+mdm))**2
    return val

# eqn 4.3: https://arxiv.org/pdf/1603.04156.pdf#page12
def calculate_si(gq, gdm, gl, mMed, mdm) :

    mn = 0.939 # GeV
    val = 6.9e-41 * (gq*gdm/0.25)**2 * (1000./mMed)**4 * (mn*mdm/(mn+mdm))**2
    #print("gq, gdm, gl, mMed, mdm = ",gq,gdm,gl,mMed,mdm)
    #print("  val =",val)
    return val

def interpolate_vertical(x1,y1,x2,y2,n) :
    ov = sorted([[y1,x1],[y2,x2]])
    new_ys = np.linspace(ov[0][0], ov[1][0], num=n)
    new_xs = np.interp(new_ys,[ov[0][0],ov[1][0]],[ov[0][1],ov[1][1]])
    vertices = list(zip(new_xs,new_ys))
    return vertices if y1 < y2 else list(reversed(vertices))

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = it.tee(iterable)
    next(b, None)
    return zip(a, b)

# Collect direct detection contours for comparison plots.
def get_dd_lines(lineinfo) :
    legend_lines = list(lineinfo.keys())
    dd_lines = []
    for name in legend_lines :
        dd_lines.append(lineinfo[name])
    return legend_lines, dd_lines

def make_plots(collider, model, contours, legend_lines, fix_couplings, couplingscan, dmhypothesis, extra_tag = "") :
    xlow = 1
    xhigh = 2000 if 'hl-lhc' in collider else 4000
    if 'vector' in model : ylow = 1e-48 if 'hl-lhc' in collider else 1e-50
    else : ylow = 1e-46 if 'hl-lhc' in collider else 1e-50
    yhigh = 1e-37 if 'hl-lhc' in collider else 1e-42    
    usepath = "plots/directdetection/"+collider
    tag_line = "{0}_{1}_{2}_{3}".format(model,couplingscan,dmhypothesis,fix_couplings)
    if extra_tag : tag_line = tag_line + "_" + extra_tag
    couplings_separate = fix_couplings.split("_")
    # At least 2 fixed couplings.
    if len(couplings_separate) > 1 :
        label_line =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),masslines[dmhypothesis],transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
    else :
        label_line = "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),analysis.capitalize(),masslines[dmhypothesis], transform_coupling(couplings_separate[0]))
    # And draw. First, version without DD experiment lines
    # Then draw the plots with DD lines on
    if 'vector' in model :
        formatted_lines = get_dd_lines(spin_independent)
        use_ylabel = "$\sigma_{SI}$ ($\chi$-nucleon) [cm$^2$]"
        drawDDPlot(contours,legend_lines, this_tag = tag_line, plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=True, transluscent=True, xhigh=xhigh, ylow=ylow, yhigh=yhigh)
        drawDDPlot(contours,legend_lines, this_tag = tag_line+"_withDD", plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=True, transluscent=True, xhigh=xhigh, ylow=ylow, yhigh=yhigh, dd_curves = formatted_lines[1], dd_legendlines = formatted_lines[0])
    else :
        use_ylabel = "$\sigma_{SD}$ ($\chi$-nucleon) [cm$^2$]"
        drawDDPlot(contours,legend_lines, this_tag = tag_line, plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=True, transluscent=True, xhigh=xhigh, ylow=ylow, yhigh=yhigh)
        formatted_lines_p = get_dd_lines(sd_proton)
        use_ylabel = "$\sigma_{SD}$ ($\chi$-proton) [cm$^2$]"
        drawDDPlot(contours,legend_lines, this_tag = tag_line+"_withDD_proton", plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=True, transluscent=True, xhigh=xhigh, ylow=ylow, yhigh=yhigh, dd_curves = formatted_lines_p[1], dd_legendlines = formatted_lines_p[0])  
        formatted_lines_n = get_dd_lines(sd_neutron)
        use_ylabel = "$\sigma_{SD}$ ($\chi$-neutron) [cm$^2$]"
        drawDDPlot(contours,legend_lines, this_tag = tag_line+"_withDD_neutron", plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=True, transluscent=True, xhigh=xhigh, ylow=ylow, yhigh=yhigh, dd_curves = formatted_lines_n[1], dd_legendlines = formatted_lines_n[0])

# Load pickle files with polygons
for collider in ['hl-lhc', 'fcc-hh'] : 
    for model in ['vector','axial'] :        
        ylabel = "$\sigma_{SD}$" if 'axial' in model else "$\sigma_{SI}$"
        ylabel = ylabel + " ($\chi$-nucleon) [cm$^2$]" # No difference between proton & neutron for SD unless comparing to other limits  
        exclusions_dd = {}      
        with open('{0}_exclusion_contours_couplingDMmass_{1}.pkl'.format(model,collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)

            # Make new dict re-sorted by what coupling the limit is with respect to
            # and convert the contours to DD curves.
            for signature in ['dijet','monojet','dilepton'] :
                signature_contours = loaded_polygons[signature]
                for key,exclusions in signature_contours.items() :
                    if not exclusions : continue
                    tokens = key.split("_")
                    couplingscan = "_".join([tokens[2],tokens[3]])
                    dmhypothesis = "_".join([tokens[4],tokens[5]])
                    couplings = "_".join([tokens[6],tokens[7]])
                    
                    # Get everything required to convert to DD contour
                    converted_contours = []
                    for contour in exclusions :
                        original_vertices = list(contour.exterior.coords)
                        use_vertices = []
                        # Reconnect loop by adding first one back
                        for (x, y), (xnext,ynext) in pairwise(original_vertices+[original_vertices[0]]):
                            use_vertices.append((x, y))
                            if (y > 0 or ynext > 0) and ynext < 50 :
                                new_vertices = interpolate_vertical(x,y,xnext,ynext,100)
                                use_vertices = use_vertices + new_vertices
                        contour_DD_raw = []
                        raw_graph = ROOT.TGraph()                                
                        for (x, y) in use_vertices :
                            mdm = x
                            mmed = test_mass_scenarios(dmhypothesis,mdm)
                            #print("Got mmed=",mmed,"for mdm",mdm, "hypothesis",dmhypothesis)
                            if "gq" in couplingscan :
                                gq = y
                                gdm = eval(tokens[6].replace("gdm","")) if "gdm" in tokens[6] else eval(tokens[7].replace("gdm",""))
                                gl = eval(tokens[6].replace("gl","")) if "gl" in tokens[6] else eval(tokens[7].replace("gl",""))
                            elif "gdm" in couplingscan :
                                gq = eval(tokens[6].replace("gq","")) if "gq" in tokens[6] else eval(tokens[7].replace("gq",""))
                                gdm = y
                                gl = eval(tokens[6].replace("gl","")) if "gl" in tokens[6] else eval(tokens[7].replace("gl",""))
                            elif "gl" in couplingscan :
                                gq = eval(tokens[6].replace("gq","")) if "gq" in tokens[6] else eval(tokens[7].replace("gq",""))
                                gdm = eval(tokens[6].replace("gdm","")) if "gdm" in tokens[6] else eval(tokens[7].replace("gdm",""))
                                gl = y
                            else : 
                                print("Error!!!")
                                exit(1)
                            # Collect cross section limit
                            if model=='axial' :
                                sigma = calculate_sd(gq, gdm, gl, mmed, mdm)
                            else :
                                sigma = calculate_si(gq, gdm, gl, mmed, mdm)
                            contour_DD_raw.append((mdm, sigma))
                            raw_graph.AddPoint(mdm, sigma)
                        contour_DD = shapely_pol(contour_DD_raw)
                        converted_contours.append(contour_DD)
                    # Store in an orderly fashion
                    if not couplingscan in exclusions_dd.keys() :
                        exclusions_dd[couplingscan] = {}
                    if not dmhypothesis in exclusions_dd[couplingscan].keys() :
                        exclusions_dd[couplingscan][dmhypothesis] = {}
                    if not couplings in exclusions_dd[couplingscan][dmhypothesis].keys() : 
                        exclusions_dd[couplingscan][dmhypothesis][couplings] = {}
                    if contour : exclusions_dd[couplingscan][dmhypothesis][couplings][signature] = converted_contours

            # And now we loop.
            for couplingscan in exclusions_dd.keys() :
                all_couplingsets = []
                for dmhypothesis in exclusions_dd[couplingscan].keys() :
                    for coupling_set, analysis_list in exclusions_dd[couplingscan][dmhypothesis].items() :
                        if not coupling_set in all_couplingsets : all_couplingsets.append(coupling_set)
                        # First set of plots: all analyses contributing to a particular limit scenario
                        legend_lines = []
                        contours_list = []
                        for name, contours in analysis_list.items() :
                            if contours : 
                                legend_lines.append(name)
                                contours_list.append(contours)
                        make_plots(collider, model, contours_list, legend_lines, coupling_set, couplingscan, dmhypothesis.split("_")[0])
            
                    # No need for multiple hypotheses on one plot here: we have only one, so let's just do it.
                    # Instead let's overlay other couplings.
                    contours_list_dijet = {}
                    contours_list_monojet = {}
                    contours_list_dilepton = {}
                    contours_list_merged = {}
                    for coupling_set in all_couplingsets :
                        couplings_separate = coupling_set.split("_")
                        if coupling_set not in exclusions_dd[couplingscan][dmhypothesis].keys ():
                            continue
                        analysis_list = exclusions_dd[couplingscan][dmhypothesis][coupling_set]
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
                            ordered_couplings = list(dictionary[fixedcoupling].keys())
                            ordered_couplings.sort(reverse=True)
                            for variedcoupling in ordered_couplings :
                                contour = dictionary[fixedcoupling][variedcoupling]
                                overlays.append(contour)
                                legend_lines.append(get_legend_line(variedcoupling))
                            make_plots(collider, model, overlays, legend_lines, fixedcoupling, couplingscan, dmhypothesis.split("_")[0], extra_tag = analysis)

            
        # Save contours for reuse
        with open("{0}_DD_contours_fixedMassRatio_{1}.pkl".format(model,collider), "wb") as poly_file:
            pickle.dump(exclusions_dd, poly_file, pickle.HIGHEST_PROTOCOL)    