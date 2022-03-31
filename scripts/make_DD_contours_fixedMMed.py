import numpy as np
from shapely.geometry import Polygon as shapely_pol
import pickle
from basic_plotter import *
import ROOT
import itertools as it

test_gq = [0.25,0.15,0.1,0.05,0.01]
test_gdm = [1.0,0.8,0.6,0.4,0.2]
test_gl = [0.1,0.05,0.01,0.0]

# Define an option for each here
def test_mass_scenarios(scenario,mymdm) : 
  if scenario == "DPLike_fixedMMed" : return 3*mymdm

masslines = {
  "DPLike_fixedMMed" : r"m$_{\chi}$ = m$_{\rm med}/3$",
}

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

test_coupling_scenarios = {
  "gq_lim" : {
    "test_gq" : np.logspace(np.log10(0.001),0,101),
    "test_gdm" : [0.0, 0.1, 0.2, 0.5, 1.0],
    "test_gl" : [0.0],
  },
  "gdm_lim" : {
    "test_gq" : [0.01, 0.05, 0.1, 0.15, 0.25],
    "test_gdm" : np.logspace(np.log10(0.001),0,101),
    "test_gl" : [0.0]
  },
  "gl_lim" : {
    "test_gq" : [0.01, 0.1, 0.25],
    "test_gdm" : [0.0, 1.0],
    "test_gl" : np.logspace(np.log10(0.001),0,101),
  }
}

# eqn 4.10 https://arxiv.org/pdf/1603.04156.pdf#page12
def calculate_sd(gq, gdm, gl, mMed, mdm) :

    mn = 0.939 # GeV
    val = 2.4e-42 * (gq*gdm/0.25)**2 * (1000./mMed)**4 * (mn*mdm/(mn+mdm))**2
    return val

# eqn 4.3: https://arxiv.org/pdf/1603.04156.pdf#page12
def calculate_si(gq, gdm, gl, mMed, mdm) :

    mn = 0.939 # GeV
    val = 6.9e-41 * (gq*gdm/0.25)**2 * (1000./mMed)**4 * (mn*mdm/(mn+mdm))**2
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

print("Testing point:")
mchi=400.
mmed = 1200.
gq = 0.042
gl = 0.0
gchi = 1.0
print("SI for selected point is",calculate_si(gq, gchi, gl, mmed, mchi))

# Load pickle files with polygons
for collider in ['hl-lhc', 'fcc-hh'] :
    xlow = 1
    xhigh = 2000 if 'hl-lhc' in collider else 4000
    ylow = 1e-46 if 'hl-lhc' in collider else 1e-50
    yhigh = 1e-37 if 'hl-lhc' in collider else 1e-42    
    for model in ['vector','axial'] :        
        print("Starting model",model)
        ylabel = "$\sigma_{SD}$" if 'axial' in model else "$\sigma_{SI}$"
        ylabel = ylabel + " ($\chi$-nucleon) [cm$^2$]" # No difference between proton & neutron for SD unless comparing to other limits        
        with open('{0}_exclusion_contours_couplingDMmass_{1}.pkl'.format(model,collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)

            # Make new dict re-sorted by what coupling the limit is with respect to
            # and convert the contours to DD curves.
            exclusions_dd = {}
            for signature in ['dijet','monojet','dilepton'] :
                signature_contours = loaded_polygons[signature]
                for key,exclusions in signature_contours.items() :
                    if not exclusions : continue
                    tokens = key.split("_")
                    couplingscan = "_".join([tokens[2],tokens[3]])
                    dmhypothesis = "_".join([tokens[4],tokens[5]])
                    couplings = "_".join([tokens[6],tokens[7]])
                    
                    print(tokens,signature)
                    print(dmhypothesis)
                    # Get everything required to convert to DD contour
                    converted_contours = []
                    print(exclusions)
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
                        couplings_separate = coupling_set.split("_")
                        tag = "{0}_{1}_{2}_{3}".format(model,couplingscan,dmhypothesis.split("_")[0],coupling_set)
                        label_line =  "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),masslines[dmhypothesis],transform_coupling(couplings_separate[0]),transform_coupling(couplings_separate[1]))
                        legend_lines = []
                        contours_list = []
                        for name, contours in analysis_list.items() :
                            if contours : 
                                legend_lines.append(name)
                                contours_list.append(contours)
                        drawDDPlot(contours_list, legend_lines, this_tag = tag, plot_path = "plots/directdetection/"+collider, addText = label_line,ylabel=ylabel,is_scaling=True, transluscent=True,xhigh=xhigh, ylow=ylow, yhigh=yhigh)
            
                    # No need for multiple hypotheses on one plot here: we have only one, so let's just do it.
                    # Instead let's overlay other couplings.
                    contours_list_dijet = {}
                    contours_list_monojet = {}
                    contours_list_dilepton = {}
                    contours_list_merged = {}
                    for coupling_set in all_couplingsets :
                        couplings_separate = coupling_set.split("_")
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
                            tag = "{0}_{1}_{2}_{3}_{4}".format(model,couplingscan,dmhypothesis, fixedcoupling, analysis)
                            ordered_couplings = list(dictionary[fixedcoupling].keys())
                            ordered_couplings.sort(reverse=True)
                            for variedcoupling in ordered_couplings :
                                contour = dictionary[fixedcoupling][variedcoupling]
                                overlays.append(contour)
                                legend_lines.append(get_legend_line(variedcoupling))
                            label_line = "{0}\n{1}, {2}\n{3}, {4}".format(("Axial-vector" if 'axial' in model else "Vector"),collider.upper(),analysis.capitalize(),masslines[dmhypothesis], transform_coupling(fixedcoupling))
                            drawDDPlot(overlays, legend_lines, this_tag = tag, plot_path = "plots/directdetection/"+collider, addText = label_line,ylabel=ylabel,is_scaling=True, transluscent=True,xhigh=xhigh, ylow=ylow, yhigh=yhigh)
