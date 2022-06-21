import numpy as np
from shapely.geometry import Polygon as shapely_pol
import pickle
from basic_plotter import *
import ROOT
import itertools as it
import sys, os
sys.path.insert(1, '../inputs/directdetection')
import collect_dd

test_gq = [0.25,0.1,0.05,0.02]#,0.01]
test_gdm = [1.0,0.2,0.1,0.05]
test_gl = [0.1,0.05,0.01,0.0]

sd_proton = collect_dd.get_sd_proton()
sd_neutron = collect_dd.get_sd_neutron()
spin_independent = collect_dd.get_spin_independent()

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

# Collect direct detection contours for comparison plots.
def get_dd_lines(lineinfo) :
    legend_lines = list(lineinfo.keys())
    dd_lines = []
    for name in legend_lines :
        dd_lines.append(lineinfo[name])
    return legend_lines, dd_lines

def make_plots(collider, model, contours, legend_lines, fix_couplings, extra_tag = "") :
    xlow = 1
    xhigh = 2000 if 'hl-lhc' in collider else 4000
    if 'vector' in model : ylow = 1e-48 if 'hl-lhc' in collider else 1e-50
    else : ylow = 1e-46 if 'hl-lhc' in collider else 1e-50
    yhigh = 1e-37 if 'hl-lhc' in collider else 1e-42    
    usepath = "plots/directdetection/"+collider
    formatters = {"gq" : "q", "gdm" : "\chi", "gl" : "l"}
    treat_as_scaling = False
    # At least 2 fixed couplings.
    if len(fix_couplings.keys()) > 2 :
        label_line =  "{0}\n{7}, g$_{5}$={2}\ng$_{4}$={1}, g$_{6}$={3}".format(("Axial-vector" if 'axial' in model else "Vector"),fix_couplings["gq"],fix_couplings["gdm"],fix_couplings["gl"],"q","\chi","l",collider.upper())
        tag_line = model+"_gq{0}_gdm{1}_gl{2}".format(fix_couplings["gq"],fix_couplings["gdm"],fix_couplings["gl"])
    else :
        treat_as_scaling = True
        usecouplings = list(fix_couplings.keys())
        useformats = []
        vals = []
        for coupling in usecouplings :
            useformats.append(formatters[coupling])
            vals.append(fix_couplings[coupling])
        label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}".format(("Axial-vector" if 'axial' in model else "Vector"),vals[0],vals[1],useformats[0],useformats[1],collider.upper())
        tag_line = model+"_{0}{1}_{2}{3}".format(usecouplings[0],vals[0],usecouplings[1],vals[1])
    if extra_tag : tag_line = tag_line + "_" + extra_tag
    # And draw. First, version without DD experiment lines
    # Then draw the plots with DD lines on
    if 'vector' in model :
        formatted_lines = get_dd_lines(spin_independent)
        use_ylabel = "$\sigma_{SI}$ ($\chi$-nucleon) [cm$^2$]"
        drawDDPlot(contours,legend_lines, this_tag = tag_line, plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=treat_as_scaling, transluscent=treat_as_scaling, xhigh=xhigh, ylow=ylow, yhigh=yhigh)
        drawDDPlot(contours,legend_lines, this_tag = tag_line+"_withDD", plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=treat_as_scaling, transluscent=treat_as_scaling, xhigh=xhigh, ylow=ylow, yhigh=yhigh, dd_curves = formatted_lines[1], dd_legendlines = formatted_lines[0])
    else :
        use_ylabel = "$\sigma_{SD}$ ($\chi$-nucleon) [cm$^2$]"
        drawDDPlot(contours,legend_lines, this_tag = tag_line, plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=treat_as_scaling, transluscent=treat_as_scaling, xhigh=xhigh, ylow=ylow, yhigh=yhigh)
        formatted_lines_p = get_dd_lines(sd_proton)
        use_ylabel = "$\sigma_{SD}$ ($\chi$-proton) [cm$^2$]"
        drawDDPlot(contours,legend_lines, this_tag = tag_line+"_withDD_proton", plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=treat_as_scaling, transluscent=treat_as_scaling, xhigh=xhigh, ylow=ylow, yhigh=yhigh, dd_curves = formatted_lines_p[1], dd_legendlines = formatted_lines_p[0])  
        formatted_lines_n = get_dd_lines(sd_neutron)
        use_ylabel = "$\sigma_{SD}$ ($\chi$-neutron) [cm$^2$]"
        drawDDPlot(contours,legend_lines, this_tag = tag_line+"_withDD_neutron", plot_path = usepath, addText = label_line, ylabel=use_ylabel, is_scaling=treat_as_scaling, transluscent=treat_as_scaling, xhigh=xhigh, ylow=ylow, yhigh=yhigh, dd_curves = formatted_lines_n[1], dd_legendlines = formatted_lines_n[0])

# Load pickle files with polygons
for collider in ['hl-lhc', 'fcc-hh'] :
    for model in ['vector','axial'] :

        # Limits with fixed couplings
        with open('{0}_exclusion_contours_{1}.pkl'.format(model,collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)

            # Grid of plots:
            exclusions_dd = {'dijet' : {},'monojet' : {},'dilepton' : {}}
            exclusions_separate_dd = {'dijet' : {},'monojet' : {},'dilepton' : {}}
            for gdm in test_gdm :
                for gq in test_gq :
                    contours_list_couplingscan = []
                    legend_lines_couplingscan = []                
                    for gl in test_gl :
                        contours_list = []
                        legend_lines = []
                        for signature in ['dijet','monojet','dilepton'] :
                            # e.g. no coupling to leptons, skip:
                            if (gq, gdm, gl) not in loaded_polygons[signature].keys() :
                                continue
                            exclusions = loaded_polygons[signature][(gq, gdm, gl)]
                            if not exclusions : continue
                            inner_contours = []
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
                                    if model=='axial' :
                                        sigma = calculate_sd(gq, gdm, gl, x, y)
                                    else :
                                        sigma = calculate_si(gq, gdm, gl, x, y)
                                    contour_DD_raw.append((y, sigma))
                                    raw_graph.AddPoint(y, sigma)
                                contour_DD = shapely_pol(contour_DD_raw)
                                inner_contours.append(contour_DD)
                            contours_list.append(inner_contours)
                            legend_lines.append(signature.capitalize())
                            exclusions_dd[signature][(gq, gdm, gl)] = contours_list
                            exclusions_separate_dd[signature][(gq, gdm, gl)] = inner_contours
                        # First set of plots: 3 contours, one plot for every coupling combo
                        make_plots(collider, model, contours_list, legend_lines, {"gq" : gq, "gdm" : gdm, "gl" : gl})
                        full_polygons = merge_exclusions(contours_list)
                        if not full_polygons : continue 
                        contours_list_couplingscan.append(full_polygons)
                        legend_lines_couplingscan.append("g$_{0}$={1}".format("l",gl))
                    # Second set of plots: merge all contours; fix gq and vary gl.
                    # Note this is not meaningful where we don't have dilepton projections - skip then.
                    make_plots(collider, model, contours_list_couplingscan, legend_lines_couplingscan, {"gq" : gq, "gdm" : gdm})
                    # Could do signature only, fixed gq and varying gl, 
                    # but don't think we need it right now. Would go here.

                # Need second set of plots with gl fixed instead:
                for gl in test_gl :
                    contours_list_couplingscan = []
                    legend_lines_couplingscan = []
                    for gq in test_gq :
                        contours_list = []
                        for signature in ['dijet','monojet','dilepton'] :
                            # e.g. no coupling to leptons, skip:
                            if (gq, gdm, gl) not in exclusions_dd[signature].keys() :
                                continue
                            exclusions = exclusions_dd[signature][(gq, gdm, gl)]
                            contours_list+=exclusions
                        full_polygons = merge_exclusions(contours_list)
                        if not full_polygons : continue
                        contours_list_couplingscan.append(full_polygons)
                        legend_lines_couplingscan.append("g$_{0}$={1}".format("q",gq))
                    make_plots(collider, model, contours_list_couplingscan, legend_lines_couplingscan, {"gl" : gl, "gdm" : gdm})
                
                    # A version overlaying all monojet and overlaying all dijet, but not combining
                    for signature in ['dijet','monojet','dilepton'] :
                        sub_contours_list = []
                        sub_legends_list = []
                        for gq in test_gq :
                            if (gq, gdm, gl) not in exclusions_separate_dd[signature].keys() : continue
                            sub_contours_list.append(exclusions_separate_dd[signature][(gq, gdm, gl)])
                            sub_legends_list.append("g$_{0}$={1}".format("q",gq))
                        make_plots(collider, model, sub_contours_list, sub_legends_list, {"gl" : gl, "gdm" : gdm}, extra_tag = signature)
            # And now third set of plots with gq and gl fixed:
            for gq in test_gq :
                for gl in test_gl :
                    contours_list_couplingscan = []
                    legend_lines_couplingscan = []
                    for gdm in test_gdm :
                        contours_list = []
                        for signature in ['dijet','monojet','dilepton'] :
                            # e.g. no coupling to leptons, skip:
                            if (gq, gdm, gl) not in exclusions_dd[signature].keys() :
                                continue                            
                            exclusions = exclusions_dd[signature][(gq, gdm, gl)]
                            contours_list+=exclusions
                        if all(not i for i in contours_list) : continue
                        full_polygons = merge_exclusions(contours_list)
                        contours_list_couplingscan.append(full_polygons)
                        legend_lines_couplingscan.append("g$_{0}$={1}".format("\chi",gdm))
                    make_plots(collider, model, contours_list_couplingscan, legend_lines_couplingscan, {"gq" : gq, "gl" : gl})
                    # A version overlaying all monojet and overlaying all dijet, but not combining
                    for signature in ['dijet','monojet','dilepton'] :
                        sub_contours_list = []
                        sub_legends_list = []
                        for gdm in test_gdm :
                            if (gq, gdm, gl) not in exclusions_separate_dd[signature].keys() : continue
                            sub_contours_list.append(exclusions_separate_dd[signature][(gq, gdm, gl)])
                            sub_legends_list.append("g$_{0}$={1}".format("\chi",gdm))                   
                        make_plots(collider, model, sub_contours_list,sub_legends_list, {"gq" : gq, "gl" : gl}, extra_tag = signature)