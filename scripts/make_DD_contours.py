import numpy as np
from shapely.geometry import Polygon as shapely_pol
import pickle
from basic_plotter import *
import ROOT
import itertools as it

test_gq = [0.25,0.15,0.1,0.05,0.01]
test_gl = [0.1,0.05,0.01,0.0]

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

# Load pickle files with polygons
for collider in ['hl-lhc', 'fcc-hh'] :
    for model in ['vector','axial'] :
        with open('{0}_exclusion_contours_{1}.pkl'.format(model,collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)

            # Grid of plots:
            exclusions_dd = {'dijet' : {},'monojet' : {},'dilepton' : {}}
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
                        print("Starting coupling",gq,1.0,gl)
                        exclusions = loaded_polygons[signature][(gq, 1.0, gl)]
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
                                    sigma = calculate_sd(gq, 1.0, gl, x, y)
                                else :
                                    sigma = calculate_si(gq, 1.0, gl, x, y)
                                contour_DD_raw.append((y, sigma))
                                raw_graph.AddPoint(y, sigma)
                            contour_DD = shapely_pol(contour_DD_raw)
                            inner_contours.append(contour_DD)
                        contours_list.append(inner_contours)
                        legend_lines.append(signature)
                        exclusions_dd[signature][(gq, 1.0, gl)] = contours_list
                    # First set of plots: 3 contours, one plot for every coupling combo
                    label_line =  "{0}\n{7}, g$_{5}$={2}\ng$_{4}$={1}, g$_{6}$={3}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,1.0,gl,"q","\chi","l",collider.upper())
                    print(contours_list)
                    drawDDPlot(contours_list, legend_lines, this_tag = model+"_gq{0}_gdm1.0_gl{1}".format(gq, gl), plot_path = "plots/directdetection/"+collider, addText=label_line, xhigh=2000, yhigh=1e-37)
                    full_polygons = merge_exclusions(contours_list)
                    contours_list_couplingscan.append(full_polygons)
                    legend_lines_couplingscan.append("g$_{0}$={1}".format("l",gl))
                # Second set of plots: merge all contours; fix gq and vary gl.
                # Note this is not meaningful where we don't have dilepton projections - skip then.
                label_line = "{0}, {5}\ng$_{3}$={1}, g$_{4}$={2}".format(("Axial-vector" if 'axial' in model else "Vector"),gq,1.0,"q","\chi",collider.upper())
                drawDDPlot(contours_list_couplingscan,legend_lines_couplingscan, this_tag = model+"_gq{0}_gdm1.0".format(gq), plot_path = "plots/massmass/"+collider, addText = label_line,is_scaling=True, xhigh=2000, yhigh=1e-37)                        