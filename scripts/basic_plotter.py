import os
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.colors as cols
import numpy as np

from shapely.geometry import Polygon as shapely_pol
from shapely.geometry import MultiPoint
from shapely.ops import unary_union
from shapely import validation
from matplotlib.path import Path
from matplotlib.patches import Polygon
from matplotlib.colors import ColorConverter,ListedColormap, LinearSegmentedColormap,rgb2hex
from matplotlib import cm

import colorsys

# This is just how many points we get - should be enough for anything
colourmap = cm.get_cmap('Blues_r', 200)

def get_aspect_ratio(ax) :
  ratio = 1.0
  xleft, xright = ax.get_xlim()
  ybottom, ytop = ax.get_ylim()
  return abs((xright-xleft)/(ybottom-ytop))*ratio

def scale_lightness(matplotlib_col, scale_l):
    rgb = ColorConverter.to_rgb(matplotlib_col)
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(ColorConverter.to_rgb(c1))
    c2=np.array(ColorConverter.to_rgb(c2))
    return cols.to_hex((1-mix)*c1 + mix*c2)

def get_colours(npoints):
    # Ignore the last 20% of the scale to keep from approaching white too closely
    breaks = [float(i)/(npoints*1.2) for i in range(npoints)]
    colours = [rgb2hex(colourmap(bb)) for bb in breaks]
    return colours

def get_contours(xvals, yvals, zvals) :

  # Don't try to plot if it's empty
  if xvals.size == 0 : return {0 : None}

  # This separates into two sets of filled spaces:
  # spaces under 1 (that is, between 1 and 0)
  # and spaces over 1.
  cntr = plt.tricontourf(xvals, yvals, zvals, levels=[0,1.0,1e5], cmap="RdBu_r")

  # Extract them and organise them for returning.
  contours = {}

  # Loop through all polygons that have a specific intensity level
  for level,col in enumerate(cntr.collections) :

    # There can be multiple contours at each level
    sub_contours = []
    for contour_path in col.get_paths() :

        # Each path can also contain what should be
        # mulitple contours: split by MOVETO codes
        # where necessary
        vertices = contour_path.vertices
        codes = contour_path.codes
        idx = np.where(codes==Path.MOVETO)[0]
        vertex_segs = np.split(vertices,idx)[1:]

        # Now create and save polygons!
        for sub_segment in vertex_segs :
          #multi_shape = MultiPoint(sub_segment)          
          new_shape = shapely_pol(sub_segment)

          # Check validity, fail if not valid
          if not new_shape.is_valid :
            print("Invalid polygon from extracted contour! Buffering....")
            print(sub_segment)
            bval=0.01
            while (not new_shape.is_valid) and bval < 0.1 :
              new_shape = new_shape.buffer(bval) # smoothing out invalidities?
              bval+=0.1
              if new_shape.is_valid : break
            
            if not new_shape.is_valid :
              print("Not fixed by buffering.")
              print(validation.explain_validity(new_shape))
              exit(1)

          sub_contours.append(new_shape)

    contours[level] = sub_contours
  
  return contours  

def merge_exclusions(contour_list) :

    # If there's only one, do nothing.
    if len(contour_list) == 1 :
      return contour_list[0]
    elif len(contour_list) == 0 :
      return None

    # No longer matters which input each contour is from:
    # flatten the input to a simple list of polygons
    contour_list = [item for sublist in contour_list for item in sublist]

    # Merge
    merged = unary_union(contour_list)

    # Split multipolygon back into list of polygons, if needed
    poly_list = []
    if merged.geom_type == 'MultiPolygon':
      poly_list = [poly for poly in merged.geoms]
    else :
      poly_list = [merged]

    return poly_list 

def drawContourPlotRough(grid_list, addPoints = False, this_tag = "default", plot_path = "plots", addText = "", xlow = None, xhigh=None, ylow = None, yhigh=None, vsCoupling=False) :

  # Check output
  if not os.path.exists(plot_path) :
    os.makedirs(plot_path)

  # Object for plotting
  fig,ax=plt.subplots(1,1)

  #levels = range(26)  # Levels must be increasing.
  levels = [0.001,0.005,0.01,0.05,0.1,0.5,1,5,10]
  usexlow = xlow if xlow else 0
  usexhigh = xhigh if xhigh else 7500
  useylow = ylow if ylow else 0
  useyhigh = yhigh if yhigh else 1200
  ax.set_xlim(xlow, usexhigh)
  ax.set_ylim(ylow, useyhigh)
  plt.rc('font',size=17)
  ratio = get_aspect_ratio(ax)
  ax.set_aspect(ratio)

  ax.set_xlabel(r"m$_{\rm med}$ [GeV]")
  if vsCoupling :
    ax.set_ylabel("Coupling")
    #ax.set_yscale('log')
  else :
    ax.set_ylabel("m$_{\chi}$ [GeV]")

  # Add text
  if addText :
    plt.figtext(0.2,0.75,addText,size=14)

  contour_list = []
  plot_empty = True
  for grid in grid_list :
    # Don't try to plot if it's empty
    if grid[0].size == 0 : continue
    else : plot_empty = False

    colors = get_colours(len(levels))
    cp = ax.tricontourf(grid[0], grid[1], grid[2], levels=levels, colors=colors)

    # Want points under contour, if adding them.
    if addPoints :
      # Separate into two populations: excluded and non excluded.
      xexcl,yexcl = [],[]
      xnon,ynon = [],[]
      for x,y,z in zip(grid[0], grid[1], grid[2]) :
        if z < 1. : 
          xexcl.append(x)
          yexcl.append(y)
        else :
          xnon.append(x)
          ynon.append(y)
      ax.scatter(xnon,ynon,color='red', marker='o',facecolors='none',linewidths=2)
      ax.scatter(xexcl,yexcl,color='white', marker='o',facecolors='none',linewidths=2)

    # Now add exclusion contour (if not doing official - harder to see with both)
    excontour = ax.tricontour(grid[0], grid[1], grid[2],levels=[1],colors=['r'],linewidths=[2])
    for path in excontour.collections[0].get_paths():
        this_contour = [path.vertices[:,0],path.vertices[:,1]]
        contour_list.append(this_contour)

  if plot_empty : return

  fig.colorbar(cp)

  if not vsCoupling :
    print("Making plot",plot_path+'/massmass_{0}.pdf'.format(this_tag))
    plt.savefig(plot_path+'/massmass_{0}.pdf'.format(this_tag),bbox_inches='tight')
  else :
    print("Making plot",plot_path+'/couplingmass_{0}.pdf'.format(this_tag))
    plt.savefig(plot_path+'/couplingmass_{0}.pdf'.format(this_tag),bbox_inches='tight')  

  plt.close(fig)

  return contour_list

def drawMassMassPlot(contour_groups, legend_lines, this_tag = "default", plot_path = "plots", addText = "",is_scaling=False, transluscent=False,xlow=None,xhigh=None, ylow=None,yhigh=None) :

    # Check output
    if not os.path.exists(plot_path) :
        os.makedirs(plot_path)

    # Object for plotting
    fig,ax=plt.subplots(1,1)

    usexlow = xlow if xlow else 0
    usexhigh = xhigh if xhigh else 7500
    useylow = ylow if ylow else 0
    useyhigh = yhigh if yhigh else 1200
    ax.set_xlim(usexlow, usexhigh)
    ax.set_ylim(useylow, useyhigh)
    plt.rc('font',size=16)
    ratio = get_aspect_ratio(ax)
    ax.set_aspect(ratio)

    ax.set_xlabel(r"m$_{\rm med}$ [GeV]")
    ax.set_ylabel("m$_{\chi}$ [GeV]")   

    if addText :
        if addText.count('\n') == 1 :
            plt.figtext(0.23,0.77,addText,size=14)
        elif addText.count('\n') == 2 :
            plt.figtext(0.23,0.72,addText,size=14)
        else : plt.figtext(0.23,0.82,addText,size=14)

    # Need 3 cute colours
    if is_scaling :
        ncols = len(contour_groups)
        # old version: all blues. Great for overlaid but hard to tell apart if transluscent.
        #fill_colours = [scale_lightness('cornflowerblue',0.5+i/(ncols-1) for i in range(ncols)]
        if ncols < 2 :
          fill_colours = [colorFader('cornflowerblue','turquoise',0.5)] # other good option: 'lightgreen'
        else :
          fill_colours = [colorFader('cornflowerblue','turquoise',i/(ncols-1)) for i in range(ncols)]
        if transluscent :
          fill_colours = [ColorConverter.to_rgba(col, alpha=0.5) for col in fill_colours]
          line_colours = ['black' for i in fill_colours]
        else :
          line_colours = fill_colours          
        line_width = 1
    else :
        colours_raw = ['cornflowerblue','turquoise','mediumorchid']
        fillOpacity = 0.5
        fill_colours = [ColorConverter.to_rgba(col, alpha=fillOpacity) for col in colours_raw]
        line_colours = colours_raw
        line_width = 2
    for contour_group,label_line,face_col,line_col in zip(contour_groups,legend_lines,fill_colours,line_colours) :
        for index, contour in enumerate(contour_group) :
            if index == 0 :
                patch = Polygon(list(contour.exterior.coords), facecolor=face_col, edgecolor=line_col, zorder=2, label=label_line,linewidth=line_width) #alpha=fillOpacity, 
            else :
                patch = Polygon(list(contour.exterior.coords), facecolor=face_col, edgecolor=line_col, zorder=2, label="_",linewidth=line_width) # alpha=fillOpacity,
            ax.add_patch(patch)
    ax.legend(fontsize=14,loc='upper right')

    #plt.savefig(plot_path+'/massmass_{0}.eps'.format(this_tag),bbox_inches='tight')
    plt.savefig(plot_path+'/massmass_{0}.pdf'.format(this_tag),bbox_inches='tight')

    plt.close(fig)    

def drawCouplingMassPlot(contour_groups, legend_lines, this_tag = "default", plot_path = "plots",xlabel="",ylabel="coupling", addText = "",is_scaling=False, transluscent=False,xlow=None, xhigh=None, ylow=None,yhigh=None, use_colourscheme=False) :

    # Check output
    if not os.path.exists(plot_path) :
        os.makedirs(plot_path)

    # Object for plotting
    fig,ax=plt.subplots(1,1)

    usexlow = xlow if xlow else 0
    usexhigh = xhigh if xhigh else 7500
    useylow = ylow if ylow else 0.003
    useyhigh = yhigh if yhigh else 0.5
    ax.set_xlim(usexlow, usexhigh)
    ax.set_ylim(useylow, useyhigh)

    ax.set_yscale('log') 
    plt.rc('font',size=16)
    #ratio = get_aspect_ratio(ax)
    #ax.set_aspect(ratio)

    if not xlabel :
      ax.set_xlabel(r"m$_{\rm med}$ [GeV]")
    else :
      ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)   

    if addText :
        if addText.count('\n') == 1 :
            plt.figtext(0.23,0.77,addText,size=14)
        elif addText.count('\n') == 2 :
            plt.figtext(0.23,0.72,addText,size=14)
        else : plt.figtext(0.23,0.82,addText,size=14)

    # Need 3 cute colours
    if is_scaling :
        ncols = len(contour_groups)
        if ncols < 2 :
          fill_colours = [colorFader('cornflowerblue','turquoise',0.5)] # other good option: 'lightgreen'
        else :
          fill_colours = [colorFader('cornflowerblue','turquoise',i/(ncols-1)) for i in range(ncols)]
        if transluscent :
          fill_colours = [ColorConverter.to_rgba(col, alpha=0.5) for col in fill_colours]
          line_colours = ['black' for i in fill_colours]
        else :
          line_colours = fill_colours          
        line_width = 1
    else :
        colours_raw = ['cornflowerblue','turquoise','mediumorchid']
        fillOpacity = 0.5
        fill_colours = [ColorConverter.to_rgba(col, alpha=fillOpacity) for col in colours_raw]
        line_colours = colours_raw
        line_width = 2
    # Going to do a colour scheme for clarity
    for group_index, (contour_group,label_line) in enumerate(zip(contour_groups,legend_lines)) :
        if use_colourscheme :
          if "dijet" in label_line : icol = 0
          elif "mono" in label_line : icol = 1
          else : icol = 2
        else :
          icol = group_index
        for index, contour in enumerate(contour_group) :
            if index == 0 :
                patch = Polygon(list(contour.exterior.coords), facecolor=fill_colours[icol], edgecolor=line_colours[icol], zorder=2, label=label_line,linewidth=line_width) #alpha=fillOpacity, 
            else :
                patch = Polygon(list(contour.exterior.coords), facecolor=fill_colours[icol], edgecolor=line_colours[icol], zorder=2, label="_",linewidth=line_width) # alpha=fillOpacity,
            ax.add_patch(patch)
    ax.legend(fontsize=14,loc='upper right')

    print("Creating plot",plot_path+'/couplingmass_{0}.pdf'.format(this_tag))
    plt.savefig(plot_path+'/couplingmass_{0}.pdf'.format(this_tag),bbox_inches='tight')

    plt.close(fig)        

def drawDDPlot(contour_groups, legend_lines, this_tag = "default", plot_path = "plots", addText = "",ylabel="\sigma",is_scaling=False, transluscent=False, xlow=None, xhigh=None, ylow=None,yhigh=None) :

    # Check output
    if not os.path.exists(plot_path) :
        os.makedirs(plot_path)

    # Object for plotting
    fig,ax=plt.subplots(1,1)

    usexhigh = xhigh if xhigh else 2000
    useylow = ylow if ylow else 1e-46
    useyhigh = yhigh if yhigh else 1e-37
    ax.set_xlim(1, usexhigh)
    ax.set_ylim(useylow,useyhigh)
    plt.rc('font',size=16)
    ax.set_xscale('log')
    ax.set_yscale('log')
    #plt.rc('font',size=16)

    ax.set_xlabel("m$_{\chi}$ [GeV]", fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    if addText :
        # For low corner: 0.15/0.19/0.23
        if addText.count('\n') == 1 :
            plt.figtext(0.16,0.77,addText,size=14)
        elif addText.count('\n') == 2 :
            plt.figtext(0.16,0.72,addText,size=14)
        else : plt.figtext(0.16,0.82,addText,size=14)

    # Need 3+ cute colours
    if is_scaling :
        ncols = len(contour_groups)
        # old version: all blues. Great for overlaid but hard to tell apart if transluscent.
        #fill_colours = [scale_lightness('cornflowerblue',0.5+i/(ncols-1) for i in range(ncols)]
        if ncols < 2 :
          fill_colours = [colorFader('cornflowerblue','turquoise',0.5)] # other good option: 'lightgreen'
        else :
          fill_colours = [colorFader('cornflowerblue','turquoise',i/(ncols-1)) for i in range(ncols)]
        if transluscent :
          fill_colours = [ColorConverter.to_rgba(col, alpha=0.5) for col in fill_colours]
          line_colours = ['black' for i in fill_colours]
        else :
          line_colours = fill_colours          
        line_width = 1
    else :
        colours_raw = ['cornflowerblue','turquoise','mediumorchid']
        fillOpacity = 0.5
        fill_colours = [ColorConverter.to_rgba(col, alpha=fillOpacity) for col in colours_raw]
        line_colours = colours_raw
        line_width = 2
    for contour_group,label_line,face_col,line_col in zip(contour_groups,legend_lines,fill_colours,line_colours) :
        for index, contour in enumerate(contour_group) :
            if len(list(contour.exterior.coords)) == 0 : continue
            if index == 0 :
                patch = Polygon(list(contour.exterior.coords), facecolor=face_col, edgecolor=line_col, zorder=2, label=label_line,linewidth=line_width) 
            else :
                patch = Polygon(list(contour.exterior.coords), facecolor=face_col, edgecolor=line_col, zorder=2, label="_",linewidth=line_width)
            ax.add_patch(patch)
    ax.legend(fontsize=14,loc='upper right')

    plt.savefig(plot_path+'/directdetection_{0}.pdf'.format(this_tag),bbox_inches='tight')

    plt.close(fig)    
            

