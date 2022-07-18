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

def get_aspect_ratio(ax, isloglog=False, islinearlog=False) :
  ratio = 1.0
  xleft, xright = ax.get_xlim()
  ybottom, ytop = ax.get_ylim()
  if isloglog :
    nUnitsX = np.log10(xright) - np.log10(xleft)
    nUnitsY = np.log10(ytop) - np.log10(ybottom)
    return nUnitsX/nUnitsY
  elif islinearlog :
    nUnitsY = np.log10(ytop) - np.log10(ybottom)
    nUnitsX = (xright-xleft)
    return abs(nUnitsX/nUnitsY)
  else :
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

def get_gradient_fill(axes,contour,use_colour,logx,logy) :
  (xmin, ymin, xmax, ymax) = contour.bounds
  z = np.empty((1, 100, 4), dtype=float)
  rgb = ColorConverter.to_rgb(use_colour)
  # Just basically array of 100 steps in gradient with rgb
  z[:,:,:3] = rgb
  z[:,:,-1] = np.linspace(0, 0.9, 100)[None,:] # Go to 1 for true solid
  # Aspect ratio nonsense
  if logx and logy :
    aspect = get_aspect_ratio(axes,isloglog=True)
  elif logy and not logx :
    aspect = get_aspect_ratio(axes,islinearlog=True)
  else :
    aspect = get_aspect_ratio(axes)
  im = axes.imshow(z, aspect=aspect, extent=[xmin, xmax, ymin, ymax],
                   origin='lower')
  return im

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

def coreDrawFunction(axes,contour_groups,legend_lines,addText = "",is_scaling=False, transluscent=False,xlow=None,xhigh=None, ylow=None,yhigh=None, use_colourscheme=False, dash_contour = [],gradient_fill=[],dashed_lines=[], dashed_legends=[],text_spot = 0,logx=False,logy=False) :

    # Check input
    if not dash_contour : 
      dash_contour = [False for i in contour_groups]
    elif len(dash_contour) != len(contour_groups) :
      print("Error: dash_contour list must have an entry for each contour!")
      exit(1)
    if not gradient_fill :
      gradient_fill = [False for i in contour_groups]
    elif len(gradient_fill) != len(contour_groups) :
      print("Error: gradient_fill list must have an entry for each contour!")    
      exit(1)

    if addText :
        text_x = 0.25 if not text_spot else text_spot[0]
        text_y = 0.82 if not text_spot else text_spot[1]
        text_ploty = text_y - 0.05*addText.count('\n')
        plt.figtext(text_x,text_ploty,addText,size=14)

    # Need cute colours
    if is_scaling :
        ncols = len(contour_groups)
        ncols_dashed = len(dashed_lines)
        # old version: all blues. Great for overlaid but hard to tell apart if transluscent.
        #fill_colours = [scale_lightness('cornflowerblue',0.5+i/(ncols-1) for i in range(ncols)]
        if ncols < 2 :
          colours_raw = [colorFader('cornflowerblue','turquoise',0.5)] # other good option: 'lightgreen'
        else :
          colours_raw = [colorFader('cornflowerblue','turquoise',i/(ncols-1)) for i in range(ncols)]
        # Additional dashed line colours: want contrast if we're doing scaling
        if ncols_dashed < 2 :
          colours_dashed = [colorFader('crimson','gold',0.5)]
        else :
          colours_dashed = [colorFader('crimson','gold',i/(ncols-1)) for i in range(ncols)]
        if transluscent :
          fill_colours = [ColorConverter.to_rgba(col, alpha=0.5) for col in colours_raw]
          line_colours = ['black' for i in fill_colours]
        else :
          line_colours = colours_raw
        line_width = 1
    else :
        if len(contour_groups) < 4 or use_colourscheme:
          colours_raw = ['cornflowerblue','turquoise','mediumorchid']
          colours_dashed = ['royalblue','darkcyan','darkorchid']
        else :
          from matplotlib.pyplot import cm
          colours_raw = cm.rainbow(np.linspace(0, 1, len(contour_groups)))
          colours_dashed = colours_raw
        fillOpacity = 0.5
        fill_colours = [ColorConverter.to_rgba(col, alpha=fillOpacity) for col in colours_raw]
        line_colours = colours_raw
        line_width = 2
    # Going to do a colour scheme for clarity
    for group_index, (contour_group,label_line, dashed, gradient) in enumerate(zip(contour_groups,legend_lines,dash_contour,gradient_fill)) :
        if use_colourscheme :
          if "dijet" in label_line or "Dijet" in label_line : icol = 0
          elif "mono" in label_line or "Mono" in label_line : icol = 1
          else : icol = 2
        else :
          icol = group_index
        for index, contour in enumerate(contour_group) :
          if index == 0 : uselabel = label_line
          else : uselabel = "_"
          if dashed :
            #usefacecolor = 'none'
            uselinecolor = colours_dashed[icol]
            #uselinecolor = line_colours[icol]
            linestyle='dashed'
          else :
            uselinecolor = line_colours[icol]
            linestyle='solid'
          if gradient :
            im = get_gradient_fill(axes,contour,fill_colours[icol],logx,logy) # colours_raw
            usefacecolor = 'none'
            # Need to mess with legend with an un-plotted patch.
            # This is hacky but whatever.
            if index == 0 :
              patch_forleg = Polygon([[-2,0],[-1,0],[-1,-1]], facecolor=fill_colours[icol], edgecolor=uselinecolor, label=uselabel,linewidth=line_width,linestyle=linestyle)
              axes.add_patch(patch_forleg)
              uselabel = "_"
          else :
            usefacecolor = fill_colours[icol]
          # Drawing here
          patch = Polygon(list(contour.exterior.coords), facecolor=usefacecolor, edgecolor=uselinecolor, zorder=2, label=uselabel,linewidth=line_width,linestyle=linestyle) 
          axes.add_patch(patch)
          if gradient : im.set_clip_path(patch)
      
    if dashed_lines :
      for i,(newline,label) in enumerate(zip(dashed_lines,dashed_legends)) :
        #plt.plot(newline[0],newline[1], color=dd_colours[i]) # for an actual line
        for j,line in enumerate(newline) :
          patch = Polygon(list(line.exterior.coords),facecolor='none',edgecolor=colours_dashed[i],label=(label if j==0 else "_"),zorder=2,linewidth=2,linestyle='dashed')
          axes.add_patch(patch)

    return

def drawMassMassPlot(contour_groups, legend_lines, this_tag = "default", plot_path = "plots", addText = "",is_scaling=False, transluscent=False,xlow=None,xhigh=None, ylow=None,yhigh=None, use_colourscheme=False, dash_contour = [],gradient_fill=[],dashed_lines=[],dashed_legends=[],text_spot = 0) :

    # Check input
    if not dash_contour : 
      dash_contour = [False for i in contour_groups]
    elif len(dash_contour) != len(contour_groups) :
      print("Error: dash_contour list must have an entry for each contour!")
      exit(1)
    if not gradient_fill :
      gradient_fill = [False for i in contour_groups]
    elif len(gradient_fill) != len(contour_groups) :
      print("Error: gradient_fill list must have an entry for each contour!")    
      exit(1)

    # Check output
    if not os.path.exists(plot_path) :
        os.makedirs(plot_path)

    # Object for plotting
    plt.rc('font',size=16)
    fig,ax=plt.subplots(1,1)

    usexlow = xlow if xlow else 0
    usexhigh = xhigh if xhigh else 7500
    useylow = ylow if ylow else 0
    useyhigh = yhigh if yhigh else 1200
    ax.set_xlim(usexlow, usexhigh)
    ax.set_ylim(useylow, useyhigh)
    ratio = get_aspect_ratio(ax)
    ax.set_aspect(ratio)

    ax.set_xlabel(r"m$_{\rm med}$ [GeV]")
    ax.set_ylabel("m$_{\chi}$ [GeV]")   

    coreDrawFunction(ax,contour_groups,legend_lines,addText,is_scaling, transluscent,xlow,xhigh, ylow,yhigh, use_colourscheme, dash_contour,gradient_fill,dashed_lines,dashed_legends,text_spot)
    ax.legend(fontsize=14,loc='upper right')

    plt.savefig(plot_path+'/{0}.pdf'.format(this_tag),bbox_inches='tight')

    plt.close(fig)    

def drawCouplingMassPlot(contour_groups, legend_lines, this_tag = "default", plot_path = "plots",xlabel="",ylabel="coupling", addText = "",is_scaling=False, transluscent=False,xlow=None, xhigh=None, ylow=None,yhigh=None, use_colourscheme=False, dash_contour = [], gradient_fill=[],dashed_lines=[], dashed_legends=[],text_spot = 0,legend_spot='upper right') :

    # Check input
    if not dash_contour : 
      dash_contour = [False for i in contour_groups]
    elif len(dash_contour) != len(contour_groups) :
      print("Error: dash_contour list must have an entry for each contour!")
      exit(1)
    if not gradient_fill :
      gradient_fill = [False for i in contour_groups]
    elif len(gradient_fill) != len(contour_groups) :
      print("Error: gradient_fill list must have an entry for each contour!")    
      exit(1)

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
    ratio = get_aspect_ratio(ax, islinearlog=True)
    ax.set_aspect(ratio)

    if not xlabel :
      ax.set_xlabel(r"m$_{\rm med}$ [GeV]")
    else :
      ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)   

    if 'outside' in legend_spot : text_spot = [0.84, 0.83]

    coreDrawFunction(ax,contour_groups,legend_lines,addText,is_scaling, transluscent,xlow,xhigh, ylow,yhigh, use_colourscheme, dash_contour,gradient_fill,dashed_lines,dashed_legends,text_spot,logx=False,logy=True)

    if 'outside' in legend_spot :
      leg_y = text_spot[1]-0.05*(addText.count('\n')-2)
      leg = ax.legend(fontsize=14,bbox_to_anchor=(1.02,leg_y),loc="upper left")
      leg.get_frame().set_linewidth(0.0)
    else :
      ax.legend(fontsize=14,loc=legend_spot)

    plt.savefig(plot_path+'/{0}.pdf'.format(this_tag),bbox_inches='tight')

    plt.close(fig)        

def drawDDPlot(contour_groups, legend_lines, this_tag = "default", plot_path = "plots", addText = "",ylabel="\sigma",is_scaling=False, transluscent=False, xlow=None, xhigh=None, ylow=None,yhigh=None, dd_curves = None, dd_legendlines = None,dashed_lines=[], dashed_legends=[]) :

    # Check output
    if not os.path.exists(plot_path) :
        os.makedirs(plot_path)

    # Object for plotting
    fig,ax=plt.subplots(1,1)

    usexhigh = xhigh if xhigh else 2000
    useylow = ylow if ylow else 1e-46
    useyhigh = yhigh if yhigh else 1e-37
    ax.set_xscale('log')
    ax.set_yscale('log')    
    ax.set_xlim(1, usexhigh)
    ax.set_ylim(useylow,useyhigh)
    ratio = get_aspect_ratio(ax, isloglog=True)
    ax.set_aspect(ratio)
    plt.rc('font',size=16)

    ax.set_xlabel("m$_{\chi}$ [GeV]", fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    text_spot = [0.84, 0.83]
    coreDrawFunction(ax,contour_groups,legend_lines,addText,is_scaling, transluscent,xlow,xhigh, ylow,yhigh, use_colourscheme=False, dash_contour=[],gradient_fill=[],dashed_lines=dashed_lines,dashed_legends=dashed_legends,text_spot=text_spot,logx=True,logy=True)

    leg_y = text_spot[1]-0.05*(addText.count('\n')-2)
    leg = ax.legend(fontsize=14,bbox_to_anchor=(1.02,leg_y),loc="upper left")
    leg.get_frame().set_linewidth(0.0)

    plt.savefig(plot_path+'/{0}.pdf'.format(this_tag),bbox_inches='tight')

    plt.close(fig)    
            

