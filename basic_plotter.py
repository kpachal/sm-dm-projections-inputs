import os
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import numpy as np

from shapely.geometry import Polygon
from shapely.ops import cascaded_union
from shapely import validation
from matplotlib.path import Path
from matplotlib.patches import Polygon

def get_aspect_ratio(ax) :
  ratio = 1.0
  xleft, xright = ax.get_xlim()
  ybottom, ytop = ax.get_ylim()
  return abs((xright-xleft)/(ybottom-ytop))*ratio

# def getContours(xvals, yvals, zvals) :
#   fig,ax=plt.subplots()
#   excontour = ax.tricontour(xvals, yvals, zvals, levels=[1])
#   contour_list = []
#   for path in excontour.collections[0].get_paths():
#     this_contour = [path.vertices[:,0],path.vertices[:,1]]
#     contour_list.append(this_contour)
#   plt.close(fig)
#   return contour_list

def get_contours(xvals, yvals, zvals) :

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
          new_shape = Polygon(sub_segment)

          # Check validity, fail if not valid
          if not new_shape.is_valid :
            print("Invalid polygon from extracted contour!")
            print(validation.explain_validity(new_shape))
            exit(1)

          sub_contours.append(new_shape)

    contours[level] = sub_contours
  
  return contours  

def merge_exclusions(contour_list) :

    # If there's only one, do nothing.
    if len(contour_list) < 2 :
      return contour_list

    # No longer matters which input each contour is from:
    # flatten the input to a simple list of polygons
    contour_list = [item for sublist in contour_list for item in sublist]

    # Merge
    merged = cascaded_union(contour_list)

    # Split multipolygon back into list of polygons, if needed
    if merged.geom_type == 'MultiPolygon':
      poly_list = [poly for poly in merged]
    else :
      poly_list = [merged]

    return poly_list 

def drawContourPlotRough(grid_list, addPoints = False, this_tag = "default", plot_path = "plots", addText = "") :

  # Check output
  if not os.path.exists(plot_path) :
    os.makedirs(plot_path)

  # Object for plotting
  fig,ax=plt.subplots(1,1)

  levels = range(26)  # Levels must be increasing.
  ax.set_xlim(0, 7500)
  ax.set_ylim(0, 1200)
  plt.rc('font',size=17)
  ratio = get_aspect_ratio(ax)
  ax.set_aspect(ratio)

  ax.set_xlabel("m$_{ZA}$ [GeV]")
  ax.set_ylabel("m$_{\chi}$ [GeV]")   

  # Add text
  if addText :
    plt.figtext(0.2,0.75,addText,size=14)

  contour_list = []
  for grid in grid_list :
    cp = ax.tricontourf(grid[0], grid[1], grid[2], levels=levels, cmap='Blues_r')

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

  fig.colorbar(cp)

  plt.savefig(plot_path+'/massmass_{0}.eps'.format(this_tag),bbox_inches='tight')
  plt.savefig(plot_path+'/massmass_{0}.pdf'.format(this_tag),bbox_inches='tight')

  plt.close(fig)

  return contour_list

def drawMassMassPlot(contour_groups, legend_lines, this_tag = "default", plot_path = "plots", addText = "") :

    # Check output
    if not os.path.exists(plot_path) :
        os.makedirs(plot_path)

    # Object for plotting
    fig,ax=plt.subplots(1,1)

    ax.set_xlim(0, 7500)
    ax.set_ylim(0, 1200)
    plt.rc('font',size=16)
    ratio = get_aspect_ratio(ax)
    ax.set_aspect(ratio)

    ax.set_xlabel("m$_{ZA}$ [GeV]")
    ax.set_ylabel("m$_{\chi}$ [GeV]")   

    if addText :
        plt.figtext(0.23,0.77,addText,size=14)

    # Need 3 cute colours
    colours = ['cornflowerblue','turquoise','mediumorchid']    
    for contour_group,label_line,col in zip(contour_groups,legend_lines,colours) :
        for contour in contour_group :
            patch = Polygon(list(contour.exterior.coords), facecolor=col, edgecolor=col, alpha=0.5, zorder=2, label=label_line)
            ax.add_patch(patch)
    ax.legend(fontsize=14)

    plt.savefig(plot_path+'/massmass_{0}.eps'.format(this_tag),bbox_inches='tight')
    plt.savefig(plot_path+'/massmass_{0}.pdf'.format(this_tag),bbox_inches='tight')

    plt.close(fig)    

            

