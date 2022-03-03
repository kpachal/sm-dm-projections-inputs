import os
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches

# Object for plotting
fig,ax=plt.subplots()

def get_aspect_ratio(ax) :
  ratio = 1.0
  xleft, xright = ax.get_xlim()
  ybottom, ytop = ax.get_ylim()
  return abs((xright-xleft)/(ybottom-ytop))*ratio

def getContours(xvals, yvals, zvals) :
  #plt.figure(1).clf()
  #fig,ax=plt.subplots(1,1)
  excontour = ax.tricontour(xvals, yvals, zvals, levels=[1])
  contour_list = []
  for path in excontour.collections[0].get_paths():
    this_contour = [path.vertices[:,0],path.vertices[:,1]]
    contour_list.append(this_contour)
  plt.show()
  #plt.cla()
  return contour_list

def drawContourPlot(grid_list, addPoints = False, this_tag = "default", plot_path = "plots", addText = "") :

  # Check output
  if not os.path.exists(plot_path) :
    os.makedirs(plot_path)

  # Object for plotting
  #fig,ax=plt.subplots(1,1)
  plt.cla()

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
    print("Doing grid", grid)
    cp = ax.tricontourf(grid[0], grid[1], grid[2], levels=levels, cmap='Blues_r')

    #plt.show()

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

  #fig.colorbar(cp)

  plt.savefig(plot_path+'/massmass_{0}.eps'.format(this_tag),bbox_inches='tight')
  plt.savefig(plot_path+'/massmass_{0}.pdf'.format(this_tag),bbox_inches='tight')

  plt.cla()

  return contour_list
