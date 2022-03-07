import pickle
import ROOT

test_gq = [0.25,0.15,0.1,0.05,0.01]
test_gl = [0.1,0.05,0.01,0.0]

# Load pickle files with polygons
for collider in ['hl-lhc'] : # 'fcc-hh'
    for model in ['vector','axial'] :
        with open('{0}_exclusion_contours_{1}.pkl'.format(model,collider), "rb") as poly_file:
            loaded_polygons = pickle.load(poly_file)

            outfile_name = '{0}_exclusion_contours_{1}.root'.format(model,collider)
            outfile = ROOT.TFile.Open(outfile_name,"RECREATE")
            outfile.cd()

            # Grid of plots: 
            for gq in test_gq :
                for gl in test_gl :
                    for signature in ['dijet','monojet','dilepton'] :
                        # e.g. no coupling to leptons, skip:
                        if (gq, 1.0, gl) not in loaded_polygons[signature].keys() :
                            continue
                        print("Starting coupling",gq,1.0,gl)
                        exclusions = loaded_polygons[signature][(gq, 1.0, gl)]
                        contours_list = []
                        for contour in exclusions :
                          newgraph = ROOT.TGraph()
                          for x, y in list(contour.exterior.coords) :
                            newgraph.AddPoint(x,y)
                          newgraph.Print()
                          contours_list.append(newgraph)
                        if len(contours_list)==0 : continue
                        base_name = "contour_{0}_gq{1}_gdm1.0_gl{2}".format(signature,gq,gl)
                        base_name = base_name.replace(".","p")
                        if len(contours_list)== 1 :
                          contours_list[0].Write(base_name)
                        else :
                          for i,igraph in enumerate(contours_list) :
                            igraph_name = base_name+"_{0}".format(i)
                            igraph.Write(igraph_name)
            
            outfile.Close()

                          