import model_file_handling
import sys
import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

#get the model
args = sys.argv[1:]
mode = None
for arg in args:
    if arg=="-help" or arg=="--help":
        print("\ncommand line call for this program:\n\tpython3 edge_span.py -name <model file name> -folder <path of folder of model>")
        
    if arg=="-folder":
        mode = "folder"
    elif arg=="-name":
        mode = "name"
    elif mode=="folder":
        model_folder = arg
        mode = None
    elif mode=="name":
        model_name = arg
        mode = None

reordering_names = ["Fiedler_reordered_", "", "scrambled3_"]
#reordering_names = ["", "Fiedler_reordered_", "scrambled_","scrambled2_","scrambled3_"]
colors = ["blue", "black", "red"]
plt.rcParams.update({'font.size': 22})
plt.xscale("log")
plt.yscale("log")


bin_width = 50
print(bin_width)

histogram_string = ""
for r_index, reordering_name in enumerate(reordering_names):
    #extract the info
    name = reordering_name + model_name
    print(name)
    nodes, face_data, extra_info = model_file_handling.extract_model(name, model_folder)       #extract the model from the file
    nodes, face_data = model_file_handling.prune_to_make_fully_connected(nodes, face_data)

    #get the edge spans
    edge_spans = []
    for i,node in enumerate(nodes):
        #idea: add an edge to edge spans if the other index > the current index.
        # 4 Node example:
        #    1:    1-2, 1-3, 1-4
        #    2:    2-3, 2-4
        #    3:    3-4
        for other_node in node.connections:
            potential_span = other_node - i
            if potential_span > 0:
                edge_spans.append(potential_span)
    print("number of edges", len(edge_spans))
        
    #make a histogram
    histogram_x = []
    histogram = []
    max_bin_index = int(max(edge_spans)/bin_width) + 1
    for i in range(max_bin_index):
        histogram_x.append(bin_width*i)
        histogram.append(0)
    print("len histogram", len(histogram))
    for edge_span in edge_spans:
        bin_index_to_use = int(edge_span/bin_width)
        #histogram[min(max_bin_index-1, bin_index_to_use)] += 1
        #if bin_index_to_use < max_bin_index:
        histogram[bin_index_to_use] += 1

    plt.plot(histogram_x, histogram, label=name, linewidth=2.5, color=colors[r_index])
    
    edge_spans_array = np.array(edge_spans)
    histogram_string += name+": "+str(np.mean(edge_spans_array))+" +/- "+str(np.std(edge_spans_array))+"\n"
print(histogram_string)
#ax = plt.gca()
#axins = zoomed_inset_axes(ax, zoom=0.25, loc='upper right')
plt.xlabel("Edge span")
plt.ylabel("Number of edges")
plt.legend()
fig = plt.gcf()
fig.set_figheight(9)
fig.set_figwidth(12)
fig.set_dpi(120)
#plt.show()
fig.tight_layout()
plt.savefig(model_name)#, bbox_inches='tight')

