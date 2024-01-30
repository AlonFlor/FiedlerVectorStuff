import model_file_handling
import numpy as np
from scipy import sparse
import sys
import pyamg


def construct_Laplacian_for_nodes(nodes):
    '''construct graph Laplacian for vertices'''
    N = len(nodes)
    L = sparse.lil_matrix((N,N), dtype=float)
    for i in np.arange(N):
        node = nodes[i]
        node_sum = 0.
        #adjusted_location = node.location*np.array([0.9,1.,0.95])
        for j in node.connections:
            #neighbor = nodes[j]
            #adjusted_neighbor_location = neighbor.location*np.array([0.9,1.,0.95])
            L[i,j] = -1.#/np.linalg.norm(node.location - neighbor.location)#(adjusted_location - adjusted_neighbor_location)#
            node_sum -= L[i,j]
        L[i,i] = node_sum
    return sparse.csr_matrix(L)

def construct_Laplacian_for_faces(face_data, num_nodes):
    '''construct graph Laplacian for faces'''
    M = len(face_data)
    L = sparse.lil_matrix((M,M), dtype=float)

    #get connections between faces
    face_groups = {}
    for i in np.arange(num_nodes):
        face_groups[i] = []
    for i in np.arange(M):
        face = face_data[i]
        for node_index in face:
            face_groups[node_index].append(i)
    face_connections = {}
    for i in np.arange(M):
        face = face_data[i]
        face_connections[i] = set()
        for node_index in face:
            neighbor_list = face_groups[node_index]
            for neighbor in neighbor_list:
                face_connections[i].add(neighbor)
        face_connections[i].remove(i)

    #build the Laplacian matrix
    for i in np.arange(M):
        face_sum = 0.
        for j in face_connections[i]:
            L[i,j] = -1.
            face_sum -= L[i, j]
        L[i, i] = face_sum
    return sparse.csr_matrix(L)



def pyamg_solve(L, initial_guess=None):
    K=2

    # create the AMG hierarchy
    ml = pyamg.smoothed_aggregation_solver(L)

    # initial approximation to the K eigenvectors
    X = np.random.rand(L.shape[0], K) - 0.5
    X[:,0] = (1./L.shape[0])*np.ones(L.shape[0])

    #X[:,1] = 2.*np.array([i for i in range(L.shape[0])])
    #X[:,1] = (1./L.shape[0])*rankings_estimate
    if initial_guess is not None:
        X[:, 1] = initial_guess

    # preconditioner based on ml
    M = ml.aspreconditioner()

    return sparse.linalg.lobpcg(L, X, M=M, tol=1e-8, largest=False, maxiter=400)

def get_Fiedler_vector_reordering(data, num_nodes=None):
    # construct combinatorial Laplacian
    print()
    print("get Laplacian"+("" if num_nodes==None else " for faces"))
    if num_nodes==None:
        L = construct_Laplacian_for_nodes(data)
    else:
        L = construct_Laplacian_for_faces(data, num_nodes)
    print("extract Fiedler vector"+("" if num_nodes==None else " for faces"))
    #extract Fiedler vector
    W, V = pyamg_solve(L)
    Fiedler_vector = V[:,1]
    print("Fiedler_vector",Fiedler_vector,"\ne-value:",W[1])

    #get ordering from Fiedler vector
    resort_indices = np.argsort(Fiedler_vector)
    resort_ranks = np.zeros_like(Fiedler_vector)
    resort_ranks[resort_indices] = np.arange(len(Fiedler_vector))
    resort_ranks_normed = resort_ranks / resort_ranks.shape[0]
    resort_ranks = resort_ranks.astype(int)

    return resort_ranks, resort_ranks_normed

def get_deterministic_scrambling_reordering(data, region_size):
    number_of_regions = int(np.ceil(len(data) / region_size))

    resort_ranks = []
    for i in np.arange(region_size):
        for j in np.arange(number_of_regions):
            new_rank = region_size*j + i
            if new_rank < len(data):
                resort_ranks.append(new_rank)
            else:
                break

    resort_ranks = np.array(resort_ranks)
    resort_ranks_normed = resort_ranks / resort_ranks.shape[0]
    resort_ranks = resort_ranks.astype(int)

    return resort_ranks, resort_ranks_normed

def get_scrambled_reordering(data):
    resort_ranks = np.array([i for i in np.arange(len(data))])
    np.random.shuffle(resort_ranks)

    resort_ranks_normed = resort_ranks / resort_ranks.shape[0]
    resort_ranks = resort_ranks.astype(int)

    return resort_ranks, resort_ranks_normed

def get_partially_scrambled_reordering(data, chunk_size):
    resort_chunk_ranks = np.array([i for i in np.arange(int(np.ceil(len(data) / chunk_size)))])
    np.random.shuffle(resort_chunk_ranks)

    resort_ranks = []
    for i in np.arange(len(resort_chunk_ranks)):
        for j in np.arange(chunk_size):
            candidate = chunk_size * resort_chunk_ranks[i] + j
            if candidate >= len(data):
                break
            resort_ranks.append(candidate)
    resort_ranks = np.array(resort_ranks)

    resort_ranks_normed = resort_ranks / resort_ranks.shape[0]
    resort_ranks = resort_ranks.astype(int)

    return resort_ranks, resort_ranks_normed


#start code execution here
args = sys.argv[1:]
mode = None
model_folder = None
model_name = None
rearrange_faces = False
command = "Fiedler"
size_of_chunks = 1
for arg in args:
    if arg=="-help" or arg=="--help":
        print("\ncommand line call for this program:\n\tpython3 main.py -name <model file name> -folder <path of folder of model> -command <command>")
        print("Can optionally include the -f tag to reorder the faces.")
        print()
        print("Available commands:\n\tFiedler: reorder the vertices of the model by its Fiedler vector."+
              "\n\tsame: keep the order of the vertices of the model the same.\n\tscramble: randomly rearrange the vertices of the model."+
              "\n\tdet_scramble: scramble the vertices in a deterministic manner.\n")
        print("Can give a -chunk_size <size of chunks> flag, which only applies when command is set to scramble or det_scramble."+
              "It controls the chunk size of the scrambled shape. Recommended for det_scramble to work, since the default value is 1.")
        exit(0)
    if arg=="-folder":
        mode = "folder"
    elif arg=="-name":
        mode = "name"
    elif arg=="-command":
        mode = "command"
    elif arg=="-chunk_size":
        mode = "chunks"
    elif mode=="folder":
        model_folder = arg
        mode = None
    elif mode=="name":
        model_name = arg
        mode = None
    elif mode=="command":
        command = arg
        mode = None
    elif mode=="chunks":
        size_of_chunks = int(arg)
        mode = None
    elif arg=="-f":
        rearrange_faces = True

if model_name==None:
    print("Need a model name")
    exit(1)
if command != "Fiedler" and command != "scramble" and command != "det_scramble" and command !="same":
    print("Invalid command. Must be Fielder, scramble, det_scramble, or same. If left out, default is Fiedler.")
    exit(1)
if (command == "scramble" or command == "det_scramble") and size_of_chunks <= 0:
    print("Invalid size of chunks. Must be a positive integer.")
    exit(1)

nodes, face_data, extra_info = model_file_handling.extract_model(model_name, model_folder)       #extract the model from the file
nodes, face_data = model_file_handling.prune_to_make_fully_connected(nodes, face_data)

resort_ranks = None
resort_ranks_normed = None
print(command+("" if rearrange_faces is False else " + rearrange faces"))
file_prefix = ""
if command=="Fiedler":
    resort_ranks, resort_ranks_normed = get_Fiedler_vector_reordering(nodes)
    file_prefix = "Fiedler_reordered_"
elif command=="scramble":
    if size_of_chunks == 1:
        resort_ranks, resort_ranks_normed = get_scrambled_reordering(nodes)
        file_prefix = "scrambled_"
    else:
        resort_ranks, resort_ranks_normed = get_partially_scrambled_reordering(nodes, size_of_chunks)
        file_prefix = f"scrambled_{size_of_chunks}_node_chunks_"
elif command=="det_scramble":
    resort_ranks, resort_ranks_normed = get_deterministic_scrambling_reordering(nodes, size_of_chunks)
    file_prefix = f"det_scrambled_region_size_{size_of_chunks}_"
elif command=="same":
    resort_ranks = np.array([i for i in np.arange(len(nodes))])
    resort_ranks_normed = resort_ranks / resort_ranks.shape[0]

face_resort_ranks = None
if rearrange_faces:
    file_prefix = "vf_" + file_prefix
    if command=="scramble":
        if size_of_chunks == 1:
            face_resort_ranks, face_resort_ranks_normed = get_scrambled_reordering(face_data)
        else:
            face_resort_ranks, face_resort_ranks_normed = get_partially_scrambled_reordering(face_data, size_of_chunks)
    if command=="det_scramble":
        face_resort_ranks, face_resort_ranks_normed = get_deterministic_scrambling_reordering(face_data, size_of_chunks)
    elif command=="Fiedler":
        face_resort_ranks, face_resort_ranks_normed = get_Fiedler_vector_reordering(face_data, len(nodes))


#print out data
if("." not in model_name[-4:]):
    #Dealing with .node and .ele files
    model_file_handling.write_color_PLY_file("color_gradient_"+model_name+".ply", nodes, face_data, resort_ranks_normed)
    model_file_handling.write_color_PLY_file("contour_"+model_name+".ply", nodes, face_data, resort_ranks_normed, contour=True)
    model_file_handling.write_reordered_node_ele_files(file_prefix+model_name, nodes, face_data, resort_ranks, face_rankings=face_resort_ranks)
else:
    #print color gradient and contour
    model_file_handling.write_color_PLY_file("color_gradient_"+model_name[:-4]+".ply", nodes, face_data, resort_ranks_normed)
    model_file_handling.write_color_PLY_file("contour_"+model_name[:-4]+".ply", nodes, face_data, resort_ranks_normed, contour=True)

    #print reordered model
    model_name_to_check = model_name.lower()[-4:]
    if model_name_to_check == ".ply":
        model_file_handling.write_reordered_PLY_file(file_prefix+model_name, nodes, face_data, resort_ranks, face_rankings=face_resort_ranks)
    elif model_name_to_check == ".veg":
        model_file_handling.write_reordered_VEG_file(file_prefix+model_name, nodes, face_data, resort_ranks, extra_info, face_rankings=face_resort_ranks)
    elif model_name_to_check == ".obj":
        model_file_handling.write_reordered_OBJ_file(file_prefix+model_name, nodes, face_data, resort_ranks, extra_info, face_rankings=face_resort_ranks)

