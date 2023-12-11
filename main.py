import model_file_handling
import numpy as np
from scipy import sparse
import sys
import pyamg

def bfs_traverse(nodes, start_node=None):
    queue = []
    visited = np.zeros((len(nodes)))
    current_score = 0.

    #get the start node furthest from the center
    if start_node is None:
        #find center of model
        center_location = np.zeros((3))
        for node in nodes:
            center_location += node.location
        center_location /= len(nodes)

        #get the start node
        start_node = 0
        max_dist = 0.
        for i in np.arange(len(nodes)):
            node = nodes[i]
            dist = np.linalg.norm(node.location - center_location)
            if dist > max_dist:
                max_dist = dist
                start_node = i

    queue.append(start_node)
    visited[start_node] = 1.
    while queue:
        current = nodes[queue.pop(0)]
        for j in current.connections:
            if visited[j] == 0.:
                queue.append(j)
                visited[j] = 1.
                current_score += 1.
    return visited

def prune_to_make_fully_connected(nodes, face_data):
    '''remove disconnected regions and decrement all indices after them'''

    print("removing disconnected regions and decrementing all indices after the indices of the removed nodes")
    to_delete = []
    decrement_stops = []
    #get unconnected regions and the indices of their nodes
    visited = bfs_traverse(nodes)
    tries = 0
    while np.sum(visited)<0.5*len(nodes):
        #assume at least half of the graph is connected. If this fails after 100 tries, exit
        if tries > 100:
            print("Graph appears to be split in half or into lots of little pieces")
            exit(1)
        visited = bfs_traverse(nodes, np.random.randint(0,len(nodes)))
        tries += 1
    for i in np.arange(len(nodes)):
        node = nodes[i]
        if visited[i] == 0:
            to_delete.append(node)
            decrement_stops.append(i)
    #get the altered indices by decrementing indices after the indices of nodes in unconnected regions
    altered_indices = [i for i in np.arange(len(nodes))]
    for i in np.arange(len(altered_indices)):
        dec_num = 0
        for j in decrement_stops:
            if i>j:
                dec_num +=1
        altered_indices[i] -= dec_num
    #put altered indices in nodes' connections
    for node in nodes:
        for i in np.arange(len(node.connections)):
            node.connections[i] = altered_indices[node.connections[i]]
    #mark faces to in isolated regions to delete
    faces_to_delete = []
    for i in np.arange(len(face_data)):
        delete_this = False
        for j in np.arange(len(face_data[i])):
            if nodes[face_data[i][j]] in to_delete:
                delete_this = True
                break
        if delete_this:
            faces_to_delete.append(face_data[i])
    #put altered indices in face data
    for i in np.arange(len(face_data)):
        for j in np.arange(len(face_data[i])):
            face_data[i][j] = altered_indices[face_data[i][j]]
    #delete unconnected nodes
    for node in to_delete:
        nodes.remove(node)
    #delete faces
    for face in faces_to_delete:
        face_data.remove(face)
    print(f"Done. Removed {len(to_delete)} nodes and {len(faces_to_delete)} faces.")

    return nodes, face_data


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
              "\n\tsame: keep the order of the vertices of the model the same.\n\tscramble: randomly rearrange the vertices of the model.\n")
        print("Can give a -chunk_size <size of chunks> flag, which only applies when command is set to scramble. It controls the chunk size of the scrambled shape.")
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
if command != "Fiedler" and command != "scramble" and command !="same":
    print("Invalid command. Must be Fielder, scramble, or same. If left out, default is Fiedler.")
    exit(1)
if command == "scramble" and size_of_chunks <= 0:
    print("Invalid size of chunks. Must be a positive integer.")
    exit(1)

nodes, face_data, extra_info = model_file_handling.extract_model(model_name, model_folder)       #extract the model from the file
nodes, face_data = prune_to_make_fully_connected(nodes, face_data)

resort_ranks = None
resort_ranks_normed = None
print(command)
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
elif command=="same":
    resort_ranks = np.array([i for i in np.arange(len(nodes))])
    resort_ranks_normed = resort_ranks / resort_ranks.shape[0]

face_resort_ranks = None
if rearrange_faces:
    if command=="scramble":
        if size_of_chunks == 1:
            face_resort_ranks, face_resort_ranks_normed = get_scrambled_reordering(face_data)
        else:
            face_resort_ranks, face_resort_ranks_normed = get_partially_scrambled_reordering(face_data, size_of_chunks)
    elif command=="Fiedler":
        face_resort_ranks, face_resort_ranks_normed = get_Fiedler_vector_reordering(face_data, len(nodes))


#print out data
model_file_handling.write_color_PLY_file("color_gradient_"+model_name[:-4]+".ply", nodes, face_data, resort_ranks_normed)
model_file_handling.write_color_PLY_file("contour_"+model_name[:-4]+".ply", nodes, face_data, resort_ranks_normed, contour=True)

model_name_to_check = model_name.lower()[-4:]

if model_name_to_check == ".ply":
    model_file_handling.write_reordered_PLY_file(file_prefix+model_name, nodes, face_data, resort_ranks, face_rankings=face_resort_ranks)
elif model_name_to_check == ".veg":
    model_file_handling.write_reordered_VEG_file(file_prefix+model_name, nodes, face_data, resort_ranks, extra_info, face_rankings=face_resort_ranks)
elif model_name_to_check == ".obj":
    model_file_handling.write_reordered_OBJ_file(file_prefix+model_name, nodes, face_data, resort_ranks, extra_info, face_rankings=face_resort_ranks)
