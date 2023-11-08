import numpy as np
from scipy import sparse
import time
import os

models_folder = os.path.join("..","models")

def read_PLY_file(model_name):
    '''read PLY file'''
    model_path = os.path.join(models_folder,model_name)
    model_file = open(model_path, "r")
    model_data = []
    for line in model_file:
        model_data.append(line.strip().split(" "))
    model_file.close()
    return model_data

def write_PLY_file(model_name, nodes, face_data, resort_indices):
    '''write a PLY file with colors'''

    header = f"ply\nformat ascii 1.0\nelement vertex {len(nodes)}\nproperty float x\nproperty float y\nproperty float z\n" + \
             "property uchar red\nproperty uchar green\nproperty uchar blue\n" + \
             f"element face {len(face_data)}\n" + \
             "property list uchar int vertex_indices\nend_header\n"
    vec_data_str = ""
    resort_indices_normed = resort_indices / resort_indices.shape[0]
    for i,node in enumerate(nodes):
        #color = int(np.round(255 * resort_indices_normed[i]))
        color = int(np.round(255/5 * np.round(resort_indices_normed[i]*5)))
        #color = 255 * (0 if resort_indices_normed[i]<0.5 else 1)
        vec_data_str += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + " " + str(color) + " " + str(color) + " 255\n"
    face_data_str = ""
    for face_datum in face_data:
        face_data_str += "3 " + str(face_datum[0]) + " " + str(face_datum[1]) + " " + str(face_datum[2]) + "\n"

    string_to_write = header + vec_data_str + face_data_str[:-1]

    data_file = open(model_name, "w")
    data_file.write(string_to_write)
    data_file.close()

class Node:
    def __init__(self, location):
        self.location = np.array(location)
        self.connections = []

def extract_model(model_name):
    nodes = []
    model_data = read_PLY_file(model_name)
    face_data = []

    #extract the data
    for line in model_data:
        skip = False
        try:
            float(line[0])
        except:
            skip=True
        if not skip:
            num = float(line[0])
            if int(num) - num == 0.:
                node1 = int(line[1])
                node2 = int(line[2])
                node3 = int(line[3])
                face_data.append([node1, node2, node3])
                if node2 not in nodes[node1].connections:
                    nodes[node1].connections.append(node2)
                if node3 not in nodes[node1].connections:
                    nodes[node1].connections.append(node3)
                if node1 not in nodes[node2].connections:
                    nodes[node2].connections.append(node1)
                if node3 not in nodes[node2].connections:
                    nodes[node2].connections.append(node3)
                if node1 not in nodes[node3].connections:
                    nodes[node3].connections.append(node1)
                if node2 not in nodes[node3].connections:
                    nodes[node3].connections.append(node2)
            else:
                location = [float(line[0]), float(line[1]), float(line[2])]
                nodes.append(Node(location))

    #remove disconnected nodes and decrement all indices after them

    print("removing disconnected nodes and decrementing all indices after the indices of the removed nodes")
    to_delete = []
    decrement_stops = []
    #get unconnected nodes and their indices
    for i in np.arange(len(nodes)):
        node = nodes[i]
        if len(node.connections) == 0:
            to_delete.append(node)
            decrement_stops.append(i)
    #get the altered indices by decrementing indices after the indices of unconnected nodes
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
    #put altered indices in face data
    for i in np.arange(len(face_data)):
        for j in np.arange(len(face_data[i])):
            face_data[i][j] = altered_indices[face_data[i][j]]
    #delete unconnected nodes
    for node in to_delete:
        nodes.remove(node)
    print(f"Done. Removed {len(to_delete)} nodes.")

    return nodes, face_data




def construct_Laplacian(nodes):
    '''construct graph Laplacian'''
    N = len(nodes)
    L = sparse.lil_matrix((N,N), dtype=float)
    for i in np.arange(N):
        node = nodes[i]
        node_sum = 0.
        #adjusted_location = node.location*np.array([0.9,1.,0.95])
        for j in node.connections:
            neighbor = nodes[j]
            #adjusted_neighbor_location = neighbor.location*np.array([0.9,1.,0.95])
            L[i,j] = -1./np.linalg.norm(node.location - neighbor.location)#(adjusted_location - adjusted_neighbor_location)#
            node_sum -= L[i,j]
        L[i,i] = node_sum
    return sparse.csr_matrix(L)

model_name = "bun_zipper.ply"#"dragon_vrip.ply"#
nodes, face_data = extract_model(model_name)       #extract the model from the PLY file
L = construct_Laplacian(nodes)          #construct graph Laplacian
#L_explicit = L.toarray()
#print("L-L.T\n",L_explicit-L_explicit.T)
print()
W, V = sparse.linalg.eigsh(L, k=2, which="SM")  #extract Fiedler vector
print(W, V.T)
#print(np.matmul(L_explicit,V[:,0]))
resort_indices = np.argsort(V[:,1])     #get ordering from Fiedler vector
print(resort_indices)
write_PLY_file(model_name, nodes, face_data, resort_indices)



"""print("Hello world")

solve_start = time.perf_counter_ns()
A = np.array([[3, 2, 0], [1, -1, 0], [0, 5, 1]])
b = np.array([2, 4, -1])
x = linalg.solve(A, b)
solve_end = time.perf_counter_ns()
print(x, "\n", np.dot(A,x)-b)

time_to_run_sims = (solve_end - solve_start) / 1e9
print('Time to run solve:', time_to_run_sims, 's\t\t(', time_to_run_sims/60., 'm)\t\t(', time_to_run_sims/3600., 'h)')"""