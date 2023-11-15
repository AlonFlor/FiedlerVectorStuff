import numpy as np
from scipy import sparse
import time
import os
import subprocess
import pyamg

models_folder = os.path.join("..","models")

def read_data_file(model_name):
    '''read data file'''
    model_path = os.path.join(models_folder,model_name)
    model_file = open(model_path, "r")
    model_data = []
    for line in model_file:
        model_data.append(line.strip().split(" "))
    model_file.close()
    return model_data

def write_PLY_file(model_name, nodes, face_data, vect, contour=False):
    '''write a PLY file with colors'''

    header = f"ply\nformat ascii 1.0\nelement vertex {len(nodes)}\nproperty float x\nproperty float y\nproperty float z\n" + \
             "property uchar red\nproperty uchar green\nproperty uchar blue\n" + \
             f"element face {len(face_data)}\n" + \
             "property list uchar int vertex_indices\nend_header\n"
    vec_data_str = ""
    for i,node in enumerate(nodes):
        color = int(np.round(255 * vect[i]))
        if contour:
            color = 255 if (int(np.round(255 * vect[i])) % 10 > 5) else 0
        #color = int(np.round(255 * resort_indices_normed[i]))
        #color = int(np.round(255/5 * np.round(resort_indices_normed[i]*5)))
        #color = 255 * (0 if resort_indices_normed[i]<0.5 else 1)
        vec_data_str += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + " " + str(color) + " " + str(color) + " 255\n"
    face_data_str = ""
    for face_datum in face_data:
        face_data_str += "3 " + str(face_datum[0]) + " " + str(face_datum[1]) + " " + str(face_datum[2]) + "\n"

    string_to_write = header + vec_data_str + face_data_str[:-1]

    data_file = open(model_name, "w")
    data_file.write(string_to_write)
    data_file.close()

def write_basic_PLY_file(model_name, nodes, face_data):
    '''write a PLY file with colors'''

    header = f"ply\nformat ascii 1.0\nelement vertex {len(nodes)}\nproperty float x\nproperty float y\nproperty float z\n" + \
             f"element face {len(face_data)}\n" + \
             "property list uchar int vertex_indices\nend_header\n"
    vec_data_str = ""
    for i,node in enumerate(nodes):
        #color = int(np.round(255 * resort_indices_normed[i]))
        #color = int(np.round(255/5 * np.round(resort_indices_normed[i]*5)))
        #color = 255 * (0 if resort_indices_normed[i]<0.5 else 1)
        vec_data_str += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"
    face_data_str = ""
    for face_datum in face_data:
        face_data_str += "3 " + str(face_datum[0]) + " " + str(face_datum[1]) + " " + str(face_datum[2]) + "\n"

    string_to_write = header + vec_data_str + face_data_str[:-1]

    data_file = open(model_name, "w")
    data_file.write(string_to_write)
    data_file.close()



def read_numerical_csv_file(file_path, num_type=float):
    file = open(file_path)
    data_raw = []
    for line in file:
        data_raw.append(line.strip().split(","))
    data = []
    for line in data_raw[1:]:
        line_data = []
        for item in line:
            line_data.append(num_type(item))
        data.append(line_data)
    file.close()
    return np.array(data)


'''Note: bfs rankings got an almost-good result with the lobpcg for the medium dragon.'''
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



class Node:
    def __init__(self, location):
        self.location = np.array(location)
        self.connections = []

def extract_model(model_name):
    nodes = []
    model_data = read_data_file(model_name)
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

    #remove disconnected regions and decrement all indices after them

    print("removing disconnected areas and decrementing all indices after the indices of the removed nodes")
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
    #get the altered indices by decrementing indices after the indices of nodes in unconnected areas
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
            L[i,j] = -1.#/np.linalg.norm(node.location - neighbor.location)#(adjusted_location - adjusted_neighbor_location)#
            node_sum -= L[i,j]
        L[i,i] = node_sum
    return sparse.csr_matrix(L)


def get_Laplacian_libigl(nodes, face_data):
    #set up the nodes and faces for libigl
    vectors_string = ""
    for node in nodes:
        vectors_string += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"
    data_file = open("vectors.txt", "w")
    data_file.write(vectors_string)
    data_file.close()
    faces_string = ""
    for face_data_list in face_data:
        faces_string += str(face_data_list[0]) + " " + str(face_data_list[1]) + " " + str(face_data_list[2]) + "\n"
    data_file = open("faces.txt", "w")
    data_file.write(faces_string)
    data_file.close()

    #call the libgl code
    command_str = "./generate_cotangent_Laplacian"
    process = subprocess.Popen(command_str, shell=True, stdout=subprocess.PIPE)
    process.wait()

    #read the csv file returned
    data = read_numerical_csv_file("cotangent_Laplacian.csv")

    #generate Laplacian
    N = len(nodes)
    L = sparse.lil_matrix((N, N), dtype=float)
    for i in np.arange(len(data)):
        row,col,value = data[i]
        row=int(row)
        col=int(col)
        L[row, col] = value
    return sparse.csr_matrix(L)


def interpolation_matrix_for_coarsening(L, M):
    '''Return a nxm interpolation matrix that is the coarsening of the graph whose Laplacian is L.'''
    weights_positive = 1.
    if L[0,0] > 0.:
        weights_positive = -1.

    N = L.shape[0]
    threshold = 0.05
    threshold_incr = 0.05

    #get representatives
    R = np.zeros((N))
    R_list = []
    iters = 0
    ordering = np.arange(N)
    np.random.shuffle(ordering)
    while True:
        node_index = ordering[iters]

        r_node_num = 0.
        r_node_denom = 0.
        row = np.array(L.getrow(node_index).tocoo().col)
        for i in row:
            if node_index==i:
                continue
            weight = weights_positive*L[node_index,i]
            r_node_denom += weight
            if R[i] == 1:
                r_node_num += weight
        r_node = r_node_num/r_node_denom
        if r_node < threshold:
            R[node_index] = 1.
            R_list.append(node_index)
            if len(R_list) == M:
                break

        #every n nodes, increment the threshold
        iters += 1
        if iters % N == 0:
            threshold += threshold_incr
            iters = 0

    '''region_indices = np.zeros((N))
    distances_to_center = np.zeros((N))
    visited = np.zeros((N))
    region_queues = []
    for i in np.arange(1):#M):
        region_queue = []
        rep_index = R_list[i]
        row = np.array(L.getrow(rep_index).tocoo().col)
        for j in row:
            region_queue.append(j)
            region_indices[j] = i+1
            distances_to_center[j] = 1.
            visited[j] = 1.
        distances_to_center[rep_index]=0.
        region_queues.append(region_queue)
    while region_queues:
        to_delete = []
        for i in np.arange(len(region_queues)):
            region_queue = region_queues[i]
            next = region_queue.pop(0)
            row = np.array(L.getrow(next).tocoo().col)

            print(len(row))
            if visited[next]!=1.:
                print("ha")
                exit(1)
            for j in row:
                if j!=next:
                    if visited[j] != 0.:
                        if region_indices[j] == 0:
                            print("haha")
                            exit(1)
                        if region_indices[j] != i+1:
                            if distances_to_center[j] > distances_to_center[next] + 1.:
                                region_queue.append(j)
                                region_indices[j] = i+1
                                distances_to_center[j] = distances_to_center[next] + 1.
                    else:
                        region_queue.append(j)
                        region_indices[j] = i+1
                        distances_to_center[j] = distances_to_center[next] + 1.
                    visited[j] = 1.

            print(len(region_queue),end=" ")
            if len(region_queue)==0:
                to_delete.append(region_queue)
        print()
        for region_queue in to_delete:
            region_queues.remove(region_queue)
        print("\t\t",len(region_queues),"\t\t\t",np.sum(visited),"\t",N)
    print("\n\n", len(region_queues), "\t\t\t", np.sum(visited), "\t", N)

    for i in np.arange(N):
        if visited[i]==0:
            row = np.array(L.getrow(i).tocoo().col)
            print("visited",i,"neighbors:",end=" ")
            for j in row:
                if j!=i:
                    print(j,visited[j],end="\t\t")
            print()
    exit()

    #build the interpolation matrix
    print("going to build it",N,M)
    A = sparse.lil_matrix((N,M), dtype=float)
    for i in np.arange(N):
        region_index = int(region_indices[i])-1
        print(region_index)
        A[i,region_index]=1.'''

    #build the interpolation matrix
    print("going to build it",N,M)
    A = sparse.lil_matrix((N,M), dtype=float)
    for i in np.arange(M):
        rep_index = R_list[i]
        #weights from non-representatives to representatives
        row = np.array(L.getrow(rep_index).tocoo().col)
        print(i,len(row))
        for j in row:
            if R[j] != 0.:
                continue
            j_row = np.array(L.getrow(j).tocoo().col)
            denom = 0.
            for k in j_row:
                denom += weights_positive * k * R[k]
            A[j,i] = weights_positive * L[rep_index, j] / denom
        #weights between representatives
        for j in np.arange(M):
            if j == i:
                continue
            A[rep_index,j] = 1.
    print("built it")

    return sparse.csr_matrix(A)


'''Note: bfs rankings got an almost-good result with the lobpcg for the medium dragon even before I discovered that it is a disconnected model.'''
def bfs_rankings(nodes):
    queue = []
    visited = np.zeros((len(nodes)))
    scores = np.zeros((len(nodes)))
    current_score = 0.

    #find center of model
    center_location = np.zeros((3))
    for node in nodes:
        center_location += node.location
    center_location /= len(nodes)

    #get the start node furthest from the center
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
                scores[j] = current_score + 0.
    return scores


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


def multigrid_solve(L, target_size, multiplier):
    if L.shape[0] <= target_size:
        return pyamg_solve(L)

    interpolatation_matrix = interpolation_matrix_for_coarsening(L, int(np.ceil(L.shape[0]*multiplier)))
    print("interpolated")

    #note: this is way slower than it should be, due to MxM step in creating interpolation matrix, and the steps to make coarse L.

    coarse_L = interpolatation_matrix.T @ (L @ interpolatation_matrix.tocsc()).tocsc()
    guess = multigrid_solve(coarse_L, target_size, multiplier)[1][:,1]

    #refinement
    de_interpolated_guess = interpolatation_matrix @ guess
    return pyamg_solve(L, de_interpolated_guess)


model_name = "dragon_vrip_connected.ply"#"dragon_vrip_res2_connected.ply"#"dragon_vrip_res2_connected.ply"#"Armadillo_digital.ply"#"bun_connected.ply"#
nodes, face_data = extract_model(model_name)       #extract the model from the PLY file

#write_basic_PLY_file(model_name, nodes, face_data)
#exit()

#L = get_Laplacian_libigl(nodes, face_data)         #cotangent Laplacian
L = construct_Laplacian(nodes)          #construct combinatorial Laplacian
#print(L)
#L_explicit = L.toarray()
#print("L-L.T\n",L_explicit-L_explicit.T)
print()
print("extract Fiedler vector")
#W, V = sparse.linalg.eigsh(L, k=2, which="SM")  #extract Fiedler vector
W, V = pyamg_solve(L)
#W, V = multigrid_solve(L, 10000, 0.05)

#print(W, "\n", V.T)
#print(np.matmul(L_explicit,V[:,0]))
Fiedler_vector = V[:,1]

#data = read_numerical_csv_file("Fiedler.csv")
#Fiedler_vector = data.reshape((data.shape[0]))

print("Fiedler_vector",Fiedler_vector)
#print("Fiedler_vector",Fiedler_vector,"\ne-value:",W[1])

#get ordering from Fiedler vector
resort_indices = np.argsort(Fiedler_vector)
resort_ranks = np.zeros_like(Fiedler_vector)
resort_ranks[resort_indices] = np.arange(len(Fiedler_vector))
resort_ranks = resort_ranks / resort_ranks.shape[0]

#print out data
write_PLY_file(model_name, nodes, face_data, resort_ranks)
write_PLY_file("contour_"+model_name, nodes, face_data, resort_ranks, contour=True)




"""print("Hello world")

solve_start = time.perf_counter_ns()
A = np.array([[3, 2, 0], [1, -1, 0], [0, 5, 1]])
b = np.array([2, 4, -1])
x = linalg.solve(A, b)
solve_end = time.perf_counter_ns()
print(x, "\n", np.dot(A,x)-b)

time_to_run_sims = (solve_end - solve_start) / 1e9
print('Time to run solve:', time_to_run_sims, 's\t\t(', time_to_run_sims/60., 'm)\t\t(', time_to_run_sims/3600., 'h)')"""