import numpy as np
import os

class Node:
    def __init__(self, location):
        self.location = np.array(location)
        self.connections = []

def read_node_ele_data_files(model_path):
    model_data_node = []
    print("getting " + model_path+".node")
    print("getting " + model_path+".ele")
    model_file = open(model_path+".node", "r")
    for line in model_file:
        model_data_node.append(line.strip().split(" "))
    model_file.close()
    model_data_ele = []
    model_file = open(model_path+".ele", "r")
    for line in model_file:
        model_data_ele.append(line.strip().split(" "))
    model_file.close()
    return [model_data_node, model_data_ele]

def read_data_file(model_path):
    '''read data file'''
    if("." not in model_path[-4:]):
        #dealing with .node and .ele files
        return read_node_ele_data_files(model_path)
    print("getting " + model_path)
    model_file = open(model_path, "r")
    model_data = []
    for line in model_file:
        model_data.append(line.strip().split(" "))
    model_file.close()
    return model_data


def extract_model_node_ele(model_data):
    nodes = []
    tet_data = []
    model_data_nodes, model_data_tets = model_data
    
    #handle nodes
    for line in model_data_nodes[1:]:
        if line[0].startswith("#"):
            continue
        numbers = []
        for component in line:
            try:
                float(component)
            except:
                continue
            numbers.append(float(component))
        location = [float(numbers[1]),float(numbers[2]),float(numbers[3])]
        nodes.append(Node(location))
        
    #handle tets
    for line in model_data_tets[1:]:
        if line[0].startswith("#"):
            continue
        numbers = []
        for component in line:
            try:
                int(component)
            except:
                continue
            numbers.append(int(component))
        node1 = int(numbers[1])
        node2 = int(numbers[2])
        node3 = int(numbers[3])
        node4 = int(numbers[4])
        tet = [node1,node2,node3,node4]
        tet_data.append(tet)
        if node2 not in nodes[node1].connections:
            nodes[node1].connections.append(node2)
        if node3 not in nodes[node1].connections:
            nodes[node1].connections.append(node3)
        if node4 not in nodes[node1].connections:
            nodes[node1].connections.append(node4)

        if node1 not in nodes[node2].connections:
            nodes[node2].connections.append(node1)
        if node3 not in nodes[node2].connections:
            nodes[node2].connections.append(node3)
        if node4 not in nodes[node2].connections:
            nodes[node2].connections.append(node4)

        if node1 not in nodes[node3].connections:
            nodes[node3].connections.append(node1)
        if node2 not in nodes[node3].connections:
            nodes[node3].connections.append(node2)
        if node4 not in nodes[node3].connections:
            nodes[node3].connections.append(node4)

        if node1 not in nodes[node4].connections:
            nodes[node4].connections.append(node1)
        if node2 not in nodes[node4].connections:
            nodes[node4].connections.append(node2)
        if node3 not in nodes[node4].connections:
            nodes[node4].connections.append(node3)
    
    return nodes, tet_data, None
    

def extract_model_PLY(model_data):
    nodes = []
    face_data = []

    #extract the data
    for line in model_data:
        skip = False
        try:
            float(line[0])
        except:
            skip=True
        if not skip:
            isFace = True
            for element in line:
                if "." in element or "-" in element:
                    isFace = False
                    break
            if isFace:
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
    return nodes, face_data, None

def extract_model_OBJ(model_data):
    nodes = []
    face_data = []
    other = []

    #extract the data
    for line in model_data:
        if line[0] == "v":
            location = [float(line[1]), float(line[2]), float(line[3])]
            nodes.append(Node(location))
        elif line[0] == "f":
            new_line = []
            #get rid of "/"
            for i in np.arange(1, len(line)):
                line_part = line[i]
                if "/" in line_part:
                    line_part = line_part[:line_part.index("/")]
                new_line.append(line_part)
            node1 = int(new_line[0]) - 1
            node2 = int(new_line[1]) - 1
            node3 = int(new_line[2]) - 1

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
            other.append(line)

    return nodes, face_data, other


def extract_model_VEG(model_data):
    nodes = []
    tet_data = []
    other = []

    mode = ""
    skip = False

    #extract the data
    count= 0
    for line in model_data:
        count+=1
        if skip:
            skip=False
            continue
        if mode != "other" and len(line[0])==0:
            continue
        if line[0]=="#":
            continue
        if line[0].startswith("*"):
            if line[0]=="*VERTICES":
                mode = "nodes"
                skip=True
            elif line[0]=="*ELEMENTS":
                mode = "faces"
            else:
                mode = "other"
                other.append(line)
            continue
        if mode == "other":
            other.append(line)
            continue

        try:
            float(line[0])
        except:
            if mode=="nodes":
                print(f"In veg file, line {count-1} should be a number but is not.")
                exit(1)
            elif mode == "faces":
                skip=True #we are in the line for TET or CUBIC, skip the line after it
                continue

        #dealing with numbers
        if mode=="nodes":
            location = [float(line[1]), float(line[2]), float(line[3])]
            nodes.append(Node(location))
        elif mode=="faces":
            node1 = int(line[1]) - 1
            node2 = int(line[2]) - 1
            node3 = int(line[3]) - 1
            node4 = int(line[4]) - 1
            tet_data.append([node1, node2, node3, node4])
            #tet_data.append([node1, node2, node3])
            #tet_data.append([node1, node2, node4])
            #tet_data.append([node2, node3, node4])
            #tet_data.append([node1, node3, node4])
            if node2 not in nodes[node1].connections:
                nodes[node1].connections.append(node2)
            if node3 not in nodes[node1].connections:
                nodes[node1].connections.append(node3)
            if node4 not in nodes[node1].connections:
                nodes[node1].connections.append(node4)

            if node1 not in nodes[node2].connections:
                nodes[node2].connections.append(node1)
            if node3 not in nodes[node2].connections:
                nodes[node2].connections.append(node3)
            if node4 not in nodes[node2].connections:
                nodes[node2].connections.append(node4)

            if node1 not in nodes[node3].connections:
                nodes[node3].connections.append(node1)
            if node2 not in nodes[node3].connections:
                nodes[node3].connections.append(node2)
            if node4 not in nodes[node3].connections:
                nodes[node3].connections.append(node4)

            if node1 not in nodes[node4].connections:
                nodes[node4].connections.append(node1)
            if node2 not in nodes[node4].connections:
                nodes[node4].connections.append(node2)
            if node3 not in nodes[node4].connections:
                nodes[node4].connections.append(node3)

    return nodes, tet_data, other


def extract_model(model_name, model_folder=None):
    model_path = model_name
    if model_folder is not None:
        model_path = os.path.join(model_folder, model_name)
    model_data = read_data_file(model_path)
    model_name_to_check = model_name.lower()
    if model_name_to_check.endswith(".ply"):
        return extract_model_PLY(model_data)
    elif model_name_to_check.endswith(".obj"):
        return extract_model_OBJ(model_data)
    elif model_name_to_check.endswith(".veg"):
        return extract_model_VEG(model_data)
    else:
        return extract_model_node_ele(model_data)

def write_color_PLY_file(model_name, nodes, face_data, rankings, contour=False):
    '''write a PLY file with colors'''

    header = f"ply\nformat ascii 1.0\nelement vertex {len(nodes)}\nproperty float x\nproperty float y\nproperty float z\n" + \
             "property uchar red\nproperty uchar green\nproperty uchar blue\n" + \
             f"element face {len(face_data)}\n" + \
             "property list uchar int vertex_indices\nend_header\n"
    vec_data_str = ""
    for i, node in enumerate(nodes):
        color = int(np.round(255 * rankings[i]))
        if contour:
            color = 255 if (int(np.round(255 * rankings[i])) % 10 > 5) else 0
        # color = int(np.round(255 * resort_indices_normed[i]))
        # color = int(np.round(255/5 * np.round(resort_indices_normed[i]*5)))
        # color = 255 * (0 if resort_indices_normed[i]<0.5 else 1)
        vec_data_str += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + " " + str(color) + " " + str(color) + " 255\n"

    face_data_str = ""
    if len(face_data[0]) == 4:
        #handle tetrahedral meshes
        for face in face_data:
            face_data_str += "4 " + str(face[0]) + " " + str(face[1]) + " " + str(face[2]) + " " + str(face[3]) + "\n"
    else:
        for face in face_data:
            face_data_str += "3 " + str(face[0]) + " " + str(face[1]) + " " + str(face[2]) + "\n"

    string_to_write = header + vec_data_str + face_data_str[:-1]

    data_file = open(model_name, "w")
    data_file.write(string_to_write)
    data_file.close()

def reorder_faces(face_data, face_resort_ranks, num_nodes):
    # make sure faces are sorted in ascending order, by vertex indices in the faces
    first_face_index = np.argmin(np.array(face_resort_ranks))
    if face_data[first_face_index][0] > 0.5*num_nodes:
        for i in np.arange(len(face_resort_ranks)):
            face_resort_ranks[i] = len(face_resort_ranks) - 1 - face_resort_ranks[i]

    #resort faces based on ranks
    faces_ranks_dict = {}
    for i in np.arange(len(face_data)):
        faces_ranks_dict[face_resort_ranks[i]] = face_data[i]
    new_face_data = []
    for i in np.arange(len(face_data)):
        face = faces_ranks_dict[i]
        new_face_data.append(face)

    return new_face_data

def write_reordered_node_ele_files(model_name, nodes, tet_data, rankings, face_rankings=None):
    '''write reordered .node and .ele files'''
    
    vec_data_str = f"{len(nodes)} 3 0 0\n"
    # write nodes
    nodes_ranks_dict = {}
    for i in np.arange(len(nodes)):
        nodes_ranks_dict[rankings[i]]=nodes[i]
    for i in np.arange(len(nodes)):
        node = nodes_ranks_dict[i]
        vec_data_str += str(i) + " " + str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"
    data_file = open(model_name+".node", "w")
    data_file.write(vec_data_str)
    data_file.close()

    # change vertices in faces
    for i in np.arange(len(tet_data)):
        for j in np.arange(len(tet_data[i])):
            tet_data[i][j] = rankings[tet_data[i][j]]

    # reorder faces
    if face_rankings is not None:
        tet_data = reorder_faces(tet_data, face_rankings, len(nodes))

    # write faces
    tet_data_str = f"{len(tet_data)}"
    #handle tetrahedral meshes
    for i in np.arange(len(tet_data)):
        tet = tet_data[i]
        tet_data_str += str(i) + " " + str(tet[0]) + " " + str(tet[1]) + " " + str(tet[2]) + " " + str(tet[3]) + "\n"
    
    data_file = open(model_name+".ele", "w")
    data_file.write(tet_data_str)
    data_file.close()

def write_reordered_PLY_file(model_name, nodes, face_data, rankings, face_rankings=None):
    '''write a reordered PLY file'''

    header = f"ply\nformat ascii 1.0\nelement vertex {len(nodes)}\nproperty float x\nproperty float y\nproperty float z\n" + \
             f"element face {len(face_data)}\n" + \
             "property list uchar int vertex_indices\nend_header\n"
    vec_data_str = ""

    # write nodes
    nodes_ranks_dict = {}
    for i in np.arange(len(nodes)):
        nodes_ranks_dict[rankings[i]]=nodes[i]
    for i in np.arange(len(nodes)):
        node = nodes_ranks_dict[i]
        vec_data_str += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"

    # change vertices in faces
    for i in np.arange(len(face_data)):
        for j in np.arange(len(face_data[i])):
            face_data[i][j] = rankings[face_data[i][j]]

    # reorder faces
    if face_rankings is not None:
        face_data = reorder_faces(face_data, face_rankings, len(nodes))

    # write faces
    face_data_str = ""
    if len(face_data[0]) == 4:
        #handle tetrahedral meshes
        for face in face_data:
            face_data_str += "4 " + str(face[0]) + " " + str(face[1]) + " " + str(face[2]) + " " + str(face[3]) + "\n"
    else:
        for face in face_data:
            face_data_str += "3 " + str(face[0]) + " " + str(face[1]) + " " + str(face[2]) + "\n"

    string_to_write = header + vec_data_str + face_data_str[:-1]

    data_file = open(model_name, "w")
    data_file.write(string_to_write)
    data_file.close()


def write_reordered_VEG_file(model_name, nodes, tet_data, rankings, other, face_rankings=None):
    '''write a reordered VEG file'''

    vec_data_str = "\n*VERTICES\n"
    vec_data_str += f"{len(nodes)} 3 0 0\n"

    # write nodes
    nodes_ranks_dict = {}
    for i in np.arange(len(nodes)):
        nodes_ranks_dict[rankings[i]]=nodes[i]
    for i in np.arange(len(nodes)):
        node = nodes_ranks_dict[i]
        vec_data_str += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"

    # change vertices in tetss
    for i in np.arange(len(tet_data)):
        for j in np.arange(4):
            tet_data[i][j] = rankings[tet_data[i][j]]

    # reorder faces
    if face_rankings is not None:
        tet_data = reorder_faces(tet_data, face_rankings, len(nodes))

    # write faces
    tet_data_str = "\n*ELEMENTS\nTET\n"
    tet_data_str += f"{len(tet_data)} 4 0\n"
    for i in np.arange(len(tet_data)):
        tet = tet_data[i]
        tet_data_str += str(i+1) + " " + str(tet[0]+1) + " " + str(tet[1]+1) + " " + str(tet[2]+1) + " " + str(tet[3]+1) + "\n"

    other_str = "\n"
    for line in other:
        for element in line:
            if len(element) > 0:
                other_str += element + " "
        other_str += "\n"

    string_to_write = vec_data_str + tet_data_str + other_str

    data_file = open(model_name, "w")
    data_file.write(string_to_write)
    data_file.close()


def write_reordered_OBJ_file(model_name, nodes, face_data, rankings, other, face_rankings=None):
    '''write a reordered VEG file'''

    vec_data_str = ""

    # write nodes
    nodes_ranks_dict = {}
    for i in np.arange(len(nodes)):
        nodes_ranks_dict[rankings[i]]=nodes[i]
    for i in np.arange(len(nodes)):
        node = nodes_ranks_dict[i]
        vec_data_str += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"

    # change vertices in faces
    for i in np.arange(len(face_data)):
        for j in np.arange(len(face_data[i])):
            face_data[i][j] = rankings[face_data[i][j]]

    # reorder faces
    if face_rankings is not None:
        face_data = reorder_faces(face_data, face_rankings, len(nodes))

    # write faces
    face_data_str = ""
    for i in np.arange(len(face_data)):
        local_face_data_str = "f"
        face = face_data[i]
        for j in np.arange(len(face)):
            local_face_data_str += " " + str(face[j]+1)+"//"+str(face[j]+1)
        local_face_data_str += "\n"
        face_data_str += local_face_data_str

    other_str = ""
    if other is not None:
        for line in other:
            for i in np.arange(len(line)):
                element = line[i]
                if len(element) > 0:
                    other_str += element
                    if i<len(line)-1:
                        other_str += " "
            other_str += "\n"

    string_to_write = vec_data_str + other_str + face_data_str

    data_file = open(model_name, "w")
    data_file.write(string_to_write)
    data_file.close()

def write_basic_PLY_file(model_name, nodes, face_data):
    '''write a PLY file with colors'''

    header = f"ply\nformat ascii 1.0\nelement vertex {len(nodes)}\nproperty float x\nproperty float y\nproperty float z\n" + \
             f"element face {len(face_data)}\n" + \
             "property list uchar int vertex_indices\nend_header\n"
    vec_data_str = ""
    for i, node in enumerate(nodes):
        # color = int(np.round(255 * resort_indices_normed[i]))
        # color = int(np.round(255/5 * np.round(resort_indices_normed[i]*5)))
        # color = 255 * (0 if resort_indices_normed[i]<0.5 else 1)
        vec_data_str += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"
    face_data_str = ""
    for face in face_data:
        face_data_str += "3 " + str(face[0]) + " " + str(face[1]) + " " + str(face[2]) + "\n"

    string_to_write = header + vec_data_str + face_data_str[:-1]

    data_file = open(model_name, "w")
    data_file.write(string_to_write)
    data_file.close()


