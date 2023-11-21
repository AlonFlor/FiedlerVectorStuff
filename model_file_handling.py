import numpy as np
import os

class Node:
    def __init__(self, location):
        self.location = np.array(location)
        self.connections = []

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
    return nodes, face_data, None

def extract_model_OBJ(model_data):
    nodes = []
    face_data = []
    other = []

    #extract the data
    for line in model_data:
        print(line)
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
    face_data = []
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
            face_data.append([node1, node2, node3, node4])
            #face_data.append([node1, node2, node3])
            #face_data.append([node1, node2, node4])
            #face_data.append([node2, node3, node4])
            #face_data.append([node1, node3, node4])
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

    return nodes, face_data, other


def extract_model(model_name):
    model_data = read_data_file(model_name)
    model_name_to_check = model_name.lower()
    if model_name_to_check.endswith(".ply"):
        return extract_model_PLY(model_data)
    elif model_name_to_check.endswith(".obj"):
        return extract_model_OBJ(model_data)
    elif model_name_to_check.endswith(".veg"):
        return extract_model_VEG(model_data)

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

def bucket_sort(data):
    buckets = {}
    max_key = 0
    for item_list in data:
        bucket_id = min(item_list)
        if bucket_id > max_key:
            max_key = bucket_id
        if bucket_id not in buckets:
            buckets[bucket_id] = [item_list]
        else:
            buckets[bucket_id].append(item_list)
    # sort buckets
    sorted_data = []
    for key in np.arange(max_key + 1):
        # print(key)
        try:
            for item_list in buckets[key]:
                sorted_data.append(item_list)
        except:
            pass
    # print(len(sorted_data))
    return sorted_data

def write_reordered_PLY_file(model_name, nodes, face_data, rankings, rearrange_faces=False):
    '''write a reordered PLY file'''

    header = f"ply\nformat ascii 1.0\nelement vertex {len(nodes)}\nproperty float x\nproperty float y\nproperty float z\n" + \
             f"element face {len(face_data)}\n" + \
             "property list uchar int vertex_indices\nend_header\n"
    vec_data_str = ""

    # write nodes
    for i in np.arange(len(nodes)):
        node_index = np.where(rankings == i)[0][0]
        node = nodes[node_index]
        vec_data_str += str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"

    # change vertices in faces
    for i in np.arange(len(face_data)):
        for j in np.arange(len(face_data[i])):
            face_data[i][j] = rankings[face_data[i][j]]

    # reorder faces
    if rearrange_faces:
        face_data = bucket_sort(face_data)

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


def write_reordered_VEG_file(model_name, nodes, face_data, rankings, other, rearrange_faces=False):
    '''write a reordered VEG file'''

    vec_data_str = "\n*VERTICES\n"
    vec_data_str += f"{len(nodes)} 3 0 0\n"

    # write nodes
    for i in np.arange(len(nodes)):
        node_index = np.where(rankings == i)[0][0]
        node = nodes[node_index]
        vec_data_str += str(i+1) + " " + str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"

    # change vertices in faces
    for i in np.arange(len(face_data)):
        for j in np.arange(4):
            face_data[i][j] = rankings[face_data[i][j]]

    # reorder faces
    if rearrange_faces:
        face_data = bucket_sort(face_data)

    # write faces
    face_data_str = "\n*ELEMENTS\nTET\n"
    face_data_str += f"{len(face_data)} 4 0\n"
    for i in np.arange(len(face_data)):
        face = face_data[i]
        face_data_str += str(i+1) + " " + str(face[0]+1) + " " + str(face[1]+1) + " " + str(face[2]+1) + " " + str(face[3]+1) + "\n"

    other_str = "\n"
    for line in other:
        for element in line:
            if len(element) > 0:
                other_str += element + " "
        other_str += "\n"

    string_to_write = vec_data_str + face_data_str + other_str

    data_file = open(model_name, "w")
    data_file.write(string_to_write)
    data_file.close()


def write_reordered_OBJ_file(model_name, nodes, face_data, rankings, other, rearrange_faces=False):
    '''write a reordered VEG file'''

    vec_data_str = ""

    # write nodes
    for i in np.arange(len(nodes)):
        node_index = np.where(rankings == i)[0][0]
        node = nodes[node_index]
        vec_data_str += "v " + str(node.location[0]) + " " + str(node.location[1]) + " " + str(node.location[2]) + "\n"

    # change vertices in faces
    for i in np.arange(len(face_data)):
        for j in np.arange(len(face_data[i])):
            face_data[i][j] = rankings[face_data[i][j]]

    # reorder faces
    if rearrange_faces:
        face_data = bucket_sort(face_data)

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


