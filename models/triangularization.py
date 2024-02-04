import os
import sys
#import numpy as np


name = sys.argv[1]
f = open(name+".ele", "r")
start=False
tets = []
for line in f:
	if line.startswith("#"):
		continue
	line=line.strip()
	if not start:
		#skip over header
		start=True
	else:
		vert_indices_str = line.split()[1:]
		vert_indices = []
		for vert_index_str in vert_indices_str:
			vert_indices.append(int(vert_index_str))
		tets.append(sorted(vert_indices))
f.close()

#get the faces on the tetrahedra and get the biggest vertex index
all_faces = {}
max_vert_index = 0
for tet in tets:
	index_1, index_2, index_3, index_4 = tet
	face_candidate1 = (index_1, index_2, index_3)
	face_candidate2 = (index_1, index_2, index_4)
	face_candidate3 = (index_1, index_3, index_4)
	face_candidate4 = (index_2, index_3, index_4)
	face_candidates = [face_candidate1, face_candidate2, face_candidate3, face_candidate4]
	for face_candidate in face_candidates:
		if face_candidate in all_faces:
			all_faces[face_candidate] += 1
		else:
			all_faces[face_candidate] = 1
			
	max_vert_index = max(index_1, max_vert_index)
	max_vert_index = max(index_2, max_vert_index)
	max_vert_index = max(index_3, max_vert_index)
	max_vert_index = max(index_4, max_vert_index)

#get unique faces: these are on the surface
faces = []
for face in all_faces:
	if all_faces[face] == 1:
		faces.append(face)
print(len(tets), "tetrahedra\t\t", len(all_faces), "faces total","\t\t",len(faces), "surface faces")

#get excluded vertices
included_vert_indices = set()
excluded_vert_indices = set()
for face in faces:
	index_1, index_2, index_3 = face
	included_vert_indices.add(index_1)
	included_vert_indices.add(index_2)
	included_vert_indices.add(index_3)
for i in range(max_vert_index+1):
	if i not in included_vert_indices:
		excluded_vert_indices.add(i)
print(len(excluded_vert_indices))

out = open(os.path.join("triangle_lists",name+"_faces.txt"), "w")
out.write(str(len(faces)) + "\n")
for face in faces:
	out.write(str(face[0]) + " " + str(face[1]) + " " + str(face[2]) + "\n")
out.close()

