#import os
import sys

input_file = sys.argv[1]
f = open(input_file, "r")
almost_write_nodes = False
almost_write_elements = False
write_nodes = False
write_elements = False
count=0

out_nodes = open(input_file[:-4] + ".node", "w")
out_elements = open(input_file[:-4] + ".ele", "w")
for line in f:
	print(count)
	line=line.strip()
	if(len(line) == 0):
		write_nodes = False
		if write_elements:
			write_elements = False
			break
	
	if("*VERTICES" in line):
		almost_write_nodes = True
		continue
	if almost_write_nodes:
		out_nodes.write(line + "\n")
		almost_write_nodes = False
		write_nodes = True
		continue
	if write_nodes:
		parts = line.split(" ")
		rest_of_line = ""
		for part in parts[1:]:
			rest_of_line += " " + part
		out_nodes.write(str(int(parts[0])-1) + rest_of_line + "\n")
	if("TET" in line):
		almost_write_elements = True
		continue
	if almost_write_elements:
		out_elements.write(line + "\n")
		almost_write_elements = False
		write_elements = True
		continue
	if write_elements:
		parts = line.split(" ")
		for i,part in enumerate(parts):
			part = parts[i]
			out_elements.write(str(int(part)-1) + (" " if i<len(parts)-1 else "\n"))
	count+=1
f.close()


out_nodes.close()
out_elements.close()

