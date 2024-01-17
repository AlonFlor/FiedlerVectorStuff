import os
import bpy
import math
import sys

project_directory="/data/local/af656/Fiedler_vector/possible_examples/ADMM_Projective_Dynamics_extended/"

def make_folder_if_not_exists(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

args = sys.argv
object_models_directory=project_directory+"output/"
model_name = args[5]
scenario_name = args[6]

animation_name = model_name + "_" + scenario_name
base_output_dir = project_directory + "animations/" + animation_name + "/"
make_folder_if_not_exists(base_output_dir)
output_dir = base_output_dir + "images/"
make_folder_if_not_exists(output_dir)


#create the material
mat = bpy.data.materials.new(name="my_color")
mat.use_nodes = True
principled = mat.node_tree.nodes["Principled BSDF"]
principled.inputs["Base Color"].default_value = (0.710, 0.396, 0.114, 1.) #Armadillo
mat.use_fake_user = True

#height_offset = 7.
i=0
while(True):
    str_i = str(i).zfill(4)
    print(i)
    num_objects_start = len(list(bpy.context.scene.objects))
    try:
        bpy.context.scene.render.filepath = output_dir + str_i
        bpy.ops.import_scene.obj(filepath = object_models_directory+animation_name+"_"+str_i+".obj")#, axis_up='Y')
    except:
        break

    '''if i>0:
        current_shape.location=(1000.,1000.,1000.)
        current_shape.keyframe_insert(data_path="location", frame=i-1)
    current_shape.location=(0.,0.,0.+height_offset)
    current_shape.keyframe_insert(data_path="location", frame=i)
    
    current_shape.location=(1000.,1000.,1000.)
    current_shape.keyframe_insert(data_path="location", frame=i+1)'''
    
    mat = bpy.data.materials.get("my_color")
    # Assign matrial to object
    if current_shape.data.materials:
        # assign to 1st material slot
        current_shape.data.materials[0] = mat
    else:
        # no slots
        current_shape.data.materials.append(mat)

    #render    
    bpy.ops.render.render(write_still=True)

    #delete
    while(len(bpy.context.scene.objects)> num_objects_start):
        bpy.ops.object.select_all(action='DESELECT')
        bpy.context.scene.objects[-1].select_set(True)
        bpy.ops.object.delete()
    bpy.ops.outliner.orphans_purge(do_local_ids=True, do_linked_ids=True, do_recursive=True)
    
    i+=1
    
#make video
os.system(f"python3 make_video.py -folder {animation_name} -fps 24")

