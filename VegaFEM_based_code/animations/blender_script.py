import os
import bpy
import math
import sys

project_directory="/common/home/af656/Desktop/Fiedler_vector/possible_examples/VegaFEM_based_code/"
object_models_directory=project_directory+"output/"

def make_folder_if_not_exists(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

#model_name = "scrambled_dragon"
#constitutive_model_name = "corotationalLinear"

args = sys.argv
model_name = args[5]
constitutive_model_name = args[6]

model_attempt_number = 0
animation_name = model_name + "_" + constitutive_model_name + "_" + str(model_attempt_number)
output_dir = project_directory + "animations/"+animation_name+"/"
make_folder_if_not_exists(output_dir)
output_dir += "images/"
make_folder_if_not_exists(output_dir)

height_offset = 7.
i=0
while(True):
    str_i = str(i).zfill(4)
    print(i)
    try:
        bpy.context.scene.render.filepath = output_dir + str_i
        bpy.ops.import_mesh.ply(filepath = object_models_directory+animation_name+"_"+str_i+".ply")#, axis_up='Y')
    except:
        break
    #cannot use axis_up parameter, so doing that manually
    object_name = bpy.context.scene.objects[list(bpy.context.scene.objects)[-1].name].name
    current_shape = bpy.context.scene.objects[object_name]
    current_shape.rotation_euler.x = math.pi/2
    print(current_shape.rotation_euler)

    '''if i>0:
        current_shape.location=(1000.,1000.,1000.)
        current_shape.keyframe_insert(data_path="location", frame=i-1)
    current_shape.location=(0.,0.,0.+height_offset)
    current_shape.keyframe_insert(data_path="location", frame=i)
    
    current_shape.location=(1000.,1000.,1000.)
    current_shape.keyframe_insert(data_path="location", frame=i+1)'''

    #render    
    bpy.ops.render.render(write_still=True)

    #delete
    bpy.ops.object.select_all(action='DESELECT')
    for o in bpy.context.scene.objects:
        if o.name == object_name:
            o.select_set(True)
        else:
            o.select_set(False)
    bpy.ops.object.delete()
    
    i+=1
    
#bpy.data.scenes[0].frame_start = 0
#bpy.data.scenes[0].frame_end = i-1

