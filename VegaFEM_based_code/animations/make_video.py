import sys
import subprocess

def make_video(frame_rate, scenario_folder, video_base_name="video", image_type_name=""):
    command_str = f"ffmpeg -framerate {frame_rate} -i {scenario_folder}/images/%04d{image_type_name}.png" \
                  f" -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {scenario_folder}/{video_base_name}{image_type_name}.mp4"
    process = subprocess.Popen(command_str, shell=True, stdout=subprocess.PIPE)
    process.wait()
    
args = sys.argv[1:]
mode = None
folder = None
fps = None
for arg in args:
    if arg=="-help" or arg=="--help":
        print("determine frame rate: -fps [int]")
        print("determine folder to make video: -folder [string]")
    if arg=="-fps":
        mode = "fps"
    elif arg=="-folder":
        mode = "folder"
    elif mode=="fps":
        fps = int(arg)
    elif mode=="folder":
        folder = arg
if folder is not None and fps is not None:
    make_video(fps, folder, video_base_name=folder)
else:
    print("Need to specify frame rate and folder.")
    print("determine frame rate: -fps [int]")
    print("determine folder to make video: -folder [string]")

