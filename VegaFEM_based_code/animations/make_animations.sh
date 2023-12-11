#!/bin/bash
blender -b generate_animations.blend -P blender_script.py scrambled_dragon corotationalLinear
blender -b generate_animations.blend -P blender_script.py scrambled_dragon linear
blender -b generate_animations.blend -P blender_script.py scrambled_dragon StVK
blender -b generate_animations.blend -P blender_script.py scrambled_dragon massSpring
blender -b generate_animations.blend -P blender_script.py dragon corotationalLinear
blender -b generate_animations.blend -P blender_script.py dragon linear
blender -b generate_animations.blend -P blender_script.py dragon StVK
blender -b generate_animations.blend -P blender_script.py dragon massSpring
blender -b generate_animations.blend -P blender_script.py reordered_dragon corotationalLinear
blender -b generate_animations.blend -P blender_script.py reordered_dragon linear
blender -b generate_animations.blend -P blender_script.py reordered_dragon StVK
blender -b generate_animations.blend -P blender_script.py reordered_dragon massSpring

python3 make_video.py -folder scrambled_dragon_corotationalLinear_0 -fps 24
python3 make_video.py -folder scrambled_dragon_linear_0 -fps 24
python3 make_video.py -folder scrambled_dragon_StVK_0 -fps 24
python3 make_video.py -folder scrambled_dragon_massSpring_0 -fps 24
python3 make_video.py -folder dragon_corotationalLinear_0 -fps 24
python3 make_video.py -folder dragon_linear_0 -fps 24
python3 make_video.py -folder dragon_StVK_0 -fps 24
python3 make_video.py -folder dragon_massSpring_0 -fps 24
python3 make_video.py -folder reordered_dragon_corotationalLinear_0 -fps 24
python3 make_video.py -folder reordered_dragon_linear_0 -fps 24
python3 make_video.py -folder reordered_dragon_StVK_0 -fps 24
python3 make_video.py -folder reordered_dragon_massSpring_0 -fps 24
