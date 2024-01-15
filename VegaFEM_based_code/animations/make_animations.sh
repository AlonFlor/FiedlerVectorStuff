#!/bin/bash
blender -b generate_animations_improved.blend -P blender_script.py dragon corotationalLinear falling
#blender -b generate_animations_improved.blend -P blender_script.py dragon linear falling
#blender -b generate_animations_improved.blend -P blender_script.py dragon StVK falling
#blender -b generate_animations_improved.blend -P blender_script.py dragon massSpring falling
#
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py dragon corotationalLinear side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py dragon linear side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py dragon StVK side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py dragon massSpring side_motion_constrained
#
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py dragon corotationalLinear falling_constrained
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py dragon linear falling_constrained
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py dragon StVK falling_constrained
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py dragon massSpring falling_constrained
#
