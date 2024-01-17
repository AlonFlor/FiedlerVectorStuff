#!/bin/bash
blender -b generate_animations_improved.blend -P blender_script.py scrambled_dragon corotationalLinear falling
#blender -b generate_animations_improved.blend -P blender_script.py dragon corotationalLinear falling
#blender -b generate_animations_improved.blend -P blender_script.py Fiedler_reordered_dragon corotationalLinear falling
#blender -b generate_animations_improved.blend -P blender_script.py vf_Fiedler_reordered_dragon corotationalLinear falling
#blender -b generate_animations_improved.blend -P blender_script.py scrambled2_dragon corotationalLinear falling
#blender -b generate_animations_improved.blend -P blender_script.py scrambled3_dragon corotationalLinear falling
#
#blender -b generate_animations_improved.blend -P blender_script.py dragon linear falling
#blender -b generate_animations_improved.blend -P blender_script.py dragon StVK falling
#blender -b generate_animations_improved.blend -P blender_script.py dragon massSpring falling
#
#
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py dragon corotationalLinear side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py dragon linear side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py dragon StVK side_motion_constrained
#
blender -b generate_animations_improved_side_motion.blend -P blender_script.py scrambled_dragon massSpring side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py dragon massSpring side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py Fiedler_reordered_dragon massSpring side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py vf_Fiedler_reordered_dragon massSpring side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py scrambled2_dragon massSpring side_motion_constrained
#blender -b generate_animations_improved_side_motion.blend -P blender_script.py scrambled3_dragon massSpring side_motion_constrained
#
#
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py dragon corotationalLinear falling_constrained
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py dragon linear falling_constrained
#
blender -b generate_animations_improved_down_motion.blend -P blender_script.py scrambled_dragon StVK falling_constrained
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py dragon StVK falling_constrained
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py Fiedler_reordered_dragon StVK falling_constrained
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py vf_Fiedler_reordered_dragon StVK falling_constrained
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py scrambled2_dragon StVK falling_constrained
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py scrambled3_dragon StVK falling_constrained
#
#blender -b generate_animations_improved_down_motion.blend -P blender_script.py dragon massSpring falling_constrained
#
#
