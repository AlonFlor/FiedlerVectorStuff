#!/bin/bash
./VegaFEM_program corotationalLinear falling 240 3 scrambled_dragon.veg 3 dragon.veg 3 Fiedler_reordered_dragon.veg 3 vf_scrambled_dragon.veg 3 vf_Fiedler_reordered_dragon.veg
./VegaFEM_program StVK falling_constrained 240 3 scrambled_dragon.veg 3 dragon.veg 3 Fiedler_reordered_dragon.veg 3 vf_scrambled_dragon.veg 3 vf_Fiedler_reordered_dragon.veg
