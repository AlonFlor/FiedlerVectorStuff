Code for Spectral reordering for faster elasticity simulations (https://dl.acm.org/doi/10.1007/s00371-024-03513-0)

To reorder shapes:
        call main.py on the command line:
                python3 main.py -name <model file name> -folder <path of folder of model> -command <command>
                
        Can optionally include the -f tag to reorder the faces.
        
        Available commands:
                Fiedler: reorder the vertices of the model by its Fiedler vector.
                same: keep the order of the vertices of the model the same.
                scramble: randomly rearrange the vertices of the model.
                det_scramble: scramble the vertices in a deterministic manner.
        Can give a -chunk_size <size of chunks> flag, which only applies when command is set to scramble or det_scramble.
        It controls the chunk size of the scrambled shape. Recommended for det_scramble to work, since the default value is 1.

Code to run the simulated examples is in the VegaFEM_based_code, ADMM_Projective_Dynamics, and ADMM_Projective_Dynamics_extended folders.
Run any of the bash script files in a Linux environment.
