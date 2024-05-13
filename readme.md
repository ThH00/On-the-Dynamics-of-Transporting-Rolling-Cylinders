The repository contains the companion codes to the paper `On the Dynamics of Transporting Rolling Cylinders` that is in preperation for submition to Nonlinear Dynamics in May 2024.

# Frobenius Integrability Criteria
The folder `frobenius-integrability-criteria` constaints a Matlab implementations of the Frobenius Criteria (symbolically and with a numerical example) for
- one cylinder (`constraint_examination_one_cylinder.m`)
- two cylinders subject to 9 constraints: rolling between each cylinder and the ground and rolling at the contact point between the two cylinder (`constraint_examination_two_cylinders_7_constraints.m`)
- two cylinders subject to 7 constraints: rolling between each cylinder and the ground and touching at the contact point between the two cylinders (`constraint_examination_two_cylinders_9_constraints.m`)

# Simulation
To run the nonsmooth generalized-alpha with frictional contact simulation:
1. Create an input file that includes in order the cylinder weight, radius, height, angle of inclination with the vertical, and the initial coordinates of the centers of mass of both cylinders in the manner of `generalized-alpha-simulation/inputs-movie.txt`
2. Indicate your input file in line 12 of `runner.py`
3. Run `generalized-alpha-simulation/runner.py` by typing `python3 runner.py` in the terminal.
Two input files are provided, `generalized-alpha-simulation/inputs_movie.txt` and `generalized-alpha-simulation/inputs-angles.txt`. The latter provides the results presented in Figure 14.

Given the inclination angle of the cylinders, the initial positions of the center of mass of each cylinder is calculated by calling the Matlab function `get-initial-configuration>get_config.m`. To obtain the initial conditions saved in `generalized-alpha-simulation/inputs-angles.txt`, run `generalized-alpha-simulation/obtaining_initial_configs.m`
