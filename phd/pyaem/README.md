# PyAEM

-- Author: Ivan S. P. Marin
-- email: ispmarin@gmail.com


AEM Crack simulator using the \textit{chi} formalism and direct integration for the coefficients, and comparison with line doublets with ellipsoid form. This simulation does both the crack and the line doublet simulations and compare both for the same problem. Note that this simulation is not solving a system with cracks and line doublets, but solving the system for cracks and them the same system for line doublets, and compares both.

PyAEM is using the packages matplotlib for plotting mad, pyaudio for the END PROGRAM warning

Input files: input_crack.dat: Defines the parameters for the crack simulation, aquifer and uniform flow. It is in the following way:

4 # n: max number of coefficients for the element. This uses the C notation, so coefficients are from 0 to n-1 18 # max_iter: maximum number of iterations for the solver 65 #deprecated 1 # k_ext, aquifer hydraulic conductivity 1 # base, aquifer base, not in use 1000 # height, aquifer max height, not in use 20 # head at the referece point 1000 # x coordinate of the reference point 1000 # y coordinate of the reference point 0.1 # uniform flow Q0 30 # angle of uniform flow, given in degrees

input_ld.dat: Defines the parameters for the crack simulation, aquifer and uniform flow. It is in the following way:

4 # n: max number of coefficients for the element. This uses the C notation, so coefficients are from 0 to n-1 100 # j: max number of coefficients for the far field expansion. This uses the C notation, so coefficients are from 0 to n-1 18 # max_iter: maximum number of iterations for the solver 65 #deprecated 1 # k_ext, aquifer hydraulic conductivity 1 # base, aquifer base, not in use 1000 # height, aquifer max height, not in use 20 # head at the referece point 1000 # x coordinate of the reference point 1000 # y coordinate of the reference point 0.1 # uniform flow Q0 30 # angle of uniform flow, given in degrees

crack.dat : Defines the crack parameters. All the cracks should be inserted in this file.

1000 1 20 20 40 20 -0 k_int ap x1 y1 x2 y2 dep

k_int: crack's hydraulic conductivity ap: crack's aperture x1, y1, x2, y2: coordinates of the beginning and the end of the crack dep: deprecated

The simulation now assumes unconfined flow, with uniform flow. The integration is being done using scipy.integrate quad for the crack and basic_int for the line doublets, and complex numbers for all the inputs.

To run the program, execute python main.py
