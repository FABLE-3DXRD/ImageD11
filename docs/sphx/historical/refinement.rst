Refinement
==========
Assuming you have managed to find some trial orientation matrices the quality of the derived parameters will be quite poor without post refinement of the experimental setup. For each grain the translation of the grain from the rotation axis perturbs the position of the diffraction spots. Also the precise orientation of the rotation axis affects the final peak positions. How to do this?

Input data
==========

| Filtered diffraction peak positions (output from peaksearch menu) for one or more scans. 
| UBI matrices for grains to be treated. 
| Starting geometry parameters (as used in transformation menu).

Output results
==============

| Refined UBI matrix and displacement (x,y,z) for at least one grain in each rotation. 
| "Observed" and calculated values of indexed reflections

FIXME: This should write about fitgrain/filtergrain/refine_em/avg_par
