# On the behavior of Lagrange multipliers in convex and non-convex infeasible interior point methods

This git repo contains code to allow one to replicate the results of the paper:

"On the behavior of Lagrange multipliers in convex and non-convex infeasible interior point methods. Gabriel Haeser, Oliver Hinder, Yinyu Ye. 2018."
The arxiv paper can be found here https://arxiv.org/pdf/1707.07327.pdf

Run
``
julia results/results.jl
``
and it will create the figures and tables used the paper inside the folders "results/figures" and "results/tables" respectively.

In order to run this file you need to install several packages including
- Ipopt
- OnePhase IPM. (links)

TODO:
- delete notebooks and clean up repo
- write README
- print results for small test problems (run IPOPT, etc seperately)
