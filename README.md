# On the behavior of Lagrange multipliers in convex and non-convex infeasible interior point methods

This git repo contains code to allow one to replicate the results of the paper: "On the behavior of Lagrange multipliers in convex and non-convex infeasible interior point methods. Gabriel Haeser, Oliver Hinder, Yinyu Ye. 2018." The arxiv paper can be found here https://arxiv.org/pdf/1707.07327.pdf.

To replicate the results run
``
julia results/results.jl
``
in the console and it will create the figures and tables used the paper inside the folders "results/figures" and "results/tables" respectively.

In order to run this code you need to install several packages including
- Ipopt (https://github.com/JuliaOpt/Ipopt.jl)
- One Phase IPM (https://github.com/ohinder/OnePhase.jl)
