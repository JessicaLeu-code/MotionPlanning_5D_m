# MotionPlanning_5D_m

This package contains motion planning examples that runs the following four solvers: PSGCFS, RRT*-CFS, CFS.

To run these examples (main files), users will need to install supporting packages. (Install other packeges suggested by MATLAB when running the examples.)  

There are three examples in this package.
- `main_FANUC.m` is a simple implementation of the Convex Feasible Set (CFS) algorithm [1] and projected stochastic gradient Convex Feasible Set (PSGCFS) in 5Dof planning.
- `main_2L.m` implements CFS and PSGCFS for 2Dof robotic arm planning.
- `RRTstar_CFS.m` implements RRT*-CFS [2] in 5Dof planning. 

![GitHub Logo](/planning.jpg)

[1] Liu, Changliu, Chung-Yen Lin, and Masayoshi Tomizuka. "The convex feasible set algorithm for real time optimization in motion planning." SIAM Journal on Control and optimization 56.4 (2018): 2712-2733.

[2] https://jessicaleu24.github.io/ACC2021.html


