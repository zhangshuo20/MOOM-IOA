# MOOM-IOA

A Multi-Objective Optimization Model based on the Input-Output Analysis for the structural adjustment of economic sectors with fine-resolution power sectors is built to coordinate contradictions caused by diverse goals among economy, environment, and society.

mop_eps.gms - This program file is to realize the multi-objective optimization model proposed by this article. The solution method of the multi-objective optimization model is the augmented ε-constraint method.

MCMD.m, TOPSIS.m - These grograms are the final solutions screening method to select the solutions of Pareto optimal frontier.

tst.xlsx - This data file contains all input data,  where the A_2017, F_2020, CI, EI, LI, output_2020, lim_energy, lim_emp, rGDP, x0, x2017 represent the input-output technogical coefficient, total final demand in 2020, carbon emissions intensity, energy consumption intenisty, employment intensity, total output in 2020, the limit amount of energy consumption, the limit amount of labour force, the growth rate of GDP, the sectoral output in 2020, and the sectoral output in 2017, respectively.
