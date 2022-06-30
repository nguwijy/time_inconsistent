# TICP solver
Authors: Jiang Yu Nguwi and Nicolas Privault.

If this code is used for research purposes, please cite as

J.Y. Nguwi and N. Privault.
A constructive approach to existence of equilibria in time-inconsistent stochastic control problems.
*SIAM Journal on Control and Optimization,*
60(2), 674--698, 2022.

TICP solver aims at constructing equilibria in time-inconsistent control problems (TICP).
The paper is available [here](doc/m.pdf).
<br/><br/>

## Explanation
All script files are organized as follows:
1. [equilibrium_discrete_continuous.py](equilibrium_discrete_continuous.py) is the main scripts solving linear-quadratic and linear-quartic TIC problems. It can be run by setting `LQ_or_quartic` (0 for linear-quadratic and 1 for linear-quartic) and `python equilibrium_discrete_continuous.py`.
2. [plot.R](plot/plot.R) reads the data files output by above (saved in [final_log](final_log)) and outputs the data files to be read by `gnuplot`. It can be run by setting `which_prob` (1 for linear-quadratic and 2 for linear-quartic) and `Rscript plot.R`.
3. [run_plot](plot/run_plot) is the bash script used to generate Figures 1-7 in the [paper](doc/m.pdf). It can be run by `bash run_plot`.

**Note:** The data files output by
[equilibrium_discrete_continuous.py](equilibrium_discrete_continuous.py)
are too large to be uploaded in GitHub.
They are available [here](https://drive.google.com/drive/folders/1PmMCkxUge7yS98uZWVdpCa1oKBTYymcP?usp=sharing).
