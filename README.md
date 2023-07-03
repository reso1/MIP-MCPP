# MIP-MCPP
This repository is the implementation of the MIP, MIP-PRH and MIP-SRH models for the rooted min-max tree cover problem and their corresponding planners for the graph-based multi-robot coverage path planning problem from the following paper:

*Mixed Integer Programming for Time-Optimal Multi-Robot Coverage Path Planning with Efficient Heuristics. Jingtao Tang and Hang Ma. * [[paper]](https://arxiv.org/abs/2306.17609)

Please cite this article if you use this code for the multi-robot coverage path planning problem.

## Installation
### 1. Python lib:
`pip install -r requirements.txt`

### 2. Gurobi lib:
> optional if you don't want to run solver.py for MIP optimization. Pre-run model solutions are provided in directory 'data/solutions'.

Please refer to [[Gurobi]](https://www.gurobi.com/) for the installation. (they have trial and academic licenses)

## Description

### 1. The RMMTC MIP Solver
> The MCPP problem is reduced to the RMMTC problem and then solved with the STC algorithm. Please refer to the paper for more details.

#### Usage
```bash
python solver.py [-h] [--solver_cfg SOLVER_CFG] [--alpha ALPHA] [--beta BETA] [--warm_start WARM_START] istc
```
- Required:
  - `istc`: the instance name stored in directory 'data/instances'.
- Optional:
  - `--solver_cfg SOLVER_CFG`: path to the Gurobi configuration file. (see 'data/cfgs' for reference)
  - `alpha ALPHA`: parameter of Parabolic Removal Heuristics (PRH). Will solve the MIP-PRH model if specified.
  - `beta BETA`: parameter of Subgraph Removal Heuristics (SRH). Will solve the MIP-SRH model if specified.
  - `--warm_start WARM_START`: type of warm-startup for the model optimization. Use 'RTC' for the original MIP model and 'MST' for MIP-PRH and MIP-SRH.

### 2. The Instance Maker
A simple routine to create random RMMTC instance.
- if map is provided, then a terrain with uniform terrain weight of 1 is generated, encoded by:
  - obstacle vertex: black pixel, rgb=(0,0,0)
  - free vertex: white pixel, rgb=(1,1,1)
  - root vertex: red pixel, rgb=(1,0,0)
- otherwise, an empty terrain with random weights will be generated.

#### Usage
```bash
python instance_maker.py [-h] [--map MAP] name
```

- Required:
  - `name`: the instance name in the format of '[grid width]x[grid height]-[Characteristics]-k[# of roots]'.
    - If no map is provided, the generated instance is a `[grid width]`x`[grid height]` empty terrain with `[# of roots]` subtrees (or robots) and randomized terrain weights.

- Optional:
  - `--map MAP`: path to the map to create the instance.

### 3. The MCPP Planner
The MCPP planners with simulation, including MFC, MSTC$^*$ and MIP (the method in this paper).

#### Usage:
```bash
python planner.py [-h] [--method METHOD] [--istc_sol_name ISTC_SOL_NAME] [--scale SCALE] [--dt DT] [--write WRITE] istc
```
- Required:
  - `istc`: the instance name stored in directory 'data/instances'.
- Optional:
  - `--method METHOD`: planner type choose from {MFC, MSTC*, MIP}.
  - `-istc_sol_name ISTC_SOL_NAME`: MIP solution name stored in the directroy 'data/solutions'. (only required when planner type is MIP)
  - `--scale SCALE`: the canvas scaling factor for visualization.
  - `--dt DT`: delta time of simulation.
  - `--write WRITE`: is writing the simulation as MP4. (ffmpeg lib is required)


## MCPP Simulation Results
- The floor-medium instance using the RMMTC solution from MIP-SRH model

![](figs/floor-medium-MIP.gif)

## License
MIP-MCPP is released under the GPL version 3. See LICENSE.txt for further details.
