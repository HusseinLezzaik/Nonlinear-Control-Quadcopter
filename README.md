# Nonlinear Control Quadcopter

## Introduction
The main objective of this project is to propose a stabilization control algorithm for a quadrotor, based on the " Nested Saturations Approach " .

## Getting Started
1.  Clone our repo: `git clone https://github.com/HusseinLezzaik/Nonlinear-Control-Quadcopter.git`
2.  Install [MATLAB](https://fr.mathworks.com/products/matlab-online.html).

## Overview of the Repository
* The MATLAB script `"main.m"` simulates the behaviour of a quadcopter 
controlled by nested saturation functions. The other two scripts are 
helper functions for quadcopter animation.

* Simulation parameters that can be changed: The first few blocks of 
   the script contain all adjustable parameters:
	- Reference Trajectory: can be `"spiral"`,`"cylinder"`,`"reference"`, or
							`"custom"`.
	- `Noise factor` on each input and on each state. Can be set to zero.
	- `Animation`: can be set on or off.
	- `Physical parameters` of the drone.
	- `Simulation time` and `sampling time`.
	- `Gains` of the control law proposed.
	- `Initial positions` and `velocities` of the quadcopter.
	- Custom `trajectory` and its derivetive, if trajectory was chosed to
	  be custom.

* After adjusting the parameters, simply run the script (`"main.m"`).

## What the script does:
	- if animation enabled, the drone and its desired trajectory are plotted
	  as animation.
	- Plots of the errors on x,y,z, and yaw angle.
	- Plots of the applied main thrust and the three torques.
	- 3D Plots of the trajectory followed, and the desired one.

## Maintainers
If you have any question, or if anything of the above is not working, don't hestitate to contact us! We are more than happy to help!
* [Hussein Lezzaik](www.husseinlezzaik.com)
* [Hasan Kassem](https://www.linkedin.com/in/hasan-kassem-02625119b/)
* [Ahmad Shour](https://www.linkedin.com/in/ahmad-shour-1531371a8/)
