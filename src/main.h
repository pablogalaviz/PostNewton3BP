#include "evolution.h"
#include "initial_data.h"
#include "output.h"
#include "utils.h"
/*

   Title: Instalation


   Download the code from:

> git clone https://github.com/pablogalaviz/PostNewtonian3BP.git



   Function: main

   Controls the main loop of the code.

   Parameters:

      help,h				- Shows a help message (type: flag | default: false).
      input_file 			- Is a <Parameter File> (type: string | default None).
      silent,s				- Starts the program in silent mode (type: flag | default: false).
      					  Output to the log file exclusively.  
      debug,d 				- Shows additional debug messages in log (type: flag | default: false).
      domain.size_x			- Numerical domain size in x direction  (type: positive float | default: 1).
      					  Domain extends from -size_z to size_x.
      domain.size_y			- Similar to domain.size_x but for y direction (type: float | default: inherits from domain.size_x). 
      domain.size_z			- Similar to domain.size_x but for z direction (type: float | default: inherits from domain.size_y). 
      grid.size_x			- Grid size in x direction  (type: integer > 10 X number of cpus | default: 33).
      					  The numerical domain is partitioned by the number of cpus.
					  The number of grid points is adjusted to the multigrid requirements.  
      grid.size_y			- Similar to grid.size_x but for y direction (type: integer | default: inherits from grid.size_x).
      grid.size_z			- Similar to grid.size_x but for z direction (type: integer | default: inherits from grid.size_y).
      grid.order			- Finite differences stencil order (type: integer, {2,4,6,8} | default: 2).
      elliptic.verbose		     	- Shows information about the multigrid solver progress (type: boolean | default: false). 
      elliptic.pre_cycles	     	- Number of iteration in the pre-relaxation phase of the V-Cycle (type: integer | default: 5).
      elliptic.post_cycles	     	- Number of iteration in the post-relaxation phase of the V-Cycle (type: integer | default: 8).
      elliptic.fine_grid_tolerance   	- The maximum tolerance of the residual in the finer level of the multigrid method (type: float | default: 1e-6). 
      elliptic.coarse_grid_tolerance 	- The maximum tolerance of the residual in the coarser level of the multigrid method (type: float | default: 1e-6).
      hyperbolic.verbose		- Shows information about the ODE solver progress (type: boolean | default: false). 
      hyperbolic.ode_method		- ODE method (type: string, {rk2,rk4,rkck,rk8pd,rk2imp,rk4imp,bsimp} | default: rkck). *missing reimplementation*
      hyperbolic.initial_dt		- Initial time-step (type: float | default 1e-6).
      hyperbolic.eps_rel		- Relative error tolerance (type: float | default 1e-6).
      hyperbolic.eps_abs		- Absolute error tolerance (type: float | default 1e-6).
      hyperbolic.final_time		- Final evolution time (type: float | default 1).
      boundary.field			- A list of all the fields. The order defines the boundary property applied to each field (type: list of string | default: None). 
      boundary.type			- A list of boundary type, (type: list of strings, {PERIODIC, ZERO, DIRICHLET, ROBIN} | default: None)
      boundary.asymptote		- A list of asymptotic values for all the robin boundary conditions. (type: float | default: None)
      boundary.power			- A list of inverse power values for all the robin boundary conditions. (type: positive float | default: None)
      boundary.dirichlet_field		- A list of parametric fields for dirichlet boundary condition. (type: string | default: None)


   Returns:

      Zero if the program terminates normally. One otherwise 

   See Also:



    Title: Parameter File

    A parameters file has the following format:

>[grid]
>
> size_x = <int>  
> size_y = <int>
> size_z = <int>
>
> order = <int>
>
>
>[domain]
>
> size_x = <float>
> size_y = <float>
> size_z = <float>
>
>    
>[problem]
>
> name = <string>  //example waveEquation 
>
>    
>[output]
>
> directory=waveData
>
> fields=u
> fields=v
> .
> .
> .
> fields=V
>
> verbose=no
> delta_time=0
> delta_time_analysis=0
>
>
>[hyperbolic]
>
> verbose=no
>
> final_time=2
>
> initial_dt=0.00983
>
>
>[boundary]
>
> field=u
> type=PERIODIC
> dirichlet_field=U
>
> field=v
> type=PERIODIC
> dirichlet_field=V
>
>

	*Note:* Every parameter is optional. The values in the file are overwritten by the command line values. 

*/
      
