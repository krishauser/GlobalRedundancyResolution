# Global Redundancy Resolution

Kris Hauser, Scott Emmons
<kris.hauser@duke.edu>, <scott.emmons@duke.edu>
Duke University

A software package for computing global redundancy resolutions for redundant robots over 2D, 3D, and 6D Cartesian workspaces.  A redundancy resolution is a map from some Cartesian workspace to a robot's configuration space.  This accompanies the paper:

> K. Hauser.  Continuous Pseudoinversion of a Multivariate Function: Application to Global Redundancy Resolution. Workshop on the Algorithmic Foundations of Robotics (WAFR), 2016. 

More information and examples can be found at:

> http://motion.pratt.duke.edu/redundancyresolution/


## Prerequisites

The main requirements are Python 2.7, the Klamp't Python API, the IKDB library, and scipy.

Instructions for installing Klamp't can be found here:

> http://motion.pratt.duke.edu/klampt/tutorial_install.html

Version 0.7+ is required. Even if you are using a binary Python install, you should still install Klampt since many of the example problems refer to the robot files in the Klampt/data folder.

To install IKDB, git clone the repository as follows:

> git clone https://github.com/krishauser/ikdb.git

And then move or add a symbolic link to the ikdb/ikdb folder inside GlobalRedundancyResolution/grr, e.g., with

> cd GlobalReduncancyResolution/grr/
> ln -s ~/ikdb/ikdb 

To test whether you have correctly installed the prerequisites, in this folder run:

> python redundancy.py atlas

You should then see the ATLAS robot and a previously computed redundancy resolution.  (The default visualization assumes that you have installed Klamp't in the user's home directory; if this is not the case you can edit the line ``` klampt_directory = os.path.expand("~/Klampt") ``` in grr/visualization.py.)


## File Structure

The main file structure is as follows

- grr/: The main Python module
  - ikdb/: A link to the IK Database package, used primarily just for defining rich IK problems.
  - utils/: Various utilities.
  - graph.py: Defines the main RedundancyResolutionGraph data structure containing and querying a redundancy resolution.
  - solver.py: The redundancy resolution solver.
  - visualization.py: A basic GUI for visualizing a redundancy resolution.  Can be subclassed or run in script mode.
- output/: The program will generate output in this folder
- problems/: Problems should be put into this folder.
  - atlas/: Sample problems for the Boston Dynamics ATLAS.
  - baxter/: Sample problems for the Rethink Robotics Baxter.
  - jaco/: Sample problems for the Kinova Jaco.
  - planar_3R/: Sample problems for a planar 3R robot.
  - planar_5R/: Sample problems for a planar 5R robot.
  - planar_10R/: Sample problems for a planar 10R robot.
  - planar_20R/: Sample problems for a planar 20R robot.
  - robonaut2/: Sample problems for the NASA Robonaut2t
- shelf/: Temporary folder for persistent data stored by the solver.  Note that no more than 1 instance of the solvers should be run at the same time.
- redundancy.py: A GUI for solving a redundancy resolution problem.


## Running the main program

The command line arguments for both redundancy.py and grr/visualization.py are as follows:

> python redundancy.py PROBLEM [JSONFILE] [OPT1=VAL1] [OPT2=VAL2] [clean]...

The PROBLEM variable is some folder inside the problems/ directory.  If JSONFILE is a .json file given as the second argument, then this is the name of the file inside problems/PROBLEM that defines the parameters of the problem.

Any options specified on the command line override the definitions in the JSON file.  Valid options are specified below.

If "clean" is specified on the command line, any previously computed resolutions are not loaded from disk.  Otherwise, if a previously computed resolution is available, it will be loaded.



## Problem JSON files

These files specify the location of the robot model, the Cartesian workspace, and the number of nodes defined in the workspace.  It is probably easiest just to start with an existing file and 

Valid items are as follows:
- filename: the URDF robot, .rob Klamp't robot, or .xml Klamp't world filename.
- Nw: the number of workspace nodes
- workspace_graph: either "staggered_grid" (default), "grid", or "random".
- link: the end effector link.
- eelocal: the local coordinates of the end effector point.
- orientation: either "free", "variable", or "fixed".
- fixedRotation: if orientation is "fixed" then this is the desired orientation, in Klamp't so3 module format.  Can also be a snipped of python code.
- domain: the domain of the Cartesian task.  If orientation is "free" or "fixed", this is 3D.  If orientation is "variable" then this is 6D.
- useCollision: if true, collisions are avoided.  Otherwise, self-collisions are allowed.
- links: the links that are movable. By default this is all links of the robot.
- Nq: the number of configurations per workspace node, used for parallel computation.
- output: the output folder name (default is of the form ```output/[PROBLEM]/[JSONFILE (without extension)]_W[Nw]```)
- output_prefix: prepends a custom prefix to the output folder.

When setting up the domain and other parameters, it is useful to launch one of the visualization files and then type 'g'.  This will display the workspace graph.  It is also useful to drag around the end effector widget to see whether the IK problem is being solved correctly.


## Solver GUI commands

The solver builds a roadmap linking a workspace grid and random samples in configuration space.  The accuracy of the redundancy resolution solve is determined primarily by the number of workspace nodes N<sub>w</sub> and number of configurations sampled per node N<sub>q</sub>.  We implement three solvers, of which the first and third are the most useful.

### Basic solver

The basic solver was defined in the WAFR paper above and simply builds the roadmap connecting all neighboring pairs of configurations.  Its running time is O(DN<sub>w</sub>N<sub>q</sub><sup>2</sup>) where D is the degree of the workspace graph.  However, it can be run in an incremental manner in which a solve with low N<sub>q</sub> can be inspected, and then augmented with more samples if desired.

To interact with the solver in the GUI, type these commands:

- 'q': Samples 10 more configurations and tests their connectivity.
- 'u': Solves for a redundancy resolution using a heuristic CSP assignment.
- 'd': Improves a redundancy resolution using randomized descent.
- 'r': Performs a random assignment and then uses randomized descent.
- 'c': Samples more configurations only along conflicting edges.

### Visibility graph solver

The visibility graph solver uses the notion of a self-motion manifold to reduce the practical impact of the N<sub>q</sub><sup>2</sup> factor to N<sub>q</sub>.  In typical cases where the self-motion manifolds are highly connected, the running time can be as low as O(N<sub>w</sub>N<sub>q</sub>+DN<sub>w</sub>).  However, this solver does not always achieve this best-case scenario, and is not incremental.  In other words, it uses a fixed sampling N<sub>q</sub> and cannot be inspected and revised.  It is recommended to set N<sub>q</sub> to 100 or more for this solver.

This solver also uses a more complex resolution technique that is more expensive and less customizable than the basic solver.

To run this solver, type these commands:

- 'v': compute the visibility graph roadmap
- 'R': resolve the visibility graph roadmap.

### Parallel visibility graph solver

The parallel visibility graph solver uses multithreading to speed up the solve.  However, like the visibility graph solver it is not incremental, and uses a fixed sampling Nq.  It is recommended to set Nq to 100 or more for this solver.

To run this solver, type these commands:

- 't': compute and resolve the visibility graph.  You will be prompted to enter the number of threads.

### Smoothing and Saving

Other miscellaneous keyboard commands include:
- 'o': smooth the resolution using 10 iterations of coordinate descent optimization.
- 's': save the resolution (and current solver progress) to disk.


## Visualization GUI commands

The following commands are used both in the visualization GUI and the solver GUI.

After a redundancy resolution has been computed, clicking and dragging on the widget will use the resolution to find the instantaneous solution.

Keyboard commands: 

- 'g': draw the workspace graph.
- 'x': validate the existing redundancy resolution, deleting any invalid edges.
- 'b': toggle drawing the external boundaries of the workspace (default off for >=3D, on for 2D).
- 'w': perform a random walk on the redundancy resolution.
- 'i': toggle "inspection" mode.
- 'm': begin saving movie frames to disk (in PPM format).
- 'M': save a 360 degree spin-around movie to disk (in PPM format).


## Using resolutions in external projects

Once a resolution is solved, you may use the RedundancyGraph object for use in your own external project. You will need to copy the problem settings JSON file and "resolution_graph.pickle" file found in the appropriate output directory to whereever you need it.  You may also find the parse_args() and make_program() functions in grr/visualization.py to be useful for setting up the appropriate IK problem.


