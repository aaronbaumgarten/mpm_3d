---------------------------------------------------------------------
MPM_3D README.txt
Version 3.0 (Nov 21 2018)
---------------------------------------------------------------------


---------------------------------------------------------------------
INSTALLATION
---------------------------------------------------------------------
Installation is completed using CMake 3.2.2 with gcc 5.2.1,
though earlier versions may be supported. Requires Eigen C++
library.

Using the command line, make a build directory:
"$ mkdir build"
"$ cd build"

Then call cmake in the build directory:
"build$ cmake <path-to-README.txt>"
"build$ make"

An executable "mpm_v3" will be created in the build directory.
This executable should not be moved, but can be called from
any subfolder or read/write accessible directory.
---------------------------------------------------------------------


---------------------------------------------------------------------
CREATING SIMULATION
---------------------------------------------------------------------
For starting a new simulation, a configuration file and
.points files are needed. Example simulations can be found in
the "examples/" directory.

In general, configuration files define the various objects and
object parameters needed to run a simulation. The format is roughly
shown here:

"
    #this is a comment
    job #this is a header
    {
        #this is how values are set
        dt = 1e-3
        t = 0
        TYPE = 2
    }

    serializer
    {
        class = "DefaultVTK"
        properties = {60} #sample rate [fps]
        int-properties = {}
        #{frame directory, save directory, name}
        str-properties = {"output","save","test2D"}
    }

    driver
    {
        class = "DefaultDriver"
        properties = {1} #stop time [s]
        int-properties = {}
        str-properties = {}
    }

    solver
    {
        class = "ExplicitUSL"
        properties = {}
        #{use_gimp, implicit_contact}
        int-properties = {1, 1}
        str-properties = {}
    }

    grid
    {
        class = "CartesianLinear"
        #{Lx, Ly}
        properties = {1.0,1.0}
        #{Nx, Ny}
        int-properties = {20,20}
        str-properties = {}
    }

    body
    {
        class = "DefaultBody"
        name = "block_0"
        point-file = "block_0.points"

        point-class = "DefaultPoints"
        point-props = {}
        point-int-props = {}
        point-str-props = {}

        node-class = "DefaultNodes"
        node-props = {}
        node-int-props = {}
        node-str-props = {}
        
        material-class = "IsotropicLinearElasticity"
        #{E, nu}
        material-props = {1e5,0.3}
        material-int-props = {}
        material-str-props = {}
        
        boundary-class = "CartesianBox"
        boundary-props = {}
        boundary-int-props = {}
        boundary-str-props = {}
    }

    body
    {
        class = "DefaultBody"
        name = "block_1"
        point-file = "block_1.points"
        
        point-class = "DefaultPoints"
        node-class = "DefaultNodes"    

        material-class = "IsotropicLinearElasticity"
        material-props = {1e5,0.3}
        
        boundary-class = "CartesianBox"
    }

    contact
    {
        name = "collision"
        class = "ContactHuang"
        properties = {0.4}
        int-properties = {}
        str-properties = {"block_0","block_1"}
    }
"

In addition to the configuration file, each body needs a ".points"
file as defined in the configuration file. These files have the
following general format (space delimited):

"
    <number-of-points>
    <mass> <volume> <position-vector> <velocity-vector> <active?>
    .      .        .                 .                 .
    .      .        .                 .                 .
    .      .        .                 .                 .
    <mass> <volume> <position-vector> <velocity-vector> <active?>
"

All values in the ".points" file are space delimited and each line
represents a different point. The number of dimensions in the above
vectors should match the dimension of the simulation specified by
the "TYPE" value passed in the configuration file.

Some useful generation files are included in "mpm_pygen/" and
examples in "examples/".

To run a simulation using a configuration file,
use the following command:
"$<path-to-mpm_v3>/mpm_v3 -c <configuration-file>"
---------------------------------------------------------------------


---------------------------------------------------------------------
CUSTOMIZING MPM_3D
---------------------------------------------------------------------
MPM_3D is highly customizeable, letting you write user defined
materials, boundary conditions, contact algorithms, solvers,
serializers, drivers, grids, bodies, nodes, and points. It is 
recommended that users look at currently implemented examples to 
more fully understand how these interact with the base code.

If the code snippet needs run-time defined parameters, these should
be clearly written in the "init()" function.

Within the "src/" directory are subdirectories for each of these
objects. When adding a new file, the user must also add it to the
relevant "CMakeLists.txt" file and the registry.


---------------------------------------------------------------------
RESTARTING SIMULATIONS
---------------------------------------------------------------------
In theory, MPM_3D, can be fully serialized. If save files are written
succesfully by the serializer, then a simulation can be restarted as
follows:

"$<path-to-mpm_v3>/mpm_v3 -r <path-to-save-folder>/default_save_file.txt"
