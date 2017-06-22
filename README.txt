---------------------------------------------------------------------
MPM_V2 README.txt
Version 2.0 (Jun 22 2017)
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

An executable "mpm_v2" will be created in the build directory.
This executable should not be moved, but can be called from
any subfolder or read/write accessible directory.
---------------------------------------------------------------------


---------------------------------------------------------------------
CREATING SIMULATION
---------------------------------------------------------------------
For starting a new simulation, a configuration file and
.point files are needed. Example simulations can be found in
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
        DIM = 2
    }

    serializer
    {
        #filepath = "src/serializers/"
        filename = "default_vtk.so"
        properties = {60} #sample rate [fps]
        int-properties = {}
        #{frame directory, save directory, name}
        str-properties = {"output","save","test2D"}
    }

    driver
    {
        #filepath = "src/drivers/"
        filename = "default_driver.so"
        properties = {1} #stop time [s]
        int-properties = {}
        str-properties = {}
    }

    solver
    {
        #filepath = "src/solvers/"
        filename = "explicit_usl.so"
        properties = {}
        int-properties = {}
        str-properties = {}
    }

    grid
    {
        #filepath = "src/grids/"
        filename = "cartesian.so"
        #{Lx, Ly}
        properties = {1.0,1.0}
        #{Nx, Ny}
        int-properties = {20,20}
        str-properties = {}
    }

    body
    {
        name = "body_0"
        point-file = "body_0.points"
        
        #material-filepath = "src/materials/"
        material-filename = "isolin.so"
        #{E, nu}
        material-props = {1e5,0.3}
        material-int-props = {}
        material-str-props = {}
        
        #boundary-filepath = "src/boundaries/"
        boundary-filename = "cartesian_box.so"
        boundary-props = {}
        boundary-int-props = {}
        boundary-str-props = {}
    }

    body
    {
        name = "body_1"
        point-file = "body_1.points"
        material-filename = "isolin.so"
        material-props = {1e5,0.3}
        boundary-filename = "cartesian_box.so"
    }

    contact
    {
        name = "collision"
        filename = "contact_huang.so"
        properties = {0.4}
        int-properties = {}
        str-properties = {"body_0","body_1"}
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
vectors should match the "DIM" value passed in the configuration
file.

Some useful generation files are included in "mpm_pygen/" and
examples in "examples/".

To run a simulation using a configuration file,
use the following command:
"$<path-to-mpm_v2>/mpm_v2 -c <configuration-file>"
---------------------------------------------------------------------


---------------------------------------------------------------------
CUSTOMIZING MPM_V2
---------------------------------------------------------------------
MPM_V2 is highly customizeable, letting you write user defined
materials, boundary conditions, contact algorithms, solvers,
serializers, drivers, and grids. It is recommended that users
look at currently implemented examples to more fully understand how
these interact with the base code.

If the code snippet needs run-time defined parameters, these should
be clearly written in the "<object>Init()" function.

Within the "src/" directory are subdirectories for each of these
obejcts. When adding a new file, the user must also add it to the
relevant "CMakeLists.txt" file.


---------------------------------------------------------------------
RESTARTING SIMULATIONS
---------------------------------------------------------------------
In theory, MPM_V2, should be fully serialized (assuming user defined
files are written appropriately). If save files are written
succesfully by the serializer, then a simulation can be restarted as
follows:

"$<path-to-mpm_v2>/mpm_v2 -r <path-to-save-folder>/default_save_file.txt"
