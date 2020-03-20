# MCAC

MCAC is a Monte-Carlo Aggregation Code.

It's purpose is to produce physical aggregates.

The process goes as follow :
 1. We initialize some spheres in a box (random position and size)
    At first each sphere is one aggregate
 3. We pick randomly an aggregate
 4. We move the aggregate in a random direction for
    * its mean free path
    * or until there is a collision
 5. In case of a collision, we merge the aggregates
 6. loop to 2. until we meet some criteria

Developers are invited to read the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines.

WARNING the documentation is outdated and is being rewriting  
WARNING the the code installation process is under refactoring, and may not work currently

Scientific papers are on their way, based on the [JCIS2020](https://gitlab.coria-cfd.fr/MCAC/MCAC/tree/JCIS2020) version

## Dependencies

We need
 * `libhdf5` for writing our results
 * `mpfr` and `gmp` for accurate computation of sphere overlapping (Optionnal)
 * `numpy`, `pandas`, `dask`, `h5py`, `scipy`, `matplotlib`, `scikit-image` and `cython` for python post-processing

The other dependencies that may be needed will be compiled with the code
 * `libxdmf` for writing our results, but this will need `boost`
 * `CGAL` and `sbl` for accurate computation of sphere overlapping (Optionnal)

Finally two other libraries will be compiled with the code:
 * `inipp` for reading input files
 * `gsl` for internal use

## Installation

Download the current version of the code

    git clone git@gitlab.coria-cfd.fr:MCAC/MCAC.git
    
Compile the code (It will automatically download and compile some other dependencies)

    mkdir build
    cd build
    cmake ..
    make -j4
    cd ..

If you want the SBL part, replace the `cmake` line with

    cmake -DWITH_SBL=ON ..

If you will use the python post processing tools, you have to install them as well

    pip install -e .

The "-e" part allows you to modify this library without the need to reinstall it
(except for the compiled part)

### Computer specific instructions

#### Newton

You will have to enable the correct compilation environment before the cmake part

    scl enable devtoolset-9 bash

## Update

Just go to the MCAC folder and use git to obtain the last version

    git pull
    
And recompile the compiled part that may have change

    cd build
    make -j4
    cd ..

If you experience troubles, you can alway destroy the `build` dir and restart the compilation from scratch 

And optionnaly recompile the python post-processing tools

    pip install -e .

## Usage

The binaries are produced in the bin folder (release and debug).
The code need an input file containing:

    1000        N           []                  Initial number of spheres
    1000        FV          [ppm]               Volume fraction
    30          Dpm         [nm]                Initial diameter of spheres
    1           Mode        []                  Type of dispersion of the initial diameter of spheres 1=Normal,0=LogNormal
    1.25        sigmaDpm    [sigmageo] or [nm]  Dispersion of the initial diameter of spheres (depend on the dispersion type)
    1                       []                  Physical model (0 or 1)
    1700        T           [K]                 Temperature
    101300      P           [Pa]                Pressure
    1800        Rho         [kg/m3]             Density
    1.4         kfe         []                  Fractal law linear                             (to be converged)
    1.8         dfe         []                  Fractal law power law                          (to be converged)
    1e-3        coeffB      []                  Surface Growth linear                          (stick to zero)
    2           xsurfgrowth []                  Surface Growth power law
    0                       []                  Temporal variation of physical parameters       (stick to zero)
    1           Aggmin      []                  Final number of Aggregates
    -1          WaitLimit   []                  Stop if WaitLimit iterations without collision (negative value to deactivate)
    60          CPULimit    [s]                 Stop when computation time is longer           (negative value to deactivate)
    10          DeltaSauve  []                  Number of time step per file
    output_dir  output_dir  []                  Output directory

Every element of the first column is mandatory in this order.  
Everything else is just comments and are optional.  
The second column is the name of the parameter.  
The third column is the unit if the parameter.

You then juste have to use the executable with this parameters file

    bin/release/MCAC params.txt
    
Again, on Newton you have to activate the correct environment before

    scl enable devtoolset-9 bash

## Output
### Terminal

On the screen, it will output after every collision:
- `NAgg` is the current number of aggregate, at first it is the number of spheres and decrease one by one
- `Time` is the physical time (at least statisticly)
- `CPU` is the computational time (in seconds)
- `after n it` is the nuber of time iteration without collision that precedes this one

It will also produce some **files.h5** and **files.xmf**

### Paraview

The xmf files are natively readable by paraview
To see the spheres:
0. apply
1. add a **glyph** filter
2. glyph type: **sphere**
3. advanced : glyph radius: **1**
4. scaling mode: **scalar**
5. scale factor: **1**
6. apply
7. reset camera (in order to auto scale)

### Python

A library is available in the pymcac folder
It is installable with (in the root folder of MCAC)

    pip install -e .
