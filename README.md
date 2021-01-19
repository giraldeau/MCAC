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
WARNING validation process and unit testing is under construction

## Scientific papers

  * **José Morán, Alexandre Poux, Jérôme Yon,** Impact of the competition between aggregation and surface growth on the morphology of soot particles formed in an ethylene laminar premixed flame  
    *Journal of Aerosol Science, 2020.* ([10.1016/j.jaerosci.2020.105690](https://doi.org/10.1016/j.jaerosci.2020.105690)).
 * **J. Morán, J. Yon, A. Poux, F. Corbin, F.-X. Ouf, et al..** Monte Carlo Aggregation Code (MCAC) Part 2: Application to soot agglomeration, highlighting the importance of primary particles.  
   *Journal of Colloid and Interface Science, Elsevier, 2020, 575, pp.274-285. ⟨[10.1016/j.jcis.2020.04.085](https://dx.doi.org/10.1016/j.jcis.2020.04.085)⟩.* ⟨[hal-02563051](https://hal.archives-ouvertes.fr/hal-02563051)⟩
 * **J. Morán, J. Yon, A. Poux.** Monte Carlo Aggregation Code (MCAC) Part 1: Fundamentals.  
   *Journal of Colloid and Interface Science, Elsevier, 2020, 569, pp.184-194.* ⟨[10.1016/j.jcis.2020.02.039](https://dx.doi.org/10.1016/j.jcis.2020.02.039)⟩. ⟨[hal-02494461](https://hal.archives-ouvertes.fr/hal-02494461)⟩

The corresponding code versions are [tagged](https://gitlab.coria-cfd.fr/MCAC/MCAC/-/tags).

## Dependencies

We need
 * `libhdf5` for writing our results
 * `mpfr` and `gmp` for accurate computation of sphere overlapping (Optionnal)
 * `python3.6+` for python post-processing (Optionnal)

The other dependencies that may be needed will be compiled with the code
 * `libxdmf` for writing our results, if not found it will be compiled this will need `boost`
 * `CGAL` and `SBL` for accurate computation of sphere overlapping (Optionnal)  
    (`SBL` being a large repository, you can clone it separatly and give the path to cmake with `-DSBL_GIT_REP=path`)
 * multiple python packages for python post-processing (Optionnal)

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

## Update

Just go to the MCAC folder and use git to obtain the last version

    git pull
    
And recompile the compiled part that may have change

    cd build
    make -j4
    cd ..

If you experience troubles, you can alway destroy the `build` dir and restart the compilation from scratch 

If you have updated your python version, you may have to remove the folder venv

## Usage

### MCAC

The binaries are produced in the bin folder (Release or Debug).
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

### Python

A library is available in the pymcac folder.  
It is prepared when you compile the code in a virtual environment  
In order to use `pymcac` you need to activate the virtual environment

    source venv/bin/activate
    
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

## Computer specific instructions

### Newton

You need the correct compilation environment for the compilation and execution of MCAC  
You can activate it with

    scl enable devtoolset-9 bash

If you need SBL, you'll have to add this option to cmake

     cmake .. -DBoost_INCLUDE_DIR=/usr/include/boost169

### Myria

You'll need the following module

    module load cmake python3/3.6.1 libxml2 compilers/gnu/7.3.0 
    
You'll need to export the use of gcc

    export CC=gcc
    export CXX=g++
    
You'll need to add the following option to cmake

     cmake .. -DCMAKE_PREFIX_PATH="${LIBXML2_ROOT};/soft/library/hdf5-1.8.18-gnu-serial/" 

