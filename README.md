#EIT cooling of trapped ion crystals

This is a compiled semiclassical molecular-dynamics code in C++ capable of simulating electromagnetically-induced transparency (EIT) cooling of large ion crystals confined in a Penning trap.  Penning trap parameters such as the magnetic field magnitude, axial confinement frequency, and rotating wall strength and frequency can be specified by the user.  Similarly, ion mass and charge can be specified.  While this code has only been tested using particles with the same charge and mass, it should be able to handle the general case of unique charges and masses with minimal, if any, modifications.  To accelerate the calculation of Coulomb forces between ions, the code incorporates the fast multipole method using the FMM3D library. Therefore, the time to compute Coulomb forces scales approximately linearly with ion number when N is sufficiently large, making this code especially useful in studying large crystals. 

## Installation
To compile the code, it is necessary to install the FMM3D library. The Doppler cooling code has been benchmarked using v1.0.1 of the FMM3D library, although it should be possible to use more recent versions.  When including the FMM3D v1.0.1 C headers from C++ code, we found it necessary to add the following compatibility guard to FMM3D/c/utils.h:

```c
#ifdef __cplusplus
#define complex _Complex
#endif
```

Download v1.0.1 of the FMM3D library 

    git clone --branch v1.0.1 --depth 1 https://github.com/flatironinstitute/FMM3D.git FMM3D_v1.0.1

and install using 

    make install PREFIX=/path/to/FMM3D_v1.0.1_install FAST_KER=ON

Further details on the installation of FMM3D can be found at https://fmm3d.readthedocs.io/en/latest/

To compile the Doppler cooling code, first

    cd ion-crystal-eit-cooling

then 

    make

### FMM3D paths

By default, the Makefile assumes

    FMM3D_DIR=$HOME/FMM3D_v1.0.1
    FMM3D_LIB=$HOME/FMM3D_v1.0.1_install

These can be overridden at build time:

    make FMM3D_DIR=/path/to/FMM3D \
         FMM3D_LIB=/path/to/directory/containing/libfmm3d.so

At runtime, ensure that the directory containing `libfmm3d.so`
is included in `LD_LIBRARY_PATH`, for example:

    export LD_LIBRARY_PATH=/path/to/directory/containing/libfmm3d.so:$LD_LIBRARY_PATH


## Usage
Jupyter notebook files for pre- and post-processing are provided in the examples folder. Before running the Python examples, set FMM3D_PYTHON_PATH in the notebook ./examples/initialize_crystal (or ./examples/initialize_crystal_3D) to the location of your local FMM3D/python directory.

The files './examples/initialize_crystal.ipynb' and './examples/initialize_crystal_3D.ipynb' produce ion crystals whose cooling can then be simulated.  Use the former to initialize 2D (planar) crystals and the latter to initialize 3D crystals. After the user specifies particle, trap, and laser parameters, an equilibrium ion configuration is then computed and the crystal is excited to a specified initial temperature. An input text file is generated with the simulation parameters and initial crystal state.  

Several example scripts to run the EIT cooling simulation are provided in the folder ./examples/run_scripts and can be run from the repository root.  The input files required to run this script are generated in initialize_crystal.ipynb and initialize_crystal_3D.ipynb.  In the file 'run_example_100i.sh' and 'run_example_100i_3D.sh', the executable name is followed by the name of the input file and the desired name of the output file.  The code is also capable of starting a cooling simulation from the end state of a previous run.  This resume feature is shown in 'run_example_100i_3D_restart.sh', which reads the same input file, but also takes the previous output files as inputs.  

After the simulation has finished, the output files can be read using './examples/read_output.ipynb' and './examples/read_output_density_matrix.ipynb', which read the classical coordinates and the internal state density matrices, respectively.

We note that the direct Coulomb calculation can be disabled and replaced with the fast multipole method by uncommenting out 'coulomb_force_fmm();'  and commenting 'coulomb_force();' in the function 'run_trap_potential(double kz)' (both function calls appear twice).   For small crystals, the direct Coulomb force evaluation is faster than the FMM.

## License

This project is distributed under an MIT license. See the `LICENSE` file for the full license text.

Third-party software included in or used by this project remains subject to its own license terms.

### Third-party software

* **dSFMT** — The repository includes source code from the dSFMT pseudorandom-number generator under `external/dSFMT/`. These files remain under the original dSFMT BSD license. The corresponding license notice is included with the bundled dSFMT source. https://github.com/MersenneTwister-Lab/dSFMT
* **FMM3D / FMM3DPy** — Long-range Coulomb interactions may be evaluated using FMM3D, developed by the Flatiron Institute FMM3D development team. FMM3D is distributed under the Apache License 2.0. https://github.com/flatironinstitute/fmm3d
* **NumPy, SciPy, Matplotlib, and Jupyter** — These Python packages are external dependencies used for initialization, analysis, visualization, and example workflows.

## Acknowledgements

Long-range Coulomb calculations use FMM3D / FMM3DPy, developed by the Flatiron Institute FMM3D development team. Users of this repository are encouraged to cite the FMM3D project and associated publications when appropriate.

The Python analysis and initialization tools make use of NumPy, SciPy, Matplotlib, and Jupyter.

## Dependencies

### C/C++

* GCC / G++ with C++17 support
* FMM3D

### Python

* Python 3
* NumPy
* SciPy
* Matplotlib
* Jupyter
* FMM3DPy
  
## Support
If you encounter any issues or have questions, please email john.zaris@colorado.edu.





