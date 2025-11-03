# CMake Build System for Pegasus 2.6.1

This CMake configuration converts your original Makefile into a modern CMake build system.

## Quick Start

### Basic Build

```bash
# Create a build directory
mkdir build
cd build

# Configure
cmake ..

# Build
cmake --build .

# Or use make
make
```

### Custom Configuration

You can override compiler settings and flags:

```bash
cmake .. \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DFCFLAGS="-O2 -g" \
  -DCPPFLAGS="-DLINUX" \
  -DLDFLAGS="-L/usr/local/lib" \
  -DEXE="my_program.x"
```

## Configuration Options

The following CMake variables can be set:

- `CMAKE_Fortran_COMPILER`: Fortran compiler (e.g., gfortran, ifort)
- `CPP`: C preprocessor (default: cpp)
- `CPPFLAGS`: Preprocessor flags
- `FCFLAGS`: Fortran compiler flags
- `LDFLAGS`: Linker flags
- `EXE`: Executable name (default: pegasus.x)

## Using Configuration Files

If you have your original `cambox_config.make.in` and `cambox_config.cpp.in` files:

1. You can convert them to a CMake cache initialization file:

```cmake
# Create a file named config.cmake
set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "")
set(FCFLAGS "-O2 -march=native" CACHE STRING "")
set(CPPFLAGS "-DLINUX -DDEBUG" CACHE STRING "")
set(LDFLAGS "-lnetcdf -lnetcdff" CACHE STRING "")
```

2. Then use it:
```bash
cmake .. -C config.cmake
```

## Build Targets

### Main Target
- `make` or `cmake --build .`: Build the executable

### Clean Targets
- `make clean`: Standard CMake clean
- `make clean-all`: Clean everything including output directory
- `make cleanmod`: Clean only module files (equivalent to `make cleanmod` in original Makefile)

## Output Directory

The executable and intermediate files are placed in `../runUPD` relative to your source directory, matching the original Makefile behavior.

## Source File Organization

The build system expects source files organized as:

- **OBJ1**: Shared modules (shr_kind_mod, shr_const_mod)
- **OBJ2**: Core modules (cam_logfile, spmd_utils, etc.)
- **OBJ3**: Modal aerosol data
- **OBJ4**: Radiation constituents
- **OBJ5**: MOSAIC and modal aerosol modules
- **OBJ9**: Driver and main program

## Preprocessing

The original Makefile uses CPP preprocessing for `.F` and `.F90` files. Currently, this CMake configuration relies on CMake's native Fortran preprocessing support. If you need custom preprocessing:

1. Ensure your Fortran compiler supports preprocessing (most modern compilers do)
2. Add appropriate flags (e.g., `-cpp` for gfortran, `-fpp` for Intel)

Example:
```bash
cmake .. -DFCFLAGS="-cpp -O2"
```

## Troubleshooting

### Missing Configuration Files
If you get warnings about missing `cambox_config.make.in` or `cambox_config.cpp.in`, you can:
1. Create a `config.cmake` file with your settings
2. Pass flags directly via command line
3. Edit `CMakeLists.txt` to set default values

### Module Dependencies
The CMakeLists.txt maintains the dependency order from your Makefile:
OBJ1 → OBJ2 → OBJ3 → OBJ4 → OBJ5 → OBJ9

### Missing Source Files
If you get warnings about missing source files, check that:
1. The file extensions match (.F, .f, .F90, .f90)
2. Files are in the same directory as CMakeLists.txt
3. File names match exactly (case-sensitive on Linux)

## Advantages of CMake

- **Parallel builds**: Use `cmake --build . -j8` for faster compilation
- **Out-of-source builds**: Keep source directory clean
- **Better dependency tracking**: Automatic recompilation when needed
- **Cross-platform**: Works on Linux, macOS, Windows
- **IDE integration**: Generate project files for Visual Studio, Xcode, etc.
- **Modern tooling**: Better integration with debuggers and profilers

## Converting Your Config Files

If you share your `cambox_config.make.in` and `cambox_config.cpp.in` files, I can help create a proper CMake configuration script that imports those settings automatically.

## Example Full Build

```bash
# From the directory containing CMakeLists.txt
mkdir build
cd build

# Configure with your preferred compiler and flags
cmake .. \
  -DCMAKE_Fortran_COMPILER=ifort \
  -DFCFLAGS="-O3 -xHost -fp-model precise" \
  -DCPPFLAGS="-DLINUX -DMAM_CONFIG=3" \
  -DLDFLAGS="-lnetcdff -lnetcdf"

# Build with 8 parallel jobs
cmake --build . -j8

# The executable will be in ../runUPD/pegasus.x
```

## Installation

To install the executable to a system location:

```bash
cmake --install . --prefix /usr/local
```

This installs the executable to `/usr/local/bin/`.
