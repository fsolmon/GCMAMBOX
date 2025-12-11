cmake_minimum_required(VERSION 3.18)
project(pegasus VERSION 2.6.1 LANGUAGES Fortran C)

# Set the Fortran standard
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Output directory
set(OUT_DIR "${CMAKE_SOURCE_DIR}/../modbin/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUT_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${OUT_DIR})

# Create output directory if it doesn't exist
file(MAKE_DIRECTORY ${OUT_DIR})

# Include configuration files (if they exist)
# These should define FC, CPP, CPPFLAGS, FCFLAGS, LDFLAGS, EXE
if(EXISTS "${CMAKE_SOURCE_DIR}/config.cmake")
    message(STATUS "Found config.cmake")
endif()

# Compiler settings
# You can override these with: cmake -DCMAKE_Fortran_COMPILER=gfortran ...
if(NOT CMAKE_Fortran_COMPILER_ID)
    message(WARNING "Fortran compiler not explicitly set")
endif()

#-----------------------------------------------------------------------------
# Find NetCDF Libraries
#-----------------------------------------------------------------------------
# Option 1: Use pkg-config (recommended)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(NETCDF IMPORTED_TARGET netcdf-fortran netcdf)
    if(NETCDF_FOUND)
        message(STATUS "Found NetCDF via pkg-config")
        message(STATUS "  NetCDF include dirs: ${NETCDF_INCLUDE_DIRS}")
        message(STATUS "  NetCDF library dirs: ${NETCDF_LIBRARY_DIRS}")
        message(STATUS "  NetCDF libraries: ${NETCDF_LIBRARIES}")
    endif()
endif()

# Option 2: Manual search if pkg-config fails
if(NOT NETCDF_FOUND)
    message(STATUS "Searching for NetCDF manually...")

    # Set potential search paths
    set(NETCDF_SEARCH_PATHS
        /home/fsolmon/miniconda3
        /usr
        /usr/local
        $ENV{NETCDF_DIR}
        $ENV{NETCDF_ROOT}
	$ENV{NETCDF_INCLUDE}
	$ENV{NETCDF_LIB}
        $ENV{NETCDF_FORTRAN_HOME}
    )

    # Find NetCDF C library
    find_library(NETCDF_C_LIB
        NAMES netcdf
        PATHS ${NETCDF_SEARCH_PATHS}
        PATH_SUFFIXES lib lib64
    )

    # Find NetCDF Fortran library
    find_library(NETCDF_F_LIB
        NAMES netcdff
        PATHS ${NETCDF_SEARCH_PATHS}
        PATH_SUFFIXES lib lib64
    )

    # Find NetCDF include directory
    find_path(NETCDF_INCLUDE_DIR
        NAMES netcdf.mod netcdf.inc
        PATHS ${NETCDF_SEARCH_PATHS}
        PATH_SUFFIXES include Include
    )

    if(NETCDF_C_LIB AND NETCDF_F_LIB AND NETCDF_INCLUDE_DIR)
        set(NETCDF_FOUND TRUE)
        set(NETCDF_LIBRARIES ${NETCDF_F_LIB} ${NETCDF_C_LIB})
        set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
        message(STATUS "Found NetCDF manually")
        message(STATUS "  NetCDF C library: ${NETCDF_C_LIB}")
        message(STATUS "  NetCDF Fortran library: ${NETCDF_F_LIB}")
        message(STATUS "  NetCDF include dir: ${NETCDF_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "NetCDF libraries not found. Please install NetCDF or set NETCDF_DIR environment variable.")
    endif()
endif()

# Preprocessor settings
# Set your CPP and CPPFLAGS here or via command line
set(CPP "cpp" CACHE STRING "C preprocessor")
set(CPPFLAGS "-P -traditional -DMOSAIC_SPECIES -DCAMBOX_ACTIVATE_THIS -DCAMBOX_DEACTIVATE_THIS -DMODAL_AERO -DMODAL_AERO_4MODE_MOM -DRAIN_EVAP_TO_COARSE_AERO -DUSE_NC4" CACHE STRING "C preprocessor flags")

# Fortran compiler flags
if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
    set(FCFLAGS "-g -O0 -fno-range-check -ffree-line-length-none -fallow-invalid-boz" CACHE STRING "Fortran compiler flags")
else
    set(FCFLAGS "-g -O0 -fno-range-check -ffree-line-length-none " CACHE STRING "Fortran compiler flags")
endif()

# Build the full Fortran flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FCFLAGS} ${CPPFLAGS}")

# Executable name
set(EXE "gcmambox" CACHE STRING "Executable name")

#-----------------------------------------------------------------------------
# Custom preprocessing function for Fortran files
#-----------------------------------------------------------------------------
function(preprocess_fortran input_file output_file)
    get_filename_component(file_name ${input_file} NAME_WE)
    get_filename_component(file_ext ${input_file} EXT)

    if(file_ext STREQUAL ".F" OR file_ext STREQUAL ".F90")
        # Determine output extension
        if(file_ext STREQUAL ".F")
            set(out_ext ".f")
        else()
            set(out_ext ".f90")
        endif()

        set(preprocessed_file "${OUT_DIR}/${file_name}${out_ext}")

        add_custom_command(
            OUTPUT ${preprocessed_file}
            COMMAND ${CPP} ${CPPFLAGS} ${input_file} ${preprocessed_file}
            DEPENDS ${input_file}
            COMMENT "Preprocessing ${input_file}"
        )

        set(${output_file} ${preprocessed_file} PARENT_SCOPE)
    else()
        set(${output_file} ${input_file} PARENT_SCOPE)
    endif()
endfunction()

#-----------------------------------------------------------------------------
# Source file lists
# normally the OBJ1-OBJ5 are sources that are used in the MAM lib in GeosCore
# OBJ 9 is the driver
#-----------------------------------------------------------------------------

# OBJ1: Shared modules
set(OBJ1_SOURCES
    precision_mod.F90
    physconstants.F90
    mam_utils.F90
)

# OBJ2: declaration and utilities from CAM env.
set(OBJ2_SOURCES
    constituents.F90
    radconstants.F90
    physconst.F90
    buffer.F90
    physics_buffer.F90
    physics_types.F90
    phys_control.F90
    chem_mods.F90
    wv_saturation.F90
    error_function.F90
)

# OBJ3: Modal aerosol data
set(OBJ3_SOURCES
    modal_aero_data.F90
    rad_constituents.F90
)

# OBJ5: MOSAIC and MAM aerosol modules
set(OBJ5_SOURCES
    module_data_mosaic_aero.F90
    module_mosaic_support.F90
    module_data_mosaic_asecthp.F90
    module_mosaic_ext.F90
    module_mosaic_astem.F90
    module_mosaic_init_aerpar.F90
    module_data_mosaic_main.F90
    module_mosaic_lsode.F90
    module_mosaic_box_aerchem.F90
    modal_aero_rename.F90
    modal_aero_newnuc.F90
    modal_aero_gasaerexch.F90
    modal_aero_coag.F90
    modal_aero_amicphys.F90
    module_mosaic_cam_init.F90
    modal_aero_wateruptake.F90
    modal_aero_calcsize.F90
    modal_aero_initialize_data.F90
)

# OBJ9: Driver and utilitie for box model
set(OBJ9_SOURCES
    gaschem_simple.F90
    cloudchem_simple.F90
    driver.F90
    main.F90
)

# Combine all sources
set(ALL_SOURCES
    ${OBJ1_SOURCES}
    ${OBJ2_SOURCES}
    ${OBJ3_SOURCES}
    ${OBJ5_SOURCES}
    ${OBJ9_SOURCES}
)

# Check if source files exist (with common extensions)
set(VERIFIED_SOURCES)
foreach(src_file ${ALL_SOURCES})
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${src_file}")
        list(APPEND VERIFIED_SOURCES ${src_file})
    else()
        get_filename_component(base_name ${src_file} NAME_WE)
        get_filename_component(ext ${src_file} EXT)

        # Try lowercase extension
        string(TOLOWER ${ext} ext_lower)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${base_name}${ext_lower}")
            list(APPEND VERIFIED_SOURCES ${base_name}${ext_lower})
        else()
            message(WARNING "Source file not found: ${src_file}")
        endif()
    endif()
endforeach()

#-----------------------------------------------------------------------------
# Build the executable
#-----------------------------------------------------------------------------
add_executable(${EXE} ${VERIFIED_SOURCES})

# Target properties
set_target_properties(${EXE} PROPERTIES
    Fortran_MODULE_DIRECTORY ${OUT_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${OUT_DIR}
)

#-----------------------------------------------------------------------------
# Link NetCDF libraries to the executable
#-----------------------------------------------------------------------------
if(PkgConfig_FOUND AND TARGET PkgConfig::NETCDF)
    # If using pkg-config, link with the imported target
    target_link_libraries(${EXE} PRIVATE PkgConfig::NETCDF)
else()
    # Manual linking
    target_include_directories(${EXE} PRIVATE ${NETCDF_INCLUDE_DIRS})
    target_link_libraries(${EXE} PRIVATE ${NETCDF_LIBRARIES})
endif()

# Clean targets
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${OUT_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${OUT_DIR}
    COMMENT "Cleaning all build artifacts"
)

add_custom_target(cleanmod
    COMMAND ${CMAKE_COMMAND} -E chdir ${OUT_DIR} ${CMAKE_COMMAND} -E remove *.x *.mod module*.o module*.f90 module*.f
    COMMENT "Cleaning module files"
)

#-----------------------------------------------------------------------------
# Installation (optional)
#-----------------------------------------------------------------------------
install(TARGETS ${EXE}
    RUNTIME DESTINATION bin
)

# Print configuration summary
message(STATUS "=============================================")
message(STATUS "MAMBOX v${PROJECT_VERSION} Configuration")
message(STATUS "=============================================")
message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran Flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "CPP: ${CPP}")
message(STATUS "CPP Flags: ${CPPFLAGS}")
message(STATUS "NetCDF Found: ${NETCDF_FOUND}")
if(NETCDF_FOUND)
    message(STATUS "NetCDF Include Dirs: ${NETCDF_INCLUDE_DIRS}")
    message(STATUS "NetCDF Libraries: ${NETCDF_LIBRARIES}")
endif()
message(STATUS "Output Directory: ${OUT_DIR}")
message(STATUS "Executable: ${EXE}")
message(STATUS "=============================================")
