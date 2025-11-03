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

# Preprocessor settings
# Set your CPP and CPPFLAGS here or via command line
set(CPP "cpp" CACHE STRING "C preprocessor")
set(CPPFLAGS "-P -traditional -DMOSAIC_SPECIES  -DCAMBOX_ACTIVATE_THIS -DCAMBOX_DEACTIVATE_THIS -DMODAL_AERO -DMODAL_AERO_4MODE_MOM -DRAIN_EVAP_TO_COARSE_AERO  -DPCNST=54 -DPCOLS=1 -DPVER=1 -DNBC=1 -DNPOA=1 -DNSOA=1"  CACHE STRING "C preprocessor flags")

# Fortran compiler flags
# Set your FCFLAGS here or via command line
set(FCFLAGS "-g -O0  -fno-range-check -ffree-line-length-none -I/usr/local/netcdf-4.6.1-4.4.4-gnu-7.4.1/include/" CACHE STRING "Fortran compiler flags")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FCFLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FCFLAGS} ${CPPFLAGS}")
# Linker flags
set(LDFLAGS "-L/usr/local/netcdf-4.6.1-4.4.4-gnu-7.4.1/lib -lnetcdff -lnetcdf" CACHE STRING "Linker flags")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LDFLAGS}")

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
#-----------------------------------------------------------------------------

# OBJ1: Shared modules
set(OBJ1_SOURCES
    shr_kind_mod.F90
    shr_const_mod.F90
    precision_mod.F90
)

# OBJ2: Core modules
set(OBJ2_SOURCES
    cam_logfile.F90
    spmd_utils.F90
    abortutils.F90
    cam_abortutils.F90
    ppgrid.F90
    pmgrid.F90
    constituents.F90
    radconstants.F90
    physconst.F90
    cam_history.F90
    units.F90
    ref_pres.F90
    dyn_grid.F90
    buffer.F90
    physics_buffer.F90
    physics_types.F90
    phys_control.F90
    mo_constants.F90
    chem_mods.F90
    mo_tracname.F90
    mo_chem_utls.F90
    gffgch.F90
    wv_saturation.F90
    error_function.F90
    shr_spfn_mod.F90
    infnan.F90
    seasalt_model.F90
    time_manager.F90
    aerodep_flx.F90
    modal_aero_convproc.F90
    modal_aero_deposition.F90
)

# OBJ3: Modal aerosol data
set(OBJ3_SOURCES
    modal_aero_data.F90
)

# OBJ4: Radiation constituents
set(OBJ4_SOURCES
    rad_constituents.F90
)

# OBJ5: MOSAIC and modal aerosol modules
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

# OBJ9: Driver and main
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
    ${OBJ4_SOURCES}
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

# Set up dependencies between object groups
# OBJ1 -> OBJ2 -> OBJ3 -> OBJ4 -> OBJ5 -> OBJ9

# Target properties
set_target_properties(${EXE} PROPERTIES
    Fortran_MODULE_DIRECTORY ${OUT_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${OUT_DIR}
)

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
message(STATUS "Linker Flags: ${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "Output Directory: ${OUT_DIR}")
message(STATUS "Executable: ${EXE}")
message(STATUS "=============================================")
