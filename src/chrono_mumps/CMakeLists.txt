#=============================================================================
# CMake configuration file for the Chrono MUMPS module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================

option(ENABLE_MODULE_MUMPS "Enable the Chrono MUMPS module" OFF)

if(NOT ENABLE_MODULE_MUMPS)
	# mark_as_advanced(FORCE MUMPS_DIR)
  return()
endif()

message(STATUS "\n==== Chrono Mumps module ====\n")

enable_language(Fortran)

# mark_as_advanced(CLEAR MUMPS_DIR)

find_package(MUMPS CONFIG REQUIRED)
			 
if(NOT MUMPS_d_FOUND)
    message(FATAL_ERROR "MUMPS library has been found, but not the double precision version. Cannot continue.")
    return()
endif()

message(STATUS "MUMPS version:          ${MUMPS_VERSION}")

# ------------------------------------------------------------------------------
# List all files in the Chrono mumps module
# ------------------------------------------------------------------------------

set(ChronoEngine_MUMPS_HEADERS
  ChApiMumps.h
  ChSolverMumps.h
  ChMumpsEngine.h
)

set(ChronoEngine_MUMPS_SOURCES
  ChSolverMumps.cpp
  ChMumpsEngine.cpp
)

source_group("" FILES ${ChronoEngine_MUMPS_HEADERS} ${ChronoEngine_MUMPS_SOURCES})

# ------------------------------------------------------------------------------
# Add the ChronoEngine_mumps library
# ------------------------------------------------------------------------------

# add_definitions("-D_OPENMP_NOFORCE_MANIFEST")

add_library(ChronoEngine_mumps ${ChronoEngine_MUMPS_SOURCES} ${ChronoEngine_MUMPS_HEADERS})
target_compile_definitions(ChronoEngine_mumps PUBLIC CH_API_COMPILE_MUMPS)


set_target_properties(ChronoEngine_mumps PROPERTIES
                      COMPILE_FLAGS "${CH_CXX_FLAGS} ${CH_MUMPS_COMPILE_FLAGS}"
                      LINK_FLAGS "${CH_LINKERFLAG_LIB} ${CH_MUMPS_LINK_FLAGS}")

target_link_libraries(ChronoEngine_mumps ChronoEngine MUMPS::MUMPS)		  
                            
install(TARGETS ChronoEngine_mumps
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(FILES ${ChronoEngine_MUMPS_HEADERS} 
        DESTINATION include/chrono_mumps)