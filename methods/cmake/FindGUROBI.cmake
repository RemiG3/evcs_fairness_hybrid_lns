#### Taken from http://www.openflipper.org/svnrepo/CoMISo/trunk/CoMISo/cmake/FindGUROBI.cmake


# - Try to find GUROBI
# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi

if (GUROBI_INCLUDE_DIR)
  # in cache already
  set(GUROBI_FOUND TRUE)
  set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
  set(GUROBI_LIBRARIES "${GUROBI_LIBRARY};${GUROBI_CXX_LIBRARY}" )
else (GUROBI_INCLUDE_DIR)



# Try to get GUROBI_HOME from CMake cache variable first, then environment variable
if(NOT GUROBI_HOME)
    set(GUROBI_HOME $ENV{GUROBI_HOME})
endif()

# If GUROBI_HOME is not set, try common installation paths
if(NOT GUROBI_HOME)
    # Common installation paths for different platforms
    if(WIN32)
        # Windows common paths
        file(GLOB GUROBI_PATHS "C:/gurobi*/win64" "C:/Program Files/gurobi*/win64")
    elseif(APPLE)
        # macOS common paths
        file(GLOB GUROBI_PATHS "/Library/gurobi*/macos_universal2" "/opt/gurobi*/macos_universal2" "$ENV{HOME}/gurobi*/macos_universal2")
    else()
        # Linux common paths
        file(GLOB GUROBI_PATHS "/opt/gurobi*/linux64" "$ENV{HOME}/gurobi*/linux64" "/usr/local/gurobi*/linux64")
    endif()
    
    # Use the first found path (usually the newest version)
    if(GUROBI_PATHS)
        list(GET GUROBI_PATHS 0 GUROBI_HOME)
        message(STATUS "Found Gurobi installation at: ${GUROBI_HOME}")
    endif()
endif()

# If still not found, show helpful error message
if(NOT GUROBI_HOME)
    message(FATAL_ERROR "GUROBI_HOME not found. Please either:\n"
                        "  1. Set the GUROBI_HOME environment variable to your Gurobi installation directory, or\n"
                        "  2. Install Gurobi in a standard location:\n"
                        "     - Linux: /opt/gurobi*/linux64 or ~/gurobi*/linux64\n"
                        "     - macOS: /Library/gurobi*/macos_universal2\n"
                        "     - Windows: C:/gurobi*/win64")
endif()

message(STATUS "Using GUROBI_HOME: ${GUROBI_HOME}")




find_path(GUROBI_INCLUDE_DIR
          NAMES gurobi_c++.h
          PATHS "${GUROBI_HOME}/include"
                "${GUROBI_HOME}/lib"
          )

find_library( GUROBI_LIBRARY
              NAMES gurobi110
                    gurobi1102
              PATHS "${GUROBI_HOME}/lib"
              )

find_library( GUROBI_CXX_LIBRARY
              NAMES gurobi_c++
              PATHS "${GUROBI_HOME}/lib"
              )

set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}" )

# use c++ headers as default
# set(GUROBI_COMPILER_FLAGS "-DIL_STD" CACHE STRING "Gurobi Compiler Flags")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GUROBI  DEFAULT_MSG
                                  GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_INCLUDE_DIR)

mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY)

endif(GUROBI_INCLUDE_DIR)
