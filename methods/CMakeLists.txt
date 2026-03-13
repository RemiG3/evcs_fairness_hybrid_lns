cmake_minimum_required(VERSION 3.5)

project(evcs_optimizer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(evcs_optimizer main.cpp lns.cpp lns.h mip.cpp mip.h utils.h utils.cpp)

# Let FindGUROBI.cmake handle all the path detection
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(GUROBI REQUIRED)

# Use the variables set by FindGUROBI.cmake
target_include_directories(evcs_optimizer PRIVATE ${GUROBI_INCLUDE_DIRS})
target_link_libraries(evcs_optimizer ${GUROBI_LIBRARIES} stdc++fs)



