# CMake configuration for Geant4 project: g4_pigknee

cmake_minimum_required(VERSION 3.16)

# Project name and language
project(pigknee CXX)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the Geant4 package and its components
find_package(Geant4 REQUIRED)

# Include Geant4's headers
include(${Geant4_USE_FILE})

# Add your source files to be compiled
file(GLOB sources "src/*.cc")
add_executable(pigknee ${sources})

# --- Add project's include directory ---
# This MUST come AFTER add_executable()
target_include_directories(pigknee PUBLIC include)

# Link the executable against the Geant4 libraries
target_link_libraries(pigknee ${Geant4_LIBRARIES})

# --- Installation ---
# This section defines where the executable should be installed.
install(TARGETS pigknee DESTINATION bin)
