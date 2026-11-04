cmake_minimum_required(VERSION 3.13.4) # cmake version
project(cp2-feynman-ppiep C) # the project (NAME LANGUAGE)

set(CMAKE_C_STANDARD 11) # c standard

# include directories
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
find_package(fftw3 REQUIRED)
include_directories(${FFTW3_INCLUDE_DIRS})
link_directories(/usr/local/Cellar/fftw/3.3.10/lib) # apparently needed for mac cellar installations?

add_subdirectory(modules)
add_subdirectory(tests)
include_directories(${PROJECT_SOURCE_DIR}/modules)



# executables (NAME SOURCE SOURCE ...)
add_executable(groundstate groundstate.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
    ${PROJECT_SOURCE_DIR}/modules/conjugategradient.c
    ${PROJECT_SOURCE_DIR}/modules/geometry.c
    ${PROJECT_SOURCE_DIR}/modules/hamiltonian.c
    ${PROJECT_SOURCE_DIR}/modules/powermethod.c
    ${PROJECT_SOURCE_DIR}/modules/vectorutils.c
)
target_link_libraries(groundstate m)

add_executable(groundstate_v2 groundstate_v2.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
    ${PROJECT_SOURCE_DIR}/modules/conjugategradient.c
    ${PROJECT_SOURCE_DIR}/modules/geometry.c
    ${PROJECT_SOURCE_DIR}/modules/hamiltonian.c
    ${PROJECT_SOURCE_DIR}/modules/observables.c
    ${PROJECT_SOURCE_DIR}/modules/powermethod.c
    ${PROJECT_SOURCE_DIR}/modules/vectorutils.c
)
target_link_libraries(groundstate_v2 m)



add_executable(evolution evolution.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
    ${PROJECT_SOURCE_DIR}/modules/conjugategradient.c
    ${PROJECT_SOURCE_DIR}/modules/geometry.c
    ${PROJECT_SOURCE_DIR}/modules/hamiltonian.c
    ${PROJECT_SOURCE_DIR}/modules/observables.c
    ${PROJECT_SOURCE_DIR}/modules/powermethod.c
    ${PROJECT_SOURCE_DIR}/modules/timeevolution.c
    ${PROJECT_SOURCE_DIR}/modules/vectorutils.c
)
target_link_libraries(evolution m -lfftw3 -lm)

add_executable(evolution_v2 evolution_v2.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
    ${PROJECT_SOURCE_DIR}/modules/conjugategradient.c
    ${PROJECT_SOURCE_DIR}/modules/geometry.c
    ${PROJECT_SOURCE_DIR}/modules/hamiltonian.c
    ${PROJECT_SOURCE_DIR}/modules/observables.c
    ${PROJECT_SOURCE_DIR}/modules/powermethod.c
    ${PROJECT_SOURCE_DIR}/modules/timeevolution.c
    ${PROJECT_SOURCE_DIR}/modules/vectorutils.c
)
target_link_libraries(evolution_v2 m -lfftw3 -lm)
