include_directories(${FFTW_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/modules)
link_directories(/usr/local/Cellar/fftw/3.3.10/lib)



add_executable(cfg cfg.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
)
target_link_libraries(cfg m)
add_test(NAME cfg COMMAND cfg)



add_executable(powermethod powermethod.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
    ${PROJECT_SOURCE_DIR}/modules/powermethod.c
    ${PROJECT_SOURCE_DIR}/modules/vectorutils.c
)
target_link_libraries(powermethod m)
add_test(NAME powermethod COMMAND powermethod)



add_executable(conjugategradient conjugategradient.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
    ${PROJECT_SOURCE_DIR}/modules/conjugategradient.c
    ${PROJECT_SOURCE_DIR}/modules/vectorutils.c
)
target_link_libraries(conjugategradient m)
add_test(NAME conjugategradient COMMAND conjugategradient)



add_executable(geometry geometry.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
    ${PROJECT_SOURCE_DIR}/modules/geometry.c
)
target_link_libraries(geometry m)
add_test(NAME geometry COMMAND geometry)



add_executable(hamiltonian hamiltonian.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
    ${PROJECT_SOURCE_DIR}/modules/hamiltonian.c
    ${PROJECT_SOURCE_DIR}/modules/vectorutils.c
    ${PROJECT_SOURCE_DIR}/modules/geometry.c
)
target_link_libraries(hamiltonian m)
add_test(NAME hamiltonian COMMAND hamiltonian)



add_executable(timeevolution timeevolution.c
    ${PROJECT_SOURCE_DIR}/modules/cfg.c
    ${PROJECT_SOURCE_DIR}/modules/conjugategradient.c
    ${PROJECT_SOURCE_DIR}/modules/hamiltonian.c
    ${PROJECT_SOURCE_DIR}/modules/vectorutils.c
    ${PROJECT_SOURCE_DIR}/modules/geometry.c
    ${PROJECT_SOURCE_DIR}/modules/timeevolution.c
    ${PROJECT_SOURCE_DIR}/modules/observables.c
)
target_link_libraries(timeevolution m -lfftw3 -lm)
add_test(NAME timeevolution COMMAND timeevolution)
