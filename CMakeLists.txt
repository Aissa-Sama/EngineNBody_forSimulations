cmake_minimum_required(VERSION 3.15)
project(nbody_sim)

# ============================================================================
# CONFIGURACIÓN BÁSICA DEL LENGUAJE
# ============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================================
# 0. FLAGS DE DETERMINISMO IEEE 754
# ============================================================================
if(MSVC)
    add_compile_options(/fp:precise)
    add_compile_options(/fp:except-)
else()
    add_compile_options(-ffp-contract=off)
    add_compile_options(-fno-unsafe-math-optimizations)
    add_compile_options(-fno-associative-math)
endif()

# ============================================================================
# 1. LIBRERÍA ESTÁTICA
# ============================================================================
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}/analysis
    ${CMAKE_CURRENT_SOURCE_DIR}/integrators
    ${CMAKE_CURRENT_SOURCE_DIR}/io
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/ks
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/chain
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/hierarchy
)

file(GLOB LIB_SOURCES
    "analysis/*.cpp"
    "integrators/*.cpp"
    "io/*.cpp"
    "regularization/ks/*.cpp"
    "regularization/chain/*.cpp"
    "initial_conditions.cpp"
    "core/units.cpp"
    "core/nbody_system.cpp"
    "regularization/hierarchy/*.cpp"
    "regularization/chain/chain3_bs_integrator.cpp"
)

add_library(nbody_lib STATIC ${LIB_SOURCES})

# ============================================================================
# 2. EJECUTABLE PRINCIPAL
# ============================================================================
add_executable(nbody_sim main.cpp)
target_link_libraries(nbody_sim PRIVATE nbody_lib)

# ============================================================================
# 3. EJECUTABLES DE TEST
# ============================================================================

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ks.cpp)
    add_executable(test_ks tests/test_ks.cpp)
    target_link_libraries(test_ks PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain_ks.cpp)
    add_executable(test_chain_ks tests/test_chain_ks.cpp)
    target_link_libraries(test_chain_ks PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain3_state.cpp)
    add_executable(test_chain3_state tests/test_chain3_state.cpp)
    target_link_libraries(test_chain3_state PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain3_integration.cpp)
    add_executable(test_chain3_integration tests/test_chain3_integration.cpp)
    target_link_libraries(test_chain3_integration PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain3_diagnostic.cpp)
    add_executable(test_chain3_diagnostic tests/test_chain3_diagnostic.cpp)
    target_link_libraries(test_chain3_diagnostic PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hybrid_chain3.cpp)
    add_executable(test_hybrid_chain3 tests/test_hybrid_chain3.cpp)
    target_link_libraries(test_hybrid_chain3 PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hierarchical.cpp)
    add_executable(test_hierarchical tests/test_hierarchical.cpp)
    target_link_libraries(test_hierarchical PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_archain3.cpp)
    add_executable(test_archain3 tests/test_archain3.cpp)
    target_link_libraries(test_archain3 PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain3_bs.cpp)
    add_executable(test_chain3_bs tests/test_chain3_bs.cpp)
    target_link_libraries(test_chain3_bs PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_bs_vs_rk4.cpp)
    add_executable(test_bs_vs_rk4 tests/test_bs_vs_rk4.cpp)
    target_link_libraries(test_bs_vs_rk4 PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hierarchical_pn.cpp)
    add_executable(test_hierarchical_pn tests/test_hierarchical_pn.cpp)
    target_link_libraries(test_hierarchical_pn PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_archain_n_ks.cpp)
    add_executable(test_archain_n_ks tests/test_archain_n_ks.cpp)
    target_link_libraries(test_archain_n_ks PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_gw_observables.cpp)
    add_executable(test_gw_observables tests/test_gw_observables.cpp)
    target_link_libraries(test_gw_observables PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_physical_benchmarks.cpp)
    add_executable(test_physical_benchmarks tests/test_physical_benchmarks.cpp)
    target_link_libraries(test_physical_benchmarks PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_archain3_bs.cpp)
    add_executable(test_archain3_bs tests/test_archain3_bs.cpp)
    target_link_libraries(test_archain3_bs PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_validation_benchmarks.cpp)
    add_executable(test_validation_benchmarks tests/test_validation_benchmarks.cpp)
    target_link_libraries(test_validation_benchmarks PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_archain_n.cpp)
    add_executable(test_archain_n tests/test_archain_n.cpp)
    target_link_libraries(test_archain_n PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_archain_n_pn.cpp)
    add_executable(test_archain_n_pn tests/test_archain_n_pn.cpp)
    target_link_libraries(test_archain_n_pn PRIVATE nbody_lib)
endif()

# ============================================================================
# 4. PYTHON BINDING — nbody_core.pyd   (Fase 6D)
#    Usa pybind11 local: extern/pybind11 (git clone)
# ============================================================================
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/CMakeLists.txt
   AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python/nbody_bindings.cpp)

    add_subdirectory(extern/pybind11)

    pybind11_add_module(nbody_core python/nbody_bindings.cpp)

    target_include_directories(nbody_core PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/core
        ${CMAKE_CURRENT_SOURCE_DIR}/analysis
        ${CMAKE_CURRENT_SOURCE_DIR}/integrators
        ${CMAKE_CURRENT_SOURCE_DIR}/io
        ${CMAKE_CURRENT_SOURCE_DIR}/regularization
        ${CMAKE_CURRENT_SOURCE_DIR}/regularization/ks
        ${CMAKE_CURRENT_SOURCE_DIR}/regularization/chain
        ${CMAKE_CURRENT_SOURCE_DIR}/regularization/hierarchy
    )

    target_link_libraries(nbody_core PRIVATE nbody_lib)

    if(MSVC)
        target_compile_options(nbody_core PRIVATE /W4 /utf-8 /wd4100 /wd4505)
    endif()

    message(STATUS "Target nbody_core.pyd configurado (pybind11 local)")

else()
    message(STATUS "extern/pybind11 o python/nbody_bindings.cpp no encontrado — binding omitido")
endif()

# ============================================================================
# 5. OPCIONES DE ADVERTENCIAS
# ============================================================================
set(ALL_TEST_TARGETS
    test_ks test_chain_ks test_chain3_state test_chain3_integration
    test_chain3_diagnostic test_hybrid_chain3 test_hierarchical
    test_archain3 test_chain3_bs test_bs_vs_rk4
    test_archain3_bs test_validation_benchmarks test_archain_n
    test_archain_n_pn test_hierarchical_pn
    test_archain_n_ks test_gw_observables
    test_physical_benchmarks
)

if(MSVC)
    target_compile_options(nbody_lib PRIVATE /W4)
    target_compile_options(nbody_sim PRIVATE /W4 /utf-8)
    foreach(_t ${ALL_TEST_TARGETS})
        if(TARGET ${_t})
            target_compile_options(${_t} PRIVATE /W4 /utf-8)
        endif()
    endforeach()
else()
    target_compile_options(nbody_lib PRIVATE -Wall -Wextra)
    target_compile_options(nbody_sim PRIVATE -Wall -Wextra)
    foreach(_t ${ALL_TEST_TARGETS})
        if(TARGET ${_t})
            target_compile_options(${_t} PRIVATE -Wall -Wextra)
        endif()
    endforeach()
endif()
