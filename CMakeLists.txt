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
# Garantizan resultados bit-a-bit idénticos entre runs con el mismo
# compilador y hardware. Ver REPRODUCIBILITY_04_03.md para justificación.
#
# NUNCA usar:
#   MSVC:  /fp:fast         — permite reordenar FP, rompe determinismo
#   GCC:   -ffast-math      — activa todas las optimizaciones FP inseguras
#          -ffinite-math-only — asume sin NaN/Inf → resultados incorrectos
#                              en encuentros cercanos (singularidad 1/r)
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
# 1. LIBRERÍA ESTÁTICA CON TODO EL CÓDIGO (excepto main)
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

# NOTA: el GLOB captura todos los .cpp de cada subdirectorio,
# incluyendo archain3_bs_integrator.cpp y chain3_bs_integrator.cpp.
file(GLOB LIB_SOURCES
    "analysis/*.cpp"
    "integrators/*.cpp"
    "io/*.cpp"
    "regularization/ks/*.cpp"
    "regularization/chain/*.cpp"
    "regularization/hierarchy/*.cpp"
    "initial_conditions.cpp"
    "core/units.cpp"
    "core/nbody_system.cpp"
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
foreach(_t
    test_ks
    test_chain_ks
    test_chain3_state
    test_chain3_integration
    test_chain3_diagnostic
    test_hybrid_chain3
    test_hierarchical
    test_archain3
    test_archain3_bs
    test_chain3_bs
    test_bs_vs_rk4
    test_validation_benchmarks
    test_archain_n
    test_archain_n_pn
)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/${_t}.cpp)
        add_executable(${_t} tests/${_t}.cpp)
        target_link_libraries(${_t} PRIVATE nbody_lib)
    endif()
endforeach()

# ============================================================================
# 4. FLAGS DE ADVERTENCIAS (separadas de los flags de determinismo)
# ============================================================================
set(ALL_TEST_TARGETS
    test_ks test_chain_ks test_chain3_state test_chain3_integration
    test_chain3_diagnostic test_hybrid_chain3 test_hierarchical
    test_archain3 test_archain3_bs test_chain3_bs test_bs_vs_rk4
    test_validation_benchmarks test_archain_n test_archain_n_pn
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

# ============================================================================
# COMPILAR:  cmake --build .
# EJECUTAR:  .\test_archain3.exe  /  .\test_hierarchical.exe  etc.
# ============================================================================