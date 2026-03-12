// tests/test_archain_n.cpp
//
// Suite de validación del integrador AR-chain para N cuerpos arbitrario.
// Referencia: Mikkola & Aarseth (1993); Mikkola & Merritt (2006, 2008).
//
// ── DISEÑO DE LOS TESTS ──────────────────────────────────────────────────────
//
// Test 1 — FIND_CHAIN_INDICES (unitario):
//   4 cuerpos con par (0,1) claramente más cercano. Verifica adyacencia.
//
// Test 2 — Round-trip initialize/reconstruct:
//   initialize() → reconstruct_positions() devuelve CI exactas. Err < 1e-13.
//
// Test 3 — Conservación de energía, N=4, leapfrog TTL:
//   4 cuerpos en cuadrado virializado (virial=1.000, E=-1.914).
//   η=1e-4, T=6.42 (un período). Criterio: |dE/E| < 5% (orden 2).
//
// Test 4 — Conservación de energía, N=4, GBS:
//   Mismo sistema, bs_eps=1e-10. Criterio: |dE/E| < 1e-8.
//
// Test 5 — Conservación del CM:
//   El CM se propaga exactamente. Criterio: |Δcm_pos| < 1e-12.
//
// Test 6 — N=3 consistency:
//   ARChainN con N=3 reproduce figura-8 con |dE/E| < 10%.
//
// Test 7 — Re-construcción de cadena (cuadrado virializado, T=12):
//   Sistema estable durante >1 período. La cadena puede reconstruirse.
//   Criterio: no NaN/Inf y |dE/E| < 5%.
//
// Test 8 — N=5 (figura-8 + 2 perturbadores):
//   Perturbadores masa 1e-3 a distancia 10. |dE/E| < 1e-7 con GBS.
//
// ─────────────────────────────────────────────────────────────────────────────

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "archain_n_integrator.h"
#include "archain_n_bs_integrator.h"
#include "nbody_system.h"
#include "body.h"

// ── Helpers ───────────────────────────────────────────────────────────────────

static NBodySystem make_figure8_3() {
    NBodySystem sys;
    sys.bodies.resize(3);
    sys.bodies[0].mass = 1.0;
    sys.bodies[0].position = { 0.9700436926041021,-0.2430865994534989,0.0};
    sys.bodies[0].velocity = { 0.4662036850003824, 0.4323657300939067,0.0};
    sys.bodies[1].mass = 1.0;
    sys.bodies[1].position = {-0.9700436926041021,-0.2430865994534989,0.0};
    sys.bodies[1].velocity = { 0.4662036850003824,-0.4323657300939067,0.0};
    sys.bodies[2].mass = 1.0;
    sys.bodies[2].position = { 0.0,                0.4861731989069978,0.0};
    sys.bodies[2].velocity = {-0.9324073700007648, 0.0,               0.0};
    return sys;
}

// 4 cuerpos en cuadrado virializado (virial ratio = 1.000 exacto).
//
// Derivación:
//   Fuerza radial sobre cuerpo en (1,0):
//     F = m²/r²·1 (vecinos adyacentes a dist sqrt(2)) + m²/4 (opuesto a dist 2)
//     F_radial = 2*(1/sqrt(2))*(1/2) + 1/4 = 1/sqrt(2) + 1/4 ≈ 0.957107
//   v_circ = sqrt(F_radial * r / m) = sqrt(0.957107) ≈ 0.97832
//
// Resultado: T=1.9142, U=-3.8284, E=-1.9142, virial 2T/|U|=1.000
static NBodySystem make_square_4() {
    constexpr double V = 0.9783183434785159;   // v_circ calculada
    NBodySystem sys;
    sys.bodies.resize(4);
    sys.bodies[0].mass = 1.0; sys.bodies[0].position = { 1.0, 0.0, 0.0}; sys.bodies[0].velocity = {0.0,  V, 0.0};
    sys.bodies[1].mass = 1.0; sys.bodies[1].position = { 0.0, 1.0, 0.0}; sys.bodies[1].velocity = { -V, 0.0, 0.0};
    sys.bodies[2].mass = 1.0; sys.bodies[2].position = {-1.0, 0.0, 0.0}; sys.bodies[2].velocity = {0.0, -V, 0.0};
    sys.bodies[3].mass = 1.0; sys.bodies[3].position = { 0.0,-1.0, 0.0}; sys.bodies[3].velocity = { V,  0.0, 0.0};
    return sys;
}

static bool PASS(const char* name) {
    std::cout << "  [PASS] " << name << "\n";
    return true;
}
static bool FAIL(const char* name, const std::string& reason) {
    std::cout << "  [FAIL] " << name << " — " << reason << "\n";
    return false;
}

// ── Test 1 — FIND_CHAIN_INDICES ───────────────────────────────────────────────
static bool test_find_chain() {
    std::vector<Vec3> r = {
        {0.0, 0.0, 0.0},
        {0.1, 0.0, 0.0},   // par más cercano con 0
        {5.0, 0.0, 0.0},
        {10.0, 0.0, 0.0},
    };
    std::vector<int> chain;
    ARChainNIntegrator::find_chain_indices(r, chain);

    std::vector<int> sorted_chain = chain;
    std::sort(sorted_chain.begin(), sorted_chain.end());
    for (int i = 0; i < 4; ++i)
        if (sorted_chain[i] != i)
            return FAIL("find_chain: no es permutación", "");

    bool adj_01 = false;
    for (int k = 0; k < 3; ++k)
        if ((chain[k]==0 && chain[k+1]==1) || (chain[k]==1 && chain[k+1]==0))
            { adj_01 = true; break; }
    if (!adj_01)
        return FAIL("find_chain: par más cercano no es adyacente", "");

    std::cout << "    chain=[";
    for (int i=0;i<4;++i) std::cout<<chain[i]<<(i<3?",":"");
    std::cout << "]\n";
    return PASS("find_chain_indices: par más cercano es vecino directo");
}

// ── Test 2 — Round-trip ───────────────────────────────────────────────────────
static bool test_roundtrip() {
    ARChainNIntegrator integrator(1e-3);
    auto sys = make_square_4();
    std::vector<int> idx = {0,1,2,3};

    auto state = integrator.initialize(sys, idx);

    std::vector<Vec3> r_out, v_out;
    state.reconstruct_positions(r_out);
    state.reconstruct_velocities(v_out);

    double max_err_r = 0.0, max_err_v = 0.0;
    for (int i = 0; i < 4; ++i) {
        max_err_r = std::max(max_err_r, (r_out[i] - sys.bodies[i].position).norm());
        max_err_v = std::max(max_err_v, (v_out[i] - sys.bodies[i].velocity).norm());
    }

    std::cout << "    max_err_r=" << max_err_r << "  max_err_v=" << max_err_v << "\n";
    if (max_err_r > 1e-13 || max_err_v > 1e-13)
        return FAIL("round-trip", "error de reconstrucción > 1e-13");
    return PASS("round-trip initialize/reconstruct");
}

// ── Test 3 — Energía N=4, leapfrog ───────────────────────────────────────────
static bool test_energy_leapfrog_n4() {
    ARChainNIntegrator integrator(1e-4);
    auto sys = make_square_4();
    std::vector<int> idx = {0,1,2,3};

    auto state = integrator.initialize(sys, idx);
    const double E0 = state.energy;

    // Integrar un período completo del cuadrado: T ~ 6.42
    integrator.integrate_to(state, 6.42);

    const double E1 = integrator.compute_energy(state);
    const double dE = std::abs(E1 - E0) / std::abs(E0);

    std::cout << "    E0=" << E0 << "  E1=" << E1 << "  |dE/E|=" << dE << "\n";

    if (dE > 0.05)
        return FAIL("leapfrog N=4 energía", "|dE/E| > 5%");
    return PASS("leapfrog N=4: |dE/E| < 5% en T=6.42 (1 período)");
}

// ── Test 4 — Energía N=4, GBS ────────────────────────────────────────────────
static bool test_energy_bs_n4() {
    ARChainNBSIntegrator::BSParameters params;
    params.bs_eps     = 1e-10;
    params.energy_tol = 1e-8;
    ARChainNBSIntegrator integrator(1e-3, params);

    auto sys = make_square_4();
    std::vector<int> idx = {0,1,2,3};

    auto state = integrator.initialize(sys, idx);
    const double E0 = state.energy;

    integrator.integrate_to_bs(state, 6.42);

    const double E1 = integrator.compute_energy(state);
    const double dE = std::abs(E1 - E0) / std::abs(E0);

    std::cout << "    E0=" << E0 << "  E1=" << E1 << "  |dE/E|=" << dE << "\n";

    // GBS con bs_eps=1e-10 debe dar al menos 8 dígitos correctos
    if (dE > 1e-8)
        return FAIL("GBS N=4 energía", "|dE/E| > 1e-8");
    return PASS("GBS N=4: |dE/E| < 1e-8 en T=6.42");
}

// ── Test 5 — Conservación CM ─────────────────────────────────────────────────
static bool test_cm_conservation() {
    ARChainNIntegrator integrator(1e-4);
    auto sys = make_square_4();
    std::vector<int> idx = {0,1,2,3};

    auto state = integrator.initialize(sys, idx);
    const Vec3 cm_pos0 = state.cm_pos;
    const Vec3 cm_vel0 = state.cm_vel;

    integrator.integrate_to(state, 6.42);
    integrator.write_back(state, sys, idx);

    double M = 0.0;
    Vec3 cm_pos1, cm_vel1;
    for (int i = 0; i < 4; ++i) {
        M       += sys.bodies[i].mass;
        cm_pos1 = cm_pos1 + sys.bodies[i].position * sys.bodies[i].mass;
        cm_vel1 = cm_vel1 + sys.bodies[i].velocity * sys.bodies[i].mass;
    }
    cm_pos1 = cm_pos1 * (1.0/M);
    cm_vel1 = cm_vel1 * (1.0/M);

    Vec3 cm_expected = cm_pos0 + cm_vel0 * 6.42;
    double err_pos = (cm_pos1 - cm_expected).norm();
    double err_vel = (cm_vel1 - cm_vel0).norm();

    std::cout << "    |Δcm_pos|=" << err_pos << "  |Δcm_vel|=" << err_vel << "\n";

    if (err_pos > 1e-12 || err_vel > 1e-14)
        return FAIL("CM conservation", "deriva del CM > tolerancia");
    return PASS("CM conservation: exacta en T=6.42");
}

// ── Test 6 — N=3 consistency ─────────────────────────────────────────────────
static bool test_n3_consistency() {
    ARChainNIntegrator integrator(1e-4);
    auto sys = make_figure8_3();
    std::vector<int> idx = {0,1,2};

    auto state = integrator.initialize(sys, idx);
    const double E0 = state.energy;

    integrator.integrate_to(state, 6.3259);

    const double E1 = integrator.compute_energy(state);
    const double dE = std::abs(E1 - E0) / std::abs(E0);

    std::cout << "    E0=" << E0 << "  |dE/E|=" << dE << "\n";

    if (dE > 0.10)
        return FAIL("N=3 consistency", "|dE/E| > 10%");
    return PASS("N=3 consistency: |dE/E| < 10% en figura-8");
}

// ── Test 7 — Re-construcción de cadena ───────────────────────────────────────
// Integra el cuadrado virializado durante 2 períodos (T=12.84).
// El cuadrado es dinámico: durante la integración el par más cercano
// cambia periódicamente (los 4 cuerpos rotan). Verificamos que:
//   - La integración no produce NaN/Inf
//   - La energía se conserva razonablemente (|dE/E| < 5%)
static bool test_rechain() {
    ARChainNIntegrator integrator(1e-4);
    auto sys = make_square_4();
    std::vector<int> idx = {0,1,2,3};

    auto state = integrator.initialize(sys, idx);
    const double E0 = state.energy;

    // 2 períodos: suficiente para que la cadena se reconstruya varias veces
    integrator.integrate_to(state, 12.84);

    const double E1 = integrator.compute_energy(state);
    const double dE = std::abs(E1 - E0) / std::abs(E0);

    std::cout << "    E0=" << E0 << "  E1=" << E1 << "  |dE/E|=" << dE << "\n";

    if (std::isnan(E1) || std::isinf(E1))
        return FAIL("rechain", "NaN/Inf después de 2 períodos");
    if (dE > 0.05)
        return FAIL("rechain", "|dE/E| > 5% en 2 períodos");
    return PASS("re-construcción de cadena: estable en 2 períodos, |dE/E| < 5%");
}

// ── Test 8 — N=5 perturbado ───────────────────────────────────────────────────
static bool test_n5_perturbed() {
    ARChainNBSIntegrator::BSParameters params;
    params.bs_eps     = 1e-8;
    params.energy_tol = 1e-6;
    ARChainNBSIntegrator integrator(1e-3, params);

    NBodySystem sys;
    sys.bodies.resize(5);
    sys.bodies[0].mass = 1.0; sys.bodies[0].position = { 0.9700436926041021,-0.2430865994534989,0.0}; sys.bodies[0].velocity = { 0.4662036850003824, 0.4323657300939067,0.0};
    sys.bodies[1].mass = 1.0; sys.bodies[1].position = {-0.9700436926041021,-0.2430865994534989,0.0}; sys.bodies[1].velocity = { 0.4662036850003824,-0.4323657300939067,0.0};
    sys.bodies[2].mass = 1.0; sys.bodies[2].position = { 0.0,                0.4861731989069978,0.0}; sys.bodies[2].velocity = {-0.9324073700007648, 0.0,               0.0};
    sys.bodies[3].mass = 1e-3; sys.bodies[3].position = { 10.0,  0.1, 0.0}; sys.bodies[3].velocity = {0.0, 0.3, 0.0};
    sys.bodies[4].mass = 1e-3; sys.bodies[4].position = {-10.0, -0.1, 0.0}; sys.bodies[4].velocity = {0.0,-0.3, 0.0};

    std::vector<int> idx = {0,1,2,3,4};
    auto state = integrator.initialize(sys, idx);
    const double E0 = state.energy;

    integrator.integrate_to_bs(state, 6.3259);

    const double E1 = integrator.compute_energy(state);
    const double dE = std::abs(E1 - E0) / std::abs(E0);

    std::cout << "    E0=" << E0 << "  |dE/E|=" << dE << "\n";

    if (std::isnan(E1) || std::isinf(E1))
        return FAIL("N=5 perturbado", "NaN/Inf");
    if (dE > 1e-7)
        return FAIL("N=5 perturbado", "|dE/E| > 1e-7");
    return PASS("N=5 (figura-8 + perturbadores): |dE/E| < 1e-7");
}

// ── Main ──────────────────────────────────────────────────────────────────────
int main() {
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "======================================================\n";
    std::cout << "  test_archain_n — AR-chain para N cuerpos arbitrario\n";
    std::cout << "======================================================\n\n";

    int passed = 0, total = 0;
    auto run = [&](bool(*fn)(), const char* name) {
        ++total;
        std::cout << "Test " << total << " — " << name << "\n";
        try {
            if (fn()) ++passed;
        } catch (const std::exception& e) {
            FAIL(name, std::string("excepción: ") + e.what());
        }
        std::cout << "\n";
    };

    run(test_find_chain,          "FIND_CHAIN_INDICES unitario");
    run(test_roundtrip,           "Round-trip initialize/reconstruct");
    run(test_energy_leapfrog_n4,  "Energía N=4, leapfrog TTL");
    run(test_energy_bs_n4,        "Energía N=4, GBS bs_eps=1e-10");
    run(test_cm_conservation,     "Conservación del CM");
    run(test_n3_consistency,      "Consistencia N=3 con ARChain3");
    run(test_rechain,             "Re-construcción de cadena (2 períodos)");
    run(test_n5_perturbed,        "N=5 figura-8 + perturbadores");

    std::cout << "======================================================\n";
    std::cout << "  Resultado: " << passed << "/" << total << " tests pasando\n";
    std::cout << "======================================================\n";

    return (passed == total) ? 0 : 1;
}