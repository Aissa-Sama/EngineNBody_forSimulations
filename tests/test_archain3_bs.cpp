// tests/test_archain3_bs.cpp
//
// Tests del integrador Bulirsch-Stoer sobre AR-chain TTL.
// Referencia: Gragg (1965); Bulirsch & Stoer (1966);
//             Mikkola & Aarseth (2002), Celest. Mech. 84, 343.
//
// ── DISEÑO DE LOS TESTS ──────────────────────────────────────────────────────
//
// Test A — Figura-8, precisión alta:
//   bs_eps = 1e-10, T = 6.3259. Criterio: dE_rel < 1e-9.
//   Compara con el leapfrog ord. 2 (3.32% con eta=1e-4).
//   Ganancia mínima esperada: 6 órdenes de magnitud.
//
// Test B — Convergencia GBS (tabla de tolerancias):
//   Verifica que al pedir más precisión (eps=1e-6,1e-8,1e-10) el error
//   la sigue. Confirma el escalado exponencial propio del GBS.
//
// Test C — Pitágoras (masas 3,4,5), hasta t=20:
//   No blowup numérico. El encuentro en t≈16 debe resolverse sin NaN/Inf.
//   Con BS se espera más precisión que el leapfrog, pero el criterio es
//   simplemente estabilidad (sistema caótico → error físico O(1) es normal).
//
// Test D — Figura-8 × 10 períodos (sin deriva secular):
//   BS mantiene dE < 1e-8 uniformemente en 63.259 u.t.
//   El leapfrog simpléctico oscila ±3% — BS debe ser radicalmente más plano.
//
// Test E — Conservación del CM:
//   El GBS no puede degradar la conservación del CM respecto al leapfrog.
//   Criterio: |Δcm_pos| < 1e-14, |Δmomento| < 1e-14.
//
// Test F — Tabla comparativa BS vs leapfrog:
//   Mismo η, misma figura-8. Muestra la ganancia en precisión y el ratio
//   de coste (pasos físicos totales BS vs leapfrog).
//
// ─────────────────────────────────────────────────────────────────────────────

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include "archain3_bs_integrator.h"
#include "nbody_system.h"
#include "body.h"

// ── Helpers ───────────────────────────────────────────────────────────────────

static NBodySystem make_figure8() {
    // Figura-8 de Chenciner & Montgomery (2000).
    // CI de alta precisión (Simó 2002, 15 dígitos significativos).
    NBodySystem sys;
    sys.bodies.resize(3);
    sys.bodies[0].mass     =  1.0;
    sys.bodies[0].position = {  0.9700436926041021, -0.2430865994534989, 0.0 };
    sys.bodies[0].velocity = {  0.4662036850003824,  0.4323657300939067, 0.0 };
    sys.bodies[1].mass     =  1.0;
    sys.bodies[1].position = { -0.9700436926041021, -0.2430865994534989, 0.0 };
    sys.bodies[1].velocity = {  0.4662036850003824, -0.4323657300939067, 0.0 };
    sys.bodies[2].mass     =  1.0;
    sys.bodies[2].position = {  0.0,                 0.4861731989069978, 0.0 };
    sys.bodies[2].velocity = { -0.9324073700007648,  0.0,                0.0 };
    return sys;
}

static NBodySystem make_pythagorean() {
    // Problema de Pitágoras: masas 3, 4, 5 en triángulo rectángulo.
    // Szebehely & Peters (1967). Sistema caótico — eyección en t≈70.
    NBodySystem sys;
    sys.bodies.resize(3);
    sys.bodies[0].mass     = 3.0;
    sys.bodies[0].position = { 1.0,  3.0, 0.0 };
    sys.bodies[0].velocity = { 0.0,  0.0, 0.0 };
    sys.bodies[1].mass     = 4.0;
    sys.bodies[1].position = {-2.0, -1.0, 0.0 };
    sys.bodies[1].velocity = { 0.0,  0.0, 0.0 };
    sys.bodies[2].mass     = 5.0;
    sys.bodies[2].position = { 1.0, -1.0, 0.0 };
    sys.bodies[2].velocity = { 0.0,  0.0, 0.0 };
    return sys;
}

// ── Test A: figura-8, precisión alta ─────────────────────────────────────────

static bool test_A_figure8_high_precision() {
    std::cout << "\n── Test A: figura-8, precisión alta (bs_eps = 1e-10) ───────\n";
    std::cout << "   Criterio: dE_rel < 1e-9  (leapfrog ord.2: 3.32% con eta=1e-4)\n\n";

    NBodySystem sys = make_figure8();

    ARChain3BSIntegrator::BSParameters p;
    p.bs_eps     = 1e-10;
    p.verbose    = true;
    p.initial_ds = 1e-3;

    ARChain3BSIntegrator bs(1e-3, p);
    ARChain3State state = bs.initialize(sys, 0, 1, 2);
    const double E0 = state.energy;

    std::cout << std::fixed << std::setprecision(12)
              << "  E₀ = " << E0 << "\n\n";

    const double checkpoints[] = { 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.3259 };

    std::cout << std::left
              << std::setw(10) << "  t"
              << std::setw(16) << "|dE/E₀|"
              << "sep_min\n"
              << "  " << std::string(42, '-') << "\n";

    double max_dE = 0.0;
    for (double t_tgt : checkpoints) {
        bs.integrate_to_bs(state, t_tgt);
        const double dE  = bs.energy_error(state, E0);
        const double sep = bs.compute_sep_min(state);
        max_dE = std::max(max_dE, dE);
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  " << std::setw(8) << t_tgt;
        std::cout << std::scientific << std::setprecision(3)
                  << std::setw(16) << dE;
        std::cout << std::fixed << std::setprecision(6) << sep << "\n";
    }

    std::cout << "\n  max|dE/E₀| = " << std::scientific << max_dE;
    bool ok = (max_dE < 1e-9);
    std::cout << (ok ? "  ✅  (< 1e-9)" : "  ❌  (>= 1e-9)") << "\n";
    std::cout << "  Resultado: " << (ok ? "PASADO ✅" : "FALLADO ❌") << "\n";
    return ok;
}

// ── Test B: convergencia GBS (tabla de tolerancias) ──────────────────────────

static bool test_B_convergence_table() {
    std::cout << "\n── Test B: tabla de convergencia GBS ───────────────────────\n";
    std::cout << "   Criterio: error ≈ bs_eps (factor ≤ 10 permitido)\n\n";

    const double eps_list[] = { 1e-4, 1e-6, 1e-8, 1e-10 };
    const double tol_mult   = 10.0;

    std::cout << std::left
              << std::setw(14) << "  bs_eps"
              << std::setw(16) << "dE_rel"
              << std::setw(12) << "ratio"
              << "estado\n"
              << "  " << std::string(50, '-') << "\n";

    bool all_ok = true;
    double prev_dE = -1.0;

    for (double eps : eps_list) {
        NBodySystem sys = make_figure8();
        ARChain3BSIntegrator::BSParameters p;
        p.bs_eps     = eps;
        p.verbose    = false;
        p.initial_ds = 1e-3;

        ARChain3BSIntegrator bs(1e-3, p);
        ARChain3State state = bs.initialize(sys, 0, 1, 2);
        const double E0 = state.energy;

        bs.integrate_to_bs(state, 6.3259);
        const double dE = bs.energy_error(state, E0);

        bool ok = (dE <= eps * tol_mult || dE < 1e-14);
        std::cout << "  " << std::scientific << std::setprecision(1)
                  << std::setw(14) << eps
                  << std::setw(16) << dE;

        if (prev_dE > 0.0) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(1)
                      << (prev_dE / dE);
        } else {
            std::cout << std::setw(12) << "—";
        }
        std::cout << (ok ? " ✅" : " ❌") << "\n";
        if (!ok) all_ok = false;
        prev_dE = dE;
    }

    std::cout << "  Resultado: " << (all_ok ? "PASADO ✅" : "FALLADO ❌") << "\n";
    return all_ok;
}

// ── Test C: Pitágoras — estabilidad numérica ──────────────────────────────────

static bool test_C_pythagorean_stability() {
    std::cout << "\n── Test C: Pitágoras (masas 3,4,5) — estabilidad ──────────\n";
    std::cout << "   Criterio: no NaN/Inf hasta t=20\n\n";

    NBodySystem sys = make_pythagorean();

    ARChain3BSIntegrator::BSParameters p;
    p.bs_eps     = 1e-8;
    p.verbose    = false;
    p.initial_ds = 1e-3;

    ARChain3BSIntegrator bs(1e-3, p);
    ARChain3State state = bs.initialize(sys, 0, 1, 2);
    const double E0 = state.energy;

    std::cout << std::fixed << std::setprecision(8)
              << "  E₀ = " << E0 << "\n\n";

    const double targets[] = { 5.0, 10.0, 15.0, 20.0 };
    bool ok = true;

    for (double t : targets) {
        try {
            bs.integrate_to_bs(state, t);
            const double E_now = bs.compute_energy(state);
            const double dE    = bs.energy_error(state, E0);
            const double sep   = bs.compute_sep_min(state);
            bool step_ok = std::isfinite(E_now) && std::isfinite(sep) && sep > 0.0;
            std::cout << "  t=" << std::setw(5) << t
                      << "  sep=" << std::fixed << std::setprecision(4)
                      << std::setw(10) << sep
                      << std::scientific << std::setprecision(2)
                      << "  |dE/E|=" << dE
                      << (step_ok ? "  ✅ finito" : "  ❌ NaN/Inf") << "\n";
            if (!step_ok) ok = false;
        } catch (const std::exception& e) {
            std::cout << "  t=" << t << "  EXCEPCIÓN: " << e.what() << " ❌\n";
            ok = false;
            break;
        }
    }

    std::cout << "  Resultado: " << (ok ? "PASADO ✅" : "FALLADO ❌") << "\n";
    return ok;
}

// ── Test D: figura-8 × 10 períodos — sin deriva secular ─────────────────────

static bool test_D_long_integration() {
    std::cout << "\n── Test D: figura-8 × 10 períodos — sin deriva secular ─────\n";
    std::cout << "   Criterio: max|dE/E| < 1e-8 en 10 períodos (63.259 u.t.)\n\n";

    NBodySystem sys = make_figure8();

    ARChain3BSIntegrator::BSParameters p;
    p.bs_eps     = 1e-10;
    p.verbose    = false;
    p.initial_ds = 1e-3;

    ARChain3BSIntegrator bs(1e-3, p);
    ARChain3State state = bs.initialize(sys, 0, 1, 2);
    const double E0 = state.energy;

    double max_dE = 0.0;

    for (int i = 1; i <= 10; ++i) {
        const double t_tgt = i * 6.3259;
        bs.integrate_to_bs(state, t_tgt);
        const double dE = bs.energy_error(state, E0);
        max_dE = std::max(max_dE, dE);

        std::cout << "  período " << std::setw(2) << i
                  << "  t=" << std::fixed << std::setprecision(3)
                  << std::setw(8) << state.t_phys
                  << std::scientific << std::setprecision(3)
                  << "  |dE/E|=" << dE << "\n";
    }

    std::cout << "\n  max|dE/E| en 10 períodos = "
              << std::scientific << max_dE;
    bool ok = (max_dE < 1e-8) && std::isfinite(max_dE);
    std::cout << (ok ? "  ✅  (< 1e-8)" : "  ❌") << "\n";
    std::cout << "  Resultado: " << (ok ? "PASADO ✅" : "FALLADO ❌") << "\n";
    return ok;
}

// ── Test E: conservación del CM ───────────────────────────────────────────────

static bool test_E_cm_conservation() {
    std::cout << "\n── Test E: conservación del CM ─────────────────────────────\n";
    std::cout << "   Criterio: |Δcm_pos| < 1e-14, |Δmomento| < 1e-14\n\n";

    NBodySystem sys = make_figure8();

    ARChain3BSIntegrator::BSParameters p;
    p.bs_eps     = 1e-10;
    p.verbose    = false;
    p.initial_ds = 1e-3;

    ARChain3BSIntegrator bs(1e-3, p);
    ARChain3State state = bs.initialize(sys, 0, 1, 2);

    const Vec3   cm0  = state.cm_pos;
    const Vec3   pcm0 = state.cm_vel * state.total_mass();

    bs.integrate_to_bs(state, 6.3259);

    const Vec3   pcm_f = state.cm_vel * state.total_mass();
    const double dcm   = (state.cm_pos - cm0).norm();
    const double dp    = (pcm_f - pcm0).norm();

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  |Δcm_pos|  = " << dcm << (dcm < 1e-14 ? "  ✅" : "  ❌") << "\n";
    std::cout << "  |Δmomento| = " << dp  << (dp  < 1e-14 ? "  ✅" : "  ❌") << "\n";

    bool ok = (dcm < 1e-14 && dp < 1e-14);
    std::cout << "  Resultado: " << (ok ? "PASADO ✅" : "FALLADO ❌") << "\n";
    return ok;
}

// ── Test F: tabla comparativa BS vs leapfrog ─────────────────────────────────

static bool test_F_bs_vs_leapfrog() {
    std::cout << "\n── Test F: comparativa BS vs leapfrog TTL (figura-8) ───────\n";
    std::cout << "   Criterio: dE_BS < dE_leapfrog para cada eta\n\n";

    const double eta_list[] = { 1e-3, 1e-4 };

    std::cout << std::left
              << std::setw(8)  << "  eta"
              << std::setw(16) << "dE_leapfrog"
              << std::setw(16) << "dE_BS(1e-10)"
              << std::setw(10) << "ganancia"
              << "estado\n"
              << "  " << std::string(60, '-') << "\n";

    bool all_ok = true;

    for (double eta : eta_list) {
        // ── Leapfrog ord. 2 ──────────────────────────────────────────────
        {
            NBodySystem sys = make_figure8();
            ARChain3Integrator lf(eta);
            ARChain3State st = lf.initialize(sys, 0, 1, 2);
            const double E0  = st.energy;
            lf.integrate_to(st, 6.3259);
            const double dE_lf = std::abs(lf.compute_energy(st) - E0)
                               / std::abs(E0);

            // ── BS ────────────────────────────────────────────────────────
            NBodySystem sys2 = make_figure8();
            ARChain3BSIntegrator::BSParameters p;
            p.bs_eps  = 1e-10;
            p.verbose = false;
            ARChain3BSIntegrator bs(eta, p);
            ARChain3State st2 = bs.initialize(sys2, 0, 1, 2);
            const double E02 = st2.energy;
            bs.integrate_to_bs(st2, 6.3259);
            const double dE_bs = bs.energy_error(st2, E02);

            const double ganancia = (dE_bs > 1e-16) ? dE_lf / dE_bs : 1e16;
            bool ok = (dE_bs < dE_lf);

            std::cout << "  " << std::scientific << std::setprecision(1)
                      << std::setw(8)  << eta
                      << std::setw(16) << dE_lf
                      << std::setw(16) << dE_bs
                      << std::setw(10) << std::setprecision(1)
                      << std::fixed    << ganancia
                      << (ok ? " ✅" : " ❌") << "\n";
            if (!ok) all_ok = false;
        }
    }

    std::cout << "\n  Resultado: " << (all_ok ? "PASADO ✅" : "FALLADO ❌") << "\n";
    return all_ok;
}

// ── main ──────────────────────────────────────────────────────────────────────

int main() {
    std::cout << "════════════════════════════════════════════════════════════════\n";
    std::cout << "  TEST AR-chain BS — Extrapolación Gragg-Bulirsch-Stoer sobre TTL\n";
    std::cout << "  Referencia: Mikkola & Aarseth (2002), Celest. Mech. 84, 343\n";
    std::cout << "  Método base: Mikkola & Tanikawa (1999), MNRAS 310, 745\n";
    std::cout << "════════════════════════════════════════════════════════════════\n";

    int passed = 0, total = 0;

    auto run = [&](bool (*test)(), const char* name) {
        ++total;
        try {
            if (test()) ++passed;
        } catch (const std::exception& e) {
            std::cout << "  EXCEPCIÓN en " << name << ": " << e.what() << " ❌\n";
        }
    };

    run(test_A_figure8_high_precision,  "A: figura-8 bs_eps=1e-10");
    run(test_B_convergence_table,       "B: convergencia GBS");
    run(test_C_pythagorean_stability,   "C: Pitágoras estabilidad");
    run(test_D_long_integration,        "D: × 10 períodos");
    run(test_E_cm_conservation,         "E: CM conservado");
    run(test_F_bs_vs_leapfrog,          "F: BS vs leapfrog");

    std::cout << "\n════════════════════════════════════════════════════════════════\n";
    std::cout << "  Resultado final: " << passed << "/" << total << " tests pasados";
    std::cout << (passed == total ? "  ✅\n" : "  ❌\n");

    if (passed < total) {
        std::cout << "\n  Nota: Los tests están calibrados para el GBS sobre TTL\n";
        std::cout << "  con bs_eps = 1e-10. Referencia: Mikkola & Aarseth (2002).\n";
    }

    std::cout << "════════════════════════════════════════════════════════════════\n";
    return (passed == total) ? 0 : 1;
}