// tests/test_hierarchical.cpp
// FASE 7A: Añadido Test E — cuádruple cercano → GROUP_AR_CHAIN (N=4)
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <functional>
#include "nbody_system.h"
#include "hierarchical_integrator.h"
#include "hierarchy_builder.h"
#include "leapfrog_integrator.h"

static double total_energy(const NBodySystem& sys) {
    double T = 0, U = 0;
    const int N = sys.bodies.size();
    for (int i = 0; i < N; ++i) {
        T += 0.5 * sys.bodies[i].mass * sys.bodies[i].velocity.norm2();
        for (int j = i+1; j < N; ++j) {
            Vec3 r = sys.bodies[j].position - sys.bodies[i].position;
            U -= sys.G * sys.bodies[i].mass * sys.bodies[j].mass / r.norm();
        }
    }
    return T + U;
}

// ============================================================================
// TEST A: Binaria aislada → PAIR_KS
// ============================================================================
bool test_A() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST A: Binaria aislada (debe producir PAIR_KS)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys; sys.G = 1.0;
    sys.bodies.push_back({{-0.5, 0, 0}, {0, -0.5, 0}, 1.0});
    sys.bodies.push_back({{ 0.5, 0, 0}, {0,  0.5, 0}, 1.0});

    HierarchyBuilder::Params p; p.r_ks_threshold = 2.0;
    HierarchyBuilder hb(p);
    auto tree = hb.build(sys);

    bool has_pair_ks = false;
    std::function<void(const HierarchyNode&)> check = [&](const HierarchyNode& n) {
        if (n.type == HierarchyNode::Type::PAIR_KS) has_pair_ks = true;
        for (const auto& c : n.children) check(*c);
    };
    check(*tree);

    std::cout << "  Arbol: " << (has_pair_ks ? "PAIR_KS detectado ✓" : "PAIR_KS NO detectado ✗") << "\n";

    HierarchyBuilder::Params bp; bp.r_ks_threshold = 2.0;
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 2.0, 1e-4, bp);

    double E0 = total_energy(sys);
    std::vector<bool> used(2, false);
    for (int s = 0; s < 500; ++s) integrator.step(sys, 0.01, used);

    double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
    std::cout << "  |dE/E0| = " << dE << "\n";

    bool ok = has_pair_ks && (dE < 0.01);
    std::cout << "RESULTADO TEST A: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

// ============================================================================
// TEST B: Triple cercano → TRIPLE_CHAIN
// ============================================================================
bool test_B() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST B: Triple cercano (debe producir TRIPLE_CHAIN)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys; sys.G = 1.0;
    sys.bodies.push_back({{ 0.970, -0.243, 0}, { 0.466,  0.433, 0}, 1.0});
    sys.bodies.push_back({{-0.970, -0.243, 0}, { 0.466, -0.433, 0}, 1.0});
    sys.bodies.push_back({{ 0.000,  0.970, 0}, {-0.932,  0.000, 0}, 1.0});

    HierarchyBuilder::Params p; p.r_ks_threshold = 3.0; p.strong_coupling_eta = 5.0;
    HierarchyBuilder hb(p);
    auto tree = hb.build(sys);

    bool has_chain = false;
    std::function<void(const HierarchyNode&)> check = [&](const HierarchyNode& n) {
        if (n.type == HierarchyNode::Type::TRIPLE_CHAIN) has_chain = true;
        for (const auto& c : n.children) check(*c);
    };
    check(*tree);

    std::cout << "  Arbol: " << (has_chain ? "TRIPLE_CHAIN detectado ✓" : "TRIPLE_CHAIN NO detectado ✗") << "\n";

    HierarchyBuilder::Params bp; bp.r_ks_threshold = 3.0; bp.strong_coupling_eta = 5.0;
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 3.0, 1e-4, bp);

    double E0 = total_energy(sys);
    std::vector<bool> used(3, false);
    double max_dE = 0.0, r_max = 0.0;
    for (int s = 0; s < 300; ++s) {
        integrator.step(sys, 0.005, used);
        double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
        if (dE > max_dE) max_dE = dE;
        for (const auto& b : sys.bodies) r_max = std::max(r_max, b.position.norm());
    }

    std::cout << "  Error max energia: " << max_dE << "\n";
    std::cout << "  Radio maximo:      " << r_max << "\n";

    bool ok = has_chain && (r_max < 20.0) && (max_dE < 5.0);
    std::cout << "RESULTADO TEST B: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

// ============================================================================
// TEST C: Binaria + estrella lejana → PAIR_KS + LEAF
// ============================================================================
bool test_C() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST C: Binaria + estrella lejana (PAIR_KS + LEAF)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys; sys.G = 1.0;
    sys.bodies.push_back({{-0.5, 0, 0}, {0, -0.5, 0}, 1.0});
    sys.bodies.push_back({{ 0.5, 0, 0}, {0,  0.5, 0}, 1.0});
    sys.bodies.push_back({{50.0, 0, 0}, {0,  0.1, 0}, 0.01});

    HierarchyBuilder::Params bp; bp.r_ks_threshold = 2.0;
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 2.0, 1e-4, bp);

    double E0 = total_energy(sys);
    std::vector<bool> used(3, false);
    for (int s = 0; s < 300; ++s) integrator.step(sys, 0.01, used);

    double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
    double r_field = sys.bodies[2].position.norm();

    std::cout << "  |dE/E0| = " << dE << "\n";
    std::cout << "  Estrella lejana en r = " << r_field << " (debe ser ~50)\n";

    bool ok = (dE < 0.05) && (std::abs(r_field - 50.0) < 5.0);
    std::cout << "RESULTADO TEST C: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

// ============================================================================
// TEST D: Triple + estrella lejana → TRIPLE_CHAIN + LEAF
// ============================================================================
bool test_D() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST D: Triple + estrella lejana (TRIPLE_CHAIN + LEAF)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys; sys.G = 1.0;
    sys.bodies.push_back({{ 0.970, -0.243, 0}, { 0.466,  0.433, 0}, 1.0});
    sys.bodies.push_back({{-0.970, -0.243, 0}, { 0.466, -0.433, 0}, 1.0});
    sys.bodies.push_back({{ 0.000,  0.970, 0}, {-0.932,  0.000, 0}, 1.0});
    sys.bodies.push_back({{100.0,   0.000, 0}, { 0.000,  0.050, 0}, 0.01});

    HierarchyBuilder::Params bp; bp.r_ks_threshold = 3.0; bp.strong_coupling_eta = 5.0;
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 3.0, 1e-4, bp);

    double E0 = total_energy(sys);
    std::vector<bool> used(4, false);
    double max_dE = 0.0;
    for (int s = 0; s < 200; ++s) {
        integrator.step(sys, 0.005, used);
        double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
        if (dE > max_dE) max_dE = dE;
    }

    double r_field = sys.bodies[3].position.norm();
    std::cout << "  Error max energia: " << max_dE << "\n";
    std::cout << "  Estrella lejana en r = " << r_field << " (debe ser ~100)\n";

    bool ok = (max_dE < 5.0) && (std::abs(r_field - 100.0) < 10.0);
    std::cout << "RESULTADO TEST D: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

// ============================================================================
// TEST E: Cuádruple cercano → GROUP_AR_CHAIN (N=4)   [FASE 7A]
//
// Sistema: cuadrado virializado de 4 masas iguales en (±1,0) y (0,±1).
// Separaciones: adyacentes √2 ≈ 1.41, diagonales 2.0 — todas < r_ks=3.0.
// Energía colectiva < 0 → grupo ligado → debe detectarse como GROUP_AR_CHAIN.
//
// Criterios:
//   - Árbol produce GROUP_AR_CHAIN con 4 índices     (detección correcta)
//   - |ΔE/E0| < 5% en 300 pasos dt=0.005            (integración estable)
//   - r_max < 20                                      (sin explosión)
// ============================================================================
bool test_E() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST E: Cuadruple cercano → GROUP_AR_CHAIN (N=4) [Fase 7A]\n";
    std::cout << std::string(60,'=') << "\n";

    // Cuadrado virializado: virial ratio = 1.000 exacto
    // T=1.9142, U=-3.8284, E=-1.9142
    constexpr double V = 0.9783183434785159;
    NBodySystem sys; sys.G = 1.0;
    sys.bodies.push_back({{ 1.0, 0.0, 0}, {0.0,  V, 0}, 1.0});
    sys.bodies.push_back({{ 0.0, 1.0, 0}, { -V, 0.0, 0}, 1.0});
    sys.bodies.push_back({{-1.0, 0.0, 0}, {0.0, -V, 0}, 1.0});
    sys.bodies.push_back({{ 0.0,-1.0, 0}, {  V, 0.0, 0}, 1.0});

    // r_ks_threshold=3.0: sep_max=2.0 < 3.0 → todos los pares dentro del umbral
    // ar_chain_threshold=0.5: sep_min=√2≈1.41 > 0.5, pero N=4 → GROUP_AR_CHAIN igualmente
    HierarchyBuilder::Params bp;
    bp.r_ks_threshold    = 3.0;
    bp.ar_chain_threshold = 0.5;
    bp.ar_chain_eta      = 1e-3;

    // Verificar detección en el árbol
    {
        HierarchyBuilder hb(bp);
        auto tree = hb.build(sys);

        bool has_group4 = false;
        std::function<void(const HierarchyNode&)> check =
            [&](const HierarchyNode& n) {
                if (n.type == HierarchyNode::Type::GROUP_AR_CHAIN
                    && (int)n.body_indices.size() == 4)
                    has_group4 = true;
                for (const auto& c : n.children) check(*c);
            };
        check(*tree);

        std::cout << "  Arbol: "
                  << (has_group4
                      ? "GROUP_AR_CHAIN(N=4) detectado ✓"
                      : "GROUP_AR_CHAIN(N=4) NO detectado ✗") << "\n";

        if (!has_group4) {
            std::cout << "RESULTADO TEST E: ✗ FALLADO\n";
            return false;
        }
    }

    // Integrar con HierarchicalIntegrator
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 3.0, 1e-4, bp);

    const double E0 = total_energy(sys);
    std::cout << "  E0 = " << std::setprecision(8) << E0 << "\n";

    std::vector<bool> used(4, false);
    double max_dE = 0.0, r_max = 0.0;

    for (int s = 0; s < 300; ++s) {
        integrator.step(sys, 0.005, used);
        double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
        if (dE > max_dE) max_dE = dE;
        for (const auto& b : sys.bodies)
            r_max = std::max(r_max, b.position.norm());
    }

    std::cout << "  Error max energia: " << max_dE << "\n";
    std::cout << "  Radio maximo:      " << r_max  << "\n";
    std::cout << "  ¿Sin explosion?    " << (r_max < 20.0 ? "✓" : "✗") << "\n";
    std::cout << "  ¿Energia OK?       " << (max_dE < 0.05 ? "✓" : "✗")
              << "  (criterio <5%)\n";

    bool ok = (r_max < 20.0) && (max_dE < 0.05);
    std::cout << "RESULTADO TEST E: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

// ============================================================================
// MAIN
// ============================================================================
int main() {
    std::cout << std::string(72,'#') << "\n";
    std::cout << "TESTS: HierarchicalIntegrator (Ruta B)\n";
    std::cout << std::string(72,'#') << "\n";

    bool a = test_A();
    bool b = test_B();
    bool c = test_C();
    bool d = test_D();
    bool e = test_E();

    std::cout << "\n" << std::string(72,'=') << "\n";
    std::cout << "RESUMEN:\n";
    std::cout << "  A (binaria → PAIR_KS):              " << (a ? "✓ PASADO" : "✗ FALLADO") << "\n";
    std::cout << "  B (triple → TRIPLE_CHAIN):          " << (b ? "✓ PASADO" : "✗ FALLADO") << "\n";
    std::cout << "  C (binaria + campo):                " << (c ? "✓ PASADO" : "✗ FALLADO") << "\n";
    std::cout << "  D (triple + campo):                 " << (d ? "✓ PASADO" : "✗ FALLADO") << "\n";
    std::cout << "  E (cuadruple → GROUP_AR_CHAIN N=4): " << (e ? "✓ PASADO" : "✗ FALLADO") << "\n";
    std::cout << std::string(72,'=') << "\n";

    if (a && b && c && d && e) {
        std::cout << "✓ Todos los tests pasaron. Fase 7A completa.\n";
        return 0;
    }
    std::cout << "✗ Algunos tests fallaron.\n";
    return 1;
}
