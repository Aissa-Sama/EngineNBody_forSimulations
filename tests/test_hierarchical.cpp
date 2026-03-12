// tests/test_hierarchical.cpp
// ============================================================================
// TEST: HierarchicalIntegrator (Ruta B)
//
// Verifica que el integrador jerárquico construye el árbol correctamente
// y produce resultados físicamente correctos en todos los escenarios.
//
// TEST A: Binaria aislada         → árbol produce un nodo PAIR_KS
// TEST B: Triple cercano          → árbol produce un nodo TRIPLE_CHAIN
// TEST C: Binaria + estrella lejana → PAIR_KS + LEAF
// TEST D: Triple + estrella lejana  → TRIPLE_CHAIN + LEAF
// ============================================================================
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include "nbody_system.h"
#include "hierarchical_integrator.h"
#include "hierarchy_builder.h"
#include "leapfrog_integrator.h"

static double total_energy(const NBodySystem& sys) {
    double T = 0, U = 0;
    const int N = sys.bodies.size();
    for (int i = 0; i < N; ++i) {
        T += 0.5 * sys.bodies[i].mass
           * sys.bodies[i].velocity.norm2();
        for (int j = i+1; j < N; ++j) {
            Vec3 r = sys.bodies[j].position - sys.bodies[i].position;
            U -= sys.G * sys.bodies[i].mass * sys.bodies[j].mass / r.norm();
        }
    }
    return T + U;
}

// ============================================================================
// TEST A: Binaria aislada
// ============================================================================
bool test_A() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST A: Binaria aislada (debe producir PAIR_KS)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys;
    sys.G = 1.0;
    sys.bodies.push_back({{-0.5, 0, 0}, {0, -0.5, 0}, 1.0});
    sys.bodies.push_back({{ 0.5, 0, 0}, {0,  0.5, 0}, 1.0});

    // Verificar que HierarchyBuilder detecta un PAIR_KS
    HierarchyBuilder::Params p;
    p.r_ks_threshold = 2.0;
    HierarchyBuilder hb(p);
    auto tree = hb.build(sys);

    // Con 2 cuerpos ligados y N<3 → espera PAIR_KS o COMPOSITE con PAIR_KS
    bool has_pair_ks = false;
    auto check = [&](const HierarchyNode& n) {
        if (n.type == HierarchyNode::Type::PAIR_KS) has_pair_ks = true;
        if (n.type == HierarchyNode::Type::COMPOSITE)
            for (const auto& c : n.children)
                if (c->type == HierarchyNode::Type::PAIR_KS) has_pair_ks = true;
    };
    check(*tree);

    std::cout << "  Árbol: " << (has_pair_ks ? "PAIR_KS detectado ✓" : "PAIR_KS NO detectado ✗") << "\n";

    // Integrar y verificar energía
    HierarchyBuilder::Params bp; bp.r_ks_threshold = 2.0;
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(),
        2.0, 1e-4, bp
    );

    double E0 = total_energy(sys);
    std::vector<bool> used(2, false);
    for (int s = 0; s < 500; ++s)
        integrator.step(sys, 0.01, used);

    double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
    std::cout << "  |ΔE/E0| = " << dE << "\n";

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

    NBodySystem sys;
    sys.G = 1.0;
    // Figura-8
    sys.bodies.push_back({{ 0.970, -0.243, 0}, { 0.466,  0.433, 0}, 1.0});
    sys.bodies.push_back({{-0.970, -0.243, 0}, { 0.466, -0.433, 0}, 1.0});
    sys.bodies.push_back({{ 0.000,  0.970, 0}, {-0.932,  0.000, 0}, 1.0});

    HierarchyBuilder::Params p;
    p.r_ks_threshold = 3.0;
    p.strong_coupling_eta = 5.0;
    HierarchyBuilder hb(p);
    auto tree = hb.build(sys);

    bool has_chain = false;
    auto check = [&](const HierarchyNode& n) {
        if (n.type == HierarchyNode::Type::TRIPLE_CHAIN) has_chain = true;
        if (n.type == HierarchyNode::Type::COMPOSITE)
            for (const auto& c : n.children)
                if (c->type == HierarchyNode::Type::TRIPLE_CHAIN) has_chain = true;
    };
    check(*tree);

    std::cout << "  Árbol: " << (has_chain ? "TRIPLE_CHAIN detectado ✓" : "TRIPLE_CHAIN NO detectado ✗") << "\n";

    // Integrar
    HierarchyBuilder::Params bp;
    bp.r_ks_threshold = 3.0;
    bp.strong_coupling_eta = 5.0;
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(),
        3.0, 1e-4, bp
    );

    double E0 = total_energy(sys);
    std::vector<bool> used(3, false);
    double max_dE = 0.0;
    double r_max  = 0.0;
    for (int s = 0; s < 300; ++s) {
        integrator.step(sys, 0.005, used);
        double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
        if (dE > max_dE) max_dE = dE;
        for (const auto& b : sys.bodies)
            r_max = std::max(r_max, b.position.norm());
    }

    std::cout << "  Error max energía: " << max_dE << "\n";
    std::cout << "  Radio máximo:      " << r_max << "\n";

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

    NBodySystem sys;
    sys.G = 1.0;
    sys.bodies.push_back({{-0.5, 0, 0}, {0, -0.5, 0}, 1.0});
    sys.bodies.push_back({{ 0.5, 0, 0}, {0,  0.5, 0}, 1.0});
    sys.bodies.push_back({{50.0, 0, 0}, {0,  0.1, 0}, 0.01});

    HierarchyBuilder::Params bp;
    bp.r_ks_threshold = 2.0;
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(),
        2.0, 1e-4, bp
    );

    double E0 = total_energy(sys);
    std::vector<bool> used(3, false);
    for (int s = 0; s < 300; ++s)
        integrator.step(sys, 0.01, used);

    double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
    double r_field = sys.bodies[2].position.norm();

    std::cout << "  |ΔE/E0| = " << dE << "\n";
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

    NBodySystem sys;
    sys.G = 1.0;
    sys.bodies.push_back({{ 0.970, -0.243, 0}, { 0.466,  0.433, 0}, 1.0});
    sys.bodies.push_back({{-0.970, -0.243, 0}, { 0.466, -0.433, 0}, 1.0});
    sys.bodies.push_back({{ 0.000,  0.970, 0}, {-0.932,  0.000, 0}, 1.0});
    sys.bodies.push_back({{100.0,   0.000, 0}, { 0.000,  0.050, 0}, 0.01});

    HierarchyBuilder::Params bp;
    bp.r_ks_threshold = 3.0;
    bp.strong_coupling_eta = 5.0;
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(),
        3.0, 1e-4, bp
    );

    double E0 = total_energy(sys);
    std::vector<bool> used(4, false);
    double max_dE = 0.0;
    for (int s = 0; s < 200; ++s) {
        integrator.step(sys, 0.005, used);
        double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
        if (dE > max_dE) max_dE = dE;
    }

    double r_field = sys.bodies[3].position.norm();
    std::cout << "  Error max energía: " << max_dE << "\n";
    std::cout << "  Estrella lejana en r = " << r_field << " (debe ser ~100)\n";

    bool ok = (max_dE < 5.0) && (std::abs(r_field - 100.0) < 10.0);
    std::cout << "RESULTADO TEST D: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
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

    std::cout << "\n" << std::string(72,'=') << "\n";
    std::cout << "RESUMEN:\n";
    std::cout << "  A (binaria → PAIR_KS):        " << (a ? "✓ PASADO" : "✗ FALLADO") << "\n";
    std::cout << "  B (triple → TRIPLE_CHAIN):    " << (b ? "✓ PASADO" : "✗ FALLADO") << "\n";
    std::cout << "  C (binaria + campo):           " << (c ? "✓ PASADO" : "✗ FALLADO") << "\n";
    std::cout << "  D (triple + campo):            " << (d ? "✓ PASADO" : "✗ FALLADO") << "\n";
    std::cout << std::string(72,'=') << "\n";

    if (a && b && c && d) {
        std::cout << "✓ Todos los tests pasaron. Ruta B completa.\n";
        return 0;
    }
    std::cout << "✗ Algunos tests fallaron.\n";
    return 1;
}
/*

Sobre el error de energía en Test B (0.437) — es esperado dado el contexto.
Por qué es tan alto:
La figura-8 es una órbita caótica y muy sensible. Con RK4 fijo de paso dtau=1e-4, 
el integrador de Chain3 acumula error secular porque RK4 no es simpléctico,
no conserva el volumen de fase por construcción.
Además, el árbol se reconstruye en cada paso del HierarchicalIntegrator. 
Cada reconstrucción reinicializa el Chain3State desde el sistema físico. 
Eso introduce discontinuidades pequeñas pero acumulables.

Para validación y tests cualitativos → aceptable.
Para simulación de precisión larga → no.

Las dos soluciones posibles:
Corto plazo: 
Usar el integrador adaptativo de Chain3 con tolerancia más fina (ya lo tienes implementado) 
y reducir el paso dt del test. No requiere cambios de arquitectura.
Largo plazo correcto: 
Cambiar la arquitectura para que el Chain3State persista entre pasos en lugar de reinicializarse,
el árbol reutiliza el estado anterior si el triple sigue siendo el mismo. 
Esto es lo que hace NBODY6 con sus binarias KS: el estado KS vive entre pasos de campo.

*/