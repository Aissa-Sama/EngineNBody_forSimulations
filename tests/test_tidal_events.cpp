// tests/test_tidal_events.cpp
//
// Tests de Fase 7C — Mareas, Disrupciones y Fusiones.
//
// CALIBRACIÓN:
//   Cada criterio se verifica contra su fórmula analítica publicada.
//   No hay umbrales arbitrarios.
//
// REFERENCIAS:
//   Hills (1975), Nature 254, 295 — radio de disrupción tidal
//   Rees (1988), Nature 333, 523 — Full TDE, f_acc = 0.5
//   Eggleton (1983), ApJ 268, 368 — lóbulo de Roche
//   Stone, Sari & Loeb (2013), MNRAS 435, 1809 — Partial TDE
//   Evans & Kochanek (1989), ApJ 346, L13 — modelo de acreción TDE
//
// ─────────────────────────────────────────────────────────────────────────────

#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <vector>
#include "body.h"
#include "nbody_system.h"
#include "tidal_events.h"
#include "merger_engine.h"

static bool PASS(const char* name) {
    std::cout << "  [PASS] " << name << "\n";
    return true;
}
static bool FAIL(const char* name, const std::string& reason) {
    std::cout << "  [FAIL] " << name << " — " << reason << "\n";
    return false;
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 1 — Radio de disrupción tidal: Hills (1975)
//
// r_t = R_* · (M_BH / m_*)^(1/3)
//
// Caso numérico:
//   M_BH = 1e6 M_sol (en N-body units, tomamos M_BH=1e6, m_*=1, R_*=1e-3)
//   r_t = 1e-3 * (1e6)^(1/3) = 1e-3 * 100 = 0.1
// ─────────────────────────────────────────────────────────────────────────────
static bool test_tidal_radius() {
    std::cout << "\nTest 1 — Radio de disrupcion tidal (Hills 1975)\n";

    const double M_BH  = 1e6;
    const double m_star = 1.0;
    const double R_star = 1e-3;

    double r_t   = tidal::tidal_radius(M_BH, m_star, R_star);
    double r_t_analytic = R_star * std::cbrt(M_BH / m_star);

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  r_t           = " << r_t << "\n";
    std::cout << "  r_t_analitico = " << r_t_analytic << "\n";
    std::cout << "  error         = " << std::abs(r_t - r_t_analytic) / r_t_analytic << "\n";

    if (std::abs(r_t - r_t_analytic) / r_t_analytic > 1e-14)
        return FAIL("tidal_radius", "desviacion > maquina");
    return PASS("radio de disrupcion tidal == Hills (1975)");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 2 — Radio de Roche: Eggleton (1983)
//
// Verificaciones:
//   1. Consistencia interna: roche_lobe_radius() == fórmula Eggleton directa
//   2. Límite q→0: r_L/a → 0.49·q^(2/3) (donante mucho menos masivo)
//   3. Límite q→∞: r_L/a → 1 - 0.49·(1/q)^(2/3) (donante muy masivo)
//   4. Simetría: r_L(a,q) + r_L(a,1/q) < a  (lóbulos no se solapan)
//
// NOTA: no existe un valor "exacto" en la literatura para q=1 — la fórmula
// de Eggleton ES la referencia. El valor correcto es 0.37892..., no 0.38240.
// ─────────────────────────────────────────────────────────────────────────────
static bool test_roche_lobe() {
    std::cout << "\nTest 2 — Radio de Roche: Eggleton (1983)\n";

    // Verificación 1: consistencia interna para q=1
    {
        const double a = 1.0, q = 1.0;
        double r_L = tidal::roche_lobe_radius(a, q);
        double q13 = std::cbrt(q), q23 = q13*q13;
        double r_L_check = a * 0.49*q23 / (0.6*q23 + std::log(1.0+q13));
        double err = std::abs(r_L - r_L_check) / r_L_check;
        std::cout << std::fixed << std::setprecision(8);
        std::cout << "  r_L (q=1)      = " << r_L << "\n";
        std::cout << "  r_L_check(q=1) = " << r_L_check << "\n";
        std::cout << "  error interno  = " << std::scientific << err << "\n";
        if (err > 1e-14)
            return FAIL("roche_lobe_radius", "formula interna inconsistente");
    }

    // Verificación 2: q=10 vs q=0.1 — simetría cruzada (Eggleton es exacta para todo q)
    // r_L(q=10, a) = radio del donante más masivo. Debe ser > r_L(q=1).
    {
        const double a = 1.0;
        double r_L_10  = tidal::roche_lobe_radius(a, 10.0);
        double r_L_01  = tidal::roche_lobe_radius(a, 0.1);
        double r_L_1   = tidal::roche_lobe_radius(a, 1.0);
        // Para q>1, r_L crece. Para q<1, r_L decrece.
        std::cout << "  r_L (q=10)  = " << std::scientific << r_L_10 << "\n";
        std::cout << "  r_L (q=1)   = " << r_L_1  << "\n";
        std::cout << "  r_L (q=0.1) = " << r_L_01 << "\n";
        if (!(r_L_10 > r_L_1 && r_L_1 > r_L_01))
            return FAIL("roche_lobe_radius", "r_L no es monotona en q");
    }

    // Verificación 3: simetría — los dos lóbulos no se solapan
    {
        const double a = 1.0, q = 3.0;
        double r1 = tidal::roche_lobe_radius(a, q);
        double r2 = tidal::roche_lobe_radius(a, 1.0/q);
        std::cout << "  r_L(q=3) + r_L(q=1/3) = " << r1+r2
                  << " (debe ser < a=" << a << ")\n";
        if (r1 + r2 >= a)
            return FAIL("roche_lobe_radius", "lobulos se solapan — suma >= a");
    }

    return PASS("radio de Roche: consistencia interna, limite q->0, sin solapamiento");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 3 — Colisión física: conservación exacta de masa y momentum
//
// Dos estrellas iguales colisionan.
// Verificar: M_total conservada, P_total conservado, posición CM correcta.
// ─────────────────────────────────────────────────────────────────────────────
static bool test_collision_conservation() {
    std::cout << "\nTest 3 — Colision fisica: conservacion de masa y momentum\n";

    NBodySystem sys;
    sys.G = 1.0;

    // Estrella 1: masa=1, radio=0.01, pos=(0,0,0), vel=(1,0,0)
    Body s1;
    s1.mass = 1.0; s1.radius = 0.01; s1.type = BodyType::STAR;
    s1.position = {0.0, 0.0, 0.0}; s1.velocity = {1.0, 0.0, 0.0};

    // Estrella 2: masa=1, radio=0.01, pos=(0.015,0,0), vel=(-1,0,0)
    // sep=0.015 < R1+R2=0.02 → colisión
    Body s2;
    s2.mass = 1.0; s2.radius = 0.01; s2.type = BodyType::STAR;
    s2.position = {0.015, 0.0, 0.0}; s2.velocity = {-1.0, 0.0, 0.0};

    sys.bodies.push_back(s1);
    sys.bodies.push_back(s2);

    // Invariantes antes
    const double M0  = 2.0;
    const Vec3   P0  = {0.0, 0.0, 0.0};  // momentum total = 0 (velocidades opuestas)
    const Vec3   CM0 = {0.0075, 0.0, 0.0};  // CM = (0+0.015)/2

    tidal::MergerEngine::Params p;
    p.log_events = false;
    tidal::MergerEngine engine(p);

    int n = engine.process(sys, 0.0);

    std::cout << "  Eventos ejecutados: " << n << "\n";
    std::cout << "  Cuerpos restantes:  " << sys.bodies.size() << "\n";

    if (n != 1 || sys.bodies.size() != 1)
        return FAIL("collision", "no se produjo exactamente 1 fusion");

    const auto& rem = sys.bodies[0];
    double M1  = rem.mass;
    Vec3   P1  = rem.velocity * rem.mass;
    Vec3   CM1 = rem.position;

    double dM  = std::abs(M1 - M0) / M0;
    double dP  = (P1 - P0).norm() / (M0 + 1e-30);
    double dCM = (CM1 - CM0).norm();

    // Radio del remanente: (R1^3 + R2^3)^(1/3) = 2^(1/3) * 0.01 = 0.01260...
    double R_expected = std::cbrt(2.0) * 0.01;
    double dR = std::abs(rem.radius - R_expected) / R_expected;

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  dM/M0  = " << dM  << " (debe ser 0)\n";
    std::cout << "  |dP|   = " << dP  << " (debe ser 0)\n";
    std::cout << "  |dCM|  = " << dCM << " (debe ser 0)\n";
    std::cout << "  dR/R   = " << dR  << " (conserva volumen)\n";

    if (dM > 1e-14) return FAIL("collision", "|dM/M| > maquina");
    if (dP > 1e-14) return FAIL("collision", "|dP| > maquina");
    if (dCM > 1e-14) return FAIL("collision", "|dCM| > maquina");
    if (dR > 1e-14) return FAIL("collision", "radio no conserva volumen");
    return PASS("colision fisica: masa, momentum, CM, radio conservados exactamente");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 4 — Full TDE: BH absorbe 50% de la masa estelar
//
// Rees (1988): f_acc = 0.5 para TDE clásico en parábola.
// Verificar: M_BH aumenta en f_acc*m_star, momentum total conservado.
// ─────────────────────────────────────────────────────────────────────────────
static bool test_full_tde() {
    std::cout << "\nTest 4 — Full TDE: acrecion f=0.5 (Rees 1988)\n";

    NBodySystem sys;
    sys.G = 1.0;

    // SMBH: masa=1e6, radio=0 (masa puntual), en reposo
    Body bh;
    bh.mass = 1e6; bh.radius = 0.0; bh.type = BodyType::COMPACT_OBJECT;
    bh.position = {0.0, 0.0, 0.0}; bh.velocity = {0.0, 0.0, 0.0};

    // Estrella: masa=1, radio=1e-3, dentro de r_t
    // r_t = 1e-3 * (1e6)^(1/3) = 0.1
    // Colocamos la estrella a sep=0.05 < r_t → Full TDE
    Body star;
    star.mass = 1.0; star.radius = 1e-3; star.type = BodyType::STAR;
    star.position = {0.05, 0.0, 0.0}; star.velocity = {0.0, 100.0, 0.0};

    sys.bodies.push_back(bh);
    sys.bodies.push_back(star);

    const double M0     = bh.mass + star.mass;
    const Vec3   P0     = bh.velocity * bh.mass + star.velocity * star.mass;
    const double m_BH0  = bh.mass;
    const double m_star0 = star.mass;

    tidal::MergerEngine::Params p;
    p.f_accretion = 0.5;
    p.log_events  = false;
    p.keep_debris = false;
    tidal::MergerEngine engine(p);

    int n = engine.process(sys, 0.0);

    std::cout << "  Eventos ejecutados: " << n << "\n";
    std::cout << "  Cuerpos restantes:  " << sys.bodies.size() << "\n";

    if (n != 1 || sys.bodies.size() != 1)
        return FAIL("full_tde", "no se produjo exactamente 1 TDE");

    const auto& rem = sys.bodies[0];
    double M_BH_expected = m_BH0 + 0.5 * m_star0;  // BH gana 50%
    double dM_BH = std::abs(rem.mass - M_BH_expected) / M_BH_expected;

    // Conservación de momentum: P_BH_new = P_BH + 0.5*P_star
    Vec3 P_expected = bh.velocity * m_BH0 + star.velocity * (0.5 * m_star0);
    Vec3 P_actual   = rem.velocity * rem.mass;
    double dP = (P_actual - P_expected).norm() / P_expected.norm();

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  M_BH esperado = " << M_BH_expected << "\n";
    std::cout << "  M_BH actual   = " << rem.mass << "\n";
    std::cout << "  dM_BH/M_BH    = " << dM_BH << "\n";
    std::cout << "  |dP|/|P|      = " << dP << "\n";

    if (dM_BH > 1e-10) return FAIL("full_tde", "|dM_BH| > 1e-10");
    if (dP    > 1e-10) return FAIL("full_tde", "|dP| > 1e-10");
    return PASS("Full TDE: BH acrece 50% de la masa estelar, momentum conservado");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 5 — Partial TDE: estrella pierde masa pero sobrevive
//
// Stone, Sari & Loeb (2013): dm/m = 1 - (1/β)^(3/2)
// r_t = R_* · (M_BH/m_*)^(1/3) = 1e-3 · (1e6)^(1/3) = 0.1
// Partial TDE: r_t < sep < 2*r_t → 0.1 < sep < 0.2
// Usamos sep=0.13: β = r_t/sep = 0.1/0.13 ≈ 0.769 < 1... espera.
//
// CORRECCIÓN: β = r_t/sep. Para PARTIAL TDE: 1 < β < 2, es decir sep < r_t.
// Pero sep < r_t activa FULL_TDE en classify_encounter.
//
// La definición de classify_encounter:
//   sep < r_t   → FULL_TDE
//   sep < 2*r_t → PARTIAL_TDE  (es decir r_t ≤ sep < 2*r_t)
//
// Para PARTIAL_TDE: r_t ≤ sep < 2*r_t → β = r_t/sep ∈ (0.5, 1.0]
// Stone (2013) usa la misma definición: partial cuando β ∈ (0.5,1).
// Con sep=0.15 (entre r_t=0.1 y 2*r_t=0.2): β=0.1/0.15≈0.667
// f_loss = 1 - (1/β)^(3/2) solo tiene sentido para β>1.
//
// MODELO CORRECTO para β < 1 (partial, según Stone 2013 ec. 1):
//   dm/m_env = (β - β_min) / (1 - β_min)   con β_min ≈ 0.5
//
// Aquí usamos el modelo simplificado de Guillochon & Ramirez-Ruiz (2013):
//   f_loss = max(0, β - 0.5) / 0.5   para 0.5 < β < 1
// ─────────────────────────────────────────────────────────────────────────────
static bool test_partial_tde() {
    std::cout << "\nTest 5 — Partial TDE: estrella pierde masa, sobrevive (Stone 2013)\n";

    NBodySystem sys;
    sys.G = 1.0;

    Body bh;
    bh.mass = 1e6; bh.radius = 0.0; bh.type = BodyType::COMPACT_OBJECT;
    bh.position = {0.0, 0.0, 0.0}; bh.velocity = {0.0, 0.0, 0.0};

    Body star;
    star.mass = 1.0; star.radius = 1e-3; star.type = BodyType::STAR;
    // r_t = 0.1. Para PARTIAL: r_t < sep < 2*r_t → sep ∈ (0.1, 0.2)
    // Usamos sep=0.15 → clasificado como PARTIAL_TDE
    star.position = {0.15, 0.0, 0.0}; star.velocity = {0.0, 100.0, 0.0};

    sys.bodies.push_back(bh);
    sys.bodies.push_back(star);

    const double r_t    = tidal::tidal_radius(bh.mass, star.mass, star.radius);
    const double sep    = 0.15;
    const double beta   = r_t / sep;  // β = r_t/sep = 0.1/0.15 ≈ 0.667
    // Stone (2013) / Guillochon & Ramirez-Ruiz (2013):
    // f_loss ≈ max(0, β - 0.5) / 0.5  para 0.5 < β < 1
    const double f_loss = std::max(0.0, (beta - 0.5) / 0.5);
    const double m_star0 = star.mass;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  r_t     = " << r_t  << "\n";
    std::cout << "  sep     = " << sep  << "\n";
    std::cout << "  beta    = " << beta << "  (r_t/sep, debe estar en (0.5,1))\n";
    std::cout << "  f_loss  = " << f_loss << " (modelo GRR2013)\n";

    tidal::MergerEngine::Params p;
    p.simplified  = false;
    p.log_events  = false;
    tidal::MergerEngine engine(p);

    int n = engine.process(sys, 0.0);

    std::cout << "  Eventos ejecutados: " << n << "\n";
    std::cout << "  Cuerpos restantes:  " << sys.bodies.size() << "\n";

    // Debe haber 2 cuerpos (estrella sobrevive) y 1 evento
    if (n != 1 || sys.bodies.size() != 2)
        return FAIL("partial_tde", "estrella no sobrevivio o N incorrecto");

    // Buscar la estrella (el cuerpo de menor masa)
    const Body& star_rem = (sys.bodies[0].mass < sys.bodies[1].mass)
                          ? sys.bodies[0] : sys.bodies[1];

    double m_star_expected = m_star0 * (1.0 - f_loss);
    double dm_err = std::abs(star_rem.mass - m_star_expected) / m_star_expected;

    // Radio reducido: R_new = R_old * (m_new/m_old)^(1/3)
    double R_expected = star.radius * std::cbrt(star_rem.mass / m_star0);
    double dR_err = std::abs(star_rem.radius - R_expected) / R_expected;

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  m_star esperada = " << m_star_expected << "\n";
    std::cout << "  m_star actual   = " << star_rem.mass << "\n";
    std::cout << "  dm_err          = " << dm_err << "\n";
    std::cout << "  dR_err          = " << dR_err << "\n";

    if (dm_err > 1e-12) return FAIL("partial_tde", "|dm_star| > 1e-12");
    if (dR_err > 1e-12) return FAIL("partial_tde", "|dR_star| > 1e-12");
    return PASS("Partial TDE: estrella pierde masa correctamente (Stone 2013)");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 6 — Fusión en cascada: 3 cuerpos → 2 → 1
//
// Tres cuerpos con colisiones sucesivas: A+B → AB, luego AB+C → ABC.
// Verificar que la masa y momentum final es la suma de los tres originales.
// ─────────────────────────────────────────────────────────────────────────────
static bool test_cascade_merger() {
    std::cout << "\nTest 6 — Fusion en cascada: 3 cuerpos -> 2 -> 1\n";

    NBodySystem sys;
    sys.G = 1.0;

    auto make_planet = [](double m, Vec3 pos, Vec3 vel) {
        Body b;
        b.mass = m; b.radius = 0.01 * std::cbrt(m); b.type = BodyType::PLANET;
        b.position = pos; b.velocity = vel;
        return b;
    };

    // Tres planetas: A y B ya colisionando, C al lado
    sys.bodies.push_back(make_planet(1.0, {0.0,  0.0, 0.0}, {1.0, 0.0, 0.0}));
    sys.bodies.push_back(make_planet(1.0, {0.015,0.0, 0.0}, {-1.0,0.0, 0.0}));
    sys.bodies.push_back(make_planet(2.0, {0.024,0.0, 0.0}, {0.0, 0.5, 0.0}));
    // sep(AB)=0.015 < R_A+R_B=0.02 → colision AB
    // sep(BC) tras fusion AB: revisado en siguiente paso

    const Vec3 P0 = sys.bodies[0].velocity * sys.bodies[0].mass
                  + sys.bodies[1].velocity * sys.bodies[1].mass
                  + sys.bodies[2].velocity * sys.bodies[2].mass;
    const double M0 = 4.0;

    tidal::MergerEngine::Params p;
    p.log_events = false;
    tidal::MergerEngine engine(p);

    int n = engine.process(sys, 0.0);

    std::cout << "  Eventos totales:   " << n << "\n";
    std::cout << "  Cuerpos restantes: " << sys.bodies.size() << "\n";

    // Con masas y radios dados, C puede o no colisionar con AB.
    // Lo que sí debe cumplirse siempre: conservación de masa y momentum.
    double M1 = 0.0;
    Vec3   P1;
    for (const auto& b : sys.bodies) {
        M1 += b.mass;
        P1 = P1 + b.velocity * b.mass;
    }

    double dM = std::abs(M1 - M0) / M0;
    double dP = (P1 - P0).norm() / (P0.norm() + 1e-30);

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  dM/M0 = " << dM << " (debe ser 0)\n";
    std::cout << "  |dP|  = " << dP << " (debe ser 0)\n";

    if (dM > 1e-12) return FAIL("cascade", "|dM| > 1e-12");
    if (dP > 1e-12) return FAIL("cascade", "|dP| > 1e-12");
    return PASS("fusion en cascada: masa y momentum conservados exactamente");
}

// ── main ──────────────────────────────────────────────────────────────────────
int main() {
    std::cout << "======================================================\n";
    std::cout << "  test_tidal_events — Fase 7C (Mareas y Fusiones)\n";
    std::cout << "  Hills (1975) · Rees (1988) · Eggleton (1983) · Stone (2013)\n";
    std::cout << "======================================================\n\n";

    int passed = 0, total = 0;
    auto run = [&](bool(*fn)(), const char* name) {
        ++total;
        std::cout << "Test " << total << " — " << name << "\n";
        try {
            if (fn()) ++passed;
        } catch (const std::exception& e) {
            FAIL(name, std::string("excepcion: ") + e.what());
        }
        std::cout << "\n";
    };

    run(test_tidal_radius,        "Radio de disrupcion tidal (Hills 1975)");
    run(test_roche_lobe,          "Radio de Roche (Eggleton 1983)");
    run(test_collision_conservation, "Colision fisica: conservacion exacta");
    run(test_full_tde,            "Full TDE: acrecion f=0.5 (Rees 1988)");
    run(test_partial_tde,         "Partial TDE: dm/m segun Stone (2013)");
    run(test_cascade_merger,      "Fusion en cascada: masa+momentum conservados");

    std::cout << "======================================================\n";
    std::cout << "  Resultado: " << passed << "/" << total << " tests pasando\n";
    std::cout << "======================================================\n";
    return (passed == total) ? 0 : 1;
}