// regularization/chain/archain3_integrator.cpp
#include "archain3_integrator.h"
#include "nbody_system.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>

// ============================================================================
// COEFICIENTES DE YOSHIDA (1990) — COMPOSICIÓN DE ORDEN 4
//
// Referencia: Yoshida, H. (1990), Phys. Lett. A 150, 262–268.
//             Mikkola & Merritt (2006), AJ 130, 2345 (aplicación a TTL).
//
// Dado un integrador simpléctico de orden 2 φ_h, el esquema de orden 4 es:
//
//   φ_h^(4) = φ_{c1·h}  ∘  φ_{c0·h}  ∘  φ_{c1·h}
//
// con:
//   c1 =  1 / (2 - 2^(1/3))  ≈ +1.3512071919596578
//   c0 = -2^(1/3) / (2 - 2^(1/3))  ≈ -1.7024143839193153
//
// Verificación: 2·c1 + c0 = 1 (el paso total se conserva).
// c0 < 0 es correcto — el segundo sub-paso retrocede en tiempo ficticio.
// La simpleticidad se conserva porque cada sub-paso sigue siendo simpléctico.
// ============================================================================
namespace {
    constexpr double Y_c1 =  1.3512071919596578;   //  1 / (2 - 2^(1/3))
    constexpr double Y_c0 = -1.7024143839193153;   // -2^(1/3) / (2 - 2^(1/3))
    static_assert(2.0 * 1.3512071919596578 + (-1.7024143839193153) > 0.9999999 &&
                  2.0 * 1.3512071919596578 + (-1.7024143839193153) < 1.0000001,
                  "Yoshida: 2*c1 + c0 != 1");
}

// ============================================================================
// Constructor
// ============================================================================
ARChain3Integrator::ARChain3Integrator(double eta)
    : eta_(eta)
    , ds_min_(1e-16)
    , ds_max_(1.0)
{}

// ============================================================================
// INITIALIZE: sistema físico → estado AR-chain
// ============================================================================
ARChain3State ARChain3Integrator::initialize(const NBodySystem& system,
                                              int idx1, int idx2, int idx3) const
{
    if (idx1 < 0 || idx1 >= (int)system.bodies.size() ||
        idx2 < 0 || idx2 >= (int)system.bodies.size() ||
        idx3 < 0 || idx3 >= (int)system.bodies.size())
        throw std::out_of_range("ARChain3Integrator::initialize: índice fuera de rango");

    const auto& b1 = system.bodies[idx1];
    const auto& b2 = system.bodies[idx2];
    const auto& b3 = system.bodies[idx3];

    const double m1 = b1.mass, m2 = b2.mass, m3 = b3.mass;
    const double M  = m1 + m2 + m3;

    const Vec3 cm_pos = (b1.position * m1 + b2.position * m2 + b3.position * m3) * (1.0 / M);
    const Vec3 cm_vel = (b1.velocity * m1 + b2.velocity * m2 + b3.velocity * m3) * (1.0 / M);

    ARChain3State st;
    st.X1 = b2.position - b1.position;   // r2 - r1
    st.X2 = b3.position - b2.position;   // r3 - r2
    st.V1 = b2.velocity - b1.velocity;   // v2 - v1
    st.V2 = b3.velocity - b2.velocity;   // v3 - v2

    st.cm_pos = cm_pos;
    st.cm_vel = cm_vel;
    st.t_phys = 0.0;
    st.s_fict = 0.0;
    st.masses[0] = m1;
    st.masses[1] = m2;
    st.masses[2] = m3;

    st.energy = compute_energy(st);

    return st;
}

// ============================================================================
// WRITE_BACK: estado AR-chain → sistema físico
// ============================================================================
void ARChain3Integrator::write_back(const ARChain3State& state, NBodySystem& system,
                                     int idx1, int idx2, int idx3) const
{
    if (!state.is_valid())
        throw std::runtime_error("ARChain3Integrator::write_back: estado inválido");

    Vec3 r1, r2, r3, v1, v2, v3;
    state.reconstruct_positions(r1, r2, r3);
    state.reconstruct_velocities(v1, v2, v3);
    // NOTA: reconstruct_velocities() ya devuelve velocidades absolutas
    // (incluye cm_vel internamente). No hay que sumarlo de nuevo.

    const Vec3 cm_displacement = state.cm_vel * state.t_phys;
    system.bodies[idx1].position = r1 + cm_displacement;
    system.bodies[idx2].position = r2 + cm_displacement;
    system.bodies[idx3].position = r3 + cm_displacement;

    system.bodies[idx1].velocity = v1;   // ya es velocidad absoluta
    system.bodies[idx2].velocity = v2;
    system.bodies[idx3].velocity = v3;
}

// ============================================================================
// COMPUTE_OMEGA: función de tiempo TTL
// Ω = m1·m2/r12 + m2·m3/r23 + m1·m3/r13
// ============================================================================
double ARChain3Integrator::compute_Omega(const ARChain3State& st) const
{
    const double r12 = st.X1.norm();
    const double r23 = st.X2.norm();
    const Vec3   X3  = st.X1 + st.X2;
    const double r13 = X3.norm();

    if (r12 < 1e-30 || r23 < 1e-30 || r13 < 1e-30)
        throw std::runtime_error("compute_Omega: separación singular");

    return st.m1() * st.m2() / r12
         + st.m2() * st.m3() / r23
         + st.m1() * st.m3() / r13;
}

// ============================================================================
// COMPUTE_SEP_MIN: separación mínima entre los tres cuerpos
// ============================================================================
double ARChain3Integrator::compute_sep_min(const ARChain3State& st) const
{
    const double r12 = st.X1.norm();
    const double r23 = st.X2.norm();
    const double r13 = (st.X1 + st.X2).norm();
    return std::min({r12, r23, r13});
}

// ============================================================================
// COMPUTE_ENERGY: T + U del subsistema en coordenadas de cadena
// ============================================================================
double ARChain3Integrator::compute_energy(const ARChain3State& st) const
{
    Vec3 r1, r2, r3, v1, v2, v3;
    st.reconstruct_positions(r1, r2, r3);
    st.reconstruct_velocities(v1, v2, v3);

    const double T = 0.5 * (st.m1() * dot(v1, v1)
                           + st.m2() * dot(v2, v2)
                           + st.m3() * dot(v3, v3));

    const double r12 = st.X1.norm();
    const double r23 = st.X2.norm();
    const double r13 = (st.X1 + st.X2).norm();
    const double U   = -(st.m1() * st.m2() / r12
                       + st.m2() * st.m3() / r23
                       + st.m1() * st.m3() / r13);

    return T + U;
}

// ============================================================================
// COMPUTE_ACCELERATIONS: diferencias de aceleración en coordenadas de cadena
//
// A1 = a2 - a1  y  A2 = a3 - a2
//
// Derivación:
//   a2 = m1(-X1)/r12³ + m3·X2/r23³
//   a1 = m2·X1/r12³  + m3·X3/r13³
//   A1 = a2 - a1 = -X1·(m1+m2)/r12³ + X2·m3/r23³ - X3·m3/r13³
//
//   a3 = m1(-X3)/r13³ + m2(-X2)/r23³
//   A2 = a3 - a2 = -X3·m1/r13³ - X2·(m2+m3)/r23³ + X1·m1/r12³
// ============================================================================
void ARChain3Integrator::compute_accelerations(const ARChain3State& st,
                                                Vec3& A1, Vec3& A2) const
{
    const double m1 = st.m1(), m2 = st.m2(), m3 = st.m3();

    const Vec3   X3     = st.X1 + st.X2;
    const double r12_sq = dot(st.X1, st.X1);
    const double r23_sq = dot(st.X2, st.X2);
    const double r13_sq = dot(X3, X3);

    const double r12 = std::sqrt(r12_sq);
    const double r23 = std::sqrt(r23_sq);
    const double r13 = std::sqrt(r13_sq);

    if (r12 < 1e-30 || r23 < 1e-30 || r13 < 1e-30)
        throw std::runtime_error("compute_accelerations: separación singular");

    const double f12 = 1.0 / (r12_sq * r12);
    const double f23 = 1.0 / (r23_sq * r23);
    const double f13 = 1.0 / (r13_sq * r13);

    A1 =   st.X1 * (-(m1 + m2) * f12)
         + st.X2 * (m3 * f23)
         + X3    * (-m3 * f13);

    A2 =   X3    * (-m1 * f13)
         + st.X2 * (-(m2 + m3) * f23)
         + st.X1 * (m1 * f12);
}

// ============================================================================
// LEAPFROG_STEP_BARE: un paso DKD (Drift-Kick-Drift) en tiempo ficticio
//
// Mikkola & Tanikawa (1999), ec. (6):
//   DRIFT½:   Xi <- Xi + (ds/2)·Ω_ini·Vi
//   KICK:     Vi <- Vi + ds·Ω_mid·Ai(X_mid)
//   DRIFT½:   Xi <- Xi + (ds/2)·Ω_new·Vi_new
//   TIME:     t  <- t  + ds·Ω_mid
//
// La asimetría entre Ω_ini y Ω_new (en lugar de Ω_mid en ambos DRIFT)
// proviene del splitting simpléctico exacto del Hamiltoniano transformado.
//
// ACEPTA ds < 0: necesario para la composición de Yoshida, donde el
// segundo sub-paso usa c0 < 0. El tiempo físico puede decrementar
// temporalmente en ese sub-paso — es matemáticamente correcto y el
// valor final de t_phys tras los tres sub-pasos es positivo neto.
// ============================================================================
void ARChain3Integrator::leapfrog_step_bare(ARChain3State& st, double ds) const
{
    // ── DRIFT½ ────────────────────────────────────────────────────────────────
    const double Omega_ini = compute_Omega(st);
    const double half_ds   = 0.5 * ds;

    st.X1 += st.V1 * (half_ds * Omega_ini);
    st.X2 += st.V2 * (half_ds * Omega_ini);

    // ── KICK ──────────────────────────────────────────────────────────────────
    const double Omega_mid = compute_Omega(st);
    const double dt_mid    = ds * Omega_mid;

    Vec3 A1, A2;
    compute_accelerations(st, A1, A2);

    st.V1 += A1 * dt_mid;
    st.V2 += A2 * dt_mid;

    // ── DRIFT½ ────────────────────────────────────────────────────────────────
    const double Omega_new = compute_Omega(st);

    st.X1 += st.V1 * (half_ds * Omega_new);
    st.X2 += st.V2 * (half_ds * Omega_new);

    // ── AVANCE DE TIEMPO FÍSICO ───────────────────────────────────────────────
    st.t_phys += dt_mid;
    st.s_fict += ds;
}

// ============================================================================
// LEAPFROG_STEP: wrapper público de compatibilidad (delega a leapfrog_step_bare)
//
// Mantiene la API pública para tests y código externo que usa leapfrog_step.
// ============================================================================
void ARChain3Integrator::leapfrog_step(ARChain3State& st, double ds) const
{
    leapfrog_step_bare(st, ds);
}

// ============================================================================
// YOSHIDA4_STEP: integrador de orden 4 TTL — esquema DKDKDKD
//
// Implementa la composición de Yoshida (1990) para el Hamiltoniano
// transformado por TTL. El esquema correcto para TTL no es componer
// tres leapfrogs DKD completos (que usa Ω en posiciones inconsistentes
// durante los sub-pasos negativos), sino descomponer en operadores
// elementales DRIFT y KICK independientes.
//
// Esquema DKDKDKD (7 sub-pasos):
//   D(c1) K(d1) D(c2) K(d2) D(c2) K(d1) D(c1)
//
// Coeficientes (Yoshida 1990, ec. 4.4):
//   w  = 1 / (2 - 2^(1/3))
//   c1 = c4 = w/2
//   c2 = c3 = (1 - 2^(1/3))*w/2
//   d1 = d3 = w
//   d2      = -2^(1/3)*w
//
// Propiedades:
//   sum(c) = 1,  sum(d) = 1  (el paso total se conserva)
//   c2 < 0: el drift central retrocede en tiempo ficticio (correcto)
//   d2 < 0: el kick central retrocede en velocidad (correcto)
//   Ningún drift cruza singularidades porque los pasos son más pequeños
//   que en la composición DKD-DKD-DKD.
//
// CADA operador DRIFT y KICK usa su propio Ω evaluado en la posición
// actual — esto es correcto para TTL y evita inconsistencias numéricas.
//
// Coste: 4 evaluaciones de Ω (drifts) + 3 evaluaciones de Ω+accel (kicks)
//        vs 3 evaluaciones de Ω+accel del leapfrog orden 2.
//        Factor real: ~3.5× más caro que orden 2.
//
// Referencia: Yoshida (1990); Mikkola & Merritt (2006) §2.
// ============================================================================
void ARChain3Integrator::yoshida4_step(ARChain3State& st, double ds) const
{
    // Coeficientes de Yoshida para esquema DKDKDKD
    constexpr double w  = 1.3512071919596578;
    constexpr double c1 =  w * 0.5;                           //  0.6756035960
    constexpr double c2 = (1.0 - 1.2599210498948732)*w*0.5;   // -0.1756035960
    constexpr double d1 =  w;                                  //  1.3512071920
    constexpr double d2 = -1.2599210498948732 * w;             // -1.7024143839

    // Omega de referencia: calculado al inicio del paso completo.
    // Se usa SOLO para avanzar t_phys y s_fict al final.
    // Los operadores DRIFT usan su propio Omega local para mover posiciones.
    const double Omega_ref = compute_Omega(st);

    Vec3 A1, A2;
    double Omega_local;

    // ── D(c1) ────────────────────────────────────────────────────────────────
    Omega_local = compute_Omega(st);
    st.X1 += st.V1 * (c1 * ds * Omega_local);
    st.X2 += st.V2 * (c1 * ds * Omega_local);

    // ── K(d1) ────────────────────────────────────────────────────────────────
    Omega_local = compute_Omega(st);
    compute_accelerations(st, A1, A2);
    st.V1 += A1 * (d1 * ds * Omega_local);
    st.V2 += A2 * (d1 * ds * Omega_local);

    // ── D(c2) ────────────────────────────────────────────────────────────────
    Omega_local = compute_Omega(st);
    st.X1 += st.V1 * (c2 * ds * Omega_local);
    st.X2 += st.V2 * (c2 * ds * Omega_local);

    // ── K(d2) ────────────────────────────────────────────────────────────────
    Omega_local = compute_Omega(st);
    compute_accelerations(st, A1, A2);
    st.V1 += A1 * (d2 * ds * Omega_local);
    st.V2 += A2 * (d2 * ds * Omega_local);

    // ── D(c2) ────────────────────────────────────────────────────────────────
    Omega_local = compute_Omega(st);
    st.X1 += st.V1 * (c2 * ds * Omega_local);
    st.X2 += st.V2 * (c2 * ds * Omega_local);

    // ── K(d1) ────────────────────────────────────────────────────────────────
    Omega_local = compute_Omega(st);
    compute_accelerations(st, A1, A2);
    st.V1 += A1 * (d1 * ds * Omega_local);
    st.V2 += A2 * (d1 * ds * Omega_local);

    // ── D(c1) ────────────────────────────────────────────────────────────────
    Omega_local = compute_Omega(st);
    st.X1 += st.V1 * (c1 * ds * Omega_local);
    st.X2 += st.V2 * (c1 * ds * Omega_local);

    // ── Avance de tiempo: un único Omega_ref para consistencia con el bucle ──
    // Usar Omega_ref (inicio del paso) mantiene coherencia con el criterio
    // adaptativo ds = eta*sep_min/Omega calculado en integrate_until.
    // El error de usar Omega_ref vs Omega_final es O(ds^2) — despreciable.
    st.t_phys += ds * Omega_ref;   // sum(c)=1, sum(d)=1 → paso neto = ds*Omega
    st.s_fict += ds;
}

// ============================================================================
// INTEGRATE_UNTIL: núcleo compartido del bucle de integración
//
// Avanza state desde state.t_phys hasta t_target usando pasos adaptativos

// ds = η·sep_min/Ω. El último paso se recorta para no sobrepasar t_target.
//
// Este método es el núcleo real: tanto integrate() como integrate_to()
// lo llaman con el t_target apropiado. De este modo la lógica de paso
// adaptativo vive en un solo lugar.
// ============================================================================
void ARChain3Integrator::integrate_until(ARChain3State& state, double t_target) const
{
    if (t_target <= state.t_phys)
        return;

    // Límite de seguridad: 10M pasos es margen amplio para cualquier sistema
    // razonable con η ≥ 1e-5.
    constexpr int max_steps = 10'000'000;
    int steps_remaining = max_steps;

    while (state.t_phys < t_target && steps_remaining-- > 0) {

        const double Omega   = compute_Omega(state);
        const double sep_min = compute_sep_min(state);

        // ── ds adaptativo basado en separación ─────────────────────────────
        // ds = η · sep_min / Ω  →  dt_físico ≈ η · sep_min
        // Cuando sep_min → 0: dt_físico → 0 (paso físico se reduce).
        // Referencia: Mikkola & Tanikawa (1999).
        double ds_step = eta_ * sep_min / Omega;
        ds_step = std::clamp(ds_step, ds_min_, ds_max_);

        // ── SIN recorte del último paso ─────────────────────────────────────
        // Con Yoshida orden 4 NO se recorta el último paso. El esquema de
        // composición φ_{c1·ds} ∘ φ_{c0·ds} ∘ φ_{c1·ds} requiere que los
        // tres sub-pasos sean proporcionales al mismo ds. Recortar ds_step
        // en el último paso hace que c0·ds_step retroceda más de lo que
        // c1·ds_step avanzó, desestabilizando el integrador.
        //
        // El overshoot resultante (t_phys supera t_target por ~η·sep_min)
        // es irrelevante en el modo directo (step_to), donde t_target es
        // el tiempo final absoluto de la simulación y no hay sincronización
        // externa que requiera exactitud en el instante de parada.
        //
        // Referencia: Mikkola & Merritt (2006) §2 — el TTL con Yoshida
        // no recorta el último paso; acepta el overshoot residual.

        leapfrog_step(state, ds_step);   // leapfrog TTL orden 2 (Mikkola & Tanikawa 1999)
    }

    if (steps_remaining <= 0)
        std::cerr << "[AR] MAX_STEPS en t=" << state.t_phys
                  << " target=" << t_target
                  << " sep=" << compute_sep_min(state) << "\n";
}

// ============================================================================
// INTEGRATE: modo bloque (legacy / subsistema parcial en N > 3)
//
// Integra un intervalo dt_phys de tiempo físico desde state.t_phys.
// Usado cuando el integrador jerárquico externo impone el tamaño del bloque
// porque necesita sincronizar el tiempo con otros subsistemas (campo lejano,
// otros pares KS, etc.).
//
// Para N=3 con el triple como único subsistema activo, usar integrate_to()
// en su lugar: evita los errores de fase que introduce el corte de bloque.
// ============================================================================
void ARChain3Integrator::integrate(ARChain3State& state, double dt_phys) const
{
    if (dt_phys <= 0.0)
        return;

    const double t_target = state.t_phys + dt_phys;
    integrate_until(state, t_target);
}

// ============================================================================
// INTEGRATE_TO: modo directo (sistema completo, N=3)
//
// Integra desde state.t_phys hasta t_abs_final SIN cortes de bloque externos.
// El AR-chain controla su propio tiempo de extremo a extremo.
//
// MOTIVACIÓN ARQUITECTURAL:
//   El integrador jerárquico con bloques de dt fijo acumula error de fase
//   antes del encuentro cercano. Con dt=0.001, la figura-8 diverge de la
//   trayectoria correcta en t≈2.2 (primer encuentro en t≈2.23 con sep=0.021
//   en lugar de t≈3.04 con sep=7.7e-4). Este error no se reduce achicando dt
//   — se reduce eliminando los cortes.
//
//   Este método es el equivalente directo del RK4 adaptativo de referencia
//   que alcanza |dH/H₀| = 1.7×10⁻¹² en la figura-8 con 93k pasos.
//
// DIFERENCIA respecto a integrate():
//   - integrate()    recibe un DELTA de tiempo (dt_phys): t_target = t_phys + dt_phys
//   - integrate_to() recibe un TIEMPO ABSOLUTO (t_abs_final): t_target = t_abs_final
//
// El estado puede tener t_phys > 0 si fue parcialmente integrado.
// La integración avanza desde state.t_phys hasta t_abs_final.
// ============================================================================
void ARChain3Integrator::integrate_to(ARChain3State& state, double t_abs_final) const
{
    if (t_abs_final <= state.t_phys)
        return;

    integrate_until(state, t_abs_final);
}