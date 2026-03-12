// regularization/chain/chain3_integrator.cpp
// ============================================================================
// CHAIN REGULARIZATION PARA 3 CUERPOS — Mikkola & Aarseth (1993)
//
// CONVENCIÓN T_MIKKOLA:
//   T = D1 + D2 = A1·W1 + A2·W2   (SIN factor 0.5)
//   state.energy = T_Mikkola + U   (constante del movimiento regularizado)
//   L = T_Mikkola - U > 0,  g = 1/L
// ============================================================================

#include "chain3_integrator.h"
#include "chain_ks.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

// ============================================================================
// FUNCIONES AUXILIARES LOCALES
// ============================================================================

// Restricción bilineal Q·P = 0
static void enforce_bilinear(Vec4& Q, Vec4& P) {
    double q2 = Q.norm2();
    if (q2 < 1e-30) return;
    double v = (Q.x*P.x + Q.y*P.y + Q.z*P.z + Q.w*P.w) / q2;
    P.x -= v * Q.x;
    P.y -= v * Q.y;
    P.z -= v * Q.z;
    P.w -= v * Q.w;
}

// ============================================================================
// CONSTRUCTOR
// ============================================================================
Chain3Integrator::Chain3Integrator(double internal_dt)
    : internal_dt(internal_dt) {}

// ============================================================================
// INICIALIZACIÓN: Físico → Cadena (VERSIÓN CORREGIDA)
// ============================================================================
Chain3State Chain3Integrator::initialize(const NBodySystem& system,
                                         int idx1, int idx2, int idx3) {
    if (idx1 < 0 || idx1 >= (int)system.bodies.size() ||
        idx2 < 0 || idx2 >= (int)system.bodies.size() ||
        idx3 < 0 || idx3 >= (int)system.bodies.size())
        throw std::out_of_range("initialize: índice fuera de rango");

    const auto& b1 = system.bodies[idx1];
    const auto& b2 = system.bodies[idx2];
    const auto& b3 = system.bodies[idx3];

    double m1 = b1.mass, m2 = b2.mass, m3 = b3.mass;
    double M = m1 + m2 + m3;

    Vec3 r1 = b1.position, r2 = b2.position, r3 = b3.position;
    Vec3 v1 = b1.velocity, v2 = b2.velocity, v3 = b3.velocity;

    // Centro de masa
    Vec3 r0 = (m1 * r1 + m2 * r2 + m3 * r3) / M;
    Vec3 v0 = (m1 * v1 + m2 * v2 + m3 * v3) / M;

    // Momentos en el marco del CM
    Vec3 p1 = m1 * (v1 - v0);
    Vec3 p2 = m2 * (v2 - v0);
    Vec3 p3 = m3 * (v3 - v0);

    // Vectores de la cadena (eslabones)
    Vec3 R1 = r2 - r1;   // eslabón 1: de cuerpo 1 a cuerpo 2
    Vec3 R2 = r3 - r2;   // eslabón 2: de cuerpo 2 a cuerpo 3

    // Transformación a variables KS (W1 = -p1, W2 = p3)
    Vec3 W1 = -p1;
    Vec3 W2 =  p3;

    // Transformación a coordenadas KS
    Vec4 Q1 = r_to_ks(R1);
    Vec4 Q2 = r_to_ks(R2);
    Vec4 P1 = levi_civita_transpose_multiply(Q1, W1) * 2.0;
    Vec4 P2 = levi_civita_transpose_multiply(Q2, W2) * 2.0;

    // ── OPCIÓN B: Proyección de gauge Q·P = 0 ─────────────────────────────
    // Proyectar P sobre el subespacio perpendicular a Q (variedad Q·P = 0).
    // Esto satisface el constraint bilineal de la transformación KS, que
    // es necesario para que el integrador BS (extrapolación de Richardson)
    // converja correctamente a tiempos largos.
    //
    // La proyección elimina la componente de W paralela a R:
    //   W_gauge = W_físico - (R·W / |R|²) · R_hat
    //
    // Almacenamos los escalares de corrección para recuperar W_físico
    // en write_back():
    //   gauge_s1 = Q1·P1 / |Q1|²   (antes de proyectar)
    //   gauge_s2 = Q2·P2 / |Q2|²
    double q1sq = Q1.norm2();
    double q2sq = Q2.norm2();
    double alpha1 = (Q1.x*P1.x + Q1.y*P1.y + Q1.z*P1.z + Q1.w*P1.w) / q1sq;
    double alpha2 = (Q2.x*P2.x + Q2.y*P2.y + Q2.z*P2.z + Q2.w*P2.w) / q2sq;
    P1 = P1 - Q1 * alpha1;   // ahora Q1·P1 = 0
    P2 = P2 - Q2 * alpha2;   // ahora Q2·P2 = 0

    // Recalcular W en gauge (ks_to_w con P proyectado)
    Vec3 W1g = ks_to_w(Q1, P1);
    Vec3 W2g = ks_to_w(Q2, P2);

    // Matriz T y vectores A (Mikkola) — con W en gauge
    double T11 = 1.0 / m1 + 1.0 / m2;
    double T12 = -1.0 / m2;
    double T22 = 1.0 / m2 + 1.0 / m3;

    Vec3 A1 = T11 * W1g + T12 * W2g;
    Vec3 A2 = T12 * W1g + T22 * W2g;

    // Energía cinética de Mikkola en gauge (sin factor 1/2)
    double T_mikkola = dot(A1, W1g) + dot(A2, W2g);

    // Energía potencial
    double d12 = R1.norm();
    double d23 = R2.norm();
    double d13 = (r3 - r1).norm();
    double U = -(m1 * m2 / d12 + m2 * m3 / d23 + m1 * m3 / d13);

    // Construir estado final
    Chain3State state;
    state.Q1 = Q1;
    state.P1 = P1;
    state.Q2 = Q2;
    state.P2 = P2;
    state.energy = T_mikkola + U;   // Hamiltoniano regularizado en gauge
    state.tau = 0.0;
    state.cm_pos = r0;
    state.cm_vel = v0;
    state.cm_time = 0.0;
    state.masses[0] = m1;
    state.masses[1] = m2;
    state.masses[2] = m3;
    // Correcciones de gauge: alpha = Q·P_original/|Q|²
    // Necesarias en write_back para recuperar W_físico desde W_gauge.
    state.gauge_s1 = alpha1;
    state.gauge_s2 = alpha2;

    return state;
}

// ============================================================================
// RECONSTRUCCIÓN: Cadena → Físico
// ============================================================================
void Chain3Integrator::write_back(const Chain3State& state, NBodySystem& system,
                                   int idx1, int idx2, int idx3) {
    if (!state.is_valid())
        throw std::runtime_error("write_back: estado inválido");

    double m1 = state.m1(), m2 = state.m2(), m3 = state.m3();
    double M = m1 + m2 + m3;

    Vec3 R1 = ks_to_r(state.Q1);
    Vec3 R2 = ks_to_r(state.Q2);
    Vec3 W1 = ks_to_w(state.Q1, state.P1);
    Vec3 W2 = ks_to_w(state.Q2, state.P2);

    // Corrección de gauge: recuperar la componente de W paralela a R
    // que fue eliminada en initialize() al proyectar P sobre la variedad Q·P=0.
    //
    // alpha = Q·P_original/|Q|² = 2·R·W_físico/|R|  (almacenado en gauge_s)
    // ΔW = W_físico - W_gauge = (R·W/|R|²)·R = (alpha/2) · R_hat
    //                                         = alpha · R / (2·|R|)
    //
    // gauge_s = 0 cuando el estado se inicializó sin proyección (RK4).
    if (state.gauge_s1 != 0.0) {
        double r1_norm = R1.norm();
        if (r1_norm > 1e-14)
            W1 = W1 + R1 * (state.gauge_s1 / (2.0 * r1_norm));
    }
    if (state.gauge_s2 != 0.0) {
        double r2_norm = R2.norm();
        if (r2_norm > 1e-14)
            W2 = W2 + R2 * (state.gauge_s2 / (2.0 * r2_norm));
    }

    // Momentos en el marco del CM
    Vec3 p1 = -W1;
    Vec3 p2 =  W1 - W2;
    Vec3 p3 =  W2;

    // Posiciones relativas al origen arbitrario
    Vec3 q1(0, 0, 0);
    Vec3 q2 = q1 + R1;
    Vec3 q3 = q2 + R2;

    // Centro de masa de las posiciones relativas
    Vec3 qcm = (m1 * q1 + m2 * q2 + m3 * q3) / M;

    // Posiciones absolutas
    Vec3 r1 = q1 - qcm + state.cm_pos;
    Vec3 r2 = q2 - qcm + state.cm_pos;
    Vec3 r3 = q3 - qcm + state.cm_pos;

    // Velocidades absolutas
    Vec3 v1 = p1 / m1 + state.cm_vel;
    Vec3 v2 = p2 / m2 + state.cm_vel;
    Vec3 v3 = p3 / m3 + state.cm_vel;

    system.bodies[idx1].position = r1;
    system.bodies[idx1].velocity = v1;
    system.bodies[idx2].position = r2;
    system.bodies[idx2].velocity = v2;
    system.bodies[idx3].position = r3;
    system.bodies[idx3].velocity = v3;
}

// ============================================================================
// DERIVADAS — Mikkola & Aarseth (1993) ec. 51-68 (CON PROTECCIÓN)
// ============================================================================
void Chain3Integrator::compute_derivatives(const Chain3State& state,
                                           const NBodySystem&,
                                           Chain3State& deriv) const {
    const double m1 = state.m1(), m2 = state.m2(), m3 = state.m3();

    const Vec3 R1 = ks_to_r(state.Q1);
    const Vec3 R2 = ks_to_r(state.Q2);
    const Vec3 R3 = R1 + R2;

    const Vec3 W1 = ks_to_w(state.Q1, state.P1);
    const Vec3 W2 = ks_to_w(state.Q2, state.P2);

    const double d12 = R1.norm();
    const double d23 = R2.norm();
    const double d13 = R3.norm();

    // --- PROTECCIÓN CONTRA SINGULARIDADES ---
    if (d12 < 1e-12 || d23 < 1e-12 || d13 < 1e-12) {
        throw std::runtime_error("compute_derivatives: singularidad - distancias demasiado pequeñas");
    }

    const double T11 = 1.0 / m1 + 1.0 / m2;
    const double T12 = -1.0 / m2;
    const double T22 = 1.0 / m2 + 1.0 / m3;

    const Vec3 A1 = T11 * W1 + T12 * W2;
    const Vec3 A2 = T12 * W1 + T22 * W2;

    const double D1 = dot(A1, W1);
    const double D2 = dot(A2, W2);
    const double T_kin = D1 + D2;

    const double U = -(m1 * m2 / d12 + m2 * m3 / d23 + m1 * m3 / d13);
    const double L = T_kin - U;
    const double g = 1.0 / L;

    const double H_reg = T_kin + U;
    const double Gamma = g * (H_reg - state.energy);
    const double Gamma_T = (1.0 - Gamma) * g;
    const double Gamma_U = -(1.0 + Gamma) * g;

    const double c12 = m1 * m2 / (d12 * d12 * d12);
    const double c23 = m2 * m3 / (d23 * d23 * d23);
    const double c13 = m1 * m3 / (d13 * d13 * d13);

    const Vec3 dU_dR1 = c12 * R1 + c13 * R3;
    const Vec3 dU_dR2 = c23 * R2 + c13 * R3;

    const double Q1_2 = state.Q1.norm2();
    const double Q2_2 = state.Q2.norm2();

    const Vec4 dU_dQ1 = levi_civita_transpose_multiply(state.Q1, dU_dR1) * 2.0;
    const Vec4 dU_dQ2 = levi_civita_transpose_multiply(state.Q2, dU_dR2) * 2.0;

    const Vec4 T_P1 = levi_civita_transpose_multiply(state.Q1, A1) * (0.5 / Q1_2);
    const Vec4 T_P2 = levi_civita_transpose_multiply(state.Q2, A2) * (0.5 / Q2_2);

    const Vec4& P1 = state.P1;
    const Vec4& P2 = state.P2;

    // Términos L(P)·A (necesarios para T_Q)
    const Vec4 LP1A1(
        P1.x * A1.x - P1.y * A1.y - P1.z * A1.z,
        P1.y * A1.x + P1.x * A1.y - P1.w * A1.z,
        P1.z * A1.x + P1.w * A1.y + P1.x * A1.z,
        P1.w * A1.x - P1.z * A1.y + P1.y * A1.z
    );

    const Vec4 LP2A2(
        P2.x * A2.x - P2.y * A2.y - P2.z * A2.z,
        P2.y * A2.x + P2.x * A2.y - P2.w * A2.z,
        P2.z * A2.x + P2.w * A2.y + P2.x * A2.z,
        P2.w * A2.x - P2.z * A2.y + P2.y * A2.z
    );

    const Vec4 T_Q1 = (LP1A1 - 4.0 * D1 * state.Q1) * (0.5 / Q1_2);
    const Vec4 T_Q2 = (LP2A2 - 4.0 * D2 * state.Q2) * (0.5 / Q2_2);
    // La derivada de T respecto a Q1 se verifica analítica y numéricamente como:
    // dT1/dQ1_k = [A·(dL/dQ_k)·P] / (2|Q1|²) - 2·D1·Q1_k / |Q1|²

    deriv.Q1 = Gamma_T * T_P1;
    deriv.Q2 = Gamma_T * T_P2;
    deriv.P1 = -(Gamma_T * T_Q1 + Gamma_U * dU_dQ1);
    deriv.P2 = -(Gamma_T * T_Q2 + Gamma_U * dU_dQ2);

    deriv.energy = 0.0;
    deriv.cm_pos = Vec3(0, 0, 0);
    deriv.cm_vel = Vec3(0, 0, 0);
    deriv.cm_time = g;
    deriv.tau = 1.0;
    deriv.masses[0] = deriv.masses[1] = deriv.masses[2] = 0.0;
}

// ============================================================================
// PASO RK4
// ============================================================================
void Chain3Integrator::rk4_step(Chain3State& state, double dtau,
                                 const NBodySystem& full_system) {
    Chain3State k1, k2, k3, k4, tmp;

    compute_derivatives(state, full_system, k1);

    tmp = state;
    tmp.Q1 = state.Q1 + k1.Q1 * (dtau / 2);
    tmp.P1 = state.P1 + k1.P1 * (dtau / 2);
    tmp.Q2 = state.Q2 + k1.Q2 * (dtau / 2);
    tmp.P2 = state.P2 + k1.P2 * (dtau / 2);
    tmp.cm_time = state.cm_time + k1.cm_time * (dtau / 2);
    compute_derivatives(tmp, full_system, k2);

    tmp = state;
    tmp.Q1 = state.Q1 + k2.Q1 * (dtau / 2);
    tmp.P1 = state.P1 + k2.P1 * (dtau / 2);
    tmp.Q2 = state.Q2 + k2.Q2 * (dtau / 2);
    tmp.P2 = state.P2 + k2.P2 * (dtau / 2);
    tmp.cm_time = state.cm_time + k2.cm_time * (dtau / 2);
    compute_derivatives(tmp, full_system, k3);

    tmp = state;
    tmp.Q1 = state.Q1 + k3.Q1 * dtau;
    tmp.P1 = state.P1 + k3.P1 * dtau;
    tmp.Q2 = state.Q2 + k3.Q2 * dtau;
    tmp.P2 = state.P2 + k3.P2 * dtau;
    tmp.cm_time = state.cm_time + k3.cm_time * dtau;
    compute_derivatives(tmp, full_system, k4);

    const double w = dtau / 6.0;
    state.Q1 += w * (k1.Q1 + 2.0 * k2.Q1 + 2.0 * k3.Q1 + k4.Q1);
    state.P1 += w * (k1.P1 + 2.0 * k2.P1 + 2.0 * k3.P1 + k4.P1);
    state.Q2 += w * (k1.Q2 + 2.0 * k2.Q2 + 2.0 * k3.Q2 + k4.Q2);
    state.P2 += w * (k1.P2 + 2.0 * k2.P2 + 2.0 * k3.P2 + k4.P2);

    const double dt_phys = w * (k1.cm_time + 2.0 * k2.cm_time + 2.0 * k3.cm_time + k4.cm_time);
    state.cm_time += dt_phys;
    state.cm_pos += state.cm_vel * dt_phys;
    state.tau += dtau;

    if (!std::isfinite(state.Q1.norm()) || !std::isfinite(state.Q2.norm()))
        throw std::runtime_error("rk4_step: NaN/Inf detectado");
}

// ============================================================================
// ESTIMACIÓN DE ERROR (paso doble)
// ============================================================================
double Chain3Integrator::estimate_error(const Chain3State& s1, const Chain3State& s2) {
    return ((s2.Q1 - s1.Q1).norm() + (s2.P1 - s1.P1).norm() +
            (s2.Q2 - s1.Q2).norm() + (s2.P2 - s1.P2).norm()) / 4.0;
}

// ============================================================================
// INTEGRACIÓN ADAPTATIVA
// ============================================================================
void Chain3Integrator::integrate(Chain3State& state,
                                 double t_target,
                                 double& t_achieved,
                                 const IntegrationParams& params,
                                 const NBodySystem& full_system) {
    double dtau = internal_dt;

    while (state.cm_time < t_target - 1e-12) {
        // Paso de prueba con dtau
        Chain3State trial = state;
        bool step_ok = true;
        try {
            rk4_step(trial, dtau, full_system);
        } catch (const std::runtime_error& e) {
            // Capturar excepción de singularidad y reducir paso
            step_ok = false;
        } catch (...) {
            step_ok = false;
        }

        if (!step_ok) {
            dtau = std::max(dtau * 0.5, params.min_dtau);
            if (dtau == params.min_dtau) {
                // Forzar avance mínimo
                try { rk4_step(state, dtau, full_system); } catch (...) {}
                break;
            }
            continue;
        }

        // Paso doble con dtau/2 para estimar error
        Chain3State half = state;
        try {
            rk4_step(half, dtau * 0.5, full_system);
            rk4_step(half, dtau * 0.5, full_system);
        } catch (...) {
            dtau = std::max(dtau * 0.5, params.min_dtau);
            continue;
        }

        double err = estimate_error(trial, half);
        double tol = params.abs_tol;

        if (err <= tol || dtau <= params.min_dtau) {
            // Aceptar (usamos half, más preciso)
            state = half;
            if (err > 0.0) {
                double factor = 0.9 * std::pow(tol / err, 0.2);
                dtau = std::clamp(dtau * factor, params.min_dtau, params.max_dtau);
            } else {
                dtau = std::min(dtau * 2.0, params.max_dtau);
            }
        } else {
            double factor = 0.9 * std::pow(tol / err, 0.25);
            dtau = std::max(dtau * factor, params.min_dtau);
        }
    }

    t_achieved = state.cm_time;
}

// ============================================================================
// MÉTODOS AUXILIARES
// ============================================================================
void Chain3Integrator::rk4_step_with_error(Chain3State& state, double dtau,
                                           const NBodySystem& fs,
                                           Chain3State& err) {
    Chain3State y0 = state;
    rk4_step(state, dtau, fs);
    err.Q1 = state.Q1 - y0.Q1;
    err.P1 = state.P1 - y0.P1;
    err.Q2 = state.Q2 - y0.Q2;
    err.P2 = state.P2 - y0.P2;
}

void Chain3Integrator::update_time(Chain3State&, double) {}

void Chain3Integrator::integrate_simple(Chain3State& s, double dt,
                                        const NBodySystem& fs) {
    double t_target = s.cm_time + dt;
    double t_achieved;
    IntegrationParams p;
    p.abs_tol = 1e-8;
    p.min_dtau = 1e-8;
    p.max_dtau = 1e-1;
    integrate(s, t_target, t_achieved, p, fs);
}