// regularization/chain/archain_n_pn_integrator.cpp
#include "archain_n_pn_integrator.h"
#include <cmath>
#include <stdexcept>

// ============================================================================
// CONSTRUCTOR
// ============================================================================
ARChainNPNIntegrator::ARChainNPNIntegrator(double eta,
                                           double c_speed,
                                           int    pn_order)
    : ARChainNIntegrator(eta)
    , c_(c_speed)
    , pn_order_(pn_order)
{
    if (c_ <= 0.0)
        throw std::invalid_argument("ARChainNPNIntegrator: c_speed debe ser > 0");
}

// ============================================================================
// COMPUTE_ACCELERATIONS — override con correcciones PN
//
// Las aceleraciones absolutas de cada cuerpo son:
//   a_total[i] = a_N[i] + a_PN[i]
//
// Las diferencias de cadena son:
//   A[k] = a_total[chain[k+1]] - a_total[chain[k]]
//
// IMPLEMENTACIÓN:
//   1. Reconstruir velocidades absolutas (necesario para PN)
//   2. Calcular a_N[i] para cada cuerpo (igual que la base)
//   3. Calcular a_PN[i] = suma de correcciones PN
//   4. A[k] = (a_N + a_PN)[chain[k+1]] - (a_N + a_PN)[chain[k]]
//
// NOTA sobre r_abs:
//   r_abs ya está reconstruido en orden de CADENA: r_abs[chain[k]] es la
//   posición del k-ésimo cuerpo de la cadena. Pero la convención en
//   ARChainNIntegrator es que r_abs[i] = posición del cuerpo con índice i
//   en el sistema local (0..N-1 del AR-chain, no del sistema global).
//   Las aceleraciones absolutas a[i] usan el mismo índice local.
// ============================================================================
void ARChainNPNIntegrator::compute_accelerations(
    const ARChainNState& state,
    const std::vector<Vec3>& r_abs,
    std::vector<Vec3>& A_out) const
{
    const int N = state.N;

    // ── Velocidades absolutas ──────────────────────────────────────────────
    // Necesarias para las correcciones PN (dependen de velocidades)
    std::vector<Vec3> v_abs;
    state.reconstruct_velocities(v_abs);

    // ── Masas en orden de índice LOCAL (0..N-1) ───────────────────────────
    // state.masses[i] = masa del cuerpo con índice local i
    // r_abs[chain[k]] = posición del k-ésimo eslabón de la cadena

    // Construir r[], v[], m[] en el mismo orden (índice local 0..N-1)
    // r_abs ya está en ese orden (r_abs[i] = posición del cuerpo i)
    // v_abs está en orden LOCAL también

    // ── Aceleraciones Newtonianas absolutas ───────────────────────────────
    std::vector<Vec3> a_newton(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            Vec3 dr = r_abs[j] - r_abs[i];
            double r2 = dot(dr, dr);
            double r  = std::sqrt(r2);
            if (r < 1e-30)
                throw std::runtime_error("ARChainNPN: separación singular (Newton)");
            a_newton[i] = a_newton[i] + dr * (state.masses[j] / (r2 * r));
        }
    }

    // ── Correcciones post-Newtonianas ─────────────────────────────────────
    std::vector<Vec3> a_pn(N);
    PNForces::compute_pn_corrections(r_abs, v_abs, state.masses,
                                     c_, pn_order_, a_pn);

    // ── Diferencias de aceleración en coordenadas de cadena ───────────────
    // A[k] = a_total[chain[k+1]] - a_total[chain[k]]
    // donde los índices chain[k] son locales (0..N-1)
    A_out.resize(N - 1);
    for (int k = 0; k < N - 1; ++k) {
        const int ci  = state.chain[k];
        const int ci1 = state.chain[k+1];
        A_out[k] = (a_newton[ci1] + a_pn[ci1]) - (a_newton[ci] + a_pn[ci]);
    }
}

// ============================================================================
// COMPUTE_ENERGY_PN1 — energía conservativa incluyendo corrección PN1
//
// H_PN1 = T_N + T_PN1 + U_N + U_PN1
//
// donde (en gauge armónica):
//   T_PN1 = −(3/8c²) × Σᵢ mᵢ v⁴ᵢ − (1/2c²) × Σᵢ<ⱼ mᵢmⱼ/rᵢⱼ × (vᵢ−vⱼ)²
//           (corrección cinética)
//   U_PN1 = −(1/2c²) × Σᵢ<ⱼ mᵢmⱼ/rᵢⱼ × (4vᵢ·vⱼ − 3v²ᵢ − 3v²ⱼ + ...)
//           (corrección potencial — incluye velocidades)
//
// Referencia: Blanchet & Iyer (2003) ec. (145); Soffel (1989) §3.4.
//
// NOTA: la energía PN1 no es un invariante exacto del leapfrog — solo lo
// es del MMP + GBS. Se proporciona como diagnóstico.
// ============================================================================
double ARChainNPNIntegrator::compute_energy_pn1(const ARChainNState& state) const
{
    std::vector<Vec3> r_abs, v_abs;
    state.reconstruct_positions(r_abs);
    state.reconstruct_velocities(v_abs);

    const int N = state.N;
    const double inv_c2 = 1.0 / (c_ * c_);

    // ── Energía Newtoniana ─────────────────────────────────────────────────
    double T_N = 0.0, U_N = 0.0;
    for (int i = 0; i < N; ++i)
        T_N += 0.5 * state.masses[i] * dot(v_abs[i], v_abs[i]);
    for (int i = 0; i < N; ++i)
        for (int j = i+1; j < N; ++j) {
            double rij = (r_abs[i] - r_abs[j]).norm();
            U_N -= state.masses[i] * state.masses[j] / rij;
        }

    // ── Corrección cinética PN1: −(3/8c²) × Σᵢ mᵢ v⁴ᵢ ───────────────────
    double T_PN1 = 0.0;
    for (int i = 0; i < N; ++i) {
        double v2 = dot(v_abs[i], v_abs[i]);
        T_PN1 -= (3.0/8.0) * state.masses[i] * v2 * v2 * inv_c2;
    }

    // ── Corrección potencial PN1 ──────────────────────────────────────────
    double U_PN1 = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
            Vec3   dr   = r_abs[i] - r_abs[j];
            double rij  = dr.norm();
            double v2i  = dot(v_abs[i], v_abs[i]);
            double v2j  = dot(v_abs[j], v_abs[j]);
            double vivj = dot(v_abs[i], v_abs[j]);
            Vec3   nij  = dr * (1.0/rij);
            double vi_n = dot(v_abs[i], nij);
            double vj_n = dot(v_abs[j], nij);

            // Forma de Blanchet & Iyer (2003)
            U_PN1 += state.masses[i] * state.masses[j] / rij * inv_c2
                   * (-1.5*v2i - 1.5*v2j + 3.5*vivj
                      - 0.5*vi_n*vj_n
                      - 0.5*state.masses[i]/rij - 0.5*state.masses[j]/rij);
        }
    }

    return T_N + U_N + T_PN1 + U_PN1;
}
