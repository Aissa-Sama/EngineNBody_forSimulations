// regularization/chain/archain_n_integrator.cpp
#include "archain_n_integrator.h"
#include "nbody_system.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <limits>

// ============================================================================
// CONSTRUCTOR
// ============================================================================
ARChainNIntegrator::ARChainNIntegrator(double eta)
    : eta_(eta), ds_min_(1e-16), ds_max_(1.0)
{}

// ============================================================================
// FIND_CHAIN_INDICES — algoritmo greedy de Mikkola & Aarseth (1993)
//
// Construye una permutación chain[0..N-1] tal que los cuerpos más cercanos
// son vecinos directos en la cadena. Esto minimiza el error de redondeo
// en las separaciones críticas (eslabones directos = diferencias exactas).
//
// ALGORITMO:
//   1. Encontrar el par (i*,j*) con sep mínima → primer eslabón.
//   2. Marcar ambos como 'en_cadena'. chain = [i*, j*].
//   3. Mientras queden cuerpos libres:
//      - Calcular dist(extremo_izq, libre_k) para todo k libre.
//      - Calcular dist(extremo_der, libre_k) para todo k libre.
//      - Escoger la extensión (izq o der) que minimiza la distancia.
//      - Prepend (izq) o append (der) el cuerpo seleccionado.
// ============================================================================
void ARChainNIntegrator::find_chain_indices(const std::vector<Vec3>& r,
                                             std::vector<int>& chain)
{
    const int N = static_cast<int>(r.size());
    chain.clear();
    chain.reserve(N);

    std::vector<bool> in_chain(N, false);

    // ── Paso 1: par más cercano ──────────────────────────────────────────────
    double sep_min = std::numeric_limits<double>::max();
    int best_i = 0, best_j = 1;
    for (int i = 0; i < N; ++i)
        for (int j = i+1; j < N; ++j) {
            double d = (r[i] - r[j]).norm();
            if (d < sep_min) { sep_min = d; best_i = i; best_j = j; }
        }

    chain.push_back(best_i);
    chain.push_back(best_j);
    in_chain[best_i] = true;
    in_chain[best_j] = true;

    // ── Pasos siguientes: extensión greedy ──────────────────────────────────
    for (int step = 2; step < N; ++step) {
        const Vec3& left  = r[chain.front()];
        const Vec3& right = r[chain.back()];

        double best_dist_L = std::numeric_limits<double>::max();
        double best_dist_R = std::numeric_limits<double>::max();
        int    best_k_L = -1, best_k_R = -1;

        for (int k = 0; k < N; ++k) {
            if (in_chain[k]) continue;
            double dL = (r[k] - left).norm();
            double dR = (r[k] - right).norm();
            if (dL < best_dist_L) { best_dist_L = dL; best_k_L = k; }
            if (dR < best_dist_R) { best_dist_R = dR; best_k_R = k; }
        }

        if (best_dist_L <= best_dist_R) {
            chain.insert(chain.begin(), best_k_L);
            in_chain[best_k_L] = true;
        } else {
            chain.push_back(best_k_R);
            in_chain[best_k_R] = true;
        }
    }
}

// ============================================================================
// REBUILD_CHAIN — re-construye X[], W[] desde posiciones absolutas
// ============================================================================
void ARChainNIntegrator::rebuild_chain(ARChainNState& state,
                                        const std::vector<Vec3>& r_abs,
                                        const std::vector<Vec3>& v_abs)
{
    find_chain_indices(r_abs, state.chain);
    for (int k = 0; k < state.N - 1; ++k) {
        state.X[k] = r_abs[state.chain[k+1]] - r_abs[state.chain[k]];
        state.W[k] = v_abs[state.chain[k+1]] - v_abs[state.chain[k]];
    }
}

// ============================================================================
// INITIALIZE
// ============================================================================
ARChainNState ARChainNIntegrator::initialize(const NBodySystem& system,
                                              const std::vector<int>& indices) const
{
    const int N = static_cast<int>(indices.size());
    if (N < 2)
        throw std::invalid_argument("ARChainNIntegrator::initialize: N < 2");
    for (int idx : indices)
        if (idx < 0 || idx >= (int)system.bodies.size())
            throw std::out_of_range("ARChainNIntegrator::initialize: índice fuera de rango");

    ARChainNState st(N);

    // Masas y CM
    double M = 0.0;
    Vec3 cm_pos, cm_vel;
    for (int i = 0; i < N; ++i) {
        const auto& b = system.bodies[indices[i]];
        st.masses[i] = b.mass;
        M += b.mass;
    }
    for (int i = 0; i < N; ++i) {
        const auto& b = system.bodies[indices[i]];
        cm_pos = cm_pos + b.position * (b.mass / M);
        cm_vel = cm_vel + b.velocity * (b.mass / M);
    }
    st.cm_pos = cm_pos;
    st.cm_vel = cm_vel;
    st.t_phys = 0.0;
    st.s_fict = 0.0;

    // Posiciones y velocidades absolutas (en orden de indices[])
    std::vector<Vec3> r_abs(N), v_abs(N);
    for (int i = 0; i < N; ++i) {
        r_abs[i] = system.bodies[indices[i]].position;
        v_abs[i] = system.bodies[indices[i]].velocity;
    }

    // Construir cadena óptima y calcular X[], W[]
    rebuild_chain(st, r_abs, v_abs);

    st.energy = compute_energy(st);

    return st;
}

// ============================================================================
// WRITE_BACK
// ============================================================================
void ARChainNIntegrator::write_back(const ARChainNState& state,
                                     NBodySystem& system,
                                     const std::vector<int>& indices) const
{
    if (!state.is_valid())
        throw std::runtime_error("ARChainNIntegrator::write_back: estado inválido");

    std::vector<Vec3> r_abs, v_abs;
    state.reconstruct_positions(r_abs);
    state.reconstruct_velocities(v_abs);

    const Vec3 cm_disp = state.cm_vel * state.t_phys;
    for (int i = 0; i < state.N; ++i) {
        system.bodies[indices[i]].position = r_abs[i] + cm_disp;
        system.bodies[indices[i]].velocity = v_abs[i];
    }
}

// ============================================================================
// COMPUTE_OMEGA — Ω = Σᵢ<ⱼ mᵢmⱼ / rᵢⱼ para todos los N(N-1)/2 pares
//
// Las separaciones entre cuerpos NO adyacentes en la cadena se calculan
// como suma de eslabones intermedios: R(i,j) = X[i] + ... + X[j-1].
// Esto evita cancelación catastrófica para los pares más cercanos
// (que son vecinos directos en la cadena, i.e., j = i+1).
// ============================================================================
double ARChainNIntegrator::compute_Omega(const ARChainNState& st) const
{
    const int N = st.N;
    double omega = 0.0;

    for (int i = 0; i < N-1; ++i) {
        Vec3 rij;
        for (int j = i+1; j < N; ++j) {
            rij = rij + st.X[j-1];                // R(i,j) acumulado
            double r = rij.norm();
            if (r < 1e-30)
                throw std::runtime_error("compute_Omega: separación singular");
            omega += st.masses[st.chain[i]] * st.masses[st.chain[j]] / r;
        }
    }
    return omega;
}

// ============================================================================
// COMPUTE_SEP_MIN
// ============================================================================
double ARChainNIntegrator::compute_sep_min(const ARChainNState& st) const
{
    double sep_min = std::numeric_limits<double>::max();
    // Eslabones directos (máxima precisión)
    for (int k = 0; k < st.N - 1; ++k) {
        double d = st.X[k].norm();
        if (d < sep_min) sep_min = d;
    }
    // Pares no adyacentes (pueden ser más cercanos si la cadena no es óptima)
    for (int i = 0; i < st.N - 2; ++i) {
        Vec3 rij;
        for (int j = i + 2; j < st.N; ++j) {
            rij = rij + st.X[j-1];
            double d = rij.norm();
            if (d < sep_min) sep_min = d;
        }
    }
    return sep_min;
}

// ============================================================================
// COMPUTE_ENERGY — T + U
// ============================================================================
double ARChainNIntegrator::compute_energy(const ARChainNState& st) const
{
    std::vector<Vec3> r_abs, v_abs;
    st.reconstruct_positions(r_abs);
    st.reconstruct_velocities(v_abs);

    double T = 0.0;
    for (int i = 0; i < st.N; ++i)
        T += 0.5 * st.masses[i] * dot(v_abs[i], v_abs[i]);

    double U = 0.0;
    for (int i = 0; i < st.N - 1; ++i) {
        Vec3 rij;
        for (int j = i + 1; j < st.N; ++j) {
            rij = rij + st.X[j-1];
            // NOTA: los índices en X[] corresponden a la cadena, pero para
            // la energía necesitamos las masas en el orden de la cadena.
            // chain[i] y chain[j] son los índices globales (dentro de indices[]).
            U -= st.masses[st.chain[i]] * st.masses[st.chain[j]] / rij.norm();
        }
    }

    return T + U;
}

// ============================================================================
// COMPUTE_ACCELERATIONS — diferencias de aceleración en coordenadas de cadena
//
// A[k] = a[chain[k+1]] - a[chain[k]]   para k = 0..N-2
//
// La aceleración absoluta del i-ésimo cuerpo en la cadena es:
//   a[chain[i]] = Σ_{j≠i} m[chain[j]] * (r[chain[j]] - r[chain[i]]) / r³ᵢⱼ
//
// Para calcular A[k] = a[chain[k+1]] - a[chain[k]], se puede reformular
// directamente usando las separaciones de cadena R(i,j) para evitar
// reconstruir posiciones absolutas.
//
// IMPLEMENTACIÓN:
//   1. Precalcular todas las separaciones R(i,j) = chain_sep(i,j) para i<j.
//   2. Calcular aceleraciones absolutas de cada cuerpo.
//   3. A[k] = a[k+1] - a[k].
//
// Nota: para N pequeño (≤ 20), la reconstrucción de posiciones absolutas
// y el cálculo directo es equivalente en precisión a usar chain_sep(),
// porque los pares críticos (vecinos directos) ya están bien representados
// por X[k]. Se usa reconstruct_positions() para simplicidad y correctitud.
// ============================================================================
void ARChainNIntegrator::compute_accelerations(const ARChainNState& st,
                                                const std::vector<Vec3>& r_abs,
                                                std::vector<Vec3>& A_out) const
{
    const int N = st.N;

    // Aceleraciones absolutas de cada cuerpo en orden de cadena
    std::vector<Vec3> a(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            Vec3 dr = r_abs[st.chain[j]] - r_abs[st.chain[i]];
            double r2 = dot(dr, dr);
            double r  = std::sqrt(r2);
            if (r < 1e-30)
                throw std::runtime_error("compute_accelerations: separación singular");
            a[i] = a[i] + dr * (st.masses[st.chain[j]] / (r2 * r));
        }
    }

    // Diferencias de aceleración en coordenadas de cadena
    A_out.resize(N - 1);
    for (int k = 0; k < N - 1; ++k)
        A_out[k] = a[k+1] - a[k];
}

// ============================================================================
// LEAPFROG_STEP_BARE — un paso DKD en tiempo ficticio. Acepta ds < 0.
//
// Mismo esquema DKD que ARChain3Integrator::leapfrog_step_bare,
// generalizado a N eslabones:
//
//   DRIFT½:   X[k] += W[k] * (ds/2 * Ω_ini)      para k = 0..N-2
//   KICK:     W[k] += A[k] * (ds   * Ω_mid)
//   DRIFT½:   X[k] += W[k] * (ds/2 * Ω_new)
//   TIME:     t    += ds * Ω_mid
// ============================================================================
void ARChainNIntegrator::leapfrog_step_bare(ARChainNState& st, double ds) const
{
    const int N = st.N;

    // ── DRIFT½ ────────────────────────────────────────────────────────────────
    const double Omega_ini = compute_Omega(st);
    const double half_ds   = 0.5 * ds;
    const double drift_ini = half_ds * Omega_ini;

    for (int k = 0; k < N-1; ++k)
        st.X[k] = st.X[k] + st.W[k] * drift_ini;

    // ── KICK ──────────────────────────────────────────────────────────────────
    const double Omega_mid = compute_Omega(st);
    const double dt_mid    = ds * Omega_mid;

    // Reconstruir posiciones para compute_accelerations
    std::vector<Vec3> r_abs;
    st.reconstruct_positions(r_abs);

    std::vector<Vec3> A;
    compute_accelerations(st, r_abs, A);

    for (int k = 0; k < N-1; ++k)
        st.W[k] = st.W[k] + A[k] * dt_mid;

    // ── DRIFT½ ────────────────────────────────────────────────────────────────
    const double Omega_new  = compute_Omega(st);
    const double drift_new  = half_ds * Omega_new;

    for (int k = 0; k < N-1; ++k)
        st.X[k] = st.X[k] + st.W[k] * drift_new;

    // ── AVANCE DE TIEMPO ──────────────────────────────────────────────────────
    st.t_phys += dt_mid;
    st.s_fict += ds;
}

void ARChainNIntegrator::leapfrog_step(ARChainNState& st, double ds) const
{
    leapfrog_step_bare(st, ds);
}

// ============================================================================
// INTEGRATE_UNTIL — núcleo del bucle de integración
//
// Idéntico al de ARChain3Integrator, pero con re-construcción de cadena
// cuando la configuración cambia significativamente.
// ============================================================================
void ARChainNIntegrator::integrate_until(ARChainNState& state,
                                          double t_target) const
{
    if (t_target <= state.t_phys) return;

    constexpr int max_steps = 10'000'000;
    int steps_remaining = max_steps;

    while (state.t_phys < t_target && steps_remaining-- > 0) {

        // ── Re-construcción de cadena si es necesario ──────────────────────
        {
            std::vector<Vec3> r_abs;
            state.reconstruct_positions(r_abs);
            if (state.needs_rechain(r_abs)) {
                std::vector<Vec3> v_abs;
                state.reconstruct_velocities(v_abs);
                rebuild_chain(state, r_abs, v_abs);
            }
        }

        const double Omega   = compute_Omega(state);
        const double sep_min = compute_sep_min(state);

        double ds_step = eta_ * sep_min / Omega;
        ds_step = std::clamp(ds_step, ds_min_, ds_max_);

        leapfrog_step(state, ds_step);
    }

    if (steps_remaining <= 0)
        std::cerr << "[AR-N] MAX_STEPS en t=" << state.t_phys
                  << " target=" << t_target
                  << " sep=" << compute_sep_min(state) << "\n";
}

// ============================================================================
// INTEGRATE / INTEGRATE_TO
// ============================================================================
void ARChainNIntegrator::integrate(ARChainNState& state, double dt_phys) const
{
    if (dt_phys <= 0.0) return;
    integrate_until(state, state.t_phys + dt_phys);
}

void ARChainNIntegrator::integrate_to(ARChainNState& state,
                                       double t_abs_final) const
{
    if (t_abs_final <= state.t_phys) return;
    integrate_until(state, t_abs_final);
}