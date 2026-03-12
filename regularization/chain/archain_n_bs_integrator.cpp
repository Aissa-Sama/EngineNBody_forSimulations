// regularization/chain/archain_n_bs_integrator.cpp
#include "archain_n_bs_integrator.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// ============================================================================
// CONSTRUCTOR
// ============================================================================
ARChainNBSIntegrator::ARChainNBSIntegrator(double eta,
                                           const BSParameters& params)
    : ARChainNIntegrator(eta)
    , params_(params)
{}

// ============================================================================
// MODIFIED MIDPOINT STEP (Gragg 1965) — generalizado a N cuerpos
//
// Idéntico al de ARChain3BSIntegrator pero opera sobre vectores X[], W[]
// de longitud N-1 en lugar de (X1,X2,V1,V2) individuales.
//
// Ecuaciones de movimiento en tiempo ficticio τ:
//   dX[k]/dτ = Ω(X) · W[k]
//   dW[k]/dτ = Ω(X) · A[k](X)
//   dt_phys/dτ = Ω(X)
// ============================================================================
void ARChainNBSIntegrator::modified_midpoint_step(
    const ARChainNState& y0,
    double               ds_total,
    int                  n,
    ARChainNState&       result) const
{
    const int N  = y0.N;
    const double h = ds_total / static_cast<double>(n);

    ARChainNState z_prev = y0;
    ARChainNState z_curr = y0;

    // ── Primer paso: Euler ─────────────────────────────────────────────────
    const double Omega0 = compute_Omega(z_prev);
    std::vector<Vec3> r_abs0;
    z_prev.reconstruct_positions(r_abs0);
    std::vector<Vec3> A0;
    compute_accelerations(z_prev, r_abs0, A0);

    const double dt0 = h * Omega0;
    for (int k = 0; k < N-1; ++k) {
        z_curr.X[k]  = z_prev.X[k] + z_prev.W[k] * dt0;
        z_curr.W[k]  = z_prev.W[k] + A0[k]        * dt0;
    }
    z_curr.t_phys = z_prev.t_phys + dt0;
    z_curr.s_fict = z_prev.s_fict + h;

    // ── Pasos medios: m = 1 .. n-1 ────────────────────────────────────────
    for (int m = 1; m < n; ++m) {
        const double Omega_m = compute_Omega(z_curr);
        std::vector<Vec3> r_abs_m;
        z_curr.reconstruct_positions(r_abs_m);
        std::vector<Vec3> A_m;
        compute_accelerations(z_curr, r_abs_m, A_m);

        const double dt_m = 2.0 * h * Omega_m;

        ARChainNState z_next = z_curr;
        for (int k = 0; k < N-1; ++k) {
            z_next.X[k] = z_prev.X[k] + z_curr.W[k] * dt_m;
            z_next.W[k] = z_prev.W[k] + A_m[k]       * dt_m;
        }
        z_next.t_phys = z_prev.t_phys + dt_m;
        z_next.s_fict = z_curr.s_fict + h;

        z_prev = z_curr;
        z_curr = z_next;
    }

    // ── Paso final suavizante ──────────────────────────────────────────────
    const double Omega_n = compute_Omega(z_curr);
    std::vector<Vec3> r_abs_n;
    z_curr.reconstruct_positions(r_abs_n);
    std::vector<Vec3> A_n;
    compute_accelerations(z_curr, r_abs_n, A_n);

    const double dt_n = h * Omega_n;

    result = y0;
    for (int k = 0; k < N-1; ++k) {
        result.X[k] = 0.5 * (z_curr.X[k] + z_prev.X[k] + z_curr.W[k] * dt_n);
        result.W[k] = 0.5 * (z_curr.W[k] + z_prev.W[k] + A_n[k]       * dt_n);
    }
    result.t_phys = 0.5 * (z_curr.t_phys + z_prev.t_phys + dt_n);
    result.s_fict = y0.s_fict + ds_total;
}

// ============================================================================
// EXTRAPOLACIÓN POLINÓMICA DE RICHARDSON (Neville-Aitken)
//
// Idéntica a ARChain3BSIntegrator pero opera sobre vectores X[], W[].
// La tabla T[k] es un ARChainNState completo.
// ============================================================================
void ARChainNBSIntegrator::polynomial_extrapolation(
    const std::vector<ARChainNState>& raw,
    const std::vector<double>&        step_sizes,
    int                               k_current,
    ARChainNState&                    extrap_out) const
{
    const int N = raw[0].N;
    std::vector<ARChainNState> T(raw.begin(), raw.begin() + k_current + 1);

    for (int m = 1; m <= k_current; ++m) {
        for (int k = k_current; k >= m; --k) {
            const double ratio  = step_sizes[k - m] / step_sizes[k];
            double       factor = ratio * ratio - 1.0;
            if (std::abs(factor) < 1e-15) factor = 1e-15;
            const double inv_f = 1.0 / factor;

            for (int j = 0; j < N-1; ++j) {
                T[k].X[j] = T[k].X[j] + (T[k].X[j] - T[k-1].X[j]) * inv_f;
                T[k].W[j] = T[k].W[j] + (T[k].W[j] - T[k-1].W[j]) * inv_f;
            }
            T[k].t_phys = T[k].t_phys
                        + (T[k].t_phys - T[k-1].t_phys) * inv_f;
        }
    }

    extrap_out = T[k_current];
}

// ============================================================================
// ESTIMATE_ERROR — norma escalada de diferencias entre dos extrapolaciones
// ============================================================================
double ARChainNBSIntegrator::estimate_error(
    const ARChainNState& hi,
    const ARChainNState& lo) const
{
    const int N = hi.N;
    double err = 0.0;
    int    cnt = 0;

    for (int k = 0; k < N-1; ++k) {
        const double sX = std::max(hi.X[k].norm(), 1e-10);
        const double sW = std::max(hi.W[k].norm(), 1e-10);
        err += (hi.X[k] - lo.X[k]).norm() / sX;
        err += (hi.W[k] - lo.W[k]).norm() / sW;
        cnt += 2;
    }
    return (cnt > 0) ? err / cnt : 0.0;
}

// ============================================================================
// STEP_SCALE_FACTOR
// ============================================================================
double ARChainNBSIntegrator::step_scale_factor(double error, int k_opt) const
{
    if (error < 1e-30) return 4.0;
    const double exponent = 1.0 / (2.0 * k_opt + 1.0);
    const double scale    = params_.safety
                          * std::pow(params_.bs_eps / error, exponent);
    return std::clamp(scale, 0.1, 4.0);
}

// ============================================================================
// ENERGY_ERROR
// ============================================================================
double ARChainNBSIntegrator::energy_error(const ARChainNState& state,
                                           double E_ref) const
{
    return std::abs(compute_energy(state) - E_ref)
         / (std::abs(E_ref) + 1e-30);
}

// ============================================================================
// BS_STEP
// ============================================================================
bool ARChainNBSIntegrator::bs_step(
    ARChainNState& state,
    double         ds_total,
    double&        error_out,
    double&        ds_next_out)
{
    const ARChainNState y_save = state;

    std::vector<ARChainNState> raw;
    std::vector<double>        step_sizes;
    raw.reserve(K_MAX);
    step_sizes.reserve(K_MAX);

    ARChainNState extrap_prev;
    bool  has_prev  = false;
    error_out       = 1e30;

    for (int k = 0; k < params_.k_max && k < K_MAX; ++k) {
        const int    n    = n_seq_[k];
        const double ds_k = ds_total / static_cast<double>(n);

        ARChainNState r(y_save.N);
        try {
            modified_midpoint_step(y_save, ds_total, n, r);
        } catch (const std::exception&) {
            ds_next_out = ds_total * 0.25;
            return false;
        }
        raw.push_back(r);
        step_sizes.push_back(ds_k);

        if (k >= 1) {
            ARChainNState extrap_curr(y_save.N);
            polynomial_extrapolation(raw, step_sizes, k, extrap_curr);
            // Preservar metadata que no se extrapola
            extrap_curr.chain  = y_save.chain;
            extrap_curr.masses = y_save.masses;
            extrap_curr.cm_pos = y_save.cm_pos;
            extrap_curr.cm_vel = y_save.cm_vel;
            extrap_curr.energy = y_save.energy;
            extrap_curr.N      = y_save.N;

            if (has_prev) {
                error_out = estimate_error(extrap_curr, extrap_prev);
                const double scaled = error_out / params_.bs_eps;

                if (scaled <= 1.0) {
                    state       = extrap_curr;
                    ds_next_out = ds_total * step_scale_factor(error_out, k);
                    return true;
                }
                if (scaled > 1000.0 && k >= 3) break;
            }
            extrap_prev = extrap_curr;
            has_prev    = true;
        }
    }

    if (has_prev) {
        state       = extrap_prev;
        ds_next_out = ds_total * 0.5;
        return true;
    }

    ds_next_out = ds_total * 0.5;
    return false;
}

// ============================================================================
// INTEGRATE_TO_BS — bucle principal
// ============================================================================
void ARChainNBSIntegrator::integrate_to_bs(ARChainNState& state,
                                            double         t_abs_final)
{
    if (t_abs_final <= state.t_phys) return;

    const double E0      = state.energy;
    double       ds      = params_.initial_ds;
    int          n_steps = 0, n_reject = 0;
    int          n_consecutive_reject = 0;

    while (state.t_phys < t_abs_final - 1e-14) {

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

        const double t_left = t_abs_final - state.t_phys;
        const double Omega   = compute_Omega(state);
        const double sep_min = compute_sep_min(state);

        double ds_eta  = eta_ * sep_min / Omega;
        ds_eta = std::clamp(ds_eta, params_.min_ds, params_.max_ds);

        const double ds_max_t = t_left / (Omega * 0.5 + 1e-30);
        ds = std::min({ds, ds_eta, ds_max_t});
        ds = std::clamp(ds, params_.min_ds, params_.max_ds);

        ARChainNState trial = state;
        double error = 0.0, ds_next = ds;

        bool accepted = false;
        try {
            accepted = bs_step(trial, ds, error, ds_next);
            n_consecutive_reject = 0;
        } catch (const std::runtime_error&) {
            ds *= 0.25;
            ds = std::max(ds, params_.min_ds);
            ++n_reject; ++n_consecutive_reject;
            if (n_consecutive_reject > 30)
                throw std::runtime_error("ARChainNBS: singularidad persistente");
            continue;
        }

        if (!accepted) {
            ds = std::clamp(ds_next, params_.min_ds, params_.max_ds);
            ++n_reject; ++n_consecutive_reject;
            if (n_consecutive_reject > 30)
                throw std::runtime_error("ARChainNBS: no convergencia persistente");
            continue;
        }
        n_consecutive_reject = 0;

        if (trial.t_phys > t_abs_final + 1e-12) {
            ds *= 0.5;
            ds = std::max(ds, params_.min_ds);
            ++n_reject;
            continue;
        }

        const double dE = energy_error(trial, E0);
        if (dE > params_.energy_tol
            && ds > params_.min_ds * 4.0
            && n_reject < params_.max_steps / 2) {
            ds = std::max(ds * 0.5, params_.min_ds);
            ++n_reject;
            continue;
        }

        state = trial;
        ds    = std::clamp(ds_next, params_.min_ds, params_.max_ds);
        ++n_steps;

        if (n_steps > params_.max_steps)
            throw std::runtime_error("ARChainNBS: máximo de pasos excedido");
    }

    if (params_.verbose)
        std::cout << "[AR-N-BS] " << n_steps << " pasos, "
                  << n_reject << " rechazados"
                  << "  t_final=" << state.t_phys << "\n";
}