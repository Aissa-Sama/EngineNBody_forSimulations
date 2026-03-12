// regularization/chain/archain_n_bs_integrator.h
#pragma once
#include "archain_n_integrator.h"
#include <array>
#include <vector>

/**
 * @brief Integrador Bulirsch-Stoer para AR-chain de N cuerpos arbitrario.
 *
 * Exacta generalización de ARChain3BSIntegrator a N cuerpos.
 * El esquema GBS es idéntico — solo cambia el integrador base (ARChainN
 * en lugar de ARChain3).
 *
 * La extrapolación opera sobre el vector de estado generalizado:
 *   - X[0..N-2], W[0..N-2] (coordenadas de cadena)
 *   - t_phys
 *
 * REFERENCIAS:
 *   Gragg (1965), Bulirsch & Stoer (1966) — método GBS.
 *   Mikkola & Aarseth (2002), Cel. Mech. 84, 343 — GBS sobre TTL.
 *   Mikkola & Merritt (2008), AJ 135, 2398 — AR-CHAIN N cuerpos.
 */
class ARChainNBSIntegrator : public ARChainNIntegrator {
public:

    struct BSParameters {
        double bs_eps      = 1e-10;
        int    k_max       = 8;
        double initial_ds  = 1e-3;
        double min_ds      = 1e-14;
        double max_ds      = 0.5;
        double safety      = 0.9;
        double energy_tol  = 0.01;
        int    max_steps   = 1000000;
        bool   verbose     = false;
    };

    explicit ARChainNBSIntegrator(double eta, const BSParameters& params);
    explicit ARChainNBSIntegrator(double eta = 1e-3)
        : ARChainNBSIntegrator(eta, BSParameters{}) {}

    // ── Interfaz principal ───────────────────────────────────────────────────

    void integrate_to_bs(ARChainNState& state, double t_abs_final);

    bool bs_step(ARChainNState& state,
                 double         ds_total,
                 double&        error_out,
                 double&        ds_next_out);

    double energy_error(const ARChainNState& state, double E_ref) const;

private:
    BSParameters params_;

    static constexpr int K_MAX = 8;
    static constexpr std::array<int, K_MAX> n_seq_ = {2, 4, 6, 8, 10, 12, 14, 16};

    void modified_midpoint_step(const ARChainNState& y0,
                                double               ds_total,
                                int                  n,
                                ARChainNState&       result) const;

    void polynomial_extrapolation(const std::vector<ARChainNState>& raw,
                                  const std::vector<double>&        step_sizes,
                                  int                               k_current,
                                  ARChainNState&                    extrap_out) const;

    double estimate_error(const ARChainNState& high_order,
                          const ARChainNState& low_order) const;

    double step_scale_factor(double error, int k_opt) const;
};