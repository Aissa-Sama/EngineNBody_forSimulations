// regularization/chain/archain_n_pn_bs_integrator.h
#pragma once
#include "archain_n_pn_integrator.h"
#include "archain_n_bs_integrator.h"

/**
 * @brief Integrador Bulirsch-Stoer sobre AR-chain con correcciones PN.
 *
 * Combina ARChainNPNIntegrator (aceleraciones PN) con ARChainNBSIntegrator
 * (extrapolación GBS). El MMP generalizado (Mikkola & Merritt 2006) maneja
 * correctamente los términos dependientes de velocidad (PN2.5 disipativo)
 * porque no asume estructura Hamiltoniana.
 *
 * JERARQUÍA DE HERENCIA:
 *   ARChainNIntegrator            — base: leapfrog TTL puro
 *   ├── ARChainNPNIntegrator      — override compute_accelerations (PN)
 *   └── ARChainNBSIntegrator      — GBS sobre TTL
 *       └── ARChainNPNBSIntegrator — GBS + PN (hereda de ambos vía PN)
 *
 * DISEÑO:
 *   ARChainNPNBSIntegrator hereda de ARChainNPNIntegrator (para tener
 *   compute_accelerations con PN) y reimplementa el MMP y la extrapolación
 *   del mismo modo que ARChainNBSIntegrator, pero usando la clase PN como base.
 *
 * NOTA: La herencia múltiple diamond (ambas ramas convergen en ARChainNIntegrator)
 *   se resuelve con herencia virtual en C++. Para evitar complejidad, se usa
 *   composición: ARChainNPNBSIntegrator hereda solo de ARChainNPNIntegrator
 *   y duplica la lógica GBS del integrador BS base.
 */
class ARChainNPNBSIntegrator : public ARChainNPNIntegrator {
public:

    struct BSParameters {
        double bs_eps     = 1e-10;
        int    k_max      = 8;
        double initial_ds = 1e-3;
        double min_ds     = 1e-14;
        double max_ds     = 0.5;
        double safety     = 0.9;
        double energy_tol = 0.01;
        int    max_steps  = 1000000;
        bool   verbose    = false;
    };

    explicit ARChainNPNBSIntegrator(double eta, double c_speed, int pn_order,
                                    const BSParameters& params);

    explicit ARChainNPNBSIntegrator(double eta = 1e-3,
                                    double c_speed = 1e4,
                                    int    pn_order = 1)
        : ARChainNPNBSIntegrator(eta, c_speed, pn_order, BSParameters{}) {}

    void integrate_to_bs(ARChainNState& state, double t_abs_final);

    bool bs_step(ARChainNState& state,
                 double ds_total,
                 double& error_out,
                 double& ds_next_out);

    double energy_error(const ARChainNState& state, double E_ref) const;

private:
    BSParameters params_;

    static constexpr int K_MAX = 8;
    static constexpr std::array<int, K_MAX> n_seq_ = {2, 4, 6, 8, 10, 12, 14, 16};

    void modified_midpoint_step(const ARChainNState& y0,
                                double ds_total, int n,
                                ARChainNState& result) const;

    void polynomial_extrapolation(const std::vector<ARChainNState>& raw,
                                  const std::vector<double>& step_sizes,
                                  int k_current,
                                  ARChainNState& extrap_out) const;

    double estimate_error(const ARChainNState& hi,
                          const ARChainNState& lo) const;

    double step_scale_factor(double error, int k_opt) const;
};
