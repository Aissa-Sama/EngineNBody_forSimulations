// integrators/hierarchical_integrator.h
// FASE 6A: ARChainNPNBSIntegrator con activación automática PN.
// FASE 6C: step_to() — integra hasta t_final con callback por paso.
// FASE 7A: GROUP_AR_CHAIN — N >= 2 cuerpos en AR-chain (cuadruples, quintuples).
#pragma once
#define _USE_MATH_DEFINES
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <tuple>
#include "integrator.h"
#include "hierarchy_node.h"
#include "hierarchy_builder.h"
#include "chain3_integrator.h"
#include "archain3_integrator.h"
#include "archain3_state.h"
#include "ks_perturbed_integrator.h"
#include "ks_integrator.h"
#include "regime_logger.h"
#include "nbody_system.h"

#include "../regularization/chain/archain_n_state.h"
#include "../regularization/chain/archain_n_pn_bs_integrator.h"
// FASE 7A
#include "../regularization/chain/archain_n_ks_integrator.h"

class HierarchicalIntegrator : public Integrator {
public:
    HierarchicalIntegrator(
        std::unique_ptr<Integrator> far_integrator,
        double r_ks_threshold,
        double ks_internal_dt,
        const HierarchyBuilder::Params& builder_params = HierarchyBuilder::Params{},
        RegimeLogger* logger = nullptr
    );

    void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) override;

    void step_to(
        NBodySystem& system,
        double t_final,
        double dt_hint,
        std::function<void(double)> on_step = nullptr
    );

    void step_to(
        NBodySystem& system,
        double t_final,
        std::function<void(double)> on_step
    );

    const HierarchyNode* last_tree() const { return last_root.get(); }

    bool pn_active_for(const std::vector<int>& indices) const;
    int  pn_active_count() const;

private:
    using TripleKey = std::tuple<int,int,int>;

    static TripleKey make_triple_key(int i, int j, int k) {
        int a = i, b = j, c = k;
        if (a > b) std::swap(a, b);
        if (b > c) std::swap(b, c);
        if (a > b) std::swap(a, b);
        return {a, b, c};
    }

    // ── Despachadores ────────────────────────────────────────────────────────
    void integrate_node(
        HierarchyNode& node, NBodySystem& system, double dt,
        std::vector<bool>& in_subsystem);

    void integrate_leaf(
        HierarchyNode& node, NBodySystem& system, double dt,
        std::vector<bool>& in_subsystem);

    void integrate_pair_ks(
        HierarchyNode& node, NBodySystem& system, double dt);

    void integrate_triple_chain(
        HierarchyNode& node, NBodySystem& system, double dt);

    void integrate_triple_ar_chain(
        HierarchyNode& node, NBodySystem& system, double dt);

    // FASE 7A: N >= 2 cuerpos en GROUP_AR_CHAIN
    void integrate_group_ar_chain(
        HierarchyNode& node, NBodySystem& system, double dt);

    void integrate_pn_group(
        HierarchyNode& node, NBodySystem& system, double dt);

    void integrate_composite(
        HierarchyNode& node, NBodySystem& system, double dt,
        std::vector<bool>& in_subsystem);

    void collect_active_ar_keys(
        const HierarchyNode& node, std::vector<TripleKey>& keys) const;

    // ── Miembros ─────────────────────────────────────────────────────────────
    std::unique_ptr<Integrator>    far;
    KSIntegrator                   ks_simple;
    KSPerturbedIntegrator          ks_perturbed;
    Chain3Integrator               chain3;
    ARChain3Integrator             ar_chain_;
    HierarchyBuilder               builder;
    double                         tidal_threshold;
    RegimeLogger*                  logger;
    std::size_t                    step_counter = 0;
    double                         dt_hint_;

    std::unique_ptr<HierarchyNode> last_root;
    std::map<int, ARChain3State>   ar_chain_states_;

    // FASE 7A: integrador KS generalizado para GROUP_AR_CHAIN
    ARChainNKSIntegrator                         ar_chain_ks_;
    std::map<std::string, ARChainNState>         ar_chain_n_states_;

    // FASE 6A: integrador PN
    std::unique_ptr<ARChainNPNBSIntegrator>      pn_integrator_;
    std::map<std::string, ARChainNState>         pn_states_;
    std::map<std::string, bool>                  pn_cache_;
};
