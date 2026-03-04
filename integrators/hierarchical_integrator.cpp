#include "hierarchical_integrator.h"
#include "binary_state.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

HierarchicalIntegrator::HierarchicalIntegrator(
    std::unique_ptr<Integrator> far_integrator,
    double r_ks_threshold,
    double ks_internal_dt,
    const HierarchyBuilder::Params& builder_params,
    RegimeLogger* logger_
)
    : far(std::move(far_integrator))
    , ks_simple(ks_internal_dt)
    , ks_perturbed(ks_internal_dt)
    , chain3(1e-4)
    , ar_chain_(builder_params.ar_chain_eta)
    , builder(builder_params)
    , tidal_threshold(builder_params.tidal_threshold)
    , logger(logger_)
{
    (void)r_ks_threshold;
}

void HierarchicalIntegrator::step(
    NBodySystem& system,
    double dt,
    const std::vector<bool>& 
) {
    const int N = static_cast<int>(system.bodies.size());
    std::vector<bool> in_subsystem(N, false);

    last_root = builder.build(system);
    {
        std::vector<TripleKey> active_keys;
        collect_active_ar_keys(*last_root, active_keys);

        for (auto it = ar_chain_states_.begin(); it != ar_chain_states_.end(); ) {
            bool found = std::find(active_keys.begin(), active_keys.end(),
                                   it->first) != active_keys.end();
            it = found ? std::next(it) : ar_chain_states_.erase(it);
        }
    }
    integrate_node(*last_root, system, dt, in_subsystem);
    far->step(system, dt, in_subsystem);
    ++step_counter;
}

void HierarchicalIntegrator::collect_active_ar_keys(
    const HierarchyNode& node,
    std::vector<TripleKey>& keys
) const {
    if (node.type == HierarchyNode::Type::TRIPLE_AR_CHAIN) {
        keys.push_back(make_key(
            node.body_indices[0],
            node.body_indices[1],
            node.body_indices[2]
        ));
    }
    for (const auto& child : node.children)
        collect_active_ar_keys(*child, keys);
}

void HierarchicalIntegrator::integrate_node(
    HierarchyNode& node,
    NBodySystem& system,
    double dt,
    std::vector<bool>& in_subsystem
) {
    switch (node.type) {
        case HierarchyNode::Type::LEAF:
            integrate_leaf(node, system, dt, in_subsystem);
            break;
        case HierarchyNode::Type::PAIR_KS:
            integrate_pair_ks(node, system, dt);
            for (int idx : node.body_indices) in_subsystem[idx] = true;
            break;
        case HierarchyNode::Type::TRIPLE_CHAIN:
            integrate_triple_chain(node, system, dt);
            for (int idx : node.body_indices) in_subsystem[idx] = true;
            break;
        case HierarchyNode::Type::TRIPLE_AR_CHAIN:
            integrate_triple_ar_chain(node, system, dt);
            for (int idx : node.body_indices) in_subsystem[idx] = true;
            break;
        case HierarchyNode::Type::COMPOSITE:
            integrate_composite(node, system, dt, in_subsystem);
            break;
    }
}

void HierarchicalIntegrator::integrate_leaf(
    HierarchyNode& node,
    NBodySystem& ,
    double ,
    std::vector<bool>& 
) {
    (void)node;
}

void HierarchicalIntegrator::integrate_pair_ks(
    HierarchyNode& node,
    NBodySystem& system,
    double dt
) {
    const int i = node.body_indices[0];
    const int j = node.body_indices[1];
    if (logger) {
        logger->log({
            static_cast<int>(step_counter), i, j,
            node.tidal_parameter < tidal_threshold
                ? "ENTER_KS_SIMPLE" : "ENTER_KS_PERTURBED"
        });
    }

    BinaryState state(system.bodies[i], system.bodies[j]);
    if (node.tidal_parameter < tidal_threshold)
        ks_simple.integrate(state, dt);
    else
        ks_perturbed.integrate_perturbed(state, dt, system, i, j);
    state.write_back(system.bodies[i], system.bodies[j]);
}

void HierarchicalIntegrator::integrate_triple_chain(
    HierarchyNode& node,
    NBodySystem& system,
    double dt
) {
    const int i = node.body_indices[0];
    const int j = node.body_indices[1];
    const int k = node.body_indices[2];
    if (logger)
        logger->log({static_cast<int>(step_counter), i, j, "ENTER_CHAIN3"});
    Chain3State state = chain3.initialize(system, i, j, k);
    double t_target   = state.cm_time + dt;
    double t_achieved = 0.0;
    IntegrationParams params;
    params.abs_tol  = 1e-10;
    params.min_dtau = 1e-8;
    params.max_dtau = 1e-1;
    chain3.integrate(state, t_target, t_achieved, params, system);
    chain3.write_back(state, system, i, j, k);
    if (logger)
        logger->log({static_cast<int>(step_counter), i, j, "EXIT_CHAIN3"});
}

void HierarchicalIntegrator::integrate_triple_ar_chain(
    HierarchyNode& node,
    NBodySystem& system,
    double dt
) {
    const int i = node.body_indices[0];
    const int j = node.body_indices[1];
    const int k = node.body_indices[2];
    const TripleKey key = make_key(i, j, k);
    if (logger)
        logger->log({static_cast<int>(step_counter), i, j, "ENTER_AR_CHAIN"});
    auto it = ar_chain_states_.find(key);
    ARChain3State state;
    if (it != ar_chain_states_.end()) {
        state = it->second;
        state.cm_pos = state.cm_pos + state.cm_vel * state.t_phys;
        state.t_phys = 0.0;
    } else {
        state = ar_chain_.initialize(system, i, j, k);
    }
    ar_chain_.integrate(state, dt);
    ar_chain_states_[key] = state;
    ar_chain_.write_back(state, system, i, j, k);
    if (logger)
        logger->log({static_cast<int>(step_counter), i, j, "EXIT_AR_CHAIN"});
}

void HierarchicalIntegrator::integrate_composite(
    HierarchyNode& node,
    NBodySystem& system,
    double dt,
    std::vector<bool>& in_subsystem
) {
    for (auto& child : node.children)
        integrate_node(*child, system, dt, in_subsystem);
}
