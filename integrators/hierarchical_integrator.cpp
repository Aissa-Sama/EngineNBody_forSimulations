// integrators/hierarchical_integrator.cpp
// FASE 6A: integrate_pn_group() + lógica de despacho PN.
// FASE 6C: step_to() con callback.
// FASE 7A: integrate_group_ar_chain() — N >= 2 cuerpos en GROUP_AR_CHAIN.
#define _USE_MATH_DEFINES
#include "hierarchical_integrator.h"
#include "binary_state.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// ============================================================================
// CONSTRUCTOR
// ============================================================================
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
    , dt_hint_(ks_internal_dt * 100.0)
    , ar_chain_ks_(builder_params.ar_chain_eta)   // FASE 7A
{
    (void)r_ks_threshold;

    if (builder_params.pn.enabled) {
        ARChainNPNBSIntegrator::BSParameters bs;
        bs.bs_eps      = builder_params.pn.bs_eps;
        bs.energy_tol  = 1.0;
        bs.max_steps   = 500000;

        pn_integrator_ = std::make_unique<ARChainNPNBSIntegrator>(
            builder_params.pn.eta_pn,
            builder_params.pn.c_speed,
            builder_params.pn.pn_order,
            bs
        );
    }
}

// ============================================================================
// PASO PRINCIPAL
// ============================================================================
void HierarchicalIntegrator::step(
    NBodySystem& system, double dt, const std::vector<bool>& /*external_mask*/)
{
    const int N = static_cast<int>(system.bodies.size());
    std::vector<bool> in_subsystem(N, false);
    last_root = builder.build(system, &pn_cache_);
    integrate_node(*last_root, system, dt, in_subsystem);
    far->step(system, dt, in_subsystem);
    ++step_counter;
}

// ============================================================================
// STEP_TO (con dt_hint explícito)
// ============================================================================
void HierarchicalIntegrator::step_to(
    NBodySystem& system, double t_final, double dt_hint,
    std::function<void(double)> on_step)
{
    if (t_final <= 0.0) return;
    if (dt_hint <= 0.0)
        throw std::invalid_argument("HierarchicalIntegrator::step_to: dt_hint debe ser > 0");

    const double tol = 1e-14 * std::abs(t_final);
    double t_current = 0.0;
    const int N = static_cast<int>(system.bodies.size());

    while (t_current < t_final - tol) {
        const double dt = std::min(dt_hint, t_final - t_current);
        std::vector<bool> in_subsystem(N, false);
        last_root = builder.build(system, &pn_cache_);
        integrate_node(*last_root, system, dt, in_subsystem);
        far->step(system, dt, in_subsystem);
        ++step_counter;
        t_current += dt;
        if (on_step) on_step(t_current);
    }
}

// ============================================================================
// STEP_TO (compatibilidad — sin dt_hint)
// ============================================================================
void HierarchicalIntegrator::step_to(
    NBodySystem& system, double t_final,
    std::function<void(double)> on_step)
{
    step_to(system, t_final, dt_hint_, std::move(on_step));
}

// ============================================================================
// INSPECCIÓN FASE 6A
// ============================================================================
bool HierarchicalIntegrator::pn_active_for(const std::vector<int>& indices) const {
    std::string key = HierarchyBuilder::make_group_key(
        std::vector<int>(indices.begin(), indices.end()));
    auto it = pn_cache_.find(key);
    return (it != pn_cache_.end()) ? it->second : false;
}

int HierarchicalIntegrator::pn_active_count() const {
    int count = 0;
    for (const auto& [k, v] : pn_cache_)
        if (v) ++count;
    return count;
}

// ============================================================================
// DESPACHADOR RECURSIVO
// ============================================================================
void HierarchicalIntegrator::integrate_node(
    HierarchyNode& node, NBodySystem& system, double dt,
    std::vector<bool>& in_subsystem)
{
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
        case HierarchyNode::Type::GROUP_AR_CHAIN:   // FASE 7A (incluye TRIPLE_AR_CHAIN)
            integrate_group_ar_chain(node, system, dt);
            for (int idx : node.body_indices) in_subsystem[idx] = true;
            break;
        case HierarchyNode::Type::COMPOSITE:
            integrate_composite(node, system, dt, in_subsystem);
            break;
    }
}

// ============================================================================
// LEAF
// ============================================================================
void HierarchicalIntegrator::integrate_leaf(
    HierarchyNode& node, NBodySystem&, double, std::vector<bool>&)
{
    (void)node;
}

// ============================================================================
// PAIR_KS
// ============================================================================
void HierarchicalIntegrator::integrate_pair_ks(
    HierarchyNode& node, NBodySystem& system, double dt)
{
    if (node.pn_active && pn_integrator_) {
        if (logger) logger->log({static_cast<int>(step_counter),
            node.body_indices[0], node.body_indices[1], "ENTER_PAIR_PN"});
        integrate_pn_group(node, system, dt);
        if (logger) logger->log({static_cast<int>(step_counter),
            node.body_indices[0], node.body_indices[1], "EXIT_PAIR_PN"});
        return;
    }

    const int i = node.body_indices[0];
    const int j = node.body_indices[1];

    if (logger) logger->log({static_cast<int>(step_counter), i, j,
        node.tidal_parameter < tidal_threshold
            ? "ENTER_KS_SIMPLE" : "ENTER_KS_PERTURBED"});

    BinaryState state(system.bodies[i], system.bodies[j]);

    if (node.tidal_parameter < tidal_threshold)
        ks_simple.integrate(state, dt);
    else
        ks_perturbed.integrate_perturbed(state, dt, system, i, j);

    state.write_back(system.bodies[i], system.bodies[j]);
}

// ============================================================================
// TRIPLE_CHAIN (legacy)
// ============================================================================
void HierarchicalIntegrator::integrate_triple_chain(
    HierarchyNode& node, NBodySystem& system, double dt)
{
    const int i = node.body_indices[0];
    const int j = node.body_indices[1];
    const int k = node.body_indices[2];

    if (logger) logger->log({static_cast<int>(step_counter), i, j, "ENTER_CHAIN3"});

    Chain3State state = chain3.initialize(system, i, j, k);
    double t_target = state.cm_time + dt, t_achieved = 0.0;
    IntegrationParams params;
    params.abs_tol = 1e-10; params.min_dtau = 1e-8; params.max_dtau = 1e-1;
    chain3.integrate(state, t_target, t_achieved, params, system);
    chain3.write_back(state, system, i, j, k);

    if (logger) logger->log({static_cast<int>(step_counter), i, j, "EXIT_CHAIN3"});
}

// ============================================================================
// TRIPLE_AR_CHAIN (legacy — ARChain3, mantenido para compatibilidad)
// ============================================================================
void HierarchicalIntegrator::integrate_triple_ar_chain(
    HierarchyNode& node, NBodySystem& system, double dt)
{
    if (node.pn_active && pn_integrator_) {
        if (logger) logger->log({static_cast<int>(step_counter),
            node.body_indices[0], node.body_indices[1], "ENTER_TRIPLE_PN"});
        integrate_pn_group(node, system, dt);
        if (logger) logger->log({static_cast<int>(step_counter),
            node.body_indices[0], node.body_indices[1], "EXIT_TRIPLE_PN"});
        return;
    }

    const int i = node.body_indices[0];
    const int j = node.body_indices[1];
    const int k = node.body_indices[2];

    if (logger) logger->log({static_cast<int>(step_counter), i, j, "ENTER_AR_CHAIN"});

    int a = i, b = j, c = k;
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    const int key = a * 10000 + b * 100 + c;

    ARChain3State state;
    auto it = ar_chain_states_.find(key);
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

    if (logger) logger->log({static_cast<int>(step_counter), i, j, "EXIT_AR_CHAIN"});
}

// ============================================================================
// FASE 7A — GROUP_AR_CHAIN
//
// Integra N >= 2 cuerpos con ARChainNKSIntegrator.
// Clave canonica: indices ordenados separados por '_' (e.g. "0_1_2_3").
// Estado persistido en ar_chain_n_states_ entre pasos.
//
// PATRON CM:
//   Al retomar el estado del bloque anterior, absorber el desplazamiento
//   de CM acumulado (cm_vel * t_phys) y resetear t_phys = 0.
//   Mismo patron que integrate_triple_ar_chain().
// ============================================================================
void HierarchicalIntegrator::integrate_group_ar_chain(
    HierarchyNode& node, NBodySystem& system, double dt)
{
    // FASE 6A: si PN activo, delegar
    if (node.pn_active && pn_integrator_) {
        if (logger) logger->log({static_cast<int>(step_counter),
            node.body_indices[0], node.body_indices[1], "ENTER_GROUP_PN"});
        integrate_pn_group(node, system, dt);
        if (logger) logger->log({static_cast<int>(step_counter),
            node.body_indices[0], node.body_indices[1], "EXIT_GROUP_PN"});
        return;
    }

    const std::vector<int>& indices = node.body_indices;

    // Clave canónica: índices ordenados
    std::string key = HierarchyBuilder::make_group_key(
        std::vector<int>(indices.begin(), indices.end()));

    if (logger) logger->log({static_cast<int>(step_counter),
        indices[0], indices.back(), "ENTER_GROUP_AR_CHAIN"});

    ARChainNState state;
    auto it = ar_chain_n_states_.find(key);
    if (it != ar_chain_n_states_.end()) {
        state = it->second;
        // Absorber desplazamiento CM del bloque anterior
        state.cm_pos = state.cm_pos + state.cm_vel * state.t_phys;
        state.t_phys = 0.0;
    } else {
        state = ar_chain_ks_.initialize(system, indices);
    }

    ar_chain_ks_.integrate_to_ks(state, state.t_phys + dt);

    ar_chain_n_states_[key] = state;
    ar_chain_ks_.write_back(state, system, indices);

    if (logger) logger->log({static_cast<int>(step_counter),
        indices[0], indices.back(), "EXIT_GROUP_AR_CHAIN"});
}

// ============================================================================
// FASE 6A — INTEGRATE_PN_GROUP
// ============================================================================
void HierarchicalIntegrator::integrate_pn_group(
    HierarchyNode& node, NBodySystem& system, double dt)
{
    if (!pn_integrator_) return;

    const std::string key = HierarchyBuilder::make_group_key(node.body_indices);

    ARChainNState state;
    auto it = pn_states_.find(key);
    if (it != pn_states_.end()) {
        state = it->second;
        state.cm_pos = state.cm_pos + state.cm_vel * state.t_phys;
        state.t_phys = 0.0;
    } else {
        state = pn_integrator_->initialize(system, node.body_indices);
    }

    pn_integrator_->integrate_to_bs(state, state.t_phys + dt);
    pn_states_[key] = state;
    pn_integrator_->write_back(state, system, node.body_indices);
}

// ============================================================================
// COMPOSITE
// ============================================================================
void HierarchicalIntegrator::integrate_composite(
    HierarchyNode& node, NBodySystem& system, double dt,
    std::vector<bool>& in_subsystem)
{
    for (auto& child : node.children)
        integrate_node(*child, system, dt, in_subsystem);
}

// ============================================================================
// COLLECT_ACTIVE_AR_KEYS (legacy)
// ============================================================================
void HierarchicalIntegrator::collect_active_ar_keys(
    const HierarchyNode& node, std::vector<TripleKey>& keys) const
{
    if (node.type == HierarchyNode::Type::TRIPLE_AR_CHAIN
        && node.body_indices.size() == 3) {
        keys.push_back(make_triple_key(
            node.body_indices[0], node.body_indices[1], node.body_indices[2]));
    }
    for (const auto& child : node.children)
        collect_active_ar_keys(*child, keys);
}
