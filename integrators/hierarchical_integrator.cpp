// integrators/hierarchical_integrator.cpp
#include "hierarchical_integrator.h"
#include "binary_state.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <functional>

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
{
    (void)r_ks_threshold;
}

// ============================================================================
// IS_PURE_AR_TRIPLE: detecta si el árbol es un único TRIPLE_AR_CHAIN
//
// Devuelve true si la raíz es exactamente un nodo TRIPLE_AR_CHAIN con
// tres body_indices. Rellena i, j, k con esos índices.
//
// Casos que NO son triple puro y devuelven false:
//   - Raíz COMPOSITE (hay campo lejano o pares KS adicionales)
//   - Raíz TRIPLE_CHAIN (separación moderada, no AR-chain)
//   - Cualquier otro tipo de nodo
// ============================================================================
bool HierarchicalIntegrator::is_pure_ar_triple(const HierarchyNode& root,
                                                int& i, int& j, int& k) const
{
    if (root.type != HierarchyNode::Type::TRIPLE_AR_CHAIN)
        return false;
    if (root.body_indices.size() != 3)
        return false;

    i = root.body_indices[0];
    j = root.body_indices[1];
    k = root.body_indices[2];
    return true;
}

// ============================================================================
// STEP: modo bloque (uso general, N > 3 o subsistemas mixtos)
// ============================================================================
void HierarchicalIntegrator::step(
    NBodySystem& system,
    double dt,
    const std::vector<bool>& /*external_mask*/
) {
    const int N = static_cast<int>(system.bodies.size());
    std::vector<bool> in_subsystem(N, false);

    last_root = builder.build(system);
    integrate_node(*last_root, system, dt, in_subsystem);
    far->step(system, dt, in_subsystem);

    ++step_counter;
}

// ============================================================================
// STEP_TO: modo directo (sistema completo AR-chain, N=3)
//
// Integra de t=0 a t_final sin bloques de dt externos.
//
// ALGORITMO:
//   1. Construir el árbol para verificar que el sistema es un triple puro.
//   2. Si es TRIPLE_AR_CHAIN puro:
//      a. Inicializar el estado desde el sistema físico (t_phys=0).
//      b. Llamar integrate_to(state, t_final) — el AR-chain controla
//         su propio tiempo sin interrupciones de ningún tipo.
//      c. Escribir el estado final al sistema físico.
//   3. Si NO es triple puro (campo lejano activo, pares KS, etc.):
//      Caer al modo step() con dt heurístico. Garantiza que step_to()
//      nunca falla silenciosamente en sistemas inesperados.
//
// POR QUÉ ESTE MODO ES MÁS PRECISO:
//   En el modo step(), cada bloque de dt introduce un pequeño error de
//   fase. Con 633 bloques (dt=0.01) o 6326 (dt=0.001), ese error se
//   acumula y altera la trayectoria antes del encuentro cercano, haciendo
//   que el cruce ocurra en t≈2.23 con sep=0.021 en lugar del correcto
//   t≈3.04 con sep=7.7e-4. En este modo, no hay bloques: el AR-chain
//   integra directamente con su propio control adaptativo de paso.
//   Referencia: NBody_Maestro_v4 §5 (diagnóstico arquitectural).
// ============================================================================
void HierarchicalIntegrator::step_to(
    NBodySystem& system,
    double t_final,
    std::function<void(double)> logger_cb
) {
    // ── 1. Construir árbol y detectar modo ───────────────────────────────────
    last_root = builder.build(system);

    int i, j, k;
    if (!is_pure_ar_triple(*last_root, i, j, k)) {
        // El sistema no es un triple puro — caer al modo step con dt heurístico.
        // El caller debería usar step() directamente para casos mixtos.
        const double dt_fallback = t_final / 1000.0;
        double t_cur = 0.0;
        while (t_cur < t_final) {
            const double dt = std::min(dt_fallback, t_final - t_cur);
            step(system, dt, {});
            t_cur += dt;
            if (logger_cb) logger_cb(t_cur);
        }
        return;
    }

    // ── 2. Triple puro: inicializar estado ───────────────────────────────────
    if (logger)
        logger->log({static_cast<int>(step_counter), i, j, "ENTER_AR_CHAIN_DIRECT"});

    ARChain3State state = ar_chain_.initialize(system, i, j, k);
    // t_phys = 0.0 tras initialize(); integrate_to avanza hasta t_final absoluto.

    // Clave canónica para el caché (consistente con modo bloque)
    int a = i, b = j, c = k;
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    const int key = a * 10000 + b * 100 + c;

    // ── 3. Integración directa sin cortes ────────────────────────────────────
    ar_chain_.integrate_to(state, t_final);

    // ── 4. Escribir estado final al sistema físico ───────────────────────────
    ar_chain_.write_back(state, system, i, j, k);

    // Persistir estado final en el caché (diagnóstico/reinicio)
    ar_chain_states_[key] = state;

    if (logger)
        logger->log({static_cast<int>(step_counter), i, j, "EXIT_AR_CHAIN_DIRECT"});

    ++step_counter;

    // ── 5. Callback de logging (t_final) ─────────────────────────────────────
    if (logger_cb) logger_cb(t_final);
}

// ============================================================================
// DESPACHADOR RECURSIVO
// ============================================================================
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

// ============================================================================
// LEAF: marcar como libre (far lo integrará)
// ============================================================================
void HierarchicalIntegrator::integrate_leaf(
    HierarchyNode& node,
    NBodySystem& /*system*/,
    double /*dt*/,
    std::vector<bool>& /*in_subsystem*/
) {
    (void)node;
}

// ============================================================================
// PAIR_KS: KS simple o perturbado según tidal_parameter
// ============================================================================
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

// ============================================================================
// TRIPLE_CHAIN: Chain3Integrator (KS-chain, separación moderada)
// ============================================================================
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

// ============================================================================
// TRIPLE_AR_CHAIN: ARChain3Integrator en modo bloque
//
// Modo bloque: llamado desde step() con un dt fijo.
// Para el modo directo (sin bloques), usar step_to().
//
// PERSISTENCIA DEL ESTADO ENTRE BLOQUES:
//   El árbol se reconstruye en cada paso, destruyendo los nodos.
//   El estado se guarda en ar_chain_states_ indexado por clave canónica.
//
// MANEJO DEL CM ENTRE BLOQUES:
//   Al retomar el estado del bloque anterior, t_phys tiene el valor
//   final de ese bloque. Antes de integrar el nuevo bloque:
//     1. Absorber el desplazamiento del CM: cm_pos += cm_vel * t_phys
//     2. Resetear el reloj interno:         t_phys  = 0
// ============================================================================
void HierarchicalIntegrator::integrate_triple_ar_chain(
    HierarchyNode& node,
    NBodySystem& system,
    double dt
) {
    const int i = node.body_indices[0];
    const int j = node.body_indices[1];
    const int k = node.body_indices[2];

    if (logger)
        logger->log({static_cast<int>(step_counter), i, j, "ENTER_AR_CHAIN"});

    // ── Clave normalizada ────────────────────────────────────────────────────
    int a = i, b = j, c = k;
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    const int key = a * 10000 + b * 100 + c;

    ARChain3State state;

    auto it = ar_chain_states_.find(key);
    if (it != ar_chain_states_.end()) {
        state = it->second;
        // Absorber desplazamiento del CM del bloque anterior y resetear reloj
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

// ============================================================================
// COMPOSITE: integrar hijos recursivamente
// ============================================================================
void HierarchicalIntegrator::integrate_composite(
    HierarchyNode& node,
    NBodySystem& system,
    double dt,
    std::vector<bool>& in_subsystem
) {
    for (auto& child : node.children)
        integrate_node(*child, system, dt, in_subsystem);
}