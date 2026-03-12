// integrators/hierarchical_integrator.h
#pragma once
#include <memory>
#include <vector>
#include <map>
#include <tuple>
#include <functional>
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

// ============================================================================
// HIERARCHICAL INTEGRATOR — Ruta B
//
// Cada paso:
//   1. HierarchyBuilder::build() → árbol del instante actual
//   2. Recorrer el árbol e integrar cada nodo con su método especializado:
//        LEAF             → far (Leapfrog/RK4)
//        PAIR_KS          → KS simple o KS perturbado (tidal_parameter)
//        TRIPLE_CHAIN     → Chain3Integrator  (sep ≥ ar_chain_threshold)
//        TRIPLE_AR_CHAIN  → ARChain3Integrator (sep < ar_chain_threshold)
//        COMPOSITE        → integrar hijos recursivamente
//
// ── DOS MODOS DE OPERACIÓN ───────────────────────────────────────────────────
//
// 1. step(system, dt) — modo bloque (uso general, N > 3 o subsistemas mixtos):
//    Avanza el sistema un intervalo dt. El integrador exterior sincroniza
//    el tiempo entre subsistemas (campo lejano, pares KS, triples AR-chain).
//    Cada subsistema se integra en un bloque de tamaño dt.
//
// 2. step_to(system, t_final) — modo directo (sistema completo AR-chain):
//    Usado cuando el árbol entero es un único TRIPLE_AR_CHAIN.
//    En lugar de dividir la integración en bloques de dt, delega directamente
//    a ARChain3Integrator::integrate_to(), que integra de t=0 a t=t_final
//    sin interrupciones. El AR-chain controla su propio tiempo.
//
//    CUÁNDO USAR step_to():
//      - N=3 con ar_thresh suficientemente alto (e.g. 99.0) para que el
//        triple esté siempre activo durante toda la simulación.
//      - La figura-8, el problema de Pitágoras, y cualquier sistema donde
//        los tres cuerpos permanecen fuertemente acoplados todo el tiempo.
//
//    POR QUÉ ES MÁS PRECISO:
//      Los cortes de bloque fijo (modo step) acumulan error de fase antes
//      del encuentro cercano, alterando la trayectoria. Con dt=0.001 la
//      figura-8 diverge en t≈2.23 (sep=0.021) en lugar del encuentro real
//      en t≈3.04 (sep=7.7e-4). Este error no se reduce reduciendo dt.
//      Referencia: diagnóstico arquitectural en NBody_Maestro_v4 (§5).
//
// ── PERSISTENCIA DEL ESTADO AR-CHAIN ─────────────────────────────────────────
//   El árbol se reconstruye desde cero en cada paso — los nodos son efímeros.
//   El ARChain3State NO puede guardarse en el nodo porque se destruye con él.
//
//   Solución: ar_chain_states_ es un std::map<int, ARChain3State>
//   indexado por clave canónica a*10000+b*100+c (con a<b<c).
//   El estado persiste en el integrador entre pasos. No se limpia entre
//   pasos para sobrevivir oscilaciones del builder entre regímenes.
// ============================================================================
class HierarchicalIntegrator : public Integrator {
public:
    HierarchicalIntegrator(
        std::unique_ptr<Integrator> far_integrator,
        double r_ks_threshold,
        double ks_internal_dt,
        const HierarchyBuilder::Params& builder_params = HierarchyBuilder::Params{},
        RegimeLogger* logger = nullptr
    );

    // ── Interfaz principal ───────────────────────────────────────────────────

    /**
     * @brief Modo bloque: avanza el sistema un intervalo dt.
     *
     * Uso general para N > 3 o cuando el sistema tiene subsistemas mixtos
     * (campo lejano + pares KS + triples AR-chain).
     */
    void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) override;

    /**
     * @brief Modo directo: integra de t=0 a t_final sin bloques externos.
     *
     * Detecta si el árbol entero es un único TRIPLE_AR_CHAIN y delega a
     * ARChain3Integrator::integrate_to(). Si el sistema no es un triple
     * puro, cae al modo step() normal.
     *
     * Uso: simulaciones donde N=3 y el triple está siempre activo
     * (ar_thresh alto). Evita el error de fase de los bloques dt.
     *
     * @param system    Sistema físico de entrada/salida.
     * @param t_final   Tiempo físico absoluto final.
     * @param logger_cb Callback opcional invocado tras cada paso de observables
     *                  (para logging externo). Firma: void(double t).
     *                  Si es nullptr no se invoca nada.
     */
    void step_to(
        NBodySystem& system,
        double t_final,
        std::function<void(double)> logger_cb = nullptr
    );

    /// Acceso al árbol del último paso (útil para tests y diagnóstico)
    const HierarchyNode* last_tree() const { return last_root.get(); }

private:
    // -----------------------------------------------------------------------
    // CLAVE CANÓNICA PARA EL CACHÉ
    // Los índices se ordenan a < b < c para que el mismo triple detectado
    // en cualquier permutación recupere el mismo estado.
    // -----------------------------------------------------------------------
    using TripleKey = std::tuple<int,int,int>;

    static TripleKey make_key(int i, int j, int k) {
        int a = i, b = j, c = k;
        if (a > b) std::swap(a, b);
        if (b > c) std::swap(b, c);
        if (a > b) std::swap(a, b);
        return {a, b, c};
    }

    // -----------------------------------------------------------------------
    // INTEGRACIÓN POR TIPO DE NODO
    // -----------------------------------------------------------------------

    void integrate_node(
        HierarchyNode& node,
        NBodySystem& system,
        double dt,
        std::vector<bool>& in_subsystem
    );

    void integrate_leaf(
        HierarchyNode& node,
        NBodySystem& system,
        double dt,
        std::vector<bool>& in_subsystem
    );

    void integrate_pair_ks(
        HierarchyNode& node,
        NBodySystem& system,
        double dt
    );

    void integrate_triple_chain(
        HierarchyNode& node,
        NBodySystem& system,
        double dt
    );

    void integrate_triple_ar_chain(
        HierarchyNode& node,
        NBodySystem& system,
        double dt
    );

    void integrate_composite(
        HierarchyNode& node,
        NBodySystem& system,
        double dt,
        std::vector<bool>& in_subsystem
    );

    // -----------------------------------------------------------------------
    // DETECCIÓN DEL MODO DIRECTO
    // Devuelve true si el árbol raíz es un único TRIPLE_AR_CHAIN con
    // los índices (i,j,k), y rellena esos índices.
    // -----------------------------------------------------------------------
    bool is_pure_ar_triple(const HierarchyNode& root,
                           int& i, int& j, int& k) const;

    // -----------------------------------------------------------------------
    // MIEMBROS
    // -----------------------------------------------------------------------
    std::unique_ptr<Integrator>   far;
    KSIntegrator                  ks_simple;
    KSPerturbedIntegrator         ks_perturbed;
    Chain3Integrator              chain3;
    ARChain3Integrator            ar_chain_;
    HierarchyBuilder              builder;
    double                        tidal_threshold;
    RegimeLogger*                 logger;
    std::size_t                   step_counter = 0;

    std::unique_ptr<HierarchyNode> last_root;

    // Mapa de persistencia de estado AR-chain entre pasos (modo bloque).
    // Clave: triple normalizado (a<b<c) como a*10000+b*100+c.
    // No se limpia entre pasos para sobrevivir oscilaciones del builder.
    std::map<int, ARChain3State> ar_chain_states_;
};