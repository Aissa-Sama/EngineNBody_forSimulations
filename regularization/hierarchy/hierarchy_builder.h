// regularization/hierarchy/hierarchy_builder.h
#pragma once
#include <vector>
#include <tuple>
#include <memory>
#include "hierarchy_node.h"
#include "nbody_system.h"

// ============================================================================
// HIERARCHY BUILDER
//
// Construye el árbol de jerarquía a partir del estado físico del sistema.
//
// Algoritmo (Mikkola & Aarseth):
//   1. find_close_pairs()        → pares con r_ij < r_threshold
//   2. filter_bound()            → pares con E_bind < 0
//   3. find_triples()            → dos pares que comparten un cuerpo
//   4. classify_triple()         → AR-chain si sep_min < umbral,
//                                   KS-chain si acoplamiento moderado
//   5. compute_tidal_parameter() → decide KS simple vs KS perturbado
//   6. Construir árbol con los nodos resultantes
//
// DECISIÓN AR-chain vs KS-chain:
//   sep_min < ar_chain_threshold → TRIPLE_AR_CHAIN
//   sep_min ≥ ar_chain_threshold → TRIPLE_CHAIN
//
//   Justificación del umbral 0.1:
//     Experimentos numéricos (Nota Técnica Interna, Marzo 2026) muestran que
//     KS-chain falla para sep_min < ~0.15. Un umbral de 0.1 activa AR-chain
//     solo cuando el encuentro es genuinamente cercano, evitando el costo
//     de AR-chain cuando KS-chain es suficiente.
// ============================================================================
class HierarchyBuilder {
public:
    // -----------------------------------------------------------------------
    // PARÁMETROS DE CONSTRUCCIÓN
    // -----------------------------------------------------------------------
    struct Params {
        double r_ks_threshold      = 1.0;  ///< Radio para pares cercanos
        double tidal_threshold     = 0.1;  ///< ε < thr → KS simple; else → KS perturbado
        double strong_coupling_eta = 3.0;  ///< r_3cm < eta * r_binary → triple fuerte

        // Umbral para elegir AR-chain vs KS-chain en triples
        // sep_min < ar_chain_threshold → AR-chain TTL (robusto para encuentros reales)
        // sep_min ≥ ar_chain_threshold → KS-chain     (más rápido, separación moderada)
        double ar_chain_threshold  = 0.5;  ///< u.a. — calibrado en experimentos (usado 0.1 anteriormente)

        // Parámetro η del ARChain3Integrator: ds = η · sep_min / Ω
        // η = 1e-3 → ~8k pasos/período, dE ~ 0.2 (rápido)
        // η = 1e-4 → ~90k pasos/período, dE ~ 0.03 (preciso)
        double ar_chain_eta        = 1e-3;
    };

    explicit HierarchyBuilder(const Params& p = Params{}) : params(p) {}

    // -----------------------------------------------------------------------
    // INTERFAZ PRINCIPAL
    // -----------------------------------------------------------------------

    /// Construye el árbol completo a partir del sistema.
    std::unique_ptr<HierarchyNode> build(const NBodySystem& system) const;

    // -----------------------------------------------------------------------
    // MÉTODOS PÚBLICOS (visibles para tests)
    // -----------------------------------------------------------------------

    std::vector<std::pair<int,int>> find_close_pairs(
        const NBodySystem& system) const;

    double compute_binding_energy(
        const NBodySystem& system, int i, int j) const;

    double compute_tidal_parameter(
        const NBodySystem& system, int i, int j) const;

    std::vector<std::pair<int,int>> select_bound_pairs(
        const NBodySystem& system) const;

    std::vector<std::tuple<int,int,int>> find_triples(
        const std::vector<std::pair<int,int>>& pairs) const;

    bool triple_is_strongly_coupled(
        const NBodySystem& system, int i, int j, int k) const;

    /// Separación mínima entre los tres cuerpos del triple.
    double compute_sep_min(
        const NBodySystem& system, int i, int j, int k) const;

private:
    Params params;
};