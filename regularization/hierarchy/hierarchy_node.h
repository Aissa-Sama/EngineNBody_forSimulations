// regularization/hierarchy/hierarchy_node.h
// FASE 6A: pn_active para activación automática de PN.
// FASE 7A: TRIPLE_AR_CHAIN fusionado en GROUP_AR_CHAIN (N >= 2 cuerpos).
#pragma once
#include <vector>
#include <memory>
#include <variant>
#include "vec3.h"
#include "chain_state.h"
#include "ks_state.h"

// ============================================================================
// NODO DEL ÁRBOL DE JERARQUÍA
//
// FASE 7A — GROUP_AR_CHAIN:
//   Reemplaza TRIPLE_AR_CHAIN. Soporta N >= 2 cuerpos en AR-chain TTL+GBS.
//   body_indices.size() == 2 → par  (antes PAIR_KS en régimen cercano)
//   body_indices.size() == 3 → triple (antes TRIPLE_AR_CHAIN)
//   body_indices.size() >= 4 → cuádruple, quíntuples, etc. (nuevo en Fase 7A)
//
//   El integrador despacha siempre a ARChainNKSIntegrator con los N índices.
//   No hay código especial por tamaño de grupo.
//
// FASE 6A — pn_active:
//   El builder marca pn_active=true cuando sep_min < r_PN.
//   El integrador elige ARChainNPNBSIntegrator cuando pn_active=true.
// ============================================================================
struct HierarchyNode {

    enum class Type {
        LEAF,           ///< Cuerpo individual → integrador de campo
        PAIR_KS,        ///< Par ligado (sep moderada) → KS simple o perturbado
        TRIPLE_CHAIN,   ///< Triple (sep moderada) → Chain3 KS (legacy)
        GROUP_AR_CHAIN, ///< N >= 2 cuerpos cercanos → AR-chain TTL+GBS+KS
        TRIPLE_AR_CHAIN = GROUP_AR_CHAIN,  ///< Alias legacy (Fase 4/6) — usar GROUP_AR_CHAIN
        COMPOSITE       ///< Contiene sub-nodos
    };

    Type type = Type::LEAF;

    std::vector<int> body_indices;
    std::vector<std::unique_ptr<HierarchyNode>> children;

    double binding_energy  = 0.0;
    double tidal_parameter = 0.0;
    bool   pn_active       = false;

    // Estado regularizado — GROUP_AR_CHAIN usa ARChainNState (gestionado
    // en ar_chain_states_ del integrador, no aquí)
    std::variant<
        std::monostate,
        KSState,
        Chain3State
    > regularized_state;

    // ── Factories ──────────────────────────────────────────────────────────

    static std::unique_ptr<HierarchyNode> make_leaf(int i) {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::LEAF;
        n->body_indices = {i};
        return n;
    }

    static std::unique_ptr<HierarchyNode> make_pair_ks(
        int i, int j, double binding_e, double tidal, bool pn = false)
    {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::PAIR_KS;
        n->body_indices = {i, j};
        n->binding_energy  = binding_e;
        n->tidal_parameter = tidal;
        n->pn_active       = pn;
        return n;
    }

    static std::unique_ptr<HierarchyNode> make_triple_chain(
        int i, int j, int k, double binding_e, bool pn = false)
    {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::TRIPLE_CHAIN;
        n->body_indices = {i, j, k};
        n->binding_energy = binding_e;
        n->pn_active      = pn;
        return n;
    }

    /// Factory principal de Fase 7A.
    /// Acepta cualquier numero de indices (N >= 2).
    /// Reemplaza make_triple_ar_chain de Fase 4/6.
    static std::unique_ptr<HierarchyNode> make_group_ar_chain(
        std::vector<int> indices, double binding_e, bool pn = false)
    {
        auto n = std::make_unique<HierarchyNode>();
        n->type           = Type::GROUP_AR_CHAIN;
        n->body_indices   = std::move(indices);
        n->binding_energy = binding_e;
        n->pn_active      = pn;
        return n;
    }

    /// Alias de compatibilidad para codigo existente que llama make_triple_ar_chain.
    static std::unique_ptr<HierarchyNode> make_triple_ar_chain(
        int i, int j, int k, double binding_e, bool pn = false)
    {
        return make_group_ar_chain({i, j, k}, binding_e, pn);
    }

    static std::unique_ptr<HierarchyNode> make_composite(
        std::vector<std::unique_ptr<HierarchyNode>> ch)
    {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::COMPOSITE;
        for (auto& c : ch)
            for (int idx : c->body_indices)
                n->body_indices.push_back(idx);
        n->children = std::move(ch);
        return n;
    }

    // ── Utilidades ─────────────────────────────────────────────────────────

    bool is_regularized() const {
        return type == Type::PAIR_KS
            || type == Type::TRIPLE_CHAIN
            || type == Type::GROUP_AR_CHAIN;
    }

    int size() const { return static_cast<int>(body_indices.size()); }
};
