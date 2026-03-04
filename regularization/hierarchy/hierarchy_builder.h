#pragma once
#include <vector>
#include <tuple>
#include <memory>
#include "hierarchy_node.h"
#include "nbody_system.h"

class HierarchyBuilder {
public:
    struct Params {
        double r_ks_threshold      = 1.0;  
        double tidal_threshold     = 0.1;  
        double strong_coupling_eta = 3.0;  
        double ar_chain_threshold  = 0.1;  
        double ar_chain_eta        = 1e-3;
    };

    explicit HierarchyBuilder(const Params& p = Params{}) : params(p) {}

    std::unique_ptr<HierarchyNode> build(const NBodySystem& system) const;

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

    double compute_sep_min(
        const NBodySystem& system, int i, int j, int k) const;

private:
    Params params;
};
