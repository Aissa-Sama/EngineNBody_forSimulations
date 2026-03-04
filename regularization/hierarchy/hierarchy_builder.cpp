#include "hierarchy_builder.h"
#include <algorithm>
#include <cmath>
#include <limits>

std::vector<std::pair<int,int>>
HierarchyBuilder::find_close_pairs(const NBodySystem& system) const {
    std::vector<std::pair<int,int>> pairs;
    const int N = static_cast<int>(system.bodies.size());
    for (int i = 0; i < N; ++i)
        for (int j = i+1; j < N; ++j) {
            Vec3 dr = system.bodies[j].position - system.bodies[i].position;
            if (norm(dr) < params.r_ks_threshold)
                pairs.push_back({i, j});
        }
    return pairs;
}

double HierarchyBuilder::compute_binding_energy(
    const NBodySystem& system, int i, int j) const
{
    const auto& a = system.bodies[i];
    const auto& b = system.bodies[j];
    Vec3 r = b.position - a.position;
    Vec3 v = b.velocity - a.velocity;
    double mu = (a.mass * b.mass) / (a.mass + b.mass);
    return 0.5 * mu * dot(v, v)
         - system.G * a.mass * b.mass / norm(r);
}

double HierarchyBuilder::compute_tidal_parameter(
    const NBodySystem& system, int i, int j) const
{
    const int N = static_cast<int>(system.bodies.size());
    const auto& a = system.bodies[i];
    const auto& b = system.bodies[j];

    Vec3 r_int = b.position - a.position;
    double d_int = norm(r_int);
    double F_int = system.G * a.mass * b.mass / (d_int * d_int);
    if (F_int < 1e-30) return 1e30;

    Vec3 cm = (a.mass * a.position + b.mass * b.position) / (a.mass + b.mass);

    Vec3 F_ext_vec{0, 0, 0};
    for (int k = 0; k < N; ++k) {
        if (k == i || k == j) continue;
        Vec3 r_ext = system.bodies[k].position - cm;
        double d_ext = norm(r_ext);
        double f = system.G * (a.mass + b.mass) * system.bodies[k].mass
                   / (d_ext * d_ext);
        F_ext_vec = F_ext_vec + (f / d_ext) * r_ext;
    }

    return norm(F_ext_vec) / F_int;
}

std::vector<std::pair<int,int>>
HierarchyBuilder::select_bound_pairs(const NBodySystem& system) const {
    auto close = find_close_pairs(system);

    struct Candidate { int i, j; double binding_energy; };
    std::vector<Candidate> candidates;
    for (auto [i, j] : close) {
        double E = compute_binding_energy(system, i, j);
        if (E < 0.0) candidates.push_back({i, j, E});
    }

    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate& a, const Candidate& b){
            return a.binding_energy < b.binding_energy;
        });

    const int N = static_cast<int>(system.bodies.size());
    std::vector<bool> used(N, false);
    std::vector<std::pair<int,int>> selected;
    for (const auto& c : candidates) {
        if (!used[c.i] && !used[c.j]) {
            used[c.i] = used[c.j] = true;
            selected.push_back({c.i, c.j});
        }
    }
    return selected;
}

std::vector<std::tuple<int,int,int>>
HierarchyBuilder::find_triples(
    const std::vector<std::pair<int,int>>& pairs) const
{
    std::vector<std::tuple<int,int,int>> triples;
    const int M = static_cast<int>(pairs.size());
    for (int a = 0; a < M; ++a) {
        for (int b = a+1; b < M; ++b) {
            auto [ai, aj] = pairs[a];
            auto [bi, bj] = pairs[b];

            if      (ai == bi) triples.push_back({aj, ai, bj});
            else if (ai == bj) triples.push_back({aj, ai, bi});
            else if (aj == bi) triples.push_back({ai, aj, bj});
            else if (aj == bj) triples.push_back({ai, aj, bi});
        }
    }
    return triples;
}

bool HierarchyBuilder::triple_is_strongly_coupled(
    const NBodySystem& system, int i, int j, int k) const
{
    const auto& bi = system.bodies[i];
    const auto& bj = system.bodies[j];
    const auto& bk = system.bodies[k];

    double r_binary = norm(bj.position - bi.position);
    Vec3 cm_ij = (bi.mass * bi.position + bj.mass * bj.position)
                 / (bi.mass + bj.mass);
    double r_3cm = norm(bk.position - cm_ij);

    return r_3cm < params.strong_coupling_eta * r_binary;
}

double HierarchyBuilder::compute_sep_min(
    const NBodySystem& system, int i, int j, int k) const
{
    double r12 = norm(system.bodies[j].position - system.bodies[i].position);
    double r23 = norm(system.bodies[k].position - system.bodies[j].position);
    double r13 = norm(system.bodies[k].position - system.bodies[i].position);
    return std::min({r12, r23, r13});
}

std::unique_ptr<HierarchyNode>
HierarchyBuilder::build(const NBodySystem& system) const {
    const int N = static_cast<int>(system.bodies.size());
    std::vector<bool> used(N, false);
    std::vector<std::unique_ptr<HierarchyNode>> nodes;

    auto close = find_close_pairs(system);

    struct Candidate { int i, j; double E; };
    std::vector<Candidate> all_bound;
    for (auto [i, j] : close) {
        double E = compute_binding_energy(system, i, j);
        if (E < 0.0) all_bound.push_back({i, j, E});
    }
    std::vector<std::pair<int,int>> bound_pairs;
    for (const auto& c : all_bound)
        bound_pairs.push_back({c.i, c.j});

    auto triples = find_triples(bound_pairs);

    for (const auto& [ti, tj, tk] : triples) {
        if (used[ti] || used[tj] || used[tk]) continue;

        if (!triple_is_strongly_coupled(system, ti, tj, tk)) continue;

        const auto& a = system.bodies[ti];
        const auto& b = system.bodies[tj];
        const auto& c = system.bodies[tk];
        Vec3 vcm = (a.mass * a.velocity + b.mass * b.velocity + c.mass * c.velocity)
                   / (a.mass + b.mass + c.mass);
        double T = 0.5 * a.mass * dot(a.velocity - vcm, a.velocity - vcm)
                 + 0.5 * b.mass * dot(b.velocity - vcm, b.velocity - vcm)
                 + 0.5 * c.mass * dot(c.velocity - vcm, c.velocity - vcm);
        double U = -system.G * a.mass * b.mass / norm(b.position - a.position)
                  - system.G * b.mass * c.mass / norm(c.position - b.position)
                  - system.G * a.mass * c.mass / norm(c.position - a.position);
        double binding_e = T + U;

        double sep_min = compute_sep_min(system, ti, tj, tk);

        if (sep_min < params.ar_chain_threshold) {
            nodes.push_back(
                HierarchyNode::make_triple_ar_chain(ti, tj, tk, binding_e));
        } else {
            nodes.push_back(
                HierarchyNode::make_triple_chain(ti, tj, tk, binding_e));
        }
        used[ti] = used[tj] = used[tk] = true;
    }

    std::sort(all_bound.begin(), all_bound.end(),
        [](const Candidate& a, const Candidate& b){ return a.E < b.E; });

    for (const auto& c : all_bound) {
        if (used[c.i] || used[c.j]) continue;
        double eps = compute_tidal_parameter(system, c.i, c.j);
        nodes.push_back(HierarchyNode::make_pair_ks(c.i, c.j, c.E, eps));
        used[c.i] = used[c.j] = true;
    }

    for (int i = 0; i < N; ++i)
        if (!used[i])
            nodes.push_back(HierarchyNode::make_leaf(i));

    if (nodes.size() == 1) return std::move(nodes[0]);
    return HierarchyNode::make_composite(std::move(nodes));
}
