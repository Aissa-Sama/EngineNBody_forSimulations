// regularization/ks/ks_perturbed_integrator.cpp
#include "ks_perturbed_integrator.h"
#include "ks_transform.h"
#include <algorithm>
#include <iostream>
#include <iomanip>

KSPerturbedIntegrator::KSPerturbedIntegrator(double internal_dt)
    : KSIntegrator(internal_dt)
{
}

std::vector<KSPerturbedIntegrator::Perturber>
KSPerturbedIntegrator::find_perturbers(
    const NBodySystem& system,
    int idx_i, int idx_j,
    double r_binary
) const {
    std::vector<Perturber> perturbers;

    double M = system.bodies[idx_i].mass + system.bodies[idx_j].mass;
    Vec3 cm_pos = (system.bodies[idx_i].position * system.bodies[idx_i].mass +
                   system.bodies[idx_j].position * system.bodies[idx_j].mass) / M;

    for (size_t k = 0; k < system.bodies.size(); ++k) {
        if (k == static_cast<size_t>(idx_i) || k == static_cast<size_t>(idx_j)) continue;

        Vec3 r_pert = system.bodies[k].position - cm_pos;
        double dist = norm(r_pert);

        if (dist < 50.0 * r_binary) {
            Perturber p;
            p.index = static_cast<int>(k);
            p.mass = system.bodies[k].mass;
            p.position = r_pert;
            p.velocity = system.bodies[k].velocity -
                         (system.bodies[idx_i].velocity + system.bodies[idx_j].velocity) / 2.0;
            p.distance = dist;
            perturbers.push_back(p);
        }
    }

    std::sort(perturbers.begin(), perturbers.end(),
              [](const Perturber& a, const Perturber& b) { return a.distance < b.distance; });

    return perturbers;
}

Vec3 KSPerturbedIntegrator::compute_perturbation_force(
    const BinaryState& binary,
    const std::vector<Perturber>& perturbers,
    const NBodySystem& system
) const {
    Vec3 F_total{0, 0, 0};

    double M = binary.total_mass();
    Vec3 r_rel = binary.relative_position();

    for (const auto& p : perturbers) {
        Vec3 r1 = - (binary.mass2() / M) * r_rel;
        Vec3 r2 =   (binary.mass1() / M) * r_rel;

        Vec3 dr1 = p.position - r1;
        Vec3 dr2 = p.position - r2;

        double d1 = norm(dr1) + 1e-10;
        double d2 = norm(dr2) + 1e-10;

        Vec3 a1 = (system.G * p.mass / (d1 * d1 * d1)) * dr1;
        Vec3 a2 = (system.G * p.mass / (d2 * d2 * d2)) * dr2;

        Vec3 F_marea = binary.reduced_mass() * (a1 - a2);
        F_total = F_total + F_marea;
    }

    return F_total;
}

void KSPerturbedIntegrator::force_to_ks(
    const Vec3& F_phys,
    const double u[4],
    double F_ks[4]
) const {
    double L[3][4];
    levi_civita_matrix(u, L);

    double u_norm_sq = ks_radius(u);

    for (int k = 0; k < 4; ++k) {
        F_ks[k] = (u_norm_sq / 4.0) * (  // reducido para evitar explosión
            L[0][k] * F_phys.x +
            L[1][k] * F_phys.y +
            L[2][k] * F_phys.z
        );
    }
}

void KSPerturbedIntegrator::leapfrog_step_perturbed(
    KSState& ks,
    const double F_ks[4],
    double dtau
) const {
    // Guardar energía antes del paso
    double h_before = ks.energy;
    
    for (int i = 0; i < 4; ++i) {
        ks.w[i] -= 0.5 * dtau * ks.energy * ks.u[i];
    }

    for (int i = 0; i < 4; ++i) {
        ks.u[i] += dtau * ks.w[i];
    }

    for (int i = 0; i < 4; ++i) {
        ks.w[i] -= 0.5 * dtau * ks.energy * ks.u[i];
        ks.w[i] += dtau * F_ks[i];
    }

    double r = ks_radius(ks.u);
    double work = 0.0;
    for (int i = 0; i < 4; ++i) {
        work += F_ks[i] * ks.w[i];
    }
    
    // dH/dτ = -r·(F_pert·v): aplicar sin limitador.
    // Si dh es inestable, la solución correcta es reducir dtau, no recortar la energía.
    double dh = -r * work * dtau;
    
    ks.energy += dh;
}

void KSPerturbedIntegrator::integrate_perturbed(
    BinaryState& binary,
    double dt_global,
    const NBodySystem& system,
    int idx_i, int idx_j
) {
    if (dt_global <= 0.0) return;

    KSState ks = to_ks(binary);

    double r0 = norm(binary.relative_position());
    if (r0 < 1e-12) {
        throw std::runtime_error("KS perturbed: r0 nulo");
    }

    auto perturbers = find_perturbers(system, idx_i, idx_j, r0);

    double F_ks[4] = {0,0,0,0};
    bool has_pert = !perturbers.empty();

    if (!has_pert) {
        integrate(binary, dt_global);
        return;
    }

    Vec3 F_phys = compute_perturbation_force(binary, perturbers, system);
    force_to_ks(F_phys, ks.u, F_ks);

    // CORRECCIÓN: relación correcta de la transformación KS: dt = 2·r·dτ → dτ = dt/(2r)
    double dtau_total = dt_global / (2.0 * r0);

    int n_steps = std::max(1000,static_cast<int>(std::ceil(dtau_total / internal_dt)));
    if (has_pert) n_steps *= 5;  // más precisión para evitar explosión
    double dtau = dtau_total / n_steps;

    double h0 = ks.energy;

    for (int i = 0; i < n_steps; ++i) {
        leapfrog_step_perturbed(ks, F_ks, dtau);

        if (i % 50 == 0) {
            double unorm = ks_radius(ks.u);
            if (unorm > 2000.0 * r0 || unorm < 0.0005 * r0) {
                std::cerr << "[KS Pert] Abandono por inestabilidad (|u|^2 = " << unorm << ")\n";
                break;
            }
        }
    }

    from_ks(ks, binary);
    binary.advance_cm(dt_global);

/*
    double dh = ks.energy - h0;
    if (std::abs(dh) > 1e-8) {
        std::cout << "[KS Pert] Δh total = " << dh << " (pert = " << perturbers.size() << ")\n";
    }
*/
    // dh diagnóstico eliminado — era debug de desarrollo
    (void)(ks.energy - h0);

}
