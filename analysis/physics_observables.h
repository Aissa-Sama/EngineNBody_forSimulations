// analysis/physics_observables.h
// ============================================================================
// OBSERVABLES FÍSICOS — FUNCIONES PURAS
// ============================================================================
//
// Funciones libres que calculan invariantes físicos a partir del estado
// del sistema. Son independientes de cualquier integrador, logger o
// estructura interna del motor.
//
// Principios de diseño:
//   - Todas son funciones puras: mismo input → mismo output, sin efectos
//     secundarios, sin acceso a estado global.
//   - Toman const NBodySystem& — no modifican el sistema.
//   - No dependen de HierarchicalIntegrator, HybridIntegrator ni ningún
//     otro componente del motor. Solo de NBodySystem y Vec3.
//   - Adecuadas para logging estructurado, validación y benchmarks.
//
// Uso típico:
//   auto obs = compute_observables(system, t, E0, P0, L0);
//   std::cout << obs.energy_error << "\n";
//   write_csv_row(file, obs);
//
// ============================================================================
#pragma once
#include "nbody_system.h"
#include "vec3.h"
#include <cmath>
#include <string>
#include <limits>

namespace phys {

// ── Struct de observables ────────────────────────────────────────────────────

/// Snapshot completo de todos los invariantes físicos en un instante t.
/// Todos los campos son valores escalares o Vec3 — sin punteros ni
/// referencias a estado externo.
struct Observables {
    // Tiempo físico
    double t = 0.0;

    // Energía
    double E          = 0.0;   ///< Energía total: T + V
    double E_kinetic  = 0.0;   ///< Energía cinética T = Σ ½mᵢvᵢ²
    double E_potential= 0.0;   ///< Energía potencial V = -Σ Gmᵢmⱼ/rᵢⱼ
    double dE_rel     = 0.0;   ///< |E - E0| / |E0|  (error relativo)

    // Momento lineal
    Vec3   P       = {0,0,0};  ///< Momento total Σ mᵢvᵢ
    double dP_rel  = 0.0;      ///< |P - P0| / |P0|  (0 si |P0|=0)

    // Momento angular
    Vec3   L       = {0,0,0};  ///< Momento angular Σ mᵢ(rᵢ×vᵢ)
    double dL_rel  = 0.0;      ///< |L - L0| / |L0|  (0 si |L0|=0)

    // Centro de masa
    Vec3   cm_pos  = {0,0,0};  ///< Posición del CM = Σ mᵢrᵢ / M
    Vec3   cm_vel  = {0,0,0};  ///< Velocidad del CM = Σ mᵢvᵢ / M

    // Separaciones
    double sep_min = 0.0;      ///< Separación mínima entre cualquier par
    double sep_max = 0.0;      ///< Separación máxima entre cualquier par

    // Virial
    double virial  = 0.0;      ///< -2T/V (equilibrio virial = 1)

    // Metadata
    std::string regime = "";   ///< Régimen activo en este paso (opcional)
};

// ── Observables individuales ─────────────────────────────────────────────────

/// Energía cinética total: T = Σ ½ mᵢ |vᵢ|²
inline double kinetic_energy(const NBodySystem& sys) {
    double T = 0.0;
    for (const auto& b : sys.bodies)
        T += 0.5 * b.mass * b.velocity.norm2();
    return T;
}

/// Energía potencial gravitacional: V = -Σᵢ<ⱼ G mᵢ mⱼ / |rᵢⱼ|
/// Usa softening ε=0 por defecto. Pasar ε>0 para sistemas con singularidades.
inline double potential_energy(const NBodySystem& sys, double softening = 0.0) {
    double V = 0.0;
    const size_t N = sys.bodies.size();
    const double eps2 = softening * softening;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            Vec3 r = sys.bodies[j].position - sys.bodies[i].position;
            double d = std::sqrt(r.norm2() + eps2);
            V -= sys.G * sys.bodies[i].mass * sys.bodies[j].mass / d;
        }
    }
    return V;
}

/// Energía total: E = T + V
inline double total_energy(const NBodySystem& sys, double softening = 0.0) {
    return kinetic_energy(sys) + potential_energy(sys, softening);
}

/// Momento lineal total: P = Σ mᵢ vᵢ
inline Vec3 total_momentum(const NBodySystem& sys) {
    Vec3 P = {0, 0, 0};
    for (const auto& b : sys.bodies)
        P = P + b.mass * b.velocity;
    return P;
}

/// Momento angular total respecto al origen: L = Σ mᵢ (rᵢ × vᵢ)
inline Vec3 total_angular_momentum(const NBodySystem& sys) {
    Vec3 L = {0, 0, 0};
    for (const auto& b : sys.bodies) {
        const Vec3& r = b.position;
        const Vec3& v = b.velocity;
        L = L + b.mass * Vec3{
            r.y * v.z - r.z * v.y,
            r.z * v.x - r.x * v.z,
            r.x * v.y - r.y * v.x
        };
    }
    return L;
}

/// Masa total del sistema
inline double total_mass(const NBodySystem& sys) {
    double M = 0.0;
    for (const auto& b : sys.bodies) M += b.mass;
    return M;
}

/// Posición del centro de masa: r_cm = Σ mᵢ rᵢ / M
inline Vec3 center_of_mass_position(const NBodySystem& sys) {
    Vec3 cm = {0, 0, 0};
    double M = 0.0;
    for (const auto& b : sys.bodies) {
        cm = cm + b.mass * b.position;
        M += b.mass;
    }
    return (M > 0.0) ? cm / M : cm;
}

/// Velocidad del centro de masa: v_cm = Σ mᵢ vᵢ / M
inline Vec3 center_of_mass_velocity(const NBodySystem& sys) {
    Vec3 vcm = {0, 0, 0};
    double M = 0.0;
    for (const auto& b : sys.bodies) {
        vcm = vcm + b.mass * b.velocity;
        M += b.mass;
    }
    return (M > 0.0) ? vcm / M : vcm;
}

/// Separación mínima entre cualquier par de cuerpos
inline double min_separation(const NBodySystem& sys) {
    double sep = std::numeric_limits<double>::max();
    const size_t N = sys.bodies.size();
    for (size_t i = 0; i < N; ++i)
        for (size_t j = i + 1; j < N; ++j) {
            Vec3 r = sys.bodies[j].position - sys.bodies[i].position;
            sep = std::min(sep, std::sqrt(r.norm2()));
        }
    return (sep == std::numeric_limits<double>::max()) ? 0.0 : sep;
}

/// Separación máxima entre cualquier par de cuerpos
inline double max_separation(const NBodySystem& sys) {
    double sep = 0.0;
    const size_t N = sys.bodies.size();
    for (size_t i = 0; i < N; ++i)
        for (size_t j = i + 1; j < N; ++j) {
            Vec3 r = sys.bodies[j].position - sys.bodies[i].position;
            sep = std::max(sep, std::sqrt(r.norm2()));
        }
    return sep;
}

/// Ratio virial: -2T/V  (valor esperado en equilibrio: 1.0)
inline double virial_ratio(const NBodySystem& sys) {
    double V = potential_energy(sys);
    if (V == 0.0) return 0.0;
    return -2.0 * kinetic_energy(sys) / V;
}

// ── Errores relativos ─────────────────────────────────────────────────────────

/// |E - E0| / |E0|   —  0 si E0 == 0
inline double relative_energy_error(double E, double E0) {
    if (E0 == 0.0) return std::abs(E);
    return std::abs(E - E0) / std::abs(E0);
}

/// |v - v0| / |v0|  para Vec3  —  0 si |v0| == 0
inline double relative_vector_error(const Vec3& v, const Vec3& v0) {
    double n0 = std::sqrt(v0.norm2());
    if (n0 == 0.0) return std::sqrt((v - v0).norm2());
    return std::sqrt((v - v0).norm2()) / n0;
}

// ── Snapshot completo ─────────────────────────────────────────────────────────

/// Calcula todos los observables en un paso.
/// E0, P0, L0: valores de referencia (t=0) para calcular errores relativos.
/// regime: string opcional con el nombre del método activo (para logging).
inline Observables compute_observables(
    const NBodySystem& sys,
    double t,
    double E0,
    const Vec3& P0,
    const Vec3& L0,
    const std::string& regime = "",
    double softening = 0.0)
{
    Observables o;
    o.t           = t;
    o.E_kinetic   = kinetic_energy(sys);
    o.E_potential = potential_energy(sys, softening);
    o.E           = o.E_kinetic + o.E_potential;
    o.dE_rel      = relative_energy_error(o.E, E0);
    o.P           = total_momentum(sys);
    o.dP_rel      = relative_vector_error(o.P, P0);
    o.L           = total_angular_momentum(sys);
    o.dL_rel      = relative_vector_error(o.L, L0);
    o.cm_pos      = center_of_mass_position(sys);
    o.cm_vel      = center_of_mass_velocity(sys);
    o.sep_min     = min_separation(sys);
    o.sep_max     = max_separation(sys);
    o.virial      = virial_ratio(sys);
    o.regime      = regime;
    return o;
}

// ── Cabecera CSV ──────────────────────────────────────────────────────────────

/// Cabecera CSV para logging estructurado.
inline const char* csv_header() {
    return "t,E,E_kin,E_pot,dE_rel,Px,Py,Pz,dP_rel,Lx,Ly,Lz,dL_rel,"
           "cm_x,cm_y,cm_z,vcm_x,vcm_y,vcm_z,sep_min,sep_max,virial,regime";
}

/// Escribe una fila CSV con los observables (sin salto de línea al final).
inline std::string csv_row(const Observables& o) {
    char buf[1024];
    std::snprintf(buf, sizeof(buf),
        "%.10e,%.10e,%.10e,%.10e,%.6e,"
        "%.10e,%.10e,%.10e,%.6e,"
        "%.10e,%.10e,%.10e,%.6e,"
        "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,"
        "%.6e,%.6e,%.6e,%s",
        o.t, o.E, o.E_kinetic, o.E_potential, o.dE_rel,
        o.P.x, o.P.y, o.P.z, o.dP_rel,
        o.L.x, o.L.y, o.L.z, o.dL_rel,
        o.cm_pos.x, o.cm_pos.y, o.cm_pos.z,
        o.cm_vel.x, o.cm_vel.y, o.cm_vel.z,
        o.sep_min, o.sep_max, o.virial,
        o.regime.c_str()
    );
    return std::string(buf);
}

} // namespace phys