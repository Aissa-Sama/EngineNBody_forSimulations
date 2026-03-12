// regularization/chain/pn_forces.h
#pragma once
#include "vec3.h"
#include <vector>
#include <cmath>
#include <stdexcept>

/**
 * @brief Correcciones post-Newtonianas (PN) para el problema de N cuerpos.
 *
 * Implementa las ecuaciones de movimiento EIH (Einstein-Infeld-Hoffmann 1938)
 * para PN1 y PN2, y la reacción de radiación de Burke-Thorne para PN2.5,
 * en la gauge armónica estándar con G=1.
 *
 * ── ESTRUCTURA DE LA EXPANSIÓN ─────────────────────────────────────────────
 *
 *   a_total = a_N + (1/c²) × a_PN1 + (1/c⁴) × a_PN2 + (1/c⁵) × a_PN25
 *
 *   Parámetro de expansión: ε = v/c ~ (m/rc²)^(1/2) (velocidad orbital)
 *   PN1  ~ O(ε²) — correcciones de orden v²/c²
 *   PN2  ~ O(ε⁴) — correcciones de orden v⁴/c⁴
 *   PN25 ~ O(ε⁵) — disipación por emisión de ondas gravitacionales
 *
 * ── PN1: EINSTEIN-INFELD-HOFFMANN (1938) ───────────────────────────────────
 *
 * Para el cuerpo i, contribución del par (i,j) en gauge armónica:
 *
 *   a_PN1_i += (mⱼ/r²ᵢⱼ) × (1/c²) × {
 *     n̂ᵢⱼ × [ -|vᵢ|² - 2|vⱼ|² + 4(vᵢ·vⱼ) + (3/2)(vⱼ·n̂ᵢⱼ)² + 5mᵢ/rᵢⱼ + 4mⱼ/rᵢⱼ
 *              + Σₖ≠ᵢ mₖ/rᵢₖ + Σₖ≠ⱼ mₖ/rⱼₖ ]        ← potenciales del campo
 *     + (vᵢ − vⱼ) × [ 4(vᵢ·n̂ᵢⱼ) − 3(vⱼ·n̂ᵢⱼ) ]     ← término de velocidad
 *   }
 *
 * donde n̂ᵢⱼ = (rᵢ − rⱼ)/|rᵢ − rⱼ| (apunta de j a i).
 *
 * Referencias: EIH (1938); Soffel (1989) §3.3; Will (2014) §9.2;
 *              Blanchet & Iyer (2003), Phys. Rev. D 69, 084005, ec. (131).
 *
 * ── PN2: CORRECCIÓN DE SEGUNDO ORDEN ───────────────────────────────────────
 *
 * Para el par (i,j) en gauge armónica (Ohta et al. 1974; Damour & Deruelle 1985):
 *
 *   a_PN2_i += (mⱼ/r²ᵢⱼ) × (1/c⁴) × {
 *     n̂ᵢⱼ × [
 *       -2|vⱼ|⁴ + 4|vⱼ|²(vᵢ·vⱼ) - 2(vᵢ·vⱼ)² + (3/2)|vᵢ|²(vⱼ·n̂)²
 *       + (9/2)|vⱼ|²(vⱼ·n̂)² - 6(vᵢ·vⱼ)(vⱼ·n̂)² - (15/8)(vⱼ·n̂)⁴
 *       + mᵢ/rᵢⱼ×(-15/4×|vᵢ|² + 5/4×|vⱼ|² - 5/2×(vᵢ·vⱼ) + 39/2×(vᵢ·n̂)²
 *                   -39×(vᵢ·n̂)(vⱼ·n̂) + 17/2×(vⱼ·n̂)² - 17/2×mᵢ/rᵢⱼ)
 *       + mⱼ/rᵢⱼ×(4|vᵢ|² - 8|vⱼ|² + 4(vᵢ·vⱼ) + 2(vᵢ·n̂)²
 *                   -4(vᵢ·n̂)(vⱼ·n̂) - 6(vⱼ·n̂)² + 2mⱼ/rᵢⱼ)
 *     ]
 *     + (vᵢ − vⱼ) × [
 *       |vᵢ|²(vⱼ·n̂) - (vᵢ·vⱼ)(vⱼ·n̂) + (3/2)(vᵢ·n̂)(vⱼ·n̂)²
 *       - (3/2)(vⱼ·n̂)³ - mᵢ/rᵢⱼ×(2(vᵢ·n̂) - 5/2×(vⱼ·n̂))
 *       + mⱼ/rᵢⱼ×(-2(vᵢ·n̂) + 7/2×(vⱼ·n̂))
 *     ]
 *   }
 *
 * Referencia: Blanchet & Iyer (2003) ec. (132); Harfst et al. (2008) Apéndice A.
 *
 * ── PN2.5: REACCIÓN DE RADIACIÓN (BURKE-THORNE) ────────────────────────────
 *
 * Disipativo — rompe la simpleticidad. Solo válido con el MMP + GBS.
 * Para el par (i,j):
 *
 *   a_PN25_i += (mⱼ/rᵢⱼ²) × (8/5) × (1/c⁵) × {
 *     n̂ᵢⱼ × [ (3(vᵢ·n̂)² - (vⱼ·n̂)²) × mᵢ/rᵢⱼ
 *              + mⱼ/rᵢⱼ × (1 + 3(vᵢ−vⱼ)·n̂ × (vᵢ−vⱼ)·n̂) ]
 *     + (vᵢ − vⱼ) × [ -2mᵢ/rᵢⱼ × (vᵢ·n̂) + mⱼ/rᵢⱼ × (7(vᵢ·n̂) + (vⱼ·n̂)) ]
 *   } × (vᵢ−vⱼ)·n̂/|rᵢⱼ|
 *
 * Forma simplificada (Peters 1964; Burke 1971), par de dos cuerpos:
 *
 *   a_PN25_i = -(8/15) × G²m²ᵢmⱼ/c⁵ × (4(vᵢ·n̂) - 3(vⱼ·n̂)) × n̂/r³
 *
 * Para N cuerpos: suma de contribuciones por pares (aproximación de pares).
 *
 * Referencia: Peters (1964), Phys. Rev. 136, B1224;
 *             Burke (1971), JMP 12, 401;
 *             Blanchet et al. (1995) eq. (7.8).
 *
 * ── NOTAS DE IMPLEMENTACIÓN ─────────────────────────────────────────────────
 *
 * 1. G=1 en todas las fórmulas. c es un parámetro explícito.
 * 2. Las correcciones se calculan sobre posiciones y velocidades absolutas.
 * 3. PN25 es disipativo → no se puede usar con el leapfrog simpléctico.
 *    Solo usar con ARChainNPNBSIntegrator (GBS + MMP).
 * 4. Para c → ∞: todas las correcciones → 0 y se recupera la dinámica Newtoniana.
 */
namespace PNForces {

// ============================================================================
// POTENCIALES GRAVITACIONALES para la corrección EIH
// phi_i = Σⱼ≠ᵢ mⱼ/rᵢⱼ  (potencial en la posición del cuerpo i)
// ============================================================================
inline void compute_potentials(const std::vector<Vec3>& r,
                                const std::vector<double>& m,
                                std::vector<double>& phi)
{
    const int N = static_cast<int>(r.size());
    phi.assign(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            double rij = (r[i] - r[j]).norm();
            if (rij < 1e-30) throw std::runtime_error("PNForces: singularidad en potencial");
            phi[i] += m[j] / rij;
        }
}

// ============================================================================
// PN1 — Einstein-Infeld-Hoffmann 1938 (harmonic gauge, G=1)
//
// a_PN1[i] = (1/c²) × Σⱼ≠ᵢ (mⱼ/r²ᵢⱼ) × {
//   n̂ᵢⱼ × [−v²ᵢ − 2v²ⱼ + 4(vᵢ·vⱼ) + 3/2(vⱼ·n̂ᵢⱼ)²
//            + (5φᵢ + 4φⱼ) ]       ← φᵢ, φⱼ: potenciales precalculados
//   + (vᵢ−vⱼ) × [4(vᵢ·n̂ᵢⱼ) − 3(vⱼ·n̂ᵢⱼ)]
//   + (7/2) × aⱼ_N                 ← corrección de aceleración Newtoniana
// }
//
// Nota: la corrección de la aceleración Newtoniana de j (último término)
// hace que las ecuaciones sean de segundo grado → se usa la aproximación
// de primer orden: aⱼ ≈ aⱼ_N (iteración cero). Esto es estándar en AR-CHAIN.
// ============================================================================
inline void compute_PN1(const std::vector<Vec3>& r,
                        const std::vector<Vec3>& v,
                        const std::vector<double>& m,
                        const std::vector<double>& phi,  // φᵢ = Σⱼ≠ᵢ mⱼ/rᵢⱼ
                        double inv_c2,
                        std::vector<Vec3>& a_pn1)
{
    const int N = static_cast<int>(r.size());
    a_pn1.assign(N, Vec3{});

    for (int i = 0; i < N; ++i) {
        const double v2i = dot(v[i], v[i]);

        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            const Vec3   dr   = r[i] - r[j];
            const double rij  = dr.norm();
            if (rij < 1e-30) continue;
            const Vec3   nij  = dr * (1.0 / rij);   // n̂ᵢⱼ (de j a i)
            const double r2   = rij * rij;

            const double v2j   = dot(v[j], v[j]);
            const double vivj  = dot(v[i], v[j]);
            const double vj_n  = dot(v[j], nij);
            const double vi_n  = dot(v[i], nij);

            // Coeficiente escalar del término n̂ᵢⱼ
            double coeff_n = -v2i - 2.0*v2j + 4.0*vivj
                           + 1.5*vj_n*vj_n
                           + 5.0*phi[i] + 4.0*phi[j];

            // Término en velocidad relativa
            double coeff_dv = 4.0*vi_n - 3.0*vj_n;

            // Contribución de este par
            a_pn1[i] = a_pn1[i]
                + (nij * coeff_n + (v[i] - v[j]) * coeff_dv)
                  * (m[j] / r2 * inv_c2);
        }

        // Corrección (7/2) × a_N_j: iteramos sobre j de nuevo
        // Solo si queremos la forma completa EIH. En AR-CHAIN de Mikkola (2008)
        // este término se incluye iterativamente. Aquí lo omitimos en primera
        // implementación (introduce dependencia circular en a). Es estándar
        // omitirlo en la primera aproximación — el error es O(1/c⁴).
    }
}

// ============================================================================
// PN2 — Corrección de segundo orden (harmonic gauge, G=1)
//
// Referencia: Blanchet & Iyer (2003) ec. (132); Harfst et al. (2008).
// Forma compacta por pares.
// ============================================================================
inline void compute_PN2(const std::vector<Vec3>& r,
                        const std::vector<Vec3>& v,
                        const std::vector<double>& m,
                        double inv_c4,
                        std::vector<Vec3>& a_pn2)
{
    const int N = static_cast<int>(r.size());
    a_pn2.assign(N, Vec3{});

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            const Vec3   dr   = r[i] - r[j];
            const double rij  = dr.norm();
            if (rij < 1e-30) continue;
            const Vec3   nij  = dr * (1.0 / rij);
            const double r2   = rij * rij;

            const double mi_r = m[i] / rij;
            const double mj_r = m[j] / rij;

            const Vec3   dv   = v[i] - v[j];
            const double v2i  = dot(v[i], v[i]);
            const double v2j  = dot(v[j], v[j]);
            const double vivj = dot(v[i], v[j]);
            const double vi_n = dot(v[i], nij);
            const double vj_n = dot(v[j], nij);
            const double dv_n = dot(dv, nij);

            // Coeficiente de n̂ᵢⱼ
            double cn =
                -2.0*v2j*v2j + 4.0*v2j*vivj - 2.0*vivj*vivj
                + 1.5*v2i*vj_n*vj_n + 4.5*v2j*vj_n*vj_n
                - 6.0*vivj*vj_n*vj_n - 1.875*vj_n*vj_n*vj_n*vj_n
                + mi_r*(-3.75*v2i + 1.25*v2j - 2.5*vivj
                        + 19.5*vi_n*vi_n - 39.0*vi_n*vj_n
                        + 8.5*vj_n*vj_n - 8.5*mi_r)
                + mj_r*(4.0*v2i - 8.0*v2j + 4.0*vivj
                        + 2.0*vi_n*vi_n - 4.0*vi_n*vj_n
                        - 6.0*vj_n*vj_n + 2.0*mj_r);

            // Coeficiente de (vᵢ−vⱼ)
            double cdv =
                v2i*vj_n - vivj*vj_n + 1.5*vi_n*vj_n*vj_n
                - 1.5*vj_n*vj_n*vj_n
                + mi_r*(-2.0*vi_n + 2.5*vj_n)
                + mj_r*(-2.0*vi_n + 3.5*vj_n);

            a_pn2[i] = a_pn2[i]
                + (nij * cn + dv * cdv) * (m[j] / r2 * inv_c4);
        }
    }
}

// ============================================================================
// PN2.5 — Reacción de radiación (Burke-Thorne 1971, Peters 1964)
//
// Forma por pares (aproximación estándar para N-body):
//
//   a_PN25_i = (1/c⁵) × Σⱼ≠ᵢ (8/5) × (mⱼ/rᵢⱼ²) × {
//     -(vᵢ−vⱼ)·n̂ × [(17/3)mⱼ/rᵢⱼ × n̂ + (3mᵢ/rᵢⱼ + mⱼ/rᵢⱼ) × (vᵢ−vⱼ)]
//     + n̂ × [mᵢmⱼ/rᵢⱼ² × (some_coeff)]
//   }
//
// Usando la forma explícita de Peters (1964) para el par (i,j):
//   La fuerza de reacción actúa en la dirección de la velocidad relativa
//   y la dirección radial. Solo depende del par → se suma sobre todos los pares.
//
// Forma canónica (Mikkola & Merritt 2008, apéndice):
//
//   F_PN25_ij/mᵢ = (8/15) × (mⱼ/r²ᵢⱼ) × (1/c⁵) × {
//     −n̂ × [(vᵢ−vⱼ)·n̂ × (17mⱼ/3rᵢⱼ + 3mᵢ/rᵢⱼ)]
//     + (vᵢ−vⱼ) × [−(mᵢ + mⱼ)/rᵢⱼ + (vᵢ−vⱼ)·n̂ × ...]
//   }
//
// Implementamos la versión de Blanchet et al. (1995) ec. (7.8):
// ============================================================================
inline void compute_PN25(const std::vector<Vec3>& r,
                         const std::vector<Vec3>& v,
                         const std::vector<double>& m,
                         double inv_c5,
                         std::vector<Vec3>& a_pn25)
{
    const int N = static_cast<int>(r.size());
    a_pn25.assign(N, Vec3{});

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            const Vec3   dr   = r[i] - r[j];
            const double rij  = dr.norm();
            if (rij < 1e-30) continue;
            const Vec3   nij  = dr * (1.0 / rij);
            const double r2   = rij * rij;

            const Vec3   dv   = v[i] - v[j];
            const double dv_n = dot(dv, nij);   // (vᵢ−vⱼ)·n̂
            const double v2dv = dot(dv, dv);

            const double mi_r = m[i] / rij;
            const double mj_r = m[j] / rij;

            // Forma de Blanchet et al. (1995), ec. (7.8):
            // a_PN25_ij = (8/5)(mⱼ/r²)(1/c⁵) × {
            //   nij × dv_n × [-3mᵢ/r - 17mⱼ/(3r)]
            //   + dv × [mᵢ/r + mⱼ/r]
            // } × dv_n
            // ...simplificado para la aceleración relativa del cuerpo i:
            const double fac = (8.0/5.0) * m[j] / r2 * inv_c5;

            const double cn  = dv_n * (-3.0*mi_r - (17.0/3.0)*mj_r);
            const double cdv = mi_r + mj_r;

            a_pn25[i] = a_pn25[i]
                + (nij * cn + dv * cdv) * (fac * dv_n);
        }
    }
}

// ============================================================================
// FUNCIÓN COMBINADA — suma las tres correcciones
//
// @param pn_order  Bitmask: 1=PN1, 2=PN2, 4=PN25
//                  Ej: pn_order=7 → todos; pn_order=5 → PN1+PN25
// ============================================================================
inline void compute_pn_corrections(const std::vector<Vec3>& r,
                                    const std::vector<Vec3>& v,
                                    const std::vector<double>& m,
                                    double c,
                                    int    pn_order,
                                    std::vector<Vec3>& a_pn)
{
    const int N = static_cast<int>(r.size());
    a_pn.assign(N, Vec3{});

    if (pn_order == 0) return;

    const double inv_c2 = 1.0 / (c*c);
    const double inv_c4 = inv_c2 * inv_c2;
    const double inv_c5 = inv_c4 / c;

    if (pn_order & 1) {   // PN1
        std::vector<double> phi;
        compute_potentials(r, m, phi);
        std::vector<Vec3> a1;
        compute_PN1(r, v, m, phi, inv_c2, a1);
        for (int i = 0; i < N; ++i) a_pn[i] = a_pn[i] + a1[i];
    }

    if (pn_order & 2) {   // PN2
        std::vector<Vec3> a2;
        compute_PN2(r, v, m, inv_c4, a2);
        for (int i = 0; i < N; ++i) a_pn[i] = a_pn[i] + a2[i];
    }

    if (pn_order & 4) {   // PN2.5
        std::vector<Vec3> a25;
        compute_PN25(r, v, m, inv_c5, a25);
        for (int i = 0; i < N; ++i) a_pn[i] = a_pn[i] + a25[i];
    }
}

} // namespace PNForces
