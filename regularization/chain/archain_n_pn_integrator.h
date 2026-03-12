// regularization/chain/archain_n_pn_integrator.h
#pragma once
#include "archain_n_integrator.h"
#include "pn_forces.h"

/**
 * @brief Integrador AR-chain TTL para N cuerpos con correcciones post-Newtonianas.
 *
 * Extiende ARChainNIntegrator añadiendo correcciones PN1, PN2 y PN2.5 a las
 * ecuaciones de movimiento. Las correcciones se añaden en compute_accelerations()
 * sobreescribiendo el método virtual.
 *
 * ── COMPATIBILIDAD CON GBS ────────────────────────────────────────────────
 *
 * PN1 y PN2 son conservativas → el leapfrog TTL sigue siendo simpléctico
 * (aproximadamente) si se incluyen en el KICK. Sin embargo, la precisión
 * de orden 2 del leapfrog no es suficiente para estos términos en sistemas
 * relativistas. Se recomienda usar ARChainNPNBSIntegrator (GBS sobre este).
 *
 * PN2.5 es DISIPATIVO → el leapfrog simpléctico NO es correcto. Solo usar
 * con el integrador GBS (MMP + extrapolación de Richardson). Mikkola &
 * Merritt (2008) §4 discuten explícitamente que el MMP maneja correctamente
 * los términos dependientes de velocidad (PN25) porque no asume Hamiltoniano.
 *
 * ── USO DE UNIDADES ──────────────────────────────────────────────────────
 *
 * G=1 en todo el código. c es un parámetro explícito en unidades N-body.
 * Para conversión física:
 *   c_SI = 2.998×10⁸ m/s
 *   En unidades N-body (M_sol, pc, km/s):
 *     c_Nbody ≈ 2.998×10⁵ km/s / v_unit
 *   Para v_unit ~ 1 km/s: c_Nbody ≈ 3×10⁵
 *   Para binarias compactas (v~0.1c): elegir c~10 para hacer PN visible
 *
 * ── ESTIMACIÓN DEL PARÁMETRO PN ─────────────────────────────────────────
 *
 * La corrección PN1 es relevante cuando:
 *   ε_PN = m/(r·c²) ~ v²/c² ≥ 10⁻⁶  (1 ppm de la Newtoniana)
 *
 * Para separaciones sub-mpc con masas ~10⁶ M_sol:
 *   ε_PN ~ 10⁻⁴ a 10⁻² → PN1 y PN2 son importantes.
 *
 * PN2.5 produce inspiral orbital:
 *   τ_merge ~ (5/256) × c⁵ × a⁴ / (G³ m₁ m₂ (m₁+m₂))  [Peters 1964]
 *
 * REFERENCIAS:
 *   Mikkola & Merritt (2008), AJ 135, 2398 — AR-CHAIN con PN.
 *   Harfst et al. (2008), MNRAS 389, 2 — implementación PN en AR-CHAIN.
 *   Einstein, Infeld & Hoffmann (1938) — ecuaciones EIH.
 *   Burke (1971), JMP 12, 401 — reacción de radiación.
 *   Peters (1964), Phys. Rev. 136, B1224 — tiempo de merge PN2.5.
 */
class ARChainNPNIntegrator : public ARChainNIntegrator {
public:

    /**
     * @param eta       Parámetro de paso TTL.
     * @param c_speed   Velocidad de la luz en unidades N-body (G=1).
     * @param pn_order  Bitmask de órdenes PN activos:
     *                    1 = PN1 (EIH, O(1/c²))
     *                    2 = PN2 (O(1/c⁴))
     *                    4 = PN25 (reacción radiación O(1/c⁵)) — solo con GBS
     *                    3 = PN1+PN2,   7 = PN1+PN2+PN25
     */
    explicit ARChainNPNIntegrator(double eta    = 1e-3,
                                  double c_speed = 1e4,
                                  int    pn_order = 1);

    double c_speed() const { return c_; }
    int    pn_order() const { return pn_order_; }

    // ── Observables (públicos para tests) ───────────────────────────────────

    /** @brief Energía post-Newtoniana PN1 (aproximada, conservativa). */
    double compute_energy_pn1(const ARChainNState& state) const;

protected:
    double c_;
    int    pn_order_;

    /**
     * @brief Override: añade correcciones PN a las aceleraciones Newtonianas.
     *
     * A[k] = A_N[k] + A_PN[k]
     * donde A_PN[k] = a_pn[chain[k+1]] - a_pn[chain[k]]
     */
    void compute_accelerations(const ARChainNState& state,
                                const std::vector<Vec3>& r_abs,
                                std::vector<Vec3>& A_out) const;
};
