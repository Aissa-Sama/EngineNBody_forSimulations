// regularization/chain/archain3_state.h
#pragma once
#include "vec3.h"

/**
 * @brief Estado del integrador AR-chain (regularización algorítmica) para 3 cuerpos.
 *
 * A diferencia de Chain3State (variables KS), este estado trabaja directamente
 * en coordenadas físicas. No hay transformación de coordenadas, no hay gauge,
 * no hay constraint Q·P.
 *
 * Referencia: Mikkola & Tanikawa (1999), MNRAS 310, 745-749.
 *
 * REPRESENTACIÓN EN CADENA:
 *   Se almacenan las DIFERENCIAS de posición y velocidad (coordenadas de cadena)
 *   para reducir el error de cancelación en encuentros muy cercanos:
 *
 *     X1 = r2 - r1   (separación eslabón 1)
 *     X2 = r3 - r2   (separación eslabón 2)
 *     V1 = v2 - v1   (velocidad relativa eslabón 1)
 *     V2 = v3 - v2   (velocidad relativa eslabón 2)
 *
 *   Las posiciones absolutas se reconstruyen desde X1, X2 y el CM.
 *
 * FUNCIÓN DE TIEMPO (TTL):
 *   El tiempo ficticio s se relaciona con el tiempo físico t mediante:
 *     dt/ds = Ω(r)   donde   Ω = Σᵢ<ⱼ mᵢmⱼ / rᵢⱼ
 *   Cuando rᵢⱼ → 0, Ω → ∞, lo que hace que dt → 0 para ds fijo
 *   (el paso físico se vuelve pequeño automáticamente en encuentros cercanos).
 */
struct ARChain3State {
    // ── Coordenadas de cadena ────────────────────────────────────────────────
    Vec3 X1;   ///< r2 - r1  (separación físicia, eslabón 1)
    Vec3 X2;   ///< r3 - r2  (separación física, eslabón 2)
    Vec3 V1;   ///< v2 - v1  (velocidad relativa, eslabón 1)
    Vec3 V2;   ///< v3 - v2  (velocidad relativa, eslabón 2)

    // ── Centro de masa ───────────────────────────────────────────────────────
    Vec3 cm_pos;   ///< Posición del CM (constante, se propaga linealmente)
    Vec3 cm_vel;   ///< Velocidad del CM (constante)

    // ── Tiempo ───────────────────────────────────────────────────────────────
    double t_phys;   ///< Tiempo físico acumulado
    double s_fict;   ///< Tiempo ficticio acumulado

    // ── Energía (invariante de movimiento) ──────────────────────────────────
    double energy;   ///< Energía total del subsistema de 3 cuerpos (T + U, G=1)

    // ── Masas ────────────────────────────────────────────────────────────────
    double masses[3];   ///< masses[0]=m1, masses[1]=m2, masses[2]=m3

    // ── Constructor por defecto ──────────────────────────────────────────────
    ARChain3State()
        : X1(), X2(), V1(), V2()
        , cm_pos(), cm_vel()
        , t_phys(0.0), s_fict(0.0)
        , energy(0.0)
    {
        masses[0] = masses[1] = masses[2] = 0.0;
    }

    // ── Accesores de masa ────────────────────────────────────────────────────
    double m1() const { return masses[0]; }
    double m2() const { return masses[1]; }
    double m3() const { return masses[2]; }
    double total_mass() const { return masses[0] + masses[1] + masses[2]; }

    // ── Validación ───────────────────────────────────────────────────────────
    bool is_valid() const {
        return (m1() > 0.0 && m2() > 0.0 && m3() > 0.0);
    }

    // ── Reconstrucción de posiciones absolutas ───────────────────────────────
    /**
     * @brief Reconstruye r1, r2, r3 desde X1, X2 y cm_pos.
     *
     * Dado que X1 = r2 - r1 y X2 = r3 - r2, y el CM es:
     *   r_cm = (m1·r1 + m2·r2 + m3·r3) / M
     *
     * Se obtiene r1 desde:
     *   r_cm = (m1·r1 + m2·(r1+X1) + m3·(r1+X1+X2)) / M
     *         = r1 + (m2·X1 + m3·(X1+X2)) / M
     *
     *   → r1 = r_cm - (m2·X1 + m3·(X1+X2)) / M
     */
    void reconstruct_positions(Vec3& r1, Vec3& r2, Vec3& r3) const {
        double M = total_mass();
        Vec3 X3 = X1 + X2;   // r3 - r1
        r1 = cm_pos - (X1 * m2() + X3 * m3()) * (1.0 / M);
        r2 = r1 + X1;
        r3 = r2 + X2;
    }

    /**
     * @brief Reconstruye v1, v2, v3 desde V1, V2 y cm_vel.
     * Misma lógica que reconstruct_positions pero para velocidades.
     */
    void reconstruct_velocities(Vec3& v1, Vec3& v2, Vec3& v3) const {
        double M = total_mass();
        Vec3 V3 = V1 + V2;   // v3 - v1
        v1 = cm_vel - (V1 * m2() + V3 * m3()) * (1.0 / M);
        v2 = v1 + V1;
        v3 = v2 + V2;
    }
};