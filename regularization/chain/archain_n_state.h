// regularization/chain/archain_n_state.h
#pragma once
#include "vec3.h"
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <cmath>

/**
 * @brief Estado del integrador AR-chain para N cuerpos arbitrario.
 *
 * Generalización directa de ARChain3State a N cuerpos.
 *
 * REPRESENTACIÓN EN CADENA:
 *   La cadena está definida por un vector de índices chain_indices[N]
 *   que es una permutación de {0,...,N-1}. Los N-1 eslabones son:
 *
 *     X[k] = r[chain[k+1]] - r[chain[k]]   k = 0..N-2
 *     W[k] = v[chain[k+1]] - v[chain[k]]
 *
 *   La cadena se construye de forma que los cuerpos más cercanos sean
 *   vecinos directos (algoritmo greedy de Mikkola & Aarseth 1993).
 *   Esto minimiza el error de redondeo en las separaciones críticas.
 *
 * SEPARACIONES CRUZADAS:
 *   La separación entre cuerpos no adyacentes en la cadena se obtiene
 *   sumando eslabones intermedios:
 *     R(i,j) = X[i] + X[i+1] + ... + X[j-1]   para i < j
 *
 * FUNCIÓN DE TIEMPO (TTL):
 *   dt/ds = Ω(r)   donde   Ω = Σᵢ<ⱼ mᵢmⱼ / rᵢⱼ
 *
 * REFERENCIAS:
 *   Mikkola & Aarseth (1993), Cel. Mech. Dyn. Astron. 57, 439.
 *   Mikkola & Merritt (2006), MNRAS 372, 219.
 */
struct ARChainNState {
    int N = 0;  ///< Número de cuerpos

    // ── Coordenadas de cadena (N-1 vectores) ────────────────────────────────
    std::vector<Vec3> X;  ///< X[k] = r[chain[k+1]] - r[chain[k]]
    std::vector<Vec3> W;  ///< W[k] = v[chain[k+1]] - v[chain[k]]

    // ── Índices de la cadena ─────────────────────────────────────────────────
    std::vector<int> chain;  ///< chain[k] = índice del k-ésimo cuerpo en la cadena

    // ── Centro de masa ───────────────────────────────────────────────────────
    Vec3 cm_pos;
    Vec3 cm_vel;

    // ── Tiempo ───────────────────────────────────────────────────────────────
    double t_phys = 0.0;
    double s_fict = 0.0;

    // ── Energía ──────────────────────────────────────────────────────────────
    double energy = 0.0;

    // ── Masas ────────────────────────────────────────────────────────────────
    std::vector<double> masses;  ///< masses[i] = masa del i-ésimo cuerpo

    // ── Constructor ──────────────────────────────────────────────────────────
    ARChainNState() = default;

    explicit ARChainNState(int n)
        : N(n), X(n-1), W(n-1), chain(n), masses(n, 0.0)
    {}

    // ── Accesores ────────────────────────────────────────────────────────────
    double total_mass() const {
        double M = 0.0;
        for (double m : masses) M += m;
        return M;
    }

    bool is_valid() const {
        if (N < 2) return false;
        for (double m : masses) if (m <= 0.0) return false;
        return true;
    }

    // ── Separación acumulada entre posiciones i y j en la cadena (i < j) ───
    /**
     * @brief Suma de eslabones X[i] + ... + X[j-1] = r[chain[j]] - r[chain[i]]
     *
     * Esta operación es la clave de la estructura de cadena: la separación
     * relativa entre dos cuerpos se expresa como suma de eslabones sin
     * cancelación catastrófica.
     */
    Vec3 chain_sep(int i, int j) const {
        Vec3 r;
        for (int k = i; k < j; ++k) r = r + X[k];
        return r;
    }

    // ── Reconstrucción de posiciones absolutas ───────────────────────────────
    /**
     * @brief Reconstruye r[0..N-1] (en orden de chain) desde X[] y cm_pos.
     *
     * r[chain[0]] se fija para que el CM sea exactamente cm_pos.
     * Los demás se calculan acumulando eslabones.
     *
     * @param r_out  Vector de N posiciones en orden ORIGINAL (no de cadena).
     */
    void reconstruct_positions(std::vector<Vec3>& r_out) const {
        r_out.resize(N);
        const double M = total_mass();

        // r[chain[0]] = cm_pos - Σ_{k=1}^{N-1} (m[chain[k]] * R_0k) / M
        // donde R_0k = Σ_{j=0}^{k-1} X[j]
        Vec3 correction;
        Vec3 accum;
        for (int k = 1; k < N; ++k) {
            accum = accum + X[k-1];                       // accum = R_{0,k}
            correction = correction + accum * masses[chain[k]];
        }
        r_out[chain[0]] = cm_pos - correction * (1.0 / M);

        // r[chain[k]] = r[chain[k-1]] + X[k-1]
        for (int k = 1; k < N; ++k) {
            r_out[chain[k]] = r_out[chain[k-1]] + X[k-1];
        }
    }

    /**
     * @brief Reconstruye v[0..N-1] (en orden ORIGINAL) desde W[] y cm_vel.
     */
    void reconstruct_velocities(std::vector<Vec3>& v_out) const {
        v_out.resize(N);
        const double M = total_mass();

        Vec3 correction;
        Vec3 accum;
        for (int k = 1; k < N; ++k) {
            accum = accum + W[k-1];
            correction = correction + accum * masses[chain[k]];
        }
        v_out[chain[0]] = cm_vel - correction * (1.0 / M);

        for (int k = 1; k < N; ++k) {
            v_out[chain[k]] = v_out[chain[k-1]] + W[k-1];
        }
    }

    // ── Necesita re-construcción de cadena ──────────────────────────────────
    /**
     * @brief Comprueba si la cadena actual es subóptima y debe reconstruirse.
     *
     * La cadena debe reconstruirse cuando el par más cercano en la simulación
     * ya no es un vecino directo. Criterio conservador:
     *   sep_global_min < 0.5 * |X[k_min]|
     *
     * donde sep_global_min es la menor separación entre cualquier par,
     * y X[k_min] es el eslabón más corto de la cadena actual.
     *
     * @param r_abs  Posiciones absolutas reconstruidas.
     */
    bool needs_rechain(const std::vector<Vec3>& r_abs) const {
        // Separación mínima global (todos los pares)
        double sep_global = 1e30;
        for (int i = 0; i < N; ++i)
            for (int j = i+1; j < N; ++j) {
                double d = (r_abs[i] - r_abs[j]).norm();
                if (d < sep_global) sep_global = d;
            }

        // Eslabón más corto de la cadena actual
        double link_min = 1e30;
        for (int k = 0; k < N-1; ++k) {
            double d = X[k].norm();
            if (d < link_min) link_min = d;
        }

        return sep_global < 0.5 * link_min;
    }
};