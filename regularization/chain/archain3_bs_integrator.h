// regularization/chain/archain3_bs_integrator.h
#pragma once
#include "archain3_integrator.h"
#include <array>
#include <vector>

// ============================================================================
// INTEGRADOR BULIRSCH-STOER PARA AR-CHAIN TTL (coordenadas físicas)
//
// Implementa la extrapolación Gragg-Bulirsch-Stoer (GBS) sobre el leapfrog
// DKD del Time-Transformed Leapfrog (TTL) de Mikkola & Tanikawa (1999).
//
// PRINCIPIO:
//   El leapfrog TTL (ARChain3Integrator::leapfrog_step_bare) es time-symmetric:
//   aplicar el paso con −ds reproduce exactamente la condición inicial.
//   Esta simetría garantiza que el error de un único paso contiene únicamente
//   potencias PARES de h = ds/n. Por lo tanto, el método de paso medio
//   modificado (Gragg 1965) aplicado sobre él produce una expansión de error:
//
//     y_n − y(s+S) = α₁h² + α₂h⁴ + α₃h⁶ + ...
//
//   La extrapolación de Richardson (Neville-Aitken) elimina los términos de
//   dos en dos. Con K etapas y secuencia n = {2,4,6,...,2K} se alcanza
//   orden 2K. Con K_MAX = 8 → orden 16 teórico (límite práctico: precisión
//   de doble aritmética, ~12-14 dígitos correctos).
//
// DIFERENCIA respecto a Chain3BSIntegrator (base KS):
//   Chain3BSIntegrator opera sobre Chain3State (variables cuaternión KS Q,P).
//   Este integrador opera sobre ARChain3State (variables físicas X,V) y usa
//   directamente compute_Omega() + compute_accelerations() de la clase base.
//   No hay transformación de coordenadas ni constraint de gauge Q·P.
//
// VENTAJA clave sobre el leapfrog TTL puro (orden 2):
//   Leapfrog η=1e-4:  dE_rel ≈ 3.32%
//   GBS bs_eps=1e-10: dE_rel < 1e-10  (ganancia de 7+ órdenes de magnitud)
//   Coste adicional:  ~8-12× más evaluaciones por paso físico (pero pasos
//                     adaptativos mucho más grandes → coste total razonable)
//
// RECORTE DEL ÚLTIMO PASO:
//   A diferencia del leapfrog TTL puro (que no recorta el último paso por
//   razones de coherencia con la composición de Yoshida), este integrador
//   SÍ puede recortar el último paso ds para terminar exactamente en t_target.
//   El MMP funciona correctamente con cualquier ds > 0, incluyendo valores
//   muy pequeños cerca del objetivo.
//
// REFERENCIAS:
//   Gragg (1965), SIAM J. Numer. Anal. 2, 384 — MMP, expansión error en h²
//   Bulirsch & Stoer (1966), Numer. Math. 8, 1 — Extrapolación GBS
//   Mikkola & Tanikawa (1999), MNRAS 310, 745 — Leapfrog TTL base
//   Mikkola & Aarseth (2002), Celest. Mech. 84, 343 — GBS sobre TTL
//   Press et al. (2007), Numerical Recipes, cap. 17.3
// ============================================================================

class ARChain3BSIntegrator : public ARChain3Integrator {
public:

    // ── Parámetros de control ────────────────────────────────────────────────
    struct BSParameters {
        double bs_eps      = 1e-10;   ///< Tolerancia de error GBS (relativa)
        int    k_max       = 8;       ///< Niveles de extrapolación (≤ K_MAX)
        double initial_ds  = 1e-3;    ///< Paso ficticio inicial sugerido
        double min_ds      = 1e-14;   ///< Paso mínimo (evita bucles infinitos)
        double max_ds      = 0.5;     ///< Paso máximo (evita sub-pasos enormes)
        double safety      = 0.9;     ///< Factor de seguridad para escala de paso
        double energy_tol  = 0.01;    ///< Rechaza paso si |dE/E| > energy_tol
        int    max_steps   = 1000000; ///< Límite de pasos BS totales
        bool   verbose     = false;   ///< Diagnóstico por stdout
    };

    // ── Constructor ──────────────────────────────────────────────────────────

    /**
     * @param eta    Parámetro de paso del leapfrog TTL base (hereda de
     *               ARChain3Integrator). Usado para escalar el ds inicial:
     *               ds = eta * sep_min / Omega.
     * @param params Parámetros BS.
     */
    explicit ARChain3BSIntegrator(double eta, const BSParameters& params);

    /** @brief Constructor de conveniencia con parámetros por defecto. */
    explicit ARChain3BSIntegrator(double eta = 1e-3)
        : ARChain3BSIntegrator(eta, BSParameters{}) {}

    // ── Interfaz principal ───────────────────────────────────────────────────

    /**
     * @brief Integra desde state.t_phys hasta t_abs_final con control GBS.
     *
     * Reemplaza integrate_to() del leapfrog base con precisión de orden alto.
     * La tolerancia de error se controla con params.bs_eps.
     *
     * @param state        Estado de entrada/salida (ARChain3State).
     * @param t_abs_final  Tiempo físico absoluto objetivo.
     */
    void integrate_to_bs(ARChain3State& state, double t_abs_final);

    /**
     * @brief Un paso GBS completo en tiempo ficticio ds_total.
     *
     * Intenta avanzar state en ds_total. Si el error GBS converge a la
     * tolerancia bs_eps, acepta el paso y devuelve true. Si no converge
     * en K_MAX etapas, acepta el mejor resultado con paso reducido.
     *
     * @param state         Estado entrada/salida (modificado solo si ok).
     * @param ds_total      Paso ficticio τ a intentar.
     * @param error_out     Error GBS estimado (siempre actualizado).
     * @param ds_next_out   Sugerencia para el siguiente paso.
     * @return true si el paso fue aceptado con error < bs_eps.
     */
    bool bs_step(ARChain3State& state,
                 double         ds_total,
                 double&        error_out,
                 double&        ds_next_out);

    // ── Diagnóstico (público para tests) ────────────────────────────────────

    /** @brief Error relativo de energía respecto a E_ref. */
    double energy_error(const ARChain3State& state, double E_ref) const;

private:
    BSParameters params_;

    // Secuencia de sub-divisiones n_k = 2, 4, 6, 8, 10, 12, 14, 16
    // (Stoer & Bulirsch 1980, cap. 7; Numerical Recipes §17.3)
    static constexpr int K_MAX = 8;
    static constexpr std::array<int, K_MAX> n_seq_ = {2, 4, 6, 8, 10, 12, 14, 16};

    // ── Paso medio modificado (Gragg 1965) ───────────────────────────────────

    /**
     * @brief Avanza y0 en ds_total usando n sub-pasos de tamaño h = ds_total/n.
     *
     * El esquema del paso medio modificado es:
     *   z₀ = y0
     *   z₁ = z₀ + h·f(z₀)                       ← Euler inicial
     *   z_{m+1} = z_{m-1} + 2h·f(z_m)            ← pasos medios
     *   result = ½·(z_n + z_{n-1} + h·f(z_n))    ← paso final suavizante
     *
     * Para ARChain3State, "f(z)" no es una derivada única sino el operador
     * leapfrog DKD completo. Sin embargo, para el MMP se necesita el
     * integrador de PRIMER orden (Euler), no el leapfrog. Las ecuaciones
     * de movimiento en tiempo ficticio τ son:
     *   dX/dτ = Ω(X)·V
     *   dV/dτ = Ω(X)·A(X)
     *   dt_phys/dτ = Ω(X)
     *
     * Se evalúa Ω en el estado ACTUAL de cada sub-paso (coherente con
     * la simetría temporal del método).
     *
     * @param y0       Estado inicial.
     * @param ds_total Intervalo total de tiempo ficticio.
     * @param n        Número de sub-pasos.
     * @param result   Estado resultado.
     */
    void modified_midpoint_step(const ARChain3State& y0,
                                double               ds_total,
                                int                  n,
                                ARChain3State&       result) const;

    // ── Extrapolación polinómica de Richardson (Neville-Aitken) ─────────────

    /**
     * @brief Extrapolación de Richardson usando la expansión en h².
     *
     * Construye la tabla T[k][m] iterativamente:
     *   T[k][m] = T[k][m-1] + (T[k][m-1]-T[k-1][m-1]) / ((h_{k-m}/h_k)²-1)
     *
     * El denominador usa h² (no h) porque la expansión de error del MMP
     * sobre un integrador time-symmetric contiene únicamente potencias pares.
     * Esto duplica el orden ganado por etapa: K etapas → orden 2K.
     *
     * @param raw        Resultados del MMP para n = n₀, ..., nₖ.
     * @param step_sizes Tamaños de sub-paso h_k = ds_total/n_k.
     * @param k_current  Nivel actual (se usan raw[0..k_current]).
     * @param extrap_out Estado extrapolado de orden 2*(k_current+1).
     */
    void polynomial_extrapolation(const std::vector<ARChain3State>& raw,
                                  const std::vector<double>&        step_sizes,
                                  int                               k_current,
                                  ARChain3State&                    extrap_out) const;

    // ── Estimación de error ──────────────────────────────────────────────────

    /**
     * @brief Error relativo entre dos extrapolaciones consecutivas.
     *
     * Mide la norma máxima escalada de las diferencias en X1,V1,X2,V2.
     * La escala evita divisiones por cero cerca del origen.
     */
    double estimate_error(const ARChain3State& high_order,
                          const ARChain3State& low_order) const;

    // ── Control del paso ─────────────────────────────────────────────────────

    /**
     * @brief Factor de escala para el siguiente ds, dado el error actual.
     *
     * Basado en la fórmula estándar GBS (Numerical Recipes §17.3):
     *   scale = safety · (1/error)^(1/(2*k+1))
     * Limitado al intervalo [0.1, 4.0].
     */
    double step_scale_factor(double error, int k_opt) const;
};