// regularization/chain/archain_n_integrator.h
#pragma once
#include "archain_n_state.h"
#include "nbody_system.h"

/**
 * @brief Integrador AR-chain TTL para N cuerpos arbitrario.
 *
 * Generalización de ARChain3Integrator a N cuerpos usando:
 *   - Coordenadas de cadena con índices óptimos (FIND_CHAIN_INDICES)
 *   - Leapfrog DKD en tiempo ficticio s con función de tiempo Ω
 *   - Re-construcción automática de cadena cuando cambia la configuración
 *   - Extrapolación GBS disponible vía ARChainNBSIntegrator (subclase)
 *
 * FORMALISMO TTL PARA N ARBITRARIO (Mikkola & Merritt 2006, §4, caso β=1):
 *
 *   Ω = Σᵢ<ⱼ mᵢmⱼ / rᵢⱼ   (todos los N(N-1)/2 pares)
 *
 *   K-paso:  drₖ/ds = vₖ / Ω   →   Xₖ += Wₖ · (h·Ω)
 *   D-paso:  dvₖ/ds = aₖ / Ω   →   Wₖ += Aₖ · (h·Ω)
 *   Tiempo:  dt/ds = 1/Ω        →   t  += h·Ω
 *
 *   donde Aₖ = a[chain[k+1]] - a[chain[k]] son las diferencias de
 *   aceleración en coordenadas de cadena.
 *
 * REFERENCIAS:
 *   Mikkola & Aarseth (1993), Cel. Mech. Dyn. Astron. 57, 439.
 *   Mikkola & Merritt (2006), MNRAS 372, 219.
 *   Mikkola & Merritt (2008), AJ 135, 2398.
 */
class ARChainNIntegrator {
public:
    explicit ARChainNIntegrator(double eta = 1e-3);
    virtual ~ARChainNIntegrator() = default;

    // ── Interfaz principal ───────────────────────────────────────────────────

    /**
     * @brief Inicializa el estado desde el sistema físico.
     *
     * Construye la cadena óptima (FIND_CHAIN_INDICES), calcula X[], W[],
     * cm_pos, cm_vel, masas y energía inicial.
     *
     * @param system   Sistema N-body.
     * @param indices  Índices de los cuerpos a incluir en el AR-chain.
     */
    ARChainNState initialize(const NBodySystem& system,
                             const std::vector<int>& indices) const;

    /**
     * @brief Integra desde state.t_phys hasta t_abs_final.
     *
     * Modo directo sin cortes externos. El AR-chain controla su propio tiempo.
     * Re-construye la cadena automáticamente si la configuración cambia.
     */
    void integrate_to(ARChainNState& state, double t_abs_final) const;

    /**
     * @brief Integra un intervalo dt_phys (modo bloque, legacy).
     */
    void integrate(ARChainNState& state, double dt_phys) const;

    /**
     * @brief Escribe posiciones y velocidades de vuelta al sistema físico.
     */
    void write_back(const ARChainNState& state, NBodySystem& system,
                    const std::vector<int>& indices) const;

    // ── Observables (públicos para tests y subclases) ────────────────────────

    /** @brief Ω = Σᵢ<ⱼ mᵢmⱼ/rᵢⱼ usando separaciones de cadena. */
    double compute_Omega(const ARChainNState& state) const;

    /** @brief Energía total T + U. */
    double compute_energy(const ARChainNState& state) const;

    /** @brief Separación mínima entre cualquier par de cuerpos. */
    double compute_sep_min(const ARChainNState& state) const;

    // ── Algoritmo de cadena ──────────────────────────────────────────────────

    /**
     * @brief Construye la permutación óptima de la cadena (greedy, O(N²)).
     *
     * Algoritmo de Mikkola & Aarseth (1993):
     *   1. Empieza con el par más cercano.
     *   2. En cada paso, extiende por el extremo más cercano al nodo libre
     *      más próximo.
     *   3. Garantiza que los cuerpos más cercanos son vecinos directos.
     *
     * @param r     Posiciones absolutas de los N cuerpos.
     * @param chain Salida: permutación óptima de índices {0..N-1}.
     */
    static void find_chain_indices(const std::vector<Vec3>& r,
                                   std::vector<int>& chain);

    /**
     * @brief Re-construye la cadena desde posiciones absolutas actuales.
     *
     * Llama a find_chain_indices y actualiza X[], W[] en state.
     */
    static void rebuild_chain(ARChainNState& state,
                               const std::vector<Vec3>& r_abs,
                               const std::vector<Vec3>& v_abs);

protected:
    double eta_;
    double ds_min_;
    double ds_max_;

    /**
     * @brief Diferencias de aceleración en coordenadas de cadena.
     *
     * A[k] = a[chain[k+1]] - a[chain[k]]
     *
     * Las aceleraciones absolutas se calculan desde las separaciones
     * reconstruidas por chain_sep(). La estructura de cadena garantiza
     * que los pares más cercanos (los más sensibles numéricamente) se
     * calculan directamente desde X[k] sin acumulación de error.
     *
     * @param state   Estado AR-chain.
     * @param r_abs   Posiciones absolutas (precalculadas para eficiencia).
     * @param A_out   Salida: N-1 diferencias de aceleración.
     */
    void compute_accelerations(const ARChainNState& state,
                                const std::vector<Vec3>& r_abs,
                                std::vector<Vec3>& A_out) const;

    /**
     * @brief Un paso DKD (leapfrog) en tiempo ficticio. Acepta ds < 0.
     */
    void leapfrog_step_bare(ARChainNState& state, double ds) const;

    void leapfrog_step(ARChainNState& state, double ds) const;

    /**
     * @brief Núcleo del bucle de integración.
     *
     * Avanza state desde state.t_phys hasta t_target con pasos adaptativos
     * ds = η·sep_min/Ω. Re-construye la cadena si needs_rechain().
     */
    void integrate_until(ARChainNState& state, double t_target) const;
};