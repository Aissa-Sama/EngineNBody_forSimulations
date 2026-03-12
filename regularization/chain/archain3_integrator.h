// regularization/chain/archain3_integrator.h
#pragma once
#include "archain3_state.h"
#include "nbody_system.h"

/**
 * @brief Integrador AR-chain (regularización algorítmica) para 3 cuerpos.
 *
 * Implementa el método TTL (Time-Transformed Leapfrog) de
 * Mikkola & Tanikawa (1999), MNRAS 310, 745-749.
 *
 * ── PRINCIPIO DEL TTL ────────────────────────────────────────────────────────
 * El tiempo ficticio s se relaciona con el tiempo físico t mediante:
 *   dt/ds = Ω(r)   donde   Ω = Σᵢ<ⱼ mᵢmⱼ / rᵢⱼ
 *
 * El leapfrog resultante es simpléctico: conserva una energía sombra
 * acotada para todo t. La energía física H oscila con amplitud O(ds²)
 * pero sin drift secular.
 *
 * ── CONTROL DE PASO ADAPTATIVO ───────────────────────────────────────────────
 * El parámetro central es eta (η), que controla el paso en tiempo ficticio:
 *
 *   ds_step = η · sep_min / Ω
 *
 * donde sep_min = min(r₁₂, r₂₃, r₁₃) es la separación mínima actual.
 *
 * COMPORTAMIENTO:
 *   • Lejos del encuentro (sep_min grande):  ds_step grande → pasos eficientes
 *   • Cerca del encuentro (sep_min → 0):     ds_step → 0   → paso físico se reduce
 *
 * El paso físico resultante es: dt = ds_step · Ω = η · sep_min
 * Esto replica el comportamiento del RK4 adaptativo con dt ∝ sep_min,
 * que es la estrategia que permite cruzar encuentros cercanos con precisión.
 *
 * ── DOS MODOS DE INTEGRACIÓN ─────────────────────────────────────────────────
 *
 * 1. integrate(state, dt_phys) — modo bloque (legacy / N > 3):
 *    Integra un intervalo dt_phys de tiempo físico. Para uso dentro del
 *    integrador jerárquico cuando el AR-chain es un subsistema parcial y
 *    el integrador exterior necesita sincronizar el tiempo.
 *
 * 2. integrate_to(state, t_abs_final) — modo directo (N=3, sistema completo):
 *    Integra desde state.t_phys hasta t_abs_final SIN cortes externos.
 *    El AR-chain controla su propio tiempo de extremo a extremo.
 *    Este es el modo correcto cuando el triple es el único subsistema activo:
 *    elimina los errores de fase introducidos por cortes de dt fijo que
 *    acumulan error antes del encuentro cercano y alteran la trayectoria.
 *
 *    Referencia: diagnóstico arquitectural documentado en
 *    "NBody_Maestro_v4" (Marzo 2026) — §5: el error residual dE=3–7%
 *    no es numérico sino arquitectural; proviene de los cortes de dt.
 *
 * ── ELECCIÓN DE η ────────────────────────────────────────────────────────────
 *   η = 1e-3  →  ~8k pasos para figura-8,  max|ΔE/E₀| ~ 0.2  (rápido)
 *   η = 1e-4  →  ~90k pasos para figura-8, max|ΔE/E₀| ~ 0.03 (preciso)
 *   η = 1e-5  →  ~900k pasos,              max|ΔE/E₀| ~ 0.003 (muy preciso)
 */
class ARChain3Integrator {
public:
    /**
     * @param eta   Parámetro de control de paso. El paso en tiempo ficticio
     *              se calcula como ds = eta * sep_min / Omega.
     *              Valor típico: 1e-3 (rápido) a 1e-4 (preciso).
     */
    explicit ARChain3Integrator(double eta = 1e-3);

    // ── Interfaz principal ───────────────────────────────────────────────────

    /**
     * @brief Construye el estado AR-chain desde el sistema físico.
     *
     * Inicializa X1, X2, V1, V2, cm_pos, cm_vel, masas y energía.
     * t_phys y s_fict se ponen a 0.
     */
    ARChain3State initialize(const NBodySystem& system,
                             int idx1, int idx2, int idx3) const;

    /**
     * @brief Modo bloque: integra dt_phys de tiempo físico desde state.t_phys.
     *
     * Uso: integrador jerárquico con subsistemas parciales (N > 3).
     * El integrador exterior impone el intervalo; el AR-chain integra
     * adaptativamente dentro de ese intervalo.
     *
     * @param state      Estado de entrada/salida.
     * @param dt_phys    Intervalo de tiempo físico a integrar (delta).
     */
    void integrate(ARChain3State& state, double dt_phys) const;

    /**
     * @brief Modo directo: integra desde state.t_phys hasta t_abs_final.
     *
     * Uso: sistema completo donde el triple es el único subsistema activo.
     * El AR-chain controla su propio tiempo sin cortes externos impuestos.
     *
     * DIFERENCIA CRÍTICA respecto a integrate():
     *   - integrate()    recibe un DELTA de tiempo (dt_phys)
     *   - integrate_to() recibe un TIEMPO ABSOLUTO (t_abs_final)
     *
     * El estado puede tener t_phys > 0 si fue inicializado previamente.
     * La integración avanza desde state.t_phys hasta t_abs_final.
     *
     * @param state         Estado de entrada/salida.
     * @param t_abs_final   Tiempo físico absoluto objetivo.
     */
    void integrate_to(ARChain3State& state, double t_abs_final) const;

    /**
     * @brief Escribe las posiciones y velocidades de vuelta al sistema físico.
     *
     * Las posiciones incluyen el desplazamiento del CM: r_i + cm_vel * t_phys.
     * Las velocidades son absolutas (incluyen cm_vel).
     */
    void write_back(const ARChain3State& state, NBodySystem& system,
                    int idx1, int idx2, int idx3) const;

    // ── Métodos públicos para tests ──────────────────────────────────────────

    /** @brief Función de tiempo Ω = Σᵢ<ⱼ mᵢmⱼ/rᵢⱼ. */
    double compute_Omega(const ARChain3State& state) const;

    /** @brief Energía total T + U calculada desde el estado actual. */
    double compute_energy(const ARChain3State& state) const;

    /** @brief Separación mínima entre los tres cuerpos. */
    double compute_sep_min(const ARChain3State& state) const;

protected:
    double eta_;      ///< Parámetro de control de paso
    double ds_min_;   ///< Paso mínimo (evita división por cero)
    double ds_max_;   ///< Paso máximo (evita pasos enormes cuando sep>>1)

    void compute_accelerations(const ARChain3State& state,
                               Vec3& A1, Vec3& A2) const;

    /**
     * @brief Núcleo DKD: un paso leapfrog en tiempo ficticio. Acepta ds < 0.
     *
     * Versión "bare" usada internamente por yoshida4_step (que necesita
     * sub-pasos con coeficiente c0 < 0). No hace comprobación de signo.
     */
    void leapfrog_step_bare(ARChain3State& state, double ds) const;

    /**
     * @brief Wrapper público de leapfrog_step_bare (compatibilidad con tests).
     */
    void leapfrog_step(ARChain3State& state, double ds) const;

    /**
     * @brief Un paso de orden 4 mediante composición de Yoshida (1990).
     *
     * Aplica tres sub-pasos leapfrog con coeficientes c1, c0, c1:
     *   φ_ds^(4) = φ_{c1·ds} ∘ φ_{c0·ds} ∘ φ_{c1·ds}
     *
     * Coste: 3× leapfrog orden 2. Ganancia: error O(η^4) vs O(η^2).
     * Referencia: Yoshida (1990); Mikkola & Merritt (2006).
     */
    void yoshida4_step(ARChain3State& state, double ds) const;

    /**
     * @brief Nucleo del bucle de integracion (compartido por integrate e integrate_to).
     *
     * Avanza el estado desde state.t_phys hasta t_target usando pasos
     * adaptativos ds = η·sep_min/Ω. Recorta el ultimo paso para no
     * sobrepasar t_target. Usa yoshida4_step internamente (orden 4).
     *
     * @param state     Estado de entrada/salida.
     * @param t_target  Tiempo fisico absoluto hasta el que integrar.
     */
    void integrate_until(ARChain3State& state, double t_target) const;
};