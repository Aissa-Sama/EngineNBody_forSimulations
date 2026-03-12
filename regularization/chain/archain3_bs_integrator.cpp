// regularization/chain/archain3_bs_integrator.cpp
#include "archain3_bs_integrator.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// ============================================================================
// CONSTRUCTOR
// ============================================================================

ARChain3BSIntegrator::ARChain3BSIntegrator(double eta,
                                           const BSParameters& params)
    : ARChain3Integrator(eta)
    , params_(params)
{}

// ============================================================================
// MODIFIED MIDPOINT STEP (Gragg 1965)
//
// Avanza y0 en ds_total usando n sub-pasos uniformes de tamaño h = ds_total/n.
//
// NOTA SOBRE LAS ECUACIONES DE MOVIMIENTO EN TIEMPO FICTICIO:
//   Las ecuaciones de movimiento del TTL en tiempo ficticio τ son:
//     dXᵢ/dτ = Ω(X)·Vᵢ
//     dVᵢ/dτ = Ω(X)·Aᵢ(X)
//     dt_phys/dτ = Ω(X)
//
//   Omega Ω = Σ mᵢmⱼ/rᵢⱼ DEPENDE de las posiciones X, por lo que actúa
//   como una función de escalado no lineal. El paso medio modificado evalúa
//   Ω en la posición z_curr.X del sub-paso actual — esto es correcto y
//   mantiene la simetría temporal del método (a diferencia de usar Ω_fijo).
//
// NOTA SOBRE EL PASO EULER INICIAL vs LEAPFROG:
//   El MMP usa un Euler explícito (orden 1) como integrador base, NO el
//   leapfrog. Esto es correcto: la potencia del GBS viene de la extrapolación,
//   no del integrador base. El leapfrog se usa solo en integrate_to() del
//   padre para integración directa sin GBS.
//
// NOTA SOBRE t_phys y s_fict:
//   t_phys se propaga acumulando Ω·ds en cada sub-paso (= la derivada
//   dt_phys/dτ = Ω integrada con paso h). La acumulación en el paso medio
//   es consistente con la del resultado final: result.t_phys es la suma
//   de todos los incrementos h·Ω evaluados en cada sub-paso.
// ============================================================================

void ARChain3BSIntegrator::modified_midpoint_step(
    const ARChain3State& y0,
    double               ds_total,
    int                  n,
    ARChain3State&       result) const
{
    const double h = ds_total / static_cast<double>(n);

    ARChain3State z_prev = y0;
    ARChain3State z_curr = y0;

    // ── Primer paso: Euler ─────────────────────────────────────────────────
    //   z₁ = z₀ + h · f(z₀)
    const double Omega0 = compute_Omega(z_prev);
    Vec3 A1_0, A2_0;
    compute_accelerations(z_prev, A1_0, A2_0);

    const double dt0 = h * Omega0;
    z_curr.X1     = z_prev.X1 + z_prev.V1 * dt0;
    z_curr.X2     = z_prev.X2 + z_prev.V2 * dt0;
    z_curr.V1     = z_prev.V1 + A1_0 * dt0;
    z_curr.V2     = z_prev.V2 + A2_0 * dt0;
    z_curr.t_phys = z_prev.t_phys + dt0;
    z_curr.s_fict = z_prev.s_fict + h;
    // cm_pos, cm_vel, masas, energy: heredadas de y0 (invariantes)

    // ── Pasos medios: m = 1 .. n-1 ────────────────────────────────────────
    //   z_{m+1} = z_{m-1} + 2h · f(z_m)
    for (int m = 1; m < n; ++m) {
        const double Omega_m = compute_Omega(z_curr);
        Vec3 A1_m, A2_m;
        compute_accelerations(z_curr, A1_m, A2_m);

        const double dt_m = 2.0 * h * Omega_m;

        ARChain3State z_next = z_curr;     // hereda cm_pos, cm_vel, masas
        z_next.X1     = z_prev.X1 + z_curr.V1 * dt_m;
        z_next.X2     = z_prev.X2 + z_curr.V2 * dt_m;
        z_next.V1     = z_prev.V1 + A1_m * dt_m;
        z_next.V2     = z_prev.V2 + A2_m * dt_m;
        z_next.t_phys = z_prev.t_phys + dt_m;
        z_next.s_fict = z_curr.s_fict + h;

        z_prev = z_curr;
        z_curr = z_next;
    }

    // ── Paso final suavizante ──────────────────────────────────────────────
    //   result = ½·(z_n + z_{n-1} + h·f(z_n))
    const double Omega_n = compute_Omega(z_curr);
    Vec3 A1_n, A2_n;
    compute_accelerations(z_curr, A1_n, A2_n);

    const double dt_n = h * Omega_n;

    result        = y0;     // hereda cm_pos, cm_vel, masas, energy
    result.X1     = 0.5 * (z_curr.X1 + z_prev.X1 + z_curr.V1 * dt_n);
    result.X2     = 0.5 * (z_curr.X2 + z_prev.X2 + z_curr.V2 * dt_n);
    result.V1     = 0.5 * (z_curr.V1 + z_prev.V1 + A1_n * dt_n);
    result.V2     = 0.5 * (z_curr.V2 + z_prev.V2 + A2_n * dt_n);
    result.t_phys = 0.5 * (z_curr.t_phys + z_prev.t_phys + dt_n);
    result.s_fict = y0.s_fict + ds_total;
}

// ============================================================================
// EXTRAPOLACIÓN POLINÓMICA DE RICHARDSON (Neville-Aitken)
//
// Construye la tabla de Richardson T[k][m] usando h² como variable de
// extrapolación (no h), aprovechando que la expansión de error del MMP
// sobre un integrador time-symmetric contiene únicamente potencias pares.
//
// Fórmula de recurrencia:
//   T[k][m] = T[k][m-1] + (T[k][m-1] - T[k-1][m-1]) / ((h_{k-m}/h_k)² - 1)
//
// El denominador (ratio² - 1), con ratio = step_sizes[k-m] / step_sizes[k],
// cancela el término de orden 2m en la expansión de error. Cada barrido
// de m elimina un orden adicional.
//
// IMPLEMENTACIÓN IN-PLACE:
//   Se trabaja con un vector de trabajo T[] que se actualiza de derecha
//   a izquierda para evitar usar los valores ya actualizados en el mismo
//   barrido. Al final, T[k_current] contiene la extrapolación de orden
//   2*(k_current+1).
// ============================================================================

void ARChain3BSIntegrator::polynomial_extrapolation(
    const std::vector<ARChain3State>& raw,
    const std::vector<double>&        step_sizes,
    int                               k_current,
    ARChain3State&                    extrap_out) const
{
    // Copia de trabajo: T[j] = estimación del nivel j
    std::vector<ARChain3State> T(raw.begin(),
                                 raw.begin() + k_current + 1);

    for (int m = 1; m <= k_current; ++m) {
        for (int k = k_current; k >= m; --k) {
            const double ratio  = step_sizes[k - m] / step_sizes[k];
            double       factor = ratio * ratio - 1.0;
            if (std::abs(factor) < 1e-15) factor = 1e-15;

            // Extrapolación campo a campo del estado
            T[k].X1 = T[k].X1 + (T[k].X1 - T[k-1].X1) * (1.0 / factor);
            T[k].X2 = T[k].X2 + (T[k].X2 - T[k-1].X2) * (1.0 / factor);
            T[k].V1 = T[k].V1 + (T[k].V1 - T[k-1].V1) * (1.0 / factor);
            T[k].V2 = T[k].V2 + (T[k].V2 - T[k-1].V2) * (1.0 / factor);
            T[k].t_phys = T[k].t_phys
                        + (T[k].t_phys - T[k-1].t_phys) * (1.0 / factor);
        }
    }

    extrap_out = T[k_current];
}

// ============================================================================
// ESTIMACIÓN DE ERROR
//
// Error relativo entre dos extrapolaciones consecutivas (niveles k y k-1).
// Usa la norma euclídea de las diferencias escaladas por el máximo de cada
// componente para medir el error relativo de forma homogénea.
//
// Se promedian las cuatro componentes (X1, V1, X2, V2) para obtener una
// medida escalar del error que refleja simultáneamente el error en posición
// y en velocidad.
// ============================================================================

double ARChain3BSIntegrator::estimate_error(
    const ARChain3State& hi,
    const ARChain3State& lo) const
{
    const double sX1 = std::max(hi.X1.norm(), 1e-10);
    const double sV1 = std::max(hi.V1.norm(), 1e-10);
    const double sX2 = std::max(hi.X2.norm(), 1e-10);
    const double sV2 = std::max(hi.V2.norm(), 1e-10);

    return ((hi.X1 - lo.X1).norm() / sX1
          + (hi.V1 - lo.V1).norm() / sV1
          + (hi.X2 - lo.X2).norm() / sX2
          + (hi.V2 - lo.V2).norm() / sV2) * 0.25;
}

// ============================================================================
// FACTOR DE ESCALA DEL PASO
//
// Basado en la fórmula estándar GBS (Numerical Recipes §17.3):
//   scale = safety · (bs_eps / error)^(1/(2·k+1))
//
// El exponente 1/(2k+1) refleja que el método de orden 2k tiene error ∝ h^(2k).
// Un factor de seguridad de 0.9 reduce el paso sugerido ligeramente por debajo
// del óptimo teórico para aumentar la tasa de aceptación en el siguiente intento.
//
// Límites: [0.1, 4.0] — evita pasos irracionales tras errores extremos.
// ============================================================================

double ARChain3BSIntegrator::step_scale_factor(double error, int k_opt) const
{
    if (error < 1e-30) return 4.0;
    const double exponent = 1.0 / (2.0 * k_opt + 1.0);
    const double scale    = params_.safety
                          * std::pow(params_.bs_eps / error, exponent);
    return std::clamp(scale, 0.1, 4.0);
}

// ============================================================================
// ENERGY_ERROR — diagnóstico
// ============================================================================

double ARChain3BSIntegrator::energy_error(const ARChain3State& state,
                                          double E_ref) const
{
    return std::abs(compute_energy(state) - E_ref)
         / (std::abs(E_ref) + 1e-30);
}

// ============================================================================
// UN PASO GBS COMPLETO
//
// Intenta avanzar state en ds_total aplicando el MMP con n₁,n₂,...,nₖ
// sub-pasos y extrapolando al límite h→0 mediante la tabla de Richardson.
//
// FLUJO:
//   Para k = 0, 1, ..., K_MAX-1:
//     1. Aplica MMP con n_seq_[k] sub-pasos → raw[k]
//     2. Actualiza la tabla de Richardson → extrap_curr
//     3. Si k >= 2: estima el error entre extrap_curr y extrap_prev
//     4. Si error < bs_eps: ACEPTA el paso
//     5. Si error > 1000*bs_eps y k >= 3: ABANDONA (paso demasiado grande)
//   Si no converge: acepta la mejor estimación con ds_next reducido.
//
// CONTROL DE PASO SUGERIDO (ds_next_out):
//   Si el paso se acepta: ds_next = ds_total * step_scale_factor(error, k)
//   Si el paso falla:     ds_next = ds_total * 0.5 (reducción conservadora)
//
// NOTA: el estado solo se modifica si el paso es ACEPTADO.
// ============================================================================

bool ARChain3BSIntegrator::bs_step(
    ARChain3State& state,
    double         ds_total,
    double&        error_out,
    double&        ds_next_out)
{
    const ARChain3State y_save = state;

    std::vector<ARChain3State> raw;
    std::vector<double>        step_sizes;
    raw.reserve(K_MAX);
    step_sizes.reserve(K_MAX);

    ARChain3State extrap_prev;
    bool  has_prev  = false;
    error_out       = 1e30;

    for (int k = 0; k < params_.k_max && k < K_MAX; ++k) {
        const int    n      = n_seq_[k];
        const double ds_k   = ds_total / static_cast<double>(n);

        // ── Paso medio modificado ──────────────────────────────────────────
        ARChain3State r;
        try {
            modified_midpoint_step(y_save, ds_total, n, r);
        } catch (const std::exception&) {
            // Singularidad numérica: el paso es demasiado grande
            ds_next_out = ds_total * 0.25;
            return false;
        }
        raw.push_back(r);
        step_sizes.push_back(ds_k);

        if (k >= 1) {
            // ── Extrapolación de Richardson ────────────────────────────────
            ARChain3State extrap_curr;
            polynomial_extrapolation(raw, step_sizes, k, extrap_curr);

            if (has_prev) {
                // ── Estimación y criterio de error ────────────────────────
                error_out = estimate_error(extrap_curr, extrap_prev);
                const double scaled = error_out / params_.bs_eps;

                if (scaled <= 1.0) {
                    // ✅ Convergió: aceptar
                    state       = extrap_curr;
                    ds_next_out = ds_total * step_scale_factor(error_out, k);
                    return true;
                }

                // Abandono temprano si el error es muy grande
                if (scaled > 1000.0 && k >= 3) break;
            }
            extrap_prev = extrap_curr;
            has_prev    = true;
        }
    }

    // ── No convergió en K_MAX etapas: acepta la mejor estimación ──────────
    // Esto ocurre raramente (ds demasiado grande o encuentro muy cercano).
    // El paso se acepta igualmente para no bloquear la integración, pero
    // ds_next_out se reduce para el siguiente intento.
    if (has_prev) {
        state       = extrap_prev;
        ds_next_out = ds_total * 0.5;
        return true;
    }

    ds_next_out = ds_total * 0.5;
    return false;
}

// ============================================================================
// INTEGRATE_TO_BS — bucle de integración principal
//
// Avanza state desde state.t_phys hasta t_abs_final usando pasos GBS
// adaptativos. El control de paso opera en tiempo ficticio τ (ds), pero
// se monitoriza t_phys para la condición de parada.
//
// ESTIMACIÓN DEL ds INICIAL:
//   Igual que el leapfrog TTL base: ds = η · sep_min / Ω
//   Esto garantiza que el primer paso GBS tenga un tamaño razonablemente
//   bueno para el encuentro más cercano en ese instante.
//
// RECORTE DEL ÚLTIMO PASO:
//   A diferencia del leapfrog TTL puro, aquí SÍ recortamos el ds para
//   terminar exactamente en t_abs_final. El MMP funciona correctamente
//   con pasos pequeños y arbitrarios. La estimación del ds_final se hace
//   de forma conservadora usando t_left / (Ω_actual * safety_factor).
//
// RECHAZO POR ENERGÍA:
//   Si |ΔE/E| > params_.energy_tol después de un paso aceptado, el paso
//   se rechaza y ds se reduce. Esto previene la acumulación de error en
//   la extrapolación cerca de singularidades.
// ============================================================================

void ARChain3BSIntegrator::integrate_to_bs(ARChain3State& state,
                                           double         t_abs_final)
{
    if (t_abs_final <= state.t_phys) return;

    const double E0       = state.energy;
    double       ds       = params_.initial_ds;
    double       max_dE   = 0.0;
    int          n_steps  = 0;
    int          n_reject = 0;
    int          n_consecutive_reject = 0;

    while (state.t_phys < t_abs_final - 1e-14) {

        // ── Tiempo restante y estimación de ds ────────────────────────────
        const double t_left = t_abs_final - state.t_phys;

        // ds adaptativo basado en separación (igual que el leapfrog base)
        double Omega   = compute_Omega(state);
        double sep_min = compute_sep_min(state);
        double ds_eta  = eta_ * sep_min / Omega;
        ds_eta = std::clamp(ds_eta, params_.min_ds, params_.max_ds);

        // No iniciar un paso que claramente sobrepasaría t_final por mucho
        // Estimamos: dt_fisico_esperado ≈ ds * Omega
        // Para no sobrepasar: ds ≤ t_left / Omega * factor
        const double ds_max_t = t_left / (Omega * 0.5 + 1e-30);
        ds = std::min({ds, ds_eta, ds_max_t});
        ds = std::clamp(ds, params_.min_ds, params_.max_ds);

        // ── Intento de paso GBS ───────────────────────────────────────────
        ARChain3State trial = state;
        double error = 0.0, ds_next = ds;

        bool accepted = false;
        try {
            accepted = bs_step(trial, ds, error, ds_next);
            n_consecutive_reject = 0;
        } catch (const std::runtime_error&) {
            ds *= 0.25;
            ds = std::max(ds, params_.min_ds);
            ++n_reject;
            ++n_consecutive_reject;
            if (n_consecutive_reject > 30)
                throw std::runtime_error(
                    "ARChain3BS: singularidad persistente");
            continue;
        }

        if (!accepted) {
            ds = std::clamp(ds_next, params_.min_ds, params_.max_ds);
            ++n_reject;
            ++n_consecutive_reject;
            if (n_consecutive_reject > 30)
                throw std::runtime_error(
                    "ARChain3BS: no convergencia persistente");
            continue;
        }
        n_consecutive_reject = 0;

        // ── Verificar que no sobrepasamos t_final ─────────────────────────
        if (trial.t_phys > t_abs_final + 1e-12) {
            // El paso físico superó el objetivo: reducir ds y reintentar
            ds *= 0.5;
            ds = std::max(ds, params_.min_ds);
            ++n_reject;
            continue;
        }

        // ── Control de energía ────────────────────────────────────────────
        const double dE = energy_error(trial, E0);
        if (dE > params_.energy_tol
            && ds > params_.min_ds * 4.0
            && n_reject < params_.max_steps / 2) {
            ds = std::max(ds * 0.5, params_.min_ds);
            ++n_reject;
            continue;
        }

        // ── Aceptar ───────────────────────────────────────────────────────
        state = trial;
        ds    = std::clamp(ds_next, params_.min_ds, params_.max_ds);
        ++n_steps;
        if (dE > max_dE) max_dE = dE;

        if (n_steps > params_.max_steps)
            throw std::runtime_error(
                "ARChain3BS: maximo de pasos excedido");
    }

    if (params_.verbose) {
        std::cout << "[AR-BS] " << n_steps  << " pasos, "
                  << n_reject << " rechazados"
                  << "  max|dE/E| = " << max_dE
                  << "  t_final = " << state.t_phys
                  << "\n";
    }
}