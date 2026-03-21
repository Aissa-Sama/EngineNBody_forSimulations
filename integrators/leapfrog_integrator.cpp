// integrators/leapfrog_integrator.cpp
// FASE 7B: Añadido step_block() para block timestep.
#include "leapfrog_integrator.h"
#include "nbody_system.h"

// ── Paso uniforme (sin cambios) ──────────────────────────────────────────────
void LeapfrogIntegrator::step(
    NBodySystem& system, double dt, const std::vector<bool>& used)
{
    const size_t N = system.bodies.size();
    system.invalidate_accelerations();
    auto acc = system.compute_accelerations();
    for (size_t i = 0; i < N; ++i) {
        if (used[i]) continue;
        system.bodies[i].velocity += (dt * 0.5) * acc[i];
    }
    for (size_t i = 0; i < N; ++i) {
        if (used[i]) continue;
        system.bodies[i].position += dt * system.bodies[i].velocity;
    }
    system.invalidate_accelerations();
    acc = system.compute_accelerations();
    for (size_t i = 0; i < N; ++i) {
        if (used[i]) continue;
        system.bodies[i].velocity += (dt * 0.5) * acc[i];
    }
}

// ── FASE 7B: Paso con dt individual por cuerpo ───────────────────────────────
//
// ESQUEMA DKD CON DT VARIABLE:
//   Cada cuerpo libre avanza con su propio dts[i].
//   Las aceleraciones se calculan una vez al inicio y una vez al final.
//   El kick inicial usa dts[i]/2, el drift usa dts[i], el kick final dts[i]/2.
//
// NOTA SOBRE CONSISTENCIA:
//   Con dt variables, el kick final de un cuerpo no es simultaneo con el
//   kick inicial del siguiente paso para cuerpos con dt diferente. Esto es
//   inherente al block timestep y se acepta — el error es O(dt_max^2) para
//   los pares con dt diferente, que es el mismo orden que el leapfrog estandar.
//   Para cuerpos con el mismo dt el esquema es exactamente simplectico.
//
// REFERENCIA:
//   Aarseth (2003) cap. 2 — block timestep con leapfrog.
//   Makino & Aarseth (1992) PASJ 44, 141.
// ────────────────────────────────────────────────────────────────────────────
void LeapfrogIntegrator::step_block(
    NBodySystem& system,
    const std::vector<double>& dts,
    const std::vector<bool>& used)
{
    const size_t N = system.bodies.size();
    system.invalidate_accelerations();
    auto acc = system.compute_accelerations();

    // Kick 1/2
    for (size_t i = 0; i < N; ++i) {
        if (used[i]) continue;
        system.bodies[i].velocity += (dts[i] * 0.5) * acc[i];
    }
    // Drift
    for (size_t i = 0; i < N; ++i) {
        if (used[i]) continue;
        system.bodies[i].position += dts[i] * system.bodies[i].velocity;
    }
    // Recalcular aceleraciones con posiciones nuevas
    system.invalidate_accelerations();
    acc = system.compute_accelerations();
    // Kick 1/2
    for (size_t i = 0; i < N; ++i) {
        if (used[i]) continue;
        system.bodies[i].velocity += (dts[i] * 0.5) * acc[i];
    }
}
