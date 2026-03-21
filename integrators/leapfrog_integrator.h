// integrators/leapfrog_integrator.h
// FASE 7B: Añadida sobrecarga step_block() para block timestep.
#pragma once
#include <vector>
#include "integrator.h"
#include "nbody_system.h"

class LeapfrogIntegrator : public Integrator {
public:
    /// Paso uniforme — todos los cuerpos libres avanzan dt.
    void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) override;

    /// FASE 7B: Paso con dt individual por cuerpo.
    /// dts[i] = paso de tiempo del cuerpo i.
    /// used[i] = true → saltar (ya integrado por subsistema regularizado).
    /// Cada cuerpo avanza con su propio dts[i].
    void step_block(
        NBodySystem& system,
        const std::vector<double>& dts,
        const std::vector<bool>& used
    );
};
