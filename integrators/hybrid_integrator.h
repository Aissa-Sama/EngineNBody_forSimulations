// integrators/hybrid_integrator.h
#pragma once

#include <memory>
#include <vector>
#include <limits>

#include "integrator.h"
#include "nbody_system.h"
#include "vec3.h"
#include "binary_state.h"
#include "ks_perturbed_integrator.h"
#include "chain3_integrator.h"
#include "regime_logger.h"

// ============================================================================
// STRUCTS DE DATOS
// ============================================================================

/// Candidato a triple detectado por detect_triple()
struct TripleCandidate {
    int    i, j, k;          ///< Índices de los tres cuerpos (orden de cadena)
    double binding_energy;   ///< Energía de ligadura estimada del triple
};

/// Par ligado detectado por detect_all_bound_pairs()
struct BinaryPair {
    int    i, j;             ///< Índices de los dos cuerpos
    double binding_energy;   ///< Energía de ligadura (negativa si ligado)
};

// ============================================================================
// HYBRID INTEGRATOR — Ruta A (legacy)
//
// Jerarquía de decisión en cada paso:
//   1. ¿Hay un triple fuertemente acoplado?  → Chain3Integrator
//   2. ¿Hay pares ligados?                   → KS perturbado
//   3. El resto                              → integrador de campo (far)
//
// Esta clase es el integrador "Ruta A" del proyecto. Para simulaciones de
// precisión se prefiere HierarchicalIntegrator (Ruta B). HybridIntegrator
// se mantiene como referencia histórica y para test_hybrid_chain3.
// ============================================================================
class HybridIntegrator : public Integrator {
public:
    /// Constructor
    /// @param far_integrator  Integrador de campo lejano (e.g. VelocityVerlet)
    /// @param r_close         Radio de detección de encuentros cercanos
    /// @param ks_dt           Paso interno del integrador KS perturbado
    /// @param logger          Logger de regímenes (puede ser nullptr)
    HybridIntegrator(
        std::unique_ptr<Integrator> far_integrator,
        double r_close_,
        double ks_dt,
        RegimeLogger* logger_ = nullptr
    );

    /// Paso principal: detecta subsistemas y delega al integrador apropiado
    void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) override;

private:
    // -----------------------------------------------------------------------
    // DETECCIÓN
    // -----------------------------------------------------------------------

    /// ¿Tres cuerpos i,j,k están mutuamente dentro de r_close?
    bool is_close_triple(
        const NBodySystem& system,
        int i, int j, int k
    ) const;

    /// ¿Un tercer cuerpo invade la zona KS de la binaria (bi, bj)?
    /// Si es así, devuelve true y escribe el índice en third_idx.
    bool third_body_too_close(
        const NBodySystem& system,
        int bi, int bj,
        int& third_idx
    ) const;

    /// Detecta el triple más fuertemente acoplado del sistema.
    /// Devuelve true si existe; llena `triple` con los índices y energía.
    bool detect_triple(
        const NBodySystem& system,
        TripleCandidate& triple
    ) const;

    /// ¿El par (i,j) es una binaria ligada (E < 0) dentro de r_close?
    bool is_bound_binary(
        const NBodySystem& system,
        int i, int j,
        double& binding_energy
    ) const;

    /// Detecta todos los pares ligados del sistema.
    std::vector<BinaryPair> detect_all_bound_pairs(
        const NBodySystem& system
    ) const;

    /// Selecciona los pares óptimos (sin solapamiento) ordenados por energía.
    std::vector<BinaryPair> select_optimal_binaries(
        std::vector<BinaryPair>& candidates,
        std::vector<bool>& used,
        const NBodySystem& system
    ) const;

    // -----------------------------------------------------------------------
    // INTEGRACIÓN POR SUBSISTEMA
    // -----------------------------------------------------------------------

    /// Integra un triple con Chain3Integrator
    void handle_triple(
        NBodySystem& system,
        const TripleCandidate& triple,
        double dt
    );

    /// Integra una binaria con KSPerturbedIntegrator
    void handle_binary_ks(
        NBodySystem& system,
        const BinaryPair& bin,
        double dt
    );

    /// Calcula fuerzas externas sobre binarias (perturbaciones de campo)
    void compute_external_forces(
        const NBodySystem& system,
        const std::vector<bool>& in_binary,
        std::vector<Vec3>& acc_binary,
        std::vector<Vec3>& acc_field
    ) const;

    // -----------------------------------------------------------------------
    // MIEMBROS
    // -----------------------------------------------------------------------
    std::unique_ptr<Integrator> far;          ///< Integrador de campo lejano
    KSPerturbedIntegrator       ks_perturbed; ///< KS con perturbaciones
    Chain3Integrator            chain3;       ///< Integrador de triple
    double                      r_close;      ///< Radio de detección
    RegimeLogger*               logger;       ///< Logger (puede ser nullptr)
    size_t                      step_counter = 0; ///< Contador de pasos
};