// python/nbody_bindings.cpp
// ============================================================================
// PYBIND11 BINDINGS — nbody_core
//
// Expone a Python las clases y funciones principales del motor N-body:
//   Vec3, Body, NBodySystem, InitialConditions,
//   HierarchicalIntegrator, LeapfrogIntegrator,
//   phys::Observables, phys::compute_observables
// ============================================================================
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>

#include "vec3.h"
#include "body.h"
#include "nbody_system.h"
#include "initial_conditions.h"
#include "leapfrog_integrator.h"
#include "hierarchical_integrator.h"
#include "hierarchy_builder.h"
#include "physics_observables.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(nbody_core, m) {
    m.doc() = "nbody_core — Motor N-body gravitacional de alta precisión (C++17)";

    // ========================================================================
    // Vec3
    // ========================================================================
    py::class_<Vec3>(m, "Vec3")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z)
        .def("norm",  &Vec3::norm,  "Norma euclidiana |v|")
        .def("norm2", &Vec3::norm2, "Norma al cuadrado |v|²")
        .def("dot",   &Vec3::dot,   "Producto punto", py::arg("other"))
        // Operadores aritméticos
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())
        .def(-py::self)
        .def("__repr__", [](const Vec3& v) {
            return "Vec3(" + std::to_string(v.x) + ", "
                           + std::to_string(v.y) + ", "
                           + std::to_string(v.z) + ")";
        });

    // Funciones libres de Vec3
    m.def("dot",   &dot,   "Producto punto de dos Vec3", py::arg("a"), py::arg("b"));
    m.def("cross", &cross, "Producto vectorial de dos Vec3", py::arg("a"), py::arg("b"));

    // ========================================================================
    // Body
    // ========================================================================
    py::class_<Body>(m, "Body")
        .def(py::init<>())
        .def_readwrite("position", &Body::position)
        .def_readwrite("velocity", &Body::velocity)
        .def_readwrite("mass",     &Body::mass)
        .def("__repr__", [](const Body& b) {
            return "Body(m=" + std::to_string(b.mass)
                + ", r=(" + std::to_string(b.position.x) + ", "
                + std::to_string(b.position.y) + ", "
                + std::to_string(b.position.z) + "))";
        });

    // ========================================================================
    // NBodySystem
    // ========================================================================
    py::class_<NBodySystem>(m, "NBodySystem")
        .def(py::init<>())
        .def_readwrite("bodies", &NBodySystem::bodies)
        .def_readwrite("G",      &NBodySystem::G)
        .def("compute_accelerations",    &NBodySystem::compute_accelerations)
        .def("invalidate_accelerations", &NBodySystem::invalidate_accelerations)
        .def("kinetic_energy",           &NBodySystem::kinetic_energy)
        .def("potential_energy",         &NBodySystem::potential_energy)
        .def("total_energy",             &NBodySystem::total_energy)
        .def("total_momentum",           &NBodySystem::total_momentum)
        .def("total_angular_momentum",   &NBodySystem::total_angular_momentum)
        .def("__repr__", [](const NBodySystem& sys) {
            return "NBodySystem(N=" + std::to_string(sys.bodies.size())
                + ", G=" + std::to_string(sys.G) + ")";
        });

    // ========================================================================
    // InitialConditions (métodos estáticos)
    // ========================================================================
    py::class_<InitialConditions>(m, "InitialConditions")
        .def_static("solar_system",      &InitialConditions::solar_system,
            "Sistema solar simplificado")
        .def_static("binary_with_field", &InitialConditions::binary_with_field,
            "Binaria + estrella de campo")
        .def_static("plummer_cluster",   &InitialConditions::plummer_cluster,
            "Cúmulo Plummer", py::arg("n_bodies"))
        .def_static("figure_eight",      &InitialConditions::figure_eight,
            "Órbita figura-8 (3 cuerpos, estable)")
        .def_static("random_system",     &InitialConditions::random_system,
            "Sistema aleatorio", py::arg("n_bodies"), py::arg("radius"))
        .def_static("kepler_binary",     &InitialConditions::kepler_binary,
            "Binaria Kepleriana", py::arg("a") = 1.0, py::arg("e") = 0.0);

    // ========================================================================
    // HierarchyBuilder::Params  (para configurar el integrador)
    // ========================================================================
    py::class_<HierarchyBuilder::Params::PNParams>(m, "PNParams")
        .def(py::init<>())
        .def_readwrite("enabled",            &HierarchyBuilder::Params::PNParams::enabled)
        .def_readwrite("c_speed",            &HierarchyBuilder::Params::PNParams::c_speed)
        .def_readwrite("pn_order",           &HierarchyBuilder::Params::PNParams::pn_order)
        .def_readwrite("activation_epsilon", &HierarchyBuilder::Params::PNParams::activation_epsilon);

    py::class_<HierarchyBuilder::Params>(m, "HierarchyParams")
        .def(py::init<>())
        .def_readwrite("r_ks_threshold",      &HierarchyBuilder::Params::r_ks_threshold)
        .def_readwrite("tidal_threshold",     &HierarchyBuilder::Params::tidal_threshold)
        .def_readwrite("strong_coupling_eta", &HierarchyBuilder::Params::strong_coupling_eta)
        .def_readwrite("ar_chain_threshold",  &HierarchyBuilder::Params::ar_chain_threshold)
        .def_readwrite("ar_chain_eta",        &HierarchyBuilder::Params::ar_chain_eta)
        .def_readwrite("pn",                  &HierarchyBuilder::Params::pn);

    // ========================================================================
    // HierarchicalIntegrator
    // ========================================================================
    py::class_<HierarchicalIntegrator>(m, "HierarchicalIntegrator")
        .def(py::init([](double r_ks, double ks_dt,
                         const HierarchyBuilder::Params& params) {
            return std::make_unique<HierarchicalIntegrator>(
                std::make_unique<LeapfrogIntegrator>(),
                r_ks, ks_dt, params, nullptr
            );
        }),
            py::arg("r_ks_threshold") = 1.0,
            py::arg("ks_internal_dt") = 1e-4,
            py::arg("params") = HierarchyBuilder::Params{},
            "Crea un HierarchicalIntegrator con LeapfrogIntegrator como far-field"
        )
        .def("step", [](HierarchicalIntegrator& self,
                         NBodySystem& system, double dt) {
            std::vector<bool> used(system.bodies.size(), false);
            self.step(system, dt, used);
        }, py::arg("system"), py::arg("dt"),
           "Realiza un paso de integración")
        .def("step_to",
            static_cast<void (HierarchicalIntegrator::*)(
                NBodySystem&, double, double, std::function<void(double)>
            )>(&HierarchicalIntegrator::step_to),
            py::arg("system"),
            py::arg("t_final"),
            py::arg("dt_hint"),
            py::arg("on_step") = nullptr,
            "Integra hasta t_final con pasos de tamaño dt_hint")
        .def("pn_active_count", &HierarchicalIntegrator::pn_active_count,
            "Número de grupos con PN activo en el último paso");

    // ========================================================================
    // Physics Observables (namespace phys)
    // ========================================================================
    py::class_<phys::Observables>(m, "Observables")
        .def(py::init<>())
        .def_readonly("t",           &phys::Observables::t)
        .def_readonly("E",           &phys::Observables::E)
        .def_readonly("E_kinetic",   &phys::Observables::E_kinetic)
        .def_readonly("E_potential", &phys::Observables::E_potential)
        .def_readonly("dE_rel",      &phys::Observables::dE_rel)
        .def_readonly("P",           &phys::Observables::P)
        .def_readonly("dP_rel",      &phys::Observables::dP_rel)
        .def_readonly("L",           &phys::Observables::L)
        .def_readonly("dL_rel",      &phys::Observables::dL_rel)
        .def_readonly("cm_pos",      &phys::Observables::cm_pos)
        .def_readonly("cm_vel",      &phys::Observables::cm_vel)
        .def_readonly("sep_min",     &phys::Observables::sep_min)
        .def_readonly("sep_max",     &phys::Observables::sep_max)
        .def_readonly("virial",      &phys::Observables::virial)
        .def_readonly("regime",      &phys::Observables::regime);

    // Funciones libres del namespace phys
    m.def("kinetic_energy",     &phys::kinetic_energy,     py::arg("system"));
    m.def("potential_energy",   &phys::potential_energy,    py::arg("system"), py::arg("softening") = 0.0);
    m.def("total_energy",       &phys::total_energy,       py::arg("system"), py::arg("softening") = 0.0);
    m.def("total_momentum",     &phys::total_momentum,     py::arg("system"));
    m.def("total_angular_momentum", &phys::total_angular_momentum, py::arg("system"));
    m.def("total_mass",         &phys::total_mass,         py::arg("system"));
    m.def("min_separation",     &phys::min_separation,     py::arg("system"));
    m.def("max_separation",     &phys::max_separation,     py::arg("system"));
    m.def("virial_ratio",       &phys::virial_ratio,       py::arg("system"));
    m.def("compute_observables", &phys::compute_observables,
        py::arg("system"),
        py::arg("t"),
        py::arg("E0"),
        py::arg("P0"),
        py::arg("L0"),
        py::arg("regime") = "",
        py::arg("softening") = 0.0,
        "Calcula todos los observables físicos en un instante t");
}
