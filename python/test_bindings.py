"""Test rápido de nbody_core Python bindings."""
import nbody_core as nc

# Vec3
v = nc.Vec3(1.0, 2.0, 3.0)
print(f"Vec3: {v}")
print(f"  norm = {v.norm():.6f}")
print(f"  norm2 = {v.norm2():.6f}")
print(f"  dot(v,v) = {v.dot(v):.6f}")

# Body
b = nc.Body()
b.mass = 1.0
b.position = nc.Vec3(0, 0, 0)
b.velocity = nc.Vec3(1, 0, 0)
print(f"\nBody: {b}")

# InitialConditions - figure eight
sys = nc.InitialConditions.figure_eight()
print(f"\nFigure-8: {sys}")
print(f"  N = {len(sys.bodies)}")
print(f"  E0 = {nc.total_energy(sys):.15e}")
print(f"  |P| = {nc.total_momentum(sys).norm():.6e}")
print(f"  sep_min = {nc.min_separation(sys):.6f}")

# Kepler binary
kb = nc.InitialConditions.kepler_binary(1.0, 0.0)
print(f"\nKepler binary: N={len(kb.bodies)}, E={nc.total_energy(kb):.10e}")

# HierarchicalIntegrator
integrator = nc.HierarchicalIntegrator(r_ks_threshold=1.0, ks_internal_dt=1e-4)
E0 = nc.total_energy(sys)
integrator.step(sys, 0.01)
E1 = nc.total_energy(sys)
print(f"\nDespués de 1 paso (dt=0.01):")
print(f"  E = {E1:.15e}")
print(f"  dE/E0 = {abs(E1-E0)/abs(E0):.6e}")

print("\n=== TODOS LOS TESTS PASARON ===")
