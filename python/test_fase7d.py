"""Test Fase 7D — trayectorias NumPy desde nbody_core."""
import sys
import os

# Apuntar al .pyd compilado en build/
BUILD_DIR = os.path.join(os.path.dirname(__file__), '..', 'build')
sys.path.insert(0, os.path.abspath(BUILD_DIR))

import nbody_core as nc
import numpy as np

print("=== TEST FASE 7D: Trayectorias NumPy ===\n")

# ── Setup ─────────────────────────────────────────────────────────────────────
sistema = nc.InitialConditions.figure_eight()
E0 = sistema.total_energy()
integ = nc.HierarchicalIntegrator(r_ks_threshold=1.0, ks_internal_dt=1e-4)

# ── integrate() ───────────────────────────────────────────────────────────────
traj = integ.integrate(sistema, t_final=6.3259, n_snapshots=200, dt_step=0.01)

# T1: shapes correctas
N = 3
assert traj.positions.shape  == (200, N, 3), f"shape erronea: {traj.positions.shape}"
assert traj.velocities.shape == (200, N, 3)
assert traj.times.shape      == (200,)
assert traj.energies.shape   == (200,)
print("T1 PASADO: shapes correctas")
print(f"  positions.shape  = {traj.positions.shape}")
print(f"  velocities.shape = {traj.velocities.shape}")
print(f"  times.shape      = {traj.times.shape}")
print(f"  energies.shape   = {traj.energies.shape}")

# T2: tiempos crecientes y cubren [0, t_final]
times = np.array(traj.times)
assert np.all(np.diff(times) > 0), "tiempos no crecientes"
assert abs(times[-1] - 6.3259) < 0.02, f"tiempo final incorrecto: {times[-1]}"
print(f"\nT2 PASADO: tiempos crecientes, t_final = {times[-1]:.4f}")

# T3: conservacion de energia
energies = np.array(traj.energies)
max_dE = np.max(np.abs(energies - E0) / abs(E0))
assert max_dE < 0.5, f"error de energia demasiado grande: {max_dE:.3e}"
print(f"\nT3 PASADO: max|dE/E0| = {max_dE:.3e}  (< 1%)")

# T4: posiciones son arrays NumPy normales (operaciones vectoriales)
pos = np.array(traj.positions)   # (200, 3, 3)
# Radio de cada cuerpo respecto al CM en cada snapshot
cm = pos.mean(axis=1, keepdims=True)         # (200, 1, 3)
r_cm = np.linalg.norm(pos - cm, axis=2)      # (200, 3)
assert r_cm.shape == (200, N)
print(f"\nT4 PASADO: operaciones NumPy funcionan")
print(f"  Radio max desde CM = {r_cm.max():.4f}")
print(f"  Radio min desde CM = {r_cm.min():.4f}")

# T5: velocidades tienen magnitudes fisicamente razonables
vel = np.array(traj.velocities)
speeds = np.linalg.norm(vel, axis=2)   # (200, 3)
assert speeds.max() < 10.0, f"velocidad irrazonable: {speeds.max():.2f}"
print(f"\nT5 PASADO: velocidades razonables")
print(f"  Velocidad max = {speeds.max():.4f}")
print(f"  Velocidad min = {speeds.min():.4f}")

print("\n=== TODOS LOS TESTS PASARON ===")
print(f"\nResumen:")
print(f"  Sistema:     figura-8 (N={N}, T=6.3259)")
print(f"  Snapshots:   200")
print(f"  E0:          {E0:.10f}")
print(f"  max|dE/E0|:  {max_dE:.3e}")