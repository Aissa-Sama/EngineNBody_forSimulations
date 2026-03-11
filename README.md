# nbody_core

Motor N-body gravitacional de alta precisión en C++17, construido sobre métodos del estado del arte en dinámica estelar computacional. El objetivo es reproducir fielmente los algoritmos usados en los códigos de referencia de la comunidad (ARCHAIN, NBODY6/7) con arquitectura limpia, tests verificados y documentación rigurosa.

---

## Resultados clave

| Método | Sistema | `dE/E` | Referencia |
|---|---|---|---|
| Leapfrog TTL ord. 2, η=1e-3 | Figura-8 | 19.3% | Mikkola & Tanikawa (1999) |
| Leapfrog TTL ord. 2, η=1e-4 | Figura-8 | 3.32% | ídem |
| **GBS sobre TTL, bs_eps=1e-10** | **Figura-8** | **1.78×10⁻¹³** | **Mikkola & Aarseth (2002)** |
| GBS sobre TTL, ×10 períodos | Figura-8 | 2.5×10⁻¹¹ | estable, sin deriva secular |

La ganancia de GBS respecto al leapfrog puro es de **~10¹³×** en precisión de energía, a un coste de ~10× en evaluaciones de fuerza por paso físico.

---

## Compilar

Requiere: MSVC 19.x (o GCC/Clang), CMake ≥ 3.16, Ninja.

```powershell
cmake -S . -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl
cmake --build build --config Release
```

---

## Tests

```powershell
cd build
.\test_archain3.exe       # Leapfrog TTL — 6/6
.\test_archain3_bs.exe    # GBS sobre TTL — 6/6
.\test_hierarchical.exe   # Integrador jerárquico — 4/4
```

Todos los tests tienen criterios calibrados con referencias publicadas (Simó 2002, figura-8 de 15 dígitos significativos). No hay umbrales arbitrarios.

---

## Ejecutar una simulación

```powershell
cd build

# Figura-8, un período, modo directo (sin cortes de dt externos)
.\nbody_sim.exe --scenario figure8 --t-final 6.3259 --eta 1e-4 --direct-mode --output run1

# Tres cuerpos random, modo jerárquico
.\nbody_sim.exe --scenario random --t-final 50.0 --dt 0.01 --eta 1e-3 --output run2
```

Salida: archivos `<output>_positions.csv` y `<output>_energy.csv`.

Opciones disponibles:

| Opción | Descripción | Default |
|---|---|---|
| `--scenario` | `binary`, `field`, `figure8`, `solar`, `random`, `plummer` | `binary` |
| `--t-final` | Tiempo físico total | `10.0` |
| `--dt` | Paso de campo lejano | `0.01` |
| `--eta` | Parámetro η del AR-chain TTL | `1e-3` |
| `--r-ks` | Umbral de regularización KS (u.a.) | `1.0` |
| `--ar-thresh` | Umbral de activación AR-chain (u.a.) | `0.1` |
| `--direct-mode` | Sin bloques `dt` externos; AR-chain controla su propio tiempo | off |
| `--output` | Prefijo de archivos CSV | `sim` |
| `--softening` | Suavizado gravitacional | `0.0` |

---

## Arquitectura

```
nbody_core/
├── core/                          # Vec3, Vec4, Body, NBodySystem
├── integrators/
│   ├── hierarchical_integrator    # Integrador principal: detecta régimen y delega
│   ├── leapfrog_integrator        # Campo lejano
│   └── rk45_integrator            # Alternativa adaptativa
├── regularization/
│   ├── hierarchy/
│   │   └── hierarchy_builder      # Construye el árbol jerárquico en cada paso
│   ├── ks/
│   │   ├── ks_integrator          # Regularización Kustaanheimo-Stiefel (par)
│   │   └── ks_perturbed_integrator
│   └── chain/
│       ├── archain3_integrator    # AR-chain TTL — leapfrog DKD orden 2
│       ├── archain3_bs_integrator # AR-chain TTL — extrapolación GBS orden ~16
│       ├── chain3_integrator      # Cadena KS clásica
│       └── chain3_bs_integrator   # Cadena KS con GBS
├── analysis/                      # Detección de binarias, observables físicos
├── io/                            # Logger de regímenes, snapshots CSV
├── validation/                    # Herramientas de validación (Viviani, etc.)
└── tests/                         # Suite de tests verificados
```

### Flujo de integración (modo directo, N=3)

```
main.cpp
  └── HierarchicalIntegrator::step_to(system, t_final)
        └── HierarchyBuilder::build() → TRIPLE_AR_CHAIN
              └── ARChain3BSIntegrator::integrate_to_bs(state, t_final)
                    ├── bs_step() × N                     ← GBS adaptativo
                    │     ├── modified_midpoint_step()    ← Gragg 1965
                    │     └── polynomial_extrapolation()  ← Neville-Aitken
                    └── write_back() → NBodySystem
```

---

## Métodos implementados

### AR-chain TTL (Mikkola & Tanikawa 1999)

Regularización algorítmica en coordenadas físicas de cadena. La función de tiempo TTL es:

```
dt/ds = Ω(r)    donde    Ω = Σᵢ<ⱼ mᵢmⱼ / rᵢⱼ
```

Cuando `rᵢⱼ → 0`, `Ω → ∞`, haciendo que el paso físico `dt → 0` automáticamente en encuentros cercanos. El leapfrog DKD resultante es simpléctico: conserva una energía sombra acotada para todo `t`, sin deriva secular.

El paso adaptativo es `ds = η · sep_min / Ω`, que produce `dt = η · sep_min` — proporcional a la separación mínima actual.

### Extrapolación GBS sobre TTL (Mikkola & Aarseth 2002)

El leapfrog TTL es time-symmetric (DKD). Por el teorema de Gragg (1965), su expansión de error contiene solo potencias pares de `h`. La extrapolación de Richardson por lo tanto sube el orden de dos en dos por etapa. Con la secuencia de subdivisiones `n = {2, 4, 6, 8, 10, 12, 14, 16}` y K=8 etapas se alcanza orden 16 teórico — en la práctica, el límite es la aritmética de doble precisión (~14 dígitos).

### Regularización KS (Kustaanheimo & Stiefel 1965)

Para pares binarios, las variables de posición se transforman a cuaterniones 4D mediante la transformación KS. La singularidad `1/r` desaparece en las ecuaciones de movimiento transformadas. Se usa para encuentros binarios dentro del integrador jerárquico.

### Integrador jerárquico

En cada paso, `HierarchyBuilder` analiza el sistema y construye un árbol jerárquico:

- **PAIR_KS**: par binario cercano → regularización KS
- **TRIPLE_AR_CHAIN**: triple cercano → AR-chain TTL (+ GBS)
- **LEAF**: cuerpo aislado en campo lejano → leapfrog

El integrador delega a cada subsistema según su régimen, permitiendo mezclar métodos dentro de una misma simulación.

---

## Referencias

- Mikkola & Tanikawa (1999), *MNRAS* 310, 745 — AR-chain TTL
- Mikkola & Aarseth (2002), *Celest. Mech. Dyn. Astron.* 84, 343 — GBS sobre TTL
- Mikkola & Merritt (2006), *AJ* 130, 2345 — punto medio generalizado (fuerzas PN, Fase 5)
- Gragg (1965), *SIAM J. Numer. Anal.* 2, 384 — expansión error en `h²`
- Bulirsch & Stoer (1966), *Numer. Math.* 8, 1 — extrapolación GBS
- Kustaanheimo & Stiefel (1965), *J. Reine Angew. Math.* 218, 204 — regularización KS
- Chenciner & Montgomery (2000), *Ann. Math.* 152, 881 — figura-8
- Simó (2002) — condiciones iniciales de la figura-8, 15 dígitos significativos
- Szebehely & Peters (1967) — problema de Pitágoras (masas 3, 4, 5)

---

## Hoja de ruta

- [x] **Fase 1** — Rediseño arquitectural: `step_to()` sin bloques `dt` externos
- [x] **Fase 2** — Extrapolación GBS sobre TTL: `ARChain3BSIntegrator`
- [ ] **Fase 3** — Validación con benchmarks publicados (perihelios, eyección Pitágoras completa)
- [ ] **Fase 4** — Extensión a N > 3 en el AR-chain
- [ ] **Fase 5** — Física adicional: mareas, correcciones post-Newtonianas (Mikkola & Merritt 2006)
