[INFO] Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 17 entities
Info    : 7904 nodes
Info    : 15806 elements
Info    : Done reading 'mesh.msh'
[INFO] Mesh successfully loaded from Gmsh file.
[INFO] Mesh topology dimension d=2
[INFO] 
Available volume tags (dx):
[INFO]   Tag ID: 1
[INFO] 
Unique tags found in facet data: [2 3 4 5 6]
[INFO] Label map loaded from geometry:
[INFO]   uo2          → 1
[INFO]   ymin         → 2
[INFO]   xmin         → 3
[INFO]   xmax         → 4
[INFO]   ymax         → 5
[INFO]   cavity       → 6
[INFO]   Lz = 0.000 m
[INFO]   Lx = 0.000 m, Ly = 0.000 m
[INFO]   area = 2.700e-09 m², perimeter = 2.100e-04 m
[INFO] === Mesh summary ===
[INFO]   Topology dim: 2
[INFO]   Facet dim: 1
[INFO]   Num cells: 15094
[INFO]   Cell tags: {np.int32(1)}
[INFO]   Facet tags: {np.int32(2), np.int32(3), np.int32(4), np.int32(5), np.int32(6)}
[INFO]   Geometry type: rect


***

Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
Author: Giovanni Zullo
Version: 0.1.0 (2025)

***



## Description

Z3ST is an open-source framework for the thermo-mechanical modelling
of materials. Built on FEniCSx, it supports transient simulations,
complex geometries, and user-defined boundary conditions.


### Config initializer

  → Geometry            : geometry.yaml
  → Mesh                : mesh.msh
  → Boundary conditions : boundary_conditions.yaml
  → Time steps          : 201
  → Regime              : 2d
  → Models active       :
      thermal    → OFF
      mechanical → ON
      damage     → ON
      cluster    → OFF
      plasticity → OFF
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, triangle, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, triangle, 1, gll_warped, unset, False, float64, []), (2,)))
Scalar function space (V_d): FunctionSpace(<Mesh #0>, Basix element (P, triangle, 1, gll_warped, unset, False, float64, []))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, triangle, 0, gll_warped, unset, True, float64, []))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 0.9
  → Displacement : 0.7
  → Damage       : 0.4
  Adaptive relaxation enabled
  → relax_growth  : 1.1
  → relax_shrink : 0.7
  → relax_min  : 0.05
  → relax_max : 1.0


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  order               : 1
  linear_solver       : direct_mumps
  rtol                : 1e-06
  stag_tol            : 0.0001
  convergence         : rel_norm
  debug               : False
DamageModel initializer
Options loaded from input.yaml:
  type                : AT1
  split               : amor
  solver              : linear
  linear_solver       : direct_mumps
  rtol                : 1e-06
  stag_tol            : 0.0001
  convergence         : rel_norm
  lc                  : 2e-06
  hybrid_constraint   : True
[spine.load_materials]
Material loaded: uo2
  → k defined as constant: 5.0
  → Gc not defined for uo2
  - Material 'uo2': Gc (AT1) from sigma_c = 4.22e+08 Pa
  → constitutive model: lame
  E               → 358000000000.0 (float)
  G               → 145528455284.55286 (float)
  Gc              → 2.6467374301675974 (float)
  T_initial       → 1023.15 (float)
  T_ref           → 1023.15 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 220987654320.98764 (float)
  constitutive_mode → lame (str)
  cp              → 280.0 (float)
  k               → 5.0 (float)
  lmbda           → 123968684131.28575 (float)
  name            → UO2 (str)
  nu              → 0.23 (float)
  rho             → 10970.0 (float)
  sigma_c         → 421500000.0 (float)
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

Initializing the damage field...



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Clamp_y mechanical BC on 'uo2' → 0.0 (first step) at region 'ymin'
  **[INFO]** Clamp_x mechanical BC on 'uo2' → 0.0 (first step) at region 'xmin'
  **[INFO]** Clamp_x mechanical BC on 'uo2' → 0.0 (first step) at region 'xmax'
  **[INFO]** Neumann mechanical BC on 'uo2' → cavity: 0.0 Pa (list loaded)

Setting damage boundary conditions...
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/201: t = 0.00e+00 s | LHR = 0.00e+00 W/m

  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → 0.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.70

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.40
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 0.0000e+00 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 0.0000e+00 J


## Step 02/201: t = 5.00e-03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -4000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.70

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.40
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -4000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.000e-01
  [adaptive] relax_u=0.77

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.28
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -4000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.900e-02
  [adaptive] relax_u=0.85

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.20
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -4000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.505e-02
  [adaptive] relax_u=0.93

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.14
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -4000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.215e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.10
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -4000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.090e-04
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -4000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.984e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9061e-09 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.9061e-09 J


## Step 03/201: t = 1.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -8000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.000e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -8000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.350e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1624e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.1624e-08 J


## Step 04/201: t = 1.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -12000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.333e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -12000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.857e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1155e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.1155e-08 J


## Step 05/201: t = 2.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -16000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.500e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -16000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.045e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2650e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2650e-07 J


## Step 06/201: t = 2.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -20000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.000e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -20000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.239e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9765e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.9765e-07 J


## Step 07/201: t = 3.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -24000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.667e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -24000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.381e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8462e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8462e-07 J


## Step 08/201: t = 3.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -28000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.429e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -28000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.079e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8740e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.8740e-07 J


## Step 09/201: t = 4.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -32000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.250e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -32000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.807e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0599e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.0599e-07 J


## Step 10/201: t = 4.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -36000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.111e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -36000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.978e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4040e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.4040e-07 J


## Step 11/201: t = 5.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -40000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.000e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -40000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.375e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9061e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.9061e-07 J


## Step 12/201: t = 5.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -44000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.091e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -44000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.625e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.5664e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.5664e-07 J


## Step 13/201: t = 6.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -48000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.333e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -48000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.548e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1385e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1385e-06 J


## Step 14/201: t = 6.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -52000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.692e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -52000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.042e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3361e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3361e-06 J


## Step 15/201: t = 7.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -56000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.143e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -56000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.953e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5496e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5496e-06 J


## Step 16/201: t = 7.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -60000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.667e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -60000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.044e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7789e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.7789e-06 J


## Step 17/201: t = 8.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -64000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.250e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -64000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.898e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0240e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0240e-06 J


## Step 18/201: t = 8.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -68000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.882e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -68000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.273e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2849e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2849e-06 J


## Step 19/201: t = 9.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -72000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.556e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -72000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.932e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5616e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5616e-06 J


## Step 20/201: t = 9.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -76000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.263e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -76000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.533e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8541e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8541e-06 J


## Step 21/201: t = 1.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -80000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.000e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -80000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.871e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1624e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.1624e-06 J


## Step 22/201: t = 1.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -84000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.762e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -84000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.447e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4866e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.4866e-06 J


## Step 23/201: t = 1.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -88000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.545e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -88000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.614e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8266e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.8266e-06 J


## Step 24/201: t = 1.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -92000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.348e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -92000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.236e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1823e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.1823e-06 J


## Step 25/201: t = 1.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -96000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.167e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -96000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.408e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5539e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.5539e-06 J


## Step 26/201: t = 1.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -100000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.000e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -100000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.376e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9413e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.9413e-06 J


## Step 27/201: t = 1.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -104000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.846e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -104000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.311e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3445e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.3445e-06 J


## Step 28/201: t = 1.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -108000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.704e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -108000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.736e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7636e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.7636e-06 J


## Step 29/201: t = 1.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -112000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.571e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -112000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.467e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1984e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.1984e-06 J


## Step 30/201: t = 1.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -116000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.448e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -116000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.416e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6490e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.6490e-06 J


## Step 31/201: t = 1.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -120000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.333e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -120000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.601e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1155e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.1155e-06 J


## Step 32/201: t = 1.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -124000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.226e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -124000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.5978e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.5978e-06 J


## Step 33/201: t = 1.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -128000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.125e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -128000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.634e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0959e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.0959e-06 J


## Step 34/201: t = 1.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -132000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.030e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -132000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.705e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.6098e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.6098e-06 J


## Step 35/201: t = 1.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -136000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.941e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -136000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.273e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.1395e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.1395e-06 J


## Step 36/201: t = 1.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -140000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.857e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -140000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.018e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6850e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.6850e-06 J


## Step 37/201: t = 1.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -144000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.778e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -144000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.940e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0246e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0246e-05 J


## Step 38/201: t = 1.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -148000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.703e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -148000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.447e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0823e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0823e-05 J


## Step 39/201: t = 1.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -152000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.632e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -152000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.567e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1416e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1416e-05 J


## Step 40/201: t = 1.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -156000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.564e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -156000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.157e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2025e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2025e-05 J


## Step 41/201: t = 2.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -160000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.500e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -160000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.235e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2650e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2650e-05 J


## Step 42/201: t = 2.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 41 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -164000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.439e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -164000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.741e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3290e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3290e-05 J


## Step 43/201: t = 2.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 42 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -168000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.381e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -168000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.420e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3946e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3946e-05 J


## Step 44/201: t = 2.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 43 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -172000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.326e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -172000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.419e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4618e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4618e-05 J


## Step 45/201: t = 2.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 44 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -176000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.273e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -176000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.757e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5306e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5306e-05 J


## Step 46/201: t = 2.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 45 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -180000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.222e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -180000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6010e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.6010e-05 J


## Step 47/201: t = 2.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 46 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -184000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.174e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -184000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.785e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6729e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.6729e-05 J


## Step 48/201: t = 2.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 47 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -188000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.128e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -188000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.356e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7465e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.7465e-05 J


## Step 49/201: t = 2.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 48 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -192000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.083e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -192000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.953e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8216e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.8216e-05 J


## Step 50/201: t = 2.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 49 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -196000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.041e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -196000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.435e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8983e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.8983e-05 J


## Step 51/201: t = 2.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 50 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -200000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.000e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -200000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.449e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9765e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.9765e-05 J


## Step 52/201: t = 2.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 51 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -204000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.961e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -204000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.760e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0564e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0564e-05 J


## Step 53/201: t = 2.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 52 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -208000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.923e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -208000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.848e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1378e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.1378e-05 J


## Step 54/201: t = 2.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 53 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -212000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.887e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -212000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.420e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2208e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2208e-05 J


## Step 55/201: t = 2.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 54 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -216000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.852e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -216000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.327e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3054e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.3054e-05 J


## Step 56/201: t = 2.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 55 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -220000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.818e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -220000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.385e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3916e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.3916e-05 J


## Step 57/201: t = 2.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 56 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -224000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.786e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -224000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.269e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4794e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.4794e-05 J


## Step 58/201: t = 2.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 57 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -228000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.754e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -228000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.230e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5687e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5687e-05 J


## Step 59/201: t = 2.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 58 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -232000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.724e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -232000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.438e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6596e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.6596e-05 J


## Step 60/201: t = 2.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 59 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -236000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.695e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -236000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.405e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7521e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.7521e-05 J


## Step 61/201: t = 3.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 60 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -240000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.667e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -240000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.603e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8462e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8462e-05 J


## Step 62/201: t = 3.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 61 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -244000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.639e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -244000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.128e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9419e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.9419e-05 J


## Step 63/201: t = 3.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 62 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -248000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.613e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -248000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0391e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.0391e-05 J


## Step 64/201: t = 3.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 63 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -252000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.587e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -252000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1379e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.1379e-05 J


## Step 65/201: t = 3.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 64 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -256000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.562e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -256000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.950e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2383e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.2383e-05 J


## Step 66/201: t = 3.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 65 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -260000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.538e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -260000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.244e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3403e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.3403e-05 J


## Step 67/201: t = 3.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 66 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -264000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.515e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -264000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.551e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4439e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.4439e-05 J


## Step 68/201: t = 3.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 67 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -268000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.493e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -268000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.234e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5491e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.5491e-05 J


## Step 69/201: t = 3.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 68 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -272000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.471e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -272000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.060e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6558e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.6558e-05 J


## Step 70/201: t = 3.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 69 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -276000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.449e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -276000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.635e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7641e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.7641e-05 J


## Step 71/201: t = 3.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 70 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -280000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.429e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -280000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.490e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8740e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.8740e-05 J


## Step 72/201: t = 3.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 71 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -284000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.408e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -284000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.990e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9855e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.9855e-05 J


## Step 73/201: t = 3.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 72 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -288000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.389e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -288000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.731e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0985e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.0985e-05 J


## Step 74/201: t = 3.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 73 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -292000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.370e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -292000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.865e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2132e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.2132e-05 J


## Step 75/201: t = 3.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 74 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -296000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.351e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -296000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.560e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3294e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.3294e-05 J


## Step 76/201: t = 3.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 75 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -300000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.333e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -300000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.807e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4472e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.4472e-05 J


## Step 77/201: t = 3.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 76 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -304000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.316e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -304000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5666e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.5666e-05 J


## Step 78/201: t = 3.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 77 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -308000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.299e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -308000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.622e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6875e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.6875e-05 J


## Step 79/201: t = 3.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 78 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -312000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.282e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -312000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.491e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8101e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.8101e-05 J


## Step 80/201: t = 3.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 79 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -316000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.266e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -316000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.218e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9342e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.9342e-05 J


## Step 81/201: t = 4.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 80 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -320000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.250e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -320000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.489e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0599e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.0599e-05 J


## Step 82/201: t = 4.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 81 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -324000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.235e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -324000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1872e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.1872e-05 J


## Step 83/201: t = 4.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 82 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -328000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.220e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -328000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3161e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.3161e-05 J


## Step 84/201: t = 4.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 83 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -332000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.205e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -332000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.188e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4465e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.4465e-05 J


## Step 85/201: t = 4.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 84 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -336000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.190e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -336000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.138e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5786e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.5786e-05 J


## Step 86/201: t = 4.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 85 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -340000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.176e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -340000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.261e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7122e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.7122e-05 J


## Step 87/201: t = 4.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 86 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -344000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.163e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -344000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.648e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.8474e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.8474e-05 J


## Step 88/201: t = 4.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 87 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -348000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.149e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -348000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.897e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9841e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.9841e-05 J


## Step 89/201: t = 4.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 88 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -352000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.136e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -352000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.582e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1225e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.1225e-05 J


## Step 90/201: t = 4.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 89 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -356000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.124e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -356000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.011e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2624e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.2624e-05 J


## Step 91/201: t = 4.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 90 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -360000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.111e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -360000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.006e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4040e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.4040e-05 J


## Step 92/201: t = 4.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 91 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -364000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.099e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -364000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.794e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.5471e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.5471e-05 J


## Step 93/201: t = 4.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 92 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -368000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.087e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -368000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.290e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6917e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.6917e-05 J


## Step 94/201: t = 4.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 93 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -372000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.075e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -372000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.877e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8380e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.8380e-05 J


## Step 95/201: t = 4.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 94 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -376000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.064e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -376000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.347e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9858e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.9858e-05 J


## Step 96/201: t = 4.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 95 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -380000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.053e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -380000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.305e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1353e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.1353e-05 J


## Step 97/201: t = 4.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 96 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -384000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.042e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -384000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.423e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.2863e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.2863e-05 J


## Step 98/201: t = 4.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 97 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -388000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.031e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -388000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.029e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4389e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.4389e-05 J


## Step 99/201: t = 4.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 98 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -392000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.020e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -392000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.896e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.5930e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.5930e-05 J


## Step 100/201: t = 4.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 99 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -396000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.010e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -396000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.106e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.7488e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.7488e-05 J


## Step 101/201: t = 5.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 100 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.000e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 1.610e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.349e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.500e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.530e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.365e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.928e-01
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 1.599e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.714e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.529e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.374e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.120e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.912e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.435e-03

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.727e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.264e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 1.492e-03

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.574e-01
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 1.542e-03

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.707e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.715e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.081e-03

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.733e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.017e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.130e-03

Convergence check


#### Iteration 10/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.728e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.294e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 1.175e-03

Convergence check


#### Iteration 11/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.649e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.538e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 1.214e-03

Convergence check


#### Iteration 12/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.793e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.740e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 1.246e-03

Convergence check


#### Iteration 13/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.701e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.891e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.271e-03

Convergence check


#### Iteration 14/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.628e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.079e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 8.179e-04

Convergence check


#### Iteration 15/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.352e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.272e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 8.490e-04

Convergence check


#### Iteration 16/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.980e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.439e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 8.760e-04

Convergence check


#### Iteration 17/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.380e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.575e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 8.979e-04

Convergence check


#### Iteration 18/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.163e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.673e-01
  [adaptive] relax_D=0.09
  |ΔD|_∞ = 9.135e-04

Convergence check


#### Iteration 19/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.721e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.725e-01
  [adaptive] relax_D=0.10
  |ΔD|_∞ = 9.220e-04

Convergence check


#### Iteration 20/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.337e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.726e-01
  [adaptive] relax_D=0.11
  |ΔD|_∞ = 9.221e-04

Convergence check


#### Iteration 21/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.451e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.669e-01
  [adaptive] relax_D=0.12
  |ΔD|_∞ = 9.130e-04

Convergence check


#### Iteration 22/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.452e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.551e-01
  [adaptive] relax_D=0.13
  |ΔD|_∞ = 8.940e-04

Convergence check


#### Iteration 23/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.351e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.368e-01
  [adaptive] relax_D=0.15
  |ΔD|_∞ = 8.646e-04

Convergence check


#### Iteration 24/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.626e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.120e-01
  [adaptive] relax_D=0.16
  |ΔD|_∞ = 8.246e-04

Convergence check


#### Iteration 25/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.625e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.809e-01
  [adaptive] relax_D=0.18
  |ΔD|_∞ = 7.745e-04

Convergence check


#### Iteration 26/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.528e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.439e-01
  [adaptive] relax_D=0.19
  |ΔD|_∞ = 7.149e-04

Convergence check


#### Iteration 27/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.539e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.019e-01
  [adaptive] relax_D=0.21
  |ΔD|_∞ = 6.473e-04

Convergence check


#### Iteration 28/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.717e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.561e-01
  [adaptive] relax_D=0.24
  |ΔD|_∞ = 5.735e-04

Convergence check


#### Iteration 29/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.745e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.078e-01
  [adaptive] relax_D=0.26
  |ΔD|_∞ = 4.958e-04

Convergence check


#### Iteration 30/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.445e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.589e-01
  [adaptive] relax_D=0.28
  |ΔD|_∞ = 4.169e-04

Convergence check


#### Iteration 31/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.961e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.110e-01
  [adaptive] relax_D=0.31
  |ΔD|_∞ = 3.398e-04

Convergence check


#### Iteration 32/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.526e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.660e-01
  [adaptive] relax_D=0.34
  |ΔD|_∞ = 2.673e-04

Convergence check


#### Iteration 33/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.724e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.254e-01
  [adaptive] relax_D=0.38
  |ΔD|_∞ = 2.019e-04

Convergence check


#### Iteration 34/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.735e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.036e-02
  [adaptive] relax_D=0.42
  |ΔD|_∞ = 1.455e-04

Convergence check


#### Iteration 35/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.787e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.170e-02
  [adaptive] relax_D=0.46
  |ΔD|_∞ = 9.937e-05

Convergence check


#### Iteration 36/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.352e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.956e-02
  [adaptive] relax_D=0.50
  |ΔD|_∞ = 6.371e-05

Convergence check


#### Iteration 37/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.993e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.355e-02
  [adaptive] relax_D=0.56
  |ΔD|_∞ = 3.793e-05

Convergence check


#### Iteration 38/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.992e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.283e-02
  [adaptive] relax_D=0.61
  |ΔD|_∞ = 2.066e-05

Convergence check


#### Iteration 39/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.350e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.277e-03
  [adaptive] relax_D=0.67
  |ΔD|_∞ = 1.011e-05

Convergence check


#### Iteration 40/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.721e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.688e-03
  [adaptive] relax_D=0.74
  |ΔD|_∞ = 4.329e-06

Convergence check


#### Iteration 41/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.231e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.704e-04
  [adaptive] relax_D=0.81
  |ΔD|_∞ = 1.563e-06

Convergence check


#### Iteration 42/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.539e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.786e-04
  [adaptive] relax_D=0.89
  |ΔD|_∞ = 4.487e-07

Convergence check


#### Iteration 43/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -400000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.729e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.735e-05
  [adaptive] relax_D=0.98
  |ΔD|_∞ = 9.236e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 43 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9041e-05 J
  → Fracture energy : 6.0932e-09 J
  → Total energy    : 7.9047e-05 J


## Step 102/201: t = 5.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 101 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.013e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.144e-01
  [adaptive] relax_D=0.98
  |ΔD|_∞ = 4.902e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.714e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.337e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.047e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.753e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.231e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.343e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.994e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.138e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.577e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0584e-05 J
  → Fracture energy : 4.0180e-08 J
  → Total energy    : 8.0624e-05 J


## Step 103/201: t = 5.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 102 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -408000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.089e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.710e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.304e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -408000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.195e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.996e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.388e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.2174e-05 J
  → Fracture energy : 1.2125e-07 J
  → Total energy    : 8.2295e-05 J


## Step 104/201: t = 5.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 103 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -412000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.223e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.791e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.933e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -412000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.589e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.540e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.3831e-05 J
  → Fracture energy : 2.8374e-07 J
  → Total energy    : 8.4115e-05 J


## Step 105/201: t = 5.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 104 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -416000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.456e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.229e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.013e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -416000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.116e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.727e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.327e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5583e-05 J
  → Fracture energy : 5.8544e-07 J
  → Total energy    : 8.6168e-05 J


## Step 106/201: t = 5.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 105 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -420000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.868e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.850e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.308e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -420000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.162e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.729e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7472e-05 J
  → Fracture energy : 1.1156e-06 J
  → Total energy    : 8.8588e-05 J


## Step 107/201: t = 5.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 106 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -424000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.615e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.547e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.739e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -424000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.770e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.682e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.9614e-05 J
  → Fracture energy : 2.0629e-06 J
  → Total energy    : 9.1677e-05 J


## Step 108/201: t = 5.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 107 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -428000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.984e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.310e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.387e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -428000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.695e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.782e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2324e-05 J
  → Fracture energy : 3.6351e-06 J
  → Total energy    : 9.5959e-05 J


## Step 109/201: t = 5.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 108 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -432000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.491e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.091e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.257e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -432000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.300e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.409e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.759e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6583e-05 J
  → Fracture energy : 5.9864e-06 J
  → Total energy    : 1.0257e-04 J


## Step 110/201: t = 5.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 109 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -436000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.522e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.922e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.187e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -436000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.831e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.767e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0378e-04 J
  → Fracture energy : 9.3849e-06 J
  → Total energy    : 1.1316e-04 J


## Step 111/201: t = 5.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 110 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -440000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.995e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.783e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.827e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -440000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.798e-20
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.421e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9685e-03 J
  → Fracture energy : 1.4727e-05 J
  → Total energy    : 6.9833e-03 J


## Step 112/201: t = 5.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 111 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -444000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.548e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.755e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.289e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -444000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.1952e-02 J
  → Fracture energy : 2.3709e-05 J
  → Total energy    : 9.1976e-02 J


## Step 113/201: t = 5.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 112 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -448000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.055e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.059e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.885e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -448000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.385e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3230e-01 J
  → Fracture energy : 4.1944e-05 J
  → Total energy    : 3.3235e-01 J


## Step 114/201: t = 5.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 113 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -452000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.639e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.279e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.361e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -452000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.265e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.126e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1405e-01 J
  → Fracture energy : 8.2762e-05 J
  → Total energy    : 8.1413e-01 J


## Step 115/201: t = 5.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 114 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -456000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.618e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.641e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.632e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -456000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.323e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6518e+00 J
  → Fracture energy : 1.2356e-04 J
  → Total energy    : 1.6519e+00 J


## Step 116/201: t = 5.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 115 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -460000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.455e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.494e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.638e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -460000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.896e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.419e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7711e+01 J
  → Fracture energy : 1.4416e-04 J
  → Total energy    : 1.7712e+01 J


## Step 117/201: t = 5.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 116 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -464000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.041e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.330e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.384e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -464000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.064e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0352e+01 J
  → Fracture energy : 1.5609e-04 J
  → Total energy    : 3.0352e+01 J


## Step 118/201: t = 5.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 117 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -468000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.848e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.211e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.165e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -468000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.054e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.061e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7986e+01 J
  → Fracture energy : 1.6901e-04 J
  → Total energy    : 3.7987e+01 J


## Step 119/201: t = 5.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 118 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -472000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.132e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.167e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.712e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -472000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.328e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.472e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2338e+01 J
  → Fracture energy : 1.7985e-04 J
  → Total energy    : 4.2338e+01 J


## Step 120/201: t = 5.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 119 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -476000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.824e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.142e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.249e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -476000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.276e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.848e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6802e+01 J
  → Fracture energy : 1.8680e-04 J
  → Total energy    : 4.6802e+01 J


## Step 121/201: t = 6.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 120 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -480000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.236e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.218e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.523e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -480000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.076e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.393e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0910e+01 J
  → Fracture energy : 1.9520e-04 J
  → Total energy    : 5.0910e+01 J


## Step 122/201: t = 6.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 121 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -484000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.441e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.318e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.798e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -484000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.818e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.976e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4887e+01 J
  → Fracture energy : 2.0352e-04 J
  → Total energy    : 5.4888e+01 J


## Step 123/201: t = 6.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 122 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -488000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.533e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.098e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.433e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -488000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.718e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.991e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9169e+01 J
  → Fracture energy : 2.1195e-04 J
  → Total energy    : 5.9169e+01 J


## Step 124/201: t = 6.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 123 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -492000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.761e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.261e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.973e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -492000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.951e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.939e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2795e+01 J
  → Fracture energy : 2.1837e-04 J
  → Total energy    : 6.2796e+01 J


## Step 125/201: t = 6.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 124 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -496000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.789e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.857e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.761e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -496000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.967e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.7810e+01 J
  → Fracture energy : 2.2321e-04 J
  → Total energy    : 6.7811e+01 J


## Step 126/201: t = 6.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 125 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -500000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.694e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.650e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.389e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -500000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.133e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1739e+01 J
  → Fracture energy : 2.2515e-04 J
  → Total energy    : 7.1739e+01 J


## Step 127/201: t = 6.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 126 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -504000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.757e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.466e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.821e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -504000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.323e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.420e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4401e+01 J
  → Fracture energy : 2.2916e-04 J
  → Total energy    : 7.4401e+01 J


## Step 128/201: t = 6.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 127 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -508000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.962e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.215e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.640e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -508000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.338e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.958e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.7280e+01 J
  → Fracture energy : 2.3486e-04 J
  → Total energy    : 7.7281e+01 J


## Step 129/201: t = 6.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 128 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -512000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.071e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.033e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.302e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -512000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.281e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.198e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9715e+01 J
  → Fracture energy : 2.3882e-04 J
  → Total energy    : 7.9715e+01 J


## Step 130/201: t = 6.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 129 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -516000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.259e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.005e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.376e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -516000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.406e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.296e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.2642e+01 J
  → Fracture energy : 2.4432e-04 J
  → Total energy    : 8.2643e+01 J


## Step 131/201: t = 6.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 130 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -520000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.807e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.847e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.898e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -520000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.208e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.227e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5045e+01 J
  → Fracture energy : 2.5225e-04 J
  → Total energy    : 8.5045e+01 J


## Step 132/201: t = 6.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 131 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -524000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.440e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.583e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.288e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -524000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.744e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.352e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7310e+01 J
  → Fracture energy : 2.5542e-04 J
  → Total energy    : 8.7311e+01 J


## Step 133/201: t = 6.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 132 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -528000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.555e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.344e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.394e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -528000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.118e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.776e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0928e+01 J
  → Fracture energy : 2.6602e-04 J
  → Total energy    : 9.0928e+01 J


## Step 134/201: t = 6.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 133 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -532000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.869e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.245e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.137e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -532000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.325e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.627e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3077e+01 J
  → Fracture energy : 2.8217e-04 J
  → Total energy    : 9.3078e+01 J


## Step 135/201: t = 6.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 134 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -536000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.508e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.175e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.287e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -536000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.246e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.672e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.788e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6337e+01 J
  → Fracture energy : 2.9314e-04 J
  → Total energy    : 9.6337e+01 J


## Step 136/201: t = 6.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 135 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -540000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.406e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.242e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.041e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -540000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.034e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.640e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.998e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0242e+02 J
  → Fracture energy : 2.9793e-04 J
  → Total energy    : 1.0242e+02 J


## Step 137/201: t = 6.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 136 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -544000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.706e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.635e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.404e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -544000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.333e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.828e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.998e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0879e+02 J
  → Fracture energy : 3.0095e-04 J
  → Total energy    : 1.0879e+02 J


## Step 138/201: t = 6.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 137 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -548000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.919e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.314e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.980e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -548000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.798e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.683e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1298e+02 J
  → Fracture energy : 3.0240e-04 J
  → Total energy    : 1.1298e+02 J


## Step 139/201: t = 6.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 138 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -552000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.383e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.359e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.190e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -552000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.797e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.815e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1562e+02 J
  → Fracture energy : 3.0227e-04 J
  → Total energy    : 1.1562e+02 J


## Step 140/201: t = 6.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 139 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -556000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.163e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.648e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.227e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -556000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.992e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1775e+02 J
  → Fracture energy : 3.0243e-04 J
  → Total energy    : 1.1775e+02 J


## Step 141/201: t = 7.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 140 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -560000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.143e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.379e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.525e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -560000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.168e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.459e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1945e+02 J
  → Fracture energy : 3.0397e-04 J
  → Total energy    : 1.1945e+02 J


## Step 142/201: t = 7.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 141 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -564000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.093e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.510e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.491e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -564000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.865e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.477e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2116e+02 J
  → Fracture energy : 3.0609e-04 J
  → Total energy    : 1.2116e+02 J


## Step 143/201: t = 7.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 142 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -568000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.541e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.494e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.960e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -568000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.805e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.780e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2289e+02 J
  → Fracture energy : 3.0999e-04 J
  → Total energy    : 1.2289e+02 J


## Step 144/201: t = 7.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 143 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -572000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.706e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.318e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.749e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -572000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.128e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.579e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2527e+02 J
  → Fracture energy : 3.1470e-04 J
  → Total energy    : 1.2527e+02 J


## Step 145/201: t = 7.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 144 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -576000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.698e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.362e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.381e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -576000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.825e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.886e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2748e+02 J
  → Fracture energy : 3.1825e-04 J
  → Total energy    : 1.2748e+02 J


## Step 146/201: t = 7.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 145 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -580000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.152e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.628e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.565e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -580000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.103e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.012e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3078e+02 J
  → Fracture energy : 3.2091e-04 J
  → Total energy    : 1.3078e+02 J


## Step 147/201: t = 7.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 146 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -584000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.584e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.338e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.094e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -584000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.752e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.528e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.101e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3344e+02 J
  → Fracture energy : 3.2280e-04 J
  → Total energy    : 1.3344e+02 J


## Step 148/201: t = 7.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 147 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -588000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.050e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.439e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.467e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -588000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.452e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3642e+02 J
  → Fracture energy : 3.2462e-04 J
  → Total energy    : 1.3642e+02 J


## Step 149/201: t = 7.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 148 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.817e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.939e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.258e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.019e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.395e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3830e+02 J
  → Fracture energy : 3.2916e-04 J
  → Total energy    : 1.3830e+02 J


## Step 150/201: t = 7.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 149 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -596000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.836e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.879e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.883e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -596000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.552e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.699e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4017e+02 J
  → Fracture energy : 3.3851e-04 J
  → Total energy    : 1.4017e+02 J


## Step 151/201: t = 7.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 150 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -600000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.802e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.294e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.450e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -600000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.743e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.214e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4241e+02 J
  → Fracture energy : 3.4535e-04 J
  → Total energy    : 1.4241e+02 J


## Step 152/201: t = 7.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 151 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -604000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.784e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.924e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.208e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -604000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.977e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.107e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4600e+02 J
  → Fracture energy : 3.4939e-04 J
  → Total energy    : 1.4600e+02 J


## Step 153/201: t = 7.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 152 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -608000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.875e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.812e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.450e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -608000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.007e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.210e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5037e+02 J
  → Fracture energy : 3.5285e-04 J
  → Total energy    : 1.5037e+02 J


## Step 154/201: t = 7.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 153 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -612000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.840e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.380e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.453e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -612000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.201e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.655e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5340e+02 J
  → Fracture energy : 3.5744e-04 J
  → Total energy    : 1.5340e+02 J


## Step 155/201: t = 7.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 154 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -616000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.967e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.378e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.129e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -616000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.296e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.582e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5573e+02 J
  → Fracture energy : 3.6229e-04 J
  → Total energy    : 1.5573e+02 J


## Step 156/201: t = 7.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 155 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -620000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.088e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.250e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.677e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -620000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.891e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.254e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.982e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5914e+02 J
  → Fracture energy : 3.7099e-04 J
  → Total energy    : 1.5914e+02 J


## Step 157/201: t = 7.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 156 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -624000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.304e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.813e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.194e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -624000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.533e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.322e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6141e+02 J
  → Fracture energy : 3.7367e-04 J
  → Total energy    : 1.6141e+02 J


## Step 158/201: t = 7.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 157 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -628000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.926e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.143e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.761e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -628000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.654e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.084e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6775e+02 J
  → Fracture energy : 3.7396e-04 J
  → Total energy    : 1.6775e+02 J


## Step 159/201: t = 7.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 158 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -632000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.337e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.910e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.995e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -632000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.980e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.354e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.313e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6990e+02 J
  → Fracture energy : 3.7648e-04 J
  → Total energy    : 1.6990e+02 J


## Step 160/201: t = 7.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 159 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -636000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.290e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.799e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.588e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -636000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.102e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.099e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.109e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7205e+02 J
  → Fracture energy : 3.8200e-04 J
  → Total energy    : 1.7205e+02 J


## Step 161/201: t = 8.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 160 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -640000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.530e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.377e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.126e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -640000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.698e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7423e+02 J
  → Fracture energy : 3.8672e-04 J
  → Total energy    : 1.7423e+02 J


## Step 162/201: t = 8.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 161 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -644000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.737e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.227e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.129e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -644000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.220e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.039e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7771e+02 J
  → Fracture energy : 3.8942e-04 J
  → Total energy    : 1.7771e+02 J


## Step 163/201: t = 8.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 162 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -648000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.722e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.975e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.271e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -648000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.091e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.809e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8084e+02 J
  → Fracture energy : 3.9122e-04 J
  → Total energy    : 1.8084e+02 J


## Step 164/201: t = 8.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 163 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -652000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.552e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.749e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.422e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -652000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.893e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.656e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8376e+02 J
  → Fracture energy : 3.9668e-04 J
  → Total energy    : 1.8376e+02 J


## Step 165/201: t = 8.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 164 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -656000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.170e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.837e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.983e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -656000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.405e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8602e+02 J
  → Fracture energy : 3.9985e-04 J
  → Total energy    : 1.8602e+02 J


## Step 166/201: t = 8.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 165 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -660000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.274e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.157e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.537e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -660000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.184e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.138e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8985e+02 J
  → Fracture energy : 4.0423e-04 J
  → Total energy    : 1.8985e+02 J


## Step 167/201: t = 8.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 166 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -664000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.349e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.898e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.549e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -664000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.508e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.728e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9217e+02 J
  → Fracture energy : 4.0457e-04 J
  → Total energy    : 1.9217e+02 J


## Step 168/201: t = 8.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 167 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -668000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.228e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.659e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.460e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -668000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.764e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.469e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.998e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9648e+02 J
  → Fracture energy : 4.0456e-04 J
  → Total energy    : 1.9648e+02 J


## Step 169/201: t = 8.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 168 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -672000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.952e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.724e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.915e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -672000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.997e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.540e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9884e+02 J
  → Fracture energy : 4.0458e-04 J
  → Total energy    : 1.9884e+02 J


## Step 170/201: t = 8.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 169 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -676000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.917e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.990e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.751e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -676000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.360e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.959e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0121e+02 J
  → Fracture energy : 4.0491e-04 J
  → Total energy    : 2.0121e+02 J


## Step 171/201: t = 8.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 170 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -680000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.882e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.456e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.920e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -680000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.613e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.349e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0360e+02 J
  → Fracture energy : 4.0717e-04 J
  → Total energy    : 2.0360e+02 J


## Step 172/201: t = 8.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 171 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -684000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.849e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.876e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.410e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -684000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.972e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.344e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0600e+02 J
  → Fracture energy : 4.1087e-04 J
  → Total energy    : 2.0600e+02 J


## Step 173/201: t = 8.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 172 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -688000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.235e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.061e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.814e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -688000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.429e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.178e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0842e+02 J
  → Fracture energy : 4.1112e-04 J
  → Total energy    : 2.0842e+02 J


## Step 174/201: t = 8.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 173 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -692000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.634e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.067e-09
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.384e-08

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -692000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.571e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.945e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1257e+02 J
  → Fracture energy : 4.1112e-04 J
  → Total energy    : 2.1257e+02 J


## Step 175/201: t = 8.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 174 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -696000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.747e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.626e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.121e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -696000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.677e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.695e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1503e+02 J
  → Fracture energy : 4.1112e-04 J
  → Total energy    : 2.1503e+02 J


## Step 176/201: t = 8.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 175 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -700000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.714e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.912e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.562e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -700000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.088e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.954e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1751e+02 J
  → Fracture energy : 4.1111e-04 J
  → Total energy    : 2.1751e+02 J


## Step 177/201: t = 8.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 176 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -704000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.682e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.614e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.299e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -704000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.920e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.059e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2000e+02 J
  → Fracture energy : 4.1107e-04 J
  → Total energy    : 2.2000e+02 J


## Step 178/201: t = 8.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 177 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -708000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.650e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.188e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.204e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -708000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.826e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.665e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2251e+02 J
  → Fracture energy : 4.1099e-04 J
  → Total energy    : 2.2251e+02 J


## Step 179/201: t = 8.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 178 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -712000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.618e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.078e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.828e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -712000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.520e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.279e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2503e+02 J
  → Fracture energy : 4.1110e-04 J
  → Total energy    : 2.2503e+02 J


## Step 180/201: t = 8.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 179 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -716000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.587e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.226e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.718e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -716000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.669e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.099e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2757e+02 J
  → Fracture energy : 4.1279e-04 J
  → Total energy    : 2.2757e+02 J


## Step 181/201: t = 9.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 180 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -720000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.556e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.726e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.278e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -720000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.235e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.127e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3012e+02 J
  → Fracture energy : 4.1590e-04 J
  → Total energy    : 2.3012e+02 J


## Step 182/201: t = 9.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 181 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -724000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.699e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.275e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.002e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -724000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.223e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.061e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3268e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.3268e+02 J


## Step 183/201: t = 9.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 182 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -728000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.808e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.874e-09
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.649e-07

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -728000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.129e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.315e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3642e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.3642e+02 J


## Step 184/201: t = 9.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 183 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -732000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.464e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.478e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.648e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -732000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.365e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.267e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3902e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.3902e+02 J


## Step 185/201: t = 9.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 184 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -736000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.435e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.351e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.569e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -736000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.219e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.315e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4164e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.4164e+02 J


## Step 186/201: t = 9.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 185 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -740000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.405e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.232e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.494e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -740000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.216e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.843e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4427e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.4427e+02 J


## Step 187/201: t = 9.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 186 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -744000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.376e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.115e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.421e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -744000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.940e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.829e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4692e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.4692e+02 J


## Step 188/201: t = 9.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 187 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -748000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.348e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.998e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.350e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -748000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.229e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.076e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4958e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.4958e+02 J


## Step 189/201: t = 9.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 188 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -752000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.319e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.889e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.283e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -752000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.996e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.367e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5226e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.5226e+02 J


## Step 190/201: t = 9.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 189 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -756000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.291e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.781e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.214e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -756000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.012e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.319e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5495e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.5495e+02 J


## Step 191/201: t = 9.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 190 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -760000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.263e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.679e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.152e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -760000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.623e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.395e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5766e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.5766e+02 J


## Step 192/201: t = 9.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 191 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -764000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.236e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.574e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.088e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -764000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.528e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.560e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6038e+02 J
  → Fracture energy : 4.1614e-04 J
  → Total energy    : 2.6038e+02 J


## Step 193/201: t = 9.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 192 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -768000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.208e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.135e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.985e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -768000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.334e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.213e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.947e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6311e+02 J
  → Fracture energy : 4.1613e-04 J
  → Total energy    : 2.6311e+02 J


## Step 194/201: t = 9.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 193 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -772000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.181e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.530e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.106e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -772000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.104e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.773e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6586e+02 J
  → Fracture energy : 4.1612e-04 J
  → Total energy    : 2.6586e+02 J


## Step 195/201: t = 9.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 194 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -776000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.155e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.200e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.335e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -776000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.399e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.344e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6862e+02 J
  → Fracture energy : 4.1612e-04 J
  → Total energy    : 2.6862e+02 J


## Step 196/201: t = 9.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 195 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -780000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.128e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.254e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.104e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -780000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.261e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.542e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7140e+02 J
  → Fracture energy : 4.1628e-04 J
  → Total energy    : 2.7140e+02 J


## Step 197/201: t = 9.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 196 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -784000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.102e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.938e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.705e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -784000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.084e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.494e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7419e+02 J
  → Fracture energy : 4.1687e-04 J
  → Total energy    : 2.7419e+02 J


## Step 198/201: t = 9.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 197 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -788000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.076e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.781e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.207e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -788000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.783e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.571e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7699e+02 J
  → Fracture energy : 4.1836e-04 J
  → Total energy    : 2.7699e+02 J


## Step 199/201: t = 9.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 198 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -792000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.051e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.374e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.015e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -792000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.109e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.860e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7981e+02 J
  → Fracture energy : 4.2217e-04 J
  → Total energy    : 2.7981e+02 J


## Step 200/201: t = 9.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 199 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -796000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.025e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.490e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.580e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -796000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.576e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.004e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8264e+02 J
  → Fracture energy : 4.2807e-04 J
  → Total energy    : 2.8264e+02 J


## Step 201/201: t = 1.00e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 200 | dt = 5.00e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -800000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.080e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.960e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.205e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -800000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.062e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8549e+02 J
  → Fracture energy : 4.3021e-04 J
  → Total energy    : 2.8549e+02 J

Simulation completed in 152.86 s
Total time steps solved: 201
