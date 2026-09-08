[INFO] Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 17 entities
Info    : 1149 nodes
Info    : 2296 elements
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
[INFO]   Num cells: 2080
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
  ||Δu||/||u|| = 1.382e-17
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
  → Elastic energy  : 7.8436e-09 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.8436e-09 J


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
  ||Δu||/||u|| = 5.342e-17
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
  → Elastic energy  : 3.1374e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.1374e-08 J


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
  ||Δu||/||u|| = 7.106e-17
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
  → Elastic energy  : 7.0593e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.0593e-08 J


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
  ||Δu||/||u|| = 4.265e-17
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
  → Elastic energy  : 1.2550e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2550e-07 J


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
  ||Δu||/||u|| = 3.553e-17
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
  → Elastic energy  : 1.9609e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.9609e-07 J


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
  → Elastic energy  : 2.8237e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8237e-07 J


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
  ||Δu||/||u|| = 7.250e-17
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
  → Elastic energy  : 3.8434e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.8434e-07 J


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
  ||Δu||/||u|| = 3.350e-17
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
  → Elastic energy  : 5.0199e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.0199e-07 J


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
  ||Δu||/||u|| = 8.284e-17
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
  → Elastic energy  : 6.3533e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.3533e-07 J


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
  ||Δu||/||u|| = 1.332e-16
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
  → Elastic energy  : 7.8436e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.8436e-07 J


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
  ||Δu||/||u|| = 1.201e-18
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
  → Elastic energy  : 9.4908e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.4908e-07 J


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
  ||Δu||/||u|| = 8.206e-17
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
  → Elastic energy  : 1.1295e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1295e-06 J


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
  ||Δu||/||u|| = 1.063e-16
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
  → Elastic energy  : 1.3256e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3256e-06 J


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
  ||Δu||/||u|| = 1.355e-16
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
  → Elastic energy  : 1.5373e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5373e-06 J


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
  ||Δu||/||u|| = 2.520e-17
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
  → Elastic energy  : 1.7648e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.7648e-06 J


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
  ||Δu||/||u|| = 7.887e-18
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
  → Elastic energy  : 2.0080e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0080e-06 J


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
  ||Δu||/||u|| = 1.864e-16
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
  → Elastic energy  : 2.2668e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2668e-06 J


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
  ||Δu||/||u|| = 4.153e-17
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
  → Elastic energy  : 2.5413e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5413e-06 J


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
  ||Δu||/||u|| = 1.470e-17
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
  → Elastic energy  : 2.8315e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8315e-06 J


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
  ||Δu||/||u|| = 1.819e-18
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
  → Elastic energy  : 3.1374e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.1374e-06 J


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
  ||Δu||/||u|| = 8.662e-17
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
  → Elastic energy  : 3.4590e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.4590e-06 J


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
  ||Δu||/||u|| = 1.736e-16
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
  → Elastic energy  : 3.7963e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.7963e-06 J


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
  ||Δu||/||u|| = 4.979e-17
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
  → Elastic energy  : 4.1493e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.1493e-06 J


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
  → Elastic energy  : 4.5179e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.5179e-06 J


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
  ||Δu||/||u|| = 3.983e-17
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
  → Elastic energy  : 4.9023e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.9023e-06 J


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
  ||Δu||/||u|| = 1.998e-16
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
  → Elastic energy  : 5.3023e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.3023e-06 J


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
  ||Δu||/||u|| = 5.419e-18
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
  → Elastic energy  : 5.7180e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.7180e-06 J


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
  → Elastic energy  : 6.1494e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.1494e-06 J


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
  ||Δu||/||u|| = 7.141e-17
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
  → Elastic energy  : 6.5965e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.5965e-06 J


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
  ||Δu||/||u|| = 8.853e-17
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
  → Elastic energy  : 7.0593e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.0593e-06 J


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
  ||Δu||/||u|| = 1.393e-16
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
  → Elastic energy  : 7.5377e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.5377e-06 J


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
  ||Δu||/||u|| = 5.633e-18
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
  → Elastic energy  : 8.0319e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.0319e-06 J


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
  ||Δu||/||u|| = 8.123e-17
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
  → Elastic energy  : 8.5417e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.5417e-06 J


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
  ||Δu||/||u|| = 7.309e-17
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
  → Elastic energy  : 9.0672e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.0672e-06 J


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
  → Elastic energy  : 9.6084e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.6084e-06 J


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
  ||Δu||/||u|| = 9.856e-17
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
  → Elastic energy  : 1.0165e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0165e-05 J


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
  ||Δu||/||u|| = 1.285e-16
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
  → Elastic energy  : 1.0738e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0738e-05 J


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
  ||Δu||/||u|| = 2.407e-17
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
  → Elastic energy  : 1.1326e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1326e-05 J


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
  ||Δu||/||u|| = 6.745e-17
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
  → Elastic energy  : 1.1930e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1930e-05 J


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
  ||Δu||/||u|| = 1.821e-17
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
  → Elastic energy  : 1.2550e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2550e-05 J


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
  → Elastic energy  : 1.3185e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3185e-05 J


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
  ||Δu||/||u|| = 3.804e-17
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
  → Elastic energy  : 1.3836e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3836e-05 J


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
  ||Δu||/||u|| = 1.041e-16
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
  → Elastic energy  : 1.4503e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4503e-05 J


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
  ||Δu||/||u|| = 1.140e-16
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
  → Elastic energy  : 1.5185e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5185e-05 J


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
  ||Δu||/||u|| = 1.189e-17
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
  → Elastic energy  : 1.5883e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5883e-05 J


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
  ||Δu||/||u|| = 7.283e-17
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
  → Elastic energy  : 1.6597e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.6597e-05 J


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
  → Elastic energy  : 1.7327e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.7327e-05 J


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
  ||Δu||/||u|| = 2.070e-17
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
  → Elastic energy  : 1.8072e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.8072e-05 J


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
  ||Δu||/||u|| = 6.442e-18
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
  → Elastic energy  : 1.8833e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.8833e-05 J


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
  ||Δu||/||u|| = 1.654e-17
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
  → Elastic energy  : 1.9609e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.9609e-05 J


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
  ||Δu||/||u|| = 1.255e-17
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
  → Elastic energy  : 2.0401e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0401e-05 J


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
  → Elastic energy  : 2.1209e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.1209e-05 J


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
  → Elastic energy  : 2.2033e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2033e-05 J


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
  ||Δu||/||u|| = 2.766e-17
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
  → Elastic energy  : 2.2872e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2872e-05 J


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
  ||Δu||/||u|| = 3.710e-17
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
  → Elastic energy  : 2.3727e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.3727e-05 J


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
  ||Δu||/||u|| = 1.286e-16
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
  → Elastic energy  : 2.4598e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.4598e-05 J


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
  ||Δu||/||u|| = 9.196e-18
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
  → Elastic energy  : 2.5484e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5484e-05 J


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
  → Elastic energy  : 2.6386e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.6386e-05 J


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
  → Elastic energy  : 2.7304e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.7304e-05 J


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
  ||Δu||/||u|| = 9.498e-17
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
  → Elastic energy  : 2.8237e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8237e-05 J


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
  → Elastic energy  : 2.9186e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.9186e-05 J


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
  ||Δu||/||u|| = 1.927e-16
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
  → Elastic energy  : 3.0151e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.0151e-05 J


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
  ||Δu||/||u|| = 5.202e-17
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
  → Elastic energy  : 3.1131e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.1131e-05 J


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
  ||Δu||/||u|| = 8.430e-18
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
  → Elastic energy  : 3.2127e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.2127e-05 J


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
  ||Δu||/||u|| = 1.287e-17
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
  → Elastic energy  : 3.3139e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.3139e-05 J


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
  ||Δu||/||u|| = 7.476e-17
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
  → Elastic energy  : 3.4167e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.4167e-05 J


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
  ||Δu||/||u|| = 1.503e-16
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
  → Elastic energy  : 3.5210e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.5210e-05 J


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
  ||Δu||/||u|| = 7.609e-17
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
  → Elastic energy  : 3.6269e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.6269e-05 J


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
  ||Δu||/||u|| = 9.124e-18
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
  → Elastic energy  : 3.7343e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.7343e-05 J


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
  → Elastic energy  : 3.8434e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.8434e-05 J


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
  ||Δu||/||u|| = 1.120e-16
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
  → Elastic energy  : 3.9540e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.9540e-05 J


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
  ||Δu||/||u|| = 9.362e-17
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
  → Elastic energy  : 4.0661e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.0661e-05 J


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
  ||Δu||/||u|| = 2.578e-16
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
  → Elastic energy  : 4.1799e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.1799e-05 J


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
  ||Δu||/||u|| = 1.969e-17
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
  → Elastic energy  : 4.2952e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.2952e-05 J


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
  → Elastic energy  : 4.4120e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.4120e-05 J


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
  ||Δu||/||u|| = 1.111e-16
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
  → Elastic energy  : 4.5305e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.5305e-05 J


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
  → Elastic energy  : 4.6505e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.6505e-05 J


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
  ||Δu||/||u|| = 7.422e-17
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
  → Elastic energy  : 4.7721e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.7721e-05 J


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
  ||Δu||/||u|| = 5.869e-17
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
  → Elastic energy  : 4.8952e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.8952e-05 J


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
  ||Δu||/||u|| = 1.128e-17
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
  → Elastic energy  : 5.0199e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.0199e-05 J


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
  ||Δu||/||u|| = 5.643e-17
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
  → Elastic energy  : 5.1462e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.1462e-05 J


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
  → Elastic energy  : 5.2741e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.2741e-05 J


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
  ||Δu||/||u|| = 8.410e-17
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
  → Elastic energy  : 5.4035e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.4035e-05 J


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
  ||Δu||/||u|| = 8.778e-17
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
  → Elastic energy  : 5.5345e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.5345e-05 J


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
  → Elastic energy  : 5.6670e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.6670e-05 J


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
  ||Δu||/||u|| = 1.082e-16
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
  → Elastic energy  : 5.8011e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.8011e-05 J


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
  ||Δu||/||u|| = 7.464e-17
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
  → Elastic energy  : 5.9368e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.9368e-05 J


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
  ||Δu||/||u|| = 2.489e-17
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
  → Elastic energy  : 6.0741e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.0741e-05 J


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
  → Elastic energy  : 6.2129e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.2129e-05 J


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
  ||Δu||/||u|| = 1.356e-17
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
  → Elastic energy  : 6.3533e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.3533e-05 J


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
  ||Δu||/||u|| = 1.269e-16
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
  → Elastic energy  : 6.4953e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.4953e-05 J


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
  ||Δu||/||u|| = 4.980e-17
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
  → Elastic energy  : 6.6388e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.6388e-05 J


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
  ||Δu||/||u|| = 1.426e-16
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
  → Elastic energy  : 6.7839e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.7839e-05 J


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
  ||Δu||/||u|| = 9.382e-17
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
  → Elastic energy  : 6.9306e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.9306e-05 J


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
  ||Δu||/||u|| = 3.237e-18
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
  → Elastic energy  : 7.0789e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.0789e-05 J


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
  ||Δu||/||u|| = 2.059e-17
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
  → Elastic energy  : 7.2287e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.2287e-05 J


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
  → Elastic energy  : 7.3801e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.3801e-05 J


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
  ||Δu||/||u|| = 1.218e-16
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
  → Elastic energy  : 7.5330e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.5330e-05 J


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
  ||Δu||/||u|| = 1.299e-19
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
  → Elastic energy  : 7.6875e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.6875e-05 J


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
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

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
  ||Δu||/||u|| = 8.588e-19
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
  → Elastic energy  : 7.8436e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.8436e-05 J


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
  ||Δu||/||u|| = 9.901e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 7.768e-04

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
  ||Δu||/||u|| = 1.344e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.500e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 7.379e-04

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
  ||Δu||/||u|| = 9.858e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.927e-01
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 7.712e-04

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
  ||Δu||/||u|| = 2.862e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.529e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 6.625e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.817e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.912e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 6.923e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.773e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.264e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 7.197e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.001e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.574e-01
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 7.437e-04

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.358e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.715e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 5.216e-04

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.815e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.017e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 5.451e-04

Convergence check


#### Iteration 10/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.815e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.294e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 5.666e-04

Convergence check


#### Iteration 11/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.621e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.538e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 5.855e-04

Convergence check


#### Iteration 12/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.006e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.740e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 6.012e-04

Convergence check


#### Iteration 13/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.161e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.891e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 6.129e-04

Convergence check


#### Iteration 14/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.164e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.079e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 3.945e-04

Convergence check


#### Iteration 15/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.394e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.272e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 4.095e-04

Convergence check


#### Iteration 16/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.815e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.439e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 4.225e-04

Convergence check


#### Iteration 17/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.725e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.575e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 4.331e-04

Convergence check


#### Iteration 18/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.819e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.673e-01
  [adaptive] relax_D=0.09
  |ΔD|_∞ = 4.406e-04

Convergence check


#### Iteration 19/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.817e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.725e-01
  [adaptive] relax_D=0.10
  |ΔD|_∞ = 4.447e-04

Convergence check


#### Iteration 20/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.057e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.726e-01
  [adaptive] relax_D=0.11
  |ΔD|_∞ = 4.448e-04

Convergence check


#### Iteration 21/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.057e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.669e-01
  [adaptive] relax_D=0.12
  |ΔD|_∞ = 4.404e-04

Convergence check


#### Iteration 22/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.606e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.551e-01
  [adaptive] relax_D=0.13
  |ΔD|_∞ = 4.312e-04

Convergence check


#### Iteration 23/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.993e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.368e-01
  [adaptive] relax_D=0.15
  |ΔD|_∞ = 4.170e-04

Convergence check


#### Iteration 24/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.770e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.120e-01
  [adaptive] relax_D=0.16
  |ΔD|_∞ = 3.977e-04

Convergence check


#### Iteration 25/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.620e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.809e-01
  [adaptive] relax_D=0.18
  |ΔD|_∞ = 3.736e-04

Convergence check


#### Iteration 26/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.006e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.439e-01
  [adaptive] relax_D=0.19
  |ΔD|_∞ = 3.448e-04

Convergence check


#### Iteration 27/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.606e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.019e-01
  [adaptive] relax_D=0.21
  |ΔD|_∞ = 3.122e-04

Convergence check


#### Iteration 28/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.561e-01
  [adaptive] relax_D=0.24
  |ΔD|_∞ = 2.766e-04

Convergence check


#### Iteration 29/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.606e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.078e-01
  [adaptive] relax_D=0.26
  |ΔD|_∞ = 2.391e-04

Convergence check


#### Iteration 30/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.589e-01
  [adaptive] relax_D=0.28
  |ΔD|_∞ = 2.011e-04

Convergence check


#### Iteration 31/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.817e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.110e-01
  [adaptive] relax_D=0.31
  |ΔD|_∞ = 1.639e-04

Convergence check


#### Iteration 32/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.337e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.660e-01
  [adaptive] relax_D=0.34
  |ΔD|_∞ = 1.289e-04

Convergence check


#### Iteration 33/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.337e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.254e-01
  [adaptive] relax_D=0.38
  |ΔD|_∞ = 9.738e-05

Convergence check


#### Iteration 34/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.036e-02
  [adaptive] relax_D=0.42
  |ΔD|_∞ = 7.019e-05

Convergence check


#### Iteration 35/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.815e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.170e-02
  [adaptive] relax_D=0.46
  |ΔD|_∞ = 4.793e-05

Convergence check


#### Iteration 36/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.862e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.956e-02
  [adaptive] relax_D=0.50
  |ΔD|_∞ = 3.073e-05

Convergence check


#### Iteration 37/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.394e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.355e-02
  [adaptive] relax_D=0.56
  |ΔD|_∞ = 1.829e-05

Convergence check


#### Iteration 38/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.394e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.283e-02
  [adaptive] relax_D=0.61
  |ΔD|_∞ = 9.966e-06

Convergence check


#### Iteration 39/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.277e-03
  [adaptive] relax_D=0.67
  |ΔD|_∞ = 4.876e-06

Convergence check


#### Iteration 40/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.394e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.688e-03
  [adaptive] relax_D=0.74
  |ΔD|_∞ = 2.088e-06

Convergence check


#### Iteration 41/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.704e-04
  [adaptive] relax_D=0.81
  |ΔD|_∞ = 7.538e-07

Convergence check


#### Iteration 42/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.786e-04
  [adaptive] relax_D=0.89
  |ΔD|_∞ = 2.164e-07

Convergence check


#### Iteration 43/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -404000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.735e-05
  [adaptive] relax_D=0.98
  |ΔD|_∞ = 4.455e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 43 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9990e-05 J
  → Fracture energy : 1.9669e-09 J
  → Total energy    : 7.9992e-05 J


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
  ||Δu||/||u|| = 1.007e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.021e-01
  [adaptive] relax_D=0.98
  |ΔD|_∞ = 4.737e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.317e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.776e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -408000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.014e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.197e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.298e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -408000000.0 Pa
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

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1533e-05 J
  → Fracture energy : 3.7243e-08 J
  → Total energy    : 8.1570e-05 J


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
  ||Δu||/||u|| = 1.086e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.868e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.929e-02

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
  ||Δu||/||u|| = 2.448e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.667e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.290e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.3092e-05 J
  → Fracture energy : 1.0870e-07 J
  → Total energy    : 8.3200e-05 J


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
  ||Δu||/||u|| = 1.266e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.175e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.100e-02

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
  → Elastic energy  : 8.4781e-05 J
  → Fracture energy : 2.9999e-07 J
  → Total energy    : 8.5081e-05 J


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
  ||Δu||/||u|| = 1.515e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.739e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.037e-01

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
  ||Δu||/||u|| = 4.804e-19
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
  → Elastic energy  : 8.6552e-05 J
  → Fracture energy : 6.3781e-07 J
  → Total energy    : 8.7189e-05 J


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
  ||Δu||/||u|| = 1.983e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.525e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.374e-01

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
  ||Δu||/||u|| = 1.176e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.106e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.8490e-05 J
  → Fracture energy : 1.2866e-06 J
  → Total energy    : 8.9776e-05 J


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
  ||Δu||/||u|| = 2.825e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.474e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.878e-01

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
  ||Δu||/||u|| = 9.692e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.275e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0710e-05 J
  → Fracture energy : 2.4056e-06 J
  → Total energy    : 9.3115e-05 J


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
  ||Δu||/||u|| = 4.389e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.389e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.581e-01

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
  ||Δu||/||u|| = 8.125e-18
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
  → Elastic energy  : 9.3644e-05 J
  → Fracture energy : 4.2532e-06 J
  → Total energy    : 9.7897e-05 J


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
  ||Δu||/||u|| = 7.237e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.284e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.506e-01

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
  ||Δu||/||u|| = 2.399e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.654e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8477e-05 J
  → Fracture energy : 6.9683e-06 J
  → Total energy    : 1.0545e-04 J


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
  ||Δu||/||u|| = 1.957e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.172e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.391e-01

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
  ||Δu||/||u|| = 3.041e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.640e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.483e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0662e-04 J
  → Fracture energy : 1.1014e-05 J
  → Total energy    : 1.1764e-04 J


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
  ||Δu||/||u|| = 9.996e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.008e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.021e-01

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
  ||Δu||/||u|| = 3.822e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.041e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.759e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.5747e-03 J
  → Fracture energy : 1.7563e-05 J
  → Total energy    : 6.5923e-03 J


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
  ||Δu||/||u|| = 9.191e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.973e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.525e-01

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
  ||Δu||/||u|| = 1.250e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.873e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3207e-01 J
  → Fracture energy : 2.8910e-05 J
  → Total energy    : 1.3210e-01 J


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
  ||Δu||/||u|| = 7.907e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.559e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.277e-01

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
  ||Δu||/||u|| = 1.177e-20
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.492e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7361e-01 J
  → Fracture energy : 5.5776e-05 J
  → Total energy    : 4.7366e-01 J


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
  ||Δu||/||u|| = 6.332e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.553e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.920e-01

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
  ||Δu||/||u|| = 9.532e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.852e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1007e+00 J
  → Fracture energy : 1.0013e-04 J
  → Total energy    : 1.1008e+00 J


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
  ||Δu||/||u|| = 5.479e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.763e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.350e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.348e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1285e+00 J
  → Fracture energy : 1.3746e-04 J
  → Total energy    : 2.1286e+00 J


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
  ||Δu||/||u|| = 9.732e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.532e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.789e-01

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
  ||ΔD||/||D|| = 6.108e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4587e+01 J
  → Fracture energy : 1.5846e-04 J
  → Total energy    : 2.4587e+01 J


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
  ||Δu||/||u|| = 2.777e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.404e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.505e-01

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
  ||Δu||/||u|| = 1.017e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.184e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4996e+01 J
  → Fracture energy : 1.7520e-04 J
  → Total energy    : 3.4997e+01 J


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
  ||Δu||/||u|| = 1.424e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.281e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.444e-01

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
  ||Δu||/||u|| = 6.017e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.170e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1084e+01 J
  → Fracture energy : 1.8361e-04 J
  → Total energy    : 4.1084e+01 J


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
  ||Δu||/||u|| = 1.198e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.268e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.908e-01

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
  ||Δu||/||u|| = 8.487e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.843e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7158e+01 J
  → Fracture energy : 1.9254e-04 J
  → Total energy    : 4.7158e+01 J


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
  ||Δu||/||u|| = 6.878e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.327e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.240e-01

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
  ||Δu||/||u|| = 7.338e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.814e-11
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.822e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0579e+01 J
  → Fracture energy : 2.0332e-04 J
  → Total energy    : 5.0580e+01 J


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
  ||Δu||/||u|| = 6.174e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.399e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.382e-01

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
  ||Δu||/||u|| = 3.167e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.361e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3964e+01 J
  → Fracture energy : 2.1343e-04 J
  → Total energy    : 5.3964e+01 J


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
  ||Δu||/||u|| = 8.302e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.038e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.235e-01

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
  ||Δu||/||u|| = 2.661e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.634e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.682e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.8857e+01 J
  → Fracture energy : 2.2284e-04 J
  → Total energy    : 5.8857e+01 J


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
  ||Δu||/||u|| = 6.337e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.665e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.790e-01

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
  ||Δu||/||u|| = 8.729e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.078e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2975e+01 J
  → Fracture energy : 2.3614e-04 J
  → Total energy    : 6.2976e+01 J


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
  ||Δu||/||u|| = 5.802e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.316e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.415e-01

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
  ||Δu||/||u|| = 2.612e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.961e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.160e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.7390e+01 J
  → Fracture energy : 2.4730e-04 J
  → Total energy    : 6.7390e+01 J


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
  ||Δu||/||u|| = 8.588e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.063e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.292e-01

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
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3303e+01 J
  → Fracture energy : 2.5392e-04 J
  → Total energy    : 7.3303e+01 J


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
  ||Δu||/||u|| = 7.821e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.100e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.455e-01

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
  ||Δu||/||u|| = 1.964e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.168e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.8944e+01 J
  → Fracture energy : 2.6145e-04 J
  → Total energy    : 7.8945e+01 J


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
  ||Δu||/||u|| = 4.946e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.277e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.571e-01

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
  → Elastic energy  : 8.3017e+01 J
  → Fracture energy : 2.6449e-04 J
  → Total energy    : 8.3018e+01 J


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
  ||Δu||/||u|| = 6.508e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.818e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.155e-01

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
  ||Δu||/||u|| = 3.189e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.059e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7891e+01 J
  → Fracture energy : 2.6802e-04 J
  → Total energy    : 8.7891e+01 J


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
  ||Δu||/||u|| = 3.672e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.750e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.902e-01

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
  ||Δu||/||u|| = 8.791e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.256e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0536e+01 J
  → Fracture energy : 2.7321e-04 J
  → Total energy    : 9.0536e+01 J


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
  ||Δu||/||u|| = 2.387e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.823e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.036e-01

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
  → Elastic energy  : 9.3049e+01 J
  → Fracture energy : 2.7605e-04 J
  → Total energy    : 9.3049e+01 J


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
  ||Δu||/||u|| = 3.448e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.516e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.319e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.576e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6124e+01 J
  → Fracture energy : 2.7667e-04 J
  → Total energy    : 9.6125e+01 J


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
  ||Δu||/||u|| = 1.771e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.421e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.225e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.757e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8258e+01 J
  → Fracture energy : 2.7833e-04 J
  → Total energy    : 9.8258e+01 J


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
  ||Δu||/||u|| = 1.208e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.936e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.330e-01

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
  ||Δu||/||u|| = 6.323e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.788e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.9927e+01 J
  → Fracture energy : 2.8274e-04 J
  → Total energy    : 9.9927e+01 J


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
  ||Δu||/||u|| = 7.962e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.829e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.851e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.856e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0150e+02 J
  → Fracture energy : 2.8650e-04 J
  → Total energy    : 1.0150e+02 J


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
  ||Δu||/||u|| = 2.153e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.089e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.121e-01

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
  ||Δu||/||u|| = 2.041e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.597e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0394e+02 J
  → Fracture energy : 2.8832e-04 J
  → Total energy    : 1.0394e+02 J


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
  ||Δu||/||u|| = 2.084e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.058e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.079e-01

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
  ||Δu||/||u|| = 1.232e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.171e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0607e+02 J
  → Fracture energy : 2.8863e-04 J
  → Total energy    : 1.0607e+02 J


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
  ||Δu||/||u|| = 2.894e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.750e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.332e-01

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
  → Elastic energy  : 1.0887e+02 J
  → Fracture energy : 2.9214e-04 J
  → Total energy    : 1.0887e+02 J


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
  ||Δu||/||u|| = 7.252e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.008e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.228e-01

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
  ||Δu||/||u|| = 2.822e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.826e-11
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.113e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1046e+02 J
  → Fracture energy : 2.9710e-04 J
  → Total energy    : 1.1046e+02 J


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
  ||Δu||/||u|| = 2.254e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.414e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.135e-01

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
  ||Δu||/||u|| = 1.253e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.394e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1286e+02 J
  → Fracture energy : 3.0208e-04 J
  → Total energy    : 1.1286e+02 J


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
  ||Δu||/||u|| = 1.089e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.679e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.799e-01

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
  ||Δu||/||u|| = 8.189e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.761e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1459e+02 J
  → Fracture energy : 3.1004e-04 J
  → Total energy    : 1.1459e+02 J


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
  ||Δu||/||u|| = 2.582e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.324e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.522e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.591e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1775e+02 J
  → Fracture energy : 3.2098e-04 J
  → Total energy    : 1.1775e+02 J


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
  ||Δu||/||u|| = 2.967e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.261e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.618e-01

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
  ||Δu||/||u|| = 1.119e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.340e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2116e+02 J
  → Fracture energy : 3.2910e-04 J
  → Total energy    : 1.2116e+02 J


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
  ||Δu||/||u|| = 3.031e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.519e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.631e-01

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
  ||Δu||/||u|| = 1.019e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.934e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2520e+02 J
  → Fracture energy : 3.3680e-04 J
  → Total energy    : 1.2520e+02 J


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
  ||Δu||/||u|| = 2.523e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.191e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.354e-01

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
  ||Δu||/||u|| = 1.340e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.866e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2852e+02 J
  → Fracture energy : 3.4473e-04 J
  → Total energy    : 1.2852e+02 J


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
  ||Δu||/||u|| = 2.105e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.461e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.227e-01

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
  ||Δu||/||u|| = 2.201e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.374e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3164e+02 J
  → Fracture energy : 3.4819e-04 J
  → Total energy    : 1.3164e+02 J


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
  ||Δu||/||u|| = 6.606e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.208e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.269e-01

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
  ||Δu||/||u|| = 2.043e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.501e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3877e+02 J
  → Fracture energy : 3.5143e-04 J
  → Total energy    : 1.3877e+02 J


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
  ||Δu||/||u|| = 1.857e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.254e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.703e-01

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
  ||ΔD||/||D|| = 7.093e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4120e+02 J
  → Fracture energy : 3.5462e-04 J
  → Total energy    : 1.4120e+02 J


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
  ||Δu||/||u|| = 1.014e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.073e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.328e-01

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
  ||Δu||/||u|| = 1.498e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.531e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4357e+02 J
  → Fracture energy : 3.5653e-04 J
  → Total energy    : 1.4357e+02 J


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
  ||Δu||/||u|| = 2.029e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.452e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.549e-01

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
  ||Δu||/||u|| = 7.248e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.498e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4632e+02 J
  → Fracture energy : 3.6093e-04 J
  → Total energy    : 1.4632e+02 J


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
  ||Δu||/||u|| = 6.693e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.622e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.213e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.278e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4829e+02 J
  → Fracture energy : 3.6203e-04 J
  → Total energy    : 1.4829e+02 J


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
  ||Δu||/||u|| = 2.940e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.507e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.098e-01

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
  ||Δu||/||u|| = 2.911e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.545e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5213e+02 J
  → Fracture energy : 3.6200e-04 J
  → Total energy    : 1.5213e+02 J


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
  ||Δu||/||u|| = 6.580e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.176e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.561e-01

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
  ||Δu||/||u|| = 3.899e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.240e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5415e+02 J
  → Fracture energy : 3.6225e-04 J
  → Total energy    : 1.5415e+02 J


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
  ||Δu||/||u|| = 6.536e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.386e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.343e-01

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
  ||Δu||/||u|| = 1.119e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.243e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5618e+02 J
  → Fracture energy : 3.6371e-04 J
  → Total energy    : 1.5619e+02 J


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
  ||Δu||/||u|| = 6.495e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.319e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.980e-01

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
  ||Δu||/||u|| = 1.143e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.580e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5823e+02 J
  → Fracture energy : 3.6481e-04 J
  → Total energy    : 1.5823e+02 J


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
  ||Δu||/||u|| = 1.381e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.535e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.735e-03

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
  ||Δu||/||u|| = 2.355e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.555e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6068e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.6068e+02 J


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
  ||Δu||/||u|| = 7.937e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.334e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.285e-11

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
  ||Δu||/||u|| = 3.582e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.429e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6301e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.6301e+02 J


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
  ||Δu||/||u|| = 6.369e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.698e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.599e-12

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
  ||Δu||/||u|| = 4.968e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.471e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6510e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.6510e+02 J


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
  ||Δu||/||u|| = 6.329e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.638e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.479e-12

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
  ||Δu||/||u|| = 4.873e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.267e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6721e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.6721e+02 J


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
  ||Δu||/||u|| = 6.289e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.578e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.359e-12

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
  → Elastic energy  : 1.6934e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.6934e+02 J


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
  ||Δu||/||u|| = 6.250e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.523e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.245e-12

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
  ||Δu||/||u|| = 2.795e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.228e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7147e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.7147e+02 J


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
  ||Δu||/||u|| = 6.211e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.469e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.136e-12

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
  ||Δu||/||u|| = 8.041e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.185e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7362e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.7362e+02 J


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
  ||Δu||/||u|| = 6.173e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.417e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.029e-12

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
  ||Δu||/||u|| = 3.268e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.282e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7579e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.7579e+02 J


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
  ||Δu||/||u|| = 6.135e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.367e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.927e-12

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
  → Elastic energy  : 1.7796e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.7796e+02 J


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
  ||Δu||/||u|| = 6.098e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.319e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.829e-12

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
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8015e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.8015e+02 J


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
  ||Δu||/||u|| = 6.061e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.273e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.735e-12

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
  ||Δu||/||u|| = 5.226e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.702e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8236e+02 J
  → Fracture energy : 3.6482e-04 J
  → Total energy    : 1.8236e+02 J


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
  ||Δu||/||u|| = 6.024e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.277e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.654e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.771e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8457e+02 J
  → Fracture energy : 3.6481e-04 J
  → Total energy    : 1.8457e+02 J


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
  ||Δu||/||u|| = 5.988e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.857e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.820e-02

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
  ||Δu||/||u|| = 1.013e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.262e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8680e+02 J
  → Fracture energy : 3.6481e-04 J
  → Total energy    : 1.8680e+02 J


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
  ||ΔD||/||D|| = 4.880e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.529e-01

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
  ||Δu||/||u|| = 5.605e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.381e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8905e+02 J
  → Fracture energy : 3.6502e-04 J
  → Total energy    : 1.8905e+02 J


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
  ||ΔD||/||D|| = 1.052e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.296e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.874e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9131e+02 J
  → Fracture energy : 3.6661e-04 J
  → Total energy    : 1.9131e+02 J


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
  ||Δu||/||u|| = 5.883e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.246e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.858e-01

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
  ||Δu||/||u|| = 1.827e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.265e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9358e+02 J
  → Fracture energy : 3.7037e-04 J
  → Total energy    : 1.9358e+02 J


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
  ||Δu||/||u|| = 5.868e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.223e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.531e-01

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
  ||Δu||/||u|| = 2.310e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.497e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9586e+02 J
  → Fracture energy : 3.7116e-04 J
  → Total energy    : 1.9586e+02 J


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
  ||Δu||/||u|| = 1.672e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.844e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.105e-01

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
  ||Δu||/||u|| = 1.354e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.133e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9915e+02 J
  → Fracture energy : 3.7340e-04 J
  → Total energy    : 1.9915e+02 J


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
  ||Δu||/||u|| = 5.782e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.472e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.267e-01

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
  ||Δu||/||u|| = 2.050e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.202e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0147e+02 J
  → Fracture energy : 3.7784e-04 J
  → Total energy    : 2.0147e+02 J


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
  ||Δu||/||u|| = 1.362e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.546e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.511e-01

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
  → Elastic energy  : 2.0443e+02 J
  → Fracture energy : 3.8146e-04 J
  → Total energy    : 2.0443e+02 J


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
  ||Δu||/||u|| = 1.312e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.687e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.904e-01

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
  ||Δu||/||u|| = 4.951e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.976e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0753e+02 J
  → Fracture energy : 3.8658e-04 J
  → Total energy    : 2.0753e+02 J


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
  ||Δu||/||u|| = 6.247e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.117e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.504e-01

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
  ||Δu||/||u|| = 4.582e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.028e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1007e+02 J
  → Fracture energy : 3.9199e-04 J
  → Total energy    : 2.1008e+02 J


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
  ||Δu||/||u|| = 1.490e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.852e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.219e-01

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
  ||Δu||/||u|| = 3.836e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.576e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1381e+02 J
  → Fracture energy : 3.9692e-04 J
  → Total energy    : 2.1381e+02 J


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
  ||Δu||/||u|| = 1.695e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.746e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.586e-01

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
  ||Δu||/||u|| = 6.635e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.099e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1764e+02 J
  → Fracture energy : 4.0369e-04 J
  → Total energy    : 2.1764e+02 J


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
  ||Δu||/||u|| = 6.042e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.489e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.003e-01

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
  ||Δu||/||u|| = 9.040e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.088e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2011e+02 J
  → Fracture energy : 4.1085e-04 J
  → Total energy    : 2.2011e+02 J


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
  ||Δu||/||u|| = 1.835e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.261e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.167e-01

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
  ||Δu||/||u|| = 1.408e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.052e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2423e+02 J
  → Fracture energy : 4.2003e-04 J
  → Total energy    : 2.2423e+02 J


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
  ||Δu||/||u|| = 1.919e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.154e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.085e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.271e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2825e+02 J
  → Fracture energy : 4.2363e-04 J
  → Total energy    : 2.2825e+02 J


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
  ||Δu||/||u|| = 3.488e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.274e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.001e-01

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
  ||Δu||/||u|| = 7.106e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.523e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3456e+02 J
  → Fracture energy : 4.2573e-04 J
  → Total energy    : 2.3456e+02 J


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
  ||Δu||/||u|| = 9.111e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.189e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.617e-01

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
  ||Δu||/||u|| = 2.880e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.319e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3817e+02 J
  → Fracture energy : 4.2951e-04 J
  → Total energy    : 2.3817e+02 J


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
  ||Δu||/||u|| = 5.574e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.510e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.355e-01

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
  ||Δu||/||u|| = 1.216e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.048e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4078e+02 J
  → Fracture energy : 4.3034e-04 J
  → Total energy    : 2.4078e+02 J


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
  ||Δu||/||u|| = 2.003e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.340e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.224e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.057e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4479e+02 J
  → Fracture energy : 4.3321e-04 J
  → Total energy    : 2.4479e+02 J


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
  ||Δu||/||u|| = 5.378e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.646e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.095e-01

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
  ||Δu||/||u|| = 3.246e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.575e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4745e+02 J
  → Fracture energy : 4.3557e-04 J
  → Total energy    : 2.4745e+02 J


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
  ||Δu||/||u|| = 1.347e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.995e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.205e-03

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
  ||Δu||/||u|| = 9.244e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.921e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5069e+02 J
  → Fracture energy : 4.3559e-04 J
  → Total energy    : 2.5069e+02 J


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
  ||Δu||/||u|| = 8.424e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.523e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-10

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
  ||Δu||/||u|| = 2.764e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.883e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5402e+02 J
  → Fracture energy : 4.3559e-04 J
  → Total energy    : 2.5402e+02 J


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
  ||ΔD||/||D|| = 1.116e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.632e-12

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
  ||Δu||/||u|| = 1.292e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.247e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5673e+02 J
  → Fracture energy : 4.3559e-04 J
  → Total energy    : 2.5673e+02 J


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
  ||ΔD||/||D|| = 1.083e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.557e-12

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
  ||Δu||/||u|| = 2.698e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.839e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5945e+02 J
  → Fracture energy : 4.3559e-04 J
  → Total energy    : 2.5945e+02 J


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
  ||ΔD||/||D|| = 1.052e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.484e-12

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.874e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6219e+02 J
  → Fracture energy : 4.3559e-04 J
  → Total energy    : 2.6219e+02 J


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
  ||ΔD||/||D|| = 1.021e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.413e-12

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
  ||Δu||/||u|| = 1.128e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.612e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6494e+02 J
  → Fracture energy : 4.3559e-04 J
  → Total energy    : 2.6494e+02 J


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
  ||ΔD||/||D|| = 2.954e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.314e-03

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
  ||Δu||/||u|| = 8.923e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.835e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6771e+02 J
  → Fracture energy : 4.3559e-04 J
  → Total energy    : 2.6771e+02 J


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
  ||ΔD||/||D|| = 1.325e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.177e-02

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
  ||Δu||/||u|| = 1.527e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.552e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7049e+02 J
  → Fracture energy : 4.3560e-04 J
  → Total energy    : 2.7049e+02 J


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
  ||ΔD||/||D|| = 4.097e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.292e-01

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
  ||Δu||/||u|| = 2.463e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.462e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7329e+02 J
  → Fracture energy : 4.3577e-04 J
  → Total energy    : 2.7329e+02 J


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
  ||ΔD||/||D|| = 1.018e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.211e-01

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
  ||Δu||/||u|| = 1.739e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.896e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7610e+02 J
  → Fracture energy : 4.3720e-04 J
  → Total energy    : 2.7610e+02 J


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
  ||ΔD||/||D|| = 1.349e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.254e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.363e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7892e+02 J
  → Fracture energy : 4.4130e-04 J
  → Total energy    : 2.7892e+02 J


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
  ||Δu||/||u|| = 5.065e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.324e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.330e-02

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
  ||Δu||/||u|| = 4.333e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.887e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8176e+02 J
  → Fracture energy : 4.4226e-04 J
  → Total energy    : 2.8176e+02 J


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
  ||Δu||/||u|| = 1.400e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.351e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.230e-02

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
  ||Δu||/||u|| = 3.408e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.859e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8588e+02 J
  → Fracture energy : 4.4237e-04 J
  → Total energy    : 2.8588e+02 J


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
  ||Δu||/||u|| = 5.000e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.368e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.620e-01

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
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8876e+02 J
  → Fracture energy : 4.4293e-04 J
  → Total energy    : 2.8877e+02 J

Simulation completed in 45.39 s
Total time steps solved: 201
