[INFO] Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 17 entities
Info    : 11614 nodes
Info    : 23226 elements
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
[INFO]   area = 1.200e-08 m², perimeter = 5.200e-04 m
[INFO] === Mesh summary ===
[INFO]   Topology dim: 2
[INFO]   Facet dim: 1
[INFO]   Num cells: 21774
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
  ||Δu||/||u|| = 0.000e+00
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
  → Elastic energy  : 7.0523e-09 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.0523e-09 J


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
  → Elastic energy  : 2.8209e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8209e-08 J


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
  → Elastic energy  : 6.3470e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.3470e-08 J


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
  → Elastic energy  : 1.1284e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1284e-07 J


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
  → Elastic energy  : 1.7631e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.7631e-07 J


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
  → Elastic energy  : 2.5388e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5388e-07 J


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
  ||Δu||/||u|| = 3.854e-17
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
  → Elastic energy  : 3.4556e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.4556e-07 J


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
  → Elastic energy  : 4.5135e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.5135e-07 J


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
  → Elastic energy  : 5.7123e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.7123e-07 J


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
  → Elastic energy  : 7.0523e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.0523e-07 J


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
  → Elastic energy  : 8.5333e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.5333e-07 J


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
  → Elastic energy  : 1.0155e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0155e-06 J


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
  ||Δu||/||u|| = 5.868e-18
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
  → Elastic energy  : 1.1918e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1918e-06 J


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
  ||Δu||/||u|| = 6.556e-17
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
  → Elastic energy  : 1.3822e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3822e-06 J


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
  → Elastic energy  : 1.5868e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5868e-06 J


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
  → Elastic energy  : 1.8054e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.8054e-06 J


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
  → Elastic energy  : 2.0381e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0381e-06 J


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
  → Elastic energy  : 2.2849e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2849e-06 J


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
  → Elastic energy  : 2.5459e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5459e-06 J


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
  ||Δu||/||u|| = 1.902e-17
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
  → Elastic energy  : 2.8209e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8209e-06 J


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
  → Elastic energy  : 3.1101e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.1101e-06 J


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
  → Elastic energy  : 3.4133e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.4133e-06 J


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
  → Elastic energy  : 3.7307e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.7307e-06 J


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
  → Elastic energy  : 4.0621e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.0621e-06 J


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
  → Elastic energy  : 4.4077e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.4077e-06 J


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
  → Elastic energy  : 4.7673e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.7673e-06 J


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
  → Elastic energy  : 5.1411e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.1411e-06 J


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
  → Elastic energy  : 5.5290e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.5290e-06 J


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
  → Elastic energy  : 5.9310e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.9310e-06 J


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
  → Elastic energy  : 6.3470e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.3470e-06 J


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
  → Elastic energy  : 6.7772e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.7772e-06 J


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
  → Elastic energy  : 7.2215e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.2215e-06 J


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
  ||Δu||/||u|| = 3.422e-17
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
  → Elastic energy  : 7.6799e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.6799e-06 J


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
  → Elastic energy  : 8.1524e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.1524e-06 J


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
  ||Δu||/||u|| = 5.667e-17
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
  → Elastic energy  : 8.6390e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.6390e-06 J


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
  → Elastic energy  : 9.1397e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.1397e-06 J


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
  → Elastic energy  : 9.6546e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.6546e-06 J


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
  → Elastic energy  : 1.0183e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0183e-05 J


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
  ||Δu||/||u|| = 7.163e-17
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
  → Elastic energy  : 1.0727e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0727e-05 J


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
  ||Δu||/||u|| = 1.902e-17
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
  → Elastic energy  : 1.1284e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1284e-05 J


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
  ||Δu||/||u|| = 1.988e-18
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
  → Elastic energy  : 1.1855e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1855e-05 J


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
  → Elastic energy  : 1.2440e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2440e-05 J


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
  → Elastic energy  : 1.3040e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3040e-05 J


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
  → Elastic energy  : 1.3653e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3653e-05 J


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
  → Elastic energy  : 1.4281e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4281e-05 J


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
  → Elastic energy  : 1.4923e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4923e-05 J


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
  → Elastic energy  : 1.5578e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5578e-05 J


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
  ||Δu||/||u|| = 1.196e-16
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
  → Elastic energy  : 1.6248e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.6248e-05 J


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
  → Elastic energy  : 1.6933e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.6933e-05 J


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
  → Elastic energy  : 1.7631e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.7631e-05 J


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
  ||Δu||/||u|| = 5.254e-17
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
  → Elastic energy  : 1.8343e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.8343e-05 J


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
  → Elastic energy  : 1.9069e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.9069e-05 J


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
  ||Δu||/||u|| = 2.326e-17
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
  → Elastic energy  : 1.9810e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.9810e-05 J


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
  → Elastic energy  : 2.0564e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0564e-05 J


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
  ||Δu||/||u|| = 3.504e-17
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
  → Elastic energy  : 2.1333e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.1333e-05 J


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
  → Elastic energy  : 2.2116e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2116e-05 J


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
  → Elastic energy  : 2.2913e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2913e-05 J


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
  → Elastic energy  : 2.3724e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.3724e-05 J


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
  → Elastic energy  : 2.4549e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.4549e-05 J


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
  → Elastic energy  : 2.5388e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5388e-05 J


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
  → Elastic energy  : 2.6242e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.6242e-05 J


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
  → Elastic energy  : 2.7109e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.7109e-05 J


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
  ||Δu||/||u|| = 2.390e-17
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
  → Elastic energy  : 2.7990e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.7990e-05 J


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
  → Elastic energy  : 2.8886e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8886e-05 J


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
  → Elastic energy  : 2.9796e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.9796e-05 J


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
  → Elastic energy  : 3.0720e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.0720e-05 J


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
  → Elastic energy  : 3.1658e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.1658e-05 J


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
  → Elastic energy  : 3.2610e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.2610e-05 J


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
  → Elastic energy  : 3.3576e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.3576e-05 J


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
  → Elastic energy  : 3.4556e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.4556e-05 J


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
  → Elastic energy  : 3.5551e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.5551e-05 J


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
  → Elastic energy  : 3.6559e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.6559e-05 J


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
  → Elastic energy  : 3.7582e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.7582e-05 J


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
  → Elastic energy  : 3.8618e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.8618e-05 J


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
  ||Δu||/||u|| = 4.898e-17
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
  → Elastic energy  : 3.9669e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.9669e-05 J


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
  → Elastic energy  : 4.0734e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.0734e-05 J


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
  → Elastic energy  : 4.1813e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.1813e-05 J


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
  → Elastic energy  : 4.2906e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.2906e-05 J


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
  → Elastic energy  : 4.4013e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.4013e-05 J


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
  → Elastic energy  : 4.5135e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.5135e-05 J


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
  → Elastic energy  : 4.6270e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.6270e-05 J


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
  → Elastic energy  : 4.7419e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.7419e-05 J


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
  → Elastic energy  : 4.8583e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.8583e-05 J


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
  → Elastic energy  : 4.9761e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.9761e-05 J


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
  → Elastic energy  : 5.0953e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.0953e-05 J


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
  → Elastic energy  : 5.2159e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.2159e-05 J


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
  → Elastic energy  : 5.3379e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.3379e-05 J


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
  → Elastic energy  : 5.4613e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.4613e-05 J


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
  → Elastic energy  : 5.5861e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.5861e-05 J


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
  → Elastic energy  : 5.7123e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.7123e-05 J


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
  → Elastic energy  : 5.8400e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.8400e-05 J


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
  ||Δu||/||u|| = 3.049e-17
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
  → Elastic energy  : 5.9690e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.9690e-05 J


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
  → Elastic energy  : 6.0995e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.0995e-05 J


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
  → Elastic energy  : 6.2314e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.2314e-05 J


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
  → Elastic energy  : 6.3647e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.3647e-05 J


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
  → Elastic energy  : 6.4994e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.4994e-05 J


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
  → Elastic energy  : 6.6355e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.6355e-05 J


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
  → Elastic energy  : 6.7730e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.7730e-05 J


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
  → Elastic energy  : 6.9119e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.9119e-05 J


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
  → Elastic energy  : 7.0523e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.0523e-05 J


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
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

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
  → Elastic energy  : 7.1940e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.1940e-05 J


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
  ||Δu||/||u|| = 9.804e-03
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
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3372e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.3372e-05 J


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
  ||Δu||/||u|| = 9.709e-03
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
  **[INFO]** Updating traction on region 6 → -412000000.0 Pa
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
  → Elastic energy  : 7.4818e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.4818e-05 J


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
  ||Δu||/||u|| = 9.615e-03
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
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6277e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.6277e-05 J


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
  ||Δu||/||u|| = 9.524e-03
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
  **[INFO]** Updating traction on region 6 → -420000000.0 Pa
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
  → Elastic energy  : 7.7751e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.7751e-05 J


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
  ||Δu||/||u|| = 9.434e-03
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
  **[INFO]** Updating traction on region 6 → -424000000.0 Pa
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
  → Elastic energy  : 7.9239e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.9239e-05 J


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
  ||Δu||/||u|| = 9.346e-03
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
  **[INFO]** Updating traction on region 6 → -428000000.0 Pa
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
  → Elastic energy  : 8.0741e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.0741e-05 J


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
  ||Δu||/||u|| = 9.259e-03
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
  **[INFO]** Updating traction on region 6 → -432000000.0 Pa
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
  → Elastic energy  : 8.2258e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.2258e-05 J


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
  ||Δu||/||u|| = 9.174e-03
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
  **[INFO]** Updating traction on region 6 → -436000000.0 Pa
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
  → Elastic energy  : 8.3788e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.3788e-05 J


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
  ||Δu||/||u|| = 9.091e-03
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
  **[INFO]** Updating traction on region 6 → -440000000.0 Pa
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
  → Elastic energy  : 8.5333e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.5333e-05 J


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
  ||Δu||/||u|| = 9.009e-03
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
  **[INFO]** Updating traction on region 6 → -444000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.176e-16
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
  → Elastic energy  : 8.6891e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.6891e-05 J


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
  ||Δu||/||u|| = 8.929e-03
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
  **[INFO]** Updating traction on region 6 → -448000000.0 Pa
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
  → Elastic energy  : 8.8464e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.8464e-05 J


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
  ||Δu||/||u|| = 8.850e-03
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
  **[INFO]** Updating traction on region 6 → -452000000.0 Pa
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
  → Elastic energy  : 9.0050e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.0050e-05 J


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
  ||Δu||/||u|| = 8.772e-03
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
  **[INFO]** Updating traction on region 6 → -456000000.0 Pa
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
  → Elastic energy  : 9.1651e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.1651e-05 J


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
  ||Δu||/||u|| = 8.696e-03
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
  **[INFO]** Updating traction on region 6 → -460000000.0 Pa
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
  → Elastic energy  : 9.3266e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.3266e-05 J


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
  ||Δu||/||u|| = 8.621e-03
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
  **[INFO]** Updating traction on region 6 → -464000000.0 Pa
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
  → Elastic energy  : 9.4895e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.4895e-05 J


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
  ||Δu||/||u|| = 8.547e-03
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
  **[INFO]** Updating traction on region 6 → -468000000.0 Pa
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
  → Elastic energy  : 9.6539e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.6539e-05 J


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
  ||Δu||/||u|| = 8.475e-03
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
  **[INFO]** Updating traction on region 6 → -472000000.0 Pa
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
  → Elastic energy  : 9.8196e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.8196e-05 J


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
  ||Δu||/||u|| = 8.403e-03
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
  **[INFO]** Updating traction on region 6 → -476000000.0 Pa
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
  → Elastic energy  : 9.9867e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.9867e-05 J


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
  ||Δu||/||u|| = 8.333e-03
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
  **[INFO]** Updating traction on region 6 → -480000000.0 Pa
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
  → Elastic energy  : 1.0155e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0155e-04 J


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
  ||Δu||/||u|| = 8.264e-03
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
  **[INFO]** Updating traction on region 6 → -484000000.0 Pa
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
  → Elastic energy  : 1.0325e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0325e-04 J


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
  ||Δu||/||u|| = 8.197e-03
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
  **[INFO]** Updating traction on region 6 → -488000000.0 Pa
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
  → Elastic energy  : 1.0497e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0497e-04 J


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
  ||Δu||/||u|| = 8.130e-03
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
  **[INFO]** Updating traction on region 6 → -492000000.0 Pa
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
  → Elastic energy  : 1.0669e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0669e-04 J


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
  ||Δu||/||u|| = 8.065e-03
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
  **[INFO]** Updating traction on region 6 → -496000000.0 Pa
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
  → Elastic energy  : 1.0844e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0844e-04 J


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
  ||Δu||/||u|| = 8.000e-03
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
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1019e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1019e-04 J


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
  ||Δu||/||u|| = 7.937e-03
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
  **[INFO]** Updating traction on region 6 → -504000000.0 Pa
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
  → Elastic energy  : 1.1196e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1196e-04 J


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
  ||Δu||/||u|| = 7.874e-03
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
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1375e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1375e-04 J


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
  ||Δu||/||u|| = 7.812e-03
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
  **[INFO]** Updating traction on region 6 → -512000000.0 Pa
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
  → Elastic energy  : 1.1554e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1554e-04 J


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
  ||Δu||/||u|| = 7.752e-03
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
  **[INFO]** Updating traction on region 6 → -516000000.0 Pa
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
  → Elastic energy  : 1.1736e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1736e-04 J


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
  ||Δu||/||u|| = 7.692e-03
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
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1918e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1918e-04 J


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
  ||Δu||/||u|| = 7.634e-03
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
  **[INFO]** Updating traction on region 6 → -524000000.0 Pa
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
  → Elastic energy  : 1.2102e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2102e-04 J


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
  ||Δu||/||u|| = 7.576e-03
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
  **[INFO]** Updating traction on region 6 → -528000000.0 Pa
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
  → Elastic energy  : 1.2288e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2288e-04 J


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
  ||Δu||/||u|| = 7.519e-03
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
  **[INFO]** Updating traction on region 6 → -532000000.0 Pa
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
  → Elastic energy  : 1.2475e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2475e-04 J


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
  ||Δu||/||u|| = 7.463e-03
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
  **[INFO]** Updating traction on region 6 → -536000000.0 Pa
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
  → Elastic energy  : 1.2663e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2663e-04 J


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
  ||Δu||/||u|| = 7.407e-03
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
  **[INFO]** Updating traction on region 6 → -540000000.0 Pa
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
  → Elastic energy  : 1.2853e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2853e-04 J


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
  ||Δu||/||u|| = 7.353e-03
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
  **[INFO]** Updating traction on region 6 → -544000000.0 Pa
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
  → Elastic energy  : 1.3044e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3044e-04 J


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
  ||Δu||/||u|| = 7.299e-03
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
  **[INFO]** Updating traction on region 6 → -548000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.930e-17
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
  → Elastic energy  : 1.3236e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3236e-04 J


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
  ||Δu||/||u|| = 7.246e-03
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
  **[INFO]** Updating traction on region 6 → -552000000.0 Pa
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
  → Elastic energy  : 1.3430e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3430e-04 J


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
  ||Δu||/||u|| = 7.194e-03
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
  **[INFO]** Updating traction on region 6 → -556000000.0 Pa
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
  → Elastic energy  : 1.3626e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3626e-04 J


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
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

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
  → Elastic energy  : 1.3822e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3822e-04 J


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
  ||Δu||/||u|| = 7.092e-03
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
  **[INFO]** Updating traction on region 6 → -564000000.0 Pa
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
  → Elastic energy  : 1.4021e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4021e-04 J


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
  ||Δu||/||u|| = 7.042e-03
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
  **[INFO]** Updating traction on region 6 → -568000000.0 Pa
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
  → Elastic energy  : 1.4220e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4220e-04 J


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
  ||Δu||/||u|| = 6.993e-03
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
  **[INFO]** Updating traction on region 6 → -572000000.0 Pa
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
  → Elastic energy  : 1.4421e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4421e-04 J


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
  ||Δu||/||u|| = 6.944e-03
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
  **[INFO]** Updating traction on region 6 → -576000000.0 Pa
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
  → Elastic energy  : 1.4624e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4624e-04 J


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
  ||Δu||/||u|| = 6.897e-03
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
  **[INFO]** Updating traction on region 6 → -580000000.0 Pa
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
  → Elastic energy  : 1.4827e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4827e-04 J


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
  ||Δu||/||u|| = 6.849e-03
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
  **[INFO]** Updating traction on region 6 → -584000000.0 Pa
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
  → Elastic energy  : 1.5033e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5033e-04 J


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
  ||Δu||/||u|| = 6.803e-03
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
  **[INFO]** Updating traction on region 6 → -588000000.0 Pa
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
  → Elastic energy  : 1.5239e-04 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5239e-04 J


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
  ||Δu||/||u|| = 6.757e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 9.878e-04

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.500e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 9.384e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.928e-01
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 9.807e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.529e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 8.425e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.912e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 8.804e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.264e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 9.152e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
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
  |ΔD|_∞ = 9.458e-04

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.715e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 6.633e-04

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.017e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 6.931e-04

Convergence check


#### Iteration 10/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.294e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 7.205e-04

Convergence check


#### Iteration 11/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.538e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 7.446e-04

Convergence check


#### Iteration 12/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.740e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 7.646e-04

Convergence check


#### Iteration 13/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.891e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 7.795e-04

Convergence check


#### Iteration 14/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.079e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 5.017e-04

Convergence check


#### Iteration 15/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.361e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.272e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 5.208e-04

Convergence check


#### Iteration 16/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.361e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.439e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 5.373e-04

Convergence check


#### Iteration 17/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.575e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 5.507e-04

Convergence check


#### Iteration 18/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.673e-01
  [adaptive] relax_D=0.09
  |ΔD|_∞ = 5.604e-04

Convergence check


#### Iteration 19/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.725e-01
  [adaptive] relax_D=0.10
  |ΔD|_∞ = 5.655e-04

Convergence check


#### Iteration 20/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.726e-01
  [adaptive] relax_D=0.11
  |ΔD|_∞ = 5.656e-04

Convergence check


#### Iteration 21/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.669e-01
  [adaptive] relax_D=0.12
  |ΔD|_∞ = 5.600e-04

Convergence check


#### Iteration 22/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.551e-01
  [adaptive] relax_D=0.13
  |ΔD|_∞ = 5.484e-04

Convergence check


#### Iteration 23/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.368e-01
  [adaptive] relax_D=0.15
  |ΔD|_∞ = 5.303e-04

Convergence check


#### Iteration 24/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.120e-01
  [adaptive] relax_D=0.16
  |ΔD|_∞ = 5.058e-04

Convergence check


#### Iteration 25/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.809e-01
  [adaptive] relax_D=0.18
  |ΔD|_∞ = 4.750e-04

Convergence check


#### Iteration 26/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.439e-01
  [adaptive] relax_D=0.19
  |ΔD|_∞ = 4.385e-04

Convergence check


#### Iteration 27/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.019e-01
  [adaptive] relax_D=0.21
  |ΔD|_∞ = 3.970e-04

Convergence check


#### Iteration 28/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
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
  |ΔD|_∞ = 3.517e-04

Convergence check


#### Iteration 29/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.078e-01
  [adaptive] relax_D=0.26
  |ΔD|_∞ = 3.041e-04

Convergence check


#### Iteration 30/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
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
  |ΔD|_∞ = 2.557e-04

Convergence check


#### Iteration 31/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.110e-01
  [adaptive] relax_D=0.31
  |ΔD|_∞ = 2.085e-04

Convergence check


#### Iteration 32/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.660e-01
  [adaptive] relax_D=0.34
  |ΔD|_∞ = 1.640e-04

Convergence check


#### Iteration 33/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.254e-01
  [adaptive] relax_D=0.38
  |ΔD|_∞ = 1.238e-04

Convergence check


#### Iteration 34/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
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
  |ΔD|_∞ = 8.926e-05

Convergence check


#### Iteration 35/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.170e-02
  [adaptive] relax_D=0.46
  |ΔD|_∞ = 6.095e-05

Convergence check


#### Iteration 36/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.956e-02
  [adaptive] relax_D=0.50
  |ΔD|_∞ = 3.908e-05

Convergence check


#### Iteration 37/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.355e-02
  [adaptive] relax_D=0.56
  |ΔD|_∞ = 2.326e-05

Convergence check


#### Iteration 38/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.283e-02
  [adaptive] relax_D=0.61
  |ΔD|_∞ = 1.267e-05

Convergence check


#### Iteration 39/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
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
  |ΔD|_∞ = 6.201e-06

Convergence check


#### Iteration 40/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.688e-03
  [adaptive] relax_D=0.74
  |ΔD|_∞ = 2.655e-06

Convergence check


#### Iteration 41/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
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
  |ΔD|_∞ = 9.586e-07

Convergence check


#### Iteration 42/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
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
  |ΔD|_∞ = 2.752e-07

Convergence check


#### Iteration 43/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -592000000.0 Pa
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
  |ΔD|_∞ = 5.665e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 43 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5446e-04 J
  → Fracture energy : 2.3449e-09 J
  → Total energy    : 1.5446e-04 J


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
  ||Δu||/||u|| = 6.813e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.981e-01
  [adaptive] relax_D=0.98
  |ΔD|_∞ = 3.911e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.310e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.420e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -596000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.186e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.071e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating traction on region 6 → -596000000.0 Pa
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
  → Elastic energy  : 1.5651e-04 J
  → Fracture energy : 2.4235e-08 J
  → Total energy    : 1.5654e-04 J


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
  ||Δu||/||u|| = 7.140e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.350e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.050e-02

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
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5859e-04 J
  → Fracture energy : 7.4147e-08 J
  → Total energy    : 1.5867e-04 J


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
  ||Δu||/||u|| = 7.819e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.428e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.462e-02

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
  → Elastic energy  : 1.6073e-04 J
  → Fracture energy : 2.0675e-07 J
  → Total energy    : 1.6094e-04 J


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
  ||Δu||/||u|| = 9.005e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.154e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.298e-02

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
  → Elastic energy  : 1.6295e-04 J
  → Fracture energy : 4.6692e-07 J
  → Total energy    : 1.6341e-04 J


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
  ||Δu||/||u|| = 1.112e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.910e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.116e-01

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
  → Elastic energy  : 1.6527e-04 J
  → Fracture energy : 9.2454e-07 J
  → Total energy    : 1.6619e-04 J


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
  ||Δu||/||u|| = 1.497e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.783e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.515e-01

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
  → Elastic energy  : 1.6777e-04 J
  → Fracture energy : 1.7949e-06 J
  → Total energy    : 1.6957e-04 J


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
  ||Δu||/||u|| = 2.180e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.711e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.129e-01

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
  ||Δu||/||u|| = 1.172e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.570e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7058e-04 J
  → Fracture energy : 3.3242e-06 J
  → Total energy    : 1.7390e-04 J


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
  ||Δu||/||u|| = 3.410e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.634e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.960e-01

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
  → Elastic energy  : 1.7409e-04 J
  → Fracture energy : 5.9409e-06 J
  → Total energy    : 1.8003e-04 J


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
  ||Δu||/||u|| = 5.883e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.445e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.764e-01

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
  → Elastic energy  : 1.7916e-04 J
  → Fracture energy : 1.0189e-05 J
  → Total energy    : 1.8935e-04 J


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
  ||Δu||/||u|| = 7.046e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.096e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.397e-01

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
  ||Δu||/||u|| = 2.927e-17
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
  → Elastic energy  : 1.8615e-04 J
  → Fracture energy : 1.7091e-05 J
  → Total energy    : 2.0324e-04 J


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
  ||Δu||/||u|| = 9.993e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.682e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.661e-01

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
  ||ΔD||/||D|| = 1.451e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3267e-01 J
  → Fracture energy : 2.7922e-05 J
  → Total energy    : 1.3270e-01 J


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
  ||Δu||/||u|| = 9.282e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.287e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.692e-01

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
  ||Δu||/||u|| = 1.319e-21
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.555e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0049e+00 J
  → Fracture energy : 4.4157e-05 J
  → Total energy    : 1.0050e+00 J


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
  ||Δu||/||u|| = 8.152e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.068e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.559e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.020e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7147e+00 J
  → Fracture energy : 6.6201e-05 J
  → Total energy    : 3.7148e+00 J


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
  ||Δu||/||u|| = 6.978e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.474e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.757e-01

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
  ||Δu||/||u|| = 1.100e-21
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.619e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4057e+00 J
  → Fracture energy : 9.2442e-05 J
  → Total energy    : 9.4058e+00 J


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
  ||Δu||/||u|| = 5.438e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.326e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.612e-01

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
  ||ΔD||/||D|| = 1.090e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8002e+01 J
  → Fracture energy : 1.1905e-04 J
  → Total energy    : 1.8002e+01 J


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
  ||Δu||/||u|| = 4.941e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.118e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.940e-01

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
  ||ΔD||/||D|| = 1.416e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3239e+01 J
  → Fracture energy : 1.4866e-04 J
  → Total energy    : 3.3239e+01 J


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
  ||Δu||/||u|| = 3.438e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.771e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.897e-01

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
  → Elastic energy  : 4.6319e+01 J
  → Fracture energy : 1.8048e-04 J
  → Total energy    : 4.6319e+01 J


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
  ||Δu||/||u|| = 2.842e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.799e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.965e-01

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
  ||ΔD||/||D|| = 1.334e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1068e+01 J
  → Fracture energy : 2.1000e-04 J
  → Total energy    : 6.1068e+01 J


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
  ||Δu||/||u|| = 1.862e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.431e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.451e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.134e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1772e+01 J
  → Fracture energy : 2.3894e-04 J
  → Total energy    : 7.1773e+01 J


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
  ||Δu||/||u|| = 1.340e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.632e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.214e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.321e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1269e+01 J
  → Fracture energy : 2.7002e-04 J
  → Total energy    : 8.1269e+01 J


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
  ||Δu||/||u|| = 1.111e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.984e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.537e-01

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
  ||ΔD||/||D|| = 8.612e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.9388e+01 J
  → Fracture energy : 3.0436e-04 J
  → Total energy    : 8.9388e+01 J


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
  ||Δu||/||u|| = 1.123e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.866e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.935e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.210e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8574e+01 J
  → Fracture energy : 3.2592e-04 J
  → Total energy    : 9.8574e+01 J


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
  ||Δu||/||u|| = 1.160e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.928e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.498e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.650e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1047e+02 J
  → Fracture energy : 3.4779e-04 J
  → Total energy    : 1.1047e+02 J


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
  ||Δu||/||u|| = 2.500e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.915e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.051e-01

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
  ||Δu||/||u|| = 7.487e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.679e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1318e+02 J
  → Fracture energy : 3.6179e-04 J
  → Total energy    : 1.1318e+02 J


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
  ||Δu||/||u|| = 7.532e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.491e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.426e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.581e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2106e+02 J
  → Fracture energy : 3.8060e-04 J
  → Total energy    : 1.2106e+02 J


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
  ||Δu||/||u|| = 1.095e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.991e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.890e-01

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
  ||ΔD||/||D|| = 1.593e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2290e+02 J
  → Fracture energy : 4.0280e-04 J
  → Total energy    : 1.2290e+02 J


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
  ||Δu||/||u|| = 5.253e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.280e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.724e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.738e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2857e+02 J
  → Fracture energy : 4.1808e-04 J
  → Total energy    : 1.2857e+02 J


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
  ||Δu||/||u|| = 5.693e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.007e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.209e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.055e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3480e+02 J
  → Fracture energy : 4.3250e-04 J
  → Total energy    : 1.3480e+02 J


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
  ||Δu||/||u|| = 1.727e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.440e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.692e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.359e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3747e+02 J
  → Fracture energy : 4.4134e-04 J
  → Total energy    : 1.3747e+02 J


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
  ||Δu||/||u|| = 3.766e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.396e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.972e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.398e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4244e+02 J
  → Fracture energy : 4.4805e-04 J
  → Total energy    : 1.4245e+02 J


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
  ||Δu||/||u|| = 2.067e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.243e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.672e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.329e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4538e+02 J
  → Fracture energy : 4.5355e-04 J
  → Total energy    : 1.4538e+02 J


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
  ||Δu||/||u|| = 2.809e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.331e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.605e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.002e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4882e+02 J
  → Fracture energy : 4.5963e-04 J
  → Total energy    : 1.4882e+02 J


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
  ||Δu||/||u|| = 1.066e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.545e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.221e-01

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
  ||ΔD||/||D|| = 1.465e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5069e+02 J
  → Fracture energy : 4.6387e-04 J
  → Total energy    : 1.5069e+02 J


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
  ||Δu||/||u|| = 1.337e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.915e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.063e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.365e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5337e+02 J
  → Fracture energy : 4.6596e-04 J
  → Total energy    : 1.5337e+02 J


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
  ||Δu||/||u|| = 2.329e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.217e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.251e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.501e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5662e+02 J
  → Fracture energy : 4.7079e-04 J
  → Total energy    : 1.5662e+02 J


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
  ||Δu||/||u|| = 5.442e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.058e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.374e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.406e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5834e+02 J
  → Fracture energy : 4.7321e-04 J
  → Total energy    : 1.5834e+02 J


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
  ||Δu||/||u|| = 2.257e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.093e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.035e-01

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
  ||ΔD||/||D|| = 6.790e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6152e+02 J
  → Fracture energy : 4.7549e-04 J
  → Total energy    : 1.6152e+02 J


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
  ||Δu||/||u|| = 5.420e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.053e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.335e-01

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
  ||Δu||/||u|| = 5.884e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.629e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6328e+02 J
  → Fracture energy : 4.7732e-04 J
  → Total energy    : 1.6328e+02 J


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
  ||Δu||/||u|| = 1.225e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.983e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.990e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.076e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6566e+02 J
  → Fracture energy : 4.7903e-04 J
  → Total energy    : 1.6566e+02 J


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
  ||Δu||/||u|| = 5.876e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.065e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.018e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.710e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6744e+02 J
  → Fracture energy : 4.7904e-04 J
  → Total energy    : 1.6744e+02 J


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
  ||Δu||/||u|| = 1.192e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.554e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.729e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.628e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6990e+02 J
  → Fracture energy : 4.7903e-04 J
  → Total energy    : 1.6990e+02 J


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
  ||Δu||/||u|| = 5.269e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.584e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.399e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.621e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7170e+02 J
  → Fracture energy : 4.7896e-04 J
  → Total energy    : 1.7170e+02 J


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
  ||Δu||/||u|| = 5.240e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.132e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.829e-02

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
  ||ΔD||/||D|| = 2.389e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7351e+02 J
  → Fracture energy : 4.7887e-04 J
  → Total energy    : 1.7351e+02 J


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
  ||Δu||/||u|| = 5.213e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.332e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.334e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.452e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7533e+02 J
  → Fracture energy : 4.7890e-04 J
  → Total energy    : 1.7533e+02 J


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
  ||Δu||/||u|| = 5.184e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.797e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.930e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.192e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7717e+02 J
  → Fracture energy : 4.7913e-04 J
  → Total energy    : 1.7717e+02 J


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
  ||Δu||/||u|| = 5.159e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.618e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.469e-02

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.683e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7901e+02 J
  → Fracture energy : 4.7967e-04 J
  → Total energy    : 1.7901e+02 J


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
  ||Δu||/||u|| = 5.100e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.350e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.234e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.090e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8086e+02 J
  → Fracture energy : 4.8100e-04 J
  → Total energy    : 1.8086e+02 J


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
  ||Δu||/||u|| = 5.106e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.185e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.282e-01

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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.416e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8272e+02 J
  → Fracture energy : 4.8470e-04 J
  → Total energy    : 1.8272e+02 J


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
  ||Δu||/||u|| = 5.075e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.664e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.623e-01

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
  ||ΔD||/||D|| = 1.626e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8459e+02 J
  → Fracture energy : 4.9451e-04 J
  → Total energy    : 1.8459e+02 J


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
  ||Δu||/||u|| = 5.056e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.868e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.036e-01

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
  ||Δu||/||u|| = 1.247e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.161e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8647e+02 J
  → Fracture energy : 5.0850e-04 J
  → Total energy    : 1.8647e+02 J


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
  ||Δu||/||u|| = 6.947e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.047e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.717e-01

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
  ||Δu||/||u|| = 1.001e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.907e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8843e+02 J
  → Fracture energy : 5.2252e-04 J
  → Total energy    : 1.8843e+02 J


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
  ||Δu||/||u|| = 1.886e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.915e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.351e-01

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
  → Elastic energy  : 1.9245e+02 J
  → Fracture energy : 5.3480e-04 J
  → Total energy    : 1.9245e+02 J

Simulation completed in 190.58 s
Total time steps solved: 201
