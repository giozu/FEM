[INFO] Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 13 entities
Info    : 3635 nodes
Info    : 7268 elements
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
[INFO]   area = 3.600e-09 m², perimeter = 2.400e-04 m
[INFO] === Mesh summary ===
[INFO]   Topology dim: 2
[INFO]   Facet dim: 1
[INFO]   Num cells: 6790
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
  **[INFO]** Updating traction on region 6 → -6000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -6000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -6000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -6000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -6000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -6000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -6000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.977e-16
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
  → Elastic energy  : 8.0009e-09 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.0009e-09 J


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
  **[INFO]** Updating traction on region 6 → -12000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -12000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.924e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2004e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.2004e-08 J


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
  **[INFO]** Updating traction on region 6 → -18000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -18000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.586e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.2008e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.2008e-08 J


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
  **[INFO]** Updating traction on region 6 → -24000000.0 Pa
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
  → Elastic energy  : 1.2801e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2801e-07 J


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
  **[INFO]** Updating traction on region 6 → -30000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -30000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.860e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0002e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0002e-07 J


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
  **[INFO]** Updating traction on region 6 → -36000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -36000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.162e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8803e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8803e-07 J


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
  **[INFO]** Updating traction on region 6 → -42000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -42000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.368e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9204e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.9204e-07 J


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
  **[INFO]** Updating traction on region 6 → -48000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -48000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.454e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1206e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.1206e-07 J


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
  **[INFO]** Updating traction on region 6 → -54000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -54000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.066e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4807e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.4807e-07 J


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
  **[INFO]** Updating traction on region 6 → -60000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -60000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.410e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0009e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.0009e-07 J


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
  **[INFO]** Updating traction on region 6 → -66000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -66000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.690e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6811e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.6811e-07 J


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
  **[INFO]** Updating traction on region 6 → -72000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -72000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.586e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1521e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1521e-06 J


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
  **[INFO]** Updating traction on region 6 → -78000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -78000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.039e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3521e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3521e-06 J


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
  **[INFO]** Updating traction on region 6 → -84000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -84000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.970e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5682e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5682e-06 J


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
  **[INFO]** Updating traction on region 6 → -90000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -90000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.151e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8002e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.8002e-06 J


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
  **[INFO]** Updating traction on region 6 → -96000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -96000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.362e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0482e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0482e-06 J


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
  **[INFO]** Updating traction on region 6 → -102000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -102000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.143e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3123e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.3123e-06 J


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
  **[INFO]** Updating traction on region 6 → -108000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -108000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.061e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5923e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5923e-06 J


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
  **[INFO]** Updating traction on region 6 → -114000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -114000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.407e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8883e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8883e-06 J


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
  **[INFO]** Updating traction on region 6 → -120000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -120000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.268e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2004e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.2004e-06 J


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
  **[INFO]** Updating traction on region 6 → -126000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -126000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.547e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5284e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.5284e-06 J


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
  **[INFO]** Updating traction on region 6 → -132000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -132000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.154e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8724e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.8724e-06 J


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
  **[INFO]** Updating traction on region 6 → -138000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -138000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.724e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2325e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.2325e-06 J


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
  **[INFO]** Updating traction on region 6 → -144000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -144000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.232e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6085e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.6085e-06 J


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
  **[INFO]** Updating traction on region 6 → -150000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -150000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.113e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0006e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.0006e-06 J


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
  **[INFO]** Updating traction on region 6 → -156000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -156000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.226e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4086e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.4086e-06 J


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
  **[INFO]** Updating traction on region 6 → -162000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -162000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.460e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.8326e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.8326e-06 J


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
  **[INFO]** Updating traction on region 6 → -168000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -168000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.998e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2727e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.2727e-06 J


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
  **[INFO]** Updating traction on region 6 → -174000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -174000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.637e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.7287e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.7287e-06 J


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
  **[INFO]** Updating traction on region 6 → -180000000.0 Pa
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
  → Elastic energy  : 7.2008e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.2008e-06 J


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
  **[INFO]** Updating traction on region 6 → -186000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -186000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.553e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6889e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.6889e-06 J


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
  **[INFO]** Updating traction on region 6 → -192000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -192000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.401e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1929e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.1929e-06 J


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
  **[INFO]** Updating traction on region 6 → -198000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -198000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.939e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7130e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.7130e-06 J


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
  **[INFO]** Updating traction on region 6 → -204000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -204000000.0 Pa
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
  → Elastic energy  : 9.2490e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.2490e-06 J


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
  **[INFO]** Updating traction on region 6 → -210000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -210000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.820e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8011e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.8011e-06 J


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
  **[INFO]** Updating traction on region 6 → -216000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -216000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.035e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0369e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0369e-05 J


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
  **[INFO]** Updating traction on region 6 → -222000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -222000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.660e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0953e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0953e-05 J


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
  **[INFO]** Updating traction on region 6 → -228000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -228000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.462e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1553e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1553e-05 J


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
  **[INFO]** Updating traction on region 6 → -234000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -234000000.0 Pa
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
  → Elastic energy  : 1.2169e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2169e-05 J


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
  **[INFO]** Updating traction on region 6 → -240000000.0 Pa
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
  → Elastic energy  : 1.2801e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2801e-05 J


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
  **[INFO]** Updating traction on region 6 → -246000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -246000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.828e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3449e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3449e-05 J


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
  **[INFO]** Updating traction on region 6 → -252000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -252000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.553e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4114e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4114e-05 J


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
  **[INFO]** Updating traction on region 6 → -258000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -258000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.430e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4794e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4794e-05 J


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
  **[INFO]** Updating traction on region 6 → -264000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -264000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.916e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5490e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5490e-05 J


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
  **[INFO]** Updating traction on region 6 → -270000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -270000000.0 Pa
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
  → Elastic energy  : 1.6202e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.6202e-05 J


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
  **[INFO]** Updating traction on region 6 → -276000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -276000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.429e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6930e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.6930e-05 J


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
  **[INFO]** Updating traction on region 6 → -282000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -282000000.0 Pa
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
  → Elastic energy  : 1.7674e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.7674e-05 J


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
  **[INFO]** Updating traction on region 6 → -288000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -288000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.162e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8434e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.8434e-05 J


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
  **[INFO]** Updating traction on region 6 → -294000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -294000000.0 Pa
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
  → Elastic energy  : 1.9210e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.9210e-05 J


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
  **[INFO]** Updating traction on region 6 → -300000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -300000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.174e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0002e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0002e-05 J


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
  **[INFO]** Updating traction on region 6 → -306000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -306000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.265e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0810e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0810e-05 J


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
  **[INFO]** Updating traction on region 6 → -312000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -312000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.944e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1634e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.1634e-05 J


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
  **[INFO]** Updating traction on region 6 → -318000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -318000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.490e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2474e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2474e-05 J


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
  **[INFO]** Updating traction on region 6 → -324000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -324000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.842e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3331e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.3331e-05 J


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
  **[INFO]** Updating traction on region 6 → -330000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -330000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.952e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4203e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.4203e-05 J


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
  **[INFO]** Updating traction on region 6 → -336000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -336000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.920e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5091e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5091e-05 J


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
  **[INFO]** Updating traction on region 6 → -342000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -342000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.836e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5995e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5995e-05 J


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
  **[INFO]** Updating traction on region 6 → -348000000.0 Pa
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
  → Elastic energy  : 2.6915e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.6915e-05 J


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
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
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
  → Elastic energy  : 2.7851e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.7851e-05 J


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
  **[INFO]** Updating traction on region 6 → -360000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -360000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.054e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8803e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8803e-05 J


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
  **[INFO]** Updating traction on region 6 → -366000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -366000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.521e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9771e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.9771e-05 J


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
  **[INFO]** Updating traction on region 6 → -372000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -372000000.0 Pa
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
  → Elastic energy  : 3.0755e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.0755e-05 J


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
  **[INFO]** Updating traction on region 6 → -378000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -378000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.545e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1756e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.1756e-05 J


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
  **[INFO]** Updating traction on region 6 → -384000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -384000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.508e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2772e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.2772e-05 J


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
  **[INFO]** Updating traction on region 6 → -390000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -390000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.197e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3804e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.3804e-05 J


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
  **[INFO]** Updating traction on region 6 → -396000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -396000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.077e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4852e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.4852e-05 J


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
  **[INFO]** Updating traction on region 6 → -402000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -402000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.918e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5916e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.5916e-05 J


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
  **[INFO]** Updating traction on region 6 → -408000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -408000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.424e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6996e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.6996e-05 J


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
  **[INFO]** Updating traction on region 6 → -414000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -414000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.400e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8092e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.8092e-05 J


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
  **[INFO]** Updating traction on region 6 → -420000000.0 Pa
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
  → Elastic energy  : 3.9204e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.9204e-05 J


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
  **[INFO]** Updating traction on region 6 → -426000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -426000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.980e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0332e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.0332e-05 J


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
  **[INFO]** Updating traction on region 6 → -432000000.0 Pa
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
  → Elastic energy  : 4.1477e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.1477e-05 J


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
  **[INFO]** Updating traction on region 6 → -438000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -438000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.778e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2637e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.2637e-05 J


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
  **[INFO]** Updating traction on region 6 → -444000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -444000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.660e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3813e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.3813e-05 J


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
  **[INFO]** Updating traction on region 6 → -450000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -450000000.0 Pa
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
  → Elastic energy  : 4.5005e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.5005e-05 J


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
  **[INFO]** Updating traction on region 6 → -456000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -456000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.169e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6213e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.6213e-05 J


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
  **[INFO]** Updating traction on region 6 → -462000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -462000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.461e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7437e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.7437e-05 J


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
  **[INFO]** Updating traction on region 6 → -468000000.0 Pa
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
  → Elastic energy  : 4.8677e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.8677e-05 J


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
  **[INFO]** Updating traction on region 6 → -474000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -474000000.0 Pa
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
  → Elastic energy  : 4.9934e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.9934e-05 J


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
  **[INFO]** Updating traction on region 6 → -480000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -480000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.978e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1206e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.1206e-05 J


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
  **[INFO]** Updating traction on region 6 → -486000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -486000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.794e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2494e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.2494e-05 J


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
  **[INFO]** Updating traction on region 6 → -492000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -492000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.523e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3798e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.3798e-05 J


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
  **[INFO]** Updating traction on region 6 → -498000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -498000000.0 Pa
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
  → Elastic energy  : 5.5118e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.5118e-05 J


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
  **[INFO]** Updating traction on region 6 → -504000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -504000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.271e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.6454e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.6454e-05 J


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
  **[INFO]** Updating traction on region 6 → -510000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -510000000.0 Pa
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
  → Elastic energy  : 5.7806e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.7806e-05 J


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
  **[INFO]** Updating traction on region 6 → -516000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -516000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.931e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9175e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.9175e-05 J


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
  **[INFO]** Updating traction on region 6 → -522000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -522000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.141e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.0559e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.0559e-05 J


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
  **[INFO]** Updating traction on region 6 → -528000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -528000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.780e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1959e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.1959e-05 J


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
  **[INFO]** Updating traction on region 6 → -534000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -534000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.618e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.3375e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.3375e-05 J


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
  **[INFO]** Updating traction on region 6 → -540000000.0 Pa
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
  **[INFO]** Updating traction on region 6 → -540000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.337e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4807e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.4807e-05 J


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
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.099e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 2.193e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 2.084e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.927e-01
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 2.177e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.870e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.955e-03

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 2.032e-03

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.449e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.574e-01
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 2.100e-03

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.449e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.715e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.473e-03

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.539e-03

Convergence check


#### Iteration 10/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.600e-03

Convergence check


#### Iteration 11/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.653e-03

Convergence check


#### Iteration 12/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.698e-03

Convergence check


#### Iteration 13/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.731e-03

Convergence check


#### Iteration 14/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.114e-03

Convergence check


#### Iteration 15/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.272e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 1.156e-03

Convergence check


#### Iteration 16/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.439e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 1.193e-03

Convergence check


#### Iteration 17/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.223e-03

Convergence check


#### Iteration 18/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.244e-03

Convergence check


#### Iteration 19/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.256e-03

Convergence check


#### Iteration 20/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.256e-03

Convergence check


#### Iteration 21/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.361e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.669e-01
  [adaptive] relax_D=0.12
  |ΔD|_∞ = 1.243e-03

Convergence check


#### Iteration 22/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.361e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.551e-01
  [adaptive] relax_D=0.13
  |ΔD|_∞ = 1.217e-03

Convergence check


#### Iteration 23/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.177e-03

Convergence check


#### Iteration 24/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.123e-03

Convergence check


#### Iteration 25/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.055e-03

Convergence check


#### Iteration 26/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 9.736e-04

Convergence check


#### Iteration 27/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 8.815e-04

Convergence check


#### Iteration 28/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 7.809e-04

Convergence check


#### Iteration 29/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.895e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.078e-01
  [adaptive] relax_D=0.26
  |ΔD|_∞ = 6.752e-04

Convergence check


#### Iteration 30/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.895e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.589e-01
  [adaptive] relax_D=0.28
  |ΔD|_∞ = 5.678e-04

Convergence check


#### Iteration 31/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 4.628e-04

Convergence check


#### Iteration 32/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.462e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.660e-01
  [adaptive] relax_D=0.34
  |ΔD|_∞ = 3.640e-04

Convergence check


#### Iteration 33/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.462e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.254e-01
  [adaptive] relax_D=0.38
  |ΔD|_∞ = 2.749e-04

Convergence check


#### Iteration 34/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.449e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.036e-02
  [adaptive] relax_D=0.42
  |ΔD|_∞ = 1.982e-04

Convergence check


#### Iteration 35/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.449e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.170e-02
  [adaptive] relax_D=0.46
  |ΔD|_∞ = 1.353e-04

Convergence check


#### Iteration 36/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.436e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.956e-02
  [adaptive] relax_D=0.50
  |ΔD|_∞ = 8.677e-05

Convergence check


#### Iteration 37/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.124e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.355e-02
  [adaptive] relax_D=0.56
  |ΔD|_∞ = 5.165e-05

Convergence check


#### Iteration 38/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.097e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.283e-02
  [adaptive] relax_D=0.61
  |ΔD|_∞ = 2.814e-05

Convergence check


#### Iteration 39/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.909e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.277e-03
  [adaptive] relax_D=0.67
  |ΔD|_∞ = 1.377e-05

Convergence check


#### Iteration 40/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 5.895e-06

Convergence check


#### Iteration 41/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 2.128e-06

Convergence check


#### Iteration 42/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 6.111e-07

Convergence check


#### Iteration 43/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -546000000.0 Pa
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
  |ΔD|_∞ = 1.258e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 43 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6227e-05 J
  → Fracture energy : 7.8771e-09 J
  → Total energy    : 6.6235e-05 J


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
  **[INFO]** Updating traction on region 6 → -552000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.132e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.581e-01
  [adaptive] relax_D=0.98
  |ΔD|_∞ = 6.392e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -552000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.694e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.244e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.049e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -552000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.694e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.077e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.751e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -552000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.959e-17
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
  → Elastic energy  : 6.7636e-05 J
  → Fracture energy : 4.1055e-08 J
  → Total energy    : 6.7677e-05 J


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
  **[INFO]** Updating traction on region 6 → -558000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.264e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.364e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.900e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -558000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.699e-17
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
  → Elastic energy  : 6.9099e-05 J
  → Fracture energy : 1.4306e-07 J
  → Total energy    : 6.9242e-05 J


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
  **[INFO]** Updating traction on region 6 → -564000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.498e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.777e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.192e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -564000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.312e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.189e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.140e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.0610e-05 J
  → Fracture energy : 3.6160e-07 J
  → Total energy    : 7.0972e-05 J


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
  **[INFO]** Updating traction on region 6 → -570000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.963e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.468e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.654e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -570000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.546e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.559e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.360e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.2206e-05 J
  → Fracture energy : 7.7989e-07 J
  → Total energy    : 7.2986e-05 J


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
  **[INFO]** Updating traction on region 6 → -576000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.842e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.219e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.254e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -576000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.129e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.375e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.608e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3990e-05 J
  → Fracture energy : 1.5749e-06 J
  → Total energy    : 7.5565e-05 J


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
  **[INFO]** Updating traction on region 6 → -582000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.433e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.900e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.161e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -582000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.675e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.793e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.955e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6244e-05 J
  → Fracture energy : 2.9369e-06 J
  → Total energy    : 7.9181e-05 J


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
  **[INFO]** Updating traction on region 6 → -588000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.578e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.473e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.151e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
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
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9642e-05 J
  → Fracture energy : 5.0426e-06 J
  → Total energy    : 8.4685e-05 J


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
  **[INFO]** Updating traction on region 6 → -594000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.335e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.957e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.742e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -594000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.444e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.766e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.482e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5179e-05 J
  → Fracture energy : 8.2207e-06 J
  → Total energy    : 9.3399e-05 J


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
  **[INFO]** Updating traction on region 6 → -600000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.990e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.513e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.848e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -600000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.829e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0602e-02 J
  → Fracture energy : 1.3138e-05 J
  → Total energy    : 9.0615e-02 J


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
  **[INFO]** Updating traction on region 6 → -606000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.096e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.198e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.851e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -606000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.194e-20
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.033e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5908e-01 J
  → Fracture energy : 2.0218e-05 J
  → Total energy    : 5.5910e-01 J


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
  **[INFO]** Updating traction on region 6 → -612000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.786e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.989e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.792e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -612000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.298e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.866e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7779e+00 J
  → Fracture energy : 3.0016e-05 J
  → Total energy    : 1.7780e+00 J


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
  **[INFO]** Updating traction on region 6 → -618000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.208e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.894e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.780e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -618000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.883e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.096e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8466e+00 J
  → Fracture energy : 4.2713e-05 J
  → Total energy    : 3.8467e+00 J


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
  **[INFO]** Updating traction on region 6 → -624000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.983e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.528e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.176e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -624000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.282e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.509e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8699e+00 J
  → Fracture energy : 5.7308e-05 J
  → Total energy    : 6.8699e+00 J


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
  **[INFO]** Updating traction on region 6 → -630000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.248e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.414e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.880e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -630000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.161e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.445e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1209e+01 J
  → Fracture energy : 7.2423e-05 J
  → Total energy    : 1.1210e+01 J


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
  **[INFO]** Updating traction on region 6 → -636000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.760e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.129e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.183e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -636000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.344e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.089e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7393e+01 J
  → Fracture energy : 8.8700e-05 J
  → Total energy    : 1.7394e+01 J


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
  **[INFO]** Updating traction on region 6 → -642000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.039e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.040e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.752e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -642000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.368e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3823e+01 J
  → Fracture energy : 1.0590e-04 J
  → Total energy    : 2.3823e+01 J


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
  **[INFO]** Updating traction on region 6 → -648000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.633e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.872e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.846e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -648000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.964e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.826e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0099e+01 J
  → Fracture energy : 1.2188e-04 J
  → Total energy    : 3.0099e+01 J


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
  **[INFO]** Updating traction on region 6 → -654000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.096e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.959e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.607e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -654000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.108e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.760e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5978e+01 J
  → Fracture energy : 1.3563e-04 J
  → Total energy    : 3.5978e+01 J


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
  **[INFO]** Updating traction on region 6 → -660000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.345e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.859e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.138e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -660000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.681e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1842e+01 J
  → Fracture energy : 1.5201e-04 J
  → Total energy    : 4.1842e+01 J


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
  **[INFO]** Updating traction on region 6 → -666000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.165e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.786e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.373e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -666000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.247e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.807e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5365e+01 J
  → Fracture energy : 1.6785e-04 J
  → Total energy    : 4.5365e+01 J


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
  **[INFO]** Updating traction on region 6 → -672000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.043e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.741e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.763e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -672000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.402e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9767e+01 J
  → Fracture energy : 1.7733e-04 J
  → Total energy    : 4.9767e+01 J


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
  **[INFO]** Updating traction on region 6 → -678000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.637e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.253e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.279e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -678000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.436e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.726e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4729e+01 J
  → Fracture energy : 1.8787e-04 J
  → Total energy    : 5.4729e+01 J


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
  **[INFO]** Updating traction on region 6 → -684000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.177e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.660e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.108e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -684000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.121e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.463e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.6884e+01 J
  → Fracture energy : 1.9780e-04 J
  → Total energy    : 5.6884e+01 J


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
  **[INFO]** Updating traction on region 6 → -690000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.769e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.964e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.244e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -690000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.108e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.044e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.0527e+01 J
  → Fracture energy : 2.0968e-04 J
  → Total energy    : 6.0527e+01 J


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
  **[INFO]** Updating traction on region 6 → -696000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.852e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.560e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.293e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -696000000.0 Pa
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
  → Elastic energy  : 6.3149e+01 J
  → Fracture energy : 2.1797e-04 J
  → Total energy    : 6.3150e+01 J


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
  **[INFO]** Updating traction on region 6 → -702000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.619e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.984e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.611e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -702000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.532e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.220e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.7514e+01 J
  → Fracture energy : 2.2538e-04 J
  → Total energy    : 6.7514e+01 J


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
  **[INFO]** Updating traction on region 6 → -708000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.265e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.326e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.773e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -708000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.428e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.045e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8799e+01 J
  → Fracture energy : 2.3225e-04 J
  → Total energy    : 6.8799e+01 J


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
  **[INFO]** Updating traction on region 6 → -714000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.972e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.267e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.233e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -714000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.468e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.782e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.0877e+01 J
  → Fracture energy : 2.3742e-04 J
  → Total energy    : 7.0877e+01 J


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
  **[INFO]** Updating traction on region 6 → -720000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.414e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.252e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.804e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -720000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.163e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.561e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3046e+01 J
  → Fracture energy : 2.4233e-04 J
  → Total energy    : 7.3046e+01 J


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
  **[INFO]** Updating traction on region 6 → -726000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.457e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.555e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.977e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -726000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.295e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.255e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4288e+01 J
  → Fracture energy : 2.4442e-04 J
  → Total energy    : 7.4289e+01 J


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
  **[INFO]** Updating traction on region 6 → -732000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.705e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.741e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.095e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -732000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.831e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.104e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6179e+01 J
  → Fracture energy : 2.4748e-04 J
  → Total energy    : 7.6179e+01 J


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
  **[INFO]** Updating traction on region 6 → -738000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.737e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.562e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.512e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -738000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.176e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.325e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.7467e+01 J
  → Fracture energy : 2.4889e-04 J
  → Total energy    : 7.7467e+01 J


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
  **[INFO]** Updating traction on region 6 → -744000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.102e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.563e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.812e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -744000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.771e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.696e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.8914e+01 J
  → Fracture energy : 2.5004e-04 J
  → Total energy    : 7.8915e+01 J


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
  **[INFO]** Updating traction on region 6 → -750000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.064e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.705e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.091e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -750000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.294e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.864e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0200e+01 J
  → Fracture energy : 2.5251e-04 J
  → Total energy    : 8.0200e+01 J


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
  **[INFO]** Updating traction on region 6 → -756000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.956e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.056e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.583e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -756000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.400e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.757e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1490e+01 J
  → Fracture energy : 2.5782e-04 J
  → Total energy    : 8.1490e+01 J


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
  **[INFO]** Updating traction on region 6 → -762000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.908e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.581e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.244e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -762000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.250e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.929e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.2792e+01 J
  → Fracture energy : 2.6731e-04 J
  → Total energy    : 8.2792e+01 J


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
  **[INFO]** Updating traction on region 6 → -768000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.946e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.909e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.880e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -768000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.340e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.664e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.4156e+01 J
  → Fracture energy : 2.8074e-04 J
  → Total energy    : 8.4156e+01 J


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
  **[INFO]** Updating traction on region 6 → -774000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.130e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.174e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.362e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -774000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.931e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.770e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5626e+01 J
  → Fracture energy : 2.9053e-04 J
  → Total energy    : 8.5626e+01 J


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
  **[INFO]** Updating traction on region 6 → -780000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.054e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.794e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.140e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -780000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.021e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.609e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0314e+01 J
  → Fracture energy : 3.0113e-04 J
  → Total energy    : 9.0314e+01 J


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
  **[INFO]** Updating traction on region 6 → -786000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.187e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.892e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.220e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -786000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.246e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.177e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2867e+01 J
  → Fracture energy : 3.1054e-04 J
  → Total energy    : 9.2867e+01 J


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
  **[INFO]** Updating traction on region 6 → -792000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.705e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.163e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.615e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -792000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.833e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6140e+01 J
  → Fracture energy : 3.1684e-04 J
  → Total energy    : 9.6140e+01 J


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
  **[INFO]** Updating traction on region 6 → -798000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.637e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.578e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.094e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -798000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.961e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8415e+01 J
  → Fracture energy : 3.2347e-04 J
  → Total energy    : 9.8415e+01 J


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
  **[INFO]** Updating traction on region 6 → -804000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.675e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.878e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.253e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -804000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.980e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.828e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0001e+02 J
  → Fracture energy : 3.3161e-04 J
  → Total energy    : 1.0001e+02 J


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
  **[INFO]** Updating traction on region 6 → -810000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.350e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.698e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.525e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -810000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.485e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.532e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0152e+02 J
  → Fracture energy : 3.3686e-04 J
  → Total energy    : 1.0152e+02 J


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
  **[INFO]** Updating traction on region 6 → -816000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.796e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.590e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.504e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -816000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.747e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.447e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0405e+02 J
  → Fracture energy : 3.4193e-04 J
  → Total energy    : 1.0405e+02 J


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
  **[INFO]** Updating traction on region 6 → -822000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.091e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.087e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.239e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -822000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.467e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0565e+02 J
  → Fracture energy : 3.4806e-04 J
  → Total energy    : 1.0565e+02 J


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
  **[INFO]** Updating traction on region 6 → -828000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.325e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.017e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.560e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -828000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.003e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.092e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0721e+02 J
  → Fracture energy : 3.5354e-04 J
  → Total energy    : 1.0721e+02 J


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
  **[INFO]** Updating traction on region 6 → -834000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.764e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.990e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.087e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -834000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.205e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.402e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0885e+02 J
  → Fracture energy : 3.5911e-04 J
  → Total energy    : 1.0885e+02 J


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
  **[INFO]** Updating traction on region 6 → -840000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.319e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.675e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.736e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -840000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.819e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.368e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1055e+02 J
  → Fracture energy : 3.6583e-04 J
  → Total energy    : 1.1055e+02 J


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
  **[INFO]** Updating traction on region 6 → -846000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.547e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.416e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.926e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -846000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.511e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.781e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1215e+02 J
  → Fracture energy : 3.7355e-04 J
  → Total energy    : 1.1215e+02 J


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
  **[INFO]** Updating traction on region 6 → -852000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.296e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.333e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.713e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -852000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.491e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1377e+02 J
  → Fracture energy : 3.8331e-04 J
  → Total energy    : 1.1377e+02 J


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
  **[INFO]** Updating traction on region 6 → -858000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.320e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.417e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.387e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -858000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.766e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1540e+02 J
  → Fracture energy : 3.9486e-04 J
  → Total energy    : 1.1540e+02 J


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
  **[INFO]** Updating traction on region 6 → -864000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.948e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.103e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.637e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -864000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.177e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.612e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1712e+02 J
  → Fracture energy : 4.1376e-04 J
  → Total energy    : 1.1712e+02 J


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
  **[INFO]** Updating traction on region 6 → -870000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.982e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.796e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.394e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -870000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.211e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.303e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1876e+02 J
  → Fracture energy : 4.4576e-04 J
  → Total energy    : 1.1876e+02 J


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
  **[INFO]** Updating traction on region 6 → -876000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.890e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.226e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.614e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -876000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.415e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.507e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.786e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2041e+02 J
  → Fracture energy : 5.0654e-04 J
  → Total energy    : 1.2041e+02 J


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
  **[INFO]** Updating traction on region 6 → -882000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.829e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.540e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.502e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -882000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.881e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.310e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2207e+02 J
  → Fracture energy : 6.8678e-04 J
  → Total energy    : 1.2207e+02 J


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
  **[INFO]** Updating traction on region 6 → -888000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.286e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.702e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.980e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -888000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.183e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.876e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2376e+02 J
  → Fracture energy : 1.0132e-03 J
  → Total energy    : 1.2376e+02 J


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
  **[INFO]** Updating traction on region 6 → -894000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.912e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.036e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.135e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -894000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.091e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.676e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3372e+02 J
  → Fracture energy : 1.1257e-03 J
  → Total energy    : 1.3372e+02 J


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
  **[INFO]** Updating traction on region 6 → -900000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.649e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.039e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.327e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -900000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.974e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.806e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.181e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5265e+02 J
  → Fracture energy : 1.1857e-03 J
  → Total energy    : 1.5266e+02 J


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
  **[INFO]** Updating traction on region 6 → -906000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.039e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.081e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.746e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -906000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.186e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.535e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.735e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6628e+02 J
  → Fracture energy : 1.2190e-03 J
  → Total energy    : 1.6628e+02 J


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
  **[INFO]** Updating traction on region 6 → -912000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.856e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.282e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.466e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -912000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.146e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.959e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.496e-13

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7266e+02 J
  → Fracture energy : 1.2511e-03 J
  → Total energy    : 1.7266e+02 J


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
  **[INFO]** Updating traction on region 6 → -918000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.588e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.068e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.182e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -918000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.800e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.788e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.653e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7721e+02 J
  → Fracture energy : 1.2700e-03 J
  → Total energy    : 1.7721e+02 J


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
  **[INFO]** Updating traction on region 6 → -924000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.322e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.772e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.610e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -924000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.875e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.944e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.007e-13

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8269e+02 J
  → Fracture energy : 1.2748e-03 J
  → Total energy    : 1.8269e+02 J


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
  **[INFO]** Updating traction on region 6 → -930000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.686e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.690e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.725e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -930000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.088e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.449e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.785e-13

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8656e+02 J
  → Fracture energy : 1.2769e-03 J
  → Total energy    : 1.8656e+02 J


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
  **[INFO]** Updating traction on region 6 → -936000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.044e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.339e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.439e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -936000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.654e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.475e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8925e+02 J
  → Fracture energy : 1.2814e-03 J
  → Total energy    : 1.8925e+02 J


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
  **[INFO]** Updating traction on region 6 → -942000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.920e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.201e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.308e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -942000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.579e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.119e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9173e+02 J
  → Fracture energy : 1.2847e-03 J
  → Total energy    : 1.9173e+02 J


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
  **[INFO]** Updating traction on region 6 → -948000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.522e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.060e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.361e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -948000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.062e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.017e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9435e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 1.9435e+02 J


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
  **[INFO]** Updating traction on region 6 → -954000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.922e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.909e-11
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.234e-09

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -954000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.237e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9728e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 1.9728e+02 J


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
  **[INFO]** Updating traction on region 6 → -960000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.250e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.094e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.026e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -960000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.138e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.051e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9977e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 1.9977e+02 J


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
  **[INFO]** Updating traction on region 6 → -966000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.211e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.017e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.967e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -966000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.222e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0228e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 2.0228e+02 J


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
  **[INFO]** Updating traction on region 6 → -972000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.173e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.934e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.909e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -972000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.370e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.423e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0480e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 2.0480e+02 J


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
  **[INFO]** Updating traction on region 6 → -978000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.135e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.851e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.853e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -978000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.785e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.154e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0733e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 2.0733e+02 J


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
  **[INFO]** Updating traction on region 6 → -984000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.098e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.768e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.799e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -984000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.520e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.974e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0988e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 2.0989e+02 J


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
  **[INFO]** Updating traction on region 6 → -990000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.061e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.686e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.745e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -990000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.809e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1245e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 2.1245e+02 J


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
  **[INFO]** Updating traction on region 6 → -996000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.024e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.608e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.693e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -996000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.045e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.919e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1503e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 2.1504e+02 J


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
  **[INFO]** Updating traction on region 6 → -1002000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.988e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.534e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.643e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1002000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.690e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.123e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.109e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1763e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 2.1763e+02 J


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
  **[INFO]** Updating traction on region 6 → -1008000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.952e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.239e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.929e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1008000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.910e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2025e+02 J
  → Fracture energy : 1.2848e-03 J
  → Total energy    : 2.2025e+02 J


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
  **[INFO]** Updating traction on region 6 → -1014000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.917e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.955e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.164e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1014000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.068e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2288e+02 J
  → Fracture energy : 1.2849e-03 J
  → Total energy    : 2.2288e+02 J


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
  **[INFO]** Updating traction on region 6 → -1020000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.882e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.205e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.291e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1020000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.321e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2552e+02 J
  → Fracture energy : 1.2853e-03 J
  → Total energy    : 2.2552e+02 J


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
  **[INFO]** Updating traction on region 6 → -1026000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.848e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.003e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.193e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1026000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.987e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.817e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2818e+02 J
  → Fracture energy : 1.2864e-03 J
  → Total energy    : 2.2818e+02 J


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
  **[INFO]** Updating traction on region 6 → -1032000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.814e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.592e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.544e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1032000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.662e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.466e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3086e+02 J
  → Fracture energy : 1.2882e-03 J
  → Total energy    : 2.3086e+02 J


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
  **[INFO]** Updating traction on region 6 → -1038000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.784e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.929e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.340e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1038000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.437e-16
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
  → Elastic energy  : 2.3355e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.3355e+02 J


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
  **[INFO]** Updating traction on region 6 → -1044000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.104e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.181e-07
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.033e-06

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1044000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.928e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.505e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3634e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.3634e+02 J


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
  **[INFO]** Updating traction on region 6 → -1050000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.714e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.249e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.830e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1050000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.815e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.499e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3907e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.3907e+02 J


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
  **[INFO]** Updating traction on region 6 → -1056000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.682e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.235e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.796e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1056000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.744e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.366e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4181e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.4181e+02 J


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
  **[INFO]** Updating traction on region 6 → -1062000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.650e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.213e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.764e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1062000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.277e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4456e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.4456e+02 J


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
  **[INFO]** Updating traction on region 6 → -1068000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.618e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.182e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.732e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1068000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.327e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.444e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4733e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.4733e+02 J


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
  **[INFO]** Updating traction on region 6 → -1074000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.587e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.131e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.701e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1074000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.363e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.028e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5012e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.5012e+02 J


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
  **[INFO]** Updating traction on region 6 → -1080000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.556e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.077e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.671e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1080000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.139e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.171e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5292e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.5292e+02 J


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
  **[INFO]** Updating traction on region 6 → -1086000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.525e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.024e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.641e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1086000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.400e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.050e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5574e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.5574e+02 J


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
  **[INFO]** Updating traction on region 6 → -1092000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.495e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.973e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.612e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1092000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.245e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.401e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5857e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.5857e+02 J


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
  **[INFO]** Updating traction on region 6 → -1098000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.464e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.922e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.584e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1098000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.718e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.467e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6142e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.6142e+02 J


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
  **[INFO]** Updating traction on region 6 → -1104000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.435e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.870e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.556e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1104000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.067e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6429e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.6429e+02 J


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
  **[INFO]** Updating traction on region 6 → -1110000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.405e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.820e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.529e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1110000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.721e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6717e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.6717e+02 J


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
  **[INFO]** Updating traction on region 6 → -1116000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.376e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.771e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.503e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1116000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.233e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.183e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7006e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.7007e+02 J


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
  **[INFO]** Updating traction on region 6 → -1122000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.348e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.723e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.477e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1122000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.644e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.813e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7298e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.7298e+02 J


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
  **[INFO]** Updating traction on region 6 → -1128000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.319e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.676e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.451e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1128000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.359e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.889e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7590e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.7590e+02 J


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
  **[INFO]** Updating traction on region 6 → -1134000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.291e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.630e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.427e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1134000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.510e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.688e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7885e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.7885e+02 J


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
  **[INFO]** Updating traction on region 6 → -1140000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.263e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.585e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.402e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1140000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.388e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.105e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8180e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.8181e+02 J


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
  **[INFO]** Updating traction on region 6 → -1146000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.236e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.541e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.379e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1146000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.838e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.118e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8478e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.8478e+02 J


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
  **[INFO]** Updating traction on region 6 → -1152000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.208e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.498e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.355e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1152000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.225e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.828e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8777e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.8777e+02 J


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
  **[INFO]** Updating traction on region 6 → -1158000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.181e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.456e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.333e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1158000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.576e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.904e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9077e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.9078e+02 J


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
  **[INFO]** Updating traction on region 6 → -1164000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.155e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.415e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.310e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1164000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.160e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.888e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9379e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.9380e+02 J


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
  **[INFO]** Updating traction on region 6 → -1170000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.128e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.374e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.288e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1170000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.656e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9683e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.9683e+02 J


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
  **[INFO]** Updating traction on region 6 → -1176000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.102e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.335e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.267e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1176000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.281e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.954e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9988e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 2.9989e+02 J


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
  **[INFO]** Updating traction on region 6 → -1182000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.076e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.296e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.246e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1182000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.761e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.866e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0295e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 3.0295e+02 J


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
  **[INFO]** Updating traction on region 6 → -1188000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.051e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.258e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.226e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1188000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.494e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0604e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 3.0604e+02 J


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
  **[INFO]** Updating traction on region 6 → -1194000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.025e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.221e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.206e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1194000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.019e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.979e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0913e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 3.0914e+02 J


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
  **[INFO]** Updating traction on region 6 → -1200000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.000e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.185e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.186e-10

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1200000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=2.65e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.150e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1225e+02 J
  → Fracture energy : 1.2885e-03 J
  → Total energy    : 3.1225e+02 J

Simulation completed in 153.53 s
Total time steps solved: 201
