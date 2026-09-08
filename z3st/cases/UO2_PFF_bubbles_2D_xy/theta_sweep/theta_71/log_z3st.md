Info    : Reading 'mesh.msh'...
Info    : 13 entities
Info    : 1179 nodes
Info    : 2356 elements
Info    : Done reading 'mesh.msh'


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
  lc                  : 1e-06
  hybrid_constraint   : True
[spine.load_materials]
Material loaded: uo2
  → k defined as constant: 5.0
  → Gc not defined for uo2
  - Material 'uo2': Gc (AT1) from sigma_c = 4.22e+08 Pa
  → constitutive model: lame
  E               → 358000000000.0 (float)
  G               → 145528455284.55286 (float)
  Gc              → 1.3233687150837987 (float)
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0831e-09 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.0831e-09 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 3.724e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6332e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.6332e-08 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.539e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1748e-08 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.1748e-08 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.219e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4533e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4533e-07 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2708e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2708e-07 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2699e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.2699e-07 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 2.292e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4507e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.4507e-07 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.8132e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.8132e-07 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3573e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.3573e-07 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 2.229e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0831e-07 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.0831e-07 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 9.429e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0991e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0991e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.539e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3080e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3080e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5350e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5350e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7803e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.7803e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0437e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0437e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.219e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3253e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.3253e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.306e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6250e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.6250e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.928e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9429e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.9429e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2790e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.2790e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6332e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.6332e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 3.926e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0056e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.0056e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3962e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.3962e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8050e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 4.8050e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2319e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.2319e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 6.880e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.6769e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 5.6769e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.616e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1402e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.1402e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6216e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 6.6216e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 2.292e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1212e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.1212e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6389e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 7.6389e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1748e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.1748e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 2.420e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7289e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 8.7289e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.219e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3011e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.3011e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 8.404e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8915e-06 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 9.8915e-06 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0500e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.0500e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1127e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1127e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.923e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1772e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.1772e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 3.349e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2435e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.2435e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 3.953e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3116e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3116e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 3.776e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3815e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.3815e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 2.229e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4533e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.4533e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 4.355e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5269e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.5269e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 3.926e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6023e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.6023e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 2.197e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6795e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.6795e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 9.429e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7585e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.7585e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8393e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.8393e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.220e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9220e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 1.9220e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0065e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0065e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0927e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.0927e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1809e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.1809e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.500e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2708e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.2708e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 1.106e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3625e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.3625e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 2.867e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4561e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.4561e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5514e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.5514e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6486e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.6486e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 3.384e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7476e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.7476e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8485e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.8485e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  ||Δu||/||u|| = 2.265e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9511e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 2.9511e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0556e-05 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 3.0556e-05 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 2.435e-03

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.500e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 2.313e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.928e-01
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 2.417e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.529e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 2.076e-03

Convergence check


#### Iteration 5/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.912e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 2.170e-03

Convergence check


#### Iteration 6/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.264e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 2.256e-03

Convergence check


#### Iteration 7/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.574e-01
  [adaptive] relax_D=0.05
  |ΔD|_∞ = 2.331e-03

Convergence check


#### Iteration 8/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.715e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.635e-03

Convergence check


#### Iteration 9/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.017e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.708e-03

Convergence check


#### Iteration 10/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.294e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 1.776e-03

Convergence check


#### Iteration 11/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.538e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 1.835e-03

Convergence check


#### Iteration 12/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.740e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 1.884e-03

Convergence check


#### Iteration 13/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.891e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.921e-03

Convergence check


#### Iteration 14/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.079e-01
  [adaptive] relax_D=0.06
  |ΔD|_∞ = 1.237e-03

Convergence check


#### Iteration 15/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.272e-01
  [adaptive] relax_D=0.07
  |ΔD|_∞ = 1.283e-03

Convergence check


#### Iteration 16/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.439e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 1.324e-03

Convergence check


#### Iteration 17/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.575e-01
  [adaptive] relax_D=0.08
  |ΔD|_∞ = 1.357e-03

Convergence check


#### Iteration 18/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.673e-01
  [adaptive] relax_D=0.09
  |ΔD|_∞ = 1.381e-03

Convergence check


#### Iteration 19/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.725e-01
  [adaptive] relax_D=0.10
  |ΔD|_∞ = 1.394e-03

Convergence check


#### Iteration 20/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.726e-01
  [adaptive] relax_D=0.11
  |ΔD|_∞ = 1.394e-03

Convergence check


#### Iteration 21/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.669e-01
  [adaptive] relax_D=0.12
  |ΔD|_∞ = 1.380e-03

Convergence check


#### Iteration 22/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.551e-01
  [adaptive] relax_D=0.13
  |ΔD|_∞ = 1.352e-03

Convergence check


#### Iteration 23/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.368e-01
  [adaptive] relax_D=0.15
  |ΔD|_∞ = 1.307e-03

Convergence check


#### Iteration 24/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.120e-01
  [adaptive] relax_D=0.16
  |ΔD|_∞ = 1.247e-03

Convergence check


#### Iteration 25/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.809e-01
  [adaptive] relax_D=0.18
  |ΔD|_∞ = 1.171e-03

Convergence check


#### Iteration 26/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.439e-01
  [adaptive] relax_D=0.19
  |ΔD|_∞ = 1.081e-03

Convergence check


#### Iteration 27/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.019e-01
  [adaptive] relax_D=0.21
  |ΔD|_∞ = 9.786e-04

Convergence check


#### Iteration 28/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.561e-01
  [adaptive] relax_D=0.24
  |ΔD|_∞ = 8.670e-04

Convergence check


#### Iteration 29/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.078e-01
  [adaptive] relax_D=0.26
  |ΔD|_∞ = 7.495e-04

Convergence check


#### Iteration 30/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.589e-01
  [adaptive] relax_D=0.28
  |ΔD|_∞ = 6.303e-04

Convergence check


#### Iteration 31/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.110e-01
  [adaptive] relax_D=0.31
  |ΔD|_∞ = 5.138e-04

Convergence check


#### Iteration 32/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.660e-01
  [adaptive] relax_D=0.34
  |ΔD|_∞ = 4.041e-04

Convergence check


#### Iteration 33/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.254e-01
  [adaptive] relax_D=0.38
  |ΔD|_∞ = 3.052e-04

Convergence check


#### Iteration 34/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.036e-02
  [adaptive] relax_D=0.42
  |ΔD|_∞ = 2.200e-04

Convergence check


#### Iteration 35/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.170e-02
  [adaptive] relax_D=0.46
  |ΔD|_∞ = 1.502e-04

Convergence check


#### Iteration 36/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.956e-02
  [adaptive] relax_D=0.50
  |ΔD|_∞ = 9.632e-05

Convergence check


#### Iteration 37/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -354000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.266e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.355e-02
  [adaptive] relax_D=0.56
  |ΔD|_∞ = 5.734e-05

Convergence check


#### Iteration 38/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.283e-02
  [adaptive] relax_D=0.61
  |ΔD|_∞ = 3.124e-05

Convergence check


#### Iteration 39/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.277e-03
  [adaptive] relax_D=0.67
  |ΔD|_∞ = 1.528e-05

Convergence check


#### Iteration 40/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.688e-03
  [adaptive] relax_D=0.74
  |ΔD|_∞ = 6.544e-06

Convergence check


#### Iteration 41/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.704e-04
  [adaptive] relax_D=0.81
  |ΔD|_∞ = 2.363e-06

Convergence check


#### Iteration 42/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.786e-04
  [adaptive] relax_D=0.89
  |ΔD|_∞ = 6.783e-07

Convergence check


#### Iteration 43/200


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.735e-05
  [adaptive] relax_D=0.98
  |ΔD|_∞ = 1.396e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 43 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1602e-05 J
  → Fracture energy : 4.6386e-09 J
  → Total energy    : 3.1606e-05 J


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
  ||Δu||/||u|| = 1.711e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.450e-01
  [adaptive] relax_D=0.98
  |ΔD|_∞ = 9.527e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -360000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.223e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.564e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -360000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.265e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.041e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.610e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -360000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.963e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2658e-05 J
  → Fracture energy : 3.3197e-08 J
  → Total energy    : 3.2691e-05 J


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
  ||Δu||/||u|| = 1.817e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.322e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.240e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -366000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3699e-05 J
  → Fracture energy : 1.1449e-07 J
  → Total energy    : 3.3813e-05 J


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
  ||Δu||/||u|| = 2.278e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.583e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.918e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -372000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.192e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.374e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.585e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4792e-05 J
  → Fracture energy : 3.3077e-07 J
  → Total energy    : 3.5123e-05 J


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
  ||Δu||/||u|| = 3.274e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.077e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.671e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -378000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.634e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5976e-05 J
  → Fracture energy : 8.0092e-07 J
  → Total energy    : 3.6777e-05 J


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
  ||Δu||/||u|| = 5.278e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.551e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.529e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -384000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7388e-05 J
  → Fracture energy : 1.6324e-06 J
  → Total energy    : 3.9021e-05 J


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
  ||Δu||/||u|| = 1.413e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.884e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.357e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -390000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9356e-05 J
  → Fracture energy : 2.8516e-06 J
  → Total energy    : 4.2208e-05 J


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
  ||Δu||/||u|| = 9.995e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.385e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.587e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -396000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.853e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.106e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0059e-03 J
  → Fracture energy : 4.5105e-06 J
  → Total energy    : 1.0104e-03 J


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
  ||Δu||/||u|| = 9.701e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.072e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.713e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -402000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4277e-01 J
  → Fracture energy : 6.7711e-06 J
  → Total energy    : 1.4278e-01 J


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
  ||Δu||/||u|| = 7.007e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.878e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.656e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -408000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.090e-21
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.044e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9733e-01 J
  → Fracture energy : 9.6538e-06 J
  → Total energy    : 3.9734e-01 J


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
  ||Δu||/||u|| = 5.428e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.708e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.788e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -414000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4317e-01 J
  → Fracture energy : 1.3090e-05 J
  → Total energy    : 7.4318e-01 J


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
  ||Δu||/||u|| = 4.684e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.499e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.620e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -420000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.004e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.111e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3028e+00 J
  → Fracture energy : 1.7122e-05 J
  → Total energy    : 1.3028e+00 J


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
  ||Δu||/||u|| = 3.972e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.287e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.546e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -426000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9985e+00 J
  → Fracture energy : 2.1805e-05 J
  → Total energy    : 1.9985e+00 J


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
  ||Δu||/||u|| = 2.912e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.114e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.476e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6850e+00 J
  → Fracture energy : 2.7267e-05 J
  → Total energy    : 2.6850e+00 J


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
  ||Δu||/||u|| = 3.157e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.950e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.414e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -438000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.183e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.288e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7744e+00 J
  → Fracture energy : 3.3356e-05 J
  → Total energy    : 3.7745e+00 J


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
  ||Δu||/||u|| = 2.066e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.840e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.344e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -444000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.462e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.186e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6482e+00 J
  → Fracture energy : 3.9864e-05 J
  → Total energy    : 4.6482e+00 J


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
  ||Δu||/||u|| = 2.559e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.746e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.544e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -450000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.298e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.238e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1417e+00 J
  → Fracture energy : 4.7775e-05 J
  → Total energy    : 6.1418e+00 J


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
  ||Δu||/||u|| = 1.602e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.670e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.300e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -456000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1132e+00 J
  → Fracture energy : 5.6101e-05 J
  → Total energy    : 7.1132e+00 J


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
  ||Δu||/||u|| = 1.962e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.575e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.796e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -462000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.680e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.561e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.8155e+00 J
  → Fracture energy : 6.5430e-05 J
  → Total energy    : 8.8155e+00 J


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
  ||Δu||/||u|| = 1.555e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.554e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.022e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -468000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.860e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.970e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0269e+01 J
  → Fracture energy : 7.5823e-05 J
  → Total energy    : 1.0269e+01 J


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
  ||Δu||/||u|| = 1.325e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.590e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.142e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -474000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.207e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1851e+01 J
  → Fracture energy : 8.7989e-05 J
  → Total energy    : 1.1851e+01 J


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
  ||Δu||/||u|| = 1.365e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.427e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.289e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -480000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.454e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.244e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3573e+01 J
  → Fracture energy : 9.8110e-05 J
  → Total energy    : 1.3573e+01 J


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
  ||Δu||/||u|| = 2.090e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.087e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.851e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -486000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.618e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.122e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6990e+01 J
  → Fracture energy : 1.0852e-04 J
  → Total energy    : 1.6990e+01 J


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
  ||Δu||/||u|| = 1.292e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.717e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.076e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -492000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.679e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.193e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9527e+01 J
  → Fracture energy : 1.2131e-04 J
  → Total energy    : 1.9527e+01 J


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
  ||Δu||/||u|| = 9.448e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.114e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.808e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -498000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.253e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.231e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1453e+01 J
  → Fracture energy : 1.3393e-04 J
  → Total energy    : 2.1453e+01 J


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
  ||Δu||/||u|| = 1.078e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.543e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.769e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -504000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3904e+01 J
  → Fracture energy : 1.4641e-04 J
  → Total energy    : 2.3904e+01 J


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
  ||Δu||/||u|| = 8.980e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.668e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.261e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -510000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.104e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.221e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6370e+01 J
  → Fracture energy : 1.5729e-04 J
  → Total energy    : 2.6370e+01 J


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
  ||Δu||/||u|| = 8.344e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.469e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.730e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -516000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.524e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.603e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8660e+01 J
  → Fracture energy : 1.6759e-04 J
  → Total energy    : 2.8660e+01 J


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
  ||Δu||/||u|| = 7.734e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.116e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.506e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -522000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1145e+01 J
  → Fracture energy : 1.7793e-04 J
  → Total energy    : 3.1145e+01 J


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
  ||Δu||/||u|| = 9.560e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.691e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.781e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -528000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3936e+01 J
  → Fracture energy : 1.8892e-04 J
  → Total energy    : 3.3936e+01 J


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
  ||Δu||/||u|| = 4.960e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.889e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.598e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -534000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.286e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.757e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5436e+01 J
  → Fracture energy : 1.9688e-04 J
  → Total energy    : 3.5436e+01 J


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
  ||Δu||/||u|| = 1.185e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.389e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.840e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -540000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9197e+01 J
  → Fracture energy : 2.0431e-04 J
  → Total energy    : 3.9197e+01 J


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
  ||Δu||/||u|| = 2.438e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.887e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.048e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0211e+01 J
  → Fracture energy : 2.1224e-04 J
  → Total energy    : 4.0211e+01 J


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
  ||Δu||/||u|| = 4.072e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.230e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.000e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -552000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.484e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.294e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.901e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1791e+01 J
  → Fracture energy : 2.2311e-04 J
  → Total energy    : 4.1791e+01 J


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
  ||Δu||/||u|| = 2.040e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.592e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.407e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -558000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.749e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.888e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2786e+01 J
  → Fracture energy : 2.3595e-04 J
  → Total energy    : 4.2787e+01 J


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
  ||Δu||/||u|| = 4.901e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.013e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.586e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -564000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.367e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.159e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4766e+01 J
  → Fracture energy : 2.4936e-04 J
  → Total energy    : 4.4766e+01 J


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
  ||Δu||/||u|| = 3.062e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.669e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.636e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -570000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.381e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.674e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6027e+01 J
  → Fracture energy : 2.5876e-04 J
  → Total energy    : 4.6027e+01 J


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
  ||Δu||/||u|| = 9.446e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.357e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.944e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -576000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0096e+01 J
  → Fracture energy : 2.6528e-04 J
  → Total energy    : 5.0096e+01 J


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
  ||Δu||/||u|| = 5.005e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.957e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.896e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -582000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2574e+01 J
  → Fracture energy : 2.7341e-04 J
  → Total energy    : 5.2574e+01 J


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
  ||Δu||/||u|| = 2.067e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.824e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.713e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4236e+01 J
  → Fracture energy : 2.8147e-04 J
  → Total energy    : 5.4236e+01 J


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
  ||Δu||/||u|| = 1.460e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.629e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.742e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -594000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.277e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.140e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5391e+01 J
  → Fracture energy : 2.8641e-04 J
  → Total energy    : 5.5392e+01 J


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
  ||Δu||/||u|| = 6.037e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.608e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.237e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.8326e+01 J
  → Fracture energy : 2.9284e-04 J
  → Total energy    : 5.8327e+01 J


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
  ||Δu||/||u|| = 1.180e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.665e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.937e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -606000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9520e+01 J
  → Fracture energy : 2.9921e-04 J
  → Total energy    : 5.9521e+01 J


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
  ||Δu||/||u|| = 2.051e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.985e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.851e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -612000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.205e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.956e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1111e+01 J
  → Fracture energy : 3.0829e-04 J
  → Total energy    : 6.1112e+01 J


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
  ||Δu||/||u|| = 1.251e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.567e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.844e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -618000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2475e+01 J
  → Fracture energy : 3.1556e-04 J
  → Total energy    : 6.2475e+01 J


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
  ||Δu||/||u|| = 5.962e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.900e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.590e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -624000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.127e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.114e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.5218e+01 J
  → Fracture energy : 3.2105e-04 J
  → Total energy    : 6.5219e+01 J


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
  ||Δu||/||u|| = 1.404e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.378e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.591e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -630000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6798e+01 J
  → Fracture energy : 3.2847e-04 J
  → Total energy    : 6.6798e+01 J


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
  ||Δu||/||u|| = 9.841e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.156e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.759e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -636000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8082e+01 J
  → Fracture energy : 3.3805e-04 J
  → Total energy    : 6.8083e+01 J


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
  ||Δu||/||u|| = 1.102e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.635e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.387e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9395e+01 J
  → Fracture energy : 3.5167e-04 J
  → Total energy    : 6.9396e+01 J


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
  ||Δu||/||u|| = 1.000e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.909e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.725e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -648000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.0713e+01 J
  → Fracture energy : 3.7171e-04 J
  → Total energy    : 7.0714e+01 J


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
  ||Δu||/||u|| = 1.030e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.129e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.676e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -654000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.2045e+01 J
  → Fracture energy : 4.0885e-04 J
  → Total energy    : 7.2046e+01 J


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
  ||Δu||/||u|| = 1.889e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.399e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.247e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3424e+01 J
  → Fracture energy : 4.9216e-04 J
  → Total energy    : 7.3424e+01 J


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
  ||Δu||/||u|| = 1.176e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.729e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.329e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -666000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.171e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.217e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.948e-13

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4784e+01 J
  → Fracture energy : 6.5835e-04 J
  → Total energy    : 7.4785e+01 J


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
  ||Δu||/||u|| = 1.403e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.369e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.269e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0589e+01 J
  → Fracture energy : 7.6375e-04 J
  → Total energy    : 8.0590e+01 J


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
  ||Δu||/||u|| = 2.864e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.393e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.658e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -678000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.386e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.793e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.2779e+01 J
  → Fracture energy : 8.3791e-04 J
  → Total energy    : 8.2780e+01 J


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
  ||Δu||/||u|| = 2.523e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.171e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.274e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -684000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.638e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.522e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.039e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6069e+01 J
  → Fracture energy : 8.7663e-04 J
  → Total energy    : 9.6070e+01 J


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
  ||Δu||/||u|| = 9.717e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.156e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.615e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -690000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.780e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.094e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.803e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0389e+02 J
  → Fracture energy : 8.9261e-04 J
  → Total energy    : 1.0389e+02 J


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
  ||Δu||/||u|| = 7.117e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.266e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.276e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -696000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.091e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.330e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1083e+02 J
  → Fracture energy : 9.0386e-04 J
  → Total energy    : 1.1083e+02 J


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
  ||Δu||/||u|| = 2.122e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.317e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.633e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -702000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.186e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.518e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.429e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1366e+02 J
  → Fracture energy : 9.0965e-04 J
  → Total energy    : 1.1366e+02 J


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
  ||Δu||/||u|| = 1.563e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.893e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.320e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -708000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.270e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.547e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1616e+02 J
  → Fracture energy : 9.1188e-04 J
  → Total energy    : 1.1616e+02 J


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
  ||Δu||/||u|| = 1.366e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.525e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.375e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -714000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.841e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.127e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.632e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1854e+02 J
  → Fracture energy : 9.1771e-04 J
  → Total energy    : 1.1854e+02 J


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
  ||Δu||/||u|| = 8.334e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.549e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.962e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -720000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2054e+02 J
  → Fracture energy : 9.2469e-04 J
  → Total energy    : 1.2054e+02 J


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
  ||Δu||/||u|| = 9.989e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.085e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.153e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -726000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.667e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.105e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2262e+02 J
  → Fracture energy : 9.3168e-04 J
  → Total energy    : 1.2262e+02 J


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
  ||Δu||/||u|| = 1.269e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.473e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.431e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -732000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2500e+02 J
  → Fracture energy : 9.3639e-04 J
  → Total energy    : 1.2500e+02 J


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
  ||Δu||/||u|| = 1.925e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.627e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.174e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -738000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.202e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.282e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2787e+02 J
  → Fracture energy : 9.4154e-04 J
  → Total energy    : 1.2787e+02 J


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
  ||Δu||/||u|| = 8.289e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.177e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.796e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -744000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.031e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2996e+02 J
  → Fracture energy : 9.4466e-04 J
  → Total energy    : 1.2996e+02 J


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
  ||Δu||/||u|| = 2.228e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.403e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.531e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -750000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.375e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3261e+02 J
  → Fracture energy : 9.4504e-04 J
  → Total energy    : 1.3261e+02 J


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
  ||Δu||/||u|| = 1.101e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.687e-07
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.449e-06

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -756000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.764e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.784e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3501e+02 J
  → Fracture energy : 9.4504e-04 J
  → Total energy    : 1.3501e+02 J


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
  ||Δu||/||u|| = 7.874e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.328e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.009e-11

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -762000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.135e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.304e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3717e+02 J
  → Fracture energy : 9.4504e-04 J
  → Total energy    : 1.3717e+02 J


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
  ||Δu||/||u|| = 7.812e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.325e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.683e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -768000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.825e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.427e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3933e+02 J
  → Fracture energy : 9.4504e-04 J
  → Total energy    : 1.3933e+02 J


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
  ||Δu||/||u|| = 7.752e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.134e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.302e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -774000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.334e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.823e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4152e+02 J
  → Fracture energy : 9.4504e-04 J
  → Total energy    : 1.4152e+02 J


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
  ||Δu||/||u|| = 7.692e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.151e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.718e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -780000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4372e+02 J
  → Fracture energy : 9.4506e-04 J
  → Total energy    : 1.4372e+02 J


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
  ||Δu||/||u|| = 7.634e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.026e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.315e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -786000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4594e+02 J
  → Fracture energy : 9.4525e-04 J
  → Total energy    : 1.4594e+02 J


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
  ||Δu||/||u|| = 7.576e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.916e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.265e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -792000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.774e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.492e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.820e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4818e+02 J
  → Fracture energy : 9.4601e-04 J
  → Total energy    : 1.4818e+02 J


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
  ||Δu||/||u|| = 7.519e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.040e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.359e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.031e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5043e+02 J
  → Fracture energy : 9.4828e-04 J
  → Total energy    : 1.5043e+02 J


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
  ||Δu||/||u|| = 7.463e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.044e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.374e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -804000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.031e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5270e+02 J
  → Fracture energy : 9.5095e-04 J
  → Total energy    : 1.5270e+02 J


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
  ||Δu||/||u|| = 7.501e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.216e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.332e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -810000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.644e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.063e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5499e+02 J
  → Fracture energy : 9.5115e-04 J
  → Total energy    : 1.5499e+02 J


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
  ||Δu||/||u|| = 1.168e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.023e-08
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.305e-07

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -816000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.038e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.819e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5763e+02 J
  → Fracture energy : 9.5115e-04 J
  → Total energy    : 1.5763e+02 J


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
  ||Δu||/||u|| = 7.299e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.667e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.620e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -822000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.395e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.629e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5995e+02 J
  → Fracture energy : 9.5119e-04 J
  → Total energy    : 1.5995e+02 J


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
  ||Δu||/||u|| = 7.246e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.298e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.196e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -828000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6230e+02 J
  → Fracture energy : 9.5142e-04 J
  → Total energy    : 1.6230e+02 J


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
  ||Δu||/||u|| = 7.194e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.944e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.275e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -834000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.168e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.483e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6466e+02 J
  → Fracture energy : 9.5213e-04 J
  → Total energy    : 1.6466e+02 J


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
  ||Δu||/||u|| = 7.143e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.555e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.087e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -840000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.435e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6703e+02 J
  → Fracture energy : 9.5414e-04 J
  → Total energy    : 1.6703e+02 J


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
  ||Δu||/||u|| = 7.092e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.320e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.232e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -846000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6943e+02 J
  → Fracture energy : 9.5763e-04 J
  → Total energy    : 1.6943e+02 J


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
  ||Δu||/||u|| = 7.045e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.951e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.698e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7184e+02 J
  → Fracture energy : 9.5934e-04 J
  → Total energy    : 1.7184e+02 J


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
  ||Δu||/||u|| = 9.449e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.159e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.738e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7457e+02 J
  → Fracture energy : 9.6179e-04 J
  → Total energy    : 1.7458e+02 J


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
  ||Δu||/||u|| = 6.945e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.164e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.650e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -864000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7702e+02 J
  → Fracture energy : 9.6490e-04 J
  → Total energy    : 1.7703e+02 J


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
  ||Δu||/||u|| = 6.948e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.825e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.458e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -870000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7949e+02 J
  → Fracture energy : 9.6647e-04 J
  → Total energy    : 1.7949e+02 J


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
  ||Δu||/||u|| = 1.094e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.355e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.357e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -876000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8228e+02 J
  → Fracture energy : 9.6914e-04 J
  → Total energy    : 1.8228e+02 J


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
  ||Δu||/||u|| = 6.803e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.646e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.893e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -882000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.157e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.327e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8478e+02 J
  → Fracture energy : 9.7133e-04 J
  → Total energy    : 1.8478e+02 J


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
  ||Δu||/||u|| = 1.240e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.150e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.724e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -888000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.825e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.303e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8753e+02 J
  → Fracture energy : 9.7386e-04 J
  → Total energy    : 1.8753e+02 J


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
  ||Δu||/||u|| = 7.406e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.054e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.412e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -894000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.822e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.290e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.683e-13

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9017e+02 J
  → Fracture energy : 9.7658e-04 J
  → Total energy    : 1.9017e+02 J


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
  ||Δu||/||u|| = 6.734e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.721e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.512e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -900000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.442e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.714e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9273e+02 J
  → Fracture energy : 9.7761e-04 J
  → Total energy    : 1.9273e+02 J


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
  ||Δu||/||u|| = 9.672e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.065e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.935e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -906000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.939e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9558e+02 J
  → Fracture energy : 9.7947e-04 J
  → Total energy    : 1.9558e+02 J


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
  ||Δu||/||u|| = 6.579e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.438e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.369e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -912000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.143e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.141e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9818e+02 J
  → Fracture energy : 9.8359e-04 J
  → Total energy    : 1.9818e+02 J


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
  ||Δu||/||u|| = 6.538e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.795e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.696e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -918000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0080e+02 J
  → Fracture energy : 9.8814e-04 J
  → Total energy    : 2.0080e+02 J


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
  ||Δu||/||u|| = 8.961e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.234e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.668e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -924000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.425e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0367e+02 J
  → Fracture energy : 9.9021e-04 J
  → Total energy    : 2.0367e+02 J


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
  ||Δu||/||u|| = 1.331e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.156e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.890e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -930000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0663e+02 J
  → Fracture energy : 9.9025e-04 J
  → Total energy    : 2.0663e+02 J


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
  ||Δu||/||u|| = 7.906e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.733e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.678e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -936000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.675e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.722e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0958e+02 J
  → Fracture energy : 9.9039e-04 J
  → Total energy    : 2.0959e+02 J


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
  ||Δu||/||u|| = 6.369e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.508e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.141e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -942000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1228e+02 J
  → Fracture energy : 9.9101e-04 J
  → Total energy    : 2.1228e+02 J


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
  ||Δu||/||u|| = 6.329e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.733e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.196e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -948000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.295e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.586e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.382e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1499e+02 J
  → Fracture energy : 9.9306e-04 J
  → Total energy    : 2.1499e+02 J


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
  ||Δu||/||u|| = 6.289e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.533e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.943e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1772e+02 J
  → Fracture energy : 9.9786e-04 J
  → Total energy    : 2.1772e+02 J


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
  ||Δu||/||u|| = 6.251e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.403e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.355e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -960000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.508e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.018e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.271e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2047e+02 J
  → Fracture energy : 1.0022e-03 J
  → Total energy    : 2.2047e+02 J


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
  ||Δu||/||u|| = 8.079e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.497e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.117e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.844e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2345e+02 J
  → Fracture energy : 1.0036e-03 J
  → Total energy    : 2.2346e+02 J


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
  ||Δu||/||u|| = 9.577e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.479e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.751e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -972000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.316e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.724e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2651e+02 J
  → Fracture energy : 1.0052e-03 J
  → Total energy    : 2.2652e+02 J


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
  ||Δu||/||u|| = 6.136e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.416e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.596e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -978000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.506e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.299e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.018e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2932e+02 J
  → Fracture energy : 1.0086e-03 J
  → Total energy    : 2.2932e+02 J


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
  ||Δu||/||u|| = 6.099e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.902e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.531e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -984000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.518e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.162e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3214e+02 J
  → Fracture energy : 1.0101e-03 J
  → Total energy    : 2.3214e+02 J


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
  ||Δu||/||u|| = 8.239e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.187e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.178e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3528e+02 J
  → Fracture energy : 1.0123e-03 J
  → Total energy    : 2.3528e+02 J


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
  ||Δu||/||u|| = 6.049e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.737e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.431e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -996000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3815e+02 J
  → Fracture energy : 1.0176e-03 J
  → Total energy    : 2.3815e+02 J


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
  ||Δu||/||u|| = 5.989e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.459e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.227e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1002000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.893e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.337e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4103e+02 J
  → Fracture energy : 1.0226e-03 J
  → Total energy    : 2.4103e+02 J


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
  ||Δu||/||u|| = 8.093e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.005e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.510e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1008000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.897e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.065e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4413e+02 J
  → Fracture energy : 1.0232e-03 J
  → Total energy    : 2.4413e+02 J


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
  ||Δu||/||u|| = 8.333e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.367e-07
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.744e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1014000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.583e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.028e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4748e+02 J
  → Fracture energy : 1.0232e-03 J
  → Total energy    : 2.4748e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.089e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.292e-12

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1020000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.661e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.099e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5042e+02 J
  → Fracture energy : 1.0232e-03 J
  → Total energy    : 2.5042e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.592e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.957e-11

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1026000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5337e+02 J
  → Fracture energy : 1.0232e-03 J
  → Total energy    : 2.5337e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.414e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.126e-11

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1032000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5634e+02 J
  → Fracture energy : 1.0232e-03 J
  → Total energy    : 2.5634e+02 J


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
  ||Δu||/||u|| = 5.780e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.144e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.346e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1038000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5933e+02 J
  → Fracture energy : 1.0232e-03 J
  → Total energy    : 2.5933e+02 J


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
  ||Δu||/||u|| = 5.747e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.626e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.281e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1044000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.375e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.419e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6234e+02 J
  → Fracture energy : 1.0235e-03 J
  → Total energy    : 2.6234e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.784e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1050000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6536e+02 J
  → Fracture energy : 1.0244e-03 J
  → Total energy    : 2.6536e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.104e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.588e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1056000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.204e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.230e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.050e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6840e+02 J
  → Fracture energy : 1.0267e-03 J
  → Total energy    : 2.6840e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.174e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.815e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7146e+02 J
  → Fracture energy : 1.0297e-03 J
  → Total energy    : 2.7146e+02 J


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
  ||Δu||/||u|| = 5.634e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.628e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.117e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1068000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7454e+02 J
  → Fracture energy : 1.0307e-03 J
  → Total energy    : 2.7454e+02 J


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
  ||Δu||/||u|| = 7.824e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.740e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.019e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1074000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.312e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.952e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7791e+02 J
  → Fracture energy : 1.0331e-03 J
  → Total energy    : 2.7791e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.819e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.261e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1080000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.565e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8102e+02 J
  → Fracture energy : 1.0389e-03 J
  → Total energy    : 2.8102e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.372e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.560e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1086000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.784e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.534e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8415e+02 J
  → Fracture energy : 1.0441e-03 J
  → Total energy    : 2.8415e+02 J


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
  ||Δu||/||u|| = 5.939e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.504e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.999e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1092000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.613e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.694e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.853e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8732e+02 J
  → Fracture energy : 1.0463e-03 J
  → Total energy    : 2.8732e+02 J


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
  ||Δu||/||u|| = 1.216e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.747e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.015e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1098000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9114e+02 J
  → Fracture energy : 1.0518e-03 J
  → Total energy    : 2.9114e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.028e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.220e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9433e+02 J
  → Fracture energy : 1.0602e-03 J
  → Total energy    : 2.9433e+02 J


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
  ||Δu||/||u|| = 5.457e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.740e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.550e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1110000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 5.468e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.814e-13
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.169e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9754e+02 J
  → Fracture energy : 1.0649e-03 J
  → Total energy    : 2.9754e+02 J


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
  ||Δu||/||u|| = 1.557e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.164e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.192e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1116000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0157e+02 J
  → Fracture energy : 1.0679e-03 J
  → Total energy    : 3.0157e+02 J


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
  ||Δu||/||u|| = 1.194e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.732e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.039e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1122000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0540e+02 J
  → Fracture energy : 1.0726e-03 J
  → Total energy    : 3.0540e+02 J


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
  ||Δu||/||u|| = 5.385e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.396e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.317e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1128000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0868e+02 J
  → Fracture energy : 1.0765e-03 J
  → Total energy    : 3.0869e+02 J


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
  ||Δu||/||u|| = 6.852e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 5.092e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.609e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1134000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.227e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1213e+02 J
  → Fracture energy : 1.0778e-03 J
  → Total energy    : 3.1213e+02 J


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
  ||Δu||/||u|| = 9.341e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.297e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.033e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1140000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1583e+02 J
  → Fracture energy : 1.0797e-03 J
  → Total energy    : 3.1583e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.348e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.399e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1146000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1916e+02 J
  → Fracture energy : 1.0831e-03 J
  → Total energy    : 3.1916e+02 J


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
  ||Δu||/||u|| = 5.211e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.069e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.578e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1152000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.149e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.096e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2251e+02 J
  → Fracture energy : 1.0841e-03 J
  → Total energy    : 3.2252e+02 J


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
  ||Δu||/||u|| = 7.572e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 9.024e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.945e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1158000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2620e+02 J
  → Fracture energy : 1.0843e-03 J
  → Total energy    : 3.2620e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.021e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.858e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1164000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.804e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2959e+02 J
  → Fracture energy : 1.0849e-03 J
  → Total energy    : 3.2959e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 8.121e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.650e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1170000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 8.154e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.953e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3299e+02 J
  → Fracture energy : 1.0866e-03 J
  → Total energy    : 3.3300e+02 J


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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.337e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.365e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1176000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.325e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 6.527e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.130e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3642e+02 J
  → Fracture energy : 1.0899e-03 J
  → Total energy    : 3.3642e+02 J


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
  ||Δu||/||u|| = 5.077e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 4.503e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.407e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1182000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.922e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 2.813e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.184e-13

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3986e+02 J
  → Fracture energy : 1.0914e-03 J
  → Total energy    : 3.3986e+02 J


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
  ||Δu||/||u|| = 6.695e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.946e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.592e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1188000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.354e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.139e-11
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.524e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4359e+02 J
  → Fracture energy : 1.0931e-03 J
  → Total energy    : 3.4360e+02 J


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
  ||Δu||/||u|| = 5.041e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 1.468e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.776e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 6 → -1194000000.0 Pa
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.487e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.231e-12
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.055e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4708e+02 J
  → Fracture energy : 1.0968e-03 J
  → Total energy    : 3.4708e+02 J


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
  ||Δu||/||u|| = 5.001e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 7.760e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.164e-01

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
  - Material 'uo2': AT1 solve. Gc=1.32e+00, sigma_c=4.22e+08
  ||ΔD||/||D|| = 3.399e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5058e+02 J
  → Fracture energy : 1.0995e-03 J
  → Total energy    : 3.5058e+02 J

Simulation completed in 56.28 s
Total time steps solved: 201
