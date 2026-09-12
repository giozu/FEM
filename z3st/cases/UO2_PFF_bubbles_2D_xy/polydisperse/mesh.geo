// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  2D periodic row of POLYDISPERSE lenticular cavities along a grain
//  boundary (y=0): bubble A (Rp_1) and bubble B (Rp_2) alternate,
//  Lx apart center-to-center, EVERY consecutive pair (A-then-B or
//  B-then-A) at the same spacing Lx. Same quarter-model structure as
//  periodic_row/mesh.geo, just with two different Rp instead of one.
//
//  SYMMETRY (corrected from an earlier, wrong rejection of this
//  approach): an alternating A-B-A-B row with UNIFORM center-to-center
//  spacing Lx *does* have a mirror plane through every bubble center,
//  regardless of size. Proof sketch: unfolding the cell below across
//  xmin (bubble A's own center) reproduces a full A; unfolding across
//  xmax (bubble B's own center) reproduces a full B; repeating this
//  tiling gives A(0), B(Lx), A(2Lx), B(3Lx), ... -- uniform spacing Lx
//  everywhere, hence a SINGLE ligament value throughout the infinite
//  row, not two. (What is NOT valid is a construction with two
//  DIFFERENT gaps around each bubble -- that would break the
//  per-bubble mirror symmetry; this is not that.)
//
//  Lx = A-center-to-B-center distance. Ligament = Lx - Rp_1 - Rp_2.
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

// --- Parameters (overridable via -setnumber) ---
If(!Exists(Rp_1))
  Rp_1 = 11.2e-6;      // bubble A projected radius (m), centered on xmin
EndIf
If(!Exists(Rp_2))
  Rp_2 = 5.6e-6;       // bubble B projected radius (m), centered on xmax
EndIf
If(!Exists(theta_deg))
  theta_deg = 48.1;    // semi-dihedral angle (deg), same for both bubbles
EndIf
If(!Exists(Lx))
  Lx = 39.4e-6;        // center(A)-to-center(B) spacing (m); gives a
                        // 22.6e-6 ligament, matching periodic_row/Lx_45
                        // (monodisperse) for direct comparison
EndIf
If(!Exists(Ly))
  Ly = 60e-6;
EndIf
If(!Exists(h_cavity))
  h_cavity = 0.15e-6;  // unchanged calibration. No per-bubble curvature
                        // criterion applies (the tip is a wedge corner,
                        // not an ellipse -- see README correction note);
                        // the only real constraint is h <= lc/4 = 0.5e-6,
                        // satisfied here for both Rp_1 and Rp_2.
EndIf
If(!Exists(h_plate))
  h_plate = 4.0e-6;
EndIf

// Rectangular plate (quart-de-symetrie: one bubble center on each of
// xmin, xmax; ymin is the grain-boundary symmetry plane)
Rectangle(1) = {0, 0, 0, Lx, Ly};

//----------Lenticular shape (same macro as the single-bubble case)------------------
Macro MakeLenticularbubble
  p1 = newp; Point(p1) = {xc - Rp, 0, 0};
  p2 = newp; Point(p2) = {xc + Rp, 0, 0};
  p5 = newp; Point(p5) = {xc, yc, 0};
  p6 = newp; Point(p6) = {xc, -yc, 0};

  c1 = newc; Circle(c1) = {p1, p5, p2};
  c2 = newc; Circle(c2) = {p2, p6, p1};

  ll = newll; Curve Loop(ll) = {c1, c2};
  s_ = news;  Plane Surface(s_) = {ll};

  lens_surface = s_;
Return

// Bubble A (Rp_1), half-lentille centered at x=0 (cut by xmin)
Rp = Rp_1;
ay_1 = Rp * Tan(theta_deg * Pi / 360.0);
yc = (ay_1*ay_1 - Rp*Rp) / (2*ay_1);
xc = 0;
Call MakeLenticularbubble;
lens_A = lens_surface;

// Bubble B (Rp_2), half-lentille centered at x=Lx (cut by xmax)
Rp = Rp_2;
ay_2 = Rp * Tan(theta_deg * Pi / 360.0);
yc = (ay_2*ay_2 - Rp*Rp) / (2*ay_2);
xc = Lx;
Call MakeLenticularbubble;
lens_B = lens_surface;

// Subtract both half-lentilles from the rectangle in one boolean call
BooleanDifference{ Surface{1}; Delete; }{ Surface{lens_A}; Surface{lens_B}; Delete; }

//---------------------------------------------------------
// --- Groupes Physiques ---
Physical Surface("uo2") = {1};

// eps must stay well above OpenCASCADE's built-in confusion tolerance
// (Precision::Confusion() = 1e-7); it also sets the minimum viable
// ligament: Lx - Rp_1 - Rp_2 must exceed ~2*eps = 1e-6 m for the
// cavity/ligament bounding boxes below to stay cleanly separated.
eps = 0.5e-6;

// Single ligament: from bubble A's tip (x=Rp_1) to bubble B's tip
// (x=Lx-Rp_2). Uniform spacing Lx means this is the ONLY ligament
// value in the infinite row (see symmetry note above).
c_ymin() = Curve In BoundingBox{Rp_1-eps, -eps, -eps, Lx-Rp_2+eps, eps, eps};

// xmin (bubble A's own symmetry plane) spans y=ay_1 to Ly;
// xmax (bubble B's own symmetry plane) spans y=ay_2 to Ly -- DIFFERENT
// height than xmin since ay_2 != ay_1 (Rp_2 != Rp_1).
c_xmin() = Curve In BoundingBox{-eps, ay_1-eps, -eps, eps, Ly+eps, eps};
c_xmax() = Curve In BoundingBox{Lx-eps, ay_2-eps, -eps, Lx+eps, Ly+eps, eps};

c_ymax() = Curve In BoundingBox{-eps, Ly-eps, -eps, Lx+eps, Ly+eps, eps};

// cavity_A and cavity_B are kept as SEPARATE physical groups (so
// different pressures can be applied later if wanted). NOTE: no
// combined "cavity" group is created on top of these -- dolfinx's
// mesh reader (dolfinx.io.gmsh.read_from_msh, what z3st actually
// calls) raises "All cells are expected to be tagged once, found
// duplicates" if a curve belongs to more than one physical group of
// the same dimension (verified directly). A uniform pressure on both
// bubbles must instead be expressed as two boundary_conditions.yaml
// entries (region: cavity_A / region: cavity_B) with identical
// traction values -- functionally equivalent, no meshing trick needed.
c_cavity_A() = Curve In BoundingBox{-eps, -eps, -eps, Rp_1+eps, ay_1+eps, eps};
c_cavity_B() = Curve In BoundingBox{Lx-Rp_2-eps, -eps, -eps, Lx+eps, ay_2+eps, eps};

Printf("Boundary group ymin: %g curve(s) found", #c_ymin());
Printf("Boundary group xmin: %g curve(s) found", #c_xmin());
Printf("Boundary group xmax: %g curve(s) found", #c_xmax());
Printf("Boundary group ymax: %g curve(s) found", #c_ymax());
Printf("Boundary group cavity_A (bubble A, Rp_1): %g curve(s) found", #c_cavity_A());
Printf("Boundary group cavity_B (bubble B, Rp_2): %g curve(s) found", #c_cavity_B());
Printf("Ligament length (Lx - Rp_1 - Rp_2) = %g m", Lx - Rp_1 - Rp_2);

Physical Curve("ymin") = {c_ymin()};
Physical Curve("xmin") = {c_xmin()};
Physical Curve("xmax") = {c_xmax()};
Physical Curve("ymax") = {c_ymax()};
Physical Curve("cavity_A") = {c_cavity_A()};
Physical Curve("cavity_B") = {c_cavity_B()};

// Mesh Refinement
Field[1] = Distance;
Field[1].CurvesList = {c_cavity_A(), c_cavity_B(), c_ymin()};
Field[1].NumPointsPerCurve = 400;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_cavity;
Field[2].SizeMax = h_plate;
Field[2].DistMin = 0.5e-6;
Field[2].DistMax = 10.0e-6;

Background Field = 2;

// Mesh generation options
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6;
