// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  2D periodic row of lenticular cavities along a grain boundary (y=0)
//  Unit cell = quarter-model: one half-lentille centered at x=0 (cut
//  by the xmin symmetry plane) and one half-lentille centered at
//  x=Lx (cut by the xmax symmetry plane). Lx IS the center-to-center
//  directory: this is NOT two isolated bubbles -- xmax is a mirror
//  plane, not a free surface.
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

// --- Parameters (overridable via -setnumber) ---
If(!Exists(Rp))
  Rp = 11.2e-6;        // Projected radius of the cavity (m)
EndIf
If(!Exists(theta_deg))
  theta_deg = 48.1;    // Semi-dihedral angle (deg)
EndIf
If(!Exists(Lx))
  Lx = 4.500000e-05;          // center-to-center bubble spacing (m) = quarter-cell width
EndIf
If(!Exists(Ly))
  Ly = 60e-6;
EndIf
If(!Exists(h_cavity))
  h_cavity = 0.15e-6;  // retained from the angle sweep; refinement proportional to defect size, h <= lc/4 = 0.5e-6 (see README correction note)
EndIf
If(!Exists(h_plate))
  h_plate = 4.0e-6;    // coarse mesh size
EndIf

// Rectangular plate (Quart de symetrie : coin inferieur gauche en 0,0)
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

ay = Rp * Tan(theta_deg * Pi / 360.0);
yc = (ay*ay - Rp*Rp) / (2*ay);

// Left half-bubble, centered at x=0 (cut by xmin)
xc = 0;
Call MakeLenticularbubble;
lens_left = lens_surface;

// Right half-bubble, centered at x=Lx (cut by xmax)
xc = Lx;
Call MakeLenticularbubble;
lens_right = lens_surface;

// Subtract both half-lentilles from the rectangle in one boolean call
BooleanDifference{ Surface{1}; Delete; }{ Surface{lens_left}; Surface{lens_right}; Delete; }

//---------------------------------------------------------
// --- Groupes Physiques ---
Physical Surface("uo2") = {1};

// Sélection dynamique des frontières par boîte englobante pour éviter les erreurs d'IDs post-opération booléenne
// eps must stay well above OpenCASCADE's built-in confusion tolerance
// (Precision::Confusion() = 1e-7); it also sets the minimum viable
// ligament: the ligament (Lx - 2*Rp) must exceed ~2*eps = 1e-6 m for
// the cavity/ligament bounding boxes below to stay cleanly separated.
eps = 0.5e-6;

// Ligament: from the left bubble's tip (x=Rp) to the right bubble's
// tip (x=Lx-Rp) -- NOT to the domain edge anymore.
c_ymin() = Curve In BoundingBox{Rp-eps, -eps, -eps, Lx-Rp+eps, eps, eps};

// xmin and xmax are now BOTH symmetry planes cut by a half-bubble:
// each spans from the local tip height y=ay up to the top y=Ly.
c_xmin() = Curve In BoundingBox{-eps, ay-eps, -eps, eps, Ly+eps, eps};
c_xmax() = Curve In BoundingBox{Lx-eps, ay-eps, -eps, Lx+eps, Ly+eps, eps};

c_ymax() = Curve In BoundingBox{-eps, Ly-eps, -eps, Lx+eps, Ly+eps, eps};

// Cavity = union of BOTH half-lentille arcs (kept as a single physical
// group so a uniform pressure BC applies to both at once).
c_cavity_left()  = Curve In BoundingBox{-eps, -eps, -eps, Rp+eps, ay+eps, eps};
c_cavity_right() = Curve In BoundingBox{Lx-Rp-eps, -eps, -eps, Lx+eps, ay+eps, eps};
c_cavity() = {c_cavity_left(), c_cavity_right()};

Printf("Boundary group ymin: %g curve(s) found", #c_ymin());
Printf("Boundary group xmin: %g curve(s) found", #c_xmin());
Printf("Boundary group xmax: %g curve(s) found", #c_xmax());
Printf("Boundary group ymax: %g curve(s) found", #c_ymax());
Printf("Boundary group cavity_left: %g curve(s) found", #c_cavity_left());
Printf("Boundary group cavity_right: %g curve(s) found", #c_cavity_right());
Printf("Boundary group cavity (total): %g curve(s) found", #c_cavity());
Printf("Ligament length (Lx - 2*Rp) = %g m", Lx - 2*Rp);

Physical Curve("ymin") = {c_ymin()};
Physical Curve("xmin") = {c_xmin()};
Physical Curve("xmax") = {c_xmax()};
Physical Curve("ymax") = {c_ymax()};
Physical Curve("cavity") = {c_cavity()};

// Mesh Refinement -- the Distance field now covers the cavity arcs AND
// the full ligament segment between the two tips (x=Rp to x=Lx-Rp),
// not the old single-bubble "tip to far domain edge" span.
Field[1] = Distance;
Field[1].CurvesList = {c_cavity(), c_ymin()};
Field[1].NumPointsPerCurve = 400;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_cavity;
Field[2].SizeMax = h_plate;
Field[2].DistMin = 0.5e-6;   // Taille minimale jusqu'à 0.5 µm des courbes
Field[2].DistMax = 10.0e-6;  // Transition progressive vers la taille grossière sur 10 µm

Background Field = 2;

// Mesh generation options
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6; // Frontal-Delaunay for better quality in 2D
