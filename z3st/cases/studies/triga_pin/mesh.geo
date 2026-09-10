// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 2D axisymmetric (r-z) TRIGA fuel element, 103/104 series
//
//  Geometry from the OpenMC model of the Pavia TRIGA Mark II
//  (MP-OFELIA, Tutorials/openmc/TRIGA): central Zr rod 0.285 cm,
//  U-ZrH meat to 1.82 cm, SS-304 cladding to 1.88 cm, active
//  height 8.81 -> 46.91 cm.
//
//  The central Zr rod is not meshed: it is a hole with an adiabatic
//  inner boundary. [TBC] whether the rod's conduction path matters.
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

r_1_i = 0.00285;    // Zr rod radius = fuel inner radius (m)
r_1_o = 0.01820;    // fuel outer radius (m)
r_2_i = 0.01820;    // clad inner radius (m)
r_2_o = 0.01880;    // clad outer radius (m)
h     = 0.38100;    // active fuel height (m)

n_r1 = 21;          // radial divisions, fuel
n_r2 = 5;           // radial divisions, clad
n_z  = 41;          // axial divisions (shared)

// --. fuel : rectangle [r_1_i, r_1_o] x [0, h] --..
Point(1) = {r_1_i, 0, 0};
Point(2) = {r_1_o, 0, 0};
Point(3) = {r_1_o, h, 0};
Point(4) = {r_1_i, h, 0};
Line(1) = {1, 2};   // bottom_1  (z = 0)
Line(2) = {2, 3};   // lateral_1 (r = r_1_o)
Line(3) = {3, 4};   // top_1     (z = h)
Line(4) = {4, 1};   // inner_1   (r = r_1_i, Zr rod surface)
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --. clad : rectangle [r_2_i, r_2_o] x [0, h] --..
Point(5) = {r_2_i, 0, 0};
Point(6) = {r_2_o, 0, 0};
Point(7) = {r_2_o, h, 0};
Point(8) = {r_2_i, h, 0};
Line(5) = {5, 6};   // bottom_2 (z = 0)
Line(6) = {6, 7};   // outer_2  (r = r_2_o)
Line(7) = {7, 8};   // top_2    (z = h)
Line(8) = {8, 5};   // inner_2  (r = r_2_i)
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

// --- structured quad mesh ---
Transfinite Line {1, 3} = n_r1;
Transfinite Line {5, 7} = n_r2;
Transfinite Line {2, 4, 6, 8} = n_z;
Transfinite Surface {1};
Transfinite Surface {2};
Recombine Surface {1, 2};

// --- physical groups (ids must match geometry.yaml) ---
Physical Surface("fuel", 1) = {1};
Physical Surface("clad", 2) = {2};

Physical Curve("bottom_1", 1)  = {1};
Physical Curve("bottom_2", 2)  = {5};
Physical Curve("lateral_1", 3) = {2};
Physical Curve("outer_2", 4)   = {6};
Physical Curve("inner_2", 5)   = {8};
Physical Curve("top_2", 6)     = {7};
Physical Curve("top_1", 7)     = {3};
Physical Curve("inner_1", 8)   = {4};

Mesh.ElementOrder = 1;
