// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for the cohesive strength-surface test
//  (Vicentini et al. 2026, Fig. 21)
//
//  Square block loaded by opposing displacements on opposing edges, so the
//  stress state stays homogeneous until nucleation. The element count is odd
//  so that no node sits on the centre lines.
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

L  = 1.0e-3;    // edge length (m)
nx = 41;        // nodes per side -> 40 cells, h = ell

Point(1) = {-L/2, -L/2, 0, 1.0};
Point(2) = { L/2, -L/2, 0, 1.0};
Point(3) = { L/2,  L/2, 0, 1.0};
Point(4) = {-L/2,  L/2, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Structured triangles (no Recombine).
Transfinite Curve {1, 2, 3, 4} = nx Using Progression 1;
Transfinite Surface {1};

Physical Curve("bottom") = {1};
Physical Curve("right")  = {2};
Physical Curve("top")    = {3};
Physical Curve("left")   = {4};
Physical Surface("block") = {1};
