# Define a polytope somehow

DeclareGlobalVariable( "Simplex1", "line segment - simplex of dimension 1" );

DeclareGlobalVariable( "Simplex2", "triangle - simplex of dimension 2" );

DeclareGlobalVariable( "Simplex3", "tetrahedron - simplex of dimension 2" );

DeclareGlobalVariable( "3sq", "triple square" );

DeclareGlobalVariable( "T2", "2-dimensional torus, 4 symmetries: 2 translations, x -> -x in both coordinates together, rotation" );

DeclareGlobalVariable( "Cube", "cube" );

DeclareGlobalVariable( "S3cubic", "sphere S^3 as two cubes" );

DeclareGlobalVariable( "Sphere1", "circle, with two vertices and symmetries" );

DeclareGlobalVariable( "Sphere2", "sphere S^2 made of two triangles" );

DeclareGlobalVariable( "Sphere4", "sphere S^4 made of two 4-simplices" );

DeclareGlobalVariable( "Disk3", "disk D^3 with two triangles as boundary and one vertex, namely A, inside" );

DeclareGlobalVariable( "3cubes", "non-manifold: 3 cubes all glued at the same boundary" );

DeclareGlobalVariable( "Bigon", "bigon" );

DeclareGlobalVariable( "3pillow", "3-d pillow" );

DeclareGlobalVariable( "3pillow1", "3-d pillow - with another order on vertices" );

DeclareGlobalVariable( "Lantern", "lantern - D^3 with 3 faces" );

DeclareGlobalVariable( "DoubleLantern", "double lantern (interior and exterior) with 4 edges" );

DeclareGlobalVariable( "MobiusBand", "Mobius band made of two squares" );

DeclareGlobalVariable( "PoincareSphere", "Poincare sphere");
DeclareGlobalVariable( "s2s2twisted", "Twisted product of 2-sphere");
DeclareGlobalVariable( "cp2", "Complex projective plane" );

############################################################# now functions

DeclareGlobalFunction( "Lens" ); # 3-d lens spaces; usage: Lens(7,1)

DeclareGlobalFunction( "PolPrint" ); # print all the faces of simplitial complex in terms of names of vertices

DeclareGlobalFunction( "ballAB" ); # a ball of dimension n and just two vertices A and B

DeclareGlobalFunction( "sphereAB" ); # a sphere of dimension n and just two vertices A and B

DeclareGlobalFunction( "KummerSurface" ); # KummerSurface() makes the Kummer surface

DeclareGlobalFunction( "TorTwist" ); # Вычисляется n раз скученный тор по одной матрице
DeclareGlobalFunction( "projectivePlane" ); # Вещественные проективные пространства

