DeclareGlobalFunction( "PolBnd" );
# Creates an index of boundary faces of face [d,fn] of complex s
# <result>[i] --- index of (i-1)-dimensional faces of s which are in the boundary of [d,fn]
# Input data: polytope, d, fn

##  DeclareGlobalFunction( "PolFaceVertices" );
##  # function returning set of vertices (as numbers) bounding given face of complex
##  # s for complex, d for dimension of the face, f for number of the face

DeclareGlobalFunction( "PolCheckComb" );
# check if a face of given complex is combinatorial complex
# that is, whether every its subface of lower dimension is uniquely
# determined by its vertices (and dimension)
# input: complex p, face dimension d and face number fn

DeclareGlobalFunction( "PolProduct" );
# Cartesian product of two polytopes

DeclareGlobalFunction( "PolProductSyms" );
# Cartesian product of two polytopes with symmetries of multipliers transferred to it.
# First go the symmetries of the first multiplier, then - the second

DeclareGlobalFunction( "PolProductSymsDict" );
# Cartesian product of two polytopes with symmetries of multipliers transferred to it.
# First go the symmetries of the first multiplier, then - the second.
# Also returns the face dictionary.

DeclareGlobalFunction( "PolTriangulate" );
# triangulating a polytope ( = ball complex)

DeclareGlobalFunction( "OrientTriangulated" );
# Function computes consistent orientation on simplices of greatest dimension
# on a given triangulated complex s.
# Returns array of -1,1-s which correspond to the orientation of simplices 
# of greatest dimension of s

DeclareGlobalFunction( "PolInnerFaces" );
# build index of inner faces of given polytope complex
# (!) <returned value>[i] --- set of inner faces of dimension (i-1)
# Any face is outer if it has at most 1 adjacent face of higher dimension
# or if it lies in the boundary of such a face.
# And inner faces are not outer faces.

DeclareGlobalFunction( "PolDoubleCone" );
# Make a double cone with vertices V1 and V2 over the given polytope p

DeclareGlobalFunction( "PolCylinder" );
# Maybe not needed in view of the existence of PolProduct(Syms), but anyhow:
# triangulates a polytope M and creates a triangulated cylinder over it;
# in the list for faces of any dimensions (including vertices),
# first go two identical copies of triangulated M -   M x {0}  and  M x {1}
# Also, <result>.l[i+1] is the number of i-faces in the triangulated M

DeclareGlobalFunction( "MaxTree" );
# finds a maximal tree in the 1-skeleton of a polytope as a list of edges

DeclareGlobalFunction( "CellOrient" );
# Argument: ball complex  p .
# provides inductively some orientations for cells of dimensions 1..(dim p) 
# as consistent orientations of their faces.

DeclareGlobalFunction( "PolOrient" );
# Argument: ball complex  p of dimension >1.
# If  p  is orientable, gives a consistent orientation of n-faces (n= dim p),
# otherwise returns fail

DeclareGlobalFunction( "Slovo" );
# Arguments: two lists, both corresponding to a certain 2-face.
# The first list consists of its edges represented in the standard way as sets of their vertices.
# The second list consists of group elements ascribed to these edges.
# Output: a word ("slovo" in Russian) corresponding to the 2-face

DeclareGlobalFunction( "FundGroup" );
# Computes the fundamental group of the given polytope

DeclareGlobalFunction( "PolMinusFace" );
# cuts out a neighborhood of a face with given address from polytope
# Arguments: p, [d, fn]

##  DeclareGlobalFunction( "PolMinusPol" );
# cuts out a neighborhood of a subpolytope from polytope
# Arguments: p - polytope, sp=rec(vertices:=[...],faces:=[...]) - addresses of cut out faces
# sp.faces[i] is the list of numbers of i-faces to be removed.
# output: new polytope новый политоп
# dependencies: PolMinusFace.g

DeclareGlobalFunction( "DelFace" );
# an auxiliary function: cuts out a face with address [k, m] from polytope  p,
# shifting back all necessary numbers. This returns a polytope-like structure which
# may not be ball complex, but it works well if employed skillfully.
# Arguments: p, address

DeclareGlobalFunction( "PolSimplify1" );
# symplifies a polytope by merging two k-faces when possible, starting from k=n
# and down to k=1. If actual simplification is achieved, try this again and so on.
# Probably, there will be "PolSymplify2" and so on in the future.

DeclareGlobalFunction( "PolFactorInvolution" );
# PolFactorInvolution := function( p, invol )
# p is the polytope with symmetries
# invol is such a list of some of its symmetries (repetitions possible)
# that it is known that the product of symmetries in s is an involution
# returns the factored polytope


