DeclareGlobalFunction( "Matrices4dAffineTriangulated" );
# Returns a structure with matrices f2,f3_full,f4_full,f5 and supplementary information
# for given polytope complex of dimension 4. Arguments:
# p is the initial polytope (not necessarily triangulated),
# num is a logical (I mean Boolean) variable: if num is true, the coordinates z[i] are numbers,
# f is the function determining the coordinates: z[i]=f(i). So, if num is false f does not affect the result.
# (!) Does not work if some matrices are empty

DeclareGlobalFunction( "f3E4" );
# Arguments: 4-simplex  four_simplex = [i1,i2,i3,i4,i5] ,  its 2-face  two_face  and edge  edge .
# Computes the (partial) derivative of the squared sine of dihedral angle at two_face wrt squared edge length

DeclareGlobalFunction( "Matrices4dEuclid4" );
# return structure with matrices f2_full,f3_FULL,f4_full,f5
# constructed using Euclidean 4-dimensional geometry
# and supplementary information
# for given 4-manifold with one-component boundary. Arguments:
# p is the initial polytope (not necessarily triangulated),
# coord is the list [w,x,y,z], where w,x,y,z are long enough lists of rational numbers
# and such that no five vertices lie in a 3-plane




