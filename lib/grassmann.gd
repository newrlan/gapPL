DeclareGlobalFunction( "GrassmannSum" );
# GrassmannSum := function( x,y )
# x and y are records representing elements of a Grassmann algebra
# RecNames are "monomials" and "coeffs"
# "monomials" is a list of monomials, each element being a product of Grassmann generators
# if this element is [3,4], then it is (generator no. 3)*(generator no. 4)
# and these elements *must* be sets!
# "coeffs" is the list of coefficients at corresponding products

DeclareGlobalFunction( "GrassmannMonomialsProduct" );
# GrassmannMonomialsProduct := function( x1,x2 )
# x1 = [<monomial>,<coefficient>], x2 similarly
# returns their product written in a similar way

DeclareGlobalFunction( "GrassmannProduct" );
# GrassmannProduct := function( x,y )
# x and y are records representing elements of Grassmann algebra
# returns x*y

DeclareGlobalFunction( "BerezinIntegral" );
# BerezinIntegral := function( x, a )
# x is a record representing Grassmann algebra element
# a is a (number of) Grassmann variable
# returns the record representing the integral

DeclareGlobalFunction( "BerezinMultipleIntegral" );
# BerezinMultipleIntegral := function( x, l )
# x is a record representing Grassmann algebra element
# l is a list of (numbers of) Grassmann variables
# integration goes first in the first variable, then the second etc.
# returns the record representing the integral

