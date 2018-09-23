InstallGlobalFunction( GrassmannSum, 

function( x,y )

# x and y are records representing elements of a Grassmann algebra
# RecNames are "monomials" and "coeffs"
# "monomials" is a list of monomials, each element being a product of Grassmann generators
# if this element is [3,4], then it is (generator no. 3)*(generator no. 4)
# and these elements *must* be sets!
# "coeffs" is the list of coefficients at corresponding products

local i,j, z;

z := StructuralCopy(x);

for i in [1..Length(y.monomials)] do
  if y.monomials[i] in z.monomials then
    j := Position(z.monomials, y.monomials[i]);
    z.coeffs[j] := z.coeffs[j] + y.coeffs[i];
  else Add( z.monomials, y.monomials[i] );
    Add( z.coeffs, y.coeffs[i] );
  fi;
od;

return z;

end );

########################################################################################


InstallGlobalFunction( GrassmannMonomialsProduct,

function( x1,x2 )

# x1 = [<monomial>,<coefficient>], x2 similarly
# returns their product written in a similar way

local m1,m2,x,c,u,s;

m1 := x1[1];
m2 := x2[1];

if Intersection(m1,m2) <> [] then
  x := [ [], 0 ];
else c := Concatenation( m1,m2 );
  u := UnionSet( m1,m2 );
  s := SignPerm( PermListList( c,u ) );
  x := [ u, s*x1[2]*x2[2] ];
fi;

return(x);

end );

########################################################################################

InstallGlobalFunction( GrassmannProduct,

function( x,y )

# x and y are records representing elements of Grassmann algebra
# returns x*y

local z,i,j,x1,y1,z1;

z := rec( monomials := [], coeffs := [] );

for i in [1..Length(x.monomials)] do
  x1 := [ x.monomials[i], x.coeffs[i] ];
  for j in [1..Length(y.monomials)] do
    y1 := [ y.monomials[j], y.coeffs[j] ];
    z1 := GrassmannMonomialsProduct( x1,y1 );
    z := GrassmannSum( z, rec( monomials:=[z1[1]], coeffs:=[z1[2]] ) );
  od;
od;

return z;

end );

########################################################################################


InstallGlobalFunction( BerezinIntegral,

function( x, a )

# x is a record representing Grassmann algebra element
# a is a (number of) Grassmann variable
# returns the record representing the integral

local i,int,l,m,s;

int:= rec(monomials:=[],coeffs:=[]);

for i in [1..Length(x.monomials)] do
  if a in x.monomials[i] then
    l := Difference( x.monomials[i], [a] );
    m := Concatenation( l, [a] );
    s := SignPerm( PermListList( x.monomials[i], m ) );
    Add( int.monomials, l );
    Add( int.coeffs, s*x.coeffs[i] );
  fi;
od;

return(int);

end );

########################################################################################


InstallGlobalFunction( BerezinMultipleIntegral,

function( x, l )

# x is a record representing Grassmann algebra element
# l is a list of (numbers of) Grassmann variables
# integration goes first in the first variable, then the second etc.
# returns the record representing the integral

local i,b;

b := StructuralCopy( x );

for i in [1..Length(l)] do
  b := BerezinIntegral( b,l[i] );
od;  

return(b);

end );

########################################################################################



