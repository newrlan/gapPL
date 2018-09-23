InstallGlobalFunction( Matrices4dAffineTriangulated,

# return structure with matrices f2,f3_full,f4_full,f5
# and supplementary information
# for given polytope complex of dimension 4

# Does not work if some matrices are empty

function(p, num, f)

# p is the initial polytope (not necessarily triangulated)
# num is a logical (I mean Boolean) variable:
      # if num is true, the z[i] below are numbers
# f is the function determining these numbers: z[i]=f(i)
      # so, if num is false f does not affect the result

  local s         # triangulated p
    ,F            # structure to be returned
    ,inner_faces  # list of inner faces (inner_faces[i] - set of 
                  # inner faces of dimension (i-1)
    ,ind_f2       # index of all 2-faces, inner faces first
    ,ind_f2_li    # position of last inner 2-face in index of all 2-faces
    ,ind_f3       # index of all 3-faces, inner faces first
    ,ind_f3_li    # position of last inner 3-face in index of all 3-faces
    ,faces_ind    # backward index of the faces
                  # faces_ind[d+1][i] - index of face of dimension d 
                  # in the list of "our" faces of corresponding dimension
    ,f2           # 
    ,f3_full      #
    ,f4_full      #
    ,f5           # matrices that we are looking for
    ,i,j,k,l      # counters
    ,nl           # 
    ,l_v,j_v      # and local stuff
    ,z            # array of zetas 
    ,s4_orient    # orientations of 4-simplices
    ,xxx          # temporary Rational indeterminate
    ;

 # First triangulate, if needed
  s := PolTriangulate(p);

  # Get index of inner faces
  inner_faces := PolInnerFaces(s);

  # Enumeration of inner vertices is done automatically in inner_faces.
  # The same is true for inner 2-faces and inner 3-faces

  # Enumerate all 2-faces (inner 2-faces first)
  ind_f2 := StructuralCopy(inner_faces[3]);
  ind_f2_li := Length(inner_faces[3]);
  Append(ind_f2, Difference([1..Length(s.faces[2])], ind_f2));

  # Enumerate all 3-faces (inner 3-faces first)
  ind_f3 := StructuralCopy(inner_faces[4]);
  ind_f3_li := Length(inner_faces[4]);
  Append(ind_f3, Difference([1..Length(s.faces[3])],ind_f3));

  # build backward index
  faces_ind := List([1..5], i->[]);
  faces_ind[1] := List([1..Length(s.vertices)], i->0);
  for i in [2..5] do
    faces_ind[i] := List([1..Length(s.faces[i-1])], j->0);
  od;
  for i in [1,2] do
    for j in [1..Length(inner_faces[i])] do
      faces_ind[i][inner_faces[i][j]] := j;
    od;
  od;
  for i in [1..Length(ind_f2)] do
    faces_ind[3][ind_f2[i]] := i;
  od;
  for i in [1..Length(ind_f3)] do
    faces_ind[4][ind_f3[i]] := i;
  od;

  # assign xxx
  if num 
    then xxx := 1;
    else xxx := Indeterminate(Rationals, "xxx");
  fi;
  # init matrices
  f2 := NullMat(Length(ind_f2), Length(inner_faces[1]));
  f3_full := NullMat(2*Length(ind_f3), Length(ind_f2));
  f4_full := NullMat(3*Length(s.faces[4]),  2*Length(ind_f3));
  f5 := NullMat(Length(inner_faces[1]), 3*Length(s.faces[4]));

  # Now: fill out elements of our matrices that are non-zero

  # First of all, prepare the zetas
  if num
    then z := List([1..Length(s.vertices)], f);
    else z := List([1..Length(s.vertices)], 
      i->Indeterminate(Rationals,Concatenation("z",String(i))));
  fi;

  # f2: 
  nl := DimensionsMat(f2)[1];
  # loop over rows
  for l in [1..nl] do
    # get sorted list of 3 vertices bounding l-th (not necessary) inner 2-face
    l_v := PolBnd(s,2,ind_f2[l])[1];
    Sort(l_v);
    i := faces_ind[1][l_v[1]];
    j := faces_ind[1][l_v[2]];
    k := faces_ind[1][l_v[3]];
    # fill out some elements of the row f2[l]
    if i<>0 then
      f2[l][i] := 1/(z[l_v[1]]-z[l_v[2]]) - 1/(z[l_v[1]]-z[l_v[3]]);
    fi;
    if j<>0 then
      f2[l][j] := - 1/(z[l_v[1]]-z[l_v[2]]);
    fi;
    if k<>0 then
      f2[l][k] := 1/(z[l_v[1]]-z[l_v[3]]);
    fi;
  od;

  # f3_full:
  nl := Length(ind_f3);
  # loop over all 3-faces
  for l in [1..nl] do
    # get sorted list of 4 vertices bounding l-th 3-face
    l_v := PolBnd(s,3,ind_f3[l])[1];
    Sort(l_v);
    # loop over all (not only inner) 2-faces which bound this 3-face
    for k in s.faces[3][ind_f3[l]] do
      j := faces_ind[3][k];
      # get sorted vertices of 2-face
      j_v := PolBnd(s,2,ind_f2[j])[1];
      Sort(j_v);
      # the conditional statement below seems to be superfluous
      if IsSubset(l_v,j_v) then
        # Fill out the respective f3 elements
        if j_v = l_v{[1,2,3]} then
          f3_full[2*l-1][j] := 1;
          f3_full[2*l  ][j] :=   (z[l_v[3]]-z[l_v[1]])/(z[l_v[2]]-z[l_v[3]]);
        elif j_v = l_v{[1,2,4]} then
          f3_full[2*l-1][j] := -1;
          f3_full[2*l  ][j] := - (z[l_v[4]]-z[l_v[1]])/(z[l_v[2]]-z[l_v[4]]);
        elif j_v = l_v{[1,3,4]} then
          f3_full[2*l-1][j] := 1;
        elif j_v = l_v{[2,3,4]} then
          f3_full[2*l  ][j] := -1;
        fi;
      fi;
    od;
  od;

  # f4_full
  # loop over all 4-faces
  for l in [1..Length(s.faces[4])] do
    # get sorted list of vertices
    l_v := PolBnd(s,4,l)[1];
    Sort(l_v);
    # loop over 3-faces bounding this 4-face
    for j in s.faces[4][l] do
      # get sorted list of vertices of 3-face
      j_v := PolBnd(s,3,j)[1];
      Sort(j_v);
      k := faces_ind[4][j];
      # Fill out respective f4_full elements
      if j_v = l_v{[1,2,3,4]} then
        f4_full[3*l-2][2*k-1] := 1;
        f4_full[3*l-1][2*k-1] := 0;
        f4_full[3*l  ][2*k-1] := - (z[l_v[1]]-z[l_v[4]])/(z[l_v[3]]-z[l_v[4]]);
        f4_full[3*l-2][2*k  ] := 0;
        f4_full[3*l-1][2*k  ] := 1;
        f4_full[3*l  ][2*k  ] := - (z[l_v[2]]-z[l_v[4]])/(z[l_v[3]]-z[l_v[4]]);
      elif j_v = l_v{[1,2,3,5]} then
        f4_full[3*l-2][2*k-1] := - 1;
        f4_full[3*l-1][2*k-1] := 0;
        f4_full[3*l  ][2*k-1] :=   (z[l_v[1]]-z[l_v[5]])/(z[l_v[3]]-z[l_v[5]]);
        f4_full[3*l-2][2*k  ] := 0;
        f4_full[3*l-1][2*k  ] := - 1;
        f4_full[3*l  ][2*k  ] :=   (z[l_v[2]]-z[l_v[5]])/(z[l_v[3]]-z[l_v[5]]);
      elif j_v = l_v{[1,2,4,5]} then
        f4_full[3*l-2][2*k-1] := 1;
        f4_full[3*l-1][2*k-1] := 0;
        f4_full[3*l  ][2*k-1] := 0;
        f4_full[3*l-2][2*k  ] := 0;
        f4_full[3*l-1][2*k  ] := 1;
        f4_full[3*l  ][2*k  ] := 0;
      elif j_v = l_v{[1,3,4,5]} then
        f4_full[3*l-2][2*k-1] := - 1;
        f4_full[3*l-1][2*k-1] := 0;
        f4_full[3*l  ][2*k-1] := 0;
        f4_full[3*l-2][2*k  ] := 0;
        f4_full[3*l-1][2*k  ] := 0;
        f4_full[3*l  ][2*k  ] := - 1;
      elif j_v = l_v{[2,3,4,5]} then
        f4_full[3*l-2][2*k-1] := 0;
        f4_full[3*l-1][2*k-1] := 1;
        f4_full[3*l  ][2*k-1] := 0;
        f4_full[3*l-2][2*k  ] := 0;
        f4_full[3*l-1][2*k  ] := 0;
        f4_full[3*l  ][2*k  ] := 1;
      fi;
   od;
  od;

  # f5
  # make consistent orientation on 4-simplices
  s4_orient := OrientTriangulated(s);
  if s4_orient=[] then
    Print("Non-orientable manifold, aborting!\n");
  fi;
  # loop over inner vertices
  for i in [1..Length(inner_faces[1])] do
    # loop over all 4-simplices having inner_faces[1][i] as vertex
    for j in [1..Length(s.faces[4])] do
      j_v := FaceComp(s,[4,j]).0;
      if inner_faces[1][i] in j_v then
        k := Position(j_v, inner_faces[1][i]);
        if k=1 then
          f5[i][3*j-2] := s4_orient[j];
        elif k=2 then
          f5[i][3*j-1] := s4_orient[j];
        elif k=3 then
          f5[i][3*j  ] := s4_orient[j];
        elif k=4 then
          f5[i][3*j-2] := -s4_orient[j] * (z[j_v[1]]-z[j_v[5]])/(z[j_v[4]]-z[j_v[5]]);
          f5[i][3*j-1] := -s4_orient[j] * (z[j_v[2]]-z[j_v[5]])/(z[j_v[4]]-z[j_v[5]]);
          f5[i][3*j  ] := -s4_orient[j] * (z[j_v[3]]-z[j_v[5]])/(z[j_v[4]]-z[j_v[5]]);
        elif k=5 then
          f5[i][3*j-2] :=  s4_orient[j] * (z[j_v[1]]-z[j_v[4]])/(z[j_v[4]]-z[j_v[5]]);
          f5[i][3*j-1] :=  s4_orient[j] * (z[j_v[2]]-z[j_v[4]])/(z[j_v[4]]-z[j_v[5]]);
          f5[i][3*j  ] :=  s4_orient[j] * (z[j_v[3]]-z[j_v[4]])/(z[j_v[4]]-z[j_v[5]]);
        fi;
      fi;
    od;
  od;

  F := rec (
    f2 := f2*xxx/xxx,
    f3_full := f3_full*xxx/xxx,
    f4_full := f4_full*xxx/xxx,
    f5 := f5*xxx/xxx,
    last_inner_2face := ind_f2_li,
    last_inner_3face := ind_f3_li,
    z := z
  );

  return F;
end );

###############################################################################################################


InstallGlobalFunction( f3E4,
# считаем, что дан 4-симплекс four_simplex = [i1,i2,i3,i4,i5] и его двумерная грань two_face и ребро edge
# вычисляем производную квадрата синуса двугранного угла при two_face по квадрату длины edge'а
function( four_simplex, two_face, edge, coord )

local w, x, y, z # coordinate arrays
      , L12, L13, L14, L15, L23, L24, L25, L34, L35, L45 # squared lengths
      , w1, x1, y1, z1,  w2, x2, y2, z2,  w3, x3, y3, z3,  w4, x4, y4, z4,  w5, x5, y5, z5 # coordinates
      , n  # число вершин, общих для 2-грани и ребра
      , l  # список пяти вершин в нужном для данного случая порядке
      , CM123 # the Cayley-Menger determinant for triangle i1,i2,i3
      , W12345 # 24 * the (oriented) volume of 4-simplex i1,i2,i3,i4,i5
      , A, B # two auxiliary matrices
      , v1234, v1235 # two auxiliary vectors
      , r0, r1, r2  # three variants for the derivative of squared sine w.r.t. squared length
      , r # what we choose from r0, r1, r2
      , d  # result to return
;

# координаты
w := coord[1]; x := coord[2]; y := coord[3]; z := coord[4];

n := Length( Intersection( two_face, edge ) );

if n = 0 then l := Concatenation( two_face, edge );
  elif n = 1 then l := Concatenation( Intersection( two_face, edge ), Difference( two_face, edge ), 
                                Difference( four_simplex, UnionSet( two_face, edge ) ), Difference( edge, two_face ) );
  elif n = 2 then l := Concatenation( edge, Difference( two_face, edge ), Difference( four_simplex, two_face ) );
fi;

w1:=w[l[1]]; x1:=x[l[1]]; y1:=y[l[1]]; z1:=z[l[1]];
w2:=w[l[2]]; x2:=x[l[2]]; y2:=y[l[2]]; z2:=z[l[2]];
w3:=w[l[3]]; x3:=x[l[3]]; y3:=y[l[3]]; z3:=z[l[3]];
w4:=w[l[4]]; x4:=x[l[4]]; y4:=y[l[4]]; z4:=z[l[4]];
w5:=w[l[5]]; x5:=x[l[5]]; y5:=y[l[5]]; z5:=z[l[5]];

# так что теперь можно считать, что мы берём 4-симплекс 12345 
# и дифференцируем, в данный момент, квадрат синуса двугранного угла 123
# (1) по квадрату длины 45
# (2) по квадрату длины 15
# (3) по квадрату длины 12

# вначале пишем выражения квадратов длин через координаты

L12:= (w2-w1)^2+(x2-x1)^2+(y2-y1)^2+(z2-z1)^2; L13:= (w3-w1)^2+(x3-x1)^2+(y3-y1)^2+(z3-z1)^2; 
L14:= (w4-w1)^2+(x4-x1)^2+(y4-y1)^2+(z4-z1)^2; L15:= (w5-w1)^2+(x5-x1)^2+(y5-y1)^2+(z5-z1)^2; 
L23:= (w3-w2)^2+(x3-x2)^2+(y3-y2)^2+(z3-z2)^2; 
L24:= (w4-w2)^2+(x4-x2)^2+(y4-y2)^2+(z4-z2)^2; L25:= (w5-w2)^2+(x5-x2)^2+(y5-y2)^2+(z5-z2)^2; 
L34:= (w4-w3)^2+(x4-x3)^2+(y4-y3)^2+(z4-z3)^2; L35:= (w5-w3)^2+(x5-x3)^2+(y5-y3)^2+(z5-z3)^2; 
L45:= (w5-w4)^2+(x5-x4)^2+(y5-y4)^2+(z5-z4)^2;

# производную от квадрата синуса по квадрату длины, пока она ещё выражена через квадраты длин, факторизуем
# и в таком виде и используем

#####     случай (1) - по квадрату длины 45     #####

r0 := -(L23^2-2*L13*L23-2*L12*L23+L13^2-2*L12*L13+L12^2)
       *(L23^2*L45-2*L13*L23*L45-2*L12*L23*L45+L13^2*L45-2*L12*L13*L45
                  +L12^2*L45-2*L12*L34*L35-L23*L24*L35+L13*L24*L35+L12*L24*L35
                  +L14*L23*L35+L12*L23*L35-L13*L14*L35+L12*L14*L35+L12*L13*L35
                  -L12^2*L35-L23*L25*L34+L13*L25*L34+L12*L25*L34+L15*L23*L34
                  +L12*L23*L34-L13*L15*L34+L12*L15*L34+L12*L13*L34-L12^2*L34
                  -2*L13*L24*L25+L14*L23*L25+L13*L23*L25+L13*L14*L25
                  -L12*L14*L25-L13^2*L25+L12*L13*L25+L15*L23*L24+L13*L23*L24
                  +L13*L15*L24-L12*L15*L24-L13^2*L24+L12*L13*L24-L15*L23^2
                  -L14*L23^2-2*L14*L15*L23+L13*L15*L23+L12*L15*L23+L13*L14*L23
                  +L12*L14*L23-2*L12*L13*L23)
       /(2*(L12*L34^2+L23*L24*L34-L13*L24*L34-L12*L24*L34-L14*L23*L34
                     -L12*L23*L34+L13*L14*L34-L12*L14*L34-L12*L13*L34
                     +L12^2*L34+L13*L24^2-L14*L23*L24-L13*L23*L24-L13*L14*L24
                     +L12*L14*L24+L13^2*L24-L12*L13*L24+L14*L23^2+L14^2*L23
                     -L13*L14*L23-L12*L14*L23+L12*L13*L23)
          *(L12*L35^2+L23*L25*L35-L13*L25*L35-L12*L25*L35-L15*L23*L35
                     -L12*L23*L35+L13*L15*L35-L12*L15*L35-L12*L13*L35
                     +L12^2*L35+L13*L25^2-L15*L23*L25-L13*L23*L25-L13*L15*L25
                     +L12*L15*L25+L13^2*L25-L12*L13*L25+L15*L23^2+L15^2*L23
                     -L13*L15*L23-L12*L15*L23+L12*L13*L23))
  ;

#####     случай (2) - по квадрату длины 15     #####

r1 := -(L23^2-2*L13*L23-2*L12*L23+L13^2-2*L12*L13+L12^2)
       *(L23^2*L45-2*L13*L23*L45-2*L12*L23*L45+L13^2*L45-2*L12*L13*L45
                  +L12^2*L45-2*L12*L34*L35-L23*L24*L35+L13*L24*L35+L12*L24*L35
                  +L14*L23*L35+L12*L23*L35-L13*L14*L35+L12*L14*L35+L12*L13*L35
                  -L12^2*L35-L23*L25*L34+L13*L25*L34+L12*L25*L34+L15*L23*L34
                  +L12*L23*L34-L13*L15*L34+L12*L15*L34+L12*L13*L34-L12^2*L34
                  -2*L13*L24*L25+L14*L23*L25+L13*L23*L25+L13*L14*L25
                  -L12*L14*L25-L13^2*L25+L12*L13*L25+L15*L23*L24+L13*L23*L24
                  +L13*L15*L24-L12*L15*L24-L13^2*L24+L12*L13*L24-L15*L23^2
                  -L14*L23^2-2*L14*L15*L23+L13*L15*L23+L12*L15*L23+L13*L14*L23
                  +L12*L14*L23-2*L12*L13*L23)
       *(L23*L35*L45-L13*L35*L45+L12*L35*L45+L23*L25*L45+L13*L25*L45
                    -L12*L25*L45-L23^2*L45-2*L15*L23*L45+L13*L23*L45
                    +L12*L23*L45-L24*L35^2+L14*L35^2-L12*L35^2+L25*L34*L35
                    -L15*L34*L35+L12*L34*L35+L24*L25*L35-2*L23*L25*L35
                    -2*L14*L25*L35+L13*L25*L35+L12*L25*L35+L23*L24*L35
                    +L15*L24*L35+L13*L24*L35-2*L12*L24*L35+L15*L23*L35
                    -2*L14*L23*L35+L12*L23*L35-L25^2*L34+L23*L25*L34
                    +L15*L25*L34-2*L13*L25*L34+L12*L25*L34+L15*L23*L34
                    -L12*L23*L34+L14*L25^2-L13*L25^2-L15*L24*L25+L13*L24*L25
                    +L15*L23*L25-2*L14*L23*L25+L13*L23*L25+L15*L23*L24
                    -L13*L23*L24-L15*L23^2+L14*L23^2)
       /(4*(L12*L34^2+L23*L24*L34-L13*L24*L34-L12*L24*L34-L14*L23*L34
                     -L12*L23*L34+L13*L14*L34-L12*L14*L34-L12*L13*L34
                     +L12^2*L34+L13*L24^2-L14*L23*L24-L13*L23*L24-L13*L14*L24
                     +L12*L14*L24+L13^2*L24-L12*L13*L24+L14*L23^2+L14^2*L23
                     -L13*L14*L23-L12*L14*L23+L12*L13*L23)
          *(L12*L35^2+L23*L25*L35-L13*L25*L35-L12*L25*L35-L15*L23*L35
                     -L12*L23*L35+L13*L15*L35-L12*L15*L35-L12*L13*L35
                     +L12^2*L35+L13*L25^2-L15*L23*L25-L13*L23*L25-L13*L15*L25
                     +L12*L15*L25+L13^2*L25-L12*L13*L25+L15*L23^2+L15^2*L23
                     -L13*L15*L23-L12*L15*L23+L12*L13*L23)
           ^2)
  ;


#####     случай (3) - по квадрату длины 12     #####

r2 := (L23^2*L45-2*L13*L23*L45-2*L12*L23*L45+L13^2*L45-2*L12*L13*L45+L12^2*L45
                -2*L12*L34*L35-L23*L24*L35+L13*L24*L35+L12*L24*L35+L14*L23*L35
                +L12*L23*L35-L13*L14*L35+L12*L14*L35+L12*L13*L35-L12^2*L35
                -L23*L25*L34+L13*L25*L34+L12*L25*L34+L15*L23*L34+L12*L23*L34
                -L13*L15*L34+L12*L15*L34+L12*L13*L34-L12^2*L34-2*L13*L24*L25
                +L14*L23*L25+L13*L23*L25+L13*L14*L25-L12*L14*L25-L13^2*L25
                +L12*L13*L25+L15*L23*L24+L13*L23*L24+L13*L15*L24-L12*L15*L24
                -L13^2*L24+L12*L13*L24-L15*L23^2-L14*L23^2-2*L14*L15*L23
                +L13*L15*L23+L12*L15*L23+L13*L14*L23+L12*L14*L23
                -2*L12*L13*L23)
       *(2*L12*L23^2*L34^2*L35^2*L45-4*L12*L13*L23*L34^2*L35^2*L45
                                    +2*L12*L13^2*L34^2*L35^2*L45
                                    -2*L12^3*L34^2*L35^2*L45
                                    +L23^3*L24*L34*L35^2*L45
                                    -3*L13*L23^2*L24*L34*L35^2*L45
                                    +3*L13^2*L23*L24*L34*L35^2*L45
                                    +4*L12*L13*L23*L24*L34*L35^2*L45
                                    -3*L12^2*L23*L24*L34*L35^2*L45
                                    -L13^3*L24*L34*L35^2*L45
                                    -4*L12*L13^2*L24*L34*L35^2*L45
                                    +3*L12^2*L13*L24*L34*L35^2*L45
                                    +2*L12^3*L24*L34*L35^2*L45
                                    -L14*L23^3*L34*L35^2*L45
                                    -2*L12*L23^3*L34*L35^2*L45
                                    +3*L13*L14*L23^2*L34*L35^2*L45
                                    -4*L12*L14*L23^2*L34*L35^2*L45
                                    +2*L12*L13*L23^2*L34*L35^2*L45
                                    +3*L12^2*L23^2*L34*L35^2*L45
                                    -3*L13^2*L14*L23*L34*L35^2*L45
                                    +4*L12*L13*L14*L23*L34*L35^2*L45
                                    +3*L12^2*L14*L23*L34*L35^2*L45
                                    +2*L12*L13^2*L23*L34*L35^2*L45
                                    -6*L12^2*L13*L23*L34*L35^2*L45
                                    +L13^3*L14*L34*L35^2*L45
                                    -3*L12^2*L13*L14*L34*L35^2*L45
                                    +2*L12^3*L14*L34*L35^2*L45
                                    -2*L12*L13^3*L34*L35^2*L45
                                    +3*L12^2*L13^2*L34*L35^2*L45
                                    -L12^4*L34*L35^2*L45
                                    +L13*L23^2*L24^2*L35^2*L45
                                    -2*L13^2*L23*L24^2*L35^2*L45
                                    +2*L12*L13*L23*L24^2*L35^2*L45
                                    +L13^3*L24^2*L35^2*L45
                                    +2*L12*L13^2*L24^2*L35^2*L45
                                    -3*L12^2*L13*L24^2*L35^2*L45
                                    -L14*L23^3*L24*L35^2*L45
                                    -L13*L23^3*L24*L35^2*L45
                                    +L13*L14*L23^2*L24*L35^2*L45
                                    +3*L13^2*L23^2*L24*L35^2*L45
                                    -4*L12*L13*L23^2*L24*L35^2*L45
                                    +L13^2*L14*L23*L24*L35^2*L45
                                    -8*L12*L13*L14*L23*L24*L35^2*L45
                                    +3*L12^2*L14*L23*L24*L35^2*L45
                                    -3*L13^3*L23*L24*L35^2*L45
                                    +4*L12*L13^2*L23*L24*L35^2*L45
                                    +3*L12^2*L13*L23*L24*L35^2*L45
                                    -L13^3*L14*L24*L35^2*L45
                                    +3*L12^2*L13*L14*L24*L35^2*L45
                                    -2*L12^3*L14*L24*L35^2*L45
                                    +L13^4*L24*L35^2*L45
                                    -3*L12^2*L13^2*L24*L35^2*L45
                                    +2*L12^3*L13*L24*L35^2*L45
                                    +L14*L23^4*L35^2*L45+L14^2*L23^3*L35^2*L45
                                    -3*L13*L14*L23^3*L35^2*L45
                                    +2*L12*L13*L23^3*L35^2*L45
                                    -2*L13*L14^2*L23^2*L35^2*L45
                                    +2*L12*L14^2*L23^2*L35^2*L45
                                    +3*L13^2*L14*L23^2*L35^2*L45
                                    +4*L12*L13*L14*L23^2*L35^2*L45
                                    -3*L12^2*L14*L23^2*L35^2*L45
                                    -4*L12*L13^2*L23^2*L35^2*L45
                                    +L13^2*L14^2*L23*L35^2*L45
                                    +2*L12*L13*L14^2*L23*L35^2*L45
                                    -3*L12^2*L14^2*L23*L35^2*L45
                                    -L13^3*L14*L23*L35^2*L45
                                    -4*L12*L13^2*L14*L23*L35^2*L45
                                    +3*L12^2*L13*L14*L23*L35^2*L45
                                    +2*L12^3*L14*L23*L35^2*L45
                                    +2*L12*L13^3*L23*L35^2*L45
                                    -2*L12^3*L13*L23*L35^2*L45
                                    +L23^3*L25*L34^2*L35*L45
                                    -3*L13*L23^2*L25*L34^2*L35*L45
                                    +3*L13^2*L23*L25*L34^2*L35*L45
                                    +4*L12*L13*L23*L25*L34^2*L35*L45
                                    -3*L12^2*L23*L25*L34^2*L35*L45
                                    -L13^3*L25*L34^2*L35*L45
                                    -4*L12*L13^2*L25*L34^2*L35*L45
                                    +3*L12^2*L13*L25*L34^2*L35*L45
                                    +2*L12^3*L25*L34^2*L35*L45
                                    -L15*L23^3*L34^2*L35*L45
                                    -2*L12*L23^3*L34^2*L35*L45
                                    +3*L13*L15*L23^2*L34^2*L35*L45
                                    -4*L12*L15*L23^2*L34^2*L35*L45
                                    +2*L12*L13*L23^2*L34^2*L35*L45
                                    +3*L12^2*L23^2*L34^2*L35*L45
                                    -3*L13^2*L15*L23*L34^2*L35*L45
                                    +4*L12*L13*L15*L23*L34^2*L35*L45
                                    +3*L12^2*L15*L23*L34^2*L35*L45
                                    +2*L12*L13^2*L23*L34^2*L35*L45
                                    -6*L12^2*L13*L23*L34^2*L35*L45
                                    +L13^3*L15*L34^2*L35*L45
                                    -3*L12^2*L13*L15*L34^2*L35*L45
                                    +2*L12^3*L15*L34^2*L35*L45
                                    -2*L12*L13^3*L34^2*L35*L45
                                    +3*L12^2*L13^2*L34^2*L35*L45
                                    -L12^4*L34^2*L35*L45
                                    +2*L23^3*L24*L25*L34*L35*L45
                                    +2*L13*L23^2*L24*L25*L34*L35*L45
                                    -6*L12*L23^2*L24*L25*L34*L35*L45
                                    -10*L13^2*L23*L24*L25*L34*L35*L45
                                    +4*L12*L13*L23*L24*L25*L34*L35*L45
                                    +6*L12^2*L23*L24*L25*L34*L35*L45
                                    +6*L13^3*L24*L25*L34*L35*L45
                                    +2*L12*L13^2*L24*L25*L34*L35*L45
                                    -6*L12^2*L13*L24*L25*L34*L35*L45
                                    -2*L12^3*L24*L25*L34*L35*L45
                                    -L23^4*L25*L34*L35*L45
                                    -4*L14*L23^3*L25*L34*L35*L45
                                    +2*L13*L23^3*L25*L34*L35*L45
                                    +2*L12*L23^3*L25*L34*L35*L45
                                    +4*L13*L14*L23^2*L25*L34*L35*L45
                                    +6*L12*L14*L23^2*L25*L34*L35*L45
                                    -10*L12*L13*L23^2*L25*L34*L35*L45
                                    +4*L13^2*L14*L23*L25*L34*L35*L45
                                    -12*L12*L13*L14*L23*L25*L34*L35*L45
                                    -2*L13^3*L23*L25*L34*L35*L45
                                    +6*L12*L13^2*L23*L25*L34*L35*L45
                                    +6*L12^2*L13*L23*L25*L34*L35*L45
                                    -2*L12^3*L23*L25*L34*L35*L45
                                    -4*L13^3*L14*L25*L34*L35*L45
                                    +6*L12*L13^2*L14*L25*L34*L35*L45
                                    -2*L12^3*L14*L25*L34*L35*L45
                                    +L13^4*L25*L34*L35*L45
                                    +2*L12*L13^3*L25*L34*L35*L45
                                    -6*L12^2*L13^2*L25*L34*L35*L45
                                    +2*L12^3*L13*L25*L34*L35*L45
                                    +L12^4*L25*L34*L35*L45
                                    -L23^4*L24*L34*L35*L45
                                    -4*L15*L23^3*L24*L34*L35*L45
                                    +2*L13*L23^3*L24*L34*L35*L45
                                    +2*L12*L23^3*L24*L34*L35*L45
                                    +4*L13*L15*L23^2*L24*L34*L35*L45
                                    +6*L12*L15*L23^2*L24*L34*L35*L45
                                    -10*L12*L13*L23^2*L24*L34*L35*L45
                                    +4*L13^2*L15*L23*L24*L34*L35*L45
                                    -12*L12*L13*L15*L23*L24*L34*L35*L45
                                    -2*L13^3*L23*L24*L34*L35*L45
                                    +6*L12*L13^2*L23*L24*L34*L35*L45
                                    +6*L12^2*L13*L23*L24*L34*L35*L45
                                    -2*L12^3*L23*L24*L34*L35*L45
                                    -4*L13^3*L15*L24*L34*L35*L45
                                    +6*L12*L13^2*L15*L24*L34*L35*L45
                                    -2*L12^3*L15*L24*L34*L35*L45
                                    +L13^4*L24*L34*L35*L45
                                    +2*L12*L13^3*L24*L34*L35*L45
                                    -6*L12^2*L13^2*L24*L34*L35*L45
                                    +2*L12^3*L13*L24*L34*L35*L45
                                    +L12^4*L24*L34*L35*L45
                                    +L15*L23^4*L34*L35*L45
                                    +L14*L23^4*L34*L35*L45
                                    +2*L12*L23^4*L34*L35*L45
                                    +6*L14*L15*L23^3*L34*L35*L45
                                    -2*L13*L15*L23^3*L34*L35*L45
                                    +2*L12*L15*L23^3*L34*L35*L45
                                    -2*L13*L14*L23^3*L34*L35*L45
                                    +2*L12*L14*L23^3*L34*L35*L45
                                    -6*L12^2*L23^3*L34*L35*L45
                                    -10*L13*L14*L15*L23^2*L34*L35*L45
                                    +2*L12*L14*L15*L23^2*L34*L35*L45
                                    +6*L12*L13*L15*L23^2*L34*L35*L45
                                    -6*L12^2*L15*L23^2*L34*L35*L45
                                    +6*L12*L13*L14*L23^2*L34*L35*L45
                                    -6*L12^2*L14*L23^2*L34*L35*L45
                                    -4*L12*L13^2*L23^2*L34*L35*L45
                                    +6*L12^2*L13*L23^2*L34*L35*L45
                                    +6*L12^3*L23^2*L34*L35*L45
                                    +2*L13^2*L14*L15*L23*L34*L35*L45
                                    +4*L12*L13*L14*L15*L23*L34*L35*L45
                                    -6*L12^2*L14*L15*L23*L34*L35*L45
                                    +2*L13^3*L15*L23*L34*L35*L45
                                    -10*L12*L13^2*L15*L23*L34*L35*L45
                                    +6*L12^2*L13*L15*L23*L34*L35*L45
                                    +2*L12^3*L15*L23*L34*L35*L45
                                    +2*L13^3*L14*L23*L34*L35*L45
                                    -10*L12*L13^2*L14*L23*L34*L35*L45
                                    +6*L12^2*L13*L14*L23*L34*L35*L45
                                    +2*L12^3*L14*L23*L34*L35*L45
                                    +6*L12^2*L13^2*L23*L34*L35*L45
                                    -4*L12^3*L13*L23*L34*L35*L45
                                    -2*L12^4*L23*L34*L35*L45
                                    +2*L13^3*L14*L15*L34*L35*L45
                                    -6*L12*L13^2*L14*L15*L34*L35*L45
                                    +6*L12^2*L13*L14*L15*L34*L35*L45
                                    -2*L12^3*L14*L15*L34*L35*L45
                                    -L13^4*L15*L34*L35*L45
                                    +2*L12*L13^3*L15*L34*L35*L45
                                    -2*L12^3*L13*L15*L34*L35*L45
                                    +L12^4*L15*L34*L35*L45
                                    -L13^4*L14*L34*L35*L45
                                    +2*L12*L13^3*L14*L34*L35*L45
                                    -2*L12^3*L13*L14*L34*L35*L45
                                    +L12^4*L14*L34*L35*L45
                                    +2*L12*L13^4*L34*L35*L45
                                    -6*L12^2*L13^3*L34*L35*L45
                                    +6*L12^3*L13^2*L34*L35*L45
                                    -2*L12^4*L13*L34*L35*L45
                                    +3*L13*L23^2*L24^2*L25*L35*L45
                                    +2*L13^2*L23*L24^2*L25*L35*L45
                                    -6*L12*L13*L23*L24^2*L25*L35*L45
                                    -5*L13^3*L24^2*L25*L35*L45
                                    +2*L12*L13^2*L24^2*L25*L35*L45
                                    +3*L12^2*L13*L24^2*L25*L35*L45
                                    -2*L14*L23^3*L24*L25*L35*L45
                                    -4*L13*L23^3*L24*L25*L35*L45
                                    -8*L13*L14*L23^2*L24*L25*L35*L45
                                    +6*L12*L14*L23^2*L24*L25*L35*L45
                                    +4*L13^2*L23^2*L24*L25*L35*L45
                                    +6*L12*L13*L23^2*L24*L25*L35*L45
                                    +6*L13^2*L14*L23*L24*L25*L35*L45
                                    +8*L12*L13*L14*L23*L24*L25*L35*L45
                                    -6*L12^2*L14*L23*L24*L25*L35*L45
                                    +4*L13^3*L23*L24*L25*L35*L45
                                    -12*L12*L13^2*L23*L24*L25*L35*L45
                                    +4*L13^3*L14*L24*L25*L35*L45
                                    -6*L12*L13^2*L14*L24*L25*L35*L45
                                    +2*L12^3*L14*L24*L25*L35*L45
                                    -4*L13^4*L24*L25*L35*L45
                                    +6*L12*L13^3*L24*L25*L35*L45
                                    -2*L12^3*L13*L24*L25*L35*L45
                                    +2*L14*L23^4*L25*L35*L45
                                    +L13*L23^4*L25*L35*L45
                                    +3*L14^2*L23^3*L25*L35*L45
                                    +2*L13*L14*L23^3*L25*L35*L45
                                    -6*L12*L14*L23^3*L25*L35*L45
                                    -3*L13^2*L23^3*L25*L35*L45
                                    +2*L13*L14^2*L23^2*L25*L35*L45
                                    -6*L12*L14^2*L23^2*L25*L35*L45
                                    -10*L13^2*L14*L23^2*L25*L35*L45
                                    +4*L12*L13*L14*L23^2*L25*L35*L45
                                    +6*L12^2*L14*L23^2*L25*L35*L45
                                    +3*L13^3*L23^2*L25*L35*L45
                                    +4*L12*L13^2*L23^2*L25*L35*L45
                                    -3*L12^2*L13*L23^2*L25*L35*L45
                                    -5*L13^2*L14^2*L23*L25*L35*L45
                                    +2*L12*L13*L14^2*L23*L25*L35*L45
                                    +3*L12^2*L14^2*L23*L25*L35*L45
                                    +6*L13^3*L14*L23*L25*L35*L45
                                    +2*L12*L13^2*L14*L23*L25*L35*L45
                                    -6*L12^2*L13*L14*L23*L25*L35*L45
                                    -2*L12^3*L14*L23*L25*L35*L45
                                    -L13^4*L23*L25*L35*L45
                                    -4*L12*L13^3*L23*L25*L35*L45
                                    +3*L12^2*L13^2*L23*L25*L35*L45
                                    +2*L12^3*L13*L23*L25*L35*L45
                                    -L13*L23^3*L24^2*L35*L45
                                    -5*L13*L15*L23^2*L24^2*L35*L45
                                    +L13^2*L23^2*L24^2*L35*L45
                                    +2*L13^2*L15*L23*L24^2*L35*L45
                                    +2*L12*L13*L15*L23*L24^2*L35*L45
                                    +L13^3*L23*L24^2*L35*L45
                                    -8*L12*L13^2*L23*L24^2*L35*L45
                                    +3*L12^2*L13*L23*L24^2*L35*L45
                                    +3*L13^3*L15*L24^2*L35*L45
                                    -6*L12*L13^2*L15*L24^2*L35*L45
                                    +3*L12^2*L13*L15*L24^2*L35*L45
                                    -L13^4*L24^2*L35*L45
                                    +3*L12^2*L13^2*L24^2*L35*L45
                                    -2*L12^3*L13*L24^2*L35*L45
                                    +L14*L23^4*L24*L35*L45
                                    +L13*L23^4*L24*L35*L45
                                    +4*L14*L15*L23^3*L24*L35*L45
                                    +6*L13*L15*L23^3*L24*L35*L45
                                    -2*L12*L14*L23^3*L24*L35*L45
                                    -2*L13^2*L23^3*L24*L35*L45
                                    +2*L12*L13*L23^3*L24*L35*L45
                                    +6*L13*L14*L15*L23^2*L24*L35*L45
                                    -6*L12*L14*L15*L23^2*L24*L35*L45
                                    -10*L13^2*L15*L23^2*L24*L35*L45
                                    +2*L12*L13*L15*L23^2*L24*L35*L45
                                    -2*L13^2*L14*L23^2*L24*L35*L45
                                    +10*L12*L13*L14*L23^2*L24*L35*L45
                                    +6*L12*L13^2*L23^2*L24*L35*L45
                                    -6*L12^2*L13*L23^2*L24*L35*L45
                                    -8*L13^2*L14*L15*L23*L24*L35*L45
                                    +8*L12*L13*L14*L15*L23*L24*L35*L45
                                    +2*L13^3*L15*L23*L24*L35*L45
                                    +4*L12*L13^2*L15*L23*L24*L35*L45
                                    -6*L12^2*L13*L15*L23*L24*L35*L45
                                    +10*L12*L13^2*L14*L23*L24*L35*L45
                                    -12*L12^2*L13*L14*L23*L24*L35*L45
                                    +2*L12^3*L14*L23*L24*L35*L45
                                    +2*L13^4*L23*L24*L35*L45
                                    -10*L12*L13^3*L23*L24*L35*L45
                                    +6*L12^2*L13^2*L23*L24*L35*L45
                                    +2*L12^3*L13*L23*L24*L35*L45
                                    -2*L13^3*L14*L15*L24*L35*L45
                                    +6*L12*L13^2*L14*L15*L24*L35*L45
                                    -6*L12^2*L13*L14*L15*L24*L35*L45
                                    +2*L12^3*L14*L15*L24*L35*L45
                                    +2*L13^4*L15*L24*L35*L45
                                    -6*L12*L13^3*L15*L24*L35*L45
                                    +6*L12^2*L13^2*L15*L24*L35*L45
                                    -2*L12^3*L13*L15*L24*L35*L45
                                    +L13^4*L14*L24*L35*L45
                                    -2*L12*L13^3*L14*L24*L35*L45
                                    +2*L12^3*L13*L14*L24*L35*L45
                                    -L12^4*L14*L24*L35*L45-L13^5*L24*L35*L45
                                    +2*L12*L13^4*L24*L35*L45
                                    -2*L12^3*L13^2*L24*L35*L45
                                    +L12^4*L13*L24*L35*L45-L14*L23^5*L35*L45
                                    -4*L14*L15*L23^4*L35*L45
                                    -L13*L15*L23^4*L35*L45-L14^2*L23^4*L35*L45
                                    +2*L13*L14*L23^4*L35*L45
                                    +2*L12*L14*L23^4*L35*L45
                                    -2*L12*L13*L23^4*L35*L45
                                    -5*L14^2*L15*L23^3*L35*L45
                                    +4*L13*L14*L15*L23^3*L35*L45
                                    +6*L12*L14*L15*L23^3*L35*L45
                                    +3*L13^2*L15*L23^3*L35*L45
                                    -4*L12*L13*L15*L23^3*L35*L45
                                    +L13*L14^2*L23^3*L35*L45
                                    -10*L12*L13*L14*L23^3*L35*L45
                                    +2*L12*L13^2*L23^3*L35*L45
                                    +3*L12^2*L13*L23^3*L35*L45
                                    +2*L13*L14^2*L15*L23^2*L35*L45
                                    +2*L12*L14^2*L15*L23^2*L35*L45
                                    +4*L13^2*L14*L15*L23^2*L35*L45
                                    -12*L12*L13*L14*L15*L23^2*L35*L45
                                    -3*L13^3*L15*L23^2*L35*L45
                                    +4*L12*L13^2*L15*L23^2*L35*L45
                                    +3*L12^2*L13*L15*L23^2*L35*L45
                                    +L13^2*L14^2*L23^2*L35*L45
                                    -8*L12*L13*L14^2*L23^2*L35*L45
                                    +3*L12^2*L14^2*L23^2*L35*L45
                                    -2*L13^3*L14*L23^2*L35*L45
                                    +6*L12*L13^2*L14*L23^2*L35*L45
                                    +6*L12^2*L13*L14*L23^2*L35*L45
                                    -2*L12^3*L14*L23^2*L35*L45
                                    +2*L12*L13^3*L23^2*L35*L45
                                    -6*L12^2*L13^2*L23^2*L35*L45
                                    +3*L13^2*L14^2*L15*L23*L35*L45
                                    -6*L12*L13*L14^2*L15*L23*L35*L45
                                    +3*L12^2*L14^2*L15*L23*L35*L45
                                    -4*L13^3*L14*L15*L23*L35*L45
                                    +6*L12*L13^2*L14*L15*L23*L35*L45
                                    -2*L12^3*L14*L15*L23*L35*L45
                                    +L13^4*L15*L23*L35*L45
                                    -3*L12^2*L13^2*L15*L23*L35*L45
                                    +2*L12^3*L13*L15*L23*L35*L45
                                    -L13^3*L14^2*L23*L35*L45
                                    +3*L12^2*L13*L14^2*L23*L35*L45
                                    -2*L12^3*L14^2*L23*L35*L45
                                    +L13^4*L14*L23*L35*L45
                                    +2*L12*L13^3*L14*L23*L35*L45
                                    -6*L12^2*L13^2*L14*L23*L35*L45
                                    +2*L12^3*L13*L14*L23*L35*L45
                                    +L12^4*L14*L23*L35*L45
                                    -2*L12*L13^4*L23*L35*L45
                                    +3*L12^2*L13^3*L23*L35*L45
                                    -L12^4*L13*L23*L35*L45
                                    +L13*L23^2*L25^2*L34^2*L45
                                    -2*L13^2*L23*L25^2*L34^2*L45
                                    +2*L12*L13*L23*L25^2*L34^2*L45
                                    +L13^3*L25^2*L34^2*L45
                                    +2*L12*L13^2*L25^2*L34^2*L45
                                    -3*L12^2*L13*L25^2*L34^2*L45
                                    -L15*L23^3*L25*L34^2*L45
                                    -L13*L23^3*L25*L34^2*L45
                                    +L13*L15*L23^2*L25*L34^2*L45
                                    +3*L13^2*L23^2*L25*L34^2*L45
                                    -4*L12*L13*L23^2*L25*L34^2*L45
                                    +L13^2*L15*L23*L25*L34^2*L45
                                    -8*L12*L13*L15*L23*L25*L34^2*L45
                                    +3*L12^2*L15*L23*L25*L34^2*L45
                                    -3*L13^3*L23*L25*L34^2*L45
                                    +4*L12*L13^2*L23*L25*L34^2*L45
                                    +3*L12^2*L13*L23*L25*L34^2*L45
                                    -L13^3*L15*L25*L34^2*L45
                                    +3*L12^2*L13*L15*L25*L34^2*L45
                                    -2*L12^3*L15*L25*L34^2*L45
                                    +L13^4*L25*L34^2*L45
                                    -3*L12^2*L13^2*L25*L34^2*L45
                                    +2*L12^3*L13*L25*L34^2*L45
                                    +L15*L23^4*L34^2*L45+L15^2*L23^3*L34^2*L45
                                    -3*L13*L15*L23^3*L34^2*L45
                                    +2*L12*L13*L23^3*L34^2*L45
                                    -2*L13*L15^2*L23^2*L34^2*L45
                                    +2*L12*L15^2*L23^2*L34^2*L45
                                    +3*L13^2*L15*L23^2*L34^2*L45
                                    +4*L12*L13*L15*L23^2*L34^2*L45
                                    -3*L12^2*L15*L23^2*L34^2*L45
                                    -4*L12*L13^2*L23^2*L34^2*L45
                                    +L13^2*L15^2*L23*L34^2*L45
                                    +2*L12*L13*L15^2*L23*L34^2*L45
                                    -3*L12^2*L15^2*L23*L34^2*L45
                                    -L13^3*L15*L23*L34^2*L45
                                    -4*L12*L13^2*L15*L23*L34^2*L45
                                    +3*L12^2*L13*L15*L23*L34^2*L45
                                    +2*L12^3*L15*L23*L34^2*L45
                                    +2*L12*L13^3*L23*L34^2*L45
                                    -2*L12^3*L13*L23*L34^2*L45
                                    +3*L13*L23^2*L24*L25^2*L34*L45
                                    +2*L13^2*L23*L24*L25^2*L34*L45
                                    -6*L12*L13*L23*L24*L25^2*L34*L45
                                    -5*L13^3*L24*L25^2*L34*L45
                                    +2*L12*L13^2*L24*L25^2*L34*L45
                                    +3*L12^2*L13*L24*L25^2*L34*L45
                                    -L13*L23^3*L25^2*L34*L45
                                    -5*L13*L14*L23^2*L25^2*L34*L45
                                    +L13^2*L23^2*L25^2*L34*L45
                                    +2*L13^2*L14*L23*L25^2*L34*L45
                                    +2*L12*L13*L14*L23*L25^2*L34*L45
                                    +L13^3*L23*L25^2*L34*L45
                                    -8*L12*L13^2*L23*L25^2*L34*L45
                                    +3*L12^2*L13*L23*L25^2*L34*L45
                                    +3*L13^3*L14*L25^2*L34*L45
                                    -6*L12*L13^2*L14*L25^2*L34*L45
                                    +3*L12^2*L13*L14*L25^2*L34*L45
                                    -L13^4*L25^2*L34*L45
                                    +3*L12^2*L13^2*L25^2*L34*L45
                                    -2*L12^3*L13*L25^2*L34*L45
                                    -2*L15*L23^3*L24*L25*L34*L45
                                    -4*L13*L23^3*L24*L25*L34*L45
                                    -8*L13*L15*L23^2*L24*L25*L34*L45
                                    +6*L12*L15*L23^2*L24*L25*L34*L45
                                    +4*L13^2*L23^2*L24*L25*L34*L45
                                    +6*L12*L13*L23^2*L24*L25*L34*L45
                                    +6*L13^2*L15*L23*L24*L25*L34*L45
                                    +8*L12*L13*L15*L23*L24*L25*L34*L45
                                    -6*L12^2*L15*L23*L24*L25*L34*L45
                                    +4*L13^3*L23*L24*L25*L34*L45
                                    -12*L12*L13^2*L23*L24*L25*L34*L45
                                    +4*L13^3*L15*L24*L25*L34*L45
                                    -6*L12*L13^2*L15*L24*L25*L34*L45
                                    +2*L12^3*L15*L24*L25*L34*L45
                                    -4*L13^4*L24*L25*L34*L45
                                    +6*L12*L13^3*L24*L25*L34*L45
                                    -2*L12^3*L13*L24*L25*L34*L45
                                    +L15*L23^4*L25*L34*L45
                                    +L13*L23^4*L25*L34*L45
                                    +4*L14*L15*L23^3*L25*L34*L45
                                    -2*L12*L15*L23^3*L25*L34*L45
                                    +6*L13*L14*L23^3*L25*L34*L45
                                    -2*L13^2*L23^3*L25*L34*L45
                                    +2*L12*L13*L23^3*L25*L34*L45
                                    +6*L13*L14*L15*L23^2*L25*L34*L45
                                    -6*L12*L14*L15*L23^2*L25*L34*L45
                                    -2*L13^2*L15*L23^2*L25*L34*L45
                                    +10*L12*L13*L15*L23^2*L25*L34*L45
                                    -10*L13^2*L14*L23^2*L25*L34*L45
                                    +2*L12*L13*L14*L23^2*L25*L34*L45
                                    +6*L12*L13^2*L23^2*L25*L34*L45
                                    -6*L12^2*L13*L23^2*L25*L34*L45
                                    -8*L13^2*L14*L15*L23*L25*L34*L45
                                    +8*L12*L13*L14*L15*L23*L25*L34*L45
                                    +10*L12*L13^2*L15*L23*L25*L34*L45
                                    -12*L12^2*L13*L15*L23*L25*L34*L45
                                    +2*L12^3*L15*L23*L25*L34*L45
                                    +2*L13^3*L14*L23*L25*L34*L45
                                    +4*L12*L13^2*L14*L23*L25*L34*L45
                                    -6*L12^2*L13*L14*L23*L25*L34*L45
                                    +2*L13^4*L23*L25*L34*L45
                                    -10*L12*L13^3*L23*L25*L34*L45
                                    +6*L12^2*L13^2*L23*L25*L34*L45
                                    +2*L12^3*L13*L23*L25*L34*L45
                                    -2*L13^3*L14*L15*L25*L34*L45
                                    +6*L12*L13^2*L14*L15*L25*L34*L45
                                    -6*L12^2*L13*L14*L15*L25*L34*L45
                                    +2*L12^3*L14*L15*L25*L34*L45
                                    +L13^4*L15*L25*L34*L45
                                    -2*L12*L13^3*L15*L25*L34*L45
                                    +2*L12^3*L13*L15*L25*L34*L45
                                    -L12^4*L15*L25*L34*L45
                                    +2*L13^4*L14*L25*L34*L45
                                    -6*L12*L13^3*L14*L25*L34*L45
                                    +6*L12^2*L13^2*L14*L25*L34*L45
                                    -2*L12^3*L13*L14*L25*L34*L45
                                    -L13^5*L25*L34*L45+2*L12*L13^4*L25*L34*L45
                                    -2*L12^3*L13^2*L25*L34*L45
                                    +L12^4*L13*L25*L34*L45
                                    +2*L15*L23^4*L24*L34*L45
                                    +L13*L23^4*L24*L34*L45
                                    +3*L15^2*L23^3*L24*L34*L45
                                    +2*L13*L15*L23^3*L24*L34*L45
                                    -6*L12*L15*L23^3*L24*L34*L45
                                    -3*L13^2*L23^3*L24*L34*L45
                                    +2*L13*L15^2*L23^2*L24*L34*L45
                                    -6*L12*L15^2*L23^2*L24*L34*L45
                                    -10*L13^2*L15*L23^2*L24*L34*L45
                                    +4*L12*L13*L15*L23^2*L24*L34*L45
                                    +6*L12^2*L15*L23^2*L24*L34*L45
                                    +3*L13^3*L23^2*L24*L34*L45
                                    +4*L12*L13^2*L23^2*L24*L34*L45
                                    -3*L12^2*L13*L23^2*L24*L34*L45
                                    -5*L13^2*L15^2*L23*L24*L34*L45
                                    +2*L12*L13*L15^2*L23*L24*L34*L45
                                    +3*L12^2*L15^2*L23*L24*L34*L45
                                    +6*L13^3*L15*L23*L24*L34*L45
                                    +2*L12*L13^2*L15*L23*L24*L34*L45
                                    -6*L12^2*L13*L15*L23*L24*L34*L45
                                    -2*L12^3*L15*L23*L24*L34*L45
                                    -L13^4*L23*L24*L34*L45
                                    -4*L12*L13^3*L23*L24*L34*L45
                                    +3*L12^2*L13^2*L23*L24*L34*L45
                                    +2*L12^3*L13*L23*L24*L34*L45
                                    -L15*L23^5*L34*L45-L15^2*L23^4*L34*L45
                                    -4*L14*L15*L23^4*L34*L45
                                    +2*L13*L15*L23^4*L34*L45
                                    +2*L12*L15*L23^4*L34*L45
                                    -L13*L14*L23^4*L34*L45
                                    -2*L12*L13*L23^4*L34*L45
                                    -5*L14*L15^2*L23^3*L34*L45
                                    +L13*L15^2*L23^3*L34*L45
                                    +4*L13*L14*L15*L23^3*L34*L45
                                    +6*L12*L14*L15*L23^3*L34*L45
                                    -10*L12*L13*L15*L23^3*L34*L45
                                    +3*L13^2*L14*L23^3*L34*L45
                                    -4*L12*L13*L14*L23^3*L34*L45
                                    +2*L12*L13^2*L23^3*L34*L45
                                    +3*L12^2*L13*L23^3*L34*L45
                                    +2*L13*L14*L15^2*L23^2*L34*L45
                                    +2*L12*L14*L15^2*L23^2*L34*L45
                                    +L13^2*L15^2*L23^2*L34*L45
                                    -8*L12*L13*L15^2*L23^2*L34*L45
                                    +3*L12^2*L15^2*L23^2*L34*L45
                                    +4*L13^2*L14*L15*L23^2*L34*L45
                                    -12*L12*L13*L14*L15*L23^2*L34*L45
                                    -2*L13^3*L15*L23^2*L34*L45
                                    +6*L12*L13^2*L15*L23^2*L34*L45
                                    +6*L12^2*L13*L15*L23^2*L34*L45
                                    -2*L12^3*L15*L23^2*L34*L45
                                    -3*L13^3*L14*L23^2*L34*L45
                                    +4*L12*L13^2*L14*L23^2*L34*L45
                                    +3*L12^2*L13*L14*L23^2*L34*L45
                                    +2*L12*L13^3*L23^2*L34*L45
                                    -6*L12^2*L13^2*L23^2*L34*L45
                                    +3*L13^2*L14*L15^2*L23*L34*L45
                                    -6*L12*L13*L14*L15^2*L23*L34*L45
                                    +3*L12^2*L14*L15^2*L23*L34*L45
                                    -L13^3*L15^2*L23*L34*L45
                                    +3*L12^2*L13*L15^2*L23*L34*L45
                                    -2*L12^3*L15^2*L23*L34*L45
                                    -4*L13^3*L14*L15*L23*L34*L45
                                    +6*L12*L13^2*L14*L15*L23*L34*L45
                                    -2*L12^3*L14*L15*L23*L34*L45
                                    +L13^4*L15*L23*L34*L45
                                    +2*L12*L13^3*L15*L23*L34*L45
                                    -6*L12^2*L13^2*L15*L23*L34*L45
                                    +2*L12^3*L13*L15*L23*L34*L45
                                    +L12^4*L15*L23*L34*L45
                                    +L13^4*L14*L23*L34*L45
                                    -3*L12^2*L13^2*L14*L23*L34*L45
                                    +2*L12^3*L13*L14*L23*L34*L45
                                    -2*L12*L13^4*L23*L34*L45
                                    +3*L12^2*L13^3*L23*L34*L45
                                    -L12^4*L13*L23*L34*L45
                                    +4*L13^2*L23*L24^2*L25^2*L45
                                    +4*L13^3*L24^2*L25^2*L45
                                    -4*L12*L13^2*L24^2*L25^2*L45
                                    -3*L13*L14*L23^2*L24*L25^2*L45
                                    -5*L13^2*L23^2*L24*L25^2*L45
                                    -10*L13^2*L14*L23*L24*L25^2*L45
                                    +6*L12*L13*L14*L23*L24*L25^2*L45
                                    +2*L13^3*L23*L24*L25^2*L45
                                    +2*L12*L13^2*L23*L24*L25^2*L45
                                    -3*L13^3*L14*L24*L25^2*L45
                                    +6*L12*L13^2*L14*L24*L25^2*L45
                                    -3*L12^2*L13*L14*L24*L25^2*L45
                                    +3*L13^4*L24*L25^2*L45
                                    -6*L12*L13^3*L24*L25^2*L45
                                    +3*L12^2*L13^2*L24*L25^2*L45
                                    +3*L13*L14*L23^3*L25^2*L45
                                    +L13^2*L23^3*L25^2*L45
                                    +4*L13*L14^2*L23^2*L25^2*L45
                                    +2*L13^2*L14*L23^2*L25^2*L45
                                    -6*L12*L13*L14*L23^2*L25^2*L45
                                    -2*L13^3*L23^2*L25^2*L45
                                    +2*L12*L13^2*L23^2*L25^2*L45
                                    +4*L13^2*L14^2*L23*L25^2*L45
                                    -4*L12*L13*L14^2*L23*L25^2*L45
                                    -5*L13^3*L14*L23*L25^2*L45
                                    +2*L12*L13^2*L14*L23*L25^2*L45
                                    +3*L12^2*L13*L14*L23*L25^2*L45
                                    +L13^4*L23*L25^2*L45
                                    +2*L12*L13^3*L23*L25^2*L45
                                    -3*L12^2*L13^2*L23*L25^2*L45
                                    -3*L13*L15*L23^2*L24^2*L25*L45
                                    -5*L13^2*L23^2*L24^2*L25*L45
                                    -10*L13^2*L15*L23*L24^2*L25*L45
                                    +6*L12*L13*L15*L23*L24^2*L25*L45
                                    +2*L13^3*L23*L24^2*L25*L45
                                    +2*L12*L13^2*L23*L24^2*L25*L45
                                    -3*L13^3*L15*L24^2*L25*L45
                                    +6*L12*L13^2*L15*L24^2*L25*L45
                                    -3*L12^2*L13*L15*L24^2*L25*L45
                                    +3*L13^4*L24^2*L25*L45
                                    -6*L12*L13^3*L24^2*L25*L45
                                    +3*L12^2*L13^2*L24^2*L25*L45
                                    +2*L14*L15*L23^3*L24*L25*L45
                                    +4*L13*L15*L23^3*L24*L25*L45
                                    +4*L13*L14*L23^3*L24*L25*L45
                                    +6*L13^2*L23^3*L24*L25*L45
                                    +14*L13*L14*L15*L23^2*L24*L25*L45
                                    -6*L12*L14*L15*L23^2*L24*L25*L45
                                    +6*L13^2*L15*L23^2*L24*L25*L45
                                    -6*L12*L13*L15*L23^2*L24*L25*L45
                                    +6*L13^2*L14*L23^2*L24*L25*L45
                                    -6*L12*L13*L14*L23^2*L24*L25*L45
                                    -10*L13^3*L23^2*L24*L25*L45
                                    +2*L12*L13^2*L23^2*L24*L25*L45
                                    +14*L13^2*L14*L15*L23*L24*L25*L45
                                    -20*L12*L13*L14*L15*L23*L24*L25*L45
                                    +6*L12^2*L14*L15*L23*L24*L25*L45
                                    -8*L13^3*L15*L23*L24*L25*L45
                                    +8*L12*L13^2*L15*L23*L24*L25*L45
                                    -8*L13^3*L14*L23*L24*L25*L45
                                    +8*L12*L13^2*L14*L23*L24*L25*L45
                                    +2*L13^4*L23*L24*L25*L45
                                    +4*L12*L13^3*L23*L24*L25*L45
                                    -6*L12^2*L13^2*L23*L24*L25*L45
                                    +2*L13^3*L14*L15*L24*L25*L45
                                    -6*L12*L13^2*L14*L15*L24*L25*L45
                                    +6*L12^2*L13*L14*L15*L24*L25*L45
                                    -2*L12^3*L14*L15*L24*L25*L45
                                    -2*L13^4*L15*L24*L25*L45
                                    +6*L12*L13^3*L15*L24*L25*L45
                                    -6*L12^2*L13^2*L15*L24*L25*L45
                                    +2*L12^3*L13*L15*L24*L25*L45
                                    -2*L13^4*L14*L24*L25*L45
                                    +6*L12*L13^3*L14*L24*L25*L45
                                    -6*L12^2*L13^2*L14*L24*L25*L45
                                    +2*L12^3*L13*L14*L24*L25*L45
                                    +2*L13^5*L24*L25*L45
                                    -6*L12*L13^4*L24*L25*L45
                                    +6*L12^2*L13^3*L24*L25*L45
                                    -2*L12^3*L13^2*L24*L25*L45
                                    -2*L14*L15*L23^4*L25*L45
                                    -L13*L15*L23^4*L25*L45
                                    -4*L13*L14*L23^4*L25*L45
                                    -L13^2*L23^4*L25*L45
                                    -3*L14^2*L15*L23^3*L25*L45
                                    -8*L13*L14*L15*L23^3*L25*L45
                                    +6*L12*L14*L15*L23^3*L25*L45
                                    +L13^2*L15*L23^3*L25*L45
                                    -5*L13*L14^2*L23^3*L25*L45
                                    +4*L13^2*L14*L23^3*L25*L45
                                    +6*L12*L13*L14*L23^3*L25*L45
                                    +3*L13^3*L23^3*L25*L45
                                    -4*L12*L13^2*L23^3*L25*L45
                                    -10*L13*L14^2*L15*L23^2*L25*L45
                                    +6*L12*L14^2*L15*L23^2*L25*L45
                                    +6*L13^2*L14*L15*L23^2*L25*L45
                                    +8*L12*L13*L14*L15*L23^2*L25*L45
                                    -6*L12^2*L14*L15*L23^2*L25*L45
                                    +L13^3*L15*L23^2*L25*L45
                                    -8*L12*L13^2*L15*L23^2*L25*L45
                                    +3*L12^2*L13*L15*L23^2*L25*L45
                                    +2*L13^2*L14^2*L23^2*L25*L45
                                    +2*L12*L13*L14^2*L23^2*L25*L45
                                    +4*L13^3*L14*L23^2*L25*L45
                                    -12*L12*L13^2*L14*L23^2*L25*L45
                                    -3*L13^4*L23^2*L25*L45
                                    +4*L12*L13^3*L23^2*L25*L45
                                    +3*L12^2*L13^2*L23^2*L25*L45
                                    -3*L13^2*L14^2*L15*L23*L25*L45
                                    +6*L12*L13*L14^2*L15*L23*L25*L45
                                    -3*L12^2*L14^2*L15*L23*L25*L45
                                    +4*L13^3*L14*L15*L23*L25*L45
                                    -6*L12*L13^2*L14*L15*L23*L25*L45
                                    +2*L12^3*L14*L15*L23*L25*L45
                                    -L13^4*L15*L23*L25*L45
                                    +3*L12^2*L13^2*L15*L23*L25*L45
                                    -2*L12^3*L13*L15*L23*L25*L45
                                    +3*L13^3*L14^2*L23*L25*L45
                                    -6*L12*L13^2*L14^2*L23*L25*L45
                                    +3*L12^2*L13*L14^2*L23*L25*L45
                                    -4*L13^4*L14*L23*L25*L45
                                    +6*L12*L13^3*L14*L23*L25*L45
                                    -2*L12^3*L13*L14*L23*L25*L45
                                    +L13^5*L23*L25*L45
                                    -3*L12^2*L13^3*L23*L25*L45
                                    +2*L12^3*L13^2*L23*L25*L45
                                    +3*L13*L15*L23^3*L24^2*L45
                                    +L13^2*L23^3*L24^2*L45
                                    +4*L13*L15^2*L23^2*L24^2*L45
                                    +2*L13^2*L15*L23^2*L24^2*L45
                                    -6*L12*L13*L15*L23^2*L24^2*L45
                                    -2*L13^3*L23^2*L24^2*L45
                                    +2*L12*L13^2*L23^2*L24^2*L45
                                    +4*L13^2*L15^2*L23*L24^2*L45
                                    -4*L12*L13*L15^2*L23*L24^2*L45
                                    -5*L13^3*L15*L23*L24^2*L45
                                    +2*L12*L13^2*L15*L23*L24^2*L45
                                    +3*L12^2*L13*L15*L23*L24^2*L45
                                    +L13^4*L23*L24^2*L45
                                    +2*L12*L13^3*L23*L24^2*L45
                                    -3*L12^2*L13^2*L23*L24^2*L45
                                    -2*L14*L15*L23^4*L24*L45
                                    -4*L13*L15*L23^4*L24*L45
                                    -L13*L14*L23^4*L24*L45-L13^2*L23^4*L24*L45
                                    -3*L14*L15^2*L23^3*L24*L45
                                    -5*L13*L15^2*L23^3*L24*L45
                                    -8*L13*L14*L15*L23^3*L24*L45
                                    +6*L12*L14*L15*L23^3*L24*L45
                                    +4*L13^2*L15*L23^3*L24*L45
                                    +6*L12*L13*L15*L23^3*L24*L45
                                    +L13^2*L14*L23^3*L24*L45
                                    +3*L13^3*L23^3*L24*L45
                                    -4*L12*L13^2*L23^3*L24*L45
                                    -10*L13*L14*L15^2*L23^2*L24*L45
                                    +6*L12*L14*L15^2*L23^2*L24*L45
                                    +2*L13^2*L15^2*L23^2*L24*L45
                                    +2*L12*L13*L15^2*L23^2*L24*L45
                                    +6*L13^2*L14*L15*L23^2*L24*L45
                                    +8*L12*L13*L14*L15*L23^2*L24*L45
                                    -6*L12^2*L14*L15*L23^2*L24*L45
                                    +4*L13^3*L15*L23^2*L24*L45
                                    -12*L12*L13^2*L15*L23^2*L24*L45
                                    +L13^3*L14*L23^2*L24*L45
                                    -8*L12*L13^2*L14*L23^2*L24*L45
                                    +3*L12^2*L13*L14*L23^2*L24*L45
                                    -3*L13^4*L23^2*L24*L45
                                    +4*L12*L13^3*L23^2*L24*L45
                                    +3*L12^2*L13^2*L23^2*L24*L45
                                    -3*L13^2*L14*L15^2*L23*L24*L45
                                    +6*L12*L13*L14*L15^2*L23*L24*L45
                                    -3*L12^2*L14*L15^2*L23*L24*L45
                                    +3*L13^3*L15^2*L23*L24*L45
                                    -6*L12*L13^2*L15^2*L23*L24*L45
                                    +3*L12^2*L13*L15^2*L23*L24*L45
                                    +4*L13^3*L14*L15*L23*L24*L45
                                    -6*L12*L13^2*L14*L15*L23*L24*L45
                                    +2*L12^3*L14*L15*L23*L24*L45
                                    -4*L13^4*L15*L23*L24*L45
                                    +6*L12*L13^3*L15*L23*L24*L45
                                    -2*L12^3*L13*L15*L23*L24*L45
                                    -L13^4*L14*L23*L24*L45
                                    +3*L12^2*L13^2*L14*L23*L24*L45
                                    -2*L12^3*L13*L14*L23*L24*L45
                                    +L13^5*L23*L24*L45
                                    -3*L12^2*L13^3*L23*L24*L45
                                    +2*L12^3*L13^2*L23*L24*L45
                                    +2*L14*L15*L23^5*L45+L13*L15*L23^5*L45
                                    +L13*L14*L23^5*L45+3*L14*L15^2*L23^4*L45
                                    +L13*L15^2*L23^4*L45+3*L14^2*L15*L23^4*L45
                                    +2*L13*L14*L15*L23^4*L45
                                    -6*L12*L14*L15*L23^4*L45
                                    -3*L13^2*L15*L23^4*L45+L13*L14^2*L23^4*L45
                                    -3*L13^2*L14*L23^4*L45
                                    +2*L12*L13^2*L23^4*L45
                                    +4*L14^2*L15^2*L23^3*L45
                                    +2*L13*L14*L15^2*L23^3*L45
                                    -6*L12*L14*L15^2*L23^3*L45
                                    -2*L13^2*L15^2*L23^3*L45
                                    +2*L12*L13*L15^2*L23^3*L45
                                    +2*L13*L14^2*L15*L23^3*L45
                                    -6*L12*L14^2*L15*L23^3*L45
                                    -10*L13^2*L14*L15*L23^3*L45
                                    +4*L12*L13*L14*L15*L23^3*L45
                                    +6*L12^2*L14*L15*L23^3*L45
                                    +3*L13^3*L15*L23^3*L45
                                    +4*L12*L13^2*L15*L23^3*L45
                                    -3*L12^2*L13*L15*L23^3*L45
                                    -2*L13^2*L14^2*L23^3*L45
                                    +2*L12*L13*L14^2*L23^3*L45
                                    +3*L13^3*L14*L23^3*L45
                                    +4*L12*L13^2*L14*L23^3*L45
                                    -3*L12^2*L13*L14*L23^3*L45
                                    -4*L12*L13^3*L23^3*L45
                                    +4*L13*L14^2*L15^2*L23^2*L45
                                    -4*L12*L14^2*L15^2*L23^2*L45
                                    -5*L13^2*L14*L15^2*L23^2*L45
                                    +2*L12*L13*L14*L15^2*L23^2*L45
                                    +3*L12^2*L14*L15^2*L23^2*L45
                                    +L13^3*L15^2*L23^2*L45
                                    +2*L12*L13^2*L15^2*L23^2*L45
                                    -3*L12^2*L13*L15^2*L23^2*L45
                                    -5*L13^2*L14^2*L15*L23^2*L45
                                    +2*L12*L13*L14^2*L15*L23^2*L45
                                    +3*L12^2*L14^2*L15*L23^2*L45
                                    +6*L13^3*L14*L15*L23^2*L45
                                    +2*L12*L13^2*L14*L15*L23^2*L45
                                    -6*L12^2*L13*L14*L15*L23^2*L45
                                    -2*L12^3*L14*L15*L23^2*L45
                                    -L13^4*L15*L23^2*L45
                                    -4*L12*L13^3*L15*L23^2*L45
                                    +3*L12^2*L13^2*L15*L23^2*L45
                                    +2*L12^3*L13*L15*L23^2*L45
                                    +L13^3*L14^2*L23^2*L45
                                    +2*L12*L13^2*L14^2*L23^2*L45
                                    -3*L12^2*L13*L14^2*L23^2*L45
                                    -L13^4*L14*L23^2*L45
                                    -4*L12*L13^3*L14*L23^2*L45
                                    +3*L12^2*L13^2*L14*L23^2*L45
                                    +2*L12^3*L13*L14*L23^2*L45
                                    +2*L12*L13^4*L23^2*L45
                                    -2*L12^3*L13^2*L23^2*L45
                                    -L23^2*L24^2*L34*L35^3
                                    +2*L13*L23*L24^2*L34*L35^3
                                    +L12*L23*L24^2*L34*L35^3
                                    -L13^2*L24^2*L34*L35^3
                                    +L12*L13*L24^2*L34*L35^3
                                    +2*L14*L23^2*L24*L34*L35^3
                                    +L12*L23^2*L24*L34*L35^3
                                    -4*L13*L14*L23*L24*L34*L35^3
                                    -2*L12*L14*L23*L24*L34*L35^3
                                    -2*L12*L13*L23*L24*L34*L35^3
                                    +2*L13^2*L14*L24*L34*L35^3
                                    -2*L12*L13*L14*L24*L34*L35^3
                                    +L12*L13^2*L24*L34*L35^3
                                    -L12^3*L24*L34*L35^3-L14^2*L23^2*L34*L35^3
                                    +L12*L14*L23^2*L34*L35^3
                                    +2*L13*L14^2*L23*L34*L35^3
                                    +L12*L14^2*L23*L34*L35^3
                                    -2*L12*L13*L14*L23*L34*L35^3
                                    -L12^3*L23*L34*L35^3-L13^2*L14^2*L34*L35^3
                                    +L12*L13*L14^2*L34*L35^3
                                    +L12*L13^2*L14*L34*L35^3
                                    -L12^3*L14*L34*L35^3-L12^3*L13*L34*L35^3
                                    +L12^4*L34*L35^3-L13*L23*L24^3*L35^3
                                    +L13^2*L24^3*L35^3-L12*L13*L24^3*L35^3
                                    +L14*L23^2*L24^2*L35^3
                                    +L13*L23^2*L24^2*L35^3
                                    +L13*L14*L23*L24^2*L35^3
                                    -L12*L14*L23*L24^2*L35^3
                                    -2*L13^2*L23*L24^2*L35^3
                                    +2*L12*L13*L23*L24^2*L35^3
                                    -2*L13^2*L14*L24^2*L35^3
                                    +2*L12*L13*L14*L24^2*L35^3
                                    +L13^3*L24^2*L35^3-4*L12*L13^2*L24^2*L35^3
                                    +3*L12^2*L13*L24^2*L35^3
                                    -L14*L23^3*L24*L35^3
                                    -2*L14^2*L23^2*L24*L35^3
                                    +L13*L14*L23^2*L24*L35^3
                                    +2*L12*L14*L23^2*L24*L35^3
                                    -L12*L13*L23^2*L24*L35^3
                                    +L13*L14^2*L23*L24*L35^3
                                    +2*L12*L14^2*L23*L24*L35^3
                                    +L13^2*L14*L23*L24*L35^3
                                    -3*L12^2*L14*L23*L24*L35^3
                                    +2*L12*L13^2*L23*L24*L35^3
                                    -3*L12^2*L13*L23*L24*L35^3
                                    +L13^2*L14^2*L24*L35^3
                                    -L12*L13*L14^2*L24*L35^3
                                    -L13^3*L14*L24*L35^3
                                    +2*L12*L13^2*L14*L24*L35^3
                                    -3*L12^2*L13*L14*L24*L35^3
                                    +2*L12^3*L14*L24*L35^3-L12*L13^3*L24*L35^3
                                    +3*L12^2*L13^2*L24*L35^3
                                    -2*L12^3*L13*L24*L35^3+L14^2*L23^3*L35^3
                                    -L12*L14*L23^3*L35^3+L14^3*L23^2*L35^3
                                    -2*L13*L14^2*L23^2*L35^3
                                    -4*L12*L14^2*L23^2*L35^3
                                    +2*L12*L13*L14*L23^2*L35^3
                                    +3*L12^2*L14*L23^2*L35^3
                                    -L13*L14^3*L23*L35^3-L12*L14^3*L23*L35^3
                                    +L13^2*L14^2*L23*L35^3
                                    +2*L12*L13*L14^2*L23*L35^3
                                    +3*L12^2*L14^2*L23*L35^3
                                    -L12*L13^2*L14*L23*L35^3
                                    -3*L12^2*L13*L14*L23*L35^3
                                    -2*L12^3*L14*L23*L35^3
                                    +2*L12^3*L13*L23*L35^3
                                    +2*L23^2*L24*L25*L34^2*L35^2
                                    -4*L13*L23*L24*L25*L34^2*L35^2
                                    -2*L12*L23*L24*L25*L34^2*L35^2
                                    +2*L13^2*L24*L25*L34^2*L35^2
                                    -2*L12*L13*L24*L25*L34^2*L35^2
                                    -2*L14*L23^2*L25*L34^2*L35^2
                                    -L12*L23^2*L25*L34^2*L35^2
                                    +4*L13*L14*L23*L25*L34^2*L35^2
                                    +2*L12*L14*L23*L25*L34^2*L35^2
                                    +2*L12*L13*L23*L25*L34^2*L35^2
                                    -2*L13^2*L14*L25*L34^2*L35^2
                                    +2*L12*L13*L14*L25*L34^2*L35^2
                                    -L12*L13^2*L25*L34^2*L35^2
                                    +L12^3*L25*L34^2*L35^2
                                    -2*L15*L23^2*L24*L34^2*L35^2
                                    -L12*L23^2*L24*L34^2*L35^2
                                    +4*L13*L15*L23*L24*L34^2*L35^2
                                    +2*L12*L15*L23*L24*L34^2*L35^2
                                    +2*L12*L13*L23*L24*L34^2*L35^2
                                    -2*L13^2*L15*L24*L34^2*L35^2
                                    +2*L12*L13*L15*L24*L34^2*L35^2
                                    -L12*L13^2*L24*L34^2*L35^2
                                    +L12^3*L24*L34^2*L35^2
                                    +2*L14*L15*L23^2*L34^2*L35^2
                                    -L12*L15*L23^2*L34^2*L35^2
                                    -L12*L14*L23^2*L34^2*L35^2
                                    -4*L13*L14*L15*L23*L34^2*L35^2
                                    -2*L12*L14*L15*L23*L34^2*L35^2
                                    +2*L12*L13*L15*L23*L34^2*L35^2
                                    +2*L12*L13*L14*L23*L34^2*L35^2
                                    +2*L12^3*L23*L34^2*L35^2
                                    +2*L13^2*L14*L15*L34^2*L35^2
                                    -2*L12*L13*L14*L15*L34^2*L35^2
                                    -L12*L13^2*L15*L34^2*L35^2
                                    +L12^3*L15*L34^2*L35^2
                                    -L12*L13^2*L14*L34^2*L35^2
                                    +L12^3*L14*L34^2*L35^2
                                    +2*L12^3*L13*L34^2*L35^2
                                    -2*L12^4*L34^2*L35^2
                                    +L13*L23*L24^2*L25*L34*L35^2
                                    -L13^2*L24^2*L25*L34*L35^2
                                    +L12*L13*L24^2*L25*L34*L35^2
                                    -L23^3*L24*L25*L34*L35^2
                                    -4*L14*L23^2*L24*L25*L34*L35^2
                                    -L13*L23^2*L24*L25*L34*L35^2
                                    +3*L12*L23^2*L24*L25*L34*L35^2
                                    +6*L13*L14*L23*L24*L25*L34*L35^2
                                    +4*L12*L14*L23*L24*L25*L34*L35^2
                                    +5*L13^2*L23*L24*L25*L34*L35^2
                                    -3*L12^2*L23*L24*L25*L34*L35^2
                                    -2*L13^2*L14*L24*L25*L34*L35^2
                                    +2*L12*L13*L14*L24*L25*L34*L35^2
                                    -3*L13^3*L24*L25*L34*L35^2
                                    +5*L12*L13^2*L24*L25*L34*L35^2
                                    -3*L12^2*L13*L24*L25*L34*L35^2
                                    +L12^3*L24*L25*L34*L35^2
                                    +4*L14*L23^3*L25*L34*L35^2
                                    +L12*L23^3*L25*L34*L35^2
                                    +4*L14^2*L23^2*L25*L34*L35^2
                                    -8*L13*L14*L23^2*L25*L34*L35^2
                                    -4*L12*L14*L23^2*L25*L34*L35^2
                                    -L12*L13*L23^2*L25*L34*L35^2
                                    -3*L12^2*L23^2*L25*L34*L35^2
                                    -7*L13*L14^2*L23*L25*L34*L35^2
                                    -4*L12*L14^2*L23*L25*L34*L35^2
                                    +4*L13^2*L14*L23*L25*L34*L35^2
                                    +2*L12*L13*L14*L23*L25*L34*L35^2
                                    -L12*L13^2*L23*L25*L34*L35^2
                                    +3*L12^2*L13*L23*L25*L34*L35^2
                                    +3*L12^3*L23*L25*L34*L35^2
                                    +3*L13^2*L14^2*L25*L34*L35^2
                                    -3*L12*L13*L14^2*L25*L34*L35^2
                                    -6*L12*L13^2*L14*L25*L34*L35^2
                                    +6*L12^2*L13*L14*L25*L34*L35^2
                                    +L12*L13^3*L25*L34*L35^2
                                    -L12^4*L25*L34*L35^2+L23^3*L24^2*L34*L35^2
                                    +3*L15*L23^2*L24^2*L34*L35^2
                                    -3*L12*L23^2*L24^2*L34*L35^2
                                    -7*L13*L15*L23*L24^2*L34*L35^2
                                    -3*L12*L15*L23*L24^2*L34*L35^2
                                    -3*L13^2*L23*L24^2*L34*L35^2
                                    -2*L12*L13*L23*L24^2*L34*L35^2
                                    +3*L12^2*L23*L24^2*L34*L35^2
                                    +4*L13^2*L15*L24^2*L34*L35^2
                                    -4*L12*L13*L15*L24^2*L34*L35^2
                                    +2*L13^3*L24^2*L34*L35^2
                                    -L12*L13^2*L24^2*L34*L35^2
                                    -L12^3*L24^2*L34*L35^2
                                    -3*L14*L23^3*L24*L34*L35^2
                                    -L12*L23^3*L24*L34*L35^2
                                    -2*L14*L15*L23^2*L24*L34*L35^2
                                    +4*L13*L15*L23^2*L24*L34*L35^2
                                    -6*L12*L15*L23^2*L24*L34*L35^2
                                    +3*L13*L14*L23^2*L24*L34*L35^2
                                    +8*L12*L14*L23^2*L24*L34*L35^2
                                    +2*L12*L13*L23^2*L24*L34*L35^2
                                    +3*L12^2*L23^2*L24*L34*L35^2
                                    +6*L13*L14*L15*L23*L24*L34*L35^2
                                    +2*L12*L14*L15*L23*L24*L34*L35^2
                                    -8*L13^2*L15*L23*L24*L34*L35^2
                                    +2*L12*L13*L15*L23*L24*L34*L35^2
                                    +6*L12^2*L15*L23*L24*L34*L35^2
                                    +3*L13^2*L14*L23*L24*L34*L35^2
                                    -4*L12*L13*L14*L23*L24*L34*L35^2
                                    -3*L12^2*L14*L23*L24*L34*L35^2
                                    -L12*L13^2*L23*L24*L34*L35^2
                                    -3*L12^3*L23*L24*L34*L35^2
                                    -4*L13^2*L14*L15*L24*L34*L35^2
                                    +4*L12*L13*L14*L15*L24*L34*L35^2
                                    +4*L13^3*L15*L24*L34*L35^2
                                    -4*L12*L13^2*L15*L24*L34*L35^2
                                    -3*L13^3*L14*L24*L34*L35^2
                                    +8*L12*L13^2*L14*L24*L34*L35^2
                                    -3*L12^2*L13*L14*L24*L34*L35^2
                                    -2*L12^3*L14*L24*L34*L35^2
                                    -3*L12^2*L13^2*L24*L34*L35^2
                                    +2*L12^3*L13*L24*L34*L35^2
                                    +L12^4*L24*L34*L35^2
                                    -3*L14*L15*L23^3*L34*L35^2
                                    +L12*L15*L23^3*L34*L35^2
                                    +2*L14^2*L23^3*L34*L35^2
                                    -L14^2*L15*L23^2*L34*L35^2
                                    +5*L13*L14*L15*L23^2*L34*L35^2
                                    +5*L12*L14*L15*L23^2*L34*L35^2
                                    -L12*L13*L15*L23^2*L34*L35^2
                                    -3*L13*L14^2*L23^2*L34*L35^2
                                    -L12*L14^2*L23^2*L34*L35^2
                                    -L12*L13*L14*L23^2*L34*L35^2
                                    -3*L12^2*L14*L23^2*L34*L35^2
                                    +L13*L14^2*L15*L23*L34*L35^2
                                    +L12*L14^2*L15*L23*L34*L35^2
                                    -L13^2*L14*L15*L23*L34*L35^2
                                    -3*L12^2*L14*L15*L23*L34*L35^2
                                    -L12*L13^2*L15*L23*L34*L35^2
                                    +3*L12^2*L13*L15*L23*L34*L35^2
                                    -2*L12*L13*L14^2*L23*L34*L35^2
                                    +2*L12*L13^2*L14*L23*L34*L35^2
                                    +2*L12^3*L14*L23*L34*L35^2
                                    -2*L12^3*L13*L23*L34*L35^2
                                    -L13^3*L14*L15*L34*L35^2
                                    +3*L12*L13^2*L14*L15*L34*L35^2
                                    -3*L12^2*L13*L14*L15*L34*L35^2
                                    +L12^3*L14*L15*L34*L35^2
                                    +L12*L13^3*L15*L34*L35^2
                                    -3*L12^2*L13^2*L15*L34*L35^2
                                    +3*L12^3*L13*L15*L34*L35^2
                                    -L12^4*L15*L34*L35^2+L13^3*L14^2*L34*L35^2
                                    -3*L12*L13^2*L14^2*L34*L35^2
                                    +3*L12^2*L13*L14^2*L34*L35^2
                                    -L12^3*L14^2*L34*L35^2
                                    -L12*L13^3*L14*L34*L35^2
                                    +3*L12^2*L13^2*L14*L34*L35^2
                                    -3*L12^3*L13*L14*L34*L35^2
                                    +L12^4*L14*L34*L35^2
                                    -L13*L23*L24^3*L25*L35^2
                                    -L13^2*L24^3*L25*L35^2
                                    +L12*L13*L24^3*L25*L35^2
                                    +2*L13*L14*L23*L24^2*L25*L35^2
                                    -L13^2*L23*L24^2*L25*L35^2
                                    +3*L12*L13*L23*L24^2*L25*L35^2
                                    +4*L13^2*L14*L24^2*L25*L35^2
                                    -4*L12*L13*L14*L24^2*L25*L35^2
                                    +L13^3*L24^2*L25*L35^2
                                    +2*L12*L13^2*L24^2*L25*L35^2
                                    -3*L12^2*L13*L24^2*L25*L35^2
                                    +2*L14*L23^3*L24*L25*L35^2
                                    +L13*L23^3*L24*L25*L35^2
                                    +2*L14^2*L23^2*L24*L25*L35^2
                                    -6*L12*L14*L23^2*L24*L25*L35^2
                                    -L13^2*L23^2*L24*L25*L35^2
                                    -3*L12*L13*L23^2*L24*L25*L35^2
                                    -5*L13*L14^2*L23*L24*L25*L35^2
                                    -2*L12*L14^2*L23*L24*L25*L35^2
                                    -4*L13^2*L14*L23*L24*L25*L35^2
                                    +4*L12*L13*L14*L23*L24*L25*L35^2
                                    +6*L12^2*L14*L23*L24*L25*L35^2
                                    -L13^3*L23*L24*L25*L35^2
                                    +2*L12*L13^2*L23*L24*L25*L35^2
                                    -3*L13^2*L14^2*L24*L25*L35^2
                                    +3*L12*L13*L14^2*L24*L25*L35^2
                                    +2*L13^3*L14*L24*L25*L35^2
                                    -2*L12^3*L14*L24*L25*L35^2
                                    +L13^4*L24*L25*L35^2
                                    -3*L12*L13^3*L24*L25*L35^2
                                    +2*L12^3*L13*L24*L25*L35^2
                                    -2*L14*L23^4*L25*L35^2
                                    -5*L14^2*L23^3*L25*L35^2
                                    +4*L13*L14*L23^3*L25*L35^2
                                    +6*L12*L14*L23^3*L25*L35^2
                                    -L12*L13*L23^3*L25*L35^2
                                    -2*L14^3*L23^2*L25*L35^2
                                    +9*L13*L14^2*L23^2*L25*L35^2
                                    +8*L12*L14^2*L23^2*L25*L35^2
                                    -2*L13^2*L14*L23^2*L25*L35^2
                                    -12*L12*L13*L14*L23^2*L25*L35^2
                                    -6*L12^2*L14*L23^2*L25*L35^2
                                    +2*L12*L13^2*L23^2*L25*L35^2
                                    +3*L12^2*L13*L23^2*L25*L35^2
                                    +4*L13*L14^3*L23*L25*L35^2
                                    +2*L12*L14^3*L23*L25*L35^2
                                    -4*L13^2*L14^2*L23*L25*L35^2
                                    -11*L12*L13*L14^2*L23*L25*L35^2
                                    -3*L12^2*L14^2*L23*L25*L35^2
                                    +10*L12*L13^2*L14*L23*L25*L35^2
                                    +6*L12^2*L13*L14*L23*L25*L35^2
                                    +2*L12^3*L14*L23*L25*L35^2
                                    -L12*L13^3*L23*L25*L35^2
                                    -3*L12^2*L13^2*L23*L25*L35^2
                                    -2*L12^3*L13*L23*L25*L35^2
                                    +L13*L23^2*L24^3*L35^2
                                    +4*L13*L15*L23*L24^3*L35^2
                                    +L13^2*L23*L24^3*L35^2
                                    -L12*L13*L23*L24^3*L35^2
                                    -2*L13^2*L15*L24^3*L35^2
                                    +2*L12*L13*L15*L24^3*L35^2
                                    -2*L13^3*L24^3*L35^2
                                    +2*L12*L13^2*L24^3*L35^2
                                    -L14*L23^3*L24^2*L35^2
                                    -L13*L23^3*L24^2*L35^2
                                    -3*L14*L15*L23^2*L24^2*L35^2
                                    -4*L13*L15*L23^2*L24^2*L35^2
                                    -4*L13*L14*L23^2*L24^2*L35^2
                                    +3*L12*L14*L23^2*L24^2*L35^2
                                    -5*L13*L14*L15*L23*L24^2*L35^2
                                    +3*L12*L14*L15*L23*L24^2*L35^2
                                    +9*L13^2*L15*L23*L24^2*L35^2
                                    -11*L12*L13*L15*L23*L24^2*L35^2
                                    +2*L13^2*L14*L23*L24^2*L35^2
                                    +3*L12*L13*L14*L23*L24^2*L35^2
                                    -3*L12^2*L14*L23*L24^2*L35^2
                                    +3*L13^3*L23*L24^2*L35^2
                                    -L12*L13^2*L23*L24^2*L35^2
                                    +2*L13^2*L14*L15*L24^2*L35^2
                                    -2*L12*L13*L14*L15*L24^2*L35^2
                                    -5*L13^3*L15*L24^2*L35^2
                                    +8*L12*L13^2*L15*L24^2*L35^2
                                    -3*L12^2*L13*L15*L24^2*L35^2
                                    +3*L13^3*L14*L24^2*L35^2
                                    -7*L12*L13^2*L14*L24^2*L35^2
                                    +3*L12^2*L13*L14*L24^2*L35^2
                                    +L12^3*L14*L24^2*L35^2-2*L13^4*L24^2*L35^2
                                    +7*L12*L13^3*L24^2*L35^2
                                    -6*L12^2*L13^2*L24^2*L35^2
                                    +L12^3*L13*L24^2*L35^2+L14*L23^4*L24*L35^2
                                    +2*L14*L15*L23^3*L24*L35^2
                                    +3*L14^2*L23^3*L24*L35^2
                                    +2*L13*L14*L23^3*L24*L35^2
                                    -4*L12*L14*L23^3*L24*L35^2
                                    +L12*L13*L23^3*L24*L35^2
                                    +4*L14^2*L15*L23^2*L24*L35^2
                                    -4*L13*L14*L15*L23^2*L24*L35^2
                                    -2*L13^2*L15*L23^2*L24*L35^2
                                    +10*L12*L13*L15*L23^2*L24*L35^2
                                    +2*L13*L14^2*L23^2*L24*L35^2
                                    -7*L12*L14^2*L23^2*L24*L35^2
                                    -6*L13^2*L14*L23^2*L24*L35^2
                                    -2*L12*L13*L14*L23^2*L24*L35^2
                                    +6*L12^2*L14*L23^2*L24*L35^2
                                    -L12*L13^2*L23^2*L24*L35^2
                                    +2*L13*L14^2*L15*L23*L24*L35^2
                                    -4*L12*L14^2*L15*L23*L24*L35^2
                                    +4*L12*L13*L14*L15*L23*L24*L35^2
                                    +4*L13^3*L15*L23*L24*L35^2
                                    -12*L12*L13^2*L15*L23*L24*L35^2
                                    +6*L12^2*L13*L15*L23*L24*L35^2
                                    -4*L13^2*L14^2*L23*L24*L35^2
                                    +3*L12*L13*L14^2*L23*L24*L35^2
                                    +3*L12^2*L14^2*L23*L24*L35^2
                                    +2*L13^3*L14*L23*L24*L35^2
                                    -2*L12*L13^2*L14*L23*L24*L35^2
                                    -4*L12^3*L14*L23*L24*L35^2
                                    -L12*L13^3*L23*L24*L35^2
                                    +3*L12^2*L13^2*L23*L24*L35^2
                                    +2*L13^3*L14*L15*L24*L35^2
                                    -6*L12*L13^2*L14*L15*L24*L35^2
                                    +6*L12^2*L13*L14*L15*L24*L35^2
                                    -2*L12^3*L14*L15*L24*L35^2
                                    -2*L13^4*L15*L24*L35^2
                                    +6*L12*L13^3*L15*L24*L35^2
                                    -6*L12^2*L13^2*L15*L24*L35^2
                                    +2*L12^3*L13*L15*L24*L35^2
                                    -L13^3*L14^2*L24*L35^2
                                    +3*L12*L13^2*L14^2*L24*L35^2
                                    -3*L12^2*L13*L14^2*L24*L35^2
                                    +L12^3*L14^2*L24*L35^2+L13^4*L14*L24*L35^2
                                    -4*L12*L13^3*L14*L24*L35^2
                                    +6*L12^2*L13^2*L14*L24*L35^2
                                    -4*L12^3*L13*L14*L24*L35^2
                                    +L12^4*L14*L24*L35^2+L12*L13^4*L24*L35^2
                                    -3*L12^2*L13^3*L24*L35^2
                                    +3*L12^3*L13^2*L24*L35^2
                                    -L12^4*L13*L24*L35^2+L14*L15*L23^4*L35^2
                                    -2*L14^2*L23^4*L35^2+L12*L14*L23^4*L35^2
                                    +L14^2*L15*L23^3*L35^2
                                    -L13*L14*L15*L23^3*L35^2
                                    -3*L12*L14*L15*L23^3*L35^2
                                    -L12*L13*L15*L23^3*L35^2
                                    -2*L14^3*L23^3*L35^2
                                    +3*L13*L14^2*L23^3*L35^2
                                    +7*L12*L14^2*L23^3*L35^2
                                    -L12*L13*L14*L23^3*L35^2
                                    -3*L12^2*L14*L23^3*L35^2
                                    -L14^3*L15*L23^2*L35^2
                                    -L13*L14^2*L15*L23^2*L35^2
                                    +2*L12*L14^2*L15*L23^2*L35^2
                                    -L13^2*L14*L15*L23^2*L35^2
                                    +2*L12*L13*L14*L15*L23^2*L35^2
                                    +2*L12*L13^2*L15*L23^2*L35^2
                                    -3*L12^2*L13*L15*L23^2*L35^2
                                    +L13*L14^3*L23^2*L35^2
                                    +2*L12*L14^3*L23^2*L35^2
                                    -L12*L13*L14^2*L23^2*L35^2
                                    -6*L12^2*L14^2*L23^2*L35^2
                                    -L12*L13^2*L14*L23^2*L35^2
                                    +3*L12^2*L13*L14*L23^2*L35^2
                                    +3*L12^3*L14*L23^2*L35^2
                                    -L12^3*L13*L23^2*L35^2
                                    -L13*L14^3*L15*L23*L35^2
                                    +L12*L14^3*L15*L23*L35^2
                                    +3*L12*L13*L14^2*L15*L23*L35^2
                                    -3*L12^2*L14^2*L15*L23*L35^2
                                    +L13^3*L14*L15*L23*L35^2
                                    -3*L12*L13^2*L14*L15*L23*L35^2
                                    +2*L12^3*L14*L15*L23*L35^2
                                    -L12*L13^3*L15*L23*L35^2
                                    +3*L12^2*L13^2*L15*L23*L35^2
                                    -2*L12^3*L13*L15*L23*L35^2
                                    +L13^2*L14^3*L23*L35^2
                                    -L12*L13*L14^3*L23*L35^2
                                    -L13^3*L14^2*L23*L35^2
                                    +L12^3*L14^2*L23*L35^2
                                    +L12*L13^3*L14*L23*L35^2
                                    -L12^4*L14*L23*L35^2-L12^3*L13^2*L23*L35^2
                                    +L12^4*L13*L23*L35^2-L23^2*L25^2*L34^3*L35
                                    +2*L13*L23*L25^2*L34^3*L35
                                    +L12*L23*L25^2*L34^3*L35
                                    -L13^2*L25^2*L34^3*L35
                                    +L12*L13*L25^2*L34^3*L35
                                    +2*L15*L23^2*L25*L34^3*L35
                                    +L12*L23^2*L25*L34^3*L35
                                    -4*L13*L15*L23*L25*L34^3*L35
                                    -2*L12*L15*L23*L25*L34^3*L35
                                    -2*L12*L13*L23*L25*L34^3*L35
                                    +2*L13^2*L15*L25*L34^3*L35
                                    -2*L12*L13*L15*L25*L34^3*L35
                                    +L12*L13^2*L25*L34^3*L35
                                    -L12^3*L25*L34^3*L35-L15^2*L23^2*L34^3*L35
                                    +L12*L15*L23^2*L34^3*L35
                                    +2*L13*L15^2*L23*L34^3*L35
                                    +L12*L15^2*L23*L34^3*L35
                                    -2*L12*L13*L15*L23*L34^3*L35
                                    -L12^3*L23*L34^3*L35-L13^2*L15^2*L34^3*L35
                                    +L12*L13*L15^2*L34^3*L35
                                    +L12*L13^2*L15*L34^3*L35
                                    -L12^3*L15*L34^3*L35-L12^3*L13*L34^3*L35
                                    +L12^4*L34^3*L35
                                    +L13*L23*L24*L25^2*L34^2*L35
                                    -L13^2*L24*L25^2*L34^2*L35
                                    +L12*L13*L24*L25^2*L34^2*L35
                                    +L23^3*L25^2*L34^2*L35
                                    +3*L14*L23^2*L25^2*L34^2*L35
                                    -3*L12*L23^2*L25^2*L34^2*L35
                                    -7*L13*L14*L23*L25^2*L34^2*L35
                                    -3*L12*L14*L23*L25^2*L34^2*L35
                                    -3*L13^2*L23*L25^2*L34^2*L35
                                    -2*L12*L13*L23*L25^2*L34^2*L35
                                    +3*L12^2*L23*L25^2*L34^2*L35
                                    +4*L13^2*L14*L25^2*L34^2*L35
                                    -4*L12*L13*L14*L25^2*L34^2*L35
                                    +2*L13^3*L25^2*L34^2*L35
                                    -L12*L13^2*L25^2*L34^2*L35
                                    -L12^3*L25^2*L34^2*L35
                                    -L23^3*L24*L25*L34^2*L35
                                    -4*L15*L23^2*L24*L25*L34^2*L35
                                    -L13*L23^2*L24*L25*L34^2*L35
                                    +3*L12*L23^2*L24*L25*L34^2*L35
                                    +6*L13*L15*L23*L24*L25*L34^2*L35
                                    +4*L12*L15*L23*L24*L25*L34^2*L35
                                    +5*L13^2*L23*L24*L25*L34^2*L35
                                    -3*L12^2*L23*L24*L25*L34^2*L35
                                    -2*L13^2*L15*L24*L25*L34^2*L35
                                    +2*L12*L13*L15*L24*L25*L34^2*L35
                                    -3*L13^3*L24*L25*L34^2*L35
                                    +5*L12*L13^2*L24*L25*L34^2*L35
                                    -3*L12^2*L13*L24*L25*L34^2*L35
                                    +L12^3*L24*L25*L34^2*L35
                                    -3*L15*L23^3*L25*L34^2*L35
                                    -L12*L23^3*L25*L34^2*L35
                                    -2*L14*L15*L23^2*L25*L34^2*L35
                                    +3*L13*L15*L23^2*L25*L34^2*L35
                                    +8*L12*L15*L23^2*L25*L34^2*L35
                                    +4*L13*L14*L23^2*L25*L34^2*L35
                                    -6*L12*L14*L23^2*L25*L34^2*L35
                                    +2*L12*L13*L23^2*L25*L34^2*L35
                                    +3*L12^2*L23^2*L25*L34^2*L35
                                    +6*L13*L14*L15*L23*L25*L34^2*L35
                                    +2*L12*L14*L15*L23*L25*L34^2*L35
                                    +3*L13^2*L15*L23*L25*L34^2*L35
                                    -4*L12*L13*L15*L23*L25*L34^2*L35
                                    -3*L12^2*L15*L23*L25*L34^2*L35
                                    -8*L13^2*L14*L23*L25*L34^2*L35
                                    +2*L12*L13*L14*L23*L25*L34^2*L35
                                    +6*L12^2*L14*L23*L25*L34^2*L35
                                    -L12*L13^2*L23*L25*L34^2*L35
                                    -3*L12^3*L23*L25*L34^2*L35
                                    -4*L13^2*L14*L15*L25*L34^2*L35
                                    +4*L12*L13*L14*L15*L25*L34^2*L35
                                    -3*L13^3*L15*L25*L34^2*L35
                                    +8*L12*L13^2*L15*L25*L34^2*L35
                                    -3*L12^2*L13*L15*L25*L34^2*L35
                                    -2*L12^3*L15*L25*L34^2*L35
                                    +4*L13^3*L14*L25*L34^2*L35
                                    -4*L12*L13^2*L14*L25*L34^2*L35
                                    -3*L12^2*L13^2*L25*L34^2*L35
                                    +2*L12^3*L13*L25*L34^2*L35
                                    +L12^4*L25*L34^2*L35
                                    +4*L15*L23^3*L24*L34^2*L35
                                    +L12*L23^3*L24*L34^2*L35
                                    +4*L15^2*L23^2*L24*L34^2*L35
                                    -8*L13*L15*L23^2*L24*L34^2*L35
                                    -4*L12*L15*L23^2*L24*L34^2*L35
                                    -L12*L13*L23^2*L24*L34^2*L35
                                    -3*L12^2*L23^2*L24*L34^2*L35
                                    -7*L13*L15^2*L23*L24*L34^2*L35
                                    -4*L12*L15^2*L23*L24*L34^2*L35
                                    +4*L13^2*L15*L23*L24*L34^2*L35
                                    +2*L12*L13*L15*L23*L24*L34^2*L35
                                    -L12*L13^2*L23*L24*L34^2*L35
                                    +3*L12^2*L13*L23*L24*L34^2*L35
                                    +3*L12^3*L23*L24*L34^2*L35
                                    +3*L13^2*L15^2*L24*L34^2*L35
                                    -3*L12*L13*L15^2*L24*L34^2*L35
                                    -6*L12*L13^2*L15*L24*L34^2*L35
                                    +6*L12^2*L13*L15*L24*L34^2*L35
                                    +L12*L13^3*L24*L34^2*L35
                                    -L12^4*L24*L34^2*L35
                                    +2*L15^2*L23^3*L34^2*L35
                                    -3*L14*L15*L23^3*L34^2*L35
                                    +L12*L14*L23^3*L34^2*L35
                                    -L14*L15^2*L23^2*L34^2*L35
                                    -3*L13*L15^2*L23^2*L34^2*L35
                                    -L12*L15^2*L23^2*L34^2*L35
                                    +5*L13*L14*L15*L23^2*L34^2*L35
                                    +5*L12*L14*L15*L23^2*L34^2*L35
                                    -L12*L13*L15*L23^2*L34^2*L35
                                    -3*L12^2*L15*L23^2*L34^2*L35
                                    -L12*L13*L14*L23^2*L34^2*L35
                                    +L13*L14*L15^2*L23*L34^2*L35
                                    +L12*L14*L15^2*L23*L34^2*L35
                                    -2*L12*L13*L15^2*L23*L34^2*L35
                                    -L13^2*L14*L15*L23*L34^2*L35
                                    -3*L12^2*L14*L15*L23*L34^2*L35
                                    +2*L12*L13^2*L15*L23*L34^2*L35
                                    +2*L12^3*L15*L23*L34^2*L35
                                    -L12*L13^2*L14*L23*L34^2*L35
                                    +3*L12^2*L13*L14*L23*L34^2*L35
                                    -2*L12^3*L13*L23*L34^2*L35
                                    +L13^3*L15^2*L34^2*L35
                                    -3*L12*L13^2*L15^2*L34^2*L35
                                    +3*L12^2*L13*L15^2*L34^2*L35
                                    -L12^3*L15^2*L34^2*L35
                                    -L13^3*L14*L15*L34^2*L35
                                    +3*L12*L13^2*L14*L15*L34^2*L35
                                    -3*L12^2*L13*L14*L15*L34^2*L35
                                    +L12^3*L14*L15*L34^2*L35
                                    -L12*L13^3*L15*L34^2*L35
                                    +3*L12^2*L13^2*L15*L34^2*L35
                                    -3*L12^3*L13*L15*L34^2*L35
                                    +L12^4*L15*L34^2*L35
                                    +L12*L13^3*L14*L34^2*L35
                                    -3*L12^2*L13^2*L14*L34^2*L35
                                    +3*L12^3*L13*L14*L34^2*L35
                                    -L12^4*L14*L34^2*L35
                                    +2*L13*L23*L24^2*L25^2*L34*L35
                                    +2*L13^2*L24^2*L25^2*L34*L35
                                    -2*L12*L13*L24^2*L25^2*L34*L35
                                    -L13*L23^2*L24*L25^2*L34*L35
                                    -6*L13*L14*L23*L24*L25^2*L34*L35
                                    -2*L12*L13*L23*L24*L25^2*L34*L35
                                    -2*L13^2*L14*L24*L25^2*L34*L35
                                    +2*L12*L13*L14*L24*L25^2*L34*L35
                                    +L13^3*L24*L25^2*L34*L35
                                    -4*L12*L13^2*L24*L25^2*L34*L35
                                    +3*L12^2*L13*L24*L25^2*L34*L35
                                    -L14*L23^3*L25^2*L34*L35
                                    -2*L13*L23^3*L25^2*L34*L35
                                    -3*L14^2*L23^2*L25^2*L34*L35
                                    +6*L13*L14*L23^2*L25^2*L34*L35
                                    +3*L12*L14*L23^2*L25^2*L34*L35
                                    +3*L13^2*L23^2*L25^2*L34*L35
                                    +7*L12*L13*L23^2*L25^2*L34*L35
                                    +10*L13*L14^2*L23*L25^2*L34*L35
                                    +3*L12*L14^2*L23*L25^2*L34*L35
                                    -7*L13^2*L14*L23*L25^2*L34*L35
                                    +4*L12*L13*L14*L23*L25^2*L34*L35
                                    -3*L12^2*L14*L23*L25^2*L34*L35
                                    -L12*L13^2*L23*L25^2*L34*L35
                                    -6*L12^2*L13*L23*L25^2*L34*L35
                                    -3*L13^2*L14^2*L25^2*L34*L35
                                    +3*L12*L13*L14^2*L25^2*L34*L35
                                    +2*L13^3*L14*L25^2*L34*L35
                                    +3*L12*L13^2*L14*L25^2*L34*L35
                                    -6*L12^2*L13*L14*L25^2*L34*L35
                                    +L12^3*L14*L25^2*L34*L35
                                    -L13^4*L25^2*L34*L35
                                    +L12^3*L13*L25^2*L34*L35
                                    -L13*L23^2*L24^2*L25*L34*L35
                                    -6*L13*L15*L23*L24^2*L25*L34*L35
                                    -2*L12*L13*L23*L24^2*L25*L34*L35
                                    -2*L13^2*L15*L24^2*L25*L34*L35
                                    +2*L12*L13*L15*L24^2*L25*L34*L35
                                    +L13^3*L24^2*L25*L34*L35
                                    -4*L12*L13^2*L24^2*L25*L34*L35
                                    +3*L12^2*L13*L24^2*L25*L34*L35
                                    +4*L13*L23^3*L24*L25*L34*L35
                                    +8*L14*L15*L23^2*L24*L25*L34*L35
                                    +2*L13*L15*L23^2*L24*L25*L34*L35
                                    +2*L13*L14*L23^2*L24*L25*L34*L35
                                    -4*L13^2*L23^2*L24*L25*L34*L35
                                    -8*L12*L13*L23^2*L24*L25*L34*L35
                                    -8*L12*L14*L15*L23*L24*L25*L34*L35
                                    -4*L13^3*L23*L24*L25*L34*L35
                                    +12*L12^2*L13*L23*L24*L25*L34*L35
                                    +8*L13^2*L14*L15*L24*L25*L34*L35
                                    -8*L12*L13*L14*L15*L24*L25*L34*L35
                                    -2*L13^3*L15*L24*L25*L34*L35
                                    -4*L12*L13^2*L15*L24*L25*L34*L35
                                    +6*L12^2*L13*L15*L24*L25*L34*L35
                                    -2*L13^3*L14*L24*L25*L34*L35
                                    -4*L12*L13^2*L14*L24*L25*L34*L35
                                    +6*L12^2*L13*L14*L24*L25*L34*L35
                                    +4*L13^4*L24*L25*L34*L35
                                    -8*L12*L13^3*L24*L25*L34*L35
                                    +12*L12^2*L13^2*L24*L25*L34*L35
                                    -8*L12^3*L13*L24*L25*L34*L35
                                    +L15*L23^4*L25*L34*L35
                                    -2*L14*L15*L23^3*L25*L34*L35
                                    +2*L13*L15*L23^3*L25*L34*L35
                                    -4*L12*L15*L23^3*L25*L34*L35
                                    +2*L14^2*L23^3*L25*L34*L35
                                    -8*L13*L14*L23^3*L25*L34*L35
                                    +2*L12*L14*L23^3*L25*L34*L35
                                    -2*L14^2*L15*L23^2*L25*L34*L35
                                    -4*L12*L14*L15*L23^2*L25*L34*L35
                                    -6*L13^2*L15*L23^2*L25*L34*L35
                                    -2*L12*L13*L15*L23^2*L25*L34*L35
                                    +6*L12^2*L15*L23^2*L25*L34*L35
                                    -7*L13*L14^2*L23^2*L25*L34*L35
                                    +3*L12*L14^2*L23^2*L25*L34*L35
                                    +16*L13^2*L14*L23^2*L25*L34*L35
                                    +6*L12*L13*L14*L23^2*L25*L34*L35
                                    -6*L12^2*L14*L23^2*L25*L34*L35
                                    -L12*L13^2*L23^2*L25*L34*L35
                                    -3*L12^2*L13*L23^2*L25*L34*L35
                                    -6*L13*L14^2*L15*L23*L25*L34*L35
                                    +2*L12*L14^2*L15*L23*L25*L34*L35
                                    +2*L13^2*L14*L15*L23*L25*L34*L35
                                    +6*L12^2*L14*L15*L23*L25*L34*L35
                                    +2*L13^3*L15*L23*L25*L34*L35
                                    -2*L12*L13^2*L15*L23*L25*L34*L35
                                    -4*L12^3*L15*L23*L25*L34*L35
                                    +6*L13^2*L14^2*L23*L25*L34*L35
                                    +4*L12*L13*L14^2*L23*L25*L34*L35
                                    -6*L12^2*L14^2*L23*L25*L34*L35
                                    -8*L13^3*L14*L23*L25*L34*L35
                                    +6*L12*L13^2*L14*L23*L25*L34*L35
                                    -12*L12^2*L13*L14*L23*L25*L34*L35
                                    +6*L12^3*L14*L23*L25*L34*L35
                                    +2*L12*L13^3*L23*L25*L34*L35
                                    +2*L12^3*L13*L23*L25*L34*L35
                                    +L13^4*L15*L25*L34*L35
                                    -4*L12*L13^3*L15*L25*L34*L35
                                    +6*L12^2*L13^2*L15*L25*L34*L35
                                    -4*L12^3*L13*L15*L25*L34*L35
                                    +L12^4*L15*L25*L34*L35
                                    -L13^3*L14^2*L25*L34*L35
                                    +3*L12*L13^2*L14^2*L25*L34*L35
                                    -3*L12^2*L13*L14^2*L25*L34*L35
                                    +L12^3*L14^2*L25*L34*L35
                                    +2*L12*L13^3*L14*L25*L34*L35
                                    -6*L12^2*L13^2*L14*L25*L34*L35
                                    +6*L12^3*L13*L14*L25*L34*L35
                                    -2*L12^4*L14*L25*L34*L35
                                    -L12*L13^4*L25*L34*L35
                                    +3*L12^2*L13^3*L25*L34*L35
                                    -3*L12^3*L13^2*L25*L34*L35
                                    +L12^4*L13*L25*L34*L35
                                    -L15*L23^3*L24^2*L34*L35
                                    -2*L13*L23^3*L24^2*L34*L35
                                    -3*L15^2*L23^2*L24^2*L34*L35
                                    +6*L13*L15*L23^2*L24^2*L34*L35
                                    +3*L12*L15*L23^2*L24^2*L34*L35
                                    +3*L13^2*L23^2*L24^2*L34*L35
                                    +7*L12*L13*L23^2*L24^2*L34*L35
                                    +10*L13*L15^2*L23*L24^2*L34*L35
                                    +3*L12*L15^2*L23*L24^2*L34*L35
                                    -7*L13^2*L15*L23*L24^2*L34*L35
                                    +4*L12*L13*L15*L23*L24^2*L34*L35
                                    -3*L12^2*L15*L23*L24^2*L34*L35
                                    -L12*L13^2*L23*L24^2*L34*L35
                                    -6*L12^2*L13*L23*L24^2*L34*L35
                                    -3*L13^2*L15^2*L24^2*L34*L35
                                    +3*L12*L13*L15^2*L24^2*L34*L35
                                    +2*L13^3*L15*L24^2*L34*L35
                                    +3*L12*L13^2*L15*L24^2*L34*L35
                                    -6*L12^2*L13*L15*L24^2*L34*L35
                                    +L12^3*L15*L24^2*L34*L35
                                    -L13^4*L24^2*L34*L35
                                    +L12^3*L13*L24^2*L34*L35
                                    +L14*L23^4*L24*L34*L35
                                    +2*L15^2*L23^3*L24*L34*L35
                                    -2*L14*L15*L23^3*L24*L34*L35
                                    -8*L13*L15*L23^3*L24*L34*L35
                                    +2*L12*L15*L23^3*L24*L34*L35
                                    +2*L13*L14*L23^3*L24*L34*L35
                                    -4*L12*L14*L23^3*L24*L34*L35
                                    -2*L14*L15^2*L23^2*L24*L34*L35
                                    -7*L13*L15^2*L23^2*L24*L34*L35
                                    +3*L12*L15^2*L23^2*L24*L34*L35
                                    -4*L12*L14*L15*L23^2*L24*L34*L35
                                    +16*L13^2*L15*L23^2*L24*L34*L35
                                    +6*L12*L13*L15*L23^2*L24*L34*L35
                                    -6*L12^2*L15*L23^2*L24*L34*L35
                                    -6*L13^2*L14*L23^2*L24*L34*L35
                                    -2*L12*L13*L14*L23^2*L24*L34*L35
                                    +6*L12^2*L14*L23^2*L24*L34*L35
                                    -L12*L13^2*L23^2*L24*L34*L35
                                    -3*L12^2*L13*L23^2*L24*L34*L35
                                    -6*L13*L14*L15^2*L23*L24*L34*L35
                                    +2*L12*L14*L15^2*L23*L24*L34*L35
                                    +6*L13^2*L15^2*L23*L24*L34*L35
                                    +4*L12*L13*L15^2*L23*L24*L34*L35
                                    -6*L12^2*L15^2*L23*L24*L34*L35
                                    +2*L13^2*L14*L15*L23*L24*L34*L35
                                    +6*L12^2*L14*L15*L23*L24*L34*L35
                                    -8*L13^3*L15*L23*L24*L34*L35
                                    +6*L12*L13^2*L15*L23*L24*L34*L35
                                    -12*L12^2*L13*L15*L23*L24*L34*L35
                                    +6*L12^3*L15*L23*L24*L34*L35
                                    +2*L13^3*L14*L23*L24*L34*L35
                                    -2*L12*L13^2*L14*L23*L24*L34*L35
                                    -4*L12^3*L14*L23*L24*L34*L35
                                    +2*L12*L13^3*L23*L24*L34*L35
                                    +2*L12^3*L13*L23*L24*L34*L35
                                    -L13^3*L15^2*L24*L34*L35
                                    +3*L12*L13^2*L15^2*L24*L34*L35
                                    -3*L12^2*L13*L15^2*L24*L34*L35
                                    +L12^3*L15^2*L24*L34*L35
                                    +2*L12*L13^3*L15*L24*L34*L35
                                    -6*L12^2*L13^2*L15*L24*L34*L35
                                    +6*L12^3*L13*L15*L24*L34*L35
                                    -2*L12^4*L15*L24*L34*L35
                                    +L13^4*L14*L24*L34*L35
                                    -4*L12*L13^3*L14*L24*L34*L35
                                    +6*L12^2*L13^2*L14*L24*L34*L35
                                    -4*L12^3*L13*L14*L24*L34*L35
                                    +L12^4*L14*L24*L34*L35
                                    -L12*L13^4*L24*L34*L35
                                    +3*L12^2*L13^3*L24*L34*L35
                                    -3*L12^3*L13^2*L24*L34*L35
                                    +L12^4*L13*L24*L34*L35-L15^2*L23^4*L34*L35
                                    +4*L14*L15*L23^4*L34*L35
                                    -L12*L15*L23^4*L34*L35-L14^2*L23^4*L34*L35
                                    -L12*L14*L23^4*L34*L35
                                    +L14*L15^2*L23^3*L34*L35
                                    +L14^2*L15*L23^3*L34*L35
                                    -4*L13*L14*L15*L23^3*L34*L35
                                    -8*L12*L14*L15*L23^3*L34*L35
                                    +2*L12*L13*L15*L23^3*L34*L35
                                    +3*L12^2*L15*L23^3*L34*L35
                                    +2*L12*L13*L14*L23^3*L34*L35
                                    +3*L12^2*L14*L23^3*L34*L35
                                    +2*L14^2*L15^2*L23^2*L34*L35
                                    -4*L12*L14*L15^2*L23^2*L34*L35
                                    +3*L13^2*L15^2*L23^2*L34*L35
                                    -L12*L13*L15^2*L23^2*L34*L35
                                    -4*L12*L14^2*L15*L23^2*L34*L35
                                    -4*L13^2*L14*L15*L23^2*L34*L35
                                    +12*L12^2*L14*L15*L23^2*L34*L35
                                    -L12*L13^2*L15*L23^2*L34*L35
                                    -3*L12^3*L15*L23^2*L34*L35
                                    +3*L13^2*L14^2*L23^2*L34*L35
                                    -L12*L13*L14^2*L23^2*L34*L35
                                    -L12*L13^2*L14*L23^2*L34*L35
                                    -3*L12^3*L14*L23^2*L34*L35
                                    +2*L12^3*L13*L23^2*L34*L35
                                    +2*L13*L14^2*L15^2*L23*L34*L35
                                    -2*L12*L14^2*L15^2*L23*L34*L35
                                    -L13^2*L14*L15^2*L23*L34*L35
                                    -2*L12*L13*L14*L15^2*L23*L34*L35
                                    +3*L12^2*L14*L15^2*L23*L34*L35
                                    -2*L13^3*L15^2*L23*L34*L35
                                    +7*L12*L13^2*L15^2*L23*L34*L35
                                    -6*L12^2*L13*L15^2*L23*L34*L35
                                    +L12^3*L15^2*L23*L34*L35
                                    -L13^2*L14^2*L15*L23*L34*L35
                                    -2*L12*L13*L14^2*L15*L23*L34*L35
                                    +3*L12^2*L14^2*L15*L23*L34*L35
                                    +4*L13^3*L14*L15*L23*L34*L35
                                    -8*L12*L13^2*L14*L15*L23*L34*L35
                                    +12*L12^2*L13*L14*L15*L23*L34*L35
                                    -8*L12^3*L14*L15*L23*L34*L35
                                    -3*L12^2*L13^2*L15*L23*L34*L35
                                    +2*L12^3*L13*L15*L23*L34*L35
                                    +L12^4*L15*L23*L34*L35
                                    -2*L13^3*L14^2*L23*L34*L35
                                    +7*L12*L13^2*L14^2*L23*L34*L35
                                    -6*L12^2*L13*L14^2*L23*L34*L35
                                    +L12^3*L14^2*L23*L34*L35
                                    -3*L12^2*L13^2*L14*L23*L34*L35
                                    +2*L12^3*L13*L14*L23*L34*L35
                                    +L12^4*L14*L23*L34*L35
                                    +2*L12^3*L13^2*L23*L34*L35
                                    -2*L12^4*L13*L23*L34*L35
                                    -2*L13*L14*L23*L24^2*L25^2*L35
                                    -2*L13^2*L23*L24^2*L25^2*L35
                                    -2*L13^2*L14*L24^2*L25^2*L35
                                    +2*L12*L13*L14*L24^2*L25^2*L35
                                    -2*L13^3*L24^2*L25^2*L35
                                    +2*L12*L13^2*L24^2*L25^2*L35
                                    +7*L13*L14*L23^2*L24*L25^2*L35
                                    +L13^2*L23^2*L24*L25^2*L35
                                    +5*L13*L14^2*L23*L24*L25^2*L35
                                    +4*L13^2*L14*L23*L24*L25^2*L35
                                    -10*L12*L13*L14*L23*L24*L25^2*L35
                                    -L13^3*L23*L24*L25^2*L35
                                    +2*L12*L13^2*L23*L24*L25^2*L35
                                    +3*L13^2*L14^2*L24*L25^2*L35
                                    -3*L12*L13*L14^2*L24*L25^2*L35
                                    -3*L13^3*L14*L24*L25^2*L35
                                    +3*L12^2*L13*L14*L24*L25^2*L35
                                    +3*L12*L13^3*L24*L25^2*L35
                                    -3*L12^2*L13^2*L24*L25^2*L35
                                    -5*L13*L14*L23^3*L25^2*L35
                                    +L13^2*L23^3*L25^2*L35
                                    +L14^3*L23^2*L25^2*L35
                                    -12*L13*L14^2*L23^2*L25^2*L35
                                    +9*L13^2*L14*L23^2*L25^2*L35
                                    +8*L12*L13*L14*L23^2*L25^2*L35
                                    -2*L13^3*L23^2*L25^2*L35
                                    -4*L12*L13^2*L23^2*L25^2*L35
                                    -5*L13*L14^3*L23*L25^2*L35
                                    -L12*L14^3*L23*L25^2*L35
                                    +8*L13^2*L14^2*L23*L25^2*L35
                                    +10*L12*L13*L14^2*L23*L25^2*L35
                                    -4*L13^3*L14*L23*L25^2*L35
                                    -11*L12*L13^2*L14*L23*L25^2*L35
                                    -3*L12^2*L13*L14*L23*L25^2*L35
                                    +L13^4*L23*L25^2*L35
                                    +2*L12*L13^3*L23*L25^2*L35
                                    +3*L12^2*L13^2*L23*L25^2*L35
                                    +2*L13*L15*L23*L24^3*L25*L35
                                    +2*L13^2*L23*L24^3*L25*L35
                                    +2*L13^2*L15*L24^3*L25*L35
                                    -2*L12*L13*L15*L24^3*L25*L35
                                    +2*L13^3*L24^3*L25*L35
                                    -2*L12*L13^2*L24^3*L25*L35
                                    -3*L13*L15*L23^2*L24^2*L25*L35
                                    -2*L13*L14*L23^2*L24^2*L25*L35
                                    +L13^2*L23^2*L24^2*L25*L35
                                    +4*L13^2*L15*L23*L24^2*L25*L35
                                    -8*L13^2*L14*L23*L24^2*L25*L35
                                    +8*L12*L13*L14*L23*L24^2*L25*L35
                                    -4*L12*L13^2*L23*L24^2*L25*L35
                                    -4*L13^2*L14*L15*L24^2*L25*L35
                                    +4*L12*L13*L14*L15*L24^2*L25*L35
                                    +7*L13^3*L15*L24^2*L25*L35
                                    -10*L12*L13^2*L15*L24^2*L25*L35
                                    +3*L12^2*L13*L15*L24^2*L25*L35
                                    -2*L13^3*L14*L24^2*L25*L35
                                    +8*L12*L13^2*L14*L24^2*L25*L35
                                    -6*L12^2*L13*L14*L24^2*L25*L35
                                    -L13^4*L24^2*L25*L35
                                    -2*L12*L13^3*L24^2*L25*L35
                                    +3*L12^2*L13^2*L24^2*L25*L35
                                    -2*L14*L15*L23^3*L24*L25*L35
                                    +2*L13*L15*L23^3*L24*L25*L35
                                    +L14^2*L23^3*L24*L25*L35
                                    -2*L13*L14*L23^3*L24*L25*L35
                                    -3*L13^2*L23^3*L24*L25*L35
                                    -4*L14^2*L15*L23^2*L24*L25*L35
                                    -6*L13*L14*L15*L23^2*L24*L25*L35
                                    +6*L12*L14*L15*L23^2*L24*L25*L35
                                    -4*L13^2*L15*L23^2*L24*L25*L35
                                    +5*L13*L14^2*L23^2*L24*L25*L35
                                    -3*L12*L14^2*L23^2*L24*L25*L35
                                    -4*L12*L13*L14*L23^2*L24*L25*L35
                                    +5*L13^3*L23^2*L24*L25*L35
                                    +5*L12*L13^2*L23^2*L24*L25*L35
                                    +4*L12*L14^2*L15*L23*L24*L25*L35
                                    -6*L13^2*L14*L15*L23*L24*L25*L35
                                    +4*L12*L13*L14*L15*L23*L24*L25*L35
                                    -6*L12^2*L14*L15*L23*L24*L25*L35
                                    +4*L12*L13^2*L15*L23*L24*L25*L35
                                    +5*L13^2*L14^2*L23*L24*L25*L35
                                    -12*L12*L13*L14^2*L23*L24*L25*L35
                                    +3*L12^2*L14^2*L23*L24*L25*L35
                                    +2*L13^3*L14*L23*L24*L25*L35
                                    +6*L12^2*L13*L14*L23*L24*L25*L35
                                    -L13^4*L23*L24*L25*L35
                                    -3*L12^2*L13^2*L23*L24*L25*L35
                                    -2*L13^3*L14*L15*L24*L25*L35
                                    +6*L12*L13^2*L14*L15*L24*L25*L35
                                    -6*L12^2*L13*L14*L15*L24*L25*L35
                                    +2*L12^3*L14*L15*L24*L25*L35
                                    +2*L13^4*L15*L24*L25*L35
                                    -6*L12*L13^3*L15*L24*L25*L35
                                    +6*L12^2*L13^2*L15*L24*L25*L35
                                    -2*L12^3*L13*L15*L24*L25*L35
                                    +L13^3*L14^2*L24*L25*L35
                                    -3*L12*L13^2*L14^2*L24*L25*L35
                                    +3*L12^2*L13*L14^2*L24*L25*L35
                                    -L12^3*L14^2*L24*L25*L35-L13^5*L24*L25*L35
                                    +3*L12*L13^4*L24*L25*L35
                                    -3*L12^2*L13^3*L24*L25*L35
                                    +L12^3*L13^2*L24*L25*L35
                                    +2*L14*L15*L23^4*L25*L35
                                    -L13*L15*L23^4*L25*L35-L14^2*L23^4*L25*L35
                                    +4*L13*L14*L23^4*L25*L35
                                    +7*L14^2*L15*L23^3*L25*L35
                                    -6*L12*L14*L15*L23^3*L25*L35
                                    +L13^2*L15*L23^3*L25*L35
                                    +2*L12*L13*L15*L23^3*L25*L35
                                    -2*L14^3*L23^3*L25*L35
                                    +6*L13*L14^2*L23^3*L25*L35
                                    +3*L12*L14^2*L23^3*L25*L35
                                    -8*L13^2*L14*L23^3*L25*L35
                                    -4*L12*L13*L14*L23^3*L25*L35
                                    +L12*L13^2*L23^3*L25*L35
                                    +2*L14^3*L15*L23^2*L25*L35
                                    +4*L13*L14^2*L15*L23^2*L25*L35
                                    -10*L12*L14^2*L15*L23^2*L25*L35
                                    -4*L13^2*L14*L15*L23^2*L25*L35
                                    +4*L12*L13*L14*L15*L23^2*L25*L35
                                    +6*L12^2*L14*L15*L23^2*L25*L35
                                    +L13^3*L15*L23^2*L25*L35
                                    -3*L12^2*L13*L15*L23^2*L25*L35
                                    +2*L12*L14^3*L23^2*L25*L35
                                    -7*L13^2*L14^2*L23^2*L25*L35
                                    +4*L12*L13*L14^2*L23^2*L25*L35
                                    -3*L12^2*L14^2*L23^2*L25*L35
                                    +4*L13^3*L14*L23^2*L25*L35
                                    +2*L12*L13^2*L14*L23^2*L25*L35
                                    -2*L12*L13^3*L23^2*L25*L35
                                    +2*L13*L14^3*L15*L23*L25*L35
                                    -2*L12*L14^3*L15*L23*L25*L35
                                    -3*L13^2*L14^2*L15*L23*L25*L35
                                    +3*L12^2*L14^2*L15*L23*L25*L35
                                    +2*L13^3*L14*L15*L23*L25*L35
                                    -2*L12^3*L14*L15*L23*L25*L35
                                    -L13^4*L15*L23*L25*L35
                                    +2*L12*L13^3*L15*L23*L25*L35
                                    -3*L12^2*L13^2*L15*L23*L25*L35
                                    +2*L12^3*L13*L15*L23*L25*L35
                                    -2*L13^2*L14^3*L23*L25*L35
                                    +2*L12*L13*L14^3*L23*L25*L35
                                    +2*L13^3*L14^2*L23*L25*L35
                                    +3*L12*L13^2*L14^2*L23*L25*L35
                                    -6*L12^2*L13*L14^2*L23*L25*L35
                                    +L12^3*L14^2*L23*L25*L35
                                    -6*L12*L13^3*L14*L23*L25*L35
                                    +6*L12^2*L13^2*L14*L23*L25*L35
                                    +L12*L13^4*L23*L25*L35
                                    -L12^3*L13^2*L23*L25*L35
                                    -2*L13*L15*L23^2*L24^3*L35
                                    -2*L13^2*L23^2*L24^3*L35
                                    -5*L13*L15^2*L23*L24^3*L35
                                    +2*L12*L13*L15*L23*L24^3*L35
                                    +L13^3*L23*L24^3*L35
                                    +2*L12*L13^2*L23*L24^3*L35
                                    +L13^2*L15^2*L24^3*L35
                                    -L12*L13*L15^2*L24^3*L35
                                    -2*L13^3*L15*L24^3*L35
                                    +2*L12*L13^2*L15*L24^3*L35+L13^4*L24^3*L35
                                    -L12*L13^3*L24^3*L35
                                    +L14*L15*L23^3*L24^2*L35
                                    +2*L13*L15*L23^3*L24^2*L35
                                    +3*L13*L14*L23^3*L24^2*L35
                                    +2*L13^2*L23^3*L24^2*L35
                                    +3*L14*L15^2*L23^2*L24^2*L35
                                    +8*L13*L15^2*L23^2*L24^2*L35
                                    +5*L13*L14*L15*L23^2*L24^2*L35
                                    -3*L12*L14*L15*L23^2*L24^2*L35
                                    -7*L13^2*L15*L23^2*L24^2*L35
                                    +3*L12*L13*L15*L23^2*L24^2*L35
                                    +2*L13^2*L14*L23^2*L24^2*L35
                                    -7*L12*L13*L14*L23^2*L24^2*L35
                                    -3*L13^3*L23^2*L24^2*L35
                                    -L12*L13^2*L23^2*L24^2*L35
                                    +5*L13*L14*L15^2*L23*L24^2*L35
                                    -3*L12*L14*L15^2*L23*L24^2*L35
                                    -12*L13^2*L15^2*L23*L24^2*L35
                                    +10*L12*L13*L15^2*L23*L24^2*L35
                                    +5*L13^2*L14*L15*L23*L24^2*L35
                                    -12*L12*L13*L14*L15*L23*L24^2*L35
                                    +3*L12^2*L14*L15*L23*L24^2*L35
                                    +6*L13^3*L15*L23*L24^2*L35
                                    +4*L12*L13^2*L15*L23*L24^2*L35
                                    -6*L12^2*L13*L15*L23*L24^2*L35
                                    -4*L13^3*L14*L23*L24^2*L35
                                    +3*L12*L13^2*L14*L23*L24^2*L35
                                    +3*L12^2*L13*L14*L23*L24^2*L35
                                    -2*L12*L13^3*L23*L24^2*L35
                                    +L13^3*L14*L15*L24^2*L35
                                    -3*L12*L13^2*L14*L15*L24^2*L35
                                    +3*L12^2*L13*L14*L15*L24^2*L35
                                    -L12^3*L14*L15*L24^2*L35
                                    -L13^4*L15*L24^2*L35
                                    +3*L12*L13^3*L15*L24^2*L35
                                    -3*L12^2*L13^2*L15*L24^2*L35
                                    +L12^3*L13*L15*L24^2*L35
                                    -L13^4*L14*L24^2*L35
                                    +3*L12*L13^3*L14*L24^2*L35
                                    -3*L12^2*L13^2*L14*L24^2*L35
                                    +L12^3*L13*L14*L24^2*L35+L13^5*L24^2*L35
                                    -3*L12*L13^4*L24^2*L35
                                    +3*L12^2*L13^3*L24^2*L35
                                    -L12^3*L13^2*L24^2*L35-L14^2*L23^4*L24*L35
                                    -3*L13*L14*L23^4*L24*L35
                                    -3*L14*L15^2*L23^3*L24*L35
                                    -4*L13*L15^2*L23^3*L24*L35
                                    -2*L14^2*L15*L23^3*L24*L35
                                    +2*L13*L14*L15*L23^3*L24*L35
                                    +4*L13^2*L15*L23^3*L24*L35
                                    -6*L12*L13*L15*L23^3*L24*L35
                                    -4*L13*L14^2*L23^3*L24*L35
                                    +3*L12*L14^2*L23^3*L24*L35
                                    +3*L13^2*L14*L23^3*L24*L35
                                    +8*L12*L13*L14*L23^3*L24*L35
                                    -L12*L13^2*L23^3*L24*L35
                                    -2*L14^2*L15^2*L23^2*L24*L35
                                    +4*L13*L14*L15^2*L23^2*L24*L35
                                    +9*L13^2*L15^2*L23^2*L24*L35
                                    -11*L12*L13*L15^2*L23^2*L24*L35
                                    -8*L13*L14^2*L15*L23^2*L24*L35
                                    +8*L12*L14^2*L15*L23^2*L24*L35
                                    -8*L13^3*L15*L23^2*L24*L35
                                    +2*L12*L13^2*L15*L23^2*L24*L35
                                    +6*L12^2*L13*L15*L23^2*L24*L35
                                    +2*L13^2*L14^2*L23^2*L24*L35
                                    +3*L12*L13*L14^2*L23^2*L24*L35
                                    -3*L12^2*L14^2*L23^2*L24*L35
                                    +3*L13^3*L14*L23^2*L24*L35
                                    -4*L12*L13^2*L14*L23^2*L24*L35
                                    -3*L12^2*L13*L14*L23^2*L24*L35
                                    +2*L12*L13^3*L23^2*L24*L35
                                    -2*L13*L14^2*L15^2*L23*L24*L35
                                    +2*L12*L14^2*L15^2*L23*L24*L35
                                    +7*L13^2*L14*L15^2*L23*L24*L35
                                    -10*L12*L13*L14*L15^2*L23*L24*L35
                                    +3*L12^2*L14*L15^2*L23*L24*L35
                                    -5*L13^3*L15^2*L23*L24*L35
                                    +8*L12*L13^2*L15^2*L23*L24*L35
                                    -3*L12^2*L13*L15^2*L23*L24*L35
                                    -2*L13^2*L14^2*L15*L23*L24*L35
                                    +8*L12*L13*L14^2*L15*L23*L24*L35
                                    -6*L12^2*L14^2*L15*L23*L24*L35
                                    -2*L13^3*L14*L15*L23*L24*L35
                                    -4*L12*L13^2*L14*L15*L23*L24*L35
                                    +6*L12^2*L13*L14*L15*L23*L24*L35
                                    +4*L13^4*L15*L23*L24*L35
                                    -4*L12*L13^3*L15*L23*L24*L35
                                    +3*L13^3*L14^2*L23*L24*L35
                                    -7*L12*L13^2*L14^2*L23*L24*L35
                                    +3*L12^2*L13*L14^2*L23*L24*L35
                                    +L12^3*L14^2*L23*L24*L35
                                    -3*L13^4*L14*L23*L24*L35
                                    +8*L12*L13^3*L14*L23*L24*L35
                                    -3*L12^2*L13^2*L14*L23*L24*L35
                                    -2*L12^3*L13*L14*L23*L24*L35
                                    -L12*L13^4*L23*L24*L35
                                    +L12^3*L13^2*L23*L24*L35-L14*L15*L23^5*L35
                                    +L14^2*L23^5*L35+L13*L15^2*L23^4*L35
                                    -L14^2*L15*L23^4*L35-L13*L14*L15*L23^4*L35
                                    +3*L12*L14*L15*L23^4*L35
                                    +L12*L13*L15*L23^4*L35+L14^3*L23^4*L35
                                    -3*L12*L14^2*L23^4*L35
                                    -L12*L13*L14*L23^4*L35
                                    -2*L14^2*L15^2*L23^3*L35
                                    -L13*L14*L15^2*L23^3*L35
                                    +3*L12*L14*L15^2*L23^3*L35
                                    -2*L13^2*L15^2*L23^3*L35
                                    +2*L12*L13*L15^2*L23^3*L35
                                    +2*L14^3*L15*L23^3*L35
                                    -2*L12*L14^2*L15*L23^3*L35
                                    +5*L13^2*L14*L15*L23^3*L35
                                    -3*L12^2*L14*L15*L23^3*L35
                                    -2*L12*L13^2*L15*L23^3*L35
                                    +L13*L14^3*L23^3*L35-L12*L14^3*L23^3*L35
                                    -3*L13^2*L14^2*L23^3*L35
                                    -2*L12*L13*L14^2*L23^3*L35
                                    +3*L12^2*L14^2*L23^3*L35
                                    +2*L12*L13^2*L14*L23^3*L35
                                    -2*L13*L14^2*L15^2*L23^2*L35
                                    +2*L12*L14^2*L15^2*L23^2*L35
                                    +L13^2*L14*L15^2*L23^2*L35
                                    +2*L12*L13*L14*L15^2*L23^2*L35
                                    -3*L12^2*L14*L15^2*L23^2*L35
                                    +L13^3*L15^2*L23^2*L35
                                    -4*L12*L13^2*L15^2*L23^2*L35
                                    +3*L12^2*L13*L15^2*L23^2*L35
                                    +2*L13*L14^3*L15*L23^2*L35
                                    -2*L12*L14^3*L15*L23^2*L35
                                    +L13^2*L14^2*L15*L23^2*L35
                                    -4*L12*L13*L14^2*L15*L23^2*L35
                                    +3*L12^2*L14^2*L15*L23^2*L35
                                    -3*L13^3*L14*L15*L23^2*L35
                                    +5*L12*L13^2*L14*L15*L23^2*L35
                                    -3*L12^2*L13*L14*L15*L23^2*L35
                                    +L12^3*L14*L15*L23^2*L35
                                    +L12*L13^3*L15*L23^2*L35
                                    -L12^3*L13*L15*L23^2*L35
                                    -2*L13^2*L14^3*L23^2*L35
                                    +2*L12*L13*L14^3*L23^2*L35
                                    +2*L13^3*L14^2*L23^2*L35
                                    -L12*L13^2*L14^2*L23^2*L35
                                    -L12^3*L14^2*L23^2*L35
                                    -L12*L13^3*L14*L23^2*L35
                                    +L12^3*L13*L14*L23^2*L35
                                    -L13*L23*L25^3*L34^3+L13^2*L25^3*L34^3
                                    -L12*L13*L25^3*L34^3+L15*L23^2*L25^2*L34^3
                                    +L13*L23^2*L25^2*L34^3
                                    +L13*L15*L23*L25^2*L34^3
                                    -L12*L15*L23*L25^2*L34^3
                                    -2*L13^2*L23*L25^2*L34^3
                                    +2*L12*L13*L23*L25^2*L34^3
                                    -2*L13^2*L15*L25^2*L34^3
                                    +2*L12*L13*L15*L25^2*L34^3
                                    +L13^3*L25^2*L34^3-4*L12*L13^2*L25^2*L34^3
                                    +3*L12^2*L13*L25^2*L34^3
                                    -L15*L23^3*L25*L34^3
                                    -2*L15^2*L23^2*L25*L34^3
                                    +L13*L15*L23^2*L25*L34^3
                                    +2*L12*L15*L23^2*L25*L34^3
                                    -L12*L13*L23^2*L25*L34^3
                                    +L13*L15^2*L23*L25*L34^3
                                    +2*L12*L15^2*L23*L25*L34^3
                                    +L13^2*L15*L23*L25*L34^3
                                    -3*L12^2*L15*L23*L25*L34^3
                                    +2*L12*L13^2*L23*L25*L34^3
                                    -3*L12^2*L13*L23*L25*L34^3
                                    +L13^2*L15^2*L25*L34^3
                                    -L12*L13*L15^2*L25*L34^3
                                    -L13^3*L15*L25*L34^3
                                    +2*L12*L13^2*L15*L25*L34^3
                                    -3*L12^2*L13*L15*L25*L34^3
                                    +2*L12^3*L15*L25*L34^3-L12*L13^3*L25*L34^3
                                    +3*L12^2*L13^2*L25*L34^3
                                    -2*L12^3*L13*L25*L34^3+L15^2*L23^3*L34^3
                                    -L12*L15*L23^3*L34^3+L15^3*L23^2*L34^3
                                    -2*L13*L15^2*L23^2*L34^3
                                    -4*L12*L15^2*L23^2*L34^3
                                    +2*L12*L13*L15*L23^2*L34^3
                                    +3*L12^2*L15*L23^2*L34^3
                                    -L13*L15^3*L23*L34^3-L12*L15^3*L23*L34^3
                                    +L13^2*L15^2*L23*L34^3
                                    +2*L12*L13*L15^2*L23*L34^3
                                    +3*L12^2*L15^2*L23*L34^3
                                    -L12*L13^2*L15*L23*L34^3
                                    -3*L12^2*L13*L15*L23*L34^3
                                    -2*L12^3*L15*L23*L34^3
                                    +2*L12^3*L13*L23*L34^3
                                    -L13*L23*L24*L25^3*L34^2
                                    -L13^2*L24*L25^3*L34^2
                                    +L12*L13*L24*L25^3*L34^2
                                    +L13*L23^2*L25^3*L34^2
                                    +4*L13*L14*L23*L25^3*L34^2
                                    +L13^2*L23*L25^3*L34^2
                                    -L12*L13*L23*L25^3*L34^2
                                    -2*L13^2*L14*L25^3*L34^2
                                    +2*L12*L13*L14*L25^3*L34^2
                                    -2*L13^3*L25^3*L34^2
                                    +2*L12*L13^2*L25^3*L34^2
                                    +2*L13*L15*L23*L24*L25^2*L34^2
                                    -L13^2*L23*L24*L25^2*L34^2
                                    +3*L12*L13*L23*L24*L25^2*L34^2
                                    +4*L13^2*L15*L24*L25^2*L34^2
                                    -4*L12*L13*L15*L24*L25^2*L34^2
                                    +L13^3*L24*L25^2*L34^2
                                    +2*L12*L13^2*L24*L25^2*L34^2
                                    -3*L12^2*L13*L24*L25^2*L34^2
                                    -L15*L23^3*L25^2*L34^2
                                    -L13*L23^3*L25^2*L34^2
                                    -3*L14*L15*L23^2*L25^2*L34^2
                                    -4*L13*L15*L23^2*L25^2*L34^2
                                    +3*L12*L15*L23^2*L25^2*L34^2
                                    -4*L13*L14*L23^2*L25^2*L34^2
                                    -5*L13*L14*L15*L23*L25^2*L34^2
                                    +3*L12*L14*L15*L23*L25^2*L34^2
                                    +2*L13^2*L15*L23*L25^2*L34^2
                                    +3*L12*L13*L15*L23*L25^2*L34^2
                                    -3*L12^2*L15*L23*L25^2*L34^2
                                    +9*L13^2*L14*L23*L25^2*L34^2
                                    -11*L12*L13*L14*L23*L25^2*L34^2
                                    +3*L13^3*L23*L25^2*L34^2
                                    -L12*L13^2*L23*L25^2*L34^2
                                    +2*L13^2*L14*L15*L25^2*L34^2
                                    -2*L12*L13*L14*L15*L25^2*L34^2
                                    +3*L13^3*L15*L25^2*L34^2
                                    -7*L12*L13^2*L15*L25^2*L34^2
                                    +3*L12^2*L13*L15*L25^2*L34^2
                                    +L12^3*L15*L25^2*L34^2
                                    -5*L13^3*L14*L25^2*L34^2
                                    +8*L12*L13^2*L14*L25^2*L34^2
                                    -3*L12^2*L13*L14*L25^2*L34^2
                                    -2*L13^4*L25^2*L34^2
                                    +7*L12*L13^3*L25^2*L34^2
                                    -6*L12^2*L13^2*L25^2*L34^2
                                    +L12^3*L13*L25^2*L34^2
                                    +2*L15*L23^3*L24*L25*L34^2
                                    +L13*L23^3*L24*L25*L34^2
                                    +2*L15^2*L23^2*L24*L25*L34^2
                                    -6*L12*L15*L23^2*L24*L25*L34^2
                                    -L13^2*L23^2*L24*L25*L34^2
                                    -3*L12*L13*L23^2*L24*L25*L34^2
                                    -5*L13*L15^2*L23*L24*L25*L34^2
                                    -2*L12*L15^2*L23*L24*L25*L34^2
                                    -4*L13^2*L15*L23*L24*L25*L34^2
                                    +4*L12*L13*L15*L23*L24*L25*L34^2
                                    +6*L12^2*L15*L23*L24*L25*L34^2
                                    -L13^3*L23*L24*L25*L34^2
                                    +2*L12*L13^2*L23*L24*L25*L34^2
                                    -3*L13^2*L15^2*L24*L25*L34^2
                                    +3*L12*L13*L15^2*L24*L25*L34^2
                                    +2*L13^3*L15*L24*L25*L34^2
                                    -2*L12^3*L15*L24*L25*L34^2
                                    +L13^4*L24*L25*L34^2
                                    -3*L12*L13^3*L24*L25*L34^2
                                    +2*L12^3*L13*L24*L25*L34^2
                                    +L15*L23^4*L25*L34^2
                                    +3*L15^2*L23^3*L25*L34^2
                                    +2*L14*L15*L23^3*L25*L34^2
                                    +2*L13*L15*L23^3*L25*L34^2
                                    -4*L12*L15*L23^3*L25*L34^2
                                    +L12*L13*L23^3*L25*L34^2
                                    +4*L14*L15^2*L23^2*L25*L34^2
                                    +2*L13*L15^2*L23^2*L25*L34^2
                                    -7*L12*L15^2*L23^2*L25*L34^2
                                    -4*L13*L14*L15*L23^2*L25*L34^2
                                    -6*L13^2*L15*L23^2*L25*L34^2
                                    -2*L12*L13*L15*L23^2*L25*L34^2
                                    +6*L12^2*L15*L23^2*L25*L34^2
                                    -2*L13^2*L14*L23^2*L25*L34^2
                                    +10*L12*L13*L14*L23^2*L25*L34^2
                                    -L12*L13^2*L23^2*L25*L34^2
                                    +2*L13*L14*L15^2*L23*L25*L34^2
                                    -4*L12*L14*L15^2*L23*L25*L34^2
                                    -4*L13^2*L15^2*L23*L25*L34^2
                                    +3*L12*L13*L15^2*L23*L25*L34^2
                                    +3*L12^2*L15^2*L23*L25*L34^2
                                    +4*L12*L13*L14*L15*L23*L25*L34^2
                                    +2*L13^3*L15*L23*L25*L34^2
                                    -2*L12*L13^2*L15*L23*L25*L34^2
                                    -4*L12^3*L15*L23*L25*L34^2
                                    +4*L13^3*L14*L23*L25*L34^2
                                    -12*L12*L13^2*L14*L23*L25*L34^2
                                    +6*L12^2*L13*L14*L23*L25*L34^2
                                    -L12*L13^3*L23*L25*L34^2
                                    +3*L12^2*L13^2*L23*L25*L34^2
                                    -L13^3*L15^2*L25*L34^2
                                    +3*L12*L13^2*L15^2*L25*L34^2
                                    -3*L12^2*L13*L15^2*L25*L34^2
                                    +L12^3*L15^2*L25*L34^2
                                    +2*L13^3*L14*L15*L25*L34^2
                                    -6*L12*L13^2*L14*L15*L25*L34^2
                                    +6*L12^2*L13*L14*L15*L25*L34^2
                                    -2*L12^3*L14*L15*L25*L34^2
                                    +L13^4*L15*L25*L34^2
                                    -4*L12*L13^3*L15*L25*L34^2
                                    +6*L12^2*L13^2*L15*L25*L34^2
                                    -4*L12^3*L13*L15*L25*L34^2
                                    +L12^4*L15*L25*L34^2-2*L13^4*L14*L25*L34^2
                                    +6*L12*L13^3*L14*L25*L34^2
                                    -6*L12^2*L13^2*L14*L25*L34^2
                                    +2*L12^3*L13*L14*L25*L34^2
                                    +L12*L13^4*L25*L34^2
                                    -3*L12^2*L13^3*L25*L34^2
                                    +3*L12^3*L13^2*L25*L34^2
                                    -L12^4*L13*L25*L34^2-2*L15*L23^4*L24*L34^2
                                    -5*L15^2*L23^3*L24*L34^2
                                    +4*L13*L15*L23^3*L24*L34^2
                                    +6*L12*L15*L23^3*L24*L34^2
                                    -L12*L13*L23^3*L24*L34^2
                                    -2*L15^3*L23^2*L24*L34^2
                                    +9*L13*L15^2*L23^2*L24*L34^2
                                    +8*L12*L15^2*L23^2*L24*L34^2
                                    -2*L13^2*L15*L23^2*L24*L34^2
                                    -12*L12*L13*L15*L23^2*L24*L34^2
                                    -6*L12^2*L15*L23^2*L24*L34^2
                                    +2*L12*L13^2*L23^2*L24*L34^2
                                    +3*L12^2*L13*L23^2*L24*L34^2
                                    +4*L13*L15^3*L23*L24*L34^2
                                    +2*L12*L15^3*L23*L24*L34^2
                                    -4*L13^2*L15^2*L23*L24*L34^2
                                    -11*L12*L13*L15^2*L23*L24*L34^2
                                    -3*L12^2*L15^2*L23*L24*L34^2
                                    +10*L12*L13^2*L15*L23*L24*L34^2
                                    +6*L12^2*L13*L15*L23*L24*L34^2
                                    +2*L12^3*L15*L23*L24*L34^2
                                    -L12*L13^3*L23*L24*L34^2
                                    -3*L12^2*L13^2*L23*L24*L34^2
                                    -2*L12^3*L13*L23*L24*L34^2
                                    -2*L15^2*L23^4*L34^2+L14*L15*L23^4*L34^2
                                    +L12*L15*L23^4*L34^2-2*L15^3*L23^3*L34^2
                                    +L14*L15^2*L23^3*L34^2
                                    +3*L13*L15^2*L23^3*L34^2
                                    +7*L12*L15^2*L23^3*L34^2
                                    -L13*L14*L15*L23^3*L34^2
                                    -3*L12*L14*L15*L23^3*L34^2
                                    -L12*L13*L15*L23^3*L34^2
                                    -3*L12^2*L15*L23^3*L34^2
                                    -L12*L13*L14*L23^3*L34^2
                                    -L14*L15^3*L23^2*L34^2
                                    +L13*L15^3*L23^2*L34^2
                                    +2*L12*L15^3*L23^2*L34^2
                                    -L13*L14*L15^2*L23^2*L34^2
                                    +2*L12*L14*L15^2*L23^2*L34^2
                                    -L12*L13*L15^2*L23^2*L34^2
                                    -6*L12^2*L15^2*L23^2*L34^2
                                    -L13^2*L14*L15*L23^2*L34^2
                                    +2*L12*L13*L14*L15*L23^2*L34^2
                                    -L12*L13^2*L15*L23^2*L34^2
                                    +3*L12^2*L13*L15*L23^2*L34^2
                                    +3*L12^3*L15*L23^2*L34^2
                                    +2*L12*L13^2*L14*L23^2*L34^2
                                    -3*L12^2*L13*L14*L23^2*L34^2
                                    -L12^3*L13*L23^2*L34^2
                                    -L13*L14*L15^3*L23*L34^2
                                    +L12*L14*L15^3*L23*L34^2
                                    +L13^2*L15^3*L23*L34^2
                                    -L12*L13*L15^3*L23*L34^2
                                    +3*L12*L13*L14*L15^2*L23*L34^2
                                    -3*L12^2*L14*L15^2*L23*L34^2
                                    -L13^3*L15^2*L23*L34^2
                                    +L12^3*L15^2*L23*L34^2
                                    +L13^3*L14*L15*L23*L34^2
                                    -3*L12*L13^2*L14*L15*L23*L34^2
                                    +2*L12^3*L14*L15*L23*L34^2
                                    +L12*L13^3*L15*L23*L34^2
                                    -L12^4*L15*L23*L34^2
                                    -L12*L13^3*L14*L23*L34^2
                                    +3*L12^2*L13^2*L14*L23*L34^2
                                    -2*L12^3*L13*L14*L23*L34^2
                                    -L12^3*L13^2*L23*L34^2+L12^4*L13*L23*L34^2
                                    +2*L13*L14*L23*L24*L25^3*L34
                                    +2*L13^2*L23*L24*L25^3*L34
                                    +2*L13^2*L14*L24*L25^3*L34
                                    -2*L12*L13*L14*L24*L25^3*L34
                                    +2*L13^3*L24*L25^3*L34
                                    -2*L12*L13^2*L24*L25^3*L34
                                    -2*L13*L14*L23^2*L25^3*L34
                                    -2*L13^2*L23^2*L25^3*L34
                                    -5*L13*L14^2*L23*L25^3*L34
                                    +2*L12*L13*L14*L23*L25^3*L34
                                    +L13^3*L23*L25^3*L34
                                    +2*L12*L13^2*L23*L25^3*L34
                                    +L13^2*L14^2*L25^3*L34
                                    -L12*L13*L14^2*L25^3*L34
                                    -2*L13^3*L14*L25^3*L34
                                    +2*L12*L13^2*L14*L25^3*L34+L13^4*L25^3*L34
                                    -L12*L13^3*L25^3*L34
                                    -2*L13*L15*L23*L24^2*L25^2*L34
                                    -2*L13^2*L23*L24^2*L25^2*L34
                                    -2*L13^2*L15*L24^2*L25^2*L34
                                    +2*L12*L13*L15*L24^2*L25^2*L34
                                    -2*L13^3*L24^2*L25^2*L34
                                    +2*L12*L13^2*L24^2*L25^2*L34
                                    -2*L13*L15*L23^2*L24*L25^2*L34
                                    -3*L13*L14*L23^2*L24*L25^2*L34
                                    +L13^2*L23^2*L24*L25^2*L34
                                    -8*L13^2*L15*L23*L24*L25^2*L34
                                    +8*L12*L13*L15*L23*L24*L25^2*L34
                                    +4*L13^2*L14*L23*L24*L25^2*L34
                                    -4*L12*L13^2*L23*L24*L25^2*L34
                                    -4*L13^2*L14*L15*L24*L25^2*L34
                                    +4*L12*L13*L14*L15*L24*L25^2*L34
                                    -2*L13^3*L15*L24*L25^2*L34
                                    +8*L12*L13^2*L15*L24*L25^2*L34
                                    -6*L12^2*L13*L15*L24*L25^2*L34
                                    +7*L13^3*L14*L24*L25^2*L34
                                    -10*L12*L13^2*L14*L24*L25^2*L34
                                    +3*L12^2*L13*L14*L24*L25^2*L34
                                    -L13^4*L24*L25^2*L34
                                    -2*L12*L13^3*L24*L25^2*L34
                                    +3*L12^2*L13^2*L24*L25^2*L34
                                    +L14*L15*L23^3*L25^2*L34
                                    +3*L13*L15*L23^3*L25^2*L34
                                    +2*L13*L14*L23^3*L25^2*L34
                                    +2*L13^2*L23^3*L25^2*L34
                                    +3*L14^2*L15*L23^2*L25^2*L34
                                    +5*L13*L14*L15*L23^2*L25^2*L34
                                    -3*L12*L14*L15*L23^2*L25^2*L34
                                    +2*L13^2*L15*L23^2*L25^2*L34
                                    -7*L12*L13*L15*L23^2*L25^2*L34
                                    +8*L13*L14^2*L23^2*L25^2*L34
                                    -7*L13^2*L14*L23^2*L25^2*L34
                                    +3*L12*L13*L14*L23^2*L25^2*L34
                                    -3*L13^3*L23^2*L25^2*L34
                                    -L12*L13^2*L23^2*L25^2*L34
                                    +5*L13*L14^2*L15*L23*L25^2*L34
                                    -3*L12*L14^2*L15*L23*L25^2*L34
                                    +5*L13^2*L14*L15*L23*L25^2*L34
                                    -12*L12*L13*L14*L15*L23*L25^2*L34
                                    +3*L12^2*L14*L15*L23*L25^2*L34
                                    -4*L13^3*L15*L23*L25^2*L34
                                    +3*L12*L13^2*L15*L23*L25^2*L34
                                    +3*L12^2*L13*L15*L23*L25^2*L34
                                    -12*L13^2*L14^2*L23*L25^2*L34
                                    +10*L12*L13*L14^2*L23*L25^2*L34
                                    +6*L13^3*L14*L23*L25^2*L34
                                    +4*L12*L13^2*L14*L23*L25^2*L34
                                    -6*L12^2*L13*L14*L23*L25^2*L34
                                    -2*L12*L13^3*L23*L25^2*L34
                                    +L13^3*L14*L15*L25^2*L34
                                    -3*L12*L13^2*L14*L15*L25^2*L34
                                    +3*L12^2*L13*L14*L15*L25^2*L34
                                    -L12^3*L14*L15*L25^2*L34
                                    -L13^4*L15*L25^2*L34
                                    +3*L12*L13^3*L15*L25^2*L34
                                    -3*L12^2*L13^2*L15*L25^2*L34
                                    +L12^3*L13*L15*L25^2*L34
                                    -L13^4*L14*L25^2*L34
                                    +3*L12*L13^3*L14*L25^2*L34
                                    -3*L12^2*L13^2*L14*L25^2*L34
                                    +L12^3*L13*L14*L25^2*L34+L13^5*L25^2*L34
                                    -3*L12*L13^4*L25^2*L34
                                    +3*L12^2*L13^3*L25^2*L34
                                    -L12^3*L13^2*L25^2*L34
                                    +7*L13*L15*L23^2*L24^2*L25*L34
                                    +L13^2*L23^2*L24^2*L25*L34
                                    +5*L13*L15^2*L23*L24^2*L25*L34
                                    +4*L13^2*L15*L23*L24^2*L25*L34
                                    -10*L12*L13*L15*L23*L24^2*L25*L34
                                    -L13^3*L23*L24^2*L25*L34
                                    +2*L12*L13^2*L23*L24^2*L25*L34
                                    +3*L13^2*L15^2*L24^2*L25*L34
                                    -3*L12*L13*L15^2*L24^2*L25*L34
                                    -3*L13^3*L15*L24^2*L25*L34
                                    +3*L12^2*L13*L15*L24^2*L25*L34
                                    +3*L12*L13^3*L24^2*L25*L34
                                    -3*L12^2*L13^2*L24^2*L25*L34
                                    +L15^2*L23^3*L24*L25*L34
                                    -2*L14*L15*L23^3*L24*L25*L34
                                    -2*L13*L15*L23^3*L24*L25*L34
                                    +2*L13*L14*L23^3*L24*L25*L34
                                    -3*L13^2*L23^3*L24*L25*L34
                                    -4*L14*L15^2*L23^2*L24*L25*L34
                                    +5*L13*L15^2*L23^2*L24*L25*L34
                                    -3*L12*L15^2*L23^2*L24*L25*L34
                                    -6*L13*L14*L15*L23^2*L24*L25*L34
                                    +6*L12*L14*L15*L23^2*L24*L25*L34
                                    -4*L12*L13*L15*L23^2*L24*L25*L34
                                    -4*L13^2*L14*L23^2*L24*L25*L34
                                    +5*L13^3*L23^2*L24*L25*L34
                                    +5*L12*L13^2*L23^2*L24*L25*L34
                                    +4*L12*L14*L15^2*L23*L24*L25*L34
                                    +5*L13^2*L15^2*L23*L24*L25*L34
                                    -12*L12*L13*L15^2*L23*L24*L25*L34
                                    +3*L12^2*L15^2*L23*L24*L25*L34
                                    -6*L13^2*L14*L15*L23*L24*L25*L34
                                    +4*L12*L13*L14*L15*L23*L24*L25*L34
                                    -6*L12^2*L14*L15*L23*L24*L25*L34
                                    +2*L13^3*L15*L23*L24*L25*L34
                                    +6*L12^2*L13*L15*L23*L24*L25*L34
                                    +4*L12*L13^2*L14*L23*L24*L25*L34
                                    -L13^4*L23*L24*L25*L34
                                    -3*L12^2*L13^2*L23*L24*L25*L34
                                    +L13^3*L15^2*L24*L25*L34
                                    -3*L12*L13^2*L15^2*L24*L25*L34
                                    +3*L12^2*L13*L15^2*L24*L25*L34
                                    -L12^3*L15^2*L24*L25*L34
                                    -2*L13^3*L14*L15*L24*L25*L34
                                    +6*L12*L13^2*L14*L15*L24*L25*L34
                                    -6*L12^2*L13*L14*L15*L24*L25*L34
                                    +2*L12^3*L14*L15*L24*L25*L34
                                    +2*L13^4*L14*L24*L25*L34
                                    -6*L12*L13^3*L14*L24*L25*L34
                                    +6*L12^2*L13^2*L14*L24*L25*L34
                                    -2*L12^3*L13*L14*L24*L25*L34
                                    -L13^5*L24*L25*L34+3*L12*L13^4*L24*L25*L34
                                    -3*L12^2*L13^3*L24*L25*L34
                                    +L12^3*L13^2*L24*L25*L34
                                    -L15^2*L23^4*L25*L34
                                    -3*L13*L15*L23^4*L25*L34
                                    -2*L14*L15^2*L23^3*L25*L34
                                    -4*L13*L15^2*L23^3*L25*L34
                                    +3*L12*L15^2*L23^3*L25*L34
                                    -3*L14^2*L15*L23^3*L25*L34
                                    +2*L13*L14*L15*L23^3*L25*L34
                                    +3*L13^2*L15*L23^3*L25*L34
                                    +8*L12*L13*L15*L23^3*L25*L34
                                    -4*L13*L14^2*L23^3*L25*L34
                                    +4*L13^2*L14*L23^3*L25*L34
                                    -6*L12*L13*L14*L23^3*L25*L34
                                    -L12*L13^2*L23^3*L25*L34
                                    -2*L14^2*L15^2*L23^2*L25*L34
                                    -8*L13*L14*L15^2*L23^2*L25*L34
                                    +8*L12*L14*L15^2*L23^2*L25*L34
                                    +2*L13^2*L15^2*L23^2*L25*L34
                                    +3*L12*L13*L15^2*L23^2*L25*L34
                                    -3*L12^2*L15^2*L23^2*L25*L34
                                    +4*L13*L14^2*L15*L23^2*L25*L34
                                    +3*L13^3*L15*L23^2*L25*L34
                                    -4*L12*L13^2*L15*L23^2*L25*L34
                                    -3*L12^2*L13*L15*L23^2*L25*L34
                                    +9*L13^2*L14^2*L23^2*L25*L34
                                    -11*L12*L13*L14^2*L23^2*L25*L34
                                    -8*L13^3*L14*L23^2*L25*L34
                                    +2*L12*L13^2*L14*L23^2*L25*L34
                                    +6*L12^2*L13*L14*L23^2*L25*L34
                                    +2*L12*L13^3*L23^2*L25*L34
                                    -2*L13*L14^2*L15^2*L23*L25*L34
                                    +2*L12*L14^2*L15^2*L23*L25*L34
                                    -2*L13^2*L14*L15^2*L23*L25*L34
                                    +8*L12*L13*L14*L15^2*L23*L25*L34
                                    -6*L12^2*L14*L15^2*L23*L25*L34
                                    +3*L13^3*L15^2*L23*L25*L34
                                    -7*L12*L13^2*L15^2*L23*L25*L34
                                    +3*L12^2*L13*L15^2*L23*L25*L34
                                    +L12^3*L15^2*L23*L25*L34
                                    +7*L13^2*L14^2*L15*L23*L25*L34
                                    -10*L12*L13*L14^2*L15*L23*L25*L34
                                    +3*L12^2*L14^2*L15*L23*L25*L34
                                    -2*L13^3*L14*L15*L23*L25*L34
                                    -4*L12*L13^2*L14*L15*L23*L25*L34
                                    +6*L12^2*L13*L14*L15*L23*L25*L34
                                    -3*L13^4*L15*L23*L25*L34
                                    +8*L12*L13^3*L15*L23*L25*L34
                                    -3*L12^2*L13^2*L15*L23*L25*L34
                                    -2*L12^3*L13*L15*L23*L25*L34
                                    -5*L13^3*L14^2*L23*L25*L34
                                    +8*L12*L13^2*L14^2*L23*L25*L34
                                    -3*L12^2*L13*L14^2*L23*L25*L34
                                    +4*L13^4*L14*L23*L25*L34
                                    -4*L12*L13^3*L14*L23*L25*L34
                                    -L12*L13^4*L23*L25*L34
                                    +L12^3*L13^2*L23*L25*L34
                                    -5*L13*L15*L23^3*L24^2*L34
                                    +L13^2*L23^3*L24^2*L34
                                    +L15^3*L23^2*L24^2*L34
                                    -12*L13*L15^2*L23^2*L24^2*L34
                                    +9*L13^2*L15*L23^2*L24^2*L34
                                    +8*L12*L13*L15*L23^2*L24^2*L34
                                    -2*L13^3*L23^2*L24^2*L34
                                    -4*L12*L13^2*L23^2*L24^2*L34
                                    -5*L13*L15^3*L23*L24^2*L34
                                    -L12*L15^3*L23*L24^2*L34
                                    +8*L13^2*L15^2*L23*L24^2*L34
                                    +10*L12*L13*L15^2*L23*L24^2*L34
                                    -4*L13^3*L15*L23*L24^2*L34
                                    -11*L12*L13^2*L15*L23*L24^2*L34
                                    -3*L12^2*L13*L15*L23*L24^2*L34
                                    +L13^4*L23*L24^2*L34
                                    +2*L12*L13^3*L23*L24^2*L34
                                    +3*L12^2*L13^2*L23*L24^2*L34
                                    -L15^2*L23^4*L24*L34
                                    +2*L14*L15*L23^4*L24*L34
                                    +4*L13*L15*L23^4*L24*L34
                                    -L13*L14*L23^4*L24*L34
                                    -2*L15^3*L23^3*L24*L34
                                    +7*L14*L15^2*L23^3*L24*L34
                                    +6*L13*L15^2*L23^3*L24*L34
                                    +3*L12*L15^2*L23^3*L24*L34
                                    -6*L12*L14*L15*L23^3*L24*L34
                                    -8*L13^2*L15*L23^3*L24*L34
                                    -4*L12*L13*L15*L23^3*L24*L34
                                    +L13^2*L14*L23^3*L24*L34
                                    +2*L12*L13*L14*L23^3*L24*L34
                                    +L12*L13^2*L23^3*L24*L34
                                    +2*L14*L15^3*L23^2*L24*L34
                                    +2*L12*L15^3*L23^2*L24*L34
                                    +4*L13*L14*L15^2*L23^2*L24*L34
                                    -10*L12*L14*L15^2*L23^2*L24*L34
                                    -7*L13^2*L15^2*L23^2*L24*L34
                                    +4*L12*L13*L15^2*L23^2*L24*L34
                                    -3*L12^2*L15^2*L23^2*L24*L34
                                    -4*L13^2*L14*L15*L23^2*L24*L34
                                    +4*L12*L13*L14*L15*L23^2*L24*L34
                                    +6*L12^2*L14*L15*L23^2*L24*L34
                                    +4*L13^3*L15*L23^2*L24*L34
                                    +2*L12*L13^2*L15*L23^2*L24*L34
                                    +L13^3*L14*L23^2*L24*L34
                                    -3*L12^2*L13*L14*L23^2*L24*L34
                                    -2*L12*L13^3*L23^2*L24*L34
                                    +2*L13*L14*L15^3*L23*L24*L34
                                    -2*L12*L14*L15^3*L23*L24*L34
                                    -2*L13^2*L15^3*L23*L24*L34
                                    +2*L12*L13*L15^3*L23*L24*L34
                                    -3*L13^2*L14*L15^2*L23*L24*L34
                                    +3*L12^2*L14*L15^2*L23*L24*L34
                                    +2*L13^3*L15^2*L23*L24*L34
                                    +3*L12*L13^2*L15^2*L23*L24*L34
                                    -6*L12^2*L13*L15^2*L23*L24*L34
                                    +L12^3*L15^2*L23*L24*L34
                                    +2*L13^3*L14*L15*L23*L24*L34
                                    -2*L12^3*L14*L15*L23*L24*L34
                                    -6*L12*L13^3*L15*L23*L24*L34
                                    +6*L12^2*L13^2*L15*L23*L24*L34
                                    -L13^4*L14*L23*L24*L34
                                    +2*L12*L13^3*L14*L23*L24*L34
                                    -3*L12^2*L13^2*L14*L23*L24*L34
                                    +2*L12^3*L13*L14*L23*L24*L34
                                    +L12*L13^4*L23*L24*L34
                                    -L12^3*L13^2*L23*L24*L34+L15^2*L23^5*L34
                                    -L14*L15*L23^5*L34+L15^3*L23^4*L34
                                    -L14*L15^2*L23^4*L34-3*L12*L15^2*L23^4*L34
                                    -L13*L14*L15*L23^4*L34
                                    +3*L12*L14*L15*L23^4*L34
                                    -L12*L13*L15*L23^4*L34+L13*L14^2*L23^4*L34
                                    +L12*L13*L14*L23^4*L34
                                    +2*L14*L15^3*L23^3*L34+L13*L15^3*L23^3*L34
                                    -L12*L15^3*L23^3*L34
                                    -2*L14^2*L15^2*L23^3*L34
                                    -2*L12*L14*L15^2*L23^3*L34
                                    -3*L13^2*L15^2*L23^3*L34
                                    -2*L12*L13*L15^2*L23^3*L34
                                    +3*L12^2*L15^2*L23^3*L34
                                    -L13*L14^2*L15*L23^3*L34
                                    +3*L12*L14^2*L15*L23^3*L34
                                    +5*L13^2*L14*L15*L23^3*L34
                                    -3*L12^2*L14*L15*L23^3*L34
                                    +2*L12*L13^2*L15*L23^3*L34
                                    -2*L13^2*L14^2*L23^3*L34
                                    +2*L12*L13*L14^2*L23^3*L34
                                    -2*L12*L13^2*L14*L23^3*L34
                                    +2*L13*L14*L15^3*L23^2*L34
                                    -2*L12*L14*L15^3*L23^2*L34
                                    -2*L13^2*L15^3*L23^2*L34
                                    +2*L12*L13*L15^3*L23^2*L34
                                    -2*L13*L14^2*L15^2*L23^2*L34
                                    +2*L12*L14^2*L15^2*L23^2*L34
                                    +L13^2*L14*L15^2*L23^2*L34
                                    -4*L12*L13*L14*L15^2*L23^2*L34
                                    +3*L12^2*L14*L15^2*L23^2*L34
                                    +2*L13^3*L15^2*L23^2*L34
                                    -L12*L13^2*L15^2*L23^2*L34
                                    -L12^3*L15^2*L23^2*L34
                                    +L13^2*L14^2*L15*L23^2*L34
                                    +2*L12*L13*L14^2*L15*L23^2*L34
                                    -3*L12^2*L14^2*L15*L23^2*L34
                                    -3*L13^3*L14*L15*L23^2*L34
                                    +5*L12*L13^2*L14*L15*L23^2*L34
                                    -3*L12^2*L13*L14*L15*L23^2*L34
                                    +L12^3*L14*L15*L23^2*L34
                                    -L12*L13^3*L15*L23^2*L34
                                    +L12^3*L13*L15*L23^2*L34
                                    +L13^3*L14^2*L23^2*L34
                                    -4*L12*L13^2*L14^2*L23^2*L34
                                    +3*L12^2*L13*L14^2*L23^2*L34
                                    +L12*L13^3*L14*L23^2*L34
                                    -L12^3*L13*L14*L23^2*L34
                                    -L13*L14^2*L23*L24*L25^3
                                    +2*L13^2*L14*L23*L24*L25^3
                                    -L13^3*L23*L24*L25^3-L13^2*L14^2*L24*L25^3
                                    +L12*L13*L14^2*L24*L25^3
                                    +2*L13^3*L14*L24*L25^3
                                    -2*L12*L13^2*L14*L24*L25^3-L13^4*L24*L25^3
                                    +L12*L13^3*L24*L25^3+L13*L14^2*L23^2*L25^3
                                    -2*L13^2*L14*L23^2*L25^3+L13^3*L23^2*L25^3
                                    +2*L13*L14^3*L23*L25^3
                                    -5*L13^2*L14^2*L23*L25^3
                                    -L12*L13*L14^2*L23*L25^3
                                    +4*L13^3*L14*L23*L25^3
                                    +2*L12*L13^2*L14*L23*L25^3-L13^4*L23*L25^3
                                    -L12*L13^3*L23*L25^3
                                    +2*L13*L14*L15*L23*L24^2*L25^2
                                    -2*L13^2*L15*L23*L24^2*L25^2
                                    -2*L13^2*L14*L23*L24^2*L25^2
                                    +2*L13^3*L23*L24^2*L25^2
                                    +2*L13^2*L14*L15*L24^2*L25^2
                                    -2*L12*L13*L14*L15*L24^2*L25^2
                                    -2*L13^3*L15*L24^2*L25^2
                                    +2*L12*L13^2*L15*L24^2*L25^2
                                    -2*L13^3*L14*L24^2*L25^2
                                    +2*L12*L13^2*L14*L24^2*L25^2
                                    +2*L13^4*L24^2*L25^2
                                    -2*L12*L13^3*L24^2*L25^2
                                    -4*L13*L14*L15*L23^2*L24*L25^2
                                    +4*L13^2*L15*L23^2*L24*L25^2
                                    +3*L13*L14^2*L23^2*L24*L25^2
                                    -2*L13^2*L14*L23^2*L24*L25^2
                                    -L13^3*L23^2*L24*L25^2
                                    -2*L13*L14^2*L15*L23*L24*L25^2
                                    +4*L12*L13*L14*L15*L23*L24*L25^2
                                    +2*L13^3*L15*L23*L24*L25^2
                                    -4*L12*L13^2*L15*L23*L24*L25^2
                                    +5*L13^2*L14^2*L23*L24*L25^2
                                    -3*L12*L13*L14^2*L23*L24*L25^2
                                    -6*L13^3*L14*L23*L24*L25^2
                                    +2*L12*L13^2*L14*L23*L24*L25^2
                                    +L13^4*L23*L24*L25^2
                                    +L12*L13^3*L23*L24*L25^2
                                    +2*L13*L14*L15*L23^3*L25^2
                                    -2*L13^2*L15*L23^3*L25^2
                                    -3*L13*L14^2*L23^3*L25^2
                                    +4*L13^2*L14*L23^3*L25^2-L13^3*L23^3*L25^2
                                    -L14^3*L15*L23^2*L25^2
                                    +5*L13*L14^2*L15*L23^2*L25^2
                                    -5*L13^2*L14*L15*L23^2*L25^2
                                    -2*L12*L13*L14*L15*L23^2*L25^2
                                    +L13^3*L15*L23^2*L25^2
                                    +2*L12*L13^2*L15*L23^2*L25^2
                                    -5*L13*L14^3*L23^2*L25^2
                                    +10*L13^2*L14^2*L23^2*L25^2
                                    +3*L12*L13*L14^2*L23^2*L25^2
                                    -7*L13^3*L14*L23^2*L25^2
                                    -4*L12*L13^2*L14*L23^2*L25^2
                                    +2*L13^4*L23^2*L25^2+L12*L13^3*L23^2*L25^2
                                    -L13*L14^3*L15*L23*L25^2
                                    +L12*L14^3*L15*L23*L25^2
                                    +3*L13^2*L14^2*L15*L23*L25^2
                                    -3*L12*L13*L14^2*L15*L23*L25^2
                                    -3*L13^3*L14*L15*L23*L25^2
                                    +3*L12*L13^2*L14*L15*L23*L25^2
                                    +L13^4*L15*L23*L25^2
                                    -L12*L13^3*L15*L23*L25^2
                                    +L13^2*L14^3*L23*L25^2
                                    -L12*L13*L14^3*L23*L25^2
                                    -3*L13^3*L14^2*L23*L25^2
                                    +3*L12*L13^2*L14^2*L23*L25^2
                                    +3*L13^4*L14*L23*L25^2
                                    -3*L12*L13^3*L14*L23*L25^2-L13^5*L23*L25^2
                                    +L12*L13^4*L23*L25^2
                                    -L13*L15^2*L23*L24^3*L25
                                    +2*L13^2*L15*L23*L24^3*L25
                                    -L13^3*L23*L24^3*L25-L13^2*L15^2*L24^3*L25
                                    +L12*L13*L15^2*L24^3*L25
                                    +2*L13^3*L15*L24^3*L25
                                    -2*L12*L13^2*L15*L24^3*L25-L13^4*L24^3*L25
                                    +L12*L13^3*L24^3*L25
                                    +3*L13*L15^2*L23^2*L24^2*L25
                                    -4*L13*L14*L15*L23^2*L24^2*L25
                                    -2*L13^2*L15*L23^2*L24^2*L25
                                    +4*L13^2*L14*L23^2*L24^2*L25
                                    -L13^3*L23^2*L24^2*L25
                                    -2*L13*L14*L15^2*L23*L24^2*L25
                                    +5*L13^2*L15^2*L23*L24^2*L25
                                    -3*L12*L13*L15^2*L23*L24^2*L25
                                    +4*L12*L13*L14*L15*L23*L24^2*L25
                                    -6*L13^3*L15*L23*L24^2*L25
                                    +2*L12*L13^2*L15*L23*L24^2*L25
                                    +2*L13^3*L14*L23*L24^2*L25
                                    -4*L12*L13^2*L14*L23*L24^2*L25
                                    +L13^4*L23*L24^2*L25
                                    +L12*L13^3*L23*L24^2*L25
                                    -3*L13*L15^2*L23^3*L24*L25
                                    +8*L13*L14*L15*L23^3*L24*L25
                                    -2*L13^2*L15*L23^3*L24*L25
                                    -3*L13*L14^2*L23^3*L24*L25
                                    -2*L13^2*L14*L23^3*L24*L25
                                    +2*L13^3*L23^3*L24*L25
                                    +2*L14^2*L15^2*L23^2*L24*L25
                                    -5*L13^2*L15^2*L23^2*L24*L25
                                    +3*L12*L13*L15^2*L23^2*L24*L25
                                    -8*L12*L13*L14*L15*L23^2*L24*L25
                                    +6*L13^3*L15*L23^2*L24*L25
                                    +2*L12*L13^2*L15*L23^2*L24*L25
                                    -5*L13^2*L14^2*L23^2*L24*L25
                                    +3*L12*L13*L14^2*L23^2*L24*L25
                                    +6*L13^3*L14*L23^2*L24*L25
                                    +2*L12*L13^2*L14*L23^2*L24*L25
                                    -4*L13^4*L23^2*L24*L25
                                    -2*L12*L13^3*L23^2*L24*L25
                                    +2*L13*L14^2*L15^2*L23*L24*L25
                                    -2*L12*L14^2*L15^2*L23*L24*L25
                                    -4*L13^2*L14*L15^2*L23*L24*L25
                                    +4*L12*L13*L14*L15^2*L23*L24*L25
                                    +2*L13^3*L15^2*L23*L24*L25
                                    -2*L12*L13^2*L15^2*L23*L24*L25
                                    -4*L13^2*L14^2*L15*L23*L24*L25
                                    +4*L12*L13*L14^2*L15*L23*L24*L25
                                    +8*L13^3*L14*L15*L23*L24*L25
                                    -8*L12*L13^2*L14*L15*L23*L24*L25
                                    -4*L13^4*L15*L23*L24*L25
                                    +4*L12*L13^3*L15*L23*L24*L25
                                    +2*L13^3*L14^2*L23*L24*L25
                                    -2*L12*L13^2*L14^2*L23*L24*L25
                                    -4*L13^4*L14*L23*L24*L25
                                    +4*L12*L13^3*L14*L23*L24*L25
                                    +2*L13^5*L23*L24*L25
                                    -2*L12*L13^4*L23*L24*L25
                                    +L13*L15^2*L23^4*L25
                                    -4*L13*L14*L15*L23^4*L25
                                    +2*L13^2*L15*L23^4*L25
                                    +3*L13*L14^2*L23^4*L25
                                    -2*L13^2*L14*L23^4*L25
                                    -2*L14^2*L15^2*L23^3*L25
                                    +2*L13*L14*L15^2*L23^3*L25
                                    +L13^2*L15^2*L23^3*L25
                                    -L12*L13*L15^2*L23^3*L25
                                    +2*L14^3*L15*L23^3*L25
                                    -6*L13*L14^2*L15*L23^3*L25
                                    +6*L13^2*L14*L15*L23^3*L25
                                    +4*L12*L13*L14*L15*L23^3*L25
                                    -4*L13^3*L15*L23^3*L25
                                    -2*L12*L13^2*L15*L23^3*L25
                                    +4*L13*L14^3*L23^3*L25
                                    -7*L13^2*L14^2*L23^3*L25
                                    -3*L12*L13*L14^2*L23^3*L25
                                    +4*L13^3*L14*L23^3*L25
                                    +2*L12*L13^2*L14*L23^3*L25
                                    -2*L13*L14^2*L15^2*L23^2*L25
                                    +2*L12*L14^2*L15^2*L23^2*L25
                                    +4*L13^2*L14*L15^2*L23^2*L25
                                    -4*L12*L13*L14*L15^2*L23^2*L25
                                    -2*L13^3*L15^2*L23^2*L25
                                    +2*L12*L13^2*L15^2*L23^2*L25
                                    +2*L13*L14^3*L15*L23^2*L25
                                    -2*L12*L14^3*L15*L23^2*L25
                                    -2*L13^2*L14^2*L15*L23^2*L25
                                    +2*L12*L13*L14^2*L15*L23^2*L25
                                    -2*L13^3*L14*L15*L23^2*L25
                                    +2*L12*L13^2*L14*L15*L23^2*L25
                                    +2*L13^4*L15*L23^2*L25
                                    -2*L12*L13^3*L15*L23^2*L25
                                    -2*L13^2*L14^3*L23^2*L25
                                    +2*L12*L13*L14^3*L23^2*L25
                                    +4*L13^3*L14^2*L23^2*L25
                                    -4*L12*L13^2*L14^2*L23^2*L25
                                    -2*L13^4*L14*L23^2*L25
                                    +2*L12*L13^3*L14*L23^2*L25
                                    +L13*L15^2*L23^2*L24^3
                                    -2*L13^2*L15*L23^2*L24^3+L13^3*L23^2*L24^3
                                    +2*L13*L15^3*L23*L24^3
                                    -5*L13^2*L15^2*L23*L24^3
                                    -L12*L13*L15^2*L23*L24^3
                                    +4*L13^3*L15*L23*L24^3
                                    +2*L12*L13^2*L15*L23*L24^3-L13^4*L23*L24^3
                                    -L12*L13^3*L23*L24^3
                                    -3*L13*L15^2*L23^3*L24^2
                                    +2*L13*L14*L15*L23^3*L24^2
                                    +4*L13^2*L15*L23^3*L24^2
                                    -2*L13^2*L14*L23^3*L24^2-L13^3*L23^3*L24^2
                                    -L14*L15^3*L23^2*L24^2
                                    -5*L13*L15^3*L23^2*L24^2
                                    +5*L13*L14*L15^2*L23^2*L24^2
                                    +10*L13^2*L15^2*L23^2*L24^2
                                    +3*L12*L13*L15^2*L23^2*L24^2
                                    -5*L13^2*L14*L15*L23^2*L24^2
                                    -2*L12*L13*L14*L15*L23^2*L24^2
                                    -7*L13^3*L15*L23^2*L24^2
                                    -4*L12*L13^2*L15*L23^2*L24^2
                                    +L13^3*L14*L23^2*L24^2
                                    +2*L12*L13^2*L14*L23^2*L24^2
                                    +2*L13^4*L23^2*L24^2+L12*L13^3*L23^2*L24^2
                                    -L13*L14*L15^3*L23*L24^2
                                    +L12*L14*L15^3*L23*L24^2
                                    +L13^2*L15^3*L23*L24^2
                                    -L12*L13*L15^3*L23*L24^2
                                    +3*L13^2*L14*L15^2*L23*L24^2
                                    -3*L12*L13*L14*L15^2*L23*L24^2
                                    -3*L13^3*L15^2*L23*L24^2
                                    +3*L12*L13^2*L15^2*L23*L24^2
                                    -3*L13^3*L14*L15*L23*L24^2
                                    +3*L12*L13^2*L14*L15*L23*L24^2
                                    +3*L13^4*L15*L23*L24^2
                                    -3*L12*L13^3*L15*L23*L24^2
                                    +L13^4*L14*L23*L24^2
                                    -L12*L13^3*L14*L23*L24^2-L13^5*L23*L24^2
                                    +L12*L13^4*L23*L24^2+3*L13*L15^2*L23^4*L24
                                    -4*L13*L14*L15*L23^4*L24
                                    -2*L13^2*L15*L23^4*L24+L13*L14^2*L23^4*L24
                                    +2*L13^2*L14*L23^4*L24
                                    +2*L14*L15^3*L23^3*L24
                                    +4*L13*L15^3*L23^3*L24
                                    -2*L14^2*L15^2*L23^3*L24
                                    -6*L13*L14*L15^2*L23^3*L24
                                    -7*L13^2*L15^2*L23^3*L24
                                    -3*L12*L13*L15^2*L23^3*L24
                                    +2*L13*L14^2*L15*L23^3*L24
                                    +6*L13^2*L14*L15*L23^3*L24
                                    +4*L12*L13*L14*L15*L23^3*L24
                                    +4*L13^3*L15*L23^3*L24
                                    +2*L12*L13^2*L15*L23^3*L24
                                    +L13^2*L14^2*L23^3*L24
                                    -L12*L13*L14^2*L23^3*L24
                                    -4*L13^3*L14*L23^3*L24
                                    -2*L12*L13^2*L14*L23^3*L24
                                    +2*L13*L14*L15^3*L23^2*L24
                                    -2*L12*L14*L15^3*L23^2*L24
                                    -2*L13^2*L15^3*L23^2*L24
                                    +2*L12*L13*L15^3*L23^2*L24
                                    -2*L13*L14^2*L15^2*L23^2*L24
                                    +2*L12*L14^2*L15^2*L23^2*L24
                                    -2*L13^2*L14*L15^2*L23^2*L24
                                    +2*L12*L13*L14*L15^2*L23^2*L24
                                    +4*L13^3*L15^2*L23^2*L24
                                    -4*L12*L13^2*L15^2*L23^2*L24
                                    +4*L13^2*L14^2*L15*L23^2*L24
                                    -4*L12*L13*L14^2*L15*L23^2*L24
                                    -2*L13^3*L14*L15*L23^2*L24
                                    +2*L12*L13^2*L14*L15*L23^2*L24
                                    -2*L13^4*L15*L23^2*L24
                                    +2*L12*L13^3*L15*L23^2*L24
                                    -2*L13^3*L14^2*L23^2*L24
                                    +2*L12*L13^2*L14^2*L23^2*L24
                                    +2*L13^4*L14*L23^2*L24
                                    -2*L12*L13^3*L14*L23^2*L24-L13*L15^2*L23^5
                                    +2*L13*L14*L15*L23^5-L13*L14^2*L23^5
                                    -L14*L15^3*L23^4-L13*L15^3*L23^4
                                    +2*L14^2*L15^2*L23^4+L13*L14*L15^2*L23^4
                                    +2*L13^2*L15^2*L23^4+L12*L13*L15^2*L23^4
                                    -L14^3*L15*L23^4+L13*L14^2*L15*L23^4
                                    -4*L13^2*L14*L15*L23^4
                                    -2*L12*L13*L14*L15*L23^4-L13*L14^3*L23^4
                                    +2*L13^2*L14^2*L23^4+L12*L13*L14^2*L23^4
                                    -L13*L14*L15^3*L23^3+L12*L14*L15^3*L23^3
                                    +L13^2*L15^3*L23^3-L12*L13*L15^3*L23^3
                                    +2*L13*L14^2*L15^2*L23^3
                                    -2*L12*L14^2*L15^2*L23^3
                                    -L13^2*L14*L15^2*L23^3
                                    +L12*L13*L14*L15^2*L23^3-L13^3*L15^2*L23^3
                                    +L12*L13^2*L15^2*L23^3-L13*L14^3*L15*L23^3
                                    +L12*L14^3*L15*L23^3-L13^2*L14^2*L15*L23^3
                                    +L12*L13*L14^2*L15*L23^3
                                    +2*L13^3*L14*L15*L23^3
                                    -2*L12*L13^2*L14*L15*L23^3
                                    +L13^2*L14^3*L23^3-L12*L13*L14^3*L23^3
                                    -L13^3*L14^2*L23^3+L12*L13^2*L14^2*L23^3)
       /(4*(L12*L34^2+L23*L24*L34-L13*L24*L34-L12*L24*L34-L14*L23*L34
                     -L12*L23*L34+L13*L14*L34-L12*L14*L34-L12*L13*L34
                     +L12^2*L34+L13*L24^2-L14*L23*L24-L13*L23*L24-L13*L14*L24
                     +L12*L14*L24+L13^2*L24-L12*L13*L24+L14*L23^2+L14^2*L23
                     -L13*L14*L23-L12*L14*L23+L12*L13*L23)
           ^2
          *(L12*L35^2+L23*L25*L35-L13*L25*L35-L12*L25*L35-L15*L23*L35
                     -L12*L23*L35+L13*L15*L35-L12*L15*L35-L12*L13*L35
                     +L12^2*L35+L13*L25^2-L15*L23*L25-L13*L23*L25-L13*L15*L25
                     +L12*L15*L25+L13^2*L25-L12*L13*L25+L15*L23^2+L15^2*L23
                     -L13*L15*L23-L12*L15*L23+L12*L13*L23)
           ^2)
  ;

if n=0 then r:=r0;
  elif n=1 then r:=r1;
  elif n=2 then r:=r2;
fi;

# теперь переходим к искомой величине d, умножив которую на ориентацию, получим матричный элемент

A := [ [w2-w1, x2-x1, y2-y1, z2-z1], [w3-w1, x3-x1, y3-y1, z3-z1], [w4-w1, x4-x1, y4-y1, z4-z1] ];
v1234 := [ Determinant(A{[1,2,3]}{[1,2,3]}), -Determinant(A{[1,2,3]}{[1,2,4]}),
  Determinant(A{[1,2,3]}{[1,3,4]}), -Determinant(A{[1,2,3]}{[2,3,4]}) ]; 

B := [ [w2-w1, x2-x1, y2-y1, z2-z1], [w3-w1, x3-x1, y3-y1, z3-z1], [w5-w1, x5-x1, y5-y1, z5-z1] ]; 
v1235 := [ Determinant(B{[1,2,3]}{[1,2,3]}), -Determinant(B{[1,2,3]}{[1,2,4]}),
  Determinant(B{[1,2,3]}{[1,3,4]}), -Determinant(B{[1,2,3]}{[2,3,4]}) ]; 

CM123 := Determinant([ [0,1,1,1], [1,0,L12,L13], [1,L12,0,L23], [1,L13,L23,0] ]);

W12345 := Determinant([ [w2-w1, x2-x1, y2-y1, z2-z1], [w3-w1, x3-x1, y3-y1, z3-z1],
  [w4-w1, x4-x1, y4-y1, z4-z1], [w5-w1, x5-x1, y5-y1, z5-z1] ]);

d := r * v1234^2 * v1235^2 / ( CM123 * W12345 * (v1234*v1235) )
  * SignPerm( PermListList( four_simplex, l ) ) ;


return d;

end );


###################################################################################################################


InstallGlobalFunction( Matrices4dEuclid4,
# return structure with matrices f2_full,f3_FULL,f4_full,f5
# constructed using Euclidean 4-dimensional geometry
# and supplementary information
# for given 4-manifold with one-component boundary

function(p, coord)

# p is the initial polytope (not necessarily triangulated)
# coord is the list [w,x,y,z], where w,x,y,z are long enough lists of rational numbers
      # and such that no five vertices lie in a 3-plane

  local s         # triangulated p
    ,F            # structure to be returned
    ,inner_faces  # list of inner faces (inner_faces[i] - set of 
                  # inner faces of dimension (i-1)
    ,ind_f1       # index of all edges, inner edges first
    ,ind_f1_li    # position of last inner edge in index of all edges
    ,ind_f2       # index of all 2-faces, inner faces first
    ,ind_f2_li    # position of last inner 2-face in index of all 2-faces
    ,faces_ind    # backward index of the faces
                  # faces_ind[d+1][i] - index of face of dimension d 
                  # in the list of "our" faces of corresponding dimension
    ,f2_full      # 
    ,f3_FULL      #
    ,f4_full      #
    ,f5           # matrices that we are looking for
    ,i,j,k,l      # 
    ,k1,k2        # counters
    ,nl           # 
    ,l_v,j_v      #
    ,v, Lij       #
    ,j_A,j_B,j_C  #
    ,l1_v,l2_v    # and local stuff
    ,w,x,y,z      # arrays of coordinates
    ,s4_orient    # orientations of 4-simplices
    ,edge         #
    ,two_face     #
    ,four_simplex # the stuff needed for f3_FULL
    ;

  # First triangulate, if needed
  s := PolTriangulate(p);

  # Get index of inner faces
  inner_faces := PolInnerFaces(s);
  # remember about PolInnerFaces: <returned value>[i] --- set of inner faces of dimension (i-1)

  # Enumeration of inner vertices is done automatically in inner_faces.
  # The same is true for inner edges and inner 2-faces

  # we need vertices (inner only), edges and 2-faces

  # Enumerate all edges, inner edges first
  ind_f1 := StructuralCopy(inner_faces[2]);
  ind_f1_li := Length(inner_faces[2]);
  Append(ind_f1, Difference([1..Length(s.faces[1])], ind_f1));

  # Enumerate all 2-faces (inner 2-faces first)
  ind_f2 := StructuralCopy(inner_faces[3]);
  ind_f2_li := Length(inner_faces[3]);
  Append(ind_f2, Difference([1..Length(s.faces[2])], ind_f2));

  # build backward index
  # as we remember, the simplex dimension, here again, is the first index minus one, 
  # e.g., faces_ind[1] contains numbers of inner vertices
  faces_ind := List([1..5], i->[]);
  faces_ind[1] := List([1..Length(s.vertices)], i->0);
  for i in [2..5] do
    faces_ind[i] := List([1..Length(s.faces[i-1])], j->0);
  od;

  for j in [1..Length(inner_faces[1])] do
    faces_ind[1][inner_faces[1][j]] := j;
  od; ##### внутренней вершине сопоставляется её номер по матрицам,
      ##### т.е. среди только внутренних; граничным - 0

  for i in [1..Length(ind_f1)] do
    faces_ind[2][ind_f1[i]] := i;
  od;
  for i in [1..Length(ind_f2)] do
    faces_ind[3][ind_f2[i]] := i;
  od;

  # init matrices
  f2_full := NullMat(Length(ind_f1), 4*Length(inner_faces[1]));
  f3_FULL := NullMat(Length(ind_f2), Length(ind_f1));
  f4_full := NullMat(4*Length(inner_faces[2]),  Length(ind_f2));
  f5 := NullMat(6*Length(inner_faces[1]), 4*Length(inner_faces[2]));

  # Now: fill out elements of our matrices that are non-zero

  # First of all, get the coordinates:
  w := coord[1]; x := coord[2]; y := coord[3]; z := coord[4]; 

  # f2_full: 
  nl := DimensionsMat(f2_full)[1];
  # loop over rows
  for l in [1..nl] do
    # get sorted list of two vertices bounding l-th edge
    l_v := PolBnd(s,1,ind_f1[l])[1]; ##### помним: ind_f1[l] это исходный номер ребра, которому соответствует l-я строка
    Sort(l_v);
    i := faces_ind[1][l_v[1]]; # i это номер по матрицам вершины с исходным номером l_v[1]
    j := faces_ind[1][l_v[2]]; # j аналогично
    # fill out some elements of the row f2[l]
    if i<>0 then
      f2_full[l][4*i-3] := w[l_v[1]]-w[l_v[2]];
      f2_full[l][4*i-2] := x[l_v[1]]-x[l_v[2]];
      f2_full[l][4*i-1] := y[l_v[1]]-y[l_v[2]];
      f2_full[l][4*i] := z[l_v[1]]-z[l_v[2]];
    fi;
    if j<>0 then
      f2_full[l][4*j-3] := w[l_v[2]]-w[l_v[1]];
      f2_full[l][4*j-2] := x[l_v[2]]-x[l_v[1]];
      f2_full[l][4*j-1] := y[l_v[2]]-y[l_v[1]];
      f2_full[l][4*j] := z[l_v[2]]-z[l_v[1]];
    fi;
  od;

  # f3_FULL:
  # make consistent orientation on 4-simplices
  s4_orient := OrientTriangulated(s);
  if s4_orient=[] then
    Print("Non-orientable manifold, aborting!\n");
  fi;
  # loop over all 4-simplices
  for l in [1..Length(s.faces[4])] do
    four_simplex := PolBnd(s,4,l)[1]; Sort(four_simplex);
    # get sorted lists of edges and 2-faces in the boundary of l-th 4-simplex
    l1_v := PolBnd(s,4,l)[2]; Sort(l1_v);
    l2_v := PolBnd(s,4,l)[3]; Sort(l2_v);
    for i in l1_v do
      edge := PolBnd(s,1,i)[1]; Sort(edge);
      for j in l2_v do
      two_face := PolBnd(s,2,j)[1]; Sort(two_face);
      ##### делаем цикл по i in l1_v и j in l2_v: делаем edge и two_face как (sorted) списки вершин 
      ##### и *прибавляем* к матричному элементу соответствующую производную
      ##### а точнее, нужно взять от i и j faces_ind[2] и faces_ind[3] соответственно
      f3_FULL[ faces_ind[3][j] ][ faces_ind[2][i] ] := f3_FULL[ faces_ind[3][j] ][ faces_ind[2][i] ] 
       + s4_orient[l]*f3E4( four_simplex, two_face, edge, coord );
      od;
    od;
  od;

  # f4_full
  # loop over all 2-faces
  for l in [1..Length(s.faces[2])] do
    # get sorted list of vertices # сортируем по традиции
    l_v := PolBnd(s,2,l)[1]; Sort(l_v);
    # loop over edges bounding this 2-face
    for j in s.faces[2][l] do
      # get sorted list of vertices
      v := s.faces[1][j]; 
      # get the number of 2-face's vertex outside the edge
      j_A := Difference( l_v, v )[1]; # потому что Difference состоит из одного элемента
      # и два номера вершин, принадлежащих ребру
      j_B := v[1]; j_C := v[2];
      k1 := faces_ind[2][j]; # номер данного ребра по матрицам
      k2 := faces_ind[3][l]; # номер данной 2-грани по матрицам
      # Fill out respective f4_full elements
        if k1 <= ind_f1_li then # чтобы не вылезти за пределы матрицы
          f4_full[4*k1-3][k2] := -(w[j_B]*z[j_C]^2-w[j_A]*z[j_C]^2-z[j_B]*w[j_C]*z[j_C]
                       +z[j_A]*w[j_C]*z[j_C]-w[j_B]*z[j_B]*z[j_C]
                       +2*w[j_A]*z[j_B]*z[j_C]-z[j_A]*w[j_B]*z[j_C]
                       +w[j_B]*y[j_C]^2-w[j_A]*y[j_C]^2-y[j_B]*w[j_C]*y[j_C]
                       +y[j_A]*w[j_C]*y[j_C]-w[j_B]*y[j_B]*y[j_C]
                       +2*w[j_A]*y[j_B]*y[j_C]-y[j_A]*w[j_B]*y[j_C]
                       +w[j_B]*x[j_C]^2-w[j_A]*x[j_C]^2-x[j_B]*w[j_C]*x[j_C]
                       +x[j_A]*w[j_C]*x[j_C]-w[j_B]*x[j_B]*x[j_C]
                       +2*w[j_A]*x[j_B]*x[j_C]-x[j_A]*w[j_B]*x[j_C]
                       +z[j_B]^2*w[j_C]-z[j_A]*z[j_B]*w[j_C]+y[j_B]^2*w[j_C]
                       -y[j_A]*y[j_B]*w[j_C]+x[j_B]^2*w[j_C]
                       -x[j_A]*x[j_B]*w[j_C]-w[j_A]*z[j_B]^2
                       +z[j_A]*w[j_B]*z[j_B]-w[j_A]*y[j_B]^2
                       +y[j_A]*w[j_B]*y[j_B]-w[j_A]*x[j_B]^2
                       +x[j_A]*w[j_B]*x[j_B]) /4; # производная по координате дубль-вэ точки j_A
          f4_full[4*k1-2][k2] := -(x[j_B]*z[j_C]^2-x[j_A]*z[j_C]^2-z[j_B]*x[j_C]*z[j_C]
                        +z[j_A]*x[j_C]*z[j_C]-x[j_B]*z[j_B]*z[j_C]
                        +2*x[j_A]*z[j_B]*z[j_C]-z[j_A]*x[j_B]*z[j_C]
                        +x[j_B]*y[j_C]^2-x[j_A]*y[j_C]^2-y[j_B]*x[j_C]*y[j_C]
                        +y[j_A]*x[j_C]*y[j_C]-x[j_B]*y[j_B]*y[j_C]
                        +2*x[j_A]*y[j_B]*y[j_C]-y[j_A]*x[j_B]*y[j_C]
                        -w[j_B]*w[j_C]*x[j_C]+w[j_A]*w[j_C]*x[j_C]
                        +z[j_B]^2*x[j_C]-z[j_A]*z[j_B]*x[j_C]+y[j_B]^2*x[j_C]
                        -y[j_A]*y[j_B]*x[j_C]+w[j_B]^2*x[j_C]
                        -w[j_A]*w[j_B]*x[j_C]+x[j_B]*w[j_C]^2-x[j_A]*w[j_C]^2
                        -w[j_B]*x[j_B]*w[j_C]-w[j_A]*x[j_B]*w[j_C]
                        +2*x[j_A]*w[j_B]*w[j_C]-x[j_A]*z[j_B]^2
                        +z[j_A]*x[j_B]*z[j_B]-x[j_A]*y[j_B]^2
                        +y[j_A]*x[j_B]*y[j_B]+w[j_A]*w[j_B]*x[j_B]
                        -x[j_A]*w[j_B]^2) /4; # производная по координате икс точки j_A
          f4_full[4*k1-1][k2] := -(y[j_B]*z[j_C]^2-y[j_A]*z[j_C]^2-z[j_B]*y[j_C]*z[j_C]
                        +z[j_A]*y[j_C]*z[j_C]-y[j_B]*z[j_B]*z[j_C]
                        +2*y[j_A]*z[j_B]*z[j_C]-z[j_A]*y[j_B]*z[j_C]
                        -x[j_B]*x[j_C]*y[j_C]+x[j_A]*x[j_C]*y[j_C]
                        -w[j_B]*w[j_C]*y[j_C]+w[j_A]*w[j_C]*y[j_C]
                        +z[j_B]^2*y[j_C]-z[j_A]*z[j_B]*y[j_C]+x[j_B]^2*y[j_C]
                        -x[j_A]*x[j_B]*y[j_C]+w[j_B]^2*y[j_C]
                        -w[j_A]*w[j_B]*y[j_C]+y[j_B]*x[j_C]^2-y[j_A]*x[j_C]^2
                        -x[j_B]*y[j_B]*x[j_C]-x[j_A]*y[j_B]*x[j_C]
                        +2*y[j_A]*x[j_B]*x[j_C]+y[j_B]*w[j_C]^2
                        -y[j_A]*w[j_C]^2-w[j_B]*y[j_B]*w[j_C]
                        -w[j_A]*y[j_B]*w[j_C]+2*y[j_A]*w[j_B]*w[j_C]
                        -y[j_A]*z[j_B]^2+z[j_A]*y[j_B]*z[j_B]
                        +x[j_A]*x[j_B]*y[j_B]+w[j_A]*w[j_B]*y[j_B]
                        -y[j_A]*x[j_B]^2-y[j_A]*w[j_B]^2) /4; # производная по координате игрек точки j_A
          f4_full[4*k1][k2] := (y[j_B]*y[j_C]*z[j_C]-y[j_A]*y[j_C]*z[j_C]+x[j_B]*x[j_C]*z[j_C]
                            -x[j_A]*x[j_C]*z[j_C]+w[j_B]*w[j_C]*z[j_C]
                            -w[j_A]*w[j_C]*z[j_C]-y[j_B]^2*z[j_C]
                            +y[j_A]*y[j_B]*z[j_C]-x[j_B]^2*z[j_C]
                            +x[j_A]*x[j_B]*z[j_C]-w[j_B]^2*z[j_C]
                            +w[j_A]*w[j_B]*z[j_C]-z[j_B]*y[j_C]^2
                            +z[j_A]*y[j_C]^2+y[j_B]*z[j_B]*y[j_C]
                            +y[j_A]*z[j_B]*y[j_C]-2*z[j_A]*y[j_B]*y[j_C]
                            -z[j_B]*x[j_C]^2+z[j_A]*x[j_C]^2
                            +x[j_B]*z[j_B]*x[j_C]+x[j_A]*z[j_B]*x[j_C]
                            -2*z[j_A]*x[j_B]*x[j_C]-z[j_B]*w[j_C]^2
                            +z[j_A]*w[j_C]^2+w[j_B]*z[j_B]*w[j_C]
                            +w[j_A]*z[j_B]*w[j_C]-2*z[j_A]*w[j_B]*w[j_C]
                            -y[j_A]*y[j_B]*z[j_B]-x[j_A]*x[j_B]*z[j_B]
                            -w[j_A]*w[j_B]*z[j_B]+z[j_A]*y[j_B]^2
                            +z[j_A]*x[j_B]^2+z[j_A]*w[j_B]^2) /4; # производная по координате зет точки j_A
        fi;
    od;
  od;

  # f5
  # loop over inner edges
  for l in [1..Length(inner_faces[2])] do
    l_v := s.faces[1][ inner_faces[2][l] ]; # список двух вершин данного ребра - исходные номера
    i := faces_ind[1][l_v[1]]; # i это номер по матрицам вершины с исходным номером l_v[1]
    j := faces_ind[1][l_v[2]]; # j аналогично
    # the edge length:
    Lij := (w[l_v[2]]-w[l_v[1]])^2 + (x[l_v[2]]-x[l_v[1]])^2 + (y[l_v[2]]-y[l_v[1]])^2 + (z[l_v[2]]-z[l_v[1]])^2 ;
    # fill out some elements
    if i<>0 then
      f5[6*i-5][4*l-3] := -(x[l_v[2]]-x[l_v[1]])/Lij;
      f5[6*i-5][4*l-2] := (w[l_v[2]]-w[l_v[1]])/Lij;
      f5[6*i-5][4*l-1] := 0;
      f5[6*i-5][4*l] := 0;
      f5[6*i-4][4*l-3] := -(y[l_v[2]]-y[l_v[1]])/Lij;
      f5[6*i-4][4*l-2] := 0;
      f5[6*i-4][4*l-1] := (w[l_v[2]]-w[l_v[1]])/Lij;
      f5[6*i-4][4*l] := 0;
      f5[6*i-3][4*l-3] := -(z[l_v[2]]-z[l_v[1]])/Lij;
      f5[6*i-3][4*l-2] := 0;
      f5[6*i-3][4*l-1] := 0;
      f5[6*i-3][4*l] := (w[l_v[2]]-w[l_v[1]])/Lij;
      f5[6*i-2][4*l-3] := 0;
      f5[6*i-2][4*l-2] := -(y[l_v[2]]-y[l_v[1]])/Lij;
      f5[6*i-2][4*l-1] := (x[l_v[2]]-x[l_v[1]])/Lij;
      f5[6*i-2][4*l] := 0;
      f5[6*i-1][4*l-3] := 0;
      f5[6*i-1][4*l-2] := -(z[l_v[2]]-z[l_v[1]])/Lij;
      f5[6*i-1][4*l-1] := 0;
      f5[6*i-1][4*l] := (x[l_v[2]]-x[l_v[1]])/Lij;
      f5[6*i][4*l-3] := 0;
      f5[6*i][4*l-2] := 0;
      f5[6*i][4*l-1] := -(z[l_v[2]]-z[l_v[1]])/Lij;
      f5[6*i][4*l] := (y[l_v[2]]-y[l_v[1]])/Lij;
    fi;
    if j<>0 then
      f5[6*j-5][4*l-3] := (x[l_v[2]]-x[l_v[1]])/Lij;
      f5[6*j-5][4*l-2] := -(w[l_v[2]]-w[l_v[1]])/Lij;
      f5[6*j-5][4*l-1] := 0;
      f5[6*j-5][4*l] := 0;
      f5[6*j-4][4*l-3] := (y[l_v[2]]-y[l_v[1]])/Lij;
      f5[6*j-4][4*l-2] := 0;
      f5[6*j-4][4*l-1] := -(w[l_v[2]]-w[l_v[1]])/Lij;
      f5[6*j-4][4*l] := 0;
      f5[6*j-3][4*l-3] := (z[l_v[2]]-z[l_v[1]])/Lij;
      f5[6*j-3][4*l-2] := 0;
      f5[6*j-3][4*l-1] := 0;
      f5[6*j-3][4*l] := -(w[l_v[2]]-w[l_v[1]])/Lij;
      f5[6*j-2][4*l-3] := 0;
      f5[6*j-2][4*l-2] := (y[l_v[2]]-y[l_v[1]])/Lij;
      f5[6*j-2][4*l-1] := -(x[l_v[2]]-x[l_v[1]])/Lij;
      f5[6*j-2][4*l] := 0;
      f5[6*j-1][4*l-3] := 0;
      f5[6*j-1][4*l-2] := (z[l_v[2]]-z[l_v[1]])/Lij;
      f5[6*j-1][4*l-1] := 0;
      f5[6*j-1][4*l] := -(x[l_v[2]]-x[l_v[1]])/Lij;
      f5[6*j][4*l-3] := 0;
      f5[6*j][4*l-2] := 0;
      f5[6*j][4*l-1] := (z[l_v[2]]-z[l_v[1]])/Lij;
      f5[6*j][4*l] := -(y[l_v[2]]-y[l_v[1]])/Lij;
    fi;
  od;


  F := rec (
    f2_full := f2_full,
    f3_FULL := f3_FULL,
    f4_full := f4_full,
    f5 := f5,
    last_inner_edge := ind_f1_li,
    last_inner_2face := ind_f2_li
  );

  return F;
end );


###########################################################################################################



