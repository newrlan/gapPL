Print("We gonna eat your brain\n");

# Define polytopes somehow

InstallValue(
Simplex1, rec(
  vertices := List(["A", "B"]),
  faces := List([
    # 1-faces
    List([Set([1,2])])
  ])  
));

InstallValue(
Simplex2, rec( 
  vertices := List(["A", "B", "C"]),
  faces := List([
    # 1-faces
    List([Set([1,2]), Set([2,3]), Set([3,1])]),
    # 2-faces
    List([Set([1,2,3])])
  ])
));

InstallValue(
Simplex3, rec( 
  vertices := List(["A", "B", "C","D"]),
  faces := List([
    # 1-faces
    List([Set([1,2]), Set([2,3]), Set([3,1]), 
          Set([4,1]), Set([4,2]), Set([4,3])]),
    # 2-faces
    List([Set([1,2,3]), Set([1,5,4]), Set([3,4,6]), Set([2,5,6])]),
    # 3-faces
    List([[1,2,3,4]])
  ])
));

InstallValue(
3sq, rec(
  vertices := List(["A", "B", "C", "D"]),
  faces := List([
    # 1-faces
    List([Set([1,2]), Set([2,3]), Set([3,4]), Set([4,1])]),
    # 2-faces
    List([[1,2,3,4],[1,2,3,4],[1,2,3,4]])
  ])
));

InstallValue(
T2, rec(
  vertices := List(["A", "B", "C", "D"]),
  faces := List([
    # 1-faces
    List([Set([1,2]), Set([2,3]), Set([3,4]), Set([4,1]),
          Set([1,2]), Set([2,3]), Set([3,4]), Set([4,1])]),
    # 2-faces
    List([Set([1,2,3,4]),Set([5,6,7,8]),Set([2,5,4,7]),Set([3,6,1,8])])
  ]),
  syms := [ [(1,2)(3,4), (1,5)(2,4)(3,7)(6,8), (1,3)(2,4)], # translation
            [(1,4)(2,3), (1,3)(2,6)(4,8)(5,7), (1,4)(2,3)], # another translation
            [(), (1,5)(2,6)(3,7)(4,8), (1,2)(3,4)], #  x -> -x (both coordinates;
                                                    # at the moment, no need to do it for each separately) 
            [(2,4), (1,4,5,8)(2,7,6,3), (1,3,2,4)] # rotation
          ]
));

InstallValue(
Cube, rec(
  vertices := List(["A", "B", "C", "D", "A'", "B'", "C'", "D'"]),
  faces := List([
    # 1-faces
    List([Set([1,2]), Set([2,3]), Set([3,4]), Set([4,1]),
          Set([1,5]), Set([2,6]), Set([3,7]), Set([4,8]),
          Set([5,6]), Set([6,7]), Set([7,8]), Set([8,5])]),
    # 2-faces
    List([Set([1,6,9,5]),Set([6,2,7,10]),Set([7,11,8,3]),Set([8,12,5,4]),
          Set([1,2,3,4]),Set([9,10,11,12])]),
    # 3-faces
    List([Set([1,2,3,4,5,6])])
  ])
));


InstallValue(
S3cubic, rec(
  vertices := List(["A", "B", "C", "D", "A'", "B'", "C'", "D'"]),
  faces := List([
    # 1-faces
    List([Set([1,2]), Set([2,3]), Set([3,4]), Set([4,1]),
          Set([1,5]), Set([2,6]), Set([3,7]), Set([4,8]),
          Set([5,6]), Set([6,7]), Set([7,8]), Set([8,5])]),
    # 2-faces
    List([Set([1,6,9,5]),Set([6,2,7,10]),Set([7,11,8,3]),Set([8,12,5,4]),
          Set([1,2,3,4]),Set([9,10,11,12])]),
    # 3-faces
    List([Set([1,2,3,4,5,6]),Set([1,2,3,4,5,6])])
  ])
));

InstallValue(
Sphere1, rec(
  vertices := ["A","B"],
  faces := [
    # 1-faces
    [[1,2], [1,2]]
  ],
  syms := [ [(1,2),(1,2)], [(),(1,2)] ]
));

InstallValue(
Sphere2, rec(
  vertices := ["A", "B", "C"],
  faces := [
    # 1-faces
    [ [1,2], [1,3], [2,3] ],
    # 2-faces
    [ [1,2,3], [1,2,3] ],
  ],
  syms := [ [(1,2,3),(1,3,2),()], [(),(),(1,2)], [(2,3),(1,2),()] ]
));
    

InstallValue(
Sphere4, rec(
  vertices := ["A","B","C","D","E"],
  faces := [
    # 1-faces
    [ [1,2], [1,3], [1,4], [1,5], [2,3], [2,4], [2,5], [3,4], [3,5], [4,5] ],
    # 2-faces
    [ [1,2,5], [1,3,6], [1,4,7], [2,3,8], [2,4,9], [3,4,10], [5,6,8], [5,7,9], [6,7,10], [8,9,10] ],
    # their vertices:
    #  1,2,3    1,2,4    1,2,5    1,3,4    1,3,5    1,4,5     2,3,4    2,3,5    2,4,5     3,4,5
    # 3-faces
    [ [1,2,4,7], [1,3,5,8], [2,3,6,9], [4,5,6,10], [7,8,9,10] ],
    # their vertices:
    #  1,2,3,4    1,2,3,5    1,2,4,5    1,3,4,5     2,3,4,5
    # 4-faces
    [ [1,2,3,4,5], [1,2,3,4,5] ]
  ]
));

InstallValue(
Disk3, rec(
  vertices := List(["A","B","C","D"]),
  faces := List([
    # 1-faces
    [ [1,2], [1,3], [1,4], [2,3], [3,4], [2,4] ],
    # 2-faces
    [ [1,2,4], [1,3,6], [2,3,5], [4,5,6], [4,5,6] ],
    # 3-faces
    [ [1,2,3,4], [1,2,3,5] ]
  ])
));

InstallValue(
3cubes, rec(
  vertices := List(["A", "B", "C", "D", "A'", "B'", "C'", "D'"]),
  faces := [
    # 1-faces
    [Set([1,2]), Set([2,3]), Set([3,4]), Set([4,1]),
          Set([1,5]), Set([2,6]), Set([3,7]), Set([4,8]),
          Set([5,6]), Set([6,7]), Set([7,8]), Set([8,5])],
    # 2-faces
    [Set([1,6,9,5]),Set([6,2,7,10]),Set([7,11,8,3]),Set([8,12,5,4]),
          Set([1,2,3,4]),Set([9,10,11,12])],
    # 3-faces
    [[1,2,3,4,5,6],[1,2,3,4,5,6],[1,2,3,4,5,6]]
    ]
));

InstallValue(
Bigon, rec(
  vertices := ["A","B"],
  faces := [
    # 1-faces
    [ [1,2], [1,2] ],
    # 2-faces
    [ [1,2] ]
  ]
));

InstallValue(
3pillow, rec(
  vertices := ["A", "B", "C", "D"],
  faces := [
    # 1-faces
    [ [1,2], [2,4], [1,3], [3,4], [1,4], [1,4] ],
    # 2-faces
    [ [3,4,5], [3,4,6], [1,2,5], [1,2,6] ],
    # 3-faces
    [ [1,2,3,4] ]
  ]
));

InstallValue(
3pillow1, rec(
  vertices := ["A", "B", "C", "D"],
  faces := [
    # 1-faces
    [ [1,2], [1,3], [2,3], [2,3], [2,4], [3,4] ],
    # 2-faces
    [ [1,2,3], [1,2,4], [3,5,6], [4,5,6] ],
    # 3-faces
    [ [1,2,3,4] ]
  ]
));

InstallValue(
Lantern, rec(
  vertices := ["A", "B"],
  faces := [
    # 1-faces
    [ [1,2], [1,2], [1,2] ],
    # 2-faces
    [ [1,2], [2,3], [1,3] ],
    # 3-faces
    [ [1,2,3] ]
  ]
));

InstallValue(
DoubleLantern, rec(
  vertices := ["A", "B"],
  faces := [ 
    # 1-faces
    [ [1,2], [1,2], [1,2], [1,2] ],
    # 2-faces
    [ [1,2], [2,3], [3,4], [1,4] ],
    # 3-faces
    [ [1,2,3,4], [1,2,3,4] ]
  ]
));

InstallValue(
MobiusBand, rec(
  vertices := ["A", "B", "C", "D"],
  faces := [
    # 1-faces
    [ [1,2], [3,4], [2,3], [2,4], [1,4], [1,3] ],
    # 2-faces
    [ [1,2,3,5], [1,2,4,6] ]
  ]
));

InstallValue(PoincareSphere, rec( 
      # вершина "O" находится в центре додекаэдра
      vertices := ["A", "B", "C", "D", "E", "O"],
      faces := [
        # 1-faces
          # сначала рёбра, идущие в центр, к вершине "O"
          [ [1,6], [5,6], [2,6], [4,6], [3,6],
            [2,6], [3,6], [4,6], [1,6], [3,6], [5,6], [1,6], [2,6], [5,6], [4,6],
            [1,6], [5,6], [2,6], [4,6], [3,6],
          # теперь рёбра, идущие по поверхности додекаэдра
            [1,5], [2,5], [2,4], [3,4], [1,3], [1,2], [4,5], [2,3], [1,4], [3,5]  ],
        # 2-faces
          # сначала 2-клетки с вершиной "O", это треугольники
          [ [1,2,21], [2,3,22], [3,4,23], [4,5,24], [1,5,25],
            [1,6,26], [2,8,27], [3,10,28], [4,12,29], [5,14,30],
            [6,7,28], [7,8,24], [8,9,29], [9,10,25], [10,11,30], 
            [11,12,21], [12,13,26], [13,14,22], [14,15,27], [6,15,23],
            [15,16,29], [7,17,30], [9,18,26], [11,19,27], [13,20,28],
            [16,17,21], [17,18,22], [18,19,23], [19,20,24], [16,20,25],
          # теперь 2-клетки на поверхности додекаэдра, это пятиугольники
            SortedList([21,22,23,24,25]),  # клетка 31: [AE,BE,BD,CD,AC]
            SortedList([25,30,27,23,26]),  # клетка 32: [AC,CE,DE,BD,AB]
            SortedList([26,28,24,27,21]),  # клетка 33: [AB,BC,CD,DE,AE]
            SortedList([25,28,22,27,29]),  # клетка 34: [AC,BC,BE,DE,AD]
            SortedList([29,23,28,30,21]),  # клетка 35: [AD,BD,BC,CE,AE]
            SortedList([26,22,30,24,29])   # клетка 36: [AB,BE,CE,CD,AD]
          ],
        # 3-faces
          [ [1,2,3,4,5,31], [5,6,10,19,20,32], [1,6,7,11,12,33], 
                  [2,7,8,13,14,34], [3,8,9,15,16,35], [4,9,10,17,18,36],
            [26,27,28,29,30,31], [14,15,23,24,28,32], [16,17,24,25,29,33],
                  [18,19,21,25,30,34], [11,20,21,22,26,35], [12,13,22,23,27,36]
          ]
        ] )
);

InstallValue(s2s2twisted, rec( 
      vertices := ["A", "B", "C", "D"],
      faces := [
        # 1-faces
          [ [1,2], [1,2], [3,4], [3,4], [1,3], [1,3], [2,4], [2,4], [2,3], [2,3], [1,4], [1,4] ],
        # 2-faces
         # сначала те 2-грани, которые входят в "базовый" тор 2*2
          [ [1,8,11], [3,6,11], [3,7,9], [1,5,9], [2,6,10], [4,8,10], [4,5,12], [2,7,12],
         # теперь две 2-грани для "верхнего" полнотория и две для "нижнего"
            [5,6], [7,8],   [5,6], [7,8],
         # теперь полноторие с меридианами ABA и CDC
            [1,2], [3,4],   # эти 2-грани имеют номера 13 и 14
         # а теперь полноторие с меридианами ADA и CBC
            [11,12], [9,10]    # эти 2-грани имеют номера 15 и 16,
                               # а в тетрадке были тоже 13 и 14
          ],
        # 3-faces
         # для верхнего полнотория
          [ [1,2,3,4,9,10], [5,6,7,8,9,10],
         # для нижнего полнотория
            [1,2,3,4,11,12], [5,6,7,8,11,12],
         # для полнотория с меридианами ABA и CDC
            [1,2,5,6,13,14], [3,4,7,8,13,14],   # эти 3-грани имеют номера 5 и 6
         # для полнотория с меридианами ADA и CBC
            [1,4,6,7,15,16], [2,3,5,8,15,16]   # а эти 3-грани имеют номера 7 и 8,
                                               # хотя в тетрадке они были 5 и 6
                                               # и, разумеется, в них 13,14 -> 15,16
          ],
        # 4-faces
         # вокруг полнотория с меридианами ABA и CDC
          [ [1,2,5,6], [3,4,5,6],
         # вокруг полнотория с меридианами ADA и CBC
            [1,2,7,8], [3,4,7,8]   # здесь по сравнению с тетрадкой 5,6 -> 7,8
          ]
          ] )
);

# Комплексная проективная плоскость.
InstallValue(cp2, rec(
  vertices := ["A", "B", "C", "D"],
  faces := [
    # 1-faces
      [ [1,2], [1,2], [3,4], [3,4], [1,3], [1,3], [2,4], [2,4], [2,3], [2,3], [1,4], [1,4] ],
    # 2-faces
     # сначала те 2-грани, которые входят в "базовый" тор 2*2
      [ [1,8,11], [3,6,11], [3,7,9], [1,5,9], [2,6,10], [4,8,10], [4,5,12], [2,7,12],
     # теперь две 2-грани для "верхнего" полнотория   ### и "нижнего" нет
        [5,6], [7,8],
     # теперь полноторие с меридианами ABA и CDC — "левое"
        [1,2], [3,4],   # эти 2-грани имеют номера 11 и 12
     # а теперь полноторие с меридианами ADA и CBC — "правое"
        [11,12], [9,10]    # эти 2-грани имеют номера 13 и 14
      ],
    # 3-faces
     # для верхнего полнотория
      [ [1,2,3,4,9,10], [5,6,7,8,9,10],
        ### нижнего полнотория нет
     # для "левого" полнотория с меридианами ABA и CDC
        [1,2,5,6,11,12], [3,4,7,8,11,12],   # эти 3-грани имеют номера 3 и 4
     # для "правого" полнотория с меридианами ADA и CBC
        [1,4,6,7,13,14], [2,3,5,8,13,14]   # а эти 3-грани имеют номера 5 и 6
      ],
    # 4-faces
     # вокруг верхнего полнотория
      [ [1,2,3,4], [1,2,5,6],
     # снизу от левого и правого — всего один шарик!
        [3,4,5,6] 
      ]
  ] )
);



##### Now functions - but KummerSurface will go *after* the function KummerSurface_ !

###############################################################################################
###############################################################################################

InstallGlobalFunction( Lens,
function (p,q)
    local i, lens_vertices, 
      edges_AC, edges_BC, edges_AD, edges_BD, edges_AB, edges_CD,
      faces_ACD_upper, faces_ACD_lower, faces_BCD_upper, faces_BCD_lower, 
        faces_ABC_left, faces_ABC_right, faces_ABD_left, faces_ABD_right,
      tetra_upper_left, tetra_upper_right, tetra_lower_left, tetra_lower_right,
      all_faces
    ;

    lens_vertices := ["A", "B", "C", "D"];

    edges_AC := List( [1..p], i -> [1,3] );
    edges_BC := List( [1..p], i -> [2,3] );
    edges_AD := List( [1..p], i -> [1,4] );
    edges_BD := List( [1..p], i -> [2,4] );
    edges_AB := [[1,2], [1,2]]; # first left, then right
    edges_CD := [[3,4], [3,4]]; # first upper, then lower
    # and then all edges will be assembled in a single list in this very order

    faces_ACD_upper := List( [1..p], i -> [i, 2*p+i, 4*p+3] );
    faces_ACD_lower := List( [1..p], i -> [1+(q+i-1) mod p, 2*p+i, 4*p+4] );
    faces_BCD_upper := List( [1..p], i -> [p+i, 3*p+i, 4*p+3] );
    faces_BCD_lower := List( [1..p], i -> [p+ 1+(q+i-1) mod p, 3*p+i, 4*p+4] );
    faces_ABC_left := List( [1..p], i -> [i, p+i, 4*p+1] );
    faces_ABC_right := List( [1..p], i -> [1+i mod p, p+i, 4*p+2] );
    faces_ABD_left := List( [1..p], i -> [2*p+i, 3*p+i, 4*p+1] );
    faces_ABD_right := List( [1..p], i -> [2*p +1+i mod p, 3*p+i, 4*p+2] );
    # and then all 2-faces will be assembled in a single list in this very order

    tetra_upper_left := List( [1..p], i -> [i, 2*p+i, 4*p+i, 6*p+i] );
    tetra_upper_right := List( [1..p], i -> [1+i mod p, 2*p+i, 5*p+i, 7*p+i] );
    tetra_lower_left := List( [1..p], i -> [p+i, 3*p+i, 4*p+ 1+(q+i-1) mod p, 6*p+i] );
    tetra_lower_right := List( [1..p], i -> [p+ 1+i mod p, 3*p+i, 5*p+ 1+(q+i-1) mod p, 7*p+i] );

    all_faces := [ 
      Concatenation( edges_AC, edges_BC, edges_AD, edges_BD, edges_AB, edges_CD ),
      Concatenation( faces_ACD_upper, faces_ACD_lower, faces_BCD_upper, faces_BCD_lower, 
        faces_ABC_left, faces_ABC_right, faces_ABD_left, faces_ABD_right ),
      Concatenation( tetra_upper_left, tetra_upper_right, tetra_lower_left, tetra_lower_right ) 
    ];

    return rec ( vertices := lens_vertices, faces := all_faces );
end );



###############################################################################################

InstallGlobalFunction( PolPrint,
# print all the faces of simplitial complex in terms of names of vertices
# (in a bit awkward fashion, but who cares)
function (s) 
  local d,d1,f,f1,f2,i,j,k,fs;

  Print("Vertices: ", s.vertices, "\n");
  
  # loop over dimension of faces (actual dimensions are d-1)
  for d in [1..Length(s.faces)] do
    # loop over faces of given dimension
    for k in [1..Length(s.faces[d])] do
      f := s.faces[d][k];
      f1 := f;
      # recursively substitute faces
      for d1 in [d-1,d-2..1] do
        f2 := [];
	for i in [1..Length(f1)] do
	  UniteSet(f2, s.faces[d1][f1[i]]); 
	od;
	f1 := f2;
      od;
      # finally, substitute faces names
      fs := [];
      for i in [1..Length(f1)] do
        UniteSet(fs, [s.vertices[f1[i]]]);
      od;
      Print("Face ", d, " #", k, " (", f, ") is ", fs, "\n");
    od;
  od;
end );


###############################################################################################

InstallGlobalFunction( ballAB,
# a ball of dimension n and just two vertices A and B
function(n)
  local l, m;
    l := List( [1..n-1], x -> [ [1,2], [1,2] ] );
    m := Concatenation( l, [[[1,2]]] );

  return rec(
    vertices := ["A", "B"],
    faces := m
  );
end );

###############################################################################################

InstallGlobalFunction( sphereAB,
# a sphere of dimension n and just two vertices A and B
function(n)
  local l;
    l := List( [1..n], x -> [ [1,2], [1,2] ] );

  return rec(
    vertices := ["A", "B"],
    faces := l
  );
end );

###############################################################################################

InstallGlobalFunction( KummerSurface, function()
    local s,isyms,sym,i0,i1,i2,lv1,lv2,i,j,k,pl,4l,i_,j_,k_,i_1,i_2,j_1,j_2,
    k_1,k_2,l1,l2,k3,k2_1,k2_2,l,x;

    s := PolProductSymsDict( T2, T2 );

    # numbers of vertices and faces in T2 x T2
    i0 := Length( s.vertices );  # the number 16
    i1 := Length( s.faces[1] );  # the number 64
    i2 := Length( s.faces[2] );

    # symmetries of T2 x T2
    isyms := StructuralCopy(s.syms);


    ### begin vertices:
    lv1 := StructuralCopy(s.vertices);
    lv2 := StructuralCopy(s.vertices);

    for i in [1..i0] do Add(lv1[i],"rho"); od;

    for i in [1..i0] do Add(lv2[i],"sigma"); od;

    s.vertices := Concatenation( lv1, lv2 );

    for i in [1..Length(isyms)] do
      pl := Concatenation( Permuted([1..i0], isyms[i][1]), Permuted([i0+1..2*i0], isyms[i][1]) );
      s.syms[i][1] := PermListList( pl,[1..2*i0] );
    od;
    ### end vertices


    ### begin edges

    # changing boundary vertices for edges of type "edge x vertex"
    for i in [1..Length(T2.faces[1])] do
      for j in [1..Length(T2.vertices)] do
        k := LookupDictionary(s.fd, [[1,i],[0,j]]);
        s.faces[1][k] := s.faces[1][k] + i0; # adds i0 to both the edge beginning and end
      od;
    od;

    # creating meridians
    # the numbers of meridians corresponding to initial vertex j will be
      # i1+j, i1+i0+j, i1+2*i0+j, i1+3*i0+j
    for i in [1..4] do
      for j in [1..i0] do
        Add( s.faces[1], [j,i0+j] );
      od;
    od;

    # creating the action of symmetries on meridians
    for i in [1..Length(s.syms)] do
      # first write down what is due to action on vertices
      4l := List( [1..4], j-> Permuted([i1+(j-1)*i0+1..i1+j*i0], isyms[i][1]) );
      # then action on meridians within each vertex
      if i=3 then 4l := Permuted(4l, (1,3)(2,4)); # reflection in the first T2
      elif i=4 then 4l := Permuted(4l, (1,2,3,4)); # rotation in the first T2
      elif i=7 then 4l := Permuted(4l, (1,3)(2,4)); # reflection in the second T2
      elif i=8 then 4l := Permuted(4l, (1,2,3,4)^-1); # rotation in the second T2
      fi;
      4l := Concatenation( 4l );
      pl := Concatenation( Permuted([1..i1], isyms[i][2]), 4l );
      s.syms[i][2] := PermListList( pl,[1..i1+4*i0] );
    od;
    ### end edges


    ### begin 2-faces

    # finding a vertex in the first 2-face of type "edge x edge"
    k_ := LookupDictionary(s.fd, [[1,1],[1,1]]); # the number 17
    i_ := s.faces[2][k][1]; # the number 1
    j_ := s.faces[1][i][1]; # again the number 1

    # putting a meridian corresponding to this vertex into the boundary of the 2-face
      # and the meridians obtained from this by symmetries
        # into the boundaries of the 2-faces obtained by the same symmetries
    for i_1 in [0..1] do for i_2 in [0..1] do # first translation in the first and second T2
      for j_1 in [0..1] do for j_2 in [0..1] do # second translation in the first and second T2
        for k_1 in [0..3] do for k_2 in [0..3] do
          l2 := Permuted( [1..i2], isyms[1][3]^i_1 * isyms[5][3]^i_2 * 
           isyms[2][3]^j_1 * isyms[6][3]^j_2 * isyms[4][3]^k_1 * isyms[8][3]^k_2 );
          l1 := Permuted( [1..Length(s.faces[1])], s.syms[1][2]^i_1 * s.syms[5][2]^i_2 * 
           s.syms[2][2]^j_1 * s.syms[6][2]^j_2 * s.syms[4][2]^k_1 * s.syms[8][2]^k_2 );
          AddSet( s.faces[2][l2[k_]], l1[i1+j_] ); # или здесь позиции? наверно, всё равно
        od; od;
      od; od;
    od; od;

    # creating sectors
    # the numbers of sectors corresponding to initial vertex j will be
      # i2+j, i2+i0+j, i2+2*i0+j, i2+3*i0+j
        # also remember that the numbers of meridians corresponding to initial vertex j are
         # i1+j, i1+i0+j, i1+2*i0+j, i1+3*i0+j
    for i in [0..3] do
      for j in [1..i0] do
        Add( s.faces[2], Set( [i1+i*i0+j, i1+((i+1) mod 4)*i0+j] ) );
      od;
    od;

    # creating the action of symmetries on sectors
    for i in [1..Length(s.syms)] do
      # first write down what is due to action on vertices
      4l := List( [1..4], j-> Permuted([i2+(j-1)*i0+1..i2+j*i0], isyms[i][1]) );
      # then action on meridians within each vertex
      if i=3 then 4l := Permuted(4l, (1,3)(2,4)); # reflection in the first T2
      elif i=4 then 4l := Permuted(4l, (1,2,3,4)); # rotation in the first T2
      elif i=7 then 4l := Permuted(4l, (1,3)(2,4)); # reflection in the second T2
      elif i=8 then 4l := Permuted(4l, (1,2,3,4)^-1); # rotation in the second T2
      fi;
      4l := Concatenation( 4l );
      pl := Concatenation( Permuted([1..i2], isyms[i][3]), 4l );
      s.syms[i][3] := PermListList( pl,[1..i2+4*i0] );
    od;
    ### end 2-faces



    ### begin 3-faces

    # "2-face x edge" - adding sector(s)
    for i in [1..Length(T2.faces[2])] do
      for k_1 in [1..Length(T2.faces[2][i])] do
        for k_2 in [1..Length(T2.faces[2][i])] do
          l1 := T2.faces[2][i][k_1]; # a number of edge entering in 2-face i
          l2 := T2.faces[2][i][k_2]; # also a number of edge entering in 2-face i
          if k_1<>k_2 and Intersection(T2.faces[1][l1], T2.faces[1][l2]) <> [] then
            for j in [1..Length(T2.faces[1])] do
              k3 := LookupDictionary(s.fd, [[2,i],[1,j]]);
              k2_1 := LookupDictionary(s.fd, [[1,l1],[1,j]]);
              k2_2 := LookupDictionary(s.fd, [[1,l2],[1,j]]);
              for l in [i2+1..Length(s.faces[2])] do
                if Intersection( s.faces[2][l], s.faces[2][k2_1] ) <> []
                  and Intersection( s.faces[2][l], s.faces[2][k2_2] ) <> [] then
                  AddSet( s.faces[3][k3], l );
                fi;
              od;
            od;
          fi;
        od;
      od;
    od;

    # "edge x 2-face" - adding sector(s)
    for i in [1..Length(T2.faces[2])] do
      for k_1 in [1..Length(T2.faces[2][i])] do
        for k_2 in [1..Length(T2.faces[2][i])] do
          l1 := T2.faces[2][i][k_1]; # a number of edge entering in 2-face i
          l2 := T2.faces[2][i][k_2]; # also a number of edge entering in 2-face i
          if k_1<>k_2 and Intersection(T2.faces[1][l1], T2.faces[1][l2]) <> [] then
            for j in [1..Length(T2.faces[1])] do
              k3 := LookupDictionary(s.fd, [[1,j],[2,i]]);
              k2_1 := LookupDictionary(s.fd, [[1,j],[1,l1]]);
              k2_2 := LookupDictionary(s.fd, [[1,j],[1,l2]]);
              for l in [i2+1..Length(s.faces[2])] do
                if Intersection( s.faces[2][l], s.faces[2][k2_1] ) <> []
                  and Intersection( s.faces[2][l], s.faces[2][k2_2] ) <> [] then
                  AddSet( s.faces[3][k3], l );
                fi;
              od;
            od;
          fi;
        od;
      od;
    od;
    ### end 3-faces


    ### symmetry

    x := PolFactorInvolution( s, [3,7] );

    Unbind(x.syms);
    Unbind(x.fd);

    return(x);

end );


###############################################################################################
#            <ManSection><Func Name="TorTwist" Arg="n" />
#                <Description>
#                    Создаем n раз скрученный тор, под действием матрицы
#                    <M> \begin{bmatrix}
#                           1 & -1 & 0 \\
#                           0 & 1  & 0 \\
#                           0 & 0  & 1
#                        \end{bmatrix}^n
#                    </M>. При положительных целых <M>n</M>.
#                </Description>
#            </ManSection>

InstallGlobalFunction(TorTwist, function(n)
    local t2, I, tower1, tower2, l2, subpol1, c, pol, subpol2, verts, k, subpol11, subpol12, tower, subpol21, subpol22, i;

    # T2    2-тор (поставляется вместе с нашей библиотекой функций)
    # Что бы осуществить склейку необходимо обеспечить наличие хотя бы одной
    # комбинаторной 2-клетки (т.е. однозначно определяющейся набором своих
    # вершин), для этого добавляем дополнительную вершину на одном из ребер.
    t2 := DivideFace(T2, [1, 1], "K");
    # t2 := DivideFace(t2, [1, 5], "L");
    # t2 := DivideFace(t2, [1, 2], "N");
    # t2 := DivideFace(t2, [1, 6], "M");
    I := ballAB(1);
    I.vertices := [0, 1];

    tower1 := PolProduct(I, t2);
    tower2 := StructuralCopy(tower1);
    l2 := Length(tower1.faces[2]);
    # Теперь на нижнем торе цилиндра нужно создать необходимое разбиение на
    # треугольники. По построению функции PolProduct индексы клеток нулевого
    # основания совпадают с индексами клеток с первого политопа в произведении.
    subpol1 := [];
    tower1 := DivideFace(tower1, [2, 1], [1, 3]);
    tower1 := DivideFace(tower1, [2, 2], [1, 3]);
    tower1 := DivideFace(tower1, [2, 3], [2, 4]);
    tower1 := DivideFace(tower1, [2, 4], [2, 4]);

    tower2 := DivideFace(tower2, [2, 1], [2, 4]);
    tower2 := DivideFace(tower2, [2, 2], [2, 4]);
    tower2 := DivideFace(tower2, [2, 3], [1, 3]);
    tower2 := DivideFace(tower2, [2, 4], [1, 3]);

    # По построению функции DivideFace, клетка разбивается на две, причем одна
    # из новых клеток заменяет существующую клетку, а вторая дописывается в
    # конец списка клеток нужной размерности. Новая добавленная клетка меньшей
    # размерности записывается в конец соответствующего списка.
    subpol1 := [1, 2, 3, 4, l2+1, l2+2, l2+3, l2+4];

    # Осталось сделать склейку двух башень. 
    # Переименуем во втором торе вершины  С <-> D.
    c := tower2.vertices[3];
    tower2.vertices[3] := tower2.vertices[4];
    tower2.vertices[4] := c;

    l2 := Length(tower1.faces[2]);
    pol := FreeUnionPol(tower1, tower2);
    for i in [1 .. Length(pol.vertices)] do
        if pol.vertices[i][2][1] = 0 then
            pol.vertices[i] := pol.vertices[i][2];
        fi;
    od;
    Print("A\n");
    pol := VerticesRullGluePol(pol, subpol1, l2 + subpol1, 2);
    Print("B\n");

    subpol1 := [];
    subpol2 := [];
    for i in [1 .. Length(pol.faces[2])] do
        verts := FaceComp(pol, [2, i]).0;
        verts := List(pol.vertices{verts}, x -> x[1]);
        if Set(verts) = [1] then
            Add(subpol1, i);
        elif Set(verts) = [2] then
            Add(subpol2, i);
        fi;
    od;
    # По построению сейчас получается взаимнооднозначное соответствие 2-клеток в
    # списках subpol1 и subpol2.
    # pol := GlueIndenticalSubpolitops(pol, subpol1, subpol2, 2);

    k := 1; # проделанное количество скруток
    l2 := Length(pol.faces[2]);
    subpol11 := StructuralCopy(subpol1);
    subpol12 := StructuralCopy(subpol2);
    tower := StructuralCopy(pol);
    while k < n do
        pol := FreeUnionPol(pol, tower);
        # Теперь необходимо пересчитать индексы subpol1 и subpol2 на дубликате.
        subpol21 := subpol1 + l2;
        subpol22 := subpol2 + l2;
        pol := GlueIndenticalSubpolitops(pol, subpol12, subpol21, 2);
        k := k + 1;
        subpol12 := StructuralCopy(subpol22 - 4);
        l2 := Length(pol.faces[2]);
        # l2 := 2 * l2 - 4;
    od;

    pol := GlueIndenticalSubpolitops(pol, subpol11, subpol12, 2);

    return pol;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="projectivePlane" Arg="n" />
#                <Description>
#                    Вещественная проективная плоскость размрености <M>n</M>
#                </Description>
#            </ManSection>

InstallGlobalFunction(projectivePlane, function(n)
    local names, vertices, simplex, delta, i, vector, pol, gluepair, pos, xs, x, sim, pair;


    names := Tuples([-1, 1], n+1);
    vertices := [];
    simplex := [];
    for xs in names do
        delta := [];
        i := 1;
        for x in xs do
            vector := List([1 .. n+1], x -> 0);
            vector[i] := xs[i];
            Add(vertices, vector);
            Add(delta, vector);
            i := i + 1;
        od;
        Add(simplex, delta);
    od;

    pol := FromSimplexToPolytope(simplex);
    gluepair := [];
    i := 1;
    for sim in simplex do
        pos := Position(simplex, - sim);
        Add(gluepair, Set([i, pos]));
        i := i + 1;
    od;
    gluepair := Set(gluepair);
    pol.vertices := List(pol.vertices, xs -> List(xs, x -> AbsoluteValue(x)));
    for pair in gluepair do
        pol := VerticesRullGlueFace(pol, pair, n);
    od;

    return pol;
end);
