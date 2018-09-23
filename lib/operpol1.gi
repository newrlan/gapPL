InstallGlobalFunction( PolBnd,
# Creates an index of boundary faces of face [d,fn] of complex s
# <result>[i] --- index of (i-1)-dimensional faces of s which are in the boundary of [d,fn].
# Input data: polytope, d, fn

#function (s, d, fn)	
function (s,adr)
 
  local i, d1, f_b, d, fn;

  d	:=adr[1];
  fn	:=adr[2];

  f_b := List([1..d], x->Set([]) );
  f_b[d] := StructuralCopy(s.faces[d][fn]);
  for d1 in [d-1,d-2..1] do
    for i in [1..Length(f_b[d1+1])] do
      UniteSet(f_b[d1], s.faces[d1][f_b[d1+1][i]]);
    od;
  od;

  return f_b;

end );

################################################################################

InstallGlobalFunction( PolCheckComb,
# check if a face of given complex is combinatorial complex
# that is, whether every its subface of lower dimension is uniquely
# determined by its vertices (and dimension)

# input: complex p, face dimension d and face number fn

function (p, adr)

  local f_b         # index of boundary faces 
    , d1
    , f1n
    , faces_in_vertices
    , faces_in_vertices_n, d, fn
    ;

	d:=adr[1];
	fn:=adr[2];
  f_b := PolBnd(p,[d,fn]);

  # check that every subface of p which is in the index f_b
  # is uniquely determined by its vertices (and dimension)

  # in a loop over dimensions
  for d1 in [2..Length(f_b)] do
    faces_in_vertices := [];
    faces_in_vertices_n := 0;
    # in a loop over faces
    for f1n in f_b[d1] do
      AddSet(faces_in_vertices, FaceComp(p,[d1-1,f1n]).0);
      faces_in_vertices_n := faces_in_vertices_n + 1;
      if Length(faces_in_vertices)<>faces_in_vertices_n then
        return false;
      fi;
    od;
  od;

  return true;

end );

################################################################################

InstallGlobalFunction( PolProduct,
# Cartesian product of two polytopes

function (s1, s2)

  local d1,d2,i,i1,i2,i1_,i2_,s,faces_dict, f;

  # init new complex
  s := rec(
    vertices := List([]),
    faces := List([1..(Length(s1.faces)+Length(s2.faces))], x->[])
  );
 
  # build the correspondense table for 2 faces (one from s1, another from s2) with their product
  # (extensively)
  faces_dict := NewDictionary([[1,1],[1,1]], true);  
  for d1 in [0..Length(s1.faces)] do
    for d2 in [0..Length(s2.faces)] do
      if d1=0 then
        i1_ := Length(s1.vertices);
      else
        i1_ := Length(s1.faces[d1]);
      fi;
      if d2=0 then
        i2_ := Length(s2.vertices);
      else 
        i2_ := Length(s2.faces[d2]);
      fi;
      for i1 in [1..i1_] do
        for i2 in [1..i2_] do
	  if d1=0 and d2=0 then
	    Add(s.vertices,0);
	    AddDictionary(faces_dict, [[d1,i1],[d2,i2]], Length(s.vertices)); 
	  else
	    Add(s.faces[d1+d2],0);
	    AddDictionary(faces_dict, [[d1,i1],[d2,i2]], Length(s.faces[d1+d2])); 
	  fi;
	od;
      od;
    od;
  od;

  # form new set of vertices
  for i1 in [1..Length(s1.vertices)] do
    for i2 in [1..Length(s2.vertices)] do
      s.vertices[LookupDictionary(faces_dict,[[0,i1],[0,i2]])] := 
        StructuralCopy([s1.vertices[i1], s2.vertices[i2]]);
    od;
  od;
 
  # compute the faces
  for d1 in [0..Length(s1.faces)] do
    for d2 in [0..Length(s2.faces)] do
      if d1=0 then
        i1_ := Length(s1.vertices);
      else
        i1_ := Length(s1.faces[d1]);
      fi;
      if d2=0 then
        i2_ := Length(s2.vertices);
      else 
        i2_ := Length(s2.faces[d2]);
      fi;
      if d1>0 or d2>0 then
        for i1 in [1..i1_] do
          for i2 in [1..i2_] do
	    f := Set([]);
	    # boundaries of s1 by whole s2
	    if d1>0 then
	      for i in s1.faces[d1][i1] do
	        AddSet(f, LookupDictionary(faces_dict, [[d1-1,i],[d2,i2]])); 
	      od;
	    fi;
	    # whole s1 by boundaries of s2
	    if d2>0 then
	      for i in s2.faces[d2][i2] do
	        AddSet(f, LookupDictionary(faces_dict, [[d1,i1],[d2-1,i]])); 
	      od;
	    fi;
	    # write down new face
	    s.faces[d1+d2][LookupDictionary(faces_dict, [[d1,i1],[d2,i2]])] := f;
	  od;
        od;
      fi;
    od; 
  od;

  return s;
end );

#############################################################################

InstallGlobalFunction( PolProductSyms,
# Cartesian product of two polytopes with symmetries of multipliers transferred to it.
# First go the symmetries of the first multiplier, then - the second

function (s1, s2)

  local d1,d2,i,i1,i2,i1_,i2_,s,faces_dict, f,
##### new locals and other stuff added by me - I.K.
l_, k, s1syms, s2syms, sym_l1, sym_l2, sym_l, s1sym_l, s2sym_l, j1, j2;

# The list of generating symmetries is arranged according to the following principle:
# syms [ symmetry_number ][ face_dimension+1 ] = permutation_of_faces

if IsBound(s1.syms) then s1syms := ShallowCopy(s1.syms); else s1syms := []; fi;
if IsBound(s2.syms) then s2syms := ShallowCopy(s2.syms); else s2syms := []; fi;
#####


  # init new complex
  s := rec(
    vertices := List([]),
    faces := List([1..(Length(s1.faces)+Length(s2.faces))], x->[]),
#####
    syms := List( [1..Length(s1syms)+Length(s2syms)], i -> List( [1..Length(s1.faces)+Length(s2.faces)+1], j -> [] ) ),
#####
  );
 
  # build the correspondense table for 2 faces (one from s1, another from s2) with their product
  # (extensively)
  faces_dict := NewDictionary([[1,1],[1,1]], true);  
  for d1 in [0..Length(s1.faces)] do
    for d2 in [0..Length(s2.faces)] do
      if d1=0 then
        i1_ := Length(s1.vertices);
      else
        i1_ := Length(s1.faces[d1]);
      fi;
      if d2=0 then
        i2_ := Length(s2.vertices);
      else 
        i2_ := Length(s2.faces[d2]);
      fi;
      for i1 in [1..i1_] do
        for i2 in [1..i2_] do
	  if d1=0 and d2=0 then
	    Add(s.vertices,0);
	    AddDictionary(faces_dict, [[d1,i1],[d2,i2]], Length(s.vertices)); 
	  else
	    Add(s.faces[d1+d2],0);
	    AddDictionary(faces_dict, [[d1,i1],[d2,i2]], Length(s.faces[d1+d2])); 
	  fi;
	od;
      od;
    od;
  od;

  # form new set of vertices
  for i1 in [1..Length(s1.vertices)] do
    for i2 in [1..Length(s2.vertices)] do
      s.vertices[LookupDictionary(faces_dict,[[0,i1],[0,i2]])] := 
        StructuralCopy([s1.vertices[i1], s2.vertices[i2]]);
    od;
  od;
 
  # compute the faces
  for d1 in [0..Length(s1.faces)] do
    for d2 in [0..Length(s2.faces)] do
      if d1=0 then
        i1_ := Length(s1.vertices);
      else
        i1_ := Length(s1.faces[d1]);
      fi;
      if d2=0 then
        i2_ := Length(s2.vertices);
      else 
        i2_ := Length(s2.faces[d2]);
      fi;
      if d1>0 or d2>0 then
        for i1 in [1..i1_] do
          for i2 in [1..i2_] do
	    f := Set([]);
	    # boundary faces of  s1  multiplied by the whole  s2
	    if d1>0 then
	      for i in s1.faces[d1][i1] do
	        AddSet(f, LookupDictionary(faces_dict, [[d1-1,i],[d2,i2]])); 
	      od;
	    fi;
	    # the whole  s1  multiplied by boundary faces of  s2
	    if d2>0 then
	      for i in s2.faces[d2][i2] do
	        AddSet(f, LookupDictionary(faces_dict, [[d1,i1],[d2-1,i]])); 
	      od;
	    fi;
	    # write down new face
	    s.faces[d1+d2][LookupDictionary(faces_dict, [[d1,i1],[d2,i2]])] := f;
	  od;
        od;
      fi;
    od; 
  od;




#####

# preparing permuted lists of vertices/faces for complexes  s1  and  s2,
# for each face dimension  d1  or  d2  and for each face number  i
# These lists are organized according to the following principle:
# sym_l [ symmetry_number ][ face_dimension+1 ] = permuted_list_of_faces

s1sym_l := List( [1..Length(s1syms)], i -> List( [1..Length(s1.faces)+1], j -> [] ) );
s2sym_l := List( [1..Length(s2syms)], i -> List( [1..Length(s2.faces)+1], j -> [] ) );

for i in [1..Length(s1syms)] do
  for d1 in [0..Length(s1.faces)] do
    if d1=0 then
      i1_ := Length(s1.vertices);
    else
      i1_ := Length(s1.faces[d1]);
    fi;
    s1sym_l[i][d1+1] := Permuted( [1..i1_], s1syms[i][d1+1] );
  od;
od;

for i in [1..Length(s2syms)] do
  for d2 in [0..Length(s2.faces)] do
    if d2=0 then
      i2_ := Length(s2.vertices);
    else
      i2_ := Length(s2.faces[d2]);
    fi;
    s2sym_l[i][d2+1] := Permuted( [1..i2_], s2syms[i][d2+1] );
  od;
od;


# now preparing permuted lists of vertices/faces for complex  s,
# using dictionary  faces_dict

sym_l1 := List( [1..Length(s1syms)], i -> List( [1..Length(s1.faces)+Length(s2.faces)+1], j -> [] ) );
sym_l2 := List( [1..Length(s2syms)], i -> List( [1..Length(s1.faces)+Length(s2.faces)+1], j -> [] ) );

  for d1 in [0..Length(s1.faces)] do
    for d2 in [0..Length(s2.faces)] do
      if d1=0 then
        i1_ := Length(s1.vertices);
      else
        i1_ := Length(s1.faces[d1]);
      fi;
      if d2=0 then
        i2_ := Length(s2.vertices);
      else 
        i2_ := Length(s2.faces[d2]);
      fi;
        for i1 in [1..i1_] do
          for i2 in [1..i2_] do
            l_ := LookupDictionary( faces_dict, [[d1,i1],[d2,i2]] );
            for i in [1..Length(s1syms)] do
              j1 := s1sym_l[i][d1+1][i1];
              sym_l1[i][d1+d2+1][l_] := LookupDictionary( faces_dict, [[d1,j1],[d2,i2]] ); # here  i  numbers the symmetries of  s1
            od;
            for i in [1..Length(s2syms)] do
              j2 := s2sym_l[i][d2+1][i2];
              sym_l2[i][d1+d2+1][l_] := LookupDictionary( faces_dict, [[d1,i1],[d2,j2]] ); # here  i  numbers the symmetries of  s2
            od;
          od;
        od;
    od;
  od;

sym_l := Concatenation( sym_l1, sym_l2 );


# making the final lists of permutations

for i in [1..Length(sym_l)] do
  for k in [1..Length(sym_l[i])] do
s.syms[i][k] := PermListList( [1..Length(sym_l[i][k])], sym_l[i][k] );
  od;
od;

#####

  return s;
end );


################################################################################

InstallGlobalFunction( PolProductSymsDict,

# Cartesian product of two polytopes with symmetries of multipliers transferred to it.
# First go the symmetries of the first multiplier, then - the second.
# Also returns the face dictionary.

function (s1, s2)

  local d1,d2,i,i1,i2,i1_,i2_,s,faces_dict, f,
##### new locals and other stuff added by me - I.K.
l_, k, s1syms, s2syms, sym_l1, sym_l2, sym_l, s1sym_l, s2sym_l, j1, j2;

# The list of generating symmetries is arranged according to the following principle:
# syms [ symmetry_number ][ face_dimension+1 ] = permutation_of_faces

if IsBound(s1.syms) then s1syms := ShallowCopy(s1.syms); else s1syms := []; fi;
if IsBound(s2.syms) then s2syms := ShallowCopy(s2.syms); else s2syms := []; fi;
#####


  # init new complex
  s := rec(
    vertices := List([]),
    faces := List([1..(Length(s1.faces)+Length(s2.faces))], x->[]),
#####
    syms := List( [1..Length(s1syms)+Length(s2syms)], i -> List( [1..Length(s1.faces)+Length(s2.faces)+1], j -> [] ) ),
#####
  );
 
  # build the correspondense table for 2 faces (one from s1, another from s2) with their product
  # (extensively)
  faces_dict := NewDictionary([[1,1],[1,1]], true);  
  for d1 in [0..Length(s1.faces)] do
    for d2 in [0..Length(s2.faces)] do
      if d1=0 then
        i1_ := Length(s1.vertices);
      else
        i1_ := Length(s1.faces[d1]);
      fi;
      if d2=0 then
        i2_ := Length(s2.vertices);
      else 
        i2_ := Length(s2.faces[d2]);
      fi;
      for i1 in [1..i1_] do
        for i2 in [1..i2_] do
	  if d1=0 and d2=0 then
	    Add(s.vertices,0);
	    AddDictionary(faces_dict, [[d1,i1],[d2,i2]], Length(s.vertices)); 
	  else
	    Add(s.faces[d1+d2],0);
	    AddDictionary(faces_dict, [[d1,i1],[d2,i2]], Length(s.faces[d1+d2])); 
	  fi;
	od;
      od;
    od;
  od;

  # form new set of vertices
  for i1 in [1..Length(s1.vertices)] do
    for i2 in [1..Length(s2.vertices)] do
      s.vertices[LookupDictionary(faces_dict,[[0,i1],[0,i2]])] := 
        StructuralCopy([s1.vertices[i1], s2.vertices[i2]]);
    od;
  od;
 
  # compute the faces
  for d1 in [0..Length(s1.faces)] do
    for d2 in [0..Length(s2.faces)] do
      if d1=0 then
        i1_ := Length(s1.vertices);
      else
        i1_ := Length(s1.faces[d1]);
      fi;
      if d2=0 then
        i2_ := Length(s2.vertices);
      else 
        i2_ := Length(s2.faces[d2]);
      fi;
      if d1>0 or d2>0 then
        for i1 in [1..i1_] do
          for i2 in [1..i2_] do
	    f := Set([]);
	    # boundary faces of  s1  multiplied by the whole  s2
	    if d1>0 then
	      for i in s1.faces[d1][i1] do
	        AddSet(f, LookupDictionary(faces_dict, [[d1-1,i],[d2,i2]])); 
	      od;
	    fi;
	    # the whole  s1  multiplied by boundary faces of  s2
	    if d2>0 then
	      for i in s2.faces[d2][i2] do
	        AddSet(f, LookupDictionary(faces_dict, [[d1,i1],[d2-1,i]])); 
	      od;
	    fi;
	    # write down new face
	    s.faces[d1+d2][LookupDictionary(faces_dict, [[d1,i1],[d2,i2]])] := f;
	  od;
        od;
      fi;
    od; 
  od;

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# preparing permuted lists of vertices/faces for complexes  s1  and  s2,
# for each face dimension  d1  or  d2  and for each face number  i
# These lists are organized according to the following principle:
# sym_l [ symmetry_number ][ face_dimension+1 ] = permuted_list_of_faces

s1sym_l := List( [1..Length(s1syms)], i -> List( [1..Length(s1.faces)+1], j -> [] ) );
s2sym_l := List( [1..Length(s2syms)], i -> List( [1..Length(s2.faces)+1], j -> [] ) );

for i in [1..Length(s1syms)] do
  for d1 in [0..Length(s1.faces)] do
    if d1=0 then
      i1_ := Length(s1.vertices);
    else
      i1_ := Length(s1.faces[d1]);
    fi;
    s1sym_l[i][d1+1] := Permuted( [1..i1_], s1syms[i][d1+1] );
  od;
od;

for i in [1..Length(s2syms)] do
  for d2 in [0..Length(s2.faces)] do
    if d2=0 then
      i2_ := Length(s2.vertices);
    else
      i2_ := Length(s2.faces[d2]);
    fi;
    s2sym_l[i][d2+1] := Permuted( [1..i2_], s2syms[i][d2+1] );
  od;
od;


# now preparing permuted lists of vertices/faces for complex  s,
# using dictionary  faces_dict

sym_l1 := List( [1..Length(s1syms)], i -> List( [1..Length(s1.faces)+Length(s2.faces)+1], j -> [] ) );
sym_l2 := List( [1..Length(s2syms)], i -> List( [1..Length(s1.faces)+Length(s2.faces)+1], j -> [] ) );

  for d1 in [0..Length(s1.faces)] do
    for d2 in [0..Length(s2.faces)] do
      if d1=0 then
        i1_ := Length(s1.vertices);
      else
        i1_ := Length(s1.faces[d1]);
      fi;
      if d2=0 then
        i2_ := Length(s2.vertices);
      else 
        i2_ := Length(s2.faces[d2]);
      fi;
        for i1 in [1..i1_] do
          for i2 in [1..i2_] do
            l_ := LookupDictionary( faces_dict, [[d1,i1],[d2,i2]] );
            for i in [1..Length(s1syms)] do
              j1 := s1sym_l[i][d1+1][i1];
              sym_l1[i][d1+d2+1][l_] := LookupDictionary( faces_dict, [[d1,j1],[d2,i2]] ); 
              # here  i  numbers the symmetries of  s1
            od;
            for i in [1..Length(s2syms)] do
              j2 := s2sym_l[i][d2+1][i2];
              sym_l2[i][d1+d2+1][l_] := LookupDictionary( faces_dict, [[d1,i1],[d2,j2]] ); 
              # here  i  numbers the symmetries of  s2
            od;
          od;
        od;
    od;
  od;

sym_l := Concatenation( sym_l1, sym_l2 );


# making the final lists of permutations

for i in [1..Length(sym_l)] do
  for k in [1..Length(sym_l[i])] do
s.syms[i][k] := PermListList( [1..Length(sym_l[i][k])], sym_l[i][k] );
  od;
od;
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  # and now adding the dictionary to the record
  s.fd := faces_dict;

  return s;
end );



################################################################################


InstallGlobalFunction( PolTriangulate,
# triangulating a polytope ( = ball complex)
function (p) 

  local s, 
    d, d1, fn, f, f_b, v,
    new_faces, f_found,
    n_a_f,
    i, j, k,
    f_new_ind,
    n_f, newV, newF;

  s := StructuralCopy(p);

  # loop over dimensions of s (we could start at 2 BTW)
  for d in [1..Length(s.faces)] do
    # loop over faces of dimension d
    for fn in [1..Length(s.faces[d])] do
      f := s.faces[d][fn];
      
      # 0. Find mimimal vertex
      v := Minimum(FaceComp(s,[ d, fn]).0);

      # compute n_f --- number of faces "against" (i.e., not touching) v
      n_f := 0;
      if d=1 then 
        n_f := Length(f)-1;
      else
        for i in f do
          if not (v in FaceComp(s,[d-1,i]).0) then
            n_f := n_f + 1;
          fi;
        od;
      fi;  

      # if [d,fn] is not combinatorial complex, add a new vertex in the center
      if (not PolCheckComb(s,[d,fn])) then
        
        # add new vertex
        i := 1;
        while Concatenation("V", String(i)) in s.vertices do
          i := i + 1;
        od;
        Add(s.vertices, Concatenation("V", String(i)));
        newV := Length(s.vertices);

        # "connect" existing boundary faces of f to a new vertex
        # starting with dimension 0
        f_b := PolBnd(s,[d,fn]);
        # to every element from f_b there corresponds a new face of a bigger dimension having newV as vertex;
        # we will use a sort of index (face of s -> new face):
        f_new_ind := [List(s.vertices, x->0)];
        for i in [1..d] do
          Add(f_new_ind, List([1..Length(s.faces[i])], x->0));
        od;

        for i in [1..Length(f_b)] do
          for j in [1..Length(f_b[i])] do
            if i=1 then
              # simple case - add a line segment
              Add(s.faces[i], Set([f_b[i][j],newV]));
            else
              # init index of facets of new face
              newF := [f_b[i][j]];
              # find those facets among already existing faces of s
              for k in s.faces[i-1][f_b[i][j]] do
                AddSet(newF, f_new_ind[i-1][k]);
              od;
              # add new face when i>1
              Add(s.faces[i], StructuralCopy(newF));
            fi;
            # register new face in the index
            f_new_ind[i][f_b[i][j]] := Length(s.faces[i]);
          od;
        od;

        # Erase [d,fn] from the history
        # First, find out the faces f has been divided into
        newF := Difference(Set(f_new_ind[d]),[0]);
	s.faces[d][fn] := [];
	if Length(s.faces)>d then 
	  for i in [1..Length(s.faces[d+1])] do
	    if fn in s.faces[d+1][i] then
	      RemoveSet(s.faces[d+1][i], fn);
	      UniteSet(s.faces[d+1][i], newF);
	    fi;
	  od;
        fi;
      
      # if [d,fn] IS combinatorial, and n_f (number of faces with dimension d-1 not touching v) is greater than 1, use
      # "usual" triangulation
      elif n_f > 1 then
        # Normal triangulation of a big face

	# 1. Create f_b --- an index of boundaries of subcomplex [d,fn]
	# f_b[i] --- index of (i-1)-dimensional faces of s which are boundaries of [d,fn]
	f_b := PolBnd(s,[d,fn]);
	Append(f_b, [[]]);
	
	# 2. Create (initially empty) index of newly added faces
	new_faces := List([1..d+1],x->[]);

	# 3. In the loop over lower dimensions
	for d1 in [0..d-1] do
	  # In the loop over faces-boundaries of [d,fn] of dimension d1 not touching v
	  for i in [1..Length(f_b[d1+1])] do
	    if (d1=0 and f_b[d1+1][i]<>v) or 
	       (d1>0 and  not ( v in FaceComp(s,[d1,f_b[d1+1][i]]).0)) then
	      # look for face of dimension d1+1, having [d1, f_b[d1+1][i]] and v as boundaries
	      f_found := false;
	      # among boundaries of [d,fn] and among newly added faces
	      for j in Union(f_b[d1+2], new_faces[d1+1]) do
	      	if (f_b[d1+1][i] in s.faces[d1+1][j])
		  and (v in FaceComp(s,[d1+1,j]).0) then

		  # found!
		  f_found := true;
		  break;
		fi;
	      od;
	      
	      if f_found=false then
	        # If there's no such face, then add it. All needed boundaried should already be already added to s
	        n_a_f := [f_b[d1+1][i]];
		if d1=0 then
		  AddSet(n_a_f, v);
	        else 
		  # loop over boundaries of [d1,f_b[d1+1][i]]
	          for k in s.faces[d1][f_b[d1+1][i]] do
	            # look for face of dimension d1 having [d1-1, k] and [0,v] as boundaries 
		    for j in Union(f_b[d1+1], new_faces[d1]) do
	      	      if (k in s.faces[d1][j]) and 
		         (v in FaceComp(s,[d1,j]).0) then
		        # found!
		        AddSet(n_a_f, j);
		      fi;
	            od;
	          od;
		fi;
		# add n_a_f (newly added face) to s
		Append(s.faces[d1+1], [n_a_f]);
		# register n_a_f in index of newly added faces
		AddSet(new_faces[d1+1], Length(s.faces[d1+1]));
	      fi;

	    fi;
	  od;
	od;

	# 4. Erase [d,fn] from the history
	# Remove(s.faces[d],fn);
	s.faces[d][fn] := [];
	if Length(s.faces)>d then 
	  for i in [1..Length(s.faces[d+1])] do
	    if fn in s.faces[d+1][i] then
	      #RemoveSet(s.faces[d+1][i], fn); # какая то огреха с mutable\immutable над множетсвами
	      s.faces[d+1][i]:=Difference(s.faces[d+1][i], [fn]);
	      UniteSet(s.faces[d+1][i], new_faces[d]);
	    fi;
	  od;
        fi;
      fi;
    od;
  od;

  # Cleanup - remove deleted faces
  for d in [1..Length(s.faces)] do
    fn := 1;
    while fn<=Length(s.faces[d]) do
      if s.faces[d][fn]=[] then
        Remove(s.faces[d],fn);
	if d<Length(s.faces) then
	  for i in [1..Length(s.faces[d+1])] do
	    for j in [1..Length(s.faces[d+1][i])] do
	      if s.faces[d+1][i][j]>fn then
	        s.faces[d+1][i][j] := s.faces[d+1][i][j]-1;
	      fi;
	    od;
	  od;
	fi;
      else
        fn := fn + 1;
      fi;
    od;
  od;

  return s;

end );


################################################################################

InstallGlobalFunction( OrientTriangulated,
# Function computes consistent orientation on simplices of greatest dimension
# on a given triangulated complex s

# Returns array of -1,1-s which correspond to the orientation of simplices 
# of greatest dimension of s

# probably quite efficiently, but there is still space for further optimization

function(s)
  local orient,d,
    con, f_con_sort,
    i,j,k,i_k,j_k,
    i_v,j_v,k_v,
    i_s,j_s;

  # dimension
  d := Length(s.faces);
  # init the array of orientations
  orient := List([1..Length(s.faces[d])], x->0);

  # create an index of connectors
  # "connector" is a list of 
  # [is_satisfied, i, j, k, i_k, j_k]
  # where is_satisfied is 0 if connector is not satisfied,
  # and 1 otherwise,
  # i,j are indices of d-faces, k is index of their common (d-1)-face
  # i_k and j_k are signs with which k is included into boundaries of
  # i and j
  con := [];
  for i in [1..Length(s.faces[d])] do
    for j in [(i+1)..Length(s.faces[d])] do
      for k in Intersection(s.faces[d][i],s.faces[d][j]) do
        # vertices of involved simplices
        i_v := FaceComp(s,[d,i]).0;
        j_v := FaceComp(s,[d,j]).0;
        k_v := FaceComp(s,[d-1,k]).0;
        # posisions of removed vertices in i-th and j-th d-faces minus 1 modulo 2 
        i_s := RemInt(Position(i_v,Difference(i_v,k_v)[1])-1,2);
        j_s := RemInt(Position(j_v,Difference(j_v,k_v)[1])-1,2);
        Add(con, [0,i,j,k,(-1)^i_s, (-1)^j_s]);
      od;
    od;
  od;

  # now work on the list of "connectors" trying to satisfy the all
  orient[1] := 1;

  f_con_sort := function(x,y) 
         return x[1]<y[1] 
                or 
                ( x[1]=y[1] 
                  and (AbsInt(orient[x[2]]) + AbsInt(orient[x[3]]) > 
                       AbsInt(orient[y[2]]) + AbsInt(orient[y[3]]))
                );
  end;
 

  while true do
    # sort connectors making not satisfied and having at least one face oriented
    # go first
    Sort(con, f_con_sort);
    # check if all connectors are satisfied
    if con[1][1]=1 then
      break;
    fi;
    # try to satisfy the very first connector
    if orient[con[1][2]]=0 and orient[con[1][3]]=0 then
      orient[con[1][2]] := 1;
      if con[1][5]=con[1][6] then
        orient[con[1][3]] := -1;
      else
        orient[con[1][3]] := 1;
      fi;
      con[1][1] := 1;
    elif orient[con[1][2]]<>0 and orient[con[1][3]]=0 then
      if con[1][5]=con[1][6] then
        orient[con[1][3]] := - orient[con[1][2]];
      else
        orient[con[1][3]] := orient[con[1][2]];
      fi;
      con[1][1] := 1;
    elif orient[con[1][2]]=0 and orient[con[1][3]]<>0 then
      if con[1][5]=con[1][6] then
        orient[con[1][2]] := - orient[con[1][3]];
      else
        orient[con[1][2]] := orient[con[1][3]];
      fi;
      con[1][1] := 1;
    else # both are not zero
      if (con[1][5]=con[1][6] and orient[con[1][2]]=-orient[con[1][3]])
         or
         (con[1][5]<>con[1][6] and orient[con[1][2]]=orient[con[1][3]])
         then
        con[1][1] := 1;
      else
        # complex is not orientable at all
        return [];
      fi;
    fi;
  od;

  return orient;
end );


################################################################################


InstallGlobalFunction( PolInnerFaces,
# build index of inner faces of given polytope complex
# (!) <returned value>[i] --- set of inner faces of dimension (i-1)

# Any face is outer if it has at most 1 adjacent face of higher dimension
# or if it lies in the boundary of such a face.
# And inner faces are not outer faces.

function(p)

  local i,j,k,ni,
    outer_faces,
    d,d1,d2,a,ib,
    inner_faces;

  # dimension of p
  d := Length(p.faces);

  # outer_faces[i] - set of outer faces of dimension i-1
  outer_faces := List([1..d], x->Set([]));

  
  # look for isolated vertices
  outer_faces[1] := Difference([1..Length(p.vertices)], Union(p.faces[1]));

  # loop over dimension
  for d1 in [0..(d-1)] do
    # loop over faces of dimension d1
    if d1=0 then 
      ni := Length(p.vertices); 
    else
      ni := Length(p.faces[d1]);
    fi;
    for i in [1..ni] do
      # check if the face has at most 1 adjacent face of dimension d1+1
      a := 0;
      for j in p.faces[d1+1] do
        if i in j then
          a := a + 1;
        fi;
      od;
      if a<=1 then
        # Print("Adding face [", d1,", ", i, "]\n");
        # add [d1,i] and all its bounaries to the list of outer faces
        AddSet(outer_faces[d1+1],i);
        # get index of bounaries of [d1,i]
        if d1>0 then
          ib := PolBnd(p,[d1,i]);
          # add all those to outer_faces
          for k in [1..Length(ib)] do
            UniteSet(outer_faces[k],ib[k]);
          od;
        fi;
      fi;
    od;
  od;

  # finally, compute inner faces as not-outer faces
  inner_faces := List([1..d], x->[]);
  for i in [1..d] do
    if i=1 then
      j := Length(p.vertices);
    else
      j := Length(p.faces[i-1]);
    fi;
    inner_faces[i] := Difference([1..j], outer_faces[i]);
  od;

  return inner_faces;
end );

################################################################################

InstallGlobalFunction( PolDoubleCone,
# Make a double cone with vertices V1 and V2 over the given polytope p

function (p)	
 
  local pdc, d, j, k, n, nk1, nk2, v, pfk1;

  d := Length( p.faces );
  v := Length( p.vertices );

  n := [];
  for k in [1..d] do
    n[k] := Length( p.faces[k] );
  od;

  # a rough piece to be made into the double cone
  pdc := StructuralCopy( p );

  # V1 and V2 will be our double cone vertices
  # they go *after* the already existing vertices
  Append( pdc.vertices, ["V1", "V2"] );

  # appending one more dimension to the new complex
  Append( pdc.faces, [[]] );

  # below nk1 and nk2 are, essentially, n[k-1] and n[k-2], respectively,
  # with `n[0]' being the number v of vertices and `n[-1]' being just 1
  # because there exists just one cell of dimension -1  :)

  # similarly, pfk1 is, essentially, p.faces[k-1],
  # with `p.faces[0]' being [ [1], [1], ... , [1] ]

  for k in [1..d+1] do
    if k=1 then nk1:=v; nk2:=1; pfk1:=List( [1..v], x->[1] );
      elif k=2 then nk1:=n[1]; nk2:=v; pfk1:=p.faces[1];
      else nk1:=n[k-1]; nk2:=n[k-2]; pfk1:=p.faces[k-1];
    fi;
    
    for j in [1..nk1] do
     Append( pdc.faces[k], [ Concatenation( [j], List( pfk1[j], x->x+nk1 ) ) ] );
    od;
    
    for j in [1..nk1] do
     Append( pdc.faces[k], [ Concatenation( [j], List( pfk1[j], x->x+nk1+nk2 ) ) ] );
    od;
  od;
    
  return pdc;

end );

################################################################################

InstallGlobalFunction( MaxTree,
# finds a maximal tree in the 1-skeleton of a polytope as a list of edges

function(p)
         local l, n, # l-список ребер, n-количество вершин
               maxd, # максимальное дерево
               zv,   # звезда множества вершин
               vzv,  # вершины звезды
               vd,   # вершины дерева
               v,    # вершины графа которые еще не входят в дерево
               i, t,
               s,    # счетчики
               rebra; # номера ребер входящих в дерево
          l:=p.faces[1];
          n:= Length(p.vertices); # изъятие значений из триангуляции
          maxd:=[];
          zv:=[];
          vzv:=[];
          vd:=[1];
          v:=[2..n];
          rebra:=[];
   while Length(v)<>0 do
         for t in vd do
             for i in v do
                 if Set([t,i]) in l then
                    Add(zv, Set([t,i])); # создание звезды
                    Add(vzv,i); # создание вершин звезды
                 fi;
             od;
          od;
          Add(vd,vzv[1]); # добавление к проверенным вершинам
          v:=Difference(v,vd);
          Add(maxd,zv[1]);
          zv:=[];
          vzv:=[]; # обнуление параметров звезды
   od;
   for s in maxd do     # сопоставление ребру из maxd номера ребра из l
         Add(rebra,Position(l,s));
   od;
return rebra;
end );

################################################################################

InstallGlobalFunction( CellOrient,
# ОПИСАНИЕ:
# программа которая считает ориентации с которыми входят клетки размерности (k-1) в k-грани.
# входные данные: p-политоп.
# выходные данные: список ориентаций граней
# зависимости: нет.

function(p)
     local n,s,orient,i,k,ind,del,a,t,j,d,da,dj;

n:=Length(p.faces);

# 1) ориентируем 1-грани стандартным способом, выбрав направление от вершины
#    с меньшим индексом к вершине с большим индексом.

s:=1;
orient:=[[]];
for i in p.faces[1] do
    p.faces[1][s]:=Set(i);
#    Print("CellOrient orders vertices in edges, which is probably unnecessary","\n");
    orient[1][s]:=[-1,+1];
    s:=s+1;
od;

# 2) ориентируем грани размерности k=2,...,n.

for k in [2..n] do
    orient[k]:=[];
    s:=1;
    for i in p.faces[k] do
        orient[k][s]:=[];
        orient[k][s][1]:=-1; # здесь был +1 и было написано:
                             # "в ориентации данной клетки мы можем выбрать, произвольно, и -1".
                             # Я и поменял на -1, чтобы красивее получалось для ballAB и sphereAB - И.К.
        ind:=[2..Length(i)];
        a:=[1];
        t:=1;
        while ind<>[] do
            del:=[];
            for j in ind do
                d:=Intersection(p.faces[k-1][i[a[t]]],p.faces[k-1][i[j]]);
                if d<>[] then
                   d:=d[1];   # выбираем одну из списка пересечения
                   da:=Position(p.faces[k-1][i[a[t]]],d);# определяем ее позицию в первой (k-1)-мерной грани 
                   dj:=Position(p.faces[k-1][i[j]],d); # определяем ее позицию в j-той (k-1)-мерной грани
                   Add(del,j); # эти индексы мы потом удалим, т.к. их ориентации мы уже можем вычислить
                   Add(a,j);# 
                   orient[k][s][j]:=(-1)*(orient[k-1][i[j]][dj])*(orient[k-1][i[a[t]]][da])*(orient[k][s][a[t]]); 
                    # считаем ориентацию грани
                fi; # исходя из 
            od;
            ind:=Difference(ind,del);
            t:=t+1;
         od;
         s:=s+1;
    od;
od;


return orient;
end );


################################################################################

InstallGlobalFunction( PolOrient,
# ОПИСАНИЕ:
# программа которая вычисляет является ли многообразие ориентируемым
# входные данные: p-политоп.
# выходные данные: если да - получаем ориентацию n-граней,
#                  если нет - fail
# dependencies: CellOrient

function(p)
     local n,e,w,E,m,k,del,int,b,ab,i,s,j;

n:=Length(p.faces);
e:= CellOrient(p);

# 1) посчитаем предварительно ориентации каждой n-грани 
w:=[2..Length(p.faces[n])];
E:=[1];
m:=[1];
k:=1;
while w<>[] do
   del:=[];
   for i in w do
      int:=Intersection(p.faces[n][m[k]],p.faces[n][i]);
      if int <>[] then
         b:=int[1];
         E[i]:=(-1)*e[n][m[k]][Position(p.faces[n][m[k]],b)]*e[n][i][Position(p.faces[n][i],b)]*E[m[k]];
         Add(m,i);
         Add(del,i);
      fi;
   od;
   w:=Difference(w,del);
   k:=k+1;
od;

# 2) проверим, не возникает ли противоречий
i:=1;
while (i < Length(p.faces[n-1])+1) and (E<>[]) do
   ab:=[];
   s:=1;
   for j in p.faces[n] do
      if i in j then Add(ab,s); fi;
      s:=s+1;
   od;
   if Length(ab)=2 then
     if e[n][ab[1]][Position(p.faces[n][ab[1]],i)]*E[ab[1]] <> 
              (-1)*e[n][ab[2]][Position(p.faces[n][ab[2]],i)]*E[ab[2]]
         then E:=[];
      fi;
   fi;
   i:=i+1;
od;

if E=[] then 
  E := fail;
fi;

return E;

end );


################################################################################

InstallGlobalFunction( Slovo,
function( A, B )    # A  is a list of edges surrounding the 2-face, each given as the set of its 2 vertices
                    # B  is the list of corresponding group elements
       local S,     # последовательность вершин
             Ac,    # копия списка А, чтобы не повредить сам список А
             Sor,   # ориентации граней
             i,k,   # локальные переменные
             oo,    # порядок рёбер в направлении обхода
             Aoo,   # перестановка
             w,     # перенос порядка на список В
             word,  # слово 
             ls;    # длина списка S
       
       Sor:=[1];
       Ac:=ShallowCopy(A);
       S:=ShallowCopy(Ac[1]);
       Remove(Ac,1);

    while Length(Ac)<>0 do
          ls:=Length(S);            # в списке А определяем позицию списка, в
       for i in Ac do               # котором содержится нужный нам последний
            if S[ls] in i then      # элемент из списка S.
               k:=Position(Ac,i);
            fi;
       od;
       if S[ls]=Ac[k][1] then       # создание списка S (список номеров вершин
           Add(S,Ac[k][2]);         # по порядку в направлении следования ориентации)
           Add(Sor,1);              # и списка Sor (ориентаций граней рёбер вершинами).
                        else
           Add(S,Ac[k][1]);
           Add(Sor, -1);
       fi;
       Remove(Ac,k);
    od;
       oo:=[];     # список ребер по порядку по направлению ориентации.
    for i in [1..ls] do 
           Add(oo,Set([S[i],S[i+1]]));  # создание списка оо
    od;
       Aoo:= PermListList(oo,A); # считаем перестановку из А в оо
       w:=Permuted(B,Aoo); # действуем найденной перестановкой на список В
   word:=(w[1]^Sor[1]);
   for i in [2..ls] do
       word:=word*(w[i]^Sor[i]);
   od;
   return word;
end );

################################################################################

InstallGlobalFunction( FundGroup,
function(p)
           local 1gr, 2gr, # 1-мерные и 2-мерные грани, соответственно
                 n, f, genf, m, A, B,
                 SG, i, md, SMD, fs,    # i - свободная переменная
                 FG, FGO, g,
                 P;


1gr:=p.faces[1];

  if Length(p.faces)<2 then
          2gr:=[];
        else 2gr:=p.faces[2];
  fi;

n:=Length(1gr);
f:=FreeGroup(n); 
genf:=GeneratorsOfGroup(f);
A:=[];
B:=[];
m:=Length(2gr);
SG:=[];
 
  
     for i in [1..m] do
              A[i]:=1gr{2gr[i]};  # рёбра, входящие в 2-грань i и заданные номерами вершин
              B[i]:=genf{2gr[i]}; # элементы группы на этих рёбрах
              SG[i]:=Slovo(A[i],B[i]);  # перевод граней в соотношения группы
     od;
 

  md:=MaxTree(p); 
  SMD:=genf{md};      # образующие входящие в максимальное дерево
  fs:=Union(SMD,SG);  # список отношений фундаментальной группы
  g:=f/fs;
  FG:=PresentationFpGroup(g);          # фундаментальная группа в общем виде
  TzInitGeneratorImages(FG);    # определяем старые генераторы
  TzGo(FG);      # упрощение фундаментальной группы
  
                 # imoldgen- обазы старых образующих выраженные через  новые
  P:=rec(FG:=FpGroupPresentation(FG), imoldgen:=TzImagesOldGens(FG));
return P;

 #    RelatorsOfFpGroup(_) команда вывода соотношений группы
 #    FpGroupPresentation(FG) команда перевода из представления группы в формат группы
  
end );

################################################################################

InstallGlobalFunction( PolMinusFace,
# ОПИСАНИЕ:
# функция которая из политопа p вырезает грань с адресом m
# входные данные: p - политоп, m - адрес вырезаемой грани
# выходные данные: новый политоп
# зависимости: PolBnd

function(p,m)
  local g,n,ind,i,st,d,s,pos,k,K,st_old,j,new,newk,gr;

g := StructuralCopy(p);
g.syms := [];

# 0) проверяем, имеет ли грань максимальную размерность
     # если да, до просто убираем её    
n:=Length(p.faces);
if n=m[1] then

       ind:=[1..Length(p.faces[n])];
       Remove( ind, m[2] );
       g.faces[n]:=g.faces[n]{ind};
 
else
 
# 1) вычислим в какие (m[1]+1)- грани входит грань m
  st:=[];
  s:=1;
  for i in g.faces[m[1]+1] do
    if m[2] in i then
       Add(st,s);
    fi;
    s:=s+1;
  od;

# если ни в одну, то убираем её из списка m[1]-граней 
  # и также уменьшаем на единицу номера следующих за ней m-граней в (m[1]+1)-гранях
  # из-за этого в функции PolMinusPol мы будем перебирать грани в обратном порядке
    #  - тогда их номера будут не испорчены

  if st=[] then
    if m[1]=0 then
      Remove(g.vertices, m[2]);
    else 
      Remove(g.faces[m[1]], m[2]);
    fi;
    for i in [1..Length(g.faces[m[1]+1])] do
      for j in [1..Length(g.faces[m[1]+1][i])] do
        if g.faces[m[1]+1][i][j] > m[2] then
          g.faces[m[1]+1][i][j] := g.faces[m[1]+1][i][j] - 1;
        fi;
      od;
    od;
  else

# 2) продублируем грань m.
    d:=Length(st);

# грань m входит в d различных (k+1)-мерных граней. 
# Создадим d-1 дубликат грани m, дубликатом под номером d объявим саму грань.
    new:=[m[2]];
    s:=1;

    if m[1]=0 then # если вырезаем вершину, добавляем еще (d-1) новых вершин, присваиваем им новые имена.
      K:=Length(g.vertices);
        for i in st{[1..d-1]} do
          Add(g.vertices,K+s); # добавление новой вершины 
          pos:=Position(g.faces[1][i],m[2]); # вычисление местоположения вершины m в i-м ребре
          Add(new,K+s); # индексы новых вершин.
          g.faces[1][i][pos]:=K+s; # замена старого индекса на новый.
          s:=s+1;
          Sort( g.faces[1][i] ); # задаем естественную ориентацию ребра
        od;

    else
      K:=Length(g.faces[m[1]]);
      for i in st{[1..d-1]} do
        Add(g.faces[m[1]],StructuralCopy(g.faces[m[1]][m[2]])); # добавление копий грани m
        pos:=Position(g.faces[m[1]+1][i],m[2]); # вычисление местоположения грани m в i-ой (m[1]+1)-грани
        Add(new,K+s); # индексы новых граней.
        g.faces[m[1]+1][i][pos]:=K+s; # замена старого индекса на новый.
        Sort( g.faces[m[1]+1][i] ); # ordering
        s:=s+1;
      od;
    fi;

# !!

    st_old:=st;
    for j in [m[1]+2..n] do # индукция по всем большим размерностям
      st:=[];
      s:=1; # счетчик индексов
      for i in g.faces[j] do  # определяем в какие j-грани входят (j-1)-грани из st_old (посчитано на предыдущем шаге)
        if Intersection(i,st_old)<>[] then
           Add(st,s);
        fi;
        s:=s+1;
      od;
# теперь нужно посчитать какие (j-1) грани добавятся при текущем шаге
# и какие индексы надо им приписать, включая в состав j-граней
      K:=Length(g.faces[j-1]);
      s:=1;
      newk:=[];
      for i in st do 
        gr:=Intersection(PolBnd(g,[j,i])[j-1],new);  # находим, какие новые (j-1)-грани добавятся
        Add(newk,K+s); # учитываем индекс созданной грани.
        Add(g.faces[j-1],gr); # добавляем новую грань в список (j-1)-граней
        Add(g.faces[j][i],K+s); # добавляем индекс новой грани в список образующих текущей грани 
        s:=s+1;
      od;
      new:=newk;
      st_old:=st;
    od;

  fi;
fi;

Unbind(g.syms); # unbinds the symmetries, if there were any
Unbind(g.fd); # unbinds the face dictionary, if any and if it was called "fd"

return g; 
end );

################################################################################

InstallGlobalFunction( DelFace,
# ОПИСАНИЕ:
# программа которая удаляет из списка грань с адресом [k,m] (k-размерность грани, m-индекс грани)
# входные данные: p - политоп, address - адрес грани.
# выходные данные: r - политоп без грани.
# зависимости: нет.

function(p,address)
     local k,m,Nk,s,t,i,j,r ;

r:=StructuralCopy(p);
# 1) extract  k  and  m
k:=address[1];
m:=address[2];

# 2) удаляем k-грань c индексом m
if k=0 then
   Nk:=[1..Length(r.vertices)];
   Nk:=Difference(Nk,[m]);
   r.vertices:=r.vertices{Nk};
 else
   Nk:=[1..Length(r.faces[k])];
   Nk:=Difference(Nk,[m]);
   r.faces[k]:=r.faces[k]{Nk};
fi;

# 3) удаляем ссылки на удаленную грань, эти ссылки могут быть только в (k+1)-гранях.
if k < Length(r.faces) then
   s:=1;
   for i in [1..Length(r.faces[k+1])] do
      r.faces[k+1][i]:=Difference(r.faces[k+1][i],[m]); # удаляем ссылку на грань
      t:=1;
      for j in r.faces[k+1][i] do # пересчитываем индексы
          if  j<m then
              r.faces[k+1][s][t]:=j;
            else
              r.faces[k+1][s][t]:=j-1;
          fi;
          t:=t+1;
      od;
      s:=s+1;
   od;
fi;

return r;
end );

################################################################################

InstallGlobalFunction( PolSimplify1,
# merge k-faces, k = n, n-1, ..., 1
# and if we actually simplified p, we do this cycle in  k  again
# until the marker  m  below shows that we did what we could with this method

function(p)
     local g,n,m,ch,k,ta,a,t,v1,v2,va,s,i,j,l,sa,len;

n := Length(p.faces);
g := StructuralCopy(p);
m := 1; # marker

while m>0 do
   m:=0;

   for k in [n,n-1..1] do
#   a) по всем (k-1)-граням
      sa:=0;
      if k=1 
        then len:=Length(g.vertices);
        else len:=Length(g.faces[k-1]);
      fi;

      for ta in [1..len] do

#   б) считаем, в каких k-гранях встречается текущая (k-1)-грань
         ch:=[];
         s:=1;
         for t in g.faces[k] do
            if (ta-sa) in t then
               Add(ch,s);
            fi;
            s:=s+1;
         od;

#   в) если текущая (k-1)-грань встречается только в двух k-гранях
         if Length(ch)=2 then
          if g.faces[k][ch[1]]<>g.faces[k][ch[2]] then # 
            v1:=PolBnd(g,[k,ch[1]]);
            v2:=PolBnd(g,[k,ch[2]]);
            va:=PolBnd(g,[k-1,(ta-sa)]);

#   г) если исследуемые k-грани пересекаются только по текущей грани

            l:=List([1..k], x->[]);
            for j in [1..k] do
              l[j] := IntersectionSet(v1[j],v2[j]); 
            od; # сделали набор общих граничных (j-1)-клеток для [k,ch[1]] и [k,ch[2]]

            # теперь смотрим, что они в размерностях 0..(k-2) такие же, как у клетки [k-1,(ta-sa)],
            # а в размерности (k-1) это просто клетка [k-1,(ta-sa)]
            if l{[1..k-1]}=va and l[k]=[ta-sa] then
               UniteSet(g.faces[k][ch[1]],StructuralCopy(g.faces[k][ch[2]])); # объединили клетки в одну
               if k < n then
                 for i in [1..Length(g.faces[k+1])] do # в те (k+1)-грани, в которых содержится ch[2], добавляем ch[1]
                   if ch[2] in g.faces[k+1][i] then
                     Add(g.faces[k+1][i],ch[1]);     # т.е. добавили ссылку на объединенную грань.
                   fi;
                 od;
               fi;
               g:=DelFace(g,[k,ch[2]]); # удалим k-грань с индексом ch[2]
               g:=DelFace(g,[k-1,ta-sa]); # удалим (k-1)-грань с индексом (ta-sa)
               sa:=sa+1; # пересчитали смещение индекса ta
            fi;
          fi;
        fi;
      od; #  завершили весь цикл для (k-1)-граней 
      m:=m+sa; # учитываем количество смещений за цикл
   od; # завершение цикла по k
od;

return g;
end );

################################################################################

InstallGlobalFunction( PolFactorInvolution,

function( p, invol )

# p is the polytope with symmetries
# invol is such a list of some of its symmetries (repetitions possible)
# that it is known that the product of symmetries in s is an involution

# returns the factored polytope

local s,sym,sp,sp_,lv,i,j,k,n;

s := StructuralCopy(p);

sym := List( [1..Length(s.faces)+1 ], 
  i -> Product( [1..Length(invol)], j->s.syms[ invol[j] ][i] ) );

sp := []; # here will be stored vertex numbers to be deleted
lv := Permuted([1..Length(s.vertices)], sym[1]);
for i in [1..Length(s.vertices)] do
  if i<lv[i] then 
    Add(sp,lv[i]);
    for j in [1..Length(s.faces[1])] do
      if lv[i] in s.faces[1][j] then
        s.faces[1][j] := Difference(s.faces[1][j],[lv[i]]);
        AddSet(s.faces[1][j],i);
      fi;
    od;
  fi;
od;
sp_ := Difference([1..Length(s.vertices)], sp); # list of numbers of remaining vertices
s.vertices := s.vertices{sp_};
for i in [1..Length(s.faces[1])] do
  for j in [1..Length(s.faces[1][i])] do
    s.faces[1][i][j] := Position( sp_, s.faces[1][i][j] );
  od;
od;

for k in [1..Length(s.faces)-1] do
  sp := [];
  lv := Permuted([1..Length(s.faces[k])], sym[k+1]);
  for i in [1..Length(s.faces[k])] do
    if i<lv[i] then 
      Add(sp,lv[i]);
      for j in [1..Length(s.faces[k+1])] do
        if lv[i] in s.faces[k+1][j] then
          s.faces[k+1][j] := Difference(s.faces[k+1][j],[lv[i]]);
          AddSet(s.faces[k+1][j],i);
        fi;
      od;
    fi;
  od;
  sp_ := Difference([1..Length(s.faces[k])], sp); # list of numbers of remaining k-faces
  s.faces[k] := s.faces[k]{sp_};
  for i in [1..Length(s.faces[k+1])] do
    for j in [1..Length(s.faces[k+1][i])] do
      s.faces[k+1][i][j] := Position( sp_, s.faces[k+1][i][j] );
    od;
  od;
od;

n := Length( s.faces );
sp := [];
lv := Permuted([1..Length(s.faces[n])], sym[n+1]);
for i in [1..Length(s.faces[n])] do
  if i<lv[i] then 
    Add(sp,lv[i]);
  fi;
od;
sp_ := Difference([1..Length(s.faces[n])], sp);
s.faces[n] := s.faces[n]{sp_};

Unbind( s.syms );
Unbind( s.fd );

return(s);

end );

################################################################################

