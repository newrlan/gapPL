###############################################################################

#			<ManSection><Func Name="Pfaffian" Arg="mat" />
#				<Description>
#					Вычисляется пфаффиан кососимметрической матрицы.
#					<Example>
#gap> PrintArray(mat);
#[ [   0,  -3,  -1,   1 ],
#  [   3,   0,  -1,   2 ],
#  [   1,   1,   0,   3 ],
#  [  -1,  -2,  -3,   0 ] ]
#gap> Pfaffian(mat);
#-8
#					</Example>
#				</Description>
#			</ManSection>

InstallGlobalFunction(Pfaffian, function(mat1)
	  local lm,pf,s,a,perm,mat,j,i;

	mat:=MutableCopyMat(mat1);
	lm:=Length(mat);
	pf:=1;
	if lm mod 2 = 1 then 
		pf:=0;
		lm:=1;
	fi;

	while lm>4 do
		s:=1;
		repeat # в первой строке ищем не нулевой элемент
			s:=s+1;
			a:=mat[1][s];
		until not IsZero(a) or s = lm;
		if not IsZero(a) then 
			perm:=[1 .. lm];
			perm[2]:=s;
			perm[s]:=2;
			mat:=mat{perm}{perm};
			if s>2 then 
			 pf:=(-1)*pf;
			fi;
			for j in [1,2] do
				for i in [3 .. lm] do
					mat[i]:=mat[i] - mat[i][2]/mat[1][2] * mat[1];
				od;
				mat:=MutableCopyMat(TransposedMat(mat));
			od;
			pf:=pf*mat[1][2];
			mat:=mat{[3..lm]}{[3..lm]};
			lm:=lm-2;
		else 
			  pf:=0;
		fi;
	od;
	if lm=4 then 
		pf:=pf*(mat[1][2]*mat[3][4] -mat[1][3]*mat[2][4] +mat[1][4]*mat[2][3]);
	elif lm=2 then
		pf:=pf*mat[1][2];
	fi;

return pf;
end);
