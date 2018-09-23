
###############################################################################

# ОПИСАНИЕ
# "Выделение связной компоненты" если множество пересекается хотя бы с одни
# множеством из компоненты, то оно тоже принадлежит компоненте. Выводится первая
# же попавшаяся компонента.
# ЗАМЕЧАНИЕ:
# входные данные: list1 - список списков (подразумевается, что его элементами являются списки)
# 
# выходные данные: выводятся индексы элементов из списка list1.
#
# зависимости:
#



 InstallGlobalFunction( ConnectedSubset, function(list1)
     local ver, s,ind,del,i,list;


list:=StructuralCopy(list1);
ver:=list[1];
s:=1;
ind:=[1];
del:=[2,3..Length(list)];
while s<Length(ver)+1 do
   for i in del do
      if ver[s] in list[i] then
         Append(ver, list[i]);
         Add(ind,i);
      fi; 
   od;
   del:=Difference(del,ind);
   s:=s+1;
od;


return ind;
end );

# ПОЯСНЕНИЕ:
#

###############################################################################

InstallGlobalFunction( SortCircle, function(list)
	local	n,sort,s,t;

	n:=Length(list);
	sort:=[Remove(list,1)];
	s:=1;
	while s<n do
		t:=0;
		repeat t:=t+1;
		until not IsEmpty(Intersection(sort[s],list[t]));
		Add(sort,Remove(list,t));
		s:=s+1;
	od;

return sort;
end );

################################################################################
#						NOTE: По сути дела я искал линейный граф, описание
#						функции надо сделать через как раз через это

#			<ManSection><Func Name="LineOrdering" Arg="list" />
#				<Description>
#
#					Пусть имеется некоторое разбиение отрезка на интервалы,
#					которое записано в список <M>list.</M> Каждый интервал
#					представлен либо парой своих элементов - концов данного
#					интервала. На таком множестве можно ввести линейный порядок.
#					Функция возвращает список индексов, который задает линейных
#					порядок внутри списка <M>list</M> под именем <M>.order</M> и
#					ориентации каждого пары под именем
#					<M>.orient</M>					
#					<Example>
#gap> list:=[ [ 2, 3 ], [ 4, 1 ], [ 2, 4 ], [ 15, 1 ], [ 15, 5 ] ];;
#gap> LineOrdering(list);
#rec( order := [ 5, 4, 2, 3, 1 ], orient := [ 1, -1, 1, -1, 1 ] )
#					</Example>
#
#					Список list может представлять из себя цикл.
#					<Example>
#gap> circ:=[[1,2],[3,4],[4,2],[1,3]];;
#gap> LineOrdering(circ);
#rec( order := [ 3, 1, 4, 2 ], orient := [ 1, -1, 1, 1 ] )
#					</Example>
#				</Description>
#			</ManSection>

InstallGlobalFunction(LineOrdering, function(list0)
	local	list, n, order, ab, ind, s, i, index, verify, t, a, b, orient;

	list:=StructuralCopy(list0);
	n:=Length(list);
	index:=[1..n];
	order:=[ Remove(index) ];
	ab:=Remove(list);
	a:=Minimum(ab);
	b:=Maximum(ab);

	s:=StructuralCopy(n)-1;
	verify:= not IsEmpty(list);
	t:=0;
	while verify do
		ab:=Remove(list,1);
		i:=Remove(index,1);
		if a in ab then
			a:=Difference(ab,[a]);
			if not IsEmpty(a) then
				a:=a[1];
			fi;
			Add(order,i,1);
			t:=0;
		elif b in ab then
			b:=Difference(ab, [b]);
			if not IsEmpty(b) then
				b:=b[1];
			fi;
			Add(order,i);
			t:=0;
		else
			Add(list,ab);
			Add(index,i);
			t:=t+1;
		fi;

		verify:=not IsEmpty(list);
		verify:=verify and (t<s+1);
		if t=0 then s:=Length(list); fi;
	od;

	orient:=[1];
	list:=list0{order};
	for i in [2 .. Length(order)] do
		if list[i-1][2]=list[i][1] then
			Add(orient, orient[i-1]);
		else
			Add(orient, -orient[i-1]);
		fi;
	od;

return rec(order:=order, orient:=orient);
end);
