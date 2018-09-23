###############################################################################

# ОПИСАНИЕ
# Стягивание k-клетки которая состоит из двух вершин, дву 1-клеток, двух
# 2-клеток и т.д., до (k-1)-клеток (назовем такие клетки минимальными). k-клетка стягивается на одну из своих (k-1)-клеток
# ЗАМЕЧАНИЕ: Клетка стягивается на одну из ее (k-1)-клеток. Проверку того, что
# стягивание возможно оставляем за пользователем.
# входные данные: pol1 - политоп
#                 adr=[размерность, индекс] - адрес минимальной клетки
# выходные данные: политоп
# зависимости:

InstallGlobalFunction(ContractMiniFace, function(pol1,adr)
     local pol,ind,s,i;

	 pol:=StructuralCopy(pol1);
	if Length(pol.faces[adr[1]][adr[2]])=2 then
##  		pol:=StructuralCopy(pol1);
		pol:=DelFace(pol,adr);
		ind:=pol1.faces[adr[1]][adr[2]]; # узнаем (k-1)-клетки которые нужно склеить
		ind:=Set(ind); # вообще это множество должно уже быть сортированным

		# Во все клетки размерности adr[1] добавляем клетку ind[2].
		s:=1;
		for i in pol.faces[adr[1]] do
		  	if ind[2] in i then 
				Add(i,ind[1]);
				pol.faces[adr[1]][s]:=Set(i);
   			fi; 
   		s:=s+1;
		od;

		pol:=DelFace(pol,[adr[1]-1, ind[2]]); # удаление лишней клетки ind[2]
	else
		Print("It isn't a minimal face.\n");
		s:=s.stop;
	fi;


return pol;
end);

###############################################################################

# ОПИСАНИЕ
# Дробим выбранную k-клетку (k-1)-клеткой. 
# ЗАМЕЧАНИЕ: Проверку того, что данную клетку можно подразбить выбранным образом
# оставляем пользователю.
# входные данные: pol - политоп
#                 adr - адрес k-клетки которую будем дробить
#                 nabor - индексы (k-2)-клеток на которые будет натянута дробящая (k-1)-клетка
#                         если дробится 1-клетка, то объект nabor это имя новой вершины
# выходные данные: политоп
# зависимости:

InstallGlobalFunction(DivideFace, function(pol1,adr,nabor)
    local interes,pos,ind,l,s,i,pol, v;

    pol:=StructuralCopy(pol1);

    if adr[1]>1 then 
    #-------------------------------------------------------------------------

        # 1) выбираем все (k-1)-клетки разбиваемой
        interes:=pol.faces[adr[1]][adr[2]];
        # 2) смотрим из каких (k-2)-клеток они состоят
        interes:=pol.faces[adr[1]-1]{interes};
        # 3) из этих наборов выкидываем (k-2)-клетки на которые будут натянута новая (k-1)-клетка
        interes:=List(interes, i->Difference(i, nabor));
        # 4) разделяем полученный список множеств на два списка связных наборов
        pos:=ConnectedSubset(interes);
        ind:=StructuralCopy(pol.faces[adr[1]][adr[2]]);
        pos:=ind{pos}; # выделели (k-1)-клетки которые будут образовывать одну k-клетку
        # 5) добавляем дробящую (k-1)-клетку
        Add(pol.faces[adr[1]-1],Set(nabor));
        l:=Length(pol.faces[adr[1]-1]); # количество (k-1)-клеток
        # 6) заменяем раздрабливаемую клетку на половинку
        Add(pos, l); 
        pol.faces[adr[1]][adr[2]]:=pos;
        # 7) создаем вторую часть дробленной клетки
        ind:=Difference(ind,pos);
        Add(ind, l);
        # 8) добавляем вторую полвинку дробленной клетки
        Add(pol.faces[adr[1]], ind);
        # 9) учтем в (k+1)-клетках что мы подразбили клетку
        l:=Length(pol.faces[adr[1]]);
        if adr[1]<Length(pol.faces) then
           s:=1;
           for i in pol.faces[adr[1]+1] do
              if adr[2] in i then
                 Add(pol.faces[adr[1]+1][s],l);
              fi;
              s:=s+1;
           od;
        fi;

    else
    #--------------------------------------------------------------------------------
    #                         если дробится ребро

       v:=StructuralCopy(pol.faces[1][adr[2]][2]);
       l:=Length(pol.vertices);
       Add(pol.vertices,nabor); # добавили имя новой вершины
       l:=l+1;
       # дробим ребро
       pol.faces[1][adr[2]][2]:=StructuralCopy(l);
       Add(pol.faces[1],[v,l]);
       # отражаем дробление ребра на 2-клетках
       if Length(pol.faces)>1 then
            l:=Length(pol.faces[1]);
            s:=1;
            for i in pol.faces[2] do
                if adr[2] in i then
                    Add(pol.faces[2][s], l);
                fi;
                s:=s+1;
            od;
        fi;
    fi;

return pol;
end);


###############################################################################

# ОПИСАНИЕ
#  Программа сортирует политоп так, что бы первыми в списках .faces[k] шли
#  клетки которые лежат на границе.
# ЗАМЕЧАНИЕ: Если \alpha это i-ая k-клетка на гарнице, то эта клетка cтавится на место i, а та клетка которая имела индекс i занимает освободившееся место.
# входные данные: pol1 - политоп
# выходные данные:
# зависимости:

 InstallGlobalFunction(FirstBoundary, function(pol1)
      local pol,n,bord,sostB,sost,i,j,ver,s,k;


pol:=StructuralCopy(pol1);
n:=Length(pol.faces);

bord:=PolBoundary(pol);
sostB:=List(bord, i-> PolBnd(pol, [n-1, i]));

sost:=[];
for i in [1 .. n-1] do
   sost[i]:=[];
   for j in sostB do
      Append(sost[i], j[i]);
   od;
   sost[i]:=Set(sost[i]);
od; # вычислили все клетки на границе

ver:=Remove(sost, 1); # в одельный список выделили набор вершин лежащих на границе
s:=1;
for i in ver do
   if s = i then ;
   else
      pol:=PermFaces(pol,(s,i),0);
   fi;
   s:=s+1;
od; # теперь в списке pol.vertices сперва идут вершины которые лежат на границе

Add(sost, bord);

for k in [1 .. n-1] do
   s:=1;
   for i in sost[k] do
      if s = i then ;
      else
         pol:=PermFaces(pol,(s,i),k);
      fi;
      s:=s+1;
   od;
od;


return pol;
end);


###############################################################################

# ОПИСАНИЕ
#	Производим объединение двух политопов.
# ЗАМЕЧАНИЕ:
# входные данные:	pol1
#			pol2
# выходные данные:
#			объект типа политоп
# зависимости:

 InstallGlobalFunction(FreeUnionPol, function(pol1,pol2)
		local	pol,pol2ver,lp,dim2,i;


#Создаем новые вершины
pol:=StructuralCopy(pol1);
pol.vertices:=List(pol.vertices,x->[1,x]);
pol2ver:=List(pol2.vertices, x->[2,x]);
Append(pol.vertices,pol2ver);

lp:=Length(pol1.vertices); # length previously
dim2:=Length(pol2.faces);

for i in [1 .. dim2] do
	Append(pol.faces[i], pol2.faces[i]+lp);
	lp:=Length(pol1.faces[i]);
od;

return pol;
end);


###############################################################################

# ОПИСАНИЕ
#  Проверка является ли объект политопом
# ЗАМЕЧАНИЕ: Полная проверка на то, что все клетки являются дисками пока не
# представляется возможной. Тем не менее для размерностей от 1 до 3(4?) можно
# утверждать, что проверка является точной.
# входные данные: pol - политоп
# выходные данные: true\false
# зависимости:

#		  И Д Е Я
# 1) В записи обязательно должны присутствовать поля "vertices" и "faces".
# 2) Все клетки должны быть натянуты(ссылаться) на клетки которы были описаны ранее (кроме вершин).
# 3) Все клетки должны быть шарами\дисками. Вкачестве проверки этого факта взята характеристика Эйлера.
# 4) Не должно существовать клеток размерности меньшей dim, которые не входят ни в какие клетки максимальной размерности.


 InstallGlobalFunction(IsPolytope, function(pol)
     local	name_space,verify,kolichestvo,dim,n,ostov,
     		odnorodnost,invEuler,sostav,length,i;

# 1)---------------------------------------------------------------------------
# произовдится проверка наличия соответствующих обязательных полей

verify:=IsRecord(pol);
if verify then
	name_space:=RecNames(pol);
	verify:=("vertices" in name_space and "faces" in name_space);
else
	Print("Record havn't fild .vertices or .faces\n");
fi;

# 2)---------------------------------------------------------------------------

if verify then
	kolichestvo:=Length(pol.vertices);
	dim:=Length(pol.faces);
	n:=1;
	repeat
		ostov:=Set(Concatenation(pol.faces[n]));
		verify:=(Length(ostov) = kolichestvo);
		kolichestvo:=Length(pol.faces[n]);
		n:=n+1;
	until n>dim or (verify=false);
fi;
if verify then
else
	Print("В политопе имеются не корректные ссылки\n");
fi;

# 3)---------------------------------------------------------------------------

if verify then
	n:=1;
	repeat
		length:=Length(pol.faces[n]);
		
		i:=1;
		repeat
			sostav:=FaceComp(pol,[n,i]);
			invEuler:=List([0 .. n-1], x->(-1)^x *Length(sostav.(x)));
			invEuler:=Sum(invEuler);
			verify:=(invEuler = 1 - (-1)^n);
			i:=i+1;
		until 
			(i>length) or (verify=false);

		n:=n+1;
	until
		(n>dim) or (verify=false);
fi;

return verify;
end);

# ПРОВЕРКА:
# Проводилась на политопе S3cubic и на производных объекта полученных
# добавлением\удалением некоторых клеток. Так же попробовали провести проверку
# для Trefoil

###############################################################################

# ОПИСАНИЕ
#  Осуществляем заданную перестановку k-клеток политопа
# ЗАМЕЧАНИЕ: 
# входные данные: pol - политоп
#                 perm - перестановка 
#                 k - размерность клеток, которые перестанавливаем
# выходные данные:
# зависимости:

 InstallGlobalFunction(PermFaces, function(pol1, perm, k)
     local pol,lf,ind,lp,i;


pol:=StructuralCopy(pol1);

if k>0 then 
   pol.faces[k]:=Permuted(pol.faces[k], perm);
   lf:=Length(pol.faces[k]);
else
   pol.vertices:=Permuted(pol.vertices, perm);
   lf:=Length(pol.vertices);
fi;

if k<Length(pol.faces) then 
   # в списке ind на месте i  стоит индекс который получит i-ая клетка после перестановки
   ind:=ListPerm(perm);
   lp:=Length(ind);
   lf:=lf-lp;
   Append(ind,lp+[1 .. lf]);

   pol.faces[k+1]:=List(pol.faces[k+1], i->Set(ind{i}));
fi;


return pol;
end);


###############################################################################

# ОПИСАНИЕ:
# Упрощенный алгоритм поиска граничных (n-1)-клеток.
# ЗАМЕЧАНИЕ:
# входные данные:pol
# выходные данные:
# зависимости: НЕТ

 InstallGlobalFunction(PolBoundary, function(pol)
     local n,l,s,i,j,bondary,t;

n:=Length(pol.faces);
l:=Length(pol.faces[n-1]);
s:=[1..l]*0;

for i in pol.faces[n] do
   for j in i do
      s[j]:=s[j]+1;
   od;
od;

bondary:=[];
t:=1;
for i in s do
   if i=1 then 
      Add(bondary, t);
   fi;
   t:=t+1;
od;

return bondary;
end);


###############################################################################

# ОПИСАНИЕ:
# Пусть для некоторого симплекса \delta политопа pol имеется такая окрестность,
# что если удалить сам рассматриваемый симплекс из нее вместе с его гарницей, то
# мы получим несколько не пересекающихся дисков. То в таком случае данный
# симплекс \delta можно вырезать из политопа, более экономичнм способом.

# ЗАМЕЧАНИЕ: Проверку того, что такая окрестность для данного симплекса
# существует, оставляем за пользователем.
# входные данные: pol - политоп
#                 adr - адрес клетки которую удаляем
# выходные данные: политоп

 InstallGlobalFunction(PolMinusFaceDoublingMethod, function(pol1,adr)
         local pol,n,clas,dim,pos,lc,star,clasters,ldim,ind,i,name;

pol:=StructuralCopy(pol1);
n:=Length(pol.faces); # размерность политопа

dim:=adr[1];
pos:=adr[2];

if dim = n then 
	Remove(pol.faces[n],pos);
else
	# (dim+1)-клетки звезды объединяем в кластеры по n-клеткам звезды
	star:=StarFace(pol,adr);
	clasters:=List(star.(n), i -> FaceComp(pol,[n,i]).(dim+1));
	clasters:=List(clasters, i -> Intersection(i,star.(dim+1)));

	# список clasters должен распасться на не пересекающиеся списки, что и
	# соответствует распаду окрестности симплекса
	lc:=Length(clasters);
	clas:=ConnectedSubset(clasters); 	# нашли первый класс, в котором в
										# качестве дубликата будет сама клетка
	ind:=Difference([1..lc],clas);
	clasters:=clasters{ind}; 		# что осталось в кластерах

	lc:=Length(clasters);
	if dim=0 then 
		ldim:=Length(pol.vertices);
	else
		ldim:=Length(pol.faces[dim]);
	fi;
	while lc > 0 do
		ldim:=ldim+1;
		clas:=ConnectedSubset(clasters);
		ind:=Difference([1..lc],clas);
		clas:=Set(Concatenation(clasters{clas})); 
		for i in clas do	# замена исходного симплекса на дубликат
			pol.faces[dim+1][i]:=Difference(pol.faces[dim+1][i],[pos]);
			Add(pol.faces[dim+1][i],ldim);
		od;
		if dim=0 then	# нужно добавить новую вершину
			name:=StructuralCopy(ldim)-1;
			repeat
				name:=name+1;
			until 
				(name in pol.vertices)=false;

			Add(pol.vertices, name);
		else
			Add(pol.faces[dim],pol.faces[dim][pos]);
		fi;
		clasters:=clasters{ind};
		lc:=Length(ind);
	od;

fi;


return pol;
end);

# ПРОВЕРКА:
# на кластере из двух тетраэдров трехмерного движения Пахрена. Было произведено
# разделение кластера на два не пересекающихся тетраэдра.

###############################################################################

# ОПИСАНИЕ
# ЗАМЕЧАНИЕ:
# входные данные:	pol - политоп
#			setoffaces - набор индексов клеток для которых хотим найти гарницу
# 			dim -  размерность этих клеток
# выходные данные:	список индексов (dim-1)-клеток лежащих на гарнице исследуемых клеток
# зависимости:

 InstallGlobalFunction(SetOfFacesBoundary, function(pol,setoffaces,dim)
		local	nabor,setnabor,ind,t,i,s,j;


nabor:=pol.faces[dim]{setoffaces};
setnabor:=Set(Concatenation(nabor));
ind:=[];
t:=1;
for i in setnabor do
	s:=0;
	for j in nabor do
		if i in j then
			s:=s+1;
		fi;
	od;
	if s=1 then
		Add(ind,t);
	fi;
	t:=t+1;
od;


return setnabor{ind};
end);

# ПРОВЕРКА:
# проводилась с помощью программы ZeifertSurface, т.е. проверялось, что границей поверхности Зеферта действительно является заданный узел.

###############################################################################

# ОПИСАНИЕ
#  Находим все клети которые содержатся в рассматриваемой клетке adr.
# ЗАМЕЧАНИЕ: Саму клетку также включаем в свой состав
# входные данные:	pol - политоп
# 			adr - адрес клетки
# выходные данные: запись списков индексов входящих в состав клеток по размерностям
# зависимости: НЕТ

 InstallGlobalFunction(FaceComp, function(pol,adr)
	local	dim,pos,sostav,i;



dim:=adr[1];
pos:=adr[2];

sostav:=rec();
sostav.(dim):=[pos];

for i in [dim-1, dim-2 .. 0] do
	sostav.(i):=pol.faces[i+1]{sostav.(i+1)};
	sostav.(i):=Set(Concatenation(sostav.(i)));
od;

return sostav;
end);


###############################################################################

# ОПИСАНИЕ
# Функция создания триангулированной n-сферы и -шара (n-симплекса)
# ЗАМЕЧАНИЕ: 
# входные данные: n - размерность
# выходные данные: p - политоп
# зависимости: нет

# триангулированная n-сфера
 InstallGlobalFunction(sphereTriangul, function(n)
      local p,k,i,j,vertgr,ind ;


p:=rec(vertices:=[],faces:=[]);
p.vertices:=[1..n+2];
p.faces[1]:=Combinations(p.vertices,2);
vertgr:=[];
vertgr[1]:=p.faces[1];
for k in [2..n] do
   vertgr[k]:=Combinations(p.vertices,k+1);
   p.faces[k]:=[];
   for i in vertgr[k] do
      ind:=[];
      for j in Combinations(i,k) do
         Add(ind,Position(vertgr[k-1],j));
      od;
      Add(p.faces[k],ind);
   od;   
od;

return p; 
end);

# n-симплекс (триангулированный n-диск)
 InstallGlobalFunction(ballTriangul, function(n)
      local p;

if n=1 then
   p:=rec(vertices:=[1,2],faces:=[[[1,2]]]);
else
   p:=sphereTriangul(n-1);
   p.faces[n]:=[[1..n+1]];
fi;

return p;
end);

###############################################################################

# ОПИСАНИЕ:
# функция которая вычисляет звезду грани adr в политопе pol
# В качестве определения звезды берем следующее:
#	ЗВЕЗДОЙ для данной клетки назовем набор клеток большей размерности в которых содержится иследуемая нами.
# По таком определению сама клетка НЕ ВХОДИТ в свою звезду
# входные данные:	pol - политоп
#			adr - адрес грани, для которой строится звезда ( adr=[размерность,индекс])
# выходные данные: запись (record)  по размерностям клеток
# зависимости:

 InstallGlobalFunction(StarFace, function(pol,adr)
	local	dim,pos,star,n,s,ind,i,j;



dim:=adr[1];
pos:=adr[2];
star:=rec();
n:=Length(pol.faces); # размерность политопа

s:=1;
ind:=[];
for i in pol.faces[dim+1] do
	if pos in i then
		Add(ind,s);
	fi;
	s:=s+1;
od;
star.(dim+1):=StructuralCopy(ind);

for j in [dim+2 .. n] do
	s:=1;
	ind:=[];
	for i in pol.faces[j] do
		if IsEmpty(Intersection(i,star.(j-1))) then ;
		else
			Add(ind,s);
		fi;
		s:=s+1;
	od;
	star.(j):=StructuralCopy(ind);
od;


return star;
end);

# ПРОВЕРКА:
# на кластерах тетраэдров образующих трехмерные движения Пахнера.

###############################################################################

# ОПИСАНИЕ
#  Выделяем подполитоп из данного политопа.
# ЗАМЕЧАНИЕ:
# входные данные: pol - политоп
#                 ind - индексы клеток высших размерностей из которых состоит подполитоп (высших размерностей подполитопа)
#                 dim - размерность подполитопа
# выходные данные:
# зависимости:

 InstallGlobalFunction(SubPolytope, function(pol1, ind, dim)
     local Isost,pol,sost,k,i,ver,s,l;


# Основная идея заклучается в следующем: 
# Все клетки принадлежашиее клеткам ind перегоняем в начало списков. В таком случае первые (для каждой размерности разное количество) клетки в данной размерности будут образовывать интересующий нас подполитоп и нам не нужно будет возиться с определением того, какой индекс будет иметь та или иная клетка в новом политопе.

##  pol:=StructuralCopy(pol1);
pol:=rec();
pol.vertices:=StructuralCopy(pol1.vertices);
pol.faces:=StructuralCopy(pol1.faces);

# выделяем состав каждой отдельной клетки
Isost:=List(ind, i -> PolBnd(pol, [dim , i]));

# создаем общий список клеток размерности k лежащий на выделяемом подполитопе
sost:=[];
for k in [1 .. dim] do
   sost[k] := [];
   for i in Isost do
      Append(sost[k], i[k]);
   od;
   sost[k]:=Set(sost[k]);
od;

# список вершин выделим отдельно
ver:=Remove(sost,1);
# добавляем индексы самих клеток
#Add(sost,Set(ind));
# NOTE: В данном случае список ind должен быть сортированным

s:=1;
for i in ver do
	if i<>s then
		pol:=PermFaces(pol,(s,i),0);
	fi;
	s:=s+1;
od;

l:=Length(ver);
pol.vertices:=pol.vertices{[1 .. l]};

for k in [1 .. dim-1] do
	s:=1;
	for i in sost[k] do
		if i<>s then
			pol:=PermFaces(pol, (s,i),k);
		fi;
		s:=s+1;
	od;
	l:=Length(sost[k]);
	pol.faces[k]:=pol.faces[k]{[1 .. l]};
od;

pol.faces:=pol.faces{[1 .. dim]};
pol.faces[dim]:=pol.faces[dim]{ind};


return pol;
end);

###############################################################################

# ОПИСАНИЕ
# Вычисляем данные о движении Пахнера. Это задание симплекса до и после
# движения, представление этих данных в виде политопов, вычисление индексов в
# матрице fk которые соответствуют внутренним граням.
# ЗАМЕЧАНИЕ:
# входные данные: n - размерность пространства над которым вычисляется движение
#                 k - движение пахнера из k -> n+2-k (k симплексов --- левая часть, n+2-k --- правая) k=1,...,n+1
# выходные данные: информация о соответствующем движении Пахнера
#              .l - данные для левой части движения
#              .r - данные для правой части движения
#                 .pol - политоп
#                 .sim - симплекс
#                 .vnut - индексы внутренних граней
# зависимости: SimPol, 
#

 InstallGlobalFunction(dataPachner, function(n,k)
     local comb,date,e,pol,s,vnut,i,t,j,x, l,real,future,new_pos,old_pos;


# 1) создаем левый и правый симплексы движения пахнера k -> n+2-k
comb:=Combinations([1..n+2], n+1);
#sim:=rec(1:=comb{[1..k]}, 2:=comb{[k+1..n+2]}); # левый и правый симплексы

# подготовка данных
date:=rec(1:=rec(sim:=comb{[1..k]}),2:=rec(sim:=comb{[k+1..n+2]})); # создали левый и правый симплексы

# 2) Вычисляем внутренние грани (номера по которым идет интегрирование)
for e in [1,2] do
   # if e=1, then we calcul for l.h.s.
   # if e=2, then we calcul for r.h.s.
   pol:=FromSimplexToPolytope(date.(e).sim); # переводим данные из формата "симплекс" в формат "политоп"

   # 3) Вычисляем внутренние симплексы.
   s:=1;
   vnut:=[];
   for i in pol.faces[n-1] do
      t:=0;
      for j in pol.faces[n] do
         if s in j then t:=t+1; fi;
      od;
      if t=2 then Add(vnut,s); fi;
      s:=s+1;
   od; # вычислили индексы внутренних (n-1)-клеток в политопе pol
   # теперь эти индексы нужно перевести в индексы которым они соответсвуют столбцам в матрице (саму матрицу можно создать программой MatrixFkTriangulPol).

# Будет удобно если внутренние (n-1)-грани будут стоять на последних местах, а
# для одинаковых внешних (n-1)-граней слева и справа натянутых на одни и теже
# вершины индексы совпадали. Граням как внутренним так и внешним присваивается
# естественный порядок по упорядоваченному множеству вершин на которые грань
# натянута.
l:=Length(pol.faces[n-1]);
real:=List([1..l],x->PolBnd(pol,[n-1,x])[1]);
future:=Difference([1..l],vnut);
future:=StructuralCopy(real{future});
future:=Set(future);
Append(future,StructuralCopy(Set(real{vnut}))); # теперь список future содержит желаемый порядок

new_pos:=List(future,x->Position(real,x)); # старые индексы (n-1)-граней на новых местах
old_pos:=List(real,x->Position(future,x)); # новые индексы (n-1)-граней
pol.faces[n-1]:=pol.faces[n-1]{new_pos}; # пересортировали грани
pol.faces[n]:=List(pol.faces[n],x->Set(old_pos{x})); # преобразовываем индексы n-граней
vnut:=old_pos{vnut};

   # формируем выходные данные
   date.(e).pol:=StructuralCopy(pol);
   date.(e).vnut:=StructuralCopy(vnut);

od;

#redate:=rec(l:=date.1, r:=date.2);

return rec(l:=date.1, r:=date.2);
end);


###############################################################################

# ОПИСАНИЕ
# программа переводящая симплекс в политоп 
# ЗАМЕЧАНИЕ: Реализованный алгоритм возможно будет слишком затратным для больших триангуляций.
# входные данные: simp - симплекс
# выходные данные: p - политоп
# зависимости:

 InstallGlobalFunction(FromSimplexToPolytope, function(simp)
     local n,p,cs,vr,lv,s,j,i,x,tek;


# 1) подготовка данных
n:=Length(simp[1])-1; # размерность
p:=rec(vertices:=[], faces:=[]);
cs:=StructuralCopy(simp);

# 2) запуск основного цикла
while n>1 do
   p.faces[n]:=[];
    vr:=[]; lv:=0;
    s:=1;
   for j in cs do
      tek:=Combinations(j,n); # тут все комбинации будут упорядоченными
      p.faces[n][s]:=[];
      for i in tek do
         if i in vr then Add(p.faces[n][s],Position(vr,i)); # s:=s.stop;
         else Add(vr,i); lv:=lv+1; Add(p.faces[n][s],lv);
         fi;
      od;
      p.faces[n][s]:=Set(p.faces[n][s]);
      s:=s+1;
   od;
   cs:=vr;
   n:=n-1;
od;
p.vertices:=Union(simp);
# именами вершин могут не оказаться цифры, по этому нужно сделать подсчет индексов
p.faces[1]:=[];
for i in cs do
   tek:=List(i, x->Position(p.vertices,x));
   tek:=Set(tek);
   Add(p.faces[1],tek);
od;

return p;
end);

################################################################################
# version PL-2.3
################################################################################


# <ManSection><Func Name="GlueFaces" Arg="pol, faces, dim" />
# 	<Description>
# 		Склеить две клетки в многообразие у которых общая граница. Если склейку
# 		указанных дисков провести невозможно функция вернет fail. Переменная
# 		<C>faces</C> список из двух индексов клеток которые необходимо склеить,
# 		<C>dim</C> размерность этих клеток.
#       <!--
#       Мы говорим, что невозможно склеивать только в том случае, если нарушется
#       свойство комбинаторности политопа (когда существуют клетки косающиеся
#       сами себя). Если бы существовали клетки в которых бы нарушилась
#       комбинаторность (клетка стала косаться сама себя), то эти клетки имели
#       бы размерность `d+1` (мы склеиваем `d`-клетки). Из этого следует, что
#       существует такая `(d+k)`-клетка которая содержит обе склеиваемые клетки,
#       следовательно либо эта клетка содержит только эти две клетки в своей
#       границе, либо в ней нарушается комбинаторность, так как возникает
#       условие склеивания по границе.
#
#       Ну и конечно проверить, что скливаемые клетки имеют одинаковые границы
#       очень легко.
#       -->
#	</Description>
# </ManSection>

InstallGlobalFunction(GlueFaces, function(pol1, pos, dim)
	local	ab, N, ind, pol;

    N := Length(pol1.faces);
    pol := fail;
    # TODO: Тут на каждом этапе можно положить сообщение об ошибке, с рассказом
    # о том, почему такую склейку нельзя проделать.
    if pos[1] = pos[2] then
        # Получается, что мы пытаемся приклеить клетку к самой себе.
    elif dim > 0 and not Set(pol1.faces[dim][pos[1]]) = Set(pol1.faces[dim][pos[2]]) then
        # Проверили, что у клеток одна и таже граница.
    elif dim < N and Set(pos) in List(pol1.faces[dim+1], x -> Set(x)) then
        # Проверили, что несуществует клетки с границей из склеиваемых клеток.
    else
        pol := StructuralCopy(pol1);
        ab := Set(pos);
        if dim < N then
            Add(pol.faces[dim + 1], ab);
        else
            Add(pol.faces, [ab]);
        fi;
        ind := Length(pol.faces[dim + 1]);
        pol := ContractMiniFace(pol, [dim + 1, ind]);
        pol.faces := pol.faces{[1 .. N]};
    fi;

	# pol:=StructuralCopy(pol1);
    # if pos[1]=pos[2] then
    # else
    #     # идяе такая: добавляем "минимальную" клетку граница которой состоит из
    #     # face1 и face2, которые мы хотим склеить. Так как у этих клеток
    #     # одинаковая граница, то эта операция корректна и мы получаем клетку
    #     # размерности dim+1. После этого мы стягиваем минимальную клетку на
    #     # своию границу, такая функция нами уже была организована.
    #     ab:=Set([pos[1],pos[2]]);
    #     N := Length(pol.faces);
    #     if dim < N then
    #         Add(pol.faces[dim+1], ab);
    #     else
    #         Add(pol.faces, [ ab ]);
    #     fi;
    #     ind:=Length(pol.faces[dim+1]);
    #     pol:=ContractMiniFace(pol,[dim+1,ind]);
    #     pol.faces:=pol.faces{[1..N]};
    # fi;

return pol;
end);

################################################################################

# <ManSection><Func Name="ImageInPolProduct" Arg="pol1, pol2, d1xd2" />
# 	<Description>
# 	Фукнция указывает клеку которую образует декартово произведение двух
# 	клеток [adr1,pos1] и [adr2,pos2] из политопов pol1 и pol2,
# 	соответственно. На вход функции можно положить либо сами политопы которые
# 	участвовали в произведении, либо списки количества клеток каждой размерности
# 	для каждого из политопов.
#		<Example>
#gap> s1:=sphereAB(1);;
#gap> t2:=PolProduct(s1,s1);;
#gap> ImageInPolProduct(s1,s1,[[1,1],[1,2]]);
#gap> [2,2]
#gap> len:=rec(0:=Lenght(s1.vertices), 1:=Length(s1.faces[1]));
#gap> ImageInPolProduct(len,len,[[1,1],[1,2]]);
#gap> [2,2]
#		</Example>
#	</Description>
# </ManSection>

InstallGlobalFunction(ImageInPolProduct, function(pol1,pol2,d1xd2)
	local	d1, d2, k1, k2, dimpol1, dimpol2, l, i, adr, n;



	d1:=d1xd2[1][1];	# размерность первого диска
	d2:=d1xd2[2][1];	# размерность второго диска
	k1:=d1xd2[1][2];	# позиция первого диска в pol1
	k2:=d1xd2[2][2];	# позиция второго диска в pol2

	# Создаем возможность выбора данных.
	l:=rec(1:=rec(), 2:=rec());
	if ("vertices" in RecNames(pol1)) and ("faces" in RecNames(pol1)) then
		dimpol1:=Length(pol1.faces);
		l.1.0:=Length(pol1.vertices);
		for i in [1 .. d1+d2] do
			if i<= dimpol1 then
				l.1.(i):=Length(pol1.faces[i]);
			fi;
		od;
	else
		l.1:=pol1;
	fi;
	if ("vertices" in RecNames(pol2)) and ("faces" in RecNames(pol2)) then
		dimpol2:=Length(pol2.faces);
		l.2.0:=Length(pol2.vertices);
		for i in [1 .. d1+d2] do
			if i<= dimpol2 then
				l.2.(i):=Length(pol2.faces[i]);
			fi;
		od;
	else
		l.2:=pol2;
	fi;
	for i in [Length(RecNames(l.1)) .. d1+d2] do
		l.1.(i):=0;
	od;
	for i in [Length(RecNames(l.2)) .. d1+d2] do
		l.2.(i):=0;
	od;

	adr:=0;
	n:=d1+d2;
	for i in [0 .. d1-1] do
		adr:=adr + l.1.(i) * l.2.(n-i);
	od;
	adr:=adr + (k1 - 1)*l.2.(d2) + k2;

	
return [n, adr];
end);

################################################################################

# <ManSection><Func Name="VerticesRullGlueFace" Arg="pol,para,dim" />
# 	<Description>
# 		Производится склейка двух клеток политопа. Клетки индексы клеток
# 		помещаются в список para, размерность dim клеток указывается отдельно.
#		<Example>
#		</Example>
#	</Description>
# </ManSection>

InstallGlobalFunction(VerticesRullGlueFace, function(pol0, ab, dim)
	local	pol, sost_a, sost_b, namespace, set, pairs, kl, name, pos, kl2, p,
	d;

	pol:=StructuralCopy(pol0);
	sost_a:=FaceComp(pol,[dim,ab[1]]);
	sost_b:=FaceComp(pol,[dim,ab[2]]);
	namespace:=pol.vertices{sost_b.0};
	
	#--- нулевая итерация ------------------------------------------------------
	set:=Difference(sost_a.0,sost_b.0);
	pairs:=[];
	for kl in set do
		name:=pol.vertices[kl];
		pos:=Position(namespace,name);
		kl2:=Remove(sost_b.0,pos);
		Remove(namespace,pos);
		Add(pairs, Set([kl,kl2]));
	od;
	Sort(pairs, function(u,v) return u[2]>v[2]; end);
	for p in pairs do
		pol:=GlueFaces(pol,p,0);
	od;

	#--- запуск цикла ----------------------------------------------------------
	for d in [1 .. dim-1] do
		namespace:=pol.faces[d]{sost_b.(d)};
		set:=Difference(sost_a.(d),sost_b.(d));
		pairs:=[];
		for kl in set do
			name:=pol.faces[d][kl];
			pos:=Position(namespace,name);
			kl2:=Remove(sost_b.(d),pos);
			Remove(namespace,pos);
			Add(pairs, Set([kl,kl2]));
		od;
		Sort(pairs, function(u,v) return u[2]>v[2]; end);
		for p in pairs do
			pol:=GlueFaces(pol,p,d);
		od;
	od;

	#--- последняя итерация ----------------------------------------------------
	pol:=GlueFaces(pol,ab,dim);


return pol;
end);

################################################################################

# <ManSection><Func Name="VerticesRullGluePol" Arg="pol,subpol1,subpol2,dim" />
# 	<Description>
# 	В политопе указывается два набора связных подполитопов, которые необходимо
# 	склеить. Правила склеики подполитопов указываются в именах вершин (вершины с
# 	одинаковыми именами склеиваются в одну), далее все это индуцируется на
# 	клетки большей размерности. Для работы данного алгоритма необходимо, что бы
# 	у подполитопов были одинаковые pl-разбиения, а также, что бы в этом
# 	pl-разбиении содержалась хоть одна n-клетка натянутая хотябы на (n+1)
# 	вершин.
#		<Example>
#		</Example>
#	</Description>
# </ManSection>

# TODO: учесть смену индексации в списке sp2 при проведении склейки

InstallGlobalFunction(VerticesRullGluePol, function(pol1,subpol1,subpol2,dim)
	local	pol, sp1, sp2, s, l, sostav, verify, pairs, kl, name, namespace,
	pos, i2, i, para, normal;

	pol:=StructuralCopy(pol1);
	sp1:=Difference(subpol1,subpol2);
	sp2:=Difference(subpol2,subpol1);

	#--- нахождение комбинаторной клетки ---------------------------------------
    # Под комбинаторной клеткой мы подразумеваем такую клетку, которая может
    # быть однозначно идентифицирована в многообразии теми вершинами на которые
    # она натянута.
	s:=1;
	verify:=false;
    # Удаляем все клетки, которые имеют дубликаты.
	normal:=List(sp1, x -> FaceComp(pol, [dim,x]).0);
	normal:=List(normal, x -> Length(Positions(normal, x)));
	normal:=Positions(normal, 1);
	normal:=sp1{normal};
	l:=Length(normal);
	while s<=l and (not verify) do
        # Проверяем, что все подклетки также являются комбинаторными.
		sostav:=FaceComp(pol,[dim,normal[s]]);
		verify:=List([1..dim-1],x->pol.faces[x]{sostav.(x)});
		verify:=List([1..dim-1],x->Length(Set(verify[x]))=Length(verify[x]));
		verify:=Set(verify);
		verify:=(verify=[true]);
		s:=s+1;
	od;
	s:=s-1;
	if not verify then
		Print("Sorry, all faces is not combinatorial.\n");
		# break;
        #s := s.stop;    # TODO необходимо посмотреть как в GAP делается
                        # безболезненная остановка функций.
	fi;
	# На клетках sp1 необходимо содать упорядочение таким образом, что бы каждая
	# следующая клетка в подмногообразии имела общую грань хотябы с одной
	# клеткой которая уже была склеена.
	s:=Remove(sp1,s);
	Add(sp1,s,1);

	#--- запуск основного цикала -----------------------------------------------
    # Производим поиск пар клеток подполитопов которые должны быть склеены.
	pairs:=[];
	namespace:=List(sp2, x -> FaceComp(pol,[dim,x]).0);
	namespace:=List(namespace, x -> pol.vertices{x});
    namespace := List(namespace, x -> Set(x));      # добавлено 17 июля 2017
	for i in sp1 do
		kl:=[dim,i];
		name:=FaceComp(pol,kl).0;
		name:=pol.vertices{name};
        name := Set(name);                          # добавлено 17 июля 2017
		pos:=Position(namespace,name);
		i2:=Remove(sp2,pos);
		Add(pairs,Set([i,i2]));
		Remove(namespace,pos);
	od;
	Sort(pairs, function(u,v) return u[2]>v[2]; end);
	for para in pairs do
		pol:=VerticesRullGlueFace(pol,para,dim);
	od;

return pol;
end);

################################################################################

# <ManSection><Func Name="PreimageInPolProduct" Arg="pol1, pol2, imageface" />
# 	<Description>
# 		По заданному образу в декартовом произведении политопов pol1 и pol2
# 		указываем из каких клеток была составлена данная клетка imageface.
#		<Example>
#		</Example>
#	</Description>
# </ManSection>

InstallGlobalFunction(PreimageInPolProduct, function(pol1,pol2, imageface)
	local	dim, k, pol, l, verify, d, ind, b, a, i;

	dim:=imageface[1];
	l:=rec();
	#--- различные возможности -------------------------------------------------
	k:=1;
	for pol in [pol1,pol2] do
		if IsSubset(RecNames(pol),["vertices","faces"]) then
			l.(k):=rec();
			l.(k).0:=Length(pol.vertices);
			for i in [1 .. Length(pol.faces)] do
				l.(k).(i):=Length(pol.faces[i]);
			od;
		elif "0" in RecNames(pol) then
			l.(k):=pol;
		else
			Print("Data isn't correct in function PreimageInPolProduct\n");
			break;
		fi;
		for i in [Length(RecNames(l.(k))) .. dim] do
			l.(k).(i):=0;
		od;
		k:=k+1;
	od;

	#--- запуск цикла ----------------------------------------------------------
	verify:=true;
	d:=-1;
	ind:=StructuralCopy(imageface[2]);
	while verify do
		d:=d+1;
		ind:=ind - l.1.(d)*l.2.(dim-d);
		verify := (ind > 0);
	od;
	b:= ind + l.1.(d)*l.2.(dim-d);
	a:=1;
	verify:=true;
	while verify do
		b:= b - l.2.(dim-d);
		verify:= b > 0;
		a:=a+1;
	od;
	a:=a-1;
	b:=b + l.2.(dim-d);

return [[d,a],[dim-d,b]];
end);

################################################################################

# <ManSection><Func Name="EulerNumber" Arg="pol" />
# 	<Description>
# 		функция вычисляет число Эйлера для шарового комплекса. Данная функция
# 		полиморфна и способна принимать для вычислений данные типа политоп
# 		(IsPolytope), именнованные списки по размерностям в котором содержится
# 		структура политопа и именованный список по размерностям элементами
# 		которого могут выступать количества клеток определенной размерности. 
#		<Example>
# gap> EulerNumber(T2);
# 0
# gap> EulerNumber(sphereAB(4));
# 2
# gap> EulerNumber(sphereAB(3));
# 0
#		</Example>
#	</Description>
# </ManSection>

# зависимости:
# Read("~/ProgGAP/***.g");

InstallGlobalFunction(EulerNumber, function(data)
	local	namespace, ln, l, i, en;

	namespace:=RecNames(data);
	ln:=Length(namespace);
	if "vertices" in namespace then
		l:=rec(0:=Length(data.vertices));
		for i in [1 .. Length(data.faces)] do
			l.(i):=Length(data.faces[i]);
		od;
	elif "0" in namespace then
		l:=rec();
		if IsList(data.0) then
			l.0:=Length(data.0);
			for i in [1..ln] do
				if String(i) in namespace then
					l.(i):=Length(data.(i));
				fi;
			od;
		elif IsInt(data.0) then
			l.0:=data.0;
			for i in [1 .. ln] do
				if String(i) in namespace then
					l.(i):=data.(i);
				fi;
			od;
		fi;
	fi;

	en:=0;
	for i in RecNames(l) do
		en:=en + (-1)^(Int(i))*l.(i);
	od;

return en;
end);

################################################################################

# <ManSection><Func Name="UnionFaces" Arg="pol, kl1, kl2" />
# 	<Description>
# 		на вход функции посылаются две клетки одной и той же размерности для
# 		объединения их в одну клетку. По сути данной функцией реализуется
# 		обратная операция к функции <M>DivideFace</M> разбивающей клетку на две
# 		части.
# 		<Example>
# gap> UnionFaces(p3,[3,1],[3,2]);
# rec(
#   faces :=
#     [ [ [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], [ 1, 4 ], [ 2, 4 ], [ 3, 4 ], [ 1, 5 ],
#           [ 2, 5 ], [ 3, 5 ], [ 4, 5 ] ],
#       [ [ 2, 4, 6 ], [ 2, 7, 9 ], [ 4, 7, 10 ], [ 3, 5, 6 ], [ 3, 8, 9 ],
#           [ 5, 8, 10 ], [ 1, 4, 5 ], [ 1, 7, 8 ] ],
#       [ [ 1, 2, 4, 5, 7, 8 ], [ 3, 6, 7, 8 ] ] ],
#   vertices := [ 1, 2, 3, 4, 5 ] )
# 		</Example>
#
# 		Объединение двух клеток <M>D^k_1</M> и <M>D^k_2</M> в одну в политопе
# 		<M>pol</M> можно провести только в том случае если звезда пересечения
# 		этих клеток состоит только из <M>D^k_1</M> и <M>D^k_2</M>.
#
#		Для проведения объединения функция проверяет, что
# 		данные клетки пересекаются по одному диску. Если после объединения
# 		указанных клеток комплекс перестанет быть шаровым разбиением, то функция
# 		не будет проводить объединение, на выход будет подан начальный политоп и
# 		информационная строка с пояснением почему функция отказывается работать.
#		<Example>
# gap> UnionFaces(T2,[2,1],[2,3]);;
# This faces intersected on some balls or not intersected.
# I cannot union the faces in a polytope.
#		</Example>
#	</Description>
# </ManSection>


InstallGlobalFunction(UnionFaces, function(pol0, adr1, adr2)
	local	pol, dim, sost1, sost2, i, int, disk, generalsost, max, min, l,
	verify, ssilki;

	
	pol:=StructuralCopy(pol0);
	verify:=true;
	# первый критерий корректности совпадение размерностей шаров
	if adr1[1] = adr2[1] then
		dim:=adr1[1];
	else
##  		Print("Dimensions of both faces must be equal. \n");
		verify:=false;
	fi;
	sost1:=FaceComp(pol,adr1);
	sost2:=FaceComp(pol,adr2);
	int:=rec();
	for i in [0 .. dim-1] do
		# так как клетки не пересекаются по размерности dim, то мы ее исключим
		int.(i):=Intersection(sost1.(i), sost2.(i));
	od;

	# второй критерий корректности единственность в int.(dim-1)
	if Length(int.(dim-1)) = 1 then
		disk:=[dim-1, int.(dim-1)[1]];
	else
		# Пересечение указанных клеток не является диском
##  		Print("The intersection of the cells isn't a ball.\n");
		verify:=false;
		disk:=[];
	fi;

	# проверка на звезду
	if verify then
		ssilki:=Concatenation(pol.faces[dim]);
		ssilki:=Positions(ssilki,disk[2]);
		if Length(ssilki)=2 then
			generalsost:=FaceComp(pol,disk);
		else
			# звезда клетки diski имеет более чем две (dim+1)-клетки
##  			Print("The star the face ", disk, " have more then two (", dim,
##  			")-balls.\n");
			verify:=false;
			generalsost:=rec();
		fi;
	fi;

	# объединение двух клеток можно проводить если списки ind и generalsost
	# совпадают, это обусловлено тем, что если это так, то выбранные клетки
	# пересекаются только по одному однородному пространству, а для дисков это
	# может быть только диск размерность меньше.
	if verify and generalsost=int then
		pol:=DelFace(pol,disk);
		max:=Maximum(adr1[2],adr2[2]);
		min:=Minimum(adr1[2],adr2[2]);
		UniteSet(pol.faces[dim][min], pol.faces[dim][max]);
		# нужно во всех клетках размерности dim+1 перенаправить ссылки
		if Length(pol.faces) >= dim then
			l:=0;
		else
			l:=pol.faces[dim+1];
		fi;
		for i in [1 .. l] do
			if max in pol.faces[dim+1][i] then
				pol.faces[dim+1][i]:=Difference(pol.faces[dim+1][i],[max]);
				UniteSet(pol.faces[dim+1][i],[min]);
			fi;
		od;
		pol:=DelFace(pol,[dim,max]);
##  		Print("\t\t were unite \n");
	else
##  		Print("I cannot union the faces in a polytope. \n");
	fi;

return pol;
end);

################################################################################

# <ManSection><Func Name="wasDelFace" Arg="pol, adr" />
# 	<Description>
# 		функция корректирует сопуствующую информацию прикрепленную к политопу
# 		<M>pol</M> которая должна была измениться после удаления одной из
# 		клеток. На вход функции подается политоп и адрес той клетки которая
# 		была удалена. Напомним, что после удаления клетки индексы больших клеток
# 		данной размерности понижаются на едининцу, это изменение индексации
# 		должно быть отображено в той информации которая сопутствует данному
# 		шаровому комплексу.
#
#		Если из какой-либо сопуствующей информации была удалена клетка, то будет
#		выведено соответствующее сообщение, но изменения будут проведены.
#		
#		При обработке информации по 2-узле в политопе будет, в случае если
#		удаляется 2-клетка, будет выведено соответствующее сообщение, индекс
#		2-клетки будет удален из списка .2knot.sheets, если на данную 2-клетку
#		есть ссылка в .2knot.dpoints.(1kl), то соответсвующая позиция будет
#		очищена
#		<Example>
#		</Example>
#	</Description>
# </ManSection>

InstallGlobalFunction(wasDelFace, function(pol0,adr)
	local	namespace, pol, i, kl, names1kl, newdpoints, 1kl, r, list;

	namespace:=RecNames(pol0);
	pol:=StructuralCopy(pol0);
	
#--- корректрировка 2-узла -----------------------------------------------------
	if "2knot" in namespace then
		if adr[1] = 2 then
			i:=1;
			r:=0;
			for kl in pol.2knot.sheets do
				if kl > adr[2] then
					pol.2knot.sheets[i] := kl - 1;
				elif kl = adr[2] then
					Print("Attention.\t The 2-face of the 2-knot was delite.\n");
					r:=StructuralCopy(i);
				else
					pol.2knot.sheets[i] := kl;
				fi;
				i := i + 1;
			od;
			if r > 0 then
				Remove(pol.2knot.sheets,r);
			fi;

			for 1kl in RecNames(pol.2knot.dpoints) do
				i := 1;
				list:=[];
				for kl in pol.2knot.dpoints.(1kl) do
					if kl > adr[2] then
						list[i] := kl - 1;
					elif kl < adr[2] then
						list[i] := kl;
					fi;
					i := i + 1;
				od;
				pol.2knot.dpoints.(1kl):=list;
			od;
		fi;

		if adr[1] = 1 then
			names1kl:=List(RecNames(pol.2knot.dpoints), x->Int(x));
			i:=1;
			newdpoints:=rec();
			for 1kl in names1kl do
				if 1kl < adr[2] then
					newdpoints.(1kl):=pol.2knot.dpoints.(1kl);
				elif 1kl > adr[2] then
					newdpoints.(1kl-1):=pol.2knot.dpoints.(1kl);
				elif 1kl = adr[2] then
					Print("Attention.\t The 1-face of the 2-knot was delite.\n");
					Unbind(pol.(1kl));
				fi;
			od;
			pol.2knot.dpoints:=newdpoints;
		fi;
	fi;

return pol;
end);

################################################################################

# <ManSection><Func Name="PolMinusPol" Arg="pol, subpol, dim" />
# 	<Description>
# 		Функция вырезает подполитоп <M>subpol</M> из политопа <M>pol<M>. По
# 		возможности выбирается такой способ вырезания, который будет более
# 		экономичным.
#		<Example>
#		</Example>
#	</Description>
# </ManSection>

# зависимости:
# Read("~/ProgGAP/***.g");



InstallGlobalFunction(PolMinusPol, function(pol1,subpol,dim)
	local	pol, sostav, d, bigdim, kl;

	pol:=StructuralCopy(pol1);
	sostav:=rec();
	sostav.(dim):=subpol;
	for d in [dim, dim-1 .. 1] do
		sostav.(d-1):=Set(Concatenation(pol.faces[d]{sostav.(d)}));
	od;
	bigdim:=Length(pol.faces);
	if bigdim = dim+1 then
		for d in [dim, dim-1 .. 0] do
			for kl in sostav.(d) do
				pol:=PolMinusFaceDoublingMethod(pol,[d,kl]);
			od;
		od;
	else
		for d in [dim, dim-1 .. 0] do
			for kl in sostav.(d) do
				pol:=PolMinusFace(pol,[d,kl]);
			od;
		od;
	fi;

return pol;
end);

################################################################################

# <ManSection><Func Name="PolSimplify" Arg="pol" />
# 	<Description>
# 		проводит упрощение политопа <M>pol</M> фукцией UnionFaces. Данная
# 		функция перебирает все возможности начиная с клеток максимальной
# 		размерности. Вновь появившиеся возможности функция не исследует. По
# 		этому функцию можно запускать несколько раз если целью стоит максимально
# 		упростить политоп.
#		<Example>
# gap> a:=ballTriangul(3);;
# gap> PolSimplify(a);
# rec( faces := [ [ [ 1, 2 ], [ 1, 2 ] ], [ [ 1, 2 ], [ 1, 2 ] ], [ [ 1, 2 ] ]
#      ], vertices := [ 1, 2 ] )
#		</Example>
#	</Description>
# </ManSection>

# зависимости:
# Read("~/ProgGAP/***.g");



InstallGlobalFunction(PolSimplify, function(pol0)
	local	pol, dim, d, l, ssilki, chastots, wherecan, kl, pos;

	pol:=StructuralCopy(pol0);
	dim:=Length(pol.faces);

	for d in [dim-1, dim-2 .. 1] do
		l:=Length(pol.faces[d]);
		ssilki:=Concatenation(pol.faces[d+1]);
		chastots:=List([1 .. l], x -> Length(Positions(ssilki,x)));
		wherecan:=Positions(chastots, 2);
		Sort(wherecan, function(x,y) return x > y; end);
		for kl in wherecan do
			pos:=List([1 .. Length(pol.faces[d+1])], 
		 		x -> (kl in pol.faces[d+1][x]));
			pos:=Positions(pos,true);
			pol:=UnionFaces(pol,[d+1,pos[1]],[d+1,pos[2]]);
		od;
	od;

	l:=Length(pol.vertices);
	ssilki:=Concatenation(pol.faces[1]);
	chastots:=List([1 .. l], x -> Length(Positions(ssilki, x)));
	wherecan := Positions(chastots, 2);
	Sort(wherecan, function(x,y) return x > y; end);
	for kl in wherecan do
		pos:=List([1 .. Length(pol.faces[1])], x->(kl in pol.faces[1][x]));
		pos:=Positions(pos,true);
		pol:=UnionFaces(pol,[1,pos[1]],[1,pos[2]]);
	od;


return pol;
end);


################################################################################

# <ManSection><Func Name="ParallelSimplify" Arg="pol,subpol,dim" />
# 	<Description>
# 		В политопе <M>pol</M> производится упрощение с помощью функции
# 		UnionFaces, параллельно с этим упрощается подполитоп <M>subpol</M>
# 		размерности <M>dim.</M> Информация о клетках подполитопа помещается в
# 		список .subpol. Функция работает только для подполитопов чья размерность
# 		меньше размерности объемлющего пространства.
#
# 		Если вложенное подмногообразие <M>A</M> той же размерности, что и
# 		политоп <M>M</M>, можно провести параллельное упрощение c пересечением
# 		<M>C = A \cap \overline{(M\A)}.</M>
#		<Example>
#		</Example>
#	</Description>
# </ManSection>


InstallGlobalFunction(ParallelSimplify, function(pol0,subpol,dim)
	local	bigdim, sostav, d, l, may, kl, where, kl1, kl2, pol, paras,
	verify, elem, sost, i, newl, ind, zapret, wasdel;

	pol:=StructuralCopy(pol0);
	bigdim:=Length(pol.faces);

	d:=bigdim-1;
	while d > dim do
		l:=Length(pol.faces[d]);
		may:=Concatenation(pol.faces[d+1]);
		may:=List([1..l], x -> Length(Positions(may,x)));
		may:=Positions(may,2);
		while not IsEmpty(may) do
			kl:=Remove(may);
			where:=List(pol.faces[d+1], x -> kl in x);
			where:=Positions(where, true);
			kl1:=[d+1, where[1]];
			kl2:=[d+1, where[2]];
			pol:=UnionFaces(pol,kl1,kl2);
		od;
		d:=d-1;
	od;

	zapret:=StructuralCopy(subpol);
	l:=Length(pol.faces[d]);
	ind:=[1 .. l];
	may:=Concatenation(pol.faces[d+1]);
	may:=List(ind, x -> Length(Positions(may,x)));
	may:=Positions(may,2);
	may:=Difference(may, zapret);
	wasdel:=[];
	while not IsEmpty(may) do
		kl:=Remove(may);
		where:=List(pol.faces[d+1], x -> kl in x);
		where:=Positions(where, true);
		kl1:=[d+1, where[1]];
		kl2:=[d+1, where[2]];
		pol:=UnionFaces(pol,kl1,kl2);
		newl:=Length(pol.faces[d]);
		if newl < l then
			Add(wasdel,kl);
			l:=StructuralCopy(newl);
		fi;
	od;
	while not IsEmpty(wasdel) do
		i:=Remove(wasdel);
		Remove(ind,i);
	od;
	zapret:=List(zapret, x -> Position(ind,x));
	
	sost:=Set(Concatenation(pol.faces[d]{zapret}));
	d:=d-1;
	while d > -1 do
		if d = 0 then
			l:=Length(pol.vertices);
		else
			l:=Length(pol.faces[d]);
		fi;
		ind:=[1 .. l];
		may:=Concatenation(pol.faces[d+1]);
		may:=List(ind, x -> Length(Positions(may,x)));
		may:=Positions(may,2);
		while not IsEmpty(may) do
			kl:=Remove(may);
			where:=List(pol.faces[d+1], x -> kl in x);
			where:=Positions(where, true);
			kl1:=[d+1, where[1]];
			kl2:=[d+1, where[2]];
			verify:=Intersection(where, sost);
			verify:= not (Length(verify)=1);
			if verify then
				pol:=UnionFaces(pol,kl1,kl2);
				newl:=Length(pol.faces[d+1]);
				if l > newl then
					l:=StructuralCopy(newl);
					i:=1;
					for elem in sost do
						if elem > where[2] then
							sost[i]:=sost[i]-1;
						elif elem = where[2] then
							sost[i]:=where[1];
						fi;
						i:=i+1;
					od;
				fi;
			fi;
		od;
		if d > 0 then
			sost:=Set(Concatenation(pol.faces[d]{Set(sost)}));
		fi;
		d:=d-1;
	od;
	
	pol.subpol:=zapret;

	# Допустим клетка kl принадлежит подмногообразию, эту клетку можно удалять
	# из шарового комплекса и две клетки на размернсоть больше из звезды
	# Star(kl) не принадлежат подподмногообразию. Данный случай возможен только
	# водной единственной размерности, когда dim(kl) = dim. Во всех остальных
	# случаях в Star(kl) должны быть клетки из подмногообразия, следовательно
	# мощьность |Star(kl)| > 2 и данное удаление неразрешенно.

return pol;
end);


################################################################################

# <ManSection><Func Name="LengthPol" Arg="pol" />
# 	<Description>
# 		В <M>resul.d</M> указывается мощьность <M>d</M>-мерного остова в
# 		политопе <M>pol.</M>
#		<Example>
# gap> LengthPol(T2);
# total16
# rec( 0 := 4, 1 := 8, 2 := 4 )
#		</Example>
#	</Description>
# </ManSection>


InstallGlobalFunction(LengthPol, function(pol)
	local	l, dim, i; # total;
	
	l:=rec();
	l.0:=Length(pol.vertices);
	dim:=Length(pol.faces);
	for i in [1 .. dim] do
		l.(i):=Length(pol.faces[i]);
	od;
	
##  	total:=List([0 .. dim], x -> l.(x));
##  	Print("total: ", Sum(total),"\n");

return l;
end);

################################################################################

#			<ManSection><Func Name="ConnectedSum" Arg="N,M,e" />
#				<Description>
#					Фукнция создает связную сумму двух политопов <M>N</M> и
#					<M>M</M> одинаковой размерности. Параметр <M>e</M>
#					определяет будут ли склеены многообразия с сохранием
#					ориентаций или ориентация одного из многообразий (<M>M</M>)
#					будет изменена на противоположную.
#					<Example>
#gap> s2:=sphereAB(2);;
#gap> ConnectedSum(s2,s2);
#rec(
#  faces := [ [ [ 1, 2 ], [ 1, 2 ], [ 1, 2 ], [ 1, 2 ] ],
#      [ [ 1, 3 ], [ 1, 2 ], [ 3, 4 ], [ 2, 4 ] ] ],
#  vertices := [ [ 1, "A" ], [ 1, "B" ] ] )
#					</Example>
#					Если попробовать нарисовать этот пример, то получится д
#				</Description>
#			</ManSection>

InstallGlobalFunction(ConnectedSum, function(pol10, pol20, e)
    local pol1, pol2, n, a, b, bord1, bord1kl_a, bord2, bord1kl_b, sp1, sp2, l,
    kl, pol, vertname, ind, sost1, sost2, pairs, cellorient, orient, cell_id,
    pairs_orient, pos, i, j, cell;

	pol1:=StructuralCopy(pol10);
	pol2:=StructuralCopy(pol20);
	if Length(pol1.faces) = Length(pol2.faces) then
		n:=Length(pol1.faces);
	else
		Print("Dimensions of this polytopes must be equal.\n");
	fi;
	#
	# Еще более простой идеей является создать минимальную (n-1)-клетку в
	# каждом многообразии и уже тут и проводить вырезания. Думаю такая идея
	# является наиболее оптимальной и экономичной в данном случае. Плюс ко всему
	# это дает решение задачи о построении цилиндра над двумя сферами. Возмоно
	# даже самое оптимальное построение с точки зрения шаровых комплексов.
	
	# Нам не обязательно искать или строить минимальную клетку, поскольку для
	# нас имеет важность минимальность границы, то есть разбиения сферы.
	# Поступим так. Вырежем окрестность некоторого ребра в двух разбиениях.
	# Получившаяся сфера будет состоять из двуугольных гиппершаров. Лишние
	# клетки этого разбиения можно будет стягивать. Данное стягивание можно
	# продолжить до тех пор пока не останется 2 клеки в границе (это самое
	# минимальное разбиение сферы), Очевидно, что операция стягивания не будет
	# выводить из категории шаровых комплексов. 
	
	a:=Length(pol1.faces[1]);
	b:=Length(pol2.faces[1]);
	bord1:=PolBoundary(pol1);
	bord1kl_a:=List(bord1, x -> FaceComp(pol1,[n-1,x]).1);
	bord1kl_a:=Set(Concatenation(bord1kl_a));
	bord2:=PolBoundary(pol2);
	bord1kl_b:=List(bord2, x -> FaceComp(pol2,[n-1,x]).1);
	bord1kl_b:=Set(Concatenation(bord1kl_b));
	while a in bord1 do a:=a-1; od;
	while b in bord2 do b:=b-1; od;
	pol1:=PolMinusFace(pol1,[1,a]);
	pol2:=PolMinusFace(pol2,[1,b]);
	sp1:=Difference(PolBoundary(pol1), bord1);
	sp2:=Difference(PolBoundary(pol2), bord2);

	l:=Length(sp1);
	while l > 2 do
		kl:=Remove(sp1);
		pol1:=ContractMiniFace(pol1,[n-1, kl]);
		l:=l-1;
	od;
	l:=Length(sp2);
	while l > 2 do
		kl:=Remove(sp2);
		pol2:=ContractMiniFace(pol2,[n-1, kl]);
		l:=l-1;
	od;

	l:=Length(pol1.faces[n-1]);
	sp2:=sp2 + l;
	pol:=FreeUnionPol(pol1,pol2);
	vertname:=FaceComp(pol,[n-1,sp1[1]]).0;
	vertname:=pol.vertices{vertname};
	ind:=FaceComp(pol,[n-1,sp2[1]]).0;
	pol.vertices[ind[1]]:=vertname[1];
	pol.vertices[ind[2]]:=vertname[2];
	# по построению списки sp1 и sp2 упорядочены по возрастанию, так же будут
	# упорядочены и списки sostX.(i)
	sost1:=FaceComp(pol,[n-1,sp1[2]]);
	sost2:=FaceComp(pol,[n-1,sp2[2]]);
	sost1.(n-1):=sp1;
	sost2.(n-1):=sp2;

    # Оба многообразия должны содержать ориентации клеток высшей размерности,
    # если они ориентируемы. Ориентации должны находиться в поле orient. Если
    # соответствующего поля нет, или хоть в одно из многообразий orient=fail
    # (т.е. многообразие неориентируемо), то будет создана произвольная связная
    # сумма.
	for i in [0 .. n-2] do
		for j in [2,1] do
			pol:=GlueFaces(pol,[sost1.(i)[j],sost2.(i)[j]],i);
		od;
	od;
    # Осталось приклеить правльным образом друг к другу (n-1)-клетки.
    # Соответствие для склейки будет выбрано на основе ориентаций. Если хоть
    # одно из многообразий не ориентируемо или его ориентация не задана, то
    # будет реализована случайная связная сумма (без учета ориентации склейки).

    if IsList(pol1.orient) and IsList(pol2.orient) then
        pairs := StructuralCopy(sp1);
        Append(pairs, sp2);
        cellorient := CellOrient(pol1)[n];
        Append(cellorient, CellOrient(pol2)[n]);
        orient := StructuralCopy(pol1.orient);
        Append(orient, e * pol2.orient);
        cell_id := 0;
        pairs_orient := [];
        for cell in cellorient do
            cell_id := cell_id + 1;
            for j in [1 .. 4] do
                pos := Position(pol.faces[n][cell_id], pairs[j]);
                if not pos = fail then
                    pairs_orient[j] := orient[cell_id] * cell[pos];
                fi;
            od;
        od;
        Print(pairs_orient{[1,2]}); # Проверочная часть, эти пары должны быть
        Print(pairs_orient{[3,4]}); # равны либо [-1,1], либо [1,-1].

        if pairs_orient{[1,2]} = pairs_orient{[3,4]} then
            pol := GlueFaces(pol, sp1, sp2{[2,1]});
            pol.orient := orient;
        elif pairs_orient{[2,1]} = pairs_orient{[3,4]} then
            pol := GlueFaces(pol, sp1, sp2);
            pol.orient := orient;
        else
            Print("WARNING: Something wrong with induce orient on board.\n");
            pol := fail;
        fi;
    else
        pol := GlueFaces(pol, sp1, sp2);
    fi;

return pol;
end);

################################################################################
# version PL-2.8
################################################################################

#            <ManSection><Func Name="PolCanonicalOrder" Arg="pol" />
#                <Description>
#                    Функция создает каноническое упорядочение клеток внутри
#                    политопа. Каноническим упорядочением мы называем порядок
#                    внесенный вершинами клеток, то есть клетка d1 меньше клетки d2,
#                    если множество вершин d1 меньше множества вершин d2 в
#                    лексикографическом смысле.
#                </Description>
#            </ManSection>

InstallGlobalFunction(PolCanonicalOrder, function(pol0)
    local pol, n, l, newind_CanonicalOrder, vertcell, neword, newind, s, d, i;

    pol := StructuralCopy(pol0);
    n := Length(pol.faces);     # размерность многообразия
    l := LengthPol(pol);
    newind_CanonicalOrder := [];

    vertcell := List([1 .. l.0], x -> [x]);
    for d in [1 .. n-1] do
        neword := [1 .. l.(d)];
        vertcell := List(pol.faces[d], x -> Set(Concatenation(vertcell{x})));
        SortParallel(vertcell, neword);
        newind := [];
        s := 1;
        for i in neword do
            newind[i] := s;
            s := s + 1;
        od;
        pol.faces[d] := pol.faces[d]{neword};
        pol.faces[d+1] := List(pol.faces[d+1], x -> newind{x});
        Add(newind_CanonicalOrder, newind);
    od;

    neword := [1 .. l.(n)];
    vertcell := List(pol.faces[n], x -> vertcell{x});
    SortParallel(vertcell, neword);
    pol.faces[n] := pol.faces[n]{neword};
    newind := [];
    s := 1;
    for i in neword do
        newind[i] := s;
        s := s + 1;
    od;
    Add(newind_CanonicalOrder, newind);
    if Length(RecNames(pol)) > 2 then
        pol.newind_CanonicalOrder := newind_CanonicalOrder;
    fi;

    pol.faces := List( pol.faces, flist -> List( flist , y -> Set(y) ) );

    return pol;
end);

################################################################################

#            <ManSection><Func Name="GlueIndenticalSubpolitops" Arg="pol, sub1, sub2, dim" />
#                <Description>
#                    Склеить два подполитопа которые имею идентичные разбиеня на
#                    клетки.
#                </Description>
#            </ManSection>

InstallGlobalFunction(GlueIndenticalSubpolitops, function(pol0, sub1, sub2, dim)
    local pol, sost1, sost2, ind1, ind2, l, d, pair;

    pol := StructuralCopy(pol0);
    sost1 := List(sub1, x -> FaceComp(pol, [dim, x]));
    sost2 := List(sub2, x -> FaceComp(pol, [dim, x]));
    for d in [0 .. dim] do
        ind1 := Union(List(sost1, x -> x.(d)));
        ind2 := Union(List(sost2, x -> x.(d)));
        l := Length(ind1);
        while l > 0 do
            # WARNING считаем что старший индекс затирается при склейке
            pair := [ind1[l], ind2[l]];
            pol := GlueFaces(pol, Set(pair), d);
            l := l - 1;
        od;
    od;

    return pol;
end);

################################################################################
################################################################################
