###############################################################################
# удаляем свободные петли в диаграмме узла, если есть

InstallGlobalFunction( Reidemeister10Everywhere, function(knot1)
	local knot, newname, name, l, ind, i, s ;

	knot:=StructuralCopy(knot1);
	newname:=List(knot.kod, i->i[1]);
	name:=List(knot.orient, i->[1]);
	l:=Length(knot.kod);
	ind:=[];
	for i in [1 .. l-1] do
		if newname[i]=newname[i+1] then
			Add(ind, i);
			Add(ind, i+1); # мы нашли петлю
			s:=Position(name, newname[i]);
			Remove(knot.orient,s); # удаляем информацию о петле
			Remove(name,s);
   		fi;
	od;
	Sort(ind, function(u,w) return u>w; end);
	for i in ind do # удаляем сами петли
   		Remove(knot.kod, i);
	od;

return knot;
end );

###############################################################################
# ОПИСАНИЕ
#  Создаем диаграму узла, удобную для построение паралели на дополнении узла,
#  имеющей нулевой коэффициент зацепления с данным узлом.
# ЗАМЕЧАНИЕ:
# входные данные: knot1 - диаграмма узла
# выходные данные:
# зависимости:

InstallGlobalFunction( ZeroLinkFromKnot, function(knot1)
	local knot,newname,link,abs,sign,l,s,i,a,b,c,name,j;
	
	knot:=Reidemeister10Everywhere(knot1);
    name := List(knot.orient, i -> i[1]);   # имена двойных точек (05.12.2017)
	# считаем линк зацепления по узлу
	link:=Sum(List(knot.orient, i->i[2]));
	abs:=AbsInt(link);
	sign:=SignInt(link);
	l:=Length(name)+1;
	s:=Length(knot.kod);

	for i in [1 .. abs] do
        # нужно создать имена новых вершин
        newname:=[];
        for j in [1 .. 3] do
            while l in name do
                l:=l+1;
            od;
            Add(newname, l);
            l:=l+1;
        od;

        a:=newname[1];
        b:=newname[2];
        c:=newname[3];

        # добавляем "зацепленную" петлю
        Add(knot.kod, [ StructuralCopy(c), -1],s);
        Add(knot.kod, [ StructuralCopy(b), -1],s);
        Add(knot.kod, [ StructuralCopy(a), -1],s);

        Add(knot.kod, [ StructuralCopy(c), 1],s);
        Add(knot.kod, [ StructuralCopy(b), 1],s);
        Add(knot.kod, [ StructuralCopy(a), 1],s);

        Add(knot.orient, [ StructuralCopy(a), -sign]);
        Add(knot.orient, [ StructuralCopy(b),  sign]);
        Add(knot.orient, [ StructuralCopy(c), -sign]);

        s:=s-1; 
        # Данный счетчик корректен так, как количество ребер в два раза превышает
        # количество двойных точек, а коэффициент зацепления не может быть больше
        # чем количество двойных точек.
    od;


return knot;
end );

###############################################################################

# ОПИСАНИЕ
# Возможно улучшение результата --- можно переписать указание поверхности
# Зейферта на тех 2-клетках которые уже есть, без добавления новой нижней.
# Может быть вообще получится упростить данный участок кода.

# Создаем поверхность Зейферта с заданной границей в виде узла K

# ЗАМЕЧАНИЕ:
# входные данные: K - диаграмма узла
# выходные данные: pol. - политоп
#                   .vertices
#                   .faces
#                   .knot - индексы 1-клеток по которым проходит узел
#                   .zeifert - индексы 2-клеток на которые натянута поверхность Зейферта
# зависимости:

 InstallGlobalFunction( ZeifertSurface, function(K)
     local d2,n,i,
           int,out,s,t,j, rebra,zcircles,okr,na4alo,konec,vot,l,
	   under,l2,ostov,pol,l3, underZ,tyt,ind,
	   grd2,lu2, underD,ost,gran,diski,nad,
	   1kl,surface,max,glubina,uzlovie,k, del, vertikal,
	   v,before,l0,vn,uzl,newname,name,ver;


# 1) построение граффа диаграммы

d2:=rec(vertices:=[], faces:=[[]]);
d2.vertices:=List(K.orient, i->i[1]);
n:=Length(K.kod); # это количество 1-клеток граффа
d2.faces[1]:=List([1..n-1], i->[K.kod[i][1], K.kod[i+1][1]]);
Add(d2.faces[1],[K.kod[n][1],K.kod[1][1]]); 
d2.faces[1]:=List(d2.faces[1],i->[Position(d2.vertices,i[1]), Position(d2.vertices, i[2])]);

# 2) определяем окружности Зейферта

#  а) для каждой точки определяем входящие и исходящие ребра
int:=[];
out:=[];
s:=1;
for i in [1 .. n/2] do
   int[s]:=[];
   out[s]:=[];
   t:=1;
   for j in d2.faces[1] do
      if i = j[1] then # то это исходящее ребро
         Add(out[s],t);
      elif i = j[2] then # то это входящее ребро
         Add(int[s],t);
      fi;
      t:=t+1;
   od;
   s:=s+1;
od;

# б) непосредственное выделение окружностей Зейферта
rebra:=[1 .. n]; # индексы всех ребер
zcircles:=[]; # окружности Зейферта
l:=StructuralCopy(n);
while l>0 do
   okr:=[rebra[1]];
   Remove(rebra,1);
   l:=l-1;
   s:=1;
   na4alo:=d2.faces[1][okr[1]][1];
   konec:=d2.faces[1][okr[s]][2];
   while na4alo <> konec do
      # ищем ребра исходящие из вершины konec, среди них выбираем то ребро чей индекс не равен okr[s]+1
      if okr[s] = n then # учитываем цикличность
         vot:=Difference(out[konec],[1])[1];
      else
         vot:=Difference(out[konec],[okr[s]+1])[1];
      fi;
      Add(okr, vot);
      s:=s+1;
      l:=l-1;
      konec:=d2.faces[1][okr[s]][2];
      rebra:=Difference(rebra,[vot]);
   od;
   # нашли замкнутую ориентируему окружность которая по построению является окружностью Зейферта
   Add(zcircles, StructuralCopy(okr));
od; #..........................................................................

# 3) вкладываем заданный узел в S^3
pol:=KnotInS3(K);
#    удаляем "заглушку" (3-диск который приклеивался, что бы создать полноценное S^3)
l3:=Length(pol.faces[3]);
Remove(pol.faces[3],l3);

# 4) Выделяем нижнее основание вложенного узла (under)

under:=[]; # список 2-клеток лежащих на нижнем основании (основании уровня "0")
l2:=Length(pol.faces[2]);
# узнаем на какие врешины натянута каждая 2-клетка
ostov:=List([1..l2], i->PolBnd(pol,[2,i])[1]);
# переводим индекс каждой вершины в информацию о том на каком основании она находится
ostov:=List(ostov, i->Set(pol.vertices{i}[2]));

s:=1;
for i in ostov do
   if i=["0"] then # все вершины этой клетки лежат на нижнем основании
      Add(under,s);
   fi;
   s:=s+1;
od;
# При таком выделении 2-клеток в списке under могут оказаться вертикальные
# 2-клетки, состоящие всего из двух вершин. Лишние 2-клетки из списка under
# можно удалить проверив лежат ли на ней ребра из узла
s:=1;
ind:=[];
for i in under do
   tyt:=Intersection(pol.faces[2][i], pol.knot);
   if IsEmpty(tyt) then
      Add(ind,s);
   fi;
   s:=s+1;
od;
under:=under{ind}; # выбрали клетки которые действительно лежат на нижнем основании

# 5) теперь нужно спроектировать окружности Зейферта на выделенное основание
# 1-клетки на основаниях имеют порядковую нумерацию соответствующую порядковой
# нумерации в списке pol.knot

ostov:=List(under,i->StructuralCopy(pol.faces[2][i]));
underZ:=Set(Concatenation(ostov));
# проекции окружностей Зейферта на нижнее основание
underZ:=List(zcircles, i->underZ{i});

# 6) выделяем диски которые ограничивают окружностями Зейферта

#  а) выделяем границу диска основания
grd2:=List(under, i->StructuralCopy(pol.faces[2][i])); 
lu2:=Length(under); # количество 2-клеток на нижнем основании
for i in [2 .. lu2] do
   tyt:=Difference(grd2[i],grd2[1]);
   grd2[1]:=Difference(grd2[1],grd2[i]);
   Append(grd2[1],tyt);
od;
grd2:=grd2[1]; # Список 1-клеток на границе диска основания
#ostov:=List(under,i->StructuralCopy(pol.faces[2][i]));

# б) Для каждой проекции окружности Зейферта находим диски, которые они
# ограничиывают. Идея заключается в следующем: последовательно выкидываем
# 2-клетки примыкающие к краю, пока не оголим границу исследуемой окружности,
# все что останется и будет нужным списком дисков.
underD:=[];
for okr in underZ do
   diski:=StructuralCopy(under); # текущий список 2-клеток основания 
   gran:=StructuralCopy(grd2); # текущий список 1-клеток на границе основания
   ost:=StructuralCopy(ostov); # текущие 1-клетки основания

   while (Set(okr)=Set(gran)) = false do
      s:=0;
      repeat  
         s:=s+1;
	 # смотрим исть ли у данной клетки ребро на крае диска
	 t:=Intersection(ost[s], gran);
	 # смотрим лежит ли это ребро на окружности (их может быть несколько)
         tyt:=Intersection(t,okr); #s:=s.top;
	 if IsEmpty(t) then
            tyt:=[0]; # учитываем, что эта 2-клетка внутренняя
	 fi; 
      until
         IsEmpty(tyt);

      tyt:=Remove(ost,s);
      Remove(diski,s);
      Append(gran,tyt);
      gran:=Difference(gran,t);
   od;

   Add(underD, StructuralCopy(diski));
od;
#------------------------------------------------------------------------------

# выделяем поверхность Зейферта

# выдлеяем вешины поверхности Зейферта над узлом
vertikal:=[];
surface:=[]; # список 2-клеток принадлежащих поверхности Зейферта
ind:=[];
1kl:=Set(Concatenation(underZ)); # набор 1-клеток на нижнем основании
s:=1;
for i in pol.faces[2] do
   v:=Intersection(i,pol.knot);
   if IsEmpty(v) then ;
   else
      if IsEmpty(Intersection(i,1kl)) then ;
      else
         Add(vertikal,s);
	 Add(ind,v[1]);
      fi;
   fi;
   s:=s+1;
od;
SortParallel(ind, vertikal);

# Сперва работаем с дисками которые лежат "выше" всех
s:=1;
ind:=[];
for tyt in underD do
   if Length(tyt)=1 then
      Append(surface, vertikal{zcircles[s]});
      Add(surface,tyt[1]); # список tyt это одноэлементный список
   else
      Add(ind,s);
   fi;
   s:=s+1;
od;

# работаем с оставшимися дисками
underD:=underD{ind};
underZ:=underZ{ind};
zcircles:=zcircles{ind};
# считаем уровни для оставшися дисков
nad:=[];
for i in underD do
   ind:=[];
   s:=1;
   for j in underD do
      if Intersection(i,j)=Set(i) then
         Add(ind,s);
      fi;
      s:=s+1;
   od;
   Add(nad,StructuralCopy(ind));
od;

glubina:=List(nad, i->Length(i));
max:=Maximum(glubina); # считам максимальную глубину вхождения дисков
# в списке ind[i] содержатся позиции дисков в которых содежится i-ый диск,
# имеются в виду позиции которые они имеют в списке underD
l2:=Length(pol.faces[2]);
l3:=Length(pol.faces[3]);
while max>0 do
   s:=1;
   for i in glubina do
      if i=max then
         l2:=l2+1;
	 Append(surface, vertikal{zcircles[s]});
	 Add(surface, StructuralCopy(l2));
	 # создаем "днище" которое приклеиваем к рассматриваемому диску (указания на диск которые сейчас создадим мы уже внесли в список surface)
         Add(pol.faces[2], StructuralCopy(underZ[s])); # добавили клетку по границе диска
	 Add(pol.faces[3], StructuralCopy(underD[s])); # создали новую 3-клетку
	 l3:=l3+1;
	 Add(pol.faces[3][l3], StructuralCopy(l2)); # указали 2-клетку закрывающую созданную 3-клетку

	 # теперь подкорректируем списки underD, укажем что некоторые клетки мы заменили новой
	 for j in nad[s] do
            underD[j]:=Difference(underD[j], underD[s]);
	    Add(underD[j],l2);
	 od;
      else 
         Add(ind,s);
      fi;
      s:=s+1;
   od;
   glubina:=glubina{ind};
   underD:=underD{ind};
   underZ:=underZ{ind};
   zcircles:=zcircles{ind};
   nad:=nad{ind};
   max:=max-1;
od;

# Заклеиваем заглушку
Add(pol.faces[3],PolBoundary(pol));




pol.zeifert:=surface;
return pol;
end );

###############################################################################

# ОПИСАНИЕ
# Построение простой границы многообразия Зейферта, которая состоит из двух отрезков (1-симплексов) 
# ЗАМЕЧАНИЕ:
# входные данные: knot - узел для которого стоим поверхность Зейферта
# выходные данные: pol - политом с указанной поверхностью Зейферта
# зависимости:

# ИДЕЯ: 
# 1) Выделяем порядок вершин и порядок ребер узла. 
# 2) Каждое очередное новое ребро подразбиваем, кадую очередную новую 2-клетку так же.
# 3) Делаем разрез по новым клеткам. Вклеиваем 3-диск с необходимым нам подразбиением.
# 4) Добавлем новую клетку к поверхности.



 InstallGlobalFunction( ZeifertSurfaceWithSimplyBoundary, function(knot)
     local 	pol,ver,i,vnut_1kl_zeif,kolvo_ver,previous_ver,l1,l2,new_ver,
     		star_z,order,pos,s,order_1kl,order_2kl,j,
		isol_ver,subpol,kolco,new_3kl,l3;



pol:=ZeifertSurface(knot);

# 1)---------------------------------------------------------------------------

ver:=[]; # порядок обхода вершин, соответствующий порядку обхода узла
for i in knot.kod do
	if i[2]=1 then
		s:=[i[1],"1"];
		Add(ver,Position(pol.vertices,s));
	else
		s:=[i[1],"0"];
		Add(ver,Position(pol.vertices,s));				# ++
	fi;
od;

vnut_1kl_zeif:=pol.faces[2]{pol.zeifert};
vnut_1kl_zeif:=Set(Concatenation(vnut_1kl_zeif));
vnut_1kl_zeif:=Difference(vnut_1kl_zeif,pol.knot);

# 2)---------------------------------------------------------------------------

kolvo_ver:=Length(ver);
previous_ver:=StructuralCopy(ver[1]);
l1:=Length(pol.faces[1]);
l2:=Length(pol.faces[2]);
new_ver:=StructuralCopy(kolvo_ver);

       #                                 # for checking
					#if IsPolytope(pol) then ;
					#else
						#Print("\n\n Oops! Line: 68. \n\n");
						#s:=s.stop;
					#fi;

for i in [2 .. kolvo_ver-1] do

		star_z:=StarFace(pol,[0,ver[i]]); # звезда вершины на поверхности Зейферта
		star_z.1:=Intersection(vnut_1kl_zeif,star_z.1);
		star_z.2:=Intersection(pol.zeifert,star_z.2);

	# Порядок встречи 1- и 2-клеток если идти рядом с границей поверхности.
	order:=StructuralCopy(pol.faces[2]{star_z.2});
	order:=List(order, x->Intersection(x, star_z.1));
	pos:=0;
	repeat pos:=pos+1; until pol.knot[i-1] in pol.faces[2][star_z.2[pos]];
	s:=Remove(star_z.2, pos);	# параллельное изменение списков pol.zeifert и order
	Add(star_z.2,s,1);
	s:=Remove(order,pos);
	Add(order,s,1);								# ++

	s:=ConnectedSubset(order);
	order:=order{s};
	order_1kl:=[order[1][1]];
	order_2kl:=star_z.2{s};
	for j in order do
		s:=Difference(j,order_1kl);
		Append(order_1kl,s);
	od;
	
	s:=1;
	for j in order_1kl do

			new_ver:=new_ver+1;
			pol:=DivideFace(pol,[1,j],new_ver);
			l1:=l1+1;
			Add(vnut_1kl_zeif,l1);

		pol:=DivideFace(pol,[2,order_2kl[s]],[previous_ver,new_ver]);
		previous_ver:=StructuralCopy(new_ver);
		s:=s+1;
		l2:=l2+1;
		l1:=l1+1;
		Add(pol.zeifert,l2);
		Add(vnut_1kl_zeif,l1);
	od;
	
					# for checking
       #                                 if IsPolytope(pol) then ;
					#else
						#Print("\n\n Oops! Line: 117. \n\n");
						#s:=s.stop;
					#fi;
od;

pos:=0;
repeat pos:=pos+1; until pol.knot[kolvo_ver-1] in pol.faces[2][pol.zeifert[pos]];
pol:=DivideFace(pol,[2,pol.zeifert[pos]],[ver[kolvo_ver],previous_ver]);
l2:=l2+1;
l1:=l1+1;
Add(pol.zeifert,l2);
Add(vnut_1kl_zeif,l1);

					# for checking
					#if IsPolytope(pol) then ;
					#else
						#Print("\n\n Oops! Line: 133. \n\n");
						#s:=s.stop;
					#fi;
# 3)---------------------------------------------------------------------------

# Сперва нужно разобраться какие 2-клетки добавились при обходе по краю. 
isol_ver:=ver{[2 .. kolvo_ver-1]};
subpol:=rec(vertices:=[], faces:=[]);
subpol.faces[2]:=[];
for i in pol.zeifert do
	s:=FaceComp(pol,[2,i]).0;
	s:=Intersection(s,isol_ver);
	if IsEmpty(s) then ;
	else
		Add(subpol.faces[2],i);
	fi;
od;
subpol.faces[1]:=Set(Concatenation(pol.faces[2]{subpol.faces[2]}));

kolco:=SetOfFacesBoundary(pol,subpol.faces[2],2);

new_3kl:=StructuralCopy(subpol.faces[2]);
for i in subpol.faces[2] do
	pol:=PolMinusFaceDoublingMethod(pol,[2,i]);
	l2:=l2+1;
	Add(new_3kl, l2);
od;

for i in Difference(subpol.faces[1],kolco) do
	pol:=PolMinusFaceDoublingMethod(pol,[1,i]);
od;

					# for checking
					#if IsPolytope(pol) then ;
					#else
						#Print("\n\n Oops! Line: 168. \n\n");
						#s:=s.stop;
					#fi;
# 4)---------------------------------------------------------------------------

Add(pol.faces[3],new_3kl);
l3:=Length(pol.faces[3]);

pol:=DivideFace(pol,[3,l3],kolco);
l2:=l2+1;

pol:=DivideFace(pol,[2,l2],pol.faces[1][pol.knot[kolvo_ver]]);
	
					# for checking
					#if IsPolytope(pol) then ;
					#else
						#Print("\n\n Oops! Line: 184. \n\n");
						#s:=s.stop;
					#fi;
#-------------------------------------------------(формирование выходных данных)

pol.knot:=[pol.knot[kolvo_ver], Length(pol.faces[1])];
pol.zeifert:=Difference(pol.zeifert,subpol.faces[2]);

	j:=pol.faces[2][l2];
#	second:=pol.faces[2][l2+1];
	s:=pol.faces[2]{pol.zeifert};
	if IsEmpty(Intersection(j,s)) then
		Add(pol.zeifert,l2+1);
	else
		Add(pol.zeifert,l2);
	fi;



for i in [1,2,3] do
	s:=1;
	for j in pol.faces[i] do
		pol.faces[i][s]:=Set(pol.faces[i][s]);
		s:=s+1;
	od;
od;

return pol;
end );

# ПРОВЕРКА:
# Проверка проводилась на узлах: трилистнике, восьмерке, неузле и узле 7_7. 
# В качестве проверочного фактора были взяты:
#	1) полученная запись должна быть политопом;
#	2) количество 2-клеток поверхности Зейферта до и после упрощения границы должно остаться неизменным;
#	3) граница поверхности Зейферта должна состоять из двух отрезков.

###############################################################################



# ОПИСАНИЕ
# Строим дополение узла в сфере S^3. Дополнительно указываем два метридиана, две
# параллели и коэффициент зацепления параллей (обе параллели будут иметь один и
# тот же коэффициент зацепления). Сводим границу полученного многообразия к
# четырем прямоугольникам (граница is тор).

#------------------------------------------------------------------------------#
#                        нарисуй рисунки к пояснениям алгоритма
#------------------------------------------------------------------------------#

# ЗАМЕЧАНИЕ: 
# входные данные: knot - диаграмма узла
#
# выходные данные:
#
# зависимости:
#

#


 InstallGlobalFunction( ComplementOfKnot, function(knot1)
     local knot, petly,link,sign,a,l,vn,newname,kosan,b,c,
           disk, ostov,s,vertical,i,pod,nad,pol,l3,l2,kriska,t,ind,j,sost,
           kasanie,l0,sost_nad,sost_pod,v,
           namek,drob,name,l1, antiv,eta,
           meridian,lp,ku4a,parallel, para,
           ver,bord,levo,prav;

# 1) погружаем узел в диск
disk:=KnotInS3(knot); # построили вложение узла в сферу S^3

# 2) Находим вертикальные 2-клетки над узлом и под узлом.
#  а) Для этого мы используем следующую идею: Так как все ребра узла проходят по
#  вертикальным 2-клеткам, и под каждой 1-клеткой на основании находится ровно
#  одно ребро, то каждая вертикальная 2-клетка лежит либо над либо под ребром
#  узла. Выделить вертикальные 2-клетки мы сможем по этому принцыпу (наличия в
#  ее составе ребра узла).
s:=1;
l2:=Length(disk.faces[2]);
ind:=[1..l2];
vertical:=[];
for i in disk.knot do
   t:=0;
   j:=1;
   vertical[s]:=[];
   while t<2 do
      if i in disk.faces[2][ind[j]] then
         Add(vertical[s], ind[j]);
         t:=t+1;
      fi;
      j:=j+1;
   od;
   s:=s+1;
od;

#  б) Все вертикальные 2-клетки могут быть либо треугольными, либо двуугольными
#  и четрырехугольными. Если посмотреть на вершины на которые натянуты 2-клетки
#  будет видно, что если клетка состоит из трех ветршин и она "под", то она
#  будет иметь две вершины отмеченные "0"  и одну вершину отмеченну "1" (и
#  наоборот). Если же клетка 4-угольная, то противоположная ей 2-клетка,
#  содержащая это же ребро узла будет двуугольником.
nad:=[];
pod:=[];
for i in vertical do
   sost:=List(i, x -> PolBnd(disk, [2,x])[1]);
   sost:=List(sost, x -> disk.vertices{x}[2]); # смотрим какими индексами помечены вершины (верхн,низ)
   s:=Length(sost[1]);
   if s=2 then
      if "0" in sost[1] then 
         Add(pod, i[1]);
         Add(nad, i[2]);
      else
         Add(nad, i[1]);
         Add(pod, i[2]);
      fi;
   elif s=3 then 
      Sort(sost[1]);
      if sost[1] = ["0", "0", "1"] then
         Add(pod, i[1]);
         Add(nad, i[2]);
      else 
         Add(nad, i[1]);
         Add(pod, i[2]);
      fi;
   elif s=4 then
      if "0" in sost[2] then 
         Add(pod, i[2]);
         Add(nad, i[1]);
      else
         Add(nad, i[2]);
         Add(pod, i[1]);
      fi;
   fi;
od;

# 3) Для дальнейших пострений мы предполагаем, что каждое ребро узла входит
# только в две 2-клетки и соответсвенно две 3-клетки. Что обеспечивается за счет
# способа построение погружения узла.

# Сейчас, будем вырезать все 1-клетки по которым проходит узел. Функция
# вырезания PolMinusFace устроена так, что при вырезании клетки (елси ее
# размерность не равна размерности многообарзия) она НЕ повреждает индексы
# остальных клеток которые существовали до вырезания, а новые клетки
# образованные при этой процедуре она добавляет в конец списков (замечание: для
# того, что бы осуществить сохранение индексов, одна из новых клеток вставляется
# на место вырезаемой).
# Параллельно, можно проверять сколько добавилось новых 2-клеток. По построению
# количество новых 2-клеток при каждом вырезании должно быть равно двум.
for i in disk.knot do
	pol:=PolMinusFace(disk, [1,i]);
od;

# 4) Для каждой вершины выделим 2-клетки которые ее касаются, если вершина
# принадлежит верхнему основанию, то 2-клетки будем выделять из списка
# вертикальних "под", и наоборот, если вершина принадлежит нижнему основанию, то
# 2-клетки будем выделять из списка вертикальных "над".
kasanie:=[];
l0:=Length(disk.vertices);
sost_nad:=List(nad, i -> PolBnd(pol,[2,i])[1]); # вершины на которые натянуты эти клетки
sost_pod:=List(pod, i -> PolBnd(pol,[2,i])[1]); # вершины на которые натянуты эти клетки
for v in [1..l0] do
   kasanie[v]:=[];
   if pol.vertices[v][2] = "1" then # индекс "1" означает верх, верхнее основание
      s:=0;
      t:=1;
      while s<2 do
         if v in sost_pod[t] then
            Add(kasanie[v],t);
            s:=s+1; # мы нашли одну 2-клетку которая касается данной вершины (всего их 2 по построению)
         fi;
         t:=t+1;
      od;
      # Сейчас индексы с списках kasanie это позиции в списках sost_pod. Нам нужны идексы 2-клеток
      kasanie[v]:=StructuralCopy(pod{kasanie[v]});
   else # в данном случае вершина принадлежит нижнему основанию
      s:=0;
      t:=1;
      while s<2 do
         if v in sost_nad[t] then
            Add(kasanie[v],t);
            s:=s+1;
         fi;
         t:=t+1;
      od;
      # Сейчас индексы с списках kasanie[v] это позиции в списках sost_nad. Нам нужны идексы 2-клеток 
      kasanie[v]:=StructuralCopy(nad{kasanie[v]});
   fi;
od; 



# 4) Теперь вырезаем вершины, но особым образом. (Лучше всего идея показана на
# рисунках). Идея заключается в том, что используется на алгоритм вырезания
# PolMinusFace, а для того что бы вырезать окрестность точки, сама эта тока
# дублируется, и некоторые ребра в качестве своей точки будут иметь
# дублированную. Так же добавляется еще пара ребер, которые будут принадлежать
# некоторым 2-клеткам, вместо вырезаемой точки. При таком вырезании список
# 3-клеток не затрагивается.
for v in [1..l0] do

   # a) выделим ребра которые имеют в себе исследуемую вершину
   s:=1;
   namek:=[];
   for i in pol.faces[1] do
      if v in i then
         Add(namek,s);
      fi;
      s:=s+1;
   od;
   
   # б) Находим индексы 1-клеток которые будут иметь дубликат данной вершины.
   # Все такие 1-клетки находятся на 2-клетках которые мы назвали касательными к
   # данной вершине
   drob:=List(kasanie[v], i-> StructuralCopy(pol.faces[2][i])); # 
   drob:=List(drob, i-> Intersection(i, namek));
   drob:=Concatenation(drob);
   drob:=Set(drob); # в этом списке будут содержаться повторяющиеся индексы
   
   # г) придумываем имя новой вершине.
   l0:=Length(pol.vertices); 
   name:=l0+1; # по построению такого имени вершины еще не было
   
   # д) замещение вершины на дубликат
   Add(pol.vertices,name);
   for i in drob do # по построению список drob должен быть длины 3
      pol.faces[1][i] := Difference(pol.faces[1][i],[v]);
      Add(pol.faces[1][i], name); # name больше всех остальных индексов
   od;
   
   # e) создаем две 1-клетки натянутые на исследуемую вершину и ее дубль
   l1:=Length(pol.faces[1]);
   Append(pol.faces[1], [ [v,name],[v,name] ]);
   
   # ё) Нам нужно узнать какие 2-клетки теперь будут иметь в себе новые
   # созданные нами 1-клетки.
   # Можно представить, что два ребра узла проходящие через вершину которую мы
   # пытаемся удалить (пусть и несколько хитрым способом) лежат на некоторой
   # вертикальной плоскости. Тогда все 2-клетки лежащие с одной стороны от этой
   # плоскости будут содеражать одно и тоже добавленное ребро (1-клетку),
   # 2-клетки с другой строны от этой плоскости будут содерать вторую копию от
   # этого ребра.
   antiv:=[disk.vertices[v][1]]; # находим противопложную вершину для данной
   if disk.vertices[v][2]="1" then
      Add(antiv,"0");
   else
      Add(antiv,"1");
   fi; # сейчас мы только нашли ее имя
   antiv:= Position(disk.vertices, antiv); # нашли индекс противопложной вершины
   
   # находим 2-клетки с той и этой строны в которых надо добавить ребра
   l1:=l1+1; # теперь это индекс первого добавленного ребра
   for i in kasanie[antiv] do 
	   # 3-клетки с одной сторны должны содержать одну(какую-то) клетку которая
	   # указана в списке kasinie для этой вершины.
      # Т.о. мы разделили клетки с "той" и "этой" строны. Теперь нужно найти
	  # 2-клетки с соответсвующей стороны в которые нужно добавить это 1-клетку.
	  # Одну такую 2-клетку мы уже знаем это клетка с индексом i, но мы не будем
	  # делать этого сразу, а осуществим это в цикле, позже.
      # непосредственно находим 3-клетки с "этой" стороны
      eta:=[];
      s:=0;
      t:=1;
      while s<2 do
         if i in pol.faces[3][t] then
            Add(eta,t);
            s:=s+1;
         fi;
         t:=t+1;
      od;
      
      # одно и тоже ребро должно быть добавленно к тем 2-клеткам с "этой" строны
	  # которые имеют и вершину v и вершину name.
      eta:=List(eta, i -> pol.faces[3][i]);
      eta:=Set(Concatenation(eta)); # список индексов 2-клеток с "этой" стороны
      s:=0;
      t:=1;
      while s<3 do 
		  # Такое условие выбрано по тому, что 2-клеток в которые добавится
		  # новое ребро с "этой" стороны будет всего 3.
         ostov:=PolBnd(pol,[2,eta[t]])[1]; # указали вершины на которые натянута данная 2-клетка
         if (v in ostov) and (name in ostov) then
            Add(pol.faces[2][eta[t]], l1);
            s:=s+1;
         fi;
         t:=t+1;
      od;
      
      l1:=l1+1; # 
   od;

od; #++++++++++++++++++++++++++++++ (были проверены: 1-остов и граница, вручную) +++++++++++++++++++++++++++++++++++++++++
#s:=s.stop;

# 4) Упрощение PL-разбиения
# Заметим, что если стянуть все вертикальные 2-клетки которые состоят только из
# двух ребер на одно из них, то ребро узла (или уже ребро параллели) просто
# окажется на основании диска. Такое стягивание не повредит политоп (т.е. после
# осущетсвления этой операции, все клетки остануться шарами и само многообразие
# не изменится). Так же можно стянуть все минимальные 2-клетки (состоящии только
# из двух ребер) на одно из из собственных ребер.

l2:=Length(pol.faces[2]);
s:=1;
while s<l2+1 do
   if Length(pol.faces[2][s])=2 then
      pol:=ContractMiniFace(pol, [2,s]);
      l2:=l2-1;
	  # если клетка минимальна, то мы ее удалили из списка 2-клеток
	  # (соответсвенно, список уменьшился, s менять не нужно)
   else
      s:=s+1; # если клетка не минимальная то переходим к следующему индексу
   fi;
od;

# 5) Выделение параллелей и меридианов.
# Все 1-клетки лежащие на меридианах находятся в конце списка (их количество равно количеству вершин)
l0:=Length(pol.vertices);
meridian:=[];
l1:=Length(pol.faces[1]);
meridian:=l1 +1 -[1 ..  l0]; # выбрали последние l0 индексов из всех индексов 1-клеток (пока еще этот список мы не разделяли)
ku4a:=PolBoundary(pol);
ku4a:=List(ku4a, i-> StructuralCopy(pol.faces[2][i]));
ku4a:=Set(Concatenation(ku4a));
ku4a:=Difference(ku4a, meridian); # список 1-клеток которые будут образовывать параллели

# Теперь мы можем указать две паралели. 
ostov:=List(ku4a, i-> pol.faces[1][i]);
parallel:=[];
parallel[1]:=ConnectedSubset(ostov);
parallel[1]:=ku4a{parallel[1]};
parallel[2]:=Difference(ku4a, parallel[1]);

## Мы можем выделить два меридиана: в списке meridian они группируются по два (по построению). 
#meridian[1]:=StructuralCopy([meridian[1], meridian[2]]);
#meridian[2]:=StructuralCopy([meridian[3], meridian[4]]);
#meridian:=meridian{[1,2]};

# Выбирем меридиан удобным для нас способом. Выбираем меридиан на вершинах узла
# между которыми расположено только одно ребро выбранное для параллели.
v:=disk.faces[1][disk.knot[1]];
# меридианы будем брать на этих вершинах
sost:=List(meridian, i->pol.faces[1][i]);
a:=[]; b:=[];
s:=0; t:=1;
while s<4 do
   if v[1] in sost[t] then
      Add(a,meridian[t]);
      s:=s+1;
   elif v[2] in sost[t] then
      Add(b, meridian[t]);
      s:=s+1;
   fi;
   t:=t+1;
od;
meridian:=StructuralCopy([a,b]);
#++++++++++++++++++++++++(ручная проверка параллелей и меридиан)++++++++++++++++++++++++++++++++

pol:=rec(vertices:=pol.vertices, faces:=pol.faces, parallel:=parallel, meridian:=meridian);


#-----------------------------------------------------------------------------------------------#
#------------------- Упрощаем границу до четырех прямоугольников -------------------------------#

# А) Подготовка к вклеиванию нового участка трубки для тора.

# Определяем верхние и нижние вершины меридианов.
ver:=List(pol.meridian{[1,2]}[1], i->pol.faces[1][i]); 
# из каждого меридиана взяли по одному ребру и выяснили на какие вершины он
# натянут. Так как все списки pol.faces[1][i] отсортированны по возрастанию, то
# по построению первым в этом списке будет идти индекс вершины построенной при
# погружении узла в 3-сферу, второй индекс это индекс дубликата этой вершины.
for i in [1,2] do
   if pol.vertices[ver[i][1]][2]="0" then # нулем были помечены вершины на нижнем основании
      ver[i]:=ver[i]{[2,1]}; 
   fi;
od; # сделали так, чтоб в списках ver первыми шли индексы вершин сверху

# Когда мы указывали меридианы, мы сделали так, что этими меридианами ребра
# образующие параллели сверху\снизу разбились на два участка. Один из этих
# участков содержит только одно ребро, второй все остальные ребра. Обозначим
# отдельными символами индексы ребер лежащих на первых участках.
s:=Position(pol.faces[1], Set(ver{[1,2]}[1])); # нашли индекс ребра натянутого на верхние вершины меридианов (по построению он единственный)
t:=Position(pol.faces[1], Set(ver{[1,2]}[2])); # нашли индекс ребра натянутого на нижние вершины меридианов
if s in pol.parallel[1] then
   parallel:=pol.parallel;
else
   parallel:=pol.parallel{[2,1]};
fi; # первым в списке parallel будут идти индексы ребер верхней параллели

# выделяем списки индексов на втором участке
parallel[1]:=Difference(parallel[1],[s]);
parallel[2]:=Difference(parallel[2],[t]);

bord:=PolBoundary(pol);
sost:=List(bord, i -> pol.faces[2][i]{[3,4]}); # выбрали только те ребра которые лежат на меридианах (по построению они последние в списках)
ind:=ConnectedSubset(sost);
levo:=bord{ind};             # индексы 2-клеток границы слева от узла
prav:=Difference(bord, levo);# индексы 2-клеток границы справа от узла
# лево и право выбраны условно

# определим "левое" и "правое" ребро выделенного меридиана
sost:=sost{ind}; # список sost состоит только из индексов ребер по которым проходят меридианы (не обязательно выделенные нами)
sost:=(Concatenation(sost));
for i in [1,2] do
   if pol.meridian[i][1] in sost then
   else
      pol.meridian[i]:=pol.meridian[i]{[2,1]};
   fi;
od; # сделали так, что бы первыми в списках меридианов шли индексы ребер слева

# из списка levo выкидываем индекс 2-клетки которые лежат на первом участке
sost:=List(levo, i->pol.faces[2][i]);
i:=0; j:=1;
while i<1 do
   if s in sost[j] then
      Remove(levo, j);
      i:=1;
   fi;
   j:=j+1;
od;
# из списка prav выкидываем индекс 2-клетки которые лежат на первом участке
sost:=List(prav, i->pol.faces[2][i]);
i:=0; j:=1;
while i<1 do
   if s in sost[j] then
      Remove(prav, j);
      i:=1;
   fi;
   j:=j+1;
od;
# теперь списки levo и pravo списки 2-клеток на оставшемся участке тора

# Б) Вклеиваине нового участка трубки тора

# Добавляем новое ребро между верхними\нижними вершинами меридианов
l1:=Length(pol.faces[1]);
Add(pol.faces[1], Set(ver{[1,2]}[1])); # верхнее ребро, его индекс l1+1
Add(pol.faces[1], Set(ver{[1,2]}[2])); # ниждее ребро, его индекс l1+2

# Создаем 2-клетки заменяющие 2-клетки слева и справа
l2:=Length(pol.faces[2]);
for i in [1,2] do
   j:=pol.meridian{[1,2]}[i]; # 1-клетки меридиана слева\справа
   Append(j, l1+[1,2]); # добавляем к ним связующие 1-клетки для меридианов (которые мы добавили)
   Add(pol.faces[2], StructuralCopy(j));
od; # заменяющая левая 2-клетка будет иметь индекс l2+1, правая --- l2+2

# Создаем 2-клетки сверху и снизу, на добавленных 1-клетках
Add(parallel[1], l1+1); # сверху
Add(parallel[2], l1+2); # снизу

Add(pol.faces[2], parallel[1]); # ее индекс l2+3
Add(pol.faces[2], parallel[2]); # ее индекс l2+4

# Создаем 3-клетки
# левую:
   Append(levo, [l2+1, l2+3, l2+4]);
   Add(pol.faces[3], levo);
# правую:
   Append(prav, [l2+2, l2+3, l2+4]);
   Add(pol.faces[3], prav);
   
# В) указываем новые параллели
pol.parallel:=[ [s, l1+1], # верхняя параллель
                [t, l1+2] ]; # нижняя параллель



for i in [1..3] do
   pol.faces[i]:=List(pol.faces[i], x->Set(x));
od;

return pol;
end );

###############################################################################

# ОПИСАНИЕ
#  Построение триангуляции многообразия полученного при вырезании трубчатой
#  окрестности узла. В качестве дополнительной информации указываются две
#  параллели и два меридиана.

# ЗАМЕЧАНИЕ: В качестве параллели выделена линия на трубчатой окрестности узла,
# имеющая нулевой коэффициент зацепления с этим узлом.
# входные данные: knot - узел
# 
# выходные данные:
#
# зависимости:
#



 
 InstallGlobalFunction( TriangulateComplementOfKnot, function(knot,orient)
    local pol,bord,a,b,s,i,sost,ver,j,perm,newpor,verxparal,vermerid;

    # Вырезание погружение узла из S^3.
    pol:=ComplementOfKnot(knot); 

    # вычисляем границу
    bord:=PolBoundary(pol); 



    #         схема разбиения границы на данный момент 
    #  +---------------------------+---------------------------+
    # 1|                          3|                          1|
    #  |                           |                           |
    #  |            c              |             d             |      vnes
    #  |                           |                           |
    #  |                           |                           |
    #  |                           |                           |
    #  |---------------------------+---------------------------| <-- верхняя параллель (условно)
    # 2|                          4|                          2|
    #  |                           |                           |
    #  |            a              |             b             |      vnut
    #  |                           |                           | 
    #  |                           |                           |
    #  |                           |                           |
    #  +---------------------------+---------------------------+
    # 1                           3                           1     1,2,3,4 --- vertetices name


    # Наша задача сейчас триангулизировать гарницу многообразия. Ее PL-разбиение
    # (разбиение границы) представлено на схеме выше. Для нас нужна триангуляция
    # тора двумя не пресекающимимся друг с другом группами линий, которые на данной
    # схеме образуют диагонали прямоугольников. Если на прямоугольнике а задать
    # диагональ (14), то такая же диагональ будет на прямоугольнике d, на
    # прямоугольниках c и b будет диагональ (23). Но на прямоугольнике a можно
    # выбрать другую диагональ и уже на остальных прямоугольниках диагонали будут
    # выбраны в соответсвии с нашим выбором.

    # Не произвольно на диаграмме (в нашем способе задания) определено некоторое
    # напраление по узлу. В соответсвии с этим можно было сказать какая клетка лежит
    # слева от узла (по направлению) какая справа, чем мы и пользовались при
    # построении погружения и вырезания узла. В программе CopmlementOfKnot в списках
    # meridian[i] первыми шли индексы 1-клеток которые лежат слева от узла по
    # направлению обхода. Первая же параллель в списке parallel это параллель над
    # узлом (вторая соответсвенно под узлом).


    #------------------------------------------------------------------------------#
    # К чему бы можно было привязать выбор набора диагоналей на клетках? Пока ответ
    # на этот вопрос не ясен. Поэтому здесь по алгоритму выбирается произвольная
    # клетка на которой добавляесят диагональ (14) (если orient = +1) или диагональ
    # (23) (если orient = -1)
    #------------------------------------------------------------------------------#

    # Сперва переместим границу в начало списков pol.faces и находим образы
    # меридианов и параллелей в новой индексации. Новый индекс у граничного ребра i,
    # будет его порядковый номер среди граничных ребер.
    bord:=PolBoundary(pol);
    sost:=List(bord, i->pol.faces[2][i]);
    sost:=Set(Concatenation(sost));

    pol:=FirstBoundary(pol);

    for i in [1,2] do
       for j in [1,2] do
          pol.meridian[i][j]:=Position(sost, pol.meridian[i][j]);
          pol.parallel[i][j]:=Position(sost, pol.parallel[i][j]);
       od;
    od;
    # в списке bord тоже изменились индексы 2-клеток на границе, т.к. у нас заведомо всего четыре прямоугольника, то
    bord:=[1..4];

    # Условно, вершинам на одно из параллелей (пусть для определенности это будет нижняя параллель) присвоим имена 1 и 3 (произвольно) и на другой параллели вершинам присвоим именя 2 и 4.
    verxparal:=pol.faces[1][pol.parallel[1][1]]; # смотрим на какие вершины натянута верхняя параллель
    vermerid:=List(pol.meridian, i-> pol.faces[1][i[1]]); # смотрим на какие вершины натянуты меридианы

    newpor:=[];
    newpor[2]:=verxparal[1];
    newpor[4]:=verxparal[2];

       if newpor[2] in vermerid[1] then 
          newpor[1]:=Difference(vermerid[1], newpor)[1];
          newpor[3]:=Difference(vermerid[2], newpor)[1];
       else
          newpor[3]:=Difference(vermerid[1], newpor)[1];
          newpor[1]:=Difference(vermerid[2], newpor)[1];
       fi;
    # newpor это порядок вершин который нам нужно задать.
    # Переставим вершины.
    perm:=PermListList(newpor,[1 .. 4]);
    pol:=PermFaces(pol,perm,0);
    # для удобства обращения переименуем первые 4 вершины их порядковыми номерами
    for i in [1 .. 4] do
       pol.vertices[i]:=i;
    od;

    # выбираем первую же попавшуюся 2-клетку и добавляем на ней диагональ, в пару к ищем еще одну 2-клетку которая не имеет общих ребер с первой клеткой
    a:=[bord[1]]; # первая группа прямоугольников на которых добавляется диагональ на одинаковые вершины
    b:=[];        # вторая группа прямоугольников на которых добавляется диагональ на одинаковые вершины
    for i in bord{[2..4]} do
       if IsEmpty(Intersection(pol.faces[2][a[1]], pol.faces[2][i])) then
          Add(a,i);
       else
          Add(b,i);
       fi;
    od; 

    # непосредтсвенное добавление диагоналей
    if orient = +1 then
       pol:=Diagonal2(pol,a[1],[1,4]);
       pol:=Diagonal2(pol,a[2],[1,4]);
       pol:=Diagonal2(pol,b[1],[2,3]);
       pol:=Diagonal2(pol,b[2],[2,3]);
    elif orient = -1 then
       pol:=Diagonal2(pol,b[1],[1,4]);
       pol:=Diagonal2(pol,b[2],[1,4]);
       pol:=Diagonal2(pol,a[1],[2,3]);
       pol:=Diagonal2(pol,a[2],[2,3]);
    else
       Print(" Function input mast be (knot, +-1). Please, chek youre the input. \n");
       # break;
       s := s.stop;     # TODO: Надо посомтреть как безболезненно прервать функцию в GAP
    fi; # триангулизовали границу
     
    # триангуляция остального политопа
    pol:=PolTriangulate(pol);



return pol;
end );

###############################################################################
#### PL-2.3
###############################################################################

# <ManSection><Func Name="SingularitySet2Knot" Arg="pol" />
# 	<Description>
# 		Singularity set --- граф двойных точек диаграммы 2-узла. В графе
# 		перечислены все 1-клетки объемлющего политопа <M>pol</M> которые
# 		содеражт описываемый граф. Так же в соответствии каждой вершине
# 		сопоставлен список. Данный список строится из звезды <M>Star(v)</M> в
# 		описываемом графе упорядоченной таким образом. Первая пара элементов
# 		принадлежит верхней линии двойных точек в этой вершине (получаемая на
# 		пересечении верхнего и среднего листов), вторая пара элементов
# 		принадлежит средней линии двойных точек (пересечение верхнего и нижнего
# 		листов) и третья пара принадлежит нижней линии двойных точек
# 		(пересечение среднего и нижнего листов). Таким образом соответсвующие
# 		списки для тройной точки будут состоять из шести элементов, для двойной
# 		точки из двух и для точки ветвления из одного.
#		<Example>
# gap> SingularitySet2Knot(TurnKnot(Trefoil,2));
# 	...
# rec( graf := [ 7, 4, 5, 6, 19, 20, 21, 22, 17, 18, 23, 8, 9, 10, 12, 45, 32,
#       31, 25, 26, 51, 54, 55, 64, 112, 113, 114, 115, 134, 136, 137, 138,
#       131, 130, 129, 119, 109, 118, 127, 106, 124, 105, 123, 104, 122, 102,
#       99, 98, 79, 97, 78, 53, 94, 52, 59, 93, 92, 57, 73, 56, 72, 90, 70, 69,
#       68, 67, 66, 65 ],
#   order := rec( 1 := [ 137, 92 ], 10 := [ 12, 102 ], 11 := [ 10, 12 ],
#       13 := [ 4 ], 14 := [ 22, 17, 19, 20, 104, 97 ],
#       15 := [ 19, 18, 21, 22, 105, 98 ], 16 := [ 20, 21, 18, 23, 106, 99 ],
#       19 := [ 25 ], 2 := [ 138, 94 ], 20 := [ 32, 31, 25, 26, 109, 102 ],
#       21 := [ 23, 26 ], 23 := [ 17, 32 ], 24 := [ 31, 112 ],
#       25 := [ 113, 104 ], 26 := [ 114, 105 ], 27 := [ 115, 106 ],
#       30 := [ 45, 109 ], 34 := [ 45, 112 ], 36 := [ 113, 118 ],
#       37 := [ 115, 119 ], 40 := [ 51, 56, 54, 53, 118, 122 ],
#       41 := [ 53, 52, 55, 56, 114, 123 ], 42 := [ 54, 55, 52, 57, 119, 124 ],
#       45 := [ 127, 59 ], 46 := [ 59, 57 ], 48 := [ 51 ],
#       49 := [ 64, 69, 67, 66, 129, 122 ], 5 := [ 4, 9, 7, 6, 97, 92 ],
#       50 := [ 66, 65, 69, 68, 130, 123 ], 51 := [ 68, 67, 70, 65, 131, 124 ],
#       54 := [ 72 ], 55 := [ 79, 78, 73, 72, 134, 127 ], 56 := [ 73, 70 ],
#       58 := [ 64, 79 ], 59 := [ 136, 78 ], 6 := [ 5, 6, 8, 9, 98, 93 ],
#       60 := [ 137, 129 ], 61 := [ 130, 93 ], 62 := [ 138, 131 ],
#       64 := [ 134, 90 ], 67 := [ 136, 90 ], 7 := [ 7, 8, 5, 10, 99, 94 ] ) )
#		</Example>
#	</Description>
# </ManSection>


InstallGlobalFunction(SingularitySet2Knot, function(pol)
	local	points, graf, chastots, ver, v , pos, 1kl, porjadok, k, 2kl, u, d,
	order;

	points:=TripleDoubleBranchPoints(pol);
	graf:=List(RecNames(pol.2knot.dpoints), x -> Int(x));
	chastots:=Concatenation(pol.faces[1]{graf});
	ver:=Set(chastots);
	order:=rec();
	for v in ver do
		pos:=Positions(chastots, v);
		pos:=pos - 1;
		pos:=pos - (pos mod 2);
		pos:=pos / 2;
		pos:=pos + 1;
		1kl:=graf{pos};
		porjadok:=[];
		if Length(1kl) = 6 then
			porjadok:=[];
			for k in 1kl do
				2kl:=pol.2knot.dpoints.(k);
				u:=Intersection(points.triple.(v).u, 2kl);
				d:=Intersection(points.triple.(v).d, 2kl);
				if IsEmpty(u) then
					Add(porjadok, 3);
				elif IsEmpty(d) then
					Add(porjadok, 1);
				else
					Add(porjadok, 2);
				fi;
			od;
			SortParallel(porjadok, 1kl);
			order.(v):=1kl;
		else
			order.(v):=1kl;
		fi;
	od;

return rec(order:=order, graf:=graf);
end);

################################################################################

# <ManSection><Func Name="PolSimplifyWith2Knot" Arg="pol" />
# 	<Description>
# 		упрощает политоп содержащий 2-узел. В качестве упрощающей функции была
# 		выбрана функция UnionFaces, которая объединяет два шара в политопе.
#
# 		Так же как и функция PolSimplify данная функция не проверяет возможны ли
# 		дальнейшие упрощение политопа <M>pol.</M> Для этой проверки необходимо
# 		еще раз запустить эту функцию на вновь полученных данных.
#		<Example>
#		
#		</Example>
#	</Description>
# </ManSection>

InstallGlobalFunction(PolSimplifyWith2Knot, function(pol0)
	local	pol, l2, ssilki, chastota, canwash, kl, where, pos, l3, l1,
	mayornot, para, i, p, l2new, name, 1klknot, l1new, l0, 1kl, 2kl, list, may,
	verify, 0kl;

	#--- объединение 3-клеток --------------------------------------------------
	pol:=StructuralCopy(pol0);
	l2:=Length(pol.faces[2]);
	canwash:=[1 .. l2];
	canwash:=Difference(canwash, pol.2knot.sheets);
	canwash:=Difference(canwash, PolBoundary(pol));
	Sort(canwash, function(x,y) return x > y; end);
	for kl in canwash do
		where:=List(pol.faces[3], x -> kl in x);
		pos:=Positions(where, true);
		2kl:=Intersection(pol.faces[3]{pos});
		l3:=Length(pol.faces[3]);
		pol:=UnionFaces(pol,[3,pos[1]], [3,pos[2]]);
		if l3 > Length(pol.faces[3]) then
			pol:=wasDelFace(pol, [2,2kl[1]] );
		fi;
	od;

	#--- объединение 2-клеток -------------------------------------------------

	# В политопе содержащем диаграмму 2-узла можно проводить объединение двух
	# 2-клеток в том случае, если это объединение можно провести в самом
	# шаровом комплексе и если обе клетки либо принадлежат, либо не принадлежат
	# диаграмме узла. В таком случае, если в политопе указана диаграмма
	# двумерной заузленной поверхности без края, то для всех пар клеток
	# подозрилтельных на объединение возможно только два случая, либо обе эти
	# клетки принадлежат диаграмме 2-узла, либо обе они не принадлежат этой
	# диаграмме.

	l1:=Length(pol.faces[1]);
	ssilki:=Concatenation(pol.faces[2]);
	chastota:=List([1 .. l1], x -> Positions(ssilki, x));
	chastota:=List(chastota, x -> Length(x));
	canwash:=Positions(chastota, 2);
	l2:=Length(pol.faces[2]);
	while not IsEmpty(canwash) do
		1kl:=Remove(canwash);
		pos:=List([1 .. l2], x -> 1kl in pol.faces[2][x]);
		pos:=Positions(pos,true);
		may:=Intersection(pos,pol.2knot.sheets);
		if Length(may) in [0,2] then
			pol:=UnionFaces(pol,[2,pos[1]],[2,pos[2]]);
			l2new:=Length(pol.faces[2]);
			if l2 > l2new then
				pol:=wasDelFace(pol,[1,1kl]);
				pol:=wasDelFace(pol,[2,pos[2]]);
				if pos[1] in pol.2knot.sheets then
					for name in RecNames(pol.2knot.dpoints) do
						list:=pol.2knot.dpoints.(name);
						for i in [1 .. 4] do
							if not IsBound(list[i]) then
								list[i]:=pos[1];
							fi;
						od;
					od;
				fi;
				l2:=l2new;
			fi;
		fi;
	od;

	#--- объединение 1-клеток --------------------------------------------------
	l0:=Length(pol.vertices);
	ssilki:=Concatenation(pol.faces[1]);
	chastota:=List([1 .. l0], x -> Positions(ssilki, x));
	chastota:=List(chastota, x -> Length(x));
	canwash:=Positions(chastota, 2);
	1klknot:=Set(Concatenation(pol.faces[2]{pol.2knot.sheets}));
	i:=1;
	l1 := Length(pol.faces[1]);
	while not IsEmpty(canwash) do
		0kl:=Remove(canwash);
		pos:=List(pol.faces[1], x -> (0kl in x) );
		pos:=Positions(pos, true);
		may:=Intersection(pos, 1klknot);
		if Length(may) in [0,2] then
			pol:=UnionFaces(pol,[1,pos[1]],[1,pos[2]]);
			l1new:=Length(pol.faces[1]);
			if l1 > l1new then
				pol:=wasDelFace(pol,[1,pos[2]]);
			fi;
		fi;
	od;

return pol;
end);

################################################################################

# <ManSection><Func Name="Knot1OnSphere2" Arg="knot" />
# 	<Description>
#		по диаграмме узла создается двумерная сфера <M>S^2,</M> в разбиении
#		которой указана данная диаграмма. Узел содержится в прикрепленном
#		именованном списке .1knot, который построен в стиле задания диаграмм
#		двумерных заузленных поверхностей. 
#		<Example>
#gap> Knot1OnSphere2(Figure8);
#rec(
#  1knot :=
#    rec(
#      dpoints := rec( 1 := [ 8, 1, 4, 3 ], 2 := [ 4, 5, 8, 7 ],
#          3 := [ 2, 3, 5, 6 ], 4 := [ 6, 7, 1, 2 ] ), sheets := [ 1 .. 8 ] ),
#  faces :=
#    [ [ [ 1, 4 ], [ 3, 4 ], [ 1, 3 ], [ 1, 2 ], [ 2, 3 ], [ 3, 4 ], [ 2, 4 ],
#          [ 1, 2 ] ],
#      [ [ 1, 3, 6 ], [ 1, 4, 7 ], [ 2, 5, 7 ], [ 2, 6 ], [ 3, 5, 8 ], [ 4, 8 ]
#         ] ], vertices := [ "a", "b", "c", "d" ] )
#		</Example>
#	</Description>
# </ManSection>

# зависимости: 



InstallGlobalFunction(Knot1OnSphere2, function(knot)
	local	graph, l, i, sections, rebra, ab,
		name, ver, down, rebro;

graph:=rec(vertices:=List(knot.orient, x -> x[1]), faces:=[[]]);

rebra:=List(knot.kod, x->x[1]);
Add(rebra, rebra[1]);
l:=Length(rebra);
while l > 1 do
	ab:=rebra{[1,2]};
	ab:=List(ab, x -> Position(graph.vertices, x));
	Add(graph.faces[1], Set(ab));
	Remove(rebra,1);
	l:=l-1;
od;

#...............................................................................

l:=Length(graph.faces[1]);
graph.faces[2]:=[];
for i in [1 .. l] do
	sections:=SectionLR(knot,i);
	Add(graph.faces[2], Set(sections.sectionL));
	Add(graph.faces[2], Set(sections.sectionR));
od;
graph.faces[2]:=Set(graph.faces[2]);
# После осуществления цикла в списке graph.faces[2] содержатся лишние клетки.
# Так как граф узла планарный, то это значит что в нем не содержится ни одной
# дублированный клетки (т.е. натянутой на те же самые ребра).

# Мы хотим указать как именно узел нужно располагать на сфере S^2. Для этого
# прикрепим дополнительные данные следующей структуры:
# 	.sheets 		- список ребер по которым проходит узел
# 	.dpoints.ver 	- список ребер для двойной точки ver, первые два
#			 	элемента этого списка соответствуют ребрам выше,
#			 	последние два --- ниже (это список строго из 4-х
#			 	элементов).
graph.1knot:=rec(sheets:=[1 .. l], dpoints:=rec());
for i in [1 .. l] do
	if knot.kod[i][2] = 1 then
		name:=knot.kod[i][1];
		ver:=Position(graph.vertices, name);
		if i=1 then rebro:=l; else rebro:=i-1; fi;
		down:=LeftRight(knot,rebro,name);
		graph.1knot.dpoints.(ver):=[rebro,i, down.left, down.right];
	fi;
od;





return graph;
end);

################################################################################

# <ManSection><Func Name="TurnKnot" Arg="knot, number" />
# 	<Description>
# 		создается диаграмма 2-узла вложенная в трехмерную сферу <M>S^3</M> с
# 		помощью алгоритма SpunTwist на основании диаграммы одномерного узла
# 		<M>knot</M>. Число <M>number<M> задает количество оборотов данной
# 		диаграммы при осуществлении twist-движения, при этом отрицательный знак
# 		данного числа задает обращение узла в противоположном направлении, при
# 		этом если <M>number</M> указать равным нулю, тогда twist-оборотов в
# 		диаграмме не будет и мы получим простую spun-диаграмму.
#		<Example>
#gap> TurnKnot(Trefoil,0);
#I'm trying simplify a polytope.
#
#  ...
#
# All good!
#
#rec(
#  2knot :=
#    rec( dpoints := rec( 11 := [ 12, 7, 9, 10 ], 12 := [ 8, 9, 11, 12 ],
#          13 := [ 10, 11, 13, 8 ], 14 := [ 18, 7, 15, 16 ],
#          15 := [ 14, 15, 17, 18 ], 16 := [ 16, 17, 13, 14 ] ),
#      sheets := [ 7, 8, 14, 9, 15, 10, 16, 11, 17, 12, 18, 13 ] ),
#  faces :=
#    [ [ [ 2, 3 ], [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], [ 1, 2 ], [ 5, 6 ], [ 4, 5 ],
#          [ 4, 6 ], [ 5, 6 ], [ 4, 5 ], [ 1, 4 ], [ 2, 5 ], [ 3, 6 ],
#          [ 1, 4 ], [ 2, 5 ], [ 3, 6 ] ],
#      [ [ 1, 3, 5 ], [ 1, 4 ], [ 2, 5 ], [ 6, 8, 10 ], [ 6, 9 ], [ 7, 10 ],
#          [ 11, 14 ], [ 1, 6, 12, 13 ], [ 2, 7, 11, 12 ], [ 3, 8, 11, 13 ],
#          [ 4, 9, 12, 13 ], [ 5, 10, 11, 12 ], [ 13, 16 ], [ 1, 6, 15, 16 ],
#          [ 2, 7, 14, 15 ], [ 3, 8, 14, 16 ], [ 4, 9, 15, 16 ],
#          [ 5, 10, 14, 15 ] ],
#      [ [ 7, 10, 13, 16 ], [ 1, 4, 8, 10, 12 ], [ 2, 5, 8, 11 ],
#          [ 3, 6, 9, 12 ], [ 1, 4, 14, 16, 18 ], [ 2, 5, 14, 17 ],
#          [ 3, 6, 15, 18 ] ] ],
#  vertices := [ [ 1, "a" ], [ 1, "b" ], [ 1, "c" ], [ 2, "a" ], [ 2, "b" ],
#      [ 2, "c" ] ] )
#		</Example>
#	</Description>
# </ManSection>

# TODO: после выбора осевой клетки и дробления в точках отвечающих за северный и
# южный полюса, т.к. они являются двойными точками узла произошла смена индексов
# клеток входящих в отношение выше\ниже. Эту информацию надо исправить.

InstallGlobalFunction(TurnKnot, function(knot,number)
	local	fk, bord, rebro, rebra_knot, name, verknot, i, verify, l, l3, l2,
	l1, l0, os, p, 2kl_os, vnut_2kl_knot, new1kl, north_r, south_r, north,
	south, point, cepi, 0kl, 1kl, 2kl, 3kl, kl, kl1, kl2, p2p3, p5p7, new2kl,
	sost_new2kl, k, true_ord, perm, rebra, sost_rebra, r, kadr, n, s1, ls1, lfk,
	pol, 2sheets, 2dpoints, level_sheets, order, nk, lsheets, orient, star, v,
	tor, granica, set_kl, dim, set, kolco1, kolco2, kolco, 1kl_tor, 1kl_2knot,
	lines, subsets, 1kl_kolco, 2kl_kolco, ind, d1, d2, plus, diski, para, osa,
	osb, a, b, list, bok, sost2kl, set1kl, chastota, j, image, vertical,
	1kl1,1kl2,0kl2, ord, pos, massa, newmassa, images, netak, s;

	#--- подготовка данных -----------------------------------------------------

	fk:=Knot1OnSphere2(knot);
	bord:=Remove(fk.faces[2],1);
	rebro:=Remove(bord,1);
	rebra_knot:=Difference([1 .. Length(fk.faces[1])], [rebro]);
	name:=Length(knot.orient);
	verknot:=List(knot.orient, x -> x[1]);
	for i in [1,2] do
		repeat
			name:=name+1;
			verify:=not (name in verknot);
		until verify;
		fk:=DivideFace(fk,[1,rebro],name);
	od;
	l1:=Length(fk.faces[1]);
	l0:=Length(fk.vertices);
	os:=[-1,0] + l0;
	if fk.faces[1][l1] = os then	# создаем осевое ребро
		os:=StructuralCopy(l1);
	elif fk.faces[1][l1-1] = os then
		os:=StructuralCopy(l1 - 1);
	else
		os:=StructuralCopy(rebro);
	fi;
	p:=fk.faces[1][os];
	Append(fk.1knot.sheets, [l1-1,l1]);
	fk.1knot.sheets:=Difference(fk.1knot.sheets,[os]);
	if number = 0 then
		# исправление ошибки с усиками
		for name in RecNames(fk.1knot.dpoints) do
			order:=fk.1knot.dpoints.(name);
			v:=Int(name);
			star:=StarFace(fk,[0,v]).1;
			star:=Intersection(star, fk.1knot.sheets);
			netak:=Difference(order,star);
			if Length(netak)=0 then
			else
				netak:=netak[1];
				pos:=Position(order,netak);
				a:=Difference(star,order)[1];
				order[pos]:=a;
				fk.1knot.dpoints.(name):=order;
			fi;
		od;

		s1:=rec( vertices:=[1,2], faces:=[ [ [1,2],[1,2] ] ]);
		ls1:=rec(0:=2, 1:=2);
		lfk:=rec(	0:=Length(fk.vertices),
					1:=Length(fk.faces[1]),
		 			2:=Length(fk.faces[2])		);
		pol:=PolProduct(s1,fk);
		2sheets:=[];
		2dpoints:=rec();
		for i in fk.1knot.sheets do
			images:=List([1,2],x->ImageInPolProduct(ls1,lfk,[[1,x],[1,i]])[2]);
			Append(2sheets,images);
		od;
		for v in RecNames(fk.1knot.dpoints) do
			for s in [1,2] do
				name:=ImageInPolProduct(ls1,lfk,[[1,s],[0,Int(v)]])[2];
				ord:=List(fk.1knot.dpoints.(v),
					  x -> ImageInPolProduct(ls1,lfk,[[1,s],[1,x]])[2]);
				2dpoints.(name):=StructuralCopy(ord);
			od;
		od;

		# создаине 3-диска из полнотория
		north:=pol.faces[1][os][1];
		south:=pol.faces[1][os][2];
		for v in pol.faces[1][os] do
			kl:=List([1,2],x->ImageInPolProduct(ls1,lfk,[[1,x],[0,v]])[2]);
			Add(pol.faces[2],kl);
		od;
		new2kl:=[-1,0] + Length(pol.faces[2]);
		Append(2sheets, new2kl);
		kl:=List([1,2],x->ImageInPolProduct(ls1,lfk,[[1,x],[1,os]])[2]);
		Append(kl, new2kl );
		Add(pol.faces[3], kl);
		
	else
	#--- создание 2-диска с необходимым разбиенеем -----------------------------
		# находим индексы внутренних 2-клеток узла
		2kl_os:=StarFace(fk,[1,os]).2[1];
		vnut_2kl_knot:=Difference([1 .. Length(fk.faces[2])],[2kl_os]);
		# Ниже, интересующее нас разбиение мы создаем вручную.
		new1kl:=[rebro, l1-1, l1];
		new1kl:=Difference(new1kl,[os]);
		if p[1] in fk.faces[1][new1kl[1]] then
			north_r:=new1kl[1];
			south_r:=new1kl[2];
		else
			north_r:=new1kl[2];
			south_r:=new1kl[1];
		fi;
		north:=Difference(fk.faces[1][north_r],[p[1]])[1];
		south:=Difference(fk.faces[1][south_r],[p[2]])[1];
		# тут можно устроить проверку поскольку списки north и south должны быть
		# одноэлементными.
		for i in [1 .. 3] do
			repeat
				name:=name+1;
				verify:= not (name in verknot);
			until verify;
			fk:=DivideFace(fk,[1,north_r],name);
			l1:=l1+1;
			Add(new1kl, l1);
		od;
		for i in [1,2] do
			repeat
				name:=name+1;
				verify:= not (name in verknot);
			until verify;
			fk:=DivideFace(fk,[1,south_r],name);
			l1:=l1+1;
			Add(new1kl,l1);
		od;
		point:=p[1];
		cepi:=[];
		# При создании новых точек, было созданно 7 клеток (осевую клетку мы
		# не рассматриваем). Эти 7 клеток распадаются на две части, условно мы
		# их обозначаем как северная часть и южная. В первую часть 1-клеток
		# входят четыре ребра, но нас итересуют только первые 3 ребра от точки
		# p[1], для того что бы восстановить порядок созданных вершин.
		while Length(new1kl)>4 do
			1kl:=Remove(new1kl,1);
			if point in fk.faces[1][1kl] then
				point:=Difference(fk.faces[1][1kl],[point])[1];
				Add(p, point);
				Add(cepi, 1kl);
			else
				Add(new1kl, 1kl);
			fi;
		od;
		# Аналогично с южной стороны для восстановления порядка вершин, нас
		# будет интересовать составленная из первых двух ребер. В настоящий
		# момент в списке new1kl должно было остаться 4 ребра.
		point:=south;
		while Length(new1kl)>2 do
			1kl:=Remove(new1kl,1);
			if point in fk.faces[1][1kl] then
				point:=Difference(fk.faces[1][1kl],[point])[1];
				Add(p,point);
				Add(cepi,1kl);
			else
				Add(new1kl,1kl);
			fi;
		od;

		Append(bord, cepi);
		Append(bord, new1kl);
		bord:=Difference(bord,[cepi[1]]);
		Add(fk.faces[1],Set(p{[2,3]}));
		p2p3:=Length(fk.faces[1]);
		Add(bord,p2p3);
		Add(fk.faces[2],bord);
		l2:=Length(fk.faces[2]);
		fk:=DivideFace(fk,[2,l2],Set(p{[5,7]}));
		p5p7:=p2p3+1;
		l0:=Length(fk.vertices);
		for 1kl in [p2p3, p5p7] do
			repeat
				name:=name+1;
				verify:= not (name in verknot);
			until verify;
			fk:=DivideFace(fk,[1,1kl],name);
			l0:=l0+1;
			Add(p,l0);
		od;

		l2:=Length(fk.faces[2]);
		new2kl:=Difference([1 .. l2], vnut_2kl_knot); 
		sost_new2kl:=List(new2kl, x -> FaceComp(fk,[2,x]).0);
		for 1kl in [ p{[8,9]}, p{[4,8]}, p{[4,6]}, p{[1,7]} ] do
			1kl:=Set(1kl);
			i:=0;
			repeat
				i:=i+1;
				verify:=Intersection(sost_new2kl[i],1kl) = 1kl;
			until verify;
			2kl:=new2kl[i];
			fk:=DivideFace(fk,[2,2kl],1kl);
			l2:=l2+1;
			Add(new2kl,l2);
			sost_new2kl[i]:=FaceComp(fk,[2,2kl]).0;
			Add(sost_new2kl,FaceComp(fk,[2,l2]).0);
		od;
		# Указываем индексы созданных 1-клеток и 2-клеток. Они необходимы для
		# того что бы в дальнейшем по ним пустить диаграмму узла. В данном
		# случае тут создается искусственное упорядочение.
		k:=[];
		sost_new2kl:=List(sost_new2kl, x -> Intersection(x,p));
		true_ord:=[p{[1,3,4,6,7]}, p{[4,5,6]}, p{[5,6,7,9]}, p{[2,7,8,9]},
								p{[4,5,8,9]}, p{[3,4,8]}, p{[1,2,7]}];
		true_ord:=List(true_ord, x -> Set(x));
		perm:=PermListList(true_ord, sost_new2kl);
		k:=Permuted(new2kl, perm);

		rebra:=Difference([1 .. Length(fk.faces[1])], rebra_knot);
		rebra:=Difference(rebra, [os]);
		sost_rebra:=List(rebra, x -> Intersection(fk.faces[1][x],p));
		true_ord:=[p{[1,3]}, p{[3,4]}, p{[4,5]}, p{[5]}, p{[6]}, p{[6,7]},
					p{[2,7]}, p{[3,8]}, p{[4,8]}, p{[5,9]}, p{[8,9]}, p{[2,8]},
					p{[7,9]}, p{[1,7]}, p{[4,6]}];
		true_ord:=List(true_ord, x -> Set(x));
		perm:=PermListList(true_ord, sost_rebra);
		r:=Permuted(rebra, perm);
		Add(fk.faces[1],Set(p{[1,2]}));
		Add(r,Length(fk.faces[1]));
		Add(fk.faces[2],Set(r{[1,8,12,16]}));
		Add(k,Length(fk.faces[2]));
		os:=fk.faces[1]{fk.faces[2][k[7]]};
		os:=Position(os,Set(p{[1,2]}));
		os:=fk.faces[2][k[7]][os];
		Add(r,os);	# ребро r[17] осевое
		# корректируем dpoints в северной и южной точках
		i:=4;
		for v in [north, south] do
			ord:=fk.1knot.dpoints.(v);
			pos:=fail;
			list:=r{[17,1,2,3,4,5,6,7]};
			j:=1;
			while pos=fail do
				pos:=Position(ord,list[j]);
				j:=j+1;
			od;
			ord[pos]:=r[i];
			i:=i+1;
		od;
	
	#--- создание полнотория с указанным 2-узлом -------------------------------
		kadr:=4;
		n:=AbsInt(number);
		s1:=rec(vertices:=[0 .. n*kadr-1], faces:=[[]]);
		s1.faces[1]:=List([1 .. n*kadr-1],x-> [x,x+1]);
		Add(s1.faces[1],[1,n*kadr]);
		ls1:=rec(0:=n*kadr, 1:=n*kadr, 2:=0, 3:=0);
		lfk:=rec(	0:=Length(fk.vertices),
					1:=Length(fk.faces[1]),
		 			2:=Length(fk.faces[2]),
					3:=0	);
		pol:=PolProduct(s1,fk);
		2sheets:=[];
		2dpoints:=rec();
		# Считаем образ узла на вертикальных клетках (т.е среди 2-клеток
		# созданных при умножении на ребра окружности s1).
		vertical:=List([1 .. 4], x -> StructuralCopy(fk.1knot));
		Append(vertical[1].sheets,r{[1 .. 7]});
		Append(vertical[2].sheets,r{[1,2,3,4,5,15,9,12]});
		Append(vertical[3].sheets,r{[14,13,11,8,2,3,4,5,15,9,12]});
		Append(vertical[4].sheets,r{[14,13,10,4,5,15,9,12]});
		if number>0 then
			perm:=();
		else
			perm:=(4,1)(2,3);
		fi;
		vertical[2].dpoints.(p[4]):=Permuted(r{[15, 9, 2, 3]},perm);
		vertical[3].dpoints.(p[4]):=Permuted(r{[15, 9, 2, 3]},perm);
		vertical[3].dpoints.(p[8]):=Permuted(r{[ 9,12, 8,11]},perm);
		for i in [1 .. kadr] do
			list:=[0 .. n-1]*kadr + i;
			for j in list do
				1kl1:=[1,j];
				for s in Set(vertical[i].sheets) do
					1kl2:=[1,s];
					Add(2sheets,ImageInPolProduct(ls1,lfk,[1kl1,1kl2])[2]);
				od;
				for v in RecNames(vertical[i].dpoints) do
					v:=Int(v);
					0kl2:=[0,v];
					name:=ImageInPolProduct(ls1,lfk,[1kl1,0kl2])[2];
					2dpoints.(name):=[];
					for s in vertical[i].dpoints.(v) do
						1kl2:=[1,s];
						Add(2dpoints.(name),
								ImageInPolProduct(ls1,lfk,[1kl1,1kl2])[2]);
					od;
				od;
			od;
		od;		
		# Считаем клетки 2-узла на горизонтальных уровнях (т.е.
		# уровни получающихся при умножении 2-клеток диска fk на вершины
		# отрезка окружности s1).
		level_sheets:=[];
		Add(level_sheets,Union(vnut_2kl_knot, k{[2,3,4,5]}));
		Add(level_sheets,Union(vnut_2kl_knot, k{[1,2,3,5,6]}));
		Add(level_sheets, k{[5,6]});
		for nk in [0 .. n-1] do
			for l in [2,3,4] do
				i:=nk*kadr+l;
				lsheets:=level_sheets[l-1];
				0kl:=[0,i];
				for kl in lsheets do
					2kl:=[2,kl];
					Add(2sheets,ImageInPolProduct(ls1,lfk,[0kl,2kl])[2]);
				od;
			od;
		od;
		# По построению звезда каждого горизонтального ребра в политопе pol
		# состоит ровно из 4-х 2-клеток. Если они все являются клетками узла, то
		# следовательно рассматриваемое ребро является ребром двойных точек.
		# Так же по построению горизонтальные клетки имеют меньший индекс чем
		# вертикальные клетки (это можно понять разобравшись с программой
		# ImageInPolProduct). Это означает что, если ребро является ребром
		# двойных точек, то первые 2-е 2-клетки зведы являются горизонтальными
		# клетками, вторые две 2-клетки --- вертикальными. Воспользуемся
		# описанной выше информацией для того, что бы найти все горизонтальные
		# ребра двойных точек и указать отношение выше\ниже для них.
		orient:=[number>0, number>0,number<0, number<0];
		sost2kl:=pol.faces[2]{2sheets};
		1kl:=Concatenation(sost2kl);
		set1kl:=Set(1kl);
		chastota:=List(set1kl, x -> Length(Positions(1kl,x)));
		1kl:=set1kl{Positions(chastota,4)};
		1kl:=Difference(1kl,List(RecNames(2dpoints),x->Int(x)));
		for i in 1kl do
			order:=[];
			s:=0;
			j:=1;
			while s<4 do
				if i in sost2kl[j] then
					Add(order, 2sheets[j]);
					s:=s+1;
				fi;
				j:=j+1;
			od;
			order:=Set(order);
			# По построению список 1kl состоит только из горизонтальных
			# 1-клеток, т.к. все вертикальные двойные ребра мы выкинули.
			image:=PreimageInPolProduct(ls1,lfk,[1,i]);
			name:=((image[1][2] - 1) mod kadr) + 1;
			if orient[name] then
				2dpoints.(i):=order;
			else
				2dpoints.(i):=order{[4,3,2,1]};
			fi;
		od;

	#--- запуск цикла вырезаний ------------------------------------------------
		# В каждой вершине которая является 1 mod 4 сделаем разре вокруг
		# внутренних клеток узла на этих уровнях. В итоге получим тор T^2 в
		# качестве границы. Склеим этот тор в кольцо так, что бы клетки
		# составляющие листы 2-узла совпали.
		l:=1;
		l1:=Length(pol.faces[1]);
		l2:=Length(pol.faces[2]);
		l3:=Length(pol.faces[3]);
		for l in [0 .. n-1]*4 + 1 do
			#--- вырезание тора ------------------------------------------------
			tor:=List(k,x->ImageInPolProduct(ls1,lfk,[[0,l],[2,x]])[2]);
			granica:=SetOfFacesBoundary(pol,tor,2);
			for i in [1..Length(tor)] do
				l2:=l2+1;
				Add(tor,l2);
			od;
			set_kl:=[k,r{[1..15]},p{[3..9]}];
			dim:=2;
			for set in set_kl do
				set:=List(set,x->ImageInPolProduct(ls1,lfk,[[0,l],[dim,x]])[2]);
				for kl in set do
					pol:=PolMinusFaceDoublingMethod(pol,[dim,kl]);
				od;
				dim:=dim-1;
			od;
			l1:=Length(pol.faces[1]);

			#--- упрощение тора ------------------------------------------------
			# Мы получили тор, теперь его 2-клетки нужно разбить на две
			# логических части по меридианам
			set:=pol.faces[2]{tor};
			set:=List(set, x -> Difference(x,granica));
			kolco1:=ConnectedSubset(set);
			kolco2:=Difference([1..Length(set)],kolco1);
			kolco1:=tor{kolco1};
			kolco2:=tor{kolco2};
			# Нужно восстановить линии которые индуцирует 2-узел на этом торе.
			1kl_tor:=Set(Concatenation(pol.faces[2]{tor}));
			1kl_2knot:=Set(Concatenation(pol.faces[2]{2sheets}));
			lines:=Intersection(1kl_tor, 1kl_2knot);
			# Упрощение разбиений на кольцах тора.
			subsets:=[];
			3kl:=[];
			for kolco in [kolco1,kolco2] do
				# линии узла на данном кольце
				1kl_kolco:=Set(Concatenation(pol.faces[2]{kolco}));
				1kl_kolco:=Intersection(lines,1kl_kolco);
				1kl_kolco:=Difference(1kl_kolco,granica);
				ind:=ConnectedSubset(pol.faces[1]{1kl_kolco});
				d1:=[1kl_kolco{ind}];
				ind:=Difference([1..Length(1kl_kolco)],ind);
				Add(d1,1kl_kolco{ind});
				plus:=[];
				for diski in d1 do
					bord:=SetOfFacesBoundary(pol,diski,1);
					Add(pol.faces[1],bord);
					l1:=l1+1;
					Add(diski,l1);
					Add(pol.faces[2],diski);
					l2:=l2+1;
					Add(plus,l2);
				od;

				2kl_kolco:=List(pol.faces[2]{kolco},x->Difference(x,1kl_kolco));
				ind:=ConnectedSubset(2kl_kolco);
				d2:=[kolco{ind}];
				ind:=Difference([1..Length(kolco)],ind);
				Add(d2,kolco{ind});
				for diski in d2 do
					Append(diski,plus);
					bord:=SetOfFacesBoundary(pol,diski,2);
					Add(pol.faces[2],bord);
					l2:=l2+1;
					Add(diski,l2);
					Add(pol.faces[3],diski);
					l3:=l3+1;
				od;
				Add(subsets,[-1,0]+l2);
				Add(3kl,[-1,0]+l3);
			od;

			#--- склеивание тора, указание листов узла -------------------------
			pol:=VerticesRullGluePol(pol,subsets[1],subsets[2],2);
			l2:=l2-2;
			for i in [1,2] do
				para:=pol.faces[3]{3kl[i]};
				para:=Intersection(para[1],para[2]);
				Append(2sheets,para);
				# При создании нужной нам границы нужно указать те 2-клетки
				# которые являюстя клетками 2-узла. Алгоритм специально построен
				# так, что на данном уровне не возникает горизонтальных ребер
				# двойных точек. Клетки 2-узла это те клетки которые получаются
				# на пересечении 3-клеток упрощающих разбиение тора.
			od;

		od;
	#--- превращаем полноторие в 3-диск ----------------------------------------
		# находим образ осевого ребра и его вершин
		os:=[1,r[17]];	# осевое ребро
		osa:=[0,p[1]];	# вершина оси
		osb:=[0,p[2]];	# вершина оси
		list:=List([1 .. kadr*n], x -> [1,x]);
		a:=List(list,x->ImageInPolProduct(ls1,lfk,[x,osa])[2]);
		b:=List(list,x->ImageInPolProduct(ls1,lfk,[x,osb])[2]);
		bok:=List(list,x->ImageInPolProduct(ls1,lfk,[x,os])[2]);
		Add(pol.faces[2],a);
		l2:=l2+1;
		Add(2sheets,l2);
		Add(bok,l2);

		Add(pol.faces[2],b);
		l2:=l2+1;
		Add(2sheets,l2);
		Add(bok,l2);

		Add(pol.faces[3],bok);
	fi;

	pol.2knot:=rec();
	pol.2knot.sheets:=2sheets;
	pol.2knot.dpoints:=2dpoints;
	bord:=PolBoundary(pol);
	Add(pol.faces[3], bord);

	#--- урощение полученного политопа -----------------------------------------
	Print("I'm trying simplify a polytope.\n\n");
	massa:=Length(pol.vertices)+Sum(List([1..3], x->Length(pol.faces[x])));
	a:=StructuralCopy(massa);
	verify:=false;
	while not verify do
		pol:=PolSimplifyWith2Knot(pol);
		newmassa:=Length(pol.vertices)+Sum(List([1..3],
			 x->Length(pol.faces[x])));
		verify := massa = newmassa;
		massa:=newmassa;
	od;


##  	Print("\n All good!\n\n");

return pol;
end);

################################################################################

# <ManSection><Func Name="TripleDoubleBranchPoints" Arg="pol" />
# 	<Description>
# 		Для шарового разбиения многообразия, внутри которого указана диаграмма
# 		заузленной поверхнсоти вычисляются те вершины политопа, которые
# 		являются тройными точками, двойными точками или точками ветвления. По
# 		возможности для каждой точки указываются 2-клетки диаграммы
# 		группированные по принадлежности различным листам. Функция выводит
# 		именованный список с полями .triple, .double и .branch. В именованном
# 		списке result.triple каждой вершине <M>v</M> под списком поля .u можно
# 		узнать 2-клетки звезды на верхнем листе в тройной точке <M>v</M>, под
# 		списком .m --- 2-клетки на среднем листе и под списком .d --- на
# 		нижнем. Для двойной точки указываются только списки .u и .d. Для точки
# 		ветвления указан только список 2-клеток узла лежащих в звезде этой
# 		вершины.
#		<Example>
# gap> pol:=TurnKnot(Trefoil,1);;
# ...
# gap> TripleDoubleBranchPoints(pol);
#rec( branch := rec( 13 := [ 2, 6, 22, 30 ], 19 := [ 13, 16, 37, 44 ] ),
#  double :=
#    rec( 1 := rec( d := [ 24, 25, 41, 57 ], u := [ 22, 57, 51, 27 ] ),
#      10 := rec( d := [ 28, 37, 38 ], u := [ 2, 30, 7 ] ),
#      11 := rec( d := [ 28, 36, 38 ], u := [ 2, 6, 7 ] ),
#      2 := rec( d := [ 23, 28, 55, 52, 55 ], u := [ 25, 41, 26 ] ),
#      21 := rec( d := [ 36, 38, 43 ], u := [ 8, 12, 14 ] ),
#      23 := rec( d := [ 30, 39, 48 ], u := [ 8, 12, 16 ] ),
#      24 := rec( d := [ 7, 46, 47 ], u := [ 13, 44, 14 ] ),
#      25 := rec( d := [ 24, 40, 41 ], u := [ 27, 39, 51 ] ),
#      26 := rec( d := [ 26, 42, 27 ], u := [ 23, 24, 40 ] ),
#      27 := rec( d := [ 23, 43, 52 ], u := [ 26, 42, 41 ] ),
#      30 := rec( d := [ 47, 48, 51, 57 ], u := [ 19, 44, 20, 43 ] ),
#      34 := rec( d := [ 46, 54, 57, 47 ], u := [ 14, 20, 19, 44 ] ) ),
#  triple :=
#    rec(
#      14 := rec( d := [ 32, 33, 40, 41 ], m := [ 27, 30, 35, 39 ],
#          u := [ 8, 9, 11, 12 ] ),
#      15 := rec( d := [ 27, 34, 35, 42 ], m := [ 23, 31, 32, 40 ],
#          u := [ 9, 10, 11, 12 ] ),
#      16 := rec( d := [ 23, 31, 36, 43 ], m := [ 33, 34, 41, 42 ],
#          u := [ 8, 9, 10, 12 ] ),
#      20 := rec( d := [ 37, 38, 43, 44 ], m := [ 7, 30, 47, 48 ],
#          u := [ 8, 13, 14, 16 ] ),
#      5 := rec( d := [ 24, 25, 32, 33 ], m := [ 22, 27, 30, 35 ],
#          u := [ 2, 3, 5, 6 ] ),
#      6 := rec( d := [ 26, 27, 34, 35 ], m := [ 23, 24, 31, 32 ],
#          u := [ 3, 4, 5, 6 ] ),
#      7 := rec( d := [ 23, 28, 31, 36 ], m := [ 25, 26, 33, 34 ],
#          u := [ 2, 3, 4, 6 ] ) ) )
#		</Example>
#	</Description>
# </ManSection>

InstallGlobalFunction(TripleDoubleBranchPoints, function(pol)
	local	dpoints, double_rebras, verss, ver, kratnost, triple, branch, star,
	inner1kl, components, spiski, i, ind, sheets, point, sp, level, 1kl, ab,
	cd, v, uitni, dpoint, dp;

	# Построим зведу вершины внутри диаграммы 2-узла. Далее,
	# 2-клетки звезды распадаются на 12 компонент связности, если не
	# учитывать при вычислении связности двойные ребра. Объединение компонент
	# связности проводим по соотношениям на двойных ребрах, как и отношение
	# высот.

	sheets:=pol.2knot.sheets;
	dpoints:=RecNames(pol.2knot.dpoints);
	double_rebras:=List(dpoints, x -> Int(x));
	verss:=List(double_rebras, x -> pol.faces[1][x]);
	verss:=Concatenation(verss);
	ver:=Set(verss);

	kratnost:=List(ver, v -> Length(Positions(verss,v)) );
	triple:=ver{Positions(kratnost,6)};
	branch:=ver{Positions(kratnost,1)};

	#--- тройные точки ---------------------------------------------------------
	point:=rec();
	for v in triple do
		star:=StarFace(pol,[0,v]);
		star.2:=Intersection(star.2, sheets);
		inner1kl:=Difference(star.1, double_rebras);
		components:=List(pol.faces[2]{star.2},x->Intersection(x,inner1kl));
		# те клетки которые не имеют ни одного внутреннего ребра являются
		# самостоятельными компонентами.

		# те клетки которые не имеют ни одного внутреннего ребра являются
		# самостоятельными компонентами.
		spiski:=[];
		while not IsEmpty(components) do
				ind:=ConnectedSubset(components);
				Add(spiski,star.2{ind});
				while not IsEmpty(ind) do
					i:=Remove(ind);
					Remove(components, i);
					Remove(star.2, i);
				od;
		od;
		star.1:=Intersection(star.1, double_rebras);
		point.(v):=rec(u:=[],m:=[],d:=[]);
		for sp in spiski do
			level:=[];
			for 1kl in star.1 do
				ab:=pol.2knot.dpoints.(1kl){[1,2]};
				cd:=pol.2knot.dpoints.(1kl){[3,4]};
				if not IsEmpty(Intersection(ab,sp)) then
					Add(level,1);
				elif not IsEmpty(Intersection(cd,sp)) then
					Add(level,-1);
				fi;
			od;
			# Если строка level равна [1,1] то рассматриваемая компонетна
			# принадлежит верхнему листу, если [-1,1] или [1,-1] - то среднему,
			# и если [-1,-1] - то нижнему.
			level:=Set(level);
			if level=[1] then
				Append(point.(v).u, sp);
			elif level=[-1] then
				Append(point.(v).d, sp);
			else
				Append(point.(v).m, sp);
			fi;
		od;
	od;

	#--- точки ветвления -------------------------------------------------------
	uitni:=rec();
	for v in branch do
		star:=StarFace(pol,[0,v]).2;
		star:=Intersection(star,sheets);
		uitni.(v):=StructuralCopy(star);
	od;

	#--- двойные точки ---------------------------------------------------------
	dpoint:=Set(Concatenation(pol.faces[1]{double_rebras}));
	dpoint:=Difference(dpoint,triple);
	dpoint:=Difference(dpoint,branch);
	dp:=rec();
	for v in dpoint do
		star:=StarFace(pol,[0,v]);
		star.2:=Intersection(star.2, sheets);
		inner1kl:=Difference(star.1, double_rebras);
		components:=List(pol.faces[2]{star.2},x->Intersection(x,inner1kl));
		spiski:=[];
		while not IsEmpty(components) do
				ind:=ConnectedSubset(components);
				Add(spiski,star.2{ind});
				while not IsEmpty(ind) do
					i:=Remove(ind);
					Remove(components,i);
					Remove(star.2,i);
				od;
		od;
		star.1:=Intersection(star.1, double_rebras);
		dp.(v):=rec(u:=[],d:=[]);
		1kl:=star.1[1];
		ab:=pol.2knot.dpoints.(1kl){[1,2]};
		for sp in spiski do
			level:=Intersection(ab,sp);
			if IsEmpty(level) then
				Append(dp.(v).d, sp);
			else
				Append(dp.(v).u, sp);
			fi;
		od;
	od;



		
return rec(triple:=point, double:=dp, branch:=uitni);
end);
################################################################################
# <ManSection><Func Name="2KnotInS4" Arg="pol" />
# 	<Description>
# 		По трехмерному политопу <M>pol</M> в котором содержится диаграмма
# 		заузленной двумерной поверхности <M>T</M> создается вложение этой
# 		поверхности <M>T</M> в четырехмерную сферу <M>S^4.</M>
#		<Example>
#		</Example>
#	</Description>
# </ManSection>
#
# NOTE: Тройные вершины помещаются в точки 1, 1/2, 2, двойные вершины в точки 1
# и 2, обычные вершины узла и точки ветвления в точку 1.

InstallGlobalFunction(2KnotInS4, function(s3)
	local	I, pol, ls3, singularity, tdbp, triple, double, branch, all1kl,
	allver, imagesver, v, l0, name, im1, images1kl, l2, l1, kl, bord1, bord2,
	im2, perm, 2knot, bord, im3, 2kl, alternativa, pos, 1kl, new1kl, list,
	heits, ver, ab, 1kls, sost;

	I:=ballAB(1);
	I.vertices:=[0,1];
	pol:=PolProduct(I,s3);
	ls3:=LengthPol(s3);
	singularity:=SingularitySet2Knot(s3);
	tdbp:=TripleDoubleBranchPoints(s3);
	triple:=List(RecNames(tdbp.triple), x -> Int(x));
	double:=List(RecNames(tdbp.double), x -> Int(x));
	branch:=List(RecNames(tdbp.branch), x -> Int(x));
	all1kl:=Set(Concatenation(s3.faces[2]{s3.2knot.sheets}));
	allver:=Set(Concatenation(s3.faces[1]{all1kl}));
	imagesver:=[];
	for v in allver do
		imagesver[v]:=[ImageInPolProduct(I,ls3,[[0,1],[0,v]])[2]];
	od;
	for v in double do
		Add(imagesver[v],ImageInPolProduct(I,ls3,[[0,2],[0,v]])[2]);
	od;
	l0:=Length(pol.vertices);
	for v in triple do
		Add(imagesver[v],ImageInPolProduct(I,ls3,[[0,2],[0,v]])[2]);
		name:=[1/2, s3.vertices[v]];
		im1:=ImageInPolProduct(I,ls3,[[1,1],[0,v]]);
		pol:=DivideFace(pol,im1,name);
		l0:=l0 + 1;
		Add(imagesver[v], l0, 2);
	od;

	#--- поднимаем ребра 2-узла ------------------------------------------------
	images1kl:=[];
	l2:=Length(pol.faces[2]);
	l1:=Length(pol.faces[1]);
	for kl in singularity.graf do 
		ver:=pol.faces[1][kl];
		bord1:=[];
		bord2:=[];
		for v in ver do
			if v in triple then
				pos:=Position(singularity.order.(v), kl);
				if pos < 3 then
					Add(bord1, imagesver[v][3]);
					Add(bord2, imagesver[v][2]);
				elif pos < 5 then
					Add(bord1, imagesver[v][3]);
					Add(bord2, imagesver[v][1]);
				else
					Add(bord1, imagesver[v][2]);
					Add(bord2, imagesver[v][1]);
				fi;
			elif v in branch then
				Add(bord1, imagesver[v][1]);
				Add(bord2, imagesver[v][1]);
			else
				Add(bord1, imagesver[v][2]);
				Add(bord2, imagesver[v][1]);
			fi;
		od;
		# по построению bord1 > bord2
		im2:=ImageInPolProduct(I,ls3,[[1,1],[1,kl]]);
		sost:=pol.faces[1]{pol.faces[2][im2[2]]};
		images1kl[kl]:=[];
		for bord in [Set(bord1), Set(bord2)] do
			if bord in sost then
				pos:=Position(sost, bord);
				Add(images1kl[kl], pol.faces[2][im2[2]][pos]);
			else
				pol:=DivideFace(pol,im2,bord);
				l2:=l2+1;
				l1:=l1+1;
				Add(images1kl[kl], l1);
				# Вертикальная 2-клетка была раздроблена на две части. Что бы в
				# дальнейшем знать индексацию клеток сделаем так, чтобы нижняя
				# 2-клетка имела индекс изначальной клетки, поскольку ее мы
				# будем разбивать снова (если это только первый шаг).
				im1:=ImageInPolProduct(I,ls3,[[0,1],[1,kl]])[2];
				if im1 in pol.faces[2][im2[2]] then
					perm:=();
				else
					perm:=(im2[2],l2);
				fi;
				pol:=PermFaces(pol, perm, 2);
			fi;
		od;
	od;
	# Для остальных ребер образ будем вычислять по ходу поднятия 2-клеток,
	# причем очевидно, что полностью внутренние ребра диаграммы поднимаются
	# однозначно в точку 1. 
	
	#--- поднимаем грани 2-узла ------------------------------------------------
	l2:=Length(pol.faces[2]);
	2knot:=[];
	for 2kl in s3.2knot.sheets do
		bord:=[];
		1kls:=s3.faces[2][2kl];
		alternativa:=[[],[]];
		for 1kl in 1kls do
			if 1kl in singularity.graf then
				ab:=s3.2knot.dpoints.(1kl){[1,2]};
				if 2kl in ab then
					Add(bord, images1kl[1kl][1]); # верхняя
				else
					Add(bord, images1kl[1kl][2]); # нижняя
				fi;
			else
				ver:=s3.faces[1][1kl];
				heits:=[];
				for v in ver do
					if v in triple then
						if 2kl in tdbp.triple.(v).u then
							Add(heits,3);
						elif 2kl in tdbp.triple.(v).d then
							Add(heits,1);
						else
							Add(heits,2);
						fi;
					elif v in double then
						if 2kl in tdbp.double.(v).u then
							Add(heits,2);
						else
							Add(heits,1);
						fi;
					else
						# Во всех остальных случаях (внутренняя точка и точка
						# ветвления) вершина имеет единственный образ в
						# поднятии.
						Add(heits,1);
					fi;
				od;
				new1kl:=[imagesver[ver[1]][heits[1]],
			 				imagesver[ver[2]][heits[2]]];
			 	new1kl:=Set(new1kl);
				im2:=ImageInPolProduct(I,ls3,[[1,1],[1,1kl]])[2];
				list:=pol.faces[1]{pol.faces[2][im2]};
				if new1kl in list then
					pos:=Position(list, new1kl);
					Add(bord, pol.faces[2][im2][pos]);
				else
					pol:=DivideFace(pol,[2,im2], new1kl);
					l1:=l1+1;
					l2:=l2+1;
					Add(bord, l1);
				fi;
			fi;
			Add(alternativa[1], ImageInPolProduct(I,ls3,[[0,1],[1,1kl]])[2]);
			Add(alternativa[2], ImageInPolProduct(I,ls3,[[0,2],[1,1kl]])[2]);
		od;
		alternativa:=List(alternativa, x -> Set(x));
		bord:=Set(bord);
		if bord in alternativa then
			pos:=Position(alternativa, bord);
			Add(2knot, ImageInPolProduct(I,ls3,[[0,pos],[2,2kl]])[2]);
		else
			im3:=ImageInPolProduct(I,ls3,[[1,1],[2,2kl]]);
			pol:=DivideFace(pol,im3, bord);
			l2:=l2+1;
			Add(2knot, l2);
		fi;
	od;

	pol.2knot:=2knot;
	# на настоящий момент pol это I x S3 сделаем из нее полноценную 4-сферу
	bord:=PolBoundary(pol);
	ab:=ConnectedSubset(pol.faces[3]{bord});
	ab:=bord{ab};
	bord:=Difference(bord,ab);
	Add(pol.faces[4],Set(ab));
	Add(pol.faces[4],bord);

return pol;
end);

################################################################################

# <ManSection><Func Name="IsDiagrammOf2Knot" Arg="pol" />
# 	<Description>
# 		Проверяет диаграмму 2-узла вложенную в трехмерное многообразие. Для
# 		этого проверяются: 1) корректность всех ссылок на клетки политопа, 2)
# 		граф двойных точек, 3) отсутствие точек самокасания.
#		<Example>
# gap> pol:=TurnKnot(Figure8,-1);;
# gap> IsDiagrammOf2Knot(pol);
# true
#		</Example>
#	</Description>
# </ManSection>


InstallGlobalFunction(IsDiagrammOf2Knot, function(pol)
	local	1kl, names, vert, chastota, verify, knot, l, v, max, l1, l2, name,
	uslovie1, uslovie2, star, uslovie3;

	# Для того, что бы проверить на диаграммность можно проверить следующее
	# 0) все ссылки в прикрепленных данных .2knot действительны
	# 1) граф двойных точек имеет вершины валентности 1,2,6
	# 2) зведы всех 1-клеток заузленной поверхности содержат 1,2 или 4 2-клетки
	# 3) Строится прообраз двумерной заузленной поверхности.
	# 4) Проверяем, что в получившейся поверхности нет самокасаний.

	#--- проверка что все ссылки в .2knot действительные -----------------------
	max:=Maximum(pol.2knot.sheets);
	l2:=Length(pol.faces[2]);
	if max <= l2 then
		verify:=true;
	else
		verify:=false;
		Print("Ошибка в списке pol.2knot.sheets, не существующая клетка\n");
	fi;
	1kl:=List(RecNames(pol.2knot.dpoints), x -> Int(x));
	max:=Maximum(1kl);
	l1:=Length(pol.faces[1]);
	if (max <= l1) and verify then
	else
		verify:=false;
		Print("Ошибка в именованном списке pol.2knot.dpoints, ");
		Print("не сущствующее вдойное ребро.\n");
	fi;

	for name in RecNames(pol.2knot.dpoints) do
		uslovie1 := (Length(pol.2knot.dpoints.(name)) = 4);
		uslovie2 := IsSubset(pol.2knot.sheets, pol.2knot.dpoints.(name));
		uslovie3 := Set(List(pol.2knot.dpoints.(name), 
			  x -> (Int(name) in pol.faces[2][x])))  = [true];
		if uslovie1 and uslovie2 and uslovie3 then
		else
			verify:=false;
			Print("Ошибка в построении списка pol.2knot.dpoints.", Int(name),
			"\n");
		fi;
	od;

	#--- проверка валентности вершин графа -------------------------------------
	if verify then
		names:=Concatenation(pol.faces[1]{1kl});
		vert:=Set(names);
		chastota:=List(vert, x -> Length(Positions(names,x)));
		chastota:=Set(chastota);
		if IsSubset([1,2,6],chastota) then
			verify:=true;
		else
			verify:=false;
			Print("Нарушена валентность вершин на графе двойных точек.\n");
		fi;
	fi;

	if verify then
		knot:=SubPolytope(pol,pol.2knot.sheets,2);
		for 1kl in [1 .. Length(knot.faces[1])] do
			l:=Length(StarFace(knot,[1,1kl]).2);
			verify:=verify and (l in [1,2,4]);
		od;
	fi;
	if not verify then
		Print("В каком-то ребре диаграммы сходятся три или более листа.\n");
	fi;

	if verify then
		knot:=SurfaceOf2Knot(pol);
		for v in [1 .. Length(knot.vertices)] do
			star:=StarFace(knot,[0,v]).2;
			star:=knot.faces[2]{star};
			l:=Length(ConnectedSubset(star));
			if l = Length(star) then
			else
				verify:=false;
			fi;
		od;
		verify:=(knot = SurfaceOf2Knot(pol));
	fi;
	if not verify then
		Print("В диаграмме есть недопустимые точки (самокасание).\n");
	fi;

return verify;
end);

###############################################################################
#### PL-2.4
###############################################################################

#			<ManSection><Func Name="KnotGroup" Arg="knot" />
#				<Description>
#					Вычисляется фундаментальная группа узла, которая
#					определяется как фундаментальная группа дополнения узла в
#					трехмерной сфере <M>S^3.</M> Для создания группы
#					используются соотношения Виртингера.
#					<Example>
#					</Example>
#				</Description>
#			</ManSection>

InstallGlobalFunction(KnotGroup, function(knot)
	local	n, genf, lgenf, i, gr, generators, relators, namepoint, slovo, fg,
	lr, r;

	# создание свободной группы
	n:=Length(knot.kod);
	r:=n/2;
	genf:=[1];
	lgenf:=1;
	for i in [2 .. n-1] do
		if knot.kod[i][2] = -1 then
			lgenf:=lgenf+1;
		fi;
		Add(genf,lgenf);
	od;
	if knot.kod[n][2] = -1 then
		Add(genf, 1);
	else
		Add(genf,lgenf);
	fi;
	gr:=FreeGroup(r);
	generators:=GeneratorsOfGroup(gr);

	# создание соотношений группы
	relators:=[];
	for i in [1 .. n] do
		if knot.kod[i][2] = 1 then
			namepoint:=knot.kod[i][1];
			lr:=LeftRight(knot,i,namepoint);
			Add(relators, [ i, lr.right, i, lr.left ]);
			# степени у соответствующих образующих [1,1,-1,-1]
		fi;
	od;

	for i in [1 .. r] do
		slovo:=genf{relators[i]};
		slovo:=generators{slovo};
		slovo[3]:=slovo[3]^(-1);
		slovo[4]:=slovo[4]^(-1);
		relators[i]:=Product(slovo);
	od;

	# факторизация и упрощение
	gr:=gr/relators;
	fg:=PresentationFpGroup(gr);
	TzGo(fg);

return FpGroupPresentation(fg);
end);

################################################################################

#			<ManSection><Func Name="TorusKnot" Arg="q,p" />
#				Создается диаграмма торического узла <M>(q,p),</M> если <M>q</M>
#				и <M>p</M> взаимнопростые, если это не так, то будет создано
#				соответствующее зацепление. Параметр <M>q > 0</M> соответствует
#				количеству нитей, а параметр <M>p</M> соответствует количеству
#				оборотов. 
#				<Description>
#					<Example>
#					</Example>
#				</Description>
#			</ManSection>

InstallGlobalFunction(TorusKnot, function(q,p1)
	local	t, name, i, ind, next, l, d, h, exit, p, ij, j, orient, kod, sgn;
	
	p:=AbsInt(p1);
	t:=List([1 .. p], i -> []);
	name:=1;
	for i in [1 .. p] do
		t[i][1]:=[];
		for j in [2 .. q] do
			Add(t[i][1], [name,1]);
			t[i][j]:=[[name, -1]];
			name:=name+1;
		od;
	od;

	# Вычислим все траектории по которым будет гулять торический узел
	ind:=[ [0,0] ];
	next:= [1, -1] mod [p,q];
	l:=1;
	while not next = ind[1] do
		Add(ind, next);
		l:=l+1;
		next:=(ind[l] + [1,-1]) mod [p,q];
	od;
	sgn:=SignInt(p1);
	orient:=List([1 .. (q-1)*p], i -> [i, sgn]);
	if l = p*q then # тогда это узел
		kod:=[];
		for ij in ind do
			i:=ij[1]+1;
			j:=ij[2]+1;
			Append(kod, t[i][j]);
		od;
	else # тогда это зацепление
		d:=p*q/l;
		h:=0;
		kod:=rec();
		while h < d do
			h:=h+1;
			kod.(h):=[];
			for ij in ind do
				i:=ij[1]+1;
				j:=ij[2]+1;
				Append(kod.(h), t[i][j]);
			od;
			ind:=List(ind, x -> [x[1], (x[2]+1) mod q]);
		od;
	fi;

return rec(kod:=kod, orient:=orient);
end);

################################################################################

#			<ManSection><Func Name="SurfaceOf2Knot" Arg="M3" />
#				<Description>
#					Для заузленной поверхности указанной внутри некоторого
#					3-многообразия <M>M3</M> создается прообраз этой поверхности
#					с указанием прообразов двойных ребер, тройных точек и точек
#					ветвления. Прообраз создается как pl-комплекс, к которому
#					прикреплена дополнительная информация содеражащаяся в
#					именнованном списке <M>.preimage.</M> <M>Preimage</M>
#					является списком дублированных прообразов. В него входят
#					список <M>.1</M> и <M>.0.</M> Список
#					<M>.1</M> содержит прообразы для каждого двойного
#					ребра, первым элементом пары является ребро-прообраз лежащее
#					на нижнем листе, вторым, соответственно ребро-прообраз
#					лежащее на верхнем листе образа в диаграмме. Список
#					<M>.0</M> содержит состоит из списков длины 1 и 3,
#					которые является списками прообразов точек ветвления и
#					тройных точек, соответственно. Причем, для тройных точек
#					прообраз тройной точки, лежащий на верхнем листе будет
#					третьим в списке, на среднем - вторым и на нижнем,
#					соответственно, первым.
#				</Description>
#			</ManSection>

InstallGlobalFunction(SurfaceOf2Knot, function(s3)
	local	pol, preimage, 2kl, sost, dpoints, l1, vert, r, ind, order,
	sost1kl, i, heits, l0, v, triple, list, s, j, k, kl, image;

	pol:=SubPolytope(s3,s3.2knot.sheets,2);
	# После такого выделения внутренний порядок в подмножестве клеток, которые
	# образуют заузленную диаграмму сохраняется. Сейчас нужно добавить дубликаты
	# для двойных ребер и тройных точек.
	preimage:=rec(points:=[], rebras:=[]);
	image:=rec();

	#--- создание прообраза ----------------------------------------------------
	2kl:=s3.2knot.sheets;
	sost:=List(2kl, x -> FaceComp(s3,[2,x]).1);
	sost:=Set(Concatenation(sost));
	#... дублирование ребер ....................................................
	dpoints:=List(RecNames(s3.2knot.dpoints), x -> Int(x));
	image.rebras:=dpoints;
	l1:=Length(sost);
	vert:=[];
	for r in dpoints do
		ind:=Position(sost, r);
		order:=List(s3.2knot.dpoints.(r){[1,2]}, x -> Position(2kl, x));
		# Мы нашли двойное ребро на поверхности, теперь его надо продублировать
		# и верхнюю часть поверхности пустить по дубликату. В графе будем
		# указывать именно пару ребер на поверхности, которые образуют одно
		# двойное ребро как пару, первым элементом которой идет прообраз нижнего
		# ребра, вторым прообраз верхнего ребра.
		sost1kl:=pol.faces[1][ind];
		Add(pol.faces[1],sost1kl);
		Append(vert, sost1kl);
		l1:=l1+1;
		for i in order do
			kl:=Difference(pol.faces[2][i], [ind]);
			Add(kl, l1);
			pol.faces[2][i]:=StructuralCopy(kl);
		od;
		Add(preimage.rebras, [ind, l1]);
	od;
	#... дублирование тройных точек ............................................
	heits:=[];
	for r in preimage.rebras do
		for i in [1,2] do
			heits[r[i]]:=i;
		od;
	od;
	image.points:=[];
	vert:=Set(vert);
	l0:=Length(pol.vertices);
	for v in vert do
		pol:=PolMinusFaceDoublingMethod(pol,[0,v]);
		if IsBound(pol.vertices[l0+2]) then # это была тройная точка
			triple:=[v, l0+1, l0+2];
			Add(image.points,v);
			# Получили три прообраза тройной точки, нам надо упорядочить их
			# относительно того каким листам они принадлежат нижнему, среднему
			# или верхнему.
			list:=List(triple, x -> StarFace(pol,[0,x]).1);
			order:=[];
			for i in [1,2,3] do
				s:=[];
				for j in list[i] do
					if IsBound(heits[j]) then
						Add(s,heits[j]);
					fi;
				od;
				s:=Set(s);
				if s = [2] then
					order[3]:=i;
				elif s = [1] then
					order[1]:=i;
				else
					order[2]:=i;
				fi;
			od;
			triple:=triple{order};
			Add(preimage.points, triple);
			l0:=l0+2;
		elif IsBound(pol.vertices[l0+1]) then # это была двойная точка
			l0:=l0+1;
		else	# это может быть только точка ветвления
			Add(preimage.points, [v]);
			Add(image.points, v);
		fi;
		# Список preimage.points будет состоять только из списков длины 1 и 3,
		# которые отвечают точкам ветвления и тройным точкам, соответственно.
		# Причем для прообразов тройных точек, порядок соответствует высоте
		# листа в тройной точке, чем выше порядок, тем выше лист.
	od;
	
	list:=Set(Concatenation(s3.faces[2]{s3.2knot.sheets}));
	list:=Set(Concatenation(s3.faces[1]{list}));
	image.points:=list{image.points};
	pol.images:=rec(0:=image.points, 1:=image.rebras);
	pol.preimages:=rec(0:=preimage.points, 1:=preimage.rebras);

return pol;
end);
# Было проверенно, что прообраз является многообразием. Для spun трилистника и
# для 2-twist трилистника было показано, что эти прообразы сферы. Было проверено
# количество тройных точек и точек ветвления, так же для одной из тройных точек
# проверялось корректность упорядочения по принадлежности высотам.

################################################################################

#			<ManSection><Func Name="OrientBrockenDiagramm" Arg="s3" />
#				<Description>
#					Функция построит разорванную диаграмму для ориентируемой
#					линейно связной заузленной поверхности. 
#
#					Для ориентируемого 2-узла вычисляется информация о данной
#					диаграмме в которую входят разорванная диаграмма 2-узла и
#					граф двойных. На выход подается именованный список. 
#
#					Список rez.manifold содержит прообраз двумерной заузленной
#					поверхности.
#					
#					В списках rez.images и rez.preimages содержатся образ и
#					прообраз двойных линий и сингулярных точек (тройных и точек
#					ветвления), соответственно. Считается, что на двойных ребрах
#					введена нумерация которая отображается в списках
#					rez.images.1 и rez.preimages.1, то есть rez.images.1[i] и
#					rez.preimages.1[i] это образ и прообраз i-того двойного
#					ребра. На тройных точках и на точках ветвления введена общая
#					нумерация которая так же отображается в списках
#					rez.preimages.0  и rez.images.0.
#					
#					Граф сингулярных точек заузленной поверности состоит из
#					набора точек соответствующих тройным точкам и точкам
#					ветвления, при этом дуги графа могут проходить по нескольким
#					ребрам шарового разбиения. На дугах графа двойных точек
#					введена нумерация. Список rez.colines i-тому двойному ребру
#					сопоставляет дугу которая графа двойных точек которая
#					проходит по данному ребру.
#
#					В rez.cosheets для i-той 2-клетки заузленной поверхности
#					указывает какому листу разорванной диаграммы эта 2-клетка
#					принадлежит.
#
#					Так как по условиею прообраз заузленной поверхности является
#					ориентируемым многообразием, то дуги графа двойных точек
#					можно ориентировать, а так же присвоить ориентации вершинам
#					графа. Эти ориентации содержатся в списках rez.coorient.dim,
#					где dim=0 или 1. Заметим, что оринетации в rez.coorient.1
#					сопоставляются не по дугам и по двойным ребрам.
#
#					Список rez.cofaces.0 на i-том месте содержит списко дуг
#					графа двойных точек которые содержат i-тую вершину графа.
#					Причем если вершина является тройной точкой, то
#					соответствующий список состоит из шести элементов. Первая
#					тройка списка rez.cofaces.0[i] это дуги входящие в i-тую
#					тройную точку, оставшаяся тройка это дуги которые исходят из
#					указанной вершины.
#
#					Список rez.cofaces.1 для каждой дуги графа двойных точек
#					указывает список из трех элементов в которым первым
#					указывается номер верхнего листа, вторым индекс нижнего
#					листа по направлению нормали нормали верхнего листа, третьим
#					- индекс нижнего листа против направления номрали верхенго
#					листа.
#					<Example>
#					</Example>
#				</Description>
#			</ManSection>

InstallGlobalFunction(OrientBrockenDiagramm, function(s3)
	local	s2, s3_orient, s2_orient, cell_orient, dpoints, l2, ldp, downimage,
	2cells, indexis, cosheets, list, i, s, grafpoints, sostrebra, colines,
	sheets, ll, line_ind, r, triple, plate, 2kl, star_r, 3kletki, 3kl, set,
	normal, p, part, 3star_2kl, j, wedo, wewilldo, ind, para, ident, verify,
	ab, cd, sost, pos, info, dp_orient, direction, cofaces, a, b, enter_exit,
	l, data, order, v, 1kl, 1kls, star, height, tdbp, coorient, gamma;

	s2:=SurfaceOf2Knot(s3);
	s3_orient:=PolOrient(s3);
	cell_orient:=CellOrient(s3);
	dpoints:=s2.images.1;
	sheets:=StructuralCopy(s3.2knot.sheets);

	#--- разбиение на листы и дуги --------------------------------------------
	l2:=Length(s2.faces[2]);
	ldp:=Length(dpoints);
	#... разбиение на листы ...................................................
	downimage:=List(s2.preimages.1, x -> x[1]);
	2cells:=StructuralCopy(s2.faces[2]);
	2cells:=List(2cells, x -> Difference(x, downimage));
	indexis:=[1 .. l2];
	cosheets:=[];
	coorient:=rec(0:=[], 1:=[]);
	s:=1;
	while not IsEmpty(2cells) do
		list:=Set(ConnectedSubset(2cells));
		while not IsEmpty(list) do
			i:=Remove(list);
			Remove(2cells,i);
			i:=Remove(indexis, i);
			cosheets[i]:=StructuralCopy(s);
		od;
		s:=s+1;
	od;

	#... разбиение на дуги ....................................................
	grafpoints:=StructuralCopy(s2.images.0);
	sostrebra:=s3.faces[1]{dpoints};
	sostrebra:=List(sostrebra, x -> Difference(x, grafpoints));
	indexis:=[1 .. ldp];
	colines:=[];
	s:=1;
	while not IsEmpty(sostrebra) do
		list:=Set(ConnectedSubset(sostrebra));
		while not IsEmpty(list) do
			i:=Remove(list);
			Remove(sostrebra, i);
			i:=Remove(indexis, i);
			colines[i]:=StructuralCopy(s);
		od;
		s:=s+1;
	od;
	# NOTE: до сих проверенно

	#--- ориентирование листов, дуг и тройных точек ---------------------------
	#... ориентирование диаграммы узла ........................................
	sost:=StructuralCopy(s3.faces[2]{sheets});
	s2_orient:=[1];
	wedo:=[1];
	wewilldo:=[2 .. Length(sheets)];
	while not IsEmpty(wewilldo) do
		j:=Remove(wewilldo,1);
		list:=List(wedo, x -> Intersection(sost[x], sost[j]));
		pos:=List(list, x -> IsEmpty(x));
		pos:=Positions(pos, false);
		ident:=wedo{pos};
		verify:=not IsEmpty(ident);
		para:=[];
		while verify do
			ind:=Remove(ident);
			set:=Remove(list, Remove(pos));
			while not IsEmpty(set) and IsEmpty(para) do
				r:=Remove(set);
				para:=sheets{[j, ind]};
				if r in dpoints then
					ab:=s3.2knot.dpoints.(r){[1,2]};
					cd:=s3.2knot.dpoints.(r){[3,4]};
					if IsSubset(ab, para) or IsSubset(cd, para) then
					else
						para:=[];
					fi;
				fi;
			od;
			verify:=not IsEmpty(ident);
			if not IsEmpty(para) then
				verify:=false;
			fi;
		od;
		if not IsEmpty(para) then
			s2_orient[j]:= - s2_orient[ind];
			for i in para do
				p:=Position(s3.faces[2][i], r);
				s2_orient[j]:=s2_orient[j] * cell_orient[2][i][p];
			od;
			Add(wedo, StructuralCopy(j));
		else
			Add(wewilldo, j);
		fi;
	od;

	#..........................................................................

	ll:=StructuralCopy(s)-1;	# количество ребер графа двойных точек
	cofaces:=rec(1:=[], 0:=[]);
	dp_orient:=[];
	enter_exit:=List(s2.images.0, x -> [ [], [] ]);
	for line_ind in [1 .. ll] do
		indexis:=Positions(colines, line_ind);

	#... создание троек листов в ребрах .......................................
		r:=dpoints[indexis[1]];
		triple:=StructuralCopy(s3.2knot.dpoints.(r));
		plate:=StructuralCopy(triple{[1,2]});
		2kl:=Remove(triple,1);
		star_r:=StarFace(s3, [1,r]);
		3kletki:=s3.faces[3]{star_r.3};
		3kletki:=List(3kletki, x -> Difference(x, plate));
		# По построению, 3-клетки из star_r содержатся только либо по
		# направлению нормали к верхнему листу (клеткам plate), либо против
		# этого направления.
		3kletki:=List(3kletki, x -> Intersection(x, star_r.2));
		part:=ConnectedSubset(3kletki);
		3star_2kl:=StarFace(s3,[2,2kl]).3;
		3kl:=Intersection(3star_2kl, star_r.3{part})[1];
		normal:=s3_orient[3kl];
		p:=Position(s3.faces[3][3kl], 2kl);
		normal:=normal * cell_orient[3][3kl][p];
		p:=Position(sheets, 2kl);
		normal:=normal * s2_orient[p];
		set:=Union(3kletki{part});
		if (normal = 1) = (triple[2] in set) then
			triple:=triple{[1,3,2]};
		fi;
		triple:=List(triple, x -> Position(sheets, x));
		Add(cofaces.1, cosheets{triple});

	#... ориентирование ребер .................................................
		list:=dpoints{indexis};
		list:=s3.faces[1]{list};
		info:=LineOrdering(list);
		2kl:=sheets[triple[3]];
		p:=Position(s3.faces[2][2kl], r);
		direction:=StructuralCopy(cell_orient[2][2kl][p]);
		p:=Position(sheets, 2kl);
		direction:=direction * s2_orient[p];
		p:=Position(info.order,1);
		if not direction = info.orient[p] then
			info.orient:= - info.orient;
		fi;
		s:=1;
		for i in info.order do
			j:=indexis[i];
			coorient.1[j]:=info.orient[s];
			s:=s+1;
		od;

	#... вход-выход по тройным точкам .........................................
		list:=list{info.order};
	   if not IsEmpty(s2.images.0) then
		if Length(info.order)=1 and not IsEmpty(s2.images.0) then
			if info.orient[1]=1 then
				ab:=list[1];
			else
				ab:=list[1]{[2,1]};
			fi;
		else
			l:=Length(info.order);
			a:=Difference(list[1],list[2])[1];
			b:=Difference(list[l], list[l-1])[1];
			ab:=[];
			if not (a in s2.images.0) then
				# Если точка a не принадлежит тройным точкам, это должно
				# означать, что точка b тоже не принадлежит тройным точкам, а
				# следовательно a=b и данная линия является окружностью. Если на
				# данной окружности лежит тройная точка, то она должна быть
				# единственной (по потстроению). Для данной тройной точки,
				# данное ребро является как входящим так и исходящим.
				a:=Union(list);
				a:=Intersection(a, s2.images.0);
				if not IsEmpty(a) then
					ab:=a{[1,1]};
				fi;
			else
				if (list[1][1]=a) = (info.orient[1]=1) then
					ab:=[a,b];
				else
					ab:=[b,a];
				fi;
			fi;
		fi;
		for i in [1,2] do
			p:=Position(s2.images.0, ab[i]);
			Add(enter_exit[p][i],line_ind);
		od;
	   fi;
		# По построению в списке enter_exit первая группа это исходящие, вторая
		# группа это входящие двойные линии.
	#..........................................................................
		# сбор информации о ориентации двойных линий
		list:=indexis{info.order};
		while not IsEmpty(list) do
			i:=Remove(list);
			dp_orient[i]:=Remove(info.orient);
		od;	
	od;
	
	#... упорядочение по высотам в тройных точках .............................
	s:=1;
	tdbp:=TripleDoubleBranchPoints(s3).triple;
	for v in s2.images.0 do
		if Length(s2.preimages.0[s])=3 then
			star:=StarFace(s3,[0,v]);
			height:=[];
			1kls:=List(["d","m","u"],
				x->Concatenation(s3.faces[2]{tdbp.(v).(x)}));
##  			for j in [[1,2],[1,3],[2,3]] do # [md, ud, um, md, ud, um]
			for j in [[2,3],[1,3],[1,2]] do # [um, ud, md, um, ud, md]
				1kl:=Intersection(1kls{j});
				1kl:=Intersection(1kl, star.1);
				1kl:=List(1kl, x -> Position(dpoints, x));
				1kl:=colines{1kl};
				if 1kl[1] in enter_exit[s][1] and 1kl[2] in enter_exit[s][2] then
				else
					1kl:=1kl{[2,1]};
				fi;
				Append(height,1kl);
			od;
			Add(cofaces.0,height{[2,4,6,1,3,5]});
			# Список cofaces.0 устроен следующим образом. Первые три элемента
			# это индексы входящих в вершину двойных линий (в соответствии с
			# направлениями на них), вторая тройка числе это, соответственно,
			# исходящие.
	#... ориентирование тройных точек .........................................
			# Для тройной точки построить нормаль к верхнему листу. Направление
			# построенной нормали сравнить с направлением соответствующей
			# двойной линии (пересечение нижнего и среднего листов). Если
			# направления совпадают, то тройная точка положительна.
			star:=StarFace(s3,[0,v]);
			sost:=s3.faces[3]{star.3};
			sost:=List(sost, x -> Intersection(x, star.2));
			sost:=List(sost, x -> Difference(x, tdbp.(v).u));
			ind:=ConnectedSubset(sost);
			3kl:=star.3[1];
			2kl:=Intersection(s3.faces[3][3kl], tdbp.(v).u)[1];
			normal:=s3_orient[3kl];
			p:=Position(s3.faces[3][3kl], 2kl);
			normal:=normal*cell_orient[3][3kl][p];
			p:=Position(sheets, 2kl);
			normal:=normal*s2_orient[p];
			if normal=-1 then
				ind:=Difference([1..Length(sost)], ind);
			fi;
			sost:=Set(Concatenation(sost{ind}));
			1kl:=Intersection(s3.faces[2]{sost});
			1kl:=Intersection(1kl, star.1);
			1kl:=1kl[1];
			normal:=StructuralCopy(s3.faces[1][1kl]);
			if normal[1]=v then
				normal:=1;
			else
				normal:=-1;
			fi;
			p:=Position(dpoints, 1kl);
			Add(coorient.0, coorient.1[p] * normal);
		else
			Add(cofaces.0,Concatenation(enter_exit[s]));
			if IsEmpty(enter_exit[s][1]) then
				Add(coorient.0,1);
			else
				Add(coorient.0,-1);
			fi;
		fi;
		s:=s+1;
	od;

	data:=rec(	manifold:=rec(vertices:=s2.vertices, faces:=s2.faces),
				images:=s2.images,
				preimages:=s2.preimages,
				colines:=colines,
				cosheets:=cosheets,
				cofaces:=cofaces,
				coorient:=coorient);

return data;
end);

# ПРОВЕРКА: 
# 1) Проверка количественных характеристик для различных диаграмм по
# Twist-узлам (количество листов).
# 2) Проверялась тривиальность раскрашивающего квандла для 1-twist диаграмм
# 3) Проверялась согласованность построения троек в списках .cofaces[1]
# 4) Проверялось количество входных и выходных линий в тройных точках
# 5) #NOTE: Ориентация троных точек не проверена
