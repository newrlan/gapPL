# Формат данных RatFunc.
# Полином <M>f</M> можно разложить на множители <M>f = a f_1^{k_1} \dots
# f_n^{k_n}</M>. Это дает возможность задать функцию <M>f</M> в электронном виде
# более экономично. А именно функция будет записана так
# [a, f_1, k_1, ... , f_n, k_n].
# Как видно первый элемент списка является коэффициентом функции при старшей
# степени, затем каждый четный элемент списка объявляет множители функции, а
# каждый нечетный кратности множителей стоящих перед ним. Такую запись полинома
# будем называть форматом RatFunc.

# В формате RatFunc рациональная функция <M>r</M> представляется в виде
# двухэлементного списка. Первый и второй элементы этого списка это числитель и
# знаменатель рациональной функции <M>r</M>. Если вместо числителя(знаменателя) стоит
# пустое множество, считаем, что числитель(знаменатель) равен единице.

#------------------------------------------------------------------------------
#			<ManSection><Func Name="ConvertPolynomeToRatFunc" Arg="f" />
#				<Description>
#					преобразует полином <M>f</M> в формат RatFunc
#				</Description>
#			</ManSection>

InstallGlobalFunction( ConvertPolynomeToRatFunc,function(f0)
	local coef, list, f;

	if IsPolynomial(f0) then
		coef:=LeadingCoefficient(f0);
		if IsZero(coef) then
			list:=[0];
		else
			f:=f0/coef;
			if IsOne(f) then
				f:=[];
			else
				f:=Factors(f);
			fi;
			list:=Collected(f);
			list:=Concatenation(list);
			Add(list,coef,1);
		fi;
	else
		list:=[f0];
	fi;

	return list;
end);

#------------------------------------------------------------------------------
#			<ManSection><Func Name="ConvertToRatFunc" Arg="f" />
#				<Description>
#					преобразует рациональную функцию <M>f</M> в формат RatFunc.
#				</Description>
#			</ManSection>

InstallGlobalFunction( ConvertToRatFunc,function(f)
	local numerator, denominator;

	numerator:=NumeratorOfRationalFunction(f);
	numerator:=ConvertPolynomeToRatFunc(numerator);
	denominator:=DenominatorOfRationalFunction(f);
	denominator:=ConvertPolynomeToRatFunc(denominator);

	return [numerator, denominator];
end);

#------------------------------------------------------------------------------
#			<ManSection><Func Name="GcdPolynomial" Arg="f,g" />
#				<Description>
#					вычисляется наибольший общий делитель полиномов <M>f</M> и
#					<M>g</M>. Входящие данные могут быть как в формате RatFunc,
#					так и в виде полиномов, ответ будет создан только в формате
#					RatFunc.
#				</Description>
#				<Description>
#					В вычислении наибольшего общего делителя не участвуют
#					коэффициенты фукнций <M>f</M> и <M>g</M>, их следует
#					вычислять отдельно.
#				</Description>
#			</ManSection>

InstallGlobalFunction( GcdPolynomial,function(f0,g0)
	local f, g, gcd, gk, deg, pos;

	if IsPolynomial(f0) then
		f:=ConvertPolynomeToRatFunc(f0);
	else
		f:=StructuralCopy(f0);
	fi;
	if IsPolynomial(g0) then
		g:=ConvertPolynomeToRatFunc(g0);
	else
		g:=StructuralCopy(g0);
	fi;

	gcd:=[1];
	Remove(g,1);
	Remove(f,1);
	while (not IsEmpty(g)) and (not IsEmpty(f)) do
		gk:=Remove(g,1);
		deg:=Remove(g,1);
		pos:=Position(f, gk);
		if pos = fail then
		else
			Remove(f, pos);
			Add(gcd, gk);
			Add(gcd, Minimum(deg, Remove(f, pos)));
		fi;
	od;

	return gcd;
end);

#------------------------------------------------------------------------------
#			<ManSection><Func Name="LcmPolynomial" Arg="f,g" />
#				<Description>
#					вычисляется наименьшее общее кратное полиномов <M>f</M> и
#					<M>g</M>. Входящие данные могут быть как в формате RatFunc,
#					так и в виде полиномов, ответ будет создан в формате
#					RatFunc.
#				</Description>
#				<Description>
#					В вычислении наименьшего общего кратного не участвуют
#					коэффициенты фукнций <M>f</M> и <M>g</M>, их следует
#					вычислять отдельно.
#				</Description>
#			</ManSection>

InstallGlobalFunction( LcmPolynomial,function(f0,g0)
	local f, g, lcm, gk, deg, pos;

	if IsPolynomial(f0) then
		f:=ConvertPolynomeToRatFunc(f0);
	else
		f:=StructuralCopy(f0);
	fi;
	if IsPolynomial(g0) then
		g:=ConvertPolynomeToRatFunc(g0);
	else
		g:=StructuralCopy(g0);
	fi;

	lcm:=[1];
	Remove(g,1);
	Remove(f,1);
	while (not IsEmpty(g)) and (not IsEmpty(f)) do
		gk:=Remove(g,1);
		deg:=Remove(g,1);
		pos:=Position(f, gk);
		if pos = fail then
		else
			Remove(f, pos);
			Add(lcm, gk);
			Add(lcm, Maximum(deg, Remove(f, pos)));
		fi;
	od;
	Append(lcm, g);
	Append(lcm, f);

	return lcm;
end);

#------------------------------------------------------------------------------
#			<ManSection><Func Name="SimplifyRatFunc" Arg="f" />
#				<Description>
#					упрощение рациональной функции.
#					<Example>
# gap> f:=(x^2-y^2)/(x+y)^2;;
# gap> f1:=ConvertToRatFunc(f);
# [ [ 1, x-y, 1, x+y, 1 ], [ 1, x+y, 2 ] ]
# gap> SimplifyRatFunc(f1);
# [ [ 1, x-y, 1 ], [ 1, x+y, 1 ] ]
# gap> g:=[[2,x-y^2,0],[]];;
# gap> SimplifyRatFunc(g);
# [ [ 2 ], [ 1 ] ]
#					</Example>
#				</Description>
#			</ManSection>

InstallGlobalFunction( SimplifyRatFunc,function(f0)
	local f, l, s, i, j, nechet, chet, gameover, ind, g;

	if IsList(f0) then
		f:=StructuralCopy(f0);
	else
		f:=ConvertToRatFunc(f0);
	fi;
	if IsEmpty(f[1]) then f[1]:=[1]; fi;
	if IsEmpty(f[2]) then f[2]:=[1]; fi;

	# 1) проверяем наличие одинаковых множителей, но записанных по отдельности
	for g in f do
		l:=Length(g); # по построению это нечетное число
		s:=2;
		for i in [l-1, l-3 .. s+2] do
			if g[i] = g[s] then
				Remove(g,i);
				g[s+1]:=g[s+1] + g[i];
				Remove(g,i);
				l:=l-2;
			fi;
		od;
		s:=s+2;
	od;

	# 2) проводим сокращение числителя и знаменателя
	# сейчас в l записано значение Length(f[2])
	i:=2;
	while i < l do
		j:=Position(f[1], f[2][i]);
		if j = fail then
		else
			if f[2][i+1] > f[1][j+1] then # остается в знаменателе
				Remove(f[1], j);
				f[2][i+1]:=f[2][i+1] - Remove(f[1],j);
			elif f[2][i+1] < f[1][j+1] then # остается в числителе
				Remove(f[2], i);
				f[1][j+1]:=f[1][j+1] - Remove(f[2],i);
				l:=l-2;
			else
				Remove(f[2],i);
				Remove(f[2],i);
				Remove(f[1],j);
				Remove(f[1],j);
				l:=l-2;
			fi;
		fi;
		i:=i+2;
	od;

	# 3) ищем нули и единицы
	for g in f do
		l:=Length(g);
		nechet:=List([1..(l-1)/2], x -> 2*x+1);
		for i in nechet do
			if IsZero(g[i]) then
				if not IsZero(g[i-1]) then
					g[i-1]:=1;
				else
					Print("Function f has a 0^0 as a multiplier.\n");
					Print("I don't know what I have to do.\n");
					break;
				fi;
			fi;
		od;
		chet:=nechet - 1;
		gameover:=false;
		ind:=[];
		for i in chet do
			if IsZero(g[i]) then
				g[1]:=0;
				gameover:=true;
				ind:=[];
			elif IsOne(g[i]) then
				Add(ind, i);
			fi;
		od;
		while not IsEmpty(ind) do
			i:=Remove(ind);
			Remove(g,i);
			Remove(g,i);
		od;
		if gameover then
			g:=[0];
		fi;
	od;
	if IsZero(f[1][1]) then f[1]:=[0]; fi;
	if IsZero(f[2][1]) then f[2]:=[0]; fi;

	return f;
end);

#------------------------------------------------------------------------------
#			<ManSection><Func Name="ProdRatFunc" Arg="f,g" />
#				<Description>
#					функция произведения двух рациональных функций формата
#					RatFunc
#				</Description>
#			</ManSection>

InstallGlobalFunction( ProdRatFunc,function(f0,g0)
	local f, g, i;

	f:=StructuralCopy(f0);
	g:=StructuralCopy(g0);
	for i in [1,2] do
		if f[i] = [] then f[i]:=[1]; fi;
		if g[i] = [] then g[i]:=[1]; fi;
		f[i][1]:=f[i][1]*g[i][1];
		Remove(g[i],1);
		Append(f[i], g[i]);
	od;
	f:=SimplifyRatFunc(f);

	return f;
end);

#------------------------------------------------------------------------------
#			<ManSection><Func Name="ConvertFromRatFuncToPolynom" Arg="f" />
#				<Description>
#					преобразует полином <M>f</M> из формата RatFunc в обычный формат
#				</Description>
#			</ManSection>

InstallGlobalFunction( ConvertFromRatFuncToPolynom,function(f0)
	local f, g, f1, deg;

	f:=StructuralCopy(f0);
	g:=Remove(f,1);
	while not IsEmpty(f) do
		f1:=Remove(f,1);
		deg:=Remove(f,1);
		g:=g*f1^deg;
	od;

	return g;
end);

#------------------------------------------------------------------------------
#			<ManSection><Func Name="ConvertFromRatFunc" Arg="f" />
#				<Description>
#					преобразует рациональную функцию <M>f</M> из формата RatFunc в
#					обычный формат.
#				</Description>
#			</ManSection>

InstallGlobalFunction( ConvertFromRatFunc,function(f)
	local numerator, denominator, g;

	numerator:=ConvertFromRatFuncToPolynom(f[1]);
	denominator:=ConvertFromRatFuncToPolynom(f[2]);
	g:=numerator/denominator;

	return g;
end);


#------------------------------------------------------------------------------
#			<ManSection><Func Name="SumRatFunc" Arg="f,g" />
#				<Description>
#					Суммируем рациональные функции заданные в формате RatFunc
#				</Description>
#			</ManSection>

InstallGlobalFunction( SumRatFunc,function(f,g)
	local a, b, denominator, left, right, gcd, result;

	a:=[[],f[2]];
	b:=[[],g[2]];
	denominator:=ProdRatFunc(a,b)[2];
	
	a:=[f[1],[]];
	b:=[g[2],[]];
	left:=ProdRatFunc(a,b);

	a:=[g[1],[]];
	b:=[f[2],[]];
	right:=ProdRatFunc(a,b);

	gcd:=GcdPolynomial(left[1], right[1]);
	result:=[gcd, denominator];
	gcd:=[[],gcd];
	left:=ProdRatFunc(left, gcd);
	left:=ConvertFromRatFunc(left);
	right:=ProdRatFunc(right, gcd);
	right:=ConvertFromRatFunc(right);
	g:=left + right;
	g:=ConvertPolynomeToRatFunc(g);
	g:=[g,[]];
	result:=ProdRatFunc(g, result);

	return result;
end);

#------------------------------------------------------------------------------
#			<ManSection><Func Name="DerivativePolynomRatFunc" Arg="f,x" />
#				<Description>
#					вычисляется производная полинома заданной в формате RatFunc
#				</Description>
#			</ManSection>

InstallGlobalFunction( DerivativePolynomRatFunc,function(f,x)
	local l, nechet, chet, gcd, derivat, shablon, g, d, i;

	# Можно было бы сразу перевести полином из формата RatFunc в обычный, но это
	# не продуктивно, поскольку производная может иметь часть множителей
	# изначальной функции. Внесение общих множителей осложнит задачу
	# факторизации.
	
	l:=Length(f);
	nechet:=List([1..(l-1)/2], i -> 2*i+1);
	chet:=List([1..(l-1)/2], i -> 2*i);
	gcd:=StructuralCopy(f);
	for i in nechet do
		gcd[i]:=gcd[i]-1;
	od;

	derivat:=0;
	shablon:=StructuralCopy(f);
	shablon[1]:=1;
	for i in nechet do
		shablon[i]:=1;
	od;

	for i in chet do
		g:=StructuralCopy(shablon);
		d:=Remove(g,i);
		Remove(g,i);
		d:=Derivative(d,x);
		d:=ConvertPolynomeToRatFunc(d);
		g:=ProdRatFunc([d,[]],[g,[]]);
		g:=ConvertFromRatFuncToPolynom(g[1]);
		derivat:=derivat + f[i+1]*g;
	od;
	derivat:=ConvertPolynomeToRatFunc(derivat);
	derivat:=ProdRatFunc([derivat,[]], [gcd,[]]);

	return derivat[1];
end);


#------------------------------------------------------------------------------
#			<ManSection><Func Name="DerivativeRatFunc" Arg="f,x" />
#				<Description>
#					вычисляется производная рациональной функции заданной в формате RatFunc
#				</Description>
#			</ManSection>

InstallGlobalFunction( DerivativeRatFunc,function(f,x)
	local denominator, dn, dd, left, right, numerator, l, g, i;

	denominator:=[[],f[2]];

	dn:=DerivativePolynomRatFunc(f[1],x);
	dd:=DerivativePolynomRatFunc(f[2],x);
	left:=ProdRatFunc([dn,[]], [f[2],[]]);
	right:=ProdRatFunc([dd,[]],[f[1],[]]);
	right[1][1]:= - right[1][1];
	numerator:=SumRatFunc(left, right);
	denominator:=StructuralCopy(f[2]);
	denominator[1]:=denominator[1]^2;
	l:=Length(denominator);
	for i in [1 .. (l-1)/2] do
		denominator[2*i +1] := 2*denominator[2*i+1];
	od;
	denominator:=[[],denominator];
	g:=ProdRatFunc(numerator, denominator);

	return g;
end);

#------------------------------------------------------------------------------
#			<ManSection><Func Name="JacobiMatRatFunc" Arg="listf,listx" />
#				<Description>
#					Для системы рациональных фукнций listf вычисляется матрица
#					якоби по переменным listx. Формат записи данных в список
#					listf может быть смешанным, то есть в список могут входить
#					как обычные функции, так и функции заданные в формате
#					RatFunc.
#				</Description>
#			</ManSection>

InstallGlobalFunction( JacobiMatRatFunc,function(listf,listx)
	local fk, mat, f;

	fk:=[];
	for f in listf do
		if IsRationalFunction(f) then
			Add(fk, ConvertToRatFunc(f));
		else
			Add(fk, f);
		fi;
	od;
	mat:=List(fk,f->List(listx,x->DerivativeRatFunc(f,x)));

	return mat;
end);

