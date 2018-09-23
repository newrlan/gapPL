## Работа начата 28 июня 2016
# Разреженная матрица сжимается по строкам, то есть, если в строке матрицы есть
# ненулевой элемент он выкидывается. Структура разреженных матриц:
# sm := ref(
#           indices := [ [Int] ];
#           elements := [ [a] ];
#       );
# mat[i][j] = 0 если sm.indices[i] не содержит j
# mat[i][j] = sm.elements[i][p], где p позиция элемента j в списке sm.indices[i]

#------------------------------------------------------------------------------
#            <ManSection><Func Name="toSparceMatrix" Arg="mat" />
#                <Description>
#                    Первести матрицу в формат разреженной.
#                </Description>
#            </ManSection>
InstallGlobalFunction( toSM,
function(mat)
  local sm, inds, elem, i, str, e;

    sm := rec( indices := []
             , elements := []
             , dimensions := MutableCopyMat(DimensionsMat(mat))
             );
    for str in mat do
        inds := [];
        elem := [];
        i := 1;
        for e in str do
            if not IsZero(e) then
                Add(inds, i);
                Add(elem, e);
            fi;
            i := i + 1;
        od;
        Add(sm.indices, inds);
        Add(sm.elements, elem);
    od;

    return sm;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="fromSparceMatrix" Arg="sm" />
#                <Description>
#                    Создать матрицу по разреженной.
#                </Description>
#            </ManSection>
InstallGlobalFunction( fromSM,
function(sm)
  local mat, str, xs, s, inds, col, l;

    mat := [];
    str := 1;
    l := sm.dimensions[2];
    for inds in sm.indices do
        xs := List([1 .. l], x -> 0);
        s := 1;
        for col in inds do
            xs[col] := sm.elements[str][s];
            s := s + 1;
        od;
        str := str + 1;
        Add(mat, xs);
    od;

    return mat;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="SMtranspose" Arg="mat" />
#                <Description>
#                    Транспонирование разреженной матрицы.
#                </Description>
#            </ManSection>
InstallGlobalFunction( SMtranspose, 
function(mat)
	local n, tmat, elements, i, j, string, ind;

    n := List(mat.indices, x -> Maximum(x));
    n := Maximum( n );

    tmat := rec( indices := List([1 .. n], x -> []),
                 elements := List([1 .. n], x -> []) );
    i := 0;
    for string in mat.indices do
        i := i + 1;
        j := 0;
        for ind in string do
            j := j + 1;
            Add(tmat.indices[ind], i);
            Add(tmat.elements[ind], mat.elements[i][j]);
        od;
    od;

    tmat.dimensions := mat.dimensions{[2,1]};

    return tmat;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="SMprodT" Arg="mat1, mat2" />
#                <Description>
#                    Умножение разреженной матрицы mat1 на
#                    транспонированную к разреженной матрице mat2. То есть
#                    будет вычисленно <M>mat1 * mat2^T</M>.
#                </Description>
#            </ManSection>

InstallGlobalFunction( SMprodT, 
function(mat1, mat2)
	local pmat, i, j, elem, p, q, string1, string2, ind;

    pmat := rec( indices := [], elements := [] );
    i := 0;
    for string1 in mat1.indices do
        i := i + 1;
        Add(pmat.indices,[]);
        Add(pmat.elements,[]);
        j := 0;
        for string2 in mat2.indices do
            j := j + 1;
            elem := 0;
            p := 0;
            for ind in string1 do
                p := p + 1;
                if ind in string2 then 
                    q := Position(string2, ind);
                    elem := elem + mat1.elements[i][p]*mat2.elements[j][q];
                fi;
            od;
            if not IsZero(elem) then
                Add(pmat.indices[i], j);
                Add(pmat.elements[i], elem);
            fi;
        od;
    od;

    pmat.dimensions := [mat1.dimensions[1], mat2.dimensions[1]];

    return pmat;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="SMprod" Arg="mat1, mat2" />
#                <Description>
#                    Умножение двух разреженных матриц.
#                </Description>
#            </ManSection>

InstallGlobalFunction( SMprod, 
function(mat1, mat2)
	local mat3, mat0;

    mat3 := SMtranspose(mat2);
    mat0 := SMprodT(mat1,mat3);

    return mat0;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="IsZeroSM" Arg="mat" />
#                <Description>
#                    Провека является ли матрица mat нулевой.
#                </Description>
#            </ManSection>

InstallGlobalFunction( IsZeroSM, 
function(mat)
	local answer;

    if Set(mat.indices) = [[]] or Set(mat.indices) = [] then
        answer := true;
    else
        answer := false;
    fi;

    return answer;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="sm_string_prod" Arg="mat, ind, coef" />
#                <Description>
#                    Домножение строки ind на неулевой элемент coef
#                </Description>
#            </ManSection>

InstallGlobalFunction( sm_string_prod,
function(mat, ind, coef)

    mat.elements[ind] := mat.elements[ind] * coef;

    return ;
end);
#------------------------------------------------------------------------------
#            <ManSection><Func Name="sm_string_add" Arg="mat, ind1, ind2, coef" />
#                <Description>
#                    Сложение строк с коэффициентом. К строке ind1 прибавляется
#                    строка ind2 умноженная на coeff
#                </Description>
#            </ManSection>

InstallGlobalFunction( sm_string_add,
function(mat, ind1, ind2, coef)
    local xs, ys, x_ind, y_ind, inds, elem, x, y, a;

    xs := StructuralCopy(mat.indices[ind1]);
    ys := StructuralCopy(mat.indices[ind2]);
    x_ind := 1;
    y_ind := 1;
    inds := [];
    elem := [];
    while not (IsEmpty(xs) or IsEmpty(ys)) do
        x := xs[1];
        y := ys[1];
        if x = y then
            a := mat.elements[ind1][x_ind] + coef * mat.elements[ind2][y_ind];
            if not IsZero(a) then
                Add(elem, a);
                Add(inds, x);
            fi;
            Remove(xs, 1);
            Remove(ys, 1);
            x_ind := x_ind + 1;
            y_ind := y_ind + 1;
        elif x < y then
            Add(inds, x);
            Add(elem, mat.elements[ind1][x_ind]);
            Remove(xs, 1);
            x_ind := x_ind + 1;
        else
            Add(inds, y);
            Add(elem, coef * mat.elements[ind2][y_ind]);
            Remove(ys, 1);
            y_ind := y_ind + 1;
        fi;
    od;

    for x in xs do
        Add(inds, x);
        Add(elem, mat.elements[ind1][x_ind]);
        x_ind := x_ind + 1;
    od;
    for y in ys do
        Add(inds, y);
        Add(elem, coef * mat.elements[ind2][y_ind]);
        y_ind := y_ind + 1;
    od;

    mat.indices[ind1] := inds;
    mat.elements[ind1] := elem;


    return ;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="sm_string_permut" Arg="mat, ind1, ind2" />
#                <Description>
#                    Перестановка строк.
#                </Description>
#            </ManSection>
InstallGlobalFunction( sm_string_permut,
function(mat, ind1, ind2)
  local perm;

    perm := (ind1, ind2);
    mat.indices := Permuted(mat.indices, perm);
    mat.elements := Permuted(mat.elements, perm);

    return ;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="sm_column_prod" Arg="mat, ind, coef" />
#                <Description>
#                    Домножение столбца ind на неулевой элемент coef
#                </Description>
#            </ManSection>
InstallGlobalFunction( sm_column_prod ,
function(mat, ind1, coef)
  local i, j, xs;

    i := 1;
    for xs in mat.indices do
        j := Position(xs, ind1);
        if not j = fail then
            mat.elements[i][j] := mat.elements[i][j] * coef;
        fi;
        i := i + 1;
    od;

    return;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="sm_column_permut" Arg="mat, ind1, ind2" />
#                <Description>
#                    Перестановка столбцов.
#                </Description>
#            </ManSection>
InstallGlobalFunction( sm_column_permut ,
function(mat, ind1, ind2)
  local a, b, so, i;

    if not ind1 = ind2 then
        for i in [1 .. mat.dimensions[1]] do
            a := Position(mat.indices[i], ind1);
            b := Position(mat.indices[i], ind2);
            so := false;
            if not a = fail then
                mat.indices[i][a] := ind2;
                so := true;
            fi;
            if not b = fail then
                mat.indices[i][b] := ind1;
                so := true;
            fi;
            if so then
                SortParallel(mat.indices[i], mat.elements[i]);
            fi;
        od;
    fi;

    return;
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="sm_column_add" Arg="mat, ind1, ind2, coef" />
#                <Description>
#                    Сложение столбцов с коэффициентом. К столбцу ind1 прибавляется
#                    столбец ind2 умноженный на coef.
#                </Description>
#            </ManSection>
InstallGlobalFunction( sm_column_add ,
function(mat, ind1, ind2, coef)
    local a, b, c, i;

    if not ind1 = ind2 then

        for i in [1 .. mat.dimensions[1]] do
            a := Position(mat.indices[i], ind1);
            b := Position(mat.indices[i], ind2);
            if not b = fail then
                c := coef * mat.elements[i][b];
                if not a = fail then
                    c := c + mat.elements[i][a];
                    if IsZero(c) then
                        Remove(mat.elements[i],a);
                        Remove(mat.indices[i],a);
                    else
                        mat.elements[i][a] := c;
                    fi;
                    #mat.elements[i][a] := mat.elements[i][a] + c;
                else;
                    while b > 1 and mat.indices[i][b-1] > a do
                        b := b - 1;
                    od;
                    Add(mat.elements[i], c, b);
                    Add(mat.indices[i], ind1, b);
                fi;
            fi;
        od;

    fi;

    return ;
end);

#-------------------------------------------------------------------------------
# Пробую написать первую итерацию реализации алгоритма приведения к нормальной
# форме.


#------------------------------------------------------------------------------
#            <ManSection><Func Name="sm_abs_argmin" Arg="M" />
#                <Description>
#                    Адрес минимального по модулю не нулевого элемента в
#                    матирце M.
#                </Description>
#            </ManSection>
InstallGlobalFunction( sm_abs_argmin ,
function(M)
    local m, ij, i, j, me, es, e;

    i := 1;
    while i <= M.dimensions[1] and IsEmpty(M.elements[i]) do
        i := i + 1;
    od;
    if i <= M.dimensions[1] then
        m := AbsInt( M.elements[i][1] );
        ij := [i,1];
        for es in M.elements{[i .. M.dimensions[1]]} do
            j := 1;
            for e in es do
                me := AbsInt(e);
                if me < m then
                    ij := [i, j];
                    m := StructuralCopy(me);
                fi;
                j := j + 1;
            od;
            if IsOne(m) then break; fi;
            i := i + 1;
        od;
    else
        ij := fail;
    fi;

    return ij;
end);


#------------------------------------------------------------------------------
# Description: sm_inter_SNF
# Проводится одна итерация алгоритма приведения матрицы к нормальной форме
# Смита. Элемент ij будет перенесен в ячейку 11, все остальные элементы в первой
# строке и первом столбце будут обнулены. В течении итерации элемент в ячейке 11
# может быть изменен на меньший по модулю.

# Arguments:
# M     разреженная матирца
# ij    адрес элемента который хотим поставить в ячейку 11

# Results: В Н И М А Н И Е : Матрица M будет изменена. 
# Возврачает словарь с двумя полями:
#   string: список действий над строками матрицы M, которые были проведены.
# Если string содержит пару [a,b], то это перестановка строк a и b. Если
# список string содержит триплет [a,b,c], то это сложение к строке a строки b
# домноженной на число c.
#   column: список действий над столбцами матрицы M, которые были проведены.
# Если column содержит пару [a,b], то это перестановка столбцов a и b. Если
# список column содержит триплет [a,b,c], то это сложение к столбцу a столбца b
# домноженного на число c.

InstallGlobalFunction( sm_iter_SNF ,
function(M, ij)
    local column, string, i, j, rep, str, s, r, d, l, xs;

    column := [];
    string := [];

    # поиск минимального элемента ..............................................
        i := ij[1];
        j := M.indices[i][ij[2]];
        if not IsOne(i) then 
            sm_string_permut(M, 1, i);
            Add(string, [1, i]);
        fi;
        if not IsOne(j) then
            sm_column_permut(M, 1, j); 
            Add(column, [1, j]);
        fi;

    # редуцирование столбца\строки .............................................
    #   Логика редуцирования следущаня:
    #   1) редуцировать полностью первый столбец
    #   2) редуцировать ненулевой элемент в строке (конечно, исключая элемент 11)
    #   3) Если элемент не редуцирован, переставить столбец на первое место
    #   4) Возвращаться к шагу 1, до тех пор пока первая строк не обнулится.

    # По построению сейчас элемент M_11 самый минимальный по абсолютному
    # значению во всей матрице.
    rep := true;
    while rep do
        str := [];
        s := 1;
        for xs in M.indices do
            if (not IsEmpty(xs)) and IsOne(xs[1]) then
                Add(str, s);
            fi;
            s := s + 1;
        od;
        Remove(str, 1);
        while not IsEmpty(str) do
            s := str[1];
            r := M.elements[s][1] mod M.elements[1][1];
            d := (M.elements[s][1] - r) / M.elements[1][1];
            sm_string_add(M, s, 1, -d);
            Add(string, [s, 1, -d]);
            if IsZero(r) then
                Remove(str, 1);
            else
                sm_string_permut(M, 1, s);
                Add(string, [1,s]);
            fi;
        od;

        rep := true;
        l := Length(M.indices[1]);
        while rep and (l > 1) do
            s := M.indices[1][2];
            r := M.elements[1][2] mod M.elements[1][1];
            d := (M.elements[1][2] - r) / M.elements[1][1];
            sm_column_add(M, s, 1, -d);
            Add(column, [s, 1, -d]);
            if IsZero(r) then
                l := l - 1;
            else
                sm_column_permut(M, s, 1);
                Add(column, [s, 1]);
                rep := false;
            fi;
        od;
        rep := not rep;
    od;

    return rec(column := column, string := string);
end);

#------------------------------------------------------------------------------
#            <ManSection><Func Name="IdenticalMat" Arg="n" />
#                <Description>
#                    Создаем единичную матрицу размера <M>n\times n</M>.
#                </Description>
#            </ManSection>
InstallGlobalFunction( IdenticalMatSM,
function(n)
    local E;

    E := rec( dimensions := [n, n]
            , indices := List([1 .. n], x -> [x])
            , elements := List([1 .. n], x -> [1])
            ) ;
    return E;
end);
