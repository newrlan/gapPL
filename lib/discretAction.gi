#------------------------------------------------------------------------------
#            <ManSection><Func Name="getLsystem" Arg="pol, action" />
#                <Description>
#                   По действию на стандартном пентахоре (12345) вычисляется
#                   действие на всем многообразии <C>pol</C>. 
#                </Description>
#                <Description>
#                   Действие на стандартном пентахоре (12345) задается во
#                   входном параметре <C>action</C> и представляет из себя
#                   список троек индексов тетраэдров <M>1:1234, 2:1235, 3:1245,
#                   4:1345, 5:2345</M>.
#                </Description>
#                </Description>
#                   Результат будет содержать поля <C>dim</C> --- количество
#                   переменных в многочлене действия, <C>monoms</C> --- мономы
#                   многочлена действия, <C>coeffs</C> --- коэффициенты каждого
#                   из монома, <C>rough_invariant</C> --- грубый инвариант.
#                   <Example>
#                   <!-- TODO создать пример на 4-торе -->
#                   </Example>
#                <Description>
#            </ManSection>

InstallGlobalFunction(getLsystem, function(pol1, func)
  local pol, characteristic, field, l, st3, matL, mat, c3, e, bord, ind, s, V, W, f, L, BL, num_coloring, matBL, rebra_mat, r, subind, i, RV, U, BU, lu, monoms, coeffs, j, res, m, xs, pos, monomes, rough_invariant, pair, c4, t, monom, c;

    if not Length(pol1.faces) = 4 then
        Print("Dimenstion of the polytope isn't equal to 4.\n");
    fi;
    pol := PolTriangulate(pol1);
    pol := PolCanonicalOrder(pol);

    #--- Размерность ядра и матрица отображения из ядра в прообраз --------------
    if not IsEmpty(func.coeffs) then
        characteristic := Characteristic(func.coeffs[1]);
        field := GF(characteristic);
    else
        field := Rationals;
    fi;
    l := LengthPol(pol);
    st3 := List([1 .. l.3], x -> StarFace(pol, [3, x]).4);
    matL := [ [ 0, 0,-1, 1, 0]
            , [ 2,-1,-2, 3,-1]
            , [ 1, 0,-2, 2,-1]
            , [-1, 1, 0,-1, 0]
            , [-2, 1, 1,-2, 0]
            ];
    matL := matL * One(field);
    mat := NullMat(l.3, l.3);
    c3 := 0;
    for pair in st3 do
        c3 := c3 + 1;
        e := 1;
        for c4 in pair do
            bord := pol.faces[4][c4];
            ind := Position(bord, c3);
            s := 1;
            for i in bord do
                mat[c3][i] := mat[c3][i] + e * matL[ind][s];
                s := s + 1;
            od;
            e := - 1;
        od;
    od;
    mat := mat * One(field);
    mat := TransposedMat(mat);

    V := field^l.3;
    W := field^l.3;
    f := LeftModuleHomomorphismByMatrix(Basis(V), mat, Basis(W));
    L := KernelOfAdditiveGeneralMapping(f);

    BL := Basis(L);
    num_coloring := Length(BL) - l.4 / 2 - 2 * l.0;
    # Print(Length(BL) - l.4 / 2 - 2 * l.0, "\tValue of colorings.\n");
    matBL := List(BL, x -> x);
    matBL := List(matBL, x -> List(x, y -> y));

    #--- Факторизация пространства L -------------------------------------------

    # Построение векторов по ребрам
    rebra_mat := NullMat(l.1, l.3);
    for t in [1 .. l.3] do
        r := FaceComp(pol, [3,t]).1;
        rebra_mat[r[1]][t] :=-1;
        rebra_mat[r[2]][t] := 2;
        rebra_mat[r[3]][t] :=-1;
        rebra_mat[r[4]][t] :=-1;
        rebra_mat[r[5]][t] := 0;
        rebra_mat[r[6]][t] := 1;
    od;
    rebra_mat := rebra_mat * One(field);
    rebra_mat := TriangulizedMat(rebra_mat);
    subind := [];
    i := Length(rebra_mat);
    while IsZero(rebra_mat[i]) do
        Remove(rebra_mat);
        i := i - 1;
    od;             # выделяем базис в подпространстве реберных операторов
    RV := VectorSpace(field, rebra_mat);
    U := L / RV;    # факторизация пространства L по ядерному пространству ребер
    BU := Basis(U);
    # Print("Length of BU basis:\t", Length(BU), "\n");
    lu := NaturalHomomorphismBySubspace(L, RV);
    mat := List(BU, vect -> PreImagesRepresentative(lu, vect));

    monoms := [];   # мономы составленные по триангуляции
    coeffs := [];   # коэффициенты соответствующих мономов из списка 
    # Действие на каждом пентахоре домножается на согласованную ориентацию
    # пентахора внутри триангуляции. значит ли это, что частоты для k и -k
    # должны совпадать?
    if characteristic = 2 then
        e := List(pol.faces[4], x -> 1);
    else;
        e := PolOrient(pol);
    fi;
    j := 0; # индекс пентахора
    for c4 in pol.faces[4] do
        j := j + 1;
        i := 0;
        for monom in func.monoms do
            i := i + 1;
            Add(monoms, c4{monom});
            Add(coeffs, func.coeffs[i] * e[j]);
        od;
    od;
    coeffs := coeffs * One(field);

    # WARNING: Надо проверить, что для вычислений достаточно указать только
    # характеристику поля над которым работаем. Дело в том, что в коде часто
    # возникают умножения. Характеристики может оказаться просто не достаточно.

    #--- Перевод триплетов в ядро ----------------------------------------------
    mat := TransposedMatMutable(mat);
    res := [];
    res := rec(monoms := [], coeffs := []);
    if not IsEmpty(mat) then
        m := 0; # индекс монома
        for monom in monoms do
            m := m + 1;
            xs := linearReplacement(mat{monom}, coeffs[m]);
            # По построению сейчас xs = rec(monosm, coeffs) причем поле coeffs
            # не содержит ни одного нуля.
            i := 1;
            for c in xs.coeffs do
                pos := Position(res.monoms, xs.monoms[i]);
                if not pos = fail then
                    res.coeffs[pos] := res.coeffs[pos] + c;
                    if IsZero(res.coeffs[pos]) then
                        Remove(res.coeffs, pos);
                        Remove(res.monoms, pos);
                    fi;
                else;
                    Add(res.monoms, xs.monoms[i]);
                    Add(res.coeffs, c);
                fi;
                i := i + 1;
            od;
        od;
    fi;

    # Проводим переиндексацию оставшихся перменных после сокращений
    ind := Set(Concatenation(res.monoms));
    res.monoms := List(res.monoms, abc -> List(abc, x -> Position(ind, x)));

    return rec( dim := Length(ind), 
                monomes := res.monoms,
                coeffs := res.coeffs, 
                rough_invariant := num_coloring );
end);

#-------------------------------------------------------------------------------
# Вычисляются чистоты по многочлену дискретного дейсвия над указанным полем
# Галуа FG(2^k). На вход функции поступает дискретное действие (в формате
# выдаваемом функцией getLsystem_char2_deg3) и степень двойки для поля Галуа.
# Результатом будет список частот каждого элемента поля (порядок следования
# элементов задается <C>List(field, x -> x)</C>.

InstallGlobalFunction(lookover, function(func, k)
    local vol, field, elems, frec, action, m, pos, vect, abc;

    if not IsEmpty(func.coeffs) then
        vol := (Characteristic(func.coeffs[1]))^k;
        field := GF(vol);
        elems := List(field, x -> x);
        frec := List(elems, x -> 0);
    else
        vol := 0;
        elems := [];
    fi;
    if not IsEmpty(func.monomes) then
        for vect in IteratorOfTuples(elems, func.dim) do
            action := 0;
            m := 0;
            for abc in func.monomes do
                m := m + 1;
                action := action + func.coeffs[m] * Product(vect{abc});
            od;
            pos := Position(elems, action);
            frec[pos] := frec[pos] + 1;
        od;
    else
        frec := [];
    fi;

    return frec/Sum(frec);
end);

#-------------------------------------------------------------------------------
# Структура данных содержащая нетривиальные действия для различных
# характеристик. Все действия представляют из себя однородные многочлены над
# некоторыми конечными полями. Список <C>actionlib.char.deg</C> содержит
# некоторые примеры действий степени <C>deg</C> для характеристики <C>char</C>.
InstallValue(actionlib, rec(
    2 := rec(
         3 := [ rec( coeffs := [1,1,1,1,1] * Z(2)
                   , monoms := [ [1,2,4], [1,3,4], [1,3,5], [2,3,5], [2,4,5] ]
                   )
              , rec( coeffs := [1,1,1,1,1,1,1,1,1] * Z(2)
                   , monoms := [ [4,5,5], [3,5,5], [1,5,5], [3,3,4], [2,2,4]
                               , [3,3,3], [1,3,3], [2,2,3], [1,2,2] ]
                   )
              ]
            )
));
