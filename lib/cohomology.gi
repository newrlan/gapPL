#-------------------------------------------------------------------------------
# редукция пары матриц A и B, A*B = 0
# Матрица A будет приведена в нормальную форму Смита. Соответствующие изменения
# будут отображены на базисах Ea и Eb, и на матрице B. Результатом будет словарь
# со следующими полями:
# dimensions    размерность пространства (= dim Ea)
# torsion       коэффициенты кручений группы (ко)гомологий
# elements 
# indices
#   elements[i] и indices[i] описывают как выражается i-ый базисный элементы
#   когомологии через базис пространства. indices[i] указывает базисные элементы
#   пространства образуют элемент когомологии, elements[i] указывает с какими
#   коэффициентами эти базисные элементы пространства должны быть включены.
# 
# i-ый базисный вектор когомологии имеет кручение torsion[i]. Если torsion[i] не
# существует, а базисный элемент есть, то его кручение равно нулю.

InstallGlobalFunction( hom_pair_reduction,
function(A, Ea, B, Eb)
    local hom, ij, n, res, x, l, s;

    hom := rec( dimensions := [0, Ea.dimensions[2]]
              , indices    := []
              , elements   := []
              , torsion    := []
              ) ;

    ij := sm_abs_argmin(A);
    n := 0;
    while not ij = fail do
        res := sm_iter_SNF(A, ij);
        #.......................................................................
        # Производим преобразование базиса Ea в соответсвии с изменениями в
        # матрице A
        while not IsEmpty(res.string) do
            x := Remove(res.string);
            l := Length(x);
            if l = 2 then
                sm_column_permut(Ea, x[2], x[1]);
            elif l = 3 then
                sm_column_add(Ea, x[2], x[1], x[3]);
            else
                Print("Something wrong in res.string at function hom_pair_reduction.\n");
            fi;
        od;

        # Убираем знак главного элемента (ести требуется)
        if A.elements[1][1] < 0 then
            sm_string_prod(A, 1, -1);
            sm_string_prod(Ea, 1, -1);
        fi;

        # Если "имеет смысл" сохраняем кручение и вектор из базиса
        if not IsOne(A.elements[1][1]) then
            Add(hom.torsion, A.elements[1][1]);
            Add(hom.indices, Remove(Ea.indices, 1));
            Add(hom.elements, Remove(Ea.elements, 1));
            # TODO тут можно проверить числа Бети.
        else
            Remove(Ea.indices, 1);
            Remove(Ea.elements, 1);
        fi;
        Ea.dimensions[1] := Ea.dimensions[1] - 1;

        Remove(A.indices, 1);
        Remove(A.elements, 1);
        A.indices := A.indices - 1;
        A.dimensions := A.dimensions - 1;

        ## TODO надо написать удаление пустых строк, чтобы не возиться с ними.
        ## Удаляем нулевые строки из матрицы A
        #s := 1;
        #while s < A.dimensions[1] + 1 do
        #    if IsEmpty(A.indices[s]) then
        #        Remove(A.indices, s);
        #        Remove(A.elements, s);
        #        A.dimensions[1] := A.dimensions[1] - 1;

        #        Remove(Ea.indices, s);
        #        Remove(Ea.elements, s);
        #        Ea.dimensions[1] := Ea.dimensions[1] - 1;
        #    else
        #        s := s + 1;
        #    fi;
        #od;
        ## TODO нехватает учета в hom и torsion

        #.......................................................................
        # Производим преобразования базиса Eb в соответствии с изменениями в
        # матрице A.
        while not IsEmpty(res.column) do
            x := Remove(res.column, 1);
            l := Length(x);
            if l = 2 then
                sm_string_permut(B, x[1], x[2]);
                sm_string_permut(Eb, x[1], x[2]);
            elif l = 3 then
                sm_string_add(B, x[2], x[1], -x[3]);
                sm_string_add(Eb, x[2], x[1], x[3]);
            else
                Print("Something wrong in res.column at function hom_pair_reduction.\n");
            fi;
        od;

        if IsEmpty(Remove(B.indices, 1)) then
            Remove(B.elements, 1);
            B.dimensions[1] := B.dimensions[1] - 1;

            Remove(Eb.elements, 1);
            Remove(Eb.indices, 1);
            Eb.dimensions[1] := Eb.dimensions[1] - 1;
        else # ошибка: всегда должна получаться нулевая строка
            Print("First string of B isn't zero!\n");
            n := n.stop; # TODO: сделать более аккуратную остановку
        fi;

        ij := sm_abs_argmin(A);

        n := n + 1;
    od;

    Append(hom.elements, Ea.elements);
    Append(hom.indices, Ea.indices);
    hom.dimensions[1] := hom.dimensions[1] + Ea.dimensions[1];

    return hom;
end);

#-------------------------------------------------------------------------------
# Вычисляются группы (ко)гомологий.
# Результатом будет канонически упорядоченный триангулированный политоп и
# структура cohom по каждой размерности. cohom.(n) имеет поля: dimensions,
# torsion, indices, elements --- которые подробно описаны в функции
# hom_pair_reduction.

InstallGlobalFunction( IntCoHomology,
function(pol0)
    local pol, n, l, elem, A, elements, indices, Ea, B, Eb;

    pol := PolTriangulate(pol0);
    pol := PolCanonicalOrder(pol);
    pol.cohom := rec();
    n := Length(pol.faces);     # размерность матрицы
    l := LengthPol(pol);

    elem := List([n+1, n .. 1], x -> (-1)^x);
    A := rec( dimensions := [l.(n), l.(n-1)]
            , elements := List(pol.faces[n], x -> StructuralCopy(elem))
            , indices := StructuralCopy(pol.faces[n])
            ) ;
    Ea := IdenticalMatSM(l.(n));

    while n > 1 do
        n := n - 1;
        Remove(elem, 1);
        B := rec( dimensions := [l.(n), l.(n-1)]
                , elements := List(pol.faces[n], x -> StructuralCopy(elem))
                , indices := StructuralCopy(pol.faces[n])
                ) ;
        Eb := IdenticalMatSM(l.(n));

        # Проверка, XXX этот участок надо будет удалить
        if not IsZeroSM(SMprod(A, B)) then Print("Mistake. A * B not zero!\n"); fi;

        pol.cohom.(n+1) := hom_pair_reduction(A, Ea, B, Eb);

        A := StructuralCopy(B);
        Ea := StructuralCopy(Eb);
    od;

    Eb := IdenticalMatSM(A.dimensions[2]);
    B := rec( dimensions := [A.dimensions[2], 1]
            , elements := List([1 .. A.dimensions[2]], x -> [])
            , indices := List([1 .. A.dimensions[2]], x -> [])
            ) ;
    pol.cohom.1 := hom_pair_reduction(A, Ea, B, Eb);
    pol.cohom.0 := Eb;
    pol.cohom.0.torsion := [];

    for n in [0 .. Length(pol.faces)] do
        pol.cohom.(n).beti := pol.cohom.(n).dimensions[1];
        Unbind(pol.cohom.(n).dimensions);
    od;

    return pol;
end);
