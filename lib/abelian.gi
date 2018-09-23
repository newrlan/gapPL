#------------------------------------------------------------------------------
#            <ManSection><Func Name="linearReplacement" Arg="mat, c" />
#                <Description>
#                   Положим имеется некоторый моном над коммутативным кольцом.
#                   Переменные с - коэффициент монома, mat - матрица линейной
#                   замены для монома (количество строк равно степени монома, в
#                   строке i указывается линейная замена для i-го множителя
#                   монома). Функция возвращает мономы в новых переменных и
#                   список их коэффициентов (мономы с нулевыми коэффициентами не
#                   включаются).
#                </Description>
#            </ManSection>

InstallGlobalFunction(linearReplacement, function(mat, c)
  local res, dim, coef, pos, monom;

    # TODO: Хороший прирост производительности на этом шаге можно будет
    # получить, если внутри функции перевести матрицу mat в формат разреженных
    # матриц.

    res := rec(monoms := [], coeffs := []);
    dim := DimensionsMat(mat);
    # dim[1]    отвечает за степень монома
    # dim[2]    количество новых переменных

    # Для экономии памяти сохраняем только те мономы которые имеют ненулевой
    # коэффициент.
    for monom in IteratorOfTuples([1 .. dim[2]], dim[1]) do
        coef := List([1 .. dim[1]], x -> mat[x][monom[x]]);
        coef := c * Product(coef);
        # можно проще
        # coef := c * Product(mat{[1 .. dim[2]]}{monom});
        if not IsZero(coef) then
            Sort(monom);
            pos := Position(res.monoms, monom);
            if not pos = fail then
                res.coeffs[pos] := res.coeffs[pos] + coef;
                if IsZero(res.coeffs[pos]) then
                    Remove(res.monoms, pos);
                    Remove(res.coeffs, pos);
                fi;
            else;
                Add(res.monoms, monom);
                Add(res.coeffs, coef);
            fi;
        fi;
    od;
    # TODO: Интересно можно ли именами словаря сделать списки. Идея ---
    # запихнуть сами мономы в эти имена, а значениями сделать коэффициенты. Если
    # коэффициенты получаются нулевыми то удалять. Если это окажется эффективной
    # идеей, то лучше пересобрать все именно так. Потенциальное преимущество это
    # время обращения к мономам.

    return res;
end);
