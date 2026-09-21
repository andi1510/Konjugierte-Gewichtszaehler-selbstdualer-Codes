// Laplace-Operator und harmonische Projektion zu Abschnitt 5.1.
// Vor der Rechnung CCWE_Genus2.m laden.


// Laplace-Operator

function ConjugateLaplace(f)
    // Variablen 1,...,4 sind x00,...,x11.
    // Variablen 5,...,8 sind y00,...,y11.
    return &+[
        Derivative(Derivative(f, i), i+4)
        : i in [1..4]
    ];
end function;


// Hilfsfunktionen fuer Koeffizientenvektoren

function AllMonomials(polynomials)
    monomials := [];

    for f in polynomials do
        for m in Monomials(f) do
            if m notin monomials then
                Append(~monomials, m);
            end if;
        end for;
    end for;

    return monomials;
end function;


function PolynomialCoefficientVector(f, monomials)
    Q := BaseRing(Parent(f));

    return Vector(
        Q,
        [ MonomialCoefficient(f, m) : m in monomials ]
    );
end function;


// Rang des von Polynomen erzeugten Raumes

function PolynomialSpanRank(polynomials)
    if #polynomials eq 0 then
        return 0;
    end if;

    Q := BaseRing(Parent(polynomials[1]));
    monomials := AllMonomials(polynomials);

    if #monomials eq 0 then
        return 0;
    end if;

    entries := &cat[
        [ MonomialCoefficient(f, m) : m in monomials ]
        : f in polynomials
    ];

    M := Matrix(Q, #polynomials, #monomials, entries);

    return Rank(M);
end function;


// LGS fuer den harmonischen Anteil im Grad (6,6)

// Ansatz: h = F - a*S^6 - b*S^2*W44 mit Delta(h) = 0.
// Die Koeffizienten von Delta(S^6) und Delta(S^2*W44)
// bilden die beiden Zeilen der Koeffizientenmatrix.
// Gesucht ist (a,b) mit (a,b)*coefficientMatrix = rightSide.

function HarmonicProjection66(F, S, W44)
    Q := BaseRing(Parent(F));

    deltaF := ConjugateLaplace(F);
    deltaS6 := ConjugateLaplace(S^6);
    deltaS2W := ConjugateLaplace(S^2 * W44);

    monomials := AllMonomials([
        deltaF,
        deltaS6,
        deltaS2W
    ]);

    coefficientMatrix := Matrix(
        Q,
        2,
        #monomials,
        [ MonomialCoefficient(deltaS6, m) : m in monomials ]
        cat
        [ MonomialCoefficient(deltaS2W, m) : m in monomials ]
    );

    rightSide := PolynomialCoefficientVector(deltaF, monomials);

    isConsistent, solution, nullspace :=
        IsConsistent(coefficientMatrix, rightSide);

    if not isConsistent then
        error "Das lineare Gleichungssystem fuer a und b ist nicht loesbar.";
    end if;

    // Rang 2 bedeutet, dass a und b eindeutig bestimmt sind.
    if Rank(coefficientMatrix) ne 2 then
        error "Die Koeffizienten a und b sind nicht eindeutig bestimmt.";
    end if;

    a := solution[1];
    b := solution[2];

    h := F - a*S^6 - b*S^2*W44;

    if ConjugateLaplace(h) ne 0 then
        error "Interner Fehler: Der berechnete Anteil ist nicht harmonisch.";
    end if;

    return h, a, b;
end function;
