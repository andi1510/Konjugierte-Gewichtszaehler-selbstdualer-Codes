///////////////////////////////////////////////////////////////////////////
// Laplace-Operator und Loesung des linearen Gleichungssystems aus
// Abschnitt 5.1.
//
// Voraussetzung:
//   Zuvor muss CCWE_Genus2.m geladen worden sein. Dadurch sind der
//   Polynomring P und die Variablen
//
//     x00,x01,x10,x11,y00,y01,y10,y11
//
//   definiert.
//
// Der Laplace-Operator ist
//
//   Delta = Sum_{v in F_2^2} d^2 / (d x_v d bar{x}_v).
//
// In unserem Polynomring entsprechen y00,...,y11 den konjugierten
// Variablen.
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// 1. Laplace-Operator
///////////////////////////////////////////////////////////////////////////

function ConjugateLaplace(f)
    // Variablen 1,...,4 sind x00,...,x11.
    // Variablen 5,...,8 sind y00,...,y11.
    return &+[
        Derivative(Derivative(f, i), i+4)
        : i in [1..4]
    ];
end function;


///////////////////////////////////////////////////////////////////////////
// 2. Hilfsfunktionen fuer Koeffizientenvektoren
///////////////////////////////////////////////////////////////////////////

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


///////////////////////////////////////////////////////////////////////////
// 3. Rang des von Polynomen erzeugten Raumes
///////////////////////////////////////////////////////////////////////////

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


///////////////////////////////////////////////////////////////////////////
// 4. LGS fuer den harmonischen Anteil im Grad (6,6)
///////////////////////////////////////////////////////////////////////////

/*
    Nach Abschnitt 5.1 wird

        h = F - a*S^6 - b*S^2*W44

    angesetzt.

    Die Bedingung Delta(h)=0 liefert

        Delta(F)
          = a*Delta(S^6)
          + b*Delta(S^2*W44).

    Durch Vergleich aller Monomkoeffizienten entsteht ein lineares
    Gleichungssystem fuer a und b.

    Magmas IsConsistent(A,w) loest Systeme in der Form

        v*A = w.

    Deshalb werden die Koeffizienten von Delta(S^6) und
    Delta(S^2*W44) als die beiden ZEILEN der Matrix A gespeichert.
    Der gesuchte Vektor v ist dann genau (a,b).
*/
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
