// Ausgabe der Rechnung aus Abschnitt 5.1

procedure PrintHarmonicCalculation(F, number, S, W44)

    Q := Rationals();

    printf "\n";
    printf "============================================================\n";
    printf "F_%o\n", number;
    printf "============================================================\n\n";


    // Ausgangspolynom

    printf "Polynom F_%o:\n\n", number;
    printf "%o\n\n", F;


    // Laplace-Operator

    deltaF := ConjugateLaplace(F);
    deltaS6 := ConjugateLaplace(S^6);
    deltaS2W44 := ConjugateLaplace(S^2 * W44);

    printf "Laplace-Operator auf F_%o:\n\n", number;
    printf "Delta(F_%o) =\n", number;
    printf "%o\n\n", deltaF;

    printf "Delta(S^6) =\n";
    printf "%o\n\n", deltaS6;

    printf "Delta(S^2 * W_(4,4)) =\n";
    printf "%o\n\n", deltaS2W44;


    // Koeffizientenvergleich fuer die harmonische Projektion.

    monomials := AllMonomials([
        deltaF,
        deltaS6,
        deltaS2W44
    ]);

    A := Matrix(
        Q,
        2,
        #monomials,
        [ MonomialCoefficient(deltaS6, m) : m in monomials ]
        cat
        [ MonomialCoefficient(deltaS2W44, m) : m in monomials ]
    );

    rhs := PolynomialCoefficientVector(deltaF, monomials);


    printf "Lineares Gleichungssystem:\n\n";

    printf "Gesucht: (a_%o, b_%o)\n\n", number, number;

    printf "(a_%o, b_%o) * A_%o = r_%o\n\n",
        number, number, number, number;

    printf "Koeffizientenmatrix A_%o:\n\n", number;
    printf "%o\n\n", A;

    printf "Rechte Seite r_%o:\n\n", number;
    printf "%o\n\n", rhs;


    // LGS loesen

    isConsistent, solution, nullspace :=
        IsConsistent(A, rhs);

    printf "LGS loesbar: %o\n\n",
        isConsistent select "ja" else "nein";

    if not isConsistent then
        printf "Keine Loesung gefunden.\n";
        return;
    end if;

    a := solution[1];
    b := solution[2];

    printf "Loesung:\n\n";
    printf "a_%o = %o\n", number, a;
    printf "b_%o = %o\n\n", number, b;


    // Harmonisches Polynom

    h := F - a*S^6 - b*S^2*W44;

    printf "Harmonisches Polynom:\n\n";

    printf "h_%o = F_%o - (%o) S^6 - (%o) S^2 W_(4,4)\n\n",
        number, number, a, b;

    printf "h_%o =\n", number;
    printf "%o\n\n", h;


    // Kontrolle

    deltaH := ConjugateLaplace(h);

    printf "Kontrolle:\n\n";
    printf "Delta(h_%o) =\n", number;
    printf "%o\n\n", deltaH;

    printf "h_%o ist harmonisch: %o\n",
        number,
        (deltaH eq 0) select "ja" else "nein";

end procedure;
