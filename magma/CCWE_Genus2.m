// Vollstaendiger konjugierter Gewichtzaehler des Genus 2.
// Berechnungen zu Abschnitt 5.1 fuer binaere Codes vom Typ (N1,N2).
//
// x00,...,x11 gehoeren zum linken Block.
// y00,...,y11 stehen fuer die konjugierten Variablen und werden
// fuer die Rechnung als unabhaengige Variablen behandelt.
// Summiert wird ueber alle geordneten Paare aus C x C.


// Polynomring

Q := Rationals();

P<x00,x01,x10,x11,y00,y01,y10,y11> :=
    PolynomialRing(Q, 8);

XVariables := [ x00, x01, x10, x11 ];
YVariables := [ y00, y01, y10, y11 ];


// Index eines Vektors aus F_2^2

function PairIndex(a, b)
    // Reihenfolge der Paare: 00, 01, 10, 11.

    if a eq 0 and b eq 0 then
        return 1;
    elif a eq 0 and b eq 1 then
        return 2;
    elif a eq 1 and b eq 0 then
        return 3;
    else
        return 4;
    end if;
end function;


// Genus-2-CCWE

function CCWEGenus2(C, N1, N2)
    if Length(C) ne N1 + N2 then
        error "Die Codelaenge stimmt nicht mit N1+N2 ueberein.";
    end if;

    F := BaseRing(GeneratorMatrix(C));

    if Characteristic(F) ne 2 or #F ne 2 then
        error "CCWEGenus2 ist hier nur fuer binaere Codes ueber GF(2) implementiert.";
    end if;

    result := P!0;

    // Nach der Definition von ccwe(C^(2)) laufen wir ueber alle
    // geordneten Paare von Codewoertern (c,d) in C x C.
    for c in C do
        for d in C do
            monomial := P!1;

            // Erste N1 Koordinaten: x-Variablen.
            for j in [1..N1] do
                index := PairIndex(c[j], d[j]);
                monomial *:= XVariables[index];
            end for;

            // Letzte N2 Koordinaten: konjugierte Variablen.
            for j in [N1+1..N1+N2] do
                index := PairIndex(c[j], d[j]);
                monomial *:= YVariables[index];
            end for;

            result +:= monomial;
        end for;
    end for;

    return result;
end function;


// Aeussere Summe unter Beibehaltung der beiden Bloecke

// Koordinatenreihenfolge der Summe:
// links(C), links(D) | rechts(C), rechts(D).

function TypeDirectSum(C, N1, N2, D, M1, M2)
    if Length(C) ne N1 + N2 then
        error "Falsche Blocklaengen fuer C.";
    end if;

    if Length(D) ne M1 + M2 then
        error "Falsche Blocklaengen fuer D.";
    end if;

    F := BaseRing(GeneratorMatrix(C));

    if BaseRing(GeneratorMatrix(D)) ne F then
        error "C und D muessen ueber demselben Koerper definiert sein.";
    end if;

    GC := GeneratorMatrix(C);
    GD := GeneratorMatrix(D);

    rows := [];

    // Zeilen aus C.
    for i in [1..Nrows(GC)] do
        row :=
            [ GC[i,j] : j in [1..N1] ] cat
            [ F!0 : j in [1..M1] ] cat
            [ GC[i,N1+j] : j in [1..N2] ] cat
            [ F!0 : j in [1..M2] ];

        rows cat:= row;
    end for;

    // Zeilen aus D.
    for i in [1..Nrows(GD)] do
        row :=
            [ F!0 : j in [1..N1] ] cat
            [ GD[i,j] : j in [1..M1] ] cat
            [ F!0 : j in [1..N2] ] cat
            [ GD[i,M1+j] : j in [1..M2] ];

        rows cat:= row;
    end for;

    G := Matrix(
        F,
        Nrows(GC) + Nrows(GD),
        N1 + M1 + N2 + M2,
        rows
    );

    return LinearCode(G);
end function;


// 5. Kontrolle der Typ-II-Eigenschaft

function IsSignedDoublyEven(C, N1, N2)
    for c in C do
        signedWeight :=
            #[ j : j in [1..N1] | c[j] ne 0 ]
            -
            #[ j : j in [N1+1..N1+N2] | c[j] ne 0 ];

        if signedWeight mod 4 ne 0 then
            return false;
        end if;
    end for;

    return true;
end function;


procedure CheckExampleCode(C, N1, N2, name)
    printf "%o: Laenge (%o,%o), Dimension %o\n",
        name, N1, N2, Dimension(C);

    printf "  selbstdual: %o\n",
        IsSelfDual(C) select "ja" else "nein";

    printf "  doppelt-gerade bzgl. signiertem Gewicht: %o\n",
        IsSignedDoublyEven(C, N1, N2) select "ja" else "nein";
end procedure;


// Beispielcodes nach Bannai--Oura--Zhao.

F2 := GF(2);

// g_(1,1)
g11 := Matrix(F2, 1, 2, [
    1,1
]);
C11 := LinearCode(g11);


// g_(4,4)
g44 := Matrix(F2, 4, 8, [
    1,0,0,1,  0,1,1,0,
    0,1,0,1,  0,1,0,1,
    0,0,1,1,  0,0,1,1,
    0,0,0,0,  1,1,1,1
]);
C44 := LinearCode(g44);


// g^a_(6,6)
ga66 := Matrix(F2, 6, 12, [
    1,0,0,0,0,0,  1,1,1,1,1,0,
    0,1,0,0,0,0,  1,1,1,1,0,1,
    0,0,1,0,0,0,  1,1,1,0,1,1,
    0,0,0,1,0,0,  1,1,0,1,1,1,
    0,0,0,0,1,0,  1,0,1,1,1,1,
    0,0,0,0,0,1,  0,1,1,1,1,1
]);
Ca66 := LinearCode(ga66);


// g^b_(6,6)
gb66 := Matrix(F2, 6, 12, [
    1,0,0,0,1,0,  0,0,0,0,1,1,
    0,1,0,0,0,1,  0,0,0,0,1,1,
    0,0,1,0,1,1,  0,0,1,1,1,0,
    0,0,0,1,1,1,  0,0,1,1,0,1,
    0,0,0,0,0,0,  1,0,1,0,1,1,
    0,0,0,0,0,0,  0,1,0,1,1,1
]);
Cb66 := LinearCode(gb66);


// C1 = g_(1,1)^(oplus 6)
C1 := C11;
leftLength := 1;
rightLength := 1;

for i in [2..6] do
    C1 := TypeDirectSum(
        C1, leftLength, rightLength,
        C11, 1, 1
    );

    leftLength +:= 1;
    rightLength +:= 1;
end for;


// C2 = g_(4,4) oplus g_(1,1)^(oplus 2)
C11x2 := TypeDirectSum(C11, 1, 1, C11, 1, 1);
C2 := TypeDirectSum(C44, 4, 4, C11x2, 2, 2);

// C3 und C4 sind die beiden unzerlegbaren Beispielcodes.
C3 := Ca66;
C4 := Cb66;


// Die in Abschnitt 5.1 verwendeten Polynome

// Fuer C_(1,1) ist der Genus-2-CCWE genau
// S = Sum_v x_v * bar{x}_v.
S := CCWEGenus2(C11, 1, 1);

// W_(4,4) aus der Masterarbeit.
W44 := CCWEGenus2(C44, 4, 4);

// Die vier Basisinvarianten vom Grad (6,6).
F1 := CCWEGenus2(C1, 6, 6);
F2poly := CCWEGenus2(C2, 6, 6);
F3 := CCWEGenus2(C3, 6, 6);
F4 := CCWEGenus2(C4, 6, 6);
