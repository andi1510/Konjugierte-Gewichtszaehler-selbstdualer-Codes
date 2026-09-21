///////////////////////////////////////////////////////////////////////////
// Strukturanalyse fuer binaere Typ-II-Codes vom Typ (N1,N2), N1 != N2
//
// Die Analyse ist als Ergaenzung zum Kneser-Programm gedacht.
// Erwartet wird eine Liste "Klassen" der berechneten Vertreter.
//
// Geprueft werden:
//   1. allgemeine aeussere orthogonale Zerlegung in Koordinatenfaktoren
//   2. balancierte Faktoren (a,a)
//   3. einseitige Typ-II-Faktoren (0,8r) bzw. (8r,0)
//   4. unbalancierte Faktoren (a,b), a != b
//   5. linke und rechte Kerne K_L, K_R
//   6. Identitaeten pi_L(C)=K_L^perp und pi_R(C)=K_R^perp
//   7. Quotientendimensionen der allgemeinen Verklebung
//   8. volle Projektion auf den kleineren Block
//   9. verallgemeinerte Darstellung bei voller Projektion
//
// Hinweis:
//   Der triviale Code T_(N,N) und Double(K) sind als Gesamtkonstruktionen
//   nur fuer gleich grosse Bloecke definiert. Sie koennen aber als
//   balancierte Faktoren einer aeusseren Summe auftreten.
///////////////////////////////////////////////////////////////////////////


// Hilfsfunktionen

function MatrixColumns(A, cols)
    return Submatrix(A, [1..Nrows(A)], cols);
end function;

// Vergleicht die Dimensionen und Eintraege zweier Matrizen.

function MatricesEqualEntrywise(A, B)
    if Nrows(A) ne Nrows(B) or Ncols(A) ne Ncols(B) then
        return false;
    end if;

    for i in [1..Nrows(A)] do
        for j in [1..Ncols(A)] do
            if A[i,j] ne B[i,j] then
                return false;
            end if;
        end for;
    end for;

    return true;
end function;


function IsZeroMatrixEntrywise(A)
    F := BaseRing(A);

    for i in [1..Nrows(A)] do
        for j in [1..Ncols(A)] do
            if A[i,j] ne F!0 then
                return false;
            end if;
        end for;
    end for;

    return true;
end function;


// Vergleicht die Codes ueber den Rang ihrer gemeinsam
// angeordneten Generatormatrizen.

function CodesEqualByRowSpace(C, D)
    if Length(C) ne Length(D) then
        return false;
    end if;

    if Dimension(C) ne Dimension(D) then
        return false;
    end if;

    if Dimension(C) eq 0 then
        return true;
    end if;

    GC := GeneratorMatrix(C);
    GD := GeneratorMatrix(D);

    if Ncols(GC) ne Ncols(GD) then
        return false;
    end if;

    return Rank(VerticalJoin(GC, GD)) eq Dimension(C);
end function;


function ProjectionCodeOnPositions(C, pos)
    G := GeneratorMatrix(C);
    GP := MatrixColumns(G, pos);
    return LinearCode(GP);
end function;


function SupportedCodeOnPositions(C, pos)
    // Gesucht sind Codewoerter, die ausserhalb von pos verschwinden.
    // Fuer ein Codewort a*G bedeutet das a*GOutside = 0.
    // Die zugehoerigen Koeffizienten bilden den linken Nullraum.

    G := GeneratorMatrix(C);
    F := BaseRing(G);
    n := Ncols(G);
    k := Nrows(G);

    keep := [ Integers()!j : j in pos ];
    outside := [ j : j in [1..n] | j notin keep ];

    GKeep := MatrixColumns(G, keep);

    // Falls alle Koordinaten behalten werden, ist keine
    // Nullraumbedingung notwendig. Wir erzeugen den Code trotzdem
    // ueber den ausgewaehlten Spaltenblock, damit auch die Reihenfolge
    // in keep respektiert wird.
    if #outside eq 0 then
        if k eq 0 then
            return ZeroCode(F, #keep);
        end if;
        return LinearCode(GKeep);
    end if;

    GOutside := MatrixColumns(G, outside);
    KCoeffs := NullspaceMatrix(GOutside);

    // Keine nichttriviale Linearkombination verschwindet ausserhalb
    // der gewuenschten Positionen.
    if Nrows(KCoeffs) eq 0 then
        return ZeroCode(F, #keep);
    end if;

    // Die zugehoerigen Codewoerter auf die behaltenen Koordinaten
    // einschraenken.
    GSupported := KCoeffs * GKeep;

    return LinearCode(GSupported);
end function;


function SignedWeight(c, a, b)
    weightLeft := #[ j : j in [1..a] | c[j] ne 0 ];
    weightRight := #[ j : j in [a+1..a+b] | c[j] ne 0 ];
    return weightLeft - weightRight;
end function;


function IsSignedDoublyEven(C, a, b)
    for c in C do
        if SignedWeight(c, a, b) mod 4 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;


// Gewoehnliche Koordinatenkomponenten

// Bestimmt die Koordinatenkomponenten fuer die direkte Summenzerlegung.

function CoordinateComponents(C)
    G := GeneratorMatrix(C);
    E := EchelonForm(G);

    k := Dimension(C);
    n := Length(C);
    pivotColumns := [];

    for i in [1..k] do
        pivot := 0;

        for j in [1..n] do
            if E[i,j] ne 0 then
                pivot := j;
                break;
            end if;
        end for;

        if pivot eq 0 then
            error "In der Zeilenstufenform trat eine Nullzeile auf.";
        end if;

        Append(~pivotColumns, pivot);
    end for;

    pivotMatrix := MatrixColumns(G, pivotColumns);

    if Rank(pivotMatrix) ne k then
        error "Die bestimmten Pivotspalten sind nicht linear unabhaengig.";
    end if;

    systematicMatrix := pivotMatrix^(-1) * G;
    adjacency := [ [] : j in [1..n] ];

    for i in [1..k] do
        pivot := pivotColumns[i];

        for j in [1..n] do
            if j ne pivot and systematicMatrix[i,j] ne 0 then
                Append(~adjacency[pivot], j);
                Append(~adjacency[j], pivot);
            end if;
        end for;
    end for;

    visited := [ false : j in [1..n] ];
    components := [];

    for start in [1..n] do
        if not visited[start] then
            queue := [ start ];
            visited[start] := true;
            component := [];
            queuePosition := 1;

            while queuePosition le #queue do
                u := queue[queuePosition];
                queuePosition +:= 1;
                Append(~component, u);

                for v in adjacency[u] do
                    if not visited[v] then
                        visited[v] := true;
                        Append(~queue, v);
                    end if;
                end for;
            end while;

            Sort(~component);
            Append(~components, component);
        end if;
    end for;

    return components;
end function;


//    Der Faktor wird so angeordnet, dass zuerst seine linken und danach
//    seine rechten Koordinaten stehen.

function ComponentFactor(C, positions, N1)
    leftPositions := [ j : j in positions | j le N1 ];
    rightPositions := [ j : j in positions | j gt N1 ];

    Sort(~leftPositions);
    Sort(~rightPositions);

    orderedPositions := leftPositions cat rightPositions;
    factor := SupportedCodeOnPositions(C, orderedPositions);

    return factor, #leftPositions, #rightPositions;
end function;


// Kerne und allgemeine Quotientenverklebung

function GeneralKernelCodes(C, N1, N2)
    if Length(C) ne N1 + N2 then
        error "Die Blocklaengen stimmen nicht mit der Codelaenge ueberein.";
    end if;

    KLeft := SupportedCodeOnPositions(C, [1..N1]);
    KRight := SupportedCodeOnPositions(C, [N1+1..N1+N2]);

    return KLeft, KRight;
end function;


function GeneralGluingData(C, N1, N2)
    KLeft, KRight := GeneralKernelCodes(C, N1, N2);

    projectionLeft := ProjectionCodeOnPositions(C, [1..N1]);
    projectionRight := ProjectionCodeOnPositions(C, [N1+1..N1+N2]);

    quotientDimensionLeft :=
        Dimension(Dual(KLeft)) - Dimension(KLeft);

    quotientDimensionRight :=
        Dimension(Dual(KRight)) - Dimension(KRight);

    projectionIdentityLeft := CodesEqualByRowSpace(projectionLeft, Dual(KLeft));
    projectionIdentityRight := CodesEqualByRowSpace(projectionRight, Dual(KRight));

    fullLeftProjection := Dimension(projectionLeft) eq N1;
    fullRightProjection := Dimension(projectionRight) eq N2;

    return
        KLeft,
        KRight,
        projectionLeft,
        projectionRight,
        quotientDimensionLeft,
        quotientDimensionRight,
        projectionIdentityLeft,
        projectionIdentityRight,
        fullLeftProjection,
        fullRightProjection;
end function;


// Verallgemeinerte volle Projektion


//    Falls die linke Projektion voll ist, besitzt der Code nach
//    Zeilenoperationen eine Darstellung
//
//        [ I_N1 | A ]
//        [   0  | K ]
//
//    mit K = K_R und dim(K)=(N2-N1)/2.
//
//    Die Funktion liefert A und K_R.

function GeneralLeftProjectionForm(C, N1, N2)
    F := BaseRing(GeneratorMatrix(C));
    G := GeneratorMatrix(C);
    k := Nrows(G);

    // Bloecke direkt mit Magmas Submatrix extrahieren.
    GLeft  := Submatrix(G, [1..k], [1..N1]);
    GRight := Submatrix(G, [1..k], [N1+1..N1+N2]);

    // Volle linke Projektion <=> Rang(GLeft)=N1.
    if Rank(GLeft) ne N1 then
        return false,
               ZeroMatrix(F, 0, 0),
               ZeroCode(F, N2);
    end if;

    // Liftungen der Standardbasis:
    // U * GLeft = I_N1.
    U := Solution(GLeft, IdentityMatrix(F, N1));

    if Nrows(U) ne N1 or Ncols(U) ne k then
        error "Unerwartete Dimension der Loesungsmatrix fuer die linke Projektion.";
    end if;

    if not MatricesEqualEntrywise(U * GLeft, IdentityMatrix(F, N1)) then
        error "Die berechnete Loesungsmatrix U erfuellt U*GLeft = I nicht.";
    end if;

    // Kernkoeffizienten der linken Projektion:
    // N * GLeft = 0.
    N := NullspaceMatrix(GLeft);

    if Ncols(N) ne k then
        error "Unerwartete Spaltenzahl der Kernmatrix der linken Projektion.";
    end if;

    if Nrows(N) ne k - N1 then
        error "Die Kerndimension der linken Projektion stimmt nicht mit k-N1 ueberein.";
    end if;

    if not IsZeroMatrixEntrywise(N * GLeft) then
        error "Die berechneten Kernkoeffizienten verschwinden nicht auf dem linken Block.";
    end if;

    // U und N muessen zusammen eine Basis des Koeffizientenraums bilden.
    T := VerticalJoin(U, N);
    if Nrows(T) ne k or Ncols(T) ne k or Rank(T) ne k then
        error "Die aus Liftungen und Kernbasis gebildete Basistransformation ist nicht invertierbar.";
    end if;

    // Rechte Bloecke der Normalform:
    // A = U*GRight und GK = N*GRight.
    A := U * GRight;
    GK := N * GRight;

    if Nrows(GK) eq 0 then
        KRight := ZeroCode(F, N2);
    else
        KRight := LinearCode(GK);
    end if;

    // Direkte mathematische Kontrollen der Normalform:
    // [U;N] * [GLeft|GRight] = [I|A; 0|GK].
    if not MatricesEqualEntrywise(T * GLeft,
        VerticalJoin(IdentityMatrix(F, N1), ZeroMatrix(F, k-N1, N1))) then
        error "Interner Fehler: Die linke Seite der Normalform ist nicht [I;0].";
    end if;

    // Fuer einen selbstorthogonalen Code muessen diese Beziehungen gelten.
    if not MatricesEqualEntrywise(A * Transpose(A), IdentityMatrix(F, N1)) then
        error "Die Normalform wurde konstruiert, aber A*A^tr ist nicht I.";
    end if;

    if Nrows(GK) gt 0 and not IsZeroMatrixEntrywise(A * Transpose(GK)) then
        error "Die Normalform wurde konstruiert, aber A*G_KR^tr ist nicht 0.";
    end if;

    return true, A, KRight;
end function;

// Spiegelbildliche Version fuer N1 > N2. 
function GeneralRightProjectionForm(C, N1, N2)
    F := BaseRing(GeneratorMatrix(C));
    G := GeneratorMatrix(C);

    if Rank(MatrixColumns(G, [N1+1..N1+N2])) ne N2 then
        return false,
               ZeroMatrix(F, 0, 0),
               ZeroCode(F, N1);
    end if;

    // Fuer die rechte Variante ordnen wir die Bloecke temporaer um.
    GSwapped := HorizontalJoin(
        MatrixColumns(G, [N1+1..N1+N2]),
        MatrixColumns(G, [1..N1])
    );

    CSwapped := LinearCode(GSwapped);
    ok, B, KLeft := GeneralLeftProjectionForm(CSwapped, N2, N1);

    return ok, B, KLeft;
end function;


// Analyse einer einzelnen Klasse

procedure AnalyseUnequalConstructions(C, N1, N2, classNumber)
    printf "\n";
    printf "==================================================\n";
    printf "Klasse %o, Typ (%o,%o)\n", classNumber, N1, N2;
    printf "==================================================\n";

    if N1 eq N2 then
        error "Dieses Programm ist fuer N1 != N2 gedacht.";
    end if;

    if Length(C) ne N1 + N2 then
        error "Die Codelaenge stimmt nicht mit N1+N2 ueberein.";
    end if;

    F := BaseRing(GeneratorMatrix(C));

    if Characteristic(F) ne 2 or #F ne 2 then
        error "Der Algorithmus ist nur fuer binaere Codes ueber GF(2) geschrieben.";
    end if;

    printf "Dimension: %o\n", Dimension(C);
    printf "Selbstdual (binaer): %o\n",
        CodesEqualByRowSpace(C, Dual(C)) select "ja" else "nein";
    printf "Doppelt-gerade bzgl. signiertem Gewicht: %o\n",
        IsSignedDoublyEven(C, N1, N2) select "ja" else "nein";

    printf "\nKlassische Konstruktionen aus dem Fall (N,N):\n";
    printf "  Trivialer Gesamtcode T_(N,N): nicht anwendbar fuer N1 != N2\n";
    printf "  Gewoehnliches Double(K): nicht anwendbar als Gesamtkonstruktion\n";
    printf "  Hinweis: balancierte Summanden koennen selbst trivial oder Doubling-Codes sein.\n";

    // Koordinatenzerlegung

    components := CoordinateComponents(C);
    isDecomposable := #components gt 1;
    numberOneSided := 0;
    numberBalanced := 0;
    numberUnbalanced := 0;

    printf "\nAllgemeine aeussere orthogonale Zerlegung: %o\n",
        isDecomposable select "ja" else "nein";
    printf "Anzahl der Koordinatenfaktoren: %o\n", #components;

    for i in [1..#components] do
        factor, a, b := ComponentFactor(C, components[i], N1);
        factorSelfDual := CodesEqualByRowSpace(factor, Dual(factor));
        factorDoublyEven := IsSignedDoublyEven(factor, a, b);

        printf "  Faktor %o: Typ (%o,%o), Dimension %o\n",
            i, a, b, Dimension(factor);
        printf "    selbstdual: %o\n",
            factorSelfDual select "ja" else "nein";
        printf "    doppelt-gerade: %o\n",
            factorDoublyEven select "ja" else "nein";
        printf "    Differenz a-b = %o\n", a-b;

        if a eq b then
            numberBalanced +:= 1;
            if a eq 1 then
                printf "    Typ: balancierter Faktor T_(1,1)\n";
            else
                printf "    Typ: balancierter Faktor; kann mit der (N,N)-Analyse weiter untersucht werden\n";
            end if;
        elif a eq 0 or b eq 0 then
            numberOneSided +:= 1;
            printf "    Typ: einseitiger Typ-II-Faktor\n";
            printf "    gewoehnliche Koordinatenkomponenten dieses Faktors: %o\n",
                #CoordinateComponents(factor);
        else
            numberUnbalanced +:= 1;
            printf "    Typ: unbalancierter Faktor / unbalancierter Kern\n";
        end if;
    end for;

    // Quotientenverklebung

    KLeft,
    KRight,
    projectionLeft,
    projectionRight,
    qLeft,
    qRight,
    projectionIdentityLeft,
    projectionIdentityRight,
    fullLeftProjection,
    fullRightProjection := GeneralGluingData(C, N1, N2);

    printf "\nAllgemeine Quotientenverklebung:\n";
    printf "  dim K_L = %o\n", Dimension(KLeft);
    printf "  dim K_R = %o\n", Dimension(KRight);
    printf "  Erwartete Kerndifferenz (N2-N1)/2 = %o\n", (N2-N1) div 2;
    printf "  Tatsaechliche Kerndifferenz dim(K_R)-dim(K_L) = %o\n",
        Dimension(KRight)-Dimension(KLeft);
    printf "  pi_L(C) = K_L^perp: %o\n",
        projectionIdentityLeft select "ja" else "nein";
    printf "  pi_R(C) = K_R^perp: %o\n",
        projectionIdentityRight select "ja" else "nein";
    printf "  dim(K_L^perp/K_L) = %o\n", qLeft;
    printf "  dim(K_R^perp/K_R) = %o\n", qRight;
    printf "  Quotientendimensionen stimmen ueberein: %o\n",
        (qLeft eq qRight) select "ja" else "nein";

    printf "  Volle linke Projektion: %o\n",
        fullLeftProjection select "ja" else "nein";
    printf "  Volle rechte Projektion: %o\n",
        fullRightProjection select "ja" else "nein";

    // Verallgemeinerte Projektionsform

    if N1 lt N2 and fullLeftProjection then
        ok, A, K := GeneralLeftProjectionForm(C, N1, N2);

        if ok then
            printf "\nVerallgemeinerte volle Projektion auf den kleineren linken Block:\n";
            printf "  C besitzt eine Darstellung [ I_%o | A ; 0 | K_R ].\n", N1;
            printf "  Matrixgroesse A: %o x %o\n", Nrows(A), Ncols(A);
            printf "  dim K_R = %o\n", Dimension(K);
            printf "  A*A^tr = I_%o: %o\n", N1,
                MatricesEqualEntrywise(A*Transpose(A), IdentityMatrix(F,N1))
                select "ja" else "nein";

            if Dimension(K) gt 0 then
                GK := GeneratorMatrix(K);
                printf "  A*K_R^tr = 0: %o\n",
                    IsZeroMatrixEntrywise(A*Transpose(GK))
                    select "ja" else "nein";
                printf "  K_R selbstorthogonal: %o\n",
                    IsSelfOrthogonal(K) select "ja" else "nein";
            end if;
        end if;
    elif N2 lt N1 and fullRightProjection then
        ok, B, K := GeneralRightProjectionForm(C, N1, N2);

        if ok then
            printf "\nVerallgemeinerte volle Projektion auf den kleineren rechten Block:\n";
            printf "  Nach Vertauschung der Bloecke besitzt C eine Darstellung [ I_%o | B ; 0 | K_L ].\n", N2;
            printf "  Matrixgroesse B: %o x %o\n", Nrows(B), Ncols(B);
            printf "  dim K_L = %o\n", Dimension(K);
        end if;
    end if;

    // Zusammenfassung

    printf "\nZusammenfassung:\n";

    if isDecomposable then
        printf "  allgemeine aeussere orthogonale Summe\n";
    else
        printf "  koordinatenweise unzerlegbarer Code\n";
    end if;

    if numberOneSided gt 0 then
        printf "  enthaelt %o einseitige(n) Typ-II-Faktor(en)\n", numberOneSided;
    end if;

    if numberBalanced gt 0 then
        printf "  enthaelt %o balancierte(n) Faktor(en)\n", numberBalanced;
    end if;

    if numberUnbalanced gt 0 then
        printf "  enthaelt %o unbalancierte(n) Faktor(en)\n", numberUnbalanced;
    end if;

    if (N1 lt N2 and fullLeftProjection) or
       (N2 lt N1 and fullRightProjection) then
        printf "  besitzt volle Projektion auf den kleineren Block\n";
    else
        printf "  keine volle Projektion auf den kleineren Block\n";
    end if;

    printf "  allgemeine Beschreibung durch (K_L,K_R,psi) mit Quotientendimension %o\n",
        qLeft;
end procedure;


// Analyse einer gesamten Vertreterliste

procedure AnalyseAllUnequalConstructions(codeRepresentatives, N1, N2)
    if N1 eq N2 then
        error "Dieses Programm ist fuer N1 != N2 gedacht.";
    end if;

    printf "==================================================\n";
    printf "Strukturanalyse fuer Typ (%o,%o)\n", N1, N2;
    printf "Anzahl der Vertreter: %o\n", #codeRepresentatives;
    printf "==================================================\n";

    numberDecomposable := 0;
    numberIndecomposable := 0;
    numberWithOneSidedFactor := 0;
    numberWithFullSmallProjection := 0;
    quotientDimensions := [];

    for i in [1..#codeRepresentatives] do
        C := codeRepresentatives[i];
        components := CoordinateComponents(C);

        if #components gt 1 then
            numberDecomposable +:= 1;
        else
            numberIndecomposable +:= 1;
        end if;

        hasOneSided := false;
        for positions in components do
            factor, a, b := ComponentFactor(C, positions, N1);
            if a eq 0 or b eq 0 then
                hasOneSided := true;
            end if;
        end for;

        if hasOneSided then
            numberWithOneSidedFactor +:= 1;
        end if;

        KLeft,
        KRight,
        projectionLeft,
        projectionRight,
        qLeft,
        qRight,
        projectionIdentityLeft,
        projectionIdentityRight,
        fullLeftProjection,
        fullRightProjection := GeneralGluingData(C, N1, N2);

        Append(~quotientDimensions, qLeft);

        if (N1 lt N2 and fullLeftProjection) or
           (N2 lt N1 and fullRightProjection) then
            numberWithFullSmallProjection +:= 1;
        end if;

        AnalyseUnequalConstructions(C, N1, N2, i);
    end for;

    printf "\n";
    printf "==================================================\n";
    printf "GESAMTAUSWERTUNG fuer (%o,%o)\n", N1, N2;
    printf "==================================================\n";
    printf "Klassen insgesamt: %o\n", #codeRepresentatives;
    printf "Allgemein aeusserlich zerlegbar: %o\n", numberDecomposable;
    printf "Koordinatenweise unzerlegbar: %o\n", numberIndecomposable;
    printf "Mit mindestens einem einseitigen Typ-II-Faktor: %o\n",
        numberWithOneSidedFactor;
    printf "Mit voller Projektion auf den kleineren Block: %o\n",
        numberWithFullSmallProjection;

    qValues := Setseq(Seqset(quotientDimensions));
    Sort(~qValues);

    printf "Quotientendimensionen:\n";
    for q in qValues do
        frequency := #[ x : x in quotientDimensions | x eq q ];
        printf "  q = %o: %o Klassen\n", q, frequency;
    end for;
end procedure;


procedure WriteUnequalConstructionAnalysis(
    codeRepresentatives,
    N1,
    N2,
    fileName
)
    SetOutputFile(fileName : Overwrite := true);

    try
        AnalyseAllUnequalConstructions(codeRepresentatives, N1, N2);
        UnsetOutputFile();
    catch e
        UnsetOutputFile();
        print "FEHLER waehrend der Strukturanalyse:";
        print e`Object;
        error e;
    end try;
end procedure;


// AUFRUF
//
// Am Ende des Kneser-Programms, nachdem "Klassen", N1 und N2
// berechnet wurden:
//
// WriteUnequalConstructionAnalysis(
//     Klassen,
//     N1,
//     N2,
//     Sprintf("Konstruktionen_ungleich_%o_%o.txt", N1, N2)
// );
