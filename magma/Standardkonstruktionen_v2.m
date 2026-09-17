///////////////////////////////////////////////////////////////////////////
// Standardkonstruktionen fuer binaere Codes vom Typ (N,N)
//
// Geprueft werden:
//   1. trivialer Code T_(N,N)
//   2. aeussere orthogonale Summen
//   3. Doubling-Konstruktion Double(K)
//   4. volle linke bzw. rechte Projektion
//
// Die Aequivalenz wird blockweise unter S_N x S_N geprueft.
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// 1. Grundlegende Hilfsfunktionen
///////////////////////////////////////////////////////////////////////////

function MatrixColumns(A, cols)
    F := BaseRing(A);

    if #cols eq 0 then
        return ZeroMatrix(F, Nrows(A), 0);
    end if;

    return Matrix(
        F,
        Nrows(A),
        #cols,
        [ A[i,j] : i in [1..Nrows(A)], j in cols ]
    );
end function;


/*
    Projektion eines Codes auf die Koordinaten in pos.
*/
function ProjectionCodeOnPositions(C, pos)
    G := GeneratorMatrix(C);
    GP := MatrixColumns(G, pos);

    return LinearCode(GP);
end function;


/*
    Bestimmt alle Codewoerter, deren Traeger vollstaendig in pos liegt,
    und entfernt danach die uebrigen Koordinaten.

    Fuer die kleinen Dimensionen in den betrachteten Beispielen wird
    dies direkt ueber die Codewoerter bestimmt. Das vermeidet
    Orientierungsfragen bei rechteckigen Nullraum-Berechnungen.
*/
function SupportedCodeOnPositions(C, pos)
    G := GeneratorMatrix(C);
    F := BaseRing(G);
    n := Length(C);

    outside := [ j : j in [1..n] | j notin pos ];

    // Fuer die hier betrachteten Codes (Dimension hoechstens ca. 11)
    // ist die direkte Codewortpruefung sehr klein und besonders robust.
    // Wir sammeln genau die Codewoerter, die ausserhalb von pos nur
    // Nullen besitzen, und beschraenken sie anschliessend auf pos.
    supportedRows := [];

    for c in C do
        isSupported := true;

        for j in outside do
            if c[j] ne F!0 then
                isSupported := false;
                break;
            end if;
        end for;

        if isSupported then
            row := [ c[j] : j in pos ];

            // Die Nullzeile muss nicht als Erzeuger gespeichert werden.
            if exists{ a : a in row | a ne F!0 } then
                Append(~supportedRows, row);
            end if;
        end if;
    end for;

    if #supportedRows eq 0 then
        return ZeroCode(F, #pos);
    end if;

    M := Matrix(F, #supportedRows, #pos, &cat supportedRows);
    return LinearCode(M);
end function;


function BinaryBiweightDistribution(C, N1, N2)
    distribution := [ 0 : i in [1..(N1 + 1)*(N2 + 1)] ];

    for c in C do
        weightLeft := #[ j : j in [1..N1] | c[j] ne 0 ];
        weightRight := #[ j : j in [N1+1..N1+N2] | c[j] ne 0 ];

        index := weightLeft*(N2 + 1) + weightRight + 1;
        distribution[index] +:= 1;
    end for;

    return distribution;
end function;


///////////////////////////////////////////////////////////////////////////
// 2. Blockerhaltende Aequivalenz unter S_N1 x S_N2
///////////////////////////////////////////////////////////////////////////

/*
    Markierter Inzidenzgraph des Codes.

    Label 1: linke Koordinaten
    Label 2: rechte Koordinaten
    Label 3: Codewoerter

    Ein Graphisomorphismus muss dadurch linke Koordinaten auf linke,
    rechte auf rechte und Codewoerter auf Codewoerter abbilden.
*/
function BlockIncidenceGraph(C, N1, N2)
    n := N1 + N2;

    if Length(C) ne n then
        error "Die Codelaenge stimmt nicht mit N1+N2 ueberein.";
    end if;

    F := BaseRing(GeneratorMatrix(C));
    words := [ c : c in C ];
    numberOfVertices := n + #words;

    // Genau wie im funktionierenden Kneser-Programm wird der Graph
    // ueber eine Nachbarfolge der Laenge numberOfVertices konstruiert.
    neighbours := [ { Integers() | } : i in [1..numberOfVertices] ];

    for t in [1..#words] do
        c := words[t];
        wordVertex := n + t;

        for j in [1..n] do
            if c[j] eq F!1 then
                Include(~neighbours[j], wordVertex);
                Include(~neighbours[wordVertex], j);
            end if;
        end for;
    end for;

    graph := Graph< numberOfVertices | neighbours : SparseRep := true >;

    labels :=
        [ "L" : j in [1..N1] ] cat
        [ "R" : j in [1..N2] ] cat
        [ "W" : j in [1..#words] ];

    AssignLabels(VertexSet(graph), labels);

    return graph;
end function;


function AreBlockEquivalent(C, D, N1, N2)
    if Length(C) ne Length(D) then
        return false;
    end if;

    if Dimension(C) ne Dimension(D) then
        return false;
    end if;

    if #C ne #D then
        return false;
    end if;

    if WeightDistribution(C) ne WeightDistribution(D) then
        return false;
    end if;

    if BinaryBiweightDistribution(C, N1, N2) ne
       BinaryBiweightDistribution(D, N1, N2) then
        return false;
    end if;

    graphC := BlockIncidenceGraph(C, N1, N2);
    graphD := BlockIncidenceGraph(D, N1, N2);

    isomorphic := IsIsomorphic(graphC, graphD);

    return isomorphic;
end function;


///////////////////////////////////////////////////////////////////////////
// 3. Gewoehnliche Koordinatenkomponenten des Codes
///////////////////////////////////////////////////////////////////////////

/*
    Berechnet die Zusammenhangskomponenten des Spaltenmatroids einer
    Generatormatrix. Diese Komponenten beschreiben die gewoehnliche
    direkte Summenzerlegung des Codes nach Koordinaten.
*/
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


function UnionOfComponents(components, indices)
    positions := [];

    for index in indices do
        positions cat:= components[index];
    end for;

    Sort(~positions);

    return positions;
end function;


function IsBalancedPositionSet(positions, N1)
    numberLeft := #[ j : j in positions | j le N1 ];
    numberRight := #[ j : j in positions | j gt N1 ];

    return numberLeft eq numberRight;
end function;


/*
    Fasst die gewoehnlichen Koordinatenkomponenten zu minimalen
    balancierten Gruppen zusammen. Jede Gruppe beschreibt einen Faktor
    vom Typ (m,m).
*/
function TypeFactorSupports(C, N1, N2)
    if N1 + N2 ne Length(C) then
        error "Die Blocklaengen stimmen nicht mit der Codelaenge ueberein.";
    end if;

    if N1 ne N2 then
        error "Diese Zerlegung ist fuer Codes vom Typ (N,N) formuliert.";
    end if;

    components := CoordinateComponents(C);
    remaining := [ 1..#components ];
    factorSupports := [];

    while #remaining gt 0 do
        first := remaining[1];
        others := [ remaining[i] : i in [2..#remaining] ];

        best := remaining;
        bestSize := #remaining;

        for mask in [0..2^(#others)-1] do
            candidate := [ first ];

            for t in [1..#others] do
                bit := (mask div 2^(t-1)) mod 2;

                if bit eq 1 then
                    Append(~candidate, others[t]);
                end if;
            end for;

            positions := UnionOfComponents(components, candidate);

            if IsBalancedPositionSet(positions, N1) and
               #candidate lt bestSize then
                best := candidate;
                bestSize := #candidate;
            end if;
        end for;

        factorPositions := UnionOfComponents(components, best);

        if not IsBalancedPositionSet(factorPositions, N1) then
            error "Es konnte keine balancierte Faktorzerlegung bestimmt werden.";
        end if;

        Append(~factorSupports, factorPositions);
        remaining := [ i : i in remaining | i notin best ];
    end while;

    return factorSupports;
end function;


function TypeFactorCode(C, positions, N1)
    leftPositions := [ j : j in positions | j le N1 ];
    rightPositions := [ j : j in positions | j gt N1 ];

    Sort(~leftPositions);
    Sort(~rightPositions);

    orderedPositions := leftPositions cat rightPositions;
    factor := SupportedCodeOnPositions(C, orderedPositions);

    return factor, #leftPositions;
end function;


///////////////////////////////////////////////////////////////////////////
// 4. Trivialer Code T_(N,N)
///////////////////////////////////////////////////////////////////////////

function TrivialTypeCode(F, N)
    identity := IdentityMatrix(F, N);
    generatorMatrix := HorizontalJoin(identity, identity);

    return LinearCode(generatorMatrix);
end function;


function IsTrivialTypeCode(C, N)
    F := BaseRing(GeneratorMatrix(C));
    factors := TypeFactorSupports(C, N, N);

    if #factors ne N then
        return false;
    end if;

    trivialOne := LinearCode(Matrix(F, 1, 2, [ 1, 1 ]));

    for positions in factors do
        factor, m := TypeFactorCode(C, positions, N);

        if m ne 1 then
            return false;
        end if;

        if factor ne trivialOne then
            return false;
        end if;
    end for;

    return true;
end function;


///////////////////////////////////////////////////////////////////////////
// 5. Doubling-Konstruktion
///////////////////////////////////////////////////////////////////////////

function DoubleCode(K)
    return PlotkinSum(Dual(K), K);
end function;


/*
    K ist der linke Kern
        K = {x | (x,0) in C}.

    Ist C unter S_N x S_N aequivalent zu Double(K0), dann ist der
    linke Kern eine Koordinatenpermutation von K0. Daher wird C mit
    Double(K) blockweise verglichen.
*/
function IsDoublingCode(C, N)
    if Length(C) ne 2*N then
        error "Der Code muss die Laenge 2*N besitzen.";
    end if;

    leftPositions := [ 1..N ];
    K := SupportedCodeOnPositions(C, leftPositions);
    candidate := DoubleCode(K);

    isDoubling := AreBlockEquivalent(C, candidate, N, N);

    return isDoubling, K;
end function;


///////////////////////////////////////////////////////////////////////////
// 6. Volle Projektionen und Darstellung (I_N | A)
///////////////////////////////////////////////////////////////////////////

function LeftProjectionMatrix(C, N)
    F := BaseRing(GeneratorMatrix(C));
    G := GeneratorMatrix(C);

    if Length(C) ne 2*N or Dimension(C) ne N then
        return false, ZeroMatrix(F, 0, 0);
    end if;

    GLeft := MatrixColumns(G, [1..N]);
    GRight := MatrixColumns(G, [N+1..2*N]);

    if Rank(GLeft) ne N then
        return false, ZeroMatrix(F, 0, 0);
    end if;

    A := GLeft^(-1) * GRight;

    return true, A;
end function;


function RightProjectionMatrix(C, N)
    F := BaseRing(GeneratorMatrix(C));
    G := GeneratorMatrix(C);

    if Length(C) ne 2*N or Dimension(C) ne N then
        return false, ZeroMatrix(F, 0, 0);
    end if;

    GLeft := MatrixColumns(G, [1..N]);
    GRight := MatrixColumns(G, [N+1..2*N]);

    if Rank(GRight) ne N then
        return false, ZeroMatrix(F, 0, 0);
    end if;

    B := GRight^(-1) * GLeft;

    return true, B;
end function;


function IsPermutationMatrixBinary(A)
    if Nrows(A) ne Ncols(A) then
        return false;
    end if;

    n := Nrows(A);

    for i in [1..n] do
        rowWeight := #[ j : j in [1..n] | A[i,j] ne 0 ];

        if rowWeight ne 1 then
            return false;
        end if;
    end for;

    for j in [1..n] do
        columnWeight := #[ i : i in [1..n] | A[i,j] ne 0 ];

        if columnWeight ne 1 then
            return false;
        end if;
    end for;

    return true;
end function;


///////////////////////////////////////////////////////////////////////////
// 7. Eigenschaften des Ausgangscodes K
///////////////////////////////////////////////////////////////////////////

function IsDoublyEvenBinaryCode(K)
    for v in K do
        if Weight(v) mod 4 ne 0 then
            return false;
        end if;
    end for;

    return true;
end function;


function SupportSizeOfCode(K)
    G := GeneratorMatrix(K);
    supportSize := 0;

    for j in [1..Length(K)] do
        nonzeroColumn := exists{ i : i in [1..Nrows(G)] | G[i,j] ne 0 };

        if nonzeroColumn then
            supportSize +:= 1;
        end if;
    end for;

    return supportSize;
end function;



///////////////////////////////////////////////////////////////////////////
// 8. Verklebungsanalyse: linke und rechte Kerne
///////////////////////////////////////////////////////////////////////////

/*
    Fuer C <= F_2^N x F_2^N definieren wir

        K_L = { x | (x,0) in C },
        K_R = { y | (0,y) in C }.

    Die Funktion SupportedCodeOnPositions liefert dabei bereits die
    auf N Koordinaten eingeschraenkten Codes.
*/
function LeftRightKernelCodes(C, N)
    if Length(C) ne 2*N then
        error "Der Code muss die Laenge 2*N besitzen.";
    end if;

    KLeft := SupportedCodeOnPositions(C, [1..N]);
    KRight := SupportedCodeOnPositions(C, [N+1..2*N]);

    return KLeft, KRight;
end function;


/*
    Gewoehnliche Permutationsaequivalenz zweier Codes gleicher Laenge.

    Der Inzidenzgraph besitzt zwei Knotentypen:
      "C" = Koordinaten,
      "W" = Codewoerter.

    Ein labelerhaltender Graphisomorphismus entspricht genau einer
    Koordinatenpermutation, welche die beiden Codes ineinander ueberfuehrt.
*/
function OneBlockCodeGraph(C)
    n := Length(C);
    F := BaseRing(GeneratorMatrix(C));
    words := [ c : c in C ];
    numberOfVertices := n + #words;

    neighbours := [ { Integers() | } : i in [1..numberOfVertices] ];

    for t in [1..#words] do
        c := words[t];
        wordVertex := n + t;

        for j in [1..n] do
            if c[j] eq F!1 then
                Include(~neighbours[j], wordVertex);
                Include(~neighbours[wordVertex], j);
            end if;
        end for;
    end for;

    graph := Graph< numberOfVertices | neighbours : SparseRep := true >;

    labels :=
        [ "C" : j in [1..n] ] cat
        [ "W" : j in [1..#words] ];

    AssignLabels(VertexSet(graph), labels);

    return graph;
end function;


function ArePermutationEquivalentCodes(C, D)
    if Length(C) ne Length(D) then
        return false;
    end if;

    if Dimension(C) ne Dimension(D) then
        return false;
    end if;

    if WeightDistribution(C) ne WeightDistribution(D) then
        return false;
    end if;

    graphC := OneBlockCodeGraph(C);
    graphD := OneBlockCodeGraph(D);

    return IsIsomorphic(graphC, graphD);
end function;


/*
    Liefert die grundlegenden Daten der Quotientenverklebung.

    Fuer einen selbstdualen Code gilt theoretisch
        pi_L(C) = K_L^perp,
        pi_R(C) = K_R^perp.

    Die Dimension der Quotienten ist dann
        dim(K_L^perp/K_L) = N - 2 dim(K_L)
    und analog rechts.
*/
function GluingData(C, N)
    KLeft, KRight := LeftRightKernelCodes(C, N);

    projectionLeft := ProjectionCodeOnPositions(C, [1..N]);
    projectionRight := ProjectionCodeOnPositions(C, [N+1..2*N]);

    kernelsEquivalent :=
        ArePermutationEquivalentCodes(KLeft, KRight);

    quotientDimensionLeft :=
        Dimension(Dual(KLeft)) - Dimension(KLeft);

    quotientDimensionRight :=
        Dimension(Dual(KRight)) - Dimension(KRight);

    projectionIdentityLeft :=
        projectionLeft eq Dual(KLeft);

    projectionIdentityRight :=
        projectionRight eq Dual(KRight);

    return
        KLeft,
        KRight,
        kernelsEquivalent,
        quotientDimensionLeft,
        quotientDimensionRight,
        projectionIdentityLeft,
        projectionIdentityRight;
end function;


///////////////////////////////////////////////////////////////////////////
// 9. Analyse einer einzelnen Klasse
///////////////////////////////////////////////////////////////////////////

procedure AnalyseStandardConstructions(C, N1, N2, classNumber)
    printf "\n";
    printf "==================================================\n";
    printf "Klasse %o\n", classNumber;
    printf "==================================================\n";

    if Length(C) ne N1 + N2 then
        error "Die Codelaenge stimmt nicht mit N1+N2 ueberein.";
    end if;

    F := BaseRing(GeneratorMatrix(C));

    if Characteristic(F) ne 2 or #F ne 2 then
        error "Der Algorithmus ist nur fuer binaere Codes ueber GF(2) geschrieben.";
    end if;

    if N1 ne N2 then
        printf "Die Standardkonstruktionen werden nur fuer Codes vom Typ (N,N) ausgewertet.\n";
        return;
    end if;

    N := N1;

    ///////////////////////////////////////////////////////////////
    // Trivialer Code
    ///////////////////////////////////////////////////////////////

    isTrivial := IsTrivialTypeCode(C, N);

    printf "Trivialer Code T_(%o,%o): %o\n",
        N, N, isTrivial select "ja" else "nein";

    ///////////////////////////////////////////////////////////////
    // Aeussere orthogonale Zerlegung
    ///////////////////////////////////////////////////////////////

    factorSupports := TypeFactorSupports(C, N, N);
    isDecomposable := #factorSupports gt 1;

    printf "Aeusserlich orthogonal zerlegbar: %o\n",
        isDecomposable select "ja" else "nein";
    printf "Anzahl der Faktoren vom Typ (m,m): %o\n", #factorSupports;

    for factorNumber in [1..#factorSupports] do
        positions := factorSupports[factorNumber];
        leftPositions := [ j : j in positions | j le N ];
        rightPositionsLocal := [ j-N : j in positions | j gt N ];
        factor, m := TypeFactorCode(C, positions, N);

        printf "  Faktor %o: Laenge (%o,%o), Dimension %o\n",
            factorNumber, m, m, Dimension(factor);
        printf "    linke Koordinaten: %o\n", leftPositions;
        printf "    rechte Koordinaten: %o\n", rightPositionsLocal;

        if m eq 1 then
            printf "    Standardbezeichnung: T_(1,1)\n";
        elif m eq 4 then
            printf "    Standardbezeichnung: g_(4,4)\n";
        else
            printf "    Standardbezeichnung: kein Name aus dem Paper\n";
        end if;

        printf "    Generatormatrix des Faktors:\n";
        printf "%o\n", GeneratorMatrix(factor);
    end for;

    ///////////////////////////////////////////////////////////////
    // Doubling
    ///////////////////////////////////////////////////////////////

    isDoubling, K := IsDoublingCode(C, N);

    printf "Doubling-Konstruktion: %o\n",
        isDoubling select "ja" else "nein";

    if isDoubling then
        printf "  C ist unter S_%o x S_%o aequivalent zu Double(K).\n", N, N;
        printf "  Dimension von K: %o\n", Dimension(K);
        printf "  Traegergroesse von K: %o\n", SupportSizeOfCode(K);
        printf "  K ist selbstorthogonal: %o\n",
            IsSelfOrthogonal(K) select "ja" else "nein";
        printf "  K ist doppelt-gerade: %o\n",
            IsDoublyEvenBinaryCode(K) select "ja" else "nein";
        printf "  Gewichtsverteilung von K:\n";
        printf "  %o\n", WeightDistribution(K);
        printf "  Generatormatrix von K:\n";
        printf "%o\n", GeneratorMatrix(K);
    end if;

    ///////////////////////////////////////////////////////////////
    // Projektionen
    ///////////////////////////////////////////////////////////////

    hasLeftProjection, A := LeftProjectionMatrix(C, N);

    printf "Volle linke Projektion: %o\n",
        hasLeftProjection select "ja" else "nein";

    if hasLeftProjection then
        identity := IdentityMatrix(F, N);
        one := Vector(F, [ 1 : i in [1..N] ]);

        printf "  Darstellung: (I_%o | A)\n", N;
        printf "  A*A^tr = I: %o\n",
            (A * Transpose(A) eq identity) select "ja" else "nein";
        printf "  1*A = 1: %o\n",
            (one * A eq one) select "ja" else "nein";
        printf "  A ist eine Permutationsmatrix: %o\n",
            IsPermutationMatrixBinary(A) select "ja" else "nein";
        printf "  Matrix A:\n";
        printf "%o\n", A;
    end if;

    hasRightProjection, B := RightProjectionMatrix(C, N);

    printf "Volle rechte Projektion: %o\n",
        hasRightProjection select "ja" else "nein";

    if hasRightProjection then
        printf "  Darstellung: (B | I_%o)\n", N;
        printf "  Matrix B:\n";
        printf "%o\n", B;
    end if;

    ///////////////////////////////////////////////////////////////
    // Verklebungsanalyse
    ///////////////////////////////////////////////////////////////

    KLeft, KRight, kernelsEquivalent,
    quotientDimensionLeft, quotientDimensionRight,
    projectionIdentityLeft, projectionIdentityRight :=
        GluingData(C, N);

    printf "Verklebungsdaten:\n";
    printf "  Dimension K_L: %o\n", Dimension(KLeft);
    printf "  Dimension K_R: %o\n", Dimension(KRight);
    printf "  Gewichtsverteilung K_L: %o\n", WeightDistribution(KLeft);
    printf "  Gewichtsverteilung K_R: %o\n", WeightDistribution(KRight);
    printf "  K_L und K_R permutationsaequivalent: %o\n",
        kernelsEquivalent select "ja" else "nein";
    printf "  Quotientendimension links: %o\n", quotientDimensionLeft;
    printf "  Quotientendimension rechts: %o\n", quotientDimensionRight;
    printf "  pi_L(C) = K_L^perp: %o\n",
        projectionIdentityLeft select "ja" else "nein";
    printf "  pi_R(C) = K_R^perp: %o\n",
        projectionIdentityRight select "ja" else "nein";

    if Dimension(KLeft) gt 0 then
        printf "  Generatormatrix K_L:\n";
        printf "%o\n", GeneratorMatrix(KLeft);
    end if;

    if Dimension(KRight) gt 0 then
        printf "  Generatormatrix K_R:\n";
        printf "%o\n", GeneratorMatrix(KRight);
    end if;

    ///////////////////////////////////////////////////////////////
    // Zusammenfassung
    ///////////////////////////////////////////////////////////////

    printf "Zusammenfassung:\n";

    if isTrivial then
        printf "  T_(%o,%o) = Double(0)\n", N, N;
    elif isDecomposable and isDoubling then
        printf "  aeussere orthogonale Summe und Doubling-Code\n";
    elif isDecomposable then
        printf "  aeussere orthogonale Summe, kein Doubling erkannt\n";
    elif isDoubling then
        printf "  unzerlegbarer Doubling-Code\n";
    elif hasLeftProjection or hasRightProjection then
        printf "  unzerlegbarer Code mit voller Projektion\n";
    elif kernelsEquivalent then
        printf "  Twisted-Doubling-Kandidat:\n";
        printf "  K_L und K_R sind permutationsaequivalent, aber C ist kein gewoehnliches Double(K).\n";
        printf "  Der Code wird durch eine nichttriviale Verklebung der Quotienten beschrieben.\n";
    else
        printf "  asymmetrische Verklebung:\n";
        printf "  K_L und K_R sind nicht permutationsaequivalent.\n";
        printf "  Der Code wird durch eine Verklebung verschiedener Quotientencodes beschrieben.\n";
    end if;
end procedure;


///////////////////////////////////////////////////////////////////////////
// 10. Analyse einer gesamten Vertreterliste
///////////////////////////////////////////////////////////////////////////

procedure AnalyseAllStandardConstructions(codeRepresentatives, N1, N2)
    printf "==================================================\n";
    printf "Analyse der Standardkonstruktionen fuer (%o,%o)\n", N1, N2;
    printf "Anzahl der Vertreter: %o\n", #codeRepresentatives;
    printf "==================================================\n";

    for i in [1..#codeRepresentatives] do
        AnalyseStandardConstructions(codeRepresentatives[i], N1, N2, i);
    end for;
end procedure;


procedure WriteStandardConstructionAnalysis(
    codeRepresentatives,
    N1,
    N2,
    fileName
)
    SetOutputFile(fileName : Overwrite := true);

    try
        AnalyseAllStandardConstructions(codeRepresentatives, N1, N2);
        UnsetOutputFile();
    catch e
        // Wichtig: Datei auch bei einem Fehler schliessen, damit bereits
        // erzeugte Ausgabe geschrieben wird und die Konsole wieder sichtbar ist.
        UnsetOutputFile();
        print "FEHLER waehrend der Standardkonstruktionsanalyse:";
        print e`Object;
        error e;
    end try;
end procedure;


///////////////////////////////////////////////////////////////////////////
// 11. AUFRUF
///////////////////////////////////////////////////////////////////////////

// Ausgabe auf dem Bildschirm:
// AnalyseAllStandardConstructions(Klassen, N1, N2);

// Ausgabe in eine Datei:
// WriteStandardConstructionAnalysis(
//     Klassen,
//     N1,
//     N2,
//     Sprintf("Standardkonstruktionen_%o_%o.txt", N1, N2)
// );
