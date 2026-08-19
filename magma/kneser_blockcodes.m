///////////////////////////////////////////////////////////////////////////
// Sicherer Kneser-Algorithmus fuer binaere Typ-II-Codes der Laenge (N1,N2)
// Exakte Hyperebenenorbits unter der wirklichen Codewirkung des Stabilisators
//
// Gewicht gemaess Bannai--Oura--Zhao:
// wt(c) = Summe(linker Block) - Summe(rechter Block).
///////////////////////////////////////////////////////////////////////////

// Aufruf:
// Genus, Stabilisatoren := KneserBlockCodes(4);
// Genus, Stabilisatoren := KneserBlockCodesGeneral(1, 9);

KneserBlockCodesGeneral := function(N1, N2)

Laenge := N1 + N2;
DimensionCode := Laenge div 2;
p := 2;

if N1 lt 0 or N2 lt 0 or Laenge eq 0 then
    error "N1 und N2 muessen nichtnegativ sein und duerfen nicht beide null sein.";
end if;

if (Laenge mod 2) ne 0 then
    error "Fuer einen selbstdualen Code muss N1 + N2 gerade sein.";
end if;

if ((N1-N2) mod 8) ne 0 then
    error "Ein binaerer Typ-II-Code erfordert N1 - N2 = 0 modulo 8.";
end if;

F := GF(p);

//////////////////////////////
// Dualitaet und signiertes Gewicht
//////////////////////////////

SignedWeight := function(v)
    links := 0;
    rechts := 0;

    for i in [1..N1] do
        links +:= Integers()!v[i];
    end for;

    for i in [1..N2] do
        rechts +:= Integers()!v[N1+i];
    end for;

    return links - rechts;
end function;

IsSelfDual := function(C)
    // Bei Laenge 2n ist ein selbstorthogonaler n-dimensionaler Code
    // automatisch selbstdual.
    if Dimension(C) ne DimensionCode then
        return false;
    end if;
    return C eq Dual(C);
end function;

IsDoublyEven := function(C)
    // Fuer einen selbstorthogonalen binaeren Code genuegt der Test an
    // einer Basis: Die Kreuzterme sind modulo 4 gleich null.
    g := GeneratorMatrix(C);

    for i in [1..Nrows(g)] do
        if (SignedWeight(g[i]) mod 4) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

//////////////////////////////
// Gefaerbter Inzidenzgraph
//////////////////////////////

// Der Graph besitzt drei durch Labels getrennte Knotentypen:
//   L = linke Koordinaten, R = rechte Koordinaten, W = Codewoerter.
// Ein labelerhaltender Graphisomorphismus ist daher genau eine
// Codeisomorphie unter S_n x S_n.
CodeGraph := function(C)
    Woerter := [w : w in C];
    AnzahlWoerter := #Woerter;
    AnzahlKnoten := Laenge + AnzahlWoerter;
    Nachbarn := [{ Integers() | } : i in [1..AnzahlKnoten]];

    for a in [1..AnzahlWoerter] do
        WortKnoten := Laenge + a;
        w := Woerter[a];

        for i in [1..Laenge] do
            if w[i] eq F!1 then
                Include(~Nachbarn[i], WortKnoten);
                Include(~Nachbarn[WortKnoten], i);
            end if;
        end for;
    end for;

    X := Graph< AnzahlKnoten | Nachbarn : SparseRep := true >;
    Labels := ["L" : i in [1..N1]] cat
              ["R" : i in [1..N2]] cat
              ["W" : i in [1..AnzahlWoerter]];
    AssignLabels(VertexSet(X), Labels);

    return X;
end function;

IsIsomorphicUnderBlocks := function(X1, X2)
    return IsIsomorphic(X1, X2);
end function;

BlockStabilizer := function(C)
    // Vollstaendige Koordinaten-Automorphismengruppe per Backtracking.
    A := PermutationGroup(C);

    // Nur Automorphismen, die den linken Block mengenweise festhalten.
    LinkerBlock := {1..N1};
    return Stabilizer(A, LinkerBlock);
end function;

//////////////////////////////
// Automatischer Startcode
//////////////////////////////

// Zuerst min(N1,N2) Kopien von g_(1,1). Die uebrigen Koordinaten
// liegen vollstaendig in einem Block und werden mit E8-Bloecken gefuellt.
Startmatrix := ZeroMatrix(F, DimensionCode, Laenge);
Gemeinsam := Minimum(N1, N2);
Zeile := 1;

for i in [1..Gemeinsam] do
    Startmatrix[Zeile][i] := 1;
    Startmatrix[Zeile][N1+i] := 1;
    Zeile +:= 1;
end for;

E8 := Matrix(F, 4, 8, [
    1,0,0,0,0,1,1,1,
    0,1,0,0,1,0,1,1,
    0,0,1,0,1,1,0,1,
    0,0,0,1,1,1,1,0
]);

Differenz := Abs(N1-N2);

if Differenz gt 0 then
    if N1 gt N2 then
        RestStart := Gemeinsam + 1;
    else
        RestStart := N1 + Gemeinsam + 1;
    end if;

    for b in [0..(Differenz div 8)-1] do
        for i in [1..4] do
            for j in [1..8] do
                Startmatrix[Zeile+i-1][RestStart+8*b+j-1] := E8[i][j];
            end for;
        end for;
        Zeile +:= 4;
    end for;
end if;

Startcode := LinearCode(Startmatrix);

if not IsSelfDual(Startcode) then
    error "Der Startcode ist nicht selbstdual.";
end if;

if not IsDoublyEven(Startcode) then
    error "Der Startcode ist bezueglich des signierten Gewichts nicht doppelt-gerade.";
end if;

Genus := [Startcode];
CodeGraphen := [CodeGraph(Startcode)];
Stabilisatoren := [];

//////////////////////////////
// Kneser-Algorithmus
//////////////////////////////

k := 1;

while k le #Genus do
    Caktuell := Genus[k];
    genmat := GeneratorMatrix(Caktuell);

    // Stabilisator fuer Orbitberechnung und spaetere Ausgabe
    G1 := BlockStabilizer(Caktuell);
    Append(~Stabilisatoren, G1);

    //////////////////////////
    // Alle Hyperebenen von Caktuell
    //////////////////////////

    // Ueber GF(2) bestimmt jeder von null verschiedene Vektor a genau
    // eine Hyperebene ker(a). Zunaechst werden garantiert alle maximalen
    // Teilcodes von Caktuell konstruiert.
    Vinfo := VectorSpace(F, DimensionCode);
    Hyperebenen := [];
    HyperebenenCodes := [];

    for a in Vinfo do
        if a eq Vinfo!0 then
            continue;
        end if;

        A := Matrix(F, DimensionCode, 1, Eltseq(a));
        B := BasisMatrix(Kernel(A));
        Append(~Hyperebenen, B);
        Append(~HyperebenenCodes, LinearCode(B*genmat));
    end for;

    // Sichere Orbitreduktion: G1 wirkt hier direkt durch
    // Koordinatenpermutationen auf den Hyperebenen-Codes. Damit wird keine
    // kuenstliche Matrixdarstellung verwendet. Aus jedem tatsaechlichen
    // G1-Orbit wird genau ein Vertreter ausgewaehlt.
    Verwendet := [false : i in [1..#Hyperebenen]];
    HyperebenenVertreter := [];

    for i in [1..#Hyperebenen] do
        if Verwendet[i] then
            continue;
        end if;

        Append(~HyperebenenVertreter, Hyperebenen[i]);
        Orb := HyperebenenCodes[i]^G1;

        for j in [i..#Hyperebenen] do
            if not Verwendet[j] and HyperebenenCodes[j] in Orb then
                Verwendet[j] := true;
            end if;
        end for;
    end for;

    //////////////////////////
    // Kneser-Nachbarn
    //////////////////////////

    for B in HyperebenenVertreter do
        SC := B * genmat;
        C0Code := LinearCode(SC);
        CC := Dual(C0Code);

        Bneu := [];
        Raum := sub<CC | C0Code>;

        for i in [1..Dimension(CC)] do
            if not CC.i in Raum then
                Append(~Bneu, CC.i);
                Raum := sub<CC | Raum, CC.i>;
            end if;
        end for;

        // Nach Kapitel 4 ist C0 eine Hyperebene eines selbstdualen Codes.
        // Daher muss C0^perp/C0 genau Dimension 2 besitzen. Ein Abweichen
        // weist auf einen Programm- oder Eingabefehler hin und darf nicht
        // durch das Ueberspringen dieser Hyperebene verborgen werden.
        if #Bneu ne 2 then
            error "Der Quotient C0^perp/C0 hat nicht Dimension 2.";
        end if;

        // Der Quotient C0Code^perp/C0Code hat Dimension 2. Ueber GF(2)
        // besitzt er genau drei eindimensionale Unterraeume. Alle drei
        // muessen untersucht werden; einer davon liefert Caktuell selbst,
        // die beiden anderen liefern die beiden moeglichen Nachbarn.
        Erweiterungsvektoren := [
            Bneu[1],
            Bneu[2],
            Bneu[1] + Bneu[2]
        ];

        for b in Erweiterungsvektoren do
            C1 := sub<CC | C0Code, b>;

            if C1 eq Caktuell then
                continue;
            end if;

            if IsSelfDual(C1) and IsDoublyEven(C1) then
                Neu := true;
                X1 := CodeGraph(C1);

                for i in [1..#Genus] do
                    if IsIsomorphicUnderBlocks(X1, CodeGraphen[i]) then
                        Neu := false;
                        break;
                    end if;
                end for;

                if Neu then
                    Append(~Genus, C1);
                    Append(~CodeGraphen, X1);
                end if;
            end if;
        end for;
    end for;

    // Genau eine Statusmeldung, nachdem die aktuelle Klasse vollstaendig
    // bearbeitet wurde.
    printf "Aktuelle Klasse %o | Gefundene Klassen %o\n", k, #Genus;
    k +:= 1;
end while;

//////////////////////////////
// Ausgabe
//////////////////////////////

print "==================================================";
printf "Laenge: (%o,%o)\n", N1, N2;
printf "Anzahl der Aequivalenzklassen: %o\n", #Genus;
print "==================================================";

for i in [1..#Genus] do
    printf "\nKlasse %o\n", i;
    print "Vertreter (Generatormatrix):";
    print GeneratorMatrix(Genus[i]);

    printf "Stabilisator in S_%o x S_%o:\n", N1, N2;
    print Stabilisatoren[i];
    printf "Ordnung des Stabilisators: %o\n", #Stabilisatoren[i];
    printf "Abstrakter Gruppentyp: %o\n", GroupName(Stabilisatoren[i]);

    if CanIdentifyGroup(#Stabilisatoren[i]) then
        printf "SmallGroup-Kennung: %o\n", IdentifyGroup(Stabilisatoren[i]);
    end if;

    print "Erzeuger des Stabilisators:";
    print Generators(Stabilisatoren[i]);
end for;

//////////////////////////////
// Haeufigkeiten der Stabilisatoren nach SmallGroup-Kennung
//////////////////////////////

SmallGroupKennungen := [];
SmallGroupHaeufigkeiten := [];
SmallGroupVertreter := [];
NichtIdentifizierbar := [];

for i in [1..#Stabilisatoren] do
    S := Stabilisatoren[i];

    if CanIdentifyGroup(#S) then
        Kennung := IdentifyGroup(S);
        PositionKennung := Position(SmallGroupKennungen, Kennung);

        if PositionKennung eq 0 then
            Append(~SmallGroupKennungen, Kennung);
            Append(~SmallGroupHaeufigkeiten, 1);
            // Der Stabilisator der ersten zugehoerigen Klasse dient als
            // konkreter Vertreter dieses abstrakten Isomorphietyps.
            Append(~SmallGroupVertreter, i);
        else
            SmallGroupHaeufigkeiten[PositionKennung] +:= 1;
        end if;
    else
        Append(~NichtIdentifizierbar, i);
    end if;
end for;

print "\nHaeufigkeiten der Stabilisatortypen nach SmallGroup-Kennung:";

for j in [1..#SmallGroupKennungen] do
    VertreterIndex := SmallGroupVertreter[j];
    S := Stabilisatoren[VertreterIndex];

    print "--------------------------------------------------";
    printf "SmallGroup-Kennung: %o\n", SmallGroupKennungen[j];
    printf "Abstrakter Gruppentyp: %o\n", GroupName(S);
    printf "Ordnung: %o\n", #S;
    printf "Haeufigkeit: %o\n", SmallGroupHaeufigkeiten[j];
    printf "Vertreter des Stabilisator-Isomorphietyps (Stabilisator von Klasse %o):\n",
        VertreterIndex;
    print S;
    print "Erzeuger dieses Stabilisatorvertreters:";
    print Generators(S);
end for;

if #NichtIdentifizierbar gt 0 then
    print "--------------------------------------------------";
    print "Nicht durch die SmallGroups-Datenbank identifizierbare Stabilisatoren:";

    for i in NichtIdentifizierbar do
        printf "Klasse %o: Ordnung %o, GroupName %o\n",
            i, #Stabilisatoren[i], GroupName(Stabilisatoren[i]);
    end for;
end if;

//////////////////////////////
// Sicherheitscheck
//////////////////////////////

print "\nIsomorphie-Sicherheitscheck:";
Fehler := false;

for i in [1..#Genus] do
    for j in [i+1..#Genus] do
        if IsIsomorphicUnderBlocks(CodeGraphen[i], CodeGraphen[j]) then
            printf "FEHLER: Klassen %o und %o sind isomorph.\n", i, j;
            Fehler := true;
        end if;
    end for;
end for;

if not Fehler then
    print "Keine doppelten Klassen gefunden.";
end if;

//////////////////////////////
// Massformel fuer den binaeren Typ 2_II und Laenge (N1,N2)
//////////////////////////////

Q := Rationals();
n := (N1 + N2) div 2;

// Anzahl t_(N1,N2) aller gelabelten Codes des Typs 2_II^(N1,N2):
// t_(N1,N2) = Produkt_{j=0}^{n-2} (2^j + 1),
// wobei n = (N1+N2)/2 ist.
tN1N2 := 1;
if n ge 2 then
    for j in [0..n-2] do
        tN1N2 *:= 2^j + 1;
    end for;
end if;

ErwarteteMasse :=
    (Q!tN1N2)/(Q!(Factorial(N1) * Factorial(N2)));
BerechneteMasse := Q!0;

for S in Stabilisatoren do
    BerechneteMasse +:= (Q!1)/(Q!#S);
end for;

print "\nMassformel-Kontrolle:";
printf "Anzahl t_(%o,%o) aller gelabelten Codes: %o\n",
    N1, N2, tN1N2;
printf "Erwartete Masse t_(N1,N2)/(N1!*N2!): %o\n",
    ErwarteteMasse;
printf "Berechnete Masse Summe 1/|Aut(C)|: %o\n",
    BerechneteMasse;

if BerechneteMasse eq ErwarteteMasse then
    print "Massformel erfuellt: Die Klassenliste ist vollstaendig.";
else
    print "WARNUNG: Massformel nicht erfuellt.";
    print "Es fehlen Klassen oder ein Stabilisator wurde nicht korrekt bestimmt.";
end if;

return Genus, Stabilisatoren;

end function;

// Bequemer Spezialfall fuer die diagonale Laenge (n,n).
KneserBlockCodes := function(n)
    return KneserBlockCodesGeneral(n, n);
end function;
