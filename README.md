# Conjugate Weight Enumerators of Self-Dual Codes

This repository contains the computational material accompanying my Master's thesis

**“Konjugierte Gewichtszähler selbstdualer Codes”**

The thesis studies conjugate complete weight enumerators of self-dual binary codes, their invariant-theoretic properties, and related classification problems.

The computations were carried out using **MAGMA**.

## Contents

The repository contains MAGMA programs and computational results used in the thesis. In particular, it includes code for

* computing and classifying self-dual Type II codes,
* applying Kneser's neighbor method,
* determining equivalence classes under the action of
  ($S_{N_1} \times S_{N_2}$),
* computing stabilizer groups,
* checking the corresponding mass formula,
* analyzing standard constructions of self-dual codes,
* and computing conjugate complete weight enumerators.

## Repository Structure

```text
magma/
    MAGMA programs used for the computations

Results/
    Computational results and output files

README.md
```

The individual MAGMA files contain additional comments describing their purpose and usage.

## Programs and Corresponding Results

All paths below are relative to the repository root.

| Program | Purpose | Thesis reference | Example result file |
|---|---|---|---|
| `magma/kneser_blockcodes.m` | Classification, stabilizers, and mass-formula checks | Chapters 4 and 5 | `Results/Kneser_N_N/Ergebnisse_11_11.txt` |
| `magma/Standardkonstruktionen_v2.m` | Construction analysis for equal block lengths | Chapter 5 | `Results/Constructions_N_N/Standardkonstruktionen_11_11.txt` |
| `magma/Standardkonstruktionen_ungleich.m` | Structural analysis for unequal block lengths | Chapter 5 | `Results/Constructions_N1_N2/Konstruktionen_ungleich_6_14.txt` |
| `magma/CCWE_Genus2.m` | Conjugate complete weight enumerators | Section 5.1 | `Results/Harmonic_Polynomials/Beispiel_5_1_Ausgabe.txt` |
| `magma/Harmonische_Projektion.m` | Laplace operator and harmonic invariants | Section 5.1 | `Results/Harmonic_Polynomials/Beispiel_5_1_Ausgabe.txt` |

## Requirements

The computations require the computer algebra system

**MAGMA**

The programs were developed and tested with MAGMA during the preparation of the Master's thesis.

## Running the Programs

Run the following commands in a MAGMA session with the repository root as the working directory. This is the directory containing `magma/`, `Results/`, and `README.md`.

### 1. Small Classification Example

The case $(N_1,N_2)=(6,6)$ provides a small example for checking the classification program. The Type II condition is understood as defined in the thesis.

```magma
load "kneser_blockcodes.m";

Klassen, Stabilisatoren := KneserBlockCodes(6);
```

**Expected result:** four equivalence classes under $S_6 \times S_6$.

The output contains representatives, stabilizer information, and the mass-formula check. The computed mass should agree with the theoretical value.

**Result file:** `Results/Kneser_N_N/Ergebnisse_6_6.txt`

### 2. Construction Analysis for Equal Block Lengths

The following example analyzes the codes of length $(11,11)$ considered in Chapter 5.

```magma
load "Standardkonstruktionen_v2.m";


// Prepare the representatives in the format required by the analysis.
Klassen, Stabilisatoren := KneserBlockCodes(11);

// Analyze the representatives
AnalyseAllStandardConstructions(Klassen, 11, 11);
```

The output records the construction or structural description assigned to each analyzed class.

**Result file:** `Results/Constructions_N_N/Standardkonstruktionen_11_11.txt`

### 3. Weight Enumerators and Harmonic Invariants

This example reproduces the polynomial computations in Section 5.1 for genus $m=2$ and bidegree $(6,6)$.

The files CCWE_Genus2.m and Harmonische_Projektion.m define the conjugate complete weight enumerators, the Laplace operator, and the auxiliary functions required for the linear systems. The file PrintHarmonicCalculation.m contains the procedure used to write the individual calculations to the output file.

The computation is started with

```magma
load "CCWE_Genus2.m";
load "Harmonische_Projektion.m";
load "PrintHarmonicCalculation.m";

PrintHarmonicCalculation(F1, 1, S, W44);
PrintHarmonicCalculation(F2poly, 2, S, W44);
PrintHarmonicCalculation(F3, 3, S, W44);
PrintHarmonicCalculation(F4, 4, S, W44);
```

The output documents

- the polynomials $F_1,\ldots,F_4$,
- their images under the Laplace operator,
- the linear system used to determine the harmonic invariants,
- and its solutions.

**Expected result:** the harmonic parts of $F_1$ and $F_2$ vanish. The harmonic parts of $F_3$ and $F_4$ are linearly independent and are annihilated by the Laplace operator.

**Result file:** `Results/Harmonic_Polynomials/Beispiel_5_1_Ausgabe.txt`

### 4. Structural Analysis for Unequal Block Lengths

The following example analyzes the codes of length $(6,14)$ considered in Chapter 5.

```magma
load "Standardkonstruktionen_ungleich.m";

// Prepare the representatives in the format required by the analysis.
Klassen, Stabilisatoren := KneserBlockCodesGeneral(6,14);

// Analyze the representatives 
AnalyseAllUnequalConstructions(Klassen, 6, 14);
```

The output records decompositions, projections, projection kernels, and quotient-space descriptions.

**Result file:** `Results/Constructions_N1_N2/Konstruktionen_ungleich_6_14.txt`

## Reproducibility

The files in this repository are intended to make the computational results presented in the thesis reproducible.

For some of the larger parameters, in particular the classification computations using Kneser's neighbor method, the running time can be substantial.

## Master's Thesis

**Title:** *Konjugierte Gewichtszähler selbstdualer Codes*
**Author:** Andreas Hild
**Year:** 2026

The notation and mathematical background used in the programs are explained in the thesis.

