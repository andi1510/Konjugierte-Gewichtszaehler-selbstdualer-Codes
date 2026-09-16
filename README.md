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
| `magma/kneser_blockcodes.m` | Classification, stabilizers, and mass-formula checks | Chapters 4 and 5 | `Results/Harmonic_Polynomials/Beispiel_5_1_Ausgabe.txt` |
| `magma/Standardkonstruktionen_v2.m` | Construction analysis for equal block lengths | Chapter 5 | `Results/Constructions_N_N/Standardkonstruktionen_11_11.txt` |
| `magma/Stamdardlpmstruktionen_ungleich.m` | Structural analysis for unequal block lengths | Chapter 5 | `Results/Constructions_N1_N2/Konstruktionen_ungleich_6_14.txt` |
| `magma/CCWE_Genus2.m` | Conjugate complete weight enumerators | Section 5.1 | `Results/Harmonic_Polynomials/Beispiel_5_1_Ausgabe.txt` |
| `magma/Harmonische_Projektion.m` | Laplace operator and harmonic invariants | Section 5.1 | `Results/Kneser_N_N/Ergebnisse_11_11.txt` |

## Requirements

The computations require the computer algebra system

**MAGMA**

The programs were developed and tested with MAGMA during the preparation of the Master's thesis.

## Reproducibility

The files in this repository are intended to make the computational results presented in the thesis reproducible.

For some of the larger parameters, in particular the classification computations using Kneser's neighbor method, the running time can be substantial.

## Master's Thesis

**Title:** *Konjugierte Gewichtszähler selbstdualer Codes*
**Author:** Andreas Hild
**Year:** 2026

The notation and mathematical background used in the programs are explained in the thesis.

