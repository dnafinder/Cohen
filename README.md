[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/Cohen&file=kappa.m)

# Cohen – Cohen's kappa in MATLAB (kappa.m)

## 📌 Overview

This repository provides a MATLAB implementation of **Cohen's kappa** coefficient for inter-rater agreement, with support for:

- Unweighted kappa (nominal categories)
- Linearly weighted kappa (ordinal categories)
- Quadratically weighted kappa (ordinal categories)
- Custom user-defined weight matrices

The main function is:

- kappa.m

It computes Cohen's kappa, its standard error, confidence interval, and several useful descriptive statistics and diagnostic measures.

---

## ✨ Features

- Supports **unweighted, linear, quadratic, and custom** weighting schemes
- Accepts any **square contingency table** of nonnegative integer frequencies
- Returns:
  - Kappa estimate
  - Confidence interval for kappa
  - Detailed statistics in a structured output
- Implements **Landis & Koch** qualitative classification of agreement
- Performs **z-test** for the null hypothesis of purely accidental agreement
- Optional **textual output** for quick inspection in the Command Window
- Robust input validation and basic numerical stability checks

---

## 📥 Function syntax

Basic calls:

- kappa(X)
- kappa(X, W)
- kappa(X, W, ALPHA)
- kappa(X, W, ALPHA, DISPLAY)
- [K, CI] = kappa(...)
- [K, CI, STATS] = kappa(...)

Where:

- X  
  Square M-by-M matrix of **nonnegative integer frequencies**, representing the cross-classification of two raters (or methods) on M mutually exclusive categories.

- W  
  Weight specification:
  - Scalar:
    - 0 → unweighted (default)
    - 1 → linearly weighted
    - 2 → quadratically weighted
  - Matrix:
    - M-by-M custom weight matrix (for M categories).

- ALPHA  
  Significance level for the confidence interval (default 0.05).

- DISPLAY  
  Logical flag (true/false).  
  - true (default) → prints a formatted report in the Command Window  
  - false → silent mode (only outputs are returned)

Outputs:

- K  
  Cohen's kappa estimate.

- CI  
  1-by-2 vector with the lower and upper bounds of the confidence interval for kappa.

- STATS  
  Structure containing auxiliary results:
  - po                     – observed agreement
  - pe                     – expected agreement by chance
  - trueAgreement          – po - pe
  - residualAgreement      – 1 - pe
  - kappa                  – kappa estimate
  - se                     – standard error of kappa
  - ci                     – confidence interval for kappa
  - alpha                  – significance level
  - z                      – z statistic
  - p                      – p-value
  - kappaMax               – maximum possible kappa given the margins
  - kappaRatio             – kappa / kappaMax
  - weightsType            – 'unweighted', 'linear', 'quadratic', or 'custom'
  - weightMatrix           – weight matrix actually used
  - n                      – total number of observations
  - m                      – number of categories
  - maxObservableAgreement – maximum observable agreement (pom)
  - landisKochClass        – qualitative classification of agreement

---

## 📊 Example

Simple example with a 3×3 contingency table:

x = [88 14 18; ...
     10 40 10; ...
      2  6 12];

% Unweighted kappa, with printed report
kappa(x);

% Linear weights, 95% CI, silent mode with full stats
[k, ci, stats] = kappa(x, 1, 0.05, false);

The first call prints a complete report (observed and expected agreement, kappa, standard error, confidence interval, Landis & Koch class, z statistic, p-value, etc.). The second call suppresses printing and returns all key quantities programmatically.

---

## 🧮 Weighting schemes

- Unweighted (W = 0)  
  Only exact agreement (diagonal) is given full credit (weight = 1); off-diagonal cells have weight 0.

- Linear weights (W = 1)  
  Weights decrease linearly as categories become more distant, assuming an ordinal structure:
  - For categories i and j, the weight is:
    1 - |i - j| / (M - 1)

- Quadratic weights (W = 2)  
  Weights decrease quadratically with distance:
  - For categories i and j, the weight is:
    1 - ((i - j) / (M - 1))^2

- Custom weights (W = M×M matrix)  
  User provides a full M-by-M matrix of weights, typically with:
  - Values in [0, 1]
  - Diagonal elements equal to 1 (perfect agreement)
  This allows for arbitrary weighting structures tailored to specific problems.

---

## 🧾 Landis & Koch benchmarks

The function classifies the resulting kappa according to the widely cited Landis and Koch (1977) benchmarks:

- k < 0           → Poor agreement
- 0.00–0.20       → Slight agreement
- 0.21–0.40       → Fair agreement
- 0.41–0.60       → Moderate agreement
- 0.61–0.80       → Substantial agreement
- 0.81–1.00       → Perfect agreement

The qualitative class is returned in STATS.landisKochClass and printed in the report when DISPLAY is true.

---

## 📚 References

- Cohen J. (1960). A coefficient of agreement for nominal scales. Educational and Psychological Measurement, 20(1), 37–46.
- Landis JR, Koch GG. (1977). The measurement of observer agreement for categorical data. Biometrics, 33(1), 159–174.
- Cardillo G. (2007). Cohen's kappa: compute the Cohen's kappa ratio on a square matrix. Available from:  
  https://github.com/dnafinder/Cohen

---

## 📁 Repository structure

Main files:

- kappa.m  
  Core function implementing Cohen's kappa (unweighted, weighted, and custom weights).

Additional files (if present) may include:

- Examples or demo scripts
- Test scripts
- Documentation or auxiliary utilities

---

## 🚀 Getting started

1. Clone or download the repository:

   - GitHub: https://github.com/dnafinder/Cohen

2. Add the folder containing kappa.m to your MATLAB path, for example:

   - In MATLAB: Home → Set Path → Add Folder
   - Or programmatically with addpath

3. Call kappa from the Command Window or from your own scripts/functions:

   - kappa(X)
   - [k, ci, stats] = kappa(X, 1, 0.01, false);

---

## 🧪 Testing and validation

To validate your installation and understand the output:

- Run the example in the help section of kappa.m
- Compare with published examples or with other implementations (e.g., statistical software) for known contingency tables

If you use custom weight matrices, always verify that your weights are appropriate for the ordinal or quasi-ordinal structure of your categories.

---

## ⚖️ License and citation

Please refer to the license file provided in the repository (if present) for detailed licensing terms.

If you use this code in scientific work, software, or teaching material, an appropriate citation is:

Cardillo G. (2007). Cohen's kappa: compute the Cohen's kappa ratio on a square matrix. Available from:  
https://github.com/dnafinder/Cohen
