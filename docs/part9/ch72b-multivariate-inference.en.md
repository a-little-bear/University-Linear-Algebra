# Chapter 72B: Multivariate Statistical Inference

<div class="context-flow" markdown>

**Prerequisites**: Matrix Distributions (Ch72A) · Positive Definite Matrices (Ch16) · Generalized Inverses (Ch33) · Projections (Ch07)

**Chapter Outline**: From Univariate to Multivariate → Definition of the Multivariate Normal (MVN) Distribution and its Quadratic Density → Parameter Estimation: Maximum Likelihood Estimation (MLE) of Mean Vectors and Covariance Matrices → Hypothesis Testing: Hotelling’s $T^2$ Statistic (Matrix version of the $t$-test) → Wilks’s Lambda ($\Lambda$) Distribution and Likelihood Ratio Tests → Linear Algebraic Essence of Discriminant Analysis (LDA) → Canonical Correlation Analysis (CCA) as a Generalized Eigenvalue Problem → Applications: Multi-indicator Difference Testing in Clinical Medicine, Psychological Measurement, and Multivariate Control Charts in Quality Control

**Extension**: Multivariate statistical inference is the supreme application of linear algebra in "verifying the truth"; it transforms complex hypothesis tests into interval estimations of operator spectra and matrix quadratic forms. It proves that the reliability of scientific conclusions depends on the geometric stability of statistical operators in high-dimensional space—the mathematical judge of modern empirical research.

</div>

In statistics, when we observe multiple variables simultaneously (e.g., height, weight, blood pressure), analyzing each individually ignores their correlations. **Multivariate Statistical Inference** achieves a unified analysis by integrating observations into vectors and matrices. Using matrix-based statistics like **Hotelling’s $T^2$** and **Wilks’s Lambda**, we can determine whether group differences are essential within a probabilistic framework. This chapter introduces this algebraic inference framework serving as the bedrock of scientific evidence.

---

## 72B.1 Multivariate Normal Distribution and MLE

!!! definition "Definition 72B.1 (MVN Distribution)"
    A vector $\mathbf{x} \in \mathbb{R}^p$ follows a multivariate normal distribution $N(\boldsymbol{\mu}, \Sigma)$ if its probability density function is:
    $$f(\mathbf{x}) = \frac{1}{(2\pi)^{p/2} |\Sigma|^{1/2}} \exp \left( -\frac{1}{2} (\mathbf{x}-\boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x}-\boldsymbol{\mu}) \right)$$
    where the exponent is a **positive-definite quadratic form**, defining hyper-ellipsoidal contours of equal probability.

---

## 72B.2 Hotelling’s $T^2$ Statistic

!!! theorem "Theorem 72B.1 (Mean Testing)"
    To test $H_0: \boldsymbol{\mu} = \boldsymbol{\mu}_0$, we construct the $T^2$ statistic:
    $$T^2 = n (\bar{\mathbf{x}} - \boldsymbol{\mu}_0)^T S^{-1} (\bar{\mathbf{x}} - \boldsymbol{\mu}_0)$$
    where $S$ is the sample covariance matrix.
    **Algebraic Essence**: This is the square of the **Mahalanobis Distance**, quantifying the statistical distance of the sample center from the target.

---

## 72B.3 Canonical Correlation Analysis (CCA)

!!! technique "Technique: Generalized Eigenvalue Problem"
    CCA seeks the maximum correlation between two sets of variables $X$ and $Y$. This is equivalent to solving a generalized eigenvalue problem involving cross-covariance matrices:
    $$\Sigma_{XY} \Sigma_{YY}^{-1} \Sigma_{YX} \mathbf{a} = \rho^2 \Sigma_{XX} \mathbf{a}$$
    This demonstrates how linear algebra uncovers deep consistency across different dimensions through operator composition.

---

## Exercises

**1. [Basics] Let $\mathbf{x} \sim N(\boldsymbol{\mu}, \Sigma)$. Find the distribution of the linear transform $\mathbf{y} = A\mathbf{x} + \mathbf{b}$.**

??? success "Solution"
    **Using Linearity of Expectation and Variance:**
    1. $E[\mathbf{y}] = A E[\mathbf{x}] + \mathbf{b} = A\boldsymbol{\mu} + \mathbf{b}$.
    2. $\operatorname{Var}(\mathbf{y}) = A \operatorname{Var}(\mathbf{x}) A^T = A\Sigma A^T$.
    **Conclusion**: $\mathbf{y} \sim N(A\boldsymbol{\mu} + \mathbf{b}, A\Sigma A^T)$. This proves normality is preserved under affine transformations.

**2. [MLE] Prove that the MLE $\hat{\boldsymbol{\mu}}$ for the multivariate mean is the sample mean $\bar{\mathbf{x}}$.**

??? success "Solution"
    **Matrix Differentiation:**
    1. The log-likelihood involves the term $-\sum (\mathbf{x}_i - \boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x}_i - \boldsymbol{\mu})$.
    2. Differentiate with respect to $\boldsymbol{\mu}$ (using Ch47A formulas): $\sum 2\Sigma^{-1} (\mathbf{x}_i - \boldsymbol{\mu}) = 0$.
    3. Since $\Sigma^{-1}$ is non-singular, $\sum \mathbf{x}_i - n\boldsymbol{\mu} = 0$.
    **Conclusion**: $\hat{\boldsymbol{\mu}} = \frac{1}{n} \sum \mathbf{x}_i$.

**3. [Calculation] For $n=100, p=2$, if the sample mean difference $\mathbf{d} = (1, 1)^T$ and covariance $S = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$, compute $T^2$.**

??? success "Solution"
    **Steps:**
    1. $S^{-1} = \frac{1}{0.75} \begin{pmatrix} 1 & -0.5 \\ -0.5 & 1 \end{pmatrix} = \begin{pmatrix} 4/3 & -2/3 \\ -2/3 & 4/3 \end{pmatrix}$.
    2. $\mathbf{d}^T S^{-1} \mathbf{d} = (1, 1) \begin{pmatrix} 2/3 \\ 2/3 \end{pmatrix} = 4/3$.
    3. $T^2 = 100 \cdot (4/3) \approx 133.3$.
    **Conclusion**: Since $T^2$ is large, we reject the null hypothesis and conclude a significant difference in means.

**4. [Wishart] What distribution does the sample covariance matrix $S$ follow?**

??? success "Solution"
    **Conclusion: The Wishart Distribution.**
    Specifically, $(n-1)S \sim W_p(n-1, \Sigma)$. This is the natural matrix-space generalization of the scalar Chi-squared distribution, forming the sampling distribution foundation for multivariate inference.

**5. [LDA] What equation defines the optimal direction $\mathbf{w}$ in Fisher Linear Discriminant Analysis?**

??? success "Solution"
    **Conclusion: The generalized eigenvalue equation $S_B \mathbf{w} = \lambda S_W \mathbf{w}$.**
    Where $S_B$ is the between-class scatter and $S_W$ is the within-class scatter. The goal is to maximize the Rayleigh quotient $J(w) = \frac{w^T S_B w}{w^T S_W w}$, a classic problem of finding the direction of maximum scaling.

**6. [Properties] Prove that the contours of an MVN density are hyper-ellipsoids.**

??? success "Solution"
    **Reasoning:**
    1. Constant density requires the exponent $(\mathbf{x}-\boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x}-\boldsymbol{\mu}) = c$.
    2. Since $\Sigma$ is PD, $\Sigma^{-1}$ is PD.
    3. In analytic geometry, the set defined by a PD quadratic form equal to a constant is a **hyper-ellipsoid** with axes determined by the eigenvectors.

**7. [Independence] What does it imply if the covariance matrix $\Sigma$ is diagonal?**

??? success "Solution"
    **Conclusion: The variables are mutually independent.**
    For the normal distribution, being uncorrelated (zero covariance) is equivalent to independence. Matrix-wise, this means the joint density factors into the product of marginal densities.

**8. [Limits] Can the $T^2$ statistic be computed if $p > n$?**

??? success "Solution"
    **Conclusion: No (not directly).**
    **Reasoning**: If the sample size is less than the number of variables, the sample covariance matrix $S$ is singular ($\operatorname{rank} \le n-1 < p$). Thus $S^{-1}$ does not exist. One must use **generalized inverses** (Ch33) or regularization.

**9. [Wilks] What is Wilks’s Lambda ($\Lambda$)?**

??? success "Solution"
    **Definition:**
    $\Lambda = \frac{|S_W|}{|S_W + S_B|}$. It is the core metric in Multivariate ANOVA (MANOVA).
    **Algebraic meaning**: It is the ratio of the "unexplained residual volume" to the "total deviation volume." A smaller $\Lambda$ implies stronger group separation.

**10. [Application] Briefly state how Principal Component Regression (PCR) solves multi-collinearity.**

??? success "Solution"
    1. If features $X$ are multi-collinear, $X^T X$ is near singular.
    2. Extract the first $k$ principal components $Z = U_k \Sigma_k$ via SVD.
    3. Perform regression on $Z$.
    **Conclusion**: Since the columns of $Z$ are orthogonal, the new gram matrix $Z^T Z = \Sigma_k^2$ is diagonal and perfectly conditioned, eliminating numerical instability.

## Chapter Summary

Multivariate statistical inference is the "rule of truth" for empirical science in linear algebra:

1.  **Geometric Evidence**: It quantifies statistical differences as distances (Hotelling) and volume ratios (Wilks) in vector space, providing rigorous geometric criteria for scientific observations.
2.  **Operator Insight**: Through generalized eigenvalue problems, CCA and LDA reveal core correlation modes hidden behind massive data noise, showcasing the feature extraction power of linear algebra.
3.  **Completeness of Distributions**: From MVN to Wishart matrices, this chapter establishes the algebraic order of the high-dimensional stochastic world, proving that statistical inference is an extension of operator properties under probability measures.
