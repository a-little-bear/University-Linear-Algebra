# Chapter 72B: Multivariate Statistical Inference

<div class="context-flow" markdown>

**Prerequisites**: Matrix Distributions (Ch72A) · Matrix Analysis (Ch14) · Quadratic Forms (Ch09) · Linear Equations (Ch01)

**Chapter Outline**: From Univariate to Multivariate Inference → Hotelling's $T^2$ Distribution (Matrix Version of Student's t) → Wilks' Lambda Distribution (Ratio of Determinants) → Multivariate Analysis of Variance (MANOVA) → Likelihood Ratio Tests (LRT) → Testing Equality of Covariance Matrices → Canonical Correlation Analysis (CCA) → Multivariate Regression Inference → Challenges of High-dimensional Inference ($p > n$) → Applications: Social Sciences, Psychometrics, and Medical Clinical Trials

**Extension**: Multivariate inference is the highest application of linear algebra in decision theory; by investigating the weighted sums or products of eigenvalues, it determines whether observed differences in high-dimensional space stem from true effects or random noise.

</div>

After establishing the matrix-based description of data (Ch72A), the core task of statistics becomes: how to make rigorous inferences based on observed matrix samples. **Multivariate Statistical Inference** utilizes the trace, determinant, and eigenvalues of matrices to construct test statistics. It answers questions such as "Do two groups of multidimensional data differ significantly?" or "Is there a latent link between two sets of variables?"

---

## 72B.1 Hotelling's $T^2$ Test

!!! definition "Definition 72B.1 (Hotelling's $T^2$ Statistic)"
    To test whether a mean vector $\boldsymbol{\mu}$ equals a hypothesized $\boldsymbol{\mu}_0$, we define the statistic:
    $$T^2 = n (\bar{\mathbf{x}} - \boldsymbol{\mu}_0)^T \mathbf{S}^{-1} (\bar{\mathbf{x}} - \boldsymbol{\mu}_0)$$
    where $\mathbf{S}$ is the sample covariance matrix.
    **Algebraic Essence**: This is the square of the Mahalanobis distance, which utilizes the inverse of the covariance matrix to normalize fluctuations across different dimensions.

---

## 72B.2 Multivariate Analysis of Variance (MANOVA)

!!! technique "Matrix Decomposition Perspective"
    In MANOVA, we decompose the total sum of squares and cross-products matrix $\mathbf{T}$ into the within-group error matrix $\mathbf{E}$ and the between-group effect matrix $\mathbf{H}$:
    $$\mathbf{T} = \mathbf{H} + \mathbf{E}$$
    **Statistics**:
    - **Wilks' Lambda**: $\Lambda = \frac{\det(\mathbf{E})}{\det(\mathbf{H} + \mathbf{E})}$. Its distribution is determined by the product of eigenvalues.
    - **Pillai’s Trace**: Based on $\operatorname{tr}(\mathbf{H}(\mathbf{H}+\mathbf{E})^{-1})$.

---

## 72B.3 Canonical Correlation Analysis (CCA)

!!! definition "Definition 72B.2 (CCA)"
    Given two sets of variables $\mathbf{x}$ and $\mathbf{y}$, find projection vectors $\mathbf{a}, \mathbf{b}$ such that the correlation between the linear combinations $u = \mathbf{a}^T \mathbf{x}$ and $v = \mathbf{b}^T \mathbf{y}$ is maximized.
    **Solution**: This is equivalent to solving a generalized eigenvalue problem involving the cross-covariance matrices.

---

## 72B.4 Challenges of High-dimensional Inference ($p > n$)

!!! warning "Curse of Dimensionality"
    When the number of variables $p$ exceeds the sample size $n$, the sample covariance matrix $\mathbf{S}$ is **singular** (not invertible), causing the $T^2$ statistic to fail.
    **Strategies**: Use ridge regularization (shrinkage estimators) or projection-based non-parametric methods.

---

## Exercises

1.  **[Basics] What does Hotelling's $T^2$ reduce to when the number of variables $p=1$?**
    ??? success "Solution"
        It reduces to the square of the univariate $t$-statistic: $t^2 = \frac{n(\bar{x}-\mu_0)^2}{s^2}$.

2.  **[Determinant] Why does Wilks' Lambda use the ratio of determinants?**
    ??? success "Solution"
        The determinant represents "generalized variance" (volume) in high-dimensional space. The ratio reflects the proportion of residual volume to total fluctuation volume; a smaller ratio indicates more significant differences between groups.

3.  **[Calculation] If $\mathbf{H} = \operatorname{diag}(10, 0)$ and $\mathbf{E} = \operatorname{diag}(1, 1)$, calculate Wilks' Lambda.**
    ??? success "Solution"
        $\det(\mathbf{E}) = 1$. $\det(\mathbf{H}+\mathbf{E}) = \det(\operatorname{diag}(11, 1)) = 11$.
        $\Lambda = 1/11 \approx 0.09$.

4.  **[Eigenvalue] Prove: Wilks' Lambda can be expressed as $\prod \frac{1}{1+\lambda_i}$, where $\lambda_i$ are the eigenvalues of $\mathbf{E}^{-1}\mathbf{H}$.**
    ??? success "Solution"
        $\frac{\det(\mathbf{E})}{\det(\mathbf{H}+\mathbf{E})} = \det(\mathbf{E}(\mathbf{H}+\mathbf{E})^{-1}) = \det(I + \mathbf{E}^{-1}\mathbf{H})^{-1}$. The result follows from the property of determinants.

5.  **[Application] In what scenario is Pillai's trace preferred over Wilks' Lambda?**
    ??? success "Solution"
        Pillai's trace is generally more robust when the assumption of equal covariance matrices (homoscedasticity) is violated.

6.  **[CCA] The maximum canonical correlation coefficient in CCA is related to which matrix property?**
    ??? success "Solution"
        It is related to the singular values of the cross-covariance matrix between the two groups (after variance normalization).

7.  **[Invariance] Prove that the $T^2$ statistic is invariant under linear coordinate transformations.**
    ??? success "Solution"
        If $\mathbf{x} \to \mathbf{Ax}$, then the mean becomes $\mathbf{A}\bar{\mathbf{x}}$ and the covariance becomes $\mathbf{ASA}^T$. Substituting these into the formula, the $\mathbf{A}$ matrices cancel out.

8.  **[Detection] What information does Box's M test utilize to check for the equality of two covariance matrices?**
    ??? success "Solution"
        It uses the weighted difference between the log-determinants of the sample covariance matrices.

9.  **[Regression] In multivariate regression, what does the trace of the error matrix $\operatorname{tr}(\mathbf{E})$ correspond to?**
    ??? success "Solution"
        It corresponds to the total sum of squared residuals (Total SSE) across all dependent variables.

10. **[Limits] What distribution does the $T^2$ distribution approach in large samples?**

   ??? success "Solution"
        It approaches a $\chi^2$ distribution with $p$ degrees of freedom.

## Chapter Summary

Multivariate statistical inference is the ultimate algebraic judgment in empirical science:

1.  **Algebraization of Distance**: Hotelling's $T^2$ proves that through matrix inversion, we can correct complex anisotropic fluctuations into standard statistical distances, establishing benchmarks for multi-dimensional difference testing.
2.  **Competition of Volumes**: Wilks' Lambda reduces complex group comparisons to a "game of determinants" (volumes), revealing the algebraic share of explanatory power within multi-dimensional space.
3.  **Deconstruction of Correlation**: CCA demonstrates how to utilize generalized eigenvalue problems to extract resonant linear signals from two seemingly chaotic data streams, providing a mathematical scalpel for understanding coupling mechanisms in complex systems.
