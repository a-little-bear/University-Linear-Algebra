# Chapter 72B: Multivariate Statistical Inference

<div class="context-flow" markdown>

**Prerequisites**: Matrix Calculus (Ch47) · Wishart Distribution (Ch72A) · QR Decomposition (Ch10) · SVD (Ch11)

**Chapter Outline**: Parameter Estimation (MLE) → Hotelling’s $T^2$ Test → Multivariate Analysis of Variance (MANOVA) → Canonical Correlation Analysis (CCA) → Structural Equation Modeling (SEM) → High-dimensional Covariance Matrix Estimation → Shrinkage Estimators (Ledoit-Wolf)

**Extension**: Multivariate statistical inference applies the spectral theory of random matrices to decision-making under uncertainty; CCA is the algebraic bridge for finding correlations between two sets of multidimensional variables.

</div>

Multivariate statistical inference utilizes matrix algebra to perform hypothesis testing and parameter estimation on high-dimensional data. This field converts the comparison of group differences into the comparison of matrix spectra (eigenvalues) and the dependencies between variables into the geometry of subspace angles.

---

## 72B.1 Core Estimators and Testing Theory

!!! definition "Definition 72B.1 (MLE of Mean and Covariance)"
    For $\mathcal{N}_p(\mu, \Sigma)$, given $n$ samples, the Maximum Likelihood Estimates are:
    $$\hat{\mu} = \bar{x}, \quad \hat{\Sigma} = \frac{1}{n} \sum (x_i - \bar{x})(x_i - \bar{x})^T$$
    These are the orthogonal projections of the sample cloud onto the parameter space.

!!! theorem "Theorem 72B.3 (Wilks' Lambda and MANOVA)"
    In MANOVA, the test statistic for comparing group means is Wilks' Lambda:
    $$\Lambda = \frac{\det(E)}{\det(E + H)}$$
    where $E$ is the error (within-group) sum of squares matrix and $H$ is the hypothesis (between-group) sum of squares matrix. This is a function of the eigenvalues of $E^{-1}H$.

---

## Exercises

1. **[Group Comparison] Why does MANOVA use matrix determinants instead of just summing variances?**
   ??? success "Solution"
       The determinant $\det(E)$ represents the "generalized variance" or volume of the sample cluster. MANOVA accounts for the correlation between variables; summing variances would ignore the covariance structure and lead to incorrect inference when variables are non-orthogonal.

2. **[Hotelling's T-squared] Show that Hotelling’s $T^2$ is the multivariate generalization of the squared t-statistic.**
   ??? success "Solution"
       $T^2 = n(\bar{x}-\mu)^T S^{-1}(\bar{x}-\mu)$. In the scalar case ($p=1$), this reduces to $n(\bar{x}-\mu)^2 / s^2 = t^2$. The use of the matrix inverse $S^{-1}$ (Mahalanobis distance) standardizes the distance across all correlated dimensions.

3. **[CCA Geometry] In Canonical Correlation Analysis (CCA), how do the canonical correlations relate to the angles between two subspaces?**
   ??? success "Solution"
       The canonical correlations are the cosines of the **principal angles** between the subspace spanned by variables $X$ and the subspace spanned by variables $Y$. Maximizing the correlation is equivalent to finding the vectors in each subspace that are closest in an angular sense.

4. **[Calculation] Given $E = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ and $H = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$. Calculate Wilks' $\Lambda$.**
   ??? success "Solution"
       $E+H = \begin{pmatrix} 3 & 1 \\ 1 & 2 \end{pmatrix}$.
       $\det(E) = 4-1 = 3$; $\det(E+H) = 6-1 = 5$.
       $\Lambda = 3/5 = 0.6$. A smaller $\Lambda$ would indicate more evidence against the null hypothesis.

5. **[Discriminant Analysis] Relate Fisher’s Linear Discriminant to the generalized eigenvalue problem.**
   ??? success "Solution"
       Fisher seeks a projection vector $w$ that maximizes $J(w) = \frac{w^T H w}{w^T E w}$. Differentiating leads to the generalized eigenvalue equation $Hw = \lambda Ew$. The optimal projection is the eigenvector corresponding to the largest eigenvalue of $E^{-1}H$.

6. **[Shrinkage Estimation] Why is the MLE estimator $\hat{\Sigma}$ poorly conditioned in high dimensions ($p \approx n$)?**
   ??? success "Solution"
       As $p/n$ increases, the sample eigenvalues become more dispersed (Marchenko-Pastur law). The smallest eigenvalues are biased towards zero, making $S^{-1}$ unstable. Shrinkage estimators like Ledoit-Wolf use a convex combination $S(\lambda) = (1-\lambda)S + \lambda I$ to pull eigenvalues away from the boundaries.

7. **[SVD and Regression] Describe the role of SVD in Principal Component Regression (PCR).**
   ??? success "Solution"
       PCR first performs SVD on the design matrix $X = U\Sigma V^T$, then regresses the response $y$ on the first $k$ columns of $U$. This eliminates collinearity by using the orthogonal principal components as predictors.

8. **[Invariance] Is Hotelling’s $T^2$ statistic invariant under non-singular linear transformations $x \mapsto Ax + b$?**
   ??? success "Solution"
       Yes. Substituting the transformed mean and covariance into the $T^2$ formula results in the cancellation of the matrices $A$ and $A^{-1}$, proving that the test result is independent of the choice of coordinate system.

9. **[Partial Correlation] Explain the relationship between the precision matrix $\Omega = \Sigma^{-1}$ and partial correlations.**
   ??? success "Solution"
       The $(i,j)$ entry of the inverse covariance matrix, after normalization, is equal to the negative of the partial correlation between variables $i$ and $j$ given all other variables. Zero entries in $\Omega$ indicate conditional independence (Gaussian graphical models).

10. **[Spectral Hypothesis] How does the distribution of the largest eigenvalue $\lambda_{\max}(E^{-1}H)$ relate to the Roy’s Largest Root test?**
    ??? success "Solution"
        Roy’s test uses only the largest eigenvalue as the test statistic. It is the most powerful test when the group differences are concentrated along a single dimension (rank-1 alternative), reflecting the sensitivity of the spectral radius to structured perturbations.

## Chapter Summary

This chapter applies matrix theory to the logic of statistical decision-making:

1. **Spectral Comparison**: Established group difference testing as the analysis of relative matrix spectra ($E^{-1}H$).
2. **Subspace Alignment**: Formulated CCA as the problem of minimizing angles between variable subspaces via SVD.
3. **Regularized Estimation**: Introduced shrinkage methods to solve the numerical instability of sample matrices in high-dimensional settings.
4. **Distance Geometry**: Utilized Mahalanobis metrics and quadratic forms to define robust multivariate hypothesis tests.
