# Chapter 72A: Matrix-valued Distributions

<div class="context-flow" markdown>

**Prerequisites**: Matrix Calculus (Ch47) · Kronecker Product (Ch19) · Random Matrices (Ch23) · Probability Theory

**Chapter Outline**: Multivariate Normal Distribution → Matrix Normal Distribution → Wishart Distribution (Sample Covariance) → Matrix Beta and Gamma Distributions → Jacobian of Matrix Transformations → Matrix Variate T-distribution → Characteristic Functions of Matrix Distributions

**Extension**: The Wishart distribution is the algebraic foundation of Multivariate Statistical Analysis (MANOVA, PCA); matrix distributions characterize the uncertainty of high-dimensional estimators.

</div>

Matrix-valued distributions extend random variables from scalars and vectors to matrices. This field utilizes the Kronecker product to describe the internal correlation structure of matrices, mapping the statistical properties of multi-dimensional data to the properties of random linear operators.

---

## 72A.1 Matrix Normal and Wishart Distributions

!!! definition "Definition 72A.1 (Matrix Normal Distribution)"
    A random matrix $X \in \mathbb{R}^{n \times p}$ follows a **Matrix Normal Distribution** $\mathcal{MN}_{n,p}(M, U, V)$ if:
    $$\operatorname{vec}(X) \sim \mathcal{N}_{np}(\operatorname{vec}(M), V \otimes U)$$
    where $U$ captures correlations between rows and $V$ captures correlations between columns.

!!! theorem "Theorem 72A.3 (Wishart Distribution and Sample Covariance)"
    Let $X_1, \dots, X_n$ be i.i.d. samples from $\mathcal{N}_p(0, \Sigma)$. Then the matrix $S = \sum_{i=1}^n X_i X_i^T$ follows a **Wishart Distribution** $W_p(n, \Sigma)$.

---

## Exercises

1. **[Matrix Normal] Explain why the covariance of a Matrix Normal distribution is represented by a Kronecker product $V \otimes U$.**
   ??? success "Solution"
       The Kronecker product $V \otimes U$ encodes a separable correlation structure: $U$ represents the covariance between the $n$ rows (e.g., temporal correlation), and $V$ represents the covariance between the $p$ columns (e.g., spatial or variable correlation). This drastically reduces the number of parameters compared to a general $np \times np$ covariance matrix.

2. **[Wishart] Prove: If $S \sim W_p(n, \Sigma)$, then its trace $\operatorname{tr}(S)$ is a sum of squared independent normal variables when $\Sigma=I$.**
   ??? success "Solution"
       $\operatorname{tr}(S) = \operatorname{tr}(\sum X_i X_i^T) = \sum X_i^T X_i = \sum_{i,j} X_{ij}^2$. If $\Sigma=I$, the $X_{ij}$ are i.i.d. standard normals, so $\operatorname{tr}(S)$ follows a Chi-squared distribution with $np$ degrees of freedom.

3. **[Jacobian] Calculate the Jacobian of the matrix transformation $Y = AXB$.**
   ??? success "Solution"
       Using the differential $dY = A(dX)B \implies \operatorname{vec}(dY) = (B^T \otimes A) \operatorname{vec}(dX)$. The Jacobian is the determinant of the representation matrix: $|B^T \otimes A| = (\det B)^n (\det A)^p$.

4. **[Determinant Expectation] Find the expected value of $\det(S)$ for $S \sim W_p(n, I)$.**
   ??? success "Solution"
       The determinant of a Wishart matrix relates to the product of independent Chi-squared variables. $\mathbb{E}[\det S] = \prod_{i=0}^{p-1} (n-i)$. This reflects the volume evolution of the sample cluster in $p$-dimensional space.

5. **[Inversion] Define the Inverse Wishart distribution and its role in Bayesian statistics.**
   ??? success "Solution"
       If $S \sim W_p(n, \Sigma)$, then $S^{-1}$ follows an Inverse Wishart distribution. It is the conjugate prior for the covariance matrix of a multivariate normal distribution, enabling efficient posterior updates in Bayesian inference.

6. **[Beta Distribution] Describe the Matrix Beta distribution in terms of two independent Wishart matrices.**
   ??? success "Solution"
       Let $S_1 \sim W_p(n_1, \Sigma)$ and $S_2 \sim W_p(n_2, \Sigma)$. The matrix $U = (S_1+S_2)^{-1/2} S_1 (S_1+S_2)^{-1/2}$ follows a Matrix Beta distribution. It generalizes the ratio of Chi-squared variables to the matrix domain.

7. **[Bartlett Decomposition] Explain the Bartlett decomposition of a Wishart matrix.**
   ??? success "Solution"
       A Wishart matrix $S \sim W_p(n, I)$ can be factored as $S = T T^T$, where $T$ is a lower triangular matrix with $T_{ii}^2 \sim \chi_{n-i+1}^2$ and $T_{ij} \sim \mathcal{N}(0, 1)$ for $i > j$. This provides an efficient way to simulate Wishart samples.

8. **[Characteristic Function] Write the form of the characteristic function for the Matrix Normal distribution.**
   ??? success "Solution"
       $\phi_X(Z) = \exp(i \operatorname{tr}(Z^T M) - \frac{1}{2} \operatorname{tr}(Z^T U Z V))$. The quadratic form in the trace captures the aggregated variance structure across all matrix entries.

9. **[Singular Wishart] When does a Wishart matrix become singular?**
   ??? success "Solution"
       A Wishart matrix $W_p(n, \Sigma)$ is singular with probability 1 if the number of samples $n$ is less than the dimension $p$. In this case, the sample covariance matrix is rank-deficient and does not possess a density with respect to the Lebesgue measure on the PSD cone.

10. **[Entropy] How does the entropy of a Matrix Normal distribution relate to the determinants of $U$ and $V$?**
    ??? success "Solution"
        The entropy is proportional to $\log \det(V \otimes U) = n \log \det V + p \log \det U$. This shows that the information content (uncertainty) of the matrix is the sum of the uncertainties in its row and column structures.

## Chapter Summary

This chapter explores the statistical distribution theory of matrix variables:

1. **Structured Uncertainty**: Defined Matrix Normal distributions using Kronecker products to capture row and column correlations.
2. **Covariance Dynamics**: Established the Wishart distribution as the fundamental model for sample covariance matrices.
3. **Algebraic Geometry of Randomness**: Utilized matrix Jacobians to derive the densities of transformed random matrices.
4. **Bayesian Conjugacy**: Linked Inverse Wishart and Matrix Beta distributions to the estimation of high-dimensional uncertainty.
