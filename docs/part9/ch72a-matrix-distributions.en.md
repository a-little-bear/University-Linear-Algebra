# Chapter 72A: Matrix Probability Distributions

<div class="context-flow" markdown>

**Prerequisites**: Inner Product (Ch8) · Eigenvalues (Ch6) · Matrix Exponential (Ch13) · Random Matrices (Ch23)

**Chapter Outline**: Random Vectors and Covariance Matrices → Multivariate Normal Distribution (Application of Positive Definite Matrices) → Wishart Distribution (Foundation of Random Matrices) → Matrix Normal Distribution → Role of Kronecker Product in Probability → Matrix Form of Characteristic Functions

**Extension**: Matrix distribution theory is the mathematical floor for financial portfolio analysis, wireless communication channel modeling, and large-scale Gaussian process regression.

</div>

Matrix probability theory studies the second-moment structures between multidimensional random variables and their distribution characteristics in matrix spaces. By expressing probability density functions as functionals of matrix invariants (trace, determinant, quadratic forms), linear algebra provides rigorous parametric means for describing high-dimensional stochastic associations.

---

## 72A.1 Covariance Structure and Multivariate Normal Distribution

!!! definition "Definition 72A.1 (Covariance Operator)"
    Let $\mathbf{X} \in \mathbb{R}^p$ be a random vector. Its covariance matrix $\Sigma = \mathbb{E}[(\mathbf{X}-\mu)(\mathbf{X}-\mu)^T]$ is a symmetric semi-positive definite matrix. $\Sigma$ encodes the linear correlations between dimensions and the distribution of variance.

!!! theorem "Theorem 72A.1 (Wishart Random Matrix)"
    Let $\mathbf{X}_1, \dots, \mathbf{X}_n$ be independent random vectors following $N_p(0, \Sigma)$. Then the scatter matrix $S = \sum \mathbf{X}_i \mathbf{X}_i^T$ follows a Wishart distribution with parameters $(n, \Sigma)$. It is the algebraic foundation for sample covariance matrix analysis.

---

## Exercises

1. **[Positive Semi-definiteness] Prove: For any random vector, its covariance matrix $\Sigma$ must be positive semi-definite.**
   ??? success "Solution"
       For any constant vector $\mathbf{a} \in \mathbb{R}^p$, the variance of the scalar random variable $\mathbf{a}^T \mathbf{X}$ is $\operatorname{Var}(\mathbf{a}^T \mathbf{X}) = \mathbf{a}^T \Sigma \mathbf{a}$. Since the variance of a real random variable must be non-negative, this quadratic form is non-negative for all $\mathbf{a}$, thus $\Sigma \succeq 0$.

2. **[Linear Transformation] Let $\mathbf{Y} = A \mathbf{X} + \mathbf{b}$, where $\mathbf{X} \sim (\mu, \Sigma)$. Derive the expression for the covariance matrix of $\mathbf{Y}$.**
   ??? success "Solution"
       $\Sigma_Y = \mathbb{E}[A(\mathbf{X}-\mu)(A(\mathbf{X}-\mu))^T] = A \mathbb{E}[(\mathbf{X}-\mu)(\mathbf{X}-\mu)^T] A^T = A \Sigma A^T$. This demonstrates the evolution rules of covariance under affine transformations.

3. **[Independence Determination] Prove: In a multivariate normal distribution, the components are independent if and only if the covariance matrix $\Sigma$ is diagonal.**
   ??? success "Solution"
       The density function of a multivariate normal distribution contains the term $\exp(-\frac{1}{2}(\mathbf{x}-\mu)^T \Sigma^{-1} (\mathbf{x}-\mu))$. If $\Sigma$ is diagonal, the quadratic form decomposes into a sum of squares of individual components, and the density function degenerates into a product of univariate normal densities, satisfying the definition of independence.

4. **[Mahalanobis Distance] Calculate the inverse of the matrix $A = \begin{pmatrix} 1 & 0.5 \ 0.5 & 1 \end{pmatrix}$ and explain its metric role in the Mahalanobis distance $\sqrt{\mathbf{x}^T \Sigma^{-1} \mathbf{x}}$.**
   ??? success "Solution"
       $\Sigma^{-1} = \frac{4}{3} \begin{pmatrix} 1 & -0.5 \ -0.5 & 1 \end{pmatrix}$. The inverse matrix acts as a weighting adjuster: in directions with strong correlation (corresponding to large eigenvalues of $\Sigma$), the inverse matrix shrinks the Euclidean distance in that direction, achieving normalization of data scale and correlation.

5. **[Wishart Expectation] Prove that the expectation of the sample scatter matrix $S$ is $\mathbb{E}[S] = n \Sigma$.**
   ??? success "Solution"
       $\mathbb{E}[\sum \mathbf{X}_i \mathbf{X}_i^T] = \sum \mathbb{E}[\mathbf{X}_i \mathbf{X}_i^T]$. Since $\mathbb{E}[\mathbf{X}_i]=0$, we have $\mathbb{E}[\mathbf{X}_i \mathbf{X}_i^T] = \Sigma$. Summing $n$ times yields $n \Sigma$.

6. **[Precision Matrix] Define the precision matrix $\Omega = \Sigma^{-1}$ and explain its components $\omega_{ij}$'s relationship with partial correlation coefficients.**
   ??? success "Solution"
       $\Omega$ reflects conditional correlations. In a Gaussian graphical model, $\omega_{ij} = 0$ implies that variables $i$ and $j$ are conditionally independent given all other variables. Its standardized components are the partial correlation coefficients.

7. **[Eigenstructure] Analyze the equivalence between the eigenvalue decomposition of a covariance matrix and the direction of maximum variance in Principal Component Analysis (PCA).**
   ??? success "Solution"
       The eigenvector corresponding to the largest eigenvalue is the solution to $\max_{\mathbf{v}^T \mathbf{v}=1} \mathbf{v}^T \Sigma \mathbf{v}$. Geometrically, this corresponds to the principal axis direction where the random point cloud is most dispersed.

8. **[Information Entropy] Prove: The differential entropy $H$ of a multivariate normal distribution is linearly related to the logarithm of $\det(\Sigma)$.**
   ??? success "Solution"
       $H(\mathbf{X}) = \frac{p}{2}(1 + \ln(2\pi)) + \frac{1}{2} \ln \det(\Sigma)$. The determinant quantifies the volume occupied by the random variable in space, thereby reflecting the information uncertainty of the system.

9. **[Matrix Normality] Explain the parameter compression significance of the Kronecker product $V \otimes U$ for the covariance structure in the matrix normal distribution $\mathcal{MN}_{n 	imes p}(M, U, V)$.**
   ??? success "Solution"
       It assumes that the row-wise correlation $U$ and column-wise correlation $V$ are decoupled. This simplifies a full covariance matrix with $n^2 p^2$ elements to only $(n^2+p^2)$ parameters, greatly mitigating the curse of dimensionality in estimation.

10. **[Characteristic Function] Write the matrix form of the characteristic function $\phi(\mathbf{t}) = \mathbb{E}[e^{j \mathbf{t}^T \mathbf{X}}]$ for a multivariate normal distribution and analyze the quadratic form in its exponent.**
    ??? success "Solution"
        $\phi(\mathbf{t}) = \exp(j \mathbf{t}^T \mu - \frac{1}{2} \mathbf{t}^T \Sigma \mathbf{t})$. The quadratic form $\mathbf{t}^T \Sigma \mathbf{t}$ in the exponent determines the convergence properties of the probability density in the frequency domain and is the mapping of the covariance structure under the Fourier transform.

## Chapter Summary

This chapter discusses matrix parameterization modeling in stochastic analysis:

1. **Second-Moment Representation**: Established symmetric semi-positive definite matrices as a universal tool for describing spatial associations of random vectors.
2. **Generative Distributions**: Established the Wishart distribution as the matrix algebraic standard for sample statistic analysis.
3. **Structured Covariance**: Revealed conditional dependence and independence laws among high-dimensional data using Kronecker products and inverse matrices (precision matrices).
4. **Invariant Metrics**: Quantified entropy increase and energy distribution characteristics of stochastic processes through determinants and traces.
