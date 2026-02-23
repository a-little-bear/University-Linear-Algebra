# Chapter 72B: Matrix Methods in Multivariate Statistical Inference

<div class="context-flow" markdown>

**Prerequisites**: Singular Value Decomposition (Ch11) · Eigenvalues (Ch6) · Covariance (Ch72A) · Projections (Ch5)

**Chapter Outline**: Matrix Derivative Forms of Maximum Likelihood Estimation → Geometric Projection Perspective of Linear Regression → Generalized Least Squares (GLS) → Discriminant Analysis (LDA) → Canonical Correlation Analysis (CCA) → Factor Analysis → Foundations of Structural Equation Modeling (SEM)

**Extension**: Multivariate statistical inference is the theoretical root of modern machine learning (from parametric regression to deep learning weight optimization).

</div>

Statistical inference aims to estimate true operator parameters from finite samples. Linear algebra solves the geometric uniqueness problem of least squares estimation through orthogonal projection theory and establishes optimal criteria for multidimensional data classification and correlation analysis using generalized eigenvalue theory.

---

## 72B.1 Orthogonal Projections and Estimation Theory

!!! definition "Definition 72B.1 (Least Squares Operator)"
    Given an observation vector $\mathbf{y} \in \mathbb{R}^n$ and a design matrix $X \in \mathbb{R}^{n 	imes p}$. The least squares estimate of linear regression coefficients is $\hat{\beta} = (X^T X)^{-1} X^T \mathbf{y}$. Geometrically, $\hat{\mathbf{y}} = X\hat{\beta}$ is the orthogonal projection of $\mathbf{y}$ onto the column space of $X$.

!!! theorem "Theorem 72B.1 (Generalized Eigenvalue Problem in Discriminant Analysis)"
    In Fisher Linear Discriminant Analysis (LDA), the optimal projection direction $\mathbf{w}$ is determined by the following generalized eigenvalue problem:
    $$S_b \mathbf{w} = \lambda S_w \mathbf{w}$$
    where $S_b$ is the between-class scatter matrix and $S_w$ is the within-class scatter matrix. The eigenvector corresponding to the largest eigenvalue maximizes the ratio of between-class separation to within-class dispersion (Rayleigh quotient).

---

## Exercises

1. **[Normal Equations] From the geometric condition that the residual vector $\mathbf{e} = \mathbf{y} - X\beta$ is orthogonal to the column space of $X$, derive the normal equations $X^T X \beta = X^T \mathbf{y}$ for the least squares method.**
   ??? success "Solution"
       Orthogonality requires $X^T \mathbf{e} = 0$. Substituting $\mathbf{e}$ yields $X^T (\mathbf{y} - X\beta) = 0 \implies X^T \mathbf{y} - X^T X \beta = 0$. Rearranging gives the normal equations.

2. **[Projection Matrix] Prove that the projection matrix $H = X(X^T X)^{-1} X^T$ is symmetric and idempotent ($H^2 = H$), and state why its eigenvalues can only be 0 or 1.**
   ??? success "Solution"
       $H^T = (X(X^T X)^{-1} X^T)^T = X((X^T X)^{-1})^T X^T = H$. $H^2 = X(X^T X)^{-1} X^T X(X^T X)^{-1} X^T = X I (X^T X)^{-1} X^T = H$. The eigenvalues $\lambda$ of an idempotent operator satisfy $\lambda^2 = \lambda$, hence $\lambda \in \{0, 1\}$.

3. **[Calculation] Given a design matrix $X = \begin{pmatrix} 1 & 1 \ 1 & 2 \end{pmatrix}$ and observation $\mathbf{y} = \begin{pmatrix} 2 \ 3 \end{pmatrix}$. Solve for $\hat{\beta}$ and verify if the predicted value $\hat{\mathbf{y}}$ lies in the column space of $X$.**
   ??? success "Solution"
       $X^T X = \begin{pmatrix} 2 & 3 \ 3 & 5 \end{pmatrix}$, its inverse is $\begin{pmatrix} 5 & -3 \ -3 & 2 \end{pmatrix}$. $X^T \mathbf{y} = \begin{pmatrix} 5 \ 8 \end{pmatrix}$.
       $\hat{\beta} = \begin{pmatrix} 5 & -3 \ -3 & 2 \end{pmatrix} \begin{pmatrix} 5 \ 8 \end{pmatrix} = \begin{pmatrix} 1 \ 1 \end{pmatrix}$. $\hat{\mathbf{y}} = X \hat{\beta} = [2, 3]^T = \mathbf{y}$ (here, since $X$ is full rank and $n=p$, the projection is the identity operator).

4. **[CCA Analysis] Describe how Canonical Correlation Analysis (CCA) transforms the problem of maximizing the correlation between two sets of random variables into an SVD problem of the product of two covariance matrix inverse roots.**
   ??? success "Solution"
       The goal is to maximize $\operatorname{corr}(\mathbf{a}^T \mathbf{X}, \mathbf{b}^T \mathbf{Y})$. This is algebraically equivalent to analyzing the singular values of the operator $\Sigma_{XX}^{-1/2} \Sigma_{XY} \Sigma_{YY}^{-1/2}$. The singular vector pairs are the linear combination coefficients that achieve maximum correlation.

5. **[Factor Analysis] Analyze the algebraic division of labor between the low-rank matrix $\Lambda \Lambda^T$ and the residual diagonal matrix $\Psi$ in capturing information in the factor analysis model $\Sigma = \Lambda \Lambda^T + \Psi$.**
   ??? success "Solution"
       $\Lambda \Lambda^T$ captures the shared variance among variables (common factors) via a rank-$k$ ($k \ll p$) approximation, corresponding to the correlation structure; $\Psi$ characterizes the unique variance (noise) independent to each variable. This is a stochastic form of structured low-rank decomposition.

6. **[Gauss-Markov] Explain the relationship between the "minimum variance" property of the least-squares estimator in the Gauss-Markov theorem and the eigenvalues of the $(X^T X)$ inverse matrix.**
   ??? success "Solution"
       The covariance of the estimator $\hat{\beta}$ is $\sigma^2 (X^T X)^{-1}$. Minimum variance means that in specific directions, the estimation uncertainty is bounded by the spectral distribution of $(X^T X)$. Orthogonality of the design matrix ($X^T X = I$) achieves optimal balance in variance distribution.

7. **[Ridge Regression] Prove that the ridge regression estimator $\hat{\beta}_{ridge} = (X^T X + \lambda I)^{-1} X^T \mathbf{y}$ is a singular value shrinkage of the original least squares solution.**
   ??? success "Solution"
       Using SVD expansion $X = U \Sigma V^T$. The original solution contains $1/\sigma_i$, which ridge regression replaces with $\sigma_i / (\sigma_i^2 + \lambda)$. As $\sigma_i 	o 0$, this term approaches 0 instead of infinity, thus achieving regularization suppression of numerical instability.

8. **[Generalized Least Squares] Derive the algebraic form of the Generalized Least Squares (GLS) solution $\hat{\beta} = (X^T V^{-1} X)^{-1} X^T V^{-1} \mathbf{y}$ when the error covariance matrix is $V$.**
   ??? success "Solution"
       Pre-process the original observations (Whitening) via a linear transformation $L^{-1}$ (where $LL^T = V$), turning errors into white noise. Apply standard least squares in the new space and transform back to the original space to obtain the GLS expression.

9. **[LDA Eigenvectors] Explain why in LDA, we typically only focus on the first $C-1$ characteristic directions (where $C$ is the number of classes).**
   ??? success "Solution"
       The between-class scatter matrix $S_b$ is defined as the weighted sum of outer products of $C$ mean vectors. Since these $C$ vectors are constrained by the total mean, their linearly independent rank is at most $C-1$. Thus, $S_b$ has only $C-1$ non-zero eigenvalues.

10. **[Statistical Coherence] Analyze the explosion behavior of $(X^T X)^{-1}$ elements when strong linear correlation (multicollinearity) exists among the column vectors of the design matrix $X$.**
    ??? success "Solution"
        Multicollinearity causes $X^T X$ to approach singularity, with its smallest eigenvalue $\lambda_{min} 	o 0$. The inverse matrix terms contain $1/\lambda_{min}$, leading to maximized estimation variance. Algebraically, this manifests as the estimator parameters exhibiting extreme sensitivity to observational perturbations.

## Chapter Summary

This chapter discusses the algebraic implementations of core algorithms in multivariate statistical inference:

1. **Geometric Estimation**: Established the optimal criteria and computational paradigms for linear regression through orthogonal projections.
2. **Dimension Compression**: Demonstrated the application of SVD and eigenvalue decomposition in extracting core data structures (Factor Analysis, CCA).
3. **Classification and Discrimination**: Utilized generalized eigenvalue theory to establish optimal hyperplane criteria for linear space partitioning.
4. **Regularization Control**: Explored stability compensation mechanisms for matrix inversion, establishing robust links from theoretical estimation to numerical computation.
