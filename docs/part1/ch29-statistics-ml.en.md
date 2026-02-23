# Chapter 29: Applications of Linear Algebra in Statistics and Machine Learning

<div class="context-flow" markdown>

**Prerequisites**: SVD (Ch11) · Matrix Analysis (Ch14) · Positive Definite Matrices (Ch16) · Optimization (Ch25)

**Chapter Outline**: Data Matrix and Centering → Covariance Matrix and Eigendecomposition → Principal Component Analysis (PCA) → Linear Discriminant Analysis (LDA) → Least Squares and Projection Operators → Ridge Regression and Tikhonov Regularization → Weight Matrices and Backpropagation in Neural Networks → Kernel Methods and PD Kernels → Introduction to Dimensionality Reduction

**Extension**: Machine learning is geometry in high-dimensional spaces; almost all mainstream algorithms can be viewed as performing matrix decompositions or projections under specific constraints.

</div>

Statistics and machine learning study how to extract structure from data. From the perspective of linear algebra, feature extraction is subspace projection, dimensionality reduction is low-rank approximation, and model training is optimization of large-scale matrix parameters.

---

## 29.1 Data Representation and Core Algorithms

!!! definition "Definition 29.1 (Sample Covariance Matrix)"
    Let $X$ be an $n 	imes d$ centered data matrix. Its sample covariance matrix is $S = \frac{1}{n-1} X^T X$, which is symmetric and positive semi-definite.

!!! theorem "Theorem 29.3 (PCA and Eigendecomposition)"
    The first $k$ principal directions of PCA correspond to the eigenvectors associated with the $k$ largest eigenvalues of the covariance matrix $S$.

---

## Exercises

1. **[Centering] Given data matrix $X = \begin{pmatrix} 1 & 2 \ 3 & 4 \end{pmatrix}$, compute its centered version.**
   ??? success "Solution"
       Mean vector $\mu = [(1+3)/2, (2+4)/2] = [2, 3]$.
       Centered matrix $X_c = X - \mathbf{1}\mu = \begin{pmatrix} 1-2 & 2-3 \ 3-2 & 4-3 \end{pmatrix} = \begin{pmatrix} -1 & -1 \ 1 & 1 \end{pmatrix}$.

2. **[PCA and Variance] Why does PCA choose the directions with the largest eigenvalues?**
   ??? success "Solution"
       Eigenvalues represent the variance of the data along their corresponding eigenvectors. To retain as much info as possible (diversity) during dimensionality reduction, we project onto directions with maximum variance.

3. **[Least Squares] Prove: The least squares solution $\hat{\beta} = (X^T X)^{-1} X^T y$ is the coefficient of projection of $y$ onto the column space of $X$.**
   ??? success "Solution"
       Projection vector $p = X\hat{\beta} = X(X^T X)^{-1} X^T y = P_X y$. Here $P_X = X(X^T X)^{-1} X^T$ is exactly the orthogonal projection matrix onto the column space. Minimizing error is equivalent to finding the projection.

4. **[Ridge Regression] In ridge regression, how does the regularization term $\lambda I$ affect the eigenvalues of $X^T X$?**
   ??? success "Solution"
       It adds $\lambda$ to each eigenvalue $\lambda_i$. This improves the condition number, solving instability issues when the data matrix is nearly rank-deficient, thereby enhancing model robustness.

5. **[LDA vs PCA] Briefly describe the difference between PCA and LDA in subspace selection.**
   ??? success "Solution"
       - **PCA**: Unsupervised, seeks directions that maximize **total variance**.
       - **LDA**: Supervised, seeks directions that maximize **between-class variance** while minimizing **within-class variance**.

6. **[Weight Matrix] In a fully connected layer $y = \sigma(Wx + b)$, what role does the weight matrix $W$ play?**
   ??? success "Solution"
       $W$ performs a linear mapping (rotation, scaling, and dim reduction/increase) from input space to output space. It is the core feature transformation operator learned by the system.

7. **[SVD and Collaborative Filtering] How is SVD used in recommendation systems for user-rating matrices?**
   ??? success "Solution"
       By truncated SVD, decompose the sparse matrix $R$ into $U \Sigma V^T$, where $U$ represents latent user preferences and $V$ represents latent item features. Low-rank approximation fills missing entries to predict ratings for unobserved items.

8. **[Kernel Trick] Why is Mercer's Theorem important for machine learning?**
   ??? success "Solution"
       Mercer's Theorem guarantees that if a kernel function satisfies positive definiteness, it corresponds to an inner product in some high-dimensional space. This allows computing high-dim projections (like SVM) in low-dim space without constructing high-dim vectors.

9. **[Calculation] Find the principal directions for $X = \begin{pmatrix} 1 & 1 \ 1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       $S \propto X^T X = \begin{pmatrix} 2 & 2 \ 2 & 2 \end{pmatrix}$. Eigenvalues are 4 and 0. The eigenvector for 4 is $[1/\sqrt{2}, 1/\sqrt{2}]^T$. This is the first principal component.

10. **[Singular Values and Stability] How does the distribution of singular values of $X$ reflect its quality?**
    ??? success "Solution"
        If singular values decay rapidly, the data has a strong low-rank structure, making it suitable for dim reduction. If all are large and close, the dimensions are nearly independent, and dim reduction will lose significant information.

## Chapter Summary

Machine learning is the modern frontier of linear algebra:

1. **Projection is Core**: All estimates and approximations are essentially finding the closest point in a subspace.
2. **Spectral Insight**: Eigenvalues reveal the primary contradictions in data, determining the line between signal and noise.
3. **Regularization as Constraint**: By manually adjusting eigenvalues, we strike an algebraic balance between precision and stability.
