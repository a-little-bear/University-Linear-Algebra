# Chapter 29: Linear Algebra in Statistics and Machine Learning

<div class="context-flow" markdown>

**Prerequisites**: Matrix Multiplication (Ch2) · Orthogonality (Ch7) · SVD (Ch11) · Matrix Calculus (Ch47A) · Probability

**Chapter Outline**: Data Matrix $X$ → Covariance and Correlation Matrices → Linear Regression and Least Squares → Normal Equations → Principal Component Analysis (PCA) → Dimension Reduction → Ridge and Lasso Regularization → Linear Discriminant Analysis (LDA) → Kernel Trick and Support Vector Machines

**Extension**: Machine learning is essentially the search for an optimal linear operator in high-dimensional feature spaces; matrices are the structures that store and transform data.

</div>

Linear algebra is the "assembly language" of modern data science. In statistics and machine learning, data is represented as a matrix $X$, and finding patterns involves decomposing this matrix. **Principal Component Analysis (PCA)** uses the spectral decomposition of the covariance matrix to find the directions of maximum variance, allowing for massive dimensionality reduction. **Linear Regression** uses projections to find the best-fitting hyperplane through a data cloud. This chapter explores how the core theorems of linear algebra solve the most important problems in learning from data.

---

## 29.1 Data Matrices and Projections

!!! definition "Definition 29.1 (Covariance Matrix)"
    For a zero-centered data matrix $X \in \mathbb{R}^{n \times p}$ (n samples, p features), the sample covariance matrix is $\Sigma = \frac{1}{n-1} X^T X$. It is always symmetric and positive semi-definite.

!!! theorem "Theorem 29.1 (The Least Squares Solution)"
    The solution to $\min \|y - X\beta\|^2$ is given by the **Normal Equations**:
    $$X^T X \beta = X^T y$$
    Geometrically, $\hat{y} = X\hat{\beta}$ is the orthogonal projection of $y$ onto the column space of $X$.

---

## Exercises

1. **[Fundamentals] Compute the covariance matrix for $X = \begin{pmatrix} 1 & 2 \\ 2 & 1 \\ 3 & 3 \end{pmatrix}$ (center the data first).**
   ??? success "Solution"
       Means are $\bar{x}_1 = 2, \bar{x}_2 = 2$. Centered $X = \begin{pmatrix} -1 & 0 \\ 0 & -1 \\ 1 & 1 \end{pmatrix}$. $X^T X = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$. $\Sigma = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$.

2. **[PCA] How do you find the first principal component of a data set?**
   ??? success "Solution"
       It is the eigenvector corresponding to the largest eigenvalue of the covariance matrix $\Sigma$. It is the direction that maximizes the variance of the projected data.

3. **[SVD] Relate PCA to the Singular Value Decomposition.**
   ??? success "Solution"
       If $X = U \Sigma V^T$ is the SVD of the centered data, the columns of $V$ are the principal components (eigenvectors of $X^T X$) and the columns of $U \Sigma$ are the projected coordinates (principal scores).

4. **[Regularization] What matrix is added to $X^T X$ in Ridge Regression?**
   ??? success "Solution"
       The term $\lambda I$. This "regularizes" the problem, making the matrix $X^T X + \lambda I$ non-singular even if $X$ is rank-deficient, preventing overfitting.

5. **[Normal Equations] When does $X^T X \beta = X^T y$ have a unique solution?**
   ??? success "Solution"
       When $X^T X$ is invertible, which occurs iff $X$ has full column rank.

6. **[Hat Matrix] Define the "Hat Matrix" $H = X(X^T X)^{-1} X^T$.**
   ??? success "Solution"
       It is the orthogonal projection matrix onto the column space of $X$. It maps the actual outcomes $y$ to the predicted values $\hat{y} = Hy$.

7. **[Dimension] If we keep the first $k$ principal components, what percentage of variance is retained?**
   ??? success "Solution"
       $\frac{\sum_{i=1}^k \lambda_i}{\sum_{j=1}^p \lambda_j}$. This cumulative variance plot (Scree plot) helps choose the optimal reduced dimension.

8. **[Kernel Trick] Briefly describe the "Kernel Trick" in linear terms.**
   ??? success "Solution"
       Mapping data to a much higher-dimensional space where it becomes linearly separable. The matrix $\Phi(X)^T \Phi(X)$ is replaced by a kernel matrix $K_{ij} = k(x_i, x_j)$ without explicitly computing the mapping.

9. **[Multicollinearity] Why is high correlation between features (multicollinearity) a problem?**
   ??? success "Solution"
       It makes $X^T X$ near-singular, leading to high variance in the estimated coefficients $\beta$ (high condition number). SVD or PCA can solve this by decorrelating the features.

10. **[Trace] Show that the total variance of a dataset is $\operatorname{tr}(\Sigma)$.**
    ??? success "Solution"
        Total variance is the sum of the variances of each feature: $\sum \sigma_i^2$. Since $\Sigma_{ii} = \sigma_i^2$, the sum is the trace. By the trace-eigenvalue identity, it is also the sum of all spectral variance $\sum \lambda_j$.

## Chapter Summary

This chapter establishes the algebraic engine of modern data analytics:

1. **Second-Moment Modeling**: Formulated the covariance matrix as the definitive capture of data correlation and variance.
2. **Geometric Regression**: Integrated orthogonal projections into the solution of linear predictive models.
3. **Spectral Compression**: Developed PCA as the optimal linear tool for discovering latent structures and reducing data dimensionality.
4. **Regularized Stability**: Linked matrix inversion to model robustness, providing the mathematical basis for handling ill-posed learning problems.
