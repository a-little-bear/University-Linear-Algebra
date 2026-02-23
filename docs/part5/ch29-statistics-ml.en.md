# Chapter 29: Linear Algebra in Statistics and Machine Learning

<div class="context-flow" markdown>

**Prerequisites**: Singular Value Decomposition (Ch11) · Matrix Calculus (Ch47A) · Projections and Least Squares (Ch07)

**Chapter Outline**: Data as Matrices → Data Preprocessing: Centering and Scaling → Matrix Expression of Statistics: The Covariance Matrix $\Sigma$ → Core Dimensionality Reduction: PCA derived via SVD → Regression Analysis: The Projection View of OLS → Linear Classification: Hyperplane Geometry in SVM and Logistic Regression → Regularization: Algebraic Meaning of the Ridge Parameter $\lambda$ → Applications: Feature Extraction, Denoising, and Latent Semantic Analysis (LSA)

**Extension**: Linear algebra is the "underlying architecture" of machine learning; it transforms tedious data points into geometric distributions in high-dimensional space. It proves that the process of learning is essentially finding optimal projection subspaces and separating hyperplanes—the mathematical common term for understanding all deep learning algorithms.

</div>

In machine learning, every sample is a vector, and every feature is a dimension. **Linear Algebra** weaves these scattered data points into matrices, allowing us to use geometric intuition to find patterns. Whether compressing information via **PCA** or predicting trends via **Least Squares**, the underlying logic is always finding an optimal projection in vector space. This chapter introduces the algebraic framework serving as the engine for AI.

---

## 29.1 Data Matrices and Covariance

!!! definition "Definition 29.1 (Centered Data Matrix $X$)"
    Given $n$ samples and $d$ features, the matrix $X \in \mathbb{R}^{n \times d}$ obtained by subtracting the mean of each column is the centered matrix.
    The **Covariance Matrix** is: $\Sigma = \frac{1}{n-1} X^T X$.
    It describes the linear correlation between features.

---

## 29.2 Principal Component Analysis (PCA)

!!! technique "Technique: SVD-based PCA"
    To reduce dimensionality while preserving maximum variance:
    1.  Perform SVD on the centered matrix $X$: $X = U \Sigma V^T$.
    2.  The columns of the right singular matrix $V$ are the principal component directions.
    3.  The eigenvalues (squares of singular values) represent the data variance along these directions.
    **Significance**: PCA is essentially finding the "principal axes" of the data distribution.

---

## 29.3 Linear Regression and Regularization

!!! theorem "Theorem 29.1 (Normal Equation)"
    The least squares solution $\mathbf{w} = (X^T X)^{-1} X^T \mathbf{y}$ is the orthogonal projection of the target vector $\mathbf{y}$ onto the column space of $X$. When $X^T X$ is ill-conditioned, we introduce Ridge Regression: $\mathbf{w} = (X^T X + \lambda I)^{-1} X^T \mathbf{y}$.

---

## Exercises

**1. [Basics] If feature 1 has variance 4, feature 2 has variance 9, and their covariance is 3, write the $2 \times 2$ covariance matrix.**

??? success "Solution"
    **Construction:**
    1. Diagonals are variances: $a_{11}=4, a_{22}=9$.
    2. Off-diagonals are covariances: $a_{12}=a_{21}=3$.
    **Conclusion**: $\Sigma = \begin{pmatrix} 4 & 3 \\ 3 & 9 \end{pmatrix}$.

**2. [PCA] In PCA, why do we keep the directions with the largest eigenvalues?**

??? success "Solution"
    **Algebraic Explanation:**
    1. The eigenvalue represents the **variance** of the data along that direction.
    2. Greater variance implies higher discriminability and information content.
    3. Discarding small eigenvalue directions is equivalent to removing low-energy noise or unimportant details, achieving optimal information compression (see Ch11 Best Low-Rank Approximation).

**3. [Calculation] Given data points $(1, 1)^T$ and $(-1, -1)^T$, find the covariance matrix.**

??? success "Solution"
    **Steps:**
    1. Mean is $(0, 0)^T$. Data is already centered.
    2. Matrix $X = \begin{pmatrix} 1 & 1 \\ -1 & -1 \end{pmatrix}$.
    3. $X^T X = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ -1 & -1 \end{pmatrix} = \begin{pmatrix} 2 & 2 \\ 2 & 2 \end{pmatrix}$.
    4. $\Sigma = \frac{1}{2-1} \begin{pmatrix} 2 & 2 \\ 2 & 2 \end{pmatrix} = \begin{pmatrix} 2 & 2 \\ 2 & 2 \end{pmatrix}$.

**4. [Geometry] What does the decision boundary $w^T x + b = 0$ represent geometrically?**

??? success "Solution"
    **Conclusion: A Hyperplane.**
    The vector $w$ is the **normal vector** of the hyperplane, determining its orientation. The bias $b$ determines the displacement from the origin. Classification is essentially determining which side of the hyperplane a point $x$ lies on.

**5. [Logistic Regression] Briefly explain the origin of matrix multiplication in the gradient update for logistic regression.**

??? success "Solution"
    The gradient of the loss function with respect to $w$ can be written as $\nabla J(w) = X^T (\mathbf{p} - \mathbf{y})$. Here $X^T$ is the transpose of the feature matrix, and $(\mathbf{p} - \mathbf{y})$ is the error vector. This matrix-vector product computes the update for all parameters simultaneously, serving as the basis for parallelized training.

**6. [Regularization] Why does adding $\lambda I$ in Ridge regression solve the non-invertibility of $X^T X$?**

??? success "Solution"
    **Algebraic Reason:**
    1. $X^T X$ is positive semi-definite, so eigenvalues are $\ge 0$. If columns are dependent, the smallest eigenvalue is 0.
    2. The eigenvalues of $X^T X + \lambda I$ are $\lambda_i + \lambda$.
    3. As long as $\lambda > 0$, all eigenvalues are strictly positive.
    **Conclusion**: The matrix becomes strictly positive definite, guaranteeing the existence of the inverse and numerical stability.

**7. [Application] What is a "Fully Connected Layer" in a neural network?**

??? success "Solution"
    A fully connected layer is essentially a **matrix multiplication** $y = Wx + b$. Each row of matrix $W$ represents the weight vector of one neuron. Linear algebra proves that every layer performs a linear transformation (rotation, scaling, shear) of the space.

**8. [Calculation] If the singular values of data matrix $X$ are $\{10, 5, 0.1\}$, find the proportion of variance explained by the first two principal components.**

??? success "Solution"
    **Steps:**
    1. Total variance is proportional to the sum of squared singular values: $10^2 + 5^2 + 0.1^2 = 125.01$.
    2. Variance contribution of the first two: $100 + 25 = 125$.
    3. Ratio: $125 / 125.01 \approx 99.99\%$.
    **Conclusion**: The first two principal components capture nearly all the information.

**9. [Convexity] Prove the loss function of logistic regression is convex.**

??? success "Solution"
    **Proof Key:**
    The Hessian of the loss function takes the form $X^T D X$, where $D$ is a diagonal matrix with entries $p_i(1-p_i) > 0$. Since this is a Gram-matrix form, the Hessian is positive semi-definite. According to convexity criteria (Ch64B), the optimization problem has a unique global minimum.

**10. [Application] Describe the role of SVD in recommendation systems (LSA).**

??? success "Solution"
    By performing truncated SVD on the "User-Item" interaction matrix, $X \approx U_k \Sigma_k V_k^T$. The columns of $U_k$ extract latent user preferences, while columns of $V_k$ extract latent item features. This low-rank decomposition automatically merges similar features, solving the data sparsity problem.

## Chapter Summary

Linear algebra is the "geometric skeleton" of modern intelligent algorithms:

1.  **Compression of Space**: Through PCA and SVD, linear algebra proves that high-dimensional data often collapses onto low-dimensional linear manifolds, establishing the mathematical standard for feature extraction.
2.  **Projection of Errors**: From OLS to Logistic Regression, the learning process is unified as finding an algebraic solution that is closest to the target vector (projection) or separates classes most cleanly (hyperplane).
3.  **Acceleration of Computation**: Matrix-based symbolic expressions not only simplify derivations but also leverage high-performance libraries like BLAS/LAPACK to turn training into lightning-fast parallel matrix operations, supporting today's massive AI models.
