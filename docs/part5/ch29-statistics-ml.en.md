# Chapter 29  Linear Algebra in Statistics and Machine Learning

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues/SVD (Ch6, 11) · positive definite matrices (Ch16) · least squares (Ch25) · matrix decompositions (Ch10)

**Chapter arc**: Covariance matrix (symmetric PSD) → PCA (eigendecomposition/SVD) → linear regression (normal equations/QR) → LDA (generalized eigenvalue) → kernel methods (Gram matrix) → SVM (quadratic programming) → neural networks (matrix chain multiplication/gradients) → recommender systems (low-rank matrix factorization)

</div>

Statistics and machine learning represent one of the most important application domains of linear algebra. From multivariate covariance analysis to deep learning backpropagation, linear algebra provides a unified mathematical language and efficient computational tools. This chapter systematically presents the central role of linear algebra in statistics and machine learning, covering multivariate statistics, principal component analysis, linear regression and regularization, linear discriminant analysis, kernel methods, support vector machines, neural networks, and matrix factorization in recommender systems.

---

## 29.1 Multivariate Statistics and the Covariance Matrix

<div class="context-flow" markdown>

**Linear algebra perspective**: $n$ samples, $p$ features → data matrix $X \in \mathbb{R}^{n \times p}$ → covariance matrix $\Sigma = \frac{1}{n-1}X_c^T X_c$ (positive semidefinite, Ch16) → eigenvalues = variance in each direction → Mahalanobis distance = $\Sigma^{-1}$-weighted inner product

**Link**: Ch8 inner product spaces · Ch16 positive definite matrices

</div>

The foundation of multivariate statistics is organizing data into matrices and characterizing the correlation structure among variables through the covariance matrix.

!!! definition "Definition 29.1 (Data matrix and centering)"
    Given $n$ samples, each with $p$ features, the **data matrix** $X \in \mathbb{R}^{n \times p}$ has its $i$-th row $\mathbf{x}_i^T$ as the feature vector of sample $i$. The sample mean is

    $$
    \bar{\mathbf{x}} = \frac{1}{n}\sum_{i=1}^n \mathbf{x}_i = \frac{1}{n}X^T \mathbf{1}_n.
    $$

    The **centered data matrix** is defined as

    $$
    X_c = X - \mathbf{1}_n \bar{\mathbf{x}}^T = \left(I_n - \frac{1}{n}\mathbf{1}_n\mathbf{1}_n^T\right)X = H X,
    $$

    where $H = I_n - \frac{1}{n}\mathbf{1}_n\mathbf{1}_n^T$ is the **centering matrix**, an orthogonal projection matrix ($H^2 = H$, $H^T = H$).

!!! definition "Definition 29.2 (Sample covariance matrix)"
    The **sample covariance matrix** is defined as

    $$
    S = \frac{1}{n-1}X_c^T X_c \in \mathbb{R}^{p \times p}.
    $$

    $S$ is symmetric positive semidefinite. When $n > p$ and the data do not lie in any hyperplane, $S$ is positive definite. The $(i,j)$-entry of $S$ is the sample covariance between features $i$ and $j$:

    $$
    S_{ij} = \frac{1}{n-1}\sum_{k=1}^n (x_{ki} - \bar{x}_i)(x_{kj} - \bar{x}_j).
    $$

!!! theorem "Theorem 29.1 (Spectral properties of the covariance matrix)"
    Let $S \in \mathbb{R}^{p \times p}$ be the sample covariance matrix with eigendecomposition

    $$
    S = Q \Lambda Q^T = \sum_{j=1}^p \lambda_j \mathbf{q}_j \mathbf{q}_j^T,
    $$

    where $\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_p \ge 0$. Then:

    1. The variance of the data in direction $\mathbf{q}_j$ is $\lambda_j$.
    2. **Total variance** $= \operatorname{tr}(S) = \sum_{j=1}^p \lambda_j$.
    3. **Generalized variance** $= \det(S) = \prod_{j=1}^p \lambda_j$.
    4. $\operatorname{rank}(S) = \operatorname{rank}(X_c) \le \min(n-1, p)$.

??? proof "Proof"
    For a direction $\mathbf{v}$ ($\|\mathbf{v}\| = 1$), the projection of the data onto this direction is $X_c \mathbf{v}$, with variance

    $$
    \frac{1}{n-1}\|X_c \mathbf{v}\|^2 = \frac{1}{n-1}\mathbf{v}^T X_c^T X_c \mathbf{v} = \mathbf{v}^T S \mathbf{v}.
    $$

    When $\mathbf{v} = \mathbf{q}_j$, $\mathbf{q}_j^T S \mathbf{q}_j = \lambda_j$. Total variance $\sum_j \operatorname{Var}(X_c \mathbf{e}_j) = \sum_j S_{jj} = \operatorname{tr}(S) = \sum_j \lambda_j$. Since $S = \frac{1}{n-1}X_c^T X_c$, $\operatorname{rank}(S) = \operatorname{rank}(X_c)$. Centering removes one degree of freedom, so $\operatorname{rank}(X_c) \le \min(n-1, p)$. $\blacksquare$

!!! definition "Definition 29.3 (Mahalanobis distance)"
    Let $S$ be a positive definite covariance matrix. The **Mahalanobis distance** between points $\mathbf{x}$ and $\mathbf{y}$ is

    $$
    d_M(\mathbf{x}, \mathbf{y}) = \sqrt{(\mathbf{x} - \mathbf{y})^T S^{-1} (\mathbf{x} - \mathbf{y})}.
    $$

    This is equivalent to first whitening the data via $\mathbf{z} = S^{-1/2}\mathbf{x}$ and then computing Euclidean distance. The Mahalanobis distance accounts for feature scales and correlations.

!!! example "Example 29.1"
    **Computing and analyzing a covariance matrix.** Let three 2D samples be $\mathbf{x}_1 = (1,2)^T$, $\mathbf{x}_2 = (3,4)^T$, $\mathbf{x}_3 = (5,6)^T$. The data matrix and mean:

    $$
    X = \begin{pmatrix}1&2\\3&4\\5&6\end{pmatrix}, \quad \bar{\mathbf{x}} = (3, 4)^T.
    $$

    Centered data:

    $$
    X_c = \begin{pmatrix}-2&-2\\0&0\\2&2\end{pmatrix}, \quad S = \frac{1}{2}X_c^T X_c = \frac{1}{2}\begin{pmatrix}8&8\\8&8\end{pmatrix} = \begin{pmatrix}4&4\\4&4\end{pmatrix}.
    $$

    Eigenvalues: $\lambda_1 = 8$, $\lambda_2 = 0$, with eigenvector $\mathbf{q}_1 = \frac{1}{\sqrt{2}}(1,1)^T$. The data lies entirely along the $(1,1)^T$ direction, consistent with $\operatorname{rank}(S) = 1$.

!!! example "Example 29.2"
    **Geometric meaning of Mahalanobis distance.** Let

    $$
    S = \begin{pmatrix}4&2\\2&3\end{pmatrix}, \quad S^{-1} = \frac{1}{8}\begin{pmatrix}3&-2\\-2&4\end{pmatrix}.
    $$

    The Mahalanobis distance between $\mathbf{a} = (0,0)^T$ and $\mathbf{b} = (2,1)^T$:

    $$
    d_M^2 = (2,1)\frac{1}{8}\begin{pmatrix}3&-2\\-2&4\end{pmatrix}\begin{pmatrix}2\\1\end{pmatrix} = \frac{1}{8}(2,1)\begin{pmatrix}4\\0\end{pmatrix} = 1.
    $$

    The Euclidean distance is $\sqrt{5} \approx 2.24$. The Mahalanobis distance is smaller because the direction $(2,1)^T$ aligns with a direction of larger data variance.

---

## 29.2 Principal Component Analysis (PCA)

<div class="context-flow" markdown>

**Core idea**: PCA = eigendecompose the covariance matrix = SVD of centered data → top $k$ principal components = $k$ orthogonal directions of maximum variance = best rank-$k$ approximation

**Link**: Ch11 SVD · Ch6 eigenvalues · Eckart-Young theorem

</div>

Principal component analysis is the most fundamental dimensionality reduction method in multivariate statistics, with its mathematical essence being eigenvalue decomposition or singular value decomposition.

!!! definition "Definition 29.4 (Principal components)"
    Let the centered data matrix be $X_c \in \mathbb{R}^{n \times p}$ and the sample covariance matrix $S = \frac{1}{n-1}X_c^T X_c$ have eigendecomposition $S = Q\Lambda Q^T$. The $j$-th **principal component** is

    $$
    \mathbf{z}_j = X_c \mathbf{q}_j \in \mathbb{R}^n, \quad j = 1, \ldots, p,
    $$

    where $\mathbf{q}_j$ is the $j$-th eigenvector of $S$ (ordered by decreasing eigenvalue). $\mathbf{q}_j$ is called the $j$-th **loading vector**. Retaining the top $k$ principal components gives the reduced representation $Z_k = X_c Q_k \in \mathbb{R}^{n \times k}$.

!!! theorem "Theorem 29.2 (Optimality of PCA)"
    The top $k$ principal component directions $\mathbf{q}_1, \ldots, \mathbf{q}_k$ solve the optimization problem:

    $$
    \max_{\substack{V \in \mathbb{R}^{p \times k} \\ V^T V = I_k}} \operatorname{tr}(V^T S V).
    $$

    Equivalently, among all rank-$k$ matrices, $\hat{X}_c = X_c Q_k Q_k^T$ minimizes the reconstruction error $\|X_c - \hat{X}_c\|_F^2$. The proportion of variance explained by the top $k$ components is

    $$
    \frac{\sum_{j=1}^k \lambda_j}{\sum_{j=1}^p \lambda_j}.
    $$

??? proof "Proof"
    Let $X_c = U\Sigma Q^T$ be the SVD, where $\sigma_j^2 = (n-1)\lambda_j$. For any $V$ with $V^T V = I_k$,

    $$
    \operatorname{tr}(V^T S V) = \frac{1}{n-1}\operatorname{tr}(V^T X_c^T X_c V) = \frac{1}{n-1}\|X_c V\|_F^2.
    $$

    By the Eckart-Young theorem (Ch11), $X_c Q_k Q_k^T$ is the best rank-$k$ approximation of $X_c$. Projecting onto the subspace spanned by the $k$ columns of $Q_k$ retains the maximum Frobenius norm, equivalent to maximizing $\operatorname{tr}(V^T S V)$ at $V = Q_k$. The reconstruction error is

    $$
    \|X_c - X_c Q_k Q_k^T\|_F^2 = \sum_{j=k+1}^p \sigma_j^2 = (n-1)\sum_{j=k+1}^p \lambda_j. \quad \blacksquare
    $$

!!! theorem "Theorem 29.3 (Computing PCA via SVD)"
    Let the centered data matrix have SVD $X_c = U\Sigma V^T$. Then:

    1. The principal component directions are the columns of $V$.
    2. The principal component scores are $Z = U\Sigma$, i.e., the $j$-th column of $Z$ is $\sigma_j \mathbf{u}_j$.
    3. The covariance eigenvalues are $\lambda_j = \sigma_j^2/(n-1)$.

    When $n \gg p$, directly eigendecomposing $S$ ($O(p^3)$) is more efficient; when $p \gg n$, first computing the eigendecomposition of $X_c X_c^T$ ($O(n^3)$) or using truncated SVD is more efficient.

??? proof "Proof"
    From $X_c = U\Sigma V^T$,

    $$
    S = \frac{1}{n-1}X_c^T X_c = \frac{1}{n-1}V\Sigma^T U^T U\Sigma V^T = V\left(\frac{\Sigma^2}{n-1}\right)V^T.
    $$

    This is exactly the eigendecomposition of $S$ with eigenvalues $\lambda_j = \sigma_j^2/(n-1)$ and eigenvectors as the columns of $V$. The principal component scores are $Z = X_c V = U\Sigma V^T V = U\Sigma$. $\blacksquare$

!!! example "Example 29.3"
    **PCA for dimensionality reduction.** Let the centered data matrix for 4 three-dimensional data points be

    $$
    X_c = \begin{pmatrix}2&1&0\\-1&0&1\\0&-1&-1\\-1&0&0\end{pmatrix}.
    $$

    The covariance matrix is $S = \frac{1}{3}X_c^T X_c = \frac{1}{3}\begin{pmatrix}6&2&-1\\2&2&-1\\-1&-1&2\end{pmatrix}$. The eigenvalues (approximately) are $\lambda_1 \approx 2.54$, $\lambda_2 \approx 0.98$, $\lambda_3 \approx 0.15$. The proportion of variance explained by the first two components is approximately $(2.54 + 0.98)/3.67 \approx 96\%$, so $k = 2$ principal components effectively represent the original 3D data.

---

## 29.3 Linear Regression and Regularization

<div class="context-flow" markdown>

**Linear algebra core**: $\hat{\boldsymbol{\beta}} = (X^TX)^{-1}X^T\mathbf{y}$ (normal equations) = project $\mathbf{y}$ onto $\operatorname{col}(X)$ → ridge regression = spectral filtering ($\sigma_i^2 \to \sigma_i^2 + \lambda$) → LASSO = $\ell_1$ sparsity

**Link**: Ch25 least squares · Ch11 SVD · Ch7 orthogonal projection

</div>

Linear regression is the most fundamental statistical modeling tool, with its solution and analysis entirely relying on linear algebra.

!!! definition "Definition 29.5 (Linear regression model)"
    The **linear regression model** assumes a relationship between the response variable $\mathbf{y} \in \mathbb{R}^n$ and the design matrix $X \in \mathbb{R}^{n \times p}$:

    $$
    \mathbf{y} = X\boldsymbol{\beta} + \boldsymbol{\varepsilon},
    $$

    where $\boldsymbol{\beta} \in \mathbb{R}^p$ is the regression coefficient vector and $\boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0}, \sigma^2 I_n)$ is the noise. The **ordinary least squares** (OLS) estimate is

    $$
    \hat{\boldsymbol{\beta}}_{\text{OLS}} = \arg\min_{\boldsymbol{\beta}} \|\mathbf{y} - X\boldsymbol{\beta}\|_2^2 = (X^T X)^{-1}X^T \mathbf{y}.
    $$

!!! theorem "Theorem 29.4 (Geometric and spectral interpretation of OLS)"
    Let $X = U\Sigma V^T$ be the SVD with $\operatorname{rank}(X) = r$. Then:

    1. **Geometric interpretation**: $\hat{\mathbf{y}} = X\hat{\boldsymbol{\beta}}_{\text{OLS}} = H\mathbf{y}$, where $H = X(X^TX)^{-1}X^T = U_r U_r^T$ is the orthogonal projection matrix onto $\operatorname{col}(X)$ (the hat matrix).
    2. **Spectral interpretation**: $\hat{\boldsymbol{\beta}}_{\text{OLS}} = \sum_{j=1}^r \frac{\mathbf{u}_j^T \mathbf{y}}{\sigma_j}\mathbf{v}_j$.
    3. **Variance**: $\operatorname{Var}(\hat{\boldsymbol{\beta}}_{\text{OLS}}) = \sigma^2(X^TX)^{-1} = \sigma^2 V \Sigma^{-2} V^T$, so small singular values $\sigma_j$ lead to large estimation variance.

??? proof "Proof"
    The normal equations give $X^TX\hat{\boldsymbol{\beta}} = X^T\mathbf{y}$. Substituting the SVD: $(V\Sigma^2 V^T)\hat{\boldsymbol{\beta}} = V\Sigma U^T\mathbf{y}$, so $\hat{\boldsymbol{\beta}} = V\Sigma^{-1}U^T\mathbf{y} = \sum_j \frac{\mathbf{u}_j^T\mathbf{y}}{\sigma_j}\mathbf{v}_j$. The hat matrix is $H = X(X^TX)^{-1}X^T = U\Sigma V^T \cdot V\Sigma^{-2}V^T \cdot V\Sigma U^T = UU^T$ (taking the first $r$ columns). The variance follows from $\hat{\boldsymbol{\beta}} = (X^TX)^{-1}X^T\mathbf{y}$ and $\operatorname{Var}(\mathbf{y}) = \sigma^2 I$, giving $\operatorname{Var}(\hat{\boldsymbol{\beta}}) = \sigma^2(X^TX)^{-1}X^T X (X^TX)^{-1} = \sigma^2(X^TX)^{-1}$. $\blacksquare$

!!! definition "Definition 29.6 (Ridge regression and LASSO)"
    **Ridge regression** adds an $\ell_2$ penalty to OLS:

    $$
    \hat{\boldsymbol{\beta}}_{\text{ridge}} = \arg\min_{\boldsymbol{\beta}} \|\mathbf{y} - X\boldsymbol{\beta}\|_2^2 + \lambda\|\boldsymbol{\beta}\|_2^2 = (X^TX + \lambda I)^{-1}X^T\mathbf{y}.
    $$

    **LASSO** (Least Absolute Shrinkage and Selection Operator) adds an $\ell_1$ penalty:

    $$
    \hat{\boldsymbol{\beta}}_{\text{LASSO}} = \arg\min_{\boldsymbol{\beta}} \frac{1}{2}\|\mathbf{y} - X\boldsymbol{\beta}\|_2^2 + \lambda\|\boldsymbol{\beta}\|_1.
    $$

    LASSO tends to produce sparse solutions (some $\beta_j = 0$), thereby simultaneously performing feature selection.

!!! theorem "Theorem 29.5 (Spectral filtering interpretation of ridge regression)"
    Let $X = U\Sigma V^T$. The ridge regression solution can be written as

    $$
    \hat{\boldsymbol{\beta}}_{\text{ridge}} = \sum_{j=1}^r \frac{\sigma_j^2}{\sigma_j^2 + \lambda} \cdot \frac{\mathbf{u}_j^T \mathbf{y}}{\sigma_j} \mathbf{v}_j.
    $$

    The factor $\frac{\sigma_j^2}{\sigma_j^2 + \lambda}$ is a **spectral filter**: it barely attenuates large singular values ($\sigma_j^2 \gg \lambda$) but strongly suppresses small singular values ($\sigma_j^2 \ll \lambda$). This is the key mechanism for handling multicollinearity.

??? proof "Proof"
    $(X^TX + \lambda I)^{-1} = (V\Sigma^2 V^T + \lambda VV^T)^{-1} = V(\Sigma^2 + \lambda I)^{-1}V^T$ (assuming $V$ is square or handled accordingly). Thus

    $$
    \hat{\boldsymbol{\beta}}_{\text{ridge}} = V(\Sigma^2 + \lambda I)^{-1}\Sigma U^T\mathbf{y} = \sum_{j=1}^r \frac{\sigma_j}{\sigma_j^2+\lambda}(\mathbf{u}_j^T\mathbf{y})\mathbf{v}_j = \sum_{j=1}^r \frac{\sigma_j^2}{\sigma_j^2+\lambda}\cdot\frac{\mathbf{u}_j^T\mathbf{y}}{\sigma_j}\mathbf{v}_j. \quad \blacksquare
    $$

!!! example "Example 29.4"
    **Multicollinearity and ridge regression.** Let the design matrix be

    $$
    X = \begin{pmatrix}1&1\\1&1.01\\1&0.99\end{pmatrix}, \quad \mathbf{y} = \begin{pmatrix}2\\2.5\\1.5\end{pmatrix}.
    $$

    The two columns of $X$ are nearly collinear: $\sigma_1 \approx 2.45$, $\sigma_2 \approx 0.01$. The OLS estimate is extremely sensitive to noise (variance $\propto 1/\sigma_2^2 \approx 10^4$). With $\lambda = 0.1$, ridge regression suppresses the small singular value direction by the factor $\frac{0.01^2}{0.01^2 + 0.1} \approx 0.001$, yielding a stable estimate.

---

## 29.4 Linear Discriminant Analysis (LDA)

<div class="context-flow" markdown>

**Core idea**: Fisher's discriminant = find direction $\mathbf{w}$ maximizing between-class scatter / within-class scatter → generalized eigenvalue problem $S_B\mathbf{w} = \lambda S_W\mathbf{w}$

**Link**: Ch6 generalized eigenvalues · Ch7 Rayleigh quotient

</div>

Linear discriminant analysis finds optimal classification projection directions by solving a generalized eigenvalue problem.

!!! definition "Definition 29.7 (Within-class and between-class scatter matrices)"
    Given $C$ classes with $n_c$ samples in class $c$, class means $\boldsymbol{\mu}_c$, and overall mean $\boldsymbol{\mu}$, the **within-class scatter matrix** and **between-class scatter matrix** are

    $$
    S_W = \sum_{c=1}^C \sum_{\mathbf{x} \in \text{class } c} (\mathbf{x} - \boldsymbol{\mu}_c)(\mathbf{x} - \boldsymbol{\mu}_c)^T, \quad S_B = \sum_{c=1}^C n_c (\boldsymbol{\mu}_c - \boldsymbol{\mu})(\boldsymbol{\mu}_c - \boldsymbol{\mu})^T.
    $$

    The total scatter matrix is $S_T = S_W + S_B$. Note that $\operatorname{rank}(S_B) \le C - 1$.

!!! theorem "Theorem 29.6 (Fisher's discriminant criterion)"
    Fisher's linear discriminant seeks the projection direction $\mathbf{w}$ that maximizes the **Fisher criterion function**:

    $$
    J(\mathbf{w}) = \frac{\mathbf{w}^T S_B \mathbf{w}}{\mathbf{w}^T S_W \mathbf{w}}.
    $$

    The optimal $\mathbf{w}$ is the eigenvector corresponding to the largest eigenvalue of the generalized eigenvalue problem $S_B \mathbf{w} = \lambda S_W \mathbf{w}$. When $S_W$ is invertible, this is equivalent to the eigenvector of $S_W^{-1}S_B$ for its largest eigenvalue.

    For multi-class problems ($C > 2$), the top $C - 1$ generalized eigenvectors yield the optimal $(C-1)$-dimensional discriminant subspace.

??? proof "Proof"
    $J(\mathbf{w})$ is a generalized Rayleigh quotient. Setting its gradient to zero:

    $$
    \frac{\partial}{\partial \mathbf{w}}\frac{\mathbf{w}^T S_B\mathbf{w}}{\mathbf{w}^T S_W\mathbf{w}} = 0 \implies S_B\mathbf{w} = J(\mathbf{w}) S_W \mathbf{w}.
    $$

    This is the generalized eigenvalue problem $S_B\mathbf{w} = \lambda S_W\mathbf{w}$, and the maximum of $J(\mathbf{w})$ equals the largest generalized eigenvalue $\lambda_1$. Since $\operatorname{rank}(S_B) \le C - 1$, there are at most $C - 1$ nonzero generalized eigenvalues. $\blacksquare$

!!! theorem "Theorem 29.7 (Closed-form solution for two-class LDA)"
    For binary classification ($C = 2$), $S_B = n_1 n_2/(n_1+n_2) \cdot (\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)^T$ (rank 1). The optimal discriminant direction is

    $$
    \mathbf{w}^* \propto S_W^{-1}(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2).
    $$

    Classification rule: project a new sample $\mathbf{x}$ onto $\mathbf{w}^*$ and compare the projection to the threshold $\frac{1}{2}\mathbf{w}^{*T}(\boldsymbol{\mu}_1 + \boldsymbol{\mu}_2)$.

??? proof "Proof"
    Since $S_B\mathbf{w} \propto (\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)^T\mathbf{w}$, $S_B\mathbf{w}$ always points in the direction $(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)$. From $S_B\mathbf{w} = \lambda S_W\mathbf{w}$, we get $\mathbf{w} \propto S_W^{-1}(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)$. $\blacksquare$

!!! example "Example 29.5"
    **Two-class LDA computation.** Let two classes of 2D data have: class 1 samples $(0,0)^T$, $(2,2)^T$ (mean $\boldsymbol{\mu}_1 = (1,1)^T$); class 2 samples $(4,2)^T$, $(6,4)^T$ (mean $\boldsymbol{\mu}_2 = (5,3)^T$).

    $$
    S_W = \begin{pmatrix}1&1\\1&1\end{pmatrix} + \begin{pmatrix}1&1\\1&1\end{pmatrix} = \begin{pmatrix}2&2\\2&2\end{pmatrix}.
    $$

    $S_W$ is singular, so we add regularization $S_W + 0.01I$. $\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2 = (-4,-2)^T$. The optimal direction is $\mathbf{w}^* \propto (S_W + 0.01I)^{-1}(-4,-2)^T$. After projection, the separation between class means is maximized.

---

## 29.5 Kernel Methods and Feature Spaces

<div class="context-flow" markdown>

**Core idea**: Nonlinear map $\phi: \mathbb{R}^p \to \mathcal{H}$ (high/infinite-dimensional) linearizes data → kernel trick: $k(\mathbf{x},\mathbf{y}) = \langle\phi(\mathbf{x}),\phi(\mathbf{y})\rangle$ avoids computing $\phi$ explicitly → Gram matrix $K$ replaces $X^TX$

**Link**: Ch8 inner product spaces · Ch16 positive definite matrices · Mercer's theorem

</div>

Kernel methods extend linear methods to nonlinear problems through implicit high-dimensional mappings, with their mathematical foundation in positive definite kernels and reproducing kernel Hilbert spaces.

!!! definition "Definition 29.8 (Positive definite kernel and Gram matrix)"
    A function $k: \mathbb{R}^p \times \mathbb{R}^p \to \mathbb{R}$ is called a **positive definite kernel** if for any $n$ points $\mathbf{x}_1, \ldots, \mathbf{x}_n$, the **Gram matrix**

    $$
    K = [k(\mathbf{x}_i, \mathbf{x}_j)]_{n \times n}
    $$

    is symmetric positive semidefinite. Common kernel functions include:

    - **Linear kernel**: $k(\mathbf{x}, \mathbf{y}) = \mathbf{x}^T\mathbf{y}$.
    - **Polynomial kernel**: $k(\mathbf{x}, \mathbf{y}) = (\mathbf{x}^T\mathbf{y} + c)^d$, $c \ge 0$, $d \in \mathbb{N}$.
    - **Gaussian (RBF) kernel**: $k(\mathbf{x}, \mathbf{y}) = \exp\!\left(-\frac{\|\mathbf{x}-\mathbf{y}\|^2}{2\sigma^2}\right)$.

!!! theorem "Theorem 29.8 (Mercer's theorem, finite-dimensional version)"
    A symmetric function $k: \mathbb{R}^p \times \mathbb{R}^p \to \mathbb{R}$ is a positive definite kernel if and only if there exists a mapping $\phi: \mathbb{R}^p \to \mathcal{H}$ ($\mathcal{H}$ being some Hilbert space) such that

    $$
    k(\mathbf{x}, \mathbf{y}) = \langle \phi(\mathbf{x}), \phi(\mathbf{y}) \rangle_{\mathcal{H}}.
    $$

    Equivalently, $K$ is PSD $\Leftrightarrow$ $K$ can be factored as $K = \Phi\Phi^T$, where the $i$-th row of $\Phi$ is $\phi(\mathbf{x}_i)^T$.

??? proof "Proof"
    (Sufficiency) If $k(\mathbf{x},\mathbf{y}) = \langle\phi(\mathbf{x}),\phi(\mathbf{y})\rangle$, then for any $\mathbf{c} \in \mathbb{R}^n$, $\mathbf{c}^T K\mathbf{c} = \sum_{i,j}c_ic_j\langle\phi(\mathbf{x}_i),\phi(\mathbf{x}_j)\rangle = \|\sum_i c_i\phi(\mathbf{x}_i)\|^2 \ge 0$, so $K$ is PSD.

    (Necessity) Since $K$ is PSD, it has eigendecomposition $K = Q\Lambda Q^T$. Taking $\Phi = Q\Lambda^{1/2}$ gives $K = \Phi\Phi^T$. More generally, one can construct a reproducing kernel Hilbert space $\mathcal{H}$ and set $\phi(\mathbf{x}) = k(\cdot, \mathbf{x})$. $\blacksquare$

!!! definition "Definition 29.9 (Kernel PCA)"
    **Kernel PCA** performs PCA in the feature space $\mathcal{H}$. Let the centered Gram matrix be

    $$
    \tilde{K} = HKH, \quad H = I_n - \frac{1}{n}\mathbf{1}_n\mathbf{1}_n^T.
    $$

    Eigendecomposing $\tilde{K} = Q\Lambda Q^T$, the kernel principal components are $\boldsymbol{\alpha}_j = \frac{1}{\sqrt{\lambda_j}}\mathbf{q}_j$. The projection of a new sample $\mathbf{x}$ onto the $j$-th kernel principal component is

    $$
    z_j = \sum_{i=1}^n \alpha_{ji} \tilde{k}(\mathbf{x}_i, \mathbf{x}),
    $$

    where $\tilde{k}$ is the centered kernel function.

!!! example "Example 29.6"
    **Feature map of a polynomial kernel.** For the quadratic kernel $k(\mathbf{x},\mathbf{y}) = (\mathbf{x}^T\mathbf{y})^2$ ($\mathbf{x}, \mathbf{y} \in \mathbb{R}^2$), the explicit feature map is

    $$
    \phi(x_1, x_2) = (x_1^2,\; \sqrt{2}\,x_1 x_2,\; x_2^2)^T.
    $$

    Verification: $\phi(\mathbf{x})^T\phi(\mathbf{y}) = x_1^2 y_1^2 + 2x_1 x_2 y_1 y_2 + x_2^2 y_2^2 = (x_1 y_1 + x_2 y_2)^2 = (\mathbf{x}^T\mathbf{y})^2$. The power of the kernel trick is that even when $\phi$ maps to an extremely high-dimensional space (the Gaussian kernel corresponds to infinite dimensions), computation only requires evaluating the kernel function $k(\mathbf{x},\mathbf{y})$.

---

## 29.6 Support Vector Machines

<div class="context-flow" markdown>

**Core idea**: Maximum-margin classifier = constrained optimization (QP) → dual problem involves only inner products $\mathbf{x}_i^T\mathbf{x}_j$ → kernelization → support vectors = samples with $\alpha_i > 0$

**Link**: Ch25 optimization · kernel methods (29.5)

</div>

Support vector machines are one of the most successful applications of kernel methods, with training reducing to a convex quadratic programming problem.

!!! definition "Definition 29.10 (Hard-margin SVM)"
    Given linearly separable training data $\{(\mathbf{x}_i, y_i)\}_{i=1}^n$, $y_i \in \{+1, -1\}$, the **hard-margin SVM** solves

    $$
    \min_{\mathbf{w}, b} \frac{1}{2}\|\mathbf{w}\|^2, \quad \text{s.t.} \quad y_i(\mathbf{w}^T\mathbf{x}_i + b) \ge 1, \quad i = 1, \ldots, n.
    $$

    The geometric margin is $2/\|\mathbf{w}\|$. The **soft-margin SVM** introduces slack variables $\xi_i \ge 0$:

    $$
    \min_{\mathbf{w}, b, \boldsymbol{\xi}} \frac{1}{2}\|\mathbf{w}\|^2 + C\sum_{i=1}^n \xi_i, \quad \text{s.t.} \quad y_i(\mathbf{w}^T\mathbf{x}_i + b) \ge 1 - \xi_i.
    $$

!!! theorem "Theorem 29.9 (SVM dual problem)"
    The Lagrange dual of the SVM is

    $$
    \max_{\boldsymbol{\alpha}} \sum_{i=1}^n \alpha_i - \frac{1}{2}\sum_{i,j=1}^n \alpha_i \alpha_j y_i y_j \mathbf{x}_i^T\mathbf{x}_j, \quad \text{s.t.} \quad 0 \le \alpha_i \le C, \quad \sum_{i=1}^n \alpha_i y_i = 0.
    $$

    In matrix form:

    $$
    \max_{\boldsymbol{\alpha}} \mathbf{1}^T\boldsymbol{\alpha} - \frac{1}{2}\boldsymbol{\alpha}^T Q\boldsymbol{\alpha}, \quad Q_{ij} = y_i y_j \mathbf{x}_i^T\mathbf{x}_j,
    $$

    where $Q = \operatorname{diag}(\mathbf{y}) \cdot X X^T \cdot \operatorname{diag}(\mathbf{y})$ is positive semidefinite. Kernelized SVM simply replaces $\mathbf{x}_i^T\mathbf{x}_j$ with $k(\mathbf{x}_i, \mathbf{x}_j)$.

??? proof "Proof"
    Introducing Lagrange multipliers $\alpha_i \ge 0$, the Lagrangian is

    $$
    L = \frac{1}{2}\|\mathbf{w}\|^2 - \sum_i \alpha_i[y_i(\mathbf{w}^T\mathbf{x}_i + b) - 1].
    $$

    Setting the gradient with respect to $\mathbf{w}$ to zero: $\mathbf{w} = \sum_i \alpha_i y_i \mathbf{x}_i$. Setting the gradient with respect to $b$ to zero: $\sum_i \alpha_i y_i = 0$. Substituting back into $L$:

    $$
    L = \sum_i \alpha_i - \frac{1}{2}\sum_{i,j}\alpha_i\alpha_j y_i y_j \mathbf{x}_i^T\mathbf{x}_j.
    $$

    $Q$ is PSD since $Q = D_y X X^T D_y$ where $D_y = \operatorname{diag}(\mathbf{y})$, so $\mathbf{c}^T Q\mathbf{c} = \|X^T D_y\mathbf{c}\|^2 \ge 0$. $\blacksquare$

!!! example "Example 29.7"
    **Simple SVM example.** Two classes: class $+1$: $(2,1)^T$, $(3,2)^T$; class $-1$: $(0,0)^T$, $(1,1)^T$. The dual $Q$ matrix:

    $$
    Q_{ij} = y_i y_j \mathbf{x}_i^T\mathbf{x}_j, \quad X X^T = \begin{pmatrix}5&8&0&3\\8&13&0&5\\0&0&0&0\\3&5&0&2\end{pmatrix}.
    $$

    By the KKT conditions, support vectors ($\alpha_i > 0$) lie on the margin boundary. For this dataset, the optimal hyperplane separates the two classes with maximum margin. After solving the dual QP, $\mathbf{w} = \sum_i \alpha_i y_i \mathbf{x}_i$ and the decision function is $f(\mathbf{x}) = \operatorname{sign}(\mathbf{w}^T\mathbf{x} + b)$.

---

## 29.7 Linear Algebra in Neural Networks

<div class="context-flow" markdown>

**Core idea**: Forward propagation = matrix chain multiplication + elementwise nonlinearities → backpropagation = Jacobian matrix multiplication in the chain rule → condition number of weight matrices affects training stability

**Link**: Ch2 matrix multiplication · Ch6 eigenvalues · Ch15 norms

</div>

The core operation at each layer of a deep neural network is matrix multiplication. Understanding this is essential for efficient training and analysis.

!!! definition "Definition 29.11 (Matrix representation of fully-connected neural networks)"
    An $L$-layer **fully-connected neural network** can be expressed as the composite mapping

    $$
    f(\mathbf{x}) = \sigma_L(W_L \sigma_{L-1}(W_{L-1} \cdots \sigma_1(W_1 \mathbf{x} + \mathbf{b}_1) \cdots + \mathbf{b}_{L-1}) + \mathbf{b}_L),
    $$

    where $W_\ell \in \mathbb{R}^{d_\ell \times d_{\ell-1}}$ is the weight matrix of layer $\ell$, $\mathbf{b}_\ell \in \mathbb{R}^{d_\ell}$ is the bias, and $\sigma_\ell$ is the elementwise activation function. For a mini-batch input $X \in \mathbb{R}^{d_0 \times m}$ ($m$ samples), the forward computation at layer $\ell$ is

    $$
    Z_\ell = W_\ell A_{\ell-1} + \mathbf{b}_\ell \mathbf{1}^T, \quad A_\ell = \sigma_\ell(Z_\ell),
    $$

    where $A_0 = X$. The core operation is the matrix multiplication $W_\ell A_{\ell-1}$.

!!! theorem "Theorem 29.10 (Matrix form of backpropagation)"
    Let the gradient of the loss $\mathcal{L}$ with respect to the pre-activation of layer $\ell$ be $\delta_\ell = \frac{\partial \mathcal{L}}{\partial Z_\ell} \in \mathbb{R}^{d_\ell \times m}$. The backpropagation recurrence is

    $$
    \delta_\ell = (W_{\ell+1}^T \delta_{\ell+1}) \odot \sigma_\ell'(Z_\ell),
    $$

    where $\odot$ is the elementwise (Hadamard) product. The weight gradient is

    $$
    \frac{\partial \mathcal{L}}{\partial W_\ell} = \frac{1}{m}\delta_\ell A_{\ell-1}^T.
    $$

    The core operations remain matrix multiplications: $W_{\ell+1}^T \delta_{\ell+1}$ (dimensions $d_\ell \times m$) and $\delta_\ell A_{\ell-1}^T$ (dimensions $d_\ell \times d_{\ell-1}$).

??? proof "Proof"
    By the chain rule, $\frac{\partial \mathcal{L}}{\partial Z_\ell} = \frac{\partial A_\ell}{\partial Z_\ell} \cdot \frac{\partial Z_{\ell+1}}{\partial A_\ell} \cdot \frac{\partial \mathcal{L}}{\partial Z_{\ell+1}}$. Since $A_\ell = \sigma_\ell(Z_\ell)$ is elementwise, $\frac{\partial A_\ell}{\partial Z_\ell}$ is diagonal (in the elementwise sense), contributing $\sigma_\ell'(Z_\ell)$. From $Z_{\ell+1} = W_{\ell+1}A_\ell + \mathbf{b}_{\ell+1}\mathbf{1}^T$, we get $\frac{\partial Z_{\ell+1}}{\partial A_\ell} = W_{\ell+1}$, giving $\delta_\ell = (W_{\ell+1}^T\delta_{\ell+1})\odot\sigma_\ell'(Z_\ell)$. For the weights, $\frac{\partial \mathcal{L}}{\partial W_\ell} = \delta_\ell \frac{\partial Z_\ell}{\partial W_\ell}$, and from $Z_\ell = W_\ell A_{\ell-1}$, we get $\frac{\partial \mathcal{L}}{\partial W_\ell} = \frac{1}{m}\delta_\ell A_{\ell-1}^T$. $\blacksquare$

!!! theorem "Theorem 29.11 (Gradient vanishing/exploding and singular values)"
    For linear networks ($\sigma_\ell = \text{id}$), gradient propagation from layer $L$ to layer $\ell$ involves the matrix product $\prod_{k=\ell+1}^L W_k^T$. Let $\sigma_{\max}(W_k)$ and $\sigma_{\min}(W_k)$ denote the largest and smallest singular values of $W_k$. Then:

    1. If $\sigma_{\max}(W_k) > 1$ for most $k$, gradients may **explode exponentially**.
    2. If $\sigma_{\max}(W_k) < 1$ for most $k$, gradients may **vanish exponentially**.
    3. When $W_k$ is orthogonal (all singular values equal 1), gradient norms are preserved — this is the theoretical basis for orthogonal initialization and orthogonal-constrained training.

??? proof "Proof"
    For linear networks, $\delta_\ell = \prod_{k=\ell+1}^L W_k^T \cdot \delta_L$. By submultiplicativity of singular values, $\|\delta_\ell\| \le \prod_{k=\ell+1}^L \sigma_{\max}(W_k) \|\delta_L\|$. Similarly, $\|\delta_\ell\| \ge \prod_{k=\ell+1}^L \sigma_{\min}(W_k) \|\delta_L\|$. When all $W_k$ are orthogonal, $\sigma_{\max} = \sigma_{\min} = 1$, so $\|\delta_\ell\| = \|\delta_L\|$. $\blacksquare$

!!! example "Example 29.8"
    **Gradient propagation analysis.** Consider a 5-layer linear network with weight matrices $W_k \in \mathbb{R}^{100 \times 100}$.

    **Case 1**: $W_k = 1.1 I$. After 4 layers, gradients are amplified by $1.1^4 \approx 1.46$. After 50 layers, amplification is $1.1^{50} \approx 117$ — gradient explosion.

    **Case 2**: $W_k = 0.9 I$. After 50 layers, gradients decay to $0.9^{50} \approx 0.005$ — gradient vanishing.

    **Case 3**: $W_k$ is a random orthogonal matrix. Regardless of depth, $\|\delta_\ell\| = \|\delta_L\|$, perfectly preserving gradients. This explains why orthogonal initialization (e.g., Saxe et al., 2014) significantly improves deep network training.

---

## 29.8 Matrix Factorization in Recommender Systems

<div class="context-flow" markdown>

**Core idea**: User-item rating matrix $R$ is sparse and incomplete → low-rank approximation $R \approx UV^T$ → NMF yields interpretable nonnegative factors → regularization prevents overfitting

**Link**: Ch11 SVD · Ch10 matrix decompositions · Ch17 nonnegative matrices

</div>

Recommender systems represent one of the most successful industrial applications of matrix factorization, with the core idea being collaborative filtering through the low-rank structure of rating matrices.

!!! definition "Definition 29.12 (Matrix factorization model)"
    Let $R \in \mathbb{R}^{m \times n}$ be the **user-item rating matrix**, where $R_{ij}$ is user $i$'s rating for item $j$ (with many missing values). The **matrix factorization model** seeks a low-rank approximation

    $$
    R \approx U V^T, \quad U \in \mathbb{R}^{m \times k}, \quad V \in \mathbb{R}^{n \times k},
    $$

    where $k \ll \min(m,n)$ is the number of **latent factors**. The $i$-th row $\mathbf{u}_i^T$ of $U$ is user $i$'s latent vector, the $j$-th row $\mathbf{v}_j^T$ of $V$ is item $j$'s latent vector, and the predicted rating is $\hat{R}_{ij} = \mathbf{u}_i^T\mathbf{v}_j$.

!!! definition "Definition 29.13 (Regularized matrix factorization)"
    Regularized matrix factorization solves the optimization problem:

    $$
    \min_{U, V} \sum_{(i,j) \in \Omega} (R_{ij} - \mathbf{u}_i^T\mathbf{v}_j)^2 + \lambda(\|U\|_F^2 + \|V\|_F^2),
    $$

    where $\Omega$ is the index set of observed ratings and $\lambda > 0$ is the regularization parameter. A common algorithm is **Alternating Least Squares** (ALS):

    - Fix $V$, update each user: $\mathbf{u}_i = (V_{\Omega_i}^T V_{\Omega_i} + \lambda I)^{-1}V_{\Omega_i}^T \mathbf{r}_i$.
    - Fix $U$, update each item: $\mathbf{v}_j = (U_{\Omega_j}^T U_{\Omega_j} + \lambda I)^{-1}U_{\Omega_j}^T \mathbf{r}_j$.

    Here $\Omega_i$ is the set of items rated by user $i$, and $V_{\Omega_i}$ is the submatrix of corresponding rows.

!!! theorem "Theorem 29.12 (SVD and matrix completion)"
    If the complete rating matrix $R$ has exactly rank $k$, and the observation set $\Omega$ satisfies certain conditions (e.g., uniform random sampling with $|\Omega| \ge C \cdot (m+n)k \log^2(m+n)$), then $R$ can be exactly recovered via nuclear norm minimization:

    $$
    \min_{M} \|M\|_* \quad \text{s.t.} \quad M_{ij} = R_{ij}, \quad (i,j) \in \Omega,
    $$

    where $\|M\|_* = \sum_j \sigma_j(M)$ is the nuclear norm (sum of singular values), the convex relaxation of the rank function.

??? proof "Proof"
    (Proof sketch) The nuclear norm is the convex envelope of the rank function on $\{M : \|M\|_{\text{op}} \le 1\}$. When $R$ satisfies the incoherence condition and $|\Omega|$ is large enough, the solution to nuclear norm minimization is unique and equals $R$. This result was proved by Candes and Recht (2009), with the key tools being matrix concentration inequalities and dual certificate construction. The complete proof requires tools from random matrix theory (Ch23). $\blacksquare$

!!! definition "Definition 29.14 (Nonnegative matrix factorization, NMF)"
    **Nonnegative Matrix Factorization** (NMF) requires nonnegative factor matrices:

    $$
    R \approx WH, \quad W \in \mathbb{R}_{\ge 0}^{m \times k}, \quad H \in \mathbb{R}_{\ge 0}^{k \times n}.
    $$

    The nonnegativity constraint gives the factorization a "parts-based" interpretability: each latent factor represents a "topic" or "pattern," and users/items are represented as nonnegative weighted combinations of these topics. The common multiplicative update rules are

    $$
    W \leftarrow W \odot \frac{RH^T}{WHH^T}, \quad H \leftarrow H \odot \frac{W^TR}{W^TWH},
    $$

    where division is elementwise.

!!! theorem "Theorem 29.13 (Convergence of NMF multiplicative updates)"
    The Lee-Seung multiplicative update rules ensure that the objective $\|R - WH\|_F^2$ is monotonically nonincreasing. Specifically, define the auxiliary function

    $$
    G(h, h') = F(h') + F'(h')(h - h') + \frac{(WHH^T)_{ij}}{h'_{ij}}(h - h')^2.
    $$

    Then $G(h, h') \ge F(h)$ and $G(h', h') = F(h')$. The multiplicative update is equivalent to $h^{(t+1)} = \arg\min_h G(h, h^{(t)})$, and the auxiliary function property guarantees $F(h^{(t+1)}) \le F(h^{(t)})$.

??? proof "Proof"
    Differentiating $\|R - WH\|_F^2$ with respect to $H_{ij}$: $\frac{\partial F}{\partial H_{ij}} = -2(W^TR)_{ij} + 2(W^TWH)_{ij}$. The multiplicative update can be written as

    $$
    H_{ij} \leftarrow H_{ij} \cdot \frac{(W^TR)_{ij}}{(W^TWH)_{ij}}.
    $$

    Define the auxiliary function (a quadratic upper bound) tangent at $h = h'$. Minimizing the auxiliary function yields the update rule. Since each step does not increase the auxiliary function value, and the auxiliary function is an upper bound on the objective, the objective is monotonically nonincreasing. A bounded decreasing sequence converges. $\blacksquare$

!!! example "Example 29.9"
    **Simple matrix factorization example.** Let 3 users rate 4 movies (? denotes missing):

    $$
    R = \begin{pmatrix}5&3&?&1\\4&?&?&1\\1&1&?&5\end{pmatrix}.
    $$

    Set latent factors $k = 2$. Initialize $U \in \mathbb{R}^{3 \times 2}$, $V \in \mathbb{R}^{4 \times 2}$. After ALS iterations, suppose convergence to

    $$
    U \approx \begin{pmatrix}2.1&0.3\\1.8&0.2\\0.5&2.0\end{pmatrix}, \quad V \approx \begin{pmatrix}2.3&0.2\\1.4&0.1\\0.8&0.5\\0.1&2.4\end{pmatrix}.
    $$

    Predicted missing ratings: $\hat{R}_{13} = \mathbf{u}_1^T\mathbf{v}_3 \approx 2.1 \times 0.8 + 0.3 \times 0.5 = 1.83$, $\hat{R}_{22} = 1.8 \times 1.4 + 0.2 \times 0.1 = 2.54$. The latent factors are interpretable: the first dimension might correspond to "action movie preference" and the second to "art-house movie preference."
