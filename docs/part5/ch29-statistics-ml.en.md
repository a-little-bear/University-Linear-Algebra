# Chapter 29: Linear Algebra in Statistics and Machine Learning

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues & SVD (Ch06, 11) · Positive Definite Matrices (Ch16) · Least Squares (Ch07) · Matrix Decompositions (Ch10)

**Chapter Outline**: Multivariate Statistics & Covariance Matrices → Principal Component Analysis (PCA) → Linear Regression & Normal Equations → Regularization (Ridge and LASSO) → Linear Discriminant Analysis (LDA) → Kernel Methods & Gram Matrices → Support Vector Machines (SVM) → Neural Networks (Forward/Backward Pass) → Matrix Factorization in Recommender Systems

**Extension**: Linear algebra is the "operating system" of machine learning; from the covariance analysis of classical statistics to the backpropagation of deep learning, matrices are the primary objects of computation.

</div>

Machine learning is essentially large-scale, high-dimensional geometry. Every data sample is a vector, every dataset is a matrix, and every algorithm is a series of matrix transformations. Whether it's reducing the dimensions of a massive dataset using SVD or calculating the gradients of a neural network using chain rules, linear algebra provides the unified language and efficient computational engine. This chapter explores how these mathematical structures power modern AI.

---

## 29.1 Foundations: Covariance and PCA

<div class="context-flow" markdown>

**Matrix Perspective**: $n$ samples, $p$ features → Data matrix $X \in \mathbb{R}^{n \times p}$. The **Sample Covariance Matrix** $\Sigma$ captures the relationships between features.

</div>

!!! definition "Definition 29.1 (Covariance Matrix)"
    For a centered data matrix $X_c = X - \mathbf{1}\bar{x}^T$, the sample covariance matrix is:
    $$\Sigma = \frac{1}{n-1} X_c^T X_c$$
    $\Sigma$ is symmetric positive semi-definite. Its diagonal entries are the variances of the features.

!!! technique "Principal Component Analysis (PCA)"
    PCA seeks to project data onto a lower-dimensional subspace while preserving the maximum variance.
    - The directions of maximum variance are the **eigenvectors** of $\Sigma$.
    - This is equivalent to taking the **SVD** of the centered data matrix $X_c$.

---

## 29.2 Linear Regression and Regularization

!!! theorem "Theorem 29.1 (The Normal Equation)"
    The Ordinary Least Squares (OLS) solution to $X\beta = y$ is:
    $$\hat{\beta} = (X^T X)^{-1} X^T y$$
    This projects $y$ onto the column space of $X$.

!!! definition "Definition 29.2 (Ridge and LASSO)"
    - **Ridge Regression (L2)**: $\min \|y - X\beta\|^2 + \lambda \|\beta\|^2$. The solution is $(X^T X + \lambda I)^{-1} X^T y$, which is more stable for ill-conditioned $X$.
    - **LASSO (L1)**: $\min \|y - X\beta\|^2 + \lambda \|\beta\|_1$. This promotes sparsity, setting some coefficients exactly to zero.

---

## 29.3 Kernel Methods and SVMs

<div class="context-flow" markdown>

**The Kernel Trick**: Non-linear mapping $\phi(x)$ into a high-dimensional space can be handled implicitly using the inner product $K(x, y) = \langle \phi(x), \phi(y) \rangle$.

</div>

!!! theorem "Theorem 29.2 (Mercer's Theorem)"
    A symmetric function $K(x, y)$ can be represented as an inner product in some Hilbert space if and only if its **Gram matrix** $K_{ij} = K(x_i, x_j)$ is positive semi-definite for any set of points.

---

## 29.4 Neural Networks and Backpropagation

!!! technique "Matrix Calculus in Deep Learning"
    - **Forward Pass**: A layer is a matrix-vector product followed by a non-linear activation: $a^{(l)} = \sigma(W^{(l)} a^{(l-1)} + b^{(l)})$.
    - **Backward Pass**: Gradient updates involve the **Jacobian** of the loss function. The error signal flows back through the chain rule as a sequence of matrix-transpose multiplications: $\delta^{(l)} = ((W^{(l+1)})^T \delta^{(l+1)}) \odot \sigma'(z^{(l)})$.

---

## Exercises


****
??? success "Solution"
     Because $\Sigma = \frac{1}{n-1} X_c^T X_c$. For any $v$, $v^T \Sigma v = \frac{1}{n-1} \|X_c v\|^2 \ge 0$.


****
??? success "Solution"
     It is the direction of the major axis of the hyper-ellipsoid formed by the data distribution.


****
??? success "Solution"
     $X^T X$ becomes singular (not invertible). One must use the Moore-Penrose pseudoinverse $X^+$ or add regularization (Ridge).


****
??? success "Solution"
    Its Gram matrix is $X X^T$ (for data matrix $X$), which is always PSD.


****
??? success "Solution"
     $10 \times 100$.


****
??? success "Solution"
     If the singular values of the weight matrices are all significantly less than 1, the product of Jacobians in backpropagation will shrink exponentially, causing the gradient to vanish.


****
??? success "Solution"
     Original: $10^{11}$ entries. Factorized: $(10^6 + 10^5) \times 100 \approx 1.1 \times 10^8$. Compression is roughly 1000x.


****
??? success "Solution"
     It solves a generalized eigenvalue problem $S_B w = \lambda S_W w$, where $S_B$ is the between-class scatter and $S_W$ is the within-class scatter.


****
??? success "Solution"
     $d = \sqrt{(x-\mu)^T \Sigma^{-1} (x-\mu)}$. It scales distances by the inverse of the data's variance.

****
??? success "Solution"
    ## Chapter Summary

Statistics and Machine Learning are the practical realizations of high-dimensional matrix analysis:


****: PCA and SVD show that high-dimensional data often live on lower-dimensional "manifolds," and linear algebra is the tool to find them.

****: Identified OLS and SVM as specific linear systems or quadratic programs, proving that training AI is a form of matrix optimization.

****: Explained Neural Networks as matrix chains, linking the stability of deep learning to the spectral properties of weight matrices.

****: Used Kernel methods to prove that linear algebra can solve non-linear problems by working in the right feature space.
