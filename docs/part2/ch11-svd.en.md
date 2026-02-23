# Chapter 11: Singular Value Decomposition (SVD)

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Orthogonality (Ch7) · Matrix Decompositions (Ch10)

**Chapter Outline**: Definition $A = U \Sigma V^*$ → Geometric Meaning (Hyperellipsoidal Mapping) → Relation between Singular Values and Eigenvalues → Truncated SVD and Low-Rank Approximation (Eckart-Young Theorem) → SVD for Generalized Inverses → Applications (Image Compression, PCA, Collaborative Filtering)

**Extension**: SVD is the most versatile decomposition in linear algebra, valid for rectangular and singular matrices alike, and is often called the "bedrock of linear algebra."

</div>

While eigendecomposition explores the intrinsic structure of square matrices, the **Singular Value Decomposition (SVD)** provides a comprehensive analysis of operators of any shape. It reveals that any linear mapping can be viewed as: rotation → anisotropic scaling → rotation. In the age of data science, SVD is the ultimate tool for extracting core features from high-dimensional data.

---

## 11.1 Definitions and Core Theorems

!!! definition "Definition 11.1 (SVD)"
    For any $m \times n$ complex matrix $A$, there exist an $m \times m$ unitary matrix $U$ and an $n \times n$ unitary matrix $V$ such that:
    $$A = U \Sigma V^*$$
    where $\Sigma$ is a diagonal-like matrix with non-negative real entries $\sigma_1 \ge \sigma_2 \dots \ge \sigma_r > 0$ on the diagonal.

!!! theorem "Theorem 11.3 (Eckart-Young Theorem)"
    The optimal rank-$k$ approximation of $A$ in the Frobenius norm is given by the truncated SVD:
    $$A_k = \sum_{i=1}^k \sigma_i u_i v_i^*$$

---

## Exercises

1. **[Basic Calculation] Calculate the singular values of $A = \begin{pmatrix} 3 & 0 \\ 0 & -2 \end{pmatrix}$.**
   ??? success "Solution"
       Singular values are the square roots of the eigenvalues of $A^* A$.
       $A^* A = \begin{pmatrix} 9 & 0 \\ 0 & 4 \end{pmatrix}$. Eigenvalues are 9 and 4.
       Thus $\sigma_1 = 3, \sigma_2 = 2$. Note singular values are always non-negative.

2. **[SVD and Rank] If $A$ has 3 non-zero singular values, what is its rank?**
   ??? success "Solution"
       The rank is 3. The number of non-zero singular values is exactly equal to the rank of the matrix.

3. **[Geometry] What shape does the linear transformation $A = U \Sigma V^*$ map the unit circle to?**
   ??? success "Solution"
       It maps it to a **hyperellipse**. The lengths of the semi-axes are determined by $\sigma_i$, and the directions are determined by the left singular vectors $u_i$.

4. **[Spectral Norm] Prove that the spectral norm $\|A\|_2$ equals its maximum singular value $\sigma_{\max}$.**
   ??? success "Solution"
       $\|A\|_2 = \max \frac{\|Ax\|}{\|x\|} = \max \sqrt{\frac{x^* A^* A x}{x^* x}}$.
       By the Rayleigh quotient, the maximum is $\sqrt{\lambda_{\max}(A^* A)} = \sigma_{\max}$.

5. **[Condition Number] Define the condition number $\kappa(A)$ using singular values.**
   ??? success "Solution"
       $\kappa(A) = \sigma_{\max} / \sigma_{\min}$. It measures how much the matrix amplifies errors in the inversion process.

6. **[Low-Rank Approx] If $A$ has singular values $100, 50, 1, 0.1$, what is the Frobenius norm error of its rank-2 approximation?**
   ??? success "Solution"
       The error is the square root of the sum of squares of the discarded singular values: $\sqrt{1^2 + 0.1^2} \approx 1.005$.

7. **[Eigenvalue Comparison] If $A$ is a positive definite symmetric matrix, how do its singular values relate to its eigenvalues?**
   ??? success "Solution"
       For a positive definite symmetric matrix, the singular values are exactly equal to the eigenvalues.

8. **[Calculation] Find the SVD of $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$.**
   ??? success "Solution"
       $A^* A = (2)$, eigenvalue 2, so $\sigma_1 = \sqrt{2}$.
       $V = (1)$. $u_1 = Av_1/\sigma_1 = \begin{pmatrix} 1/\sqrt{2} \\ 1/\sqrt{2} \end{pmatrix}$.
       $A = \begin{pmatrix} 1/\sqrt{2} \\ 1/\sqrt{2} \end{pmatrix} \begin{pmatrix} \sqrt{2} \end{pmatrix} (1)$.

9. **[Application] How does keeping the top $k$ singular values achieve compression in image processing?**
   ??? success "Solution"
       The original image requires $m \times n$ storage. Keeping the top $k$ singular values and their $u, v$ vectors requires only $k(m+n+1)$ space. When $k \ll \min(m,n)$, the compression ratio is massive.

10. **[Pseudoinverse] Write the Moore-Penrose pseudoinverse $A^\dagger$ using SVD.**
    ??? success "Solution"
        $A^\dagger = V \Sigma^\dagger U^*$, where $\Sigma^\dagger$ is obtained by inverting the non-zero singular values and transposing.

## Chapter Summary

SVD is the crown jewel of applied mathematics:

1. **Universality**: It breaks the limitations of square matrices, providing a unified view for all linear systems.
2. **Energy Concentration**: Through the distribution of singular values, it reveals the dominant components in data.
3. **Robustness**: As the core of low-rank approximation, SVD is the algebraic standard for noise filtering and info compression.
