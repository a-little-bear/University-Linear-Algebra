# Chapter 11: Singular Value Decomposition

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Orthogonality (Ch7) · Matrix Factorization (Ch10) · Rank (Ch2)

**Chapter Outline**: Definition of SVD ($A = U \Sigma V^T$) → Existence and Uniqueness → Singular Values vs. Eigenvalues → Geometric Interpretation (Rotating and Stretching) → Low-Rank Approximation (Eckart-Young Theorem) → Matrix Norms and SVD → Image Compression → Principal Component Analysis (PCA) Overview → Pseudoinverse via SVD

**Extension**: SVD is the "Swiss Army Knife" of linear algebra; it provides the definitive answers to the rank, condition number, and fundamental subspaces of any matrix.

</div>

The **Singular Value Decomposition** (SVD) is the most powerful and versatile matrix factorization. While diagonalization ($A = PDP^{-1}$) only applies to specific square matrices, SVD ($A = U \Sigma V^T$) exists for *any* matrix, rectangular or square. It decomposes any linear map into a rotation ($V^T$), followed by a coordinate-wise scaling ($\Sigma$), followed by another rotation ($U$). This decomposition reveals the "energy" of the transformation in the singular values and allows for the optimal approximation of high-dimensional data by lower-dimensional structures.

---

## 11.1 The SVD Framework

!!! definition "Definition 11.1 (Singular Value Decomposition)"
    For any $m \times n$ matrix $A$, the SVD is $A = U \Sigma V^T$, where:
    - $U$ is an $m \times m$ orthogonal matrix (left singular vectors);
    - $V$ is an $n \times n$ orthogonal matrix (right singular vectors);
    - $\Sigma$ is an $m \times n$ diagonal matrix with non-negative entries $\sigma_1 \ge \sigma_2 \ge \dots \ge 0$ (singular values).

!!! theorem "Theorem 11.1 (Eckart-Young Theorem)"
    The best rank-$k$ approximation of $A$ in the Frobenius norm is obtained by keeping the $k$ largest singular values and setting the rest to zero:
    $$A_k = \sum_{i=1}^k \sigma_i u_i v_i^T$$

---

## Exercises

1. **[Fundamentals] Compute the singular values of $A = \begin{pmatrix} 3 & 0 \\ 0 & -2 \end{pmatrix}$.**
   ??? success "Solution"
       Singular values are the square roots of the eigenvalues of $A^T A = \begin{pmatrix} 9 & 0 \\ 0 & 4 \end{pmatrix}$. They are $\sigma_1 = 3, \sigma_2 = 2$. Note that singular values are always non-negative.

2. **[Geometry] Describe the image of the unit circle under $A = U \Sigma V^T$.**
   ??? success "Solution"
       It is an ellipse in $\mathbb{R}^m$ whose semi-axes have lengths $\sigma_i$ and are oriented along the directions of $u_i$.

3. **[Rank] How does SVD identify the rank of a matrix?**
   ??? success "Solution"
       The rank of $A$ is exactly the number of non-zero singular values. In practice, we count the number of singular values above a small numerical threshold.

4. **[Eigenvalues] Relate the singular values of $A$ to the eigenvalues of $A^T A$ and $A A^T$.**
   ??? success "Solution"
       $\sigma_i(A) = \sqrt{\lambda_i(A^T A)} = \sqrt{\lambda_i(A A^T)}$. The right singular vectors are eigenvectors of $A^T A$, and left singular vectors are eigenvectors of $A A^T$.

5. **[Norm] Find the spectral norm $\|A\|_2$ using SVD.**
   ??? success "Solution"
       $\|A\|_2 = \max \sigma_i = \sigma_1$. The spectral norm is determined by the most dominant scaling factor of the operator.

6. **[Pseudo-inverse] Compute the Moore-Penrose inverse $A^\dagger$ using SVD.**
   ??? success "Solution"
       $A^\dagger = V \Sigma^\dagger U^T$, where $\Sigma^\dagger$ is the diagonal matrix formed by taking the reciprocal of non-zero singular values.

7. **[Compression] How is SVD used in image compression?**
   ??? success "Solution"
       By treating an image as a matrix and keeping only the first $k$ terms of the SVD. Since the singular values of natural images decay rapidly, a small $k$ can capture most of the visual information.

8. **[Condition Number] Define the condition number of a matrix in terms of singular values.**
   ??? success "Solution"
       $\kappa(A) = \sigma_{\max} / \sigma_{\min}$. It measures the sensitivity of the system $Ax = b$ to perturbations.

9. **[Symmetry] If $A$ is symmetric and positive definite, how do its singular values relate to its eigenvalues?**
   ??? success "Solution"
       In this case, singular values and eigenvalues are identical: $\sigma_i = \lambda_i$. If $A$ is symmetric but not positive definite, $\sigma_i = |\lambda_i|$.

10. **[PCA] What is the connection between SVD and Principal Component Analysis?**
    ??? success "Solution"
        PCA is performed by taking the SVD of the zero-centered data matrix $X$. The right singular vectors are the principal components, and the singular values are proportional to the standard deviation along those components.

## Chapter Summary

This chapter establishes the definitive structural decomposition of any linear operator:

1. **Spectral Generality**: Positioned SVD as the universal counterpart to diagonalization, applicable to all rectangular matrices.
2. **Optimal Compression**: Utilized the Eckart-Young theorem to establish SVD as the foundation for high-fidelity data approximation.
3. **Geometric Analysis**: Linked singular values to the distortion of unit spheres into ellipses, revealing the geometric power of the operator.
4. **Numerical Benchmarking**: Formulated rank, condition number, and norm analysis through the lens of singular value spectra.
