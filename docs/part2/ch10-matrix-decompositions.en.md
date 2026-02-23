# Chapter 10: Matrix Decompositions

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Gaussian Elimination (Ch1) · Eigenvalues (Ch6) · Orthogonality (Ch7)

**Chapter Outline**: LU Decomposition (Matrix-based Elimination) → Cholesky Decomposition (Positive Definite Matrices) → QR Decomposition (Orthogonalization) → Eigendecomposition → Schur Decomposition (Unitary Similarity) → Polar Decomposition → Numerical Stability of Decompositions

**Extension**: Matrix decomposition is the soul of numerical linear algebra, decomposing a "black box" operator into components with excellent physical or geometric properties.

</div>

Matrix decomposition is the engineering art of handling large-scale computational problems. Instead of operating directly on a complex matrix $A$, it is better to break it down into simple triangular, diagonal, or orthogonal matrices. Different decompositions correspond to different application perspectives: LU for solving equations, QR for stability, and eigendecomposition for the fundamental behavior of the operator.

---

## 10.1 Core Decomposition Models

!!! definition "Definition 10.1 (LU Decomposition)"
    If a square matrix $A$ can be written as $A = LU$, where $L$ is lower triangular and $U$ is upper triangular, it is called LU decomposition. This corresponds to Gaussian elimination without row swaps.

!!! theorem "Theorem 10.3 (Cholesky Decomposition)"
    Every symmetric positive definite matrix $A$ can be uniquely factored as $A = LL^T$, where $L$ is a lower triangular matrix with positive diagonal entries.

---

## Exercises

1. **[LU Decomposition] Perform LU decomposition on $A = \begin{pmatrix} 2 & 1 \\ 4 & 7 \end{pmatrix}$.**
   ??? success "Solution"
       $L = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}$, $U = \begin{pmatrix} 2 & 1 \\ 0 & 5 \end{pmatrix}$.
       Check: $LU = \begin{pmatrix} 2 & 1 \\ 4 & 2+5 \end{pmatrix} = A$.

2. **[Cholesky] Compute the Cholesky decomposition of $A = \begin{pmatrix} 4 & 12 \\ 12 & 37 \end{pmatrix}$.**
   ??? success "Solution"
       Let $L = \begin{pmatrix} l_{11} & 0 \\ l_{21} & l_{22} \end{pmatrix}$.
       $l_{11} = \sqrt{4} = 2$.
       $l_{21} = 12/2 = 6$.
       $l_{22} = \sqrt{37 - 6^2} = 1$.
       Thus $L = \begin{pmatrix} 2 & 0 \\ 6 & 1 \end{pmatrix}$.

3. **[QR Meaning] In the QR decomposition $A=QR$, what geometric property do the columns of $Q$ have?**
   ??? success "Solution"
       The columns of $Q$ form an orthonormal basis for the column space of $A$. They are obtained by performing Gram-Schmidt orthogonalization on the columns of $A$.

4. **[Eigendecomposition] If $A$ is diagonalizable, write its eigendecomposition form.**
   ??? success "Solution"
       $A = PDP^{-1}$, where $D$ is the diagonal matrix of eigenvalues and the columns of $P$ are the corresponding eigenvectors.

5. **[Schur Decomposition] Prove that the diagonal entries of the Schur decomposition $A = U T U^*$ are the eigenvalues of $A$.**
   ??? success "Solution"
       Since $A$ is similar to the upper triangular matrix $T$, they share the same characteristic polynomial. The eigenvalues of an upper triangular matrix are its diagonal entries.

6. **[Polar Decomposition] In the polar decomposition $A = UP$, what do $U$ and $P$ represent?**
   ??? success "Solution"
       $U$ is an isometry (rotation or reflection) and $P$ is a positive semi-definite Hermitian matrix (stretch). This is analogous to the polar form of a complex number $z = re^{i\theta}$.

7. **[Existence] Does every matrix have an LU decomposition?**
   ??? success "Solution"
       No. An LU decomposition exists only if all leading principal minors are non-zero. If row swaps are needed during elimination, a PLU decomposition (with permutation matrix $P$) is required.

8. **[Computational Cost] Is solving $Ax=b$ via LU decomposition faster than direct inversion?**
   ??? success "Solution"
       Yes. LU decomposition has complexity $O(n^3/3)$, and subsequent back-substitution takes only $O(n^2)$. Direct inversion is generally more computationally expensive and numerically less stable.

9. **[Orthogonal Similarity] What is special about the eigendecomposition if $A$ is symmetric?**
   ??? success "Solution"
       It can be written as $A = QDQ^T$, where $Q$ is an orthogonal matrix. This means symmetric matrices are orthogonally diagonalizable.

10. **[Application] Why is QR decomposition often used for least squares in computer vision?**
    ??? success "Solution"
        QR decomposition is numerically more stable (lower condition number) than directly using the normal equations $A^T A x = A^T b$, significantly reducing rounding errors in floating-point calculations.

## Chapter Summary

Matrix decomposition is the "anatomy" of linear algebra:

1. **Structure Revelation**: Reducing complex operators to basic geometric actions (rotation, projection, stretching).
2. **Computational Optimization**: Lowering $O(n^3)$ complexity to $O(n^2)$ for repeated use through triangulation or diagonalization.
3. **Theoretical Unity**: Different decomposition theorems (e.g., Schur, Spectral Theorem) delineate the boundaries of matrix theory.
