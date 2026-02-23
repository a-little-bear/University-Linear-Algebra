# Chapter 10: Matrix Decompositions

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch01) · Matrix Algebra (Ch02) · Orthogonality (Ch07) · Eigenvalues (Ch06)

**Chapter Outline**: Motivation for Decompositions (Divide and Conquer & Numerical Stability) → LU Decomposition (The Matrix Expression of Gaussian Elimination) → Cholesky Decomposition (The Privilege of Symmetric Positive Definite Matrices) → QR Decomposition (Algorithmic Implementation of Orthogonalization) → Schur Decomposition (Unitary Similarity and Spectral Theory) → LDL^T and LDU Decompositions → Polar Decomposition → Rank Factorization → Numerical Algorithm Comparison: Householder Transforms vs. Givens Rotations → Applications: Industrial-grade Linear Solvers and Efficient Least Squares Solutions

**Extension**: Matrix decomposition is the soul of Numerical Linear Algebra (Ch22); it proves that complex operators can be "disassembled" into components with well-defined properties (such as triangular or orthogonal). Almost all scientific computing software (LAPACK, MATLAB) is built around these decomposition algorithms.

</div>

Matrix decomposition (also known as matrix factorization) is the process of expressing a matrix as a product of several matrices with specific structures. This not only simplifies theoretical proofs but is also key to improving numerical efficiency and enhancing stability. This chapter systematically introduces the most commonly used decomposition methods in linear algebra, exploring the trade-offs between computational cost and precision.

---

## 10.1 LU Decomposition: The Algebra of Elimination

!!! definition "Definition 10.1 (LU Decomposition)"
    Factors a square matrix $A$ into the product of a lower triangular matrix $L$ and an upper triangular matrix $U$: $A = LU$.
    - **Existence**: If all leading principal minors of $A$ are non-zero, then $A$ has a unique LU decomposition (with $L$ having ones on the diagonal).
    - **Significance**: It transforms solving $Ax = b$ into solving two extremely simple triangular systems: $Ly = b$ (forward substitution) and $Ux = y$ (backward substitution).

---

## 10.2 Cholesky Decomposition: Efficiency for SPD Matrices

!!! theorem "Theorem 10.1 (Cholesky Decomposition)"
    If $A$ is a real symmetric positive definite matrix, there exists a unique lower triangular matrix $L$ with positive diagonal entries such that:
    $$A = LL^T$$
    **Advantages**: It requires only half the computational effort of LU and is numerically extremely stable, requiring no pivoting.

---

## 10.3 QR Decomposition: The Orthogonality Standard

!!! theorem "Theorem 10.2 (QR Decomposition)"
    Every matrix $A$ with linearly independent columns can be factored as $A = QR$.
    - $Q$ is an orthogonal (or unitary) matrix.
    - $R$ is an upper triangular invertible matrix.
    **Implementation**: Beyond the Gram-Schmidt process, industrial applications favor the **Householder Transformation** for its superior ability to maintain orthogonality.

---

## 10.4 Schur Decomposition: The Universal Spectral Form

!!! theorem "Theorem 10.3 (Schur Decomposition)"
    Every complex square matrix $A$ is unitarily similar to an upper triangular matrix $T$:
    $$U^* A U = T$$
    where the diagonal entries of $T$ are the eigenvalues of $A$.
    **Significance**: This proves that any operator can be represented in a "layered" form in an appropriate basis, serving as the most general form of eigenvalue theory.

---

## Exercises

**1. [LU] Find the LU decomposition of $A = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. Perform elimination: $R_2 \gets R_2 - 2R_1$.
    2. Record the multiplier in $L$: $l_{21} = 2$.
    3. The elimination result is $U$: $U = \begin{pmatrix} 1 & 2 \\ 0 & 5-4 \end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$.
    **Conclusion**: $L = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}, U = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$.

**2. [Cholesky] Find the Cholesky decomposition for the matrix $A$ from the previous problem.**

??? success "Solution"
    **Analysis:**
    Since the LU decomposition resulted in $U = L^T$, the matrix is positive definite.
    **Conclusion**: $L = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}$. Verify: $LL^T = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix} \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix} = A$.

**3. [QR] What is the QR decomposition of $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$?**

??? success "Solution"
    **Observation:**
    1. The matrix is already upper triangular.
    2. Check columns: $\mathbf{a}_1 = (1, 0)^T, \mathbf{a}_2 = (1, 1)^T$. They already satisfy the structural characteristics of an orthonormal basis setup (first is a unit vector).
    **Conclusion**: Since $A$ is upper triangular and normalized in the first column, we can take $Q=I$ and $R=A$ (if we don't strictly require $Q$ to be $I$, G-S yields this directly).

**4. [Schur] Prove that the trace of $T$ in a Schur decomposition equals the trace of $A$.**

??? success "Solution"
    **Proof:**
    1. Trace has the cyclic property: $\operatorname{tr}(XYZ) = \operatorname{tr}(ZXY)$.
    2. $\operatorname{tr}(T) = \operatorname{tr}(U^* A U) = \operatorname{tr}(A U U^*) = \operatorname{tr}(A I) = \operatorname{tr}(A)$.
    **Significance**: This explains why the sum of eigenvalues equals the sum of diagonal entries.

**5. [Existence] Give an example of a non-singular matrix that has no LU decomposition.**

??? success "Solution"
    **Example:**
    $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.
    **Analysis**: The first leading principal minor is 0. During elimination, we cannot use the first row to clear the element in the second row because the pivot is zero. Although invertible, it requires row swapping (PLU decomposition).

**6. [Pivoting] What is PLU decomposition and what problem does it solve?**

??? success "Solution"
    **Definition**: It uses a permutation matrix $P$ to swap rows such that $PA = LU$.
    **Function**:
    1. Solves the existence problem caused by zero pivots.
    2. Mitigates numerical instability caused by very small pivots (see Ch15). This is the default algorithm for linear solves in standard libraries.

**7. [Polar] In the polar decomposition $A = UP$, what properties does $P$ satisfy and what is its physical meaning?**

??? success "Solution"
    **Properties**: $P$ is positive semi-definite (positive definite if $A$ is invertible).
    **Physical Meaning**: Analogous to the polar form of complex numbers $z = re^{i\theta}$.
    - $U$ (unitary) represents a **pure rotation**.
    - $P$ (positive definite) represents a **pure stretch** along principal axes.
    This decomposition is vital in continuum mechanics for describing deformation.

**8. [Rank] If $\operatorname{rank}(A) = r$, what is the benefit of factoring it as $A = FG$ ($F$ is $m \times r, G$ is $r \times n$)?**

??? success "Solution"
    **Value:**
    1. **Storage Compression**: For large low-rank matrices, storing $F$ and $G$ requires $(m+n)r$ entries, which is much less than $mn$.
    2. **Acceleration**: Computing $Ax$ becomes $F(Gx)$, reducing complexity from $O(mn)$ to $O((m+n)r)$.

**9. [Stability] Why is Householder preferred over Gram-Schmidt?**

??? success "Solution"
    **Analysis:**
    1. Classical Gram-Schmidt (CGS) rapidly loses orthogonality between vectors due to floating-point rounding errors.
    2. Householder transformations use reflection matrices where orthogonality is structurally guaranteed by the matrix products, resulting in minimal accumulated error.
    **Conclusion**: Householder is the gold standard for robust QR decomposition.

**10. [Application] How do decompositions speed up solving $Ax=b$?**

??? success "Solution"
    **Process:**
    1. Perform a one-time high-cost decomposition (e.g., $O(n^3)$ for LU).
    2. For different $b$ vectors, only $O(n^2)$ forward and backward substitutions are needed.
    3. In engineering (like Finite Element Analysis), the matrix $A$ is often fixed while the load $b$ changes; this pre-factorization strategy greatly improves simulation throughput.

## Chapter Summary

Matrix decompositions are the algebraic art of "simplifying the complex":

1.  **Structural Mapping**: LU, QR, and other factorizations map general matrices into subgroups with desirable properties (triangular, orthogonal), establishing shortcuts for calculation.
2.  **Numerical Foundation**: Decomposition is not just a tool for proofs but the bulwark of numerical stability analysis, distinguishing "robust" from "fragile" algorithms.
3.  **Information Extraction**: Schur and Rank decompositions demonstrate how essential information like eigenvalues and rank can be "squeezed" into the diagonal or a few rows through specific product forms.
