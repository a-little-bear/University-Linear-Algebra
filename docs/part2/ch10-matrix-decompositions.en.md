# Chapter 10: Matrix Decompositions

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch01) · Matrix Algebra (Ch02) · Orthogonality (Ch07) · Eigenvalues (Ch06)

**Chapter Outline**: Motivation for Decompositions → LU Decomposition (Gaussian Elimination) → Cholesky Decomposition (Symmetric Positive Definite) → QR Decomposition (Orthogonal Projections) → Schur Decomposition (Unitary Similarity) → LDL^T and LDU Decompositions → Polar Decomposition → Rank Factorization → Numerical Algorithms and Stability Analysis

**Extension**: Matrix decompositions are the soul of Numerical Linear Algebra (Ch22); almost all industrial-grade linear algebra libraries (LAPACK, BLAS) are built around these decomposition algorithms.

</div>

Matrix decomposition (also known as matrix factorization) is the process of expressing a matrix as a product of several matrices with specific structures (e.g., triangular, orthogonal, diagonal). This not only simplifies theoretical proofs but is also the key to improving numerical efficiency and enhancing stability. This chapter systematically introduces the most commonly used decomposition methods in linear algebra.

---

## 10.1 LU Decomposition

!!! definition "Definition 10.1 (LU Decomposition)"
    Factors a square matrix $A$ into the product of a lower triangular matrix $L$ and an upper triangular matrix $U$: $A = LU$.
    - **Existence**: If all leading principal minors of $A$ are non-zero, then $A$ has a unique LU decomposition (with $L$ having ones on the diagonal).
    - **Application**: Solving $Ax = b$ becomes solving $Ly = b$ and $Ux = y$, greatly increasing efficiency for systems with multiple right-hand sides.

---

## 10.2 Cholesky Decomposition

!!! theorem "Theorem 10.1 (Cholesky Decomposition)"
    If $A$ is a real symmetric positive definite matrix, there exists a unique lower triangular matrix $L$ with positive diagonal entries such that:
    $$A = LL^T$$
    **Significance**: Cholesky decomposition is more stable than LU and requires only half the computational effort.

---

## 10.3 QR Decomposition

!!! theorem "Theorem 10.2 (QR Decomposition)"
    Every matrix $A$ with linearly independent columns can be factored as $A = QR$.
    - $Q$ is an orthogonal (or unitary) matrix.
    - $R$ is an upper triangular invertible matrix.
    **Implementation**: Beyond the Gram-Schmidt process, it can be realized via Householder reflections or Givens rotations, which offer better numerical stability.

---

## 10.4 Schur Decomposition

!!! theorem "Theorem 10.3 (Schur Decomposition)"
    Every complex square matrix $A$ is unitarily similar to an upper triangular matrix $T$:
    $$U^* A U = T$$
    where $U$ is unitary and the diagonal entries of $T$ are the eigenvalues of $A$.
    **Extension**: This is a generalization of the Spectral Theorem, proving the universality of unitary transformations in revealing spectral structures.

---

## Exercises

1. **[LU] Find the LU decomposition of $A = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}$.**

   ??? success "Solution"
       $L = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}, U = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$.

2. **[Cholesky] Find the Cholesky decomposition for the same matrix $A$.**

   ??? success "Solution"
       Since $U = L^T$ in the previous problem, $L = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}$ and $A = LL^T$.

3. **[QR] What is the QR decomposition of $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$?**

   ??? success "Solution"
       Since it is already upper triangular and its columns are orthogonal (though not normalized), $Q=I$ and $R=A$ (if we don't require normalization) or $Q=I, R=A$ is fine.

4. **[Schur] Prove that the trace of $T$ in a Schur decomposition equals the trace of $A$.**

   ??? success "Solution"
       Trace is a similarity invariant: $\operatorname{tr}(T) = \operatorname{tr}(U^* A U) = \operatorname{tr}(A U U^*) = \operatorname{tr}(A)$.

5. **[Existence] Give an example of a non-singular matrix that has no LU decomposition.**

   ??? success "Solution"
       $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$. The first leading principal minor is 0, making initial elimination impossible.

6. **[Pivoting] What is PLU decomposition?**

   ??? success "Solution"
       It uses a permutation matrix $P$ to swap rows such that $PA = LU$. This handles zero minors and improves numerical stability.

7. **[Polar] In the polar decomposition $A = UP$, what property does $P$ satisfy?**

   ??? success "Solution"
       $P$ is a positive semi-definite matrix (geometrically representing scaling), and $U$ is orthogonal (representing rotation).

8. **[Rank] If $\operatorname{rank}(A) = r$, what is the benefit of factoring it as $A = FG$ ($F$ is $m \times r, G$ is $r \times n$)?**

   ??? success "Solution"
       It compresses a large matrix into two skinny matrices, saving storage and simplifying operations.

9. **[Stability] Why is Householder preferred over Gram-Schmidt?**

   ??? success "Solution"
       Householder uses orthogonal reflections, avoiding the loss of orthogonality that occurs in Gram-Schmidt due to rounding errors.

10. **[Application] How do decompositions speed up solving $Ax=b$?**

   ??? success "Solution"
        By performing a one-time high-cost decomposition (like LU at $O(n^3)$), subsequent solves are reduced to low-cost triangular substitutions ($O(n^2)$).

## Chapter Summary

Matrix decompositions are the algebraic art of "simplifying the complex":

1.  **Structural Mapping**: Decompositions like LU and QR map general matrices into subgroups with excellent properties (triangular, orthogonal), establishing shortcuts for calculation.
2.  **Numerical Foundation**: Factorization is not just a tool for theoretical proofs but also the bulwark of numerical stability analysis, distinguishing "robust" from "fragile" algorithms.
3.  **Information Extraction**: Schur and Rank decompositions demonstrate how to "squeeze" essential information like eigenvalues and rank into the diagonal or a few rows through specific product forms.
