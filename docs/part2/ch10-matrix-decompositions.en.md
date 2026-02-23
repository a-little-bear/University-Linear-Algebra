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

****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ## Chapter Summary

Matrix decompositions are the algebraic art of "simplifying the complex":


****: Decompositions like LU and QR map general matrices into subgroups with excellent properties (triangular, orthogonal), establishing shortcuts for calculation.

****: Factorization is not just a tool for theoretical proofs but also the bulwark of numerical stability analysis, distinguishing "robust" from "fragile" algorithms.

****: Schur and Rank decompositions demonstrate how to "squeeze" essential information like eigenvalues and rank into the diagonal or a few rows through specific product forms.
