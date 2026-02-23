# Chapter 12: Jordan Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Decompositions (Ch10) · Polynomial Algebra (Ch00)

**Chapter Outline**: Limitations of Diagonalization (Defective Matrices) → Definition of Jordan Blocks → Generalized Eigenvectors and Jordan Chains → Existence and Uniqueness of Jordan Canonical Form (JCF) → Relationship between Minimal and Characteristic Polynomials → Steps to Determine JCF (Rank Methods, Weyr Characteristic) → Jordan Analysis of Matrix Powers and Series → Numerical Instability

**Extension**: The Jordan Canonical Form is the ultimate representation under similarity transformations; it is indispensable for the theoretical analysis of linear differential equations (Ch26) and matrix functions (Ch13).

</div>

Not all square matrices can be diagonalized. When the geometric multiplicity of an eigenvalue is less than its algebraic multiplicity, the matrix is called "defective." The Jordan Canonical Form (JCF) provides the structure closest to diagonal for such matrices. It not only reveals the deep structure of linear operators but also serves as the most important theoretical tool in matrix analysis.

---

## 12.1 Jordan Blocks and Generalized Eigenvectors

!!! definition "Definition 12.1 (Jordan Block)"
    A **Jordan block** $J_k(\lambda)$ of order $k$ is a square matrix with $\lambda$ on the main diagonal, 1s on the super-diagonal, and 0s elsewhere:
    $$J_k(\lambda) = \begin{pmatrix} \lambda & 1 & 0 & \cdots & 0 \\ 0 & \lambda & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & 0 & \lambda & 1 \\ 0 & 0 & 0 & 0 & \lambda \end{pmatrix}$$

!!! definition "Definition 12.2 (Generalized Eigenvector)"
    A vector $\mathbf{v}$ is a **generalized eigenvector of rank $k$** corresponding to $\lambda$ if $(A - \lambda I)^k \mathbf{v} = \mathbf{0}$ but $(A - \lambda I)^{k-1} \mathbf{v} \neq \mathbf{0}$.

---

## 12.2 The Jordan Canonical Form Theorem

!!! theorem "Theorem 12.1 (Jordan Canonical Form Theorem)"
    Every complex square matrix $A$ is similar to a Jordan Canonical Form $J$, and $J$ is unique up to the permutation of its blocks:
    $$P^{-1} A P = J = \operatorname{diag}(J_{k_1}(\lambda_1), J_{k_2}(\lambda_2), \ldots, J_{k_m}(\lambda_m))$$
    - Each Jordan block corresponds to one linearly independent eigenvector (geometric multiplicity).
    - The sum of the orders of all Jordan blocks for a given eigenvalue equals its algebraic multiplicity.

---

## 12.3 The Minimal Polynomial

!!! definition "Definition 12.3 (Minimal Polynomial)"
    The **minimal polynomial** $m(\lambda)$ is the monic polynomial of lowest degree such that $m(A) = O$.
    **Properties**:
    1.  $m(\lambda)$ divides the characteristic polynomial $p(\lambda)$.
    2.  $A$ is diagonalizable $\iff$ $m(\lambda)$ has no repeated roots.
    3.  The multiplicity of $\lambda_i$ in $m(\lambda)$ is equal to the size of the largest Jordan block corresponding to that eigenvalue.

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

The Jordan Canonical Form is the final verdict on the structure of square matrices:

1.  **Completing the Defective**: It perfectly fills the gap in eigenvectors for non-diagonalizable matrices by introducing the "1" step structure (Jordan block).
2.  **Polynomial Depth**: The correspondence between the minimal polynomial and JCF block sizes reveals the geometric depth of a matrix as a root of a polynomial.
3.  **Structural Uniqueness**: JCF establishes the classification standard for matrix similarity classes, providing the most precise theoretical framework for analyzing matrix functions and power series convergence.
