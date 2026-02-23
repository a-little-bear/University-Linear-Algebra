# Chapter 06: Eigenvalues and Eigenvectors

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch03) · Linear Transformations (Ch05)

**Chapter Outline**: Definition of Eigenvalues and Eigenvectors → The Characteristic Equation and Polynomial → Algebraic and Geometric Multiplicities → Similarity Transformations → Diagonalization Criteria → Spectra of Special Matrices (Symmetric, Triangular) → Cayley-Hamilton Theorem → Matrix Powers and Stability

**Extension**: Eigenvalue analysis is key to understanding dynamical system stability, Google's PageRank algorithm, and quantum energy levels; it is the prerequisite for the Jordan Canonical Form (Ch12).

</div>

Eigenvalues and eigenvectors reveal the most essential "invariant directions" of a linear transformation. A complex matrix acting on a specific vector may manifest simply as a scaling operation. Finding these scaling factors and their corresponding directions is the core method for simplifying matrix arithmetic and analyzing long-term system behavior.

---

## 06.1 Definitions and Characteristic Equations

!!! definition "Definition 06.1 (Eigenvalues and Eigenvectors)"
    Let $A$ be an $n \times n$ matrix. If there exists a non-zero vector $\mathbf{v}$ and a scalar $\lambda$ such that:
    $$A\mathbf{v} = \lambda\mathbf{v}$$
    then $\lambda$ is an **eigenvalue** of $A$, and $\mathbf{v}$ is the corresponding **eigenvector**.

!!! definition "Definition 06.2 (Characteristic Equation)"
    The equation $\det(A - \lambda I) = 0$ is called the **characteristic equation** of $A$. The left side $p(\lambda) = \det(A - \lambda I)$ is a polynomial of degree $n$, known as the **characteristic polynomial**.

---

## 06.2 Multiplicity and Eigenspaces

!!! definition "Definition 06.3 (Algebraic and Geometric Multiplicity)"
    1.  **Algebraic Multiplicity**: The multiplicity of $\lambda_i$ as a root of the characteristic equation.
    2.  **Geometric Multiplicity**: The dimension of the eigenspace $E_{\lambda_i} = \ker(A - \lambda_i I)$, i.e., the maximum number of linearly independent eigenvectors.
    **Property**: Geometric Multiplicity $\leq$ Algebraic Multiplicity.

---

## 06.3 Diagonalization

!!! theorem "Theorem 06.1 (Diagonalization Theorem)"
    An $n \times n$ matrix $A$ is diagonalizable $\iff$ $A$ has $n$ linearly independent eigenvectors $\iff$ for every eigenvalue, its geometric multiplicity equals its algebraic multiplicity.
    The diagonal form is: $P^{-1}AP = \Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$.

---

## 06.4 Cayley-Hamilton Theorem

!!! theorem "Theorem 06.2 (Cayley-Hamilton Theorem)"
    Every square matrix satisfies its own characteristic equation. That is, if $p(\lambda) = \det(A - \lambda I)$, then $p(A) = O$.
    **Application**: This theorem can be used to efficiently calculate high powers of a matrix and find the inverse $A^{-1}$.

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

Eigen-analysis is the deepest deconstruction of a square matrix:

1.  **Coordinate Invariance**: Eigenvalues are intrinsic properties of a matrix, independent of the basis, making them ideal carriers for physical laws.
2.  **Deconstruction and Reconstruction**: The process of diagonalization is essentially finding a perfect basis that makes the action of a linear operator pure and independent.
3.  **Polynomial Constraints**: The Cayley-Hamilton theorem establishes the ultimate link between matrix arithmetic and polynomial algebra, revealing a certain self-recursive property of square matrices.
