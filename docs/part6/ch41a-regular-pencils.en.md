# Chapter 41A: Matrix Pencils and Regular Pencil Theory

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Jordan Form (Ch12) · Smith Normal Form (Ch13B) · Schur Decomposition (Ch10)

**Chapter Outline**: Definition of Matrix Pencils → Regular vs. Singular Pencils → Generalized Eigenvalues (Homogeneous Coordinates) → Strict Equivalence → Weierstrass Canonical Form → Deflating Subspaces → Generalized Schur Decomposition (QZ Decomposition) → QZ Algorithm → Definite Hermite Pencils → Polynomial Eigenvalue Problems (PEP) & Linearization

**Extension**: Regular pencil theory is the mathematical basis for generalized eigenvalue solvers in LAPACK (`xGGES`); the QZ algorithm is the core of MATLAB's `eig(A, B)`.

</div>

The standard eigenvalue problem $Ax = \lambda x$ studies the pencil $A - \lambda I$. Replacing $I$ with a general matrix $B$ leads to the **Generalized Eigenvalue Problem** $Ax = \lambda Bx$ and the **Matrix Pencil** $A - \lambda B$. If $\det(A - \lambda B)$ is not identically zero, the pencil is **regular** and possesses a Weierstrass Canonical Form—a structure that separates finite eigenvalues from infinite ones.

---

## 41A.1 Core Concepts

!!! definition "Definition 41A.2 (Regular vs. Singular)"
    A square pencil $A - \lambda B$ is **regular** if $\det(A - \lambda B) 
ot\equiv 0$. Otherwise, it is **singular**.

!!! theorem "Theorem 41A.3 (Weierstrass Canonical Form)"
    For a regular pencil $A - \lambda B$, there exist non-singular $P, Q$ such that
    $$P(A - \lambda B)Q = \begin{pmatrix} J - \lambda I & 0 \ 0 & I - \lambda N \end{pmatrix},$$
    where $J$ is in Jordan form (finite eigenvalues) and $N$ is nilpotent (infinite eigenvalues).

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

This chapter extends eigenvalue theory to matrix pairs:

1. **Regularity Framework**: Established the criteria for well-posed generalized eigenvalue problems.
2. **Canonical Decomposition**: Developed the Weierstrass form to classify pencils under strict equivalence.
3. **Unitary Algorithms**: Introduced the QZ decomposition and algorithm as the numerically stable path for spectral computation.
4. **Polynomial Mapping**: Formulated linearizations to bridge high-degree matrix polynomials and linear pencils.
