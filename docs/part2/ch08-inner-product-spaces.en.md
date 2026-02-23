# Chapter 08: Inner Product Spaces

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Orthogonality Basics (Ch07)

**Chapter Outline**: Definition of Abstract Inner Products → Axioms → Cauchy-Schwarz Inequality → Generalization of Norms and Distance → Adjoint Operators $T^*$ → Normal, Self-adjoint (Hermitian), and Unitary Operators → The Spectral Theorem (Spectral Decomposition) → Positive Operators → Introduction to Completeness and Hilbert Spaces

**Extension**: Inner product spaces introduce continuity and geometry to linear spaces, serving as a finite-dimensional microcosm of functional analysis; the spectral theorem is the mathematical core of "observables" in quantum mechanics.

</div>

Inner product spaces are a reinforcement of vector spaces. By defining an abstract inner product, we can talk not only about linear combinations but also about lengths, angles, and approximations. It is the mathematical foundation for modern signal processing, quantum mechanics, and optimization theory.

---

## 08.1 Definition and Core Properties

!!! definition "Definition 08.1 (Inner Product)"
    A vector space $V$ over a field $F$ is equipped with an inner product $\langle \cdot, \cdot \rangle$ if it satisfies:
    1.  **Positivity**: $\langle v, v \rangle \ge 0$, and equals 0 if and only if $v=0$.
    2.  **Conjugate Symmetry**: $\langle u, v \rangle = \overline{\langle v, u \rangle}$.
    3.  **Additivity and Homogeneity (Linearity in first slot)**: $\langle au+bv, w \rangle = a\langle u, w \rangle + b\langle v, w \rangle$.

!!! theorem "Theorem 08.1 (Cauchy-Schwarz Inequality)"
    For any vectors $u, v$ in an inner product space:
    $$|\langle u, v \rangle| \le \|u\| \|v\|$$
    Equality holds if and only if $u$ and $v$ are linearly dependent.

---

## 08.2 Adjoint Operators

!!! definition "Definition 08.2 (Adjoint Operator $T^*$)"
    Let $T$ be a linear operator on $V$. If there exists an operator $T^*$ such that for all $u, v \in V$:
    $$\langle Tu, v \rangle = \langle u, T^* v \rangle$$
    then $T^*$ is called the **adjoint** of $T$. In an orthonormal basis, its matrix representation is the conjugate transpose $A^*$.

---

## 08.3 Operator Classification and the Spectral Theorem

!!! definition "Definition 08.3 (Special Operator Classes)"
    1.  **Self-adjoint (Hermitian)**: $T^* = T$. Its eigenvalues must be real.
    2.  **Unitary**: $T^* T = TT^* = I$. Preserves inner products and norms.
    3.  **Normal**: $T^* T = TT^*$.

!!! theorem "Theorem 08.2 (The Spectral Theorem)"
    An operator $T$ is diagonalizable in an orthonormal basis $\iff$ $T$ is normal.
    In particular, if $T$ is self-adjoint, all its eigenvalues are real, and there exists an orthonormal basis of eigenvectors.

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

Inner product spaces establish the metrics of linear algebra:

1.  **Geometrization**: The introduction of norms, angles, and projections gives abstract spaces physical intuition and approximation power.
2.  **Operator Symmetry**: The classification of self-adjoint, normal, and unitary operators forms the core of spectral theory, revealing the most stable structures of linear operators.
3.  **Approximation Essence**: Orthogonal projection is the mathematical endpoint for solving all linear approximation and error-minimization problems.
