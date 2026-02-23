# Chapter 09: Quadratic and Bilinear Forms

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Eigenvalues and Diagonalization (Ch06)

**Chapter Outline**: Definition of Bilinear Forms → Symmetric Bilinear Forms → Definition of Quadratic Forms $\mathbf{x}^T A \mathbf{x}$ → Canonical and Normal Forms → Congruence Transformations → Sylvester's Law of Inertia → Definiteness Criteria (Hurwitz Criterion) → Geometric Interpretation (Quadric Surfaces) → Applications in Optimization (Hessian Matrices)

**Extension**: Quadratic forms are the algebraic core of multi-variable extrema analysis, kinetic energy expressions in classical mechanics, and optimization theory; they are the starting point for studying Symplectic Geometry and Complex Manifolds.

</div>

Quadratic forms study scalar functions that take the form of quadratic homogeneous polynomials. By linear substitution (congruence transformations), we can eliminate complex cross-products to reveal the essential geometric features of the function. This chapter establishes standard methods for determining the "sign" of a quadratic form.

---

## 09.1 Bilinear and Quadratic Forms

!!! definition "Definition 09.1 (Bilinear and Quadratic Forms)"
    1.  **Bilinear Form**: A mapping $B: V \times V \to F$ that is linear in each argument.
    2.  **Quadratic Form**: A function $Q(\mathbf{x}) = B(\mathbf{x}, \mathbf{x})$ derived from a symmetric bilinear form. In coordinates, it is $Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}$, where $A$ is a symmetric matrix.

---

## 09.2 Congruence and the Law of Inertia

!!! definition "Definition 09.2 (Congruence Transformation)"
    Two matrices $A$ and $B$ are **congruent** if there exists a non-singular matrix $P$ such that $B = P^T A P$. Congruence preserves symmetry.

!!! theorem "Theorem 09.1 (Sylvester's Law of Inertia)"
    Over the field of real numbers, when a quadratic form is reduced to canonical form via congruence, the number of positive coefficients ($p$) and negative coefficients ($q$) are invariants.

---

## 09.3 Criteria for Definiteness

!!! definition "Definition 09.3 (Classification of Definiteness)"
    1.  **Positive Definite**: $Q(\mathbf{x}) > 0$ for all $\mathbf{x} \neq \mathbf{0}$.
    2.  **Positive Semi-definite**: $Q(\mathbf{x}) \ge 0$.
    3.  **Indefinite**: Takes both positive and negative values.

!!! theorem "Theorem 09.2 (Hurwitz Criterion / Leading Principal Minors)"
    A real symmetric matrix $A$ is positive definite if and only if all its leading principal minors are positive.

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

Quadratic forms bridge polynomial algebra and spatial geometry:

1.  **Structural Simplification**: Congruence transformations are the blades that eliminate cross-terms and simplify coordinate systems, revealing the core inertia of symmetric operators.
2.  **Sign Dynamics**: The Law of Inertia establishes the deepest topological invariants of a quadratic function, which remain unchanged regardless of coordinate rotation or scaling.
3.  **Core of Definiteness**: As the algebraic criterion for energy, distance, and stability analysis, the theory of definiteness is the vital junction connecting linear algebra with analysis and physics.
