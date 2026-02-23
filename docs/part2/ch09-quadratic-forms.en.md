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

1. **[Fundamentals] Write the symmetric matrix representation for $f(x, y) = x^2 + 4xy + 3y^2$.**
   ??? success "Solution"
       $A = \begin{pmatrix} 1 & 2 \\ 2 & 3 \end{pmatrix}$. Note the cross-term coefficient 4 is split between $a_{12}$ and $a_{21}$.

2. **[Completion] Use the method of completing the square to find the canonical form of $Q(x, y) = x^2 + 4xy$.**
   ??? success "Solution"
       $Q = (x+2y)^2 - 4y^2$. Let $u = x+2y, v = y$, the form is $u^2 - 4v^2$.

3. **[Congruence] If $A$ and $B$ are congruent, do they have the same eigenvalues?**
   ??? success "Solution"
       Not necessarily. Congruence $P^T AP$ is different from similarity $P^{-1} AP$. Congruence preserves inertia, not eigenvalues.

4. **[Definiteness] Use principal minors to determine if $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ is positive definite.**
   ??? success "Solution"
       $D_1 = 2 > 0, D_2 = 4 - 1 = 3 > 0$. All leading principal minors are positive, so it is positive definite.

5. **[Inertia] Find the inertia and signature of $Q = x^2 - y^2$.**
   ??? success "Solution"
       Positive inertia $p=1$, negative inertia $q=1$. The signature is $p-q = 0$.

6. **[Geometry] What curves do $x^2 + y^2 = 1$ and $x^2 - y^2 = 1$ represent?**
   ??? success "Solution"
       The first represents an ellipse (circle), the second represents a hyperbola.

7. **[Rayleigh Quotient] What is the maximum value of $R(\mathbf{x}) = \frac{\mathbf{x}^T A \mathbf{x}}{\mathbf{x}^T \mathbf{x}}$?**
   ??? success "Solution"
       The largest eigenvalue $\lambda_{\max}$.

8. **[Skew-symmetric] Prove: for any skew-symmetric matrix $A$, $\mathbf{x}^T A \mathbf{x}$ is always 0.**
   ??? success "Solution"
       $\mathbf{x}^T A \mathbf{x} = (\mathbf{x}^T A \mathbf{x})^T = \mathbf{x}^T A^T \mathbf{x} = \mathbf{x}^T (-A) \mathbf{x} \implies 2\mathbf{x}^T A \mathbf{x} = 0$.

9. **[Normal Form] Transform $Q = 2x^2 + 2y^2$ into normal form.**
   ??? success "Solution"
       Let $u = \sqrt{2}x, v = \sqrt{2}y$, the normal form is $u^2 + v^2$.

10. **[Application] How is a quadratic form used to identify a minimum in optimization?**
    ??? success "Solution"
        If the Hessian matrix at a critical point is positive definite, the point is a local minimum.

## Chapter Summary

Quadratic forms bridge polynomial algebra and spatial geometry:

1.  **Structural Simplification**: Congruence transformations are the blades that eliminate cross-terms and simplify coordinate systems, revealing the core inertia of symmetric operators.
2.  **Sign Dynamics**: The Law of Inertia establishes the deepest topological invariants of a quadratic function, which remain unchanged regardless of coordinate rotation or scaling.
3.  **Core of Definiteness**: As the algebraic criterion for energy, distance, and stability analysis, the theory of definiteness is the vital junction connecting linear algebra with analysis and physics.
