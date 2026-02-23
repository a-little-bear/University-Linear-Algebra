# Chapter 09: Quadratic Forms and Inertia

<div class="context-flow" markdown>

**Prerequisites**: Matrix Symmetry (Ch2) · Eigenvalues (Ch6) · Inner Product (Ch8) · Change of Basis (Ch4)

**Chapter Outline**: Definition of Quadratic Forms → Matrix Representation ($x^T A x$) → Congruence and Change of Variables → Sylvester's Law of Inertia → Definiteness (Positive, Negative, Indefinite) → Principal Axes Theorem (Diagonalization) → Quadratic Forms and Conic Sections → Signature of a Form

**Extension**: Quadratic forms describe the "energy surface" of a system; they are the second-order approximations used in optimization and general relativity.

</div>

A **quadratic form** is a scalar-valued function $q(x) = x^T A x$ where $A$ is a symmetric matrix. Unlike linear forms, quadratic forms capture the curvature and curvature-like behavior of functions. By choosing the right coordinate system (the "principal axes"), we can simplify a quadratic form into a weighted sum of squares. This chapter details the classification of forms into positive definite, negative definite, or indefinite categories and establishes **Sylvester's Law of Inertia**, which ensures that the fundamental "shape" of the form is invariant under non-singular coordinate changes.

---

## 09.1 Matrix Representation and Definiteness

!!! definition "Definition 09.1 (Quadratic Form)"
    A quadratic form on $\mathbb{R}^n$ is a function $q(\mathbf{x})$ that can be written as:
    $$q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x} = \sum_{i,j} a_{ij} x_i x_j$$
    where $A$ is a symmetric matrix.

!!! theorem "Theorem 09.1 (Sylvester's Law of Inertia)"
    The number of positive, negative, and zero eigenvalues of a symmetric matrix $A$ is invariant under congruence transformations ($P^T A P$ for non-singular $P$). This triplet $(p, n, z)$ is called the **inertia** of the matrix.

---

## Exercises

1. **[Fundamentals] Write the matrix representation of $q(x, y) = x^2 + 4xy + 3y^2$.**
   ??? success "Solution"
       The diagonal entries are the coefficients of squares, and the off-diagonals are half the coefficient of the cross-term. $A = \begin{pmatrix} 1 & 2 \\ 2 & 3 \end{pmatrix}$.

2. **[Principal Axes] Find the principal axes transformation for $q(x, y) = 5x^2 + 5y^2$.**
   ??? success "Solution"
       The matrix is $5I$, which is already diagonal. The principal axes are the standard $x$ and $y$ axes. The form represents a circular bowl.

3. **[Definiteness] Is $q(x, y) = x^2 + 2xy + y^2$ positive definite?**
   ??? success "Solution"
       $q(x, y) = (x+y)^2$. It is $\ge 0$ for all $(x, y)$, but it is 0 for non-zero vectors like $(1, -1)$. Thus it is **positive semi-definite** but not positive definite.

4. **[Eigenvalues] Relate the definiteness of $A$ to its eigenvalues.**
   ??? success "Solution"
       - **Positive Definite**: All $\lambda_i > 0$.
       - **Positive Semi-definite**: All $\lambda_i \ge 0$.
       - **Indefinite**: At least one $\lambda > 0$ and one $\lambda < 0$.

5. **[Congruence] Show that $A$ and $P^T A P$ represent the same quadratic form in different coordinates.**
   ??? success "Solution"
       Let $\mathbf{x} = P \mathbf{y}$. Then $\mathbf{x}^T A \mathbf{x} = (P\mathbf{y})^T A (P\mathbf{y}) = \mathbf{y}^T (P^T A P) \mathbf{y}$. Coordinate changes transform the matrix by congruence.

6. **[Sylvester] Find the inertia of $A = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$.**
   ??? success "Solution"
       The eigenvalues are 1 and -1. The inertia is $(1, 1, 0)$. The form is indefinite (a saddle shape).

7. **[Signature] Define the signature of a quadratic form.**
   ??? success "Solution"
       The difference between the number of positive and negative eigenvalues: $s = p - n$. Together with the rank, it uniquely identifies the form under congruence.

8. **[Calculation] Check if $\begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ is positive definite using the leading principal minors.**
   ??? success "Solution"
       Minors: $D_1 = 2 > 0$ and $D_2 = \det A = 3 > 0$. Since all leading principal minors are positive, the matrix is positive definite.

9. **[Geometry] What geometric shape is represented by $x^2 - y^2 = 1$?**
   ??? success "Solution"
       A hyperbola. In quadratic forms, indefinite matrices correspond to hyperbolic conic sections.

10. **[Energy] In physics, why is the kinetic energy $T = \frac{1}{2} \dot{q}^T M \dot{q}$ always a positive definite quadratic form?**
    ??? success "Solution"
        Because kinetic energy must be positive for any non-zero motion. This requires the mass matrix $M$ to be positive definite, ensuring stable and realistic physical behavior.

## Chapter Summary

This chapter explores the second-order structure of linear operators:

1. **Functional Mapping**: Formalized quadratic forms as scalar fields governed by symmetric matrices.
2. **Spectral Canonicalization**: Utilized the Principal Axes Theorem to decouple variables and reveal the core curvature of the form.
3. **Inertia Invariance**: Established Sylvester's Law as the fundamental preservation law for the "shape" of space under transformation.
4. **Classification Theory**: Developed definiteness criteria to characterize the stability and extremum properties of quadratic surfaces.
