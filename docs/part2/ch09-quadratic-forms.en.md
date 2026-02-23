# Chapter 09: Quadratic and Bilinear Forms

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Eigenvalues and Diagonalization (Ch06) · Vector Spaces (Ch04)

**Chapter Outline**: Definition of Bilinear Forms → Symmetric Bilinear Forms → Definition of Quadratic Forms $\mathbf{x}^T A \mathbf{x}$ → Canonical and Normal Forms → Completing the Square and Elementary Transformations → Congruence Transformations → Sylvester's Law of Inertia (Inertia Indices) → Definiteness Criteria (Hurwitz Criterion) → Geometric Interpretation: Classification and Rotation of Quadric Surfaces → Applications: Extrema Analysis of Multivariate Functions (Hessian Matrices)

**Extension**: Quadratic forms are the algebraic characterization of local shapes of multivariate functions; they are the algebraic core of kinetic energy expressions in classical mechanics and optimization theory; they provide the foundation for studying manifold curvature and gravitational field equations.

</div>

Quadratic forms study scalar functions that take the form of homogeneous polynomials of degree two. By linear substitution (congruence transformations), we can eliminate complex cross-product terms, thereby revealing the essential geometric features of the function. This chapter establishes standard methods for determining the "sign" of a quadratic form and maps algebraic properties perfectly to the morphology of spatial surfaces.

---

## 9.1 Bilinear and Quadratic Forms

!!! definition "Definition 09.1 (Bilinear and Quadratic Forms)"
    1.  **Bilinear Form**: A mapping $B: V \times V \to F$ that is linear in each argument separately.
    2.  **Quadratic Form**: A function $Q(\mathbf{x}) = B(\mathbf{x}, \mathbf{x})$ derived from a symmetric bilinear form.
    In a given basis, it can be written as:
    $$Q(\mathbf{x}) = \sum_{i=1}^n \sum_{j=1}^n a_{ij} x_i x_j = \mathbf{x}^T A \mathbf{x}$$
    where $A$ can always be chosen as a **symmetric matrix**.

---

## 9.2 Congruence and Canonical Forms

!!! definition "Definition 09.2 (Congruence)"
    Two matrices $A$ and $B$ are **congruent** if there exists a non-singular matrix $P$ such that $B = P^T A P$.
    **Geometric Insight**: Congruence corresponds to a linear change of coordinates and does not change the set of values attained by the quadratic form.

!!! theorem "Theorem 09.1 (Sylvester's Law of Inertia)"
    When a quadratic form over $\mathbb{R}$ is reduced to canonical form $\sum d_i y_i^2$ via congruence, the number of positive coefficients $p$ (positive inertia index) and negative coefficients $q$ (negative inertia index) are invariants.
    - **Rank**: $r = p + q$.
    - **Signature**: $s = p - q$.

---

## 9.3 Criteria for Definiteness

!!! definition "Definition 09.3 (Definiteness Classification)"
    1.  **Positive Definite**: $Q(\mathbf{x}) > 0$ for all $\mathbf{x} \neq \mathbf{0}$. Eigenvalues are all positive.
    2.  **Positive Semi-definite**: $Q(\mathbf{x}) \ge 0$ for all $\mathbf{x}$. Eigenvalues are non-negative.
    3.  **Indefinite**: Takes both positive and negative values. Eigenvalues have mixed signs.

!!! theorem "Theorem 09.2 (Hurwitz Criterion / Leading Principal Minors)"
    A real symmetric matrix $A$ is positive definite iff all its leading principal minors are strictly positive.

---

## Exercises

**1. [Basics] Write the symmetric matrix representation for $f(x, y) = x^2 + 4xy + 3y^2$.**

??? success "Solution"
    **Steps:**
    1. Place squared coefficients on the diagonal: $a_{11}=1, a_{22}=3$.
    2. Split the cross-term coefficient 4 equally: $a_{12}=2, a_{21}=2$.
    **Matrix Representation**:
    $A = \begin{pmatrix} 1 & 2 \\ 2 & 3 \end{pmatrix}$.
    Verification: $(x \ y) \begin{pmatrix} 1 & 2 \\ 2 & 3 \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix} = x(x+2y) + y(2x+3y) = x^2 + 4xy + 3y^2$.

**2. [Completion] Use the method of completing the square to find the canonical form of $Q(x, y) = x^2 + 4xy$.**

??? success "Solution"
    **Steps:**
    1. Complete the square: $Q = (x^2 + 4xy + 4y^2) - 4y^2$.
    2. Write as a full square: $Q = (x+2y)^2 - 4y^2$.
    3. Let $u = x+2y, v = y$.
    **Conclusion**: The canonical form is $u^2 - 4v^2$. Its inertia indices are $p=1, q=1$.

**3. [Congruence] If $A$ and $B$ are congruent, must they have the same eigenvalues?**

??? success "Solution"
    **Conclusion**: Not necessarily.
    **Analysis**:
    - **Similarity** $P^{-1}AP$ preserves eigenvalues.
    - **Congruence** $P^T AP$ preserves **inertia indices** (the number of positive/negative eigenvalues) but not the specific values.
    - For example, $\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$ and $\begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$ are congruent (let $P=\sqrt{2}I$), but have different eigenvalues.

**4. [Definiteness] Use principal minors to determine if $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ is positive definite.**

??? success "Solution"
    **Calculation:**
    1. 1st minor: $D_1 = 2 > 0$.
    2. 2nd minor: $D_2 = 2 \cdot 2 - 1 \cdot 1 = 3 > 0$.
    **Conclusion**: Since all leading principal minors are positive, the matrix is positive definite.

**5. [Inertia] Find the inertia and signature of $Q = x^2 - y^2$.**

??? success "Solution"
    **Analysis:**
    1. The form is already canonical: $1x^2 + (-1)y^2$.
    2. Number of positive coefficients $p=1$.
    3. Number of negative coefficients $q=1$.
    4. Signature $s = p - q = 0$.
    **Conclusion**: The inertia is $(1, 1)$, and the signature is 0.

**6. [Geometry] What curves do $x^2 + y^2 = 1$ and $x^2 - y^2 = 1$ represent?**

??? success "Solution"
    **Classification:**
    1. $x^2 + y^2 = 1$: Represents a circle (or ellipse). The corresponding quadratic form is **positive definite**, resulting in closed level sets.
    2. $x^2 - y^2 = 1$: Represents a hyperbola. The corresponding quadratic form is **indefinite**, resulting in open level sets.
    This shows how the sign character of a quadratic form dictates the topology of its geometric level sets.

**7. [Rayleigh] What is the maximum value of the Rayleigh Quotient $R(\mathbf{x}) = \frac{\mathbf{x}^T A \mathbf{x}}{\mathbf{x}^T \mathbf{x}}$?**

??? success "Solution"
    **Theorem:**
    The maximum value of the Rayleigh Quotient is the **largest eigenvalue** $\lambda_{\max}$ of $A$, and the minimum is $\lambda_{\min}$.
    This builds a bridge between geometric "stretching" and the algebraic spectrum.

**8. [Skew-symmetry] Prove that for any skew-symmetric matrix $A$, the induced quadratic form $\mathbf{x}^T A \mathbf{x}$ is identically 0.**

??? success "Solution"
    **Proof:**
    1. A quadratic form is a scalar, so it equals its transpose: $\mathbf{x}^T A \mathbf{x} = (\mathbf{x}^T A \mathbf{x})^T$.
    2. Expand: $= \mathbf{x}^T A^T \mathbf{x}$.
    3. Use $A^T = -A$: $= \mathbf{x}^T (-A) \mathbf{x} = -\mathbf{x}^T A \mathbf{x}$.
    4. A number equal to its negative must be 0.
    **Conclusion**: This is why we only consider symmetric matrices when studying quadratic forms.

**9. [Normal Form] Transform $Q = 2x^2 + 2y^2$ into normal form.**

??? success "Solution"
    **Steps:**
    1. It's already in canonical form.
    2. To get the normal form (coefficients in $\{1, -1, 0\}$), scale the variables.
    3. Let $u = \sqrt{2}x, v = \sqrt{2}y$.
    **Conclusion**: The normal form is $u^2 + v^2$.

**10. [Application] How is a quadratic form used to determine the nature of critical points in optimization?**

??? success "Solution"
    **Analysis:**
    Using the **Hessian matrix** $H$ of a multivariate function:
    1. If $H$ at a critical point is **positive definite**, the point is a local minimum.
    2. If $H$ is **negative definite**, it is a local maximum.
    3. If $H$ is **indefinite**, it is a saddle point.
    This illustrates the core role of quadratic forms in analyzing the curvature of continuous space.

## Chapter Summary

Quadratic forms bridge polynomial algebra and spatial geometry:

1.  **Structural Simplification**: Congruence transformations are the blades that eliminate cross-terms and simplify coordinate systems, revealing the core inertia of symmetric operators.
2.  **Sign Dynamics**: The Law of Inertia establishes the deepest topological invariants of a quadratic function, which remain unchanged regardless of coordinate rotation or scaling.
3.  **Core of Definiteness**: As the algebraic criterion for energy, distance, and stability analysis, the theory of definiteness is the vital junction connecting linear algebra with analysis and physics.
