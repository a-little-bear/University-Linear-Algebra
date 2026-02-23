# Chapter 09: Quadratic Forms and Bilinear Forms

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Orthogonality (Ch7) · Eigenvalues (Ch6)

**Chapter Outline**: Definition of Bilinear Forms → Symmetric Bilinear Forms → Definition of Quadratic Forms $x^T Ax$ → Canonical Forms (Completing the Square, Congruence) → Law of Inertia (Sylvester) → Normal Forms → Positivity Tests (Hurwitz/Sylvester Criterion) → Geometric Meaning (Classification of Quadrics)

**Extension**: Quadratic forms are the algebraic heart of extremum analysis for multivariable functions, kinetic energy in classical mechanics, and optimization theory.

</div>

The study of quadratic forms focuses on scalar functions that are quadratic polynomials of variables. Through linear substitution (congruence transformation), we can eliminate complex cross-terms to reveal the essential geometric features of the function.

---

## 09.1 Quadratic Forms and Congruence

!!! definition "Definition 09.1 (Quadratic Form)"
    A quadratic form in $n$ variables is a homogeneous polynomial $Q(x) = x^T A x$, where $A$ is a symmetric matrix.

!!! theorem "Theorem 09.3 (Sylvester's Law of Inertia)"
    For real symmetric matrices, the number of positive entries (positive inertia index) and negative entries (negative inertia index) in the normal form are invariants under non-singular linear substitutions.

---

## Exercises

1. **[Fundamentals] Write the symmetric matrix representation for the quadratic form $f(x, y) = x^2 + 4xy + 3y^2$.**
   ??? success "Solution"
       $A = \begin{pmatrix} 1 & 2 \\ 2 & 3 \end{pmatrix}$. Note the cross-term $4xy$ is split equally between $a_{12}$ and $a_{21}$.

2. **[Completing the Square] Transform $Q(x, y) = x^2 + 4xy$ into a canonical form by completing the square.**
   ??? success "Solution"
       $Q = (x+2y)^2 - 4y^2$. Let $u = x+2y, v = y$, then $Q = u^2 - 4v^2$.

3. **[Congruence] If $A$ and $B$ are congruent, must they have the same eigenvalues?**
   ??? success "Solution"
       No. Similarity transformations preserve eigenvalues, but congruence transformations $B = P^T A P$ only preserve the count of positive and negative inertia indices. For example, $\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$ is congruent to $\begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$, but they have different eigenvalues.

4. **[Positivity Test] Use principal minors to determine if $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ is positive definite.**
   ??? success "Solution"
       $D_1 = 2 > 0$.
       $D_2 = \det A = 4-1 = 3 > 0$.
       Since all leading principal minors are positive, $A$ is positive definite.

5. **[Inertia Index] Find the inertia index for $Q = x^2 - y^2$.**
   ??? success "Solution"
       Positive inertia index $p = 1$, negative inertia index $q = 1$. The signature is $p-q = 0$.

6. **[Geometric Classification] What shape does $x^2 + y^2 = 1$ represent in $\mathbb{R}^2$? How about $x^2 - y^2 = 1$?**
   ??? success "Solution"
       - $x^2 + y^2 = 1$ represents a circle (ellipse).
       - $x^2 - y^2 = 1$ represents a hyperbola.

7. **[Rayleigh Quotient] What is the maximum value of the Rayleigh Quotient $R(x) = \frac{x^T A x}{x^T x}$?**
   ??? success "Solution"
       The maximum value is the largest eigenvalue $\lambda_{\max}$ of $A$, attained at the corresponding eigenvector.

8. **[Skew-symmetric] Prove: For any skew-symmetric matrix $A$ ($A^T = -A$), the quadratic form $x^T A x$ is identically 0.**
   ??? success "Solution"
       Since $x^T A x$ is a scalar, $(x^T A x)^T = x^T A^T x = x^T (-A) x = -x^T A x$.
       A number equal to its own negative must be 0.

9. **[Normal Form] Transform $Q = 2x^2 + 2y^2$ into normal form (coefficients only $1, -1, 0$).**
   ??? success "Solution"
       Let $u = \sqrt{2}x, v = \sqrt{2}y$, then $Q = u^2 + v^2$. The normal form is $\sum y_i^2$.

10. **[Application] How are quadratic forms used to determine extremum points in optimization?**
    ??? success "Solution"
        By checking the Hessian matrix (second derivative matrix) of the objective function at a critical point. If the Hessian is a positive definite quadratic form, the point is a local minimum; if negative definite, a local maximum; if indefinite, a saddle point.

## Chapter Summary

Quadratic forms bridge polynomial algebra and spatial geometry:

1. **Structural Simplification**: Congruence transformations are the algebraic tools for eliminating cross-terms and simplifying coordinate systems.
2. **Signature Dynamics**: The Law of Inertia reveals the deepest topological properties of quadratic functions.
3. **Positivity**: As the core of energy, distance, and stability analysis, positive definiteness is the most applied part of quadratic form theory.
