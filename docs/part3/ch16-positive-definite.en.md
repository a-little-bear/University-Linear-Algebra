# Chapter 16: Positive Definite Matrices

<div class="context-flow" markdown>

**Prerequisites**: Quadratic Forms (Ch9) · Eigenvalues (Ch6) · Matrix Factorization (Ch10) · Inner Product (Ch8)

**Chapter Outline**: Definition of Positive Definiteness → Equivalent Characterizations (Eigenvalues, Pivots, Determinants, Cholesky) → Leading Principal Minors (Sylvester's Criterion) → Geometry of Positive Definite Matrices (Ellipsoids) → The Cone of Positive Definite Matrices → Positive Semi-definiteness → Sums and Inverses of PD Matrices → Generalized Eigenvalue Problems → Applications in Optimization

**Extension**: Positive definite matrices are the "positive numbers" of matrix algebra; they define energy surfaces and metrics in optimization and statistics.

</div>

Positive definite (PD) matrices are the most important class of symmetric matrices. They arise in every field of science, from the stiffness of structures in engineering to the covariance of random variables in statistics. A matrix is PD if its associated quadratic form $x^T Ax$ is strictly positive for all non-zero $x$. This property ensures that the matrix is invertible, has real positive eigenvalues, and defines a "bowl-shaped" energy surface. This chapter explores the multiple ways to identify PD matrices and their role in defining distances and solving optimization problems.

---

## 16.1 Characterizations of Positive Definiteness

!!! definition "Definition 16.1 (Positive Definite)"
    A symmetric matrix $A$ is **positive definite** (PD) if $x^T Ax > 0$ for all $x \in \mathbb{R}^n \setminus \{0\}$. If $x^T Ax \ge 0$, it is **positive semi-definite** (PSD).

!!! theorem "Theorem 16.1 (Equivalent Conditions for PD)"
    For a symmetric matrix $A$, the following are equivalent:
    1. $A$ is positive definite.
    2. All eigenvalues of $A$ are strictly positive ($\lambda_i > 0$).
    3. All leading principal minors are strictly positive (Sylvester's Criterion).
    4. All pivots in Gaussian elimination (without row swaps) are strictly positive.
    5. There exists a unique Cholesky factorization $A = LL^T$ with $l_{ii} > 0$.

---

## Exercises

1. **[Fundamentals] Is $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$ positive definite?**
   ??? success "Solution"
       Leading principal minors: $D_1 = 2 > 0$ and $D_2 = \det A = 3 > 0$. Yes, it is positive definite.

2. **[Eigenvalues] What is the smallest possible eigenvalue of a PD matrix?**
   ??? success "Solution"
       It must be strictly greater than 0. The minimum eigenvalue $\lambda_{\min}$ is the global minimum of the Rayleigh quotient $x^T Ax / \|x\|^2$.

3. **[Geometry] Describe the set $\{ x : x^T Ax = 1 \}$ for a $2 \times 2$ PD matrix.**
   ??? success "Solution"
       It is an ellipse centered at the origin. The semi-axes lengths are $1/\sqrt{\lambda_i}$, oriented along the eigenvectors of $A$.

4. **[Sums] Show that if $A$ and $B$ are PD, then $A+B$ is PD.**
   ??? success "Solution"
       $x^T(A+B)x = x^T Ax + x^T Bx$. Since $x^T Ax > 0$ and $x^T Bx > 0$ for non-zero $x$, their sum is also strictly positive. The set of PD matrices forms a **convex cone**.

5. **[Inversion] If $A$ is PD, is $A^{-1}$ also PD?**
   ??? success "Solution"
       Yes. The eigenvalues of $A^{-1}$ are $1/\lambda_i$. Since $\lambda_i > 0$, we have $1/\lambda_i > 0$. Thus $A^{-1}$ is PD.

6. **[Determinant] Prove that $\det A > 0$ for any PD matrix.**
   ??? success "Solution"
       Since $\det A = \prod \lambda_i$ and all $\lambda_i > 0$, the product must be positive. This is a necessary (but not sufficient) condition for PD.

7. **[Cholesky] Why is Cholesky factorization a good test for positive definiteness?**
   ??? success "Solution"
       It is computationally efficient ($O(n^3/3)$) and numerically stable. If the algorithm attempts to take the square root of a non-positive number, the matrix is not PD.

8. **[Inner Product] Show that $A$ is PD iff $\langle x, y \rangle_A = x^T Ay$ defines a valid inner product.**
   ??? success "Solution"
       Linearity and symmetry are inherited from matrix properties. Positivity $\langle x, x \rangle_A > 0$ for $x \neq 0$ is exactly the definition of positive definiteness. Every PD matrix defines a geometry.

9. **[Semi-definite] How do the leading principal minors of a positive semi-definite matrix behave?**
   ??? success "Solution"
       They are non-negative ($D_i \ge 0$). Note: unlike the PD case, $D_i \ge 0$ for all *leading* minors is not enough for PSD; *all* principal minors must be non-negative.

10. **[Rank] What is the rank of an $n \times n$ PD matrix?**
    ??? success "Solution"
        It is always full rank ($n$). Since all eigenvalues are non-zero, the kernel is $\{0\}$.

## Chapter Summary

This chapter explores the "positive scalars" of the matrix world:

1. **Analytical Unity**: Integrated eigenvalues, pivots, and minors into a unified theory of positivity.
2. **Computational Landmark**: Developed the Cholesky factorization as the definitive algorithm for PD systems.
3. **Geometric Ellipsoids**: Linked PD matrices to the curvature of space, providing the metrics used in high-dimensional optimization.
4. **Conical Structure**: Positioned the set of PD matrices as a convex cone, establishing the foundation for modern Semidefinite Programming.
