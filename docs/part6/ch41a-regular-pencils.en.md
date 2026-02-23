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

1. **[Regularity] Determine if $A - \lambda B$ is regular for $A = \begin{pmatrix} 1 & 0 \ 0 & 0 \end{pmatrix}, B = \begin{pmatrix} 0 & 0 \ 0 & 1 \end{pmatrix}$. Find eigenvalues.**

   ??? success "Solution"
       $\det(A - \lambda B) = \det \begin{pmatrix} 1 & 0 \ 0 & -\lambda \end{pmatrix} = -\lambda$. Since it's not identically zero, it's regular. Finite eigenvalue $\lambda=0$; one infinite eigenvalue.

2. **[Homogeneous] Why are homogeneous coordinates $(\alpha, \beta)$ used for generalized eigenvalues?**

   ??? success "Solution"
       They provide a unified treatment for finite ($\lambda = \alpha/\beta$) and infinite ($\beta=0$) eigenvalues, mapping the spectrum to the projective line $\mathbb{P}^1(\mathbb{C})$.

3. **[Equivalence] Define strict equivalence of pencils.**

   ??? success "Solution"
       $(A_1, B_1) \sim_s (A_2, B_2)$ if there exist non-singular $P, Q$ such that $PA_1Q = A_2$ and $PB_1Q = B_2$.

4. **[QZ Decomposition] What is the generalized Schur decomposition (QZ decomposition)?**

   ??? success "Solution"
       For any $A, B$, there exist unitary $Q, Z$ such that $Q^*AZ$ and $Q^*BZ$ are both upper triangular. The ratios of diagonal entries give the generalized eigenvalues.

5. **[Deflating Subspaces] Distinguish between invariant and deflating subspaces.**

   ??? success "Solution"
       An invariant subspace satisfies $A\mathcal{V} \subseteq \mathcal{V}$. A deflating subspace for $(A, B)$ satisfies $\dim(A\mathcal{V} + B\mathcal{V}) \le \dim \mathcal{V}$.

6. **[Hermite Pencils] When are all generalized eigenvalues of a Hermite pencil $A - \lambda B$ real?**

   ??? success "Solution"
       This is guaranteed if the pencil is **definite**, meaning there exists a linear combination $\alpha A + \beta B$ that is positive definite.

7. **[Condition Number] Define the condition number of a simple generalized eigenvalue $\lambda_0$.**

   ??? success "Solution"
       $\kappa(\lambda_0) = \frac{\|y\| \|x\|}{|y^* B x|}$, where $x, y$ are right and left eigenvectors. Sensitivity increases as $B$ becomes singular or the pencil approaches singularity.

8. **[Linearization] How can a quadratic eigenvalue problem $(\lambda^2 M + \lambda C + K)x = 0$ be solved?**

   ??? success "Solution"
       By transforming it into a linear pencil of twice the size (linearization), typically using the companion form: $\lambda \begin{pmatrix} M & 0 \ 0 & I \end{pmatrix} - \begin{pmatrix} -C & -K \ I & 0 \end{pmatrix}$.

9. **[Symmetry Preservation] What is a structure-preserving linearization?**

   ??? success "Solution"
       A linearization that preserves the symmetry or Hermite properties of the original matrix polynomial (e.g., if $M, C, K$ are symmetric).

10. **[Weierstrass to Jordan] How does the Weierstrass form reduce to the Jordan form when $B = I$?**

   ??? success "Solution"
        If $B=I$, then $s=0$ (no infinite part), and the form is simply $P(A - \lambda I)Q = J - \lambda I$, which is the standard Jordan form problem.

## Chapter Summary

This chapter extends eigenvalue theory to matrix pairs:

1. **Regularity Framework**: Established the criteria for well-posed generalized eigenvalue problems.
2. **Canonical Decomposition**: Developed the Weierstrass form to classify pencils under strict equivalence.
3. **Unitary Algorithms**: Introduced the QZ decomposition and algorithm as the numerically stable path for spectral computation.
4. **Polynomial Mapping**: Formulated linearizations to bridge high-degree matrix polynomials and linear pencils.
