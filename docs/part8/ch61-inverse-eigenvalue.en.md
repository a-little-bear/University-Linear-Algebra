# Chapter 61: Inverse Eigenvalue Problems

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Positive Definite Matrices (Ch16) · Nonnegative Matrices (Ch17) · Jordan Form (Ch12)

**Chapter Outline**: Framework of IEP → Symmetric IEP → Nonnegative IEP (NIEP) → Jacobi IEP → Toeplitz IEP → Stochastic IEP → Numerical Methods → Open Problems

**Extension**: Inverse eigenvalue problems are critical in structural engineering (mass-spring model refinement), control theory (pole placement), and molecular spectroscopy.

</div>

Inverse Eigenvalue Problems (IEP) involve constructing a matrix that satisfies certain structural constraints while possessing a prescribed spectrum. Unlike the direct problem (finding eigenvalues of a given matrix), IEPs deal with existence, uniqueness, and construction from the "effect" back to the "cause."

---

## 61.1 General Framework and Symmetry

!!! definition "Definition 61.1 (Inverse Eigenvalue Problem)"
    Given a set of scalars $\Lambda = \{\lambda_1, \dots, \lambda_n\}$ and a class of matrices $\mathcal{S}$, find $A \in \mathcal{S}$ such that $\sigma(A) = \Lambda$.

!!! theorem "Theorem 61.2 (Schur-Horn Theorem)"
    There exists a real symmetric matrix with eigenvalues $\lambda$ and diagonal entries $d$ if and only if $d$ is **majorized** by $\lambda$ ($d \prec \lambda$).

---

## Exercises

1. **[Concept] Contrast the logical direction of IEP with the standard eigenvalue problem.**
   ??? success "Solution"
       The standard problem is "Matrix $\to$ Spectrum," which is always solvable and the solution is unique. IEP is "Spectrum $\to$ Matrix (with constraints)," which involves complex existence conditions and often lacks a unique solution, acting as a form of "reverse engineering."

2. **[Construction] Construct a $2 \times 2$ symmetric matrix with eigenvalues 5 and 1.**
   ??? success "Solution"
       The simplest solution is the diagonal matrix $\operatorname{diag}(5, 1)$. A non-trivial one can be obtained via rotation: $\begin{pmatrix} 3 & 2 \\ 2 & 3 \end{pmatrix}$, which has characteristic equation $\lambda^2 - 6\lambda + 5 = 0$.

3. **[Schur-Horn] Determine if there exists a symmetric matrix with spectrum $\{4, 2\}$ and diagonal $\{3, 3\}$.**
   ??? success "Solution"
       Check majorization: $3 \le 4$ and $3+3=4+2=6$. The conditions are satisfied. Thus, such a matrix exists. An example is $\begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$.

4. **[NIEP] Why is the set $\{3, 1, 1\}$ unlikely to be the spectrum of a non-diagonal nonnegative matrix?**
   ??? success "Solution"
       By the Perron-Frobenius theorem, the Perron root must dominate the spectrum. While 3 is the Perron root here, the trace $\operatorname{tr}(A) = 5$ would be large. Usually, non-diagonal nonnegative matrices require some eigenvalues to have negative real parts or complex phases to "balance" the positive trace required for non-negativity.

5. **[Jacobi] State the interlacing property required for the Jacobi IEP.**
   ??? success "Solution"
       Given eigenvalues $\{\lambda_i\}$ of $J_n$ and $\{\mu_j\}$ of its $(n-1)$-order principal submatrix, a Jacobi matrix exists if and only if the eigenvalues are strictly interlaced: $\lambda_1 < \mu_1 < \lambda_2 < \mu_2 < \dots < \mu_{n-1} < \lambda_n$.

6. **[Calculation] Write the companion matrix for the spectrum $\{2, -1, -1\}$. Is it nonnegative?**
   ??? success "Solution"
       Characteristic polynomial $p(\lambda) = (\lambda-2)(\lambda+1)^2 = \lambda^3 - 3\lambda - 2$. The companion matrix is $\begin{pmatrix} 0 & 0 & 2 \\ 1 & 0 & 3 \\ 0 & 1 & 0 \end{pmatrix}$. All entries are nonnegative, so it is a valid NIEP solution.

7. **[Weyl-Horn] What is the relationship between the product of eigenvalues and the product of singular values?**
   ??? success "Solution"
       According to the Weyl product inequalities, $\prod_{i=1}^k |\lambda_i| \le \prod_{i=1}^k \sigma_i$ for all $k = 1, \dots, n$, with equality at $k=n$ (matching the absolute value of the determinant).

8. **[Applications] Describe how IEP is used in mass-spring system design.**
   ??? success "Solution"
       Engineers measure the natural frequencies (eigenvalues) of a structure and use IEP algorithms to back-calculate the required stiffness coefficients or mass distributions (entries of the stiffness/mass matrices) to achieve those frequencies.

9. **[Numerical] Describe the "Lift-and-Project" approach for solving structural IEPs.**
   ??? success "Solution"
       1. **Lift**: Find a matrix with the correct spectrum (e.g., $Q\Lambda Q^T$). 2. **Project**: Force the matrix into the constraint set $\mathcal{S}$ (e.g., zeroing specific entries). 3. **Iterate**: Repeat these steps until the matrix satisfies both spectral and structural constraints.

10. **[Stochastic IEP] What are the necessary conditions for a set to be the spectrum of a stochastic matrix?**
    ??? success "Solution"
        The spectrum must satisfy: (1) The Perron root is 1. (2) All eigenvalues satisfy $|\lambda_i| \le 1$. (3) Non-real eigenvalues appear in conjugate pairs. (4) The trace conditions $\sum \lambda_i^k \ge 0$ must hold for all $k$.

## Chapter Summary

This chapter explores the inverse mapping from spectral data to matrix structure:

1. **Existence Criteria**: Detailed the majorization conditions for diagonal entries and the Perron-Frobenius requirements for nonnegativity.
2. **Structural IEPs**: Analyzed Jacobi and Toeplitz structures, highlighting the role of interlacing eigenvalues in parameters reconstruction.
3. **Invariance Relations**: Established the link between eigenvalues and singular values via Weyl-Horn theory.
4. **Computational Frameworks**: Introduced numerical methods like homotopy continuation and projection for solving large-scale structured inverse problems.
