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
       The standard problem is "Matrix $	o$ Spectrum," which is always solvable. IEP is "Spectrum $	o$ Matrix (with constraints)," which involves complex existence conditions and often lacks a unique solution.

2. **[Construction] Construct a $2 	imes 2$ symmetric matrix with eigenvalues 5 and 1.**
   ??? success "Solution"
       The simplest solution is $\operatorname{diag}(5, 1)$. A non-trivial one can be obtained via rotation: $\begin{pmatrix} 3 & 2 \ 2 & 3 \end{pmatrix}$.

3. **[Schur-Horn] Determine if there exists a symmetric matrix with spectrum $\{4, 2\}$ and diagonal $\{3, 3\}$.**
   ??? success "Solution"
       Yes. Check majorization: $3 \le 4$ and $3+3=4+2=6$. The condition holds. An example is $\begin{pmatrix} 3 & 1 \ 1 & 3 \end{pmatrix}$.

4. **[NIEP] Why is the set $\{3, 1, 1\}$ unlikely to be the spectrum of a nonnegative matrix without a large diagonal?**
   ??? success "Solution"
       By the Perron-Frobenius theorem, the Perron root must dominate the spectrum. While $\{3, 1, 1\}$ satisfies this, the trace condition $\operatorname{tr}(A) = 5$ and the requirement for non-diagonal entries to "shift" eigenvalues suggests a purely positive real spectrum is restricted for non-diagonal nonnegative matrices.

5. **[Jacobi] State the interlacing property required for the Jacobi IEP.**
   ??? success "Solution"
       Given eigenvalues $\lambda$ of $J_n$ and $\mu$ of its $(n-1)$-order principal submatrix, a Jacobi matrix exists if and only if the eigenvalues are strictly interlaced: $\lambda_1 < \mu_1 < \lambda_2 < \dots < \mu_{n-1} < \lambda_n$.

6. **[Calculation] Write the companion matrix for the spectrum $\{2, -1, -1\}$. Is it nonnegative?**
   ??? success "Solution"
       Characteristic polynomial $p(\lambda) = (\lambda-2)(\lambda+1)^2 = \lambda^3 - 3\lambda - 2$. Companion matrix: $\begin{pmatrix} 0 & 0 & 2 \ 1 & 0 & 3 \ 0 & 1 & 0 \end{pmatrix}$. It is nonnegative.

7. **[Weyl-Horn] What is the relationship between the product of eigenvalues and the product of singular values?**
   ??? success "Solution"
       According to the Weyl product inequalities, $\prod_{i=1}^k |\lambda_i| \le \prod_{i=1}^k s_i$ for all $k$, with equality at $k=n$.

8. **[Applications] Describe how IEP is used in mass-spring system design.**
   ??? success "Solution"
       Engineers measure the natural frequencies (eigenvalues) of a structure and use IEP to back-calculate the required stiffness coefficients or mass distributions to achieve desired vibration characteristics.

9. **[Numerical] Describe the "Lift-and-Project" approach for IEP.**
   ??? success "Solution"
       1. Lift: Find a matrix with the correct spectrum (e.g., $Q\Lambda Q^T$). 2. Project: Enforce the structural constraint (e.g., set specific entries to zero). 3. Iterate until the matrix satisfies both conditions.

10. **[Stochastic IEP] What are the necessary conditions for a set to be the spectrum of a stochastic matrix?**
    ??? success "Solution"
        1. The Perron root is 1. 2. All eigenvalues satisfy $|\lambda_i| \le 1$. 3. Non-real eigenvalues appear in conjugate pairs. 4. The trace conditions $\sum \lambda_i^k \ge 0$ hold.

## Chapter Summary

This chapter explores the inverse mapping from spectral data to matrix structure:

1. **Existence Criteria**: Detailed the majorization conditions for diagonal entries and the Perron-Frobenius requirements for nonnegativity.
2. **Structural IEPs**: Analyzed Jacobi and Toeplitz structures, highlighting the role of interlacing eigenvalues in reconstruction.
3. **Invariance Relations**: Established the link between eigenvalues and singular values via Weyl-Horn theory.
4. **Computational Frameworks**: Introduced numerical methods like homotopy continuation and projection for solving large-scale structured inverse problems.
