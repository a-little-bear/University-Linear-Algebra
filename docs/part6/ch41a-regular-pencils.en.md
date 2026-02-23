# Chapter 41A: Regular Matrix Pencils

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · $\lambda$-matrices (Ch13B) · Generalized Eigenvalue Problems

**Chapter Outline**: From Single Matrices to Matrix Pairs → Definition of Matrix Pencils $A - \lambda B$ → Regular vs. Singular Pencils → Generalized Eigenvalues and Eigenvectors → The Generalized Characteristic Equation $\det(A - \lambda B) = 0$ → Weierstrass Canonical Form (Regular Case) → Handling Infinite Eigenvalues → Applications: Structural Dynamics and State-Space Analysis of Generalized Linear Systems

**Extension**: Matrix pencil theory is the mathematical foundation for studying generalized linear dynamical systems (such as Differential-Algebraic Equations, DAEs); it generalizes spectral theory from a single operator $A$ to the relative evolution of two operators $A$ and $B$, which is key to stability analysis in complex systems.

</div>

In classical eigenvalue problems, we study $A\mathbf{x} = \lambda \mathbf{x}$. However, in many engineering problems (such as finite element analysis), equations take the form $A\mathbf{x} = \lambda B\mathbf{x}$. **Matrix Pencils** $A - \lambda B$ are the tools for describing such relative spectral relationships. When the characteristic equation is not identically zero, the pencil is termed **Regular**. This chapter introduces the decomposition theory of regular matrix pencils and their dominant role in dynamical systems.

---

## 41A.1 Basic Concepts of Matrix Pencils

!!! definition "Definition 41A.1 (Matrix Pencil)"
    Given two $m \times n$ matrices $A$ and $B$, the set $\{ A - \lambda B : \lambda \in \mathbb{C} \}$ is called a **Matrix Pencil**.

!!! definition "Definition 41A.2 (Regular Pencil)"
    If $A$ and $B$ are square matrices of the same order and the characteristic polynomial $p(\lambda) = \det(A - \lambda B)$ is not identically zero, the pencil is **Regular**. Otherwise, it is **Singular**.

---

## 41A.2 Generalized Eigenvalues and Canonical Forms

!!! definition "Definition 41A.3 (Generalized Eigenvalues)"
    Scalars $\lambda$ satisfying $\det(A - \lambda B) = 0$ are **finite generalized eigenvalues**.
    If $\det(B) = 0$, the pencil may also possess **infinite eigenvalues** $\lambda = \infty$.

!!! theorem "Theorem 41A.1 (Weierstrass Canonical Form)"
    For a regular pencil $A - \lambda B$, there exist non-singular matrices $P$ and $Q$ such that:
    $$P(A - \lambda B)Q = \operatorname{diag}(J - \lambda I, I - \lambda N)$$
    - $J$ is in Jordan form, corresponding to finite eigenvalues.
    - $N$ is a nilpotent Jordan matrix, corresponding to infinite eigenvalues.

---

## 41A.3 Application in Dynamics

!!! technique "Application: Vibration Analysis"
    In mechanical vibrations, the governing equation is $M \ddot{x} + K x = 0$. Assuming a solution $x = e^{i\omega t} v$ leads to the generalized eigenvalue problem $Kv = \omega^2 Mv$. Here $K$ is the stiffness matrix and $M$ is the mass matrix.

---

## Exercises

**1. [Calculation] Find the characteristic polynomial of the pencil $A - \lambda B$ where $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}, B = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. Write the difference matrix: $A - \lambda B = \begin{pmatrix} 1-\lambda & 0 \\ 0 & 2 \end{pmatrix}$.
    2. Compute the determinant: $\det(A - \lambda B) = (1-\lambda) \cdot 2 = 2 - 2\lambda$.
    **Conclusion**: The characteristic polynomial is $2 - 2\lambda$.

**2. [Eigenvalues] Find the generalized eigenvalues (including $\infty$) for the previous problem.**

??? success "Solution"
    **Analysis:**
    1. Finite eigenvalues: Set $2 - 2\lambda = 0 \implies \lambda = 1$.
    2. Infinite eigenvalues: Since $\det(B) = 0$, an infinite eigenvalue exists.
    3. Degree check: Since $n=2$ but the polynomial degree is 1, the missing root is $\infty$.
    **Conclusion**: The generalized eigenvalues are $\{1, \infty\}$.

**3. [Regularity] Determine if $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, B = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ is a regular pencil.**

??? success "Solution"
    **Calculation:**
    $A - \lambda B = \begin{pmatrix} 1 & -\lambda \\ 0 & 0 \end{pmatrix}$.
    $\det(A - \lambda B) = 1 \cdot 0 - (-\lambda) \cdot 0 = 0$.
    **Conclusion**: Since the determinant is identically zero, the pencil is **Singular**.

**4. [Relation] If $B$ is invertible, how are generalized eigenvalues related to standard ones?**

??? success "Solution"
    **Conclusion:**
    If $B$ is invertible, $A\mathbf{x} = \lambda B\mathbf{x}$ is equivalent to **$(B^{-1}A)\mathbf{x} = \lambda \mathbf{x}$**.
    In this case, all eigenvalues are finite and correspond to the standard eigenvalues of $B^{-1}A$.

**5. [Infinite] What is the physical meaning of an infinite eigenvalue in a matrix pencil?**

??? success "Solution"
    **Explanation:**
    In Differential-Algebraic Equations (DAEs), infinite eigenvalues typically correspond to **algebraic constraints**. They represent variables with infinitely fast response times (instantaneous response) or indicate a degeneracy in the system's dynamic order. They are described by the nilpotent part $N$ in the Weierstrass form.

**6. [Calculation] Find the generalized eigenvalues of $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ relative to $B = I$.**

??? success "Solution"
    **Calculation:**
    Since $B=I$, this reduces to the standard eigenvalue problem for $A$.
    $\det(A - \lambda I) = \lambda^2 - 1 = 0 \implies \lambda = \pm 1$.

**7. [Weierstrass] In the block $I - \lambda N$ of the Weierstrass form, what properties does $N$ satisfy?**

??? success "Solution"
    **Property:**
    $N$ is a **nilpotent matrix** ($N^k = O$ for some $k$). Its non-zero entries are on the super-diagonal, representing the Jordan chains associated with the infinite eigenvalue.

**8. [Stability] If the real parts of all finite generalized eigenvalues are negative and there are no infinite eigenvalues, is the system stable?**

??? success "Solution"
    **Yes.**
    This ensures that the differential part of the evolution decays and there are no impulsive terms or order mismatches caused by algebraic constraints.

**9. [Diagonalization] Under what condition can two matrices $A$ and $B$ be simultaneously diagonalized?**

??? success "Solution"
    **Conclusion:**
    If $A$ and $B$ are both Hermitian and one of them (usually $B$) is positive definite, they can be simultaneously diagonalized. This is the basis for "modal decomposition" in mechanical vibrations.

**10. [Application] How are generalized eigenvalues used to find transmission zeros in control theory?**

??? success "Solution"
    **Connection:**
    Transmission zeros are defined as the complex values $s$ that make the system matrix pencil (the Rosenbrock matrix) lose rank. This is essentially searching for generalized eigenvalues within a specific matrix pencil structure.

## Chapter Summary

Regular matrix pencil theory is the ultimate map for generalized linear systems:

1.  **Relativity of Spectra**: It elevates eigenvalues from an inherent property of one operator to a measure of interference between two, establishing a universal framework for relative evolution in physical systems.
2.  **Parsing Infinity**: By introducing infinite eigenvalues and nilpotent structures, pencil theory perfectly captures abrupt changes, constraints, and singularities in continuous systems.
3.  **Standardization of Structure**: The Weierstrass Canonical Form provides the ultimate reference for coordinate transformations in complex DAEs, enabling complete decoupling of dynamical modes.
