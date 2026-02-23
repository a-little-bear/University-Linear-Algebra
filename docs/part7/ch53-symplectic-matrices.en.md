# Chapter 53: Symplectic and Hamiltonian Matrices

<div class="context-flow" markdown>

**Prerequisites**: Dual Spaces (Ch13A) · Matrix Groups (Ch55) · Positive Definite Matrices (Ch16)

**Chapter Outline**: Symplectic Forms & the Standard Antisymmetric Matrix $J$ → Definition and Properties of Symplectic Matrices → The Symplectic Group $Sp(2n, \mathbb{R})$ → Spectral Symmetry (Eigenvalues in Reciprocal Pairs) → Hamiltonian Matrices & Lie Algebras → Exponential Mapping & Cayley Transforms → Williamson's Theorem (Symplectic Diagonalization) → Applications: Hamiltonian Mechanics (Conservation of Phase Space), Symplectic Integrators, and Quantum Gaussian States

**Extension**: Symplectic matrices are the algebraic language for "Hamiltonian evolution" in physics; they preserve the area (or volume) of phase space and are indispensable for studying conservative systems, robotics, and quantum computing.

</div>

While orthogonal transformations preserve Euclidean distance, there is a class of matrices that preserves a more subtle "area structure." These are known as **Symplectic Matrices**. They are deeply connected to Hamiltonian mechanics, revealing the intrinsic constraints of dynamical systems in phase space. This chapter establishes the symplectic form, explores the duality between symplectic and Hamiltonian matrices, and highlights their role in modern numerical computation.

---

## 53.1 Symplectic Forms and the Standard Operator $J$

!!! definition "Definition 53.1 (Standard Symplectic Matrix $J$)"
    In $2n$-dimensional space, we define the standard skew-symmetric matrix $J$:
    $$J = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}$$
    It satisfies $J^2 = -I$ and $J^T = -J = J^{-1}$.

!!! definition "Definition 53.2 (Symplectic Matrix)"
    A matrix $M \in M_{2n}(\mathbb{R})$ is **Symplectic** if it satisfies:
    $$M^T J M = J$$
    This means $M$ preserves the symplectic bilinear form $\omega(u, v) = u^T J v$.

---

## 53.2 The Symplectic Group and Spectral Properties

!!! theorem "Theorem 53.1 (Properties of Symplectic Matrices)"
    1.  **Determinant**: All symplectic matrices have $\det(M) = 1$.
    2.  **Spectral Symmetry**: If $\lambda$ is an eigenvalue of $M$, then $1/\lambda, \bar{\lambda}, \text{ and } 1/\bar{\lambda}$ are also eigenvalues.
    3.  **Group Structure**: The set of $2n \times 2n$ symplectic matrices forms the **Symplectic Group** $Sp(2n, \mathbb{R})$ under multiplication.

---

## 53.3 Hamiltonian Matrices

!!! definition "Definition 53.3 (Hamiltonian Matrix)"
    A matrix $H$ is **Hamiltonian** if it satisfies:
    $$(JH)^T = JH$$
    This is equivalent to $H^T J + JH = 0$.
    **Relationship**: Hamiltonian matrices form the Lie algebra $\mathfrak{sp}(2n)$ associated with the symplectic group.

---

## 53.4 Williamson's Theorem

!!! theorem "Theorem 53.2 (Williamson's Theorem)"
    For any $2n \times 2n$ real positive definite matrix $A \succ 0$, there exists a symplectic matrix $M$ such that:
    $$M^T A M = \operatorname{diag}(d_1, \ldots, d_n, d_1, \ldots, d_n)$$
    The values $d_i > 0$ are known as the **Symplectic Eigenvalues** of $A$.

---

## Exercises

1.  **[Calculation] Verify if $J = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ is itself a symplectic matrix.**
    ??? success "Solution"
        $J^T J J = (-J)J^2 = (-J)(-I) = J$. The condition is satisfied.

2.  **[Determinant] Prove: If $M$ is symplectic, it must be non-singular.**
    ??? success "Solution"
        From $M^T J M = J$, taking the determinant yields $\det(M)^2 \det(J) = \det(J)$. Since $\det(J)=1 \neq 0$, it follows that $\det(M)^2=1 \implies \det(M) = \pm 1$. (A deeper proof shows $\det M$ must be +1).

3.  **[Spectral] Explain why a symplectic matrix cannot have 2 as its only eigenvalue.**
    ??? success "Solution"
        Because eigenvalues must appear in reciprocal pairs; if 2 is an eigenvalue, $1/2 = 0.5$ must also be an eigenvalue.

4.  **[Hamiltonian] When is $\begin{pmatrix} A & B \\ C & -A^T \end{pmatrix}$ a Hamiltonian matrix (assume $B, C$ are symmetric)?**
    ??? success "Solution"
        Always. By checking $H^T J + JH = 0$ or observing the symmetry of $JH = \begin{pmatrix} C & -A^T \\ -A & -B \end{pmatrix}$.

5.  **[Exp] Prove: If $H$ is Hamiltonian, then $e^H$ is symplectic.**
    ??? success "Solution"
        $(e^H)^T J e^H = e^{H^T} J e^H = J e^{-H} e^H = J$. This utilizes the identity $H^T J = -JH \implies H^T = -JHJ^{-1}$.

6.  **[Williamson] Find the symplectic eigenvalues of $A = \operatorname{diag}(2, 2)$.**
    ??? success "Solution"
        Since it is already in the form required by Williamson's theorem, the symplectic eigenvalue is $d_1 = 2$.

7.  **[Inverse] Prove that the inverse of a symplectic matrix is also symplectic.**
    ??? success "Solution"
        Invert both sides of $M^T J M = J$ to get $M^{-1} J^{-1} (M^T)^{-1} = J^{-1}$. Substituting $J^{-1} = -J$ yields $(M^{-1})^T J M^{-1} = J$.

8.  **[Numerical] What is a Symplectic Integrator?**
    ??? success "Solution"
        A numerical algorithm for solving Hamiltonian equations where each iteration matrix is symplectic, ensuring the numerical solution preserves system energy and phase space volume over long periods.

9.  **[Intersection] What is the intersection of the symplectic group and the unitary group?**
    ??? success "Solution"
        The intersection is a subgroup isomorphic to $U(n)$.

****

??? success "Solution"
    

## Chapter Summary

Symplectic and Hamiltonian matrices define the algebraic skeleton of conservative systems:

1.  **Geometric Invariants**: Symplectic matrices are operators that preserve "symplectic area" in differential geometry, revealing rigid constraints in physical evolution.
2.  **Spectral Harmony**: The four-fold symmetry of eigenvalues reflects the deep connection between energy conservation and time-reversal symmetry in Hamiltonian systems.
3.  **Numerical Guardians**: Via Williamson's theorem and symplectic algorithms, linear algebra provides the ultimate guarantee of "structure preservation" for simulations of complex mechanical systems, preventing physical errors caused by numerical dissipation.
