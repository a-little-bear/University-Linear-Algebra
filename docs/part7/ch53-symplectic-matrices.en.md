# Chapter 53A: Symplectic and Hamiltonian Matrices

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces (Ch08) · Matrix Equations (Ch20) · Classical Mechanics Basics

**Chapter Outline**: From Euclidean to Symplectic Geometry → The Symplectic Form and Standard Skew-symmetric Matrix $J$ → Definition and Group Structure of Symplectic Matrices ($Sp(2n)$) → Spectral Properties of Symplectic Matrices (Eigenvalue Pairing) → Definition and Lie Algebra Structure of Hamiltonian Matrices → The Link between Hamiltonian Matrices and Riccati Equations → Symplectic Decompositions (Symplectic QR and Schur) → Applications: Volume Preservation in Phase Space, State-Costate Equations in Optimal Control, and Quantum Optics

**Extension**: Symplectic matrices are operators that describe "area-preserving" transformations; they reveal the inherent algebraic symmetries of physical systems under energy conservation, serving as the advanced algebraic language connecting classical mechanics, PDEs, and control theory.

</div>

In linear algebra, we are familiar with orthogonal matrices that preserve inner products. In physics, particularly Hamiltonian mechanics, it is more important to preserve a "skew-symmetric inner product." **Symplectic Matrices** are the algebraic characterization of such operators. they describe the conservation of area in phase space and lead to a class of matrices with perfectly symmetric eigenvalue structures: **Hamiltonian Matrices**. This chapter explores these special matrices that profoundly influence physical evolution and optimal control.

---

## 53A.1 Symplectic Form and Matrices

!!! definition "Definition 53A.1 (Standard Symplectic Matrix $J$)"
    Define the $2n \times 2n$ block matrix $J$ as:
    $$J = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}$$
    satisfying $J^2 = -I$ and $J^T = -J$. It defines the symplectic inner product $[x, y] = x^T J y$.

!!! definition "Definition 53A.2 (Symplectic Matrix)"
    A square matrix $M \in M_{2n}$ is **Symplectic** if it preserves the symplectic form:
    $$M^T J M = J$$
    The set of all such matrices forms the **Symplectic Group** $Sp(2n)$.

---

## 53A.2 Hamiltonian Matrices

!!! definition "Definition 53A.3 (Hamiltonian Matrix)"
    A square matrix $H \in M_{2n}$ is **Hamiltonian** if $JH$ is symmetric:
    $$(JH)^T = JH \iff H^T J + JH = 0$$
    **Property**: If $M(t)$ is a symplectic trajectory with $M(0)=I$, its derivative $\dot{M}(0)$ is necessarily a Hamiltonian matrix.

---

## 53A.3 Spectral Properties and Pairing

!!! theorem "Theorem 53A.1 (Eigenvalue Pairing)"
    1.  **Symplectic**: If $\lambda$ is an eigenvalue, then $1/\lambda$ must also be an eigenvalue.
    2.  **Hamiltonian**: If $\lambda$ is an eigenvalue, then $-\lambda$ must also be an eigenvalue.
    **Physical Meaning**: This pairing reflects the symmetric balance between stable and unstable modes in conservative physical systems.

---

## Exercises

**1. [Basics] Calculate the determinant of a $2 \times 2$ symplectic matrix.**

??? success "Solution"
    **Conclusion: 1.**
    **Proof**: From $M^T J M = J$, taking the determinant of both sides gives:
    $\det(M^T) \det(J) \det(M) = \det(J)$. Since $\det(J) = 1$ for $2n \times 2n$, we have $(\det M)^2 = 1 \implies \det M = \pm 1$. More advanced theory (Pfaffians) proves the determinant must be strictly +1, corresponding to volume preservation in phase space.

**2. [Hamiltonian] Determine if $H = \begin{pmatrix} A & B \\ C & -A^T \end{pmatrix}$ is Hamiltonian (with $B, C$ symmetric).**

??? success "Solution"
    **Verification Steps:**
    1. $JH = \begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix} \begin{pmatrix} A & B \\ C & -A^T \end{pmatrix} = \begin{pmatrix} C & -A^T \\ -A & -B \end{pmatrix}$.
    2. Transpose: $(JH)^T = \begin{pmatrix} C^T & -A^T \\ -A & -B^T \end{pmatrix}$.
    3. Since $B, C$ are symmetric, $C^T=C$ and $B^T=B$.
    **Conclusion**: $(JH)^T = JH$, so this block form is the standard construction for a Hamiltonian matrix.

**3. [Basics] Verify if $J = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ itself is a symplectic matrix.**

??? success "Solution"
    **Calculation:**
    $J^T J J = (-J) J J = (-J)(-I) = J$.
    **Conclusion**: Yes. In 2D, $J$ is both symplectic and Hamiltonian.

**4. [Spectral] If a Hamiltonian matrix has an eigenvalue $2+3i$, find the other three eigenvalues that must exist.**

??? success "Solution"
    **Using Pairing Properties:**
    1. Being real, it must have the conjugate: $2-3i$.
    2. Being Hamiltonian, it must have the negative: $-(2+3i) = -2-3i$.
    3. The conjugate of the negative: $-2+3i$.
    **Conclusion**: Eigenvalues appear in quartets $\{ \pm \lambda, \pm \bar{\lambda} \}$.

**5. [Area Preservation] Prove that a symplectic operator in $\mathbb{R}^2$ preserves the area of a parallelogram.**

??? success "Solution"
    **Proof:**
    1. Area is given by the symplectic form $x^T J y$.
    2. After transformation: $(Mx)^T J (My) = x^T (M^T J M) y$.
    3. By definition of a symplectic matrix, this equals $x^T J y$.
    **Conclusion**: The signed area is invariant under symplectic transforms.

**6. [Riccati] How are Hamiltonian matrices used to solve Algebraic Riccati Equations (ARE)?**

??? success "Solution"
    **Algebraic Link:**
    The solution $X$ to an ARE corresponds to the **stable invariant subspace** of its associated Hamiltonian matrix $H$. By finding the eigenvectors of $H$ and performing block combinations, the positive definite solution can be constructed directly, bypassing non-linear iterations.

**7. [Inverse] Prove: If $M$ is symplectic, then its inverse $M^{-1}$ is also symplectic.**

??? success "Solution"
    **Proof:**
    1. Given $M^T J M = J$.
    2. Left-multiply by $(M^T)^{-1}$ and right-multiply by $M^{-1}$:
    3. $J = (M^T)^{-1} J M^{-1} = (M^{-1})^T J M^{-1}$.
    **Conclusion**: This satisfies the definition of a symplectic matrix.

**8. [Lie Algebra] Prove: The commutator $[H_1, H_2]$ of two Hamiltonian matrices is Hamiltonian.**

??? success "Solution"
    **Reasoning:**
    This is because Hamiltonian matrices form the **Lie Algebra** $\mathfrak{sp}(2n)$ of the symplectic group $Sp(2n)$. A Lie algebra is by definition closed under the commutator bracket.

**9. [Basics] Write the inverse of the $4 \times 4$ standard symplectic matrix $J$.**

??? success "Solution"
    **Conclusion: $-J$.**
    **Reasoning**: Since $J^2 = -I$, $J(-J) = I$. In physics, this corresponds to reversing the orientation of the symplectic rotation.

**10. [Application] Briefly state the advantage of "Symplectic Integrators" in orbital calculations.**

??? success "Solution"
    **Reasoning:**
    1. Traditional numerical methods (like Runge-Kutta) accumulate energy errors over time, causing planetary orbits to drift.
    2. Symplectic integrators strictly enforce $M^T J M = J$ at every step.
    3. This ensures the Hamiltonian (total energy) of the system remains almost constant, preserving the closure and stability of orbits over vast time scales.

## Chapter Summary

Symplectic and Hamiltonian matrices reveal the deep geometric order of dynamical systems:

1.  **Guardians of Area**: Symplectic matrices ensure the incompressibility of flow in phase space, serving as the algebraic root of Liouville's theorem in statistical mechanics.
2.  **Symmetry of Spectra**: The $\pm \lambda$ pairing of Hamiltonian eigenvalues perfectly characterizes local stability near equilibrium in conservative systems, providing a natural basis for identifying stable manifolds.
3.  **Bridge to Control**: Through the link with Riccati equations, symplectic algebra enables the leap from abstract physical symmetry to engineering optimal control solutions, serving as a powerful pillar of modern systems science.
