# Chapter 63A: Simultaneous Triangularization

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Schur Decomposition (Ch10) · Matrix Groups and Lie Algebras (Ch55)

**Chapter Outline**: From Single Operators to Operator Families → Definition of Simultaneous Triangularization → Core Criterion: Commuting Families → Proof of Existence of Common Eigenvectors → Necessary and Sufficient Conditions for Simultaneous Diagonalization → Lie Theoretic Background: Engel’s and Lie’s Theorems → Consistency of Structure in Linear Systems → Applications: Common Observables in Quantum Mechanics (Commutators), Decoupling of Dynamical Systems, and Structure Analysis of Matrix Algebras

**Extension**: Simultaneous triangularization studies how a set of matrices can be "ordered" through a single coordinate transformation; it reveals that commutativity is not just an algebraic property but a geometric "consistency of direction." It is the core algebraic framework for understanding systems with multi-variable interactions.

</div>

When dealing with multiple matrices, the ideal scenario is to find a single basis where all matrices simultaneously exhibit simple structures (such as upper triangular or diagonal). **Simultaneous Triangularization** studies the conditions under which a set of operators shares a common eigen-structure. This theory proves that as long as operators "understand" each other (i.e., they commute), they can be simplified together. This chapter introduces this fundamental law at the heart of quantum mechanics and matrix analysis.

---

## 63A.1 Common Eigenvectors

!!! theorem "Theorem 63A.1 (Existence of Common Eigenvectors)"
    Let $\mathcal{F}$ be a family of commuting linear operators on $V$. If each operator has an eigenvalue in the field (e.g., over $\mathbb{C}$), then there exists a **common eigenvector** for all operators in $\mathcal{F}$.
    **Physical Intuition**: Commuting physical observables can be measured simultaneously with precision.

---

## 63A.2 Simultaneous Triangularization and Diagonalization

!!! definition "Definition 63A.1 (Simultaneous Triangularization)"
    A family of matrices $\{A_i\}$ is simultaneously triangularizable if there exists a single non-singular matrix $P$ such that $P^{-1} A_i P$ is upper triangular for all $i$.

!!! theorem "Theorem 63A.2 (Criteria)"
    1.  **Triangularization**: A family of complex matrices is simultaneously triangularizable iff the matrices are **commuting** (or more generally, they span a solvable Lie algebra).
    2.  **Diagonalization**: A family of matrices is simultaneously diagonalizable iff each matrix is individually diagonalizable and the family is pairwise commuting.

---

## 63A.3 The Lie Theoretic Perspective

!!! note "Engel’s and Lie’s Theorems"
    - **Lie’s Theorem**: Any finite-dimensional representation of a solvable Lie algebra over $\mathbb{C}$ has a common eigenvector (and thus can be simultaneously triangularized).
    - This extends the scope from simple commutativity to operators with a hierarchical structure (solvable algebras).

---

## Exercises

**1. [Basics] Determine if $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ and $B = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ are simultaneously triangularizable.**

??? success "Solution"
    **Steps:**
    1. Check commutativity: $AB = \begin{pmatrix} 2 & 3 \\ 0 & 2 \end{pmatrix}$ and $BA = \begin{pmatrix} 2 & 3 \\ 0 & 2 \end{pmatrix}$.
    2. Since $AB=BA$ and both are already upper triangular in the standard basis.
    **Conclusion**: Yes, they are already simultaneously upper triangularized.

**2. [Calculation] Find a common eigenvector for $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$ and $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.**

??? success "Solution"
    **Analysis:**
    1. Check the commutator: $[A, B] = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \neq O$.
    2. Since they do not commute, they do not necessarily share an eigenvector.
    3. $A$'s eigenvectors are $e_1, e_2$.
    4. $B$'s eigenvectors are $(1, 1)^T$ and $(1, -1)^T$.
    **Conclusion**: They share no common eigenvectors.

**3. [Diagonalization] Prove: If $A$ and $B$ are simultaneously diagonalizable, they must commute.**

??? success "Solution"
    **Proof:**
    1. Let $P^{-1}AP = D_1$ and $P^{-1}BP = D_2$ where $D_1, D_2$ are diagonal.
    2. Diagonal matrices always commute: $D_1 D_2 = D_2 D_1$.
    3. Substituting: $(P^{-1}AP)(P^{-1}BP) = (P^{-1}BP)(P^{-1}AP)$.
    4. Cancel $P$ and $P^{-1}$: $P^{-1}ABP = P^{-1}BAP \implies AB = BA$.

**4. [Property] Over $\mathbb{C}$, do two commuting matrices always simultaneously diagonalize?**

??? success "Solution"
    **Conclusion: Not necessarily.**
    **Example**: $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ and $B = I$. They commute, but $A$ is not diagonalizable. Thus, they can only be **simultaneously triangularized**, not diagonalized.

**5. [Dimension] If a family of $n \times n$ matrices commutes, what is the minimum dimension of their common eigenspace?**

??? success "Solution"
    **Conclusion: At least 1.**
    This is the starting point of the theory: a commuting family over an algebraically closed field always shares at least one invariant line (eigen-direction).

**6. [Normal Matrices] For a set of normal matrices ($AA^*=A^*A$), what does simultaneous triangularization imply?**

??? success "Solution"
    **Conclusion: It implies simultaneous diagonalization.**
    **Reasoning**: An upper triangular matrix that is also normal must be diagonal (Schur decomposition properties). Thus, commuting normal matrices are unitarily simultaneously diagonalizable.

**7. [Application] In quantum mechanics, what does the operator condition $[A, B] = 0$ physically represent?**

??? success "Solution"
    **Physical Context:**
    It represents that the two observables are **compatible**. According to simultaneous diagonalization, there exists a common basis of eigenstates. This means both quantities can be measured precisely at the same time without being limited by the Uncertainty Principle.

**8. [Calculation] Find a basis that simultaneously triangularizes $\begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}$ and $\begin{pmatrix} 4 & 5 \\ 0 & 6 \end{pmatrix}$.**

??? success "Solution"
    Since both matrices are already upper triangular, the standard basis $\{e_1, e_2\}$ is the basis that triangularizes them. The transformation matrix is $P = I$.

**9. [Lie’s Theorem] Briefly state the intuitive meaning of "solvability" in Lie's Theorem.**

??? success "Solution"
    In operator algebras, solvability means that while operators don't necessarily commute, they follow a "hierarchical structure" (commutators fall into smaller ideals). This hierarchy ensures that we can find common eigenspaces by peeling back layers of the algebra.

**10. [Stability] If a state equation is described by two commuting matrices $A$ and $B$, what is special about the solution $e^{At}e^{Bt}$?**

??? success "Solution"
    **Conclusion:**
    Due to commutativity, $e^{At}e^{Bt} = e^{(A+B)t}$. This means the two evolution processes decouple and superimpose linearly. the complex dynamics is equivalent to the evolution of the sum of the operators, greatly simplifying multivariate control analysis.

## Chapter Summary

Simultaneous triangularization is the ultimate algebraic expression of operator synergy:

1.  **The Commutativity Dividend**: It establishes commutativity as the algebraic guarantee of structural consistency, proving that compatible operators can be simplified within the same geometric frame.
2.  **Core of Quantum and Classical**: From Heisenberg commutators to classic vibration modes, simultaneous diagonalization provides the only valid path for decomposing coupled systems into independent components.
3.  **Hierarchy of Structure**: Through the extension of Lie's Theorem, the theory reveals that even in the absence of perfect commutativity, hierarchical dependencies can induce an ordered triangular structure—a cornerstone of modern representation theory.
