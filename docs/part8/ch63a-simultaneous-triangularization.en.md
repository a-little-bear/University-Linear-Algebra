# Chapter 63A: Simultaneous Triangularization

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Lie Algebras (Ch55B) · Schur Decomposition (Ch10) · Matrix Groups (Ch55)

**Chapter Outline**: Definition of Simultaneous Triangularization (ST) → Commutativity and the ST Property → Existence of Common Eigenvectors → McCoy’s Theorem (ST for Commuting Families) → Engel’s Theorem (Nilpotent Lie Algebras) → Lie-Kolchin Theorem (Solvable Groups) → Criteria for Simultaneous Diagonalization (Commuting Normal Matrices) → Applications: Shared Observables in Quantum Mechanics, Decoupling Analysis of Linear Systems, and Representation Theory

**Extension**: Simultaneous triangularization investigates how a family of matrices can "harmoniously coexist"; it reveals the hidden hierarchical structure within operator families and is key to understanding how non-commutative algebras can degenerate into quasi-commutative structures through coordination mechanisms.

</div>

In studying individual matrices, Schur's theorem guarantees that every matrix can be triangularized. But what if we face a set of matrices $\{A_1, A_2, \ldots\}$? Does there exist a single basis that simultaneously transforms all of them into upper triangular form? This is the problem of **Simultaneous Triangularization** (ST). it directly pertains to whether multiple physical quantities can be measured simultaneously and whether multiple control loops can be decoupled.

---

## 63A.1 Commutativity and Common Eigenvectors

!!! definition "Definition 63A.1 (Simultaneous Triangularization)"
    A family of matrices $\mathcal{A} \subset M_n(\mathbb{C})$ is **simultaneously triangularizable** if there exists a unitary matrix $U$ such that $U^* A U$ is upper triangular for every $A \in \mathcal{A}$.

!!! theorem "Theorem 63A.1 (Common Eigenvectors)"
    If a family of matrices $\mathcal{A}$ is commutative (i.e., $AB = BA$ for all $A, B \in \mathcal{A}$), then they share at least one common eigenvector.
    **Significance**: This is the logical starting point for ST, allowing the space to be compressed dimension by dimension via induction.

---

## 63A.2 McCoy’s and Engel’s Theorems

!!! theorem "Theorem 63A.2 (McCoy’s Theorem)"
    A finite set of matrices $A_1, \ldots, A_k$ is simultaneously triangularizable iff for every polynomial $p(x_1, \ldots, x_k)$, the matrix $p(A_1, \ldots, A_k)(A_i A_j - A_j A_i)$ is nilpotent.
    **Corollary**: Any set of commuting matrices is always simultaneously triangularizable.

!!! theorem "Theorem 63A.3 (Engel’s Theorem)"
    If every element of a matrix Lie algebra $\mathfrak{g}$ is nilpotent, then $\mathfrak{g}$ is simultaneously strictly upper triangularizable (all diagonal entries are 0).

---

## 63A.3 Simultaneous Diagonalization

!!! theorem "Theorem 63A.4 (Diagonalization Criterion)"
    A family of diagonalizable matrices $\mathcal{A}$ can be simultaneously diagonalized iff they are pairwise commutative.
    **Quantum Application**: In physics, this implies that the corresponding observable operators share a complete set of eigenstates and can be measured accurately at the same time.

---

## Exercises

1.  **[Basics] Determine if $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ and $B = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ are simultaneously triangularizable.**
    ??? success "Solution"
        Yes. They are already upper triangular in the standard basis. Also, one can verify that $AB = BA$.

2.  **[Commutativity] Give an example of two matrices that do not commute but are simultaneously triangularizable.**
    ??? success "Solution"
        $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}, B = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$. They do not commute, yet both are upper triangular in the standard basis.

3.  **[Common Vector] Prove: If $AB=BA$ and $\lambda$ is a simple eigenvalue of $A$ with eigenvector $x$, then $x$ is also an eigenvector of $B$.**
    ??? success "Solution"
        $A(Bx) = B(Ax) = B(\lambda x) = \lambda (Bx)$. Since $\lambda$ is a simple root, the eigenspace is 1-dimensional. Thus $Bx$ must be proportional to $x$, meaning $Bx = \mu x$.

4.  **[Diagonalization] Can two non-normal matrices be simultaneously diagonalized if they commute?**
    ??? success "Solution"
        Not necessarily. For simultaneous diagonalization, each matrix must first be diagonalizable. If one is defective (like a Jordan block), it cannot be diagonalized at all.

5.  **[Dimension] If an $n \times n$ matrix family is simultaneously triangularizable, what is the maximum dimension of the algebra they span?**
    ??? success "Solution"
        $n(n+1)/2$, which is the dimension of the space of all upper triangular matrices.

6.  **[Lie-Kolchin] What algebraic structure does the Lie-Kolchin theorem apply to?**
    ??? success "Solution"
        It applies to linear algebraic groups, stating that any solvable linear algebraic group over an algebraically closed field is simultaneously triangularizable.

7.  **[Property] Prove: If a family is simultaneously triangularizable, their commutators $[A, B]$ must be nilpotent.**
    ??? success "Solution"
        The commutator of two upper triangular matrices is strictly upper triangular (zeros on the diagonal). Any strictly upper triangular matrix is nilpotent.

8.  **[Calculation] Can $\begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$ and $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ be simultaneously triangularized?**
    ??? success "Solution"
        No. Both are normal matrices (and symmetric with distinct eigenvalues). If they were simultaneously triangularizable, they would be simultaneously diagonalizable, which requires commutativity. However, they do not commute.

9.  **[Control] How does simultaneous triangularization help in system decoupling?**
    ??? success "Solution"
        If the state transition matrices can be simultaneously triangularized, a unified coordinate transform can turn a multivariable system into a sequence of unidirectionally coupled subsystems, simplifying controller design.

10. **[Rank] Is the existence of a common eigenvector related to the rank of the matrices?**

   ??? success "Solution"
        Not directly. It depends on the commutation structure of the operators and the stability of their actions on the space.

## Chapter Summary

Simultaneous triangularization establishes the collective order of matrix sets:

1.  **From Individual to Group**: It generalizes Schur decomposition from single operators to operator families, revealing the possibility of coordinating multiple transformations within a single physical context.
2.  **Depth of Commutativity**: It proves that commutativity is a necessary and sufficient condition for "structural harmony" (for diagonalization) and a sufficient condition for triangularization, serving as a core principle in modern operator algebra.
3.  **Discovery of Hierarchy**: The triangular structure essentially defines a filtration, proving that even amidst complex interactions, there often exists a common, step-by-step sequence of invariant subspaces.
