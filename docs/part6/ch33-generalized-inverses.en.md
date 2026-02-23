# Chapter 33: Generalized Inverses

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Singular Value Decomposition (Ch11) · Linear Equations (Ch01)

**Chapter Outline**: Limitations of the Standard Inverse → Definition of the Moore-Penrose Inverse ($A^+$) via the 4 Penrose Conditions → Existence and Uniqueness → SVD-based Computation → Least Squares and Minimum Norm Solutions → Drazin Inverse & Group Inverse (Handling Nilpotent Structures) → Inverses of Block Matrices → Applications: Singular Differential Equations and Ill-conditioned Systems

**Extension**: Generalized inverses break the shackles of "square and full rank" for inversion, serving as the final mathematical judgment for handling all linear uncertainties (no solution or infinite solutions); it is the cornerstone of linear model estimation (Gauss-Markov Theorem) in statistics.

</div>

In classical matrix algebra, only non-singular square matrices have inverses. However, in engineering, statistics, and control theory, we frequently encounter rectangular or rank-deficient matrices. **Generalized Inverses** break this restriction, defining a "inverse" in some sense for every matrix. Among them, the Moore-Penrose inverse stands out for its perfect behavior regarding solution minimization.

---

## 33.1 The Moore-Penrose Inverse $A^+$

!!! definition "Definition 33.1 (Penrose Conditions)"
    For a matrix $A \in M_{m \times n}(\mathbb{C})$, the **Moore-Penrose Inverse** is the unique matrix $A^+$ satisfying the following four conditions:
    1.  $AA^+A = A$ (Inner consistency)
    2.  $A^+AA^+ = A^+$ (Outer consistency)
    3.  $(AA^+)^* = AA^+$ (Left hermiticity)
    4.  $(A^+A)^* = A^+A$ (Right hermiticity)

!!! theorem "Theorem 33.1 (Existence and Uniqueness)"
    For any matrix $A$, a matrix $A^+$ satisfying the four Penrose conditions **exists and is unique**.

---

## 33.2 Computation and Least Squares

!!! technique "Computation: SVD Approach"
    If $A = U \Sigma V^*$, then $A^+ = V \Sigma^+ U^*$, where $\Sigma^+$ is obtained by taking the reciprocal of the non-zero diagonal entries of $\Sigma$ and transposing.

!!! theorem "Theorem 33.2 (Best Approximation Properties)"
    For the linear system $Ax = b$:
    1.  The vector $\hat{x} = A^+ b$ minimizes $\|Ax - b\|_2$ (the Least Squares solution).
    2.  Among all least squares solutions, $\hat{x} = A^+ b$ has the minimum norm $\|\hat{x}\|_2$ (the Minimum Norm solution).

---

## 33.3 The Drazin Inverse

!!! definition "Definition 33.2 (Drazin Inverse)"
    For a square matrix $A$, the **Drazin Inverse** $A^D$ satisfies:
    1.  $A^k A^D A = A^k$ (where $k$ is the index of $A$, the smallest power such that $\operatorname{rank}(A^k)$ stabilizes)
    2.  $A^D A A^D = A^D$
    3.  $AA^D = A^D A$
    **Application**: The Drazin inverse is highly effective in handling singular differential equations and convergence analysis of Markov chains.

---

## Exercises


****
??? success "Solution"
     $A^+ = \begin{pmatrix} 0.5 & 0 \\ 0 & 0 \end{pmatrix}$. Verify: $AA^+ = \operatorname{diag}(1, 0)$, which is Hermitian and satisfies the conditions.


****
??? success "Solution"
     Substitute into the four Penrose conditions. For example, $A^+ A = (A^* A)^{-1} A^* A = I$, which is clearly Hermitian and satisfies the other requirements.


****
??? success "Solution"
     Since $A$ satisfies the four Penrose conditions for $A^+$ (with roles swapped), and $A^{++}$ is unique, they must be equal.


****
??? success "Solution"
     From $AA^+A=A$, $\operatorname{rank}(A) \le \operatorname{rank}(A^+)$. From $A^+AA^+=A^+$, $\operatorname{rank}(A^+) \le \operatorname{rank}(A)$. Thus the ranks are equal.


****
??? success "Solution"
     It represents the orthogonal projection matrix onto the column space $C(A)$.


****
??? success "Solution"
     If consistent, $b \in C(A)$. Then $A(A^+ b) = (AA^+)b = b$, as the projection of a vector already in the space is the vector itself.


****
??? success "Solution"
     A Group Inverse is the Drazin inverse when the index $k=1$. It satisfies $AA^{\#}A = A, A^{\#}AA^{\#} = A^{\#}$, and $AA^{\#} = A^{\#}A$.


****
??? success "Solution"
     $1/2, 1, 0$.


****
??? success "Solution"
     Verify the Penrose conditions for $A^*$ and $(A^+)^*$; the symmetry of the conditions ensures the result.

****
??? success "Solution"
    ## Chapter Summary

The generalized inverse is the ultimate cure for linear pathology:


****: The Moore-Penrose inverse, through four symmetric consistency axioms, finds a unique and optimally behaving "pseudo-inverse" for every linear operator, unifying the theories of invertible and singular matrices.

****: Its dual optimization properties regarding least squares and minimum norms make it the standard mathematical tool for solving inconsistent and underdetermined systems.

****: The introduction of Drazin and Group inverses demonstrates how to use generalized inverses to handle the nilpotent parts and complex index structures of matrices, providing algebraic leverage for analyzing singular dynamical systems.
