# Chapter 34: Schur Complement

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Positive Definite Matrices (Ch16) · Determinants (Ch03)

**Chapter Outline**: Motivation for Block Elimination → Definition of the Schur Complement → Determinant Decomposition Formula → Inversion Formula for Block Matrices (Banachiewicz Formula) → Criteria for Positive Definiteness → Inertia Additivity (Haynsworth Formula) → Applications in Matrix Equations → Statistical Significance (Partial Correlation and Conditional Variance)

**Extension**: The Schur complement is the "scalpel" of block matrix computation; it is not only the core of divide-and-conquer solvers for large linear systems but also the mathematical key to understanding Gaussian processes and Kernel methods (Ch29) in modern probability theory.

</div>

In handling large-scale systems, we frequently partition matrices into sub-blocks. The **Schur Complement** is the critical intermediate structure obtained through local elimination. It provides not only explicit expressions for the determinant and inverse of a block matrix but also deeply reveals the correlations between different components of the matrix.

---

## 34.1 Definition and Determinant Formula

!!! definition "Definition 34.1 (Schur Complement)"
    Let $M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$ be a partitioned matrix. If $A$ is invertible, the **Schur Complement** of $A$ in $M$ is defined as:
    $$S = D - C A^{-1} B$$

!!! theorem "Theorem 34.1 (Determinant Decomposition)"
    The determinant of the block matrix $M$ can be factored as:
    $$\det(M) = \det(A) \det(D - C A^{-1} B)$$
    This implies that the "volume" of the entire system equals the volume of the top-left subsystem multiplied by the volume of its Schur complement.

---

## 34.2 Inversion Formulas for Block Matrices

!!! technique "Technique: Banachiewicz Formula"
    If both $A$ and the Schur complement $S$ are invertible, the inverse of $M$ is:
    $$M^{-1} = \begin{pmatrix} A^{-1} + A^{-1} B S^{-1} C A^{-1} & -A^{-1} B S^{-1} \\ -S^{-1} C A^{-1} & S^{-1} \end{pmatrix}$$
    This is the numerical cornerstone for solving large structured systems.

---

## 34.3 Positive Definiteness and Inertia

!!! theorem "Theorem 34.2 (Positive Definiteness Criterion)"
    Let $M = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$ be symmetric. Then:
    $$M \succ 0 \iff A \succ 0 \text{ and } C - B^T A^{-1} B \succ 0$$

!!! theorem "Theorem 34.3 (Haynsworth Inertia Formula)"
    The inertia (number of positive, negative, and zero eigenvalues) of $M$ is the sum of the inertia of $A$ and the inertia of its Schur complement $S$:
    $$\operatorname{In}(M) = \operatorname{In}(A) + \operatorname{In}(D - C A^{-1} B)$$

---

## Exercises

1.  **[Calculation] Find the Schur complement of $\begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}$.**
    ??? success "Solution"
        Take $A=(1)$. $S = 5 - 2(1)^{-1}2 = 1$.

2.  **[Determinant] Use the Schur complement to calculate the determinant of the matrix from the previous exercise.**
    ??? success "Solution"
        $\det(M) = \det(A)\det(S) = 1 \cdot 1 = 1$.

3.  **[Property] If $B=0$, what is the Schur complement of $A$ in $M$?**
    ??? success "Solution"
        $S = D - C A^{-1} 0 = D$. When blocks are decoupled, the Schur complement is simply the original diagonal block.

4.  **[Symmetry] Prove that if $M$ is symmetric, its Schur complement is also symmetric.**
    ??? success "Solution"
        $S^T = (D - C A^{-1} B)^T = D^T - B^T (A^{-1})^T C^T = D - C A^{-1} B = S$ (for symmetric $M$, $C = B^T$ and $D^T = D$).

5.  **[Rank] Prove $\operatorname{rank}(M) = \operatorname{rank}(A) + \operatorname{rank}(S)$.**
    ??? success "Solution"
        Using block elimination: $\begin{pmatrix} I & 0 \\ -CA^{-1} & I \end{pmatrix} \begin{pmatrix} A & B \\ C & D \end{pmatrix} \begin{pmatrix} I & -A^{-1}B \\ 0 & I \end{pmatrix} = \begin{pmatrix} A & 0 \\ 0 & S \end{pmatrix}$. Since the matrices on the left and right are non-singular, the rank is preserved.

6.  **[Application] In statistics, if $M$ is a covariance matrix, what does $S$ represent?**
    ??? success "Solution"
        $S$ represents the **conditional covariance** of the set of variables $D$ given the set of variables $A$.

7.  **[Dual] Write the Schur complement of $D$ in $M$.**
    ??? success "Solution"
        $S_D = A - B D^{-1} C$ (assuming $D$ is invertible).

8.  **[Singular Values] Is there a simple additive relationship between the singular values of $M$ and those of $A$ and $S$?**
    ??? success "Solution"
        No. Singular values do not possess the simple block-additivity properties that eigenvalues or inertia indices do.

9.  **[Inverse] If $M$ is a block lower triangular matrix ($B=0$), find its inverse.**
    ??? success "Solution"
        $\begin{pmatrix} A^{-1} & 0 \\ -D^{-1} C A^{-1} & D^{-1} \end{pmatrix}$.

****

??? success "Solution"
    

## Chapter Summary

The Schur complement is the core syntax of block algebra:

1.  **Dimensional Folding**: It demonstrates how to compress the complexity of a high-dimensional system into a low-dimensional remainder using local invertibility—the mathematical prerequisite for distributed algorithms.
2.  **Transmission of Positivity**: The Schur complement formula establishes quantitative links between the global stability of a block matrix and the local stability of its subsystems and interaction terms.
3.  **Statistical Bridge**: In probability theory, the Schur complement reveals the algebraic essence of the conditioning process in Gaussian distributions, serving as the ultimate tool for handling "remaining correlation" between variables.
