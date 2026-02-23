# Chapter 38A: M-Matrices and Z-Matrices

<div class="context-flow" markdown>

**Prerequisites**: Matrix Inversion (Ch2) · Eigenvalues (Ch6) · Nonnegative Matrices (Ch17) · Stability (Ch36)

**Chapter Outline**: Definition of Z-Matrices → Definition of M-Matrices → Characterizations of M-Matrices (Positivity of Inverse, Principal Minors, etc.) → Property of $A^{-1} \ge 0$ → M-Matrices and Stability → Singular M-Matrices → Comparison Theorem → Applications in Markov Chains and Economics

**Extension**: M-matrices are the fundamental tool for analyzing systems where negative off-diagonal interactions are balanced by strong positive self-regulation, ensuring the positivity of the solution.

</div>

**M-matrices** are a class of matrices whose off-diagonal entries are non-positive and whose principal minors are positive. They arise in contexts where "competition" or "consumption" (negative off-diagonals) is dominated by "production" or "supply" (positive diagonals). The defining characteristic of an M-matrix is that its **inverse is nonnegative** ($A^{-1} \ge 0$), ensuring that positive inputs always produce positive outputs. This makes them indispensable in the Leontief input-output model (Ch69) and the analysis of Markov chains.

---

## 38A.1 Z-Matrices and M-Matrices

!!! definition "Definition 38A.1 (Z-Matrix)"
    A square matrix $A$ is a **Z-matrix** if $a_{ij} \le 0$ for all $i 
eq j$.

!!! definition "Definition 38A.2 (M-Matrix)"
    A Z-matrix $A$ is an **M-matrix** if it satisfies any of the following equivalent conditions:
    1. $A$ is invertible and $A^{-1} \ge 0$ (element-wise).
    2. All principal minors of $A$ are positive (Hawkins-Simon conditions).
    3. All eigenvalues of $A$ have positive real parts.
    4. There exists a vector $x > 0$ such that $Ax > 0$.

---

## Exercises

1. **[Fundamentals] Is $A = \begin{pmatrix} 2 & -1 \ -1 & 2 \end{pmatrix}$ an M-matrix?**
   ??? success "Solution"
       Yes. It is a Z-matrix (off-diagonals $\le 0$). Principal minors: $2 > 0$ and $\det A = 3 > 0$. Since all principal minors are positive, it is an M-matrix. Note $A^{-1} = \frac{1}{3} \begin{pmatrix} 2 & 1 \ 1 & 2 \end{pmatrix} \ge 0$.

2. **[Inverse Positivity] Prove that if $A = sI - B$ with $B \ge 0$ and $s > ho(B)$, then $A$ is an M-matrix.**
   ??? success "Solution"
       Using the Neumann series: $A^{-1} = (sI - B)^{-1} = \frac{1}{s}(I - B/s)^{-1} = \frac{1}{s} \sum_{k=0}^\infty (B/s)^k$. Since $B \ge 0$ and $s > 0$, every term in the sum is nonnegative. The series converges because $ho(B/s) < 1$. Thus $A^{-1} \ge 0$.

3. **[Stability] Relate M-matrices to Hurwitz stability.**
   ??? success "Solution"
       A Z-matrix $A$ is an M-matrix iff $-A$ is Hurwitz stable. This implies that dynamical systems of the form $\dot{x} = -Ax$ (where $A$ is an M-matrix) always converge to the origin.

4. **[Economics] How does the M-matrix property relate to the Leontief input-output model?**
   ??? success "Solution"
       In the model $(I-C)x = d$, the matrix $I-C$ must be an M-matrix to ensure that for any positive demand $d$, the required production $x = (I-C)^{-1}d$ is non-negative and finite.

5. **[Minors] Show that any principal submatrix of an M-matrix is also an M-matrix.**
   ??? success "Solution"
       Since all principal minors of an M-matrix are positive, all principal minors of a principal submatrix (which are a subset of the original minors) are also positive. Thus the submatrix satisfies the M-matrix definition.

6. **[Diagonal Dominance] Prove that a Z-matrix with strictly positive diagonal dominance is an M-matrix.**
   ??? success "Solution"
       If $a_{ii} > \sum_{j 
eq i} |a_{ij}|$, then by Gershgorin's Circle Theorem, all eigenvalues have positive real parts. A Z-matrix with eigenvalues in the right half-plane is an M-matrix.

7. **[Monotonicity] Prove the Comparison Theorem: If $A$ is an M-matrix and $B \ge A$ is a Z-matrix, then $B$ is an M-matrix.**
   ??? success "Solution"
       $B \ge A \implies I - B \le I - A$ (roughly). More rigorously, since $A$ is an M-matrix, there exists $x > 0$ such that $Ax > 0$. Then $Bx \ge Ax > 0$. By the $Ax > 0$ characterization, $B$ is an M-matrix.

8. **[Singular M-matrices] Define a singular M-matrix and its property.**
   ??? success "Solution"
       A singular M-matrix $A$ is a Z-matrix where $\operatorname{Re}(\lambda_i) \ge 0$ and $ho(B) = s$. It arises in Markov chains where $I-P$ is a singular M-matrix (each column/row sum is zero).

9. **[Hadamard] Show that the Hadamard product of two M-matrices is not necessarily an M-matrix.**
   ??? success "Solution"
       While $A, B$ are M-matrices, $A \circ B$ is a Z-matrix, but its inverse might not be non-negative. However, if $A, B$ are M-matrices, then $A \circ B^{-1}$ and similar forms have specific positivity properties.

10. **[Markov Chains] Why is $L = I - P$ (where $P$ is a transition matrix) related to M-matrices?**
    ??? success "Solution"
        $L$ is a Z-matrix with zero row sums. It is a singular M-matrix. Its properties (like the existence of a positive null-vector) are fundamental to the theory of stationary distributions.

## Chapter Summary

This chapter establishes the theory of matrices with non-negative inverses:

1. **Sign Pattern Analysis**: Defined M-matrices as Z-matrices whose diagonal dominance ensures positivity.
2. **Equivalent Characterizations**: Detailed the multiple paths to verifying the M-matrix property, from minors to eigenvalues.
3. **Productive Balance**: Explored the role of M-matrices in guaranteeing non-negative solutions in economics and probability.
4. **Stability Link**: Demonstrated how the M-matrix structure provides a robust guarantee for the stability of positive systems.
