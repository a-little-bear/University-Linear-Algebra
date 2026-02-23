# Chapter 38A: M-matrices and Z-matrices

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices (Ch17) · Matrix Analysis (Ch14) · Positive Definite Matrices (Ch16)

**Chapter Outline**: Definition of Z-matrices (Non-positive Off-diagonals) → Definition of M-matrices → 10+ Equivalent Characterizations (Inverse-positivity, Principal Minors, Positive Vector Criterion) → Matrix Splitting & Iterative Convergence → Varga Comparison Theorem → Singular M-matrices → Applications: Numerical Solution of PDEs (Finite Difference Matrices), Leontief Economic Models, Markov Chain Steady States

**Extension**: M-matrices bridge the gap between non-negative and positive definite matrices; they provide the algebraic criterion for whether numerical schemes (like finite differences) satisfy the "Maximum Principle" and are core conditions for guaranteeing the absolute convergence of iterative methods (like Gauss-Seidel).

</div>

In numerical analysis and economic modeling, a special class of matrices arises: their off-diagonal entries are all non-positive, yet they possess exceptional positivity properties (such as having a strictly positive inverse). These are known as **M-matrices**. They are not only the bedrock for determining the convergence of iterative solvers but also the algebraic abstraction of diffusion phenomena in physics and input-output balance in economics.

---

## 38A.1 Definitions of Z-matrices and M-matrices

!!! definition "Definition 38A.1 (Z-matrix)"
    A matrix $A$ is a **Z-matrix** if all its off-diagonal entries are non-positive:
    $$a_{ij} \le 0, \quad \forall i \neq j$$

!!! definition "Definition 38A.2 (M-matrix)"
    A Z-matrix $A$ is an **M-matrix** if it can be expressed as:
    $$A = sI - B, \quad B \ge 0, \quad s \ge \rho(B)$$
    where $\rho(B)$ is the spectral radius of $B$. If $s > \rho(B)$, $A$ is a **non-singular M-matrix**.

---

## 38A.2 Core Characterizations

!!! theorem "Theorem 38A.1 (Equivalent Conditions for Non-singular M-matrices)"
    For a Z-matrix $A$, the following are equivalent to $A$ being a non-singular M-matrix:
    1.  **Inverse-positivity**: $A$ is invertible and $A^{-1} \ge 0$.
    2.  **Principal Minors**: All leading principal minors of $A$ are positive.
    3.  **Positive Vector**: There exists a positive vector $\mathbf{x} > 0$ such that $A\mathbf{x} > 0$.
    4.  **Spectral Property**: The real parts of all eigenvalues of $A$ are positive.

---

## 38A.3 Iterative Convergence and Comparison

!!! theorem "Theorem 38A.2 (Convergence of Matrix Splittings)"
    Let $A = M - N$ be a regular splitting of $A$ ($M$ is invertible, $M^{-1} \ge 0$, and $N \ge 0$). If $A$ is an M-matrix, then the spectral radius $\rho(M^{-1}N) < 1$.
    **Significance**: This guarantees that for systems with M-matrices, Jacobi and Gauss-Seidel iterations always converge.

---

## Exercises

1.  **[Criteria] Determine if $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$ is an M-matrix.**
    ??? success "Solution"
        It is a Z-matrix, and its principal minors $D_1=2, D_2=3$ are both positive. Thus, it is a non-singular M-matrix.

2.  **[Inverse] Prove: If $A$ is a non-singular M-matrix, then $Ax=b$ with $b \ge 0$ implies $x \ge 0$.**
    ??? success "Solution"
        Since $x = A^{-1}b$ and $A^{-1} \ge 0$ with $b \ge 0$, the property of non-negative matrices ensures $x \ge 0$.

3.  **[Dominance] Prove that a Z-matrix with positive diagonals and strict diagonal dominance is an M-matrix.**
    ??? success "Solution"
        Take $\mathbf{x} = \mathbf{1}$. Diagonal dominance implies $A\mathbf{1} > 0$, satisfying the positive vector criterion.

4.  **[Diagonals] Can the diagonal entries of an M-matrix be negative?**
    ??? success "Solution"
        No. If $A$ is an M-matrix, the principal minor criterion requires $a_{ii} > 0$ (or, from $s > \rho(B)$, $s$ must be greater than any diagonal entry of $B$).

5.  **[Economics] In the Leontief model $(I-A)x = d$, why is $I-A$ usually an M-matrix?**
    ??? success "Solution"
        The consumption matrix $A$ is non-negative ($A \ge 0$), and since the economic system must produce a surplus, its spectral radius $\rho(A) < 1$. This matches the definition of an M-matrix.

6.  **[Determinant] Prove: If $A$ is an M-matrix, its determinant $\det(A) > 0$.**
    ??? success "Solution"
        By the principal minor criterion, all leading principal minors are positive, including $\det(A)$.

7.  **[Comparison] If $A$ and $B$ are Z-matrices and $A \le B$ entry-wise, if $A$ is an M-matrix, is $B$ also an M-matrix?**
    ??? success "Solution"
        Yes. This is the Comparison Theorem: increasing diagonal entries or decreasing the absolute value of off-diagonal entries strengthens the M-matrix property.

8.  **[Singular] Give an example of a singular M-matrix.**
    ??? success "Solution"
        $A = \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$. Here $s = \rho(B) = 1$.

9.  **[Jacobi] Prove the spectral radius of the Jacobi iteration matrix for a diagonally dominant M-matrix is $\rho(J) < 1$.**
    ??? success "Solution"
        Using $\rho(A) \le \|A\|_\infty$, for $J = D^{-1}(L+U)$, the row sums are $\sum_{j \neq i} |a_{ij}|/a_{ii} < 1$ due to dominance.

10. **[PDEs] Is the matrix obtained from numerically solving the Laplace equation $\Delta u = f$ an M-matrix?**

   ??? success "Solution"
        Yes. The five-point stencil matrix obtained via finite differences is a typical strictly diagonally dominant M-matrix, ensuring the numerical solution satisfies the Maximum Principle (no spurious oscillations).

## Chapter Summary

M-matrices represent the intersection of numerical stability and physical reality:

1.  **Inverse Positivity**: The most profound trait of M-matrices is inverse-positivity, ensuring that positive inputs (stimuli) lead to positive outputs (responses)—a fundamental requirement for models to remain physically logical.
2.  **Anchor of Convergence**: In handling large-scale linear systems, the M-matrix structure is the ultimate defense against divergence in iterative algorithms (such as discrete PDE solvers).
3.  **Structural Dominance**: Through comparison theorems, M-matrices provide a powerful lever for estimating operator norms and eigenvalue ranges, transforming complex operator comparisons into simple entry-wise checks.
