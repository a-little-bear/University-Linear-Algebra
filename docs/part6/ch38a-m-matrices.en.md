# Chapter 38A: M-matrices and Z-matrices

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices & Perron-Frobenius (Ch17) · Matrix Analysis (Ch14) · Positive Definite Matrices (Ch16)

**Chapter Outline**: Definition of Z-matrices (Non-positive Off-diagonals) → Various Definitions of M-matrices → Core Equivalent Characterizations (Inverse-positivity, Principal Minors, Positive Vector) → Matrix Splitting Theory → The Role of M-matrices in Iterative Convergence → Comparison Theorems & Norm Bounds → Singular M-matrices → Applications: Stability of Discrete PDEs, Leontief Models in Economics, and Ecosystem Balance Analysis

**Extension**: M-matrices bridge the gap between non-negative and positive definite matrices; they provide the algebraic criterion for whether numerical schemes (like finite difference methods) satisfy the "Maximum Principle," and are core conditions for guaranteeing the absolute convergence of iterative solvers for large systems.

</div>

In numerical analysis and economic modeling, a class of matrices frequently appears with non-positive off-diagonal entries but excellent overall properties. **M-matrices** (named after Minkowski) are the algebraic abstraction of these structures. Their most significant trait is having a non-negative inverse, which ensures that physical systems respond positively to positive stimuli. This chapter explores over 10 equivalent criteria for M-matrices and their foundational role in computational science.

---

## 38A.1 Definitions of Z-matrices and M-matrices

!!! definition "Definition 38A.1 (Z-matrix)"
    A square matrix $A$ is a **Z-matrix** if all its off-diagonal entries are non-positive: $a_{ij} \le 0$ for all $i 
eq j$.

!!! definition "Definition 38A.2 (M-matrix)"
    A Z-matrix $A$ is an **M-matrix** if it can be expressed as:
    $$A = sI - B, \quad B \ge 0, \quad s \ge ho(B)$$
    where $ho(B)$ is the spectral radius of the non-negative matrix $B$. If $s > ho(B)$, $A$ is a **non-singular M-matrix**.

---

## 38A.2 Core Characterizations

!!! theorem "Theorem 38A.1 (Identifying Non-singular M-matrices)"
    For a Z-matrix $A$, the following conditions are equivalent:
    1.  **Inverse-positivity**: $A$ is invertible and $A^{-1} \ge 0$ (entry-wise non-negative).
    2.  **Eigenvalues**: The real parts of all eigenvalues of $A$ are positive.
    3.  **Principal Minors**: All leading principal minors of $A$ are strictly positive.
    4.  **Positive Vector**: There exists a vector $\mathbf{x} > 0$ such that $A\mathbf{x} > 0$.

---

## 38A.3 Matrix Splitting and Iteration

!!! technique "Technique: Convergence Criterion"
    Let $A = M - N$. If $M$ is a non-singular M-matrix and $N \ge 0$, then the splitting is convergent (i.e., $ho(M^{-1}N) < 1$).
    **Significance**: This provides the most robust algebraic proof for the convergence of iterative methods like Jacobi and Gauss-Seidel for M-matrix systems.

---

## Exercises

**1. [Criteria] Determine if $A = \begin{pmatrix} 2 & -1 \ -1 & 2 \end{pmatrix}$ is an M-matrix.**

??? success "Solution"
    **Steps:**
    1. Check Z-matrix condition: Off-diagonal entries are -1, both $\le 0$. Yes.
    2. Calculate principal minors: $D_1 = 2 > 0, D_2 = 4 - 1 = 3 > 0$.
    3. Since all leading principal minors are positive and it is a Z-matrix.
    **Conclusion**: It is a non-singular M-matrix.

**2. [Inverse Positivity] If $A$ is a non-singular M-matrix, prove that for $Ax=b$, if $b \ge 0$, then $x \ge 0$.**

??? success "Solution"
    **Proof:**
    1. The solution is $x = A^{-1}b$.
    2. From M-matrix properties, $A^{-1} \ge 0$ (Inverse-positivity).
    3. Multiplying a non-negative matrix by a non-negative vector yields a non-negative result.
    **Conclusion**: $x \ge 0$. This physically guarantees that "positive input results in positive response."

**3. [Dominance] Prove that a Z-matrix with positive diagonal entries and strict diagonal dominance is an M-matrix.**

??? success "Solution"
    **Proof:**
    1. Let $\mathbf{x} = (1, 1, \ldots, 1)^T$.
    2. $(A\mathbf{x})_i = a_{ii} + \sum_{j 
eq i} a_{ij} = a_{ii} - \sum_{j 
eq i} |a_{ij}|$.
    3. Due to strict diagonal dominance, this value is greater than 0 for all $i$.
    4. Thus, there exists a positive vector $\mathbf{x} > 0$ such that $A\mathbf{x} > 0$.
    **Conclusion**: It satisfies the positive vector criterion and is therefore an M-matrix.

**4. [Spectral] Can the eigenvalues of an M-matrix lie on the imaginary axis?**

??? success "Solution"
    **Conclusion:**
    **No** (for the non-singular case).
    **Reasoning**: All eigenvalues of a non-singular M-matrix have real parts strictly greater than 0. This restricts the eigenvalues to the open right half-plane, ensuring the asymptotic stability of the system.

**5. [Economics] In the Leontief model $(I-A)x = d$, why is $I-A$ usually an M-matrix?**

??? success "Solution"
    **Explanation:**
    1. $A$ is the consumption matrix with non-negative entries.
    2. Realistic economic systems must be "profitable," implying $ho(A) < 1$.
    3. Following the definition with $s=1 > ho(A)$, $I-A$ fits the form of a non-singular M-matrix.
    This structure ensures the system produces a positive net output to satisfy external demand.

**6. [Comparison] If $A, B$ are Z-matrices and $A \le B$ (entry-wise), if $A$ is an M-matrix, prove $B$ is also one.**

??? success "Solution"
    **Logic:**
    1. Consider the positive vector $\mathbf{x}$ such that $A\mathbf{x} > 0$.
    2. $B\mathbf{x} = A\mathbf{x} + (B-A)\mathbf{x}$.
    3. Since $B-A \ge 0$ and $\mathbf{x} > 0$, $(B-A)\mathbf{x} \ge 0$.
    4. Thus $B\mathbf{x} \ge A\mathbf{x} > 0$.
    **Conclusion**: $B$ inherits the positive vector criterion and is also an M-matrix.

**7. [Properties] Prove: The diagonal entries of an M-matrix must be positive.**

??? success "Solution"
    **Proof:**
    By the principal minor criterion, $a_{ii}$ is a 1st-order leading principal minor and must be strictly positive. Physically, this represents that the "self-return" strength at each node must exceed the negative coupling from other nodes.

**8. [Singularity] Give an example of a singular M-matrix.**

??? success "Solution"
    **Example:** $A = \begin{pmatrix} 1 & -1 \ -1 & 1 \end{pmatrix}$.
    **Analysis**: It is a Z-matrix, and the spectral radius of the non-negative part $B = \begin{pmatrix} 0 & 1 \ 1 & 0 \end{pmatrix}$ is exactly $ho = 1$. Here $s=1=ho$, and the determinant is 0.

**9. [Iteration] In Jacobi iteration for $Ax=b$, what does it mean for $A$ to be an M-matrix?**

??? success "Solution"
    **Conclusion:**
    It means the Jacobi iteration matrix $B = D^{-1}(L+U)$ is non-negative. Its convergence is governed by $ho(B)$, and for M-matrices, this convergence is typically monotonic and analytically guaranteed.

**10. [Application] Why is it desirable for the discretization matrix of a diffusion equation to be an M-matrix?**

??? success "Solution"
    **Reasoning:**
    Diffusion processes follow the **Maximum Principle** (concentrations do not spontaneously create local extrema). If the discretized matrix is not an M-matrix, the numerical solution might produce "spurious oscillations" or negative values at steep gradients, violating physical laws. The M-matrix structure guarantees the **fidelity** of the numerical scheme.

## Chapter Summary

M-matrices represent the intersection of numerical stability and physical reality:

1.  **Inverse Positivity**: The most profound trait of M-matrices is inverse-positivity, ensuring that positive inputs lead to positive outputs—a fundamental requirement for linear systems to remain logical.
2.  **Anchor of Convergence**: In handling large-scale systems, the M-matrix structure is the ultimate defense against divergence in iterative algorithms.
3.  **Structural Dominance**: Through comparison theorems, M-matrices provide a powerful lever for estimating operator norms and eigenvalue ranges, transforming complex operator comparisons into simple entry-wise checks.
