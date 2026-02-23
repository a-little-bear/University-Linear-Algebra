# Chapter 12: Jordan Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Decompositions (Ch10) · Polynomial Algebra (Ch00) · Spectral Theory (Ch13)

**Chapter Outline**: Limitations of Diagonalization (Defective Matrices) → Definition of Jordan Blocks → Generalized Eigenvectors and Jordan Chains → Root Subspace Decomposition → Existence and Uniqueness of Jordan Canonical Form (JCF) → Relationship between Minimal Polynomials and JCF → Steps to Determine JCF (Rank Methods, Weyr Characteristic) → Jordan Analysis of Matrix Powers and Series → Numerical Instability

**Extension**: The Jordan Canonical Form is the ultimate representation under similarity transformations; it perfectly reveals the "quasi-diagonal" structure of non-diagonalizable linear operators and is the necessary path for the theoretical derivation of systems of linear differential equations (Ch26) and matrix functions (Ch13).

</div>

Not all square matrices can be diagonalized. When the geometric multiplicity of an eigenvalue is less than its algebraic multiplicity, the matrix is called "defective." The **Jordan Canonical Form** (JCF) provides the structure closest to diagonal for such matrices. By introducing the "1" step structure (Jordan blocks), it quantifies the degree of space degeneracy as the length of algebraic chains. This chapter dives into the final verdict on matrix structure: the JCF.

---

## 12.1 Jordan Blocks and Generalized Eigenvectors

!!! definition "Definition 12.1 (Jordan Block)"
    A **Jordan block** $J_k(\lambda)$ of order $k$ is a square matrix with $\lambda$ on the main diagonal, 1s on the super-diagonal, and 0s elsewhere:
    $$J_k(\lambda) = \begin{pmatrix} \lambda & 1 & 0 & \cdots & 0 \\ 0 & \lambda & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & 0 & \lambda & 1 \\ 0 & 0 & 0 & 0 & \lambda \end{pmatrix}$$

!!! definition "Definition 12.2 (Generalized Eigenvectors and Jordan Chains)"
    A vector $\mathbf{v}_k$ is a **generalized eigenvector of rank $k$** corresponding to $\lambda$ if $(A - \lambda I)^k \mathbf{v}_k = \mathbf{0}$ but $(A - \lambda I)^{k-1} \mathbf{v}_k \neq \mathbf{0}$.
    The sequence $\{ \mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k \}$ forms a **Jordan chain**, where $\mathbf{v}_1$ is a standard eigenvector.

---

## 12.2 The Jordan Canonical Form Theorem

!!! theorem "Theorem 12.1 (JCF Existence and Uniqueness)"
    Every complex square matrix $A$ is similar to a Jordan Canonical Form $J$, and $J$ is unique up to the permutation of its blocks:
    $$P^{-1} A P = J = \operatorname{diag}(J_{k_1}(\lambda_1), J_{k_2}(\lambda_2), \ldots, J_{k_m}(\lambda_m))$$
    - Each Jordan block corresponds to one linearly independent eigenvector (geometric multiplicity).
    - The sum of the sizes of all Jordan blocks for a given eigenvalue equals its algebraic multiplicity.

---

## 12.3 Minimal Polynomials and JCF

!!! theorem "Theorem 12.2 (Minimal Polynomial Criterion)"
    The **minimal polynomial** $m(\lambda)$ is the monic polynomial of lowest degree such that $m(A) = O$.
    - **Diagonalization**: $A$ is diagonalizable $\iff$ $m(\lambda)$ has no repeated roots.
    - **Block Size**: The multiplicity of $\lambda_i$ in $m(\lambda)$ equals the **size of the largest Jordan block** for that eigenvalue in the JCF.

---

## Exercises

**1. [Jordan Block] Find the square of $J_2(5)$ and its eigenvalues.**

??? success "Solution"
    **Steps:**
    1. $J_2(5) = \begin{pmatrix} 5 & 1 \\ 0 & 5 \end{pmatrix}$.
    2. Multiplication: $\begin{pmatrix} 5 & 1 \\ 0 & 5 \end{pmatrix} \begin{pmatrix} 5 & 1 \\ 0 & 5 \end{pmatrix} = \begin{pmatrix} 25 & 10 \\ 0 & 25 \end{pmatrix}$.
    **Eigenvalue Analysis**:
    Since the result is upper triangular, the eigenvalue is $25$ with algebraic multiplicity 2. Note the result is still a modified Jordan block (non-diagonal).

**2. [Diagonalization] Determine if $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ is diagonalizable.**

??? success "Solution"
    **Logic:**
    1. Eigenvalue is 2 with algebraic multiplicity $\alpha = 2$.
    2. Calculate geometric multiplicity $\gamma = \dim E_2$:
       $A - 2I = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.
       The rank is 1, so $\gamma = 2 - 1 = 1$.
    3. Since $\gamma < \alpha$ (insufficient eigenvectors), the matrix is **not diagonalizable**. It is actually a $2 \times 2$ Jordan block.

**3. [JCF Determination] If a $3 \times 3$ matrix $A$ has all eigenvalues equal to 0 and $\operatorname{rank}(A)=1$, what is its JCF?**

??? success "Solution"
    **Process:**
    1. Total algebraic multiplicity is 3.
    2. Geometric multiplicity $\gamma = n - \operatorname{rank}(A) = 3 - 1 = 2$.
    3. $\gamma = 2$ means there are **2 Jordan blocks** in the JCF.
    4. Partition the total size 3 into 2 positive integers: only $2 + 1$ is possible.
    **Conclusion**: The JCF is $\operatorname{diag}(J_2(0), J_1(0))$.

**4. [Minimal Polynomial] Given $J = \operatorname{diag}(J_3(2), J_2(2))$, find its characteristic and minimal polynomials.**

??? success "Solution"
    **Calculation:**
    1. **Characteristic polynomial** $p(\lambda)$: Sum of block sizes. $(\lambda - 2)^{3+2} = (\lambda - 2)^5$.
    2. **Minimal polynomial** $m(\lambda)$: Size of the largest block. The largest block for $\lambda=2$ has size 3.
    **Conclusion**: $m(\lambda) = (\lambda - 2)^3$.

**5. [JCF Possibility] If the characteristic polynomial is $(\lambda-1)^4$ and the minimal polynomial is $(\lambda-1)^2$, find all possible JCFs.**

??? success "Solution"
    **Constraints:**
    1. Total size of blocks is 4.
    2. Maximum block size must be exactly 2.
    **Possible Combinations:**
    - Option A: $2 + 2$. JCF is $\operatorname{diag}(J_2(1), J_2(1))$.
    - Option B: $2 + 1 + 1$. JCF is $\operatorname{diag}(J_2(1), J_1(1), J_1(1))$.
    **Conclusion**: There are two non-similar possible structures.

**6. [Nilpotency] Describe $J_k(0)^n$ for $n \ge k$.**

??? success "Solution"
    **Conclusion:**
    The result is the **zero matrix**.
    **Reasoning**: $J_k(0)$ is strictly upper triangular. Each power shifts the super-diagonal 1s further to the top-right. After $k$ steps, all 1s fall off the boundary. This proves that Jordan blocks with zero eigenvalues are nilpotent.

**7. [Rank Method] How do you determine the number of blocks of size 1 for eigenvalue $\lambda$?**

??? success "Solution"
    **Derivation:**
    Let $n_k$ be the number of blocks of size $k$. Based on rank properties:
    $n_1 = \operatorname{rank}(A-\lambda I)^2 - 2\operatorname{rank}(A-\lambda I) + \operatorname{rank}(A-\lambda I)^0$.
    More generally, the second-order difference of the ranks of powers of $(A-\lambda I)$ yields the count of each block size.

**8. [Spaces] What is the difference between a generalized eigenspace and a standard eigenspace?**

??? success "Solution"
    **Contrast:**
    - **Eigenspace**: Set of vectors scaled by the operator ($(A-\lambda I)v=0$).
    - **Generalized Eigenspace**: Set of vectors eventually annihilated by repeated application of the shifted operator ($(A-\lambda I)^k v=0$).
    In JCF theory, the generalized eigenspace contains entire Jordan chains, filling the dimensional gap in non-diagonalizable matrices.

**9. [Uniqueness] If $A$ and $B$ have the same JCF, are they similar?**

??? success "Solution"
    **Conclusion:**
    **Yes**.
    **Reasoning**: JCF is a complete invariant under similarity. If two matrices are similar to the same canonical form, they are similar to each other by transitivity. JCF provides a perfect classification of matrix similarity classes.

**10. [Numerical] Why is JCF rarely computed directly in numerical software?**

??? success "Solution"
    **Stability Analysis:**
    1. JCF is extremely sensitive to perturbations.
    2. Example: $\begin{pmatrix} 0 & 1 \\ \epsilon & 0 \end{pmatrix}$ has distinct eigenvalues $\pm\sqrt{\epsilon}$ and is diagonalizable for any $\epsilon \neq 0$.
    3. At $\epsilon = 0$, it suddenly becomes a single $2 \times 2$ Jordan block.
    **Conclusion**: This discontinuity makes JCF unstable under floating-point arithmetic. **Schur decomposition** or **SVD** are typically used instead.

## Chapter Summary

The Jordan Canonical Form is the final verdict on the structure of square matrices:

1.  **Completing the Defective**: It fills the gap in eigenvectors via Jordan chains, establishing the root subspace structure of the space.
2.  **Polynomial Depth**: The correspondence between the minimal polynomial and JCF block sizes reveals the geometric depth of a matrix as a root of a polynomial.
3.  **Structural Uniqueness**: JCF establishes the classification standard for similarity classes, providing the most precise theoretical framework for analyzing matrix functions and power series.
