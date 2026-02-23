# Chapter 22: Numerical Linear Algebra

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch01) · Matrix Decompositions (Ch10) · Norms and Perturbation (Ch15)

**Chapter Outline**: Floating-Point Arithmetic & Rounding Errors → Algorithmic Complexity → Condition Numbers vs. Numerical Stability → Forward and Backward Stability → Direct Solvers (LU, Cholesky, QR) → Iterative Solvers (Jacobi, Gauss-Seidel, SOR) → Krylov Subspace Methods (Conjugate Gradient, GMRES) → Preconditioning → Eigenvalue Algorithms (Power Method, QR Algorithm) → Sparse Matrix Techniques

**Extension**: Numerical Linear Algebra is the engine under the hood of almost all scientific simulations, weather forecasting, and modern deep learning frameworks; it bridges the gap between pure mathematical theory and finite-precision hardware.

</div>

Theoretical linear algebra assumes exact arithmetic over fields like $\mathbb{R}$ or $\mathbb{C}$. In practice, however, computers operate using finite-precision floating-point numbers, leading to inevitable rounding errors. **Numerical Linear Algebra** is the study of algorithms that are not only computationally efficient but also robust against these errors. A theoretically correct algorithm may fail catastrophically in practice if it is numerically unstable. This chapter explores the tradeoff between complexity, precision, and stability in large-scale computation.

---

## 22.1 Stability and Complexity

<div class="context-flow" markdown>

**The Fundamental Tradeoff**: Fast algorithms (like Strassen's) often trade numerical stability for speed. Large matrices require iterative solutions to avoid the $O(n^3)$ cost of direct methods.

</div>

!!! definition "Definition 22.1 (Backward Stability)"
    An algorithm $\tilde{f}(x)$ for a problem $f(x)$ is **backward stable** if for every $x$, $\tilde{f}(x) = f(x + \Delta x)$ for some small $\Delta x$.
    In other words, a backward stable algorithm gives the *exact* solution to a *slightly perturbed* problem.

!!! theorem "Theorem 22.1 (Complexity Benchmarks)"
    - **Matrix Multiplication**: $O(n^3)$ (standard), $O(n^{2.81})$ (Strassen).
    - **Gaussian Elimination / LU**: $O(n^3)$.
    - **Back-substitution**: $O(n^2)$.
    - **QR Decomposition**: $O(n^3)$ (but with a higher constant than LU).

---

## 22.2 Direct Solvers: LU, Cholesky, and QR

!!! technique "Pivoting in LU Decomposition"
    Gaussian elimination without pivoting is numerically unstable because dividing by a small pivot magnifies errors. **Partial Pivoting** ($PA = LU$) ensures that the multipliers are always $\le 1$, providing a stable (though not guaranteed backward stable) result for most matrices.

!!! theorem "Theorem 22.2 (Stability of Cholesky)"
    For symmetric positive definite (SPD) matrices, Cholesky decomposition $A = LL^T$ is **unconditionally backward stable**. It is also twice as fast as LU and does not require pivoting.

!!! theorem "Theorem 22.3 (Householder QR)"
    The Householder transformation-based QR decomposition is the gold standard for numerical stability in orthogonalization, outperforming the classical Gram-Schmidt process which suffers from loss of orthogonality due to rounding.

---

## 22.3 Iterative Solvers and Krylov Subspaces

<div class="context-flow" markdown>

**Motivation**: For extremely large or sparse matrices (e.g., $n > 10^6$), $O(n^3)$ direct methods are impossible. Iterative methods refine an initial guess $x_0$ to approach the true solution.

</div>

!!! definition "Definition 22.2 (Krylov Subspace)"
    The $k$-th Krylov subspace associated with $A$ and $b$ is:
    $$\mathcal{K}_k(A, b) = \operatorname{span}\{b, Ab, A^2b, \ldots, A^{k-1}b\}$$
    Modern solvers like **CG** and **GMRES** find the "best" approximation of $x$ within these subspaces.

!!! theorem "Theorem 22.4 (Conjugate Gradient Convergence)"
    For an SPD matrix $A$, the Conjugate Gradient (CG) method converges in at most $n$ steps (in exact arithmetic). The rate of convergence depends on the condition number $\kappa(A)$:
    $$\|x_k - x\|_A \le 2 \left( \frac{\sqrt{\kappa(A)} - 1}{\sqrt{\kappa(A)} + 1} \right)^k \|x_0 - x\|_A$$

---

## 22.4 Preconditioning

!!! technique "Preconditioning"
    To accelerate iterative solvers, we solve $M^{-1}Ax = M^{-1}b$ instead of $Ax=b$. If $M \approx A$, then $M^{-1}A \approx I$, which has a condition number near 1, leading to near-instant convergence. Common choices include Jacobi, Incomplete LU (ILU), and Sparse Approximate Inverse (SPAI).

---

## 22.5 Eigenvalue Algorithms

!!! algorithm "Algorithm 22.1 (The QR Algorithm)"
    The most widely used method for finding all eigenvalues of a dense matrix:
    1.  Start with $A_0 = A$.
    2.  Iterate: Factor $A_k = Q_k R_k$, then set $A_{k+1} = R_k Q_k$.
    3.  $A_k$ converges to an upper triangular matrix (the Schur form), with eigenvalues on the diagonal.
    *Note: Shifting techniques are used to accelerate convergence to cubic rates.*

---

## Exercises

1.  **[Complexity] Calculate the total number of floating-point operations (FLOPs) for a standard dot product of two $n$-vectors.**
    ??? success "Solution"
        $n$ multiplications and $n-1$ additions, total $2n-1 \approx 2n$ FLOPs.

2.  **[Stability] Why is $1.000001 - 1.000000$ dangerous in numerical computing?**
    ??? success "Solution"
        This is **catastrophic cancellation**. The result has only one significant digit, losing the precision of the original numbers. In matrix operations, this often happens during elimination without pivoting.

3.  **[Condition] Find the condition number of $\begin{pmatrix} 100 & 0 \\ 0 & 0.01 \end{pmatrix}$ relative to the 2-norm.**
    ??? success "Solution"
        $\sigma_{\max}=100, \sigma_{\min}=0.01 \implies \kappa = 100 / 0.01 = 10,000$.

4.  **[Iterative] Under what condition is the Jacobi iteration $x_{k+1} = D^{-1}(b - (L+U)x_k)$ guaranteed to converge?**
    ??? success "Solution"
        If $A$ is **strictly diagonally dominant**. This ensures the spectral radius of the iteration matrix is less than 1.

5.  **[Krylov] Define the Krylov subspace $\mathcal{K}_3(A, b)$.**
    ??? success "Solution"
        $\operatorname{span}\{b, Ab, A^2b\}$.

6.  **[Preconditioning] Why is the identity matrix $I$ a bad preconditioner?**
    ??? success "Solution"
        $I^{-1}A = A$. It does not change the condition number or the spectrum of the problem, thus providing no acceleration.

7.  **[QR Algorithm] Prove that $A_{k+1} = R_k Q_k$ is similar to $A_k$.**
    ??? success "Solution"
        $A_k = Q_k R_k \implies R_k = Q_k^{-1} A_k$.
        Substituting into $A_{k+1}$ gives $A_{k+1} = (Q_k^{-1} A_k) Q_k = Q_k^T A_k Q_k$.
        Since $A_{k+1}$ is a unitary similarity transform of $A_k$, they share the same eigenvalues.

8.  **[Sparse] What is "fill-in" in sparse LU decomposition?**
    ??? success "Solution"
        Fill-in occurs when an entry that was zero in the original matrix becomes non-zero during the elimination process. Minimizing fill-in via row reordering (e.g., Cuthill-McKee) is vital for memory efficiency.

9.  **[Power Method] If the two largest eigenvalues of $A$ are $\lambda_1 = 10$ and $\lambda_2 = 9.9$, why will the Power Method converge slowly?**
    ??? success "Solution"
        The convergence rate is governed by the ratio $|\lambda_2 / \lambda_1|$. Here $9.9/10 = 0.99$, which is very close to 1, meaning the error term vanishes very slowly.

****

??? success "Solution"
    

## Chapter Summary

Numerical Linear Algebra is the physics of mathematical computation:

1.  **Error Awareness**: Recognized that precision is a finite resource and that algorithms must be designed to contain, not amplify, rounding noise.
2.  **Structural Exploitation**: Showed how SPD, sparse, and banded structures allow for algorithms that bypass the $O(n^3)$ complexity wall.
3.  **Iterative Dominance**: Established Krylov subspace methods as the primary tools for modern high-dimensional engineering and data science.
4.  **Stability First**: Validated that backward stability is the essential requirement for any algorithm intended for real-world execution.
