# Chapter 22: Numerical Linear Algebra

<div class="context-flow" markdown>

**Prerequisites**: Matrix Decompositions (Ch10) · Norms and Perturbation (Ch15) · Systems of Equations (Ch1)

**Chapter Outline**: Rounding Errors and Floating-Point Arithmetic → Computational Complexity ($O(n^3)$) → Stability of Direct Methods (LU, QR, Cholesky) → Iterative Methods (Jacobi, Gauss-Seidel, SOR) → Krylov Subspace Methods (GMRES, Conjugate Gradient CG) → Preconditioning Techniques → Sparse Matrix Computation

**Extension**: Numerical linear algebra is the underlying engine for all high-performance scientific computing (weather forecasting, finite element analysis, neural network training).

</div>

Numerical linear algebra studies how to efficiently and stably solve linear algebra problems on computers with finite precision. It focuses not on theoretical solutions, but on **execution time**, **memory efficiency**, and **error control**.

---

## 22.1 Stability and Complexity

!!! definition "Definition 22.1 (Numerical Stability)"
    An algorithm is **backward stable** if, in the presence of rounding errors, the solution it produces is exactly the solution to a slightly perturbed problem.

!!! theorem "Theorem 22.3 (Convergence of Conjugate Gradient)"
    For a symmetric positive definite matrix $A$, the CG method converges in at most $n$ steps (ignoring rounding errors). Its convergence rate depends on the condition number $\sqrt{\kappa(A)}$.

---

## Exercises

1. **[Complexity] What is the theoretical time complexity of multiplying two $n$-th order matrices? How does the Strassen algorithm optimize it?**
   ??? success "Solution"
       Standard multiplication is $O(n^3)$. Strassen's algorithm reduces the number of multiplications through clever partitioning, lowering complexity to approximately $O(n^{2.81})$. For large-scale applications, even lower-complexity algorithms (like Coppersmith-Winograd) exist.

2. **[Rounding Error] Why does calculating $1.000001 - 1.000000$ lead to precision loss?**
   ??? success "Solution"
       This is **Catastrophic Cancellation**. Subtracting two close numbers cancels the high-order significant digits, drastically reducing the number of valid digits in the result and causing relative error to skyrocket.

3. **[Iterative Methods] State a sufficient condition for the Jacobi method to converge.**
   ??? success "Solution"
       The matrix $A$ is **strictly diagonally dominant** (i.e., $|a_{ii}| > \sum_{j 
eq i} |a_{ij}|$). In this case, the spectral radius of the iteration matrix satisfies $ho(B) < 1$.

4. **[Preconditioning] Why is the preconditioner $M$ chosen to be an approximation of $A$?**
   ??? success "Solution"
       Preconditioning aims to solve $M^{-1}Ax = M^{-1}b$. If $M \approx A$, then $M^{-1}A \approx I$. Since the identity matrix has a condition number of 1 and concentrated eigenvalues, this greatly speeds up the convergence of Krylov methods like CG or GMRES.

5. **[QR Algorithm] Briefly describe the steps of the QR algorithm for computing all eigenvalues of a matrix.**
   ??? success "Solution"
       1. Set $A_0 = A$.
       2. Perform QR decomposition: $A_k = Q_k R_k$.
       3. Form the next generation: $A_{k+1} = R_k Q_k$.
       Repeat until $A_k$ converges to upper triangular form (the diagonal entries are the eigenvalues). Shifts are often used to accelerate convergence.

6. **[Calculation] Perform one Jacobi iteration on $A = \begin{pmatrix} 10 & 1 \ 1 & 10 \end{pmatrix}$ with $x_0=0$ and $b=(11, 11)^T$.**
   ??? success "Solution"
       $x_1^{(1)} = (11 - 1(0))/10 = 1.1$.
       $x_1^{(2)} = (11 - 1(0))/10 = 1.1$.
       Thus $x_1 = (1.1, 1.1)^T$.

7. **[Krylov Subspace] Define the $n$-th order Krylov subspace $\mathcal{K}_n(A, b)$.**
   ??? success "Solution"
       $\mathcal{K}_n(A, b) = \operatorname{span}\{b, Ab, A^2b, \dots, A^{n-1}b\}$. This is the basis where iterative methods (like Lanczos or Arnoldi) search for the solution.

8. **[Sparse Matrices] For large sparse matrices, why is computing $A^{-1}$ directly avoided?**
   ??? success "Solution"
       The inverse of a sparse matrix is typically **dense** (Fill-in phenomenon). Storing $A^{-1}$ would cause memory overflow, and the computational cost is far higher than using direct or iterative solvers.

9. **[Householder] What is the advantage of Householder transformations over Gram-Schmidt in QR decomposition?**
   ??? success "Solution"
       Householder transformations possess superior **numerical orthogonality**. While Modified Gram-Schmidt (MGS) is better than standard GS, it can still lose orthogonality in ill-conditioned cases, whereas Householder is unconditionally backward stable.

10. **[Application] Why is Stochastic Gradient Descent (SGD) preferred over Newton's method in deep learning?**
    ??? success "Solution"
        Newton's method requires storing and inverting the Hessian matrix, which has space complexity $O(n^2)$ and time complexity $O(n^3)$. For millions of dimensions, this cost is unacceptable. SGD only requires $O(n)$ complexity.

## Chapter Summary

Numerical linear algebra is the physical boundary of computation:

1. **Error is a First-Class Citizen**: The quality of an algorithm is determined not by analytical correctness, but by the rate of error accumulation.
2. **Structured Computation**: Specialized algorithms for sparse, symmetric, or banded structures are the cornerstones of high-performance computing.
3. **The Art of Iteration**: In ultra-large systems, replacing global exact calculations with finite local updates is the only way forward.
