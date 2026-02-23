# Chapter 22: Numerical Linear Algebra

<div class="context-flow" markdown>

**Prerequisites**: Matrix Factorization (Ch10) · Norms and Perturbations (Ch15) · Matrix Analysis (Ch14)

**Chapter Outline**: From Theoretical to Numerical Solutions → Error Analysis Basics: Absolute vs. Relative Error → Floating Point Operations (FLOPs) and Rounding Errors → Backward Error Analysis → Direct Methods for Linear Systems (LU with Pivoting, QR) → Eigenvalue Algorithms: Power Method, Inverse Iteration, and QR Iteration → Sparse Matrix Techniques → Iterative Method Foundations: Jacobi, Gauss-Seidel, and SOR → Applications: Large-scale Simulations on Supercomputers and LAPACK Design Principles

**Extension**: Numerical Linear Algebra is "linear algebra that works"; it studies not the mathematical existence of $Ax=b$, but how to obtain the closest truth using finite-precision computers in $O(n^3)$ time or less—the engine behind all modern simulation software.

</div>

In theoretical mathematics, we assume real numbers have infinite precision. In a computer, however, all numbers are represented by finite bits (e.g., 64 bits). **Numerical Linear Algebra** specializes in designing algorithms that are both efficient and **robust** under the constraints of finite precision. A good numerical algorithm must not only be fast but also reliable when faced with ill-conditioned inputs. This chapter introduces the underlying algorithms and error analysis techniques that build the foundation of modern scientific computing.

---

## 22.1 Error Analysis and FLOPs

!!! definition "Definition 22.1 (Relative Error)"
    If $\hat{x}$ is an approximation and $x$ is the true value, the relative error is $\frac{\|\hat{x} - x\|}{\|x\|}$.
    **Backward Error Analysis**: Instead of asking how much the result deviates from the truth, we ask for which "nearby problem" our result is an exact solution. If this "nearby problem" is extremely close to the original, the algorithm is stable.

---

## 22.2 Eigenvalue Computation: QR Iteration

!!! theorem "Theorem 22.1 (Limits of Eigenvalue Algorithms)"
    The Abel-Ruffini theorem proves there is no general formula for roots of polynomials of degree 5 or higher. Consequently, all eigenvalue algorithms for $n \ge 5$ are essentially **iterative**.
    The **QR Iteration** is the most versatile eigenvalue algorithm:
    1.  $A_k = Q_k R_k$
    2.  $A_{k+1} = R_k Q_k$
    As the iteration proceeds, $A_k$ typically converges to an upper triangular form with eigenvalues on the diagonal.

---

## 22.3 Iterative Solvers for Linear Systems

!!! technique "Core Iterative Methods"
    1.  **Jacobi Iteration**: Based on the split $x_{k+1} = D^{-1}(b - (L+U)x_k)$.
    2.  **Gauss-Seidel Iteration**: Uses the latest computed components, typically converging faster than Jacobi.
    3.  **Successive Over-Relaxation (SOR)**: Accelerates convergence by introducing a relaxation parameter $\omega$.

---

## Exercises

**1. [Basics] What are FLOPs? How many FLOPs are required for an $n$-dimensional dot product?**

??? success "Solution"
    **Definition and Calculation:**
    1. **FLOPs**: Floating Point Operations.
    2. For a dot product $\sum a_i b_i$:
       - Requires $n$ multiplications.
       - Requires $n-1$ additions.
    **Conclusion**: Approximately $2n - 1$ operations, denoted as $O(n)$.

**2. [Complexity] Explain why an $O(n^3)$ algorithm is unacceptable for $n=10^6$.**

??? success "Solution"
    **Scale Estimation:**
    1. For $n = 10^6$, $n^3 = 10^{18}$ operations.
    2. Assume a computer performing $10^9$ operations per second (1 GFLOPS).
    3. Time required: $10^{18} / 10^9 = 10^9$ seconds.
    4. $10^9$ seconds is approximately 31.7 years.
    **Conclusion**: This illustrates why numerical algebra must strive for $O(n \log n)$ or more efficient algorithms for large-scale data.

**3. [Stability] Why is "Pivoting" necessary in LU decomposition?**

??? success "Solution"
    **Reasoning:**
    If a pivot $a_{ii}$ is close to zero, calculating the multiplier $l_{ji} = a_{ji}/a_{ii}$ results in a massive number. In subsequent subtractions, small rounding errors are multiplied by this huge factor, completely destroying the accuracy of the $U$ block. Pivoting (row swapping) ensures the denominator is the largest possible element in the column, keeping error amplification to a minimum.

**4. [QR Iteration] If $A = QR$, why does $RQ$ have the same eigenvalues as $A$?**

??? success "Solution"
    **Proof:**
    1. Given $A = QR$.
    2. Since $Q$ is orthogonal (unitary), $R = Q^* A$.
    3. $RQ = (Q^* A) Q = Q^* A Q$.
    **Conclusion**: $RQ$ is **unitarily similar** to $A$. Similarity transformations preserve eigenvalues.

**5. [Convergence] What is the necessary and sufficient condition for Jacobi iteration to converge?**

??? success "Solution"
    **Conclusion: The spectral radius $\rho(D^{-1}(L+U)) < 1$.**
    In practice, if the matrix $A$ is **strictly diagonally dominant**, Jacobi iteration is guaranteed to converge.

**6. [Sparse Matrices] How much memory is needed to store a $10^6 \times 10^6$ tridiagonal matrix?**

??? success "Solution"
    **Calculation:**
    1. A tridiagonal matrix has entries only on the main diagonal and the two adjacent diagonals.
    2. Number of non-zero elements $\approx 3n = 3 \times 10^6$.
    3. Assuming 64-bit doubles (8 bytes): $3 \times 10^6 \times 8 \approx 24$ MB.
    **Contrast**: Storing it as a full matrix requires $(10^6)^2 \times 8 = 8$ TB.
    **Conclusion**: Sparse storage is vital for large-scale computation.

**7. [Backward Error] What is "Backward Stability"?**

??? success "Solution"
    An algorithm is backward stable if the solution $\hat{y}$ it produces for input $x$ is the exact solution to a "nearby" problem with input $\hat{x}$ (where $\hat{x}$ is extremely close to $x$). This means the error produced by the algorithm is no worse than the noise already present in the input data.

**8. [Power Method] Briefly state the purpose and limitation of the Power Method.**

??? success "Solution"
    **Purpose**: To find the **dominant eigenvalue** (maximum magnitude) of a matrix.
    **Limitation**:
    1. It only finds the dominant eigenvalue.
    2. Convergence speed depends on the ratio $|\lambda_2 / \lambda_1|$; if they are close, it converges very slowly.

**9. [Conditioning] If $\kappa(A) = 10^{12}$, how many digits of the solution $x$ are likely correct in 16-digit double precision?**

??? success "Solution"
    **Conclusion: Approximately 4 digits.**
    **Rule of Thumb**: Significant digits $\approx$ Machine precision - $\log_{10} \kappa(A)$.
    Calculation: $16 - 12 = 4$. This indicates the problem is highly sensitive.

**10. [Application] Describe the relationship between LAPACK and BLAS.**

??? success "Solution"
    **Hierarchical Relationship:**
    1. **BLAS** (Basic Linear Algebra Subprograms): The low-level core responsible for simple vector addition and matrix multiplication, optimized for specific hardware architectures.
    2. **LAPACK** (Linear Algebra Package): The higher-level library that calls BLAS primitives to implement complex factorizations like LU, Eigen-decomposition, and SVD.
    This layering ensures numerical software achieves high-level abstraction without sacrificing performance.

## Chapter Summary

Numerical linear algebra is the "last mile" in turning mathematical truth into physical simulation:

1.  **Dominance of Stability**: It reveals that the merit of an algorithm lies not in the elegance of its form, but in its ability to suppress rounding errors, establishing backward error analysis as a scientific standard.
2.  **Efficiency Limits**: Through precise control of FLOPs and memory layouts, numerical algebra breaks the $O(n^3)$ barrier, making it possible to handle million-dimensional linear systems.
3.  **Philosophy of Iteration**: Facing the non-solvability of higher-order equations, iterative methods bridge the gap between finite steps and infinite precision, serving as the ultimate means for describing continuous dynamical processes.
