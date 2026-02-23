# Chapter 22: Numerical Linear Algebra and Stability

<div class="context-flow" markdown>

**Prerequisites**: Gaussian Elimination (Ch1) · Condition Number (Ch15) · Matrix Decompositions (Ch10) · Floating Point Arithmetic

**Chapter Outline**: Precision and Rounding Errors → Forward and Backward Stability → Complexity Analysis of Matrix Operations → Condition of Problems vs. Stability of Algorithms → Direct Solvers (LU, Cholesky, QR) → Iterative Solvers (Conjugate Gradient, GMRES) → Preconditioning → Eigenvalue Algorithms (QR Algorithm, Lanczos)

**Extension**: Numerical linear algebra is the science of computing efficiently and accurately on real hardware; it is the engine of large-scale scientific simulation and machine learning.

</div>

Theoretical linear algebra assumes exact arithmetic, but real-world computation is limited by finite precision and rounding errors. **Numerical linear algebra** studies the design and analysis of algorithms that are efficient ($O(n^3)$ or better) and **backward stable** (meaning the computed solution is the exact solution to a slightly perturbed problem). This chapter explores the tradeoff between direct and iterative solvers and introduces the core algorithms used to handle matrices with millions of entries.

---

## 22.1 Floating Point and Algorithm Stability

!!! definition "Definition 22.1 (Backward Stability)"
    An algorithm $\tilde{f}(x)$ for a problem $f(x)$ is **backward stable** if for any $x$, $\tilde{f}(x) = f(x + \delta x)$ for some small $\delta x$ relative to $x$. This ensures that errors are no worse than those inherent in the input data.

!!! theorem "Theorem 22.1 (Conditioning vs. Stability)"
    The error in a computed solution $\tilde{y}$ satisfies:
    $$\frac{\|\tilde{y} - y\|}{\|y\|} \lesssim \kappa(\text{problem}) \times \text{error}(\text{algorithm})$$
    Total error is the product of the problem's sensitivity (condition number) and the algorithm's numerical noise.

---

## Exercises

1. **[Fundamentals] Compute the FLOP count for multiplying two $n \times n$ matrices.**
   ??? success "Solution"
       $n^3$ multiplications and $n^2(n-1)$ additions. The complexity is $O(n^3)$.

2. **[Stability] Why is Gaussian elimination without pivoting considered unstable?**
   ??? success "Solution"
       Because a small pivot $a_{ii}$ leads to massive multipliers $a_{ji}/a_{ii}$, which amplify rounding errors in subsequent rows. Partial pivoting ensures multipliers are $\le 1$.

3. **[Direct vs. Iterative] When is an iterative solver (like Conjugate Gradient) preferred over a direct solver (like LU)?**
   ??? success "Solution"
       For large, sparse matrices where $n$ is too large for $O(n^3)$ work, or where $A$ is not explicitly stored but available as a function $x \mapsto Ax$.

4. **[Calculation] Find the condition number of $\begin{pmatrix} 1 & 1 \\ 1 & 1.0001 \end{pmatrix}$.**
   ??? success "Solution"
       The matrix is near-singular ($\det \approx 0.0001$). $\kappa \approx 10^4$. Solving $Ax=b$ with this matrix will lose approximately 4 digits of precision.

5. **[Preconditioning] What is the goal of preconditioning in $Mx = b$?**
   ??? success "Solution"
       To solve $P^{-1}Ax = P^{-1}b$ instead, where $P \approx A$. A good preconditioner $P$ reduces the condition number, accelerating the convergence of iterative methods.

6. **[Complexity] Compare the complexity of solving $Ax=b$ via Cramer's Rule versus Gaussian Elimination.**
   ??? success "Solution"
       Cramer's Rule is $O(n \cdot n!)$, which is computationally impossible for $n > 20$. Gaussian Elimination is $O(n^3)$. This gap highlights the importance of algorithmic design.

7. **[Eigenvalues] Describe the QR Algorithm for finding eigenvalues.**
   ??? success "Solution"
       It iteratively factors $A_k = Q_k R_k$ and then sets $A_{k+1} = R_k Q_k$. Under mild conditions, $A_k$ converges to a triangular matrix whose diagonal contains the eigenvalues.

8. **[Fill-in] In sparse matrix solvers, what is "fill-in"?**
   ??? success "Solution"
       The creation of non-zero entries in positions that were originally zero during factorization. Minimizing fill-in (via reordering) is key to sparse LU and Cholesky efficiency.

9. **[Precision] How does the machine epsilon $\epsilon_{\text{mach}}$ limit accuracy?**
   ??? success "Solution"
       $\epsilon_{\text{mach}}$ is the smallest number such that $1 + \epsilon > 1$. It defines the base noise level of all calculations. For double precision, $\epsilon \approx 10^{-16}$.

10. **[Spectral Gap] How does the spectral gap $\lambda_1 / \lambda_2$ affect the power method?**
    ??? success "Solution"
        The error decays as $(\lambda_2/\lambda_1)^k$. A larger gap ensures faster convergence to the dominant eigenvector.

## Chapter Summary

This chapter explores the practical limits of matrix computation:

1. **Error Calculus**: Defined forward and backward stability as the measures of algorithmic reliability.
2. **Structural Complexity**: Benchmarked the $O(n^3)$ cost of matrix operations and explored the path to sub-cubic alternatives.
3. **Iterative Convergence**: Established preconditioning and spectral gaps as the governing factors for high-dimensional solvers.
4. **Numerical Robustness**: Linked the condition number of a problem to the ultimate precision attainable on real hardware.
