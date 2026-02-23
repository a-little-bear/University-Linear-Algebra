# Chapter 17: Nonnegative Matrices and Perron-Frobenius Theory

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Matrix Analysis (Ch14) · Graph Theory Basics (Ch27)

**Chapter Outline**: Definitions of Nonnegative and Positive Matrices → Irreducibility → Perron-Frobenius Theorem → Properties of the Perron Root $\rho(A)$ → Stochastic Matrices → Exponent of a Matrix → Power Method → Applications (PageRank, Population Models, Markov Chains)

**Extension**: Perron-Frobenius theory reveals the "dominant growth mode" of a system, serving as the link between algebraic structure and long-term behavior.

</div>

Nonnegative matrices are ubiquitous in probability, economics, and biology. Perron-Frobenius theory is the crown of nonnegative matrix analysis, ensuring that under certain connectivity conditions, a unique, positive principal eigenvector exists.

---

## 17.1 Perron-Frobenius Theorem

!!! definition "Definition 17.1 (Irreducible Matrix)"
    A nonnegative matrix $A$ is **irreducible** if its associated directed graph is strongly connected (a path exists between any two vertices).

!!! theorem "Theorem 17.3 (Perron-Frobenius Theorem)"
    If $A$ is an irreducible nonnegative matrix, then:
    1. The spectral radius $\lambda = \rho(A)$ is an eigenvalue of $A$ (the Perron root).
    2. There exists a unique (normalized) positive eigenvector $v > 0$ such that $Av = \lambda v$.
    3. $\lambda$ is a simple eigenvalue.

---

## Exercises

1. **[Fundamentals] Is $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ irreducible? Is it primitive?**
   ??? success "Solution"
       - It is irreducible because there are paths $1 \to 2$ and $2 \to 1$.
       - It is not primitive because it is periodic. $A^2 = I, A^3 = A, \dots$ never becomes a strictly positive matrix. Its eigenvalues are $1, -1$, both with modulus 1.

2. **[Perron Root] If each row sum of a nonnegative matrix $A$ is $s$, prove $\rho(A) = s$.**
   ??? success "Solution"
       Let $\mathbf{1} = (1, \dots, 1)^T$. The row-sum condition implies $A \mathbf{1} = s \mathbf{1}$. Thus $s$ is an eigenvalue. Since $A \ge 0$ and all eigenvalues are bounded by the maximum row sum (Gershgorin), $\rho(A) = s$.

3. **[Monotonicity] If $0 \le A \le B$, prove $\rho(A) \le \rho(B)$.**
   ??? success "Solution"
       Since $A, B$ are nonnegative, $A^k \le B^k$ for all $k$. By Gelfand's formula $\rho(A) = \lim \|A^k\|^{1/k} \le \lim \|B^k\|^{1/k} = \rho(B)$. This reflects the monotonicity of the spectral radius for nonnegative matrices.

4. **[Primitivity] Determine if $A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}$ is primitive.**
   ??? success "Solution"
       Calculate $A^2 = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix} > 0$. Since a power of the matrix is strictly positive, $A$ is primitive.

5. **[Collatz-Wielandt] Use row sums to estimate the spectral radius of $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.**
   ??? success "Solution"
       Row sums are 3 and 7. Thus $3 \le \rho(A) \le 7$. The exact value is $\frac{5+\sqrt{13}}{2} \approx 5.30$.

6. **[Application] In PageRank, why is a damping factor (adding a multiple of an all-ones matrix) used?**
   ??? success "Solution"
       To make the matrix **primitive**. This ensures the power method converges to a unique positive stationary distribution (the rank vector), eliminating convergence issues caused by isolated nodes or dead loops in the original graph.

7. **[Spectral Distribution] If $A > 0$, what condition do the other eigenvalues' moduli satisfy relative to the Perron root?**
   ??? success "Solution"
       For a primitive matrix (especially $A>0$), all other eigenvalues satisfy $|\lambda_i| < \rho(A)$ ($i \ge 2$).

8. **[Stochastic Matrix] Prove that the spectral radius of a stochastic matrix (row sums = 1) is 1.**
   ??? success "Solution"
       $A \mathbf{1} = 1 \mathbf{1}$ shows 1 is an eigenvalue. By Gershgorin, all eigenvalues have modulus $\le 1$. Thus $\rho(A) = 1$.

9. **[Irreducibility Test] Determine the irreducibility of $\begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       Reducible. Node 1 cannot reach node 2 (only $2 \to 1$ exists). The matrix is lower triangular, reflecting its reducibility.

10. **[Limit Behavior] If $A$ is primitive and $\rho(A)=1$, prove $\lim_{k \to \infty} A^k = v w^T$, where $v, w$ are the Perron eigenvectors.**
    ??? success "Solution"
        Since all other eigenvalues have modulus $< 1$, in the spectral decomposition, only the term corresponding to $\lambda=1$ (the rank-1 projection) survives in the limit.

## Chapter Summary

Nonnegative matrix theory interweaves analysis and combinatorics:

1. **Growth Dominance**: Spectral radius is no longer just an abstract value but the actual expansion rate of the system.
2. **Positivity Guarantee**: Irreducibility is the structural prerequisite for ensuring every component in the system can evolve positively.
3. **Inevitability of Convergence**: Primitivity establishes the unique endpoint for the evolution of discrete dynamical systems toward a steady state.
