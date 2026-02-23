# Chapter 17: Non-negative Matrices and Perron-Frobenius Theory

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Analysis (Ch14) · Graph Theory Basics (Ch27)

**Chapter Outline**: Motivation for Non-negative Matrices (Probability, Population, Money) → Matrix Irreducibility → Primitive Matrices → Perron’s Theorem (Positive Matrices) → Frobenius’s Theorem (Irreducible Non-negative Matrices) → Algebraic and Geometric Dominance of the Spectral Radius $\rho(A)$ → The Perron Vector (Steady State) → Collatz-Wielandt Formula → Applications: Google’s PageRank Algorithm, Leontief Economic Equilibrium, and Leslie Population Models

**Extension**: Perron-Frobenius theory is the ultimate tool for describing "positive feedback" and "stable flow"; it reveals that as long as a system is strongly connected and non-negative, its long-term evolution will converge to a physically meaningful positive steady state—the heart of discrete dynamical system analysis.

</div>

In physical, biological, and economic models, most variables (such as population, probability, and price) are logically non-negative. **Perron-Frobenius Theory** is the crown jewel for studying the evolution of such systems. It breaks the convention that eigenvalues might be complex or negative, proving that under non-negative constraints, a matrix must possess a "dominant" positive eigenvalue and an associated positive eigenvector. This chapter establishes the profound framework linking topological connectivity to algebraic stability.

---

## 17.1 Irreducibility and Graph Criteria

!!! definition "Definition 17.1 (Non-negative and Positive Matrices)"
    1.  **Non-negative Matrix $A \ge 0$**: All entries $a_{ij} \ge 0$.
    2.  **Positive Matrix $A > 0$**: All entries $a_{ij} > 0$.

!!! definition "Definition 17.2 (Irreducible Matrix)"
    A non-negative matrix $A$ is **irreducible** if there is no permutation matrix $P$ such that $P^T A P$ is in block upper triangular form.
    **Graph Equivalent**: The associated directed graph of the matrix is **strongly connected** (i.e., there is a path between any two nodes).

---

## 17.2 The Perron-Frobenius Theorem

!!! theorem "Theorem 17.1 (Theorem for Irreducible Non-negative Matrices)"
    Let $A \ge 0$ be irreducible, and let $r = \rho(A)$ be its spectral radius. Then:
    1.  **Positivity**: $r$ is an eigenvalue of $A$ and $r > 0$.
    2.  **Simplicity**: $r$ is a simple root of the characteristic equation (algebraic multiplicity 1).
    3.  **Positive Eigenvector**: The eigenvector $\mathbf{v}$ associated with $r$ can be chosen strictly positive ($v_i > 0$).
    4.  **Monotonicity**: $r$ increases strictly as any entry of the matrix increases.

!!! theorem "Theorem 17.2 (Dominance of Primitive Matrices)"
    If $A$ is **primitive** (i.e., there exists $k$ such that $A^k > 0$), then for any other eigenvalue $\lambda$, $|\lambda| < r$.
    **Physical Meaning**: This means that regardless of the initial state, the long-term evolution of the system is governed by $r$ and its associated positive vector (steady-state distribution).

---

## 17.3 Stochastic Matrices and Computation

!!! technique "Technique: Collatz-Wielandt Formula"
    The spectral radius $r$ can be found via the following extremal problem:
    $$r = \max_{\mathbf{x} > 0} \min_{i} \frac{(A\mathbf{x})_i}{x_i}$$
    This provides a numerical path to estimate the dominant eigenvalue without solving the characteristic equation explicitly.

---

## Exercises

**1. [Basics] Determine if $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ is positive, non-negative, and irreducible.**

??? success "Solution"
    **Determination:**
    1. **Non-negativity**: Entries are 0 or 1, so it is **non-negative**.
    2. **Positivity**: Contains a 0, so it is **not positive**.
    3. **Irreducibility**: The associated graph is $1 \leftrightarrow 2$, which is strongly connected. Alternatively, $A+A^2 = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} > 0$. Thus it is **irreducible**.

**2. [Eigenvalues] Find the eigenvalues of $A$ from the previous problem and verify the Perron-Frobenius conclusions.**

??? success "Solution"
    **Calculation:**
    Characteristic equation: $\lambda^2 - 1 = 0 \implies \lambda = \pm 1$.
    1. Spectral radius $r = \rho(A) = 1$.
    2. $r=1$ is indeed an eigenvalue and is real.
    3. The eigenvector for $r=1$ is $(1, 1)^T$, which is strictly positive.
    Note: Since the matrix is not primitive (it is cyclic), there exists another eigenvalue $-1$ with modulus equal to $r$.

**3. [Irreducibility] Determine if $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ is irreducible.**

??? success "Solution"
    **Graph Observation:**
    1. Node 1 has a self-loop and an edge to node 2.
    2. Node 2 has no edge back to node 1 ($a_{21}=0$).
    3. The graph is not strongly connected.
    **Conclusion**: It is a reducible matrix. It is already in block upper triangular form.

**4. [Stochastic] Prove that the spectral radius of a row-stochastic matrix (row sums = 1) is exactly 1.**

??? success "Solution"
    **Proof:**
    1. Let $\mathbf{1} = (1, 1, \ldots, 1)^T$.
    2. The $i$-th component of $A\mathbf{1}$ is the sum of the $i$-th row.
    3. By definition, $A\mathbf{1} = 1 \cdot \mathbf{1}$.
    4. Thus 1 is an eigenvalue, so $\rho(A) \ge 1$.
    5. Since $\|A\|_\infty = \max (\text{row sums}) = 1$, and $\rho(A) \le \|A\|$, we have $\rho(A) \le 1$.
    **Conclusion**: $\rho(A) = 1$.

**5. [Primitivity] Give an example of an irreducible but non-primitive matrix.**

??? success "Solution"
    **Example: The permutation matrix** $P = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.
    **Analysis**: It is irreducible (strongly connected), but its powers $P, I, P, I \ldots$ always contain zeros and can never become strictly positive. Such matrices are called **cyclic**.

**6. [Perron Vector] In PageRank, why is a strictly positive perturbation (damping factor) added to the link matrix?**

??? success "Solution"
    **Algebraic Explanation:**
    1. The raw link graph may not be strongly connected (reducible) or may be cyclic (non-primitive).
    2. This could lead to a non-unique dominant eigenvalue or lack of convergence in power iterations.
    3. Adding the perturbation (the Google matrix) ensures the matrix is a **primitive positive matrix**.
    4. By Perron's theorem, the dominant eigenvalue 1 is then unique and strictly dominant, ensuring convergence to a unique ranking from any starting vector.

**7. [Properties] Prove: If $A \ge 0$ is irreducible, then $(I + A)^{n-1} > 0$.**

??? success "Solution"
    **Logic:**
    1. $(I+A)^{n-1} = \sum \binom{n-1}{k} A^k$.
    2. In graph terms, a non-zero entry in $A^k$ at $(i,j)$ means there exists a path of length $k$.
    3. Irreducibility means there is a path of length $\le n-1$ between any two nodes.
    4. For every pair $(i,j)$, at least one $A^k$ in the sum will have a positive entry.
    **Conclusion**: The sum is strictly positive, showing that self-feedback turns an irreducible structure into a primitive one.

**8. [Comparison] If $0 \le A \le B$ (entry-wise), prove $\rho(A) \le \rho(B)$.**

??? success "Solution"
    **Derivation:**
    1. Use the Collatz-Wielandt formula: $\rho(A) = \inf_{\mathbf{x} > 0} \max_i \frac{(A\mathbf{x})_i}{x_i}$.
    2. Since $A \le B$ and $\mathbf{x} > 0$, we have $A\mathbf{x} \le B\mathbf{x}$.
    3. For any $\mathbf{x}$, the growth ratio under $A$ is less than or equal to that under $B$.
    **Conclusion**: The spectral radius of a non-negative matrix is a monotonic non-decreasing function of its entries.

**9. [Application] What is the steady state of the Leontief production equation?**

??? success "Solution"
    **Explanation:**
    In a closed model, $Ax = x$, meaning output exactly matches consumption. This corresponds to the Perron eigenvector of the consumption matrix $A$. Perron-Frobenius theory ensures that this equilibrium output vector $\mathbf{x}$ is strictly positive, which is the only physically meaningful result in economics.

**10. [Limits] Let $A > 0$ be a positive matrix with $\rho(A) = 1$. Prove $\lim_{k \to \infty} A^k = \mathbf{v}\mathbf{w}^T$.**

??? success "Solution"
    **Proof:**
    1. Since $A > 0$, it is primitive, and 1 is the unique eigenvalue of maximum modulus.
    2. Write the Jordan decomposition: $A = \mathbf{v}\mathbf{w}^T + \sum_{j \ge 2} \lambda_j J_j$.
    3. Take powers: $A^k = (\mathbf{v}\mathbf{w}^T)^k + \sum \lambda_j^k J_j^k = \mathbf{v}\mathbf{w}^T + \sum \lambda_j^k J_j^k$.
    4. Since all other $|\lambda_j| < 1$, as $k \to \infty$, these terms vanish.
    **Conclusion**: The power of the matrix collapses to a rank-1 matrix formed by the outer product of the principal eigenvectors.

## Chapter Summary

Non-negative matrix theory establishes a harmonic resonance between physical laws and algebraic structures:

1.  **Preservation of Positivity**: Perron-Frobenius theory establishes the positivity of the principal eigenvector, providing unique and valid steady-state solutions for problems in probability and resource allocation.
2.  **Structural Connectivity**: Irreducibility links matrix algebra to graph topology, proving that "global circulation" is a prerequisite for a system to converge to a unique stable state.
3.  **Computational Convergence**: Primitivity ensures that long-term behavior is non-oscillatory, revealing the final destination of information flow in complex networks.
