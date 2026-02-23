# Chapter 71: Markov Chains

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices & Perron-Frobenius (Ch17) · Probability Theory · Matrix Analysis (Ch14)

**Chapter Outline**: Stochastic Flow of States → Markov Property and the Transition Matrix $P$ → Row-stochastic Matrices → Stationary Distribution $\boldsymbol{\pi}$ as an Eigenvector → Classification of States: Absorbing, Transient, and Ergodic → Key Theorem: Fundamental Convergence Theorem (Dominance of Eigenvalue 1) → Mixing Time and the Spectral Gap → Fundamental Matrix of Absorbing Markov Chains → Applications: Random Walk interpretation of Google PageRank, Win-rate Prediction in Games, and MCMC Sampling Algorithms (Metropolis-Hastings)

**Extension**: Markov chains are the "probabilistic dynamization" of linear algebra; they encapsulate the future of a system entirely within the current probability vector. They prove that long-term stability is essentially the convergence of operator eigenspaces—the algebraic foundation for modern stochastic simulation and reinforcement learning.

</div>

In many systems, the future state depends only on the present, not on the past. This "memoryless" property is known as the **Markov Property**. Linear algebra provides the perfect tool for describing this stochastic flow: the **Transition Matrix**. By studying the spectral properties of this matrix, we can foresee whether a system, after countless random fluctuations, will converge to a definite equilibrium. This chapter explores the algebraic system connecting probabilistic logic with operator iteration.

---

## 71.1 Transition Matrices and Stationary Distributions

!!! definition "Definition 71.1 (Transition Matrix $P$)"
    For a system with $n$ states, the entry $p_{ij}$ of $P$ represents the probability of transitioning from state $i$ to state $j$.
    **Property**: $P$ is row-stochastic, meaning the sum of each row is exactly 1.

!!! definition "Definition 71.2 (Stationary Distribution $\boldsymbol{\pi}$)"
    A probability vector $\boldsymbol{\pi}$ is a **stationary distribution** if $\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T$.
    **Algebraic Essence**: It is the left eigenvector of $P$ corresponding to the eigenvalue 1.

---

## 71.2 Convergence and Spectral Gaps

!!! theorem "Theorem 71.1 (Fundamental Convergence Theorem)"
    If a Markov chain is irreducible and aperiodic (i.e., $P$ is primitive), then for any initial distribution $\mathbf{x}_0$:
    $$\lim_{k \to \infty} \mathbf{x}_0^T P^k = \boldsymbol{\pi}^T$$
    **Speed**: The rate of convergence is determined by the second-largest eigenvalue magnitude $|\lambda_2|$. The value $1 - |\lambda_2|$ is the **Spectral Gap**.

---

## 71.3 Absorbing Markov Chains

!!! technique "Technique: The Fundamental Matrix $N$"
    For a chain with absorbing states (states that once entered, cannot be left), the transition matrix can be partitioned as $\begin{pmatrix} Q & R \\ 0 & I \end{pmatrix}$. The **Fundamental Matrix** is $N = (I - Q)^{-1}$.
    The entry $N_{ij}$ represents the expected number of times the system stays in transient state $j$ before absorption, given it started in state $i$.

---

## Exercises

**1. [Basics] Determine if $\begin{pmatrix} 0.5 & 0.5 \\ 0.2 & 0.8 \end{pmatrix}$ is a valid transition matrix.**

??? success "Solution"
    **Verification:**
    1. All entries are non-negative.
    2. Row 1 sum: $0.5 + 0.5 = 1$.
    3. Row 2 sum: $0.2 + 0.8 = 1$.
    **Conclusion**: Yes, it is a valid row-stochastic matrix.

**2. [Calculation] Find the stationary distribution $\boldsymbol{\pi}$ for the matrix in the previous problem.**

??? success "Solution"
    **Solve $\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T$:**
    1. Let $\boldsymbol{\pi} = (p, 1-p)^T$.
    2. Component 1 equation: $0.5p + 0.2(1-p) = p$.
    3. $0.2 - 0.2p = 0.5p \implies 0.2 = 0.7p \implies p = 2/7$.
    **Conclusion**: $\boldsymbol{\pi} = (2/7, 5/7)^T$.

**3. [Property] Prove that every row-stochastic matrix has an eigenvalue of 1.**

??? success "Solution"
    **Proof:**
    1. Let $\mathbf{1} = (1, 1, \ldots, 1)^T$.
    2. $P\mathbf{1}$ results in a vector where each component is the sum of a row of $P$.
    3. Since row sums are 1, $P\mathbf{1} = \mathbf{1}$.
    **Conclusion**: 1 is an eigenvalue corresponding to the all-ones eigenvector.

**4. [Aperiodicity] Determine the convergence of $P = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.**

??? success "Solution"
    **Determination:**
    1. This is a periodic matrix with period 2.
    2. Eigenvalues are $\{1, -1\}$.
    3. Since there is an eigenvalue with modulus 1 other than 1 itself (namely -1), the sequence $P^k$ oscillates and does not converge to a single distribution.
    **Conclusion**: The chain does not converge to a steady state.

**5. [Absorbing] Given states $\{1, 2, 3\}$ where 3 is absorbing and $Q = \begin{pmatrix} 0.5 & 0.1 \\ 0.1 & 0.5 \end{pmatrix}$, find $N$.**

??? success "Solution"
    **Steps:**
    1. $I - Q = \begin{pmatrix} 0.5 & -0.1 \\ -0.1 & 0.5 \end{pmatrix}$.
    2. $\det = 0.25 - 0.01 = 0.24$.
    3. $N = (I-Q)^{-1} = \frac{1}{0.24} \begin{pmatrix} 0.5 & 0.1 \\ 0.1 & 0.5 \end{pmatrix} \approx \begin{pmatrix} 2.08 & 0.42 \\ 0.42 & 2.08 \end{pmatrix}$.
    **Meaning**: Starting in state 1, the system expects to stay in state 1 for 2.08 steps and state 2 for 0.42 steps before being absorbed by state 3.

**6. [Convergence Rate] If $|\lambda_2| = 0.9$, roughly how many iterations are needed to reduce error to $1/e$?**

??? success "Solution"
    **Estimation:**
    Error decays as $|\lambda_2|^k$.
    $0.9^k \approx e^{-1} \implies k \ln(0.9) \approx -1 \implies k \approx 1/0.1 = 10$.
    **Conclusion**: Approximately 10 transitions. A larger spectral gap $1-|\lambda_2|$ leads to faster convergence.

**7. [Balance] What is the "Detailed Balance" condition?**

??? success "Solution"
    **Definition:**
    $\pi_i p_{ij} = \pi_j p_{ji}$.
    If this holds, the Markov chain is **reversible**. This condition guarantees the existence of a stationary distribution and is the core principle for constructing transition probabilities in MCMC algorithms like the Metropolis algorithm.

**8. [Calculation] For a $2 \times 2$ matrix $P$ with rows $(1-\alpha, \alpha)$ and $(\beta, 1-\beta)$, find the eigenvalues.**

??? success "Solution"
    **Conclusion: The eigenvalues are 1 and $1-\alpha-\beta$.**
    This illustrates how the spectral gap is directly determined by the "staying probabilities" of the states.

**9. [Application] Briefly state the algebraic meaning of Markov chains in the "infinite monkey" experiment.**

??? success "Solution"
    Text generation can be modeled as a Markov process. The transition matrix encodes the statistical features of a language (e.g., 'u' follows 'q' with high probability). By studying the matrix powers, one can calculate the expected "waiting time" to generate a specific word like "HAMLET," which is essentially solving for the expected steps in an absorbing chain.

**10. [Application] What is "Burn-in" in MCMC?**

??? success "Solution"
    Because the initial distribution $\mathbf{x}_0$ may be far from the stationary $\boldsymbol{\pi}$, early samples do not reflect the target distribution. Linear algebra shows that after $k \gg 1/(1-|\lambda_2|)$ iterations, the state vector enters the stable subspace dominated by eigenvalue 1. This discarded initial phase is the "Burn-in" period.

## Chapter Summary

Markov chains are the ultimate framework for describing stochastic processes in linear algebra:

1.  **Deterministic Probability**: They prove that macroscopic stationary distributions are the inevitable result of spectral structures, elevating randomness to algebraic necessity.
2.  **Physical Meaning of Spectra**: The spectral gap $1-|\lambda_2|$ establishes the speed at which a system "forgets the past," serving as the benchmark for evaluating the efficiency of all iterative sampling algorithms.
3.  **Evolution of Structure**: From residence times in absorbing chains to mixing times in ergodic ones, matrix analysis provides a unified engine for understanding queuing, search, and diffusion in the real world.
