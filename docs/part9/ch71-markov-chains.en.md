# Chapter 71: Markov Chains and Stochastic Matrices

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices & Perron-Frobenius (Ch17) · Matrix Analysis (Ch14) · Probability Theory

**Chapter Outline**: Definition of Discrete-Time Markov Chains (DTMC) → The Transition Probability Matrix $P$ → State Classification: Irreducible, Aperiodic, Recurrent, and Transient → Stationary Distributions ($\pi P = \pi$) → Convergence Theorem: Perron-Frobenius in Stochastic Processes → Absorbing Markov Chains & the Fundamental Matrix $(I-Q)^{-1}$ → Mean First Passage Time → Applications: Google PageRank, Queueing Theory, and Financial Modeling

**Extension**: Markov chains are the algebraic template for describing "memoryless" evolution systems; they prove that the long-term trend of complex systems depends solely on the spectral structure of the transition matrix, serving as the mathematical heart of modern Monte Carlo methods (MCMC).

</div>

The future depends only on the present, not on the past. This simple physical intuition evolved into the most powerful model in probability theory: the **Markov Chain**. By arranging the transition probabilities of a system into a matrix, we transform the evolution of a stochastic process into problems of matrix powers and eigenvalue analysis. This chapter demonstrates how linear algebra provides deterministic asymptotic predictions for an uncertain world.

---

## 71.1 Transition Matrices and Distribution Evolution

!!! definition "Definition 71.1 (Transition Matrix)"
    For a system with $n$ states, the **Transition Matrix** $P$ has entries $p_{ij}$ representing the probability of moving from state $i$ to state $j$ in one time step.
    - **Stochasticity**: $P$ is a **row-stochastic** matrix (each row sums to 1).
    - **Evolution**: If the current probability distribution is the row vector $\pi_k$, the distribution at the next step is $\pi_{k+1} = \pi_k P$.

---

## 71.2 Stationary Distributions and Spectral Analysis

!!! definition "Definition 71.2 (Stationary Distribution)"
    A probability vector $\pi$ satisfying $\pi P = \pi$ and $\sum \pi_i = 1$ is a **Stationary Distribution**.
    **Algebraic Essence**: $\pi^T$ is the eigenvector of $P^T$ associated with the eigenvalue 1.

!!! theorem "Theorem 71.1 (Convergence Theorem)"
    If $P$ is irreducible and aperiodic (a primitive matrix), then for any initial distribution $\pi_0$:
    $$\lim_{k \to \infty} \pi_0 P^k = \pi$$
    **Rate**: The convergence speed is governed by the second largest eigenvalue $|\lambda_2|$ (the spectral gap).

---

## 71.3 Absorbing Markov Chains

!!! technique "Absorbing States"
    A chain is **Absorbing** if some states cannot be left once entered. Its transition matrix can be partitioned as $P = \begin{pmatrix} Q & R \\ 0 & I \end{pmatrix}$.
    - **Fundamental Matrix**: $N = (I - Q)^{-1} = I + Q + Q^2 + \cdots$.
    - **Meaning**: $N_{ij}$ represents the expected number of times the system stays in transient state $j$ before being absorbed, starting from state $i$.

---

## Exercises

1.  **[Basics] Determine if $P = \begin{pmatrix} 0.5 & 0.5 \\ 0.2 & 0.8 \end{pmatrix}$ is a stochastic matrix.**
    ??? success "Solution"
        Yes. Each entry is non-negative and each row sums to 1.

2.  **[Stationary] Find the stationary distribution for the matrix $P$ above.**
    ??? success "Solution"
        Solve $(\pi_1, \pi_2) \begin{pmatrix} 0.5 & 0.5 \\ 0.2 & 0.8 \end{pmatrix} = (\pi_1, \pi_2)$.
        This gives $0.5\pi_1 + 0.2\pi_2 = \pi_1 \implies 0.2\pi_2 = 0.5\pi_1 \implies \pi_2 = 2.5\pi_1$.
        With $\pi_1 + \pi_2 = 1$, we get $\pi = (2/7, 5/7)$.

3.  **[Spectrum] Prove that the spectral radius of any row-stochastic matrix is exactly 1.**
    ??? success "Solution"
        Since $P \mathbf{1} = \mathbf{1}$, 1 is an eigenvalue. Since $\|P\|_\infty = 1$, the spectral radius $\rho(P) \le \|P\|_\infty = 1$. Thus, $\rho(P) = 1$.

4.  **[Periodicity] Does the system with $P = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ converge?**
    ??? success "Solution"
        No. It is a periodic matrix (cycle of length 2). While it has a stationary distribution $(0.5, 0.5)$, the state will oscillate between $(1, 0)$ and $(0, 1)$ indefinitely.

5.  **[Fundamental] In an absorbing chain with $Q = (0.5)$, find the fundamental matrix $N$.**
    ??? success "Solution"
        $N = (1 - 0.5)^{-1} = 2$. This means the system stays in the transient state twice on average before absorption.

6.  **[Application] Why is PageRank a Markov chain?**
    ??? success "Solution"
        Web pages are states and links are transitions. A random surfer's behavior defines the transition matrix, and the rankings are the components of its stationary distribution.

7.  **[Uniqueness] Prove: If $P$ is a strictly positive matrix, the stationary distribution is unique.**
    ??? success "Solution"
        By the Perron-Frobenius theorem, a positive matrix has a unique dominant eigenvalue 1 with a strictly positive eigenvector.

8.  **[Limit] Find $P^k$ as $k \to \infty$ for $P = \begin{pmatrix} 1 & 0 \\ 0.5 & 0.5 \end{pmatrix}$.**
    ??? success "Solution"
        State 1 is absorbing. $Q=(0.5), R=(0.5)$. $P^\infty = \begin{pmatrix} 1 & 0 \\ 1 & 0 \end{pmatrix}$. All units are eventually absorbed into state 1.

9.  **[Convergence] Does a larger spectral gap $1 - |\lambda_2|$ imply faster or slower convergence?**
    ??? success "Solution"
        Faster convergence. The error terms decay as $|\lambda_2|^k$.

10. **[Doubly Stochastic] What is the stationary distribution of a doubly stochastic matrix?**

   ??? success "Solution"
        The uniform distribution $\pi = (1/n, \ldots, 1/n)$ because $\mathbf{1} P = \mathbf{1}$.

## Chapter Summary

Markov chains represent the ultimate resonance of probability and algebra:

1.  **Limits of Determinism**: They prove that the long-term behavior of a random process is not random, but a deterministic result pre-encoded in the transition matrix's Perron vector.
2.  **Structural Stability**: Through state classification and irreducibility analysis, the theory links the rate of "memory loss" in a system to its network topology (spectral gap).
3.  **Engineering of Computation**: From inverting absorbing states to PageRank's power iteration, Markov chains transform profound statistical philosophy into efficient linear algebra algorithms, establishing standard protocols for handling uncertain systems.
