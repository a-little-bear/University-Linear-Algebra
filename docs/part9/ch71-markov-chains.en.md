# Chapter 71: Markov Chains

<div class="context-flow" markdown>

**Prerequisites**: Matrix Multiplication (Ch2) · Eigenvalues (Ch6) · Non-negative Matrices (Ch17) · Stochastic Matrices

**Chapter Outline**: Definition of Markov Chains → Transition Probability Matrix → Chapman-Kolmogorov Equations → Classification of States (Irreducible, Aperiodic) → Stationary Distribution → Convergence to Equilibrium → Spectral Gap and Mixing Time → Absorbing Markov Chains → Application in PageRank

**Extension**: Markov chains are the mathematical engine for Google Search (PageRank) and the foundation for MCMC sampling in modern statistical computation.

</div>

Markov chains describe the evolution of stochastic systems where the "future is independent of the past, given the present." This memoryless property is perfectly captured by the linear mapping of probability vectors through **Stochastic Matrices**.

---

## 71.1 Stochastic Matrices and Equilibrium

!!! definition "Definition 71.1 (Stochastic Matrix)"
    A square matrix $P$ is **row-stochastic** if all $p_{ij} \ge 0$ and $\sum_j p_{ij} = 1$ for all $i$. It maps probability vectors to probability vectors.

!!! theorem "Theorem 71.1 (Convergence to Stationary Distribution)"
    If a Markov chain is irreducible and aperiodic (primitive matrix), there exists a unique stationary distribution $\pi$ such that $\pi P = \pi$. Furthermore, for any initial distribution $x_0$, $\lim_{k \to \infty} x_0 P^k = \pi$.

---

## Exercises

1. **[Stochasticity] Prove: If $P$ is a stochastic matrix, then $\lambda=1$ is always an eigenvalue of $P$.**
   ??? success "Solution"
       Let $\mathbf{1}$ be the all-ones column vector. The row-sum condition $\sum_j p_{ij} = 1$ is equivalent to $P \mathbf{1} = 1 \cdot \mathbf{1}$. Thus, $\mathbf{1}$ is a right eigenvector of $P$ for $\lambda=1$. Since $P$ and $P^T$ share the same spectrum, $1 \in \sigma(P)$.

2. **[Stationary Distribution] Find the stationary distribution for $P = \begin{pmatrix} 0.7 & 0.3 \\ 0.4 & 0.6 \end{pmatrix}$.**
   ??? success "Solution"
       Solve $\pi P = \pi$ with $\pi_1 + \pi_2 = 1$.
       $0.7\pi_1 + 0.4\pi_2 = \pi_1 \implies 0.4\pi_2 = 0.3\pi_1 \implies \pi_1 = \frac{4}{3}\pi_2$.
       $\frac{4}{3}\pi_2 + \pi_2 = 1 \implies \frac{7}{3}\pi_2 = 1 \implies \pi_2 = 3/7, \pi_1 = 4/7$.
       $\pi = [4/7, 3/7]$.

3. **[Mixing Time] Explain how the spectral gap $1 - |\lambda_2|$ determines the convergence speed of a Markov chain.**
   ??? success "Solution"
       The error after $k$ steps is roughly $\|x_k - \pi\| \approx C |\lambda_2|^k$. A larger spectral gap implies a smaller $|\lambda_2|$, resulting in faster geometric convergence to the equilibrium state.

4. **[Irreducibility] What is the graph-theoretic condition for a Markov chain to be irreducible?**
   ??? success "Solution"
       The directed graph associated with the transition matrix $P$ must be strongly connected, meaning there is a path of non-zero probability between any two states $i$ and $j$.

5. **[PageRank] Why is the Google PageRank matrix $G = \alpha P + (1-\alpha) \frac{1}{n} J$ guaranteed to have a unique stationary distribution?**
   ??? success "Solution"
       The addition of the term $(1-\alpha)\frac{1}{n}J$ (teleportation) makes the matrix $G$ strictly positive ($G \gg 0$). According to the Perron-Frobenius theorem, a strictly positive matrix is primitive, ensuring a unique strictly positive eigenvector for $\lambda=1$.

6. **[Absorbing States] Define an absorbing state and the structure of an absorbing Markov chain matrix.**
   ??? success "Solution"
       A state $i$ is absorbing if $p_{ii} = 1$. The matrix has the block form $P = \begin{pmatrix} Q & R \\ 0 & I \end{pmatrix}$, where $Q$ corresponds to the transitions between transient states.

7. **[Fundamental Matrix] How is the matrix $(I-Q)^{-1}$ used in absorbing Markov chains?**
   ??? success "Solution"
       $(I-Q)^{-1} = I + Q + Q^2 + \dots$ is the fundamental matrix. The entry $n_{ij}$ represents the expected number of times the chain is in transient state $j$ given it started in state $i$ before absorption.

8. **[Double Stochasticity] Show that if $P$ is doubly stochastic, the uniform distribution $\pi = [\frac{1}{n}, \dots, \frac{1}{n}]$ is a stationary distribution.**
   ??? success "Solution"
       If $P$ is doubly stochastic, then $\mathbf{1}^T P = \mathbf{1}^T$. Dividing by $n$ gives $[\frac{1}{n} \dots \frac{1}{n}] P = [\frac{1}{n} \dots \frac{1}{n}]$. Thus, the uniform distribution is invariant.

9. **[Reversibility] State the detailed balance equation and its implication.**
   ??? success "Solution"
       $\pi_i p_{ij} = \pi_j p_{ji}$. If this holds, the chain is reversible. Reversibility implies that the transition matrix is similar to a symmetric matrix, and thus all its eigenvalues are real.

10. **[Aperiodicity] Provide an example of an irreducible but periodic Markov chain.**
    ??? success "Solution"
        $P = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$. The chain toggles between states 1 and 2. It is irreducible (can reach any state), but periodic with period 2 (can only return to state 1 in an even number of steps).

## Chapter Summary

This chapter explores the linear algebraic properties of stochastic processes:

1. **State Transitions**: Formulated the evolution of probability distributions as matrix-vector multiplications.
2. **Equilibrium Theory**: Utilized the Perron-Frobenius theorem to establish the existence and uniqueness of stationary distributions.
3. **Spectral Dynamics**: Linked the convergence rate (mixing time) to the spectral gap of the transition operator.
4. **Network Ranking**: Demonstrated the power of Markov chains in large-scale information retrieval (PageRank) and transient state analysis.
