# Chapter 71: Markov Chains and Stochastic Processes

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices (Ch17) · Eigenvalues (Ch6) · Basic Probability

**Chapter Outline**: State Space and Transition Matrices → Probability Distribution Vectors → Chapman-Kolmogorov Equations → Steady-state Distribution → Ergodicity → Periodicity and Classification → Absorbing Markov Chains → Fundamental Matrix → PageRank Algorithm

**Extension**: Markov chains are the core algorithmic models for PageRank search engines, financial series modeling, and modern Natural Language Processing (HMM).

</div>

Markov chains utilize non-negative matrices to describe stochastic evolution on discrete state spaces. Their core feature is that the future state of the system depends only on the current operator's action and is independent of its historical trajectory. This property allows the spectral theory of linear algebra to fully characterize the long-term equilibrium behavior of stochastic systems.

---

## 71.1 Transition Probability Matrix and Spectral Structure

!!! definition "Definition 71.1 (Row Stochastic Matrix)"
    The element $p_{ij}$ of a transition matrix $P$ represents the probability of transitioning from state $i$ to state $j$. It satisfies $p_{ij} \ge 0$ and $\sum_j p_{ij} = 1$. The evolution of the probability distribution row vector $\pi_k^T$ follows the linear mapping $\pi_{k+1}^T = \pi_k^T P$.

!!! theorem "Theorem 71.1 (Ergodic Theorem and Unique Equilibrium)"
    If the transition matrix $P$ is irreducible and aperiodic, then according to the Perron-Frobenius theorem, its spectral radius $ho(P) = 1$ is a simple eigenvalue. There exists a unique steady-state distribution $\pi^*$ satisfying $\pi^T P = \pi^T$, and the system converges to this steady state from any initial distribution.

---

## Exercises

1. **[Unity Eigenvalue] Prove: For any row stochastic matrix $P$, 1 is necessarily an eigenvalue.**
   ??? success "Solution"
       Let $\mathbf{1}$ be the all-ones vector. From the definition $\sum_j p_{ij} = 1$, we have $P \mathbf{1} = \mathbf{1} \cdot 1$. This shows that $\mathbf{1}$ is a right eigenvector of $P$ corresponding to the eigenvalue 1.

2. **[Convergence Rate] Prove: The rate at which a Markov chain converges to its steady-state distribution is determined by the second largest eigenvalue $|\lambda_2|$. Define the spectral gap of the system.**
   ??? success "Solution"
       From the spectral decomposition $\pi_k^T = \pi^* + \sum_{i=2}^n c_i \lambda_i^k v_i^T$. As $k 	o \infty$, the error term is dominated by $|\lambda_2|^k$. The spectral gap $1 - |\lambda_2|$ measures how quickly the system approaches equilibrium.

3. **[Calculation] Given the transition matrix $P = \begin{pmatrix} 0.7 & 0.3 \ 0.2 & 0.8 \end{pmatrix}$. Solve for the stationary probability distribution vector.**
   ??? success "Solution"
       Solve the system $\begin{pmatrix} \pi_1 & \pi_2 \end{pmatrix} \begin{pmatrix} -0.3 & 0.3 \ 0.2 & -0.2 \end{pmatrix} = 0$ and $\pi_1 + \pi_2 = 1$.
       From $-0.3\pi_1 + 0.2\pi_2 = 0$, we have $3\pi_1 = 2\pi_2$. The solution is $\pi^* = [0.4, 0.6]^T$.

4. **[Periodicity Determination] Given $P = \begin{pmatrix} 0 & 1 \ 1 & 0 \end{pmatrix}$. Analyze the eigenvalues of this matrix and the asymptotic behavior of the distribution vector.**
   ??? success "Solution"
       The eigenvalues are $\pm 1$. The matrix is periodic (with period 2). The distribution vector will oscillate between two states indefinitely and will not converge to a fixed value, even though a solution to the stationary equation exists.

5. **[PageRank Analysis] Explain how the "damping factor" in Google's PageRank algorithm converts a web link matrix into a primitive matrix.**
   ??? success "Solution"
       The original link matrix may have dangling nodes (zero rows) or cycles. The convex combination $G = d P + (1-d) \frac{1}{n} J$ ensures all elements of $G$ are strictly positive. This eliminates periodicity and guarantees irreducibility, thereby ensuring the unique convergence of the power iteration method.

6. **[Absorbing Chains] Prove: In an absorbing Markov chain, the probability of being absorbed from any state is 1. Define its fundamental matrix $N = (I-Q)^{-1}$.**
   ??? success "Solution"
       Since $Q$ describes transitions between non-absorbing states, its spectral radius $ho(Q) < 1$. The series $I+Q+Q^2+\dots$ converges to $(I-Q)^{-1}$. This represents that the probability of staying in non-absorbing states decays exponentially to zero over time.

7. **[Detailed Balance] Define the reversibility of a Markov chain and explain its relationship to the symmetry of the matrix $\Delta P$ (where $\Delta$ is the steady-state diagonal matrix).**
   ??? success "Solution"
       The detailed balance condition is $\pi_i p_{ij} = \pi_j p_{ji}$. This is equivalent to the diagonal matrix $\Delta_\pi$ satisfying $\Delta_\pi P = P^T \Delta_\pi$, meaning $P$ is a self-adjoint operator under the $\pi$-weighted inner product.

8. **[Mixing Time] Briefly describe the relationship between mixing time and eigenvalues other than the spectral radius.**
   ??? success "Solution"
       Mixing time is the number of steps required for the distribution to reach within $\epsilon$ total variation distance from the steady state. It is inversely proportional to the spectral gap $\gamma = 1 - \lambda_2$. A smaller $\lambda_2$ implies that the random walk can traverse the state space more rapidly.

9. **[Graph Theory Link] Prove: The transition matrix $P$ is irreducible if and only if its associated directed graph is strongly connected.**
   ??? success "Solution"
       Irreducibility means $(I+P)^{n-1} > 0$. This is equivalent to the existence of at least one path of finite length between any two nodes in the graph, which satisfies the definition of strong connectivity.

10. **[Ergodicity] Explain how ergodicity is reflected in matrix language through the projection onto the principal eigenspace.**
    ??? success "Solution"
        Ergodicity implies $P^k 	o \mathbf{1} \pi^T$. This indicates that the matrix power operator, in the long-run evolution, degenerates into a rank-1 projection operator that maps any initial distribution onto the one-dimensional invariant subspace belonging to the eigenvalue 1.

## Chapter Summary

This chapter discusses the linear operator theory describing stochastic evolutionary behavior:

1. **Probability Transition Framework**: Established transition matrices as the core operator model for stochastic state updates.
2. **Equilibrium Theory**: Utilized the Perron-Frobenius theorem to establish the existence, uniqueness, and convergence of steady-state distributions.
3. **Structural Classification**: Distinguished between different Markov systems such as irreducible, aperiodic, and absorbing chains.
4. **Algorithmic Mapping**: Demonstrated how modern algorithms like PageRank transform combinatorial structure analysis into large-scale matrix eigenvector problems.
