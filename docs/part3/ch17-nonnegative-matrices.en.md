# Chapter 17: Nonnegative Matrices and Perron-Frobenius Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Eigenvalues (Ch6) · Graph Theory (Ch27) · Convergence (Ch14)

**Chapter Outline**: Definition of Nonnegative and Positive Matrices → Perron's Theorem for Positive Matrices → Reducibility and Irreducibility → Frobenius's Theorem for Irreducible Matrices → Primitive Matrices → Stochastic Matrices → Convergence of Markov Chains → Google PageRank Algorithm

**Extension**: Nonnegative matrices are the mathematical language of probability, population dynamics, and economic networks; the Perron-Frobenius theorem is the "Spectral Theorem" for non-symmetric positive operators.

</div>

Nonnegative matrices—those with all entries $a_{ij} \ge 0$—govern systems where quantities cannot be negative, such as probability distributions, population counts, or economic values. The **Perron-Frobenius theorem** is one of the most elegant results in matrix theory: it guarantees that such a matrix always has a unique, positive "dominant" eigenvalue and a corresponding positive eigenvector. This theorem provides the rigorous justification for the convergence of Markov chains and the ranking logic of search engines like Google.

---

## 17.1 Perron-Frobenius Theory

!!! definition "Definition 17.1 (Positive and Irreducible)"
    A matrix $A$ is **positive** ($A > 0$) if all $a_{ij} > 0$.
    A nonnegative matrix $A$ is **irreducible** if its associated directed graph is strongly connected.

!!! theorem "Theorem 17.1 (Perron-Frobenius Theorem)"
    If $A$ is an irreducible nonnegative matrix, then:
    1. The spectral radius $\rho(A)$ is an eigenvalue (the **Perron root**).
    2. $\rho(A) > 0$ and it is a simple eigenvalue.
    3. There exists a strictly positive eigenvector $v > 0$ such that $Av = \rho(A)v$.
    4. $\rho(A)$ is the unique eigenvalue with a positive eigenvector.

---

## Exercises

1. **[Fundamentals] Is $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ irreducible? Is it positive?**
   ??? success "Solution"
       Irreducible: Yes, the graph $1 \leftrightarrow 2$ is strongly connected. Positive: No, it has zeros on the diagonal.

2. **[Perron Root] Find the Perron root of $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       The eigenvalues are $\{2, 0\}$. The Perron root is 2. The corresponding eigenvector is $(1, 1)^T > 0$.

3. **[Stochastic] Define a row-stochastic matrix and its Perron root.**
   ??? success "Solution"
       A matrix $P \ge 0$ where each row sum is 1. Its Perron root is always $\lambda=1$, with the all-ones vector $\mathbf{1}$ as the right eigenvector.

4. **[Primitivity] What is a primitive matrix?**
   ??? success "Solution"
       An irreducible matrix $A$ such that $A^k > 0$ for some $k$. Primitive matrices have a unique eigenvalue of maximum modulus, ensuring convergence to a steady state.

5. **[PageRank] How does PageRank use Perron-Frobenius?**
   ??? success "Solution"
       PageRank models the web as a graph and its link structure as a stochastic matrix. The PageRank scores are the entries of the principal eigenvector (Perron vector), which exists and is unique by the theorem.

6. **[Reducibility] What happens if $A$ is reducible?**
   ??? success "Solution"
       The matrix can be permuted to a block upper-triangular form. The Perron root still exists, but the corresponding eigenvector might only be non-negative (not strictly positive), and the root might not be simple.

7. **[Collatz-Wielandt] State the Collatz-Wielandt formula for the Perron root.**
   ??? success "Solution"
       $\rho(A) = \max_{x > 0} \min_i \frac{(Ax)_i}{x_i}$. This provides a way to bound the Perron root using any positive vector.

8. **[Cycle] If $A$ is irreducible and has period $h > 1$, what can you say about its eigenvalues?**
   ??? success "Solution"
       The spectrum is invariant under rotation by $2\pi/h$ in the complex plane. There are $h$ eigenvalues on the spectral circle $|z| = \rho(A)$.

9. **[Population] In a Leslie matrix model, what does the Perron root represent?**
   ??? success "Solution"
       The long-term growth rate of the population. If $\rho(L) > 1$, the population grows; if $\rho(L) < 1$, it faces extinction.

10. **[M-matrices] Relate nonnegative matrices to M-matrices (Ch38A).**
    ??? success "Solution"
        If $A \ge 0$, then $M = sI - A$ is an M-matrix for any $s > \rho(A)$. M-matrices are the "inverse-positive" counterparts to nonnegative matrices.

## Chapter Summary

This chapter explores the spectral properties of operators that preserve the positive orthant:

1. **Spectral Positivity**: Established the Perron root as the unique dominant real eigenvalue for non-negative systems.
2. **Structural Connectivity**: Linked the irreducibility of a matrix to the strong connectivity of its underlying graph.
3. **Equilibrium Logic**: Formulated the theory of steady states for stochastic processes and search algorithms.
4. **Dominant Invariants**: Developed the Perron vector as the definitive descriptor of long-term ratios in linear dynamical systems.
