# Chapter 65B: Combinatorial Structures

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch3) · Graph Theory (Ch27) · Eigenvalues (Ch6) · Sign Patterns (Ch65A) · Permanents (Ch40a)

**Chapter Outline**: Minimum Rank of a Graph → Colin de Verdière Parameters → Zero Forcing Number → Exponent of a Primitive Matrix → Scrambling Matrices and Ergodicity → Hadamard Matrices → Permanent Bounds → Tournament Matrices → Equitable Partitions and Quotient Matrices → Sinkhorn-Knopp Scaling

**Extension**: Minimum rank problems link quantum information (independent sets of graph states) to communication complexity; Hadamard matrices are central to error-correcting codes and compressed sensing.

</div>

This chapter examines the combinatorial properties of matrices, specifically how graph theory and discrete mathematics constrain algebraic properties. We explore the minimum rank of a graph, the zero forcing number, and discrete structures with extreme properties like Hadamard matrices.

---

## 65B.1 Spectral Parameters of Graphs

!!! definition "Definition 65B.1 (Minimum Rank of a Graph)"
    $\operatorname{mr}(G) = \min\{\operatorname{rank}(A) : A^T=A, \operatorname{sgn}(A) 	ext{ matches } G\}$. This parameter measures the limit of information compression while maintaining a specific connectivity structure.

!!! theorem "Theorem 65B.2 ($Z(G) \ge M(G)$)"
    The zero forcing number $Z(G)$ is a pure combinatorial upper bound on the maximum nullity $M(G)$ of a graph. This provides a method to estimate the rank of a matrix via a simple graph-coloring game.

---

## 65B.2 Discrete Matrix Structures

!!! theorem "Theorem 65B.7 (Hadamard Determinant Bound)"
    For a real $n 	imes n$ matrix with entries $|a_{ij}| \le 1$, $|\det A| \le n^{n/2}$. Equality is attained if and only if $A$ is a Hadamard matrix.

!!! theorem "Theorem 65B.12 (Sinkhorn-Knopp Theorem)"
    A nonnegative matrix has total support if and only if it can be scaled by diagonal matrices to be doubly stochastic.

---

## Exercises

1. **[Fundamentals] Compute the minimum rank $\operatorname{mr}(P_3)$ of the path graph on 3 vertices.**
   ??? success "Solution"
       $P_3$ corresponds to $3 	imes 3$ symmetric tridiagonal matrices. By eigenvalue interlacing, if the off-diagonal entries are non-zero, the rank is at least $n-1 = 2$. Thus $\operatorname{mr}(P_3) = 2$.

2. **[Zero Forcing] Show that the zero forcing number of the cycle graph $C_4$ is 2.**
   ??? success "Solution"
       Let the vertices be $\{1, 2, 3, 4\}$. Color $\{1, 2\}$ black. Since 1 has only 2 and 4 as neighbors and 2 is already black, 1 forces 4 to be black. Similarly, 4 then forces 3 to be black. Total blackening requires 2 initial points.

3. **[Colin de Verdière] Why does $\mu(G) \le 3$ characterize planar graphs?**
   ??? success "Solution"
       This was Colin de Verdière's deep discovery. He constructed spectral parameters such that topological properties like planarity were mapped to the maximum multiplicity of the second-smallest eigenvalue. For $K_5$ (non-planar), this multiplicity is 4.

4. **[Hadamard] Prove that a $3 	imes 3$ Hadamard matrix does not exist.**
   ??? success "Solution"
       For $n > 2$, a Hadamard matrix must have order $n \equiv 0 \pmod 4$. Since 3 is not a multiple of 4, no such matrix exists.

5. **[Permanent] Compute the permanent $\operatorname{perm}(J_2)$ of the $2 	imes 2$ all-ones matrix.**
   ??? success "Solution"
       $\operatorname{perm}(J_2) = 1 \cdot 1 + 1 \cdot 1 = 2$.

6. **[Tournaments] What is a Tournament Matrix? What is unique about its diagonal?**
   ??? success "Solution"
       A tournament matrix is an orientation of a complete graph. $t_{ij}=1$ means $i$ beats $j$. Since no player can beat themselves, the diagonal is always zero.

7. **[Score Sequences] Determine if the sequence $(1, 1, 1)$ is a score sequence for a tournament of order 3.**
   ??? success "Solution"
       Check Landau's conditions: $1 \ge \binom{1}{2}=0$; $1+1=2 \ge \binom{2}{2}=1$; $1+1+1=3 = \binom{3}{2}=3$. All conditions are met. This is the score sequence of a 3-cycle tournament.

8. **[Wielandt] For an $n=3$ primitive matrix, what is the maximum possible exponent $\gamma(A)$?**
   ??? success "Solution"
       According to the Wielandt bound: $(n-1)^2 + 1 = (3-1)^2 + 1 = 5$. This ensures $A^5$ is strictly positive.

9. **[Sinkhorn] What is the purpose of the Sinkhorn-Knopp algorithm?**
   ??? success "Solution"
       It balances a nonnegative matrix into a doubly stochastic matrix by iteratively normalizing rows and columns. It is used in optimization, economics, and attention mechanisms in deep learning.

10. **[Quotient Matrices] How do equitable partitions relate to the spectrum of a graph?**
    ??? success "Solution"
        If a graph has an equitable partition, the spectrum of the associated quotient matrix is a subset of the spectrum of the graph's adjacency matrix. This allows for reducing high-dimensional spectral problems to smaller, more manageable matrices.

## Chapter Summary

This chapter explores the mapping between graph topology and matrix algebra:

1. **Complexity Metrics**: Defined minimum rank and zero forcing number as measures of structural complexity in matrices.
2. **Topological Invariants**: Utilized Colin de Verdière parameters to express pure topological properties (like planarity) through eigenvalue multiplicities.
3. **Discrete Extremal Structures**: Examined Hadamard matrices and their role in reaching the maximum possible determinant for bounded entries.
4. **Stochastic Balance**: Investigated the exponent of primitive matrices and the Sinkhorn-Knopp scaling theory, providing tools for analyzing convergence in discrete systems.
