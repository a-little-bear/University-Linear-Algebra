# Chapter 40A: Permanents

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch03) · Doubly Stochastic Matrices (Ch31)

**Chapter Outline**: Definition of the Permanent → Formal Comparison with the Determinant → Computational Complexity (#P-completeness) → Ryser’s Formula → Permanents of Doubly Stochastic Matrices → van der Waerden Conjecture (Theorem) → Minc’s Conjecture (Upper Bounds) → Applications: Counting Perfect Matchings in Bipartite Graphs and the Dimer Problem in Statistical Mechanics

**Extension**: The permanent is one of the most significant matrix functions in combinatorics; although it differs from the determinant by only a single "sign," this difference leads to a massive computational gap between polynomial time and exponential time.

</div>

If the determinant represents the "volume" of linear algebra, then the **Permanent** is its "counter." Formally striking in its similarity to the determinant, the permanent omits the alternating sign factor of permutations. This minor change renders the permanent devoid of properties like multiplicative consistency or transpose invariance (for non-square matrices), yet transforms it into the ultimate tool for solving coloring and matching problems in graph theory.

---

## 40A.1 Definition and Formal Comparison

!!! definition "Definition 40A.1 (The Permanent)"
    The **permanent** of an $n \times n$ matrix $A$ is defined as:
    $$\operatorname{perm}(A) = \sum_{\sigma \in S_n} a_{1,\sigma(1)} a_{2,\sigma(2)} \cdots a_{n,\sigma(n)}$$
    Contrast with: $\det(A) = \sum \operatorname{sgn}(\sigma) a_{1,\sigma(1)} \cdots a_{n,\sigma(n)}$.

!!! note "Computational Complexity"
    While the determinant can be solved in $O(n^3)$ time, computing the permanent is proven to be **#P-complete** (Valiant's Theorem). This implies that no polynomial-time algorithm is currently known.

---

## 40A.2 The van der Waerden Theorem

!!! theorem "Theorem 40A.1 (van der Waerden Theorem)"
    For an $n \times n$ doubly stochastic matrix $P$ (non-negative with row and column sums equal to 1):
    $$\operatorname{perm}(P) \ge \frac{n!}{n^n}$$
    Equality holds if and only if $P = J_n/n$ (the all-ones matrix divided by $n$). This result reveals the lower bound for the permanent of doubly stochastic matrices.

---

## 40A.3 Computational Formula: Ryser’s Formula

!!! algorithm "Algorithm 40A.1 (Ryser’s Formula)"
    To compute the permanent faster than enumerating all $n!$ permutations, Ryser provided a formula based on the principle of inclusion-exclusion:
    $$\operatorname{perm}(A) = (-1)^n \sum_{S \subseteq \{1,\ldots,n\}} (-1)^{|S|} \prod_{i=1}^n \left( \sum_{j \in S} a_{ij} \right)$$
    **Complexity**: $O(n 2^n)$. This is still exponential but much faster than naive summation for medium-sized matrices.

---

## Exercises


****
??? success "Solution"
     $1(4) + 2(3) = 4 + 6 = 10$. (Note: the determinant is -2).


****
??? success "Solution"
     It is the product of the diagonal entries. In this case, the permanent equals the determinant.


****
??? success "Solution"
     It represents the total number of **perfect matchings** in that bipartite graph.


****
??? success "Solution"
     $\operatorname{perm} = 0.5 \cdot 0.5 + 0.5 \cdot 0.5 = 0.5$. The bound is $2!/2^2 = 2/4 = 0.5$. Equality holds.


****
??? success "Solution"
     Each term in the sum contains exactly one entry from each row. Thus, the sum maintains linearity relative to any single row vector.


****
??? success "Solution"
     $\operatorname{perm}(A) \le \prod (r_i !)^{1/r_i}$.


****
??? success "Solution"
     Generally no. This is the most significant algebraic deficiency of the permanent compared to the determinant.


****
??? success "Solution"
     No. For example, $\operatorname{perm} \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = -1$.


****
??? success "Solution"
     $\operatorname{perm}(PAQ) = \operatorname{perm}(A)$ for any permutation matrices $P, Q$.

****
??? success "Solution"
    ## Chapter Summary

The permanent serves as the bridge between algebra and combinatorial hardness:


****: By removing the alternating signs of the determinant, the permanent loses most of its elegant algebraic properties (like the multiplicative rule), directly resulting in its "catastrophic" computational complexity.

****: By encoding combinatorial constraints into (0,1)-matrices, the permanent becomes the master key for solving counting problems involving permutations, matchings, and colorings.

****: The van der Waerden theorem demonstrates how the "evenness" of matrix entries (doubly stochasticity) strictly limits the richness of its topological structures (the size of the permanent).
