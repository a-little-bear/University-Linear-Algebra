# Chapter 40A: Permanents

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch03) · Doubly Stochastic Matrices (Ch31) · Combinatorics Basics

**Chapter Outline**: From Determinants to Permanents (Dropping the Sign) → Definition of the Permanent → Formal Contrast and Property Differences with Determinants → Computational Complexity: #P-completeness and Valiant’s Theorem → Ryser’s Inclusion-Exclusion Formula (Exponential-time Algorithm) → Permanents of Doubly Stochastic Matrices: The van der Waerden Theorem (Lower Bound) → Minc’s Conjecture (Upper Bound) → Applications: Counting Perfect Matchings in Bipartite Graphs, Statistical Mechanics, and Quantum Optics (Boson Sampling)

**Extension**: The permanent is the bridge connecting algebra to hard combinatorial problems; formally it differs from the determinant by only a "sign," but this tiny difference causes a massive chasm in computational complexity (polynomial vs. exponential). It is the ultimate subject for understanding the relationship between operators and counting.

</div>

If the determinant is the "volume" of linear algebra, then the **permanent** is its "counter." It is strikingly similar in form to the determinant but omits the sign of the permutations. This minor modification makes the permanent lose transpose-invariance (for non-square matrices) and the multiplicative property, yet it makes it the ultimate tool for solving perfect matching problems in graph theory. This chapter explores this field where algebraic form meets combinatorial hardness.

---

## 40A.1 Definition and Formal Comparison

!!! definition "Definition 40A.1 (Permanent)"
    The **permanent** of an $n \times n$ matrix $A$ is defined as:
    $$\operatorname{perm}(A) = \sum_{\sigma \in S_n} a_{1,\sigma(1)} a_{2,\sigma(2)} \cdots a_{n,\sigma(n)}$$
    Contrast: $\det(A) = \sum_{\sigma \in S_n} \operatorname{sgn}(\sigma) a_{1,\sigma(1)} \cdots a_{n,\sigma(n)}$.

!!! note "Property Differences"
    1.  **Non-multiplicative**: $\operatorname{perm}(AB) \neq \operatorname{perm}(A)\operatorname{perm}(B)$ (This is the permanent's most significant algebraic flaw).
    2.  **Multilinear**: Like the determinant, the permanent is linear with respect to each row (or column).
    3.  **Hardness**: Computing the permanent is **#P-complete**, meaning no polynomial-time algorithm is known.

---

## 40A.2 The van der Waerden Theorem and Bounds

!!! theorem "Theorem 40A.1 (van der Waerden Theorem)"
    For any $n \times n$ doubly stochastic matrix $P$:
    $$\operatorname{perm}(P) \ge \frac{n!}{n^n}$$
    Equality holds iff $P = J_n/n$ (the matrix with all entries equal to $1/n$).

---

## 40A.3 Computational Technique: Ryser’s Formula

!!! algorithm "Algorithm 40A.1 (Ryser’s Formula)"
    Based on the principle of inclusion-exclusion, the permanent can be expressed as:
    $$\operatorname{perm}(A) = (-1)^n \sum_{S \subseteq \{1,\ldots,n\}} (-1)^{|S|} \prod_{i=1}^n \left( \sum_{j \in S} a_{ij} \right)$$
    **Significance**: This reduces the complexity from $O(n!)$ to $O(n 2^n)$, serving as the standard method for medium-sized matrices.

---

## Exercises

**1. [Calculation] Find the permanent of $\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. Enumerate permutations: $\sigma_1 = (1, 2), \sigma_2 = (2, 1)$.
    2. Term 1: $a_{11}a_{22} = 1 \cdot 4 = 4$.
    3. Term 2: $a_{12}a_{21} = 2 \cdot 3 = 6$.
    4. Sum: $4 + 6 = 10$.
    **Conclusion**: $\operatorname{perm} = 10$. (Note: the determinant is -2).

**2. [Diagonal] For a diagonal matrix, how does the permanent relate to the determinant?**

??? success "Solution"
    **Conclusion:**
    For a diagonal matrix, **they are equal**.
    **Reasoning**: In both definitions, the only non-zero term corresponds to the identity permutation $\sigma=id$. Since $\operatorname{sgn}(id) = 1$, both results yield the product of the diagonal elements.

**3. [Graph Theory] If $A$ is the (0,1)-adjacency matrix of a bipartite graph, what is the meaning of $\operatorname{perm}(A)$?**

??? success "Solution"
    **Conclusion:**
    It represents the **total number of perfect matchings** in the bipartite graph.
    **Reasoning**: Each non-zero term in the permanent expansion corresponds to picking one edge for each vertex such that no vertices conflict. A product of 1s indicates a valid matching.

**4. [van der Waerden] Calculate the permanent of $J_2/2 = \begin{pmatrix} 0.5 & 0.5 \\ 0.5 & 0.5 \end{pmatrix}$ and verify the lower bound.**

??? success "Solution"
    **Calculation:**
    $\operatorname{perm} = 0.5 \cdot 0.5 + 0.5 \cdot 0.5 = 0.25 + 0.25 = 0.5$.
    **Lower Bound**: $n!/n^n = 2!/2^2 = 2/4 = 0.5$.
    **Conclusion**: Equality holds, consistent with the van der Waerden theorem.

**5. [Multilinearity] Prove: Multiplying the first row of $A$ by scalar $k$ multiplies the permanent by $k$.**

??? success "Solution"
    **Proof:**
    By definition, every product in the sum for $\operatorname{perm}(A)$ contains exactly one element from the first row ($a_{1,\sigma(1)}$). If $a_{1,j} \to k a_{1,j}$, every product gains a factor of $k$. By distributivity, the total sum is multiplied by $k$.

**6. [Upper Bound] What is Minc's Theorem?**

??? success "Solution"
    **Description:**
    For a (0,1)-matrix $A$ with row sums $r_i$, the permanent is bounded by:
    $$\operatorname{perm}(A) \le \prod_{i=1}^n (r_i !)^{1/r_i}$$
    This limits the maximum number of matchings a matrix with fixed sparsity can contain.

**7. [Powers] Prove for $A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$ that $\operatorname{perm}(A^2) \neq (\operatorname{perm} A)^2$.**

??? success "Solution"
    **Calculation:**
    1. $\operatorname{perm} A = 1+1 = 2$, so $(\operatorname{perm} A)^2 = 4$.
    2. $A^2 = \begin{pmatrix} 2 & 2 \\ 2 & 2 \end{pmatrix}$.
    3. $\operatorname{perm}(A^2) = 2 \cdot 2 + 2 \cdot 2 = 8$.
    **Conclusion**: $8 \neq 4$. This demonstrates the lack of a multiplicative law for permanents.

**8. [Skew-symmetry] Is the permanent of an odd-order skew-symmetric matrix necessarily zero?**

??? success "Solution"
    **Conclusion: No.**
    **Reasoning**: While the **determinant** of an odd-order skew-symmetric matrix is always zero, the permanent lacks alternating signs to cancel terms. For example, for $\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$, $\det = 1$ while $\operatorname{perm} = -1$.

**9. [Invariants] Under which matrix operations is the permanent invariant?**

??? success "Solution"
    **Conclusion:**
    1. **Permutation**: $A \to P A Q$ (where $P, Q$ are permutation matrices).
    2. **Transposition**: $\operatorname{perm}(A) = \operatorname{perm}(A^T)$ for square matrices.
    Note: It is **not** invariant under row additions (elimination).

**10. [Application] What is "Boson Sampling" in quantum computing?**

??? success "Solution"
    **Connection:**
    Boson Sampling is a computational task proving that the probability distribution of photons in an interference experiment is determined by the **permanents** of submatrices of the transition matrix. Since computing permanents is #P-complete, this task provides a path to demonstrating "Quantum Supremacy."

## Chapter Summary

The permanent is the link connecting algebra to hard combinatorial problems:

1.  **Cost of Symmetry**: By removing the alternating signs of the determinant, the permanent loses most desirable algebraic properties, leading directly to "catastrophic" computational complexity.
2.  **Combinatorial Counter**: By encoding constraints into (0,1)-matrices, the permanent becomes the master key for counting permutations and matchings, proving that algebraic forms can perfectly encapsulate discrete logic.
3.  **Boundaries of Calculation**: Ryser’s formula and the van der Waerden theorem define the operational limits of permanents from algorithmic and theoretical perspectives, making them central objects in modern complexity theory.
