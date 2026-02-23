# Chapter 40A: Permanents

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch3) · Combinatorial Structures (Ch65B) · Doubly Stochastic Matrices (Ch64A)

**Chapter Outline**: Definition of the Permanent → Contrast with Determinants → Permanents and Matchings → Properties of Permanents → Van der Waerden Conjecture → Egorychev-Falikman Theorem (Minimum of Permanent) → Upper Bounds (Minc's Conjecture) → Computational Complexity (#P-complete) → Permanent of Nonnegative Matrices

</div>

The **permanent** of a matrix is a polynomial in its entries that resembles the determinant but omits the sign factor of permutations. Despite this seemingly small difference, the permanent is vastly more difficult to compute than the determinant. While the determinant counts oriented volumes and can be computed in polynomial time, the permanent counts perfect matchings in bipartite graphs and is **#P-complete**. This chapter explores the combinatorial significance of the permanent and its extremum properties on the set of doubly stochastic matrices.

---

## 40A.1 Definition and Combinatorial Role

!!! definition "Definition 40A.1 (Permanent)"
    For an $n \times n$ matrix $A$, the permanent is defined as:
    $$\operatorname{perm}(A) = \sum_{\sigma \in S_n} \prod_{i=1}^n a_{i, \sigma(i)}$$

!!! theorem "Theorem 40A.1 (Egorychev-Falikman Theorem, 1981)"
    For any $n \times n$ doubly stochastic matrix $A$:
    $$\operatorname{perm}(A) \ge \frac{n!}{n^n}$$
    Equality holds if and only if $A = \frac{1}{n} J_n$ (the all-ones matrix scaled). This resolved the long-standing **van der Waerden conjecture**.

---

## Exercises

1. **[Fundamentals] Compute the permanent of $A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       $\operatorname{perm}(A) = 1 \cdot 1 + 1 \cdot 1 = 2$. Note $\det A = 0$.

2. **[Matchings] Show that for a (0,1)-matrix $A$, $\operatorname{perm}(A)$ counts the number of perfect matchings in the associated bipartite graph.**
   ??? success "Solution"
       Each permutation $\sigma$ represents a potential matching. The product $\prod a_{i, \sigma(i)}$ is 1 if all selected edges exist in the graph and 0 otherwise. Summing over all permutations yields the total count of perfect matchings.

3. **[Complexity] Why is computing the permanent harder than the determinant?**
   ??? success "Solution"
       The determinant is invariant under Gaussian elimination ($O(n^3)$). The permanent is not. Valiant (1979) proved that computing the permanent is #P-complete, meaning it is as hard as counting the solutions to an NP-complete problem.

4. **[Linearity] Is the permanent a multilinear function of the rows/columns?**
   ??? success "Solution"
       Yes. Like the determinant, the permanent is multilinear. However, it is **symmetric** rather than alternating: swapping two rows does not change the value of the permanent.

5. **[Van der Waerden] Use the Egorychev-Falikman theorem to bound the permanent of a $3 \times 3$ doubly stochastic matrix.**
   ??? success "Solution"
       $\operatorname{perm}(A) \ge \frac{3!}{3^3} = \frac{6}{27} \approx 0.222$.

6. **[Upper Bound] State Minc's Conjecture (Bregman's Theorem) for (0,1)-matrices.**
   ??? success "Solution"
       For a (0,1)-matrix $A$ with row sums $r_i$, $\operatorname{perm}(A) \le \prod (r_i!)^{1/r_i}$. This provides a sharp upper bound on the number of perfect matchings.

7. **[Positive Definiteness] Does $A \succeq 0$ imply $\operatorname{perm}(A) \ge \det A$?**
   ??? success "Solution"
       Yes. This is the **Schur permanent inequality**. For positive semi-definite matrices, the permanent is an upper bound on the determinant.

8. **[Invariance] What group of matrices $M$ satisfies $\operatorname{perm}(MAN) = \operatorname{perm}(A)$?**
   ??? success "Solution"
       Only the products of permutation matrices and diagonal matrices (monomial matrices) with appropriate scaling. This shows the permanent has far fewer symmetries than the determinant.

9. **[Glynn's Formula] Briefly describe a method to estimate the permanent.**
   ??? success "Solution"
       Since exact computation is hard, randomized algorithms are used. Glynn's formula provides an unbiased estimator based on sums over random vectors, while MCMC methods provide approximate counting.

10. **[Identity] Relate $\operatorname{perm}(J_n)$ to the number of permutations.**
    ??? success "Solution"
        $\operatorname{perm}(J_n) = n!$, as there are $n!$ permutations each contributing exactly 1 to the sum.

## Chapter Summary

This chapter explores the combinatorial analogue of the determinant:

1. **Matching Calculus**: Established the permanent as the counting function for bipartite matchings.
2. **Minimum Entropy**: Detailed the resolution of the van der Waerden conjecture, identifying the uniform distribution as the global minimizer.
3. **Hardness vs. Volume**: Analyzed the computational gap between permanents and determinants, highlighting the lack of elimination invariance.
4. **Spectral Constraints**: Investigated inequalities for PSD matrices, linking permanental bounds to spectral invariants.
