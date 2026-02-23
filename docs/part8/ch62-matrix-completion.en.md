# Chapter 62: Matrix Completion

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · SVD (Ch11) · Optimization (Ch25) · Graph Theory (Ch27)

**Chapter Outline**: General Framework → Positive Definite Completion (Chordal Graphs) → Low-rank Completion (Netflix Problem) → Nuclear Norm Minimization → Exact Recovery Conditions → Euclidean Distance Matrix Completion → Algorithms → Noisy Completion → Matrix Sensing

**Extension**: Low-rank matrix completion is the mathematical foundation for recommendation systems (Netflix, Amazon); PD completion is vital in spatial statistics and sparse covariance estimation.

</div>

The matrix completion problem asks to recover the missing entries of a matrix given a subset of its entries, such that the completed matrix satisfies a specific property (e.g., positive definiteness, low rank, or distance metric). This problem intersects graph theory, convex optimization, and random matrix theory.

---

## 62.1 Framework and Criteria

!!! definition "Definition 62.2 (Matrix Completion Problem)"
    Given a partial matrix $A_\Omega$ with known entries on index set $\Omega$, find $X$ such that $X_{ij} = A_{ij}$ for all $(i,j) \in \Omega$ and $X \in \mathcal{P}$ (where $\mathcal{P}$ is a property class).

!!! theorem "Theorem 62.1 (Grone's Theorem, 1984)"
    Every positive definite partial matrix with a given sparsity pattern graph $G$ can be completed to a positive definite matrix if and only if $G$ is a **chordal graph**.

---

## Exercises

1. **[Concept] Why can the Netflix recommendation problem be modeled as a low-rank matrix completion problem?**
   ??? success "Solution"
       The user-movie rating matrix is extremely sparse. Since users' preferences and movie characteristics are typically influenced by a few latent factors (e.g., genre, director), the full matrix is approximately low-rank. Matrix completion seeks to predict missing ratings by filling the matrix under this low-rank constraint.

2. **[PD Completion] Explain the significance of "Chordal Graphs" in Grone's Theorem.**
   ??? success "Solution"
       Chordal graphs ensure that local positive definiteness (on all cliques) implies global consistency. If the graph has non-chordal cycles of length $\ge 4$, one can construct cases where all specified sub-blocks are positive definite, yet no global positive definite completion exists due to non-local constraints.

3. **[Nuclear Norm] Define the nuclear norm and its relation to the rank of a matrix.**
   ??? success "Solution"
       The nuclear norm $\|X\|_*$ is the sum of the singular values. It is the convex envelope (the tightest convex lower bound) of the rank function on the unit spectral norm ball. In optimization, we minimize the nuclear norm as a convex proxy for the non-convex rank minimization problem.

4. **[Incoherence] Why is "Incoherence" critical for exact recovery in matrix completion?**
   ??? success "Solution"
       If a matrix is highly coherent (e.g., $e_1 e_1^T$), its information is concentrated in a single entry. Random sampling might miss this entry entirely, making recovery impossible. Incoherence ensures the information is spread across all entries, allowing for global inference from local samples.

5. **[Calculation] For a $1000 \times 1000$ matrix of rank 10, calculate the number of degrees of freedom. How many samples are typically required?**
   ??? success "Solution"
       Degrees of freedom: $r(m+n-r) = 10(1000+1000-10) = 19,900$. Information theory and the Candes-Recht theorem suggest $O(rn \log n)$ samples are needed for exact recovery with high probability.

6. **[Algorithm] Outline the Singular Value Thresholding (SVT) algorithm.**
   ??? success "Solution"
       1. **Project**: Compute the residual on the known entries. 2. **Update**: Add a step based on the residual to the current estimate. 3. **Threshold**: Apply the soft-thresholding operator to the singular values (SVD followed by shrinking). 4. **Repeat** until convergence.

7. **[EDM] How is Euclidean Distance Matrix (EDM) completion related to Semidefinite Programming (SDP)?**
   ??? success "Solution"
       Using the relation $D_{ij} = G_{ii} + G_{jj} - 2G_{ij}$, completing a distance matrix $D$ is equivalent to completing a Gram matrix $G \succeq 0$ under linear constraints. This is a standard SDP problem used in sensor network localization.

8. **[Noise] What does the error bound $O(\sigma\sqrt{r/p})$ imply in noisy matrix completion?**
   ??? success "Solution"
       It shows that the recovery error (RMSE) scales linearly with the noise level $\sigma$ and decreases with the square root of the sampling rate $p$. This demonstrates that nuclear norm minimization is stable under noisy observations.

9. **[RIP] Contrast Matrix RIP with the entry-wise incoherence condition.**
   ??? success "Solution"
       Matrix Restricted Isometry Property (RIP) is a property of the sensing operator (mapping matrices to measurements) that ensures it preserves the energy of low-rank matrices. Entry-wise sampling does not satisfy RIP in the general sense, which is why incoherence is required as a separate assumption.

10. **[Max-Det] Describe the "Maximum Determinant Completion" for PD matrices.**
    ??? success "Solution"
        Among all positive definite completions of a chordal partial matrix, there is a unique one that maximizes the determinant. This completion is characterized by having an inverse matrix whose entries are zero on the missing index set, linking it to Markov random fields.

## Chapter Summary

This chapter details the mathematical techniques for inferring global structure from local observations:

1. **Topological Constraints**: Established chordality as the graph-theoretic requirement for guaranteed positive definite completion.
2. **Convex Relaxations**: Demonstrated nuclear norm minimization as the primary tool for low-rank recovery.
3. **Recovery Theory**: Analyzed the interplay between incoherence and sampling density in ensuring exact reconstruction.
4. **Practical Solvers**: Introduced SVT and alternating minimization algorithms capable of handling massive sparse datasets.
