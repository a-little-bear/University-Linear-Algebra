# Chapter 58: Non-negative Matrix Factorization (NMF)

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices (Ch17) · Matrix Decompositions (Ch10) · SVD (Ch11) · Optimization (Ch25)

**Chapter Outline**: Definition of NMF → The "Parts-based" Motivation → Fundamental Differences between NMF and PCA/SVD (Interpretability) → Non-convexity and NP-hardness → Geometric Interpretation: Simplicial Cones → Core Algorithms: Multiplicative Update Rules (Lee & Seung) and Alternating Least Squares (ALS) → Regularized NMF (Sparsity and Manifold Constraints) → Applications: Image Feature Extraction, Topic Modeling, and Bioinformatics (Gene Clustering)

**Extension**: NMF is one of the most important dimensionality reduction tools in data mining; by enforcing non-negativity constraints, it achieves an automated deconstruction of complex data into its constituent parts, providing the mathematical foundation for computers to understand the logic of "parts forming a whole."

</div>

In traditional matrix decompositions like SVD, basis vectors and coefficients can be negative, which often lacks physical meaning when processing images or text. **Non-negative Matrix Factorization** (NMF) enforces that all components must be non-negative, ensuring that the results have a natural "parts-based" interpretability. For instance, a human face can be automatically decomposed into components like eyes, nose, and mouth. This chapter explores this decomposition technique, which presents both non-convex challenges and immense practical utility.

---

## 58.1 Definition and Motivation

!!! definition "Definition 58.1 (NMF)"
    Given a non-negative matrix $V \in \mathbb{R}^{m \times n}$ ($V \ge 0$), NMF seeks to find two non-negative matrices $W \in \mathbb{R}^{m \times k}$ and $H \in \mathbb{R}^{k \times n}$ such that:
    $$V \approx WH$$
    where $k \ll \min(m, n)$ is the number of bases (the rank).

!!! intuition "Interpretability: Local Representations"
    Because subtraction is not allowed, each column of $V$ must be a purely additive combination of the columns of $W$. This forces the columns of $W$ to represent **fundamental local features** within the data, while the rows of $H$ represent the **activation intensities** of those features.

---

## 58.2 Geometric Interpretation

!!! technique "Geometry: Simplicial Cones"
    The essence of NMF is finding a minimal **non-negative simplicial cone** that encloses all data points in $V$. Each column of data is restricted to the convex cone spanned by the columns of $W$. This contrasts sharply with PCA, which seeks a linear subspace along directions of maximum variance.

---

## 58.3 Core Algorithms

!!! algorithm "Algorithm 58.1 (Multiplicative Update Rules)"
    The classic algorithm proposed by Lee & Seung ensures that the objective function $\|V-WH\|_F$ is non-increasing through the following iterations:
    $$H_{aj} \leftarrow H_{aj} \frac{(W^T V)_{aj}}{(W^T WH)_{aj}}, \quad W_{ia} \leftarrow W_{ia} \frac{(VH^T)_{ia}}{(WHH^T)_{ia}}$$
    **Advantages**: Extremely simple to implement and automatically preserves non-negativity.
    **Disadvantages**: Converges to local minima and may get stuck in saddle points.

---

## 58.4 Variants and Extensions

!!! technique "Sparse NMF"
    By adding an $L_1$ regularization term $\lambda \|H\|_1$ to the objective, one can force the basis vectors to be more "pure" and sparse, further enhancing the discriminative power of the features.

---

## Exercises

1.  **[Basics] Prove: If $V = WH$ is an NMF, then for any positive diagonal matrix $D$, $V = (WD)(D^{-1}H)$ is also an NMF.**
    ??? success "Solution"
        $(WD)(D^{-1}H) = W(DD^{-1})H = WH = V$. Since both $D$ and its inverse have positive diagonals, non-negativity is preserved. This demonstrates the scale-indeterminacy of NMF.

2.  **[Contrast] What is the primary difference between NMF and PCA?**
    ??? success "Solution"
        PCA allows negative values and seeks to maximize variance (often leading to global, non-interpretable features). NMF enforces non-negativity and seeks interpretable local parts.

3.  **[NP-Hard] Why is NMF NP-hard to solve exactly?**
    ??? success "Solution"
        NMF is a non-convex optimization problem and is equivalent to finding a specific nesting of simplices, for which no polynomial-time global solver exists in general dimensions.

4.  **[Calculation] Using multiplicative updates, if an entry is initialized to 0, will it ever change?**
    ??? success "Solution"
        No. In multiplicative updates, 0 times any factor remains 0. Therefore, initialization is critical for NMF performance.

5.  **[Topic Modeling] In text mining, what do $W$ and $H$ represent?**
    ??? success "Solution"
        Each column of $W$ represents a "topic" (a distribution over words); each column of $H$ represents a document's "membership" or weight in each topic.

6.  **[Rank] Is the NMF rank $k$ always equal to the algebraic rank of $V$?**
    ??? success "Solution"
        Not necessarily. The non-negative rank is typically greater than or equal to the algebraic rank.

7.  **[Measure] Beyond the Frobenius norm, what other distance metric is common in NMF?**
    ??? success "Solution"
        **KL Divergence** (Kullback-Leibler) is frequently used, especially for count data and probability distributions.

8.  **[Uniqueness] Give a simple example where NMF is not unique.**
    ??? success "Solution"
        The identity matrix $I = I \cdot I$. Due to scaling and potential rotations that preserve non-negativity, the decomposition is generally not unique.

9.  **[Initialization] Why is all-zero initialization bad for NMF?**
    ??? success "Solution"
        Zeros are stationary points for the gradient; the iteration will never move and no features will be extracted. Randomized positive values or SVD-based initializations (like NNDSVD) are used instead.

****

??? success "Solution"
    

## Chapter Summary

NMF achieves a harmony between linear algebra and human cognition:

1.  **Philosophy of Constraint**: Proved that mathematical constraints (non-negativity) are not just limitations but sources of "meaning" (local features), establishing the benchmark for interpretability in reduction algorithms.
2.  **Challenge of Non-convexity**: Demonstrated how even simple linear products evolve into complex non-convex landscapes when signs are restricted, driving the development of alternating optimization algorithms.
3.  **Pattern Deconstruction**: From image recognition to genomics, NMF serves as a universal "pattern discovery" tool, proving that complex real-world data are often the superposition of a few pure, non-negative atomic components.
