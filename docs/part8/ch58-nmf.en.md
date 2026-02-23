# Chapter 58: Non-negative Matrix Factorization (NMF)

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices (Ch17) · Singular Value Decomposition (Ch11) · Matrix Analysis (Ch14)

**Chapter Outline**: Motivation for NMF (The Non-negative Essence of Data) → Definition of Non-negative Matrix Factorization → The Concept of Non-negative Rank $rank_{nn}(V)$ → Differences from SVD: From Global Orthogonality to Local Non-negativity → Core Trait: Part-based Representation → Algorithmic Implementation: Multiplicative Update Rules → Objective Functions: Frobenius Distance vs. KL Divergence → Applications: Text Mining (Topic Modeling), Facial Recognition (Feature Detection), Hyperspectral Unmixing, and Recommender Systems

**Extension**: NMF is the product of linear algebra's compromise with "interpretability"; by sacrificing orthogonality for non-negativity, it reveals how complex data is composed of simple, physically meaningful components through "addition" rather than "cancellation." It is a cornerstone of modern unsupervised learning.

</div>

In traditional SVD, basis vectors may contain negative values, which often lack intuitive meaning when processing image pixels or word counts. **Non-negative Matrix Factorization** (NMF) imposes a strong non-negativity constraint, forcing the system to synthesize data only through "addition." This "add-only" logic unexpectedly results in a **part-based** representation, allowing us to automatically extract meaningful "components" from raw data. This chapter explores this decomposition technique that balances computational challenge with explanatory power.

---

## 58.1 Definition and Non-negative Rank

!!! definition "Definition 58.1 (NMF)"
    Given a non-negative matrix $V \in \mathbb{R}^{m \times n}$ ($V \ge 0$), find non-negative matrices $W \in \mathbb{R}^{m \times k}$ and $H \in \mathbb{R}^{k \times n}$ such that:
    $$V \approx WH, \quad W, H \ge 0$$
    - $W$ is usually called the **Basis Matrix** (features).
    - $H$ is usually called the **Coefficient Matrix** (weights).

!!! definition "Definition 58.2 (Non-negative Rank)"
    The smallest dimension $k$ for which $V = WH$ holds exactly with $W, H \ge 0$ is the **non-negative rank** of $V$, denoted $rank_{nn}(V)$.
    **Property**: $rank(V) \le rank_{nn}(V)$. Computing the non-negative rank is NP-hard.

---

## 58.2 Core Trait: Part-based Representation

!!! technique "Explanation: Why extract "parts"?"
    Because entries of $W$ and $H$ are non-negative, a data point in $V$ can only be obtained through the **linear superposition** of basis vectors. To reconstruct complex patterns, the algorithm tends to make the basis vectors in $W$ sparse, representing distinct parts of an object (like eyes or a nose in a face). In SVD, basis vectors can cancel each other out (positive and negative), leading to overlapping global patterns that often lack physical meaning.

---

## 58.3 Algorithm: Multiplicative Update

!!! algorithm "Algorithm 58.1 (Lee-Seung Multiplicative Update)"
    For the Frobenius norm objective, the update rules are:
    $$H_{aj} \leftarrow H_{aj} \frac{(W^T V)_{aj}}{(W^T WH)_{aj}}, \quad W_{ia} \leftarrow W_{ia} \frac{(VH^T)_{ia}}{(WHH^T)_{ia}}$$
    **Trait**: As long as the initial values are positive, the update automatically maintains non-negativity and decreases the objective function at each step.

---

## Exercises

**1. [Basics] Find the non-negative rank of $V = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$.**

??? success "Solution"
    **Analysis:**
    1. The standard rank is 2.
    2. Since it is a diagonal matrix with 1s, we cannot synthesize these two orthogonal directions using only one non-negative column vector.
    **Conclusion**: The non-negative rank is 2. For the identity matrix, the non-negative rank equals the standard rank.

**2. [Comparison] Briefly state the main difference between NMF and SVD regarding basis vectors.**

??? success "Solution"
    **Comparison:**
    - **SVD**: Basis vectors are mutually orthogonal and contain both positive and negative values. Ideal for variance explanation but often lacks interpretability (e.g., negative pixel values).
    - **NMF**: Basis vectors are not necessarily orthogonal but must be non-negative. This produces "local" features with strong physical interpretability.

**3. [Calculation] Find the non-negative rank and an NMF factorization for $V = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.**

??? success "Solution"
    **Construction:**
    1. Standard rank is 1.
    2. $V = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \end{pmatrix}$.
    3. Factors are non-negative.
    **Conclusion**: The non-negative rank is 1.

**4. [Existence] Does every non-negative matrix have an NMF where $k = rank(V)$?**

??? success "Solution"
    **Conclusion: Not necessarily.**
    This is one of the deepest questions in NMF theory. There exist matrices (e.g., specific $5 \times 5$ constructions) with standard rank 3, but where every non-negative factorization requires $k \ge 4$ or $5$. This "rank gap" reflects the algebraic rigidity introduced by non-negativity constraints.

**5. [Objective] Name another common loss function for NMF besides the Frobenius norm.**

??? success "Solution"
    **Conclusion: Kullback-Leibler (KL) Divergence.**
    $$D(V || WH) = \sum (V_{ij} \log \frac{V_{ij}}{(WH)_{ij}} - V_{ij} + (WH)_{ij})$$
    This performs better when handling data with Poisson noise, such as count data or word frequencies in text.

**6. [Uniqueness] Is the NMF decomposition unique?**

??? success "Solution"
    **Conclusion: Usually not.**
    **Reasoning**: If $V = WH$, then for any diagonal matrix $D \succ 0$, $V = (WD)(D^{-1}H)$. There can also be complex rotational uncertainties. To obtain a unique solution, additional constraints such as **sparsity** are usually imposed.

**7. [Text Mining] In text mining, what does the basis matrix $W$ represent?**

??? success "Solution"
    **Explanation:**
    In a term-document matrix, each column of $W$ represents a **Topic**. Words with large values in a column are the core keywords for that topic. The coefficient matrix $H$ represents the distribution weights of each document over these topics.

**8. [Sparsity] Why do NMF results tend to be sparse?**

??? success "Solution"
    **Algebraic Intuition:**
    Since elements can only be added, to reconstruct a sparse data matrix, $W$ and $H$ must contain many zeros. The non-negativity constraint pushes the solution toward the "boundary" of the non-negative orthant (the axes), naturally inducing sparsity.

**9. [Complexity] Why is NMF a non-convex optimization problem?**

??? success "Solution"
    **Reasoning:**
    While the objective function is convex in $W$ (given $H$) and convex in $H$ (given $W$), it is **jointly non-convex** in $(W, H)$. This means algorithms may get stuck in local minima, and the choice of initialization is critical.

**10. [Application] Briefly describe NMF in a recommender system (e.g., movie ratings).**

??? success "Solution"
    1. $V$ is the user-movie rating matrix.
    2. $W$ represents user "latent preferences" (e.g., liking for sci-fi, action, romance).
    3. $H$ represents movie "latent attributes" (e.g., how much a movie belongs to those genres).
    4. Predicted ratings are the inner products of these non-negative vectors. Non-negativity ensures that preferences and attributes accumulate, which is intuitive.

## Chapter Summary

Non-negative matrix factorization represents a revolution of "interpretability" in linear algebra:

1.  **Reward of Constraints**: It proves that by imposing seemingly restrictive non-negativity, we gain structured decompositions with real-world physical meaning beyond pure mathematical optimization.
2.  **Part and Whole**: NMF successfully realizes the "part-based" decomposition, degrading complex signals into atomic components and establishing a new paradigm for feature extraction.
3.  **Computational Challenge**: The hardness of non-negative rank and joint non-convexity show that "structural interpretation" comes at a cost, driving the evolution of randomized and alternating optimization algorithms.
