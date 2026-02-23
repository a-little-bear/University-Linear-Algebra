# Chapter 58: Nonnegative Matrix Factorization

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · SVD (Ch11) · Nonnegative Matrices (Ch17) · Optimization (Ch25)

**Chapter Outline**: NMF Problem Definition → Contrast with SVD/PCA → Multiplicative Update Rules (Lee-Seung) → Alternating Nonnegative Least Squares (ANLS) → Uniqueness Conditions → Nonnegative Rank → Sparse and Regularized NMF → Applications

**Extension**: NMF is a fundamental tool in topic modeling (text mining), source separation (audio processing), and hyperspectral unmixing.

</div>

Nonnegative Matrix Factorization (NMF) seeks to decompose a nonnegative matrix $V \in \mathbb{R}_{\ge 0}^{m \times n}$ into two nonnegative matrices $W \in \mathbb{R}_{\ge 0}^{m \times r}$ and $H \in \mathbb{R}_{\ge 0}^{r \times n}$ such that $V \approx WH$. Unlike SVD, the nonnegativity constraint ensures a "parts-based" representation, where the data is reconstructed through additive combinations of basis vectors without subtractive cancellations.

---

## 58.1 Problem Formulation and Complexity

!!! definition "Definition 58.1 (NMF)"
    Given $V \ge 0$ and a rank $r$, NMF solves:
    $$\min_{W, H \ge 0} \|V - WH\|_F^2$$
    or other divergence measures like Kullback-Leibler (KL) divergence.

!!! theorem "Theorem 58.1 (Complexity)"
    Determining if a nonnegative matrix $V$ has a rank-$r$ nonnegative factorization $V=WH$ is **NP-hard** (Vavasis, 2009).

---

## Exercises

1. **[Parts-based Representation] Contrast NMF with PCA in terms of their representational logic.**
   ??? success "Solution"
       NMF constraints $W, H \ge 0$, forcing basis vectors to be additive components (e.g., local features). PCA allows negative weights, leading to "global" patterns where local features are created through subtractive cancellations of whole-body features.

2. **[Rank Relation] Prove that for any nonnegative matrix $V$, $\operatorname{rank}(V) \le \operatorname{rank}_+(V)$.**
   ??? success "Solution"
       Let $\operatorname{rank}_+(V) = r$. By definition, there exist $W \in \mathbb{R}^{m \times r}$ and $H \in \mathbb{R}^{r \times n}$ such that $V=WH$. Since $\operatorname{rank}(V) = \operatorname{rank}(WH) \le \min(\operatorname{rank}(W), \operatorname{rank}(H)) \le r$, the linear rank is always bounded above by the nonnegative rank.

3. **[Calculation] For $V = \begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}$, provide a rank-1 exact NMF.**
   ??? success "Solution"
       The columns are proportional. Let $W = \begin{pmatrix} 1 \\ 3 \end{pmatrix}$ and $H = \begin{pmatrix} 1 & 2 \end{pmatrix}$. Then $WH = V$ and $W, H \ge 0$.

4. **[Ambiguity] Discuss the scaling ambiguity in NMF.**
   ??? success "Solution"
       If $V=WH$, then $V=(WD)(D^{-1}H)$ for any diagonal matrix $D \succ 0$. This means we can scale the basis vectors as long as we compensate by inversely scaling the coefficients.

5. **[Algorithms] Describe the core idea of the Lee-Seung multiplicative update rules.**
   ??? success "Solution"
       The rule $H \leftarrow H \odot \frac{W^T V}{W^T WH}$ is derived from gradient descent with an adaptive step size that ensures the sign of the variables never changes, maintaining nonnegativity throughout the iteration.

6. **[Nonnegative Rank] Provide a $4 \times 4$ example where $\operatorname{rank}_+(V) > \operatorname{rank}(V)$.**
   ??? success "Solution"
       The matrix $V = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 1 & 0 & 1 & 0 \\ 0 & 1 & 0 & 1 \\ 0 & 0 & 1 & 1 \end{pmatrix}$ has linear rank 3 but nonnegative rank 4. This gap reflects the higher geometric complexity of the nonnegative cone compared to linear subspaces.

7. **[Sparsity] How does $L_1$ regularization affect NMF?**
   ??? success "Solution"
       It forces many entries in $H$ to zero, ensuring each document or image is represented by only a few "parts," which significantly enhances the interpretability of the topics or features.

8. **[Orthogonal NMF] Show that NMF with $HH^T = I$ is equivalent to K-means clustering.**
   ??? success "Solution"
       $H \ge 0$ and $HH^T = I$ force each column of $H$ to have exactly one non-zero entry. This partitions the data points into disjoint clusters, where $W$ columns act as cluster centroids.

9. **[Applications] Interpret $W$ and $H$ in topic modeling.**
   ??? success "Solution"
       In a term-document matrix $V$, $W$ contains topics (word distributions) and $H$ contains document topic-mixtures (the degree to which each document belongs to a topic).

10. **[Convexity] Is the NMF objective function convex?**
    ??? success "Solution"
        No. While it is convex in $W$ when $H$ is fixed, and vice versa, it is not jointly convex in $(W, H)$. This leads to multiple local minima and sensitivity to initialization.

## Chapter Summary

This chapter explores the theory and practice of Nonnegative Matrix Factorization:

1. **Representation Logic**: Established nonnegativity as the driver for additive, parts-based data representations.
2. **Algorithmic Evolution**: Detailed the multiplicative update rules and Alternating Nonnegative Least Squares (ANLS) frameworks.
3. **Structural Properties**: Analyzed the gap between linear rank and nonnegative rank, and the conditions for uniqueness (separability).
4. **Data Semantics**: Demonstrated NMF's superior interpretability in text mining, audio separation, and image analysis compared to traditional subspace methods.
