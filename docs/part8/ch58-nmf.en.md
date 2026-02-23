# Chapter 58: Nonnegative Matrix Factorization

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · SVD (Ch11) · Nonnegative Matrices (Ch17) · Optimization (Ch25)

**Chapter Outline**: NMF Problem Definition → Contrast with SVD/PCA → Multiplicative Update Rules (Lee-Seung) → Alternating Nonnegative Least Squares (ANLS) → Uniqueness Conditions → Nonnegative Rank → Sparse and Regularized NMF → Applications

**Extension**: NMF is a fundamental tool in topic modeling (text mining), source separation (audio processing), and hyperspectral unmixing.

</div>

Nonnegative Matrix Factorization (NMF) seeks to decompose a nonnegative matrix $V \in \mathbb{R}_{\ge 0}^{m 	imes n}$ into two nonnegative matrices $W \in \mathbb{R}_{\ge 0}^{m 	imes r}$ and $H \in \mathbb{R}_{\ge 0}^{r 	imes n}$ such that $V \approx WH$. Unlike SVD, the nonnegativity constraint ensures a "parts-based" representation, where the data is reconstructed through additive combinations of basis vectors without subtractive cancellations.

---

## 58.1 Problem Formulation and Complexity

!!! definition "Definition 58.1 (NMF)"
    Given $V \ge 0$ and a rank $r$, NMF solves:
    $$\min_{W, H \ge 0} \|V - WH\|_F^2$$
    or other divergence measures like Kullback-Leibler (KL) divergence.

!!! theorem "Theorem 58.1 (Complexity)"
    Determining if a nonnegative matrix $V$ has a rank-$r$ nonnegative factorization $V=WH$ is **NP-hard** (Vavasis, 2009). This contrasts with the polynomial-time solvability of standard low-rank approximation via SVD.

---

## Exercises

1. **[Parts-based Representation] Contrast NMF with PCA in terms of their representational logic.**
   ??? success "Solution"
       NMF constraints $W, H \ge 0$, which forces the basis vectors to be purely additive components (e.g., facial features like eyes or noses). PCA allows negative weights, leading to "eigenfaces" that represent global patterns and require subtractive cancellations to represent local features.

2. **[Rank Relation] Prove that for any nonnegative matrix $V$, $\operatorname{rank}(V) \le \operatorname{rank}_+(V)$.**
   ??? success "Solution"
       Let $\operatorname{rank}_+(V) = r$. By definition, there exist $W \in \mathbb{R}^{m 	imes r}$ and $H \in \mathbb{R}^{r 	imes n}$ such that $V=WH$. From basic linear algebra, $\operatorname{rank}(V) = \operatorname{rank}(WH) \le \min(\operatorname{rank}(W), \operatorname{rank}(H)) \le r$. Thus, the nonnegative rank is always an upper bound on the linear rank.

3. **[Calculation] For $V = \begin{pmatrix} 1 & 2 \ 3 & 6 \end{pmatrix}$, provide a rank-1 exact NMF.**
   ??? success "Solution"
       The columns are proportional. Let $W = \begin{pmatrix} 1 \ 3 \end{pmatrix}$ and $H = \begin{pmatrix} 1 & 2 \end{pmatrix}$. $WH = V$ and $W, H \ge 0$.

4. **[Uniqueness] Discuss the scaling and permutation ambiguities in NMF.**
   ??? success "Solution"
       If $V=WH$, then $V=(WD)(D^{-1}H)$ for any diagonal matrix $D \succ 0$. Similarly, $V=(WP)(P^T H)$ for any permutation matrix $P$. These represent the inherent non-uniqueness of the factorization.

5. **[Algorithms] Describe the core idea of the Lee-Seung multiplicative update rules.**
   ??? success "Solution"
       The rule $H \leftarrow H \odot \frac{W^T V}{W^T WH}$ is a gradient descent variant with an adaptive step size that cancels the additive update terms, ensuring that variables remain nonnegative throughout the iteration provided the initial state is positive.

6. **[Nonnegative Rank] Provide an example where $\operatorname{rank}_+(V) > \operatorname{rank}(V)$.**
   ??? success "Solution"
       The $4 	imes 4$ "Euclidean distance" pattern matrix $V = \begin{pmatrix} 1 & 1 & 0 & 0 \ 1 & 0 & 1 & 0 \ 0 & 1 & 0 & 1 \ 0 & 0 & 1 & 1 \end{pmatrix}$ has linear rank 3 but nonnegative rank 4. This gap reflects the geometric complexity of the nonnegative cone.

7. **[Sparsity] How does $L_1$ regularization on $H$ affect the NMF result?**
   ??? success "Solution"
       It promotes sparsity in the coefficient matrix $H$, meaning each data point is represented as a combination of only a few basis vectors. This enhances the interpretability of the parts-based decomposition.

8. **[Orthogonal NMF] Show that NMF with the constraint $HH^T = I$ is equivalent to K-means clustering.**
   ??? success "Solution"
       The combination of $H \ge 0$ and $HH^T = I$ forces each column of $H$ to have exactly one non-zero entry. This effectively assigns each data point to exactly one basis vector (centroid), which is the definition of hard clustering.

9. **[Applications] Interpret $W$ and $H$ in the context of topic modeling.**
   ??? success "Solution"
       In a term-document matrix $V$, $W$ represents the "term-topic" matrix (columns are topics defined by word distributions), and $H$ represents the "topic-document" matrix (columns are documents defined by topic mixtures).

10. **[Stationary Points] Do the multiplicative updates guarantee convergence to a global minimum?**
    ??? success "Solution"
        No. The NMF objective is non-convex in both $W$ and $H$ simultaneously. The updates only guarantee that the objective function is non-increasing and that the limit point is a stationary point (satisfying KKT conditions).

## Chapter Summary

This chapter explores the theory and practice of Nonnegative Matrix Factorization:

1. **Representation Logic**: Established nonnegativity as the driver for additive, parts-based data representations.
2. **Algorithmic Evolution**: Detailed the multiplicative update rules and Alternating Nonnegative Least Squares (ANLS) frameworks.
3. **Structural Properties**: Analyzed the gap between linear rank and nonnegative rank, and the conditions for uniqueness (separability).
4. **Data Semantics**: Demonstrated NMF's superior interpretability in text mining, audio separation, and image analysis compared to traditional subspace methods.
