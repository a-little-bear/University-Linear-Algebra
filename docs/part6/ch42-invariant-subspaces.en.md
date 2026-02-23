# Chapter 42: Invariant Subspaces and Perturbations

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch4) · Linear Transformations (Ch5) · Eigenvalues (Ch6) · Jordan Form (Ch12) · Norms & Perturbation (Ch15)

**Chapter Outline**: Definition of Invariant Subspaces → Lattice of Subspaces → Hyperinvariant Subspaces → Reducing Subspaces → Canonical Angles between Subspaces → Davis-Kahan sin Θ Theorem → Spectral Projections → Wedin sin Θ Theorem → Stewart tan Θ Theorem → Rosenblum's Theorem → Invariant Subspace Problem

**Extension**: The Davis-Kahan theorem is the theoretical foundation for PCA perturbation analysis in modern statistics and machine learning.

</div>

Invariant subspaces are fundamental to understanding linear operators. By finding subspaces that are "closed" under $T$, we decompose the operator into simpler components. This chapter explores the algebraic structure of invariant subspaces and the geometric sensitivity of these spaces under matrix perturbations—a core topic in numerical linear algebra and data science.

---

## 42.1 Definitions and Basic Properties

!!! definition "Definition 42.1 ($T$-Invariant Subspace)"
    A subspace $\mathcal{M} \subseteq V$ is **$T$-invariant** if $T(\mathcal{M}) \subseteq \mathcal{M}$, i.e., $Tv \in \mathcal{M}$ for all $v \in \mathcal{M}$.

!!! theorem "Theorem 42.9 (Davis-Kahan sin Θ Theorem)"
    Let $A, 	ilde{A}$ be Hermitian. Let $\mathcal{V}, 	ilde{\mathcal{V}}$ be invariant subspaces associated with separated spectral components. If $\delta$ is the spectral gap, then
    $$\|\sin \Theta(\mathcal{V}, 	ilde{\mathcal{V}})\|_F \leq \frac{\|	ilde{A}-A\|_F}{\delta}.$$

---

## Exercises

1. **[Fundamentals] Prove that the kernel $\ker(T)$ and the image $\operatorname{Im}(T)$ are always $T$-invariant subspaces.**

   ??? success "Solution"
       - If $v \in \ker(T)$, then $Tv = 0 \in \ker(T)$. Thus $\ker(T)$ is invariant.
       - If $v \in \operatorname{Im}(T)$, then $Tv \in \operatorname{Im}(T)$ by the definition of the image. Thus $\operatorname{Im}(T)$ is invariant.

2. **[Eigenvectors] Show that every one-dimensional invariant subspace of $T$ is spanned by an eigenvector.**

   ??? success "Solution"
       Let $\mathcal{M} = \operatorname{span}\{v\}$ with $v 
e 0$. For $\mathcal{M}$ to be $T$-invariant, $Tv$ must be in $\mathcal{M}$. Thus $Tv = \lambda v$ for some scalar $\lambda$, which is the definition of an eigenvector.

3. **[Block Form] If $A$ has an invariant subspace $\mathcal{V}$ of dimension $k$, show that $A$ is similar to a block upper triangular matrix.**

   ??? success "Solution"
       Choose a basis $\{v_1, \dots, v_k\}$ for $\mathcal{V}$ and extend it to a basis for the whole space. In this basis, the first $k$ columns of the matrix representation will have zeros below the $k$-th row, resulting in a block upper triangular form $\begin{pmatrix} B & C \ 0 & D \end{pmatrix}$.

4. **[Hyperinvariant] Define a hyperinvariant subspace. Is every invariant subspace hyperinvariant?**

   ??? success "Solution"
       A subspace is hyperinvariant if it is invariant under every operator that commutes with $T$. Not every invariant subspace is hyperinvariant. For $T=I$, every subspace is invariant, but only $\{0\}$ and $V$ are hyperinvariant.

5. **[Reducing Subspaces] Prove: $\mathcal{M}$ reduces $T$ if and only if both $\mathcal{M}$ and $\mathcal{M}^\perp$ are $T$-invariant.**

   ??? success "Solution"
       This is the definition of a reducing subspace. In terms of projections, it is equivalent to $T$ commuting with the orthogonal projection onto $\mathcal{M}$.

6. **[Canonical Angles] Find the canonical angle between $\mathcal{F} = \operatorname{span}\{e_1\}$ and $\mathcal{G} = \operatorname{span}\{(1, 1)^T\}$ in $\mathbb{R}^2$.**

   ??? success "Solution"
       The cosine of the angle is the singular value of $P_{\mathcal{F}} P_{\mathcal{G}}$. $\cos 	heta = \frac{|(1,0) \cdot (1,1)|}{1 \cdot \sqrt{2}} = \frac{1}{\sqrt{2}}$. Thus $	heta = \pi/4$.

7. **[sin Θ Bound] Using Davis-Kahan, estimate the deviation of the principal eigenvector of $A = \operatorname{diag}(10, 1)$ if it is perturbed by $E = \begin{pmatrix} 0 & 0.1 \ 0.1 & 0 \end{pmatrix}$.**

   ??? success "Solution"
       Spectral gap $\delta = 10 - 1 = 9$. $\|E\|_2 = 0.1$. $\sin 	heta \le 0.1 / 9 \approx 0.011$.

8. **[Spectral Projection] Describe the spectral projection associated with an isolated eigenvalue $\lambda$.**

   ??? success "Solution"
       It is given by the contour integral $P = \frac{1}{2\pi i} \oint_\Gamma (zI - A)^{-1} dz$ where $\Gamma$ encloses only $\lambda$. For diagonalizable matrices, this is the sum of projections onto the eigenspaces of $\lambda$.

9. **[Wedin's Theorem] When do we use Wedin's sin Θ theorem instead of Davis-Kahan?**

   ??? success "Solution"
       Wedin's theorem is used for the perturbation of **singular subspaces** (SVD), whereas Davis-Kahan is for **eigenspaces** of Hermitian matrices.

10. **[Rosenblum's Theorem] Under what condition does the Sylvester equation $AX - XB = C$ have a unique solution?**

   ??? success "Solution"
        It has a unique solution if and only if $A$ and $B$ have no common eigenvalues: $\sigma(A) \cap \sigma(B) = \emptyset$.

## Chapter Summary

This chapter explores the stability and structure of linear operators:

1. **Subspace Dynamics**: Defined invariant, hyperinvariant, and reducing subspaces as the atoms of operator decomposition.
2. **Geometric Metrics**: Introduced canonical angles to quantify the distance between linear subspaces.
3. **Perturbation Robustness**: Decoded the Davis-Kahan and Wedin theorems, linking spectral gaps to the stability of eigenvectors and singular vectors.
4. **Analytical Tools**: Utilized spectral projections and Rosenblum's theorem to provide a rigorous framework for subspace analysis.
