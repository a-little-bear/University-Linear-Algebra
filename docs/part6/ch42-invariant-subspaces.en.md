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

****
??? success "Solution"
         - If $v \in \operatorname{Im}(T)$, then $Tv \in \operatorname{Im}(T)$ by the definition of the image. Thus $\operatorname{Im}(T)$ is invariant.

****
??? success "Solution"
    
e 0$. For $\mathcal{M}$ to be $T$-invariant, $Tv$ must be in $\mathcal{M}$. Thus $Tv = \lambda v$ for some scalar $\lambda$, which is the definition of an eigenvector.

****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ## Chapter Summary

This chapter explores the stability and structure of linear operators:


****: Defined invariant, hyperinvariant, and reducing subspaces as the atoms of operator decomposition.

****: Introduced canonical angles to quantify the distance between linear subspaces.

****: Decoded the Davis-Kahan and Wedin theorems, linking spectral gaps to the stability of eigenvectors and singular vectors.

****: Utilized spectral projections and Rosenblum's theorem to provide a rigorous framework for subspace analysis.
