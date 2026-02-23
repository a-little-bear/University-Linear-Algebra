# Chapter 21: Multilinear Algebra and Tensors

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch4) · Dual Spaces (Ch13a) · Kronecker Product (Ch19) · Exterior Algebra (Ch49)

**Chapter Outline**: Dual Spaces and Multilinear Forms → Tensor Product of Vector Spaces $V \otimes W$ → General Tensors $T_{q}^p$ → Covariant and Contravariant Indices → Tensor Contraction → Inner Products on Tensor Spaces → Symmetric and Antisymmetric Tensors → Tensor Rank and Decomposition (CP/Tucker) → Applications in Physics and Data Science

**Extension**: Multilinear algebra is the mathematical language of general relativity (curvature tensors) and modern high-dimensional data analysis (tensor networks).

</div>

Linear algebra typically deals with mappings between two spaces (matrices). **Multilinear algebra** extends this to relationships involving multiple spaces simultaneously. A **tensor** is the universal object that represents multilinear maps. By introducing the **tensor product** $\otimes$, we can describe interactions that cannot be captured by simple vectors, such as the coupling between particles or the multi-way correlations in big data. This chapter formalizes the "index calculus" of tensors and introduces the concept of **tensor rank**, which is significantly more complex than matrix rank.

---

## 21.1 Dual Spaces and Tensors

!!! definition "Definition 21.1 (Multilinear Map)"
    A map $f: V_1 \times \dots \times V_k \to W$ is **multilinear** if it is linear in each argument separately when the others are fixed.

!!! definition "Definition 21.2 (Tensor Product)"
    The tensor product $V \otimes W$ is the unique vector space (up to isomorphism) that converts every multilinear map $V \times W \to U$ into a unique linear map $V \otimes W \to U$.

---

## Exercises

1. **[Fundamentals] Compute the dimension of $V \otimes W$ if $\dim V = m$ and $\dim W = n$.**
   ??? success "Solution"
       $\dim(V \otimes W) = mn$. A basis is given by $\{v_i \otimes w_j\}$ where $\{v_i\}$ and $\{w_j\}$ are bases for $V$ and $W$ respectively.

2. **[Contraction] Define the contraction of a $(1, 1)$-tensor $T \in V \otimes V^*$.**
   ??? success "Solution"
       The contraction is the linear map that sends $v \otimes \phi$ to the scalar $\phi(v)$. For a matrix (representing a $(1,1)$-tensor), the contraction is exactly the **trace**.

3. **[Dual Space] Show that $(V \otimes W)^* \cong V^* \otimes W^*$.**
   ??? success "Solution"
       This isomorphism identifies a multilinear form on $V \times W$ with a tensor in the product of the dual spaces.

4. **[Tensor Rank] Define the rank of a tensor $T \in V \otimes W \otimes U$.**
   ??? success "Solution"
       The smallest $r$ such that $T = \sum_{i=1}^r v_i \otimes w_i \otimes u_i$. Unlike matrices, the rank of a 3-way tensor can exceed its dimensions and is NP-hard to compute.

5. **[Symmetry] Distinguish between the symmetric space $S^k(V)$ and the alternating space $\Lambda^k(V)$.**
   ??? success "Solution"
       $S^k(V)$ consists of tensors invariant under index permutations (e.g., $u \otimes v + v \otimes u$). $\Lambda^k(V)$ consists of tensors that change sign under swaps (e.g., $u \otimes v - v \otimes u$).

6. **[Bilinear Forms] Relate $V^* \otimes V^*$ to the space of bilinear forms on $V$.**
   ??? success "Solution"
       Every bilinear form $B(u, v)$ can be uniquely written as a tensor $\sum b_{ij} \epsilon^i \otimes \epsilon^j$ where $\epsilon^i$ are dual basis elements.

7. **[Universal Property] Explain the importance of the universal property of the tensor product.**
   ??? success "Solution"
       It allows for the systematic "linearization" of multilinear problems. Any operation involving products of variables can be handled by standard linear algebra in the larger tensor space.

8. **[Physics] In physics, what is a contravariant vs. a covariant tensor?**
   ??? success "Solution"
       Contravariant tensors (vectors) transform like the basis change matrix $P$. Covariant tensors (dual vectors/forms) transform like $(P^{-1})^T$. Tensors of type $(p, q)$ have $p$ contravariant and $q$ covariant indices.

9. **[Decomposition] What is the CP decomposition (CANDECOMP/PARAFAC)?**
   ??? success "Solution"
       A generalization of SVD to tensors where a tensor is approximated by a sum of rank-1 tensors. It is used for blind source separation and data mining.

10. **[Inner Product] Define the induced inner product on $V \otimes W$.**
    ??? success "Solution"
        $\langle v_1 \otimes w_1, v_2 \otimes w_2 \rangle = \langle v_1, v_2 \rangle_V \langle w_1, w_2 \rangle_W$. This makes the tensor product of Hilbert spaces a Hilbert space.

## Chapter Summary

This chapter establishes the calculus of high-dimensional interactions:

1. **Multilinear Grounding**: Defined tensors as the universal objects representing linear interactions across multiple spaces.
2. **Structural Composition**: Developed the tensor product as the mechanism for expanding the state space of linear systems.
3. **Index Calculus**: Formulated contractions and trace operators as the fundamental reductions of multilinear forms.
4. **Rank and Complexity**: Analyzed tensor decomposition, highlighting the transition from the simple rank of matrices to the complex landscapes of high-order tensors.
