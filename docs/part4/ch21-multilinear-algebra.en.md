# Chapter 21: Multilinear Algebra and Tensors

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Linear Transformations (Ch05) · Dual Spaces (Ch13A)

**Chapter Outline**: From Linear to Multilinear → Definition of Multilinear Mappings → Axiomatic Construction and Universal Property of the Tensor Product ($\otimes$) → Component Representation and Index Notation → Dimension of Tensor Spaces $V \otimes W$ → Tensor Rank and CP Decomposition → Tucker Decomposition → Correspondence between Operators and Tensors → Applications: Multivariate Statistics, Quantum Entanglement, and High-dimensional Data Compression

**Extension**: Tensors are the natural generalization of linear algebra to high-dimensional spaces; they break the "row and column" binary limit of matrices, introducing algebraic structures with multi-axis coupling. They are the underlying mathematical skeleton for modern Deep Learning (TensorFlow) and General Relativity.

</div>

Linear algebra primarily studies vectors (1st-order tensors) and matrices (2nd-order tensors). However, in modern physics and data science, we need to handle quantities with multiple indices (3rd-order or higher). **Multilinear Algebra** extends linear transformation theory to arbitrary dimensions through the **Tensor Product**. This chapter introduces the universal property of tensors and their central role in describing complex high-dimensional interactions.

---

## 21.1 Definition of the Tensor Product

!!! definition "Definition 21.1 (Tensor Product)"
    The **tensor product** $V \otimes W$ of vector spaces $V$ and $W$ is a new vector space equipped with a bilinear map $\otimes: V \times W \to V \otimes W$ satisfying the **Universal Property**: every bilinear map from $V \times W$ to a space $U$ induces a unique linear map from $V \otimes W$ to $U$.

!!! note "Dimension Property"
    $\dim(V \otimes W) = \dim(V) \cdot \dim(W)$.

---

## 21.2 Order and Rank of Tensors

!!! definition "Definition 21.2 (Order and Rank)"
    1.  **Order (Mode)**: The number of indices. A matrix is a 2nd-order tensor.
    2.  **Tensor Rank**: The minimum number of pure tensors (Rank-1 tensors, $v_1 \otimes v_2 \otimes \cdots$) required to sum to the given tensor.
    **Warning**: Unlike matrix rank, calculating the rank of higher-order tensors is NP-hard.

---

## 21.3 Tensor Decompositions

!!! technique "Core Decomposition Models"
    1.  **CP Decomposition (CANDECOMP/PARAFAC)**: Expresses a tensor as a sum of Rank-1 tensors.
    2.  **Tucker Decomposition**: Expresses a tensor as a product of a core tensor and factor matrices for each mode (analogous to a higher-order SVD).

---

## Exercises

**1. [Basics] If $\dim V = 2$ and $\dim W = 3$, find the dimension of $V \otimes W$ and list its basis vectors.**

??? success "Solution"
    **Calculation:**
    1. $\dim(V \otimes W) = 2 \cdot 3 = 6$.
    2. Let the basis of $V$ be $\{e_1, e_2\}$ and $W$ be $\{f_1, f_2, f_3\}$.
    **Basis**: $\{e_1 \otimes f_1, e_1 \otimes f_2, e_1 \otimes f_3, e_2 \otimes f_1, e_2 \otimes f_2, e_2 \otimes f_3\}$.

**2. [Property] Prove: $(v_1 + v_2) \otimes w = v_1 \otimes w + v_2 \otimes w$.**

??? success "Solution"
    **Reasoning:**
    This is directly guaranteed by the **bilinearity** requirement in the tensor product definition. The operator $\otimes$ is linear with respect to each of its input arguments.

**3. [Matrix View] Prove: Any $m \times n$ matrix can be viewed as a tensor in $V \otimes W^*$.**

??? success "Solution"
    **Algebraic Mapping:**
    1. Linear transformations $T: W \to V$ form the space $\mathcal{L}(W, V)$.
    2. There is a canonical isomorphism $\mathcal{L}(W, V) \cong V \otimes W^*$.
    3. A pure tensor $v \otimes f$ corresponds to the operator $T(w) = f(w)v$, which is a rank-1 matrix.
    **Conclusion**: Matrix addition and scaling are perfectly consistent with tensor properties.

**4. [Calculation] Given $v = (1, 2)^T$ and $w = (0, 3)^T$, write the tensor product $v \otimes w$ as a matrix.**

??? success "Solution"
    **Steps:**
    In coordinate spaces, the tensor product corresponds to the outer product $v w^T$.
    $v \otimes w = \begin{pmatrix} 1 \cdot 0 & 1 \cdot 3 \\ 2 \cdot 0 & 2 \cdot 3 \end{pmatrix} = \begin{pmatrix} 0 & 3 \\ 0 & 6 \end{pmatrix}$.

**5. [Rank] Determine the rank of the tensor $\begin{pmatrix} 0 & 3 \\ 0 & 6 \end{pmatrix}$.**

??? success "Solution"
    **Conclusion: Rank 1.**
    **Reasoning**: This matrix is explicitly generated as the tensor product of two vectors. By definition, a quantity representable by a single pure tensor has rank 1.

**6. [Higher Order] How many entries are in a $3 \times 3 \times 3$ tensor?**

??? success "Solution"
    **Calculation:**
    $3 \cdot 3 \cdot 3 = 27$.
    The number of elements grows exponentially with the order, a phenomenon known as the "Curse of Dimensionality."

**7. [Index Notation] Write the tensor contraction $C_{ik} = \sum_j A_{ij} B_{jk}$ using Einstein summation convention.**

??? success "Solution"
    **Conclusion: $C_{ik} = A_{ij} B^{jk}$.**
    In Einstein notation, indices repeated in upper and lower positions imply summation, significantly simplifying multilinear expressions.

**8. [CP Decomposition] Briefly state the role of CP decomposition in data analysis.**

??? success "Solution"
    **Explanation:**
    CP decomposition aims to find the underlying independent components of data. For instance, in a (Time $\times$ Frequency $\times$ Space) signal, CP decomposition breaks down the complex tensor into "characteristic timelines," "spectral signatures," and "spatial maps," enabling feature extraction across modes.

**9. [Uniqueness] How does CP decomposition of high-order tensors differ from matrix SVD regarding uniqueness?**

??? success "Solution"
    **Key Advantage:**
    Under mild conditions, the CP decomposition of a high-order tensor is **unique** (up to permutation and scaling of factors). Unlike matrix factorizations, which have infinite rotational possibilities without constraints, tensor decompositions enable automatic signal separation (blind source separation) without human intervention.

**10. [Application] Why can't quantum entangled states be simply factored?**

??? success "Solution"
    **Physical Insight:**
    1. A multi-particle state $|\psi\rangle$ lives in the tensor product of Hilbert spaces $H_1 \otimes H_2$.
    2. If $|\psi\rangle$ can be written as $v_1 \otimes v_2$ (rank-1 tensor), it is a separable state (no entanglement).
    3. **Entangled states** are tensors with rank $> 1$. Because they cannot be decomposed into independent components, measuring one particle instantly affects the description of the other—a direct physical manifestation of tensor superposition.

## Chapter Summary

Multilinear algebra pushes linear thinking into all-dimensional space:

1.  **Universality of Construction**: The tensor product establishes the only valid way for multi-axis data interaction through its universal property, providing a template for describing stress, electromagnetic fields, and quantum states.
2.  **Evolution of Rank**: Moving from matrix to tensor rank reveals the complexity of "structural simplification" in high-dimensional spaces while introducing desirable traits like decomposition uniqueness.
3.  **Weaving Information**: Tensor theory proves that complex system properties are often hidden in the coupling between indices, making it an indispensable algebraic foundation for modern big data analysis and high-energy physics.
