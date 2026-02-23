# Chapter 21: Multilinear Algebra and Tensors

<div class="context-flow" markdown>

**Prerequisites**: Kronecker Product (Ch19) · Vector Spaces (Ch4) · Dual Spaces (Ch13a)

**Chapter Outline**: Definition of Multilinear Mappings → Tensor Product $V \otimes W$ → Tensor Rank and Decompositions (CP Decomposition, Tucker Decomposition) → Mode-$n$ Product → Tensor Contraction → Arithmetic Tensors (Symmetric and Antisymmetric) → Applications (High-dimensional Data Compression, Quantum Entanglement, Chemometrics)

**Extension**: Tensors are the natural higher-order generalization of matrices; multilinear algebra pushes matrix-based binary relations into the complex realm of multivariate relations.

</div>

Multilinear algebra studies the linear interactions between multiple vector spaces. A tensor is more than just a multi-dimensional array; it is the algebraic representation of a multilinear mapping. In the era of deep learning and quantum computing, tensor decomposition has become a core technology for tackling the "curse of dimensionality."

---

## 21.1 Tensor Products and Decompositions

!!! definition "Definition 21.1 (Tensor Product Space)"
    The tensor product $V \otimes W$ of vector spaces $V$ and $W$ is the vector space that satisfies the universal property: any bilinear map $B: V 	imes W 	o U$ factors uniquely through a linear map $	ilde{B}: V \otimes W 	o U$.

!!! theorem "Theorem 21.3 (CP Decomposition)"
    Any third-order tensor $\mathcal{X} \in \mathbb{R}^{I 	imes J 	imes K}$ can be decomposed as a sum of a finite number of rank-1 tensors:
    $$\mathcal{X} = \sum_{r=1}^R a_r \circ b_r \circ c_r$$
    The smallest such $R$ is called the **CP rank** of the tensor.

---

## Exercises

1. **[Dimension] If $\dim V = 3$ and $\dim W = 4$, find $\dim(V \otimes W)$.**
   ??? success "Solution"
       $\dim(V \otimes W) = \dim V \cdot \dim W = 3 \cdot 4 = 12$.

2. **[Tensor Rank] Give an example of a second-order tensor (matrix) with rank 1.**
   ??? success "Solution"
       An outer product $v \circ w = v w^T$. For example, $\begin{pmatrix} 1 \ 2 \end{pmatrix} \begin{pmatrix} 3 & 4 \end{pmatrix} = \begin{pmatrix} 3 & 4 \ 6 & 8 \end{pmatrix}$. Its rank is 1.

3. **[Mode-n Product] Describe the meaning of the mode-1 product $\mathcal{X} 	imes_1 M$ of a tensor $\mathcal{X}$ and matrix $M$.**
   ??? success "Solution"
       It represents applying the linear transformation $M$ along the first dimension (rows) of the tensor. If you view the tensor as a stack of matrix slices, every column in every slice is transformed by $M$.

4. **[Calculation] If $\mathcal{X} = a \circ b \circ c$, find its vectorized form $\operatorname{vec}(\mathcal{X})$.**
   ??? success "Solution"
       $\operatorname{vec}(\mathcal{X}) = c \otimes b \otimes a$ (following the standard lexicographical convention).

5. **[Tucker Decomposition] What is the role of the Core Tensor in Tucker decomposition?**
   ??? success "Solution"
       The core tensor captures the interactions between basis vectors of different dimensions. It is analogous to the singular value matrix in SVD but allows off-diagonal entries, enabling descriptions of complex cross-dimensional correlations.

6. **[Contraction] Explain the full contraction of two second-order tensors (matrices) $A, B$.**
   ??? success "Solution"
       This corresponds to the Frobenius inner product $\langle A, B angle = \sum A_{ij} B_{ij} = \operatorname{tr}(A^T B)$. Contraction operations reduce the order of tensors by summing over identical indices.

7. **[Symmetric Tensors] Provide an example of a third-order symmetric tensor.**
   ??? success "Solution"
       The third-order moment $\mathbb{E}[X \otimes X \otimes X]$ of a multivariate random variable is symmetric because $\mathbb{E}[X_i X_j X_k]$ is invariant under permutations of indices $i, j, k$.

8. **[Rank Gap] Why is the CP rank of a tensor harder to compute than the matrix rank?**
   ??? success "Solution"
       Matrix rank can be solved exactly in polynomial time via SVD. Determining tensor CP rank is NP-hard, and the rank might differ over $\mathbb{R}$ and $\mathbb{C}$. Furthermore, tensor rank can exhibit "non-closure" phenomena.

9. **[Application] Why are tensor models preferred over matrix models in hyperspectral image processing?**
   ??? success "Solution"
       Images have space (2D) and spectrum (1D) dimensions. Matrix models require forced flattening (vectorization), which destroys spatial locality. Tensor models directly preserve space-spectral structural features, leading to more effective feature extraction and denoising.

10. **[Antisymmetric] Prove that in 3D space, the dimension of the third-order fully antisymmetric tensor space is 1.**
    ??? success "Solution"
        Antisymmetric tensors correspond to the exterior algebra $\Lambda^3(\mathbb{R}^3)$. Its dimension is $\binom{3}{3} = 1$. The basis element is the Levi-Civita symbol $\epsilon_{ijk}$.

## Chapter Summary

Multilinear algebra extends the boundaries of linear thinking:

1. **Dimensional Leap**: Shifted from linear "lines" and "planes" to nested "volume" structures.
2. **Structural Decomposition**: CP and Tucker decompositions provide the algebraic blueprint for mining hidden correlations in high-dimensional data.
3. **Geometric Rigidity**: The uniqueness of tensor rank (e.g., CP decomposition is unique under weak conditions) gives tensors a massive advantage over matrices in fields like blind source separation.
