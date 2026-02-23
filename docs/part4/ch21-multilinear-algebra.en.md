# Chapter 21: Multilinear Algebra and Tensors

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Linear Transformations (Ch05) · Inner Product Spaces (Ch08)

**Chapter Outline**: Dual Spaces & Linear Functionals → Multilinear Maps → Tensor Product of Vector Spaces → Universality → Tensor Product of Linear Maps → Tensors of Type $(r, s)$ → Einstein Summation Convention → Symmetric and Antisymmetric Tensors → Exterior Algebra (Wedge Product) → Determinants as Exterior Products → Tensor Rank → CP and Tucker Decompositions → Applications in Physics and Data Science

**Extension**: Multilinear algebra is the mathematical language of General Relativity (Riemann Curvature), Quantum Information (Entanglement), and modern High-order Data Analysis; Exterior algebra provides the foundation for Differential Forms and modern Geometry.

</div>

Linear algebra primarily deals with mappings between two spaces (matrices). **Multilinear Algebra** extends these concepts to relationships involving multiple spaces simultaneously. The central object of this field is the **Tensor**, a universal construction that generalizes vectors and matrices. Tensors allow for the description of complex interactions—such as the coupling of particles in quantum mechanics or multi-way correlations in big data—that cannot be captured by simple linear arrays. This chapter systematically constructs the theory of tensors from dual spaces to high-order decompositions.

---

## 21.1 Dual Spaces and Linear Functionals

<div class="context-flow" markdown>

**Motivation**: A vector space $V$ is a set of "arrows." Its dual $V^*$ is the set of "measurement devices" that act on these arrows to produce scalars.

</div>

!!! definition "Definition 21.1 (Linear Functional)"
    Let $V$ be a vector space over $F$. A **linear functional** is a linear mapping $f: V \to F$.
    The set of all such functionals forms the **dual space** $V^* = \operatorname{Hom}(V, F)$.

!!! theorem "Theorem 21.1 (Dual Basis)"
    If $V$ has a basis $\{e_1, \ldots, e_n\}$, then $V^*$ has a unique **dual basis** $\{e^1, \ldots, e^n\}$ defined by the property:
    $$e^i(e_j) = \delta^i_j = \begin{cases} 1 & i=j \\ 0 & i \neq j \end{cases}$$
    Consequently, $\dim V^* = \dim V$ for finite-dimensional spaces.

!!! theorem "Theorem 21.2 (Natural Isomorphism to Double Dual)"
    There is a **canonical isomorphism** $\Phi: V \to V^{**}$ defined by $\Phi(v)(f) = f(v)$. This map is independent of the choice of basis, meaning we can view vectors as operators acting on measurements.

---

## 21.2 Multilinear Mappings

!!! definition "Definition 21.2 (Multilinear Map)"
    A map $f: V_1 \times V_2 \times \cdots \times V_k \to W$ is **$k$-linear** (or multilinear) if it is linear in each argument separately.
    - If $k=2$, it is a **bilinear map**.
    - Examples include the dot product, cross product (trilinear in $\mathbb{R}^3$), and the determinant (multilinear in the rows of a matrix).

---

## 21.3 The Tensor Product

<div class="context-flow" markdown>

**The Core Idea**: The tensor product $V \otimes W$ is the "largest" space that linearizes all bilinear maps from $V \times W$. It is the universal solution to the problem of turning products of vectors into linear combinations.

</div>

!!! definition "Definition 21.3 (Universal Property of Tensor Products)"
    The **tensor product** of vector spaces $V$ and $W$ is a vector space $V \otimes W$ together with a bilinear map $\otimes: V \times W \to V \otimes W$ such that for any bilinear map $B: V \times W \to U$, there exists a unique linear map $\tilde{B}: V \otimes W \to U$ making the diagram commute:
    $$B(v, w) = \tilde{B}(v \otimes w)$$

!!! theorem "Theorem 21.3 (Basis and Dimension)"
    If $\{e_i\}$ is a basis for $V$ and $\{f_j\}$ is a basis for $W$, then $\{e_i \otimes f_j\}$ is a basis for $V \otimes W$.
    Thus, $\dim(V \otimes W) = \dim V \cdot \dim W$.

---

## 21.4 Tensors of Type $(r, s)$

!!! definition "Definition 21.4 (General Tensors)"
    A tensor of type $(r, s)$ is a multilinear functional:
    $$T: \underbrace{V^* \times \cdots \times V^*}_{r} \times \underbrace{V \times \cdots \times V}_{s} \to F$$
    - Type (1, 0) tensors are **vectors** in $V$.
    - Type (0, 1) tensors are **covectors** in $V^*$.
    - Type (1, 1) tensors are **linear operators** in $\operatorname{End}(V) \cong V \otimes V^*$.

!!! notation "Einstein Summation Convention"
    In tensor calculus, summation is implied over repeated upper and lower indices:
    $$T^i_j v^j \equiv \sum_{j=1}^n T^i_j v^j$$
    Upper indices denote **contravariance** (vectors), and lower indices denote **covariance** (functionals).

---

## 21.5 Exterior Algebra and Wedge Products

!!! definition "Definition 21.5 (Wedge Product)"
    The **wedge product** $\mathbf{u} \wedge \mathbf{v}$ is an antisymmetric tensor product:
    $$\mathbf{u} \wedge \mathbf{v} = \mathbf{u} \otimes \mathbf{v} - \mathbf{v} \otimes \mathbf{u}$$
    It represents the oriented area element spanned by $\mathbf{u}$ and $\mathbf{v}$.

!!! theorem "Theorem 21.4 (Determinants and Exterior Products)"
    The determinant of an $n \times n$ matrix $A$ with columns $\mathbf{a}_1, \ldots, \mathbf{a}_n$ is the unique scalar such that:
    $$\mathbf{a}_1 \wedge \mathbf{a}_2 \wedge \cdots \wedge \mathbf{a}_n = \det(A) (\mathbf{e}_1 \wedge \mathbf{e}_2 \wedge \cdots \wedge \mathbf{e}_n)$$
    This provides the geometric definition of the determinant as a volume scaling factor.

---

## 21.6 Tensor Decompositions

<div class="context-flow" markdown>

**Challenge**: Unlike matrices (rank-2 tensors), high-order tensors do not have a unique "best" rank decomposition. Tensor rank calculation is NP-hard.

</div>

!!! definition "Definition 21.6 (CP Decomposition)"
    The **Canonical Polyadic (CP)** decomposition expresses a tensor as a sum of a minimum number of rank-1 tensors:
    $$\mathcal{T} = \sum_{r=1}^R \mathbf{a}_r \otimes \mathbf{b}_r \otimes \mathbf{c}_r$$
    The smallest $R$ is the **tensor rank**.

!!! definition "Definition 21.7 (Tucker Decomposition)"
    The **Tucker decomposition** generalizes SVD to tensors by extracting a "core tensor" and factor matrices for each dimension:
    $$\mathcal{T} \approx \mathcal{G} \times_1 U \times_2 V \times_3 W$$

---

## Exercises

1.  **[Dual] Find the dual basis for $\{(1, 1), (0, 1)\}$ in $\mathbb{R}^2$.**
    ??? success "Solution"
        Let $\{f^1, f^2\}$ be the dual basis. $f^1(1, 1)=1, f^1(0, 1)=0 \implies f^1(x, y) = x$.
        $f^2(1, 1)=0, f^2(0, 1)=1 \implies f^2(x, y) = y - x$.

2.  **[Dimension] If $\dim V = 3$ and $\dim W = 4$, what is the dimension of the space of bilinear forms on $V \times W$?**
    ??? success "Solution"
        The space of bilinear forms is isomorphic to $V^* \otimes W^*$. The dimension is $3 \times 4 = 12$.

3.  **[Simple Tensors] Is every element in $V \otimes W$ a "simple tensor" $v \otimes w$?**
    ??? success "Solution"
        No. Most elements are linear combinations of simple tensors. In quantum mechanics, non-simple tensors represent **entangled states**.

4.  **[Operators] Explain the isomorphism $\operatorname{End}(V) \cong V \otimes V^*$.**
    ??? success "Solution"
        A simple tensor $v \otimes f$ acts as an operator: $(v \otimes f)(u) = f(u)v$. This is a rank-1 operator. Sums of these represent all linear operators.

5.  **[Einstein] Expand the expression $A^i_j B^j_k$ using standard summation.**
    ??? success "Solution"
        $\sum_{j=1}^n A_{ij} B_{jk}$. This is the standard definition of matrix multiplication $(AB)_{ik}$.

6.  **[Wedge] Prove that $v \wedge v = 0$ for any vector $v$.**
    ??? success "Solution"
        By antisymmetry, $v \wedge v = - (v \wedge v)$. Thus $2(v \wedge v) = 0 \implies v \wedge v = 0$. Geometrically, a vector spans zero area with itself.

7.  **[Exterior Dim] What is the dimension of $\Lambda^k(V)$ if $\dim V = n$?**
    ??? success "Solution"
        The dimension is the binomial coefficient $\binom{n}{k}$.

8.  **[Rank] How does tensor rank differ from matrix rank?**
    ??? success "Solution"
        Matrix rank is easy to compute (SVD/Gaussian). Tensor rank is NP-hard to compute, can change over different fields (Real vs. Complex), and can exceed the dimensions of the tensor.

9.  **[Contraction] Define the contraction of a $(1, 1)$ tensor $T^i_j$.**
    ??? success "Solution"
        The contraction is the scalar $T^i_i$ (summing over $i=j$), which is exactly the **trace** of the corresponding matrix.

10. **[Physics] Why is the Stress Tensor in mechanics a type (0, 2) tensor?**

   ??? success "Solution"
        Stress measures force (a covector/functional) acting on an area element (another vector/functional property). It maps two vectors to a scalar energy density.

## Chapter Summary

This chapter elevates linear mapping to the universal language of tensors:

1.  **Measurement Duality**: Dual spaces formalize the distinction between vectors ("arrows") and functionals ("measurements"), enabling rigorous index calculus.
2.  **Linearization of Products**: The tensor product provides the mathematical machinery to treat non-linear interactions as linear objects in a higher-dimensional space.
3.  **Geometric Exterior**: Exterior algebra captures the concept of oriented volume and antisymmetry, providing the coordinate-free foundation for determinants and differential forms.
4.  **High-order Complexity**: Tensor decompositions extend SVD to multi-way data, revealing the profound computational challenges inherent in high-dimensional structures.
