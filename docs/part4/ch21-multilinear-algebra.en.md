# Chapter 21  Multilinear Algebra and Tensors

<div class="context-flow" markdown>

**Prerequisites**: Linear maps (Ch3) · Determinants (Ch5) · Inner product spaces (Ch7) · **Arc**: Dual spaces → Tensor products (universal property) → Exterior algebra (the true home of determinants) → Tensor decompositions (high-dimensional data)
**Essence**: Linear algebra generalizes from "linearity on one vector space" to "multilinearity on multiple spaces" — dimensions multiply rather than add

</div>

Multilinear algebra is the natural generalization of linear algebra that studies multilinear relationships among multiple vector spaces. Tensors, the central objects of multilinear algebra, not only provide the fundamental language for differential geometry and general relativity, but also play an increasingly important role in machine learning, quantum computing, and data science. This chapter starts from dual spaces and systematically develops the basic theories of multilinear maps, tensor products, and exterior algebra, and introduces modern applications such as tensor decompositions.

## 21.1 Dual Spaces

<div class="context-flow" markdown>

**Motivation**: Linear functionals $V \to \mathbb{F}$ on a vector space $V$ themselves form a space $V^*$ → Dual basis and original basis "pair dually" → $V \cong V^{**}$ (canonical isomorphism, no basis needed)

</div>

### Linear Functionals and Dual Spaces

!!! definition "Definition 21.1 (Linear functional)"
    Let $V$ be a vector space over the field $\mathbb{F}$ ($\mathbb{F} = \mathbb{R}$ or $\mathbb{C}$). A **linear functional** is a linear map $f: V \to \mathbb{F}$, i.e., satisfying:

    $$f(\alpha \mathbf{u} + \beta \mathbf{v}) = \alpha f(\mathbf{u}) + \beta f(\mathbf{v}), \quad \forall \mathbf{u}, \mathbf{v} \in V, \ \alpha, \beta \in \mathbb{F}.$$

!!! definition "Definition 21.2 (Dual space)"
    Let $V$ be a vector space over $\mathbb{F}$. The **dual space** of $V$ is defined as the set of all linear functionals on $V$, denoted $V^*$:

    $$V^* = \text{Hom}(V, \mathbb{F}) = \{ f: V \to \mathbb{F} \mid f \text{ is a linear map} \}.$$

    $V^*$ forms a vector space over $\mathbb{F}$ under pointwise addition and scalar multiplication.

### Dual Basis

!!! definition "Definition 21.3 (Dual basis)"
    Let $V$ be an $n$-dimensional vector space with basis $\{\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\}$. The **dual basis** $\{\mathbf{e}^1, \mathbf{e}^2, \ldots, \mathbf{e}^n\}$ is a set of linear functionals in $V^*$ satisfying:

    $$\mathbf{e}^i(\mathbf{e}_j) = \delta^i_j = \begin{cases} 1, & i = j, \\ 0, & i \neq j. \end{cases}$$

    where $\delta^i_j$ is the Kronecker delta.

!!! theorem "Theorem 21.1 (Existence and uniqueness of dual basis)"
    Let $V$ be a finite-dimensional vector space with basis $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$. Then the dual basis $\{\mathbf{e}^1, \ldots, \mathbf{e}^n\}$ exists and is unique, and forms a basis of $V^*$. In particular, $\dim V^* = \dim V$.

??? proof "Proof"
    **Existence**: For each $i = 1, \ldots, n$, define $\mathbf{e}^i: V \to \mathbb{F}$ as follows. For any $\mathbf{v} = \sum_{j=1}^n v_j \mathbf{e}_j \in V$, let

    $$\mathbf{e}^i(\mathbf{v}) = v_i.$$

    Clearly $\mathbf{e}^i$ is linear, and $\mathbf{e}^i(\mathbf{e}_j) = \delta^i_j$.

    **Uniqueness**: If $f \in V^*$ satisfies $f(\mathbf{e}_j) = \delta^i_j$, then for any $\mathbf{v} = \sum_j v_j \mathbf{e}_j$,

    $$f(\mathbf{v}) = \sum_j v_j f(\mathbf{e}_j) = \sum_j v_j \delta^i_j = v_i = \mathbf{e}^i(\mathbf{v}).$$

    Hence $f = \mathbf{e}^i$.

    **$\{\mathbf{e}^1, \ldots, \mathbf{e}^n\}$ is a basis of $V^*$**:

    - *Linear independence*: Suppose $\sum_i c_i \mathbf{e}^i = 0$. Applying to $\mathbf{e}_j$ gives $c_j = 0$, $\forall j$.
    - *Spanning $V^*$*: Let $f \in V^*$, set $c_i = f(\mathbf{e}_i)$. Then $f$ and $\sum_i c_i \mathbf{e}^i$ agree on basis vectors, so by linearity they are equal, i.e., $f = \sum_i c_i \mathbf{e}^i$.

!!! example "Example 21.1"
    Let $V = \mathbb{R}^3$ with the standard basis $\{\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3\}$. The dual basis is $\{\mathbf{e}^1, \mathbf{e}^2, \mathbf{e}^3\}$, where

    $$\mathbf{e}^i(x_1, x_2, x_3) = x_i.$$

    That is, $\mathbf{e}^i$ is the projection function that "extracts the $i$-th coordinate."

    For example, $\mathbf{e}^2(3, -1, 5) = -1$.

### Double Dual Space and Canonical Isomorphism

<div class="context-flow" markdown>

**Key insight**: The isomorphism $V \to V^*$ depends on the choice of basis; the isomorphism $V \to V^{**}$ given by $\Phi(\mathbf{v})(f) = f(\mathbf{v})$ **does not depend on any basis** — this is what "canonical" means, foreshadowing category-theoretic ideas

</div>

!!! definition "Definition 21.4 (Double dual space)"
    The **double dual space** of $V$ is defined as $(V^*)^*$, denoted $V^{**}$. Elements of $V^{**}$ are linear functionals on $V^*$.

!!! theorem "Theorem 21.2 (Canonical isomorphism)"
    Let $V$ be a finite-dimensional vector space. Define the map $\Phi: V \to V^{**}$,

    $$\Phi(\mathbf{v})(f) = f(\mathbf{v}), \quad \forall \mathbf{v} \in V, \ f \in V^*.$$

    Then $\Phi$ is a **canonical isomorphism**, i.e., $\Phi$ is a vector space isomorphism that does not depend on the choice of basis.

??? proof "Proof"
    **$\Phi$ is a linear map**: For any $\alpha, \beta \in \mathbb{F}$, $\mathbf{u}, \mathbf{v} \in V$, $f \in V^*$,

    $$\Phi(\alpha \mathbf{u} + \beta \mathbf{v})(f) = f(\alpha \mathbf{u} + \beta \mathbf{v}) = \alpha f(\mathbf{u}) + \beta f(\mathbf{v}) = (\alpha \Phi(\mathbf{u}) + \beta \Phi(\mathbf{v}))(f).$$

    **$\Phi$ is injective**: Suppose $\Phi(\mathbf{v}) = 0$, i.e., $f(\mathbf{v}) = 0$ for all $f \in V^*$. Taking the dual basis $\mathbf{e}^i$, we get $\mathbf{e}^i(\mathbf{v}) = v_i = 0$, so $\mathbf{v} = \mathbf{0}$.

    **$\Phi$ is an isomorphism**: Since $\dim V = \dim V^{**}$ ($V$ is finite-dimensional) and $\Phi$ is injective, $\Phi$ is an isomorphism.

    **Canonicality**: The definition of $\Phi$ involves no choice of basis, hence it is "canonical." In the categorical sense, $\Phi$ constitutes a natural transformation from the identity functor to the double dual functor.

!!! note "Note"
    When $V$ is infinite-dimensional, $V$ and $V^*$ are no longer isomorphic ($V^*$ may be "larger"), and the canonical map $\Phi: V \to V^{**}$ is still injective but no longer surjective. This distinction is particularly important in functional analysis, leading to the concept of reflexive spaces.

!!! example "Example 21.2"
    Let $V = \mathbb{R}^2$ with basis $\{\mathbf{e}_1, \mathbf{e}_2\}$. The dual basis is $\{\mathbf{e}^1, \mathbf{e}^2\}$. For $\mathbf{v} = 3\mathbf{e}_1 + 2\mathbf{e}_2$,

    $$\Phi(\mathbf{v})(\mathbf{e}^1) = \mathbf{e}^1(\mathbf{v}) = 3, \quad \Phi(\mathbf{v})(\mathbf{e}^2) = \mathbf{e}^2(\mathbf{v}) = 2.$$

    Therefore $\Phi(\mathbf{v}) = 3(\mathbf{e}^1)^* + 2(\mathbf{e}^2)^*$, where $\{(\mathbf{e}^1)^*, (\mathbf{e}^2)^*\}$ is the dual basis of $\{\mathbf{e}^1, \mathbf{e}^2\}$ in $V^{**}$.

## 21.2 Multilinear Maps

<div class="context-flow" markdown>

**Transition**: Dual space = 1-linear form → Generalize to $k$-linear forms · Inner products, determinants, and matrix multiplication are all instances of multilinear maps
**Core**: A $k$-linear map is uniquely determined by $n_1 n_2 \cdots n_k$ values on basis elements → Dimensions multiply

</div>

### Bilinear Maps

!!! definition "Definition 21.5 (Bilinear map)"
    Let $V, W, U$ be vector spaces over $\mathbb{F}$. A map $B: V \times W \to U$ is called a **bilinear map** if it is linear in each variable:

    - Fixing $\mathbf{w} \in W$, the map $\mathbf{v} \mapsto B(\mathbf{v}, \mathbf{w})$ is linear;
    - Fixing $\mathbf{v} \in V$, the map $\mathbf{w} \mapsto B(\mathbf{v}, \mathbf{w})$ is linear.

    That is, for any $\alpha, \beta \in \mathbb{F}$:

    $$B(\alpha \mathbf{v}_1 + \beta \mathbf{v}_2, \mathbf{w}) = \alpha B(\mathbf{v}_1, \mathbf{w}) + \beta B(\mathbf{v}_2, \mathbf{w}),$$

    $$B(\mathbf{v}, \alpha \mathbf{w}_1 + \beta \mathbf{w}_2) = \alpha B(\mathbf{v}, \mathbf{w}_1) + \beta B(\mathbf{v}, \mathbf{w}_2).$$

!!! example "Example 21.3"
    **The inner product is a bilinear map.** Let $V$ be a real inner product space. Then the inner product $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{R}$ is a bilinear map.

    **Matrix multiplication is a bilinear map.** The map $B: \mathbb{R}^{m \times n} \times \mathbb{R}^{n \times p} \to \mathbb{R}^{m \times p}$, $B(A, B) = AB$ is bilinear.

    **The determinant as a function of row vectors.** Viewing an $n \times n$ matrix as $n$ row vectors $(\mathbf{r}_1, \ldots, \mathbf{r}_n)$, $\det(\mathbf{r}_1, \ldots, \mathbf{r}_n)$ is linear in each $\mathbf{r}_i$ (i.e., the determinant is a multilinear map).

### Multilinear Maps

!!! definition "Definition 21.6 (Multilinear map)"
    Let $V_1, V_2, \ldots, V_k, U$ be vector spaces over $\mathbb{F}$. A map $f: V_1 \times V_2 \times \cdots \times V_k \to U$ is called a **$k$-linear map** (or multilinear map) if, with the other $k-1$ variables fixed, $f$ is linear in each variable.

    The set of all $k$-linear maps from $V_1 \times \cdots \times V_k$ to $U$ is denoted $\mathcal{L}^k(V_1, \ldots, V_k; U)$. When $U = \mathbb{F}$, these are called **$k$-linear forms**.

!!! theorem "Theorem 21.3 (Multilinear maps are determined by values on bases)"
    Let $V_i$ be finite-dimensional vector spaces, $\dim V_i = n_i$, with basis $\{\mathbf{e}^{(i)}_1, \ldots, \mathbf{e}^{(i)}_{n_i}\}$ for $V_i$ ($i = 1, \ldots, k$). Then a $k$-linear map $f: V_1 \times \cdots \times V_k \to U$ is completely determined by its values on all combinations of basis vectors

    $$f(\mathbf{e}^{(1)}_{j_1}, \mathbf{e}^{(2)}_{j_2}, \ldots, \mathbf{e}^{(k)}_{j_k}), \quad 1 \leq j_i \leq n_i$$

    There are $n_1 n_2 \cdots n_k$ such values, which can be freely specified. Hence

    $$\dim \mathcal{L}^k(V_1, \ldots, V_k; U) = n_1 n_2 \cdots n_k \cdot \dim U.$$

??? proof "Proof"
    Let $\mathbf{v}_i = \sum_{j_i=1}^{n_i} v^{j_i}_i \mathbf{e}^{(i)}_{j_i} \in V_i$. By multilinearity,

    $$f(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \sum_{j_1=1}^{n_1} \cdots \sum_{j_k=1}^{n_k} v^{j_1}_1 \cdots v^{j_k}_k \, f(\mathbf{e}^{(1)}_{j_1}, \ldots, \mathbf{e}^{(k)}_{j_k}).$$

    This shows that $f$ is completely determined by the values $f(\mathbf{e}^{(1)}_{j_1}, \ldots, \mathbf{e}^{(k)}_{j_k})$.

    Conversely, for any given values $\mathbf{u}_{j_1 \cdots j_k} \in U$, one can define $f$ using the formula above, and direct verification shows that $f$ is $k$-linear.

    When $U = \mathbb{F}$, $\dim \mathcal{L}^k(V_1, \ldots, V_k; \mathbb{F}) = n_1 n_2 \cdots n_k$.

## 21.3 Tensor Products

<div class="context-flow" markdown>

**Core idea**: The tensor product linearizes "bilinear problems" — **Universal property**: Any bilinear $B: V \times W \to U$ factors uniquely as $V \otimes W \xrightarrow{\tilde{B}} U$
**Dimension**: $\dim(V \otimes W) = \dim V \cdot \dim W$ · $V^* \otimes W \cong \operatorname{Hom}(V, W)$ unifies linear maps and tensors

</div>

### Definition via Universal Property

!!! definition "Definition 21.7 (Universal property of tensor products)"
    Let $V, W$ be vector spaces over $\mathbb{F}$. The **tensor product** of $V$ and $W$ is a vector space $V \otimes W$ together with a bilinear map $\otimes: V \times W \to V \otimes W$ such that the following universal property holds:

    For any vector space $U$ and any bilinear map $B: V \times W \to U$, there exists a unique linear map $\tilde{B}: V \otimes W \to U$ such that $B = \tilde{B} \circ \otimes$, i.e., the following diagram commutes:

    $$V \times W \xrightarrow{\otimes} V \otimes W$$

    $$\searrow^{B} \quad \downarrow^{\tilde{B}}$$

    $$\qquad \quad U$$

    Elements $\mathbf{v} \otimes \mathbf{w}$ ($\mathbf{v} \in V, \mathbf{w} \in W$) are called **simple tensors** (or decomposable tensors).

<div class="context-flow" markdown>

**Insight**: The universal property guarantees that the tensor product is unique up to isomorphism — the construction (quotient space) is merely an existence proof; uniqueness comes from the abstract categorical argument

</div>

!!! theorem "Theorem 21.4 (Existence and uniqueness of tensor products)"
    For any vector spaces $V, W$, the tensor product $V \otimes W$ exists and is unique up to isomorphism.

??? proof "Proof"
    **Existence (construction)**: Let $F(V \times W)$ be the free vector space with basis consisting of all elements $(\mathbf{v}, \mathbf{w})$ of $V \times W$. Let $R$ be the subspace generated by elements of the form:

    - $(\alpha \mathbf{v}_1 + \beta \mathbf{v}_2, \mathbf{w}) - \alpha(\mathbf{v}_1, \mathbf{w}) - \beta(\mathbf{v}_2, \mathbf{w})$
    - $(\mathbf{v}, \alpha \mathbf{w}_1 + \beta \mathbf{w}_2) - \alpha(\mathbf{v}, \mathbf{w}_1) - \beta(\mathbf{v}, \mathbf{w}_2)$

    Define $V \otimes W = F(V \times W) / R$, and let $\mathbf{v} \otimes \mathbf{w}$ be the image of $(\mathbf{v}, \mathbf{w})$ in the quotient space. Then the map $\otimes: V \times W \to V \otimes W$ is bilinear and satisfies the universal property.

    **Uniqueness**: Suppose $(T_1, \otimes_1)$ and $(T_2, \otimes_2)$ both satisfy the universal property. By the universal property of $T_1$, for the bilinear map $\otimes_2$, there exists a unique linear map $\varphi: T_1 \to T_2$ such that $\otimes_2 = \varphi \circ \otimes_1$. Similarly, there exists $\psi: T_2 \to T_1$ such that $\otimes_1 = \psi \circ \otimes_2$. Then $\psi \circ \varphi \circ \otimes_1 = \otimes_1$, and by the uniqueness part of the universal property, $\psi \circ \varphi = \text{id}_{T_1}$. Similarly $\varphi \circ \psi = \text{id}_{T_2}$. Hence $T_1 \cong T_2$.

### Basis and Dimension of Tensor Products

!!! theorem "Theorem 21.5 (Basis of tensor products)"
    Let $\{\mathbf{e}_1, \ldots, \mathbf{e}_m\}$ be a basis of $V$ and $\{\mathbf{f}_1, \ldots, \mathbf{f}_n\}$ be a basis of $W$. Then

    $$\{\mathbf{e}_i \otimes \mathbf{f}_j : 1 \leq i \leq m, \ 1 \leq j \leq n\}$$

    is a basis of $V \otimes W$. In particular,

    $$\dim(V \otimes W) = \dim V \cdot \dim W.$$

??? proof "Proof"
    **Spanning**: $V \otimes W$ is generated by linear combinations of simple tensors $\mathbf{v} \otimes \mathbf{w}$. Let $\mathbf{v} = \sum_i a_i \mathbf{e}_i$, $\mathbf{w} = \sum_j b_j \mathbf{f}_j$, then

    $$\mathbf{v} \otimes \mathbf{w} = \left(\sum_i a_i \mathbf{e}_i\right) \otimes \left(\sum_j b_j \mathbf{f}_j\right) = \sum_{i,j} a_i b_j \, \mathbf{e}_i \otimes \mathbf{f}_j.$$

    **Linear independence**: Suppose $\sum_{i,j} c_{ij} \, \mathbf{e}_i \otimes \mathbf{f}_j = 0$. For each pair $(p, q)$, take linear functionals $\mathbf{e}^p \in V^*$, $\mathbf{f}^q \in W^*$, and define the bilinear map $B_{pq}(\mathbf{v}, \mathbf{w}) = \mathbf{e}^p(\mathbf{v}) \mathbf{f}^q(\mathbf{w})$. By the universal property, there exists a linear map $\tilde{B}_{pq}: V \otimes W \to \mathbb{F}$ such that $\tilde{B}_{pq}(\mathbf{e}_i \otimes \mathbf{f}_j) = \delta^p_i \delta^q_j$. Applying $\tilde{B}_{pq}$ to both sides gives $c_{pq} = 0$.

!!! example "Example 21.4"
    Let $V = \mathbb{R}^2$, $W = \mathbb{R}^3$. Then $V \otimes W$ is a $6$-dimensional space. Taking the standard bases $\{\mathbf{e}_1, \mathbf{e}_2\}$ of $V$ and $\{\mathbf{f}_1, \mathbf{f}_2, \mathbf{f}_3\}$ of $W$, the basis of $V \otimes W$ is:

    $$\{\mathbf{e}_1 \otimes \mathbf{f}_1, \mathbf{e}_1 \otimes \mathbf{f}_2, \mathbf{e}_1 \otimes \mathbf{f}_3, \mathbf{e}_2 \otimes \mathbf{f}_1, \mathbf{e}_2 \otimes \mathbf{f}_2, \mathbf{e}_2 \otimes \mathbf{f}_3\}.$$

    $V \otimes W$ is isomorphic to $\mathbb{R}^{2 \times 3}$ (the space of $2 \times 3$ matrices): $\mathbf{v} \otimes \mathbf{w} \mapsto \mathbf{v}\mathbf{w}^T$.

    For example, $\begin{pmatrix} 1 \\ 2 \end{pmatrix} \otimes \begin{pmatrix} 3 \\ 0 \\ -1 \end{pmatrix} \mapsto \begin{pmatrix} 3 & 0 & -1 \\ 6 & 0 & -2 \end{pmatrix}$.

!!! proposition "Proposition 21.1 (Basic properties of tensor products)"
    Tensor products satisfy the following properties:

    1. **Associativity**: $(V \otimes W) \otimes U \cong V \otimes (W \otimes U)$;
    2. **Commutativity**: $V \otimes W \cong W \otimes V$;
    3. **Distributivity**: $V \otimes (W \oplus U) \cong (V \otimes W) \oplus (V \otimes U)$;
    4. **Tensor product with the scalar field**: $V \otimes \mathbb{F} \cong V$;
    5. **Relation with dual spaces**: $V^* \otimes W \cong \text{Hom}(V, W)$.

??? proof "Proof"
    We prove property 5. Define the map $\Phi: V^* \otimes W \to \text{Hom}(V, W)$, defined on simple tensors as

    $$\Phi(f \otimes \mathbf{w})(\mathbf{v}) = f(\mathbf{v}) \mathbf{w}, \quad f \in V^*, \ \mathbf{w} \in W, \ \mathbf{v} \in V.$$

    By the universal property of tensor products, the map $(f, \mathbf{w}) \mapsto [\ \mathbf{v} \mapsto f(\mathbf{v})\mathbf{w}\ ]$ is bilinear, so $\Phi$ exists and is linear.

    Take a basis $\{\mathbf{e}_i\}$ of $V$, a basis $\{\mathbf{f}_j\}$ of $W$, and the dual basis $\{\mathbf{e}^i\}$. Then $\Phi(\mathbf{e}^i \otimes \mathbf{f}_j)$ is the linear map sending $\mathbf{e}_i$ to $\mathbf{f}_j$ and all other basis vectors to $\mathbf{0}$, which is exactly the standard basis of $\text{Hom}(V, W)$. Therefore $\Phi$ maps basis to basis and is an isomorphism.

## 21.4 Component Representation of Tensors

<div class="context-flow" markdown>

**Transition**: Abstract tensor products → Coordinate representation · Type $(r,s)$ tensors = $r$ contravariant + $s$ covariant indices → **Einstein convention** simplifies notation → Transformation law reflects "tensors are geometric objects independent of basis"

</div>

### Index Notation

In physics and engineering, tensors are usually represented using **index notation**. Let $V$ be an $n$-dimensional vector space with basis $\{\mathbf{e}_i\}$ and dual basis $\{\mathbf{e}^i\}$.

!!! definition "Definition 21.8 (Type $(r, s)$ tensor)"
    A **type $(r, s)$ tensor** is a multilinear map

    $$T: \underbrace{V^* \times \cdots \times V^*}_{r} \times \underbrace{V \times \cdots \times V}_{s} \to \mathbb{F}.$$

    In a given basis, the **components** of $T$ are

    $$T^{i_1 \cdots i_r}_{j_1 \cdots j_s} = T(\mathbf{e}^{i_1}, \ldots, \mathbf{e}^{i_r}, \mathbf{e}_{j_1}, \ldots, \mathbf{e}_{j_s}).$$

    Upper indices are called **contravariant indices** and lower indices are called **covariant indices**.

### Einstein Summation Convention

!!! definition "Definition 21.9 (Einstein summation convention)"
    The **Einstein summation convention** stipulates that when an index appears simultaneously as an upper index and a lower index in an expression (called a **dummy index**), summation over that index is implied. For example:

    $$a^i b_i \equiv \sum_{i=1}^n a^i b_i, \quad T^i_{\ j} v^j \equiv \sum_{j=1}^n T^i_{\ j} v^j.$$

!!! example "Example 21.5"
    Using the Einstein convention, here are some common tensor operations:

    1. **Components of a vector**: $\mathbf{v} = v^i \mathbf{e}_i$ (summed over $i$).
    2. **Action of a linear functional**: $f(\mathbf{v}) = f_i v^i$.
    3. **Matrix times vector**: $(A\mathbf{v})^i = A^i_{\ j} v^j$.
    4. **Matrix multiplication**: $(AB)^i_{\ k} = A^i_{\ j} B^j_{\ k}$.
    5. **Trace**: $\text{tr}(A) = A^i_{\ i}$.

!!! theorem "Theorem 21.6 (Transformation law for tensor components)"
    Let $\{\mathbf{e}_i\}$ and $\{\tilde{\mathbf{e}}_i\}$ be two bases of $V$ with change-of-basis matrix $P$: $\tilde{\mathbf{e}}_j = P^i_{\ j} \mathbf{e}_i$. Then the transformation law for the components of a type $(r, s)$ tensor $T$ is:

    $$\tilde{T}^{i_1 \cdots i_r}_{j_1 \cdots j_s} = (P^{-1})^{i_1}_{\ k_1} \cdots (P^{-1})^{i_r}_{\ k_r} \, P^{l_1}_{\ j_1} \cdots P^{l_s}_{\ j_s} \, T^{k_1 \cdots k_r}_{l_1 \cdots l_s}.$$

??? proof "Proof"
    The dual basis transforms as $\tilde{\mathbf{e}}^i = (P^{-1})^i_{\ j} \mathbf{e}^j$ (derived from $\tilde{\mathbf{e}}^i(\tilde{\mathbf{e}}_k) = \delta^i_k$). Therefore

    $$\tilde{T}^{i_1 \cdots i_r}_{j_1 \cdots j_s} = T(\tilde{\mathbf{e}}^{i_1}, \ldots, \tilde{\mathbf{e}}^{i_r}, \tilde{\mathbf{e}}_{j_1}, \ldots, \tilde{\mathbf{e}}_{j_s}).$$

    Substituting $\tilde{\mathbf{e}}^{i_\alpha} = (P^{-1})^{i_\alpha}_{\ k_\alpha} \mathbf{e}^{k_\alpha}$ and $\tilde{\mathbf{e}}_{j_\beta} = P^{l_\beta}_{\ j_\beta} \mathbf{e}_{l_\beta}$, the result follows by multilinearity of $T$.

## 21.5 Symmetric and Antisymmetric Tensors

<div class="context-flow" markdown>

**Branch**: $V^{\otimes k}$ decomposes under the action of the symmetric group → **Symmetric** part (Sym) + **Antisymmetric** part (Alt) · Alt is a projection operator, $\text{Alt}^2 = \text{Alt}$ → Leads to exterior algebra

</div>

!!! definition "Definition 21.10 (Symmetric and antisymmetric tensors)"
    Let $T \in V^{\otimes k} = \underbrace{V \otimes \cdots \otimes V}_{k}$ be a $k$-th order tensor.

    - $T$ is **symmetric** if for every permutation $\sigma \in S_k$, $T(\mathbf{v}_{\sigma(1)}, \ldots, \mathbf{v}_{\sigma(k)}) = T(\mathbf{v}_1, \ldots, \mathbf{v}_k)$; in components, $T_{i_{\sigma(1)} \cdots i_{\sigma(k)}} = T_{i_1 \cdots i_k}$.
    - $T$ is **antisymmetric** (alternating) if $T(\mathbf{v}_{\sigma(1)}, \ldots, \mathbf{v}_{\sigma(k)}) = \text{sgn}(\sigma) \, T(\mathbf{v}_1, \ldots, \mathbf{v}_k)$; in components, $T_{i_{\sigma(1)} \cdots i_{\sigma(k)}} = \text{sgn}(\sigma) \, T_{i_1 \cdots i_k}$.

!!! definition "Definition 21.11 (Symmetrization and antisymmetrization operators)"
    The **symmetrization operator** $\text{Sym}$ and **antisymmetrization operator** $\text{Alt}$ are defined as:

    $$\text{Sym}(T)(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \frac{1}{k!} \sum_{\sigma \in S_k} T(\mathbf{v}_{\sigma(1)}, \ldots, \mathbf{v}_{\sigma(k)}),$$

    $$\text{Alt}(T)(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \frac{1}{k!} \sum_{\sigma \in S_k} \text{sgn}(\sigma) \, T(\mathbf{v}_{\sigma(1)}, \ldots, \mathbf{v}_{\sigma(k)}).$$

!!! theorem "Theorem 21.7 (Properties of symmetrization and antisymmetrization)"
    1. $\text{Sym}$ and $\text{Alt}$ are projection operators from $V^{\otimes k}$ to itself (i.e., $\text{Sym}^2 = \text{Sym}$, $\text{Alt}^2 = \text{Alt}$).
    2. $\text{Sym}(T)$ is a symmetric tensor; $\text{Alt}(T)$ is an antisymmetric tensor.
    3. $T$ is symmetric if and only if $\text{Sym}(T) = T$; $T$ is antisymmetric if and only if $\text{Alt}(T) = T$.

??? proof "Proof"
    We prove that $\text{Alt}$ is a projection operator. Let $S = \text{Alt}(T)$, then

    $$\text{Alt}(S)(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \frac{1}{k!} \sum_{\tau \in S_k} \text{sgn}(\tau) S(\mathbf{v}_{\tau(1)}, \ldots, \mathbf{v}_{\tau(k)}).$$

    Since $S = \text{Alt}(T)$ is antisymmetric, $S(\mathbf{v}_{\tau(1)}, \ldots, \mathbf{v}_{\tau(k)}) = \text{sgn}(\tau) S(\mathbf{v}_1, \ldots, \mathbf{v}_k)$. Hence

    $$\text{Alt}(S)(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \frac{1}{k!} \sum_{\tau \in S_k} \text{sgn}(\tau)^2 S(\mathbf{v}_1, \ldots, \mathbf{v}_k) = S(\mathbf{v}_1, \ldots, \mathbf{v}_k).$$

    That is, $\text{Alt}^2 = \text{Alt}$. The remaining properties are proved similarly.

!!! example "Example 21.6"
    Let $T: \mathbb{R}^2 \times \mathbb{R}^2 \to \mathbb{R}$ be defined by $T(\mathbf{u}, \mathbf{v}) = u_1 v_2$. Then:

    $$\text{Sym}(T)(\mathbf{u}, \mathbf{v}) = \frac{1}{2}(u_1 v_2 + u_2 v_1),$$

    $$\text{Alt}(T)(\mathbf{u}, \mathbf{v}) = \frac{1}{2}(u_1 v_2 - u_2 v_1).$$

    Note that $\text{Alt}(T)(\mathbf{u}, \mathbf{v})$ is exactly $\frac{1}{2} \det \begin{pmatrix} u_1 & u_2 \\ v_1 & v_2 \end{pmatrix}$.

## 21.6 Exterior Algebra

<div class="context-flow" markdown>

**Core**: Wedge product $\wedge$ = antisymmetrized tensor product → $\dim \Lambda^k(V) = \binom{n}{k}$ → $\Lambda^n(V)$ is one-dimensional → **The determinant is the coordinate of the wedge product of $n$ vectors in $\Lambda^n$**
**Link**: The algebraic properties of determinants from Ch5 (multilinearity + antisymmetry) receive geometric interpretation here · Cross product = dual of $\Lambda^2$ in $\mathbb{R}^3$

</div>

### Wedge Product and Exterior Powers

!!! definition "Definition 21.12 (Exterior/wedge product)"
    Let $V$ be an $n$-dimensional vector space. For $\omega \in \Lambda^p(V^*)$, $\eta \in \Lambda^q(V^*)$ (where $\Lambda^k(V^*)$ is the space of $k$-th order antisymmetric tensors on $V^*$), the **wedge product** $\omega \wedge \eta \in \Lambda^{p+q}(V^*)$ is defined as:

    $$\omega \wedge \eta = \frac{(p+q)!}{p! \, q!} \, \text{Alt}(\omega \otimes \eta).$$

    For the vector space $V$ itself, the **exterior power** $\Lambda^k(V)$ is defined as the subspace of antisymmetric tensors in $V^{\otimes k}$. For $\mathbf{v}_1, \ldots, \mathbf{v}_k \in V$, define

    $$\mathbf{v}_1 \wedge \mathbf{v}_2 \wedge \cdots \wedge \mathbf{v}_k = \sum_{\sigma \in S_k} \text{sgn}(\sigma) \, \mathbf{v}_{\sigma(1)} \otimes \cdots \otimes \mathbf{v}_{\sigma(k)}.$$

!!! theorem "Theorem 21.8 (Properties of the wedge product)"
    The wedge product satisfies the following properties:

    1. **Bilinearity**: $(\alpha \omega_1 + \beta \omega_2) \wedge \eta = \alpha (\omega_1 \wedge \eta) + \beta (\omega_2 \wedge \eta)$;
    2. **Associativity**: $(\omega \wedge \eta) \wedge \zeta = \omega \wedge (\eta \wedge \zeta)$;
    3. **Graded commutativity**: $\omega \wedge \eta = (-1)^{pq} \eta \wedge \omega$, where $\omega \in \Lambda^p$, $\eta \in \Lambda^q$;
    4. **Antisymmetry**: $\mathbf{v} \wedge \mathbf{v} = 0$, $\forall \mathbf{v} \in V$.

??? proof "Proof"
    We prove property 3. Let $\sigma_0$ be the permutation sending $(1, \ldots, p, p+1, \ldots, p+q)$ to $(p+1, \ldots, p+q, 1, \ldots, p)$. This permutation requires $pq$ transpositions (moving $q$ elements past $p$ elements one by one), so $\text{sgn}(\sigma_0) = (-1)^{pq}$.

    By properties of antisymmetric tensors:

    $$(\omega \wedge \eta)(\mathbf{v}_{\sigma_0(1)}, \ldots, \mathbf{v}_{\sigma_0(p+q)}) = (-1)^{pq} (\omega \wedge \eta)(\mathbf{v}_1, \ldots, \mathbf{v}_{p+q}).$$

    And $(\mathbf{v}_{\sigma_0(1)}, \ldots, \mathbf{v}_{\sigma_0(p+q)}) = (\mathbf{v}_{p+1}, \ldots, \mathbf{v}_{p+q}, \mathbf{v}_1, \ldots, \mathbf{v}_p)$; by the definition of the wedge product, this gives exactly $(\eta \wedge \omega)(\mathbf{v}_1, \ldots, \mathbf{v}_{p+q})$.

!!! theorem "Theorem 21.9 (Dimension of exterior powers)"
    Let $\dim V = n$. Then

    $$\dim \Lambda^k(V) = \binom{n}{k}.$$

    In particular, $\Lambda^0(V) = \mathbb{F}$, $\Lambda^1(V) = V$, $\Lambda^n(V)$ is $1$-dimensional, and $\Lambda^k(V) = \{0\}$ for $k > n$.

??? proof "Proof"
    Let $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ be a basis of $V$. Then

    $$\{\mathbf{e}_{i_1} \wedge \mathbf{e}_{i_2} \wedge \cdots \wedge \mathbf{e}_{i_k} : 1 \leq i_1 < i_2 < \cdots < i_k \leq n\}$$

    forms a basis of $\Lambda^k(V)$ (linear independence and spanning can be verified directly). The number of such basis vectors is $\binom{n}{k}$.

    When $k > n$, in $\mathbf{e}_{i_1} \wedge \cdots \wedge \mathbf{e}_{i_k}$ there must be repeated basis vectors, and by antisymmetry $\mathbf{e}_{i_1} \wedge \cdots \wedge \mathbf{e}_{i_k} = 0$.

!!! definition "Definition 21.13 (Exterior algebra)"
    The **exterior algebra** on $V$ is defined as

    $$\Lambda(V) = \bigoplus_{k=0}^{n} \Lambda^k(V),$$

    equipped with the wedge product operation. $\Lambda(V)$ is a graded anti-commutative algebra with total dimension $\sum_{k=0}^n \binom{n}{k} = 2^n$.

### Determinant as Exterior Product

<div class="context-flow" markdown>

**Insight**: $\mathbf{v}_1 \wedge \cdots \wedge \mathbf{v}_n = \det(A) \, \mathbf{e}_1 \wedge \cdots \wedge \mathbf{e}_n$ — the Leibniz formula for determinants is no longer "defined out of nowhere," but is the natural consequence of multilinearity + antisymmetry of the wedge product

</div>

!!! theorem "Theorem 21.10 (Exterior product interpretation of determinants)"
    Let $\mathbf{v}_1, \ldots, \mathbf{v}_n \in V$ ($\dim V = n$), with $\mathbf{v}_i = \sum_j a_{ji} \mathbf{e}_j$ (i.e., $a_{ji}$ is the $j$-th component of $\mathbf{v}_i$). Then

    $$\mathbf{v}_1 \wedge \mathbf{v}_2 \wedge \cdots \wedge \mathbf{v}_n = \det(A) \, \mathbf{e}_1 \wedge \mathbf{e}_2 \wedge \cdots \wedge \mathbf{e}_n,$$

    where $A = (a_{ij})$ is the matrix with $\mathbf{v}_j$ as columns.

??? proof "Proof"
    By multilinearity and antisymmetry of the wedge product:

    $$\mathbf{v}_1 \wedge \cdots \wedge \mathbf{v}_n = \left(\sum_{i_1} a_{i_1 1} \mathbf{e}_{i_1}\right) \wedge \cdots \wedge \left(\sum_{i_n} a_{i_n n} \mathbf{e}_{i_n}\right)$$

    $$= \sum_{i_1, \ldots, i_n} a_{i_1 1} \cdots a_{i_n n} \, \mathbf{e}_{i_1} \wedge \cdots \wedge \mathbf{e}_{i_n}.$$

    Since $\mathbf{e}_{i_1} \wedge \cdots \wedge \mathbf{e}_{i_n} = 0$ whenever $(i_1, \ldots, i_n)$ contains repeated indices, the nonzero terms occur only when $(i_1, \ldots, i_n)$ is a permutation $\sigma$ of $(1, \ldots, n)$. In that case $\mathbf{e}_{\sigma(1)} \wedge \cdots \wedge \mathbf{e}_{\sigma(n)} = \text{sgn}(\sigma) \, \mathbf{e}_1 \wedge \cdots \wedge \mathbf{e}_n$. Hence

    $$\mathbf{v}_1 \wedge \cdots \wedge \mathbf{v}_n = \left(\sum_{\sigma \in S_n} \text{sgn}(\sigma) \, a_{\sigma(1),1} \cdots a_{\sigma(n),n}\right) \mathbf{e}_1 \wedge \cdots \wedge \mathbf{e}_n = \det(A) \, \mathbf{e}_1 \wedge \cdots \wedge \mathbf{e}_n.$$

!!! example "Example 21.7"
    Let $V = \mathbb{R}^3$, $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 2 \\ 0 \end{pmatrix}$, $\mathbf{v}_2 = \begin{pmatrix} 0 \\ 1 \\ 3 \end{pmatrix}$. Then

    $$\mathbf{v}_1 \wedge \mathbf{v}_2 = (1 \cdot 1 - 2 \cdot 0)(\mathbf{e}_1 \wedge \mathbf{e}_2) + (1 \cdot 3 - 0 \cdot 0)(\mathbf{e}_1 \wedge \mathbf{e}_3) + (2 \cdot 3 - 0 \cdot 1)(\mathbf{e}_2 \wedge \mathbf{e}_3)$$

    $$= \mathbf{e}_1 \wedge \mathbf{e}_2 + 3 \, \mathbf{e}_1 \wedge \mathbf{e}_3 + 6 \, \mathbf{e}_2 \wedge \mathbf{e}_3.$$

    This corresponds to $\mathbf{v}_1 \times \mathbf{v}_2 = (6, -3, 1)^T$ (the connection between the exterior product and the cross product).

## 21.7 Tensor Decompositions

<div class="context-flow" markdown>

**Transition**: Theory → Applications · Matrix SVD (Ch8) generalized to higher orders → **CP decomposition** (sum of rank-ones) / **Tucker decomposition** (core tensor + factor matrices)
**Warning**: Tensor rank $\neq$ matrix rank — computing it is NP-hard, real rank $\neq$ complex rank · Links to Ch25 low-rank approximation / recommender systems

</div>

Tensor decomposition is the technique of representing higher-order tensors as sums of simple tensors, with broad applications in signal processing, machine learning, and chemometrics.

### CP Decomposition

!!! definition "Definition 21.14 (CP decomposition)"
    The **CP decomposition** (Canonical Polyadic decomposition, also called CANDECOMP/PARAFAC decomposition) represents a $k$-th order tensor $\mathcal{T} \in \mathbb{R}^{n_1 \times n_2 \times \cdots \times n_k}$ as a sum of rank-one tensors:

    $$\mathcal{T} = \sum_{r=1}^{R} \lambda_r \, \mathbf{a}^{(1)}_r \otimes \mathbf{a}^{(2)}_r \otimes \cdots \otimes \mathbf{a}^{(k)}_r,$$

    where $\lambda_r \in \mathbb{R}$, $\mathbf{a}^{(i)}_r \in \mathbb{R}^{n_i}$. The smallest positive integer $R$ for which this holds is called the **tensor rank** of $\mathcal{T}$.

!!! theorem "Theorem 21.11 (Differences between tensor rank and matrix rank)"
    Unlike matrix rank, tensor rank has the following properties:

    1. Tensor rank can exceed the matrix rank of any mode unfolding;
    2. Over the real and complex fields, the rank of the same tensor may differ;
    3. Computing tensor rank is NP-hard.

??? proof "Proof"
    A classic counterexample for property 1: consider the $2 \times 2 \times 2$ tensor $\mathcal{T}$ with components $\mathcal{T}_{111} = \mathcal{T}_{222} = 1$ and all other components zero. One can verify that all three mode unfoldings have matrix rank $2$, and $\mathcal{T}$ has tensor rank $2$ (over the reals).

    An example for property 2: there exist tensors with real rank $3$ but complex rank $2$. For instance, certain $2 \times 2 \times 2$ tensors have real rank $3$ but can be decomposed as a sum of $2$ rank-one tensors over $\mathbb{C}$.

    Property 3 was rigorously proved by Hillar and Lim (2013).

### Tucker Decomposition

!!! definition "Definition 21.15 (Tucker decomposition)"
    The **Tucker decomposition** represents a $k$-th order tensor as the product of a core tensor and mode factor matrices:

    $$\mathcal{T} = \mathcal{G} \times_1 U^{(1)} \times_2 U^{(2)} \times_3 \cdots \times_k U^{(k)},$$

    where $\mathcal{G} \in \mathbb{R}^{r_1 \times r_2 \times \cdots \times r_k}$ is the core tensor, $U^{(i)} \in \mathbb{R}^{n_i \times r_i}$ are factor matrices, and $\times_i$ denotes the mode-$i$ matrix product.

    In component form:

    $$\mathcal{T}_{i_1 i_2 \cdots i_k} = \sum_{j_1=1}^{r_1} \cdots \sum_{j_k=1}^{r_k} \mathcal{G}_{j_1 j_2 \cdots j_k} \, U^{(1)}_{i_1 j_1} U^{(2)}_{i_2 j_2} \cdots U^{(k)}_{i_k j_k}.$$

!!! example "Example 21.8"
    Consider a $3 \times 4 \times 2$ third-order tensor $\mathcal{T}$ with CP decomposition

    $$\mathcal{T} = \mathbf{a}_1 \otimes \mathbf{b}_1 \otimes \mathbf{c}_1 + \mathbf{a}_2 \otimes \mathbf{b}_2 \otimes \mathbf{c}_2,$$

    where $\mathbf{a}_r \in \mathbb{R}^3$, $\mathbf{b}_r \in \mathbb{R}^4$, $\mathbf{c}_r \in \mathbb{R}^2$. In component form,

    $$\mathcal{T}_{ijk} = a_{1i} b_{1j} c_{1k} + a_{2i} b_{2j} c_{2k}, \quad 1 \leq i \leq 3, \ 1 \leq j \leq 4, \ 1 \leq k \leq 2.$$

    This is a tensor of rank $2$. A Tucker decomposition can use $U^{(1)} \in \mathbb{R}^{3 \times 2}$, $U^{(2)} \in \mathbb{R}^{4 \times 2}$, $U^{(3)} \in \mathbb{R}^{2 \times 2}$, and core tensor $\mathcal{G} \in \mathbb{R}^{2 \times 2 \times 2}$.

## 21.8 Applications of Tensors

<div class="context-flow" markdown>

**Applications**: Recommender systems (user $\otimes$ item $\otimes$ time) · Quantum computing ($(\mathbb{C}^2)^{\otimes n}$, entanglement = tensor rank > 1) → Tensor networks and MPS are the computational language of quantum many-body physics

</div>

### Multidimensional Data Analysis

Tensor methods are increasingly widespread in data science. Multidimensional data (such as user-item-time recommendation data, color image sequences, EEG signals, etc.) are naturally represented by tensors.

!!! example "Example 21.9"
    **Tensor decomposition in recommender systems**: Suppose there are $m$ users, $n$ items, and $p$ time points. The rating of user $i$ for item $j$ at time $k$ forms a third-order tensor $\mathcal{T} \in \mathbb{R}^{m \times n \times p}$. The CP decomposition represents it as

    $$\mathcal{T} \approx \sum_{r=1}^{R} \mathbf{u}_r \otimes \mathbf{v}_r \otimes \mathbf{w}_r,$$

    where $\mathbf{u}_r \in \mathbb{R}^m$ encodes user features, $\mathbf{v}_r \in \mathbb{R}^n$ encodes item features, and $\mathbf{w}_r \in \mathbb{R}^p$ encodes temporal patterns. Low-rank approximation naturally performs data completion and dimensionality reduction.

### Tensors in Quantum Information

!!! example "Example 21.10"
    **Quantum entangled states.** In quantum computing, the state space of $n$ qubits is $(\mathbb{C}^2)^{\otimes n}$, i.e., the tensor product of $n$ copies of $\mathbb{C}^2$, with dimension $2^n$. A general quantum state $|\psi\rangle \in (\mathbb{C}^2)^{\otimes n}$ can be written as

    $$|\psi\rangle = \sum_{i_1, \ldots, i_n \in \{0,1\}} c_{i_1 \cdots i_n} |i_1\rangle \otimes \cdots \otimes |i_n\rangle.$$

    When $|\psi\rangle$ cannot be written in the form $|\phi_1\rangle \otimes \cdots \otimes |\phi_n\rangle$, the state is called **entangled**. For example, the Bell state

    $$|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$$

    is an entangled state. Determining entanglement is closely related to tensor rank.

!!! proposition "Proposition 21.2 (Entanglement and tensor rank)"
    A pure state $|\psi\rangle \in \mathcal{H}_A \otimes \mathcal{H}_B$ is separable if and only if its tensor rank is $1$, i.e., it can be written as $|\psi\rangle = |\phi_A\rangle \otimes |\phi_B\rangle$. Entangled states correspond to tensor rank greater than $1$.

??? proof "Proof"
    This is a direct consequence of the definition of tensor products. $|\psi\rangle$ being separable means $|\psi\rangle = |\phi_A\rangle \otimes |\phi_B\rangle$, which is exactly a simple (rank-one) tensor. Conversely, if $|\psi\rangle$ has tensor rank $1$, then by definition it can be written as a simple tensor, i.e., a separable state.
