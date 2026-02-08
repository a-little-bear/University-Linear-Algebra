# Chapter 5  Linear Transformations

<div class="context-flow" markdown>

**Prerequisites**: Chapter 4 vector spaces · basis and dimension · coordinates

**Chapter arc**: definition of linear transformation → kernel and image → rank-nullity theorem → matrix representation → change of basis/similarity → composition and inverse → isomorphism → invariant subspaces

**Further connections**：Linear transformations are fundamental in differential equations (linear operator $D = d/dx$), quantum mechanics (operator algebras), and signal processing (linear filters); infinite-dimensional operator theory (spectral theory, compact operators) is central to functional analysis

</div>

A linear transformation is one of the most central concepts in linear algebra, unifying the algebraic structure and geometric intuition of vector spaces. In previous chapters, we have seen that matrix multiplication $A\mathbf{x}$ maps one vector to another; this chapter takes the abstract viewpoint and recognizes that the essence of such a mapping — preserving addition and scalar multiplication — is the truly critical property. We will systematically study the definition and properties of linear transformations, kernel and image, the rank-nullity theorem, matrix representation, change of basis and similar matrices, composition and inverse, isomorphism, and invariant subspaces, laying a solid foundation for the subsequent theories of eigenvalues and inner product spaces.

---

## 5.1 Definition of Linear Transformations

<div class="context-flow" markdown>

**Core axiom**: $T(c\mathbf{u}+d\mathbf{v}) = cT(\mathbf{u})+dT(\mathbf{v})$ → preserves linear combinations → rotation, projection, differentiation are all linear transformations; translation is not → $T$ is **uniquely determined** by its values on a basis (Theorem 5.2)

</div>

!!! definition "Definition 5.1 (Linear transformation)"
    Let $V$ and $W$ be vector spaces over a field $\mathbb{F}$. A mapping $T: V \to W$ is called a **linear transformation** from $V$ to $W$ if for all $\mathbf{u}, \mathbf{v} \in V$ and all scalars $c \in \mathbb{F}$:

    1. **Additivity**: $T(\mathbf{u} + \mathbf{v}) = T(\mathbf{u}) + T(\mathbf{v})$
    2. **Homogeneity**: $T(c\mathbf{v}) = cT(\mathbf{v})$

    When $W = V$, $T$ is also called a **linear operator** on $V$.

!!! theorem "Theorem 5.1 (Equivalent conditions for linear transformations)"
    Let $T: V \to W$ be a mapping. Then the following conditions are equivalent:

    1. $T$ is a linear transformation.
    2. For all $\mathbf{u}, \mathbf{v} \in V$ and $c, d \in \mathbb{F}$, $T(c\mathbf{u} + d\mathbf{v}) = cT(\mathbf{u}) + dT(\mathbf{v})$.
    3. For any finite collection of vectors $\mathbf{v}_1, \ldots, \mathbf{v}_k \in V$ and scalars $c_1, \ldots, c_k \in \mathbb{F}$,

    $$T\left(\sum_{i=1}^{k} c_i \mathbf{v}_i\right) = \sum_{i=1}^{k} c_i T(\mathbf{v}_i)$$

??? proof "Proof"
    $(1) \Rightarrow (2)$: By additivity, $T(c\mathbf{u} + d\mathbf{v}) = T(c\mathbf{u}) + T(d\mathbf{v})$; then by homogeneity, $= cT(\mathbf{u}) + dT(\mathbf{v})$.

    $(2) \Rightarrow (3)$: By induction on $k$. When $k = 1$, setting $d = 0$ gives $T(c_1 \mathbf{v}_1) = c_1 T(\mathbf{v}_1)$. Assuming the result holds for $k-1$,

    $$T\left(\sum_{i=1}^{k} c_i \mathbf{v}_i\right) = T\left(\sum_{i=1}^{k-1} c_i \mathbf{v}_i + c_k \mathbf{v}_k\right) = T\left(\sum_{i=1}^{k-1} c_i \mathbf{v}_i\right) + c_k T(\mathbf{v}_k) = \sum_{i=1}^{k} c_i T(\mathbf{v}_i)$$

    $(3) \Rightarrow (1)$: Setting $k = 2$, $c_1 = c_2 = 1$ yields additivity; setting $k = 1$ yields homogeneity. $\blacksquare$

!!! proposition "Proposition 5.1 (Basic properties of linear transformations)"
    Let $T: V \to W$ be a linear transformation. Then:

    1. $T(\mathbf{0}_V) = \mathbf{0}_W$
    2. $T(-\mathbf{v}) = -T(\mathbf{v})$
    3. $T(\mathbf{u} - \mathbf{v}) = T(\mathbf{u}) - T(\mathbf{v})$

??? proof "Proof"
    1. $T(\mathbf{0}) = T(0 \cdot \mathbf{v}) = 0 \cdot T(\mathbf{v}) = \mathbf{0}$.
    2. $T(-\mathbf{v}) = T((-1)\mathbf{v}) = (-1)T(\mathbf{v}) = -T(\mathbf{v})$.
    3. $T(\mathbf{u} - \mathbf{v}) = T(\mathbf{u} + (-\mathbf{v})) = T(\mathbf{u}) + T(-\mathbf{v}) = T(\mathbf{u}) - T(\mathbf{v})$. $\blacksquare$

!!! example "Example 5.1"
    **Rotation transformation**: $T: \mathbb{R}^2 \to \mathbb{R}^2$ rotates vectors in the plane counterclockwise by angle $\theta$:

    $$T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} x\cos\theta - y\sin\theta \\ x\sin\theta + y\cos\theta \end{pmatrix}$$

    One can verify that $T$ satisfies additivity and homogeneity, so it is a linear transformation. Geometrically, it rigidly rotates the entire plane by angle $\theta$.

!!! example "Example 5.2"
    **Orthogonal projection**: $T: \mathbb{R}^3 \to \mathbb{R}^3$ projects vectors in space onto the $xy$-plane:

    $$T\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} x \\ y \\ 0 \end{pmatrix}$$

    Verification: $T(c\mathbf{u} + d\mathbf{v}) = cT(\mathbf{u}) + dT(\mathbf{v})$, so $T$ is a linear transformation.

!!! example "Example 5.3"
    **Differentiation operator**: Let $V = C^1(\mathbb{R})$ be the space of all continuously differentiable functions, and $W = C(\mathbb{R})$ the space of all continuous functions. Define $D: V \to W$ by $D(f) = f'$. By linearity of the derivative:

    $$D(cf + dg) = (cf + dg)' = cf' + dg' = cD(f) + dD(g)$$

    Therefore the differentiation operator $D$ is a linear transformation.

!!! example "Example 5.4"
    **Counterexample of a nonlinear map**: The mapping $T: \mathbb{R}^2 \to \mathbb{R}^2$, $T(\mathbf{v}) = \mathbf{v} + \begin{pmatrix} 1 \\ 0 \end{pmatrix}$ (translation) is not a linear transformation, since $T(\mathbf{0}) = \begin{pmatrix} 1 \\ 0 \end{pmatrix} \neq \mathbf{0}$.

!!! theorem "Theorem 5.2 (A linear transformation is uniquely determined by its values on a basis)"
    Let $V$ be a finite-dimensional vector space with basis $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$, and let $W$ be any vector space. For any given $\mathbf{w}_1, \ldots, \mathbf{w}_n \in W$, there exists a unique linear transformation $T: V \to W$ such that

    $$T(\mathbf{v}_i) = \mathbf{w}_i, \quad i = 1, 2, \ldots, n$$

??? proof "Proof"
    **Existence**: Every $\mathbf{v} \in V$ can be uniquely written as $\mathbf{v} = c_1 \mathbf{v}_1 + \cdots + c_n \mathbf{v}_n$. Define

    $$T(\mathbf{v}) = c_1 \mathbf{w}_1 + \cdots + c_n \mathbf{w}_n$$

    Verification of linearity: let $\mathbf{u} = a_1 \mathbf{v}_1 + \cdots + a_n \mathbf{v}_n$, $\mathbf{v} = b_1 \mathbf{v}_1 + \cdots + b_n \mathbf{v}_n$. Then

    $$T(\alpha\mathbf{u} + \beta\mathbf{v}) = T\left(\sum_i (\alpha a_i + \beta b_i)\mathbf{v}_i\right) = \sum_i (\alpha a_i + \beta b_i)\mathbf{w}_i = \alpha \sum_i a_i \mathbf{w}_i + \beta \sum_i b_i \mathbf{w}_i = \alpha T(\mathbf{u}) + \beta T(\mathbf{v})$$

    **Uniqueness**: If $S: V \to W$ also satisfies $S(\mathbf{v}_i) = \mathbf{w}_i$, then for any $\mathbf{v} = \sum c_i \mathbf{v}_i$,

    $$S(\mathbf{v}) = \sum c_i S(\mathbf{v}_i) = \sum c_i \mathbf{w}_i = T(\mathbf{v})$$

    Hence $S = T$. $\blacksquare$

---

## 5.2 Kernel and Image

<div class="context-flow" markdown>

**Information loss and coverage**: $\ker(T)$ = vectors mapped to $\mathbf{0}$ (generalization of Chapter 1 homogeneous solution sets) → $\operatorname{im}(T)$ = range → injective $\Leftrightarrow$ $\ker(T)=\{\mathbf{0}\}$

</div>

!!! definition "Definition 5.2 (Kernel)"
    Let $T: V \to W$ be a linear transformation. The **kernel** of $T$, also called the **null space**, is defined as

    $$\ker(T) = \{\mathbf{v} \in V : T(\mathbf{v}) = \mathbf{0}\}$$

!!! definition "Definition 5.3 (Image)"
    Let $T: V \to W$ be a linear transformation. The **image** of $T$, also called the **range**, is defined as

    $$\operatorname{im}(T) = \{T(\mathbf{v}) : \mathbf{v} \in V\} = \{w \in W : \exists \mathbf{v} \in V,\; T(\mathbf{v}) = \mathbf{w}\}$$

!!! theorem "Theorem 5.3 (The kernel and image are subspaces)"
    Let $T: V \to W$ be a linear transformation. Then:

    1. $\ker(T)$ is a subspace of $V$.
    2. $\operatorname{im}(T)$ is a subspace of $W$.

??? proof "Proof"
    **(1)** Verification that $\ker(T)$ is a subspace:

    - Nonempty: By Proposition 5.1, $T(\mathbf{0}) = \mathbf{0}$, so $\mathbf{0} \in \ker(T)$.
    - Closure: Let $\mathbf{u}, \mathbf{v} \in \ker(T)$, $c \in \mathbb{F}$. Then $T(\mathbf{u} + \mathbf{v}) = T(\mathbf{u}) + T(\mathbf{v}) = \mathbf{0} + \mathbf{0} = \mathbf{0}$, so $\mathbf{u} + \mathbf{v} \in \ker(T)$. $T(c\mathbf{u}) = cT(\mathbf{u}) = c\mathbf{0} = \mathbf{0}$, so $c\mathbf{u} \in \ker(T)$.

    **(2)** Verification that $\operatorname{im}(T)$ is a subspace:

    - Nonempty: $T(\mathbf{0}) = \mathbf{0} \in \operatorname{im}(T)$.
    - Closure: Let $\mathbf{w}_1 = T(\mathbf{v}_1), \mathbf{w}_2 = T(\mathbf{v}_2) \in \operatorname{im}(T)$, $c \in \mathbb{F}$. Then $\mathbf{w}_1 + \mathbf{w}_2 = T(\mathbf{v}_1) + T(\mathbf{v}_2) = T(\mathbf{v}_1 + \mathbf{v}_2) \in \operatorname{im}(T)$. $c\mathbf{w}_1 = cT(\mathbf{v}_1) = T(c\mathbf{v}_1) \in \operatorname{im}(T)$. $\blacksquare$

!!! definition "Definition 5.4 (Nullity and rank)"
    Let $T: V \to W$ be a linear transformation between finite-dimensional vector spaces.

    - The **nullity** of $T$ is defined as $\operatorname{nullity}(T) = \dim(\ker(T))$
    - The **rank** of $T$ is defined as $\operatorname{rank}(T) = \dim(\operatorname{im}(T))$

!!! theorem "Theorem 5.4 (Characterizations of injectivity and surjectivity)"
    Let $T: V \to W$ be a linear transformation. Then:

    1. $T$ is **injective** if and only if $\ker(T) = \{\mathbf{0}\}$.
    2. $T$ is **surjective** if and only if $\operatorname{im}(T) = W$.

??? proof "Proof"
    **(1)** $(\Rightarrow)$ If $T$ is injective and $\mathbf{v} \in \ker(T)$, then $T(\mathbf{v}) = \mathbf{0} = T(\mathbf{0})$, so by injectivity $\mathbf{v} = \mathbf{0}$.

    $(\Leftarrow)$ If $\ker(T) = \{\mathbf{0}\}$ and $T(\mathbf{u}) = T(\mathbf{v})$, then $T(\mathbf{u} - \mathbf{v}) = \mathbf{0}$, i.e., $\mathbf{u} - \mathbf{v} \in \ker(T) = \{\mathbf{0}\}$, so $\mathbf{u} = \mathbf{v}$.

    **(2)** This follows directly from the definition of image. $\blacksquare$

!!! example "Example 5.5"
    Let $T: \mathbb{R}^3 \to \mathbb{R}^2$ be defined by $T\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} x - y \\ 2y + z \end{pmatrix}$.

    Finding $\ker(T)$: solve $T(\mathbf{v}) = \mathbf{0}$, i.e., $x - y = 0$, $2y + z = 0$. Setting $y = t$, we get $x = t$, $z = -2t$, so

    $$\ker(T) = \left\{ t\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix} : t \in \mathbb{R} \right\}$$

    $\operatorname{nullity}(T) = 1$.

    Finding $\operatorname{im}(T)$: $T(\mathbf{v}) = x\begin{pmatrix} 1 \\ 0 \end{pmatrix} + y\begin{pmatrix} -1 \\ 2 \end{pmatrix} + z\begin{pmatrix} 0 \\ 1 \end{pmatrix}$. Since $\begin{pmatrix} 1 \\ 0 \end{pmatrix}$ and $\begin{pmatrix} -1 \\ 2 \end{pmatrix}$ are linearly independent, $\operatorname{im}(T) = \mathbb{R}^2$ and $\operatorname{rank}(T) = 2$.

    Verification: $\operatorname{nullity}(T) + \operatorname{rank}(T) = 1 + 2 = 3 = \dim(\mathbb{R}^3)$.

---

## 5.3 The Rank-Nullity Theorem

<div class="context-flow" markdown>

**The fundamental equation of finite-dimensional linear algebra**: $\dim V = \dim\ker T + \dim\operatorname{im} T$ → corollary: between spaces of equal dimension, injective $\Leftrightarrow$ surjective $\Leftrightarrow$ bijective

</div>

!!! theorem "Theorem 5.5 (Rank-Nullity Theorem)"
    Let $V$ be a finite-dimensional vector space and $T: V \to W$ a linear transformation. Then

    $$\dim(V) = \operatorname{rank}(T) + \operatorname{nullity}(T) = \dim(\operatorname{im}(T)) + \dim(\ker(T))$$

??? proof "Proof"
    Let $\dim(V) = n$ and $\dim(\ker(T)) = r$. Let $\{\mathbf{u}_1, \ldots, \mathbf{u}_r\}$ be a basis of $\ker(T)$. Extend it to a basis $\{\mathbf{u}_1, \ldots, \mathbf{u}_r, \mathbf{v}_1, \ldots, \mathbf{v}_s\}$ of $V$, where $r + s = n$.

    We show that $\{T(\mathbf{v}_1), \ldots, T(\mathbf{v}_s)\}$ is a basis of $\operatorname{im}(T)$.

    **Spanning**: For any $\mathbf{w} \in \operatorname{im}(T)$, there exists $\mathbf{v} = a_1 \mathbf{u}_1 + \cdots + a_r \mathbf{u}_r + b_1 \mathbf{v}_1 + \cdots + b_s \mathbf{v}_s$ such that $T(\mathbf{v}) = \mathbf{w}$. Since $T(\mathbf{u}_i) = \mathbf{0}$,

    $$\mathbf{w} = T(\mathbf{v}) = b_1 T(\mathbf{v}_1) + \cdots + b_s T(\mathbf{v}_s)$$

    **Linear independence**: Suppose $c_1 T(\mathbf{v}_1) + \cdots + c_s T(\mathbf{v}_s) = \mathbf{0}$. Then $T(c_1 \mathbf{v}_1 + \cdots + c_s \mathbf{v}_s) = \mathbf{0}$, so $c_1 \mathbf{v}_1 + \cdots + c_s \mathbf{v}_s \in \ker(T)$. Therefore there exist scalars $d_1, \ldots, d_r$ such that

    $$c_1 \mathbf{v}_1 + \cdots + c_s \mathbf{v}_s = d_1 \mathbf{u}_1 + \cdots + d_r \mathbf{u}_r$$

    i.e., $d_1 \mathbf{u}_1 + \cdots + d_r \mathbf{u}_r - c_1 \mathbf{v}_1 - \cdots - c_s \mathbf{v}_s = \mathbf{0}$. Since $\{\mathbf{u}_1, \ldots, \mathbf{u}_r, \mathbf{v}_1, \ldots, \mathbf{v}_s\}$ is a basis of $V$ (linearly independent), all coefficients are zero; in particular $c_1 = \cdots = c_s = 0$.

    Therefore $\dim(\operatorname{im}(T)) = s = n - r = \dim(V) - \dim(\ker(T))$. $\blacksquare$

!!! corollary "Corollary 5.1"
    Let $T: V \to W$ be a linear transformation with $\dim(V) = n$ and $\dim(W) = m$. Then:

    1. $\operatorname{rank}(T) \leq \min(n, m)$
    2. If $n = m$, then $T$ is injective $\Leftrightarrow$ $T$ is surjective $\Leftrightarrow$ $T$ is bijective

??? proof "Proof"
    1. $\operatorname{rank}(T) = \dim(\operatorname{im}(T)) \leq \dim(W) = m$, and $\operatorname{rank}(T) = n - \operatorname{nullity}(T) \leq n$.

    2. When $n = m$: $T$ is injective $\Leftrightarrow$ $\operatorname{nullity}(T) = 0$ $\Leftrightarrow$ $\operatorname{rank}(T) = n = m$ $\Leftrightarrow$ $\operatorname{im}(T) = W$ $\Leftrightarrow$ $T$ is surjective. $\blacksquare$

!!! example "Example 5.6"
    Let $T: \mathbb{P}_3 \to \mathbb{P}_2$ be the differentiation operator $T(p) = p'$. Then $\ker(T)$ is the set of constant polynomials (i.e., $\mathbb{P}_0$), so $\dim(\ker(T)) = 1$. By the rank-nullity theorem:

    $$\operatorname{rank}(T) = \dim(\mathbb{P}_3) - \operatorname{nullity}(T) = 4 - 1 = 3 = \dim(\mathbb{P}_2)$$

    Hence $T$ is surjective.

---

## 5.4 Matrix Representation of Linear Transformations

<div class="context-flow" markdown>

**Abstract → computation**: Choosing bases for $V$ and $W$ makes $T$ completely equivalent to an $m \times n$ matrix → $[T(\mathbf{v})]_\mathcal{C} = [T]_\mathcal{B}^\mathcal{C} \cdot [\mathbf{v}]_\mathcal{B}$ → matrix multiplication is the coordinatization of linear transformations

</div>

!!! definition "Definition 5.5 (Matrix representation of a linear transformation)"
    Let $V$ be an $n$-dimensional vector space, $W$ an $m$-dimensional vector space, $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ a basis of $V$, and $\mathcal{C} = \{\mathbf{w}_1, \ldots, \mathbf{w}_m\}$ a basis of $W$. Let $T: V \to W$ be a linear transformation.

    For each $j = 1, \ldots, n$, express $T(\mathbf{v}_j)$ in terms of basis $\mathcal{C}$:

    $$T(\mathbf{v}_j) = a_{1j}\mathbf{w}_1 + a_{2j}\mathbf{w}_2 + \cdots + a_{mj}\mathbf{w}_m$$

    The $m \times n$ matrix

    $$[T]_{\mathcal{B}}^{\mathcal{C}} = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}$$

    is called the **matrix representation** of $T$ with respect to bases $\mathcal{B}$ and $\mathcal{C}$.

!!! theorem "Theorem 5.6 (Matrix representation and coordinate vectors)"
    With the notation of Definition 5.5, if $\mathbf{v} \in V$ has coordinate vector $[\mathbf{v}]_{\mathcal{B}} \in \mathbb{F}^n$ with respect to basis $\mathcal{B}$, then

    $$[T(\mathbf{v})]_{\mathcal{C}} = [T]_{\mathcal{B}}^{\mathcal{C}} \cdot [\mathbf{v}]_{\mathcal{B}}$$

    That is, the action of a linear transformation on coordinate vectors is equivalent to matrix multiplication.

??? proof "Proof"
    Let $\mathbf{v} = c_1 \mathbf{v}_1 + \cdots + c_n \mathbf{v}_n$, so $[\mathbf{v}]_{\mathcal{B}} = (c_1, \ldots, c_n)^T$.

    $$T(\mathbf{v}) = \sum_{j=1}^n c_j T(\mathbf{v}_j) = \sum_{j=1}^n c_j \left(\sum_{i=1}^m a_{ij} \mathbf{w}_i\right) = \sum_{i=1}^m \left(\sum_{j=1}^n a_{ij} c_j\right) \mathbf{w}_i$$

    Therefore the $i$-th component of $[T(\mathbf{v})]_{\mathcal{C}}$ is $\sum_{j=1}^n a_{ij} c_j$, which is precisely the $i$-th component of the product $[T]_{\mathcal{B}}^{\mathcal{C}} \cdot [\mathbf{v}]_{\mathcal{B}}$. $\blacksquare$

!!! example "Example 5.7"
    Let $T: \mathbb{P}_2 \to \mathbb{P}_1$ be the differentiation operator $T(p) = p'$. Take the basis $\mathcal{B} = \{1, t, t^2\}$ of $\mathbb{P}_2$ and the basis $\mathcal{C} = \{1, t\}$ of $\mathbb{P}_1$.

    $$T(1) = 0 = 0 \cdot 1 + 0 \cdot t, \quad T(t) = 1 = 1 \cdot 1 + 0 \cdot t, \quad T(t^2) = 2t = 0 \cdot 1 + 2 \cdot t$$

    Hence

    $$[T]_{\mathcal{B}}^{\mathcal{C}} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$$

    Verification: $p(t) = 3 + 2t - t^2$, $[p]_{\mathcal{B}} = (3, 2, -1)^T$.

    $$[T]_{\mathcal{B}}^{\mathcal{C}} [p]_{\mathcal{B}} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} 3 \\ 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 2 \\ -2 \end{pmatrix}$$

    That is, $T(p) = p' = 2 - 2t$, whose coordinate vector with respect to $\mathcal{C}$ is $(2, -2)^T$. Correct.

---

## 5.5 Change of Basis and Similar Matrices

<div class="context-flow" markdown>

**Same transformation, different matrices**: change of basis → $[T]_{\mathcal{B}'} = P^{-1}[T]_\mathcal{B}P$ → **similar matrices** are different "language descriptions" of the same linear operator → finding a good basis = simplifying the matrix (→ Chapter 6 diagonalization)

</div>

!!! definition "Definition 5.6 (Transition matrix)"
    Let $V$ be an $n$-dimensional vector space, and $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ and $\mathcal{B}' = \{\mathbf{v}_1', \ldots, \mathbf{v}_n'\}$ two bases of $V$. The **transition matrix** (change of basis matrix) $P$ from $\mathcal{B}$ to $\mathcal{B}'$ is the $n \times n$ matrix whose $j$-th column is the coordinate vector of $\mathbf{v}_j$ with respect to $\mathcal{B}'$:

    $$P = \Big([\mathbf{v}_1]_{\mathcal{B}'} \;\; [\mathbf{v}_2]_{\mathcal{B}'} \;\; \cdots \;\; [\mathbf{v}_n]_{\mathcal{B}'}\Big)$$

    Then for any $\mathbf{v} \in V$, $[\mathbf{v}]_{\mathcal{B}'} = P [\mathbf{v}]_{\mathcal{B}}$. (Note: some textbooks define the direction in reverse; be mindful of the convention used.)

!!! definition "Definition 5.7 (Similar matrices)"
    Let $A, B$ be $n \times n$ matrices. If there exists an invertible matrix $P$ such that

    $$B = P^{-1}AP$$

    then $A$ and $B$ are called **similar**, denoted $A \sim B$.

!!! theorem "Theorem 5.7 (Change of basis formula)"
    Let $T: V \to V$ be a linear operator, $\mathcal{B}$ and $\mathcal{B}'$ two bases of $V$, and $P$ the transition matrix from $\mathcal{B}$ to $\mathcal{B}'$. Then

    $$[T]_{\mathcal{B}'} = P^{-1} [T]_{\mathcal{B}} P$$

    That is, the matrix representations of $T$ with respect to different bases are similar.

??? proof "Proof"
    Let $A = [T]_{\mathcal{B}}$, $B = [T]_{\mathcal{B}'}$. For any $\mathbf{v} \in V$,

    $$[T(\mathbf{v})]_{\mathcal{B}} = A [\mathbf{v}]_{\mathcal{B}}$$

    Since $[\mathbf{v}]_{\mathcal{B}'} = P[\mathbf{v}]_{\mathcal{B}}$, i.e., $[\mathbf{v}]_{\mathcal{B}} = P^{-1}[\mathbf{v}]_{\mathcal{B}'}$, we have

    $$[T(\mathbf{v})]_{\mathcal{B}'} = P [T(\mathbf{v})]_{\mathcal{B}} = P A [\mathbf{v}]_{\mathcal{B}} = P A P^{-1} [\mathbf{v}]_{\mathcal{B}'}$$

    By uniqueness, $B = PAP^{-1}$.

    **Note**: The form of the formula depends on the convention for the transition matrix. If the convention $[\mathbf{v}]_{\mathcal{B}} = P[\mathbf{v}]_{\mathcal{B}'}$ is adopted instead, then $B = P^{-1}AP$. Different textbooks use different conventions, but the essence is the same. $\blacksquare$

!!! theorem "Theorem 5.8 (Similarity is an equivalence relation)"
    The similarity relation on matrices satisfies:

    1. **Reflexivity**: $A \sim A$ (take $P = I$)
    2. **Symmetry**: If $A \sim B$, then $B \sim A$
    3. **Transitivity**: If $A \sim B$ and $B \sim C$, then $A \sim C$

??? proof "Proof"
    1. $A = I^{-1}AI$.
    2. If $B = P^{-1}AP$, then $A = PBP^{-1} = (P^{-1})^{-1}B(P^{-1})$; take $Q = P^{-1}$.
    3. If $B = P^{-1}AP$ and $C = Q^{-1}BQ$, then $C = Q^{-1}P^{-1}APQ = (PQ)^{-1}A(PQ)$. $\blacksquare$

!!! example "Example 5.8"
    Let $T: \mathbb{R}^2 \to \mathbb{R}^2$, $T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2x + y \\ 3y \end{pmatrix}$.

    With the standard basis $\mathcal{B} = \{\mathbf{e}_1, \mathbf{e}_2\}$, $[T]_{\mathcal{B}} = A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}$.

    Take the new basis $\mathcal{B}' = \left\{\begin{pmatrix} 1 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 \\ 1 \end{pmatrix}\right\}$. The transition matrix is $P = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$, $P^{-1} = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$.

    $$[T]_{\mathcal{B}'} = P^{-1}AP = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$$

    Under the new basis, the matrix of $T$ becomes diagonal, indicating that $\mathcal{B}'$ is a basis consisting of eigenvectors of $T$.

---

## 5.6 Composition and Inverse of Linear Transformations

<div class="context-flow" markdown>

**Composition ↔ matrix multiplication**: $[S \circ T] = [S][T]$ → invertible transformation $\Leftrightarrow$ bijection $\Leftrightarrow$ invertible matrix (Chapter 2)

</div>

!!! theorem "Theorem 5.9 (The composition of linear transformations is a linear transformation)"
    Let $T: U \to V$ and $S: V \to W$ be linear transformations. Then the composition $S \circ T: U \to W$ is also a linear transformation, and

    $$[S \circ T]_{\mathcal{B}}^{\mathcal{D}} = [S]_{\mathcal{C}}^{\mathcal{D}} \cdot [T]_{\mathcal{B}}^{\mathcal{C}}$$

    where $\mathcal{B}, \mathcal{C}, \mathcal{D}$ are bases of $U, V, W$ respectively.

??? proof "Proof"
    **Linearity**: For any $\mathbf{u}_1, \mathbf{u}_2 \in U$ and $c \in \mathbb{F}$,

    $$(S \circ T)(c\mathbf{u}_1 + \mathbf{u}_2) = S(T(c\mathbf{u}_1 + \mathbf{u}_2)) = S(cT(\mathbf{u}_1) + T(\mathbf{u}_2)) = cS(T(\mathbf{u}_1)) + S(T(\mathbf{u}_2)) = c(S \circ T)(\mathbf{u}_1) + (S \circ T)(\mathbf{u}_2)$$

    **Matrix equation**: For any $\mathbf{u} \in U$,

    $$[(S \circ T)(\mathbf{u})]_{\mathcal{D}} = [S(T(\mathbf{u}))]_{\mathcal{D}} = [S]_{\mathcal{C}}^{\mathcal{D}} [T(\mathbf{u})]_{\mathcal{C}} = [S]_{\mathcal{C}}^{\mathcal{D}} [T]_{\mathcal{B}}^{\mathcal{C}} [\mathbf{u}]_{\mathcal{B}}$$

    By uniqueness, $[S \circ T]_{\mathcal{B}}^{\mathcal{D}} = [S]_{\mathcal{C}}^{\mathcal{D}} [T]_{\mathcal{B}}^{\mathcal{C}}$. $\blacksquare$

!!! definition "Definition 5.8 (Invertible linear transformation)"
    A linear transformation $T: V \to W$ is called **invertible** if there exists a linear transformation $T^{-1}: W \to V$ such that

    $$T^{-1} \circ T = I_V, \quad T \circ T^{-1} = I_W$$

    where $I_V, I_W$ are the identity transformations on $V, W$ respectively.

!!! theorem "Theorem 5.10 (Equivalent conditions for invertibility)"
    Let $T: V \to W$ be a linear transformation with $\dim(V) = \dim(W) = n$. Then the following conditions are equivalent:

    1. $T$ is invertible.
    2. $T$ is bijective.
    3. $\ker(T) = \{\mathbf{0}\}$.
    4. $\operatorname{im}(T) = W$.
    5. $T$ maps a basis of $V$ to a basis of $W$.
    6. $[T]_{\mathcal{B}}^{\mathcal{C}}$ is an invertible matrix (for any bases $\mathcal{B}, \mathcal{C}$).

??? proof "Proof"
    $(1) \Leftrightarrow (2)$: Invertible if and only if bijective — this is a basic fact from the theory of mappings.

    $(2) \Leftrightarrow (3) \Leftrightarrow (4)$: By Theorem 5.4 and Corollary 5.1, when $\dim(V) = \dim(W)$, injectivity, surjectivity, and bijectivity are equivalent.

    $(2) \Rightarrow (5)$: Let $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ be a basis of $V$. By surjectivity, $\{T(\mathbf{v}_1), \ldots, T(\mathbf{v}_n)\}$ spans $W$. By injectivity, if $\sum c_i T(\mathbf{v}_i) = \mathbf{0}$, then $T(\sum c_i \mathbf{v}_i) = \mathbf{0}$, so $\sum c_i \mathbf{v}_i = \mathbf{0}$, hence $c_i = 0$.

    $(5) \Rightarrow (4)$: If $T$ maps a basis to a basis, then $\operatorname{im}(T) \supseteq \operatorname{span}\{T(\mathbf{v}_1), \ldots, T(\mathbf{v}_n)\} = W$.

    $(3) \Leftrightarrow (6)$: $\ker(T) = \{\mathbf{0}\}$ if and only if $[T]_{\mathcal{B}}^{\mathcal{C}} \mathbf{x} = \mathbf{0}$ has only the trivial solution, i.e., the matrix is invertible. $\blacksquare$

!!! example "Example 5.9"
    Let $T: \mathbb{R}^2 \to \mathbb{R}^2$, $T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2x - y \\ x + 3y \end{pmatrix}$. The standard matrix is $A = \begin{pmatrix} 2 & -1 \\ 1 & 3 \end{pmatrix}$.

    $\det(A) = 6 - (-1) = 7 \neq 0$, so $T$ is invertible.

    $$A^{-1} = \frac{1}{7}\begin{pmatrix} 3 & 1 \\ -1 & 2 \end{pmatrix}$$

    Therefore $T^{-1}\begin{pmatrix} x \\ y \end{pmatrix} = \frac{1}{7}\begin{pmatrix} 3x + y \\ -x + 2y \end{pmatrix}$.

---

## 5.7 Isomorphism

<div class="context-flow" markdown>

**Dimension determines everything**: $V \cong W \Leftrightarrow \dim V = \dim W$ → every $n$-dimensional space is isomorphic to $\mathbb{R}^n$ → this is why matrices can solve all finite-dimensional linear algebra problems

</div>

!!! definition "Definition 5.9 (Isomorphism)"
    If there exists an invertible linear transformation $T: V \to W$, then $V$ and $W$ are called **isomorphic**, denoted $V \cong W$. In this case $T$ is called an **isomorphism** from $V$ to $W$.

!!! theorem "Theorem 5.11 (Classification theorem for finite-dimensional vector spaces)"
    Let $V, W$ be finite-dimensional vector spaces over a field $\mathbb{F}$. Then

    $$V \cong W \quad \Longleftrightarrow \quad \dim(V) = \dim(W)$$

    In particular, every $n$-dimensional vector space is isomorphic to $\mathbb{F}^n$.

??? proof "Proof"
    $(\Rightarrow)$ Let $T: V \to W$ be an isomorphism and $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ a basis of $V$. By condition (5) of Theorem 5.10, $\{T(\mathbf{v}_1), \ldots, T(\mathbf{v}_n)\}$ is a basis of $W$, so $\dim(W) = n = \dim(V)$.

    $(\Leftarrow)$ Let $\dim(V) = \dim(W) = n$. Take a basis $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ of $V$ and a basis $\mathcal{C} = \{\mathbf{w}_1, \ldots, \mathbf{w}_n\}$ of $W$. By Theorem 5.2, there exists a unique linear transformation $T$ with $T(\mathbf{v}_i) = \mathbf{w}_i$. Since $T$ maps a basis to a basis, by Theorem 5.10, $T$ is invertible.

    In particular, the coordinate map $\phi_{\mathcal{B}}: V \to \mathbb{F}^n$, $\phi_{\mathcal{B}}(\mathbf{v}) = [\mathbf{v}]_{\mathcal{B}}$ is an isomorphism. $\blacksquare$

!!! note "Note"
    The isomorphism theorem states that in the finite-dimensional setting, the "essence" of a vector space is entirely determined by its dimension. All $n$-dimensional vector spaces — whether $\mathbb{R}^n$, $\mathbb{P}_{n-1}$, or $\mathbb{R}^{m \times k}$ (when $mk = n$) — are "the same" from the viewpoint of linear algebra. This is a powerful aspect of the theory and explains why we can always use matrices and vectors to handle abstract problems.

!!! example "Example 5.10"
    $\mathbb{R}^{2 \times 2}$ (the space of $2 \times 2$ real matrices) is isomorphic to $\mathbb{R}^4$. A concrete isomorphism is

    $$T\begin{pmatrix} a & b \\ c & d \end{pmatrix} = \begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}$$

    This map "flattens" a matrix into a column vector (vectorization).

!!! theorem "Theorem 5.12 (Isomorphism preserves linear structure)"
    Let $T: V \to W$ be an isomorphism. Then:

    1. $T$ maps linearly independent sets to linearly independent sets.
    2. $T$ maps spanning sets to spanning sets.
    3. $T$ maps bases of $V$ to bases of $W$.
    4. For any subspace $U \leq V$, $T(U)$ is a subspace of $W$ and $\dim(T(U)) = \dim(U)$.

??? proof "Proof"
    1. Let $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ be linearly independent. If $\sum c_i T(\mathbf{v}_i) = \mathbf{0}$, then $T(\sum c_i \mathbf{v}_i) = \mathbf{0}$. By injectivity of $T$, $\sum c_i \mathbf{v}_i = \mathbf{0}$, so $c_i = 0$.

    2. Let $\operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\} = V$. For any $\mathbf{w} \in W$, by surjectivity of $T$, there exists $\mathbf{v} = \sum c_i \mathbf{v}_i \in V$ with $T(\mathbf{v}) = \mathbf{w}$, i.e., $\mathbf{w} = \sum c_i T(\mathbf{v}_i)$.

    3. Follows immediately from 1 and 2.

    4. The subspace property of $T(U)$ is guaranteed by the linearity of $T$. The restriction of $T$ to $U$ remains injective, so $\ker(T|_U) = \{\mathbf{0}\}$, and by the rank-nullity theorem, $\dim(T(U)) = \dim(U)$. $\blacksquare$

---

## 5.8 Invariant Subspaces

<div class="context-flow" markdown>

**Theoretical basis for diagonalization**: $T(U) \subseteq U$ → $T$ can be "studied independently" on an invariant subspace → direct sum decomposition $V = U_1 \oplus U_2$ corresponds to a block diagonal matrix → one-dimensional invariant subspace = eigenvector (→ Chapter 6)

</div>

!!! definition "Definition 5.10 (Invariant subspace)"
    Let $T: V \to V$ be a linear operator and $U$ a subspace of $V$. If $T(U) \subseteq U$ (i.e., for every $\mathbf{u} \in U$, $T(\mathbf{u}) \in U$), then $U$ is called an **invariant subspace** of $T$.

!!! proposition "Proposition 5.2 (Trivial invariant subspaces)"
    Let $T: V \to V$ be a linear operator. Then:

    1. $\{\mathbf{0}\}$ and $V$ are invariant subspaces of $T$.
    2. $\ker(T)$ is an invariant subspace of $T$.
    3. $\operatorname{im}(T)$ is an invariant subspace of $T$ (note: $\operatorname{im}(T) \subseteq V$ since $T$ is an operator).

??? proof "Proof"
    1. Clearly $T(\{\mathbf{0}\}) = \{\mathbf{0}\} \subseteq \{\mathbf{0}\}$ and $T(V) \subseteq V$.
    2. If $\mathbf{u} \in \ker(T)$, then $T(\mathbf{u}) = \mathbf{0} \in \ker(T)$.
    3. If $\mathbf{w} \in \operatorname{im}(T)$, then $\mathbf{w} = T(\mathbf{v})$ for some $\mathbf{v} \in V$, and $T(\mathbf{w}) = T(T(\mathbf{v})) \in \operatorname{im}(T)$. $\blacksquare$

!!! theorem "Theorem 5.13 (Invariant subspaces and block diagonal matrices)"
    Let $T: V \to V$ be a linear operator with $\dim(V) = n$. If $V = U_1 \oplus U_2$ (direct sum decomposition), $\dim(U_1) = k$, $\dim(U_2) = n - k$, and both $U_1, U_2$ are invariant subspaces of $T$, then there exists a basis $\mathcal{B}$ of $V$ such that

    $$[T]_{\mathcal{B}} = \begin{pmatrix} A_1 & 0 \\ 0 & A_2 \end{pmatrix}$$

    where $A_1$ is a $k \times k$ matrix (the matrix representation of $T|_{U_1}$) and $A_2$ is an $(n-k) \times (n-k)$ matrix (the matrix representation of $T|_{U_2}$).

??? proof "Proof"
    Take a basis $\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$ of $U_1$ and a basis $\{\mathbf{u}_{k+1}, \ldots, \mathbf{u}_n\}$ of $U_2$; together they form a basis $\mathcal{B}$ of $V$.

    Since $U_1$ is an invariant subspace, for $j \leq k$, $T(\mathbf{u}_j) \in U_1$, so $T(\mathbf{u}_j)$ can be expressed as a linear combination of $\mathbf{u}_1, \ldots, \mathbf{u}_k$ with zero coefficients for $\mathbf{u}_{k+1}, \ldots, \mathbf{u}_n$.

    Similarly, for $j > k$, $T(\mathbf{u}_j) \in U_2$, with zero coefficients for $\mathbf{u}_1, \ldots, \mathbf{u}_k$.

    Therefore the matrix has block diagonal form. $\blacksquare$

!!! theorem "Theorem 5.14 (Eigenvalues and invariant subspaces)"
    Let $T: V \to V$ be a linear operator and $\lambda$ a scalar. Then the following are equivalent:

    1. $\lambda$ is an eigenvalue of $T$ (i.e., there exists a nonzero $\mathbf{v}$ with $T(\mathbf{v}) = \lambda\mathbf{v}$).
    2. $T$ has a one-dimensional invariant subspace $U = \operatorname{span}\{\mathbf{v}\}$ on which $T|_U = \lambda \cdot I_U$.

??? proof "Proof"
    $(1) \Rightarrow (2)$: Suppose $T(\mathbf{v}) = \lambda\mathbf{v}$ with $\mathbf{v} \neq \mathbf{0}$. Let $U = \operatorname{span}\{\mathbf{v}\}$. For any $c\mathbf{v} \in U$, $T(c\mathbf{v}) = cT(\mathbf{v}) = c\lambda\mathbf{v} = \lambda(c\mathbf{v}) \in U$.

    $(2) \Rightarrow (1)$: Let $U = \operatorname{span}\{\mathbf{v}\}$ be a one-dimensional invariant subspace. Then $T(\mathbf{v}) \in U$, i.e., $T(\mathbf{v}) = \lambda\mathbf{v}$ for some scalar $\lambda$. Since $\mathbf{v} \neq \mathbf{0}$, $\lambda$ is an eigenvalue of $T$. $\blacksquare$

!!! example "Example 5.11"
    Let $V = \mathbb{R}^2$ and $T$ be reflection about the $x$-axis: $T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} x \\ -y \end{pmatrix}$.

    - $U_1 = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 0 \end{pmatrix}\right\}$ (the $x$-axis) is an invariant subspace of $T$, with $T|_{U_1} = I$.
    - $U_2 = \operatorname{span}\left\{\begin{pmatrix} 0 \\ 1 \end{pmatrix}\right\}$ (the $y$-axis) is an invariant subspace of $T$, with $T|_{U_2} = -I$.
    - $V = U_1 \oplus U_2$, and with respect to the basis $\{\mathbf{e}_1, \mathbf{e}_2\}$, $[T] = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$, a diagonal matrix.

!!! example "Example 5.12"
    Let $T: \mathbb{R}^2 \to \mathbb{R}^2$ be counterclockwise rotation by $\pi/2$: $T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} -y \\ x \end{pmatrix}$.

    $T$ has no one-dimensional invariant subspace in $\mathbb{R}^2$ (since a $90°$ rotation does not map any line to itself), so $T$ has no eigenvalue over $\mathbb{R}$. Over $\mathbb{C}$, however, the eigenvalues of $T$ are $\pm i$.

    This shows that the existence of invariant subspaces (and eigenvalues) may depend on the choice of the underlying field.

---

## Chapter Summary

This chapter studied, from an abstract perspective, maps between vector spaces that preserve linear structure — linear transformations. The core topics include:

1. **Definition of linear transformations**: mappings that preserve addition and scalar multiplication, uniquely determined by their values on a basis.
2. **Kernel and image**: characterize how much information a linear transformation "loses" and "covers," respectively.
3. **Rank-nullity theorem**: $\dim(V) = \operatorname{rank}(T) + \operatorname{nullity}(T)$, one of the most fundamental equations in finite-dimensional linear algebra.
4. **Matrix representation**: a linear transformation with respect to given bases is equivalent to matrix multiplication, bridging abstract theory and computation.
5. **Change of basis and similarity**: the matrices of the same linear operator with respect to different bases are related by $B = P^{-1}AP$.
6. **Composition and inverse**: composition corresponds to matrix multiplication; invertible if and only if bijective.
7. **Isomorphism**: finite-dimensional vector spaces are completely classified by dimension.
8. **Invariant subspaces**: provide the theoretical foundation for diagonalization and block structure of matrices, leading directly to eigenvalue theory.
