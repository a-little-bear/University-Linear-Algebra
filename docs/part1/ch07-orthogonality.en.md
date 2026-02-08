# Chapter 7  Orthogonality and Least Squares

<div class="context-flow" markdown>

**Prerequisites**: Chapter 4 subspaces/four fundamental subspaces · Chapter 6 real symmetric matrices/spectral theorem · **Chapter arc**: inner product/norm → Cauchy-Schwarz → orthogonal sets → orthogonal complements → orthogonal projection → Gram-Schmidt → orthogonal matrices → least squares → QR decomposition
**One-line essence**: The inner product endows a vector space with **geometric structure** (length, angle, perpendicularity) — orthogonal projection is the algebraic expression of "best approximation," and QR decomposition is its algorithmic realization

</div>

Orthogonality is the core concept where geometric intuition and algebraic theory meet in linear algebra. In $\mathbb{R}^n$, the inner product gives vectors the notions of "length" and "angle," enabling us to discuss orthogonality, projection, best approximation, and related problems. This chapter proceeds from inner products and norms to systematically study orthogonal vectors, orthogonal complements, orthogonal projection, the Gram-Schmidt orthogonalization process, orthogonal matrices, the method of least squares, and QR decomposition. These tools carry deep theoretical significance and find extremely broad applications in numerical linear algebra, statistics, signal processing, and data science.

---

## 7.1 Inner Products and Norms in $\mathbb{R}^n$

<div class="context-flow" markdown>

**Injecting geometric structure**: Chapters 1-6 are purely algebraic (addition/multiplication/determinants) → the inner product $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^T\mathbf{v}$ introduces length and angle → the Cauchy-Schwarz inequality is the cornerstone of everything

</div>

!!! definition "Definition 7.1 (Inner product)"
    The **inner product** (also called **dot product**) of vectors $\mathbf{u} = (u_1, \ldots, u_n)^T$ and $\mathbf{v} = (v_1, \ldots, v_n)^T$ in $\mathbb{R}^n$ is defined as

    $$\langle \mathbf{u}, \mathbf{v} \rangle = \mathbf{u} \cdot \mathbf{v} = \mathbf{u}^T\mathbf{v} = \sum_{i=1}^n u_i v_i$$

!!! definition "Definition 7.2 (Norm)"
    The **norm** (or **length**) of a vector $\mathbf{v} \in \mathbb{R}^n$ is defined as

    $$\|\mathbf{v}\| = \sqrt{\langle \mathbf{v}, \mathbf{v} \rangle} = \sqrt{\mathbf{v}^T\mathbf{v}} = \sqrt{\sum_{i=1}^n v_i^2}$$

    If $\|\mathbf{v}\| = 1$, then $\mathbf{v}$ is called a **unit vector**. The process of converting a nonzero vector $\mathbf{v}$ to a unit vector in the same direction is called **normalization**: $\hat{\mathbf{v}} = \frac{\mathbf{v}}{\|\mathbf{v}\|}$.

!!! definition "Definition 7.3 (Distance)"
    The **distance** between two points $\mathbf{u}$ and $\mathbf{v}$ in $\mathbb{R}^n$ is defined as

    $$d(\mathbf{u}, \mathbf{v}) = \|\mathbf{u} - \mathbf{v}\|$$

!!! proposition "Proposition 7.1 (Basic properties of the inner product)"
    For any $\mathbf{u}, \mathbf{v}, \mathbf{w} \in \mathbb{R}^n$ and $c \in \mathbb{R}$:

    1. $\langle \mathbf{u}, \mathbf{v} \rangle = \langle \mathbf{v}, \mathbf{u} \rangle$ (symmetry)
    2. $\langle \mathbf{u} + \mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{u}, \mathbf{w} \rangle + \langle \mathbf{v}, \mathbf{w} \rangle$ (additivity)
    3. $\langle c\mathbf{u}, \mathbf{v} \rangle = c\langle \mathbf{u}, \mathbf{v} \rangle$ (homogeneity)
    4. $\langle \mathbf{v}, \mathbf{v} \rangle \geq 0$, with equality if and only if $\mathbf{v} = \mathbf{0}$ (positive definiteness)

!!! theorem "Theorem 7.1 (Cauchy-Schwarz inequality)"
    For any $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$,

    $$|\langle \mathbf{u}, \mathbf{v} \rangle| \leq \|\mathbf{u}\| \cdot \|\mathbf{v}\|$$

    Equality holds if and only if $\mathbf{u}$ and $\mathbf{v}$ are linearly dependent (i.e., one is a scalar multiple of the other).

??? proof "Proof"
    If $\mathbf{v} = \mathbf{0}$, both sides are $0$ and the inequality holds. Assume $\mathbf{v} \neq \mathbf{0}$.

    For any $t \in \mathbb{R}$, consider

    $$0 \leq \|\mathbf{u} - t\mathbf{v}\|^2 = \langle \mathbf{u} - t\mathbf{v}, \mathbf{u} - t\mathbf{v} \rangle = \|\mathbf{u}\|^2 - 2t\langle \mathbf{u}, \mathbf{v} \rangle + t^2\|\mathbf{v}\|^2$$

    This is a quadratic function in $t$, minimized at $t = \frac{\langle \mathbf{u}, \mathbf{v} \rangle}{\|\mathbf{v}\|^2}$. Substituting:

    $$0 \leq \|\mathbf{u}\|^2 - \frac{\langle \mathbf{u}, \mathbf{v} \rangle^2}{\|\mathbf{v}\|^2}$$

    That is, $\langle \mathbf{u}, \mathbf{v} \rangle^2 \leq \|\mathbf{u}\|^2 \|\mathbf{v}\|^2$. Taking square roots gives $|\langle \mathbf{u}, \mathbf{v} \rangle| \leq \|\mathbf{u}\| \cdot \|\mathbf{v}\|$.

    Equality holds $\Leftrightarrow$ $\|\mathbf{u} - t\mathbf{v}\|^2 = 0$ $\Leftrightarrow$ $\mathbf{u} = t\mathbf{v}$ $\Leftrightarrow$ $\mathbf{u}, \mathbf{v}$ are linearly dependent. $\blacksquare$

!!! theorem "Theorem 7.2 (Triangle inequality)"
    For any $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$,

    $$\|\mathbf{u} + \mathbf{v}\| \leq \|\mathbf{u}\| + \|\mathbf{v}\|$$

??? proof "Proof"
    $$\|\mathbf{u} + \mathbf{v}\|^2 = \langle \mathbf{u} + \mathbf{v}, \mathbf{u} + \mathbf{v} \rangle = \|\mathbf{u}\|^2 + 2\langle \mathbf{u}, \mathbf{v} \rangle + \|\mathbf{v}\|^2$$

    By the Cauchy-Schwarz inequality:

    $$\leq \|\mathbf{u}\|^2 + 2\|\mathbf{u}\|\|\mathbf{v}\| + \|\mathbf{v}\|^2 = (\|\mathbf{u}\| + \|\mathbf{v}\|)^2$$

    Taking square roots (both sides are nonnegative) yields the result. $\blacksquare$

!!! example "Example 7.1"
    Let $\mathbf{u} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$, $\mathbf{v} = \begin{pmatrix} 4 \\ -1 \\ 2 \end{pmatrix}$.

    $\langle \mathbf{u}, \mathbf{v} \rangle = 1 \cdot 4 + 2 \cdot (-1) + 3 \cdot 2 = 4 - 2 + 6 = 8$.

    $\|\mathbf{u}\| = \sqrt{1 + 4 + 9} = \sqrt{14}$, $\|\mathbf{v}\| = \sqrt{16 + 1 + 4} = \sqrt{21}$.

    Verifying Cauchy-Schwarz: $|8| = 8 \leq \sqrt{14} \cdot \sqrt{21} = \sqrt{294} \approx 17.15$. Holds.

    The angle between the two vectors: $\cos\theta = \frac{\langle \mathbf{u}, \mathbf{v} \rangle}{\|\mathbf{u}\|\|\mathbf{v}\|} = \frac{8}{\sqrt{294}}$, $\theta = \arccos\frac{8}{\sqrt{294}} \approx 62.2°$.

---

## 7.2 Orthogonal Vectors and Orthogonal Sets

<div class="context-flow" markdown>

**Orthogonal = "most independent"**: $\langle\mathbf{u},\mathbf{v}\rangle=0$ → orthogonal sets are automatically linearly independent (an inner-product upgrade of the Chapter 4 proof technique) → orthogonal bases simplify the coordinate formula to inner-product computations

</div>

!!! definition "Definition 7.4 (Orthogonal)"
    Two vectors $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$ are called **orthogonal** if $\langle \mathbf{u}, \mathbf{v} \rangle = 0$, denoted $\mathbf{u} \perp \mathbf{v}$.

!!! definition "Definition 7.5 (Orthogonal set and orthogonal basis)"
    A set of vectors $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ is called an **orthogonal set** if any two distinct vectors in it are orthogonal:

    $$\langle \mathbf{v}_i, \mathbf{v}_j \rangle = 0, \quad \forall\, i \neq j$$

    If furthermore each vector is a unit vector ($\|\mathbf{v}_i\| = 1$), the set is called an **orthonormal set**.

    If an orthogonal set (or orthonormal set) also forms a basis of some subspace, it is called an **orthogonal basis** (or **orthonormal basis**, respectively).

!!! theorem "Theorem 7.3 (An orthogonal set is linearly independent)"
    If $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ is an orthogonal set of nonzero vectors, then the set is linearly independent.

??? proof "Proof"
    Suppose $c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_k\mathbf{v}_k = \mathbf{0}$. For any $j$, take the inner product of both sides with $\mathbf{v}_j$:

    $$0 = \langle \mathbf{0}, \mathbf{v}_j \rangle = \left\langle \sum_{i=1}^k c_i\mathbf{v}_i, \mathbf{v}_j \right\rangle = \sum_{i=1}^k c_i \langle \mathbf{v}_i, \mathbf{v}_j \rangle = c_j \langle \mathbf{v}_j, \mathbf{v}_j \rangle = c_j \|\mathbf{v}_j\|^2$$

    Since $\mathbf{v}_j \neq \mathbf{0}$, $\|\mathbf{v}_j\|^2 > 0$, so $c_j = 0$. $\blacksquare$

!!! theorem "Theorem 7.4 (Coordinate formula with respect to an orthogonal basis)"
    Let $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ be an orthogonal basis of $V$. Then for any $\mathbf{w} \in V$,

    $$\mathbf{w} = \sum_{i=1}^n \frac{\langle \mathbf{w}, \mathbf{v}_i \rangle}{\langle \mathbf{v}_i, \mathbf{v}_i \rangle}\mathbf{v}_i = \sum_{i=1}^n \frac{\langle \mathbf{w}, \mathbf{v}_i \rangle}{\|\mathbf{v}_i\|^2}\mathbf{v}_i$$

    If the basis is orthonormal, the formula simplifies to $\mathbf{w} = \sum_{i=1}^n \langle \mathbf{w}, \mathbf{v}_i \rangle \mathbf{v}_i$.

??? proof "Proof"
    Let $\mathbf{w} = c_1\mathbf{v}_1 + \cdots + c_n\mathbf{v}_n$. Taking the inner product of both sides with $\mathbf{v}_j$:

    $$\langle \mathbf{w}, \mathbf{v}_j \rangle = c_j \langle \mathbf{v}_j, \mathbf{v}_j \rangle = c_j\|\mathbf{v}_j\|^2$$

    Hence $c_j = \frac{\langle \mathbf{w}, \mathbf{v}_j \rangle}{\|\mathbf{v}_j\|^2}$. $\blacksquare$

!!! example "Example 7.2"
    In $\mathbb{R}^3$, let $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$, $\mathbf{v}_2 = \begin{pmatrix} -1 \\ 1 \\ 0 \end{pmatrix}$, $\mathbf{v}_3 = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$.

    Checking orthogonality: $\langle \mathbf{v}_1, \mathbf{v}_2 \rangle = -1 + 1 + 0 = 0$, $\langle \mathbf{v}_1, \mathbf{v}_3 \rangle = 0$, $\langle \mathbf{v}_2, \mathbf{v}_3 \rangle = 0$.

    So $\{\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3\}$ is an orthogonal basis of $\mathbb{R}^3$.

    Expressing $\mathbf{w} = \begin{pmatrix} 3 \\ 5 \\ 7 \end{pmatrix}$ in this orthogonal basis:

    $$c_1 = \frac{\langle \mathbf{w}, \mathbf{v}_1 \rangle}{\|\mathbf{v}_1\|^2} = \frac{3 + 5}{2} = 4, \quad c_2 = \frac{\langle \mathbf{w}, \mathbf{v}_2 \rangle}{\|\mathbf{v}_2\|^2} = \frac{-3 + 5}{2} = 1, \quad c_3 = \frac{\langle \mathbf{w}, \mathbf{v}_3 \rangle}{\|\mathbf{v}_3\|^2} = \frac{7}{1} = 7$$

    Verification: $4\mathbf{v}_1 + 1\mathbf{v}_2 + 7\mathbf{v}_3 = \begin{pmatrix} 4 - 1 \\ 4 + 1 \\ 7 \end{pmatrix} = \begin{pmatrix} 3 \\ 5 \\ 7 \end{pmatrix} = \mathbf{w}$. Correct.

---

## 7.3 Orthogonal Complements

<div class="context-flow" markdown>

**Orthogonal decomposition**: $\mathbb{R}^n = W \oplus W^\perp$ (the geometric realization of Chapter 4 direct sums) → the row space and null space are orthogonal complements $\operatorname{Row}(A)^\perp = \operatorname{Null}(A)$ → the complete geometric picture of the four fundamental subspaces

</div>

!!! definition "Definition 7.6 (Orthogonal complement)"
    Let $W$ be a subspace of $\mathbb{R}^n$. The **orthogonal complement** of $W$ is defined as

    $$W^\perp = \{\mathbf{v} \in \mathbb{R}^n : \langle \mathbf{v}, \mathbf{w} \rangle = 0, \;\forall\, \mathbf{w} \in W\}$$

!!! theorem "Theorem 7.5 (The orthogonal complement is a subspace)"
    $W^\perp$ is a subspace of $\mathbb{R}^n$.

??? proof "Proof"
    $\mathbf{0} \in W^\perp$ (since $\langle \mathbf{0}, \mathbf{w} \rangle = 0$ for all $\mathbf{w}$).

    Let $\mathbf{u}, \mathbf{v} \in W^\perp$, $c \in \mathbb{R}$. For any $\mathbf{w} \in W$:

    $\langle \mathbf{u} + \mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{u}, \mathbf{w} \rangle + \langle \mathbf{v}, \mathbf{w} \rangle = 0 + 0 = 0$, so $\mathbf{u} + \mathbf{v} \in W^\perp$.

    $\langle c\mathbf{u}, \mathbf{w} \rangle = c\langle \mathbf{u}, \mathbf{w} \rangle = 0$, so $c\mathbf{u} \in W^\perp$. $\blacksquare$

!!! theorem "Theorem 7.6 (Orthogonal decomposition theorem)"
    Let $W$ be a subspace of $\mathbb{R}^n$. Then

    $$\mathbb{R}^n = W \oplus W^\perp$$

    That is, every vector $\mathbf{v} \in \mathbb{R}^n$ can be uniquely decomposed as

    $$\mathbf{v} = \mathbf{w} + \mathbf{w}^\perp, \quad \mathbf{w} \in W, \;\; \mathbf{w}^\perp \in W^\perp$$

??? proof "Proof"
    **Existence**: Take an orthonormal basis $\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$ of $W$. Define

    $$\mathbf{w} = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{u}_i \rangle \mathbf{u}_i \in W, \quad \mathbf{w}^\perp = \mathbf{v} - \mathbf{w}$$

    We verify that $\mathbf{w}^\perp \in W^\perp$: for any $j$,

    $$\langle \mathbf{w}^\perp, \mathbf{u}_j \rangle = \langle \mathbf{v} - \mathbf{w}, \mathbf{u}_j \rangle = \langle \mathbf{v}, \mathbf{u}_j \rangle - \sum_{i=1}^k \langle \mathbf{v}, \mathbf{u}_i \rangle \langle \mathbf{u}_i, \mathbf{u}_j \rangle = \langle \mathbf{v}, \mathbf{u}_j \rangle - \langle \mathbf{v}, \mathbf{u}_j \rangle = 0$$

    Therefore $\mathbf{w}^\perp$ is orthogonal to the basis of $W$, i.e., $\mathbf{w}^\perp \in W^\perp$.

    **Uniqueness**: If $\mathbf{v} = \mathbf{w}_1 + \mathbf{w}_1^\perp = \mathbf{w}_2 + \mathbf{w}_2^\perp$, then $\mathbf{w}_1 - \mathbf{w}_2 = \mathbf{w}_2^\perp - \mathbf{w}_1^\perp$. The left side belongs to $W$ and the right side to $W^\perp$, so $\mathbf{w}_1 - \mathbf{w}_2 \in W \cap W^\perp$. For $\mathbf{x} \in W \cap W^\perp$, $\langle \mathbf{x}, \mathbf{x} \rangle = 0$, so $\mathbf{x} = \mathbf{0}$. Hence $\mathbf{w}_1 = \mathbf{w}_2$ and $\mathbf{w}_1^\perp = \mathbf{w}_2^\perp$.

    **Dimension relation**: $\dim(W) + \dim(W^\perp) = n$. $\blacksquare$

!!! corollary "Corollary 7.1"
    $(W^\perp)^\perp = W$.

??? proof "Proof"
    Clearly $W \subseteq (W^\perp)^\perp$ (every vector in $W$ is orthogonal to every vector in $W^\perp$). Moreover, $\dim((W^\perp)^\perp) = n - \dim(W^\perp) = n - (n - \dim W) = \dim W$. A subspace containment with equal dimensions implies equality. $\blacksquare$

!!! theorem "Theorem 7.7 (Orthogonal relationship between row space and null space)"
    Let $A$ be an $m \times n$ matrix. Then:

    1. $\operatorname{Row}(A)^\perp = \operatorname{Null}(A)$, i.e., the orthogonal complement of the row space is the null space.
    2. $\operatorname{Col}(A)^\perp = \operatorname{Null}(A^T)$, i.e., the orthogonal complement of the column space is the left null space.

??? proof "Proof"
    **1.** $\mathbf{x} \in \operatorname{Null}(A)$ $\Leftrightarrow$ $A\mathbf{x} = \mathbf{0}$ $\Leftrightarrow$ the inner product of each row of $A$ with $\mathbf{x}$ is zero $\Leftrightarrow$ $\mathbf{x}$ is orthogonal to every row of $A$ $\Leftrightarrow$ $\mathbf{x} \in \operatorname{Row}(A)^\perp$.

    **2.** Applying 1 to $A^T$: $\operatorname{Row}(A^T)^\perp = \operatorname{Null}(A^T)$, i.e., $\operatorname{Col}(A)^\perp = \operatorname{Null}(A^T)$. $\blacksquare$

!!! example "Example 7.3"
    Let $W = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}, \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}\right\} \subseteq \mathbb{R}^3$. Find $W^\perp$.

    $W^\perp = \operatorname{Null}(A)$ where $A = \begin{pmatrix} 1 & 1 & 1 \\ 1 & -1 & 0 \end{pmatrix}$.

    Row reduction: $\begin{pmatrix} 1 & 1 & 1 \\ 0 & -2 & -1 \end{pmatrix} \to \begin{pmatrix} 1 & 0 & 1/2 \\ 0 & 1 & 1/2 \end{pmatrix}$.

    $x_3 = t$ is free, $x_1 = -t/2$, $x_2 = -t/2$.

    $$W^\perp = \operatorname{span}\left\{\begin{pmatrix} -1 \\ -1 \\ 2 \end{pmatrix}\right\} = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}\right\}$$

    Verification: $\dim(W) + \dim(W^\perp) = 2 + 1 = 3 = \dim(\mathbb{R}^3)$.

---

## 7.4 Orthogonal Projection

<div class="context-flow" markdown>

**Best approximation**: $\operatorname{proj}_W(\mathbf{v})$ is the closest point in $W$ to $\mathbf{v}$ → projection matrix $P = A(A^TA)^{-1}A^T$ satisfies $P^2=P$, $P^T=P$ → the geometric foundation of least squares

</div>

!!! definition "Definition 7.7 (Orthogonal projection)"
    Let $W$ be a subspace of $\mathbb{R}^n$. The **orthogonal projection** of vector $\mathbf{v}$ onto $W$ is the $\mathbf{w}$ component in the orthogonal decomposition $\mathbf{v} = \mathbf{w} + \mathbf{w}^\perp$, denoted

    $$\operatorname{proj}_W(\mathbf{v}) = \mathbf{w}$$

!!! theorem "Theorem 7.8 (Orthogonal projection formula)"
    Let $\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$ be an orthogonal basis of $W$. Then

    $$\operatorname{proj}_W(\mathbf{v}) = \sum_{i=1}^k \frac{\langle \mathbf{v}, \mathbf{u}_i \rangle}{\langle \mathbf{u}_i, \mathbf{u}_i \rangle}\mathbf{u}_i$$

    If $\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$ is an orthonormal basis, then

    $$\operatorname{proj}_W(\mathbf{v}) = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{u}_i \rangle \mathbf{u}_i$$

??? proof "Proof"
    This follows directly from the construction in the proof of Theorem 7.6. The key verification is that $\mathbf{v} - \operatorname{proj}_W(\mathbf{v}) \in W^\perp$, which was already established in that proof. $\blacksquare$

!!! theorem "Theorem 7.9 (Best approximation theorem)"
    Let $W$ be a subspace of $\mathbb{R}^n$ and $\mathbf{v} \in \mathbb{R}^n$. Then $\operatorname{proj}_W(\mathbf{v})$ is the point in $W$ closest to $\mathbf{v}$. That is, for any $\mathbf{w} \in W$,

    $$\|\mathbf{v} - \operatorname{proj}_W(\mathbf{v})\| \leq \|\mathbf{v} - \mathbf{w}\|$$

    Equality holds if and only if $\mathbf{w} = \operatorname{proj}_W(\mathbf{v})$.

??? proof "Proof"
    Let $\hat{\mathbf{v}} = \operatorname{proj}_W(\mathbf{v})$. For any $\mathbf{w} \in W$,

    $$\|\mathbf{v} - \mathbf{w}\|^2 = \|(\mathbf{v} - \hat{\mathbf{v}}) + (\hat{\mathbf{v}} - \mathbf{w})\|^2$$

    Note that $\mathbf{v} - \hat{\mathbf{v}} \in W^\perp$ and $\hat{\mathbf{v}} - \mathbf{w} \in W$ are orthogonal. By the Pythagorean theorem:

    $$= \|\mathbf{v} - \hat{\mathbf{v}}\|^2 + \|\hat{\mathbf{v}} - \mathbf{w}\|^2 \geq \|\mathbf{v} - \hat{\mathbf{v}}\|^2$$

    Equality holds $\Leftrightarrow$ $\|\hat{\mathbf{v}} - \mathbf{w}\| = 0$ $\Leftrightarrow$ $\mathbf{w} = \hat{\mathbf{v}}$. $\blacksquare$

!!! definition "Definition 7.8 (Projection matrix)"
    Let the columns of matrix $A$ span the subspace $W = \operatorname{Col}(A)$. If the columns of $A$ are linearly independent, then the orthogonal projection matrix from $\mathbb{R}^n$ onto $W$ is

    $$P = A(A^TA)^{-1}A^T$$

    This matrix satisfies $\operatorname{proj}_W(\mathbf{v}) = P\mathbf{v}$.

!!! proposition "Proposition 7.2 (Properties of the projection matrix)"
    The projection matrix $P = A(A^TA)^{-1}A^T$ satisfies:

    1. $P^2 = P$ (idempotence)
    2. $P^T = P$ (symmetry)
    3. $\operatorname{rank}(P) = \operatorname{rank}(A) = \dim(W)$

??? proof "Proof"
    1. $P^2 = A(A^TA)^{-1}A^T \cdot A(A^TA)^{-1}A^T = A(A^TA)^{-1}(A^TA)(A^TA)^{-1}A^T = A(A^TA)^{-1}A^T = P$.

    2. $P^T = (A(A^TA)^{-1}A^T)^T = A((A^TA)^{-1})^T A^T = A((A^TA)^T)^{-1}A^T = A(A^TA)^{-1}A^T = P$ (since $A^TA$ is symmetric, so is its inverse).

    3. The column space of $P$ equals $\operatorname{Col}(A)$, so $\operatorname{rank}(P) = \operatorname{rank}(A)$. $\blacksquare$

!!! example "Example 7.4"
    Project $\mathbf{v} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$ onto $W = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}\right\}$.

    Note that these two vectors are not orthogonal (their inner product is $0 + 0 + 1 = 1 \neq 0$), so we use the projection matrix formula.

    $A = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{pmatrix}$, $A^TA = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$, $(A^TA)^{-1} = \frac{1}{3}\begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$.

    $$\hat{\mathbf{v}} = A(A^TA)^{-1}A^T\mathbf{v} = A \cdot \frac{1}{3}\begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}\begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$$

    $A^T\mathbf{v} = \begin{pmatrix} 4 \\ 5 \end{pmatrix}$, $(A^TA)^{-1}A^T\mathbf{v} = \frac{1}{3}\begin{pmatrix} 3 \\ 6 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$.

    $$\hat{\mathbf{v}} = A\begin{pmatrix} 1 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$$

    Coincidentally, $\mathbf{v}$ itself lies in $W$! Verification: $\mathbf{v} = 1 \cdot \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix} + 2 \cdot \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$.

---

## 7.5 Gram-Schmidt Orthogonalization

<div class="context-flow" markdown>

**Algorithm: arbitrary basis → orthogonal basis**: successively subtract projection components along existing directions → $\mathbf{v}_j = \mathbf{x}_j - \operatorname{proj}_{W_{j-1}}(\mathbf{x}_j)$ → directly produces QR decomposition

</div>

!!! theorem "Theorem 7.10 (Gram-Schmidt orthogonalization process)"
    Let $\{\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_k\}$ be a basis of a subspace $W$. Define:

    $$\mathbf{v}_1 = \mathbf{x}_1$$

    $$\mathbf{v}_2 = \mathbf{x}_2 - \frac{\langle \mathbf{x}_2, \mathbf{v}_1 \rangle}{\langle \mathbf{v}_1, \mathbf{v}_1 \rangle}\mathbf{v}_1$$

    $$\mathbf{v}_3 = \mathbf{x}_3 - \frac{\langle \mathbf{x}_3, \mathbf{v}_1 \rangle}{\langle \mathbf{v}_1, \mathbf{v}_1 \rangle}\mathbf{v}_1 - \frac{\langle \mathbf{x}_3, \mathbf{v}_2 \rangle}{\langle \mathbf{v}_2, \mathbf{v}_2 \rangle}\mathbf{v}_2$$

    In general,

    $$\mathbf{v}_j = \mathbf{x}_j - \sum_{i=1}^{j-1} \frac{\langle \mathbf{x}_j, \mathbf{v}_i \rangle}{\langle \mathbf{v}_i, \mathbf{v}_i \rangle}\mathbf{v}_i, \quad j = 2, 3, \ldots, k$$

    Then $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\}$ is an orthogonal basis of $W$. Further normalizing $\mathbf{q}_i = \frac{\mathbf{v}_i}{\|\mathbf{v}_i\|}$ yields an orthonormal basis.

    Moreover, for each $j$, $\operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_j\} = \operatorname{span}\{\mathbf{x}_1, \ldots, \mathbf{x}_j\}$.

??? proof "Proof"
    By induction on $j$.

    **Base case**: $j = 1$: $\{\mathbf{v}_1\} = \{\mathbf{x}_1\}$ is a single-element set, trivially orthogonal.

    **Inductive step**: Assume $\{\mathbf{v}_1, \ldots, \mathbf{v}_{j-1}\}$ is pairwise orthogonal and $\operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_{j-1}\} = \operatorname{span}\{\mathbf{x}_1, \ldots, \mathbf{x}_{j-1}\}$.

    By definition, $\mathbf{v}_j = \mathbf{x}_j - \operatorname{proj}_{W_{j-1}}(\mathbf{x}_j)$, where $W_{j-1} = \operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_{j-1}\}$.

    Verifying orthogonality: for any $\ell < j$,

    $$\langle \mathbf{v}_j, \mathbf{v}_\ell \rangle = \langle \mathbf{x}_j, \mathbf{v}_\ell \rangle - \sum_{i=1}^{j-1}\frac{\langle \mathbf{x}_j, \mathbf{v}_i \rangle}{\|\mathbf{v}_i\|^2}\langle \mathbf{v}_i, \mathbf{v}_\ell \rangle = \langle \mathbf{x}_j, \mathbf{v}_\ell \rangle - \frac{\langle \mathbf{x}_j, \mathbf{v}_\ell \rangle}{\|\mathbf{v}_\ell\|^2}\|\mathbf{v}_\ell\|^2 = 0$$

    The last step uses $\langle \mathbf{v}_i, \mathbf{v}_\ell \rangle = 0$ ($i \neq \ell$).

    $\mathbf{v}_j \neq \mathbf{0}$: if $\mathbf{v}_j = \mathbf{0}$, then $\mathbf{x}_j = \operatorname{proj}_{W_{j-1}}(\mathbf{x}_j) \in W_{j-1} = \operatorname{span}\{\mathbf{x}_1, \ldots, \mathbf{x}_{j-1}\}$, contradicting the linear independence of $\{\mathbf{x}_1, \ldots, \mathbf{x}_k\}$.

    Equal spans: $\mathbf{v}_j \in \operatorname{span}\{\mathbf{x}_1, \ldots, \mathbf{x}_j\}$ (by definition), and $\mathbf{x}_j \in \operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_j\}$ (by solving from the definition), so the two spans are equal. $\blacksquare$

!!! example "Example 7.5"
    Apply Gram-Schmidt orthogonalization to the basis $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$, $\mathbf{x}_2 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}$, $\mathbf{x}_3 = \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}$ of $\mathbb{R}^3$.

    **Step 1**: $\mathbf{v}_1 = \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$.

    **Step 2**:

    $$\frac{\langle \mathbf{x}_2, \mathbf{v}_1 \rangle}{\|\mathbf{v}_1\|^2} = \frac{1 \cdot 1 + 0 \cdot 1 + 1 \cdot 0}{1 + 1 + 0} = \frac{1}{2}$$

    $$\mathbf{v}_2 = \mathbf{x}_2 - \frac{1}{2}\mathbf{v}_1 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix} - \frac{1}{2}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix}$$

    Verification: $\langle \mathbf{v}_1, \mathbf{v}_2 \rangle = 1/2 - 1/2 + 0 = 0$.

    **Step 3**:

    $$\frac{\langle \mathbf{x}_3, \mathbf{v}_1 \rangle}{\|\mathbf{v}_1\|^2} = \frac{0 + 1 + 0}{2} = \frac{1}{2}, \quad \frac{\langle \mathbf{x}_3, \mathbf{v}_2 \rangle}{\|\mathbf{v}_2\|^2} = \frac{0 - 1/2 + 1}{1/4 + 1/4 + 1} = \frac{1/2}{3/2} = \frac{1}{3}$$

    $$\mathbf{v}_3 = \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix} - \frac{1}{2}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} - \frac{1}{3}\begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix} = \begin{pmatrix} -2/3 \\ 2/3 \\ 2/3 \end{pmatrix}$$

    Verification: $\langle \mathbf{v}_1, \mathbf{v}_3 \rangle = -2/3 + 2/3 + 0 = 0$, $\langle \mathbf{v}_2, \mathbf{v}_3 \rangle = -1/3 - 1/3 + 2/3 = 0$.

    Normalization:

    $$\mathbf{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad \mathbf{q}_2 = \frac{1}{\sqrt{3/2}}\begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix} = \sqrt{\frac{2}{3}}\begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix}, \quad \mathbf{q}_3 = \frac{1}{\sqrt{4/3}}\begin{pmatrix} -2/3 \\ 2/3 \\ 2/3 \end{pmatrix}$$

---

## 7.6 Orthogonal Matrices

<div class="context-flow" markdown>

**Matrices that preserve geometric structure**: $Q^TQ=I$ → columns are orthonormal → preserves inner products/norms (isometry) → $\det Q = \pm 1$ (rotation or reflection) → the $Q$ in the spectral theorem $A=QDQ^T$ from Chapter 6

</div>

!!! definition "Definition 7.9 (Orthogonal matrix)"
    A real square matrix $Q$ is called an **orthogonal matrix** if

    $$Q^TQ = QQ^T = I$$

    That is, $Q^{-1} = Q^T$.

!!! theorem "Theorem 7.11 (Equivalent characterizations of orthogonal matrices)"
    Let $Q$ be an $n \times n$ real matrix. The following conditions are equivalent:

    1. $Q$ is an orthogonal matrix.
    2. The columns of $Q$ form an orthonormal basis of $\mathbb{R}^n$.
    3. The rows of $Q$ form an orthonormal basis of $\mathbb{R}^n$.
    4. For any $\mathbf{x} \in \mathbb{R}^n$, $\|Q\mathbf{x}\| = \|\mathbf{x}\|$ (norm-preserving / isometric).
    5. For any $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$, $\langle Q\mathbf{x}, Q\mathbf{y} \rangle = \langle \mathbf{x}, \mathbf{y} \rangle$ (inner-product-preserving).

??? proof "Proof"
    $(1) \Leftrightarrow (2)$: $Q^TQ = I$ means $(Q^TQ)_{ij} = \mathbf{q}_i^T\mathbf{q}_j = \delta_{ij}$ (Kronecker delta), i.e., the columns are orthonormal.

    $(1) \Leftrightarrow (3)$: $QQ^T = I$ means the rows are orthonormal.

    $(1) \Rightarrow (5)$: $\langle Q\mathbf{x}, Q\mathbf{y} \rangle = (Q\mathbf{x})^T(Q\mathbf{y}) = \mathbf{x}^TQ^TQ\mathbf{y} = \mathbf{x}^T\mathbf{y} = \langle \mathbf{x}, \mathbf{y} \rangle$.

    $(5) \Rightarrow (4)$: Setting $\mathbf{y} = \mathbf{x}$, $\|Q\mathbf{x}\|^2 = \langle Q\mathbf{x}, Q\mathbf{x} \rangle = \langle \mathbf{x}, \mathbf{x} \rangle = \|\mathbf{x}\|^2$.

    $(4) \Rightarrow (1)$: From $\|Q\mathbf{x}\|^2 = \|\mathbf{x}\|^2$ for all $\mathbf{x}$, we get $\mathbf{x}^TQ^TQ\mathbf{x} = \mathbf{x}^T\mathbf{x}$ for all $\mathbf{x}$, so $\mathbf{x}^T(Q^TQ - I)\mathbf{x} = 0$ for all $\mathbf{x}$. Since $Q^TQ - I$ is symmetric, this implies $Q^TQ - I = O$, i.e., $Q^TQ = I$. $\blacksquare$

!!! theorem "Theorem 7.12 (Properties of orthogonal matrices)"
    Let $Q$ be an orthogonal matrix. Then:

    1. $\det(Q) = \pm 1$
    2. $Q^{-1} = Q^T$ is also an orthogonal matrix
    3. If $Q_1, Q_2$ are both orthogonal matrices, then $Q_1Q_2$ is also an orthogonal matrix
    4. The eigenvalues of $Q$ have modulus $1$ (i.e., $|\lambda| = 1$)

??? proof "Proof"
    1. $\det(Q^TQ) = \det(I) = 1$. $\det(Q^T)\det(Q) = (\det Q)^2 = 1$, so $\det Q = \pm 1$.

    2. $(Q^T)^T Q^T = QQ^T = I$.

    3. $(Q_1Q_2)^T(Q_1Q_2) = Q_2^TQ_1^TQ_1Q_2 = Q_2^TQ_2 = I$.

    4. Let $Q\mathbf{v} = \lambda\mathbf{v}$ ($\mathbf{v} \neq \mathbf{0}$, $\lambda \in \mathbb{C}$). Taking the conjugate-transpose norm: $\|\mathbf{v}\|^2 = \overline{\mathbf{v}}^T\mathbf{v}$, $\|Q\mathbf{v}\|^2 = \overline{(Q\mathbf{v})}^T(Q\mathbf{v}) = |\lambda|^2\overline{\mathbf{v}}^T\mathbf{v} = |\lambda|^2\|\mathbf{v}\|^2$. By the isometric property $\|Q\mathbf{v}\| = \|\mathbf{v}\|$, we get $|\lambda| = 1$. $\blacksquare$

!!! example "Example 7.6"
    The rotation matrix $Q = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$ is an orthogonal matrix.

    Verification: $Q^TQ = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$.

    $\det Q = \cos^2\theta + \sin^2\theta = 1 > 0$, so this is a **rotation** (an orientation-preserving orthogonal transformation).

    The reflection matrix $H = \begin{pmatrix} \cos 2\alpha & \sin 2\alpha \\ \sin 2\alpha & -\cos 2\alpha \end{pmatrix}$ (reflection about the line at angle $\alpha$) is also an orthogonal matrix, but $\det H = -1$ (orientation-reversing).

---

## 7.7 Least Squares

<div class="context-flow" markdown>

**The optimal substitute when $A\mathbf{x}=\mathbf{b}$ has no solution**: project $\mathbf{b}$ onto $\operatorname{Col}(A)$ → residual $\perp$ column space → normal equations $A^TA\hat{\mathbf{x}}=A^T\mathbf{b}$ → back to Chapter 1: solve a **system of equations**

</div>

When the linear system $A\mathbf{x} = \mathbf{b}$ has no solution (i.e., $\mathbf{b} \notin \operatorname{Col}(A)$), we seek the approximate solution that minimizes the error — this is the core idea of the least squares method.

!!! definition "Definition 7.10 (Least squares solution)"
    Let $A$ be an $m \times n$ matrix and $\mathbf{b} \in \mathbb{R}^m$. A vector $\hat{\mathbf{x}} \in \mathbb{R}^n$ is called a **least squares solution** of $A\mathbf{x} = \mathbf{b}$ if

    $$\|A\hat{\mathbf{x}} - \mathbf{b}\| \leq \|A\mathbf{x} - \mathbf{b}\|, \quad \forall\, \mathbf{x} \in \mathbb{R}^n$$

    That is, $\hat{\mathbf{x}}$ minimizes the residual $\|\mathbf{b} - A\mathbf{x}\|$.

!!! theorem "Theorem 7.13 (Normal equations)"
    $\hat{\mathbf{x}}$ is a least squares solution of $A\mathbf{x} = \mathbf{b}$ if and only if $\hat{\mathbf{x}}$ satisfies the **normal equations**:

    $$A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$$

    If the columns of $A$ are linearly independent, then $A^TA$ is invertible and the least squares solution is unique:

    $$\hat{\mathbf{x}} = (A^TA)^{-1}A^T\mathbf{b}$$

??? proof "Proof"
    **Geometric argument**: $A\hat{\mathbf{x}}$ is the orthogonal projection of $\mathbf{b}$ onto $\operatorname{Col}(A)$. By the best approximation theorem (Theorem 7.9), $\|A\hat{\mathbf{x}} - \mathbf{b}\|$ is minimized if and only if $A\hat{\mathbf{x}} = \operatorname{proj}_{\operatorname{Col}(A)}(\mathbf{b})$, i.e.,

    $$\mathbf{b} - A\hat{\mathbf{x}} \in \operatorname{Col}(A)^\perp = \operatorname{Null}(A^T)$$

    That is, $A^T(\mathbf{b} - A\hat{\mathbf{x}}) = \mathbf{0}$, which simplifies to $A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$.

    **Uniqueness**: If the columns of $A$ are linearly independent, then $\operatorname{Null}(A) = \{\mathbf{0}\}$. From $A^TA\mathbf{x} = \mathbf{0}$ $\Rightarrow$ $\mathbf{x}^TA^TA\mathbf{x} = 0$ $\Rightarrow$ $\|A\mathbf{x}\|^2 = 0$ $\Rightarrow$ $A\mathbf{x} = \mathbf{0}$ $\Rightarrow$ $\mathbf{x} = \mathbf{0}$, so $A^TA$ is invertible. $\blacksquare$

!!! example "Example 7.7"
    Use least squares to fit the line $y = c_0 + c_1 x$ to the data points $(1, 1), (2, 3), (3, 2), (4, 5)$.

    Setting up the matrix equation $A\mathbf{c} = \mathbf{b}$:

    $$A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \\ 1 & 4 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 1 \\ 3 \\ 2 \\ 5 \end{pmatrix}, \quad \mathbf{c} = \begin{pmatrix} c_0 \\ c_1 \end{pmatrix}$$

    Normal equations $A^TA\mathbf{c} = A^T\mathbf{b}$:

    $$A^TA = \begin{pmatrix} 4 & 10 \\ 10 & 30 \end{pmatrix}, \quad A^T\mathbf{b} = \begin{pmatrix} 11 \\ 33 \end{pmatrix}$$

    Solving: $\begin{pmatrix} 4 & 10 \\ 10 & 30 \end{pmatrix}\begin{pmatrix} c_0 \\ c_1 \end{pmatrix} = \begin{pmatrix} 11 \\ 33 \end{pmatrix}$.

    $\det(A^TA) = 120 - 100 = 20$.

    $$(A^TA)^{-1} = \frac{1}{20}\begin{pmatrix} 30 & -10 \\ -10 & 4 \end{pmatrix}$$

    $$\hat{\mathbf{c}} = \frac{1}{20}\begin{pmatrix} 30 & -10 \\ -10 & 4 \end{pmatrix}\begin{pmatrix} 11 \\ 33 \end{pmatrix} = \frac{1}{20}\begin{pmatrix} 330 - 330 \\ -110 + 132 \end{pmatrix} = \frac{1}{20}\begin{pmatrix} 0 \\ 22 \end{pmatrix} = \begin{pmatrix} 0 \\ 1.1 \end{pmatrix}$$

    The least squares fit line is $y = 1.1x$.

!!! example "Example 7.8"
    **Geometric interpretation of least squares**: When $A\mathbf{x} = \mathbf{b}$ has no solution, $\mathbf{b}$ does not lie in $\operatorname{Col}(A)$. The least squares solution gives $A\hat{\mathbf{x}} = \operatorname{proj}_{\operatorname{Col}(A)}(\mathbf{b})$, i.e., it "projects" $\mathbf{b}$ onto the column space of $A$. The residual vector $\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}}$ is perpendicular to $\operatorname{Col}(A)$, which is precisely the meaning of the normal equations $A^T\mathbf{r} = \mathbf{0}$.

    Let $A = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$, $\mathbf{b} = \begin{pmatrix} 1 \\ 2 \\ 4 \end{pmatrix}$.

    $A^TA = 3$, $A^T\mathbf{b} = 7$, $\hat{x} = 7/3$.

    $A\hat{x} = \begin{pmatrix} 7/3 \\ 7/3 \\ 7/3 \end{pmatrix}$, which is the average of $1, 2, 4$, namely $7/3$ — the least squares solution for a "constant fit."

---

## 7.8 QR Decomposition

<div class="context-flow" markdown>

**Gram-Schmidt in matrix form**: $A = QR$ (orthogonal $\times$ upper triangular) → least squares solution simplifies to $R\hat{\mathbf{x}}=Q^T\mathbf{b}$ (back substitution in upper triangular system) → numerically more stable than the normal equations

</div>

!!! definition "Definition 7.11 (QR decomposition)"
    Let $A$ be an $m \times n$ matrix ($m \geq n$) with linearly independent columns. The **QR decomposition** (QR factorization) of $A$ writes $A$ as

    $$A = QR$$

    where $Q$ is an $m \times n$ matrix whose columns form an orthonormal set ($Q^TQ = I_n$), and $R$ is an $n \times n$ upper triangular matrix with positive diagonal entries.

!!! theorem "Theorem 7.14 (Existence and uniqueness of QR decomposition)"
    Let $A$ be an $m \times n$ matrix ($m \geq n$) with linearly independent columns. Then $A$ has a unique QR decomposition $A = QR$, where the columns of $Q$ are orthonormal and $R$ is upper triangular with positive diagonal entries.

??? proof "Proof"
    **Existence**: Apply Gram-Schmidt orthogonalization to the columns $\{\mathbf{a}_1, \ldots, \mathbf{a}_n\}$ of $A$ to obtain the orthonormal set $\{\mathbf{q}_1, \ldots, \mathbf{q}_n\}$.

    By the properties of the Gram-Schmidt process, for each $j$:

    $$\operatorname{span}\{\mathbf{q}_1, \ldots, \mathbf{q}_j\} = \operatorname{span}\{\mathbf{a}_1, \ldots, \mathbf{a}_j\}$$

    Therefore $\mathbf{a}_j \in \operatorname{span}\{\mathbf{q}_1, \ldots, \mathbf{q}_j\}$ and can be written as

    $$\mathbf{a}_j = r_{1j}\mathbf{q}_1 + r_{2j}\mathbf{q}_2 + \cdots + r_{jj}\mathbf{q}_j$$

    where $r_{ij} = \langle \mathbf{a}_j, \mathbf{q}_i \rangle$ ($i \leq j$) and $r_{ij} = 0$ ($i > j$).

    Setting $Q = (\mathbf{q}_1 \;\; \cdots \;\; \mathbf{q}_n)$ and $R = (r_{ij})$, we get $A = QR$ with $R$ upper triangular.

    The diagonal entries satisfy $r_{jj} = \langle \mathbf{a}_j, \mathbf{q}_j \rangle$. From the Gram-Schmidt construction, $r_{jj} = \|\mathbf{v}_j\| > 0$ (where $\mathbf{v}_j$ is the vector after orthogonalization but before normalization).

    **Uniqueness**: Suppose $A = Q_1R_1 = Q_2R_2$ are two QR decompositions. Then $Q_2^TQ_1 = R_2R_1^{-1}$. The left side has orthonormal columns ($Q_2^TQ_1$ is an orthogonal matrix), and the right side is upper triangular. A matrix that is both orthogonal and upper triangular must be diagonal with diagonal entries of modulus $1$. Since $R_1, R_2$ have positive diagonal entries, $R_2R_1^{-1}$ also has positive diagonal entries, so the diagonal entries are $1$, giving $Q_1 = Q_2$ and $R_1 = R_2$. $\blacksquare$

!!! theorem "Theorem 7.15 (Least squares via QR decomposition)"
    If $A = QR$ (QR decomposition), then the least squares solution of $A\mathbf{x} = \mathbf{b}$ is

    $$\hat{\mathbf{x}} = R^{-1}Q^T\mathbf{b}$$

    That is, one only needs to solve the upper triangular system $R\hat{\mathbf{x}} = Q^T\mathbf{b}$.

??? proof "Proof"
    From the normal equations $A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$, substituting $A = QR$:

    $$(QR)^T(QR)\hat{\mathbf{x}} = (QR)^T\mathbf{b}$$
    $$R^TQ^TQR\hat{\mathbf{x}} = R^TQ^T\mathbf{b}$$
    $$R^TR\hat{\mathbf{x}} = R^TQ^T\mathbf{b}$$

    Since $R$ is invertible (positive diagonal entries), so is $R^T$. Multiplying on the left by $(R^T)^{-1}$:

    $$R\hat{\mathbf{x}} = Q^T\mathbf{b}$$

    This is an upper triangular system that can be efficiently solved by back substitution. $\blacksquare$

!!! note "Note"
    QR decomposition is numerically more stable than the normal equations. The normal equations require computing $A^TA$, which squares the condition number and amplifies rounding errors. QR decomposition avoids this issue. Therefore, in practical numerical computation, QR decomposition is the preferred method for solving least squares problems.

!!! example "Example 7.9"
    Find the QR decomposition of $A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \\ 0 & 1 \end{pmatrix}$.

    $\mathbf{a}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$, $\mathbf{a}_2 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}$.

    **Gram-Schmidt**:

    $\mathbf{v}_1 = \mathbf{a}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$, $\|\mathbf{v}_1\| = \sqrt{2}$.

    $\mathbf{v}_2 = \mathbf{a}_2 - \frac{\langle \mathbf{a}_2, \mathbf{v}_1 \rangle}{\|\mathbf{v}_1\|^2}\mathbf{v}_1 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix} - \frac{1}{2}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix}$, $\|\mathbf{v}_2\| = \sqrt{1/4 + 1/4 + 1} = \sqrt{3/2}$.

    **Normalization**:

    $$\mathbf{q}_1 = \frac{\mathbf{v}_1}{\|\mathbf{v}_1\|} = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad \mathbf{q}_2 = \frac{\mathbf{v}_2}{\|\mathbf{v}_2\|} = \sqrt{\frac{2}{3}}\begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix} = \frac{1}{\sqrt{6}}\begin{pmatrix} 1 \\ -1 \\ 2 \end{pmatrix}$$

    **Computing $R$**:

    $$r_{11} = \langle \mathbf{a}_1, \mathbf{q}_1 \rangle = \frac{1}{\sqrt{2}}(1 + 1) = \sqrt{2}$$

    $$r_{12} = \langle \mathbf{a}_2, \mathbf{q}_1 \rangle = \frac{1}{\sqrt{2}}(1 + 0) = \frac{1}{\sqrt{2}}$$

    $$r_{22} = \langle \mathbf{a}_2, \mathbf{q}_2 \rangle = \frac{1}{\sqrt{6}}(1 + 0 + 2) = \frac{3}{\sqrt{6}} = \sqrt{\frac{3}{2}}$$

    Therefore:

    $$Q = \begin{pmatrix} \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{6}} \\ \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{6}} \\ 0 & \frac{2}{\sqrt{6}} \end{pmatrix}, \quad R = \begin{pmatrix} \sqrt{2} & \frac{1}{\sqrt{2}} \\ 0 & \sqrt{\frac{3}{2}} \end{pmatrix}$$

    Verifying $A = QR$:

    $$QR = \begin{pmatrix} \frac{1}{\sqrt{2}}\sqrt{2} + \frac{1}{\sqrt{6}} \cdot 0 & \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}} + \frac{1}{\sqrt{6}}\sqrt{\frac{3}{2}} \\ \frac{1}{\sqrt{2}}\sqrt{2} - \frac{1}{\sqrt{6}} \cdot 0 & \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}} - \frac{1}{\sqrt{6}}\sqrt{\frac{3}{2}} \\ 0 + \frac{2}{\sqrt{6}} \cdot 0 & 0 + \frac{2}{\sqrt{6}}\sqrt{\frac{3}{2}} \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 0 \\ 0 & 1 \end{pmatrix}$$

    Verified.

!!! example "Example 7.10"
    Use QR decomposition to solve the least squares problem from Example 7.7.

    $A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \\ 1 & 4 \end{pmatrix}$. Applying Gram-Schmidt to the columns:

    $\mathbf{a}_1 = \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}$, $\|\mathbf{a}_1\| = 2$, $\mathbf{q}_1 = \frac{1}{2}\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}$.

    $\mathbf{v}_2 = \begin{pmatrix} 1 \\ 2 \\ 3 \\ 4 \end{pmatrix} - \frac{\langle \mathbf{a}_2, \mathbf{q}_1 \rangle}{\|\mathbf{q}_1\|^2}\mathbf{q}_1 \cdot \|\mathbf{q}_1\|$. First compute $\langle \mathbf{a}_2, \mathbf{a}_1 \rangle / \|\mathbf{a}_1\|^2 = 10/4 = 5/2$.

    $\mathbf{v}_2 = \begin{pmatrix} 1 \\ 2 \\ 3 \\ 4 \end{pmatrix} - \frac{5}{2}\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} -3/2 \\ -1/2 \\ 1/2 \\ 3/2 \end{pmatrix}$, $\|\mathbf{v}_2\| = \sqrt{9/4 + 1/4 + 1/4 + 9/4} = \sqrt{5}$.

    $\mathbf{q}_2 = \frac{1}{\sqrt{5}}\begin{pmatrix} -3/2 \\ -1/2 \\ 1/2 \\ 3/2 \end{pmatrix}$.

    $R = \begin{pmatrix} 2 & 5/2 \cdot 2/1 \\ 0 & \sqrt{5} \end{pmatrix}$. More precisely: $r_{11} = \|\mathbf{a}_1\| = 2$, $r_{12} = \langle \mathbf{a}_2, \mathbf{q}_1 \rangle = 10/2 = 5$, $r_{22} = \|\mathbf{v}_2\| = \sqrt{5}$.

    $R = \begin{pmatrix} 2 & 5 \\ 0 & \sqrt{5} \end{pmatrix}$.

    $Q^T\mathbf{b} = \begin{pmatrix} \mathbf{q}_1^T\mathbf{b} \\ \mathbf{q}_2^T\mathbf{b} \end{pmatrix} = \begin{pmatrix} (1+3+2+5)/2 \\ (-3/2-1/2\cdot 3+1/2\cdot 2+3/2\cdot 5)/\sqrt{5} \end{pmatrix} = \begin{pmatrix} 11/2 \\ (-3/2-3/2+1+15/2)/\sqrt{5} \end{pmatrix}$

    Computing $\mathbf{q}_2^T\mathbf{b}$: $\frac{1}{\sqrt{5}}(-3/2 \cdot 1 + (-1/2) \cdot 3 + 1/2 \cdot 2 + 3/2 \cdot 5) = \frac{1}{\sqrt{5}}(-3/2 - 3/2 + 1 + 15/2) = \frac{1}{\sqrt{5}} \cdot \frac{-3-3+2+15}{2} = \frac{11}{2\sqrt{5}}$.

    Solving $R\hat{\mathbf{c}} = Q^T\mathbf{b}$:

    $$\begin{pmatrix} 2 & 5 \\ 0 & \sqrt{5} \end{pmatrix}\begin{pmatrix} c_0 \\ c_1 \end{pmatrix} = \begin{pmatrix} 11/2 \\ 11/(2\sqrt{5}) \end{pmatrix}$$

    From the second row: $\sqrt{5}\, c_1 = \frac{11}{2\sqrt{5}}$, so $c_1 = \frac{11}{10} = 1.1$.

    From the first row: $2c_0 + 5 \cdot 1.1 = 5.5$, so $2c_0 = 0$, $c_0 = 0$.

    Result: $\hat{\mathbf{c}} = \begin{pmatrix} 0 \\ 1.1 \end{pmatrix}$, consistent with Example 7.7.

---

## Chapter Summary

This chapter systematically studied the theory of orthogonality in $\mathbb{R}^n$ and its applications. The main topics include:

1. **Inner products and norms**: the inner product endows a vector space with geometric structure; the Cauchy-Schwarz inequality and triangle inequality are the most fundamental inequalities.
2. **Orthogonal sets**: orthogonal sets of nonzero vectors are automatically linearly independent, and orthogonal bases make coordinate computation extremely convenient.
3. **Orthogonal complements**: $\mathbb{R}^n = W \oplus W^\perp$ guarantees the existence and uniqueness of the orthogonal decomposition.
4. **Orthogonal projection**: the projection formula $P = A(A^TA)^{-1}A^T$ gives the best approximation.
5. **Gram-Schmidt orthogonalization**: an algorithm that converts any basis into an orthogonal basis (or orthonormal basis).
6. **Orthogonal matrices**: matrices satisfying $Q^TQ = I$ preserve inner products and norms, representing rigid transformations (rotations or reflections).
7. **Least squares**: the normal equations $A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$ provide the optimal approximate solution to an overdetermined system.
8. **QR decomposition**: $A = QR$ decomposes a matrix into an orthogonal part and an upper triangular part, forming the basis of numerically stable least squares algorithms.
