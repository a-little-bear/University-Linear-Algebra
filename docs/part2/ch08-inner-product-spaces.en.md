# Chapter 8  Inner Product Spaces

<div class="context-flow" markdown>

**Prerequisites**: Ch7 Orthogonality in $\mathbb{R}^n$

**Chapter arc**: Inner product axioms → Norm/Cauchy-Schwarz → Orthogonal complement/Projection → Adjoint operators → Normal/Self-adjoint operators → **Spectral theorem**

**Further connections**：Inner product spaces generalize to Hilbert spaces ($L^2$ spaces, Sobolev spaces) in infinite dimensions; quantum mechanical state spaces are Hilbert spaces; Fourier analysis is essentially orthogonal expansion in $L^2$; reproducing kernel Hilbert spaces (RKHS) underpin kernel methods in machine learning

</div>

An inner product space is one of the most geometrically intuitive structures in linear algebra. By introducing an inner product operation on a vector space, we can rigorously define geometric concepts such as length, angle, and orthogonality, thereby deeply fusing algebraic structure with geometric intuition. Starting from the axiomatic definition of inner products, this chapter systematically establishes the theory of norms and distances, studies orthogonality and orthogonal decompositions, delves into the properties of adjoint, normal, and self-adjoint operators, and ultimately derives one of the most profound results in linear algebra — the spectral theorem.

---

## 8.1 Definition of Inner Products

<div class="context-flow" markdown>

**Axiomatize** the dot product $\mathbf{x}^T\mathbf{y}$ in $\mathbb{R}^n$: positive definiteness + symmetry (or conjugate symmetry) + linearity → applicable to function spaces and other general spaces

</div>

In Euclidean space $\mathbb{R}^n$, the dot product $\mathbf{x} \cdot \mathbf{y} = \sum_{i=1}^n x_i y_i$ endows vectors with the concepts of length and angle. The inner product is a generalization of the dot product, introducing this geometric structure into general vector spaces.

### Real Inner Product Spaces

!!! definition "Definition 8.1 (Real inner product)"
    Let $V$ be a vector space over $\mathbb{R}$. A **real inner product** on $V$ is a map $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{R}$ satisfying the following four axioms: for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and $\alpha \in \mathbb{R}$,

    1. **Positive definiteness**: $\langle \mathbf{v}, \mathbf{v} \rangle \geq 0$, and $\langle \mathbf{v}, \mathbf{v} \rangle = 0$ if and only if $\mathbf{v} = \mathbf{0}$;
    2. **Symmetry**: $\langle \mathbf{u}, \mathbf{v} \rangle = \langle \mathbf{v}, \mathbf{u} \rangle$;
    3. **Homogeneity in the first argument**: $\langle \alpha\mathbf{u}, \mathbf{v} \rangle = \alpha \langle \mathbf{u}, \mathbf{v} \rangle$;
    4. **Additivity in the first argument**: $\langle \mathbf{u} + \mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{u}, \mathbf{w} \rangle + \langle \mathbf{v}, \mathbf{w} \rangle$.

    A real vector space equipped with a real inner product $(V, \langle \cdot, \cdot \rangle)$ is called a **real inner product space**.

!!! note "Note"
    Axioms 3 and 4 together are equivalent to the inner product being linear in the first argument: $\langle \alpha\mathbf{u} + \beta\mathbf{v}, \mathbf{w} \rangle = \alpha\langle \mathbf{u}, \mathbf{w} \rangle + \beta\langle \mathbf{v}, \mathbf{w} \rangle$. By symmetry, the real inner product is also linear in the second argument, i.e., the real inner product is **bilinear**.

### Complex Inner Product Spaces

!!! definition "Definition 8.2 (Complex inner product)"
    Let $V$ be a vector space over $\mathbb{C}$. A **complex inner product** (or Hermitian inner product) on $V$ is a map $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{C}$ satisfying the following four axioms: for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and $\alpha \in \mathbb{C}$,

    1. **Positive definiteness**: $\langle \mathbf{v}, \mathbf{v} \rangle \geq 0$ (i.e., $\langle \mathbf{v}, \mathbf{v} \rangle \in \mathbb{R}$ and is nonnegative), and $\langle \mathbf{v}, \mathbf{v} \rangle = 0$ if and only if $\mathbf{v} = \mathbf{0}$;
    2. **Conjugate symmetry**: $\langle \mathbf{u}, \mathbf{v} \rangle = \overline{\langle \mathbf{v}, \mathbf{u} \rangle}$;
    3. **Homogeneity in the first argument**: $\langle \alpha\mathbf{u}, \mathbf{v} \rangle = \alpha \langle \mathbf{u}, \mathbf{v} \rangle$;
    4. **Additivity in the first argument**: $\langle \mathbf{u} + \mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{u}, \mathbf{w} \rangle + \langle \mathbf{v}, \mathbf{w} \rangle$.

    A complex vector space equipped with a complex inner product $(V, \langle \cdot, \cdot \rangle)$ is called a **complex inner product space**.

!!! note "Note"
    In a complex inner product space, the inner product is linear in the first argument but **conjugate linear** in the second argument: $\langle \mathbf{u}, \alpha\mathbf{v} \rangle = \bar{\alpha}\langle \mathbf{u}, \mathbf{v} \rangle$. This property is called **sesquilinearity**. Note that some references (e.g., physics literature) adopt the convention of linearity in the second argument; here we follow the more common mathematical convention of linearity in the first argument.

!!! definition "Definition 8.3 (Inner product space)"
    An **inner product space** is a vector space equipped with an inner product $(V, \langle \cdot, \cdot \rangle)$. When the base field is $\mathbb{R}$ it is called a real inner product space, and when the base field is $\mathbb{C}$ it is called a complex inner product space. A finite-dimensional real inner product space is called a **Euclidean space**, and a finite-dimensional complex inner product space is called a **unitary space**.

!!! example "Example 8.1"
    **Standard inner product on $\mathbb{R}^n$.** For $\mathbf{x} = (x_1, \ldots, x_n)^T, \mathbf{y} = (y_1, \ldots, y_n)^T \in \mathbb{R}^n$, define

    $$\langle \mathbf{x}, \mathbf{y} \rangle = \mathbf{x}^T\mathbf{y} = \sum_{i=1}^n x_i y_i$$

    This is the most common inner product, namely the standard dot product.

!!! example "Example 8.2"
    **Standard inner product on $\mathbb{C}^n$.** For $\mathbf{x} = (x_1, \ldots, x_n)^T, \mathbf{y} = (y_1, \ldots, y_n)^T \in \mathbb{C}^n$, define

    $$\langle \mathbf{x}, \mathbf{y} \rangle = \mathbf{x}^H\mathbf{y} = \sum_{i=1}^n x_i \overline{y_i}$$

    where $\mathbf{x}^H = \overline{\mathbf{x}}^T$ is the conjugate transpose.

!!! example "Example 8.3"
    **Inner product on the space of continuous functions.** On the real vector space $C[a, b]$ (continuous real-valued functions on the interval $[a, b]$), define

    $$\langle f, g \rangle = \int_a^b f(t)g(t) \, dt$$

    One can verify that this satisfies the four axioms of a real inner product. In particular, positive definiteness is guaranteed by the fact that if $f$ is continuous and $\int_a^b f(t)^2 \, dt = 0$, then $f \equiv 0$.

!!! example "Example 8.4"
    **Weighted inner product.** Let $w_1, w_2, \ldots, w_n > 0$ be positive real numbers (called weights). For $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$, define

    $$\langle \mathbf{x}, \mathbf{y} \rangle_w = \sum_{i=1}^n w_i x_i y_i$$

    This is an inner product on $\mathbb{R}^n$, called the weighted inner product.

!!! proposition "Proposition 8.1"
    Let $(V, \langle \cdot, \cdot \rangle)$ be an inner product space. For any $\mathbf{u}, \mathbf{v} \in V$,

    $$\langle \mathbf{u}, \mathbf{0} \rangle = \langle \mathbf{0}, \mathbf{v} \rangle = 0$$

??? proof "Proof"
    By linearity, $\langle \mathbf{u}, \mathbf{0} \rangle = \langle \mathbf{u}, 0 \cdot \mathbf{v} \rangle = \bar{0} \langle \mathbf{u}, \mathbf{v} \rangle = 0$. Similarly, $\langle \mathbf{0}, \mathbf{v} \rangle = \langle 0 \cdot \mathbf{u}, \mathbf{v} \rangle = 0 \cdot \langle \mathbf{u}, \mathbf{v} \rangle = 0$. $\blacksquare$

---

## 8.2 Norms and Distances

<div class="context-flow" markdown>

Inner product → Norm $\|\mathbf{v}\| = \sqrt{\langle \mathbf{v},\mathbf{v}\rangle}$ → **Cauchy-Schwarz** (the most fundamental inequality in inner product spaces) → Triangle inequality → Distance

Key insight: The parallelogram identity characterizes "which norms come from an inner product"

</div>

An inner product naturally induces a notion of "length" (norm) and "distance" between vectors.

!!! definition "Definition 8.4 (Induced norm)"
    Let $(V, \langle \cdot, \cdot \rangle)$ be an inner product space. The **induced norm** is defined as

    $$\|\mathbf{v}\| = \sqrt{\langle \mathbf{v}, \mathbf{v} \rangle}, \quad \forall \mathbf{v} \in V$$

    The **distance** between vectors $\mathbf{v}$ and $\mathbf{w}$ is defined as $d(\mathbf{v}, \mathbf{w}) = \|\mathbf{v} - \mathbf{w}\|$.

!!! definition "Definition 8.5 (Unit vector)"
    If $\|\mathbf{v}\| = 1$, then $\mathbf{v}$ is called a **unit vector**. For any nonzero vector $\mathbf{v}$, the vector $\hat{\mathbf{v}} = \frac{\mathbf{v}}{\|\mathbf{v}\|}$ is a unit vector, and the process of transforming $\mathbf{v}$ into $\hat{\mathbf{v}}$ is called **normalization**.

### Cauchy-Schwarz Inequality

<div class="context-flow" markdown>

**Insight**: The core of the Cauchy-Schwarz proof is choosing the optimal $t$ in $\|\mathbf{u} - t\mathbf{v}\|^2 \ge 0$ — essentially the optimality of **projection**

</div>

!!! theorem "Theorem 8.1 (Cauchy-Schwarz inequality)"
    Let $(V, \langle \cdot, \cdot \rangle)$ be an inner product space. For any $\mathbf{u}, \mathbf{v} \in V$,

    $$|\langle \mathbf{u}, \mathbf{v} \rangle| \leq \|\mathbf{u}\| \cdot \|\mathbf{v}\|$$

    Equality holds if and only if $\mathbf{u}$ and $\mathbf{v}$ are linearly dependent.

??? proof "Proof"
    If $\mathbf{v} = \mathbf{0}$, both sides equal $0$ and the inequality holds trivially.

    Assume $\mathbf{v} \neq \mathbf{0}$. For any scalar $t$, by positive definiteness,

    $$0 \leq \|\mathbf{u} - t\mathbf{v}\|^2 = \langle \mathbf{u} - t\mathbf{v}, \mathbf{u} - t\mathbf{v} \rangle = \|\mathbf{u}\|^2 - t\langle \mathbf{v}, \mathbf{u}\rangle - \bar{t}\langle \mathbf{u}, \mathbf{v}\rangle + |t|^2\|\mathbf{v}\|^2$$

    Taking $t = \dfrac{\langle \mathbf{u}, \mathbf{v}\rangle}{\|\mathbf{v}\|^2}$ and substituting gives

    $$0 \leq \|\mathbf{u}\|^2 - \frac{|\langle \mathbf{u}, \mathbf{v}\rangle|^2}{\|\mathbf{v}\|^2}$$

    Rearranging yields $|\langle \mathbf{u}, \mathbf{v}\rangle|^2 \leq \|\mathbf{u}\|^2 \|\mathbf{v}\|^2$, and taking square roots gives the desired inequality.

    Equality holds if and only if $\|\mathbf{u} - t\mathbf{v}\|^2 = 0$, i.e., $\mathbf{u} = t\mathbf{v}$, in which case $\mathbf{u}, \mathbf{v}$ are linearly dependent. $\blacksquare$

!!! theorem "Theorem 8.2 (Triangle inequality)"
    Let $(V, \langle \cdot, \cdot \rangle)$ be an inner product space. For any $\mathbf{u}, \mathbf{v} \in V$,

    $$\|\mathbf{u} + \mathbf{v}\| \leq \|\mathbf{u}\| + \|\mathbf{v}\|$$

??? proof "Proof"

    $$\|\mathbf{u} + \mathbf{v}\|^2 = \langle \mathbf{u}+\mathbf{v}, \mathbf{u}+\mathbf{v}\rangle = \|\mathbf{u}\|^2 + \langle \mathbf{u}, \mathbf{v}\rangle + \langle \mathbf{v}, \mathbf{u}\rangle + \|\mathbf{v}\|^2$$

    $$= \|\mathbf{u}\|^2 + 2\operatorname{Re}\langle \mathbf{u}, \mathbf{v}\rangle + \|\mathbf{v}\|^2 \leq \|\mathbf{u}\|^2 + 2|\langle \mathbf{u}, \mathbf{v}\rangle| + \|\mathbf{v}\|^2$$

    By the Cauchy-Schwarz inequality, $|\langle \mathbf{u}, \mathbf{v}\rangle| \leq \|\mathbf{u}\|\|\mathbf{v}\|$, so

    $$\|\mathbf{u}+\mathbf{v}\|^2 \leq \|\mathbf{u}\|^2 + 2\|\mathbf{u}\|\|\mathbf{v}\| + \|\mathbf{v}\|^2 = (\|\mathbf{u}\| + \|\mathbf{v}\|)^2$$

    Taking square roots gives the result. $\blacksquare$

!!! theorem "Theorem 8.3 (Parallelogram identity)"
    Let $(V, \langle \cdot, \cdot \rangle)$ be an inner product space. For any $\mathbf{u}, \mathbf{v} \in V$,

    $$\|\mathbf{u} + \mathbf{v}\|^2 + \|\mathbf{u} - \mathbf{v}\|^2 = 2(\|\mathbf{u}\|^2 + \|\mathbf{v}\|^2)$$

??? proof "Proof"
    Direct computation:

    $$\|\mathbf{u}+\mathbf{v}\|^2 = \|\mathbf{u}\|^2 + 2\operatorname{Re}\langle \mathbf{u}, \mathbf{v}\rangle + \|\mathbf{v}\|^2$$

    $$\|\mathbf{u}-\mathbf{v}\|^2 = \|\mathbf{u}\|^2 - 2\operatorname{Re}\langle \mathbf{u}, \mathbf{v}\rangle + \|\mathbf{v}\|^2$$

    Adding the two, the middle terms cancel, giving $\|\mathbf{u}+\mathbf{v}\|^2 + \|\mathbf{u}-\mathbf{v}\|^2 = 2\|\mathbf{u}\|^2 + 2\|\mathbf{v}\|^2$. $\blacksquare$

!!! note "Note"
    The geometric meaning of the parallelogram identity is: the sum of the squares of the two diagonals of a parallelogram equals the sum of the squares of all four sides. Conversely, a norm satisfying the parallelogram identity is a necessary and sufficient condition for it to be induced by some inner product (Jordan-von Neumann theorem).

!!! theorem "Theorem 8.4 (Polarization identity)"
    In a real inner product space, the inner product can be recovered from the norm:

    $$\langle \mathbf{u}, \mathbf{v} \rangle = \frac{1}{4}\left(\|\mathbf{u}+\mathbf{v}\|^2 - \|\mathbf{u}-\mathbf{v}\|^2\right)$$

    In a complex inner product space:

    $$\langle \mathbf{u}, \mathbf{v} \rangle = \frac{1}{4}\sum_{k=0}^{3} i^k \|\mathbf{u} + i^k\mathbf{v}\|^2 = \frac{1}{4}\left(\|\mathbf{u}+\mathbf{v}\|^2 - \|\mathbf{u}-\mathbf{v}\|^2 + i\|\mathbf{u}+i\mathbf{v}\|^2 - i\|\mathbf{u}-i\mathbf{v}\|^2\right)$$

??? proof "Proof"
    For the real case: expand $\|\mathbf{u}+\mathbf{v}\|^2 - \|\mathbf{u}-\mathbf{v}\|^2$:

    $$\|\mathbf{u}+\mathbf{v}\|^2 - \|\mathbf{u}-\mathbf{v}\|^2 = (\|\mathbf{u}\|^2 + 2\langle \mathbf{u},\mathbf{v}\rangle + \|\mathbf{v}\|^2) - (\|\mathbf{u}\|^2 - 2\langle \mathbf{u},\mathbf{v}\rangle + \|\mathbf{v}\|^2) = 4\langle \mathbf{u},\mathbf{v}\rangle$$

    The complex case is similar: the weighted sum of the four expansions combines the real and imaginary parts into exactly $4\langle \mathbf{u},\mathbf{v}\rangle$. $\blacksquare$

!!! example "Example 8.5"
    In $\mathbb{R}^3$ with the standard inner product, let $\mathbf{u} = (1, 2, 3)^T$, $\mathbf{v} = (4, -1, 2)^T$. Then

    $$\langle \mathbf{u}, \mathbf{v} \rangle = 1 \cdot 4 + 2 \cdot (-1) + 3 \cdot 2 = 8$$

    $$\|\mathbf{u}\| = \sqrt{1+4+9} = \sqrt{14}, \quad \|\mathbf{v}\| = \sqrt{16+1+4} = \sqrt{21}$$

    Verifying the Cauchy-Schwarz inequality: $|8| = 8 \leq \sqrt{14} \cdot \sqrt{21} = \sqrt{294} \approx 17.15$, which holds.

    The angle $\theta$ between $\mathbf{u}$ and $\mathbf{v}$ satisfies $\cos\theta = \dfrac{\langle \mathbf{u}, \mathbf{v}\rangle}{\|\mathbf{u}\|\|\mathbf{v}\|} = \dfrac{8}{\sqrt{294}}$.

---

## 8.3 Orthogonality and Orthogonal Complements

<div class="context-flow" markdown>

$\langle \mathbf{u},\mathbf{v}\rangle = 0$ defines orthogonality → Orthogonal sets are automatically linearly independent → **Orthogonal complement** $W^\perp$ gives $V = W \oplus W^\perp$ (direct sum decomposition)

</div>

Orthogonality is one of the most important concepts in inner product spaces, generalizing the geometric concept of "perpendicularity."

!!! definition "Definition 8.6 (Orthogonality)"
    Let $(V, \langle \cdot, \cdot \rangle)$ be an inner product space.

    - Two vectors $\mathbf{u}, \mathbf{v} \in V$ are called **orthogonal** if $\langle \mathbf{u}, \mathbf{v} \rangle = 0$, denoted $\mathbf{u} \perp \mathbf{v}$.
    - A vector $\mathbf{v}$ is orthogonal to a set $S \subseteq V$ if $\mathbf{v}$ is orthogonal to every vector in $S$.
    - A set $S$ is called an **orthogonal set** if any two distinct vectors in $S$ are orthogonal.
    - A set $S$ is called an **orthonormal set** if $S$ is an orthogonal set and every vector is a unit vector.

!!! definition "Definition 8.7 (Orthogonal complement)"
    Let $W$ be a subspace of the inner product space $V$. The **orthogonal complement** of $W$ is defined as

    $$W^\perp = \{\mathbf{v} \in V : \langle \mathbf{v}, \mathbf{w} \rangle = 0, \, \forall \mathbf{w} \in W\}$$

!!! theorem "Theorem 8.5 (Properties of the orthogonal complement)"
    Let $W$ be a subspace of a finite-dimensional inner product space $V$. Then

    1. $W^\perp$ is a subspace of $V$;
    2. $W \cap W^\perp = \{\mathbf{0}\}$;
    3. $(W^\perp)^\perp = W$;
    4. $\dim W + \dim W^\perp = \dim V$.

??? proof "Proof"
    (1) If $\mathbf{u}, \mathbf{v} \in W^\perp$ and $\alpha, \beta$ are scalars, then for any $\mathbf{w} \in W$:

    $$\langle \alpha\mathbf{u} + \beta\mathbf{v}, \mathbf{w}\rangle = \alpha\langle \mathbf{u}, \mathbf{w}\rangle + \beta\langle \mathbf{v}, \mathbf{w}\rangle = 0$$

    So $\alpha\mathbf{u}+\beta\mathbf{v} \in W^\perp$, hence $W^\perp$ is a subspace.

    (2) If $\mathbf{v} \in W \cap W^\perp$, then $\mathbf{v} \in W$ and $\langle \mathbf{v}, \mathbf{w}\rangle = 0$ for all $\mathbf{w} \in W$. In particular, taking $\mathbf{w} = \mathbf{v}$ gives $\langle \mathbf{v}, \mathbf{v}\rangle = 0$, so by positive definiteness $\mathbf{v} = \mathbf{0}$.

    (3) and (4) will be proved after the orthogonal projection theorem. $\blacksquare$

!!! theorem "Theorem 8.6 (Orthogonal direct sum decomposition)"
    Let $V$ be a finite-dimensional inner product space and $W$ a subspace of $V$. Then

    $$V = W \oplus W^\perp$$

    That is, any $\mathbf{v} \in V$ can be uniquely written as $\mathbf{v} = \mathbf{w} + \mathbf{w}'$, where $\mathbf{w} \in W$ and $\mathbf{w}' \in W^\perp$.

??? proof "Proof"
    Let $\{\mathbf{e}_1, \ldots, \mathbf{e}_k\}$ be an orthonormal basis of $W$. For any $\mathbf{v} \in V$, let

    $$\mathbf{w} = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i \rangle \mathbf{e}_i, \quad \mathbf{w}' = \mathbf{v} - \mathbf{w}$$

    Clearly $\mathbf{w} \in W$. To verify $\mathbf{w}' \in W^\perp$: for any $j = 1, \ldots, k$,

    $$\langle \mathbf{w}', \mathbf{e}_j \rangle = \langle \mathbf{v} - \mathbf{w}, \mathbf{e}_j\rangle = \langle \mathbf{v}, \mathbf{e}_j\rangle - \sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i\rangle\langle \mathbf{e}_i, \mathbf{e}_j\rangle = \langle \mathbf{v}, \mathbf{e}_j\rangle - \langle \mathbf{v}, \mathbf{e}_j\rangle = 0$$

    So $\mathbf{w}' \perp \mathbf{e}_j$ for all $j$, hence $\mathbf{w}' \in W^\perp$. Thus $\mathbf{v} = \mathbf{w} + \mathbf{w}'$ is the desired decomposition.

    Uniqueness is guaranteed by $W \cap W^\perp = \{\mathbf{0}\}$: if $\mathbf{v} = \mathbf{w}_1 + \mathbf{w}_1' = \mathbf{w}_2 + \mathbf{w}_2'$, then $\mathbf{w}_1 - \mathbf{w}_2 = \mathbf{w}_2' - \mathbf{w}_1' \in W \cap W^\perp = \{\mathbf{0}\}$. $\blacksquare$

!!! proposition "Proposition 8.2 (Linear independence of orthogonal sets)"
    Let $S = \{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ be an orthogonal set of nonzero vectors in an inner product space $V$. Then $S$ is linearly independent.

??? proof "Proof"
    Suppose $\alpha_1\mathbf{v}_1 + \cdots + \alpha_k\mathbf{v}_k = \mathbf{0}$. For any $j$, take the inner product with $\mathbf{v}_j$:

    $$0 = \left\langle \sum_{i=1}^k \alpha_i\mathbf{v}_i, \mathbf{v}_j \right\rangle = \sum_{i=1}^k \alpha_i \langle \mathbf{v}_i, \mathbf{v}_j\rangle = \alpha_j \|\mathbf{v}_j\|^2$$

    Since $\mathbf{v}_j \neq \mathbf{0}$, we have $\|\mathbf{v}_j\|^2 > 0$, so $\alpha_j = 0$. $\blacksquare$

!!! example "Example 8.6"
    In $\mathbb{R}^3$, let $W = \operatorname{span}\{(1, 1, 0)^T, (0, 1, 1)^T\}$. Find $W^\perp$.

    Let $\mathbf{v} = (x, y, z)^T \in W^\perp$, then

    $$\langle \mathbf{v}, (1,1,0)^T\rangle = x + y = 0, \quad \langle \mathbf{v}, (0,1,1)^T\rangle = y + z = 0$$

    Solving gives $y = -x$, $z = -y = x$. Hence $W^\perp = \operatorname{span}\{(1, -1, 1)^T\}$.

    Verification: $\dim W + \dim W^\perp = 2 + 1 = 3 = \dim \mathbb{R}^3$.

---

## 8.4 Orthogonal Projection Theorem

<div class="context-flow" markdown>

Orthogonal direct sum $V = W \oplus W^\perp$ → Projection $\operatorname{proj}_W$ is the **best approximation** (minimizes distance) → The projection operator is both idempotent $P^2 = P$ and self-adjoint $P^* = P$

</div>

!!! definition "Definition 8.8 (Orthogonal projection)"
    Let $V$ be a finite-dimensional inner product space and $W$ a subspace of $V$. By the orthogonal direct sum decomposition $V = W \oplus W^\perp$, for any $\mathbf{v} \in V$, define the **orthogonal projection** $\operatorname{proj}_W: V \to V$ as

    $$\operatorname{proj}_W(\mathbf{v}) = \mathbf{w}$$

    where $\mathbf{v} = \mathbf{w} + \mathbf{w}'$, $\mathbf{w} \in W$, $\mathbf{w}' \in W^\perp$.

    If $\{\mathbf{e}_1, \ldots, \mathbf{e}_k\}$ is an orthonormal basis of $W$, then

    $$\operatorname{proj}_W(\mathbf{v}) = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i \rangle \mathbf{e}_i$$

<div class="context-flow" markdown>

**Insight**: Best approximation = orthogonal projection; the core of the proof is the Pythagorean theorem $\|\mathbf{v}-\mathbf{w}\|^2 = \|\mathbf{v}-\hat{\mathbf{v}}\|^2 + \|\hat{\mathbf{v}}-\mathbf{w}\|^2$ → leads directly to Ch11 least squares

</div>

!!! theorem "Theorem 8.7 (Best approximation theorem)"
    Let $V$ be an inner product space and $W$ a finite-dimensional subspace of $V$. For any $\mathbf{v} \in V$, $\operatorname{proj}_W(\mathbf{v})$ is the closest vector in $W$ to $\mathbf{v}$. That is, for any $\mathbf{w} \in W$,

    $$\|\mathbf{v} - \operatorname{proj}_W(\mathbf{v})\| \leq \|\mathbf{v} - \mathbf{w}\|$$

    Equality holds if and only if $\mathbf{w} = \operatorname{proj}_W(\mathbf{v})$.

??? proof "Proof"
    Let $\hat{\mathbf{v}} = \operatorname{proj}_W(\mathbf{v})$. For any $\mathbf{w} \in W$,

    $$\|\mathbf{v} - \mathbf{w}\|^2 = \|(\mathbf{v} - \hat{\mathbf{v}}) + (\hat{\mathbf{v}} - \mathbf{w})\|^2$$

    Note that $\mathbf{v} - \hat{\mathbf{v}} \in W^\perp$ and $\hat{\mathbf{v}} - \mathbf{w} \in W$, so they are orthogonal. By the Pythagorean theorem:

    $$\|\mathbf{v} - \mathbf{w}\|^2 = \|\mathbf{v} - \hat{\mathbf{v}}\|^2 + \|\hat{\mathbf{v}} - \mathbf{w}\|^2 \geq \|\mathbf{v} - \hat{\mathbf{v}}\|^2$$

    Equality holds if and only if $\|\hat{\mathbf{v}} - \mathbf{w}\|^2 = 0$, i.e., $\mathbf{w} = \hat{\mathbf{v}}$. $\blacksquare$

!!! proposition "Proposition 8.3 (Properties of the orthogonal projection operator)"
    Let $P = \operatorname{proj}_W$ be the orthogonal projection operator. Then

    1. $P$ is a linear operator;
    2. $P^2 = P$ (idempotent);
    3. $\operatorname{Im}(P) = W$, $\operatorname{Ker}(P) = W^\perp$;
    4. $P^* = P$ (self-adjoint);
    5. For any $\mathbf{v} \in V$, $\|\mathbf{v}\|^2 = \|P\mathbf{v}\|^2 + \|\mathbf{v} - P\mathbf{v}\|^2$.

??? proof "Proof"
    (1) Let $\mathbf{v}_1 = \mathbf{w}_1 + \mathbf{w}_1'$, $\mathbf{v}_2 = \mathbf{w}_2 + \mathbf{w}_2'$, then $\alpha\mathbf{v}_1 + \beta\mathbf{v}_2 = (\alpha\mathbf{w}_1 + \beta\mathbf{w}_2) + (\alpha\mathbf{w}_1' + \beta\mathbf{w}_2')$, so $P(\alpha\mathbf{v}_1 + \beta\mathbf{v}_2) = \alpha\mathbf{w}_1 + \beta\mathbf{w}_2 = \alpha P\mathbf{v}_1 + \beta P\mathbf{v}_2$.

    (2) For $\mathbf{v} = \mathbf{w} + \mathbf{w}'$, $P\mathbf{v} = \mathbf{w} \in W$, so $P(P\mathbf{v}) = P\mathbf{w} = \mathbf{w} = P\mathbf{v}$.

    (3) Obvious.

    (4) Let $\mathbf{u} = \mathbf{w}_1 + \mathbf{w}_1'$, $\mathbf{v} = \mathbf{w}_2 + \mathbf{w}_2'$. Then $\langle P\mathbf{u}, \mathbf{v}\rangle = \langle \mathbf{w}_1, \mathbf{w}_2 + \mathbf{w}_2'\rangle = \langle \mathbf{w}_1, \mathbf{w}_2\rangle$ and $\langle \mathbf{u}, P\mathbf{v}\rangle = \langle \mathbf{w}_1 + \mathbf{w}_1', \mathbf{w}_2\rangle = \langle \mathbf{w}_1, \mathbf{w}_2\rangle$.

    (5) Since $\mathbf{v} = P\mathbf{v} + (\mathbf{v} - P\mathbf{v})$ and the two are orthogonal, the result follows from the Pythagorean theorem. $\blacksquare$

!!! example "Example 8.7"
    Let $W = \operatorname{span}\left\{\mathbf{e}_1 = \frac{1}{\sqrt{2}}(1, 1, 0)^T, \, \mathbf{e}_2 = \frac{1}{\sqrt{6}}(1, -1, 2)^T\right\}$ be a subspace of $\mathbb{R}^3$ ($\mathbf{e}_1, \mathbf{e}_2$ are already an orthonormal basis). Find the orthogonal projection of $\mathbf{v} = (1, 2, 3)^T$ onto $W$.

    $$\operatorname{proj}_W(\mathbf{v}) = \langle \mathbf{v}, \mathbf{e}_1\rangle \mathbf{e}_1 + \langle \mathbf{v}, \mathbf{e}_2\rangle \mathbf{e}_2$$

    $$\langle \mathbf{v}, \mathbf{e}_1\rangle = \frac{1}{\sqrt{2}}(1+2) = \frac{3}{\sqrt{2}}$$

    $$\langle \mathbf{v}, \mathbf{e}_2\rangle = \frac{1}{\sqrt{6}}(1-2+6) = \frac{5}{\sqrt{6}}$$

    $$\operatorname{proj}_W(\mathbf{v}) = \frac{3}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\\0\end{pmatrix} + \frac{5}{\sqrt{6}} \cdot \frac{1}{\sqrt{6}}\begin{pmatrix}1\\-1\\2\end{pmatrix} = \frac{3}{2}\begin{pmatrix}1\\1\\0\end{pmatrix} + \frac{5}{6}\begin{pmatrix}1\\-1\\2\end{pmatrix} = \begin{pmatrix}7/3\\2/3\\5/3\end{pmatrix}$$

---

## 8.5 Bessel Inequality and Parseval's Identity

<div class="context-flow" markdown>

Projection only captures "partial energy" → **Bessel**: $\sum |\langle \mathbf{v},\mathbf{e}_i\rangle|^2 \le \|\mathbf{v}\|^2$; equality when the set is extended to a basis → **Parseval** → leads to Fourier analysis

</div>

!!! theorem "Theorem 8.8 (Bessel inequality)"
    Let $(V, \langle \cdot, \cdot \rangle)$ be an inner product space and $\{\mathbf{e}_1, \ldots, \mathbf{e}_k\}$ an orthonormal set in $V$. Then for any $\mathbf{v} \in V$,

    $$\sum_{i=1}^k |\langle \mathbf{v}, \mathbf{e}_i \rangle|^2 \leq \|\mathbf{v}\|^2$$

??? proof "Proof"
    Let $W = \operatorname{span}\{\mathbf{e}_1, \ldots, \mathbf{e}_k\}$, $\hat{\mathbf{v}} = \operatorname{proj}_W(\mathbf{v}) = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i\rangle\mathbf{e}_i$.

    By the Pythagorean theorem:

    $$\|\mathbf{v}\|^2 = \|\hat{\mathbf{v}}\|^2 + \|\mathbf{v} - \hat{\mathbf{v}}\|^2 \geq \|\hat{\mathbf{v}}\|^2$$

    And $\|\hat{\mathbf{v}}\|^2 = \left\|\sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i\rangle \mathbf{e}_i\right\|^2 = \sum_{i=1}^k |\langle \mathbf{v}, \mathbf{e}_i\rangle|^2$ (using orthonormality).

    Hence $\sum_{i=1}^k |\langle \mathbf{v}, \mathbf{e}_i\rangle|^2 \leq \|\mathbf{v}\|^2$. $\blacksquare$

!!! theorem "Theorem 8.9 (Parseval's identity)"
    Let $V$ be a finite-dimensional inner product space and $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ an orthonormal basis of $V$. Then for any $\mathbf{v} \in V$,

    $$\sum_{i=1}^n |\langle \mathbf{v}, \mathbf{e}_i \rangle|^2 = \|\mathbf{v}\|^2$$

    More generally, for any $\mathbf{u}, \mathbf{v} \in V$,

    $$\langle \mathbf{u}, \mathbf{v} \rangle = \sum_{i=1}^n \langle \mathbf{u}, \mathbf{e}_i \rangle \overline{\langle \mathbf{v}, \mathbf{e}_i \rangle}$$

??? proof "Proof"
    Since $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ is an orthonormal basis of $V$, $\mathbf{v} = \sum_{i=1}^n \langle \mathbf{v}, \mathbf{e}_i\rangle \mathbf{e}_i$, so

    $$\|\mathbf{v}\|^2 = \left\langle \sum_{i=1}^n \langle \mathbf{v}, \mathbf{e}_i\rangle \mathbf{e}_i, \sum_{j=1}^n \langle \mathbf{v}, \mathbf{e}_j\rangle \mathbf{e}_j \right\rangle = \sum_{i=1}^n \sum_{j=1}^n \langle \mathbf{v}, \mathbf{e}_i\rangle \overline{\langle \mathbf{v}, \mathbf{e}_j\rangle}\langle \mathbf{e}_i, \mathbf{e}_j\rangle = \sum_{i=1}^n |\langle \mathbf{v}, \mathbf{e}_i\rangle|^2$$

    The general identity (the general form of Parseval's identity) is proved similarly. $\blacksquare$

!!! note "Note"
    Parseval's identity can be understood as: the sum of the squared magnitudes of the coordinate components of a vector in an orthonormal basis equals the squared norm of the vector. Bessel's inequality says that if the orthonormal set is not a basis (i.e., does not span the entire space), then the "energy" (squared norm) cannot be fully captured.

!!! example "Example 8.8"
    In $\mathbb{R}^3$ with the standard orthonormal basis $\{\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3\}$, let $\mathbf{v} = (3, 4, 0)^T$.

    Parseval's identity: $|\langle \mathbf{v}, \mathbf{e}_1\rangle|^2 + |\langle \mathbf{v}, \mathbf{e}_2\rangle|^2 + |\langle \mathbf{v}, \mathbf{e}_3\rangle|^2 = 9 + 16 + 0 = 25 = \|\mathbf{v}\|^2$.

    Using only the first two basis vectors (Bessel inequality): $|\langle \mathbf{v}, \mathbf{e}_1\rangle|^2 + |\langle \mathbf{v}, \mathbf{e}_2\rangle|^2 = 25 \leq 25 = \|\mathbf{v}\|^2$. In this example equality holds because $\mathbf{v}$ itself lies in the subspace spanned by the first two components.

---

## 8.6 Adjoint Operators

<div class="context-flow" markdown>

"Transfer" the inner product to the other side: $\langle T\mathbf{v},\mathbf{w}\rangle = \langle \mathbf{v},T^*\mathbf{w}\rangle$ → In an orthonormal basis, the matrix of $T^*$ is $A^H$ → connects to the decomposition theory in Ch10

</div>

!!! definition "Definition 8.9 (Adjoint operator)"
    Let $V$ be a finite-dimensional inner product space and $T: V \to V$ a linear operator. The **adjoint operator** $T^*: V \to V$ of $T$ is the unique linear operator satisfying

    $$\langle T\mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{v}, T^*\mathbf{w} \rangle, \quad \forall \mathbf{v}, \mathbf{w} \in V$$

!!! theorem "Theorem 8.10 (Existence and uniqueness of the adjoint)"
    Let $V$ be a finite-dimensional inner product space and $T: V \to V$ a linear operator. Then $T^*$ exists and is unique.

??? proof "Proof"
    **Existence:** Let $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ be an orthonormal basis of $V$. The matrix of $T$ in this basis is $A = (a_{ij})$, where $a_{ij} = \langle T\mathbf{e}_j, \mathbf{e}_i\rangle$.

    Define $T^*$ as the operator represented by the matrix $B = (b_{ij})$, where $b_{ij} = \overline{a_{ji}}$, i.e., $B = A^H$ (conjugate transpose). Then

    $$\langle T\mathbf{v}, \mathbf{w}\rangle = (A\mathbf{v})^H \mathbf{w} = \mathbf{v}^H A^H \mathbf{w} = \mathbf{v}^H (B\mathbf{w}) = \langle \mathbf{v}, T^*\mathbf{w}\rangle$$

    **Uniqueness:** If $S$ also satisfies $\langle T\mathbf{v}, \mathbf{w}\rangle = \langle \mathbf{v}, S\mathbf{w}\rangle$, then $\langle \mathbf{v}, T^*\mathbf{w} - S\mathbf{w}\rangle = 0$ for all $\mathbf{v}$. Taking $\mathbf{v} = T^*\mathbf{w} - S\mathbf{w}$ gives $T^*\mathbf{w} = S\mathbf{w}$ for all $\mathbf{w}$. $\blacksquare$

!!! theorem "Theorem 8.11 (Properties of the adjoint operator)"
    Let $V$ be a finite-dimensional inner product space and $S, T: V \to V$ linear operators, $\alpha$ a scalar. Then

    1. $(S + T)^* = S^* + T^*$
    2. $(\alpha T)^* = \bar{\alpha} T^*$
    3. $(ST)^* = T^* S^*$
    4. $(T^*)^* = T$
    5. $\operatorname{Ker}(T^*) = (\operatorname{Im}(T))^\perp$
    6. $\operatorname{Im}(T^*) = (\operatorname{Ker}(T))^\perp$

??? proof "Proof"
    We prove property 5; the others are similar.

    $\mathbf{w} \in \operatorname{Ker}(T^*)$ $\Leftrightarrow$ $T^*\mathbf{w} = \mathbf{0}$ $\Leftrightarrow$ $\langle \mathbf{v}, T^*\mathbf{w}\rangle = 0, \, \forall \mathbf{v} \in V$ $\Leftrightarrow$ $\langle T\mathbf{v}, \mathbf{w}\rangle = 0, \, \forall \mathbf{v} \in V$ $\Leftrightarrow$ $\mathbf{w} \in (\operatorname{Im}(T))^\perp$. $\blacksquare$

!!! note "Note"
    In an orthonormal basis, when the matrix of $T$ is $A$, the matrix of $T^*$ is $A^H = \overline{A}^T$ (conjugate transpose). For real inner product spaces, the matrix of $T^*$ is simply $A^T$.

!!! example "Example 8.9"
    Let $T: \mathbb{C}^2 \to \mathbb{C}^2$ have matrix in the standard basis

    $$A = \begin{pmatrix} 1+i & 2 \\ 3i & 4-i \end{pmatrix}$$

    Then the matrix of $T^*$ is

    $$A^H = \overline{A}^T = \begin{pmatrix} 1-i & -3i \\ 2 & 4+i \end{pmatrix}$$

!!! example "Example 8.10"
    Let $V = C[0,1]$ with inner product $\langle f, g\rangle = \int_0^1 f(t)\overline{g(t)} \, dt$. Define the linear operator $T: V \to V$ by $(Tf)(t) = tf(t)$ (multiplication by the variable). Then $T^* = T$, because

    $$\langle Tf, g\rangle = \int_0^1 tf(t)\overline{g(t)} \, dt = \int_0^1 f(t)\overline{tg(t)} \, dt = \langle f, Tg\rangle$$

    where the last step uses the fact that $t$ is real.

---

## 8.7 Normal Operators

<div class="context-flow" markdown>

$TT^* = T^*T$ (commutes with its adjoint) → Eigenvectors for distinct eigenvalues are automatically orthogonal → Normal = necessary and sufficient condition for **unitarily diagonalizable** (complex case)

</div>

!!! definition "Definition 8.10 (Normal operator)"
    Let $V$ be a finite-dimensional inner product space and $T: V \to V$ a linear operator. If

    $$TT^* = T^*T$$

    then $T$ is called a **normal operator**. Correspondingly, if an $n \times n$ matrix $A$ satisfies $AA^H = A^HA$, then $A$ is called a **normal matrix**.

!!! note "Note"
    Self-adjoint operators ($T^* = T$), unitary operators ($T^*T = I$), and skew-adjoint operators ($T^* = -T$) are all special cases of normal operators.

!!! theorem "Theorem 8.12 (Equivalent conditions for normal operators)"
    Let $V$ be a finite-dimensional inner product space and $T: V \to V$ a linear operator. The following conditions are equivalent:

    1. $T$ is normal;
    2. $\|T\mathbf{v}\| = \|T^*\mathbf{v}\|$ for all $\mathbf{v} \in V$;
    3. $T - \lambda I$ is normal for every scalar $\lambda$;
    4. If $T\mathbf{v} = \lambda\mathbf{v}$, then $T^*\mathbf{v} = \bar{\lambda}\mathbf{v}$;
    5. Eigenvectors of $T$ corresponding to distinct eigenvalues are orthogonal.

??? proof "Proof"
    **(1)$\Rightarrow$(2):** $\|T\mathbf{v}\|^2 = \langle T\mathbf{v}, T\mathbf{v}\rangle = \langle T^*T\mathbf{v}, \mathbf{v}\rangle = \langle TT^*\mathbf{v}, \mathbf{v}\rangle = \langle T^*\mathbf{v}, T^*\mathbf{v}\rangle = \|T^*\mathbf{v}\|^2$.

    **(2)$\Rightarrow$(1):** Condition (2) implies $\langle (T^*T - TT^*)\mathbf{v}, \mathbf{v}\rangle = 0$ for all $\mathbf{v}$. Since $T^*T - TT^*$ is self-adjoint, this implies $T^*T - TT^* = 0$.

    **(1)$\Rightarrow$(3):** $(T-\lambda I)^* = T^* - \bar{\lambda}I$. Direct verification:

    $$(T-\lambda I)(T-\lambda I)^* = TT^* - \lambda T^* - \bar{\lambda}T + |\lambda|^2 I$$

    $$(T-\lambda I)^*(T-\lambda I) = T^*T - \bar{\lambda}T - \lambda T^* + |\lambda|^2 I$$

    By $TT^* = T^*T$, the two are equal.

    **(3)$\Rightarrow$(4):** If $T\mathbf{v} = \lambda\mathbf{v}$, i.e., $(T-\lambda I)\mathbf{v} = \mathbf{0}$. By (3), $T-\lambda I$ is normal, and by (2), $\|(T-\lambda I)^*\mathbf{v}\| = \|(T-\lambda I)\mathbf{v}\| = 0$, i.e., $(T^*-\bar{\lambda}I)\mathbf{v} = \mathbf{0}$, so $T^*\mathbf{v} = \bar{\lambda}\mathbf{v}$.

    **(4)$\Rightarrow$(5):** Let $T\mathbf{v}_1 = \lambda_1\mathbf{v}_1$, $T\mathbf{v}_2 = \lambda_2\mathbf{v}_2$, $\lambda_1 \neq \lambda_2$. By (4), $T^*\mathbf{v}_1 = \bar{\lambda}_1\mathbf{v}_1$. Then

    $$\lambda_1\langle \mathbf{v}_1, \mathbf{v}_2\rangle = \langle T\mathbf{v}_1, \mathbf{v}_2\rangle = \langle \mathbf{v}_1, T^*\mathbf{v}_2\rangle = \langle \mathbf{v}_1, \bar{\lambda}_2\mathbf{v}_2\rangle = \lambda_2\langle \mathbf{v}_1, \mathbf{v}_2\rangle$$

    So $(\lambda_1 - \lambda_2)\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$. Since $\lambda_1 \neq \lambda_2$, we get $\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$. $\blacksquare$

!!! example "Example 8.11"
    The matrix $A = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$ is a normal matrix, since

    $$AA^T = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$$

    $$A^TA = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$$

    But $A$ is not symmetric ($A \neq A^T$), nor is it orthogonal ($A^TA \neq I$).

---

## 8.8 Self-Adjoint Operators

<div class="context-flow" markdown>

The most important special case of normal operators: $T^* = T$ → Eigenvalues are **all real** → Eigenvectors for distinct eigenvalues are orthogonal → leads directly to the spectral theorem

</div>

!!! definition "Definition 8.11 (Self-adjoint operator)"
    Let $V$ be a finite-dimensional inner product space and $T: V \to V$ a linear operator. If $T^* = T$, i.e.,

    $$\langle T\mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{v}, T\mathbf{w} \rangle, \quad \forall \mathbf{v}, \mathbf{w} \in V$$

    then $T$ is called a **self-adjoint operator**. In real inner product spaces it is also called a **symmetric operator**. Correspondingly, a matrix satisfying $A^H = A$ is called a **Hermitian matrix** (in the real case, a symmetric matrix $A^T = A$).

!!! theorem "Theorem 8.13 (Eigenvalues of a self-adjoint operator are all real)"
    Let $V$ be a finite-dimensional inner product space (real or complex) and $T: V \to V$ a self-adjoint operator. Then all eigenvalues of $T$ are real.

??? proof "Proof"
    Let $T\mathbf{v} = \lambda\mathbf{v}$, $\mathbf{v} \neq \mathbf{0}$. Then

    $$\lambda\|\mathbf{v}\|^2 = \lambda\langle \mathbf{v}, \mathbf{v}\rangle = \langle \lambda\mathbf{v}, \mathbf{v}\rangle = \langle T\mathbf{v}, \mathbf{v}\rangle = \langle \mathbf{v}, T\mathbf{v}\rangle = \langle \mathbf{v}, \lambda\mathbf{v}\rangle = \bar{\lambda}\langle \mathbf{v}, \mathbf{v}\rangle = \bar{\lambda}\|\mathbf{v}\|^2$$

    Since $\|\mathbf{v}\|^2 > 0$, we get $\lambda = \bar{\lambda}$, i.e., $\lambda \in \mathbb{R}$. $\blacksquare$

!!! theorem "Theorem 8.14 (Eigenvectors of a self-adjoint operator for distinct eigenvalues are orthogonal)"
    Let $T$ be a self-adjoint operator, $\lambda_1 \neq \lambda_2$ two distinct eigenvalues of $T$, and $\mathbf{v}_1, \mathbf{v}_2$ corresponding eigenvectors. Then $\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$.

??? proof "Proof"
    By Theorem 8.12 property (5) (self-adjoint operators are a special case of normal operators), or by direct proof:

    $$\lambda_1\langle \mathbf{v}_1, \mathbf{v}_2\rangle = \langle T\mathbf{v}_1, \mathbf{v}_2\rangle = \langle \mathbf{v}_1, T\mathbf{v}_2\rangle = \lambda_2\langle \mathbf{v}_1, \mathbf{v}_2\rangle$$

    where the last step uses $\lambda_2 \in \mathbb{R}$ (so $\bar{\lambda}_2 = \lambda_2$). Therefore $(\lambda_1 - \lambda_2)\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$, and since $\lambda_1 \neq \lambda_2$, we get $\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$. $\blacksquare$

!!! proposition "Proposition 8.4"
    Let $V$ be a finite-dimensional **complex** inner product space and $T: V \to V$ a linear operator. If $\langle T\mathbf{v}, \mathbf{v}\rangle = 0$ for all $\mathbf{v} \in V$, then $T = 0$.

??? proof "Proof"
    For any $\mathbf{u}, \mathbf{v} \in V$, using the polarization identity:

    $$\langle T\mathbf{u}, \mathbf{v}\rangle = \frac{1}{4}\sum_{k=0}^3 i^k \langle T(\mathbf{u}+i^k\mathbf{v}), \mathbf{u}+i^k\mathbf{v}\rangle = 0$$

    Taking $\mathbf{v} = T\mathbf{u}$ gives $\|T\mathbf{u}\|^2 = 0$, so $T\mathbf{u} = \mathbf{0}$ for all $\mathbf{u}$. $\blacksquare$

!!! note "Note"
    Proposition 8.4 does not hold in real inner product spaces. For example, $T: \mathbb{R}^2 \to \mathbb{R}^2$, $T(x,y)^T = (-y, x)^T$ (rotation by 90 degrees), satisfies $\langle T\mathbf{v}, \mathbf{v}\rangle = 0$ for all $\mathbf{v}$, but $T \neq 0$.

!!! example "Example 8.12"
    The matrix $A = \begin{pmatrix} 2 & 1+i \\ 1-i & 3 \end{pmatrix}$ is a Hermitian matrix since $A^H = A$. Its characteristic polynomial is

    $$\det(A - \lambda I) = (2-\lambda)(3-\lambda) - |1+i|^2 = \lambda^2 - 5\lambda + 4 = (\lambda - 1)(\lambda - 4)$$

    The eigenvalues are $\lambda_1 = 1, \lambda_2 = 4$, both real.

---

## 8.9 Spectral Theorem

<div class="context-flow" markdown>

**Chapter climax**: Inner product + self-adjoint/normal → There exists an **orthonormal eigenbasis** → $A = Q\Lambda Q^T$ (real symmetric) / $A = U\Lambda U^H$ (normal)

→ Ch9 uses this to diagonalize quadratic forms, Ch10 uses this for spectral decomposition, Ch11 SVD is essentially "spectral theorem on both sides"

</div>

The spectral theorem is one of the most important and profound theorems in linear algebra. It states that normal operators (including self-adjoint operators) can be diagonalized with respect to an orthonormal basis.

### Real Spectral Theorem

<div class="context-flow" markdown>

**Insight**: The key to the inductive proof of the spectral theorem — after finding one eigenvector, its orthogonal complement is an **invariant subspace**, allowing recursive dimension reduction

</div>

!!! theorem "Theorem 8.15 (Spectral theorem for real symmetric matrices)"
    Let $A$ be an $n \times n$ real symmetric matrix (i.e., $A^T = A$). Then

    1. All eigenvalues of $A$ are real;
    2. Eigenvectors of $A$ corresponding to distinct eigenvalues are orthogonal;
    3. There exists an orthogonal matrix $Q$ (i.e., $Q^TQ = I$) such that $A = Q\Lambda Q^T$, where $\Lambda$ is a diagonal matrix.

    Equivalently, $\mathbb{R}^n$ has an orthonormal basis consisting of eigenvectors of $A$.

??? proof "Proof"
    Properties 1 and 2 have been proved in Theorems 8.13 and 8.14. We now prove property 3.

    By induction on $n$. The case $n=1$ is trivial. Assume the result holds for $(n-1) \times (n-1)$ matrices.

    **Key step:** First we need to show that the real symmetric matrix $A$ has a real eigenvalue. The characteristic polynomial $\det(A - \lambda I)$ is a degree $n$ real polynomial. Over $\mathbb{C}$ it has $n$ roots (counting multiplicity), and by Theorem 8.13 all roots are real.

    Let $\lambda_1$ be a (real) eigenvalue of $A$ and $\mathbf{q}_1$ a corresponding unit eigenvector. Extend $\mathbf{q}_1$ to an orthonormal basis $\{\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n\}$ of $\mathbb{R}^n$, and let $Q_1 = [\mathbf{q}_1 | \mathbf{q}_2 | \cdots | \mathbf{q}_n]$. Then

    $$Q_1^T A Q_1 = \begin{pmatrix} \lambda_1 & \mathbf{b}^T \\ \mathbf{0} & A_1 \end{pmatrix}$$

    By the symmetry of $A$, $Q_1^TAQ_1$ is also symmetric, so $\mathbf{b} = \mathbf{0}$ and $A_1$ is an $(n-1)\times(n-1)$ real symmetric matrix.

    By the induction hypothesis, there exists an $(n-1) \times (n-1)$ orthogonal matrix $P_1$ such that $P_1^TA_1P_1 = \operatorname{diag}(\lambda_2, \ldots, \lambda_n)$. Let

    $$P = \begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & P_1 \end{pmatrix}$$

    Then $Q = Q_1P$ is an orthogonal matrix and $Q^TAQ = \operatorname{diag}(\lambda_1, \lambda_2, \ldots, \lambda_n)$. $\blacksquare$

### Complex Spectral Theorem

!!! theorem "Theorem 8.16 (Spectral theorem for Hermitian matrices)"
    Let $A$ be an $n \times n$ Hermitian matrix (i.e., $A^H = A$). Then there exists a unitary matrix $U$ (i.e., $U^HU = I$) such that

    $$A = U\Lambda U^H$$

    where $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$ and all $\lambda_i$ are real.

??? proof "Proof"
    The proof is completely analogous to the real symmetric case, replacing "orthogonal" by "unitary" and "transpose" by "conjugate transpose." $\blacksquare$

!!! theorem "Theorem 8.17 (Spectral theorem for normal operators — complex case)"
    Let $V$ be a finite-dimensional **complex** inner product space and $T: V \to V$ a linear operator. Then $T$ is normal if and only if there exists an orthonormal basis of $V$ with respect to which the matrix of $T$ is diagonal.

    Equivalently, an $n \times n$ complex matrix $A$ is normal if and only if there exists a unitary matrix $U$ such that

    $$A = U\Lambda U^H$$

    where $\Lambda$ is a diagonal matrix (with possibly complex diagonal entries).

??? proof "Proof"
    **Sufficiency:** If $A = UDU^H$ ($D$ diagonal), then $A^H = UD^HU^H$, so

    $$AA^H = UDD^HU^H = UD^HDU^H = A^HA$$

    (since diagonal matrices always commute.)

    **Necessity:** By induction on $n$. The case $n = 1$ is trivial. Assume the result holds for $n-1$.

    By the fundamental theorem of algebra, $T$ has an eigenvalue $\lambda_1$ over $\mathbb{C}$. Let $\mathbf{e}_1$ be a corresponding unit eigenvector.

    Let $W = \operatorname{span}\{\mathbf{e}_1\}$. By Theorem 8.12 property (4), $T^*\mathbf{e}_1 = \bar{\lambda}_1\mathbf{e}_1$, so $W$ is $T^*$-invariant.

    We show $W^\perp$ is $T$-invariant: if $\mathbf{v} \in W^\perp$, then for any $\mathbf{w} \in W$,

    $$\langle T\mathbf{v}, \mathbf{w}\rangle = \langle \mathbf{v}, T^*\mathbf{w}\rangle$$

    Since $W$ is $T^*$-invariant, $T^*\mathbf{w} \in W$, so $\langle \mathbf{v}, T^*\mathbf{w}\rangle = 0$. Thus $T\mathbf{v} \in W^\perp$.

    $T|_{W^\perp}$ is a normal operator on $W^\perp$ (normality is inherited). By the induction hypothesis, there exists an orthonormal basis $\{\mathbf{e}_2, \ldots, \mathbf{e}_n\}$ of $W^\perp$ such that $T|_{W^\perp}$ is diagonal in this basis. Then $\{\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\}$ is an orthonormal basis of $V$ in which $T$ is diagonal. $\blacksquare$

!!! note "Note"
    The **real spectral theorem** requires additional conditions. A general real normal matrix cannot necessarily be orthogonally diagonalized over the reals (e.g., rotation matrices have complex eigenvalues). A real normal matrix can be orthogonally similar to a block diagonal form: diagonal blocks are $1\times 1$ real numbers or $2\times 2$ rotation-scaling blocks $\begin{pmatrix} a & -b \\ b & a \end{pmatrix}$ ($b \neq 0$), corresponding to conjugate complex eigenvalue pairs $a \pm bi$.

!!! example "Example 8.13"
    Orthogonally diagonalize the real symmetric matrix $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$.

    Characteristic polynomial: $\det(A - \lambda I) = (2-\lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda-1)(\lambda-3)$.

    $\lambda_1 = 1$: $(A-I)\mathbf{x} = \mathbf{0}$ gives $\mathbf{v}_1 = \frac{1}{\sqrt{2}}(1, -1)^T$.

    $\lambda_2 = 3$: $(A-3I)\mathbf{x} = \mathbf{0}$ gives $\mathbf{v}_2 = \frac{1}{\sqrt{2}}(1, 1)^T$.

    Verify orthogonality: $\langle \mathbf{v}_1, \mathbf{v}_2\rangle = \frac{1}{2}(1-1) = 0$.

    Let $Q = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$, then $A = Q\begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}Q^T$.

!!! example "Example 8.14"
    Let $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$. $A$ is a real normal matrix ($AA^T = A^TA = I$), but its eigenvalues are $\pm i$, which are not real. Therefore $A$ cannot be orthogonally diagonalized over the reals, but it can be unitarily diagonalized:

    $$A = U\begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}U^H, \quad U = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ i & -i \end{pmatrix}$$

---

## 8.10 Unitary and Orthogonal Operators

<div class="context-flow" markdown>

$T^*T = I$ $\leftrightarrow$ preserves inner product/norm/orthonormal bases → Eigenvalues satisfy $|\lambda|=1$ → Orthogonal matrices have $\det = \pm 1$ (rotation or reflection) → The "rotation part" in Ch10 QR decomposition / Ch11 SVD

</div>

!!! definition "Definition 8.12 (Unitary and orthogonal operators)"
    Let $V$ be a finite-dimensional inner product space and $T: V \to V$ a linear operator.

    - If $V$ is a complex inner product space and $T^*T = TT^* = I$ (i.e., $T^* = T^{-1}$), then $T$ is called a **unitary operator**.
    - If $V$ is a real inner product space and $T^*T = TT^* = I$, then $T$ is called an **orthogonal operator**.

!!! definition "Definition 8.13 (Unitary and orthogonal matrices)"
    - An $n \times n$ complex matrix $U$ satisfying $U^HU = I$ (equivalently $UU^H = I$) is called a **unitary matrix**.
    - An $n \times n$ real matrix $Q$ satisfying $Q^TQ = I$ (equivalently $QQ^T = I$) is called an **orthogonal matrix**.

!!! theorem "Theorem 8.18 (Equivalent conditions for unitary/orthogonal operators)"
    Let $V$ be a finite-dimensional inner product space and $T: V \to V$ a linear operator. The following conditions are equivalent:

    1. $T$ is a unitary (orthogonal) operator;
    2. $T$ preserves the inner product: $\langle T\mathbf{u}, T\mathbf{v}\rangle = \langle \mathbf{u}, \mathbf{v}\rangle$, $\forall \mathbf{u}, \mathbf{v} \in V$;
    3. $T$ preserves the norm (isometry): $\|T\mathbf{v}\| = \|\mathbf{v}\|$, $\forall \mathbf{v} \in V$;
    4. $T$ maps orthonormal bases to orthonormal bases;
    5. The columns of $T$ form an orthonormal set.

??? proof "Proof"
    **(1)$\Rightarrow$(2):** $\langle T\mathbf{u}, T\mathbf{v}\rangle = \langle \mathbf{u}, T^*T\mathbf{v}\rangle = \langle \mathbf{u}, I\mathbf{v}\rangle = \langle \mathbf{u}, \mathbf{v}\rangle$.

    **(2)$\Rightarrow$(3):** Take $\mathbf{u} = \mathbf{v}$.

    **(3)$\Rightarrow$(1):** By the polarization identity, preserving the norm implies preserving the inner product. That is, $\langle T^*T\mathbf{v}, \mathbf{v}\rangle = \langle \mathbf{v}, \mathbf{v}\rangle$ for all $\mathbf{v}$, so $\langle (T^*T - I)\mathbf{v}, \mathbf{v}\rangle = 0$. Since $T^*T - I$ is self-adjoint, by reasoning similar to Proposition 8.4 (which holds for self-adjoint operators), $T^*T = I$.

    **(2)$\Leftrightarrow$(4):** Let $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ be an orthonormal basis. Then $\langle T\mathbf{e}_i, T\mathbf{e}_j\rangle = \langle \mathbf{e}_i, \mathbf{e}_j\rangle = \delta_{ij}$, so $\{T\mathbf{e}_1, \ldots, T\mathbf{e}_n\}$ is an orthonormal basis. The converse also holds.

    **(4)$\Leftrightarrow$(5):** Directly from the relationship between matrix columns and images of basis vectors. $\blacksquare$

!!! theorem "Theorem 8.19 (Eigenvalues of unitary/orthogonal operators)"
    Let $T$ be a unitary (orthogonal) operator. Then all eigenvalues of $T$ have modulus equal to $1$, i.e., $|\lambda| = 1$.

??? proof "Proof"
    Let $T\mathbf{v} = \lambda\mathbf{v}$, $\mathbf{v} \neq \mathbf{0}$. By norm preservation:

    $$\|\mathbf{v}\| = \|T\mathbf{v}\| = \|\lambda\mathbf{v}\| = |\lambda|\|\mathbf{v}\|$$

    Since $\|\mathbf{v}\| > 0$, we get $|\lambda| = 1$. $\blacksquare$

!!! proposition "Proposition 8.5 (Determinant of orthogonal matrices)"
    If $Q$ is an orthogonal matrix, then $\det(Q) = \pm 1$. If $U$ is a unitary matrix, then $|\det(U)| = 1$.

??? proof "Proof"
    $Q^TQ = I$ $\Rightarrow$ $\det(Q^T)\det(Q) = 1$ $\Rightarrow$ $(\det Q)^2 = 1$ $\Rightarrow$ $\det Q = \pm 1$.

    For unitary matrices: $U^HU = I$ $\Rightarrow$ $\overline{\det U} \cdot \det U = 1$ $\Rightarrow$ $|\det U|^2 = 1$ $\Rightarrow$ $|\det U| = 1$. $\blacksquare$

!!! example "Example 8.15"
    **General form of $2 \times 2$ orthogonal matrices.** Every $2 \times 2$ orthogonal matrix is either a rotation matrix

    $$R_\theta = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}, \quad \det R_\theta = 1$$

    or a reflection matrix

    $$S_\theta = \begin{pmatrix} \cos\theta & \sin\theta \\ \sin\theta & -\cos\theta \end{pmatrix}, \quad \det S_\theta = -1$$

    Orthogonal matrices with $\det Q = 1$ are called **special orthogonal matrices**, and they form the special orthogonal group $SO(n)$.
