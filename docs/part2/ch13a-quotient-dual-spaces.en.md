# Chapter 13A  Quotient Spaces and Dual Spaces

<div class="context-flow" markdown>

**Prerequisites**: Vector spaces and subspaces (Chapter 4) · Linear maps (Chapter 5) · Dimension formula (Chapter 5)

**Chapter arc**: Equivalence classes and quotient spaces → Dimension formula and isomorphism theorems → Linear functionals and dual spaces → Dual basis → Annihilators → Transpose maps → Double dual and canonical isomorphism → Applications in finite dimensions

**Further connections**：Quotient spaces are widely used in algebraic topology (construction of homology groups) and functional analysis (quotient Banach spaces); dual spaces underlie distribution theory (Schwartz distributions = dual of test function spaces) and weak topologies; in algebraic geometry, duality manifests as Serre duality

</div>

Quotient spaces and dual spaces are two deep and elegant constructions in linear algebra. The quotient space simplifies the structure of a vector space by "collapsing a subspace to a point," embodying the algebraic ideas of equivalence relations and congruence. The dual space examines a vector space from the perspective of "measurement," organizing linear functionals themselves into a new vector space. Their interplay — annihilators, transpose maps, the double dual — reveals the deep symmetry underlying linear algebra.

---

## 13A.1 Quotient Space Construction

<div class="context-flow" markdown>

**Core question**: Given a subspace $W$, can we treat vectors in $V$ that are "equivalent modulo $W$" as a single element? → Cosets → The quotient space $V/W$

</div>

### Cosets and equivalence relations

!!! definition "Definition 13A.1 (Coset)"
    Let $V$ be a vector space over $\mathbb{F}$ and $W$ a subspace of $V$. For $v \in V$, the set

    $$
    v + W = \{v + w : w \in W\}
    $$

    is called the **coset** (or **affine subset**) of $v$ with respect to $W$. The element $v$ is called a **representative** of the coset.

!!! theorem "Theorem 13A.1 (Equivalent characterizations of cosets)"
    Let $u, v \in V$ and $W$ be a subspace of $V$. The following are equivalent:

    1. $u + W = v + W$;
    2. $u - v \in W$;
    3. $u \in v + W$.

??? proof "Proof"
    $(1) \Rightarrow (3)$: $u = u + 0 \in u + W = v + W$.

    $(3) \Rightarrow (2)$: $u \in v + W$ means $u = v + w$ for some $w \in W$, so $u - v = w \in W$.

    $(2) \Rightarrow (1)$: Let $u - v = w_0 \in W$. For any $u + w \in u + W$, we have $u + w = v + w_0 + w \in v + W$ (since $w_0 + w \in W$). Similarly $v + W \subseteq u + W$. $\blacksquare$

!!! note "Note"
    The essence of cosets is an equivalence relation: define $u \sim v$ if and only if $u - v \in W$. Then $\sim$ is an equivalence relation on $V$, and $v + W$ is precisely the equivalence class of $v$.

!!! definition "Definition 13A.2 (Quotient space)"
    Let $W$ be a subspace of $V$. The set of all cosets

    $$
    V/W = \{v + W : v \in V\}
    $$

    equipped with the operations:

    - **Addition**: $(u + W) + (v + W) = (u + v) + W$;
    - **Scalar multiplication**: $\lambda(v + W) = \lambda v + W$, $\lambda \in \mathbb{F}$,

    forms a vector space over $\mathbb{F}$, called the **quotient space** of $V$ by $W$. The zero element of $V/W$ is $0 + W = W$.

!!! theorem "Theorem 13A.2 (Well-definedness of quotient space operations)"
    The addition and scalar multiplication on $V/W$ are independent of the choice of representatives.

??? proof "Proof"
    Suppose $u + W = u' + W$ and $v + W = v' + W$. Then $u - u' \in W$ and $v - v' \in W$.

    **Addition**: $(u + v) - (u' + v') = (u - u') + (v - v') \in W$, so $(u + v) + W = (u' + v') + W$.

    **Scalar multiplication**: $\lambda u - \lambda u' = \lambda(u - u') \in W$, so $\lambda u + W = \lambda u' + W$. $\blacksquare$

!!! example "Example 13A.1"
    Let $V = \mathbb{R}^3$ and $W = \{(x, y, 0) : x, y \in \mathbb{R}\}$ (the $xy$-plane).

    The coset $(a, b, c) + W = \{(x, y, c) : x, y \in \mathbb{R}\}$, i.e., the horizontal plane at height $c$. Two cosets are equal if and only if their third coordinates agree. Thus $V/W \cong \mathbb{R}$, with the isomorphism $(a, b, c) + W \mapsto c$.

!!! example "Example 13A.2"
    Let $V = \mathbb{R}^2$ and $W = \{(t, t) : t \in \mathbb{R}\}$ (the line through the origin with slope $1$).

    The coset $(a, b) + W = \{(a + t, b + t) : t \in \mathbb{R}\}$, a family of lines with slope $1$. Two cosets coincide if and only if $b - a$ is the same. Hence $V/W \cong \mathbb{R}$.

### The quotient map

!!! definition "Definition 13A.3 (Quotient map)"
    The **quotient map** $\pi: V \to V/W$ is defined by

    $$
    \pi(v) = v + W.
    $$

    It is a surjective linear map with $\ker \pi = W$.

!!! theorem "Theorem 13A.3 (Properties of the quotient map)"
    The quotient map $\pi: V \to V/W$ satisfies:

    1. $\pi$ is a linear map;
    2. $\pi$ is surjective;
    3. $\ker \pi = W$.

??? proof "Proof"
    1. $\pi(u + v) = (u + v) + W = (u + W) + (v + W) = \pi(u) + \pi(v)$; $\pi(\lambda v) = \lambda v + W = \lambda(v + W) = \lambda \pi(v)$.

    2. For any $v + W \in V/W$, $\pi(v) = v + W$.

    3. $\pi(v) = 0_{V/W} = W$ if and only if $v + W = W$ if and only if $v \in W$. $\blacksquare$

---

## 13A.2 Dimension and Isomorphism Theorems

<div class="context-flow" markdown>

**Core question**: What is the dimension of $V/W$? → Dimension formula → The first isomorphism theorem — one of the most important structural results in linear algebra

</div>

!!! theorem "Theorem 13A.4 (Dimension of the quotient space)"
    Let $V$ be a finite-dimensional vector space and $W$ a subspace of $V$. Then

    $$
    \dim(V/W) = \dim V - \dim W.
    $$

??? proof "Proof"
    The quotient map $\pi: V \to V/W$ is surjective with $\ker \pi = W$. By the rank-nullity theorem:

    $$
    \dim V = \dim \ker \pi + \dim \operatorname{im} \pi = \dim W + \dim(V/W).
    $$

    Therefore $\dim(V/W) = \dim V - \dim W$. $\blacksquare$

!!! note "Note"
    $\dim(V/W)$ is also called the **codimension** of $W$ in $V$, denoted $\operatorname{codim} W$.

!!! theorem "Theorem 13A.5 (First isomorphism theorem)"
    Let $T: V \to U$ be a linear map. Then there exists a unique isomorphism

    $$
    \widetilde{T}: V/\ker T \xrightarrow{\;\sim\;} \operatorname{im} T,
    $$

    such that $\widetilde{T}(v + \ker T) = T(v)$. That is, $V/\ker T \cong \operatorname{im} T$.

??? proof "Proof"
    **Well-definedness**: If $v + \ker T = v' + \ker T$, then $v - v' \in \ker T$, so $T(v) = T(v')$.

    **Linearity**: $\widetilde{T}((u + \ker T) + (v + \ker T)) = \widetilde{T}((u+v) + \ker T) = T(u+v) = T(u) + T(v) = \widetilde{T}(u + \ker T) + \widetilde{T}(v + \ker T)$. Scalar multiplication is similar.

    **Injectivity**: $\widetilde{T}(v + \ker T) = 0$ means $T(v) = 0$, i.e., $v \in \ker T$, so $v + \ker T = 0 + \ker T$.

    **Surjectivity**: For any $T(v) \in \operatorname{im} T$, $\widetilde{T}(v + \ker T) = T(v)$.

    **Uniqueness**: Determined uniquely by $\widetilde{T} \circ \pi = T$. $\blacksquare$

!!! example "Example 13A.3"
    Let $T: \mathbb{R}^3 \to \mathbb{R}^2$ be defined by $T(x, y, z) = (x + y, y + z)$.

    $\ker T = \{(x, y, z) : x + y = 0,\; y + z = 0\} = \{(t, -t, t) : t \in \mathbb{R}\}$, of dimension $1$.

    $\operatorname{im} T = \mathbb{R}^2$ (one can verify $T$ is surjective).

    The first isomorphism theorem gives $\mathbb{R}^3/\ker T \cong \mathbb{R}^2$, and $\dim(\mathbb{R}^3/\ker T) = 3 - 1 = 2 = \dim \mathbb{R}^2$.

!!! theorem "Theorem 13A.6 (Second isomorphism theorem)"
    Let $U, W$ be subspaces of $V$. Then

    $$
    (U + W)/W \cong U/(U \cap W).
    $$

    In particular, $\dim(U + W) - \dim W = \dim U - \dim(U \cap W)$.

??? proof "Proof"
    Define $\varphi: U \to (U + W)/W$ by $\varphi(u) = u + W$. This is a linear map.

    **Surjectivity**: For $(u + w) + W \in (U+W)/W$, we have $(u + w) + W = u + W = \varphi(u)$.

    **Kernel**: $\varphi(u) = W$ if and only if $u \in W$, i.e., $u \in U \cap W$. So $\ker \varphi = U \cap W$.

    By the first isomorphism theorem, $U/(U \cap W) \cong (U + W)/W$. $\blacksquare$

!!! theorem "Theorem 13A.7 (Third isomorphism theorem)"
    Let $W \subseteq U$ be subspaces of $V$. Then

    $$
    (V/W)/(U/W) \cong V/U.
    $$

??? proof "Proof"
    Define $\varphi: V/W \to V/U$ by $\varphi(v + W) = v + U$.

    **Well-definedness**: If $v + W = v' + W$, then $v - v' \in W \subseteq U$, so $v + U = v' + U$.

    **Linearity** and **surjectivity** are immediate.

    **Kernel**: $\varphi(v + W) = U$ if and only if $v \in U$, i.e., $v + W \in U/W$. So $\ker \varphi = U/W$.

    By the first isomorphism theorem, $(V/W)/(U/W) \cong V/U$. $\blacksquare$

---

## 13A.3 Linear Functionals and the Dual Space

<div class="context-flow" markdown>

**Core question**: What is special about linear maps $V \to \mathbb{F}$? → Linear functionals → All linear functionals form the dual space $V^*$

</div>

!!! definition "Definition 13A.4 (Linear functional)"
    Let $V$ be a vector space over $\mathbb{F}$. A linear map $\varphi: V \to \mathbb{F}$ is called a **linear functional** on $V$.

!!! example "Example 13A.4"
    The following are all linear functionals:

    1. $\varphi: \mathbb{R}^n \to \mathbb{R}$, $\varphi(x_1, \ldots, x_n) = a_1 x_1 + \cdots + a_n x_n$ (a linear combination of coordinates).
    2. $\varphi: C[0,1] \to \mathbb{R}$, $\varphi(f) = \int_0^1 f(t)\, dt$ (definite integration).
    3. $\varphi: \mathbb{F}[x]_{\leq n} \to \mathbb{F}$, $\varphi(p) = p(c)$ (evaluation at the point $c$).

!!! definition "Definition 13A.5 (Dual space)"
    The set of all linear functionals on $V$, equipped with pointwise addition and scalar multiplication, forms a vector space

    $$
    V^* = \mathcal{L}(V, \mathbb{F}) = \operatorname{Hom}(V, \mathbb{F}),
    $$

    called the **dual space** (or **algebraic dual**) of $V$.

!!! theorem "Theorem 13A.8 (Dimension of the dual space)"
    If $V$ is finite-dimensional, then

    $$
    \dim V^* = \dim V.
    $$

??? proof "Proof"
    $V^* = \mathcal{L}(V, \mathbb{F})$. By the dimension formula for spaces of linear maps, $\dim \mathcal{L}(V, \mathbb{F}) = \dim V \cdot \dim \mathbb{F} = \dim V \cdot 1 = \dim V$. $\blacksquare$

!!! theorem "Theorem 13A.9 (Kernel of a linear functional)"
    Let $\varphi \in V^*$ be a nonzero linear functional. Then $\ker \varphi$ is a subspace of $V$ of codimension $1$ (a hyperplane):

    $$
    \dim \ker \varphi = \dim V - 1.
    $$

    Conversely, every codimension-$1$ subspace of $V$ is the kernel of some nonzero linear functional.

??? proof "Proof"
    $\varphi \neq 0$ implies $\operatorname{im} \varphi = \mathbb{F}$ (since $\operatorname{im} \varphi$ is a subspace of $\mathbb{F}$, and a nonzero subspace of $\mathbb{F}$ must be $\mathbb{F}$ itself). By the rank-nullity theorem, $\dim \ker \varphi = \dim V - \dim \operatorname{im} \varphi = \dim V - 1$.

    Conversely, let $\dim W = \dim V - 1$. Choose $v_0 \in V \setminus W$, so that $V = W \oplus \operatorname{span}(v_0)$. Define $\varphi(w + cv_0) = c$. Then $\varphi$ is a nonzero linear functional with $\ker \varphi = W$. $\blacksquare$

!!! example "Example 13A.5"
    In $\mathbb{R}^3$, the linear functional $\varphi(x, y, z) = 2x - y + 3z$ has kernel

    $$
    \ker \varphi = \{(x, y, z) : 2x - y + 3z = 0\},
    $$

    a plane through the origin (a $2$-dimensional subspace, codimension $1$).

---

## 13A.4 The Dual Basis

<div class="context-flow" markdown>

**Core question**: How does a basis of $V$ determine a basis of $V^*$? → Definition and construction of the dual basis → Dual basis vectors as coordinate extraction maps

</div>

!!! definition "Definition 13A.6 (Dual basis)"
    Let $\mathcal{B} = \{e_1, e_2, \ldots, e_n\}$ be a basis for the finite-dimensional vector space $V$. The **dual basis** $\mathcal{B}^* = \{e_1^*, e_2^*, \ldots, e_n^*\}$ of $V^*$ is defined by the condition

    $$
    e_i^*(e_j) = \delta_{ij} = \begin{cases} 1 & \text{if } i = j, \\ 0 & \text{if } i \neq j. \end{cases}
    $$

!!! theorem "Theorem 13A.10 (Existence and uniqueness of the dual basis)"
    The dual basis $\mathcal{B}^*$ exists, is unique, and forms a basis for $V^*$.

??? proof "Proof"
    **Existence and uniqueness**: For each $i$, $e_i^*$ is uniquely determined by its values on the basis $\{e_1, \ldots, e_n\}$. Defining $e_i^*(e_j) = \delta_{ij}$ and extending linearly gives a well-defined linear functional.

    **It is a basis**: Suppose $\sum c_i e_i^* = 0$ (the zero functional). Evaluating at $e_j$: $\sum c_i e_i^*(e_j) = c_j = 0$. So $e_1^*, \ldots, e_n^*$ are linearly independent. Since $\dim V^* = n$, they form a basis. $\blacksquare$

!!! theorem "Theorem 13A.11 (Coordinate interpretation of the dual basis)"
    Let $v = \sum_{i=1}^n a_i e_i \in V$. Then

    $$
    e_i^*(v) = a_i.
    $$

    That is, $e_i^*$ extracts the $i$-th coordinate of $v$ with respect to the basis $\mathcal{B}$.

    Furthermore, any $\varphi \in V^*$ can be written as $\varphi = \sum_{i=1}^n \varphi(e_i)\, e_i^*$.

??? proof "Proof"
    $e_i^*(v) = e_i^*\!\left(\sum_j a_j e_j\right) = \sum_j a_j\, e_i^*(e_j) = \sum_j a_j\, \delta_{ij} = a_i$.

    For the second claim, let $\psi = \sum_i \varphi(e_i)\, e_i^*$. On basis vectors: $\psi(e_j) = \sum_i \varphi(e_i)\, \delta_{ij} = \varphi(e_j)$. Hence $\psi = \varphi$. $\blacksquare$

!!! example "Example 13A.6"
    $V = \mathbb{R}^3$ with the standard basis $\{e_1, e_2, e_3\}$. The dual basis is

    $$
    e_1^*(x, y, z) = x, \quad e_2^*(x, y, z) = y, \quad e_3^*(x, y, z) = z.
    $$

    These are the three coordinate functions.

!!! example "Example 13A.7"
    $V = \mathbb{R}^2$ with basis $\{v_1, v_2\} = \{(1, 1), (1, -1)\}$. Find the dual basis.

    Let $v_1^*(x, y) = ax + by$. From $v_1^*(1, 1) = 1$ and $v_1^*(1, -1) = 0$: $a + b = 1$, $a - b = 0$, giving $a = b = \frac{1}{2}$.

    Let $v_2^*(x, y) = cx + dy$. From $v_2^*(1, 1) = 0$ and $v_2^*(1, -1) = 1$: $c + d = 0$, $c - d = 1$, giving $c = \frac{1}{2}$, $d = -\frac{1}{2}$.

    Therefore $v_1^*(x, y) = \frac{x + y}{2}$ and $v_2^*(x, y) = \frac{x - y}{2}$.

---

## 13A.5 Annihilators

<div class="context-flow" markdown>

**Core question**: Which functionals in $V^*$ vanish on a given subspace $U$? → The annihilator $U^0$ → Dimension formula for annihilators → Annihilators establish a duality between subspaces of $V$ and subspaces of $V^*$

</div>

!!! definition "Definition 13A.7 (Annihilator)"
    Let $U$ be a subspace of $V$. The **annihilator** of $U$ is

    $$
    U^0 = \{\varphi \in V^* : \varphi(u) = 0 \text{ for all } u \in U\}.
    $$

    $U^0$ is a subspace of $V^*$.

!!! theorem "Theorem 13A.12 (Dimension of the annihilator)"
    Let $V$ be finite-dimensional and $U$ a subspace of $V$. Then

    $$
    \dim U^0 = \dim V - \dim U.
    $$

??? proof "Proof"
    Choose a basis $\{e_1, \ldots, e_k\}$ for $U$ and extend it to a basis $\{e_1, \ldots, e_k, e_{k+1}, \ldots, e_n\}$ for $V$. Let $\{e_1^*, \ldots, e_n^*\}$ be the dual basis.

    Claim: $U^0 = \operatorname{span}\{e_{k+1}^*, \ldots, e_n^*\}$.

    On one hand, for $i > k$ and $j \leq k$: $e_i^*(e_j) = 0$, so $e_i^* \in U^0$.

    On the other hand, let $\varphi = \sum_{i=1}^n c_i e_i^* \in U^0$. For $j \leq k$: $\varphi(e_j) = c_j = 0$. So $\varphi = \sum_{i=k+1}^n c_i e_i^*$.

    Hence $U^0 = \operatorname{span}\{e_{k+1}^*, \ldots, e_n^*\}$ and $\dim U^0 = n - k = \dim V - \dim U$. $\blacksquare$

!!! theorem "Theorem 13A.13 (Duality properties of annihilators)"
    Let $V$ be finite-dimensional and $U, W$ subspaces of $V$. Then:

    1. $U \subseteq W \Leftrightarrow W^0 \subseteq U^0$.
    2. $(U + W)^0 = U^0 \cap W^0$.
    3. $(U \cap W)^0 = U^0 + W^0$.
    4. $(U^0)^0 = U$ (under the canonical isomorphism $V \cong V^{**}$).

??? proof "Proof"
    1. If $U \subseteq W$ and $\varphi \in W^0$, then $\varphi$ vanishes on $W$, hence on $U$, so $\varphi \in U^0$. The converse follows by dimension counting.

    2. $\varphi \in (U + W)^0$ iff $\varphi(u + w) = 0$ for all $u \in U, w \in W$, iff $\varphi(u) = 0$ for all $u \in U$ and $\varphi(w) = 0$ for all $w \in W$, iff $\varphi \in U^0 \cap W^0$.

    3. By (2), $(U^0 \cap W^0)^0 = ((U+W)^0)^0 = U + W$ (using (4)). Taking annihilators again and using (4), or by a direct dimension count:

    $$
    \dim(U^0 + W^0) = \dim U^0 + \dim W^0 - \dim(U^0 \cap W^0).
    $$

    By (2), $\dim(U^0 \cap W^0) = \dim(U + W)^0 = n - \dim(U+W)$. Substituting and simplifying gives $\dim(U^0 + W^0) = n - \dim(U \cap W) = \dim(U \cap W)^0$. The inclusion relation then yields equality.

    4. See the discussion of the double dual in Section 13A.7. $\blacksquare$

!!! example "Example 13A.8"
    $V = \mathbb{R}^3$, $U = \operatorname{span}\{(1, 0, 1), (0, 1, 1)\}$. Find $U^0$.

    Let $\varphi(x, y, z) = ax + by + cz \in U^0$. From $\varphi(1, 0, 1) = a + c = 0$ and $\varphi(0, 1, 1) = b + c = 0$: $a = -c$, $b = -c$. Taking $c = 1$: $\varphi(x, y, z) = -x - y + z$.

    $$
    U^0 = \operatorname{span}\{(-1, -1, 1)\}.
    $$

    Check: $\dim U^0 = 1 = 3 - 2 = \dim V - \dim U$.

!!! example "Example 13A.9"
    Quotient spaces and annihilators: the dual space $(V/U)^*$ is naturally isomorphic to $U^0$.

    Define $\Phi: U^0 \to (V/U)^*$ by $\Phi(\varphi)(v + U) = \varphi(v)$. This is well-defined ($\varphi \in U^0$ ensures independence of the representative) and is an isomorphism. This follows from the equality of dimensions ($\dim U^0 = \dim V - \dim U = \dim(V/U) = \dim(V/U)^*$) together with injectivity.

---

## 13A.6 The Transpose (Dual) Map

<div class="context-flow" markdown>

**Core question**: How does a linear map $T: V \to W$ induce a map from $W^*$ to $V^*$? → The transpose map $T^*$ → Transpose maps reverse direction → Relationship to matrix transpose

</div>

!!! definition "Definition 13A.8 (Transpose map)"
    Let $T: V \to W$ be a linear map. The **transpose** (or **dual map**) $T^*: W^* \to V^*$ is defined by

    $$
    T^*(\psi) = \psi \circ T, \quad \psi \in W^*.
    $$

    That is, for $v \in V$, $(T^*\psi)(v) = \psi(T(v))$.

!!! note "Note"
    Notice the reversal of direction: $T: V \to W$ induces $T^*: W^* \to V^*$ (a contravariant functor).

!!! theorem "Theorem 13A.14 (Properties of the transpose map)"
    Let $S, T: V \to W$ and $R: W \to U$ be linear maps, and $\lambda \in \mathbb{F}$. Then:

    1. $T^*$ is a linear map.
    2. $(S + T)^* = S^* + T^*$.
    3. $(\lambda T)^* = \lambda T^*$.
    4. $(R \circ T)^* = T^* \circ R^*$ (note the reversal of order).
    5. $(\operatorname{id}_V)^* = \operatorname{id}_{V^*}$.

??? proof "Proof"
    1. $T^*(\psi_1 + \psi_2) = (\psi_1 + \psi_2) \circ T = \psi_1 \circ T + \psi_2 \circ T = T^*\psi_1 + T^*\psi_2$. Scalar multiplication is similar.

    4. For $\varphi \in U^*$: $(R \circ T)^*(\varphi) = \varphi \circ (R \circ T) = (\varphi \circ R) \circ T = T^*(R^*(\varphi)) = (T^* \circ R^*)(\varphi)$.

    The remaining parts are similar. $\blacksquare$

!!! theorem "Theorem 13A.15 (Kernel and image of the transpose)"
    Let $T: V \to W$ be a linear map between finite-dimensional spaces. Then:

    1. $\ker T^* = (\operatorname{im} T)^0$.
    2. $\operatorname{im} T^* = (\ker T)^0$.
    3. $\operatorname{rank} T^* = \operatorname{rank} T$.

??? proof "Proof"
    1. $\psi \in \ker T^*$ iff $T^*\psi = 0$, i.e., $\psi \circ T = 0$, i.e., $\psi(T(v)) = 0$ for all $v \in V$, i.e., $\psi$ vanishes on $\operatorname{im} T$, i.e., $\psi \in (\operatorname{im} T)^0$.

    2. First, $\operatorname{im} T^* \subseteq (\ker T)^0$: if $\varphi = T^*\psi$, then for $v \in \ker T$, $\varphi(v) = \psi(T(v)) = \psi(0) = 0$, so $\varphi \in (\ker T)^0$.

    Dimension count: $\dim \operatorname{im} T^* = \dim W^* - \dim \ker T^* = \dim W - \dim(\operatorname{im} T)^0 = \dim W - (\dim W - \dim \operatorname{im} T) = \operatorname{rank} T$. And $\dim(\ker T)^0 = \dim V - \dim \ker T = \operatorname{rank} T$. Since the dimensions are equal and we have the inclusion, equality holds.

    3. From (2), $\operatorname{rank} T^* = \dim \operatorname{im} T^* = \operatorname{rank} T$. $\blacksquare$

!!! theorem "Theorem 13A.16 (Transpose map and matrix transpose)"
    Let $\mathcal{B}$ be a basis for $V$, $\mathcal{C}$ a basis for $W$, and $\mathcal{B}^*, \mathcal{C}^*$ the corresponding dual bases. If the matrix of $T: V \to W$ with respect to $\mathcal{B}, \mathcal{C}$ is $A$, then the matrix of $T^*: W^* \to V^*$ with respect to $\mathcal{C}^*, \mathcal{B}^*$ is $A^T$ (the transpose of $A$).

??? proof "Proof"
    Let $A = (a_{ij})$, so that $T(e_j) = \sum_i a_{ij} f_i$. Then

    $$
    (T^* f_i^*)(e_j) = f_i^*(T(e_j)) = f_i^*\!\left(\sum_k a_{kj} f_k\right) = a_{ij}.
    $$

    On the other hand, if $T^* f_i^* = \sum_j b_{ji}\, e_j^*$, then $(T^* f_i^*)(e_j) = b_{ji}$.

    Therefore $b_{ji} = a_{ij}$, so the matrix of $T^*$ is $A^T$. $\blacksquare$

!!! example "Example 13A.10"
    Let $T: \mathbb{R}^2 \to \mathbb{R}^3$, $T(x, y) = (x + y, 2x, y)$. With respect to the standard bases, the matrix of $T$ is

    $$
    A = \begin{pmatrix} 1 & 1 \\ 2 & 0 \\ 0 & 1 \end{pmatrix}.
    $$

    The matrix of $T^*: (\mathbb{R}^3)^* \to (\mathbb{R}^2)^*$ with respect to the dual standard bases is

    $$
    A^T = \begin{pmatrix} 1 & 2 & 0 \\ 1 & 0 & 1 \end{pmatrix}.
    $$

    Verification: $T^*(\varphi)(x, y) = \varphi(x + y, 2x, y)$. If $\varphi(a, b, c) = \alpha a + \beta b + \gamma c$, then $T^*\varphi(x, y) = \alpha(x+y) + 2\beta x + \gamma y = (\alpha + 2\beta)x + (\alpha + \gamma)y$, consistent with the matrix multiplication by $A^T$.

---

## 13A.7 The Double Dual and Canonical Isomorphism

<div class="context-flow" markdown>

**Core question**: What is the relationship between $V^{**} = (V^*)^*$ and $V$? → Canonical isomorphism → "Canonical" means independent of the choice of basis

</div>

!!! definition "Definition 13A.9 (Evaluation map)"
    For $v \in V$, define $\hat{v}: V^* \to \mathbb{F}$ by

    $$
    \hat{v}(\varphi) = \varphi(v), \quad \varphi \in V^*.
    $$

    Then $\hat{v}$ is a linear functional on $V^*$, i.e., $\hat{v} \in V^{**}$.

!!! theorem "Theorem 13A.17 (Canonical isomorphism)"
    The map $\iota: V \to V^{**}$ defined by $\iota(v) = \hat{v}$ is an **injective linear map**. When $V$ is finite-dimensional, $\iota$ is an **isomorphism**, called the **canonical isomorphism**.

??? proof "Proof"
    **Linearity**: $\widehat{u + v}(\varphi) = \varphi(u + v) = \varphi(u) + \varphi(v) = \hat{u}(\varphi) + \hat{v}(\varphi)$, so $\widehat{u+v} = \hat{u} + \hat{v}$. Scalar multiplication is similar.

    **Injectivity**: Suppose $\iota(v) = 0$, i.e., $\hat{v} = 0$. Then $\varphi(v) = 0$ for all $\varphi \in V^*$. If $v \neq 0$, extend $v$ to a basis $\{v, e_2, \ldots, e_n\}$ of $V$ and let $v^*$ be the first dual basis vector. Then $v^*(v) = 1 \neq 0$, a contradiction. Hence $v = 0$.

    **Isomorphism in finite dimensions**: $\dim V^{**} = \dim V^* = \dim V$. An injective linear map between spaces of equal dimension is automatically an isomorphism. $\blacksquare$

!!! note "Note"
    The word "canonical" has a precise mathematical meaning (natural transformation in category theory): the definition of $\iota$ does not depend on a choice of basis and is compatible with linear maps. Concretely, if $T: V \to W$, then $T^{**} \circ \iota_V = \iota_W \circ T$ (a commutative diagram).

    By contrast, although $V$ and $V^*$ are isomorphic in finite dimensions (having equal dimension), there is no canonical isomorphism between them — any isomorphism requires a choice of basis or the introduction of an inner product.

!!! theorem "Theorem 13A.18 (Dual of the dual basis)"
    Let $\{e_1, \ldots, e_n\}$ be a basis for $V$ and $\{e_1^*, \ldots, e_n^*\}$ the dual basis. Under the canonical isomorphism $\iota: V \to V^{**}$,

    $$
    \iota(e_i) = e_i^{**},
    $$

    where $\{e_1^{**}, \ldots, e_n^{**}\}$ is the dual basis of $\{e_1^*, \ldots, e_n^*\}$. That is, the basis of $V$ maps under the canonical isomorphism to the dual of the dual basis.

??? proof "Proof"
    We need to verify $\hat{e}_i(e_j^*) = \delta_{ij}$. By definition, $\hat{e}_i(e_j^*) = e_j^*(e_i) = \delta_{ji} = \delta_{ij}$. Therefore $\hat{e}_i$ is precisely the $i$-th element of the dual basis of $\{e_1^*, \ldots, e_n^*\}$. $\blacksquare$

!!! example "Example 13A.11"
    Application of the canonical isomorphism: proving annihilator property (4) — $(U^0)^0 = U$.

    Under the canonical isomorphism $\iota: V \to V^{**}$, $(U^0)^0 = \{\Phi \in V^{**} : \Phi(\varphi) = 0 \text{ for all } \varphi \in U^0\}$.

    If $v \in U$, then for any $\varphi \in U^0$, $\hat{v}(\varphi) = \varphi(v) = 0$, so $\iota(U) \subseteq (U^0)^0$.

    Dimension count: $\dim(U^0)^0 = \dim V^* - \dim U^0 = \dim V - (\dim V - \dim U) = \dim U = \dim \iota(U)$.

    Hence $\iota(U) = (U^0)^0$, i.e., under the canonical isomorphism, $(U^0)^0$ corresponds to $U$.

---

## 13A.8 Applications in Finite Dimensions

<div class="context-flow" markdown>

**Core question**: How does dual space theory apply to concrete linear algebra problems? → Dual description of linear systems → Dual perspective on change of basis → Dual characterization of rank

</div>

### Dual description of linear systems

!!! theorem "Theorem 13A.19 (Solution space and annihilators)"
    The solution space $S$ of the homogeneous system $A\mathbf{x} = \mathbf{0}$ ($A$ an $m \times n$ matrix) satisfies

    $$
    S = \{v \in \mathbb{F}^n : \varphi_i(v) = 0, \; i = 1, \ldots, m\},
    $$

    where $\varphi_i \in (\mathbb{F}^n)^*$ is the linear functional defined by the $i$-th row of $A$. Setting $W = \operatorname{span}\{\varphi_1, \ldots, \varphi_m\} \subseteq (\mathbb{F}^n)^*$, we have $S = W^0$ (viewing $W$ as a subspace of $(\mathbb{F}^n)^*$ and $W^0$ as the corresponding subspace of $\mathbb{F}^n$ under the canonical isomorphism).

    Therefore $\dim S = n - \dim W = n - \operatorname{rank} A$.

??? proof "Proof"
    $v \in S$ iff $Av = 0$, i.e., the inner product of each row of $A$ with $v$ is zero, i.e., $\varphi_i(v) = 0$ for all $i$. Vanishing on all $\varphi_i$ is equivalent to vanishing on $\operatorname{span}\{\varphi_i\}$.

    Through the canonical isomorphism $\iota: \mathbb{F}^n \to (\mathbb{F}^n)^{**}$, $S$ is precisely the annihilator of $W$ in $(\mathbb{F}^n)^{**}$ pulled back to $\mathbb{F}^n$ via $\iota$, i.e., $S = W^0$ (in the generalized sense). Hence $\dim S = n - \dim W = n - \operatorname{rank} A$. $\blacksquare$

### Dual perspective on change of basis

!!! theorem "Theorem 13A.20 (Change-of-basis formula for dual bases)"
    Let $\mathcal{B} = \{e_1, \ldots, e_n\}$ and $\mathcal{B}' = \{e_1', \ldots, e_n'\}$ be two bases for $V$ with transition matrix $P$ (i.e., $e_j' = \sum_i p_{ij} e_i$). Then the transition matrix for the dual bases is $(P^{-1})^T$:

    $$
    e_j'^* = \sum_i \left[(P^{-1})^T\right]_{ij} e_i^* = \sum_i (P^{-1})_{ji}\, e_i^*.
    $$

??? proof "Proof"
    Let $e_j'^* = \sum_i q_{ij}\, e_i^*$. From $e_j'^*(e_k') = \delta_{jk}$:

    $$
    \delta_{jk} = e_j'^*(e_k') = \sum_i q_{ij}\, e_i^*\!\left(\sum_l p_{lk}\, e_l\right) = \sum_i q_{ij}\, p_{ik} = (Q^T P)_{jk}.
    $$

    Hence $Q^T P = I$, i.e., $Q^T = P^{-1}$, so $Q = (P^{-1})^T$. $\blacksquare$

### Dual characterization of rank

!!! theorem "Theorem 13A.21 (Dual proof that row rank equals column rank)"
    The row rank of an $m \times n$ matrix $A$ equals its column rank.

??? proof "Proof"
    View $A$ as a linear map $T: \mathbb{F}^n \to \mathbb{F}^m$. The column rank is $\dim \operatorname{im} T = \operatorname{rank} T$.

    The row rank is $\dim \operatorname{span}\{\text{rows of } A\}$. The rows of $A$ define linear functionals $\varphi_1, \ldots, \varphi_m$ in $(\mathbb{F}^n)^*$. The image of $A^T: (\mathbb{F}^m)^* \to (\mathbb{F}^n)^*$ is precisely $\operatorname{span}\{\varphi_1, \ldots, \varphi_m\}$ (since $T^*(f_i^*) = f_i^* \circ T = \varphi_i$, where $\{f_1^*, \ldots, f_m^*\}$ is the standard basis of $(\mathbb{F}^m)^*$). Hence the row rank equals $\operatorname{rank} T^*$.

    By Theorem 13A.15 (3), $\operatorname{rank} T^* = \operatorname{rank} T$. $\blacksquare$

!!! example "Example 13A.12"
    Dual theory in Lagrange interpolation.

    Let $V = \mathbb{F}[x]_{\leq n}$ (the space of polynomials of degree $\leq n$), and choose $n + 1$ distinct points $c_0, c_1, \ldots, c_n \in \mathbb{F}$. Define the evaluation functionals $\varepsilon_i: V \to \mathbb{F}$ by $\varepsilon_i(p) = p(c_i)$.

    $\{\varepsilon_0, \varepsilon_1, \ldots, \varepsilon_n\}$ is a basis for $V^*$ (since if $\sum a_i \varepsilon_i = 0$, then the polynomial $\sum a_i \ell_i(x) = 0$, where $\ell_i$ are the Lagrange basis polynomials, so $a_i = 0$).

    The dual basis of this set (under $V^{**} \cong V$) is precisely the **Lagrange basis polynomials** $\ell_0, \ell_1, \ldots, \ell_n$, since $\varepsilon_i(\ell_j) = \ell_j(c_i) = \delta_{ij}$.

    This is the dual-space interpretation of the polynomial interpolation from Chapter 0.

!!! example "Example 13A.13"
    Let $T: V \to V$ be a linear operator ($V$ finite-dimensional) with $T^2 = T$ ($T$ is idempotent). Prove $(T^*)^2 = T^*$.

    By Theorem 13A.14 (4): $(T^*)^2 = T^* \circ T^* = (T \circ T)^* = (T^2)^* = T^*$.

    Furthermore, $V = \ker T \oplus \operatorname{im} T$, and in the dual space $V^* = \ker T^* \oplus \operatorname{im} T^*$. By Theorem 13A.15:

    - $\ker T^* = (\operatorname{im} T)^0$;
    - $\operatorname{im} T^* = (\ker T)^0$.

    This demonstrates the perfect symmetry of the idempotent decomposition in the dual space.
