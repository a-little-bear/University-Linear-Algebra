# Chapter 4  Vector Spaces

<div class="context-flow" markdown>

**Prerequisites**: Chapter 1 solution set structure · Chapter 2 properties of $\mathbb{R}^{m \times n}$ · **Chapter arc**: 8 axioms → subspaces → linear combinations/span → linear independence → basis and dimension → coordinates/transition matrices → four fundamental subspaces → rank-nullity theorem
**One-line essence**: Abstract from $\mathbb{R}^n$ to an axiomatic framework — "dimension" becomes the sole essential invariant of a vector space (→ Chapter 5 isomorphism theorem)

</div>

A vector space is the core abstract concept of linear algebra. In the preceding chapters, we mainly discussed problems within the concrete framework of $\mathbb{R}^n$. This chapter introduces the axiomatic definition of vector spaces, making linear algebra theory applicable to objects far more general than $\mathbb{R}^n$ — matrices, polynomials, functions, etc. We will systematically discuss subspaces, linear independence, basis and dimension, establish the rank-nullity theorem, and introduce the four fundamental subspaces.

---

## 4.1 Definition and Axioms of Vector Spaces

<div class="context-flow" markdown>

**Axiomatic leap**: The 8 algebraic properties satisfied by $\mathbb{R}^{m \times n}$ in Chapter 2 → abstracted as axioms on any set → "vectors" are no longer limited to arrows or column vectors

</div>

!!! definition "Definition 4.1 (Vector space)"
    Let $\mathbb{F}$ be a field (usually $\mathbb{R}$ or $\mathbb{C}$), and $V$ a nonempty set equipped with **addition** $+: V \times V \to V$ and **scalar multiplication** $\cdot: \mathbb{F} \times V \to V$. If the following 8 axioms hold for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and all $c, d \in \mathbb{F}$, then $V$ is called a **vector space** over $\mathbb{F}$, elements of $V$ are called **vectors**, and elements of $\mathbb{F}$ are called **scalars**:

    **Addition axioms**:

    1. $\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}$ (commutativity)
    2. $(\mathbf{u} + \mathbf{v}) + \mathbf{w} = \mathbf{u} + (\mathbf{v} + \mathbf{w})$ (associativity)
    3. There exists a **zero vector** $\mathbf{0} \in V$ such that $\mathbf{v} + \mathbf{0} = \mathbf{v}$ (additive identity)
    4. For each $\mathbf{v} \in V$, there exists $-\mathbf{v} \in V$ such that $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$ (additive inverse)

    **Scalar multiplication axioms**:

    5. $c(\mathbf{u} + \mathbf{v}) = c\mathbf{u} + c\mathbf{v}$ (distributivity over vector addition)
    6. $(c + d)\mathbf{v} = c\mathbf{v} + d\mathbf{v}$ (distributivity over scalar addition)
    7. $c(d\mathbf{v}) = (cd)\mathbf{v}$ (associativity of scalar multiplication)
    8. $1\mathbf{v} = \mathbf{v}$ (scalar multiplicative identity)

!!! proposition "Proposition 4.1 (Basic consequences of vector space axioms)"
    In a vector space $V$:

    1. The zero vector is unique.
    2. The additive inverse of each vector is unique.
    3. $0\mathbf{v} = \mathbf{0}$ (scalar $0$ times any vector gives the zero vector).
    4. $c\mathbf{0} = \mathbf{0}$ (any scalar times the zero vector gives the zero vector).
    5. $(-1)\mathbf{v} = -\mathbf{v}$.

??? proof "Proof"
    3. $0\mathbf{v} = (0+0)\mathbf{v} = 0\mathbf{v} + 0\mathbf{v}$. Adding the additive inverse of $0\mathbf{v}$ to both sides gives $\mathbf{0} = 0\mathbf{v}$.

    5. $\mathbf{v} + (-1)\mathbf{v} = 1\mathbf{v} + (-1)\mathbf{v} = (1 + (-1))\mathbf{v} = 0\mathbf{v} = \mathbf{0}$. Therefore $(-1)\mathbf{v}$ is the additive inverse of $\mathbf{v}$, and by uniqueness, $(-1)\mathbf{v} = -\mathbf{v}$. $\blacksquare$

---

## 4.2 Common Examples of Vector Spaces

<div class="context-flow" markdown>

**The power of axioms**: $\mathbb{R}^n$, matrix space $\mathbb{R}^{m \times n}$, polynomial space $\mathbb{P}_n$, function space $C(\mathbb{R})$ → diverse in form but sharing the same algebraic structure

</div>

!!! example "Example 4.1"
    **$\mathbb{R}^n$ space**: The set of $n$-dimensional real column vectors $\mathbb{R}^n$, with componentwise addition and scalar multiplication, is the most basic vector space.

!!! example "Example 4.2"
    **Matrix space $\mathbb{R}^{m \times n}$**: The set of all $m \times n$ real matrices forms a vector space under matrix addition and scalar multiplication.

!!! example "Example 4.3"
    **Polynomial space $\mathbb{P}_n$**: The set of all real-coefficient polynomials of degree at most $n$ forms a vector space, with the usual polynomial addition and scalar multiplication. The zero vector is the zero polynomial.

    The set of all real-coefficient polynomials is denoted $\mathbb{P}$, which also forms a vector space (infinite-dimensional).

!!! example "Example 4.4"
    **Function space $\mathcal{F}(\mathbb{R}, \mathbb{R})$**: The set of all functions from $\mathbb{R}$ to $\mathbb{R}$, with pointwise addition $(f+g)(x) = f(x) + g(x)$ and scalar multiplication $(cf)(x) = cf(x)$, forms a vector space. The continuous function space $C(\mathbb{R})$, the differentiable function space $C^1(\mathbb{R})$, etc., are subspaces.

!!! example "Example 4.5"
    **Zero space $\{0\}$**: The set containing only the zero vector forms a vector space, called the **zero space** or **trivial vector space**, with dimension $0$.

---

## 4.3 Subspaces

<div class="context-flow" markdown>

**"Spaces within spaces"**: Closed under addition and scalar multiplication = subspace → the homogeneous solution set from Chapter 1 is the prototype → must pass through the origin (contain $\mathbf{0}$)

</div>

!!! definition "Definition 4.2 (Subspace)"
    Let $V$ be a vector space and $W \subseteq V$ a nonempty subset. If $W$ is itself a vector space under the addition and scalar multiplication of $V$, then $W$ is called a **subspace** of $V$.

!!! theorem "Theorem 4.1 (Subspace criterion)"
    Let $V$ be a vector space and $W \subseteq V$ a nonempty subset. Then $W$ is a subspace of $V$ if and only if:

    1. **Closed under addition**: If $\mathbf{u}, \mathbf{v} \in W$, then $\mathbf{u} + \mathbf{v} \in W$.
    2. **Closed under scalar multiplication**: If $\mathbf{v} \in W$, $c \in \mathbb{F}$, then $c\mathbf{v} \in W$.

    Equivalently, $W$ is a subspace if and only if: $W$ is nonempty, and for any $\mathbf{u}, \mathbf{v} \in W$, $c, d \in \mathbb{F}$, we have $c\mathbf{u} + d\mathbf{v} \in W$ (closed under linear combinations).

??? proof "Proof"
    $(\Rightarrow)$ A vector space is obviously closed under addition and scalar multiplication.

    $(\Leftarrow)$ We need to verify 8 axioms. Commutativity, associativity, distributivity, etc., are inherited from $V$. Taking $c = 0$ in condition 2 gives $\mathbf{0} = 0\mathbf{v} \in W$ (additive identity exists). Taking $c = -1$ gives $-\mathbf{v} \in W$ (additive inverse exists). $\blacksquare$

!!! example "Example 4.6"
    The following are subspaces of $\mathbb{R}^3$:

    - $\{\mathbf{0}\}$ (zero subspace).
    - A line through the origin $\{t\mathbf{v} : t \in \mathbb{R}\}$, where $\mathbf{v} \neq \mathbf{0}$.
    - A plane through the origin $\{s\mathbf{u} + t\mathbf{v} : s, t \in \mathbb{R}\}$, where $\mathbf{u}, \mathbf{v}$ are linearly independent.
    - $\mathbb{R}^3$ itself.

    Lines or planes not passing through the origin are **not** subspaces (they do not contain the zero vector).

---

## 4.4 Linear Combinations and Span

<div class="context-flow" markdown>

**Generating subspaces**: Given a set of vectors → the set of all linear combinations = the **smallest subspace** containing these vectors → $A\mathbf{x}=\mathbf{b}$ has a solution $\Leftrightarrow$ $\mathbf{b} \in \operatorname{Span}\{\text{column vectors}\}$

</div>

!!! definition "Definition 4.3 (Linear combination)"
    Let $V$ be a vector space, $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k \in V$. A vector of the form

    $$
    c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_k\mathbf{v}_k \quad (c_i \in \mathbb{F})
    $$

    is called a **linear combination** of $\mathbf{v}_1, \ldots, \mathbf{v}_k$.

!!! definition "Definition 4.4 (Span)"
    The set of all linear combinations of the vectors $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\}$ is called the **span** of this set, denoted

    $$
    \operatorname{Span}\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\} = \left\{\sum_{i=1}^k c_i\mathbf{v}_i : c_1, \ldots, c_k \in \mathbb{F}\right\}.
    $$

    If $\operatorname{Span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\} = V$, then $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ is said to **span** (or **generate**) $V$.

!!! theorem "Theorem 4.2 (The span is a subspace)"
    $\operatorname{Span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ is a subspace of $V$, and it is the smallest subspace containing $\mathbf{v}_1, \ldots, \mathbf{v}_k$.

??? proof "Proof"
    Let $W = \operatorname{Span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$. Taking all coefficients to be zero gives the zero vector, so $W$ is nonempty.

    Let $\mathbf{u} = \sum a_i\mathbf{v}_i \in W$, $\mathbf{w} = \sum b_i\mathbf{v}_i \in W$, $c, d \in \mathbb{F}$. Then

    $$
    c\mathbf{u} + d\mathbf{w} = \sum(ca_i + db_i)\mathbf{v}_i \in W.
    $$

    So $W$ is closed under linear combinations and is a subspace. Any subspace containing $\mathbf{v}_1, \ldots, \mathbf{v}_k$ must contain all their linear combinations, so $W$ is the smallest such subspace. $\blacksquare$

---

## 4.5 Linear Independence and Linear Dependence

<div class="context-flow" markdown>

**No redundancy = linear independence**: $\sum c_i \mathbf{v}_i = \mathbf{0}$ has only the trivial solution → no vector is a linear combination of the others → in $\mathbb{R}^n$ equivalent to the matrix having full column rank (Chapter 2)

</div>

!!! definition "Definition 4.5 (Linear independence and linear dependence)"
    The set $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\}$ is called **linearly independent** if the equation

    $$
    c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_k\mathbf{v}_k = \mathbf{0}
    $$

    has only the trivial solution $c_1 = c_2 = \cdots = c_k = 0$.

    Otherwise, the set is called **linearly dependent**, i.e., there exist scalars $c_1, \ldots, c_k$, not all zero, such that $\sum c_i\mathbf{v}_i = \mathbf{0}$.

!!! theorem "Theorem 4.3 (Equivalent condition for linear dependence)"
    The set $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ ($k \ge 2$) is linearly dependent if and only if at least one vector can be expressed as a linear combination of the others.

??? proof "Proof"
    $(\Rightarrow)$ If $\sum c_i\mathbf{v}_i = \mathbf{0}$ with some $c_j \neq 0$, then $\mathbf{v}_j = -\frac{1}{c_j}\sum_{i \neq j} c_i\mathbf{v}_i$.

    $(\Leftarrow)$ If $\mathbf{v}_j = \sum_{i \neq j} a_i\mathbf{v}_i$, then $\sum_{i \neq j} a_i\mathbf{v}_i - \mathbf{v}_j = \mathbf{0}$, where the coefficient of $\mathbf{v}_j$ is $-1 \neq 0$. $\blacksquare$

!!! proposition "Proposition 4.2 (Properties of linear independence)"
    1. A single nonzero vector $\{\mathbf{v}\}$ ($\mathbf{v} \neq \mathbf{0}$) is linearly independent.
    2. A set containing the zero vector is always linearly dependent.
    3. Any subset of a linearly independent set is linearly independent.
    4. If $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ is linearly independent and $\mathbf{v}_{k+1} \notin \operatorname{Span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$, then $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{v}_{k+1}\}$ is also linearly independent.

!!! example "Example 4.7"
    In $\mathbb{R}^3$, determine whether $\{(1,0,1)^T,\; (0,1,1)^T,\; (1,1,0)^T\}$ is linearly independent.

    Form the matrix and find the row echelon form:

    $$
    \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 1 & 1 & 0 \end{pmatrix}
    \xrightarrow{R_3 - R_1}
    \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 0 & 1 & -1 \end{pmatrix}
    \xrightarrow{R_3 - R_2}
    \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & -2 \end{pmatrix}
    $$

    There are $3$ pivots, so the homogeneous equation has only the trivial solution, and the three vectors are linearly independent.

---

## 4.6 Basis and Dimension

<div class="context-flow" markdown>

**Basis = linearly independent + spanning**: Every vector is uniquely expressed as a linear combination of basis vectors → all bases have the same size (well-definedness of dimension) → $\dim(V)$ is the ultimate invariant of a vector space

</div>

!!! definition "Definition 4.6 (Basis)"
    Let $V$ be a vector space. The set $\mathcal{B} = \{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n\}$ is called a **basis** of $V$ if:

    1. $\mathcal{B}$ is linearly independent.
    2. $\operatorname{Span}(\mathcal{B}) = V$ ($\mathcal{B}$ spans $V$).

!!! example "Example 4.8"
    The **standard basis** of $\mathbb{R}^n$ is $\{\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\}$, where $\mathbf{e}_i$ has $1$ in the $i$-th component and $0$ elsewhere.

    The standard basis of $\mathbb{P}_2$ is $\{1, x, x^2\}$.

    The standard basis of $\mathbb{R}^{2 \times 2}$ is $\left\{\begin{pmatrix}1&0\\0&0\end{pmatrix}, \begin{pmatrix}0&1\\0&0\end{pmatrix}, \begin{pmatrix}0&0\\1&0\end{pmatrix}, \begin{pmatrix}0&0\\0&1\end{pmatrix}\right\}$.

!!! theorem "Theorem 4.4 (Unique representation property of a basis)"
    $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ is a basis of $V$ if and only if every vector in $V$ can be uniquely expressed as a linear combination of vectors in $\mathcal{B}$.

??? proof "Proof"
    $(\Rightarrow)$ **Existence** is guaranteed by the spanning property. **Uniqueness**: if $\mathbf{v} = \sum a_i\mathbf{v}_i = \sum b_i\mathbf{v}_i$, then $\sum(a_i - b_i)\mathbf{v}_i = \mathbf{0}$. By linear independence, $a_i - b_i = 0$, i.e., $a_i = b_i$.

    $(\Leftarrow)$ Unique representation implies spanning (every vector has at least one representation) and linear independence (the unique representation of $\mathbf{0}$ has all-zero coefficients). $\blacksquare$

!!! theorem "Theorem 4.5 (Well-definedness of dimension)"
    Any two bases of a finite-dimensional vector space $V$ contain the same number of vectors.

??? proof "Proof"
    Let $\mathcal{B}_1 = \{\mathbf{u}_1, \ldots, \mathbf{u}_m\}$ and $\mathcal{B}_2 = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ both be bases of $V$.

    Since $\mathcal{B}_1$ spans $V$ and $\mathcal{B}_2$ is linearly independent, by the Steinitz replacement lemma, $n \le m$.

    Symmetrically, since $\mathcal{B}_2$ spans $V$ and $\mathcal{B}_1$ is linearly independent, $m \le n$.

    Therefore $m = n$. $\blacksquare$

!!! definition "Definition 4.7 (Dimension)"
    The number of vectors in a basis of a finite-dimensional vector space $V$ is called the **dimension** of $V$, denoted $\dim(V)$. The dimension of the zero space $\{\mathbf{0}\}$ is defined to be $0$.

!!! theorem "Theorem 4.6 (Basis criteria)"
    Let $\dim(V) = n$. Then:

    1. Any set of more than $n$ vectors in $V$ is linearly dependent.
    2. Any set of fewer than $n$ vectors cannot span $V$.
    3. Any $n$ linearly independent vectors in $V$ form a basis of $V$.
    4. Any $n$ vectors in $V$ that span $V$ form a basis of $V$.

---

## 4.7 Coordinates and Change of Basis

<div class="context-flow" markdown>

**Bridge from abstract to concrete**: Choosing a basis → each vector corresponds to a unique coordinate vector $[\mathbf{x}]_\mathcal{B} \in \mathbb{R}^n$ → change of basis = transition matrix $P$ → Chapter 5: change of basis $\to$ similar matrices

</div>

!!! definition "Definition 4.8 (Coordinate vector)"
    Let $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ be an ordered basis of vector space $V$, and $\mathbf{x} \in V$ be uniquely expressed as

    $$
    \mathbf{x} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_n\mathbf{v}_n.
    $$

    Then the column vector $[\mathbf{x}]_\mathcal{B} = (c_1, c_2, \ldots, c_n)^T \in \mathbb{R}^n$ is called the **coordinate vector** of $\mathbf{x}$ with respect to basis $\mathcal{B}$.

!!! definition "Definition 4.9 (Transition matrix)"
    Let $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ and $\mathcal{B}' = \{\mathbf{w}_1, \ldots, \mathbf{w}_n\}$ be two bases of $V$. Express each vector in $\mathcal{B}'$ using $\mathcal{B}$:

    $$
    \mathbf{w}_j = \sum_{i=1}^n p_{ij}\mathbf{v}_i, \quad j = 1, \ldots, n.
    $$

    The matrix $P = (p_{ij})$ is called the **transition matrix** (change-of-basis matrix) from basis $\mathcal{B}'$ to basis $\mathcal{B}$, satisfying

    $$
    [\mathbf{x}]_\mathcal{B} = P [\mathbf{x}]_{\mathcal{B}'}.
    $$

!!! theorem "Theorem 4.7 (Transition matrices are invertible)"
    The transition matrix $P$ is invertible, and $P^{-1}$ is the transition matrix from $\mathcal{B}$ to $\mathcal{B}'$.

!!! example "Example 4.9"
    In $\mathbb{R}^2$, let $\mathcal{B} = \{\mathbf{e}_1, \mathbf{e}_2\}$ (standard basis), $\mathcal{B}' = \{(1,1)^T, (1,-1)^T\}$.

    Then $\mathbf{w}_1 = 1\cdot\mathbf{e}_1 + 1\cdot\mathbf{e}_2$, $\mathbf{w}_2 = 1\cdot\mathbf{e}_1 + (-1)\cdot\mathbf{e}_2$, and the transition matrix is

    $$
    P = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}.
    $$

    For vector $\mathbf{x} = (3, 1)^T$, its coordinates in $\mathcal{B}'$ are $[\mathbf{x}]_{\mathcal{B}'} = P^{-1}[\mathbf{x}]_\mathcal{B} = \frac{1}{-2}\begin{pmatrix} -1 & -1 \\ -1 & 1 \end{pmatrix}\begin{pmatrix} 3 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$.

    Verification: $2(1,1)^T + 1(1,-1)^T = (3, 1)^T$.

---

## 4.8 Row Space, Column Space, and Null Space

<div class="context-flow" markdown>

**Four fundamental subspaces of a matrix**: $\operatorname{Col}(A)$, $\operatorname{Row}(A)$, $\operatorname{Null}(A)$, $\operatorname{Null}(A^T)$ → dimensions uniformly controlled by $\operatorname{rank}(A) = r$ → orthogonal complement relations in Chapter 7

</div>

!!! definition "Definition 4.10 (Four fundamental subspaces)"
    Let $A$ be an $m \times n$ matrix. Define the following four subspaces:

    1. **Column space**: $\operatorname{Col}(A) = \{A\mathbf{x} : \mathbf{x} \in \mathbb{R}^n\} \subseteq \mathbb{R}^m$, i.e., the span of the column vectors of $A$.
    2. **Row space**: $\operatorname{Row}(A) = \operatorname{Col}(A^T) \subseteq \mathbb{R}^n$, i.e., the span of the row vectors of $A$.
    3. **Null space / kernel**: $\operatorname{Null}(A) = \{\mathbf{x} \in \mathbb{R}^n : A\mathbf{x} = \mathbf{0}\} \subseteq \mathbb{R}^n$.
    4. **Left null space**: $\operatorname{Null}(A^T) = \{\mathbf{y} \in \mathbb{R}^m : A^T\mathbf{y} = \mathbf{0}\} \subseteq \mathbb{R}^m$.

!!! theorem "Theorem 4.8 (Dimensions of the four subspaces)"
    Let $A$ be an $m \times n$ matrix with $\operatorname{rank}(A) = r$. Then:

    1. $\dim(\operatorname{Col}(A)) = r$
    2. $\dim(\operatorname{Row}(A)) = r$
    3. $\dim(\operatorname{Null}(A)) = n - r$
    4. $\dim(\operatorname{Null}(A^T)) = m - r$

!!! theorem "Theorem 4.9 (Orthogonal complement relations)"
    1. $\operatorname{Row}(A) \oplus \operatorname{Null}(A) = \mathbb{R}^n$, i.e., the row space and null space are orthogonal complements.
    2. $\operatorname{Col}(A) \oplus \operatorname{Null}(A^T) = \mathbb{R}^m$, i.e., the column space and left null space are orthogonal complements.

!!! example "Example 4.10"
    Let $A = \begin{pmatrix} 1 & 2 & 1 & 3 \\ 2 & 4 & 3 & 7 \\ 1 & 2 & 2 & 4 \end{pmatrix}$.

    Reduce to RREF:

    $$
    R = \begin{pmatrix} 1 & 2 & 0 & 2 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 0 \end{pmatrix}
    $$

    $\operatorname{rank}(A) = 2$.

    - **Column space** $\operatorname{Col}(A)$: Pivot columns are columns 1 and 3; basis $\{(1,2,1)^T,\; (1,3,2)^T\}$, dimension $2$.
    - **Row space** $\operatorname{Row}(A)$: The nonzero rows of $R$ form a basis: $\{(1,2,0,2),\; (0,0,1,1)\}$, dimension $2$.
    - **Null space** $\operatorname{Null}(A)$: Free variables $x_2 = s$, $x_4 = t$; basis $\{(-2,1,0,0)^T,\; (-2,0,-1,1)^T\}$, dimension $2$.
    - **Left null space** $\operatorname{Null}(A^T)$: dimension $3 - 2 = 1$.

---

## 4.9 Rank-Nullity Theorem (Dimension Formula)

<div class="context-flow" markdown>

**Conservation of dimension**: $\operatorname{rank}(A) + \operatorname{nullity}(A) = n$ → the "preserved dimensions" plus the "lost dimensions" of a linear map equal the dimension of the domain → this is the matrix version of the rank-nullity theorem from Chapter 5

</div>

!!! theorem "Theorem 4.10 (Rank-nullity theorem)"
    Let $A$ be an $m \times n$ matrix. Then

    $$
    \operatorname{rank}(A) + \operatorname{nullity}(A) = n,
    $$

    where $\operatorname{nullity}(A) = \dim(\operatorname{Null}(A))$ is the **nullity** of $A$.

    More generally, if $T: V \to W$ is a linear transformation between finite-dimensional vector spaces, then

    $$
    \dim(\ker T) + \dim(\operatorname{im} T) = \dim(V).
    $$

??? proof "Proof"
    Let $\dim(V) = n$, $\dim(\ker T) = k$. Take a basis $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ of $\ker T$ and extend it to a basis $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{v}_{k+1}, \ldots, \mathbf{v}_n\}$ of $V$.

    **Claim**: $\{T(\mathbf{v}_{k+1}), \ldots, T(\mathbf{v}_n)\}$ is a basis of $\operatorname{im}(T)$.

    **Spanning**: For any $\mathbf{w} \in \operatorname{im}(T)$, there exists $\mathbf{v} = \sum_{i=1}^n c_i\mathbf{v}_i \in V$ with $T(\mathbf{v}) = \mathbf{w}$. Then

    $$
    \mathbf{w} = T\left(\sum_{i=1}^n c_i\mathbf{v}_i\right) = \sum_{i=1}^k c_i T(\mathbf{v}_i) + \sum_{i=k+1}^n c_i T(\mathbf{v}_i) = \sum_{i=k+1}^n c_i T(\mathbf{v}_i),
    $$

    since $T(\mathbf{v}_i) = \mathbf{0}$ for $i = 1, \ldots, k$.

    **Linear independence**: If $\sum_{i=k+1}^n c_i T(\mathbf{v}_i) = \mathbf{0}$, then $T\left(\sum_{i=k+1}^n c_i\mathbf{v}_i\right) = \mathbf{0}$, so $\sum_{i=k+1}^n c_i\mathbf{v}_i \in \ker T$. Thus there exist scalars $d_1, \ldots, d_k$ with $\sum_{i=k+1}^n c_i\mathbf{v}_i = \sum_{i=1}^k d_i\mathbf{v}_i$, i.e., $\sum_{i=1}^k d_i\mathbf{v}_i - \sum_{i=k+1}^n c_i\mathbf{v}_i = \mathbf{0}$. By linear independence of $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$, all coefficients are zero; in particular $c_{k+1} = \cdots = c_n = 0$.

    Therefore $\dim(\operatorname{im} T) = n - k = \dim(V) - \dim(\ker T)$. $\blacksquare$

!!! example "Example 4.11"
    Let $A$ be a $5 \times 8$ matrix with $\operatorname{rank}(A) = 3$. Then $\operatorname{nullity}(A) = 8 - 3 = 5$. The solution space of $A\mathbf{x} = \mathbf{0}$ is a $5$-dimensional subspace of $\mathbb{R}^8$, and its fundamental system of solutions contains $5$ vectors.

---

## 4.10 Sum and Direct Sum of Subspaces

<div class="context-flow" markdown>

**Algebraic operations on subspaces**: $W_1 + W_2$ is the smallest subspace containing both → when $W_1 \cap W_2 = \{\mathbf{0}\}$ it is a **direct sum** → dimension formula: $\dim(W_1+W_2) = \dim W_1 + \dim W_2 - \dim(W_1 \cap W_2)$

</div>

!!! definition "Definition 4.11 (Sum of subspaces)"
    Let $W_1, W_2$ be subspaces of a vector space $V$. Their **sum** is defined as

    $$
    W_1 + W_2 = \{\mathbf{w}_1 + \mathbf{w}_2 : \mathbf{w}_1 \in W_1, \mathbf{w}_2 \in W_2\}.
    $$

    $W_1 + W_2$ is also a subspace of $V$, and is the smallest subspace containing $W_1 \cup W_2$.

!!! definition "Definition 4.12 (Direct sum)"
    If $W_1 \cap W_2 = \{\mathbf{0}\}$, the sum $W_1 + W_2$ is called a **direct sum**, denoted $W_1 \oplus W_2$. In this case, every vector in $W_1 + W_2$ can be **uniquely** decomposed as $\mathbf{w}_1 + \mathbf{w}_2$ ($\mathbf{w}_i \in W_i$).

!!! theorem "Theorem 4.11 (Equivalent conditions for direct sum)"
    Let $W_1, W_2$ be subspaces of $V$. The following conditions are equivalent:

    1. $W_1 + W_2$ is a direct sum.
    2. $W_1 \cap W_2 = \{\mathbf{0}\}$.
    3. Every vector in $W_1 + W_2$ has a unique decomposition $\mathbf{w}_1 + \mathbf{w}_2$ ($\mathbf{w}_i \in W_i$).
    4. If $\mathbf{w}_1 + \mathbf{w}_2 = \mathbf{0}$ ($\mathbf{w}_i \in W_i$), then $\mathbf{w}_1 = \mathbf{w}_2 = \mathbf{0}$.

??? proof "Proof"
    $(2) \Rightarrow (3)$: If $\mathbf{v} = \mathbf{w}_1 + \mathbf{w}_2 = \mathbf{w}_1' + \mathbf{w}_2'$, then $\mathbf{w}_1 - \mathbf{w}_1' = \mathbf{w}_2' - \mathbf{w}_2 \in W_1 \cap W_2 = \{\mathbf{0}\}$, so $\mathbf{w}_1 = \mathbf{w}_1'$ and $\mathbf{w}_2 = \mathbf{w}_2'$.

    $(3) \Rightarrow (4)$: $\mathbf{0} = \mathbf{0} + \mathbf{0}$ is one decomposition; by uniqueness the result follows.

    $(4) \Rightarrow (2)$: If $\mathbf{v} \in W_1 \cap W_2$, then $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$ ($\mathbf{v} \in W_1$, $-\mathbf{v} \in W_2$), so by condition 4, $\mathbf{v} = \mathbf{0}$. $\blacksquare$

!!! theorem "Theorem 4.12 (Dimension formula)"
    Let $W_1, W_2$ be subspaces of a finite-dimensional vector space $V$. Then

    $$
    \dim(W_1 + W_2) = \dim(W_1) + \dim(W_2) - \dim(W_1 \cap W_2).
    $$

    In particular, if $W_1 + W_2$ is a direct sum, then $\dim(W_1 \oplus W_2) = \dim(W_1) + \dim(W_2)$.

??? proof "Proof"
    Let $\dim(W_1 \cap W_2) = k$, and take a basis $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ of $W_1 \cap W_2$.

    Extend this to a basis $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{u}_1, \ldots, \mathbf{u}_p\}$ of $W_1$, where $p = \dim(W_1) - k$.

    Extend to a basis $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{w}_1, \ldots, \mathbf{w}_q\}$ of $W_2$, where $q = \dim(W_2) - k$.

    One can verify that $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{u}_1, \ldots, \mathbf{u}_p, \mathbf{w}_1, \ldots, \mathbf{w}_q\}$ is a basis of $W_1 + W_2$. Therefore

    $$
    \dim(W_1 + W_2) = k + p + q = \dim(W_1) + \dim(W_2) - \dim(W_1 \cap W_2). \quad \blacksquare
    $$

!!! example "Example 4.12"
    In $\mathbb{R}^4$, let $W_1 = \operatorname{Span}\{(1,0,1,0)^T, (0,1,0,1)^T\}$, $W_2 = \operatorname{Span}\{(1,1,0,0)^T, (0,0,1,1)^T\}$.

    $\dim(W_1) = \dim(W_2) = 2$. To check whether this is a direct sum, find $W_1 \cap W_2$: set $a(1,0,1,0)^T + b(0,1,0,1)^T = c(1,1,0,0)^T + d(0,0,1,1)^T$, giving the system

    $$
    a = c, \quad b = c, \quad a = d, \quad b = d.
    $$

    Therefore $a = b = c = d$, $W_1 \cap W_2 = \operatorname{Span}\{(1,1,1,1)^T\}$, $\dim(W_1 \cap W_2) = 1$.

    $\dim(W_1 + W_2) = 2 + 2 - 1 = 3$, so this is not a direct sum.
