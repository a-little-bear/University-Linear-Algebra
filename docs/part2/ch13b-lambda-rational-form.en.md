# Chapter 13B  Lambda-Matrices and the Rational Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Ch5 Eigenvalues and characteristic polynomials · Ch7 Jordan normal form · **Chapter arc**: $\lambda$-matrices → Elementary operations → Smith normal form → Invariant factors / Elementary divisors → Companion matrix and **Rational canonical form** → Unified derivation of Jordan form → Complete invariants of matrix similarity
Essence: Jordan form depends on an algebraically closed field ($\mathbb{C}$), while the rational canonical form holds over **any field** — invariant factors are "characteristic polynomial hierarchies without factorization"

</div>

The Jordan normal form is the most refined similarity classification tool in matrix theory, but it fundamentally depends on complete factorization of the characteristic polynomial over the base field. When the base field is not algebraically closed (e.g., $\mathbb{R}$ or $\mathbb{Q}$), the Jordan normal form may not exist. The theory of $\lambda$-matrices provides a toolkit that does not rely on eigenvalue decomposition, leading to the **rational canonical form** which exists over any field and gives the ultimate answer to the matrix similarity problem.

---

## 13B.1 Polynomial Matrices

<div class="context-flow" markdown>

**Numerical matrix** $A$ → **Polynomial matrix** $A(\lambda)$: entries extend from $\mathbb{F}$ to $\mathbb{F}[\lambda]$ → The characteristic matrix $\lambda I - A$ is the most important $\lambda$-matrix

</div>

!!! definition "Definition 13B.1 (Lambda-matrix)"
    Let $\mathbb{F}$ be a field. An $m \times n$ matrix whose entries are elements of $\mathbb{F}[\lambda]$ (the ring of polynomials over $\mathbb{F}$),

    $$
    A(\lambda) = (a_{ij}(\lambda))_{m \times n}
    $$

    is called a **$\lambda$-matrix** (polynomial matrix). If the maximum degree among the entries of $A(\lambda)$ is $s$, it can be written as

    $$
    A(\lambda) = A_s\lambda^s + A_{s-1}\lambda^{s-1} + \cdots + A_1\lambda + A_0
    $$

    where $A_i \in \mathbb{F}^{m \times n}$ are numerical matrix coefficients.

!!! definition "Definition 13B.2 (Characteristic matrix)"
    For $A \in \mathbb{F}^{n \times n}$, the matrix

    $$
    \lambda I - A
    $$

    is called the **characteristic matrix** of $A$. It is the most important $\lambda$-matrix; its determinant $\det(\lambda I - A)$ is the characteristic polynomial.

!!! definition "Definition 13B.3 (Rank of a lambda-matrix)"
    The **rank** of a $\lambda$-matrix $A(\lambda)$ is defined as the largest order of a nonzero minor, denoted $\operatorname{rank} A(\lambda)$. If $A(\lambda)$ is $n \times n$ and $\operatorname{rank} A(\lambda) = n$ (i.e., $\det A(\lambda) \not\equiv 0$), then $A(\lambda)$ is called **nonsingular**.

!!! theorem "Theorem 13B.1 (Invertibility of lambda-matrices)"
    An $n \times n$ $\lambda$-matrix $A(\lambda)$ is invertible (i.e., there exists a $\lambda$-matrix $B(\lambda)$ with $A(\lambda)B(\lambda) = I$) if and only if $\det A(\lambda)$ is a nonzero constant.

??? proof "Proof"
    Sufficiency: If $\det A(\lambda) = c \neq 0$ (a constant), then the adjugate matrix $\operatorname{adj}(A(\lambda))$ is also a $\lambda$-matrix, and $A(\lambda)^{-1} = \frac{1}{c}\operatorname{adj}(A(\lambda))$.

    Necessity: If $A(\lambda)B(\lambda) = I$, taking determinants gives $\det A(\lambda) \cdot \det B(\lambda) = 1$. Since $\det A(\lambda)$ and $\det B(\lambda)$ are polynomials in $\mathbb{F}[\lambda]$ whose product is the constant $1$, both must be nonzero constants. $\blacksquare$

!!! example "Example 13B.1"
    The $\lambda$-matrix $A(\lambda) = \begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$ is nonsingular since $\det A(\lambda) = \lambda^2 \neq 0$. However, $A(\lambda)$ is not invertible because $\lambda^2$ is not a constant.

    The $\lambda$-matrix $B(\lambda) = \begin{pmatrix} 1 & \lambda \\ 0 & 1 \end{pmatrix}$ is invertible: $\det B(\lambda) = 1$, and $B(\lambda)^{-1} = \begin{pmatrix} 1 & -\lambda \\ 0 & 1 \end{pmatrix}$.

!!! example "Example 13B.2"
    The characteristic matrix of $A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}$ is

    $$
    \lambda I - A = \begin{pmatrix} \lambda - 2 & -1 \\ 0 & \lambda - 3 \end{pmatrix}
    $$

    $\det(\lambda I - A) = (\lambda - 2)(\lambda - 3) = \lambda^2 - 5\lambda + 6$.

---

## 13B.2 Elementary Operations on Lambda-Matrices

<div class="context-flow" markdown>

Elementary row/column operations on numerical matrices generalize to $\lambda$-matrices → Note that the scaling factor must be a **nonzero constant** (not a polynomial) → The three types of elementary operations preserve **equivalence**

</div>

!!! definition "Definition 13B.4 (Elementary operations on lambda-matrices)"
    The **elementary operations** on $\lambda$-matrices are of three types:

    1. **Swap** two rows (columns): $r_i \leftrightarrow r_j$;
    2. **Multiply** a row (column) by a **nonzero constant** $c \in \mathbb{F} \setminus \{0\}$: $r_i \to c \cdot r_i$;
    3. **Add a polynomial multiple** of one row (column) to another: $r_i \to r_i + p(\lambda) r_j$ ($i \neq j$, $p(\lambda) \in \mathbb{F}[\lambda]$).

    Note that in type 2, the multiplier must be a nonzero **constant**, not an arbitrary polynomial.

!!! definition "Definition 13B.5 (Equivalence)"
    If a $\lambda$-matrix $A(\lambda)$ can be transformed into $B(\lambda)$ by a finite sequence of elementary operations, then $A(\lambda)$ and $B(\lambda)$ are called **equivalent**, denoted $A(\lambda) \sim B(\lambda)$. Equivalently, there exist invertible $\lambda$-matrices $P(\lambda), Q(\lambda)$ such that

    $$
    B(\lambda) = P(\lambda) A(\lambda) Q(\lambda)
    $$

!!! theorem "Theorem 13B.2 (Equivalence is an equivalence relation)"
    Equivalence of $\lambda$-matrices satisfies reflexivity, symmetry, and transitivity.

??? proof "Proof"
    Reflexivity: take $P = Q = I$. Symmetry: if $B = PAQ$, then $A = P^{-1}BQ^{-1}$, where $P^{-1}$ and $Q^{-1}$ are invertible $\lambda$-matrices. Transitivity: if $B = P_1AQ_1$ and $C = P_2BQ_2$, then $C = (P_2P_1)A(Q_1Q_2)$. $\blacksquare$

!!! example "Example 13B.3"
    Performing elementary operations on $A(\lambda) = \begin{pmatrix} \lambda - 2 & -1 \\ 0 & \lambda - 3 \end{pmatrix}$.

    $r_1 \to r_1 + \frac{1}{\lambda-3} r_2$? This is not allowed since $\frac{1}{\lambda-3}$ is not a polynomial.

    A correct approach uses polynomial division. For example, $c_2 \to c_2 + c_1$:

    $$
    \begin{pmatrix} \lambda - 2 & \lambda - 3 \\ 0 & \lambda - 3 \end{pmatrix}
    $$

    Then $c_1 \to c_1 - c_2$: $\begin{pmatrix} 1 & \lambda - 3 \\ -(\lambda-3) & \lambda - 3 \end{pmatrix}$. Further operations reduce this to Smith normal form.

---

## 13B.3 Smith Normal Form

<div class="context-flow" markdown>

Every $\lambda$-matrix is equivalent to a unique **diagonal form** → Diagonal entries satisfy a **divisibility chain** $d_1 \mid d_2 \mid \cdots$ → This is the central theorem of $\lambda$-matrix theory

</div>

!!! definition "Definition 13B.6 (Smith normal form)"
    The **Smith normal form** of an $m \times n$ ($m \leq n$) $\lambda$-matrix $A(\lambda)$ of rank $r$ is the diagonal matrix

    $$
    S(\lambda) = \begin{pmatrix} d_1(\lambda) & & & & \\ & d_2(\lambda) & & & \\ & & \ddots & & \\ & & & d_r(\lambda) & \\ & & & & O \end{pmatrix}
    $$

    where the $d_i(\lambda)$ are monic polynomials (leading coefficient equal to 1) satisfying the **divisibility chain**:

    $$
    d_1(\lambda) \mid d_2(\lambda) \mid \cdots \mid d_r(\lambda)
    $$

!!! theorem "Theorem 13B.3 (Existence and uniqueness of Smith normal form)"
    Every nonzero $\lambda$-matrix $A(\lambda)$ is equivalent to a unique Smith normal form.

??? proof "Proof"
    **Existence**: Constructive reduction via elementary operations.

    **Step 1**: By row and column swaps, move the nonzero entry of lowest degree to position $(1,1)$.

    **Step 2**: If the $(1,1)$ entry does not divide some entry in the first row or column, use polynomial division with remainder and elementary operations to reduce the degree of the $(1,1)$ entry. Repeat until the $(1,1)$ entry divides all entries in the first row and column.

    **Step 3**: Use the $(1,1)$ entry to eliminate all other entries in the first row and column. If the $(1,1)$ entry does not divide some entry outside the first row and column, add that entry's row to the first row and return to Step 2.

    **Step 4**: Eventually the $(1,1)$ entry divides all other entries. Make it monic and call it $d_1(\lambda)$. Recursively apply the procedure to the lower-right submatrix.

    **Uniqueness**: Guaranteed by the uniqueness of determinantal divisors (Theorem 13B.8). $\blacksquare$

!!! theorem "Theorem 13B.4 (Elementary operations preserve equivalence)"
    Elementary operations do not change the Smith normal form of a $\lambda$-matrix. Equivalently, $A(\lambda) \sim B(\lambda)$ if and only if they have the same Smith normal form.

??? proof "Proof"
    Each elementary operation can be expressed as left or right multiplication by an invertible $\lambda$-matrix (an elementary matrix) whose determinant is a nonzero constant. This does not change the determinantal divisors at any order (see Section 13B.5), so the Smith normal form remains unchanged. $\blacksquare$

!!! example "Example 13B.4"
    Find the Smith normal form of $A(\lambda) = \begin{pmatrix} \lambda & \lambda^2 \\ 1 & \lambda \end{pmatrix}$.

    $r_1 \leftrightarrow r_2$: $\begin{pmatrix} 1 & \lambda \\ \lambda & \lambda^2 \end{pmatrix}$.

    $r_2 \to r_2 - \lambda r_1$: $\begin{pmatrix} 1 & \lambda \\ 0 & 0 \end{pmatrix}$.

    $c_2 \to c_2 - \lambda c_1$: $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$.

    The Smith normal form is $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$, with rank $1$.

!!! example "Example 13B.5"
    Find the Smith normal form of $A(\lambda) = \begin{pmatrix} \lambda - 1 & 0 \\ 0 & (\lambda-1)^2 \end{pmatrix}$.

    This is already diagonal with $(\lambda-1) \mid (\lambda-1)^2$, so the Smith normal form is

    $$
    S(\lambda) = \begin{pmatrix} \lambda - 1 & 0 \\ 0 & (\lambda-1)^2 \end{pmatrix}
    $$

    $d_1(\lambda) = \lambda - 1$, $d_2(\lambda) = (\lambda-1)^2$.

---

## 13B.4 Invariant Factors and Elementary Divisors

<div class="context-flow" markdown>

Diagonal entries $d_1, \ldots, d_r$ of the Smith normal form = **invariant factors** → Each invariant factor decomposes into powers of irreducible polynomials = **elementary divisors** → Elementary divisors completely determine the Jordan block structure

</div>

!!! definition "Definition 13B.7 (Invariant factors)"
    The diagonal entries $d_1(\lambda), d_2(\lambda), \ldots, d_r(\lambda)$ of the Smith normal form of a $\lambda$-matrix $A(\lambda)$ (of rank $r$) are called the **invariant factors** of $A(\lambda)$. They satisfy the divisibility chain $d_1 \mid d_2 \mid \cdots \mid d_r$.

!!! definition "Definition 13B.8 (Elementary divisors)"
    Decompose each invariant factor $d_k(\lambda)$ into a product of powers of irreducible polynomials in $\mathbb{F}[\lambda]$:

    $$
    d_k(\lambda) = p_1(\lambda)^{e_{k1}} p_2(\lambda)^{e_{k2}} \cdots p_s(\lambda)^{e_{ks}}
    $$

    where all $e_{ki} \geq 0$. The collection of all factors $p_i(\lambda)^{e_{ki}}$ with $e_{ki} \geq 1$ (grouped by irreducible polynomial) is called the set of **elementary divisors** of $A(\lambda)$.

!!! theorem "Theorem 13B.5 (Relationship between invariant factors and elementary divisors)"
    Let $A(\lambda)$ have invariant factors $d_1(\lambda), \ldots, d_r(\lambda)$, and let $p_1(\lambda), \ldots, p_s(\lambda)$ be the irreducible polynomials involved. Then

    $$
    d_k(\lambda) = \prod_{j=1}^s p_j(\lambda)^{e_{kj}}
    $$

    where $0 \leq e_{1j} \leq e_{2j} \leq \cdots \leq e_{rj}$ (the divisibility chain condition).

    Conversely, given the set of elementary divisors, the invariant factors can be uniquely recovered: for each irreducible polynomial $p_j$, assign its powers $e_{1j} \leq \cdots \leq e_{rj}$ from the last invariant factor backward.

??? proof "Proof"
    The divisibility chain $d_1 \mid d_2 \mid \cdots \mid d_r$ directly implies $e_{1j} \leq e_{2j} \leq \cdots \leq e_{rj}$.

    Reverse recovery: the last invariant factor $d_r = \prod_j p_j^{e_{rj}}$ (taking the highest power of each irreducible factor). $d_{r-1} = \prod_j p_j^{e_{r-1,j}}$ (taking the next highest powers), and so on. Uniqueness follows from unique factorization. $\blacksquare$

!!! theorem "Theorem 13B.6 (Invariant factors and characteristic polynomial)"
    Let $A \in \mathbb{F}^{n \times n}$ and let the invariant factors of $\lambda I - A$ be $d_1(\lambda), \ldots, d_n(\lambda)$. Then

    $$
    d_1(\lambda) d_2(\lambda) \cdots d_n(\lambda) = \det(\lambda I - A) = p_A(\lambda)
    $$

    That is, the product of the invariant factors equals the characteristic polynomial. The last invariant factor $d_n(\lambda)$ equals the minimal polynomial $m_A(\lambda)$.

??? proof "Proof"
    The Smith normal form $S(\lambda)$ is equivalent to $\lambda I - A$, so $\det S(\lambda) = c \cdot \det(\lambda I - A)$ for a nonzero constant $c$. Since $S(\lambda) = \operatorname{diag}(d_1, \ldots, d_n)$ and all $d_i$ are monic, $c = 1$.

    For the minimal polynomial, the proof that $d_n(\lambda) = m_A(\lambda)$ uses the relationship between invariant factors and annihilating polynomials: $d_n(\lambda)$ divides every polynomial that annihilates $A$, and $d_n(A) = 0$ (by a generalization of the Cayley-Hamilton theorem). $\blacksquare$

!!! example "Example 13B.6"
    Let $A = \begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}$. The characteristic matrix is

    $$
    \lambda I - A = \begin{pmatrix} \lambda & 0 & 0 \\ -1 & \lambda & 0 \\ 0 & -1 & \lambda \end{pmatrix}
    $$

    After elementary operations: $c_1 \to c_1 + \lambda c_2$, $c_2 \to c_2 + \lambda c_3$:

    $$
    \begin{pmatrix} \lambda^2 - 1 & \lambda^2 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & \lambda \end{pmatrix} \to \cdots \to \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & \lambda^3 \end{pmatrix}
    $$

    Invariant factors: $d_1 = 1$, $d_2 = 1$, $d_3 = \lambda^3$. Characteristic polynomial $p_A(\lambda) = \lambda^3$, minimal polynomial $m_A(\lambda) = \lambda^3$.

    Over $\mathbb{C}$, the only elementary divisor is $\lambda^3$, corresponding to the $3 \times 3$ Jordan block $J_3(0)$.

!!! example "Example 13B.7"
    Suppose the Smith normal form of $\lambda I - A$ for a matrix $A$ is

    $$
    \operatorname{diag}(1, 1, \lambda - 1, (\lambda-1)(\lambda-2)^2)
    $$

    Invariant factors: $d_1 = d_2 = 1$, $d_3 = \lambda - 1$, $d_4 = (\lambda-1)(\lambda-2)^2$.

    Characteristic polynomial: $p_A(\lambda) = (\lambda-1)^2(\lambda-2)^2$.

    Minimal polynomial: $m_A(\lambda) = (\lambda-1)(\lambda-2)^2$.

    Elementary divisors (over $\mathbb{Q}$, $\mathbb{R}$, or $\mathbb{C}$): $(\lambda-1)^1$ (from $d_3$), $(\lambda-1)^1$ (from $d_4$), $(\lambda-2)^2$ (from $d_4$).

---

## 13B.5 Determinantal Divisors

<div class="context-flow" markdown>

Determinantal divisor $D_k(\lambda)$ = GCD of all $k \times k$ minors → $d_k = D_k/D_{k-1}$ → Determinantal divisors provide a **practical method for computing invariant factors** without performing elementary operations

</div>

!!! definition "Definition 13B.9 (Determinantal divisors)"
    Let $A(\lambda)$ be a $\lambda$-matrix of rank $r$. For $k = 1, 2, \ldots, r$, the monic greatest common divisor $D_k(\lambda)$ of all $k \times k$ minors of $A(\lambda)$ is called the **$k$-th determinantal divisor** of $A(\lambda)$. By convention, $D_0(\lambda) = 1$.

!!! theorem "Theorem 13B.7 (Divisibility of determinantal divisors)"
    The determinantal divisors satisfy $D_k(\lambda) \mid D_{k+1}(\lambda)$ ($k = 1, \ldots, r-1$). More precisely, $D_k(\lambda)$ divides every $k \times k$ minor of $A(\lambda)$.

??? proof "Proof"
    Every $(k+1) \times (k+1)$ minor, expanded along some row, is an $\mathbb{F}[\lambda]$-linear combination of $k \times k$ minors. Therefore $D_k(\lambda)$ divides every $(k+1) \times (k+1)$ minor, so $D_k(\lambda) \mid D_{k+1}(\lambda)$. $\blacksquare$

!!! theorem "Theorem 13B.8 (Invariant factors from determinantal divisors)"
    The invariant factors $d_k(\lambda)$ and determinantal divisors $D_k(\lambda)$ of a $\lambda$-matrix $A(\lambda)$ are related by

    $$
    d_k(\lambda) = \frac{D_k(\lambda)}{D_{k-1}(\lambda)}, \quad k = 1, 2, \ldots, r
    $$

    where $D_0(\lambda) = 1$.

??? proof "Proof"
    Let the Smith normal form of $A(\lambda)$ be $S(\lambda) = \operatorname{diag}(d_1, \ldots, d_r, 0, \ldots, 0)$. Since elementary operations do not change determinantal divisors (invertible $\lambda$-matrices have nonzero constant determinants, and the Binet-Cauchy formula ensures the GCD of minors is preserved), $A(\lambda)$ and $S(\lambda)$ have the same determinantal divisors.

    The nonzero $k \times k$ minors of $S(\lambda)$ are exactly $\prod_{i \in I} d_i(\lambda)$ ($|I| = k$, $I \subseteq \{1, \ldots, r\}$). By the divisibility chain $d_1 \mid d_2 \mid \cdots$, the GCD is $d_1 d_2 \cdots d_k$. Therefore

    $$
    D_k(\lambda) = d_1(\lambda)d_2(\lambda)\cdots d_k(\lambda)
    $$

    and so $d_k(\lambda) = D_k/D_{k-1}$. $\blacksquare$

!!! corollary "Corollary 13B.1"
    Equivalent $\lambda$-matrices have the same determinantal divisors, invariant factors, and elementary divisors. Conversely, having the same invariant factors (or equivalently, the same determinantal divisors) guarantees equivalence of $\lambda$-matrices.

!!! example "Example 13B.8"
    Use determinantal divisors to find the invariant factors of $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$.

    Characteristic matrix $\lambda I - A = \begin{pmatrix} \lambda - 2 & -1 \\ -1 & \lambda - 2 \end{pmatrix}$.

    $D_1(\lambda)$: The $1 \times 1$ minors are $\lambda - 2, -1, -1, \lambda - 2$. $\gcd = 1$.

    $D_2(\lambda)$: $\det(\lambda I - A) = (\lambda-2)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda-1)(\lambda-3)$.

    Invariant factors: $d_1 = D_1/D_0 = 1$, $d_2 = D_2/D_1 = (\lambda-1)(\lambda-3)$.

!!! example "Example 13B.9"
    Find the determinantal divisors and invariant factors of $A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$.

    $\lambda I - A = \operatorname{diag}(\lambda-1, \lambda-1, \lambda-2)$.

    $D_1 = \gcd(\lambda-1, \lambda-1, \lambda-2) = 1$.

    $D_2$: The $2 \times 2$ minors are $(\lambda-1)^2$, $(\lambda-1)(\lambda-2)$, $(\lambda-1)(\lambda-2)$. $D_2 = \gcd = \lambda - 1$.

    $D_3 = (\lambda-1)^2(\lambda-2)$.

    Invariant factors: $d_1 = 1$, $d_2 = \lambda - 1$, $d_3 = (\lambda-1)(\lambda-2)$.

---

## 13B.6 Companion Matrix and Rational Canonical Form

<div class="context-flow" markdown>

Degree-$n$ polynomial $d(\lambda)$ $\to$ **Companion matrix** $C(d)$: both its characteristic polynomial and minimal polynomial equal $d(\lambda)$ → Matrix $A$ is similar to a direct sum of companion matrices $\Leftrightarrow$ Rational canonical form → **Valid over any field**

</div>

!!! definition "Definition 13B.10 (Companion matrix)"
    Let $d(\lambda) = \lambda^m + a_{m-1}\lambda^{m-1} + \cdots + a_1\lambda + a_0$ be a monic polynomial in $\mathbb{F}[\lambda]$. Its **companion matrix** is defined as

    $$
    C(d) = \begin{pmatrix} 0 & 0 & \cdots & 0 & -a_0 \\ 1 & 0 & \cdots & 0 & -a_1 \\ 0 & 1 & \cdots & 0 & -a_2 \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & 1 & -a_{m-1} \end{pmatrix}_{m \times m}
    $$

!!! theorem "Theorem 13B.9 (Properties of the companion matrix)"
    Let $d(\lambda) = \lambda^m + a_{m-1}\lambda^{m-1} + \cdots + a_0$ be a monic polynomial and $C = C(d)$ its companion matrix. Then

    1. The characteristic polynomial of $C$ is $p_C(\lambda) = d(\lambda)$;
    2. The minimal polynomial of $C$ is $m_C(\lambda) = d(\lambda)$;
    3. The Smith normal form of $\lambda I - C$ is $\operatorname{diag}(1, 1, \ldots, 1, d(\lambda))$.

??? proof "Proof"
    **(1)** Expand $\det(\lambda I - C)$ along the last column:

    $$
    \det(\lambda I - C) = \det\begin{pmatrix} \lambda & 0 & \cdots & 0 & a_0 \\ -1 & \lambda & \cdots & 0 & a_1 \\ 0 & -1 & \cdots & 0 & a_2 \\ \vdots & & \ddots & & \vdots \\ 0 & 0 & \cdots & -1 & \lambda + a_{m-1} \end{pmatrix}
    $$

    Expanding along the first column and using recursion gives $\det(\lambda I - C) = \lambda^m + a_{m-1}\lambda^{m-1} + \cdots + a_0 = d(\lambda)$.

    **(2)** Let $\mathbf{e}_1 = (1, 0, \ldots, 0)^T$. Direct computation shows $C\mathbf{e}_1 = \mathbf{e}_2$, $C^2\mathbf{e}_1 = \mathbf{e}_3$, ..., $C^{m-1}\mathbf{e}_1 = \mathbf{e}_m$. Therefore $\{\mathbf{e}_1, C\mathbf{e}_1, \ldots, C^{m-1}\mathbf{e}_1\}$ forms a basis of $\mathbb{F}^m$.

    If $f(C) = 0$ with $\deg f < m$, then $f(C)\mathbf{e}_1 = f_0\mathbf{e}_1 + f_1C\mathbf{e}_1 + \cdots + f_{m-1}C^{m-1}\mathbf{e}_1 = \mathbf{0}$, and by linear independence of the basis, $f \equiv 0$. So the minimal polynomial has degree at least $m$, and by the Cayley-Hamilton theorem, $m_C = p_C = d$.

    **(3)** Follows directly from (1) and (2). $\blacksquare$

!!! theorem "Theorem 13B.10 (Rational canonical form)"
    Let $A \in \mathbb{F}^{n \times n}$ and let the invariant factors of $\lambda I - A$ be $d_1(\lambda), \ldots, d_n(\lambda)$, where $d_1 = \cdots = d_{n-s} = 1$ and $\deg d_{n-s+1} \geq 1, \ldots, \deg d_n \geq 1$. Let the nontrivial invariant factors be $f_1 = d_{n-s+1}, \ldots, f_s = d_n$ ($f_1 \mid f_2 \mid \cdots \mid f_s$). Then $A$ is similar to

    $$
    R = \begin{pmatrix} C(f_1) & & \\ & C(f_2) & \\ & & \ddots & \\ & & & C(f_s) \end{pmatrix}
    $$

    This is the **rational canonical form** (also called the **Frobenius normal form**) of $A$. It exists over any field $\mathbb{F}$.

??? proof "Proof"
    $\lambda I - A$ and $\lambda I - R$ have the same Smith normal form. Indeed, $\lambda I - R = \operatorname{diag}(\lambda I - C(f_1), \ldots, \lambda I - C(f_s))$, and the Smith normal form of $\lambda I - C(f_k)$ is $\operatorname{diag}(1, \ldots, 1, f_k(\lambda))$ (Theorem 13B.9). Combining these gives the Smith normal form $\operatorname{diag}(1, \ldots, 1, f_1, \ldots, f_s)$, which matches that of $\lambda I - A$.

    Since the Smith normal form uniquely determines the equivalence class (Corollary 13B.1), $\lambda I - A \sim \lambda I - R$, meaning there exist invertible $\lambda$-matrices $P(\lambda), Q(\lambda)$ with $P(\lambda)(\lambda I - A)Q(\lambda) = \lambda I - R$.

    From this, one can prove (using a more refined argument bridging $\lambda$-matrix equivalence and numerical matrix similarity) that $A$ and $R$ are similar. $\blacksquare$

!!! theorem "Theorem 13B.11 (Uniqueness of the rational canonical form)"
    The rational canonical form is unique up to similarity. That is, if $A$ is simultaneously similar to two rational canonical forms $R_1$ and $R_2$, then $R_1 = R_2$ (with the companion matrix blocks in the same order).

??? proof "Proof"
    $A \sim R_1$ and $A \sim R_2$ imply that $\lambda I - R_1$ and $\lambda I - R_2$ have the same Smith normal form. By the Smith normal form of companion matrices (Theorem 13B.9), the nontrivial invariant factors are exactly the polynomials defining the companion matrix blocks. By uniqueness of the Smith normal form, these polynomials are identical, so $R_1 = R_2$. $\blacksquare$

!!! example "Example 13B.10"
    Suppose a $4 \times 4$ matrix $A$ has invariant factors $d_1 = 1$, $d_2 = 1$, $d_3 = \lambda - 1$, $d_4 = (\lambda-1)(\lambda-2)^2$.

    Nontrivial invariant factors: $f_1 = \lambda - 1$, $f_2 = (\lambda-1)(\lambda-2)^2 = \lambda^3 - 5\lambda^2 + 8\lambda - 4$.

    $C(f_1) = (1)$ (a $1 \times 1$ matrix).

    $$
    C(f_2) = \begin{pmatrix} 0 & 0 & 4 \\ 1 & 0 & -8 \\ 0 & 1 & 5 \end{pmatrix}
    $$

    Rational canonical form:

    $$
    R = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 4 \\ 0 & 1 & 0 & -8 \\ 0 & 0 & 1 & 5 \end{pmatrix}
    $$

!!! example "Example 13B.11"
    Find the rational canonical form of $A = \begin{pmatrix} 0 & 1 \\ -6 & 5 \end{pmatrix}$.

    $p_A(\lambda) = \lambda^2 - 5\lambda + 6 = (\lambda-2)(\lambda-3)$.

    $D_1 = \gcd(\lambda, -1, 1, 5-\lambda) = 1$, $D_2 = (\lambda-2)(\lambda-3)$.

    Invariant factors: $d_1 = 1$, $d_2 = (\lambda-2)(\lambda-3)$. Nontrivial invariant factor: $f_1 = \lambda^2 - 5\lambda + 6$.

    $$
    C(f_1) = \begin{pmatrix} 0 & -6 \\ 1 & 5 \end{pmatrix}
    $$

    Therefore $A$ is similar to $R = \begin{pmatrix} 0 & -6 \\ 1 & 5 \end{pmatrix}$. Indeed, taking $P = \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix}$, we get $P^{-1}AP = R$.

---

## 13B.7 Relationship with Jordan Normal Form

<div class="context-flow" markdown>

Over an algebraically closed field $\mathbb{C}$: invariant factors factor completely → Elementary divisors $(\lambda - \lambda_i)^{e_i}$ → Companion matrix $C((\lambda-\lambda_i)^{e_i})$ is similar to Jordan block $J_{e_i}(\lambda_i)$ → Rational canonical form becomes **Jordan normal form**

</div>

!!! theorem "Theorem 13B.12 (Companion matrix and Jordan block)"
    Over an algebraically closed field $\mathbb{F}$, the companion matrix $C((\lambda - a)^m)$ is similar to the Jordan block $J_m(a)$.

??? proof "Proof"
    Both $C((\lambda - a)^m)$ and $J_m(a)$ have characteristic polynomial and minimal polynomial equal to $(\lambda - a)^m$. Their Smith normal forms are identical: $\operatorname{diag}(1, \ldots, 1, (\lambda-a)^m)$. By uniqueness of the rational canonical form, $C((\lambda-a)^m) \sim J_m(a)$.

    Concretely, let $C = C((\lambda-a)^m)$ and $\mathbf{e}_1 = (1,0,\ldots,0)^T$. Define $\mathbf{v}_m = \mathbf{e}_1$, $\mathbf{v}_{m-1} = (C - aI)\mathbf{v}_m$, ..., $\mathbf{v}_1 = (C-aI)^{m-1}\mathbf{v}_m$. Then $P = (\mathbf{v}_1, \ldots, \mathbf{v}_m)$ satisfies $P^{-1}CP = J_m(a)$. $\blacksquare$

!!! theorem "Theorem 13B.13 (From rational canonical form to Jordan normal form)"
    Let $A \in \mathbb{C}^{n \times n}$. After decomposing invariant factors into elementary divisors, each companion matrix $C(f_k)$ in the rational canonical form can be further transformed by similarity into a direct sum of Jordan blocks:

    If $f_k(\lambda) = \prod_{i} (\lambda - \lambda_i)^{e_{ki}}$, then

    $$
    C(f_k) \sim \bigoplus_{i} J_{e_{ki}}(\lambda_i)
    $$

    Therefore the Jordan normal form of $A$ is

    $$
    J = \bigoplus_{k,i} J_{e_{ki}}(\lambda_i)
    $$

??? proof "Proof"
    For $C(f_k)$, both the characteristic and minimal polynomials equal $f_k(\lambda)$. Since $f_k$ factors into coprime factors $f_k = \prod_i (\lambda - \lambda_i)^{e_{ki}}$, by the elementary divisor theory (or directly by the structure theorem for $\mathbb{F}[\lambda]$-modules), $C(f_k)$ is similar to $\bigoplus_i C((\lambda-\lambda_i)^{e_{ki}})$.

    By Theorem 13B.12, $C((\lambda-\lambda_i)^{e_{ki}}) \sim J_{e_{ki}}(\lambda_i)$. $\blacksquare$

!!! example "Example 13B.12"
    Continuing from Example 13B.10. The invariant factors of $A$ are $d_3 = \lambda-1$ and $d_4 = (\lambda-1)(\lambda-2)^2$.

    Elementary divisors: $(\lambda-1)^1$ (from $d_3$), $(\lambda-1)^1$ (from $d_4$), $(\lambda-2)^2$ (from $d_4$).

    Rational canonical form $\to$ Jordan normal form:

    $$
    R = \begin{pmatrix} C(\lambda-1) & \\ & C((\lambda-1)(\lambda-2)^2) \end{pmatrix} \sim \begin{pmatrix} 1 & & & \\ & 1 & & \\ & & 2 & 1 \\ & & 0 & 2 \end{pmatrix} = J
    $$

!!! example "Example 13B.13"
    Consider $A \in \mathbb{R}^{4 \times 4}$ with invariant factors $d_1 = d_2 = 1$, $d_3 = \lambda^2+1$, $d_4 = (\lambda^2+1)^2$.

    Over $\mathbb{R}$, $\lambda^2 + 1$ is irreducible, so the rational canonical form is

    $$
    R = \begin{pmatrix} C(\lambda^2+1) & \\ & C((\lambda^2+1)^2) \end{pmatrix} = \begin{pmatrix} 0 & -1 & & \\ 1 & 0 & & \\ & & 0 & 0 & 0 & 1 \\ & & 1 & 0 & 0 & 0 \\ & & 0 & 1 & 0 & -2 \\ & & 0 & 0 & 1 & 0 \end{pmatrix}
    $$

    This is the "simplest form" over $\mathbb{R}$ — it cannot be further simplified because the eigenvalues $\pm i$ are not in $\mathbb{R}$.

    Over $\mathbb{C}$, $\lambda^2+1 = (\lambda-i)(\lambda+i)$, and the elementary divisors are $(\lambda-i)^1, (\lambda+i)^1, (\lambda-i)^2, (\lambda+i)^2$. The Jordan normal form is $\operatorname{diag}(i, -i, J_2(i), J_2(-i))$.

---

## 13B.8 Complete Invariants of Matrix Similarity

<div class="context-flow" markdown>

**Ultimate criterion**: $A \sim B$ $\Leftrightarrow$ $\lambda I - A$ and $\lambda I - B$ have the same Smith normal form $\Leftrightarrow$ same invariant factors $\Leftrightarrow$ same set of elementary divisors → Solves the matrix similarity problem once and for all over any field

</div>

!!! definition "Definition 13B.11 (Complete invariants)"
    A collection of invariants of a matrix is called a **complete system of invariants** for similarity if: $A \sim B$ if and only if $A$ and $B$ have the same values for all invariants in the collection.

!!! theorem "Theorem 13B.14 (Criterion for matrix similarity)"
    Let $A, B \in \mathbb{F}^{n \times n}$. The following conditions are equivalent:

    1. $A$ and $B$ are similar (there exists an invertible $P$ with $B = P^{-1}AP$);
    2. $\lambda I - A$ and $\lambda I - B$ are equivalent (as $\lambda$-matrices);
    3. $\lambda I - A$ and $\lambda I - B$ have the same Smith normal form;
    4. $\lambda I - A$ and $\lambda I - B$ have the same determinantal divisors;
    5. $\lambda I - A$ and $\lambda I - B$ have the same invariant factors;
    6. $\lambda I - A$ and $\lambda I - B$ have the same set of elementary divisors;
    7. $A$ and $B$ have the same rational canonical form.

??? proof "Proof"
    **(1)$\Rightarrow$(2):** If $B = P^{-1}AP$ ($P$ an invertible constant matrix), then

    $$
    \lambda I - B = \lambda I - P^{-1}AP = P^{-1}(\lambda I - A)P
    $$

    Both $P$ and $P^{-1}$ are invertible $\lambda$-matrices (with nonzero constant determinants), so $\lambda I - A \sim \lambda I - B$.

    **(2)$\Leftrightarrow$(3)$\Leftrightarrow$(4)$\Leftrightarrow$(5)$\Leftrightarrow$(6):** By uniqueness of the Smith normal form and the bijective correspondence among determinantal divisors, invariant factors, and elementary divisors.

    **(2)$\Rightarrow$(1):** This is the most nontrivial direction. One must show that $\lambda$-matrix equivalence implies similarity of numerical matrices. Given $P(\lambda)(\lambda I - A)Q(\lambda) = \lambda I - B$, perform polynomial division of $P(\lambda)$ and $Q(\lambda)$ by $(\lambda I - A)$ and analyze the result to construct a constant invertible matrix $P$ with $P^{-1}AP = B$.

    **(5)$\Leftrightarrow$(7):** The rational canonical form is uniquely determined by the invariant factors. $\blacksquare$

!!! theorem "Theorem 13B.15 (Hierarchy of invariants)"
    For $A \in \mathbb{F}^{n \times n}$, the following invariants provide increasing amounts of information:

    - **Characteristic polynomial** $p_A(\lambda)$: equals the product of invariant factors $\prod d_k$. Knowing only the characteristic polynomial is **not sufficient** to determine similarity;
    - **Minimal polynomial** $m_A(\lambda)$: equals the last invariant factor $d_n$. The characteristic polynomial together with the minimal polynomial is still not sufficient to determine similarity;
    - **Set of invariant factors** $\{d_1, \ldots, d_n\}$: a **complete system of invariants**, fully determining similarity.

??? proof "Proof"
    Counterexamples show the first two items are insufficient:

    $A_1 = \operatorname{diag}(1, 1, 2)$ and $A_2 = \begin{pmatrix} 1 & 1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$ have the same characteristic polynomial $(\lambda-1)^2(\lambda-2)$ but different invariant factors: $A_1$ has invariant factors $1, \lambda-1, (\lambda-1)(\lambda-2)$; $A_2$ has invariant factors $1, 1, (\lambda-1)^2(\lambda-2)$. Therefore $A_1 \not\sim A_2$.

    The two have different minimal polynomials: $m_{A_1} = (\lambda-1)(\lambda-2)$ and $m_{A_2} = (\lambda-1)^2(\lambda-2)$. But even when both the characteristic polynomial and minimal polynomial agree, the invariant factors may still differ (counterexamples can be constructed for $n \geq 4$). $\blacksquare$

!!! theorem "Theorem 13B.16 (Refinement of the Cayley-Hamilton theorem)"
    Let the invariant factors of $A$ be $d_1 \mid d_2 \mid \cdots \mid d_n$. Then

    1. $d_n(A) = 0$ (the minimal polynomial annihilates $A$);
    2. A polynomial $f(\lambda) \in \mathbb{F}[\lambda]$ satisfies $f(A) = 0$ if and only if $d_n \mid f$;
    3. $p_A(\lambda) = d_1 \cdots d_n$ and $d_n \mid p_A$ (the Cayley-Hamilton theorem $p_A(A) = 0$ is a consequence of $d_n(A) = 0$).

??? proof "Proof"
    (1) $A$ is similar to the rational canonical form $R = \bigoplus C(f_k)$. By Theorem 13B.9, $f_k(C(f_k)) = 0$. By the divisibility chain $f_k \mid f_s = d_n$, we have $d_n(C(f_k)) = 0$. Therefore $d_n(R) = 0$, and by similarity, $d_n(A) = 0$.

    (2) $f(A) = 0$ $\Leftrightarrow$ $f(R) = 0$ $\Leftrightarrow$ $f(C(f_k)) = 0$ for all $k$ $\Leftrightarrow$ $f_k \mid f$ for all $k$ $\Leftrightarrow$ $d_n \mid f$ (since $d_n = f_s$ is the largest invariant factor). $\blacksquare$

!!! example "Example 13B.14"
    Determine whether $A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 2 \end{pmatrix}$ and $B = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 2 \end{pmatrix}$ are similar.

    $A$: $D_1 = 1$, $D_2 = \gcd((\lambda-1)(\lambda-2), (\lambda-1)(\lambda-2), (\lambda-2)^2) = \lambda-2$. $D_3 = (\lambda-1)(\lambda-2)^2$.

    Invariant factors: $d_1 = 1$, $d_2 = \lambda-2$, $d_3 = (\lambda-1)(\lambda-2)$.

    $B$: $\lambda I - B = \begin{pmatrix} \lambda-1 & 0 & 0 \\ 0 & \lambda-2 & -1 \\ 0 & 0 & \lambda-2 \end{pmatrix}$.

    The $2 \times 2$ minors include: $(\lambda-1)(\lambda-2)$, $-(\lambda-1)$, $(\lambda-1)(\lambda-2)$, $(\lambda-2)^2$, and others. Since $\gcd$ includes $-(\lambda-1)$ and $(\lambda-2)^2$ which are coprime, $D_2 = 1$.

    $D_3 = (\lambda-1)(\lambda-2)^2$.

    Invariant factors of $B$: $d_1 = 1$, $d_2 = 1$, $d_3 = (\lambda-1)(\lambda-2)^2$.

    Since $A$ and $B$ have different invariant factors, $A \not\sim B$.

    Intuition: The eigenvalue $2$ of $A$ has two linearly independent eigenvectors ($A$ is diagonalizable at $\lambda=2$), while the eigenvalue $2$ of $B$ has only one linearly independent eigenvector ($B$ has a nontrivial Jordan block at $\lambda=2$).

!!! example "Example 13B.15"
    The following two $4 \times 4$ matrices have the same characteristic polynomial $(\lambda-1)^4$ and the same minimal polynomial $(\lambda-1)^2$, but are not similar:

    $$
    A_1 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix}, \quad A_2 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

    $A_1$ ($J_2(1) \oplus J_2(1)$): invariant factors $1, 1, (\lambda-1)^2, (\lambda-1)^2$.

    $A_2$ ($J_2(1) \oplus J_1(1) \oplus J_1(1)$): invariant factors $1, \lambda-1, \lambda-1, (\lambda-1)^2$.

    Both have characteristic polynomial $(\lambda-1)^4$ and minimal polynomial $(\lambda-1)^2$, but their invariant factors differ, confirming that the combination of characteristic polynomial and minimal polynomial is still not a complete system of invariants.
