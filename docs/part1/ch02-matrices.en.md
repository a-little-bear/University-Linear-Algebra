# Chapter 2  Matrices and Matrix Operations

<div class="context-flow" markdown>

**Prerequisites**: Chapter 1 augmented matrices · elementary row operations

**Chapter arc**: Matrix definition → addition/scalar multiplication → matrix multiplication (abstraction of $A\mathbf{x}=\mathbf{b}$) → transpose → inverse matrices → elementary matrices → block matrices → rank

**Further connections**：Matrices are central to quantum mechanics (observables as Hermitian matrices), graph theory (adjacency matrices), and Markov chains (transition matrices); sparse matrix storage and operations are key to large-scale scientific computing

</div>

A matrix is the core language of linear algebra. It is not only a compact representation tool for systems of linear equations, but also the fundamental carrier for describing linear transformations. This chapter systematically introduces the basic concepts and various operations on matrices — addition, scalar multiplication, matrix multiplication, transpose, and inversion — and discusses in depth important concepts such as elementary matrices, block matrices, and the rank of a matrix.

---

## 2.1 Definition and Basic Concepts of Matrices

<div class="context-flow" markdown>

**From systems to matrices**: The rectangular array of coefficients from Chapter 1 → now given an independent mathematical identity → zero matrix, identity matrix, diagonal matrix, and other special roles

</div>

!!! definition "Definition 2.1 (Matrix)"
    An $m \times n$ **matrix** is a rectangular array of elements arranged in $m$ rows and $n$ columns, written as

    $$
    A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}
    $$

    abbreviated as $A = (a_{ij})_{m \times n}$ or $A = (a_{ij})$. The entry $a_{ij}$ is located in the $i$-th row and $j$-th column. The set of all $m \times n$ real matrices is denoted $\mathbb{R}^{m \times n}$, and the set of all $m \times n$ complex matrices is denoted $\mathbb{C}^{m \times n}$.

!!! definition "Definition 2.2 (Special matrices)"
    Let $A = (a_{ij})$ be a matrix. The following are several important types of special matrices:

    1. **Zero matrix**: A matrix with all entries equal to zero, denoted $O$ or $O_{m \times n}$.
    2. **Square matrix**: A matrix with equal numbers of rows and columns, i.e., $m = n$, called a square matrix of order $n$.
    3. **Identity matrix**: A square matrix of order $n$ with all diagonal entries equal to $1$ and all other entries equal to $0$, denoted $I_n$ or $I$:

        $$
        I_n = \begin{pmatrix} 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & 1 \end{pmatrix}
        $$

    4. **Diagonal matrix**: $a_{ij} = 0$ ($i \ne j$), denoted $\operatorname{diag}(d_1, d_2, \ldots, d_n)$.
    5. **Upper triangular matrix**: $a_{ij} = 0$ ($i > j$), i.e., all entries below the main diagonal are zero.
    6. **Lower triangular matrix**: $a_{ij} = 0$ ($i < j$), i.e., all entries above the main diagonal are zero.
    7. **Symmetric matrix**: $a_{ij} = a_{ji}$ (for all $i, j$), i.e., $A = A^T$.
    8. **Skew-symmetric matrix** (antisymmetric matrix): $a_{ij} = -a_{ji}$, i.e., $A = -A^T$. In particular, the diagonal entries of a skew-symmetric matrix must be zero.

---

## 2.2 Matrix Addition and Scalar Multiplication

<div class="context-flow" markdown>

**Entry-wise operations**: Add/scalar-multiply corresponding entries → satisfies 8 axioms → $\mathbb{R}^{m \times n}$ itself forms a **vector space** (→ an important example in Chapter 4)

</div>

!!! definition "Definition 2.3 (Matrix addition and scalar multiplication)"
    Let $A = (a_{ij})$ and $B = (b_{ij})$ be matrices of the same size (both $m \times n$), and $c$ a scalar. Define:

    1. **Matrix addition**: $A + B = (a_{ij} + b_{ij})$, i.e., add corresponding entries.
    2. **Scalar multiplication**: $cA = (ca_{ij})$, i.e., multiply every entry by $c$.

!!! theorem "Theorem 2.1 (Properties of matrix addition and scalar multiplication)"
    Let $A, B, C$ be matrices of the same size, and $c, d$ scalars. Then:

    1. $A + B = B + A$ (commutativity of addition)
    2. $(A + B) + C = A + (B + C)$ (associativity of addition)
    3. $A + O = A$ (zero matrix is the additive identity)
    4. $A + (-A) = O$ (additive inverse)
    5. $c(A + B) = cA + cB$ (distributivity of scalar multiplication over matrix addition)
    6. $(c + d)A = cA + dA$ (distributivity of scalar addition over matrix)
    7. $(cd)A = c(dA)$ (associativity of scalar multiplication)
    8. $1 \cdot A = A$

??? proof "Proof"
    All these properties follow directly from the arithmetic properties of real numbers for matrix entries. For example, for property 5:

    $(c(A + B))_{ij} = c(a_{ij} + b_{ij}) = ca_{ij} + cb_{ij} = (cA)_{ij} + (cB)_{ij} = (cA + cB)_{ij}$.

    Since this holds for all $i, j$, we have $c(A + B) = cA + cB$. The rest are similar. $\blacksquare$

!!! note "Note"
    The above 8 properties show that the set of all $m \times n$ matrices forms a vector space $\mathbb{R}^{m \times n}$ under matrix addition and scalar multiplication, with dimension $mn$.

---

## 2.3 Matrix Multiplication

<div class="context-flow" markdown>

**Core operation**: $A\mathbf{x} = x_1\mathbf{a}_1 + \cdots + x_n\mathbf{a}_n$ (linear combination of columns) → matrix multiplication = composition of linear transformations (→ Chapter 5) → Note: **not commutative**

</div>

!!! definition "Definition 2.4 (Matrix multiplication)"
    Let $A = (a_{ij})$ be an $m \times p$ matrix and $B = (b_{jk})$ a $p \times n$ matrix. The **product** $C = AB$ is an $m \times n$ matrix whose $(i, k)$-entry is

    $$
    c_{ik} = \sum_{j=1}^{p} a_{ij} b_{jk} = a_{i1}b_{1k} + a_{i2}b_{2k} + \cdots + a_{ip}b_{pk}.
    $$

    That is, the $(i, k)$-entry of $C$ is the **inner product** of the $i$-th row of $A$ and the $k$-th column of $B$.

!!! note "Note"
    The matrix product $AB$ is defined only when the number of columns of $A$ equals the number of rows of $B$. If $A$ is $m \times p$ and $B$ is $p \times n$, then $AB$ is $m \times n$.

!!! theorem "Theorem 2.2 (Properties of matrix multiplication)"
    Assuming the matrix sizes make the following operations well-defined, and $c$ is a scalar:

    1. $A(BC) = (AB)C$ (associativity)
    2. $A(B + C) = AB + AC$ (left distributivity)
    3. $(A + B)C = AC + BC$ (right distributivity)
    4. $c(AB) = (cA)B = A(cB)$
    5. $I_m A = A = AI_n$ (where $A$ is $m \times n$)

??? proof "Proof"
    We prove associativity as an example. Let $A = (a_{ij})_{m \times p}$, $B = (b_{jk})_{p \times q}$, $C = (c_{kl})_{q \times n}$.

    $(A(BC))_{il} = \sum_{j=1}^p a_{ij}(BC)_{jl} = \sum_{j=1}^p a_{ij}\left(\sum_{k=1}^q b_{jk}c_{kl}\right) = \sum_{j=1}^p \sum_{k=1}^q a_{ij}b_{jk}c_{kl}$

    $((AB)C)_{il} = \sum_{k=1}^q (AB)_{ik}c_{kl} = \sum_{k=1}^q \left(\sum_{j=1}^p a_{ij}b_{jk}\right)c_{kl} = \sum_{k=1}^q \sum_{j=1}^p a_{ij}b_{jk}c_{kl}$

    Since finite sums can be reordered, the two expressions are equal. $\blacksquare$

!!! theorem "Theorem 2.3 (Matrix multiplication is not commutative)"
    Matrix multiplication is generally **not commutative**, i.e., $AB \neq BA$ (even when both products are defined).

!!! example "Example 2.1"
    Let

    $$
    A = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}, \quad B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
    $$

    Then

    $$
    AB = \begin{pmatrix} 2 & 1 \\ 1 & 0 \end{pmatrix}, \quad BA = \begin{pmatrix} 0 & 1 \\ 1 & 2 \end{pmatrix}
    $$

    Clearly $AB \neq BA$.

!!! example "Example 2.2"
    Another way matrix multiplication differs from real number multiplication: $AB = O$ does not imply $A = O$ or $B = O$. For example,

    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \quad B = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}, \quad AB = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix} = O.
    $$

!!! example "Example 2.3"
    Using matrix multiplication, the system $A\mathbf{x} = \mathbf{b}$ can be understood as: the vector $\mathbf{b}$ is a linear combination of the column vectors $\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n$ of $A$. Specifically, if $A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)$, then

    $$
    A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n.
    $$

---

## 2.4 Matrix Transpose

<div class="context-flow" markdown>

**Swapping rows and columns**: $(AB)^T = B^TA^T$ (order reversal!) → symmetric matrices $A = A^T$ will play central roles in the spectral theorem (Chapter 6) and orthogonality (Chapter 7)

</div>

!!! definition "Definition 2.5 (Transpose)"
    Let $A = (a_{ij})$ be an $m \times n$ matrix. The **transpose** of $A$, denoted $A^T$, is the $n \times m$ matrix obtained by turning the rows of $A$ into columns (or columns into rows), i.e., $(A^T)_{ij} = a_{ji}$.

!!! theorem "Theorem 2.4 (Properties of transpose)"
    Let $A, B$ be matrices of appropriate sizes, and $c$ a scalar. Then:

    1. $(A^T)^T = A$
    2. $(A + B)^T = A^T + B^T$
    3. $(cA)^T = cA^T$
    4. $(AB)^T = B^T A^T$ (note the order reversal)

??? proof "Proof"
    We focus on property 4. Let $A$ be $m \times p$ and $B$ be $p \times n$, so $AB$ is $m \times n$ and $(AB)^T$ is $n \times m$.

    $((AB)^T)_{ji} = (AB)_{ij} = \sum_{k=1}^p a_{ik}b_{kj}$

    $(B^T A^T)_{ji} = \sum_{k=1}^p (B^T)_{jk}(A^T)_{ki} = \sum_{k=1}^p b_{kj}a_{ik} = \sum_{k=1}^p a_{ik}b_{kj}$

    The two expressions are equal, so $(AB)^T = B^T A^T$. $\blacksquare$

!!! proposition "Proposition 2.1 (Properties of symmetric and skew-symmetric matrices)"
    1. If $A$ is a square matrix, then $A + A^T$ is symmetric and $A - A^T$ is skew-symmetric.
    2. Any square matrix $A$ can be uniquely decomposed as the sum of a symmetric and a skew-symmetric matrix:

    $$
    A = \frac{A + A^T}{2} + \frac{A - A^T}{2}.
    $$

??? proof "Proof"
    1. $(A + A^T)^T = A^T + (A^T)^T = A^T + A = A + A^T$, so it is symmetric. $(A - A^T)^T = A^T - A = -(A - A^T)$, so it is skew-symmetric.

    2. The formula clearly holds. Uniqueness: if $A = S + K$ ($S$ symmetric, $K$ skew-symmetric), then $A^T = S - K$, so $S = \frac{A + A^T}{2}$, $K = \frac{A - A^T}{2}$. $\blacksquare$

---

## 2.5 Inverse Matrices

<div class="context-flow" markdown>

**Invertibility = unique solution for every system**: $A$ invertible $\Leftrightarrow$ $A\mathbf{x}=\mathbf{b}$ has a unique solution for every $\mathbf{b}$ → Theorem 2.7 gives 9 equivalent conditions, threading core concepts from the first 4 chapters

</div>

!!! definition "Definition 2.6 (Invertible matrix)"
    Let $A$ be a square matrix of order $n$. If there exists a square matrix $B$ of order $n$ such that

    $$
    AB = BA = I_n,
    $$

    then $A$ is called **invertible** (or **nonsingular**), and $B$ is called the **inverse** of $A$, denoted $A^{-1}$. If $A$ is not invertible, it is called **singular**.

!!! theorem "Theorem 2.5 (Uniqueness of the inverse)"
    If a matrix $A$ is invertible, its inverse is unique.

??? proof "Proof"
    Suppose $B$ and $C$ are both inverses of $A$. Then

    $$
    B = BI = B(AC) = (BA)C = IC = C.
    $$

    Therefore $B = C$, and the inverse is unique. $\blacksquare$

!!! theorem "Theorem 2.6 (Properties of inverses)"
    Let $A, B$ be invertible matrices of order $n$, and $c$ a nonzero scalar. Then:

    1. $(A^{-1})^{-1} = A$
    2. $(AB)^{-1} = B^{-1}A^{-1}$ (note the order reversal)
    3. $(A^T)^{-1} = (A^{-1})^T$
    4. $(cA)^{-1} = \frac{1}{c}A^{-1}$

??? proof "Proof"
    Property 2 as an example:

    $(AB)(B^{-1}A^{-1}) = A(BB^{-1})A^{-1} = AIA^{-1} = AA^{-1} = I$

    Similarly, $(B^{-1}A^{-1})(AB) = I$, so $(AB)^{-1} = B^{-1}A^{-1}$. $\blacksquare$

!!! theorem "Theorem 2.7 (Equivalent conditions for invertibility)"
    Let $A$ be a square matrix of order $n$. The following conditions are equivalent:

    1. $A$ is invertible.
    2. $A$ is row equivalent to $I_n$ (the RREF of $A$ is $I_n$).
    3. $A$ has $n$ pivot positions.
    4. The homogeneous system $A\mathbf{x} = \mathbf{0}$ has only the trivial solution.
    5. For every $\mathbf{b} \in \mathbb{R}^n$, the system $A\mathbf{x} = \mathbf{b}$ has a unique solution.
    6. The columns of $A$ are linearly independent.
    7. The columns of $A$ span $\mathbb{R}^n$.
    8. $\det(A) \neq 0$.
    9. $\operatorname{rank}(A) = n$.

### Method 1 for finding inverses: Row reduction

Augment $A$ with $I_n$ to form the $n \times 2n$ augmented matrix $[A \mid I]$, then apply row operations to reduce the left half to $I_n$. If successful, the right half is $A^{-1}$:

$$
[A \mid I] \xrightarrow{\text{row operations}} [I \mid A^{-1}].
$$

!!! example "Example 2.4"
    Find the inverse of

    $$
    A = \begin{pmatrix} 1 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 3 \end{pmatrix}
    $$

    Form the augmented matrix and perform row operations:

    $$
    \left(\begin{array}{ccc|ccc} 1 & 2 & 1 & 1 & 0 & 0 \\ 2 & 5 & 3 & 0 & 1 & 0 \\ 1 & 3 & 3 & 0 & 0 & 1 \end{array}\right)
    \xrightarrow{\text{row operations}}
    \left(\begin{array}{ccc|ccc} 1 & 0 & 0 & 6 & -3 & 1 \\ 0 & 1 & 0 & -3 & 2 & -1 \\ 0 & 0 & 1 & 1 & -1 & 1 \end{array}\right)
    $$

    Therefore

    $$
    A^{-1} = \begin{pmatrix} 6 & -3 & 1 \\ -3 & 2 & -1 \\ 1 & -1 & 1 \end{pmatrix}.
    $$

### Method 2 for finding inverses: Adjugate matrix method

For an invertible matrix $A$ of order $n$, its inverse can also be computed using the adjugate matrix (discussed in detail in the chapter on determinants):

$$
A^{-1} = \frac{1}{\det(A)} \operatorname{adj}(A),
$$

where $\operatorname{adj}(A)$ is the **adjugate matrix** of $A$, i.e., the transpose of the cofactor matrix.

!!! example "Example 2.5"
    For a $2 \times 2$ matrix $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$, if $ad - bc \neq 0$, then

    $$
    A^{-1} = \frac{1}{ad - bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}.
    $$

---

## 2.6 Elementary Matrices

<div class="context-flow" markdown>

**Representing row operations as matrices**: The three types of row operations from Chapter 1 → perform one row operation on $I$ to get **elementary matrix** $E$ → left-multiplying by $E$ = performing one row operation → $A$ invertible $\Leftrightarrow$ $A$ is a product of elementary matrices

</div>

!!! definition "Definition 2.7 (Elementary matrix)"
    A matrix obtained by applying **one** elementary row operation to the identity matrix $I_n$ is called an **elementary matrix**. Corresponding to the three types of elementary row operations, there are three types:

    1. **Row interchange matrix** $E(i, j)$: Swap rows $i$ and $j$ of $I_n$.
    2. **Row scaling matrix** $E(i(c))$: Multiply row $i$ of $I_n$ by a nonzero constant $c$.
    3. **Row replacement matrix** $E(i, j(c))$: Add $c$ times row $j$ of $I_n$ to row $i$.

!!! theorem "Theorem 2.8 (Elementary matrices and row operations)"
    Applying one elementary row operation to an $m \times n$ matrix $A$ is equivalent to **left-multiplying** $A$ by the corresponding elementary matrix of order $m$.

    Similarly, applying one elementary **column** operation to $A$ is equivalent to **right-multiplying** $A$ by the corresponding elementary matrix of order $n$.

??? proof "Proof"
    We use row replacement as an example. Let $E = E(i, j(c))$. Then the $i$-th row of $E$ is $\mathbf{e}_i + c\mathbf{e}_j$ (where $\mathbf{e}_k$ is the $k$-th standard basis vector), and the other rows are the same as $I$.

    The $i$-th row of $EA$ $= (\mathbf{e}_i + c\mathbf{e}_j)A =$ the $i$-th row of $A$ $+ c \cdot$ the $j$-th row of $A$.

    The $k$-th row of $EA$ ($k \neq i$) $= \mathbf{e}_k A =$ the $k$-th row of $A$.

    This is exactly the effect of row replacement $R_i + cR_j \to R_i$. $\blacksquare$

!!! theorem "Theorem 2.9 (Elementary matrices are invertible)"
    Every elementary matrix is invertible, and its inverse is an elementary matrix of the same type:

    1. $E(i, j)^{-1} = E(i, j)$
    2. $E(i(c))^{-1} = E(i(1/c))$
    3. $E(i, j(c))^{-1} = E(i, j(-c))$

!!! theorem "Theorem 2.10 (Invertible matrices are products of elementary matrices)"
    A square matrix $A$ of order $n$ is invertible if and only if $A$ can be expressed as a product of finitely many elementary matrices.

??? proof "Proof"
    $(\Rightarrow)$ If $A$ is invertible, then $A$ is row equivalent to $I_n$, i.e., there exist elementary matrices $E_1, E_2, \ldots, E_k$ such that $E_k \cdots E_2 E_1 A = I$. Therefore $A = E_1^{-1} E_2^{-1} \cdots E_k^{-1}$, and each $E_i^{-1}$ is also an elementary matrix.

    $(\Leftarrow)$ Elementary matrices are all invertible, and the product of invertible matrices is invertible. $\blacksquare$

!!! example "Example 2.6"
    Express $A = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}$ as a product of elementary matrices.

    Reduce $A$ to $I$:

    $$
    \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix} \xrightarrow{R_2 - 3R_1} \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} \xrightarrow{R_1 - 2R_2} \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
    $$

    That is, $E(1,2(-2)) \cdot E(2,1(-3)) \cdot A = I$, so

    $$
    A = E(2,1(-3))^{-1} \cdot E(1,2(-2))^{-1} = E(2,1(3)) \cdot E(1,2(2)) = \begin{pmatrix} 1 & 0 \\ 3 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}.
    $$

---

## 2.7 Block Matrices

<div class="context-flow" markdown>

**Divide and conquer**: Operate on large matrices block by block → block diagonal matrices decompose determinant/inverse into independent subproblems → the block perspective pervades later chapters (Chapter 6 diagonalization, Chapter 5 invariant subspaces)

</div>

When matrices are large, partitioning them into smaller submatrices (**blocks**) often simplifies analysis and computation.

!!! definition "Definition 2.8 (Block matrix)"
    Partitioning a matrix $A$ with horizontal and vertical lines into several submatrices is called a **partitioning** of $A$. Each submatrix is called a **block** of $A$. For example, partitioning an $m \times n$ matrix $A$ as

    $$
    A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}
    $$

    where $A_{11}$ is $m_1 \times n_1$, $A_{12}$ is $m_1 \times n_2$, etc., with $m_1 + m_2 = m$, $n_1 + n_2 = n$.

!!! theorem "Theorem 2.11 (Block matrix multiplication)"
    If matrices $A$ and $B$ are partitioned in a compatible way (i.e., the column partition of $A$ matches the row partition of $B$), then block matrix multiplication can be performed as if the blocks were scalar entries.

    For example, if

    $$
    A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}, \quad
    B = \begin{pmatrix} B_{11} & B_{12} \\ B_{21} & B_{22} \end{pmatrix}
    $$

    are compatibly partitioned, then

    $$
    AB = \begin{pmatrix} A_{11}B_{11} + A_{12}B_{21} & A_{11}B_{12} + A_{12}B_{22} \\ A_{21}B_{11} + A_{22}B_{21} & A_{21}B_{12} + A_{22}B_{22} \end{pmatrix}.
    $$

!!! definition "Definition 2.9 (Block diagonal matrix)"
    A matrix of the form

    $$
    A = \begin{pmatrix} A_1 & & \\ & A_2 & \\ & & \ddots & \\ & & & A_k \end{pmatrix}
    $$

    is called a **block diagonal matrix**, where $A_1, A_2, \ldots, A_k$ are square matrices and all entries outside the diagonal blocks are zero. Denoted $A = \operatorname{diag}(A_1, A_2, \ldots, A_k)$.

!!! proposition "Proposition 2.2 (Properties of block diagonal matrices)"
    Let $A = \operatorname{diag}(A_1, A_2, \ldots, A_k)$. Then:

    1. $\det(A) = \det(A_1)\det(A_2)\cdots\det(A_k)$.
    2. $A$ is invertible if and only if each $A_i$ is invertible, in which case $A^{-1} = \operatorname{diag}(A_1^{-1}, A_2^{-1}, \ldots, A_k^{-1})$.

!!! example "Example 2.7"
    Let

    $$
    A = \begin{pmatrix} A_{11} & A_{12} \\ O & A_{22} \end{pmatrix}
    $$

    be a block upper triangular matrix, where $A_{11}$ and $A_{22}$ are square. If both $A_{11}$ and $A_{22}$ are invertible, then $A$ is invertible, and

    $$
    A^{-1} = \begin{pmatrix} A_{11}^{-1} & -A_{11}^{-1}A_{12}A_{22}^{-1} \\ O & A_{22}^{-1} \end{pmatrix}.
    $$

    Verification:

    $$
    \begin{pmatrix} A_{11} & A_{12} \\ O & A_{22} \end{pmatrix}\begin{pmatrix} A_{11}^{-1} & -A_{11}^{-1}A_{12}A_{22}^{-1} \\ O & A_{22}^{-1} \end{pmatrix} = \begin{pmatrix} I & A_{12}A_{22}^{-1} - A_{12}A_{22}^{-1} \\ O & I \end{pmatrix} = \begin{pmatrix} I & O \\ O & I \end{pmatrix}.
    $$

---

## 2.8 Rank of a Matrix

<div class="context-flow" markdown>

**Rank = the "effective dimension" of a matrix**: Number of nonzero rows in REF = column rank = row rank → $\operatorname{rank}(A)$ uniformly characterizes the degrees of freedom of a system → Sylvester's inequality controls the rank of products

</div>

!!! definition "Definition 2.10 (Rank)"
    The **rank** of a matrix $A$, denoted $\operatorname{rank}(A)$ or $r(A)$, is defined as the number of nonzero rows in the row echelon form of $A$, i.e., the number of pivots.

    Equivalently, $\operatorname{rank}(A)$ equals the size of a maximal linearly independent subset of the column vectors (**column rank**), which also equals the size of a maximal linearly independent subset of the row vectors (**row rank**).

!!! theorem "Theorem 2.12 (Row rank = column rank)"
    For any matrix $A$, its row rank equals its column rank.

??? proof "Proof"
    Let $A$ be an $m \times n$ matrix with $\operatorname{rank}(A) = r$. Reduce $A$ to row echelon form $R$, which has $r$ nonzero rows, so the row rank $\le r$. Since row operations do not change the row space (the space spanned by the rows), the row rank of $A$ equals the row rank of $R$, which equals $r$ (the $r$ nonzero rows of $R$ are linearly independent).

    On the other hand, row operations do not change the linear dependence relations among column vectors (since they do not change the solution set of $A\mathbf{x} = \mathbf{0}$), so $A$ and $R$ have the same column rank. $R$ has $r$ pivot columns, corresponding to $r$ linearly independent column vectors, so the column rank also equals $r$. $\blacksquare$

!!! theorem "Theorem 2.13 (Properties of rank)"
    Let $A$ be an $m \times n$ matrix and $B$ an $n \times p$ matrix. Then:

    1. $0 \le \operatorname{rank}(A) \le \min(m, n)$.
    2. $\operatorname{rank}(A) = \operatorname{rank}(A^T)$.
    3. $\operatorname{rank}(AB) \le \min(\operatorname{rank}(A), \operatorname{rank}(B))$.
    4. If $P$ is an invertible matrix of order $m$ and $Q$ is an invertible matrix of order $n$, then $\operatorname{rank}(PAQ) = \operatorname{rank}(A)$.
    5. $\operatorname{rank}(A + B) \le \operatorname{rank}(A) + \operatorname{rank}(B)$ (when $A, B$ have the same size).

!!! theorem "Theorem 2.14 (Sylvester's rank inequality)"
    Let $A$ be an $m \times n$ matrix and $B$ an $n \times p$ matrix. Then

    $$
    \operatorname{rank}(A) + \operatorname{rank}(B) - n \le \operatorname{rank}(AB).
    $$

??? proof "Proof"
    Consider the column space of $B$ and the null space of $A$. The column space of $B$ has dimension $\operatorname{rank}(B)$, and the null space of $A$ has dimension $n - \operatorname{rank}(A)$.

    The column space of $AB$ is obtained by the action of $A$ on the column space of $B$. The vectors in the column space of $B$ that are mapped to zero by $A$ are exactly the intersection of the column space of $B$ with the null space of $A$.

    By the dimension formula:

    $$
    \operatorname{rank}(AB) = \dim(A(\operatorname{Col}(B))) = \operatorname{rank}(B) - \dim(\operatorname{Col}(B) \cap \ker(A))
    $$

    $$
    \ge \operatorname{rank}(B) - \dim(\ker(A)) = \operatorname{rank}(B) - (n - \operatorname{rank}(A))
    $$

    $$
    = \operatorname{rank}(A) + \operatorname{rank}(B) - n. \quad \blacksquare
    $$

!!! example "Example 2.8"
    Find the rank of

    $$
    A = \begin{pmatrix} 1 & 2 & 3 & 4 \\ 2 & 4 & 7 & 9 \\ 1 & 2 & 4 & 5 \end{pmatrix}
    $$

    Row operations:

    $$
    \begin{pmatrix} 1 & 2 & 3 & 4 \\ 2 & 4 & 7 & 9 \\ 1 & 2 & 4 & 5 \end{pmatrix}
    \xrightarrow{R_2 - 2R_1, \; R_3 - R_1}
    \begin{pmatrix} 1 & 2 & 3 & 4 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 1 & 1 \end{pmatrix}
    \xrightarrow{R_3 - R_2}
    \begin{pmatrix} 1 & 2 & 3 & 4 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 0 \end{pmatrix}
    $$

    There are $2$ nonzero rows, so $\operatorname{rank}(A) = 2$.

!!! example "Example 2.9"
    Let $A$ be a $3 \times 5$ matrix and $B$ a $5 \times 4$ matrix, with $\operatorname{rank}(A) = 3$ and $\operatorname{rank}(B) = 4$. By Sylvester's inequality:

    $$
    \operatorname{rank}(AB) \ge 3 + 4 - 5 = 2.
    $$

    Also $\operatorname{rank}(AB) \le \min(3, 4) = 3$, so $2 \le \operatorname{rank}(AB) \le 3$.
