# Chapter 19  Kronecker Product and Vec Operator

<div class="context-flow" markdown>

**Prerequisites**: Matrix multiplication / eigenvalues (Ch6)

**Chapter arc**: $A\otimes B$ combines two spaces → **Vec operator** transforms matrix equation $AXB=C$ into $(B^T\otimes A)\operatorname{vec}(X)=\operatorname{vec}(C)$ → Ch20 Sylvester/Lyapunov solving

**Further connections**：Kronecker products are indispensable in quantum computing (multi-qubit systems), signal processing (MIMO systems), and statistics (vectorized covariance estimation); the Vec operator converts matrix equations to vector equations and is the standard tool for computing matrix derivatives

</div>

In practical applications of linear algebra, we frequently encounter situations where matrix equations need to be converted into vector equations, or where large matrices with special structure need to be constructed. The **Kronecker product** and the **Vec operator** (vectorization operator) are precisely the core tools for handling such problems. The Kronecker product provides a systematic way to "combine" linear maps on two matrix spaces, while the Vec operator "flattens" a matrix into a vector, enabling matrix equations to be converted into standard systems of linear equations via the Kronecker product.

Starting from the definition and basic properties of the Kronecker product, this chapter introduces the Vec operator and its core formula, discusses the role of the commutation matrix, then demonstrates how to use these tools to solve matrix equations, and finally introduces the Kronecker sum and its relationship with the matrix exponential.

---

## 19.1 Definition of the Kronecker product

<div class="context-flow" markdown>


</div>

!!! definition "Definition 19.1 (Kronecker product)"
    Let $A = (a_{ij})$ be an $m \times n$ matrix and $B$ be a $p \times q$ matrix. The **Kronecker product** (also called the **tensor product**) of $A$ and $B$, denoted $A \otimes B$, is defined as the $mp \times nq$ **block matrix**:

    $$
    A \otimes B = \begin{pmatrix}
    a_{11}B & a_{12}B & \cdots & a_{1n}B \\
    a_{21}B & a_{22}B & \cdots & a_{2n}B \\
    \vdots & \vdots & \ddots & \vdots \\
    a_{m1}B & a_{m2}B & \cdots & a_{mn}B
    \end{pmatrix}.
    $$

!!! note "Note"
    The Kronecker product is not commutative: in general $A \otimes B \neq B \otimes A$, but they are **permutation similar**, i.e., there exists a permutation matrix $P$ such that $B \otimes A = P(A \otimes B)P^T$. This permutation matrix is the commutation matrix discussed later.

!!! example "Example 19.1"
    Let $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$, $B = \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix}$. Then:

    $$
    A \otimes B = \begin{pmatrix}
    1 \cdot \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix} & 2 \cdot \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix} \\[6pt]
    3 \cdot \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix} & 4 \cdot \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix}
    \end{pmatrix}
    = \begin{pmatrix}
    0 & 5 & 0 & 10 \\
    6 & 7 & 12 & 14 \\
    0 & 15 & 0 & 20 \\
    18 & 21 & 24 & 28
    \end{pmatrix}.
    $$

---

## 19.2 Properties of the Kronecker product

<div class="context-flow" markdown>

**Core property**: **Mixed-product** $(A\otimes B)(C\otimes D)=(AC)\otimes(BD)$ → implies formulas for inverse, transpose, determinant · $\det(A\otimes B)=(\det A)^n(\det B)^m$

</div>

The Kronecker product possesses rich and elegant algebraic properties that make it an indispensable tool in matrix theory.

!!! theorem "Theorem 19.1 (Mixed-product property)"
    Let $A, C$ be a compatible pair of matrices for multiplication, and $B, D$ be a compatible pair for multiplication. Then:

    $$
    (A \otimes B)(C \otimes D) = (AC) \otimes (BD).
    $$

    This property is called the **mixed-product property**.

??? proof "Proof"
    Let $A$ be $m \times n$, $B$ be $p \times q$, $C$ be $n \times r$, $D$ be $q \times s$. Then $A \otimes B$ is $mp \times nq$, $C \otimes D$ is $nq \times rs$, and the product $(A \otimes B)(C \otimes D)$ is $mp \times rs$.

    The $(i,j)$-block (of size $p \times q$) of $(A \otimes B)$ is $a_{ij}B$, and the $(j,k)$-block (of size $q \times s$) of $(C \otimes D)$ is $c_{jk}D$.

    The $(i,k)$-block of the product is:

    $$
    \sum_{j=1}^{n} (a_{ij}B)(c_{jk}D) = \sum_{j=1}^{n} a_{ij}c_{jk}(BD) = \left(\sum_{j=1}^n a_{ij}c_{jk}\right)(BD) = (AC)_{ik}(BD).
    $$

    This is precisely the $(i,k)$-block of $(AC) \otimes (BD)$. Therefore $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$. $\blacksquare$

!!! theorem "Theorem 19.2 (Basic algebraic properties of the Kronecker product)"
    Let $A, B, C$ be matrices of appropriate sizes and $\alpha$ be a scalar. Then:

    1. **Associativity**: $(A \otimes B) \otimes C = A \otimes (B \otimes C)$.
    2. **Distributivity**: $A \otimes (B + C) = A \otimes B + A \otimes C$, $(A + B) \otimes C = A \otimes C + B \otimes C$.
    3. **Scalar multiplication**: $(\alpha A) \otimes B = A \otimes (\alpha B) = \alpha(A \otimes B)$.
    4. **Transpose**: $(A \otimes B)^T = A^T \otimes B^T$.
    5. **Conjugate transpose**: $(A \otimes B)^* = A^* \otimes B^*$.
    6. **Inverse**: If $A, B$ are both invertible, then $(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$.

??? proof "Proof"
    We prove property 6. By the mixed-product property:

    $$
    (A \otimes B)(A^{-1} \otimes B^{-1}) = (AA^{-1}) \otimes (BB^{-1}) = I_m \otimes I_p = I_{mp}.
    $$

    Similarly $(A^{-1} \otimes B^{-1})(A \otimes B) = I_{mp}$. Therefore $(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$.

    The other properties can be verified directly from the definition. $\blacksquare$

!!! theorem "Theorem 19.3 (Trace, determinant, and rank of the Kronecker product)"
    Let $A$ be an $m \times m$ matrix and $B$ be an $n \times n$ matrix. Then:

    1. **Trace**: $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \cdot \operatorname{tr}(B)$.
    2. **Determinant**: $\det(A \otimes B) = (\det A)^n (\det B)^m$.
    3. **Rank**: $\operatorname{rank}(A \otimes B) = \operatorname{rank}(A) \cdot \operatorname{rank}(B)$.

??? proof "Proof"
    **(1) Trace**: The diagonal blocks of $A \otimes B$ are $a_{ii}B$ ($i = 1,\ldots,m$), therefore:

    $$
    \operatorname{tr}(A \otimes B) = \sum_{i=1}^m \operatorname{tr}(a_{ii}B) = \sum_{i=1}^m a_{ii} \operatorname{tr}(B) = \operatorname{tr}(A) \cdot \operatorname{tr}(B).
    $$

    **(2) Determinant**: Using the mixed-product property and block diagonalization. Let the eigenvalues of $A$ be $\lambda_1, \ldots, \lambda_m$ and those of $B$ be $\mu_1, \ldots, \mu_n$ (counted with multiplicity). By Theorem 19.5 (to be proved later), the eigenvalues of $A \otimes B$ are $\{\lambda_i \mu_j\}$. Therefore:

    $$
    \det(A \otimes B) = \prod_{i=1}^m \prod_{j=1}^n \lambda_i \mu_j = \left(\prod_{i=1}^m \lambda_i\right)^n \left(\prod_{j=1}^n \mu_j\right)^m = (\det A)^n (\det B)^m.
    $$

    **(3) Rank**: Let $\operatorname{rank}(A) = r$ and $\operatorname{rank}(B) = s$. There exist invertible matrices $P, Q, U, V$ such that $A = P \begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} Q$, $B = U \begin{pmatrix} I_s & 0 \\ 0 & 0 \end{pmatrix} V$.

    Then $A \otimes B = (P \otimes U) \left(\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} \otimes \begin{pmatrix} I_s & 0 \\ 0 & 0 \end{pmatrix}\right) (Q \otimes V)$.

    By direct computation, the middle matrix has rank $rs$, and $P \otimes U$ and $Q \otimes V$ are invertible, so $\operatorname{rank}(A \otimes B) = rs$. $\blacksquare$

!!! example "Example 19.2"
    Let $A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}$, $B = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$.

    **Trace**: $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \cdot \operatorname{tr}(B) = 5 \times 3 = 15$.

    Direct verification: the diagonal entries of $A \otimes B$ are $2 \cdot 1, 2 \cdot 2, 3 \cdot 1, 3 \cdot 2 = 2, 4, 3, 6$, and the trace is $2+4+3+6=15$.

    **Determinant**: $\det(A \otimes B) = (\det A)^2 (\det B)^2 = 6^2 \cdot 2^2 = 36 \cdot 4 = 144$.

    **Rank**: $\operatorname{rank}(A \otimes B) = 2 \times 2 = 4$ (full rank).

!!! definition "Definition 19.2 (Kronecker power)"
    For a square matrix $A$, the $k$-th **Kronecker power** is defined as:

    $$
    A^{\otimes k} = \underbrace{A \otimes A \otimes \cdots \otimes A}_{k \text{ copies}}.
    $$

    If $A$ is $n \times n$, then $A^{\otimes k}$ is $n^k \times n^k$.

---

## 19.3 Vec operator

<div class="context-flow" markdown>

**Bridge formula**: $\operatorname{vec}(AXB)=(B^T\otimes A)\operatorname{vec}(X)$ · Special cases: $\operatorname{vec}(AX)=(I\otimes A)\operatorname{vec}(X)$, $\operatorname{vec}(XB)=(B^T\otimes I)\operatorname{vec}(X)$

</div>

The Vec operator stacks the columns of a matrix into a vector, serving as the bridge connecting matrix equations and vector equations.

!!! definition "Definition 19.3 (Vec operator)"
    Let $A = (\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n)$ be an $m \times n$ matrix, where $\mathbf{a}_j$ is the $j$-th column. The **Vec operator** (vectorization) maps $A$ to an $mn \times 1$ column vector:

    $$
    \operatorname{vec}(A) = \begin{pmatrix} \mathbf{a}_1 \\ \mathbf{a}_2 \\ \vdots \\ \mathbf{a}_n \end{pmatrix}.
    $$

    That is, the columns of $A$ are stacked from left to right.

!!! theorem "Theorem 19.4 (Core formula of the Vec operator)"
    Let $A$ be an $m \times n$ matrix, $X$ be an $n \times p$ matrix, and $B$ be a $p \times q$ matrix. Then:

    $$
    \operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X).
    $$

??? proof "Proof"
    **Method 1** (using column vectors): Let $B = (\mathbf{b}_1, \ldots, \mathbf{b}_q)$, $Y = AXB$. The $j$-th column of $Y$ is:

    $$
    \mathbf{y}_j = AX\mathbf{b}_j = A \sum_{k=1}^p b_{kj} \mathbf{x}_k = \sum_{k=1}^p b_{kj} A \mathbf{x}_k,
    $$

    where $\mathbf{x}_k$ is the $k$-th column of $X$.

    Therefore:

    $$
    \operatorname{vec}(Y) = \begin{pmatrix} \mathbf{y}_1 \\ \vdots \\ \mathbf{y}_q \end{pmatrix} = \begin{pmatrix} \sum_k b_{k1} A\mathbf{x}_k \\ \vdots \\ \sum_k b_{kq} A\mathbf{x}_k \end{pmatrix}.
    $$

    On the other hand:

    $$
    (B^T \otimes A)\operatorname{vec}(X) = \begin{pmatrix} b_{11}A & b_{21}A & \cdots & b_{p1}A \\ b_{12}A & b_{22}A & \cdots & b_{p2}A \\ \vdots & & \ddots & \vdots \\ b_{1q}A & b_{2q}A & \cdots & b_{pq}A \end{pmatrix} \begin{pmatrix} \mathbf{x}_1 \\ \mathbf{x}_2 \\ \vdots \\ \mathbf{x}_p \end{pmatrix}.
    $$

    The $j$-th block is $\sum_{k=1}^p b_{kj} A \mathbf{x}_k = \mathbf{y}_j$. The two expressions agree. $\blacksquare$

    **Method 2** (using basis matrices): Let $E_{ij}$ be the matrix with $1$ in position $(i,j)$ and $0$ elsewhere. By linearity of $\operatorname{vec}(AE_{ij}B)$, it suffices to verify for $X = E_{ij}$. The $(r,s)$-entry of $AE_{ij}B$ is $a_{ri}b_{js}$, which matches the corresponding entry of $(B^T \otimes A)$.

!!! theorem "Theorem 19.5 (Special cases of the Vec operator)"
    The following are commonly used special cases of Theorem 19.4:

    1. $\operatorname{vec}(AX) = (I \otimes A)\operatorname{vec}(X)$ (taking $B = I$);
    2. $\operatorname{vec}(XB) = (B^T \otimes I)\operatorname{vec}(X)$ (taking $A = I$);
    3. $\operatorname{vec}(A\mathbf{x}\mathbf{b}^T) = (\mathbf{b} \otimes A)\mathbf{x}$ (when $X$ is a column vector);
    4. $\operatorname{vec}(\alpha A) = \alpha \operatorname{vec}(A)$.

??? proof "Proof"
    All are direct corollaries of Theorem 19.4, obtained by taking the corresponding $A$, $B$, or $X$ to be special matrices (identity matrix, vectors, etc.). For example:

    (1) Taking $B = I_p$, we get $\operatorname{vec}(AXI) = (I^T \otimes A)\operatorname{vec}(X) = (I \otimes A)\operatorname{vec}(X)$.

    (2) Taking $A = I_m$, we get $\operatorname{vec}(IXB) = (B^T \otimes I)\operatorname{vec}(X)$. $\blacksquare$

!!! definition "Definition 19.4 (Half-vectorization operator)"
    For an $n \times n$ symmetric matrix $A$, the **half-vectorization operator** $\operatorname{vech}(A)$ stacks the lower triangular part (including the diagonal) of $A$ column by column into an $\frac{n(n+1)}{2} \times 1$ vector:

    $$
    \operatorname{vech}(A) = (a_{11}, a_{21}, \ldots, a_{n1}, a_{22}, a_{32}, \ldots, a_{nn})^T.
    $$

!!! example "Example 19.3"
    Let $X = \begin{pmatrix} x_{11} & x_{12} \\ x_{21} & x_{22} \end{pmatrix}$, $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$.

    $\operatorname{vec}(X) = (x_{11}, x_{21}, x_{12}, x_{22})^T$.

    $AX = \begin{pmatrix} ax_{11}+bx_{21} & ax_{12}+bx_{22} \\ cx_{11}+dx_{21} & cx_{12}+dx_{22} \end{pmatrix}$.

    $\operatorname{vec}(AX) = (ax_{11}+bx_{21},\; cx_{11}+dx_{21},\; ax_{12}+bx_{22},\; cx_{12}+dx_{22})^T$.

    $(I_2 \otimes A)\operatorname{vec}(X) = \begin{pmatrix} A & 0 \\ 0 & A \end{pmatrix}\begin{pmatrix} x_{11} \\ x_{21} \\ x_{12} \\ x_{22} \end{pmatrix} = \begin{pmatrix} ax_{11}+bx_{21} \\ cx_{11}+dx_{21} \\ ax_{12}+bx_{22} \\ cx_{12}+dx_{22} \end{pmatrix}$.

    The two are equal, verifying $\operatorname{vec}(AX) = (I \otimes A)\operatorname{vec}(X)$.

!!! example "Example 19.4"
    Using the Vec formula to solve $AXB = C$.

    Let $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$, $B = \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix}$, $C = \begin{pmatrix} 6 & 2 \\ 0 & 8 \end{pmatrix}$.

    Vectorization: $(B^T \otimes A)\operatorname{vec}(X) = \operatorname{vec}(C)$.

    $B^T \otimes A = \begin{pmatrix} 3A & 0 \\ 0 & A \end{pmatrix} = \begin{pmatrix} 3 & 0 & 0 & 0 \\ 0 & 6 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 2 \end{pmatrix}$.

    $\operatorname{vec}(C) = (6, 0, 2, 8)^T$.

    Solution: $\operatorname{vec}(X) = (B^T \otimes A)^{-1}\operatorname{vec}(C) = (2, 0, 2, 4)^T$.

    Therefore $X = \begin{pmatrix} 2 & 2 \\ 0 & 4 \end{pmatrix}$.

    Verification: $AXB = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}\begin{pmatrix} 2 & 2 \\ 0 & 4 \end{pmatrix}\begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 2 \\ 0 & 8 \end{pmatrix}\begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 6 & 2 \\ 0 & 8 \end{pmatrix} = C$. Correct.

---

## 19.4 Commutation matrix

<div class="context-flow" markdown>

**Role**: $K_{m,n}\operatorname{vec}(A)=\operatorname{vec}(A^T)$ · Realizes the index rearrangement $A\otimes B \leftrightarrow B\otimes A$ · $K_{n,n}^2=I$ (involution)

</div>

The Kronecker product is not commutative, but the two orderings of the Kronecker product are related through a special permutation matrix.

!!! definition "Definition 19.5 (Commutation matrix)"
    The **commutation matrix** $K_{m,n}$ is an $mn \times mn$ permutation matrix satisfying, for any $m \times n$ matrix $A$:

    $$
    K_{m,n} \operatorname{vec}(A) = \operatorname{vec}(A^T).
    $$

    Equivalently, $K_{m,n}$ can be expressed using basis matrices as:

    $$
    K_{m,n} = \sum_{i=1}^{m} \sum_{j=1}^{n} E_{ij} \otimes E_{ji},
    $$

    where $E_{ij}$ is the $m \times n$ basis matrix (with $1$ in position $(i,j)$ and $0$ elsewhere), and $E_{ji}$ is the $n \times m$ basis matrix.

!!! theorem "Theorem 19.6 (Properties of the commutation matrix)"
    The commutation matrix $K_{m,n}$ has the following properties:

    1. $K_{m,n}^T = K_{m,n}^{-1} = K_{n,m}$.
    2. $K_{m,n}(A \otimes B)K_{p,q} = B \otimes A$ (under appropriate sizes). In particular, $K_{m,n}(A \otimes B) = (B \otimes A)K_{p,q}$.
    3. $K_{n,n}$ is symmetric and orthogonal, $K_{n,n}^2 = I_{n^2}$ (involution).
    4. $K_{1,n} = K_{n,1} = I_n$.
    5. $(A \otimes B) = K_{p,m}(B \otimes A)K_{n,q}$, where $A$ is $m \times n$ and $B$ is $p \times q$.

??? proof "Proof"
    **(1)**: For any $m \times n$ matrix $A$, $K_{m,n}\operatorname{vec}(A) = \operatorname{vec}(A^T)$. Replacing $A$ by $A^T$ ($n \times m$): $K_{n,m}\operatorname{vec}(A^T) = \operatorname{vec}(A)$. Therefore $K_{n,m}K_{m,n}\operatorname{vec}(A) = \operatorname{vec}(A)$ for all $A$, so $K_{n,m}K_{m,n} = I_{mn}$, i.e., $K_{m,n}^{-1} = K_{n,m}$.

    Since $K_{m,n}$ is a permutation matrix, $K_{m,n}^T = K_{m,n}^{-1} = K_{n,m}$.

    **(2)**: For $A$ of size $m \times n$, $B$ of size $p \times q$: consider the action of $(A \otimes B)$ on $\operatorname{vec}(X)$ (adjusting sizes as needed).

    By direct verification: $K_{m,n}$ rearranges elements from column-major to row-major order in $\operatorname{vec}$. The $((i-1)p+k, (j-1)q+l)$-entry of $(A \otimes B)$ is $a_{ij}b_{kl}$, and the $((k-1)m+i, (l-1)n+j)$-entry of $(B \otimes A)$ is also $b_{kl}a_{ij}$. The commutation matrix precisely realizes the index rearrangement $((i-1)p+k) \leftrightarrow ((k-1)m+i)$. $\blacksquare$

!!! theorem "Theorem 19.7 (Relationship between Vec and transpose)"
    For any $m \times n$ matrix $A$:

    $$
    \operatorname{vec}(A^T) = K_{m,n}\operatorname{vec}(A).
    $$

    Furthermore, for the transpose of a matrix product:

    $$
    \operatorname{vec}((AXB)^T) = \operatorname{vec}(B^T X^T A^T) = (A \otimes B^T)\operatorname{vec}(X^T) = (A \otimes B^T)K_{n,p}\operatorname{vec}(X).
    $$

??? proof "Proof"
    The first equality is the definition of $K_{m,n}$. The second equality uses $\operatorname{vec}(B^T X^T A^T) = (A \otimes B^T)\operatorname{vec}(X^T)$ (by Theorem 19.4), and then $\operatorname{vec}(X^T) = K_{n,p}\operatorname{vec}(X)$. $\blacksquare$

!!! example "Example 19.5"
    Construct $K_{2,3}$. We need a $6 \times 6$ permutation matrix such that for any $2 \times 3$ matrix $A = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \end{pmatrix}$:

    $\operatorname{vec}(A) = (a_{11}, a_{21}, a_{12}, a_{22}, a_{13}, a_{23})^T$.

    $\operatorname{vec}(A^T) = (a_{11}, a_{12}, a_{13}, a_{21}, a_{22}, a_{23})^T$.

    Therefore $K_{2,3}$ maps positions $(1,2,3,4,5,6)$ to $(1,3,5,2,4,6)$:

    $$
    K_{2,3} = \begin{pmatrix}
    1 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 1 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 1 & 0 \\
    0 & 1 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 1 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 1
    \end{pmatrix}.
    $$

---

## 19.5 Kronecker product and matrix equations

<div class="context-flow" markdown>

**Application**: $\sum A_kXB_k=C$ → $(\sum B_k^T\otimes A_k)\operatorname{vec}(X)=\operatorname{vec}(C)$

**Sylvester** $AX+XB=C$: solvable $\Leftrightarrow$ $\sigma(A)\cap\sigma(-B)=\emptyset$ → Ch20

</div>

One of the most important applications of the Kronecker product and Vec operator is converting matrix equations into standard systems of linear equations.

!!! definition "Definition 19.6 (Linear matrix equation)"
    An equation of the form

    $$
    \sum_{k=1}^{K} A_k X B_k = C
    $$

    is called a **linear matrix equation**, where $A_k, B_k, C$ are known and $X$ is the unknown matrix.

!!! theorem "Theorem 19.8 (Vectorization of matrix equations)"
    The linear matrix equation $\sum_{k=1}^K A_k X B_k = C$ is equivalent to the vector equation:

    $$
    \left(\sum_{k=1}^{K} B_k^T \otimes A_k\right) \operatorname{vec}(X) = \operatorname{vec}(C).
    $$

??? proof "Proof"
    Applying Theorem 19.4 to each term $A_k X B_k$:

    $$
    \operatorname{vec}(A_k X B_k) = (B_k^T \otimes A_k)\operatorname{vec}(X).
    $$

    Taking Vec of both sides of the equation:

    $$
    \operatorname{vec}\left(\sum_{k=1}^K A_k X B_k\right) = \sum_{k=1}^K \operatorname{vec}(A_k X B_k) = \sum_{k=1}^K (B_k^T \otimes A_k)\operatorname{vec}(X) = \left(\sum_{k=1}^K B_k^T \otimes A_k\right)\operatorname{vec}(X).
    $$

    The right side is $\operatorname{vec}(C)$. $\blacksquare$

!!! theorem "Theorem 19.9 (Kronecker product form of the Sylvester equation)"
    The Sylvester equation $AX + XB = C$ (where $A$ is $m \times m$, $B$ is $n \times n$, and $X, C$ are $m \times n$) is equivalent to:

    $$
    (I_n \otimes A + B^T \otimes I_m)\operatorname{vec}(X) = \operatorname{vec}(C).
    $$

    This equation has a unique solution if and only if $I_n \otimes A + B^T \otimes I_m$ is nonsingular, i.e., $A$ and $-B$ have no common eigenvalues.

??? proof "Proof"
    By Theorem 19.8, taking $K = 2$, $A_1 = A$, $B_1 = I$, $A_2 = I$, $B_2 = B$:

    $$
    (I^T \otimes A + B^T \otimes I)\operatorname{vec}(X) = (I_n \otimes A + B^T \otimes I_m)\operatorname{vec}(X) = \operatorname{vec}(C).
    $$

    The eigenvalues of $I_n \otimes A + B^T \otimes I_m$ are $\{\lambda_i(A) + \lambda_j(B) : i = 1,\ldots,m;\; j=1,\ldots,n\}$ (see Theorem 19.12), so it is nonsingular if and only if $\lambda_i(A) + \lambda_j(B) \neq 0$ for all $i, j$, i.e., $A$ and $-B$ have no common eigenvalues. $\blacksquare$

!!! example "Example 19.6"
    Solve the Sylvester equation $AX + XB = C$, where:

    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}, \quad B = (3), \quad C = \begin{pmatrix} 4 \\ 10 \end{pmatrix}.
    $$

    Here $A$ is $2 \times 2$, $B$ is $1 \times 1$ (the scalar 3), and $X$ is $2 \times 1$.

    Vectorization: $(I_1 \otimes A + B^T \otimes I_2)\operatorname{vec}(X) = \operatorname{vec}(C)$.

    $I_1 \otimes A = A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$, $B^T \otimes I_2 = 3I_2 = \begin{pmatrix} 3 & 0 \\ 0 & 3 \end{pmatrix}$.

    Coefficient matrix: $\begin{pmatrix} 4 & 0 \\ 0 & 5 \end{pmatrix}$.

    $\operatorname{vec}(X) = \begin{pmatrix} 4 & 0 \\ 0 & 5 \end{pmatrix}^{-1}\begin{pmatrix} 4 \\ 10 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$.

    Therefore $X = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$.

    Verification: $AX + XB = \begin{pmatrix} 1 \\ 4 \end{pmatrix} + \begin{pmatrix} 3 \\ 6 \end{pmatrix} = \begin{pmatrix} 4 \\ 10 \end{pmatrix} = C$. Correct.

!!! example "Example 19.7"
    Using the Kronecker product to determine solvability of a matrix equation.

    Consider $AX - XA = C$ (where $A$ is $n \times n$). Vectorization gives $(I \otimes A - A^T \otimes I)\operatorname{vec}(X) = \operatorname{vec}(C)$.

    If $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$, the eigenvalues of the coefficient matrix are $\{\lambda_i - \lambda_j\} = \{0, -1, 1, 0\}$.

    Since zero eigenvalues exist ($\lambda_1 - \lambda_1 = 0$ and $\lambda_2 - \lambda_2 = 0$), the coefficient matrix is singular, so the equation $AX - XA = C$ does not have a solution for all $C$.

    Specifically, $C$ must satisfy $\operatorname{tr}(C) = 0$ (since $\operatorname{tr}(AX - XA) = 0$), and having zero diagonal entries is also a necessary condition (for diagonal $A$).

---

## 19.6 Eigenvalue decomposition of the Kronecker product

<div class="context-flow" markdown>

**Intuition**: $\sigma(A\otimes B)=\{\lambda_i\mu_j\}$, eigenvectors $\mathbf{u}_i\otimes\mathbf{v}_j$ · SVD also multiplicative: $\sigma_k(A\otimes B)=\{\sigma_i(A)\sigma_j(B)\}$

</div>

!!! definition "Definition 19.7 (Spectrum of the Kronecker product)"
    Let $A$ be an $m \times m$ matrix with eigenvalues $\lambda_1, \ldots, \lambda_m$, and $B$ be an $n \times n$ matrix with eigenvalues $\mu_1, \ldots, \mu_n$. Then the $mn$ eigenvalues of $A \otimes B$ are:

    $$
    \sigma(A \otimes B) = \{\lambda_i \mu_j : i = 1,\ldots,m;\; j = 1,\ldots,n\}.
    $$

!!! theorem "Theorem 19.10 (Eigenvalues and eigenvectors of the Kronecker product)"
    If $A\mathbf{u} = \lambda\mathbf{u}$ and $B\mathbf{v} = \mu\mathbf{v}$, then:

    $$
    (A \otimes B)(\mathbf{u} \otimes \mathbf{v}) = \lambda\mu(\mathbf{u} \otimes \mathbf{v}).
    $$

    That is, $\mathbf{u} \otimes \mathbf{v}$ is an eigenvector of $A \otimes B$ corresponding to eigenvalue $\lambda\mu$.

    If $A$ and $B$ are both diagonalizable, $A = P \operatorname{diag}(\lambda_1,\ldots,\lambda_m) P^{-1}$, $B = Q \operatorname{diag}(\mu_1,\ldots,\mu_n) Q^{-1}$, then:

    $$
    A \otimes B = (P \otimes Q) \operatorname{diag}(\lambda_1\mu_1, \lambda_1\mu_2, \ldots, \lambda_m\mu_n) (P \otimes Q)^{-1}.
    $$

??? proof "Proof"
    By the mixed-product property:

    $$
    (A \otimes B)(\mathbf{u} \otimes \mathbf{v}) = (A\mathbf{u}) \otimes (B\mathbf{v}) = (\lambda\mathbf{u}) \otimes (\mu\mathbf{v}) = \lambda\mu(\mathbf{u} \otimes \mathbf{v}).
    $$

    For the diagonalizable case:

    $$
    A \otimes B = (P \Lambda_A P^{-1}) \otimes (Q \Lambda_B Q^{-1}) = (P \otimes Q)(\Lambda_A \otimes \Lambda_B)(P^{-1} \otimes Q^{-1}).
    $$

    Since $(P \otimes Q)^{-1} = P^{-1} \otimes Q^{-1}$, and $\Lambda_A \otimes \Lambda_B$ is a diagonal matrix (with diagonal entries $\lambda_i \mu_j$), the conclusion follows. $\blacksquare$

!!! theorem "Theorem 19.11 (Singular value decomposition of the Kronecker product)"
    Let $A = U_A \Sigma_A V_A^*$ and $B = U_B \Sigma_B V_B^*$ be the singular value decompositions of $A$ and $B$. Then:

    $$
    A \otimes B = (U_A \otimes U_B)(\Sigma_A \otimes \Sigma_B)(V_A \otimes V_B)^*.
    $$

    Therefore the singular values of $A \otimes B$ are $\{\sigma_i(A)\sigma_j(B)\}$.

??? proof "Proof"
    By the mixed-product property:

    $$
    A \otimes B = (U_A \Sigma_A V_A^*) \otimes (U_B \Sigma_B V_B^*) = (U_A \otimes U_B)(\Sigma_A \otimes \Sigma_B)(V_A^* \otimes V_B^*).
    $$

    Since $(V_A^* \otimes V_B^*) = (V_A \otimes V_B)^*$, and $U_A \otimes U_B$ and $V_A \otimes V_B$ are both unitary (the Kronecker product of unitary matrices is unitary), and $\Sigma_A \otimes \Sigma_B$ is a nonnegative diagonal matrix, this is precisely the singular value decomposition. $\blacksquare$

!!! example "Example 19.8"
    Let $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$, $B = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix}$.

    Eigenvalues of $A$: $\lambda_1 = 2$, $\lambda_2 = 3$. Eigenvalues of $B$: $\mu_1 = 1$, $\mu_2 = 4$.

    Eigenvalues of $A \otimes B$: $\{2 \times 1, 2 \times 4, 3 \times 1, 3 \times 4\} = \{2, 8, 3, 12\}$.

    Direct computation: $A \otimes B = \begin{pmatrix} 2 & 0 & 0 & 0 \\ 0 & 8 & 0 & 0 \\ 0 & 0 & 3 & 0 \\ 0 & 0 & 0 & 12 \end{pmatrix}$, a diagonal matrix with eigenvalues $2, 8, 3, 12$. Verification successful.

---

## 19.7 Kronecker sum

<div class="context-flow" markdown>

**Chapter arc**: $A\oplus B = A\otimes I+I\otimes B$, eigenvalues $\lambda_i+\mu_j$

**Key identity**: $e^{A\oplus B}=e^A\otimes e^B$ (because $A\otimes I$ and $I\otimes B$ commute) → Ch20 Lyapunov equation

</div>

The Kronecker sum is another operation closely related to the Kronecker product, with deep connections to the matrix exponential and Lyapunov equations.

!!! definition "Definition 19.8 (Kronecker sum)"
    Let $A$ be an $m \times m$ matrix and $B$ be an $n \times n$ matrix. The **Kronecker sum** of $A$ and $B$ is defined as:

    $$
    A \oplus B = A \otimes I_n + I_m \otimes B.
    $$

    It is an $mn \times mn$ matrix.

!!! theorem "Theorem 19.12 (Eigenvalues of the Kronecker sum)"
    Let the eigenvalues of $A$ be $\lambda_1, \ldots, \lambda_m$ and those of $B$ be $\mu_1, \ldots, \mu_n$. Then the $mn$ eigenvalues of $A \oplus B$ are:

    $$
    \sigma(A \oplus B) = \{\lambda_i + \mu_j : i = 1,\ldots,m;\; j = 1,\ldots,n\}.
    $$

    The corresponding eigenvectors are $\mathbf{u}_i \otimes \mathbf{v}_j$.

??? proof "Proof"
    Let $A\mathbf{u}_i = \lambda_i \mathbf{u}_i$ and $B\mathbf{v}_j = \mu_j \mathbf{v}_j$. Then:

    $$
    (A \oplus B)(\mathbf{u}_i \otimes \mathbf{v}_j) = (A \otimes I_n + I_m \otimes B)(\mathbf{u}_i \otimes \mathbf{v}_j)
    $$

    $$
    = (A\mathbf{u}_i) \otimes (I_n\mathbf{v}_j) + (I_m\mathbf{u}_i) \otimes (B\mathbf{v}_j) = \lambda_i(\mathbf{u}_i \otimes \mathbf{v}_j) + \mu_j(\mathbf{u}_i \otimes \mathbf{v}_j) = (\lambda_i + \mu_j)(\mathbf{u}_i \otimes \mathbf{v}_j).
    $$

    When $A, B$ are diagonalizable, this gives all $mn$ eigenvalues. The general case requires the Jordan normal form, but the conclusion is the same (counting algebraic multiplicities). $\blacksquare$

!!! theorem "Theorem 19.13 (Kronecker sum and matrix exponential)"
    Let $A$ be an $m \times m$ matrix and $B$ be an $n \times n$ matrix. Then:

    $$
    e^{A \oplus B} = e^A \otimes e^B.
    $$

??? proof "Proof"
    The key observation is that $A \otimes I$ and $I \otimes B$ **commute**:

    $$
    (A \otimes I)(I \otimes B) = A \otimes B = (I \otimes B)(A \otimes I).
    $$

    Since $A \oplus B = A \otimes I + I \otimes B$ and these two matrices commute, the matrix exponential satisfies:

    $$
    e^{A \oplus B} = e^{A \otimes I + I \otimes B} = e^{A \otimes I} \cdot e^{I \otimes B}.
    $$

    Furthermore, $e^{A \otimes I} = \sum_{k=0}^{\infty} \frac{(A \otimes I)^k}{k!} = \sum_{k=0}^{\infty} \frac{A^k \otimes I}{k!} = \left(\sum_{k=0}^{\infty}\frac{A^k}{k!}\right) \otimes I = e^A \otimes I$.

    Similarly $e^{I \otimes B} = I \otimes e^B$.

    Therefore $e^{A \oplus B} = (e^A \otimes I)(I \otimes e^B) = e^A \otimes e^B$. $\blacksquare$

!!! theorem "Theorem 19.14 (Kronecker sum and Lyapunov equation)"
    The Lyapunov equation $AX + XA^T = C$ is equivalent to:

    $$
    (A \oplus A^T)\operatorname{vec}(X) = (I \otimes A + A^* \otimes I)\operatorname{vec}(X) = \operatorname{vec}(C),
    $$

    where we note that $(A^T)^T = A$, so $I \otimes A + (A^T)^T \otimes I = I \otimes A + A \otimes I$.

    More precisely: after vectorization the equation becomes $(I_n \otimes A + \bar{A} \otimes I_n)\operatorname{vec}(X) = \operatorname{vec}(C)$ (in the real case $\bar{A} = A$).

    The eigenvalues of $A \oplus A^T$ are $\lambda_i(A) + \lambda_j(A^T) = \lambda_i(A) + \lambda_j(A)$, and the equation has a unique solution if and only if $\lambda_i(A) + \lambda_j(A) \neq 0$ for all $i,j$.

??? proof "Proof"
    Taking Vec of both sides of $AX + XA^T = C$:

    $$
    \operatorname{vec}(AX) + \operatorname{vec}(XA^T) = \operatorname{vec}(C).
    $$

    By Theorem 19.5: $\operatorname{vec}(AX) = (I \otimes A)\operatorname{vec}(X)$, $\operatorname{vec}(XA^T) = ((A^T)^T \otimes I)\operatorname{vec}(X) = (A \otimes I)\operatorname{vec}(X)$.

    Therefore $(I \otimes A + A \otimes I)\operatorname{vec}(X) = (A \oplus A)\operatorname{vec}(X) = \operatorname{vec}(C)$.

    Note that here $A \oplus A = A \otimes I + I \otimes A$ (in the real case, $A^T$ transposed back gives $A$). The eigenvalue conclusion follows directly from Theorem 19.12. $\blacksquare$

!!! example "Example 19.9"
    Let $A = \begin{pmatrix} -1 & 0 \\ 0 & -2 \end{pmatrix}$, $B = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$.

    $A \oplus B = A \otimes I_2 + I_2 \otimes B = \begin{pmatrix} -1 & 0 & 0 & 0 \\ 0 & -1 & 0 & 0 \\ 0 & 0 & -2 & 0 \\ 0 & 0 & 0 & -2 \end{pmatrix} + \begin{pmatrix} 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & -1 & 0 \end{pmatrix} = \begin{pmatrix} -1 & 1 & 0 & 0 \\ -1 & -1 & 0 & 0 \\ 0 & 0 & -2 & 1 \\ 0 & 0 & -1 & -2 \end{pmatrix}$.

    Eigenvalues of $A$: $-1, -2$. Eigenvalues of $B$: $i, -i$.

    Eigenvalues of $A \oplus B$: $\{-1+i, -1-i, -2+i, -2-i\}$.

    Verification: $\begin{pmatrix} -1 & 1 \\ -1 & -1 \end{pmatrix}$ has eigenvalues $-1 \pm i$, and $\begin{pmatrix} -2 & 1 \\ -1 & -2 \end{pmatrix}$ has eigenvalues $-2 \pm i$. Correct.

!!! example "Example 19.10"
    Verify that $e^{A \oplus B} = e^A \otimes e^B$.

    Take $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$, $B = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$.

    $e^A = \begin{pmatrix} e & 0 \\ 0 & 1 \end{pmatrix}$, $e^B = \begin{pmatrix} 1 & 0 \\ 0 & e \end{pmatrix}$.

    $e^A \otimes e^B = \begin{pmatrix} e & 0 & 0 & 0 \\ 0 & e^2 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & e \end{pmatrix}$.

    $A \oplus B = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix} + \begin{pmatrix} 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 2 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$.

    $e^{A \oplus B} = \begin{pmatrix} e & 0 & 0 & 0 \\ 0 & e^2 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & e \end{pmatrix} = e^A \otimes e^B$. Verification successful.
