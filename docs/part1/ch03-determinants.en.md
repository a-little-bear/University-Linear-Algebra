# Chapter 3  Determinants

<div class="context-flow" markdown>

**Prerequisites**: Chapter 2 square matrices · invertibility · elementary matrices

**Chapter arc**: Permutations and inversions → determinant definition → properties (effects of row operations) → cofactor expansion → adjugate matrix for inversion → determinant product formula → Cramer's rule → geometric interpretation

**Further connections**：Determinants appear in differential geometry (volume forms, Jacobians for change of variables), statistical mechanics (partition functions), and algebraic topology (incidence matrix determinants and Euler characteristic)

</div>

The determinant is a scalar value unique to square matrices, encoding a wealth of information about the matrix — invertibility, the volume-scaling effect of a linear transformation, the product of eigenvalues, etc. This chapter builds a rigorous definition of the determinant starting from permutations and inversions, systematically derives its fundamental properties, introduces the cofactor expansion theorem and various computation methods, and finally proves Cramer's rule and explains the geometric meaning of determinants.

---

## 3.1 Permutations and Inversions

<div class="context-flow" markdown>

**Combinatorial foundation of determinants**: The sign (even/odd) of each of the $n!$ permutations determines the sign of each term in the determinant definition

</div>

!!! definition "Definition 3.1 (Permutation)"
    An ordered arrangement $\sigma = (j_1, j_2, \ldots, j_n)$ of $1, 2, \ldots, n$ is called a **permutation** of $\{1, 2, \ldots, n\}$. There are $n!$ permutations of $n$ elements in total.

!!! definition "Definition 3.2 (Inversion and inversion number)"
    In a permutation $\sigma = (j_1, j_2, \ldots, j_n)$, if $s < t$ but $j_s > j_t$, then $(j_s, j_t)$ is called an **inversion**. The total number of inversions in $\sigma$ is called the **inversion number** of $\sigma$, denoted $\tau(\sigma)$ or $\operatorname{inv}(\sigma)$.

    - If $\tau(\sigma)$ is even, $\sigma$ is called an **even permutation**.
    - If $\tau(\sigma)$ is odd, $\sigma$ is called an **odd permutation**.

    The **sign** (signature) of $\sigma$ is defined as $\operatorname{sgn}(\sigma) = (-1)^{\tau(\sigma)}$.

!!! example "Example 3.1"
    The permutation $(3, 1, 2)$ has inversions $(3,1)$ and $(3,2)$, inversion number $2$, is an even permutation, with sign $+1$.

    The permutation $(3, 2, 1)$ has inversions $(3,2)$, $(3,1)$, $(2,1)$, inversion number $3$, is an odd permutation, with sign $-1$.

!!! proposition "Proposition 3.1"
    Swapping (transposing) two elements in a permutation changes the parity of the inversion number. Therefore, the parity of a permutation flips after one transposition.

??? proof "Proof"
    **Adjacent transposition**: Swapping adjacent $j_k$ and $j_{k+1}$ affects only this pair's inversion relation, changing the inversion number by $\pm 1$, so the parity changes.

    **General transposition**: Swapping positions $s$ and $t$ ($s < t$) can be decomposed into $2(t-s)-1$ adjacent transpositions (moving $j_t$ to position $s$ requires $t-s$ swaps, then moving the original $j_s$ from position $s+1$ to position $t$ requires $t-s-1$ swaps), totaling $2(t-s)-1$ (an odd number) adjacent transpositions, so the parity changes. $\blacksquare$

---

## 3.2 Definition of the Determinant

<div class="context-flow" markdown>

**From permutations to definition**: $\det(A) = \sum_\sigma \operatorname{sgn}(\sigma) \prod_i a_{i\sigma(i)}$ → each term takes exactly one entry from each row and each column → sum of $n!$ terms

</div>

!!! definition "Definition 3.3 (Determinant of order $n$)"
    Let $A = (a_{ij})$ be a square matrix of order $n$. The **determinant** of $A$ is defined as

    $$
    \det(A) = \sum_{\sigma \in S_n} \operatorname{sgn}(\sigma) \cdot a_{1\sigma(1)} a_{2\sigma(2)} \cdots a_{n\sigma(n)},
    $$

    where the sum is over all $n!$ permutations $\sigma$ of $\{1, 2, \ldots, n\}$, and $\operatorname{sgn}(\sigma) = (-1)^{\tau(\sigma)}$.

    The determinant is also written as $|A|$ or

    $$
    \begin{vmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{n1} & a_{n2} & \cdots & a_{nn} \end{vmatrix}.
    $$

!!! example "Example 3.2"
    $2 \times 2$ determinant:

    $$
    \begin{vmatrix} a & b \\ c & d \end{vmatrix} = ad - bc.
    $$

    $3 \times 3$ determinant (Sarrus' rule):

    $$
    \begin{vmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{vmatrix} = a_{11}a_{22}a_{33} + a_{12}a_{23}a_{31} + a_{13}a_{21}a_{32} - a_{13}a_{22}a_{31} - a_{12}a_{21}a_{33} - a_{11}a_{23}a_{32}.
    $$

!!! example "Example 3.3"
    Compute

    $$
    \begin{vmatrix} 2 & 1 & 3 \\ 0 & -1 & 2 \\ 1 & 0 & 1 \end{vmatrix} = 2 \cdot(-1)\cdot 1 + 1\cdot 2\cdot 1 + 3\cdot 0\cdot 0 - 3\cdot(-1)\cdot 1 - 1\cdot 0\cdot 1 - 2\cdot 2\cdot 0 = -2 + 2 + 0 + 3 - 0 - 0 = 3.
    $$

---

## 3.3 Properties of the Determinant

<div class="context-flow" markdown>

**Computational foundation**: Row interchange changes sign → row scaling multiplies by $c$ → row replacement leaves unchanged → these correspond to the three elementary row operations from Chapter 1 → $\det(AB)=\det(A)\det(B)$ is the essence of the product formula

</div>

The following properties form the theoretical basis for computing determinants.

!!! theorem "Theorem 3.1 (Invariance under transpose)"
    $\det(A^T) = \det(A)$.

??? proof "Proof"

    $$
    \det(A^T) = \sum_{\sigma} \operatorname{sgn}(\sigma) a_{\sigma(1)1}a_{\sigma(2)2}\cdots a_{\sigma(n)n}.
    $$

    Let $\tau = \sigma^{-1}$; then $\sigma(i) = j$ is equivalent to $\tau(j) = i$, and $\operatorname{sgn}(\tau) = \operatorname{sgn}(\sigma)$.

    $$
    \det(A^T) = \sum_{\tau} \operatorname{sgn}(\tau) a_{1\tau(1)}a_{2\tau(2)}\cdots a_{n\tau(n)} = \det(A). \quad \blacksquare
    $$

!!! note "Note"
    By invariance under transpose, any property of the determinant with respect to rows also holds for columns, and vice versa.

!!! theorem "Theorem 3.2 (Row interchange changes sign)"
    Swapping two rows (or two columns) of a matrix changes the sign of the determinant.

??? proof "Proof"
    Swapping rows $i$ and $j$ is equivalent to composing a transposition with each permutation $\sigma$ in the determinant definition. By Proposition 3.1, the sign of each permutation is negated, so the determinant changes sign. $\blacksquare$

!!! theorem "Theorem 3.3 (Homogeneity with respect to a row)"
    If all entries of a row (or column) of a matrix are multiplied by a constant $c$, the determinant is multiplied by $c$. That is,

    $$
    \det(\ldots, c\mathbf{r}_i, \ldots) = c \cdot \det(\ldots, \mathbf{r}_i, \ldots),
    $$

    where $\mathbf{r}_i$ denotes the $i$-th row.

??? proof "Proof"
    In the definition, each term contains exactly one entry from row $i$, namely $a_{i\sigma(i)}$. Multiplying this entry by $c$ multiplies each term by $c$, so the determinant is multiplied by $c$. $\blacksquare$

!!! corollary "Corollary 3.1"
    $\det(cA) = c^n \det(A)$, where $A$ is a square matrix of order $n$.

!!! theorem "Theorem 3.4 (Additivity with respect to a row)"
    If the $i$-th row of a matrix can be written as the sum of two row vectors $\mathbf{r}_i = \mathbf{r}_i' + \mathbf{r}_i''$, then

    $$
    \det(\ldots, \mathbf{r}_i' + \mathbf{r}_i'', \ldots) = \det(\ldots, \mathbf{r}_i', \ldots) + \det(\ldots, \mathbf{r}_i'', \ldots).
    $$

??? proof "Proof"
    In the definition, the entries of row $i$ are $a_{i\sigma(i)}' + a_{i\sigma(i)}''$. Expanding gives a sum of two parts, which correspond exactly to the two determinants on the right. $\blacksquare$

!!! theorem "Theorem 3.5 (Two identical rows implies zero determinant)"
    If a matrix has two identical rows (or columns), its determinant is $0$.

??? proof "Proof"
    Suppose rows $i$ and $j$ are identical. Swapping them leaves the matrix unchanged but negates the determinant (Theorem 3.2), so $\det(A) = -\det(A)$, giving $2\det(A) = 0$, i.e., $\det(A) = 0$. $\blacksquare$

!!! theorem "Theorem 3.6 (Row replacement does not change the determinant)"
    Adding $c$ times row $j$ to row $i$ does not change the determinant.

??? proof "Proof"

    $$
    \det(\ldots, \mathbf{r}_i + c\mathbf{r}_j, \ldots, \mathbf{r}_j, \ldots) = \det(\ldots, \mathbf{r}_i, \ldots, \mathbf{r}_j, \ldots) + c \cdot \det(\ldots, \mathbf{r}_j, \ldots, \mathbf{r}_j, \ldots).
    $$

    The second term has rows $i$ and $j$ both equal to $\mathbf{r}_j$, which is zero by Theorem 3.5. $\blacksquare$

!!! theorem "Theorem 3.7 (Determinant of a triangular matrix)"
    The determinant of an upper (or lower) triangular matrix equals the product of its diagonal entries:

    $$
    \det(A) = a_{11}a_{22}\cdots a_{nn}.
    $$

??? proof "Proof"
    Consider an upper triangular matrix. In the determinant formula $\sum_\sigma \operatorname{sgn}(\sigma) a_{1\sigma(1)} \cdots a_{n\sigma(n)}$, since $a_{ij} = 0$ ($i > j$), the product $a_{1\sigma(1)}\cdots a_{n\sigma(n)} \neq 0$ requires $\sigma(i) \ge i$ for each $i$. Since $\sigma$ is a bijection, the only possibility is $\sigma(i) = i$ (the identity permutation), whose sign is $+1$. Therefore $\det(A) = a_{11}a_{22}\cdots a_{nn}$. $\blacksquare$

<div class="context-flow" markdown>

**Key insight**: $\det(AB)=\det(A)\det(B)$ makes the determinant a **homomorphism** from the matrix multiplication group to the scalar multiplication group → Chapter 6 corollary: $\det A = \prod \lambda_i$

</div>

!!! theorem "Theorem 3.8 (Multiplicativity of determinants)"
    Let $A, B$ be square matrices of order $n$. Then

    $$
    \det(AB) = \det(A)\det(B).
    $$

??? proof "Proof"
    **Case 1**: $A$ is not invertible. Then $\operatorname{rank}(A) < n$, the columns of $AB$ are all linear combinations of the columns of $A$, $\operatorname{rank}(AB) \le \operatorname{rank}(A) < n$, so $AB$ is not invertible and $\det(AB) = 0$. Also $\det(A) = 0$, so $\det(AB) = 0 = \det(A)\det(B)$.

    **Case 2**: $A$ is invertible. Then $A$ can be written as a product of elementary matrices $A = E_1 E_2 \cdots E_k$. For elementary matrices, one can directly verify $\det(EA) = \det(E)\det(A)$ (by considering the three types of elementary matrices and using Theorems 3.2, 3.3, 3.6). Applying this repeatedly gives

    $$
    \det(AB) = \det(E_1 E_2 \cdots E_k B) = \det(E_1)\det(E_2)\cdots\det(E_k)\det(B) = \det(A)\det(B). \quad \blacksquare
    $$

---

## 3.4 Minors and Cofactors

<div class="context-flow" markdown>

**Order-reduction tool**: Delete row $i$ and column $j$ to get an $(n-1)$-order minor → with sign $(-1)^{i+j}$ gives the cofactor → prepares for Laplace expansion and the adjugate matrix

</div>

!!! definition "Definition 3.4 (Minor and cofactor)"
    Let $A$ be a square matrix of order $n$. The determinant of the $(n-1)$-order submatrix obtained by deleting row $i$ and column $j$ of $A$ is called the **minor** of $a_{ij}$, denoted $M_{ij}$.

    The **cofactor** of $a_{ij}$ is defined as

    $$
    A_{ij} = (-1)^{i+j} M_{ij}.
    $$

!!! example "Example 3.4"
    Let $A = \begin{pmatrix} 2 & 1 & 3 \\ 0 & -1 & 2 \\ 1 & 0 & 1 \end{pmatrix}$.

    The minor of $a_{11} = 2$ is $M_{11} = \begin{vmatrix} -1 & 2 \\ 0 & 1 \end{vmatrix} = -1$, and the cofactor is $A_{11} = (-1)^{1+1}(-1) = -1$.

    The minor of $a_{12} = 1$ is $M_{12} = \begin{vmatrix} 0 & 2 \\ 1 & 1 \end{vmatrix} = -2$, and the cofactor is $A_{12} = (-1)^{1+2}(-2) = 2$.

---

## 3.5 Cofactor Expansion Along a Row (or Column)

<div class="context-flow" markdown>

**Recursive computation**: $\det(A) = \sum_j a_{ij}A_{ij}$ (expand along row $i$) → expanding along a different row gives zero → combining gives $A \cdot \operatorname{adj}(A) = \det(A) \cdot I$ → from this follows $A^{-1} = \frac{1}{\det A}\operatorname{adj}(A)$

</div>

!!! theorem "Theorem 3.9 (Laplace expansion theorem)"
    The determinant $\det(A)$ of order $n$ can be expanded along row $i$ as

    $$
    \det(A) = \sum_{j=1}^n a_{ij} A_{ij} = a_{i1}A_{i1} + a_{i2}A_{i2} + \cdots + a_{in}A_{in},
    $$

    or expanded along column $j$ as

    $$
    \det(A) = \sum_{i=1}^n a_{ij} A_{ij} = a_{1j}A_{1j} + a_{2j}A_{2j} + \cdots + a_{nj}A_{nj}.
    $$

??? proof "Proof"
    Expand along row $i$. In the determinant formula, group all terms by the value of $a_{i\sigma(i)}$. For $\sigma(i) = j$, move row $i$ to the first row ($i-1$ row swaps) and column $j$ to the first column ($j-1$ column swaps); the determinant is multiplied by $(-1)^{(i-1)+(j-1)} = (-1)^{i+j}$. After removing the first row and first column, the remaining part is exactly the $(n-1)$-order determinant corresponding to $M_{ij}$. Therefore

    $$
    \det(A) = \sum_{j=1}^n a_{ij} (-1)^{i+j} M_{ij} = \sum_{j=1}^n a_{ij} A_{ij}. \quad \blacksquare
    $$

!!! theorem "Theorem 3.10 (Expansion with cofactors from a different row/column gives zero)"
    The sum of the entries of one row multiplied by the corresponding cofactors of a different row is zero:

    $$
    \sum_{k=1}^n a_{ik} A_{jk} = 0 \quad (i \neq j).
    $$

??? proof "Proof"
    Construct matrix $B$: replace row $j$ of $A$ with row $i$ (leaving other rows unchanged). Then $B$ has two identical rows (rows $i$ and $j$), so $\det(B) = 0$. Expanding $\det(B)$ along row $j$, note that the entries of row $j$ of $B$ are $a_{ik}$ ($k = 1, \ldots, n$), and the cofactors of row $j$ of $B$ are the same as those of row $j$ of $A$ (since after deleting row $j$, $B$ and $A$ agree on all other rows), giving

    $$
    0 = \det(B) = \sum_{k=1}^n a_{ik} A_{jk}. \quad \blacksquare
    $$

Combining Theorems 3.9 and 3.10, we can write

$$
\sum_{k=1}^n a_{ik}A_{jk} = \delta_{ij} \det(A),
$$

where $\delta_{ij}$ is the Kronecker delta. This is equivalent to the matrix identity $A \cdot \operatorname{adj}(A) = \det(A) \cdot I$.

!!! definition "Definition 3.5 (Adjugate matrix)"
    The **adjugate matrix** (classical adjoint) of a square matrix $A$ of order $n$ is defined as the transpose of the cofactor matrix:

    $$
    \operatorname{adj}(A) = (A_{ji})_{n \times n},
    $$

    i.e., the $(i, j)$-entry of $\operatorname{adj}(A)$ is $A_{ji}$.

!!! theorem "Theorem 3.11 (Adjugate matrix and inverse)"
    Let $A$ be a square matrix of order $n$. Then

    $$
    A \cdot \operatorname{adj}(A) = \operatorname{adj}(A) \cdot A = \det(A) \cdot I.
    $$

    Therefore, if $\det(A) \neq 0$, then $A$ is invertible and

    $$
    A^{-1} = \frac{1}{\det(A)} \operatorname{adj}(A).
    $$

!!! example "Example 3.5"
    Find the inverse using the adjugate matrix. Let $A = \begin{pmatrix} 1 & 2 \\ 3 & 5 \end{pmatrix}$.

    $\det(A) = 5 - 6 = -1$. Cofactors: $A_{11} = 5$, $A_{12} = -3$, $A_{21} = -2$, $A_{22} = 1$.

    $$
    \operatorname{adj}(A) = \begin{pmatrix} 5 & -2 \\ -3 & 1 \end{pmatrix}, \quad A^{-1} = \frac{1}{-1}\begin{pmatrix} 5 & -2 \\ -3 & 1 \end{pmatrix} = \begin{pmatrix} -5 & 2 \\ 3 & -1 \end{pmatrix}.
    $$

---

## 3.6 Methods for Computing Determinants

<div class="context-flow" markdown>

**Practical methods**: Triangularization (using properties 3.2–3.6) · recurrence method (using expansion to build recurrences) · Vandermonde determinant (classical formula) · block determinants

</div>

### Triangularization Method

Use row operations to reduce the matrix to upper triangular form (noting that row interchanges change sign and row scaling multiplies by the scaling factor), then use Theorem 3.7 to compute the determinant directly.

!!! example "Example 3.6"
    Compute

    $$
    D = \begin{vmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 0 \end{vmatrix}
    $$

    $$
    \xrightarrow{R_2 - 4R_1, \; R_3 - 7R_1}
    \begin{vmatrix} 1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & -6 & -21 \end{vmatrix}
    \xrightarrow{R_3 - 2R_2}
    \begin{vmatrix} 1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & 0 & -9 \end{vmatrix} = 1 \cdot (-3) \cdot (-9) = 27.
    $$

### Recurrence Method

For determinants with special patterns, recurrence relations can be established to find the solution.

!!! example "Example 3.7"
    Let $D_n$ be the determinant of order $n$:

    $$
    D_n = \begin{vmatrix} 2 & 1 & 0 & \cdots & 0 & 0 \\ 1 & 2 & 1 & \cdots & 0 & 0 \\ 0 & 1 & 2 & \cdots & 0 & 0 \\ \vdots & & & \ddots & & \vdots \\ 0 & 0 & 0 & \cdots & 2 & 1 \\ 0 & 0 & 0 & \cdots & 1 & 2 \end{vmatrix}
    $$

    Expanding along the first row: $D_n = 2D_{n-1} - D_{n-2}$. Initial values: $D_1 = 2$, $D_2 = 3$.

    The characteristic equation $x^2 - 2x + 1 = 0$ has a repeated root $x = 1$, so the general solution is $D_n = (c_1 + c_2 n) \cdot 1^n = c_1 + c_2 n$.

    From $D_1 = 2, D_2 = 3$, we get $c_1 = 1, c_2 = 1$, so $D_n = n + 1$.

### Vandermonde Determinant

!!! theorem "Theorem 3.12 (Vandermonde determinant)"
    The **Vandermonde determinant** is

    $$
    V_n = \begin{vmatrix} 1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\ 1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 1 & x_n & x_n^2 & \cdots & x_n^{n-1} \end{vmatrix} = \prod_{1 \le i < j \le n} (x_j - x_i).
    $$

??? proof "Proof"
    By induction on $n$. For $n = 2$, $V_2 = x_2 - x_1$, which holds.

    Assume the formula holds for order $n-1$. For the order-$n$ Vandermonde determinant, starting from the last row, subtract $x_1$ times the row above from each row ($R_n - x_1 R_{n-1}$, $R_{n-1} - x_1 R_{n-2}$, ..., $R_2 - x_1 R_1$):

    $$
    V_n = \begin{vmatrix} 1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\ 0 & x_2 - x_1 & x_2(x_2 - x_1) & \cdots & x_2^{n-2}(x_2 - x_1) \\ 0 & x_3 - x_1 & x_3(x_3 - x_1) & \cdots & x_3^{n-2}(x_3 - x_1) \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & x_n - x_1 & x_n(x_n - x_1) & \cdots & x_n^{n-2}(x_n - x_1) \end{vmatrix}
    $$

    Expanding along the first column and factoring out $(x_i - x_1)$ from each row $i$ ($i = 2, \ldots, n$):

    $$
    V_n = \prod_{i=2}^n (x_i - x_1) \cdot V_{n-1}(x_2, x_3, \ldots, x_n)
    $$

    By the induction hypothesis, $V_{n-1}(x_2, \ldots, x_n) = \prod_{2 \le i < j \le n}(x_j - x_i)$, so

    $$
    V_n = \prod_{i=2}^n(x_i - x_1) \cdot \prod_{2 \le i < j \le n}(x_j - x_i) = \prod_{1 \le i < j \le n}(x_j - x_i). \quad \blacksquare
    $$

### Block Determinants

!!! theorem "Theorem 3.13 (Determinant of block triangular matrices)"
    If

    $$
    M = \begin{pmatrix} A & B \\ O & D \end{pmatrix} \quad \text{or} \quad M = \begin{pmatrix} A & O \\ C & D \end{pmatrix},
    $$

    where $A$ and $D$ are square matrices, then $\det(M) = \det(A)\det(D)$.

!!! example "Example 3.8"
    Compute the determinant of the block matrix:

    $$
    M = \begin{pmatrix} 1 & 2 & 0 & 0 \\ 3 & 4 & 0 & 0 \\ 5 & 6 & 1 & 1 \\ 7 & 8 & 0 & 2 \end{pmatrix}
    $$

    This is a block lower triangular matrix $\begin{pmatrix} A & O \\ C & D \end{pmatrix}$, where $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$, $D = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}$.

    $\det(M) = \det(A)\det(D) = (4 - 6)(2 - 0) = (-2)(2) = -4$.

---

## 3.7 Cramer's Rule

<div class="context-flow" markdown>

**Solving equations via determinants**: $x_j = \det(A_j)/\det(A)$ → theoretical significance > computational significance ($n+1$ determinants is too expensive) → in practice, Gaussian elimination from Chapter 1 is still used

</div>

!!! theorem "Theorem 3.14 (Cramer's rule)"
    Let $A$ be an invertible matrix of order $n$ (i.e., $\det(A) \neq 0$). Then the system $A\mathbf{x} = \mathbf{b}$ has a unique solution, and

    $$
    x_j = \frac{\det(A_j)}{\det(A)}, \quad j = 1, 2, \ldots, n,
    $$

    where $A_j$ is the matrix obtained by replacing the $j$-th column of $A$ with $\mathbf{b}$.

??? proof "Proof"
    Since $A$ is invertible, the system has a unique solution $\mathbf{x} = A^{-1}\mathbf{b}$. Using the adjugate matrix formula:

    $$
    x_j = (A^{-1}\mathbf{b})_j = \frac{1}{\det(A)}\sum_{i=1}^n A_{ji} b_i.
    $$

    The sum $\sum_{i=1}^n A_{ji} b_i$ is exactly the cofactor expansion along the $j$-th column of the matrix $A_j$ obtained by replacing the $j$-th column of $A$ with $\mathbf{b}$, i.e., $\det(A_j)$. Therefore $x_j = \frac{\det(A_j)}{\det(A)}$. $\blacksquare$

!!! note "Note"
    Cramer's rule has greater theoretical than computational significance. For large systems, Gaussian elimination is far more efficient than Cramer's rule (the latter requires computing $n+1$ determinants of order $n$).

!!! example "Example 3.9"
    Solve using Cramer's rule:

    $$
    \begin{cases} 2x_1 + x_2 = 5 \\ 3x_1 + 2x_2 = 8 \end{cases}
    $$

    $$
    \det(A) = \begin{vmatrix} 2 & 1 \\ 3 & 2 \end{vmatrix} = 1, \quad \det(A_1) = \begin{vmatrix} 5 & 1 \\ 8 & 2 \end{vmatrix} = 2, \quad \det(A_2) = \begin{vmatrix} 2 & 5 \\ 3 & 8 \end{vmatrix} = 1.
    $$

    $$
    x_1 = \frac{2}{1} = 2, \quad x_2 = \frac{1}{1} = 1.
    $$

---

## 3.8 Geometric Interpretation of Determinants

<div class="context-flow" markdown>

**Bridge between algebra and geometry**: $|\det A|$ = volume of the parallelepiped spanned by the column vectors → the sign of $\det A$ = orientation (right-hand / left-hand system) → a linear transformation scales area by $|\det A|$ (→ Chapter 5)

</div>

Determinants have a deep connection to area and volume in geometry.

!!! theorem "Theorem 3.15 (Determinants and area/volume)"
    Let $\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n \in \mathbb{R}^n$ be $n$ column vectors forming the matrix $A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)$.

    1. **Two dimensions**: The **signed area** of the parallelogram with adjacent sides $\mathbf{a}_1, \mathbf{a}_2 \in \mathbb{R}^2$ equals $\det(A) = \begin{vmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{vmatrix}$.

    2. **Three dimensions**: The **signed volume** of the parallelepiped with adjacent edges $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3 \in \mathbb{R}^3$ equals $\det(A)$.

    3. **General case**: The absolute value of the $n$-dimensional volume of the $n$-dimensional parallelepiped with edges $\mathbf{a}_1, \ldots, \mathbf{a}_n$ is $|\det(A)|$.

!!! note "Note"
    $\det(A) > 0$ means the set $\{\mathbf{a}_1, \ldots, \mathbf{a}_n\}$ preserves the orientation of the standard basis (right-hand system); $\det(A) < 0$ means it reverses the orientation. $\det(A) = 0$ means the vectors are linearly dependent, the "parallelepiped" degenerates to a lower-dimensional object, and the volume is zero.

!!! example "Example 3.10"
    The area of the parallelogram spanned by vectors $\mathbf{a}_1 = (3, 0)^T$ and $\mathbf{a}_2 = (1, 2)^T$ is

    $$
    |\det(A)| = \left|\begin{vmatrix} 3 & 1 \\ 0 & 2 \end{vmatrix}\right| = |6 - 0| = 6.
    $$

    Geometrically, this parallelogram has base $3$ and height $2$, with area indeed $6$.

!!! example "Example 3.11"
    A linear transformation $T: \mathbb{R}^2 \to \mathbb{R}^2$, $T(\mathbf{x}) = A\mathbf{x}$, scales the area of any region in the plane by a factor of $|\det(A)|$. For example, if $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$, then $\det(A) = 6$, meaning $T$ enlarges the area of any region by a factor of $6$.
