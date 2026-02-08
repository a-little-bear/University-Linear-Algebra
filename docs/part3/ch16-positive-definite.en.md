# Chapter 16  Positive Definite Matrices

<div class="context-flow" markdown>

**Prerequisites**: Spectral decomposition of Hermitian matrices (Ch6)

**Chapter arc**: $A\succ 0$ $\Leftrightarrow$ all eigenvalues positive $\Leftrightarrow$ Cholesky decomposable $\Leftrightarrow$ positive quadratic form → **Schur complement** is the key to block positive definiteness → Ch18 inequalities

**Further connections**：Positive definite matrices pervade statistics (covariance matrices), machine learning (kernel matrices, Fisher information), optimization (Newton's method Hessian), and elasticity (stress-strain relations); Löwner partial order has deep applications in matrix monotone functions and quantum information theory

</div>

Positive definite matrices are among the most important classes of matrices in linear algebra. They appear ubiquitously in optimization, statistics, differential equations, signal processing, and machine learning. The theory of positive definite matrices connects quadratic forms, eigenvalues, determinants, and matrix decompositions, forming a cornerstone for applications of linear algebra. This chapter systematically develops the theory of positive definite and positive semidefinite matrices, discusses the Schur complement, Lowner partial order, operational properties of positive definite matrices, and several important determinantal inequalities.

---

## 16.1 Definition and equivalent conditions

<div class="context-flow" markdown>

**Intuition**: Core logic of 8 equivalent conditions -- positive quadratic form $\leftrightarrow$ all eigenvalues positive $\leftrightarrow$ $A=C^*C$ ($C$ invertible) $\leftrightarrow$ Cholesky $\leftrightarrow$ all leading principal minors positive

</div>

!!! definition "Definition 16.1 (Positive definite matrix)"
    Let $A \in \mathbb{C}^{n \times n}$ be a Hermitian matrix ($A = A^*$). If for all nonzero vectors $\mathbf{x} \in \mathbb{C}^n$,

    $$\mathbf{x}^*A\mathbf{x} > 0,$$

    then $A$ is called a **positive definite matrix**, denoted $A \succ 0$.

!!! theorem "Theorem 16.1 (Equivalent conditions for positive definiteness)"
    Let $A \in \mathbb{C}^{n \times n}$ be a Hermitian matrix. The following conditions are equivalent:

    (1) $A$ is positive definite, i.e., $\mathbf{x}^*A\mathbf{x} > 0$ for all $\mathbf{x} \neq \mathbf{0}$;

    (2) All eigenvalues of $A$ are positive: $\lambda_1, \lambda_2, \ldots, \lambda_n > 0$;

    (3) All leading principal minors of $A$ are positive: $\det(A_k) > 0$, $k = 1, 2, \ldots, n$ (where $A_k$ is the upper-left $k \times k$ submatrix of $A$);

    (4) There exists an invertible matrix $C$ such that $A = C^*C$;

    (5) There exists an invertible lower triangular matrix $L$ (with positive diagonal entries) such that $A = LL^*$ (Cholesky decomposition);

    (6) There exists an invertible upper triangular matrix $R$ such that $A = R^*R$;

    (7) All principal minors (not just leading ones) of $A$ are positive;

    (8) There exists a positive definite matrix $B$ such that $A = B^2$ (positive definite square root).

??? proof "Proof"
    **(1) $\Rightarrow$ (2)**: Let $A\mathbf{v} = \lambda\mathbf{v}$ ($\mathbf{v} \neq \mathbf{0}$); then $\lambda = \frac{\mathbf{v}^*A\mathbf{v}}{\mathbf{v}^*\mathbf{v}} > 0$.

    **(2) $\Rightarrow$ (1)**: $A$ is Hermitian; let $A = U\Lambda U^*$ be the spectral decomposition, where $\Lambda = \operatorname{diag}(\lambda_1,\ldots,\lambda_n)$, $\lambda_i > 0$. For $\mathbf{x} \neq \mathbf{0}$, set $\mathbf{y} = U^*\mathbf{x} \neq \mathbf{0}$; then $\mathbf{x}^*A\mathbf{x} = \mathbf{y}^*\Lambda\mathbf{y} = \sum_i\lambda_i|y_i|^2 > 0$.

    **(1) $\Rightarrow$ (3)**: $A$ positive definite $\Rightarrow$ every leading principal submatrix $A_k$ is also positive definite (restriction to a subspace) $\Rightarrow$ all eigenvalues of $A_k$ are positive $\Rightarrow$ $\det(A_k) > 0$.

    **(3) $\Rightarrow$ (5)**: By induction and Gaussian elimination. For $n = 1$: $A = (a_{11})$, $a_{11} > 0$, take $L = (\sqrt{a_{11}})$. Assume the result holds for order $n-1$. Partition $A$ as:

    $$A = \begin{pmatrix} A_{n-1} & \mathbf{a} \\ \mathbf{a}^* & a_{nn} \end{pmatrix}.$$

    By the induction hypothesis, $A_{n-1} = L_{n-1}L_{n-1}^*$. Set $\mathbf{l} = L_{n-1}^{-1}\mathbf{a}$, $\ell_{nn} = \sqrt{a_{nn} - \mathbf{l}^*\mathbf{l}}$ (one can verify that the expression under the square root is positive using $\det A > 0$ and $\det A_{n-1} > 0$). Then

    $$L = \begin{pmatrix} L_{n-1} & \mathbf{0} \\ \mathbf{l}^* & \ell_{nn} \end{pmatrix}, \quad A = LL^*.$$

    **(5) $\Rightarrow$ (4)**: Take $C = L^*$.

    **(4) $\Rightarrow$ (1)**: $\mathbf{x}^*A\mathbf{x} = \mathbf{x}^*C^*C\mathbf{x} = \|C\mathbf{x}\|^2 > 0$ (since $C$ is invertible, $C\mathbf{x} \neq \mathbf{0}$).

    **(2) $\Rightarrow$ (8)**: Let $A = U\Lambda U^*$; take $B = U\Lambda^{1/2}U^*$ (where $\Lambda^{1/2} = \operatorname{diag}(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_n})$); then $B$ is positive definite and $B^2 = A$.

    **(8) $\Rightarrow$ (1)**: $\mathbf{x}^*A\mathbf{x} = \mathbf{x}^*B^2\mathbf{x} = \|B\mathbf{x}\|^2 > 0$ ($B$ is invertible).

    **(2) $\Rightarrow$ (7)**: Any principal submatrix of a positive definite matrix (corresponding to selecting certain rows and the same column indices) is also positive definite, so its determinant is positive. $\blacksquare$

!!! definition "Definition 16.2 (Cholesky decomposition)"
    Let $A$ be a positive definite Hermitian matrix. The **Cholesky decomposition** of $A$ is the unique lower triangular matrix $L$ (with positive diagonal entries) such that $A = LL^*$.

!!! example "Example 16.1"
    Perform the Cholesky decomposition of $A = \begin{pmatrix} 4 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 6 \end{pmatrix}$.

    $l_{11} = \sqrt{4} = 2$.

    $l_{21} = a_{21}/l_{11} = 2/2 = 1$, $l_{31} = a_{31}/l_{11} = 1/2$.

    $l_{22} = \sqrt{a_{22} - l_{21}^2} = \sqrt{5 - 1} = 2$.

    $l_{32} = (a_{32} - l_{31}l_{21})/l_{22} = (3 - 1/2)/2 = 5/4$.

    $l_{33} = \sqrt{a_{33} - l_{31}^2 - l_{32}^2} = \sqrt{6 - 1/4 - 25/16} = \sqrt{6 - 29/16} = \sqrt{67/16} = \frac{\sqrt{67}}{4}$.

    $$L = \begin{pmatrix} 2 & 0 & 0 \\ 1 & 2 & 0 \\ 1/2 & 5/4 & \sqrt{67}/4 \end{pmatrix}.$$

    One can verify that $LL^T = A$.

!!! example "Example 16.2"
    Verifying the equivalent conditions for positive definiteness. Let $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$.

    - **Condition (2)**: Eigenvalues $\lambda = 2 \pm 1 = 3, 1$, all positive.
    - **Condition (3)**: $\det(A_1) = 2 > 0$, $\det(A) = 4 - 1 = 3 > 0$.
    - **Condition (4)**: $A = C^TC$, where $C = \begin{pmatrix} \sqrt{2} & -1/\sqrt{2} \\ 0 & \sqrt{3/2} \end{pmatrix}$.
    - **Condition (1)**: $\mathbf{x}^TA\mathbf{x} = 2x_1^2 - 2x_1x_2 + 2x_2^2 = (x_1 - x_2)^2 + x_1^2 + x_2^2 > 0$ (for $\mathbf{x} \neq \mathbf{0}$).

---

## 16.2 Positive semidefinite matrices

<div class="context-flow" markdown>

**Chapter arc**: $A\succeq 0$ $\Leftrightarrow$ $\lambda_i\geq 0$ $\Leftrightarrow$ $A=C^*C$ · The positive semidefinite cone $\mathbb{S}_+^n$ is a convex cone → Ch18 Lowner order and optimization

</div>

!!! definition "Definition 16.3 (Positive semidefinite matrix)"
    Let $A \in \mathbb{C}^{n \times n}$ be a Hermitian matrix. If for all $\mathbf{x} \in \mathbb{C}^n$,

    $$\mathbf{x}^*A\mathbf{x} \geq 0,$$

    then $A$ is called a **positive semidefinite matrix**, denoted $A \succeq 0$.

!!! theorem "Theorem 16.2 (Equivalent conditions for positive semidefiniteness)"
    Let $A$ be an $n \times n$ Hermitian matrix. The following conditions are equivalent:

    (1) $A$ is positive semidefinite;

    (2) All eigenvalues of $A$ are non-negative: $\lambda_i \geq 0$;

    (3) There exists a matrix $C$ (not necessarily square or invertible) such that $A = C^*C$;

    (4) All principal minors of $A$ are non-negative;

    (5) There exists a positive semidefinite matrix $B$ such that $A = B^2$.

??? proof "Proof"
    The proof is similar to the positive definite case, replacing strict inequalities with non-strict ones.

    **(1) $\Rightarrow$ (2)**: If $A\mathbf{v} = \lambda\mathbf{v}$, $\mathbf{v} \neq \mathbf{0}$, then $\lambda = \frac{\mathbf{v}^*A\mathbf{v}}{\|\mathbf{v}\|^2} \geq 0$.

    **(2) $\Rightarrow$ (3)**: Let $A = U\Lambda U^*$; take $C = \Lambda^{1/2}U^*$; then $C^*C = U\Lambda U^* = A$.

    **(3) $\Rightarrow$ (1)**: $\mathbf{x}^*A\mathbf{x} = \|C\mathbf{x}\|^2 \geq 0$. $\blacksquare$

!!! definition "Definition 16.4 (Positive semidefinite cone)"
    The set of all $n \times n$ positive semidefinite matrices is denoted $\mathbb{S}_+^n$ and called the **positive semidefinite cone**. It is a convex cone in $\mathbb{S}^n$ (the space of $n \times n$ real symmetric matrices):

    - If $A, B \in \mathbb{S}_+^n$, then $A + B \in \mathbb{S}_+^n$;
    - If $A \in \mathbb{S}_+^n$ and $\alpha \geq 0$, then $\alpha A \in \mathbb{S}_+^n$.

!!! example "Example 16.3"
    The matrix $A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$ is positive semidefinite.

    Eigenvalues: $\lambda^2 - 5\lambda = 0$, $\lambda_1 = 5$, $\lambda_2 = 0$. All non-negative.

    $\mathbf{x}^TA\mathbf{x} = x_1^2 + 4x_1x_2 + 4x_2^2 = (x_1 + 2x_2)^2 \geq 0$.

    $A$ is not positive definite because taking $\mathbf{x} = (-2, 1)^T$ gives $\mathbf{x}^TA\mathbf{x} = 0$.

    $A = C^TC$, where $C = (1, 2)$ (a $1 \times 2$ matrix).

---

## 16.3 Schur complement

<div class="context-flow" markdown>

**Core tool**: $M/A = C-B^*A^{-1}B$ · Block LDL decomposition → $M\succ 0 \Leftrightarrow A\succ 0$ and $M/A\succ 0$ · Determinant formula $\det M=\det A\cdot\det(M/A)$ → Fischer inequality

</div>

The Schur complement is the core tool for analyzing positive definiteness of block matrices, with wide applications in control theory, statistics, and optimization.

!!! definition "Definition 16.5 (Schur complement)"
    Let the block matrix $M = \begin{pmatrix} A & B \\ B^* & C \end{pmatrix}$, where $A$ is invertible. The **Schur complement** of $A$ in $M$ is defined as

    $$M/A = C - B^*A^{-1}B.$$

    Similarly, if $C$ is invertible, the Schur complement of $C$ in $M$ is $M/C = A - BC^{-1}B^*$.

!!! theorem "Theorem 16.3 (Positive definiteness criterion for block matrices)"
    Let $M = \begin{pmatrix} A & B \\ B^* & C \end{pmatrix}$ be a Hermitian matrix.

    (1) If $A \succ 0$, then $M \succ 0 \Leftrightarrow M/A = C - B^*A^{-1}B \succ 0$;

    (2) If $C \succ 0$, then $M \succ 0 \Leftrightarrow M/C = A - BC^{-1}B^* \succ 0$;

    (3) $M \succeq 0$ and $A \succ 0$ $\Leftrightarrow$ $A \succ 0$ and $M/A \succeq 0$.

??? proof "Proof"
    **(1)** Perform a block LDL decomposition. Since $A \succ 0$, we can write

    $$M = \begin{pmatrix} I & O \\ B^*A^{-1} & I \end{pmatrix}\begin{pmatrix} A & O \\ O & C-B^*A^{-1}B \end{pmatrix}\begin{pmatrix} I & A^{-1}B \\ O & I \end{pmatrix}.$$

    One can directly verify that the right-hand product equals $M$.

    Let $L = \begin{pmatrix} I & O \\ B^*A^{-1} & I \end{pmatrix}$, $D = \begin{pmatrix} A & O \\ O & M/A \end{pmatrix}$. Then $M = LDL^*$.

    Since $L$ is invertible, $M$ is positive definite $\Leftrightarrow$ $D$ is positive definite $\Leftrightarrow$ $A \succ 0$ and $M/A \succ 0$. Given $A \succ 0$, we have $M \succ 0 \Leftrightarrow M/A \succ 0$.

    **(2)** Similarly, by performing block elimination on $C$. $\blacksquare$

!!! theorem "Theorem 16.4 (Determinant formula via Schur complement)"
    If $A$ is invertible, then

    $$\det\begin{pmatrix} A & B \\ C & D \end{pmatrix} = \det(A)\cdot\det(D - CA^{-1}B).$$

??? proof "Proof"
    $$\begin{pmatrix} A & B \\ C & D \end{pmatrix} = \begin{pmatrix} I & O \\ CA^{-1} & I \end{pmatrix}\begin{pmatrix} A & B \\ O & D - CA^{-1}B \end{pmatrix}.$$

    Taking determinants: $\det(M) = 1 \cdot \det(A)\det(D - CA^{-1}B)$. $\blacksquare$

!!! example "Example 16.4"
    Let $M = \begin{pmatrix} 4 & 2 & 1 \\ 2 & 3 & 1 \\ 1 & 1 & 2 \end{pmatrix}$. Determine whether $M$ is positive definite.

    Take $A = \begin{pmatrix} 4 & 2 \\ 2 & 3 \end{pmatrix}$ (positive definite since $4 > 0$, $\det = 8 > 0$), $B = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$, $C = (2)$.

    $M/A = C - B^TA^{-1}B = 2 - (1, 1)\frac{1}{8}\begin{pmatrix} 3 & -2 \\ -2 & 4 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = 2 - \frac{1}{8}(1, 1)\begin{pmatrix} 1 \\ 2 \end{pmatrix} = 2 - \frac{3}{8} = \frac{13}{8} > 0$.

    Therefore $M \succ 0$.

!!! example "Example 16.5"
    **Schur complement in linear regression.** In the least-squares problem, the coefficient matrix $X^TX$ of the normal equations is positive semidefinite. When the model is partitioned into blocks (e.g., intercept term and other variables), the Schur complement arises naturally.

    Let $X = (\mathbf{1}, X_2)$; then

    $$(X^TX)/ ((\mathbf{1}^T\mathbf{1})) = X_2^TX_2 - X_2^T\mathbf{1}(\mathbf{1}^T\mathbf{1})^{-1}\mathbf{1}^TX_2 = X_2^T(I - \frac{1}{n}\mathbf{1}\mathbf{1}^T)X_2,$$

    which is the centered $X_2^TX_2$, i.e., the sample covariance matrix (without the factor $\frac{1}{n-1}$).

---

## 16.4 Lowner partial order

<div class="context-flow" markdown>

**Chapter arc**: $A\succeq B \Leftrightarrow A-B\succeq 0$ defines a partial order on matrices · Congruence invariance + **order reversal** ($A\succeq B\succ 0 \Rightarrow B^{-1}\succeq A^{-1}$) → Ch18 matrix monotone functions

</div>

!!! definition "Definition 16.6 (Lowner partial order)"
    For Hermitian matrices $A$ and $B$, the **Lowner partial order** is defined as:

    $$A \succeq B \quad \Longleftrightarrow \quad A - B \succeq 0 \quad (\text{i.e., } A - B \text{ is positive semidefinite}),$$

    $$A \succ B \quad \Longleftrightarrow \quad A - B \succ 0 \quad (\text{i.e., } A - B \text{ is positive definite}).$$

!!! theorem "Theorem 16.5 (Properties of the Lowner partial order)"
    Let $A, B, C$ be Hermitian matrices of the same size. The Lowner partial order satisfies:

    (1) **Reflexivity**: $A \succeq A$;

    (2) **Antisymmetry**: If $A \succeq B$ and $B \succeq A$, then $A = B$;

    (3) **Transitivity**: If $A \succeq B$ and $B \succeq C$, then $A \succeq C$;

    (4) **Additivity**: If $A \succeq B$, then $A + C \succeq B + C$;

    (5) **Congruence invariance**: If $A \succeq B$, then for any matrix $P$, $P^*AP \succeq P^*BP$;

    (6) **Order reversal**: If $A \succeq B \succ 0$, then $B^{-1} \succeq A^{-1}$.

??? proof "Proof"
    (1)--(4) follow directly from properties of the positive semidefinite cone.

    **(5)** $P^*AP - P^*BP = P^*(A-B)P$. Since $A - B \succeq 0$, for any $\mathbf{x}$, $\mathbf{x}^*P^*(A-B)P\mathbf{x} = (P\mathbf{x})^*(A-B)(P\mathbf{x}) \geq 0$.

    **(6)** From $A \succeq B \succ 0$, for any $\mathbf{x} \neq \mathbf{0}$:

    More directly: $A \succeq B \Rightarrow B^{-1/2}AB^{-1/2} \succeq I \Rightarrow (B^{-1/2}AB^{-1/2})^{-1} \preceq I$ (since $f(t) = 1/t$ is operator monotone decreasing on the positive half-line) $\Rightarrow B^{1/2}A^{-1}B^{1/2} \preceq I \Rightarrow A^{-1} \preceq B^{-1}$. $\blacksquare$

!!! example "Example 16.6"
    Let $A = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$, $B = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$.

    $A - B = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$, eigenvalues $\frac{3 \pm \sqrt{5}}{2} > 0$, so $A \succ B$.

    By order reversal, $B^{-1} \succeq A^{-1}$. Verification: $B^{-1} = \begin{pmatrix} 1/2 & 0 \\ 0 & 1 \end{pmatrix}$, $A^{-1} = \frac{1}{8}\begin{pmatrix} 3 & -1 \\ -1 & 3 \end{pmatrix}$.

    $B^{-1} - A^{-1} = \begin{pmatrix} 1/2 - 3/8 & 1/8 \\ 1/8 & 1 - 3/8 \end{pmatrix} = \begin{pmatrix} 1/8 & 1/8 \\ 1/8 & 5/8 \end{pmatrix}$, determinant $= 5/64 - 1/64 = 4/64 > 0$, and $1/8 > 0$, so it is positive semidefinite (in fact positive definite).

---

## 16.5 Operational properties of positive definite matrices

<div class="context-flow" markdown>

**Chapter arc**: Sum/scalar multiplication preserve positive definiteness · Eigenvalues of $AB$ are all positive (but $AB$ is not necessarily Hermitian)

**Schur product theorem** ($A\circ B\succeq 0$)

**Kronecker product** preserves positive definiteness → Ch19

</div>

!!! theorem "Theorem 16.6 (Sum and scalar multiplication)"
    (1) If $A \succ 0$ and $B \succ 0$, then $A + B \succ 0$;

    (2) If $A \succ 0$ and $\alpha > 0$, then $\alpha A \succ 0$.

!!! theorem "Theorem 16.7 (Positive definiteness of products)"
    Let $A, B$ be positive definite Hermitian matrices. Then $AB$ has all positive real eigenvalues, but $AB$ is not necessarily Hermitian (unless $AB = BA$).

??? proof "Proof"
    $AB$ is similar to $A^{1/2}BA^{1/2}$ (since $AB = A^{1/2}(A^{1/2}B)$ and $A^{1/2}BA^{1/2}$ have the same eigenvalues). $A^{1/2}BA^{1/2}$ is a positive definite Hermitian matrix (by congruence invariance), so its eigenvalues are all positive. Hence $AB$ has all positive eigenvalues. $\blacksquare$

!!! definition "Definition 16.7 (Hadamard product)"
    Let $A = (a_{ij})$, $B = (b_{ij}) \in \mathbb{C}^{m \times n}$. Their **Hadamard product** (also called **Schur product**) is defined as

    $$A \circ B = (a_{ij}b_{ij}),$$

    i.e., entrywise multiplication.

!!! theorem "Theorem 16.8 (Schur product theorem)"
    If $A \succeq 0$ and $B \succeq 0$, then $A \circ B \succeq 0$. If furthermore $A \succ 0$ and $B \succ 0$, then $A \circ B \succ 0$.

??? proof "Proof"
    Since $B \succeq 0$, let $B = \sum_{k=1}^{r} \mathbf{b}_k\mathbf{b}_k^*$ (rank-one decomposition). Then

    $$A \circ B = A \circ \left(\sum_k \mathbf{b}_k\mathbf{b}_k^*\right) = \sum_k A \circ (\mathbf{b}_k\mathbf{b}_k^*) = \sum_k D_k A D_k^*,$$

    where $D_k = \operatorname{diag}(\mathbf{b}_k)$ (diagonal matrix with components of $\mathbf{b}_k$ on the diagonal).

    This is because $(A \circ \mathbf{b}\mathbf{b}^*)_{ij} = a_{ij}b_ib_j^* = (DAD^*)_{ij}$, where $D = \operatorname{diag}(b_1, \ldots, b_n)$.

    Each $D_kAD_k^*$ is positive semidefinite (by congruence invariance), and their sum is also positive semidefinite.

    If $A \succ 0$ and $B \succ 0$, since the diagonal entries of $B$ are all positive ($b_{ii} > 0$), at least one $D_k$ is invertible, making $D_kAD_k^*$ positive definite. $\blacksquare$

!!! theorem "Theorem 16.9 (Kronecker product preserves positive definiteness)"
    If $A \succ 0$ and $B \succ 0$, then $A \otimes B \succ 0$.

??? proof "Proof"
    Let the eigenvalues of $A$ be $\alpha_1, \ldots, \alpha_m > 0$ and those of $B$ be $\beta_1, \ldots, \beta_n > 0$. The eigenvalues of $A \otimes B$ are $\alpha_i\beta_j > 0$. Alternatively: $A \otimes B$ is Hermitian, and $(\mathbf{x} \otimes \mathbf{y})^*(A \otimes B)(\mathbf{x} \otimes \mathbf{y}) = (\mathbf{x}^*A\mathbf{x})(\mathbf{y}^*B\mathbf{y}) > 0$. General vectors are linear combinations of such tensor products, and the argument can be completed similarly. $\blacksquare$

!!! example "Example 16.7"
    Let $A = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix} \succ 0$, $B = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 2 \end{pmatrix} \succ 0$.

    **Hadamard product**: $A \circ B = \begin{pmatrix} 2 & 0.5 \\ 0.5 & 6 \end{pmatrix}$.

    Positive definiteness check: $2 > 0$, $\det = 12 - 0.25 = 11.75 > 0$, so $A \circ B \succ 0$.

    **Ordinary matrix product**: $AB = \begin{pmatrix} 2.5 & 3 \\ 2.5 & 6.5 \end{pmatrix}$, which is not symmetric. But by Theorem 16.7, all its eigenvalues are positive.

---

## 16.6 Hadamard inequality

<div class="context-flow" markdown>

**Intuition**: $\det A\leq\prod a_{ii}$, equality $\Leftrightarrow$ diagonal matrix · Two proofs: Cholesky decomposition / AM-GM inequality → special case of the Fischer inequality with $1\times 1$ blocks

</div>

!!! theorem "Theorem 16.10 (Hadamard inequality)"
    Let $A = (a_{ij}) \in \mathbb{C}^{n \times n}$ be a positive definite Hermitian matrix. Then

    $$\det(A) \leq \prod_{i=1}^{n} a_{ii},$$

    with equality if and only if $A$ is a diagonal matrix.

??? proof "Proof"
    **Proof 1 (via Cholesky decomposition)**: Let $A = LL^*$, $L = (\ell_{ij})$. Then $a_{ii} = \sum_{k=1}^{i}|\ell_{ik}|^2 \geq |\ell_{ii}|^2 = \ell_{ii}^2$.

    $\det(A) = (\det L)^2 = \left(\prod_{i=1}^{n}\ell_{ii}\right)^2 = \prod_{i=1}^{n}\ell_{ii}^2 \leq \prod_{i=1}^{n}a_{ii}$.

    Equality holds $\Leftrightarrow$ each $a_{ii} = \ell_{ii}^2$ $\Leftrightarrow$ $L$ is diagonal $\Leftrightarrow$ $A$ is diagonal.

    **Proof 2 (via the AM-GM inequality)**: Let $D = \operatorname{diag}(a_{11}, \ldots, a_{nn})$; then $D^{-1/2}AD^{-1/2}$ is a positive definite matrix with all diagonal entries equal to $1$. Let its eigenvalues be $\mu_1, \ldots, \mu_n > 0$. Since $\operatorname{tr}(D^{-1/2}AD^{-1/2}) = n$, i.e., $\sum\mu_i = n$. By the AM-GM inequality:

    $$\frac{\det(A)}{\prod a_{ii}} = \det(D^{-1/2}AD^{-1/2}) = \prod\mu_i \leq \left(\frac{\sum\mu_i}{n}\right)^n = 1.$$

    Equality holds if and only if $\mu_1 = \cdots = \mu_n = 1$, i.e., $D^{-1/2}AD^{-1/2} = I$, i.e., $A = D$. $\blacksquare$

!!! example "Example 16.8"
    Let $A = \begin{pmatrix} 4 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 6 \end{pmatrix}$ (already verified to be positive definite).

    $\det(A) = 4(30 - 9) - 2(12 - 3) + 1(6 - 5) = 84 - 18 + 1 = 67$.

    $\prod a_{ii} = 4 \times 5 \times 6 = 120$.

    $67 \leq 120$, the Hadamard inequality holds. The ratio $67/120 \approx 0.558$ reflects the degree to which $A$ deviates from a diagonal matrix.

---

## 16.7 Fischer inequality

<div class="context-flow" markdown>

**Chapter arc**: Block generalization of the Hadamard inequality · $\det A\leq\det A_{11}\cdot\det A_{22}$, proved via Schur complement · Equality $\Leftrightarrow$ off-diagonal block is zero → Ch18 determinantal inequalities

</div>

The Fischer inequality is a generalization of the Hadamard inequality, replacing diagonal entries with determinants of diagonal block submatrices.

!!! theorem "Theorem 16.11 (Fischer inequality)"
    Let the positive definite Hermitian matrix $A$ be partitioned as

    $$A = \begin{pmatrix} A_{11} & A_{12} \\ A_{12}^* & A_{22} \end{pmatrix},$$

    where $A_{11}$ is $k \times k$ and $A_{22}$ is $(n-k) \times (n-k)$. Then

    $$\det(A) \leq \det(A_{11})\cdot\det(A_{22}),$$

    with equality if and only if $A_{12} = O$ (i.e., $A$ is block diagonal).

??? proof "Proof"
    By the Schur complement determinant formula:

    $$\det(A) = \det(A_{11})\cdot\det(A_{22} - A_{12}^*A_{11}^{-1}A_{12}).$$

    By positive definiteness, the Schur complement $S = A_{22} - A_{12}^*A_{11}^{-1}A_{12} \succ 0$. Also $A_{22} - S = A_{12}^*A_{11}^{-1}A_{12} \succeq 0$, i.e., $A_{22} \succeq S$.

    Therefore $A_{22} - S \succeq 0$, and its eigenvalues are non-negative. Since $A_{22} \succeq S \succ 0$ implies $S^{-1/2}A_{22}S^{-1/2} \succeq I$, we get $\det(S^{-1/2}A_{22}S^{-1/2}) \geq 1$, i.e., $\det(A_{22}) \geq \det(S)$.

    Therefore $\det(A) = \det(A_{11})\det(S) \leq \det(A_{11})\det(A_{22})$.

    Equality holds if and only if $S = A_{22}$, i.e., $A_{12}^*A_{11}^{-1}A_{12} = O$. Since $A_{11} \succ 0$, this is equivalent to $A_{12} = O$. $\blacksquare$

!!! corollary "Corollary 16.1"
    For a positive definite matrix $A$, more generally one can partition $A$ into multiple diagonal blocks $A_{11}, A_{22}, \ldots, A_{kk}$; then

    $$\det(A) \leq \prod_{i=1}^{k}\det(A_{ii}).$$

    The Hadamard inequality is the special case $k = n$ (each block is $1 \times 1$).

!!! example "Example 16.9"
    Let $A = \begin{pmatrix} 2 & 1 & 0.5 & 0.1 \\ 1 & 3 & 0.3 & 0.2 \\ 0.5 & 0.3 & 4 & 1 \\ 0.1 & 0.2 & 1 & 5 \end{pmatrix}$, with blocks $A_{11} = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix}$, $A_{22} = \begin{pmatrix} 4 & 1 \\ 1 & 5 \end{pmatrix}$.

    $\det(A_{11}) = 5$, $\det(A_{22}) = 19$.

    The Fischer inequality gives $\det(A) \leq 5 \times 19 = 95$.

    The Hadamard inequality gives $\det(A) \leq 2 \times 3 \times 4 \times 5 = 120$.

    The Fischer inequality gives a tighter bound.

---

## 16.8 Positive definite matrices and optimization

<div class="context-flow" markdown>

**Chapter arc**: Hessian $A\succ 0$ $\Leftrightarrow$ strictly convex $\Leftrightarrow$ unique minimum · Level sets = ellipsoids (semi-axes $1/\sqrt{\lambda_i}$) · Positive definiteness test $\to$ global optimality guarantee in convex optimization

</div>

Positive definite matrices play a central role in optimization theory. The convexity of a quadratic function is determined by the positive definiteness of its Hessian matrix, and the level sets of positive definite matrices are ellipsoids.

!!! definition "Definition 16.8 (Quadratic function)"
    Let $A$ be an $n \times n$ real symmetric matrix, $\mathbf{b} \in \mathbb{R}^n$, $c \in \mathbb{R}$. A **quadratic function** is defined as

    $$f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - \mathbf{b}^T\mathbf{x} + c.$$

!!! theorem "Theorem 16.12 (Extrema of quadratic functions)"
    Let $f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - \mathbf{b}^T\mathbf{x} + c$, where $A$ is symmetric.

    (1) $f$ has a unique global minimum if and only if $A \succ 0$. In this case the minimizer is $\mathbf{x}^* = A^{-1}\mathbf{b}$ and the minimum value is $f(\mathbf{x}^*) = c - \frac{1}{2}\mathbf{b}^TA^{-1}\mathbf{b}$.

    (2) $f$ is convex if and only if $A \succeq 0$.

??? proof "Proof"
    $\nabla f(\mathbf{x}) = A\mathbf{x} - \mathbf{b}$, Hessian matrix $\nabla^2 f = A$.

    **(1)** If $A \succ 0$, setting $\nabla f = \mathbf{0}$ gives $\mathbf{x}^* = A^{-1}\mathbf{b}$. Since $\nabla^2 f = A \succ 0$, $\mathbf{x}^*$ is a strict local minimum. By convexity ($A \succ 0$ implies strict convexity), it is also the unique global minimum.

    $f(\mathbf{x}) = f(\mathbf{x}^*) + \frac{1}{2}(\mathbf{x}-\mathbf{x}^*)^TA(\mathbf{x}-\mathbf{x}^*)$ (completing the square), where

    $f(\mathbf{x}^*) = \frac{1}{2}\mathbf{b}^TA^{-1}AA^{-1}\mathbf{b} - \mathbf{b}^TA^{-1}\mathbf{b} + c = -\frac{1}{2}\mathbf{b}^TA^{-1}\mathbf{b} + c$.

    **(2)** $f$ is convex $\Leftrightarrow$ for all $\mathbf{x}$, $\nabla^2f(\mathbf{x}) = A \succeq 0$. $\blacksquare$

!!! definition "Definition 16.9 (Ellipsoid)"
    Let $A \succ 0$ and $\mathbf{c} \in \mathbb{R}^n$. The **ellipsoid** centered at $\mathbf{c}$ with "shape matrix" $A$ is defined as

    $$\mathcal{E}(A, \mathbf{c}) = \{\mathbf{x} \in \mathbb{R}^n : (\mathbf{x} - \mathbf{c})^TA(\mathbf{x} - \mathbf{c}) \leq 1\}.$$

    The semi-axis lengths of the ellipsoid are the reciprocal square roots of the eigenvalues, $1/\sqrt{\lambda_i}$, and the directions are the corresponding eigenvectors. The volume is $\operatorname{vol}(\mathcal{E}) = \frac{\omega_n}{\sqrt{\det A}}$, where $\omega_n$ is the volume of the $n$-dimensional unit ball.

!!! theorem "Theorem 16.13 (Positive definite matrices and convex optimization)"
    Consider the constrained optimization problem $\min_{\mathbf{x}} f(\mathbf{x})$, $g_i(\mathbf{x}) \leq 0$ ($i = 1, \ldots, m$). If $f$ and all $g_i$ are quadratic functions with the Hessian of $f$ being positive definite and the Hessians of $g_i$ being positive semidefinite, then this is a convex quadratic program and any local optimum is also a global optimum.

!!! example "Example 16.10"
    Find the minimum of $f(\mathbf{x}) = 3x_1^2 + 2x_1x_2 + 2x_2^2 - 4x_1 - 6x_2$.

    In matrix form: $f(\mathbf{x}) = \mathbf{x}^T\begin{pmatrix} 3 & 1 \\ 1 & 2 \end{pmatrix}\mathbf{x} - \begin{pmatrix} 4 \\ 6 \end{pmatrix}^T\mathbf{x}$.

    $A = \begin{pmatrix} 3 & 1 \\ 1 & 2 \end{pmatrix}$, eigenvalues $\frac{5 \pm \sqrt{5}}{2} > 0$, positive definite.

    $\mathbf{x}^* = A^{-1}\mathbf{b} = \frac{1}{5}\begin{pmatrix} 2 & -1 \\ -1 & 3 \end{pmatrix}\begin{pmatrix} 4 \\ 6 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 2 \\ 14 \end{pmatrix} = \begin{pmatrix} 2/5 \\ 14/5 \end{pmatrix}$.

    $f(\mathbf{x}^*) = -\frac{1}{2}\mathbf{b}^TA^{-1}\mathbf{b} = -\frac{1}{2}\cdot\frac{1}{5}(4,6)\begin{pmatrix} 2 \\ 14 \end{pmatrix} = -\frac{1}{10}(8+84) = -\frac{92}{10} = -9.2$.

!!! example "Example 16.11"
    **Geometric description of an ellipsoid.** Let $A = \begin{pmatrix} 4 & 0 \\ 0 & 1 \end{pmatrix}$, center $\mathbf{c} = \mathbf{0}$.

    Ellipsoid equation: $4x_1^2 + x_2^2 \leq 1$, i.e., $\frac{x_1^2}{(1/2)^2} + \frac{x_2^2}{1^2} \leq 1$.

    Semi-axis lengths: $1/\sqrt{4} = 1/2$ (along the $x_1$-axis) and $1/\sqrt{1} = 1$ (along the $x_2$-axis).

    Area $= \pi \times \frac{1}{2} \times 1 = \frac{\pi}{2} = \frac{\pi}{\sqrt{\det A}} = \frac{\pi}{\sqrt{4}} = \frac{\pi}{2}$.
