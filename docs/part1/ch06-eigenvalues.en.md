# Chapter 6  Eigenvalues and Eigenvectors

<div class="context-flow" markdown>

**Prerequisites**: Chapter 3 determinants · Chapter 5 similar matrices · invariant subspaces

**Chapter arc**: $A\mathbf{v}=\lambda\mathbf{v}$ → characteristic polynomial $\det(A-\lambda I)=0$ → algebraic/geometric multiplicity → eigenspaces → diagonalization conditions → spectral theorem for real symmetric matrices → Cayley-Hamilton → similarity invariants

**Further connections**：Eigenvalues directly apply to Google PageRank (Perron vector), quantum mechanics (energy eigenstates), vibration analysis (natural frequencies), and population models (Leslie matrices); the infinite-dimensional generalization is spectral theory, where the distinction between continuous and point spectra is central to quantum mechanics

</div>

Eigenvalues and eigenvectors are among the most profound and widely applied concepts in linear algebra. Intuitively, the eigenvectors of a linear transformation are those vectors whose "direction is unchanged" — the transformation only changes their length (scaling). This simple idea gives rise to a rich theory: characteristic polynomials, diagonalization, spectral theorems, and more. Eigenvalue theory holds a central place not only in pure mathematics but also finds extensive applications in quantum mechanics, vibration analysis, principal component analysis, the Google PageRank algorithm, and many other fields. This chapter systematically studies the definition, computation, properties, and applications of eigenvalues and eigenvectors.

---

## 6.1 Definition of Eigenvalues and Eigenvectors

<div class="context-flow" markdown>

**From invariant subspaces to eigenvalues**: Chapter 5 one-dimensional invariant subspace $\operatorname{span}\{\mathbf{v}\}$ → $T(\mathbf{v})=\lambda\mathbf{v}$ → when $\lambda=0$, $\mathbf{v} \in \ker A$ ($A$ is singular); the requirement $\mathbf{v} \neq \mathbf{0}$ is essential

</div>

!!! definition "Definition 6.1 (Eigenvalue and eigenvector)"
    Let $A$ be an $n \times n$ matrix (or $T: V \to V$ a linear operator). If there exist a scalar $\lambda \in \mathbb{F}$ and a **nonzero** vector $\mathbf{v} \in \mathbb{F}^n$ (or $\mathbf{v} \in V$, $\mathbf{v} \neq \mathbf{0}$) such that

    $$A\mathbf{v} = \lambda\mathbf{v}$$

    then $\lambda$ is called an **eigenvalue** of $A$, and $\mathbf{v}$ is called an **eigenvector** of $A$ corresponding to the eigenvalue $\lambda$.

!!! note "Note"
    In $A\mathbf{v} = \lambda\mathbf{v}$, the requirement $\mathbf{v} \neq \mathbf{0}$ is essential — the zero vector satisfies this equation for any $\lambda$ but provides no information. On the other hand, $\lambda = 0$ is allowed; in that case $A\mathbf{v} = \mathbf{0}$, i.e., $\mathbf{v} \in \ker(A)$. Therefore $0$ is an eigenvalue if and only if $A$ is singular.

!!! proposition "Proposition 6.1 (Geometric meaning of eigenvalues)"
    Let $\lambda$ be an eigenvalue of $A$ and $\mathbf{v}$ a corresponding eigenvector. Then:

    - If $\lambda > 1$: the direction of $\mathbf{v}$ is "stretched."
    - If $0 < \lambda < 1$: the direction of $\mathbf{v}$ is "compressed."
    - If $\lambda = 1$: the direction of $\mathbf{v}$ is unchanged ($\mathbf{v}$ is a fixed-point direction).
    - If $\lambda = 0$: the direction of $\mathbf{v}$ is mapped to zero.
    - If $\lambda < 0$: the direction of $\mathbf{v}$ is reversed and scaled.

!!! example "Example 6.1"
    Let $A = \begin{pmatrix} 3 & 1 \\ 0 & 2 \end{pmatrix}$. Take $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$. Then

    $$A\mathbf{v}_1 = \begin{pmatrix} 3 \\ 0 \end{pmatrix} = 3\begin{pmatrix} 1 \\ 0 \end{pmatrix} = 3\mathbf{v}_1$$

    So $\lambda_1 = 3$ is an eigenvalue and $\mathbf{v}_1$ is a corresponding eigenvector.

    Take $\mathbf{v}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$. Then

    $$A\mathbf{v}_2 = \begin{pmatrix} 2 \\ -2 \end{pmatrix} = 2\begin{pmatrix} 1 \\ -1 \end{pmatrix} = 2\mathbf{v}_2$$

    So $\lambda_2 = 2$ is an eigenvalue and $\mathbf{v}_2$ is a corresponding eigenvector.

---

## 6.2 The Characteristic Polynomial

<div class="context-flow" markdown>

**Reducing the eigenvalue problem to root-finding**: $A\mathbf{v}=\lambda\mathbf{v}$ → $(A-\lambda I)\mathbf{v}=\mathbf{0}$ has a nontrivial solution → $\det(A-\lambda I)=0$ (Chapter 3) → a degree-$n$ polynomial has $n$ complex roots

</div>

!!! definition "Definition 6.2 (Characteristic polynomial)"
    Let $A$ be an $n \times n$ matrix. Then

    $$p(\lambda) = \det(A - \lambda I)$$

    is called the **characteristic polynomial** of $A$. The eigenvalues of $A$ are precisely the roots of the characteristic polynomial.

!!! note "Note"
    Some textbooks define the characteristic polynomial as $\det(\lambda I - A)$; the two differ only by a sign factor of $(-1)^n$ and have the same roots. This book adopts the convention $\det(A - \lambda I)$.

!!! theorem "Theorem 6.1 (Eigenvalues and the characteristic polynomial)"
    $\lambda_0$ is an eigenvalue of $A$ if and only if $\det(A - \lambda_0 I) = 0$.

??? proof "Proof"
    $\lambda_0$ is an eigenvalue $\Leftrightarrow$ there exists a nonzero $\mathbf{v}$ with $A\mathbf{v} = \lambda_0 \mathbf{v}$ $\Leftrightarrow$ $(A - \lambda_0 I)\mathbf{v} = \mathbf{0}$ has a nontrivial solution $\Leftrightarrow$ $A - \lambda_0 I$ is singular $\Leftrightarrow$ $\det(A - \lambda_0 I) = 0$. $\blacksquare$

!!! theorem "Theorem 6.2 (Degree and leading term of the characteristic polynomial)"
    Let $A$ be an $n \times n$ matrix. Then $p(\lambda) = \det(A - \lambda I)$ is a polynomial of degree $n$ in $\lambda$:

    $$p(\lambda) = (-1)^n \lambda^n + (-1)^{n-1}(\operatorname{tr}A)\lambda^{n-1} + \cdots + \det A$$

    where $\operatorname{tr}A = a_{11} + a_{22} + \cdots + a_{nn}$ is the **trace** of $A$.

??? proof "Proof"
    When expanding the determinant $\det(A - \lambda I)$, the $\lambda^n$ term comes from the product of the main diagonal entries:

    $$(a_{11} - \lambda)(a_{22} - \lambda)\cdots(a_{nn} - \lambda) = (-\lambda)^n + (-\lambda)^{n-1}(a_{11} + \cdots + a_{nn}) + \cdots$$

    Other terms involve at most $n-2$ diagonal entry factors, contributing at most degree $n-2$ in $\lambda$. Hence the leading term is $(-1)^n \lambda^n$ and the coefficient of $\lambda^{n-1}$ is $(-1)^{n-1}\operatorname{tr}A$.

    The constant term is $p(0) = \det(A - 0 \cdot I) = \det A$. $\blacksquare$

!!! example "Example 6.2"
    Find the eigenvalues of $A = \begin{pmatrix} 4 & -2 \\ 1 & 1 \end{pmatrix}$.

    $$p(\lambda) = \det(A - \lambda I) = \det\begin{pmatrix} 4 - \lambda & -2 \\ 1 & 1 - \lambda \end{pmatrix} = (4-\lambda)(1-\lambda) + 2 = \lambda^2 - 5\lambda + 6 = (\lambda - 2)(\lambda - 3)$$

    The eigenvalues are $\lambda_1 = 2$ and $\lambda_2 = 3$.

    Verification: $\operatorname{tr}A = 4 + 1 = 5 = 2 + 3$, $\det A = 4 \cdot 1 - (-2) \cdot 1 = 6 = 2 \times 3$.

!!! example "Example 6.3"
    Find the eigenvalues of $A = \begin{pmatrix} 2 & 1 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 3 \end{pmatrix}$.

    $$p(\lambda) = \det(A - \lambda I) = (2-\lambda)(2-\lambda)(3-\lambda) = (2-\lambda)^2(3-\lambda)$$

    The eigenvalues are $\lambda_1 = 2$ (double root) and $\lambda_2 = 3$ (simple root).

---

## 6.3 Properties of Eigenvalues

<div class="context-flow" markdown>

**Tension between algebra and geometry**: $1 \le \operatorname{GM}(\lambda) \le \operatorname{AM}(\lambda)$ → diagonalizable when both are equal → $\operatorname{tr}A = \sum\lambda_i$, $\det A = \prod\lambda_i$ (connecting Chapter 3 determinants with eigenvalues)

</div>

!!! definition "Definition 6.3 (Algebraic multiplicity)"
    The multiplicity of eigenvalue $\lambda_0$ as a root of the characteristic polynomial $p(\lambda)$ is called the **algebraic multiplicity** of $\lambda_0$, denoted $\operatorname{AM}(\lambda_0)$.

!!! definition "Definition 6.4 (Geometric multiplicity)"
    The dimension of the eigenspace $E_{\lambda_0} = \ker(A - \lambda_0 I)$ corresponding to eigenvalue $\lambda_0$ is called the **geometric multiplicity** of $\lambda_0$, denoted $\operatorname{GM}(\lambda_0)$.

!!! theorem "Theorem 6.3 (Trace and determinant)"
    Let $A$ be an $n \times n$ matrix with eigenvalues (counting multiplicity) $\lambda_1, \lambda_2, \ldots, \lambda_n$. Then

    $$\operatorname{tr}A = \sum_{i=1}^n \lambda_i, \qquad \det A = \prod_{i=1}^n \lambda_i$$

??? proof "Proof"
    The characteristic polynomial can be written as

    $$p(\lambda) = (-1)^n(\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n)$$

    Expanding, the coefficient of $\lambda^{n-1}$ is $(-1)^n \cdot (-1)(\lambda_1 + \cdots + \lambda_n) = (-1)^{n-1}(\lambda_1 + \cdots + \lambda_n)$. By Theorem 6.2, this coefficient also equals $(-1)^{n-1}\operatorname{tr}A$, so $\operatorname{tr}A = \lambda_1 + \cdots + \lambda_n$.

    Setting $\lambda = 0$: $p(0) = (-1)^n(-\lambda_1)(-\lambda_2)\cdots(-\lambda_n) = \lambda_1\lambda_2\cdots\lambda_n$. Also $p(0) = \det A$. $\blacksquare$

!!! theorem "Theorem 6.4 (Relationship between geometric and algebraic multiplicity)"
    Let $\lambda_0$ be an eigenvalue of $A$. Then

    $$1 \leq \operatorname{GM}(\lambda_0) \leq \operatorname{AM}(\lambda_0)$$

??? proof "Proof"
    $\operatorname{GM}(\lambda_0) \geq 1$ because an eigenvalue must have a nonzero eigenvector, so the eigenspace has dimension at least $1$.

    Let $\operatorname{GM}(\lambda_0) = k$. Take a basis $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ of the eigenspace $E_{\lambda_0}$ and extend it to a basis $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{v}_{k+1}, \ldots, \mathbf{v}_n\}$ of $\mathbb{F}^n$.

    Let $P = (\mathbf{v}_1 \;\; \cdots \;\; \mathbf{v}_n)$. Then

    $$P^{-1}AP = \begin{pmatrix} \lambda_0 I_k & B \\ 0 & C \end{pmatrix}$$

    where $I_k$ is the $k \times k$ identity matrix. Similar matrices have the same characteristic polynomial:

    $$p(\lambda) = \det(P^{-1}AP - \lambda I) = (\lambda_0 - \lambda)^k \det(C - \lambda I)$$

    Therefore $(\lambda_0 - \lambda)^k$ divides $p(\lambda)$, i.e., $\operatorname{AM}(\lambda_0) \geq k = \operatorname{GM}(\lambda_0)$. $\blacksquare$

!!! theorem "Theorem 6.5 (Eigenvectors corresponding to distinct eigenvalues are linearly independent)"
    Let $\lambda_1, \lambda_2, \ldots, \lambda_k$ be distinct eigenvalues of $A$, and $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k$ corresponding eigenvectors. Then $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\}$ is linearly independent.

??? proof "Proof"
    By induction on $k$. When $k = 1$, $\mathbf{v}_1 \neq \mathbf{0}$ is obviously linearly independent.

    Assume the result holds for $k - 1$; consider the case of $k$. Suppose

    $$c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_k\mathbf{v}_k = \mathbf{0} \quad \cdots (*)$$

    Applying $A$ to both sides of $(*)$:

    $$c_1\lambda_1\mathbf{v}_1 + c_2\lambda_2\mathbf{v}_2 + \cdots + c_k\lambda_k\mathbf{v}_k = \mathbf{0} \quad \cdots (**)$$

    Subtracting $\lambda_k$ times $(*)$ from $(**)$:

    $$c_1(\lambda_1 - \lambda_k)\mathbf{v}_1 + c_2(\lambda_2 - \lambda_k)\mathbf{v}_2 + \cdots + c_{k-1}(\lambda_{k-1} - \lambda_k)\mathbf{v}_{k-1} = \mathbf{0}$$

    By the induction hypothesis, $\{\mathbf{v}_1, \ldots, \mathbf{v}_{k-1}\}$ is linearly independent, and $\lambda_i - \lambda_k \neq 0$ ($i < k$), so $c_i(\lambda_i - \lambda_k) = 0$, i.e., $c_i = 0$ ($i = 1, \ldots, k-1$). Substituting back into $(*)$ gives $c_k\mathbf{v}_k = \mathbf{0}$, and since $\mathbf{v}_k \neq \mathbf{0}$, $c_k = 0$. $\blacksquare$

!!! example "Example 6.4"
    Let $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$. The characteristic polynomial is $p(\lambda) = \lambda^2 + 1$.

    Over $\mathbb{R}$ there are no real roots, so $A$ has no real eigenvalues. Over $\mathbb{C}$, however, the eigenvalues are $\lambda_1 = i$ and $\lambda_2 = -i$.

    Verification: $\operatorname{tr}A = 0 = i + (-i)$, $\det A = 1 = i \cdot (-i)$.

    This shows that a real matrix may have no real eigenvalues, but over the complex field, an $n \times n$ matrix must have $n$ eigenvalues (counting multiplicity, by the fundamental theorem of algebra).

---

## 6.4 Eigenspaces

<div class="context-flow" markdown>

**Eigenspace = $\ker(A-\lambda I)$**: the solution space of the homogeneous system from Chapter 1 → its dimension = geometric multiplicity → eigenspaces for distinct eigenvalues form a **direct sum** (Chapter 4)

</div>

!!! definition "Definition 6.5 (Eigenspace)"
    Let $\lambda_0$ be an eigenvalue of $A$. Then

    $$E_{\lambda_0} = \ker(A - \lambda_0 I) = \{\mathbf{v} \in \mathbb{F}^n : A\mathbf{v} = \lambda_0\mathbf{v}\}$$

    is called the **eigenspace** of $A$ corresponding to $\lambda_0$.

!!! theorem "Theorem 6.6 (The eigenspace is a subspace)"
    $E_{\lambda_0}$ is a subspace of $\mathbb{F}^n$.

??? proof "Proof"
    $E_{\lambda_0} = \ker(A - \lambda_0 I)$ is the solution space of a homogeneous linear system. By the results of Chapter 4, it is a subspace of $\mathbb{F}^n$. Alternatively, one can verify directly: $\mathbf{0} \in E_{\lambda_0}$ (although $\mathbf{0}$ is not an eigenvector, it belongs to the eigenspace); if $\mathbf{u}, \mathbf{v} \in E_{\lambda_0}$, then $A(\mathbf{u} + \mathbf{v}) = A\mathbf{u} + A\mathbf{v} = \lambda_0\mathbf{u} + \lambda_0\mathbf{v} = \lambda_0(\mathbf{u} + \mathbf{v})$; if $c \in \mathbb{F}$, then $A(c\mathbf{v}) = cA\mathbf{v} = c\lambda_0\mathbf{v} = \lambda_0(c\mathbf{v})$. $\blacksquare$

!!! example "Example 6.5"
    Let $A = \begin{pmatrix} 5 & -2 & 0 \\ 0 & 3 & 0 \\ 4 & -2 & 3 \end{pmatrix}$. Find the eigenvalues and eigenspaces.

    **Characteristic polynomial**:

    $$p(\lambda) = \det(A - \lambda I) = \det\begin{pmatrix} 5-\lambda & -2 & 0 \\ 0 & 3-\lambda & 0 \\ 4 & -2 & 3-\lambda \end{pmatrix}$$

    Expanding along the third column: $p(\lambda) = (3 - \lambda)\det\begin{pmatrix} 5-\lambda & -2 \\ 0 & 3-\lambda \end{pmatrix} = (3-\lambda)^2(5-\lambda)$.

    The eigenvalues are $\lambda_1 = 3$ ($\operatorname{AM} = 2$) and $\lambda_2 = 5$ ($\operatorname{AM} = 1$).

    **Eigenspace for $\lambda_1 = 3$**:

    $$A - 3I = \begin{pmatrix} 2 & -2 & 0 \\ 0 & 0 & 0 \\ 4 & -2 & 0 \end{pmatrix} \xrightarrow{\text{row reduction}} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$$

    Let us compute carefully. $R_3 - 2R_1$: $\begin{pmatrix} 2 & -2 & 0 \\ 0 & 0 & 0 \\ 0 & 2 & 0 \end{pmatrix}$, then $R_1/2$: $\begin{pmatrix} 1 & -1 & 0 \\ 0 & 0 & 0 \\ 0 & 2 & 0 \end{pmatrix}$, $R_3/2$: $\begin{pmatrix} 1 & -1 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 0 \end{pmatrix}$, swap $R_2, R_3$: $\begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$, and finally $\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$.

    Hence $x_1 = x_2 = 0$, $x_3 = t$ is free. $E_3 = \operatorname{span}\left\{\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}\right\}$, $\operatorname{GM}(3) = 1 < 2 = \operatorname{AM}(3)$.

    **Eigenspace for $\lambda_2 = 5$**:

    $$A - 5I = \begin{pmatrix} 0 & -2 & 0 \\ 0 & -2 & 0 \\ 4 & -2 & -2 \end{pmatrix} \xrightarrow{\text{row reduction}} \begin{pmatrix} 1 & 0 & -\frac{1}{2} \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$$

    Hence $x_1 = \frac{1}{2}t$, $x_2 = 0$, $x_3 = t$. $E_5 = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 0 \\ 2 \end{pmatrix}\right\}$, $\operatorname{GM}(5) = 1 = \operatorname{AM}(5)$.

---

## 6.5 Diagonalization

<div class="context-flow" markdown>

**One of the central goals of linear algebra**: $A = PDP^{-1}$ → columns of $P$ are eigenvectors, diagonal entries of $D$ are eigenvalues → diagonalizable $\Leftrightarrow$ $n$ linearly independent eigenvectors $\Leftrightarrow$ $\forall i:\operatorname{GM}=\operatorname{AM}$

</div>

!!! definition "Definition 6.6 (Diagonalizable)"
    An $n \times n$ matrix $A$ is called **diagonalizable** if there exist an invertible matrix $P$ and a diagonal matrix $D$ such that

    $$A = PDP^{-1}$$

    Equivalently, $D = P^{-1}AP$, i.e., $A$ is similar to a diagonal matrix.

!!! theorem "Theorem 6.7 (Necessary and sufficient conditions for diagonalization)"
    An $n \times n$ matrix $A$ is diagonalizable if and only if one of the following equivalent conditions holds:

    1. $A$ has $n$ linearly independent eigenvectors.
    2. For each eigenvalue $\lambda_i$ of $A$, $\operatorname{GM}(\lambda_i) = \operatorname{AM}(\lambda_i)$.
    3. The direct sum of eigenspaces equals the full space: $\bigoplus_i E_{\lambda_i} = \mathbb{F}^n$.

??? proof "Proof"
    **$A$ is diagonalizable $\Leftrightarrow$ condition 1**:

    $(\Rightarrow)$ If $A = PDP^{-1}$, let $P = (\mathbf{p}_1 \;\; \cdots \;\; \mathbf{p}_n)$ and $D = \operatorname{diag}(d_1, \ldots, d_n)$. Then $AP = PD$, i.e., $A\mathbf{p}_j = d_j \mathbf{p}_j$. So $\mathbf{p}_1, \ldots, \mathbf{p}_n$ are eigenvectors of $A$, and they are linearly independent since $P$ is invertible.

    $(\Leftarrow)$ Let $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ be $n$ linearly independent eigenvectors with $A\mathbf{v}_j = \lambda_j\mathbf{v}_j$. Set $P = (\mathbf{v}_1 \;\; \cdots \;\; \mathbf{v}_n)$ and $D = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$. Then $AP = PD$, i.e., $A = PDP^{-1}$.

    **Condition 1 $\Leftrightarrow$ condition 2**: $A$ has $n$ linearly independent eigenvectors $\Leftrightarrow$ $\sum_i \operatorname{GM}(\lambda_i) = n = \sum_i \operatorname{AM}(\lambda_i)$. Since $\operatorname{GM}(\lambda_i) \leq \operatorname{AM}(\lambda_i)$ and the sums are equal, we must have $\operatorname{GM}(\lambda_i) = \operatorname{AM}(\lambda_i)$ for each $i$. $\blacksquare$

!!! corollary "Corollary 6.1"
    If an $n \times n$ matrix $A$ has $n$ distinct eigenvalues, then $A$ is diagonalizable.

??? proof "Proof"
    By Theorem 6.5, eigenvectors corresponding to distinct eigenvalues are linearly independent; $n$ distinct eigenvalues yield $n$ linearly independent eigenvectors. By condition 1 of Theorem 6.7, $A$ is diagonalizable. $\blacksquare$

!!! example "Example 6.6"
    Diagonalize $A = \begin{pmatrix} 4 & -2 \\ 1 & 1 \end{pmatrix}$.

    From Example 6.2, $\lambda_1 = 2$, $\lambda_2 = 3$.

    $\lambda_1 = 2$: $(A - 2I)\mathbf{v} = \mathbf{0}$, $\begin{pmatrix} 2 & -2 \\ 1 & -1 \end{pmatrix}\mathbf{v} = \mathbf{0}$, giving $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$.

    $\lambda_2 = 3$: $(A - 3I)\mathbf{v} = \mathbf{0}$, $\begin{pmatrix} 1 & -2 \\ 1 & -2 \end{pmatrix}\mathbf{v} = \mathbf{0}$, giving $\mathbf{v}_2 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$.

    Set $P = \begin{pmatrix} 1 & 2 \\ 1 & 1 \end{pmatrix}$, $D = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$. Then $A = PDP^{-1}$.

    **Application**: $A^k = PD^kP^{-1} = P\begin{pmatrix} 2^k & 0 \\ 0 & 3^k \end{pmatrix}P^{-1}$, enabling fast computation of high powers of the matrix.

!!! example "Example 6.7"
    $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ is not diagonalizable.

    $p(\lambda) = (2 - \lambda)^2$; the only eigenvalue is $\lambda = 2$, with $\operatorname{AM}(2) = 2$.

    $A - 2I = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$, $\operatorname{rank}(A - 2I) = 1$, so $\operatorname{GM}(2) = 2 - 1 = 1 < 2 = \operatorname{AM}(2)$.

    $A$ does not have two linearly independent eigenvectors, so it is not diagonalizable.

---

## 6.6 Eigenvalues of Real Symmetric Matrices

<div class="context-flow" markdown>

**The perfect properties of symmetric matrices**: eigenvalues are all real → eigenvectors for distinct eigenvalues are orthogonal → **spectral theorem**: $A=QDQ^T$ (orthogonal diagonalization) → connects to Chapter 7 orthogonal matrices

</div>

!!! definition "Definition 6.7 (Real symmetric matrix)"
    A real square matrix $A$ is called **symmetric** if $A^T = A$.

!!! theorem "Theorem 6.8 (Eigenvalues of a real symmetric matrix are real)"
    Let $A$ be a real symmetric matrix. Then all eigenvalues of $A$ are real.

??? proof "Proof"
    Let $\lambda$ be an eigenvalue of $A$ (possibly complex) and $\mathbf{v} \in \mathbb{C}^n$ a corresponding eigenvector. Then $A\mathbf{v} = \lambda\mathbf{v}$.

    Taking the conjugate transpose: $\overline{\mathbf{v}}^T A^T = \overline{\lambda} \overline{\mathbf{v}}^T$ (since $A$ is real, $\overline{A} = A$). Also $A^T = A$, so $\overline{\mathbf{v}}^T A = \overline{\lambda} \overline{\mathbf{v}}^T$.

    Computing $\overline{\mathbf{v}}^T A \mathbf{v}$:

    - From the left: $\overline{\mathbf{v}}^T (A\mathbf{v}) = \overline{\mathbf{v}}^T (\lambda\mathbf{v}) = \lambda (\overline{\mathbf{v}}^T \mathbf{v})$
    - From the right: $(\overline{\mathbf{v}}^T A)\mathbf{v} = \overline{\lambda} (\overline{\mathbf{v}}^T \mathbf{v})$

    Therefore $\lambda (\overline{\mathbf{v}}^T \mathbf{v}) = \overline{\lambda} (\overline{\mathbf{v}}^T \mathbf{v})$. Since $\overline{\mathbf{v}}^T \mathbf{v} = \sum |v_i|^2 > 0$ (because $\mathbf{v} \neq \mathbf{0}$), we get $\lambda = \overline{\lambda}$, i.e., $\lambda \in \mathbb{R}$. $\blacksquare$

!!! theorem "Theorem 6.9 (Eigenvectors of a real symmetric matrix for distinct eigenvalues are orthogonal)"
    Let $A$ be a real symmetric matrix and $\lambda_1 \neq \lambda_2$ two distinct eigenvalues of $A$, with $\mathbf{v}_1, \mathbf{v}_2$ the corresponding eigenvectors. Then $\mathbf{v}_1 \perp \mathbf{v}_2$ (i.e., $\mathbf{v}_1^T \mathbf{v}_2 = 0$).

??? proof "Proof"
    $\lambda_1(\mathbf{v}_1^T\mathbf{v}_2) = (\lambda_1\mathbf{v}_1)^T\mathbf{v}_2 = (A\mathbf{v}_1)^T\mathbf{v}_2 = \mathbf{v}_1^T A^T \mathbf{v}_2 = \mathbf{v}_1^T A \mathbf{v}_2 = \mathbf{v}_1^T(\lambda_2\mathbf{v}_2) = \lambda_2(\mathbf{v}_1^T\mathbf{v}_2)$

    Therefore $(\lambda_1 - \lambda_2)(\mathbf{v}_1^T\mathbf{v}_2) = 0$. Since $\lambda_1 \neq \lambda_2$, we get $\mathbf{v}_1^T\mathbf{v}_2 = 0$. $\blacksquare$

!!! theorem "Theorem 6.10 (Spectral theorem / Orthogonal diagonalization)"
    A real symmetric matrix $A$ can always be orthogonally diagonalized: there exist an orthogonal matrix $Q$ ($Q^TQ = I$) and a diagonal matrix $D$ such that

    $$A = QDQ^T$$

??? proof "Proof"
    (Outline) By the fundamental theorem of algebra, the characteristic polynomial of $A$ has $n$ roots in $\mathbb{C}$. By Theorem 6.8, all these roots are real.

    By induction on $n$. The case $n = 1$ is trivial. Assume the result holds for $n - 1$.

    Take an eigenvalue $\lambda_1$ of $A$ and a corresponding unit eigenvector $\mathbf{q}_1$. Extend $\mathbf{q}_1$ to an orthonormal basis $\{\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n\}$ of $\mathbb{R}^n$. Let $Q_1 = (\mathbf{q}_1 \;\; \cdots \;\; \mathbf{q}_n)$. Then

    $$Q_1^T A Q_1 = \begin{pmatrix} \lambda_1 & \mathbf{b}^T \\ \mathbf{0} & A_1 \end{pmatrix}$$

    Since $A$ is symmetric, $Q_1^T A Q_1$ is also symmetric, so $\mathbf{b} = \mathbf{0}$ and $A_1$ is an $(n-1) \times (n-1)$ real symmetric matrix. By the induction hypothesis, $A_1$ can be orthogonally diagonalized. From this one can construct the orthogonal diagonalization of $A$. $\blacksquare$

!!! example "Example 6.8"
    Orthogonally diagonalize $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$.

    $p(\lambda) = (2-\lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda - 1)(\lambda - 3)$.

    $\lambda_1 = 1$: $\mathbf{v}_1 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$, normalized to $\mathbf{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix}$.

    $\lambda_2 = 3$: $\mathbf{v}_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$, normalized to $\mathbf{q}_2 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}$.

    Verification: $\mathbf{q}_1^T\mathbf{q}_2 = 0$.

    $$Q = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}, \quad D = \begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}, \quad A = QDQ^T$$

---

## 6.7 The Cayley-Hamilton Theorem

<div class="context-flow" markdown>

**A matrix satisfies its own characteristic polynomial**: $p(A)=O$ → corollary: $A^{-1}$ can be expressed as a polynomial in $A$ → implies the powers of $A$ span a finite-dimensional space in the matrix algebra

</div>

!!! theorem "Theorem 6.11 (Cayley-Hamilton theorem)"
    Let $A$ be an $n \times n$ matrix and $p(\lambda) = \det(A - \lambda I)$ its characteristic polynomial. Then

    $$p(A) = O$$

    That is, replacing $\lambda$ by $A$ (with the constant term multiplied by $I$), the resulting matrix is the zero matrix.

!!! note "Note"
    The Cayley-Hamilton theorem does **not** mean simply substituting $\lambda = A$ into $\det(A - \lambda I) = 0$. The determinant is a scalar-valued function, while here we evaluate a **matrix polynomial**. The correct meaning is: if $p(\lambda) = c_0 + c_1\lambda + \cdots + c_n\lambda^n$, then $c_0 I + c_1 A + \cdots + c_n A^n = O$.

??? proof "Proof"
    Consider the matrix $B(\lambda) = \operatorname{adj}(A - \lambda I)$ (the adjugate of $A - \lambda I$). Each entry of $B(\lambda)$ is a polynomial of degree at most $n-1$ in $\lambda$, so it can be written as

    $$B(\lambda) = B_0 + B_1\lambda + B_2\lambda^2 + \cdots + B_{n-1}\lambda^{n-1}$$

    where $B_0, B_1, \ldots, B_{n-1}$ are constant matrices.

    By the property of the adjugate: $(A - \lambda I)B(\lambda) = \det(A - \lambda I) \cdot I = p(\lambda)I$.

    Let $p(\lambda) = c_0 + c_1\lambda + \cdots + c_n\lambda^n$. Expanding the left side and comparing coefficients of powers of $\lambda$:

    $$AB_0 = c_0 I$$
    $$AB_1 - B_0 = c_1 I$$
    $$AB_2 - B_1 = c_2 I$$
    $$\vdots$$
    $$AB_{n-1} - B_{n-2} = c_{n-1} I$$
    $$-B_{n-1} = c_n I$$

    Multiplying these equations on the left by $I, A, A^2, \ldots, A^{n-1}, A^n$ respectively and summing:

    $$AB_0 + A(AB_1 - B_0) + A^2(AB_2 - B_1) + \cdots + A^n(-B_{n-1})$$
    $$= c_0 I + c_1 A + c_2 A^2 + \cdots + c_n A^n = p(A)$$

    On the left side, all terms involving $B_i$ cancel pairwise (telescoping sum), yielding $O$. Therefore $p(A) = O$. $\blacksquare$

!!! example "Example 6.9"
    Verify the Cayley-Hamilton theorem for $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.

    $p(\lambda) = \det(A - \lambda I) = (1-\lambda)(4-\lambda) - 6 = \lambda^2 - 5\lambda - 2$.

    $$A^2 = \begin{pmatrix} 7 & 10 \\ 15 & 22 \end{pmatrix}$$

    $$p(A) = A^2 - 5A - 2I = \begin{pmatrix} 7 & 10 \\ 15 & 22 \end{pmatrix} - \begin{pmatrix} 5 & 10 \\ 15 & 20 \end{pmatrix} - \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$$

    Verified.

!!! corollary "Corollary 6.2"
    Let $A$ be an $n \times n$ invertible matrix with characteristic polynomial $p(\lambda) = c_0 + c_1\lambda + \cdots + c_n\lambda^n$. Since $A$ is invertible, $c_0 = p(0) = \det A \neq 0$. By the Cayley-Hamilton theorem:

    $$c_0 I + c_1 A + \cdots + c_n A^n = O$$

    Multiplying both sides by $A^{-1}$:

    $$A^{-1} = -\frac{1}{c_0}(c_1 I + c_2 A + \cdots + c_n A^{n-1})$$

    That is, $A^{-1}$ can be expressed as a polynomial in $A$.

!!! example "Example 6.10"
    Use Corollary 6.2 to find the inverse of $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.

    $p(\lambda) = \lambda^2 - 5\lambda - 2$, so $c_0 = -2$, $c_1 = -5$, $c_2 = 1$.

    $$A^{-1} = -\frac{1}{-2}(-5I + A) = \frac{1}{2}\begin{pmatrix} -4 & 2 \\ 3 & -1 \end{pmatrix} = \begin{pmatrix} -2 & 1 \\ 3/2 & -1/2 \end{pmatrix}$$

    Verification: $AA^{-1} = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} -2 & 1 \\ 3/2 & -1/2 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$. Correct.

---

## 6.8 Properties of Similar Matrices

<div class="context-flow" markdown>

**Similarity = different representations of the same linear operator** (Chapter 5) → shared: characteristic polynomial, trace, determinant, rank → but having the same eigenvalues does not guarantee similarity (the Jordan structure matters)

</div>

!!! definition "Definition 6.8 (Similarity invariant)"
    If a property or quantity is unchanged under the similarity transformation $A \mapsto P^{-1}AP$, it is called a **similarity invariant**.

!!! theorem "Theorem 6.12 (Basic properties of similar matrices)"
    If $A \sim B$ ($B = P^{-1}AP$), then:

    1. $\det A = \det B$
    2. $\operatorname{tr} A = \operatorname{tr} B$
    3. $\operatorname{rank} A = \operatorname{rank} B$
    4. $A$ and $B$ have the same characteristic polynomial (hence the same eigenvalues, counting multiplicity)
    5. $A$ is invertible if and only if $B$ is invertible
    6. $A^k \sim B^k$ for all positive integers $k$

??? proof "Proof"
    **1.** $\det B = \det(P^{-1}AP) = \det(P^{-1})\det(A)\det(P) = \frac{1}{\det P} \cdot \det A \cdot \det P = \det A$.

    **2.** Using the cyclic property of the trace $\operatorname{tr}(XY) = \operatorname{tr}(YX)$: $\operatorname{tr}B = \operatorname{tr}(P^{-1}AP) = \operatorname{tr}(APP^{-1}) = \operatorname{tr}A$.

    **3.** $P^{-1}$ and $P$ are both invertible; the similarity transformation does not change rank.

    **4.** $\det(B - \lambda I) = \det(P^{-1}AP - \lambda P^{-1}IP) = \det(P^{-1}(A - \lambda I)P) = \det(A - \lambda I)$.

    **5.** By 1, $\det A = \det B$, so $A$ is invertible $\Leftrightarrow$ $\det A \neq 0$ $\Leftrightarrow$ $\det B \neq 0$ $\Leftrightarrow$ $B$ is invertible.

    **6.** $B^k = (P^{-1}AP)^k = P^{-1}A^k P$. $\blacksquare$

!!! note "Note"
    Similarity is a weaker relation than equality — similar matrices share many important properties but are not necessarily equal. Similarity invariants include, but are not limited to: eigenvalues, characteristic polynomial, trace, determinant, rank, minimal polynomial, Jordan normal form, etc. However, **having the same eigenvalues does not necessarily imply similarity**. For example, $\begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$ and $\begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ have the same eigenvalues but are not similar.

!!! example "Example 6.11"
    Let $A = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}$, $B = \begin{pmatrix} 3 & -4 \\ 1 & -1 \end{pmatrix}$. Determine whether $A$ and $B$ are similar.

    $p_A(\lambda) = (1 - \lambda)(3 - \lambda) = \lambda^2 - 4\lambda + 3$.

    $p_B(\lambda) = (3 - \lambda)(-1 - \lambda) + 4 = \lambda^2 - 2\lambda + 1 = (\lambda - 1)^2$.

    $A$ has eigenvalues $1, 3$; $B$ has eigenvalues $1, 1$. Their characteristic polynomials differ, so $A \not\sim B$.

!!! example "Example 6.12"
    Let $A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$, $B = \begin{pmatrix} 2 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$.

    Both $A$ and $B$ are diagonal matrices with the same characteristic polynomial $(1 - \lambda)^2(2 - \lambda)$.

    Are they similar? Since diagonal matrices are similar if and only if their diagonal entries are the same multiset (order does not matter), $A$ and $B$ are indeed similar. In fact, taking the permutation matrix $P = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix}$ (swapping the 1st and 3rd rows/columns), we get $P^{-1}AP = B$.

---

## Chapter Summary

This chapter systematically studied the theory of eigenvalues and eigenvectors. The main topics include:

1. **Definition of eigenvalues and eigenvectors**: $A\mathbf{v} = \lambda\mathbf{v}$ reveals the "invariant directions" of a linear transformation.
2. **Characteristic polynomial**: eigenvalues are found via $\det(A - \lambda I) = 0$, reducing an algebraic problem to polynomial root-finding.
3. **Properties of eigenvalues**: the trace equals the sum of eigenvalues, the determinant equals the product of eigenvalues; algebraic multiplicity is no less than geometric multiplicity.
4. **Eigenspaces**: all eigenvectors belonging to the same eigenvalue (together with the zero vector) form a subspace.
5. **Diagonalization**: a matrix is diagonalizable if and only if the geometric multiplicity equals the algebraic multiplicity for each eigenvalue.
6. **Real symmetric matrices**: eigenvalues are necessarily real, eigenvectors for distinct eigenvalues are orthogonal, and orthogonal diagonalization is always possible (spectral theorem).
7. **Cayley-Hamilton theorem**: every square matrix satisfies its own characteristic polynomial.
8. **Similarity invariants**: similar matrices share the characteristic polynomial, trace, determinant, rank, and other properties.
