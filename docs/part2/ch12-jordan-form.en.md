# Chapter 12  Jordan Normal Form

<div class="context-flow" markdown>

**Prerequisites**: Ch8 Spectral theorem (diagonalizable case) · **Chapter arc**: Invariant subspaces → Generalized eigenvectors → Nilpotent matrices → Jordan blocks $J_k(\lambda) = \lambda I + N_k$ → **Jordan normal form theorem** → Minimal polynomial → Matrix powers/ODEs
Essence: The ultimate generalization of diagonalization — every matrix over the complex numbers has a unique "nearly diagonal" canonical form, whose structure is determined by the **minimal polynomial**

</div>

Not all matrices can be diagonalized. When the geometric multiplicity of an eigenvalue is less than its algebraic multiplicity, there is no basis consisting of eigenvectors, and diagonalization is impossible. However, every square matrix can be reduced to a "nearly diagonal" canonical form — the **Jordan normal form** (Jordan canonical form). The Jordan normal form is a natural generalization of diagonalization and plays a central role in matrix theory, differential equations, control theory, and other fields. This chapter builds the Jordan normal form theory step by step, starting from invariant subspaces.

---

## 12.1 Invariant Subspaces

<div class="context-flow" markdown>

$A(W) \subseteq W$ → Direct sum decomposition into invariant subspaces ↔ **Block diagonalization** → Decompose a large matrix problem into independently solvable small blocks

</div>

An invariant subspace is a fundamental tool for understanding the structure of linear transformations and the starting point of the Jordan decomposition theory.

!!! definition "Definition 12.1 (Invariant Subspace)"
    Let $A$ be an $n \times n$ matrix (or $T$ a linear transformation on a vector space $V$), and let $W$ be a subspace of $\mathbb{R}^n$ (or $V$). If $A\mathbf{w} \in W$ for every $\mathbf{w} \in W$ (i.e., $A(W) \subseteq W$), then $W$ is called an **invariant subspace** of $A$.

!!! example "Example 12.1"
    Let $A$ be an $n \times n$ matrix. The following are all invariant subspaces of $A$:

    1. $\{0\}$ and $\mathbb{R}^n$ (trivial invariant subspaces);
    2. $\ker(A) = \{\mathbf{x} : A\mathbf{x} = \mathbf{0}\}$ (null space);
    3. $\operatorname{Im}(A) = \{A\mathbf{x} : \mathbf{x} \in \mathbb{R}^n\}$ (image / column space);
    4. Eigenspace $E_\lambda = \ker(A - \lambda I)$;
    5. $\operatorname{span}\{\mathbf{v}\}$, where $\mathbf{v}$ is an eigenvector of $A$.

    **Verification of (2):** If $A\mathbf{x} = \mathbf{0}$, then $A(A\mathbf{x}) = A\mathbf{0} = \mathbf{0}$, so $A\mathbf{x} = \mathbf{0} \in \ker(A)$.

!!! definition "Definition 12.2 (Direct sum decomposition into invariant subspaces)"
    Let $W_1, W_2, \ldots, W_k$ all be invariant subspaces of $A$. If

    $$
    \mathbb{R}^n = W_1 \oplus W_2 \oplus \cdots \oplus W_k
    $$

    (i.e., every $\mathbf{x} \in \mathbb{R}^n$ can be uniquely written as $\mathbf{x} = \mathbf{w}_1 + \cdots + \mathbf{w}_k$ with $\mathbf{w}_i \in W_i$), then $\mathbb{R}^n$ is called a **direct sum decomposition** of invariant subspaces of $A$.

!!! theorem "Theorem 12.1 (Invariant subspaces and block diagonalization)"
    Let $\mathbb{R}^n = W_1 \oplus W_2 \oplus \cdots \oplus W_k$ be a direct sum decomposition of invariant subspaces of $A$. Choose a basis $\mathcal{B}_i$ for $W_i$ and combine them into $\mathcal{B} = \mathcal{B}_1 \cup \cdots \cup \mathcal{B}_k$. Let $P$ be the matrix of basis $\mathcal{B}$ (with columns being the basis vectors). Then

    $$
    P^{-1}AP = \begin{pmatrix} A_1 & & \\ & A_2 & \\ & & \ddots & \\ & & & A_k \end{pmatrix},
    $$

    where $A_i$ is the matrix representation of $A$ restricted to $W_i$ ($\dim W_i \times \dim W_i$).

??? proof "Proof"
    Let $\dim W_i = n_i$ and $\mathcal{B}_i = \{\mathbf{w}_{i,1}, \ldots, \mathbf{w}_{i,n_i}\}$. Since $W_i$ is an invariant subspace of $A$, $A\mathbf{w}_{i,j} \in W_i$, so $A\mathbf{w}_{i,j}$ can be expressed as a linear combination of vectors in $\mathcal{B}_i$.

    Let $P = [\mathbf{w}_{1,1}, \ldots, \mathbf{w}_{1,n_1}, \ldots, \mathbf{w}_{k,1}, \ldots, \mathbf{w}_{k,n_k}]$. Then

    $$
    AP = P \begin{pmatrix} A_1 & & \\ & \ddots & \\ & & A_k \end{pmatrix},
    $$

    where the $j$-th column of $A_i$ is the coordinate vector of $A\mathbf{w}_{i,j}$ with respect to the basis $\mathcal{B}_i$. Since $P$ is invertible, the result follows. $\blacksquare$

!!! theorem "Theorem 12.2 (Eigenspaces are invariant subspaces)"
    Let $\lambda$ be an eigenvalue of matrix $A$. Then the eigenspace $E_\lambda = \ker(A - \lambda I)$ is an invariant subspace of $A$. More generally, for any polynomial $p$, $\ker(p(A))$ is an invariant subspace of $A$.

??? proof "Proof"
    Let $\mathbf{v} \in E_\lambda$, i.e., $A\mathbf{v} = \lambda\mathbf{v}$. Then $A(\lambda\mathbf{v}) = \lambda(A\mathbf{v}) = \lambda^2\mathbf{v}$. But what we need to show is $A\mathbf{v} \in E_\lambda$: $(A - \lambda I)(A\mathbf{v}) = A(A\mathbf{v}) - \lambda(A\mathbf{v}) = A(\lambda\mathbf{v}) - \lambda(\lambda\mathbf{v}) = \lambda^2\mathbf{v} - \lambda^2\mathbf{v} = \mathbf{0}$. Therefore $A\mathbf{v} \in E_\lambda$.

    For the general case: let $p(A)\mathbf{v} = \mathbf{0}$. Since $A$ commutes with $p(A)$, $p(A)(A\mathbf{v}) = A(p(A)\mathbf{v}) = A\mathbf{0} = \mathbf{0}$, so $A\mathbf{v} \in \ker(p(A))$. $\blacksquare$

---

## 12.2 Generalized Eigenvectors

<div class="context-flow" markdown>

Not enough eigenvectors → relax to $(A-\lambda I)^k\mathbf{v} = \mathbf{0}$ → generalized eigenspace $G_\lambda$ has dimension = **algebraic multiplicity** → $\mathbb{C}^n = \bigoplus G_{\lambda_i}$

</div>

When the geometric multiplicity of an eigenvalue is less than its algebraic multiplicity, there are not enough eigenvectors to form a basis. Generalized eigenvectors compensate for this deficiency by relaxing the condition.

!!! definition "Definition 12.3 (Generalized Eigenvector)"
    Let $\lambda$ be an eigenvalue of matrix $A$. A nonzero vector $\mathbf{v}$ satisfying

    $$
    (A - \lambda I)^k \mathbf{v} = \mathbf{0}
    $$

    for some positive integer $k$ is called a **generalized eigenvector** of $A$ corresponding to $\lambda$ (of rank $k$, if $(A-\lambda I)^{k-1}\mathbf{v} \neq \mathbf{0}$).

!!! definition "Definition 12.4 (Generalized Eigenspace)"
    The **generalized eigenspace** of matrix $A$ corresponding to eigenvalue $\lambda$ is

    $$
    G_\lambda = \ker(A - \lambda I)^n = \{\mathbf{v} \in \mathbb{R}^n : (A - \lambda I)^n \mathbf{v} = \mathbf{0}\},
    $$

    where $n$ is the order of the matrix. Equivalently, $G_\lambda = \bigcup_{k=1}^{\infty} \ker(A - \lambda I)^k$.

!!! theorem "Theorem 12.3 (Dimension of the generalized eigenspace)"
    Let $\lambda$ be an eigenvalue of $A$ with algebraic multiplicity $m$. Then

    $$
    \dim G_\lambda = m.
    $$

    That is, the dimension of the generalized eigenspace equals the algebraic multiplicity of the eigenvalue.

??? proof "Proof"
    By the Jordan normal form theorem (Theorem 12.8), $A$ is similar to a Jordan matrix $J$. Let $A = PJP^{-1}$, then $(A - \lambda I)^n = P(J - \lambda I)^n P^{-1}$. For the Jordan matrix $J$, $\ker(J - \lambda I)^n$ is spanned by all columns corresponding to Jordan blocks associated with $\lambda$, and its dimension equals the algebraic multiplicity $m$ of $\lambda$. Since similarity transformations preserve dimension, $\dim G_\lambda = m$. $\blacksquare$

!!! theorem "Theorem 12.4 (Direct sum decomposition of generalized eigenspaces)"
    Let $A$ be an $n \times n$ matrix (over the complex numbers) with distinct eigenvalues $\lambda_1, \ldots, \lambda_s$. Then

    $$
    \mathbb{C}^n = G_{\lambda_1} \oplus G_{\lambda_2} \oplus \cdots \oplus G_{\lambda_s}.
    $$

??? proof "Proof"
    **Dimension verification:** By Theorem 12.3, $\sum \dim G_{\lambda_i} = \sum m_i = n$ (the sum of algebraic multiplicities equals $n$).

    **Directness:** We need to show that if $\mathbf{v}_1 + \cdots + \mathbf{v}_s = \mathbf{0}$ ($\mathbf{v}_i \in G_{\lambda_i}$), then $\mathbf{v}_i = \mathbf{0}$. Suppose $(A - \lambda_i I)^{m_i} \mathbf{v}_i = \mathbf{0}$. For $i \neq j$, let $p_j(x) = \prod_{i \neq j}(x - \lambda_i)^{m_i}$. Then $p_j(A)\mathbf{v}_i = \mathbf{0}$ ($i \neq j$), while $p_j(A)\mathbf{v}_j = c_j \mathbf{v}_j$ (where $c_j = \prod_{i \neq j}(\lambda_j - \lambda_i)^{m_i} \neq 0$).

    Applying $p_j(A)$ to both sides of $\mathbf{v}_1 + \cdots + \mathbf{v}_s = \mathbf{0}$: $c_j \mathbf{v}_j = \mathbf{0}$, so $\mathbf{v}_j = \mathbf{0}$. $\blacksquare$

!!! example "Example 12.2"
    Let $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$. Find the generalized eigenspace of $A$.

    **Solution:** The eigenvalue of $A$ is $\lambda = 2$ (algebraic multiplicity 2).

    $(A - 2I) = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$, $\ker(A - 2I) = \operatorname{span}\left\{\begin{pmatrix}1\\0\end{pmatrix}\right\}$, dimension 1 (geometric multiplicity).

    $(A - 2I)^2 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$, $\ker(A - 2I)^2 = \mathbb{R}^2$, dimension 2.

    Therefore $G_2 = \mathbb{R}^2$. The generalized eigenvector $\mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}$ satisfies $(A-2I)\mathbf{v}_2 = \begin{pmatrix}1\\0\end{pmatrix} \neq \mathbf{0}$, but $(A-2I)^2\mathbf{v}_2 = \mathbf{0}$, so $\mathbf{v}_2$ is a generalized eigenvector of rank 2.

---

## 12.3 Nilpotent Matrices

<div class="context-flow" markdown>

$N^k = 0$ → All eigenvalues are $0$ → $(I-N)^{-1} = I + N + \cdots + N^{k-1}$ (finite Neumann series) → Jordan block = $\lambda I + N$

</div>

Nilpotent matrices are the key to understanding the Jordan normal form. Every Jordan block can be decomposed as a scalar matrix plus a nilpotent matrix.

!!! definition "Definition 12.5 (Nilpotent Matrix)"
    If there exists a positive integer $k$ such that $N^k = 0$, then matrix $N$ is called a **nilpotent matrix**. The smallest positive integer $k$ satisfying $N^k = 0$ is called the **nilpotency index** of $N$.

!!! theorem "Theorem 12.5 (Properties of nilpotent matrices)"
    Let $N$ be an $n \times n$ nilpotent matrix with nilpotency index $k$. Then:

    1. All eigenvalues of $N$ are $0$;
    2. $\operatorname{tr}(N) = 0$, $\det(N) = 0$;
    3. $k \le n$;
    4. The characteristic polynomial of $N$ is $p(\lambda) = \lambda^n$;
    5. $I - N$ is invertible, and $(I - N)^{-1} = I + N + N^2 + \cdots + N^{k-1}$.

??? proof "Proof"
    **(1)** Let $\lambda$ be an eigenvalue of $N$ with eigenvector $\mathbf{v}$. Then $N\mathbf{v} = \lambda\mathbf{v}$ and $N^k\mathbf{v} = \lambda^k\mathbf{v}$. From $N^k = 0$, we get $\lambda^k\mathbf{v} = \mathbf{0}$. Since $\mathbf{v} \neq \mathbf{0}$, $\lambda^k = 0$, i.e., $\lambda = 0$.

    **(2)** By (1), all eigenvalues are zero, so $\operatorname{tr}(N) = 0$ (sum of eigenvalues) and $\det(N) = 0$ (product of eigenvalues).

    **(3)** Consider the chain of subspaces $\{0\} \subseteq \ker(N) \subseteq \ker(N^2) \subseteq \cdots$. If $\ker(N^j) = \ker(N^{j+1})$, then $\ker(N^i) = \ker(N^j)$ for all $i \ge j$. Since $\ker(N^k) = \mathbb{R}^n$ (because $N^k = 0$) and the dimensions strictly increase until stabilization, $k \le n$.

    **(4)** By (1), $N$ has only zero eigenvalues, so $p(\lambda) = \lambda^n$.

    **(5)** $(I - N)(I + N + \cdots + N^{k-1}) = I - N^k = I$. $\blacksquare$

!!! example "Example 12.3"
    Let $N = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$. Verify that $N$ is nilpotent and find its nilpotency index.

    **Solution:**

    $$
    N^2 = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}, \qquad N^3 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.
    $$

    $N^2 \neq 0$ but $N^3 = 0$, so the nilpotency index is $3$.

!!! theorem "Theorem 12.6 (Jordan normal form of a nilpotent matrix)"
    Let $N$ be an $n \times n$ nilpotent matrix with nilpotency index $k$. Then $N$ is similar to a block diagonal matrix

    $$
    J = \operatorname{diag}(J_{n_1}(0), J_{n_2}(0), \ldots, J_{n_t}(0)),
    $$

    where $J_{n_i}(0)$ is an $n_i \times n_i$ nilpotent Jordan block (zeros on the diagonal, ones on the superdiagonal), $n_1 \ge n_2 \ge \cdots \ge n_t \ge 1$, $\sum n_i = n$, $n_1 = k$.

??? proof "Proof"
    This is a special case of the Jordan normal form theorem (Theorem 12.8) with $\lambda = 0$. The nilpotency index $k$ equals the size of the largest Jordan block $n_1$, since $J_{n_1}(0)^{n_1} = 0$ but $J_{n_1}(0)^{n_1-1} \neq 0$. See Theorem 12.8 for the detailed proof.

---

## 12.4 Jordan Blocks

<div class="context-flow" markdown>

$J_k(\lambda) = \lambda I + N_k$: unique eigenvalue $\lambda$, algebraic multiplicity $k$, geometric multiplicity **1** → Binomial theorem $J_k(\lambda)^m = \sum \binom{m}{j}\lambda^{m-j}N_k^j$

</div>

Jordan blocks are the basic building units of the Jordan normal form.

!!! definition "Definition 12.6 (Jordan Block)"
    The $k \times k$ **Jordan block** $J_k(\lambda)$ is defined as

    $$
    J_k(\lambda) = \begin{pmatrix}
    \lambda & 1 & 0 & \cdots & 0 \\
    0 & \lambda & 1 & \cdots & 0 \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & \lambda & 1 \\
    0 & \cdots & 0 & 0 & \lambda
    \end{pmatrix} = \lambda I_k + N_k,
    $$

    where $N_k$ is the $k \times k$ **basic nilpotent matrix** (all superdiagonal entries equal to 1, the rest zero).

!!! theorem "Theorem 12.7 (Properties of Jordan blocks)"
    $J_k(\lambda)$ has the following properties:

    1. The unique eigenvalue is $\lambda$, with algebraic multiplicity $k$ and geometric multiplicity $1$;
    2. $(J_k(\lambda) - \lambda I)^k = N_k^k = 0$, but $(J_k(\lambda) - \lambda I)^{k-1} \neq 0$;
    3. $J_k(\lambda)^m = \sum_{j=0}^{k-1} \binom{m}{j} \lambda^{m-j} N_k^j$ ($m \ge k-1$);
    4. When $\lambda \neq 0$, $J_k(\lambda)$ is invertible, and $J_k(\lambda)^{-1} = \frac{1}{\lambda}(I + \frac{1}{\lambda}N_k)^{-1}$.

??? proof "Proof"
    **(1)** $\det(J_k(\lambda) - \mu I) = (\lambda - \mu)^k$, so the unique eigenvalue is $\lambda$ with algebraic multiplicity $k$. $(J_k(\lambda) - \lambda I)\mathbf{v} = N_k\mathbf{v} = \mathbf{0}$ has solution space $\operatorname{span}\{\mathbf{e}_1\}$, so the geometric multiplicity is 1.

    **(2)** $J_k(\lambda) - \lambda I = N_k$, and $N_k^k = 0$, $N_k^{k-1} \neq 0$ (the $(1,k)$ entry of $N_k^{k-1}$ is 1).

    **(3)** By the binomial theorem, $J_k(\lambda)^m = (\lambda I + N_k)^m = \sum_{j=0}^{m} \binom{m}{j} \lambda^{m-j} N_k^j$. Since $N_k^j = 0$ for $j \ge k$, the sum effectively runs only up to $k-1$.

    **(4)** $J_k(\lambda) = \lambda(I + \frac{1}{\lambda}N_k)$. By Theorem 12.5(5), $I + \frac{1}{\lambda}N_k$ is invertible. $\blacksquare$

!!! example "Example 12.4"
    Compute $J_3(2)^4$.

    **Solution:** $J_3(2) = 2I + N_3$, where $N_3 = \begin{pmatrix}0&1&0\\0&0&1\\0&0&0\end{pmatrix}$.

    By the formula:

    $$
    J_3(2)^4 = \sum_{j=0}^{2} \binom{4}{j} 2^{4-j} N_3^j = \binom{4}{0}2^4 I + \binom{4}{1}2^3 N_3 + \binom{4}{2}2^2 N_3^2
    $$

    $$
    = 16I + 32N_3 + 24N_3^2 = \begin{pmatrix}16&32&24\\0&16&32\\0&0&16\end{pmatrix}.
    $$

---

## 12.5 The Jordan Normal Form Theorem

<div class="context-flow" markdown>

**Core of the chapter**: Any complex matrix $A = PJP^{-1}$, where $J$ consists of Jordan blocks and is **unique up to permutation** → The block structure is determined by the sequence $\operatorname{rank}(A-\lambda I)^j$

</div>

!!! definition "Definition 12.7 (Jordan Matrix)"
    A **Jordan matrix** is a block diagonal matrix

    $$
    J = \operatorname{diag}(J_{k_1}(\lambda_1), J_{k_2}(\lambda_2), \ldots, J_{k_s}(\lambda_s)),
    $$

    where each $J_{k_i}(\lambda_i)$ is a Jordan block. Different blocks may have the same $\lambda$ value.

<div class="context-flow" markdown>

**Insight**: Three steps for existence — generalized eigenspace decomposition → reduce to a nilpotent problem on each $G_\lambda$ → construct **Jordan chains** to solve block by block

</div>

!!! theorem "Theorem 12.8 (Jordan Normal Form Theorem)"
    Let $A$ be an $n \times n$ complex matrix. Then there exists an invertible matrix $P$ such that

    $$
    P^{-1}AP = J = \operatorname{diag}(J_{k_1}(\lambda_1), \ldots, J_{k_s}(\lambda_s)),
    $$

    where $J$ is a Jordan matrix.

    Moreover, $J$ is unique **up to permutation of the Jordan blocks**. That is, the sizes of the Jordan blocks and their corresponding eigenvalues are uniquely determined by $A$.

??? proof "Proof"
    **Proof outline:** (The complete proof is lengthy; we give the key steps here.)

    **Existence:**

    Step 1 (Generalized eigenspace decomposition): By Theorem 12.4, $\mathbb{C}^n = G_{\lambda_1} \oplus \cdots \oplus G_{\lambda_s}$. Each $G_{\lambda_i}$ is an invariant subspace of $A$, and $A|_{G_{\lambda_i}}$ denotes the restriction of $A$ to $G_{\lambda_i}$.

    Step 2 (Reduction to the nilpotent case): On $G_{\lambda_i}$, $A|_{G_{\lambda_i}} - \lambda_i I$ is nilpotent (since $(A - \lambda_i I)^{m_i}|_{G_{\lambda_i}} = 0$, where $m_i$ is the algebraic multiplicity of $\lambda_i$).

    Step 3 (Jordan form of nilpotent matrices): Let $N$ be a $d \times d$ nilpotent matrix. By constructing Jordan chains, one can find a basis under which $N$ has Jordan form.

    Specifically, choose $\mathbf{v}$ such that $(A-\lambda I)^{k-1}\mathbf{v} \neq \mathbf{0}$ but $(A-\lambda I)^k\mathbf{v} = \mathbf{0}$. Define the Jordan chain:

    $$
    \mathbf{v}, (A-\lambda I)\mathbf{v}, (A-\lambda I)^2\mathbf{v}, \ldots, (A-\lambda I)^{k-1}\mathbf{v}.
    $$

    These $k$ vectors are linearly independent, and $A$ restricted to their span is represented by $J_k(\lambda)$.

    Repeating this process decomposes the entire generalized eigenspace into a direct sum of Jordan chains.

    **Uniqueness:**

    The Jordan block structure is determined by the following invariants: for each eigenvalue $\lambda$ and each positive integer $j$,

    $$
    r_j = \operatorname{rank}(A - \lambda I)^j
    $$

    is a similarity invariant. The number of Jordan blocks of size $k$ associated with $J_k(\lambda)$ is

    $$
    n_k = r_{k-2} - 2r_{k-1} + r_k
    $$

    (with the convention $r_{-1} = n$, $r_0 = n - \dim\ker(A-\lambda I)$), which is uniquely determined by $A$. $\blacksquare$

!!! theorem "Theorem 12.9 (Jordan normal form and matrix properties)"
    Let the Jordan normal form of $A$ be $J = \operatorname{diag}(J_{k_1}(\lambda_1), \ldots, J_{k_s}(\lambda_s))$. Then:

    1. $A$ is diagonalizable if and only if all Jordan blocks have size $1$ (i.e., $k_i = 1$);
    2. The geometric multiplicity of eigenvalue $\lambda$ equals the number of Jordan blocks corresponding to $\lambda$;
    3. The algebraic multiplicity of eigenvalue $\lambda$ equals the sum of the sizes of all Jordan blocks corresponding to $\lambda$;
    4. $\det(A) = \prod_{i=1}^s \lambda_i^{k_i}$, $\operatorname{tr}(A) = \sum_{i=1}^s k_i \lambda_i$.

??? proof "Proof"
    **(1)** $A$ is diagonalizable $\Leftrightarrow$ for each eigenvalue, geometric multiplicity = algebraic multiplicity $\Leftrightarrow$ every Jordan block is $1 \times 1$.

    **(2)** The Jordan blocks corresponding to eigenvalue $\lambda$ are $J_{k_1}(\lambda), \ldots, J_{k_t}(\lambda)$. Each $J_{k_i}(\lambda)$ has eigenspace dimension 1, so the total geometric multiplicity is $t$ (the number of blocks).

    **(3)** Each Jordan block $J_{k_i}(\lambda)$ contributes algebraic multiplicity $k_i$, so the total algebraic multiplicity is $\sum k_i$.

    **(4)** Similarity transformations preserve the determinant and trace. $\det(J) = \prod \det(J_{k_i}(\lambda_i)) = \prod \lambda_i^{k_i}$. $\operatorname{tr}(J) = \sum k_i\lambda_i$. $\blacksquare$

!!! example "Example 12.5"
    A $5 \times 5$ matrix $A$ has eigenvalues $\lambda_1 = 2$ (algebraic multiplicity 3, geometric multiplicity 2) and $\lambda_2 = -1$ (algebraic multiplicity 2, geometric multiplicity 1). Determine the Jordan normal form of $A$.

    **Solution:**

    For $\lambda_1 = 2$: algebraic multiplicity 3, geometric multiplicity 2, so there are 2 Jordan blocks with sizes summing to 3. The possibilities are $(2,1)$ (one $2\times 2$ block and one $1\times 1$ block).

    For $\lambda_2 = -1$: algebraic multiplicity 2, geometric multiplicity 1, so there is 1 Jordan block of size 2.

    The Jordan normal form is

    $$
    J = \begin{pmatrix}
    2 & 1 & & & \\
    0 & 2 & & & \\
    & & 2 & & \\
    & & & -1 & 1 \\
    & & & 0 & -1
    \end{pmatrix}.
    $$

---

## 12.6 Minimal Polynomial

<div class="context-flow" markdown>

$m_A(\lambda) = \prod(\lambda - \lambda_i)^{d_i}$, $d_i$ = largest Jordan block size → **Diagonalizable** ↔ $m_A$ has no repeated roots → Directly connects to Ch13 matrix function definition conditions

</div>

The minimal polynomial has a close connection with the Jordan normal form.

!!! definition "Definition 12.8 (Minimal Polynomial)"
    The **minimal polynomial** $m_A(\lambda)$ of matrix $A$ is the monic polynomial of lowest degree satisfying $m_A(A) = 0$.

!!! theorem "Theorem 12.10 (Cayley-Hamilton theorem and minimal polynomial)"
    Let $p_A(\lambda)$ be the characteristic polynomial and $m_A(\lambda)$ the minimal polynomial of $A$. Then:

    1. $m_A(\lambda)$ divides $p_A(\lambda)$;
    2. $m_A(\lambda)$ and $p_A(\lambda)$ have the same roots (i.e., the same set of eigenvalues);
    3. $A$ is diagonalizable if and only if $m_A(\lambda)$ has no repeated roots.

??? proof "Proof"
    **(1)** By the Cayley-Hamilton theorem, $p_A(A) = 0$. By the definition of the minimal polynomial, $m_A$ is the lowest-degree monic polynomial with $m_A(A) = 0$. Using polynomial division $p_A = q \cdot m_A + r$, we get $r(A) = p_A(A) - q(A)m_A(A) = 0$. If $r \neq 0$, then $\deg r < \deg m_A$, a contradiction. So $r = 0$ and $m_A | p_A$.

    **(2)** Let $\lambda_0$ be a root of $p_A$ with eigenvector $\mathbf{v}$. Then $m_A(A)\mathbf{v} = m_A(\lambda_0)\mathbf{v} = \mathbf{0}$. Since $\mathbf{v} \neq \mathbf{0}$, $m_A(\lambda_0) = 0$. The reverse direction follows from (1).

    **(3)** $A$ is diagonalizable $\Leftrightarrow$ all Jordan blocks are $1 \times 1$ $\Leftrightarrow$ $m_A(\lambda) = \prod(\lambda - \lambda_i)$ (no repeated roots).

    In detail, $m_A(\lambda) = \prod_{i=1}^s (\lambda - \lambda_i)^{d_i}$, where $d_i$ is the size of the largest Jordan block for $\lambda_i$. $A$ is diagonalizable $\Leftrightarrow$ all $d_i = 1$ $\Leftrightarrow$ $m_A$ has no repeated roots. $\blacksquare$

!!! theorem "Theorem 12.11 (Relationship between minimal polynomial and Jordan form)"
    Let $d_i$ be the size of the largest Jordan block corresponding to eigenvalue $\lambda_i$ in the Jordan normal form of $A$. Then

    $$
    m_A(\lambda) = \prod_{i=1}^{s} (\lambda - \lambda_i)^{d_i}.
    $$

??? proof "Proof"
    Let $q(\lambda) = \prod_{i=1}^s (\lambda - \lambda_i)^{d_i}$. We need to show $q(A) = 0$ and that $q$ is the lowest-degree monic polynomial with this property.

    $q(J) = \operatorname{diag}(q(J_{k_1}(\lambda_1)), \ldots)$. For each Jordan block $J_k(\lambda_i)$, $q(J_k(\lambda_i))$ contains the factor $(J_k(\lambda_i) - \lambda_i I)^{d_i} = N_k^{d_i}$. Since $k \le d_i$ ($d_i$ is the largest block size), $N_k^{d_i} = 0$, so $q(J_k(\lambda_i)) = 0$. Therefore $q(J) = 0$ and $q(A) = Pq(J)P^{-1} = 0$.

    If the degree of $q$ could be lower, then for some $\lambda_i$ the power of $(\lambda-\lambda_i)$ would be less than $d_i$, but this would mean the largest Jordan block $J_{d_i}(\lambda_i)$ is not annihilated, a contradiction. $\blacksquare$

!!! example "Example 12.6"
    Let $A = \begin{pmatrix}3&1&0\\0&3&0\\0&0&5\end{pmatrix}$. Find its minimal polynomial.

    **Solution:** $A$ is already in Jordan form: $J_2(3) \oplus J_1(5)$.

    Eigenvalue $\lambda_1 = 3$ (largest block size 2), $\lambda_2 = 5$ (largest block size 1).

    Minimal polynomial: $m_A(\lambda) = (\lambda - 3)^2(\lambda - 5)$.

    Verification: $p_A(\lambda) = (\lambda-3)^2(\lambda-5)$. Here the minimal polynomial equals the characteristic polynomial.

---

## 12.7 Computing the Jordan Normal Form

<div class="context-flow" markdown>

Algorithm: Characteristic polynomial → Successively compute $\dim\ker(A-\lambda I)^j$ to determine block structure → Construct Jordan chains to obtain the transition matrix $P$

</div>

This section demonstrates how to compute the Jordan normal form and the transition matrix through a systematic method and detailed examples.

!!! definition "Definition 12.9 (Jordan Chain)"
    Let $\lambda$ be an eigenvalue of $A$. If the vector sequence $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k$ satisfies

    $$
    (A - \lambda I)\mathbf{v}_1 = \mathbf{0}, \quad (A - \lambda I)\mathbf{v}_j = \mathbf{v}_{j-1}, \quad j = 2, \ldots, k,
    $$

    then $\{\mathbf{v}_k, \mathbf{v}_{k-1}, \ldots, \mathbf{v}_1\}$ is called a **Jordan chain** of length $k$. Here $\mathbf{v}_1$ is an eigenvector, and $\mathbf{v}_2, \ldots, \mathbf{v}_k$ are generalized eigenvectors.

!!! definition "Definition 12.10 (Steps for computing the Jordan normal form)"
    Given an $n \times n$ matrix $A$, the systematic method for computing its Jordan normal form is:

    1. Find the characteristic polynomial and determine the eigenvalues $\lambda_i$ and their algebraic multiplicities $m_i$;
    2. For each $\lambda_i$, successively compute the dimensions of $\ker(A-\lambda_i I)^j$ ($j = 1, 2, \ldots$) to determine the Jordan block structure;
    3. The number of Jordan blocks of size $k$ = $\dim\ker(A-\lambda I)^k - 2\dim\ker(A-\lambda I)^{k-1} + \dim\ker(A-\lambda I)^{k-2}$;
    4. Construct Jordan chains to determine the transition matrix $P$.

!!! example "Example 12.7"
    Find the Jordan normal form of the matrix $A = \begin{pmatrix} 5 & 4 & 2 & 1 \\ 0 & 1 & -1 & -1 \\ -1 & -1 & 3 & 0 \\ 1 & 1 & -1 & 2 \end{pmatrix}$.

    **Solution:**

    **Step 1:** The characteristic polynomial (computed via determinant expansion or other methods):

    $$
    p_A(\lambda) = (\lambda - 1)^2(\lambda - 4)(\lambda - 3).
    $$

    Eigenvalues: $\lambda_1 = 1$ (algebraic multiplicity 2), $\lambda_2 = 4$ (algebraic multiplicity 1), $\lambda_3 = 3$ (algebraic multiplicity 1).

    **Step 2:** For $\lambda_1 = 1$:

    $$
    A - I = \begin{pmatrix} 4&4&2&1 \\ 0&0&-1&-1 \\ -1&-1&2&0 \\ 1&1&-1&1 \end{pmatrix}.
    $$

    Computing $\operatorname{rank}(A-I)$: by row reduction, $\operatorname{rank}(A-I) = 3$, so $\dim\ker(A-I) = 1$.

    The geometric multiplicity is 1 < algebraic multiplicity 2, so there is one $2 \times 2$ Jordan block $J_2(1)$.

    For $\lambda_2 = 4$ and $\lambda_3 = 3$: both have algebraic multiplicity 1, so each has a single $1 \times 1$ block.

    **Step 3:** The Jordan normal form is

    $$
    J = \begin{pmatrix}
    1 & 1 & & \\
    0 & 1 & & \\
    & & 4 & \\
    & & & 3
    \end{pmatrix}.
    $$

!!! example "Example 12.8"
    Find the Jordan normal form and the transition matrix of $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$.

    **Solution:** $A$ is already in Jordan form $J_2(2)$.

    Eigenvalue $\lambda = 2$ (algebraic multiplicity 2, geometric multiplicity 1).

    $(A - 2I) = \begin{pmatrix}0&1\\0&0\end{pmatrix}$, $\ker(A-2I) = \operatorname{span}\left\{\begin{pmatrix}1\\0\end{pmatrix}\right\}$.

    Eigenvector $\mathbf{v}_1 = \begin{pmatrix}1\\0\end{pmatrix}$. Generalized eigenvector $\mathbf{v}_2$ satisfying $(A-2I)\mathbf{v}_2 = \mathbf{v}_1$:

    $$
    \begin{pmatrix}0&1\\0&0\end{pmatrix}\mathbf{v}_2 = \begin{pmatrix}1\\0\end{pmatrix}, \quad \Rightarrow \quad \mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}.
    $$

    Transition matrix $P = [\mathbf{v}_1, \mathbf{v}_2] = \begin{pmatrix}1&0\\0&1\end{pmatrix} = I$.

    Verification: $P^{-1}AP = A = J_2(2)$.

!!! example "Example 12.9"
    Find the Jordan normal form of the matrix $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 8 & -12 & 6 \end{pmatrix}$.

    **Solution:** The characteristic polynomial is

    $$
    p_A(\lambda) = -\lambda^3 + 6\lambda^2 - 12\lambda + 8 = -(\lambda - 2)^3.
    $$

    The unique eigenvalue is $\lambda = 2$ with algebraic multiplicity 3.

    $$
    A - 2I = \begin{pmatrix}-2&1&0\\0&-2&1\\8&-12&4\end{pmatrix}.
    $$

    Row reduction gives $\operatorname{rank}(A - 2I) = 2$, $\dim\ker(A-2I) = 1$, so the geometric multiplicity is 1.

    Rank of $(A-2I)^2$: computing

    $$
    (A-2I)^2 = \begin{pmatrix}4&-4&1\\-8&8&-2\\-16&16&-4\end{pmatrix} + \cdots
    $$

    Actual computation gives $\operatorname{rank}(A-2I)^2 = 1$, $\dim\ker(A-2I)^2 = 2$.

    $(A-2I)^3 = 0$, $\dim\ker(A-2I)^3 = 3$.

    The number of blocks of size 3 = $3 - 2 \times 2 + 1 = 0$...

    Re-applying the formula: the Jordan block count is determined by the increasing sequence $d_j = \dim\ker(A-\lambda I)^j$. $d_0 = 0, d_1 = 1, d_2 = 2, d_3 = 3$. The increments are $\Delta_1 = 1, \Delta_2 = 1, \Delta_3 = 1$.

    The number of blocks of size $\ge k$ is $\Delta_k$. So the number of blocks of size $\ge 1$ is 1, size $\ge 2$ is 1, size $\ge 3$ is 1. Therefore there is exactly 1 Jordan block of size 3.

    Jordan normal form: $J = J_3(2) = \begin{pmatrix}2&1&0\\0&2&1\\0&0&2\end{pmatrix}$.

---

## 12.8 Applications of the Jordan Form

<div class="context-flow" markdown>

$A^n = PJ^nP^{-1}$ (matrix powers) · $e^{At} = Pe^{Jt}P^{-1}$ (differential equations) → The Jordan form reduces all matrix functions to operations on Jordan blocks → Foundation for Ch13

</div>

### 12.8.1 Computing Matrix Powers

!!! theorem "Theorem 12.12 (Computing matrix powers via the Jordan form)"
    Let $A = PJP^{-1}$. Then $A^n = PJ^nP^{-1}$, where

    $$
    J^n = \operatorname{diag}(J_{k_1}(\lambda_1)^n, \ldots, J_{k_s}(\lambda_s)^n),
    $$

    and

    $$
    J_k(\lambda)^n = \begin{pmatrix}
    \lambda^n & \binom{n}{1}\lambda^{n-1} & \binom{n}{2}\lambda^{n-2} & \cdots & \binom{n}{k-1}\lambda^{n-k+1} \\
    0 & \lambda^n & \binom{n}{1}\lambda^{n-1} & \cdots & \binom{n}{k-2}\lambda^{n-k+2} \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & \lambda^n & \binom{n}{1}\lambda^{n-1} \\
    0 & \cdots & 0 & 0 & \lambda^n
    \end{pmatrix}.
    $$

??? proof "Proof"
    $J_k(\lambda)^n = (\lambda I + N_k)^n = \sum_{j=0}^{k-1}\binom{n}{j}\lambda^{n-j}N_k^j$. The $(p,q)$ entry of $N_k^j$ is $\delta_{p+j,q}$ (the $j$-th superdiagonal is all 1's). Therefore the $(p,q)$ entry of $J_k(\lambda)^n$ is:

    $$
    [J_k(\lambda)^n]_{pq} = \begin{cases} \binom{n}{q-p}\lambda^{n-q+p} & \text{if } q \ge p, \\ 0 & \text{if } q < p. \end{cases}
    $$

    Here the convention is $\binom{n}{j} = 0$ when $j > n$ or $j < 0$. $\blacksquare$

!!! example "Example 12.10"
    Let $A = \begin{pmatrix}3&1\\0&3\end{pmatrix}$. Compute $A^{100}$.

    **Solution:** $A = J_2(3)$, so

    $$
    A^{100} = J_2(3)^{100} = \begin{pmatrix} 3^{100} & 100 \cdot 3^{99} \\ 0 & 3^{100} \end{pmatrix}.
    $$

### 12.8.2 Systems of Differential Equations

!!! theorem "Theorem 12.13 (Jordan form and linear constant-coefficient ODE systems)"
    The solution of the linear constant-coefficient ODE system $\mathbf{x}'(t) = A\mathbf{x}(t)$ is $\mathbf{x}(t) = e^{At}\mathbf{x}(0)$. Using the Jordan decomposition $A = PJP^{-1}$:

    $$
    e^{At} = Pe^{Jt}P^{-1} = P \operatorname{diag}(e^{J_{k_1}(\lambda_1)t}, \ldots, e^{J_{k_s}(\lambda_s)t}) P^{-1},
    $$

    where

    $$
    e^{J_k(\lambda)t} = e^{\lambda t}\begin{pmatrix}
    1 & t & \frac{t^2}{2!} & \cdots & \frac{t^{k-1}}{(k-1)!} \\
    0 & 1 & t & \cdots & \frac{t^{k-2}}{(k-2)!} \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & 1 & t \\
    0 & \cdots & 0 & 0 & 1
    \end{pmatrix}.
    $$

??? proof "Proof"
    $e^{J_k(\lambda)t} = e^{(\lambda I + N_k)t} = e^{\lambda t I} e^{N_k t}$ (since $\lambda I$ and $N_k$ commute).

    $e^{\lambda t I} = e^{\lambda t} I$, $e^{N_k t} = \sum_{j=0}^{k-1} \frac{t^j}{j!} N_k^j$ (since $N_k^k = 0$).

    Therefore $e^{J_k(\lambda)t} = e^{\lambda t} \sum_{j=0}^{k-1} \frac{t^j}{j!} N_k^j$, with matrix entries as described above. $\blacksquare$

!!! example "Example 12.11"
    Solve the ODE system $\mathbf{x}' = \begin{pmatrix}2&1\\0&2\end{pmatrix}\mathbf{x}$ with initial condition $\mathbf{x}(0) = \begin{pmatrix}1\\3\end{pmatrix}$.

    **Solution:** $A = J_2(2)$, so

    $$
    e^{At} = e^{2t}\begin{pmatrix}1&t\\0&1\end{pmatrix}.
    $$

    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}(0) = e^{2t}\begin{pmatrix}1&t\\0&1\end{pmatrix}\begin{pmatrix}1\\3\end{pmatrix} = e^{2t}\begin{pmatrix}1+3t\\3\end{pmatrix}.
    $$

    That is, $x_1(t) = (1+3t)e^{2t}$, $x_2(t) = 3e^{2t}$.

---

## Chapter Summary

This chapter systematically introduced the Jordan normal form theory:

1. **Invariant subspaces** provide the framework for studying linear transformations block by block;
2. **Generalized eigenvectors** compensate for the deficiency of eigenvectors; $\dim G_\lambda$ equals the algebraic multiplicity;
3. **Jordan blocks** $J_k(\lambda) = \lambda I + N_k$ are the "nearly diagonal" building units;
4. The **Jordan normal form theorem** guarantees that every square matrix (over the complex numbers) is similar to a unique Jordan matrix;
5. The **minimal polynomial** reflects the size of the largest Jordan block; $A$ is diagonalizable if and only if the minimal polynomial has no repeated roots;
6. **Applications** include efficient computation of matrix powers $A^n$ and solving linear constant-coefficient ODE systems.

The Jordan normal form is one of the most profound results in matrix theory, laying the foundation for the matrix function theory that follows.
