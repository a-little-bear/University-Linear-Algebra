# Chapter 10  Matrix Decompositions

<div class="context-flow" markdown>

**Prerequisites**: Ch8 Orthogonality/Spectral theorem · Ch9 Positive definiteness · **Chapter arc**: LU (elimination) → PLU (pivoting) → Cholesky (positive definite) → QR (orthogonalization) → Schur (triangularization) → Spectral decomposition → Polar decomposition
Essence: Factoring a matrix into a product of "structurally simple factors" — each decomposition answers a different question

</div>

Matrix decomposition (matrix factorization) is the technique of expressing a matrix as a product of several special matrices, and it is a core tool in numerical linear algebra and applied mathematics. Different matrix decompositions are suited to different applications: LU decomposition for efficiently solving linear systems, QR decomposition for least squares problems, Schur decomposition for revealing deep matrix structure, and spectral decomposition for characterizing the essence of symmetric and normal matrices. This chapter systematically introduces various important matrix decomposition methods.

---

## 10.1 LU Decomposition

<div class="context-flow" markdown>

Matrix formulation of Gaussian elimination: $A = LU$ (lower triangular $\times$ upper triangular) → Perform the $O(n^3)$ factorization once, then each right-hand side $\mathbf{b}$ requires only $O(n^2)$

</div>

!!! definition "Definition 10.1 (LU decomposition)"
    Let $A$ be an $n \times n$ matrix. If there exist an $n \times n$ lower triangular matrix $L$ (with all diagonal entries equal to 1) and an upper triangular matrix $U$ such that

    $$A = LU$$

    then $A = LU$ is called the **LU decomposition** (LU factorization) of $A$. $L$ is called a **unit lower triangular matrix** and $U$ is called an **upper triangular matrix**.

!!! theorem "Theorem 10.1 (Existence and uniqueness of LU decomposition)"
    An $n \times n$ matrix $A$ has an LU decomposition if and only if all leading principal submatrices $A_k$ ($k = 1, \ldots, n-1$) are nonsingular, i.e., $\det(A_k) \neq 0$. If $A$ is invertible, then the LU decomposition is unique.

??? proof "Proof"
    **Necessity:** If $A = LU$ with $L, U$ lower and upper triangular respectively, then $A_k = L_k U_k$, where $L_k, U_k$ are the leading $k \times k$ principal submatrices of $L, U$. Since $L_k$ is unit lower triangular, $\det(L_k) = 1 \neq 0$. If $\det(A_k) = 0$, then $\det(U_k) = 0$, meaning one of the first $k$ diagonal entries of $U$ is zero. But Gaussian elimination requires nonzero pivots, a contradiction.

    **Sufficiency:** In the language of Gaussian elimination, $\det(A_k) \neq 0$ guarantees that all pivots are nonzero and no row swaps are needed. The elimination process is recorded as $A = L_1^{-1}L_2^{-1}\cdots L_{n-1}^{-1}U$, where each $L_i$ is an elementary lower triangular matrix. $L = L_1^{-1}L_2^{-1}\cdots L_{n-1}^{-1}$ is still a unit lower triangular matrix.

    **Uniqueness (when $A$ is invertible):** If $A = L_1U_1 = L_2U_2$, then $L_2^{-1}L_1 = U_2U_1^{-1}$. The left side is unit lower triangular and the right side is upper triangular, so both equal the identity matrix $I$. $\blacksquare$

!!! example "Example 10.1"
    Perform LU decomposition of the matrix $A = \begin{pmatrix} 2 & 1 & 1 \\ 4 & 3 & 3 \\ 8 & 7 & 9 \end{pmatrix}$.

    Gaussian elimination: $R_2 \leftarrow R_2 - 2R_1$, $R_3 \leftarrow R_3 - 4R_1$:

    $$\begin{pmatrix} 2 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 3 & 5 \end{pmatrix}$$

    $R_3 \leftarrow R_3 - 3R_2$:

    $$U = \begin{pmatrix} 2 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix}$$

    Multiplier matrix $L = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 4 & 3 & 1 \end{pmatrix}$.

    Verification: $LU = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 4 & 3 & 1 \end{pmatrix}\begin{pmatrix} 2 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix} = \begin{pmatrix} 2 & 1 & 1 \\ 4 & 3 & 3 \\ 8 & 7 & 9 \end{pmatrix} = A$.

!!! note "Note"
    The main application of LU decomposition is efficiently solving linear systems $A\mathbf{x} = \mathbf{b}$: first solve $L\mathbf{y} = \mathbf{b}$ (forward substitution), then solve $U\mathbf{x} = \mathbf{y}$ (back substitution). The cost is $O(n^2)$ (assuming $L, U$ are known), while the factorization itself costs $O(\frac{2}{3}n^3)$.

---

## 10.2 PLU Decomposition

<div class="context-flow" markdown>

LU requires nonzero leading principal minors → Introduce **permutation matrices** $P$ (row swaps) → $PA = LU$ holds for **any invertible matrix**

</div>

!!! definition "Definition 10.2 (Permutation matrix)"
    A **permutation matrix** is a square matrix that has exactly one $1$ in each row and each column, with all other entries $0$. It corresponds to a permutation of the rows (columns) of the identity matrix.

!!! proposition "Proposition 10.1 (Properties of permutation matrices)"
    Let $P$ be a permutation matrix. Then

    1. $P$ is an orthogonal matrix: $P^TP = PP^T = I$, i.e., $P^{-1} = P^T$;
    2. $\det(P) = \pm 1$;
    3. The product of permutation matrices is again a permutation matrix.

??? proof "Proof"
    (1) The rows (columns) of $P$ are a permutation of the standard basis vectors, hence an orthonormal set.

    (2) From $\det(P^TP) = (\det P)^2 = 1$ we get $\det P = \pm 1$.

    (3) The composition of two permutations is again a permutation. $\blacksquare$

!!! theorem "Theorem 10.2 (PLU decomposition)"
    Any $n \times n$ invertible matrix $A$ can be decomposed as

    $$PA = LU$$

    where $P$ is a permutation matrix, $L$ is a unit lower triangular matrix, and $U$ is an upper triangular matrix. Equivalently,

    $$A = P^{-1}LU = P^TLU$$

??? proof "Proof"
    During Gaussian elimination, if a pivot is zero at some step, a row swap (partial pivoting) makes the pivot nonzero. All row swap operations correspond to a permutation matrix $P$. Performing Gaussian elimination on $PA$ without row swaps gives $PA = LU$.

    More precisely: for an invertible matrix $A$, at each elimination step, choose the element with the largest absolute value in the current column as the pivot (partial pivoting strategy), swapping rows to move it to the diagonal position. The cumulative effect of all row swaps is recorded by $P$. Since $A$ is invertible, the elimination process can always be completed. $\blacksquare$

!!! example "Example 10.2"
    Perform PLU decomposition of $A = \begin{pmatrix} 0 & 2 & 1 \\ 1 & 1 & 0 \\ 2 & 0 & 3 \end{pmatrix}$.

    The first-column pivot is $0$, so a row swap is needed. Swap rows 1 and 3:

    $$P_1A = \begin{pmatrix} 2 & 0 & 3 \\ 1 & 1 & 0 \\ 0 & 2 & 1 \end{pmatrix}$$

    $R_2 \leftarrow R_2 - \frac{1}{2}R_1$:

    $$\begin{pmatrix} 2 & 0 & 3 \\ 0 & 1 & -3/2 \\ 0 & 2 & 1 \end{pmatrix}$$

    $R_3 \leftarrow R_3 - 2R_2$:

    $$U = \begin{pmatrix} 2 & 0 & 3 \\ 0 & 1 & -3/2 \\ 0 & 0 & 4 \end{pmatrix}$$

    Therefore $P = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}$, $L = \begin{pmatrix} 1 & 0 & 0 \\ 1/2 & 1 & 0 \\ 0 & 2 & 1 \end{pmatrix}$, $PA = LU$.

---

## 10.3 Cholesky Decomposition

<div class="context-flow" markdown>

Ch9 positive definite $A = C^TC$ → Take $C = L^T$ to get $A = LL^T$ → Cost is only $\frac{1}{3}n^3$ (half of LU), no pivoting needed, numerically stable

</div>

!!! definition "Definition 10.3 (Cholesky decomposition)"
    Let $A$ be an $n \times n$ real symmetric positive definite matrix. The **Cholesky decomposition** of $A$ is

    $$A = LL^T$$

    where $L$ is a lower triangular matrix with all positive diagonal entries. For complex Hermitian positive definite matrices, $A = LL^H$.

!!! theorem "Theorem 10.3 (Existence and uniqueness of Cholesky decomposition)"
    Let $A$ be an $n \times n$ real symmetric positive definite matrix. Then there exists a unique lower triangular matrix $L$ with all positive diagonal entries such that $A = LL^T$.

??? proof "Proof"
    **Existence:** By LU decomposition (positive definite matrices have all positive leading principal minors, guaranteeing LU decomposition exists): $A = L_0 U_0$. Let $D = \operatorname{diag}(u_{11}, \ldots, u_{nn})$ (diagonal of $U_0$), then $U_0 = DU_1$ where $U_1$ is unit upper triangular.

    By $A = A^T$, we have $L_0DU_1 = U_1^TDL_0^T$. By uniqueness, $U_1 = L_0^T$, so $A = L_0DL_0^T$.

    Positive definiteness ensures $u_{ii} = d_i > 0$. Let $L = L_0 D^{1/2}$ (where $D^{1/2} = \operatorname{diag}(\sqrt{d_1}, \ldots, \sqrt{d_n})$), then $A = L_0 D^{1/2}(D^{1/2})^TL_0^T = LL^T$.

    **Uniqueness:** If $A = L_1L_1^T = L_2L_2^T$, then $L_2^{-1}L_1 = L_2^T L_1^{-T}$. The left side is lower triangular, the right side is upper triangular, so both are a diagonal matrix $D$. From $L_2^{-1}L_1 = D$ and $L_2^TL_1^{-T} = D$, we get $D = D^T = D$ and $DD^T = I$. Since diagonal entries are positive, $D = I$, so $L_1 = L_2$. $\blacksquare$

!!! example "Example 10.3"
    Perform Cholesky decomposition of $A = \begin{pmatrix} 4 & 2 & -2 \\ 2 & 10 & 2 \\ -2 & 2 & 5 \end{pmatrix}$.

    Let $L = \begin{pmatrix} l_{11} & 0 & 0 \\ l_{21} & l_{22} & 0 \\ l_{31} & l_{32} & l_{33} \end{pmatrix}$. From $A = LL^T$:

    - $l_{11}^2 = 4 \Rightarrow l_{11} = 2$
    - $l_{21}l_{11} = 2 \Rightarrow l_{21} = 1$
    - $l_{31}l_{11} = -2 \Rightarrow l_{31} = -1$
    - $l_{21}^2 + l_{22}^2 = 10 \Rightarrow l_{22}^2 = 9 \Rightarrow l_{22} = 3$
    - $l_{31}l_{21} + l_{32}l_{22} = 2 \Rightarrow l_{32} = 1$
    - $l_{31}^2 + l_{32}^2 + l_{33}^2 = 5 \Rightarrow l_{33}^2 = 3 \Rightarrow l_{33} = \sqrt{3}$

    $$L = \begin{pmatrix} 2 & 0 & 0 \\ 1 & 3 & 0 \\ -1 & 1 & \sqrt{3} \end{pmatrix}$$

!!! note "Note"
    The computational cost of Cholesky decomposition is approximately $\frac{1}{3}n^3$, half that of LU decomposition. It is very numerically stable and requires no pivoting. Cholesky decomposition is also commonly used to test whether a matrix is positive definite: if the algorithm completes successfully (all diagonal entries positive), then the matrix is positive definite.

---

## 10.4 QR Decomposition

<div class="context-flow" markdown>

Orthogonal $\times$ upper triangular: $A = QR$ → Three construction methods: **Gram-Schmidt** (Ch8 orthogonalization), **Householder** reflections, **Givens** rotations → Core tool for numerical least squares

</div>

!!! definition "Definition 10.4 (QR decomposition)"
    Let $A$ be an $m \times n$ matrix ($m \geq n$). The **QR decomposition** of $A$ is

    $$A = QR$$

    where $Q$ is an $m \times m$ orthogonal (or unitary) matrix and $R$ is an $m \times n$ upper triangular matrix.

    **Thin QR** (reduced form): $A = Q_1R_1$, where $Q_1$ is $m \times n$ (with orthonormal columns, $Q_1^TQ_1 = I_n$) and $R_1$ is $n \times n$ upper triangular.

!!! theorem "Theorem 10.4 (Existence of QR decomposition)"
    Any $m \times n$ matrix $A$ ($m \geq n$) has a QR decomposition. If $A$ has full column rank and the diagonal entries of $R$ are required to be positive, then the thin QR decomposition is unique.

??? proof "Proof"
    **Existence (Gram-Schmidt method):** Let $A = [\mathbf{a}_1 | \cdots | \mathbf{a}_n]$. Apply Gram-Schmidt orthogonalization to the columns $\mathbf{a}_1, \ldots, \mathbf{a}_n$:

    $$\mathbf{u}_1 = \mathbf{a}_1, \quad \mathbf{e}_1 = \frac{\mathbf{u}_1}{\|\mathbf{u}_1\|}$$

    $$\mathbf{u}_k = \mathbf{a}_k - \sum_{j=1}^{k-1}\langle \mathbf{a}_k, \mathbf{e}_j\rangle \mathbf{e}_j, \quad \mathbf{e}_k = \frac{\mathbf{u}_k}{\|\mathbf{u}_k\|}$$

    Then $\mathbf{a}_k = \sum_{j=1}^{k}\langle \mathbf{a}_k, \mathbf{e}_j\rangle \mathbf{e}_j$; in matrix form this is $A = Q_1R_1$.

    **Uniqueness:** If $A = Q_1R_1 = Q_2R_2$ ($R_i$ with positive diagonal entries), then $Q_2^TQ_1 = R_2R_1^{-1}$. The left side has orthonormal columns; the right side is upper triangular. $Q_2^TQ_1$ is simultaneously an orthogonal matrix and an upper triangular matrix, hence a diagonal matrix with entries $\pm 1$. Since $R_1, R_2$ have positive diagonal entries, this diagonal matrix is $I$. $\blacksquare$

### Gram-Schmidt Orthogonalization

!!! definition "Definition 10.5 (Gram-Schmidt orthogonalization)"
    Let $\{\mathbf{a}_1, \ldots, \mathbf{a}_n\}$ be a linearly independent set of vectors. The **Gram-Schmidt orthogonalization** process is:

    $$\mathbf{u}_1 = \mathbf{a}_1, \quad \mathbf{e}_1 = \frac{\mathbf{u}_1}{\|\mathbf{u}_1\|}$$

    For $k = 2, \ldots, n$:

    $$\mathbf{u}_k = \mathbf{a}_k - \sum_{j=1}^{k-1} \frac{\langle \mathbf{a}_k, \mathbf{u}_j \rangle}{\langle \mathbf{u}_j, \mathbf{u}_j \rangle} \mathbf{u}_j, \quad \mathbf{e}_k = \frac{\mathbf{u}_k}{\|\mathbf{u}_k\|}$$

    Then $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ is an orthonormal set, and for each $k$, $\operatorname{span}\{\mathbf{e}_1, \ldots, \mathbf{e}_k\} = \operatorname{span}\{\mathbf{a}_1, \ldots, \mathbf{a}_k\}$.

### Householder Transformation

!!! definition "Definition 10.6 (Householder reflection)"
    Let $\mathbf{v} \in \mathbb{R}^n$ be a nonzero vector. The **Householder reflection** is defined as

    $$H = I - 2\frac{\mathbf{v}\mathbf{v}^T}{\mathbf{v}^T\mathbf{v}}$$

    $H$ is a symmetric orthogonal matrix ($H^T = H$, $H^2 = I$), and geometrically represents reflection across the hyperplane orthogonal to $\mathbf{v}$.

!!! proposition "Proposition 10.2"
    The Householder reflection matrix $H$ satisfies:

    1. $H = H^T$ (symmetric);
    2. $H^2 = I$ (involutory);
    3. $H$ is an orthogonal matrix, $\det(H) = -1$.

??? proof "Proof"
    (1) $(I - 2\frac{\mathbf{v}\mathbf{v}^T}{\mathbf{v}^T\mathbf{v}})^T = I - 2\frac{\mathbf{v}\mathbf{v}^T}{\mathbf{v}^T\mathbf{v}} = H$.

    (2) $H^2 = (I - 2\frac{\mathbf{v}\mathbf{v}^T}{\|\mathbf{v}\|^2})^2 = I - 4\frac{\mathbf{v}\mathbf{v}^T}{\|\mathbf{v}\|^2} + 4\frac{\mathbf{v}(\mathbf{v}^T\mathbf{v})\mathbf{v}^T}{\|\mathbf{v}\|^4} = I$.

    (3) From (1) and (2), $HH^T = H^2 = I$. $\det(H) = -1$ because $H$ has one eigenvalue $-1$ and the remaining $n-1$ eigenvalues are $1$. $\blacksquare$

!!! theorem "Theorem 10.5 (Householder QR)"
    For any $m \times n$ matrix $A$ ($m \geq n$), there exist Householder reflections $H_1, H_2, \ldots, H_n$ such that

    $$H_n \cdots H_2 H_1 A = R$$

    is upper triangular. Therefore $A = QR$, where $Q = H_1 H_2 \cdots H_n$.

??? proof "Proof"
    **Constructing $H_1$:** Let $\mathbf{a}_1$ be the first column of $A$. Choose $\mathbf{v}_1 = \mathbf{a}_1 + \operatorname{sign}(a_{11})\|\mathbf{a}_1\|\mathbf{e}_1$, and let $H_1 = I - 2\frac{\mathbf{v}_1\mathbf{v}_1^T}{\mathbf{v}_1^T\mathbf{v}_1}$. Then $H_1\mathbf{a}_1 = -\operatorname{sign}(a_{11})\|\mathbf{a}_1\|\mathbf{e}_1$, so the first column is reduced to having only a nonzero first entry.

    Repeat this process on the lower-right $(m-1)\times(n-1)$ submatrix of $H_1A$. Embed the $(m-1)$-dimensional Householder matrix into $m$ dimensions (padding with $1$ in the upper-left corner) to get $H_2$. Iterate for $n$ steps. $\blacksquare$

### Givens Rotation

!!! definition "Definition 10.7 (Givens rotation)"
    A **Givens rotation** $G(i, j, \theta)$ is the identity matrix modified at positions $(i,i), (i,j), (j,i), (j,j)$:

    $$g_{ii} = g_{jj} = \cos\theta, \quad g_{ij} = -\sin\theta, \quad g_{ji} = \sin\theta$$

    All other positions are unchanged. $G(i,j,\theta)$ is an orthogonal matrix representing rotation by angle $\theta$ in the $(i,j)$ plane.

!!! example "Example 10.4"
    Use the Gram-Schmidt method to compute the QR decomposition of $A = \begin{pmatrix} 1 & 1 & 0 \\ 1 & 0 & 1 \\ 0 & 1 & 1 \end{pmatrix}$.

    $\mathbf{a}_1 = (1,1,0)^T$, $\|\mathbf{a}_1\| = \sqrt{2}$, $\mathbf{e}_1 = \frac{1}{\sqrt{2}}(1,1,0)^T$.

    $\mathbf{u}_2 = \mathbf{a}_2 - \langle \mathbf{a}_2, \mathbf{e}_1\rangle \mathbf{e}_1 = (1,0,1)^T - \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}}(1,1,0)^T = (1,0,1)^T - \frac{1}{2}(1,1,0)^T = \frac{1}{2}(1,-1,2)^T$.

    $\|\mathbf{u}_2\| = \frac{1}{2}\sqrt{6}$, $\mathbf{e}_2 = \frac{1}{\sqrt{6}}(1,-1,2)^T$.

    $\mathbf{u}_3 = \mathbf{a}_3 - \langle \mathbf{a}_3, \mathbf{e}_1\rangle \mathbf{e}_1 - \langle \mathbf{a}_3, \mathbf{e}_2\rangle \mathbf{e}_2$.

    $\langle \mathbf{a}_3, \mathbf{e}_1\rangle = \frac{1}{\sqrt{2}}$, $\langle \mathbf{a}_3, \mathbf{e}_2\rangle = \frac{1}{\sqrt{6}}(0-1+2) = \frac{1}{\sqrt{6}}$.

    $\mathbf{u}_3 = (0,1,1)^T - \frac{1}{2}(1,1,0)^T - \frac{1}{6}(1,-1,2)^T = \frac{1}{3}(-1,1,1)^T \cdot \frac{3}{2}$.

    After computation, $\mathbf{u}_3 = (-\frac{2}{3}, \frac{2}{3}, \frac{2}{3})^T$, $\|\mathbf{u}_3\| = \frac{2}{\sqrt{3}}$, $\mathbf{e}_3 = \frac{1}{\sqrt{3}}(-1, 1, 1)^T$.

    $$Q = \begin{pmatrix} 1/\sqrt{2} & 1/\sqrt{6} & -1/\sqrt{3} \\ 1/\sqrt{2} & -1/\sqrt{6} & 1/\sqrt{3} \\ 0 & 2/\sqrt{6} & 1/\sqrt{3} \end{pmatrix}$$

    $$R = Q^TA = \begin{pmatrix} \sqrt{2} & 1/\sqrt{2} & 1/\sqrt{2} \\ 0 & \sqrt{6}/2 & 1/\sqrt{6} \\ 0 & 0 & 2/\sqrt{3} \end{pmatrix}$$

---

## 10.5 Schur Decomposition

<div class="context-flow" markdown>

Any square matrix can be **unitarily triangularized** $A = UTU^H$ → Diagonal entries = eigenvalues → More general than diagonalization (no diagonalizability required) → For normal matrices, $T$ reduces to diagonal

</div>

<div class="context-flow" markdown>

**Insight**: The existence proof of the Schur decomposition is structurally identical to the spectral theorem — find an eigenvector, invariant orthogonal complement, inductive dimension reduction — but without requiring normality

</div>

!!! theorem "Theorem 10.6 (Schur decomposition)"
    Let $A$ be an $n \times n$ complex matrix. Then there exists a unitary matrix $U$ such that

    $$U^H A U = T$$

    where $T$ is an upper triangular matrix whose diagonal entries are the eigenvalues of $A$ (in any specified order). Equivalently, $A = UTU^H$.

??? proof "Proof"
    By induction on $n$. The case $n = 1$ is trivial. Assume the result holds for $(n-1) \times (n-1)$ matrices.

    By the fundamental theorem of algebra, $A$ has an eigenvalue $\lambda_1$. Let $\mathbf{u}_1$ be a corresponding unit eigenvector. Extend $\mathbf{u}_1$ to an orthonormal basis $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$ of $\mathbb{C}^n$, and let $U_1 = [\mathbf{u}_1 | \cdots | \mathbf{u}_n]$. Then

    $$U_1^H A U_1 = \begin{pmatrix} \lambda_1 & \mathbf{b}^H \\ \mathbf{0} & A_1 \end{pmatrix}$$

    where $A_1$ is an $(n-1) \times (n-1)$ matrix (the first column is the zero vector because $\mathbf{u}_j^HA\mathbf{u}_1 = \lambda_1 \mathbf{u}_j^H\mathbf{u}_1 = 0$ for $j \geq 2$).

    By the induction hypothesis, there exists an $(n-1) \times (n-1)$ unitary matrix $V_1$ such that $V_1^HA_1V_1 = T_1$ is upper triangular. Let

    $$V = \begin{pmatrix} 1 & \mathbf{0}^H \\ \mathbf{0} & V_1 \end{pmatrix}$$

    Then $U = U_1V$ is unitary and $U^HAU$ is upper triangular with the eigenvalues of $A$ on the diagonal. $\blacksquare$

!!! theorem "Theorem 10.7 (Real Schur decomposition)"
    Let $A$ be an $n \times n$ real matrix. Then there exists an orthogonal matrix $Q$ such that

    $$Q^TAQ = T$$

    where $T$ is a quasi-upper triangular matrix, i.e., diagonal blocks are $1\times 1$ blocks (corresponding to real eigenvalues) or $2\times 2$ blocks $\begin{pmatrix} a & b \\ c & d \end{pmatrix}$ (corresponding to conjugate complex eigenvalue pairs).

??? proof "Proof"
    The idea is similar to the complex Schur decomposition. For real eigenvalues, the argument is the same as in the complex case. For conjugate complex eigenvalue pairs $\lambda = a \pm bi$, the corresponding conjugate complex eigenvectors $\mathbf{v} \pm i\mathbf{w}$ have real and imaginary parts $\mathbf{v}, \mathbf{w}$ that span a $2$-dimensional real invariant subspace, and the restriction of $A$ to this subspace gives the $2\times 2$ block. $\blacksquare$

!!! example "Example 10.5"
    For $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$ (the 90-degree rotation matrix), the eigenvalues are $\pm i$.

    Complex Schur decomposition: $U = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ i & -i \end{pmatrix}$, $U^HAU = \begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}$.

    Real Schur decomposition: $A$ itself is already in quasi-upper triangular form (a single $2\times 2$ block): $Q = I$, $T = A$.

---

## 10.6 Spectral Decomposition

<div class="context-flow" markdown>

$A = \sum \lambda_i P_i$ (eigenvalue $\times$ orthogonal projection sum) → Decompose the matrix into a **superposition of independent modes** → Ch11 SVD is the non-symmetric version, Ch13 $f(A) = \sum f(\lambda_i)P_i$

</div>

!!! definition "Definition 10.8 (Spectral decomposition)"
    Let $A$ be a diagonalizable matrix. If $A$ has eigenvalues $\lambda_1, \ldots, \lambda_k$ (distinct), with corresponding eigenspace projections $P_1, \ldots, P_k$, then

    $$A = \lambda_1 P_1 + \lambda_2 P_2 + \cdots + \lambda_k P_k$$

    is called the **spectral decomposition** of $A$.

<div class="context-flow" markdown>

**Insight**: In $A = \sum \lambda_i \mathbf{q}_i\mathbf{q}_i^T$, each $\mathbf{q}_i\mathbf{q}_i^T$ is a rank-one orthogonal projection — the matrix is decomposed into "frequency components," the finite-dimensional version of Fourier analysis

</div>

!!! theorem "Theorem 10.8 (Spectral decomposition of symmetric matrices)"
    Let $A$ be an $n \times n$ real symmetric matrix with eigenvalues $\lambda_1, \ldots, \lambda_n$ (repetitions allowed) and corresponding orthonormal eigenvectors $\mathbf{q}_1, \ldots, \mathbf{q}_n$. Then

    $$A = Q\Lambda Q^T = \sum_{i=1}^n \lambda_i \mathbf{q}_i \mathbf{q}_i^T$$

    Combining terms with equal eigenvalues, let the distinct eigenvalues be $\mu_1, \ldots, \mu_k$ with corresponding eigenspace projections $P_j = \sum_{\lambda_i = \mu_j} \mathbf{q}_i\mathbf{q}_i^T$. Then

    $$A = \sum_{j=1}^k \mu_j P_j, \quad \text{where} \quad \sum_{j=1}^k P_j = I, \quad P_j P_l = \delta_{jl} P_j$$

??? proof "Proof"
    By the spectral theorem, $A = Q\Lambda Q^T$. Expanding the matrix product:

    $$A = \sum_{i=1}^n \lambda_i \mathbf{q}_i \mathbf{q}_i^T$$

    Verifying the properties of $P_j$: $P_j$ is the orthogonal projection onto the eigenspace $E_{\mu_j}$, so $P_j^2 = P_j$, $P_j^T = P_j$, $P_jP_l = 0$ ($j \neq l$, since distinct eigenspaces are orthogonal), $\sum P_j = I$ (since the eigenvectors form a complete orthonormal basis). $\blacksquare$

!!! theorem "Theorem 10.9 (Spectral decomposition of normal matrices)"
    Let $A$ be an $n \times n$ normal matrix with eigenvalues $\lambda_1, \ldots, \lambda_n$ and corresponding orthonormal eigenvectors $\mathbf{u}_1, \ldots, \mathbf{u}_n$. Then

    $$A = U\Lambda U^H = \sum_{i=1}^n \lambda_i \mathbf{u}_i \mathbf{u}_i^H$$

??? proof "Proof"
    By the spectral theorem for normal matrices (Theorem 8.17), $A$ is unitarily diagonalizable as $A = U\Lambda U^H$. Expanding gives the result. $\blacksquare$

!!! example "Example 10.6"
    Compute the spectral decomposition of $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$.

    Eigenvalues $\lambda_1 = 1$, $\mathbf{q}_1 = \frac{1}{\sqrt{2}}(1,-1)^T$; $\lambda_2 = 3$, $\mathbf{q}_2 = \frac{1}{\sqrt{2}}(1,1)^T$.

    $$P_1 = \mathbf{q}_1\mathbf{q}_1^T = \frac{1}{2}\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}, \quad P_2 = \mathbf{q}_2\mathbf{q}_2^T = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$$

    $$A = 1 \cdot P_1 + 3 \cdot P_2 = \frac{1}{2}\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix} + \frac{3}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$$

    Verification: $P_1 + P_2 = I$, $P_1P_2 = 0$.

---

## 10.7 Polar Decomposition

<div class="context-flow" markdown>

$A = UP$: unitary $\times$ positive definite $\leftrightarrow$ complex polar coordinates $z = e^{i\theta}|z|$ → Computed via SVD: $U = VW^H$, $P = W\Sigma W^H$ → Separates "rotation" and "stretching"

</div>

!!! definition "Definition 10.9 (Polar decomposition)"
    Let $A$ be an $n \times n$ invertible complex matrix. The **polar decomposition** of $A$ is

    $$A = UP$$

    where $U$ is a unitary matrix ($U^HU = I$) and $P$ is a positive definite Hermitian matrix ($P^H = P$, $P > 0$). This is called the **right polar decomposition**. Similarly, $A = P'U$ ($P'$ positive definite) is the **left polar decomposition**.

!!! theorem "Theorem 10.10 (Existence and uniqueness of polar decomposition)"
    (1) Any $n \times n$ invertible matrix $A$ has a unique right polar decomposition $A = UP$ ($U$ unitary, $P$ positive definite Hermitian) and a unique left polar decomposition $A = P'U$.

    (2) More generally, any $m \times n$ matrix $A$ can be decomposed as $A = UP$, where $U$ has orthonormal columns ($U^HU = I_n$) and $P$ is a positive semidefinite Hermitian matrix. $P$ is uniquely determined ($P = \sqrt{A^HA}$), but $U$ is not unique when $A$ is rank-deficient.

??? proof "Proof"
    **Existence:** $A^HA$ is a positive semidefinite Hermitian matrix. Let $P = (A^HA)^{1/2}$ (the positive semidefinite Hermitian square root, constructible via the spectral theorem).

    If $A$ is invertible, then $A^HA$ is positive definite, so $P$ is positive definite and invertible. Let $U = AP^{-1}$. Verify that $U$ is unitary:

    $$U^HU = (AP^{-1})^H(AP^{-1}) = P^{-H}A^HAP^{-1} = P^{-1}P^2P^{-1} = I$$

    **Uniqueness (invertible case):** If $A = U_1P_1 = U_2P_2$, then $A^HA = P_1^2 = P_2^2$. By the uniqueness of the positive definite square root (see Theorem 13.7), $P_1 = P_2$. Hence $U_1 = U_2$.

    **Left polar decomposition:** $P' = (AA^H)^{1/2}$, $U = (P')^{-1}A$. $\blacksquare$

!!! note "Note"
    The name "polar decomposition" comes from the analogy with complex polar representation $z = e^{i\theta}|z|$: $U$ corresponds to "rotation" ($e^{i\theta}$) and $P$ corresponds to "stretching" ($|z|$). The polar decomposition can be conveniently computed via SVD: if $A = V\Sigma W^H$ (SVD), then $U = VW^H$, $P = W\Sigma W^H$.

!!! example "Example 10.7"
    Compute the polar decomposition of $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$.

    $A^TA = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$, eigenvalues $\frac{3 \pm \sqrt{5}}{2}$.

    Let $\lambda_1 = \frac{3+\sqrt{5}}{2}$, $\lambda_2 = \frac{3-\sqrt{5}}{2}$. $P = (A^TA)^{1/2}$, i.e., the eigenvalues of $P$ are $\sqrt{\lambda_1}, \sqrt{\lambda_2}$.

    $P = Q\operatorname{diag}(\sqrt{\lambda_1}, \sqrt{\lambda_2})Q^T$ ($Q$ consists of orthonormal eigenvectors of $A^TA$).

    $U = AP^{-1}$ is an orthogonal matrix.

    Numerically: $P \approx \begin{pmatrix} 1.0515 & 0.4851 \\ 0.4851 & 1.2965 \end{pmatrix}$, $U \approx \begin{pmatrix} 0.8507 & 0.5257 \\ -0.5257 & 0.8507 \end{pmatrix}$.
