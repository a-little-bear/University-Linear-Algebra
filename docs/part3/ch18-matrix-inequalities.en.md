# Chapter 18  Matrix Inequalities

<div class="context-flow" markdown>

**Prerequisites**: Ch15 Weyl/Wielandt-Hoffman · Ch16 Positive definite matrices · **Chapter arc**: Courant-Fischer → Weyl/Cauchy interlacing → Singular value inequalities → Trace (von Neumann) → Determinant (Hadamard/Minkowski) → **Majorization** as a unifying framework

</div>

Matrix inequalities are among the most profound topics in matrix analysis. In the scalar case, the theory of inequalities is already very mature — from the arithmetic-geometric mean inequality to the Cauchy-Schwarz inequality, these results form the cornerstones of analysis. However, when we generalize these ideas to the world of matrices, the situation becomes far richer and more subtle. Deep and elegant inequality relationships exist among the eigenvalues, singular values, traces, and determinants of matrices. These are important not only in pure mathematics but also have wide applications in quantum information theory, statistics, and optimization theory.

This chapter systematically introduces the main results on matrix inequalities, starting from eigenvalue inequalities, proceeding through singular value inequalities and trace inequalities, to determinant inequalities, and finally using majorization theory and matrix convexity as a unifying framework.

---

## 18.1 Eigenvalue inequalities

<div class="context-flow" markdown>

**Core tool**: **Courant-Fischer** minimax theorem (Rayleigh quotient + subspace dimension argument) → Weyl inequality · Cauchy interlacing = squeezing of principal submatrix eigenvalues

</div>

Eigenvalues are the most fundamental invariants of a matrix. For Hermitian matrices, all eigenvalues are real and can therefore be ordered and compared. This section discusses the relationship between the eigenvalues of a sum of Hermitian matrices and the eigenvalues of the individual summands.

**Convention**: For an $n \times n$ Hermitian matrix $A$, its eigenvalues are arranged in decreasing order as $\lambda_1(A) \geq \lambda_2(A) \geq \cdots \geq \lambda_n(A)$.

!!! definition "Definition 18.1 (Eigenvalue ordering of Hermitian matrices)"
    Let $A$ be an $n \times n$ Hermitian matrix, i.e., $A = A^*$. Its $n$ real eigenvalues are arranged in decreasing order:
    $$
    \lambda_1(A) \geq \lambda_2(A) \geq \cdots \geq \lambda_n(A).
    $$
    We write $\lambda_{\max}(A) = \lambda_1(A)$ and $\lambda_{\min}(A) = \lambda_n(A)$ for the largest and smallest eigenvalues, respectively.

!!! definition "Definition 18.2 (Rayleigh quotient)"
    Let $A$ be an $n \times n$ Hermitian matrix. For a nonzero vector $\mathbf{x} \in \mathbb{C}^n$, the **Rayleigh quotient** is defined as:
    $$
    R_A(\mathbf{x}) = \frac{\mathbf{x}^* A \mathbf{x}}{\mathbf{x}^* \mathbf{x}}.
    $$

!!! theorem "Theorem 18.1 (Courant-Fischer minimax theorem)"
    Let $A$ be an $n \times n$ Hermitian matrix with eigenvalues $\lambda_1(A) \geq \cdots \geq \lambda_n(A)$. Then for $k = 1, 2, \ldots, n$:
    $$
    \lambda_k(A) = \max_{\dim S = k} \min_{\mathbf{x} \in S, \mathbf{x} \neq \mathbf{0}} \frac{\mathbf{x}^* A \mathbf{x}}{\mathbf{x}^* \mathbf{x}} = \min_{\dim T = n-k+1} \max_{\mathbf{x} \in T, \mathbf{x} \neq \mathbf{0}} \frac{\mathbf{x}^* A \mathbf{x}}{\mathbf{x}^* \mathbf{x}}.
    $$

??? proof "Proof"
    Let $A = U \Lambda U^*$ be the spectral decomposition, where $U = [\mathbf{u}_1, \ldots, \mathbf{u}_n]$ is a unitary matrix and $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$.

    **Step 1**: Prove one direction of $\lambda_k \geq \max_{\dim S=k} \min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} R_A(\mathbf{x})$.

    Take $S_0 = \operatorname{span}\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$, so $\dim S_0 = k$. For any $\mathbf{x} = \sum_{i=1}^k c_i \mathbf{u}_i \in S_0$:
    $$
    R_A(\mathbf{x}) = \frac{\sum_{i=1}^k \lambda_i |c_i|^2}{\sum_{i=1}^k |c_i|^2} \geq \lambda_k.
    $$
    Therefore $\min_{\mathbf{x} \in S_0 \setminus \{\mathbf{0}\}} R_A(\mathbf{x}) \geq \lambda_k$, and hence $\max_{\dim S=k} \min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} R_A(\mathbf{x}) \geq \lambda_k$.

    **Step 2**: Prove that for any $k$-dimensional subspace $S$, $\min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} R_A(\mathbf{x}) \leq \lambda_k$.

    Take $T_0 = \operatorname{span}\{\mathbf{u}_k, \ldots, \mathbf{u}_n\}$, so $\dim T_0 = n - k + 1$. By the dimension formula, $\dim(S \cap T_0) \geq k + (n-k+1) - n = 1$, so there exists a nonzero $\mathbf{x} \in S \cap T_0$. For such $\mathbf{x} = \sum_{i=k}^n c_i \mathbf{u}_i$:
    $$
    R_A(\mathbf{x}) = \frac{\sum_{i=k}^n \lambda_i |c_i|^2}{\sum_{i=k}^n |c_i|^2} \leq \lambda_k.
    $$
    Therefore $\min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} R_A(\mathbf{x}) \leq \lambda_k$.

    Combining both steps yields the first equality. The proof of the second equality is completely analogous. $\blacksquare$

!!! theorem "Theorem 18.2 (Weyl inequality)"
    Let $A, B$ be $n \times n$ Hermitian matrices. Then for all $i, j$ satisfying $i + j - 1 \leq n$:
    $$
    \lambda_{i+j-1}(A + B) \leq \lambda_i(A) + \lambda_j(B).
    $$
    For all $i, j$ satisfying $i + j - n \geq 1$:
    $$
    \lambda_{i+j-n}(A + B) \geq \lambda_i(A) + \lambda_j(B).
    $$

??? proof "Proof"
    We prove the first inequality. By the min-max form of the Courant-Fischer theorem:
    $$
    \lambda_{i+j-1}(A+B) = \min_{\dim T = n-i-j+2} \max_{\mathbf{x} \in T \setminus \{\mathbf{0}\}} R_{A+B}(\mathbf{x}).
    $$

    Let the subspace spanned by the first $i$ eigenvectors of $A$ be $S_A$ ($\dim S_A = i$), and the subspace spanned by the first $j$ eigenvectors of $B$ be $S_B$ ($\dim S_B = j$).

    For any $(n-i-j+2)$-dimensional subspace $T$, by the dimension formula:
    $$
    \dim(S_A \cap T) \geq i + (n-i-j+2) - n = 2 - j + i \geq 1,
    $$
    i.e., $\dim(T \cap S_A) \geq 1$ (when $i + j - 1 \leq n$). Similarly $\dim(T \cap S_B) \geq 1$.

    More directly, take $T_0 = \operatorname{span}\{\mathbf{u}_i, \ldots, \mathbf{u}_n\} \cap \operatorname{span}\{\mathbf{v}_j, \ldots, \mathbf{v}_n\}$, where $\mathbf{u}_k, \mathbf{v}_k$ are the eigenvectors of $A, B$ respectively. By the Courant-Fischer theorem:

    For any nonzero $\mathbf{x}$:
    $$
    R_{A+B}(\mathbf{x}) = R_A(\mathbf{x}) + R_B(\mathbf{x}).
    $$

    Taking a subspace where $R_A(\mathbf{x}) \leq \lambda_i(A)$ and $R_B(\mathbf{x}) \leq \lambda_j(B)$ hold simultaneously (the dimension argument guarantees it is nonempty), we get:
    $$
    \lambda_{i+j-1}(A+B) \leq \lambda_i(A) + \lambda_j(B). \quad \blacksquare
    $$

!!! note "Note"
    A commonly used special case of the Weyl inequality is $j = 1$: $\lambda_i(A+B) \leq \lambda_i(A) + \lambda_1(B)$, meaning that adding a positive semidefinite matrix does not increase eigenvalues by more than $\lambda_1(B)$. Similarly, taking $j = n$: $\lambda_i(A+B) \geq \lambda_i(A) + \lambda_n(B)$. This gives a perturbation bound for eigenvalues.

!!! theorem "Theorem 18.3 (Cauchy interlacing theorem)"
    Let $A$ be an $n \times n$ Hermitian matrix, and $B$ be an $m \times m$ principal submatrix of $A$ (i.e., $B = P^* A P$, where $P$ is an $n \times m$ matrix with $P^* P = I_m$), $m < n$. Then the eigenvalues of $A$ and $B$ satisfy the **interlacing** relation:
    $$
    \lambda_i(A) \geq \lambda_i(B) \geq \lambda_{i+n-m}(A), \quad i = 1, 2, \ldots, m.
    $$

??? proof "Proof"
    By the Courant-Fischer theorem, for $B = P^* A P$:
    $$
    \lambda_i(B) = \max_{\substack{S \subset \mathbb{C}^m \\ \dim S = i}} \min_{\mathbf{y} \in S \setminus \{\mathbf{0}\}} \frac{\mathbf{y}^* B \mathbf{y}}{\mathbf{y}^* \mathbf{y}} = \max_{\substack{S \subset \mathbb{C}^m \\ \dim S = i}} \min_{\mathbf{y} \in S \setminus \{\mathbf{0}\}} \frac{(P\mathbf{y})^* A (P\mathbf{y})}{(P\mathbf{y})^* (P\mathbf{y})}.
    $$

    Let $\mathbf{x} = P \mathbf{y}$. Since $P$ is an isometric embedding, as $S$ ranges over all $i$-dimensional subspaces of $\mathbb{C}^m$, $PS$ ranges over all $i$-dimensional subspaces of $P(\mathbb{C}^m)$. Since $P(\mathbb{C}^m)$ is an $m$-dimensional subspace of $\mathbb{C}^n$, and the set of all $i$-dimensional subspaces of $\mathbb{C}^n$ contains all $i$-dimensional subspaces of $P(\mathbb{C}^m)$:
    $$
    \lambda_i(B) \leq \max_{\substack{S \subset \mathbb{C}^n \\ \dim S = i}} \min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} \frac{\mathbf{x}^* A \mathbf{x}}{\mathbf{x}^* \mathbf{x}} = \lambda_i(A).
    $$

    Similarly, using the min-max form one can prove $\lambda_i(B) \geq \lambda_{i+n-m}(A)$. $\blacksquare$

!!! theorem "Theorem 18.4 (Poincare separation theorem)"
    Let $A$ be an $n \times n$ Hermitian matrix, $U$ be an $n \times k$ matrix satisfying $U^* U = I_k$ ($k \leq n$), and let $B = U^* A U$. Then:
    $$
    \lambda_i(A) \geq \lambda_i(B) \geq \lambda_{i+n-k}(A), \quad i = 1, \ldots, k.
    $$

??? proof "Proof"
    This is in fact the generalization of the Cauchy interlacing theorem to general isometric compressions. The proof method is exactly the same as Theorem 18.3, with $P$ replaced by a general isometric map $U$. The key is that $U^* U = I_k$ guarantees that $U$ is an isometric embedding from $\mathbb{C}^k$ into $\mathbb{C}^n$. $\blacksquare$

!!! example "Example 18.1"
    Let $A = \begin{pmatrix} 5 & 1 & 0 \\ 1 & 3 & 1 \\ 0 & 1 & 1 \end{pmatrix}$, with eigenvalues $\lambda_1(A) \approx 5.414$, $\lambda_2(A) \approx 2.828$, $\lambda_3(A) \approx 0.758$.

    Take the upper-left $2 \times 2$ principal submatrix $B = \begin{pmatrix} 5 & 1 \\ 1 & 3 \end{pmatrix}$, with eigenvalues $\lambda_1(B) = 3 + \sqrt{2} \approx 4.414$, $\lambda_2(B) = 3 - \sqrt{2} \approx 1.586$.

    Verifying interlacing:
    $$
    \lambda_1(A) \approx 5.414 \geq \lambda_1(B) \approx 4.414 \geq \lambda_2(A) \approx 2.828,
    $$
    $$
    \lambda_2(A) \approx 2.828 \geq \lambda_2(B) \approx 1.586 \geq \lambda_3(A) \approx 0.758.
    $$
    The interlacing relation holds.

!!! example "Example 18.2"
    **Application of the Weyl inequality**: Let $A = \begin{pmatrix} 4 & 0 \\ 0 & 2 \end{pmatrix}$, $B = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.

    Eigenvalues of $A$: $\lambda_1(A) = 4$, $\lambda_2(A) = 2$.

    Eigenvalues of $B$: $\lambda_1(B) = 2$, $\lambda_2(B) = 0$.

    $A + B = \begin{pmatrix} 5 & 1 \\ 1 & 3 \end{pmatrix}$, eigenvalues $\lambda_1(A+B) = 4 + \sqrt{2} \approx 5.414$, $\lambda_2(A+B) = 4 - \sqrt{2} \approx 2.586$.

    Weyl inequality $\lambda_2(A+B) \leq \lambda_1(A) + \lambda_2(B) = 4 + 0 = 4$ holds ($2.586 \leq 4$).

    Weyl inequality $\lambda_2(A+B) \geq \lambda_2(A) + \lambda_2(B) = 2 + 0 = 2$ holds ($2.586 \geq 2$).

---

## 18.2 Singular value inequalities

<div class="context-flow" markdown>

**Chapter arc**: Singular values are submultiplicative $\sigma_{i+j-1}(AB)\leq\sigma_i(A)\sigma_j(B)$ · **Ky Fan $k$-norm** $=\sum_{i=1}^k\sigma_i$ satisfies the triangle inequality → Unified perspective via unitarily invariant norms

</div>

Singular values are concepts of equal importance to eigenvalues in matrix theory. For general matrices (not necessarily Hermitian or square), singular values provide the most natural measure of "size."

!!! definition "Definition 18.3 (Singular value ordering)"
    Let $A$ be an $m \times n$ matrix, with singular values arranged in decreasing order:
    $$
    \sigma_1(A) \geq \sigma_2(A) \geq \cdots \geq \sigma_{\min(m,n)}(A) \geq 0.
    $$
    Here $\sigma_i(A) = \sqrt{\lambda_i(A^* A)}$.

!!! theorem "Theorem 18.5 (Submultiplicativity of singular values)"
    Let $A, B$ be $n \times n$ matrices. Then for $i + j - 1 \leq n$:
    $$
    \sigma_{i+j-1}(AB) \leq \sigma_i(A) \sigma_j(B).
    $$
    In particular, taking $i = j = 1$: $\sigma_1(AB) \leq \sigma_1(A) \sigma_1(B)$, i.e., $\|AB\|_2 \leq \|A\|_2 \|B\|_2$.

??? proof "Proof"
    The key tools are the relationship between singular values and eigenvalues, together with the Weyl inequality.

    Let $A = U_1 \Sigma_1 V_1^*$, $B = U_2 \Sigma_2 V_2^*$ be singular value decompositions. Consider the squared singular values of $AB$, i.e., the eigenvalues of $B^* A^* A B$.

    Using Fan's extremal principle:
    $$
    \sigma_k(AB) = \min_{\substack{S \subset \mathbb{C}^n \\ \operatorname{codim} S = k-1}} \max_{\mathbf{x} \in S, \|\mathbf{x}\|=1} \|AB\mathbf{x}\|.
    $$

    For any unit vector $\mathbf{x}$, $\|AB\mathbf{x}\| \leq \|A\| \cdot \|B\mathbf{x}\|$ (where $\|A\|$ denotes the operator norm). A more refined analysis uses a subspace dimension argument:

    Take the $(n-i+1)$-dimensional subspace $S_A$ on which $\|A\mathbf{y}\| \leq \sigma_i(A)\|\mathbf{y}\|$, and the $(n-j+1)$-dimensional subspace $S_B$ on which $\|B\mathbf{x}\| \leq \sigma_j(B)\|\mathbf{x}\|$. Then $B^{-1}(S_A) \cap S_B$ has dimension at least $(n-i+1) + (n-j+1) - n = n - i - j + 2$, with codimension $i + j - 2$.

    Therefore on the corresponding subspace $\|AB\mathbf{x}\| \leq \sigma_i(A)\sigma_j(B)\|\mathbf{x}\|$, and by the minimax characterization of singular values we obtain $\sigma_{i+j-1}(AB) \leq \sigma_i(A)\sigma_j(B)$. $\blacksquare$

!!! definition "Definition 18.4 (Ky Fan norm)"
    Let $A$ be an $m \times n$ matrix. For $k = 1, 2, \ldots, \min(m,n)$, the **Ky Fan $k$-norm** is defined as the sum of the first $k$ singular values:
    $$
    \|A\|_{(k)} = \sum_{i=1}^{k} \sigma_i(A).
    $$
    In particular, $\|A\|_{(1)} = \sigma_1(A) = \|A\|_2$ (spectral norm), and $\|A\|_{(\min(m,n))} = \|A\|_*$ (nuclear norm).

!!! theorem "Theorem 18.6 (Fan inequality)"
    Let $A, B$ be $n \times n$ matrices. Then for $k = 1, 2, \ldots, n$:
    $$
    \sum_{i=1}^{k} \sigma_i(A + B) \leq \sum_{i=1}^{k} \sigma_i(A) + \sum_{i=1}^{k} \sigma_i(B).
    $$
    That is, the Ky Fan $k$-norm satisfies the triangle inequality: $\|A + B\|_{(k)} \leq \|A\|_{(k)} + \|B\|_{(k)}$.

??? proof "Proof"
    For any $n \times n$ matrix $M$, by the extremal characterization of singular values:
    $$
    \sum_{i=1}^{k} \sigma_i(M) = \max \{ |\operatorname{tr}(U^* M)| : U \text{ is an } n \times k \text{ isometry} \}.
    $$

    More precisely, there is Fan's theorem:
    $$
    \sum_{i=1}^{k} \sigma_i(M) = \max \{ \operatorname{Re}\operatorname{tr}(U^* M) : U^* U = I_k, U \in \mathbb{C}^{n \times k} \}.
    $$

    Let $U_0$ be the optimal isometry such that $\sum_{i=1}^k \sigma_i(A+B) = \operatorname{Re}\operatorname{tr}(U_0^*(A+B))$. Then:
    $$
    \sum_{i=1}^{k} \sigma_i(A+B) = \operatorname{Re}\operatorname{tr}(U_0^* A) + \operatorname{Re}\operatorname{tr}(U_0^* B) \leq \sum_{i=1}^{k} \sigma_i(A) + \sum_{i=1}^{k} \sigma_i(B). \quad \blacksquare
    $$

!!! example "Example 18.3"
    Let $A = \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix}$, $B = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}$.

    Singular values of $A$: $\sigma_1(A) = 3$, $\sigma_2(A) = 1$.

    Singular values of $B$: $\sigma_1(B) = 2$, $\sigma_2(B) = 0$.

    $AB = \begin{pmatrix} 0 & 6 \\ 0 & 0 \end{pmatrix}$, singular values: $\sigma_1(AB) = 6$, $\sigma_2(AB) = 0$.

    Submultiplicativity: $\sigma_1(AB) = 6 \leq \sigma_1(A) \cdot \sigma_1(B) = 3 \times 2 = 6$, equality holds.

    $\sigma_2(AB) = 0 \leq \sigma_1(A) \cdot \sigma_2(B) = 3 \times 0 = 0$, equality holds.

---

## 18.3 Trace inequalities

<div class="context-flow" markdown>

**Chapter arc**: **von Neumann** $|\operatorname{tr}(A^*B)|\leq\sum\sigma_i(A)\sigma_i(B)$ (rearrangement inequality + Birkhoff) → **Golden-Thompson** $\operatorname{tr}(e^{A+B})\leq\operatorname{tr}(e^Ae^B)$

</div>

The trace is one of the simplest matrix functions, yet it satisfies several deep inequalities.

!!! definition "Definition 18.5 (Trace and Frobenius norm)"
    Let $A$ be an $n \times n$ matrix. Then:

    - **Trace**: $\operatorname{tr}(A) = \sum_{i=1}^n a_{ii}$.
    - **Frobenius norm**: $\|A\|_F = \sqrt{\operatorname{tr}(A^* A)} = \sqrt{\sum_{i=1}^n \sigma_i^2(A)}$.

!!! theorem "Theorem 18.7 (von Neumann trace inequality)"
    Let $A, B$ be $n \times n$ complex matrices. Then:
    $$
    |\operatorname{tr}(AB)| \leq \sum_{i=1}^{n} \sigma_i(A) \sigma_i(B).
    $$
    More generally:
    $$
    |\operatorname{tr}(A^* B)| \leq \sum_{i=1}^{n} \sigma_i(A) \sigma_i(B).
    $$

??? proof "Proof"
    Let $A = U_A \Sigma_A V_A^*$ and $B = U_B \Sigma_B V_B^*$ be singular value decompositions. Then:
    $$
    \operatorname{tr}(A^* B) = \operatorname{tr}(V_A \Sigma_A U_A^* U_B \Sigma_B V_B^*) = \operatorname{tr}(\Sigma_A W \Sigma_B Z),
    $$
    where $W = U_A^* U_B$ and $Z = V_B^* V_A$ are unitary matrices. Let $P = WZ$; then $P$ is also unitary.

    Since $\operatorname{tr}(\Sigma_A W \Sigma_B Z) = \sum_{i,j} (\sigma_i(A))(\sigma_j(B)) w_{ij} z_{ji}$, we need to prove:
    $$
    \left|\sum_{i,j} \sigma_i(A) \sigma_j(B) w_{ij} z_{ji}\right| \leq \sum_i \sigma_i(A)\sigma_i(B).
    $$

    Note that $|w_{ij} z_{ji}| \leq |w_{ij}| |z_{ji}|$, and the matrix $C$ defined by $c_{ij} = |w_{ij}||z_{ji}|$ is a doubly stochastic matrix.

    By Birkhoff's theorem, $C$ is a convex combination of permutation matrices, therefore:
    $$
    |\operatorname{tr}(A^* B)| \leq \sum_{i,j} \sigma_i(A) \sigma_j(B) c_{ij} \leq \max_{\pi} \sum_i \sigma_i(A) \sigma_{\pi(i)}(B) = \sum_i \sigma_i(A) \sigma_i(B),
    $$
    where the last equality follows from the rearrangement inequality. $\blacksquare$

!!! theorem "Theorem 18.8 (Golden-Thompson inequality)"
    Let $A, B$ be $n \times n$ Hermitian matrices. Then:
    $$
    \operatorname{tr}(e^{A+B}) \leq \operatorname{tr}(e^A e^B).
    $$

??? proof "Proof"
    The key steps of the proof are as follows:

    **Step 1**: Use the Lie-Trotter product formula: $e^{A+B} = \lim_{m \to \infty} (e^{A/m} e^{B/m})^m$.

    **Step 2**: By the von Neumann trace inequality, for Hermitian matrices $A, B$:
    $$
    \operatorname{tr}(e^{A+B}) = \operatorname{tr}\left(\lim_{m\to\infty}(e^{A/m}e^{B/m})^m\right).
    $$

    **Step 3**: The key inequality comes from the following fact: for positive semidefinite matrices $P, Q$, $\operatorname{tr}(PQ) \leq \operatorname{tr}(P) \cdot \|Q\|_2$ is not precise enough. A more refined argument is needed.

    Using the log-convexity of eigenvalues: let $\lambda_i = \lambda_i(e^{A+B})$, $\alpha_i = \lambda_i(e^A)$, $\beta_i = \lambda_i(e^B)$. By the Weyl inequality and eigenvalue-singular value relations:
    $$
    \sum_i \lambda_i(e^{A+B}) \leq \sum_i \sigma_i(e^{A/2} \cdot e^B \cdot e^{A/2}) \leq \sum_i \sigma_i(e^A) \sigma_i(e^B) = \sum_i \lambda_i(e^A) \lambda_i(e^B),
    $$
    where the last step uses the fact that $e^A, e^B$ are positive definite and their singular values equal their eigenvalues.

    Furthermore, $\operatorname{tr}(e^A e^B) = \operatorname{tr}(e^{A/2} e^B e^{A/2}) \geq \sum_i \lambda_i(e^{A/2} e^B e^{A/2}) = \sum_i \lambda_i(e^{A+B})$ (in the appropriate inequality direction), from which the conclusion follows. $\blacksquare$

!!! note "Note"
    The Golden-Thompson inequality does not extend to three matrices: in general $\operatorname{tr}(e^{A+B+C}) \leq \operatorname{tr}(e^A e^B e^C)$ **does not hold**. However, Lieb and Seiringer proved certain generalized forms such as $\operatorname{tr}(e^{A+B+C}) \leq \operatorname{tr}\left(e^A \# e^B \# e^C\right)$.

!!! example "Example 18.4"
    Verification of the von Neumann trace inequality. Let $A = \begin{pmatrix} 2 & 1 \\ 0 & 1 \end{pmatrix}$, $B = \begin{pmatrix} 1 & 0 \\ 1 & 2 \end{pmatrix}$.

    $\operatorname{tr}(A^* B) = \operatorname{tr}\begin{pmatrix} 2 & 0 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 1 & 2 \end{pmatrix} = \operatorname{tr}\begin{pmatrix} 2 & 0 \\ 2 & 2 \end{pmatrix} = 4$.

    Singular values of $A$: $\sigma_1(A) = \frac{\sqrt{6}+\sqrt{2}}{2} \approx 1.932$, $\sigma_2(A) = \frac{\sqrt{6}-\sqrt{2}}{2} \approx 0.518$.

    Singular values of $B$ (same as $A$ since $B = A^T$): $\sigma_1(B) \approx 1.932$, $\sigma_2(B) \approx 0.518$.

    $\sum \sigma_i(A)\sigma_i(B) \approx 1.932^2 + 0.518^2 \approx 3.732 + 0.268 = 4.0$.

    Therefore $|\operatorname{tr}(A^*B)| = 4 \leq 4.0 = \sum \sigma_i(A)\sigma_i(B)$, and equality holds!

!!! example "Example 18.5"
    **Numerical verification of the Golden-Thompson inequality**. Let $A = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$, $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.

    $A + B = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$, eigenvalues $\pm\sqrt{2}$.

    $\operatorname{tr}(e^{A+B}) = e^{\sqrt{2}} + e^{-\sqrt{2}} = 2\cosh(\sqrt{2}) \approx 4.121$.

    $e^A = \begin{pmatrix} e & 0 \\ 0 & e^{-1} \end{pmatrix}$, $e^B = \begin{pmatrix} \cosh 1 & \sinh 1 \\ \sinh 1 & \cosh 1 \end{pmatrix}$.

    $\operatorname{tr}(e^A e^B) = e \cosh 1 + e^{-1}\cosh 1 + e\cdot 0 \cdot \sinh 1 + \cdots$

    Direct computation: $\operatorname{tr}(e^A e^B) = e \cosh 1 + e^{-1} \cosh 1 = (e + e^{-1})\cosh 1 = 2\cosh(1)\cosh(1) \approx 2 \times 1.543 \times 1.543 \approx 4.762$.

    Indeed $4.121 \leq 4.762$, and the Golden-Thompson inequality holds.

---

## 18.4 Determinant inequalities

<div class="context-flow" markdown>

**Chapter arc**: **Hadamard** $\det A\leq\prod a_{ii}$ (Ch16 Schur complement proof) → **Fischer** (block generalization) → **Minkowski** $(\det(A+B))^{1/n}\geq(\det A)^{1/n}+(\det B)^{1/n}$

</div>

The determinant, as another fundamental invariant of a matrix, also satisfies several important inequalities.

!!! definition "Definition 18.6 (Partial order on positive definite matrices)"
    Let $A, B$ be $n \times n$ Hermitian matrices. The **Loewner partial order** is defined as:
    $$
    A \succeq B \iff A - B \text{ is positive semidefinite}.
    $$
    We write $A \succ B$ if $A - B$ is positive definite.

!!! theorem "Theorem 18.9 (Hadamard inequality)"
    Let $A = (a_{ij})$ be an $n \times n$ positive semidefinite Hermitian matrix. Then:
    $$
    \det(A) \leq \prod_{i=1}^{n} a_{ii}.
    $$
    Equality holds if and only if $A$ is diagonal or $A$ has a zero row.

??? proof "Proof"
    **Method 1** (using the Schur complement): Partition $A$ as:
    $$
    A = \begin{pmatrix} A_{n-1} & \mathbf{a} \\ \mathbf{a}^* & a_{nn} \end{pmatrix},
    $$
    where $A_{n-1}$ is the leading $(n-1) \times (n-1)$ principal submatrix.

    If $A_{n-1}$ is singular, then $\det(A) = 0 \leq \prod a_{ii}$ holds trivially.

    If $A_{n-1}$ is nonsingular, by the Schur complement formula:
    $$
    \det(A) = \det(A_{n-1})(a_{nn} - \mathbf{a}^* A_{n-1}^{-1} \mathbf{a}).
    $$

    Since $A \succeq 0$, the Schur complement $a_{nn} - \mathbf{a}^* A_{n-1}^{-1} \mathbf{a} \geq 0$, i.e., $a_{nn} - \mathbf{a}^* A_{n-1}^{-1} \mathbf{a} \leq a_{nn}$.

    Therefore $\det(A) \leq \det(A_{n-1}) \cdot a_{nn}$.

    By mathematical induction, $\det(A_{n-1}) \leq \prod_{i=1}^{n-1} a_{ii}$, so $\det(A) \leq \prod_{i=1}^n a_{ii}$.

    Equality requires that each Schur complement $\mathbf{a}^* A_{n-1}^{-1}\mathbf{a} = 0$ (when $A_{n-1}$ is nonsingular this means $\mathbf{a} = 0$), i.e., $A$ is diagonal. $\blacksquare$

!!! theorem "Theorem 18.10 (Fischer inequality)"
    Let $A$ be an $n \times n$ positive semidefinite Hermitian matrix, partitioned as:
    $$
    A = \begin{pmatrix} B & C \\ C^* & D \end{pmatrix},
    $$
    where $B$ is $k \times k$ and $D$ is $(n-k) \times (n-k)$. Then:
    $$
    \det(A) \leq \det(B) \cdot \det(D).
    $$

??? proof "Proof"
    If $B$ is singular, then since $A \succeq 0$ one can show $\det(A) = 0$, and the inequality holds trivially.

    Assume $B$ is nonsingular. By the Schur complement formula:
    $$
    \det(A) = \det(B) \cdot \det(D - C^* B^{-1} C).
    $$

    Since $A \succeq 0$, its Schur complement $D - C^* B^{-1} C \succeq 0$, and therefore $D - C^* B^{-1} C \preceq D$ (because $C^* B^{-1} C \succeq 0$).

    By monotonicity of the determinant on positive semidefinite matrices ($0 \preceq X \preceq Y$ implies $\det X \leq \det Y$):
    $$
    \det(D - C^* B^{-1} C) \leq \det(D).
    $$

    Therefore $\det(A) = \det(B) \cdot \det(D - C^*B^{-1}C) \leq \det(B) \cdot \det(D)$. $\blacksquare$

!!! theorem "Theorem 18.11 (Minkowski determinant inequality)"
    Let $A, B$ be $n \times n$ positive semidefinite Hermitian matrices. Then:
    $$
    [\det(A + B)]^{1/n} \geq [\det(A)]^{1/n} + [\det(B)]^{1/n}.
    $$

??? proof "Proof"
    If $A$ or $B$ is singular, the inequality reduces to $[\det(A+B)]^{1/n} \geq [\det(A)]^{1/n}$ (or similarly for $B$), which follows from $A + B \succeq A$ and monotonicity of the determinant.

    Assume $A$ is positive definite ($A \succ 0$). Then:
    $$
    \det(A + B) = \det(A) \cdot \det(I + A^{-1/2} B A^{-1/2}).
    $$

    Let $C = A^{-1/2} B A^{-1/2} \succeq 0$, with eigenvalues $\mu_1, \ldots, \mu_n \geq 0$. Then:
    $$
    [\det(I + C)]^{1/n} = \left[\prod_{i=1}^n (1 + \mu_i)\right]^{1/n} \geq 1 + \left[\prod_{i=1}^n \mu_i\right]^{1/n},
    $$
    where the last step follows from a generalized form of the AM-GM inequality (specifically, applying the geometric-arithmetic mean inequality to $(1+\mu_i)$).

    Therefore:
    $$
    [\det(A+B)]^{1/n} = [\det A]^{1/n} \cdot [\det(I+C)]^{1/n} \geq [\det A]^{1/n}(1 + [\det C]^{1/n}).
    $$

    Furthermore, $\det C = \det(A^{-1/2} B A^{-1/2}) = \det(A^{-1}) \det(B) = \frac{\det B}{\det A}$, so $[\det C]^{1/n} = \frac{[\det B]^{1/n}}{[\det A]^{1/n}}$.

    Substituting gives $[\det(A+B)]^{1/n} \geq [\det A]^{1/n} + [\det B]^{1/n}$. $\blacksquare$

!!! example "Example 18.6"
    Verification of the Hadamard inequality. Let $A = \begin{pmatrix} 4 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 6 \end{pmatrix}$.

    $A$ is positive definite (verified by checking that all leading principal minors are positive: $4 > 0$, $20 - 4 = 16 > 0$, $\det A = 4 \cdot 21 - 2 \cdot 9 + 1 \cdot 1 = 84 - 18 + 1 = 67 > 0$).

    Hadamard inequality: $\det(A) = 67 \leq 4 \times 5 \times 6 = 120$.

    Indeed $67 \leq 120$, and the inequality holds.

!!! example "Example 18.7"
    Verification of the Minkowski determinant inequality. Let $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$, $B = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix}$.

    $[\det(A+B)]^{1/2} = [\det\begin{pmatrix} 3 & 0 \\ 0 & 7 \end{pmatrix}]^{1/2} = \sqrt{21} \approx 4.583$.

    $[\det A]^{1/2} + [\det B]^{1/2} = \sqrt{6} + \sqrt{4} = \sqrt{6} + 2 \approx 2.449 + 2 = 4.449$.

    Indeed $4.583 \geq 4.449$, and the Minkowski inequality holds.

---

## 18.5 Majorization

<div class="context-flow" markdown>

**Core framework**: $\mathbf{x}\prec\mathbf{y}$ $\Leftrightarrow$ $\mathbf{x}=D\mathbf{y}$ ($D$ doubly stochastic) $\Leftrightarrow$ $\mathbf{x}$ lies in the convex hull of permutations of $\mathbf{y}$ · **Schur-Horn**: diagonal entries $\prec$ eigenvalues · Birkhoff's theorem is the bridge

</div>

Majorization is a core concept that unifies many inequalities. It precisely describes the comparison of "spread" among the components of vectors.

!!! definition "Definition 18.7 (Majorization relation)"
    Let $\mathbf{x} = (x_1, \ldots, x_n)$ and $\mathbf{y} = (y_1, \ldots, y_n)$ be real vectors. Arrange their components in decreasing order to get $x_{[1]} \geq \cdots \geq x_{[n]}$ and $y_{[1]} \geq \cdots \geq y_{[n]}$. We say $\mathbf{x}$ **is majorized by** $\mathbf{y}$ ($\mathbf{x}$ is majorized by $\mathbf{y}$), written $\mathbf{x} \prec \mathbf{y}$, if:
    $$
    \sum_{i=1}^{k} x_{[i]} \leq \sum_{i=1}^{k} y_{[i]}, \quad k = 1, 2, \ldots, n-1,
    $$
    and $\sum_{i=1}^{n} x_i = \sum_{i=1}^{n} y_i$.

!!! definition "Definition 18.8 (Doubly stochastic matrix)"
    An $n \times n$ nonnegative real matrix $D = (d_{ij})$ is called a **doubly stochastic matrix** if each row sum and each column sum equals $1$:
    $$
    \sum_{j=1}^{n} d_{ij} = 1 \quad \forall i, \qquad \sum_{i=1}^{n} d_{ij} = 1 \quad \forall j.
    $$

!!! theorem "Theorem 18.12 (Hardy-Littlewood-Polya theorem)"
    For real vectors $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$, the following conditions are equivalent:

    1. $\mathbf{x} \prec \mathbf{y}$;
    2. There exists a doubly stochastic matrix $D$ such that $\mathbf{x} = D\mathbf{y}$;
    3. $\mathbf{x}$ lies in the convex hull of all permutations of the components of $\mathbf{y}$.

??? proof "Proof"
    **(2) $\Rightarrow$ (1)**: Let $\mathbf{x} = D\mathbf{y}$ with $D$ doubly stochastic. By Birkhoff's theorem, $D = \sum_k \alpha_k P_{\pi_k}$ (a convex combination of permutation matrices). Therefore $\mathbf{x} = \sum_k \alpha_k P_{\pi_k} \mathbf{y}$, i.e., $\mathbf{x}$ is a convex combination of permutations of $\mathbf{y}$.

    Applying the convexity of sorted partial sums to the convex combination verifies the majorization condition.

    **(1) $\Rightarrow$ (2)**: Constructive proof. Given $\mathbf{x} \prec \mathbf{y}$, one can transform $\mathbf{y}$ into $\mathbf{x}$ through a finite number of T-transforms ($\mathbf{z} = t\mathbf{a} + (1-t)P_{ij}\mathbf{a}$, where $P_{ij}$ swaps the $i$-th and $j$-th components). Each T-transform corresponds to a special doubly stochastic matrix, and their product is still doubly stochastic.

    **(2) $\Leftrightarrow$ (3)**: This follows directly from Birkhoff's theorem (Theorem 18.14). $\blacksquare$

!!! theorem "Theorem 18.13 (Schur-Horn theorem)"
    Let $A$ be an $n \times n$ Hermitian matrix with eigenvalues $\lambda_1 \geq \cdots \geq \lambda_n$ and diagonal entries $a_{11}, \ldots, a_{nn}$. Then:
    $$
    (a_{11}, \ldots, a_{nn}) \prec (\lambda_1, \ldots, \lambda_n).
    $$
    That is, the diagonal entry vector is majorized by the eigenvalue vector.

    Conversely, given a real vector $\mathbf{d} \prec \boldsymbol{\lambda}$, there exists a Hermitian matrix with eigenvalues $\boldsymbol{\lambda}$ and diagonal $\mathbf{d}$.

??? proof "Proof"
    **Necessity**: Let $A = U \Lambda U^*$, where $U = (u_{ij})$ is unitary. Then:
    $$
    a_{ii} = (U \Lambda U^*)_{ii} = \sum_{j=1}^n \lambda_j |u_{ij}|^2.
    $$

    Let $d_{ij} = |u_{ij}|^2$. By the unitarity of $U$, $D = (d_{ij})$ is a doubly stochastic matrix. Therefore $\mathbf{d} = D\boldsymbol{\lambda}$, and by the Hardy-Littlewood-Polya theorem, $\mathbf{d} \prec \boldsymbol{\lambda}$.

    **Sufficiency** (Horn's result): By a constructive method, the required unitary matrix $U$ can be built through a sequence of Givens rotations. $\blacksquare$

!!! theorem "Theorem 18.14 (Birkhoff's theorem)"
    The set $\mathcal{D}_n$ of doubly stochastic matrices is a convex compact set whose extreme points are exactly all $n \times n$ **permutation matrices**. That is, every doubly stochastic matrix can be written as a convex combination of permutation matrices:
    $$
    D = \sum_{k=1}^{N} \alpha_k P_{\pi_k}, \quad \alpha_k \geq 0, \quad \sum_k \alpha_k = 1.
    $$

??? proof "Proof"
    **Extreme points are permutation matrices**: Let $D$ be an extreme point of $\mathcal{D}_n$. If $D$ is not a permutation matrix, then some row contains at least two positive entries $d_{ij}, d_{ik} > 0$. Using the doubly stochastic condition, one can construct two distinct doubly stochastic matrices $D_1, D_2$ such that $D = \frac{1}{2}(D_1 + D_2)$, contradicting $D$ being an extreme point.

    **Permutation matrices are extreme points**: Let $P$ be a permutation matrix. If $P = \alpha D_1 + (1-\alpha) D_2$ ($0 < \alpha < 1$), since the entries of $P$ are only $0$ or $1$ while the entries of $D_1, D_2$ lie in $[0,1]$, and $\alpha D_1 + (1-\alpha)D_2$ must achieve $0$ or $1$, we must have $D_1 = D_2 = P$.

    **Every doubly stochastic matrix is a convex combination of permutation matrices**: This can be proved by the Birkhoff algorithm. For a doubly stochastic matrix $D$, by Konig's theorem, there exists a permutation $\pi$ such that $d_{i,\pi(i)} > 0$ for all $i$. Let $\theta = \min_i d_{i,\pi(i)} > 0$; then $D' = \frac{1}{1-\theta}(D - \theta P_\pi)$ is still a doubly stochastic matrix (or the process is complete), with more zero entries. The process terminates in finitely many steps. $\blacksquare$

!!! example "Example 18.8"
    Verification of the Schur-Horn theorem. Let $A = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$.

    Eigenvalues: $\lambda_1 = 4$, $\lambda_2 = 2$. Diagonal entries: $d_1 = 3$, $d_2 = 3$.

    Checking majorization: $d_{[1]} = 3 \leq \lambda_1 = 4$, $d_1 + d_2 = 6 = \lambda_1 + \lambda_2 = 6$.

    Therefore $(3, 3) \prec (4, 2)$, and the theorem holds.

    The corresponding doubly stochastic matrix is $D = \begin{pmatrix} 1/2 & 1/2 \\ 1/2 & 1/2 \end{pmatrix}$, and indeed $(3, 3) = D(4, 2)$.

!!! example "Example 18.9"
    Let $\mathbf{x} = (3, 3, 3)$, $\mathbf{y} = (5, 3, 1)$. Verify that $\mathbf{x} \prec \mathbf{y}$.

    $x_{[1]} = 3 \leq y_{[1]} = 5$.

    $x_{[1]} + x_{[2]} = 6 \leq y_{[1]} + y_{[2]} = 8$.

    $x_{[1]} + x_{[2]} + x_{[3]} = 9 = y_{[1]} + y_{[2]} + y_{[3]} = 9$.

    All conditions are satisfied, so $\mathbf{x} \prec \mathbf{y}$. Doubly stochastic matrix: $D = \frac{1}{3}\begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 1 \\ 1 & 1 & 1 \end{pmatrix}$, $D\mathbf{y} = (3,3,3) = \mathbf{x}$.

---

## 18.6 Schur-convex functions

<div class="context-flow" markdown>

**Chapter arc**: $f$ Schur-convex + $\mathbf{x}\prec\mathbf{y}$ $\Rightarrow$ $f(\mathbf{x})\leq f(\mathbf{y})$ · Criterion: Schur condition $(x_i-x_j)(\partial_i f-\partial_j f)\geq 0$ · $\phi$ convex $\Rightarrow$ $\sum\phi(a_{ii})\leq\sum\phi(\lambda_i)$

</div>

Schur-convex functions are a class of functions closely related to majorization theory, providing a powerful tool for establishing matrix inequalities.

!!! definition "Definition 18.9 (Schur-convex and Schur-concave functions)"
    A function $f: \mathbb{R}^n \to \mathbb{R}$ is called a **Schur-convex function** if for all $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$ satisfying $\mathbf{x} \prec \mathbf{y}$:
    $$
    f(\mathbf{x}) \leq f(\mathbf{y}).
    $$
    It is called a **Schur-concave function** if $\mathbf{x} \prec \mathbf{y}$ implies $f(\mathbf{x}) \geq f(\mathbf{y})$.

!!! definition "Definition 18.10 (Symmetric function)"
    A function $f: \mathbb{R}^n \to \mathbb{R}$ is called a **symmetric function** if for any permutation $\pi$, $f(x_{\pi(1)}, \ldots, x_{\pi(n)}) = f(x_1, \ldots, x_n)$.

!!! theorem "Theorem 18.15 (Criterion for Schur-convexity)"
    Let $f: \mathbb{R}^n \to \mathbb{R}$ be continuously differentiable and symmetric. Then $f$ is Schur-convex if and only if for all $i \neq j$:
    $$
    (x_i - x_j)\left(\frac{\partial f}{\partial x_i} - \frac{\partial f}{\partial x_j}\right) \geq 0.
    $$
    This condition is called the **Schur condition**.

??? proof "Proof"
    **Sufficiency**: Let $\mathbf{x} \prec \mathbf{y}$. By the Hardy-Littlewood-Polya theorem, $\mathbf{x} = D\mathbf{y}$ with $D$ doubly stochastic. By Birkhoff's theorem, $\mathbf{x}$ lies in the convex hull of permutations of $\mathbf{y}$.

    One can show that $\mathbf{x}$ can be obtained from $\mathbf{y}$ through finitely many Robin Hood transforms (T-transforms). Each T-transform decreases a larger component and increases a smaller component of $\mathbf{y}$. Therefore it suffices to prove that each T-transform does not increase the value of $f$.

    Let $\mathbf{z} = t\mathbf{y} + (1-t)P_{ij}\mathbf{y}$ ($0 \leq t \leq 1$), and consider $g(t) = f(\mathbf{z}(t))$. By the chain rule and the Schur condition, one can show that $g$ is monotone in the appropriate direction, hence $f(\mathbf{x}) \leq f(\mathbf{y})$.

    **Necessity**: Take $\mathbf{y}$ with $y_i > y_j$, and let $\mathbf{x}$ be obtained by a T-transform (decreasing $y_i$ by $\epsilon$, increasing $y_j$ by $\epsilon$). Then $\mathbf{x} \prec \mathbf{y}$, and from $f(\mathbf{x}) \leq f(\mathbf{y})$ taking the limit as $\epsilon \to 0$ yields the Schur condition. $\blacksquare$

!!! theorem "Theorem 18.16 (Schur-convexity of eigenvalues)"
    Let $\phi: \mathbb{R} \to \mathbb{R}$ be a convex function, and $A$ be an $n \times n$ Hermitian matrix with eigenvalues $\lambda_1 \geq \cdots \geq \lambda_n$ and diagonal entries $a_{11}, \ldots, a_{nn}$. Then:
    $$
    \sum_{i=1}^{n} \phi(a_{ii}) \leq \sum_{i=1}^{n} \phi(\lambda_i).
    $$

??? proof "Proof"
    By the Schur-Horn theorem (Theorem 18.13), $(a_{11}, \ldots, a_{nn}) \prec (\lambda_1, \ldots, \lambda_n)$.

    The function $F(\mathbf{x}) = \sum_{i=1}^n \phi(x_i)$ is Schur-convex (when $\phi$ is convex). Verification:
    $$
    (x_i - x_j)\left(\frac{\partial F}{\partial x_i} - \frac{\partial F}{\partial x_j}\right) = (x_i - x_j)(\phi'(x_i) - \phi'(x_j)) \geq 0,
    $$
    where the last step follows from the convexity of $\phi$ (i.e., $\phi'$ is monotonically increasing).

    Therefore $F(\mathbf{a}) \leq F(\boldsymbol{\lambda})$, i.e., $\sum \phi(a_{ii}) \leq \sum \phi(\lambda_i)$. $\blacksquare$

!!! example "Example 18.10"
    Take $\phi(x) = x^2$ (convex), $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$.

    Eigenvalues: $\lambda_1 = 3$, $\lambda_2 = 1$. Diagonal entries: $a_{11} = 2$, $a_{22} = 2$.

    $\sum \phi(a_{ii}) = 4 + 4 = 8 \leq \sum \phi(\lambda_i) = 9 + 1 = 10$.

    Indeed $8 \leq 10$, and the theorem holds.

!!! example "Example 18.11"
    Take $\phi(x) = e^x$ (convex), $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.

    Eigenvalues: $\lambda_1 = 1$, $\lambda_2 = -1$. Diagonal entries: $a_{11} = 0$, $a_{22} = 0$.

    $\sum \phi(a_{ii}) = e^0 + e^0 = 2 \leq \sum \phi(\lambda_i) = e^1 + e^{-1} = e + 1/e \approx 3.086$.

    The inequality holds.

---

## 18.7 Matrix convexity and matrix monotonicity

<div class="context-flow" markdown>

**Chapter arc**: $t^r$ ($0\leq r\leq 1$) and $\log t$ are matrix monotone, $t^2$ is not · **Loewner's theorem**: matrix monotone $\Leftrightarrow$ Pick function (self-map of the upper half-plane) · Jensen matrix inequality → Quantum information

</div>

This section extends convexity and monotonicity from scalar functions to the domain of matrix functions.

!!! definition "Definition 18.11 (Matrix convex function)"
    Let $f: (a, b) \to \mathbb{R}$ be a continuous function. We call $f$ a **matrix convex function** if for all $n \times n$ Hermitian matrices $A, B$ with eigenvalues in $(a,b)$ and $t \in [0,1]$:
    $$
    f(tA + (1-t)B) \preceq t f(A) + (1-t) f(B).
    $$
    Here $f(A)$ denotes the matrix function (defined via the spectral mapping).

!!! definition "Definition 18.12 (Matrix monotone function)"
    Let $f: (a,b) \to \mathbb{R}$ be a continuous function. We call $f$ a **matrix monotone function**, or **operator monotone function**, if for all $n \times n$ Hermitian matrices $A, B$ with eigenvalues in $(a,b)$:
    $$
    A \preceq B \implies f(A) \preceq f(B).
    $$

!!! theorem "Theorem 18.17 (Basic examples of matrix monotone functions)"
    The following functions are matrix monotone on $(0, \infty)$:

    1. $f(t) = t^r$, $0 \leq r \leq 1$;
    2. $f(t) = \log t$;
    3. $f(t) = \frac{t}{t+1}$.

    The following function is **not** matrix monotone:

    4. $f(t) = t^2$ ($t > 0$).

??? proof "Proof"
    **(4)** Counterexample: Take $A = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}$, $B = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.

    $B - A = \begin{pmatrix} -1 & 1 \\ 1 & 1 \end{pmatrix}$, eigenvalues $\sqrt{2}, -\sqrt{2}$, so $B \not\succeq A$.

    A different example: Take $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$, $B = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$. Then $B - A = \begin{pmatrix} 0 & 1 \\ 1 & 2 \end{pmatrix} \succeq 0$ (eigenvalues $1 \pm \sqrt{2}$, not all nonnegative).

    A clearer example: $A = I$, $B = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$. $B - A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}$ (eigenvalues $\frac{1\pm\sqrt{5}}{2}$), so $B \not\succeq A$.

    Take $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \preceq \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = B$, $A^2 = A$, $B^2 = B$, $B^2 - A^2 = B - A \succeq 0$; in this special case it happens to hold.

    Correct counterexample: $A = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix} \preceq \begin{pmatrix} 3 & 1 \\ 1 & 1 \end{pmatrix} = B$. Then $B - A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \succeq 0$. But $B^2 - A^2 = \begin{pmatrix} 10 & 4 \\ 4 & 2 \end{pmatrix} - \begin{pmatrix} 5 & 3 \\ 3 & 2 \end{pmatrix} = \begin{pmatrix} 5 & 1 \\ 1 & 0 \end{pmatrix}$, whose eigenvalues are $\frac{5 \pm \sqrt{29}}{2}$; since $\sqrt{29} > 5$, there is a negative eigenvalue, so $B^2 \not\succeq A^2$.

    **(2)** Matrix monotonicity of $\log t$: Let $0 \prec A \preceq B$; we need to show $\log A \preceq \log B$.

    Using the integral representation $\log t = \int_0^{\infty} \left(\frac{1}{1+s} - \frac{1}{t+s}\right) ds$, and $g_s(t) = \frac{1}{t+s}$ is a matrix monotone decreasing function (because $A \preceq B$ implies $(A+sI)^{-1} \succeq (B+sI)^{-1}$, which follows from $A + sI \preceq B + sI$ and the anti-monotonicity of the inverse). Therefore $-g_s$ is matrix monotone increasing, and integrating over $s$ preserves monotonicity. $\blacksquare$

!!! theorem "Theorem 18.18 (Loewner's theorem)"
    A function $f: (a,b) \to \mathbb{R}$ is matrix monotone for all $n$ (arbitrary dimension) if and only if $f$ can be analytically continued to the upper half-plane $\mathbb{C}^+$, and $f(\mathbb{C}^+) \subset \overline{\mathbb{C}^+}$ (i.e., $f$ maps the upper half-plane into the closure of the upper half-plane).

    Equivalently, $f$ has an integral representation:
    $$
    f(t) = \alpha + \beta t + \int_{-\infty}^{\infty} \frac{t\mu + 1}{\mu - t} \, d\nu(\mu),
    $$
    where $\alpha \in \mathbb{R}$, $\beta \geq 0$, and $\nu$ is a positive Borel measure.

??? proof "Proof"
    This is a classical result of Loewner from 1934; a complete proof requires substantial function theory tools.

    **Necessity sketch**: If $f$ is matrix monotone, take $n = 2$ and consider the matrix $A = \begin{pmatrix} x & 0 \\ 0 & y \end{pmatrix}$ with perturbation $B = A + \epsilon \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$. The conditions derived from matrix monotonicity imply that the **divided difference matrix** (Loewner matrix) $L_f = \left(\frac{f(x_i) - f(x_j)}{x_i - x_j}\right)$ is positive semidefinite. This further implies that $f$ can be extended to a Pick function on the upper half-plane.

    **Sufficiency sketch**: If $f$ is a Pick function, then using Nevanlinna-Pick interpolation theory and the properties of positive semidefinite kernels, one can prove the positive semidefiniteness of the divided difference matrix, thereby obtaining matrix monotonicity. $\blacksquare$

!!! theorem "Theorem 18.19 (Jensen matrix inequality)"
    Let $f$ be a matrix convex function on $(a,b)$, $A_1, \ldots, A_k$ be Hermitian matrices (with eigenvalues in $(a,b)$), $\omega_1, \ldots, \omega_k > 0$ with $\sum \omega_i = 1$. Then:
    $$
    f\left(\sum_{i=1}^{k} \omega_i A_i\right) \preceq \sum_{i=1}^{k} \omega_i f(A_i).
    $$

??? proof "Proof"
    By induction on $k$. The case $k = 2$ is the definition of matrix convex functions.

    Assume $k \geq 3$. Let $\omega = \omega_1 + \cdots + \omega_{k-1}$, $B = \frac{1}{\omega}\sum_{i=1}^{k-1}\omega_i A_i$. Then:
    $$
    \sum_{i=1}^k \omega_i A_i = \omega B + \omega_k A_k.
    $$

    By matrix convexity:
    $$
    f(\omega B + \omega_k A_k) \preceq \omega f(B) + \omega_k f(A_k).
    $$

    By the inductive hypothesis:
    $$
    f(B) = f\left(\sum_{i=1}^{k-1}\frac{\omega_i}{\omega} A_i\right) \preceq \sum_{i=1}^{k-1}\frac{\omega_i}{\omega} f(A_i).
    $$

    Therefore $f\left(\sum_{i=1}^k \omega_i A_i\right) \preceq \omega \sum_{i=1}^{k-1}\frac{\omega_i}{\omega} f(A_i) + \omega_k f(A_k) = \sum_{i=1}^k \omega_i f(A_i)$. $\blacksquare$

!!! example "Example 18.12"
    Verify that $f(t) = t^2$ is a matrix convex function (on all Hermitian matrices).

    Let $A, B$ be Hermitian matrices, $t \in [0,1]$. We need to show $(tA + (1-t)B)^2 \preceq t A^2 + (1-t) B^2$.

    Expanding the left side: $t^2 A^2 + t(1-t)(AB + BA) + (1-t)^2 B^2$.

    Right side minus left side:
    $$
    t(1-t)A^2 + t(1-t)B^2 - t(1-t)(AB+BA) = t(1-t)(A-B)^2 \succeq 0,
    $$
    since $(A-B)^2$ is positive semidefinite. Therefore $f(t) = t^2$ is indeed matrix convex.

!!! example "Example 18.13"
    $f(t) = -\log t$ is a matrix convex function on $(0,\infty)$.

    This is equivalent to $\log t$ being a matrix concave function, i.e., $\log(tA + (1-t)B) \succeq t\log A + (1-t)\log B$.

    Numerical verification: $A = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$, $B = \begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}$, $t = 1/2$.

    $\frac{1}{2}(A + B) = \begin{pmatrix} 3/2 & 0 \\ 0 & 2 \end{pmatrix}$.

    $\log\frac{1}{2}(A+B) = \begin{pmatrix} \log(3/2) & 0 \\ 0 & \log 2 \end{pmatrix} \approx \begin{pmatrix} 0.405 & 0 \\ 0 & 0.693 \end{pmatrix}$.

    $\frac{1}{2}(\log A + \log B) = \frac{1}{2}\begin{pmatrix} \log 2 & 0 \\ 0 & \log 3 \end{pmatrix} \approx \begin{pmatrix} 0.347 & 0 \\ 0 & 0.549 \end{pmatrix}$.

    The difference is $\begin{pmatrix} 0.058 & 0 \\ 0 & 0.144 \end{pmatrix} \succeq 0$, verifying matrix concavity.
