# Chapter 14  Matrix Analysis

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Norms (Ch15 preview)

**Chapter arc**: Topology of matrix spaces → **Spectral radius** governs convergence → **Gershgorin** locates eigenvalues → Matrix calculus → Eigenvalue–singular value bridge

**Further connections**：Spectral radius is crucial in dynamical system stability analysis, convergence of iterative methods (Jacobi, Gauss-Seidel), and population dynamics; Gershgorin's theorem provides quick eigenvalue bounds in engineering without exact computation

</div>

Matrix analysis is the intersection of linear algebra and analysis. It introduces analytic tools such as limits, continuity, and differentiation into matrix spaces, providing a powerful theoretical framework for studying asymptotic behavior, functional properties, and spectral structure of matrices. Starting from the topology of matrix spaces, this chapter discusses convergence of matrix series, characterizations of the spectral radius, eigenvalue localization (Gershgorin disc theorem), matrix calculus, and the deep connections between eigenvalues and singular values. These topics form the mathematical foundation of numerical linear algebra, control theory, and optimization.

---

## 14.1 Topology of matrix spaces

<div class="context-flow" markdown>

**Chapter arc**: Finite-dimensional $\Rightarrow$ norm equivalence $\Rightarrow$ entrywise convergence = norm convergence = complete (Banach) space

</div>

The matrix space $\mathbb{C}^{m \times n}$ (or $\mathbb{R}^{m \times n}$) can be naturally viewed as a finite-dimensional normed linear space. By introducing a norm on this space, we can discuss analytic concepts such as convergence of matrix sequences and continuity of matrix functions.

!!! definition "Definition 14.1 (Convergence of matrix sequences)"
    Let $\{A_k\}_{k=1}^{\infty}$ be a sequence of matrices in $\mathbb{C}^{m \times n}$, and $A \in \mathbb{C}^{m \times n}$. The sequence $\{A_k\}$ is said to **converge** to $A$, written $\lim_{k \to \infty} A_k = A$, if for some (and hence any) matrix norm $\|\cdot\|$ on $\mathbb{C}^{m \times n}$,

    $$\lim_{k \to \infty} \|A_k - A\| = 0.$$

!!! definition "Definition 14.2 (Entrywise convergence)"
    The matrix sequence $\{A_k\}$ **converges entrywise** to $A$ if for every $i, j$, $\lim_{k \to \infty} (A_k)_{ij} = A_{ij}$.

!!! theorem "Theorem 14.1 (Equivalence of convergence)"
    Let $\{A_k\} \subset \mathbb{C}^{m \times n}$ and $A \in \mathbb{C}^{m \times n}$. The following conditions are equivalent:

    (1) $\{A_k\}$ converges to $A$ in some matrix norm;

    (2) $\{A_k\}$ converges to $A$ in every matrix norm;

    (3) $\{A_k\}$ converges to $A$ entrywise.

??? proof "Proof"
    **(1) $\Rightarrow$ (2)**: Since $\mathbb{C}^{m \times n}$ is finite-dimensional, all norms on it are equivalent. Let $\|\cdot\|_\alpha$ and $\|\cdot\|_\beta$ be two norms; then there exist constants $c_1, c_2 > 0$ such that $c_1\|X\|_\alpha \leq \|X\|_\beta \leq c_2\|X\|_\alpha$ for all $X$. If $\|A_k - A\|_\alpha \to 0$, then $\|A_k - A\|_\beta \leq c_2\|A_k - A\|_\alpha \to 0$.

    **(2) $\Rightarrow$ (3)**: Take $\|\cdot\|$ to be the Frobenius norm; then $|(A_k)_{ij} - A_{ij}| \leq \|A_k - A\|_F \to 0$.

    **(3) $\Rightarrow$ (1)**: If entrywise convergence holds, take $\|\cdot\|_{\max} = \max_{i,j}|a_{ij}|$ (this is a norm); then $\|A_k - A\|_{\max} = \max_{i,j}|(A_k)_{ij} - A_{ij}| \to 0$. By norm equivalence, convergence holds in any norm. $\blacksquare$

!!! definition "Definition 14.3 (Normed space structure of matrix spaces)"
    $\mathbb{C}^{m \times n}$ equipped with a matrix norm $\|\cdot\|$ forms a **normed linear space**. Since $\mathbb{C}^{m \times n}$ is isomorphic to $\mathbb{C}^{mn}$ (as a vector space), it is finite-dimensional and hence **complete**, i.e., it is a **Banach space**.

!!! theorem "Theorem 14.2 (Compactness in matrix spaces)"
    A bounded closed subset of $\mathbb{C}^{m \times n}$ is compact. In particular, every bounded matrix sequence has a convergent subsequence.

??? proof "Proof"
    $\mathbb{C}^{m \times n}$ is isomorphic to $\mathbb{C}^{mn}$, which is a finite-dimensional normed space. By the Heine–Borel theorem, bounded closed subsets of finite-dimensional normed spaces are compact. A bounded sequence lies in some closed ball, and by compactness a convergent subsequence can be extracted. $\blacksquare$

!!! example "Example 14.1"
    Consider the matrix sequence $A_k = \begin{pmatrix} 1/k & 2/k \\ 0 & 1/k^2 \end{pmatrix}$.

    Entrywise, $(A_k)_{11} = 1/k \to 0$, $(A_k)_{12} = 2/k \to 0$, $(A_k)_{21} = 0 \to 0$, $(A_k)_{22} = 1/k^2 \to 0$.

    Therefore $\lim_{k \to \infty} A_k = O$ (the zero matrix).

    Verification with the Frobenius norm: $\|A_k\|_F = \sqrt{1/k^2 + 4/k^2 + 0 + 1/k^4} = \sqrt{5/k^2 + 1/k^4} \to 0$.

!!! example "Example 14.2"
    Let $A = \begin{pmatrix} 0.5 & 0.1 \\ 0 & 0.3 \end{pmatrix}$, and consider the sequence $\{A^k\}$.

    Since $A$ is upper triangular, its eigenvalues are the diagonal entries $\lambda_1 = 0.5$, $\lambda_2 = 0.3$, both satisfying $|\lambda_i| < 1$.

    Direct computation gives $A^2 = \begin{pmatrix} 0.25 & 0.08 \\ 0 & 0.09 \end{pmatrix}$, $A^3 = \begin{pmatrix} 0.125 & 0.049 \\ 0 & 0.027 \end{pmatrix}$.

    One can observe that all entries tend to $0$, so $A^k \to O$. This is consistent with the theory in Section 14.7.

---

## 14.2 Matrix series

<div class="context-flow" markdown>

**Chapter arc**: Absolute convergence (completeness) → **Neumann series** $\sum A^k = (I-A)^{-1}$ ←→ matrix generalization of geometric series · Convergence condition: $\rho(A)<1$

</div>

By analogy with scalar series, we can define matrix series and their convergence. One of the most important matrix series is the Neumann series, which gives a series representation of $(I - A)^{-1}$.

!!! definition "Definition 14.4 (Convergence of matrix series)"
    Let $\{A_k\}_{k=0}^{\infty} \subset \mathbb{C}^{n \times n}$. The matrix series $\sum_{k=0}^{\infty} A_k$ is **convergent** if its partial sum sequence $S_N = \sum_{k=0}^{N} A_k$ converges. If $\sum_{k=0}^{\infty} \|A_k\|$ converges (for some matrix norm $\|\cdot\|$), the series is said to be **absolutely convergent**.

!!! theorem "Theorem 14.3 (Absolute convergence implies convergence)"
    In $\mathbb{C}^{n \times n}$, an absolutely convergent matrix series is convergent.

??? proof "Proof"
    Suppose $\sum_{k=0}^{\infty}\|A_k\|$ converges. For $N > M$,

    $$\left\|\sum_{k=M+1}^{N} A_k\right\| \leq \sum_{k=M+1}^{N}\|A_k\| \to 0 \quad (M, N \to \infty).$$

    Thus the partial sums $\{S_N\}$ form a Cauchy sequence. By the completeness of $\mathbb{C}^{n \times n}$, the sequence converges. $\blacksquare$

!!! definition "Definition 14.5 (Neumann series)"
    Let $A \in \mathbb{C}^{n \times n}$. The series $\sum_{k=0}^{\infty} A^k$ is called the **Neumann series**.

!!! theorem "Theorem 14.4 (Neumann series convergence theorem)"
    Let $A \in \mathbb{C}^{n \times n}$. The following conditions are equivalent:

    (1) The Neumann series $\sum_{k=0}^{\infty} A^k$ converges;

    (2) $\rho(A) < 1$ (spectral radius less than $1$);

    (3) $A^k \to O$ ($k \to \infty$).

    When these conditions hold, $I - A$ is invertible and

    $$\sum_{k=0}^{\infty} A^k = (I - A)^{-1}.$$

??? proof "Proof"
    **(2) $\Rightarrow$ (3)**: See Theorem 14.10.

    **(3) $\Rightarrow$ (1)**: Let $S_N = \sum_{k=0}^{N} A^k$. Note that

    $$(I - A)S_N = S_N(I - A) = I - A^{N+1}.$$

    Since $A^{N+1} \to O$, the right-hand side tends to $I$. This shows that for sufficiently large $N$, $I - A^{N+1}$ is invertible, so $I - A$ is invertible (since $I - A$ does not depend on $N$; it is either invertible or not).

    More precisely, from $(I-A)S_N = I - A^{N+1}$, if $I - A$ is invertible then $S_N = (I-A)^{-1}(I - A^{N+1}) \to (I-A)^{-1}$.

    We now show that $I - A$ is indeed invertible: if $I - A$ were singular, then $1$ would be an eigenvalue of $A$; let $A\mathbf{v} = \mathbf{v}$ ($\mathbf{v} \neq \mathbf{0}$), then $A^k \mathbf{v} = \mathbf{v}$ for all $k$, contradicting $A^k \to O$.

    **(1) $\Rightarrow$ (2)**: If $\sum A^k$ converges, then $A^k \to O$ (necessary condition for series convergence). If there existed an eigenvalue $\lambda$ with $|\lambda| \geq 1$, say $A\mathbf{v} = \lambda\mathbf{v}$, then $A^k\mathbf{v} = \lambda^k \mathbf{v}$, $\|A^k\mathbf{v}\| = |\lambda|^k\|\mathbf{v}\| \geq \|\mathbf{v}\| > 0$, contradicting $A^k \to O$. Hence $\rho(A) < 1$. $\blacksquare$

!!! example "Example 14.3"
    Let $A = \begin{pmatrix} 0 & 1/2 \\ 1/3 & 0 \end{pmatrix}$.

    The characteristic polynomial is $\lambda^2 - 1/6 = 0$, so the eigenvalues are $\lambda = \pm 1/\sqrt{6}$. Thus $\rho(A) = 1/\sqrt{6} < 1$ and the Neumann series converges.

    $$(I - A)^{-1} = \begin{pmatrix} 1 & -1/2 \\ -1/3 & 1 \end{pmatrix}^{-1} = \frac{1}{1-1/6}\begin{pmatrix} 1 & 1/2 \\ 1/3 & 1 \end{pmatrix} = \frac{6}{5}\begin{pmatrix} 1 & 1/2 \\ 1/3 & 1 \end{pmatrix} = \begin{pmatrix} 6/5 & 3/5 \\ 2/5 & 6/5 \end{pmatrix}.$$

!!! theorem "Theorem 14.5 (Norm estimate for the Neumann series)"
    Let $\|\cdot\|$ be a submultiplicative norm on $\mathbb{C}^{n \times n}$ with $\|A\| < 1$. Then $I - A$ is invertible and

    $$\|(I - A)^{-1}\| \leq \frac{1}{1 - \|A\|}, \qquad \|(I-A)^{-1} - I\| \leq \frac{\|A\|}{1 - \|A\|}.$$

??? proof "Proof"
    Since $\|A\| < 1$ and by submultiplicativity, $\|A^k\| \leq \|A\|^k$, so $\sum_{k=0}^{\infty}\|A^k\| \leq \sum_{k=0}^{\infty}\|A\|^k = \frac{1}{1-\|A\|}$, and the Neumann series converges absolutely. Then

    $$\|(I-A)^{-1}\| = \left\|\sum_{k=0}^{\infty}A^k\right\| \leq \sum_{k=0}^{\infty}\|A^k\| \leq \frac{1}{1-\|A\|}.$$

    For the second inequality:

    $$\|(I-A)^{-1} - I\| = \left\|\sum_{k=1}^{\infty}A^k\right\| \leq \sum_{k=1}^{\infty}\|A\|^k = \frac{\|A\|}{1-\|A\|}. \quad \blacksquare$$

!!! example "Example 14.4"
    Using the Neumann series to approximately solve a linear system. Let $A = I - E$, where $E$ is a "small" perturbation matrix. If $\|E\| < 1$, then

    $$A^{-1} = (I - E)^{-1} = I + E + E^2 + \cdots$$

    Taking the first few terms as an approximation: $A^{-1} \approx I + E + E^2$.

    For example, $E = \begin{pmatrix} 0.1 & 0.05 \\ 0.02 & 0.1 \end{pmatrix}$, $\|E\|_\infty = 0.15 < 1$.

    Then $A^{-1} \approx I + E + E^2 = \begin{pmatrix} 1.111 & 0.06 \\ 0.024 & 1.111 \end{pmatrix}$ (to three decimal places).

---

## 14.3 Spectral radius

<div class="context-flow" markdown>

**Chapter arc**: $\rho(A) = \max|\lambda_i|$ → **Gelfand formula** $\rho(A)=\lim\|A^k\|^{1/k}$ bridges norms and spectra · $\rho(A)\leq\|A\|$ for any submultiplicative norm

</div>

The spectral radius is one of the most central concepts in matrix analysis. It characterizes the "size" of the eigenvalues of a matrix and determines the asymptotic behavior of matrix powers.

!!! definition "Definition 14.6 (Spectral radius)"
    Let $A \in \mathbb{C}^{n \times n}$ have eigenvalues $\lambda_1, \lambda_2, \ldots, \lambda_n$ (counting multiplicity). The **spectral radius** of $A$ is defined as

    $$\rho(A) = \max_{1 \leq i \leq n} |\lambda_i|.$$

!!! theorem "Theorem 14.6 (Spectral radius and norm)"
    Let $\|\cdot\|$ be any submultiplicative norm on $\mathbb{C}^{n \times n}$. Then

    $$\rho(A) \leq \|A\|.$$

??? proof "Proof"
    Let $\lambda$ be an eigenvalue of $A$ and $\mathbf{v}$ be a corresponding unit eigenvector (with $\|\mathbf{v}\| = 1$ in some vector norm). Consider the rank-one matrix $B = \mathbf{v}\mathbf{w}^*$, where $\mathbf{w}$ is chosen so that $\mathbf{w}^*\mathbf{v} = 1$.

    A more direct proof: let $\lambda$ be the eigenvalue of largest modulus, $A\mathbf{v} = \lambda\mathbf{v}$. Form the matrix $X = \mathbf{v}\mathbf{e}_1^T$ (where $\mathbf{e}_1$ is a standard basis vector); then $AX = \lambda\mathbf{v}\mathbf{e}_1^T = \lambda X$. By submultiplicativity, $\|AX\| \leq \|A\|\|X\|$, so $|\lambda|\|X\| \leq \|A\|\|X\|$. Since $X \neq O$, we get $|\lambda| \leq \|A\|$. Taking the maximum over all eigenvalues gives $\rho(A) \leq \|A\|$. $\blacksquare$

!!! theorem "Theorem 14.7 (Gelfand formula)"
    Let $A \in \mathbb{C}^{n \times n}$ and $\|\cdot\|$ be any submultiplicative norm. Then

    $$\rho(A) = \lim_{k \to \infty} \|A^k\|^{1/k}.$$

    This limit is independent of the choice of norm.

??? proof "Proof"
    **Lower bound**: Since $\rho(A^k) = \rho(A)^k$ (because the eigenvalues $\lambda_i$ of $A$ correspond to eigenvalues $\lambda_i^k$ of $A^k$) and by Theorem 14.6,

    $$\rho(A)^k = \rho(A^k) \leq \|A^k\|,$$

    so $\rho(A) \leq \|A^k\|^{1/k}$. Hence $\rho(A) \leq \liminf_{k\to\infty}\|A^k\|^{1/k}$.

    **Upper bound**: Let the Jordan normal form of $A$ be $A = PJP^{-1}$, where $J = \operatorname{diag}(J_1, J_2, \ldots, J_s)$. For any $\varepsilon > 0$, let $D_\varepsilon = \operatorname{diag}(1, \varepsilon, \varepsilon^2, \ldots, \varepsilon^{n-1})$; then the superdiagonal entries of $D_\varepsilon^{-1} J D_\varepsilon$ are multiplied by $\varepsilon$. By taking $\varepsilon$ sufficiently small, $\|D_\varepsilon^{-1}JD_\varepsilon\|$ can be made arbitrarily close to $\rho(A)$ in a suitable norm.

    More precisely, let $B_\varepsilon = (PD_\varepsilon)^{-1}A(PD_\varepsilon)$; then $\|B_\varepsilon\| \leq \rho(A) + \varepsilon$ (in a suitable norm). Therefore

    $$\|A^k\|^{1/k} = \|(PD_\varepsilon)B_\varepsilon^k(PD_\varepsilon)^{-1}\|^{1/k} \leq \left(\|PD_\varepsilon\|\cdot\|(PD_\varepsilon)^{-1}\|\right)^{1/k}(\rho(A)+\varepsilon).$$

    As $k \to \infty$, $\left(\|PD_\varepsilon\|\cdot\|(PD_\varepsilon)^{-1}\|\right)^{1/k} \to 1$. Hence $\limsup_{k\to\infty}\|A^k\|^{1/k} \leq \rho(A) + \varepsilon$. Since $\varepsilon$ is arbitrary, $\limsup_{k\to\infty}\|A^k\|^{1/k} \leq \rho(A)$.

    Combining the lower and upper bounds, the limit exists and equals $\rho(A)$. $\blacksquare$

!!! proposition "Proposition 14.1 (Basic properties of the spectral radius)"
    Let $A \in \mathbb{C}^{n \times n}$. Then:

    (1) $\rho(A^T) = \rho(A)$, $\rho(\bar{A}) = \rho(A)$, $\rho(A^*) = \rho(A)$;

    (2) $\rho(\alpha A) = |\alpha|\rho(A)$ for any $\alpha \in \mathbb{C}$;

    (3) $\rho(A^k) = \rho(A)^k$ for any positive integer $k$;

    (4) If $A$ is normal ($A^*A = AA^*$), then $\rho(A) = \|A\|_2$ (operator 2-norm);

    (5) For any $\varepsilon > 0$, there exists a submultiplicative norm $\|\cdot\|$ such that $\|A\| \leq \rho(A) + \varepsilon$.

!!! example "Example 14.5"
    Compute the spectral radius of $A = \begin{pmatrix} 2 & 1 \\ 0 & -3 \end{pmatrix}$.

    The eigenvalues are $\lambda_1 = 2$, $\lambda_2 = -3$, so $\rho(A) = \max\{|2|, |-3|\} = 3$.

    Verification via the Gelfand formula: $A^2 = \begin{pmatrix} 4 & -1 \\ 0 & 9 \end{pmatrix}$, $\|A^2\|_\infty = \max\{5, 9\} = 9$, $\|A^2\|_\infty^{1/2} = 3$.

    $A^4 = \begin{pmatrix} 16 & -10 \\ 0 & 81 \end{pmatrix}$, $\|A^4\|_\infty^{1/4} = 81^{1/4} = 3$. The sequence tends to $\rho(A) = 3$.

---

## 14.4 Gershgorin disc theorem

<div class="context-flow" markdown>

**Intuition**: Without solving the characteristic polynomial, one can locate eigenvalues using only matrix entries · Center = diagonal entry, radius = deleted row sum → diagonal dominance $\Rightarrow$ invertibility

</div>

The Gershgorin disc theorem is one of the most elegant and practical results in eigenvalue localization theory. It gives eigenvalue location estimates using only the matrix entries.

!!! definition "Definition 14.7 (Gershgorin disc)"
    Let $A = (a_{ij}) \in \mathbb{C}^{n \times n}$. The $i$-th **Gershgorin disc** is defined as

    $$D_i = \left\{z \in \mathbb{C} : |z - a_{ii}| \leq R_i\right\}, \quad R_i = \sum_{j \neq i} |a_{ij}|,$$

    where $R_i$ is called the **deleted row sum** of the $i$-th row.

!!! theorem "Theorem 14.8 (Gershgorin disc theorem)"
    Let $A = (a_{ij}) \in \mathbb{C}^{n \times n}$. Every eigenvalue of $A$ belongs to at least one Gershgorin disc, i.e.,

    $$\sigma(A) \subseteq \bigcup_{i=1}^{n} D_i.$$

??? proof "Proof"
    Let $\lambda$ be an eigenvalue of $A$ and $\mathbf{x} = (x_1, x_2, \ldots, x_n)^T$ be a corresponding eigenvector. Choose $p$ so that $|x_p| = \max_{1 \leq i \leq n}|x_i|$. Since $\mathbf{x} \neq \mathbf{0}$, we have $|x_p| > 0$.

    From the $p$-th component of $A\mathbf{x} = \lambda\mathbf{x}$:

    $$\sum_{j=1}^{n} a_{pj}x_j = \lambda x_p.$$

    Therefore

    $$(\lambda - a_{pp})x_p = \sum_{j \neq p} a_{pj}x_j.$$

    Taking the modulus and using the triangle inequality:

    $$|\lambda - a_{pp}| \cdot |x_p| \leq \sum_{j \neq p} |a_{pj}|\cdot|x_j| \leq |x_p| \sum_{j \neq p} |a_{pj}| = |x_p| R_p.$$

    Dividing by $|x_p| > 0$ gives $|\lambda - a_{pp}| \leq R_p$, i.e., $\lambda \in D_p$. $\blacksquare$

!!! theorem "Theorem 14.9 (Connected components of Gershgorin discs)"
    If the union of the $n$ Gershgorin discs can be partitioned into $k$ mutually disjoint connected regions, each formed by the union of $m_1, m_2, \ldots, m_k$ discs respectively ($m_1 + m_2 + \cdots + m_k = n$), then the $j$-th connected region contains exactly $m_j$ eigenvalues of $A$ (counting multiplicity).

??? proof "Proof"
    Consider the matrix family $A(t) = D + t(A - D)$, where $D = \operatorname{diag}(a_{11}, \ldots, a_{nn})$, $t \in [0, 1]$. At $t = 0$, $A(0) = D$ with eigenvalues $a_{11}, \ldots, a_{nn}$; at $t = 1$, $A(1) = A$.

    The radius of the $i$-th Gershgorin disc of $A(t)$ is $tR_i$. As $t$ increases continuously from $0$ to $1$, the discs expand continuously. Eigenvalues are continuous functions of the coefficients of the characteristic polynomial, so they vary continuously with $t$.

    At $t = 0$, the $m_j$ discs (degenerate to points) correspond to $m_j$ eigenvalues. Since eigenvalues vary continuously and the connected regions remain disjoint, eigenvalues cannot jump from one connected region to another. Hence each connected region contains exactly $m_j$ eigenvalues. $\blacksquare$

!!! corollary "Corollary 14.1 (Strictly diagonally dominant matrices are invertible)"
    If $A$ is a **strictly diagonally dominant matrix**, i.e., $|a_{ii}| > R_i = \sum_{j \neq i}|a_{ij}|$ for every $i$, then $A$ is invertible.

??? proof "Proof"
    If $A$ were singular, then $0$ would be an eigenvalue of $A$. By the Gershgorin theorem, $0 \in D_i$ for some $i$, i.e., $|a_{ii}| \leq R_i$, contradicting the strict diagonal dominance condition. $\blacksquare$

!!! example "Example 14.6"
    Let $A = \begin{pmatrix} 4 & -1 & 0 \\ 1 & 5 & -1 \\ 0 & -1 & 3 \end{pmatrix}$.

    The three Gershgorin discs are:

    - $D_1$: center $4$, radius $1$, i.e., $\{z : |z - 4| \leq 1\} = [3, 5]$;
    - $D_2$: center $5$, radius $2$, i.e., $\{z : |z - 5| \leq 2\} = [3, 7]$;
    - $D_3$: center $3$, radius $1$, i.e., $\{z : |z - 3| \leq 1\} = [2, 4]$.

    The union of the three discs is $[2, 7]$ (one connected region), so all three eigenvalues lie in $[2, 7]$.

    The actual eigenvalues are approximately $\lambda \approx 2.38, 4.18, 5.44$, all indeed in $[2, 7]$.

!!! example "Example 14.7"
    Let $A = \begin{pmatrix} 10 & 0.1 & 0.2 \\ 0.1 & 20 & 0.3 \\ 0.2 & 0.1 & 30 \end{pmatrix}$.

    The three Gershgorin discs are:

    - $D_1$: $|z - 10| \leq 0.3$, i.e., $[9.7, 10.3]$;
    - $D_2$: $|z - 20| \leq 0.4$, i.e., $[19.6, 20.4]$;
    - $D_3$: $|z - 30| \leq 0.3$, i.e., $[29.7, 30.3]$.

    The three discs are mutually disjoint, so by Theorem 14.9, each disc contains exactly one eigenvalue. In particular, $A$ has three well-separated real eigenvalues near $10, 20, 30$ respectively.

---

## 14.5 Matrix calculus

<div class="context-flow" markdown>

**Chapter arc**: Entrywise differentiation → product rule (note the order) → $\frac{d}{dt}A^{-1}=-A^{-1}\dot{A}A^{-1}$ · Key application: $\frac{d}{dt}e^{tA}=Ae^{tA}$ → Ch20 fundamental matrix of ODEs

</div>

Extending calculus to matrix-valued functions and functions of matrix variables is an important part of matrix analysis, with wide applications in optimization, control theory, and statistics.

!!! definition "Definition 14.8 (Derivative of a matrix-valued function)"
    Let $A(t) = (a_{ij}(t))$ be a matrix-valued function defined on an interval $I \subseteq \mathbb{R}$, where each $a_{ij}(t)$ is a real-valued function of $t$. The **derivative** of $A(t)$ with respect to $t$ is defined as

    $$\frac{dA}{dt} = \left(\frac{da_{ij}}{dt}\right),$$

    i.e., differentiation is performed entrywise.

!!! theorem "Theorem 14.10 (Differentiation rules for matrix-valued functions)"
    Let $A(t)$ and $B(t)$ be differentiable matrix-valued functions, and $c(t)$ a differentiable scalar function. Then:

    (1) $\frac{d}{dt}(A + B) = \frac{dA}{dt} + \frac{dB}{dt}$;

    (2) $\frac{d}{dt}(cA) = \frac{dc}{dt}A + c\frac{dA}{dt}$;

    (3) $\frac{d}{dt}(AB) = \frac{dA}{dt}B + A\frac{dB}{dt}$ (note the order);

    (4) If $A(t)$ is invertible, then $\frac{d}{dt}A^{-1} = -A^{-1}\frac{dA}{dt}A^{-1}$.

??? proof "Proof"
    (1)--(3) follow directly from entrywise differentiation.

    **(4)** From $A(t)A^{-1}(t) = I$, differentiating both sides with respect to $t$:

    $$\frac{dA}{dt}A^{-1} + A\frac{dA^{-1}}{dt} = O.$$

    Solving for $\frac{dA^{-1}}{dt} = -A^{-1}\frac{dA}{dt}A^{-1}$. $\blacksquare$

!!! definition "Definition 14.9 (Derivative of a scalar function with respect to a matrix)"
    Let $f : \mathbb{R}^{m \times n} \to \mathbb{R}$ be a scalar-valued function. The **derivative** (or **gradient**) of $f$ with respect to the matrix $X = (x_{ij})$ is defined as

    $$\frac{\partial f}{\partial X} = \left(\frac{\partial f}{\partial x_{ij}}\right) \in \mathbb{R}^{m \times n}.$$

!!! proposition "Proposition 14.2 (Common matrix derivative formulas)"
    Let $A$ be a constant matrix, $X$ a matrix variable, and $\mathbf{x}$ a vector variable. Then:

    (1) $\frac{\partial}{\partial \mathbf{x}}(\mathbf{a}^T\mathbf{x}) = \mathbf{a}$;

    (2) $\frac{\partial}{\partial \mathbf{x}}(\mathbf{x}^T A\mathbf{x}) = (A + A^T)\mathbf{x}$; if $A$ is symmetric, this equals $2A\mathbf{x}$;

    (3) $\frac{\partial}{\partial X}\operatorname{tr}(AX) = A^T$;

    (4) $\frac{\partial}{\partial X}\operatorname{tr}(X^TAX) = (A + A^T)X$;

    (5) $\frac{\partial}{\partial X}\ln\det(X) = X^{-T}$ (when $X$ is invertible).

!!! example "Example 14.8"
    Let $f(\mathbf{x}) = \mathbf{x}^TA\mathbf{x} + \mathbf{b}^T\mathbf{x} + c$, where $A$ is an $n \times n$ symmetric matrix, $\mathbf{b} \in \mathbb{R}^n$, $c \in \mathbb{R}$.

    Then $\frac{\partial f}{\partial \mathbf{x}} = 2A\mathbf{x} + \mathbf{b}$.

    Setting $\frac{\partial f}{\partial \mathbf{x}} = \mathbf{0}$ gives the critical point $\mathbf{x}^* = -\frac{1}{2}A^{-1}\mathbf{b}$ (when $A$ is invertible).

    The second derivative (Hessian matrix) is $\frac{\partial^2 f}{\partial \mathbf{x}\partial \mathbf{x}^T} = 2A$. If $A$ is positive definite, then $\mathbf{x}^*$ is a strict minimum.

!!! example "Example 14.9"
    **Derivative of the matrix exponential.** Let $A$ be a constant matrix. The matrix exponential $e^{tA} = \sum_{k=0}^{\infty}\frac{(tA)^k}{k!}$ satisfies

    $$\frac{d}{dt}e^{tA} = Ae^{tA} = e^{tA}A.$$

    This can be obtained by term-by-term differentiation:

    $$\frac{d}{dt}\sum_{k=0}^{\infty}\frac{t^k A^k}{k!} = \sum_{k=1}^{\infty}\frac{kt^{k-1}A^k}{k!} = A\sum_{k=1}^{\infty}\frac{t^{k-1}A^{k-1}}{(k-1)!} = Ae^{tA}.$$

    This is the fundamental matrix solution of the linear ODE system $\frac{d\mathbf{x}}{dt} = A\mathbf{x}$.

---

## 14.6 Eigenvalues and singular values

<div class="context-flow" markdown>

**Intuition**: $\prod_{i=1}^k|\lambda_i|\leq\prod_{i=1}^k\sigma_i$; for $k=n$ equality holds ($=|\det A|$) · Weyl singular value inequalities → Ch15 perturbation theory · Ch18 inequality network

</div>

Eigenvalues and singular values are two different but closely related sets of spectral information. Singular values are determined by the eigenvalues of $A^*A$, while eigenvalues come directly from $A$ itself. This section explores the deep connections between them.

!!! definition "Definition 14.10 (Singular values)"
    Let $A \in \mathbb{C}^{m \times n}$. The **singular values** $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_{\min(m,n)} \geq 0$ of $A$ are the non-negative square roots of the eigenvalues of $A^*A$, i.e., $\sigma_i = \sqrt{\lambda_i(A^*A)}$.

!!! theorem "Theorem 14.11 (Modulus inequality between eigenvalues and singular values)"
    Let $A \in \mathbb{C}^{n \times n}$ have eigenvalues ordered by decreasing modulus $|\lambda_1| \geq |\lambda_2| \geq \cdots \geq |\lambda_n|$, and singular values in decreasing order $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_n$. Then for each $k = 1, 2, \ldots, n$,

    $$\prod_{i=1}^{k} |\lambda_i| \leq \prod_{i=1}^{k} \sigma_i.$$

    In particular, $|\lambda_1| \leq \sigma_1$, and $\prod_{i=1}^{n}|\lambda_i| = \prod_{i=1}^{n}\sigma_i = |\det A|$.

??? proof "Proof"
    By the Schur decomposition $A = UTU^*$ ($U$ unitary, $T$ upper triangular), $A^*A = UT^*T U^*$.

    For $k = 1$: $|\lambda_1| \leq \|A\|_2 = \sigma_1$, since $\|A\|_2$ is the operator 2-norm.

    For $k = n$: $\prod_{i=1}^{n}|\lambda_i| = |\det A| = \prod_{i=1}^{n}\sigma_i$, since $|\det A|^2 = \det(A^*A) = \prod \sigma_i^2$.

    The general case requires compound matrix theory. Let $C_k(A)$ be the $k$-th compound matrix of $A$, whose eigenvalues are all products of $k$ eigenvalues $\lambda_{i_1}\cdots\lambda_{i_k}$. Then

    $$\prod_{i=1}^{k}|\lambda_i| \leq \|C_k(A)\|_2 = \sigma_1(C_k(A)) = \prod_{i=1}^{k}\sigma_i(A). \quad \blacksquare$$

!!! theorem "Theorem 14.12 (Weyl inequality -- singular value version)"
    Let $A, B \in \mathbb{C}^{m \times n}$ with singular values in decreasing order $\sigma_1(A) \geq \cdots$ and $\sigma_1(B) \geq \cdots$. Then for all $i$,

    $$|\sigma_i(A) - \sigma_i(B)| \leq \|A - B\|_2.$$

    In particular, $|\sigma_1(A) - \sigma_1(B)| \leq \|A - B\|_2$.

??? proof "Proof"
    The complete proof is given in Chapter 15. The idea is as follows: using the minimax characterization of singular values

    $$\sigma_i(A) = \min_{\dim V = n-i+1}\max_{\substack{\mathbf{x} \in V \\ \|\mathbf{x}\|=1}}\|A\mathbf{x}\|,$$

    and the triangle inequality $\|A\mathbf{x}\| \leq \|B\mathbf{x}\| + \|(A-B)\mathbf{x}\| \leq \|B\mathbf{x}\| + \|A-B\|_2$, the conclusion follows. $\blacksquare$

!!! example "Example 14.10"
    Let $A = \begin{pmatrix} 3 & 1 \\ 0 & 2 \end{pmatrix}$.

    **Eigenvalues**: $\lambda_1 = 3$, $\lambda_2 = 2$; $|\lambda_1| = 3$, $|\lambda_2| = 2$.

    **Singular values**: $A^*A = \begin{pmatrix} 9 & 3 \\ 3 & 5 \end{pmatrix}$, with eigenvalues $\frac{14 \pm \sqrt{36}}{2} = \frac{14 \pm 6}{2}$, i.e., $10$ and $4$. Singular values $\sigma_1 = \sqrt{10} \approx 3.16$, $\sigma_2 = 2$.

    Verification: $|\lambda_1| = 3 \leq 3.16 = \sigma_1$. $|\lambda_1||\lambda_2| = 6 = \sigma_1\sigma_2 = \sqrt{10}\cdot 2 = 2\sqrt{10} \approx 6.32$?

    Wait, let us re-verify. $|\det A| = |3 \cdot 2 - 1 \cdot 0| = 6$. $\sigma_1\sigma_2 = \sqrt{\det(A^*A)} = \sqrt{45 - 9} = \sqrt{36} = 6$. So $\prod|\lambda_i| = 6 = \prod\sigma_i = 6$.

    Recomputing: $A^*A = \begin{pmatrix} 9 & 3 \\ 3 & 5 \end{pmatrix}$, $\det(A^*A) = 45 - 9 = 36$, $\operatorname{tr}(A^*A) = 14$. Eigenvalues satisfy $\mu^2 - 14\mu + 36 = 0$, $\mu = 7 \pm \sqrt{13}$. So $\sigma_1 = \sqrt{7+\sqrt{13}} \approx 3.21$, $\sigma_2 = \sqrt{7-\sqrt{13}} \approx 1.87$. $\sigma_1\sigma_2 = \sqrt{36} = 6 = |\lambda_1\lambda_2|$.

    $|\lambda_1| = 3 \leq 3.21 = \sigma_1$. $|\lambda_1||\lambda_2| = 6 = \sigma_1\sigma_2$. The inequalities hold.

---

## 14.7 Limits and convergence of matrices

<div class="context-flow" markdown>

**Key result**: $A^k\to O \Leftrightarrow \rho(A)<1$ · Rate $\|A^k\|\leq Cr^k$ ($r>\rho(A)$) → Ch17 Markov chain convergence · Jacobi/Gauss--Seidel iteration convergence criteria

</div>

The convergence behavior of matrix powers $A^k$ is completely determined by the spectral radius. This result is of fundamental importance in iterative methods, Markov chains, and dynamical systems.

!!! theorem "Theorem 14.13 (Necessary and sufficient condition for convergence of matrix powers)"
    Let $A \in \mathbb{C}^{n \times n}$. Then

    $$\lim_{k \to \infty} A^k = O \quad \Longleftrightarrow \quad \rho(A) < 1.$$

??? proof "Proof"
    **($\Leftarrow$)**: Suppose $\rho(A) < 1$. By the Gelfand formula, for any $r$ with $\rho(A) < r < 1$, there exists a positive integer $N$ such that $\|A^k\|^{1/k} < r$ for $k \geq N$, i.e., $\|A^k\| < r^k$. Therefore $\|A^k\| \to 0$, i.e., $A^k \to O$.

    **($\Rightarrow$)**: Suppose $A^k \to O$. If $\rho(A) \geq 1$, then there exists an eigenvalue $\lambda$ with $|\lambda| \geq 1$. Let $A\mathbf{v} = \lambda\mathbf{v}$ ($\mathbf{v} \neq \mathbf{0}$); then $A^k\mathbf{v} = \lambda^k\mathbf{v}$, $\|A^k\mathbf{v}\| = |\lambda|^k\|\mathbf{v}\| \geq \|\mathbf{v}\| > 0$. But $A^k \to O$ implies $A^k\mathbf{v} \to \mathbf{0}$, a contradiction. $\blacksquare$

!!! theorem "Theorem 14.14 (Boundedness of matrix powers)"
    Let $A \in \mathbb{C}^{n \times n}$. The sequence $\{A^k\}_{k=0}^{\infty}$ is bounded (i.e., $\sup_k \|A^k\| < \infty$) if and only if $\rho(A) \leq 1$ and every eigenvalue of modulus $1$ is semisimple (i.e., its algebraic multiplicity equals its geometric multiplicity, or equivalently, it corresponds to $1 \times 1$ Jordan blocks).

??? proof "Proof"
    Let $A = PJP^{-1}$, where $J$ is the Jordan normal form. Then $A^k = PJ^kP^{-1}$.

    **Sufficiency**: If $\rho(A) \leq 1$ and eigenvalues of modulus $1$ are semisimple, then each Jordan block $J_i$ either corresponds to $|\lambda_i| < 1$ (so $J_i^k \to O$) or is a $1 \times 1$ block $(\lambda_i)$ with $|\lambda_i| = 1$ (so $|J_i^k| = |\lambda_i|^k = 1$). In either case $J_i^k$ is bounded, hence $J^k$ is bounded and $A^k = PJ^kP^{-1}$ is bounded.

    **Necessity**: If $\rho(A) > 1$, by the Gelfand formula $\|A^k\|^{1/k} \to \rho(A) > 1$, so $\|A^k\| \to \infty$. If there exists an eigenvalue $\lambda$ of modulus $1$ with a Jordan block of size $s \geq 2$, then the $(1,2)$ entry of $J_i^k$ is $k\lambda^{k-1}$, whose modulus is $k \to \infty$, which is unbounded. $\blacksquare$

!!! proposition "Proposition 14.3 (Rate of convergence of matrix powers)"
    If $\rho(A) < 1$, then for any $r$ with $\rho(A) < r < 1$, there exists a constant $C > 0$ such that

    $$\|A^k\| \leq Cr^k, \quad \forall k \geq 0.$$

    That is, $A^k$ tends to zero at an exponential rate determined by $\rho(A)$.

!!! example "Example 14.11"
    Determine whether the powers of the following matrices converge to the zero matrix.

    (a) $A = \begin{pmatrix} 0.5 & 0.3 \\ 0.1 & 0.4 \end{pmatrix}$.

    Characteristic polynomial: $\lambda^2 - 0.9\lambda + 0.17 = 0$, $\lambda = \frac{0.9 \pm \sqrt{0.81 - 0.68}}{2} = \frac{0.9 \pm \sqrt{0.13}}{2}$.

    $\sqrt{0.13} \approx 0.361$, so $\lambda_1 \approx 0.63$, $\lambda_2 \approx 0.27$. $\rho(A) \approx 0.63 < 1$, hence $A^k \to O$.

    (b) $B = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$.

    Eigenvalues are $\pm i$, $\rho(B) = 1$. $B$ is a normal matrix and the eigenvalues are semisimple, so $\{B^k\}$ is bounded but does not tend to zero. In fact $B^2 = -I$, $B^4 = I$, and the sequence oscillates with period $4$.

!!! example "Example 14.12"
    **Convergence of the Jacobi iteration.** Consider the linear system $A\mathbf{x} = \mathbf{b}$, with the splitting $A = D - L - U$ ($D$ is the diagonal part, $-L$ is the strictly lower triangular part, $-U$ is the strictly upper triangular part). The Jacobi iteration matrix is $T_J = D^{-1}(L + U)$.

    The Jacobi iteration $\mathbf{x}^{(k+1)} = T_J\mathbf{x}^{(k)} + D^{-1}\mathbf{b}$ converges if and only if $\rho(T_J) < 1$.

    For example, for $A = \begin{pmatrix} 4 & 1 \\ 1 & 3 \end{pmatrix}$, $T_J = \begin{pmatrix} 0 & -1/4 \\ -1/3 & 0 \end{pmatrix}$, with eigenvalues $\pm\frac{1}{2\sqrt{3}}$, $\rho(T_J) = \frac{1}{2\sqrt{3}} \approx 0.289 < 1$, so the Jacobi iteration converges.
