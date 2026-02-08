# Chapter 17  Nonnegative Matrices and Perron-Frobenius Theory

<div class="context-flow" markdown>

**Prerequisites**: Spectral radius (Ch14) · Directed graphs · **Chapter arc**: Nonnegative matrices → Irreducible = strongly connected → **Perron-Frobenius**: existence of a dominant positive eigenvalue → Markov chain convergence / PageRank

</div>

Nonnegative matrix theory is a beautiful and profound branch of linear algebra, whose core is the Perron-Frobenius theorem — it reveals that the spectral structure of nonnegative matrices possesses strikingly regular patterns. The theory originated from Perron's (1907) work on positive matrices and Frobenius's (1912) generalization to irreducible nonnegative matrices, and later became the mathematical foundation of probability theory (Markov chains), economics (input-output models), network analysis (the PageRank algorithm), and population dynamics.

---

## 17.1 Definition of nonnegative matrices

<div class="context-flow" markdown>

**Note**: Entrywise order $A\geq B$ (elementwise) $\neq$ Loewner order $A\succeq B$ (Ch16, positive semidefinite difference) · Nonnegative matrices are closed under addition and multiplication

</div>

!!! definition "Definition 17.1 (Nonnegative matrices and positive matrices)"
    Let $A = (a_{ij}) \in \mathbb{R}^{m \times n}$.

    - If $a_{ij} \geq 0$ for all $i, j$, then $A$ is called a **nonnegative matrix**, written $A \geq O$.
    - If $a_{ij} > 0$ for all $i, j$, then $A$ is called a **positive matrix**, written $A > O$.
    - Nonnegative vectors $\mathbf{x} \geq \mathbf{0}$ and positive vectors $\mathbf{x} > \mathbf{0}$ are defined analogously.

!!! definition "Definition 17.2 (Entrywise order of matrices)"
    For matrices $A, B$ of the same size, define $A \geq B$ to mean $A - B \geq O$ (i.e., $a_{ij} \geq b_{ij}$ for all $i, j$).

!!! note "Note"
    The entrywise order $A \geq B$ (elementwise comparison) and the Loewner partial order $A \succeq B$ ($A - B$ is positive semidefinite) are completely different concepts; they should not be confused.

!!! proposition "Proposition 17.1 (Basic properties of nonnegative matrices)"
    (1) If $A \geq O$ and $B \geq O$, then $A + B \geq O$;

    (2) If $A \geq O$ and $\alpha \geq 0$, then $\alpha A \geq O$;

    (3) If $A \geq O$ and $B \geq O$ with compatible dimensions, then $AB \geq O$;

    (4) If $A \geq O$, then $A^k \geq O$ for all positive integers $k$.

!!! example "Example 17.1"
    In practical applications, nonnegative matrices arise frequently:

    - **Transition probability matrices**: The transition matrix $P = (p_{ij})$ of a Markov chain, where $p_{ij} \geq 0$ and $\sum_j p_{ij} = 1$;
    - **Adjacency matrices**: The adjacency matrix $A$ of a directed graph, where $a_{ij} \in \{0, 1\}$;
    - **Input-output matrices**: The technical coefficient matrix in the Leontief economic model.

---

## 17.2 Perron's theorem (positive matrices)

<div class="context-flow" markdown>

**Core**: $A>O$ $\Rightarrow$ $\rho(A)$ is a simple algebraic root with a corresponding positive eigenvector that strictly dominates all other eigenvalues · Proof key: extremal argument on the simplex

</div>

!!! theorem "Theorem 17.1 (Perron's theorem)"
    Let $A > O$ be an $n \times n$ positive matrix ($n \geq 2$). Then:

    (1) $\rho(A) > 0$ and $\rho(A)$ is an eigenvalue of $A$ (called the **Perron root**);

    (2) $\rho(A)$ is a simple algebraic eigenvalue of $A$ (multiplicity $1$);

    (3) There exists a positive vector $\mathbf{v} > \mathbf{0}$ such that $A\mathbf{v} = \rho(A)\mathbf{v}$ (**Perron vector**);

    (4) Any other eigenvalue $\lambda$ of $A$ satisfies $|\lambda| < \rho(A)$ (strict inequality);

    (5) The Perron vector (after normalization) is the unique nonnegative eigenvector.

??? proof "Proof"
    **Existence**. Consider the compact set $\Delta = \{\mathbf{x} \in \mathbb{R}^n : \mathbf{x} \geq \mathbf{0},\; \sum x_i = 1\}$ (the simplex). Define

    $$r(\mathbf{x}) = \min_{i: x_i > 0} \frac{(A\mathbf{x})_i}{x_i}.$$

    Since $A > O$, for $\mathbf{x} \in \Delta$ ($\mathbf{x} \neq \mathbf{0}$), $A\mathbf{x} > \mathbf{0}$. On the interior $\Delta^\circ = \{\mathbf{x} > \mathbf{0}, \sum x_i = 1\}$, $r(\mathbf{x})$ is well-defined and continuous.

    Let $\rho^* = \sup_{\mathbf{x} \in \Delta} r(\mathbf{x})$. One can show that the supremum is attained at some $\mathbf{v}^* \in \Delta^\circ$, and $A\mathbf{v}^* = \rho^*\mathbf{v}^*$.

    In fact, $\rho^* = \rho(A)$. If there existed an eigenvalue $\lambda$ with $|\lambda| > \rho^*$, let $A\mathbf{w} = \lambda\mathbf{w}$ ($\mathbf{w}$ may be complex), then $A|\mathbf{w}| \geq |A\mathbf{w}| = |\lambda||\mathbf{w}|$ (entrywise), which would imply $r(|\mathbf{w}|/\||{\mathbf{w}}\||_1) \geq |\lambda| > \rho^*$, a contradiction.

    **Uniqueness (simple algebraic root)**. Suppose $\mathbf{w}$ is another linearly independent eigenvector corresponding to $\rho(A)$. Taking real and imaginary parts, one can find a real eigenvector. But a real eigenvector must have some zero or negative entries; let $A\mathbf{u} = \rho(A)\mathbf{u}$ with $\mathbf{u}$ having some non-positive entries. Since $A > O$, $\rho(A)\mathbf{u} = A\mathbf{u}$, taking absolute values yields $\rho(A)|\mathbf{u}| \leq A|\mathbf{u}|$, with strict inequality in some entries (because $\mathbf{u}$ has a mixture of positive and negative entries). This contradicts the maximality of $\rho(A)$.

    **Strict dominance**. Let $\lambda$ be another eigenvalue with $|\lambda| = \rho(A)$. Let $A\mathbf{w} = \lambda\mathbf{w}$. By the triangle inequality $\rho(A)|\mathbf{w}| = |\lambda\mathbf{w}| = |A\mathbf{w}| \leq A|\mathbf{w}|$. If strict inequality holds somewhere then $r(|\mathbf{w}|) > \rho(A)$, a contradiction. If equality holds everywhere, this requires all entries of $A\mathbf{w}$ to have consistent signs, which is only possible when $\mathbf{w}$ is a positive or negative vector, implying $\lambda = \rho(A)$ (a positive real number). $\blacksquare$

!!! example "Example 17.2"
    Let $A = \begin{pmatrix} 2 & 1 \\ 3 & 4 \end{pmatrix} > O$.

    Characteristic polynomial: $\lambda^2 - 6\lambda + 5 = 0$, $\lambda_1 = 5$, $\lambda_2 = 1$.

    $\rho(A) = 5$, the Perron root. $\lambda_2 = 1 < 5$.

    Perron vector: $A\mathbf{v} = 5\mathbf{v}$, $(A - 5I)\mathbf{v} = \mathbf{0}$, $\begin{pmatrix} -3 & 1 \\ 3 & -1 \end{pmatrix}\mathbf{v} = \mathbf{0}$, $\mathbf{v} = t\begin{pmatrix} 1 \\ 3 \end{pmatrix}$.

    Normalized Perron vector $\mathbf{v} = \begin{pmatrix} 1/4 \\ 3/4 \end{pmatrix}$ (with $\sum v_i = 1$), which is indeed a positive vector.

!!! example "Example 17.3"
    $A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix}$ is a positive matrix.

    Eigenvalues are approximately $\lambda_1 \approx 16.12$, $\lambda_2 \approx -1.12$, $\lambda_3 \approx 0$.

    $\rho(A) \approx 16.12$, a simple algebraic eigenvalue. $|\lambda_2| \approx 1.12 < 16.12$, $|\lambda_3| \approx 0 < 16.12$.

    The corresponding Perron vector (after normalization, all entries positive) is approximately $\mathbf{v} \approx (0.164, 0.387, 0.610)^T > \mathbf{0}$.

---

## 17.3 Irreducible matrices

<div class="context-flow" markdown>

**Intuition**: Irreducible $\Leftrightarrow$ directed graph is **strongly connected** $\Leftrightarrow$ $(I+A)^{n-1}>O$ · A bridge connecting linear algebra and graph theory

</div>

!!! definition "Definition 17.3 (Reducible and irreducible matrices)"
    Let $A \in \mathbb{R}^{n \times n}$ ($n \geq 2$). If there exists a permutation matrix $P$ such that

    $$P^TAP = \begin{pmatrix} B & C \\ O & D \end{pmatrix},$$

    where $B, D$ are square matrices with $B$ at least $1 \times 1$ and $D$ at least $1 \times 1$, then $A$ is called **reducible**. Otherwise $A$ is called **irreducible**.

    A $1 \times 1$ nonzero matrix is defined to be irreducible.

!!! definition "Definition 17.4 (Directed graph interpretation)"
    The **associated directed graph** $G(A)$ of an $n \times n$ nonnegative matrix $A$ has vertex set $\{1, 2, \ldots, n\}$, with a directed edge from $i$ to $j$ whenever $a_{ij} > 0$.

    $A$ is irreducible if and only if $G(A)$ is **strongly connected**, meaning that for any two vertices $i, j$, there exists a directed path from $i$ to $j$.

!!! theorem "Theorem 17.2 (Equivalent conditions for irreducibility)"
    Let $A \geq O$ be an $n \times n$ matrix ($n \geq 2$). The following conditions are equivalent:

    (1) $A$ is irreducible;

    (2) $G(A)$ is strongly connected;

    (3) $(I + A)^{n-1} > O$ (all entries are positive);

    (4) There is no index set $S \subsetneq \{1, \ldots, n\}$ ($S \neq \emptyset$) such that $a_{ij} = 0$ for all $i \in S, j \notin S$.

??? proof "Proof"
    **(1) $\Leftrightarrow$ (2)**: This is a standard result; the directed graph of a reducible matrix is not strongly connected.

    **(2) $\Leftrightarrow$ (3)**: The $(i,j)$ entry of $(I+A)^{n-1}$ is a certain combination of $\sum_{k=0}^{n-1}\binom{n-1}{k}(A^k)_{ij}$. $(I+A)^{n-1}_{ij} > 0$ if and only if there exists a path from $i$ to $j$ of length at most $n-1$ (including length-$0$ self-loops), which is equivalent to strong connectivity. $\blacksquare$

!!! example "Example 17.4"
    $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$ is irreducible.

    Associated directed graph: $1 \to 2 \to 3 \to 1$, forming a cycle, so it is strongly connected.

    $B = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}$ is reducible (already in upper triangular form; the corresponding graph has $1 \to 2$ but no $2 \to 1$).

---

## 17.4 Perron-Frobenius theorem (irreducible nonnegative matrices)

<div class="context-flow" markdown>

**Perron → Perron-Frobenius**: Generalization from positive matrices to irreducible nonnegative matrices · New phenomenon: index of imprimitivity $h$ eigenvalues uniformly distributed on the circle $|\lambda|=\rho$

</div>

!!! theorem "Theorem 17.3 (Perron-Frobenius theorem)"
    Let $A \geq O$ be an $n \times n$ irreducible nonnegative matrix. Then:

    (1) $\rho(A) > 0$ and $\rho(A)$ is an eigenvalue of $A$ (the Perron root);

    (2) The eigenvector corresponding to $\rho(A)$ can be taken to be a positive vector $\mathbf{v} > \mathbf{0}$ (the Perron vector);

    (3) The algebraic multiplicity of $\rho(A)$ as an eigenvalue is $1$;

    (4) If $A$ has $h$ eigenvalues with modulus equal to $\rho(A)$, then they are exactly

    $$\rho(A), \rho(A)e^{2\pi i/h}, \rho(A)e^{4\pi i/h}, \ldots, \rho(A)e^{2(h-1)\pi i/h},$$

    i.e., they are uniformly distributed on the circle of radius $\rho(A)$. The integer $h$ is called the **index of imprimitivity** of $A$;

    (5) If $\mathbf{x} \geq \mathbf{0}$, $\mathbf{x} \neq \mathbf{0}$, and $A\mathbf{x} \leq \lambda\mathbf{x}$ (entrywise), then $\lambda \geq \rho(A)$; if $A\mathbf{x} \geq \lambda\mathbf{x}$, then $\lambda \leq \rho(A)$.

??? proof "Proof"
    The core ideas of the proof are as follows.

    **Existence (Perron root and positive eigenvector)**. Define $\rho^* = \sup\{r \geq 0 : \exists \mathbf{x} \geq \mathbf{0}, \mathbf{x} \neq \mathbf{0}, A\mathbf{x} \geq r\mathbf{x}\}$. One can show that $\rho^* = \rho(A)$ and the supremum is attained.

    Let $\mathbf{x}^*$ be the nonnegative vector attaining the supremum. From $A\mathbf{x}^* \geq \rho^*\mathbf{x}^*$ and irreducibility, one can show that $\mathbf{x}^* > \mathbf{0}$ (because if some entry were zero, irreducibility guarantees that after several multiplications by $A$ all entries become positive). Furthermore, $A\mathbf{x}^* = \rho^*\mathbf{x}^*$ (otherwise $r$ could be increased).

    **Simple algebraic root**. Let $\rho = \rho(A)$. From the existence of a positive eigenvector $\mathbf{v} > \mathbf{0}$, one can construct the diagonal matrix $D = \operatorname{diag}(\mathbf{v})$, and $B = D^{-1}AD$ is a nonnegative matrix with row sums equal to $\rho$. If $\mathbf{w}$ is another eigenvector of $B$ corresponding to $\rho$, an argument similar to Perron's theorem shows that $\mathbf{w}$ must be a scalar multiple of $\mathbf{v}$, so the algebraic multiplicity is $1$.

    **Cyclic structure**. The greatest common divisor of the lengths of all cycles in the associated directed graph of the irreducible matrix is the index of imprimitivity $h$. By classifying the vertices according to their distance modulo $h$ from a fixed vertex, one can transform $A$ into a cyclic block form, thereby proving the rotational symmetry of the eigenvalues. $\blacksquare$

!!! example "Example 17.5"
    Let $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \\ 3 & 0 & 0 \end{pmatrix}$, which is irreducible ($1\to 2\to 3\to 1$).

    Characteristic polynomial: $\lambda^3 - 6 = 0$, $\lambda = \sqrt[3]{6}, \sqrt[3]{6}\omega, \sqrt[3]{6}\omega^2$, where $\omega = e^{2\pi i/3}$.

    $\rho(A) = \sqrt[3]{6} \approx 1.817$. All three eigenvalues have modulus $\sqrt[3]{6}$ and are uniformly distributed on the circle. Index of imprimitivity $h = 3$ (the only cycle in the graph has length $3$, so $\gcd = 3$).

    Perron vector: $A\mathbf{v} = \sqrt[3]{6}\mathbf{v}$, yielding $\mathbf{v} \propto (\sqrt[3]{36}, 2\sqrt[3]{6}, 6)^T / c > \mathbf{0}$.

---

## 17.5 Primitive matrices

<div class="context-flow" markdown>

**Chapter arc**: Primitive = irreducible + aperiodic ($h=1$) $\Leftrightarrow$ $\exists k, A^k>O$ · $(A/\rho)^k\to\mathbf{v}\mathbf{w}^T$ → Concretization of Ch14 power convergence for nonnegative matrices

</div>

!!! definition "Definition 17.5 (Primitive matrix)"
    An irreducible nonnegative matrix $A$ is called **primitive** if its index of imprimitivity $h = 1$, i.e., $\rho(A)$ is the only eigenvalue with modulus equal to $\rho(A)$.

    Equivalently, $A \geq O$ is primitive if and only if there exists a positive integer $k$ such that $A^k > O$ (all entries are positive).

!!! theorem "Theorem 17.4 (Criteria for primitivity)"
    Let $A \geq O$ be irreducible. The following conditions are equivalent:

    (1) $A$ is primitive;

    (2) The index of imprimitivity $h = 1$;

    (3) There exists $k$ such that $A^k > O$;

    (4) The greatest common divisor of all cycle lengths in the associated directed graph is $1$.

!!! theorem "Theorem 17.5 (Convergence of primitive matrices)"
    Let $A \geq O$ be an $n \times n$ primitive matrix with Perron root $\rho = \rho(A)$. Let $\mathbf{v}$ and $\mathbf{w}^T$ be the right and left Perron vectors, respectively, normalized so that $\mathbf{w}^T\mathbf{v} = 1$. Then

    $$\lim_{k \to \infty} \left(\frac{A}{\rho}\right)^k = \mathbf{v}\mathbf{w}^T.$$

??? proof "Proof"
    By the Jordan decomposition. Let the eigenvalues of $A$ be $\rho = \lambda_1 > |\lambda_2| \geq \cdots \geq |\lambda_n|$ (primitivity guarantees the strict inequality). The Jordan decomposition of $A$ is

    $$A = \rho\mathbf{v}\mathbf{w}^T + \sum_{i=2}^{n}\lambda_i P_i + N_i,$$

    where $P_i$ are projections and $N_i$ are nilpotent parts. Therefore

    $$\left(\frac{A}{\rho}\right)^k = \mathbf{v}\mathbf{w}^T + \sum_{i=2}^{n}\left(\frac{\lambda_i}{\rho}\right)^k(\cdots).$$

    Since $|\lambda_i/\rho| < 1$, the remaining terms tend to zero (the polynomial growth from nilpotent parts is overwhelmed by exponential decay). $\blacksquare$

!!! example "Example 17.6"
    $A = \begin{pmatrix} 0.5 & 0.5 \\ 0.3 & 0.7 \end{pmatrix}$ is a positive matrix, hence automatically primitive.

    Eigenvalues: $\lambda_1 = 1$, $\lambda_2 = 0.2$. $\rho(A) = 1$.

    Right Perron vector $\mathbf{v} = \begin{pmatrix} 3/8 \\ 5/8 \end{pmatrix}$ (normalized so that entries sum to $1$). Left Perron vector $\mathbf{w}^T = (1, 1)$ (since row sums are $1$). $\mathbf{w}^T\mathbf{v} = 1$.

    $A^k \to \mathbf{v}\mathbf{w}^T = \begin{pmatrix} 3/8 & 3/8 \\ 5/8 & 5/8 \end{pmatrix}$. Verification: $A^2 = \begin{pmatrix} 0.40 & 0.60 \\ 0.36 & 0.64 \end{pmatrix}$, $A^{10} \approx \begin{pmatrix} 0.375 & 0.375 \\ 0.625 & 0.625 \end{pmatrix}$ (the error comes from $0.2^{10} \approx 10^{-7}$).

    Wait, let me recalculate. $A\mathbf{v} = \mathbf{v}$, $\begin{pmatrix} 0.5 & 0.5 \\ 0.3 & 0.7 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix} = \begin{pmatrix} a \\ b \end{pmatrix}$, $0.5a + 0.5b = a$, $0.3a + 0.7b = b$, both give $a = b$... That is not right. $-0.5a + 0.5b = 0$ gives $a = b$. But $\mathbf{w}^TA = \mathbf{w}^T$, $\begin{pmatrix} w_1 & w_2 \end{pmatrix}\begin{pmatrix} 0.5 & 0.5 \\ 0.3 & 0.7 \end{pmatrix} = \begin{pmatrix} w_1 & w_2 \end{pmatrix}$, $0.5w_1 + 0.3w_2 = w_1$ gives $w_2 = \frac{5}{3}w_1$. Taking $\mathbf{w}^T = (3, 5)$, then $\mathbf{w}^T\mathbf{v} = 3+5 = 8$ (if $\mathbf{v} = (1,1)^T$). Normalizing: $\mathbf{v} = (1,1)^T$, $\mathbf{w}^T = (3/8, 5/8)$, $\mathbf{w}^T\mathbf{v} = 1$.

    $\mathbf{v}\mathbf{w}^T = \begin{pmatrix} 3/8 & 5/8 \\ 3/8 & 5/8 \end{pmatrix}$. This is the limit of $A^k$.

---

## 17.6 Stochastic matrices

<div class="context-flow" markdown>

**Chapter arc**: Row stochastic $\Rightarrow$ $\rho=1$, $\mathbf{1}$ is right eigenvector · **Birkhoff's theorem**: Doubly stochastic matrices = convex combinations of permutation matrices → Ch18 majorization

</div>

!!! definition "Definition 17.6 (Stochastic matrices)"
    Let $P = (p_{ij}) \in \mathbb{R}^{n \times n}$.

    - If $P \geq O$ and each row sums to $1$ (i.e., $P\mathbf{1} = \mathbf{1}$), then $P$ is called a **row stochastic matrix**;
    - If $P \geq O$ and each column sums to $1$ (i.e., $\mathbf{1}^TP = \mathbf{1}^T$), then $P$ is called a **column stochastic matrix**;
    - If $P$ is both row stochastic and column stochastic, it is called a **doubly stochastic matrix**.

!!! proposition "Proposition 17.2 (Spectral properties of stochastic matrices)"
    (1) The spectral radius of a row stochastic matrix is $1$, and $\mathbf{1} = (1, 1, \ldots, 1)^T$ is the right eigenvector corresponding to eigenvalue $1$;

    (2) All eigenvalues $\lambda$ satisfy $|\lambda| \leq 1$;

    (3) The left Perron vector of a doubly stochastic matrix is $\frac{1}{n}\mathbf{1}^T$.

!!! theorem "Theorem 17.6 (Birkhoff's theorem)"
    The set of $n \times n$ doubly stochastic matrices is a convex polytope (called the **Birkhoff polytope**), whose vertices are exactly all $n \times n$ permutation matrices. That is, every doubly stochastic matrix can be expressed as a convex combination of permutation matrices.

??? proof "Proof"
    **Extreme points are permutation matrices**: Let $P$ be a doubly stochastic matrix that is an extreme point of the Birkhoff polytope. If $P$ is not a permutation matrix, then $P$ has at least two entries that are neither zero nor one. Using Konig's theorem, one can find two distinct permutation matrices $Q_1, Q_2$ such that $P = \alpha Q_1 + (1-\alpha)Q_2$ ($0 < \alpha < 1$), contradicting $P$ being an extreme point.

    **Permutation matrices are extreme points**: If a permutation matrix $Q = \alpha P_1 + (1-\alpha)P_2$, since the entries of $Q$ are $0$ or $1$ and the entries of $P_1, P_2$ lie in $[0,1]$, it follows that $P_1 = P_2 = Q$.

    By the Krein-Milman theorem (or the finite-dimensional Minkowski theorem), a convex polytope equals the convex hull of its extreme points. $\blacksquare$

!!! example "Example 17.7"
    The doubly stochastic matrix $P = \begin{pmatrix} 1/3 & 1/3 & 1/3 \\ 1/3 & 1/3 & 1/3 \\ 1/3 & 1/3 & 1/3 \end{pmatrix} = \frac{1}{3}(I + Q + Q^2)$, where $Q = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$, is a convex combination of the three permutation matrices $I, Q, Q^2$ (each with weight $1/3$).

---

## 17.7 Markov chains

<div class="context-flow" markdown>

**Application**: Primitive transition matrix → unique stationary distribution $\boldsymbol{\pi}$ · $P^k\to\mathbf{1}\boldsymbol{\pi}^T$ means "forgetting the initial state" · Convergence rate determined by $|\lambda_2|$ (spectral gap)

</div>

!!! definition "Definition 17.7 (Markov chain and transition matrix)"
    Let $\{X_t\}_{t=0,1,2,\ldots}$ be a random process taking values in the finite state space $\{1, 2, \ldots, n\}$. If

    $$\Pr(X_{t+1} = j \mid X_t = i, X_{t-1}, \ldots, X_0) = \Pr(X_{t+1} = j \mid X_t = i) = p_{ij},$$

    then $\{X_t\}$ is a **(homogeneous) Markov chain**, and $P = (p_{ij})$ is its **transition matrix**. $P$ is a row stochastic matrix.

!!! definition "Definition 17.8 (Stationary distribution)"
    A row vector $\boldsymbol{\pi}^T = (\pi_1, \ldots, \pi_n)$ is called a **stationary distribution** of the Markov chain if

    $$\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T, \quad \boldsymbol{\pi} \geq \mathbf{0}, \quad \sum_i \pi_i = 1.$$

    That is, $\boldsymbol{\pi}$ is the nonnegative normalized left eigenvector of $P^T$ corresponding to eigenvalue $1$.

!!! theorem "Theorem 17.7 (Convergence theorem for Markov chains)"
    Let $P$ be an irreducible and aperiodic (i.e., primitive) transition matrix. Then:

    (1) There exists a unique stationary distribution $\boldsymbol{\pi}^T$, and $\boldsymbol{\pi} > \mathbf{0}$;

    (2) For any initial distribution $\mathbf{p}_0^T$, $\mathbf{p}_0^T P^k \to \boldsymbol{\pi}^T$ ($k \to \infty$);

    (3) $P^k \to \mathbf{1}\boldsymbol{\pi}^T$ (every row converges to $\boldsymbol{\pi}^T$).

??? proof "Proof"
    Since $P$ is a primitive row stochastic matrix, the Perron root $\rho(P) = 1$, and the corresponding right Perron vector is $\mathbf{1}$. The left Perron vector $\boldsymbol{\pi}^T$ (normalized so that $\sum\pi_i = 1$) is the stationary distribution.

    By the convergence theorem for primitive matrices (Theorem 17.5), $(P/\rho)^k = P^k \to \mathbf{1}\boldsymbol{\pi}^T$ (since $\rho = 1$).

    For any initial distribution $\mathbf{p}_0^T$, $\mathbf{p}_0^TP^k \to \mathbf{p}_0^T\mathbf{1}\boldsymbol{\pi}^T = 1 \cdot \boldsymbol{\pi}^T = \boldsymbol{\pi}^T$. $\blacksquare$

!!! example "Example 17.8"
    **Weather model**. Let state $1$ = sunny, $2$ = rainy. Transition matrix

    $$P = \begin{pmatrix} 0.9 & 0.1 \\ 0.5 & 0.5 \end{pmatrix}.$$

    $P > O$, hence primitive. The stationary distribution satisfies $\boldsymbol{\pi}^TP = \boldsymbol{\pi}^T$:

    $0.9\pi_1 + 0.5\pi_2 = \pi_1$, $0.1\pi_1 + 0.5\pi_2 = \pi_2$, $\pi_1 + \pi_2 = 1$.

    The first equation gives $\pi_2 = 0.2\pi_1$, substituting $\pi_1 + 0.2\pi_1 = 1$, $\pi_1 = 5/6$, $\pi_2 = 1/6$.

    In the long run, approximately $83.3\%$ of days are sunny and $16.7\%$ are rainy.

---

## 17.8 Spectral properties of nonnegative matrices

<div class="context-flow" markdown>

**Chapter arc**: Index of imprimitivity $h$ = $\gcd$ of cycle lengths → spectrum has $h$-fold rotational symmetry → cyclic block normal form · PageRank: $|\lambda_2|\leq\alpha$ guarantees fast convergence

</div>

!!! definition "Definition 17.9 (Index of imprimitivity)"
    The **index of imprimitivity** (or **period**) $h$ of an irreducible nonnegative matrix $A$ is defined as the greatest common divisor of the lengths of all cycles in the associated directed graph $G(A)$.

    If $h = 1$, $A$ is primitive; if $h > 1$, $A$ is **cyclic** or **imprimitive**.

!!! theorem "Theorem 17.8 (Rotational symmetry of the spectrum)"
    Let $A \geq O$ be irreducible with index of imprimitivity $h$. Then the spectrum of $A$ has $h$-fold rotational symmetry about the origin, i.e., if $\lambda$ is an eigenvalue of $A$, then $\lambda e^{2\pi i/h}$ is also an eigenvalue of $A$.

    Equivalently, the characteristic polynomial of $A$ contains only powers of $\lambda^h$, i.e., $p(\lambda) = \lambda^r q(\lambda^h)$ for some polynomial $q$.

??? proof "Proof"
    By the Perron-Frobenius theorem, the eigenvalues with modulus equal to $\rho(A)$ are $\rho(A)e^{2k\pi i/h}$ ($k = 0, 1, \ldots, h-1$).

    More generally, let $\omega = e^{2\pi i/h}$; there exists a diagonal matrix $D = \operatorname{diag}(\omega^{d_1}, \ldots, \omega^{d_n})$ ($d_i$ depends on the position of vertex $i$ in the cyclic classification) such that

    $$D^{-1}AD = \omega A.$$

    This means $A$ and $\omega A$ are similar, hence they have the same spectrum. That is, $\lambda$ is an eigenvalue $\Rightarrow$ $\omega\lambda$ is also an eigenvalue. $\blacksquare$

!!! proposition "Proposition 17.3 (Cyclic normal form)"
    Let $A \geq O$ be irreducible with index of imprimitivity $h$. Then there exists a permutation matrix $P$ such that

    $$P^TAP = \begin{pmatrix} O & A_{12} & O & \cdots & O \\ O & O & A_{23} & \cdots & O \\ \vdots & & & \ddots & \vdots \\ O & O & O & \cdots & A_{h-1,h} \\ A_{h1} & O & O & \cdots & O \end{pmatrix},$$

    i.e., $A$ can be permuted into a cyclic block form, where the nonzero blocks appear only on the "super-diagonal" and in the lower-left corner.

!!! example "Example 17.9"
    $A = \begin{pmatrix} 0 & 2 & 0 & 0 \\ 0 & 0 & 3 & 0 \\ 0 & 0 & 0 & 1 \\ 4 & 0 & 0 & 0 \end{pmatrix}$.

    Directed graph: $1 \to 2 \to 3 \to 4 \to 1$, with the only cycle having length $4$; index of imprimitivity $h = 4$.

    Characteristic polynomial: $\lambda^4 - 24 = 0$, $\lambda = \sqrt[4]{24}\cdot\omega^k$ ($k = 0,1,2,3$), where $\omega = i$.

    The four eigenvalues $\sqrt[4]{24}, i\sqrt[4]{24}, -\sqrt[4]{24}, -i\sqrt[4]{24}$ form a square in the complex plane, exhibiting $4$-fold rotational symmetry.

!!! example "Example 17.10"
    **Mathematical model of the PageRank algorithm**. Suppose the link structure among $n$ web pages is described by a directed graph, and column-normalizing the adjacency matrix yields a column stochastic matrix $H$ (the $j$-th column of $H$ represents the probability distribution of links from page $j$). The Google matrix is

    $$G = \alpha H + (1-\alpha)\frac{1}{n}\mathbf{1}\mathbf{1}^T,$$

    where $\alpha \in (0, 1)$ (typically $\alpha = 0.85$).

    $G$ is a column stochastic positive matrix ($G > O$), hence primitive. By the Perron-Frobenius theorem, $\rho(G) = 1$ is a simple algebraic root, and the corresponding positive eigenvector $\boldsymbol{\pi}$ (after normalization) gives the PageRank value of each page. The second largest eigenvalue modulus $|\lambda_2| \leq \alpha < 1$ guarantees fast convergence of the power iteration.
