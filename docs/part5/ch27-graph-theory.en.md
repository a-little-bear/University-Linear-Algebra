# Chapter 27  Linear Algebra in Graph Theory and Networks

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · spectral theorem for symmetric matrices (Ch7-8) · nonnegative matrices (Ch17) · **Arc**: Adjacency/incidence/Laplacian matrices (linear algebra representations of graphs) → graph spectrum (spectral graph theory) → Laplacian and connectivity (Kirchhoff's matrix tree theorem) → Cheeger inequality (spectral clustering) → random walks/PageRank (Markov chains) → expander graphs (Ramanujan bound) → graph coloring (eigenvalue bounds) → network flows (LP duality)
**Essence**: Combinatorial properties of graphs (connectivity, expansion, chromatic number) can be precisely characterized through matrix spectra (eigenvalues) — this is the deepest bridge between discrete mathematics and linear algebra

</div>

The intersection of graph theory and linear algebra gave birth to spectral graph theory, whose core idea is to encode the combinatorial structure of a graph as a matrix and then use spectral properties of the matrix to reveal the graph's global structure. This chapter begins with matrix representations of graphs and progresses through spectral theory, Laplacian matrices and connectivity, Cheeger inequality and graph partitioning, random walks and PageRank, expander graphs, graph coloring, and network flows with linear programming.

---

## 27.1 Matrix Representations of Graphs

<div class="context-flow" markdown>

**Three matrices**: Adjacency matrix $A$ (symmetric $\leftrightarrow$ undirected graph) → incidence matrix $B$ (vertex-edge relationships) → Laplacian $L = D - A$ (most important)
**Links**: Spectral properties of symmetric matrices from Ch7 apply directly to the adjacency matrix and Laplacian of graphs

</div>

The combinatorial information of a graph can be completely encoded as a matrix, with different matrix representations offering different advantages.

!!! definition "Definition 27.1 (Adjacency Matrix)"
    Let $G = (V, E)$ be a simple undirected graph with $n$ vertices. The **adjacency matrix** $A \in \mathbb{R}^{n \times n}$ of $G$ is defined as

    $$
    A_{ij} = \begin{cases} 1 & \text{if } \{i, j\} \in E, \\ 0 & \text{otherwise}. \end{cases}
    $$

    $A$ is symmetric ($A = A^T$) with zero diagonal (no self-loops). For weighted graphs, $A_{ij} = w_{ij} \ge 0$.

!!! definition "Definition 27.2 (Degree Matrix and Laplacian Matrix)"
    The **degree matrix** is $D = \operatorname{diag}(d_1, \ldots, d_n)$, where $d_i = \sum_j A_{ij}$. The **Laplacian matrix** is defined as

    $$
    L = D - A.
    $$

    $L$ is symmetric positive semidefinite, and $L\mathbf{1} = \mathbf{0}$ ($\mathbf{1}$ is the all-ones vector), so $0$ is always an eigenvalue of $L$. The **normalized Laplacian** is

    $$
    \mathcal{L} = D^{-1/2} L D^{-1/2} = I - D^{-1/2} A D^{-1/2}.
    $$

!!! definition "Definition 27.3 (Incidence Matrix)"
    Given an orientation of a directed graph $G$, the **incidence matrix** $B \in \mathbb{R}^{n \times m}$ ($m = |E|$) is defined as

    $$
    B_{ve} = \begin{cases} +1 & \text{if } v \text{ is the tail of edge } e, \\ -1 & \text{if } v \text{ is the head of edge } e, \\ 0 & \text{otherwise}. \end{cases}
    $$

    Key property: $L = BB^T$ (independent of the chosen orientation).

!!! theorem "Theorem 27.1 (Quadratic Form of the Laplacian)"
    For any $\mathbf{x} \in \mathbb{R}^n$,

    $$
    \mathbf{x}^T L \mathbf{x} = \sum_{\{i,j\} \in E} w_{ij}(x_i - x_j)^2.
    $$

    This directly proves $L \succeq 0$ (positive semidefinite), and $\mathbf{x}^T L \mathbf{x} = 0$ if and only if $\mathbf{x}$ is constant on each connected component.

??? proof "Proof"
    From $L = D - A$,

    $$
    \mathbf{x}^T L \mathbf{x} = \mathbf{x}^T D \mathbf{x} - \mathbf{x}^T A \mathbf{x} = \sum_i d_i x_i^2 - \sum_{\{i,j\} \in E} 2w_{ij} x_i x_j.
    $$

    Since $d_i = \sum_j w_{ij}$,

    $$
    \sum_i d_i x_i^2 = \sum_{\{i,j\} \in E} w_{ij}(x_i^2 + x_j^2).
    $$

    Therefore $\mathbf{x}^T L \mathbf{x} = \sum_{\{i,j\} \in E} w_{ij}(x_i^2 - 2x_i x_j + x_j^2) = \sum_{\{i,j\} \in E} w_{ij}(x_i - x_j)^2 \ge 0$.

    Equality holds if and only if $x_i = x_j$ for every edge $\{i, j\}$, i.e., $\mathbf{x}$ is constant on each connected component. $\blacksquare$

!!! example "Example 27.1"
    **Matrix representations of path and cycle graphs.** The path graph $P_4$ (4-vertex path $1 - 2 - 3 - 4$):

    $$
    A = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 1 & 0 & 1 & 0 \\ 0 & 1 & 0 & 1 \\ 0 & 0 & 1 & 0 \end{pmatrix}, \quad L = \begin{pmatrix} 1 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 1 \end{pmatrix}.
    $$

    The eigenvalues of $L$ are $0, 2 - \sqrt{2}, 2, 2 + \sqrt{2}$. The multiplicity of the zero eigenvalue is $1$, confirming that $P_4$ is connected.

---

## 27.2 The Spectrum of a Graph

<div class="context-flow" markdown>

**Core concept**: Graph spectrum = eigenvalue set of the adjacency matrix (or Laplacian) → the spectrum carries rich structural information about the graph (regularity, bipartiteness, diameter, etc.)
**Links**: Extremal properties of eigenvalues of symmetric matrices (Courant-Fischer) from Ch7 appear repeatedly in spectral graph theory

</div>

The spectrum of a graph — the set of eigenvalues of its matrix representation — encodes rich structural information.

!!! definition "Definition 27.4 (Graph Spectrum)"
    Let $G$ be a simple graph of order $n$. The **spectrum** of $G$ is the set of eigenvalues (with multiplicities) of the adjacency matrix $A$:

    $$
    \lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_n.
    $$

    The **Laplacian spectrum** consists of the eigenvalues of $L$: $0 = \mu_1 \le \mu_2 \le \cdots \le \mu_n$.

!!! theorem "Theorem 27.2 (Basic Properties of the Spectrum)"
    Let $G$ be a simple graph of order $n$ with adjacency spectrum $\lambda_1 \ge \cdots \ge \lambda_n$:

    1. $\sum_{i=1}^n \lambda_i = \operatorname{tr}(A) = 0$.
    2. $\sum_{i=1}^n \lambda_i^2 = \operatorname{tr}(A^2) = 2|E|$.
    3. $\sum_{i=1}^n \lambda_i^3 = \operatorname{tr}(A^3) = 6 \times (\text{number of triangles})$.
    4. $\lambda_1 \ge \bar{d}$ (average degree), with equality for regular graphs.
    5. $G$ is bipartite $\iff$ the spectrum is symmetric about $0$ ($\lambda_i = -\lambda_{n+1-i}$).

??? proof "Proof"
    (1)-(3): $\operatorname{tr}(A^k) = \sum_i \lambda_i^k$, and $(A^k)_{ii}$ counts the number of closed walks of length $k$ starting and ending at $i$. For $k = 1$: no self-loops so $(A)_{ii} = 0$. For $k = 2$: $(A^2)_{ii} = d_i$, summing gives $2|E|$. For $k = 3$: $(A^3)_{ii}$ counts twice the number of triangles through $i$, and the total sum is $6 \times$ the number of triangles.

    (4): $\lambda_1 = \max_{\|\mathbf{x}\| = 1} \mathbf{x}^T A \mathbf{x} \ge \frac{1}{n}\mathbf{1}^T A \mathbf{1} = \frac{1}{n}\sum_i d_i = \bar{d}$.

    (5): If $G$ is bipartite ($V = V_1 \cup V_2$), let $S$ be the diagonal matrix with $S_{ii} = +1$ ($i \in V_1$), $S_{ii} = -1$ ($i \in V_2$). Then $SAS = -A$, so $\lambda$ being an eigenvalue implies $-\lambda$ is also an eigenvalue. $\blacksquare$

!!! example "Example 27.2"
    **Spectra of the complete graph and Petersen graph.** For the complete graph $K_n$: $A = J - I$ ($J$ is the all-ones matrix). The eigenvalues of $J$ are $n$ (multiplicity $1$) and $0$ (multiplicity $n-1$), so the spectrum of $A$ is $n - 1$ (multiplicity $1$) and $-1$ (multiplicity $n-1$).

    The Petersen graph ($10$-vertex $3$-regular graph) has adjacency spectrum $3, 1, 1, 1, 1, 1, -2, -2, -2, -2$. $\lambda_1 = 3 = d$ (regular), and $\lambda_2 = 1$ is relatively small, indicating good expansion properties.

---

## 27.3 Laplacian Matrix and Connectivity

<div class="context-flow" markdown>

**Core result**: $\mu_2 > 0$ $\iff$ graph is connected → $\mu_2$ (Fiedler value / algebraic connectivity) measures the "strength" of connectivity → Kirchhoff's theorem: number of spanning trees $= \frac{1}{n}\mu_2\cdots\mu_n$
**Links**: Core application of determinant theory from Ch3 in Kirchhoff's theorem

</div>

The spectrum of the Laplacian matrix completely characterizes graph connectivity.

!!! definition "Definition 27.5 (Algebraic Connectivity)"
    The second smallest eigenvalue $\mu_2$ of $L$ is called the **algebraic connectivity** or **Fiedler value** of the graph. The corresponding eigenvector is called the **Fiedler vector**.

!!! theorem "Theorem 27.3 (Laplacian and Connectivity)"
    Let the eigenvalues of $L$ be $0 = \mu_1 \le \mu_2 \le \cdots \le \mu_n$:

    1. The multiplicity of $\mu_1 = 0$ equals the number of connected components $k$.
    2. $G$ is connected $\iff$ $\mu_2 > 0$.
    3. For a connected graph, $\mu_2 \le \frac{n}{n-1}\min_i d_i$.

??? proof "Proof"
    The solution space of $L\mathbf{x} = \mathbf{0}$ is spanned by the indicator vectors of the connected components. If $G$ has $k$ connected components, define $\mathbf{x}^{(j)}$ as the indicator vector of the $j$-th component, then $L\mathbf{x}^{(j)} = \mathbf{0}$ (both endpoints of every edge lie in the same component, so $(x_i - x_j)^2 = 0$). These $k$ vectors are linearly independent, so the dimension of $\ker(L)$ is $k$, i.e., the multiplicity of the zero eigenvalue is $k$.

    In particular, $k = 1$ (connected) if and only if $\mu_2 > 0$. $\blacksquare$

!!! theorem "Theorem 27.4 (Kirchhoff's Matrix Tree Theorem)"
    Let $G$ be a connected graph with Laplacian $L$. The number of spanning trees $\tau(G)$ is

    $$
    \tau(G) = \frac{1}{n}\mu_2 \mu_3 \cdots \mu_n = \frac{1}{n}\prod_{i=2}^{n} \mu_i.
    $$

    Equivalently, $\tau(G)$ equals any $(n-1) \times (n-1)$ cofactor of $L$.

??? proof "Proof"
    Define $L_i$ as the matrix obtained by deleting the $i$-th row and column from $L$. Kirchhoff's theorem states that $\tau(G) = \det(L_i)$ (for any $i$).

    From the eigenvalues of $L$, $\det(L_i)$ can be derived via the combinatorial proof of the matrix-tree theorem. On the other hand, note that

    $$
    \det(\lambda I - L) = \lambda \prod_{i=2}^n (\lambda - \mu_i),
    $$

    and the coefficient of $\lambda$ in $\det(\lambda I - L)$ relates to the sum of cofactors. Using the characteristic polynomial expansion and the Cauchy-Binet formula, one proves that $\det(L_i) = \frac{1}{n}\prod_{i=2}^n \mu_i$. $\blacksquare$

!!! example "Example 27.3"
    **Counting spanning trees of the complete graph.** The Laplacian of $K_n$ is $L = nI - J$, with eigenvalues $0$ (multiplicity $1$) and $n$ (multiplicity $n-1$). By Kirchhoff's theorem,

    $$
    \tau(K_n) = \frac{1}{n} \cdot n^{n-1} = n^{n-2}.
    $$

    This is the celebrated **Cayley's formula**. For example, $\tau(K_4) = 4^2 = 16$.

---

## 27.4 Cheeger Inequality and Graph Partitioning

<div class="context-flow" markdown>

**Core inequality**: $\frac{h^2}{2d_{\max}} \le \mu_2 \le 2h$ → the Fiedler value $\mu_2$ is equivalent to the Cheeger constant $h$ → spectral clustering algorithm: use the Fiedler vector to partition vertices into two groups
**Applications**: Image segmentation, community detection, data clustering

</div>

The Cheeger inequality connects an algebraic quantity ($\mu_2$) to a combinatorial quantity (optimal graph cut), providing the theoretical foundation for spectral clustering.

!!! definition "Definition 27.6 (Isoperimetric Constant / Cheeger Constant)"
    The **Cheeger constant** (or isoperimetric constant / conductance) of a graph $G$ is defined as

    $$
    h(G) = \min_{S \subset V,\, 0 < |S| \le n/2} \frac{|\partial(S)|}{\min(|S|, |V \setminus S|)},
    $$

    where $\partial(S) = \{\{u, v\} \in E : u \in S, v \notin S\}$ is the edge boundary of $S$. $h(G)$ measures the minimum cost of cutting the graph into two parts.

!!! theorem "Theorem 27.5 (Discrete Cheeger Inequality)"
    For a $d$-regular graph $G$, let $\mu_2$ be the second smallest eigenvalue of the normalized Laplacian $\mathcal{L}$ and $h$ be the Cheeger constant. Then

    $$
    \frac{h^2}{2} \le \mu_2 \le 2h.
    $$

    For general graphs (using $L$ and an appropriately defined $h$), analogous inequalities hold.

??? proof "Proof"
    **Upper bound** $\mu_2 \le 2h$: Take $S$ to be the set achieving $h$, and construct a test vector $\mathbf{x}$ with $x_i = |V \setminus S|$ ($i \in S$), $x_i = -|S|$ ($i \notin S$), so that $\mathbf{x} \perp \mathbf{1}$. By the Rayleigh quotient,

    $$
    \mu_2 \le \frac{\mathbf{x}^T L \mathbf{x}}{\mathbf{x}^T D \mathbf{x}} = \frac{n^2 |\partial(S)|}{d \cdot |S| \cdot |V \setminus S|} \le \frac{2|\partial(S)|}{d \cdot |S|} = \frac{2h}{d} \cdot d = 2h.
    $$

    **Lower bound** $h^2/2 \le \mu_2$: Let $\mathbf{f}$ be the eigenvector corresponding to $\mu_2$. After sorting the components of $\mathbf{f}$, a level-set analysis combined with the Cauchy-Schwarz inequality establishes the lower bound. $\blacksquare$

!!! example "Example 27.4"
    **The spectral clustering algorithm.** Given a graph $G$ (or a similarity graph of data points), spectral clustering proceeds as follows:

    1. Compute the Laplacian $L = D - A$ (or normalized Laplacian $\mathcal{L}$).
    2. Find the $k$ smallest eigenvalues of $\mathcal{L}$ and their eigenvectors $\mathbf{u}_1, \ldots, \mathbf{u}_k$.
    3. Form the matrix $U \in \mathbb{R}^{n \times k}$ and normalize each row.
    4. Apply $k$-means clustering to the row vectors of $U$.

    The Fiedler vector (eigenvector of $\mu_2$) naturally partitions vertices into two groups based on the signs of its components. For social network data, spectral clustering effectively discovers community structure.

---

## 27.5 Random Walks and PageRank

<div class="context-flow" markdown>

**Chain**: Random walk on graph → transition matrix $P = D^{-1}A$ → Perron-Frobenius theorem (Ch17) → stationary distribution → Google PageRank = stationary distribution of a damped random walk
**Links**: Direct application of nonnegative matrix theory from Ch17

</div>

Random walks on graphs transform graph-theoretic problems into stochastic matrix analysis, and PageRank is their most famous application.

!!! definition "Definition 27.7 (Random Walk on a Graph)"
    On an undirected graph $G$, a **simple random walk** moves at each step from the current vertex $i$ to a uniformly random neighbor $j$. The transition probability matrix is

    $$
    P = D^{-1}A, \quad P_{ij} = \frac{A_{ij}}{d_i}.
    $$

    $P$ is a row-stochastic matrix (each row sums to $1$).

!!! definition "Definition 27.8 (PageRank)"
    For a directed graph $G$, the **PageRank** vector $\boldsymbol{\pi}$ is the stationary distribution of a modified random walk. The Google matrix is defined as

    $$
    M = \alpha P + (1 - \alpha) \frac{1}{n} \mathbf{1}\mathbf{1}^T,
    $$

    where $\alpha \in (0, 1)$ is the damping factor (typically $\alpha = 0.85$) and $P$ is the column-stochastic link matrix. $\boldsymbol{\pi}$ satisfies $M\boldsymbol{\pi} = \boldsymbol{\pi}$, $\boldsymbol{\pi}^T \mathbf{1} = 1$.

!!! theorem "Theorem 27.6 (Existence and Uniqueness of PageRank)"
    For $0 < \alpha < 1$, the Google matrix $M$ is a strictly positive matrix (all entries $> 0$), and therefore:

    1. $M$ has a unique eigenvalue $1$ (Perron-Frobenius theorem), with corresponding eigenvector $\boldsymbol{\pi} > \mathbf{0}$ (all components positive).
    2. Power iteration $\boldsymbol{\pi}^{(k+1)} = M\boldsymbol{\pi}^{(k)}$ converges geometrically to $\boldsymbol{\pi}$ at rate $\alpha$.
    3. The second largest eigenvalue modulus $|\lambda_2| \le \alpha < 1$.

??? proof "Proof"
    Each entry of $M$ is at least $(1-\alpha)/n > 0$, so $M$ is a strictly positive matrix. By the Perron-Frobenius theorem (Ch17), $M$ has a unique largest eigenvalue $1$ (since $M$ is column-stochastic or row-stochastic), with a corresponding positive eigenvector.

    The convergence rate is controlled by the spectral gap $1 - |\lambda_2|$. Since $M = \alpha P + (1-\alpha)J/n$, where $P$ has spectral radius $1$ and $J/n$ has unique nonzero eigenvalue $1$ (multiplicity $1$), one can show that $|\lambda_2(M)| \le \alpha$. $\blacksquare$

!!! example "Example 27.5"
    **Computing PageRank for a simple network.** Consider a three-page network: page 1 links to 2 and 3, page 2 links to 3, page 3 links to 1. The column-stochastic matrix is

    $$
    P = \begin{pmatrix} 0 & 0 & 1 \\ 1/2 & 0 & 0 \\ 1/2 & 1 & 0 \end{pmatrix}.
    $$

    With $\alpha = 0.85$, $M = 0.85P + 0.05 \cdot \mathbf{1}\mathbf{1}^T$. Power iteration starting from $\boldsymbol{\pi}^{(0)} = (1/3, 1/3, 1/3)^T$ converges to $\boldsymbol{\pi} \approx (0.387, 0.214, 0.399)^T$. Page 3 has the highest PageRank (linked by pages 1 and 2), followed by page 1 (linked by page 3, which has high weight).

---

## 27.6 Expander Graphs

<div class="context-flow" markdown>

**Definition**: $d$-regular graph with $\lambda_2(A) \le d - \varepsilon$ → good expander = large spectral gap → Alon-Boppana bound: $\lambda_2 \ge 2\sqrt{d-1} - o(1)$ → Ramanujan graphs achieve this bound
**Applications**: Error-correcting codes, derandomization, network design

</div>

Expander graphs are sparse graphs with good connectivity and pseudorandom properties, whose definition and analysis essentially depend on spectral theory.

!!! definition "Definition 27.9 (Spectral Expander)"
    A $d$-regular graph $G$ is called an $(n, d, \lambda)$**-graph** if the second largest eigenvalue (in absolute value) $\lambda = \max(|\lambda_2|, |\lambda_n|) \le \lambda$. The smaller $\lambda$ is, the better the expansion. The **spectral gap** is defined as $d - \lambda$.

!!! theorem "Theorem 27.7 (Expander Mixing Lemma)"
    For an $(n, d, \lambda)$-graph $G$, for any two vertex subsets $S, T \subseteq V$, the edge count satisfies

    $$
    \left| e(S, T) - \frac{d \cdot |S| \cdot |T|}{n} \right| \le \lambda \sqrt{|S| \cdot |T|},
    $$

    where $e(S, T)$ is the number of edges between $S$ and $T$.

??? proof "Proof"
    Let $\mathbf{x} = \mathbf{1}_S$, $\mathbf{y} = \mathbf{1}_T$. Then $e(S, T) = \mathbf{x}^T A \mathbf{y}$. Decompose $\mathbf{x}$ and $\mathbf{y}$ into eigenvectors of $A$: $\mathbf{x} = \sum \hat{x}_i \mathbf{v}_i$, $\mathbf{y} = \sum \hat{y}_i \mathbf{v}_i$.

    $$
    \mathbf{x}^T A \mathbf{y} = \sum_i \lambda_i \hat{x}_i \hat{y}_i = d \hat{x}_1 \hat{y}_1 + \sum_{i \ge 2} \lambda_i \hat{x}_i \hat{y}_i.
    $$

    $\hat{x}_1 = \langle \mathbf{x}, \frac{\mathbf{1}}{\sqrt{n}} \rangle = |S|/\sqrt{n}$, similarly $\hat{y}_1 = |T|/\sqrt{n}$. The first term is $d|S||T|/n$. The second term is controlled by Cauchy-Schwarz:

    $$
    \left|\sum_{i \ge 2} \lambda_i \hat{x}_i \hat{y}_i\right| \le \lambda \sqrt{\sum_{i \ge 2} \hat{x}_i^2} \sqrt{\sum_{i \ge 2} \hat{y}_i^2} \le \lambda \|\mathbf{x}\| \|\mathbf{y}\| = \lambda\sqrt{|S| \cdot |T|}.
    $$

    $\blacksquare$

!!! theorem "Theorem 27.8 (Alon-Boppana Bound)"
    For an infinite family of $d$-regular graphs $\{G_n\}$ ($|V(G_n)| \to \infty$),

    $$
    \liminf_{n \to \infty} \lambda_2(G_n) \ge 2\sqrt{d - 1}.
    $$

    Graphs achieving this bound ($\lambda \le 2\sqrt{d-1}$) are called **Ramanujan graphs**.

??? proof "Proof"
    (Proof sketch) The argument uses the spectrum of the infinite $d$-regular tree $T_d$. The spectrum of the adjacency operator of $T_d$ is $[-2\sqrt{d-1}, 2\sqrt{d-1}]$. Finite $d$-regular graphs can "locally" approximate $T_d$ (when they have high girth), so their nontrivial eigenvalues cannot all be far from this interval. The precise proof uses the trace method and combinatorial counting. $\blacksquare$

!!! example "Example 27.6"
    **Construction of Ramanujan graphs.** The Lubotzky-Phillips-Sarnak (LPS) graphs are a classical construction of Ramanujan graphs. For primes $p \equiv 1 \pmod{4}$ and $q$ ($q \ne p$), the LPS graph is a $(p+1)$-regular Cayley graph with vertex set $\text{PSL}(2, \mathbb{F}_q)$.

    $\lambda_2 \le 2\sqrt{p}$, satisfying the Ramanujan bound. These graphs have $O(q^3)$ vertices, are $(p+1)$-regular, and have spectral gap $p + 1 - 2\sqrt{p} = (\sqrt{p} - 1)^2$, making them nearly optimal sparse expanders.

---

## 27.7 Graph Coloring and Eigenvalues

<div class="context-flow" markdown>

**Bounds**: $\chi(G) \ge 1 + \lambda_1 / |\lambda_n|$ (Hoffman bound) → eigenvalues give lower bounds on chromatic number → independence number upper bound: $\alpha(G) \le n \cdot |\lambda_n| / (\lambda_1 + |\lambda_n|)$
**Applications**: Relaxation bounds in combinatorial optimization

</div>

Eigenvalues of a graph can provide effective bounds on its chromatic number and independence number.

!!! definition "Definition 27.10 (Chromatic Number)"
    The **chromatic number** $\chi(G)$ of a graph $G$ is the minimum number of colors needed for a proper coloring (adjacent vertices receive different colors).

!!! theorem "Theorem 27.9 (Hoffman Chromatic Number Bound)"
    For a nonempty graph $G$,

    $$
    \chi(G) \ge 1 + \frac{\lambda_1}{|\lambda_n|} = 1 - \frac{\lambda_1}{\lambda_n},
    $$

    where $\lambda_1$ and $\lambda_n$ are the largest and smallest eigenvalues of $A$, respectively.

??? proof "Proof"
    Suppose $G$ has a proper $k$-coloring with color classes $C_1, \ldots, C_k$ (independent sets). For each color class $C_s$, the induced subgraph has no edges, so

    $$
    \sum_{i, j \in C_s} A_{ij} = 0, \quad \forall s.
    $$

    Let $\mathbf{x}^{(s)}$ be the (appropriately centered) indicator vector of $C_s$. Using the quadratic form of $A$ and the Rayleigh quotient, one can show that $\lambda_n \le -\lambda_1 / (k-1)$, i.e., $k \ge 1 + \lambda_1/|\lambda_n|$. $\blacksquare$

!!! theorem "Theorem 27.10 (Hoffman Independence Number Bound)"
    For a $d$-regular graph $G$, the maximum independent set size satisfies

    $$
    \alpha(G) \le \frac{n \cdot |\lambda_n|}{\lambda_1 + |\lambda_n|} = \frac{n \cdot |\lambda_n|}{d + |\lambda_n|}.
    $$

??? proof "Proof"
    Let $S$ be an independent set with $|S| = \alpha$. Set $\mathbf{x} = \mathbf{1}_S - \frac{\alpha}{n}\mathbf{1}$ ($\mathbf{x} \perp \mathbf{1}$). Since $S$ is independent, $\mathbf{1}_S^T A \mathbf{1}_S = 0$.

    $$
    \mathbf{x}^T A \mathbf{x} = -\frac{2\alpha}{n}\mathbf{1}_S^T A \mathbf{1} + \frac{\alpha^2}{n^2}\mathbf{1}^T A \mathbf{1} = -\frac{2\alpha d\alpha}{n} + \frac{\alpha^2 dn}{n^2} = -\frac{\alpha^2 d}{n}.
    $$

    Since $\mathbf{x} \perp \mathbf{1}$, $\mathbf{x}^T A \mathbf{x} \ge \lambda_n \|\mathbf{x}\|^2 = \lambda_n (\alpha - \alpha^2/n)$. Therefore

    $$
    -\frac{\alpha^2 d}{n} \ge \lambda_n \alpha (1 - \alpha/n),
    $$

    which simplifies to $\alpha \le n|\lambda_n| / (d + |\lambda_n|)$. $\blacksquare$

!!! example "Example 27.7"
    **Chromatic number of Kneser graphs.** The Kneser graph $K(n, k)$ has as vertices all $k$-element subsets of $\{1, \ldots, n\}$, with two vertices adjacent if and only if the corresponding subsets are disjoint.

    $K(5, 2)$ is the Petersen graph, with spectrum $3, 1^5, (-2)^4$. The Hoffman bound gives $\chi \ge 1 + 3/2 = 2.5$, i.e., $\chi \ge 3$. In fact $\chi(K(5, 2)) = 3$ (the Kneser conjecture proved by Lovasz), and the spectral bound is tight here.

---

## 27.8 Network Flows and Linear Programming

<div class="context-flow" markdown>

**Model**: Network flow = linearly constrained optimization on graphs → max-flow/min-cut = special case of LP duality → the incidence matrix $B$ plays a central role in flow conservation constraints
**Links**: LP theory from Ch25 made concrete on graphs

</div>

Network flow problems represent a classical intersection of graph theory and linear programming.

!!! definition "Definition 27.11 (Network Flow)"
    The **network flow problem**: Given a directed graph $G = (V, E)$, source $s$, sink $t$, and capacity function $c : E \to \mathbb{R}_{\ge 0}$. A **flow** $\mathbf{f} \in \mathbb{R}^m$ satisfies:

    1. **Capacity constraints**: $0 \le f_e \le c_e$, $\forall e \in E$.
    2. **Flow conservation**: $\sum_{e \in \delta^+(v)} f_e = \sum_{e \in \delta^-(v)} f_e$, $\forall v \neq s, t$.

    Using the incidence matrix $B$, flow conservation is written as $B\mathbf{f} = \mathbf{b}$ ($b_s = -F$, $b_t = F$, rest are $0$).

!!! theorem "Theorem 27.11 (Max-Flow Min-Cut Theorem)"
    The maximum flow value equals the minimum cut capacity:

    $$
    \max_{\mathbf{f}} F = \min_{S:\, s \in S,\, t \notin S} \sum_{e \in \delta^+(S)} c_e.
    $$

    This is the manifestation of LP duality on networks: the dual of the max-flow LP is precisely the min-cut LP.

??? proof "Proof"
    **Weak duality** ($\le$): For any flow $\mathbf{f}$ and cut $(S, \bar{S})$, $F = \sum_{e \in \delta^+(S)} f_e - \sum_{e \in \delta^-(S)} f_e \le \sum_{e \in \delta^+(S)} c_e$.

    **Strong duality** ($=$): By the LP strong duality theorem, or by the termination of the Ford-Fulkerson augmenting path algorithm. When no augmenting path from $s$ to $t$ exists, let $S$ be the set of vertices reachable from $s$ in the residual graph; then $(S, \bar{S})$ is a minimum cut, and the current flow is maximum. $\blacksquare$

!!! theorem "Theorem 27.12 (Total Unimodularity and Integer Flows)"
    The incidence matrix $B$ of a directed graph is **totally unimodular**: every square submatrix of $B$ has determinant $0$, $+1$, or $-1$. Therefore, when capacities $c_e$ are integers, the optimal solution of the max-flow LP is automatically integer.

??? proof "Proof"
    Proceed by induction on the size $k$ of a $k \times k$ submatrix $B'$ of $B$. For $k = 1$: entries of $B'$ are $0, \pm 1$. For $k > 1$: if some column of $B'$ is all zeros, the determinant is $0$. If some column has exactly one nonzero entry, expand along that column to reduce to a $(k-1)$-order submatrix. If some column has two nonzero entries ($+1$ and $-1$), the column sum is zero, and linear dependence of rows can be used for simplification, then the induction proceeds.

    Total unimodularity ensures that $B\mathbf{f} = \mathbf{b}$ with integer $\mathbf{b}$ has an integer basic feasible solution, so the LP optimum is integer. $\blacksquare$

!!! example "Example 27.8"
    **LP formulation of maximum flow.** A network has $4$ vertices $\{s, a, b, t\}$ with edges and capacities $s \to a$ (3), $s \to b$ (2), $a \to b$ (1), $a \to t$ (2), $b \to t$ (3). The LP formulation:

    $$
    \max F, \quad \text{s.t.} \quad B\mathbf{f} = F(\mathbf{e}_t - \mathbf{e}_s), \quad \mathbf{0} \le \mathbf{f} \le \mathbf{c}.
    $$

    Maximum flow $F = 4$ ($s \to a$: 3, $s \to b$: 1, $a \to b$: 1, $a \to t$: 2, $b \to t$: 2). The minimum cut is $\{s, a\}$, with cut capacity $= 1 + 2 + 1 = 4$. The zero duality gap confirms strong duality. Total unimodularity of the incidence matrix guarantees the existence of integer flows.
