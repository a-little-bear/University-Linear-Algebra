# 第 41A 章 矩阵束与正则束理论

<div class="context-flow" markdown>

**前置**：特征值与特征向量(Ch6) · Jordan 标准形(Ch12) · $\lambda$-矩阵与 Smith 标准形(Ch13B) · Schur 分解(Ch10) · Hermite 矩阵(Ch7)

**本章脉络**：矩阵束定义 $\to$ 正则束与奇异束分类 $\to$ 广义特征值（齐次坐标） $\to$ 严格等价 $\to$ Weierstrass 标准形 $\to$ Smith 形 $\to$ 收缩子空间 $\to$ 广义 Schur 分解（QZ 分解） $\to$ QZ 算法 $\to$ 对称正定束 $\to$ Hermite 束 $\to$ 扰动理论 $\to$ 多项式特征值问题与线性化 $\to$ 保结构线性化

**延伸**：正则束理论是 LAPACK 中广义特征值求解器（`xGGES`、`xGGEV`）的数学基础；QZ 算法是 MATLAB `eig(A,B)` 的核心；对称正定束在振动分析和有限元法中不可或缺；多项式特征值问题及其线性化在结构力学（阻尼振动）、流体-结构耦合、光子晶体等领域广泛出现

</div>

标准特征值问题 $Ax = \lambda x$ 研究的是矩阵束 $A - \lambda I$ 的结构。当 $I$ 替换为一般矩阵 $B$ 时，便得到**广义特征值问题** $Ax = \lambda Bx$ 和**矩阵束** $A - \lambda B$。若 $\det(A - \lambda B)$ 作为 $\lambda$ 的多项式不恒为零，则称束为**正则的**；此时存在 Weierstrass 标准形——一种精巧地将有限特征值（Jordan 块）与无穷特征值（幂零块）分离的标准形。

正则束理论不仅是矩阵论中最典雅的结构之一，也是数值线性代数中不可替代的工具。本章从矩阵束的基本概念出发，建立 Weierstrass 标准形和 Smith 标准形，引入收缩子空间的概念，发展广义 Schur 分解（QZ 分解）及其数值实现——QZ 算法，讨论具有特殊对称结构的束的性质，最后研究多项式特征值问题的线性化方法。

---

## 41A.1 矩阵束的定义

<div class="context-flow" markdown>

**核心问题**：什么是矩阵束？如何区分正则束与奇异束？

</div>

!!! definition "定义 41A.1 (矩阵束)"
    设 $A, B \in M_{m \times n}(\mathbb{C})$。**矩阵束**（matrix pencil）$L(\lambda)$ 定义为
    $$L(\lambda) = A - \lambda B,$$
    其中 $\lambda \in \mathbb{C}$ 为参数。当 $m = n$（方阵情形）时，$\det(A - \lambda B)$ 是 $\lambda$ 的至多 $n$ 次多项式（也可能恒为零）。

    更一般地，也常将矩阵束写为 $\lambda B - A$ 或 $\alpha A + \beta B$（齐次形式）。本章主要采用 $A - \lambda B$ 的记法。

!!! definition "定义 41A.2 (正则束与奇异束)"
    设 $A, B \in M_n(\mathbb{C})$（方阵）。矩阵束 $A - \lambda B$ 称为**正则的**（regular），若
    $$\det(A - \lambda B) \not\equiv 0,$$
    即 $\det(A - \lambda B)$ 作为 $\lambda$ 的多项式不恒等于零。否则称为**奇异的**（singular）。

    对于长方矩阵束（$m \ne n$），束总是奇异的，因为不存在方阵行列式。

正则性条件等价于：存在至少一个 $\lambda_0 \in \mathbb{C}$ 使得 $A - \lambda_0 B$ 非奇异。

!!! example "例 41A.1"
    **(a)** $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$，$B = I$：$\det(A - \lambda I) = (1-\lambda)(2-\lambda)$。这是标准特征值问题，正则束。

    **(b)** $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$：$\det(A - \lambda B) = \det\begin{pmatrix} 1 & 0 \\ 0 & -\lambda \end{pmatrix} = -\lambda$。正则束，有一个有限特征值 $\lambda = 0$ 和一个无穷特征值。

    **(c)** $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$：$\det(A - \lambda B) = 0$ 对所有 $\lambda$。奇异束。

    **(d)** $A = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 1 \\ 0 & 0 \\ 0 & 0 \end{pmatrix}$：$3 \times 2$ 长方束，自动为奇异束。

---

## 41A.2 广义特征值

<div class="context-flow" markdown>

**核心问题**：如何定义广义特征值？齐次坐标为何必要？

</div>

!!! definition "定义 41A.3 (广义特征值——仿射定义)"
    设 $A - \lambda B$ 为 $n \times n$ 正则束。$\lambda_0 \in \mathbb{C}$ 称为 $A - \lambda B$ 的**有限广义特征值**，若
    $$\det(A - \lambda_0 B) = 0.$$
    若 $\det(B) = 0$ 但 $\det(A - \lambda B) \not\equiv 0$，则称束具有**无穷广义特征值**。

仿射定义的不足在于无穷特征值需要特殊处理，且重数的定义不够自然。齐次坐标提供了统一的框架。

!!! definition "定义 41A.4 (广义特征值——齐次坐标)"
    设 $A, B \in M_n(\mathbb{C})$，$A - \lambda B$ 为正则束。非零对 $(\alpha, \beta) \in \mathbb{C}^2 \setminus \{(0,0)\}$ 称为**广义特征对**，若
    $$\det(\beta A - \alpha B) = 0.$$
    两个特征对 $(\alpha_1, \beta_1)$ 和 $(\alpha_2, \beta_2)$ 等价当且仅当 $(\alpha_1, \beta_1) = c(\alpha_2, \beta_2)$ 对某 $c \ne 0$，即它们在 $\mathbb{P}^1(\mathbb{C})$（复射影直线）中代表同一点。

    - **有限特征值**：$\beta \ne 0$ 时，$\lambda = \alpha/\beta$。
    - **无穷特征值**：$\beta = 0$，$\alpha \ne 0$ 时，对应 $\lambda = \infty$。

!!! theorem "定理 41A.1 (广义特征值的个数)"
    设 $A - \lambda B$ 为 $n \times n$ 正则束，$d = \deg(\det(A - \lambda B))$。则

    (a) 束有 $d$ 个有限广义特征值（计代数重数）；

    (b) 束有 $n - d$ 个无穷广义特征值（计代数重数）；

    (c) 有限与无穷广义特征值总数恰为 $n$。

??? proof "证明"
    $\det(A - \lambda B)$ 是 $\lambda$ 的至多 $n$ 次多项式。设其次数为 $d \le n$，则 $\det(A - \lambda B)$ 有恰好 $d$ 个根（计重数）。

    为计无穷特征值重数，转至齐次形式。令 $p(\alpha, \beta) = \det(\beta A - \alpha B)$，这是 $(\alpha, \beta)$ 的 $n$ 次齐次多项式。因式分解为
    $$p(\alpha, \beta) = c \prod_{i=1}^{n}(\beta_i \alpha - \alpha_i \beta),$$
    其中 $(\alpha_i, \beta_i)$ 为特征对。有 $d$ 个满足 $\beta_i \ne 0$（有限特征值），余下 $n - d$ 个满足 $\beta_i = 0$（无穷特征值）。

!!! example "例 41A.2"
    $A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$。

    $\det(A - \lambda B) = (1 - \lambda)(2 - \lambda)(3) = 3(1-\lambda)(2-\lambda)$，$d = 2$。

    有限特征值：$\lambda = 1, 2$。无穷特征值个数：$3 - 2 = 1$。

    齐次形式：$\det(\beta A - \alpha B) = \beta^3 \cdot 3 \cdot (1 - \alpha/\beta)(2 - \alpha/\beta) \cdot 1 = 3(\beta - \alpha)(2\beta - \alpha)\beta$。根为 $(\alpha, \beta) = (1,1), (2,1), (1,0)$。

---

## 41A.3 严格等价

<div class="context-flow" markdown>

**核心问题**：矩阵束的等价关系如何定义？它与矩阵相似有何异同？

</div>

!!! definition "定义 41A.5 (严格等价)"
    两个 $m \times n$ 矩阵束 $A_1 - \lambda B_1$ 和 $A_2 - \lambda B_2$ 称为**严格等价的**（strictly equivalent），若存在非奇异矩阵 $P \in GL_m(\mathbb{C})$，$Q \in GL_n(\mathbb{C})$，使得
    $$P(A_1 - \lambda B_1)Q = A_2 - \lambda B_2,$$
    即同时有 $PA_1Q = A_2$ 且 $PB_1Q = B_2$。记为 $(A_1, B_1) \sim_s (A_2, B_2)$。

严格等价是一种比矩阵相似更一般的关系。矩阵相似 $A_2 = P^{-1}A_1 P$ 是 $B = I$、$Q = P^{-1}$ 时的特殊情况。

!!! theorem "定理 41A.2 (严格等价保持的不变量)"
    若 $(A_1, B_1) \sim_s (A_2, B_2)$，则：

    (a) $A_1 - \lambda B_1$ 与 $A_2 - \lambda B_2$ 的正则性/奇异性相同；

    (b) 作为 $\mathbb{C}[\lambda]$ 上的矩阵，$A_1 - \lambda B_1$ 与 $A_2 - \lambda B_2$ 有相同的不变因子（Smith 标准形相同）；

    (c) 二者有相同的广义特征值（含重数）；

    (d) 二者有相同的 Kronecker 不变量（最小指标、Jordan 结构等）。

??? proof "证明"
    (a) $\det(A_2 - \lambda B_2) = \det(P) \det(A_1 - \lambda B_1) \det(Q)$。由于 $\det(P), \det(Q) \ne 0$，$\det(A_2 - \lambda B_2) \equiv 0$ 当且仅当 $\det(A_1 - \lambda B_1) \equiv 0$。

    (b) 在 $\mathbb{C}[\lambda]$ 上，$P$ 和 $Q$ 是常数可逆矩阵（不依赖于 $\lambda$），因此是初等行列变换的特殊情况，保持 Smith 标准形不变。

    (c) 由 (a) 的论证，$\det(A_2 - \lambda B_2) = \det(P)\det(Q) \cdot \det(A_1 - \lambda B_1)$，二者有相同的零点。

    (d) 这是 Kronecker 标准形定理的唯一性部分（见第 41B 章）。

!!! definition "定义 41A.6 (合同等价)"
    当需要保持对称或 Hermite 结构时，使用**合同等价**：$A_2 - \lambda B_2 = Q^T(A_1 - \lambda B_1)Q$（对称情况）或 $Q^*(A_1 - \lambda B_1)Q$（Hermite 情况）。合同等价是严格等价的特殊情况（$P = Q^{-T}$ 或 $P = Q^{-*}$）。

---

## 41A.4 正则束的 Weierstrass 标准形

<div class="context-flow" markdown>

**核心问题**：正则束在严格等价下的完全标准形是什么？

</div>

!!! theorem "定理 41A.3 (Weierstrass 标准形)"
    设 $A - \lambda B$ 为 $n \times n$ 正则束。则存在非奇异矩阵 $P, Q \in GL_n(\mathbb{C})$ 使得
    $$P(A - \lambda B)Q = \begin{pmatrix} J - \lambda I_r & 0 \\ 0 & I_s - \lambda N \end{pmatrix},$$
    其中：

    - $J \in M_r(\mathbb{C})$ 为 Jordan 标准形，其对角块 $J_{k_i}(\lambda_i)$ 对应有限广义特征值 $\lambda_1, \ldots, \lambda_p$；
    - $N \in M_s(\mathbb{C})$ 为**幂零** Jordan 矩阵（所有特征值为 $0$），$N = \operatorname{diag}(N_{l_1}, \ldots, N_{l_t})$，对应无穷广义特征值；
    - $r + s = n$，$r$ 为 $\det(A - \lambda B)$ 的次数。

    此标准形在块的排列顺序外是**唯一的**。

??? proof "证明"
    **步骤一（正则性的利用）**：由 $A - \lambda B$ 正则，$\det(A - \lambda B)$ 不恒为零。设 $d = \deg(\det(A - \lambda B)) = r \le n$。则存在至少一个 $\mu \in \mathbb{C}$ 使得 $\det(A - \mu B) \ne 0$，即 $A - \mu B$ 非奇异。

    **步骤二（化为标准特征值问题）**：令 $C = (A - \mu B)^{-1}B$。则广义特征值问题 $Ax = \lambda Bx$ 等价于
    $$(A - \mu B)^{-1}(A - \lambda B)x = (A - \mu B)^{-1}(A - \mu B + \mu B - \lambda B)x = (I - (\lambda - \mu)C)x = 0.$$
    令 $\lambda' = 1/(\lambda - \mu)$（当 $\lambda \ne \mu$ 时），则 $Cx = \lambda' x$ 是标准特征值问题。

    $C$ 的特征值 $\lambda'$ 与广义特征值 $\lambda$ 的关系：
    - $\lambda' = 1/(\lambda - \mu) \ne 0$：对应有限广义特征值 $\lambda = \mu + 1/\lambda'$；
    - $\lambda' = 0$：对应无穷广义特征值。

    **步骤三（Jordan 分解）**：对 $C$ 做 Jordan 分解 $C = S \begin{pmatrix} J_C & 0 \\ 0 & N_C \end{pmatrix} S^{-1}$，其中 $J_C$ 的对角块对应 $C$ 的非零特征值（有限广义特征值），$N_C$ 为幂零部分（无穷广义特征值），大小分别为 $r \times r$ 和 $s \times s$。

    **步骤四（构造变换矩阵）**：取 $Q = S$，并令 $P = (A - \mu B)^{-1} Q^{-1}$ 的适当调整。则
    $$P(A - \lambda B)Q = P(A - \mu B)Q \cdot Q^{-1}(A-\mu B)^{-1}(A - \lambda B)Q$$
    $$= \tilde{D} \cdot (I - (\lambda - \mu) \begin{pmatrix} J_C & 0 \\ 0 & N_C \end{pmatrix}).$$
    通过对各块的仿射变换 $\lambda_i = \mu + 1/\lambda'_i$ 和适当的标度变换，将有限部分化为 $J - \lambda I_r$，将幂零部分化为 $I_s - \lambda N$。

    具体地，对有限部分：$I_r - (\lambda - \mu)J_C = I_r - (\lambda - \mu)J_C$。令 $J = \mu I_r + J_C^{-1}$（对非零 $J_C$ 的各块取逆后做仿射变换），并适当调整 $P, Q$，可以得到 $J - \lambda I_r$ 的形式。

    对无穷部分：$I_s - (\lambda - \mu)N_C = I_s - \lambda N_C + \mu N_C$。由于 $N_C$ 幂零，$I_s + \mu N_C$ 可逆。左乘 $(I_s + \mu N_C)^{-1}$ 后得到 $(I_s + \mu N_C)^{-1} - \lambda (I_s + \mu N_C)^{-1}N_C$。再通过相似变换将 $(I_s + \mu N_C)^{-1}N_C$ 化为 Jordan 形 $N$，同时调整左边为 $I_s$。

    **唯一性**：$J$ 和 $N$ 的 Jordan 结构由不变因子唯一确定。有限特征值的初等因子为 $(\lambda - \lambda_i)^{k_i}$，无穷特征值的初等因子为 $\mu^{l_j}$（在 $\mu = 1/\lambda$ 变换下）。由 Jordan 标准形的唯一性即得 Weierstrass 标准形的唯一性。

!!! example "例 41A.3"
    矩阵束 $A - \lambda B$，$A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$。

    $\det(A - \lambda B) = \det\begin{pmatrix} 1 & 0 & 0 \\ 0 & -\lambda & 1 \\ 0 & 0 & -\lambda \end{pmatrix} = \lambda^2$。

    有限特征值 $\lambda = 0$（代数重数 $2$）。$\deg(\det) = 2 < n = 3$，故有 $1$ 个无穷特征值。

    Weierstrass 标准形：$r = 2$，$s = 1$。

    - $J = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$（$\lambda = 0$ 的 $2 \times 2$ Jordan 块）；
    - $N = (0)$（$1 \times 1$ 幂零块）。

    $$P(A - \lambda B)Q = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix} - \lambda \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}.$$

!!! example "例 41A.4"
    $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。

    $\det(A - \lambda B) = \det\begin{pmatrix} 2 - \lambda & 1 \\ 0 & 2 \end{pmatrix} = 2(2 - \lambda)$。$d = 1 < n = 2$。

    有限特征值 $\lambda = 2$（重数 $1$），无穷特征值 $1$ 个。

    Weierstrass 标准形：$J = (2)$，$N = (0)$。
    $$P(A - \lambda B)Q = \begin{pmatrix} 2 - \lambda & 0 \\ 0 & 1 \end{pmatrix}.$$

---

## 41A.5 正则束的 Smith 标准形

<div class="context-flow" markdown>

**核心问题**：如何通过 $\mathbb{C}[\lambda]$ 上的初等变换得到正则束的标准形？

</div>

!!! theorem "定理 41A.4 (正则束的 Smith 标准形)"
    设 $A - \lambda B$ 为 $n \times n$ 正则束。将 $A - \lambda B$ 视为 $\mathbb{C}[\lambda]$ 上的矩阵，则存在 $\mathbb{C}[\lambda]$ 上的可逆矩阵 $U(\lambda), V(\lambda)$（行列式为非零常数的多项式矩阵）使得
    $$U(\lambda)(A - \lambda B)V(\lambda) = \operatorname{diag}(d_1(\lambda), d_2(\lambda), \ldots, d_n(\lambda)),$$
    其中 $d_i(\lambda)$ 为首一多项式（或为 $1$），且 $d_i \mid d_{i+1}$（$i = 1, \ldots, n-1$）。$d_i(\lambda)$ 称为**不变因子**（invariant factors）。

!!! definition "定义 41A.7 (初等因子)"
    设正则束 $A - \lambda B$ 的不变因子为 $d_1(\lambda), \ldots, d_n(\lambda)$。将每个 $d_i(\lambda)$ 分解为不可约因子的幂：
    $$d_i(\lambda) = \prod_j (\lambda - \lambda_j)^{e_{ij}}.$$
    所有出现的 $(\lambda - \lambda_j)^{e_{ij}}$（$e_{ij} > 0$）称为有限**初等因子**（elementary divisors）。

    无穷初等因子由"反转"得到：令 $\mu = 1/\lambda$，考虑束在 $\mu = 0$ 处的结构。无穷初等因子为 $\mu^{l_j}$，其中 $l_j$ 为 Weierstrass 形中幂零块 $N$ 的各 Jordan 块大小。

!!! theorem "定理 41A.5 (Weierstrass 与 Smith 形的关系)"
    正则束 $A - \lambda B$ 的 Weierstrass 标准形中：

    (a) 有限部分 $J$ 的 Jordan 结构完全由有限初等因子 $\{(\lambda - \lambda_i)^{k_i}\}$ 决定；

    (b) 幂零部分 $N$ 的 Jordan 结构完全由无穷初等因子 $\{\mu^{l_j}\}$ 决定；

    (c) 不变因子是有限初等因子和无穷初等因子的组合编码。

!!! example "例 41A.5"
    续例 41A.3。$\det(A - \lambda B) = \lambda^2$。

    不变因子：$d_1 = 1$，$d_2 = 1$，$d_3 = \lambda^2$。

    有限初等因子：$\lambda^2$（$\lambda = 0$ 的 $2 \times 2$ Jordan 块）。

    无穷初等因子：$\mu^1$（$1 \times 1$ 幂零块）。

---

## 41A.6 收缩子空间

<div class="context-flow" markdown>

**核心问题**：什么是收缩子空间？它在矩阵束理论中扮演什么角色？

</div>

不变子空间是标准特征值问题的核心概念。对广义特征值问题，对应的概念是**收缩子空间**（deflating subspace），它为理解束的结构和设计数值算法提供了基础。

!!! definition "定义 41A.8 (收缩子空间)"
    设 $A, B \in M_n(\mathbb{C})$。$\mathbb{C}^n$ 的子空间 $\mathcal{V}$ 称为矩阵束 $A - \lambda B$ 的**右收缩子空间**（right deflating subspace），若
    $$\dim(A\mathcal{V} + B\mathcal{V}) \le \dim(\mathcal{V}).$$
    等价地，若 $V \in M_{n \times k}(\mathbb{C})$ 的列张成 $\mathcal{V}$，则存在子空间 $\mathcal{W}$（$\dim \mathcal{W} = \dim \mathcal{V}$），使得 $AV$ 和 $BV$ 的列都属于 $\mathcal{W}$。

    类似地，$\mathcal{W}$ 称为**左收缩子空间**，若 $\dim(A^*\mathcal{W} + B^*\mathcal{W}) \le \dim(\mathcal{W})$。

!!! theorem "定理 41A.6 (收缩子空间的等价刻画)"
    设 $A - \lambda B$ 为 $n \times n$ 正则束，$\mathcal{V}$ 为 $k$ 维子空间，$V \in M_{n \times k}$ 为其列基矩阵。则以下条件等价：

    (a) $\mathcal{V}$ 是右收缩子空间；

    (b) 存在矩阵 $W \in M_{n \times k}$（列满秩）和 $S, T \in M_{k \times k}$，使得 $AV = WS$，$BV = WT$；

    (c) 存在 $n \times k$ 矩阵 $W$，使得 $W^*(AV)$ 和 $W^*(BV)$ 的联合零空间为 $\{0\}$，即 $A - \lambda B$ 限制在 $\mathcal{V}$ 上后可"收缩"为 $k \times k$ 束 $S - \lambda T$。

??? proof "证明"
    $(a) \Rightarrow (b)$：由定义，$A\mathcal{V}$ 和 $B\mathcal{V}$ 都包含在某个 $k$ 维子空间 $\mathcal{W}$ 中。取 $\mathcal{W}$ 的列基矩阵 $W$，则 $AV = WS$，$BV = WT$ 对某矩阵 $S, T$ 成立。

    $(b) \Rightarrow (a)$：$A\mathcal{V} + B\mathcal{V} \subseteq \operatorname{range}(W)$，$\dim(\operatorname{range}(W)) \le k = \dim(\mathcal{V})$。

    $(b) \Leftrightarrow (c)$：取 $W$ 的列为 $\mathcal{W}$ 的正交基，则 $S = W^*AV$，$T = W^*BV$，条件即要求 $S - \lambda T$ 是 $k \times k$ 束。

!!! definition "定义 41A.9 (收缩对)"
    若 $(\mathcal{V}, \mathcal{W})$ 满足 $A\mathcal{V} \subseteq \mathcal{W}$ 且 $B\mathcal{V} \subseteq \mathcal{W}$，$\dim \mathcal{V} = \dim \mathcal{W} = k$，则称 $(\mathcal{V}, \mathcal{W})$ 为矩阵束 $A - \lambda B$ 的**收缩对**（deflating pair）。

!!! theorem "定理 41A.7 (收缩子空间与广义特征值)"
    设 $A - \lambda B$ 为正则束，$(\mathcal{V}, \mathcal{W})$ 为 $k$ 维收缩对。取列基 $V, W$，则"收缩束" $S - \lambda T = W^{-1}AV - \lambda W^{-1}BV$ 的广义特征值是原束 $A - \lambda B$ 的广义特征值的子集。

    反之，对应于广义特征值 $\lambda_1, \ldots, \lambda_k$（含可能的 $\infty$）的广义特征向量张成的子空间是收缩子空间。

!!! example "例 41A.6"
    $A = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$。

    广义特征值为 $1, 3$。取 $\mathcal{V} = \operatorname{span}\{e_1\}$（$\lambda = 1$ 的特征向量），则 $A\mathcal{V} = \operatorname{span}\{(1,0)^T\} = \mathcal{V}$，$B\mathcal{V} = \mathcal{V}$。故 $\mathcal{V}$ 是收缩子空间，收缩束为 $(1) - \lambda(1) = 1 - \lambda$，特征值为 $\lambda = 1$。

---

## 41A.7 广义 Schur 分解（QZ 分解）

<div class="context-flow" markdown>

**核心问题**：如何用酉变换同时将两个矩阵三角化？

</div>

广义 Schur 分解是标准 Schur 分解 $A = QTQ^*$ 在矩阵束中的自然推广，是数值计算广义特征值的理论基础。

!!! theorem "定理 41A.8 (广义 Schur 分解 / QZ 分解)"
    设 $A, B \in M_n(\mathbb{C})$。则存在酉矩阵 $Q, Z \in M_n(\mathbb{C})$ 使得
    $$Q^* A Z = T_A, \quad Q^* B Z = T_B,$$
    其中 $T_A, T_B$ 均为上三角矩阵。广义特征值为
    $$\lambda_i = \frac{(T_A)_{ii}}{(T_B)_{ii}} \quad \text{（当 $(T_B)_{ii} \ne 0$）}, \quad \lambda_i = \infty \quad \text{（当 $(T_B)_{ii} = 0$，$(T_A)_{ii} \ne 0$）}.$$

    更精确地，广义特征对为 $(\alpha_i, \beta_i) = ((T_A)_{ii}, (T_B)_{ii})$。

??? proof "证明"
    证明通过对矩阵阶数 $n$ 的归纳法。

    **基础情况** ($n = 1$)：$T_A = A$，$T_B = B$，$Q = Z = (1)$，结论平凡。

    **归纳步骤**：假设结论对 $n-1$ 阶矩阵成立。设 $A, B \in M_n(\mathbb{C})$。

    *情况 1*：$A - \lambda B$ 为正则束。取一个广义特征值 $\lambda_1$（有限或无穷）。

    若 $\lambda_1$ 有限，则 $\det(A - \lambda_1 B) = 0$，存在单位向量 $z_1$ 使得 $(A - \lambda_1 B)z_1 = 0$，即 $Az_1 = \lambda_1 Bz_1$。将 $z_1$ 扩充为酉矩阵 $\hat{Z} = (z_1, Z_2)$。则
    $$A\hat{Z} = (Az_1, AZ_2) = (\lambda_1 Bz_1, AZ_2).$$

    令 $w_1 = Bz_1 / \|Bz_1\|$（若 $Bz_1 \ne 0$），扩充为酉矩阵 $\hat{Q} = (w_1, Q_2)$。则
    $$\hat{Q}^* A \hat{Z} = \begin{pmatrix} w_1^* A z_1 & w_1^* A Z_2 \\ Q_2^* A z_1 & Q_2^* A Z_2 \end{pmatrix}.$$

    由 $Az_1 = \lambda_1 Bz_1$，$w_1 = Bz_1/\|Bz_1\|$，得 $Q_2^* A z_1 = \lambda_1 Q_2^* Bz_1 = \lambda_1 \|Bz_1\| Q_2^* w_1 = 0$（因为 $Q_2^* w_1 = 0$）。类似地 $Q_2^* B z_1 = 0$。

    于是
    $$\hat{Q}^* A \hat{Z} = \begin{pmatrix} \alpha_1 & * \\ 0 & A' \end{pmatrix}, \quad \hat{Q}^* B \hat{Z} = \begin{pmatrix} \beta_1 & * \\ 0 & B' \end{pmatrix},$$
    其中 $A', B' \in M_{n-1}(\mathbb{C})$。

    若 $\lambda_1 = \infty$（即 $Bz_1 = 0$），类似处理：取 $w_1 = Az_1/\|Az_1\|$。

    对 $(n-1)$ 阶矩阵对 $(A', B')$ 应用归纳假设，得酉矩阵 $Q', Z'$ 使得 $Q'^* A' Z'$ 和 $Q'^* B' Z'$ 为上三角。

    取 $Q = \hat{Q} \begin{pmatrix} 1 & 0 \\ 0 & Q' \end{pmatrix}$，$Z = \hat{Z} \begin{pmatrix} 1 & 0 \\ 0 & Z' \end{pmatrix}$，则 $Q^* A Z$ 和 $Q^* B Z$ 均为上三角。

    *情况 2*：$A - \lambda B$ 为奇异束。此时 $\det(\beta A - \alpha B) \equiv 0$ 对所有 $(\alpha, \beta)$，特别地每个 $\lambda$ 都是"广义特征值"。取任意 $z_1$ 使得 $(A - \lambda_0 B)z_1 = 0$ 对某 $\lambda_0$，按情况 1 同样的过程进行归纳。

!!! theorem "定理 41A.9 (有序广义 Schur 分解)"
    在定理 41A.8 的基础上，可以对对角元 $(\alpha_i, \beta_i)$ 按任意给定的顺序重新排列。即对任何广义特征值的排列 $\sigma$，存在酉矩阵 $Q, Z$ 使得 $(T_A)_{ii} / (T_B)_{ii}$ 按 $\sigma$ 指定的顺序排列。

    这通过一系列 $2 \times 2$ 的酉交换变换（Givens 旋转类型）实现。

---

## 41A.8 QZ 算法

<div class="context-flow" markdown>

**核心问题**：如何在数值计算中实现广义 Schur 分解？

</div>

QZ 算法是 Moler 和 Stewart 于 1973 年提出的数值算法，是 QR 算法在广义特征值问题上的自然推广。它是 LAPACK 和 MATLAB 中求解广义特征值问题的标准方法。

!!! definition "定义 41A.10 (QZ 算法)"
    **QZ 算法**分为两个阶段：

    **阶段一（预处理——化为上 Hessenberg-三角形式）**：

    1. 对 $B$ 做 QR 分解：$Q_0 B = R_B$（$R_B$ 上三角）。令 $A \leftarrow Q_0 A$，$B \leftarrow R_B$。
    2. 对 $A$ 做 Hessenberg 化：利用 Householder 变换从左边将 $A$ 化为上 Hessenberg 矩阵 $H_A$，同时用 Givens 旋转从右边恢复 $B$ 的上三角性。

    结果：$Q_1^* A Z_1 = H_A$（上 Hessenberg），$Q_1^* B Z_1 = R_B$（上三角）。

    **阶段二（隐式位移 QZ 迭代）**：

    在 Hessenberg-三角形式 $(H_A, R_B)$ 上进行迭代。每步选取位移 $\mu$（通常取 $(H_A, R_B)$ 右下角 $2 \times 2$ 子束的广义特征值），计算：

    1. 形成 $M = H_A - \mu R_B$ 的第一列 $m_1$；
    2. 构造 Householder 变换 $P_1$ 使得 $P_1 m_1 = \alpha e_1$；
    3. 从左乘 $P_1$，破坏 $H_A$ 的 Hessenberg 结构（产生一个"凸起"(bulge)）；
    4. 通过交替的左/右 Givens 旋转将凸起逐步"追赶"至矩阵底部，恢复 Hessenberg-三角形式。

    这就是**隐式双重位移**策略，类似于 QR 算法中的 Francis 位移。

!!! theorem "定理 41A.10 (QZ 算法的性质)"
    QZ 算法具有以下性质：

    (a) **复杂度**：预处理 $O(n^3)$，每步迭代 $O(n^2)$，通常 $O(n)$ 步迭代收敛，总复杂度 $O(n^3)$。

    (b) **数值后向稳定性**：计算所得的广义特征值是扰动束 $(A + \delta A) - \lambda(B + \delta B)$ 的精确广义特征值，其中 $\|\delta A\| = O(\epsilon_{\mathrm{mach}} \|A\|)$，$\|\delta B\| = O(\epsilon_{\mathrm{mach}} \|B\|)$。

    (c) **收敛性**：对几乎所有初始矩阵，使用 Wilkinson 型位移的 QZ 迭代以三次（cubic）速率收敛。子对角元 $(H_A)_{k+1,k}$ 满足
    $$|(H_A)_{k+1,k}^{(\text{new})}| = O(|(H_A)_{k+1,k}^{(\text{old})}|^3).$$

    (d) **收缩**（deflation）：当 $(H_A)_{k+1,k}$ 足够小时，将其置零，矩阵束分裂为两个更小的子问题。

!!! example "例 41A.7"
    计算束 $A - \lambda B$ 的广义特征值，$A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，$B = \begin{pmatrix} 5 & 6 \\ 7 & 8 \end{pmatrix}$。

    $\det(A - \lambda B) = (1-5\lambda)(4-8\lambda) - (2-6\lambda)(3-7\lambda)$
    $= 4 - 8\lambda - 20\lambda + 40\lambda^2 - (6 - 14\lambda - 18\lambda + 42\lambda^2)$
    $= -2\lambda^2 + 4\lambda - 2 = -2(\lambda - 1)^2$.

    广义特征值为 $\lambda = 1$（代数重数 $2$）。

    在实际 QZ 算法中，此问题经预处理化为 Hessenberg-三角形式后，仅需少数几步迭代即收敛到三角形式，对角元之比给出 $\lambda = 1$。

!!! note "注记 41A.1 (QZ 算法的实现)"
    QZ 算法在 LAPACK 中由以下子程序实现：

    - `xGGES`：计算广义 Schur 分解（三角形式），返回 $(\alpha_i, \beta_i)$；
    - `xGGEV`：计算广义特征值和特征向量；
    - `xGGHRD`：Hessenberg-三角预处理；
    - `xHGEQZ`：QZ 迭代核心。

    在 MATLAB 中，`eig(A, B)` 和 `qz(A, B)` 调用这些子程序。Python 的 `scipy.linalg.qz` 提供类似功能。

---

## 41A.9 对称正定束与同时对角化

<div class="context-flow" markdown>

**核心问题**：当矩阵束具有对称结构时，广义特征值有何特殊性质？

</div>

!!! theorem "定理 41A.11 (对称正定束的性质)"
    设 $A - \lambda B$ 为对称束（$A = A^T$，$B = B^T$），且 $B$ 正定。则：

    (a) 所有广义特征值 $\lambda_i$ 为**实数**；

    (b) 广义特征向量可选为 $B$-正交的：$v_i^T B v_j = \delta_{ij}$；

    (c) 存在非奇异矩阵 $S$ 使得 $S^T A S = \Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，$S^T B S = I$（**同时对角化**）；

    (d) 这等价于 Weierstrass 标准形退化为对角形：$J = \Lambda$，$N$ 为空（无无穷特征值，因为 $B$ 非奇异）。

??? proof "证明"
    由 $B > 0$，令 $B = LL^T$（Cholesky 分解），其中 $L$ 为下三角、对角元为正的矩阵。

    将广义问题 $Ax = \lambda Bx$ 化为标准对称特征值问题：
    $$Ax = \lambda LL^T x \implies L^{-1}A L^{-T}(L^T x) = \lambda (L^T x).$$
    令 $y = L^T x$，$\hat{A} = L^{-1}AL^{-T}$。由于 $A$ 对称，$\hat{A}$ 也对称：$\hat{A}^T = (L^{-1}AL^{-T})^T = L^{-1}A^TL^{-T} = L^{-1}AL^{-T} = \hat{A}$。

    由谱定理，$\hat{A}$ 的特征值全为实数，且存在正交矩阵 $U$ 使得 $U^T \hat{A} U = \Lambda$。

    取 $S = L^{-T}U$。则 $S^T B S = U^T L^{-1} \cdot LL^T \cdot L^{-T}U = U^T U = I$，$S^T A S = U^T L^{-1} A L^{-T} U = U^T \hat{A} U = \Lambda$。

    $B$-正交性：$S^T B S = I$ 意味着 $s_i^T B s_j = \delta_{ij}$，其中 $s_i$ 为 $S$ 的第 $i$ 列，即广义特征向量。

!!! theorem "定理 41A.12 (对称正定束的极小极大原理)"
    设 $A = A^T$，$B > 0$。设广义特征值 $\lambda_1 \le \lambda_2 \le \cdots \le \lambda_n$。则
    $$\lambda_k = \min_{\dim \mathcal{V} = k} \max_{0 \ne x \in \mathcal{V}} \frac{x^T A x}{x^T B x} = \max_{\dim \mathcal{V} = n-k+1} \min_{0 \ne x \in \mathcal{V}} \frac{x^T A x}{x^T B x}.$$
    这是 Courant-Fischer 极小极大原理在广义特征值问题上的推广。Rayleigh 商 $\rho(x) = x^T A x / x^T B x$ 取代了标准 Rayleigh 商。

!!! example "例 41A.8"
    广义特征值问题 $Ax = \lambda Bx$，
    $$A = \begin{pmatrix} 5 & 2 \\ 2 & 2 \end{pmatrix}, \quad B = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}.$$

    $\det(A - \lambda B) = (5-2\lambda)(2-\lambda) - (2-\lambda)^2 = (2-\lambda)[(5-2\lambda)-(2-\lambda)] = (2-\lambda)(3-\lambda)$。

    广义特征值 $\lambda_1 = 2, \lambda_2 = 3$，均为实数。

    $\lambda_1 = 2$：$(A-2B)x = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}x = 0$，$x_1 = (0, 1)^T$。

    $\lambda_2 = 3$：$(A-3B)x = \begin{pmatrix} -1 & -1 \\ -1 & -1 \end{pmatrix}x = 0$，$x_2 = (1, -1)^T$。

    验证 $B$-正交性：$x_1^T B x_2 = (0, 1)\begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ -1 \end{pmatrix} = (0, 1)(1, 0)^T = 0$。

---

## 41A.10 Hermite 束

<div class="context-flow" markdown>

**核心问题**：Hermite 束有哪些特殊结构性质？

</div>

!!! definition "定义 41A.11 (Hermite 束)"
    矩阵束 $A - \lambda B$ 称为 **Hermite 束**，若 $A = A^*$，$B = B^*$（即 $A, B$ 均为 Hermite 矩阵）。此时对一切实数 $\lambda$，$A - \lambda B$ 为 Hermite 矩阵。

!!! theorem "定理 41A.13 (Hermite 束的惯性)"
    设 $A - \lambda B$ 为 Hermite 束，$A, B \in M_n(\mathbb{C})$。

    (a) 若 $B > 0$（正定），则所有广义特征值为实数。

    (b) 若 $B \ge 0$（半正定），则所有**有限**广义特征值为实数。

    (c) 更一般地，设 $B$ 不定。令 $\pi(\lambda)$、$\nu(\lambda)$、$\zeta(\lambda)$ 分别为 $A - \lambda B$ 的正、负、零惯性指数。则 $\pi(\lambda)$ 和 $\nu(\lambda)$ 仅在广义特征值处发生跳变。

!!! theorem "定理 41A.14 (Hermite 束的 Crawford 数)"
    设 $A - \lambda B$ 为正则 Hermite 束。**Crawford 数**定义为
    $$\gamma(A, B) = \min_{\|x\| = 1} \sqrt{(x^* A x)^2 + (x^* B x)^2}.$$

    (a) $\gamma(A, B) > 0$ 当且仅当束是**确定的**（definite），即不存在 $\|x\| = 1$ 使得 $x^* A x = x^* B x = 0$。

    (b) 确定 Hermite 束可以通过合同变换同时对角化（类似定理 41A.11）。

    (c) Crawford 数衡量了 Hermite 束离"不确定"的距离。

??? proof "证明"
    (a) 若 $\gamma(A,B) > 0$，则对任何 $\|x\| = 1$，$(x^*Ax)^2 + (x^*Bx)^2 > 0$，即 $x^*Ax$ 和 $x^*Bx$ 不同时为零。这意味着对任何实数 $\mu$，$A - \mu B$ 的零空间中没有使得 $x^*Bx = 0$ 的向量——等价于束是确定的。

    反之，若束不确定，则存在 $\|x_0\| = 1$ 使得 $x_0^*Ax_0 = 0$，$x_0^*Bx_0 = 0$，从而 $\gamma = 0$。

    (b) 确定性意味着存在实数 $\mu$ 使得 $A - \mu B$ 正定（或负定）。以 $B' = A - \mu B > 0$ 代替 $B$，利用定理 41A.11 即得同时对角化。

!!! example "例 41A.9"
    $A = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。两者均为 Hermite。

    $\det(A - \lambda B) = 1 - \lambda^2 = (1-\lambda)(1+\lambda)$。广义特征值 $\lambda = \pm 1$，均为实数。

    检查确定性：取 $x = (1, 0)^T$，$x^*Ax = 1$，$x^*Bx = 0$。取 $x = (0, 1)^T$，$x^*Ax = -1$，$x^*Bx = 0$。取 $x = (1,1)^T/\sqrt{2}$，$x^*Ax = 0$，$x^*Bx = 1$。故存在 $x$ 使 $x^*Ax = 0$，但此时 $x^*Bx = 1 \ne 0$。Crawford 数 $\gamma > 0$，束是确定的。

!!! definition "定义 41A.12 (Hermite 正定束与半正定束)"
    Hermite 束 $A - \lambda B$（$A = A^*, B = B^*$）称为**正定束**，若存在实数 $\alpha, \beta$ 使得 $\alpha A + \beta B > 0$。称为**半正定束**，若存在 $\alpha, \beta$ 使得 $\alpha A + \beta B \ge 0$ 且 $\operatorname{rank}(\alpha A + \beta B) = n$。

---

## 41A.11 正则束的扰动理论

<div class="context-flow" markdown>

**核心问题**：广义特征值对矩阵扰动有多敏感？如何量化这种敏感性？

</div>

!!! theorem "定理 41A.15 (简单广义特征值的条件数)"
    设 $A - \lambda B$ 为正则束，$\lambda_0$ 为简单广义特征值（代数重数 $1$）。设 $x, y$ 分别为对应的右、左广义特征向量：
    $$Ax = \lambda_0 Bx, \quad y^*A = \lambda_0 y^*B.$$
    则 $\lambda_0$ 关于 $(A, B)$ 扰动的**条件数**为
    $$\kappa(\lambda_0) = \frac{\|y\| \cdot \|x\|}{|y^* B x|}.$$

    具体地，设 $\tilde{A} = A + \delta A$，$\tilde{B} = B + \delta B$，$\tilde{\lambda}$ 为 $\tilde{A} - \lambda \tilde{B}$ 对应 $\lambda_0$ 的特征值，则
    $$|\tilde{\lambda} - \lambda_0| \le \kappa(\lambda_0) (\|\delta A\| + |\lambda_0| \|\delta B\|) + O(\|\delta\|^2).$$

??? proof "证明"
    设 $\tilde{A} = A + \epsilon E$，$\tilde{B} = B + \epsilon F$，$\tilde{\lambda} = \lambda_0 + \epsilon \lambda_1 + O(\epsilon^2)$，$\tilde{x} = x + \epsilon x_1 + O(\epsilon^2)$。

    代入 $\tilde{A}\tilde{x} = \tilde{\lambda}\tilde{B}\tilde{x}$，展开至 $O(\epsilon)$：
    $$(A + \epsilon E)(x + \epsilon x_1) = (\lambda_0 + \epsilon \lambda_1)(B + \epsilon F)(x + \epsilon x_1).$$
    $O(1)$ 项：$Ax = \lambda_0 Bx$（已知）。

    $O(\epsilon)$ 项：$Ax_1 + Ex = \lambda_0 Bx_1 + \lambda_1 Bx + \lambda_0 Fx$。

    即 $(A - \lambda_0 B)x_1 = \lambda_1 Bx - Ex + \lambda_0 Fx$。

    左乘 $y^*$：$y^*(A - \lambda_0 B)x_1 = 0$（因为 $y^*(A - \lambda_0 B) = 0$），所以
    $$0 = \lambda_1 y^* Bx - y^* Ex + \lambda_0 y^* Fx.$$
    $$\lambda_1 = \frac{y^*(E - \lambda_0 F)x}{y^* Bx}.$$
    $$|\lambda_1| \le \frac{\|y\| \cdot \|E - \lambda_0 F\| \cdot \|x\|}{|y^* Bx|} = \kappa(\lambda_0) \|E - \lambda_0 F\|.$$

!!! theorem "定理 41A.16 (Bauer-Fike 型定理——广义特征值)"
    设 $A - \lambda B$ 为正则束，通过 Weierstrass 标准形 $P(A - \lambda B)Q = \operatorname{diag}(J - \lambda I, I - \lambda N)$。设 $\tilde{A} = A + \delta A$，$\tilde{B} = B + \delta B$。则对 $\tilde{A} - \lambda \tilde{B}$ 的每个有限广义特征值 $\tilde{\lambda}$，要么 $\tilde{\lambda}$ 接近原束的某个有限特征值，要么它非常大（接近无穷特征值）。具体界取决于 $\|P\|$、$\|Q\|$ 和扰动大小。

!!! note "注记 41A.2 (病态广义特征值)"
    条件数 $\kappa(\lambda_0) = \|y\|\|x\|/|y^*Bx|$ 在以下情况下很大：

    1. **$y^*Bx \approx 0$**：这可能发生在 $B$ 接近奇异时（无穷特征值"入侵"有限特征值区域）。
    2. **$\lambda_0$ 是亏损特征值**（几何重数小于代数重数）：此时 $\lambda_0$ 对扰动的敏感性为 $O(\epsilon^{1/k})$，$k$ 为 Jordan 块大小。
    3. **束接近奇异**：$\det(A - \lambda B)$ 的首项系数 $\det(B)$ 接近零时，所有有限特征值都变得病态。

!!! theorem "定理 41A.17 (Stewart-Sun 距离)"
    对两个正则束 $A_1 - \lambda B_1$ 和 $A_2 - \lambda B_2$，定义**弦距离**（chordal distance）：对广义特征值 $\lambda = \alpha/\beta$，其在 Riemann 球上的表示为 $(\alpha, \beta)/\|(\alpha, \beta)\|$，两点间的弦距离为
    $$\chi(\lambda_1, \lambda_2) = \frac{|\alpha_1 \beta_2 - \alpha_2 \beta_1|}{\|(\alpha_1, \beta_1)\| \cdot \|(\alpha_2, \beta_2)\|}.$$
    弦距离在有限和无穷特征值之间提供了统一的度量。

!!! example "例 41A.10"
    考虑 $A = I$，$B = \begin{pmatrix} 1 & 0 \\ 0 & \epsilon \end{pmatrix}$（$\epsilon$ 很小）。

    广义特征值为 $\lambda_1 = 1$，$\lambda_2 = 1/\epsilon$。当 $\epsilon \to 0$，$\lambda_2 \to \infty$（从有限"逃逸"至无穷）。

    弦距离 $\chi(\lambda_2, \infty) = \chi(1/\epsilon, \infty)$。在齐次坐标下，$\lambda_2 = (1, \epsilon)$，$\infty = (1, 0)$，$\chi = |0 - \epsilon|/(\sqrt{1+\epsilon^2} \cdot 1) \approx \epsilon$。因此弦距离自然地捕捉了"接近无穷"的行为。

---

## 41A.12 多项式特征值问题与线性化

<div class="context-flow" markdown>

**核心问题**：如何将高次多项式特征值问题转化为线性束问题？

</div>

!!! definition "定义 41A.13 (多项式特征值问题)"
    **$m$ 次多项式特征值问题**（Polynomial Eigenvalue Problem, PEP）为
    $$P(\lambda)x = 0, \quad P(\lambda) = \sum_{k=0}^{m} \lambda^k A_k,$$
    其中 $A_k \in M_n(\mathbb{C})$，$A_m \ne 0$。$P(\lambda)$ 称为**矩阵多项式**。$\lambda_0$ 为**特征值**若 $\det P(\lambda_0) = 0$。

    $m = 1$：广义特征值问题；$m = 2$：**二次特征值问题**（QEP）。

!!! theorem "定理 41A.18 (伴随线性化)"
    $m$ 次多项式特征值问题 $P(\lambda)x = 0$ 可以通过**第一伴随线性化**转化为 $mn \times mn$ 的线性束问题：
    $$(\lambda B_c - A_c)\tilde{x} = 0,$$
    其中
    $$A_c = \begin{pmatrix} -A_{m-1} & -A_{m-2} & \cdots & -A_1 & -A_0 \\ I & 0 & \cdots & 0 & 0 \\ 0 & I & \cdots & 0 & 0 \\ \vdots & & \ddots & & \vdots \\ 0 & 0 & \cdots & I & 0 \end{pmatrix}, \quad B_c = \begin{pmatrix} A_m & 0 & \cdots & 0 \\ 0 & -I & \cdots & 0 \\ \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & -I \end{pmatrix},$$
    $\tilde{x} = (\lambda^{m-1}x^T, \lambda^{m-2}x^T, \ldots, x^T)^T$。

    **第二伴随线性化**类似，转置角色互换。两种线性化保持特征值（含重数和部分结构）。

??? proof "证明"
    将 $P(\lambda)x = 0$ 展开为
    $$\lambda^m A_m x + \lambda^{m-1} A_{m-1} x + \cdots + \lambda A_1 x + A_0 x = 0.$$
    引入辅助变量 $x_k = \lambda^k x$（$k = 0, 1, \ldots, m-1$），则 $x_{k+1} = \lambda x_k$（$k = 0, \ldots, m-2$），且
    $$\lambda A_m x_{m-1} + A_{m-1} x_{m-1} + A_{m-2} x_{m-2} + \cdots + A_0 x_0 = 0.$$
    这给出 $m$ 个方程（$1$ 个来自 $P(\lambda)x = 0$，$m-1$ 个来自 $x_{k+1} = \lambda x_k$），共 $mn$ 个标量方程。

    将 $\lambda$ 提取到左边，常数部分放右边：
    $$\lambda \begin{pmatrix} A_m & 0 & \cdots & 0 \\ 0 & -I & \cdots & 0 \\ \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & -I \end{pmatrix} \begin{pmatrix} x_{m-1} \\ x_{m-2} \\ \vdots \\ x_0 \end{pmatrix} = \begin{pmatrix} -A_{m-1} & \cdots & -A_0 \\ I & \cdots & 0 \\ \vdots & \ddots & \vdots \\ 0 & \cdots & 0 \end{pmatrix} \begin{pmatrix} x_{m-1} \\ x_{m-2} \\ \vdots \\ x_0 \end{pmatrix}.$$
    即 $\lambda B_c \tilde{x} = A_c \tilde{x}$。特征值完全一致。

!!! example "例 41A.11"
    二次特征值问题（结构振动分析）：
    $$(\lambda^2 M + \lambda C + K)x = 0,$$
    其中 $M$ 为质量矩阵，$C$ 为阻尼矩阵，$K$ 为刚度矩阵（均为 $n \times n$）。

    第一伴随线性化：
    $$\lambda \begin{pmatrix} M & 0 \\ 0 & -I \end{pmatrix} \begin{pmatrix} \lambda x \\ x \end{pmatrix} = \begin{pmatrix} -C & -K \\ I & 0 \end{pmatrix} \begin{pmatrix} \lambda x \\ x \end{pmatrix}.$$
    这是 $2n \times 2n$ 的广义特征值问题，可用 QZ 算法求解。原问题有 $2n$ 个特征值（含无穷特征值——当 $M$ 奇异时）。

---

## 41A.13 保结构线性化

<div class="context-flow" markdown>

**核心问题**：如何构造保持矩阵多项式结构的线性化？

</div>

当矩阵多项式具有对称、Hermite 或其他结构时，标准伴随线性化通常破坏这些结构。Lancaster、Mackey-Mackey-Mehl-Mehrmann 等人发展了**保结构线性化**理论。

!!! definition "定义 41A.14 (线性化)"
    设 $P(\lambda) = \sum_{k=0}^m \lambda^k A_k$ 为 $n \times n$ 矩阵多项式。$mn \times mn$ 矩阵束 $L(\lambda) = \lambda X + Y$ 称为 $P(\lambda)$ 的**线性化**，若存在 $\mathbb{C}[\lambda]$ 上的可逆矩阵 $E(\lambda), F(\lambda)$ 使得
    $$E(\lambda)L(\lambda)F(\lambda) = \begin{pmatrix} P(\lambda) & 0 \\ 0 & I_{(m-1)n} \end{pmatrix}.$$

!!! definition "定义 41A.15 (向量空间 $\mathbb{L}_1$, $\mathbb{L}_2$)"
    Mackey、Mackey、Mehl 和 Mehrmann (2006) 定义了两个线性化空间：

    $$\mathbb{L}_1(P) = \{L(\lambda) = \lambda X + Y : L(\lambda)(\Lambda \otimes I_n) = v \otimes P(\lambda)\},$$
    $$\mathbb{L}_2(P) = \{L(\lambda) = \lambda X + Y : (\Lambda^T \otimes I_n)L(\lambda) = w \otimes P(\lambda)\},$$
    其中 $\Lambda = (\lambda^{m-1}, \lambda^{m-2}, \ldots, 1)^T$，$v, w \in \mathbb{C}^m$ 为参数向量。

    $\mathbb{L}_1$ 中的线性化保持右特征向量结构（可以从线性化的特征向量恢复 $P(\lambda)$ 的特征向量），$\mathbb{L}_2$ 保持左特征向量结构。

!!! theorem "定理 41A.19 (保结构线性化定理)"
    设 $P(\lambda) = \sum_{k=0}^m \lambda^k A_k$ 为 $n \times n$ 矩阵多项式。

    (a) $\mathbb{L}_1(P)$ 和 $\mathbb{L}_2(P)$ 各为 $m$ 维向量空间（参数化为 $v$ 或 $w$）。

    (b) 第一和第二伴随线性化分别属于 $\mathbb{L}_1(P)$ 和 $\mathbb{L}_2(P)$。

    (c) $\mathbb{DL}(P) = \mathbb{L}_1(P) \cap \mathbb{L}_2(P)$ 构成一个 $m$ 维子空间，其中的线性化同时保持左右特征向量结构。

    (d) 若 $P(\lambda)$ 为 $T$-偶（$A_k^T = A_k$ 对偶数 $k$，$A_k^T = -A_k$ 对奇数 $k$），则 $\mathbb{DL}(P)$ 中存在保持 $T$-偶结构的线性化，即 $L(\lambda) = \lambda X + Y$ 满足 $X^T = X$，$Y^T = -Y$。

    (e) 类似结果对 $T$-奇、$*$-偶、$*$-奇、回文（palindromic）和反回文结构成立。

!!! theorem "定理 41A.20 (Lancaster 线性化)"
    设 $P(\lambda) = \lambda^2 A_2 + \lambda A_1 + A_0$（二次矩阵多项式），$A_2$ 非奇异。若 $P$ 的特征值互不相同，则存在**可解耦合**（solvents）$S_1, S_2 \in M_n(\mathbb{C})$，使得
    $$P(\lambda) = A_2(\lambda I - S_1)(\lambda I - S_2),$$
    其中 $S_1, S_2$ 的特征值构成 $P(\lambda)$ 的 $2n$ 个特征值的一个分划。

    这给出了一种"非线性化"方法——直接将二次问题分解为两个标准特征值问题。

!!! example "例 41A.12"
    考虑 $P(\lambda) = \lambda^2 I + \lambda C + K$（$M = I$，单位质量矩阵，$n = 2$）。

    $$C = \begin{pmatrix} 0.1 & 0 \\ 0 & 0.2 \end{pmatrix}, \quad K = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix}.$$

    $P$ 为对称矩阵多项式（$C^T = C$，$K^T = K$）。

    $\mathbb{DL}(P)$ 中的保结构线性化（取 $v = (1, 0)^T$）：
    $$L(\lambda) = \lambda \begin{pmatrix} I & 0 \\ 0 & K \end{pmatrix} + \begin{pmatrix} C & K \\ -K & 0 \end{pmatrix}.$$
    $X$ 对称，$Y$ 反对称：这是 $T$-奇结构的线性化。它保持了谱对称性（特征值关于实轴对称）。

    标准伴随线性化 $\lambda \begin{pmatrix} I & 0 \\ 0 & -I \end{pmatrix} + \begin{pmatrix} C & K \\ -I & 0 \end{pmatrix}$ 则不保持任何对称结构。

---

## 41A.14 习题

!!! exercise "习题 41A"

    **基础题**

    1. 判断以下矩阵束的正则性，并求广义特征值（含无穷特征值）：

        (a) $A = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$。

        (b) $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。

        (c) $A = \begin{pmatrix} 1 & 2 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}$。

    2. 设 $A - \lambda B$ 为正则束，证明：$\lambda_0$ 为有限广义特征值当且仅当 $\mu_0 = 1/\lambda_0$ 为束 $B - \mu A$ 的有限广义特征值（假设 $\lambda_0 \ne 0$）。

    3. 求以下矩阵束的 Weierstrass 标准形（指出 $J$ 和 $N$）：

        $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$。

    4. 证明严格等价是等价关系（自反性、对称性、传递性）。

    5. 设 $A, B \in M_n(\mathbb{C})$，$B$ 非奇异。证明 $A - \lambda B$ 为正则束，且 Weierstrass 标准形中 $N$ 为空矩阵（无无穷特征值）。

    **进阶题**

    6. 设 $\mathcal{V}$ 为 $A - \lambda B$ 的 $k$ 维收缩子空间。证明 $\mathcal{V}^\perp$ 不一定是收缩子空间（给出反例），但对 Hermite 束（$A = A^*$，$B = B^*$），$B$-正交补 $\mathcal{V}^{\perp_B} = \{u : u^*Bv = 0, \forall v \in \mathcal{V}\}$ 是收缩子空间（当 $B$ 正定时）。

    7. **(QZ 算法的一步迭代)** 设 $A, B \in M_2(\mathbb{C})$，$B$ 上三角，$A$ 上 Hessenberg。选取位移 $\mu$，执行 QZ 迭代的一步。对
        $$A = \begin{pmatrix} 3 & 1 \\ 2 & 4 \end{pmatrix}, \quad B = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix},$$
        用位移 $\mu = 4$ 完成一步 QZ 迭代，验证收敛加速。

    8. **(Hermite 束)** 设 $A, B \in M_n(\mathbb{C})$ 均为 Hermite 矩阵。

        (a) 证明束的有限广义特征值全为实数的充分必要条件是存在实数 $\alpha, \beta$（不全为零）使得 $\alpha A + \beta B$ 正定。

        (b) 构造一个 Hermite 束，使其有限广义特征值不全为实数。

    9. **(扰动理论)** 设 $A - \lambda B$ 正则，$\lambda_0$ 为 $2 \times 2$ Jordan 块的亏损特征值。证明在扰动 $\delta A, \delta B$（$\|\delta\| = \epsilon$）下，$\lambda_0$ 分裂为两个特征值 $\lambda_0 \pm O(\epsilon^{1/2})$。

    10. **(二次特征值问题)** 对阻尼振动系统 $(\lambda^2 M + \lambda C + K)x = 0$，$M = I$，$C = \begin{pmatrix} 0.5 & 0 \\ 0 & 1 \end{pmatrix}$，$K = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$：

        (a) 写出第一伴随线性化和 $\mathbb{DL}(P)$ 中的一个保结构线性化；

        (b) 用手算或数值方法求所有 $4$ 个特征值；

        (c) 验证特征值关于实轴对称（因为系数矩阵均为实对称）。

    11. **(极小极大原理)** 利用定理 41A.12，证明以下交错不等式：设 $A, B \in M_n(\mathbb{R})$，$A = A^T$，$B > 0$，广义特征值 $\lambda_1 \le \cdots \le \lambda_n$。设 $\hat{A}, \hat{B}$ 分别为 $A, B$ 的 $(n-1) \times (n-1)$ 主子矩阵的对应束，广义特征值 $\hat{\lambda}_1 \le \cdots \le \hat{\lambda}_{n-1}$。则
        $$\lambda_1 \le \hat{\lambda}_1 \le \lambda_2 \le \hat{\lambda}_2 \le \cdots \le \hat{\lambda}_{n-1} \le \lambda_n.$$

    12. **(Lancaster 分解)** 对 $P(\lambda) = \lambda^2 I + \lambda \begin{pmatrix} 3 & 0 \\ 0 & 3 \end{pmatrix} + \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$，

        (a) 求所有特征值；

        (b) 找到两个 solvent $S_1, S_2$ 使得 $P(\lambda) = (\lambda I - S_1)(\lambda I - S_2)$；

        (c) 验证 $S_1 + S_2 = -C$，$S_1 S_2 = K$。

    **挑战题**

    13. **(Crawford 数的计算)** 证明 Hermite 束 $A - \lambda B$ 的 Crawford 数 $\gamma(A, B)$ 等于以下优化问题的最优值：
        $$\gamma = \min_{\theta \in [0, 2\pi)} \lambda_{\min}(\cos\theta \cdot A + \sin\theta \cdot B),$$
        其中 $\lambda_{\min}$ 表示最小特征值。

    14. **(保结构线性化空间的维数)** 设 $P(\lambda) = \sum_{k=0}^m \lambda^k A_k$，$A_k \in M_n(\mathbb{C})$。证明 $\dim \mathbb{L}_1(P) = m$，$\dim \mathbb{L}_2(P) = m$，且当 $P$ 的首项系数 $A_m$ 非奇异时，$\dim(\mathbb{L}_1 \cap \mathbb{L}_2) = m$。

    15. **(广义 Schur 分解的构造性证明)** 不使用归纳法，给出广义 Schur 分解的另一种证明：利用收缩子空间的嵌套链 $\{0\} = \mathcal{V}_0 \subset \mathcal{V}_1 \subset \cdots \subset \mathcal{V}_n = \mathbb{C}^n$，其中每个 $\mathcal{V}_k$ 是 $k$ 维收缩子空间。证明这样的链总是存在的，并从中构造 $Q, Z$。
