# 第 64 章 矩阵空间中的凸性

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 矩阵不等式(Ch18) · 优化(Ch25)

**本章脉络**：矩阵空间中的凸集 → PSD 锥 → 双随机多面体 → 相关矩阵集 → 分离定理(迹内积) → 矩阵凸函数 → 矩阵凹函数 → 极小极大定理 → 应用到 SDP

**延伸**：矩阵凸性理论是半定规划和量子信息论（量子态空间的凸结构、纠缠判定的分离定理）的数学基础；算子凸函数理论连接了矩阵分析与泛函分析

</div>

凸性是优化理论的核心概念。当我们将凸分析从 $\mathbb{R}^n$ 推广到矩阵空间 $M_n(\mathbb{R})$ 或 Hermite 矩阵空间 $S_n(\mathbb{R})$ 时，出现了许多具有丰富结构的凸集和凸函数。半正定锥 $S_n^+$ 是最重要的矩阵凸锥，它的几何性质（自对偶、正则锥）使得半定规划成为可能。Birkhoff 多面体将置换矩阵与双随机矩阵联系起来。相关矩阵集（椭圆体）在金融数学和统计学中无处不在。

本章系统研究矩阵空间中的凸性理论，从矩阵空间的线性结构出发，经由重要的凸集实例，到达矩阵凸函数理论和极小极大定理。

---

## 64.1 矩阵空间的线性结构

<div class="context-flow" markdown>

**核心问题**：如何将 $M_n$ 视为一个赋以内积结构的有限维实向量空间？

</div>

### 矩阵空间的维数

!!! definition "定义 64.1 (矩阵空间)"
    以下矩阵空间是有限维实向量空间：

    - $M_n(\mathbb{R})$：全体 $n \times n$ 实矩阵，实维数 $n^2$；
    - $M_{m \times n}(\mathbb{R})$：全体 $m \times n$ 实矩阵，实维数 $mn$；
    - $S_n(\mathbb{R}) = \{A \in M_n(\mathbb{R}) : A = A^T\}$：实对称矩阵，实维数 $n(n+1)/2$；
    - $\mathcal{A}_n(\mathbb{R}) = \{A \in M_n(\mathbb{R}) : A = -A^T\}$：反对称矩阵，实维数 $n(n-1)/2$；
    - $M_n(\mathbb{C})$：全体 $n \times n$ 复矩阵，实维数 $2n^2$；
    - $H_n(\mathbb{C}) = \{A \in M_n(\mathbb{C}) : A = A^*\}$：Hermite 矩阵，实维数 $n^2$。

### 迹内积

!!! definition "定义 64.2 (迹内积)"
    在实矩阵空间 $M_n(\mathbb{R})$ 上，**迹内积**（trace inner product）定义为
    $$
    \langle A, B \rangle = \operatorname{tr}(A^T B) = \sum_{i,j} a_{ij} b_{ij}
    $$

    在 Hermite 矩阵空间 $H_n(\mathbb{C})$ 上，迹内积定义为
    $$
    \langle A, B \rangle = \operatorname{tr}(AB) = \sum_{i,j} a_{ij} \overline{b_{ji}}
    $$
    注意对 Hermite 矩阵，$\operatorname{tr}(AB) \in \mathbb{R}$。

!!! theorem "定理 64.1 (迹内积的性质)"
    迹内积是一个内积，即满足：

    1. **双线性**（或共轭双线性）：$\langle \alpha A + \beta B, C \rangle = \alpha \langle A, C \rangle + \beta \langle B, C \rangle$；
    2. **对称性**：$\langle A, B \rangle = \langle B, A \rangle$（实情形）或 $\langle A, B \rangle = \overline{\langle B, A \rangle}$（复情形）；
    3. **正定性**：$\langle A, A \rangle = \operatorname{tr}(A^T A) = \sum_{i,j} a_{ij}^2 \geq 0$，等号成立当且仅当 $A = 0$。

    对应的范数 $\|A\|_F = \sqrt{\langle A, A \rangle} = \sqrt{\operatorname{tr}(A^T A)}$ 就是 Frobenius 范数。

!!! note "注"
    $(M_n(\mathbb{R}), \langle \cdot, \cdot \rangle)$ 同构于 $(\mathbb{R}^{n^2}, \langle \cdot, \cdot \rangle_{\text{std}})$，因此 $\mathbb{R}^{n^2}$ 中的所有凸分析结论都可以直接搬到矩阵空间上来。

### 矩阵空间的直和分解

!!! theorem "定理 64.2 (对称-反对称分解)"
    $M_n(\mathbb{R})$ 在迹内积下是对称矩阵空间和反对称矩阵空间的正交直和：
    $$
    M_n(\mathbb{R}) = S_n(\mathbb{R}) \oplus \mathcal{A}_n(\mathbb{R})
    $$
    投影算子为 $\Pi_S(A) = \frac{A + A^T}{2}$，$\Pi_A(A) = \frac{A - A^T}{2}$。

---

## 64.2 PSD 锥

<div class="context-flow" markdown>

**核心问题**：半正定矩阵集合具有什么几何结构？

</div>

### PSD 锥的定义

!!! definition "定义 64.3 (半正定锥)"
    **半正定锥**（positive semidefinite cone）定义为
    $$
    S_n^+ = \{A \in S_n(\mathbb{R}) : A \succeq 0\} = \{A \in S_n(\mathbb{R}) : x^T A x \geq 0, \; \forall x \in \mathbb{R}^n\}
    $$

    **正定锥**为其内部：
    $$
    S_n^{++} = \{A \in S_n(\mathbb{R}) : A \succ 0\} = \operatorname{int}(S_n^+)
    $$

### 正则锥

!!! definition "定义 64.4 (正则锥)"
    $\mathbb{R}^d$ 中的锥 $K$ 称为**正则锥**（proper cone），若它同时满足：

    1. **凸性**：$A, B \in K$，$\alpha, \beta \geq 0$ 蕴含 $\alpha A + \beta B \in K$；
    2. **闭性**：$K$ 是闭集；
    3. **尖性**（pointed）：$K \cap (-K) = \{0\}$；
    4. **实性**（solid）：$K$ 有非空内部，即 $\operatorname{int}(K) \neq \emptyset$。

!!! theorem "定理 64.3 (PSD 锥是正则锥)"
    $S_n^+$ 是 $S_n(\mathbb{R})$ 中的正则锥。

??? proof "证明"
    1. **凸性**：若 $A, B \succeq 0$，$\alpha, \beta \geq 0$，则对任意 $x$，$x^T(\alpha A + \beta B)x = \alpha x^TAx + \beta x^TBx \geq 0$。

    2. **闭性**：$S_n^+$ 可以写成 $\bigcap_{x \in \mathbb{R}^n} \{A \in S_n : x^TAx \geq 0\}$，是闭半空间的交，因此闭。

    3. **尖性**：若 $A \succeq 0$ 且 $-A \succeq 0$，则 $x^TAx = 0$ 对所有 $x$，故 $A = 0$。

    4. **实性**：$I \in S_n^{++}$，且 $I$ 的任意足够小的邻域中的对称矩阵仍是正定的，因此 $S_n^{++} = \operatorname{int}(S_n^+) \neq \emptyset$。

### 自对偶性

!!! definition "定义 64.5 (对偶锥)"
    锥 $K \subseteq V$ 的**对偶锥**（dual cone）定义为
    $$
    K^* = \{Y \in V : \langle X, Y \rangle \geq 0, \; \forall X \in K\}
    $$

!!! theorem "定理 64.4 (PSD 锥的自对偶性)"
    在迹内积下，$S_n^+$ 是自对偶的：$(S_n^+)^* = S_n^+$。

??? proof "证明"
    **$S_n^+ \subseteq (S_n^+)^*$**：若 $A, B \in S_n^+$，则 $\langle A, B \rangle = \operatorname{tr}(AB) \geq 0$（因为 $AB$ 的特征值非负，注意 $A, B$ 半正定时 $AB$ 与 $A^{1/2}BA^{1/2} \succeq 0$ 相似，故迹非负）。

    **$(S_n^+)^* \subseteq S_n^+$**：设 $Y \in (S_n^+)^*$，即 $\operatorname{tr}(XY) \geq 0$ 对所有 $X \succeq 0$。取 $X = vv^T$（秩 1 半正定矩阵），得 $v^T Y v = \operatorname{tr}(vv^T Y) \geq 0$ 对所有 $v$。因此 $Y \succeq 0$。

### 极端射线

!!! theorem "定理 64.5 (PSD 锥的极端射线)"
    $S_n^+$ 的极端射线（extreme ray）恰好是秩 1 半正定矩阵 $\{xx^T : x \in \mathbb{R}^n \setminus \{0\}\}$ 生成的射线。

    等价地，$S_n^+$ 的极端点（在截面 $\{A : \operatorname{tr}(A) = 1\}$ 上）是秩 1 的密度矩阵 $xx^T / \|x\|^2$。

??? proof "证明"
    设 $A \in S_n^+$ 且 $\operatorname{rank}(A) = 1$，即 $A = xx^T$。若 $A = B + C$，$B, C \in S_n^+$，则 $\operatorname{rank}(B), \operatorname{rank}(C) \leq \operatorname{rank}(A) = 1$（因为 $B, C$ 的列空间都含在 $A$ 的列空间中）。因此 $B = \alpha xx^T$，$C = \beta xx^T$，$\alpha + \beta = 1$，$\alpha, \beta \geq 0$。这说明 $A$ 生成的射线是极端的。

    反之，若 $\operatorname{rank}(A) \geq 2$，设 $A = \sum_{i=1}^r \lambda_i u_i u_i^T$（$r \geq 2$）。取 $B = \lambda_1 u_1 u_1^T$，$C = A - B$，则 $A = B + C$ 是非平凡分解。因此 $A$ 不在极端射线上。

!!! example "例 64.1"
    $S_2^+$ 中的矩阵 $A = \begin{pmatrix} a & b \\ b & c \end{pmatrix} \succeq 0$ 需满足 $a \geq 0$，$c \geq 0$，$ac - b^2 \geq 0$。

    这是 $(a, b, c)$ 空间中的一个圆锥面（二次锥），其极端射线由
    $$
    \begin{pmatrix} \cos^2\theta & \cos\theta\sin\theta \\ \cos\theta\sin\theta & \sin^2\theta \end{pmatrix} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}\begin{pmatrix} \cos\theta & \sin\theta \end{pmatrix}
    $$
    参数化。

---

## 64.3 双随机多面体

<div class="context-flow" markdown>

**核心问题**：双随机矩阵集合的凸几何结构是什么？

</div>

### Birkhoff 多面体

!!! definition "定义 64.6 (双随机矩阵)"
    $n \times n$ 矩阵 $A$ 称为**双随机**（doubly stochastic），若 $A$ 的所有元素非负，且每行、每列之和均为 $1$。全体 $n \times n$ 双随机矩阵的集合称为 **Birkhoff 多面体**，记作 $\mathcal{B}_n$。

    形式化地：
    $$
    \mathcal{B}_n = \left\{A \in M_n(\mathbb{R}) : a_{ij} \geq 0, \; \sum_j a_{ij} = 1 \; \forall i, \; \sum_i a_{ij} = 1 \; \forall j \right\}
    $$

!!! theorem "定理 64.6 (Birkhoff 定理)"
    $\mathcal{B}_n$ 是凸多面体，其顶点恰好是 $n \times n$ 置换矩阵。即每个双随机矩阵都可以写成置换矩阵的凸组合。

    因此
    $$
    \mathcal{B}_n = \operatorname{conv}\{P_\sigma : \sigma \in S_n\}
    $$
    其中 $P_\sigma$ 是置换 $\sigma$ 对应的置换矩阵，$S_n$ 是 $n$ 元对称群。

??? proof "证明"
    **顶点是置换矩阵**：设 $A$ 是 $\mathcal{B}_n$ 的顶点。$A$ 不能写成 $\mathcal{B}_n$ 中两个不同矩阵的严格凸组合。

    若 $A$ 不是置换矩阵，则存在某行有至少两个正元素，比如 $a_{i,j_1} > 0$ 和 $a_{i,j_2} > 0$。利用双随机性条件，可以找到一个"交替圈"（alternating cycle），沿此圈交替加减一个小量 $\varepsilon > 0$ 得到两个不同的双随机矩阵 $A \pm \varepsilon B$，从而 $A = \frac{1}{2}(A+\varepsilon B) + \frac{1}{2}(A - \varepsilon B)$，矛盾。

    **每个双随机矩阵是置换矩阵的凸组合**：由 Konig 定理（或等价地，Hall 定理），每个双随机矩阵的正元素支撑集包含一个完美匹配。设对应的置换矩阵为 $P$，令 $\alpha = \min_{P_{ij}=1} a_{ij} > 0$。则 $A - \alpha P$ 的行列和均为 $1 - \alpha$，因此 $\frac{1}{1-\alpha}(A - \alpha P)$ 仍是双随机矩阵（若 $\alpha < 1$）。

    通过归纳，$A$ 可以分解为至多 $n^2 - 2n + 2$ 个置换矩阵的凸组合。

### Birkhoff 多面体的维数

!!! theorem "定理 64.7 (Birkhoff 多面体的维数)"
    $\mathcal{B}_n$ 是 $M_n(\mathbb{R})$ 中维数为 $(n-1)^2$ 的凸多面体。

??? proof "证明"
    双随机矩阵的行和条件给出 $n$ 个等式约束，列和条件给出 $n$ 个等式约束。但这 $2n$ 个约束中有一个是冗余的（所有行和之和 = 所有列和之和 = 所有元素之和）。因此独立等式约束有 $2n - 1$ 个。

    $M_n(\mathbb{R})$ 的维数为 $n^2$。$\mathcal{B}_n$ 所在的仿射子空间维数为 $n^2 - (2n - 1) = (n-1)^2$。可以验证 $\mathcal{B}_n$ 在此仿射子空间中是满维的，因此 $\dim(\mathcal{B}_n) = (n-1)^2$。

!!! example "例 64.2"
    $n = 2$ 时，$\mathcal{B}_2$ 是 $(2-1)^2 = 1$ 维的，即一条线段。所有 $2 \times 2$ 双随机矩阵形如
    $$
    A = \begin{pmatrix} t & 1-t \\ 1-t & t \end{pmatrix}, \quad t \in [0, 1]
    $$
    两个顶点是 $I_2$（$t=1$）和 $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$（$t=0$）。

!!! example "例 64.3"
    $n = 3$ 时，$\mathcal{B}_3$ 是 $4$ 维多面体，有 $3! = 6$ 个顶点（置换矩阵）。

### 面格

!!! note "注"
    $\mathcal{B}_n$ 的面格（face lattice）具有丰富的组合结构，与有向图的完美匹配密切相关。一个面 $F$ 对应于一个 $n \times n$ 的 $\{0, *\}$ 模式，其中 $*$ 位置允许取正值。这与二部图的匹配理论深度交织。

---

## 64.4 相关矩阵集

<div class="context-flow" markdown>

**核心问题**：对角线元素固定为 1 的正半定矩阵集合有什么几何结构？

</div>

### 相关矩阵

!!! definition "定义 64.7 (相关矩阵)"
    $n \times n$ 实矩阵 $C$ 称为**相关矩阵**（correlation matrix），若

    1. $C$ 是半正定的：$C \succeq 0$；
    2. 对角线元素全为 $1$：$c_{ii} = 1$，$i = 1, \ldots, n$。

    全体 $n \times n$ 相关矩阵的集合记作 $\mathcal{E}_n$，也称为**椭圆体**（elliptope）。

!!! theorem "定理 64.8 (椭圆体的凸性)"
    $\mathcal{E}_n$ 是 $S_n(\mathbb{R})$ 中的凸紧集。

??? proof "证明"
    **凸性**：若 $C_1, C_2 \in \mathcal{E}_n$，$t \in [0,1]$，则 $C = tC_1 + (1-t)C_2$ 满足 $C \succeq 0$（PSD 锥的凸性）且 $c_{ii} = t \cdot 1 + (1-t) \cdot 1 = 1$。

    **有界性**：由 $c_{ii} = 1$ 和 $C \succeq 0$，柯西不等式给出 $|c_{ij}| \leq \sqrt{c_{ii} c_{jj}} = 1$。因此 $\mathcal{E}_n$ 有界。

    **闭性**：$\mathcal{E}_n = S_n^+ \cap \{C : c_{ii} = 1, \forall i\}$ 是闭集的交。

!!! example "例 64.4"
    $n = 2$ 时，$\mathcal{E}_2 = \left\{\begin{pmatrix} 1 & r \\ r & 1 \end{pmatrix} : r \in [-1, 1]\right\}$。

    这是一条线段，参数为相关系数 $r$。极端点为 $r = 1$（完全正相关）和 $r = -1$（完全负相关）。

!!! example "例 64.5"
    $n = 3$ 时，$\mathcal{E}_3$ 由满足以下条件的 $(r_{12}, r_{13}, r_{23})$ 参数化：
    $$
    \begin{pmatrix} 1 & r_{12} & r_{13} \\ r_{12} & 1 & r_{23} \\ r_{13} & r_{23} & 1 \end{pmatrix} \succeq 0
    $$
    即 $1 - r_{12}^2 - r_{13}^2 - r_{23}^2 + 2r_{12}r_{13}r_{23} \geq 0$，$|r_{ij}| \leq 1$。

    这是 $[-1,1]^3$ 中一个由行列式条件刻画的凸体，形状类似于被"膨胀"的正八面体。

### 极端点

!!! theorem "定理 64.9 (椭圆体的极端点)"
    $\mathcal{E}_n$ 的极端点恰好是秩不超过 $r$ 的相关矩阵，其中 $r(r+1)/2 \leq n$。

    特别地，单位矩阵 $I_n$ 始终是极端点。当 $n \leq 3$ 时，所有秩 1 相关矩阵（形如 $vv^T$，$v_i = \pm 1$）也是极端点。

!!! note "注"
    在金融数学中，相关矩阵描述资产收益率之间的相关关系。相关矩阵集的凸性对投资组合优化至关重要：在给定相关结构约束下的优化问题是凸优化问题。

---

## 64.5 分离定理

<div class="context-flow" markdown>

**核心问题**：如何用超平面分离矩阵空间中的凸集？

</div>

### 矩阵空间中的 Hahn-Banach 定理

!!! theorem "定理 64.10 (矩阵空间的分离定理)"
    设 $C \subseteq S_n(\mathbb{R})$ 是闭凸集，$A \notin C$。则存在 $H \in S_n(\mathbb{R})$ 和 $\alpha \in \mathbb{R}$，使得
    $$
    \operatorname{tr}(AH) > \alpha \geq \operatorname{tr}(XH), \quad \forall X \in C
    $$

    即超平面 $\{X : \operatorname{tr}(XH) = \alpha\}$ 将 $A$ 和 $C$ 严格分离。

??? proof "证明"
    由于 $(S_n(\mathbb{R}), \langle \cdot, \cdot \rangle)$ 是有限维内积空间（$\langle X, Y \rangle = \operatorname{tr}(XY)$），这是标准 Hahn-Banach 分离定理的直接推论。

    设 $X^*$ 是 $C$ 中到 $A$ 的最近点（闭凸集的最近点存在且唯一）。令 $H = A - X^*$。则对任意 $X \in C$：
    $$
    \operatorname{tr}((X - X^*)H) = \operatorname{tr}((X - X^*)(A - X^*)) \leq 0
    $$
    （最近点的变分不等式刻画），因此
    $$
    \operatorname{tr}(XH) \leq \operatorname{tr}(X^* H) < \operatorname{tr}(AH)
    $$
    取 $\alpha = \operatorname{tr}(X^* H)$ 即可。

### 对 PSD 锥的应用

!!! theorem "定理 64.11 (PSD 锥的分离)"
    若 $A \in S_n(\mathbb{R})$ 但 $A \notin S_n^+$（即 $A$ 不是半正定的），则存在 $H \succeq 0$，$H \neq 0$，使得
    $$
    \operatorname{tr}(AH) < 0
    $$

    等价地，$A$ 有负特征值当且仅当存在半正定矩阵 $H$ 使得 $\operatorname{tr}(AH) < 0$。

??? proof "证明"
    由 $S_n^+$ 的自对偶性（定理 64.4），若 $A \notin S_n^+$，则存在 $H \in S_n^+$ 使得 $\langle A, H \rangle = \operatorname{tr}(AH) < 0$。

    具体地，设 $\lambda$ 是 $A$ 的负特征值，$v$ 是对应的单位特征向量。取 $H = vv^T \succeq 0$，则
    $$
    \operatorname{tr}(AH) = \operatorname{tr}(A \cdot vv^T) = v^T A v = \lambda < 0
    $$

!!! example "例 64.6"
    设 $A = \begin{pmatrix} 1 & 3 \\ 3 & 1 \end{pmatrix}$。$A$ 的特征值为 $4$ 和 $-2$，因此 $A \notin S_2^+$。

    取 $v = \frac{1}{\sqrt{2}}(1, -1)^T$（对应特征值 $-2$），$H = vv^T = \frac{1}{2}\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$，则
    $$
    \operatorname{tr}(AH) = \frac{1}{2}\operatorname{tr}\begin{pmatrix} -2 & 2 \\ 2 & -2 \end{pmatrix} = -2 < 0
    $$

### 在量子信息中的应用

!!! note "注"
    在量子信息论中，量子态是密度矩阵 $\rho \in S_n^+$，$\operatorname{tr}(\rho) = 1$。判断一个二体量子态 $\rho_{AB}$ 是否可分离（separable）等价于判断 $\rho_{AB}$ 是否属于某个凸集。分离定理保证，每个纠缠态都可以被一个**纠缠见证**（entanglement witness）$W$ 检测出来：$\operatorname{tr}(W\rho_{AB}) < 0$。

---

## 64.6 矩阵凸函数

<div class="context-flow" markdown>

**核心问题**：哪些从对称矩阵到实数的函数是凸的？

</div>

### 定义

!!! definition "定义 64.8 (矩阵凸函数)"
    函数 $f: \mathcal{D} \to \mathbb{R}$（其中 $\mathcal{D} \subseteq S_n(\mathbb{R})$ 是凸集）称为**凸的**，若
    $$
    f(tA + (1-t)B) \leq t f(A) + (1-t) f(B), \quad \forall A, B \in \mathcal{D}, \; t \in [0,1]
    $$

### 重要例子

!!! theorem "定理 64.12 (矩阵凸函数的例子)"
    以下函数在指定的域上是凸的：

    1. **最大特征值**：$f(A) = \lambda_{\max}(A)$ 在 $S_n$ 上凸；
    2. **前 $k$ 个特征值之和**：$f(A) = \sum_{i=1}^k \lambda_i^\downarrow(A)$ 在 $S_n$ 上凸；
    3. **核范数**（迹范数）：$f(A) = \|A\|_* = \sum_i \sigma_i(A)$ 在 $M_n$ 上凸；
    4. **谱范数**：$f(A) = \|A\|_2 = \sigma_{\max}(A)$ 在 $M_n$ 上凸；
    5. **负对数行列式**：$f(A) = -\log\det(A)$ 在 $S_n^{++}$ 上凸；
    6. **迹的指数**：$f(A) = \operatorname{tr}(e^A)$ 在 $S_n$ 上凸；
    7. **矩阵幂**（$p \geq 1$）：$f(A) = \operatorname{tr}(A^p)$ 在 $S_n^+$ 上凸。

??? proof "证明"
    **$\lambda_{\max}$ 的凸性**：利用 $\lambda_{\max}(A) = \max_{\|x\|=1} x^TAx$，这是仿射函数 $x^TAx$（关于 $A$）的逐点上确界，因此凸。

    **$-\log\det$ 的凸性**：设 $A, B \in S_n^{++}$，$t \in [0,1]$，$C = tA + (1-t)B$。由 AM-GM 不等式的矩阵版本，
    $$
    \det(C) = \det(tA + (1-t)B) \geq \det(A)^t \det(B)^{1-t}
    $$
    取对数并取负号即得 $-\log\det(C) \leq -t\log\det(A) - (1-t)\log\det(B)$。

    矩阵 AM-GM 的证明：不妨设 $B = I$（通过 $B^{-1/2}$ 进行缩并）。则需证 $\det(tA + (1-t)I) \geq \det(A)^t$。设 $A$ 的特征值为 $\lambda_i > 0$，则
    $$
    \det(tA + (1-t)I) = \prod_i (t\lambda_i + 1 - t) \geq \prod_i \lambda_i^t = \det(A)^t
    $$
    其中使用了标量凹函数 $\log(t\lambda + 1-t) \geq t\log\lambda$。

!!! example "例 64.7"
    **$\lambda_{\max}$ 的凸性示例**：设 $A = \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix}$，$t = 1/2$。

    $\frac{A+B}{2} = \begin{pmatrix} 2 & 0 \\ 0 & 2.5 \end{pmatrix}$，$\lambda_{\max}\left(\frac{A+B}{2}\right) = 2.5$。

    $\frac{\lambda_{\max}(A) + \lambda_{\max}(B)}{2} = \frac{3 + 4}{2} = 3.5$。

    确实 $2.5 \leq 3.5$，验证了凸性。

### Ky Fan 范数的凸性

!!! theorem "定理 64.13 (Ky Fan 范数的凸性)"
    对任意 $k = 1, \ldots, n$，**Ky Fan $k$-范数**
    $$
    \|A\|_{(k)} = \sum_{i=1}^k \sigma_i^\downarrow(A)
    $$
    在 $M_n$ 上是凸函数（这里 $\sigma_1^\downarrow \geq \cdots \geq \sigma_n^\downarrow$ 是奇异值的降序排列）。

    特别地，$k=1$ 给出谱范数（最大奇异值），$k=n$ 给出核范数（所有奇异值之和）。

??? proof "证明"
    利用 Ky Fan 变分刻画：
    $$
    \|A\|_{(k)} = \max\{\operatorname{tr}(U^T A V) : U \in \mathbb{R}^{n \times k}, V \in \mathbb{R}^{n \times k}, U^TU = V^TV = I_k\}
    $$
    这是仿射函数 $\operatorname{tr}(U^T A V)$（关于 $A$）在紧集上的最大值，因此凸。

---

## 64.7 矩阵凹函数

<div class="context-flow" markdown>

**核心问题**：哪些矩阵函数是凹的？凹性在矩阵分析中有何意义？

</div>

### 重要例子

!!! theorem "定理 64.14 (矩阵凹函数的例子)"
    以下函数在指定的域上是凹的：

    1. **最小特征值**：$f(A) = \lambda_{\min}(A)$ 在 $S_n$ 上凹；
    2. **对数行列式**：$f(A) = \log\det(A)$ 在 $S_n^{++}$ 上凹；
    3. **矩阵幂**（$0 \leq p \leq 1$）：$f(A) = \operatorname{tr}(A^p)$ 在 $S_n^+$ 上凹；
    4. **负迹逆**：$f(A) = -\operatorname{tr}(A^{-1})$ 在 $S_n^{++}$ 上凹；
    5. **行列式的 $1/n$ 次幂**：$f(A) = (\det A)^{1/n}$ 在 $S_n^+$ 上凹。

??? proof "证明"
    **$\lambda_{\min}$ 的凹性**：$\lambda_{\min}(A) = \min_{\|x\|=1} x^TAx$，这是仿射函数的逐点下确界，因此凹。

    **$\log\det$ 的凹性**：已在定理 64.12 的证明中给出（$-\log\det$ 是凸的）。

    **$(\det A)^{1/n}$ 的凹性**：这是 Minkowski 行列式定理的推论。对 $A, B \in S_n^+$，$t \in [0,1]$，
    $$
    \det(tA + (1-t)B)^{1/n} \geq t\det(A)^{1/n} + (1-t)\det(B)^{1/n}
    $$

!!! example "例 64.8"
    **$\log\det$ 的凹性示例**：设 $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 4 & 0 \\ 0 & 1 \end{pmatrix}$，$t = 1/2$。

    $\log\det\left(\frac{A+B}{2}\right) = \log\det\begin{pmatrix} 3 & 0 \\ 0 & 2 \end{pmatrix} = \log 6 \approx 1.791$。

    $\frac{\log\det(A) + \log\det(B)}{2} = \frac{\log 6 + \log 4}{2} = \frac{\log 24}{2} \approx 1.589$。

    确实 $1.791 \geq 1.589$，验证了凹性。

### Lieb 凹性定理

!!! theorem "定理 64.15 (Lieb 凹性定理)"
    设 $K$ 是固定矩阵。映射
    $$
    A \mapsto \operatorname{tr}(K^* A^p K A^{1-p})
    $$
    对 $p \in [0, 1]$ 在 $S_n^{++}$ 上关于 $A$ 是凹的。

    更一般地，对固定的 $K$ 和 $0 \leq p, q$ 且 $p + q \leq 1$，映射
    $$
    (A, B) \mapsto \operatorname{tr}(K^* A^p K B^q)
    $$
    在 $S_n^{++} \times S_n^{++}$ 上是联合凹的（jointly concave）。

!!! note "注"
    Lieb 凹性定理（1973）是矩阵分析中最深刻的结果之一。它的推论包括：

    - **强次可加性**（Strong Subadditivity, SSA）：量子条件熵满足 $S(A|B)_\rho + S(A|C)_\rho \leq 0$（对某些特殊态）；
    - **Wigner-Yanase-Dyson 猜想**的证明；
    - **数据处理不等式**在量子信息论中的建立。

---

## 64.8 极小极大与鞍点

<div class="context-flow" markdown>

**核心问题**：在矩阵框架下，什么条件保证极小极大定理成立？

</div>

### Von Neumann 极小极大定理

!!! theorem "定理 64.16 (Von Neumann 极小极大定理)"
    设 $A \in M_{m \times n}(\mathbb{R})$。则
    $$
    \max_{x \in \Delta_m} \min_{y \in \Delta_n} x^T A y = \min_{y \in \Delta_n} \max_{x \in \Delta_m} x^T A y
    $$
    其中 $\Delta_k = \{z \in \mathbb{R}^k : z_i \geq 0, \sum z_i = 1\}$ 是 $k$ 维单纯形（概率向量集）。

    此公共值称为矩阵博弈的**值**（value of the game），记作 $v(A)$。

??? proof "证明"
    不等式 $\max\min \leq \min\max$ 对任意函数成立（弱对偶性）。关键是证明反方向。

    **证明** $\max\min \geq \min\max$：利用线性规划的对偶定理。

    玩家 I 的最优化问题：
    $$
    \max_{x \in \Delta_m} \min_{j=1,\ldots,n} (A^T x)_j = \max\{v : A^T x \geq v \mathbf{1}, \; x \in \Delta_m\}
    $$
    玩家 II 的最优化问题：
    $$
    \min_{y \in \Delta_n} \max_{i=1,\ldots,m} (Ay)_i = \min\{w : Ay \leq w \mathbf{1}, \; y \in \Delta_n\}
    $$

    这两个线性规划互为对偶，由线性规划强对偶定理，$v^* = w^*$。

!!! example "例 64.9"
    **石头-剪刀-布**的收益矩阵：
    $$
    A = \begin{pmatrix} 0 & -1 & 1 \\ 1 & 0 & -1 \\ -1 & 1 & 0 \end{pmatrix}
    $$

    由对称性，最优策略为 $x^* = y^* = (1/3, 1/3, 1/3)$，博弈值为 $v(A) = 0$。

### 矩阵鞍点问题

!!! definition "定义 64.9 (鞍点)"
    函数 $L(x, y)$ 在点 $(x^*, y^*)$ 处有**鞍点**，若
    $$
    L(x, y^*) \leq L(x^*, y^*) \leq L(x^*, y), \quad \forall x \in X, \; \forall y \in Y
    $$

!!! theorem "定理 64.17 (鞍点存在定理)"
    设 $X \subseteq \mathbb{R}^m$，$Y \subseteq \mathbb{R}^n$ 为紧凸集，$L: X \times Y \to \mathbb{R}$ 满足：

    - 对固定的 $y$，$L(\cdot, y)$ 关于 $x$ 是凹的；
    - 对固定的 $x$，$L(x, \cdot)$ 关于 $y$ 是凸的。

    则 $L$ 存在鞍点 $(x^*, y^*)$，且
    $$
    \max_x \min_y L(x, y) = \min_y \max_x L(x, y) = L(x^*, y^*)
    $$

### 在 SDP 对偶性中的应用

!!! theorem "定理 64.18 (SDP 的弱对偶性)"
    考虑半定规划（SDP）原始问题和对偶问题：

    **原始**：$\min \operatorname{tr}(CX)$，$\text{s.t.}$ $\operatorname{tr}(A_i X) = b_i$，$X \succeq 0$

    **对偶**：$\max \sum_i b_i y_i$，$\text{s.t.}$ $\sum_i y_i A_i \preceq C$

    则对偶间隙非负：原始最优值 $\geq$ 对偶最优值。

    在适当的约束规格条件（如 Slater 条件：存在严格可行的 $X \succ 0$）下，强对偶性成立，间隙为零。

??? proof "证明"
    设 $X$ 原始可行，$y$ 对偶可行。则
    $$
    \operatorname{tr}(CX) - \sum_i b_i y_i = \operatorname{tr}(CX) - \sum_i y_i \operatorname{tr}(A_i X)
    = \operatorname{tr}\left(\left(C - \sum_i y_i A_i\right) X\right)
    $$
    由对偶可行性 $C - \sum_i y_i A_i \succeq 0$，且 $X \succeq 0$，由 PSD 锥的自对偶性（定理 64.4），$\operatorname{tr}\left(\left(C - \sum_i y_i A_i\right) X\right) \geq 0$。

!!! example "例 64.10"
    **最大特征值作为 SDP**：$\lambda_{\max}(A)$ 可以表示为
    $$
    \lambda_{\max}(A) = \min\{t : tI - A \succeq 0\}
    $$
    这是一个 SDP。对偶问题为
    $$
    \lambda_{\max}(A) = \max\{\operatorname{tr}(AX) : X \succeq 0, \; \operatorname{tr}(X) = 1\}
    $$
    即在密度矩阵上最大化 $\operatorname{tr}(AX)$，最优解为 $X^* = vv^T$（$v$ 是最大特征值对应的特征向量）。

---

## 本章小结

本章研究了矩阵空间中的凸性理论。主要结果包括：

1. **矩阵空间**在迹内积下成为有限维内积空间，$\mathbb{R}^{n^2}$ 中的凸分析理论可以直接应用。

2. **PSD 锥** $S_n^+$ 是正则锥且自对偶，其极端射线由秩 1 矩阵生成。PSD 锥是半定规划的基础。

3. **Birkhoff 多面体** $\mathcal{B}_n$ 的顶点恰好是置换矩阵（Birkhoff 定理），维数为 $(n-1)^2$。

4. **椭圆体**（相关矩阵集）$\mathcal{E}_n$ 是重要的凸体，在统计和金融中有核心应用。

5. **分离定理**保证了不在闭凸集中的矩阵可以被迹超平面分离，这是 SDP 对偶性和量子纠缠见证的基础。

6. **矩阵凸/凹函数**：$\lambda_{\max}$、$-\log\det$、核范数等是凸的；$\lambda_{\min}$、$\log\det$、$(\det)^{1/n}$ 是凹的。Lieb 凹性定理是量子信息论的基石。

7. **极小极大定理**在矩阵博弈和 SDP 对偶性中有核心应用。

---

## 习题

!!! question "习题 64.1"
    证明 $S_n^+$ 的维数为 $n(n+1)/2$（即它是 $S_n(\mathbb{R})$ 中的满维锥）。

!!! question "习题 64.2"
    给出所有 $3 \times 3$ 双随机矩阵的参数化表示。

!!! question "习题 64.3"
    证明函数 $f(A) = \operatorname{tr}(A^2)$ 在 $S_n$ 上是凸的。

!!! question "习题 64.4"
    证明 $f(A) = -\operatorname{tr}(A^{-1})$ 在 $S_n^{++}$ 上是凹的。（提示：计算 Hessian 或利用凹性的定义。）

!!! question "习题 64.5"
    设 $A \in S_n(\mathbb{R})$，$A \notin S_n^+$。构造 $H \succeq 0$ 使得 $\operatorname{tr}(AH) < 0$。

!!! question "习题 64.6"
    证明 Birkhoff 多面体 $\mathcal{B}_n$ 有 $n!$ 个顶点和 $n^2$ 个面（facet）。

!!! question "习题 64.7"
    证明相关矩阵 $C \in \mathcal{E}_n$ 的非对角元素满足 $|c_{ij}| \leq 1$。

!!! question "习题 64.8"
    证明 Von Neumann 极小极大定理中的弱对偶性：$\max_x \min_y x^TAy \leq \min_y \max_x x^TAy$。

!!! question "习题 64.9"
    设 $f(A) = \lambda_1(A) + \lambda_2(A)$（最大两个特征值之和）。证明 $f$ 在 $S_n$ 上凸。

!!! question "习题 64.10"
    利用 SDP 对偶性，推导 $\lambda_{\min}(A) = \min\{\operatorname{tr}(AX) : X \succeq 0, \operatorname{tr}(X) = 1\}$。
