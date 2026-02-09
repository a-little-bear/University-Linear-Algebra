# 第 64A 章 矩阵空间中的凸集

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 矩阵不等式(Ch18) · 优化基础(Ch25) · 完全正矩阵(Ch45a) · 共正矩阵(Ch45b)

**本章脉络**：矩阵空间维数与迹内积 → PSD 锥（正则锥、自对偶性、极端射线）→ Birkhoff 多面体（Birkhoff 定理、维数）→ 相关矩阵集（椭圆体）→ 完全正锥与共正锥 → 谱面体 → 分离定理 → S-引理 → Von Neumann 极小极大定理 → SDP 对偶性 → 内点法概述

**延伸**：矩阵凸集理论是半定规划(SDP)和量子信息论（纠缠判定、量子态空间的凸结构）的数学基础

</div>

凸性是优化理论的核心概念。当凸分析从 $\mathbb{R}^n$ 推广到矩阵空间 $M_n(\mathbb{R})$ 或 Hermite 矩阵空间 $H_n(\mathbb{C})$ 时，出现了许多具有丰富结构的凸集。半正定锥 $S_n^+$ 作为最重要的矩阵凸锥，其自对偶性使得半定规划成为可能。Birkhoff 多面体连接置换矩阵与双随机矩阵。相关矩阵集（椭圆体）在金融数学和统计学中无处不在。谱面体作为 SDP 的可行集，将凸优化理论推广到矩阵领域。

本章系统研究矩阵空间中的凸集理论，从矩阵空间的线性结构出发，经由重要的凸集实例，到达分离定理、S-引理和 SDP 对偶性。

---

## 64A.1 矩阵空间的线性结构

### 矩阵空间的维数

!!! definition "定义 64A.1 (矩阵空间)"
    以下矩阵空间是有限维实向量空间：

    - $M_n(\mathbb{R})$：全体 $n \times n$ 实矩阵，实维数 $n^2$；
    - $M_{m \times n}(\mathbb{R})$：全体 $m \times n$ 实矩阵，实维数 $mn$；
    - $S_n(\mathbb{R}) = \{A \in M_n(\mathbb{R}) : A = A^T\}$：实对称矩阵，实维数 $\frac{n(n+1)}{2}$；
    - $\mathcal{A}_n(\mathbb{R}) = \{A \in M_n(\mathbb{R}) : A = -A^T\}$：反对称矩阵，实维数 $\frac{n(n-1)}{2}$；
    - $H_n(\mathbb{C}) = \{A \in M_n(\mathbb{C}) : A = A^*\}$：Hermite 矩阵，实维数 $n^2$。

### 迹内积

!!! definition "定义 64A.2 (迹内积)"
    在实矩阵空间 $M_n(\mathbb{R})$ 上，**迹内积**定义为
    $$
    \langle A, B \rangle = \operatorname{tr}(A^T B) = \sum_{i,j} a_{ij} b_{ij}
    $$

    在 Hermite 矩阵空间 $H_n(\mathbb{C})$ 上，迹内积定义为
    $$
    \langle A, B \rangle = \operatorname{tr}(AB)
    $$
    对 Hermite 矩阵，$\operatorname{tr}(AB) \in \mathbb{R}$。

!!! theorem "定理 64A.1 (迹内积的性质)"
    迹内积是一个内积，满足：

    1. **双线性**：$\langle \alpha A + \beta B, C \rangle = \alpha \langle A, C \rangle + \beta \langle B, C \rangle$；
    2. **对称性**：$\langle A, B \rangle = \langle B, A \rangle$；
    3. **正定性**：$\langle A, A \rangle = \operatorname{tr}(A^T A) = \sum_{i,j} a_{ij}^2 \geq 0$，等号当且仅当 $A = 0$。

    对应的范数 $\|A\|_F = \sqrt{\langle A, A \rangle}$ 是 Frobenius 范数。

!!! theorem "定理 64A.2 (对称-反对称正交分解)"
    $M_n(\mathbb{R})$ 在迹内积下是对称矩阵空间和反对称矩阵空间的正交直和：
    $$
    M_n(\mathbb{R}) = S_n(\mathbb{R}) \oplus^{\perp} \mathcal{A}_n(\mathbb{R})
    $$
    投影算子为 $\Pi_S(A) = \frac{A + A^T}{2}$，$\Pi_A(A) = \frac{A - A^T}{2}$。

??? proof "证明"
    设 $S \in S_n(\mathbb{R})$，$T \in \mathcal{A}_n(\mathbb{R})$。则
    $$
    \langle S, T \rangle = \operatorname{tr}(S^T T) = \operatorname{tr}(ST)
    $$
    由 $S = S^T$ 和 $T = -T^T$，有 $\operatorname{tr}(ST) = \operatorname{tr}((ST)^T) = \operatorname{tr}(T^T S^T) = \operatorname{tr}(-TS) = -\operatorname{tr}(ST)$。因此 $\operatorname{tr}(ST) = 0$，即 $S_n \perp \mathcal{A}_n$。

    维数验证：$\frac{n(n+1)}{2} + \frac{n(n-1)}{2} = n^2 = \dim M_n(\mathbb{R})$。$\blacksquare$

---

## 64A.2 PSD 锥

### 定义与正则性

!!! definition "定义 64A.3 (半正定锥)"
    **半正定锥**（PSD cone）定义为
    $$
    S_n^+ = \{A \in S_n(\mathbb{R}) : A \succeq 0\} = \{A \in S_n(\mathbb{R}) : x^T A x \geq 0, \; \forall x\}
    $$
    其内部是**正定锥**：$S_n^{++} = \{A \in S_n(\mathbb{R}) : A \succ 0\} = \operatorname{int}(S_n^+)$。

!!! definition "定义 64A.4 (正则锥)"
    $\mathbb{R}^d$ 中的锥 $K$ 称为**正则锥**（proper cone），若它满足：

    1. **凸性**：$A, B \in K$，$\alpha, \beta \geq 0$ 蕴含 $\alpha A + \beta B \in K$；
    2. **闭性**：$K$ 是闭集；
    3. **尖性**（pointed）：$K \cap (-K) = \{0\}$；
    4. **实性**（solid）：$\operatorname{int}(K) \neq \emptyset$。

!!! theorem "定理 64A.3 (PSD 锥是正则锥)"
    $S_n^+$ 是 $S_n(\mathbb{R})$ 中的正则锥。

??? proof "证明"
    1. **凸性**：若 $A, B \succeq 0$，$\alpha, \beta \geq 0$，则对任意 $x$：
    $$
    x^T(\alpha A + \beta B)x = \alpha x^TAx + \beta x^TBx \geq 0
    $$

    2. **闭性**：$S_n^+ = \bigcap_{x \in \mathbb{R}^n} \{A \in S_n : x^TAx \geq 0\}$。每个 $\{A : x^TAx \geq 0\}$ 是闭半空间（关于 $A$ 是线性不等式），因此 $S_n^+$ 是闭半空间的交集，从而闭。

    3. **尖性**：若 $A \succeq 0$ 且 $-A \succeq 0$，则 $x^TAx = 0$ 对所有 $x$。取 $x = e_i$ 得 $a_{ii} = 0$，取 $x = e_i + e_j$ 得 $a_{ij} + a_{ji} = 0$，由对称性 $a_{ij} = 0$。故 $A = 0$。

    4. **实性**：$I_n \in S_n^{++}$。对任意 $\varepsilon > 0$ 足够小和 $\|B\|_F < \varepsilon$ 的对称矩阵 $B$，$I + B$ 仍正定。因此 $S_n^{++} = \operatorname{int}(S_n^+) \neq \emptyset$。$\blacksquare$

### 自对偶性

!!! definition "定义 64A.5 (对偶锥)"
    锥 $K \subseteq V$（$V$ 为内积空间）的**对偶锥**定义为
    $$
    K^* = \{Y \in V : \langle X, Y \rangle \geq 0, \; \forall X \in K\}
    $$

!!! theorem "定理 64A.4 (PSD 锥的自对偶性)"
    在迹内积下，$S_n^+$ 是自对偶的：$(S_n^+)^* = S_n^+$。

??? proof "证明"
    **$S_n^+ \subseteq (S_n^+)^*$**：设 $A, B \in S_n^+$。需证 $\langle A, B \rangle = \operatorname{tr}(AB) \geq 0$。

    由 $A \succeq 0$，存在 $A^{1/2} \succeq 0$ 使得 $A = A^{1/2} A^{1/2}$。则
    $$
    \operatorname{tr}(AB) = \operatorname{tr}(A^{1/2} A^{1/2} B) = \operatorname{tr}(A^{1/2} B A^{1/2})
    $$
    （利用迹的循环性质）。由 $B \succeq 0$，$A^{1/2} B A^{1/2} \succeq 0$，因此 $\operatorname{tr}(A^{1/2} B A^{1/2}) \geq 0$。

    **$(S_n^+)^* \subseteq S_n^+$**：设 $Y \in (S_n^+)^*$，即 $\operatorname{tr}(XY) \geq 0$ 对所有 $X \succeq 0$ 成立。

    取 $X = vv^T$（秩 1 半正定矩阵，对任意 $v \in \mathbb{R}^n$）。则
    $$
    v^T Y v = \operatorname{tr}(vv^T Y) = \operatorname{tr}(XY) \geq 0
    $$
    因此 $Y \succeq 0$，即 $Y \in S_n^+$。$\blacksquare$

### 极端射线

!!! theorem "定理 64A.5 (PSD 锥的极端射线)"
    $S_n^+$ 的极端射线恰好是秩 1 半正定矩阵 $\{xx^T : x \in \mathbb{R}^n \setminus \{0\}\}$ 生成的射线。

??? proof "证明"
    **秩 1 矩阵生成极端射线**：设 $A = xx^T$（$x \neq 0$），$\operatorname{rank}(A) = 1$。若 $A = B + C$，$B, C \in S_n^+$，则 $B$ 和 $C$ 的列空间都包含在 $A$ 的列空间 $\operatorname{span}\{x\}$ 中。

    证明：由 $A = B + C$ 且 $B, C \succeq 0$，对任意 $y \perp x$：
    $$
    0 = y^T A y = y^T B y + y^T C y
    $$
    由 $y^T B y \geq 0$ 和 $y^T C y \geq 0$，必须 $y^T B y = y^T C y = 0$。由半正定矩阵的性质，$By = Cy = 0$ 对所有 $y \perp x$ 成立。因此 $B$ 和 $C$ 的列空间包含在 $\operatorname{span}\{x\}$ 中。

    从而 $B = \alpha xx^T$，$C = \beta xx^T$，$\alpha + \beta = 1$，$\alpha, \beta \geq 0$。这证明 $A$ 在 $S_n^+$ 的极端射线上。

    **秩 $\geq 2$ 不是极端的**：若 $\operatorname{rank}(A) \geq 2$，设 $A$ 的谱分解 $A = \sum_{i=1}^r \lambda_i u_i u_i^T$（$r \geq 2$，$\lambda_i > 0$）。取 $B = \lambda_1 u_1 u_1^T \in S_n^+$，$C = A - B \in S_n^+$。则 $A = B + C$ 是非平凡凸分解（$B$ 和 $C$ 不在 $A$ 的射线上），因此 $A$ 不在极端射线上。$\blacksquare$

!!! example "例 64A.1"
    $S_2^+$ 中的矩阵 $A = \begin{pmatrix} a & b \\ b & c \end{pmatrix} \succeq 0$ 需满足 $a \geq 0$，$c \geq 0$，$ac - b^2 \geq 0$。

    这是 $(a, b, c)$ 空间中的一个圆锥面。其极端射线由秩 1 矩阵参数化：
    $$
    \begin{pmatrix} \cos^2\theta & \cos\theta\sin\theta \\ \cos\theta\sin\theta & \sin^2\theta \end{pmatrix}, \quad \theta \in [0, \pi)
    $$

---

## 64A.3 Birkhoff 多面体

### 定义与 Birkhoff 定理

!!! definition "定义 64A.6 (双随机矩阵与 Birkhoff 多面体)"
    $n \times n$ 矩阵 $A$ 称为**双随机的**（doubly stochastic），若所有元素非负，且每行、每列之和均为 1。全体 $n \times n$ 双随机矩阵的集合称为 **Birkhoff 多面体** $\mathcal{B}_n$：
    $$
    \mathcal{B}_n = \left\{A \in M_n(\mathbb{R}) : a_{ij} \geq 0, \; \sum_j a_{ij} = 1, \; \sum_i a_{ij} = 1\right\}
    $$

!!! theorem "定理 64A.6 (Birkhoff 定理)"
    $\mathcal{B}_n$ 是凸多面体，其顶点恰好是 $n \times n$ 置换矩阵：
    $$
    \mathcal{B}_n = \operatorname{conv}\{P_\sigma : \sigma \in S_n\}
    $$

??? proof "证明"
    **顶点是置换矩阵**：设 $A$ 是 $\mathcal{B}_n$ 的顶点（极端点）。若 $A$ 不是置换矩阵，则存在某行有至少两个正元素 $a_{i,j_1}, a_{i,j_2} > 0$。

    利用双随机性条件，从 $(i, j_1)$ 出发，列 $j_1$ 中除 $a_{i,j_1}$ 外必有其他正元素（列和为 1），设 $a_{i_2,j_1} > 0$（$i_2 \neq i$）。从 $(i_2, j_1)$ 出发，行 $i_2$ 中除 $a_{i_2,j_1}$ 外必有其他正元素，依此类推。由有限性，最终回到起点，形成一个**交替圈**。

    沿交替圈，交替加减小量 $\varepsilon > 0$，得到两个不同的双随机矩阵 $A + \varepsilon B$ 和 $A - \varepsilon B$。$A = \frac{1}{2}(A+\varepsilon B) + \frac{1}{2}(A-\varepsilon B)$ 是非平凡凸组合，与 $A$ 是极端点矛盾。

    **每个双随机矩阵是置换矩阵的凸组合**：由 Konig 定理（二部图匹配定理），每个双随机矩阵 $A$ 的正元素支撑集包含一个完美匹配，对应置换矩阵 $P$。令 $\alpha = \min_{P_{ij}=1} a_{ij} > 0$，则
    $$
    A = \alpha P + (1-\alpha) \cdot \frac{A - \alpha P}{1 - \alpha}
    $$
    其中 $\frac{A - \alpha P}{1 - \alpha}$ 仍是双随机矩阵（若 $\alpha < 1$），且至少多一个零元素。

    通过至多 $n^2 - 2n + 2$ 步归纳，$A$ 分解为置换矩阵的凸组合。$\blacksquare$

### Birkhoff 多面体的维数

!!! theorem "定理 64A.7 (Birkhoff 多面体的维数)"
    $\dim(\mathcal{B}_n) = (n-1)^2$。

??? proof "证明"
    双随机矩阵满足 $n$ 个行和约束（$\sum_j a_{ij} = 1$）和 $n$ 个列和约束（$\sum_i a_{ij} = 1$）。这 $2n$ 个等式约束中有一个冗余：所有行和之和等于所有列和之和（都等于 $n$）。因此独立等式约束为 $2n - 1$ 个。

    $\mathcal{B}_n$ 所在的仿射子空间维数为 $n^2 - (2n - 1) = (n-1)^2$。可以验证 $\mathcal{B}_n$ 在此仿射子空间中是满维的（即 $\mathcal{B}_n$ 含有一个 $(n-1)^2$ 维的开球），因此 $\dim(\mathcal{B}_n) = (n-1)^2$。$\blacksquare$

!!! example "例 64A.2"
    $n = 2$：$\mathcal{B}_2$ 是 1 维线段。所有 $2 \times 2$ 双随机矩阵为
    $$
    A = \begin{pmatrix} t & 1-t \\ 1-t & t \end{pmatrix}, \quad t \in [0, 1]
    $$
    顶点为 $I_2$（$t=1$）和 $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$（$t=0$）。

!!! example "例 64A.3"
    $n = 3$：$\mathcal{B}_3$ 是 4 维多面体，有 $3! = 6$ 个顶点。

---

## 64A.4 相关矩阵集（椭圆体）

!!! definition "定义 64A.7 (相关矩阵与椭圆体)"
    $n \times n$ 实矩阵 $C$ 是**相关矩阵**（correlation matrix），若 $C \succeq 0$ 且 $c_{ii} = 1$。全体 $n \times n$ 相关矩阵的集合记作 $\mathcal{E}_n$（**椭圆体**, elliptope）。

!!! theorem "定理 64A.8 (椭圆体的凸性与紧性)"
    $\mathcal{E}_n$ 是 $S_n(\mathbb{R})$ 中的凸紧集。

??? proof "证明"
    **凸性**：若 $C_1, C_2 \in \mathcal{E}_n$，$t \in [0,1]$，则 $C = tC_1 + (1-t)C_2$ 满足 $C \succeq 0$（PSD 锥的凸性）且 $c_{ii} = t + (1-t) = 1$。

    **有界性**：由 $c_{ii} = 1$ 和 $C \succeq 0$，Cauchy-Schwarz 不等式给出 $|c_{ij}|^2 \leq c_{ii} c_{jj} = 1$，因此 $|c_{ij}| \leq 1$。$\mathcal{E}_n$ 有界。

    **闭性**：$\mathcal{E}_n = S_n^+ \cap \{C : c_{ii} = 1, \forall i\}$ 是闭集的交集。$\blacksquare$

!!! theorem "定理 64A.9 (椭圆体的极端点)"
    $\mathcal{E}_n$ 的极端点恰好是秩 $r$ 的相关矩阵，其中 $r$ 满足 $r(r+1)/2 \leq n$。特别地，$I_n$ 和所有形如 $vv^T$（$v_i = \pm 1$）的秩 1 相关矩阵是极端点。

!!! example "例 64A.4"
    $\mathcal{E}_2 = \left\{\begin{pmatrix} 1 & r \\ r & 1 \end{pmatrix} : r \in [-1, 1]\right\}$ 是一条线段。

    $\mathcal{E}_3$ 由 $(r_{12}, r_{13}, r_{23}) \in [-1,1]^3$ 中满足 $1 - r_{12}^2 - r_{13}^2 - r_{23}^2 + 2r_{12}r_{13}r_{23} \geq 0$ 的点构成。

---

## 64A.5 完全正锥与共正锥

!!! definition "定义 64A.8 (完全正矩阵锥)"
    **完全正矩阵锥**（completely positive cone）定义为
    $$
    \mathcal{CP}_n = \left\{A \in S_n(\mathbb{R}) : A = \sum_{i=1}^k b_i b_i^T, \; b_i \geq 0\right\} = \operatorname{cone}\{bb^T : b \in \mathbb{R}_+^n\}
    $$
    即可以分解为非负向量的外积之和的矩阵。

!!! definition "定义 64A.9 (共正矩阵锥)"
    **共正矩阵锥**（copositive cone）定义为
    $$
    \mathcal{COP}_n = \{A \in S_n(\mathbb{R}) : x^T A x \geq 0, \; \forall x \geq 0\}
    $$

!!! theorem "定理 64A.10 (完全正锥与共正锥的对偶)"
    在迹内积下，$\mathcal{CP}_n$ 和 $\mathcal{COP}_n$ 互为对偶锥：$\mathcal{CP}_n^* = \mathcal{COP}_n$，$\mathcal{COP}_n^* = \mathcal{CP}_n$。

??? proof "证明"
    设 $A \in \mathcal{CP}_n$，$A = \sum_i b_i b_i^T$（$b_i \geq 0$），$C \in \mathcal{COP}_n$。则
    $$
    \langle A, C \rangle = \operatorname{tr}(AC) = \sum_i \operatorname{tr}(b_i b_i^T C) = \sum_i b_i^T C b_i \geq 0
    $$
    因为 $b_i \geq 0$ 且 $C$ 共正。因此 $\mathcal{COP}_n \subseteq \mathcal{CP}_n^*$。

    反方向需要证明若 $C \notin \mathcal{COP}_n$（即存在 $b \geq 0$ 使得 $b^T C b < 0$），则 $C \notin \mathcal{CP}_n^*$。取 $A = bb^T \in \mathcal{CP}_n$，有 $\langle A, C \rangle = b^T C b < 0$，故 $C \notin \mathcal{CP}_n^*$。$\blacksquare$

!!! note "注"
    完全正矩阵和共正矩阵的更深入讨论见第 45A 章和第 45B 章。判定一个矩阵是否完全正或共正在一般情况下是 NP-困难的，这与 PSD 锥的高效可判定性形成鲜明对比。

---

## 64A.6 谱面体

!!! definition "定义 64A.10 (谱面体)"
    **谱面体**（spectrahedron）是半定规划(SDP)的可行集，即形如
    $$
    \mathcal{S} = \left\{x \in \mathbb{R}^d : A_0 + \sum_{i=1}^d x_i A_i \succeq 0\right\}
    $$
    的集合，其中 $A_0, A_1, \ldots, A_d \in S_n(\mathbb{R})$ 是给定的对称矩阵。

!!! example "例 64A.5"
    - **椭圆**：$\{(x_1, x_2) : \begin{pmatrix} 1 - x_1^2 & -x_1 x_2 \\ -x_1 x_2 & 1 - x_2^2 \end{pmatrix} \succeq 0\}$ 是 $x_1^2 + x_2^2 \leq 1$ 的圆盘。但并非所有谱面体都可以用 LMI 的标准形式表示为简单不等式。
    - **PSD 锥**本身是谱面体（取 $A_i = E_{ij}$ 为矩阵单位基）。
    - **Birkhoff 多面体**是谱面体的投影（但不一定是谱面体本身）。

!!! theorem "定理 64A.11 (谱面体的基本性质)"
    1. 每个谱面体都是凸集；
    2. 每个谱面体都是闭集（若非空）；
    3. 多面体（polyhedron）是谱面体的特例（取对角 LMI）；
    4. 谱面体的投影（称为**谱面体影子**或**投影谱面体**）不一定是谱面体，但仍然是凸半代数集。

!!! note "注"
    谱面体将多面体（线性规划的可行集）推广为包含矩阵不等式的可行集。这种推广保留了凸性，同时大大扩展了可描述的几何对象类别。

---

## 64A.7 分离定理

### Hahn-Banach 型分离

!!! theorem "定理 64A.12 (矩阵空间的分离定理)"
    设 $C \subseteq S_n(\mathbb{R})$ 是闭凸集，$A \notin C$。则存在 $H \in S_n(\mathbb{R})$ 和 $\alpha \in \mathbb{R}$，使得
    $$
    \operatorname{tr}(AH) > \alpha \geq \operatorname{tr}(XH), \quad \forall X \in C
    $$

??? proof "证明"
    由于 $(S_n(\mathbb{R}), \langle \cdot, \cdot \rangle)$ 是有限维内积空间（$\langle X, Y \rangle = \operatorname{tr}(XY)$），这是标准 Hahn-Banach 分离定理的直接推论。

    设 $X^*$ 是 $C$ 中到 $A$ 距离最近的点（闭凸集上的最近点投影存在且唯一）。令 $H = A - X^*$。

    由最近点的变分不等式刻画，对任意 $X \in C$：
    $$
    \langle X - X^*, A - X^* \rangle = \operatorname{tr}((X - X^*)H) \leq 0
    $$
    因此 $\operatorname{tr}(XH) \leq \operatorname{tr}(X^* H)$。

    同时 $\operatorname{tr}(AH) = \operatorname{tr}((A - X^*)H) + \operatorname{tr}(X^* H) = \|H\|_F^2 + \operatorname{tr}(X^* H) > \operatorname{tr}(X^* H)$。

    取 $\alpha = \operatorname{tr}(X^* H)$ 即得结论。$\blacksquare$

### PSD 锥的分离

!!! theorem "定理 64A.13 (PSD 锥的分离)"
    若 $A \in S_n(\mathbb{R})$ 不是半正定的，则存在 $H \succeq 0$，$H \neq 0$，使得 $\operatorname{tr}(AH) < 0$。

??? proof "证明"
    设 $\lambda < 0$ 是 $A$ 的负特征值，$v$ 是对应的单位特征向量。取 $H = vv^T \succeq 0$，则
    $$
    \operatorname{tr}(AH) = \operatorname{tr}(Avv^T) = v^T A v = \lambda < 0 \qquad \blacksquare
    $$

### 量子信息中的纠缠见证

!!! note "注"
    在量子信息论中，量子态是密度矩阵 $\rho \succeq 0$，$\operatorname{tr}(\rho) = 1$。判断二体态 $\rho_{AB}$ 是否可分离等价于判断 $\rho_{AB}$ 是否属于某个凸集（可分离态集合）。分离定理保证每个纠缠态都可以被**纠缠见证**（entanglement witness）$W$ 检测：$\operatorname{tr}(W\rho_{AB}) < 0$，但 $\operatorname{tr}(W\sigma) \geq 0$ 对所有可分离态 $\sigma$ 成立。

---

## 64A.8 S-引理

!!! theorem "定理 64A.14 (S-引理 / S-procedure)"
    设 $A, B \in S_n(\mathbb{R})$。若存在 $x_0$ 使得 $x_0^T B x_0 > 0$（严格可行性条件），则以下两个命题等价：

    1. 对所有 $x \in \mathbb{R}^n$，$x^T B x \geq 0 \Rightarrow x^T A x \geq 0$；
    2. 存在 $\lambda \geq 0$，使得 $A - \lambda B \succeq 0$。

??? proof "证明"
    **(2) $\Rightarrow$ (1)**：若 $A - \lambda B \succeq 0$（$\lambda \geq 0$），则对满足 $x^T B x \geq 0$ 的 $x$：
    $$
    x^T A x = x^T(A - \lambda B)x + \lambda x^T B x \geq 0 + 0 = 0
    $$

    **(1) $\Rightarrow$ (2)**：使用对偶论证。定义集合
    $$
    \mathcal{C} = \{(\operatorname{tr}(AX), \operatorname{tr}(BX)) : X \succeq 0, X \neq 0\}
    $$
    这是 $\mathbb{R}^2$ 中的凸锥（因为 $S_n^+$ 是凸锥，线性映射保持凸性）。

    条件 (1) 等价于：$\mathcal{C}$ 不包含 $\{(a, b) : a < 0, b \geq 0\}$ 中的点（即不存在 $X \succeq 0$ 使得 $\operatorname{tr}(AX) < 0$ 且 $\operatorname{tr}(BX) \geq 0$）。

    严格可行性条件（存在 $x_0$ 使得 $x_0^T B x_0 > 0$，即 $\operatorname{tr}(B \cdot x_0 x_0^T) > 0$）保证 $\mathcal{C}$ 在第二坐标方向上无界。

    由二维凸锥的分离定理，存在直线 $a + \lambda b = 0$（$\lambda \geq 0$）分离 $\mathcal{C}$ 与负半轴。这意味着 $\operatorname{tr}((A - \lambda B)X) \geq 0$ 对所有 $X \succeq 0$ 成立，由 PSD 锥的自对偶性，$A - \lambda B \succeq 0$。$\blacksquare$

!!! example "例 64A.6"
    设 $A = \begin{pmatrix} 1 & 3 \\ 3 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$。

    取 $x_0 = (1, 0)^T$，$x_0^T B x_0 = 1 > 0$。问：$x^T B x \geq 0$ 是否蕴含 $x^T A x \geq 0$？

    由 S-引理，等价于寻找 $\lambda \geq 0$ 使得 $A - \lambda B = \begin{pmatrix} 1-\lambda & 3 \\ 3 & 1+\lambda \end{pmatrix} \succeq 0$。

    需要 $(1-\lambda)(1+\lambda) - 9 \geq 0$，即 $1 - \lambda^2 \geq 9$，这不可能。因此条件 (1) 不成立——存在 $x$ 使得 $x^T B x \geq 0$ 但 $x^T A x < 0$。

---

## 64A.9 Von Neumann 极小极大定理

!!! theorem "定理 64A.15 (Von Neumann 极小极大定理)"
    设 $A \in M_{m \times n}(\mathbb{R})$。则
    $$
    \max_{x \in \Delta_m} \min_{y \in \Delta_n} x^T A y = \min_{y \in \Delta_n} \max_{x \in \Delta_m} x^T A y
    $$
    其中 $\Delta_k = \{z \in \mathbb{R}^k : z_i \geq 0, \sum z_i = 1\}$ 是概率单纯形。公共值为矩阵博弈的**值** $v(A)$。

??? proof "证明"
    **弱对偶**（$\max\min \leq \min\max$）对任意函数成立。关键是反方向。

    将两个优化问题写为线性规划：

    玩家 I：$\max v$，subject to $A^T x \geq v \mathbf{1}$，$x \in \Delta_m$。即
    $$
    \max\{v : A^T x \geq v \mathbf{1}, \; x \geq 0, \; \mathbf{1}^T x = 1\}
    $$

    玩家 II：$\min w$，subject to $Ay \leq w \mathbf{1}$，$y \in \Delta_n$。即
    $$
    \min\{w : Ay \leq w \mathbf{1}, \; y \geq 0, \; \mathbf{1}^T y = 1\}
    $$

    这两个线性规划互为对偶（可以通过标准线性规划对偶性验证）。由线性规划**强对偶定理**，原始最优值等于对偶最优值：$v^* = w^*$。$\blacksquare$

!!! example "例 64A.7"
    石头-剪刀-布的收益矩阵：
    $$
    A = \begin{pmatrix} 0 & -1 & 1 \\ 1 & 0 & -1 \\ -1 & 1 & 0 \end{pmatrix}
    $$
    由对称性，最优策略 $x^* = y^* = (1/3, 1/3, 1/3)$，博弈值 $v(A) = 0$。

---

## 64A.10 SDP 对偶性

### 弱对偶与强对偶

!!! theorem "定理 64A.16 (SDP 弱对偶性)"
    考虑 SDP 原始问题和对偶问题：

    **原始**：$p^* = \min \operatorname{tr}(CX)$，s.t. $\operatorname{tr}(A_i X) = b_i$（$i = 1, \ldots, m$），$X \succeq 0$

    **对偶**：$d^* = \max \sum_i b_i y_i$，s.t. $\sum_i y_i A_i \preceq C$

    则 $d^* \leq p^*$（弱对偶性）。若存在严格可行的 $X \succ 0$（Slater 条件），则 $d^* = p^*$（强对偶性）。

??? proof "证明"
    设 $X$ 原始可行，$y$ 对偶可行。则
    $$
    \operatorname{tr}(CX) - \sum_i b_i y_i = \operatorname{tr}(CX) - \sum_i y_i \operatorname{tr}(A_i X) = \operatorname{tr}\left(\left(C - \sum_i y_i A_i\right) X\right)
    $$
    由对偶可行性 $C - \sum_i y_i A_i \succeq 0$，且 $X \succeq 0$，由 PSD 锥的自对偶性（定理 64A.4）：
    $$
    \operatorname{tr}\left(\left(C - \sum_i y_i A_i\right) X\right) \geq 0
    $$
    因此 $\operatorname{tr}(CX) \geq \sum_i b_i y_i$。取所有可行解的最优值即得 $p^* \geq d^*$。$\blacksquare$

### 内点法概述

!!! definition "定义 64A.11 (SDP 的内点法)"
    **内点法**（interior-point method）通过在 PSD 锥的内部沿**中心路径**（central path）迭代来求解 SDP。对于原始问题，中心路径由参数化问题
    $$
    \min \operatorname{tr}(CX) - \mu \log\det(X), \quad \text{s.t. } \operatorname{tr}(A_i X) = b_i
    $$
    定义，其中 $\mu > 0$ 是障碍参数，$-\log\det(X)$ 是 PSD 锥的**自协调障碍函数**。

    随着 $\mu \to 0^+$，中心路径的极限点收敛到原始问题的最优解。

!!! theorem "定理 64A.17 (SDP 内点法的复杂度)"
    对于 $n \times n$ 的 SDP（$m$ 个约束），内点法在 $O(\sqrt{n} \log(1/\varepsilon))$ 次 Newton 迭代内达到 $\varepsilon$-最优解。每次迭代的计算复杂度主要来自求解 $m \times m$ 线性方程组和矩阵分解，总复杂度为 $O(m^2 n^2 + mn^3)$ 每步。

!!! example "例 64A.8"
    $\lambda_{\max}(A)$ 可表示为 SDP：
    $$
    \lambda_{\max}(A) = \min\{t : tI - A \succeq 0\} = \max\{\operatorname{tr}(AX) : X \succeq 0, \operatorname{tr}(X) = 1\}
    $$
    对偶最优解为 $X^* = vv^T$（$v$ 是最大特征值对应的特征向量）。

---

## 习题

!!! question "习题 64A.1"
    证明 $S_n^+$ 在 $S_n(\mathbb{R})$ 中是满维的锥（即 $\dim(S_n^+) = n(n+1)/2$）。

!!! question "习题 64A.2"
    给出所有 $3 \times 3$ 双随机矩阵的参数化表示。验证维数为 $(3-1)^2 = 4$。

!!! question "习题 64A.3"
    证明 Birkhoff 多面体 $\mathcal{B}_n$ 有 $n!$ 个顶点和 $n^2$ 个面（facet）。

!!! question "习题 64A.4"
    证明相关矩阵 $C \in \mathcal{E}_n$ 的非对角元素满足 $|c_{ij}| \leq 1$。进一步证明：$|c_{ij}| = 1$ 当且仅当 $C$ 的第 $i$ 行与第 $j$ 行（或列）之间存在线性关系。

!!! question "习题 64A.5"
    设 $A \in S_n(\mathbb{R})$，$A \notin S_n^+$。构造 $H \succeq 0$ 使得 $\operatorname{tr}(AH) < 0$，并说明 $H$ 的几何意义。

!!! question "习题 64A.6"
    证明 Von Neumann 极小极大定理中的弱对偶性：对任意函数 $f(x, y)$，$\max_x \min_y f(x, y) \leq \min_y \max_x f(x, y)$。

!!! question "习题 64A.7"
    利用 SDP 对偶性，推导 $\lambda_{\min}(A) = \min\{\operatorname{tr}(AX) : X \succeq 0, \operatorname{tr}(X) = 1\}$。

!!! question "习题 64A.8"
    证明：若 $A \in \mathcal{CP}_n$（完全正矩阵）且 $B \in \mathcal{COP}_n$（共正矩阵），则 $\operatorname{tr}(AB) \geq 0$。

!!! question "习题 64A.9"
    构造一个 $3 \times 3$ 共正但非半正定的矩阵。

!!! question "习题 64A.10"
    设 $\mathcal{S} = \{x \in \mathbb{R}^2 : \begin{pmatrix} 1 + x_1 & x_2 \\ x_2 & 1 - x_1 \end{pmatrix} \succeq 0\}$。证明 $\mathcal{S}$ 是一个圆盘，并求其半径。

!!! question "习题 64A.11"
    证明 S-引理中严格可行性条件不可省略。即给出 $A, B$ 使得条件 (1) 成立但条件 (2) 不成立的例子（其中不存在 $x_0$ 使得 $x_0^T B x_0 > 0$）。

!!! question "习题 64A.12"
    证明谱面体的投影仍然是凸集。给出一个谱面体投影不是谱面体的例子。
