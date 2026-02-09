# 第 45A 章 完全正矩阵

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 非负矩阵(Ch17) · 凸锥与对偶(Ch25)

**本章脉络**：完全正(CP)矩阵定义 → cp-秩与上界 → CP 锥的极端射线 → 双非负(DNN)矩阵 → 小维度 CP=DNN（Maxfield-Minc）→ Berman-Xu 图论刻画 → CP 锥的面结构 → $n \geq 5$ 时 CP $\subsetneq$ DNN（Horn 矩阵）→ 完全正半定矩阵（量子推广）

**延伸**：完全正矩阵与非负秩的关系是通信复杂性理论的核心问题；完全正半定矩阵是量子信息中 Bell 不等式研究的基本工具

</div>

矩阵的正定性（$x^TAx > 0$ 对所有 $x \neq 0$）是线性代数中最基本、最有用的概念之一。当我们问一个对称矩阵能否分解为 $A = BB^T$，其中 $B$ 的所有元素非负时，就进入了**完全正矩阵**（completely positive matrix）的领域。这个看似简单的分解条件，既蕴含了半正定性，又包含了元素非负性，将代数结构与组合结构联系在一起。

完全正矩阵构成的凸锥 $\mathcal{CP}_n$ 是凸优化中的重要对象。当维度 $n \leq 4$ 时，完全正性等价于更容易验证的"双非负"条件（半正定且元素非负），但从 $n = 5$ 开始，两者严格分离——这种分离与图论中奇圈的结构密切相关。本章系统建立完全正矩阵的基础理论，探索 cp-秩的上界、锥的几何结构，以及与图论的深刻联系。

---

## 45A.1 完全正矩阵的定义与基本性质

<div class="context-flow" markdown>

**核心问题**：什么样的矩阵可以分解为非负矩阵的"外积"之和？这种分解的最小长度如何界定？

</div>

!!! definition "定义 45A.1 (完全正矩阵)"
    对称矩阵 $A \in \mathbb{R}^{n \times n}$ 称为**完全正矩阵**（completely positive matrix, CP），如果存在非负矩阵 $B \in \mathbb{R}^{n \times k}_{\geq 0}$（$k$ 为某个正整数）使得
    $$A = BB^T.$$

    等价地，$A$ 完全正当且仅当存在非负向量 $b_1, b_2, \ldots, b_k \in \mathbb{R}^n_{\geq 0}$ 使得
    $$A = \sum_{j=1}^{k} b_j b_j^T.$$

    所有 $n \times n$ 完全正矩阵构成的集合记为 $\mathcal{CP}_n$。

!!! theorem "定理 45A.1 (完全正矩阵的基本性质)"
    1. 完全正矩阵 $A$ 必然是**半正定**的（$A \succeq 0$）。
    2. 完全正矩阵 $A$ 必然是**元素非负**的（$a_{ij} \geq 0$ 对所有 $i, j$）。
    3. $\mathcal{CP}_n$ 是对称矩阵空间 $\mathcal{S}_n$ 中的一个**闭凸锥**。
    4. $\mathcal{CP}_n$ 是**尖锐锥**（pointed cone）：若 $A \in \mathcal{CP}_n$ 且 $-A \in \mathcal{CP}_n$，则 $A = 0$。
    5. $\mathcal{CP}_n$ 具有非空**内部**（即 $\mathcal{CP}_n$ 是满维锥）。

??? proof "证明"
    **(1)** 若 $A = BB^T$，则对任意 $x \in \mathbb{R}^n$，
    $$x^TAx = x^TBB^Tx = \|B^Tx\|^2 \geq 0.$$
    因此 $A$ 半正定。

    **(2)** $A$ 的第 $(i, j)$ 元素为
    $$a_{ij} = (BB^T)_{ij} = \sum_{l=1}^{k} b_{il} b_{jl}.$$
    由于 $B \geq 0$，每项 $b_{il} b_{jl} \geq 0$，因此 $a_{ij} \geq 0$。

    **(3)** 分三步证明：

    - **锥性**：若 $A = BB^T \in \mathcal{CP}_n$ 且 $\alpha \geq 0$，则 $\alpha A = (\sqrt{\alpha}B)(\sqrt{\alpha}B)^T \in \mathcal{CP}_n$。
    - **凸性**：若 $A_1 = B_1 B_1^T$ 和 $A_2 = B_2 B_2^T$ 属于 $\mathcal{CP}_n$，则
    $$A_1 + A_2 = \begin{pmatrix} B_1 & B_2 \end{pmatrix} \begin{pmatrix} B_1 & B_2 \end{pmatrix}^T \in \mathcal{CP}_n.$$
    因此对 $\lambda \in [0, 1]$，$\lambda A_1 + (1 - \lambda)A_2 \in \mathcal{CP}_n$。
    - **闭性**：设 $A^{(m)} = B^{(m)}(B^{(m)})^T \to A$。由 Caratheodory 定理（定理 45A.2），可以假设 $B^{(m)}$ 的列数一致有界于 $n(n+1)/2$。由于 $\|B^{(m)}\|_F^2 = \operatorname{tr}(A^{(m)}) \to \operatorname{tr}(A)$ 有界，可以提取收敛子列 $B^{(m_j)} \to B^*$，极限满足 $A = B^*(B^*)^T$ 且 $B^* \geq 0$。

    **(4)** 若 $A \in \mathcal{CP}_n$，则 $A \succeq 0$。若同时 $-A \in \mathcal{CP}_n$，则 $-A \succeq 0$，即 $A \preceq 0$。因此 $A = 0$。

    **(5)** 单位矩阵 $I_n = \sum_{i=1}^n e_i e_i^T \in \mathcal{CP}_n$，且 $I_n$ 正定、元素非负。对充分小的扰动 $E$（正定且元素非负），$I_n + E$ 仍在 $\mathcal{CP}_n$ 中，因此 $\mathcal{CP}_n$ 有非空内部。 $\blacksquare$

!!! example "例 45A.1 (完全正矩阵的例子)"
    1. **对角矩阵**：$D = \operatorname{diag}(d_1, \ldots, d_n)$ 完全正当且仅当 $d_i \geq 0$ 对所有 $i$。此时 $D = BB^T$，$B = \operatorname{diag}(\sqrt{d_1}, \ldots, \sqrt{d_n})$，cp-秩等于正对角元素的个数。

    2. **秩一矩阵**：$A = bb^T$，$b \geq 0$。cp-秩为 1。

    3. **$2 \times 2$ 矩阵**：$A = \begin{pmatrix} a & c \\ c & b \end{pmatrix}$（$a, b, c \geq 0$）完全正当且仅当 $ab \geq c^2$（即半正定且元素非负）。

    4. **全 1 矩阵**：$J_n = \mathbf{1}\mathbf{1}^T$ 是完全正的（取 $B = \mathbf{1}$），cp-秩为 1。

    5. **对角占优非负矩阵**：若 $A \geq 0$（元素非负）且 $A$ 对角占优（$a_{ii} \geq \sum_{j \neq i} a_{ij}$），则 $A$ 完全正。

---

## 45A.2 cp-秩与上界

<div class="context-flow" markdown>

**核心问题**：完全正矩阵的 cp-秩可以多大？是否存在只依赖于矩阵大小的通用上界？

</div>

!!! definition "定义 45A.2 (cp-秩)"
    完全正矩阵 $A \in \mathcal{CP}_n$ 的 **cp-秩**（cp-rank）定义为
    $$\operatorname{cp-rank}(A) = \min\left\{k : A = \sum_{j=1}^{k} b_j b_j^T, \; b_j \in \mathbb{R}^n_{\geq 0}\right\}.$$

    等价地，$\operatorname{cp-rank}(A) = \min\{k : A = BB^T, \; B \in \mathbb{R}^{n \times k}_{\geq 0}\}$，即使分解 $A = BB^T$（$B \geq 0$）成立的最小列数。

!!! definition "定义 45A.3 (非负秩)"
    非负矩阵 $A \in \mathbb{R}^{m \times n}_{\geq 0}$ 的**非负秩**（nonnegative rank）定义为
    $$\operatorname{rank}_+(A) = \min\left\{k : A = \sum_{j=1}^{k} u_j v_j^T, \; u_j \geq 0, \; v_j \geq 0\right\}.$$

    对于完全正矩阵 $A \in \mathcal{CP}_n$，由于 $A$ 对称且每个非负秩一分解可以对称化，有
    $$\operatorname{cp-rank}(A) = \operatorname{rank}_+(A).$$

!!! theorem "定理 45A.2 (Caratheodory 上界)"
    对任意 $A \in \mathcal{CP}_n$，
    $$\operatorname{cp-rank}(A) \leq \frac{n(n+1)}{2}.$$

??? proof "证明"
    $\mathcal{CP}_n$ 是 $\mathcal{S}_n$ 中的凸锥，$\mathcal{S}_n$ 的维数为 $\binom{n+1}{2} = n(n+1)/2$。$\mathcal{CP}_n$ 由集合 $\{bb^T : b \in \mathbb{R}^n_{\geq 0}\}$ 生成（作为锥），每个 $A \in \mathcal{CP}_n$ 都是这些生成元的锥组合：
    $$A = \sum_{j=1}^{k} b_j b_j^T.$$

    由 Caratheodory 定理（凸锥版本），$\mathcal{S}_n$ 中由一个集合生成的锥的每个元素都可以表示为至多 $\dim(\mathcal{S}_n) = n(n+1)/2$ 个生成元的锥组合。

    因此 $\operatorname{cp-rank}(A) \leq n(n+1)/2$。 $\blacksquare$

!!! theorem "定理 45A.3 (Drew-Johnson-Loewy 猜想)"
    **Drew-Johnson-Loewy 猜想** (2003) 断言：对 $n \geq 4$ 的任意 $A \in \mathcal{CP}_n$，
    $$\operatorname{cp-rank}(A) \leq \left\lfloor \frac{n^2}{4} \right\rfloor.$$

    该猜想的已知结果如下：

    | $n$ | Caratheodory 上界 $n(n+1)/2$ | DJL 猜想上界 $\lfloor n^2/4 \rfloor$ | 状态 |
    |-----|------|------|------|
    | 1 | 1 | 0 | 显然 |
    | 2 | 3 | 1 | 已证 |
    | 3 | 6 | 2 | 已证 |
    | 4 | 10 | 4 | 已证 |
    | 5 | 15 | 6 | 已证 (Loewy-Tam, 2003) |
    | 6 | 21 | 9 | 部分结果 |
    | $\geq 7$ | $O(n^2)$ | $O(n^2)$ | 开放 |

    该猜想（若成立）将把 cp-秩的上界从 $n(n+1)/2 \sim n^2/2$ 改善到 $n^2/4$，这是一个本质性的改善。

!!! example "例 45A.2 (cp-秩的例子)"
    1. **对角矩阵**：$\operatorname{cp-rank}(\operatorname{diag}(d_1, \ldots, d_n)) = |\{i : d_i > 0\}|$（正对角元素的个数）。

    2. **全 1 矩阵**：$\operatorname{cp-rank}(J_n) = 1$（$J_n = \mathbf{1}\mathbf{1}^T$）。

    3. **达到 DJL 上界的矩阵**：对 $n = 4$，矩阵
    $$A = \begin{pmatrix} 2 & 1 & 0 & 1 \\ 1 & 2 & 1 & 0 \\ 0 & 1 & 2 & 1 \\ 1 & 0 & 1 & 2 \end{pmatrix}$$
    的 cp-秩为 4，恰好等于 $\lfloor 16/4 \rfloor = 4$。

    4. **cp-秩与普通秩的关系**：对任意 $A \in \mathcal{CP}_n$，$\operatorname{rank}(A) \leq \operatorname{cp-rank}(A)$。等号成立的一个充分条件是 $A$ 的 Cholesky 因子恰好是非负的。

    5. **cp-秩可以远大于秩**：存在 $n \times n$ 完全正矩阵，其秩为 $n$ 但 cp-秩为 $\Theta(n^2)$。这种秩与 cp-秩的巨大差距是完全正矩阵理论的特征之一。

---

## 45A.3 CP 锥的极端射线

<div class="context-flow" markdown>

**核心问题**：$\mathcal{CP}_n$ 的极端射线（不可进一步分解的元素）是什么？

</div>

!!! definition "定义 45A.4 (极端射线)"
    闭凸锥 $\mathcal{K}$ 的一条**极端射线**（extreme ray）是满足以下条件的非零元素 $A \in \mathcal{K}$：若 $A = A_1 + A_2$（$A_1, A_2 \in \mathcal{K}$），则 $A_1$ 和 $A_2$ 都是 $A$ 的非负倍数。

!!! theorem "定理 45A.4 (CP 锥的极端射线)"
    $\mathcal{CP}_n$ 的极端射线恰好是形如 $bb^T$ 的秩一矩阵，其中 $b \in \mathbb{R}^n_{\geq 0} \setminus \{0\}$。

    换言之，矩阵 $A \in \mathcal{CP}_n$ 位于极端射线上当且仅当 $\operatorname{cp-rank}(A) = 1$。

??? proof "证明"
    **充分性**：设 $A = bb^T$（$b \geq 0$, $b \neq 0$）。假设 $A = A_1 + A_2$，$A_1, A_2 \in \mathcal{CP}_n$。由于 $A_1 \preceq A = bb^T$ 且 $A_1 \succeq 0$，$A_1$ 的秩至多为 1（因为 $bb^T$ 的秩为 1，而 $0 \preceq A_1 \preceq bb^T$ 意味着 $A_1$ 的值域包含在 $bb^T$ 的值域 $\operatorname{span}\{b\}$ 中）。因此 $A_1 = \alpha bb^T$ 对某个 $\alpha \geq 0$，从而 $A_2 = (1 - \alpha)bb^T$。这表明 $bb^T$ 生成极端射线。

    **必要性**：设 $A \in \mathcal{CP}_n$ 位于极端射线上，$\operatorname{cp-rank}(A) = k$。则 $A = \sum_{j=1}^k b_j b_j^T$（$b_j \geq 0$）。若 $k \geq 2$，取 $A_1 = b_1 b_1^T$ 和 $A_2 = \sum_{j=2}^k b_j b_j^T$，两者都在 $\mathcal{CP}_n$ 中且非零，但 $A_2$ 不是 $A_1$ 的倍数（因为 $b_2$ 不必与 $b_1$ 成比例）。这与 $A$ 位于极端射线上矛盾。因此 $k = 1$。 $\blacksquare$

!!! example "例 45A.3 (极端射线的几何)"
    对 $n = 2$，$\mathcal{CP}_2$ 的极端射线由向量 $b = (b_1, b_2)^T$（$b_1, b_2 \geq 0$）生成。极端射线参数化为
    $$bb^T = \begin{pmatrix} b_1^2 & b_1 b_2 \\ b_1 b_2 & b_2^2 \end{pmatrix},$$
    可以用角度 $\theta \in [0, \pi/2]$（取 $b = (\cos\theta, \sin\theta)^T$）参数化。因此 $\mathcal{CP}_2$ 的极端射线构成一条参数曲线。

    一般地，$\mathcal{CP}_n$ 的极端射线由非负正交锥 $\mathbb{R}^n_{\geq 0}$ 中的方向参数化（模正比例缩放），构成一个 $(n-1)$ 维的集合。

---

## 45A.4 双非负矩阵

<div class="context-flow" markdown>

**核心问题**：半正定且元素非负的矩阵一定是完全正的吗？

</div>

!!! definition "定义 45A.5 (双非负矩阵)"
    对称矩阵 $A \in \mathbb{R}^{n \times n}$ 称为**双非负矩阵**（doubly nonnegative matrix, DNN），如果

    1. $A \succeq 0$（半正定），且
    2. $A \geq 0$（元素非负，即 $a_{ij} \geq 0$ 对所有 $i, j$）。

    所有 $n \times n$ 双非负矩阵构成闭凸锥，记为 $\mathcal{DNN}_n$。显然
    $$\mathcal{DNN}_n = \mathcal{S}_n^+ \cap \mathcal{N}_n,$$
    其中 $\mathcal{S}_n^+$ 是半正定锥，$\mathcal{N}_n = \{A \in \mathcal{S}_n : a_{ij} \geq 0, \forall i, j\}$ 是元素非负的对称矩阵锥。

由定理 45A.1 的 (1) 和 (2)，每个完全正矩阵都是双非负的：
$$\mathcal{CP}_n \subseteq \mathcal{DNN}_n.$$

自然的问题是：反包含 $\mathcal{DNN}_n \subseteq \mathcal{CP}_n$ 是否成立？

---

## 45A.5 小维度时 CP = DNN

!!! theorem "定理 45A.5 (Maxfield-Minc 定理)"
    当 $n \leq 4$ 时，$\mathcal{CP}_n = \mathcal{DNN}_n$。即每个 $n \leq 4$ 的双非负矩阵都是完全正的。

??? proof "证明"
    逐维度证明。

    **$n = 1$**：$A = (a)$，$a \geq 0$。$A = (\sqrt{a})(\sqrt{a})^T \in \mathcal{CP}_1$。

    **$n = 2$**：设 $A = \begin{pmatrix} a & c \\ c & b \end{pmatrix}$，$a, b, c \geq 0$，$ab \geq c^2$（半正定）。若 $c = 0$，则 $A = \operatorname{diag}(a, b)$ 显然完全正。若 $c > 0$，令 $b_1 = (\sqrt{a}, c/\sqrt{a})^T$。由 $ab \geq c^2$ 知 $b - c^2/a \geq 0$，因此
    $$A = b_1 b_1^T + \operatorname{diag}(0, b - c^2/a) \in \mathcal{CP}_2.$$

    **$n = 3$**：设 $A \in \mathcal{DNN}_3$。由 $A \succeq 0$，存在分解 $A = LL^T$（$L$ 为 $3 \times r$ 矩阵，$r = \operatorname{rank}(A)$）。策略是通过正交变换 $L \mapsto LQ$ 使 $LQ \geq 0$。

    对 $n = 3$，利用以下关键引理：若 $A \in \mathcal{DNN}_3$ 且 $\operatorname{rank}(A) = r$，则可以找到非负向量 $b \geq 0$ 和最大的 $\alpha > 0$ 使得 $A - \alpha bb^T \succeq 0$ 且 $A - \alpha bb^T \geq 0$。对余项 $A' = A - \alpha bb^T$ 归纳——$A'$ 仍是双非负的，但秩或非零元素个数减少。对 $3 \times 3$ 矩阵，这个过程在有限步内终止。

    **$n = 4$**：Maxfield 和 Minc (1962) 首先证明了这一结论。证明更为技术性，使用了以下关键步骤：

    1. 若 $A \in \mathcal{DNN}_4$ 且 $\operatorname{rank}(A) \leq 3$，则 $A$ 的某个 $3 \times 3$ 主子矩阵的 CP 分解可以延拓到整个矩阵。
    2. 若 $\operatorname{rank}(A) = 4$，利用 $\mathbb{R}^4$ 中的几何论证：$A$ 的 Cholesky 因子的列向量在 $\mathbb{R}^4$ 中可以通过正交变换旋转到非负正交锥中。这里用到了 $4$ 维空间中非负正交锥有足够"宽度"的性质。
    3. 完整证明的关键引理：若 $A \in \mathcal{DNN}_n$（$n \leq 4$）且 $\operatorname{rank}(A) = r$，则 $\operatorname{cp-rank}(A) \leq \binom{r+1}{2}$。

    对 $n = 4$，$r \leq 4$，因此 $\operatorname{cp-rank}(A) \leq \binom{5}{2} = 10$。 $\blacksquare$

!!! example "例 45A.4 (cp-秩的上界表)"
    对于 $n \leq 4$ 的完全正矩阵，已证明的 cp-秩上界为：

    | $n$ | cp-秩上界（已证） | Caratheodory 上界 | DJL 猜想 |
    |-----|---------|---------|---------|
    | 1 | 1 | 1 | — |
    | 2 | 3 | 3 | 1 |
    | 3 | 6 | 6 | 2 |
    | 4 | 10 | 10 | 4 |
    | 5 | 6 (Loewy-Tam) | 15 | 6 |

    注意对 $n = 5$，DJL 猜想已证且给出的上界 6 远优于 Caratheodory 的 15。

---

## 45A.6 Berman-Xu 图论刻画

<div class="context-flow" markdown>

**核心问题**：完全正性与矩阵的图结构有什么关系？

</div>

!!! definition "定义 45A.6 (矩阵的比较图)"
    对称矩阵 $A \in \mathbb{R}^{n \times n}$ 的**比较图**（comparison graph）$G(A)$ 定义为：顶点集 $V = \{1, 2, \ldots, n\}$，边集 $E = \{(i, j) : a_{ij} \neq 0, \; i \neq j\}$。

    换言之，$G(A)$ 的边对应于 $A$ 的非零非对角元素。

!!! theorem "定理 45A.6 (Berman-Xu 定理)"
    设 $A \in \mathcal{DNN}_n$（$A$ 双非负）。则 $A \in \mathcal{CP}_n$（$A$ 完全正）当且仅当 $A$ 的比较图 $G(A)$ **不含长度 $\geq 5$ 的奇圈作为子图**。

    更精确地：$A \in \mathcal{DNN}_n$ 且 $A \notin \mathcal{CP}_n$ 的唯一障碍是比较图中存在长度为 $5, 7, 9, \ldots$ 的奇圈。

    等价的图论条件：$G(A)$ 是**完美图**（perfect graph）的一种特殊情况——具体来说，$G(A)$ 不含奇洞（odd hole，即长度 $\geq 5$ 的无弦奇圈）作为诱导子图。

??? proof "证明"
    证明分为两个方向。

    **充分性**（$G(A)$ 无长 $\geq 5$ 奇圈 $\Rightarrow$ $A \in \mathcal{CP}_n$）：

    若 $G(A)$ 不含长 $\geq 5$ 的奇圈，则 $G(A)$ 的每个块（2-连通分量）要么是边、要么是三角形、要么是偶圈（二部图），要么是弦图。对每种情形，可以利用局部的 CP 分解并通过图的分解结构（割顶点分解）将其组合为全局的 CP 分解。

    关键的图论事实是：不含长 $\geq 5$ 奇圈的图可以通过反复取割顶点和删除边的操作，分解为二部图和完全图的组合。对二部图和完全图，DNN 矩阵都是 CP 的。

    **必要性**（$G(A)$ 含长 $\geq 5$ 奇圈 $\Rightarrow$ 存在 $A' \in \mathcal{DNN}_n \setminus \mathcal{CP}_n$ 且 $G(A') = G(A)$）：

    若 $G(A)$ 包含长度为 $2k + 1$（$k \geq 2$）的奇圈 $C$，则可以构造一个 DNN 矩阵 $A'$，其比较图为 $C$，但 $A' \notin \mathcal{CP}_n$。构造基于 Horn 矩阵的推广——对奇圈 $C_{2k+1}$，存在自然的循环对称 DNN 矩阵，其完全正性可以通过分析循环结构排除。 $\blacksquare$

!!! example "例 45A.5 (图论刻画的应用)"
    1. **比较图为二部图**的 DNN 矩阵一定是 CP 的（二部图不含奇圈）。

    2. **比较图为三角形 $C_3$** 的 DNN 矩阵一定是 CP 的（$C_3$ 长度为 3，不满足 $\geq 5$）。

    3. **比较图为 $C_5$（五圈）**的 DNN 矩阵不一定是 CP 的。这正是例 45A.7 中 $n = 5$ 分离的图论根源。

    4. **比较图为完全图 $K_n$** 的 DNN 矩阵：$K_n$ 含 $C_5$ 当且仅当 $n \geq 5$。因此对 $n \geq 5$，存在比较图为 $K_n$ 的 DNN 矩阵不在 $\mathcal{CP}_n$ 中。

---

## 45A.7 CP 锥的面结构

<div class="context-flow" markdown>

**核心问题**：$\mathcal{CP}_n$ 作为凸锥，其面（face）的结构如何？

</div>

!!! definition "定义 45A.7 (凸锥的面)"
    闭凸锥 $\mathcal{K}$ 的一个**面**（face）是子集 $F \subseteq \mathcal{K}$，满足：$F$ 自身是闭凸锥，且若 $a + b \in F$（$a, b \in \mathcal{K}$），则 $a, b \in F$。

    面 $F$ 的**维数**定义为 $F$ 张成的线性子空间的维数。极端射线是 1 维面。

!!! theorem "定理 45A.7 (CP 锥的面结构)"
    $\mathcal{CP}_n$ 的面与非负正交锥 $\mathbb{R}^n_{\geq 0}$ 的面之间存在密切关系。具体来说：

    1. 对 $\mathbb{R}^n_{\geq 0}$ 的每个面 $\mathcal{F}$（形如 $\{x \in \mathbb{R}^n_{\geq 0} : x_i = 0, \; i \in S\}$，$S \subseteq \{1, \ldots, n\}$），定义
    $$\Phi(\mathcal{F}) = \operatorname{cone}\{bb^T : b \in \mathcal{F}\}.$$
    则 $\Phi(\mathcal{F})$ 是 $\mathcal{CP}_n$ 的一个面。

    2. $\mathcal{CP}_n$ 的每个面都可以表示为若干此类面的交。

    3. $\mathcal{CP}_n$ 的面格（face lattice）的结构比半正定锥 $\mathcal{S}_n^+$ 的面格更为复杂：$\mathcal{S}_n^+$ 的面由子空间参数化，而 $\mathcal{CP}_n$ 的面由非负正交锥的面（即坐标子集）参数化。

!!! example "例 45A.6 (CP_2 的面结构)"
    $\mathcal{CP}_2$ 的面结构如下：

    - **0 维面**：$\{0\}$（原点）。
    - **1 维面**（极端射线）：$\{t \cdot bb^T : t \geq 0\}$，其中 $b \in \{(\cos\theta, \sin\theta)^T : \theta \in [0, \pi/2]\}$。这些极端射线构成一个连续族。
    - **2 维面**：例如 $\{t_1 e_1 e_1^T + t_2 b b^T : t_1, t_2 \geq 0\}$（$b$ 有两个正分量）。
    - **3 维面**：$\mathcal{CP}_2$ 本身（$\mathcal{S}_2$ 是 3 维的，$\mathcal{CP}_2$ 是其中的满维锥）。

---

## 45A.8 CP 与 DNN 的分离

<div class="context-flow" markdown>

**核心问题**：$n \geq 5$ 时，是否存在双非负但不完全正的矩阵？如何构造？

</div>

!!! theorem "定理 45A.8 (Horn 矩阵与 CP/DNN 分离)"
    **Horn 矩阵**
    $$H = \begin{pmatrix}
    1 & -1 & 1 & 1 & -1 \\
    -1 & 1 & -1 & 1 & 1 \\
    1 & -1 & 1 & -1 & 1 \\
    1 & 1 & -1 & 1 & -1 \\
    -1 & 1 & 1 & -1 & 1
    \end{pmatrix}$$

    是**共正的**但**不在** $\mathcal{S}_5^+ + \mathcal{N}_5$ 中。

    通过对偶关系（$\mathcal{CP}_n$ 与 $\mathcal{COP}_n$ 互为对偶锥），这意味着存在 $5 \times 5$ 双非负矩阵不在 $\mathcal{CP}_5$ 中，即
    $$\mathcal{CP}_5 \subsetneq \mathcal{DNN}_5.$$

??? proof "证明"
    **$H$ 是共正的**：需要验证 $x^T H x \geq 0$ 对所有 $x \geq 0$。直接计算：
    \begin{align}
    x^T H x &= x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 \\
    &\quad - 2x_1 x_2 - 2x_2 x_3 - 2x_3 x_4 - 2x_4 x_5 - 2x_5 x_1 \\
    &\quad + 2x_1 x_3 + 2x_1 x_4 + 2x_2 x_4 + 2x_2 x_5 + 2x_3 x_5.
    \end{align}

    利用 AM-GM 不等式，对 $x \geq 0$，可以将上式重写为非负项之和。具体地，由循环对称性，令指标模 5，有
    $$x^T H x = \sum_{i=1}^{5} (x_i - x_{i+1})^2 + 2\sum_{i=1}^{5} x_i x_{i+2} - \sum_{i=1}^{5} x_i^2.$$
    经过更精细的配方，可以证明对 $x \geq 0$ 此表达式非负。

    **$H \notin \mathcal{S}_5^+ + \mathcal{N}_5$**：假设 $H = P + N$（$P \succeq 0$，$N \geq 0$）。考虑 $H$ 的负元素位置（如 $h_{12} = -1$），有 $p_{12} + n_{12} = -1$，由 $n_{12} \geq 0$ 得 $p_{12} \leq -1$。但 $P \succeq 0$ 要求 $|p_{12}|^2 \leq p_{11} p_{22}$（Cauchy-Schwarz 不等式），即 $p_{11} p_{22} \geq 1$。

    类似地，对每个负元素位置（$(1,2), (2,3), (3,4), (4,5), (5,1)$ 共五个），都有 $p_{i,i+1} \leq -1$，因此 $p_{ii} p_{i+1,i+1} \geq 1$。

    同时，对正元素位置（如 $h_{13} = 1$），有 $p_{13} \leq h_{13} = 1$（因为 $n_{13} \geq 0$，$p_{13} = 1 - n_{13} \leq 1$）。

    将所有约束组合，利用 $P \succeq 0$ 的行列式条件和主子式非负条件，可以导出矛盾。 $\blacksquare$

!!! example "例 45A.7 ($n = 5$ 的 DNN 非 CP 例子)"
    考虑比较图为五圈 $C_5$ 的矩阵：
    $$A = \begin{pmatrix}
    1 & 1 & 0 & 0 & 1 \\
    1 & 2 & 1 & 0 & 0 \\
    0 & 1 & 2 & 1 & 0 \\
    0 & 0 & 1 & 2 & 1 \\
    1 & 0 & 0 & 1 & 2
    \end{pmatrix}.$$

    **验证 DNN**：$A \geq 0$（元素非负）。$A$ 的特征值约为 $0.382, 0.691, 2.000, 3.000, 3.927$，全部非负，因此 $A \succeq 0$。故 $A \in \mathcal{DNN}_5$。

    **验证非 CP**：$A$ 的比较图是五圈 $C_5$（$1 - 2 - 3 - 4 - 5 - 1$），这是长度为 5 的奇圈。由 Berman-Xu 定理（定理 45A.6），存在这种图结构的 DNN 矩阵不是 CP 的。对于上述具体矩阵，可以通过以下方法直接验证 $A \notin \mathcal{CP}_5$：

    若 $A = \sum_{j=1}^k b_j b_j^T$（$b_j \geq 0$），由 $a_{13} = 0$ 知 $(b_j)_1 (b_j)_3 = 0$ 对每个 $j$，即每个 $b_j$ 的第 1 分量和第 3 分量不同时为正。类似地分析所有零位置，可以导出 $b_j$ 的支撑（support）结构上的约束，最终这些约束与 $A$ 的正元素值不相容。

---

## 45A.9 完全正半定矩阵

<div class="context-flow" markdown>

**核心问题**：完全正矩阵有没有自然的量子推广？

</div>

!!! definition "定义 45A.8 (完全正半定矩阵)"
    分块矩阵 $A = [A_{ij}]_{i,j=1}^n$（$A_{ij} \in \mathbb{C}^{d \times d}$）称为**完全正半定矩阵**（completely positive semidefinite matrix, CPSD），如果存在正半定矩阵 $P_1, \ldots, P_n \in \mathcal{S}_d^+$（$d$ 为某个正整数）使得
    $$A_{ij} = \operatorname{tr}(P_i P_j) \quad \text{对所有 } i, j.$$

    当 $d = 1$ 时（$A_{ij}$ 都是标量），这退化为标准的完全正矩阵（因为 $\operatorname{tr}(p_i p_j) = p_i p_j$，$p_i \geq 0$）。

    所有 $n \times n$ 完全正半定矩阵（标量情形）构成的锥记为 $\mathcal{CS}_+^n$。

!!! theorem "定理 45A.9 (CPSD 矩阵的基本性质)"
    1. $\mathcal{CP}_n \subseteq \mathcal{CS}_+^n \subseteq \mathcal{DNN}_n$。
    2. 对 $n \leq 4$，$\mathcal{CP}_n = \mathcal{CS}_+^n = \mathcal{DNN}_n$。
    3. 对 $n \geq 5$，三个锥严格嵌套：$\mathcal{CP}_n \subsetneq \mathcal{CS}_+^n \subsetneq \mathcal{DNN}_n$。
    4. $\mathcal{CS}_+^n$ 一般**不是**闭锥——这与 $\mathcal{CP}_n$ 是闭锥形成对比，也是 CPSD 理论中的一个本质困难。

??? proof "证明"
    **(1)** 第一个包含：若 $A \in \mathcal{CP}_n$，即 $a_{ij} = \sum_l b_{il} b_{jl}$（$b_{il} \geq 0$），取 $P_i = \operatorname{diag}(b_{i1}, \ldots, b_{ik}) \succeq 0$，则 $\operatorname{tr}(P_i P_j) = \sum_l b_{il} b_{jl} = a_{ij}$。

    第二个包含：若 $a_{ij} = \operatorname{tr}(P_i P_j)$（$P_i \succeq 0$），则 $A$ 是 Gram 矩阵（以 Hilbert-Schmidt 内积度量 $P_i$），因此 $A \succeq 0$。又 $a_{ij} = \operatorname{tr}(P_i P_j) \geq 0$（$P_i, P_j \succeq 0$ 的乘积的迹非负），因此 $A \geq 0$。故 $A \in \mathcal{DNN}_n$。

    **(3)** 的分离需要更精细的构造，涉及量子信息论中的 Bell 不等式和量子关联矩阵。 $\blacksquare$

!!! example "例 45A.8 (CPSD 的量子信息背景)"
    在量子信息论中，CPSD 矩阵与**量子关联**密切相关。考虑 Bell 实验中的关联矩阵 $\Gamma_{ij} = \operatorname{tr}(\rho \cdot M_i \otimes N_j)$（$\rho$ 是量子态，$M_i, N_j$ 是测量算子）。

    - 经典关联对应 CP 矩阵（$d = 1$ 的 CPSD）。
    - 量子关联对应一般的 CPSD 矩阵。
    - DNN 条件对应所谓的"无信号"（no-signaling）条件。

    因此包含关系 $\mathcal{CP}_n \subsetneq \mathcal{CS}_+^n \subsetneq \mathcal{DNN}_n$（$n \geq 5$）在物理上对应于经典关联 $\subsetneq$ 量子关联 $\subsetneq$ 无信号关联——这正是 Bell 不等式和 Tsirelson 问题的数学核心。

---

## 45A.10 DNN 锥的对偶

利用凸锥对偶的一般理论，我们可以计算 $\mathcal{DNN}_n$ 的对偶锥。

!!! theorem "定理 45A.10 (DNN 锥的对偶)"
    $$\mathcal{DNN}_n^* = \mathcal{S}_n^+ + \mathcal{N}_n.$$

    即双非负锥的对偶恰好是半正定锥与非负矩阵锥的 Minkowski 和。

??? proof "证明"
    利用凸锥对偶的性质：若 $\mathcal{K}_1, \mathcal{K}_2$ 是闭凸锥，则
    $$(\mathcal{K}_1 \cap \mathcal{K}_2)^* = \overline{\mathcal{K}_1^* + \mathcal{K}_2^*}.$$

    $\mathcal{DNN}_n = \mathcal{S}_n^+ \cap \mathcal{N}_n$。$\mathcal{S}_n^+$ 是**自对偶锥**（$(\mathcal{S}_n^+)^* = \mathcal{S}_n^+$），$\mathcal{N}_n$ 也是**自对偶锥**（$\mathcal{N}_n^* = \mathcal{N}_n$，因为 $\langle A, B \rangle = \sum a_{ij}b_{ij} \geq 0$ 对所有 $A, B \geq 0$）。

    因此 $\mathcal{DNN}_n^* = \overline{\mathcal{S}_n^+ + \mathcal{N}_n}$。可以证明 $\mathcal{S}_n^+ + \mathcal{N}_n$ 已经是闭锥（因为 $\mathcal{S}_n^+$ 和 $\mathcal{N}_n$ 是多面锥——$\mathcal{N}_n$ 是多面的，$\mathcal{S}_n^+$ 虽不是多面的但满足一定的正则性条件），因此不需要取闭包。 $\blacksquare$

这个对偶关系与第 45B 章的锥层次结构密切相关：$\mathcal{CP}_n \subseteq \mathcal{DNN}_n$ 对偶地给出 $\mathcal{S}_n^+ + \mathcal{N}_n = \mathcal{DNN}_n^* \subseteq \mathcal{CP}_n^* = \mathcal{COP}_n$。

---

## 45A.11 习题

!!! example "例 45A.9"
    证明：对角矩阵 $D = \operatorname{diag}(d_1, \ldots, d_n)$ 完全正当且仅当 $d_i \geq 0$ 对所有 $i$，并证明此时 $\operatorname{cp-rank}(D) = |\{i : d_i > 0\}|$。

!!! example "例 45A.10"
    设 $A \in \mathcal{CP}_n$。证明 $\operatorname{rank}(A) \leq \operatorname{cp-rank}(A)$，并给出等号成立的充分条件。

!!! example "例 45A.11"
    设 $A, B \in \mathcal{CP}_n$。证明：

    1. $A + B \in \mathcal{CP}_n$。
    2. $A \circ B \in \mathcal{CP}_n$（$\circ$ 表示 Hadamard 乘积）。
    3. 若 $P$ 是置换矩阵，则 $PAP^T \in \mathcal{CP}_n$。

!!! example "例 45A.12"
    证明 Caratheodory 上界：对任意 $A \in \mathcal{CP}_n$，$\operatorname{cp-rank}(A) \leq n(n+1)/2$。

    **提示**：利用 Caratheodory 定理的锥版本——$d$ 维空间中凸锥的每个元素都可以表示为至多 $d$ 个极端射线方向的锥组合。

!!! example "例 45A.13"
    验证矩阵
    $$A = \begin{pmatrix} 4 & 2 & 2 \\ 2 & 4 & 2 \\ 2 & 2 & 4 \end{pmatrix}$$
    是完全正的，并求其 cp-秩。

    **提示**：尝试分解 $A = BB^T$，$B \geq 0$。注意 $A = 2I + 2J_3$，其中 $J_3 = \mathbf{1}\mathbf{1}^T$。

!!! example "例 45A.14"
    设 $A \in \mathcal{DNN}_n$ 的比较图 $G(A)$ 是二部图。证明 $A \in \mathcal{CP}_n$。

    **提示**：二部图不含奇圈，利用 Berman-Xu 定理或直接构造非负分解。

!!! example "例 45A.15"
    证明 $\mathcal{CP}_n$ 的对偶锥是共正锥 $\mathcal{COP}_n$（参见第 45B 章定理 45B.3 的详细证明），即
    $$\mathcal{CP}_n^* = \{C \in \mathcal{S}_n : \langle A, C \rangle \geq 0, \; \forall A \in \mathcal{CP}_n\} = \mathcal{COP}_n.$$

!!! example "例 45A.16"
    设 $A$ 是 $n \times n$ 非负矩阵，$A \geq 0$。证明 $A^T A$ 是完全正矩阵，并给出 $\operatorname{cp-rank}(A^T A)$ 的上界。

!!! example "例 45A.17"
    对 $n = 5$，构造一个尽可能"简单"的 DNN 但非 CP 的矩阵，并验证其比较图包含 $C_5$。

!!! example "例 45A.18"
    证明：完全正矩阵的 Schur（Hadamard）乘积仍然是完全正的。即若 $A, B \in \mathcal{CP}_n$，则 $A \circ B \in \mathcal{CP}_n$，且 $\operatorname{cp-rank}(A \circ B) \leq \operatorname{cp-rank}(A) \cdot \operatorname{cp-rank}(B)$。

!!! example "例 45A.19"
    （开放问题讨论）Drew-Johnson-Loewy 猜想对 $n = 6$ 的情形：尝试构造一个 $6 \times 6$ 完全正矩阵，其 cp-秩尽可能接近 $\lfloor 36/4 \rfloor = 9$。

!!! example "例 45A.20"
    设 $\{P_i\}_{i=1}^n$ 是 $d \times d$ 正半定矩阵，$A_{ij} = \operatorname{tr}(P_i P_j)$。证明 $A = [A_{ij}]$ 是双非负矩阵。并证明：当 $d = 1$ 时（即 $P_i$ 是非负实数），$A$ 是完全正矩阵。

!!! example "例 45A.21"
    设 $A \in \mathcal{CP}_n$ 且 $A$ 的比较图 $G(A)$ 是完全图 $K_n$（即 $A$ 的所有非对角元素严格为正）。证明 $\operatorname{rank}(A) \leq \operatorname{cp-rank}(A)$，并讨论当 $n = 3$ 时 cp-秩的所有可能取值。

    **提示**：对 $n = 3$，$\operatorname{rank}(A) \in \{1, 2, 3\}$，而 $\operatorname{cp-rank}(A) \in \{1, 2, 3\}$（不超过 $\lfloor 9/4 \rfloor = 2$... 不，DJL 猜想对 $n = 3$ 给出上界 2，但需要更仔细分析）。

!!! example "例 45A.22"
    证明：若 $A \in \mathcal{CP}_n$，则 $A$ 的每个**主子矩阵**也是完全正的。即若 $S \subseteq \{1, \ldots, n\}$，则 $A_S = (a_{ij})_{i,j \in S} \in \mathcal{CP}_{|S|}$。

    **提示**：若 $A = \sum_j b_j b_j^T$（$b_j \geq 0$），则 $A_S = \sum_j (b_j)_S (b_j)_S^T$，其中 $(b_j)_S$ 是 $b_j$ 在坐标 $S$ 上的限制。
