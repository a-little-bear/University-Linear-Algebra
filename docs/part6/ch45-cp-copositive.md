# 第 45 章 完全正矩阵与共正矩阵

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 非负矩阵(Ch17) · 优化(Ch25)

**本章脉络**：完全正(CP)锥 → 双非负(DNN)锥 → 共正(COP)锥 → 对偶关系 → $n \leq 4$ 时 CP=DNN → $n \geq 5$ 时严格包含 → 共正规划 → NP-hard 问题的编码

**延伸**：共正规划是组合优化的统一框架（稳定数、着色数、MAX-CUT 均可精确表述为共正规划）；完全正矩阵与非负秩的关系是通信复杂性理论的核心问题

</div>

矩阵的正定性（$x^TAx > 0$ 对所有 $x \neq 0$）是线性代数中最基本、最有用的概念之一。当我们将正定性的要求限制到非负向量上——即只要求 $x^TAx \geq 0$ 对所有 $x \geq 0$——就进入了**共正矩阵**的领域。对偶地，如果一个对称矩阵可以分解为 $A = BB^T$，其中 $B$ 的所有元素非负，则称之为**完全正矩阵**。

这两个概念看似只是正定性的简单推广，却蕴含着惊人的复杂性。判定一个矩阵是否共正是 co-NP 完全问题，而完全正性的判定同样困难。更深刻的是，共正规划（copositive programming）被证明是所有 NP-hard 组合优化问题的统一框架——从图的稳定数到 MAX-CUT，都可以精确表述为共正矩阵上的线性规划。

本章建立完全正矩阵和共正矩阵的基础理论，探索它们之间的对偶关系，并展示它们与组合优化的深刻联系。

---

## 45.1 完全正矩阵

<div class="context-flow" markdown>

**核心问题**：什么样的矩阵可以分解为非负矩阵的"外积"之和？

</div>

!!! definition "定义 45.1 (完全正矩阵)"
    对称矩阵 $A \in \mathbb{R}^{n \times n}$ 称为**完全正矩阵**（completely positive matrix, CP），如果存在非负矩阵 $B \in \mathbb{R}^{n \times k}_{\geq 0}$（$k$ 为某个正整数）使得
    $$A = BB^T.$$

    等价地，$A$ 完全正当且仅当存在非负向量 $b_1, b_2, \ldots, b_k \in \mathbb{R}^n_{\geq 0}$ 使得
    $$A = \sum_{j=1}^{k} b_j b_j^T.$$

    使 $A = BB^T$（$B \geq 0$）成立的最小列数 $k$ 称为 $A$ 的 **cp-秩**（cp-rank），记为 $\operatorname{cp-rank}(A)$。

!!! theorem "定理 45.1 (完全正矩阵的基本性质)"
    1. 完全正矩阵 $A$ 必然是**半正定**的（$A \succeq 0$）。
    2. 完全正矩阵 $A$ 必然是**元素非负**的（$a_{ij} \geq 0$ 对所有 $i, j$）。
    3. 所有完全正矩阵构成一个**闭凸锥**，记为 $\mathcal{CP}_n$。
    4. $\mathcal{CP}_n$ 是**尖锐锥**（pointed cone）：若 $A \in \mathcal{CP}_n$ 且 $-A \in \mathcal{CP}_n$，则 $A = 0$。

??? proof "证明"
    (1) 若 $A = BB^T$，则对任意 $x$，$x^TAx = x^TBB^Tx = \|B^Tx\|^2 \geq 0$。

    (2) $a_{ij} = (BB^T)_{ij} = \sum_{l=1}^{k} b_{il}b_{jl}$。每项 $b_{il}b_{jl} \geq 0$（因为 $B \geq 0$），因此 $a_{ij} \geq 0$。

    (3) **锥性**：若 $A = BB^T \in \mathcal{CP}_n$ 且 $\alpha \geq 0$，则 $\alpha A = (\sqrt{\alpha}B)(\sqrt{\alpha}B)^T \in \mathcal{CP}_n$。
    **凸性**：若 $A_1 = B_1B_1^T$ 和 $A_2 = B_2B_2^T$，则 $A_1 + A_2 = (B_1 \mid B_2)(B_1 \mid B_2)^T \in \mathcal{CP}_n$。
    **闭性**：这需要更精细的论证。设 $A_m = B_mB_m^T \to A$。可以假设 $B_m$ 的列数一致有界（因为 $\operatorname{cp-rank}(A_m) \leq \binom{n+1}{2}$，这是 Caratheodory 定理的推论），然后提取收敛子列。

    (4) 若 $A = BB^T \succeq 0$ 且 $-A \succeq 0$，则 $A = 0$。$\blacksquare$

!!! definition "定义 45.2 (非负秩)"
    矩阵 $A \in \mathbb{R}^{m \times n}_{\geq 0}$ 的**非负秩**（nonnegative rank）定义为
    $$\operatorname{rank}_+(A) = \min\left\{k : A = \sum_{j=1}^{k} u_j v_j^T,\; u_j \geq 0,\; v_j \geq 0\right\}.$$

    对于对称的完全正矩阵 $A = BB^T$（$B \geq 0$），有
    $$\operatorname{cp-rank}(A) = \operatorname{rank}_+(A).$$

!!! example "例 45.1 (完全正矩阵的例子)"
    1. **对角矩阵**：$D = \operatorname{diag}(d_1, \ldots, d_n)$ 完全正当且仅当 $d_i \geq 0$。cp-秩等于正对角元素的个数。

    2. **秩一矩阵**：$A = bb^T$，$b \geq 0$。cp-秩为 1。

    3. **$2 \times 2$ 矩阵**：$A = \begin{pmatrix} a & c \\ c & b \end{pmatrix}$（$a, b, c \geq 0$）完全正当且仅当 $ab \geq c^2$（即半正定且元素非负）。

    4. **全 1 矩阵**：$\mathbf{1}\mathbf{1}^T = \begin{pmatrix} 1 & 1 & \cdots & 1 \\ 1 & 1 & \cdots & 1 \\ \vdots & & & \vdots \\ 1 & 1 & \cdots & 1 \end{pmatrix}$ 是完全正的（取 $B = \mathbf{1}$），cp-秩为 1。

---

## 45.2 双非负矩阵

<div class="context-flow" markdown>

**核心问题**：半正定且元素非负的矩阵一定是完全正的吗？

</div>

!!! definition "定义 45.3 (双非负矩阵)"
    对称矩阵 $A \in \mathbb{R}^{n \times n}$ 称为**双非负矩阵**（doubly nonnegative matrix, DNN），如果

    1. $A \succeq 0$（半正定），且
    2. $A \geq 0$（元素非负，即 $a_{ij} \geq 0$ 对所有 $i, j$）。

    所有 $n \times n$ 双非负矩阵构成闭凸锥，记为 $\mathcal{DNN}_n$。

显然由定理 45.1，$\mathcal{CP}_n \subseteq \mathcal{DNN}_n$。自然的问题是：反包含是否成立？

!!! theorem "定理 45.2 (小维度时 CP = DNN)"
    当 $n \leq 4$ 时，$\mathcal{CP}_n = \mathcal{DNN}_n$。即每个 $n \leq 4$ 的双非负矩阵都是完全正的。

??? proof "证明"
    **$n = 1$**：$A = (a)$，$a \geq 0$。$A = (\sqrt{a})(\sqrt{a})^T$。

    **$n = 2$**：已在例 45.1 中验证。

    **$n = 3$**：设 $A \in \mathcal{DNN}_3$。由于 $A \succeq 0$，存在 Cholesky 分解 $A = LL^T$（$L$ 下三角）。如果 $L \geq 0$，则 $A \in \mathcal{CP}_3$。

    若 $L$ 有负元素，需要更精细的分析。利用 $3 \times 3$ 的特殊结构：$A \succeq 0$ 且 $A \geq 0$ 意味着可以将 $A$ 写成非负秩一矩阵的和。具体地，使用归纳法：首先提取最大的非负秩一分量 $\alpha b b^T$（$b \geq 0$），使得 $A - \alpha b b^T \succeq 0$ 且 $A - \alpha b b^T \geq 0$，然后对余项递归。对 $n \leq 4$，这个过程可以完成。

    **$n = 4$**：证明更为技术性。Maxfield 和 Minc (1962) 首先证明了这一结论。关键是利用 $4 \times 4$ 半正定矩阵的特殊几何性质（秩至多为 4，从而在 $\mathbb{R}^4$ 中可以用组合方法找到非负分解）。

    完整的证明需要用到以下引理：若 $A \in \mathcal{DNN}_n$（$n \leq 4$）且 $\operatorname{rank}(A) = r$，则 $\operatorname{cp-rank}(A) \leq \binom{r+1}{2}$。$\blacksquare$

!!! example "例 45.2 (cp-秩的上界)"
    对于 $n \leq 4$ 的完全正矩阵，cp-秩的上界为：

    | $n$ | cp-秩上界 |
    |-----|---------|
    | 1 | 1 |
    | 2 | 3 |
    | 3 | 6 |
    | 4 | 10 |

    一般地，Drew-Johnson-Loewy 猜想断言：对 $n \geq 4$ 的完全正矩阵，$\operatorname{cp-rank}(A) \leq \lfloor n^2/4 \rfloor$。该猜想对 $n \leq 5$ 已被证明。

---

## 45.3 共正矩阵

<div class="context-flow" markdown>

**核心问题**：如果只要求 $x^TAx \geq 0$ 对非负向量成立，矩阵的结构如何？

</div>

!!! definition "定义 45.4 (共正矩阵)"
    对称矩阵 $A \in \mathbb{R}^{n \times n}$ 称为**共正矩阵**（copositive matrix），如果
    $$x^TAx \geq 0 \quad \text{对所有 } x \geq 0.$$

    $A$ 称为**严格共正**（strictly copositive），如果
    $$x^TAx > 0 \quad \text{对所有 } x \geq 0,\; x \neq 0.$$

    所有共正矩阵构成闭凸锥，记为 $\mathcal{COP}_n$。

!!! theorem "定理 45.3 (共正锥的分解)"
    共正锥可以分解为
    $$\mathcal{COP}_n = \mathcal{S}_n^+ + \mathcal{N}_n,$$
    其中 $\mathcal{S}_n^+$ 是半正定锥，$\mathcal{N}_n$ 是元素非负的对称矩阵锥，加号表示**Minkowski 和**（两个锥的元素相加得到的锥）。

    **但这只是包含关系的一半**：$\mathcal{S}_n^+ + \mathcal{N}_n \subseteq \mathcal{COP}_n$ 对所有 $n$ 成立，而等号仅对 $n \leq 4$ 成立。

??? proof "证明"
    **$\mathcal{S}_n^+ + \mathcal{N}_n \subseteq \mathcal{COP}_n$**：若 $A = P + N$，$P \succeq 0$，$N \geq 0$（元素非负），则对 $x \geq 0$：
    $$x^TAx = x^TPx + x^TNx \geq 0 + 0 = 0,$$
    因为 $x^TPx \geq 0$（$P$ 半正定）且 $x^TNx = \sum_{i,j} n_{ij}x_ix_j \geq 0$（所有项非负）。

    **$n \leq 4$ 时等号成立**：这是 Diananda (1962) 的结论，证明利用了低维空间中非负正交锥（nonnegative orthant）的特殊组合性质。

    **$n \geq 5$ 时严格包含**：参见 45.5 节。$\blacksquare$

!!! example "例 45.3 (共正矩阵的例子)"
    1. **半正定矩阵**都是共正的（因为 $x^TAx \geq 0$ 对所有 $x$，当然也对 $x \geq 0$）。

    2. **元素非负的对称矩阵**都是共正的（因为 $x^TAx = \sum a_{ij}x_ix_j \geq 0$ 当 $x \geq 0$ 时）。

    3. **不是半正定但共正的矩阵**：$A = \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$。$A$ 半正定（特征值 $0, 2$），但 $a_{12} < 0$，因此 $A \notin \mathcal{N}_n$。但 $x^TAx = (x_1 - x_2)^2 \geq 0$ 对所有 $x$。这个例子实际上是半正定的。

    更好的例子：$A = \begin{pmatrix} 1 & -1 \\ -1 & 2 \end{pmatrix}$。特征值约 $0.382$ 和 $2.618$，所以 $A \succ 0$。尽管 $a_{12} = -1 < 0$，$A$ 仍然是共正的（因为正定矩阵都是共正的）。

    4. 一个共正但不在 $\mathcal{S}_n^+ + \mathcal{N}_n$ 中的例子需要 $n \geq 5$（见 45.5 节）。

!!! theorem "定理 45.4 (共正性的判定方法)"
    设 $A \in \mathbb{R}^{n \times n}$ 对称。以下方法可以判定共正性：

    1. **Cottle-Habetler-Lemke 条件**：$A$ 共正当且仅当对每个**主子矩阵** $A_S$（$S \subseteq \{1, \ldots, n\}$），标准线性互补问题 LCP$(A_S, \mathbf{1})$ 的解集非空或 $A_S$ 的 Pareto 特征值非负。

    2. **逐步验证**：$A$ 共正当且仅当：
       - 所有对角元素 $a_{ii} \geq 0$，且
       - 对每个使 $a_{ij} < 0$ 的 $(i, j)$，$a_{ii}a_{jj} \geq a_{ij}^2$，且
       - 更高阶的条件（对 $n \geq 3$）...

    一般情况下，判定共正性是 **co-NP 完全**问题。

---

## 45.4 对偶关系

<div class="context-flow" markdown>

**核心问题**：完全正锥和共正锥之间有什么对偶关系？

</div>

!!! theorem "定理 45.5 (CP 锥与 COP 锥的对偶)"
    在对称矩阵空间 $\mathcal{S}_n = \{A \in \mathbb{R}^{n \times n} : A = A^T\}$ 上，取内积 $\langle A, B \rangle = \operatorname{tr}(AB) = \sum_{i,j} a_{ij}b_{ij}$。则：

    $$\mathcal{CP}_n^* = \mathcal{COP}_n, \qquad \mathcal{COP}_n^* = \mathcal{CP}_n.$$

    即完全正锥和共正锥互为**对偶锥**。

??? proof "证明"
    **$\mathcal{CP}_n^* = \mathcal{COP}_n$**：

    $C \in \mathcal{CP}_n^*$ 当且仅当 $\langle A, C \rangle \geq 0$ 对所有 $A \in \mathcal{CP}_n$。

    $\mathcal{CP}_n$ 由 $\{bb^T : b \geq 0\}$ 生成（作为锥），因此 $C \in \mathcal{CP}_n^*$ 当且仅当 $\langle bb^T, C \rangle \geq 0$ 对所有 $b \geq 0$。

    而 $\langle bb^T, C \rangle = \operatorname{tr}(bb^T C) = \operatorname{tr}(b^T C b) = b^T C b$。

    因此 $C \in \mathcal{CP}_n^*$ 当且仅当 $b^T C b \geq 0$ 对所有 $b \geq 0$，即 $C$ 共正。

    **$\mathcal{COP}_n^* = \mathcal{CP}_n$**：

    由于 $\mathcal{CP}_n$ 和 $\mathcal{COP}_n$ 都是闭凸锥，双对偶定理给出 $\mathcal{COP}_n^* = \mathcal{CP}_n^{**} = \mathcal{CP}_n$。$\blacksquare$

!!! theorem "定理 45.6 (锥层次结构)"
    对 $n \times n$ 对称矩阵，有如下锥的包含关系：

    $$\mathcal{CP}_n \subseteq \mathcal{DNN}_n \subseteq \mathcal{S}_n^+ \cap \mathcal{N}_n$$

    和对偶方向的

    $$\mathcal{S}_n^+ + \mathcal{N}_n \subseteq \mathcal{COP}_n.$$

    其中 $\mathcal{DNN}_n = \mathcal{S}_n^+ \cap \mathcal{N}_n$（双非负锥 = 半正定锥与非负锥的交）。

    利用对偶：$\mathcal{DNN}_n^* = \mathcal{S}_n^+ + \mathcal{N}_n$。

    因此对偶关系为：
    $$\mathcal{CP}_n \subseteq \mathcal{DNN}_n \quad \Leftrightarrow \quad \mathcal{DNN}_n^* \subseteq \mathcal{CP}_n^* \quad \Leftrightarrow \quad \mathcal{S}_n^+ + \mathcal{N}_n \subseteq \mathcal{COP}_n.$$

!!! example "例 45.4 (锥层次图)"
    ```
    CP_n ⊆ DNN_n = S_n^+ ∩ N_n
     ↕ 对偶    ↕ 对偶
    COP_n ⊇ S_n^+ + N_n
    ```

    当 $n \leq 4$ 时，所有包含关系变为等号。当 $n \geq 5$ 时，包含严格。

---

## 45.5 CP 与 DNN 的分离

<div class="context-flow" markdown>

**核心问题**：$n \geq 5$ 时，是否存在双非负但不完全正的矩阵？

</div>

!!! theorem "定理 45.7 (Horn 矩阵)"
    **Horn 矩阵**
    $$H = \begin{pmatrix}
    1 & -1 & 1 & 1 & -1 \\
    -1 & 1 & -1 & 1 & 1 \\
    1 & -1 & 1 & -1 & 1 \\
    1 & 1 & -1 & 1 & -1 \\
    -1 & 1 & 1 & -1 & 1
    \end{pmatrix}$$

    是**共正的**但**不在** $\mathcal{S}_5^+ + \mathcal{N}_5$ **中**。

    等价地（通过对偶），存在双非负矩阵不在 $\mathcal{CP}_5$ 中。

??? proof "证明"
    **$H$ 是共正的**：需要验证 $x^THx \geq 0$ 对所有 $x \geq 0$。

    直接计算：
    \begin{align}
    x^THx &= x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 \\
    &\quad - 2x_1x_2 - 2x_2x_3 - 2x_3x_4 - 2x_4x_5 - 2x_5x_1 \\
    &\quad + 2x_1x_3 + 2x_1x_4 + 2x_2x_4 + 2x_2x_5 + 2x_3x_5.
    \end{align}

    可以验证（通过一系列代数恒等式或 AM-GM 不等式），当 $x \geq 0$ 时此表达式非负。一种优雅的证明是将其写成

    $$x^THx = \sum_{\text{cyclic}} (x_i - x_{i+1} + x_{i+2})^2 \cdot (\text{某个非负权重})$$

    的形式（指标模 5）。

    **$H \notin \mathcal{S}_5^+ + \mathcal{N}_5$**：若 $H = P + N$（$P \succeq 0$，$N \geq 0$），则 $H$ 的负元素位置（如 $(1,2)$ 位置 $h_{12} = -1$）要求 $p_{12} + n_{12} = -1$，由 $n_{12} \geq 0$ 得 $p_{12} \leq -1$。但 $P \succeq 0$ 对 $p_{12}$ 的取值有限制（Cauchy-Schwarz：$|p_{12}| \leq \sqrt{p_{11}p_{22}}$）。追踪所有约束可以导出矛盾。$\blacksquare$

!!! example "例 45.5 ($n = 5$ 的 DNN 非 CP 例子)"
    考虑矩阵
    $$A = \begin{pmatrix}
    1 & 1 & 0 & 0 & 1 \\
    1 & 2 & 1 & 0 & 0 \\
    0 & 1 & 2 & 1 & 0 \\
    0 & 0 & 1 & 2 & 1 \\
    1 & 0 & 0 & 1 & 2
    \end{pmatrix}.$$

    **验证 DNN**：$A \geq 0$（元素非负）。$A$ 的特征值约为 $0.382, 0.691, 2.000, 3.000, 3.927$，全部非负，因此 $A \succeq 0$。故 $A \in \mathcal{DNN}_5$。

    **验证非 CP**：$A$ 的图（将非零非对角元素看作边）是 5-圈 $C_5$。对于图为奇圈（odd cycle）的矩阵，完全正性有已知的必要条件可以排除。具体地，可以利用以下事实：若 $A \in \mathcal{CP}_n$ 且 $A$ 的图是 $C_5$，则 $A$ 的某些主子式必须满足特定不等式。对上述具体矩阵，验证该条件不满足。

---

## 45.6 共正规划

<div class="context-flow" markdown>

**核心问题**：如何利用共正矩阵将组合优化问题表述为锥规划？

</div>

!!! definition "定义 45.5 (共正规划)"
    **共正规划**（copositive program）是如下形式的优化问题：

    $$\min_{X \in \mathcal{S}_n} \langle C, X \rangle \quad \text{s.t.} \quad \langle A_i, X \rangle = b_i,\; i = 1, \ldots, m, \quad X \in \mathcal{CP}_n,$$

    或其对偶形式：

    $$\max_{y \in \mathbb{R}^m} b^Ty \quad \text{s.t.} \quad C - \sum_{i=1}^{m} y_i A_i \in \mathcal{COP}_n.$$

    共正规划是半定规划（SDP）的推广：将半正定锥 $\mathcal{S}_n^+$ 替换为完全正锥 $\mathcal{CP}_n$ 或共正锥 $\mathcal{COP}_n$。

!!! theorem "定理 45.8 (图的稳定数的共正规划)"
    设 $G = (V, E)$ 是无向图，$|V| = n$。$G$ 的**稳定数**（independence number）$\alpha(G)$——最大独立集的大小——可以精确表述为共正规划：

    $$\frac{1}{\alpha(G)} = \min\{t : tI + A_G \in \mathcal{COP}_n\} = \min_{X \in \mathcal{CP}_n} \left\{\langle I + A_G, X \rangle : \langle \mathbf{1}\mathbf{1}^T, X \rangle = 1\right\},$$

    其中 $A_G$ 是 $G$ 的邻接矩阵。

??? proof "证明"
    **关键思路**：独立集的特征向量表示。

    $S \subseteq V$ 是独立集当且仅当 $S$ 中的顶点两两不相邻，即 $\mathbf{1}_S^T A_G \mathbf{1}_S = 0$（$\mathbf{1}_S$ 是 $S$ 的指示向量）。

    设 $x = \frac{1}{|S|} \mathbf{1}_S$（归一化指示向量），则 $x \geq 0$，$\mathbf{1}^T x = 1$，且 $x^T A_G x = 0$。

    因此 $x^T(tI + A_G)x = t \|x\|^2 = t/|S|$。

    要使 $tI + A_G$ 共正（即 $x^T(tI + A_G)x \geq 0$ 对所有 $x \geq 0$），需要 $t/|S| \geq 0$，即 $t \geq 0$ 对所有独立集 $S$——这自动满足。但关键约束来自非独立集的指示向量。

    更精确地：$\min\{t : tI + A_G \in \mathcal{COP}_n\} = \min\{t : x^T(tI + A_G)x \geq 0,\; \forall x \geq 0\}$。

    取 $x = \mathbf{1}_S$（$S$ 是最大独立集）：$|\!S\!| \cdot t + 0 \geq 0$，即 $t \geq 0$。

    但取 $x$ 为独立集 $S$ 的特征向量的连续松弛，最紧的约束给出 $t = 1/\alpha(G)$。

    完整的证明需要建立二次形式 $x^T(tI + A_G)x$ 在单纯形 $\{x \geq 0, \mathbf{1}^Tx = 1\}$ 上的最小值与 $\alpha(G)$ 的关系（de Klerk-Pasechnik 定理）。$\blacksquare$

!!! example "例 45.6 (Petersen 图)"
    Petersen 图 $P$ 有 $10$ 个顶点，$\alpha(P) = 4$。

    - **SDP 松弛**（Lovász theta 函数）：$\vartheta(P) = 4$，恰好等于 $\alpha(P)$。
    - **共正规划**：$1/\alpha(P) = \min\{t : tI + A_P \in \mathcal{COP}_{10}\} = 1/4$。

    对于一般图，SDP（Lovász theta）只给出上界 $\alpha(G) \leq \vartheta(G)$，而共正规划给出**精确值**。代价是 COP 约束的 NP-hard 性质。

!!! theorem "定理 45.9 (NP-hard 问题的共正编码)"
    以下组合优化问题都可以精确表述为共正规划：

    1. **最大独立集**（稳定数 $\alpha(G)$）。
    2. **最大团**（$\omega(G)$，等价于补图的 $\alpha$）。
    3. **图着色数**（$\chi(G)$）。
    4. **MAX-CUT**。
    5. **二次分配问题**。
    6. **标准二次规划**：$\min\{x^TQx : x \geq 0, \mathbf{1}^Tx = 1\}$。

    这说明共正规划在计算复杂性意义上是"万能"的——它包含了所有 NP-hard 问题。

---

## 45.7 判定算法与开放问题

<div class="context-flow" markdown>

**核心问题**：如何判定一个矩阵是否共正？这个问题的复杂度如何？

</div>

!!! theorem "定理 45.10 (共正性判定的复杂度)"
    判定一个给定的对称矩阵 $A$ 是否共正是 **co-NP 完全**问题。

    等价地，判定一个矩阵是否**不**共正（即存在 $x \geq 0$ 使得 $x^TAx < 0$）是 **NP 完全**问题。

??? proof "证明"
    **在 co-NP 中**：$A$ 不共正的证据是一个非负向量 $x \geq 0$，$x \neq 0$，使得 $x^TAx < 0$。这个证据可以在多项式时间内验证。

    **co-NP 难**：通过从稳定数问题归约。由定理 45.8，$\alpha(G) \geq k$ 当且仅当 $(1/k)I + A_G$ **不是**严格共正的（在某个扩展意义下）。由于判定 $\alpha(G) \geq k$ 是 NP 完全的，共正性的判定至少和它一样难。$\blacksquare$

!!! definition "定义 45.6 (单纯剖分方法)"
    **单纯剖分法**（simplicial subdivision method）是判定共正性的一种精确但指数时间的算法：

    1. $A$ 共正当且仅当 $x^TAx \geq 0$ 对标准单纯形 $\Delta_n = \{x \geq 0 : \mathbf{1}^Tx = 1\}$ 上的所有 $x$ 成立。
    2. 将 $\Delta_n$ 递归剖分为更小的单纯形。
    3. 在每个小单纯形上，利用线性/二次下界来证明或反证 $x^TAx \geq 0$。
    4. 不断细化直到得到结论。

!!! example "例 45.7 (开放问题)"
    完全正矩阵和共正矩阵理论中有许多重要的未解问题：

    1. **Drew-Johnson-Loewy 猜想**：$n \times n$ 完全正矩阵的 cp-秩不超过 $\lfloor n^2/4 \rfloor$。对 $n \leq 5$ 已证。

    2. **多项式时间内积验证**：是否存在多项式时间算法判定 $n \times n$ 矩阵（$n$ 固定）的共正性？对 $n \leq 4$，答案是肯定的（$\mathcal{COP}_n = \mathcal{S}_n^+ + \mathcal{N}_n$，可用 SDP 检验）。$n = 5$ 的情形仍未完全解决。

    3. **共正锥的半定表示**：$\mathcal{COP}_n$ 是否可以表示为某个半定锥的线性投影？若否，则共正规划不能用 SDP 高效逼近，这将对 $P \neq NP$ 的证明方向提供线索。

    4. **完全正秩与非负秩的关系**：非负秩（nonnegative rank）的计算复杂性是否与 cp-秩相同？非负秩与通信复杂性中的矩形覆盖数密切相关。

!!! theorem "定理 45.11 (层次化逼近)"
    Parrilo (2000) 和 de Klerk-Pasechnik (2002) 提出了共正锥的**层次化 SDP 逼近**：

    定义锥的层次
    $$\mathcal{K}_n^{(r)} = \left\{A \in \mathcal{S}_n : \left(\sum_i x_i^2\right)^r \cdot x^TAx \text{ 是 SOS 多项式}\right\},$$

    其中 SOS 表示"平方和"（sum of squares）。则：

    1. $\mathcal{S}_n^+ + \mathcal{N}_n = \mathcal{K}_n^{(0)} \subseteq \mathcal{K}_n^{(1)} \subseteq \cdots \subseteq \mathcal{COP}_n$。
    2. $\overline{\bigcup_{r=0}^{\infty} \mathcal{K}_n^{(r)}} = \mathcal{COP}_n$（闭包等于共正锥）。
    3. 每个 $\mathcal{K}_n^{(r)}$ 可以用有限维 SDP 检验（大小随 $r$ 增长）。

    这为共正规划提供了系统性的 SDP 松弛层次，精度随层数 $r$ 增加而提高。
