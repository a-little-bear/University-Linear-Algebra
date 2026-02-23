# 第 45B 章 余正矩阵与余正规划

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 凸优化与锥规划(Ch25) · 完全正矩阵(Ch45A)

**本章脉络**：余正(COP)矩阵定义 → COP 锥分解与 $n \leq 4$ 时的等式 → 余正性判定（Cottle-Habetler-Lemke）→ CP 与 COP 的对偶关系（完整证明）→ 锥层次图 → COP 锥的极端射线与面结构 → Bomze 边界刻画 → 独立数的余正规划（de Klerk-Pasechnik）→ NP-hard 编码：稳定数、团数、着色数、MAX-CUT → co-NP 完全性 → Parrilo 的 SDP 层次逼近 → 线性互补问题 → Sturm 的 SDP 可表示性

**延伸**：余正规划是组合优化的统一框架——所有 NP-hard 组合优化问题均可精确表述为余正规划；Parrilo 的 SOS 层次提供了系统性的 SDP 逼近方法

</div>

如果将正定性的要求从"对所有向量 $x$"限制到"对所有非负向量 $x \geq 0$"，就进入了**余正矩阵**（copositive matrix）的领域。对称矩阵 $A$ 余正意味着 $x^T A x \geq 0$ 对所有 $x \geq 0$——这看似只是正定性的一个弱化，却蕴含着惊人的计算复杂性。

余正矩阵构成的闭凸锥 $\mathcal{COP}_n$ 与完全正锥 $\mathcal{CP}_n$ 互为对偶锥，这一对偶关系是本章的核心。更深刻的是，**余正规划**（copositive programming）——即以余正锥或完全正锥为约束的线性规划——被证明是所有 NP-hard 组合优化问题的统一框架。从图的独立数到 MAX-CUT，从图着色到最大团，所有这些经典难题都可以精确表述为余正规划。

本章建立余正矩阵的基本理论，深入探索余正锥的几何结构，展示组合优化问题的余正规划公式，并介绍 Parrilo 的 SDP 层次逼近方法。

---

## 45B.1 余正矩阵的定义

<div class="context-flow" markdown>

**核心问题**：如果只要求 $x^T A x \geq 0$ 对非负向量成立，矩阵需要满足什么条件？

</div>

!!! definition "定义 45B.1 (余正矩阵)"
    对称矩阵 $A \in \mathbb{R}^{n \times n}$ 称为**余正矩阵**（copositive matrix），如果
    $$x^T A x \geq 0 \quad \text{对所有 } x \geq 0.$$

    $A$ 称为**严格余正**（strictly copositive），如果
    $$x^T A x > 0 \quad \text{对所有 } x \geq 0, \; x \neq 0.$$

    所有 $n \times n$ 余正矩阵构成闭凸锥，记为 $\mathcal{COP}_n$。

!!! example "例 45B.1 (余正矩阵的例子)"
    1. **半正定矩阵**都是余正的：$A \succeq 0 \Rightarrow x^T A x \geq 0$ 对所有 $x$（不仅是 $x \geq 0$）。因此 $\mathcal{S}_n^+ \subseteq \mathcal{COP}_n$。

    2. **元素非负的对称矩阵**都是余正的：若 $a_{ij} \geq 0$ 对所有 $i, j$，则当 $x \geq 0$ 时 $x^T A x = \sum_{i,j} a_{ij} x_i x_j \geq 0$。因此 $\mathcal{N}_n \subseteq \mathcal{COP}_n$。

    3. **含负元素但余正的矩阵**：$A = \begin{pmatrix} 1 & -1 \\ -1 & 2 \end{pmatrix}$。$A$ 正定（特征值约 $0.382$ 和 $2.618$），因此余正。尽管 $a_{12} = -1 < 0$，$A$ 仍然余正。

    4. **不是半正定也不元素非负，但余正的矩阵**：这样的例子需要 $n \geq 5$（见 45B.3 节的 Horn 矩阵）。

---

## 45B.2 COP 锥的分解

<div class="context-flow" markdown>

**核心问题**：余正锥能否分解为更简单的锥的 Minkowski 和？

</div>

!!! theorem "定理 45B.1 (余正锥的分解)"
    1. 对所有 $n$，有包含关系
    $$\mathcal{S}_n^+ + \mathcal{N}_n \subseteq \mathcal{COP}_n,$$
    其中 $\mathcal{S}_n^+$ 是半正定锥，$\mathcal{N}_n$ 是元素非负的对称矩阵锥，加号表示 **Minkowski 和**。

    2. 当 $n \leq 4$ 时，等号成立：
    $$\mathcal{COP}_n = \mathcal{S}_n^+ + \mathcal{N}_n \quad (n \leq 4).$$

    3. 当 $n \geq 5$ 时，包含严格：
    $$\mathcal{S}_n^+ + \mathcal{N}_n \subsetneq \mathcal{COP}_n \quad (n \geq 5).$$

??? proof "证明"
    **(1)** 设 $A = P + N$，$P \succeq 0$，$N \geq 0$（元素非负）。对任意 $x \geq 0$：
    $$x^T A x = x^T P x + x^T N x \geq 0 + 0 = 0,$$
    因为 $x^T P x \geq 0$（$P$ 半正定）且 $x^T N x = \sum_{i,j} n_{ij} x_i x_j \geq 0$（$n_{ij} \geq 0$，$x_i x_j \geq 0$）。

    **(2)** $n \leq 4$ 时等号成立是 Diananda (1962) 的结果。证明思路如下：

    设 $A \in \mathcal{COP}_n$（$n \leq 4$）。需要证明 $A = P + N$（$P \succeq 0$，$N \geq 0$）。

    首先，$A$ 余正意味着所有对角元素 $a_{ii} \geq 0$（取 $x = e_i$）。

    对 $n = 2$：若 $a_{12} \geq 0$，取 $P = 0$，$N = A$。若 $a_{12} < 0$，由余正性 $a_{11} a_{22} \geq a_{12}^2$（否则取 $x = (\sqrt{a_{22}}, \sqrt{a_{11}})^T$ 可使 $x^T A x < 0$）。因此 $A$ 半正定，取 $P = A$，$N = 0$。

    对 $n = 3, 4$，利用低维空间中非负正交锥的特殊组合性质（$\mathbb{R}^n_{\geq 0}$ 在 $n \leq 4$ 时的特殊结构），逐步分析 $A$ 的负元素位置，将其分解为半正定部分和非负部分。

    **(3)** Horn 矩阵（定理 45A.8）提供了 $n = 5$ 的反例。 $\blacksquare$

---

## 45B.3 余正性的判定方法

<div class="context-flow" markdown>

**核心问题**：如何判定一个给定的对称矩阵是否余正？

</div>

!!! theorem "定理 45B.2 (Cottle-Habetler-Lemke 条件)"
    设 $A \in \mathbb{R}^{n \times n}$ 对称。以下条件是判定余正性的充分或必要条件：

    1. **必要条件**：
        - 所有对角元素 $a_{ii} \geq 0$。
        - 若 $a_{ij} < 0$，则 $a_{ii} > 0$ 且 $a_{jj} > 0$，且 $a_{ii} a_{jj} > a_{ij}^2$。

    2. **充分条件**（Cottle-Habetler-Lemke）：
        $A$ 余正如果对每个非空子集 $S \subseteq \{1, \ldots, n\}$，主子矩阵 $A_S$ 满足以下条件之一：
        - $A_S$ 有一个正对角元素，或
        - $A_S$ 半正定，或
        - $A_S$ 的所有元素非负。

    3. **精确刻画**：$A$ 余正当且仅当 $\min\{x^T A x : x \geq 0, \mathbf{1}^T x = 1\} \geq 0$，即标准单纯形 $\Delta_n$ 上的二次形式 $x^T A x$ 非负。

??? proof "证明"
    **(1)** 取 $x = e_i$：$e_i^T A e_i = a_{ii} \geq 0$。取 $x = \alpha e_i + \beta e_j$（$\alpha, \beta \geq 0$）：$x^T A x = \alpha^2 a_{ii} + 2\alpha\beta a_{ij} + \beta^2 a_{jj} \geq 0$。若 $a_{ij} < 0$，则需要 $a_{ii} a_{jj} \geq a_{ij}^2$（否则取 $\alpha = \sqrt{a_{jj}|a_{ij}|}$，$\beta = \sqrt{a_{ii}|a_{ij}|}$ 可使其为负）。

    **(3)** 由于 $x^T A x$ 对 $x$ 是齐次的（2 次），$x^T A x \geq 0$ 对所有 $x \geq 0$ 等价于 $x^T A x \geq 0$ 对单纯形 $\{x \geq 0 : \mathbf{1}^T x = 1\}$ 上的所有 $x$。 $\blacksquare$

!!! definition "定义 45B.2 (单纯剖分方法)"
    **单纯剖分法**（simplicial subdivision method）是判定余正性的一种精确但指数时间的算法：

    1. $A$ 余正当且仅当 $x^T A x \geq 0$ 对标准单纯形 $\Delta_n = \{x \geq 0 : \mathbf{1}^T x = 1\}$ 上的所有 $x$ 成立。
    2. 将 $\Delta_n$ 递归剖分为更小的单纯形。
    3. 在每个小单纯形上，利用线性/二次下界来证明或反证 $x^T A x \geq 0$。
    4. 不断细化直到得到结论。

---

## 45B.4 CP 与 COP 的对偶关系

<div class="context-flow" markdown>

**核心问题**：完全正锥和余正锥之间有什么对偶关系？

</div>

!!! theorem "定理 45B.3 (CP 锥与 COP 锥的对偶)"
    在对称矩阵空间 $\mathcal{S}_n = \{A \in \mathbb{R}^{n \times n} : A = A^T\}$ 上，取内积 $\langle A, B \rangle = \operatorname{tr}(AB) = \sum_{i,j} a_{ij} b_{ij}$。则完全正锥和余正锥互为**对偶锥**：

    $$\mathcal{CP}_n^* = \mathcal{COP}_n, \qquad \mathcal{COP}_n^* = \mathcal{CP}_n.$$

??? proof "证明"
    **第一个等式**：$\mathcal{CP}_n^* = \mathcal{COP}_n$。

    回忆对偶锥的定义：$\mathcal{CP}_n^* = \{C \in \mathcal{S}_n : \langle A, C \rangle \geq 0, \; \forall A \in \mathcal{CP}_n\}$。

    **$\mathcal{COP}_n \subseteq \mathcal{CP}_n^*$**：设 $C \in \mathcal{COP}_n$，$A \in \mathcal{CP}_n$。则 $A = \sum_{j=1}^k b_j b_j^T$（$b_j \geq 0$），因此
    $$\langle A, C \rangle = \operatorname{tr}(AC) = \sum_{j=1}^k \operatorname{tr}(b_j b_j^T C) = \sum_{j=1}^k b_j^T C b_j \geq 0,$$
    其中最后一步由 $C$ 余正（$b_j^T C b_j \geq 0$，因为 $b_j \geq 0$）。

    **$\mathcal{CP}_n^* \subseteq \mathcal{COP}_n$**：设 $C \in \mathcal{CP}_n^*$。$\mathcal{CP}_n$ 由 $\{bb^T : b \geq 0\}$ 生成（作为锥），因此 $C \in \mathcal{CP}_n^*$ 意味着
    $$\langle bb^T, C \rangle \geq 0 \quad \text{对所有 } b \geq 0.$$

    计算：
    $$\langle bb^T, C \rangle = \operatorname{tr}(bb^T C) = \operatorname{tr}(b^T C b) = b^T C b.$$

    因此 $b^T C b \geq 0$ 对所有 $b \geq 0$，即 $C$ 余正：$C \in \mathcal{COP}_n$。

    综合两个方向：$\mathcal{CP}_n^* = \mathcal{COP}_n$。

    **第二个等式**：$\mathcal{COP}_n^* = \mathcal{CP}_n$。

    由于 $\mathcal{CP}_n$ 是**闭凸锥**（定理 45A.1），双对偶定理（bipolar theorem）保证
    $$\mathcal{CP}_n^{**} = \mathcal{CP}_n.$$

    而 $\mathcal{CP}_n^{**} = (\mathcal{CP}_n^*)^* = \mathcal{COP}_n^*$。因此 $\mathcal{COP}_n^* = \mathcal{CP}_n$。 $\blacksquare$

!!! theorem "定理 45B.4 (锥层次结构)"
    对 $n \times n$ 对称矩阵，有如下锥的包含关系及对偶关系：

    $$\mathcal{CP}_n \subseteq \mathcal{DNN}_n = \mathcal{S}_n^+ \cap \mathcal{N}_n$$

    $$\mathcal{S}_n^+ + \mathcal{N}_n \subseteq \mathcal{COP}_n$$

    对偶关系：$\mathcal{DNN}_n^* = \mathcal{S}_n^+ + \mathcal{N}_n$（利用 $(\mathcal{K}_1 \cap \mathcal{K}_2)^* = \mathcal{K}_1^* + \mathcal{K}_2^*$）。

    因此
    $$\mathcal{CP}_n \subseteq \mathcal{DNN}_n \quad \Leftrightarrow \quad \mathcal{DNN}_n^* \subseteq \mathcal{CP}_n^* \quad \Leftrightarrow \quad \mathcal{S}_n^+ + \mathcal{N}_n \subseteq \mathcal{COP}_n.$$

    当 $n \leq 4$ 时，所有包含变为等号：
    $$\mathcal{CP}_n = \mathcal{DNN}_n, \qquad \mathcal{S}_n^+ + \mathcal{N}_n = \mathcal{COP}_n.$$

    当 $n \geq 5$ 时，所有包含严格。

!!! example "例 45B.2 (锥层次图)"
    ```
    CP_n  ⊆  DNN_n = S_n^+ ∩ N_n
     ↕ 对偶      ↕ 对偶
    COP_n ⊇  S_n^+ + N_n
    ```

    当 $n \leq 4$ 时，$\subseteq$ 和 $\supseteq$ 变为 $=$。
    当 $n \geq 5$ 时，$\subsetneq$ 和 $\supsetneq$。

    维数：$\mathcal{S}_n$ 的维数为 $n(n+1)/2$，所有上述锥在 $\mathcal{S}_n$ 中都是满维的。

---

## 45B.5 COP 锥的极端射线与面结构

<div class="context-flow" markdown>

**核心问题**：$\mathcal{COP}_n$ 的极端射线是什么？对 $n \geq 5$，是否有新类型的极端射线？

</div>

!!! theorem "定理 45B.5 (COP 锥的极端射线)"
    1. 对 $n \leq 4$，$\mathcal{COP}_n = \mathcal{S}_n^+ + \mathcal{N}_n$ 的极端射线恰好是以下两类：
        - $\mathcal{S}_n^+$ 的极端射线：形如 $vv^T$（$v \in \mathbb{R}^n$），即秩一半正定矩阵。
        - $\mathcal{N}_n$ 的极端射线：形如 $E_{ij} + E_{ji}$（$i \neq j$）或 $E_{ii}$，其中 $E_{ij}$ 是标准基矩阵。

    2. 对 $n \geq 5$，$\mathcal{COP}_n$ 有**不属于** $\mathcal{S}_n^+ + \mathcal{N}_n$ 的极端射线——称为**特殊极端射线**（exceptional extreme rays）。Horn 矩阵所在的射线就是一个例子。

    3. 对 $n = 5$，$\mathcal{COP}_5$ 的所有极端射线已被完全分类（Hildebrand, 2012）。

!!! theorem "定理 45B.6 (Bomze 的 COP 边界刻画)"
    Bomze (2009) 给出了 $\mathcal{COP}_n$ 的边界的刻画：

    对称矩阵 $A$ 位于 $\mathcal{COP}_n$ 的**边界**上当且仅当存在非零向量 $x \geq 0$ 使得 $x^T A x = 0$。

    这个简单的刻画有以下推论：

    1. $A \in \operatorname{int}(\mathcal{COP}_n)$（余正锥的内部）当且仅当 $A$ 严格余正。
    2. $A \in \partial(\mathcal{COP}_n)$（边界）当且仅当存在 $x \geq 0$（$x \neq 0$）使得 $x^T A x = 0$。
    3. 对边界上的 $A$，集合 $\{x \geq 0 : x^T A x = 0\}$ 称为 $A$ 的**零集**（zero set），它的结构决定了 $A$ 在面格中的位置。

??? proof "证明"
    **充分性**：若存在 $x \geq 0$（$x \neq 0$）使得 $x^T A x = 0$，则 $A$ 不是严格余正的。对任意 $\epsilon > 0$，取 $A_\epsilon = A - \epsilon xx^T$，则 $x^T A_\epsilon x = -\epsilon \|x\|^4 < 0$，但 $x \geq 0$，因此 $A_\epsilon \notin \mathcal{COP}_n$。这表明 $A$ 的任意邻域内都有不在 $\mathcal{COP}_n$ 中的元素，因此 $A \in \partial(\mathcal{COP}_n)$。

    **必要性**：若 $A \in \partial(\mathcal{COP}_n)$，则存在序列 $A_m \notin \mathcal{COP}_n$ 使得 $A_m \to A$。对每个 $A_m$，存在 $x_m \geq 0$（$\|x_m\| = 1$）使得 $x_m^T A_m x_m < 0$。由紧性，$x_m$ 有收敛子列 $x_{m_k} \to x^*$。则 $x^* \geq 0$，$\|x^*\| = 1$，$(x^*)^T A x^* = \lim (x_{m_k})^T A_{m_k} x_{m_k} \leq 0$。由 $A \in \mathcal{COP}_n$，$(x^*)^T A x^* \geq 0$。因此 $(x^*)^T A x^* = 0$。 $\blacksquare$

---

## 45B.6 独立数的余正规划

<div class="context-flow" markdown>

**核心问题**：如何利用余正矩阵将组合优化问题表述为锥规划？

</div>

!!! definition "定义 45B.3 (余正规划)"
    **余正规划**（copositive program）是如下形式的优化问题：

    **原始形式**（完全正锥上的规划）：
    $$\min_{X \in \mathcal{S}_n} \langle C, X \rangle \quad \text{s.t.} \quad \langle A_i, X \rangle = b_i, \; i = 1, \ldots, m, \quad X \in \mathcal{CP}_n.$$

    **对偶形式**（余正锥上的规划）：
    $$\max_{y \in \mathbb{R}^m} b^T y \quad \text{s.t.} \quad C - \sum_{i=1}^{m} y_i A_i \in \mathcal{COP}_n.$$

    余正规划是半定规划（SDP）的推广：将半正定锥 $\mathcal{S}_n^+$ 替换为 $\mathcal{CP}_n$ 或 $\mathcal{COP}_n$。

!!! theorem "定理 45B.7 (独立数的余正规划 — de Klerk-Pasechnik)"
    设 $G = (V, E)$ 是无向图，$|V| = n$，$A_G$ 是 $G$ 的邻接矩阵。$G$ 的**独立数**（independence number，也称稳定数）$\alpha(G)$——最大独立集的大小——可以精确表述为余正规划：

    $$\frac{1}{\alpha(G)} = \min\left\{t : tI + A_G \in \mathcal{COP}_n\right\}.$$

    等价地，
    $$\alpha(G) = \max\left\{\langle J, X \rangle : \langle I + A_G, X \rangle = 1, \; X \in \mathcal{CP}_n\right\},$$
    其中 $J = \mathbf{1}\mathbf{1}^T$。

??? proof "证明"
    **关键思路**：独立集的特征向量表示。

    $S \subseteq V$ 是独立集当且仅当 $S$ 中顶点两两不相邻：$\mathbf{1}_S^T A_G \mathbf{1}_S = 0$。

    **下界**：设 $S$ 是大小为 $k$ 的独立集。令 $x = \frac{1}{k}\mathbf{1}_S \geq 0$。则
    $$x^T(tI + A_G)x = t\|x\|^2 + x^T A_G x = \frac{t}{k} + 0 = \frac{t}{k}.$$
    若 $tI + A_G \in \mathcal{COP}_n$，则 $t/k \geq 0$，即 $t \geq 0$。

    更重要的是，对标准单纯形上的最优化，考虑
    $$\min\{x^T(tI + A_G)x : x \geq 0, \mathbf{1}^T x = 1\}.$$

    取 $x = \frac{1}{|S|}\mathbf{1}_S$（$S$ 为最大独立集），最小值为 $t/\alpha(G)$。

    **精确值**：$tI + A_G \in \mathcal{COP}_n$ 当且仅当 $\min\{x^T(tI + A_G)x : x \geq 0, \mathbf{1}^T x = 1\} \geq 0$。

    Motzkin-Straus 定理告诉我们，对非负正交锥上的标准单纯形，
    $$\min_{x \in \Delta_n} x^T(tI + A_G)x = t \cdot \min_{x \in \Delta_n} \|x\|^2 + \min_{x \in \Delta_n} x^T A_G x,$$
    但这不是精确的分解（两个最小值不一定在同一点取到）。

    正确的论证是：对 $t = 1/\alpha(G)$，$tI + A_G$ 的最小值在独立集的均匀分布上取到，值为 0。对 $t < 1/\alpha(G)$，存在独立集使 $x^T(tI + A_G)x < 0$。因此
    $$\min\{t : tI + A_G \in \mathcal{COP}_n\} = \frac{1}{\alpha(G)}. \quad \blacksquare$$

!!! example "例 45B.3 (Petersen 图)"
    Petersen 图 $P$ 有 10 个顶点，$\alpha(P) = 4$。

    - **SDP 松弛**（Lovasz theta 函数）：$\vartheta(P) = 4$，恰好等于 $\alpha(P)$。
    - **余正规划**：$1/\alpha(P) = \min\{t : tI + A_P \in \mathcal{COP}_{10}\} = 1/4$。

    对一般图，SDP（Lovasz theta）给出上界 $\alpha(G) \leq \vartheta(G)$，而余正规划给出**精确值**。代价是 COP 约束是 NP-hard 的。

---

## 45B.7 NP-hard 问题的余正编码

<div class="context-flow" markdown>

**核心问题**：哪些 NP-hard 组合优化问题可以表述为余正规划？如何给出具体的公式？

</div>

!!! theorem "定理 45B.8 (最大独立集的余正公式)"
    设 $G = (V, E)$，$|V| = n$，邻接矩阵 $A_G$。则
    $$\alpha(G) = \max\left\{\mathbf{1}^T X \mathbf{1} : \operatorname{tr}(X) = 1, \; (A_G)_{ij} x_{ij} = 0 \; \forall (i,j) \in E, \; X \in \mathcal{CP}_n\right\}.$$

    等价地，
    $$\alpha(G) = \max\left\{\langle J, X \rangle : \langle I + A_G, X \rangle = 1, \; X \in \mathcal{CP}_n\right\}.$$

!!! theorem "定理 45B.9 (最大团的余正公式)"
    $G$ 的**最大团数** $\omega(G)$ 等于补图 $\bar{G}$ 的独立数 $\alpha(\bar{G})$，因此
    $$\omega(G) = \alpha(\bar{G}) = \max\left\{\langle J, X \rangle : \langle I + A_{\bar{G}}, X \rangle = 1, \; X \in \mathcal{CP}_n\right\}.$$

    其中 $A_{\bar{G}} = J - I - A_G$ 是补图的邻接矩阵。

!!! theorem "定理 45B.10 (图着色数的余正公式)"
    $G$ 的**色数** $\chi(G)$（将顶点着色使相邻顶点不同色的最少颜色数）满足
    $$\chi(G) = \min\left\{k : \exists X \in \mathcal{CP}_{kn}, \; \text{满足着色约束}\right\}.$$

    更精确地，将 $kn \times kn$ 矩阵 $X$ 分为 $k \times k$ 的块结构 $X = [X_{ab}]_{a,b=1}^k$（每块 $n \times n$），则
    $$\chi(G) \leq k \quad \Leftrightarrow \quad \exists X \in \mathcal{CP}_{kn} \text{ 满足：}$$
    $$\sum_{a=1}^k (X_{aa})_{ii} = 1 \; \forall i, \qquad (X_{aa})_{ij} = 0 \; \forall (i,j) \in E, \; \forall a.$$

    直观地，$(X_{aa})_{ii}$ 表示顶点 $i$ 被着第 $a$ 种颜色的"概率"。完全正约束保证解可以解释为有效的着色方案。

!!! theorem "定理 45B.11 (MAX-CUT 的余正公式)"
    设 $G = (V, E)$ 是带权无向图，权重矩阵 $W = [w_{ij}]$（$w_{ij} \geq 0$）。**MAX-CUT 问题**是将顶点分为两组 $S$ 和 $\bar{S}$ 使得割边的权重之和最大：
    $$\operatorname{MAXCUT}(G) = \max_{S \subseteq V} \sum_{\substack{i \in S, j \notin S}} w_{ij}.$$

    MAX-CUT 可以表述为余正规划：

    **二次公式**：令 $x_i = 1$（$i \in S$）或 $x_i = 0$（$i \notin S$），则
    $$\operatorname{MAXCUT}(G) = \max\left\{\sum_{(i,j) \in E} w_{ij}(x_i - x_j)^2 : x \in \{0, 1\}^n\right\}.$$

    **余正松弛**：将 $x \in \{0, 1\}^n$ 松弛为 $X \in \mathcal{CP}_n$（利用 $x_i \in \{0, 1\} \Leftrightarrow x_i^2 = x_i$），得到

    $$\operatorname{MAXCUT}(G) = \max\left\{\frac{1}{4}\langle L, X \rangle : X_{ii} = 1 \; \forall i, \; X + J \in \mathcal{CP}_{n}\right\},$$

    其中 $L = D - W$ 是 $G$ 的 Laplacian 矩阵，$D = \operatorname{diag}(\sum_j w_{ij})$。

    更标准的形式（使用 $\pm 1$ 变量 $y_i = 2x_i - 1$）：
    $$\operatorname{MAXCUT}(G) = \max\left\{\frac{1}{4}\langle L, Y \rangle : Y_{ii} = 1 \; \forall i, \; Y \in \mathcal{COP}_n^*\right\}.$$

!!! example "例 45B.4 (标准二次规划)"
    **标准二次规划**（standard quadratic program, StQP）是在标准单纯形上的齐次二次优化：
    $$\min\{x^T Q x : x \geq 0, \; \mathbf{1}^T x = 1\}.$$

    这可以等价表述为余正规划：
    $$\min\{x^T Q x : x \geq 0, \; \mathbf{1}^T x = 1\} = \min\{\langle Q, X \rangle : X \in \mathcal{CP}_n, \; \langle J, X \rangle = 1\}.$$

    StQP 即使对不定矩阵 $Q$ 也是 NP-hard 的。这说明余正规划的 NP-hard 性质"来源于"标准单纯形上二次优化的困难性。

!!! theorem "定理 45B.12 (NP-hard 问题的统一框架)"
    以下组合优化问题都可以精确表述为余正规划：

    1. **最大独立集**（$\alpha(G)$）。
    2. **最大团**（$\omega(G)$）。
    3. **图着色数**（$\chi(G)$）。
    4. **MAX-CUT**。
    5. **二次分配问题**（QAP）。
    6. **标准二次规划**。

    这说明余正规划在计算复杂性意义上是"万能"的——它包含了所有 NP-hard 问题。余正规划的计算困难性完全集中在余正锥约束上。

---

## 45B.8 co-NP 完全性

<div class="context-flow" markdown>

**核心问题**：判定余正性的计算复杂度是什么？

</div>

!!! theorem "定理 45B.13 (余正性判定的复杂度)"
    判定一个给定的对称矩阵 $A$ 是否余正是 **co-NP 完全**问题。

    等价地，判定一个矩阵是否**不**余正（即存在 $x \geq 0$ 使得 $x^T A x < 0$）是 **NP 完全**问题。

??? proof "证明"
    **在 co-NP 中**：$A$ 不余正的证据是一个非负向量 $x \geq 0$（$x \neq 0$）使得 $x^T A x < 0$。这个证据可以在多项式时间内验证：计算 $x^T A x$ 需要 $O(n^2)$ 次运算。

    **co-NP 难**：通过从独立数问题归约。

    由定理 45B.7，$\alpha(G) \geq k$ 当且仅当 $(1/k) I + A_G$ **不是**严格余正的——即存在 $x \geq 0$（$\|x\| = 1$）使得 $x^T((1/k)I + A_G)x = 0$。

    更精确地：$(1/k)I + A_G \in \mathcal{COP}_n$ 但 $(1/k)I + A_G \notin \operatorname{int}(\mathcal{COP}_n)$ 当且仅当 $\alpha(G) \geq k$。

    设 $t^* = 1/\alpha(G)$。则 $t^* I + A_G \in \mathcal{COP}_n$ 但 $(t^* - \epsilon)I + A_G \notin \mathcal{COP}_n$（对任意 $\epsilon > 0$）。

    因此，判定 "$tI + A_G$ 是否余正"等价于判定 "$t \geq 1/\alpha(G)$"。由于判定 $\alpha(G) \geq k$（等价于 $1/\alpha(G) \leq 1/k$）是 NP 完全的（独立集问题），余正性判定至少同等困难。 $\blacksquare$

---

## 45B.9 Parrilo 的 SDP 层次逼近

<div class="context-flow" markdown>

**核心问题**：如何用可计算的方法逼近余正锥？

</div>

!!! theorem "定理 45B.14 (Parrilo 的 SOS 层次)"
    Parrilo (2000) 和 de Klerk-Pasechnik (2002) 提出了余正锥的**层次化 SDP 逼近**。

    定义锥的层次
    $$\mathcal{K}_n^{(r)} = \left\{A \in \mathcal{S}_n : \left(\sum_{i=1}^n x_i^2\right)^r \cdot x^T A x \text{ 是 SOS 多项式}\right\},$$

    其中 SOS（sum of squares）表示多项式可以写成若干多项式平方之和。则：

    1. $\mathcal{S}_n^+ + \mathcal{N}_n = \mathcal{K}_n^{(0)} \subseteq \mathcal{K}_n^{(1)} \subseteq \cdots \subseteq \mathcal{COP}_n$。
    2. $\overline{\bigcup_{r=0}^{\infty} \mathcal{K}_n^{(r)}} = \mathcal{COP}_n$（闭包等于余正锥）。
    3. 每个 $\mathcal{K}_n^{(r)}$ 可以用有限维 SDP 检验，SDP 的大小为 $\binom{n + 2r + 2}{2r + 2}$。

??? proof "证明"
    **(1)** $\mathcal{K}_n^{(0)}$ 的定义：$A \in \mathcal{K}_n^{(0)}$ 当且仅当 $x^T A x$ 是 SOS 齐次多项式。

    齐次二次 SOS 多项式 $x^T A x$ 意味着 $x^T A x = \sum_l (c_l^T x)^2$，即 $A = \sum_l c_l c_l^T$，即 $A \succeq 0$。

    但我们需要考虑 $x \geq 0$ 的约束。实际上 $\mathcal{K}_n^{(0)}$ 的正确定义考虑了非负变量：利用替换 $x_i = y_i^2$，$x^T A x = \sum_{i,j} a_{ij} y_i^2 y_j^2$，这是 SOS 当且仅当 $A \in \mathcal{S}_n^+ + \mathcal{N}_n$（Diananda 型结果）。

    **(2)** 关键引理：若 $A$ 余正，则 $(\sum x_i^2)^r \cdot x^T A x$ 在非负正交锥上严格正（对充分大的 $r$），并且可以被 SOS 逼近。闭包等式来自正多项式在非负正交锥上的 SOS 表示的 Schmudgen/Putinar 型定理。

    **(3)** 多项式 $p(x) = (\sum x_i^2)^r \cdot x^T A x$ 是 $2(r+1)$ 次齐次多项式，其 SOS 检验等价于一个关于 $A$ 的线性矩阵不等式（LMI），即 SDP 可行性问题。SDP 的大小由 $2(r+1)$ 次齐次多项式的单项式个数决定。 $\blacksquare$

!!! example "例 45B.5 (层次逼近的效果)"
    对 Horn 矩阵 $H$：

    - $r = 0$：$H \notin \mathcal{K}_5^{(0)} = \mathcal{S}_5^+ + \mathcal{N}_5$。
    - $r = 1$：$H \in \mathcal{K}_5^{(1)}$。即 $(\sum x_i^2) \cdot x^T H x$ 是 SOS 多项式。

    因此 Parrilo 层次的第一层就足以认证 Horn 矩阵的余正性。但对更复杂的矩阵，可能需要更高层次。

    一般地，对 $n \times n$ 矩阵 $A \in \mathcal{COP}_n$，在第 $r$ 层的 SDP 大小约为 $O(n^{2r+2})$（多项式大小随 $r$ 指数增长），因此层次逼近的实际可行性受 $r$ 限制。

---

## 45B.10 与线性互补问题的联系

<div class="context-flow" markdown>

**核心问题**：余正矩阵与线性互补问题有什么关系？

</div>

!!! definition "定义 45B.4 (线性互补问题)"
    给定矩阵 $M \in \mathbb{R}^{n \times n}$ 和向量 $q \in \mathbb{R}^n$，**线性互补问题** LCP$(M, q)$ 是求向量 $z \in \mathbb{R}^n$ 满足
    $$z \geq 0, \quad w = Mz + q \geq 0, \quad z^T w = 0.$$

    互补条件 $z^T w = 0$ 意味着对每个 $i$，$z_i = 0$ 或 $w_i = 0$（或两者皆为零）。

!!! theorem "定理 45B.15 (余正矩阵与 LCP)"
    1. 若 $M$ 是余正矩阵，则 LCP$(M, q)$ 对所有 $q$ 使得可行域非空时有解。

    2. 更精确地，$M$ 余正当且仅当对所有 $q$ 使得 LCP$(M, q)$ 可行，该问题有解。

    3. 余正矩阵的 LCP 与标准二次规划的 KKT 条件密切相关：$\min\{x^T M x + 2q^T x : x \geq 0\}$ 的 KKT 条件恰好是 LCP$(M, q)$（当 $M$ 对称时）。

??? proof "证明"
    **(3)** 考虑问题 $\min\{x^T M x + 2q^T x : x \geq 0\}$（$M$ 对称）。KKT 条件为：
    $$2Mx + 2q - \lambda = 0, \quad x \geq 0, \quad \lambda \geq 0, \quad x^T \lambda = 0.$$

    令 $w = \lambda/2 = Mx + q$，则 $x \geq 0$，$w \geq 0$，$x^T w = 0$，即 LCP$(M, q)$。

    若 $M$ 余正且 $q$ 使得可行域非空，则二次目标函数在可行域上有下界（因为 $x^T M x \geq 0$ 对 $x \geq 0$），因此最小值存在，KKT 条件有解。 $\blacksquare$

!!! example "例 45B.6 (LCP 的应用)"
    线性互补问题出现在以下场景中：

    1. **接触力学**：刚体间的接触力满足互补条件（接触力为零或间隙为零）。
    2. **纳什均衡**：多人博弈的均衡条件可以表述为 LCP。
    3. **期权定价**：美式期权的最优行权条件是 LCP。

    余正矩阵的 LCP 理论将这些应用与余正锥的几何结构联系起来。

---

## 45B.11 Sturm 的 SDP 可表示性

<div class="context-flow" markdown>

**核心问题**：余正锥能否表示为半定锥的投影？

</div>

!!! theorem "定理 45B.16 (Sturm 的 SDP 可表示性结果)"
    Sturm (2000) 证明了以下结果：

    1. 对**固定的** $n$，$\mathcal{COP}_n$ 是一个**谱面体的影子**（shadow of a spectrahedron）的闭包——即存在有限维 SDP 描述。但 SDP 的大小随 $n$ 超多项式增长。

    2. 不存在**多项式大小**的 SDP 表示 $\mathcal{COP}_n$（对一般 $n$），除非 $P = NP$。

    3. 更精确地：对固定 $r$，Parrilo 层次的第 $r$ 层 $\mathcal{K}_n^{(r)}$ 有确切的 SDP 表示，大小为关于 $n$ 的多项式。但每层都严格弱于 $\mathcal{COP}_n$（对 $n \geq 5$），需要无穷层的极限才能恢复 $\mathcal{COP}_n$。

    这些结果说明了为什么余正规划虽然在理论上可以统一所有 NP-hard 问题，但在实践中不能替代专门的组合优化算法。余正锥的"非 SDP 可表示性"本质上等价于 $P \neq NP$ 的假设。

---

## 45B.12 习题

!!! example "例 45B.7"
    验证以下矩阵是否余正：
    $$A = \begin{pmatrix} 1 & -1 & 1 \\ -1 & 1 & -1 \\ 1 & -1 & 1 \end{pmatrix}.$$
    **提示**：直接计算 $x^T A x$ 并分析其在 $x \geq 0$ 上的符号。

!!! example "例 45B.8"
    证明：对 $n = 2$，$\mathcal{COP}_2 = \mathcal{S}_2^+ + \mathcal{N}_2$。

    **提示**：对 $A = \begin{pmatrix} a & c \\ c & b \end{pmatrix}$ 余正，分情况讨论 $c \geq 0$ 和 $c < 0$。

!!! example "例 45B.9"
    设 $G = C_5$（五圈图），$A_G$ 是其邻接矩阵。计算 $1/\alpha(G)$ 并验证 $(1/\alpha(G))I + A_G$ 是余正的。

!!! example "例 45B.10"
    证明：若 $A \in \mathcal{COP}_n$ 且 $B \in \mathcal{COP}_n$，则 $A + B \in \mathcal{COP}_n$ 且 $\alpha A \in \mathcal{COP}_n$（$\alpha \geq 0$）。即 $\mathcal{COP}_n$ 是凸锥。

!!! example "例 45B.11"
    设 $A \in \mathcal{COP}_n$。证明 $A$ 的所有主子矩阵也是余正的。即若 $S \subseteq \{1, \ldots, n\}$，则 $A_S \in \mathcal{COP}_{|S|}$。

!!! example "例 45B.12"
    证明对偶关系 $\mathcal{DNN}_n^* = \mathcal{S}_n^+ + \mathcal{N}_n$。

    **提示**：利用凸锥对偶的性质 $(\mathcal{K}_1 \cap \mathcal{K}_2)^* = \overline{\mathcal{K}_1^* + \mathcal{K}_2^*}$，其中 $\mathcal{K}_1 = \mathcal{S}_n^+$（自对偶），$\mathcal{K}_2 = \mathcal{N}_n$（自对偶），并验证 $\mathcal{S}_n^+ + \mathcal{N}_n$ 已经是闭锥。

!!! example "例 45B.13"
    对下列图 $G$，将 $\alpha(G)$ 表述为余正规划并求值：

    1. $G = K_3$（完全图），$\alpha(K_3) = 1$。
    2. $G = C_4$（四圈），$\alpha(C_4) = 2$。
    3. $G = P_3$（路径图，3 个顶点），$\alpha(P_3) = 2$。

!!! example "例 45B.14"
    对完全图 $G = K_n$，验证 MAX-CUT 的余正公式给出 $\operatorname{MAXCUT}(K_n) = \lfloor n^2/4 \rfloor$。

    **提示**：$K_n$ 的 MAX-CUT 是将 $n$ 个顶点分为大小尽可能相等的两组。

!!! example "例 45B.15"
    证明 Parrilo 层次的第零层 $\mathcal{K}_n^{(0)} = \mathcal{S}_n^+ + \mathcal{N}_n$。

    **提示**：齐次二次多项式 $\sum a_{ij} y_i^2 y_j^2$（替换 $x_i = y_i^2$）是 SOS 当且仅当对应的矩阵属于 $\mathcal{S}_n^+ + \mathcal{N}_n$。

!!! example "例 45B.16"
    验证 Horn 矩阵 $H$ 属于 $\mathcal{K}_5^{(1)}$ 但不属于 $\mathcal{K}_5^{(0)}$。即 $(\sum x_i^2) \cdot x^T H x$ 是 SOS 多项式但 $x^T H x$ 在变量替换后不是 SOS。

!!! example "例 45B.17"
    设 $M$ 是 $3 \times 3$ 余正矩阵。证明 LCP$(M, q)$ 对 $q = (1, 1, 1)^T$ 有解，并给出一种求解方法。

!!! example "例 45B.18"
    （开放问题讨论）是否存在多项式大小的 SDP 证书来证明一个给定的 $n \times n$ 矩阵是余正的？讨论这个问题与 $P$ vs $NP$ 问题的关系。

!!! example "例 45B.19"
    设 $G$ 是 Petersen 图。将图着色数 $\chi(G) = 3$ 表述为余正规划，并解释为什么直接求解该规划是困难的。

!!! example "例 45B.20"
    证明：标准二次规划 $\min\{x^T Q x : x \geq 0, \mathbf{1}^T x = 1\}$ 可以等价表述为完全正锥上的线性规划。给出对偶问题的余正锥表述。

---

## 练习题

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

## 本章小结

余正规划构成了连续优化与离散优化之间的桥梁：

1. **受限正性**：定义了余正锥作为将半正定要求限制在非负正交锥上的推广。
2. **精确编码**：确立了 COP/CP 规划为计算机科学中最难的问题（NP-hard）提供了精确的数学公式。
3. **层次化计算**：引入了 SOS 方法，为无法直接计算的余正约束提供了单调收敛的 SDP 逼近序列。
4. **对偶锚点**：利用 CP-COP 对偶性统一了非负矩阵分解与全局最优化问题的研究。

