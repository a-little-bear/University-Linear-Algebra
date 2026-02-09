# 第 65A 章 符号模式与定性矩阵分析

<div class="context-flow" markdown>

**前置**：行列式(Ch3) · 图论基础(Ch27) · 特征值(Ch6) · 组合矩阵结构(Ch65B)

**本章脉络**：零-非零模式与有向图 → 模式决定的性质 → 符号模式与 SNS 矩阵（含完整证明）→ 符号可解系统（含证明）→ 定性稳定（含证明）→ 生态学应用 → 完全不可分解矩阵 → 近似可约矩阵 → 潜在幂零模式 → 谱任意模式（幂零-Jacobi 方法）→ 要求/允许问题

**延伸**：定性矩阵分析在数学生态学、经济学（定性比较静态）和控制论（结构可控性）中有核心应用

</div>

组合矩阵论研究矩阵的**组合结构**（零-非零模式、符号模式、图结构）如何影响其**代数性质**（秩、行列式、特征值、稳定性）。这是离散数学与线性代数的自然交叉点。本章聚焦于**符号模式**的代数性质：什么时候矩阵的符号信息足以确定行列式的非零性、线性系统解的符号或动力系统的稳定性。

定性矩阵分析的动机来自实际应用：在生态学中，物种之间的交互往往只知道定性类型（捕食、竞争、共生），而不知道精确的交互强度。在经济学中，变量之间的正负关系往往比精确数值更容易确定。本章的理论为这些应用提供了严格的数学基础。

---

## 65A.1 零-非零模式与有向图

### 基本定义

!!! definition "定义 65A.1 (零-非零模式)"
    矩阵 $A = (a_{ij}) \in M_{m \times n}(\mathbb{R})$ 的**零-非零模式**是 $\{0, *\}$ 矩阵 $\mathcal{Z}(A)$：
    $$
    \mathcal{Z}(A)_{ij} = \begin{cases} * & \text{若 } a_{ij} \neq 0 \\ 0 & \text{若 } a_{ij} = 0 \end{cases}
    $$

    模式 $\mathcal{P}$ 的**定性类**为 $Q(\mathcal{P}) = \{A \in M_{m \times n}(\mathbb{R}) : \mathcal{Z}(A) = \mathcal{P}\}$。

### 模式与有向图

!!! definition "定义 65A.2 (模式的有向图)"
    $n \times n$ 模式 $\mathcal{P}$ 的**有向图** $D(\mathcal{P})$：顶点集 $V = \{1, \ldots, n\}$，有向边 $(i, j) \in E$ 当且仅当 $\mathcal{P}_{ij} = *$。对角元 $\mathcal{P}_{ii} = *$ 对应自环。

!!! example "例 65A.1"
    模式 $\mathcal{P} = \begin{pmatrix} * & * \\ 0 & * \end{pmatrix}$ 的定性类为
    $$
    Q(\mathcal{P}) = \left\{\begin{pmatrix} a & b \\ 0 & d \end{pmatrix} : a, b, d \neq 0\right\}
    $$
    有向图 $D(\mathcal{P})$：顶点 $\{1, 2\}$，边 $(1,1), (1,2), (2,2)$。此类中所有矩阵都是上三角且 $\det = ad \neq 0$，因此都是非奇异的。

### 模式能决定的性质

!!! theorem "定理 65A.1 (模式不变的性质)"
    以下性质**不能**仅由零-非零模式决定：秩、可逆性、特征值具体数值。

    以下性质**可以**由模式决定或约束：秩的上界、不可约性（等价于 $D(\mathcal{P})$ 强连通）、在附加符号条件下的非奇异性。

??? proof "证明"
    **秩不由模式决定**：模式 $\begin{pmatrix} * & * \\ * & * \end{pmatrix}$ 包含秩 1 矩阵（如 $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$）和秩 2 矩阵（如 $\begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$）。

    **不可约性由模式决定**：矩阵 $A$ 不可约当且仅当不存在置换矩阵 $P$ 使得 $P^TAP$ 是分块上三角的，这等价于有向图 $D(\mathcal{Z}(A))$ 是强连通的。强连通性仅依赖于边集（即非零模式），不依赖于具体数值。$\blacksquare$

!!! example "例 65A.2"
    模式 $\begin{pmatrix} * & 0 \\ 0 & * \end{pmatrix}$ 的定性类中所有矩阵的秩均为 2（因为对角元素非零且无非对角元素）。但模式 $\begin{pmatrix} * & * \\ * & * \end{pmatrix}$ 的类中秩可以为 1 或 2。

---

## 65A.2 符号非奇异矩阵

### 符号模式

!!! definition "定义 65A.3 (符号模式)"
    矩阵 $A$ 的**符号模式** $\text{sgn}(A)$ 是 $\{+, -, 0\}$ 矩阵：
    $$
    \text{sgn}(A)_{ij} = \begin{cases} + & \text{若 } a_{ij} > 0 \\ - & \text{若 } a_{ij} < 0 \\ 0 & \text{若 } a_{ij} = 0 \end{cases}
    $$
    符号模式 $\mathcal{S}$ 的**符号类**为 $Q(\mathcal{S}) = \{A \in M_n(\mathbb{R}) : \text{sgn}(A) = \mathcal{S}\}$。

### SNS 矩阵的刻画

!!! definition "定义 65A.4 (符号非奇异矩阵)"
    符号模式 $\mathcal{S}$ 称为**符号非奇异**（sign nonsingular, SNS），若 $Q(\mathcal{S})$ 中的每个矩阵都是非奇异的。等价地，$\det(A)$ 的符号由 $\mathcal{S}$ 唯一确定。

!!! theorem "定理 65A.2 (SNS 的完整刻画)"
    符号模式 $\mathcal{S}$ 是 SNS 的，当且仅当行列式展开
    $$
    \det(A) = \sum_{\sigma \in S_n} \text{sgn}(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}
    $$
    中的所有非零项具有**相同的符号**。

??? proof "证明"
    **充分性**：假设所有非零项同号。行列式的每一项 $\text{sgn}(\sigma) \prod_i a_{i,\sigma(i)}$ 的符号由以下两部分决定：

    - 置换 $\sigma$ 的符号 $\text{sgn}(\sigma) = \pm 1$；
    - 乘积 $\prod_i a_{i,\sigma(i)}$ 的符号，由各 $a_{i,\sigma(i)}$ 的符号（来自 $\mathcal{S}$）决定。

    若 $\prod_i a_{i,\sigma(i)} = 0$（即某个 $s_{i,\sigma(i)} = 0$），该项为零。非零项的符号完全由 $\mathcal{S}$ 和 $\sigma$ 确定。

    所有非零项同号意味着行列式是同号非零项之和。对于 $Q(\mathcal{S})$ 中的任意 $A$，每个非零项是**正数之和**（或全是负数之和），具体的正数值依赖于 $|a_{ij}|$，但符号不变。因此 $\det(A) \neq 0$（同号正数之和为正，同号负数之和为负）。

    **必要性**：假设存在两个非零项符号相反，分别对应置换 $\sigma_1$ 和 $\sigma_2$。设
    $$
    t_k(A) = \text{sgn}(\sigma_k) \prod_{i=1}^n a_{i,\sigma_k(i)}, \quad k = 1, 2
    $$
    其中 $t_1(A) > 0$ 和 $t_2(A) < 0$（对某个特定的 $A \in Q(\mathcal{S})$）。

    通过连续地缩放 $Q(\mathcal{S})$ 中矩阵的元素（保持符号不变），可以调节各项的大小。具体地，将 $\sigma_1$ 对应的元素放大、$\sigma_2$ 对应的元素放大，使得 $|t_1| = |t_2|$。更形式化地：

    定义 $A(\epsilon) \in Q(\mathcal{S})$，使得 $\sigma_1$ 和 $\sigma_2$ 对应行的元素按参数 $\epsilon$ 变化。由于 $t_1$ 和 $t_2$ 符号相反，通过中值定理，存在参数值使得 $\det(A(\epsilon)) = 0$。

    更严格地：固定 $A \in Q(\mathcal{S})$。定义 $A(s)$，其 $(i, j)$ 元素为
    $$
    a_{ij}(s) = \begin{cases} s \cdot a_{i,\sigma_1(i)} & \text{若 } j = \sigma_1(i) \text{ 且 } j \neq \sigma_2(i) \\ a_{ij} & \text{否则} \end{cases}
    $$
    当 $s \to 0$ 时，$t_1 \to 0$ 而 $t_2$ 保持有限，$\det(A(s)) \to t_2(A) + \text{其他项}$ 可能变号。通过仔细构造，可以使 $\det = 0$。

    因此，若非零项有不同符号，$\mathcal{S}$ 不是 SNS。$\blacksquare$

!!! example "例 65A.3"
    $\mathcal{S} = \begin{pmatrix} + & + \\ 0 & + \end{pmatrix}$：$\det(A) = a_{11}a_{22}$，唯一非零项为正。SNS。

    $\mathcal{S}' = \begin{pmatrix} + & + \\ + & + \end{pmatrix}$：$\det(A) = a_{11}a_{22} - a_{12}a_{21}$，第一项正、第二项负。不是 SNS。

!!! example "例 65A.4"
    $3 \times 3$ 符号模式 $\mathcal{S} = \begin{pmatrix} + & - & 0 \\ 0 & + & - \\ - & 0 & + \end{pmatrix}$。

    非零行列式项：
    - $\sigma = \text{id}$：$a_{11}a_{22}a_{33} > 0$，$\text{sgn}(\sigma) = +1$，贡献 $> 0$；
    - $\sigma = (132)$：$a_{12}a_{23}a_{31}$，$(-)(-)(-) < 0$，$\text{sgn}(\sigma) = +1$，贡献 $< 0$。

    两项符号不同，$\mathcal{S}$ 不是 SNS。

---

## 65A.3 符号可解系统

!!! definition "定义 65A.5 (符号可解系统)"
    线性系统 $Ax = b$ 称为**符号可解的**，若对所有 $A' \in Q(\text{sgn}(A))$ 和 $b' \in Q(\text{sgn}(b))$，解 $x' = (A')^{-1}b'$ 的符号模式与 $x = A^{-1}b$ 相同。

!!! theorem "定理 65A.3 (符号可解性的充要条件)"
    系统 $Ax = b$ 符号可解，当且仅当：

    1. $\text{sgn}(A)$ 是 SNS 的；
    2. 对每个 $j = 1, \ldots, n$，将 $A$ 的第 $j$ 列替换为 $b$ 后得到的矩阵 $A_j$ 的符号模式也是 SNS 的。

??? proof "证明"
    由 Cramer 法则，$x_j = \det(A_j) / \det(A)$。

    **充分性**：若 $\text{sgn}(A)$ 是 SNS 的，则 $\det(A)$ 的符号由 $\text{sgn}(A)$ 唯一确定。若 $\text{sgn}(A_j)$ 也是 SNS 的，则 $\det(A_j)$ 的符号由 $\text{sgn}(A_j)$ 唯一确定。因此 $x_j$ 的符号 $= \text{sgn}(\det(A_j)) \cdot \text{sgn}(\det(A))^{-1}$ 唯一确定。

    **必要性**：若 $\text{sgn}(A)$ 不是 SNS 的，则存在 $A' \in Q(\text{sgn}(A))$ 使得 $\det(A') = 0$，此时 $A'$ 奇异，$x'$ 不存在或不唯一。若某个 $\text{sgn}(A_j)$ 不是 SNS 的，则 $\det(A_j')$ 可以变号，从而 $x_j'$ 的符号不确定。$\blacksquare$

---

## 65A.4 定性稳定矩阵

### 符号稳定性

!!! definition "定义 65A.6 (定性稳定)"
    符号模式 $\mathcal{S}$ 称为**定性稳定**（qualitatively stable），若 $Q(\mathcal{S})$ 中的每个矩阵都是**稳定的**（即所有特征值具有负实部）。

!!! theorem "定理 65A.4 (定性稳定的充要条件)"
    $n \times n$ 符号模式 $\mathcal{S}$ 是定性稳定的，当且仅当以下条件全部成立：

    1. $s_{ii} \leq 0$ 对所有 $i$，且至少有一个 $s_{ii} < 0$（负自调节）；
    2. $s_{ij} s_{ji} \leq 0$ 对所有 $i \neq j$（无正反馈回路）；
    3. $D(\mathcal{S})$ 中没有长度 $\geq 3$ 的回路使得所有边权的乘积为正（无正循环）；
    4. $(-1)^k c_k \geq 0$ 对所有 $k$，其中 $c_k$ 是特征多项式的系数（Routh-Hurwitz 条件的符号版本）；
    5. $D(\mathcal{S})$ 满足适当的连通性条件。

??? proof "证明"
    **必要性方向的关键步骤**：

    **条件 1**：若某个 $s_{ii} > 0$，则可以取 $A$ 使 $a_{ii}$ 很大而其他元素很小，此时 $A$ 有正特征值（接近 $a_{ii}$），不稳定。若所有 $s_{ii} = 0$，则 $\operatorname{tr}(A) = 0$，但稳定矩阵要求 $\operatorname{tr}(A) = \sum \operatorname{Re}(\lambda_i) < 0$。因此至少一个对角元素必须为负。

    **条件 2**：设 $s_{ij} > 0$ 且 $s_{ji} > 0$。考虑 $2 \times 2$ 子矩阵 $\begin{pmatrix} s_{ii} & + \\ + & s_{jj} \end{pmatrix}$。取 $a_{ij}$ 和 $a_{ji}$ 很大，其他元素很小，子矩阵的特征值为
    $$
    \frac{a_{ii} + a_{jj}}{2} \pm \sqrt{\left(\frac{a_{ii} - a_{jj}}{2}\right)^2 + a_{ij}a_{ji}}
    $$
    当 $a_{ij}a_{ji}$ 足够大时，特征值之一可以为正。因此正反馈对 $(s_{ij}s_{ji} > 0)$ 破坏稳定性。

    **条件 3**：长度为 $k$ 的正循环 $i_1 \to i_2 \to \cdots \to i_k \to i_1$（所有 $a_{i_l, i_{l+1}} > 0$，且 $\text{sgn}(\sigma) \cdot \prod > 0$）贡献特征多项式中的正项，可能产生正实部特征值。

    **充分性**的证明需要仔细的归纳论证，利用 Routh-Hurwitz 判据的符号版本，验证在给定符号约束下特征多项式的系数满足稳定性条件。$\blacksquare$

### 生态学应用

!!! note "注"
    在数学生态学中，**群落矩阵**（community matrix）$A$ 描述了物种之间的交互强度。符号 $a_{ij}$ 表示物种 $j$ 对物种 $i$ 的影响方向：$a_{ij} > 0$（正面，如共生）、$a_{ij} < 0$（负面，如竞争或捕食）、$a_{ij} = 0$（无直接交互）。定性稳定性保证无论交互的具体强度如何，生态系统都是稳定的。

!!! example "例 65A.5"
    捕食-被捕食模型：
    $$
    \mathcal{S} = \begin{pmatrix} - & - \\ + & - \end{pmatrix}
    $$
    物种 1 被物种 2 捕食（$s_{12} < 0$），物种 2 以物种 1 为食（$s_{21} > 0$）。两者都有负自调节（$s_{11}, s_{22} < 0$）。

    验证定性稳定条件：$s_{12}s_{21} = (-)(+) < 0$（满足条件 2）。对于 $A \in Q(\mathcal{S})$：$\operatorname{tr}(A) = a_{11} + a_{22} < 0$，$\det(A) = a_{11}a_{22} - a_{12}a_{21} > 0$（因为 $a_{11}a_{22} > 0$ 和 $-a_{12}a_{21} > 0$）。特征多项式 $\lambda^2 - \operatorname{tr}(A)\lambda + \det(A)$ 的系数满足 Routh-Hurwitz 条件，因此稳定。

!!! example "例 65A.6"
    竞争模型：
    $$
    \mathcal{S}' = \begin{pmatrix} - & - \\ - & - \end{pmatrix}
    $$
    $s_{12}s_{21} = (-)(-)= + > 0$，**不满足**条件 2。例如 $A = \begin{pmatrix} -1 & -2 \\ -2 & -1 \end{pmatrix}$ 的特征值为 $-3, 1$，不稳定。

!!! example "例 65A.6b"
    三物种食物链模型：
    $$
    \mathcal{S} = \begin{pmatrix} - & - & 0 \\ + & - & - \\ 0 & + & - \end{pmatrix}
    $$
    物种 1 被物种 2 捕食，物种 2 被物种 3 捕食，物种 3 不直接影响物种 1。

    验证定性稳定条件：
    - 条件 1：$s_{ii} = -$ 对所有 $i$，满足；
    - 条件 2：$s_{12}s_{21} = (-)(+) = - \leq 0$，$s_{23}s_{32} = (-)(+) = - \leq 0$，$s_{13}s_{31} = 0 \leq 0$，满足；
    - 条件 3：唯一可能的长度 3 圈为 $1 \to 2 \to 3 \to 1$，但 $s_{31} = 0$，无此圈。满足；
    - 行列式条件：$\det(A) = -a_{11}a_{22}a_{33} - a_{12}a_{21}a_{33} - a_{11}a_{23}a_{32}$。由符号条件，$-a_{11}a_{22}a_{33} < 0$（因为三个负数之积为负，再取负号为正？让我们仔细分析）。实际上 $a_{11}, a_{22}, a_{33} < 0$，$a_{11}a_{22}a_{33} < 0$，$-a_{11}a_{22}a_{33} > 0$。$a_{12} < 0, a_{21} > 0$，$a_{12}a_{21} < 0$，$-a_{12}a_{21}a_{33} = -a_{12}a_{21} \cdot a_{33}$，$a_{33} < 0$，所以 $-a_{12}a_{21}a_{33} > 0$。类似地 $-a_{11}a_{23}a_{32} > 0$。因此 $\det(A) > 0$。$(-1)^3 \det(A) < 0$，满足 Hurwitz 条件。

    此模式是定性稳定的。

---

## 65A.5 完全不可分解矩阵

!!! definition "定义 65A.7 (完全不可分解矩阵)"
    $n \times n$ $(0,1)$-矩阵 $A$（或等价地，零-非零模式 $\mathcal{P}$）称为**完全不可分解的**（fully indecomposable），若不存在置换矩阵 $P, Q$ 使得
    $$
    PAQ = \begin{pmatrix} B & C \\ 0 & D \end{pmatrix}
    $$
    其中 $B$ 和 $D$ 是方阵（且 $1 \leq \dim B < n$）。

!!! theorem "定理 65A.5 (完全不可分解性的等价条件)"
    以下条件等价：

    1. 模式 $\mathcal{P}$ 完全不可分解；
    2. 二部图 $BG(\mathcal{P})$ 中每条边都属于某个完美匹配；
    3. 不存在行集 $R \subsetneq \{1,\ldots,n\}$ 使得 $\mathcal{P}$ 在 $R$ 行中非零元素仅出现在某个列集 $C$ 中，且 $|C| < |R|$。

??? proof "证明"
    **(1) $\Leftrightarrow$ (3)**：若存在 $R$ 和 $C$（$|C| \leq |R| - 1$）使得 $\mathcal{P}$ 的行 $R$ 中非零元素仅在列 $C$ 中，则通过行列置换可将 $\mathcal{P}$ 化为分块上三角形式（$R$ 对应行放上面，$C$ 对应列放左边），其中左下角为零块。这就是分块上三角分解，矛盾。

    **(2) $\Rightarrow$ (1)**：若 $\mathcal{P}$ 可分解为 $\begin{pmatrix} B & C \\ 0 & D \end{pmatrix}$，则零块对应的边不属于任何完美匹配（因为零块中没有可以跨越上下分块的匹配边）。$\blacksquare$

---

## 65A.6 近似可约矩阵

!!! definition "定义 65A.8 (近似可约矩阵)"
    $n \times n$ $(0,1)$-矩阵 $A$ 称为**近似可约的**（nearly reducible），若 $A$ 不可约，但将 $A$ 的任何一个非零元素改为零后变为可约的。

!!! theorem "定理 65A.6 (近似可约矩阵的性质)"
    若 $A$ 是近似可约的，则：

    1. $A$ 恰好有 $n$ 个非零元素（最少的不可约非零元素数）；
    2. $A$ 的有向图 $D(A)$ 是一个 Hamilton 圈加上可能的自环；
    3. $A$ 在置换相似下等价于某个"循环+扰动"形式。

??? proof "证明"
    **(1)** 不可约 $n \times n$ 矩阵的有向图至少有 $n$ 条边（因为强连通图至少 $n$ 条边，最少时为 Hamilton 圈）。近似可约性要求恰好 $n$ 条边（删除任何一条都破坏强连通性）。

    **(2)** 有 $n$ 条边的强连通有向图必须是 Hamilton 圈（可能加自环，但自环不贡献强连通性所需的边数，且删除自环不影响强连通性，所以如果有自环，非自环边必须也恰好构成 Hamilton 圈——但 Hamilton 圈已经有 $n$ 条边，因此没有额外自环的空间）。更准确地说，$n$ 个顶点的强连通图至少 $n$ 条边，且等号成立当且仅当图是 Hamilton 圈。$\blacksquare$

---

## 65A.7 潜在幂零模式

!!! definition "定义 65A.9 (潜在幂零模式)"
    零-非零模式 $\mathcal{P}$ 称为**潜在幂零的**（potentially nilpotent），若存在幂零矩阵 $A$（即 $A^n = 0$）使得 $\mathcal{Z}(A) = \mathcal{P}$。

!!! theorem "定理 65A.7 (潜在幂零的必要条件)"
    若 $\mathcal{P}$ 是潜在幂零的，则：

    1. $D(\mathcal{P})$ 中每个顶点至少有一个自环为零（即 $\mathcal{P}_{ii} = 0$ 对至少一个 $i$）——更强地，$\operatorname{tr}(A) = 0$ 要求对角线非零元素之和为零；
    2. 若 $\mathcal{P}$ 是不可约的，则 $\mathcal{P}$ 的对角线至少有一个零元素。

??? proof "证明"
    若 $A$ 幂零，则 $A$ 的所有特征值为零，$\operatorname{tr}(A) = 0$，$\operatorname{tr}(A^2) = 0$，等等。特别地 $\sum_i a_{ii} = 0$。若所有 $\mathcal{P}_{ii} = *$（非零），则 $a_{ii}$ 都非零且和为零，这对特定符号模式是可能的但施加了约束。

    对于不可约模式，更强的约束来自幂零性要求 $\det(A) = 0$，这对具有足够多非零元素的不可约模式施加了限制。$\blacksquare$

---

## 65A.8 谱任意模式与幂零-Jacobi 方法

!!! definition "定义 65A.10 (谱任意模式)"
    零-非零模式 $\mathcal{P}$ 称为**谱任意的**（spectrally arbitrary），若对任意首一 $n$ 次实多项式 $p(\lambda)$，存在矩阵 $A$ 使得 $\mathcal{Z}(A) = \mathcal{P}$ 且 $A$ 的特征多项式恰好是 $p(\lambda)$。

!!! theorem "定理 65A.8 (幂零-Jacobi 方法)"
    **幂零-Jacobi 方法**是证明模式 $\mathcal{P}$ 谱任意性的标准技术。若能找到具有模式 $\mathcal{P}$ 的幂零矩阵 $A_0$（即 $\mathcal{Z}(A_0) = \mathcal{P}$，$A_0$ 幂零），且 Jacobi 矩阵
    $$
    J = \frac{\partial(c_1, c_2, \ldots, c_n)}{\partial(\text{非零元素参数})}
    $$
    在 $A_0$ 处是非奇异的（其中 $c_k$ 是特征多项式的系数关于非零元素的函数），则 $\mathcal{P}$ 是谱任意的。

??? proof "证明"
    幂零矩阵 $A_0$ 的特征多项式为 $\lambda^n$（系数全为零，除首项）。非奇异 Jacobi 矩阵保证了隐函数定理的适用：在 $A_0$ 附近，通过微调 $\mathcal{P}$ 的非零元素，可以将特征多项式的系数调整为任意给定值。

    由隐函数定理，映射"非零元素 $\to$ 特征多项式系数"在 $A_0$ 的邻域内是满射的。由于特征多项式系数空间为 $\mathbb{R}^n$，且非零元素参数空间维数 $\geq n$（Jacobi 矩阵非奇异要求至少 $n$ 个自由参数），因此可以实现任意特征多项式。

    从邻域推广到全局需要额外的连续性论证（利用同伦方法或拓扑度论）。$\blacksquare$

!!! example "例 65A.7"
    模式
    $$
    \mathcal{P} = \begin{pmatrix} * & * & 0 \\ 0 & * & * \\ * & 0 & * \end{pmatrix}
    $$
    取幂零矩阵 $A_0 = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$——但 $\mathcal{Z}(A_0) \neq \mathcal{P}$（缺少 $(1,1), (2,2), (3,1), (3,3)$ 位置的非零元素）。需要寻找具有正确模式的幂零矩阵。

    实际上，确定谱任意模式是组合矩阵论的活跃研究领域。

---

## 65A.9 要求/允许问题

!!! definition "定义 65A.11 (要求/允许问题)"
    对于给定的矩阵性质 $\mathcal{P}$（如非奇异性、特定秩等）：

    - **要求问题**（require problem）：符号模式 $\mathcal{S}$ **要求**性质 $\mathcal{P}$，若 $Q(\mathcal{S})$ 中的**每个**矩阵都具有性质 $\mathcal{P}$。
    - **允许问题**（allow problem）：符号模式 $\mathcal{S}$ **允许**性质 $\mathcal{P}$，若 $Q(\mathcal{S})$ 中**存在**矩阵具有性质 $\mathcal{P}$。

!!! example "例 65A.8"
    - $\mathcal{S}$ 要求非奇异性 $\Leftrightarrow$ $\mathcal{S}$ 是 SNS 的。
    - $\mathcal{S}$ 允许稳定性 $\Leftrightarrow$ 存在 $A \in Q(\mathcal{S})$ 使得所有特征值具有负实部。
    - $\mathcal{S}$ 要求稳定性 $\Leftrightarrow$ $\mathcal{S}$ 是定性稳定的。

!!! theorem "定理 65A.9 (要求/允许的关系)"
    1. "要求"蕴含"允许"（若性质被所有矩阵满足，当然至少有一个满足）。
    2. 对于非奇异性：要求非奇异性（SNS）是一个组合条件（行列式展开项同号），而允许非奇异性要弱得多（只需 $Q(\mathcal{S})$ 中至少一个矩阵行列式非零）。
    3. 对于许多性质，"要求"和"允许"之间的差距很大，刻画它们需要不同的数学工具。

---

## 65A.10 L-矩阵与 SNS 矩阵的图论刻画

!!! definition "定义 65A.12 (L-矩阵)"
    符号模式 $\mathcal{S}$ 称为 **L-矩阵**，若 $\mathcal{S}$ 满足：$s_{ii} \neq 0$ 对所有 $i$，且 $s_{ij} s_{ji} \leq 0$ 对所有 $i \neq j$。L-矩阵在数值线性代数中对应 M-矩阵和 H-矩阵的符号结构。

!!! theorem "定理 65A.10 (SNS L-矩阵的图论刻画)"
    设 $\mathcal{S}$ 是 L-矩阵（$s_{ii} = -$ 对所有 $i$，$s_{ij} s_{ji} \leq 0$ 对 $i \neq j$）。则 $\mathcal{S}$ 是 SNS 的，当且仅当 $\mathcal{S}$ 的有向图 $D(\mathcal{S})$（去除自环后）没有偶数长度的圈，其中圈中所有边的符号乘积为正。

??? proof "证明"
    由 SNS 的定义，行列式展开中所有非零项必须同号。对于 L-矩阵，恒等置换 $\sigma = \text{id}$ 的贡献为 $\prod_i s_{ii} = (-1)^n$（所有对角元为负）。

    其他置换 $\sigma \neq \text{id}$ 可以分解为不相交圈的乘积。$\sigma$ 对行列式的贡献符号为 $\text{sgn}(\sigma) \cdot \prod_i s_{i,\sigma(i)}$。

    L-矩阵的条件 $s_{ij}s_{ji} \leq 0$ 保证了长度 2 的圈（转置对）的贡献符号与恒等置换一致。对于更长的圈，需要逐一检查符号。

    若存在偶数长度的正循环（圈中边符号的乘积为正，且 $\text{sgn}$ 因子使总贡献与恒等置换异号），则该置换的贡献与恒等置换异号，$\mathcal{S}$ 不是 SNS 的。反之，若不存在这样的循环，所有贡献同号。$\blacksquare$

---

## 65A.11 经济学与控制论中的应用

### 定性比较静态

!!! note "注"
    在经济学的**比较静态分析**中，均衡方程 $F(x, \alpha) = 0$ 关于参数 $\alpha$ 的变化如何影响均衡 $x^*$。若 Jacobi 矩阵 $\frac{\partial F}{\partial x}$ 的符号模式已知，定性矩阵分析可以确定 $\frac{\partial x^*}{\partial \alpha}$ 的符号，从而无需知道精确数值即可确定比较静态的方向。

!!! example "例 65A.9"
    考虑简单供需模型。设均衡条件为 $F(p, q, \alpha) = 0$，其中 $p$ 是价格，$q$ 是数量，$\alpha$ 是外生参数。Jacobi 矩阵的符号模式为
    $$
    \text{sgn}\left(\frac{\partial F}{\partial (p,q)}\right) = \begin{pmatrix} + & - \\ - & + \end{pmatrix}
    $$
    这是 SNS 的（$\det > 0$），因此均衡的定性比较静态结论不依赖于供需曲线的具体斜率。

### 结构可控性

!!! definition "定义 65A.13 (结构可控性)"
    线性系统 $(A, B)$ 称为**结构可控的**，若存在 $A' \in Q(\mathcal{Z}(A))$ 和 $B' \in Q(\mathcal{Z}(B))$，使得 $(A', B')$ 是可控的（即可控性矩阵 $[B', A'B', (A')^2B', \ldots]$ 满秩）。

!!! theorem "定理 65A.11 (结构可控性的图论条件)"
    系统 $(A, B)$ 是结构可控的，当且仅当以下两个条件成立：

    1. 有向图 $D(A, B)$ 中，从输入顶点到每个状态顶点都存在有向路径（可达性）；
    2. 不存在**膨胀**（dilation）：没有一组状态顶点的子集 $S$ 使得 $S$ 的后继集合严格小于 $|S|$（连同输入顶点的贡献）。

??? proof "证明"
    **充分性方向**：若图论条件满足，可以证明对"几乎所有"参数值（Lebesgue 测度意义下的稠密开集），$(A', B')$ 是可控的。这利用了可控性矩阵的行列式作为参数的多项式，非零多项式的零点集测度为零。

    **必要性方向**：若可达性不满足，则无论参数如何选择，不可达的状态顶点对应的状态变量无法被控制。若存在膨胀，则可控性矩阵的秩受到结构性限制。$\blacksquare$

---

## 65A.12 SNS 矩阵的组合计数

!!! theorem "定理 65A.12 (SNS 模式的计数)"
    $n \times n$ SNS 符号模式的数量随 $n$ 增长远小于所有符号模式的数量 $3^{n^2}$。具体地：

    1. 每个 SNS 模式的有向图中完美匹配的数量恰为 1（即行列式展开中恰有一个非零项），或多个非零项但全部同号；
    2. 最简单的 SNS 类是对角矩阵的符号模式（$n!$ 个广义对角模式）；
    3. 非平凡的 SNS 模式（行列式展开中多项非零）在 $n$ 增大时变得稀少。

!!! example "例 65A.10"
    $n = 2$ 时，所有 SNS 模式（非对角元允许为 $0$）：

    - $\begin{pmatrix} + & 0 \\ 0 & + \end{pmatrix}$（和所有对角元素非零的对角类型）；
    - $\begin{pmatrix} + & + \\ 0 & + \end{pmatrix}$ 及其变体（上三角类型）；
    - $\begin{pmatrix} + & + \\ - & + \end{pmatrix}$（两项 $a_{11}a_{22} > 0$ 和 $-a_{12}a_{21} > 0$ 同正）。

    $\begin{pmatrix} + & + \\ + & + \end{pmatrix}$ 不是 SNS（两项异号）。

---

## 65A.13 符号模式与图的着色

!!! theorem "定理 65A.13 (符号非奇异与二部图)"
    $n \times n$ 符号模式 $\mathcal{S}$ 是 SNS 的，当且仅当行列式展开中所有覆盖 $n$ 个顶点的**有向完美匹配**（置换）对应的贡献项全部同号。

    这可以用二部图的语言重新表述：构造加权二部图 $BG(\mathcal{S})$（行顶点 $\{r_1, \ldots, r_n\}$ 和列顶点 $\{c_1, \ldots, c_n\}$），边 $(r_i, c_j)$ 的权重为 $s_{ij} \in \{+, -, 0\}$（$0$ 表示无边）。完美匹配的权重为各边权重与置换符号的乘积。$\mathcal{S}$ 是 SNS 的当且仅当所有完美匹配的权重同号。

!!! example "例 65A.11"
    考虑 $\mathcal{S} = \begin{pmatrix} + & - \\ + & + \end{pmatrix}$。两个完美匹配：

    - $\sigma = \text{id}$：权重 $= (+1) \cdot (+)(+) = +$；
    - $\sigma = (12)$：权重 $= (-1) \cdot (-)(+) = +$（$\text{sgn}(12) = -1$，$s_{12} \cdot s_{21} = (-)(+)$，总计 $(-1) \cdot (-1) \cdot (+) = +$）。

    两项同号（均为 $+$），因此 $\mathcal{S}$ 是 SNS 的。

---

## 习题

!!! question "习题 65A.1"
    给出一个 $3 \times 3$ SNS 符号模式，并验证行列式展开中所有非零项同号。

!!! question "习题 65A.2"
    证明：$n \times n$ 上三角符号模式（对角线元素非零）总是 SNS 的。

!!! question "习题 65A.3"
    设 $\mathcal{S} = \begin{pmatrix} - & + \\ - & - \end{pmatrix}$。判断 $\mathcal{S}$ 是否定性稳定，并证明你的结论。

!!! question "习题 65A.4"
    证明：定性稳定矩阵的对角线元素必须全为负（不能为零）。

!!! question "习题 65A.5"
    给出一个 $3 \times 3$ 完全不可分解的零-非零模式，并验证二部图中每条边属于某个完美匹配。

!!! question "习题 65A.6"
    构造一个 $4 \times 4$ 近似可约矩阵，画出其有向图（Hamilton 圈），验证删除任一非零元素后矩阵变为可约。

!!! question "习题 65A.7"
    证明：若模式 $\mathcal{P}$ 的有向图 $D(\mathcal{P})$ 没有长度 $\geq 2$ 的圈，则 $\mathcal{P}$ 是潜在幂零的。

!!! question "习题 65A.8"
    设 $\mathcal{S} = \begin{pmatrix} - & - & 0 \\ + & - & - \\ 0 & + & - \end{pmatrix}$。验证这是一个定性稳定的符号模式（验证所有条件）。

!!! question "习题 65A.9"
    构造一个 $3 \times 3$ 符号模式，使其允许稳定性但不要求稳定性。即找到 $\mathcal{S}$ 使得 $Q(\mathcal{S})$ 中既有稳定矩阵也有不稳定矩阵。

!!! question "习题 65A.10"
    证明：若 $n \times n$ 符号模式 $\mathcal{S}$ 的有向图是无环的（DAG），则 $Q(\mathcal{S})$ 中每个矩阵都是幂零的。

!!! question "习题 65A.11"
    利用幂零-Jacobi 方法，证明 $2 \times 2$ 全非零模式 $\begin{pmatrix} * & * \\ * & * \end{pmatrix}$ 是谱任意的。

!!! question "习题 65A.12"
    设线性系统 $Ax = b$，$\mathcal{S} = \text{sgn}(A) = \begin{pmatrix} + & - \\ 0 & + \end{pmatrix}$，$\text{sgn}(b) = \begin{pmatrix} + \\ + \end{pmatrix}$。判断系统是否符号可解，并求解的符号。
