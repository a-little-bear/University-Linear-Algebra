# 第 65 章 组合矩阵论

<div class="context-flow" markdown>

**前置**：行列式(Ch3) · 图论(Ch27) · 特征值(Ch6)

**本章脉络**：零-非零模式 → 符号非奇异矩阵 → 定性矩阵理论 → 最小秩问题 → 零强迫数 → Hadamard 矩阵 → 竞赛矩阵 → 矩阵模式的允许性

**延伸**：组合矩阵论是离散数学与线性代数最深刻的交叉领域；最小秩问题与量子信息（图的最小秩与量子态的独立数）和通信复杂性有深刻联系

</div>

组合矩阵论研究矩阵的**组合结构**（如零-非零模式、符号模式、图结构）如何影响其**代数性质**（如秩、行列式、特征值）。这是离散数学与线性代数最自然的交叉点：一方面，矩阵的谱论为图论提供了强大的代数工具；另一方面，图的组合结构又对矩阵的代数性质施加了深刻的约束。

本章从矩阵的零-非零模式出发，经由符号非奇异矩阵和定性矩阵分析，到达最小秩问题和零强迫数，然后讨论 Hadamard 矩阵和竞赛矩阵这两类经典组合矩阵，最后探讨矩阵模式的允许性问题。

---

## 65.1 零-非零模式

<div class="context-flow" markdown>

**核心问题**：矩阵的零-非零位置分布能够决定哪些代数性质？

</div>

### 基本定义

!!! definition "定义 65.1 (零-非零模式)"
    矩阵 $A = (a_{ij}) \in M_{m \times n}(\mathbb{R})$ 的**零-非零模式**（zero-nonzero pattern）是一个 $\{0, *\}$ 矩阵 $\mathcal{Z}(A)$，定义为
    $$
    \mathcal{Z}(A)_{ij} = \begin{cases} * & \text{若 } a_{ij} \neq 0 \\ 0 & \text{若 } a_{ij} = 0 \end{cases}
    $$

    模式 $\mathcal{P}$ 的**定性类**（qualitative class）定义为
    $$
    Q(\mathcal{P}) = \{A \in M_{m \times n}(\mathbb{R}) : \mathcal{Z}(A) = \mathcal{P}\}
    $$
    即所有具有相同零-非零位置的矩阵的集合。

!!! example "例 65.1"
    模式 $\mathcal{P} = \begin{pmatrix} * & * \\ 0 & * \end{pmatrix}$ 的定性类为
    $$
    Q(\mathcal{P}) = \left\{\begin{pmatrix} a & b \\ 0 & d \end{pmatrix} : a, b, d \neq 0\right\}
    $$
    这个类中所有矩阵都是上三角且 $\det = ad \neq 0$，因此都是非奇异的。

### 模式与有向图

!!! definition "定义 65.2 (模式的有向图)"
    $n \times n$ 模式 $\mathcal{P}$ 的**有向图**（directed graph）$D(\mathcal{P})$ 定义为：顶点集 $V = \{1, \ldots, n\}$，有向边 $(i, j) \in E$ 当且仅当 $\mathcal{P}_{ij} = *$。

    对角线元素 $\mathcal{P}_{ii} = *$ 对应自环。

!!! note "注"
    矩阵的零-非零模式与其有向图之间的对应关系是组合矩阵论的基础。矩阵的许多代数性质可以用其有向图的组合性质来刻画。

### 模式能决定的性质

!!! theorem "定理 65.1 (模式不变的性质)"
    以下性质**不能**仅由零-非零模式决定：

    1. 矩阵的秩（模式相同但秩可以不同）；
    2. 矩阵的可逆性；
    3. 特征值的具体数值。

    以下性质**可以**仅由模式决定或约束：

    1. 秩的上界（由模式决定的最大秩）；
    2. 不可约性（等价于有向图 $D(\mathcal{P})$ 是强连通的）；
    3. 在附加符号条件下的非奇异性（见 65.2 节）。

!!! example "例 65.2"
    模式 $\mathcal{P} = \begin{pmatrix} * & * \\ * & * \end{pmatrix}$ 的定性类中包含秩 1 矩阵（如 $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$）和秩 2 矩阵（如 $\begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$）。因此零-非零模式不决定秩。

    但模式 $\mathcal{P}' = \begin{pmatrix} * & 0 \\ 0 & * \end{pmatrix}$ 的定性类中所有矩阵的秩都等于 2。

---

## 65.2 符号非奇异矩阵

<div class="context-flow" markdown>

**核心问题**：什么样的符号模式保证矩阵的行列式恒不为零？

</div>

### 符号模式

!!! definition "定义 65.3 (符号模式)"
    矩阵 $A = (a_{ij})$ 的**符号模式**（sign pattern）是一个 $\{+, -, 0\}$ 矩阵 $\text{sgn}(A)$，定义为
    $$
    \text{sgn}(A)_{ij} = \begin{cases} + & \text{若 } a_{ij} > 0 \\ - & \text{若 } a_{ij} < 0 \\ 0 & \text{若 } a_{ij} = 0 \end{cases}
    $$

    符号模式 $\mathcal{S}$ 的**符号类**（sign pattern class）为
    $$
    Q(\mathcal{S}) = \{A \in M_n(\mathbb{R}) : \text{sgn}(A) = \mathcal{S}\}
    $$

### 符号非奇异

!!! definition "定义 65.4 (符号非奇异矩阵)"
    符号模式 $\mathcal{S}$ 称为**符号非奇异**（sign nonsingular, SNS），若 $Q(\mathcal{S})$ 中的每个矩阵都是非奇异的。

    等价地，$\mathcal{S}$ 是 SNS 当且仅当 $\det(A) \neq 0$ 对所有 $A \in Q(\mathcal{S})$ 成立，即行列式的符号由 $\mathcal{S}$ 唯一确定。

!!! theorem "定理 65.2 (SNS 的刻画)"
    符号模式 $\mathcal{S}$ 是符号非奇异的，当且仅当行列式展开
    $$
    \det(A) = \sum_{\sigma \in S_n} \text{sgn}(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}
    $$
    中的所有非零项具有**相同的符号**。

    等价地，$\mathcal{S}$ 是 SNS 当且仅当在有向图 $D(\mathcal{S})$ 中，所有覆盖全部 $n$ 个顶点的置换（permutation covering）所对应的行列式项（考虑符号后）都同号。

??? proof "证明"
    行列式 $\det(A) = \sum_\sigma \text{sgn}(\sigma) \prod_i a_{i,\sigma(i)}$。每一项 $\text{sgn}(\sigma) \prod_i a_{i,\sigma(i)}$ 的符号由 $\sigma$ 的奇偶性和各 $a_{i,\sigma(i)}$ 的符号共同决定——后者由符号模式 $\mathcal{S}$ 决定。

    因此，对于给定的 $\mathcal{S}$，行列式展开中每个非零项的符号是确定的。$\det(A)$ 对所有 $A \in Q(\mathcal{S})$ 不为零，当且仅当所有非零项的符号一致（从而不会相消）。

!!! example "例 65.3"
    符号模式 $\mathcal{S} = \begin{pmatrix} + & + \\ 0 & + \end{pmatrix}$ 是 SNS 的，因为 $\det(A) = a_{11}a_{22} > 0$（唯一非零项）。

    符号模式 $\mathcal{S}' = \begin{pmatrix} + & + \\ + & + \end{pmatrix}$ 不是 SNS 的，因为 $\det(A) = a_{11}a_{22} - a_{12}a_{21}$，其中 $a_{11}a_{22} > 0$ 且 $a_{12}a_{21} > 0$，两项可能相消。

!!! example "例 65.4"
    考虑 $3 \times 3$ 符号模式
    $$
    \mathcal{S} = \begin{pmatrix} + & - & 0 \\ 0 & + & - \\ - & 0 & + \end{pmatrix}
    $$

    行列式展开的非零项：
    - 恒等置换 $\sigma = (1)(2)(3)$：$a_{11}a_{22}a_{33} > 0$，$\text{sgn}(\sigma) = +1$，贡献 $> 0$；
    - 循环置换 $\sigma = (1 \; 2 \; 3)$：$a_{13}a_{21}a_{32}$，但 $a_{13} = 0$，此项为零；
    - 循环置换 $\sigma = (1 \; 3 \; 2)$：$a_{12}a_{23}a_{31}$，$\text{sgn}(\sigma) = +1$，$a_{12} < 0, a_{23} < 0, a_{31} < 0$，乘积 $< 0$，贡献 $< 0$。

    两个非零项符号不同，所以 $\mathcal{S}$ 不是 SNS。

---

## 65.3 定性矩阵分析

<div class="context-flow" markdown>

**核心问题**：什么条件下，线性系统的解的符号仅由系数矩阵的符号决定？

</div>

### 符号可解系统

!!! definition "定义 65.5 (符号可解系统)"
    线性系统 $Ax = b$ 称为**符号可解的**（sign-solvable），若对所有 $A' \in Q(\text{sgn}(A))$ 和 $b' \in Q(\text{sgn}(b))$，解 $x' = (A')^{-1}b'$ 的符号模式与 $x = A^{-1}b$ 相同。

!!! theorem "定理 65.3 (符号可解性条件)"
    系统 $Ax = b$ 是符号可解的，当且仅当以下条件成立：

    1. $\text{sgn}(A)$ 是符号非奇异的（SNS）；
    2. 对每个 $j = 1, \ldots, n$，Cramer 法则中的矩阵 $A_j$（将 $A$ 的第 $j$ 列替换为 $b$）的符号模式也是 SNS 的。

### 符号稳定矩阵

!!! definition "定义 65.6 (定性稳定矩阵)"
    符号模式 $\mathcal{S}$ 称为**定性稳定**（qualitatively stable）或**符号稳定**（sign stable），若 $Q(\mathcal{S})$ 中的每个矩阵都是稳定的（即所有特征值具有负实部）。

!!! theorem "定理 65.4 (定性稳定的充要条件)"
    $n \times n$ 符号模式 $\mathcal{S}$ 是定性稳定的，当且仅当以下条件全部成立：

    1. $s_{ii} \leq 0$ 对所有 $i$（对角线元素非正），且至少有一个 $s_{ii} < 0$；
    2. $s_{ij} s_{ji} \leq 0$ 对所有 $i \neq j$（非对角线元素的积非正——无正反馈回路）；
    3. $\mathcal{S}$ 的有向图中没有长度 $\geq 3$ 的回路使得所有边的乘积为正；
    4. 行列式条件：$(-1)^n \det(A) > 0$ 对所有 $A \in Q(\mathcal{S})$；
    5. 有向图 $D(\mathcal{S})$ 是强连通的，或满足更一般的可约条件。

!!! note "注"
    定性稳定性在数学生态学中有重要应用。生态系统中物种之间的相互作用通常只知道定性关系（捕食、竞争、共生——对应正、负、零）。如果交互矩阵的符号模式是定性稳定的，那么无论具体的交互强度如何，系统都是稳定的。

!!! example "例 65.5"
    考虑两个物种的竞争模型：
    $$
    \mathcal{S} = \begin{pmatrix} - & - \\ - & - \end{pmatrix}
    $$

    验证定性稳定条件：
    - $s_{11} = s_{22} = -$（满足条件 1）；
    - $s_{12}s_{21} = (-)(-)= + > 0$（**不满足**条件 2）。

    因此这个模式**不是**定性稳定的。实际上 $A = \begin{pmatrix} -1 & -2 \\ -2 & -1 \end{pmatrix}$ 的特征值为 $-3, 1$，有正特征值。

    而捕食-被捕食模型
    $$
    \mathcal{S}' = \begin{pmatrix} - & - \\ + & - \end{pmatrix}
    $$
    满足：$s_{12}s_{21} = (-)(+) = - < 0$（条件 2 满足）。这个模式是定性稳定的。

---

## 65.4 最小秩问题

<div class="context-flow" markdown>

**核心问题**：给定图 $G$，在所有满足零模式匹配 $G$ 的对称矩阵中，最小可能秩是多少？

</div>

### 图的最小秩

!!! definition "定义 65.7 (图的最小秩)"
    设 $G = (V, E)$ 是一个简单图，$|V| = n$。$G$ 的**最小秩**定义为
    $$
    \operatorname{mr}(G) = \min\{\operatorname{rank}(A) : A \in S_n(\mathbb{R}), \; a_{ij} \neq 0 \Leftrightarrow \{i,j\} \in E \; (i \neq j)\}
    $$

    注意对角线元素 $a_{ii}$ 不受约束（可以为任意实数）。

    **最大零度**（maximum nullity）定义为 $M(G) = n - \operatorname{mr}(G)$。

!!! note "注"
    最小秩问题将图的组合结构与矩阵的秩联系起来。直觉上，边越少的图，矩阵中非零元素越少，实现低秩的可能性越大——但精确关系远非显而易见。

!!! example "例 65.6"
    **完全图** $K_n$：$\operatorname{mr}(K_n) = 1$（取 $A = \mathbf{1}\mathbf{1}^T - I + I = \mathbf{1}\mathbf{1}^T$，秩为 1，但对角线元素为 1 而非 0；实际上取 $A = \mathbf{1}\mathbf{1}^T$，对角线为 1 亦可。更精确地，$A = J_n - (n-1)I$ 满足模式且秩可能不为 1。需要注意对角线不受约束。取 $a_{ij} = 1$ 对所有 $i, j$（包括对角线），则 $A = \mathbf{1}\mathbf{1}^T$ 秩为 1）。

    **路径图** $P_n$：$\operatorname{mr}(P_n) = n - 1$（三对角矩阵的最小秩）。

    **完全二部图** $K_{p,q}$：$\operatorname{mr}(K_{p,q}) = 2$。

### Colin de Verdiere 参数

!!! definition "定义 65.8 (Colin de Verdiere 参数)"
    图 $G$ 的 **Colin de Verdiere 参数** $\mu(G)$ 是满足以下条件的矩阵 $M \in S_n(\mathbb{R})$ 的第二小特征值的重数的最大值：

    1. $M$ 的零-非零模式与 $G$ 的邻接结构匹配（$m_{ij} < 0$ 若 $\{i,j\} \in E$，$m_{ij} = 0$ 若 $\{i,j\} \notin E$ 且 $i \neq j$）；
    2. $M$ 恰好有一个负特征值（零特征值）；
    3. 如果 $Mx = 0$，则不存在与 $M$ 模式相同的非零矩阵 $X$ 使得 $Xx = 0$（强 Arnold 假设）。

!!! theorem "定理 65.5 (Colin de Verdiere 参数的性质)"
    $\mu(G)$ 具有以下性质：

    1. $\mu(G) \leq 1$ 当且仅当 $G$ 是路径；
    2. $\mu(G) \leq 2$ 当且仅当 $G$ 是外平面图；
    3. $\mu(G) \leq 3$ 当且仅当 $G$ 是平面图；
    4. $\mu$ 在取子图运算下单调递减。

!!! note "注"
    Colin de Verdiere 参数是矩阵谱理论与拓扑图论之间最深刻的联系之一。平面性——一个纯粹拓扑的概念——竟然可以用矩阵的谱性质来完全刻画，这是组合矩阵论最优美的结果之一。

---

## 65.5 零强迫数

<div class="context-flow" markdown>

**核心问题**：如何用图上的组合过程来上界矩阵的最大零度？

</div>

### 零强迫过程

!!! definition "定义 65.9 (零强迫过程)"
    给定图 $G = (V, E)$，**零强迫过程**（zero forcing process）如下：

    1. 初始时，将 $V$ 中的某些顶点着为黑色，其余为白色；
    2. **颜色变换规则**（color change rule）：如果一个黑色顶点 $u$ 恰好有一个白色邻居 $v$，则 $v$ 变为黑色（称 $u$ 将 $v$ 强迫为黑色）；
    3. 重复应用颜色变换规则，直到不再有变化。

!!! definition "定义 65.10 (零强迫数)"
    图 $G$ 的**零强迫数**（zero forcing number）$Z(G)$ 是使得最终所有顶点变黑的初始黑色顶点集合的最小大小：
    $$
    Z(G) = \min\{|S| : S \subseteq V, \; S \text{ 通过零强迫过程可将所有顶点变黑}\}
    $$

!!! theorem "定理 65.6 (零强迫数与最大零度)"
    对任意图 $G$，
    $$
    M(G) \leq Z(G)
    $$
    即图的最大零度不超过其零强迫数。等价地，$\operatorname{mr}(G) \geq n - Z(G)$。

??? proof "证明"
    设 $A \in S_n(\mathbb{R})$ 的零模式匹配 $G$，$\operatorname{null}(A) = M(G)$。设 $\{x_1, \ldots, x_k\}$ 是 $\ker(A)$ 的一组基（$k = M(G)$）。

    取初始黑色集合 $S$，$|S| = n - k$，对应于核向量的一个"决定集"——即 $S$ 中的分量唯一确定核向量。

    关键观察：若 $Ax = 0$ 且 $x$ 在 $S$ 上的分量已知（全为零），则方程 $\sum_j a_{ij} x_j = 0$ 中，若顶点 $i \in S$（黑色）且仅有一个邻居 $v$ 不在 $S$ 中（白色），则 $a_{iv} x_v = -\sum_{j \neq v, j \sim i} a_{ij} x_j$，可以确定 $x_v$。这恰好是零强迫的颜色变换规则。

    因此 $S$ 是零强迫集，$n - M(G) \leq n - |S_{\min}|$，即 $M(G) \leq Z(G)$。

!!! example "例 65.7"
    **路径** $P_4$：顶点 $1 - 2 - 3 - 4$。

    初始着色 $\{1\}$（一个黑色顶点）：$1$ 的唯一白色邻居是 $2$，将 $2$ 变黑；$2$ 的唯一白色邻居是 $3$，将 $3$ 变黑；$3$ 的唯一白色邻居是 $4$，将 $4$ 变黑。全部变黑。

    因此 $Z(P_4) = 1$，$M(P_4) \leq 1$。实际上 $M(P_4) = 1$（三对角矩阵的零度至多为 1）。

!!! example "例 65.8"
    **圈** $C_n$：$Z(C_n) = 2$。需要初始着两个相邻顶点为黑色，然后依次沿圈传播。

    **完全图** $K_n$：$Z(K_n) = n - 1$。因为每个黑色顶点都有多个白色邻居，无法应用颜色变换规则，直到只剩一个白色顶点。

### 传播时间

!!! definition "定义 65.11 (传播时间)"
    给定零强迫集 $S$，其**传播时间**（propagation time）$\operatorname{pt}(G, S)$ 是零强迫过程完成所需的步数（每步并行应用所有可能的颜色变换）。

    图 $G$ 的**传播时间**为
    $$
    \operatorname{pt}(G) = \min_{|S| = Z(G)} \operatorname{pt}(G, S)
    $$

---

## 65.6 Hadamard 矩阵

<div class="context-flow" markdown>

**核心问题**：$n \times n$ 的 $\pm 1$ 矩阵何时能达到行列式的最大绝对值？

</div>

### 定义与基本性质

!!! definition "定义 65.12 (Hadamard 矩阵)"
    $n \times n$ 矩阵 $H$ 称为 **Hadamard 矩阵**，若 $H$ 的所有元素为 $\pm 1$，且
    $$
    H H^T = n I_n
    $$

    等价地，$H$ 的行（列）两两正交，每行（列）的 $\ell^2$ 范数为 $\sqrt{n}$。

!!! theorem "定理 65.7 (Hadamard 矩阵的必要条件)"
    若 $n \times n$ Hadamard 矩阵存在，则 $n = 1$，$n = 2$，或 $n \equiv 0 \pmod{4}$。

??? proof "证明"
    设 $H$ 是 $n \times n$ Hadamard 矩阵。不妨（通过行列置换和取负）设 $H$ 的前两行为
    $$
    h_1 = (1, 1, \ldots, 1), \quad h_2 = (\underbrace{1, \ldots, 1}_a, \underbrace{-1, \ldots, -1}_b)
    $$
    由正交性 $h_1 \cdot h_2 = a - b = 0$，故 $a = b = n/2$。

    设第三行为 $h_3$。设 $h_3$ 在前 $a$ 个位置中有 $p$ 个 $+1$ 和 $a-p$ 个 $-1$，在后 $b$ 个位置中有 $q$ 个 $+1$ 和 $b-q$ 个 $-1$。

    由 $h_1 \cdot h_3 = 0$：$p + q - (a-p) - (b-q) = 2p + 2q - n = 0$，故 $p + q = n/2$。
    由 $h_2 \cdot h_3 = 0$：$p - (a-p) - q + (b-q) = 2p - 2q = 0$，故 $p = q$。

    从而 $p = q = n/4$，要求 $n/4$ 为整数，即 $n \equiv 0 \pmod{4}$。

### Sylvester 构造

!!! theorem "定理 65.8 (Sylvester 构造)"
    对所有 $k \geq 0$，$2^k \times 2^k$ 的 Hadamard 矩阵存在。递归构造为
    $$
    H_1 = (1), \quad H_{2^k} = \begin{pmatrix} H_{2^{k-1}} & H_{2^{k-1}} \\ H_{2^{k-1}} & -H_{2^{k-1}} \end{pmatrix}
    $$

    即 $H_{2^k} = H_2 \otimes H_{2^{k-1}}$（Kronecker 积），其中 $H_2 = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$。

!!! example "例 65.9"
    $$
    H_4 = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & -1 & 1 & -1 \\ 1 & 1 & -1 & -1 \\ 1 & -1 & -1 & 1 \end{pmatrix}
    $$

    验证：$H_4 H_4^T = 4I_4$。这就是 Walsh-Hadamard 矩阵。

### Hadamard 猜想

!!! definition "定义 65.13 (Hadamard 猜想)"
    **Hadamard 猜想**断言：对所有 $n \equiv 0 \pmod{4}$ 的正整数 $n$，$n \times n$ Hadamard 矩阵存在。

!!! note "注"
    截至目前，Hadamard 猜想仍未解决。最小的未知阶数为 $n = 668$。已知的构造方法包括 Sylvester 构造、Paley 构造（利用有限域上的二次剩余）、以及各种乘法构造。

### 应用

!!! theorem "定理 65.9 (Hadamard 界)"
    对任意 $n \times n$ 实矩阵 $A = (a_{ij})$，若 $|a_{ij}| \leq 1$，则
    $$
    |\det(A)| \leq n^{n/2}
    $$
    等号成立当且仅当 $A$ 是 Hadamard 矩阵（元素为 $\pm 1$ 且行正交）。

!!! note "注"
    Hadamard 矩阵在信号处理（Walsh-Hadamard 变换）、纠错编码（Reed-Muller 码）、压缩感知（Hadamard 测量矩阵）和实验设计（正交阵列）中有广泛应用。

---

## 65.7 竞赛矩阵

<div class="context-flow" markdown>

**核心问题**：完全有向图（竞赛图）的邻接矩阵有什么特殊的代数性质？

</div>

### 定义

!!! definition "定义 65.14 (竞赛矩阵)"
    $n$ 阶**竞赛**（tournament）是完全图 $K_n$ 的一个定向，即对每对顶点 $\{i, j\}$，恰有一条有向边 $i \to j$ 或 $j \to i$。

    竞赛的**竞赛矩阵**（tournament matrix）$T = (t_{ij})$ 定义为
    $$
    t_{ij} = \begin{cases} 1 & \text{若 } i \to j \\ 0 & \text{若 } j \to i \text{ 或 } i = j \end{cases}
    $$

    注意 $T + T^T = J - I$（其中 $J$ 是全 1 矩阵），即 $t_{ij} + t_{ji} = 1$ 对 $i \neq j$。

### 分数得分

!!! definition "定义 65.15 (得分序列)"
    竞赛中顶点 $i$ 的**得分**（score）是 $s_i = \sum_{j=1}^n t_{ij}$，即 $i$ 击败的对手数。

    **得分序列**是将 $(s_1, \ldots, s_n)$ 排成非递减序后的序列。

!!! theorem "定理 65.10 (Landau 定理)"
    非负整数序列 $(s_1 \leq s_2 \leq \cdots \leq s_n)$ 是某个竞赛的得分序列，当且仅当
    $$
    \sum_{i=1}^k s_i \geq \binom{k}{2}, \quad k = 1, \ldots, n
    $$
    且 $\sum_{i=1}^n s_i = \binom{n}{2}$。

### 谱性质

!!! theorem "定理 65.11 (竞赛矩阵的特征值)"
    设 $T$ 是 $n$ 阶竞赛矩阵。则：

    1. $T$ 的特征值之和为 $\operatorname{tr}(T) = 0$；
    2. 特征值实部之和为 $\sum \operatorname{Re}(\lambda_i) = 0$；
    3. $\frac{n-1}{2}$ 是 $T$ 的一个特征值（对应特征向量 $\mathbf{1}$），因为 $(T + T^T)\mathbf{1} = (J - I)\mathbf{1} = (n-1)\mathbf{1}$，但需更精细的分析；
    4. 所有特征值的实部满足 $\operatorname{Re}(\lambda_i) \leq \frac{n-1}{2}$。

!!! example "例 65.10"
    $n = 3$ 的循环竞赛 $1 \to 2 \to 3 \to 1$：
    $$
    T = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}
    $$
    这是循环置换矩阵，特征值为 $1, \omega, \omega^2$（$\omega = e^{2\pi i/3}$）。得分序列为 $(1, 1, 1)$——完全平衡的竞赛。

---

## 65.8 矩阵模式的允许性

<div class="context-flow" markdown>

**核心问题**：给定零-非零模式，哪些特征值重数组合是允许的？

</div>

### 允许的特征值重数

!!! definition "定义 65.16 (允许的重数列表)"
    设 $\mathcal{P}$ 是 $n \times n$ 零-非零模式。正整数列表 $(m_1, m_2, \ldots, m_r)$（$\sum m_i = n$）称为被 $\mathcal{P}$ **允许的**（allowed），若存在矩阵 $A \in Q(\mathcal{P}) \cap S_n(\mathbb{R})$ 使得 $A$ 恰有 $r$ 个不同的特征值，重数分别为 $m_1, \ldots, m_r$。

!!! example "例 65.11"
    模式 $\mathcal{P} = \begin{pmatrix} * & * & 0 \\ * & * & * \\ 0 & * & * \end{pmatrix}$（路径 $P_3$ 的模式）。

    允许的重数列表：$(1, 1, 1)$（三个不同特征值——典型情形）。

    $(3)$ 即三重特征值：需要 $A = \lambda I$ 的模式，但 $\lambda I$ 的零模式要求非对角元素为零，与 $\mathcal{P}$ 矛盾。因此 $(3)$ 不被允许。

    $(1, 2)$：需要某个特征值有重数 2。对于三对角矩阵，这是可能的（选择特殊的对角和非对角元素值）。

### 逆特征值问题

!!! definition "定义 65.17 (图的逆特征值问题)"
    给定图 $G$，确定哪些**多重集**的实数可以作为某个匹配 $G$ 的零模式的对称矩阵的特征值集合。这称为图 $G$ 的**逆特征值问题**（inverse eigenvalue problem, IEP）。

!!! theorem "定理 65.12 (路径的 IEP)"
    对于路径图 $P_n$，任意 $n$ 个不同的实数 $\lambda_1 < \lambda_2 < \cdots < \lambda_n$ 都可以作为某个匹配 $P_n$ 模式的对称三对角矩阵的特征值。但不允许重复特征值。

!!! theorem "定理 65.13 (完全图的 IEP)"
    对于完全图 $K_n$，任意 $n$ 个实数的多重集（允许重复）都可以作为某个密集对称矩阵的特征值。因为 $K_n$ 对非对角元素没有零约束。

### 禁止子模式

!!! definition "定义 65.18 (禁止子模式)"
    若重数列表 $\ell$ 不被模式 $\mathcal{P}$ 允许，且存在"最小"的子模式 $\mathcal{P}'$ 使得 $\ell$ 不被包含 $\mathcal{P}'$ 的任何模式允许，则 $\mathcal{P}'$ 称为关于 $\ell$ 的**禁止子模式**。

!!! note "注"
    矩阵模式的允许性问题与图的最小秩问题密切相关。事实上，$M(G) = k$ 等价于重数列表 $(k, 1, \ldots, 1)$ 被 $G$ 的模式允许。完整地描述所有图的允许重数列表是组合矩阵论中最困难的开放问题之一。

---

## 本章小结

本章研究了组合结构与矩阵代数性质之间的深刻联系。主要内容包括：

1. **零-非零模式**定义了矩阵的组合骨架。模式决定了有向图结构，但通常不能单独决定秩和特征值。

2. **符号非奇异矩阵**的行列式符号完全由符号模式决定，等价于行列式展开中所有非零项同号。

3. **定性矩阵分析**研究仅凭符号信息可以确定的系统性质，在生态学和经济学中有重要应用。

4. **最小秩问题**将图的组合结构与矩阵秩联系起来。Colin de Verdiere 参数更是将矩阵谱与图的拓扑性质（平面性）建立了等价。

5. **零强迫数**提供了最大零度的组合上界，零强迫过程是一个简洁而有力的组合工具。

6. **Hadamard 矩阵**是达到行列式上界的 $\pm 1$ 矩阵，其存在性（Hadamard 猜想）是组合设计理论中最著名的开放问题之一。

7. **竞赛矩阵**编码了完全有向图的结构，其谱性质揭示了竞赛的代数本质。

8. **模式的允许性**研究给定零模式下可以实现的特征值重数组合，是逆特征值问题的核心。

---

## 习题

!!! question "习题 65.1"
    给出一个 $3 \times 3$ 符号非奇异矩阵的例子，并验证行列式展开中所有非零项同号。

!!! question "习题 65.2"
    证明路径图 $P_n$ 的零强迫数 $Z(P_n) = 1$。

!!! question "习题 65.3"
    构造 $8 \times 8$ Hadamard 矩阵（使用 Sylvester 方法），并验证 $HH^T = 8I$。

!!! question "习题 65.4"
    证明 $n$ 阶竞赛矩阵 $T$ 满足 $T + T^T = J_n - I_n$。

!!! question "习题 65.5"
    计算圈图 $C_5$ 的最小秩 $\operatorname{mr}(C_5)$ 和零强迫数 $Z(C_5)$。

!!! question "习题 65.6"
    证明：若 $H$ 是 $n \times n$ Hadamard 矩阵，则 $|\det(H)| = n^{n/2}$。

!!! question "习题 65.7"
    证明完全图 $K_n$ 的零强迫数为 $n - 1$。

!!! question "习题 65.8"
    设 $T$ 是 $4$ 阶竞赛矩阵，得分序列为 $(0, 1, 2, 3)$。写出 $T$ 并计算其特征值。

!!! question "习题 65.9"
    证明：定性稳定矩阵的对角线元素必须全为负（不能为零）。

!!! question "习题 65.10"
    对于星图 $K_{1,n-1}$（一个中心顶点连接 $n-1$ 个叶子），求 $\operatorname{mr}(K_{1,n-1})$ 和 $Z(K_{1,n-1})$。
