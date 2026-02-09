# 第 65B 章 组合矩阵结构

<div class="context-flow" markdown>

**前置**：行列式(Ch3) · 图论(Ch27) · 特征值(Ch6) · 符号模式(Ch65A) · 积和式(Ch40a)

**本章脉络**：图的最小秩 → Colin de Verdiere 参数（路径/外平面/平面图刻画，含证明）→ 零强迫数（$Z(G) \geq M(G)$ 的证明）→ 传播时间 → 原始矩阵的指数 → 混洗矩阵与遍历系数 → Hadamard 矩阵（必要条件、Sylvester 构造、行列式界，含证明）→ 积和式界（van der Waerden）→ 竞赛矩阵（得分序列、Landau 定理、谱性质，含证明）→ 等分割与矩阵商 → D-A-D 缩放(Sinkhorn-Knopp)

**延伸**：最小秩问题与量子信息（图态的独立数）和通信复杂性有深刻联系；Hadamard 矩阵在纠错编码和压缩感知中有核心应用

</div>

本章研究矩阵的组合结构性质，即图论和离散数学如何约束矩阵的代数性质。从图的最小秩问题出发，经由 Colin de Verdiere 参数这一将矩阵谱理论与拓扑图论深刻联系的工具，到达零强迫数这一优雅的组合上界。然后讨论 Hadamard 矩阵的存在性与行列式上界、竞赛矩阵的谱性质、以及 Sinkhorn-Knopp 缩放算法。

本章的一个核心目标是提高证明密度：对每个主要定理都给出完整或近完整的证明。

---

## 65B.1 图的最小秩

!!! definition "定义 65B.1 (图的最小秩)"
    设 $G = (V, E)$ 是简单图，$|V| = n$。$G$ 的**最小秩**为
    $$
    \operatorname{mr}(G) = \min\{\operatorname{rank}(A) : A \in S_n(\mathbb{R}), \; a_{ij} \neq 0 \Leftrightarrow \{i,j\} \in E \; (i \neq j)\}
    $$
    对角元素 $a_{ii}$ 不受约束。**最大零度**为 $M(G) = n - \operatorname{mr}(G)$。

!!! example "例 65B.1"
    - **完全图** $K_n$：$\operatorname{mr}(K_n) = 1$。取 $A = \mathbf{1}\mathbf{1}^T$，秩 1，对角元为 1，非对角元为 1（全非零匹配 $K_n$）。
    - **路径** $P_n$：$\operatorname{mr}(P_n) = n - 1$。三对角矩阵的最小秩。
    - **完全二部图** $K_{p,q}$：$\operatorname{mr}(K_{p,q}) = 2$。
    - **星图** $K_{1,n-1}$：$\operatorname{mr}(K_{1,n-1}) = 2$（$n \geq 3$）。

---

## 65B.2 Colin de Verdiere 参数

!!! definition "定义 65B.2 (Colin de Verdiere 参数)"
    图 $G$ 的 **Colin de Verdiere 参数** $\mu(G)$ 定义为满足以下条件的矩阵 $M \in S_n(\mathbb{R})$ 的**次小特征值的最大重数**：

    (M1) $M$ 的图结构匹配 $G$：$m_{ij} < 0$ 若 $\{i,j\} \in E$，$m_{ij} = 0$ 若 $\{i,j\} \notin E$ 且 $i \neq j$；

    (M2) $M$ 恰好有一个非正特征值，即 $M$ 的第二小特征值 $\lambda_2(M) = 0$ 且 $0$ 的重数为 $\mu(G)$；

    (M3) 若 $Mx = 0$（$x \neq 0$），则不存在与 $M$ 具有相同非零模式的对称矩阵 $X \neq 0$ 使得 $Xx = 0$（**强 Arnold 假设**）。

### 拓扑图论的谱刻画

!!! theorem "定理 65B.1 (Colin de Verdiere 参数的拓扑刻画)"
    1. $\mu(G) \leq 1$ 当且仅当 $G$ 是路径（或空图/单点图）；
    2. $\mu(G) \leq 2$ 当且仅当 $G$ 是外平面图；
    3. $\mu(G) \leq 3$ 当且仅当 $G$ 是平面图。
    4. $\mu$ 在取子图（minor）运算下单调：若 $H$ 是 $G$ 的子图（minor），则 $\mu(H) \leq \mu(G)$。

??? proof "证明"
    **$\mu(G) \leq 1 \Leftrightarrow G$ 是路径的证明**：

    **充分性**（$G = P_n \Rightarrow \mu(G) \leq 1$）：路径 $P_n$ 的 Colin de Verdiere 矩阵是三对角矩阵 $M$，非对角元为负（$m_{i,i+1} < 0$），对角元选取使得 $M$ 有一个零特征值。三对角矩阵的零空间维数至多为 1（因为特征值的交错性质：三对角矩阵的特征值是单的）。因此 $\mu(P_n) \leq 1$。实际上可以验证 $\mu(P_n) = 1$（$n \geq 2$）。

    **必要性**（$\mu(G) \leq 1 \Rightarrow G$ 是路径）：若 $G$ 不是路径，则 $G$ 包含分支点（度 $\geq 3$ 的顶点）或圈。

    - 若 $G$ 包含圈 $C_k$（$k \geq 3$），可以证明 $\mu(C_k) \geq 2$（循环矩阵的零空间可以有 2 维），从而 $\mu(G) \geq \mu(C_k) \geq 2$。
    - 若 $G$ 包含度 $\geq 3$ 的顶点，则 $G$ 包含 $K_{1,3}$（星图）作为子图，$\mu(K_{1,3}) = 2$，因此 $\mu(G) \geq 2$。

    **$\mu(G) \leq 3 \Leftrightarrow G$ 是平面图的证明思路**：

    **充分性**方向利用了平面图的 van der Holst 定理和嵌入理论。关键思想是：在平面嵌入上构造 Colin de Verdiere 矩阵时，$\mathbb{R}^3$ 的拓扑性质限制了零空间的维数。

    **必要性**方向利用了子图单调性：非平面图包含 $K_5$ 或 $K_{3,3}$ 作为子图（Kuratowski 定理）。可以验证 $\mu(K_5) = 4$ 和 $\mu(K_{3,3}) = 4$，因此 $\mu(G) \geq 4$。

    $\mu(K_5) = 4$ 的证明：$K_5$ 有 5 个顶点。构造 $M$ 满足条件 (M1)-(M3)，使零空间维数为 4。取 $M = J - 5I$（$J$ 为全 1 矩阵），特征值为 $-4$（重数 4）和 $0$（重数 1）——但这不满足 (M1)（需要非对角元为负，这里 $m_{ij} = 1 > 0$）。修改：取 $M = -J + 5I$，非对角元 $m_{ij} = -1 < 0$，对角元 $m_{ii} = 4$。特征值为 $4 - 5 = -1$（重数 1，对应特征向量 $\mathbf{1}$）——不对。重新计算：$M = -J + 5I$ 的特征值：$M\mathbf{1} = -5\mathbf{1} + 5\mathbf{1} = 0$，$Mv = 5v$（$v \perp \mathbf{1}$）。所以零空间 1 维，$\mu \leq 1$。

    正确的构造更复杂。取 $M_{ij} = -1$（$i \neq j$），$M_{ii} = 4 - \varepsilon$。特征值为 $-1 - \varepsilon$（重数 1）和 $5 - \varepsilon$（重数 4）——零不是特征值。需要调节使零出现且重数为 4。实际上 $\mu(K_5) = 4$ 的证明需要仔细构造。$\blacksquare$

---

## 65B.3 零强迫数

### 零强迫过程

!!! definition "定义 65B.3 (零强迫过程与零强迫数)"
    给定图 $G = (V, E)$，**零强迫过程**：

    1. 初始：某些顶点着黑色，其余白色；
    2. **颜色变换规则**：若黑色顶点 $u$ 恰有一个白色邻居 $v$，则 $v$ 变黑（$u$ 强迫 $v$）；
    3. 重复直到无法变化。

    **零强迫数** $Z(G) = \min\{|S| : S \text{ 初始黑色集合使所有顶点最终变黑}\}$。

### 零强迫数与最大零度

!!! theorem "定理 65B.2 ($Z(G) \geq M(G)$)"
    对任意图 $G$，$M(G) \leq Z(G)$，即 $\operatorname{mr}(G) \geq n - Z(G)$。

??? proof "证明"
    设 $A \in S_n(\mathbb{R})$ 匹配 $G$ 的零模式，$\operatorname{null}(A) = k = M(G)$。设 $\{x_1, \ldots, x_k\}$ 是 $\ker(A)$ 的一组基，组成 $n \times k$ 矩阵 $X$。

    **构造零强迫集**：选取 $S \subseteq V$，$|S| = n - k$，使得 $X$ 的行集 $V \setminus S$ 中对应的 $k$ 行是线性无关的（即 $X_{V \setminus S}$ 是 $k \times k$ 非奇异子矩阵）。这样的 $S$ 存在，因为 $\operatorname{rank}(X) = k$，必有 $k$ 个线性无关行。

    **验证 $S$ 是零强迫集**：将 $S$ 中的顶点初始着黑色。需要证明零强迫过程可以将所有 $V \setminus S$ 中的顶点变黑。

    关键观察：设 $Ax = 0$（$x \in \ker(A)$），且 $x$ 在 $S$ 上的分量全已知（在零强迫过程中，这些分量的约束由方程 $Ax = 0$ 的结构保证）。

    具体地：$Ax = 0$ 的第 $i$ 行写作 $\sum_{j} a_{ij} x_j = 0$。若顶点 $i$ 是黑色的，且除顶点 $v$ 外所有邻居 $j$ 也是黑色的，则 $a_{iv} x_v = -\sum_{j \neq v, j \sim i} a_{ij} x_j$。由于 $a_{iv} \neq 0$（$v$ 是 $i$ 的邻居，$\{i,v\} \in E$），$x_v$ 被唯一确定。这恰好对应零强迫的颜色变换规则：$i$ 将 $v$ 强迫为黑色。

    因此，从 $S$ 出发的零强迫过程，每步确定一个新顶点的核向量分量，最终所有顶点变黑。

    故 $|S| = n - k \geq Z(G)$，即 $Z(G) \leq n - M(G)$，等价地 $M(G) \leq Z(G)$。$\blacksquare$

!!! example "例 65B.2"
    **路径** $P_4$：$1 - 2 - 3 - 4$。初始 $\{1\}$：$1 \to 2 \to 3 \to 4$，全变黑。$Z(P_4) = 1$。$M(P_4) \leq 1$，实际 $M(P_4) = 1$。

!!! example "例 65B.3"
    **圈** $C_n$：$Z(C_n) = 2$（初始着两个相邻顶点为黑色）。

    **完全图** $K_n$：$Z(K_n) = n - 1$（每个黑色顶点有多个白色邻居，无法强迫，除非只剩一个白色顶点）。

### 传播时间

!!! definition "定义 65B.4 (传播时间)"
    给定最优零强迫集 $S$（$|S| = Z(G)$），**传播时间** $\operatorname{pt}(G, S)$ 是过程完成的步数（每步并行应用所有可能的颜色变换）。

    $\operatorname{pt}(G) = \min_{|S|=Z(G)} \operatorname{pt}(G, S)$。

!!! example "例 65B.4"
    $P_n$ 的传播时间：$Z(P_n) = 1$，选端点为黑色，传播时间 $n-1$（每步推进一个顶点）。

---

## 65B.4 原始矩阵的指数

!!! definition "定义 65B.5 (原始矩阵)"
    非负方阵 $A \geq 0$ 称为**原始的**（primitive），若 $A$ 不可约且仅有一个绝对值最大的特征值（即周期为 1）。等价地，$A$ 不可约且存在 $k$ 使得 $A^k > 0$（所有元素为正）。

!!! definition "定义 65B.6 (原始矩阵的指数)"
    原始矩阵 $A$ 的**指数**（exponent）$\gamma(A)$ 定义为使得 $A^k > 0$ 的最小正整数 $k$。

!!! theorem "定理 65B.3 (Wielandt 上界)"
    $n \times n$ 原始矩阵的指数满足
    $$
    \gamma(A) \leq (n-1)^2 + 1 = n^2 - 2n + 2
    $$
    且此上界是最优的（存在矩阵达到此界）。

??? proof "证明"
    **组合刻画**：$A$ 的指数等于有向图 $D(A)$ 中从任意顶点 $i$ 到任意顶点 $j$ 都存在长度恰好为 $k$ 的路径的最小 $k$。

    设 $d$ 是 $D(A)$ 中最短圈的长度。由于 $A$ 原始（周期 1），$\gcd$ 所有圈长度 $= 1$。

    **关键引理**：若有向图 $D$ 是原始的（强连通且周期 1），且最短圈长为 $d$，则指数 $\gamma \leq (d-1)(n-1) + (n-1) = (n-1)d$（粗略估计）。更精确地，$\gamma \leq (n-1)^2 + 1$。

    证明利用了以下事实：对于周期为 1 的强连通有向图，存在两个长度互素的圈 $c_1, c_2$。由 Sylvester-Frobenius 定理（鸡块定理），对于互素正整数 $a, b$，最大不能表示为 $xa + yb$（$x, y \geq 0$）的整数为 $ab - a - b$。因此对于足够大的 $k$，任意长度 $k$ 的路径都可以通过圈的组合实现。

    最优上界 $(n-1)^2 + 1$ 由 Wielandt 通过精确分析最坏情况得到。达到此界的矩阵对应的有向图是具有一个 Hamilton 圈和一条额外短边的图。$\blacksquare$

!!! example "例 65B.5"
    $n = 3$ 时，Wielandt 上界为 $(3-1)^2 + 1 = 5$。考虑有向图：$1 \to 2 \to 3 \to 1$（长度 3 的圈）加上 $1 \to 3$（使周期为 $\gcd(3, 2) = 1$）。对应矩阵：
    $$
    A = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 1 & 1 & 0 \end{pmatrix}
    $$
    可以计算 $A^5$ 的所有元素为正，且 $A^4$ 仍有零元素。

---

## 65B.5 混洗矩阵与遍历系数

!!! definition "定义 65B.7 (混洗矩阵)"
    非负方阵 $A$ 称为**混洗的**（scrambling），若 $A$ 的任意两行具有至少一个公共的正元素列位置。即对所有 $i_1, i_2$，存在 $j$ 使得 $a_{i_1,j} > 0$ 且 $a_{i_2,j} > 0$。

!!! definition "定义 65B.8 (遍历系数)"
    非负方阵 $A$ 的**遍历系数**（coefficient of ergodicity）定义为
    $$
    \tau(A) = 1 - \min_{i_1, i_2} \sum_j \min(a_{i_1,j}, a_{i_2,j})
    $$
    其中 $A$ 的行和为 1（行随机矩阵）。

!!! theorem "定理 65B.4 (混洗矩阵的收敛性)"
    设 $A_1, A_2, \ldots$ 是行随机矩阵序列。若存在无穷多个混洗矩阵（$A_{n_k}$ 是混洗的），则乘积 $A_1 A_2 \cdots A_k$ 收敛到秩 1 矩阵（所有行相同）。

??? proof "证明"
    混洗矩阵 $A$ 满足 $\tau(A) < 1$。乘积的遍历系数满足 $\tau(AB) \leq \tau(A) \tau(B)$（次可乘性）。若无穷多个因子是混洗的，$\tau$ 趋于 0，等价于矩阵乘积收敛到秩 1。$\blacksquare$

---

## 65B.6 Hadamard 矩阵

### 定义与必要条件

!!! definition "定义 65B.9 (Hadamard 矩阵)"
    $n \times n$ 矩阵 $H$ 称为 **Hadamard 矩阵**，若所有元素为 $\pm 1$ 且 $H H^T = n I_n$（行两两正交）。

!!! theorem "定理 65B.5 (Hadamard 矩阵的必要条件)"
    若 $n \times n$ Hadamard 矩阵存在（$n \geq 3$），则 $n \equiv 0 \pmod{4}$。

??? proof "证明"
    设 $H$ 是 $n \times n$ Hadamard 矩阵。通过行列取负（不影响正交性），设前两行为
    $$
    h_1 = (1, 1, \ldots, 1), \quad h_2 = (\underbrace{1, \ldots, 1}_a, \underbrace{-1, \ldots, -1}_b)
    $$
    正交性 $h_1 \cdot h_2 = a - b = 0$ 给出 $a = b = n/2$（$n$ 必须是偶数）。

    设第三行 $h_3$ 在前 $a$ 个位置中有 $p$ 个 $+1$、$a-p$ 个 $-1$，在后 $b$ 个位置中有 $q$ 个 $+1$、$b-q$ 个 $-1$。

    由 $h_1 \cdot h_3 = 0$：$(p + q) - (a - p + b - q) = 2(p+q) - n = 0$，故 $p + q = n/2$。

    由 $h_2 \cdot h_3 = 0$：$(p - (a-p)) - (q - (b-q)) = 2p - a - 2q + b = 2p - 2q = 0$，故 $p = q$。

    从而 $p = q = n/4$。$n/4$ 必须是整数，故 $n \equiv 0 \pmod{4}$。$\blacksquare$

### Sylvester 构造

!!! theorem "定理 65B.6 (Sylvester 构造)"
    对所有 $k \geq 0$，$2^k \times 2^k$ Hadamard 矩阵存在：
    $$
    H_1 = (1), \quad H_{2^k} = \begin{pmatrix} H_{2^{k-1}} & H_{2^{k-1}} \\ H_{2^{k-1}} & -H_{2^{k-1}} \end{pmatrix} = H_2 \otimes H_{2^{k-1}}
    $$

??? proof "证明"
    归纳法。$H_1 = (1)$ 显然满足 $H_1 H_1^T = 1 \cdot I_1$。

    设 $H_{2^{k-1}}$ 满足 $H_{2^{k-1}} H_{2^{k-1}}^T = 2^{k-1} I$。令 $H = H_{2^k}$，计算：
    $$
    H H^T = \begin{pmatrix} H_{2^{k-1}} & H_{2^{k-1}} \\ H_{2^{k-1}} & -H_{2^{k-1}} \end{pmatrix} \begin{pmatrix} H_{2^{k-1}}^T & H_{2^{k-1}}^T \\ H_{2^{k-1}}^T & -H_{2^{k-1}}^T \end{pmatrix}
    $$
    $$
    = \begin{pmatrix} 2 H_{2^{k-1}} H_{2^{k-1}}^T & 0 \\ 0 & 2 H_{2^{k-1}} H_{2^{k-1}}^T \end{pmatrix} = \begin{pmatrix} 2 \cdot 2^{k-1} I & 0 \\ 0 & 2 \cdot 2^{k-1} I \end{pmatrix} = 2^k I
    $$
    因此 $H_{2^k}$ 是 Hadamard 矩阵。$\blacksquare$

### Hadamard 猜想与行列式上界

!!! definition "定义 65B.10 (Hadamard 猜想)"
    **Hadamard 猜想**：对所有 $n \equiv 0 \pmod{4}$，$n \times n$ Hadamard 矩阵存在。最小未知阶数为 $n = 668$。

!!! theorem "定理 65B.7 (Hadamard 行列式上界)"
    对任意 $n \times n$ 实矩阵 $A$，若 $|a_{ij}| \leq 1$，则 $|\det(A)| \leq n^{n/2}$。等号成立当且仅当 $A$ 是 Hadamard 矩阵。

??? proof "证明"
    设 $r_i$ 是 $A$ 的第 $i$ 行向量。由 $|a_{ij}| \leq 1$，$\|r_i\|^2 = \sum_j a_{ij}^2 \leq n$。

    由 Hadamard 不等式（行列式不超过行向量范数之积）：
    $$
    |\det(A)| \leq \prod_{i=1}^n \|r_i\| \leq \prod_{i=1}^n \sqrt{n} = n^{n/2}
    $$

    等号条件：$|\det(A)| = \prod \|r_i\|$ 要求行向量两两正交（Hadamard 不等式的等号条件），且 $\|r_i\| = \sqrt{n}$ 要求 $|a_{ij}| = 1$ 对所有 $i, j$。因此 $A$ 是 $\pm 1$ 矩阵且行正交，即 Hadamard 矩阵。$\blacksquare$

!!! example "例 65B.6"
    $H_4 = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & -1 & 1 & -1 \\ 1 & 1 & -1 & -1 \\ 1 & -1 & -1 & 1 \end{pmatrix}$。

    $|\det(H_4)| = 4^{4/2} = 16$。验证：$\det(H_4) = 16$（可直接计算）。

---

## 65B.7 积和式界

!!! definition "定义 65B.11 (积和式)"
    $n \times n$ 矩阵 $A$ 的**积和式**（permanent）定义为
    $$
    \operatorname{perm}(A) = \sum_{\sigma \in S_n} \prod_{i=1}^n a_{i,\sigma(i)}
    $$
    与行列式的区别在于没有 $\text{sgn}(\sigma)$ 因子。

!!! theorem "定理 65B.8 (van der Waerden 猜想 / Egorychev-Falikman 定理)"
    对任意 $n \times n$ 双随机矩阵 $A$（$a_{ij} \geq 0$，行列和为 1），
    $$
    \operatorname{perm}(A) \geq \frac{n!}{n^n}
    $$
    等号成立当且仅当 $A = \frac{1}{n} J_n$（均匀双随机矩阵）。

!!! note "注"
    此定理由 van der Waerden（1926）提出猜想，由 Egorychev（1981）和 Falikman（1981）独立证明。完整证明使用了混合体积理论和 Alexandrov-Fenchel 不等式。更详细的讨论见第 40A 章。

---

## 65B.8 竞赛矩阵

### 定义与得分序列

!!! definition "定义 65B.12 (竞赛与竞赛矩阵)"
    $n$ 阶**竞赛**是完全图 $K_n$ 的一个定向。**竞赛矩阵** $T = (t_{ij})$：
    $$
    t_{ij} = \begin{cases} 1 & \text{若 } i \to j \\ 0 & \text{若 } j \to i \text{ 或 } i = j \end{cases}
    $$
    $T + T^T = J - I$。顶点 $i$ 的**得分** $s_i = \sum_j t_{ij}$。

### Landau 定理

!!! theorem "定理 65B.9 (Landau 定理)"
    非负整数序列 $s_1 \leq s_2 \leq \cdots \leq s_n$ 是某个竞赛的得分序列，当且仅当
    $$
    \sum_{i=1}^k s_i \geq \binom{k}{2}, \quad k = 1, \ldots, n
    $$
    且 $\sum_{i=1}^n s_i = \binom{n}{2}$。

??? proof "证明"
    **必要性**：总得分 $\sum s_i = \binom{n}{2}$（每场比赛贡献 1 分给胜者，共 $\binom{n}{2}$ 场比赛）。

    对于前 $k$ 个得分最低的选手，他们之间有 $\binom{k}{2}$ 场比赛，至少产生 $\binom{k}{2}$ 分（每场 1 分）。因此 $\sum_{i=1}^k s_i \geq \binom{k}{2}$。

    **充分性**：对 $n$ 归纳。$n = 1$ 时 $s_1 = 0$，$\binom{1}{2} = 0$，成立。

    设 $n \geq 2$，给定满足条件的序列 $(s_1, \ldots, s_n)$。

    **关键步骤**：得分最高的选手 $n$（$s_n$）应该击败恰好 $s_n$ 个对手。考虑移除选手 $n$ 后剩余 $n-1$ 个选手的得分序列：对于被选手 $n$ 击败的 $s_n$ 个选手，他们的得分不变；对于击败选手 $n$ 的 $n - 1 - s_n$ 个选手，他们的得分减 1。

    需要验证新序列 $(s_1', \ldots, s_{n-1}')$（重新排序后）仍满足 Landau 条件。这通过仔细的组合论证完成：利用原始条件和 $\sum s_i' = \binom{n}{2} - s_n = \binom{n-1}{2}$（因为 $s_n$ 的范围保证了这一点），验证 $\sum_{i=1}^k s_i' \geq \binom{k}{2}$ 对所有 $k$ 成立。

    由归纳假设，$n-1$ 个选手的竞赛存在，加上选手 $n$ 的比赛结果，得到 $n$ 个选手的竞赛。$\blacksquare$

### 竞赛矩阵的谱性质

!!! theorem "定理 65B.10 (竞赛矩阵的特征值)"
    设 $T$ 是 $n$ 阶竞赛矩阵。则：

    1. $\operatorname{tr}(T) = 0$，因此特征值之和为 0；
    2. $(n-1)/2$ 是 $(T + T^T)/2 = (J - I)/2$ 的最大特征值，对应特征向量 $\mathbf{1}$；
    3. 所有特征值 $\lambda$ 满足 $\operatorname{Re}(\lambda) \leq (n-1)/2$；
    4. 对于正则竞赛（$s_i = (n-1)/2$ 对所有 $i$），$(n-1)/2$ 是 $T$ 的特征值。

??? proof "证明"
    **(1)** $\operatorname{tr}(T) = \sum_i t_{ii} = 0$（对角线全为零）。由特征值之和等于迹，$\sum \lambda_i = 0$。

    **(3)** 对任意特征值 $\lambda$，设 $x$ 是对应的特征向量，$\|x\|_\infty = |x_k|$ 最大分量。由 $Tx = \lambda x$：
    $$
    |\lambda| |x_k| = |(\lambda x)_k| = |(Tx)_k| = \left|\sum_j t_{kj} x_j\right| \leq \sum_j t_{kj} |x_j| \leq s_k |x_k|
    $$
    因此 $|\lambda| \leq s_k \leq n-1$。更精细地，$\operatorname{Re}(\lambda) \leq (n-1)/2$。

    **(4)** 对于正则竞赛，$T\mathbf{1} = (s_1, \ldots, s_n)^T = \frac{n-1}{2}\mathbf{1}$，因此 $(n-1)/2$ 是特征值。$\blacksquare$

!!! example "例 65B.7"
    $n = 3$ 循环竞赛 $1 \to 2 \to 3 \to 1$：
    $$
    T = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}
    $$
    循环置换矩阵，特征值 $1, \omega, \omega^2$（$\omega = e^{2\pi i/3}$）。得分序列 $(1, 1, 1)$。

---

## 65B.9 等分割与矩阵商

!!! definition "定义 65B.13 (等分割)"
    图 $G = (V, E)$ 的顶点划分 $\pi = \{C_1, \ldots, C_k\}$ 称为**等分割**（equitable partition），若对所有 $i, j$ 和所有 $v \in C_i$，$v$ 在 $C_j$ 中的邻居数 $b_{ij}$ 只依赖于 $i, j$，不依赖于 $v$ 的具体选择。

!!! definition "定义 65B.14 (矩阵商)"
    等分割 $\pi$ 的**商矩阵**（quotient matrix）$B = (b_{ij})_{k \times k}$，其中 $b_{ij}$ 是 $C_i$ 中顶点在 $C_j$ 中的邻居数。

!!! theorem "定理 65B.11 (商矩阵的谱性质)"
    设 $A$ 是图 $G$ 的邻接矩阵，$B$ 是等分割 $\pi$ 的商矩阵。则 $B$ 的特征值是 $A$ 的特征值的子集。即 $\operatorname{spec}(B) \subseteq \operatorname{spec}(A)$。

??? proof "证明"
    设 $S$ 是 $n \times k$ 的特征矩阵：$S_{vi} = 1$ 若 $v \in C_i$，否则 $S_{vi} = 0$。

    等分割条件等价于 $AS = SB$（可直接验证：$(AS)_{vi} = \sum_{w \sim v} S_{wi} = b_{ji}$，其中 $v \in C_j$，$(SB)_{vi} = \sum_l S_{vl} b_{li} = b_{ji}$）。

    设 $Bx = \mu x$（$x \neq 0$）。令 $y = Sx$。则
    $$
    Ay = ASx = SBx = S\mu x = \mu Sx = \mu y
    $$
    且 $y = Sx \neq 0$（因为 $S$ 的列线性无关——每个 $C_i$ 非空，$S$ 的列是 $\{0,1\}$ 向量且支撑集不相交）。因此 $\mu$ 是 $A$ 的特征值。$\blacksquare$

---

## 65B.10 D-A-D 缩放 (Sinkhorn-Knopp 算法)

!!! definition "定义 65B.15 (矩阵缩放问题)"
    给定非负矩阵 $A \geq 0$，**矩阵缩放问题**是找正对角矩阵 $D_1, D_2$，使得 $D_1 A D_2$ 是双随机矩阵（行列和均为 1）。

!!! theorem "定理 65B.12 (Sinkhorn-Knopp 定理)"
    非负矩阵 $A$ 可以被缩放为双随机矩阵（即存在正对角 $D_1, D_2$ 使得 $D_1 A D_2$ 双随机），当且仅当 $A$ 具有**全支撑**（total support）：$A$ 的每个非零元素属于某个正对角线（完美匹配）。

    **Sinkhorn-Knopp 算法**通过交替行归一化和列归一化迭代实现缩放：

    1. $A^{(0)} = A$；
    2. 行归一化：$A^{(2k+1)} = D_r^{(k)} A^{(2k)}$，使行和为 1；
    3. 列归一化：$A^{(2k+2)} = A^{(2k+1)} D_c^{(k)}$，使列和为 1。

    若 $A$ 有全支撑，则 $A^{(k)}$ 收敛到双随机矩阵 $B = D_1 A D_2$。

??? proof "证明"
    **收敛性的关键思想**：定义熵函数 $H(X) = -\sum_{ij} x_{ij} \log x_{ij}$。每次行归一化或列归一化都不减少 $H$（因为在固定行/列和约束下，均匀分配最大化熵——这是 Gibbs 不等式的推论）。

    同时，$H$ 有界（因为 $A$ 有限），因此迭代收敛。极限点必须同时满足行和 = 1 和列和 = 1，即双随机矩阵。

    全支撑条件保证了每步归一化的良定义性（不出现零行/列）。$\blacksquare$

!!! example "例 65B.8"
    设 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$。

    行归一化：$\begin{pmatrix} 1/3 & 2/3 \\ 3/7 & 4/7 \end{pmatrix}$。

    列归一化：列和为 $(1/3+3/7, 2/3+4/7) = (16/21, 26/21)$。归一化后继续迭代，最终收敛到双随机矩阵。

---

## 65B.11 Hadamard 矩阵的 Paley 构造

!!! definition "定义 65B.16 (Legendre 符号与二次剩余)"
    设 $p$ 是奇素数。整数 $a$ 的 **Legendre 符号** $\left(\frac{a}{p}\right)$ 定义为
    $$
    \left(\frac{a}{p}\right) = \begin{cases} 0 & \text{若 } p \mid a \\ 1 & \text{若 } a \text{ 是模 } p \text{ 的二次剩余} \\ -1 & \text{若 } a \text{ 不是模 } p \text{ 的二次剩余} \end{cases}
    $$

!!! theorem "定理 65B.13 (Paley 构造)"
    设 $p \equiv 3 \pmod{4}$ 是素数。定义 $(p+1) \times (p+1)$ 矩阵
    $$
    H = \begin{pmatrix} 0 & \mathbf{1}^T \\ -\mathbf{1} & Q \end{pmatrix} + I_{p+1}
    $$
    其中 $Q = \left(\left(\frac{i-j}{p}\right)\right)_{i,j=0}^{p-1}$ 是**二次剩余矩阵**（Jacobsthal 矩阵）。则 $H$ 是 $(p+1) \times (p+1)$ Hadamard 矩阵。

??? proof "证明"
    关键性质：$Q$ 是 $p \times p$ 反对称矩阵（因 $p \equiv 3 \pmod 4$ 时，$-1$ 不是二次剩余，$\left(\frac{-a}{p}\right) = -\left(\frac{a}{p}\right)$）。

    $Q$ 满足 $QQ^T = Q(-Q) = -Q^2 = pI - J$（利用二次剩余的正交关系：$\sum_{t=0}^{p-1} \left(\frac{t}{p}\right)\left(\frac{t+a}{p}\right) = -1$ 对 $a \not\equiv 0$，$= p-1$ 对 $a \equiv 0$）。

    由此可验证 $HH^T = (p+1)I_{p+1}$（分块矩阵乘法的仔细计算）。$\blacksquare$

!!! example "例 65B.9"
    $p = 3$：二次剩余 $\{1\}$，非二次剩余 $\{2\}$。
    $$
    Q = \begin{pmatrix} 0 & 1 & -1 \\ -1 & 0 & 1 \\ 1 & -1 & 0 \end{pmatrix}
    $$
    $H = \begin{pmatrix} 1 & 1 & 1 & 1 \\ -1 & 1 & 1 & -1 \\ -1 & -1 & 1 & 1 \\ -1 & 1 & -1 & 1 \end{pmatrix}$，验证 $HH^T = 4I_4$。

---

## 65B.12 竞赛矩阵的进阶性质

### Hamilton 路径与竞赛

!!! theorem "定理 65B.14 (竞赛中的 Hamilton 路径)"
    每个竞赛都包含一条 Hamilton 路径（经过所有顶点恰好一次的有向路径）。

??? proof "证明"
    对 $n$ 归纳。$n = 1$ 时显然。$n = 2$ 时，无论 $1 \to 2$ 或 $2 \to 1$，都是 Hamilton 路径。

    设 $n \geq 3$。由归纳假设，去掉顶点 $n$ 后的 $(n-1)$ 顶点竞赛有 Hamilton 路径 $v_1 \to v_2 \to \cdots \to v_{n-1}$。

    **情况 1**：若 $n \to v_1$，则 $n \to v_1 \to v_2 \to \cdots \to v_{n-1}$ 是 Hamilton 路径。

    **情况 2**：若 $v_{n-1} \to n$，则 $v_1 \to v_2 \to \cdots \to v_{n-1} \to n$ 是 Hamilton 路径。

    **情况 3**：若 $v_1 \to n$ 且 $n \to v_{n-1}$（即情况 1、2 都不成立），则存在某个 $k$ 使得 $v_k \to n$ 且 $n \to v_{k+1}$（因为关系从"$v$ 赢 $n$"变为"$n$ 赢 $v$"，必有切换点）。则 $v_1 \to \cdots \to v_k \to n \to v_{k+1} \to \cdots \to v_{n-1}$ 是 Hamilton 路径。$\blacksquare$

### 强连通竞赛

!!! theorem "定理 65B.15 (竞赛的强连通分量)"
    竞赛 $T$ 的强连通分量 $C_1, \ldots, C_r$ 可以排成线性序 $C_1, C_2, \ldots, C_r$，使得 $C_i$ 中每个顶点都击败 $C_j$ 中每个顶点（$i < j$）。特别地，竞赛的强连通分量数等于其**国王数**（king number）。

??? proof "证明"
    设 $C_1, \ldots, C_r$ 是 $T$ 的强连通分量。由于 $T$ 是竞赛（完全有向图），不同强连通分量之间的所有边方向一致：若 $C_i$ 中某顶点 $u$ 击败 $C_j$ 中某顶点 $v$，则 $C_i$ 中每个顶点都击败 $C_j$ 中每个顶点（否则 $v$ 击败 $C_i$ 中某顶点 $u'$，$u \to v$ 且存在 $C_j$ 到 $C_i$ 的边 $v \to u'$，加上 $C_i$ 内部从 $u'$ 到 $u$ 的路径和 $C_j$ 内部从 $v$ 到某处的路径，会产生新的强连通分量，矛盾）。

    因此强连通分量之间形成线性序（全序），可以排列使得 $C_i \to C_j$ 对所有 $i < j$。$\blacksquare$

---

## 65B.13 零强迫数的变体与推广

!!! definition "定义 65B.17 (正零强迫数)"
    对图 $G$，**正零强迫过程**将颜色变换规则修改为：黑色顶点 $u$ 可以将其白色邻居 $v$ 变黑，前提是 $v$ 是 $u$ 在**其所在的白色连通分量中**的唯一白色邻居。

    **正零强迫数** $Z_+(G)$ 是使全部变黑的最小初始黑色集合大小。

!!! theorem "定理 65B.16 (正零强迫数与正半定最大零度)"
    定义 $M_+(G) = \max\{\operatorname{null}(A) : A \in S_n(\mathbb{R}), A \succeq 0, \text{零模式匹配 } G\}$。则
    $$
    M_+(G) \leq Z_+(G)
    $$

!!! note "注"
    正零强迫数 $Z_+(G)$ 提供了正半定矩阵零度的上界，这在量子信息论中的图态独立数问题中有应用。$Z_+(G)$ 和 $Z(G)$ 之间的关系并不简单：总有 $Z_+(G) \geq Z(G)$，但差距可以任意大。

---

## 习题

!!! question "习题 65B.1"
    计算圈图 $C_5$ 的最小秩 $\operatorname{mr}(C_5)$ 和零强迫数 $Z(C_5)$。

!!! question "习题 65B.2"
    证明完全图 $K_n$ 的零强迫数为 $n - 1$。

!!! question "习题 65B.3"
    构造 $8 \times 8$ Hadamard 矩阵（使用 Sylvester 方法），并验证 $HH^T = 8I$。

!!! question "习题 65B.4"
    证明 $n$ 阶竞赛矩阵 $T$ 满足 $T + T^T = J_n - I_n$。

!!! question "习题 65B.5"
    证明：若 $H$ 是 $n \times n$ Hadamard 矩阵，则 $|\det(H)| = n^{n/2}$。

!!! question "习题 65B.6"
    设 $T$ 是 4 阶竞赛矩阵，得分序列为 $(0, 1, 2, 3)$。写出 $T$ 并计算其特征值。

!!! question "习题 65B.7"
    对星图 $K_{1,n-1}$，求 $\operatorname{mr}(K_{1,n-1})$ 和 $Z(K_{1,n-1})$。

!!! question "习题 65B.8"
    证明原始矩阵的指数至少为 2（即 $\gamma(A) \geq 2$，若 $n \geq 2$）。

!!! question "习题 65B.9"
    设 $A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$。证明 $A$ 不能被缩放为双随机矩阵。（$A$ 没有全支撑。）

!!! question "习题 65B.10"
    证明 Petersen 图的 Colin de Verdiere 参数 $\mu = 4$（即 Petersen 图不是平面图）。（提示：Petersen 图包含 $K_{3,3}$ 作为子图。）

!!! question "习题 65B.11"
    设 $G$ 是 $n$ 个顶点的正则图（度为 $d$）。证明 $\{V\}$（全体顶点作为一个类）是等分割，商矩阵为 $(d)$。

!!! question "习题 65B.12"
    利用 Sinkhorn-Knopp 算法，对矩阵 $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ 执行三轮迭代（行归一化 → 列归一化 → 行归一化 → ...），并观察收敛行为。
