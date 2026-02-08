# 第 62 章 矩阵补全问题

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · SVD(Ch11) · 优化(Ch25) · 图论(Ch27)

**本章脉络**：矩阵补全的一般框架 → 正定补全（弦图条件） → 低秩补全（Netflix 问题） → 核范数最小化 → 精确恢复条件 → 欧氏距离矩阵补全 → 惯性补全 → 算法

**延伸**：低秩矩阵补全是推荐系统（Netflix、Amazon）和协同过滤的数学基础；正定补全在空间统计（稀疏协方差估计）和凸优化（SDP 松弛的对偶）中有重要应用

</div>

矩阵补全问题的基本设定是：给定一个矩阵的**部分条目**，将缺失的条目"补全"，使得完成后的矩阵满足某种所需的性质——如正半定性、低秩、距离矩阵等。这一看似简单的问题在过去二十年中发展为一个深刻且应用广泛的研究领域。

在理论方面，矩阵补全将图论（稀疏模式）、凸优化（核范数最小化）、随机矩阵理论（随机采样）和代数几何（秩的代数结构）交织在一起。在应用方面，低秩矩阵补全是推荐系统（Netflix 奖金问题的核心数学模型）和协同过滤的理论基石。

---

## 62.1 矩阵补全的一般框架

<div class="context-flow" markdown>

**核心问题**：什么是矩阵补全问题？有哪些主要类型？

</div>

!!! definition "定义 62.1 (部分矩阵)"
    **部分矩阵**是一个 $m \times n$ 矩阵，其中只有部分条目被指定，其余条目未知。形式化地，给定指标集 $\Omega \subset \{1,\ldots,m\} \times \{1,\ldots,n\}$ 和值 $\{a_{ij} : (i,j) \in \Omega\}$，部分矩阵 $A_\Omega$ 是在 $\Omega$ 上取值 $a_{ij}$、在 $\Omega^c$ 上取值未定的矩阵。

!!! definition "定义 62.2 (矩阵补全问题)"
    **矩阵补全问题**：给定部分矩阵 $A_\Omega$，找到完整矩阵 $X \in \mathbb{R}^{m \times n}$，使得：

    (a) $X_{ij} = a_{ij}$，$\forall (i,j) \in \Omega$（与已知条目一致）。

    (b) $X$ 满足某种所需性质 $\mathcal{P}$。

    根据性质 $\mathcal{P}$ 的不同，得到不同类型的矩阵补全问题。

!!! definition "定义 62.3 (矩阵补全的主要类型)"
    **(a) 正半定补全**：$\mathcal{P}$ = "$X$ 是正半定的"。

    **(b) 低秩补全**：$\mathcal{P}$ = "$\mathrm{rank}(X) \leq r$"。

    **(c) 欧氏距离矩阵补全**：$\mathcal{P}$ = "$X$ 是欧氏距离矩阵"。

    **(d) 非负补全**：$\mathcal{P}$ = "$X \geq 0$"。

    **(e) 全正补全**：$\mathcal{P}$ = "$X$ 是全正矩阵（所有子矩阵行列式非负）"。

!!! definition "定义 62.4 (稀疏模式图)"
    部分矩阵 $A_\Omega$ 的**稀疏模式图**（sparsity pattern graph）$G = (V, E)$ 定义为：

    - **对称情形**（$m = n$，$A$ 对称）：$V = \{1, \ldots, n\}$，$(i,j) \in E$ 当且仅当 $(i,j) \in \Omega$ 且 $i \neq j$。对角元素总假定已知。
    - **一般情形**：$G$ 是二部图，左顶点 $\{r_1, \ldots, r_m\}$，右顶点 $\{c_1, \ldots, c_n\}$，$(r_i, c_j) \in E$ 当且仅当 $(i,j) \in \Omega$。

!!! example "例 62.1"
    考虑 $4 \times 4$ 对称部分矩阵，$\Omega$ 对应的图 $G$ 为：

    $$A = \begin{pmatrix}2 & 1 & ? & ?\\1 & 3 & 2 & ?\\? & 2 & 4 & 1\\? & ? & 1 & 5\end{pmatrix}.$$

    稀疏模式图 $G$：顶点 $\{1,2,3,4\}$，边 $\{(1,2), (2,3), (3,4)\}$——这是一条路径图 $P_4$。

    问题：能否填入 $a_{13}, a_{14}, a_{24}$（及其对称位置）使得 $A$ 正半定？

---

## 62.2 正定补全

<div class="context-flow" markdown>

**核心问题**：给定对称部分矩阵，何时能补全为正半定矩阵？稀疏模式图的什么性质决定了补全的可行性？

</div>

!!! definition "定义 62.5 (弦图)"
    无向图 $G$ 称为**弦图**（chordal graph），如果 $G$ 的每个长度 $\geq 4$ 的环都有一条**弦**（连接环上两个不相邻顶点的边）。

    等价地，$G$ 不包含长度 $\geq 4$ 的**无弦环**（induced cycle）作为诱导子图。

!!! example "例 62.2"
    以下图是弦图：完全图 $K_n$；树；区间图；路径图 $P_n$。

    以下图不是弦图：长度 $\geq 4$ 的环 $C_n$（$n \geq 4$）；Petersen 图。

    4 个顶点的环 $C_4$（正方形）：$1-2-3-4-1$，没有弦（缺少 $(1,3)$ 或 $(2,4)$ 边），不是弦图。

!!! theorem "定理 62.1 (Grone 定理, 1984)"
    设 $A_\Omega$ 是 $n \times n$ 对称部分矩阵，稀疏模式图为 $G$。则以下等价：

    **(a)** 对**每一个** $A_\Omega$（只要 $A_\Omega$ 的每个全指定主子矩阵都是正半定的），都存在正半定补全 $X \succeq 0$。

    **(b)** $G$ 是**弦图**。

??? proof "证明"
    **"$(b) \Rightarrow (a)$"：弦图保证补全存在。**

    弦图的关键性质是具有**完美消除序**（perfect elimination ordering, PEO）：存在顶点的排列 $v_1, v_2, \ldots, v_n$，使得对每个 $i$，$v_i$ 在 $G[\{v_i, v_{i+1}, \ldots, v_n\}]$（诱导子图）中的邻居集构成团（完全子图）。

    **构造算法。** 按完美消除序的逆序（$v_n, v_{n-1}, \ldots, v_1$）逐步填充缺失条目。

    在第 $k$ 步，考虑顶点 $v_k$。由 PEO 的性质，$v_k$ 在 $G[\{v_k, \ldots, v_n\}]$ 中的邻居构成团 $C_k$。此时 $C_k$ 中的所有条目已被确定（或为原始数据）。需要确定 $v_k$ 与 $\{v_{k+1}, \ldots, v_n\} \setminus C_k$ 之间的条目。

    将问题化归为：已知 Schur 补的条件，选择使得 Schur 补非负的填充值。具体地，将矩阵按 $\{v_k\}$ 和其余顶点分块：

    $$\begin{pmatrix}a_{kk} & b^\top\\b & S\end{pmatrix} \succeq 0 \iff a_{kk} \geq 0 \text{ 且 } S - bb^\top / a_{kk} \succeq 0.$$

    由 $C_k$ 是团，$b$ 中已知部分对应的 $S$ 的子矩阵是已确定的正半定矩阵。可以选择 $b$ 的未知部分使得条件满足。

    **"$(a) \Rightarrow (b)$"：非弦图存在不可补全的部分矩阵。**

    如果 $G$ 不是弦图，则存在长度 $k \geq 4$ 的无弦诱导环 $C_k$。可以构造一个 $k \times k$ 对称部分矩阵，其已知条目的每个主子矩阵都是正定的，但不存在正半定补全。

    经典构造：在 $C_4$（$1-2-3-4-1$）上，取对角元为 1，已知非对角元为 $a_{12} = a_{23} = a_{34} = a_{41} = \cos\theta$（$\theta$ 接近 0）。则每个 $2 \times 2$ 已知子矩阵 $\begin{pmatrix}1 & \cos\theta\\\cos\theta & 1\end{pmatrix} \succ 0$。但要使整个 $4 \times 4$ 矩阵正半定，需要 $a_{13}$ 和 $a_{24}$ 满足非常紧的约束，当 $\theta$ 足够小时无法同时满足。$\blacksquare$

!!! theorem "定理 62.2 (最大行列式补全)"
    在 Grone 定理的条件下（$G$ 是弦图），所有正半定补全中，存在**唯一的最大行列式补全** $X^*$：

    $$X^* = \arg\max\{\det(X) : X \succeq 0,\, X_{ij} = a_{ij}\, \forall (i,j) \in \Omega \cup \{(i,i)\}\}.$$

    最大行列式补全有以下性质：$(X^*)^{-1}_{ij} = 0$，$\forall (i,j) \notin \Omega \cup \{(i,i)\}$（即 $X^*$ 的逆矩阵在 $G$ 的补图上为零）。

??? proof "证明"
    对数行列式 $\log\det(X)$ 在正定锥上是严格凹函数。在仿射约束（$X_{ij} = a_{ij}$）和半定约束（$X \succeq 0$）下，最大化凹函数有唯一解。

    KKT 条件给出：$X^{-1}_{ij} = 0$（对 $(i,j) \notin \Omega$，$i \neq j$）。这是因为 $\frac{\partial}{\partial X_{ij}}\log\det(X) = (X^{-1})_{ij}$，而未约束的变量 $X_{ij}$（$(i,j) \notin \Omega$）的梯度在最优解处为零。$\blacksquare$

!!! example "例 62.3"
    **路径图上的正定补全。** 回到例 62.1：

    $$A = \begin{pmatrix}2 & 1 & ? & ?\\1 & 3 & 2 & ?\\? & 2 & 4 & 1\\? & ? & 1 & 5\end{pmatrix}.$$

    路径图 $P_4$ 是弦图（树是弦图）。按完美消除序 $4, 3, 2, 1$：

    先确定 $a_{24}$：使用 $\{2,3,4\}$ 子矩阵。$\begin{pmatrix}3 & 2 & x\\2 & 4 & 1\\x & 1 & 5\end{pmatrix} \succeq 0$。最大行列式补全要求 $3 \times 3$ 矩阵的逆在 $(2,4)$ 位置为零。

    $\begin{pmatrix}3 & 2 & x\\2 & 4 & 1\\x & 1 & 5\end{pmatrix}^{-1}$ 的 $(1,3)$ 元素 = $\frac{2 \cdot 1 - 4x}{D}$（$D$ 为行列式）= 0，得 $x = 1/2$。

    类似地确定 $a_{13}$ 和 $a_{14}$。

---

## 62.3 低秩矩阵补全问题

<div class="context-flow" markdown>

**核心问题**：给定秩-$r$ 矩阵的一小部分条目，能否恢复整个矩阵？需要多少条目？

</div>

!!! definition "定义 62.6 (低秩矩阵补全)"
    设 $M \in \mathbb{R}^{m \times n}$ 是未知的秩-$r$ 矩阵。我们观测到 $\Omega \subset \{1,\ldots,m\} \times \{1,\ldots,n\}$ 上的条目 $\{M_{ij} : (i,j) \in \Omega\}$。**低秩矩阵补全**问题是从这些部分观测中恢复 $M$：

    $$\min_X \mathrm{rank}(X) \quad \text{s.t.} \quad X_{ij} = M_{ij},\, \forall (i,j) \in \Omega.$$

!!! theorem "定理 62.3 (自由度计数)"
    $m \times n$ 的秩-$r$ 矩阵有

    $$r(m + n - r)$$

    个自由度。因此，精确恢复至少需要 $|\Omega| \geq r(m + n - r)$ 个观测条目。

??? proof "证明"
    秩-$r$ 矩阵可以参数化为 $M = UV^\top$，$U \in \mathbb{R}^{m \times r}$，$V \in \mathbb{R}^{n \times r}$。$U$ 有 $mr$ 个参数，$V$ 有 $nr$ 个参数。但分解不唯一：$M = (UQ)(VQ)^\top$（$Q \in GL(r)$），$Q$ 有 $r^2$ 个参数。因此自由度为 $mr + nr - r^2 = r(m + n - r)$。$\blacksquare$

!!! definition "定义 62.7 (不相干性条件)"
    秩-$r$ 矩阵 $M$ 的 SVD 为 $M = U\Sigma V^\top$，$U \in \mathbb{R}^{m \times r}$，$V \in \mathbb{R}^{n \times r}$。定义不相干性参数：

    $$\mu_0 = \max\!\Bigl(\frac{m}{r}\max_i \|U^\top e_i\|^2,\, \frac{n}{r}\max_j \|V^\top e_j\|^2\Bigr).$$

    $\mu_0 \in [1, m/r]$（或 $n/r$）。$\mu_0 = 1$ 表示奇异向量均匀分布（"不相干"），$\mu_0$ 大表示能量集中在少数行/列（"相干"）。

    第二个不相干参数：

    $$\mu_1 = \frac{\sqrt{mn}}{r}\max_{i,j}|(UV^\top)_{ij}|.$$

不相干性条件排除了"不可能恢复"的病态情形。

!!! example "例 62.4"
    **不相干性的直觉。** 考虑 $m = n$，$r = 1$。

    - **最不相干**：$M = \frac{1}{n}\mathbf{1}\mathbf{1}^\top$（所有条目相同）。$\mu_0 = 1$。观测任何一个条目就能恢复整个矩阵。
    - **最相干**：$M = e_1 e_1^\top$（只有 $(1,1)$ 位置非零）。$\mu_0 = n$。如果没有观测到 $(1,1)$，完全无法恢复。

    不相干性保证了矩阵的信息分散在所有条目中，而非集中在少数位置。

---

## 62.4 核范数最小化

<div class="context-flow" markdown>

**核心问题**：如何将非凸的秩最小化问题松弛为可高效求解的凸优化问题？核范数为何是秩的良好凸代理？

</div>

!!! definition "定义 62.8 (核范数)"
    矩阵 $X$ 的**核范数**（nuclear norm，或迹范数）定义为

    $$\|X\|_* = \sum_{i=1}^{\min(m,n)} \sigma_i(X) = \mathrm{tr}\!\bigl(\sqrt{X^\top X}\bigr),$$

    其中 $\sigma_i(X)$ 是 $X$ 的奇异值。

!!! theorem "定理 62.4 (核范数是秩的凸包)"
    在单位谱范数球 $\{X : \|X\| \leq 1\}$ 上，核范数 $\|X\|_*$ 是 $\mathrm{rank}(X)$ 的**凸包**（最紧凸下界）。

    这类似于 $L_1$ 范数是 $L_0$ "范数"（非零元素个数）在单位 $L_\infty$ 球上的凸包。

??? proof "证明"
    设 $X$ 满足 $\|X\| \leq 1$，即 $\sigma_1(X) \leq 1$。则 $\|X\|_* = \sum \sigma_i \leq \mathrm{rank}(X) \cdot 1 = \mathrm{rank}(X)$（因为每个非零奇异值 $\leq 1$）。

    另一方面，当 $X$ 的所有非零奇异值恰好等于 1 时，$\|X\|_* = \mathrm{rank}(X)$。这些点是秩函数在单位谱范数球上的"极端点"。

    因此，在约束 $\|X\| \leq 1$ 下，$\|X\|_*$ 是 $\mathrm{rank}(X)$ 的最紧凸下界。$\blacksquare$

!!! definition "定义 62.9 (核范数最小化补全)"
    **核范数最小化**矩阵补全：

    $$\min_X \|X\|_* \quad \text{s.t.} \quad X_{ij} = M_{ij},\, \forall (i,j) \in \Omega.$$

    这是一个凸优化问题（核范数是凸函数，约束是仿射的），可以表示为半定规划（SDP）。

!!! theorem "定理 62.5 (核范数最小化的 SDP 表示)"
    核范数最小化等价于以下半定规划：

    $$\min_{X, W_1, W_2} \frac{1}{2}\bigl(\mathrm{tr}(W_1) + \mathrm{tr}(W_2)\bigr) \quad \text{s.t.} \quad \begin{pmatrix}W_1 & X\\X^\top & W_2\end{pmatrix} \succeq 0, \quad X_{ij} = M_{ij}\, \forall (i,j) \in \Omega.$$

??? proof "证明"
    利用 Schur 补：$\begin{pmatrix}W_1 & X\\X^\top & W_2\end{pmatrix} \succeq 0$ 要求 $W_1 \succeq XW_2^{-1}X^\top$（当 $W_2 \succ 0$）。最小化 $\mathrm{tr}(W_1) + \mathrm{tr}(W_2)$ 等价于最小化 $X$ 的奇异值之和。

    具体地，设 $X = U\Sigma V^\top$，取 $W_1 = U\Sigma U^\top$，$W_2 = V\Sigma V^\top$，则 $\mathrm{tr}(W_1) + \mathrm{tr}(W_2) = 2\sum \sigma_i = 2\|X\|_*$，且半定约束满足。可以验证这是最优解。$\blacksquare$

!!! theorem "定理 62.6 (Candes-Recht 定理, 2009)"
    设 $M \in \mathbb{R}^{n \times n}$ 是秩-$r$ 矩阵，不相干性参数为 $\mu_0$。如果 $|\Omega|$ 个条目从均匀分布中独立随机采样，且

    $$|\Omega| \geq C\,\mu_0^2\, r\, n \,(\log n)^2,$$

    则核范数最小化以概率至少 $1 - n^{-3}$ **精确恢复** $M$（即核范数最小化的解等于 $M$）。

??? proof "证明（概要）"
    证明分为几个关键步骤。

    **第一步：对偶证书。** 核范数最小化恢复 $M$ 的充分条件是存在**对偶证书** $Y \in \mathbb{R}^{n \times n}$：

    (i) $\mathcal{P}_\Omega(Y) = Y$（$Y$ 仅在 $\Omega$ 上有非零元素）；

    (ii) $\mathcal{P}_T(Y) = UV^\top$（$Y$ 在 $M$ 的切空间 $T$ 上的投影等于 $M/\|M\|_*$ 的次梯度方向）；

    (iii) $\|\mathcal{P}_{T^\perp}(Y)\| < 1$（$Y$ 在法空间上的谱范数严格小于 1）。

    **第二步：Golfing scheme 构造。** 通过随机矩阵的迭代逼近（"高尔夫球"方案），构造满足条件 (i)-(iii) 的对偶证书。

    **第三步：矩阵浓度不等式。** 利用矩阵 Bernstein 不等式（第 57 章）控制随机构造的误差，保证对偶证书以高概率满足所需条件。

    采样数 $|\Omega| \geq C\mu_0^2 rn(\log n)^2$ 的来源：$rn$ 来自自由度计数，$\mu_0^2$ 来自不相干性，$(\log n)^2$ 来自概率的 union bound。$\blacksquare$

!!! example "例 62.5"
    **Netflix 问题的参数估计。** Netflix 数据集：$m \approx 500,000$ 用户，$n \approx 17,770$ 部电影，观测到 $|\Omega| \approx 10^8$ 个评分。

    如果评分矩阵近似为秩-$r$（$r \approx 10$-$50$），自由度约 $r(m+n) \approx 5 \times 10^6$-$2.5 \times 10^7$。$|\Omega|/r(m+n) \approx 4$-$20$，满足恢复的量级要求。

---

## 62.5 精确恢复条件

<div class="context-flow" markdown>

**核心问题**：精确恢复需要什么条件？不相干性条件的具体含义是什么？信息论下界是多少？

</div>

!!! theorem "定理 62.7 (信息论下界)"
    任何算法要从随机采样的条目中恢复一般的秩-$r$ 矩阵 $M \in \mathbb{R}^{m \times n}$（$m \leq n$），至少需要

    $$|\Omega| \geq r \cdot \max(m, n)$$

    个观测。这是因为秩-$r$ 矩阵有 $r(m+n-r) \geq r \cdot \max(m,n)/2$ 个自由度。

!!! theorem "定理 62.8 (Candes-Tao 改进, 2010)"
    在 Candes-Recht 定理的基础上，Candes 和 Tao 将 $(\log n)^2$ 改进为 $\log n$：

    $$|\Omega| \geq C\,\mu_0\, r\, n\, \log n$$

    个均匀随机采样即足以保证精确恢复（以高概率）。不相干性参数从 $\mu_0^2$ 降为 $\mu_0$。

!!! theorem "定理 62.9 (不相干性的必要性)"
    不相干性条件不仅是技术性假设，而且是**本质必要的**。存在相干矩阵，即使观测了 $n^{2-\epsilon}$ 个条目也无法恢复。

    **例子。** $M = e_1 v^\top$（秩-1，但 $\mu_0 = n$）。$M$ 的信息完全集中在第一行。如果没有采样到第一行的足够多条目，无法恢复 $v$。

!!! definition "定义 62.10 (采样算子)"
    **采样算子** $\mathcal{P}_\Omega: \mathbb{R}^{m \times n} \to \mathbb{R}^{m \times n}$ 定义为

    $$(\mathcal{P}_\Omega(X))_{ij} = \begin{cases}X_{ij}, & (i,j) \in \Omega,\\0, & (i,j) \notin \Omega.\end{cases}$$

    采样算子是正交投影：$\mathcal{P}_\Omega^2 = \mathcal{P}_\Omega$，$\mathcal{P}_\Omega = \mathcal{P}_\Omega^\top$（关于 Frobenius 内积）。

!!! theorem "定理 62.10 (RIP for matrices)"
    如果 $\Omega$ 满足**矩阵受限等距性质**（matrix RIP）：对所有秩不超过 $2r$ 的矩阵 $X$，

    $$(1-\delta)\|X\|_F^2 \leq \frac{mn}{|\Omega|}\|\mathcal{P}_\Omega(X)\|_F^2 \leq (1+\delta)\|X\|_F^2,$$

    其中 $\delta < 1/2$，则核范数最小化精确恢复所有秩-$r$ 矩阵。

    然而，与压缩感知不同，均匀随机采样不直接满足矩阵 RIP（需要不相干性条件）。

---

## 62.6 欧氏距离矩阵补全

<div class="context-flow" markdown>

**核心问题**：给定部分两两距离，能否恢复完整的距离矩阵？这与传感器网络定位有什么关系？

</div>

!!! definition "定义 62.11 (欧氏距离矩阵)"
    $n \times n$ 矩阵 $D$ 是**欧氏距离矩阵**（EDM），如果存在点 $p_1, \ldots, p_n \in \mathbb{R}^d$（某个 $d$）使得

    $$D_{ij} = \|p_i - p_j\|^2, \quad \forall\, i, j.$$

!!! theorem "定理 62.11 (EDM 与 Gram 矩阵的关系)"
    设 $D$ 是 EDM，$G$ 是对应点集的 **Gram 矩阵** $G_{ij} = p_i^\top p_j$（假设质心在原点，$\sum p_i = 0$）。则

    $$D_{ij} = G_{ii} + G_{jj} - 2G_{ij}, \quad \text{即 } D = \mathrm{diag}(G)\mathbf{1}^\top + \mathbf{1}\mathrm{diag}(G)^\top - 2G.$$

    反之，$G = -\frac{1}{2}JDJ$，其中 $J = I - \frac{1}{n}\mathbf{1}\mathbf{1}^\top$ 是中心化矩阵。

    $D$ 是 EDM 当且仅当 $G = -\frac{1}{2}JDJ \succeq 0$。

??? proof "证明"
    $D_{ij} = \|p_i - p_j\|^2 = \|p_i\|^2 + \|p_j\|^2 - 2p_i^\top p_j = G_{ii} + G_{jj} - 2G_{ij}$。

    反方向：$G = -\frac{1}{2}JDJ$ 可以直接验证。$J$ 是投影到 $\mathbf{1}^\perp$ 的投影矩阵，$JDJ$ 是"双中心化"操作。如果 $D$ 是 EDM，则 $G \succeq 0$（因为 $G$ 是中心化点集的 Gram 矩阵）。$\blacksquare$

!!! definition "定义 62.12 (EDM 补全问题)"
    **欧氏距离矩阵补全**：给定部分距离 $\{D_{ij} : (i,j) \in \Omega\}$，找到完整的 EDM $D$，使其对应的点集在 $\mathbb{R}^d$（$d$ 尽量小）中。

    通过 EDM-Gram 对应，这等价于：补全 Gram 矩阵 $G$ 使得 $G \succeq 0$ 且 $\mathrm{rank}(G) \leq d$。

!!! theorem "定理 62.12 (EDM 补全与正定补全的联系)"
    EDM 补全问题可以化归为正半定补全问题（对 Gram 矩阵 $G$）：

    $$\text{EDM 补全} \xrightarrow{G = -\frac{1}{2}JDJ} \text{PSD 补全（对 } G \text{）} + \text{秩约束 rank}(G) \leq d.$$

    放松秩约束后，EDM 补全归结为 SDP（半定规划）。

!!! example "例 62.6"
    **传感器网络定位。** $n$ 个传感器分布在二维平面上。每个传感器只能测量与附近（距离 $\leq R$）传感器的距离。已知 $m$ 个"锚点"（已知位置的传感器）的坐标。

    问题：从部分距离数据恢复所有传感器的位置。

    这是 EDM 补全问题的实例：$d = 2$（或 3），$\Omega$ 由距离 $\leq R$ 的传感器对确定。锚点提供了额外的约束（固定某些 $p_i$ 的值），消除了平移和旋转的不确定性。

    SDP 松弛方法将其表述为

    $$\min \mathrm{tr}(G) \quad \text{s.t.} \quad G \succeq 0,\, G_{ii} + G_{jj} - 2G_{ij} = D_{ij}\, \forall (i,j) \in \Omega.$$

---

## 62.7 算法

<div class="context-flow" markdown>

**核心问题**：除了 SDP（计算量大），有哪些高效的矩阵补全算法？

</div>

!!! definition "定义 62.13 (交替最小化)"
    **交替最小化**（Alternating Minimization）利用低秩分解 $X = UV^\top$（$U \in \mathbb{R}^{m \times r}$，$V \in \mathbb{R}^{n \times r}$），交替优化：

    $$U^{(t+1)} = \arg\min_U \sum_{(i,j) \in \Omega} (M_{ij} - (UV^{(t)\top})_{ij})^2,$$

    $$V^{(t+1)} = \arg\min_V \sum_{(i,j) \in \Omega} (M_{ij} - (U^{(t+1)}V^\top)_{ij})^2.$$

    每步是最小二乘问题，按行/列分解后可以高效求解。

!!! theorem "定理 62.13 (交替最小化的收敛性——Jain-Netrapalli-Sanghavi, 2013)"
    在不相干性条件和足够采样（$|\Omega| \geq C\mu_0^4 r^5 n \log n$）的假设下，适当初始化的交替最小化以几何速率收敛到 $M$：

    $$\|U^{(t)}V^{(t)\top} - M\|_F \leq \rho^t \cdot \|M\|_F,$$

    其中 $\rho < 1$ 是与采样率和不相干性相关的常数。

!!! definition "定义 62.14 (奇异值阈值化, SVT)"
    **奇异值阈值化**（Singular Value Thresholding）算法求解核范数正则化问题

    $$\min_X \frac{1}{2}\|\mathcal{P}_\Omega(X - M)\|_F^2 + \tau\|X\|_*.$$

    定义**软阈值化算子**：对矩阵 $Y = U\Sigma V^\top$（SVD），

    $$\mathcal{D}_\tau(Y) = U\,\mathrm{diag}(\max(\sigma_i - \tau, 0))\,V^\top.$$

    SVT 迭代：$X^{(k+1)} = \mathcal{D}_\tau\!\bigl(X^{(k)} + \delta\,\mathcal{P}_\Omega(M - X^{(k)})\bigr)$（$\delta > 0$ 为步长）。

??? proof "证明（SVT 的最优性条件）"
    目标函数 $f(X) = \frac{1}{2}\|\mathcal{P}_\Omega(X-M)\|_F^2 + \tau\|X\|_*$ 是凸函数。其次梯度条件为

    $$0 \in \mathcal{P}_\Omega(X - M) + \tau\,\partial\|X\|_*,$$

    其中 $\partial\|X\|_*$ 是核范数的次微分。对 $X = U\Sigma V^\top$（$\sigma_i > 0$），

    $$\partial\|X\|_* = \{UV^\top + W : U^\top W = 0,\, WV = 0,\, \|W\| \leq 1\}.$$

    SVT 算法是近端梯度法（proximal gradient method）的特例，$\mathcal{D}_\tau$ 是核范数的近端算子。收敛性由近端梯度法的一般理论保证。$\blacksquare$

!!! definition "定义 62.15 (Grassmann 流形上的梯度下降)"
    **Grassmann 流形方法**将低秩矩阵参数化为 $M = U\Sigma V^\top$，其中 $U$ 和 $V$ 分别在 Grassmann 流形 $\mathrm{Gr}(r, m)$ 和 $\mathrm{Gr}(r, n)$ 上。利用流形上的梯度下降（Riemannian gradient descent）优化

    $$\min_{U, \Sigma, V} \sum_{(i,j) \in \Omega} (M_{ij} - (U\Sigma V^\top)_{ij})^2.$$

    关键优势：流形方法自然地处理了分解的非唯一性（旋转不变性），避免了交替最小化中的不稳定性。

!!! theorem "定理 62.14 (梯度下降的全局收敛——Ma-Chen-etal, 2018)"
    在适当的不相干性和采样条件下，带谱初始化的梯度下降

    $$X^{(t+1)} = X^{(t)} - \eta \cdot \nabla_{\Omega} f(X^{(t)})$$

    （$\nabla_\Omega f$ 是采样条目上的梯度）以线性速率收敛到 $M$，只需 $O(|\Omega| r)$ 的计算量每步。

!!! example "例 62.7"
    **算法效率比较。** 对 $n = 10000$，$r = 10$，$|\Omega| = 5 \times 10^6$ 的问题：

    | 算法 | 每步复杂度 | 总迭代数 | 总时间（大约） |
    |------|-----------|----------|----------------|
    | SDP (内点法) | $O(n^3)$ | — | 不可行 |
    | SVT | $O(|\Omega| r)$ | $O(n/\delta)$ | 分钟级 |
    | 交替最小化 | $O(|\Omega| r)$ | $O(\log(1/\epsilon))$ | 秒级 |
    | 梯度下降 | $O(|\Omega| r)$ | $O(\log(1/\epsilon))$ | 秒级 |

---

## 62.8 应用

<div class="context-flow" markdown>

**核心问题**：矩阵补全理论如何在推荐系统、协同过滤和其他领域中应用？

</div>

### 62.8.1 推荐系统

!!! example "例 62.8"
    **Netflix 奖金问题。** Netflix 在 2006 年发起竞赛：给定用户对电影的部分评分（1-5 星），预测缺失的评分。

    数学模型：用户-电影评分矩阵 $M \in \mathbb{R}^{m \times n}$（$m$ 用户，$n$ 部电影），$M_{ij}$ 是用户 $i$ 对电影 $j$ 的评分。低秩假设：用户的偏好和电影的特征可以用少数"隐因子"描述，$\mathrm{rank}(M) \approx r$（$r$ 远小于 $m, n$）。

    分解 $M \approx UV^\top$：$U$ 的第 $i$ 行 $u_i \in \mathbb{R}^r$ 是用户 $i$ 的"偏好向量"，$V$ 的第 $j$ 行 $v_j \in \mathbb{R}^r$ 是电影 $j$ 的"特征向量"。预测评分 $\hat{M}_{ij} = u_i^\top v_j$。

    实际系统还加入偏置项：$\hat{M}_{ij} = \mu + b_i + c_j + u_i^\top v_j$。

### 62.8.2 协同过滤

!!! example "例 62.9"
    **Amazon 购物推荐。** 用户-商品交互矩阵（1 = 购买，0 = 未购买）近似低秩。通过矩阵补全预测用户可能感兴趣的商品。

    与 Netflix 问题的区别：这里是隐式反馈（购买/未购买），而非显式评分（1-5 星）。未购买不一定表示不喜欢——可能只是没看到。需要加权矩阵补全。

### 62.8.3 运动结构恢复

!!! example "例 62.10"
    **从运动恢复结构**（Structure from Motion, SfM）。$n$ 个 3D 点在 $m$ 个相机视角下的 2D 投影坐标形成测量矩阵 $W \in \mathbb{R}^{2m \times n}$。

    在仿射相机模型下，$W$ 的秩最多为 4（3D 点坐标 + 齐次坐标）。由于遮挡（某些点在某些视角不可见），$W$ 是部分观测的。

    矩阵补全 + SVD 因子分解 = 同时恢复 3D 点坐标和相机参数。这是计算机视觉中 SfM 管道的数学基础。

### 62.8.4 功率系统状态估计

!!! example "例 62.11"
    **电力系统中的矩阵补全。** 智能电网中的相量测量单元（PMU）提供电压和电流的相量数据。由于安装成本，只有部分节点有 PMU。节点电压矩阵在稳态下近似低秩。通过矩阵补全，可以从部分 PMU 数据恢复全网状态。

---

**本章要点总结：**

1. 矩阵补全问题从部分观测条目恢复满足特定性质的完整矩阵。
2. Grone 定理：正定补全存在当且仅当稀疏模式图是弦图。
3. 低秩矩阵补全需要不相干性条件：矩阵的信息不能集中在少数行/列。
4. 核范数是秩的最佳凸代理；核范数最小化在 $O(rn\log n)$ 个随机样本下精确恢复。
5. 欧氏距离矩阵补全通过 Gram 矩阵化归为正定补全加秩约束。
6. 高效算法包括交替最小化、SVT 和流形梯度下降，规模可达百万级。
7. 核心应用包括推荐系统（Netflix）、传感器定位和运动结构恢复。
