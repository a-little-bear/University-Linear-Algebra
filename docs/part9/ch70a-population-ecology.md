# 第 70A 章 种群生态学与流行病学中的线性代数

<div class="context-flow" markdown>

**前置**：特征值与特征向量(Ch6) · 非负矩阵与 Perron-Frobenius 定理(Ch17) · 矩阵指数(Ch13) · 微分方程(Ch26) · M-矩阵(Ch17)

**本章脉络**：Leslie 矩阵(年龄结构种群) $\to$ 渐近增长率、稳定年龄分布与繁殖价值 $\to$ 灵敏度与弹性分析 $\to$ Lefkovitch 阶段结构模型 $\to$ Lotka-Volterra 系统线性化 $\to$ 多物种群落矩阵与 May 稳定性定理 $\to$ 房室模型(药代动力学) $\to$ SIR/SEIR 流行病模型 $\to$ 次代矩阵与基本再生数 $R_0$ $\to$ 群体免疫阈值

**延伸**：Leslie 矩阵模型是种群生态学和渔业管理的标准工具；灵敏度与弹性分析为保护生物学提供管理策略优先级；May 的多样性-稳定性理论是理论生态学的里程碑；房室模型是药代动力学和放射性示踪剂动力学的基础；次代矩阵方法是现代传染病数学建模的核心框架

</div>

生物学在 20 世纪下半叶经历了深刻的数学化转变。从种群动力学到流行病学，线性代数提供了不可或缺的分析工具。Leslie 矩阵将年龄结构种群的增长编码为矩阵乘法，Perron-Frobenius 定理直接给出长期增长率和稳定年龄分布。在生态学中，Lotka-Volterra 系统的稳定性分析归结为 Jacobi 矩阵的特征值问题，而 May 的随机矩阵方法揭示了生态系统多样性与稳定性之间的深刻张力。在流行病学中，基本再生数 $R_0$ 可以表示为次代矩阵的谱半径，从而将非负矩阵理论与公共卫生决策紧密联系。

本章系统展示矩阵方法如何贯穿种群生态学和流行病学——从年龄结构种群到多物种群落，从药代动力学到传染病传播。

---

## 70A.1 Leslie 矩阵与种群动力学

<div class="context-flow" markdown>

**核心问题**：一个具有年龄结构的种群，各年龄组有不同的繁殖力和存活率，其长期增长行为如何？

</div>

P. H. Leslie 于 1945 年提出了以矩阵形式描述年龄结构种群动力学的方法。这个模型至今仍是种群生态学和渔业管理的标准工具。

!!! definition "定义 70A.1 (Leslie 矩阵)"
    将种群按年龄分为 $n$ 个年龄组。设 $\mathbf{n}_t = (n_1(t), n_2(t), \ldots, n_n(t))^T$ 为第 $t$ 个时间步各年龄组的个体数。定义 **Leslie 矩阵**：

    $$L = \begin{pmatrix} f_1 & f_2 & f_3 & \cdots & f_{n-1} & f_n \\ s_1 & 0 & 0 & \cdots & 0 & 0 \\ 0 & s_2 & 0 & \cdots & 0 & 0 \\ \vdots & & \ddots & & & \vdots \\ 0 & 0 & \cdots & s_{n-1} & 0 & 0 \end{pmatrix}$$

    其中：

    - $f_i \geq 0$ 为第 $i$ 个年龄组的**繁殖力**（每个个体在一个时间步内产生的后代数）
    - $s_i \in (0, 1]$ 为从第 $i$ 个年龄组存活到第 $i+1$ 个年龄组的**存活率**

    种群动力学方程为：

    $$\mathbf{n}_{t+1} = L \mathbf{n}_t$$

    因此 $\mathbf{n}_t = L^t \mathbf{n}_0$。

!!! theorem "定理 70A.1 (Leslie 矩阵的特征方程)"
    Leslie 矩阵 $L$ 的特征方程为：

    $$\lambda^n = f_1 \lambda^{n-1} + f_2 s_1 \lambda^{n-2} + f_3 s_1 s_2 \lambda^{n-3} + \cdots + f_n s_1 s_2 \cdots s_{n-1}$$

    等价地，令 $\ell_i = s_1 s_2 \cdots s_{i-1}$（$\ell_1 = 1$）为存活到第 $i$ 个年龄组的概率，则特征方程为：

    $$1 = \sum_{i=1}^{n} f_i \ell_i \lambda^{-i}$$

??? proof "证明"
    计算 $\det(L - \lambda I)$。沿第一列展开：

    $$\det(L - \lambda I) = (f_1 - \lambda) \det(M_{11}) - s_1 \det(M_{21})$$

    其中 $M_{11}$ 和 $M_{21}$ 是相应的余子式矩阵。

    我们用归纳法证明特征多项式。设 $p_n(\lambda)$ 为 $n \times n$ Leslie 矩阵的特征多项式。

    **基础情形**：$n = 1$ 时，$L = (f_1)$，$p_1(\lambda) = f_1 - \lambda$，即 $\lambda = f_1$。

    **归纳步骤**：假设对 $n-1$ 阶 Leslie 矩阵结论成立。对 $n$ 阶情形，沿第一列展开 $\det(L - \lambda I)$：

    第一列只有两个非零元素：$(L - \lambda I)_{11} = f_1 - \lambda$ 和 $(L - \lambda I)_{21} = s_1$。

    余子式 $M_{11}$ 是一个 $(n-1) \times (n-1)$ 矩阵，其结构为带有 $f_2, \ldots, f_n$ 在第一行、$s_2, \ldots, s_{n-1}$ 在次对角线上的 Leslie 型矩阵减去 $\lambda I$。

    余子式 $M_{21}$ 的第一行为 $(f_2, f_3, \ldots, f_n)$，其余行构成下三角结构。

    通过仔细展开（沿第一列反复展开），得到：

    $$p_n(\lambda) = (-1)^n\left[\lambda^n - f_1\lambda^{n-1} - f_2 s_1 \lambda^{n-2} - \cdots - f_n s_1 s_2 \cdots s_{n-1}\right]$$

    令 $p_n(\lambda) = 0$，两端除以 $(-1)^n \lambda^n$（$\lambda \neq 0$），得：

    $$1 = \sum_{i=1}^n f_i \ell_i \lambda^{-i}$$

    其中 $\ell_i = \prod_{j=1}^{i-1} s_j$。

    $\blacksquare$

!!! theorem "定理 70A.2 (Leslie 矩阵的主特征值)"
    若 Leslie 矩阵 $L$ 至少有两个非零繁殖力 $f_i, f_j > 0$（$i < j$），且 $\gcd\{i : f_i > 0\} = 1$（即 $L$ 为本原矩阵），则 $L$ 有唯一的正主特征值 $\lambda_1 > 0$，且 $\lambda_1 > |\lambda_k|$ 对所有其他特征值 $\lambda_k$ 成立。

??? proof "证明"
    Leslie 矩阵 $L$ 是非负矩阵。我们分两步证明。

    **第一步：$L$ 不可约。** 不可约等价于其有向图强连通。Leslie 矩阵的有向图中，存在边 $i \to i+1$（对应 $s_i > 0$，$i = 1, \ldots, n-1$），以及从 $i$ 到 $1$ 的边（对应 $f_i > 0$）。只要所有 $s_i > 0$ 且至少存在一个 $f_i > 0$（特别是 $f_n > 0$ 或更一般地存在从某个年龄组回到第 1 组的路径），则图是强连通的。

    **第二步：$L$ 非周期。** 不可约非负矩阵 $L$ 的周期 $d = \gcd\{k : (L^k)_{11} > 0\}$。从状态 1 回到状态 1 的最短路径长度对应于 $f_i > 0$ 的那些 $i$ 值。因此 $d = \gcd\{i : f_i > 0\}$。条件 $\gcd\{i : f_i > 0\} = 1$ 保证 $d = 1$，即 $L$ 非周期。

    不可约且非周期的非负矩阵即为**本原矩阵**。由 Perron-Frobenius 定理，本原矩阵有唯一的正主特征值 $\lambda_1 = \rho(L)$，且严格大于所有其他特征值的模，对应的右特征向量 $\mathbf{v}_1 > 0$（分量全正）。

    $\blacksquare$

!!! example "例 70A.1"
    一个三龄组种群的 Leslie 矩阵：

    $$L = \begin{pmatrix} 0 & 4 & 3 \\ 0.5 & 0 & 0 \\ 0 & 0.25 & 0 \end{pmatrix}$$

    即第 1 龄组不繁殖（$f_1 = 0$），第 2 龄组繁殖力为 4，第 3 龄组为 3；从第 1 到第 2 龄组的存活率为 0.5，从第 2 到第 3 龄组的存活率为 0.25。

    特征方程：$\lambda^3 = 4 \times 0.5 \lambda + 3 \times 0.5 \times 0.25 = 2\lambda + 0.375$

    即 $\lambda^3 - 2\lambda - 0.375 = 0$。数值求解得主特征值 $\lambda_1 \approx 1.489$。

    这意味着种群长期按每个时间步增长 $48.9\%$ 的速率增长。$\gcd(2, 3) = 1$，故 $L$ 是本原矩阵。

---

## 70A.2 渐近增长率、稳定年龄分布与繁殖价值

<div class="context-flow" markdown>

**核心问题**：Leslie 矩阵的特征值和特征向量有什么生态学含义？

</div>

Perron-Frobenius 定理在种群生态学中有三个核心应用：确定长期增长率、稳定年龄分布和繁殖价值。

!!! theorem "定理 70A.3 (长期增长率与稳定年龄分布)"
    设 Leslie 矩阵 $L$ 本原，主特征值为 $\lambda_1$，对应的右特征向量为 $\mathbf{v}_1$（$L\mathbf{v}_1 = \lambda_1 \mathbf{v}_1$），左特征向量为 $\mathbf{w}_1$（$\mathbf{w}_1^T L = \lambda_1 \mathbf{w}_1^T$）。则：

    1. **渐近增长率**：$\lambda_1$ 是种群的渐近增长率，即 $\mathbf{n}_t \sim c \lambda_1^t \mathbf{v}_1$（$t \to \infty$）
    2. **稳定年龄分布**：$\mathbf{v}_1 / \|\mathbf{v}_1\|_1$ 是种群的稳定年龄分布——无论初始年龄结构如何，最终各年龄组的比例趋向 $\mathbf{v}_1$ 的归一化形式
    3. **繁殖价值**：$\mathbf{w}_1$ 是各年龄组的**繁殖价值**——$w_i$ 度量第 $i$ 个年龄组的一个个体对未来种群增长的贡献

??? proof "证明"
    设 $L$ 可对角化（本原矩阵在简单特征值时可对角化；即使有重复特征值，以下谱投影论证仍然成立）。$L$ 的谱分解为：

    $$L = \sum_{k=1}^{n} \lambda_k \frac{\mathbf{v}_k \mathbf{w}_k^T}{\mathbf{w}_k^T \mathbf{v}_k}$$

    其中 $\{\mathbf{v}_k\}$ 和 $\{\mathbf{w}_k\}$ 分别是右特征向量和左特征向量，满足 $\mathbf{w}_j^T \mathbf{v}_k = 0$（$j \neq k$）。

    因此：

    $$L^t = \sum_{k=1}^{n} \lambda_k^t \frac{\mathbf{v}_k \mathbf{w}_k^T}{\mathbf{w}_k^T \mathbf{v}_k}$$

    由于 $L$ 本原，$|\lambda_k / \lambda_1| < 1$ 对 $k \geq 2$，故：

    $$\frac{L^t}{\lambda_1^t} \to \frac{\mathbf{v}_1 \mathbf{w}_1^T}{\mathbf{w}_1^T \mathbf{v}_1} \quad (t \to \infty)$$

    因此：

    $$\mathbf{n}_t = L^t \mathbf{n}_0 \sim \lambda_1^t \frac{\mathbf{w}_1^T \mathbf{n}_0}{\mathbf{w}_1^T \mathbf{v}_1} \mathbf{v}_1$$

    这表明：

    **(1)** 总体增长率为 $\lambda_1$：$\|\mathbf{n}_t\|_1 \sim c \lambda_1^t$。

    **(2)** 年龄结构趋向 $\mathbf{v}_1$：$\mathbf{n}_t / \|\mathbf{n}_t\|_1 \to \mathbf{v}_1 / \|\mathbf{v}_1\|_1$。

    **(3)** 初始条件通过 $\mathbf{w}_1^T \mathbf{n}_0$ 影响总量，其中 $\mathbf{w}_1$ 的分量 $w_i$ 给出了第 $i$ 个年龄组一个个体对未来种群总量的权重。这正是 Fisher 定义的**繁殖价值**。

    收敛速度由次主特征值比 $|\lambda_2 / \lambda_1|$ 决定：该比值越小，收敛越快。

    $\blacksquare$

!!! definition "定义 70A.2 (净繁殖率 $R_0$)"
    Leslie 模型中的**净繁殖率**（Net Reproductive Rate）定义为：

    $$R_0 = \sum_{i=1}^{n} f_i \ell_i$$

    其中 $\ell_i = s_1 s_2 \cdots s_{i-1}$（$\ell_1 = 1$）为存活到第 $i$ 个年龄组的概率。

    $R_0$ 表示一个个体在其一生中期望产生的后代数。

    - $R_0 > 1$：种群增长（$\lambda_1 > 1$）
    - $R_0 = 1$：种群稳定（$\lambda_1 = 1$）
    - $R_0 < 1$：种群衰退（$\lambda_1 < 1$）

!!! theorem "定理 70A.4 ($R_0$ 与主特征值的关系)"
    令 $g(\lambda) = \sum_{i=1}^{n} f_i \ell_i \lambda^{-i}$。则 $g(1) = R_0$，$g(\lambda_1) = 1$（特征方程），且 $g$ 在 $(0, \infty)$ 上严格递减。因此：

    - $R_0 > 1 \iff \lambda_1 > 1$
    - $R_0 = 1 \iff \lambda_1 = 1$
    - $R_0 < 1 \iff \lambda_1 < 1$

??? proof "证明"
    函数 $g(\lambda) = \sum_{i=1}^{n} f_i \ell_i \lambda^{-i}$ 在 $\lambda > 0$ 上定义。

    $g'(\lambda) = \sum_{i=1}^{n} (-i) f_i \ell_i \lambda^{-i-1} < 0$

    （只要至少有一个 $f_i > 0$，即种群有繁殖），因此 $g$ 在 $(0, \infty)$ 上严格递减。

    由特征方程 $g(\lambda_1) = 1$，又 $g(1) = R_0$。

    若 $R_0 > 1$，即 $g(1) > 1 = g(\lambda_1)$，由 $g$ 严格递减知 $1 < \lambda_1$。

    同理，$R_0 = 1 \iff \lambda_1 = 1$，$R_0 < 1 \iff \lambda_1 < 1$。

    $\blacksquare$

!!! example "例 70A.2"
    对例 70A.1 中的 Leslie 矩阵，$\ell_1 = 1$，$\ell_2 = 0.5$，$\ell_3 = 0.5 \times 0.25 = 0.125$。

    $R_0 = 0 \times 1 + 4 \times 0.5 + 3 \times 0.125 = 0 + 2 + 0.375 = 2.375$

    $R_0 > 1$，与主特征值 $\lambda_1 \approx 1.489 > 1$ 一致。

    右特征向量（稳定年龄分布）：$\mathbf{v}_1 \propto (1, 0.336, 0.0564)^T$

    归一化后约为 $(0.718, 0.241, 0.041)^T$，即稳定状态下约 71.8% 的个体在第 1 龄组，24.1% 在第 2 龄组，4.1% 在第 3 龄组。

    左特征向量（繁殖价值）：$\mathbf{w}_1 \propto (1, 2.978, 4.020)^T$

    第 3 龄组个体的繁殖价值最高（约为第 1 龄组的 4 倍），因为它们即将进入高繁殖力的年龄段（第 3 龄组本身繁殖力为 3）。

---

## 70A.3 灵敏度与弹性分析

<div class="context-flow" markdown>

**核心问题**：种群增长率对生命周期参数（繁殖力、存活率）的变化有多敏感？哪些参数对管理最关键？

</div>

灵敏度分析和弹性分析是保护生物学和渔业管理中最重要的工具之一。它们通过特征值摄动理论来实现。

!!! definition "定义 70A.3 (灵敏度与弹性)"
    主特征值 $\lambda_1$ 对 Leslie 矩阵元素 $a_{ij}$ 的**灵敏度**为：

    $$\frac{\partial \lambda_1}{\partial a_{ij}} = \frac{w_i v_j}{\mathbf{w}^T \mathbf{v}}$$

    **弹性**（相对灵敏度）为：

    $$e_{ij} = \frac{a_{ij}}{\lambda_1} \cdot \frac{\partial \lambda_1}{\partial a_{ij}} = \frac{a_{ij} w_i v_j}{\lambda_1 \mathbf{w}^T \mathbf{v}}$$

    弹性矩阵的元素之和等于 1：$\sum_{ij} e_{ij} = 1$。

!!! theorem "定理 70A.5 (灵敏度公式的推导——特征值摄动)"
    设 $L$ 为本原非负矩阵，$\lambda_1$ 为其简单主特征值，$\mathbf{v}_1$ 和 $\mathbf{w}_1$ 分别为对应的右和左特征向量。则对 $L$ 的任意小摄动 $\delta L$：

    $$\delta \lambda_1 = \frac{\mathbf{w}_1^T (\delta L) \mathbf{v}_1}{\mathbf{w}_1^T \mathbf{v}_1} + O(\|\delta L\|^2)$$

    特别地，$\frac{\partial \lambda_1}{\partial a_{ij}} = \frac{w_i v_j}{\mathbf{w}^T \mathbf{v}}$。

??? proof "证明"
    设 $L(\varepsilon) = L + \varepsilon \delta L$，$\lambda_1(\varepsilon) = \lambda_1 + \varepsilon \delta\lambda_1 + O(\varepsilon^2)$，$\mathbf{v}_1(\varepsilon) = \mathbf{v}_1 + \varepsilon \delta\mathbf{v}_1 + O(\varepsilon^2)$。

    由 $L(\varepsilon)\mathbf{v}_1(\varepsilon) = \lambda_1(\varepsilon)\mathbf{v}_1(\varepsilon)$，展开一阶项：

    $$L \cdot \delta\mathbf{v}_1 + (\delta L) \cdot \mathbf{v}_1 = \lambda_1 \delta\mathbf{v}_1 + (\delta\lambda_1) \mathbf{v}_1$$

    两边左乘 $\mathbf{w}_1^T$：

    $$\mathbf{w}_1^T L \cdot \delta\mathbf{v}_1 + \mathbf{w}_1^T (\delta L) \cdot \mathbf{v}_1 = \lambda_1 \mathbf{w}_1^T \delta\mathbf{v}_1 + (\delta\lambda_1) \mathbf{w}_1^T \mathbf{v}_1$$

    由于 $\mathbf{w}_1^T L = \lambda_1 \mathbf{w}_1^T$，左边第一项 $= \lambda_1 \mathbf{w}_1^T \delta\mathbf{v}_1$，与右边第一项抵消：

    $$\mathbf{w}_1^T (\delta L) \mathbf{v}_1 = (\delta\lambda_1) \mathbf{w}_1^T \mathbf{v}_1$$

    故：

    $$\delta\lambda_1 = \frac{\mathbf{w}_1^T (\delta L) \mathbf{v}_1}{\mathbf{w}_1^T \mathbf{v}_1}$$

    取 $\delta L = E_{ij}$（第 $(i,j)$ 位为 1、其余为 0 的矩阵），得 $\mathbf{w}_1^T E_{ij} \mathbf{v}_1 = w_i v_j$，故：

    $$\frac{\partial \lambda_1}{\partial a_{ij}} = \frac{w_i v_j}{\mathbf{w}_1^T \mathbf{v}_1}$$

    $\blacksquare$

!!! theorem "定理 70A.6 (弹性之和等于 1)"
    弹性矩阵 $E = (e_{ij})$ 满足 $\sum_{i,j} e_{ij} = 1$。

??? proof "证明"
    由定义：

    $$\sum_{i,j} e_{ij} = \sum_{i,j} \frac{a_{ij} w_i v_j}{\lambda_1 \mathbf{w}^T \mathbf{v}} = \frac{1}{\lambda_1 \mathbf{w}^T \mathbf{v}} \sum_{i,j} w_i a_{ij} v_j = \frac{\mathbf{w}^T L \mathbf{v}}{\lambda_1 \mathbf{w}^T \mathbf{v}}$$

    由 $L\mathbf{v} = \lambda_1 \mathbf{v}$，得 $\mathbf{w}^T L \mathbf{v} = \lambda_1 \mathbf{w}^T \mathbf{v}$。

    故 $\sum_{i,j} e_{ij} = 1$。

    $\blacksquare$

弹性之和为 1 的性质使得弹性可以解释为各生命周期转移对种群增长率的**相对贡献比例**。

!!! example "例 70A.3"
    对例 70A.1 中的 Leslie 矩阵，$\lambda_1 \approx 1.489$，$\mathbf{v}_1 \propto (1, 0.336, 0.0564)^T$，$\mathbf{w}_1 \propto (1, 2.978, 4.020)^T$。

    灵敏度矩阵 $S_{ij} = \frac{w_i v_j}{\mathbf{w}^T \mathbf{v}}$，其中 $\mathbf{w}^T \mathbf{v} = 1 + 2.978 \times 0.336 + 4.020 \times 0.0564 \approx 2.228$。

    关键灵敏度：

    - $\partial \lambda_1 / \partial s_1 = w_2 v_1 / (\mathbf{w}^T \mathbf{v}) \approx 2.978 / 2.228 \approx 1.337$
    - $\partial \lambda_1 / \partial f_2 = w_1 v_2 / (\mathbf{w}^T \mathbf{v}) \approx 0.336 / 2.228 \approx 0.151$

    弹性：

    - $e_{21}$（$s_1$ 的弹性）$= s_1 \times 1.337 / 1.489 = 0.5 \times 1.337 / 1.489 \approx 0.449$
    - $e_{12}$（$f_2$ 的弹性）$= 4 \times 0.151 / 1.489 \approx 0.406$

    在这个例子中，第 1 龄组的存活率 $s_1$ 和第 2 龄组的繁殖力 $f_2$ 对种群增长率影响最大。

---

## 70A.4 Lefkovitch 阶段结构模型

<div class="context-flow" markdown>

**核心问题**：当生物体不能按年龄分组（如植物的幼苗、幼树、成树），而是按发育阶段分组时，如何推广 Leslie 模型？

</div>

许多生物体（特别是植物、无脊椎动物、爬行动物等）的生活史更适合用**阶段**（如幼虫、幼体、成体）而非年龄来描述。Lefkovitch（1965）提出了 Leslie 矩阵的阶段结构推广。

!!! definition "定义 70A.4 (Lefkovitch 阶段结构矩阵)"
    将种群按发育阶段分为 $n$ 个阶段。**Lefkovitch 矩阵**（或投影矩阵）$A$ 的一般形式为：

    $$A = \begin{pmatrix} p_1 + f_1 & f_2 & f_3 & \cdots & f_n \\ g_1 & p_2 & 0 & \cdots & 0 \\ 0 & g_2 & p_3 & \cdots & 0 \\ \vdots & & \ddots & \ddots & \vdots \\ 0 & 0 & \cdots & g_{n-1} & p_n \end{pmatrix}$$

    其中：

    - $f_i \geq 0$：第 $i$ 阶段的繁殖贡献
    - $p_i \in [0, 1)$：第 $i$ 阶段个体在下一时间步**留在同一阶段**的概率（滞留率）
    - $g_i \in (0, 1]$：第 $i$ 阶段个体**进入下一阶段**的概率（转移率）
    - 约束 $p_i + g_i \leq 1$（$g_n = 0$ 对最后阶段）

    与 Leslie 矩阵的关键区别：**对角线非零**——个体可以在同一阶段停留多个时间步。

!!! theorem "定理 70A.7 (Lefkovitch 模型的 Perron-Frobenius 性质)"
    若 Lefkovitch 矩阵 $A$ 满足：(1) 所有 $g_i > 0$；(2) 至少一个 $f_j > 0$；(3) $\gcd$ 条件使 $A$ 本原，则 Perron-Frobenius 定理适用：$A$ 有唯一正主特征值 $\lambda_1$，灵敏度和弹性分析公式完全相同。

??? proof "证明"
    Lefkovitch 矩阵 $A$ 是非负矩阵。条件 (1) 保证所有相邻阶段之间有正转移，条件 (2) 保证至少一个阶段有繁殖。

    **不可约性**：有向图中，$g_i > 0$ 给出边 $i \to i+1$；$f_j > 0$ 给出从阶段 $j$ 到阶段 1 的边。这保证了强连通性。

    **非周期性**：由于 $p_i > 0$ 对至少一个 $i$（或由其他路径长度组合），图中存在从某节点到自身的长度为 1 的路径（自环），因此周期为 1。事实上，只要有一个 $p_i > 0$，矩阵就是非周期的。

    因此 $A$ 本原，Perron-Frobenius 定理适用。灵敏度公式 $\partial \lambda_1 / \partial a_{ij} = w_i v_j / (\mathbf{w}^T \mathbf{v})$ 是对任意可微矩阵参数化的特征值摄动结果，不依赖于矩阵的特殊结构。

    $\blacksquare$

!!! example "例 70A.4"
    一种多年生植物的三阶段模型：种子库、幼苗、成树。

    $$A = \begin{pmatrix} 0.1 & 0 & 50 \\ 0.01 & 0.3 & 0 \\ 0 & 0.05 & 0.95 \end{pmatrix}$$

    - 种子库：$p_1 = 0.1$（10% 种子存活到下年仍为种子），$g_1 = 0.01$（1% 发芽为幼苗），成树繁殖贡献 $f_3 = 50$
    - 幼苗：$p_2 = 0.3$（30% 留在幼苗阶段），$g_2 = 0.05$（5% 长为成树）
    - 成树：$p_3 = 0.95$（95% 存活，成树寿命长）

    对角线元素 $p_1 = 0.1, p_2 = 0.3, p_3 = 0.95$ 反映了各阶段的滞留时间差异——成树平均在该阶段停留 $1/(1-0.95) = 20$ 个时间步。

    弹性分析通常发现：对长寿物种，成体存活率（$p_n$）的弹性最大；对短命多产物种，繁殖力（$f_i$）的弹性更大。

---

## 70A.5 Lotka-Volterra 与线性化

<div class="context-flow" markdown>

**核心问题**：非线性种群动力学系统（如捕食者-猎物模型）的稳定性如何通过线性化分析来判定？

</div>

!!! definition "定义 70A.5 (Lotka-Volterra 捕食者-猎物模型)"
    经典 Lotka-Volterra 系统描述猎物 ($x$) 与捕食者 ($y$) 的动态：

    $$\frac{dx}{dt} = \alpha x - \beta xy = x(\alpha - \beta y)$$

    $$\frac{dy}{dt} = \delta xy - \gamma y = y(\delta x - \gamma)$$

    其中 $\alpha, \beta, \gamma, \delta > 0$ 为正参数：$\alpha$ 为猎物内禀增长率，$\beta$ 为捕食率，$\gamma$ 为捕食者死亡率，$\delta$ 为捕食者的能量转化效率。

!!! theorem "定理 70A.8 (Lotka-Volterra 的均衡与稳定性)"
    Lotka-Volterra 系统有两个均衡点：

    1. **平凡均衡** $(0, 0)$：Jacobi 矩阵 $J_0 = \begin{pmatrix} \alpha & 0 \\ 0 & -\gamma \end{pmatrix}$，特征值 $\alpha > 0, -\gamma < 0$，是**鞍点**（不稳定）。

    2. **非平凡均衡** $(x^*, y^*) = \left(\frac{\gamma}{\delta}, \frac{\alpha}{\beta}\right)$：

    $$J^* = \begin{pmatrix} 0 & -\frac{\beta\gamma}{\delta} \\ \frac{\delta\alpha}{\beta} & 0 \end{pmatrix}$$

    特征值 $\lambda = \pm i\sqrt{\alpha\gamma}$，纯虚数，线性化给出**中心**。经典 Lotka-Volterra 系统有守恒量，实际上轨道是封闭曲线。

??? proof "证明"
    设 $f(x,y) = x(\alpha - \beta y)$，$g(x,y) = y(\delta x - \gamma)$。

    Jacobi 矩阵为：

    $$J = \begin{pmatrix} \partial f/\partial x & \partial f/\partial y \\ \partial g/\partial x & \partial g/\partial y \end{pmatrix} = \begin{pmatrix} \alpha - \beta y & -\beta x \\ \delta y & \delta x - \gamma \end{pmatrix}$$

    **在 $(0,0)$**：$J_0 = \begin{pmatrix} \alpha & 0 \\ 0 & -\gamma \end{pmatrix}$，特征值 $\alpha > 0$ 和 $-\gamma < 0$，一正一负，鞍点。

    **在 $(x^*, y^*) = (\gamma/\delta, \alpha/\beta)$**：

    $$J^* = \begin{pmatrix} \alpha - \beta(\alpha/\beta) & -\beta(\gamma/\delta) \\ \delta(\alpha/\beta) & \delta(\gamma/\delta) - \gamma \end{pmatrix} = \begin{pmatrix} 0 & -\beta\gamma/\delta \\ \delta\alpha/\beta & 0 \end{pmatrix}$$

    特征方程：$\lambda^2 - \text{tr}(J^*)\lambda + \det(J^*) = 0$，即 $\lambda^2 + \alpha\gamma = 0$。

    特征值 $\lambda = \pm i\sqrt{\alpha\gamma}$，纯虚数。$\text{tr}(J^*) = 0$，$\det(J^*) = \alpha\gamma > 0$。

    线性化给出中心（特征值在虚轴上），但线性化不能确定非线性系统中心的稳定性。然而，经典 Lotka-Volterra 系统有守恒量：

    $$H(x,y) = \delta x - \gamma \ln x + \beta y - \alpha \ln y$$

    可以验证 $dH/dt = 0$，因此轨道在 $H$ 的等值线上，是封闭曲线，均衡确实是（非线性）中心。

    $\blacksquare$

!!! definition "定义 70A.6 (竞争 Lotka-Volterra 模型)"
    两个物种竞争同一资源的模型：

    $$\frac{dx_1}{dt} = r_1 x_1 \left(1 - \frac{x_1 + \alpha_{12} x_2}{K_1}\right)$$

    $$\frac{dx_2}{dt} = r_2 x_2 \left(1 - \frac{\alpha_{21} x_1 + x_2}{K_2}\right)$$

    其中 $r_i$ 为内禀增长率，$K_i$ 为承载量，$\alpha_{12}, \alpha_{21}$ 为种间竞争系数。共存均衡 $(x_1^*, x_2^*)$ 的稳定性由 Jacobi 矩阵的特征值决定。

!!! example "例 70A.5"
    设 $\alpha = 1, \beta = 0.1, \gamma = 1.5, \delta = 0.075$。

    非平凡均衡：$x^* = 1.5/0.075 = 20$，$y^* = 1/0.1 = 10$。

    $$J^* = \begin{pmatrix} 0 & -0.1 \times 20 \\ 0.075 \times 10 & 0 \end{pmatrix} = \begin{pmatrix} 0 & -2 \\ 0.75 & 0 \end{pmatrix}$$

    特征值 $\lambda = \pm i\sqrt{1.5} \approx \pm 1.225i$。

    振荡周期 $T = 2\pi / \sqrt{1.5} \approx 5.13$ 时间单位。

---

## 70A.6 多物种群落矩阵与 May 稳定性定理

<div class="context-flow" markdown>

**核心问题**：一个由 $n$ 个物种组成的生态群落，其稳定性与物种多样性之间有什么关系？随机矩阵理论如何回答这个问题？

</div>

1972 年，Robert May 用随机矩阵理论揭示了生态学中一个反直觉的结果：在随机模型中，更复杂的生态群落反而更不稳定。

!!! definition "定义 70A.7 (群落矩阵)"
    一个 $n$ 物种生态群落的动力学在均衡点附近线性化为：

    $$\frac{d\mathbf{x}}{dt} = M \mathbf{x}$$

    其中 $\mathbf{x}$ 为各物种偏离均衡的偏差，$M \in \mathbb{R}^{n \times n}$ 为**群落矩阵**（Community Matrix），即 Jacobi 矩阵在均衡点的值。

    均衡稳定当且仅当 $M$ 的所有特征值实部为负。

!!! definition "定义 70A.8 (May 的随机群落矩阵模型)"
    May 考虑以下随机群落矩阵：

    $$M = -dI + \sigma A$$

    其中：

    - $d > 0$：每个物种的自调节强度（对角线）
    - $A$：$n \times n$ 随机矩阵，每个元素 $A_{ij}$（$i \neq j$）以概率 $C$ 从均值为 0、方差为 1 的分布中抽取，以概率 $1-C$ 为 0
    - $C \in [0,1]$：连接度（物种间交互的概率）
    - $\sigma > 0$：交互强度

!!! theorem "定理 70A.9 (May 稳定性定理)"
    在上述随机模型中，当 $n \to \infty$ 时：

    - 若 $\sigma\sqrt{nC} < d$，则群落几乎必然稳定
    - 若 $\sigma\sqrt{nC} > d$，则群落几乎必然不稳定

    即**稳定性的临界条件**为：

    $$\sigma\sqrt{nC} = d$$

    或等价地，当自调节为单位强度（$d = 1$）时：$\sigma\sqrt{nC} < 1$。

??? proof "证明"
    群落矩阵 $M = -dI + \sigma A$ 稳定当且仅当其所有特征值实部为负，即 $\sigma A$ 的所有特征值实部小于 $d$。

    关键工具是**随机矩阵理论中的圆律**（Circular Law）。对于 $n \times n$ 随机矩阵 $A$，其中非对角元素独立同分布、均值为 0、方差为 $\sigma_0^2 / n$，当 $n \to \infty$ 时，经验特征值分布趋向半径为 $\sigma_0$ 的圆盘上的均匀分布。

    在 May 的模型中，$\sigma A$ 的非零元素以概率 $C$ 出现，方差为 $\sigma^2$。等效方差为 $C\sigma^2$，有效的随机矩阵为 $n$ 维、方差 $C\sigma^2 / n$ 的矩阵乘以 $n$。

    由圆律，$\sigma A$ 的特征值渐近分布在以原点为中心、半径为 $\sigma\sqrt{nC}$ 的圆盘内。

    因此，$\sigma A$ 的最大特征值实部趋向 $\sigma\sqrt{nC}$。

    $M = -dI + \sigma A$ 的特征值 $= -d +$ ($\sigma A$ 的特征值)。所有实部为负的条件为 $\sigma\sqrt{nC} < d$。

    $\blacksquare$

May 定理的生态学含义是深刻的：**更多物种（大 $n$）、更密的连接（大 $C$）或更强的交互（大 $\sigma$）都使群落更不稳定**。这与 Elton（1958）的直觉——复杂生态系统更稳定——形成了尖锐对比，引发了理论生态学关于"多样性-稳定性"关系的持续讨论。

!!! example "例 70A.6"
    一个 $n = 100$ 物种的群落，连接度 $C = 0.1$，自调节强度 $d = 1$。

    临界交互强度：$\sigma_{\text{crit}} = 1/\sqrt{nC} = 1/\sqrt{10} \approx 0.316$。

    若平均交互强度 $\sigma = 0.2 < 0.316$，群落稳定。

    若 $\sigma = 0.5 > 0.316$，群落不稳定。

    增加物种数到 $n = 400$（$\sigma_{\text{crit}} = 1/\sqrt{40} \approx 0.158$），同样的交互强度 $\sigma = 0.2$ 就足以使群落不稳定。

---

## 70A.7 房室模型

<div class="context-flow" markdown>

**核心问题**：药物在体内的吸收、分布和消除过程如何用线性微分方程组来描述？

</div>

!!! definition "定义 70A.9 (房室模型)"
    **房室模型**将生物体分为 $n$ 个房室（如血液、组织、器官等），物质在房室间按一级动力学传输。系统方程为：

    $$\dot{\mathbf{x}}(t) = A\mathbf{x}(t) + \mathbf{u}(t)$$

    其中 $\mathbf{x}(t)$ 为各房室中的物质量，$\mathbf{u}(t)$ 为外部输入，$A$ 为**房室矩阵**，满足：

    - $a_{ij} \geq 0$ 对 $i \neq j$（物质从房室 $j$ 流向房室 $i$ 的速率常数）
    - $a_{ii} \leq 0$（从房室 $i$ 流出的总速率，包括向其他房室的转移和外排）
    - $\sum_i a_{ij} \leq 0$ 对每个 $j$（物质守恒或有外排）

!!! theorem "定理 70A.10 (房室矩阵的特征值)"
    房室矩阵 $A$ 的所有特征值具有非正实部：$\text{Re}(\lambda_k) \leq 0$。若系统是封闭的（$\sum_i a_{ij} = 0$），则 $\lambda = 0$ 是特征值。若系统是开放的（至少存在外排），则所有特征值实部严格为负。

??? proof "证明"
    $A$ 的列和 $\leq 0$（即 $\mathbf{e}^T A \leq \mathbf{0}^T$，其中 $\mathbf{e} = (1,\ldots,1)^T$）。等价地，$-A$ 的列和 $\geq 0$，且 $-A$ 的非对角元素 $\leq 0$，即 $-A$ 是 Z-矩阵。

    由 Gershgorin 圆盘定理，$A$ 的每个特征值 $\lambda$ 满足：

    $$|\lambda - a_{jj}| \leq \sum_{i \neq j} |a_{ij}| = \sum_{i \neq j} a_{ij}$$

    对某个 $j$。由于 $a_{jj} \leq 0$ 且 $|a_{jj}| \geq \sum_{i \neq j} a_{ij}$（因为 $\sum_i a_{ij} \leq 0$ 意味着 $a_{jj} + \sum_{i \neq j} a_{ij} \leq 0$，即 $|a_{jj}| = -a_{jj} \geq \sum_{i \neq j} a_{ij}$），Gershgorin 圆盘的中心在 $a_{jj} \leq 0$ 处，半径不超过 $|a_{jj}|$。因此圆盘完全在左半平面（含虚轴）。

    对于开放系统，至少存在一个 $j$ 使 $\sum_i a_{ij} < 0$（严格外排），则 $|a_{jj}| > \sum_{i \neq j} a_{ij}$，对应的 Gershgorin 圆盘严格在开左半平面内。若系统连通（$A$ 对应的图强连通），则由不可约 Gershgorin 定理，所有特征值都严格在开左半平面内。

    $\blacksquare$

!!! example "例 70A.7"
    **二房室药代动力学模型**。药物从给药部位（房室 1）吸收到血液（房室 2），并从血液中消除。

    $$A = \begin{pmatrix} -k_a & 0 \\ k_a & -k_e \end{pmatrix}$$

    其中 $k_a$ 为吸收速率常数，$k_e$ 为消除速率常数。解为：

    $$\mathbf{x}(t) = e^{At}\mathbf{x}(0)$$

    特征值为 $\lambda_1 = -k_a$，$\lambda_2 = -k_e$。设 $k_a \neq k_e$，解为：

    $$x_1(t) = x_1(0) e^{-k_a t}$$

    $$x_2(t) = \frac{k_a x_1(0)}{k_a - k_e}(e^{-k_e t} - e^{-k_a t}) + x_2(0) e^{-k_e t}$$

    血药浓度 $x_2(t)$ 先上升（吸收阶段）后下降（消除阶段），达到峰值的时间为 $t_{\max} = \frac{\ln(k_a/k_e)}{k_a - k_e}$。

!!! example "例 70A.8"
    **三房室模型**。血浆（房室 1）与快速分布组织（房室 2）和慢速分布组织（房室 3）交换，并从血浆消除。

    $$A = \begin{pmatrix} -(k_{12}+k_{13}+k_e) & k_{21} & k_{31} \\ k_{12} & -k_{21} & 0 \\ k_{13} & 0 & -k_{31} \end{pmatrix}$$

    静脉注射后，血药浓度呈三指数衰减：

    $$C(t) = A_1 e^{-\alpha t} + A_2 e^{-\beta t} + A_3 e^{-\gamma t}$$

    其中 $\alpha > \beta > \gamma > 0$ 是 $-A$ 的三个正特征值。

---

## 70A.8 SIR 模型

<div class="context-flow" markdown>

**核心问题**：最基本的传染病传播模型——SIR 模型——的数学结构和阈值条件是什么？

</div>

SIR 模型是 Kermack 和 McKendrick 于 1927 年提出的，是流行病学中最基本的房室模型。

!!! definition "定义 70A.10 (SIR 模型)"
    SIR 模型将种群分为三类：

    - $S$：易感者（Susceptible）
    - $I$：感染者（Infectious）
    - $R$：康复者（Recovered/Removed）

    动力学方程为：

    $$\frac{dS}{dt} = -\beta SI, \quad \frac{dI}{dt} = \beta SI - \gamma I, \quad \frac{dR}{dt} = \gamma I$$

    其中 $\beta > 0$ 为传播率（单位时间内一个感染者传染一个易感者的概率），$\gamma > 0$ 为恢复率（$1/\gamma$ 为平均感染期）。

!!! theorem "定理 70A.11 (SIR 模型的阈值条件)"
    定义**基本再生数** $R_0 = \beta S_0 / \gamma$，其中 $S_0$ 为初始易感人口比例。

    - 若 $R_0 > 1$：$dI/dt|_{t=0} > 0$，疫情爆发
    - 若 $R_0 \leq 1$：$dI/dt \leq 0$，疫情不会爆发

??? proof "证明"
    在无病均衡附近（$S \approx S_0$，$I \approx 0$），感染者方程线性化为：

    $$\frac{dI}{dt} \approx (\beta S_0 - \gamma) I = \gamma(R_0 - 1) I$$

    初始增长率为 $\gamma(R_0 - 1)$。

    - $R_0 > 1$：增长率为正，$I(t)$ 指数增长（疫情初期）
    - $R_0 < 1$：增长率为负，$I(t)$ 指数衰减
    - $R_0 = 1$：临界情形

    更严格地，在单方程的 SIR 框架中，线性化分析对应于 $1 \times 1$ 的"次代矩阵"$K = (\beta S_0 / \gamma)$，$R_0 = \rho(K) = \beta S_0 / \gamma$。

    $\blacksquare$

!!! theorem "定理 70A.12 (群体免疫阈值)"
    若通过疫苗接种使一定比例 $p$ 的人口免疫，则有效再生数变为 $R_{\text{eff}} = R_0(1-p)$。疫情不会爆发的条件 $R_{\text{eff}} \leq 1$ 等价于：

    $$p \geq 1 - \frac{1}{R_0}$$

    即**群体免疫阈值**为 $p_c = 1 - 1/R_0$。

??? proof "证明"
    疫苗接种使有效易感人口从 $S_0$ 减少到 $(1-p)S_0$。有效再生数：

    $$R_{\text{eff}} = \frac{\beta (1-p) S_0}{\gamma} = R_0 (1-p)$$

    要求 $R_{\text{eff}} \leq 1$：

    $$R_0(1-p) \leq 1 \implies 1 - p \leq 1/R_0 \implies p \geq 1 - 1/R_0$$

    $\blacksquare$

!!! example "例 70A.9"
    几种传染病的 $R_0$ 和群体免疫阈值：

    | 疾病 | $R_0$ | 群体免疫阈值 $p_c = 1 - 1/R_0$ |
    |------|------|------|
    | 麻疹 | 12-18 | 92%-94% |
    | 百日咳 | 12-17 | 92%-94% |
    | 天花 | 5-7 | 80%-86% |
    | COVID-19 (原始株) | 2.5-3.5 | 60%-71% |
    | 流感 | 1.5-2 | 33%-50% |

    $R_0$ 越大，需要免疫的人口比例越高。麻疹的 $R_0$ 极高，这解释了为什么麻疹疫苗接种率需要超过 95% 才能控制传播。

---

## 70A.9 SEIR 模型与次代矩阵方法

<div class="context-flow" markdown>

**核心问题**：对于更复杂的多类型传染病模型，如何系统地计算基本再生数 $R_0$？

</div>

!!! definition "定义 70A.11 (SEIR 模型)"
    SEIR 模型在 SIR 的基础上增加了潜伏期：

    - $S$：易感者（Susceptible）
    - $E$：潜伏者（Exposed，已感染但尚未具传染性）
    - $I$：感染者（Infectious）
    - $R$：康复者（Recovered）

    动力学方程为：

    $$\frac{dS}{dt} = -\beta SI, \quad \frac{dE}{dt} = \beta SI - \sigma E, \quad \frac{dI}{dt} = \sigma E - \gamma I, \quad \frac{dR}{dt} = \gamma I$$

    其中 $\beta$ 为传播率，$\sigma$ 为潜伏期倒数（$1/\sigma$ 为平均潜伏期），$\gamma$ 为恢复率。

!!! definition "定义 70A.12 (次代矩阵与基本再生数)"
    对一般的多类型传染病模型，van den Driessche 和 Watmough（2002）发展了系统的**次代矩阵方法**（Next-Generation Matrix method）。

    在无病均衡点附近，感染子系统的线性化为：

    $$\frac{d\mathbf{z}}{dt} = (F - V)\mathbf{z}$$

    其中 $\mathbf{z}$ 为感染类的偏差向量，

    - $F$：**新感染产生矩阵**——$F_{ij}$ 为状态 $j$ 中的个体单位时间产生的状态 $i$ 类新感染数
    - $V$：**转移矩阵**——描述感染类之间的转移和移除

    $F \geq 0$（非负），$V$ 是 M-矩阵（$V^{-1} \geq 0$）。

    **次代矩阵**定义为：

    $$K = FV^{-1}$$

    **基本再生数**为 $R_0 = \rho(K)$，即 $K$ 的谱半径。

    $K$ 的 $(i,j)$ 元素 $K_{ij}$ 的含义：一个处于感染状态 $j$ 的个体在其整个感染期内产生的处于感染状态 $i$ 的新感染者数。

!!! theorem "定理 70A.13 (基本再生数判据)"
    无病均衡是局部渐近稳定的当且仅当 $R_0 < 1$。当 $R_0 > 1$ 时，疾病可以侵入种群。

??? proof "证明"
    需要证明：$F - V$ 的所有特征值实部为负 $\iff$ $\rho(FV^{-1}) < 1$。

    **$V$ 是非奇异 M-矩阵**：$V$ 的非对角元素非正（转出率以负号出现），且 $V$ 的行和为正（每个感染状态最终会被移除）。因此 $V^{-1} \geq 0$。

    **关键引理**（M-矩阵理论）：设 $F \geq 0$，$V$ 为非奇异 M-矩阵。则 $s(F - V) < 0$ $\iff$ $\rho(FV^{-1}) < 1$，其中 $s(\cdot)$ 表示谱横坐标（spectral abscissa）。

    证明引理：令 $A = F - V$，$M = FV^{-1}$。

    $(\Rightarrow)$ 设 $s(A) < 0$。由于 $V$ 可逆，$A = F - V$ 的特征值与 $FV^{-1} - I$ 的特征值通过相似变换关联（$V^{-1}A = V^{-1}F - I = M - I$ 与 $A$ 相似于 $V^{-1}AV$，但更直接的论证如下）：

    设 $\lambda$ 为 $M = FV^{-1}$ 的特征值，$M\mathbf{u} = \lambda \mathbf{u}$。令 $\mathbf{v} = V^{-1}\mathbf{u}$，则 $F\mathbf{v} = \lambda V \mathbf{v}$，即 $(F - \lambda V)\mathbf{v} = 0$。

    考虑参数化矩阵 $A(\alpha) = F - \alpha V$。当 $\alpha = 1$ 时 $A(1) = F - V = A$。$A(\alpha)$ 的 Perron 根 $r(\alpha)$（作为非负矩阵 $F$ 减去 $\alpha V$ 后的谱横坐标）是 $\alpha$ 的连续递减函数。

    $r(\alpha) = 0$ 当 $\alpha = \rho(FV^{-1})$（因为此时 $F - \alpha V$ 有零特征值意味着 $\alpha^{-1}$ 是 $FV^{-1}$ 的特征值）。

    因此 $s(F - V) < 0 \iff 1 > 1/\rho(FV^{-1})$ 的意义是反的...更精确地：$r$ 在 $\alpha = 0$ 时 $= \rho(F) \geq 0$，在 $\alpha \to \infty$ 时 $\to -\infty$，在 $\alpha^* = \rho(FV^{-1})^{-1}$ 处过零（注意这里参数化需要仔细）。

    实际上，标准结论是：对 $F \geq 0$ 和非奇异 M-矩阵 $V$，$s(F-V) < 0 \iff \rho(FV^{-1}) < 1$。这在 van den Driessche-Watmough (2002) 的定理 2 中有严格证明。

    $\blacksquare$

!!! example "例 70A.10"
    对 SEIR 模型，感染类为 $(E, I)$，在无病均衡 $(S_0, 0, 0, 0)$ 附近：

    $$F = \begin{pmatrix} 0 & \beta S_0 \\ 0 & 0 \end{pmatrix}, \quad V = \begin{pmatrix} \sigma & 0 \\ -\sigma & \gamma \end{pmatrix}$$

    $$V^{-1} = \begin{pmatrix} 1/\sigma & 0 \\ 1/\gamma & 1/\gamma \end{pmatrix}$$

    $$K = FV^{-1} = \begin{pmatrix} 0 & \beta S_0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1/\sigma & 0 \\ 1/\gamma & 1/\gamma \end{pmatrix} = \begin{pmatrix} \beta S_0/\gamma & \beta S_0/\gamma \\ 0 & 0 \end{pmatrix}$$

    $R_0 = \rho(K) = \beta S_0 / \gamma$。

    这与 SIR 模型的 $R_0$ 相同——潜伏期 $1/\sigma$ 不影响 $R_0$（它影响疫情的时间尺度和感染高峰的延迟，但不影响最终是否传播）。

!!! example "例 70A.11"
    **多类型传染病模型**。考虑两种传播途径（如直接接触和环境传播）的模型，感染类为 $(I_1, I_2)$（两种感染状态）：

    $$F = \begin{pmatrix} \beta_{11} S_0 & \beta_{12} S_0 \\ \beta_{21} S_0 & \beta_{22} S_0 \end{pmatrix}, \quad V = \begin{pmatrix} \gamma_1 & 0 \\ 0 & \gamma_2 \end{pmatrix}$$

    $$K = FV^{-1} = S_0 \begin{pmatrix} \beta_{11}/\gamma_1 & \beta_{12}/\gamma_2 \\ \beta_{21}/\gamma_1 & \beta_{22}/\gamma_2 \end{pmatrix}$$

    $R_0 = \rho(K)$ 是这个 $2 \times 2$ 矩阵的谱半径。当 $K$ 不可约时，$R_0$ 不等于任何单一途径的再生数，而是两种途径的协同效应。

    例如 COVID-19 早期，若 $\beta = 0.3$/天，$\gamma = 0.1$/天，$S_0 \approx 1$，则 $R_0 = 3$，群体免疫阈值 $p_c = 1 - 1/3 \approx 67\%$。

---

## 习题

!!! exercise "习题 70A.1"
    对 Leslie 矩阵

    $$L = \begin{pmatrix} 0 & 0 & 6 & 8 \\ 0.5 & 0 & 0 & 0 \\ 0 & 0.6 & 0 & 0 \\ 0 & 0 & 0.3 & 0 \end{pmatrix}$$

    (a) 写出特征方程并求净繁殖率 $R_0$。

    (b) 判断 $L$ 是否本原。

    (c) 求主特征值 $\lambda_1$（数值解）并解释其生态含义。

!!! exercise "习题 70A.2"
    证明：对本原 Leslie 矩阵，主特征值 $\lambda_1$ 对存活率 $s_k$ 的灵敏度为：

    $$\frac{\partial \lambda_1}{\partial s_k} = \frac{w_{k+1} v_k}{\mathbf{w}^T \mathbf{v}}$$

    其中 $\mathbf{v}$ 和 $\mathbf{w}$ 分别为 $\lambda_1$ 的右和左特征向量。利用 $s_k$ 在 $L$ 中出现在 $(k+1, k)$ 位置。

!!! exercise "习题 70A.3"
    一种濒危海龟有 5 个阶段：卵/幼龟、小幼龟、大幼龟、亚成体、成体。Lefkovitch 矩阵为：

    $$A = \begin{pmatrix} 0 & 0 & 0 & 0 & 127 \\ 0.674 & 0.737 & 0 & 0 & 0 \\ 0 & 0.047 & 0.661 & 0 & 0 \\ 0 & 0 & 0.019 & 0.682 & 0 \\ 0 & 0 & 0 & 0.061 & 0.809 \end{pmatrix}$$

    (a) 求主特征值 $\lambda_1$。种群是增长还是衰退？

    (b) 计算弹性矩阵。哪个生命周期参数对 $\lambda_1$ 影响最大？

    (c) 讨论保护策略应侧重于减少幼龟死亡率还是增加产卵量。

!!! exercise "习题 70A.4"
    对竞争 Lotka-Volterra 模型，求共存均衡点 $(x_1^*, x_2^*)$ 的存在条件，并计算 Jacobi 矩阵的特征值。证明共存均衡稳定当且仅当 $\alpha_{12}\alpha_{21} < 1$。

!!! exercise "习题 70A.5"
    考虑 $n = 50$ 个物种的群落，连接度 $C = 0.2$。

    (a) 根据 May 稳定性定理，临界交互强度 $\sigma_{\text{crit}}$ 是多少？

    (b) 若将群落规模增加到 $n = 200$，但保持 $nC$ 不变（即 $C = 0.05$），$\sigma_{\text{crit}}$ 是否改变？讨论其含义。

!!! exercise "习题 70A.6"
    对三房室模型（例 70A.8），设 $k_{12} = 0.5$，$k_{13} = 0.1$，$k_{21} = 0.3$，$k_{31} = 0.05$，$k_e = 0.2$。

    (a) 写出房室矩阵 $A$ 并验证其列和 $\leq 0$。

    (b) 求 $A$ 的特征值并验证实部为负。

    (c) 画出三指数衰减的定性曲线。

!!! exercise "习题 70A.7"
    一种传染病在两个年龄组（儿童和成人）中传播，接触模式不同。$F$ 和 $V$ 矩阵为：

    $$F = \begin{pmatrix} 3 & 1 \\ 0.5 & 2 \end{pmatrix}, \quad V = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$$

    (a) 计算次代矩阵 $K$ 和基本再生数 $R_0$。

    (b) 若只对儿童接种疫苗（减少 $F$ 的第一列），需要多高的覆盖率才能使 $R_0 < 1$？

!!! exercise "习题 70A.8"
    证明群体免疫阈值公式 $p_c = 1 - 1/R_0$ 对以下推广的 SIR 模型仍然成立：种群出生率 $\mu$，自然死亡率 $\mu$（使总人口恒定），疾病死亡率 $\alpha$：

    $$\frac{dS}{dt} = \mu - \beta SI - \mu S, \quad \frac{dI}{dt} = \beta SI - (\gamma + \mu + \alpha) I, \quad \frac{dR}{dt} = \gamma I - \mu R$$

    先求 $R_0$，然后推导疫苗接种后的有效 $R_0$。

!!! exercise "习题 70A.9"
    设 Leslie 矩阵 $L$ 的主特征值为 $\lambda_1$，对应的左、右特征向量为 $\mathbf{w}_1, \mathbf{v}_1$。证明世代时间 $T$ 可以表示为：

    $$T = \frac{\sum_{i=1}^{n} i \cdot f_i \ell_i \lambda_1^{-i}}{\sum_{i=1}^{n} f_i \ell_i \lambda_1^{-i}} = \frac{\lambda_1}{\lambda_1 - 1} \cdot \frac{\ln R_0}{\ln \lambda_1}$$

    （假设 $R_0 \neq 1$）。提示：对特征方程 $g(\lambda) = 1$ 关于 $\lambda$ 求导。

!!! exercise "习题 70A.10"
    考虑一个 SEIRS 模型（康复者可以再次变为易感者），新增方程 $dR/dt = \gamma I - \delta R$，$dS/dt = -\beta SI + \delta R$。

    (a) 写出感染子系统 $(E, I)$ 的 $F$ 和 $V$ 矩阵。

    (b) 证明 $R_0$ 与 SIR/SEIR 模型相同（$R_0 = \beta S_0 / \gamma$）。解释为什么免疫丧失率 $\delta$ 不影响 $R_0$ 但影响长期动态。
