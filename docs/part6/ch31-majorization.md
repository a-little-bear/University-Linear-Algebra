# 第 31 章 Majorization 与双随机矩阵

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵不等式(Ch18)

**本章脉络**：向量 majorization 定义 → Hardy-Littlewood-Polya 定理 → Schur 凸函数 → 双随机矩阵 → Birkhoff 定理 → Schur-Horn 定理 → Lidskii 不等式 → 弱 majorization → 对数 majorization

**延伸**：Majorization 是信息论（Shannon 熵的 Schur 凹性）和量子信息（纠缠态的 LOCC 变换）的核心工具；在经济学中 Lorenz 曲线本质上就是 majorization

</div>

Majorization 是一种描述"一个向量比另一个向量更分散"的偏序关系。这一概念将特征值不等式、双随机矩阵理论和凸分析统一在一个优雅的框架之中。Hardy、Littlewood 和 Polya 在 1934 年的经典著作《Inequalities》中系统建立了这一理论，此后 majorization 成为矩阵分析中最深刻也最实用的工具之一。

本章从 majorization 的基本定义出发，经由 Hardy-Littlewood-Polya 定理建立与双随机矩阵的等价关系，进而发展 Schur 凸函数理论，最终到达 Schur-Horn 定理和 Lidskii 特征值不等式。

---

## 31.1 Majorization 的定义与等价刻画

<div class="context-flow" markdown>

**核心问题**：如何精确描述"一个向量比另一个向量更均匀"？

</div>

### 向量的降序重排

给定向量 $\boldsymbol{x} = (x_1, x_2, \ldots, x_n) \in \mathbb{R}^n$，我们用 $x_{[1]} \geq x_{[2]} \geq \cdots \geq x_{[n]}$ 表示 $\boldsymbol{x}$ 的分量按降序排列后的结果。

!!! definition "定义 31.1 (Majorization)"
    设 $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^n$。称 $\boldsymbol{x}$ 被 $\boldsymbol{y}$ **majorize**（记作 $\boldsymbol{x} \prec \boldsymbol{y}$），若满足：

    1. 对所有 $k = 1, 2, \ldots, n-1$，
    $$
    \sum_{i=1}^{k} x_{[i]} \leq \sum_{i=1}^{k} y_{[i]}
    $$
    2. 等式在 $k = n$ 时成立：
    $$
    \sum_{i=1}^{n} x_{[i]} = \sum_{i=1}^{n} y_{[i]}
    $$

    直观地说，$\boldsymbol{x} \prec \boldsymbol{y}$ 意味着 $\boldsymbol{x}$ 的分量比 $\boldsymbol{y}$ 的分量"更均匀"或"更分散"。

!!! note "注"
    Majorization 记号的方向在文献中并不统一。本书采用 $\boldsymbol{x} \prec \boldsymbol{y}$ 表示"$\boldsymbol{x}$ 被 $\boldsymbol{y}$ majorize"，即 $\boldsymbol{y}$ 的分量更集中。部分文献（如 Marshall-Olkin-Arnold 的专著）使用相反的方向。

!!! example "例 31.1"
    设 $\boldsymbol{x} = (3, 3, 3)$ 和 $\boldsymbol{y} = (5, 3, 1)$。

    验证 $\boldsymbol{x} \prec \boldsymbol{y}$：

    - $k=1$：$x_{[1]} = 3 \leq 5 = y_{[1]}$ ✓
    - $k=2$：$x_{[1]} + x_{[2]} = 6 \leq 8 = y_{[1]} + y_{[2]}$ ✓
    - $k=3$：$x_{[1]} + x_{[2]} + x_{[3]} = 9 = 9 = y_{[1]} + y_{[2]} + y_{[3]}$ ✓

    因此 $\boldsymbol{x} \prec \boldsymbol{y}$，即均匀向量 $(3,3,3)$ 被 $(5,3,1)$ majorize。

!!! example "例 31.2"
    对于任意 $\boldsymbol{x} \in \mathbb{R}^n$ 且 $\sum x_i = s$，有

    $$
    \left(\frac{s}{n}, \frac{s}{n}, \ldots, \frac{s}{n}\right) \prec \boldsymbol{x} \prec (s, 0, \ldots, 0)
    $$

    即均匀向量是"最分散"的，而将所有质量集中在一个分量上是"最集中"的。

### 部分和的几何解释

!!! definition "定义 31.2 (Lorenz 曲线)"
    给定 $\boldsymbol{x} \in \mathbb{R}^n$，其 **Lorenz 曲线** 是连接点 $(0, 0)$ 和 $(k/n, \sum_{i=1}^{k} x_{[i]} / \sum_{i=1}^n x_i)$（$k = 1, \ldots, n$）的折线。

    $\boldsymbol{x} \prec \boldsymbol{y}$ 当且仅当 $\boldsymbol{y}$ 的 Lorenz 曲线处处不低于 $\boldsymbol{x}$ 的 Lorenz 曲线（假设分量总和相等且为正）。

!!! note "注"
    Lorenz 曲线在经济学中用来度量收入不平等程度。基尼系数（Gini coefficient）正是 Lorenz 曲线与对角线之间面积的两倍。因此 majorization 理论为不平等度量提供了数学基础。

### Majorization 的等价条件

!!! theorem "定理 31.1 (Majorization 的等价刻画)"
    设 $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^n$。以下条件等价：

    1. $\boldsymbol{x} \prec \boldsymbol{y}$。
    2. 对所有凸函数 $\phi: \mathbb{R} \to \mathbb{R}$，有 $\sum_{i=1}^n \phi(x_i) \leq \sum_{i=1}^n \phi(y_i)$。
    3. 存在双随机矩阵 $D \in \mathbb{R}^{n \times n}$ 使得 $\boldsymbol{x} = D\boldsymbol{y}$。
    4. $\boldsymbol{x}$ 在 $\boldsymbol{y}$ 的所有坐标置换的凸包中，即 $\boldsymbol{x} \in \operatorname{conv}\{P\boldsymbol{y} : P \text{ 是置换矩阵}\}$。

    条件 (1) $\Leftrightarrow$ (2) 是 Schur 凸函数理论的基石，(1) $\Leftrightarrow$ (3) 是 Hardy-Littlewood-Polya 定理，(3) $\Leftrightarrow$ (4) 是 Birkhoff 定理的推论。

---

## 31.2 Hardy-Littlewood-Polya 定理

<div class="context-flow" markdown>

**核心问题**：Majorization 关系 $\boldsymbol{x} \prec \boldsymbol{y}$ 与双随机矩阵之间有什么确切联系？

</div>

### T-变换

!!! definition "定义 31.3 (T-变换)"
    一个 **T-变换**（也称 Pigou-Dalton 变换或 Robin Hood 变换）是形如

    $$
    T = (1 - t)I + tP_{ij}
    $$

    的矩阵，其中 $0 \leq t \leq 1$，$P_{ij}$ 是交换第 $i$ 和第 $j$ 个坐标的置换矩阵。

    T-变换的作用是将第 $i$ 个和第 $j$ 个分量"拉近"：若 $\boldsymbol{z} = T\boldsymbol{y}$，则

    $$
    z_i = (1-t)y_i + ty_j, \quad z_j = ty_i + (1-t)y_j
    $$

    而其余分量不变。

!!! note "注"
    T-变换是最简单的双随机矩阵（除了单位阵和置换矩阵之外）。它的直观含义是"从富人那里取一部分给穷人"，这正是 Robin Hood 变换名称的由来。

### 主定理

!!! theorem "定理 31.2 (Hardy-Littlewood-Polya, 1929)"
    设 $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^n$。则 $\boldsymbol{x} \prec \boldsymbol{y}$ 当且仅当存在双随机矩阵 $D$ 使得 $\boldsymbol{x} = D\boldsymbol{y}$。

    等价地，$\boldsymbol{x} \prec \boldsymbol{y}$ 当且仅当 $\boldsymbol{x}$ 可以由 $\boldsymbol{y}$ 经过有限次 T-变换得到。

??? proof "证明"
    **充分性**（$\boldsymbol{x} = D\boldsymbol{y} \Rightarrow \boldsymbol{x} \prec \boldsymbol{y}$）：

    设 $D$ 是双随机矩阵（每行每列之和均为 1 且元素非负）。不妨设 $\boldsymbol{y}$ 已按降序排列。

    对于 $k \leq n-1$：

    $$
    \sum_{i=1}^k x_{[i]} \leq \sum_{i=1}^k x_i' = \sum_{i=1}^k \sum_{j=1}^n d_{ij} y_j
    $$

    其中 $\boldsymbol{x}' = (x_{\sigma(1)}, \ldots, x_{\sigma(n)})$ 是使前 $k$ 个部分和最大的排列。由双随机矩阵的性质和 $y_{[1]} \geq \cdots \geq y_{[n]}$，可以用重排不等式证明：

    $$
    \sum_{i=1}^k x_{[i]} = \max_{|S|=k} \sum_{i \in S} x_i = \max_{|S|=k} \sum_{i \in S} (D\boldsymbol{y})_i \leq \sum_{i=1}^k y_{[i]}
    $$

    最后一步利用了 Ky Fan 最大化原理：对双随机矩阵 $D$，$\sum_{i \in S} \sum_j d_{ij} y_j$ 在 $|S| = k$ 时不超过 $\sum_{i=1}^k y_{[i]}$。

    而 $k = n$ 时的等式由 $D$ 的行和列和均为 1 直接保证。

    **必要性**（$\boldsymbol{x} \prec \boldsymbol{y} \Rightarrow \boldsymbol{x} = D\boldsymbol{y}$）：

    对 $n$ 用归纳法。$n = 1$ 时 $\boldsymbol{x} = \boldsymbol{y}$ 显然成立。

    设命题对 $n-1$ 成立。设 $\boldsymbol{x} \prec \boldsymbol{y}$，$\boldsymbol{x}$ 和 $\boldsymbol{y}$ 均按降序排列。

    若 $\boldsymbol{x} = \boldsymbol{y}$，取 $D = I$ 即可。否则，存在最小的 $k$ 使得 $x_{[k]} < y_{[k]}$。

    由 $\sum x_{[i]} = \sum y_{[i]}$，必存在 $j > k$ 使得 $x_{[j]} > y_{[j]}$。

    选择适当的 $t \in (0, 1]$，构造 T-变换 $T$ 将 $y_{[k]}$ 和 $y_{[j]}$ 拉近，使得 $\boldsymbol{y}' = T\boldsymbol{y}$ 满足：$y'_{[k]}$ 和 $y'_{[j]}$ 中至少有一个与 $x$ 的对应分量相等。

    这样 $\boldsymbol{x} \prec \boldsymbol{y}' \prec \boldsymbol{y}$，且 $\boldsymbol{y}'$ 与 $\boldsymbol{x}$ 在至少一个分量上相等。删去该相等分量后，降维到 $n-1$ 的情形，由归纳假设完成证明。

    总共至多需要 $n-1$ 次 T-变换。$\blacksquare$

!!! example "例 31.3"
    将 $\boldsymbol{y} = (6, 3, 0)$ 通过 T-变换变成 $\boldsymbol{x} = (4, 3, 2)$。

    **第一步**：选 $k=1, j=3$，取 $t = 1/3$。

    $$
    T_1 = \begin{pmatrix} 2/3 & 0 & 1/3 \\ 0 & 1 & 0 \\ 1/3 & 0 & 2/3 \end{pmatrix}
    $$

    $T_1 \boldsymbol{y} = (4, 3, 2) = \boldsymbol{x}$。一次 T-变换即可完成。

    验证 $T_1$ 是双随机的：每行每列之和均为 1，且所有元素非负。

### Majorization 与凸函数

!!! theorem "定理 31.3 (Schur, 1923)"
    $\boldsymbol{x} \prec \boldsymbol{y}$ 当且仅当对所有凸函数 $\phi: \mathbb{R} \to \mathbb{R}$，

    $$
    \sum_{i=1}^n \phi(x_i) \leq \sum_{i=1}^n \phi(y_i)
    $$

??? proof "证明"
    **充分性**：对每个 $k$，取 $\phi(t) = \max(t - c, 0)$（这是凸函数），令 $c$ 适当变化可以恢复部分和不等式。

    具体地，取 $c = y_{[k+1]}$ 时（假设 $\boldsymbol{y}$ 按降序排列），有

    $$
    \sum_{i=1}^n \max(x_i - c, 0) \leq \sum_{i=1}^n \max(y_i - c, 0) = \sum_{i=1}^{k} (y_{[i]} - c)
    $$

    而

    $$
    \sum_{i=1}^n \max(x_i - c, 0) \geq \sum_{i=1}^k (x_{[i]} - c)
    $$

    （只取前 $k$ 个最大的 $x_i - c$ 的正部分）。因此 $\sum_{i=1}^k x_{[i]} \leq \sum_{i=1}^k y_{[i]}$。

    **必要性**：由 Hardy-Littlewood-Polya 定理，$\boldsymbol{x} = D\boldsymbol{y}$。由 Jensen 不等式，对双随机矩阵 $D$ 的每一行：

    $$
    \phi(x_i) = \phi\left(\sum_j d_{ij} y_j\right) \leq \sum_j d_{ij} \phi(y_j)
    $$

    对 $i$ 求和：

    $$
    \sum_i \phi(x_i) \leq \sum_i \sum_j d_{ij} \phi(y_j) = \sum_j \phi(y_j) \sum_i d_{ij} = \sum_j \phi(y_j)
    $$

    最后一步用了 $D$ 的列和为 1。$\blacksquare$

---

## 31.3 Schur 凸函数与 Schur 凹函数

<div class="context-flow" markdown>

**核心问题**：什么样的多元函数保持 majorization 偏序？

</div>

### 定义

!!! definition "定义 31.4 (Schur 凸函数)"
    函数 $F: \mathbb{R}^n \to \mathbb{R}$ 称为 **Schur 凸**的（Schur-convex），若

    $$
    \boldsymbol{x} \prec \boldsymbol{y} \implies F(\boldsymbol{x}) \leq F(\boldsymbol{y})
    $$

    若不等号反向，则称 $F$ 是 **Schur 凹**的（Schur-concave）。

!!! note "注"
    Schur 凸函数不必是凸函数。例如 $F(\boldsymbol{x}) = \prod x_i$ 在正象限上是 Schur 凹的（由 AM-GM 不等式启发），但它既不是凸函数也不是凹函数。

### Schur-Ostrowski 判据

!!! theorem "定理 31.4 (Schur-Ostrowski 判据)"
    设 $F: \mathbb{R}^n \to \mathbb{R}$ 是对称函数（即对任何置换 $\sigma$，$F(x_{\sigma(1)}, \ldots, x_{\sigma(n)}) = F(x_1, \ldots, x_n)$），且连续可微。则 $F$ 是 Schur 凸的当且仅当对所有 $i \neq j$，

    $$
    (x_i - x_j)\left(\frac{\partial F}{\partial x_i} - \frac{\partial F}{\partial x_j}\right) \geq 0
    $$

??? proof "证明"
    **必要性**：设 $F$ 是 Schur 凸的。对 $i \neq j$ 和小的 $\epsilon > 0$，考虑 T-变换

    $$
    \boldsymbol{x}' = (1-\epsilon) \boldsymbol{x} + \epsilon P_{ij} \boldsymbol{x}
    $$

    则 $\boldsymbol{x}' \prec \boldsymbol{x}$，故 $F(\boldsymbol{x}') \leq F(\boldsymbol{x})$，即

    $$
    F(\boldsymbol{x}) - F(\boldsymbol{x}') \geq 0
    $$

    对 $\epsilon \to 0$ 取极限，$x'_i = x_i - \epsilon(x_i - x_j)$，$x'_j = x_j + \epsilon(x_i - x_j)$，得

    $$
    (x_i - x_j)\left(\frac{\partial F}{\partial x_i} - \frac{\partial F}{\partial x_j}\right) \geq 0
    $$

    **充分性**：由于 $F$ 是对称的，只需验证对降序排列的 $\boldsymbol{x}, \boldsymbol{y}$（$\boldsymbol{x} \prec \boldsymbol{y}$），$F(\boldsymbol{x}) \leq F(\boldsymbol{y})$。由 Hardy-Littlewood-Polya 定理，$\boldsymbol{x}$ 可由有限次 T-变换从 $\boldsymbol{y}$ 得到。因此只需证明对每次 T-变换 $\boldsymbol{z}' = T\boldsymbol{z}$（$\boldsymbol{z}' \prec \boldsymbol{z}$），$F(\boldsymbol{z}') \leq F(\boldsymbol{z})$。

    设 $T$ 作用在第 $i$ 和第 $j$ 个分量上（$z_i \geq z_j$），$z'_i = z_i - t(z_i - z_j)$，$z'_j = z_j + t(z_i - z_j)$。定义 $g(s) = F(\boldsymbol{z}(s))$，其中 $z_i(s) = z_i - s(z_i - z_j)$，$z_j(s) = z_j + s(z_i - z_j)$。则

    $$
    g'(s) = (z_i - z_j)\left(-\frac{\partial F}{\partial x_i} + \frac{\partial F}{\partial x_j}\right) \leq 0
    $$

    最后一个不等式由 Schur-Ostrowski 条件保证（当 $z_i(s) \geq z_j(s)$ 时）。因此 $g$ 单调递减，$F(\boldsymbol{z}') = g(t) \leq g(0) = F(\boldsymbol{z})$。$\blacksquare$

### 经典例子

!!! example "例 31.4 (初等对称函数)"
    第 $k$ 个初等对称函数

    $$
    e_k(\boldsymbol{x}) = \sum_{1 \leq i_1 < \cdots < i_k \leq n} x_{i_1} \cdots x_{i_k}
    $$

    是 **Schur 凹**的（在正象限上）。这可以通过 Schur-Ostrowski 判据验证：

    $$
    \frac{\partial e_k}{\partial x_i} = e_{k-1}(\boldsymbol{x}_{\hat{i}})
    $$

    其中 $\boldsymbol{x}_{\hat{i}}$ 是去掉 $x_i$ 后的向量。当 $x_i > x_j$ 时，$e_{k-1}(\boldsymbol{x}_{\hat{i}}) \leq e_{k-1}(\boldsymbol{x}_{\hat{j}})$（因为 $\boldsymbol{x}_{\hat{j}}$ 包含较大的 $x_i$ 而 $\boldsymbol{x}_{\hat{i}}$ 包含较小的 $x_j$），所以 $(x_i - x_j)(\partial e_k / \partial x_i - \partial e_k / \partial x_j) \leq 0$。

!!! example "例 31.5 (Shannon 熵)"
    Shannon 熵

    $$
    H(\boldsymbol{p}) = -\sum_{i=1}^n p_i \log p_i
    $$

    （定义在概率单纯形 $\{\boldsymbol{p} : p_i \geq 0, \sum p_i = 1\}$ 上）是 **Schur 凹**的。

    验证：$\partial H / \partial p_i = -\log p_i - 1$。当 $p_i > p_j$ 时，$-\log p_i - 1 < -\log p_j - 1$，因此

    $$
    (p_i - p_j)\left(\frac{\partial H}{\partial p_i} - \frac{\partial H}{\partial p_j}\right) = (p_i - p_j)(\log p_j - \log p_i) \leq 0
    $$

    因此 $H$ 是 Schur 凹的。这意味着：$\boldsymbol{p} \prec \boldsymbol{q} \Rightarrow H(\boldsymbol{p}) \geq H(\boldsymbol{q})$。

    更均匀的概率分布具有更高的熵，这与信息论的直觉完全一致。均匀分布 $(1/n, \ldots, 1/n)$ 是所有 $n$ 维概率向量中最被 majorize 的（最分散的），因此具有最大熵 $\log n$。

---

## 31.4 双随机矩阵

<div class="context-flow" markdown>

**核心问题**：双随机矩阵有什么结构？它们构成什么样的几何体？

</div>

### 基本定义与性质

!!! definition "定义 31.5 (双随机矩阵)"
    矩阵 $D = (d_{ij}) \in \mathbb{R}^{n \times n}$ 称为**双随机矩阵**（doubly stochastic matrix），若：

    1. $d_{ij} \geq 0$ 对所有 $i, j$ 成立（非负性）；
    2. $\sum_{j=1}^n d_{ij} = 1$ 对所有 $i$ 成立（行和为 1）；
    3. $\sum_{i=1}^n d_{ij} = 1$ 对所有 $j$ 成立（列和为 1）。

!!! note "注"
    - 每个置换矩阵都是双随机的。
    - 双随机矩阵的乘积仍然是双随机的。
    - 双随机矩阵的特征值的模不超过 1，且 1 一定是特征值（对应特征向量 $(1, 1, \ldots, 1)^T$）。
    - $n \times n$ 双随机矩阵的全体 $\Omega_n$ 是 $\mathbb{R}^{n \times n}$ 中的一个紧凸集。

!!! theorem "定理 31.5 (van der Waerden 猜想 / Egorychev-Falikman, 1981)"
    在所有 $n \times n$ 双随机矩阵中，积和式（permanent）的最小值在 $J_n = (1/n)_{n \times n}$（所有元素均为 $1/n$ 的矩阵）处取得：

    $$
    \operatorname{perm}(D) \geq \operatorname{perm}(J_n) = \frac{n!}{n^n}
    $$

    这是 van der Waerden 在 1926 年提出的猜想，由 Egorychev 和 Falikman 在 1981 年独立证明。

### Birkhoff 定理

!!! theorem "定理 31.6 (Birkhoff, 1946)"
    $n \times n$ 双随机矩阵的集合 $\Omega_n$ 恰好是所有 $n \times n$ 置换矩阵的凸包：

    $$
    \Omega_n = \operatorname{conv}\{P_\sigma : \sigma \in S_n\}
    $$

    即每个双随机矩阵 $D$ 都可以写成置换矩阵的凸组合：

    $$
    D = \sum_{k=1}^m \theta_k P_{\sigma_k}, \quad \theta_k > 0, \quad \sum_{k=1}^m \theta_k = 1
    $$

    且至多需要 $m \leq n^2 - 2n + 2$ 个置换矩阵。

??? proof "证明"
    **置换矩阵是双随机的**显然，因此 $\operatorname{conv}\{P_\sigma\} \subseteq \Omega_n$。

    反方向需要证明每个双随机矩阵是置换矩阵的凸组合。

    **第一步**：证明每个双随机矩阵都有一个"置换支撑"——即存在置换矩阵 $P$ 使得 $P$ 的非零位置是 $D$ 的非零位置的子集。

    这等价于证明：在 $D$ 的支撑图（二部图）中存在完美匹配。由 Hall 定理（婚姻定理）：对于 $\{1, \ldots, n\}$ 的任何子集 $S$，$S$ 的邻居集 $N(S) = \{j : \exists i \in S, d_{ij} > 0\}$ 满足 $|N(S)| \geq |S|$。

    验证 Hall 条件：若 $|N(S)| < |S|$，则

    $$
    |S| = \sum_{i \in S} \sum_{j=1}^n d_{ij} = \sum_{j \in N(S)} \sum_{i \in S} d_{ij} \leq \sum_{j \in N(S)} 1 = |N(S)| < |S|
    $$

    矛盾。因此 Hall 条件满足，存在完美匹配，即存在置换矩阵 $P_{\sigma_1}$ 使得 $d_{i,\sigma_1(i)} > 0$ 对所有 $i$ 成立。

    **第二步**：令 $\theta_1 = \min_i d_{i,\sigma_1(i)} > 0$。构造

    $$
    D' = \frac{D - \theta_1 P_{\sigma_1}}{1 - \theta_1}
    $$

    若 $\theta_1 = 1$，则 $D = P_{\sigma_1}$，完成。否则 $D'$ 仍是双随机矩阵（验证非负性和行列和条件），且 $D'$ 比 $D$ 至少少一个正元素。

    **第三步**：对 $D'$ 递归重复上述过程。由于每步至少消去一个正元素，而双随机矩阵至多有 $n^2$ 个正元素，过程在有限步内终止。

    更精确地，由 Caratheodory 定理，$\Omega_n$ 作为 $\mathbb{R}^{n^2}$ 中受 $2n-1$ 个独立等式约束（$n$ 个行和 + $n$ 个列和 - 1 个冗余）的凸多面体，其维数为 $(n-1)^2$，因此至多需要 $(n-1)^2 + 1 = n^2 - 2n + 2$ 个顶点的凸组合。$\blacksquare$

!!! definition "定义 31.6 (Birkhoff 多面体)"
    $\Omega_n$ 称为 **Birkhoff 多面体**，它是 $\mathbb{R}^{n^2}$ 中的一个凸多面体，具有以下性质：

    - 维数：$(n-1)^2$
    - 顶点：恰好是 $n!$ 个置换矩阵
    - 面数：$n^2$ 个面（对应于 $n^2$ 个不等式 $d_{ij} \geq 0$）

!!! example "例 31.6"
    $\Omega_2$ 是连接 $\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$ 和 $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ 的线段。

    一般的 $2 \times 2$ 双随机矩阵为

    $$
    D = \begin{pmatrix} t & 1-t \\ 1-t & t \end{pmatrix}, \quad t \in [0, 1]
    $$

    这正是两个置换矩阵的凸组合 $tI + (1-t)P_{12}$。

---

## 31.5 Schur-Horn 定理

<div class="context-flow" markdown>

**核心问题**：Hermite 矩阵的对角元素与特征值之间满足什么关系？

</div>

### 定理陈述

!!! theorem "定理 31.7 (Schur, 1923)"
    设 $A \in \mathbb{C}^{n \times n}$ 是 Hermite 矩阵，特征值为 $\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_n$，对角元素为 $a_{11}, a_{22}, \ldots, a_{nn}$。令 $\boldsymbol{d} = (a_{11}, \ldots, a_{nn})$ 和 $\boldsymbol{\lambda} = (\lambda_1, \ldots, \lambda_n)$。则

    $$
    \boldsymbol{d} \prec \boldsymbol{\lambda}
    $$

    即对角元素被特征值 majorize。

??? proof "证明"
    设 $A = U \operatorname{diag}(\lambda_1, \ldots, \lambda_n) U^*$，其中 $U$ 是酉矩阵。则

    $$
    a_{ii} = (Ae_i, e_i) = \sum_{k=1}^n \lambda_k |u_{ki}|^2
    $$

    其中 $u_{ki}$ 是 $U$ 的元素。定义矩阵 $B = (b_{ij})$，其中 $b_{ij} = |u_{ji}|^2$。由于 $U$ 是酉的：

    - $\sum_j |u_{ji}|^2 = 1$（$U$ 的第 $j$ 行的范数平方为 1）：即 $B$ 的行和为 1。
    - $\sum_j |u_{ij}|^2 = 1$（$U$ 的第 $i$ 列的范数平方为 1）：即 $B$ 的列和为 1。
    - $b_{ij} = |u_{ji}|^2 \geq 0$。

    因此 $B$ 是双随机矩阵，且 $\boldsymbol{d} = B\boldsymbol{\lambda}$。由 Hardy-Littlewood-Polya 定理，$\boldsymbol{d} \prec \boldsymbol{\lambda}$。$\blacksquare$

### Horn 的逆定理

!!! theorem "定理 31.8 (Horn, 1954)"
    设 $\boldsymbol{d}, \boldsymbol{\lambda} \in \mathbb{R}^n$ 满足 $\boldsymbol{d} \prec \boldsymbol{\lambda}$。则存在 $n \times n$ Hermite 矩阵 $A$，其特征值为 $\lambda_1, \ldots, \lambda_n$，对角元素为 $d_1, \ldots, d_n$。

??? proof "证明"
    由 Hardy-Littlewood-Polya 定理，$\boldsymbol{d} = D\boldsymbol{\lambda}$，其中 $D$ 是双随机矩阵。由 Birkhoff 定理，$D = \sum_k \theta_k P_{\sigma_k}$。

    构造：从 $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$ 出发，通过一系列 Jacobi 旋转（Givens 旋转）将对角元素调整为 $\boldsymbol{d}$。

    具体地，$\boldsymbol{d} = D\boldsymbol{\lambda}$ 意味着 $\boldsymbol{d}$ 可由 $\boldsymbol{\lambda}$ 经有限次 T-变换得到。每次 T-变换 $\boldsymbol{d}' = T_{ij}(t) \boldsymbol{d}$ 对应一次 Givens 旋转：

    $$
    A' = G_{ij}(\alpha)^* A \, G_{ij}(\alpha)
    $$

    其中 $G_{ij}(\alpha)$ 是在第 $i, j$ 平面上旋转角度 $\alpha$ 的 Givens 矩阵，$\alpha$ 选择使得新对角元素满足要求。

    由于 Givens 旋转是酉变换，$A'$ 与 $A$ 有相同的特征值。$\blacksquare$

!!! note "注"
    Schur-Horn 定理的组合形式是：设 $\mathcal{O}(\boldsymbol{\lambda})$ 为对角为 $\boldsymbol{\lambda}$ 的对角矩阵的酉轨道 $\{U \Lambda U^* : U \text{ 酉}\}$，$\pi: \mathbb{C}^{n \times n} \to \mathbb{R}^n$ 为取对角映射。则

    $$
    \pi(\mathcal{O}(\boldsymbol{\lambda})) = \{\boldsymbol{d} \in \mathbb{R}^n : \boldsymbol{d} \prec \boldsymbol{\lambda}\}
    $$

    即酉轨道在对角上的投影恰好是 majorization 多面体。这一结果在辛几何中有深远的推广（Atiyah-Guillemin-Sternberg 凸性定理）。

!!! example "例 31.7"
    设 $\boldsymbol{\lambda} = (4, 2, 0)$。哪些 $\boldsymbol{d}$ 可以作为特征值为 $4, 2, 0$ 的 Hermite 矩阵的对角？

    答：恰好是满足 $\boldsymbol{d} \prec (4, 2, 0)$ 的所有 $\boldsymbol{d}$，即

    $$
    d_{[1]} \leq 4, \quad d_{[1]} + d_{[2]} \leq 6, \quad d_1 + d_2 + d_3 = 6
    $$

    例如 $(3, 2, 1) \prec (4, 2, 0)$（验证：$3 \leq 4$，$5 \leq 6$，$6 = 6$），所以存在特征值为 $4, 2, 0$ 且对角为 $3, 2, 1$ 的 Hermite 矩阵。

    而 $(5, 1, 0) \not\prec (4, 2, 0)$（因为 $5 > 4$），所以不存在这样的矩阵。

---

## 31.6 特征值不等式与 Majorization

<div class="context-flow" markdown>

**核心问题**：矩阵和（或积）的特征值如何被分量矩阵的特征值控制？

</div>

### Lidskii 定理

!!! theorem "定理 31.9 (Lidskii, 1950)"
    设 $A, B$ 是 $n \times n$ Hermite 矩阵，特征值分别为 $\alpha_1 \geq \cdots \geq \alpha_n$ 和 $\beta_1 \geq \cdots \geq \beta_n$。设 $A + B$ 的特征值为 $\gamma_1 \geq \cdots \geq \gamma_n$。则

    $$
    \boldsymbol{\gamma} - \boldsymbol{\alpha} \prec \boldsymbol{\beta}
    $$

    即 $(\gamma_1 - \alpha_1, \ldots, \gamma_n - \alpha_n)$ 被 $(\beta_1, \ldots, \beta_n)$ majorize。

    等价地，$\boldsymbol{\gamma} \prec \boldsymbol{\alpha} + \boldsymbol{\beta}$（这里 $+$ 是分量相加，且两侧均按降序排列后进行 majorization 比较）。

??? proof "证明"
    **部分和不等式的证明**：对任意 $k$ 个指标 $1 \leq i_1 < i_2 < \cdots < i_k \leq n$，需要证明

    $$
    \sum_{j=1}^k (\gamma_{i_j} - \alpha_{i_j}) \leq \sum_{j=1}^k \beta_j
    $$

    更一般地，这是 Weyl-Lidskii 的结果。利用 Courant-Fischer 极小极大原理：

    $$
    \gamma_i = \max_{\dim V = i} \min_{\boldsymbol{x} \in V, \|\boldsymbol{x}\|=1} \boldsymbol{x}^*(A+B)\boldsymbol{x}
    $$

    设 $A$ 的特征向量为 $\boldsymbol{u}_1, \ldots, \boldsymbol{u}_n$，$A+B$ 的特征向量为 $\boldsymbol{v}_1, \ldots, \boldsymbol{v}_n$。

    对于 $\sum_{j=1}^k \gamma_j - \sum_{j=1}^k \alpha_j$ 的估计，考虑子空间 $V_k = \operatorname{span}\{\boldsymbol{v}_1, \ldots, \boldsymbol{v}_k\}$：

    $$
    \sum_{j=1}^k \gamma_j = \operatorname{tr}(P_{V_k}(A+B)) = \operatorname{tr}(P_{V_k} A) + \operatorname{tr}(P_{V_k} B)
    $$

    由 Ky Fan 最大化原理，$\operatorname{tr}(P_{V_k} A) \geq \sum_{j=1}^k \alpha_j$ 不一定成立，但

    $$
    \operatorname{tr}(P_{V_k} B) \leq \sum_{j=1}^k \beta_j
    $$

    利用更精细的分析（Fan-Pall 不等式），可以完成证明。$\blacksquare$

!!! note "注"
    Lidskii 定理可以等价地表述为：映射

    $$
    (A, B) \mapsto \lambda(A+B) - \lambda(A)
    $$

    的值域被 $\{\boldsymbol{x} : \boldsymbol{x} \prec \lambda(B)\}$ 所控制。这里 $\lambda(\cdot)$ 表示按降序排列的特征值向量。

### Weyl 不等式

!!! theorem "定理 31.10 (Weyl 不等式)"
    设 $A, B$ 是 $n \times n$ Hermite 矩阵。则对所有满足 $i + j - 1 \leq n$ 的 $i, j$：

    $$
    \gamma_{i+j-1} \leq \alpha_i + \beta_j
    $$

    对所有满足 $i + j \leq n + 1$ 的 $i, j$：

    $$
    \gamma_{i+j-1} \geq \alpha_i + \beta_j
    $$

    特别地，取 $j = 1$ 或 $j = n$：

    $$
    \alpha_i + \beta_n \leq \gamma_i \leq \alpha_i + \beta_1
    $$

!!! example "例 31.8"
    设 $A = \operatorname{diag}(5, 3, 1)$ 和 $B = \operatorname{diag}(4, 2, 0)$。

    $A + B = \operatorname{diag}(9, 5, 1)$，特征值 $\gamma = (9, 5, 1)$。

    $\boldsymbol{\gamma} - \boldsymbol{\alpha} = (4, 2, 0)$。验证 $(4, 2, 0) \prec (4, 2, 0)$：显然成立（等号情况，因为 $A$ 和 $B$ 同时对角化）。

    对于非同时对角化的情形，majorization 通常是严格的。

### 对数 Majorization

!!! definition "定义 31.7 (对数 Majorization)"
    设 $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}_+^n$（正分量向量）。称 $\boldsymbol{x}$ 被 $\boldsymbol{y}$ **对数 majorize**（log-majorize），记作 $\boldsymbol{x} \prec_{\log} \boldsymbol{y}$，若

    $$
    \prod_{i=1}^k x_{[i]} \leq \prod_{i=1}^k y_{[i]}, \quad k = 1, \ldots, n-1
    $$

    且 $\prod_{i=1}^n x_i = \prod_{i=1}^n y_i$。

    等价地，$(\log x_{[1]}, \ldots, \log x_{[n]}) \prec (\log y_{[1]}, \ldots, \log y_{[n]})$。

!!! theorem "定理 31.11 (奇异值的对数 Majorization)"
    设 $A, B \in \mathbb{C}^{n \times n}$，奇异值分别为 $\sigma_1(A) \geq \cdots \geq \sigma_n(A)$ 和 $\sigma_1(B) \geq \cdots \geq \sigma_n(B)$。则

    $$
    \sigma(AB) \prec_{\log} \sigma(A) \cdot \sigma(B)
    $$

    其中 $\sigma(A) \cdot \sigma(B) = (\sigma_1(A)\sigma_1(B), \ldots, \sigma_n(A)\sigma_n(B))$。

    即对所有 $k = 1, \ldots, n$：

    $$
    \prod_{i=1}^k \sigma_i(AB) \leq \prod_{i=1}^k \sigma_i(A) \sigma_i(B)
    $$

    $k = n$ 时取等号（$|\det(AB)| = |\det A| \cdot |\det B|$）。

!!! note "注"
    对数 majorization 蕴含通常的 majorization：若 $\boldsymbol{x} \prec_{\log} \boldsymbol{y}$ 且分量均为正，则 $\boldsymbol{x} \prec \boldsymbol{y}$（由 Schur 凸函数理论可证）。但反之不成立。

---

## 31.7 弱 Majorization 与推广

<div class="context-flow" markdown>

**核心问题**：如何推广 majorization 以处理分量总和不相等的情形？

</div>

### 弱 Majorization

!!! definition "定义 31.8 (弱 Majorization)"
    设 $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^n$。

    1. 称 $\boldsymbol{x}$ 被 $\boldsymbol{y}$ **弱 majorize**（weak majorization），记作 $\boldsymbol{x} \prec_w \boldsymbol{y}$，若对所有 $k = 1, \ldots, n$：
    $$
    \sum_{i=1}^k x_{[i]} \leq \sum_{i=1}^k y_{[i]}
    $$
    （不要求 $k = n$ 时取等号）。

    2. 称 $\boldsymbol{x}$ 被 $\boldsymbol{y}$ **弱次 majorize**（weak supermajorization），记作 $\boldsymbol{x} \prec^w \boldsymbol{y}$，若对所有 $k = 1, \ldots, n$：
    $$
    \sum_{i=1}^k x_{(i)} \geq \sum_{i=1}^k y_{(i)}
    $$
    其中 $x_{(1)} \leq x_{(2)} \leq \cdots$ 表示升序排列。

!!! theorem "定理 31.12 (弱 Majorization 的等价刻画)"
    1. $\boldsymbol{x} \prec_w \boldsymbol{y}$ 当且仅当对所有非递减凸函数 $\phi$，$\sum \phi(x_i) \leq \sum \phi(y_i)$。

    2. $\boldsymbol{x} \prec_w \boldsymbol{y}$ 当且仅当存在向量 $\boldsymbol{z} \leq \boldsymbol{x}$（分量逐个不超过）使得 $\boldsymbol{z} \prec \boldsymbol{y}$。

    3. $\boldsymbol{x} \prec_w \boldsymbol{y}$ 当且仅当存在双次随机矩阵 $D$（行和不超过 1，列和等于 1）使得 $\boldsymbol{x} = D\boldsymbol{y}$。

??? proof "证明"
    证明 (1)：

    **必要性**：与定理 31.3 类似，取 $\phi(t) = \max(t - c, 0)$（这是非递减凸函数）即可恢复部分和不等式。

    **充分性**：设 $\boldsymbol{x} \prec_w \boldsymbol{y}$。取 $\epsilon = \sum y_{[i]} - \sum x_{[i]} \geq 0$。构造 $\boldsymbol{x}' = (x_1, \ldots, x_n) + (\epsilon/n, \ldots, \epsilon/n)$（实际上取 $\boldsymbol{x}'$ 使得 $\boldsymbol{x}' \prec \boldsymbol{y}$ 且 $x'_i \geq x_i$）。

    对非递减凸函数 $\phi$，$\phi(x_i) \leq \phi(x'_i)$（因为 $x_i \leq x'_i$ 且 $\phi$ 非递减），而 $\sum \phi(x'_i) \leq \sum \phi(y_i)$（因为 $\boldsymbol{x}' \prec \boldsymbol{y}$ 且 $\phi$ 凸）。$\blacksquare$

### Ky Fan 不等式

!!! theorem "定理 31.13 (Ky Fan 不等式)"
    设 $A, B$ 是 $n \times n$ Hermite 矩阵。则对所有 $k = 1, \ldots, n$：

    $$
    \sum_{i=1}^k \lambda_i(A+B) \leq \sum_{i=1}^k \lambda_i(A) + \sum_{i=1}^k \lambda_i(B)
    $$

    这等价于 $\lambda(A+B) \prec_w \lambda(A) + \lambda(B)$（弱 majorization）。

    结合 $\operatorname{tr}(A+B) = \operatorname{tr}(A) + \operatorname{tr}(B)$（即 $k = n$ 时取等号），实际上有

    $$
    \lambda(A+B) \prec \lambda(A) + \lambda(B)
    $$

    即 Lidskii 定理的另一种表述。

!!! example "例 31.9"
    **Majorization 在量子信息中的应用**：

    在量子信息论中，纯态 $|\psi\rangle_{AB}$ 的纠缠度由约化密度矩阵的特征值（Schmidt 系数的平方）的 majorization 关系刻画。

    Nielsen 定理（1999）指出：纯态 $|\psi\rangle$ 可以通过 LOCC（局部操作与经典通信）转化为 $|\phi\rangle$，当且仅当

    $$
    \lambda(\rho_\psi) \prec \lambda(\rho_\phi)
    $$

    其中 $\rho_\psi = \operatorname{tr}_B(|\psi\rangle\langle\psi|)$ 是约化密度矩阵。

    这意味着纠缠态转换的可能性完全由 majorization 关系决定。

!!! example "例 31.10"
    **Majorization 与不等式的统一**

    许多经典不等式可以作为 majorization 的推论：

    1. **AM-GM 不等式**：由 $(1/n, \ldots, 1/n) \prec \boldsymbol{p}$（$\boldsymbol{p}$ 是概率向量），取 Schur 凹函数 $F(\boldsymbol{x}) = \prod x_i^{1/n}$ 或等价地取 $\phi(t) = -\log t$ 得到。

    2. **幂平均不等式**：$M_r(\boldsymbol{x}) \leq M_s(\boldsymbol{x})$（$r \leq s$）可通过 Schur 凸性论证。

    3. **Hadamard 不等式**：$|\det A| \leq \prod \|a_i\|$ 可通过 Schur-Horn 定理和 Oppenheim 不等式推导。

!!! note "注"
    本章建立的 majorization 理论为后续章节中的矩阵不等式、特征值扰动分析和矩阵函数理论提供了统一的框架。读者应特别注意 Schur-Horn 定理（定理 31.7-31.8）的深刻含义：它将代数问题（特征值与对角元素的关系）化为组合问题（majorization）。

## 练习题

1. **[基础] 验证向量 $\mathbf{x} = (2, 2)$ 是否被 $\mathbf{y} = (3, 1)$ 优超（majorize）。**
   ??? success "参考答案"
       1. 分量和相等：$2+2=4, 3+1=4$。✓
       2. 最大部分和：$x_{[1]}=2, y_{[1]}=3$。由于 $2 \le 3$，条件满足。
       故 $\mathbf{x} \prec \mathbf{y}$。这说明 $(2,2)$ 比 $(3,1)$ 更“均匀”。

2. **[性质] 证明：如果 $\mathbf{x} \prec \mathbf{y}$，那么对于任何常数 $c$，都有 $\mathbf{x} + (c, \dots, c) \prec \mathbf{y} + (c, \dots, c)$。**
   ??? success "参考答案"
       由于给每个分量加相同的常数不改变分量的相对排序，部分和的差值保持不变，总和的增量也相同（均为 $nc$），因此优超关系保持不变。

3. **[凸函数] 利用凸函数刻画，证明对于概率向量 $\mathbf{p}$（即 $p_i \ge 0, \sum p_i = 1$），总有 $(\frac{1}{n}, \dots, \frac{1}{n}) \prec \mathbf{p}$。**
   ??? success "参考答案"
       考虑凸函数 $\phi(t) = t^2$。均匀分布使 $\sum p_i^2$ 达到最小值 $1/n$。实际上，对于任何凸函数 $\phi$，由 Jensen 不等式：$\sum \phi(p_i) \ge n \phi(\frac{1}{n} \sum p_i) = n \phi(1/n)$。这符合 $\mathbf{x} \prec \mathbf{y}$ 的凸函数判定条件。

4. **[双随机矩阵] 证明：两个双随机矩阵的乘积仍然是双随机矩阵。**
   ??? success "参考答案"
       设 $D_1, D_2$ 为双随机矩阵。则 $(D_1 D_2) \mathbf{1} = D_1 (D_2 \mathbf{1}) = D_1 \mathbf{1} = \mathbf{1}$。同理 $\mathbf{1}^T (D_1 D_2) = \mathbf{1}^T$。且非负项相乘相加仍非负。故乘积也是双随机的。

5. **[Birkhoff定理] $3 \times 3$ 的双随机矩阵集合有多少个顶点？**
   ??? success "参考答案"
       顶点是所有的置换矩阵。$3 \times 3$ 矩阵共有 $3! = 6$ 个置换矩阵，即该 Birkhoff 多面体有 6 个顶点。

6. **[Schur-Horn] 一个 Hermite 矩阵的对角线元素全为 1，特征值为 $(2, 1, 0)$。这可能吗？**
   ??? success "参考答案"
       对角线向量 $\mathbf{d} = (1, 1, 1)$，特征值向量 $\boldsymbol{\lambda} = (2, 1, 0)$。
       检查优超关系：$1 \le 2$ (k=1), $1+1 \le 2+1$ (k=2), $1+1+1 = 2+1+0$ (k=3)。
       条件全部满足，由 Schur-Horn 定理可知，这样的矩阵确实存在。

7. **[Schur凸性] 为什么乘积函数 $f(\mathbf{x}) = \prod x_i$ 在正象限上是 Schur 凹（Concave）的？**
   ??? success "参考答案"
       因为 $-\log f(\mathbf{x}) = \sum (-\log x_i)$ 是凸函数的和（$-\log x$ 是凸的），因此是 Schur 凸的。负号反转后，原函数即为 Schur 凹。物理意义：分配越均匀，乘积（体积）越大。

8. **[熵] 证明：如果随机变量的概率分布 $\mathbf{p}$ 被 $\mathbf{q}$ 优超（$\mathbf{p} \prec \mathbf{q}$），则其信息熵满足 $H(\mathbf{p}) \ge H(\mathbf{q})$。**
   ??? success "参考答案"
       信息熵 $H(\mathbf{p}) = -\sum p_i \log p_i$ 是 Schur 凹函数。根据定义，$\mathbf{p} \prec \mathbf{q}$ 意味着 $\mathbf{p}$ 更“均匀”，因此熵更大。这符合热力学第二定律：孤立系统的演化倾向于抹平差异，增加熵。

9. **[弱优超] 什么是弱优超（Weak Majorization）？它在处理什么问题时很有用？**
   ??? success "参考答案"
       弱优超不要求总和相等，只要求部分和不等式成立（$\sum_{i=1}^k x_{[i]} \le \sum_{i=1}^k y_{[i]}$）。这在比较两个不同维度的算子，或者处理不包含所有信息的子系统谱分析时非常有用。

10. **[爱因斯坦思考题] 爱因斯坦在研究布朗运动时观察到粒子倾向于从高浓度区域向低浓度区域扩散。如果我们将浓度分布看作向量，这种扩散过程在 Majorization 理论中对应什么操作？这说明了自然界演化的什么逻辑？**
    ??? success "参考答案"
        扩散过程对应于对浓度向量施加一系列的 **T-变换（Robin Hood 变换）**，即 $\mathbf{x}_{new} = D \mathbf{x}_{old}$，其中 $D$ 是双随机矩阵。在 Majorization 理论中，这表现为 $\mathbf{x}_{new} \prec \mathbf{x}_{old}$。这说明自然界的自发演化本质上是一个“去极化”过程：系统不断通过内部的相互作用抹除极端分布，向最平庸、最均匀的平衡态回归。优超关系 $\prec$ 实际上刻画了不可逆过程的时间箭头。

## 本章小结

本章研究了刻画向量“分散程度”的深刻偏序关系——Majorization（优超理论），主要内容包括：

1. **定义与几何**：确立了部分和不等式作为衡量分量均匀性的准则，并通过 Lorenz 曲线给出了直观的几何解释。
2. **Hardy-Littlewood-Polya 定理**：揭示了优超关系与双随机矩阵作用的等价性，证明了优超等价于有限次 T-变换（均贫富过程）的复合。
3. **Schur 凸性**：引入了保持优超序的函数类，建立了对称凸函数与优超关系的内在联系，为证明各种代数不等式提供了统一步径。
4. **Birkhoff 定理**：刻画了双随机矩阵空间的极点为置换矩阵，为矩阵分析与凸几何的结合搭建了桥梁。
5. **Schur-Horn 与特征值**：证明了 Hermite 矩阵对角元总是被其特征值优超这一核心结论，将矩阵的代数性质转化为了组合性质。

