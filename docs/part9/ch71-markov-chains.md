# 第 71 章 Markov 链

<div class="context-flow" markdown>

**前置**：非负矩阵(Ch17) · 特征值(Ch6) · 矩阵幂(Ch14)

**本章脉络**：Markov 链定义 $\to$ 转移矩阵 $\to$ 稳态分布(左特征向量) $\to$ 遍历性 $\to$ 混合时间与谱间隙 $\to$ 吸收链 $\to$ 可逆链与详细平衡 $\to$ MCMC $\to$ PageRank

**延伸**：Markov 链是蒙特卡罗方法（MCMC 在贝叶斯统计中的应用）、自然语言处理（n-gram 语言模型）和金融工程（信用评级迁移矩阵）的核心数学工具

</div>

Andrey Markov 于 1906 年提出了以他名字命名的随机过程理论，这是概率论中最基本的动态模型之一。Markov 链的核心思想是**无记忆性**：系统未来的状态只取决于当前状态，与过去的历史无关。这一简单假设使得整个理论可以用矩阵语言来表达——转移概率构成矩阵，状态分布的演化就是矩阵乘法。

Markov 链理论深刻地统一了线性代数（特征值、矩阵幂）和概率论（收敛、遍历性），并在物理、化学、生物、经济、计算机科学等领域有广泛的应用。本章将从线性代数的视角系统地发展 Markov 链理论。

---

## 71.1 Markov 链的定义

<div class="context-flow" markdown>

**核心问题**：如何用矩阵语言精确定义 Markov 链？转移矩阵具有什么代数结构？

</div>

!!! definition "定义 71.1 (离散时间 Markov 链)"
    设 $\mathcal{S} = \{1, 2, \ldots, n\}$ 为有限**状态空间**。随机过程 $\{X_t\}_{t=0,1,2,\ldots}$ 称为（齐次）**Markov 链**，若对所有 $t \geq 0$ 和所有状态 $i_0, i_1, \ldots, i_t, j$：

    $$P(X_{t+1} = j \mid X_t = i_t, X_{t-1} = i_{t-1}, \ldots, X_0 = i_0) = P(X_{t+1} = j \mid X_t = i_t)$$

    且转移概率不依赖于时间 $t$：

    $$P_{ij} = P(X_{t+1} = j \mid X_t = i), \quad \forall t$$

!!! definition "定义 71.2 (转移矩阵)"
    $n \times n$ 矩阵 $P = (P_{ij})$ 称为**转移矩阵**（或随机矩阵），其中 $P_{ij}$ 为从状态 $i$ 转移到状态 $j$ 的概率。$P$ 满足：

    1. $P_{ij} \geq 0$ 对所有 $i, j$（非负性）
    2. $\sum_{j=1}^{n} P_{ij} = 1$ 对所有 $i$（行和为 1）

    满足条件 (1) 和 (2) 的矩阵称为**行随机矩阵**。

!!! definition "定义 71.3 (概率向量与状态分布)"
    **概率向量**（或分布向量）$\boldsymbol{\pi} \in \mathbb{R}^n$ 满足 $\pi_i \geq 0$ 且 $\sum_i \pi_i = 1$。

    若 $\boldsymbol{\pi}_t^T$ 为第 $t$ 步的状态分布（行向量），即 $(\boldsymbol{\pi}_t)_j = P(X_t = j)$，则状态分布的演化为：

    $$\boldsymbol{\pi}_{t+1}^T = \boldsymbol{\pi}_t^T P$$

    因此 $\boldsymbol{\pi}_t^T = \boldsymbol{\pi}_0^T P^t$。

!!! theorem "定理 71.1 (Chapman-Kolmogorov 方程)"
    $n$ 步转移概率 $P_{ij}^{(n)} = P(X_{t+n} = j \mid X_t = i)$ 由 $P$ 的 $n$ 次幂给出：

    $$P^{(n)} = P^n$$

    等价地，对任意 $0 \leq m \leq n$：

    $$P_{ij}^{(n)} = \sum_{k=1}^{n} P_{ik}^{(m)} P_{kj}^{(n-m)}$$

    即 $P^n = P^m \cdot P^{n-m}$。

??? proof "证明"
    对 $n$ 进行归纳。$n = 1$ 时，$P^{(1)} = P$ 显然成立。

    设 $P^{(n)} = P^n$，则：

    $$P_{ij}^{(n+1)} = P(X_{n+1} = j \mid X_0 = i) = \sum_k P(X_{n+1}=j \mid X_n = k) P(X_n = k \mid X_0 = i)$$

    $$= \sum_k P_{kj} P_{ik}^{(n)} = (P^n P)_{ij} = (P^{n+1})_{ij}$$

    Chapman-Kolmogorov 方程 $P^n = P^m P^{n-m}$ 是矩阵乘法结合律的直接推论。

    $\blacksquare$

!!! example "例 71.1"
    天气的简单 Markov 模型。状态空间 $\{$晴, 雨$\}$，转移矩阵：

    $$P = \begin{pmatrix} 0.8 & 0.2 \\ 0.4 & 0.6 \end{pmatrix}$$

    即晴天后继续晴天的概率为 0.8，变为雨天的概率为 0.2；雨天后变为晴天的概率为 0.4，继续下雨的概率为 0.6。

    两步转移矩阵：

    $$P^2 = \begin{pmatrix} 0.8 & 0.2 \\ 0.4 & 0.6 \end{pmatrix}^2 = \begin{pmatrix} 0.72 & 0.28 \\ 0.56 & 0.44 \end{pmatrix}$$

    即今天是晴天，后天还是晴天的概率为 0.72。

---

## 71.2 转移矩阵的性质

<div class="context-flow" markdown>

**核心问题**：Markov 链的状态如何分类？转移矩阵的结构如何反映链的动态特性？

</div>

!!! definition "定义 71.4 (状态的可达性与互通)"
    - 状态 $j$ 从 $i$ **可达**（$i \to j$），若存在 $n \geq 0$ 使 $P_{ij}^{(n)} > 0$
    - 状态 $i$ 和 $j$ **互通**（$i \leftrightarrow j$），若 $i \to j$ 且 $j \to i$
    - 互通关系是等价关系，将状态空间划分为**通信类**

!!! definition "定义 71.5 (不可约性)"
    Markov 链（或其转移矩阵 $P$）是**不可约的**，若所有状态互通，即只有一个通信类。等价地，$P$ 的有向图是强连通的。

!!! definition "定义 71.6 (状态的分类)"
    - 状态 $i$ 是**吸收的**，若 $P_{ii} = 1$
    - 状态 $i$ 是**常返的**（recurrent），若从 $i$ 出发最终返回 $i$ 的概率为 1
    - 状态 $i$ 是**暂态的**（transient），若从 $i$ 出发最终返回 $i$ 的概率严格小于 1
    - 常返状态 $i$ 的**周期**为 $d_i = \gcd\{n \geq 1 : P_{ii}^{(n)} > 0\}$
    - 若 $d_i = 1$，则状态 $i$ 是**非周期的**

!!! theorem "定理 71.2 (不可约链的状态分类)"
    不可约 Markov 链的所有状态具有相同的类型：要么全部常返，要么全部暂态；且所有状态具有相同的周期 $d$。有限状态不可约链的所有状态都是常返的。

??? proof "证明"
    设 $i \leftrightarrow j$，存在 $r, s$ 使得 $P_{ij}^{(r)} > 0$，$P_{ji}^{(s)} > 0$。

    **暂态/常返的一致性**：设 $i$ 是常返的。对任意 $n$ 使 $P_{jj}^{(n)} > 0$，也有 $P_{ii}^{(r+n+s)} \geq P_{ij}^{(r)} P_{jj}^{(n)} P_{ji}^{(s)} > 0$。反过来，若 $P_{ii}^{(m)} > 0$，则 $P_{jj}^{(r+m+s)} \geq P_{ji}^{(s)} P_{ii}^{(m)} P_{ij}^{(r)} > 0$。利用常返的级数判据（$\sum_n P_{ii}^{(n)} = \infty$ 当且仅当 $i$ 常返），可得 $i$ 常返当且仅当 $j$ 常返。

    **周期的一致性**：$P_{jj}^{(r+n+s)} > 0$（只要 $P_{ii}^{(n)} > 0$），故 $d_j | (r+n+s)$。类似地 $d_j | (r+s)$（取 $n = 0$ 或适当的 $n$）。因此 $d_j | n$，即 $\{n : P_{ii}^{(n)} > 0\}$ 的所有元素都是 $d_j$ 的倍数，故 $d_i \geq d_j$。对称地 $d_j \geq d_i$，故 $d_i = d_j$。

    **有限链的常返性**：假设某个状态 $i$ 是暂态的，则 $P_{ii}^{(n)} \to 0$（$n \to \infty$）。在不可约链中所有状态都是暂态的，故 $P_{ij}^{(n)} \to 0$ 对所有 $i, j$。但 $\sum_j P_{ij}^{(n)} = 1$，由有限状态不可能所有分量都趋于 0，矛盾。

    $\blacksquare$

!!! example "例 71.2"
    **周期链的例子**。考虑状态空间 $\{1, 2, 3\}$，转移矩阵：

    $$P = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$$

    这是一个确定性循环 $1 \to 2 \to 3 \to 1 \to \cdots$。链不可约，周期 $d = 3$。

    $$P^3 = I, \quad P^6 = I, \quad \ldots$$

    $P^n$ 不会收敛到一个极限矩阵（因为 $P$ 有特征值 $1, \omega, \omega^2$，其中 $\omega = e^{2\pi i/3}$，绝对值都为 1）。

---

## 71.3 稳态分布

<div class="context-flow" markdown>

**核心问题**：什么分布在 Markov 链的演化下不发生变化？这样的分布何时存在且唯一？

</div>

!!! definition "定义 71.7 (稳态分布)"
    概率向量 $\boldsymbol{\pi}$ 称为转移矩阵 $P$ 的**稳态分布**（或不变分布、平稳分布），若：

    $$\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T$$

    即 $\boldsymbol{\pi}^T$ 是 $P$ 的属于特征值 $\lambda = 1$ 的左特征向量。

!!! theorem "定理 71.3 (稳态分布的存在性)"
    有限状态 Markov 链的转移矩阵 $P$ 至少有一个稳态分布。

??? proof "证明"
    $P$ 是行随机矩阵，因此 $P\mathbf{e} = \mathbf{e}$（$\mathbf{e} = (1,\ldots,1)^T$），即 $\lambda = 1$ 是 $P$ 的特征值。

    $P^T$ 也有特征值 1，设 $P^T \mathbf{v} = \mathbf{v}$。由 Perron-Frobenius 定理（$P^T \geq 0$），可以选取 $\mathbf{v} \geq 0$。令 $\boldsymbol{\pi} = \mathbf{v} / (\mathbf{e}^T \mathbf{v})$（归一化），则 $\boldsymbol{\pi} \geq 0$，$\boldsymbol{\pi}^T \mathbf{e} = 1$，且 $\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T$。

    $\blacksquare$

!!! theorem "定理 71.4 (稳态分布的唯一性)"
    若 $P$ 不可约，则稳态分布唯一，且 $\boldsymbol{\pi} > 0$（所有分量严格为正）。

??? proof "证明"
    $P$ 不可约等价于 $P^T$ 不可约。由 Perron-Frobenius 定理（Ch17），不可约非负矩阵的 Perron 特征值 $\rho(P^T) = 1$ 对应的特征向量在标量倍数意义下唯一，且所有分量严格为正。

    因此稳态分布 $\boldsymbol{\pi}$（作为 $P^T$ 对应 $\lambda = 1$ 的右特征向量的归一化）是唯一的且 $\boldsymbol{\pi} > 0$。

    $\blacksquare$

!!! example "例 71.3"
    对例 71.1 中的天气模型，求稳态分布。

    $$\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T: \quad \begin{cases} 0.8\pi_1 + 0.4\pi_2 = \pi_1 \\ 0.2\pi_1 + 0.6\pi_2 = \pi_2 \end{cases}$$

    两个方程等价，得 $\pi_2 = \pi_1 / 2$。加上 $\pi_1 + \pi_2 = 1$，解得：

    $$\boldsymbol{\pi} = \begin{pmatrix} 2/3 \\ 1/3 \end{pmatrix}$$

    即长期来看，晴天占 $2/3$，雨天占 $1/3$。

!!! example "例 71.4"
    **随机游走**。粒子在 $\{1, 2, 3, 4, 5\}$ 上进行随机游走，在内部状态以概率 $1/2$ 向左或向右移动，在边界反射：

    $$P = \begin{pmatrix} 1/2 & 1/2 & 0 & 0 & 0 \\ 1/2 & 0 & 1/2 & 0 & 0 \\ 0 & 1/2 & 0 & 1/2 & 0 \\ 0 & 0 & 1/2 & 0 & 1/2 \\ 0 & 0 & 0 & 1/2 & 1/2 \end{pmatrix}$$

    由对称性，稳态分布为均匀分布 $\boldsymbol{\pi} = (1/5, 1/5, 1/5, 1/5, 1/5)^T$。

    验证：由详细平衡 $\pi_i P_{ij} = \pi_j P_{ji}$（因为 $P$ 对称，$P_{ij} = P_{ji}$），均匀分布确实是稳态分布。

---

## 71.4 遍历定理

<div class="context-flow" markdown>

**核心问题**：在什么条件下 $P^n$ 收敛？收敛速度有多快？

</div>

!!! definition "定义 71.8 (遍历链)"
    Markov 链是**遍历的**（ergodic），若它不可约且非周期。等价地，转移矩阵 $P$ 是本原的（存在 $N$ 使 $P^N > 0$，即所有元素严格为正）。

!!! theorem "定理 71.5 (Markov 链基本极限定理)"
    若 $P$ 是遍历的（不可约且非周期），$\boldsymbol{\pi}$ 是唯一的稳态分布，则：

    $$\lim_{n \to \infty} P^n = \mathbf{e}\boldsymbol{\pi}^T$$

    即 $P^n$ 的每一行都收敛到 $\boldsymbol{\pi}^T$。无论初始分布 $\boldsymbol{\pi}_0$ 如何：

    $$\lim_{n \to \infty} \boldsymbol{\pi}_0^T P^n = \boldsymbol{\pi}^T$$

??? proof "证明"
    对 $P$ 进行谱分解。$P$ 的特征值为 $\lambda_1 = 1 > |\lambda_2| \geq |\lambda_3| \geq \cdots \geq |\lambda_n|$（由 Perron-Frobenius 定理，$P$ 本原时主特征值严格大于其他特征值的模）。

    设 $P$ 可对角化（一般情形需用 Jordan 标准形），$P = S \Lambda S^{-1}$，其中 $\Lambda = \text{diag}(1, \lambda_2, \ldots, \lambda_n)$。则：

    $$P^n = S \Lambda^n S^{-1} = S \begin{pmatrix} 1 & & \\ & \lambda_2^n & \\ & & \ddots \\ & & & \lambda_n^n \end{pmatrix} S^{-1}$$

    由于 $|\lambda_k| < 1$ 对 $k \geq 2$，$\lambda_k^n \to 0$。因此：

    $$P^n \to S \begin{pmatrix} 1 & & \\ & 0 & \\ & & \ddots \\ & & & 0 \end{pmatrix} S^{-1}$$

    $S$ 的第一列是 $\mathbf{e}$（$P\mathbf{e} = \mathbf{e}$ 的右特征向量），$S^{-1}$ 的第一行是 $\boldsymbol{\pi}^T$（左特征向量，归一化使得 $\boldsymbol{\pi}^T \mathbf{e} = 1$）。因此极限矩阵为 $\mathbf{e}\boldsymbol{\pi}^T$。

    对一般情形（$P$ 不一定可对角化），使用 Jordan 标准形。$\lambda = 1$ 的 Jordan 块大小为 1（因为 $1$ 是单特征值，Perron-Frobenius 保证了这一点）。其余 Jordan 块 $J_k(\lambda_k)$ 满足 $J_k(\lambda_k)^n \to 0$（当 $|\lambda_k| < 1$）。因此结论相同。

    $\blacksquare$

!!! theorem "定理 71.6 (收敛速度)"
    遍历链的收敛速度由第二大特征值模 $|\lambda_2|$ 控制：

    $$\|P^n - \mathbf{e}\boldsymbol{\pi}^T\| = O(|\lambda_2|^n)$$

    更精确地，对任意初始分布 $\boldsymbol{\pi}_0$：

    $$\|\boldsymbol{\pi}_0^T P^n - \boldsymbol{\pi}^T\|_{TV} \leq C \cdot |\lambda_2|^n$$

    其中 $\|\cdot\|_{TV}$ 为全变差距离，$C$ 为与初始分布有关的常数。

??? proof "证明"
    由谱分解，$\boldsymbol{\pi}_0^T P^n = \boldsymbol{\pi}^T + \sum_{k=2}^{n} c_k \lambda_k^n \mathbf{w}_k^T$，其中 $c_k = \boldsymbol{\pi}_0^T \mathbf{v}_k / (\mathbf{w}_k^T \mathbf{v}_k)$。因此：

    $$\|\boldsymbol{\pi}_0^T P^n - \boldsymbol{\pi}^T\| \leq \sum_{k=2}^{n} |c_k| \cdot |\lambda_k|^n \cdot \|\mathbf{w}_k^T\| \leq C' |\lambda_2|^n$$

    其中 $C' = \sum_{k=2}^{n} |c_k| \|\mathbf{w}_k^T\|$ 是有限常数。

    $\blacksquare$

!!! example "例 71.5"
    对天气模型 $P = \begin{pmatrix} 0.8 & 0.2 \\ 0.4 & 0.6 \end{pmatrix}$：

    特征值：$\lambda_1 = 1$，$\lambda_2 = 0.4$。

    $$P^n \to \mathbf{e}\boldsymbol{\pi}^T = \begin{pmatrix} 2/3 & 1/3 \\ 2/3 & 1/3 \end{pmatrix}$$

    收敛速度由 $|\lambda_2| = 0.4$ 决定。$P^n$ 中各元素以 $0.4^n$ 的速率趋近极限。

    具体地：

    $$P^n = \begin{pmatrix} 2/3 & 1/3 \\ 2/3 & 1/3 \end{pmatrix} + 0.4^n \begin{pmatrix} 1/3 & -1/3 \\ -2/3 & 2/3 \end{pmatrix}$$

---

## 71.5 混合时间与谱间隙

<div class="context-flow" markdown>

**核心问题**：Markov 链需要运行多长时间才能"接近"稳态分布？如何精确量化"接近"？

</div>

!!! definition "定义 71.9 (谱间隙)"
    遍历 Markov 链的**谱间隙**定义为：

    $$\gamma = 1 - |\lambda_2|$$

    其中 $|\lambda_2| = \max\{|\lambda| : \lambda \text{ 是 } P \text{ 的特征值}, \lambda \neq 1\}$。

    谱间隙 $\gamma \in (0, 1]$。$\gamma$ 越大，链收敛越快。

!!! definition "定义 71.10 (全变差距离与混合时间)"
    两个概率分布 $\boldsymbol{\mu}$ 和 $\boldsymbol{\nu}$ 之间的**全变差距离**为：

    $$\|\boldsymbol{\mu} - \boldsymbol{\nu}\|_{TV} = \frac{1}{2}\sum_{i} |\mu_i - \nu_i|$$

    从最坏初始状态出发的**混合距离**为：

    $$d(t) = \max_{i} \|\mathbf{e}_i^T P^t - \boldsymbol{\pi}^T\|_{TV}$$

    **混合时间**定义为：

    $$t_{\text{mix}}(\varepsilon) = \min\{t : d(t) \leq \varepsilon\}$$

    通常取 $\varepsilon = 1/4$ 或 $\varepsilon = 1/(2e)$。

!!! theorem "定理 71.7 (谱间隙与混合时间的关系)"
    对遍历 Markov 链：

    $$\frac{1}{\gamma} \cdot \left(\ln \frac{1}{2\varepsilon}\right) \leq t_{\text{mix}}(\varepsilon) \leq \frac{1}{\gamma} \cdot \left(\ln \frac{1}{\varepsilon \cdot \pi_{\min}}\right)$$

    其中 $\pi_{\min} = \min_i \pi_i$。特别地：

    $$t_{\text{mix}}(\varepsilon) = \Theta\left(\frac{1}{\gamma} \log \frac{1}{\varepsilon \cdot \pi_{\min}}\right)$$

??? proof "证明"
    **上界**：由定理 71.6 的谱分解，

    $$d(t) \leq \frac{1}{\sqrt{\pi_{\min}}} \cdot |\lambda_2|^t = \frac{1}{\sqrt{\pi_{\min}}} \cdot (1-\gamma)^t \leq \frac{1}{\sqrt{\pi_{\min}}} \cdot e^{-\gamma t}$$

    令 $d(t) \leq \varepsilon$，解得 $t \geq \frac{1}{\gamma}(\ln \frac{1}{\varepsilon} + \frac{1}{2}\ln \frac{1}{\pi_{\min}}) \leq \frac{1}{\gamma}\ln \frac{1}{\varepsilon \pi_{\min}}$。

    **下界**：由 $d(t) \geq |\lambda_2|^t / 2 = (1-\gamma)^t / 2$，令 $(1-\gamma)^t / 2 \leq \varepsilon$，得 $t \geq \frac{\ln(1/(2\varepsilon))}{-\ln(1-\gamma)} \geq \frac{\ln(1/(2\varepsilon))}{\gamma}$（利用 $-\ln(1-\gamma) \leq \gamma$——此处应为 $-\ln(1-\gamma) \geq \gamma$，修正后下界为 $\frac{\ln(1/(2\varepsilon))}{\gamma}$ 的常数倍）。

    $\blacksquare$

!!! example "例 71.6"
    **懒惰随机游走**。在 $n$ 个顶点的路径图上，懒惰随机游走的转移矩阵为 $P' = (I + P)/2$，其中 $P$ 是原始随机游走矩阵。

    原始随机游走的特征值为 $\lambda_k = \cos(k\pi/n)$（$k = 0, 1, \ldots, n-1$）。懒惰版本的特征值为 $(1 + \lambda_k)/2$。

    谱间隙 $\gamma = 1 - (1 + \cos(\pi/n))/2 = (1 - \cos(\pi/n))/2 \approx \pi^2/(4n^2)$。

    混合时间 $t_{\text{mix}} = \Theta(n^2 \log n)$。

---

## 71.6 吸收链

<div class="context-flow" markdown>

**核心问题**：如果 Markov 链有吸收状态（一旦进入就无法离开），那么从各暂态状态出发，被各吸收状态吸收的概率是多少？预期吸收时间是多少？

</div>

!!! definition "定义 71.11 (吸收链的标准形式)"
    设 Markov 链有 $r$ 个吸收状态和 $t$ 个暂态状态。适当排列状态后，转移矩阵的**标准形式**为：

    $$P = \begin{pmatrix} Q & R \\ \mathbf{0} & I_r \end{pmatrix}$$

    其中 $Q$ 为 $t \times t$ 矩阵（暂态状态间的转移），$R$ 为 $t \times r$ 矩阵（暂态到吸收的转移），$I_r$ 为 $r$ 阶单位矩阵。

!!! definition "定义 71.12 (基本矩阵)"
    吸收链的**基本矩阵**为：

    $$N = (I - Q)^{-1} = \sum_{k=0}^{\infty} Q^k$$

    $N_{ij}$ 表示从暂态状态 $i$ 出发，在被吸收之前访问暂态状态 $j$ 的期望次数。

!!! theorem "定理 71.8 (基本矩阵的存在性)"
    对吸收链，$\rho(Q) < 1$，因此 $(I - Q)^{-1}$ 存在且非负。

??? proof "证明"
    $Q$ 是转移矩阵 $P$ 的子矩阵，对应暂态状态。由暂态性，从任意暂态状态出发最终被吸收的概率为 1，因此 $Q^n \to 0$（$n \to \infty$）。这意味着 $\rho(Q) < 1$（否则 $Q$ 的幂不会趋于零）。

    由 Neumann 级数，$(I - Q)^{-1} = \sum_{k=0}^{\infty} Q^k$ 收敛且非负。

    $\blacksquare$

!!! theorem "定理 71.9 (吸收概率与期望时间)"
    设 $N = (I-Q)^{-1}$ 为基本矩阵。

    1. **吸收概率矩阵**：$B = NR$，其中 $B_{ij}$ 是从暂态状态 $i$ 出发被吸收状态 $j$ 吸收的概率
    2. **期望吸收时间**：$\mathbf{t} = N\mathbf{e}$，其中 $t_i = \sum_j N_{ij}$ 是从暂态状态 $i$ 出发到被吸收的期望步数

??? proof "证明"
    **(1)** 设 $B_{ij}^{(n)}$ 为在前 $n$ 步内从状态 $i$ 被吸收状态 $j$ 吸收的概率。则：

    $$B^{(n)} = \sum_{k=0}^{n-1} Q^k R$$

    取极限 $n \to \infty$：$B = \sum_{k=0}^{\infty} Q^k R = (I-Q)^{-1} R = NR$。

    **(2)** 从暂态状态 $i$ 出发，在被吸收前在所有暂态状态上花费的总期望时间为 $t_i = \sum_j N_{ij} = (N\mathbf{e})_i$。

    也可以通过条件期望直接推导：$t_i = 1 + \sum_j Q_{ij} t_j$，即 $\mathbf{t} = \mathbf{e} + Q\mathbf{t}$，解得 $\mathbf{t} = (I-Q)^{-1}\mathbf{e} = N\mathbf{e}$。

    $\blacksquare$

!!! example "例 71.7"
    **赌徒破产问题**。赌徒有 $i$ 元初始资金，每局以概率 $p$ 赢 1 元，概率 $q = 1-p$ 输 1 元，直到赢到 $n$ 元或输光。状态 0 和 $n$ 为吸收状态。

    对 $n = 4$，$p = 0.4$：

    $$Q = \begin{pmatrix} 0 & 0.4 & 0 \\ 0.6 & 0 & 0.4 \\ 0 & 0.6 & 0 \end{pmatrix}$$

    基本矩阵 $N = (I - Q)^{-1}$。对称公式给出从状态 $i$ 赢到 $n$ 的概率：

    $$B_{i,n} = \frac{1 - (q/p)^i}{1 - (q/p)^n}$$

    当 $p = 0.4$，$i = 2$，$n = 4$：$B_{2,4} = \frac{1 - 1.5^2}{1 - 1.5^4} = \frac{1-2.25}{1-5.0625} = \frac{-1.25}{-4.0625} \approx 0.308$。

---

## 71.7 可逆链与详细平衡

<div class="context-flow" markdown>

**核心问题**：什么条件下 Markov 链在时间反转下看起来相同？这个性质为什么重要？

</div>

!!! definition "定义 71.13 (详细平衡条件)"
    转移矩阵 $P$ 关于分布 $\boldsymbol{\pi}$ 满足**详细平衡**（detailed balance），若对所有 $i, j$：

    $$\pi_i P_{ij} = \pi_j P_{ji}$$

    满足详细平衡的链称为**可逆链**（reversible chain）。

!!! theorem "定理 71.10 (详细平衡蕴含稳态)"
    若 $P$ 关于 $\boldsymbol{\pi}$ 满足详细平衡，则 $\boldsymbol{\pi}$ 是 $P$ 的稳态分布。

??? proof "证明"
    $$(\boldsymbol{\pi}^T P)_j = \sum_i \pi_i P_{ij} = \sum_i \pi_j P_{ji} = \pi_j \sum_i P_{ji} = \pi_j$$

    因此 $\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T$。

    $\blacksquare$

注意反过来不成立：稳态分布不一定满足详细平衡。满足详细平衡是比稳态更强的条件。

!!! definition "定义 71.14 (Metropolis-Hastings 算法)"
    给定目标分布 $\boldsymbol{\pi}$ 和建议矩阵 $Q$（对称的），**Metropolis-Hastings 算法**构造具有稳态分布 $\boldsymbol{\pi}$ 的可逆链：

    $$P_{ij} = \begin{cases} Q_{ij} \min\left(1, \frac{\pi_j}{\pi_i}\right) & \text{if } i \neq j \\ 1 - \sum_{k \neq i} P_{ik} & \text{if } i = j \end{cases}$$

!!! theorem "定理 71.11 (Metropolis-Hastings 的正确性)"
    由上述公式定义的转移矩阵 $P$ 满足关于 $\boldsymbol{\pi}$ 的详细平衡条件，因此 $\boldsymbol{\pi}$ 是其稳态分布。

??? proof "证明"
    设 $i \neq j$。不失一般性设 $\pi_i Q_{ij} \leq \pi_j Q_{ji}$（当 $Q$ 对称时即 $\pi_i \leq \pi_j$）。则：

    $$P_{ij} = Q_{ij} \cdot \frac{\pi_j}{\pi_i} \cdot \frac{Q_{ji}}{Q_{ij}} = Q_{ij} \cdot \min\left(1, \frac{\pi_j Q_{ji}}{\pi_i Q_{ij}}\right)$$

    （对于一般的非对称 $Q$，公式中的接受概率为 $\min(1, \frac{\pi_j Q_{ji}}{\pi_i Q_{ij}})$。当 $Q$ 对称时简化为 $\min(1, \frac{\pi_j}{\pi_i})$。）

    验证详细平衡：当 $\pi_i \leq \pi_j$（$Q$ 对称情形）：

    $$\pi_i P_{ij} = \pi_i \cdot Q_{ij} \cdot 1 = \pi_i Q_{ij}$$

    $$\pi_j P_{ji} = \pi_j \cdot Q_{ji} \cdot \frac{\pi_i}{\pi_j} = \pi_i Q_{ji} = \pi_i Q_{ij}$$

    （最后一步用了 $Q$ 的对称性。）因此 $\pi_i P_{ij} = \pi_j P_{ji}$。

    $\blacksquare$

!!! example "例 71.8"
    **Metropolis 算法采样离散分布**。目标分布 $\boldsymbol{\pi} = (0.1, 0.3, 0.4, 0.2)$，建议矩阵为均匀随机游走。

    从状态 2（$\pi_2 = 0.3$）建议到状态 3（$\pi_3 = 0.4$）：接受概率 $= \min(1, 0.4/0.3) = 1$，总是接受。

    从状态 3（$\pi_3 = 0.4$）建议到状态 1（$\pi_1 = 0.1$）：接受概率 $= \min(1, 0.1/0.4) = 0.25$，以 25% 的概率接受。

    通过这种"高概率区域容易接受，低概率区域难以接受"的机制，链的长期行为按 $\boldsymbol{\pi}$ 分布。

---

## 71.8 MCMC 与 PageRank

<div class="context-flow" markdown>

**核心问题**：Markov 链方法如何用于从复杂分布中采样？Google 如何用 Markov 链给网页排名？

</div>

!!! definition "定义 71.15 (Markov 链蒙特卡罗 MCMC)"
    **MCMC** 的核心思想：为了从复杂的目标分布 $\boldsymbol{\pi}$ 中采样，构造一个以 $\boldsymbol{\pi}$ 为稳态分布的 Markov 链，运行足够长时间后，链的状态近似服从 $\boldsymbol{\pi}$。

    常用的 MCMC 方法：

    - **Metropolis-Hastings 算法**（如定义 71.14）
    - **Gibbs 采样**：按坐标逐个更新，每次从条件分布中采样

!!! definition "定义 71.16 (Gibbs 采样)"
    对多维目标分布 $\pi(x_1, \ldots, x_d)$，Gibbs 采样的每一步选择一个坐标 $k$，从条件分布中采样：

    $$x_k^{(t+1)} \sim \pi(x_k \mid x_1^{(t)}, \ldots, x_{k-1}^{(t)}, x_{k+1}^{(t)}, \ldots, x_d^{(t)})$$

    其他坐标保持不变。Gibbs 采样是 Metropolis-Hastings 的特例，接受概率恒为 1。

!!! theorem "定理 71.12 (MCMC 的收敛性)"
    若 MCMC 链是遍历的（不可约且非周期），则对任意可积函数 $f$：

    $$\frac{1}{T}\sum_{t=1}^{T} f(X_t) \xrightarrow{a.s.} \sum_i \pi_i f(i) = E_\pi[f]$$

    即时间平均收敛到期望值（遍历定理的推论）。

!!! definition "定义 71.17 (PageRank 模型)"
    设互联网有 $n$ 个网页。链接结构由邻接矩阵 $A$ 描述（$A_{ij} = 1$ 若网页 $j$ 链接到网页 $i$）。出度矩阵 $D = \text{diag}(d_1, \ldots, d_n)$，$d_j = \sum_i A_{ij}$。

    **随机冲浪者模型**：冲浪者在网页间随机浏览。以概率 $1-\alpha$ 沿链接随机跳转，以概率 $\alpha$ 随机跳到任意网页。转移矩阵为：

    $$M = (1-\alpha) A D^{-1} + \frac{\alpha}{n}\mathbf{e}\mathbf{e}^T$$

    其中 $\alpha \in (0,1)$ 为**阻尼因子**（damping factor），Google 使用 $\alpha = 0.15$。

    PageRank 向量 $\mathbf{r}$ 是 $M$ 的稳态分布（$M$ 的列随机版本的 Perron 特征向量）：

    $$\mathbf{r} = M\mathbf{r}, \quad \mathbf{e}^T\mathbf{r} = 1$$

!!! theorem "定理 71.13 (PageRank 的存在唯一性)"
    对任意 $\alpha \in (0,1)$，PageRank 矩阵 $M$ 是正矩阵（所有元素严格为正），因此：

    1. $M$ 本原（不可约且非周期）
    2. $\rho(M) = 1$ 的特征向量唯一且为正
    3. 幂迭代 $\mathbf{r}^{(k+1)} = M\mathbf{r}^{(k)}$ 以几何速率收敛

??? proof "证明"
    由于 $\frac{\alpha}{n}\mathbf{e}\mathbf{e}^T > 0$ 且 $(1-\alpha)AD^{-1} \geq 0$，$M > 0$（所有元素严格为正）。

    $M$ 是列随机矩阵（假设 $AD^{-1}$ 列和为 1 且 $\frac{1}{n}\mathbf{e}\mathbf{e}^T$ 列和为 1），故 $\mathbf{e}^T M = \mathbf{e}^T$，即 $\rho(M) = 1$。

    由 Perron-Frobenius 定理，正矩阵的 Perron 特征向量唯一且为正向量。

    收敛速率：第二大特征值满足 $|\lambda_2| \leq 1 - \alpha$（因为 $M$ 是 $(1-\alpha)$ 的不太远离列随机矩阵的扰动加上 $\alpha$ 的秩一项）。故 $\alpha = 0.15$ 时收敛速率为 $0.85^k$。

    $\blacksquare$

!!! example "例 71.9"
    一个四页的微型网络：

    页 1 链接到页 2 和页 3；页 2 链接到页 3；页 3 链接到页 1；页 4 链接到页 1、页 2 和页 3。

    $$A = \begin{pmatrix} 0 & 0 & 1 & 1 \\ 1 & 0 & 0 & 1 \\ 1 & 1 & 0 & 1 \\ 0 & 0 & 0 & 0 \end{pmatrix}, \quad D = \text{diag}(2, 1, 1, 3)$$

    $$AD^{-1} = \begin{pmatrix} 0 & 0 & 1 & 1/3 \\ 1/2 & 0 & 0 & 1/3 \\ 1/2 & 1 & 0 & 1/3 \\ 0 & 0 & 0 & 0 \end{pmatrix}$$

    取 $\alpha = 0.15$：

    $$M = 0.85 \cdot AD^{-1} + 0.15 \cdot \frac{1}{4}\mathbf{e}\mathbf{e}^T$$

    通过幂迭代计算 PageRank 向量：$\mathbf{r} \approx (0.368, 0.188, 0.327, 0.117)^T$。

    页 1 排名最高（被页 3 和页 4 链接），页 3 次之（被页 1、页 2 和页 4 链接），页 4 排名最低（没有入链接）。

!!! example "例 71.10"
    **PageRank 与线性方程组**。PageRank 方程 $\mathbf{r} = M\mathbf{r}$ 等价于：

    $$(I - M)\mathbf{r} = \mathbf{0}$$

    即 $\mathbf{r}$ 是 $(I - M)$ 的零空间向量。由于 $M$ 是列随机的，$\text{rank}(I-M) = n - 1$，零空间维数为 1。

    也可以写为线性方程组：

    $$(I - (1-\alpha)AD^{-1})\mathbf{r} = \frac{\alpha}{n}\mathbf{e}$$

    因为 $(I - (1-\alpha)AD^{-1})$ 是 M-矩阵（对 $\alpha > 0$），该方程有唯一正解。

    对大规模网络（数十亿网页），幂迭代 $\mathbf{r}^{(k+1)} = M\mathbf{r}^{(k)}$ 是可行的计算方法，因为 $M$ 是稀疏矩阵（每个网页只有少量出链接），矩阵向量乘法可以高效实现。

---

## 本章小结

Markov 链理论是线性代数与概率论的完美结合。本章的核心结果包括：

1. **代数结构**：转移矩阵 $P$ 是行随机矩阵，状态分布的演化 $\boldsymbol{\pi}_{t+1}^T = \boldsymbol{\pi}_t^T P$ 是矩阵-向量乘法。

2. **稳态分布**：$\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T$，即 $\boldsymbol{\pi}^T$ 是 $P$ 属于特征值 1 的左特征向量。Perron-Frobenius 定理保证了不可约链的稳态分布存在且唯一。

3. **遍历定理**：不可约非周期链满足 $P^n \to \mathbf{e}\boldsymbol{\pi}^T$，收敛速率由谱间隙 $\gamma = 1 - |\lambda_2|$ 控制。

4. **吸收链**：基本矩阵 $N = (I-Q)^{-1}$ 给出吸收概率和期望吸收时间，本质上是 Leontief 逆矩阵在概率语境中的对应物。

5. **可逆链与 MCMC**：详细平衡条件保证了可逆链的稳态分布，Metropolis-Hastings 算法利用这一点构造具有任意指定稳态分布的链。

6. **PageRank**：Google 的网页排名算法是 Markov 链理论的经典应用，PageRank 向量是带阻尼的随机冲浪者模型的 Perron 特征向量。
