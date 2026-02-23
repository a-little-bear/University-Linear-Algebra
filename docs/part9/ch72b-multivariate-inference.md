# 第 72B 章 多元统计推断

<div class="context-flow" markdown>

**前置**：矩阵值分布 (Ch72A) · 正定矩阵 (Ch16) · 广义逆 (Ch33) · 投影 (Ch07)

**本章脉络**：从单一变量到多元联合 $\to$ 多元正态分布 (Multivariate Normal) 的定义与二次型密度 $\to$ 参数估计：均值向量与协方差矩阵的极大似然估计 (MLE) $\to$ 假设检验：Hotelling $T^2$ 统计量（$t$ 检验的矩阵化） $\to$ 威尔克斯 Lambda ($\Lambda$) 分布与似然比检验 $\to$ 判别分析 (LDA) 的线性代数本质 $\to$ 典型相关分析 (CCA) 与广义特征值问题 $\to$ 应用：临床医学的多指标差异判定、心理学测量、质量控制中的多变量控制图

**延伸**：多元统计推断是线性代数在“真理验证”中的最高应用；它将复杂的假设检验转化为对算子谱和矩阵二次型的区间估计，证明了科学结论的可靠性取决于统计算子在高维空间中的几何稳定性，是现代实证研究的数学判官

</div>

在统计学中，当我们同时观测多个指标（如身高、体重、血压）时，单独对每一项进行分析会忽略变量间的相关性。**多元统计推断**（Multivariate Statistical Inference）通过将观测值整合为向量和矩阵，实现了对多维现象的统一分析。利用 **Hotelling $T^2$** 和 **Wilks' Lambda** 等矩阵统计量，我们能在概率意义下判定两个群体是否存在本质差异。本章将介绍这一作为科学研究证据基石的代数推断框架。

---

## 72B.1 多元正态分布及其 MLE

!!! definition "定义 72B.1 (多元正态分布)"
    向量 $\mathbf{x} \in \mathbb{R}^p$ 服从多元正态分布 $N(\boldsymbol{\mu}, \Sigma)$，其概率密度函数为：
    $$f(\mathbf{x}) = \frac{1}{(2\pi)^{p/2} |\Sigma|^{1/2}} \exp \left( -\frac{1}{2} (\mathbf{x}-\boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x}-\boldsymbol{\mu}) \right)$$
    其中指数项是一个**正定二次型**，定义了高维空间中的等概率超椭球面。

---

## 72B.2 Hotelling $T^2$ 统计量

!!! theorem "定理 72B.1 (均值检验)"
    为了检验 $H_0: \boldsymbol{\mu} = \boldsymbol{\mu}_0$，构造 $T^2$ 统计量：
    $$T^2 = n (\bar{\mathbf{x}} - \boldsymbol{\mu}_0)^T S^{-1} (\bar{\mathbf{x}} - \boldsymbol{\mu}_0)$$
    其中 $S$ 是样本协方差矩阵。
    **代数本质**：这是马氏距离（Mahalanobis Distance）的平方，量化了样本中心偏离目标的统计距离。

---

## 72B.3 典型相关分析 (CCA)

!!! technique "技术：广义特征值问题"
    CCA 寻找两组变量 $X, Y$ 之间的最大相关性。
    这等价于求解涉及互协方差矩阵的广义特征值问题：
    $$\Sigma_{XY} \Sigma_{YY}^{-1} \Sigma_{YX} \mathbf{a} = \rho^2 \Sigma_{XX} \mathbf{a}$$
    这展示了线性代数如何通过算子组合挖掘不同维度间的深层一致性。

---

## 练习题

**1. [基础] 设 $\mathbf{x} \sim N(\boldsymbol{\mu}, \Sigma)$。计算线性变换 $\mathbf{y} = A\mathbf{x} + \mathbf{b}$ 的分布。**

??? success "参考答案"
    **利用期望与方差的线性性质：**
    1. $E[\mathbf{y}] = A E[\mathbf{x}] + \mathbf{b} = A\boldsymbol{\mu} + \mathbf{b}$。
    2. $\operatorname{Var}(\mathbf{y}) = A \operatorname{Var}(\mathbf{x}) A^T = A\Sigma A^T$。
    **结论**：$\mathbf{y} \sim N(A\boldsymbol{\mu} + \mathbf{b}, A\Sigma A^T)$。这证明了正态性在仿射变换下是保持的。

**2. [极大似然] 证明多元均值的 MLE $\hat{\boldsymbol{\mu}}$ 是样本均值 $\bar{\mathbf{x}}$。**

??? success "参考答案"
    **矩阵求导：**
    1. 对数似然函数包含项 $-\sum (\mathbf{x}_i - \boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x}_i - \boldsymbol{\mu})$。
    2. 对 $\boldsymbol{\mu}$ 求梯度（利用 Ch47A 二次型公式）：$\sum 2\Sigma^{-1} (\mathbf{x}_i - \boldsymbol{\mu}) = 0$。
    3. 由于 $\Sigma^{-1}$ 非奇异， $\sum \mathbf{x}_i - n\boldsymbol{\mu} = 0$。
    **结论**：$\hat{\boldsymbol{\mu}} = \frac{1}{n} \sum \mathbf{x}_i$。

**3. [计算] 若样本量 $n=100$，变量维数 $p=2$，样本均值差 $\mathbf{d} = (1, 1)^T$，协方差 $S = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$。计算 $T^2$。**

??? success "参考答案"
    **计算步骤：**
    1. $S^{-1} = \frac{1}{0.75} \begin{pmatrix} 1 & -0.5 \\ -0.5 & 1 \end{pmatrix} = \begin{pmatrix} 4/3 & -2/3 \\ -2/3 & 4/3 \end{pmatrix}$。
    2. 计算 $\mathbf{d}^T S^{-1} \mathbf{d} = (1, 1) \begin{pmatrix} 2/3 \\ 2/3 \end{pmatrix} = 4/3$。
    3. $T^2 = 100 \cdot (4/3) \approx 133.3$。
    **结论**：由于 $T^2$ 很大，远超临界值，我们拒绝原假设，判定均值存在显著差异。

**4. [Wishart] 样本协方差矩阵 $S$ 满足什么分布？**

??? success "参考答案"
    **结论：Wishart 分布。**
    具体地， $(n-1)S \sim W_p(n-1, \Sigma)$。
    这是标量卡方分布在矩阵空间的自然推广，是多元推断的抽样分布基础。

**5. [LDA] 线性判别分析（Fisher LDA）寻找的最佳方向 $w$ 满足什么方程？**

??? success "参考答案"
    **结论：广义特征值方程 $S_B \mathbf{w} = \lambda S_W \mathbf{w}$。**
    其中 $S_B$ 是类间散度阵，$S_W$ 是类内散度阵。
    我们的目标是最大化雷莱商 $J(w) = \frac{w^T S_B w}{w^T S_W w}$，这正是线性代数中寻找最大缩放方向的典型问题。

**6. [性质] 证明：多元正态分布的等密度线是超椭球面。**

??? success "参考答案"
    **理由：**
    1. 等密度要求指数项 $(\mathbf{x}-\boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x}-\boldsymbol{\mu}) = c$（常数）。
    2. 由于 $\Sigma$ 是正定的， $\Sigma^{-1}$ 也是正定的。
    3. 在解析几何中，正定二次型等于常数定义的集合正是主轴方向由特征向量确定的**超椭球**。

**7. [独立性] 若协方差矩阵 $\Sigma$ 是对角阵，这意味着什么？**

??? success "参考答案"
    **结论：各变量之间相互独立。**
    对于正态分布，不相关（协方差为 0）等价于独立。在矩阵上，这意味着联合概率密度可以完全分解为各分量密度的乘积。

**8. [计算] 判定：若 $p > n$， $T^2$ 统计量是否能计算？**

??? success "参考答案"
    **结论：不能（直接计算）。**
    **理由**：当样本量小于变量数时，样本协方差矩阵 $S$ 是奇异的（秩 $\le n-1 < p$）。此时 $S^{-1}$ 不存在。
    **对策**：需要使用**广义逆**（见 Ch33）或引入正则化项。

**9. [Wilks] 什么是 Wilks' Lambda ($\Lambda$)？**

??? success "参考答案"
    **定义：**
    $\Lambda = \frac{|S_W|}{|S_W + S_B|}$。
    它是多元方差分析 (MANOVA) 中的核心指标。
    **代数意义**：它是“未解释的残差体积”与“总偏差体积”的比值。$\Lambda$ 越小，说明自变量对组间差异的解释力越强。

**10. [应用] 简述线性代数在“主成分回归”中解决共线性问题的逻辑。**

??? success "参考答案"
    1. 原始特征 $X$ 存在多重共线性（即 $X^T X$ 接近奇异）。
    2. 利用 SVD 提取前 $k$ 个主成分 $Z = U_k \Sigma_k$。
    3. 在 $Z$ 上进行回归。
    **结论**：由于主成分的方向是正交的，新特征矩阵 $Z^T Z = \Sigma_k^2$ 是理想的对角阵，彻底消除了共线性导致的数值不稳定。

## 本章小结

多元统计推断是线性代数在实证科学中的“真理法则”：

1.  **几何的证据**：它将统计差异量化为向量空间中的距离（Hotelling）与体积比（Wilks），为科学观察提供了严密的几何判据。
2.  **算子的透视**：通过广义特征值问题，CCA 与 LDA 揭示了隐藏在海量数据噪声背后的核心关联模态，展示了线性代数强大的特征提取能力。
3.  **分布的完备性**：从多元正态到威沙特矩阵，本章确立了高维随机世界的代数秩序，证明了统计推断的本质是算子性质在概率测度下的延伸。
