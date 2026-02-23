# 第 72B 章 多元统计推断

<div class="context-flow" markdown>

**前置**：矩阵分布 (Ch72A) · 矩阵分析 (Ch14) · 二次型 (Ch09) · 线性方程组 (Ch01)

**本章脉络**：从一元推断到多元推断 $\to$ Hotelling's $T^2$ 分布（一元 t-检验的矩阵化） $\to$ Wilks' Lambda 分布（行列式比例） $\to$ 多元方差分析 (MANOVA) $\to$ 似然比检验 (LRT) $\to$ 协方差矩阵相等性检验 $\to$ 典型相关分析 (CCA) $\to$ 多元回归推断 $\to$ 高维推断 ($p > n$) 挑战 $\to$ 应用：社会科学差异性分析、心理学多指标评估、医学临床对照实验

**延伸**：多元推断是线性代数在决策论中的最高应用；它通过研究特征值的加权和或乘积，判定多维空间中的观察差异究竟是来自真实效应还是随机噪声，是现代实证研究的数理审判官

</div>

在完成了数据的矩阵化描述（Ch72A）后，统计学的核心任务变成了：如何根据观测到的矩阵样本做出严谨的推断。**多元统计推断**（Multivariate Statistical Inference）利用矩阵的迹、行列式和特征值构建检验统计量。它回答了诸如“两组多维数据是否有显著差异”或“两组变量之间是否存在潜在联系”等问题。

---

## 72B.1 Hotelling's $T^2$ 检验

!!! definition "定义 72B.1 (Hotelling's $T^2$ 统计量)"
    为了检验均值向量 $\boldsymbol{\mu}$ 是否等于 $\boldsymbol{\mu}_0$，定义统计量：
    $$T^2 = n (\bar{\mathbf{x}} - \boldsymbol{\mu}_0)^T \mathbf{S}^{-1} (\bar{\mathbf{x}} - \boldsymbol{\mu}_0)$$
    其中 $\mathbf{S}$ 是样本协方差矩阵。
    **代数本质**：这是马氏距离的平方，它利用协方差矩阵的逆来对不同维度的波动进行归一化。

---

## 72B.2 多元方差分析 (MANOVA)

!!! technique "矩阵分解视角"
    在 MANOVA 中，我们将总离差平方和矩阵 $\mathbf{T}$ 分解为组内误差矩阵 $\mathbf{E}$ 和组间效应矩阵 $\mathbf{H}$：
    $$\mathbf{T} = \mathbf{H} + \mathbf{E}$$
    **统计量**：
    - **Wilks' Lambda**：$\Lambda = \frac{\det(\mathbf{E})}{\det(\mathbf{H} + \mathbf{E})}$。其分布由 $p$ 阶特征值的乘积决定。
    - **Pillai 迹**：基于 $\operatorname{tr}(\mathbf{H}(\mathbf{H}+\mathbf{E})^{-1})$。

---

## 72B.3 典型相关分析 (CCA)

!!! definition "定义 72B.2 (CCA)"
    给定两组变量 $\mathbf{x}$ 和 $\mathbf{y}$，寻找投影向量 $\mathbf{a}, \mathbf{b}$ 使得线性组合 $u = \mathbf{a}^T \mathbf{x}$ 与 $v = \mathbf{b}^T \mathbf{y}$ 的相关系数最大。
    **求解**：这等价于求解涉及互协方差矩阵的广义特征值问题。

---

## 72B.4 高维推断的挑战 ($p > n$)

!!! warning "维数灾难"
    当变量数 $p$ 超过样本量 $n$ 时，样本协方差矩阵 $\mathbf{S}$ 是**奇异的**（不可逆），导致 $T^2$ 统计量失效。
    **对策**：使用岭正则化（Shrinkage 估计）或基于投影的非参数方法。

---

## 练习题

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    
       $\Lambda = 1/11 \approx 0.09$。

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

多元统计推断是线性代数逻辑在经验科学中的终极审判：

1.  **距离的代数化**：Hotelling's $T^2$ 证明了通过矩阵求逆，我们可以将复杂的各向异性波动纠正为标准的统计距离，确立了多维差异判定的基准。
2.  **体积的竞争**：Wilks' Lambda 将复杂的组间比较还原为行列式（体积）的博弈，揭示了效应解释力在多维空间中的代数占比。
3.  **相关性的解构**：CCA 展示了如何利用广义特征值问题，从两组看似杂乱的数据流中提取出共振的线性信号，为理解复杂系统的耦合机制提供了数学手术刀。
