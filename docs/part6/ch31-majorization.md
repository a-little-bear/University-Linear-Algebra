# 第 31 章 Majorization 与双随机矩阵

<div class="context-flow" markdown>

**前置**：特征值 (Ch06) · 凸集与凸优化 (Ch25) · 矩阵不等式 (Ch18)

**本章脉络**：优序 (Majorization) 的几何定义 $\to$ Hardy-Littlewood-Pólya 定理 $\to$ 双随机矩阵定义与性质 $\to$ Birkhoff 定理（置换矩阵的凸包） $\to$ Schur-Horn 定理（谱与对角元的纽带） $\to$ 算子单调与 Schur 凸性 $\to$ Robin 不等式 $\to$ 熵的 Majorization 性质 $\to$ 算子级数中的应用

**延伸**：Majorization 是量子信息论中纯态转换、量子信道容量以及统计学中不平等度量（Gini 系数）的统一数学框架；它量化了“分布的均匀程度”

</div>

在线性代数中，我们经常需要比较两个向量的“分散程度”或“混乱程度”。**优序理论**（Majorization）为此提供了一个强有力的数学框架。它不仅连接了矩阵的特征值与对角元素（Schur-Horn 定理），还通过**双随机矩阵**（Doubly Stochastic Matrices）将这种序关系与凸几何联系起来。本章将揭示这种隐藏在矩阵数值分布背后的深刻规律。

---

## 31.1 优序 (Majorization) 的定义

!!! definition "定义 31.1 (向量优序)"
    设 $x, y \in \mathbb{R}^n$，将其分量按非递增顺序排列为 $x_{[1]} \ge x_{(2)} \ge \cdots \ge x_{(n)}$。
    称 $y$ **优于** $x$（或 $x$ 被 $y$ 优序），记作 $x \prec y$，如果：
    1.  对 $k=1, \ldots, n-1$，有 $\sum_{i=1}^k x_{(i)} \le \sum_{i=1}^k y_{(i)}$。
    2.  总和相等：$\sum_{i=1}^n x_i = \sum_{i=1}^n y_i$。

!!! intuition "直观理解"
    $x \prec y$ 意味着 $x$ 比 $y$ “更均匀”或“更不集中”。例如，$(1/n, \ldots, 1/n) \prec x \prec (1, 0, \ldots, 0)$ 对任何和为 1 的非负向量 $x$ 都成立。

---

## 31.2 双随机矩阵与 Birkhoff 定理

!!! definition "定义 31.2 (双随机矩阵)"
    矩阵 $P \in M_n(\mathbb{R})$ 称为**双随机矩阵**，如果其元素非负且每行、每列之和均等于 1。

!!! theorem "定理 31.1 (Birkhoff-von Neumann 定理)"
    $n$ 阶双随机矩阵的全集 $\Omega_n$ 是 $n$ 阶**置换矩阵**（Permutation Matrices）的凸包。
    这意味着任何双随机矩阵都可以表示为置换矩阵的凸线性组合：$P = \sum \alpha_i P_{\sigma_i}$。

!!! theorem "定理 31.2 (Hardy-Littlewood-Pólya)"
    对于 $x, y \in \mathbb{R}^n$， $x \prec y$ 当且仅当存在双随机矩阵 $P$ 使得 $x = Py$。

---

## 31.3 Schur-Horn 定理

!!! theorem "定理 31.3 (Schur-Horn 定理)"
    设 $A$ 是 $n$ 阶 Hermite 矩阵，$\mathbf{d} = (a_{11}, \ldots, a_{nn})$ 为其对角元向量，$\boldsymbol{\lambda} = (\lambda_1, \ldots, \lambda_n)$ 为其特征值向量。则：
    $$\mathbf{d} \prec \boldsymbol{\lambda}$$
    反之，若 $\mathbf{d} \prec \boldsymbol{\lambda}$，则必存在以 $\mathbf{d}$ 为对角元、$\boldsymbol{\lambda}$ 为特征值的 Hermite 矩阵。

---

## 31.4 Schur 凸性

!!! definition "定义 31.3 (Schur 凸函数)"
    函数 $\phi: \mathbb{R}^n \to \mathbb{R}$ 称为 **Schur 凸**的，如果 $x \prec y \Rightarrow \phi(x) \le \phi(y)$。
    **例子**：$\phi(x) = \sum x_i^2$ 和 $\phi(x) = \sum x_i \log x_i$ 都是 Schur 凸的。

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

Majorization 理论确立了分布形态的序结构：


****：它为“混乱”和“平均”提供了严谨的数学刻画，证明了均匀分布是所有分布的“基石”（最小元）。

****：Birkhoff 和 Schur-Horn 定理展示了离散组合结构（置换）如何支撑起连续的凸空间，以及矩阵迹与谱之间的本质约束。

****：通过 Schur 凸性，Majorization 成为了处理熵、能量耗散和量子退相干等物理过程的自然语言。
