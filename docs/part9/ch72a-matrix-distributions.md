# 第 72A 章 矩阵值随机变量与分布

<div class="context-flow" markdown>

**前置**：随机矩阵 (Ch23) · 正定矩阵 (Ch16) · 概率论基础 · 统计学基础

**本章脉络**：矩阵值随机变量的定义 $\to$ 矩阵正态分布 ($MN_{n,p}$) $\to$ 协方差的 Kronecker 结构 $\to$ Wishart 分布（样本协方差的代数模型） $\to$ 逆 Wishart 分布（共轭先验） $\to$ 矩阵变量 T-分布 $\to$ 矩阵 Beta 与 Gamma 分布 $\to$ 应用：多元方差分析 (MANOVA)、贝叶斯多元回归、计量经济学中的结构建模

**延伸**：矩阵分布是多元统计的支柱；它将标量概率分布提升到了高维张量空间，揭示了多个变量在时间与空间维度上的联合波动规律，是处理金融市场与传感器网络数据的数学利器

</div>

在传统的统计学中，我们研究随机变量 $X$ 或随机向量 $\mathbf{x}$。但在处理具有时间序列的多元数据（如 $n$ 个时刻、 $p$ 个指标的股票收益）时，最自然的描述对象是随机矩阵。**矩阵分布**（Matrix Distributions）不仅描述了矩阵元素的整体波动，还通过特定的乘积结构刻画了变量间复杂的协方差关系。本章将介绍这一现代统计学的高级代数语言。

---

## 72A.1 矩阵正态分布

!!! definition "定义 72A.1 (矩阵正态分布)"
    随机矩阵 $X \in \mathbb{R}^{n \times p}$ 满足 **矩阵正态分布** $MN_{n,p}(M, U, V)$，如果其向量化形式满足：
    $$\operatorname{vec}(X) \sim \mathcal{N}(\operatorname{vec}(M), V \otimes U)$$
    - $M$：$n \times p$ 均值矩阵。
    - $U$：$n \times n$ 行协方差矩阵（描述样本间的相关性）。
    - $V$：$p \times p$ 列协方差矩阵（描述特征间的相关性）。

---

## 72A.2 Wishart 分布

!!! definition "定义 72A.2 (Wishart 分布)"
    设 $X_1, \ldots, X_n$ 是来自 $\mathcal{N}(0, \Sigma)$ 的独立样本。则随机矩阵 $S = \sum X_i X_i^T$ 满足 **Wishart 分布**，记作 $S \sim W_p(n, \Sigma)$。
    **地位**：Wishart 分布是多元分析中“样本协方差矩阵”的理论模型，正如 $\chi^2$ 分布是标量方差的模型。

---

## 72A.3 逆 Wishart 分布与贝叶斯

!!! definition "定义 72A.3 (逆 Wishart 分布)"
    若 $S \sim W_p(n, \Sigma)$，则 $S^{-1}$ 满足**逆 Wishart 分布**。
    **应用**：在贝叶斯统计中，它是多元正态分布协方差矩阵的**共轭先验**，极大简化了后验概率的矩阵计算。

---

## 72A.4 矩阵 T-分布

!!! technique "重尾分布"
    矩阵 T-分布是矩阵正态分布与 Wishart 尺度的混合。它比正态分布更鲁棒，能够捕捉金融数据中的“胖尾”现象（即极端事件发生频率高于正态预测）。

---

## 练习题

**1. [基础] 写出矩阵正态分布 $\operatorname{vec}(X)$ 的协方差矩阵。**

??? success "参考答案"
    **解析：**
    根据矩阵正态分布 $MN_{n,p}(M, U, V)$ 的定义，其向量化算子 $\operatorname{vec}(X)$ 满足多维正态分布。
    其对应的协方差矩阵具有特定的 **Kronecker 积** 结构：
    $$\Sigma_{\operatorname{vec}(X)} = V \otimes U$$
    其中 $V$ 是 $p \times p$ 矩阵，描述了列与列（变量间）的相关性；$U$ 是 $n \times n$ 矩阵，描述了行与行（观察值间）的相关性。这种结构反映了行与列相关性的解耦。

**2. [期望] 若 $X \sim MN(M, U, V)$，求 $E[X]$。**

??? success "参考答案"
    **推导：**
    1. 由于 $\operatorname{vec}(E[X]) = E[\operatorname{vec}(X)] = \operatorname{vec}(M)$。
    2. 向量化算子是线性的且是一一映射。
    3. 因此直接得出：$E[X] = M$。均值矩阵 $M$ 直接给出了随机矩阵每个位置的期望值。

**3. [Wishart] 证明：若 $S \sim W_p(n, \Sigma)$，则 $E[S] = n\Sigma$。**

??? success "参考答案"
    **证明过程：**
    1. 根据 Wishart 分布的构造：$S = \sum_{i=1}^n X_i X_i^T$，其中 $X_i \sim \mathcal{N}(0, \Sigma)$。
    2. 利用期望的线性性质：$E[S] = E[\sum X_i X_i^T] = \sum_{i=1}^n E[X_i X_i^T]$。
    3. 由于 $X_i$ 均值为 0，其协方差 $\Sigma = E[X_i X_i^T] - E[X_i]E[X_i^T] = E[X_i X_i^T]$。
    4. 代入得：$E[S] = \sum_{i=1}^n \Sigma = n\Sigma$。
    **物理意义**：样本协方差矩阵（未归一化）的期望是真实协方差的 $n$ 倍。

**4. [性质] 随机矩阵 $S \sim W_p(n, \Sigma)$ 什么时候是奇异的？**

??? success "参考答案"
    **代数分析：**
    1. $S$ 是 $n$ 个秩为 1 的外积阵 $X_i X_i^T$ 之和。
    2. 根据矩阵秩的不等式：$\operatorname{rank}(S) \le \sum \operatorname{rank}(X_i X_i^T) = n$。
    3. 同时，$S$ 的维度是 $p \times p$。
    4. 如果 $n < p$，则 $\operatorname{rank}(S) \le n < p$，意味着矩阵不满秩。
    **结论**：当**样本量 $n$ 小于变量维数 $p$** 时，$S$ 必然是奇异的（不可逆）。

**5. [不变性] 若 $X \sim MN(M, U, V)$，证明线性变换 $AXB$ 仍满足矩阵正态分布。**

??? success "参考答案"
    **推导：**
    1. 考虑向量化形式：$\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$。
    2. 由于 $\operatorname{vec}(X)$ 是正态分布，其线性变换 $(B^T \otimes A) \operatorname{vec}(X)$ 依然是正态分布。
    3. 均值：$(B^T \otimes A) \operatorname{vec}(M) = \operatorname{vec}(AMB)$。
    4. 协方差：$(B^T \otimes A) (V \otimes U) (B^T \otimes A)^T$。
    5. 利用 Kronecker 性质 $(M \otimes N)^T = M^T \otimes N^T$：$= (B^T V B) \otimes (A U A^T)$。
    **结论**：$AXB \sim MN(AMB, AUA^T, B^TVB)$。

**6. [Beta] 什么是矩阵值 Beta 分布？**

??? success "参考答案"
    **定义：**
    设 $S_1 \sim W_p(n_1, \Sigma)$ 和 $S_2 \sim W_p(n_2, \Sigma)$ 是独立的 Wishart 变量。
    构造矩阵 $B = (S_1 + S_2)^{-1/2} S_1 (S_1 + S_2)^{-1/2}$。
    则 $B$ 满足的分布称为 **矩阵值 Beta 分布**。
    **应用**：它在多元假设检验中用于构造似然比统计量，如 Wilks' Lambda 分布。

**7. [计算] 若 $X \in \mathbb{R}^{2 \times 2}$，且 $U=I, V=I, M=0$，求 $P(\|X\|_F^2 > t)$ 的分布类型。**

??? success "参考答案"
    **步骤：**
    1. Frobenius 范数平方 $\|X\|_F^2 = \sum_{i=1}^2 \sum_{j=1}^2 x_{ij}^2$。
    2. 由于 $U=I, V=I, M=0$，所有的 $x_{ij}$ 都是独立同分布的标准正态变量 $\mathcal{N}(0, 1)$。
    3. 共有 $2 \times 2 = 4$ 个独立变量。
    4. 独立正态变量的平方和遵循卡方分布。
    **结论**：该概率遵循**自由度为 4 的 $\chi^2$ 分布**。

**8. [贝叶斯] 为什么称逆 Wishart 是共轭先验？**

??? success "参考答案"
    **统计逻辑：**
    1. 在贝叶斯推断中，如果“似然函数 $P(Data|\Sigma)$”与“先验分布 $P(\Sigma)$”相乘后，得到的“后验分布 $P(\Sigma|Data)$”与先验属于同一类分布，则称其为共轭先验。
    2. 对于多元正态似然，协方差矩阵的逆（精度矩阵）服从 Wishart 分布，而协方差本身服从逆 Wishart。
    3. 这使得计算出的后验参数只需进行简单的矩阵加法（累加样本平方和），极大地简化了高维随机建模。

**9. [关系] 简述矩阵分布与随机矩阵理论 (RMT) 的区别。**

??? success "参考答案"
    **对比分析：**
    - **矩阵分布**：关注的是**精确的统计模型**。给定特征（如均值、协方差），研究矩阵作为整体的概率密度函数。适用于样本量有限的统计推断。
    - **RMT**：关注的是**渐近的普遍规律**。研究当维度 $n \to \infty$ 时，特征值分布的极限形态（如半圆律、MP 律）。它通常假设元素是独立同分布的，而不关注特定的均值偏移。

**10. [应用] 在信号处理中，如何利用 Wishart 分布检测信号？**

??? success "参考答案"
    **方法：**
    1. 采集传感器数据，计算样本协方差矩阵 $S$。
    2. 假设环境中只有纯噪声，则 $S$ 应服从 $W_p(n, \sigma^2 I)$。
    3. 计算 $S$ 的最大特征值 $\lambda_{\max}$。
    4. 根据 Wishart 分布的谱理论（如 Tracy-Widom 分布），计算出现该 $\lambda_{\max}$ 的概率。
    5. **判定**：若观测值远大于理论上限，则判定存在非随机信号（真实目标）。
