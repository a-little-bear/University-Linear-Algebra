# 第 72A 章 矩阵值分布

<div class="context-flow" markdown>

**前置**：矩阵微积分(Ch47) · Kronecker 积(Ch19) · 随机矩阵(Ch23) · 概率论基础

**本章脉络**：多元正态分布回顾 → 矩阵正态分布 → Wishart 分布（样本协方差阵） → 矩阵 Beta 与 Gamma 分布 → 矩阵变换的 Jacobian → 矩阵 T 分布 → 矩阵分布的特征函数

**延伸**：Wishart 分布是多元统计分析（MANOVA, PCA）的代数基石；矩阵分布刻画了高维估计量的不确定性

</div>

矩阵值分布将随机变量从标量和向量扩展到了矩阵。该领域利用 Kronecker 积描述矩阵内部的相关性结构，将多维数据的统计特性映射为随机线性算子的性质。

---

## 72A.1 矩阵正态与 Wishart 分布

!!! definition "定义 72A.1 (矩阵正态分布)"
    若随机矩阵 $X \in \mathbb{R}^{n \times p}$ 满足其向量化形式遵循：
    $$\operatorname{vec}(X) \sim \mathcal{N}_{np}(\operatorname{vec}(M), V \otimes U)$$
    其中 $U$ 刻画行间相关性，$V$ 刻画列间相关性，则称 $X$ 服从**矩阵正态分布** $\mathcal{MN}_{n,p}(M, U, V)$。

!!! theorem "定理 72A.3 (Wishart 分布与样本协方差)"
    设 $X_1, \dots, X_n$ 为来自 $\mathcal{N}_p(0, \Sigma)$ 的独立同分布样本。则矩阵 $S = \sum_{i=1}^n X_i X_i^T$ 服从参数为 $(n, \Sigma)$ 的 **Wishart 分布** $W_p(n, \Sigma)$。

---

## 练习题

1. **[矩阵正态] 解释为什么矩阵正态分布的协方差用 Kronecker 积 $V \otimes U$ 表示。**
   ??? success "参考答案"
       Kronecker 积 $V \otimes U$ 编码了一种可分离的相关结构：$U$ 代表 $n$ 个观测（行）之间的相关性（如时间相关性），而 $V$ 代表 $p$ 个变量（列）之间的相关性。相比一般的 $np \times np$ 协方差矩阵，这极大减少了参数量。

2. **[Wishart] 证明：若 $S \sim W_p(n, \Sigma)$，当 $\Sigma=I$ 时，其迹 $\operatorname{tr}(S)$ 是独立正态变量的平方和。**
   ??? success "参考答案"
       $\operatorname{tr}(S) = \operatorname{tr}(\sum X_i X_i^T) = \sum X_i^T X_i = \sum_{i,j} X_{ij}^2$。若 $\Sigma=I$，则 $X_{ij}$ 为独立标准正态变量，故 $\operatorname{tr}(S)$ 服从自由度为 $np$ 的卡方分布。

3. **[Jacobian] 计算矩阵线性变换 $Y = AXB$ 的 Jacobian。**
   ??? success "参考答案"
       利用微分 $dY = A(dX)B \implies \operatorname{vec}(dY) = (B^T \otimes A) \operatorname{vec}(dX)$。Jacobian 即为表示矩阵的行列式：$|B^T \otimes A| = (\det B)^n (\det A)^p$。

4. **[行列式期望] 求 $S \sim W_p(n, I)$ 时，$\det(S)$ 的期望值。**
   ??? success "参考答案"
       Wishart 矩阵的行列式与独立卡方变量的乘积有关。$\mathbb{E}[\det S] = \prod_{i=0}^{p-1} (n-i)$。这反映了样本点集在 $p$ 维空间中构成的平行多面体体积的演变。

5. **[逆Wishart] 定义逆 Wishart 分布及其在贝叶斯统计中的作用。**
   ??? success "参考答案"
       若 $S \sim W_p(n, \Sigma)$，则 $S^{-1}$ 服从逆 Wishart 分布。它是多元正态分布协方差矩阵的共轭先验，允许在贝叶斯推断中进行高效的后验参数更新。

6. **[Beta分布] 描述如何通过两个独立的 Wishart 矩阵构造矩阵 Beta 分布。**
   ??? success "参考答案"
       设 $S_1 \sim W_p(n_1, \Sigma)$ 且 $S_2 \sim W_p(n_2, \Sigma)$。则矩阵 $U = (S_1+S_2)^{-1/2} S_1 (S_1+S_2)^{-1/2}$ 服从矩阵 Beta 分布。它将标量卡方变量的比例推广到了矩阵域。

7. **[Bartlett分解] 解释 Wishart 矩阵的 Bartlett 分解。**
   ??? success "参考答案"
       $S \sim W_p(n, I)$ 可以分解为 $S = T T^T$，其中 $T$ 为下三角阵，其对角元 $T_{ii}^2 \sim \chi_{n-i+1}^2$，非对角元 $T_{ij} \sim \mathcal{N}(0, 1)$。这提供了模拟 Wishart 样本的高效算法。

8. **[特征函数] 写出矩阵正态分布特征函数的迹形式。**
   ??? success "参考答案"
       $\phi_X(Z) = \exp(i \operatorname{tr}(Z^T M) - \frac{1}{2} \operatorname{tr}(Z^T U Z V))$。迹中的二次项捕捉了矩阵各元素间聚合的方差结构。

9. **[奇异Wishart] 什么时候 Wishart 矩阵是奇异的？**
   ??? success "参考答案"
       当样本数 $n < p$ 时，Wishart 矩阵 $W_p(n, \Sigma)$ 以概率 1 为奇异矩阵。此时样本协方差矩阵不满秩，在正定锥上不具有相对于 Lebesgue 测度的密度函数。

10. **[熵] 矩阵正态分布的熵如何与 $U$ 和 $V$ 的行列式关联？**
    ??? success "参考答案"
        熵正比于 $\log \det(V \otimes U) = n \log \det V + p \log \det U$。这表明矩阵的信息量（不确定性）是行结构与列结构不确定性的线性叠加。

## 本章小结

本章探讨了矩阵变量的统计分布理论：

1. **结构化不确定性**：利用 Kronecker 积定义了矩阵正态分布，区分了行间与列间相关性。
2. **协方差动力学**：确立了 Wishart 分布作为样本协方差阵的基础模型。
3. **随机几何**：利用矩阵 Jacobian 推导了矩阵变换后的密度函数演变。
4. **贝叶斯共轭**：建立了逆 Wishart 与矩阵 Beta 分布在高维不确定性估算中的纽带。
