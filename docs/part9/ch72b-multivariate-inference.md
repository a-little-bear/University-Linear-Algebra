# 第 72B 章 多元统计推断中的矩阵方法

<div class="context-flow" markdown>

**前置**：奇异值分解(Ch11) · 特征值(Ch6) · 协方差(Ch72A) · 投影(Ch5)

**本章脉络**：极大似然估计的矩阵导数形式 → 线性回归的几何投影视角 → 广义最小二乘(GLS) → 判别分析(LDA) → 典型相关分析(CCA) → 因子分析 → 结构方程模型(SEM)基石

**延伸**：多元统计推断是现代机器学习（从参数回归到深度学习权重优化）的理论根源

</div>

统计推断旨在通过有限样本估计算子的真实参数。线性代数通过正交投影理论解决了最小二乘估计的几何唯一性问题，并利用广义特征值理论确立了多维数据分类与相关性分析的最优准则。

---

## 72B.1 正交投影与估计理论

!!! definition "定义 72B.1 (最小二乘算子)"
    给定观测向量 $\mathbf{y} \in \mathbb{R}^n$ 与设计矩阵 $X \in \mathbb{R}^{n \times p}$。线性回归系数的最小二乘估计为 $\hat{\beta} = (X^T X)^{-1} X^T \mathbf{y}$。在几何上，$\hat{\mathbf{y}} = X\hat{\beta}$ 是 $\mathbf{y}$ 在 $X$ 的列空间上的正交投影。

!!! theorem "定理 72B.1 (判别分析的广义特征值问题)"
    在 Fisher 线性判别分析（LDA）中，最优投影方向 $\mathbf{w}$ 由以下广义特征值问题决定：
    $$S_b \mathbf{w} = \lambda S_w \mathbf{w}$$
    其中 $S_b$ 为类间散度矩阵，$S_w$ 为类内散度矩阵。最大特征值对应的特征向量最大化了类间分离度与类内离散度的比值（Rayleigh 商）。

---

## 练习题

1. **[正规方程] 从残差向量 $\mathbf{e} = \mathbf{y} - X\beta$ 与 $X$ 的列空间正交这一几何条件，推导出最小二乘法的正规方程 $X^T X \beta = X^T \mathbf{y}$。**
   ??? success "参考答案"
       正交性要求 $X^T \mathbf{e} = 0$。代入 $\mathbf{e}$ 得 $X^T (\mathbf{y} - X\beta) = 0 \implies X^T \mathbf{y} - X^T X \beta = 0$。移项即得正规方程。

2. **[投影矩阵] 证明投影矩阵 $H = X(X^T X)^{-1} X^T$ 是对称且幂等的（$H^2 = H$），并说明其特征值只能是 0 或 1。**
   ??? success "参考答案"
       $H^T = (X(X^T X)^{-1} X^T)^T = X((X^T X)^{-1})^T X^T = H$。$H^2 = X(X^T X)^{-1} X^T X(X^T X)^{-1} X^T = X I (X^T X)^{-1} X^T = H$。幂等算子的特征值 $\lambda$ 满足 $\lambda^2 = \lambda$，故 $\lambda \in \{0, 1\}$。

3. **[计算] 给定设计矩阵 $X = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$ 和观测 $\mathbf{y} = \begin{pmatrix} 2 \\ 3 \end{pmatrix}$。求解 $\hat{\beta}$ 并验证预测值 $\hat{\mathbf{y}}$ 是否在 $X$ 的列空间内。**
   ??? success "参考答案"
       $X^T X = \begin{pmatrix} 2 & 3 \\ 3 & 5 \end{pmatrix}$，其逆为 $\begin{pmatrix} 5 & -3 \\ -3 & 2 \end{pmatrix}$。$X^T \mathbf{y} = \begin{pmatrix} 5 \\ 8 \end{pmatrix}$。
       $\hat{\beta} = \begin{pmatrix} 5 & -3 \\ -3 & 2 \end{pmatrix} \begin{pmatrix} 5 \\ 8 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。$\hat{\mathbf{y}} = X \hat{\beta} = [2, 3]^T = \mathbf{y}$（此处因 $X$ 满秩且 $n=p$，投影为单位算子）。

4. **[CCA分析] 描述典型相关分析（CCA）如何将两组随机变量的相关性最大化问题转化为两个协方差矩阵逆乘积的 SVD 问题。**
   ??? success "参考答案"
       目标是最大化 $\operatorname{corr}(\mathbf{a}^T \mathbf{X}, \mathbf{b}^T \mathbf{Y})$。代数上等价于分析算子 $\Sigma_{XX}^{-1/2} \Sigma_{XY} \Sigma_{YY}^{-1/2}$ 的奇异值。其奇异向量对即为实现最大相关的线性组合系数。

5. **[因子分析] 分析因子分析模型 $\Sigma = \Lambda \Lambda^T + \Psi$ 中，低秩阵 $\Lambda \Lambda^T$ 与残差对角阵 $\Psi$ 在信息捕捉上的代数分工。**
   ??? success "参考答案"
       $\Lambda \Lambda^T$ 通过秩-$k$（$k \ll p$）近似捕捉了变量间的共有方差（公共因子），对应于相关性结构；$\Psi$ 刻画了各变量独立的特有方差（噪声）。这是一种随机形式的结构化低秩分解。

6. **[Gauss-Markov] 说明在 Gauss-Markov 定理中，最小二乘估计具有“最小方差”属性与 $X^T X$ 逆矩阵特征值之间的关系。**
   ??? success "参考答案"
       估计量 $\hat{\beta}$ 的协方差为 $\sigma^2 (X^T X)^{-1}$。最小方差意味着在特定方向上，估计的不确定性受限于 $(X^T X)$ 的谱分布。设计矩阵的正交性（$X^T X = I$）能达到方差分布的最优平衡。

7. **[岭回归] 证明岭回归估计量 $\hat{\beta}_{ridge} = (X^T X + \lambda I)^{-1} X^T \mathbf{y}$ 是对原始最小二乘解进行的一种奇异值收缩。**
   ??? success "参考答案"
       利用 SVD 展开 $X = U \Sigma V^T$。原始解包含 $1/\sigma_i$，岭回归将其替换为 $\sigma_i / (\sigma_i^2 + \lambda)$。当 $\sigma_i \to 0$ 时，该项趋于 0 而非无穷大，从而实现了对数值不稳定性的正则化抑制。

8. **[广义最小二乘] 推导当误差协方差矩阵为 $V$ 时，广义最小二乘解 $\hat{\beta} = (X^T V^{-1} X)^{-1} X^T V^{-1} \mathbf{y}$ 的代数形式。**
   ??? success "参考答案"
       通过线性变换 $L^{-1}$（其中 $LL^T = V$）对原始观测进行预处理（Whitening），使误差变为白噪声。在新空间应用标准最小二乘，逆变换回原空间即得 GLS 表达式。

9. **[LDA特征向量] 解释为什么在 LDA 中，我们通常只关注前 $C-1$ 个特征方向（其中 $C$ 为类别数）。**
   ??? success "参考答案"
       类间散度矩阵 $S_b$ 定义为 $C$ 个均值向量的加权外积之和。由于这 $C$ 个向量受总均值约束，其线性无关的秩至多为 $C-1$。因此 $S_b$ 只有 $C-1$ 个非零特征值。

10. **[统计相干性] 分析当设计矩阵 $X$ 的列向量间存在强线性相关（多重共线性）时，$(X^T X)^{-1}$ 元素的爆炸行为。**
    ??? success "参考答案"
        共线性导致 $X^T X$ 趋于奇异，其最小特征值 $\lambda_{min} \to 0$。逆矩阵项包含 $1/\lambda_{min}$，导致估计方差极大化。这在代数上表现为估计参数对观测扰动表现出极高的敏感性。

## 本章小结

本章论述了多元统计推断中核心算法的代数实现：

1. **几何估计**：通过正交投影确立了线性回归的最优判据与计算范式。
2. **维度压缩**：展示了 SVD 与特征值分解在提取核心数据结构（因子分析、CCA）中的应用。
3. **分类与判别**：利用广义特征值理论确立了线性空间划分的最佳超平面准则。
4. **正则化控制**：探讨了矩阵逆运算的稳定性补偿机制，建立了从理论估计到数值计算的鲁棒链路。
