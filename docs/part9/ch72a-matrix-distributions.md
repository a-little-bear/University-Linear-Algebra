# 第 72A 章 矩阵概率分布

<div class="context-flow" markdown>

**前置**：内积(Ch8) · 特征值(Ch6) · 指数矩阵(Ch13) · 随机矩阵(Ch23)

**本章脉络**：随机向量与协方差矩阵 → 多元正态分布(正定矩阵应用) → Wishart 分布(随机矩阵基础) → 矩阵正态分布 → Kronecker 积在概率中的作用 → 特征函数的矩阵形式

**延伸**：矩阵分布理论是金融组合分析、无线通信信道建模以及大规模高斯过程回归的数学地基

</div>

矩阵概率理论研究多维随机变量之间的二阶矩结构及其在矩阵空间上的分布特性。通过将概率密度函数表示为关于矩阵不变量（迹、行列式、二次型）的泛函，线性代数为描述高维随机关联提供了严谨的参数化手段。

---

## 72A.1 协方差结构与多元正态分布

!!! definition "定义 72A.1 (协方差算子)"
    设随机向量 $\mathbf{X} \in \mathbb{R}^p$。其协方差矩阵 $\Sigma = \mathbb{E}[(\mathbf{X}-\mu)(\mathbf{X}-\mu)^T]$ 是对称半正定阵。$\Sigma$ 编码了各维度间的线性相关性及方差分布。

!!! theorem "定理 72A.1 (Wishart 随机矩阵)"
    设 $\mathbf{X}_1, \dots, \mathbf{X}_n$ 是独立服从 $N_p(0, \Sigma)$ 的随机向量。则散度矩阵 $S = \sum \mathbf{X}_i \mathbf{X}_i^T$ 服从参数为 $(n, \Sigma)$ 的 Wishart 分布。它是样本协方差矩阵的代数基础。

---

## 练习题

1. **[半正定性] 证明：对任意随机向量，其协方差矩阵 $\Sigma$ 必然是半正定的。**
   ??? success "参考答案"
       对于任意常向量 $\mathbf{a} \in \mathbb{R}^p$，标量随机变量 $\mathbf{a}^T \mathbf{X}$ 的方差为 $\operatorname{Var}(\mathbf{a}^T \mathbf{X}) = \mathbf{a}^T \Sigma \mathbf{a}$。由于实随机变量的方差物理上非负，该二次型对所有 $\mathbf{a}$ 均非负，故 $\Sigma \succeq 0$。

2. **[线性变换] 设 $\mathbf{Y} = A \mathbf{X} + \mathbf{b}$，其中 $\mathbf{X} \sim (\mu, \Sigma)$。推导 $\mathbf{Y}$ 的协方差矩阵表达式。**
   ??? success "参考答案"
       $\Sigma_Y = \mathbb{E}[A(\mathbf{X}-\mu)(A(\mathbf{X}-\mu))^T] = A \mathbb{E}[(\mathbf{X}-\mu)(\mathbf{X}-\mu)^T] A^T = A \Sigma A^T$。这展示了协方差在仿射变换下的演化规则。

3. **[独立性判定] 证明：在多元正态分布中，各分量相互独立的充要条件是协方差矩阵 $\Sigma$ 为对角阵。**
   ??? success "参考答案"
       多元正态分布的密度函数包含项 $\exp(-\frac{1}{2}(\mathbf{x}-\mu)^T \Sigma^{-1} (\mathbf{x}-\mu))$。若 $\Sigma$ 为对角阵，则指数项可分解为各分量的平方和，密度函数退化为一元正态密度之积，满足独立性定义。

4. **[Mahalanobis距离] 计算矩阵 $A = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$ 的逆，并解释其在 Mahalanobis 距离 $\sqrt{\mathbf{x}^T \Sigma^{-1} \mathbf{x}}$ 中的测度作用。**
   ??? success "参考答案"
       $\Sigma^{-1} = \frac{4}{3} \begin{pmatrix} 1 & -0.5 \\ -0.5 & 1 \end{pmatrix}$。逆矩阵起到权重调节作用：在相关性强的方向上（对应 $\Sigma$ 的大特征值），逆矩阵会缩减该方向的欧氏距离，实现对数据尺度与相关性的归一化。

5. **[Wishart期望] 证明样本散度矩阵 $S$ 的期望为 $\mathbb{E}[S] = n \Sigma$。**
   ??? success "参考答案"
       $\mathbb{E}[\sum \mathbf{X}_i \mathbf{X}_i^T] = \sum \mathbb{E}[\mathbf{X}_i \mathbf{X}_i^T]$。由于 $\mathbb{E}[\mathbf{X}_i]=0$，有 $\mathbb{E}[\mathbf{X}_i \mathbf{X}_i^T] = \Sigma$。累加 $n$ 次得 $n \Sigma$。

6. **[精度矩阵] 定义精度矩阵 $\Omega = \Sigma^{-1}$，并说明其分量 $\omega_{ij}$ 与偏相关系数的关系。**
   ??? success "参考答案"
       $\Omega$ 反映了条件相关性。在高斯图形模型中，$\omega_{ij} = 0$ 意味着在给定其余变量的情况下，变量 $i$ 与 $j$ 条件独立。其分量标准化后即为偏相关系数。

7. **[特征结构] 分析协方差矩阵的特征值分解与主成分分析（PCA）中方差最大化方向的等价性。**
   ??? success "参考答案"
       最大特征值对应的特征向量是 $\max_{\mathbf{v}^T \mathbf{v}=1} \mathbf{v}^T \Sigma \mathbf{v}$ 的解。这在几何上对应了随机点云散布最广的主轴方向。

8. **[信息熵] 证明：多元正态分布的微分熵 $H$ 与 $\det(\Sigma)$ 的对数成线性关系。**
   ??? success "参考答案"
       $H(\mathbf{X}) = \frac{p}{2}(1 + \ln(2\pi)) + \frac{1}{2} \ln \det(\Sigma)$。行列式量化了随机变量在空间中占据的体积，从而反映了系统的信息不确定度。

9. **[矩阵正态性] 解释矩阵正态分布 $\mathcal{MN}_{n \times p}(M, U, V)$ 中 Kronecker 积 $V \otimes U$ 对协方差结构的参数压缩意义。**
   ??? success "参考答案"
       它假设行间相关性 $U$ 与列间相关性 $V$ 是解耦的。这使得具有 $n^2 p^2$ 个元素的完整协方差矩阵被简化为仅需 $(n^2+p^2)$ 个参数，极大降低了估计的维度灾难。

10. **[特征函数] 写出多元正态分布的特征函数 $\phi(\mathbf{t}) = \mathbb{E}[e^{j \mathbf{t}^T \mathbf{X}}]$ 的矩阵形式，并分析其指数项中的二次型。**
    ??? success "参考答案"
        $\phi(\mathbf{t}) = \exp(j \mathbf{t}^T \mu - \frac{1}{2} \mathbf{t}^T \Sigma \mathbf{t})$。指数中的二次型 $\mathbf{t}^T \Sigma \mathbf{t}$ 确定了频率域内概率密度的收敛特性，是协方差结构在 Fourier 变换下的映射。

## 本章小结

本章论述了随机分析中的矩阵参数化建模：

1. **二阶矩表达**：确立了对称半正定矩阵作为描述随机向量空间关联的通用工具。
2. **生成分布**：建立了 Wishart 分布作为样本统计量分析的矩阵代数标准。
3. **结构化协方差**：利用 Kronecker 积与逆矩阵（精度矩阵）揭示了高维数据间的条件依赖与独立性规律。
4. **不变量度量**：通过行列式与迹量化了随机过程的熵增与能量分布特性。
