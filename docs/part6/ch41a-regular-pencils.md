# 第 41A 章 矩阵束与正则束理论

<div class="context-flow" markdown>

**前置**：特征值与特征向量(Ch6) · Jordan 标准形(Ch12) · $\lambda$-矩阵与 Smith 标准形(Ch13B) · Schur 分解(Ch10)

**本章脉络**：矩阵束定义 $\to$ 正则束与奇异束分类 $\to$ 广义特征值（齐次坐标） $\to$ 严格等价 $\to$ Weierstrass 标准形 $\to$ 收缩子空间 $\to$ 广义 Schur 分解（QZ 分解） $\to$ QZ 算法 $\to$ Hermite 束 $\to$ 多项式特征值问题与线性化

**延伸**：正则束理论是 LAPACK 中广义特征值求解器的数学基础；QZ 算法是 MATLAB `eig(A, B)` 的核心

</div>

标准特征值问题 $Ax = \lambda x$ 研究的是矩阵束 $A - \lambda I$。当 $I$ 替换为一般矩阵 $B$ 时，便得到**广义特征值问题** $Ax = \lambda Bx$。若 $\det(A - \lambda B)$ 不恒为零，则称束为**正则的**；此时存在 Weierstrass 标准形——一种将有限特征值与无穷特征值分离的代数结构。

---

## 41A.1 核心概念

!!! definition "定义 41A.2 (正则与奇异)"
    方阵束 $A - \lambda B$ 称为**正则的**，若 $\det(A - \lambda B) \not\equiv 0$。否则称为**奇异的**。

!!! theorem "定理 41A.3 (Weierstrass 标准形)"
    对于正则束 $A - \lambda B$，存在非奇异矩阵 $P, Q$ 使得
    $$P(A - \lambda B)Q = \begin{pmatrix} J - \lambda I & 0 \\ 0 & I - \lambda N \end{pmatrix},$$
    其中 $J$ 为 Jordan 形（对应有限特征值），$N$ 为幂零阵（对应无穷特征值）。

---

## 练习题

1. **[正则性] 判定 $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, B = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$ 是否正则，并求其特征值。**
   ??? success "参考答案"
       $\det(A - \lambda B) = \det \begin{pmatrix} 1 & 0 \\ 0 & -\lambda \end{pmatrix} = -\lambda$。不恒为零，故正则。有限特征值 $\lambda=0$；包含 1 个无穷特征值。

2. **[齐次坐标] 为什么广义特征值要使用齐次坐标 $(\alpha, \beta)$？**
   ??? success "参考答案"
       齐次坐标统一处理了有限特征值（$\lambda = \alpha/\beta$）和无穷特征值（$\beta=0$），将整个谱映射到复射影直线 $\mathbb{P}^1(\mathbb{C})$ 上，避免了无穷大的奇异性。

3. **[严格等价] 定义矩阵束的严格等价。**
   ??? success "参考答案"
       若存在非奇异阵 $P, Q$ 使得 $PA_1Q = A_2$ 且 $PB_1Q = B_2$，则称两矩阵束严格等价。

4. **[QZ分解] 什么是广义 Schur 分解（QZ 分解）？**
   ??? success "参考答案"
       对任意矩阵对 $(A, B)$，存在酉矩阵 $Q, Z$ 使得 $Q^*AZ$ 和 $Q^*BZ$ 均为上三角矩阵。对角元之比即给出广义特征值。

5. **[收缩子空间] 区分不变子空间与收缩子空间（Deflating Subspace）。**
   ??? success "参考答案"
       不变子空间满足 $A\mathcal{V} \subseteq \mathcal{V}$。矩阵对 $(A, B)$ 的收缩子空间满足 $\dim(A\mathcal{V} + B\mathcal{V}) \le \dim \mathcal{V}$。

6. **[Hermite束] 何时 Hermite 束 $A - \lambda B$ 的广义特征值全是实数？**
   ??? success "参考答案"
       当矩阵束是**确定的**（Definite）时，即存在线性组合 $\alpha A + \beta B$ 是正定矩阵时，所有特征值均为实数。

7. **[条件数] 定义简单广义特征值 $\lambda_0$ 的条件数。**
   ??? success "参考答案"
       $\kappa(\lambda_0) = \frac{\|y\| \|x\|}{|y^* B x|}$，其中 $x, y$ 分别为左右特征向量。当 $B$ 趋于奇异或矩阵束趋于奇异时，特征值变得高度敏感（病态）。

8. **[线性化] 如何求解二次特征值问题 $(\lambda^2 M + \lambda C + K)x = 0$？**
   ??? success "参考答案"
       将其转化为两倍维数的线性矩阵束（线性化），通常采用伴随形式：$\lambda \begin{pmatrix} M & 0 \\ 0 & I \end{pmatrix} - \begin{pmatrix} -C & -K \\ I & 0 \end{pmatrix}$。

9. **[保结构] 什么是保结构线性化？**
   ??? success "参考答案"
       指在将多项式特征值问题线性化时，保持原矩阵序列的对称性或 Hermite 性等物理结构不被破坏的线性化方法。

10. **[Weierstrass回归] 当 $B = I$ 时，Weierstrass 标准形如何退化为 Jordan 标准形？**
    ??? success "参考答案"
        若 $B=I$，则不存在无穷特征值部分（$s=0$），标准形简化为 $P(A - \lambda I)Q = J - \lambda I$，这正是标准矩阵相似下的 Jordan 标准形问题。

## 本章小结

本章将特征值理论扩展到了矩阵对：

1. **正则性框架**：确立了良置广义特征值问题的判定准则。
2. **规范分解**：开发了 Weierstrass 标准形，用于严格等价下的束分类。
3. **酉算法**：引入了 QZ 分解及其算法，作为谱计算的数值稳健路径。
4. **多项式映射**：通过线性化技术桥接了高阶矩阵多项式与线性矩阵束。
