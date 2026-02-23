# 第 41B 章 Kronecker 标准形与应用

<div class="context-flow" markdown>

**前置**：矩阵束与正则束理论(Ch41A) · Jordan 标准形(Ch12) · $\lambda$-矩阵与 Smith 标准形(Ch13B) · 奇异值分解(Ch11) · 控制论基础

**本章脉络**：最小指标 $\to$ 最小基理论（Forney） $\to$ Kronecker 块 $\to$ Kronecker 标准形 $\to$ 阶梯算法（Van Dooren） $\to$ Brunovsky 标准形 $\to$ 微分代数方程（DAE） $\to$ 描述子系统 $\to$ 结构化矩阵束

**延伸**：Kronecker 标准形是奇异矩阵束的最终分类结果；DAE 指标理论是电路模拟和多体动力学的数学基础

</div>

当矩阵束 $A - \lambda B$ 是奇异的——即行列式恒为零或为长方形矩阵时，Weierstrass 形不再适用。Kronecker 的理论给出了任意矩阵束在严格等价下的完全分类。其结构由**最小指标**捕捉，描述了多项式零空间的“阶梯状”结构。

---

## 41B.1 核心概念

!!! definition "定义 41B.1 (最小指标)"
    右（列）最小指标 $\varepsilon_1, \dots, \varepsilon_p$ 是多项式右零空间 $\mathcal{N}_r(\lambda) = \{x(\lambda) : (A - \lambda B)x(\lambda) = 0\}$ 的一组**最小基**中向量的次数。

!!! theorem "定理 41B.3 (Kronecker 标准形)"
    任何 $m \times n$ 矩阵束 $A - \lambda B$ 均严格等价于一个分块对角矩阵，由以下三类块组成：
    - 右 Kronecker 块 $L_\varepsilon$
    - 左 Kronecker 块 $L_\eta^T$
    - 正则部分（有限和无穷 Jordan 块）

---

## 练习题

1. **[最小指标] 求 $L(\lambda) = \begin{pmatrix} 1 & -\lambda & 0 \\ 0 & 1 & -\lambda \end{pmatrix}$ 的右最小指标。**

   ??? success "参考答案"
       求解 $L(\lambda)x(\lambda) = 0$ 得：$x_1 = \lambda x_2 = \lambda^2 x_3$。取 $x_3=1$ 得 $x(\lambda) = (\lambda^2, \lambda, 1)^T$。其最高次数为 2。故 $\varepsilon_1 = 2$。

2. **[维数] 一个秩为 2 的 $3 \times 5$ 矩阵束有多少个右最小指标和左最小指标？**

   ??? success "参考答案"
       右最小指标个数 $p = n - r = 5 - 2 = 3$。左最小指标个数 $q = m - r = 3 - 2 = 1$。满足 $p - q = n - m = 2$。

3. **[Kronecker块] 写出 $2 \times 3$ 的右 Kronecker 块 $L_2$ 的显式形式。**

   ??? success "参考答案"
       $L_2 = \begin{pmatrix} 1 & -\lambda & 0 \\ 0 & 1 & -\lambda \end{pmatrix}$。

4. **[Forney] 什么是最小基的“可预测次数性质”？**

   ??? success "参考答案"
       指基向量的任何多项式线性组合的次数，都精确等于各分量次数的最大值。这保证了最高次项不会发生意外的抵消。

5. **[阶梯算法] 为什么 Van Dooren 阶梯算法在数值上是稳定的？**

   ??? success "参考答案"
       因为它完全依赖于酉变换（SVD 或带选主元的 QR）来识别秩并提取结构指标，避免了多项式运算中的不稳定性。

6. **[Brunovsky] 将控制理论与 Kronecker 指标联系起来。**

   ??? success "参考答案"
       对于可控系统 $(A, B)$，其可控性指标 $\kappa_i$ 与矩阵束 $(\lambda I - A, -B)$ 的右最小指标 $\varepsilon_i$ 满足 $\varepsilon_i = \kappa_i - 1$。

7. **[DAE指标] 定义 DAE 方程 $E\dot{x} = Ax + f$ 的指标。**

   ??? success "参考答案"
       定义为矩阵束 $sE - A$ 的 Weierstrass 标准形中幂零块 $N$ 的幂零阶数（最大 Jordan 块的大小）。它衡量了求解 DAE 所需的求导次数。

8. **[相容初值] 为什么 DAE 限制了初值 $x(0)$ 的选择？**

   ??? success "参考答案"
       因为 DAE 的代数部分 $N\dot{z} = z + h(t)$ 意味着变量 $z(t)$ 完全由 $h$ 及其导数决定。因此 $z(0)$ 必须满足代数约束，不能任意选取。

9. **[描述子系统] 描述子系统 $(E, A, B, C, D)$ 有唯一解的条件是什么？**

   ??? success "参考答案"
       当且仅当矩阵束 $sE - A$ 是正则的（即 $\det(sE - A) \not\equiv 0$）。

10. **[谱对称性] $T$-回文束 $A - \lambda A^T$ 的谱具有什么对称性？**

   ??? success "参考答案"
        若 $\lambda_0$ 是特征值，则 $1/\lambda_0$ 也是特征值。

## 本章小结

本章为所有线性算子对提供了最终的分类方案：

1. **奇异微积分**：定义了最小指标作为多项式零空间的结构不变量。
2. **全局分解**：在 Kronecker 标准形中统一了正则与奇异行为。
3. **控制链接**：建立了矩阵结构与线性系统反馈属性之间的等价关系。
4. **动力学约束**：将理论应用于 DAE，确立了指标作为数值和分析复杂度的度量。
