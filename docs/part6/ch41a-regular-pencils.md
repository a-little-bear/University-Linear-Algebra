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

本章将特征值理论扩展到了矩阵对：


****：确立了良置广义特征值问题的判定准则。

****：开发了 Weierstrass 标准形，用于严格等价下的束分类。

****：引入了 QZ 分解及其算法，作为谱计算的数值稳健路径。

****：通过线性化技术桥接了高阶矩阵多项式与线性矩阵束。
