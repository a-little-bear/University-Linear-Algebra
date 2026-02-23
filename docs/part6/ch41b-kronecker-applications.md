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

本章为所有线性算子对提供了最终的分类方案：

1. **奇异微积分**：定义了最小指标作为多项式零空间的结构不变量。
2. **全局分解**：在 Kronecker 标准形中统一了正则与奇异行为。
3. **控制链接**：建立了矩阵结构与线性系统反馈属性之间的等价关系。
4. **动力学约束**：将理论应用于 DAE，确立了指标作为数值和分析复杂度的度量。
