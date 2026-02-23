# 第 37B 章 Vandermonde、Cauchy 矩阵与位移结构

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 行列式(Ch3) · 矩阵分解(Ch10) · 数值线性代数(Ch22) · Toeplitz/Hankel/循环矩阵(Ch37A)

**本章脉络**：Vandermonde 矩阵与多项式插值 $\to$ Cauchy 矩阵与位移秩 1 结构 $\to$ Resultant 矩阵与 Krylov 矩阵 $\to$ Sylvester/Stein 位移算子 $\to$ 各类矩阵的位移秩计算 $\to$ 广义 Schur 算法与 Gohberg-Semencul 公式 $\to$ 层次矩阵与现代推广

**延伸**：位移结构理论（Kailath-Kung-Morf, 1979）统一了所有经典结构化矩阵类的快速算法框架

</div>

在第 37A 章中，我们研究了以“沿对角线常数”为特征的结构化矩阵。本章转向 Vandermonde 和 Cauchy 矩阵，它们分别以“幂次结构”和“核函数结构”为特征。我们引入**位移结构**的统一框架，将这些矩阵纳入同一理论体系，揭示它们共同的低秩位移性质，从而推导出 $O(n^2)$ 甚至 $O(n \log^2 n)$ 的快速算法。

---

## 37B.1 核心概念

!!! definition "定义 37B.1 (Vandermonde 矩阵)"
    给定节点 $x_1, \ldots, x_n$，Vandermonde 矩阵 $V$ 定义为 $V_{ij} = x_i^{j-1}$。当 $x_i$ 两两不同时，它是非奇异的。

!!! definition "定义 37B.8 (Sylvester 位移)"
    矩阵 $A$ 关于算子 $(F, F')$ 的位移定义为 $\nabla_{F,F'}(A) = FA - AF'$。该结果矩阵的秩被称为 **位移秩**。

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

位移结构是快速线性代数的统一理论：

1. **低秩映射**：结构化矩阵的特征是在特定算子作用下具有低秩的“映射”。
2. **结构继承**：这种结构在求逆和求 Schur 补运算下是保持的。
3. **算法飞跃**：通过操作生成元，我们为广泛的工程问题打破了 $O(n^3)$ 的计算壁垒。
