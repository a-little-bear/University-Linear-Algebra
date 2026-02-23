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

1. **[Vandermonde] 计算 $V(1, 2, 3)$ 的行列式。**
   ??? success "参考答案"
       $\det V = (2-1)(3-1)(3-2) = 1 \cdot 2 \cdot 1 = 2$。

2. **[条件数] 为什么在 $[0, 1]$ 上取等距节点的 Vandermonde 矩阵是病态的？**
   ??? success "参考答案"
       其条件数随 $n$ 指数级增长。等距节点上的高次多项式插值会引发 Runge 现象，导致插值系数对数据扰动极度敏感。

3. **[Cauchy] 证明 Cauchy 矩阵 $C = (1/(x_i - y_j))$ 关于对角算子的位移秩为 1。**
   ??? success "参考答案"
       设 $D_x = \operatorname{diag}(x_i), D_y = \operatorname{diag}(y_j)$。则 $(D_x C - C D_y)_{ij} = x_i \frac{1}{x_i - y_j} - \frac{1}{x_i - y_j} y_j = \frac{x_i - y_j}{x_i - y_j} = 1$。结果是全 1 矩阵，秩为 1。

4. **[Toeplitz] 证明 Toeplitz 矩阵的逆不一定是 Toeplitz 的，但它是“类 Toeplitz”的。**
   ??? success "参考答案"
       虽然条目的对称性消失，但位移秩得以保持。若 $T$ 的位移秩为 2，则 $T^{-1}$ 相对位移算子的位移秩也为 2，这允许通过 FFT 实现 $O(n \log n)$ 的矩阵-向量乘法。

5. **[Krylov] 当 $A$ 是对角阵时，Krylov 矩阵 $K(A, b)$ 与 Vandermonde 矩阵有什么关系？**
   ??? success "参考答案"
       若 $A = \operatorname{diag}(x_i)$，则 $(A^k b)_i = x_i^k b_i$。因此 $K(A, b) = \operatorname{diag}(b) V(x)$，即行缩放的 Vandermonde 矩阵。

6. **[复杂度] 对比求解一般 $n$ 阶系统与 Toeplitz 系统的复杂度。**
   ??? success "参考答案"
       一般系统：$O(n^3)$（Gauss 消元）。Toeplitz 系统：$O(n^2)$（Levinson-Durbin）或 $O(n \log^2 n)$（超快求解器）。

7. **[Gohberg-Semencul] Gohberg-Semencul 公式的意义是什么？**
   ??? success "参考答案"
       它将 Toeplitz 矩阵的逆表示为三角 Toeplitz 矩阵乘积之差。这允许我们在不显式写出逆矩阵的情况下，利用 FFT 在 $O(n \log n)$ 内计算 $T^{-1}v$。

8. **[结式] Sylvester 结式矩阵 $\operatorname{Syl}(f, g)$ 何时是奇异的？**
   ??? success "参考答案"
       当且仅当多项式 $f(x)$ 和 $g(x)$ 有公共根时，即 $\gcd(f, g) \neq 1$。

9. **[H-矩阵] 简述层次矩阵 ($\mathcal{H}$-matrix) 的核心思想。**
   ??? success "参考答案"
       虽然整个矩阵可能满秩，但其代表“远场”相互作用的非对角块具有低秩结构。通过递归分块，可以实现 $O(n \log n)$ 的算术运算。

10. **[Schur算法] 广义 Schur 算法如何利用位移生成元？**
    ??? success "参考答案"
        它直接在低秩生成元 $(G, H)$ 上进行消元，而不是在 $n^2$ 个条目上操作，从而将计算量减少到 $O(rn^2)$，其中 $r$ 是位移秩。

## 本章小结

位移结构是快速线性代数的统一理论：

1. **低秩映射**：结构化矩阵的特征是在特定算子作用下具有低秩的“映射”。
2. **结构继承**：这种结构在求逆和求 Schur 补运算下是保持的。
3. **算法飞跃**：通过操作生成元，我们为广泛的工程问题打破了 $O(n^3)$ 的计算壁垒。
