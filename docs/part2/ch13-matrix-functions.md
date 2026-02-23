# 第 13 章 矩阵函数

<div class="context-flow" markdown>

**前置**：Jordan 标准形 (Ch12) · 矩阵分析基础 (Ch14) · 微积分幂级数

**本章脉络**：从标量函数到矩阵函数 $\to$ 矩阵函数的泰勒级数定义 $\to$ 通过 Jordan 标准形的一般定义 $\to$ Sylvester-Lagrange 插值计算法 $\to$ 核心函数：矩阵指数 ($e^A$) 与对数 ($\ln A$) $\to$ 矩阵指数的性质与运算律 $\to$ 矩阵三角函数 $\to$ 应用：线性常微分方程组 (ODE) 的解析解、量子演化算符、控制理论的状态转移矩阵

**延伸**：矩阵函数将算术代数提升到了解析代数的高度；它是连接离散矩阵结构与连续动力学演化的终极纽带，是描述一切随时间演化的线性系统的核心数学语言

</div>

矩阵函数（Matrix Functions）研究的是如何将经典的标量函数（如 $e^x, \sin x, \ln x$）作用于矩阵变量。这不仅仅是简单的逐元素运算，而是保持矩阵代数结构的算子运算。矩阵函数是处理时间演化、信号传播以及复杂系统响应的核心工具。本章将确立定义矩阵函数的三种等价路径，并深入探讨最重要的算子——矩阵指数。

---

## 13.1 定义方法：从级数到 Jordan 形

!!! definition "定义 13.1 (通过泰勒级数定义)"
    若标量函数 $f(z)$ 的泰勒级数为 $\sum_{k=0}^\infty c_k z^k$，且其收敛半径大于矩阵 $A$ 的所有特征值的模，则定义：
    $$f(A) = \sum_{k=0}^\infty c_k A^k$$

!!! definition "定义 13.2 (通过 Jordan 标准形定义)"
    若 $A = P J P^{-1}$，则 $f(A) = P f(J) P^{-1}$。
    对于每个 Jordan 块 $J_k(\lambda)$，其函数值定义为由 $f$ 及其导数构成的上三角阵：
    $$f(J_k(\lambda)) = \begin{pmatrix} f(\lambda) & f'(\lambda) & \frac{f''(\lambda)}{2!} & \cdots & \frac{f^{(k-1)}(\lambda)}{(k-1)!} \\ 0 & f(\lambda) & f'(\lambda) & \cdots & \vdots \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & \cdots & f(\lambda) & f'(\lambda) \\ 0 & 0 & \cdots & 0 & f(\lambda) \end{pmatrix}$$

---

## 13.2 矩阵指数 $e^A$

!!! theorem "定理 13.1 (矩阵指数的性质)"
    1.  **收敛性**：对任何方阵 $A$，级数 $\sum \frac{A^k}{k!}$ 始终绝对收敛。
    2.  **微分性质**：$\frac{d}{dt} e^{At} = A e^{At}$。
    3.  **乘法性质**：若 $AB = BA$，则 $e^{A+B} = e^A e^B$。
    4.  **行列式关系**：$\det(e^A) = e^{\operatorname{tr}(A)}$（Jacobi 恒等式）。

---

## 13.3 计算技术：插值法

!!! technique "技术：Sylvester-Lagrange 插值"
    若 $A$ 的特征值为 $\lambda_1, \ldots, \lambda_m$，其在最小多项式中的重数为 $n_1, \ldots, n_m$。
    存在一个次数小于 $\sum n_i$ 的多项式 $p(z)$，使得 $p(z)$ 及其导数在各特征值处与 $f(z)$ 吻合。
    则：$f(A) = p(A)$。
    **优点**：不需要求解完整的 Jordan 分解，仅需特征值和矩阵幂。

---

## 练习题

**1. [计算] 求 $e^A$，其中 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。**

??? success "参考答案"
    **步骤：**
    1. 观察矩阵性质：这是一个幂零阵，$A^2 = O$。
    2. 利用幂级数定义：$e^A = I + A + \frac{A^2}{2!} + \cdots$。
    3. 代入：$e^A = I + A + O = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。
    **结论**：$e^A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$。

**2. [性质] 证明若 $A$ 为对角阵 $\operatorname{diag}(d_1, \ldots, d_n)$，则 $f(A) = \operatorname{diag}(f(d_1), \ldots, f(d_n))$。**

??? success "参考答案"
    **证明：**
    1. 利用级数定义：$A^k = \operatorname{diag}(d_1^k, \ldots, d_n^k)$。
    2. $f(A) = \sum c_k A^k = \sum c_k \operatorname{diag}(d_1^k, \ldots, d_n^k)$。
    3. 按分量合并：$= \operatorname{diag}(\sum c_k d_1^k, \ldots, \sum c_k d_n^k)$。
    4. 每一项正好是 $f(d_i)$ 的级数定义。
    **结论**：对角阵的函数运算等价于对各分量分别求函数值。

**3. [交换性] 证明：对于任何 $f(A)$，必有 $A f(A) = f(A) A$。**

??? success "参考答案"
    **证明：**
    1. 使用多项式/级数表示：$f(A) = \sum c_k A^k$。
    2. $A f(A) = A (\sum c_k A^k) = \sum c_k A^{k+1}$。
    3. $f(A) A = (\sum c_k A^k) A = \sum c_k A^{k+1}$。
    **结论**：矩阵总是与其自身的函数值交换。这反映了算子与其自身生成的代数具有交换性。

**4. [行列式] 若 $\operatorname{tr}(A) = 0$，求 $\det(e^A)$。**

??? success "参考答案"
    **定理应用：**
    利用公式 $\det(e^A) = e^{\operatorname{tr}(A)}$。
    **计算：**
    $\det(e^A) = e^0 = 1$。
    **物理应用**：这说明如果系统演化的生成元是迹为零的（如旋转），则演化过程保持体积不变（行列式为 1）。

**5. [对数] 若 $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$，求 $\ln A$。**

??? success "参考答案"
    **方法：**
    1. 令 $A = I + N$，其中 $N = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。
    2. 利用级数：$\ln(I+N) = N - \frac{N^2}{2} + \frac{N^3}{3} - \cdots$。
    3. 由于 $N^2 = O$，级数截断。
    **结论**：$\ln A = N = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。

**6. [插值法] 使用插值法求 $f(A)$，已知 $A$ 的唯一特征值为 2，代数重数为 2。**

??? success "参考答案"
    **步骤：**
    1. 最小多项式为 $(z-2)^2$。
    2. 设插值多项式为 $p(z) = a z + b$。
    3. 匹配条件：$p(2) = f(2)$ 且 $p'(2) = f'(2)$。
    4. 解方程：$2a + b = f(2)$，$a = f'(2)$。
    5. 得 $b = f(2) - 2f'(2)$。
    **结论**：$f(A) = f'(2) A + [f(2) - 2f'(2)] I$。

**7. [ODE] 若系统满足 $\mathbf{x}' = A\mathbf{x}$，且 $\mathbf{x}(0) = \mathbf{x}_0$，其解析解是什么？**

??? success "参考答案"
    **结论：**
    $\mathbf{x}(t) = e^{At} \mathbf{x}_0$。
    **理由**：对该式求导 $\frac{d}{dt}(e^{At} \mathbf{x}_0) = A e^{At} \mathbf{x}_0 = A \mathbf{x}(t)$，且在 $t=0$ 时，$e^O = I$，满足初始条件。

**8. [幂运算] 证明 $(e^A)^k = e^{kA}$ 对任何整数 $k$ 成立。**

??? success "参考答案"
    **证明：**
    由于 $A$ 与 $A$ 总是交换的，根据乘法律：$e^A e^A = e^{A+A} = e^{2A}$。利用数学归纳法推广到 $k$。对于负数，$e^A e^{-A} = e^O = I$，故 $e^{-A} = (e^A)^{-1}$。

**9. [三角函数] 若 $A^2 = -I$，求 $\cos A$。**

??? success "参考答案"
    **级数展开：**
    $\cos A = I - \frac{A^2}{2!} + \frac{A^4}{4!} - \cdots$。
    代入 $A^2 = -I, A^4 = I \ldots$：
    $\cos A = I(1 - (-1)/2! + 1/4! - \cdots)$。
    这对应标量 $\cos(i)$ 的展开？不，注意符号：$= I(1 + 1/2! + 1/4! + \cdots) = I \cosh(1)$。
    **结论**：$\cos A = (\cosh 1) I$。

**10. [应用] 为什么在数值计算中，大矩阵的 $e^A$ 很难计算？**

??? success "参考答案"
    **难点分析：**
    1. **舍入误差**：泰勒级数在项数较多时，正负项相互抵消会产生巨大的相对误差。
    2. **缩放平方开销**：常用算法是先缩小 $A$（使 $\|A\|$ 小），计算 $e^{A/2^k}$，再通过 $2^k$ 次平方。这虽然稳定，但大矩阵的连续乘法计算开销巨大。
    3. **条件数**：如果 $A$ 接近亏损（非正规），计算过程极度敏感。

## 本章小结

矩阵函数是算子理论的解析延拓：

1.  **一致性原则**：无论是通过级数、Jordan 形还是插值定义，只要函数在谱上解析，定义就是等价的，保证了逻辑的严密性。
2.  **指数核心**：矩阵指数 $e^A$ 是线性动力系统的 DNA，它将复杂的微分演化压缩为纯粹的代数算子作用。
3.  **计算的多样性**：插值法提供了一种绕过 Jordan 分解计算函数值的捷径，揭示了矩阵性质完全由其在特征点处的局部解析行为决定。
