# 第 13 章 矩阵函数

<div class="context-flow" markdown>

**前置**：Jordan 标准形 (Ch12) · 矩阵分析基础 (Ch14) · 微积分级数理论

**本章脉络**：从标量函数到矩阵函数 $\to$ 幂级数定义 $\to$ Jordan 标准形定义法 $\to$ Sylvester-Lagrange 插值定义法 $\to$ 核心函数：矩阵指数 ($e^A$)、矩阵对数 ($\log A$)、矩阵三角函数 $\to$ 计算技术 $\to$ 性质与恒等式 $\to$ 在常微分方程组 (ODE) 中的应用

**延伸**：矩阵函数将算术代数提升到了分析代数的高度；它是控制理论、量子力学演化算符以及连续动力系统的核心数学语言

</div>

矩阵函数（Matrix Functions）是研究如何将标量函数（如 $e^x, \sin x, \log x$）作用于矩阵变量的学科。这不仅仅是简单的逐元素运算，而是保持矩阵代数结构的算子运算。矩阵函数是连接离散矩阵代数与连续物理世界的桥梁。

---

## 13.1 定义方法

!!! definition "定义 13.1 (通过 Jordan 标准形定义)"
    若 $A = P J P^{-1}$，其中 $J = \operatorname{diag}(J_1, \ldots, J_m)$。
    对每个 Jordan 块 $J_k(\lambda)$：
    $$f(J_k(\lambda)) = \begin{pmatrix} f(\lambda) & f'(\lambda) & \frac{f''(\lambda)}{2!} & \cdots & \frac{f^{(k-1)}(\lambda)}{(k-1)!} \\ 0 & f(\lambda) & f'(\lambda) & \cdots & \vdots \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & \cdots & f(\lambda) & f'(\lambda) \\ 0 & 0 & \cdots & 0 & f(\lambda) \end{pmatrix}$$
    则 $f(A) = P \operatorname{diag}(f(J_1), \ldots, f(J_m)) P^{-1}$。

!!! definition "定义 13.2 (通过幂级数定义)"
    若标量函数 $f(z)$ 的泰勒级数为 $\sum_{k=0}^\infty c_k z^k$，且其收敛半径大于 $A$ 的所有特征值的模，则定义：
    $$f(A) = \sum_{k=0}^\infty c_k A^k$$

---

## 13.2 矩阵指数 $e^A$

!!! theorem "定理 13.1 (矩阵指数的性质)"
    1.  **定义**：$e^A = I + A + \frac{A^2}{2!} + \cdots$。该级数对所有方阵 $A$ 绝对收敛。
    2.  **微分**：$\frac{d}{dt} e^{At} = A e^{At}$。这是求解 $\mathbf{x}' = A\mathbf{x}$ 的关键。
    3.  **乘法**：若 $AB = BA$，则 $e^{A+B} = e^A e^B$。
    4.  **行列式**：$\det(e^A) = e^{\operatorname{tr}(A)}$。

---

## 13.3 计算方法：插值法

!!! technique "技术：Sylvester-Lagrange 插值"
    若 $A$ 的特征值为 $\lambda_1, \ldots, \lambda_k$，对应的最大 Jordan 块阶数为 $n_1, \ldots, n_k$。
    寻找一个多项式 $q(\lambda)$ 使得它及其导数在 $\lambda_i$ 处的值与 $f(\lambda)$ 及其导数相等。则：
    $$f(A) = q(A)$$
    这避免了复杂的 Jordan 分解，仅需特征值和矩阵幂。

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

矩阵函数是算子理论的解析延拓：


****：通过级数、Jordan 形或插值定义的矩阵函数在收敛域内是等价的，保证了代数与分析的逻辑统一。

****：矩阵指数 $e^A$ 是最核心的函数，它将线性的微分演化转化为了纯粹的矩阵乘法，是处理时间演化问题的终极钥匙。

****：插值法和 Jordan 分解提供了两种互补的计算视角——前者侧重于特征值的局部解析，后者侧重于全空间的结构解构。
