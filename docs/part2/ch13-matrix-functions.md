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

1. **[计算] 求 $e^A$，其中 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。**
   ??? success "参考答案"
       $A^2 = O$，故 $e^A = I + A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$。

2. **[性质] 证明 $f(P A P^{-1}) = P f(A) P^{-1}$。**
   ??? success "参考答案"
       利用幂级数定义：$(P A P^{-1})^k = P A^k P^{-1}$。代入级数求和式，提取 $P$ 和 $P^{-1}$ 即得。

3. **[对角阵] 若 $A = \operatorname{diag}(\lambda_1, \lambda_2)$，求 $\sin A$。**
   ??? success "参考答案"
       $\sin A = \operatorname{diag}(\sin \lambda_1, \sin \lambda_2)$。

4. **[交换性] 证明 $A f(A) = f(A) A$。**
   ??? success "参考答案"
       由于 $A$ 与自身的任何幂 $A^k$ 都交换，且级数收敛，故与级数和也交换。

5. **[行列式] 若 $\operatorname{tr}(A) = 0$，求 $\det(e^A)$。**
   ??? success "参考答案"
       $\det(e^A) = e^{\operatorname{tr}(A)} = e^0 = 1$。

6. **[对数] 若 $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$，求 $\log A$。**
   ??? success "参考答案"
       令 $A = I + N$，其中 $N = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。由于 $N^2=0$，$\log(I+N) = N - N^2/2 + \cdots = N = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。

7. **[插值] 使用插值法求 $f(A)$，已知 $A$ 的唯一特征值为 2，代数重数为 2。**
   ??? success "参考答案"
       寻找 $q(\lambda) = a\lambda + b$ 满足 $q(2)=f(2)$ 且 $q'(2)=f'(2)$。解出 $a, b$ 后，$f(A) = aA + bI$。

8. **[ODE] 方程 $\mathbf{x}' = A\mathbf{x}, \mathbf{x}(0) = \mathbf{x}_0$ 的解是什么？**
   ??? success "参考答案"
       $\mathbf{x}(t) = e^{At} \mathbf{x}_0$。

9. **[三角函数] 证明 $\cos^2 A + \sin^2 A = I$。**
   ??? success "参考答案"
       利用 $\cos A = \frac{e^{iA}+e^{-iA}}{2}$ 和 $\sin A = \frac{e^{iA}-e^{-iA}}{2i}$ 展开验证即可。

10. **[应用] 为什么在处理非对角化矩阵时，函数值中会出现导数项？**
    ??? success "参考答案"
        这反映了矩阵作用的“耦合”效应。在 Jordan 块结构中，偏离对角线的 1 导致了函数展开时高阶无穷小的积累，数学上表现为 Taylor 展开的导数项。

## 本章小结

矩阵函数是算子理论的解析延拓：

1.  **一致性原则**：通过级数、Jordan 形或插值定义的矩阵函数在收敛域内是等价的，保证了代数与分析的逻辑统一。
2.  **指数核心**：矩阵指数 $e^A$ 是最核心的函数，它将线性的微分演化转化为了纯粹的矩阵乘法，是处理时间演化问题的终极钥匙。
3.  **计算多样性**：插值法和 Jordan 分解提供了两种互补的计算视角——前者侧重于特征值的局部解析，后者侧重于全空间的结构解构。
