# 第 47A 章 矩阵微积分基础

<div class="context-flow" markdown>

**前置**：矩阵运算 (Ch02) · 矩阵分析 (Ch14) · 多元微积分基础

**本章脉络**：从标量导数到矩阵导数 $\to$ 布局约定（Layout Conventions）：分子布局 vs. 分母布局 $\to$ 标量对向量的梯度（梯度向量） $\to$ 标量对矩阵的导数 $\to$ 核心公式：线性型 $a^T x$ 与二次型 $x^T Ax$ 的梯度 $\to$ 迹函数的导数（迹技巧） $\to$ 链式法则的矩阵形式 $\to$ 应用：普通最小二乘法 (OLS) 的导数推导、机器学习中的逻辑回归梯度、线性神经网络初步

**延伸**：矩阵微积分是现代最优化算法的符号引擎；它将繁琐的分量求导简化为简洁的矩阵乘法，是理解反向传播算法 (Backpropagation) 以及任何基于梯度的机器学习模型的数学必经之路

</div>

在多维空间的最优化问题中，我们经常需要寻找使标量函数（如误差、能量、成本）达到极值的方向。**矩阵微积分**（Matrix Calculus）为处理这类问题提供了一套高度紧凑的符号系统。它不仅能显著减少求导过程中的笔误，还能直接揭示最优解的代数结构。本章将确立标准的求导法则，并介绍在统计与机器学习中最常用的“迹技巧”。

---

## 47A.1 布局约定与基础导数

!!! definition "定义 47A.1 (梯度向量)"
    对于标量函数 $f(\mathbf{x})$，其关于 $n$ 维向量 $\mathbf{x}$ 的**梯度**定义为：
    - **分母布局**（常用）：$\frac{\partial f}{\partial \mathbf{x}}$ 是一个 $n \times 1$ 列向量。
    - **分子布局**：$\frac{\partial f}{\partial \mathbf{x}}$ 是一个 $1 \times n$ 行向量。
    本站默认采用**分母布局**。

!!! theorem "定理 47A.1 (线性与二次型梯度)"
    1.  $\frac{\partial (\mathbf{a}^T \mathbf{x})}{\partial \mathbf{x}} = \mathbf{a}$
    2.  $\frac{\partial (\mathbf{x}^T A \mathbf{x})}{\partial \mathbf{x}} = (A + A^T)\mathbf{x}$（若 $A$ 对称，则为 $2A\mathbf{x}$）。

---

## 47A.2 迹技巧 (Trace Tricks)

!!! technique "技术：迹函数导数"
    对于标量函数 $f(X) = \operatorname{tr}(AXB)$，其导数计算遵循以下核心公式：
    1.  $\frac{\partial \operatorname{tr}(AX)}{\partial X} = A^T$
    2.  $\frac{\partial \operatorname{tr}(X^T A X)}{\partial X} = (A + A^T)X$
    由于标量 $s$ 满足 $s = \operatorname{tr}(s)$，许多复杂的二次型导数都可以通过化为迹的形式轻松求解。

---

## 47A.3 链式法则与 Hessian

!!! definition "定义 47A.2 (Hessian 矩阵)"
    标量函数 $f$ 对向量 $\mathbf{x}$ 的二阶导数矩阵称为 **Hessian 矩阵**：
    $$\mathbf{H} = \frac{\partial^2 f}{\partial \mathbf{x} \partial \mathbf{x}^T}$$
    **意义**：Hessian 矩阵刻画了函数的曲率，其正定性决定了驻点是极小值还是极大值。

---

## 练习题

**1. [基础] 计算 $f(x) = \mathbf{x}^T \mathbf{x}$ 关于向量 $\mathbf{x}$ 的梯度。**

??? success "参考答案"
    **计算步骤：**
    1. 写成二次型形式：$f(x) = \mathbf{x}^T I \mathbf{x}$。
    2. 应用二次型梯度公式：$\frac{\partial f}{\partial \mathbf{x}} = (I + I^T)\mathbf{x} = 2I\mathbf{x}$。
    **结论**：$\nabla f = 2\mathbf{x}$。
    这与标量 $x^2$ 的导数 $2x$ 形式上完全一致。

**2. [迹运算] 求 $\frac{\partial \operatorname{tr}(A X^T)}{\partial X}$。**

??? success "参考答案"
    **步骤：**
    1. 利用迹的性质：$\operatorname{tr}(A X^T) = \operatorname{tr}(X A^T)$。
    2. 应用基本公式 $\frac{\partial \operatorname{tr}(XB)}{\partial X} = B^T$。
    3. 令 $B = A^T$。
    **结论**：导数为 $(A^T)^T = A$。

**3. [线性回归] 为最小二乘损失函数 $J(\mathbf{w}) = \|\mathbf{y} - X\mathbf{w}\|^2$ 求关于 $\mathbf{w}$ 的梯度。**

??? success "参考答案"
    **推导：**
    1. 展开范数：$J = (\mathbf{y}-X\mathbf{w})^T(\mathbf{y}-X\mathbf{w}) = \mathbf{y}^T\mathbf{y} - 2\mathbf{y}^TX\mathbf{w} + \mathbf{w}^TX^TX\mathbf{w}$。
    2. 对 $\mathbf{w}$ 求导：
       - 第一项是常数，导数为 0。
       - 第二项 $-2(X^T\mathbf{y})^T\mathbf{w}$，导数为 $-2X^T\mathbf{y}$。
       - 第三项是关于 $A=X^TX$ 的二次型，导数为 $2X^TX\mathbf{w}$。
    **结论**：$\nabla_{\mathbf{w}} J = 2X^TX\mathbf{w} - 2X^T\mathbf{y}$。
    令梯度为 0 即可得到正规方程 $X^TX\mathbf{w} = X^T\mathbf{y}$。

**4. [Hessian] 计算上题中 $J(\mathbf{w})$ 的 Hessian 矩阵。**

??? success "参考答案"
    **计算：**
    1. 对一阶梯度 $\mathbf{g} = 2X^TX\mathbf{w} - 2X^T\mathbf{y}$ 再次求导。
    2. 梯度中的第一项关于 $\mathbf{w}$ 是线性的，系数为 $2X^TX$。
    3. 第二项与 $\mathbf{w}$ 无关，导数为 0。
    **结论**：$H = 2X^TX$。由于 $X^TX$ 总是半正定的，这解释了为什么最小二乘法的损失曲面是向上开口的凸面。

**5. [行列式导数] 已知 $\frac{\partial \ln \det X}{\partial X} = (X^{-1})^T$（针对 $X \succ 0$）。求 $\frac{\partial \det X}{\partial X}$。**

??? success "参考答案"
    **利用链式法则：**
    1. 令 $y = \det X$，则 $\ln y = \ln \det X$。
    2. $\frac{\partial \ln y}{\partial X} = \frac{1}{y} \frac{\partial y}{\partial X}$。
    3. 已知 $\frac{\partial \ln y}{\partial X} = (X^{-1})^T$。
    4. 故 $\frac{\partial \det X}{\partial X} = (\det X) (X^{-1})^T$。
    **补充**：这正是伴随矩阵 $A^*$ 的转置。

**6. [验证] 若 $f = \mathbf{a}^T X \mathbf{b}$，求 $\frac{\partial f}{\partial X}$。**

??? success "参考答案"
    **迹技巧应用：**
    1. 因为 $f$ 是标量，$f = \operatorname{tr}(\mathbf{a}^T X \mathbf{b})$。
    2. 利用循环性质：$= \operatorname{tr}(\mathbf{b} \mathbf{a}^T X)$。
    3. 对 $X$ 求导得 $(\mathbf{b} \mathbf{a}^T)^T = \mathbf{a} \mathbf{b}^T$。
    **结论**：梯度是一个秩 1 矩阵。

**7. [Frobenius] 计算 $\frac{\partial \|X\|_F^2}{\partial X}$。**

??? success "参考答案"
    **计算：**
    1. $\|X\|_F^2 = \operatorname{tr}(X^T X)$。
    2. 令 $A=I$，应用 $\frac{\partial \operatorname{tr}(X^T A X)}{\partial X} = (A+A^T)X = 2X$。
    **结论**：导数为 $2X$。

**8. [分母布局] 在分母布局下，若 $y = Ax$，$y$ 对 $x$ 的导数 $\frac{\partial y}{\partial x}$ 维数是多少？**

??? success "参考答案"
    **结论：**
    是一个 **$n \times m$** 矩阵（即 $A^T$）。
    在分母布局中，如果结果是向量且变量是向量，梯度矩阵是 Jacobian 的转置。

**9. [非线性] 求 $\frac{\partial \exp(\mathbf{a}^T \mathbf{x})}{\partial \mathbf{x}}$。**

??? success "参考答案"
    **链式法则：**
    1. 设 $u = \mathbf{a}^T \mathbf{x}$。
    2. 则 $f = e^u$。
    3. $\frac{\partial f}{\partial \mathbf{x}} = \frac{df}{du} \frac{\partial u}{\partial \mathbf{x}} = e^u \cdot \mathbf{a}$。
    **结论**：$\exp(\mathbf{a}^T \mathbf{x}) \mathbf{a}$。

**10. [应用] 简述为什么矩阵微积分在深度学习的权重更新中不可或缺。**

??? success "参考答案"
    **理由：**
    神经网络包含数百万个参数，通常排列为权值矩阵 $W$。
    为了最小化损失函数 $L$，需要计算梯度 $\nabla_W L$。
    矩阵微积分允许我们直接写出如 $\nabla_W L = \delta \cdot \mathbf{x}^T$ 这种简洁的更新规则，而无需对每一个神经元的连接单独求偏导，极大提高了算法的实现效率和运行速度。

## 本章小结

矩阵微积分是将线性代数转化为优化算法的符号机器：

1.  **符号的压缩**：布局约定和迹技巧将成千上万的分量偏导数压缩为简单的矩阵外积与转置，实现了多维空间导数运算的符号化。
2.  **最优解的解析**：线性回归等经典模型的最优解，本质上都是矩阵梯度方程 $\nabla f = 0$ 的代数根，揭示了误差曲面与算子逆之间的必然联系。
3.  **智能的引擎**：作为现代梯度下降算法的底座，矩阵微积分确立了参数在感知机、逻辑回归及神经网络中流转的数学法则。
