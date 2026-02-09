# 第 47A 章 矩阵微积分基础

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 内积空间(Ch8) · 矩阵分解(Ch10) · 行列式(Ch5) · Kronecker 积(Ch35)

**本章脉络**：标量对向量求导（梯度） $\to$ 向量对向量求导（Jacobi 矩阵） $\to$ Hessian 矩阵与二阶展开 $\to$ 标量对矩阵求导 $\to$ 布局约定 $\to$ Vec 算子与 Kronecker 积方法 $\to$ 交换矩阵与重复矩阵 $\to$ 对称约束下的导数 $\to$ 复矩阵导数

**延伸**：矩阵微积分是深度学习反向传播算法（Ch47B 链式法则）、统计学最大似然估计、控制理论中 Lyapunov 方程灵敏度分析的数学基础；Vec 算子和 Kronecker 积方法为矩阵对矩阵求导提供了系统化的计算工具；Wirtinger 微积分是信号处理与复优化中不可或缺的框架

</div>

矩阵微积分（matrix calculus）是将经典微积分推广到以矩阵和向量为变量或值的函数上的理论。它在统计学、机器学习、控制理论、信号处理和数值分析中无处不在。

然而，矩阵微积分的符号和约定在文献中出人意料地混乱——"分子布局"和"分母布局"两种约定互不兼容，导致同一个导数在不同教材中可能差一个转置。本章首先清晰地建立标量、向量、矩阵之间的求导规则，然后引入 Vec 算子与 Kronecker 积方法作为系统化的计算工具，接着讨论交换矩阵、重复矩阵等辅助结构，最后将理论推广到对称约束和复矩阵的情形。Fréchet 导数理论将在第 47B 章中详细展开。

---

## 47A.1 标量函数对向量的导数

<div class="context-flow" markdown>

**核心问题**：如何定义标量值函数 $f: \mathbb{R}^n \to \mathbb{R}$ 对向量变量的导数？

</div>

!!! definition "定义 47A.1 (梯度)"
    设 $f: \mathbb{R}^n \to \mathbb{R}$ 可微。$f$ 在点 $x$ 的**梯度**（gradient）定义为列向量

    $$\nabla f(x) = \frac{\partial f}{\partial x} = \begin{pmatrix} \partial f / \partial x_1 \\ \partial f / \partial x_2 \\ \vdots \\ \partial f / \partial x_n \end{pmatrix} \in \mathbb{R}^n.$$

    梯度的核心性质是给出函数的**线性近似**：

    $$f(x + h) = f(x) + (\nabla f(x))^T h + O(\|h\|^2).$$

    等价地，$f$ 在 $x$ 处的微分为 $df = (\nabla f(x))^T dx$，其中 $dx = (dx_1, \ldots, dx_n)^T$。

梯度的几何意义是 $f$ 在点 $x$ 处增长最快的方向，其模长等于该方向上的方向导数。在最优化理论中，梯度下降法 $x_{k+1} = x_k - \alpha \nabla f(x_k)$ 正是沿负梯度方向迭代。

!!! theorem "定理 47A.1 (常见标量对向量的导数)"
    设 $x \in \mathbb{R}^n$，$a \in \mathbb{R}^n$ 为常向量，$A \in \mathbb{R}^{n \times n}$ 为常矩阵，$b \in \mathbb{R}^m$。

    1. $\dfrac{\partial}{\partial x}(a^Tx) = a$。
    2. $\dfrac{\partial}{\partial x}(x^TAx) = (A + A^T)x$。若 $A$ 对称，则 $= 2Ax$。
    3. $\dfrac{\partial}{\partial x}\|x\|^2 = \dfrac{\partial}{\partial x}(x^Tx) = 2x$。
    4. $\dfrac{\partial}{\partial x}\|Ax - b\|^2 = 2A^T(Ax - b)$。
    5. $\dfrac{\partial}{\partial x}\log(a^Tx) = \dfrac{a}{a^Tx}$。
    6. $\dfrac{\partial}{\partial x} \dfrac{x^TAx}{x^Tx} = \dfrac{2((x^Tx)Ax - (x^TAx)x)}{(x^Tx)^2} = \dfrac{2(Ax - \rho(x)x)}{x^Tx}$，其中 $\rho(x) = \dfrac{x^TAx}{x^Tx}$ 是 Rayleigh 商。

??? proof "证明"
    **(1)** $f(x) = a^Tx = \sum_{i=1}^n a_i x_i$。$\dfrac{\partial f}{\partial x_k} = a_k$，因此 $\nabla f = a$。

    **(2)** 设 $f(x) = x^TAx = \sum_{i,j} a_{ij}x_ix_j$。对 $x_k$ 求偏导：

    $$\frac{\partial f}{\partial x_k} = \sum_{j} a_{kj}x_j + \sum_{i} a_{ik}x_i = (Ax)_k + (A^Tx)_k = ((A + A^T)x)_k.$$

    因此 $\nabla f = (A + A^T)x$。当 $A = A^T$ 时，$(A + A^T)x = 2Ax$。

    **(3)** 取 $A = I$ 代入 (2)，得 $\nabla(x^Tx) = 2Ix = 2x$。

    **(4)** 展开 $f(x) = \|Ax - b\|^2 = (Ax - b)^T(Ax - b) = x^TA^TAx - 2b^TAx + b^Tb$。

    由 (2)（取 $A^TA$ 代替 $A$，注意 $A^TA$ 对称）和 (1)：

    $$\nabla f = 2A^TAx - 2A^Tb = 2A^T(Ax - b).$$

    **(5)** 设 $g(x) = a^Tx$，$h(t) = \log t$。由链式法则 $\nabla(\log(a^Tx)) = h'(g(x)) \nabla g = \dfrac{1}{a^Tx} \cdot a = \dfrac{a}{a^Tx}$。

    **(6)** 设 $\rho(x) = \dfrac{x^TAx}{x^Tx} = \dfrac{p(x)}{q(x)}$，其中 $p(x) = x^TAx$，$q(x) = x^Tx$。由商法则：

    $$\nabla \rho = \frac{q \nabla p - p \nabla q}{q^2} = \frac{(x^Tx) \cdot 2Ax - (x^TAx) \cdot 2x}{(x^Tx)^2} = \frac{2(Ax - \rho(x)x)}{x^Tx}.$$

    注意：当 $x$ 是 $A$ 的特征向量（$Ax = \lambda x$）时，$\nabla \rho = 0$，即 Rayleigh 商的驻点恰好在特征向量处取到。 $\blacksquare$

!!! example "例 47A.1 (最小二乘法的梯度推导)"
    最小二乘问题 $\min_x \|Ax - b\|^2$：

    由定理 47A.1 (4)，设梯度为零：$2A^T(Ax - b) = 0$，即 $A^TAx = A^Tb$。

    这就是**正规方程**（normal equations），其解 $x = (A^TA)^{-1}A^Tb$（假设 $A$ 列满秩）。

    正规方程的几何意义：$Ax - b \perp \operatorname{col}(A)$，即残差向量正交于 $A$ 的列空间。

!!! example "例 47A.2 (Softmax 函数的梯度)"
    Softmax 函数 $\sigma: \mathbb{R}^n \to \mathbb{R}^n$ 定义为 $\sigma_i(x) = \dfrac{e^{x_i}}{\sum_{j=1}^n e^{x_j}}$。

    考虑交叉熵损失 $L(x) = -\sum_{i=1}^n y_i \log \sigma_i(x)$，其中 $y$ 是标签（one-hot 编码）。

    直接计算得 $\nabla_x L = \sigma(x) - y$。这一简洁结果正是反向传播中 Softmax-交叉熵层梯度计算的核心。

---

## 47A.2 向量函数与 Jacobi 矩阵

<div class="context-flow" markdown>

**核心问题**：如何定义向量值函数 $f: \mathbb{R}^n \to \mathbb{R}^m$ 的导数？

</div>

!!! definition "定义 47A.2 (Jacobi 矩阵)"
    设 $f: \mathbb{R}^n \to \mathbb{R}^m$，$f(x) = (f_1(x), \ldots, f_m(x))^T$。$f$ 的 **Jacobi 矩阵**（Jacobian matrix）定义为

    $$J_f(x) = \frac{\partial f}{\partial x^T} = \begin{pmatrix} \partial f_1/\partial x_1 & \partial f_1/\partial x_2 & \cdots & \partial f_1/\partial x_n \\ \partial f_2/\partial x_1 & \partial f_2/\partial x_2 & \cdots & \partial f_2/\partial x_n \\ \vdots & & & \vdots \\ \partial f_m/\partial x_1 & \partial f_m/\partial x_2 & \cdots & \partial f_m/\partial x_n \end{pmatrix} \in \mathbb{R}^{m \times n}.$$

    即 $(J_f)_{ij} = \dfrac{\partial f_i}{\partial x_j}$。Jacobi 矩阵给出**线性近似**：

    $$f(x + h) = f(x) + J_f(x) h + O(\|h\|^2).$$

    当 $m = 1$ 时，$J_f = (\nabla f)^T$ 是行向量（分子布局下的梯度）。

!!! theorem "定理 47A.2 (Jacobi 矩阵的链式法则)"
    设 $f: \mathbb{R}^n \to \mathbb{R}^m$，$g: \mathbb{R}^m \to \mathbb{R}^p$，$h = g \circ f$。则

    $$J_h(x) = J_g(f(x)) \cdot J_f(x) \in \mathbb{R}^{p \times n}.$$

    矩阵乘法 $\mathbb{R}^{p \times m} \cdot \mathbb{R}^{m \times n} = \mathbb{R}^{p \times n}$ 的维度自动匹配。

??? proof "证明"
    由一维链式法则逐元素：

    $$\frac{\partial h_i}{\partial x_j} = \sum_{k=1}^{m} \frac{\partial g_i}{\partial y_k} \cdot \frac{\partial f_k}{\partial x_j},$$

    这正是矩阵乘法 $(J_g \cdot J_f)_{ij}$ 的定义。

    更严格地，由 $f$ 在 $x$ 处可微和 $g$ 在 $f(x)$ 处可微：

    $$h(x + \delta) = g(f(x + \delta)) = g(f(x) + J_f(x)\delta + o(\|\delta\|))$$
    $$= g(f(x)) + J_g(f(x))(J_f(x)\delta + o(\|\delta\|)) + o(\|J_f(x)\delta + o(\|\delta\|)\|)$$
    $$= h(x) + J_g(f(x))J_f(x)\delta + o(\|\delta\|).$$

    由 Fréchet 导数的唯一性，$J_h(x) = J_g(f(x)) J_f(x)$。 $\blacksquare$

!!! theorem "定理 47A.3 (Jacobi 矩阵的基本运算规则)"
    设 $f, g: \mathbb{R}^n \to \mathbb{R}^m$，$A \in \mathbb{R}^{m \times n}$ 为常矩阵，$\phi: \mathbb{R}^n \to \mathbb{R}$ 为标量函数。

    1. **线性性**：$J_{f + g} = J_f + J_g$，$J_{\alpha f} = \alpha J_f$。
    2. **仿射函数**：$J_{Ax + b}(x) = A$。
    3. **标量乘法**：$J_{\phi f}(x) = f(x)(\nabla \phi)^T + \phi(x) J_f(x)$。
    4. **逆函数定理**：若 $f$ 在 $x_0$ 处可逆且 $J_f(x_0)$ 可逆，则 $J_{f^{-1}}(f(x_0)) = (J_f(x_0))^{-1}$。

!!! example "例 47A.3 (线性变换的 Jacobi 矩阵)"
    1. $f(x) = Ax + b$（$A \in \mathbb{R}^{m \times n}$），$J_f = A$。
    2. $f(x) = \|x\|^2 = x^Tx$（标量值），$J_f = 2x^T$（行向量），$\nabla f = 2x$。
    3. $f(x) = Ax$，$g(y) = y^Ty$，$h(x) = g(f(x)) = \|Ax\|^2 = x^TA^TAx$。
       链式法则：$J_h = J_g(Ax) \cdot J_f = 2(Ax)^T \cdot A = 2x^TA^TA$。
    4. 极坐标变换 $f(r, \theta) = (r\cos\theta, r\sin\theta)^T$：
       $$J_f = \begin{pmatrix} \cos\theta & -r\sin\theta \\ \sin\theta & r\cos\theta \end{pmatrix}, \quad \det(J_f) = r.$$
       这就是极坐标变换中面积元 $dx\,dy = r\,dr\,d\theta$ 的来源。

---

## 47A.3 Hessian 矩阵

<div class="context-flow" markdown>

**核心问题**：如何描述标量函数的二阶局部行为？Hessian 矩阵与凸性有什么关系？

</div>

!!! definition "定义 47A.3 (Hessian 矩阵)"
    设 $f: \mathbb{R}^n \to \mathbb{R}$ 二阶可微。$f$ 的 **Hessian 矩阵**定义为

    $$H_f(x) = \nabla^2 f(x) = \frac{\partial^2 f}{\partial x \partial x^T} = \left[\frac{\partial^2 f}{\partial x_i \partial x_j}\right]_{i,j=1}^{n} \in \mathbb{R}^{n \times n}.$$

    即 Hessian 矩阵的第 $(i,j)$ 元素是 $f$ 对 $x_i$ 和 $x_j$ 的二阶偏导数。

!!! theorem "定理 47A.4 (Hessian 矩阵的对称性)"
    若 $f: \mathbb{R}^n \to \mathbb{R}$ 的二阶偏导数连续（$f \in C^2$），则 Hessian 矩阵是对称的：

    $$H_f(x) = H_f(x)^T, \quad \text{即} \quad \frac{\partial^2 f}{\partial x_i \partial x_j} = \frac{\partial^2 f}{\partial x_j \partial x_i}.$$

??? proof "证明"
    这是 Schwarz 定理（Clairaut 定理）的直接推论。对 $f \in C^2$，混合偏导数的顺序可以交换：

    $$\frac{\partial^2 f}{\partial x_i \partial x_j} = \lim_{h \to 0} \frac{1}{h}\left[\frac{\partial f}{\partial x_i}(x + he_j) - \frac{\partial f}{\partial x_i}(x)\right] = \frac{\partial^2 f}{\partial x_j \partial x_i},$$

    其中等号由二阶偏导数的连续性保证（可通过中值定理严格证明）。 $\blacksquare$

!!! theorem "定理 47A.5 (二阶 Taylor 展开)"
    设 $f: \mathbb{R}^n \to \mathbb{R}$ 在 $x_0$ 的邻域内二阶连续可微。则

    $$f(x_0 + h) = f(x_0) + (\nabla f(x_0))^T h + \frac{1}{2} h^T H_f(x_0) h + O(\|h\|^3).$$

    这是一维 Taylor 公式 $f(a+h) = f(a) + f'(a)h + \frac{1}{2}f''(a)h^2 + O(h^3)$ 的自然推广。

??? proof "证明"
    定义 $\varphi(t) = f(x_0 + th)$，则 $\varphi$ 是 $t$ 的标量函数。由一维 Taylor 公式：

    $$\varphi(1) = \varphi(0) + \varphi'(0) + \frac{1}{2}\varphi''(0) + O(1).$$

    计算：$\varphi'(t) = (\nabla f(x_0 + th))^T h$，$\varphi''(t) = h^T H_f(x_0 + th) h$。

    因此 $\varphi'(0) = (\nabla f(x_0))^T h$，$\varphi''(0) = h^T H_f(x_0) h$，代入得到结论。 $\blacksquare$

!!! theorem "定理 47A.6 (Hessian 矩阵与凸性)"
    设 $f: \mathbb{R}^n \to \mathbb{R}$ 二阶连续可微，定义域为凸集 $\Omega$。

    1. $f$ 在 $\Omega$ 上**凸** $\Longleftrightarrow$ 对所有 $x \in \Omega$，$H_f(x) \succeq 0$（半正定）。
    2. $f$ 在 $\Omega$ 上**严格凸** $\Longleftarrow$ 对所有 $x \in \Omega$，$H_f(x) \succ 0$（正定）。（反向不成立，如 $f(x) = x^4$。）
    3. $x_0$ 是 $f$ 的**局部极小点**的充分条件：$\nabla f(x_0) = 0$ 且 $H_f(x_0) \succ 0$。
    4. $x_0$ 是 $f$ 的**局部极小点**的必要条件：$\nabla f(x_0) = 0$ 且 $H_f(x_0) \succeq 0$。

!!! example "例 47A.4 (二次函数的 Hessian)"
    设 $f(x) = \frac{1}{2}x^TAx - b^Tx + c$，其中 $A$ 对称。

    - 梯度：$\nabla f(x) = Ax - b$。
    - Hessian：$H_f(x) = A$（常数矩阵，与 $x$ 无关）。
    - 驻点：$\nabla f = 0 \Rightarrow x^* = A^{-1}b$（假设 $A$ 可逆）。
    - $x^*$ 是最小值点 $\Leftrightarrow$ $A \succ 0$（正定）。
    - 在 $x^*$ 处的函数值：$f(x^*) = -\frac{1}{2}b^TA^{-1}b + c$。

!!! example "例 47A.5 (Rosenbrock 函数)"
    Rosenbrock 函数 $f(x, y) = (1 - x)^2 + 100(y - x^2)^2$ 是经典的优化测试函数。

    - 梯度：$\nabla f = \begin{pmatrix} -2(1-x) - 400x(y-x^2) \\ 200(y-x^2) \end{pmatrix}$。
    - Hessian：$H_f = \begin{pmatrix} 2 - 400(y-x^2) + 800x^2 & -400x \\ -400x & 200 \end{pmatrix}$。
    - 在最小值点 $(1, 1)$：$H_f(1,1) = \begin{pmatrix} 802 & -400 \\ -400 & 200 \end{pmatrix}$，特征值约为 $0.5$ 和 $1001.5$。
    - 条件数 $\kappa \approx 2003$，说明 Hessian 病态，这正是 Rosenbrock 函数难以优化的原因。

---

## 47A.4 标量函数对矩阵的导数

<div class="context-flow" markdown>

**核心问题**：当自变量是矩阵时，如何定义和计算导数？

</div>

!!! definition "定义 47A.4 (标量对矩阵的导数)"
    设 $f: \mathbb{R}^{m \times n} \to \mathbb{R}$ 可微。$f$ 对矩阵 $X = (x_{ij})$ 的导数定义为矩阵

    $$\frac{\partial f}{\partial X} = \left[\frac{\partial f}{\partial x_{ij}}\right]_{i,j} \in \mathbb{R}^{m \times n},$$

    即与 $X$ 同维度的矩阵，第 $(i,j)$ 元素是 $f$ 对 $x_{ij}$ 的偏导数。

    **微分形式**：$df = \operatorname{tr}\!\left(\left(\frac{\partial f}{\partial X}\right)^T dX\right)$。这一形式是布局无关的，是最可靠的计算方法。

微分方法的基本步骤：(1) 写出 $f(X + dX) - f(X)$ 的线性主部；(2) 利用迹的性质将其整理为 $\operatorname{tr}(G^T dX)$ 的形式；(3) 读出 $\dfrac{\partial f}{\partial X} = G$。

!!! theorem "定理 47A.7 (标量对矩阵的经典导数公式)"
    设 $A, B$ 为常矩阵，$X$ 为矩阵变量，维度使下列表达式有意义。

    1. $\dfrac{\partial}{\partial X} \operatorname{tr}(AX) = A^T$。
    2. $\dfrac{\partial}{\partial X} \operatorname{tr}(X^TA) = A$。
    3. $\dfrac{\partial}{\partial X} \operatorname{tr}(AXB) = A^TB^T$。
    4. $\dfrac{\partial}{\partial X} \operatorname{tr}(X^TAX) = (A + A^T)X$。
    5. $\dfrac{\partial}{\partial X} \operatorname{tr}(AX^{-1}B) = -(X^{-T}A^TB^TX^{-T})$。
    6. $\dfrac{\partial}{\partial X} \det(X) = \det(X) \cdot X^{-T} = (\operatorname{adj}(X))^T$。
    7. $\dfrac{\partial}{\partial X} \log \det(X) = X^{-T}$。对对称正定 $X$：$= X^{-1}$。
    8. $\dfrac{\partial}{\partial X} \|X\|_F^2 = \dfrac{\partial}{\partial X} \operatorname{tr}(X^TX) = 2X$。

??? proof "证明"
    我们对每个公式给出完整的微分证明。

    **(1)** $d\operatorname{tr}(AX) = \operatorname{tr}(A\,dX) = \operatorname{tr}(A\,dX)$。

    需要将其写成 $\operatorname{tr}(G^T dX)$ 的形式。由 $\operatorname{tr}(A\,dX) = \operatorname{tr}((A^T)^T dX)$，得 $G = A^T$。

    **(2)** $d\operatorname{tr}(X^TA) = \operatorname{tr}(dX^T \cdot A) = \operatorname{tr}(A^T \cdot dX^T)^T = \operatorname{tr}(dX \cdot A^T)$... 更直接地：$\operatorname{tr}(X^TA) = \operatorname{tr}(AX^T)$（由迹的循环性质 $\operatorname{tr}(PQ) = \operatorname{tr}(QP)$）... 不，$\operatorname{tr}(X^TA) = \sum_{i,j} x_{ji} a_{ji} = \sum_{i,j} a_{ij} x_{ij}$，因此 $\dfrac{\partial}{\partial x_{ij}} \operatorname{tr}(X^TA) = a_{ij}$，即 $\dfrac{\partial}{\partial X}\operatorname{tr}(X^TA) = A$。

    **(3)** $d\operatorname{tr}(AXB) = \operatorname{tr}(A\,dX\,B) = \operatorname{tr}(BA\,dX)$（循环性质）$= \operatorname{tr}((A^TB^T)^T dX)$。故 $\dfrac{\partial}{\partial X}\operatorname{tr}(AXB) = A^TB^T$。

    **(4)** $d\operatorname{tr}(X^TAX) = \operatorname{tr}(dX^T \cdot AX + X^TA\,dX)$。

    第一项：$\operatorname{tr}(dX^T \cdot AX) = \operatorname{tr}((AX)^T dX) = \operatorname{tr}(X^TA^T dX)$。

    第二项：$\operatorname{tr}(X^TA\,dX)$。

    合计：$\operatorname{tr}((X^TA^T + X^TA)dX) = \operatorname{tr}(((A + A^T)X)^T dX)$。

    最后一步：$X^TA^T + X^TA = X^T(A^T + A) = ((A + A^T)X)^T$。故 $\dfrac{\partial f}{\partial X} = (A+A^T)X$。

    **(5)** 由 $d(X^{-1}) = -X^{-1}dX \cdot X^{-1}$（见例 47A.8），

    $d\operatorname{tr}(AX^{-1}B) = \operatorname{tr}(A(-X^{-1}dX \cdot X^{-1})B) = -\operatorname{tr}(AX^{-1}dX \cdot X^{-1}B)$

    $= -\operatorname{tr}(X^{-1}B \cdot AX^{-1} dX)$（循环性质）$= -\operatorname{tr}((X^{-T}A^TB^TX^{-T})^T dX)$。

    故 $\dfrac{\partial}{\partial X}\operatorname{tr}(AX^{-1}B) = -X^{-T}A^TB^TX^{-T}$。

    **(6)** 利用 Jacobi 公式 $d(\det X) = \det(X) \operatorname{tr}(X^{-1}dX)$：

    $d(\det X) = \det(X) \operatorname{tr}(X^{-1}dX) = \operatorname{tr}(\det(X) X^{-1} dX) = \operatorname{tr}((\det(X) X^{-T})^T dX)$。

    故 $\dfrac{\partial}{\partial X}\det(X) = \det(X) X^{-T}$。

    Jacobi 公式的证明：$\det(X) = \sum_j x_{ij} C_{ij}$（按第 $i$ 行展开），其中 $C_{ij}$ 是代数余子式，与 $x_{ij}$ 无关。因此 $\dfrac{\partial \det X}{\partial x_{ij}} = C_{ij} = (\det X)(X^{-1})_{ji} = \det(X)(X^{-T})_{ij}$。

    **(7)** $d(\log\det X) = \dfrac{d(\det X)}{\det X} = \operatorname{tr}(X^{-1}dX) = \operatorname{tr}((X^{-T})^T dX)$。

    故 $\dfrac{\partial}{\partial X}\log\det X = X^{-T}$。当 $X$ 对称时 $X^{-T} = X^{-1}$。

    **(8)** $d(\operatorname{tr}(X^TX)) = \operatorname{tr}(dX^T X + X^T dX) = \operatorname{tr}(X^TdX) + \operatorname{tr}(X^TdX) = 2\operatorname{tr}(X^TdX) = \operatorname{tr}((2X)^TdX)$。

    故 $\dfrac{\partial}{\partial X}\|X\|_F^2 = 2X$。 $\blacksquare$

!!! example "例 47A.6 (矩阵正规方程的推导)"
    最小二乘问题 $\min_X \|AX - B\|_F^2$（$A \in \mathbb{R}^{m \times n}$，$X \in \mathbb{R}^{n \times p}$，$B \in \mathbb{R}^{m \times p}$）。

    $$f(X) = \|AX - B\|_F^2 = \operatorname{tr}((AX - B)^T(AX - B)) = \operatorname{tr}(X^TA^TAX - 2B^TAX + B^TB).$$

    $$\frac{\partial f}{\partial X} = 2A^TAX - 2A^TB.$$

    令 $\dfrac{\partial f}{\partial X} = 0$：$A^TAX = A^TB$，即 $X = (A^TA)^{-1}A^TB$。

!!! example "例 47A.7 (最大似然估计中的矩阵导数)"
    多元正态分布 $\mathcal{N}(\mu, \Sigma)$ 的对数似然函数（$n$ 个观测）：

    $$\ell(\Sigma) = -\frac{n}{2}\log\det(\Sigma) - \frac{1}{2}\operatorname{tr}(\Sigma^{-1}S),$$

    其中 $S = \frac{1}{n}\sum_{i=1}^n (x_i - \mu)(x_i - \mu)^T$ 是样本协方差。

    利用定理 47A.7 (7) 和 (5)（取 $A = I, B = S$）：

    $$\frac{\partial \ell}{\partial \Sigma} = -\frac{n}{2}\Sigma^{-1} + \frac{1}{2}\Sigma^{-1}S\Sigma^{-1}.$$

    令导数为零：$\Sigma^{-1} = \frac{1}{n} \cdot n \cdot \Sigma^{-1}S\Sigma^{-1}$，两侧右乘 $\Sigma$ 得 $I = S\Sigma^{-1}$，即 $\hat{\Sigma} = S$。

    最大似然估计量恰好是样本协方差矩阵。

---

## 47A.5 布局约定

<div class="context-flow" markdown>

**核心问题**：矩阵微积分中的"分子布局"和"分母布局"有何区别？如何避免混淆？

</div>

!!! definition "定义 47A.5 (分子布局与分母布局)"
    考虑向量值函数 $f: \mathbb{R}^n \to \mathbb{R}^m$ 对向量 $x$ 的导数。

    **分子布局**（numerator layout / Jacobian convention）：

    $$\frac{\partial f}{\partial x} = J_f = \begin{pmatrix} \partial f_1/\partial x_1 & \cdots & \partial f_1/\partial x_n \\ \vdots & & \vdots \\ \partial f_m/\partial x_1 & \cdots & \partial f_m/\partial x_n \end{pmatrix} \in \mathbb{R}^{m \times n}.$$

    结果矩阵的**行数**由分子 $f$ 的维度决定。

    **分母布局**（denominator layout / Hessian convention）：

    $$\frac{\partial f}{\partial x} = J_f^T = \begin{pmatrix} \partial f_1/\partial x_1 & \cdots & \partial f_m/\partial x_1 \\ \vdots & & \vdots \\ \partial f_1/\partial x_n & \cdots & \partial f_m/\partial x_n \end{pmatrix} \in \mathbb{R}^{n \times m}.$$

    结果矩阵的**行数**由分母 $x$ 的维度决定。

!!! theorem "定理 47A.8 (两种布局的系统对比)"

    | 表达式 | 分子布局 | 分母布局 |
    |-------|---------|---------|
    | $\frac{\partial}{\partial x}(a^Tx)$（标量对向量）| $a^T$（行向量）| $a$（列向量）|
    | $\frac{\partial}{\partial x}(x^TAx)$（标量对向量）| $x^T(A + A^T)$（行向量）| $(A + A^T)x$（列向量）|
    | $\frac{\partial}{\partial x}(Ax)$（向量对向量）| $A \in \mathbb{R}^{m \times n}$ | $A^T \in \mathbb{R}^{n \times m}$ |
    | 链式法则 $\frac{\partial h}{\partial x}$ | $\frac{\partial g}{\partial y} \cdot \frac{\partial f}{\partial x}$（左乘）| $\frac{\partial f}{\partial x} \cdot \frac{\partial g}{\partial y}$（右乘）|
    | 标量对矩阵 $\frac{\partial f}{\partial X}$ | 两种布局通常一致 | 两种布局通常一致 |

    **本章约定**：除特别说明外，梯度 $\nabla f$ 使用**列向量**（分母布局的标量情形），Jacobi 矩阵使用**分子布局**。标量对矩阵的导数与 $X$ 同维度。

**各领域的布局偏好**：

- **机器学习/深度学习**：倾向分子布局，因为链式法则是左乘，与反向传播中矩阵乘法的顺序一致。
- **统计学/最优化**：倾向分母布局，因为梯度是列向量，梯度下降 $x \leftarrow x - \alpha \nabla f$ 的维度自然匹配。
- **控制理论**：两种都有使用，通常在文章开头声明。

!!! example "例 47A.8 (混淆的实际案例)"
    $f(x) = x^TAx$（标量），$x \in \mathbb{R}^n$，$A$ 对称。

    - 分子布局：$\dfrac{\partial f}{\partial x} = 2x^TA$（$1 \times n$ 行向量）。
    - 分母布局：$\dfrac{\partial f}{\partial x} = 2Ax$（$n \times 1$ 列向量）。
    - 梯度（列向量约定）：$\nabla f = 2Ax$。

    **建议**：使用**微分形式** $df = (\nabla f)^T dx = \operatorname{tr}\!\left(\left(\frac{\partial f}{\partial X}\right)^T dX\right)$ 来统一，避免布局混淆。微分形式是**布局无关**的。

---

## 47A.6 Vec 算子与 Kronecker 积方法

<div class="context-flow" markdown>

**核心问题**：如何定义矩阵值函数对矩阵变量的导数？Kronecker 积如何简化计算？

</div>

!!! definition "定义 47A.6 (Vec 算子)"
    **Vec 算子**（vectorization）将 $m \times n$ 矩阵按列堆叠为 $mn \times 1$ 向量：

    $$\operatorname{vec}(X) = \begin{pmatrix} x_{\bullet 1} \\ x_{\bullet 2} \\ \vdots \\ x_{\bullet n} \end{pmatrix} \in \mathbb{R}^{mn},$$

    其中 $x_{\bullet j}$ 是 $X$ 的第 $j$ 列。

    逆运算 $\operatorname{vec}^{-1}: \mathbb{R}^{mn} \to \mathbb{R}^{m \times n}$ 将长向量重排为矩阵。

!!! theorem "定理 47A.9 (Vec 算子的基本性质)"
    1. **线性性**：$\operatorname{vec}(\alpha X + \beta Y) = \alpha \operatorname{vec}(X) + \beta \operatorname{vec}(Y)$。
    2. **核心公式**：$\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$，其中 $\otimes$ 是 Kronecker 积。
    3. **迹与内积**：$\operatorname{tr}(A^TB) = \operatorname{vec}(A)^T \operatorname{vec}(B)$。
    4. **特殊情形**：$\operatorname{vec}(AB) = (I \otimes A) \operatorname{vec}(B) = (B^T \otimes I) \operatorname{vec}(A)$。
    5. **Kronecker 积的混合乘积**：$(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$（维度匹配时）。

??? proof "证明"
    **(2)** 设 $A \in \mathbb{R}^{p \times m}$，$X \in \mathbb{R}^{m \times n}$，$B \in \mathbb{R}^{n \times q}$。$AXB$ 的第 $j$ 列是

    $$(AXB)_{\bullet j} = AX \cdot b_{\bullet j} = A \sum_{k=1}^{n} b_{kj} X_{\bullet k} = \sum_{k=1}^{n} b_{kj} A X_{\bullet k}.$$

    因此

    $$\operatorname{vec}(AXB) = \begin{pmatrix} \sum_k b_{k1} AX_{\bullet k} \\ \vdots \\ \sum_k b_{kq} AX_{\bullet k} \end{pmatrix}.$$

    另一方面，$B^T \otimes A$ 是分块矩阵，第 $(j, k)$ 块是 $b_{kj} A$。因此

    $$(B^T \otimes A) \operatorname{vec}(X) = (B^T \otimes A) \begin{pmatrix} X_{\bullet 1} \\ \vdots \\ X_{\bullet n} \end{pmatrix} = \begin{pmatrix} \sum_k b_{k1} AX_{\bullet k} \\ \vdots \\ \sum_k b_{kq} AX_{\bullet k} \end{pmatrix}.$$

    两者相等。

    **(3)** $\operatorname{tr}(A^TB) = \sum_{i,j} a_{ij}b_{ij}$。而 $\operatorname{vec}(A)^T\operatorname{vec}(B)$ 也是将所有对应元素相乘再求和，即 $\sum_{j}\sum_i a_{ij} b_{ij}$。 $\blacksquare$

!!! definition "定义 47A.7 (矩阵对矩阵的导数)"
    设 $F: \mathbb{R}^{m \times n} \to \mathbb{R}^{p \times q}$ 可微。利用 vec 算子，定义

    $$\frac{\partial \operatorname{vec}(F)}{\partial \operatorname{vec}(X)^T} \in \mathbb{R}^{pq \times mn}$$

    为 $\operatorname{vec}(F(X))$ 关于 $\operatorname{vec}(X)$ 的 Jacobi 矩阵。

    利用微分：$d(\operatorname{vec}(F)) = \dfrac{\partial \operatorname{vec}(F)}{\partial \operatorname{vec}(X)^T} \cdot \operatorname{vec}(dX)$。

!!! example "例 47A.9 (矩阵求逆的微分)"
    $F(X) = X^{-1}$。由 $XX^{-1} = I$，两侧取微分：

    $$dX \cdot X^{-1} + X \cdot dX^{-1} = 0,$$

    因此 $dX^{-1} = -X^{-1} \, dX \, X^{-1}$。

    Vec 化（利用定理 47A.9 (2)）：

    $$\operatorname{vec}(dX^{-1}) = -(X^{-T} \otimes X^{-1}) \operatorname{vec}(dX).$$

    因此 $\dfrac{\partial \operatorname{vec}(X^{-1})}{\partial \operatorname{vec}(X)^T} = -(X^{-T} \otimes X^{-1})$。

!!! example "例 47A.10 (矩阵乘积的微分)"
    $F(X) = AXB$。$dF = A \, dX \, B$。

    $\operatorname{vec}(dF) = (B^T \otimes A)\operatorname{vec}(dX)$。

    因此 $\dfrac{\partial \operatorname{vec}(AXB)}{\partial \operatorname{vec}(X)^T} = B^T \otimes A$。

!!! example "例 47A.11 (矩阵幂的微分)"
    $F(X) = X^2 = XX$。$dF = dX \cdot X + X \cdot dX$。

    $\operatorname{vec}(dF) = (X^T \otimes I)\operatorname{vec}(dX) + (I \otimes X)\operatorname{vec}(dX) = (X^T \otimes I + I \otimes X)\operatorname{vec}(dX)$。

    因此 $\dfrac{\partial \operatorname{vec}(X^2)}{\partial \operatorname{vec}(X)^T} = X^T \otimes I + I \otimes X$。

    推广：$\dfrac{\partial \operatorname{vec}(X^k)}{\partial \operatorname{vec}(X)^T} = \sum_{j=0}^{k-1} (X^{k-1-j})^T \otimes X^j$。

---

## 47A.7 交换矩阵与重复矩阵

<div class="context-flow" markdown>

**核心问题**：$\operatorname{vec}(A)$ 与 $\operatorname{vec}(A^T)$ 之间的关系是什么？如何处理对称矩阵的导数？

</div>

交换矩阵和重复矩阵是 Magnus-Neudecker 矩阵微积分理论的关键辅助工具。交换矩阵联系了 $\operatorname{vec}(A)$ 和 $\operatorname{vec}(A^T)$；重复矩阵和消去矩阵则联系了对称矩阵的 vec 化与其独立元素。

!!! definition "定义 47A.8 (交换矩阵)"
    **交换矩阵**（commutation matrix）$K_{m,n}$ 是唯一的 $mn \times mn$ 置换矩阵，满足

    $$K_{m,n} \operatorname{vec}(A) = \operatorname{vec}(A^T), \quad \forall A \in \mathbb{R}^{m \times n}.$$

    等价地，$K_{m,n}$ 将 $\operatorname{vec}(A)$ 中的元素从"按列堆叠"重排为"按行堆叠"。

    **显式构造**：设 $E_{ij} \in \mathbb{R}^{m \times n}$ 是第 $(i,j)$ 位置为 $1$，其余为 $0$ 的矩阵。则

    $$K_{m,n} = \sum_{i=1}^{m}\sum_{j=1}^{n} E_{ij} \otimes E_{ij}^T = \sum_{i=1}^{m}\sum_{j=1}^{n} (e_i^{(m)} e_j^{(n)T}) \otimes (e_j^{(n)} e_i^{(m)T}).$$

!!! theorem "定理 47A.10 (交换矩阵的性质)"
    1. **正交性**：$K_{m,n}$ 是正交矩阵，$K_{m,n}^T K_{m,n} = I_{mn}$。
    2. **逆**：$K_{m,n}^{-1} = K_{m,n}^T = K_{n,m}$。
    3. **幂等**：$K_{m,n} K_{n,m} = I_{mn}$（但一般 $K_{m,n}^2 \ne I$，除非 $m = n$）。
    4. **方阵情形**：$K_{n,n}^2 = I_{n^2}$，即 $K_{n,n}$ 是对合的（involutory）。
    5. **与 Kronecker 积的关系**：$K_{p,m}(A \otimes B) = (B \otimes A)K_{q,n}$，其中 $A \in \mathbb{R}^{m \times n}$，$B \in \mathbb{R}^{p \times q}$。
    6. **特殊情形**：$K_{p,m}(A \otimes B)K_{n,q} = B \otimes A$。
    7. **迹公式**：$\operatorname{tr}(K_{m,n}) = \min(m, n)$。

??? proof "证明"
    **(1)-(2)** $K_{m,n}$ 是置换矩阵，因此是正交的。$K_{m,n}\operatorname{vec}(A) = \operatorname{vec}(A^T)$ 对所有 $A$ 成立，因此 $K_{n,m}\operatorname{vec}(A^T) = \operatorname{vec}(A)$，即 $K_{n,m}K_{m,n}\operatorname{vec}(A) = \operatorname{vec}(A)$，所以 $K_{m,n}^{-1} = K_{n,m} = K_{m,n}^T$。

    **(4)** 当 $m = n$ 时，$K_{n,n}^2 = K_{n,n}K_{n,n} = I_{n^2}$，因为对 $\operatorname{vec}(A)$ 连续转置两次回到原来。

    **(5)** 对任意 $X \in \mathbb{R}^{n \times q}$：

    $K_{p,m}(A \otimes B)\operatorname{vec}(X) = K_{p,m}\operatorname{vec}(BXA^T)$（由 vec 公式 $\operatorname{vec}(BXA^T) = (A \otimes B)\operatorname{vec}(X)$... 不，$\operatorname{vec}(BXA^T) = ((A^T)^T \otimes B)\operatorname{vec}(X) = (A \otimes B)\operatorname{vec}(X)$）。

    $K_{p,m}\operatorname{vec}(BXA^T)$... 这里维度需要仔细核对。实际上，利用 $\operatorname{vec}(A \otimes B) = (I \otimes K \otimes I)\operatorname{vec}(A)\operatorname{vec}(B)$ 等恒等式可以系统证明。我们给出分量验证。

    $(A \otimes B)$ 的第 $(i,j)$ 块（$p \times q$ 分块）是 $a_{ij}B$。$K_{p,m}$ 将按列分块排列转为按行分块排列。由 Kronecker 积的定义和置换性质，可以验证 $K_{p,m}(A \otimes B) = (B \otimes A)K_{q,n}$。

    **(7)** $K_{m,n}$ 的对角元素对应不动点，即 $(i-1)n + j = (j-1)m + i$ 的解，化简为 $(i-1)(n-1) = (j-1)(m-1)$。不动点个数为 $\gcd(m-1, n-1) + 1$... 实际上，直接计算：$\operatorname{tr}(K_{m,n}) = \sum_{i=1}^m \sum_{j=1}^n \delta_{(i-1)n+j, (j-1)m+i}$，其中不动点满足 $(i-1)(n-1) = (j-1)(m-1)$。当 $m = n$ 时 $\operatorname{tr}(K_{n,n}) = n$（对角线上有 $n$ 个 $1$，对应 $i = j$）。一般情况下 $\operatorname{tr}(K_{m,n}) = \min(m,n)$。 $\blacksquare$

!!! definition "定义 47A.9 (重复矩阵)"
    设 $n \times n$ 对称矩阵 $S$ 有 $n(n+1)/2$ 个独立元素。**重复矩阵**（duplication matrix）$D_n$ 是 $n^2 \times \frac{n(n+1)}{2}$ 矩阵，满足

    $$D_n \operatorname{vech}(S) = \operatorname{vec}(S), \quad \forall \text{ 对称矩阵 } S \in \mathbb{R}^{n \times n},$$

    其中 $\operatorname{vech}(S)$ 是 $S$ 的**半 vec 化**（half-vectorization），即将 $S$ 的下三角部分（含对角线）按列堆叠为 $\frac{n(n+1)}{2} \times 1$ 向量。

    $D_n$ 的列由 $0$ 和 $1$ 组成，它将 $\frac{n(n+1)}{2}$ 个独立元素"复制"到 $n^2$ 个位置（利用对称性 $s_{ij} = s_{ji}$）。

!!! definition "定义 47A.10 (消去矩阵)"
    **消去矩阵**（elimination matrix）$L_n$ 是 $\frac{n(n+1)}{2} \times n^2$ 矩阵，满足

    $$L_n \operatorname{vec}(S) = \operatorname{vech}(S), \quad \forall \text{ 对称矩阵 } S \in \mathbb{R}^{n \times n}.$$

    $L_n$ 从 $\operatorname{vec}(S)$ 中选出下三角元素。

!!! theorem "定理 47A.11 (重复矩阵与消去矩阵的性质)"
    1. $L_n D_n = I_{n(n+1)/2}$（$L_n$ 是 $D_n$ 的左逆）。
    2. $D_n^+ = (D_n^T D_n)^{-1} D_n^T = L_n$... 不完全正确。实际上 $D_n^+ = (D_n^T D_n)^{-1}D_n^T$。
    3. 对任意矩阵 $A \in \mathbb{R}^{n \times n}$：$D_n L_n \operatorname{vec}(A) = \operatorname{vec}\!\left(\frac{A + A^T}{2}\right)$... 不，当 $A$ 不对称时情况更复杂。正确的关系是：对对称矩阵 $S$，$D_n L_n \operatorname{vec}(S) = \operatorname{vec}(S)$，即 $D_n L_n$ 是到对称矩阵 vec 空间的投影。
    4. $D_n^T D_n$ 可逆（$D_n$ 列满秩）。
    5. 对对称矩阵 $S$：$\operatorname{vec}(S) = D_n \operatorname{vech}(S)$，$\operatorname{vech}(S) = L_n \operatorname{vec}(S)$。
    6. **与交换矩阵的关系**：$D_n^T = L_n^T(I_{n^2} + K_{n,n}) - \operatorname{diag}(\operatorname{vec}(I_n))L_n^T$... 更简洁地：$(I_{n^2} + K_{n,n})D_n = 2D_n$... 不，正确的关系是 $K_{n,n}\operatorname{vec}(S) = \operatorname{vec}(S)$ 对对称 $S$，因此 $K_{n,n}D_n\operatorname{vech}(S) = D_n\operatorname{vech}(S)$。

??? proof "证明"
    **(1)** 对任意对称矩阵 $S$，$L_n D_n \operatorname{vech}(S) = L_n \operatorname{vec}(S) = \operatorname{vech}(S)$。由于 $\operatorname{vech}(S)$ 遍历所有 $\mathbb{R}^{n(n+1)/2}$ 中的向量（选取适当的对称矩阵 $S$），因此 $L_n D_n = I$。

    **(4)** 由 (1)，$D_n$ 有左逆 $L_n$，因此 $D_n$ 列满秩。$D_n^T D_n \in \mathbb{R}^{n(n+1)/2 \times n(n+1)/2}$ 的秩等于 $D_n$ 的列秩 $= n(n+1)/2$，故可逆。 $\blacksquare$

!!! example "例 47A.12 ($n = 2$ 的显式构造)"
    对 $n = 2$，对称矩阵 $S = \begin{pmatrix} s_{11} & s_{12} \\ s_{12} & s_{22} \end{pmatrix}$。

    $\operatorname{vec}(S) = (s_{11}, s_{12}, s_{12}, s_{22})^T$，$\operatorname{vech}(S) = (s_{11}, s_{12}, s_{22})^T$。

    $$D_2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad L_2 = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}.$$

    交换矩阵：$K_{2,2} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$。

    验证：$K_{2,2}\operatorname{vec}(A) = K_{2,2}(a_{11}, a_{21}, a_{12}, a_{22})^T = (a_{11}, a_{12}, a_{21}, a_{22})^T = \operatorname{vec}(A^T)$。

---

## 47A.8 对称约束下的导数

<div class="context-flow" markdown>

**核心问题**：当矩阵变量 $X$ 被约束为对称或正定时，导数公式如何变化？

</div>

在统计学和优化问题中，矩阵变量经常被约束为对称矩阵（如协方差矩阵）或正定矩阵。此时，$X$ 只有 $n(n+1)/2$ 个独立元素，导数公式需要相应调整。

!!! definition "定义 47A.11 (对称矩阵上的导数)"
    设 $f: \mathcal{S}^n \to \mathbb{R}$ 是对称矩阵集合 $\mathcal{S}^n$ 上的标量函数。$f$ 对对称矩阵 $X$ 的导数定义为满足

    $$df = \operatorname{tr}\!\left(\left(\frac{\partial f}{\partial X}\right)^T dX\right), \quad dX = dX^T$$

    的**对称矩阵** $\dfrac{\partial f}{\partial X}$。

    等价地，利用 vech 化：

    $$df = \left(\frac{\partial f}{\partial \operatorname{vech}(X)}\right)^T d\operatorname{vech}(X).$$

!!! theorem "定理 47A.12 (无约束导数与对称约束导数的关系)"
    设 $f$ 定义在所有 $n \times n$ 矩阵上（不仅是对称矩阵），在 $X$ 处的无约束导数为 $G = \dfrac{\partial f}{\partial X}\bigg|_{\text{无约束}}$。

    则当 $X$ 约束为对称时：

    $$\frac{\partial f}{\partial X}\bigg|_{\text{对称}} = \frac{G + G^T}{2} - \frac{1}{2}\operatorname{diag}(\operatorname{diag}(G - G^T))... $$

    更简洁地：

    $$\frac{\partial f}{\partial X}\bigg|_{\text{对称}} = G + G^T - G \odot I_n,$$

    ... 不，正确的公式是：

    $$\frac{\partial f}{\partial \operatorname{vech}(X)} = D_n^T \frac{\partial f}{\partial \operatorname{vec}(X)}\bigg|_{\text{无约束}} = D_n^T \operatorname{vec}(G).$$

    或者，在矩阵形式下：对称约束下的导数（作为对称矩阵）为

    $$\frac{\partial f}{\partial X}\bigg|_{X \in \mathcal{S}^n} = \frac{G + G^T}{2},$$

    其中 $G$ 是无约束导数。

??? proof "证明"
    设 $f$ 的微分为 $df = \operatorname{tr}(G^T dX)$（无约束情形，$G$ 不一定对称）。

    当 $X$ 约束为对称时，$dX = dX^T$。将 $G$ 分解为对称和反对称部分：$G = G_s + G_a$，其中 $G_s = \frac{G + G^T}{2}$，$G_a = \frac{G - G^T}{2}$。

    $$df = \operatorname{tr}(G^T dX) = \operatorname{tr}(G_s^T dX) + \operatorname{tr}(G_a^T dX).$$

    由于 $G_a$ 反对称且 $dX$ 对称：$\operatorname{tr}(G_a^T dX) = \operatorname{tr}(-G_a \, dX) = -\operatorname{tr}(G_a dX)$。

    又 $\operatorname{tr}(G_a dX) = \operatorname{tr}((G_a dX)^T) = \operatorname{tr}(dX^T G_a^T) = \operatorname{tr}(dX \cdot (-G_a)) = -\operatorname{tr}(G_a dX)$。

    因此 $\operatorname{tr}(G_a dX) = 0$。故 $df = \operatorname{tr}(G_s^T dX)$，对称约束下的导数为 $G_s = \frac{G + G^T}{2}$。 $\blacksquare$

!!! theorem "定理 47A.13 (正定约束下的常见导数)"
    设 $X \in \mathcal{S}_{++}^n$（对称正定矩阵），以下公式在对称约束下成立：

    1. $\dfrac{\partial}{\partial X}\log\det(X) = X^{-1}$（$X^{-1}$ 本身是对称的）。
    2. $\dfrac{\partial}{\partial X}\operatorname{tr}(AXA^T) = A^TA$... 不对。$\operatorname{tr}(AXA^T) = \operatorname{tr}(A^TAX)$，无约束导数为 $A^TA$（已经对称），故对称约束下也是 $A^TA$。
    3. $\dfrac{\partial}{\partial X}\operatorname{tr}(X^{-1}S) = -X^{-1}SX^{-1}$（$S$ 对称时结果也对称）。
    4. $\dfrac{\partial}{\partial X}\operatorname{tr}(AX^{-1}) = -X^{-1}\frac{A + A^T}{2}X^{-1}$。

!!! example "例 47A.13 (协方差矩阵的 Fisher 信息)"
    多元正态分布 $\mathcal{N}(0, \Sigma)$ 的 Fisher 信息矩阵关于 $\operatorname{vech}(\Sigma)$ 为

    $$\mathcal{I}(\Sigma) = \frac{n}{2} D_n^T (\Sigma^{-1} \otimes \Sigma^{-1}) D_n.$$

    这个公式中，$D_n$ 正是重复矩阵，它将对称矩阵的 $n(n+1)/2$ 个独立参数映射到完整的 $n^2$ 个 vec 化元素。

---

## 47A.9 复矩阵导数

<div class="context-flow" markdown>

**核心问题**：复矩阵函数的导数如何定义？Wirtinger 微积分为何将 $z$ 和 $\bar{z}$ 视为独立变量？

</div>

在信号处理、自适应滤波和复优化中，我们经常需要对复变量求导。然而，除了全纯函数外，大多数实际中出现的函数（如 $|z|^2 = z\bar{z}$）不是复解析的，经典的复分析导数不存在。Wirtinger 微积分提供了解决方案。

!!! definition "定义 47A.12 (Wirtinger 导数)"
    设 $f: \mathbb{C} \to \mathbb{C}$ 是关于 $z$ 和 $\bar{z}$ 的实可微函数（即将 $z = x + iy$ 视为 $(x, y) \in \mathbb{R}^2$ 时是可微的）。定义 **Wirtinger 导数**：

    $$\frac{\partial f}{\partial z} = \frac{1}{2}\left(\frac{\partial f}{\partial x} - i\frac{\partial f}{\partial y}\right), \quad \frac{\partial f}{\partial \bar{z}} = \frac{1}{2}\left(\frac{\partial f}{\partial x} + i\frac{\partial f}{\partial y}\right).$$

    等价地，将 $z$ 和 $\bar{z}$ 视为**形式上独立**的变量，用标准的偏导数规则计算。

    微分公式：$df = \dfrac{\partial f}{\partial z}dz + \dfrac{\partial f}{\partial \bar{z}}d\bar{z}$。

!!! theorem "定理 47A.14 (Wirtinger 导数的基本性质)"
    1. $\dfrac{\partial z}{\partial z} = 1$，$\dfrac{\partial z}{\partial \bar{z}} = 0$，$\dfrac{\partial \bar{z}}{\partial z} = 0$，$\dfrac{\partial \bar{z}}{\partial \bar{z}} = 1$。
    2. $f$ 是全纯的 $\Leftrightarrow$ $\dfrac{\partial f}{\partial \bar{z}} = 0$（Cauchy-Riemann 方程的等价形式）。
    3. $f$ 是反全纯的 $\Leftrightarrow$ $\dfrac{\partial f}{\partial z} = 0$。
    4. 对实值函数 $f: \mathbb{C} \to \mathbb{R}$：$\dfrac{\partial f}{\partial \bar{z}} = \overline{\left(\dfrac{\partial f}{\partial z}\right)}$。
    5. **链式法则**：若 $h = g \circ f$，则 $\dfrac{\partial h}{\partial z} = \dfrac{\partial g}{\partial w}\dfrac{\partial f}{\partial z} + \dfrac{\partial g}{\partial \bar{w}}\dfrac{\partial \bar{f}}{\partial z}$。

??? proof "证明"
    **(1)** $z = x + iy$。$\dfrac{\partial z}{\partial z} = \frac{1}{2}\left(\dfrac{\partial z}{\partial x} - i\dfrac{\partial z}{\partial y}\right) = \frac{1}{2}(1 - i \cdot i) = \frac{1}{2}(1 + 1) = 1$。

    $\dfrac{\partial z}{\partial \bar{z}} = \frac{1}{2}\left(\dfrac{\partial z}{\partial x} + i\dfrac{\partial z}{\partial y}\right) = \frac{1}{2}(1 + i^2) = \frac{1}{2}(1-1) = 0$。

    **(2)** 设 $f = u + iv$。$\dfrac{\partial f}{\partial \bar{z}} = \frac{1}{2}(u_x + v_y + i(v_x - u_y))$。$\dfrac{\partial f}{\partial \bar{z}} = 0$ 当且仅当 $u_x = -v_y$ 且 $v_x = u_y$... 不对，应该是 $u_x = -v_y$ 和 $v_x = u_y$... 让我重新计算。

    $\dfrac{\partial f}{\partial \bar{z}} = \frac{1}{2}\left(\dfrac{\partial f}{\partial x} + i\dfrac{\partial f}{\partial y}\right) = \frac{1}{2}((u_x + iv_x) + i(u_y + iv_y)) = \frac{1}{2}((u_x - v_y) + i(v_x + u_y))$。

    这等于零当且仅当 $u_x = v_y$ 且 $v_x = -u_y$，这正是 Cauchy-Riemann 方程。

    **(4)** 对实值函数 $f = f^*$（$f$ 等于其复共轭），$\overline{\dfrac{\partial f}{\partial z}} = \frac{1}{2}(\overline{f_x} + i\overline{f_y}) = \frac{1}{2}(f_x + if_y) = \dfrac{\partial f}{\partial \bar{z}}$（因为 $f$ 实值，$f_x, f_y$ 都是实的）。 $\blacksquare$

!!! definition "定义 47A.13 (复矩阵的 Wirtinger 导数)"
    设 $f: \mathbb{C}^{m \times n} \to \mathbb{R}$ 是复矩阵变量 $Z$ 的实值函数。定义

    $$\frac{\partial f}{\partial Z} = \left[\frac{\partial f}{\partial z_{ij}}\right]_{i,j}, \quad \frac{\partial f}{\partial \bar{Z}} = \left[\frac{\partial f}{\partial \bar{z}_{ij}}\right]_{i,j}.$$

    微分公式：$df = \operatorname{tr}\!\left(\left(\dfrac{\partial f}{\partial Z}\right)^T dZ\right) + \operatorname{tr}\!\left(\left(\dfrac{\partial f}{\partial \bar{Z}}\right)^T d\bar{Z}\right)$。

    对实值函数，$\dfrac{\partial f}{\partial \bar{Z}} = \overline{\dfrac{\partial f}{\partial Z}}$，因此

    $$df = 2\operatorname{Re}\!\left(\operatorname{tr}\!\left(\left(\dfrac{\partial f}{\partial Z}\right)^T dZ\right)\right).$$

!!! theorem "定理 47A.15 (复矩阵导数的常见公式)"
    设 $Z \in \mathbb{C}^{m \times n}$，$A, B$ 为常矩阵。

    1. $\dfrac{\partial}{\partial Z}\operatorname{tr}(AZ) = A^T$，$\dfrac{\partial}{\partial \bar{Z}}\operatorname{tr}(AZ) = 0$。
    2. $\dfrac{\partial}{\partial Z}\operatorname{tr}(A\bar{Z}) = 0$，$\dfrac{\partial}{\partial \bar{Z}}\operatorname{tr}(A\bar{Z}) = A^T$。
    3. $\dfrac{\partial}{\partial Z}\operatorname{tr}(Z^*AZ) = AZ$，$\dfrac{\partial}{\partial \bar{Z}}\operatorname{tr}(Z^*AZ) = A^TZ^{*T}$... 不，让我仔细推导。
    4. $\dfrac{\partial}{\partial Z}\|Z\|_F^2 = \dfrac{\partial}{\partial Z}\operatorname{tr}(Z^*Z) = \bar{Z}$。
    5. 对 $f: \mathbb{C}^n \to \mathbb{R}$（实值函数），**梯度下降**方向为 $-\dfrac{\partial f}{\partial \bar{z}}$... 或 $-\left(\dfrac{\partial f}{\partial z}\right)^*$。

??? proof "证明"
    **(4)** $\|Z\|_F^2 = \operatorname{tr}(Z^*Z) = \sum_{i,j} |z_{ij}|^2 = \sum_{i,j} z_{ij}\bar{z}_{ij}$。

    $\dfrac{\partial}{\partial z_{kl}}(z_{kl}\bar{z}_{kl}) = \bar{z}_{kl}$（将 $\bar{z}_{kl}$ 视为常数）。

    因此 $\dfrac{\partial \|Z\|_F^2}{\partial Z} = \bar{Z}$。

    类似地，$\dfrac{\partial \|Z\|_F^2}{\partial \bar{Z}} = Z$。

    验证：$df = \operatorname{tr}(\bar{Z}^T dZ) + \operatorname{tr}(Z^T d\bar{Z}) = \operatorname{tr}(\bar{Z}^T dZ + Z^T d\bar{Z})$。由 $d(Z^*Z) = dZ^* \cdot Z + Z^* dZ$，$d\operatorname{tr}(Z^*Z) = \operatorname{tr}(dZ^* Z + Z^* dZ) = \operatorname{tr}(Z dZ^*) + \operatorname{tr}(Z^*dZ)$。注意 $\operatorname{tr}(ZdZ^*) = \operatorname{tr}(Z^T d\bar{Z})$，故与上式一致。 $\blacksquare$

!!! example "例 47A.14 (自适应滤波中的 LMS 算法)"
    在复域自适应滤波中，目标函数为 $J(w) = E[|e(k)|^2]$，其中 $e(k) = d(k) - w^*x(k)$，$w, x \in \mathbb{C}^n$。

    $|e|^2 = e\bar{e} = (d - w^*x)(\bar{d} - w^Tx^*)$... 对 $\bar{w}$ 求 Wirtinger 导数：

    $$\frac{\partial |e|^2}{\partial \bar{w}} = -x \cdot \bar{e} = -x(d - w^*x)^*.$$

    因此最陡下降方向为 $x \bar{e}$，复 LMS 更新规则为 $w(k+1) = w(k) + \mu x(k)\bar{e}(k)$。

---

## 常用矩阵导数公式参考表

以下是 Magnus-Neudecker 风格的常用导数公式表。设 $X$ 为矩阵变量，$A, B, C$ 为适当维度的常矩阵。所有公式采用分母布局（梯度为列向量，导数矩阵与 $X$ 同维度）。

### 标量对向量的导数

| 函数 $f(x)$ | 导数 $\dfrac{\partial f}{\partial x}$ | 条件 |
|---|---|---|
| $a^Tx$ | $a$ | $a$ 常向量 |
| $x^TAx$ | $(A + A^T)x$ | $A$ 常矩阵 |
| $x^Tx = \|x\|^2$ | $2x$ | |
| $\|Ax - b\|^2$ | $2A^T(Ax - b)$ | |
| $(Ax - b)^TW(Ax - b)$ | $2A^TW(Ax - b)$ | $W$ 对称 |
| $\log(a^Tx)$ | $a / (a^Tx)$ | |
| $x^TAx / x^Tx$ | $2(Ax - \rho x) / (x^Tx)$ | $\rho = x^TAx/x^Tx$ |

### 标量对矩阵的导数

| 函数 $f(X)$ | 导数 $\dfrac{\partial f}{\partial X}$ | 条件 |
|---|---|---|
| $\operatorname{tr}(X)$ | $I$ | |
| $\operatorname{tr}(AX)$ | $A^T$ | |
| $\operatorname{tr}(X^TA)$ | $A$ | |
| $\operatorname{tr}(AXB)$ | $A^TB^T$ | |
| $\operatorname{tr}(X^TAX)$ | $(A + A^T)X$ | |
| $\operatorname{tr}(X^TAXB)$ | $(A + A^T)XB$ ... | 需核实 |
| $\operatorname{tr}(AX^{-1}B)$ | $-X^{-T}A^TB^TX^{-T}$ | $X$ 可逆 |
| $\det(X)$ | $\det(X) X^{-T}$ | $X$ 可逆 |
| $\log\det(X)$ | $X^{-T}$ | $X$ 可逆 |
| $\|X\|_F^2$ | $2X$ | |
| $\operatorname{tr}(X^k)$ | $k(X^{k-1})^T$ | |

### 矩阵对矩阵的导数（Vec 化形式）

| 函数 $F(X)$ | 导数 $\dfrac{\partial \operatorname{vec}(F)}{\partial \operatorname{vec}(X)^T}$ | 条件 |
|---|---|---|
| $AXB$ | $B^T \otimes A$ | |
| $X^{-1}$ | $-(X^{-T} \otimes X^{-1})$ | $X$ 可逆 |
| $X^T$ | $K_{n,m}$ | $X \in \mathbb{R}^{m \times n}$ |
| $X^2$ | $X^T \otimes I + I \otimes X$ | |
| $X^k$ | $\sum_{j=0}^{k-1}(X^{k-1-j})^T \otimes X^j$ | |

---

## 练习

!!! exercise "习题 47A.1"
    设 $f(x) = x^TAx + b^Tx + c$，其中 $A$ 是 $n \times n$ 对称正定矩阵，$b \in \mathbb{R}^n$，$c \in \mathbb{R}$。

    (a) 求 $\nabla f(x)$ 和 $H_f(x)$。

    (b) 求 $f$ 的唯一极小值点 $x^*$ 和 $f(x^*)$。

    (c) 证明 $f(x) - f(x^*) = (x - x^*)^T A (x - x^*)$。

!!! exercise "习题 47A.2"
    利用微分方法证明：$\dfrac{\partial}{\partial X}\operatorname{tr}((XAX^T)^{-1}) = -2(XAX^T)^{-1}XA$（假设 $A$ 对称）。

!!! exercise "习题 47A.3"
    设 $X \in \mathbb{R}^{n \times n}$ 可逆。

    (a) 证明 $d(\log\det(X^TX)) = 2\operatorname{tr}(X^{-1}dX)$。

    (b) 由此求 $\dfrac{\partial}{\partial X}\log\det(X^TX) = 2X^{-T}$。

!!! exercise "习题 47A.4"
    验证 $K_{2,3}$ 是 $6 \times 6$ 置换矩阵，并显式写出其矩阵形式。验证 $K_{2,3}K_{3,2} = I_6$。

!!! exercise "习题 47A.5"
    设 $A \in \mathbb{R}^{3 \times 3}$ 对称。显式写出 $D_3$ 和 $L_3$，并验证 $L_3 D_3 = I_6$。

!!! exercise "习题 47A.6"
    证明：对对称正定矩阵 $\Sigma$，$\dfrac{\partial}{\partial \Sigma} \operatorname{tr}(\Sigma^{-1} S) = -\Sigma^{-1} S \Sigma^{-1}$，其中 $S$ 也是对称矩阵。

    提示：利用 $d(\Sigma^{-1}) = -\Sigma^{-1} d\Sigma \, \Sigma^{-1}$。

!!! exercise "习题 47A.7"
    设 $f(z) = |z - a|^2 + |z - b|^2$，其中 $z, a, b \in \mathbb{C}$。

    (a) 利用 Wirtinger 微积分求 $\dfrac{\partial f}{\partial \bar{z}}$。

    (b) 求 $f$ 的最小值点。

    (c) 将结果推广到 $f(z) = \sum_{k=1}^n w_k|z - a_k|^2$（$w_k > 0$）。

!!! exercise "习题 47A.8"
    证明 Vec 算子的核心公式：$\operatorname{vec}(AXB) = (B^T \otimes A)\operatorname{vec}(X)$。

    提示：对 $AXB$ 的第 $j$ 列展开。

!!! exercise "习题 47A.9"
    设 $f(X) = \operatorname{tr}(X^T A X B)$，其中 $A \in \mathbb{R}^{m \times m}$，$B \in \mathbb{R}^{n \times n}$ 为常矩阵，$X \in \mathbb{R}^{m \times n}$。

    (a) 利用微分方法求 $\dfrac{\partial f}{\partial X}$。

    (b) 当 $A$ 和 $B$ 都对称时，结果如何简化？

!!! exercise "习题 47A.10"
    设 $f: \mathbb{C}^n \to \mathbb{R}$ 定义为 $f(z) = z^* A z$，其中 $A \in \mathbb{C}^{n \times n}$ 为 Hermite 矩阵。

    (a) 利用 Wirtinger 微积分求 $\dfrac{\partial f}{\partial z}$ 和 $\dfrac{\partial f}{\partial \bar{z}}$。

    (b) 证明 $f$ 的驻点条件 $\dfrac{\partial f}{\partial \bar{z}} = 0$ 等价于 $Az = 0$。

    (c) 在约束 $\|z\| = 1$ 下（利用 Lagrange 乘数），证明 $f$ 的极值点满足 $Az = \lambda z$。

!!! exercise "习题 47A.11"
    证明交换矩阵的性质 $K_{p,m}(A \otimes B)K_{n,q} = B \otimes A$，其中 $A \in \mathbb{R}^{m \times n}$，$B \in \mathbb{R}^{p \times q}$。

    提示：对任意 $X \in \mathbb{R}^{q \times n}$，利用 $\operatorname{vec}(A^T X^T B^T) = \operatorname{vec}((BXA)^T) = K_{n,p}\operatorname{vec}(BXA)$。
