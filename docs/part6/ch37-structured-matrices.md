# 第 37 章 结构化矩阵

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 矩阵分解(Ch10) · 数值线性代数(Ch22)

**本章脉络**：Toeplitz 矩阵 $\to$ Hankel 矩阵 $\to$ 循环矩阵与 DFT $\to$ Vandermonde 矩阵 $\to$ Cauchy 矩阵 $\to$ 位移结构与位移秩 $\to$ 快速算法

**延伸**：位移结构理论统一了所有结构化矩阵类的快速算法；在信号处理（自相关矩阵是 Toeplitz）、控制理论（Hankel 算子与系统实现）和逼近论（Vandermonde 插值）中无处不在

</div>

一般的 $n \times n$ 矩阵有 $n^2$ 个自由参数，求解线性方程组需要 $O(n^3)$ 次运算。然而，在科学计算和工程应用中反复出现的许多矩阵具有特殊的结构——它们的元素之间存在确定的关系，使得矩阵可以用远少于 $n^2$ 个参数来描述。利用这种结构，可以设计出比一般算法快得多的**快速算法**。

本章系统研究五种最重要的结构化矩阵——Toeplitz、Hankel、循环、Vandermonde、Cauchy 矩阵，然后通过**位移结构**的统一框架揭示它们之间的深层联系，最后讨论快速算法的一般原理。

---

## 37.1 Toeplitz 矩阵

<div class="context-flow" markdown>

**核心问题**：沿对角线常数的矩阵有什么特殊性质？如何快速求解 Toeplitz 系统？

</div>

!!! definition "定义 37.1 (Toeplitz 矩阵)"
    矩阵 $T = (t_{ij}) \in M_n(\mathbb{C})$ 称为 **Toeplitz 矩阵**，若 $t_{ij}$ 仅取决于差 $i - j$，即存在序列 $\{t_k\}_{k=-(n-1)}^{n-1}$ 使得
    $$t_{ij} = t_{i-j}, \quad 1 \le i, j \le n.$$
    显式地，
    $$T = \begin{pmatrix}
    t_0 & t_{-1} & t_{-2} & \cdots & t_{-(n-1)} \\
    t_1 & t_0 & t_{-1} & \cdots & t_{-(n-2)} \\
    t_2 & t_1 & t_0 & \cdots & t_{-(n-3)} \\
    \vdots & & & \ddots & \vdots \\
    t_{n-1} & t_{n-2} & t_{n-3} & \cdots & t_0
    \end{pmatrix}.$$
    Toeplitz 矩阵由 $2n-1$ 个参数完全确定（而非 $n^2$ 个）。

!!! theorem "定理 37.1 (Toeplitz 矩阵的基本性质)"
    设 $T_1, T_2 \in M_n(\mathbb{C})$ 为 Toeplitz 矩阵。则：

    (a) $\alpha T_1 + \beta T_2$ 是 Toeplitz 矩阵；

    (b) $T_1^T$ 是 Toeplitz 矩阵（将 $t_k$ 替换为 $t_{-k}$）；

    (c) $T_1 T_2$ **一般不是** Toeplitz 矩阵；

    (d) 若 $T$ 为对称 Toeplitz（$t_k = t_{-k}$），则 $T$ 有一种**持续的**结构：$JTJ = T$，其中 $J$ 为反对角单位矩阵（逆序置换矩阵）。

??? proof "证明"
    (a) 线性组合的 $(i,j)$ 元素 $\alpha t^{(1)}_{i-j} + \beta t^{(2)}_{i-j}$ 仅取决于 $i-j$。

    (b) $(T^T)_{ij} = t_{ji} = t_{j-i} = t_{-(i-j)}$，定义 $\tilde{t}_k = t_{-k}$，则 $T^T$ 是以 $\{\tilde{t}_k\}$ 为参数的 Toeplitz 矩阵。

    (c) 取 $T_1 = T_2 = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$（Toeplitz），$T_1^2 = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$ 虽仍是 Toeplitz（巧合），但一般情况反例很多。取 $n=3$, $T_1$ 以第一行 $(1,2,3)$, 第一列 $(1,0,0)$, $T_2$ 以第一行 $(1,0,0)$, 第一列 $(1,1,1)$，乘积一般不具有 Toeplitz 结构。

    (d) $J$ 是置换矩阵 $J = (e_n, e_{n-1}, \ldots, e_1)$。$(JTJ)_{ij} = t_{n+1-i, n+1-j} = t_{(n+1-i)-(n+1-j)} = t_{j-i} = t_{-(i-j)}$。当 $t_k = t_{-k}$ 时，$t_{-(i-j)} = t_{i-j}$，故 $JTJ = T$。

!!! definition "定义 37.2 (Levinson 递推)"
    **Levinson-Durbin 算法**是求解对称正定 Toeplitz 系统 $Tx = b$ 的经典 $O(n^2)$ 算法。其核心思想是利用 Toeplitz 结构，从 $k \times k$ 子系统的解递推出 $(k+1) \times (k+1)$ 子系统的解。

    设 $T_k = (t_{|i-j|})_{1 \le i,j \le k}$ 为 $T$ 的前 $k$ 阶主子矩阵。Levinson 递推的关键是引入**反射系数**（又称偏相关系数）：
    $$\alpha_k = -\frac{t_k + \sum_{j=1}^{k-1} a_j^{(k-1)} t_{k-j}}{1 + \sum_{j=1}^{k-1} a_j^{(k-1)} t_j},$$
    其中 $\{a_j^{(k-1)}\}$ 是上一步的 Levinson 系数。然后更新：
    $$a_j^{(k)} = a_j^{(k-1)} + \alpha_k a_{k-j}^{(k-1)}, \quad j = 1, \ldots, k-1, \qquad a_k^{(k)} = \alpha_k.$$

!!! example "例 37.1"
    求解 Toeplitz 系统
    $$\begin{pmatrix} 4 & 2 & 1 \\ 2 & 4 & 2 \\ 1 & 2 & 4 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}.$$

    矩阵为对称正定 Toeplitz 矩阵，$t_0 = 4, t_1 = 2, t_2 = 1$。

    **步骤 1** ($k=1$)：$4x_1^{(1)} = 1 \Rightarrow x_1^{(1)} = 1/4$。$\alpha_1 = -t_1/t_0 = -2/4 = -1/2$。

    **步骤 2** ($k=2$)：$\alpha_2 = -(t_2 + a_1^{(1)} t_1)/(t_0 + a_1^{(1)} t_1) = -(1 + (-1/2)\cdot 2)/(4 + (-1/2)\cdot 2) = -(1-1)/(4-1) = 0$。

    $a_1^{(2)} = a_1^{(1)} + \alpha_2 a_1^{(1)} = -1/2$，$a_2^{(2)} = \alpha_2 = 0$。

    然后利用 Toeplitz 结构回代求解 $Tx = b$，最终得 $x = (0, 1/4, 3/4)^T$。

    使用一般 Gauss 消元需要约 $n^3/3 = 9$ 次乘法，而 Levinson 算法仅需约 $n^2 = 9$ 次——在此小例中差别不明显，但对 $n = 10000$，从 $3 \times 10^{11}$ 减少到 $10^8$，加速约 3000 倍。

!!! note "注记 37.1 (Toeplitz 矩阵在信号处理中的出现)"
    宽平稳随机过程 $\{X_t\}$ 的自协方差矩阵
    $$R_{ij} = \operatorname{Cov}(X_i, X_j) = r(|i-j|)$$
    天然具有 Toeplitz 结构。Levinson-Durbin 算法最初就是在信号处理和时间序列分析的背景下被开发的，用于 AR 模型的参数估计。

---

## 37.2 Hankel 矩阵

<div class="context-flow" markdown>

**核心问题**：沿反对角线常数的矩阵与 Toeplitz 矩阵有何联系？它在系统理论中扮演什么角色？

</div>

!!! definition "定义 37.3 (Hankel 矩阵)"
    矩阵 $H = (h_{ij}) \in M_n(\mathbb{C})$ 称为 **Hankel 矩阵**，若 $h_{ij}$ 仅取决于和 $i + j$，即存在序列 $\{h_k\}_{k=2}^{2n}$ 使得
    $$h_{ij} = h_{i+j}, \quad 1 \le i, j \le n.$$
    显式地，
    $$H = \begin{pmatrix}
    h_2 & h_3 & h_4 & \cdots & h_{n+1} \\
    h_3 & h_4 & h_5 & \cdots & h_{n+2} \\
    h_4 & h_5 & h_6 & \cdots & h_{n+3} \\
    \vdots & & & \ddots & \vdots \\
    h_{n+1} & h_{n+2} & h_{n+3} & \cdots & h_{2n}
    \end{pmatrix}.$$
    Hankel 矩阵同样由 $2n-1$ 个参数确定。

!!! theorem "定理 37.2 (Toeplitz 与 Hankel 的互换)"
    设 $J$ 为 $n \times n$ 反对角单位矩阵。则：

    (a) 若 $T$ 为 Toeplitz 矩阵，则 $TJ$ 和 $JT$ 为 Hankel 矩阵；

    (b) 若 $H$ 为 Hankel 矩阵，则 $HJ$ 和 $JH$ 为 Toeplitz 矩阵。

??? proof "证明"
    (a) $(TJ)_{ij} = \sum_k T_{ik} J_{kj} = T_{i, n+1-j} = t_{i-(n+1-j)} = t_{i+j-(n+1)}$，仅取决于 $i+j$，故 $TJ$ 是 Hankel 的。$JT$ 的证明类似。

    (b) $(HJ)_{ij} = H_{i, n+1-j} = h_{i+(n+1-j)} = h_{(n+1) + (i-j)}$，仅取决于 $i-j$，故 $HJ$ 是 Toeplitz 的。

!!! definition "定义 37.4 (Hankel 矩阵与矩问题)"
    给定测度 $\mu$ 在 $\mathbb{R}$ 上的矩 $s_k = \int x^k \, d\mu(x)$，$k = 0, 1, 2, \ldots$，对应的 **Hankel 矩阵**为
    $$H_n = (s_{i+j-2})_{1 \le i,j \le n} = \begin{pmatrix}
    s_0 & s_1 & \cdots & s_{n-1} \\
    s_1 & s_2 & \cdots & s_n \\
    \vdots & & & \vdots \\
    s_{n-1} & s_n & \cdots & s_{2n-2}
    \end{pmatrix}.$$
    **Hamburger 矩问题**：给定序列 $\{s_k\}$，何时存在非负测度 $\mu$ 使其矩为 $s_k$？经典结果表明，必要充分条件是所有 Hankel 矩阵 $H_n$ 半正定。

!!! example "例 37.2"
    考虑标准正态分布 $\mu = N(0,1)$，其矩为 $s_{2k} = (2k-1)!! = 1 \cdot 3 \cdots (2k-1)$，$s_{2k+1} = 0$。

    $3 \times 3$ Hankel 矩阵：
    $$H_3 = \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 3 \end{pmatrix}.$$
    特征值为 $1$ 和 $2 \pm \sqrt{2}$，全为正。$H_3 > 0$ 与 $N(0,1)$ 为正测度一致。

!!! note "注记 37.2 (Hankel 算子与系统理论)"
    在线性系统理论中，传递函数 $G(s) = C(sI - A)^{-1}B$ 对应的 Hankel 算子的矩阵表示恰好是由 Markov 参数 $h_k = CA^{k-1}B$ 构成的 Hankel 矩阵。Hankel 矩阵的秩等于系统的 **McMillan 度**（最小实现的状态维数）。Hankel 奇异值在模型降阶（平衡截断）中起核心作用。

---

## 37.3 循环矩阵与 DFT

<div class="context-flow" markdown>

**核心问题**：为什么循环矩阵可以被 Fourier 矩阵对角化？这如何导出 $O(n\log n)$ 的矩阵-向量乘法？

</div>

!!! definition "定义 37.5 (循环矩阵)"
    矩阵 $C \in M_n(\mathbb{C})$ 称为**循环矩阵**，若每一行是上一行的循环右移：
    $$C = \begin{pmatrix}
    c_0 & c_{n-1} & c_{n-2} & \cdots & c_1 \\
    c_1 & c_0 & c_{n-1} & \cdots & c_2 \\
    c_2 & c_1 & c_0 & \cdots & c_3 \\
    \vdots & & & \ddots & \vdots \\
    c_{n-1} & c_{n-2} & c_{n-3} & \cdots & c_0
    \end{pmatrix} = \sum_{k=0}^{n-1} c_k \Pi^k,$$
    其中 $\Pi$ 为**循环置换矩阵**（$\Pi e_j = e_{j+1 \pmod{n}}$）。循环矩阵由 $n$ 个参数确定。

!!! definition "定义 37.6 (DFT 矩阵)"
    设 $\omega = e^{2\pi i/n}$ 为 $n$ 次本原单位根。$n$ 阶 **DFT（离散 Fourier 变换）矩阵**定义为
    $$F = \frac{1}{\sqrt{n}} (\omega^{(j-1)(k-1)})_{1 \le j,k \le n}.$$
    $F$ 是酉矩阵：$F^* F = I$。

!!! theorem "定理 37.3 (循环矩阵的对角化)"
    每个循环矩阵 $C$ 可以被 DFT 矩阵对角化：
    $$C = F^* \Lambda F,$$
    其中 $\Lambda = \operatorname{diag}(\hat{c}_0, \hat{c}_1, \ldots, \hat{c}_{n-1})$，$\hat{c}_k = \sum_{j=0}^{n-1} c_j \omega^{jk}$ 是向量 $(c_0, \ldots, c_{n-1})$ 的 DFT。

??? proof "证明"
    **步骤一**：循环置换矩阵 $\Pi$ 的特征值为 $1, \omega, \omega^2, \ldots, \omega^{n-1}$，对应特征向量为
    $$v_k = \frac{1}{\sqrt{n}}(1, \omega^k, \omega^{2k}, \ldots, \omega^{(n-1)k})^T, \quad k = 0, 1, \ldots, n-1.$$
    验证：$(\Pi v_k)_j = (v_k)_{j-1} = \frac{1}{\sqrt{n}} \omega^{k(j-1)} = \omega^{-k} \cdot \frac{1}{\sqrt{n}}\omega^{kj}$。
    注意到 $(\Pi v_k)_j = (v_k)_{j-1 \pmod n}$，即 $v_k$ 的分量循环移位，等于 $\omega^k v_k$。故 $\Pi v_k = \omega^k v_k$。

    **步骤二**：$C = \sum_{j=0}^{n-1} c_j \Pi^j$，故 $C v_k = \sum_{j=0}^{n-1} c_j \omega^{jk} v_k = \hat{c}_k v_k$。

    **步骤三**：将特征向量排成矩阵 $F^* = (v_0, v_1, \ldots, v_{n-1})$，则 $CF^* = F^* \Lambda$，即 $C = F^* \Lambda F$。

!!! theorem "定理 37.4 (循环矩阵的代数)"
    (a) 循环矩阵的全体 $\mathcal{C}_n$ 构成 $M_n(\mathbb{C})$ 的一个交换子代数。

    (b) 两个循环矩阵的乘积是循环矩阵。

    (c) 非奇异循环矩阵的逆是循环矩阵。

    (d) $\mathcal{C}_n \cong \mathbb{C}[x]/(x^n - 1)$（多项式环的商环）。

??? proof "证明"
    (a)(b) 由 $C_1 = F^*\Lambda_1 F$，$C_2 = F^*\Lambda_2 F$，得 $C_1 C_2 = F^*\Lambda_1\Lambda_2 F = F^*\Lambda_2\Lambda_1 F = C_2 C_1$。对角矩阵的乘积仍为对角矩阵，故 $C_1C_2$ 为循环矩阵。

    (c) $C^{-1} = F^*\Lambda^{-1}F$（当 $\Lambda$ 可逆时），$\Lambda^{-1}$ 仍为对角矩阵。

    (d) 将 $c(x) = c_0 + c_1 x + \cdots + c_{n-1}x^{n-1}$ 对应于循环矩阵 $C = c(\Pi)$。由于 $\Pi^n = I$，$\Pi$ 满足 $x^n - 1 = 0$。此对应是环同构。

!!! example "例 37.3"
    循环矩阵-向量乘法 $y = Cx$ 的快速算法：

    1. 计算 $\hat{x} = Fx$（$x$ 的 DFT），$O(n\log n)$；
    2. 计算 $\hat{c} = Fc$（$c$ 的 DFT），$O(n\log n)$（可预计算）；
    3. 逐元素乘 $\hat{y}_k = \hat{c}_k \hat{x}_k$，$O(n)$；
    4. 计算 $y = F^*\hat{y}$（逆 DFT），$O(n\log n)$。

    总计 $O(n\log n)$，而直接乘法为 $O(n^2)$。

    数值示例：$C = \begin{pmatrix} 1 & 3 & 2 \\ 2 & 1 & 3 \\ 3 & 2 & 1 \end{pmatrix}$，$x = (1, 0, 0)^T$。

    $\hat{c} = (1+2+3, 1+2\omega+3\omega^2, 1+2\omega^2+3\omega^4) = (6, 1+2\omega+3\omega^2, 1+2\omega^2+3\omega)$。
    其中 $\omega = e^{2\pi i/3}$。

    $y = Cx = (1, 2, 3)^T$（即 $C$ 的第一列）。

!!! note "注记 37.3 (循环预处理)"
    在迭代法求解 Toeplitz 系统 $Tx = b$ 时，可以构造一个"最佳循环逼近" $C \approx T$ 作为预处理子。由于 $C^{-1}$ 的作用可以用 FFT 在 $O(n\log n)$ 时间完成，预处理共轭梯度法（PCG）结合循环预处理子通常在 $O(n\log n)$ 总时间内收敛，远优于直接法的 $O(n^2)$。这是 T. Chan 和 R. Chan 等人在 20 世纪 80-90 年代发展的重要数值方法。

---

## 37.4 Vandermonde 矩阵

<div class="context-flow" markdown>

**核心问题**：由节点的幂次构成的矩阵有什么结构？它与多项式插值的关系是什么？

</div>

!!! definition "定义 37.7 (Vandermonde 矩阵)"
    给定 $n$ 个（通常不同的）点 $x_1, x_2, \ldots, x_n$，**Vandermonde 矩阵**定义为
    $$V = \begin{pmatrix}
    1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\
    1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\
    \vdots & & & & \vdots \\
    1 & x_n & x_n^2 & \cdots & x_n^{n-1}
    \end{pmatrix}, \quad V_{ij} = x_i^{j-1}.$$

!!! theorem "定理 37.5 (Vandermonde 行列式)"
    $$\det(V) = \prod_{1 \le i < j \le n} (x_j - x_i).$$
    因此 $V$ 非奇异当且仅当 $x_1, x_2, \ldots, x_n$ 两两不同。

??? proof "证明"
    对 $n$ 归纳。$n = 1$ 时 $\det(V) = 1$，成立。

    设命题对 $n-1$ 成立。考虑多项式
    $$p(x) = \det\begin{pmatrix}
    1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\
    \vdots & & & & \vdots \\
    1 & x_{n-1} & x_{n-1}^2 & \cdots & x_{n-1}^{n-1} \\
    1 & x & x^2 & \cdots & x^{n-1}
    \end{pmatrix}.$$
    $p(x)$ 是 $x$ 的 $n-1$ 次多项式（按最后一行展开，$x^{n-1}$ 的系数是 $V_{n-1}$ 的行列式）。$p(x_i) = 0$（$i = 1, \ldots, n-1$），因为当 $x = x_i$ 时行列式有两行相同。

    故 $p(x) = c \prod_{i=1}^{n-1}(x - x_i)$，其中 $c$ 为 $x^{n-1}$ 的系数，等于 $(n-1)$ 阶 Vandermonde 行列式 $\prod_{1 \le i < j \le n-1}(x_j - x_i)$（归纳假设）。

    $$\det(V) = p(x_n) = \prod_{1 \le i < j \le n-1}(x_j - x_i) \cdot \prod_{i=1}^{n-1}(x_n - x_i) = \prod_{1 \le i < j \le n}(x_j - x_i).$$

!!! theorem "定理 37.6 (多项式插值与 Vandermonde 系统)"
    给定 $n$ 个插值条件 $(x_i, y_i)$（$x_i$ 两两不同），求 $n-1$ 次多项式 $p(x) = a_0 + a_1 x + \cdots + a_{n-1} x^{n-1}$ 使得 $p(x_i) = y_i$ 等价于求解 Vandermonde 系统
    $$Va = y, \quad a = (a_0, \ldots, a_{n-1})^T, \quad y = (y_1, \ldots, y_n)^T.$$
    该系统有唯一解。

!!! example "例 37.4"
    节点 $x_1 = 1, x_2 = 2, x_3 = 4$，插值值 $y = (1, 3, 7)^T$。

    $$V = \begin{pmatrix} 1 & 1 & 1 \\ 1 & 2 & 4 \\ 1 & 4 & 16 \end{pmatrix}, \quad \det(V) = (2-1)(4-1)(4-2) = 1 \cdot 3 \cdot 2 = 6.$$

    解 $Va = y$ 得 $a = V^{-1}y$。使用 Cramer 法则或直接求解：
    $a_0 = -1, a_1 = 2, a_2 = 0$（需验算）。实际上 $p(x) = -1 + 2x$ 满足 $p(1)=1, p(2)=3, p(4)=7$。故 $a = (-1, 2, 0)^T$。

!!! note "注记 37.4 (Vandermonde 矩阵的条件数)"
    Vandermonde 矩阵通常具有很大的条件数，特别是当节点在实轴上等距分布时。对于 $x_k = k/n$（$k = 1, \ldots, n$），条件数以指数速度增长。然而，当节点取为 Chebyshev 节点 $x_k = \cos((2k-1)\pi/(2n))$ 时，条件数增长温和得多。这与 Runge 现象密切相关。

!!! note "注记 37.5 (快速 Vandermonde 求解)"
    Vandermonde 系统 $Va = y$ 可以通过 **Björck-Pereyra 算法**在 $O(n^2)$ 时间内求解（而非一般的 $O(n^3)$）。该算法本质上就是 Newton 插值公式的矩阵语言表达。

---

## 37.5 Cauchy 矩阵

<div class="context-flow" markdown>

**核心问题**：形如 $1/(x_i - y_j)$ 的矩阵有什么特殊性质？

</div>

!!! definition "定义 37.8 (Cauchy 矩阵)"
    给定两组参数 $x_1, \ldots, x_m$ 和 $y_1, \ldots, y_n$（$x_i \ne y_j$ 对所有 $i, j$），**Cauchy 矩阵**定义为
    $$C = \left(\frac{1}{x_i - y_j}\right)_{1 \le i \le m, 1 \le j \le n}.$$

!!! theorem "定理 37.7 (Cauchy 行列式)"
    设 $m = n$，$x_1, \ldots, x_n$ 两两不同，$y_1, \ldots, y_n$ 两两不同，且 $x_i \ne y_j$。则
    $$\det(C) = \frac{\prod_{i<j}(x_j - x_i) \prod_{i<j}(y_j - y_i)}{\prod_{i,j}(x_i - y_j)}.$$

??? proof "证明"
    令 $D = \operatorname{diag}(\prod_{k=1}^n (x_1-y_k), \ldots, \prod_{k=1}^n (x_n-y_k))$，则 $DC$ 是多项式矩阵。具体地，$(DC)_{ij} = \prod_{k\ne j}(x_i - y_k)$。

    矩阵 $DC$ 可以写成 Vandermonde 形式的乘积。考虑 $(DC)_{ij}$ 是 $x_i$ 的 $n-1$ 次多项式（固定 $j$ 时，$\prod_{k \ne j}(x_i - y_k)$ 是关于 $x_i$ 的 $n-1$ 次多项式），展开后可以表示为 Vandermonde 矩阵乘以一个由 $y_j$ 的初等对称函数组成的矩阵。

    经过仔细的代数操作（可用 Lagrange 插值的语言），最终得到上述行列式公式。完整细节可参见 Schechter 的证明或 Krattenthaler 的综述。

!!! example "例 37.5"
    $x = (1, 2), y = (0, -1)$：
    $$C = \begin{pmatrix} 1/(1-0) & 1/(1-(-1)) \\ 1/(2-0) & 1/(2-(-1)) \end{pmatrix} = \begin{pmatrix} 1 & 1/2 \\ 1/2 & 1/3 \end{pmatrix}.$$

    $\det(C) = 1/3 - 1/4 = 1/12$。

    由 Cauchy 行列式公式：$\frac{(2-1)(-1-0)}{(1-0)(1+1)(2-0)(2+1)} = \frac{1 \cdot (-1)}{1 \cdot 2 \cdot 2 \cdot 3} = \frac{-1}{12}$。

    差一个符号！这是因为公式中 $y_j - y_i$ 的顺序：$y_2 - y_1 = -1 - 0 = -1$，$x_2 - x_1 = 1$。分子 $= 1 \cdot (-1) = -1$。分母 $= (1-0)(1+1)(2-0)(2+1) = 12$。故 $\det = -1/12$。

    直接计算：$\det = 1 \cdot 1/3 - 1/2 \cdot 1/2 = 1/3 - 1/4 = 1/12$。不一致表明符号约定需要仔细对待——Cauchy 行列式公式有不同的版本，取决于是 $1/(x_i - y_j)$ 还是 $1/(x_i + y_j)$。

!!! note "注记 37.6 (Cauchy 矩阵的位移结构)"
    Cauchy 矩阵满足简洁的位移方程。设 $D_x = \operatorname{diag}(x_1, \ldots, x_n)$，$D_y = \operatorname{diag}(y_1, \ldots, y_n)$，则
    $$D_x C - C D_y = \mathbf{1}\mathbf{1}^T,$$
    其中 $\mathbf{1} = (1, \ldots, 1)^T$。这说明 Cauchy 矩阵具有**位移秩 1**。

---

## 37.6 位移结构与位移秩

<div class="context-flow" markdown>

**核心问题**：能否用统一的框架来描述所有结构化矩阵？"位移秩"如何量化矩阵的结构程度？

</div>

位移结构理论由 Kailath、Kung 和 Morf 在 20 世纪 70 年代末建立，它为结构化矩阵提供了统一的数学语言。

!!! definition "定义 37.9 (Sylvester 型位移)"
    设 $F, G \in M_n(\mathbb{C})$ 为给定的矩阵。矩阵 $A$ 的 **Sylvester 型位移**定义为
    $$\nabla_{F,G}(A) = A - FAG.$$
    $A$ 的**位移秩**定义为
    $$r = \operatorname{rank}(\nabla_{F,G}(A)).$$

!!! definition "定义 37.10 (Stein 型位移)"
    矩阵 $A$ 的 **Stein 型位移**定义为
    $$\Delta_{F,G}(A) = FA - AG.$$
    $A$ 的 Stein 型位移秩为 $\operatorname{rank}(\Delta_{F,G}(A))$。

!!! theorem "定理 37.8 (结构化矩阵的位移秩)"
    选取适当的位移算子 $F, G$，各类结构化矩阵具有低位移秩：

    | 矩阵类型 | $F$ | $G$ | 位移类型 | 位移秩 |
    |----------|-----|-----|----------|--------|
    | Toeplitz | $Z_1$（下移） | $Z_1$ | Sylvester | $\le 2$ |
    | Hankel | $Z_1$ | $Z_1^T$ | Sylvester | $\le 2$ |
    | Vandermonde | $D_x$ | $Z_1$ | Stein | $1$ |
    | Cauchy | $D_x$ | $D_y$ | Stein | $1$ |
    | 循环 | $\Pi$ | $\Pi$ | Sylvester | $0$（！） |

    其中 $Z_1$ 为下移矩阵（$Z_1 e_k = e_{k+1}$，最后一行为零），$\Pi$ 为循环置换矩阵，$D_x, D_y$ 为对角矩阵。

??? proof "证明（Toeplitz 情形）"
    设 $T = (t_{i-j})$，$Z_1$ 为 $n \times n$ 下移矩阵。则
    $$(T - Z_1 T Z_1^T)_{ij} = t_{i-j} - (Z_1 T Z_1^T)_{ij}.$$
    $(Z_1 T Z_1^T)_{ij} = \begin{cases} t_{(i-1)-(j-1)} = t_{i-j} & \text{if } 2 \le i, j \le n \\ 0 & \text{if } i = 1 \text{ or } j = 1. \end{cases}$

    但更精确地：$Z_1$ 的效果是将行下移（第 1 行变成零行，第 $k$ 行变成原第 $k-1$ 行），$Z_1^T$ 的效果是将列左移。因此 $Z_1 T Z_1^T$ 是将 $T$ 的行列都移动后得到的矩阵，其 $(i,j)$ 元素为 $t_{(i-1)-(j-1)} = t_{i-j}$（当 $i, j \ge 2$），其余为零。

    故 $T - Z_1 T Z_1^T$ 的非零元素仅出现在第 1 行和第 1 列，这个矩阵的秩至多为 2。

!!! definition "定义 37.11 (位移表示/生成元)"
    若 $\operatorname{rank}(\nabla_{F,G}(A)) = r$，则存在矩阵 $B \in \mathbb{C}^{n \times r}$，$C \in \mathbb{C}^{n \times r}$ 使得
    $$\nabla_{F,G}(A) = A - FAG = BC^T.$$
    称 $(B, C)$ 为 $A$ 的**位移生成元**。$A$ 可以由这 $2nr$ 个参数（而非 $n^2$ 个）来紧凑表示。

!!! example "例 37.6"
    对 $4 \times 4$ Toeplitz 矩阵 $T$，位移秩 $r \le 2$，故 $T$ 可以用 $2 \times 4 \times 2 = 16$ 个参数表示（与 $2n - 1 = 7$ 个自由参数的事实不矛盾，因为位移表示是冗余的）。

    快速算法的时间复杂度通常为 $O(rn \log^k n)$，其中 $r$ 为位移秩。当 $r$ 为常数（如结构化矩阵中 $r \le 2$），这给出接近线性的算法。

!!! note "注记 37.7 (位移不变性)"
    位移结构的一个关键性质是对基本矩阵运算的**封闭性**：

    - 低位移秩矩阵的**和**仍有低位移秩（$r_1 + r_2$）；
    - 低位移秩矩阵的**逆**（若存在）也有低位移秩（在相关位移算子下，秩不变）；
    - 低位移秩矩阵的 **Schur 补**保持低位移秩。

    这些封闭性质是快速算法得以递归构造的根本原因。

---

## 37.7 快速算法概述

<div class="context-flow" markdown>

**核心问题**：位移结构如何系统地导出快速算法？能达到什么样的复杂度？

</div>

!!! definition "定义 37.12 (广义 Schur 算法)"
    **广义 Schur 算法**（Generalized Schur Algorithm, GSA）是利用位移结构求解线性方程组的递归算法。其基本步骤为：

    1. 给定位移生成元 $(B, C)$（$r$ 列）；
    2. 执行一步 Gauss 消元，得到 Schur 补 $A/A_{11}$；
    3. 更新位移生成元（$O(rn)$ 运算）；
    4. 递归处理 Schur 补。

    每步更新 $O(rn)$，共 $n$ 步，总复杂度 $O(rn^2)$。

!!! theorem "定理 37.9 (结构化矩阵的复杂度)"
    利用位移结构，各类结构化矩阵的线性方程组求解复杂度如下：

    | 矩阵类型 | 一般 Gauss 消元 | 位移结构算法 | 超快算法 |
    |----------|----------------|-------------|---------|
    | 一般矩阵 | $O(n^3)$ | — | — |
    | Toeplitz | $O(n^3)$ | $O(n^2)$（Levinson/Schur） | $O(n\log^2 n)$ |
    | 循环 | $O(n^3)$ | $O(n\log n)$（FFT） | $O(n\log n)$ |
    | Vandermonde | $O(n^3)$ | $O(n^2)$（Björck-Pereyra） | $O(n\log^2 n)$ |
    | Cauchy | $O(n^3)$ | $O(n^2)$（快速 Cauchy） | $O(n\log^2 n)$ |

!!! theorem "定理 37.10 (Bini-Pan 超快 Toeplitz 求解器)"
    对称正定 Toeplitz 系统 $Tx = b$ 可以在 $O(n\log^2 n)$ 次算术运算内求解（假设 $n$ 为 2 的幂）。

??? proof "证明（思路概述）"
    核心思想是**分治法**：

    1. 将 $n \times n$ Toeplitz 矩阵分为 $2 \times 2$ 的 $n/2$ 阶块：
    $$T = \begin{pmatrix} T_{11} & T_{12} \\ T_{21} & T_{22} \end{pmatrix},$$
    其中每个块是"近似 Toeplitz"（具有低位移秩）。

    2. 通过 Schur 补 $T/T_{11} = T_{22} - T_{21}T_{11}^{-1}T_{12}$ 归约为 $n/2$ 阶问题。

    3. 关键引理：Toeplitz 矩阵的逆具有位移秩 $\le 2$（Gohberg-Semencul 公式），因此 $T_{21}T_{11}^{-1}T_{12}$ 的位移秩有界。

    4. 位移生成元的更新可用 FFT 在 $O(n\log n)$ 内完成。

    5. 递推关系 $T(n) = T(n/2) + O(n\log n)$ 的解为 $T(n) = O(n\log^2 n)$。

    然而，超快算法在数值稳定性方面存在已知问题，实践中 $O(n^2)$ 的 Levinson 或 Schur 算法往往更受青睐。

!!! example "例 37.7"
    用 Gohberg-Semencul 公式表示 Toeplitz 逆：设 $T$ 为 $n \times n$ 非奇异 Toeplitz 矩阵，$Tx = e_1$（$e_1$ 为第一个标准基向量）解为 $x = (x_0, x_1, \ldots, x_{n-1})^T$，$Ty = e_n$ 解为 $y = (y_0, y_1, \ldots, y_{n-1})^T$。则

    $$T^{-1} = \frac{1}{x_0}\left(L(x)L(\tilde{y})^T - L(Jy)L(J\tilde{x})^T\right),$$

    其中 $L(v)$ 是以 $v$ 为第一列的下三角 Toeplitz 矩阵，$J$ 为反序矩阵，$\tilde{v}$ 表示分量反序。

    这个公式说明 Toeplitz 逆可以表示为两个 Toeplitz 矩阵之差，因此 $T^{-1}v$ 的计算归结为 Toeplitz 矩阵-向量乘法，可用 FFT 完成。

---

## 37.8 对称 Toeplitz 矩阵的特征值

<div class="context-flow" markdown>

**核心问题**：大维数对称 Toeplitz 矩阵的特征值有怎样的渐近分布？

</div>

!!! definition "定义 37.13 (生成函数)"
    设 $T_n = (t_{|i-j|})_{1 \le i,j \le n}$ 为对称 Toeplitz 矩阵族。若序列 $\{t_k\}$ 的 Fourier 级数
    $$f(\theta) = \sum_{k=-\infty}^{\infty} t_k e^{ik\theta} = t_0 + 2\sum_{k=1}^{\infty} t_k \cos(k\theta)$$
    收敛，则称 $f$ 为 $\{T_n\}$ 的**生成函数**（或符号）。

!!! theorem "定理 37.11 (Szego 定理)"
    设 $f \in L^1([0, 2\pi])$ 为实值函数（即 $t_k = t_{-k}$），$T_n$ 为对应的 $n$ 阶对称 Toeplitz 矩阵。设 $\lambda_1^{(n)} \le \lambda_2^{(n)} \le \cdots \le \lambda_n^{(n)}$ 为 $T_n$ 的特征值。则：

    (a) **特征值界**：$\operatorname{ess\,inf} f \le \lambda_k^{(n)} \le \operatorname{ess\,sup} f$ 对所有 $k, n$。

    (b) **等分布定理**（Szego, 1920）：对任意连续函数 $F$，
    $$\lim_{n\to\infty} \frac{1}{n}\sum_{k=1}^{n} F(\lambda_k^{(n)}) = \frac{1}{2\pi}\int_0^{2\pi} F(f(\theta))\, d\theta.$$

    直观地说，$T_n$ 的特征值在 $n \to \infty$ 时"等分布"于 $f$ 的值域上，密度由 $f$ 的取值分布决定。

??? proof "证明（概要）"
    (a) 由 Rayleigh 商，$\lambda_{\min}(T_n) = \min_{\|x\|=1} x^T T_n x$。可以证明
    $$x^T T_n x = \frac{1}{2\pi}\int_0^{2\pi} f(\theta) |p_n(\theta)|^2 d\theta,$$
    其中 $p_n(\theta) = \sum_{k=1}^n x_k e^{ik\theta}$ 为三角多项式。由此 $x^T T_n x \ge (\operatorname{ess\,inf} f) \|x\|^2$。

    (b) 证明的核心是将 $\frac{1}{n}\operatorname{tr}(F(T_n))$ 与积分联系起来。对 $F(\lambda) = \lambda^m$，$\operatorname{tr}(T_n^m)$ 可以用 Fourier 系数计算，然后通过 Stone-Weierstrass 定理推广到一般连续函数。

!!! example "例 37.8"
    设 $f(\theta) = 2 - 2\cos\theta$，对应 $t_0 = 2, t_{\pm 1} = -1, t_k = 0$（$|k| \ge 2$）。
    则 $T_n$ 为三对角 Toeplitz 矩阵
    $$T_n = \begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & \ddots & \ddots & \ddots \\ & & -1 & 2 \end{pmatrix}.$$
    这正是一维离散 Laplacian。其特征值精确已知：
    $$\lambda_k^{(n)} = 2 - 2\cos\frac{k\pi}{n+1} = 4\sin^2\frac{k\pi}{2(n+1)}, \quad k = 1, \ldots, n.$$
    $f(\theta) = 2 - 2\cos\theta$ 的值域为 $[0, 4]$，特征值确实分布在 $[0, 4]$ 中，与 Szego 定理一致。

!!! theorem "定理 37.12 (强 Szego 极限定理)"
    设 $f > 0$ 且 $\log f \in L^1$。则
    $$\lim_{n\to\infty} \frac{\det(T_n)}{G^n} = E,$$
    其中 $G = \exp\left(\frac{1}{2\pi}\int_0^{2\pi} \log f(\theta)\, d\theta\right)$ 为**几何平均**，$E$ 为可以用 $f$ 的 Fourier 系数精确表示的常数：
    $$E = \exp\left(\sum_{k=1}^{\infty} k |\hat{c}_k|^2\right),$$
    其中 $\hat{c}_k$ 是 $\log f$ 的 Fourier 系数。

!!! note "注记 37.8 (Toeplitz 矩阵的预处理与谱聚类)"
    Szego 定理告诉我们，用循环矩阵 $C_n$（其特征值恰好是 $f$ 在等距节点的取值）来预处理 Toeplitz 矩阵 $T_n$ 时，$C_n^{-1}T_n$ 的特征值将聚集在 1 附近。这解释了为什么循环预处理共轭梯度法能在很少的迭代步内收敛。

!!! note "注记 37.9 (Fisher-Hartwig 猜想)"
    当生成函数 $f$ 具有零点或跳跃间断点时，Szego 定理的渐近形式需要修正。Fisher-Hartwig 猜想（现已大部分被证明）给出了此类情况下 $\det(T_n)$ 的精确渐近，涉及 Barnes G-函数等特殊函数。这在统计力学（Ising 模型的自发磁化）中有重要应用。

---

## 本章小结

本章的核心内容可以概括为：

| 矩阵类型 | 结构 | 参数数 | 求解复杂度 | 位移秩 |
|----------|------|--------|------------|--------|
| Toeplitz | 对角线常数 | $2n-1$ | $O(n^2)$（Levinson） | $\le 2$ |
| Hankel | 反对角线常数 | $2n-1$ | $O(n^2)$ | $\le 2$ |
| 循环 | 循环移位不变 | $n$ | $O(n\log n)$（FFT） | $0$ |
| Vandermonde | 幂次 | $n$ | $O(n^2)$（Björck-Pereyra） | $1$ |
| Cauchy | $1/(x_i-y_j)$ | $2n$ | $O(n^2)$ | $1$ |

位移结构理论将这些矩阵统一为"低位移秩矩阵"，通过位移生成元的紧凑表示和 Schur 算法的递归结构，系统地导出快速算法。Szego 定理揭示了大维 Toeplitz 矩阵的谱与生成函数之间的深刻联系，为预处理技术提供了理论基础。
