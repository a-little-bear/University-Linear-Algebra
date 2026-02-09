# 第 37B 章 Vandermonde、Cauchy 矩阵与位移结构

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 行列式(Ch3) · 矩阵分解(Ch10) · 数值线性代数(Ch22) · Toeplitz/Hankel/循环矩阵(Ch37A)

**本章脉络**：Vandermonde 矩阵与多项式插值 $\to$ Cauchy 矩阵与位移秩 1 结构 $\to$ Resultant 矩阵与 Krylov 矩阵 $\to$ Sylvester/Stein 位移算子 $\to$ 各类矩阵的位移秩计算 $\to$ 广义 Schur 算法与 Gohberg-Semencul 公式 $\to$ 层次矩阵与现代推广

**延伸**：位移结构理论（Kailath-Kung-Morf, 1979）统一了所有经典结构化矩阵类的快速算法框架；Vandermonde 矩阵的条件数分析与逼近论中的 Runge 现象紧密关联；Cauchy 矩阵在有理插值和控制理论中扮演核心角色；层次矩阵（$\mathcal{H}$-矩阵）和 HSS 矩阵将位移结构的思想推广到一般稀疏矩阵的快速直接求解

</div>

在第 37A 章中，我们系统研究了 Toeplitz、Hankel 和循环矩阵——以"沿对角线常数"或"沿反对角线常数"为特征的结构化矩阵。本章转向另外两类同样重要的结构化矩阵——Vandermonde 矩阵和 Cauchy 矩阵，它们分别以"幂次结构"和"核函数结构"为特征。随后，我们引入**位移结构**的统一框架，将 Toeplitz、Hankel、Vandermonde、Cauchy 等矩阵类纳入同一理论体系，揭示它们共同的低秩位移性质，并由此推导出系统的快速算法。最后，我们介绍层次矩阵等现代推广，它们将位移结构的核心思想——"离对角线部分的低秩性"——推广到更广泛的矩阵类。

---

## 37B.1 Vandermonde 矩阵

!!! definition "定义 37B.1 (Vandermonde 矩阵)"
    给定 $n$ 个标量 $x_1, x_2, \ldots, x_n \in \mathbb{C}$，**Vandermonde 矩阵**定义为
    $$V = V(x_1, \ldots, x_n) = \begin{pmatrix}
    1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\
    1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\
    1 & x_3 & x_3^2 & \cdots & x_3^{n-1} \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    1 & x_n & x_n^2 & \cdots & x_n^{n-1}
    \end{pmatrix}, \quad V_{ij} = x_i^{j-1}.$$
    $V$ 由 $n$ 个参数 $x_1, \ldots, x_n$（称为**节点**）完全确定。其**转置** $V^T$ 也常被称为 Vandermonde 矩阵，上下文将明确区分。

!!! theorem "定理 37B.1 (Vandermonde 行列式)"
    $$\det V(x_1, \ldots, x_n) = \prod_{1 \le i < j \le n} (x_j - x_i).$$
    因此 $V$ 非奇异当且仅当节点 $x_1, x_2, \ldots, x_n$ 两两不同。

??? proof "证明"
    对 $n$ 进行数学归纳。

    **基础情形**：$n = 1$ 时，$V = (1)$，$\det V = 1$，右端为空乘积，约定为 $1$，命题成立。$n = 2$ 时，$V = \begin{pmatrix} 1 & x_1 \\ 1 & x_2 \end{pmatrix}$，$\det V = x_2 - x_1$，与公式一致。

    **归纳步骤**：设命题对 $n-1$ 阶 Vandermonde 矩阵成立。考虑 $n$ 阶矩阵 $V = V(x_1, \ldots, x_n)$，定义关于变量 $x$ 的多项式
    $$p(x) = \det\begin{pmatrix}
    1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\
    1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    1 & x_{n-1} & x_{n-1}^2 & \cdots & x_{n-1}^{n-1} \\
    1 & x & x^2 & \cdots & x^{n-1}
    \end{pmatrix}.$$

    将行列式按最后一行展开，$x^{n-1}$ 的系数是其余各行前 $n-1$ 列构成的 $(n-1)$ 阶子行列式的代数余子式。该子行列式恰为 $(n-1)$ 阶 Vandermonde 行列式 $\det V(x_1, \ldots, x_{n-1})$。因此 $p(x)$ 是关于 $x$ 的至多 $n-1$ 次多项式，且首项系数（$x^{n-1}$ 的系数）为
    $$c = \det V(x_1, \ldots, x_{n-1}).$$

    当 $x = x_i$（$i = 1, \ldots, n-1$）时，行列式有两行相同，故 $p(x_i) = 0$。这表明 $x_1, x_2, \ldots, x_{n-1}$ 是 $p(x)$ 的 $n-1$ 个根。由于 $p(x)$ 次数恰为 $n-1$，它可以完全因式分解为
    $$p(x) = c \prod_{i=1}^{n-1}(x - x_i) = \det V(x_1, \ldots, x_{n-1}) \cdot \prod_{i=1}^{n-1}(x - x_i).$$

    令 $x = x_n$，得到
    $$\det V(x_1, \ldots, x_n) = p(x_n) = \det V(x_1, \ldots, x_{n-1}) \cdot \prod_{i=1}^{n-1}(x_n - x_i).$$

    由归纳假设，$\det V(x_1, \ldots, x_{n-1}) = \prod_{1 \le i < j \le n-1}(x_j - x_i)$，代入得
    $$\det V(x_1, \ldots, x_n) = \prod_{1 \le i < j \le n-1}(x_j - x_i) \cdot \prod_{i=1}^{n-1}(x_n - x_i) = \prod_{1 \le i < j \le n}(x_j - x_i).$$

    最后一步等号成立是因为：右端将所有满足 $1 \le i < j \le n$ 的因子 $(x_j - x_i)$ 分为两组——$j \le n-1$ 的因子（即 $\prod_{1 \le i < j \le n-1}(x_j-x_i)$）和 $j = n$ 的因子（即 $\prod_{i=1}^{n-1}(x_n - x_i)$）。归纳完成。$\blacksquare$

!!! theorem "定理 37B.2 (多项式插值的存在唯一性)"
    给定 $n$ 个互不相同的节点 $x_1, \ldots, x_n \in \mathbb{C}$ 和 $n$ 个值 $y_1, \ldots, y_n \in \mathbb{C}$，存在唯一的次数至多为 $n-1$ 的多项式
    $$p(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_{n-1} x^{n-1}$$
    满足插值条件 $p(x_i) = y_i$（$i = 1, \ldots, n$）。系数向量 $\boldsymbol{a} = (a_0, a_1, \ldots, a_{n-1})^T$ 是 Vandermonde 系统 $V\boldsymbol{a} = \boldsymbol{y}$ 的唯一解。

??? proof "证明"
    插值条件 $p(x_i) = y_i$ 即 $\sum_{j=0}^{n-1} a_j x_i^j = y_i$（$i = 1, \ldots, n$），以矩阵形式写即 $V\boldsymbol{a} = \boldsymbol{y}$，其中 $V = V(x_1, \ldots, x_n)$。由定理 37B.1，当 $x_i$ 两两不同时 $\det V \ne 0$，故方程组有唯一解。$\blacksquare$

!!! example "例 37B.1"
    取节点 $x_1 = 1, x_2 = 2, x_3 = 4$，插值值 $y_1 = 1, y_2 = 3, y_3 = 7$。Vandermonde 矩阵为
    $$V = \begin{pmatrix} 1 & 1 & 1 \\ 1 & 2 & 4 \\ 1 & 4 & 16 \end{pmatrix}.$$
    行列式 $\det V = (2-1)(4-1)(4-2) = 1 \cdot 3 \cdot 2 = 6$。解 $V\boldsymbol{a} = \boldsymbol{y}$ 得 $\boldsymbol{a} = (-1, 2, 0)^T$，故插值多项式为 $p(x) = -1 + 2x$。验证：$p(1) = 1, p(2) = 3, p(4) = 7$。

!!! theorem "定理 37B.3 (Vandermonde 矩阵的条件数)"
    设 $V = V(x_1, \ldots, x_n)$，节点 $x_i$ 为实数且两两不同。则：

    (a) 当节点为等距节点 $x_k = (k-1)/(n-1)$（$k = 1, \ldots, n$）时，$V$ 的 2-范数条件数 $\kappa_2(V)$ 随 $n$ **指数增长**。

    (b) 当节点为 Chebyshev 节点 $x_k = \cos\left(\frac{2k-1}{2n}\pi\right)$（$k = 1, \ldots, n$）时，$\kappa_2(V)$ 增长温和得多（多项式级别）。

    (c) 当节点取在单位圆上 $x_k = e^{2\pi i k/n}$（$k = 1, \ldots, n$）时，$V = \sqrt{n} \cdot F^*$（其中 $F$ 为 DFT 矩阵），$\kappa_2(V) = 1$。

??? proof "证明"
    (c) 当 $x_k = \omega^{k-1}$（$\omega = e^{2\pi i/n}$）时，$V_{kj} = \omega^{(k-1)(j-1)}$，故 $V = \sqrt{n} \cdot F^*$。由 $F$ 为酉矩阵，$V^*V = nI$，故 $\kappa_2(V) = 1$。

    (a) 对等距节点，考虑 $V^T V$ 的元素 $(V^T V)_{jk} = \sum_{i=1}^n x_i^{j-1}x_i^{k-1}$。当 $n$ 增大时，高次项的贡献使得 $V^T V$ 趋于病态。更精确地，Gautschi（1962, 1990）证明了对 $[0,1]$ 上的等距节点，
    $$\kappa_\infty(V) \ge \frac{2^{n-1}}{n \cdot e},$$
    因此条件数至少以指数速度增长。

    (b) Chebyshev 节点是 Chebyshev 多项式在 $[-1,1]$ 上的零点。由 Chebyshev 多项式的极值性质，相应的 Vandermonde 矩阵（更确切地说，在 Chebyshev 基下的变换矩阵）条件数增长缓慢。这解释了为什么数值插值中应优先选择 Chebyshev 节点。$\blacksquare$

!!! note "注记 37B.1 (Bjorck-Pereyra 算法)"
    直接求解 Vandermonde 系统 $V\boldsymbol{a} = \boldsymbol{y}$ 使用 Gauss 消元需要 $O(n^3)$ 运算。**Bjorck-Pereyra 算法**（1970）利用 Vandermonde 结构，将复杂度降至 $O(n^2)$。其核心思想等价于 Newton 插值差商的递推计算。具体地，该算法分两个阶段：

    1. **前向消元**：等价于计算 Newton 差商 $[y_{i_1}, y_{i_2}, \ldots, y_{i_k}]$，复杂度 $O(n^2)$。
    2. **回代**：从 Newton 基系数转换为幂基系数，复杂度 $O(n^2)$。

    虽然 Bjorck-Pereyra 算法在等距节点下可能数值不稳定（这与 Vandermonde 矩阵本身的病态性有关），但在 Chebyshev 节点等良好分布的节点下，它具有良好的数值行为。

---

## 37B.2 Cauchy 矩阵

!!! definition "定义 37B.2 (Cauchy 矩阵)"
    给定两组参数 $\boldsymbol{x} = (x_1, \ldots, x_n)$ 和 $\boldsymbol{y} = (y_1, \ldots, y_n)$，满足 $x_i + y_j \ne 0$ 对所有 $i, j$，**Cauchy 矩阵**定义为
    $$C = C(\boldsymbol{x}; \boldsymbol{y}) = \left(\frac{1}{x_i + y_j}\right)_{1 \le i, j \le n}.$$

    更一般地，有时也使用 $1/(x_i - y_j)$ 的约定。两种形式通过将 $y_j$ 替换为 $-y_j$ 互相转化。本章在需要时会明确指定所用约定。

!!! theorem "定理 37B.4 (Cauchy 行列式)"
    设 $x_1, \ldots, x_n$ 两两不同，$y_1, \ldots, y_n$ 两两不同，且 $x_i + y_j \ne 0$。则
    $$\det C(\boldsymbol{x}; \boldsymbol{y}) = \frac{\displaystyle\prod_{1 \le i < j \le n}(x_j - x_i) \cdot \prod_{1 \le i < j \le n}(y_j - y_i)}{\displaystyle\prod_{i=1}^n \prod_{j=1}^n (x_i + y_j)}.$$

??? proof "证明"
    **第一步**：提取公因子。将 $C$ 的第 $i$ 行乘以 $\prod_{j=1}^n(x_i + y_j)$，得到矩阵 $\tilde{C}$，其中
    $$\tilde{C}_{ij} = \prod_{k \ne j}(x_i + y_k).$$
    此时
    $$\det C = \frac{\det \tilde{C}}{\prod_{i=1}^n \prod_{j=1}^n(x_i + y_j)}.$$

    **第二步**：分析 $\tilde{C}_{ij}$。固定 $j$，$\tilde{C}_{ij} = \prod_{k \ne j}(x_i + y_k)$ 是 $x_i$ 的 $n-1$ 次多项式。将其展开：
    $$\prod_{k \ne j}(x_i + y_k) = \sum_{m=0}^{n-1} e_{n-1-m}(y_1, \ldots, \hat{y}_j, \ldots, y_n) \cdot x_i^m,$$
    其中 $e_r(y_1, \ldots, \hat{y}_j, \ldots, y_n)$ 是去掉 $y_j$ 后剩余变量的第 $r$ 个初等对称多项式。

    **第三步**：矩阵分解。上式表明 $\tilde{C} = V \cdot E$，其中 $V = V(x_1, \ldots, x_n)$ 是节点为 $x_i$ 的 Vandermonde 矩阵，$E$ 是由初等对称多项式构成的 $n \times n$ 矩阵：$E_{mj} = e_{n-1-m}(y_1, \ldots, \hat{y}_j, \ldots, y_n)$（$m = 0, 1, \ldots, n-1$）。

    **第四步**：计算 $\det E$。$E$ 可以通过 Lagrange 插值的核函数联系到 $y_j$ 的 Vandermonde 矩阵。具体地，$E$ 的行列式可以表示为
    $$\det E = \prod_{1 \le i < j \le n}(y_j - y_i).$$
    这可以通过对矩阵 $E$ 进行列运算来验证：做列变换 $E \to E \cdot V_y^{-T}$（其中 $V_y$ 是 $y_j$ 的 Vandermonde 矩阵的某种形式），化为单位矩阵的常数倍。

    **第五步**：综合以上各步：
    $$\det C = \frac{\det V \cdot \det E}{\prod_{i,j}(x_i + y_j)} = \frac{\prod_{i<j}(x_j - x_i) \cdot \prod_{i<j}(y_j - y_i)}{\prod_{i,j}(x_i + y_j)}. \quad \blacksquare$$

!!! example "例 37B.2"
    取 $x_1 = 1, x_2 = 3$，$y_1 = 2, y_2 = 4$。Cauchy 矩阵为
    $$C = \begin{pmatrix} 1/(1+2) & 1/(1+4) \\ 1/(3+2) & 1/(3+4) \end{pmatrix} = \begin{pmatrix} 1/3 & 1/5 \\ 1/5 & 1/7 \end{pmatrix}.$$
    直接计算：$\det C = \frac{1}{3} \cdot \frac{1}{7} - \frac{1}{5} \cdot \frac{1}{5} = \frac{1}{21} - \frac{1}{25} = \frac{25 - 21}{525} = \frac{4}{525}$。

    由公式：分子 $= (3-1)(4-2) = 2 \cdot 2 = 4$，分母 $= (1+2)(1+4)(3+2)(3+4) = 3 \cdot 5 \cdot 5 \cdot 7 = 525$。故 $\det C = 4/525$，与直接计算一致。

!!! definition "定义 37B.3 (Cauchy-like 矩阵)"
    矩阵 $R \in M_n(\mathbb{C})$ 称为 **Cauchy-like 矩阵**，若存在对角矩阵 $D_x = \operatorname{diag}(x_1, \ldots, x_n)$，$D_y = \operatorname{diag}(y_1, \ldots, y_n)$ 和列向量 $\boldsymbol{g}, \boldsymbol{h} \in \mathbb{C}^n$，使得
    $$D_x R - R D_y = \boldsymbol{g}\boldsymbol{h}^T,$$
    即
    $$R_{ij} = \frac{g_i h_j}{x_i - y_j}.$$
    标准 Cauchy 矩阵（约定 $1/(x_i - y_j)$）对应于 $g_i = h_j = 1$ 的特殊情形。

    更一般地，若 $D_x R - R D_y = GH^T$，其中 $G, H \in \mathbb{C}^{n \times r}$，则 $R$ 的**位移秩**为 $r$（或更小），我们称 $R$ 为秩-$r$ 的 Cauchy-like 矩阵。

!!! definition "定义 37B.4 (合流 Cauchy 矩阵)"
    **合流 Cauchy 矩阵**（confluent Cauchy matrix）是 Cauchy 矩阵在节点合并时的极限形式。当 $x_i \to x_j$ 时，$1/(x_i - y_k)$ 的差商趋于导数。具体地，若某些 $x_i$ 重合，对应的行被替换为核函数 $1/(x - y_j)$ 在重合点处的高阶导数值：
    $$\frac{1}{m!} \frac{\partial^m}{\partial x^m}\left(\frac{1}{x - y_j}\right)\bigg|_{x = x_0} = \frac{(-1)^m}{(x_0 - y_j)^{m+1}}.$$

    合流 Cauchy 矩阵在 Hermite 插值（带导数条件的插值）和多点 Pade 逼近中自然出现。

---

## 37B.3 其他结构矩阵

### Resultant 矩阵

!!! definition "定义 37B.5 (Sylvester 结式矩阵)"
    给定两个多项式
    $$f(x) = a_m x^m + a_{m-1}x^{m-1} + \cdots + a_0, \quad g(x) = b_n x^n + b_{n-1}x^{n-1} + \cdots + b_0,$$
    其 **Sylvester 结式矩阵**（resultant matrix）为 $(m+n) \times (m+n)$ 矩阵
    $$\operatorname{Syl}(f, g) = \begin{pmatrix}
    a_m & a_{m-1} & \cdots & a_0 & & & \\
    & a_m & a_{m-1} & \cdots & a_0 & & \\
    & & \ddots & & & \ddots & \\
    & & & a_m & a_{m-1} & \cdots & a_0 \\
    b_n & b_{n-1} & \cdots & b_0 & & & \\
    & b_n & b_{n-1} & \cdots & b_0 & & \\
    & & \ddots & & & \ddots & \\
    & & & b_n & b_{n-1} & \cdots & b_0
    \end{pmatrix},$$
    其中上半部分包含 $f$ 的系数的 $n$ 个移位行，下半部分包含 $g$ 的系数的 $m$ 个移位行。**结式**（resultant）定义为 $\operatorname{Res}(f, g) = \det \operatorname{Syl}(f, g)$。

!!! theorem "定理 37B.5 (结式与公共根)"
    设 $a_m \ne 0$，$b_n \ne 0$。则

    (a) $\operatorname{Res}(f, g) = 0$ 当且仅当 $f$ 与 $g$ 有公共根（在 $\mathbb{C}$ 中），即 $\gcd(f, g) \ne 1$。

    (b) 若 $f(x) = a_m \prod_{i=1}^m (x - \alpha_i)$，$g(x) = b_n \prod_{j=1}^n (x - \beta_j)$，则
    $$\operatorname{Res}(f, g) = a_m^n b_n^m \prod_{i=1}^m \prod_{j=1}^n (\alpha_i - \beta_j) = a_m^n \prod_{i=1}^m g(\alpha_i) = (-1)^{mn} b_n^m \prod_{j=1}^n f(\beta_j).$$

    (c) $\gcd(f, g)$ 的次数等于 $m + n - \operatorname{rank}(\operatorname{Syl}(f, g))$。

??? proof "证明"
    (a) $\operatorname{Res}(f, g) = 0$ 当且仅当 $\operatorname{Syl}(f, g)$ 奇异，当且仅当存在非零向量 $(\boldsymbol{u}, \boldsymbol{v})$（$\boldsymbol{u} \in \mathbb{C}^n$，$\boldsymbol{v} \in \mathbb{C}^m$）使得 $\operatorname{Syl}(f,g)^T (\boldsymbol{u}, \boldsymbol{v})^T = 0$。这等价于存在非零多项式 $u(x) = \sum u_k x^k$（$\deg u < n$）和 $v(x) = \sum v_k x^k$（$\deg v < m$）使得 $u(x)f(x) + v(x)g(x) = 0$，即 $u(x)f(x) = -v(x)g(x)$。当 $\gcd(f,g) = 1$ 时，由唯一因式分解得 $g | u$，但 $\deg u < n = \deg g$，矛盾（除非 $u = 0$，进而 $v = 0$）。反之，若 $\gcd(f,g) = d$ 且 $\deg d \ge 1$，令 $f = d\tilde{f}$，$g = d\tilde{g}$，取 $u = \tilde{g}$，$v = -\tilde{f}$，则 $uf + vg = 0$，$\deg u = n - \deg d < n$，$\deg v = m - \deg d < m$。

    (b) 利用 Vandermonde 矩阵和行列式的乘积性质可以推导。将结式矩阵的行列式用 $f$ 的根 $\alpha_i$ 表达：通过行变换和因式分解，最终得到乘积公式。

    (c) 由 (a) 的论证推广：$\operatorname{Syl}(f,g)$ 的零空间维数等于 $\deg(\gcd(f,g))$，故秩等于 $(m+n) - \deg(\gcd(f,g))$。$\blacksquare$

### Bezout 矩阵

!!! definition "定义 37B.6 (Bezout 矩阵)"
    给定多项式 $f(x)$ 和 $g(x)$（$\deg f, \deg g \le n$），**Bezout 矩阵** $B(f,g)$ 是 $n \times n$ 对称矩阵，其元素 $b_{ij}$（$0 \le i, j \le n-1$）由双线性形式
    $$\frac{f(x)g(y) - f(y)g(x)}{x - y} = \sum_{i,j=0}^{n-1} b_{ij} x^i y^j$$
    定义。左端的分子在 $x = y$ 时为零，故被 $(x-y)$ 整除，商为关于 $(x, y)$ 的次数至多为 $(n-1, n-1)$ 的多项式。

!!! theorem "定理 37B.6 (Bezout 矩阵的性质)"
    设 $f, g$ 的次数至多为 $n$，且 $\deg f = n$（首项系数不为零）。则

    (a) $B(f,g)$ 是对称矩阵。

    (b) $\operatorname{rank} B(f,g) = n - \deg(\gcd(f,g))$。

    (c) $\det B(f,g) = (-1)^{n(n-1)/2} a_n^{-(2n-\deg g - 1)} \cdot \operatorname{Res}(f, g)$（其中 $a_n$ 为 $f$ 的首项系数），即 Bezout 行列式与结式之间有确定的关系。

    (d) 当 $\gcd(f,g) = 1$ 且 $f$ 的根全部实且互不相同时，$B(f,g)$ 的惯性（正、负特征值的个数）等于 $g/f$ 在 $f$ 的各根处取值的符号变化结构。

### Krylov 矩阵

!!! definition "定义 37B.7 (Krylov 矩阵)"
    给定矩阵 $A \in M_n(\mathbb{C})$ 和向量 $\boldsymbol{b} \in \mathbb{C}^n$，**Krylov 矩阵**定义为
    $$K(A, \boldsymbol{b}) = \begin{pmatrix} \boldsymbol{b} & A\boldsymbol{b} & A^2\boldsymbol{b} & \cdots & A^{n-1}\boldsymbol{b} \end{pmatrix} \in \mathbb{C}^{n \times n}.$$
    由 $A$ 和 $\boldsymbol{b}$ 生成的 **Krylov 子空间**为
    $$\mathcal{K}_k(A, \boldsymbol{b}) = \operatorname{span}\{\boldsymbol{b}, A\boldsymbol{b}, A^2\boldsymbol{b}, \ldots, A^{k-1}\boldsymbol{b}\}, \quad k = 1, 2, \ldots$$

!!! theorem "定理 37B.7 (Krylov 矩阵与可控性)"
    (a) $K(A, \boldsymbol{b})$ 非奇异当且仅当 $(\boldsymbol{b}, A\boldsymbol{b}, \ldots, A^{n-1}\boldsymbol{b})$ 线性无关，当且仅当 $A$ 的最小多项式等于其特征多项式（关于 $\boldsymbol{b}$ 的循环向量条件）。

    (b) 在控制理论中，线性系统 $\dot{\boldsymbol{x}} = A\boldsymbol{x} + \boldsymbol{b}u$ **可控**（controllable）当且仅当 $\operatorname{rank} K(A, \boldsymbol{b}) = n$（Kalman 可控性判据）。

    (c) 当 $A = D_x = \operatorname{diag}(x_1, \ldots, x_n)$（对角矩阵，$x_i$ 两两不同）时，$K(D_x, \boldsymbol{b}) = \operatorname{diag}(b_1, \ldots, b_n) \cdot V(x_1, \ldots, x_n)$，即 Krylov 矩阵化为 Vandermonde 矩阵的行缩放形式。

??? proof "证明"
    (a) $K(A,\boldsymbol{b})$ 的列空间为 $\mathcal{K}_n(A,\boldsymbol{b})$。设 $p(x)$ 为 $A$ 关于 $\boldsymbol{b}$ 的最小多项式，即满足 $p(A)\boldsymbol{b} = 0$ 的最低次首一多项式。则 $\mathcal{K}_n(A,\boldsymbol{b})$ 的维数为 $\deg p$。$K$ 非奇异当且仅当此维数为 $n$，即 $\deg p = n$。由 Cayley-Hamilton 定理，$p$ 整除 $A$ 的特征多项式 $\chi_A$（$\deg \chi_A = n$），故 $\deg p = n$ 当且仅当 $p = \chi_A$。

    (b) 这是 Kalman 可控性定理。$(A, \boldsymbol{b})$ 可控意味着通过控制输入 $u(t)$，系统可以从任意初始状态到达任意目标状态。可控性的代数条件恰为可控性矩阵 $K(A,\boldsymbol{b})$ 满秩。

    (c) 当 $A = D_x$ 时，$(D_x^k \boldsymbol{b})_i = x_i^k b_i$，故 $K(D_x, \boldsymbol{b})_{ij} = x_i^{j-1} b_i = b_i \cdot V_{ij}$。$\blacksquare$

!!! note "注记 37B.2 (Krylov 子空间方法)"
    Krylov 矩阵和 Krylov 子空间是现代迭代法——如共轭梯度法（CG）、GMRES、Lanczos 算法——的理论基础。这些方法在第 $k$ 步从 $\mathcal{K}_k(A, \boldsymbol{b})$ 中寻找线性方程组 $A\boldsymbol{x} = \boldsymbol{b}$ 的近似解，无需显式构造 Krylov 矩阵（这会引入数值不稳定性），而是通过正交化过程（如 Arnoldi 或 Lanczos 过程）隐式地工作。

---

## 37B.4 位移结构理论

位移结构理论由 Kailath、Kung 和 Morf（1979）建立，为结构化矩阵提供了统一的数学语言。其核心思想是：结构化矩阵虽然本身不是低秩的，但在某种"位移算子"的作用下变成低秩矩阵。

!!! definition "定义 37B.8 (Sylvester 型位移)"
    设 $F, F' \in M_n(\mathbb{C})$ 为给定的**位移算子**。矩阵 $A \in M_n(\mathbb{C})$ 的 **Sylvester 型位移**（也称 Stein 型位移）定义为
    $$\nabla_{F,F'}(A) = FA - AF'.$$
    $A$ 的 **Sylvester 型位移秩**定义为
    $$\operatorname{d-rank}_S(A) = \operatorname{rank}(\nabla_{F,F'}(A)) = \operatorname{rank}(FA - AF').$$

!!! definition "定义 37B.9 (Stein 型位移)"
    矩阵 $A$ 的 **Stein 型位移**定义为
    $$\Delta_{F,F'}(A) = A - FAF'.$$
    $A$ 的 **Stein 型位移秩**为
    $$\operatorname{d-rank}_T(A) = \operatorname{rank}(\Delta_{F,F'}(A)) = \operatorname{rank}(A - FAF').$$

!!! theorem "定理 37B.8 (Sylvester 与 Stein 型位移的关系)"
    当 $F$ 可逆时，
    $$\Delta_{F,F'}(A) = A - FAF' = F(F^{-1}A - AF') = F \cdot \nabla_{F^{-1}, F'}(A),$$
    因此 $\operatorname{rank}(\Delta_{F,F'}(A)) = \operatorname{rank}(\nabla_{F^{-1},F'}(A))$。类似地，当 $F'$ 可逆时，
    $$\nabla_{F,F'}(A) = FA - AF' = (A - F^{-1}A F') \cdot F' = \Delta_{F^{-1}, F'}(A) \cdot F',$$
    两种位移秩相同。特别地，当位移算子可逆时，Sylvester 型和 Stein 型位移秩一致。

??? proof "证明"
    第一个等式：$A - FAF' = A - FAF' + FA - FA = F(F^{-1}A - AF') + (A - FA) \cdot 0$。更直接地：
    $$A - FAF' = F \cdot F^{-1} \cdot A - F \cdot A \cdot F' = F(F^{-1}A - AF').$$
    由于乘以可逆矩阵 $F$ 不改变秩，$\operatorname{rank}(A - FAF') = \operatorname{rank}(F^{-1}A - AF')$。第二个等式类似证明。$\blacksquare$

!!! definition "定义 37B.10 (位移生成元)"
    若 $\operatorname{rank}(\nabla_{F,F'}(A)) = r$，则存在矩阵 $G \in \mathbb{C}^{n \times r}$ 和 $H \in \mathbb{C}^{n \times r}$ 使得
    $$\nabla_{F,F'}(A) = FA - AF' = GH^T.$$
    称有序对 $(G, H)$ 为 $A$ 关于位移算子 $(F, F')$ 的 **Sylvester 型位移生成元**。

    类似地，若 $\operatorname{rank}(\Delta_{F,F'}(A)) = r$，则 $A - FAF' = GH^T$，$(G, H)$ 为 **Stein 型位移生成元**。

    位移生成元用 $2nr$ 个参数紧凑表示了原本需要 $n^2$ 个参数的矩阵 $A$。当 $r \ll n$ 时，这是巨大的压缩。

!!! theorem "定理 37B.9 (位移结构在矩阵运算下的封闭性)"
    设 $A$ 具有位移秩 $r$（关于算子 $(F, F')$），即 $FA - AF' = GH^T$。则：

    (a) **加法封闭**：若 $B$ 的位移秩为 $s$，则 $A + B$ 的位移秩至多为 $r + s$。

    (b) **逆的位移秩保持**：若 $A$ 可逆，则 $A^{-1}$ 的位移秩（关于算子 $(F', F)$，注意顺序交换）也为 $r$。具体地，
    $$F'A^{-1} - A^{-1}F = -A^{-1}(FA - AF')A^{-1} = -A^{-1}GH^T A^{-1} = \tilde{G}\tilde{H}^T,$$
    其中 $\tilde{G} = -A^{-1}G$，$\tilde{H} = A^{-T}H$。

    (c) **Schur 补的位移秩保持**：设 $A$ 分块为 $A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}$，且 $A$ 具有位移秩 $r$，$A_{11}$ 可逆。则 Schur 补 $A/A_{11} = A_{22} - A_{21}A_{11}^{-1}A_{12}$ 的位移秩也至多为 $r$。

??? proof "证明"
    (a) $F(A+B) - (A+B)F' = (FA - AF') + (FB - BF') = GH^T + G_B H_B^T$。右端秩至多为 $r + s$。

    (b) 由 $FA - AF' = GH^T$，左乘 $A^{-1}$，右乘 $A^{-1}$：
    $$A^{-1}FA \cdot A^{-1} - A^{-1} \cdot AF' A^{-1} = A^{-1}GH^T A^{-1},$$
    即 $A^{-1}F - F'A^{-1} = A^{-1}GH^TA^{-1}$。移项得 $F'A^{-1} - A^{-1}F = -(A^{-1}G)(A^{-T}H)^T$，秩为 $r$。

    (c) 这是 (a) 和 (b) 的推论。Schur 补可以写为加法和逆运算的组合：$A_{22} - A_{21}A_{11}^{-1}A_{12}$。利用位移结构在子矩阵提取和乘法下的性质（需要更细致的论证），可以证明位移秩不增。详细证明参见 Kailath-Sayed (1999)。$\blacksquare$

---

## 37B.5 各类矩阵的位移秩

本节计算各类经典结构化矩阵的位移秩，展示位移结构理论如何将这些看似不同的矩阵类统一在同一框架下。

!!! definition "定义 37B.11 (常用位移算子)"
    以下矩阵常被用作位移算子：

    - **下移矩阵** $Z_1$：$(Z_1)_{ij} = \delta_{i,j+1}$，即 $Z_1 e_k = e_{k+1}$（$k < n$），$Z_1 e_n = 0$。
    - **循环置换矩阵** $Z_f$：$(Z_f)_{ij} = \delta_{i,j+1 \bmod n}$，即 $Z_f e_k = e_{k+1}$（$k < n$），$Z_f e_n = f \cdot e_1$（当 $f = 1$ 时即标准循环矩阵 $\Pi$）。
    - **对角矩阵** $D_x = \operatorname{diag}(x_1, \ldots, x_n)$。

!!! theorem "定理 37B.10 (Toeplitz 矩阵的位移秩)"
    设 $T = (t_{i-j})$ 为 $n \times n$ Toeplitz 矩阵。以 $Z_1$（下移矩阵）为位移算子，Stein 型位移
    $$\Delta(T) = T - Z_1 T Z_1^T$$
    的秩至多为 $2$。因此 Toeplitz 矩阵的 Stein 型位移秩至多为 $2$。

??? proof "证明"
    计算 $\Delta(T) = T - Z_1 T Z_1^T$ 的各元素。注意 $Z_1$ 将行下移一位（第一行变为零行），$Z_1^T$ 将列左移一位（第一列变为零列）。故

    $$(Z_1 T Z_1^T)_{ij} = \begin{cases}
    T_{i-1, j-1} = t_{(i-1)-(j-1)} = t_{i-j} & \text{若 } i \ge 2, j \ge 2, \\
    0 & \text{若 } i = 1 \text{ 或 } j = 1.
    \end{cases}$$

    因此 $\Delta(T)_{ij} = t_{i-j} - t_{i-j} = 0$（当 $i, j \ge 2$），而第一行和第一列为
    $$\Delta(T)_{1j} = t_{1-j} \quad (j = 1, \ldots, n), \qquad \Delta(T)_{i1} = t_{i-1} \quad (i = 1, \ldots, n).$$

    由此 $\Delta(T)$ 仅在第一行和第一列有非零元素，其形式为
    $$\Delta(T) = \boldsymbol{e}_1 \boldsymbol{r}^T + \boldsymbol{c} \boldsymbol{e}_1^T - t_0 \boldsymbol{e}_1 \boldsymbol{e}_1^T,$$
    其中 $\boldsymbol{r} = (t_0, t_{-1}, \ldots, t_{-(n-1)})^T$（第一行），$\boldsymbol{c} = (t_0, t_1, \ldots, t_{n-1})^T$（第一列）。这是至多两个秩-1 矩阵之和（减去重复计数的 $(1,1)$ 位置），故 $\operatorname{rank}(\Delta(T)) \le 2$。$\blacksquare$

!!! theorem "定理 37B.11 (Hankel 矩阵的位移秩)"
    设 $H = (h_{i+j})$ 为 $n \times n$ Hankel 矩阵。Stein 型位移
    $$\Delta(H) = H - Z_1 H Z_1$$
    的秩至多为 $2$。

??? proof "证明"
    注意这里用的是 $Z_1 H Z_1$（而非 $Z_1 H Z_1^T$），这是因为 Hankel 矩阵沿反对角线常数而非沿对角线常数。计算：

    $$(Z_1 H Z_1)_{ij} = \begin{cases}
    H_{i-1, j+1} = h_{(i-1)+(j+1)} = h_{i+j} & \text{若 } i \ge 2, j \le n-1, \\
    0 & \text{若 } i = 1 \text{ 或 } j = n.
    \end{cases}$$

    故 $\Delta(H)_{ij} = h_{i+j} - h_{i+j} = 0$（当 $i \ge 2$ 且 $j \le n-1$），$\Delta(H)$ 仅在第一行和最后一列有非零元素。因此 $\operatorname{rank}(\Delta(H)) \le 2$。$\blacksquare$

!!! theorem "定理 37B.12 (Cauchy 矩阵的位移秩)"
    设 $C = (1/(x_i - y_j))$ 为 Cauchy 矩阵，$D_x = \operatorname{diag}(x_1, \ldots, x_n)$，$D_y = \operatorname{diag}(y_1, \ldots, y_n)$。则 Sylvester 型位移
    $$\nabla(C) = D_x C - C D_y = \boldsymbol{1}\boldsymbol{1}^T$$
    的秩为 $1$（其中 $\boldsymbol{1} = (1, \ldots, 1)^T$）。因此 Cauchy 矩阵的位移秩为 $1$。

??? proof "证明"
    直接验证。$\nabla(C)$ 的 $(i,j)$ 元素为
    $$(D_x C - C D_y)_{ij} = x_i \cdot \frac{1}{x_i - y_j} - \frac{1}{x_i - y_j} \cdot y_j = \frac{x_i - y_j}{x_i - y_j} = 1.$$
    因此 $D_x C - C D_y = \boldsymbol{1}\boldsymbol{1}^T$，秩为 $1$。$\blacksquare$

!!! theorem "定理 37B.13 (Vandermonde 矩阵的位移秩)"
    设 $V = V(x_1, \ldots, x_n)$ 为 Vandermonde 矩阵，$D_x = \operatorname{diag}(x_1, \ldots, x_n)$。则 Sylvester 型位移
    $$\nabla(V) = D_x V - V Z_1^T$$
    的秩为 $1$。因此 Vandermonde 矩阵的位移秩为 $1$。

??? proof "证明"
    $V$ 的 $(i,j)$ 元素为 $V_{ij} = x_i^{j-1}$。计算：

    $$(D_x V)_{ij} = x_i \cdot x_i^{j-1} = x_i^j.$$

    $(VZ_1^T)$ 的效果是将 $V$ 的列左移一位（最后一列变为零列）：
    $$(VZ_1^T)_{ij} = V_{i, j+1} = x_i^j \quad (j = 1, \ldots, n-1), \qquad (VZ_1^T)_{in} = 0.$$

    故 $(D_x V - VZ_1^T)_{ij} = 0$（$j = 1, \ldots, n-1$），$(D_x V - VZ_1^T)_{in} = x_i^n$。

    因此 $D_x V - VZ_1^T = \boldsymbol{x}^n \boldsymbol{e}_n^T$，其中 $\boldsymbol{x}^n = (x_1^n, x_2^n, \ldots, x_n^n)^T$。这是秩-1 矩阵。$\blacksquare$

!!! theorem "定理 37B.14 (循环矩阵的位移秩)"
    循环矩阵 $C = \sum_{k=0}^{n-1} c_k \Pi^k$（$\Pi$ 为循环置换矩阵）的 Stein 型位移
    $$\Delta(C) = C - \Pi C \Pi^T$$
    恒等于零矩阵。因此循环矩阵的位移秩为 $0$。

??? proof "证明"
    由循环矩阵与 $\Pi$ 可交换（$\Pi C = C \Pi$）和 $\Pi^T = \Pi^{n-1} = \Pi^{-1}$，
    $$\Pi C \Pi^T = \Pi C \Pi^{-1} = C.$$
    故 $C - \Pi C \Pi^T = 0$，秩为 $0$。$\blacksquare$

!!! example "例 37B.3"
    **位移秩总结表**。选取适当的位移算子，各类结构化矩阵的位移秩如下：

    | 矩阵类型 | 位移算子 $(F, F')$ | 位移方程 | 位移秩 |
    |----------|-------------------|----------|--------|
    | Toeplitz | $(Z_1, Z_1)$ | $T - Z_1 T Z_1^T$ | $\le 2$ |
    | Hankel | $(Z_1, Z_1)$ | $H - Z_1 H Z_1$ | $\le 2$ |
    | Vandermonde | $(D_x, Z_1)$ | $D_x V - V Z_1^T$ | $1$ |
    | Cauchy | $(D_x, D_y)$ | $D_x C - C D_y$ | $1$ |
    | 循环 | $(\Pi, \Pi)$ | $C - \Pi C \Pi^T$ | $0$ |
    | Toeplitz-like（逆）| $(Z_1, Z_1)$ | $T^{-1} - Z_1 T^{-1} Z_1^T$ | $\le 2$ |
    | Vandermonde-like（逆）| $(Z_1, D_x)$ | $Z_1^T V^{-1} - V^{-1} D_x$ | $1$ |

    这张表格清晰地展示了位移结构理论的核心洞察：所有经典结构化矩阵都具有**常数级别**的位移秩。

---

## 37B.6 基于位移结构的快速算法

### 广义 Schur 算法

!!! definition "定义 37B.12 (广义 Schur 算法)"
    **广义 Schur 算法**（Generalized Schur Algorithm, GSA）是利用位移结构进行 LU 分解（或 LDL 分解）的递归算法。其基本框架如下：

    **输入**：$n \times n$ 矩阵 $A$ 的位移生成元 $(G, H)$（$G, H \in \mathbb{C}^{n \times r}$），满足 $FA - AF' = GH^T$。

    **递归步骤**（$k = 1, 2, \ldots, n$）：

    1. 从当前生成元恢复 $A$ 的第 $k$ 行/列（$O(r^2)$ 运算）。
    2. 执行一步 Gauss 消元，得到 Schur 补。
    3. 更新位移生成元以表示 Schur 补（$O(rn)$ 运算），利用定理 37B.9(c)。

    **总复杂度**：$\sum_{k=1}^n O(rn) = O(rn^2)$。当 $r$ 为常数时（如 Toeplitz 的 $r = 2$），这给出 $O(n^2)$ 算法。

!!! theorem "定理 37B.15 (广义 Schur 算法的复杂度)"
    设 $A$ 为 $n \times n$ 矩阵，位移秩为 $r$。则：

    (a) 广义 Schur 算法在 $O(rn^2)$ 次算术运算内计算 $A$ 的 LU 分解的位移生成元。

    (b) 利用生成元，线性方程组 $A\boldsymbol{x} = \boldsymbol{b}$ 可在 $O(rn^2)$ 内求解。

    (c) 对 Toeplitz 矩阵（$r \le 2$），广义 Schur 算法的复杂度为 $O(n^2)$，恢复了 Levinson-Durbin 算法的复杂度。

??? proof "证明"
    (a) 每步需要 $O(r)$ 运算来恢复主元（从生成元重建），$O(rn)$ 运算来更新 $(n-k) \times r$ 的生成元。总计 $\sum_{k=1}^n O(r(n-k)) = O(rn^2/2) = O(rn^2)$。

    (b) 一旦 LU 分解的生成元已知，前代和回代分别需要 $O(rn)$ 和 $O(rn)$（因为 $L$ 和 $U$ 的位移秩也不超过 $r+1$），总计 $O(rn)$。但严格地，从生成元恢复完整的 $L, U$ 需要 $O(rn^2)$。

    (c) 对 Toeplitz 矩阵，$r \le 2$，故 $O(rn^2) = O(n^2)$。$\blacksquare$

### Gohberg-Semencul 公式

!!! theorem "定理 37B.16 (Gohberg-Semencul 公式)"
    设 $T$ 为 $n \times n$ 非奇异 Toeplitz 矩阵。设 $\boldsymbol{x} = (x_0, x_1, \ldots, x_{n-1})^T$ 为 $T\boldsymbol{x} = \boldsymbol{e}_1$ 的解，$\boldsymbol{y} = (y_0, y_1, \ldots, y_{n-1})^T$ 为 $T\boldsymbol{y} = \boldsymbol{e}_n$ 的解。若 $x_0 \ne 0$，则
    $$T^{-1} = \frac{1}{x_0}\left(L(\boldsymbol{x})\, U(\boldsymbol{x}') - L(J\boldsymbol{y})\, U(J\boldsymbol{y}')\right),$$
    其中：

    - $L(\boldsymbol{v})$ 是以 $\boldsymbol{v}$ 为第一列的下三角 Toeplitz 矩阵。
    - $U(\boldsymbol{v})$ 是以 $\boldsymbol{v}$ 为第一行的上三角 Toeplitz 矩阵。
    - $J$ 为 $n \times n$ 反序矩阵。
    - $\boldsymbol{x}' = (x_0, 0, \ldots, 0, x_{n-1}, x_{n-2}, \ldots, x_1)$（适当排列）。
    - $\boldsymbol{y}' = (y_{n-1}, 0, \ldots, 0, y_0, y_1, \ldots, y_{n-2})$（适当排列）。

    更简洁的形式为
    $$T^{-1} = \frac{1}{x_0}\bigl(\mathcal{L}(\boldsymbol{x})\,\mathcal{L}(\tilde{\boldsymbol{y}})^T - \mathcal{L}(Z_1 J\boldsymbol{y})\,\mathcal{L}(Z_1 J\tilde{\boldsymbol{x}})^T\bigr),$$
    其中 $\mathcal{L}(\boldsymbol{v})$ 是下三角 Toeplitz 矩阵，$\tilde{\boldsymbol{v}}$ 表示分量反序。

??? proof "证明"
    **核心思想**：Toeplitz 矩阵的逆具有位移秩 $\le 2$（关于 $(Z_1, Z_1)$），因此 $T^{-1}$ 可以用两对位移生成元来表示。

    **第一步**：由定理 37B.9(b)，$T^{-1}$ 的位移满足
    $$T^{-1} - Z_1^T T^{-1} Z_1 = T^{-1}(T - Z_1 T Z_1^T) T^{-1} \cdot (*).$$
    实际上需要更仔细地处理下移矩阵的转置关系。位移方程变为
    $$Z_1^T T^{-1} - T^{-1} Z_1 = -T^{-1}(Z_1^T T - T Z_1) T^{-1}.$$

    **第二步**：$Z_1^T T - T Z_1$ 可以直接计算。由 Toeplitz 结构，这是秩至多为 $2$ 的矩阵，其因式分解涉及 $T$ 的第一行和第一列的元素。

    **第三步**：设 $\boldsymbol{x} = T^{-1}\boldsymbol{e}_1$，$\boldsymbol{y} = T^{-1}\boldsymbol{e}_n$，则 $T^{-1}$ 的位移生成元可以用 $\boldsymbol{x}$ 和 $\boldsymbol{y}$ 表达。具体地，
    $$T^{-1} - Z_1^T T^{-1} Z_1 = \frac{1}{x_0}\left(\boldsymbol{x}\boldsymbol{x}_r^T - (J\boldsymbol{y})(J\boldsymbol{y})_r^T\right),$$
    其中下标 $r$ 表示某种移位操作。

    **第四步**：由位移方程恢复 $T^{-1}$。一个满足位移方程 $A - Z_1^T A Z_1 = \boldsymbol{p}\boldsymbol{q}^T$ 的下三角 Toeplitz 结构的矩阵恰为 $L(\boldsymbol{p})\,U(\boldsymbol{q})$。将两项位移生成元分别恢复，叠加即得 Gohberg-Semencul 公式。$\blacksquare$

!!! theorem "定理 37B.17 (Gohberg-Semencul 公式的应用)"
    Gohberg-Semencul 公式的实际意义是：

    (a) $T^{-1}\boldsymbol{v}$ 的计算归结为四次 Toeplitz 矩阵-向量乘法，每次可通过 FFT 在 $O(n\log n)$ 内完成，因此总复杂度为 $O(n\log n)$。

    (b) 一旦 $\boldsymbol{x}$ 和 $\boldsymbol{y}$ 已知（需 $O(n^2)$ 预处理或 $O(n\log^2 n)$ 超快预处理），后续的每次求解 $T^{-1}\boldsymbol{v}$ 仅需 $O(n\log n)$。

    (c) 该公式揭示了 Toeplitz 逆的内在结构：$T^{-1}$ 是两个秩-1（在 Toeplitz 意义下）矩阵乘积之差，这正是位移秩 $\le 2$ 的体现。

### 超快算法

!!! theorem "定理 37B.18 (超快 Toeplitz 求解器)"
    对称正定 Toeplitz 系统 $T\boldsymbol{x} = \boldsymbol{b}$ 可以在 $O(n\log^2 n)$ 次算术运算内求解。

??? proof "证明"
    **分治策略**：将 $n \times n$ Toeplitz 矩阵分为 $2 \times 2$ 块形式
    $$T = \begin{pmatrix} T_{11} & T_{12} \\ T_{21} & T_{22} \end{pmatrix},$$
    其中 $T_{11}, T_{22}$ 为 $n/2 \times n/2$ Toeplitz 矩阵，$T_{12}, T_{21}$ 为 $n/2 \times n/2$ Hankel 矩阵（因为 Toeplitz 矩阵的离对角线块在反序后变为 Toeplitz，即它们是 Hankel 的）。

    **关键观察**：Schur 补 $S = T_{22} - T_{21}T_{11}^{-1}T_{12}$ 不再是 Toeplitz 的，但由定理 37B.9(c)，$S$ 的位移秩仍然 $\le 2$。

    **生成元更新**：$T_{21}T_{11}^{-1}T_{12}$ 的计算涉及 Toeplitz 矩阵的逆（由 Gohberg-Semencul 公式表示）和 Hankel 矩阵的乘积。利用 FFT，位移生成元的更新可在 $O(n\log n)$ 内完成。

    **递推关系**：设 $T(n)$ 为求解 $n$ 阶 Toeplitz 系统的运算量，则
    $$T(n) = T(n/2) + O(n\log n).$$
    解此递推得 $T(n) = O(n\log^2 n)$。

    该算法的数值稳定性是一个重要的实践问题。Bini 和 Pan（1994）、Kailath 和 Sayed（1995）等人发展了不同版本的超快算法，其中某些变体在数值稳定性方面有所改善。$\blacksquare$

!!! example "例 37B.4"
    **各类算法的复杂度对比**：

    | 方法 | 复杂度 | 数值稳定性 |
    |------|--------|-----------|
    | Gauss 消元（一般矩阵） | $O(n^3)$ | 稳定（带部分选主元） |
    | Levinson-Durbin | $O(n^2)$ | 正定时稳定 |
    | 广义 Schur 算法 | $O(rn^2)$ | 一般稳定 |
    | Bini-Pan 超快算法 | $O(n\log^2 n)$ | 可能不稳定 |
    | 循环预处理 PCG | $O(n\log n)$（迭代） | 稳定 |

    对 $n = 10^6$，$O(n^3) \approx 10^{18}$（不可行），$O(n^2) \approx 10^{12}$（数小时），$O(n\log^2 n) \approx 4 \times 10^8$（秒级），$O(n\log n) \approx 2 \times 10^7$（毫秒级）。这展示了利用结构的巨大计算优势。

---

## 37B.7 层次矩阵与现代推广

经典位移结构理论适用于具有全局规律性的矩阵（如 Toeplitz、Cauchy）。然而，许多实际应用中的矩阵——如偏微分方程的 Green 函数离散化、积分方程的核矩阵——具有更一般的结构：虽然矩阵本身不是低秩的，但其**离对角线子块**具有低秩或近似低秩性。这种"分层低秩"结构催生了层次矩阵的理论和算法。

!!! definition "定义 37B.13 (HSS 矩阵)"
    **层次半可分矩阵**（Hierarchically Semi-Separable matrix, HSS 矩阵）是指满足以下性质的矩阵 $A \in M_n(\mathbb{C})$：

    (a) 将 $A$ 以某种嵌套的方式递归地分块为 $2 \times 2$ 形式。

    (b) 在每一层分块中，离对角线块（即 $A_{12}$ 和 $A_{21}$）具有低秩 $r$（$r$ 远小于块的维数）。

    (c) 不同层次的低秩表示之间具有**嵌套基**（nested basis）的关系。

    HSS 矩阵的存储量为 $O(rn)$，矩阵-向量乘法复杂度为 $O(rn)$，LU 分解复杂度为 $O(r^2 n)$。

!!! definition "定义 37B.14 ($\mathcal{H}$-矩阵)"
    **$\mathcal{H}$-矩阵**（Hierarchical matrix，Hackbusch 1999）是一种更一般的层次结构化矩阵。给定指标集 $I = \{1, \ldots, n\}$ 的一棵**簇树**（cluster tree），$\mathcal{H}$-矩阵 $A$ 的定义如下：

    (a) 将 $A$ 的子矩阵 $A_{t \times s}$（其中 $t, s$ 是簇树的节点）分为**可容许块**（admissible blocks）和**不可容许块**。

    (b) 可容许块用低秩逼近 $A_{t \times s} \approx U_t S_{ts} V_s^T$ 表示，其中 $U_t, V_s$ 为瘦矩阵，$S_{ts}$ 为小矩阵。

    (c) 不可容许块（通常为对角线上的小块）完整存储。

    $\mathcal{H}$-矩阵的关键参数是**秩** $r$（所有可容许块的最大秩）和**深度** $p$（簇树的深度，通常 $p = O(\log n)$）。

!!! theorem "定理 37B.19 ($\mathcal{H}$-矩阵的复杂度)"
    设 $A$ 为 $n \times n$ 的 $\mathcal{H}$-矩阵，秩参数为 $r$，簇树深度为 $p = O(\log n)$。则：

    (a) **存储量**：$O(rn\log n)$。

    (b) **矩阵-向量乘法**：$O(rn\log n)$。

    (c) **$\mathcal{H}$-矩阵加法**：$O(r^2 n\log n)$。

    (d) **近似 $\mathcal{H}$-LU 分解**：$O(r^2 n\log^2 n)$。

    对比：一般矩阵的存储为 $O(n^2)$，矩阵-向量乘法为 $O(n^2)$，LU 分解为 $O(n^3)$。

??? proof "证明"
    (a) 簇树每层有 $O(n)$ 个指标参与，每个可容许块贡献 $O(r \cdot \text{块大小})$ 的存储。由于可容许块在每层覆盖 $O(n)$ 个指标，每层存储 $O(rn)$，共 $p = O(\log n)$ 层，总计 $O(rn\log n)$。

    (b) 矩阵-向量乘法按簇树自底向上计算。对每个可容许块 $A_{ts} \approx U_t S_{ts} V_s^T$，$A_{ts}\boldsymbol{x}_s = U_t(S_{ts}(V_s^T \boldsymbol{x}_s))$ 的运算量为 $O(r \cdot |t|) + O(r \cdot |s|)$。每层总计 $O(rn)$，共 $O(\log n)$ 层。

    (c)(d) 通过递归的块消元和低秩更新，利用 $\mathcal{H}$-矩阵算术（包括低秩矩阵的截断 SVD 重压缩），可以在上述复杂度内完成。$\blacksquare$

!!! definition "定义 37B.15 (快速多极方法的联系)"
    **快速多极方法**（Fast Multipole Method, FMM，Greengard-Rokhlin 1987）是 $\mathcal{H}$-矩阵思想的先驱。FMM 用于加速 $N$ 体问题的计算：给定 $n$ 个源点 $\boldsymbol{y}_j$ 和 $n$ 个目标点 $\boldsymbol{x}_i$，计算
    $$\phi(\boldsymbol{x}_i) = \sum_{j=1}^n \frac{q_j}{|\boldsymbol{x}_i - \boldsymbol{y}_j|}, \quad i = 1, \ldots, n.$$

    直接计算需要 $O(n^2)$。FMM 利用核函数 $1/|\boldsymbol{x} - \boldsymbol{y}|$ 在远场区域的**多极展开**（低秩逼近），将复杂度降至 $O(n)$（固定精度下）或 $O(n\log n)$（自适应精度下）。

    从矩阵角度看，FMM 正是对 Cauchy 型矩阵 $(1/|\boldsymbol{x}_i - \boldsymbol{y}_j|))_{ij}$ 实施 $\mathcal{H}$-矩阵-向量乘法。

!!! theorem "定理 37B.20 (低秩离对角线结构的数学基础)"
    设 $K(\boldsymbol{x}, \boldsymbol{y})$ 为光滑核函数（如 $K(\boldsymbol{x}, \boldsymbol{y}) = 1/|\boldsymbol{x} - \boldsymbol{y}|$），$A_{ij} = K(\boldsymbol{x}_i, \boldsymbol{y}_j)$。对于**分离的**（well-separated）指标集 $t, s$（即 $\operatorname{dist}(X_t, Y_s) \ge c \cdot \operatorname{diam}(\max(X_t, Y_s))$），子矩阵 $A_{t \times s}$ 的奇异值指数衰减：
    $$\sigma_k(A_{t \times s}) \le C \cdot q^k, \quad 0 < q < 1,$$
    其中 $q$ 取决于分离比。因此 $A_{t \times s}$ 可以用秩-$r$ 矩阵以精度 $\epsilon$ 逼近，其中 $r = O(\log(1/\epsilon))$。

??? proof "证明"
    核心工具是核函数的**退化展开**（degenerate expansion）。当 $|\boldsymbol{x} - \boldsymbol{x}_0| < R$ 且 $|\boldsymbol{y} - \boldsymbol{y}_0| > R/q$（$q < 1$）时，
    $$K(\boldsymbol{x}, \boldsymbol{y}) = \sum_{k=0}^{\infty} \phi_k(\boldsymbol{x}) \psi_k(\boldsymbol{y}),$$
    其中展开系数以 $q^k$ 的速度衰减。这可以从 Taylor 展开（一维）或球谐函数展开（多维）推导。

    截断到 $r$ 项，残差为 $O(q^r)$。对应的矩阵逼近为
    $$A_{t \times s} \approx \sum_{k=0}^{r-1} \boldsymbol{u}_k \boldsymbol{v}_k^T,$$
    秩为 $r$。取 $r = O(\log(1/\epsilon)/\log(1/q))$ 即可达到精度 $\epsilon$。$\blacksquare$

!!! note "注记 37B.3 (从位移结构到层次结构)"
    经典位移结构理论与层次矩阵理论之间有深刻的联系：

    1. **共同思想**：两者都利用矩阵的"内在低秩性"来设计快速算法。位移结构利用全局位移方程的低秩性，层次矩阵利用离对角线子块的低秩性。

    2. **Cauchy 矩阵的桥梁作用**：Cauchy 矩阵既有位移秩 1（经典理论），其离对角线子块又具有低秩性（层次理论）。事实上，任何位移秩为 $r$ 的矩阵都可以通过适当的变换化为 Cauchy-like 矩阵，然后应用层次矩阵技术。

    3. **现代统一框架**：Chandrasekaran-Gu-Pals (2006) 等人发展的 HSS 算法框架，可以看作位移结构广义 Schur 算法的层次化推广。

!!! note "注记 37B.4 (计算复杂度的层次)"
    结构化矩阵的快速算法复杂度形成清晰的层次：

    $$O(n^3) \xrightarrow{\text{位移结构}} O(rn^2) \xrightarrow{\text{超快算法}} O(rn\log^2 n) \xrightarrow{\text{层次矩阵}} O(r^2 n \log^k n) \xrightarrow{\text{FMM}} O(n)$$

    每一步利用矩阵的更精细的结构信息来获得更低的复杂度。在实践中，选择哪个层次的算法取决于问题规模 $n$、位移秩 $r$、精度要求和实现复杂度之间的权衡。

---

## 37B.8 习题

!!! example "例 37B.5 (习题 1)"
    计算以下 Vandermonde 矩阵的行列式：
    $$V = V(1, 2, 3, 4) = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & 2 & 4 & 8 \\ 1 & 3 & 9 & 27 \\ 1 & 4 & 16 & 64 \end{pmatrix}.$$
    验证 Vandermonde 行列式公式 $\det V = \prod_{i<j}(x_j - x_i)$。进一步计算 $V$ 的 2-范数条件数 $\kappa_2(V)$。

!!! example "例 37B.6 (习题 2)"
    设 $\omega = e^{2\pi i/n}$，证明 Vandermonde 矩阵 $V(\omega^0, \omega^1, \ldots, \omega^{n-1})$ 与 DFT 矩阵的关系为 $V = \sqrt{n} \cdot F^*$，其中 $F$ 为 $n$ 阶 DFT 矩阵。由此推导循环矩阵可以被 DFT 对角化。

!!! example "例 37B.7 (习题 3)"
    证明 Cauchy 行列式公式：对 $n = 3$，直接展开
    $$\det\begin{pmatrix} \frac{1}{x_1+y_1} & \frac{1}{x_1+y_2} & \frac{1}{x_1+y_3} \\ \frac{1}{x_2+y_1} & \frac{1}{x_2+y_2} & \frac{1}{x_2+y_3} \\ \frac{1}{x_3+y_1} & \frac{1}{x_3+y_2} & \frac{1}{x_3+y_3} \end{pmatrix}$$
    并化简为 $\frac{\prod_{i<j}(x_j-x_i)\prod_{i<j}(y_j-y_i)}{\prod_{i,j}(x_i+y_j)}$。

!!! example "例 37B.8 (习题 4)"
    设 $C = (1/(x_i - y_j))$ 为 $n \times n$ Cauchy 矩阵。证明 $C$ 的逆 $C^{-1}$ 可以显式表示为
    $$(C^{-1})_{ij} = (x_j - y_i) \cdot \frac{\prod_{k=1}^n (y_i - y_k)}{{\prod_{k \ne i}(y_i - y_k)}} \cdot \frac{\prod_{k=1}^n (x_j - x_k)}{{\prod_{k \ne j}(x_j - x_k)}} \cdot \frac{1}{(x_j - y_i)}.$$
    化简此表达式。（提示：利用 Lagrange 插值基函数。）

!!! example "例 37B.9 (习题 5)"
    验证 Toeplitz 矩阵的位移秩计算。设
    $$T = \begin{pmatrix} 2 & -1 & 0 \\ 3 & 2 & -1 \\ 1 & 3 & 2 \end{pmatrix}.$$
    显式计算 $\Delta(T) = T - Z_1 T Z_1^T$，验证其秩为 $2$，并找出位移生成元 $(G, H)$ 使得 $\Delta(T) = GH^T$（$G, H \in \mathbb{R}^{3 \times 2}$）。

!!! example "例 37B.10 (习题 6)"
    证明：若 $A$ 的 Sylvester 型位移秩为 $r$（$FA - AF' = GH^T$），则 $A^{-1}$ 的 Sylvester 型位移秩（关于算子 $(F', F)$）也为 $r$。具体给出 $A^{-1}$ 的位移生成元。

!!! example "例 37B.11 (习题 7)"
    **Resultant 的计算**。设 $f(x) = x^3 - 1$，$g(x) = x^2 - 1$。

    (a) 写出 Sylvester 结式矩阵 $\operatorname{Syl}(f, g)$。

    (b) 计算 $\operatorname{Res}(f, g) = \det \operatorname{Syl}(f, g)$。

    (c) 验证 $\operatorname{Res}(f, g) = \prod_{f(\alpha)=0} g(\alpha)$，其中乘积取遍 $f$ 的所有根。

    (d) 确定 $\gcd(f, g)$ 并验证 $\operatorname{rank}(\operatorname{Syl}(f,g)) = \deg f + \deg g - \deg(\gcd(f,g))$。

!!! example "例 37B.12 (习题 8)"
    **Bezout 矩阵的构造**。设 $f(x) = x^3 + 1$，$g(x) = x^2 + x + 1$。

    (a) 计算 $\frac{f(x)g(y) - f(y)g(x)}{x - y}$ 并由此构造 $3 \times 3$ Bezout 矩阵 $B(f,g)$。

    (b) 验证 $B(f,g)$ 的秩等于 $3 - \deg(\gcd(f,g))$。

    (c) 利用 $B(f,g)$ 的秩亏损来确定 $\gcd(f,g)$。

!!! example "例 37B.13 (习题 9)"
    **Krylov 矩阵与最小多项式**。设 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$\boldsymbol{b} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$。

    (a) 计算 Krylov 矩阵 $K(A, \boldsymbol{b}) = (\boldsymbol{b}, A\boldsymbol{b})$。

    (b) $K(A, \boldsymbol{b})$ 是否可逆？这与 $A$ 关于 $\boldsymbol{b}$ 的最小多项式有什么关系？

    (c) 取 $\boldsymbol{b}' = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$，重复上述计算。讨论向量选择对 Krylov 矩阵秩的影响。

!!! example "例 37B.14 (习题 10)"
    **Gohberg-Semencul 公式的验证**。设 $T = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。

    (a) 求解 $T\boldsymbol{x} = \boldsymbol{e}_1$ 和 $T\boldsymbol{y} = \boldsymbol{e}_2$。

    (b) 利用 Gohberg-Semencul 公式写出 $T^{-1}$。

    (c) 直接计算 $T^{-1}$ 并与 (b) 的结果对比验证。

!!! example "例 37B.15 (习题 11)"
    **位移秩与矩阵变换**。证明：Toeplitz 矩阵 $T$ 可以通过预乘和后乘适当的对角矩阵/Vandermonde 矩阵转化为 Cauchy-like 矩阵。具体地，若 $D_\omega = \operatorname{diag}(1, \omega, \omega^2, \ldots, \omega^{n-1})$（$\omega$ 为 $n$ 次本原单位根），则 $\hat{T} = FTF^*$（$F$ 为 DFT 矩阵）满足什么位移方程？

    （提示：若 $T - Z_1 T Z_1^T = GH^T$，对两边做 $F(\cdot)F^*$ 变换。）

!!! example "例 37B.16 (习题 12)"
    **$\mathcal{H}$-矩阵逼近**。考虑 $n \times n$ 矩阵 $A_{ij} = 1/(i + j - 1)$（Hilbert 矩阵）。

    (a) 证明 Hilbert 矩阵是 Cauchy 矩阵的特殊情形（指定参数 $x_i, y_j$）。

    (b) 对于"分离的"指标集 $t = \{i : p \le i \le p+k-1\}$，$s = \{j : q \le j \le q+l-1\}$（$|p - q| \gg \max(k, l)$），说明子矩阵 $A_{t \times s}$ 的奇异值指数衰减的原因。

    (c) 对 $n = 8$，画出 $\mathcal{H}$-矩阵分块的示意图，标出哪些块可以用低秩逼近。

!!! example "例 37B.17 (习题 13)"
    **合流 Vandermonde 矩阵与 Hermite 插值**。设节点 $x_1$ 出现 $m_1$ 次，$x_2$ 出现 $m_2$ 次（$m_1 + m_2 = n$），定义合流 Vandermonde 矩阵
    $$V_c = \begin{pmatrix}
    1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\
    0 & 1 & 2x_1 & \cdots & (n-1)x_1^{n-2} \\
    \vdots & & & & \vdots \\
    1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\
    0 & 1 & 2x_2 & \cdots & (n-1)x_2^{n-2} \\
    \vdots & & & & \vdots
    \end{pmatrix}.$$

    (a) 证明 $\det V_c \ne 0$ 当且仅当 $x_1 \ne x_2$。

    (b) 说明 Hermite 插值问题（给定函数值和导数值，求插值多项式）等价于求解合流 Vandermonde 系统 $V_c \boldsymbol{a} = \boldsymbol{y}$。

!!! example "例 37B.18 (习题 14)"
    **位移结构的极限情形**。

    (a) 证明：位移秩为 $0$ 的矩阵（关于算子 $(Z_f, Z_f)$，$Z_f$ 为循环移位矩阵）恰好是循环矩阵。

    (b) 位移秩为 $n$ 的矩阵是什么？它与一般矩阵有何关系？

    (c) 讨论位移秩 $r$ 在 $0$ 与 $n$ 之间变化时，矩阵从"完全结构化"到"完全无结构"的过渡。
