# 第 37A 章 Toeplitz、Hankel 与循环矩阵

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 矩阵分解(Ch10) · 数值线性代数(Ch22) · Fourier 分析基础

**本章脉络**：Toeplitz 矩阵定义与基本性质 $\to$ Levinson-Durbin 算法 $\to$ 正定 Toeplitz 与谱分解 $\to$ Hankel 矩阵与矩问题 $\to$ 循环矩阵与 DFT 对角化 $\to$ 循环代数与环同构 $\to$ 分块 Toeplitz 矩阵 $\to$ Szegő 分布定理 $\to$ 循环预处理

**延伸**：Toeplitz 矩阵在信号处理（自相关矩阵）、时间序列分析（AR 模型）、控制理论中无处不在；循环矩阵通过 DFT 对角化为快速算法奠定基础；Szegő 定理连接了 Toeplitz 行列式的渐近与调和分析中的深刻问题

</div>

一般的 $n \times n$ 矩阵有 $n^2$ 个自由参数，求解线性方程组需要 $O(n^3)$ 次运算。然而，在信号处理、时间序列分析和偏微分方程等领域中反复出现的矩阵具有特殊的**结构**——元素之间存在确定的关系，使得矩阵可以用远少于 $n^2$ 个参数来描述。本章研究三种最基本的结构化矩阵类——Toeplitz 矩阵、Hankel 矩阵和循环矩阵，它们以"沿对角线常数"或"沿反对角线常数"的简洁模式编码了丰富的数学结构。

---

## 37A.1 Toeplitz 矩阵的定义与基本性质

!!! definition "定义 37A.1 (Toeplitz 矩阵)"
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

!!! definition "定义 37A.2 (对称 Toeplitz 矩阵与 Hermite Toeplitz 矩阵)"
    (a) 若 $t_k = t_{-k}$ 对所有 $k$，则 $T$ 为**对称 Toeplitz 矩阵**。此时 $T$ 由 $n$ 个参数 $t_0, t_1, \ldots, t_{n-1}$ 完全确定。

    (b) 若 $t_{-k} = \overline{t_k}$ 对所有 $k$，则 $T$ 为 **Hermite Toeplitz 矩阵**。此时 $T = T^*$。

!!! theorem "定理 37A.1 (Toeplitz 矩阵的基本性质)"
    设 $T_1, T_2 \in M_n(\mathbb{C})$ 为 Toeplitz 矩阵。则：

    (a) $\alpha T_1 + \beta T_2$ 是 Toeplitz 矩阵（Toeplitz 矩阵对线性组合封闭）；

    (b) $T_1^T$ 是 Toeplitz 矩阵（将 $t_k$ 替换为 $t_{-k}$）；

    (c) $T_1 T_2$ **一般不是** Toeplitz 矩阵（Toeplitz 矩阵对乘法不封闭）；

    (d) 若 $T$ 为对称 Toeplitz（$t_k = t_{-k}$），则 $T$ 具有**持续结构**（persymmetry）：$JTJ = T$，其中 $J$ 为 $n \times n$ 反对角单位矩阵（逆序置换矩阵）；

    (e) $n$ 阶 Toeplitz 矩阵全体构成 $M_n(\mathbb{C})$ 的一个 $(2n-1)$ 维线性子空间。

??? proof "证明"
    (a) 线性组合的 $(i,j)$ 元素 $\alpha t^{(1)}_{i-j} + \beta t^{(2)}_{i-j}$ 仅取决于 $i-j$。

    (b) $(T^T)_{ij} = t_{ji} = t_{j-i} = t_{-(i-j)}$，定义 $\tilde{t}_k = t_{-k}$，则 $T^T$ 是以 $\{\tilde{t}_k\}$ 为参数的 Toeplitz 矩阵。

    (c) 取 $n=3$，$T_1$ 以第一行 $(1,2,3)$、第一列 $(1,0,0)$ 定义，$T_2$ 以第一行 $(1,0,0)$、第一列 $(1,1,1)$ 定义，乘积一般不具有 Toeplitz 结构。

    (d) $J$ 是置换矩阵 $J = (e_n, e_{n-1}, \ldots, e_1)$。$(JTJ)_{ij} = t_{n+1-i, n+1-j} = t_{(n+1-i)-(n+1-j)} = t_{j-i} = t_{-(i-j)}$。当 $t_k = t_{-k}$ 时，$t_{-(i-j)} = t_{i-j}$，故 $JTJ = T$。

    (e) 由 (a) 显然，Toeplitz 矩阵对线性组合封闭。维数为 $2n-1$，因为 $T$ 由 $t_{-(n-1)}, \ldots, t_0, \ldots, t_{n-1}$ 这 $2n-1$ 个自由参数决定，且这些参数可以独立选取。

!!! example "例 37A.1"
    下列矩阵是 Toeplitz 矩阵：
    $$T_1 = \begin{pmatrix} 2 & -1 & 0 \\ 3 & 2 & -1 \\ 1 & 3 & 2 \end{pmatrix}, \quad T_2 = \begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 1 \\ 1 & 1 & 1 \end{pmatrix}.$$
    $T_1$ 不对称（$t_1 = 3 \ne t_{-1} = -1$），$T_2$ 是对称的。$T_2$ 的秩为 $1$，说明 Toeplitz 矩阵可以是奇异的。

!!! note "注记 37A.1 (Toeplitz 矩阵在信号处理中的出现)"
    宽平稳随机过程 $\{X_t\}$ 的自协方差矩阵
    $$R_{ij} = \operatorname{Cov}(X_i, X_j) = r(|i-j|)$$
    天然具有 Toeplitz 结构。Levinson-Durbin 算法最初就是在信号处理和时间序列分析的背景下被开发的，用于自回归（AR）模型的参数估计。在功率谱估计中，Toeplitz 矩阵编码了信号的二阶统计特性。

---

## 37A.2 Levinson-Durbin 算法

Levinson-Durbin 算法是求解对称正定 Toeplitz 系统最经典的快速算法，将复杂度从 $O(n^3)$ 降至 $O(n^2)$。其核心思想是利用 Toeplitz 结构，通过递推从小规模子系统的解逐步构建大规模系统的解。

!!! definition "定义 37A.3 (Yule-Walker 方程与反射系数)"
    设 $T_n = (t_{|i-j|})_{1 \le i,j \le n}$ 为 $n$ 阶对称正定 Toeplitz 矩阵。**Yule-Walker 方程**为
    $$T_k \boldsymbol{a}^{(k)} = -\boldsymbol{t}_k, \quad k = 1, 2, \ldots,$$
    其中 $\boldsymbol{t}_k = (t_1, t_2, \ldots, t_k)^T$。向量 $\boldsymbol{a}^{(k)} = (a_1^{(k)}, \ldots, a_k^{(k)})^T$ 是 $k$ 阶 AR 模型的系数。

    第 $k$ 步的**反射系数**（又称偏相关系数或 PARCOR 系数）定义为
    $$\alpha_k = a_k^{(k)}.$$

!!! theorem "定理 37A.2 (Levinson-Durbin 递推)"
    设 $T_n$ 为 $n$ 阶对称正定 Toeplitz 矩阵，$t_0 > 0$。Levinson-Durbin 算法按以下步骤递推：

    **初始化**：$\sigma_0^2 = t_0$。

    **递推**（$k = 1, 2, \ldots, n-1$）：

    (a) 计算反射系数：
    $$\alpha_k = -\frac{t_k + \sum_{j=1}^{k-1} a_j^{(k-1)} t_{k-j}}{\sigma_{k-1}^2}.$$

    (b) 更新 AR 系数：
    $$a_j^{(k)} = a_j^{(k-1)} + \alpha_k a_{k-j}^{(k-1)}, \quad j = 1, \ldots, k-1, \qquad a_k^{(k)} = \alpha_k.$$

    (c) 更新预测误差功率：
    $$\sigma_k^2 = \sigma_{k-1}^2 (1 - \alpha_k^2).$$

    算法的总运算量为 $O(n^2)$。

??? proof "证明"
    **正确性证明**：我们需要证明若 $T_{k-1}\boldsymbol{a}^{(k-1)} = -\boldsymbol{t}_{k-1}$，则按上述递推得到的 $\boldsymbol{a}^{(k)}$ 满足 $T_k\boldsymbol{a}^{(k)} = -\boldsymbol{t}_k$。

    注意到 $T_k$ 可以写成分块形式：
    $$T_k = \begin{pmatrix} T_{k-1} & \boldsymbol{t}_{k-1}' \\ (\boldsymbol{t}_{k-1}')^T & t_0 \end{pmatrix},$$
    其中 $\boldsymbol{t}_{k-1}' = (t_{k-1}, t_{k-2}, \ldots, t_1)^T = J_{k-1}\boldsymbol{t}_{k-1}$。

    定义 $\bar{\boldsymbol{a}}^{(k-1)} = J_{k-1}\boldsymbol{a}^{(k-1)}$（反转系数向量）。由对称 Toeplitz 的持续结构 $J T_{k-1} J = T_{k-1}$，有
    $$T_{k-1}\bar{\boldsymbol{a}}^{(k-1)} = -\boldsymbol{t}_{k-1}'.$$

    现在构造 $\boldsymbol{a}^{(k)} = \begin{pmatrix} \boldsymbol{a}^{(k-1)} \\ 0 \end{pmatrix} + \alpha_k \begin{pmatrix} \bar{\boldsymbol{a}}^{(k-1)} \\ 1 \end{pmatrix}$。验证 $T_k \boldsymbol{a}^{(k)} = -\boldsymbol{t}_k$ 即得递推关系中 $\alpha_k$ 的表达式。

    **复杂度分析**：第 $k$ 步需要计算 $\alpha_k$（$O(k)$ 次运算）和更新 $\boldsymbol{a}^{(k)}$（$O(k)$ 次运算）。总运算量 $\sum_{k=1}^{n-1} O(k) = O(n^2)$。

!!! theorem "定理 37A.3 (反射系数的性质)"
    设 $T_n$ 为正定 Toeplitz 矩阵。则：

    (a) $|\alpha_k| < 1$ 对所有 $k = 1, 2, \ldots, n-1$。

    (b) $\sigma_k^2 > 0$ 对所有 $k$（预测误差严格为正）。

    (c) $T_n$ 正定当且仅当 $|\alpha_k| < 1$ 对所有 $k$。

    (d) $\det(T_k) = t_0^k \prod_{j=1}^{k-1}(1-\alpha_j^2)^{k-j}$。

??? proof "证明"
    (a) 由 $\sigma_k^2 = \sigma_{k-1}^2(1-\alpha_k^2)$ 和 $\sigma_k^2 > 0$（因 $T_k$ 正定），得 $1 - \alpha_k^2 > 0$，即 $|\alpha_k| < 1$。

    (b) 由 (a) 和递推关系 $\sigma_k^2 = \sigma_{k-1}^2(1-\alpha_k^2)$ 立即得到。

    (c) 充分性：若 $|\alpha_k| < 1$ 对所有 $k$，则 $\sigma_k^2 > 0$ 对所有 $k$，而 $\sigma_k^2 = \det(T_{k+1})/\det(T_k)$，故所有顺序主子式为正，$T_n$ 正定。必要性即 (a)。

    (d) 由 $\sigma_k^2 = \det(T_{k+1})/\det(T_k)$ 递推即得。

!!! example "例 37A.2"
    求解 Toeplitz 系统
    $$\begin{pmatrix} 4 & 2 & 1 \\ 2 & 4 & 2 \\ 1 & 2 & 4 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}.$$

    矩阵为对称正定 Toeplitz 矩阵，$t_0 = 4, t_1 = 2, t_2 = 1$。

    **步骤 1** ($k=1$)：$\alpha_1 = -t_1/t_0 = -2/4 = -1/2$。$\sigma_1^2 = 4(1-1/4) = 3$。$a_1^{(1)} = -1/2$。

    **步骤 2** ($k=2$)：
    $$\alpha_2 = -\frac{t_2 + a_1^{(1)} t_1}{\sigma_1^2} = -\frac{1 + (-1/2)\cdot 2}{3} = -\frac{1-1}{3} = 0.$$

    $a_1^{(2)} = a_1^{(1)} + \alpha_2 a_1^{(1)} = -1/2 + 0 = -1/2$，$a_2^{(2)} = \alpha_2 = 0$。$\sigma_2^2 = 3(1-0) = 3$。

    然后利用 Toeplitz 结构回代求解 $Tx = b$，最终得 $x = (0, 1/4, 3/4)^T$。

    使用一般 Gauss 消元需要约 $n^3/3 = 9$ 次乘法，而 Levinson 算法仅需约 $n^2 = 9$ 次——在此小例中差别不明显，但对 $n = 10000$，从 $3 \times 10^{11}$ 减少到 $10^8$，加速约 3000 倍。

!!! note "注记 37A.2 (Schur 算法)"
    与 Levinson 算法并行的是 **Schur 算法**（也称 Schur-Cohn 算法），它直接计算 Toeplitz 矩阵的 LDL 分解的生成元而不显式求解方程组。Schur 算法同为 $O(n^2)$ 复杂度，但在并行计算中更有优势，因为其计算步骤之间的依赖性更低。

---

## 37A.3 正定 Toeplitz 矩阵与谱分解

!!! theorem "定理 37A.4 (正定 Toeplitz 矩阵的刻画)"
    设 $T_n = (t_{|i-j|})_{1 \le i,j \le n}$。以下条件等价：

    (a) $T_n$ 正定。

    (b) 存在 $L^1([0,2\pi])$ 上的非负函数 $f \ge 0$（$f$ 不几乎处处为零），使得
    $$t_k = \frac{1}{2\pi}\int_0^{2\pi} f(\theta) e^{-ik\theta}\, d\theta, \quad k = 0, \pm 1, \ldots, \pm(n-1).$$

    (c) 所有反射系数满足 $|\alpha_k| < 1$，$k = 1, \ldots, n-1$。

??? proof "证明"
    **(a) $\Rightarrow$ (b)**：这是 Herglotz 定理（有限维版本）的推论。正定 Toeplitz 矩阵的参数序列是正定序列，由 Bochner 定理，存在非负测度 $\mu$ 使得 $t_k = \int e^{-ik\theta}\, d\mu(\theta)$。当 $T_n$ 严格正定时，$\mu$ 有密度函数 $f \ge 0$。

    **(b) $\Rightarrow$ (a)**：对任意非零向量 $\boldsymbol{x} = (x_1, \ldots, x_n)^T$，
    $$\boldsymbol{x}^* T_n \boldsymbol{x} = \sum_{j,k} t_{j-k} x_j \overline{x_k} = \frac{1}{2\pi}\int_0^{2\pi} f(\theta) \left|\sum_{k=1}^n x_k e^{ik\theta}\right|^2 d\theta \ge 0.$$
    等号成立要求 $f(\theta) \cdot |p(\theta)|^2 = 0$ a.e.，当 $f$ 不几乎处处为零时，这要求 $p(\theta) = \sum_k x_k e^{ik\theta} = 0$ a.e.，从而 $\boldsymbol{x} = 0$。

    **(a) $\Leftrightarrow$ (c)**：已在定理 37A.3 中证明。

!!! definition "定义 37A.4 (谱分解)"
    设 $f(\theta) > 0$ 为正定 Toeplitz 矩阵族的生成函数，且 $\log f \in L^1$。**谱分解**（spectral factorization）是寻找函数
    $$\sigma(z) = \sigma_0 + \sigma_1 z + \sigma_2 z^2 + \cdots, \quad \sigma_0 > 0,$$
    在单位圆盘 $|z| < 1$ 内解析且无零点，使得
    $$f(\theta) = |\sigma(e^{i\theta})|^2.$$

!!! theorem "定理 37A.5 (Szegő-Kolmogorov 谱分解定理)"
    设 $f \ge 0$，$f \in L^1([0,2\pi])$。谱分解 $f(\theta) = |\sigma(e^{i\theta})|^2$（$\sigma$ 在单位圆盘内解析且无零点）存在当且仅当
    $$\frac{1}{2\pi}\int_0^{2\pi} \log f(\theta)\, d\theta > -\infty$$
    （即 $\log f \in L^1$，称为 **Szegő 条件**）。此时
    $$\sigma_0^2 = \exp\left(\frac{1}{2\pi}\int_0^{2\pi}\log f(\theta)\, d\theta\right).$$

??? proof "证明"
    **必要性**：若 $f = |\sigma|^2$，则 $\log f = 2\operatorname{Re}(\log \sigma)$。由 $\log \sigma$ 在单位圆盘内解析且 $\operatorname{Re}(\log \sigma)$ 在 $L^1$ 中，$\log f \in L^1$。由 Jensen 公式（或调和函数的均值性质），
    $$\frac{1}{2\pi}\int_0^{2\pi}\log f(\theta)\, d\theta = 2\log|\sigma(0)| = 2\log\sigma_0 > -\infty.$$

    **充分性**：定义 $\log \sigma(z) = \frac{1}{4\pi}\int_0^{2\pi}\frac{e^{i\theta}+z}{e^{i\theta}-z}\log f(\theta)\, d\theta$（Poisson 积分的复化形式）。可以验证这定义了单位圆盘内的解析函数，且 $\operatorname{Re}(\log\sigma(e^{i\theta})) = \frac{1}{2}\log f(\theta)$ a.e.，从而 $|\sigma(e^{i\theta})|^2 = f(\theta)$。

!!! note "注记 37A.3 (谱分解与线性预测)"
    谱分解在信号处理中的意义是：$\sigma_0^2$ 等于最优线性预测器的均方误差（即从过去全部历史预测当前值的最小可能误差）。Szegő 条件 $\log f \in L^1$ 恰好是过程"不完全可预测"的充要条件——当它不成立时（即 $\int \log f = -\infty$），过程是**确定性的**，可以被完美预测。

---

## 37A.4 Hankel 矩阵

!!! definition "定义 37A.5 (Hankel 矩阵)"
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
    Hankel 矩阵同样由 $2n-1$ 个参数确定。Hankel 矩阵天然是**对称的**：$h_{ij} = h_{i+j} = h_{j+i} = h_{ji}$。

!!! theorem "定理 37A.6 (Toeplitz 与 Hankel 的互换)"
    设 $J$ 为 $n \times n$ 反对角单位矩阵。则：

    (a) 若 $T$ 为 Toeplitz 矩阵，则 $TJ$ 和 $JT$ 为 Hankel 矩阵；

    (b) 若 $H$ 为 Hankel 矩阵，则 $HJ$ 和 $JH$ 为 Toeplitz 矩阵。

    换言之，Toeplitz 和 Hankel 矩阵通过左乘或右乘反序矩阵 $J$ 相互转化。

??? proof "证明"
    (a) $(TJ)_{ij} = \sum_k T_{ik} J_{kj} = T_{i, n+1-j} = t_{i-(n+1-j)} = t_{i+j-(n+1)}$，仅取决于 $i+j$，故 $TJ$ 是 Hankel 的。$JT$ 的证明类似：$(JT)_{ij} = T_{n+1-i,j} = t_{(n+1-i)-j} = t_{(n+1)-(i+j)}$，仅取决于 $i+j$。

    (b) $(HJ)_{ij} = H_{i, n+1-j} = h_{i+(n+1-j)} = h_{(n+1) + (i-j)}$，仅取决于 $i-j$，故 $HJ$ 是 Toeplitz 的。$JH$ 的证明对称。

!!! definition "定义 37A.6 (Hankel 矩阵与矩问题)"
    给定测度 $\mu$ 在 $\mathbb{R}$ 上的矩 $s_k = \int x^k \, d\mu(x)$，$k = 0, 1, 2, \ldots$，对应的 **Hankel 矩阵**为
    $$H_n = (s_{i+j-2})_{1 \le i,j \le n} = \begin{pmatrix}
    s_0 & s_1 & \cdots & s_{n-1} \\
    s_1 & s_2 & \cdots & s_n \\
    \vdots & & & \vdots \\
    s_{n-1} & s_n & \cdots & s_{2n-2}
    \end{pmatrix}.$$

!!! theorem "定理 37A.7 (Hamburger 矩问题)"
    给定实数序列 $\{s_k\}_{k=0}^{\infty}$。以下条件等价：

    (a) 存在 $\mathbb{R}$ 上的非负 Borel 测度 $\mu$ 使得 $s_k = \int_{-\infty}^{\infty} x^k\, d\mu(x)$ 对所有 $k \ge 0$。

    (b) 所有 Hankel 矩阵 $H_n = (s_{i+j-2})_{1\le i,j\le n}$ 半正定，即 $H_n \succeq 0$ 对所有 $n \ge 1$。

    进一步，$\mu$ 的支撑有恰好 $r$ 个点当且仅当 $\operatorname{rank}(H_n) = r$ 对所有充分大的 $n$。

??? proof "证明"
    **(a) $\Rightarrow$ (b)**：对任意 $\boldsymbol{c} = (c_0, \ldots, c_{n-1})^T$，
    $$\boldsymbol{c}^T H_n \boldsymbol{c} = \sum_{i,j} c_i c_j s_{i+j} = \int \left(\sum_i c_i x^i\right)^2 d\mu(x) \ge 0.$$

    **(b) $\Rightarrow$ (a)**：这是经典的 Hamburger 矩问题的解。由半正定条件，可以构造正交多项式序列和对应的 Jacobi 矩阵，然后由谱定理得到测度 $\mu$。

!!! example "例 37A.3"
    考虑标准正态分布 $\mu = N(0,1)$，其矩为 $s_{2k} = (2k-1)!! = 1 \cdot 3 \cdots (2k-1)$，$s_{2k+1} = 0$。

    $3 \times 3$ Hankel 矩阵：
    $$H_3 = \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 3 \end{pmatrix}.$$
    特征值为 $1$ 和 $2 \pm \sqrt{2}$，全为正。$H_3 > 0$ 与 $N(0,1)$ 为正测度一致。

!!! note "注记 37A.4 (Hankel 算子与系统理论)"
    在线性系统理论中，传递函数 $G(s) = C(sI - A)^{-1}B$ 对应的 Hankel 算子的矩阵表示恰好是由 Markov 参数 $h_k = CA^{k-1}B$ 构成的 Hankel 矩阵。Hankel 矩阵的秩等于系统的 **McMillan 度**（最小实现的状态维数）。Hankel 奇异值在模型降阶（平衡截断，Glover 1984）中起核心作用——AAK 定理（Adamyan-Arov-Kreĭn）精确刻画了 Hankel 算子的最佳有理逼近。

!!! theorem "定理 37A.8 (Kronecker 定理)"
    无穷 Hankel 矩阵 $H_\infty = (h_{i+j})_{i,j \ge 0}$ 具有有限秩 $r$ 当且仅当对应的 Laurent 级数
    $$\sum_{k=0}^{\infty} h_k z^{-k-1}$$
    是一个次数至多为 $r$ 的有理函数。

??? proof "证明"
    **必要性**：若 $\operatorname{rank}(H_\infty) = r$，则列 $r+1$ 是前 $r$ 列的线性组合。这给出递推关系 $h_{k+r} = \sum_{j=0}^{r-1} c_j h_{k+j}$，其特征多项式的次数为 $r$。因此 $\{h_k\}$ 满足 $r$ 阶线性递推，对应的生成函数为有理函数 $p(z)/q(z)$，$\deg q \le r$。

    **充分性**：若 Laurent 级数是 $r$ 次有理函数，则系数 $h_k$ 满足 $r$ 阶线性递推，故 $H_\infty$ 的列空间至多 $r$ 维。

---

## 37A.5 循环矩阵与 DFT 对角化

!!! definition "定义 37A.7 (循环矩阵)"
    矩阵 $C \in M_n(\mathbb{C})$ 称为**循环矩阵**，若每一行是上一行的循环右移：
    $$C = \begin{pmatrix}
    c_0 & c_{n-1} & c_{n-2} & \cdots & c_1 \\
    c_1 & c_0 & c_{n-1} & \cdots & c_2 \\
    c_2 & c_1 & c_0 & \cdots & c_3 \\
    \vdots & & & \ddots & \vdots \\
    c_{n-1} & c_{n-2} & c_{n-3} & \cdots & c_0
    \end{pmatrix} = \sum_{k=0}^{n-1} c_k \Pi^k,$$
    其中 $\Pi$ 为**循环置换矩阵**（$\Pi e_j = e_{j+1 \pmod{n}}$）。循环矩阵由 $n$ 个参数确定。循环矩阵是特殊的 Toeplitz 矩阵，附加了周期条件 $t_{k} = t_{k+n}$。

!!! definition "定义 37A.8 (DFT 矩阵)"
    设 $\omega = e^{2\pi i/n}$ 为 $n$ 次本原单位根。$n$ 阶 **DFT（离散 Fourier 变换）矩阵**定义为
    $$F = \frac{1}{\sqrt{n}} (\omega^{(j-1)(k-1)})_{1 \le j,k \le n}.$$
    $F$ 是酉矩阵：$F^* F = I$。

!!! theorem "定理 37A.9 (循环矩阵的 DFT 对角化)"
    每个循环矩阵 $C$ 可以被 DFT 矩阵对角化：
    $$C = F^* \Lambda F,$$
    其中 $\Lambda = \operatorname{diag}(\hat{c}_0, \hat{c}_1, \ldots, \hat{c}_{n-1})$，$\hat{c}_k = \sum_{j=0}^{n-1} c_j \omega^{jk}$ 是向量 $(c_0, \ldots, c_{n-1})$ 的 DFT。

??? proof "证明"
    **步骤一**：证明循环置换矩阵 $\Pi$ 的特征值和特征向量。

    $\Pi$ 满足 $\Pi^n = I$，故 $\Pi$ 的特征值是 $x^n = 1$ 的根，即 $1, \omega, \omega^2, \ldots, \omega^{n-1}$。对应的特征向量为
    $$\boldsymbol{v}_k = \frac{1}{\sqrt{n}}(1, \omega^k, \omega^{2k}, \ldots, \omega^{(n-1)k})^T, \quad k = 0, 1, \ldots, n-1.$$

    验证：$(\Pi \boldsymbol{v}_k)_j = (\boldsymbol{v}_k)_{j-1 \pmod n} = \frac{1}{\sqrt{n}}\omega^{k(j-1)}$。而 $(\omega^k \boldsymbol{v}_k)_j = \frac{1}{\sqrt{n}}\omega^{k}\omega^{kj} = \frac{1}{\sqrt{n}}\omega^{k(j+1)}$。

    修正下标约定：设 $(\boldsymbol{v}_k)_j = \frac{1}{\sqrt{n}}\omega^{kj}$（$j = 0, 1, \ldots, n-1$），则 $(\Pi\boldsymbol{v}_k)_j = (\boldsymbol{v}_k)_{j-1} = \frac{1}{\sqrt{n}}\omega^{k(j-1)} = \omega^{-k}(\boldsymbol{v}_k)_j$。但 $\Pi$ 将分量下移，故 $(\Pi\boldsymbol{v}_k)_j = (\boldsymbol{v}_k)_{j-1 \bmod n} = \frac{1}{\sqrt{n}}\omega^{k(j-1)} = \omega^{-k}(\boldsymbol{v}_k)_j$。

    因此 $\Pi \boldsymbol{v}_k = \omega^{-k}\boldsymbol{v}_k$，故 $\omega^{-k}$ 是特征值。等价地，$\Pi$ 的特征值为 $\omega^0, \omega^{-1}, \ldots, \omega^{-(n-1)}$，即 $1, \omega^{n-1}, \omega^{n-2}, \ldots, \omega$，这与 $1, \omega, \omega^2, \ldots, \omega^{n-1}$ 是同一组值（只是排列不同）。

    **步骤二**：$C = \sum_{j=0}^{n-1} c_j \Pi^j$ 是 $\Pi$ 的多项式。故 $C$ 与 $\Pi$ 同时对角化，$C$ 的特征向量也是 $\boldsymbol{v}_k$。特征值为
    $$C\boldsymbol{v}_k = \sum_{j=0}^{n-1} c_j \Pi^j \boldsymbol{v}_k = \sum_{j=0}^{n-1} c_j (\omega^{-k})^j \boldsymbol{v}_k = \sum_{j=0}^{n-1} c_j \omega^{-jk} \boldsymbol{v}_k.$$

    定义 $\hat{c}_k = \sum_{j=0}^{n-1} c_j \omega^{jk}$，则 $C\boldsymbol{v}_k = \hat{c}_{-k}\boldsymbol{v}_k = \hat{c}_{n-k}\boldsymbol{v}_k$。

    **步骤三**：将特征向量排成矩阵 $F^* = (\boldsymbol{v}_0, \boldsymbol{v}_1, \ldots, \boldsymbol{v}_{n-1})$，则 $CF^* = F^* \Lambda$，即 $C = F^* \Lambda F$。（此处 $\Lambda$ 的对角元素经过适当重排后恰为 $\hat{c}_0, \hat{c}_1, \ldots, \hat{c}_{n-1}$。） $\blacksquare$

!!! example "例 37A.4"
    循环矩阵-向量乘法 $\boldsymbol{y} = C\boldsymbol{x}$ 的快速算法：

    1. 计算 $\hat{\boldsymbol{x}} = F\boldsymbol{x}$（$\boldsymbol{x}$ 的 DFT），$O(n\log n)$；
    2. 计算 $\hat{\boldsymbol{c}} = F\boldsymbol{c}$（$\boldsymbol{c}$ 的 DFT），$O(n\log n)$（可预计算）；
    3. 逐元素乘 $\hat{y}_k = \hat{c}_k \hat{x}_k$，$O(n)$；
    4. 计算 $\boldsymbol{y} = F^*\hat{\boldsymbol{y}}$（逆 DFT），$O(n\log n)$。

    总计 $O(n\log n)$，而直接乘法为 $O(n^2)$。

    数值示例：$C = \begin{pmatrix} 1 & 3 & 2 \\ 2 & 1 & 3 \\ 3 & 2 & 1 \end{pmatrix}$，$\boldsymbol{x} = (1, 0, 0)^T$。

    $\hat{\boldsymbol{c}} = (1+2+3,\; 1+2\omega+3\omega^2,\; 1+2\omega^2+3\omega^4) = (6,\; 1+2\omega+3\omega^2,\; 1+2\omega^2+3\omega)$，其中 $\omega = e^{2\pi i/3}$。

    $\boldsymbol{y} = C\boldsymbol{x} = (1, 2, 3)^T$（即 $C$ 的第一列）。

---

## 37A.6 循环代数与环同构

!!! theorem "定理 37A.10 (循环矩阵的代数结构)"
    (a) $n$ 阶循环矩阵的全体 $\mathcal{C}_n$ 构成 $M_n(\mathbb{C})$ 的一个 $n$ 维**交换子代数**。

    (b) 两个循环矩阵的乘积是循环矩阵，且乘法可交换。

    (c) 非奇异循环矩阵的逆是循环矩阵。

    (d) 存在环同构
    $$\mathcal{C}_n \cong \mathbb{C}[x]/(x^n - 1).$$

    (e) 通过 DFT，有 $\mathbb{C}$-代数同构
    $$\mathcal{C}_n \cong \mathbb{C}^n$$
    （直积环），其中右侧的乘法是逐分量乘法。

??? proof "证明"
    (a)(b) 由 $C_1 = F^*\Lambda_1 F$，$C_2 = F^*\Lambda_2 F$，得 $C_1 C_2 = F^*\Lambda_1\Lambda_2 F = F^*\Lambda_2\Lambda_1 F = C_2 C_1$。对角矩阵的乘积仍为对角矩阵，故 $C_1C_2$ 为循环矩阵。维数为 $n$，因为循环矩阵由 $n$ 个参数 $c_0, \ldots, c_{n-1}$ 决定。

    (c) $C^{-1} = F^*\Lambda^{-1}F$（当 $\Lambda$ 可逆时），$\Lambda^{-1}$ 仍为对角矩阵。

    (d) 将 $c(x) = c_0 + c_1 x + \cdots + c_{n-1}x^{n-1}$ 对应于循环矩阵 $C = c(\Pi)$。由于 $\Pi^n = I$，$\Pi$ 满足 $x^n - 1 = 0$。映射 $c(x) \mapsto c(\Pi)$ 是满射环同态，核为 $(x^n - 1)$，故由同构定理得 $\mathcal{C}_n \cong \mathbb{C}[x]/(x^n - 1)$。

    (e) 由中国剩余定理，$\mathbb{C}[x]/(x^n-1) \cong \prod_{k=0}^{n-1}\mathbb{C}[x]/(x - \omega^k) \cong \mathbb{C}^n$。这个同构恰好由 DFT 给出：$c(x) \mapsto (c(1), c(\omega), \ldots, c(\omega^{n-1}))$。

!!! example "例 37A.5"
    $n = 4$ 时，$x^4 - 1 = (x-1)(x+1)(x-i)(x+i)$。故
    $$\mathcal{C}_4 \cong \mathbb{C}[x]/(x^4-1) \cong \mathbb{C} \times \mathbb{C} \times \mathbb{C} \times \mathbb{C}.$$
    循环矩阵 $C$ 在此同构下对应于其 4 个 DFT 特征值 $(\hat{c}_0, \hat{c}_1, \hat{c}_2, \hat{c}_3)$。$C$ 可逆当且仅当所有 $\hat{c}_k \ne 0$。

!!! note "注记 37A.5 (实循环矩阵)"
    当限制为实数域时，$\mathbb{R}[x]/(x^n-1)$ 的分解取决于 $x^n-1$ 在 $\mathbb{R}$ 上的因式分解。例如 $n$ 为偶数时，$x^n - 1 = (x-1)(x+1)\prod(x^2 - 2\cos(2\pi k/n)x + 1)$，故 $\mathcal{C}_n(\mathbb{R}) \cong \mathbb{R} \times \mathbb{R} \times \mathbb{C}^{(n-2)/2}$。

---

## 37A.7 分块 Toeplitz 矩阵与 Toeplitz 分块矩阵

!!! definition "定义 37A.9 (分块 Toeplitz 矩阵)"
    设 $T_{ij} \in M_p(\mathbb{C})$ 为 $p \times p$ 矩阵块。**分块 Toeplitz 矩阵**（Block Toeplitz Matrix, BTM）是指
    $$\mathcal{T} = (T_{i-j})_{1 \le i,j \le n} = \begin{pmatrix}
    T_0 & T_{-1} & \cdots & T_{-(n-1)} \\
    T_1 & T_0 & \cdots & T_{-(n-2)} \\
    \vdots & & \ddots & \vdots \\
    T_{n-1} & T_{n-2} & \cdots & T_0
    \end{pmatrix} \in M_{np}(\mathbb{C}),$$
    其中块 $T_{ij}$ 仅取决于差 $i - j$，但各块 $T_k$ 本身可以是任意的 $p \times p$ 矩阵。

!!! definition "定义 37A.10 (Toeplitz 分块矩阵)"
    **Toeplitz 分块矩阵**（Toeplitz Block Matrix, TBM）是指分块矩阵
    $$\mathcal{T} = (T_{ij})_{1 \le i,j \le n},$$
    其中每个块 $T_{ij}$ 本身是 Toeplitz 矩阵，但 $T_{ij}$ 不一定仅取决于 $i-j$。

!!! theorem "定理 37A.11 (分块循环矩阵的对角化)"
    设 $\mathcal{C}$ 为 $n \times n$ 分块循环矩阵（块大小 $p \times p$），即
    $$\mathcal{C} = \sum_{k=0}^{n-1} C_k \otimes E_k,$$
    其中 $E_k$ 是由 $\Pi^k$ 的位置模式确定的排列。则
    $$\mathcal{C} = (F \otimes I_p)^* \operatorname{diag}(\hat{C}_0, \hat{C}_1, \ldots, \hat{C}_{n-1})(F \otimes I_p),$$
    其中 $\hat{C}_k = \sum_{j=0}^{n-1} C_j \omega^{jk}$ 是矩阵值 DFT。

??? proof "证明"
    与标量情形类似，利用 $F \otimes I_p$ 对 $\Pi \otimes I_p$ 的对角化，以及循环矩阵是 $\Pi$ 的多项式。每个块 $\hat{C}_k$ 是 $p \times p$ 矩阵，因此分块循环矩阵的特征值问题归结为 $n$ 个 $p \times p$ 矩阵的特征值问题。

!!! note "注记 37A.6 (多通道信号处理)"
    分块 Toeplitz 矩阵自然出现在多通道（MIMO）信号处理中。若 $\{\boldsymbol{X}_t\}$ 是 $p$ 维宽平稳向量随机过程，其自协方差矩阵 $R_{ij} = \operatorname{Cov}(\boldsymbol{X}_i, \boldsymbol{X}_j) = R(i-j)$ 构成的协方差矩阵具有分块 Toeplitz 结构，其中每个块 $R(k)$ 是 $p \times p$ 矩阵。

---

## 37A.8 Szegő 分布定理

!!! definition "定义 37A.11 (生成函数/符号)"
    设 $T_n = (t_{|i-j|})_{1 \le i,j \le n}$ 为对称 Toeplitz 矩阵族。若序列 $\{t_k\}$ 的 Fourier 级数
    $$f(\theta) = \sum_{k=-\infty}^{\infty} t_k e^{ik\theta} = t_0 + 2\sum_{k=1}^{\infty} t_k \cos(k\theta)$$
    收敛（在 $L^1$ 意义下），则称 $f$ 为 $\{T_n\}$ 的**生成函数**（或**符号**）。对 Hermite Toeplitz 矩阵，$f$ 为实值函数。

!!! theorem "定理 37A.12 (Szegő 等分布定理)"
    设 $f \in L^1([0, 2\pi])$ 为实值函数，$T_n$ 为对应的 $n$ 阶对称 Toeplitz 矩阵。设 $\lambda_1^{(n)} \le \lambda_2^{(n)} \le \cdots \le \lambda_n^{(n)}$ 为 $T_n$ 的特征值。则：

    (a) **特征值界**：$\operatorname{ess\,inf} f \le \lambda_k^{(n)} \le \operatorname{ess\,sup} f$ 对所有 $k, n$。

    (b) **等分布定理**（Szegő, 1920）：对任意连续函数 $F$，
    $$\lim_{n\to\infty} \frac{1}{n}\sum_{k=1}^{n} F(\lambda_k^{(n)}) = \frac{1}{2\pi}\int_0^{2\pi} F(f(\theta))\, d\theta.$$

    直观地说，$T_n$ 的特征值在 $n \to \infty$ 时"等分布"于 $f$ 的值域上，密度由 $f$ 的取值分布决定。

??? proof "证明"
    **(a)** 由 Rayleigh 商，$\lambda_{\min}(T_n) = \min_{\|\boldsymbol{x}\|=1} \boldsymbol{x}^T T_n \boldsymbol{x}$。利用 Toeplitz 矩阵的二次型表达：
    $$\boldsymbol{x}^T T_n \boldsymbol{x} = \frac{1}{2\pi}\int_0^{2\pi} f(\theta) |p_n(\theta)|^2\, d\theta,$$
    其中 $p_n(\theta) = \sum_{k=1}^n x_k e^{ik\theta}$ 为三角多项式，$\|p_n\|_{L^2}^2 = \|\boldsymbol{x}\|^2$。由此 $\boldsymbol{x}^T T_n \boldsymbol{x} \ge (\operatorname{ess\,inf} f) \|\boldsymbol{x}\|^2$，上界类似。

    **(b)** 证明分为三步：

    **第一步**：对 $F(\lambda) = \lambda^m$（$m = 0, 1, 2, \ldots$），
    $$\frac{1}{n}\sum_{k=1}^n (\lambda_k^{(n)})^m = \frac{1}{n}\operatorname{tr}(T_n^m).$$
    利用 $T_n^m$ 的迹可以用 Fourier 系数的卷积表达，当 $n \to \infty$ 时，
    $$\frac{1}{n}\operatorname{tr}(T_n^m) \to \frac{1}{2\pi}\int_0^{2\pi} f(\theta)^m\, d\theta.$$
    这需要仔细估计 $T_n^m$ 的迹中的边界效应项（即非 Toeplitz 部分的贡献），证明它们是 $o(n)$。

    **第二步**：由上一步，等式对所有多项式 $F$ 成立。

    **第三步**：由 Stone-Weierstrass 定理，多项式在 $[m_f, M_f]$（$m_f = \operatorname{ess\,inf} f$，$M_f = \operatorname{ess\,sup} f$）上的连续函数空间中稠密。由 (a) 的特征值界和一致逼近，推广到一般连续函数 $F$。 $\blacksquare$

!!! example "例 37A.6"
    设 $f(\theta) = 2 - 2\cos\theta$，对应 $t_0 = 2, t_{\pm 1} = -1, t_k = 0$（$|k| \ge 2$）。
    则 $T_n$ 为三对角 Toeplitz 矩阵
    $$T_n = \begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & \ddots & \ddots & \ddots \\ & & -1 & 2 \end{pmatrix}.$$
    这正是一维离散 Laplacian。其特征值精确已知：
    $$\lambda_k^{(n)} = 2 - 2\cos\frac{k\pi}{n+1} = 4\sin^2\frac{k\pi}{2(n+1)}, \quad k = 1, \ldots, n.$$
    $f(\theta) = 2 - 2\cos\theta$ 的值域为 $[0, 4]$，特征值确实分布在 $[0, 4]$ 中，与 Szegő 定理一致。

    取 $F(\lambda) = \lambda$，Szegő 定理给出 $\frac{1}{n}\sum_k \lambda_k^{(n)} \to \frac{1}{2\pi}\int_0^{2\pi}(2-2\cos\theta)\,d\theta = 2$。实际上 $\frac{1}{n}\operatorname{tr}(T_n) = 2$ 对所有 $n$ 精确成立。

---

## 37A.9 强 Szegő 极限定理与 Fisher-Hartwig 猜想

!!! theorem "定理 37A.13 (强 Szegő 极限定理)"
    设 $f > 0$ 且 $\log f \in L^1$，且 $f$ 满足额外的光滑性条件（$\log f$ 的 Fourier 系数 $\hat{c}_k$ 满足 $\sum_{k=1}^{\infty} k|\hat{c}_k|^2 < \infty$）。则
    $$\lim_{n\to\infty} \frac{\det(T_n)}{G^n} = E,$$
    其中 $G = \exp\left(\frac{1}{2\pi}\int_0^{2\pi} \log f(\theta)\, d\theta\right)$ 为**几何平均**，$E$ 为可以用 $f$ 的 Fourier 系数精确表示的常数：
    $$E = \exp\left(\sum_{k=1}^{\infty} k |\hat{c}_k|^2\right),$$
    其中 $\hat{c}_k$ 是 $\log f$ 的 Fourier 系数。

??? proof "证明"
    **概要**：

    **第一步（弱 Szegő 定理）**：首先证明 $\lim_{n\to\infty} (\det T_n)^{1/n} = G$。这由 Szegő 等分布定理取 $F(\lambda) = \log\lambda$ 得到：
    $$\frac{1}{n}\log\det(T_n) = \frac{1}{n}\sum_{k=1}^n \log\lambda_k^{(n)} \to \frac{1}{2\pi}\int_0^{2\pi}\log f(\theta)\,d\theta = \log G.$$

    **第二步**：为得到更精细的渐近 $\det(T_n) \sim E \cdot G^n$，需要分析
    $$\log\det(T_n) - n\log G = \sum_{k=1}^n\log\lambda_k^{(n)} - n\log G.$$
    利用 $\det(T_n) = \prod_{k=0}^{n-1}\sigma_k^2$（其中 $\sigma_k^2$ 来自 Levinson 递推）和谱分解理论，可以证明此差趋于 $\sum_{k=1}^{\infty}k|\hat{c}_k|^2$。

    完整证明需要用到 Wiener-Hopf 因式分解和 Toeplitz 算子理论。

!!! theorem "定理 37A.14 (Fisher-Hartwig 猜想)"
    当生成函数 $f$ 具有**零点**或**跳跃间断点**时，Szegő 极限定理的渐近形式需要修正。设
    $$f(\theta) = f_0(\theta) \prod_{r=1}^R |e^{i\theta} - e^{i\theta_r}|^{2\alpha_r} e^{i\beta_r(\theta - \theta_r - \pi)},$$
    其中 $f_0$ 是光滑的正函数，$\alpha_r > -1/2$ 控制零点/奇点的阶，$\beta_r$ 控制跳跃。Fisher-Hartwig 猜想（1968，现已大部分被证明）断言
    $$\det(T_n) \sim E' \cdot G^n \cdot n^{\sum_r(\alpha_r^2 - \beta_r^2)} \cdot \prod_{r<s}|e^{i\theta_r}-e^{i\theta_s}|^{-2(\alpha_r\alpha_s - \beta_r\beta_s)},$$
    其中 $E'$ 涉及 Barnes G-函数。

!!! note "注记 37A.7 (Fisher-Hartwig 与统计力学)"
    Fisher-Hartwig 猜想在统计力学中有重要应用。二维 Ising 模型在临界温度下的自发磁化可以表示为 Toeplitz 行列式，其生成函数恰好具有 Fisher-Hartwig 型奇点。该猜想的证明（Basor, Widom, Deift-Its-Krasovsky 等人的工作）联系了随机矩阵理论、可积系统和 Riemann-Hilbert 问题。

---

## 37A.10 循环预处理

!!! definition "定义 37A.12 (最优循环逼近)"
    给定 $n$ 阶 Toeplitz 矩阵 $T$，**T. Chan 循环预处理子** $c_F(T)$ 定义为在 Frobenius 范数意义下最接近 $T$ 的循环矩阵：
    $$c_F(T) = \arg\min_{C \in \mathcal{C}_n} \|T - C\|_F.$$
    $c_F(T)$ 的第一列各分量可以显式给出：对 $k = 0, 1, \ldots, n-1$，
    $$c_k = \frac{(n-k)t_k + k t_{k-n}}{n}.$$

!!! theorem "定理 37A.15 (循环预处理的谱聚类)"
    设 $f$ 为正的连续生成函数。则预处理矩阵 $c_F(T_n)^{-1}T_n$ 的特征值在 $n \to \infty$ 时聚集（cluster）在 $1$ 附近：对任意 $\epsilon > 0$，
    $$\#\{k : |\lambda_k(c_F(T_n)^{-1}T_n) - 1| > \epsilon\} = o(n).$$

??? proof "证明"
    **思路**：$c_F(T_n)$ 的特征值为 $f$ 在等距节点 $\theta_k = 2\pi k/n$ 上的采样值，$T_n$ 的特征值渐近分布于 $f$ 的值域。因此
    $$c_F(T_n)^{-1}T_n \approx \operatorname{diag}(f(\theta_k))^{-1}\operatorname{diag}(f(\theta_k)) \approx I.$$
    严格地，$c_F(T_n) - T_n$ 的秩有界（与 $n$ 无关），由 Szegő 定理的推论得到谱聚类。

!!! note "注记 37A.8 (预处理共轭梯度法)"
    利用循环预处理子 $C = c_F(T)$，预处理共轭梯度法（PCG）求解 $Tx = b$ 的每步迭代需要：(1) 一次 Toeplitz 矩阵-向量乘（通过嵌入 $2n$ 阶循环矩阵，用 FFT 在 $O(n\log n)$ 完成）；(2) 一次循环预处理子的求解（$O(n\log n)$）。由谱聚类性质，迭代步数与 $n$ 无关（对光滑生成函数），故总复杂度为 $O(n\log n)$。这是 T. Chan 和 R. Chan 等人在 20 世纪 80-90 年代发展的重要方法。

    **Strang 循环预处理子**是另一种常用选择：取 $T$ 的中心 $n/2$ 条对角线，然后周期延拓得到循环矩阵。

---

## 37A.11 椭球波函数与结构化矩阵的完成问题

!!! definition "定义 37A.13 (椭球波函数)"
    **椭球球面波函数**（Prolate Spheroidal Wave Functions, PSWFs）是同时在时域和频域上近似集中的函数。它们是以下积分方程的特征函数：
    $$\int_{-T}^{T} \frac{\sin \Omega(t-s)}{\pi(t-s)} \psi_k(s)\, ds = \lambda_k \psi_k(t),$$
    其中 $T$ 是时间半宽，$\Omega$ 是频率半宽。

!!! theorem "定理 37A.16 (PSWFs 与 Toeplitz 矩阵的联系)"
    将上述积分方程离散化，核 $K(t,s) = \frac{\sin\Omega(t-s)}{\pi(t-s)}$ 在等距网格上的采样产生 Toeplitz 矩阵。该矩阵的特征值表现出**阶梯效应**（plunge phenomenon）：前约 $2\Omega T/\pi$ 个特征值接近 $1$，其余迅速衰减到 $0$。

    更精确地，设 $c = \Omega T$，则对 Slepian 矩阵 $B_n$（$n$ 阶离散化），
    $$\lambda_k \approx 1 \text{ 当 } k \le \lfloor 2c/\pi \rfloor, \quad \lambda_k \approx 0 \text{ 当 } k \gg 2c/\pi.$$
    过渡区间的宽度为 $O(\log n)$。

!!! note "注记 37A.9 (Slepian 的贡献)"
    Slepian, Landau 和 Pollak 在 Bell Labs 的经典系列论文（1961-1978）建立了 PSWFs 的系统理论。它们在信号处理中的意义在于量化了**时频集中度**的基本限制——一个函数不可能同时在时域和频域上任意集中（一种不确定性原理的定量形式）。PSWFs 在近年来的压缩感知和数值分析中重新焕发了活力。

!!! definition "定义 37A.14 (Toeplitz 矩阵的完成问题)"
    **正定 Toeplitz 完成问题**：给定部分 Toeplitz 数据 $t_0, t_1, \ldots, t_m$（$m < n-1$），是否存在 $t_{m+1}, \ldots, t_{n-1}$ 使得完整的 $n$ 阶 Toeplitz 矩阵 $T_n = (t_{|i-j|})$ 正定？

!!! theorem "定理 37A.17 (正定 Toeplitz 完成)"
    设 $t_0, t_1, \ldots, t_m$ 使得带状 Toeplitz 矩阵 $(t_{|i-j|})_{|i-j|\le m}$ 的所有 $(m+1) \times (m+1)$ 主子矩阵正定。则：

    (a) 存在正定的 $n$ 阶 Toeplitz 完成（对任意 $n > m$）。

    (b) 在所有正定完成中，存在唯一的**最大熵完成**（maximum entropy completion），其生成函数 $f$ 使得 $\int_0^{2\pi}\log f(\theta)\, d\theta$ 最大化。

    (c) 最大熵完成对应的生成函数 $f$ 是次数至多为 $m$ 的三角多项式的倒数。

??? proof "证明"
    **(a)** 由 Levinson 递推的反射系数解释：已知 $t_0, \ldots, t_m$ 可以确定 $\alpha_1, \ldots, \alpha_m$。只要 $|\alpha_k| < 1$（由正定性保证），可以任意选择 $|\alpha_{m+1}| < 1, \ldots, |\alpha_{n-1}| < 1$ 来扩展，对应的 $T_n$ 正定。

    **(b)(c)** 最大熵完成对应于选择 $\alpha_k = 0$（$k > m$），此时预测误差功率不再下降。对应的谱为 AR($m$) 模型的功率谱 $f(\theta) = \sigma_m^2/|1 + a_1^{(m)}e^{i\theta} + \cdots + a_m^{(m)}e^{im\theta}|^2$。

---

## 练习题

1. **[结构] 证明 $n$ 阶 Toeplitz 矩阵 $T$ 可以嵌入到一个 $2n$ 阶循环矩阵 $C$ 中。利用此性质说明如何利用 FFT 在 $O(n\log n)$ 时间内计算 $T\mathbf{x}$。**
   ??? success "参考答案"
       构造 $2n$ 阶循环矩阵 $C$，其第一列由 $T$ 的第一列、一个 0 以及 $T$ 的第一行逆序排列组成。将 $\mathbf{x}$ 补零至 $2n$ 维得到 $\tilde{\mathbf{x}}$，则 $C\tilde{\mathbf{x}}$ 的前 $n$ 个分量即为 $T\mathbf{x}$。循环矩阵乘法通过 FFT 加速。

2. **[对角化] 设 $C = \operatorname{circ}(c_0, c_1, \dots, c_{n-1})$。证明 $C$ 的特征值 $\lambda_k$ 恰好是向量 $\mathbf{c}$ 的离散傅里叶变换（DFT）。**
   ??? success "参考答案"
       定义 $\omega = e^{2\pi i/n}$。令 $\mathbf{v}_k = (1, \omega^k, \dots, \omega^{k(n-1)})^T$。直接验证 $C\mathbf{v}_k = (\sum c_j \omega^{jk}) \mathbf{v}_k$。特征向量 $\mathbf{v}_k$ 与 $\mathbf{c}$ 无关，这一性质使得循环矩阵在工程中极具价值。

3. **[Levinson] 在 Levinson-Durbin 递推中，如果某个反射系数 $|\alpha_k| = 1$，说明了矩阵 $T$ 的什么性质？**
   ??? success "参考答案"
       说明 $T_{k+1}$ 是奇异的。如果 $|\alpha_k| < 1$ 对所有 $k$ 成立，则 $T$ 是正定的。这提供了判定 Toeplitz 矩阵正定性的高效方法。

4. **[Hankel] 证明：如果 $H$ 是 Hankel 矩阵，则 $J H$ 是 Toeplitz 矩阵（其中 $J$ 是反对角单位阵）。**
   ??? success "参考答案"
       $(JH)_{ij} = H_{n-i+1, j}$。由于 $H$ 的元素取决于指标之和，即 $(n-i+1) + j = n + 1 + (j-i)$。由于此值仅取决于 $j-i$，故 $JH$ 符合 Toeplitz 矩阵的定义。

5. **[Szegő] 根据 Szegő 等分布定理，三对角矩阵 $\operatorname{tridiag}(-1, 2, -1)$ 在维数 $n \to \infty$ 时的特征值分布倾向于什么函数？**
   ??? success "参考答案"
       生成函数为 $f(\theta) = 2 - e^{i\theta} - e^{-i\theta} = 2 - 2\cos\theta$。特征值在 $[0, 4]$ 区间内按此余弦函数的形状分布。

6. **[逆矩阵] 循环矩阵的逆矩阵（若存在）一定还是循环矩阵吗？为什么？**
   ??? success "参考答案"
       是。因为循环矩阵全体在 DFT 变换下对应于对角矩阵。对角矩阵的逆（逐元素取倒数）仍是对角矩阵，逆变换回原空间后仍是循环矩阵。

7. **[卷积] 说明循环矩阵与向量的乘法 $C\mathbf{x}$ 实际上对应于序列的什么运算？**
   ??? success "参考答案"
       对应于向量 $\mathbf{c}$ 与 $\mathbf{x}$ 的**循环卷积**（Circular Convolution）。这解释了为什么 FFT 可以加速此类运算。

8. **[行列式] 利用特征值公式，求全 1 矩阵 $J_n$（它是特殊的循环矩阵）的行列式。**
   ??? success "参考答案"
       $J_n$ 的特征值为 $\hat{c}_0 = n$ 且 $\hat{c}_k = 0$ ($k \ge 1$)。由于存在零特征值，对于 $n > 1$，$\det(J_n) = 0$。

9. **[预处理] 为什么在求解大型 Toeplitz 方程组时，常用循环矩阵作为预处理子（Preconditioner）？**
   ??? success "参考答案"
       因为循环矩阵的逆极易计算（通过 FFT），且对于很多平稳过程产生的 Toeplitz 矩阵，循环矩阵能很好地逼近其谱特性，从而将特征值聚类在 1 附近，加速共轭梯度法的收敛。

10. **[爱因斯坦思考题] 爱因斯坦的广义相对论强调物理规律在坐标变换下的不变性。循环矩阵的“基无关性”（对所有循环矩阵，DFT 矩阵都是通用的特征向量基）反映了离散空间的哪种对称性？**
    ??? success "参考答案"
        反映了**离散平移对称性**（Cyclic Shift Invariance）。DFT 矩阵本质上是离散平移群的表示。这意味着在具有周期性边界条件的均匀宇宙中，无论我们站在哪个格点上，空间的动力学结构（特征模态）看起来都是完全一样的。循环矩阵是这种“空间均匀性”在代数上的最直接体现。

## 本章小结

本章研究了具有高度参数冗余的结构化矩阵，它们是连接代数与分析、信号处理的纽带：

1. **Toeplitz 与 Hankel 结构**：定义了沿对角线或反对角线恒定的矩阵，揭示了它们在描述平稳过程和矩问题中的天然优势。
2. **Levinson-Durbin 递推**：展示了如何利用结构特性将计算复杂度从 $O(n^3)$ 优化至 $O(n^2)$，这是结构化矩阵算法的典范。
3. **循环矩阵与傅里叶分析**：证明了循环矩阵可被 DFT 矩阵统一对角化，确立了卷积、滤波与线性代数之间的深刻同构关系。
4. **渐近谱理论**：通过 Szegő 定理描述了大型 Toeplitz 矩阵特征值的分布规律，将矩阵的局部性质与生成函数的全局分析联系起来。
5. **快速算法与应用**：讨论了循环预处理在 PCG 算法中的应用，以及这些结构在现代统计、预测理论和数值模拟中的基石作用。

