# 第 22 章 数值线性代数

<div class="context-flow" markdown>

**前置**：LU/QR/SVD/特征值分解(Ch5-8)

**脉络**：浮点误差是根源 → 条件数决定精度 → 直接法(LU, $O(n^3)$) vs 迭代法(Krylov, $O(n \cdot \text{nnz})$) → 稀疏结构决定算法选择

**延伸**：数值线性代数是所有大规模科学计算的引擎：有限元分析（结构工程）、计算流体力学、量子化学（Hartree-Fock 方程）、天气预报模型都依赖高效的数值线性代数算法

</div>

数值线性代数（numerical linear algebra）是科学计算的核心。在实际应用中，矩阵运算不可避免地受到浮点舍入误差的影响，算法的稳定性和效率成为关键问题。本章从浮点运算基础出发，系统介绍线性方程组的直接法与迭代法、Krylov 子空间方法、特征值的数值计算以及稀疏矩阵技术。这些方法构成了现代科学计算和工程模拟的基石。

## 22.1 浮点运算与舍入误差

<div class="context-flow" markdown>

**起点**：实数 → 有限位浮点数 · 机器精度 $\epsilon_{\text{mach}} \approx 10^{-16}$（双精度）是一切误差的原子单位 → **灾难性消去**是最常见陷阱

</div>

### IEEE 754 浮点数

!!! definition "定义 22.1 (浮点数系统)"
    一个**浮点数系统**（floating-point number system）$\mathbb{F}(\beta, t, L, U)$ 由以下参数确定：

    - $\beta$：基数（base），通常 $\beta = 2$；
    - $t$：尾数（mantissa）位数；
    - $[L, U]$：指数范围。

    该系统中的浮点数表示为

    $$x = \pm d_0.d_1 d_2 \cdots d_{t-1} \times \beta^e = \pm \beta^e \sum_{i=0}^{t-1} d_i \beta^{-i},$$

    其中 $0 \leq d_i \leq \beta - 1$，$d_0 \neq 0$（规格化），$L \leq e \leq U$。

    IEEE 754 标准定义了两种常用格式：

    | 格式 | 位数 | 尾数位 ($t$) | 指数范围 |
    |------|------|------------|---------|
    | 单精度 (float32) | 32 | 24 | $[-126, 127]$ |
    | 双精度 (float64) | 64 | 53 | $[-1022, 1023]$ |

!!! definition "定义 22.2 (机器精度)"
    **机器精度**（machine epsilon）$\epsilon_{\text{mach}}$ 定义为最小的正浮点数 $\epsilon$ 使得 $\text{fl}(1 + \epsilon) > 1$，即浮点数系统能区分 $1$ 和 $1 + \epsilon$ 的最小 $\epsilon$。等价地，

    $$\epsilon_{\text{mach}} = \frac{1}{2} \beta^{1-t}.$$

    对于 IEEE 双精度，$\epsilon_{\text{mach}} \approx 1.11 \times 10^{-16}$。

    任何实数 $x$ 的浮点表示 $\text{fl}(x)$ 满足

    $$\text{fl}(x) = x(1 + \delta), \quad |\delta| \leq \epsilon_{\text{mach}}.$$

!!! theorem "定理 22.1 (浮点运算的基本误差界)"
    设 $\oplus$ 表示浮点加法，$\ominus$ 表示浮点减法，$\otimes$ 表示浮点乘法，$\oslash$ 表示浮点除法。对任意浮点数 $a, b$，有

    $$a \circledast b = (a * b)(1 + \delta), \quad |\delta| \leq \epsilon_{\text{mach}},$$

    其中 $\circledast \in \{\oplus, \ominus, \otimes, \oslash\}$，$* \in \{+, -, \times, /\}$ 是对应的精确运算。

??? proof "证明"
    这是 IEEE 754 标准的基本保证。标准规定每一步浮点运算的结果必须是精确结果的正确舍入（round to nearest），即 $\text{fl}(a * b)$ 是最接近 $a * b$ 的浮点数。由此直接得到相对误差界 $\epsilon_{\text{mach}}$。

!!! example "例 22.1"
    **灾难性消去**（catastrophic cancellation）。考虑计算 $f(x) = \sqrt{x+1} - \sqrt{x}$，当 $x$ 很大时。

    取 $x = 10^{16}$。精确值约为 $f(x) \approx \frac{1}{2\sqrt{x}} \approx 5 \times 10^{-9}$。

    但在双精度浮点下，$\text{fl}(\sqrt{10^{16}+1}) = \text{fl}(\sqrt{x})$（因为 $10^{16}+1$ 与 $10^{16}$ 在双精度下可能是同一个浮点数），导致 $\text{fl}(f(x)) = 0$，相对误差为 $100\%$。

    改用数值稳定的等价形式 $f(x) = \frac{1}{\sqrt{x+1} + \sqrt{x}}$ 可避免此问题。

## 22.2 数值稳定性

<div class="context-flow" markdown>

**核心框架**：前向误差(输出偏差) ≤ **条件数** × 后向误差(等价输入扰动) → 后向稳定算法 + 良态问题 = 精确结果

**关键**：$\kappa(A) = \sigma_{\max}/\sigma_{\min}$（SVD, Ch8）——条件数是问题固有的，算法只控制后向误差

</div>

!!! definition "定义 22.3 (前向误差与后向误差)"
    设算法 $\hat{f}$ 是函数 $f$ 的数值近似。对输入 $x$：

    - **前向误差**（forward error）：$\|\hat{f}(x) - f(x)\|$，度量输出的误差；
    - **后向误差**（backward error）：满足 $\hat{f}(x) = f(x + \Delta x)$ 的最小 $\|\Delta x\|$，即将输出误差解释为输入扰动。

!!! definition "定义 22.4 (条件数)"
    问题 $f$ 的**条件数**（condition number）定义为

    $$\kappa = \lim_{\delta \to 0} \sup_{\|\Delta x\| \leq \delta} \frac{\|f(x + \Delta x) - f(x)\| / \|f(x)\|}{\|\Delta x\| / \|x\|}.$$

    条件数度量问题本身对输入扰动的敏感程度。条件数大的问题称为**病态的**（ill-conditioned）。

    对于线性方程组 $A\mathbf{x} = \mathbf{b}$，条件数为

    $$\kappa(A) = \|A\| \cdot \|A^{-1}\|.$$

    使用 $2$-范数时，$\kappa_2(A) = \sigma_{\max}(A) / \sigma_{\min}(A)$，其中 $\sigma_{\max}$、$\sigma_{\min}$ 分别是最大和最小奇异值。

!!! definition "定义 22.5 (数值稳定性)"
    - **前向稳定**（forward stable）：算法的前向误差很小，即 $\|\hat{f}(x) - f(x)\| / \|f(x)\| = O(\epsilon_{\text{mach}})$。
    - **后向稳定**（backward stable）：输出可以解释为略微扰动的输入的精确输出，即 $\hat{f}(x) = f(x + \Delta x)$，其中 $\|\Delta x\| / \|x\| = O(\epsilon_{\text{mach}})$。
    - **混合稳定**（mixed stable）：$\hat{f}(x) + \Delta y = f(x + \Delta x)$，其中 $\|\Delta x\|/\|x\|$ 和 $\|\Delta y\|/\|f(x)\|$ 都是 $O(\epsilon_{\text{mach}})$。

!!! theorem "定理 22.2 (后向稳定性与前向误差的关系)"
    若算法 $\hat{f}$ 是后向稳定的，则前向相对误差满足

    $$\frac{\|\hat{f}(x) - f(x)\|}{\|f(x)\|} \leq \kappa(f, x) \cdot O(\epsilon_{\text{mach}}).$$

    即后向稳定的算法对于良态问题（$\kappa$ 小）给出精确的结果。

??? proof "证明"
    由后向稳定性，$\hat{f}(x) = f(x + \Delta x)$，$\|\Delta x\|/\|x\| = O(\epsilon_{\text{mach}})$。由条件数的定义，

    $$\frac{\|f(x + \Delta x) - f(x)\|}{\|f(x)\|} \leq \kappa(f, x) \cdot \frac{\|\Delta x\|}{\|x\|} = \kappa(f, x) \cdot O(\epsilon_{\text{mach}}).$$

!!! example "例 22.2"
    **矩阵条件数的例子**。Hilbert 矩阵 $H_n$ 的元素为 $H_{ij} = 1/(i+j-1)$，是著名的病态矩阵：

    | $n$ | $\kappa_2(H_n)$ |
    |-----|----------------|
    | 5   | $4.8 \times 10^5$ |
    | 10  | $1.6 \times 10^{13}$ |
    | 15  | $3.7 \times 10^{17}$ |

    当 $n \geq 13$ 时，$\kappa_2(H_n) > 1/\epsilon_{\text{mach}}$，此时用双精度浮点求解 $H_n \mathbf{x} = \mathbf{b}$ 的结果可能完全不可靠。

## 22.3 直接法求解线性方程组

<div class="context-flow" markdown>

**直接法**：$PA = LU$（$\frac{2}{3}n^3$）→ 部分选主元保证 $|l_{ij}| \le 1$ → 后向稳定 · 增长因子 $g(n)$ 理论可达 $2^{n-1}$，实际 $O(n)$

</div>

### 带部分选主元的高斯消元

!!! definition "定义 22.6 (部分选主元)"
    在高斯消元（Gaussian elimination）的第 $k$ 步，**部分选主元**（partial pivoting）策略是在第 $k$ 列的第 $k$ 行及其以下的元素中选取绝对值最大的作为主元，然后交换行使之成为主元行。

    具体地，选取 $p = \arg\max_{k \leq i \leq n} |a^{(k)}_{ik}|$，然后交换第 $k$ 行与第 $p$ 行。

!!! theorem "定理 22.3 (LU 分解)"
    设 $A$ 是 $n \times n$ 非奇异矩阵。带部分选主元的高斯消元产生分解

    $$PA = LU,$$

    其中 $P$ 是置换矩阵，$L$ 是单位下三角矩阵（对角线元素为 $1$，且 $|l_{ij}| \leq 1$），$U$ 是上三角矩阵。

    计算量为 $\frac{2}{3}n^3 + O(n^2)$ 次浮点运算。

??? proof "证明"
    高斯消元的第 $k$ 步（$k = 1, \ldots, n-1$）：

    1. 选主元并交换行（编码为置换矩阵 $P_k$）；
    2. 对第 $i$ 行（$i > k$），计算乘数 $l_{ik} = a^{(k)}_{ik} / a^{(k)}_{kk}$，然后执行行变换 $R_i \leftarrow R_i - l_{ik} R_k$。

    由于选主元保证 $|a^{(k)}_{kk}| \geq |a^{(k)}_{ik}|$（$\forall i > k$），故 $|l_{ik}| \leq 1$。

    消元过程等价于 $M_{n-1} P_{n-1} \cdots M_1 P_1 A = U$，其中 $M_k$ 是消元矩阵。重新排列可得 $PA = LU$。

    运算量：第 $k$ 步需要 $2(n-k)^2$ 次浮点运算（含乘法和加法），总和 $\sum_{k=1}^{n-1} 2(n-k)^2 = \frac{2}{3}n^3 + O(n^2)$。

### 误差分析

!!! theorem "定理 22.4 (LU 分解的后向误差)"
    带部分选主元的 LU 分解是后向稳定的：计算得到的解 $\hat{\mathbf{x}}$ 满足

    $$(A + \Delta A)\hat{\mathbf{x}} = \mathbf{b}, \quad \|\Delta A\| \leq c \, n \, g(n) \, \epsilon_{\text{mach}} \|A\|,$$

    其中 $c$ 是一个适度常数，$g(n)$ 是**增长因子**（growth factor）：

    $$g(n) = \frac{\max_{i,j,k} |a^{(k)}_{ij}|}{\max_{i,j} |a_{ij}|}.$$

    理论上 $g(n)$ 可达 $2^{n-1}$，但在实际计算中几乎总是 $g(n) = O(n)$。

??? proof "证明"
    这是经典的 Wilkinson 后向误差分析的结果。每步消元引入的舍入误差可以等价地看成对原矩阵的扰动。带部分选主元保证 $|l_{ij}| \leq 1$，从而误差不会过度放大。增长因子 $g(n)$ 控制消元过程中矩阵元素的增长。详细的逐步分析参见 Higham (2002)。

!!! example "例 22.3"
    用带部分选主元的高斯消元求解：

    $$\begin{pmatrix} 0.001 & 1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}.$$

    精确解为 $x_1 \approx 1.001$，$x_2 \approx 0.999$。

    **不选主元**：以 $0.001$ 为主元，乘数 $l_{21} = 1000$。在有限精度下，中间步骤的大乘数导致严重的舍入误差。

    **部分选主元**：交换两行后以 $1$ 为主元，乘数 $l_{21} = 0.001$，计算稳定。

## 22.4 迭代法求解线性方程组

<div class="context-flow" markdown>

**转折**：直接法 $O(n^3)$ 对大稀疏系统不可行 → 迭代法只需矩阵-向量乘 → 收敛性由**谱半径** $\rho(G) < 1$ 决定（特征值理论, Ch6）

**层次**：Jacobi/GS（经典）→ SOR（加速）→ Krylov（最优多项式近似, §22.5）

</div>

对于大型稀疏系统，直接法的 $O(n^3)$ 复杂度过高，迭代法成为更好的选择。

### 基本迭代法

!!! definition "定义 22.7 (矩阵分裂)"
    设 $A = M - N$，其中 $M$ 可逆，则线性方程组 $A\mathbf{x} = \mathbf{b}$ 等价于

    $$\mathbf{x} = M^{-1}N\mathbf{x} + M^{-1}\mathbf{b}.$$

    由此产生迭代格式 $\mathbf{x}^{(k+1)} = M^{-1}N\mathbf{x}^{(k)} + M^{-1}\mathbf{b}$。矩阵 $G = M^{-1}N$ 称为**迭代矩阵**（iteration matrix）。

    常用的分裂方式：将 $A = D - L - U$（$D$ 为对角部分，$-L$ 为严格下三角部分，$-U$ 为严格上三角部分），则：

    | 方法 | $M$ | 迭代矩阵 $G$ |
    |------|-----|-------------|
    | **Jacobi 迭代** | $D$ | $D^{-1}(L + U)$ |
    | **Gauss-Seidel 迭代** | $D - L$ | $(D - L)^{-1}U$ |
    | **SOR（$\omega$）** | $\frac{1}{\omega}D - L$ | $(\frac{1}{\omega}D - L)^{-1}[(\frac{1}{\omega}-1)D + U]$ |

!!! theorem "定理 22.5 (迭代法收敛的充要条件)"
    迭代 $\mathbf{x}^{(k+1)} = G\mathbf{x}^{(k)} + \mathbf{c}$ 对任意初始向量 $\mathbf{x}^{(0)}$ 收敛到不动点的充要条件是

    $$\rho(G) < 1,$$

    其中 $\rho(G) = \max_i |\lambda_i(G)|$ 是 $G$ 的**谱半径**（spectral radius）。

??? proof "证明"
    设不动点为 $\mathbf{x}^* = G\mathbf{x}^* + \mathbf{c}$。令误差 $\mathbf{e}^{(k)} = \mathbf{x}^{(k)} - \mathbf{x}^*$，则 $\mathbf{e}^{(k+1)} = G\mathbf{e}^{(k)}$，故 $\mathbf{e}^{(k)} = G^k \mathbf{e}^{(0)}$。

    $G^k \to 0$（$k \to \infty$）对所有 $\mathbf{e}^{(0)}$ 成立当且仅当 $\rho(G) < 1$。

    **必要性**：若 $\rho(G) \geq 1$，取 $\mathbf{e}^{(0)}$ 为对应于模最大特征值 $\lambda$ 的特征向量，则 $\|G^k \mathbf{e}^{(0)}\| = |\lambda|^k \|\mathbf{e}^{(0)}\| \not\to 0$。

    **充分性**：若 $\rho(G) < 1$，由 Gelfand 公式 $\rho(G) = \lim_{k \to \infty} \|G^k\|^{1/k}$，对任意 $\rho(G) < r < 1$，存在 $C > 0$ 使得 $\|G^k\| \leq C r^k$，故 $\|G^k\| \to 0$。

!!! theorem "定理 22.6 (对角占优矩阵的收敛性)"
    若 $A$ 是**严格对角占优**的（strictly diagonally dominant），即

    $$|a_{ii}| > \sum_{j \neq i} |a_{ij}|, \quad \forall i,$$

    则 Jacobi 迭代和 Gauss-Seidel 迭代均收敛。

??? proof "证明"
    以 Jacobi 迭代为例。设 $G_J = D^{-1}(L + U)$。由 Gershgorin 圆盘定理，$G_J$ 的每个特征值 $\lambda$ 满足

    $$|\lambda - 0| \leq \sum_{j \neq i} \left|\frac{a_{ij}}{a_{ii}}\right| < 1$$

    （利用严格对角占优条件）。故 $\rho(G_J) < 1$。

    Gauss-Seidel 迭代的收敛证明更为技术性，可利用 Stein-Rosenberg 定理或直接构造 Lyapunov 函数。

!!! example "例 22.4"
    考虑方程组 $A\mathbf{x} = \mathbf{b}$，其中

    $$A = \begin{pmatrix} 4 & -1 & 0 \\ -1 & 4 & -1 \\ 0 & -1 & 4 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 1 \\ 5 \\ 0 \end{pmatrix}.$$

    $A$ 严格对角占优。Jacobi 迭代矩阵为

    $$G_J = \begin{pmatrix} 0 & 1/4 & 0 \\ 1/4 & 0 & 1/4 \\ 0 & 1/4 & 0 \end{pmatrix}, \quad \rho(G_J) = \frac{1}{2\sqrt{2}} \approx 0.354.$$

    Gauss-Seidel 迭代矩阵的谱半径 $\rho(G_{GS}) = \rho(G_J)^2 \approx 0.125$（对称正定情形下的 Stein-Rosenberg 关系）。因此 Gauss-Seidel 收敛更快。

## 22.5 Krylov 子空间方法

<div class="context-flow" markdown>

**核心思想**：$\mathcal{K}_k(A, \mathbf{b}) = \text{span}\{\mathbf{b}, A\mathbf{b}, \ldots, A^{k-1}\mathbf{b}\}$ = 用 $k$ 次矩阵-向量乘构建的最大信息空间 → Arnoldi 正交化 → Lanczos（对称时三对角化）

**统一框架**：CG(§22.6) = 对称正定 Krylov 最优化 · GMRES(§22.7) = 一般情形残差最小化

</div>

### Krylov 子空间

!!! definition "定义 22.8 (Krylov 子空间)"
    给定 $n \times n$ 矩阵 $A$ 和向量 $\mathbf{b}$，第 $k$ 个 **Krylov 子空间**（Krylov subspace）定义为

    $$\mathcal{K}_k(A, \mathbf{b}) = \text{span}\{\mathbf{b}, A\mathbf{b}, A^2\mathbf{b}, \ldots, A^{k-1}\mathbf{b}\}.$$

    Krylov 子空间满足嵌套关系 $\mathcal{K}_1 \subseteq \mathcal{K}_2 \subseteq \cdots$，且存在最小的 $m \leq n$ 使得 $\mathcal{K}_m = \mathcal{K}_{m+1} = \cdots$。

!!! proposition "命题 22.1 (Krylov 子空间的性质)"
    1. $\dim \mathcal{K}_k(A, \mathbf{b}) \leq k$，等号不一定成立。
    2. 设 $A$ 的最小多项式为 $p(\lambda)$，$\deg p = d$，则 $\dim \mathcal{K}_k(A, \mathbf{b}) = \min(k, d_{\mathbf{b}})$，其中 $d_{\mathbf{b}}$ 是使 $\mathcal{K}_{d_{\mathbf{b}}}(A, \mathbf{b}) = \mathcal{K}_{d_{\mathbf{b}}+1}(A, \mathbf{b})$ 的最小整数。
    3. $\mathcal{K}_k(A, \mathbf{b}) = \{p(A)\mathbf{b} : p \text{ 是次数不超过 } k-1 \text{ 的多项式}\}$。

??? proof "证明"
    第 3 条：$\mathcal{K}_k(A, \mathbf{b})$ 由 $\mathbf{b}, A\mathbf{b}, \ldots, A^{k-1}\mathbf{b}$ 张成，而 $p(A)\mathbf{b} = c_0\mathbf{b} + c_1 A\mathbf{b} + \cdots + c_{k-1}A^{k-1}\mathbf{b}$ 恰好是它们的线性组合。反之，任何这样的线性组合都可以写成 $p(A)\mathbf{b}$ 的形式。

### Arnoldi 过程

!!! definition "定义 22.9 (Arnoldi 过程)"
    **Arnoldi 过程**（Arnoldi process）是构造 $\mathcal{K}_k(A, \mathbf{b})$ 的正交基的算法：

    **输入**：矩阵 $A$，向量 $\mathbf{b}$，迭代步数 $k$。

    1. $\mathbf{q}_1 = \mathbf{b} / \|\mathbf{b}\|$
    2. **For** $j = 1, 2, \ldots, k$：
        - $\mathbf{v} = A\mathbf{q}_j$
        - **For** $i = 1, \ldots, j$：$h_{ij} = \mathbf{q}_i^T \mathbf{v}$，$\mathbf{v} = \mathbf{v} - h_{ij}\mathbf{q}_i$
        - $h_{j+1,j} = \|\mathbf{v}\|$
        - 若 $h_{j+1,j} = 0$ 则停止（Krylov 子空间不变）
        - $\mathbf{q}_{j+1} = \mathbf{v} / h_{j+1,j}$

    **输出**：正交矩阵 $Q_k = [\mathbf{q}_1, \ldots, \mathbf{q}_k]$ 和上 Hessenberg 矩阵 $H_k = (h_{ij})_{k \times k}$。

!!! theorem "定理 22.7 (Arnoldi 关系)"
    Arnoldi 过程产生的矩阵满足

    $$AQ_k = Q_k H_k + h_{k+1,k} \, \mathbf{q}_{k+1} \mathbf{e}_k^T = Q_{k+1} \tilde{H}_k,$$

    其中 $\tilde{H}_k$ 是 $(k+1) \times k$ 矩阵。等价地，$Q_k^T A Q_k = H_k$。

??? proof "证明"
    由 Arnoldi 过程的第 $j$ 步：

    $$A\mathbf{q}_j = \sum_{i=1}^{j} h_{ij} \mathbf{q}_i + h_{j+1,j} \mathbf{q}_{j+1}.$$

    将 $j = 1, \ldots, k$ 的等式合并，得 $AQ_k = Q_k H_k + h_{k+1,k} \mathbf{q}_{k+1} \mathbf{e}_k^T$。

    由 $Q_k$ 的列正交性，左乘 $Q_k^T$ 得 $Q_k^T A Q_k = H_k$（利用 $Q_k^T \mathbf{q}_{k+1} = \mathbf{0}$）。

### Lanczos 过程

!!! theorem "定理 22.8 (Lanczos 过程)"
    当 $A$ 是**对称矩阵**时，Arnoldi 过程中的 $H_k$ 退化为**三对角矩阵** $T_k$，Arnoldi 过程简化为 **Lanczos 过程**（Lanczos process）：

    $$\beta_{j+1} \mathbf{q}_{j+1} = A\mathbf{q}_j - \alpha_j \mathbf{q}_j - \beta_j \mathbf{q}_{j-1},$$

    其中 $\alpha_j = \mathbf{q}_j^T A \mathbf{q}_j$，$\beta_{j+1} = \|A\mathbf{q}_j - \alpha_j \mathbf{q}_j - \beta_j \mathbf{q}_{j-1}\|$。

    每步只需存储当前和前一个向量，且只需一次矩阵-向量乘法。

??? proof "证明"
    当 $A = A^T$ 时，$H_k = Q_k^T A Q_k$ 也是对称的。由于 $H_k$ 是上 Hessenberg 矩阵且对称，它必须是三对角矩阵。设 $\alpha_j = h_{jj}$，$\beta_j = h_{j,j-1} = h_{j-1,j}$，则 Arnoldi 递推简化为三项递推。

!!! example "例 22.5"
    对对称矩阵 $A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 3 & 1 \\ 0 & 1 & 2 \end{pmatrix}$ 和 $\mathbf{b} = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}$，执行 Lanczos 过程：

    - $\mathbf{q}_1 = (1, 0, 0)^T$，$\alpha_1 = \mathbf{q}_1^T A \mathbf{q}_1 = 2$
    - $\mathbf{r}_1 = A\mathbf{q}_1 - 2\mathbf{q}_1 = (0, 1, 0)^T$，$\beta_2 = 1$，$\mathbf{q}_2 = (0, 1, 0)^T$
    - $\alpha_2 = \mathbf{q}_2^T A \mathbf{q}_2 = 3$
    - $\mathbf{r}_2 = A\mathbf{q}_2 - 3\mathbf{q}_2 - 1 \cdot \mathbf{q}_1 = (0, 0, 1)^T$，$\beta_3 = 1$，$\mathbf{q}_3 = (0, 0, 1)^T$

    得到三对角矩阵 $T_3 = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 3 & 1 \\ 0 & 1 & 2 \end{pmatrix} = A$（完全恢复）。

## 22.6 共轭梯度法

<div class="context-flow" markdown>

**CG 的本质**：在 Krylov 子空间中求 $A$-范数最佳近似 → 收敛速度 $O(\sqrt{\kappa})$（Chebyshev 多项式最优性）→ **预条件**降低有效 $\kappa$ 是实战关键

</div>

### 算法推导

!!! definition "定义 22.10 ($A$-共轭)"
    设 $A$ 是对称正定矩阵。向量 $\mathbf{p}$ 和 $\mathbf{q}$ 称为 **$A$-共轭的**（$A$-conjugate），如果 $\mathbf{p}^T A \mathbf{q} = 0$。

**共轭梯度法**（Conjugate Gradient method, CG）求解对称正定系统 $A\mathbf{x} = \mathbf{b}$，等价于最小化二次函数

$$\phi(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T A \mathbf{x} - \mathbf{b}^T \mathbf{x}.$$

!!! theorem "定理 22.9 (共轭梯度法的最优性)"
    CG 法第 $k$ 步的迭代 $\mathbf{x}_k$ 满足

    $$\mathbf{x}_k = \arg\min_{\mathbf{x} \in \mathbf{x}_0 + \mathcal{K}_k(A, \mathbf{r}_0)} \|\mathbf{x} - \mathbf{x}^*\|_A,$$

    其中 $\mathbf{r}_0 = \mathbf{b} - A\mathbf{x}_0$ 是初始残差，$\mathbf{x}^*$ 是精确解，$\|\mathbf{v}\|_A = \sqrt{\mathbf{v}^T A \mathbf{v}}$ 是 $A$-范数。

    等价地，CG 是在 Krylov 子空间中求 $A$-范数意义下最佳近似的方法。

??? proof "证明"
    CG 法构造 $A$-共轭搜索方向 $\{\mathbf{p}_0, \mathbf{p}_1, \ldots\}$。可以证明 $\text{span}\{\mathbf{p}_0, \ldots, \mathbf{p}_{k-1}\} = \mathcal{K}_k(A, \mathbf{r}_0)$。

    由于 $\mathbf{x}_k = \mathbf{x}_0 + \sum_{i=0}^{k-1} \alpha_i \mathbf{p}_i$，且 $\alpha_i$ 的选取使得误差 $\mathbf{e}_k = \mathbf{x}_k - \mathbf{x}^*$ 满足 $\mathbf{e}_k^T A \mathbf{p}_i = 0$（$i = 0, \ldots, k-1$），这正是在 $\mathbf{x}_0 + \mathcal{K}_k$ 中关于 $A$-内积的正交投影，等价于 $A$-范数最小化。

### 收敛性分析

<div class="context-flow" markdown>

**洞察**：CG 收敛速度 $\propto \sqrt{\kappa}$ 而非 $\kappa$——因为 Chebyshev 多项式在 $[\lambda_{\min}, \lambda_{\max}]$ 上的最优逼近性质，将迭代次数从 $O(\kappa)$ 降到 $O(\sqrt{\kappa})$

</div>

!!! theorem "定理 22.10 (CG 法的收敛速度)"
    设 $A$ 的条件数 $\kappa = \kappa_2(A) = \lambda_{\max}/\lambda_{\min}$。CG 法的误差满足

    $$\|\mathbf{x}_k - \mathbf{x}^*\|_A \leq 2 \left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k \|\mathbf{x}_0 - \mathbf{x}^*\|_A.$$

    因此，将误差减少到 $\epsilon$ 倍所需的迭代步数为 $O(\sqrt{\kappa} \log(1/\epsilon))$。

??? proof "证明"
    由 CG 的最优性，

    $$\|\mathbf{e}_k\|_A = \min_{p \in \mathcal{P}_k, \, p(0)=1} \max_{\lambda \in \sigma(A)} |p(\lambda)| \cdot \|\mathbf{e}_0\|_A,$$

    其中 $\mathcal{P}_k$ 是次数不超过 $k$ 的多项式集合。选取 Chebyshev 多项式

    $$p_k(\lambda) = \frac{T_k\left(\frac{\lambda_{\max}+\lambda_{\min}-2\lambda}{\lambda_{\max}-\lambda_{\min}}\right)}{T_k\left(\frac{\lambda_{\max}+\lambda_{\min}}{\lambda_{\max}-\lambda_{\min}}\right)},$$

    利用 $T_k$ 的性质可得上述估计。

!!! example "例 22.6"
    对 $\kappa = 100$ 的对称正定系统，CG 法将误差减少到 $10^{-6}$ 大约需要

    $$k \approx \frac{\sqrt{100}}{2} \ln(2/10^{-6}) \approx 5 \times 14.5 \approx 73$$

    步迭代。使用好的预条件子（preconditioner）$M$ 将有效条件数降至 $\kappa_{\text{eff}} = 10$，则只需约 $23$ 步。

### 预条件共轭梯度法

!!! definition "定义 22.11 (预条件共轭梯度法)"
    **预条件共轭梯度法**（Preconditioned CG, PCG）是对系统 $M^{-1}A\mathbf{x} = M^{-1}\mathbf{b}$ 应用 CG 法，其中 $M$ 是**预条件子**（preconditioner），满足 $M \approx A$ 且 $M^{-1}$ 易于计算。

    PCG 算法的核心步骤为：

    1. $\mathbf{r}_k = \mathbf{b} - A\mathbf{x}_k$
    2. 求解 $M\mathbf{z}_k = \mathbf{r}_k$
    3. $\beta_k = \mathbf{r}_k^T \mathbf{z}_k / \mathbf{r}_{k-1}^T \mathbf{z}_{k-1}$
    4. $\mathbf{p}_k = \mathbf{z}_k + \beta_k \mathbf{p}_{k-1}$
    5. $\alpha_k = \mathbf{r}_k^T \mathbf{z}_k / (\mathbf{p}_k^T A \mathbf{p}_k)$
    6. $\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k$

    常用的预条件子包括：不完全 Cholesky 分解、对称 Gauss-Seidel、多重网格方法等。

## 22.7 GMRES 方法

<div class="context-flow" markdown>

**推广**：CG 限于对称正定 → GMRES 处理一般方阵 · Arnoldi 关系将大问题转化为小 Hessenberg 最小二乘 → 收敛取决于特征值聚集程度

</div>

!!! definition "定义 22.12 (GMRES 方法)"
    **GMRES**（Generalized Minimal Residual）方法用于求解一般的（可能非对称的）线性方程组 $A\mathbf{x} = \mathbf{b}$。在第 $k$ 步，GMRES 在 Krylov 子空间中寻找残差 $2$-范数最小的近似：

    $$\mathbf{x}_k = \arg\min_{\mathbf{x} \in \mathbf{x}_0 + \mathcal{K}_k(A, \mathbf{r}_0)} \|\mathbf{b} - A\mathbf{x}\|_2.$$

!!! theorem "定理 22.11 (GMRES 与 Arnoldi 的关系)"
    GMRES 第 $k$ 步等价于求解最小二乘问题

    $$\min_{\mathbf{y} \in \mathbb{R}^k} \|\beta \mathbf{e}_1 - \tilde{H}_k \mathbf{y}\|_2,$$

    其中 $\beta = \|\mathbf{r}_0\|$，$\tilde{H}_k$ 是 Arnoldi 过程产生的 $(k+1) \times k$ 上 Hessenberg 矩阵。最终 $\mathbf{x}_k = \mathbf{x}_0 + Q_k \mathbf{y}_k$。

??? proof "证明"
    由 Arnoldi 关系 $AQ_k = Q_{k+1}\tilde{H}_k$。对 $\mathbf{x} = \mathbf{x}_0 + Q_k\mathbf{y}$，

    $$\mathbf{b} - A\mathbf{x} = \mathbf{r}_0 - AQ_k\mathbf{y} = \beta Q_{k+1}\mathbf{e}_1 - Q_{k+1}\tilde{H}_k\mathbf{y} = Q_{k+1}(\beta\mathbf{e}_1 - \tilde{H}_k\mathbf{y}).$$

    由 $Q_{k+1}$ 的列正交性，

    $$\|\mathbf{b} - A\mathbf{x}\|_2 = \|\beta\mathbf{e}_1 - \tilde{H}_k\mathbf{y}\|_2.$$

    因此 GMRES 的最优化问题转化为小规模（$(k+1) \times k$）最小二乘问题。

!!! proposition "命题 22.2 (GMRES 的收敛性)"
    1. GMRES 至多 $n$ 步收敛到精确解（在精确算术下）。
    2. 若 $A$ 的特征值远离零点且聚集，则 GMRES 收敛快。具体地，若 $A$ 可对角化，$A = X \Lambda X^{-1}$，则

    $$\frac{\|\mathbf{r}_k\|}{\|\mathbf{r}_0\|} \leq \kappa(X) \min_{p \in \mathcal{P}_k, \, p(0)=1} \max_{\lambda \in \sigma(A)} |p(\lambda)|.$$

??? proof "证明"
    第 1 条：$\mathcal{K}_n(A, \mathbf{r}_0) = \mathbb{R}^n$（当 $A$ 非奇异时），故 $\mathbf{x}^* \in \mathbf{x}_0 + \mathcal{K}_n$。

    第 2 条：由 GMRES 的最优性，$\|\mathbf{r}_k\| \leq \|p_k(A)\mathbf{r}_0\|$ 对任意满足 $p_k(0) = 1$ 的 $k$ 次多项式 $p_k$。若 $A = X\Lambda X^{-1}$，则 $\|p_k(A)\| \leq \|X\| \cdot \|p_k(\Lambda)\| \cdot \|X^{-1}\| = \kappa(X) \max_\lambda |p_k(\lambda)|$。

!!! example "例 22.7"
    对上三角矩阵 $A = I + N$（$N$ 严格上三角，$N^n = 0$），GMRES 在至多 $n$ 步内精确收敛（因为 $A$ 的最小多项式次数不超过 $n$）。实际上，若 $N^m = 0$（$m < n$），则 GMRES 在 $m$ 步内收敛。

## 22.8 特征值的数值计算

<div class="context-flow" markdown>

**层次**：幂迭代（最简单, $|\lambda_2/\lambda_1|^k$）→ 反迭代+Rayleigh 商（三次收敛）→ QR 迭代（全谱, 先 Hessenberg 化再迭代）

**链接**：QR 分解(Ch8) 在这里不是分解工具，而是迭代引擎 · Ch25 Rayleigh 商优化视角

</div>

### 幂迭代法

!!! definition "定义 22.13 (幂迭代法)"
    **幂迭代法**（power iteration）用于计算矩阵 $A$ 的模最大特征值及对应的特征向量：

    1. 选取初始向量 $\mathbf{q}_0$（$\|\mathbf{q}_0\| = 1$）
    2. **For** $k = 0, 1, 2, \ldots$：
        - $\mathbf{z}_{k+1} = A\mathbf{q}_k$
        - $\mathbf{q}_{k+1} = \mathbf{z}_{k+1} / \|\mathbf{z}_{k+1}\|$
        - $\lambda_{k+1} = \mathbf{q}_{k+1}^T A \mathbf{q}_{k+1}$（Rayleigh 商）

!!! theorem "定理 22.12 (幂迭代的收敛性)"
    设 $A$ 有特征值 $|\lambda_1| > |\lambda_2| \geq \cdots \geq |\lambda_n|$，初始向量 $\mathbf{q}_0$ 在 $\lambda_1$ 对应的特征向量方向上有非零分量，则

    $$|\lambda_{k} - \lambda_1| = O\left(\left|\frac{\lambda_2}{\lambda_1}\right|^{2k}\right), \quad \text{dist}(\mathbf{q}_k, \text{span}\{\mathbf{v}_1\}) = O\left(\left|\frac{\lambda_2}{\lambda_1}\right|^k\right).$$

    使用 Rayleigh 商时，特征值的收敛速度是特征向量的两倍。

??? proof "证明"
    设 $\mathbf{q}_0 = \sum_{i=1}^n c_i \mathbf{v}_i$，$c_1 \neq 0$。则

    $$A^k \mathbf{q}_0 = \sum_{i=1}^n c_i \lambda_i^k \mathbf{v}_i = c_1 \lambda_1^k \left[\mathbf{v}_1 + \sum_{i=2}^n \frac{c_i}{c_1}\left(\frac{\lambda_i}{\lambda_1}\right)^k \mathbf{v}_i\right].$$

    由 $|\lambda_i/\lambda_1| < 1$（$i \geq 2$），当 $k \to \infty$ 时 $A^k \mathbf{q}_0 / \|A^k \mathbf{q}_0\| \to \mathbf{v}_1/\|\mathbf{v}_1\|$，收敛速度为 $|\lambda_2/\lambda_1|^k$。

    Rayleigh 商 $\lambda_k = \mathbf{q}_k^T A \mathbf{q}_k$ 的误差为 $O(\|\mathbf{q}_k - \mathbf{v}_1\|^2) = O(|\lambda_2/\lambda_1|^{2k})$。

### QR 算法

!!! definition "定义 22.14 (QR 迭代)"
    **基本 QR 迭代**（QR iteration）：

    1. 令 $A_0 = A$
    2. **For** $k = 0, 1, 2, \ldots$：
        - 计算 QR 分解 $A_k = Q_k R_k$
        - 令 $A_{k+1} = R_k Q_k$

    **带位移的 QR 迭代**：

    1. 令 $A_0 = A$
    2. **For** $k = 0, 1, 2, \ldots$：
        - 选取位移 $\mu_k$（通常为 $A_k$ 的 Wilkinson 位移）
        - 计算 QR 分解 $A_k - \mu_k I = Q_k R_k$
        - 令 $A_{k+1} = R_k Q_k + \mu_k I$

!!! theorem "定理 22.13 (QR 迭代的收敛性)"
    设 $A$ 的特征值满足 $|\lambda_1| > |\lambda_2| > \cdots > |\lambda_n|$（模互不相等）。则基本 QR 迭代中，$A_k$ 收敛到上三角矩阵，对角元素收敛到 $\lambda_1, \ldots, \lambda_n$。

    收敛速度：$A_k$ 的 $(i,j)$ 元素（$i > j$）以 $|\lambda_i/\lambda_j|^k$ 的速度趋于零。

    带 Wilkinson 位移的 QR 迭代具有**三次收敛**速度（对对称矩阵）。

??? proof "证明"
    注意 $A_{k+1} = R_k Q_k = Q_k^T (Q_k R_k) Q_k = Q_k^T A_k Q_k$，即 $A_{k+1}$ 与 $A_k$ 正交相似。设 $\hat{Q}_k = Q_0 Q_1 \cdots Q_{k-1}$，则 $A_k = \hat{Q}_k^T A \hat{Q}_k$。

    另一方面，$A^k = (Q_0 R_0)(Q_1 R_1) \cdots = \hat{Q}_k \hat{R}_k$（其中 $\hat{R}_k = R_0 R_1 \cdots R_{k-1}$）。这与 $A^k$ 的 QR 分解密切相关。利用 $A$ 的 Schur 分解 $A = USU^*$ 和 $A^k = US^kU^*$ 的 QR 分解的渐近分析，可证明 $A_k$ 的下三角部分趋于零。

### Hessenberg 化简

!!! theorem "定理 22.14 (Hessenberg 预处理)"
    任何矩阵 $A$ 可以通过 $O(\frac{2}{3}n^3)$ 次运算用正交相似变换化为上 Hessenberg 形式（upper Hessenberg form）$H = Q^T A Q$，其中 $h_{ij} = 0$（$i > j + 1$）。

    对 Hessenberg 矩阵执行一步 QR 迭代只需 $O(n^2)$ 运算（而非 $O(n^3)$），这使得 QR 算法在实际中高效可行。

??? proof "证明"
    使用 Householder 变换逐列消去。在第 $k$ 步，选取 Householder 矩阵 $P_k$ 将 $A$ 第 $k$ 列第 $k+2$ 行以下的元素消为零，执行相似变换 $A \leftarrow P_k A P_k$。保持相似性意味着从右边也要乘 $P_k$，这会在第 $k$ 列右边引入非零元素但不破坏之前各列的 Hessenberg 结构。

    对 Hessenberg 矩阵做 QR 分解只需 $n-1$ 个 Givens 旋转（每个 $O(n)$），总共 $O(n^2)$。

!!! example "例 22.8"
    将矩阵 $A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 0 \end{pmatrix}$ 先化为 Hessenberg 形式，再用 QR 迭代求特征值。

    Hessenberg 化简：使用 Householder 变换消去 $a_{31}$。取 $P_1$ 使 $(a_{21}, a_{31})^T = (4, 7)^T$ 变为 $(\pm\sqrt{65}, 0)^T$，则

    $$H = P^T A P = \begin{pmatrix} 1 & * & * \\ -\sqrt{65} & * & * \\ 0 & * & * \end{pmatrix}$$

    已是 Hessenberg 形式（$3 \times 3$ 矩阵只需一步）。之后的 QR 迭代在此 Hessenberg 矩阵上进行，每步只需 $O(n^2)$ 运算。

## 22.9 稀疏矩阵

<div class="context-flow" markdown>

**实践**：科学计算中 $n \sim 10^6$ 但 $\text{nnz} \sim O(n)$（如有限元/差分）→ COO/CSR/CSC 格式 → 矩阵-向量乘 $O(\text{nnz})$ 使 Krylov 方法可行

</div>

### 存储格式

!!! definition "定义 22.15 (稀疏矩阵存储格式)"
    **稀疏矩阵**（sparse matrix）是大部分元素为零的矩阵。设 $n \times n$ 矩阵有 $\text{nnz}$ 个非零元素，$\text{nnz} \ll n^2$。常用的存储格式：

    - **COO（Coordinate format）**：存储三个数组：行索引 `row[]`、列索引 `col[]`、值 `val[]`，每个长度为 $\text{nnz}$。
    - **CSR（Compressed Sparse Row）**：存储数组 `val[]`（长 $\text{nnz}$）、`col_idx[]`（长 $\text{nnz}$）、`row_ptr[]`（长 $n+1$）。其中 `row_ptr[i]` 到 `row_ptr[i+1]-1` 指示第 $i$ 行非零元素在 `val[]` 和 `col_idx[]` 中的位置。
    - **CSC（Compressed Sparse Column）**：类似 CSR 但按列压缩。

    | 格式 | 存储量 | 矩阵-向量乘法 | 访问 $a_{ij}$ |
    |------|-------|-------------|-------------|
    | COO  | $3 \cdot \text{nnz}$ | $O(\text{nnz})$ | $O(\text{nnz})$ |
    | CSR  | $2 \cdot \text{nnz} + n + 1$ | $O(\text{nnz})$（行优先） | $O(\log d_i)$ |
    | CSC  | $2 \cdot \text{nnz} + n + 1$ | $O(\text{nnz})$（列优先） | $O(\log d_j)$ |

    其中 $d_i$ 是第 $i$ 行的非零元素个数。

!!! example "例 22.9"
    考虑矩阵

    $$A = \begin{pmatrix} 5 & 0 & 0 & 3 \\ 0 & 8 & 0 & 0 \\ 0 & 0 & 3 & 0 \\ 0 & 6 & 0 & 1 \end{pmatrix}.$$

    **COO 格式**：

    - `val = [5, 3, 8, 3, 6, 1]`
    - `row = [0, 0, 1, 2, 3, 3]`
    - `col = [0, 3, 1, 2, 1, 3]`

    **CSR 格式**：

    - `val = [5, 3, 8, 3, 6, 1]`
    - `col_idx = [0, 3, 1, 2, 1, 3]`
    - `row_ptr = [0, 2, 3, 4, 6]`

    存储量：$\text{nnz} = 6$，COO 需要 $18$ 个数，CSR 需要 $17$ 个数，而稠密存储需要 $16$ 个数。当 $n$ 很大、$\text{nnz} \ll n^2$ 时，稀疏存储的优势显著。

!!! proposition "命题 22.3 (稀疏矩阵-向量乘法的效率)"
    CSR 格式下矩阵-向量乘法 $\mathbf{y} = A\mathbf{x}$ 的算法为：

    **For** $i = 0, \ldots, n-1$：

    $$y_i = \sum_{k=\text{row\_ptr}[i]}^{\text{row\_ptr}[i+1]-1} \text{val}[k] \cdot x_{\text{col\_idx}[k]}.$$

    计算量为 $2 \cdot \text{nnz}$ 次浮点运算，而稠密矩阵-向量乘法需要 $2n^2$ 次。对于五点差分矩阵等结构化稀疏矩阵，$\text{nnz} = O(n)$，因此稀疏方法的效率远优于稠密方法。

??? proof "证明"
    由 CSR 格式的定义，第 $i$ 行的非零元素存储在 `val[row_ptr[i]:row_ptr[i+1]]` 中，对应列索引为 `col_idx[row_ptr[i]:row_ptr[i+1]]`。矩阵-向量乘法 $y_i = \sum_j a_{ij} x_j$ 中只需遍历非零的 $a_{ij}$，因此总运算量为 $\sum_i (\text{第 } i \text{ 行非零元素数}) \times 2 = 2 \cdot \text{nnz}$。

!!! example "例 22.10"
    **五点差分格式**。在 $N \times N$ 的网格上离散 Laplace 方程 $-\Delta u = f$，得到 $n \times n$（$n = N^2$）稀疏矩阵，每行至多 $5$ 个非零元素：

    $$A = \begin{pmatrix} T & -I & & \\ -I & T & -I & \\ & \ddots & \ddots & \ddots \\ & & -I & T \end{pmatrix}, \quad T = \begin{pmatrix} 4 & -1 & & \\ -1 & 4 & -1 & \\ & \ddots & \ddots & \ddots \\ & & -1 & 4 \end{pmatrix}.$$

    此时 $\text{nnz} \leq 5n$，使用 CSR 存储和稀疏矩阵-向量乘法可以高效实现 CG 等迭代法。对 $N = 1000$（$n = 10^6$），稠密存储需要 $10^{12}$ 个浮点数（约 $8$ TB），而稀疏存储仅需约 $5 \times 10^6$ 个（约 $40$ MB）。

## 练习题

1. **[误差分析] 在浮点数运算中，为什么两个非常相近的数相减会导致严重的“灾难性相消”（Catastrophic Cancellation）？**
   ??? success "参考答案"
       因为相近的两个数在计算机中前几位有效数字相同，相减后这些有效数字变为 0，暴露出低位的舍入误差（甚至被当做有效数字放大），导致相对误差急剧增大。

2. **[条件数] 求解 $A\mathbf{x} = \mathbf{b}$ 时，矩阵的条件数 $\kappa(A)$ 很大，这会对数值求解产生什么实际影响？**
   ??? success "参考答案"
       条件数极大意味着矩阵接近奇异（病态）。此时右端向量 $\mathbf{b}$ 或系数矩阵 $A$ 上的极小舍入误差，会被条件数放大，导致计算出的解 $\hat{\mathbf{x}}$ 产生巨大的相对误差，变得毫无物理意义。

3. **[直接法] 相比于直接计算 $A^{-1}$ 然后乘 $\mathbf{b}$，为什么数值软件更倾向于先对 $A$ 做 LU 分解再回代？**
   ??? success "参考答案"
       首先计算 $A^{-1}$ 需要大约 $\frac{4}{3}n^3$ 的浮点运算，而 LU 分解加回代只需要大约 $\frac{2}{3}n^3$，效率高一倍；其次，求逆容易在病态矩阵中累积更多舍入误差，数值稳定性较差。

4. **[迭代法] 迭代法（如 Jacobi、Gauss-Seidel）能保证对任何非奇异矩阵都收敛吗？**
   ??? success "参考答案"
       不能。它们只在特定条件下保证收敛（如严格对角占优矩阵，或正定矩阵）。一般矩阵需要判断迭代矩阵的谱半径是否严格小于 1（$\rho(M) < 1$）。

5. **[Krylov子空间] 共轭梯度法（CG）为什么只适用于对称正定矩阵？**
   ??? success "参考答案"
       因为 CG 算法的推导基于寻找二次型 $f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - \mathbf{b}^T\mathbf{x}$ 的极小值。只有当 $A$ 为对称正定矩阵时，该函数的平稳点才真正对应于唯一的全局极小值，从而保证基于梯度的搜索是收敛的。

6. **[特征值] 经典的幂法（Power Method）能够求出矩阵的哪个特征值？**
   ??? success "参考答案"
       能够求出**模长最大**（主导）的特征值及其对应的特征向量。

7. **[特征值] 反幂法（Inverse Power Method）结合位移（Shift）的优势是什么？**
   ??? success "参考答案"
       反幂法作用于 $(A - \mu I)^{-1}$，其最大特征值对应于 $A$ 中最接近 $\mu$ 的特征值。这使我们能够精确地“瞄准”并快速求出矩阵中任何特定位置的特征值。

8. **[QR算法] 在实际使用 QR 迭代求解特征值之前，为什么要先将矩阵化为 Hessenberg 形式？**
   ??? success "参考答案"
       如果直接对稠密矩阵做 QR 迭代，每一步都需要 $O(n^3)$ 的计算量；而 Hessenberg 矩阵由于次下对角线以下全为 0，每步 QR 迭代仅需 $O(n^2)$ 计算量，且 Hessenberg 结构在迭代中保持不变，极大地加速了算法。

9. **[稀疏矩阵] 在机器学习和科学计算中，为什么需要 CSR 或 CSC 这样的稀疏矩阵存储格式？**
   ??? success "参考答案"
       现实中的网络、图、偏微分方程离散矩阵等，绝大部分元素都是 0。用稠密数组存储会导致内存爆炸（$O(n^2)$），而 CSR/CSC 只存储非零元素及其索引，将内存和矩阵乘法时间都降到了 $O(\text{nnz})$。

10. **[爱因斯坦思考题] 纯数学中，矩阵要么可逆，要么奇异。而在数值计算的现实世界中，这种绝对的二分法被“条件数”的连续统取代。这与物理学中从“理想模型”走向“实验测量”有何相似之处？**
    ??? success "参考答案"
        物理学中不存在绝对精确的测量（就像浮点数无法表示所有实数），任何观测都带有误差（不确定性）。纯数学的“奇异（$\det = 0$）”是一种理想的绝对状态，而在存在测量误差的真实物理世界里，我们关心的是系统面对扰动时放大的倍率（条件数）。一个高度病态的“可逆”系统，在物理观测上与“不可逆”的奇异系统其实是没有区别的。

## 本章小结

本章将线性代数从纯理论的理想殿堂引入到了充满误差与计算资源限制的现实世界，主要内容包括：

1. **数值稳定性与条件数**：确立了浮点数运算的误差模型，并利用条件数 $\kappa(A)$ 严格量化了方程求解时由于舍入误差引起的解的漂移。
2. **直接法与预处理**：回顾了 LU、Cholesky 分解在工业软件中的主导地位，并探讨了如何通过主元选取（Pivoting）维持数值稳定。
3. **迭代法体系**：介绍了从经典的 Jacobi、Gauss-Seidel 迭代到现代大规模问题首选的 Krylov 子空间方法（如针对正定阵的共轭梯度法 CG 和一般阵的 GMRES）。
4. **特征值的数值算法**：从求解单一特征值的幂法及其变体出发，深入到求解全谱特征值的工业标准——带位移的 QR 算法（结合 Hessenberg 预处理）。
5. **稀疏矩阵计算**：介绍了以 CSR/CSC 为代表的稀疏数据结构设计，说明了在大规模稀疏图或偏微分网格中，如何打破 $O(n^3)$ 的算力瓶颈。
