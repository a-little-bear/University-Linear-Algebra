# 第 25 章 线性代数在优化中的应用

<div class="context-flow" markdown>

**前置**：SVD/特征值(Ch6-8) · 正定性(Ch7) · 流形优化(Ch24)

**脉络**：LP(基=列选取) → 最小二乘(QR/SVD) → SDP(半正定锥) → 矩阵补全/压缩感知(核范数/$\ell_1$) → PCA/Rayleigh 商

**延伸**：半定规划在组合优化（MAX-CUT 的 Goemans-Williamson 近似）、控制理论（鲁棒控制）、量子信息（纠缠判定）中至关重要；压缩感知革新了 MRI 成像、雷达信号处理和天文观测

</div>

优化理论与线性代数有着深刻而广泛的联系。线性代数不仅为优化问题提供了建模语言和分析工具，其核心分解方法（如 QR 分解、SVD、特征值分解）更直接构成了许多优化算法的计算骨架。本章将从线性规划、最小二乘、半定规划、矩阵补全、压缩感知、主成分分析、低秩近似到特征值优化，系统展示线性代数在优化中的核心作用。

---

## 25.1 线性规划基础

<div class="context-flow" markdown>

**线性代数视角**：LP 最优解 = **基本可行解** = 选 $m$ 列构成可逆 $A_B$ → 单纯形法的每步 = 解线性方程组 $A_B^{-1}\mathbf{b}$ + 列交换(LU 更新)

**链接**：Ch22 LU 分解在此处直接应用

</div>

线性规划（Linear Programming, LP）是最基本的优化问题类型，其理论和算法根植于线性代数。

!!! definition "定义 25.1 (线性规划标准形)"
    **线性规划**（linear program）的标准形为

    $$
    \min_{\mathbf{x} \in \mathbb{R}^n} \mathbf{c}^T \mathbf{x}, \quad \text{s.t.} \quad A\mathbf{x} = \mathbf{b}, \quad \mathbf{x} \ge \mathbf{0},
    $$

    其中 $A \in \mathbb{R}^{m \times n}$（$m < n$）为约束矩阵，$\mathbf{b} \in \mathbb{R}^m$ 为右端向量，$\mathbf{c} \in \mathbb{R}^n$ 为目标向量。不等式 $\mathbf{x} \ge \mathbf{0}$ 表示分量非负。

!!! definition "定义 25.2 (基本可行解)"
    设 $A \in \mathbb{R}^{m \times n}$，$\operatorname{rank}(A) = m$。$A$ 的一个 **基**（basis）$B$ 是 $m$ 个列的索引集，使得 $A_B$（由这 $m$ 列组成的子矩阵）可逆。对应的 **基本可行解**（basic feasible solution, BFS）为

    $$
    \mathbf{x}_B = A_B^{-1} \mathbf{b}, \quad \mathbf{x}_N = \mathbf{0},
    $$

    其中 $N$ 为非基变量的索引集。若 $\mathbf{x}_B \ge \mathbf{0}$，则称该基为可行基。

!!! theorem "定理 25.1 (线性规划基本定理)"
    若线性规划有最优解，则存在一个基本可行解是最优解。

??? proof "证明"
    设 $\mathbf{x}^*$ 为最优解。若 $\mathbf{x}^*$ 不是基本可行解，则其正分量的个数 $k$ 小于 $m$，或者正分量对应的 $A$ 的列线性相关。

    **情形 1**：若正分量对应列线性无关且 $k < m$，可以用其他列补充为基，从而 $\mathbf{x}^*$ 是退化的基本可行解。

    **情形 2**：若正分量对应列线性相关，存在非零 $\mathbf{d}$ 使得 $A\mathbf{d} = \mathbf{0}$，$d_j = 0$（$j$ 不在正分量集中）。考虑 $\mathbf{x}^* + t\mathbf{d}$ 和 $\mathbf{x}^* - t\mathbf{d}$，两者都满足 $A\mathbf{x} = \mathbf{b}$。调整 $t$ 使得某个正分量变为零，同时保持非负性，且目标值不增。重复此过程直到正分量对应列线性无关，即得到基本可行解。$\blacksquare$

!!! theorem "定理 25.2 (单纯形法的线性代数基础)"
    单纯形法的每一步执行如下线性代数运算：

    1. **求解基本可行解**：$\mathbf{x}_B = A_B^{-1}\mathbf{b}$。
    2. **计算简约成本**：$\bar{c}_j = c_j - \mathbf{c}_B^T A_B^{-1} \mathbf{a}_j$，$j \in N$。
    3. **确定进基变量**：选择 $\bar{c}_j < 0$ 的非基变量 $j$。
    4. **计算方向**：$\mathbf{d}_B = A_B^{-1} \mathbf{a}_j$。
    5. **确定出基变量**：$\theta^* = \min_{i: d_{B_i} > 0} \frac{(x_B)_i}{d_{B_i}}$（最小比值规则）。
    6. **基交换**（pivot）：更新 $B$ 和 $A_B^{-1}$。

    核心运算是线性方程组求解 $A_B^{-1}\mathbf{b}$ 和 $A_B^{-1}\mathbf{a}_j$。实际实现中使用 LU 分解并进行列交换更新。

??? proof "证明"
    简约成本 $\bar{c}_j$ 表示非基变量 $x_j$ 增加一个单位时目标值的变化率。设当前基为 $B$，基本可行解为 $\mathbf{x}_B = A_B^{-1}\mathbf{b}$。若 $x_j$ 从 $0$ 增加到 $\theta$，为保持 $A\mathbf{x} = \mathbf{b}$，需要 $\mathbf{x}_B \to \mathbf{x}_B - \theta A_B^{-1}\mathbf{a}_j$。目标值变化为

    $$
    \Delta z = c_j \theta - \mathbf{c}_B^T (A_B^{-1}\mathbf{a}_j) \theta = \bar{c}_j \theta.
    $$

    当 $\bar{c}_j < 0$ 时，增加 $x_j$ 能降低目标值。$\theta^*$ 由非负性约束 $\mathbf{x}_B - \theta A_B^{-1}\mathbf{a}_j \ge \mathbf{0}$ 的最紧约束确定。$\blacksquare$

!!! example "例 25.1"
    **简单线性规划求解。**

    $$
    \min -x_1 - 2x_2, \quad \text{s.t.} \quad x_1 + x_2 + s_1 = 4, \quad x_1 + 3x_2 + s_2 = 6, \quad \mathbf{x}, \mathbf{s} \ge \mathbf{0}.
    $$

    初始基 $B = \{s_1, s_2\}$，$A_B = I_2$。基本可行解 $(s_1, s_2) = (4, 6)$，$z = 0$。

    简约成本 $\bar{c}_{x_1} = -1$，$\bar{c}_{x_2} = -2$。选 $x_2$ 进基，$\mathbf{d}_B = A_B^{-1}\mathbf{a}_{x_2} = (1, 3)^T$。比值：$4/1 = 4$，$6/3 = 2$，故 $\theta^* = 2$，$s_2$ 出基。新基本可行解 $(x_2, s_1) = (2, 2)$，$z = -4$。继续迭代直至所有简约成本非负。

---

## 25.2 最小二乘问题

<div class="context-flow" markdown>

**三种求解**：正规方程($A^TA$, 条件数平方) → QR 分解(稳定, Ch8) → SVD(最通用，伪逆+正则化)

**Tikhonov 正则化** = 谱过滤：大 $\sigma_i$ 保留，小 $\sigma_i$ 压制

</div>

最小二乘问题是线性代数与优化交汇的经典领域。

!!! definition "定义 25.3 (最小二乘问题)"
    给定 $A \in \mathbb{R}^{m \times n}$（$m \ge n$）和 $\mathbf{b} \in \mathbb{R}^m$，**最小二乘问题**（least squares problem）为

    $$
    \min_{\mathbf{x} \in \mathbb{R}^n} \|A\mathbf{x} - \mathbf{b}\|_2^2.
    $$

!!! definition "定义 25.4 (Tikhonov 正则化)"
    **Tikhonov 正则化**（Tikhonov regularization）或岭回归（ridge regression）问题为

    $$
    \min_{\mathbf{x} \in \mathbb{R}^n} \|A\mathbf{x} - \mathbf{b}\|_2^2 + \lambda \|\mathbf{x}\|_2^2, \quad \lambda > 0.
    $$

    正则化参数 $\lambda$ 控制解的光滑性与数据拟合之间的权衡。

!!! theorem "定理 25.3 (正规方程)"
    最小二乘问题 $\min \|A\mathbf{x} - \mathbf{b}\|_2^2$ 的解满足 **正规方程**（normal equations）：

    $$
    A^T A \mathbf{x} = A^T \mathbf{b}.
    $$

    当 $A$ 列满秩时，解唯一：$\mathbf{x}^* = (A^T A)^{-1} A^T \mathbf{b}$。Tikhonov 正则化问题的解为

    $$
    \mathbf{x}_\lambda = (A^T A + \lambda I)^{-1} A^T \mathbf{b}.
    $$

??? proof "证明"
    令 $f(\mathbf{x}) = \|A\mathbf{x} - \mathbf{b}\|_2^2 = \mathbf{x}^T A^T A \mathbf{x} - 2\mathbf{b}^T A\mathbf{x} + \|\mathbf{b}\|^2$。必要条件 $\nabla f(\mathbf{x}) = 2A^TA\mathbf{x} - 2A^T\mathbf{b} = \mathbf{0}$ 即正规方程。$A^TA$ 半正定，若 $A$ 列满秩则 $A^TA$ 正定，解唯一。

    对 Tikhonov 正则化，$g(\mathbf{x}) = f(\mathbf{x}) + \lambda\|\mathbf{x}\|^2$，$\nabla g = 2(A^TA + \lambda I)\mathbf{x} - 2A^T\mathbf{b} = \mathbf{0}$。由于 $\lambda > 0$，$A^TA + \lambda I$ 正定，解唯一。$\blacksquare$

!!! theorem "定理 25.4 (QR 分解求解最小二乘)"
    设 $A = QR$，其中 $Q \in \mathbb{R}^{m \times n}$，$Q^TQ = I_n$，$R \in \mathbb{R}^{n \times n}$ 为上三角矩阵。则最小二乘解为

    $$
    \mathbf{x}^* = R^{-1} Q^T \mathbf{b}.
    $$

    数值上，QR 分解比正规方程更稳定，因为 $\kappa(A^TA) = \kappa(A)^2$。

??? proof "证明"
    $\|A\mathbf{x} - \mathbf{b}\|^2 = \|QR\mathbf{x} - \mathbf{b}\|^2$。设 $Q_\perp$ 为 $Q$ 的正交补，使得 $[Q \ Q_\perp]$ 为正交矩阵。则

    $$
    \|QR\mathbf{x} - \mathbf{b}\|^2 = \|R\mathbf{x} - Q^T\mathbf{b}\|^2 + \|Q_\perp^T \mathbf{b}\|^2.
    $$

    第二项与 $\mathbf{x}$ 无关。最小化第一项：令 $R\mathbf{x} = Q^T\mathbf{b}$，即 $\mathbf{x} = R^{-1}Q^T\mathbf{b}$。$\blacksquare$

!!! theorem "定理 25.5 (SVD 求解最小二乘与伪逆)"
    设 $A = U\Sigma V^T$ 为 SVD，$\Sigma = \operatorname{diag}(\sigma_1, \ldots, \sigma_r, 0, \ldots, 0)$。最小二乘问题的最小范数解为

    $$
    \mathbf{x}^+ = A^+ \mathbf{b} = V \Sigma^+ U^T \mathbf{b} = \sum_{i=1}^{r} \frac{\mathbf{u}_i^T \mathbf{b}}{\sigma_i} \mathbf{v}_i,
    $$

    其中 $A^+ = V\Sigma^+ U^T$ 为 Moore-Penrose 伪逆，$\Sigma^+ = \operatorname{diag}(\sigma_1^{-1}, \ldots, \sigma_r^{-1}, 0, \ldots, 0)$。

    Tikhonov 正则化解的 SVD 表达为

    $$
    \mathbf{x}_\lambda = \sum_{i=1}^{r} \frac{\sigma_i}{\sigma_i^2 + \lambda} (\mathbf{u}_i^T \mathbf{b}) \mathbf{v}_i.
    $$

??? proof "证明"
    代入 $A = U\Sigma V^T$：

    $$
    \|A\mathbf{x} - \mathbf{b}\|^2 = \|\Sigma V^T\mathbf{x} - U^T\mathbf{b}\|^2.
    $$

    令 $\mathbf{y} = V^T\mathbf{x}$，$\mathbf{c} = U^T\mathbf{b}$，则 $\|\Sigma\mathbf{y} - \mathbf{c}\|^2 = \sum_{i=1}^r (\sigma_i y_i - c_i)^2 + \sum_{i=r+1}^m c_i^2$。最小化得 $y_i = c_i/\sigma_i$（$i \le r$），$y_i$ 任意（$i > r$）。最小范数解取 $y_i = 0$（$i > r$），即 $\mathbf{x}^+ = V\Sigma^+ U^T\mathbf{b}$。

    对 Tikhonov 正则化，$(A^TA + \lambda I)\mathbf{x} = A^T\mathbf{b}$ 变为 $(\Sigma^2 + \lambda I)\mathbf{y} = \Sigma \mathbf{c}$，故 $y_i = \sigma_i c_i / (\sigma_i^2 + \lambda)$。$\blacksquare$

!!! example "例 25.2"
    **条件数对最小二乘的影响。**

    设 $A = \begin{pmatrix} 1 & 1 \\ \epsilon & 0 \\ 0 & \epsilon \end{pmatrix}$，$\mathbf{b} = \begin{pmatrix} 2 \\ 0 \\ 0 \end{pmatrix}$。当 $\epsilon$ 很小时，$\kappa(A) \approx \sqrt{2}/\epsilon$。正规方程 $A^TA = \begin{pmatrix} 1+\epsilon^2 & 1 \\ 1 & 1+\epsilon^2 \end{pmatrix}$，$\kappa(A^TA) \approx 2/\epsilon^2$。使用 QR 分解或 SVD 可以避免条件数平方的精度损失。

!!! example "例 25.3"
    **Tikhonov 正则化的谱过滤效应。**

    SVD 表达 $\mathbf{x}_\lambda = \sum_i \frac{\sigma_i}{\sigma_i^2 + \lambda} (\mathbf{u}_i^T\mathbf{b})\mathbf{v}_i$ 表明：对于大奇异值 $\sigma_i \gg \sqrt{\lambda}$，滤波因子 $\sigma_i/(\sigma_i^2 + \lambda) \approx 1/\sigma_i$（与伪逆相同）；对于小奇异值 $\sigma_i \ll \sqrt{\lambda}$，滤波因子 $\approx \sigma_i/\lambda \approx 0$（被抑制）。这有效地截断了噪声在小奇异值方向上的放大。

---

## 25.3 半定规划（SDP）

<div class="context-flow" markdown>

**推广链**：LP($\mathbf{x} \ge 0$) → SDP($X \succeq 0$)——将非负约束推广为**半正定锥**约束 · 强对偶性 + 互补松弛 $X^*S^* = 0$

**威力**：MAX-CUT 的 Goemans-Williamson 松弛(0.878 近似比) · 核范数最小化 = SDP(§25.4)

</div>

半定规划是线性规划在矩阵空间上的自然推广。

!!! definition "定义 25.5 (半定规划)"
    **半定规划**（Semidefinite Programming, SDP）的标准形为

    $$
    \min_{X \in \operatorname{Sym}(n)} \langle C, X \rangle, \quad \text{s.t.} \quad \langle A_i, X \rangle = b_i \; (i = 1, \ldots, m), \quad X \succeq 0,
    $$

    其中 $\langle A, B \rangle = \operatorname{tr}(A^T B)$ 为矩阵内积，$C, A_i \in \operatorname{Sym}(n)$，$X \succeq 0$ 表示 $X$ 半正定。

!!! definition "定义 25.6 (SDP 对偶问题)"
    标准形 SDP 的 **对偶问题**（dual problem）为

    $$
    \max_{\mathbf{y} \in \mathbb{R}^m} \mathbf{b}^T \mathbf{y}, \quad \text{s.t.} \quad \sum_{i=1}^{m} y_i A_i + S = C, \quad S \succeq 0.
    $$

    原问题的目标值（最小化）始终不小于对偶问题的目标值（最大化）：$\langle C, X \rangle \ge \mathbf{b}^T\mathbf{y}$（弱对偶性）。

!!! theorem "定理 25.6 (SDP 强对偶性)"
    若原问题（或对偶问题）满足 **Slater 条件**（存在严格可行解 $X \succ 0$，使得 $\langle A_i, X \rangle = b_i$），则强对偶性成立：

    $$
    \min_X \langle C, X \rangle = \max_{\mathbf{y}} \mathbf{b}^T \mathbf{y},
    $$

    且对偶最优值可以达到。在最优解 $(X^*, \mathbf{y}^*, S^*)$ 处，互补松弛条件 $X^* S^* = 0$ 成立。

??? proof "证明"
    **证明思路（锥规划对偶理论）。** SDP 是锥规划的特例，约束锥为半正定锥 $\mathcal{S}_+^n$。Slater 条件保证了约束的正则性。

    弱对偶性的证明：设 $X$ 原可行，$(\mathbf{y}, S)$ 对偶可行，则

    $$
    \langle C, X \rangle - \mathbf{b}^T\mathbf{y} = \langle C, X \rangle - \sum_i y_i \langle A_i, X \rangle = \langle C - \sum_i y_i A_i, X \rangle = \langle S, X \rangle \ge 0,
    $$

    最后不等式由 $S \succeq 0$、$X \succeq 0$ 保证（$\langle S, X \rangle = \operatorname{tr}(SX) \ge 0$，因为两个半正定矩阵之积的迹非负）。

    强对偶性需要 Slater 条件和凸分析中的分离超平面定理，证明对偶间隙为零。互补松弛 $\langle S^*, X^* \rangle = 0$ 意味着 $\operatorname{tr}(S^*X^*) = 0$，由于 $S^*, X^* \succeq 0$，这等价于 $S^*X^* = 0$。$\blacksquare$

!!! theorem "定理 25.7 (LP 是 SDP 的特例)"
    线性规划 $\min \mathbf{c}^T\mathbf{x}$，s.t. $A\mathbf{x} = \mathbf{b}$，$\mathbf{x} \ge \mathbf{0}$ 可以表述为 SDP：取 $X = \operatorname{diag}(\mathbf{x})$，$C = \operatorname{diag}(\mathbf{c})$，约束 $X \succeq 0$（即 $X$ 是非负对角矩阵）。

??? proof "证明"
    令 $X = \operatorname{diag}(x_1, \ldots, x_n)$。则 $\langle C, X \rangle = \operatorname{tr}(\operatorname{diag}(\mathbf{c}) \operatorname{diag}(\mathbf{x})) = \mathbf{c}^T\mathbf{x}$。约束 $\langle A_i, X \rangle = b_i$ 可以编码 $A\mathbf{x} = \mathbf{b}$ 的每一行。$X \succeq 0$ 对对角矩阵等价于 $\mathbf{x} \ge \mathbf{0}$。$\blacksquare$

!!! example "例 25.4"
    **最大割问题的 SDP 松弛。**

    图 $G = (V, E)$ 上的最大割（MAX-CUT）问题可以表述为

    $$
    \max \frac{1}{4} \sum_{(i,j) \in E} w_{ij}(1 - x_i x_j), \quad x_i \in \{-1, +1\}.
    $$

    令 $X = \mathbf{x}\mathbf{x}^T$，则约束等价于 $X \succeq 0$，$\operatorname{rank}(X) = 1$，$X_{ii} = 1$。去掉秩一约束得到 SDP 松弛：

    $$
    \max \frac{1}{4} \langle L, X \rangle, \quad \text{s.t.} \quad X_{ii} = 1, \quad X \succeq 0,
    $$

    其中 $L$ 为图的 Laplacian 矩阵。Goemans-Williamson 定理保证此松弛的近似比至少为 $0.878$。

!!! example "例 25.5"
    **矩阵范数最小化。**

    最小化 $\|A\|_{\text{op}}$（算子范数）等价于 SDP：

    $$
    \min t, \quad \text{s.t.} \quad \begin{pmatrix} tI & A \\ A^T & tI \end{pmatrix} \succeq 0.
    $$

    这是因为 Schur 补条件 $tI - A^T(tI)^{-1}A \succeq 0$ 等价于 $t^2 I \succeq A^TA$，即 $t \ge \sigma_{\max}(A)$。

---

## 25.4 矩阵补全

<div class="context-flow" markdown>

**思路**：秩约束(NP-hard) → **核范数**松弛(秩的最紧凸包) = SDP → 非相干条件 + $O(\mu^2 r \log^2 n)$ 个观测 → 精确恢复

**链接**：Ch21 张量分解的矩阵版本 · Netflix 推荐/协同过滤的数学基础

</div>

矩阵补全（Matrix Completion）是从部分观测恢复低秩矩阵的问题。

!!! definition "定义 25.7 (矩阵补全问题)"
    给定 $M \in \mathbb{R}^{m \times n}$ 的部分元素 $\{M_{ij} : (i,j) \in \Omega\}$，**矩阵补全问题**（matrix completion problem）为

    $$
    \text{find } X \in \mathbb{R}^{m \times n}, \quad \operatorname{rank}(X) \le r, \quad X_{ij} = M_{ij} \; \forall (i,j) \in \Omega.
    $$

    由于秩约束非凸，实际中通常用核范数（nuclear norm）松弛。

!!! definition "定义 25.8 (核范数)"
    矩阵 $X \in \mathbb{R}^{m \times n}$ 的 **核范数**（nuclear norm）定义为

    $$
    \|X\|_* = \sum_{i=1}^{\min(m,n)} \sigma_i(X) = \operatorname{tr}\!\left(\sqrt{X^T X}\right),
    $$

    其中 $\sigma_i(X)$ 为奇异值。核范数是秩函数的最紧凸松弛（在算子范数球上）。

!!! theorem "定理 25.8 (核范数最小化的 SDP 表示)"
    核范数最小化

    $$
    \min_{X} \|X\|_*, \quad \text{s.t.} \quad X_{ij} = M_{ij} \; \forall (i,j) \in \Omega
    $$

    等价于 SDP：

    $$
    \min \frac{1}{2}(\operatorname{tr}(W_1) + \operatorname{tr}(W_2)), \quad \text{s.t.} \quad \begin{pmatrix} W_1 & X \\ X^T & W_2 \end{pmatrix} \succeq 0, \quad X_{ij} = M_{ij}.
    $$

??? proof "证明"
    对任意 $X$ 和半正定块矩阵 $\begin{pmatrix} W_1 & X \\ X^T & W_2 \end{pmatrix} \succeq 0$，由 Schur 补条件，$W_1 \succeq XW_2^{-1}X^T$（当 $W_2 \succ 0$ 时）。因此

    $$
    \operatorname{tr}(W_1) + \operatorname{tr}(W_2) \ge \operatorname{tr}(XW_2^{-1}X^T) + \operatorname{tr}(W_2).
    $$

    对 $W_2$ 优化，取 $W_2 = (X^TX)^{1/2}$，则

    $$
    \operatorname{tr}(X(X^TX)^{-1/2}X^T) + \operatorname{tr}((X^TX)^{1/2}) = 2\operatorname{tr}((X^TX)^{1/2}) = 2\|X\|_*.
    $$

    因此 $\frac{1}{2}(\operatorname{tr}(W_1) + \operatorname{tr}(W_2)) \ge \|X\|_*$，且等号可以达到。$\blacksquare$

!!! theorem "定理 25.9 (矩阵补全的精确恢复条件)"
    设 $M \in \mathbb{R}^{m \times n}$ 为秩 $r$ 矩阵，满足 **非相干条件**（incoherence condition）：设 $M = U\Sigma V^T$，

    $$
    \max_i \|U^T \mathbf{e}_i\|^2 \le \frac{\mu r}{m}, \quad \max_j \|V^T \mathbf{e}_j\|^2 \le \frac{\mu r}{n},
    $$

    其中 $\mu \ge 1$ 为非相干参数。若观测集 $\Omega$ 中每个元素以概率 $p$ 独立被观测，且

    $$
    p \ge C \frac{\mu^2 r \log^2(m+n)}{m \wedge n}
    $$

    （$C$ 为绝对常数），则核范数最小化以高概率精确恢复 $M$。

??? proof "证明"
    **证明思路（Candes-Recht 2009, Gross 2011）。** 需要证明 $M$ 是核范数最小化的唯一解。构造对偶证书（dual certificate）$Y$ 满足：$\mathcal{P}_\Omega(Y) = \mathcal{P}_\Omega(\text{sgn}(M))$，$\mathcal{P}_T(Y) = UV^T$（次梯度条件），$\|\mathcal{P}_{T^\perp}(Y)\|_{\text{op}} < 1$。其中 $T$ 为 $M$ 的切空间，$\mathcal{P}_\Omega$ 为观测集投影。

    对偶证书通过 **辅助随机化方法**（golfing scheme）构造：将 $\Omega$ 随机划分为多个子集 $\Omega_1, \ldots, \Omega_L$，逐步修正 $Y$ 以逼近所需性质。每一步的误差以几何速率衰减，非相干条件保证了投影算子 $\mathcal{P}_{\Omega_j}\mathcal{P}_T$ 的良好行为（RIP 性质）。$\blacksquare$

!!! example "例 25.6"
    **Netflix 推荐问题。**

    用户-电影评分矩阵 $M \in \mathbb{R}^{m \times n}$（$m$ 个用户，$n$ 部电影）通常近似低秩（用户偏好由少数潜在因子决定）。观测到的评分构成 $\Omega$。核范数最小化可以从稀疏观测中恢复完整评分矩阵，从而进行推荐。

!!! example "例 25.7"
    **低秩恢复的交替最小化。**

    实际中常用交替最小化替代核范数最小化：设 $X = LR^T$（$L \in \mathbb{R}^{m \times r}$，$R \in \mathbb{R}^{n \times r}$），交替优化

    $$
    L \leftarrow \arg\min_L \sum_{(i,j) \in \Omega} (L_i^T R_j - M_{ij})^2, \quad
    R \leftarrow \arg\min_R \sum_{(i,j) \in \Omega} (L_i^T R_j - M_{ij})^2.
    $$

    每一步为最小二乘问题，计算高效。

---

## 25.5 压缩感知

<div class="context-flow" markdown>

**核心条件**：**RIP**——测量矩阵 $A$ 在稀疏向量上近似保距 → $\delta_{2s} < \sqrt{2}-1$ 时 $\ell_1$ 最小化精确恢复 $s$-稀疏信号

**随机矩阵连接**：高斯随机 $A$ 以 $m = O(s\log(n/s))$ 行满足 RIP（Ch23 集中不等式）→ 远少于 $n$ 次测量即可重建

</div>

压缩感知（Compressed Sensing）利用信号的稀疏性从远少于 Nyquist 采样定理要求的测量中恢复信号。

!!! definition "定义 25.9 (受限等距性质 RIP)"
    矩阵 $A \in \mathbb{R}^{m \times n}$ 满足 $s$ 阶 **受限等距性质**（Restricted Isometry Property, RIP），常数为 $\delta_s \in (0, 1)$，若对所有 $s$-稀疏向量 $\mathbf{x}$（至多 $s$ 个非零分量），

    $$
    (1 - \delta_s) \|\mathbf{x}\|_2^2 \le \|A\mathbf{x}\|_2^2 \le (1 + \delta_s) \|\mathbf{x}\|_2^2.
    $$

    直觉上，RIP 要求 $A$ 在稀疏向量上近似保距。

!!! definition "定义 25.10 (基追踪)"
    **基追踪**（Basis Pursuit, BP）是通过 $\ell_1$ 范数最小化恢复稀疏信号的方法：

    $$
    \min_{\mathbf{x} \in \mathbb{R}^n} \|\mathbf{x}\|_1, \quad \text{s.t.} \quad A\mathbf{x} = \mathbf{b}.
    $$

    含噪声版本（BPDN）：$\min \|\mathbf{x}\|_1$，s.t. $\|A\mathbf{x} - \mathbf{b}\|_2 \le \epsilon$。

!!! theorem "定理 25.10 (RIP 保证精确恢复)"
    设 $A \in \mathbb{R}^{m \times n}$ 满足 $\delta_{2s} < \sqrt{2} - 1$。若 $\mathbf{x}^0$ 为 $s$-稀疏向量，$\mathbf{b} = A\mathbf{x}^0$，则基追踪的解 $\hat{\mathbf{x}}$ 唯一且 $\hat{\mathbf{x}} = \mathbf{x}^0$。

??? proof "证明"
    设 $\hat{\mathbf{x}}$ 为基追踪的解，$\mathbf{h} = \hat{\mathbf{x}} - \mathbf{x}^0$。则 $A\mathbf{h} = 0$，且 $\|\hat{\mathbf{x}}\|_1 \le \|\mathbf{x}^0\|_1$。

    设 $S$ 为 $\mathbf{x}^0$ 的支撑集（$|S| \le s$），$S^c$ 为补集。由三角不等式：

    $$
    \|\mathbf{x}^0 + \mathbf{h}\|_1 = \|\mathbf{x}^0_S + \mathbf{h}_S\|_1 + \|\mathbf{h}_{S^c}\|_1 \ge \|\mathbf{x}^0_S\|_1 - \|\mathbf{h}_S\|_1 + \|\mathbf{h}_{S^c}\|_1.
    $$

    由 $\|\hat{\mathbf{x}}\|_1 \le \|\mathbf{x}^0\|_1 = \|\mathbf{x}^0_S\|_1$，得 $\|\mathbf{h}_{S^c}\|_1 \le \|\mathbf{h}_S\|_1$（**锥约束**）。

    现在将 $S^c$ 上的分量按绝对值从大到小排列，分成大小为 $s$ 的块 $S_1, S_2, \ldots$。利用锥约束和 RIP，可以证明

    $$
    \|\mathbf{h}\|_2^2 \le \frac{2\delta_{2s}}{1 - \delta_{2s}} \|\mathbf{h}_S\|_2 \|\mathbf{h}\|_2 + \text{尾项},
    $$

    当 $\delta_{2s} < \sqrt{2} - 1$ 时，可以推导出 $\|\mathbf{h}\|_2 = 0$，即精确恢复。$\blacksquare$

!!! theorem "定理 25.11 (高斯随机矩阵满足 RIP)"
    设 $A \in \mathbb{R}^{m \times n}$，元素 $A_{ij} \sim N(0, 1/m)$ 独立同分布。若

    $$
    m \ge C \delta^{-2} s \log(n/s),
    $$

    则 $A$ 以高概率满足 $s$ 阶 RIP，常数为 $\delta_s \le \delta$。

??? proof "证明"
    **证明思路。** 固定一个 $s$-稀疏向量 $\mathbf{x}$，$\|A\mathbf{x}\|^2 = \sum_{i=1}^m (\mathbf{a}_i^T\mathbf{x})^2$，其中 $\mathbf{a}_i$ 为 $A$ 的行。$\mathbf{a}_i^T\mathbf{x} \sim N(0, \|\mathbf{x}\|^2/m)$，故 $\|A\mathbf{x}\|^2/\|\mathbf{x}\|^2$ 为 $\chi^2(m)/m$。由集中不等式：

    $$
    P\!\left(\left|\frac{\|A\mathbf{x}\|^2}{\|\mathbf{x}\|^2} - 1\right| > \delta\right) \le 2e^{-cm\delta^2}.
    $$

    对所有 $s$-稀疏向量取并界：$s$-稀疏向量的单位球可以用 $\epsilon$-网覆盖，网的大小为 $(Cn/s)^s$。取 $m \ge C'\delta^{-2}s\log(n/s)$ 使得并界的概率趋于零。$\blacksquare$

!!! example "例 25.8"
    **稀疏信号恢复。**

    设 $\mathbf{x}^0 \in \mathbb{R}^{1000}$ 有 $s = 20$ 个非零分量。取 $m = 200$ 个高斯随机测量 $\mathbf{b} = A\mathbf{x}^0$。由定理 25.11，$m = 200 \ge C \cdot 20 \cdot \log(50) \approx 78C$ 满足 RIP 条件（$C$ 为适度常数）。基追踪可以从这 $200$ 个测量中精确恢复 $1000$ 维信号。

---

## 25.6 主成分分析（PCA）

<div class="context-flow" markdown>

**SVD 即 PCA**：主成分方向 = $S = \frac{1}{n}\bar{X}^T\bar{X}$ 的特征向量 = $\bar{X}$ 的右奇异向量 → Eckart-Young 最优低秩近似

**鲁棒 PCA**：$M = L + S$（低秩+稀疏）→ $\|L\|_* + \lambda\|S\|_1$ 凸松弛 → 视频前景/背景分离 · 链接 Ch23 BBP 相变（何时能检测到信号？）

</div>

主成分分析是数据降维和特征提取的经典方法，其数学基础完全建立在线性代数之上。

!!! definition "定义 25.11 (主成分分析)"
    给定数据矩阵 $X \in \mathbb{R}^{n \times p}$（$n$ 个样本，$p$ 个特征），中心化后 $\bar{X} = X - \frac{1}{n}\mathbf{1}\mathbf{1}^TX$。**主成分分析**（Principal Component Analysis, PCA）寻找正交方向 $\mathbf{v}_1, \ldots, \mathbf{v}_k$ 以最大化投影方差：

    $$
    \mathbf{v}_1 = \arg\max_{\|\mathbf{v}\|=1} \frac{1}{n}\|\bar{X}\mathbf{v}\|^2 = \arg\max_{\|\mathbf{v}\|=1} \mathbf{v}^T S \mathbf{v},
    $$

    其中 $S = \frac{1}{n}\bar{X}^T\bar{X}$ 为样本协方差矩阵。

!!! definition "定义 25.12 (鲁棒 PCA)"
    **鲁棒 PCA**（Robust PCA）将数据矩阵分解为低秩部分和稀疏部分：

    $$
    \min_{L, S} \|L\|_* + \lambda \|S\|_1, \quad \text{s.t.} \quad L + S = M,
    $$

    其中 $\|L\|_*$ 为核范数（促进低秩），$\|S\|_1 = \sum_{ij}|S_{ij}|$（促进稀疏），$\lambda > 0$ 为权衡参数。

!!! theorem "定理 25.12 (PCA 的 SVD 推导)"
    设 $\bar{X} = U\Sigma V^T$ 为 SVD。则：

    1. 第 $k$ 个主成分方向 $\mathbf{v}_k$ 为 $V$ 的第 $k$ 列（$S$ 的第 $k$ 个特征向量）；
    2. 第 $k$ 个主成分得分为 $\bar{X}\mathbf{v}_k = \sigma_k \mathbf{u}_k$；
    3. 前 $k$ 个主成分解释的方差比例为 $\sum_{i=1}^k \sigma_i^2 / \sum_{i=1}^p \sigma_i^2$；
    4. $\bar{X}$ 的最佳秩 $k$ 近似（在 Frobenius 范数下）为 $\bar{X}_k = U_k \Sigma_k V_k^T = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$。

??? proof "证明"
    $S = \frac{1}{n}\bar{X}^T\bar{X} = \frac{1}{n}V\Sigma^2 V^T$，故 $S$ 的特征值为 $\sigma_i^2/n$，特征向量为 $V$ 的列。

    第一主成分：$\max_{\|\mathbf{v}\|=1} \mathbf{v}^T S \mathbf{v} = \sigma_1^2/n$，达到于 $\mathbf{v}_1$。

    投影方差：$\frac{1}{n}\|\bar{X}\mathbf{v}_k\|^2 = \frac{1}{n}\sigma_k^2$。

    最佳秩 $k$ 近似：由 Eckart-Young-Mirsky 定理，$\min_{\operatorname{rank}(B)\le k}\|\bar{X} - B\|_F = \sqrt{\sum_{i=k+1}^p \sigma_i^2}$，最优 $B = \bar{X}_k$。$\blacksquare$

!!! theorem "定理 25.13 (鲁棒 PCA 的精确恢复)"
    设 $M = L_0 + S_0$，$L_0$ 秩 $r$（满足非相干条件），$S_0$ 稀疏（非零元素比例 $\rho$）。取 $\lambda = 1/\sqrt{\max(m,n)}$，则在 $r$ 足够小、$\rho$ 足够小的条件下，鲁棒 PCA 的解 $(\hat{L}, \hat{S})$ 以高概率满足 $\hat{L} = L_0$，$\hat{S} = S_0$。

??? proof "证明"
    **证明思路（Candes-Li-Ma-Wright 2011）。** 构造次梯度条件的对偶证书。设 $L_0 = U\Sigma V^T$，需要构造矩阵 $W$ 满足：

    $$
    W = UV^T + D, \quad \|W\|_{\text{op}} \le 1, \quad \mathcal{P}_\Omega(W) = \lambda \operatorname{sgn}(S_0), \quad \|\mathcal{P}_{\Omega^c}(W - UV^T)\|_\infty < \lambda,
    $$

    其中 $\Omega$ 为 $S_0$ 的支撑。利用交替投影方法（类似矩阵补全中的 golfing scheme）和非相干条件的分析可以完成构造。$\blacksquare$

!!! example "例 25.9"
    **视频监控中的前景-背景分离。**

    视频序列的每一帧展平为列向量，组成矩阵 $M$。静态背景对应低秩矩阵 $L$（背景帧高度相关），运动前景对应稀疏矩阵 $S$（仅少量像素变化）。鲁棒 PCA $M = L + S$ 自动完成前景-背景分离。

---

## 25.7 低秩近似与降维

<div class="context-flow" markdown>

**两大定理**：**Eckart-Young-Mirsky**（截断 SVD = 最优低秩近似）+ **JL 引理**（随机投影 $\mathbb{R}^n \to \mathbb{R}^{O(\log N/\epsilon^2)}$ 近似保距）

**链接**：Ch23 高斯随机矩阵的集中性质是 JL 引理的证明核心

</div>

低秩近似和降维是处理高维数据的核心技术。

!!! definition "定义 25.13 (Johnson-Lindenstrauss 变换)"
    **Johnson-Lindenstrauss (JL) 变换**是从 $\mathbb{R}^n$ 到 $\mathbb{R}^m$（$m \ll n$）的随机线性映射 $\Phi: \mathbb{R}^n \to \mathbb{R}^m$，使得高维点集的两两距离在低维空间中近似保持。

!!! theorem "定理 25.14 (Johnson-Lindenstrauss 引理)"
    对任意 $\epsilon \in (0, 1)$ 和 $\mathbb{R}^n$ 中的 $N$ 个点 $\mathbf{x}_1, \ldots, \mathbf{x}_N$，存在映射 $f: \mathbb{R}^n \to \mathbb{R}^m$，$m = O(\epsilon^{-2} \log N)$，使得

    $$
    (1 - \epsilon)\|\mathbf{x}_i - \mathbf{x}_j\|_2^2 \le \|f(\mathbf{x}_i) - f(\mathbf{x}_j)\|_2^2 \le (1 + \epsilon)\|\mathbf{x}_i - \mathbf{x}_j\|_2^2
    $$

    对所有 $i, j$ 成立。特别地，$f$ 可以取为随机矩阵 $\Phi \in \mathbb{R}^{m \times n}$（$\Phi_{ij} \sim N(0, 1/m)$），$f(\mathbf{x}) = \Phi\mathbf{x}$。

??? proof "证明"
    **第一步：单个向量的集中。** 设 $\mathbf{x} \in \mathbb{R}^n$，$\|\mathbf{x}\| = 1$，$\Phi_{ij} \sim N(0, 1/m)$ 独立。则 $\|\Phi\mathbf{x}\|^2 = \sum_{i=1}^m (\boldsymbol{\phi}_i^T\mathbf{x})^2$，其中 $\boldsymbol{\phi}_i^T\mathbf{x} \sim N(0, 1/m)$。因此 $m\|\Phi\mathbf{x}\|^2 \sim \chi^2(m)$。由 $\chi^2$ 分布的集中不等式：

    $$
    P\!\left(\left|\|\Phi\mathbf{x}\|^2 - 1\right| > \epsilon\right) \le 2\exp\!\left(-\frac{m\epsilon^2}{8}\right).
    $$

    **第二步：并界。** 对 $N$ 个点的 $\binom{N}{2}$ 对差向量 $\mathbf{x}_i - \mathbf{x}_j$，取并界：

    $$
    P(\exists \, i,j : \text{距离失真} > \epsilon) \le 2\binom{N}{2}\exp\!\left(-\frac{m\epsilon^2}{8}\right) \le N^2 \exp\!\left(-\frac{m\epsilon^2}{8}\right).
    $$

    取 $m \ge \frac{16}{\epsilon^2}\ln N$，则失败概率 $\le N^2 \cdot N^{-2} = 1$（实际取更大常数使概率趋于零）。$\blacksquare$

!!! theorem "定理 25.15 (Eckart-Young-Mirsky 定理)"
    设 $A \in \mathbb{R}^{m \times n}$，$\operatorname{rank}(A) = r$，SVD 为 $A = \sum_{i=1}^r \sigma_i \mathbf{u}_i \mathbf{v}_i^T$。则对 $k < r$，

    $$
    \min_{\operatorname{rank}(B) \le k} \|A - B\|_F = \sqrt{\sum_{i=k+1}^r \sigma_i^2}, \quad
    \min_{\operatorname{rank}(B) \le k} \|A - B\|_{\text{op}} = \sigma_{k+1},
    $$

    最优解为截断 SVD $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$。

??? proof "证明"
    **Frobenius 范数情形。** 设 $B$ 为秩 $\le k$ 矩阵。$\|A - B\|_F^2 = \operatorname{tr}((A-B)^T(A-B))$。设 $B$ 的列空间为 $\mathcal{V}$（$\dim \mathcal{V} \le k$），行空间为 $\mathcal{W}$（$\dim \mathcal{W} \le k$）。

    由 Courant-Fischer 极小极大原理，$\mathcal{V}^\perp$ 是 $n - k$ 维子空间。$A$ 在 $\mathcal{V}^\perp$ 上的限制满足

    $$
    \|A - B\|_F^2 \ge \sum_{i=k+1}^r \sigma_i^2(A|_{\mathcal{V}^\perp}) \ge \sum_{i=k+1}^r \sigma_i^2(A),
    $$

    后一个不等式由奇异值的极小极大刻画得出。等号在 $B = A_k$ 时达到。

    **算子范数情形。** $\ker(B)$ 至少 $n - k$ 维。取 $\mathbf{v} \in \ker(B) \cap \operatorname{span}(\mathbf{v}_1, \ldots, \mathbf{v}_{k+1})$（维数论证保证此交集非零），$\|\mathbf{v}\| = 1$，则

    $$
    \|A - B\|_{\text{op}} \ge \|(A - B)\mathbf{v}\| = \|A\mathbf{v}\| \ge \sigma_{k+1}.
    $$

    等号在 $B = A_k$ 时达到。$\blacksquare$

!!! example "例 25.10"
    **JL 引理的应用：近似最近邻搜索。**

    在 $\mathbb{R}^{10000}$ 中有 $N = 10^6$ 个数据点。JL 引理保证可以将它们投影到 $m = O(\epsilon^{-2} \cdot 6\log 10) \approx O(\epsilon^{-2} \cdot 14)$ 维空间，保持距离的 $1 \pm \epsilon$ 因子。取 $\epsilon = 0.1$，$m \approx 1400$，显著降低了最近邻搜索的计算复杂度。

!!! example "例 25.11"
    **图像压缩中的低秩近似。**

    灰度图像可以表示为矩阵 $A \in \mathbb{R}^{m \times n}$。若 $A$ 有快速衰减的奇异值，秩 $k$ 截断 SVD $A_k$ 用 $k(m + n + 1)$ 个参数近似 $mn$ 个像素值。例如 $m = n = 1000$，$k = 50$ 时，压缩比约为 $1000000 / 100050 \approx 10:1$。

---

## 25.8 特征值优化

<div class="context-flow" markdown>

**特征值 = 极值**：Rayleigh 商 $R_A(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}/\|\mathbf{x}\|^2$ → $\lambda_{\min} \le R_A \le \lambda_{\max}$

**Courant-Fischer** 极小极大 → **Weyl 扰动不等式** $|\gamma_i - \alpha_i| \le \|B\|$

**汇聚**：Rayleigh 商迭代(三次收敛, Ch22) · 图 Laplacian 的 $\lambda_2$ = 连通性(Fiedler) · 链接 Ch24 Stiefel 流形上的特征值问题

</div>

许多优化问题的解可以通过特征值来刻画。

!!! definition "定义 25.14 (Rayleigh 商)"
    设 $A \in \operatorname{Sym}(n)$。**Rayleigh 商**（Rayleigh quotient）定义为

    $$
    R_A(\mathbf{x}) = \frac{\mathbf{x}^T A \mathbf{x}}{\mathbf{x}^T \mathbf{x}}, \quad \mathbf{x} \ne \mathbf{0}.
    $$

    Rayleigh 商是齐零次函数，因此可以限制在单位球面 $\|\mathbf{x}\| = 1$ 上。

!!! definition "定义 25.15 (广义 Rayleigh 商)"
    设 $A, B \in \operatorname{Sym}(n)$，$B \succ 0$。**广义 Rayleigh 商**（generalized Rayleigh quotient）为

    $$
    R_{A,B}(\mathbf{x}) = \frac{\mathbf{x}^T A \mathbf{x}}{\mathbf{x}^T B \mathbf{x}}, \quad \mathbf{x} \ne \mathbf{0}.
    $$

    其极值与广义特征值问题 $A\mathbf{x} = \lambda B\mathbf{x}$ 相关。

!!! theorem "定理 25.16 (Rayleigh 商的极值)"
    设 $A \in \operatorname{Sym}(n)$，特征值 $\lambda_1 \le \lambda_2 \le \cdots \le \lambda_n$，对应特征向量 $\mathbf{v}_1, \ldots, \mathbf{v}_n$。则

    $$
    \lambda_1 = \min_{\mathbf{x} \ne \mathbf{0}} R_A(\mathbf{x}), \quad \lambda_n = \max_{\mathbf{x} \ne \mathbf{0}} R_A(\mathbf{x}),
    $$

    极值分别在 $\mathbf{x} = \mathbf{v}_1$ 和 $\mathbf{x} = \mathbf{v}_n$ 处达到。

??? proof "证明"
    令 $\mathbf{x} = \sum_{i=1}^n c_i \mathbf{v}_i$，$\|\mathbf{x}\| = 1$，即 $\sum c_i^2 = 1$。则

    $$
    R_A(\mathbf{x}) = \mathbf{x}^T A \mathbf{x} = \sum_{i=1}^n \lambda_i c_i^2.
    $$

    因为 $\sum c_i^2 = 1$ 且 $\lambda_1 \le \cdots \le \lambda_n$，

    $$
    \lambda_1 = \lambda_1 \sum c_i^2 \le \sum \lambda_i c_i^2 \le \lambda_n \sum c_i^2 = \lambda_n.
    $$

    下界在 $c_1 = 1$（即 $\mathbf{x} = \mathbf{v}_1$）时达到，上界在 $c_n = 1$ 时达到。$\blacksquare$

<div class="context-flow" markdown>

**洞察**：Courant-Fischer 将特征值从"代数对象"($\det(A-\lambda I)=0$)变为"优化对象"(子空间上的极值)——由此推导 Weyl 扰动、Lidskii 不等式、图分割等一切特征值不等式

</div>

!!! theorem "定理 25.17 (Courant-Fischer 极小极大定理)"
    设 $A \in \operatorname{Sym}(n)$，特征值 $\lambda_1 \le \lambda_2 \le \cdots \le \lambda_n$。则

    $$
    \lambda_k = \min_{\dim V = k} \max_{\mathbf{x} \in V, \|\mathbf{x}\|=1} \mathbf{x}^T A \mathbf{x} = \max_{\dim W = n-k+1} \min_{\mathbf{x} \in W, \|\mathbf{x}\|=1} \mathbf{x}^T A \mathbf{x}.
    $$

??? proof "证明"
    **第一等式的证明。** 令 $V_k = \operatorname{span}(\mathbf{v}_1, \ldots, \mathbf{v}_k)$。

    **上界**：取 $V = V_k$，则 $\max_{\mathbf{x} \in V_k, \|\mathbf{x}\|=1} \mathbf{x}^TA\mathbf{x} = \lambda_k$（因为 $A$ 限制在 $V_k$ 上的最大特征值为 $\lambda_k$）。

    **下界**：设 $V$ 为任意 $k$ 维子空间。考虑 $W = \operatorname{span}(\mathbf{v}_k, \ldots, \mathbf{v}_n)$，$\dim W = n - k + 1$。由维数公式，$\dim(V \cap W) \ge k + (n-k+1) - n = 1$。取 $\mathbf{x} \in V \cap W$，$\|\mathbf{x}\| = 1$，则

    $$
    \max_{\mathbf{y} \in V, \|\mathbf{y}\|=1} \mathbf{y}^T A \mathbf{y} \ge \mathbf{x}^T A \mathbf{x} \ge \lambda_k,
    $$

    后一个不等式因为 $\mathbf{x} \in W$ 蕴含 $\mathbf{x}^T A \mathbf{x} \ge \lambda_k$（$A$ 在 $W$ 上的最小特征值为 $\lambda_k$）。

    综合得 $\lambda_k = \min_{\dim V = k} \max_{\mathbf{x} \in V} \mathbf{x}^TA\mathbf{x}$。第二等式的证明类似。$\blacksquare$

!!! theorem "定理 25.18 (Weyl 特征值扰动不等式)"
    设 $A, B \in \operatorname{Sym}(n)$，特征值分别为 $\alpha_1 \le \cdots \le \alpha_n$ 和 $\beta_1 \le \cdots \le \beta_n$，$A + B$ 的特征值为 $\gamma_1 \le \cdots \le \gamma_n$。则

    $$
    \alpha_i + \beta_1 \le \gamma_i \le \alpha_i + \beta_n, \quad i = 1, \ldots, n.
    $$

    特别地，$|\gamma_i - \alpha_i| \le \|B\|_{\text{op}}$。

??? proof "证明"
    由 Courant-Fischer 定理，

    $$
    \gamma_k = \min_{\dim V = k} \max_{\mathbf{x} \in V, \|\mathbf{x}\|=1} \mathbf{x}^T(A+B)\mathbf{x} \le \min_{\dim V = k} \max_{\mathbf{x} \in V, \|\mathbf{x}\|=1} (\mathbf{x}^TA\mathbf{x} + \beta_n) = \alpha_k + \beta_n.
    $$

    类似地，$\gamma_k \ge \alpha_k + \beta_1$。取 $B' = -B$ 可以对称地得到 $|\gamma_k - \alpha_k| \le \max(|\beta_1|, |\beta_n|) = \|B\|_{\text{op}}$。$\blacksquare$

!!! theorem "定理 25.19 (Lidskii 不等式)"
    设 $A, B \in \operatorname{Sym}(n)$，特征值分别为 $\alpha_1 \le \cdots \le \alpha_n$ 和 $\beta_1 \le \cdots \le \beta_n$。设 $C = A + B$，特征值 $\gamma_1 \le \cdots \le \gamma_n$。则对任意 $1 \le i_1 < i_2 < \cdots < i_k \le n$，

    $$
    \sum_{j=1}^{k} \gamma_{i_j} \le \sum_{j=1}^{k} \alpha_{i_j} + \sum_{j=1}^{k} \beta_{n-k+j}.
    $$

??? proof "证明"
    利用 Courant-Fischer 定理的推广。对指标集 $I = \{i_1, \ldots, i_k\}$，构造适当的子空间链 $V_1 \subset V_2 \subset \cdots \subset V_k$，使得 $\dim V_j = i_j$。利用 Courant-Fischer 刻画每个 $\gamma_{i_j}$，再通过子空间相交的维数论证，将 $\gamma_{i_j}$ 的 Rayleigh 商分解为 $A$ 和 $B$ 的贡献，从而得到不等式。完整证明可参见 Bhatia《Matrix Analysis》定理 III.4.1。$\blacksquare$

!!! example "例 25.12"
    **Rayleigh 商迭代。**

    求 $A \in \operatorname{Sym}(n)$ 的特征值和特征向量：给定初始 $\mathbf{x}_0$，$\|\mathbf{x}_0\| = 1$，迭代

    $$
    \rho_k = \mathbf{x}_k^T A \mathbf{x}_k, \quad (A - \rho_k I)\mathbf{y}_{k+1} = \mathbf{x}_k, \quad \mathbf{x}_{k+1} = \mathbf{y}_{k+1}/\|\mathbf{y}_{k+1}\|.
    $$

    Rayleigh 商迭代具有三次收敛速度：$|\rho_{k+1} - \lambda| = O(|\rho_k - \lambda|^3)$，是求特征值最快的迭代方法之一。

!!! example "例 25.13"
    **Courant-Fischer 定理在图论中的应用。**

    设 $L$ 为图 $G$ 的 Laplacian 矩阵。由 Courant-Fischer 定理，

    $$
    \lambda_2(L) = \min_{\mathbf{x} \perp \mathbf{1}, \|\mathbf{x}\|=1} \mathbf{x}^T L \mathbf{x} = \min_{\mathbf{x} \perp \mathbf{1}} \frac{\sum_{(i,j)\in E}(x_i - x_j)^2}{\sum_i x_i^2}.
    $$

    $\lambda_2(L)$（Fiedler 值）衡量了图的连通性：$\lambda_2 > 0$ 当且仅当 $G$ 连通，$\lambda_2$ 越大图的连通性越强。Fiedler 向量（$\lambda_2$ 对应的特征向量）用于图分割。

!!! example "例 25.14"
    **Weyl 不等式的应用：特征值的稳定性。**

    设 $A$ 为对称矩阵，$E$ 为对称扰动，$\|E\|_{\text{op}} = \epsilon$。由 Weyl 不等式，$A + E$ 的第 $k$ 个特征值与 $A$ 的第 $k$ 个特征值之差不超过 $\epsilon$。因此，若要求特征值精度 $10^{-6}$，只需将矩阵元素精度控制在 $O(10^{-6})$ 量级。这为数值特征值计算提供了理论保障。

---

## 本章小结

本章展示了线性代数在优化理论与算法中的核心地位：

1. **线性规划**的单纯形法本质上是线性方程组求解的迭代过程，基本可行解对应于列选取。
2. **最小二乘**问题通过正规方程、QR 分解或 SVD 求解，Tikhonov 正则化通过谱过滤抑制噪声放大。
3. **半定规划**将线性规划推广到矩阵空间，强对偶性是其理论基石。
4. **矩阵补全**利用核范数松弛将 NP 难的秩约束问题转化为凸优化。
5. **压缩感知**中，RIP 条件（线性代数的受限等距性质）保证了 $\ell_1$ 最小化的精确恢复。
6. **PCA** 通过 SVD 实现最优降维，鲁棒 PCA 分离低秩和稀疏成分。
7. **JL 引理**和 **Eckart-Young-Mirsky 定理**为维度约减提供了理论基础。
8. **Courant-Fischer 定理**和 **Weyl 不等式**是特征值优化和扰动分析的基本工具。
