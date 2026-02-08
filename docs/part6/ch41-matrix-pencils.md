# 第 41 章 矩阵束与 Kronecker 标准形

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · Jordan 标准形(Ch12) · $\lambda$-矩阵(Ch13B)

**本章脉络**：矩阵束定义 $\to$ 正则束与广义特征值 $\to$ 奇异束 $\to$ Kronecker 标准形 $\to$ QZ 算法 $\to$ 微分代数方程(DAE) $\to$ Rosenbrock 系统矩阵

**延伸**：矩阵束理论是微分代数方程（电路模拟中的指标概念）和控制理论（描述子系统）的数学基础；QZ 算法是 LAPACK 中广义特征值问题的标准求解器

</div>

标准特征值问题 $Ax = \lambda x$ 可以视为矩阵束 $A - \lambda I$ 的研究。当将 $I$ 替换为一般矩阵 $B$ 时，就得到了**广义特征值问题** $Ax = \lambda Bx$ 和**矩阵束** $A - \lambda B$。当 $B$ 奇异时，矩阵束可能是**奇异的**——它没有明确定义的特征多项式，其标准形理论远比 Jordan 标准形复杂。Kronecker 在 1890 年建立的标准形理论完整地分类了所有矩阵束，是矩阵论中最精致也是最困难的结果之一。

本章从矩阵束的基本概念出发，经过正则束的 Weierstrass 标准形和奇异束的 Kronecker 块，到达完整的 Kronecker 标准形。然后讨论数值计算中的 QZ 算法以及在微分代数方程和控制理论中的应用。

---

## 41.1 矩阵束的定义

<div class="context-flow" markdown>

**核心问题**：什么是矩阵束？正则束与奇异束有何区别？

</div>

!!! definition "定义 41.1 (矩阵束)"
    设 $A, B \in M_{m \times n}(\mathbb{C})$。**矩阵束**（matrix pencil）$L(\lambda)$ 定义为
    $$L(\lambda) = A - \lambda B,$$
    其中 $\lambda \in \mathbb{C}$ 为参数。当 $m = n$（方阵情形）时，可以考虑 $\det(A - \lambda B)$ 作为 $\lambda$ 的多项式。

!!! definition "定义 41.2 (正则束与奇异束)"
    设 $A, B \in M_n(\mathbb{C})$（方阵）。矩阵束 $A - \lambda B$ 称为**正则的**（regular），若
    $$\det(A - \lambda B) \not\equiv 0$$
    （即 $\det(A - \lambda B)$ 作为 $\lambda$ 的多项式不恒等于零）。否则称为**奇异的**（singular）。

    对于长方矩阵束（$m \ne n$），束总是奇异的（因为没有方阵行列式）。

!!! example "例 41.1"
    (a) $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$，$B = I$：$\det(A - \lambda I) = (1-\lambda)(2-\lambda)$，这是标准特征值问题，正则束。

    (b) $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$：$\det(A - \lambda B) = \det\begin{pmatrix} 1 & 0 \\ 0 & -\lambda \end{pmatrix} = -\lambda$，正则束。

    (c) $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$：$\det(A - \lambda B) = 0$ 对所有 $\lambda$，奇异束。

    (d) $A = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 1 \\ 0 & 0 \\ 0 & 0 \end{pmatrix}$：$3 \times 2$ 长方束，奇异。

!!! definition "定义 41.3 (广义特征值)"
    设 $A - \lambda B$ 为 $n \times n$ 正则束。$\lambda_0 \in \mathbb{C}$ 称为 $A - \lambda B$ 的**有限广义特征值**，若 $\det(A - \lambda_0 B) = 0$。若 $\det(B) = 0$ 但 $\det(A - \lambda B) \not\equiv 0$，则称束有**无穷广义特征值**。

    更精确地说，可以使用齐次坐标 $(\alpha, \beta)$：$(\alpha_0, \beta_0)$ 为广义特征对，若 $\det(\beta_0 A - \alpha_0 B) = 0$。有限特征值对应 $\beta_0 \ne 0$（$\lambda = \alpha_0/\beta_0$），无穷特征值对应 $\beta_0 = 0$（$\alpha_0 \ne 0$）。

!!! definition "定义 41.4 (严格等价)"
    两个 $m \times n$ 矩阵束 $A_1 - \lambda B_1$ 和 $A_2 - \lambda B_2$ 称为**严格等价的**（strictly equivalent），若存在非奇异矩阵 $P \in GL_m(\mathbb{C})$，$Q \in GL_n(\mathbb{C})$，使得
    $$P(A_1 - \lambda B_1)Q = A_2 - \lambda B_2,$$
    即 $PA_1 Q = A_2$ 且 $PB_1 Q = B_2$。

---

## 41.2 正则束与广义特征值

<div class="context-flow" markdown>

**核心问题**：正则束的标准形是什么？如何处理无穷特征值？

</div>

!!! theorem "定理 41.1 (Weierstrass 标准形)"
    设 $A - \lambda B$ 为 $n \times n$ 正则束。则存在非奇异矩阵 $P, Q \in GL_n(\mathbb{C})$ 使得
    $$P(A - \lambda B)Q = \begin{pmatrix} J - \lambda I_r & 0 \\ 0 & I_s - \lambda N \end{pmatrix},$$
    其中：

    - $J \in M_r(\mathbb{C})$ 为 Jordan 标准形，对角块对应有限广义特征值 $\lambda_1, \ldots, \lambda_p$；
    - $N \in M_s(\mathbb{C})$ 为**幂零** Jordan 矩阵（所有对角元为 0），对应无穷广义特征值；
    - $r + s = n$。

    此标准形在块的排列顺序外是唯一的。

??? proof "证明（概要）"
    **步骤一**：由 $A - \lambda B$ 正则，$\det(A - \lambda B)$ 不恒为零。设 $d = \deg(\det(A-\lambda B)) = r$，则 $A - \lambda B$ 有 $r$ 个有限特征值（计重数）。

    **步骤二**：选取 $\mu \notin \sigma(A,B)$（不是广义特征值），则 $A - \mu B$ 非奇异。令 $C = (A - \mu B)^{-1}B$。标准特征值问题 $Cx = \lambda' x$ 的特征值 $\lambda'$ 与广义特征值 $\lambda$ 的关系为 $\lambda' = 1/(\lambda - \mu)$（有限特征值）或 $\lambda' = 0$（$\lambda = \infty$）。

    **步骤三**：$C$ 的 Jordan 形为 $C = S \begin{pmatrix} J_f & 0 \\ 0 & N \end{pmatrix} S^{-1}$，其中 $J_f$ 对应非零特征值（有限广义特征值），$N$ 为幂零部分（无穷广义特征值）。

    **步骤四**：通过适当选取 $P, Q$（利用 $S$ 和 $(A-\mu B)^{-1}$），将 $A - \lambda B$ 化为 Weierstrass 形式。

!!! example "例 41.2"
    矩阵束 $A - \lambda B$，$A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$。

    $\det(A - \lambda B) = \det\begin{pmatrix} 1 & 0 & 0 \\ 0 & -\lambda & 1 \\ 0 & 0 & -\lambda \end{pmatrix} = \lambda^2$。

    有限特征值 $\lambda = 0$（代数重数 2）。$\deg(\det) = 2 < n = 3$，故有 1 个无穷特征值。

    Weierstrass 形式：$J = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$（$\lambda = 0$ 的 $2 \times 2$ Jordan 块），$N = (0)$（$1 \times 1$ 幂零块）。

    $$\text{Weierstrass 形} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix} - \lambda \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}.$$

!!! theorem "定理 41.2 (正则束的 Smith 标准形)"
    正则束 $A - \lambda B$ 的**Smith 标准形**为
    $$\operatorname{diag}(d_1(\lambda), d_2(\lambda), \ldots, d_r(\lambda), 0, \ldots, 0),$$
    其中 $d_i(\lambda)$ 为首一多项式，$d_i | d_{i+1}$，$r = \operatorname{rank}(A - \lambda B)$（作为 $\mathbb{C}(\lambda)$ 上矩阵的秩）。$d_i$ 称为**不变因子**。

---

## 41.3 奇异束与最小指标

<div class="context-flow" markdown>

**核心问题**：奇异束的结构如何描述？什么是"最小指标"？

</div>

!!! definition "定义 41.5 (列最小指标)"
    设 $L(\lambda) = A - \lambda B$ 为 $m \times n$ 矩阵束（$m \le n$ 或 $m > n$），束是奇异的。$L(\lambda)$ 的**右（列）零空间**定义为
    $$\mathcal{N}_r(\lambda) = \{x(\lambda) \in \mathbb{C}[\lambda]^n : L(\lambda) x(\lambda) = 0\}.$$
    即所有使 $L(\lambda)x(\lambda) = 0$ 的多项式向量 $x(\lambda)$。

    $\mathcal{N}_r(\lambda)$ 中非零元素的最低次数的向量称为**最小基**向量。右最小指标 $\varepsilon_1 \le \varepsilon_2 \le \cdots \le \varepsilon_p$ 是最小基中各基向量的次数。

!!! definition "定义 41.6 (行最小指标)"
    类似地，$L(\lambda)$ 的**左（行）零空间**为
    $$\mathcal{N}_l(\lambda) = \{y(\lambda)^T \in \mathbb{C}[\lambda]^{1 \times m} : y(\lambda)^T L(\lambda) = 0\}.$$
    左最小指标 $\eta_1 \le \eta_2 \le \cdots \le \eta_q$ 是左最小基中各基向量的次数。

!!! example "例 41.3"
    $L(\lambda) = \begin{pmatrix} 1 & -\lambda & 0 \\ 0 & 1 & -\lambda \end{pmatrix}$（$2 \times 3$ 束）。

    右零空间：解 $L(\lambda)x(\lambda) = 0$：$x_1 = \lambda x_2$，$x_2 = \lambda x_3$。取 $x_3 = 1$，得 $x(\lambda) = (\lambda^2, \lambda, 1)^T$。最低次数为 2，右最小指标 $\varepsilon_1 = 2$。

    左零空间：解 $y^T L = 0$：$y_1 = 0$，$-\lambda y_1 + y_2 = 0$，$-\lambda y_2 = 0$。故 $y = 0$，没有非平凡左零向量。左最小指标为空集（$q = 0$）。

!!! theorem "定理 41.3 (最小指标的维数关系)"
    设 $L(\lambda)$ 为 $m \times n$ 矩阵束，$r = \operatorname{rank}(L(\lambda))$（作为 $\mathbb{C}(\lambda)$ 上矩阵的秩）。则右最小指标的个数 $p = n - r$，左最小指标的个数 $q = m - r$。

---

## 41.4 Kronecker 标准形

<div class="context-flow" markdown>

**核心问题**：任意矩阵束在严格等价下的完全标准形是什么？

</div>

!!! definition "定义 41.7 (Kronecker 块)"
    定义以下标准块：

    **(a) 右 Kronecker 块** $L_\varepsilon$（$\varepsilon \times (\varepsilon+1)$）：
    $$L_\varepsilon = \begin{pmatrix}
    1 & -\lambda & & & \\
    & 1 & -\lambda & & \\
    & & \ddots & \ddots & \\
    & & & 1 & -\lambda
    \end{pmatrix} \in M_{\varepsilon \times (\varepsilon+1)}.$$
    特别地，$L_0 = \begin{pmatrix} 1 \end{pmatrix}$（$0 \times 1$ 没有意义——约定 $L_0$ 为 $1 \times 2$ 矩阵 $(1, -\lambda)$）。

    修正：$L_\varepsilon$ 的尺寸为 $\varepsilon \times (\varepsilon + 1)$（$\varepsilon \ge 1$），当 $\varepsilon = 0$ 时 $L_0$ 为空块（不出现）。

    **(b) 左 Kronecker 块** $L_\eta^T$（$(\eta+1) \times \eta$）：$L_\eta$ 的转置。

    **(c) Jordan 块** $J_k(\lambda_0) - \lambda I_k$：对应有限特征值 $\lambda_0$。

    **(d) 无穷 Jordan 块** $I_k - \lambda N_k$：其中 $N_k$ 为 $k \times k$ 幂零 Jordan 块。

!!! theorem "定理 41.4 (Kronecker 标准形)"
    任何 $m \times n$ 矩阵束 $A - \lambda B$ 在严格等价下可以化为**唯一的**（块的顺序除外）标准形：
    $$P(A - \lambda B)Q = \operatorname{diag}(L_{\varepsilon_1}, \ldots, L_{\varepsilon_p}, L_{\eta_1}^T, \ldots, L_{\eta_q}^T, J_{k_1}(\lambda_1) - \lambda I, \ldots, J_{k_s}(\lambda_s) - \lambda I, I_{l_1} - \lambda N_{l_1}, \ldots, I_{l_t} - \lambda N_{l_t}),$$
    其中：

    - $L_{\varepsilon_i}$：右 Kronecker 块，$\varepsilon_1 \le \cdots \le \varepsilon_p$（右最小指标）；
    - $L_{\eta_j}^T$：左 Kronecker 块，$\eta_1 \le \cdots \le \eta_q$（左最小指标）；
    - $J_{k_i}(\lambda_i) - \lambda I$：有限特征值的 Jordan 块；
    - $I_{l_j} - \lambda N_{l_j}$：无穷特征值的 Jordan 块。

    矩阵束的完全不变量为：右最小指标 $\{\varepsilon_i\}$，左最小指标 $\{\eta_j\}$，有限特征值 Jordan 结构 $\{(\lambda_i, k_i)\}$，无穷特征值 Jordan 结构 $\{l_j\}$。

??? proof "证明（概要）"
    证明分为几个阶段：

    1. **分离奇异部分和正则部分**：通过行列变换，将束化为
    $$\begin{pmatrix} L_{\text{sing}} & 0 \\ 0 & L_{\text{reg}} \end{pmatrix},$$
    其中 $L_{\text{reg}}$ 是正则束，$L_{\text{sing}}$ 是奇异束（$\det \equiv 0$ 或非方阵）。

    2. **正则部分的 Weierstrass 分解**：$L_{\text{reg}}$ 由定理 41.1 化为 Jordan 块和幂零块。

    3. **奇异部分的 Kronecker 分解**：关键步骤。利用右最小基和左最小基，将 $L_{\text{sing}}$ 化为右 Kronecker 块 $L_\varepsilon$ 和左 Kronecker 块 $L_\eta^T$ 的直和。

    唯一性的证明需要证明所有不变量（最小指标、Jordan 结构）在严格等价下不变。这可以通过分析 $L(\lambda)$ 的 Smith 形式和最小基的代数性质来完成。

!!! example "例 41.4"
    $L(\lambda) = \begin{pmatrix} 1 & -\lambda & 0 & 0 \\ 0 & 0 & 1 & -\lambda \end{pmatrix}$（$2 \times 4$ 束）。

    这是两个 $L_1$ 块的直和：$L_1 = \begin{pmatrix} 1 & -\lambda \end{pmatrix}$（$1 \times 2$）。

    Kronecker 标准形为 $\operatorname{diag}(L_1, L_1)$。右最小指标 $\varepsilon_1 = \varepsilon_2 = 1$，无左最小指标，无有限或无穷特征值。

!!! example "例 41.5"
    $L(\lambda) = \begin{pmatrix} \lambda & 1 & 0 \\ 0 & \lambda & 1 \\ 0 & 0 & 0 \end{pmatrix}$（$3 \times 3$ 束）。

    $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} -1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$。

    $\det(A - \lambda B) = \det\begin{pmatrix} \lambda & 1 & 0 \\ 0 & \lambda & 1 \\ 0 & 0 & 0 \end{pmatrix} = 0$。奇异束。

    右零空间：$\lambda x_1 + x_2 = 0$，$\lambda x_2 + x_3 = 0$。取 $x_1 = 1$，$x_2 = -\lambda$，$x_3 = \lambda^2$。最小基为 $(\lambda^2, -\lambda, 1)^T$，但需要检查是否有更低次的解。$x_3$ 任意，$x_2 = -\lambda x_3$，$x_1 = -x_2/\lambda = x_3$（仅当 $x_2/\lambda$ 为多项式时）。取 $x_3 = 1$，$x_2 = -\lambda$，$x_1 = \lambda^2$ 无误——但这是次数 2 的向量。

    左零空间：$y^T L = 0$：$\lambda y_1 = 0$，$y_1 + \lambda y_2 = 0$，$y_2 = 0$。故 $y_1 = y_2 = 0$，$y_3$ 任意。左最小指标 $\eta_1 = 0$。

    Kronecker 标准形包含一个 $L_2$ 块（右最小指标 2）和一个 $L_0^T$ 块（左最小指标 0，即 $1 \times 0$ 的空块——实际上 $L_0^T$ 表示一行额外的方程 $0 = 0$）。

---

## 41.5 广义 Schur 分解与 QZ 算法

<div class="context-flow" markdown>

**核心问题**：如何用数值稳定的算法计算正则束的广义特征值？

</div>

!!! theorem "定理 41.5 (广义 Schur 分解)"
    设 $A, B \in M_n(\mathbb{C})$。则存在酉矩阵 $Q, Z \in M_n(\mathbb{C})$ 使得
    $$QAZ = T_A, \quad QBZ = T_B,$$
    其中 $T_A, T_B$ 均为上三角矩阵。广义特征值为 $\lambda_i = (T_A)_{ii} / (T_B)_{ii}$（当 $(T_B)_{ii} \ne 0$）或 $\infty$（当 $(T_B)_{ii} = 0$）。

??? proof "证明"
    **存在性**可以通过同时 Schur 三角化来证明。取 $\mu$ 使得 $A - \mu B$ 非奇异，令 $C = (A - \mu B)^{-1}B$。对 $C$ 做 Schur 分解 $C = U T_C U^*$。然后取 $Q = U^*$，$Z = U$，将 $A$ 和 $B$ 同时变为上三角形式。

    更直接地，使用 QZ 迭代算法（见下文）。

!!! definition "定义 41.8 (QZ 算法)"
    **QZ 算法**是 Moler 和 Stewart 在 1973 年提出的计算正则束广义特征值的数值算法。其基本流程为：

    **阶段一（预处理）**：
    1. 利用 Householder 变换将 $B$ 化为上 Hessenberg 形式，同时对 $A$ 做相应变换。
    2. 进一步将 $A$ 化为上 Hessenberg 形式（同时保持 $B$ 的上三角性——通过 Givens 旋转）。

    **阶段二（迭代）**：
    QZ 迭代是 QR 迭代的推广。每步选取位移 $\mu$，计算
    $$QR = A - \mu B \quad \text{(QR 分解)},$$
    然后更新 $A \leftarrow R Q + \mu B Q$（等价操作），$B$ 也做相应更新以保持上三角。

    收敛后，$A$ 和 $B$ 同时为上三角，对角元之比给出广义特征值。

!!! theorem "定理 41.6 (QZ 算法的复杂度与稳定性)"
    QZ 算法的主要性质：

    (a) 预处理阶段：$O(n^3)$ 次运算。

    (b) 迭代阶段：每步 $O(n^2)$ 次运算，通常收敛很快（$O(n)$ 步迭代）。

    (c) 总复杂度：$O(n^3)$。

    (d) QZ 算法是**数值后向稳定的**：计算出的广义特征值是略微扰动后的束 $(A + \delta A) - \lambda(B + \delta B)$ 的精确广义特征值，其中 $\|\delta A\|, \|\delta B\| = O(\epsilon_{\text{mach}} \|A\|)$。

!!! example "例 41.6"
    计算束 $A - \lambda B$ 的广义特征值，$A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，$B = \begin{pmatrix} 5 & 6 \\ 7 & 8 \end{pmatrix}$。

    $\det(A - \lambda B) = (1-5\lambda)(4-8\lambda) - (2-6\lambda)(3-7\lambda)$
    $= 4 - 8\lambda - 20\lambda + 40\lambda^2 - (6 - 14\lambda - 18\lambda + 42\lambda^2)$
    $= 4 - 28\lambda + 40\lambda^2 - 6 + 32\lambda - 42\lambda^2$
    $= -2\lambda^2 + 4\lambda - 2 = -2(\lambda^2 - 2\lambda + 1) = -2(\lambda - 1)^2$.

    广义特征值为 $\lambda = 1$（重数 2）。

!!! note "注记 41.1 (QZ 算法的实现)"
    QZ 算法在 LAPACK 中由子程序 `xGGES`（Schur 形式）和 `xGGEV`（特征值和特征向量）实现。在 MATLAB 中，函数 `eig(A, B)` 调用 QZ 算法。该算法是工程计算中求解广义特征值问题的标准工具。

---

## 41.6 微分代数方程

<div class="context-flow" markdown>

**核心问题**：矩阵束理论如何应用于微分代数方程（DAE）的分析？

</div>

!!! definition "定义 41.9 (线性常系数 DAE)"
    **微分代数方程**（Differential-Algebraic Equation, DAE）是形如
    $$E\dot{x}(t) = Ax(t) + f(t) \tag{41.1}$$
    的方程，其中 $E, A \in M_n(\mathbb{R})$，$E$ 可能是**奇异的**。当 $E$ 非奇异时，(41.1) 化为标准 ODE $\dot{x} = E^{-1}Ax + E^{-1}f$。当 $E$ 奇异时，(41.1) 包含微分方程和代数约束的混合。

!!! definition "定义 41.10 (DAE 的指标)"
    线性常系数 DAE (41.1) 的**（Kronecker）指标**定义为矩阵束 $E - \lambda A$（注意顺序！常写为 $\lambda E - A$ 或 $sE - A$）的 Weierstrass 标准形中幂零块 $N$ 的幂零指标（即 $N$ 的最大 Jordan 块大小）。

    - 指标 0：$E$ 非奇异（纯 ODE）。
    - 指标 1：$N$ 为对角矩阵（即所有幂零 Jordan 块大小为 1）。
    - 指标 $\nu$：$N$ 的最大 Jordan 块大小为 $\nu$。

!!! theorem "定理 41.7 (DAE 的可解性)"
    考虑 DAE $E\dot{x} = Ax + f$，矩阵束 $sE - A$ 正则。设 Weierstrass 标准形为
    $$P(sE - A)Q = \begin{pmatrix} sI - J & 0 \\ 0 & sN - I \end{pmatrix}.$$
    令 $Qx = \begin{pmatrix} y \\ z \end{pmatrix}$，$P^{-1}f = \begin{pmatrix} g \\ h \end{pmatrix}$。则 DAE 等价于

    $$\dot{y} = Jy + g(t) \tag{ODE 部分}$$
    $$N\dot{z} = z + h(t) \tag{代数部分}$$

    代数部分的解为
    $$z(t) = -\sum_{k=0}^{\nu-1} N^k h^{(k)}(t),$$
    其中 $\nu$ 为 $N$ 的幂零指标，$h^{(k)}$ 为 $h$ 的 $k$ 阶导数。

    因此 DAE 的可解性要求 $f$（从而 $h$）至少 $\nu - 1$ 次可微。**指标越高，对数据的光滑性要求越高。**

??? proof "证明"
    ODE 部分 $\dot{y} = Jy + g$ 的解由标准 ODE 理论给出。

    代数部分 $N\dot{z} = z + h$，即 $z = N\dot{z} - h = N(N\ddot{z} - \dot{h}) - h = N^2\ddot{z} - N\dot{h} - h$。迭代 $\nu$ 次（$N^\nu = 0$），得
    $$z = -h - N\dot{h} - N^2\ddot{h} - \cdots - N^{\nu-1}h^{(\nu-1)}.$$

!!! example "例 41.7"
    电路方程（RC 电路）：$C\dot{v} = -Gv + Bu$，其中 $C$ 为电容矩阵（可能奇异——当电路包含纯电阻支路时），$G$ 为电导矩阵，$v$ 为节点电压。

    若 $C$ 奇异，这就是一个 DAE。例如
    $$\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} \dot{v}_1 \\ \dot{v}_2 \end{pmatrix} = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} v_1 \\ v_2 \end{pmatrix} + \begin{pmatrix} u(t) \\ 0 \end{pmatrix}.$$

    $E = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$A = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}$。

    第二个方程 $0 = v_1 - v_2$ 是代数约束（$v_1 = v_2$）。代入第一个方程 $\dot{v}_1 = -v_1 + v_1 + u = u$，即 $\dot{v}_1 = u$。

    矩阵束 $sE - A = \begin{pmatrix} s+1 & -1 \\ -1 & 1 \end{pmatrix}$，$\det = s + 1 - 1 = s$。有限特征值 $s = 0$，无穷特征值（因为 $\deg(\det) = 1 < n = 2$，有 1 个无穷特征值）。指标 1。

!!! note "注记 41.2 (指标与数值困难)"
    DAE 的指标直接影响数值方法的选择和难度：

    - 指标 1 的 DAE 可以用隐式方法（如 BDF）直接求解，与标准 ODE 类似。
    - 指标 2 的 DAE 需要特殊处理（如指标归约或投影方法）。
    - 指标 3 或更高的 DAE 在数值上极具挑战性，容易出现漂移（drift off constraint）现象。

    电路模拟软件 SPICE 中的方程通常具有指标 1 或 2。多体动力学方程（如机器人或车辆动力学）通常具有指标 3。

---

## 41.7 应用

<div class="context-flow" markdown>

**核心问题**：矩阵束理论在控制理论和工程中有哪些重要应用？

</div>

!!! definition "定义 41.11 (描述子系统)"
    **描述子系统**（descriptor system）是形如
    $$E\dot{x}(t) = Ax(t) + Bu(t), \quad y(t) = Cx(t) + Du(t)$$
    的线性系统，其中 $E$ 可能奇异。当 $E$ 非奇异时退化为标准状态空间描述。描述子系统的传递函数为
    $$G(s) = C(sE - A)^{-1}B + D.$$

!!! theorem "定理 41.8 (描述子系统的正则性)"
    描述子系统 $(E, A, B, C, D)$ 有唯一解当且仅当矩阵束 $sE - A$ 是**正则的**。此时传递函数 $G(s)$ 是有理的（rational proper + polynomial part）。

!!! theorem "定理 41.9 (多项式特征值问题的线性化)"
    $m$ 次多项式特征值问题
    $$P(\lambda)x = (\lambda^m A_m + \lambda^{m-1}A_{m-1} + \cdots + A_0)x = 0$$
    可以通过**伴随线性化**（companion linearization）转化为线性矩阵束问题：
    $$\begin{pmatrix} A_m & & \\ & I & \\ & & \ddots & \\ & & & I \end{pmatrix}\lambda - \begin{pmatrix} -A_{m-1} & -A_{m-2} & \cdots & -A_0 \\ I & 0 & \cdots & 0 \\ & \ddots & & \vdots \\ & & I & 0 \end{pmatrix}.$$

    这将 $n \times n$ 的 $m$ 次问题转化为 $mn \times mn$ 的线性束问题，然后可用 QZ 算法求解。

!!! example "例 41.8"
    二次特征值问题（在结构振动分析中常见）：
    $$(\lambda^2 M + \lambda C + K)x = 0,$$
    其中 $M$ 为质量矩阵，$C$ 为阻尼矩阵，$K$ 为刚度矩阵。

    标准线性化为
    $$\begin{pmatrix} M & 0 \\ 0 & I \end{pmatrix}\lambda\begin{pmatrix} v \\ \lambda v \end{pmatrix} = \begin{pmatrix} -C & -K \\ I & 0 \end{pmatrix}\begin{pmatrix} v \\ \lambda v \end{pmatrix},$$
    即 $\lambda B \tilde{x} = A \tilde{x}$ 的形式，$2n \times 2n$ 的广义特征值问题。

!!! note "注记 41.3 (Rosenbrock 系统矩阵)"
    在控制理论中，**Rosenbrock 系统矩阵**定义为
    $$S(s) = \begin{pmatrix} sE - A & -B \\ C & D \end{pmatrix}.$$
    $S(s)$ 是一个矩阵束（或矩阵多项式），其 Smith 标准形包含了系统的零点和极点信息。系统的**传输零点**（transmission zeros）是 $S(s)$ 的不变零点（使 $S(s)$ 的秩下降的 $s$ 值）。

!!! note "注记 41.4 (结构化矩阵束)"
    许多应用中的矩阵束具有特殊结构：

    - **对称束** $A - \lambda B$（$A, B$ 对称）：出现在振动分析中。可用对称 QZ 变体保持对称性。
    - **Hermite 束**（$A, B$ Hermite，$B$ 半正定）：特征值全为实数。
    - **回文束** $A - \lambda A^T$：出现在控制理论的最优控制问题中。
    - **Hamilton 束**：出现在 Riccati 方程的求解中。

!!! theorem "定理 41.10 (对称正定束的性质)"
    设 $A - \lambda B$ 为对称束，$B$ 正定。则：

    (a) 所有广义特征值 $\lambda_i$ 为实数。

    (b) 广义特征向量可选为 $B$-正交的：$v_i^T B v_j = \delta_{ij}$。

    (c) 存在非奇异矩阵 $S$ 使得 $S^T A S = \Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，$S^T B S = I$。

??? proof "证明"
    由 $B > 0$，令 $B = LL^T$（Cholesky 分解），将广义问题 $Ax = \lambda Bx$ 化为标准对称特征值问题 $L^{-1}AL^{-T}y = \lambda y$（$y = L^T x$）。$L^{-1}AL^{-T}$ 对称，故特征值全为实数。将标准特征向量 $\{y_i\}$ 变回 $\{x_i = L^{-T}y_i\}$，即得 $B$-正交的广义特征向量。

!!! example "例 41.9"
    广义特征值问题 $Ax = \lambda Bx$，
    $$A = \begin{pmatrix} 5 & 2 \\ 2 & 2 \end{pmatrix}, \quad B = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}.$$

    $\det(A - \lambda B) = (5-2\lambda)(2-\lambda) - (2-\lambda)^2 = (2-\lambda)(5-2\lambda - 2+\lambda) = (2-\lambda)(3-\lambda)$。

    广义特征值 $\lambda_1 = 2, \lambda_2 = 3$，均为实数（$B$ 正定）。

    $\lambda_1 = 2$：$(A-2B)x = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}x = 0$，$x_1 = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$。
    $\lambda_2 = 3$：$(A-3B)x = \begin{pmatrix} -1 & -1 \\ -1 & -1 \end{pmatrix}x = 0$，$x_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$。

    验证 $B$-正交：$x_1^T B x_2 = (0, 1)\begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ -1 \end{pmatrix} = (0,1)(1, 0)^T = 0$。正交性成立。

!!! note "注记 41.5 (矩阵束的扰动理论)"
    矩阵束的特征值对扰动的敏感性比标准特征值问题更为复杂。对正则束 $A - \lambda B$，若 $\lambda_0$ 为简单广义特征值，左右广义特征向量为 $y, x$（$Ax = \lambda_0 Bx$，$y^*A = \lambda_0 y^*B$），则
    $$\lambda_0 \text{ 的条件数} = \frac{\|y\| \|x\|}{|y^* B x|}.$$
    当 $y^* Bx \approx 0$ 时，广义特征值对扰动非常敏感。特别地，当 $B$ 接近奇异时，某些广义特征值可能变得极度病态。

    对奇异束，Kronecker 标准形本身在一般扰动下是不稳定的（最小指标可以改变），这使得奇异束的数值计算需要特别小心。Van Dooren 等人发展了用 SVD 等方法数值稳定地提取 Kronecker 结构的算法。

---

## 本章小结

| 概念 | 要点 |
|------|------|
| 矩阵束 $A - \lambda B$ | 方阵：正则（$\det \not\equiv 0$）或奇异 |
| Weierstrass 形 | 正则束 = Jordan 块（有限特征值）+ 幂零块（无穷特征值） |
| 最小指标 | 奇异束的核心不变量，描述多项式零空间的结构 |
| Kronecker 标准形 | $L_\varepsilon$（右块）+ $L_\eta^T$（左块）+ Jordan + 幂零：完全分类 |
| QZ 算法 | $O(n^3)$，后向稳定，LAPACK 标准实现 |
| DAE | $E\dot{x} = Ax$；指标 = 幂零块大小；指标越高越难求解 |
| 描述子系统 | 控制理论中正则束的应用 |

Kronecker 标准形是矩阵束理论的顶峰，它将正则束的 Weierstrass 理论和奇异束的最小指标理论统一在一个完整的框架中。从数值计算（QZ 算法）到微分代数方程（指标理论）再到控制理论（描述子系统），矩阵束的语言提供了分析线性系统结构的强大工具。
