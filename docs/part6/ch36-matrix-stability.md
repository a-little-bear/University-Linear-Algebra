# 第 36 章 矩阵稳定性与惯性

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵分析(Ch14) · 矩阵方程(Ch20)

**本章脉络**：稳定矩阵(Hurwitz) $\to$ Routh-Hurwitz 准则 $\to$ Schur-Cohn 准则 $\to$ Lyapunov 稳定性定理 $\to$ 惯性三元组 $\to$ Sylvester 惯性定理 $\to$ Ostrowski-Schneider $\to$ D-稳定性 $\to$ Kharitonov 定理 $\to$ Kreiss 矩阵定理

**延伸**：矩阵稳定性理论是控制系统设计（状态反馈镇定）和生态学（Lotka-Volterra 系统稳定性分析）的数学基础

</div>

线性常微分方程组 $\dot{x}(t) = Ax(t)$ 的解的渐近行为完全由矩阵 $A$ 的特征值决定。当所有特征值位于复平面的左半平面时，所有解都指数衰减至零——此时称系统**稳定**。这一简洁的判据将动力系统的定性理论归结为矩阵的**谱定位**问题。本章系统发展矩阵稳定性的代数判据、Lyapunov 方法以及惯性理论，揭示矩阵谱的正负分布在合同变换下的不变性。

从 Hurwitz 在 19 世纪末提出的特征多项式系数判据，到 Lyapunov 在 1892 年博士论文中建立的矩阵方程方法，再到 Sylvester、Ostrowski、Schneider 等人发展的惯性理论，矩阵稳定性分析已成为线性代数与动力系统交叉领域的核心课题。

---

## 36.1 稳定矩阵的定义

<div class="context-flow" markdown>

**核心问题**：什么条件下线性系统 $\dot{x} = Ax$ 的所有解趋于零？离散系统 $x_{k+1} = Ax_k$ 呢？

</div>

稳定性概念因连续时间与离散时间系统的不同而有两种基本形式。

!!! definition "定义 36.1 (Hurwitz 稳定)"
    设 $A \in M_n(\mathbb{C})$。若 $A$ 的所有特征值 $\lambda$ 满足
    $$\operatorname{Re}(\lambda) < 0,$$
    则称 $A$ 为 **Hurwitz 稳定矩阵**（或连续时间稳定矩阵）。全体 $n$ 阶 Hurwitz 稳定矩阵的集合记为 $\mathcal{H}_n$。

!!! definition "定义 36.2 (Schur 稳定)"
    设 $A \in M_n(\mathbb{C})$。若 $A$ 的所有特征值 $\lambda$ 满足
    $$|\lambda| < 1,$$
    则称 $A$ 为 **Schur 稳定矩阵**（或离散时间稳定矩阵）。全体 $n$ 阶 Schur 稳定矩阵的集合记为 $\mathcal{S}_n$。

!!! theorem "定理 36.1 (连续系统稳定性)"
    设 $A \in M_n(\mathbb{C})$，考虑初值问题
    $$\dot{x}(t) = Ax(t), \quad x(0) = x_0.$$
    则以下条件等价：

    (a) $A$ 是 Hurwitz 稳定的；

    (b) 对所有 $x_0 \in \mathbb{C}^n$，有 $\lim_{t\to+\infty} e^{tA} x_0 = 0$；

    (c) $\lim_{t\to+\infty} \|e^{tA}\| = 0$（对任意矩阵范数）；

    (d) 存在常数 $C > 0, \alpha > 0$ 使得 $\|e^{tA}\| \le C e^{-\alpha t}$ 对所有 $t \ge 0$ 成立。

??? proof "证明"
    **(a) $\Rightarrow$ (d)**：设 $A$ 的 Jordan 标准形为 $A = PJP^{-1}$，其中 $J = \operatorname{diag}(J_1, \ldots, J_k)$，$J_i = \lambda_i I_{n_i} + N_i$ 为 Jordan 块。则
    $$e^{tA} = P e^{tJ} P^{-1} = P \operatorname{diag}(e^{tJ_1}, \ldots, e^{tJ_k}) P^{-1}.$$
    每个 Jordan 块的矩阵指数为
    $$e^{tJ_i} = e^{t\lambda_i} \begin{pmatrix} 1 & t & \frac{t^2}{2!} & \cdots \\ 0 & 1 & t & \cdots \\ \vdots & & \ddots & \\ 0 & \cdots & 0 & 1 \end{pmatrix}.$$
    当 $\operatorname{Re}(\lambda_i) < 0$ 时，$|e^{t\lambda_i}| = e^{t\operatorname{Re}(\lambda_i)} \to 0$，且多项式增长被指数衰减压制。
    取 $\alpha = \min_i (-\operatorname{Re}(\lambda_i)) / 2 > 0$，则存在 $C > 0$ 使得 $\|e^{tA}\| \le C e^{-\alpha t}$。

    **(d) $\Rightarrow$ (c) $\Rightarrow$ (b)** 显然。

    **(b) $\Rightarrow$ (a)**：反证法。若存在特征值 $\lambda$ 满足 $\operatorname{Re}(\lambda) \ge 0$，取对应特征向量 $v$，则 $e^{tA} v = e^{t\lambda} v$，从而 $\|e^{tA} v\| = |e^{t\lambda}| \|v\| = e^{t\operatorname{Re}(\lambda)} \|v\| \ge \|v\| > 0$，矛盾。

!!! theorem "定理 36.2 (离散系统稳定性)"
    设 $A \in M_n(\mathbb{C})$，考虑迭代 $x_{k+1} = Ax_k$。则以下等价：

    (a) $A$ 是 Schur 稳定的；

    (b) 对所有 $x_0$，$\lim_{k\to\infty} A^k x_0 = 0$；

    (c) $\lim_{k\to\infty} \|A^k\| = 0$；

    (d) $\rho(A) < 1$（谱半径小于 1）。

??? proof "证明"
    (a) 与 (d) 等价是谱半径的定义。(d) $\Rightarrow$ (c)：由 Gelfand 公式 $\rho(A) = \lim_{k\to\infty} \|A^k\|^{1/k}$，若 $\rho(A) < 1$，则对 $\rho(A) < r < 1$，当 $k$ 充分大时 $\|A^k\| \le r^k \to 0$。其余方向类似定理 36.1 的证明。

!!! example "例 36.1"
    考虑矩阵 $A = \begin{pmatrix} -1 & 2 \\ 0 & -3 \end{pmatrix}$。
    特征值为 $\lambda_1 = -1, \lambda_2 = -3$，均有负实部，故 $A$ 是 Hurwitz 稳定的。

    矩阵 $B = \begin{pmatrix} 0.5 & 0.1 \\ -0.2 & 0.3 \end{pmatrix}$ 的特征值为 $\lambda = 0.4 \pm 0.1\sqrt{2}\,i$，模为 $\sqrt{0.16 + 0.02} = \sqrt{0.18} < 1$，故 $B$ 是 Schur 稳定的。

!!! note "注记 36.1 (Cayley 变换与两种稳定性的联系)"
    Cayley 变换 $B = (I + A)(I - A)^{-1}$ 将左半平面映射到单位圆盘内部。因此 $A$ 是 Hurwitz 稳定的当且仅当 $B = (I + A)(I - A)^{-1}$ 是 Schur 稳定的（假设 $1 \notin \sigma(A)$）。这一变换在控制理论中称为**双线性变换**，它连接了连续时间与离散时间系统的稳定性分析。

---

## 36.2 Routh-Hurwitz 准则

<div class="context-flow" markdown>

**核心问题**：如何仅从特征多项式的系数判定 Hurwitz 稳定性，而不必求出特征值？

</div>

对于实系数矩阵，特征多项式的系数携带了足够的信息来判定所有根是否位于左半平面。

!!! definition "定义 36.3 (Hurwitz 矩阵)"
    设实系数多项式 $p(\lambda) = \lambda^n + a_1\lambda^{n-1} + \cdots + a_{n-1}\lambda + a_n$。定义 **Hurwitz 矩阵**为
    $$H = \begin{pmatrix}
    a_1 & a_3 & a_5 & \cdots \\
    1   & a_2 & a_4 & \cdots \\
    0   & a_1 & a_3 & \cdots \\
    0   & 1   & a_2 & \cdots \\
    \vdots & & & \ddots
    \end{pmatrix} \in \mathbb{R}^{n \times n},$$
    其中 $a_j = 0$（当 $j > n$ 时）。记 $\Delta_k = \det(H[1:k, 1:k])$ 为 $H$ 的 $k$ 阶顺序主子式（称为 **Hurwitz 行列式**）。

!!! theorem "定理 36.3 (Routh-Hurwitz 准则)"
    设 $p(\lambda) = \lambda^n + a_1\lambda^{n-1} + \cdots + a_n$ 为实系数首一多项式。$p(\lambda)$ 的所有根具有负实部（即对应矩阵 Hurwitz 稳定）当且仅当所有 Hurwitz 行列式为正：
    $$\Delta_1 > 0, \quad \Delta_2 > 0, \quad \ldots, \quad \Delta_n > 0.$$

??? proof "证明（概要）"
    证明的核心思想基于 Hermite 关于实根判别的二次型方法以及 Cauchy 指标理论。

    **必要性方向**：若所有根 $\lambda_i$ 满足 $\operatorname{Re}(\lambda_i) < 0$，则可以证明 $a_k > 0$（所有系数为正，因为根的实部均为负的首一多项式的系数交替为正——但由于共轭复数根的存在，所有系数实际上为正）。进一步，通过构造与多项式根的实部相关的 Hermite 型矩阵，可以证明所有 $\Delta_k > 0$。

    **充分性方向**：采用连续性论证。考虑多项式族 $p_t(\lambda) = \lambda^n + t a_1 \lambda^{n-1} + \cdots + t^n a_n$，$t \in [0, 1]$。当 $t = 0$ 时，所有根为零。随着 $t$ 从 $0$ 增至 $1$，根连续变化。若所有 $\Delta_k > 0$，则可以证明在此过程中没有根穿过虚轴（否则某个 $\Delta_k$ 将经过零），因此最终 $p_1 = p$ 的所有根仍然在左半平面。

    完整的严格证明需要利用 Cauchy 的辐角原理和 Sturm 链理论，细节可参见 Gantmacher 的经典著作《矩阵论》第 XV 章。

!!! example "例 36.2"
    判定多项式 $p(\lambda) = \lambda^3 + 2\lambda^2 + 3\lambda + 4$ 是否 Hurwitz 稳定。

    系数：$a_1 = 2, a_2 = 3, a_3 = 4$。

    Hurwitz 矩阵为
    $$H = \begin{pmatrix} 2 & 4 & 0 \\ 1 & 3 & 0 \\ 0 & 2 & 4 \end{pmatrix}.$$

    $$\Delta_1 = 2 > 0, \quad \Delta_2 = \det\begin{pmatrix} 2 & 4 \\ 1 & 3 \end{pmatrix} = 6 - 4 = 2 > 0, \quad \Delta_3 = 4 \cdot \Delta_2 = 8 > 0.$$

    所有 Hurwitz 行列式为正，故 $p(\lambda)$ Hurwitz 稳定。

!!! example "例 36.3"
    判定 $p(\lambda) = \lambda^3 + \lambda^2 + \lambda + 2$ 是否 Hurwitz 稳定。

    系数：$a_1 = 1, a_2 = 1, a_3 = 2$。

    $$\Delta_1 = 1 > 0, \quad \Delta_2 = \det\begin{pmatrix} 1 & 2 \\ 1 & 1 \end{pmatrix} = 1 - 2 = -1 < 0.$$

    $\Delta_2 < 0$，故 $p(\lambda)$ 不是 Hurwitz 稳定的。

!!! note "注记 36.2 (低阶情形的简化条件)"
    对于低阶多项式，Routh-Hurwitz 条件可大幅简化：

    - **$n = 2$**：$\lambda^2 + a_1\lambda + a_2$ Hurwitz 稳定 $\iff$ $a_1 > 0, a_2 > 0$。
    - **$n = 3$**：$\lambda^3 + a_1\lambda^2 + a_2\lambda + a_3$ Hurwitz 稳定 $\iff$ $a_1 > 0, a_3 > 0, a_1 a_2 > a_3$。
    - **$n = 4$**：$\lambda^4 + a_1\lambda^3 + a_2\lambda^2 + a_3\lambda + a_4$ Hurwitz 稳定 $\iff$ $a_1 > 0, a_4 > 0, a_1 a_2 a_3 > a_3^2 + a_1^2 a_4$。

!!! definition "定义 36.4 (Routh 表)"
    **Routh 表**是 Routh-Hurwitz 准则的等价算法形式。对多项式 $p(\lambda) = a_0\lambda^n + a_1\lambda^{n-1} + \cdots + a_n$，构造如下表格：

    | 行   | 元素 |
    |------|------|
    | $s^n$ | $a_0, a_2, a_4, \ldots$ |
    | $s^{n-1}$ | $a_1, a_3, a_5, \ldots$ |
    | $s^{n-2}$ | $b_1, b_2, b_3, \ldots$ |
    | $s^{n-3}$ | $c_1, c_2, c_3, \ldots$ |
    | $\vdots$ | $\vdots$ |

    其中 $b_1 = \frac{a_1 a_2 - a_0 a_3}{a_1}$，$b_2 = \frac{a_1 a_4 - a_0 a_5}{a_1}$，以此类推。$p(\lambda)$ Hurwitz 稳定当且仅当 Routh 表第一列的元素全部同号。

---

## 36.3 Schur-Cohn 准则

<div class="context-flow" markdown>

**核心问题**：如何仅从多项式系数判定所有根是否位于单位圆盘内部？这是离散时间系统稳定性的代数判据。

</div>

Schur-Cohn 准则是 Routh-Hurwitz 准则在离散时间系统中的对应物。它判定多项式的所有根是否在单位圆盘 $|z| < 1$ 内。

### Schur-Cohn 矩阵

!!! definition "定义 36.5 (Schur-Cohn 矩阵)"
    设实系数多项式 $p(z) = a_0 z^n + a_1 z^{n-1} + \cdots + a_n$，$a_0 \neq 0$。定义**逆序多项式** $p^*(z) = z^n \overline{p(1/\bar{z})} = \bar{a}_n z^n + \bar{a}_{n-1} z^{n-1} + \cdots + \bar{a}_0$。

    对实系数多项式，$p^*(z) = a_n z^n + a_{n-1} z^{n-1} + \cdots + a_0$。

    **Schur-Cohn 矩阵**定义为 $n \times n$ 矩阵
    $$C = B^*B - A^*A$$
    其中

    $$
    A = \begin{pmatrix} a_0 & 0 & \cdots & 0 \\ a_1 & a_0 & \cdots & 0 \\ \vdots & & \ddots & \vdots \\ a_{n-1} & a_{n-2} & \cdots & a_0 \end{pmatrix}, \quad B = \begin{pmatrix} a_n & a_{n-1} & \cdots & a_1 \\ 0 & a_n & \cdots & a_2 \\ \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & a_n \end{pmatrix}.
    $$

!!! theorem "定理 36.3A (Schur-Cohn 准则)"
    设 $p(z) = a_0 z^n + a_1 z^{n-1} + \cdots + a_n$ 为实系数多项式，$a_0 > 0$。$p(z)$ 的所有根满足 $|z| < 1$（即 Schur 稳定）当且仅当 Schur-Cohn 矩阵 $C = B^*B - A^*A$ 正定。

    等价地，$p(z)$ Schur 稳定当且仅当 $C$ 的所有顺序主子式为正。

### Jury 稳定性判据

!!! definition "定义 36.6 (Jury 表)"
    对多项式 $p(z) = a_0 z^n + a_1 z^{n-1} + \cdots + a_n$（$a_0 > 0$），**Jury 表**按如下方式构造。

    **第一行**：$a_0, a_1, a_2, \ldots, a_n$

    **第二行**：$a_n, a_{n-1}, a_{n-2}, \ldots, a_0$（逆序）

    **第三行**：$b_0, b_1, b_2, \ldots, b_{n-1}$，其中

    $$b_k = \det\begin{pmatrix} a_0 & a_{n-k} \\ a_n & a_k \end{pmatrix} = a_0 a_k - a_n a_{n-k}, \quad k = 0, 1, \ldots, n-1$$

    **第四行**：$b_{n-1}, b_{n-2}, \ldots, b_0$（逆序）

    继续此过程：每两行产生新的两行，行长度减一，直到剩下一行长度为 1。设第 $(2i+1)$ 行的第一个元素（首元素）为 $r_i$。

!!! theorem "定理 36.3B (Jury 稳定性判据)"
    设 $p(z) = a_0 z^n + a_1 z^{n-1} + \cdots + a_n$，$a_0 > 0$。$p(z)$ 的所有根满足 $|z| < 1$ 当且仅当以下条件全部成立：

    1. $p(1) > 0$
    2. $(-1)^n p(-1) > 0$
    3. Jury 表中所有奇数行的首元素为正：$r_0 = a_0 > 0$，$r_1 = b_0 > 0$，$r_2 = c_0 > 0$，$\ldots$

    条件 3 等价于 $|a_n| < a_0$（第一级），以及后续各级的首元素均为正。

!!! example "例 36.2A"
    判定 $p(z) = z^3 - 0.5z^2 + 0.2z + 0.1$ 是否 Schur 稳定。

    系数：$a_0 = 1, a_1 = -0.5, a_2 = 0.2, a_3 = 0.1$。

    **条件 1**：$p(1) = 1 - 0.5 + 0.2 + 0.1 = 0.8 > 0$。✓

    **条件 2**：$(-1)^3 p(-1) = -(- 1 - 0.5 - 0.2 + 0.1) = -(-1.6) = 1.6 > 0$。✓

    **条件 3**（构造 Jury 表）：

    第一行：$1, -0.5, 0.2, 0.1$

    第二行：$0.1, 0.2, -0.5, 1$

    $|a_3/a_0| = 0.1 < 1$。✓

    $b_k = a_0 a_k - a_3 a_{3-k}$：

    $b_0 = 1 \cdot 1 - 0.1 \cdot 0.1 = 0.99$

    $b_1 = 1 \cdot (-0.5) - 0.1 \cdot 0.2 = -0.52$

    $b_2 = 1 \cdot 0.2 - 0.1 \cdot (-0.5) = 0.25$

    第三行：$0.99, -0.52, 0.25$

    第四行：$0.25, -0.52, 0.99$

    $c_0 = 0.99 \cdot 0.99 - 0.25 \cdot 0.25 = 0.9801 - 0.0625 = 0.9176 > 0$。✓

    所有条件满足，故 $p(z)$ Schur 稳定。

!!! note "注记 36.2A (Schur-Cohn 与 Routh-Hurwitz 的关系)"
    Schur-Cohn 准则与 Routh-Hurwitz 准则可通过双线性变换 $z = \frac{1+s}{1-s}$（或 $s = \frac{z-1}{z+1}$）相互转化。将 $p(z)$ 中的 $z$ 替换为 $\frac{1+s}{1-s}$，得到关于 $s$ 的多项式 $q(s)$，则 $p$ 的根在单位圆盘内当且仅当 $q$ 的根在左半平面。因此原则上可以用 Routh-Hurwitz 准则来判定 Schur 稳定性，但直接使用 Jury 表通常更简便。

---

## 36.4 Lyapunov 稳定性定理

<div class="context-flow" markdown>

**核心问题**：能否通过求解一个矩阵方程来判定稳定性，而不必计算特征值或特征多项式？

</div>

Lyapunov 在 1892 年的博士论文中提出了一个深刻的方法：通过构造"能量函数"（二次型 $V(x) = x^* P x$）来判定线性系统的稳定性。

!!! theorem "定理 36.4 (Lyapunov 稳定性定理)"
    设 $A \in M_n(\mathbb{C})$。则 $A$ 是 Hurwitz 稳定的当且仅当对每个（或某个）Hermite 正定矩阵 $Q > 0$，Lyapunov 方程
    $$A^* P + PA = -Q \tag{36.1}$$
    有唯一的 Hermite 正定解 $P > 0$。

??? proof "证明"
    **充分性**：设存在 $P > 0$ 和 $Q > 0$ 满足 $A^*P + PA = -Q$。定义 $V(x) = x^* P x$，则沿系统 $\dot{x} = Ax$ 的轨线，
    $$\dot{V} = \dot{x}^* P x + x^* P \dot{x} = x^* A^* P x + x^* P A x = x^*(A^*P + PA)x = -x^* Q x < 0$$
    对所有 $x \ne 0$ 成立。因此 $V$ 是严格递减的 Lyapunov 函数。

    由于 $V(x(t)) \ge 0$ 且严格递减，$V(x(t)) \to 0$，从而 $x(t) \to 0$。由定理 36.1，$A$ 是 Hurwitz 稳定的。

    **必要性**：设 $A$ 是 Hurwitz 稳定的。对给定 $Q > 0$，定义
    $$P = \int_0^{\infty} e^{tA^*} Q\, e^{tA}\, dt. \tag{36.2}$$
    由于 $A$ Hurwitz 稳定，$\|e^{tA}\| \le Ce^{-\alpha t}$，因此被积函数 $\|e^{tA^*} Q\, e^{tA}\| \le C^2 \|Q\| e^{-2\alpha t}$，积分收敛。

    **$P$ 是 Hermite 的**：$P^* = \int_0^{\infty} (e^{tA^*} Q\, e^{tA})^* dt = \int_0^{\infty} e^{tA^*} Q^*\, e^{tA}\, dt = P$（因为 $Q = Q^*$）。

    **$P > 0$**：对 $x \ne 0$，$x^* P x = \int_0^{\infty} \|Q^{1/2} e^{tA} x\|^2\, dt > 0$（被积函数连续、非负，且 $t=0$ 时为 $\|Q^{1/2}x\|^2 > 0$）。

    **$P$ 满足方程 (36.1)**：
    $$A^*P + PA = \int_0^{\infty} \left(A^* e^{tA^*} Q\, e^{tA} + e^{tA^*} Q\, e^{tA} A\right) dt = \int_0^{\infty} \frac{d}{dt}\left(e^{tA^*} Q\, e^{tA}\right) dt.$$
    $$= \left[e^{tA^*} Q\, e^{tA}\right]_0^{\infty} = 0 - Q = -Q.$$

    **唯一性**：设 $P_1, P_2$ 都满足 $A^*P_i + P_i A = -Q$，则 $A^*(P_1 - P_2) + (P_1 - P_2)A = 0$。设 $\Delta = P_1 - P_2$，则
    $$\frac{d}{dt}(e^{tA^*} \Delta\, e^{tA}) = e^{tA^*}(A^*\Delta + \Delta A)e^{tA} = 0,$$
    故 $e^{tA^*} \Delta\, e^{tA} = \Delta$ 对所有 $t$ 成立。令 $t \to \infty$，由稳定性得 $\Delta = 0$。

!!! example "例 36.4"
    判定 $A = \begin{pmatrix} -1 & 1 \\ 0 & -2 \end{pmatrix}$ 是否 Hurwitz 稳定。

    取 $Q = I$，解 Lyapunov 方程 $A^T P + PA = -I$。设 $P = \begin{pmatrix} p_{11} & p_{12} \\ p_{12} & p_{22} \end{pmatrix}$，则

    $$\begin{pmatrix} -1 & 0 \\ 1 & -2 \end{pmatrix}\begin{pmatrix} p_{11} & p_{12} \\ p_{12} & p_{22} \end{pmatrix} + \begin{pmatrix} p_{11} & p_{12} \\ p_{12} & p_{22} \end{pmatrix}\begin{pmatrix} -1 & 1 \\ 0 & -2 \end{pmatrix} = -I.$$

    展开得：
    $$-2p_{11} = -1 \implies p_{11} = 1/2,$$
    $$-3p_{12} + p_{11} = 0 \implies p_{12} = 1/6,$$
    $$-4p_{22} + 2p_{12} = -1 \implies p_{22} = 1/3.$$

    $P = \begin{pmatrix} 1/2 & 1/6 \\ 1/6 & 1/3 \end{pmatrix}$。验证正定性：$p_{11} = 1/2 > 0$，$\det(P) = 1/6 - 1/36 = 5/36 > 0$。故 $P > 0$，确认 $A$ 是 Hurwitz 稳定的。

!!! theorem "定理 36.5 (离散 Lyapunov 定理)"
    $A \in M_n(\mathbb{C})$ 是 Schur 稳定的当且仅当对每个 $Q > 0$，**离散 Lyapunov 方程**（又称 Stein 方程）
    $$A^* P A - P = -Q$$
    有唯一的正定解 $P > 0$。

??? proof "证明"
    **必要性**：$A$ Schur 稳定时，令 $P = \sum_{k=0}^{\infty} (A^*)^k Q A^k$，此级数因 $\rho(A) < 1$ 而绝对收敛。可以直接验证 $A^*PA - P = -Q$ 以及 $P > 0$。

    **充分性**：若存在 $P > 0$ 使得 $A^*PA - P = -Q < 0$，则对 $V(x) = x^*Px$ 有 $V(Ax) - V(x) = x^*(A^*PA - P)x = -x^*Qx < 0$，即 $V$ 沿离散轨道严格递减，推出 $A^k x \to 0$。

!!! note "注记 36.3 (Lyapunov 方程的数值求解)"
    Lyapunov 方程 $A^*P + PA = -Q$ 可通过 Bartels-Stewart 算法在 $O(n^3)$ 时间内求解：先将 $A$ Schur 分解为上三角形式 $A = UTU^*$，然后将方程转化为三角系统逐列求解。该算法在 MATLAB 中由函数 `lyap` 实现。

---

## 36.5 惯性与 Sylvester 惯性定理

<div class="context-flow" markdown>

**核心问题**：Hermite 矩阵的正、负、零特征值个数在合同变换下是否不变？

</div>

!!! definition "定义 36.7 (矩阵的惯性)"
    设 $A \in M_n(\mathbb{C})$。$A$ 的**惯性**定义为三元组
    $$\operatorname{In}(A) = (\pi(A), \nu(A), \delta(A)),$$
    其中

    - $\pi(A)$ = 具有正实部的特征值个数（计重数），
    - $\nu(A)$ = 具有负实部的特征值个数（计重数），
    - $\delta(A)$ = 具有零实部的特征值个数（计重数）。

    显然 $\pi(A) + \nu(A) + \delta(A) = n$。

!!! definition "定义 36.8 (合同)"
    设 $A, B \in M_n(\mathbb{C})$。若存在非奇异矩阵 $S \in M_n(\mathbb{C})$ 使得
    $$B = S^* A S,$$
    则称 $A$ 与 $B$ **（*-）合同**。若 $A, B$ 为实矩阵且 $S$ 为实矩阵，$B = S^T A S$，则称 $A$ 与 $B$ **（转置）合同**。

!!! theorem "定理 36.6 (Sylvester 惯性定理)"
    设 $H \in M_n(\mathbb{C})$ 为 Hermite 矩阵，$S \in M_n(\mathbb{C})$ 为非奇异矩阵。则
    $$\operatorname{In}(S^* H S) = \operatorname{In}(H).$$
    即 Hermite 矩阵的惯性在合同变换下不变。

??? proof "证明"
    **证法一（连续路径论证）**：设 $S = S(1)$，构造从 $I$ 到 $S$ 的连续路径 $S(t)$（$t \in [0,1]$），使得 $S(t)$ 对所有 $t$ 非奇异。（因为 $GL_n(\mathbb{C})$ 是连通的，这样的路径存在。）

    矩阵族 $H(t) = S(t)^* H S(t)$ 从 $H(0) = H$ 连续变化到 $H(1) = S^*HS$。在此过程中，特征值连续变化。

    关键观察：若 $H$ 的某个特征值 $\lambda_i(t)$ 在某时刻 $t_0$ 变为零，则 $\det(H(t_0)) = 0$，即 $\det(S(t_0)^* H S(t_0)) = |\det(S(t_0))|^2 \det(H) = 0$。但 $\det(S(t_0)) \ne 0$，故 $\det(H) = 0$，矛盾（除非 $H$ 本身奇异）。

    对一般情形（$H$ 可能奇异），可以将 $H$ 分解为正定部分、负定部分和零部分，分别在不变子空间上讨论。

    **证法二（子空间维数论证）**：不妨设 $H = \operatorname{diag}(I_p, -I_q, 0_r)$（$p = \pi, q = \nu, r = \delta$），这总可以通过合同变换达到。设 $B = S^* H S$。对 $B$ 同样做此操作，设 $B$ 合同于 $\operatorname{diag}(I_{p'}, -I_{q'}, 0_{r'})$。

    利用 Witt 定理或子空间维数论证：$H$ 定义的二次型 $x^*Hx$ 在 $p$ 维子空间 $V_+$ 上正定。合同变换 $x \mapsto Sx$ 将 $V_+$ 映为 $p$ 维子空间 $S^{-1}V_+$，且 $B$ 在 $S^{-1}V_+$ 上正定。类似地对负定子空间讨论，得 $p' = p, q' = q, r' = r$。

!!! example "例 36.5"
    设 $H = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$，$S = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$。

    $$S^T H S = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}\begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 2 & 1 \\ -1 & -1 \end{pmatrix} = \begin{pmatrix} 3 & 1 \\ 1 & 0 \end{pmatrix}.$$

    $\operatorname{In}(H) = (1, 1, 0)$。$S^THS$ 的特征值：$\lambda^2 - 3\lambda - 1 = 0$，$\lambda = \frac{3 \pm \sqrt{13}}{2}$，一个正、一个负。$\operatorname{In}(S^THS) = (1, 1, 0) = \operatorname{In}(H)$。验证了惯性定理。

!!! theorem "定理 36.7 (Hermite 矩阵的惯性与合同分类)"
    两个 Hermite 矩阵 $H_1, H_2 \in M_n(\mathbb{C})$ 合同当且仅当 $\operatorname{In}(H_1) = \operatorname{In}(H_2)$。

??? proof "证明"
    **必要性**由 Sylvester 惯性定理直接给出。

    **充分性**：设 $\operatorname{In}(H_1) = \operatorname{In}(H_2) = (p, q, r)$。则 $H_1, H_2$ 均合同于 $\operatorname{diag}(I_p, -I_q, 0_r)$，设 $S_1^* H_1 S_1 = S_2^* H_2 S_2 = D$。则 $H_1 = (S_1^{-*}) D S_1^{-1}$，$H_2 = (S_2^{-*}) D S_2^{-1}$。故 $H_2 = (S_2^{-*} S_1^*) H_1 (S_1 S_2^{-1})^{} = T^* H_1 T$，其中 $T = S_1 S_2^{-1}$。

!!! note "注记 36.4 (实对称情形)"
    对实对称矩阵，Sylvester 惯性定理说的是：实对称矩阵在实合同（$B = S^T A S$，$S$ 实非奇异）下的完全不变量恰好是惯性三元组 $(p, q, r)$。这等价于说，实二次型在实可逆线性替换下的标准形由正、负项的个数完全确定。

---

## 36.6 Ostrowski-Schneider 惯性定理

<div class="context-flow" markdown>

**核心问题**：对于一般的（非 Hermite）矩阵，是否存在类似 Lyapunov 定理的惯性判据？

</div>

Ostrowski 和 Schneider 在 1962 年将 Lyapunov 稳定性定理推广为惯性定理，揭示了 Lyapunov 方程的解与矩阵惯性之间的精确联系。

!!! theorem "定理 36.8 (Ostrowski-Schneider 惯性定理)"
    设 $A \in M_n(\mathbb{C})$，$H \in M_n(\mathbb{C})$ 为 Hermite 矩阵。若
    $$AH + HA^* = K$$
    其中 $K > 0$（正定），则 $\delta(A) = 0$（即 $A$ 没有纯虚特征值）且
    $$\operatorname{In}(H) = (\pi(A), \nu(A), 0).$$

??? proof "证明"
    **步骤一：$\delta(A) = 0$**。反证。设 $A$ 有纯虚特征值 $\lambda = i\omega$（$\omega \in \mathbb{R}$），对应特征向量 $v$。则
    $$v^* K v = v^*(AH + HA^*)v = v^*AHv + v^*HA^*v = (i\omega)(v^*Hv) + (v^*Hv)(-i\omega) = 0.$$
    但 $K > 0$ 意味着 $v^*Kv > 0$，矛盾。故 $\delta(A) = 0$。

    **步骤二：$\operatorname{In}(H) = (\pi(A), \nu(A), 0)$**。

    先证 $\delta(H) = 0$。若 $Hv = 0$ 对某 $v \ne 0$，则 $v^*Kv = v^*(AH + HA^*)v = 0$，与 $K > 0$ 矛盾。故 $H$ 非奇异。

    定义 $A(t) = tA + (1-t) \cdot \frac{1}{2}H^{-1}K$，$t \in [0, 1]$。
    则 $A(t)H + HA(t)^* = t(AH + HA^*) + (1-t)K = tK + (1-t)K = K > 0$。

    由步骤一的论证，$A(t)$ 对所有 $t \in [0,1]$ 都没有纯虚特征值。因此 $\pi(A(t))$ 和 $\nu(A(t))$ 在 $t$ 变化过程中保持不变（特征值不能穿过虚轴）。

    当 $t = 0$ 时，$A(0) = \frac{1}{2}H^{-1}K$。由于 $K > 0$，$H^{-1}K$ 的特征值与 $K^{1/2}H^{-1}K^{1/2}$ 的特征值相同（相似），而 $K^{1/2}H^{-1}K^{1/2}$ 是 Hermite 的，故特征值全为实数。$H^{-1}K$ 的正特征值个数等于 $H^{-1}$（从而 $H$）的正特征值个数 $\pi(H)$（这可由同时合同对角化论证），负特征值个数等于 $\nu(H)$。

    因此 $\pi(A) = \pi(A(0)) = \pi(H)$，$\nu(A) = \nu(A(0)) = \nu(H)$。

!!! example "例 36.6"
    设 $A = \begin{pmatrix} -1 & 2 \\ 0 & 3 \end{pmatrix}$。$A$ 的特征值为 $-1, 3$，故 $\operatorname{In}(A) = (1, 1, 0)$。

    解方程 $AH + HA^T = K$，取 $K = I$：

    $$\begin{pmatrix} -1 & 2 \\ 0 & 3 \end{pmatrix}\begin{pmatrix} h_{11} & h_{12} \\ h_{12} & h_{22} \end{pmatrix} + \begin{pmatrix} h_{11} & h_{12} \\ h_{12} & h_{22} \end{pmatrix}\begin{pmatrix} -1 & 0 \\ 2 & 3 \end{pmatrix} = I.$$

    展开各分量：

    $(1,1)$ 分量：$-h_{11} + 2h_{12} + (-h_{11} + 2h_{12}) = -2h_{11} + 4h_{12} = 1$。

    $(1,2)$ 分量：$-h_{12} + 2h_{22} + 3h_{12} = 2h_{12} + 2h_{22} = 0$，故 $h_{12} = -h_{22}$。

    $(2,2)$ 分量：$3h_{22} + 3h_{22} = 6h_{22} = 1$，故 $h_{22} = 1/6$，$h_{12} = -1/6$。

    $(1,1)$ 分量代入：$-2h_{11} + 4(-1/6) = 1$，$-2h_{11} = 1 + 2/3 = 5/3$，$h_{11} = -5/6$。

    $$H = \begin{pmatrix} -5/6 & -1/6 \\ -1/6 & 1/6 \end{pmatrix}.$$

    $\det(H) = (-5/6)(1/6) - (-1/6)^2 = -5/36 - 1/36 = -6/36 = -1/6 < 0$。

    因此 $H$ 有一个正特征值和一个负特征值：$\operatorname{In}(H) = (1, 1, 0) = \operatorname{In}(A)$。定理成立。

---

## 36.7 D-稳定性与全稳定性

<div class="context-flow" markdown>

**核心问题**：矩阵 $A$ 在与正对角矩阵相乘后是否仍然稳定？这一性质在生态学和经济学模型中为何重要？

</div>

在生态学的 Lotka-Volterra 模型和经济学的 Leontief 模型中，系统矩阵往往会被某个正对角矩阵缩放。D-稳定性保证了无论这种缩放如何改变，系统都保持稳定。

!!! definition "定义 36.9 (D-稳定性)"
    矩阵 $A \in M_n(\mathbb{R})$ 称为 **D-稳定的**，若对所有正对角矩阵 $D = \operatorname{diag}(d_1, \ldots, d_n)$（$d_i > 0$），$DA$ 都是 Hurwitz 稳定的。

!!! definition "定义 36.10 (Volterra 乘子与对角稳定)"
    正对角矩阵 $D > 0$ 称为 $A$ 的 **Volterra 乘子**（或 Lyapunov 缩放因子），若 $DA + A^T D < 0$（即 $DA + A^T D$ 负定）。

    若存在这样的 $D$，则称 $A$ 是**对角稳定的**（diagonally stable）。

!!! definition "定义 36.11 (全稳定性)"
    矩阵 $A \in M_n(\mathbb{R})$ 称为**全稳定的**（totally stable），若对所有正对角矩阵 $D_1, D_2 > 0$，$D_1 A D_2$ 都是 Hurwitz 稳定的。

!!! definition "定义 36.12 (加法 D-稳定)"
    矩阵 $A \in M_n(\mathbb{R})$ 称为**加法 D-稳定的**（additively D-stable），若对所有非负对角矩阵 $D \geq 0$，$A - D$ 都是 Hurwitz 稳定的。

!!! theorem "定理 36.9 (D-稳定性的充分条件)"
    若存在正对角矩阵 $D > 0$ 使得 $DA + A^T D < 0$（即 $A$ 对角稳定），则 $A$ 是 D-稳定的。

??? proof "证明"
    设 $\tilde{D}$ 为任意正对角矩阵。需要证明 $\tilde{D}A$ 是 Hurwitz 稳定的。

    取 $P = D\tilde{D}^{-1}$，这也是正对角矩阵。计算

    $$(\tilde{D}A)^T P + P(\tilde{D}A) = A^T\tilde{D}P + P\tilde{D}A = A^T D + DA < 0$$

    其中 $\tilde{D}P = \tilde{D} \cdot D\tilde{D}^{-1} = D$。由 Lyapunov 定理，$\tilde{D}A$ 是 Hurwitz 稳定的。

!!! theorem "定理 36.10 (2×2 矩阵 D-稳定性的完全刻画)"
    实矩阵 $A = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}$ 是 D-稳定的当且仅当以下条件同时成立：

    (a) $a_{11} < 0$，$a_{22} < 0$；

    (b) $\det(A) > 0$。

??? proof "证明"
    **必要性**：取 $D = I$，则 $A$ 本身必须 Hurwitz 稳定，这要求 $\operatorname{tr}(A) < 0$ 且 $\det(A) > 0$。

    对任意 $D = \operatorname{diag}(d_1, d_2)$（$d_i > 0$），$\operatorname{tr}(DA) = d_1 a_{11} + d_2 a_{22} < 0$ 必须对所有 $d_1, d_2 > 0$ 成立。若 $a_{11} \geq 0$，取 $d_1 \to \infty$，$d_2$ 有限，则 $\operatorname{tr}(DA) \to +\infty$，不可能始终为负。故 $a_{11} < 0$。类似地 $a_{22} < 0$。

    **充分性**：设 $a_{11} < 0, a_{22} < 0, \det(A) > 0$。对任意 $D = \operatorname{diag}(d_1, d_2)$（$d_i > 0$），$\operatorname{tr}(DA) = d_1 a_{11} + d_2 a_{22} < 0$，$\det(DA) = d_1 d_2 \det(A) > 0$。由特征值满足 $\lambda^2 - \operatorname{tr}(DA)\lambda + \det(DA) = 0$，两根之和为负、积为正，故两根均有负实部。

### 3x3 矩阵的 D-稳定性

!!! theorem "定理 36.11 (Cross 判据, 1978)"
    设 $A \in M_3(\mathbb{R})$ 是 Hurwitz 稳定的，特征多项式为 $p(\lambda) = \lambda^3 + c_1\lambda^2 + c_2\lambda + c_3$。则 $A$ 是 D-稳定的当且仅当以下条件全部成立：

    (a) $a_{ii} < 0$，$i = 1, 2, 3$；

    (b) 每个 $2 \times 2$ 主子矩阵的行列式为正：$a_{ii}a_{jj} - a_{ij}a_{ji} > 0$（$1 \leq i < j \leq 3$）；

    (c) $\det(A) < 0$（等价于 $c_3 > 0$）；

    (d) $c_1 c_2 > c_3$（Routh-Hurwitz 条件的加强）；

    (e) 对所有正对角矩阵 $D$，$\operatorname{tr}(\operatorname{adj}(DA)) > 0$，这归结为关于 $a_{ij}$ 的一组多项式不等式。

!!! note "注记 36.5 (高维情形的困难)"
    对 $n \ge 4$，D-稳定性的刻画变得极为复杂，至今没有简洁的充要条件。已知的结果包括：

    - 若 $A$ 对角稳定（存在正对角 $D$ 使 $DA + A^T D < 0$），则 $A$ 是 D-稳定的（定理 36.9）。但反之不然（$n \geq 3$ 时）。
    - 符号稳定矩阵（仅由符号模式决定稳定性的矩阵）一定是 D-稳定的。
    - D-稳定性蕴含加法 D-稳定性。

!!! theorem "定理 36.12 (对角稳定与 D-稳定的关系)"
    (a) 若 $A$ 对角稳定，则 $A$ 是 D-稳定的。

    (b) 对 $n = 2$，D-稳定等价于对角稳定。

    (c) 对 $n \geq 3$，存在 D-稳定但不对角稳定的矩阵。

??? proof "证明"
    **(a)** 已在定理 36.9 中证明。

    **(b)** 设 $A = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}$ 是 D-稳定的，则 $a_{11} < 0$，$a_{22} < 0$，$\det(A) > 0$。取 $D = \operatorname{diag}(d_1, d_2)$，其中 $d_1/d_2 = -a_{22}/a_{11} > 0$（由 $a_{11}, a_{22} < 0$）。则

    $$DA + A^T D = \begin{pmatrix} 2d_1 a_{11} & d_1 a_{12} + d_2 a_{21} \\ d_1 a_{12} + d_2 a_{21} & 2d_2 a_{22} \end{pmatrix}$$

    对角元素为负。行列式为 $4d_1 d_2 a_{11}a_{22} - (d_1 a_{12} + d_2 a_{21})^2$。当 $\det(A) > 0$ 时（即 $a_{11}a_{22} > a_{12}a_{21}$），可以选择适当的 $d_1/d_2$ 使行列式为正，从而 $DA + A^T D < 0$。

!!! example "例 36.7 (Lotka-Volterra 竞争模型)"
    两物种竞争的 Lotka-Volterra 模型在平衡点处的线性化矩阵为
    $$A = \begin{pmatrix} -r_1 x_1^* / K_1 & -r_1 x_1^* \alpha_{12} / K_1 \\ -r_2 x_2^* \alpha_{21} / K_2 & -r_2 x_2^* / K_2 \end{pmatrix}$$
    其中 $r_i, x_i^*, K_i > 0$ 为增长率、平衡种群密度和环境容纳量，$\alpha_{ij} > 0$ 为种间竞争系数。

    $a_{11} < 0, a_{22} < 0$ 自动成立。D-稳定性（此处等价于 $\det(A) > 0$）要求
    $$1 - \alpha_{12}\alpha_{21} > 0,$$
    即种间竞争系数之积小于 1。这是两物种共存平衡点稳定的著名条件。

    D-稳定性在这里的生物学意义是：即使两物种的内禀增长率 $r_i$ 和平衡密度 $x_i^*$ 发生变化（相当于对 $A$ 进行对角缩放），共存平衡点仍然稳定。

!!! note "注记 36.6 (经济学中的应用)"
    在一般均衡经济学中，Arrow 和 McManus (1958) 引入了 D-稳定性来分析市场调节过程的稳定性。Walrasian tatonnement 过程的线性化系统具有形如 $DA$ 的系统矩阵，其中 $D$ 是由调节速度决定的正对角矩阵。D-稳定性保证了无论调节速度如何选择，价格调节过程都会收敛到均衡。

---

## 36.8 Kharitonov 定理

<div class="context-flow" markdown>

**核心问题**：当多项式的系数存在不确定性（只知道区间范围）时，如何判定这族多项式是否都 Hurwitz 稳定？

</div>

Kharitonov 定理是鲁棒稳定性理论中最重要的结果之一。它将无穷族多项式的稳定性判定归结为仅检验四个"极端"多项式。

### 区间多项式

!!! definition "定义 36.13 (区间多项式族)"
    设 $\underline{a}_k \leq \bar{a}_k$（$k = 0, 1, \ldots, n$）为给定的实数区间端点。**区间多项式族**定义为

    $$
    \mathcal{P} = \left\{p(\lambda) = \sum_{k=0}^{n} a_k \lambda^k : \underline{a}_k \leq a_k \leq \bar{a}_k, \; k = 0, \ldots, n\right\}
    $$

    即系数在给定矩形区域内变化的所有多项式。

!!! theorem "定理 36.13 (Kharitonov 定理, 1978)"
    设区间多项式族 $\mathcal{P}$ 中的所有多项式都是 $n$ 次的（即 $\underline{a}_n > 0$）。则 $\mathcal{P}$ 中的所有多项式都 Hurwitz 稳定当且仅当以下四个 **Kharitonov 多项式**都 Hurwitz 稳定：

    $$
    K_1(\lambda) = \underline{a}_0 + \underline{a}_1 \lambda + \bar{a}_2 \lambda^2 + \bar{a}_3 \lambda^3 + \underline{a}_4 \lambda^4 + \underline{a}_5 \lambda^5 + \cdots
    $$

    $$
    K_2(\lambda) = \underline{a}_0 + \bar{a}_1 \lambda + \bar{a}_2 \lambda^2 + \underline{a}_3 \lambda^3 + \underline{a}_4 \lambda^4 + \bar{a}_5 \lambda^5 + \cdots
    $$

    $$
    K_3(\lambda) = \bar{a}_0 + \underline{a}_1 \lambda + \underline{a}_2 \lambda^2 + \bar{a}_3 \lambda^3 + \bar{a}_4 \lambda^4 + \underline{a}_5 \lambda^5 + \cdots
    $$

    $$
    K_4(\lambda) = \bar{a}_0 + \bar{a}_1 \lambda + \underline{a}_2 \lambda^2 + \underline{a}_3 \lambda^3 + \bar{a}_4 \lambda^4 + \bar{a}_5 \lambda^5 + \cdots
    $$

    系数选取规则为：将系数按 $a_0, a_1, a_2, a_3, a_4, a_5, \ldots$ 排列，按"下下上上"或"上上下下"等四种模式循环选取上界 $\bar{a}_k$ 或下界 $\underline{a}_k$。

??? proof "证明（概要）"
    证明的关键是利用 Hermite-Biehler 定理。

    设 $p(\lambda) = U(\omega^2) + \lambda V(\omega^2)$（其中 $\lambda = j\omega$），分离实部和虚部后得到两个关于 $\omega^2$ 的实多项式 $U$ 和 $V$。Hermite-Biehler 定理断言 $p$ Hurwitz 稳定当且仅当 $U$ 和 $V$ 的根都是实的、为负的，并且互相交错。

    对区间多项式族，$U$ 和 $V$ 的系数分别独立地在区间内变化。Kharitonov 的关键洞察是：$U$ 的稳定性（根为负实数）仅取决于 $U$ 的系数的某两种极端组合，$V$ 类似。将 $U$ 和 $V$ 的极端组合配对，恰好得到四个 Kharitonov 多项式。

    这个定理的深刻之处在于：它将检验无穷多个多项式的稳定性归结为检验恰好四个多项式，与多项式的次数 $n$ 和区间的大小无关。$\blacksquare$

!!! example "例 36.8"
    考虑区间多项式族

    $$p(\lambda) = a_3\lambda^3 + a_2\lambda^2 + a_1\lambda + a_0$$

    其中 $a_3 \in [1, 2]$，$a_2 \in [3, 5]$，$a_1 \in [2, 4]$，$a_0 \in [1, 3]$。

    四个 Kharitonov 多项式为：

    $K_1(\lambda) = 1 + 2\lambda + 5\lambda^2 + 2\lambda^3$

    $K_2(\lambda) = 1 + 4\lambda + 5\lambda^2 + 1\lambda^3$

    $K_3(\lambda) = 3 + 2\lambda + 3\lambda^2 + 2\lambda^3$

    $K_4(\lambda) = 3 + 4\lambda + 3\lambda^2 + 1\lambda^3$

    对每个 $K_i$ 应用 Routh-Hurwitz 准则。以 $K_1(\lambda) = 2\lambda^3 + 5\lambda^2 + 2\lambda + 1$ 为例：

    $a_1 = 5 > 0$，$a_3 = 1 > 0$，$a_1 a_2 = 5 \cdot 2 = 10 > 1 = a_3$。✓

    逐一验证四个多项式均 Hurwitz 稳定，即可断定整个区间族都 Hurwitz 稳定。

!!! note "注记 36.7 (Kharitonov 定理的推广)"
    Kharitonov 定理有许多重要的推广：

    - **多面体不确定性**：当系数的不确定性集合是多面体（而非矩形）时，需要检验的极端多项式数量可能增加，但边缘定理（Edge Theorem, Bartlett-Hollot-Lin 1988）将检验归结为多面体的边上。
    - **离散时间版本**：类似的结果也存在于 Schur 稳定性判定中，但直接的四多项式结论不再成立，需要检验多面体的暴露边。
    - **矩阵情形**：对区间矩阵族的 Hurwitz 稳定性，问题变为 NP-hard，Kharitonov 型的有限检验不再可能。

---

## 36.9 Kreiss 矩阵定理

<div class="context-flow" markdown>

**核心问题**：矩阵的幂 $A^k$ 的增长率如何由谱性质来界定？

</div>

Kreiss 矩阵定理给出了矩阵幂的一致有界性与预解式增长之间的等价关系，是数值分析中判定离散迭代稳定性的核心工具。

### Kreiss 常数

!!! definition "定义 36.14 (Kreiss 常数)"
    设 $A \in M_n(\mathbb{C})$。$A$ 的 **Kreiss 常数**定义为

    $$
    \mathcal{K}(A) = \sup_{|z| > 1} (|z| - 1) \|(zI - A)^{-1}\|
    $$

    其中 $\|\cdot\|$ 为谱范数（2-范数）。

!!! theorem "定理 36.14 (Kreiss 矩阵定理)"
    设 $A \in M_n(\mathbb{C})$。以下条件等价：

    (a) $\sup_{k \geq 0} \|A^k\| < \infty$（幂一致有界）；

    (b) $\mathcal{K}(A) < \infty$。

    更精确地，存在仅依赖于 $n$ 的常数使得

    $$
    \mathcal{K}(A) \leq \sup_{k \geq 0} \|A^k\| \leq e n \cdot \mathcal{K}(A)
    $$

    左侧不等式对所有 $n$ 成立。右侧不等式中的 $en$ 是最优阶：存在 $n \times n$ 矩阵 $A$ 使得 $\sup_k \|A^k\| \geq c \cdot n \cdot \mathcal{K}(A)$。

??? proof "证明（概要）"
    **左侧不等式** $\mathcal{K}(A) \leq \sup_k \|A^k\|$：对 $|z| > 1$，

    $$(zI - A)^{-1} = \frac{1}{z}\sum_{k=0}^{\infty} \left(\frac{A}{z}\right)^k = \sum_{k=0}^{\infty} \frac{A^k}{z^{k+1}}$$

    故 $\|(zI - A)^{-1}\| \leq \sum_{k=0}^{\infty} \frac{\|A^k\|}{|z|^{k+1}} \leq \frac{M}{|z| - 1}$（其中 $M = \sup_k \|A^k\|$，利用等比级数求和）。因此 $(|z| - 1)\|(zI - A)^{-1}\| \leq M$。

    **右侧不等式**的证明更为精细，需要利用矩阵幂的 Cauchy 积分表示

    $$A^k = \frac{1}{2\pi i} \oint_{|z| = r} z^k (zI - A)^{-1} dz$$

    选取适当的围道半径 $r > 1$，并对不同的 $r$ 优化估计。最优界 $en$ 的证明由 Spijker (1991) 和 LeVeque-Trefethen (1997) 给出。$\blacksquare$

!!! example "例 36.9"
    设 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$（$2 \times 2$ 幂零矩阵）。

    $A^0 = I$，$A^1 = A$，$A^k = 0$（$k \geq 2$）。$\sup_k \|A^k\| = 1$。

    $(zI - A)^{-1} = \begin{pmatrix} 1/z & 1/z^2 \\ 0 & 1/z \end{pmatrix}$，$\|(zI - A)^{-1}\| \leq \frac{1}{|z|} + \frac{1}{|z|^2}$。

    $\mathcal{K}(A) = \sup_{|z|>1} (|z|-1)\|(zI-A)^{-1}\| \leq \sup_{|z|>1} (|z|-1)\left(\frac{1}{|z|} + \frac{1}{|z|^2}\right) \leq 1$。

    Kreiss 定理的界：$1 \leq 1 \leq 2e \cdot 1$。✓

!!! note "注记 36.8 (Kreiss 定理的应用)"
    Kreiss 矩阵定理在数值分析中有重要应用：

    - **有限差分格式的稳定性**：离散化偏微分方程后，时间推进算子 $A$ 的幂 $A^k$ 的增长决定了数值格式的稳定性。Kreiss 定理将此问题化为预解式的估计。
    - **谱半径条件的不足**：$\rho(A) \leq 1$ 不足以保证 $\|A^k\|$ 有界（Jordan 块的情形）。Kreiss 条件 $\mathcal{K}(A) < \infty$ 给出了正确的充要条件。
    - **瞬态增长**：即使 $\rho(A) < 1$（渐近稳定），$\|A^k\|$ 可能先增大再衰减。Kreiss 常数定量地刻画了这种瞬态增长的上界。

---

## 36.10 本章小结

本章的核心脉络可以概括为：

| 概念 | 要点 |
|------|------|
| Hurwitz 稳定 | 所有特征值 $\operatorname{Re}(\lambda) < 0$；连续系统 $\dot{x}=Ax$ 渐近稳定 |
| Schur 稳定 | 所有特征值 $|\lambda| < 1$；离散系统 $x_{k+1}=Ax_k$ 渐近稳定 |
| Routh-Hurwitz | 从特征多项式系数判定 Hurwitz 稳定性；无需求根 |
| Schur-Cohn / Jury | 从特征多项式系数判定 Schur 稳定性；离散时间的代数判据 |
| Lyapunov 定理 | 稳定 $\iff$ Lyapunov 方程有正定解；构造性判据 |
| 惯性 | 三元组 $(\pi, \nu, \delta)$；合同不变量 |
| Sylvester 定理 | Hermite 矩阵惯性在合同变换下不变 |
| Ostrowski-Schneider | $AH+HA^*>0 \Rightarrow \operatorname{In}(H) = (\pi(A), \nu(A), 0)$ |
| D-稳定性 | $DA$ 对所有 $D>0$ 稳定；生态/经济学中的鲁棒稳定性 |
| Kharitonov 定理 | 区间多项式族的 Hurwitz 稳定性由四个极端多项式决定 |
| Kreiss 矩阵定理 | 矩阵幂一致有界 $\iff$ Kreiss 常数有限；$en$ 阶最优估计 |
