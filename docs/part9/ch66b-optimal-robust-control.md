# 第 66B 章 最优控制与鲁棒控制

<div class="context-flow" markdown>

**前置**：能控性与能观性(Ch66A) · Gramian 与 Lyapunov 方程(Ch66A) · 矩阵方程(Ch20) · 矩阵稳定性(Ch36) · Schur 补(Ch3)

**本章脉络**：Lyapunov 稳定性 → 能稳性与能检性 → 极点配置 → Ackermann 公式 → 观测器设计 → 分离原理 → LQR(连续与离散) → Riccati 方程 → LQR 鲁棒性 → Kalman 滤波 / LQG → H∞ 控制 → 线性矩阵不等式(LMI)

**延伸**：H∞ 控制推广到非线性 $H_\infty$ 控制（Hamilton-Jacobi 不等式）；LMI 方法是现代鲁棒控制的统一计算框架；自适应控制和强化学习中的线性代数结构

</div>

本章从 Lyapunov 稳定性理论出发，系统发展控制系统的设计方法。Lyapunov 矩阵方程将稳定性判定归结为正定矩阵的存在性问题。极点配置定理保证能控系统可以通过线性反馈任意配置闭环极点，Luenberger 观测器利用能观性从输出重构状态，分离原理允许控制器和观测器独立设计。线性二次调节器（LQR）将最优控制归结为代数 Riccati 方程，具有理想的鲁棒性裕度。Kalman 滤波是 LQR 的对偶问题，LQG 将两者结合。H∞ 控制处理模型不确定性下的最坏情况优化，线性矩阵不等式（LMI）提供统一的计算框架。

---

## 66B.1 Lyapunov 稳定性

<div class="context-flow" markdown>

**核心问题**：如何用矩阵方程判定线性系统 $\dot{x} = Ax$ 的渐近稳定性？

</div>

### 连续时间 Lyapunov 稳定性

!!! definition "定义 66B.1 (渐近稳定性)"
    系统 $\dot{x} = Ax$ 称为**渐近稳定的**（asymptotically stable），若对任意初始状态 $x_0$，$\lim_{t \to \infty} x(t) = \lim_{t \to \infty} e^{At}x_0 = 0$。

    这等价于 $A$ 的所有特征值具有严格负实部，即 $\operatorname{Re}(\lambda_i(A)) < 0$ 对所有 $i$。此时称 $A$ 为 **Hurwitz 矩阵**。

!!! theorem "定理 66B.1 (连续时间 Lyapunov 稳定性定理)"
    以下条件等价：

    1. $A$ 是 Hurwitz 矩阵（$\operatorname{Re}(\lambda_i) < 0$，$\forall i$）。
    2. 对任意正定矩阵 $Q \succ 0$，连续 Lyapunov 方程
    $$
    A^T P + PA = -Q
    $$
    有唯一解 $P$，且 $P \succ 0$。
    3. 存在 $P \succ 0$ 使得 $A^TP + PA \prec 0$。

??? proof "证明"
    **$(1) \Rightarrow (2)$**：若 $A$ 是 Hurwitz 矩阵，则 $e^{At} \to 0$（$t \to \infty$）。定义
    $$
    P = \int_0^{\infty} e^{A^T t} Q e^{At} \, dt
    $$

    此积分收敛（因为 $\|e^{At}\| \leq Me^{-\alpha t}$ 对某 $\alpha > 0$）。

    **$P \succ 0$ 的证明**：对任意非零 $x$，
    $$
    x^T P x = \int_0^{\infty} \|Q^{1/2} e^{At} x\|^2 \, dt > 0
    $$
    因为 $Q^{1/2} e^{At} x$ 是连续函数且在 $t = 0$ 处非零，所以积分严格正。

    **$P$ 满足 Lyapunov 方程的证明**：
    $$
    A^T P + PA = \int_0^{\infty} \left(A^T e^{A^T t} Q e^{At} + e^{A^T t} Q e^{At} A\right) dt = \int_0^{\infty} \frac{d}{dt}\left(e^{A^T t} Q e^{At}\right) dt
    $$
    $$
    = \left[e^{A^T t} Q e^{At}\right]_0^{\infty} = 0 - Q = -Q
    $$

    **唯一性**：若 $P_1, P_2$ 都满足方程，则 $A^T(P_1 - P_2) + (P_1 - P_2)A = 0$。设 $\Delta = P_1 - P_2$，则 $e^{A^Tt}\Delta e^{At}$ 对 $t$ 求导为零，所以 $e^{A^Tt}\Delta e^{At} = \Delta$。令 $t \to \infty$，左端趋于零，故 $\Delta = 0$。

    **$(2) \Rightarrow (3)$**：取 $Q = I$ 即可。

    **$(3) \Rightarrow (1)$**：设 $A^TP + PA = -Q_0 \prec 0$。定义 Lyapunov 函数 $V(x) = x^T P x$（$V > 0$ 对 $x \neq 0$）。则
    $$
    \dot{V} = \dot{x}^T P x + x^T P \dot{x} = x^T(A^TP + PA)x = -x^T Q_0 x < 0, \quad \forall x \neq 0
    $$

    因此 $V(x(t))$ 严格递减且有下界 $0$，所以 $x(t) \to 0$。

    更精确地，若 $A$ 有特征值 $\lambda$ 满足 $\operatorname{Re}(\lambda) \geq 0$，取对应特征向量 $v$（$Av = \lambda v$），则
    $$
    v^*(A^TP + PA)v = (\bar{\lambda} + \lambda) v^* P v = 2\operatorname{Re}(\lambda) \cdot v^*Pv \geq 0
    $$
    这与 $A^TP + PA \prec 0$ 矛盾。$\blacksquare$

### 离散时间 Lyapunov 稳定性

!!! definition "定义 66B.2 (离散 Schur 稳定性)"
    矩阵 $A$ 称为 **Schur 稳定的**，若其所有特征值的模严格小于 1：$|\lambda_i(A)| < 1$，$\forall i$。此时离散系统 $x_{k+1} = Ax_k$ 渐近稳定。

!!! theorem "定理 66B.2 (离散时间 Lyapunov 稳定性定理)"
    以下条件等价：

    1. $A$ 是 Schur 稳定的（$|\lambda_i| < 1$，$\forall i$）。
    2. 对任意 $Q \succ 0$，离散 Lyapunov 方程
    $$
    A^T P A - P = -Q
    $$
    有唯一正定解 $P \succ 0$。
    3. 存在 $P \succ 0$ 使得 $A^TPA - P \prec 0$。

??? proof "证明"
    **$(1) \Rightarrow (2)$**：若 $A$ Schur 稳定，定义
    $$
    P = \sum_{k=0}^{\infty} (A^T)^k Q A^k
    $$

    此级数收敛（因为 $\|A^k\| \to 0$）。

    **$P \succ 0$**：$x^T P x = \sum_{k=0}^{\infty} \|Q^{1/2} A^k x\|^2 > 0$（$x \neq 0$ 时 $k=0$ 项非零）。

    **满足方程**：
    $$
    A^TPA - P = \sum_{k=0}^{\infty} (A^T)^{k+1} Q A^{k+1} - \sum_{k=0}^{\infty} (A^T)^k Q A^k = -Q + \lim_{K \to \infty} (A^T)^{K+1} Q A^{K+1} = -Q
    $$

    **$(3) \Rightarrow (1)$**：定义 $V(x) = x^T P x$。则
    $$
    V(x_{k+1}) - V(x_k) = x_k^T(A^TPA - P)x_k < 0, \quad \forall x_k \neq 0
    $$

    因此 $V(x_k)$ 严格递减，$x_k \to 0$。特征值模必须小于 1，否则取对应特征向量可导出矛盾（与连续情形类似）。$\blacksquare$

---

## 66B.2 能稳性与能检性

<div class="context-flow" markdown>

**核心问题**：当系统不完全能控或能观时，LQR/LQG 设计还能适用吗？

</div>

!!! definition "定义 66B.3 (能稳性)"
    系统 $(A, B)$ 称为**能稳的**（stabilizable），若存在状态反馈 $K$ 使得 $A - BK$ 是 Hurwitz 矩阵。

    等价条件（PBH 形式）：对 $A$ 的每个**不稳定**特征值 $\lambda$（$\operatorname{Re}(\lambda) \geq 0$），
    $$
    \operatorname{rank}\begin{pmatrix} A - \lambda I & B \end{pmatrix} = n
    $$

!!! definition "定义 66B.4 (能检性)"
    系统 $(A, C)$ 称为**能检的**（detectable），若存在观测器增益 $L$ 使得 $A - LC$ 是 Hurwitz 矩阵。

    等价条件：对 $A$ 的每个不稳定特征值 $\lambda$，
    $$
    \operatorname{rank}\begin{pmatrix} A - \lambda I \\ C \end{pmatrix} = n
    $$

!!! note "注"
    能稳性和能检性是比能控性和能观性**更弱**的条件：

    - 能控 $\Rightarrow$ 能稳（反之不然）
    - 能观 $\Rightarrow$ 能检（反之不然）

    直觉上，能稳性只要求不稳定模态可被控制（稳定模态自然衰减），能检性只要求不稳定模态可被观测。

    对偶关系：$(A, B)$ 能稳当且仅当 $(A^T, B^T)$ 能检。

!!! example "例 66B.1"
    $A = \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$。

    能控性矩阵 $\mathcal{C} = \begin{pmatrix} 0 & 0 \\ 1 & 1 \end{pmatrix}$，$\operatorname{rank} = 1 < 2$，不能控。

    但不稳定特征值 $\lambda = 1$：$\operatorname{rank}(A - I, B) = \operatorname{rank}\begin{pmatrix} -2 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix} = 2$。

    因此系统能稳——不稳定模态可控，稳定模态 $\lambda = -1$ 自然衰减，不需要控制。

---

## 66B.3 极点配置

<div class="context-flow" markdown>

**核心问题**：通过线性状态反馈，能否任意指定闭环系统的特征值？

</div>

### 状态反馈

!!! definition "定义 66B.5 (线性状态反馈)"
    **线性状态反馈**（linear state feedback）是控制律
    $$
    u(t) = -Kx(t) + r(t)
    $$
    其中 $K \in \mathbb{R}^{m \times n}$ 是**反馈增益矩阵**，$r(t)$ 是参考输入。

    闭环系统为 $\dot{x} = (A - BK)x + Br$，闭环系统矩阵为 $A_{cl} = A - BK$。

### 极点配置定理

!!! theorem "定理 66B.3 (极点配置定理)"
    若 $(A, B)$ 能控，则对任意关于实轴对称的 $n$ 个复数集合 $\{\lambda_1, \ldots, \lambda_n\}$，存在实矩阵 $K \in \mathbb{R}^{m \times n}$，使得 $A - BK$ 的特征值恰好为 $\{\lambda_1, \ldots, \lambda_n\}$。

??? proof "证明"
    **单输入情形**（$m = 1$，$B = b$）：

    将系统化为能控标准型 $(\bar{A}, \bar{b})$（Ch66A 定义 66A.6）。在标准型中，取 $\bar{K} = (\bar{k}_0, \bar{k}_1, \ldots, \bar{k}_{n-1})$，则
    $$
    \bar{A} - \bar{b}\bar{K} = \begin{pmatrix} 0 & 1 & 0 & \cdots & 0 \\ 0 & 0 & 1 & \cdots & 0 \\ \vdots & & & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & 1 \\ -(a_0+\bar{k}_0) & -(a_1+\bar{k}_1) & \cdots & \cdots & -(a_{n-1}+\bar{k}_{n-1}) \end{pmatrix}
    $$

    其特征多项式为 $\lambda^n + (a_{n-1}+\bar{k}_{n-1})\lambda^{n-1} + \cdots + (a_0 + \bar{k}_0)$。

    要使特征值为 $\lambda_1, \ldots, \lambda_n$，令特征多项式等于 $\prod_{i=1}^n (\lambda - \lambda_i) = \lambda^n + \alpha_{n-1}\lambda^{n-1} + \cdots + \alpha_0$。

    解出 $\bar{k}_i = \alpha_i - a_i$，$i = 0, \ldots, n-1$。然后 $K = \bar{K} T^{-1}$（$T$ 是标准型变换矩阵）。

    **多输入情形**：需要利用 Heymann 引理，找到向量 $f$ 使得 $(A + Bf, B\hat{e}_j)$ 对某 $j$ 能控（$\hat{e}_j$ 是标准基向量），从而归约为单输入问题。$\blacksquare$

### Ackermann 公式

!!! theorem "定理 66B.4 (Ackermann 公式)"
    对于单输入能控系统 $(A, b)$，使闭环特征多项式为 $\alpha(\lambda) = \prod_{i=1}^n (\lambda - \lambda_i)$ 的反馈增益为
    $$
    K = e_n^T \mathcal{C}^{-1} \alpha(A)
    $$
    其中 $e_n = (0, \ldots, 0, 1)^T$，$\mathcal{C} = (b, Ab, \ldots, A^{n-1}b)$ 是能控性矩阵，$\alpha(A) = A^n + \alpha_{n-1}A^{n-1} + \cdots + \alpha_0 I$。

!!! example "例 66B.2"
    设 $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$。期望闭环极点为 $\lambda_1 = -5, \lambda_2 = -5$。

    $\mathcal{C} = \begin{pmatrix} 0 & 1 \\ 1 & -3 \end{pmatrix}$，$\det(\mathcal{C}) = -1 \neq 0$，能控。

    期望特征多项式：$\alpha(\lambda) = (\lambda + 5)^2 = \lambda^2 + 10\lambda + 25$。

    原特征多项式：$\chi_A(\lambda) = \lambda^2 + 3\lambda + 2$。

    取 $K = (k_1, k_2)$，使 $(A - BK)$ 的特征多项式为 $\lambda^2 + (3+k_2)\lambda + (2+k_1) = \lambda^2 + 10\lambda + 25$。

    解出 $k_2 = 7$，$k_1 = 23$。因此 $K = (23, 7)$。

---

## 66B.4 状态观测器

<div class="context-flow" markdown>

**核心问题**：当状态不能直接测量时，如何从输出重构状态？

</div>

### Luenberger 观测器

!!! definition "定义 66B.6 (Luenberger 观测器)"
    **Luenberger 观测器**是估计系统状态的动态系统：
    $$
    \dot{\hat{x}}(t) = A\hat{x}(t) + Bu(t) + L(y(t) - C\hat{x}(t))
    $$
    其中 $\hat{x}(t)$ 是状态估计，$L \in \mathbb{R}^{n \times p}$ 是**观测器增益**。

    定义估计误差 $e(t) = x(t) - \hat{x}(t)$，则
    $$
    \dot{e}(t) = (A - LC)e(t)
    $$

    误差动力学由矩阵 $A - LC$ 完全决定。

!!! theorem "定理 66B.5 (观测器极点配置)"
    若 $(A, C)$ 能观，则对任意期望的（关于实轴对称的）极点集合，存在 $L$ 使得 $A - LC$ 的特征值为给定值。

??? proof "证明"
    由对偶原理（Ch66A 定理 66A.6），$(A, C)$ 能观等价于 $(A^T, C^T)$ 能控。由极点配置定理（定理 66B.3），存在 $\bar{K}$ 使得 $A^T - C^T\bar{K}$ 具有期望特征值。取 $L = \bar{K}^T$，则 $A - LC = (A^T - C^T L^T)^T$ 具有相同的特征值（转置不改变特征值）。$\blacksquare$

### 分离原理

!!! theorem "定理 66B.6 (分离原理)"
    若使用基于观测器的状态反馈 $u = -K\hat{x} + r$（$\hat{x}$ 来自 Luenberger 观测器），则闭环系统的特征值为 $\sigma(A - BK) \cup \sigma(A - LC)$，即控制器极点与观测器极点的并集。

    因此，控制器和观测器可以**独立设计**。

??? proof "证明"
    闭环系统的增广状态为 $(x, e)^T$（$e = x - \hat{x}$），动力学为
    $$
    \dot{x} = Ax + B(-K\hat{x} + r) = Ax - BK(x - e) + Br = (A - BK)x + BKe + Br
    $$
    $$
    \dot{e} = (A - LC)e
    $$

    写成矩阵形式：
    $$
    \begin{pmatrix} \dot{x} \\ \dot{e} \end{pmatrix} = \begin{pmatrix} A - BK & BK \\ 0 & A - LC \end{pmatrix} \begin{pmatrix} x \\ e \end{pmatrix} + \begin{pmatrix} B \\ 0 \end{pmatrix} r
    $$

    系统矩阵是上三角分块矩阵，其特征值为对角块的特征值的并集：$\sigma(A-BK) \cup \sigma(A-LC)$。$\blacksquare$

!!! example "例 66B.3"
    对 $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$，$C = (1, 0)$。

    **控制器设计**：期望闭环极点 $-5, -5$，由例 66B.2 得 $K = (23, 7)$。

    **观测器设计**：期望观测器极点 $-10, -10$。$(A^T, C^T)$ 能控（已验证 $(A, C)$ 能观），求 $L$ 使 $A - LC$ 特征值为 $-10, -10$。

    $A - LC = \begin{pmatrix} -l_1 & 1 \\ -2-l_2 & -3 \end{pmatrix}$。特征多项式为 $\lambda^2 + (3+l_1)\lambda + (3l_1 + 2 + l_2 - 1) = \lambda^2 + 20\lambda + 100$。

    解出 $l_1 = 17$，$l_2 = 3 \cdot 17 + 2 - 1 - 100 = -48$...

    让我们重新计算：特征多项式 $= \lambda^2 + (3+l_1)\lambda + (3l_1 + 2 + l_2)$。令 $3 + l_1 = 20$，$3l_1 + 2 + l_2 = 100$。得 $l_1 = 17$，$l_2 = 100 - 51 - 2 = 47$。因此 $L = (17, 47)^T$。

---

## 66B.5 线性二次调节器 LQR

<div class="context-flow" markdown>

**核心问题**：如何设计最优的线性反馈控制律？

</div>

### 连续时间 LQR

!!! definition "定义 66B.7 (连续时间 LQR 问题)"
    **线性二次调节器**（Linear Quadratic Regulator, LQR）问题是：对系统 $\dot{x} = Ax + Bu$，求控制 $u(\cdot)$ 最小化性能指标
    $$
    J = \int_0^{\infty} \left(x(t)^T Q x(t) + u(t)^T R u(t)\right) dt
    $$
    其中 $Q \succeq 0$（状态惩罚矩阵），$R \succ 0$（控制惩罚矩阵）。

    $Q$ 惩罚状态偏离原点，$R$ 惩罚控制能量。$Q/R$ 的相对大小决定了"快速回零"与"节省能量"之间的权衡。

### 连续代数 Riccati 方程

!!! theorem "定理 66B.7 (连续 LQR 最优解)"
    若 $(A, B)$ 能稳且 $(A, Q^{1/2})$ 能检，则 LQR 问题的最优控制为线性状态反馈
    $$
    u^*(t) = -K^* x(t), \quad K^* = R^{-1} B^T P
    $$
    其中 $P \succeq 0$ 是**连续代数 Riccati 方程**（Continuous Algebraic Riccati Equation, CARE）的唯一正半定稳定化解：
    $$
    A^T P + PA - PBR^{-1}B^T P + Q = 0
    $$

    当 $(A, B)$ 能控且 $(A, Q^{1/2})$ 能观时，$P \succ 0$（严格正定）。

    最优代价为 $J^* = x_0^T P x_0$。闭环系统 $A - BK^*$ 是渐近稳定的。

??? proof "证明"
    **方法：动态规划与 Hamilton-Jacobi-Bellman 方程**

    设最优代价函数为 $V(x)$。猜测 $V(x) = x^T P x$（二次形式）。Hamilton-Jacobi-Bellman（HJB）方程为
    $$
    0 = \min_u \left\{x^T Qx + u^T Ru + \nabla V^T (Ax + Bu)\right\}
    $$

    代入 $V(x) = x^T Px$，$\nabla V = 2Px$：
    $$
    0 = \min_u \left\{x^T Qx + u^T Ru + 2x^T P(Ax + Bu)\right\}
    $$

    对 $u$ 求导并令其为零：$2Ru + 2B^T Px = 0$，得
    $$
    u^* = -R^{-1}B^T Px
    $$

    代回 HJB 方程：
    $$
    0 = x^T Qx + x^T PBR^{-1}R \cdot R^{-1}B^T Px + 2x^T PAx - 2x^T PBR^{-1}B^T Px
    $$
    $$
    = x^T Qx - x^T PBR^{-1}B^T Px + x^T PAx + x^T A^T Px
    $$
    $$
    = x^T(Q + PA + A^T P - PBR^{-1}B^T P)x
    $$

    因为对所有 $x$ 成立，得到 CARE：$A^T P + PA - PBR^{-1}B^T P + Q = 0$。

    **闭环稳定性**：以 $V(x) = x^TPx$ 为 Lyapunov 函数。沿闭环轨线：
    $$
    \dot{V} = x^T(A^T_{cl} P + PA_{cl})x
    $$
    其中 $A_{cl} = A - BK^*$。由 Riccati 方程可得
    $$
    A_{cl}^T P + PA_{cl} = A^TP + PA - 2PBR^{-1}B^TP = -Q - PBR^{-1}B^TP \preceq -Q
    $$

    当 $(A, Q^{1/2})$ 能检时，$\dot{V} = 0$ 蕴含 $Q^{1/2}x = 0$ 且 $B^TPx = 0$，由不变性原理可得 $x \to 0$。$\blacksquare$

### 离散时间 LQR

!!! definition "定义 66B.8 (离散时间 LQR)"
    对离散系统 $x_{k+1} = Ax_k + Bu_k$，最小化
    $$
    J = \sum_{k=0}^{\infty} \left(x_k^T Q x_k + u_k^T R u_k\right)
    $$

!!! theorem "定理 66B.8 (离散 LQR 最优解)"
    若 $(A, B)$ 能稳且 $(A, Q^{1/2})$ 能检，最优控制为
    $$
    u_k^* = -K^* x_k, \quad K^* = (R + B^T P B)^{-1} B^T P A
    $$

    其中 $P \succeq 0$ 是**离散代数 Riccati 方程**（Discrete ARE, DARE）的唯一正半定稳定化解：
    $$
    P = A^T P A - A^T P B(R + B^T P B)^{-1} B^T P A + Q
    $$

!!! example "例 66B.4"
    一维系统 $\dot{x} = ax + bu$，$Q = q > 0$，$R = r > 0$。

    CARE 变为标量方程 $2aP - P^2 b^2/r + q = 0$，即 $\frac{b^2}{r}P^2 - 2aP - q = 0$。

    解出 $P = \frac{ar + \sqrt{a^2r^2 + b^2rq}}{b^2} > 0$（取正根）。

    最优增益 $K^* = \frac{b}{r}P = \frac{a + \sqrt{a^2 + b^2q/r}}{b}$。

    闭环极点 $\lambda_{cl} = a - bK^* = -\sqrt{a^2 + b^2q/r} < 0$（稳定）。

---

## 66B.6 LQR 的鲁棒性

<div class="context-flow" markdown>

**核心问题**：LQR 最优控制律对系统参数的变化有多大的容忍度？

</div>

!!! theorem "定理 66B.9 (LQR 的增益和相位裕度)"
    对于连续时间 SISO LQR 系统（$m = 1$），在开环断开控制器处，系统具有以下鲁棒性裕度：

    1. **增益裕度**：$[\frac{1}{2}, +\infty)$，即控制增益可以减半到无穷大而保持稳定；
    2. **相位裕度**：至少 $60°$。

    这些裕度是 LQR 的内在性质，无需额外设计。

??? proof "证明"
    **关键工具：Kalman 恒等式（回差恒等式）**

    定义回差矩阵 $\Phi(j\omega) = I + K^*(j\omega I - A)^{-1}B$。Kalman 恒等式（有时称为 Kalman 频域不等式）为
    $$
    \Phi(j\omega)^* R \, \Phi(j\omega) = R + B^T(-j\omega I - A^T)^{-1} Q (j\omega I - A)^{-1} B
    $$

    由于右端为 $R$ 加上一个正半定矩阵，得
    $$
    \Phi(j\omega)^* R \, \Phi(j\omega) \succeq R
    $$

    对 SISO 情形（$R = r > 0$），这给出 $|\Phi(j\omega)|^2 \geq 1$，即 $|1 + L(j\omega)| \geq 1$（$L(j\omega) = K^*(j\omega I - A)^{-1}B$ 是开环传递函数）。

    $|1 + L(j\omega)| \geq 1$ 意味着 Nyquist 曲线不进入以 $(-1, 0)$ 为圆心、半径为 1 的圆。由此可推出增益裕度 $[\frac{1}{2}, +\infty)$ 和相位裕度至少 $60°$。$\blacksquare$

!!! note "注"
    MIMO LQR 也有类似的鲁棒性结果，但形式更复杂。值得注意的是，LQG（LQR + Kalman 滤波）**不**继承这些鲁棒性裕度——这是促使 $H_\infty$ 控制发展的重要动机之一。

---

## 66B.7 Kalman 滤波与 LQG

<div class="context-flow" markdown>

**核心问题**：当系统受随机噪声干扰且测量有噪声时，如何最优地估计状态？

</div>

### Kalman 滤波器

!!! definition "定义 66B.9 (随机线性系统)"
    考虑受噪声干扰的系统：
    $$
    \begin{aligned}
    \dot{x}(t) &= Ax(t) + Bu(t) + G w(t) \\
    y(t) &= Cx(t) + v(t)
    \end{aligned}
    $$

    其中 $w(t) \sim \mathcal{N}(0, W)$ 是**过程噪声**（白噪声），$v(t) \sim \mathcal{N}(0, V)$ 是**测量噪声**，$W \succeq 0$，$V \succ 0$，$w$ 和 $v$ 互不相关。

!!! theorem "定理 66B.10 (连续 Kalman 滤波器)"
    最优状态估计器（在最小均方误差意义下）为
    $$
    \dot{\hat{x}}(t) = A\hat{x}(t) + Bu(t) + L_f(y(t) - C\hat{x}(t))
    $$

    其中 Kalman 增益为
    $$
    L_f = P_f C^T V^{-1}
    $$

    $P_f \succeq 0$ 是以下 **Riccati 方程**的唯一正半定稳定化解：
    $$
    AP_f + P_f A^T - P_f C^T V^{-1} C P_f + GWG^T = 0
    $$

    $P_f$ 是稳态估计误差协方差：$P_f = \lim_{t \to \infty} \mathbb{E}[(x(t) - \hat{x}(t))(x(t) - \hat{x}(t))^T]$。

!!! theorem "定理 66B.11 (LQR-Kalman 对偶性)"
    Kalman 滤波的 Riccati 方程与 LQR 的 Riccati 方程具有精确的对偶关系：

    | LQR | Kalman 滤波 |
    |-----|-------------|
    | $A^TP + PA - PBR^{-1}B^TP + Q = 0$ | $AP_f + P_fA^T - P_fC^TV^{-1}CP_f + GWG^T = 0$ |
    | $K = R^{-1}B^TP$ | $L_f = P_fC^TV^{-1}$ |
    | $(A, B)$ 能稳 | $(A, C)$ 能检 |
    | $(A, Q^{1/2})$ 能检 | $(A, (GWG^T)^{1/2})$ 能稳 |

    Kalman 滤波等价于对偶系统 $(A^T, C^T, (GWG^T)^{1/2}, V^{1/2})$ 的 LQR 问题。

### LQG 控制

!!! definition "定义 66B.10 (LQG 控制)"
    **线性二次高斯**（Linear Quadratic Gaussian, LQG）控制将 LQR 和 Kalman 滤波结合：

    1. 用 Kalman 滤波器估计状态 $\hat{x}(t)$；
    2. 用 LQR 增益 $K^*$ 构成控制律 $u(t) = -K^* \hat{x}(t)$。

    由分离原理，闭环系统的极点为 $\sigma(A - BK^*) \cup \sigma(A - L_f C)$。

!!! note "注"
    LQG 控制器是最优的（在随机意义下），但它**不具有** LQR 那样的保证鲁棒性裕度。Doyle (1978) 的著名反例表明，LQG 闭环系统的增益裕度和相位裕度可以任意差。这一事实是 $H_\infty$ 鲁棒控制理论发展的重要动机。

---

## 66B.8 H∞ 控制

<div class="context-flow" markdown>

**核心问题**：如何在模型不确定性存在时设计保证最坏情况性能的控制器？

</div>

### H∞ 范数

!!! definition "定义 66B.11 ($H_\infty$ 范数)"
    稳定传递函数 $G(s)$ 的 **$H_\infty$ 范数**定义为
    $$
    \|G\|_{H_\infty} = \sup_{\omega \in \mathbb{R}} \bar{\sigma}(G(j\omega))
    $$

    其中 $\bar{\sigma}(\cdot)$ 表示最大奇异值。

    物理意义：$\|G\|_{H_\infty}$ 是从输入到输出的**最大增益**——在所有频率和所有输入方向上，输出能量与输入能量之比的上确界。

    等价地，$\|G\|_{H_\infty} = \sup_{u \neq 0} \frac{\|y\|_{L_2}}{\|u\|_{L_2}}$（$L_2$ 诱导范数）。

!!! theorem "定理 66B.12 ($H_\infty$ 范数与 Riccati 方程)"
    给定稳定系统 $(A, B, C, D)$ 和标量 $\gamma > 0$，则 $\|G\|_{H_\infty} < \gamma$ 当且仅当：

    1. $D^TD \prec \gamma^2 I$；
    2. 以下 Riccati 方程有正半定解 $X \succeq 0$：
    $$
    A^TX + XA + C^TC + X({\gamma^{-2}}BB^T - D_\perp)X = 0
    $$
    其中 $D_\perp$ 是与 $D$ 相关的修正项（$D = 0$ 时即为 $\gamma^{-2}BB^T$）。

### 小增益定理

!!! theorem "定理 66B.13 (小增益定理)"
    设 $G(s)$ 是稳定的传递函数，$\Delta(s)$ 是未知但稳定的不确定性，满足 $\|\Delta\|_{H_\infty} \leq 1$。则反馈互联系统
    $$
    y = G(s) (\Delta(s) y + u)
    $$
    对所有满足 $\|\Delta\|_{H_\infty} \leq 1$ 的 $\Delta$ 都内部稳定，当且仅当
    $$
    \|G\|_{H_\infty} < 1
    $$

!!! note "注"
    小增益定理是鲁棒控制的基石。直觉上，若正向通道的增益 $\|G\|$ 与反馈通道的增益 $\|\Delta\|$ 的乘积小于 1，则闭环稳定。

### H∞ 最优控制问题

!!! definition "定义 66B.12 ($H_\infty$ 最优控制问题)"
    考虑广义被控对象
    $$
    \begin{aligned}
    \dot{x} &= Ax + B_1 w + B_2 u \\
    z &= C_1 x + D_{11} w + D_{12} u \\
    y &= C_2 x + D_{21} w + D_{22} u
    \end{aligned}
    $$

    其中 $w$ 是外部干扰，$u$ 是控制输入，$z$ 是性能输出，$y$ 是可测输出。

    **$H_\infty$ 最优控制问题**：求控制器 $u = K(s) y$ 使得闭环从 $w$ 到 $z$ 的传递函数 $T_{zw}(s)$ 满足
    $$
    \|T_{zw}\|_{H_\infty} < \gamma
    $$
    并最小化 $\gamma$。

!!! theorem "定理 66B.14 ($H_\infty$ 控制的 Riccati 方程方法)"
    在标准假设下（$D_{11} = 0$，$D_{22} = 0$，$D_{12}^T D_{12} = I$，$D_{21} D_{21}^T = I$ 等），存在使 $\|T_{zw}\|_{H_\infty} < \gamma$ 的输出反馈控制器，当且仅当以下两个 Riccati 方程均有正半定稳定化解：

    **控制 Riccati 方程**：
    $$
    A^T X + XA + C_1^T C_1 + X(\gamma^{-2} B_1 B_1^T - B_2 B_2^T)X = 0, \quad X \succeq 0
    $$

    **滤波 Riccati 方程**：
    $$
    AY + YA^T + B_1 B_1^T + Y(\gamma^{-2} C_1^T C_1 - C_2^T C_2)Y = 0, \quad Y \succeq 0
    $$

    **耦合条件**：$\rho(XY) < \gamma^2$（$\rho$ 为谱半径）。

    满足条件时，一个可行的控制器（中心控制器）为
    $$
    \dot{\hat{x}} = A\hat{x} + B_2 u + (I - \gamma^{-2}YX)^{-1}(YC_2^T + B_1 D_{21}^T)(y - C_2\hat{x})
    + \gamma^{-2}(I - \gamma^{-2}YX)^{-1}YC_1^TC_1\hat{x}
    $$
    $$
    u = -B_2^T X \hat{x}
    $$

!!! note "注"
    $H_\infty$ 控制与 LQG 的关键区别：

    - LQG 假设扰动是已知统计特性的随机过程，最小化**平均**性能；
    - $H_\infty$ 控制对扰动不做统计假设，最小化**最坏情况**性能；
    - $H_\infty$ 控制可以看作一个**博弈论**问题：控制器和扰动是对弈的双方。

    当 $\gamma \to \infty$ 时，$H_\infty$ 控制退化为 LQG 控制。

---

## 66B.9 线性矩阵不等式

<div class="context-flow" markdown>

**核心问题**：能否将各种控制问题统一为凸优化问题？

</div>

### LMI 基础

!!! definition "定义 66B.13 (线性矩阵不等式)"
    **线性矩阵不等式**（Linear Matrix Inequality, LMI）是形如
    $$
    F(x) = F_0 + x_1 F_1 + x_2 F_2 + \cdots + x_k F_k \succeq 0
    $$
    的约束，其中 $x = (x_1, \ldots, x_k)$ 是决策变量，$F_i$ 是给定的对称矩阵，$\succeq 0$ 表示正半定。

    LMI 的可行域是**凸集**，因此 LMI 优化问题可以用半定规划（SDP）高效求解。

### Schur 补与 LMI

!!! theorem "定理 66B.15 (Schur 补引理)"
    设 $M = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$ 为对称矩阵，$C \succ 0$。则
    $$
    M \succeq 0 \quad \Longleftrightarrow \quad A - BC^{-1}B^T \succeq 0
    $$

    Schur 补是将**非线性矩阵不等式**转化为**线性矩阵不等式**的核心工具。

!!! example "例 66B.5"
    **Lyapunov 稳定性的 LMI 形式**：判断 $A$ 是否 Hurwitz 等价于求 $P \succ 0$ 使得
    $$
    A^TP + PA \prec 0
    $$

    这是关于 $P$ 的 LMI（$P$ 的元素是决策变量）。

!!! example "例 66B.6"
    **$H_\infty$ 范数约束的 LMI 形式**（有界实引理）：$\|G\|_{H_\infty} < \gamma$（$D = 0$）等价于存在 $P \succ 0$ 使得
    $$
    \begin{pmatrix} A^TP + PA + C^TC & PB \\ B^TP & -\gamma^2 I \end{pmatrix} \prec 0
    $$

    这是关于 $P$ 的 LMI（$\gamma$ 固定时）。

!!! example "例 66B.7"
    **状态反馈 $H_\infty$ 设计的 LMI 形式**：求 $K$ 使 $A - BK$ 稳定且 $\|T_{zw}\|_{H_\infty} < \gamma$。

    令 $Y = P^{-1}$，$Z = KY$（变量替换），问题化为关于 $(Y, Z)$ 的 LMI：
    $$
    \begin{pmatrix} AY + YA^T - B_2 Z - Z^T B_2^T + \gamma^{-2} B_1 B_1^T & (C_1 Y - D_{12} Z)^T \\ C_1 Y - D_{12} Z & -I \end{pmatrix} \prec 0
    $$
    $Y \succ 0$。求解后 $K = ZY^{-1}$。

!!! note "注"
    LMI 方法的威力在于：许多看似不同的控制问题——稳定性、$H_2$ 性能、$H_\infty$ 性能、极点区域约束、多目标优化——都可以统一表示为 LMI 约束，用半定规划一并求解。Schur 补是将非线性约束"提升"为 LMI 的关键桥梁。

---

## 本章小结

本章展示了控制系统设计中线性代数的核心作用。主要结果包括：

1. **Lyapunov 稳定性定理**将矩阵稳定性判定归结为 Lyapunov 方程 $A^TP + PA = -Q$ 的正定解的存在性（连续时间），或 $A^TPA - P = -Q$（离散时间）。

2. **能稳性和能检性**是能控性/能观性的弱化，是 LQR/LQG 设计的最小条件。

3. **极点配置定理**保证能控系统可通过状态反馈任意指定闭环极点，**Ackermann 公式**给出单输入情形的显式解。

4. **Luenberger 观测器**利用能观性设计状态估计器，**分离原理**允许控制器和观测器独立设计，闭环极点为两者的并集。

5. **LQR** 将最优控制归结为代数 Riccati 方程，最优解是线性反馈。连续 LQR 具有至少 $60°$ 相位裕度和 $[\frac{1}{2}, +\infty)$ 增益裕度。

6. **离散 LQR** 通过离散 Riccati 方程求解，结构与连续情形对偶。

7. **Kalman 滤波**是 LQR 的对偶问题，**LQG** 将两者结合但不继承 LQR 的鲁棒性裕度。

8. **$H_\infty$ 控制**通过最小化最坏情况增益处理模型不确定性，小增益定理是其理论基础，两个耦合 Riccati 方程给出控制器存在的充要条件。

9. **LMI** 通过 Schur 补等工具将控制问题统一为凸优化，是现代鲁棒控制的计算框架。

---

## 习题

!!! question "习题 66B.1"
    设 $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$。

    (a) 判断 $A$ 是否 Hurwitz 稳定。

    (b) 取 $Q = I$，求解连续 Lyapunov 方程 $A^TP + PA = -I$。讨论解的正定性。

!!! question "习题 66B.2"
    设 $A = \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$。验证系统不能控但能稳。设计 $K$ 使 $A - BK$ Hurwitz 稳定。

!!! question "习题 66B.3"
    设 $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$。设计状态反馈 $u = -Kx$ 使闭环极点为 $-1 \pm j$。

!!! question "习题 66B.4"
    对一维系统 $\dot{x} = 2x + u$，$Q = 1$，$R = 1$，求解 LQR 问题（求 $P$、$K^*$ 和闭环极点）。

!!! question "习题 66B.5"
    证明 LQR 闭环系统 $A - BK^*$ 是渐近稳定的。（提示：用 $V(x) = x^TPx$ 作为 Lyapunov 函数，利用 Riccati 方程。）

!!! question "习题 66B.6"
    证明：若 $(A, B)$ 不能控，则存在 $A$ 的特征值不能通过状态反馈改变。更精确地，不能控的模态对应的特征值在 $\sigma(A - BK)$ 中固定不变。

!!! question "习题 66B.7"
    设计 Luenberger 观测器，使观测器误差动力学 $\dot{e} = (A - LC)e$ 的特征值为 $-10, -10$，其中 $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}$，$C = (1, 0)$。

!!! question "习题 66B.8"
    对于离散时间系统 $x_{k+1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}x_k + \begin{pmatrix} 0 \\ 1 \end{pmatrix}u_k$，$Q = I$，$R = 1$：

    (a) 验证系统的能控性。

    (b) 写出离散代数 Riccati 方程（DARE）。

    (c) 尝试迭代求解 DARE（从 $P_0 = 0$ 开始）。

!!! question "习题 66B.9"
    证明 Kalman 滤波的 Riccati 方程可以通过对偶变换从 LQR 的 Riccati 方程得到。具体地，说明对偶关系中各矩阵的对应。

!!! question "习题 66B.10"
    考虑系统 $G(s) = \frac{1}{s-1}$（不稳定被控对象）。

    (a) 设计比例控制器 $u = -ky$ 使闭环稳定。求 $k$ 的范围。

    (b) 若存在乘性不确定性 $G_\Delta(s) = G(s)(1 + \Delta(s))$，$\|\Delta\|_{H_\infty} \leq 0.5$，讨论闭环的鲁棒稳定性条件。

!!! question "习题 66B.11"
    将 Lyapunov 稳定性条件 $A^TP + PA \prec 0$，$P \succ 0$ 写成标准 LMI 形式 $F_0 + \sum_i x_i F_i \succ 0$，其中 $x_i$ 是 $P$ 的独立元素。对 $2 \times 2$ 的情形显式写出。

!!! question "习题 66B.12"
    利用 Schur 补引理证明：对称矩阵 $\begin{pmatrix} P & Q \\ Q^T & R \end{pmatrix} \succeq 0$ 且 $R \succ 0$ 等价于 $P - QR^{-1}Q^T \succeq 0$。

!!! question "习题 66B.13"
    设 $A$ Hurwitz 稳定，$G(s) = C(sI - A)^{-1}B$（$D = 0$）。利用有界实引理证明 $\|G\|_{H_\infty} < \gamma$ 等价于存在 $P \succ 0$ 使得
    $$
    A^TP + PA + \gamma^{-2}PBB^TP + C^TC \prec 0
    $$
    并说明这与 Riccati 方程方法的关系。

!!! question "习题 66B.14"
    （综合题）考虑倒立摆系统，线性化后 $A = \begin{pmatrix} 0 & 1 \\ g/l & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ -1/(ml) \end{pmatrix}$，$C = (1, 0)$。取 $g = 9.8$，$l = 1$，$m = 1$。

    (a) 验证系统不稳定但能控能观。

    (b) 用 LQR 设计控制器（$Q = \operatorname{diag}(10, 1)$，$R = 1$），求 Riccati 方程的解 $P$ 和最优增益 $K^*$。

    (c) 设计 Kalman 滤波器（$W = 1$，$V = 0.1$），求滤波增益 $L_f$。

    (d) 讨论 LQG 闭环系统的极点。
