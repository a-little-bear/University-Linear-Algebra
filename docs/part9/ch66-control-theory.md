# 第 66 章 线性代数在控制理论中的应用

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵指数(Ch13) · 矩阵方程(Ch20) · 矩阵稳定性(Ch36)

**本章脉络**：状态空间模型 → 能控性(Kalman 秩条件) → 能观性 → 对偶性 → Kalman 分解 → 极点配置 → 状态观测器 → 最优控制(LQR) → H∞ 控制初步

**延伸**：控制理论的线性代数框架直接推广到非线性系统的局部分析（线性化）、分布参数系统（无穷维算子半群）和随机控制（Kalman 滤波的矩阵 Riccati 方程）

</div>

控制理论是线性代数最成功的应用领域之一。1960 年，Kalman 建立的状态空间方法彻底革新了控制系统的分析与设计，将矩阵理论置于控制工程的核心。在这一框架中，系统的动力学由矩阵 $(A, B, C, D)$ 完全描述，系统的基本性质——能控性、能观性、稳定性——都可以用矩阵的秩条件和特征值来刻画。极点配置定理保证了能控系统可以通过线性反馈任意配置闭环极点，而线性二次调节器（LQR）将最优控制问题归结为代数 Riccati 方程。

本章系统地展示这些经典控制理论结果背后的线性代数结构。

---

## 66.1 状态空间模型

<div class="context-flow" markdown>

**核心问题**：如何用矩阵语言统一描述线性动力系统？

</div>

### 连续时间模型

!!! definition "定义 66.1 (连续时间线性时不变系统)"
    **连续时间线性时不变**（continuous-time linear time-invariant, CT-LTI）系统由以下方程描述：
    $$
    \begin{aligned}
    \dot{x}(t) &= Ax(t) + Bu(t) \\
    y(t) &= Cx(t) + Du(t)
    \end{aligned}
    $$

    其中：

    - $x(t) \in \mathbb{R}^n$ 是**状态向量**（state vector）；
    - $u(t) \in \mathbb{R}^m$ 是**控制输入**（control input）；
    - $y(t) \in \mathbb{R}^p$ 是**输出**（output）；
    - $A \in \mathbb{R}^{n \times n}$ 是**系统矩阵**（system matrix）；
    - $B \in \mathbb{R}^{n \times m}$ 是**输入矩阵**（input matrix）；
    - $C \in \mathbb{R}^{p \times n}$ 是**输出矩阵**（output matrix）；
    - $D \in \mathbb{R}^{p \times m}$ 是**直接传递矩阵**（feedthrough matrix）。

### 解的公式

!!! theorem "定理 66.1 (状态方程的解)"
    给定初始条件 $x(0) = x_0$，系统 $\dot{x} = Ax + Bu$ 的解为
    $$
    x(t) = e^{At} x_0 + \int_0^t e^{A(t-\tau)} B u(\tau) \, d\tau
    $$

    相应的输出为
    $$
    y(t) = C e^{At} x_0 + \int_0^t C e^{A(t-\tau)} B u(\tau) \, d\tau + D u(t)
    $$

??? proof "证明"
    定义 $z(t) = e^{-At} x(t)$。则
    $$
    \dot{z}(t) = -A e^{-At} x(t) + e^{-At} \dot{x}(t) = -A e^{-At} x + e^{-At}(Ax + Bu) = e^{-At} Bu(t)
    $$
    积分得 $z(t) = z(0) + \int_0^t e^{-A\tau} Bu(\tau) d\tau = x_0 + \int_0^t e^{-A\tau} Bu(\tau) d\tau$。
    乘以 $e^{At}$ 即得结果。

### 传递函数

!!! definition "定义 66.2 (传递函数)"
    系统 $(A, B, C, D)$ 的**传递函数**（transfer function）是 $p \times m$ 有理矩阵
    $$
    G(s) = C(sI - A)^{-1}B + D
    $$

    传递函数将输入的 Laplace 变换 $U(s)$ 映射到输出的 Laplace 变换 $Y(s)$：$Y(s) = G(s) U(s)$（假设零初始条件）。

!!! example "例 66.1"
    **弹簧-质量-阻尼系统**：$m\ddot{q} + c\dot{q} + kq = f(t)$。

    取状态 $x = (q, \dot{q})^T$，输入 $u = f$，输出 $y = q$：
    $$
    A = \begin{pmatrix} 0 & 1 \\ -k/m & -c/m \end{pmatrix}, \;
    B = \begin{pmatrix} 0 \\ 1/m \end{pmatrix}, \;
    C = \begin{pmatrix} 1 & 0 \end{pmatrix}, \;
    D = 0
    $$

    传递函数：$G(s) = \frac{1/m}{s^2 + (c/m)s + k/m} = \frac{1}{ms^2 + cs + k}$。

### 离散时间模型

!!! definition "定义 66.3 (离散时间线性时不变系统)"
    **离散时间 LTI 系统**为
    $$
    \begin{aligned}
    x_{k+1} &= A x_k + B u_k \\
    y_k &= C x_k + D u_k
    \end{aligned}
    $$

    解为 $x_k = A^k x_0 + \sum_{j=0}^{k-1} A^{k-1-j} B u_j$。

---

## 66.2 能控性

<div class="context-flow" markdown>

**核心问题**：什么条件下，系统的状态可以通过控制输入从任意初态转移到任意终态？

</div>

### 能控性矩阵

!!! definition "定义 66.4 (能控性)"
    系统 $(A, B)$ 称为**能控的**（controllable），若对任意初始状态 $x_0$ 和目标状态 $x_f$，存在有限时间 $T > 0$ 和控制输入 $u(\cdot)$，使得系统从 $x(0) = x_0$ 到达 $x(T) = x_f$。

!!! definition "定义 66.5 (能控性矩阵)"
    **能控性矩阵**（controllability matrix）定义为
    $$
    \mathcal{C} = \begin{pmatrix} B & AB & A^2B & \cdots & A^{n-1}B \end{pmatrix} \in \mathbb{R}^{n \times nm}
    $$

!!! theorem "定理 66.2 (Kalman 能控性秩条件)"
    系统 $(A, B)$ 能控，当且仅当能控性矩阵满秩：
    $$
    \operatorname{rank}(\mathcal{C}) = \operatorname{rank}\begin{pmatrix} B & AB & A^2B & \cdots & A^{n-1}B \end{pmatrix} = n
    $$

??? proof "证明"
    **充分性**：设 $\operatorname{rank}(\mathcal{C}) = n$。对于连续时间系统，可达集（从原点出发在时间 $T$ 内可以到达的状态集合）为
    $$
    \mathcal{R}_T = \left\{\int_0^T e^{A(T-\tau)} Bu(\tau) d\tau : u \in L^2([0,T]; \mathbb{R}^m)\right\}
    $$

    可以证明 $\mathcal{R}_T = \operatorname{Im}(\mathcal{C})$。因此 $\operatorname{rank}(\mathcal{C}) = n$ 蕴含 $\mathcal{R}_T = \mathbb{R}^n$。

    具体地，$e^{A(T-\tau)}B$ 的列空间包含在 $\operatorname{Im}(\mathcal{C})$ 中（因为 $e^{At} = \sum_{k=0}^{\infty} \frac{(At)^k}{k!}$，而由 Cayley-Hamilton 定理 $A^k$ 对 $k \geq n$ 可以用 $I, A, \ldots, A^{n-1}$ 线性表示）。

    **必要性**：若 $\operatorname{rank}(\mathcal{C}) < n$，则存在非零 $v$ 使得 $v^T A^k B = 0$ 对所有 $k = 0, 1, \ldots, n-1$。由 Cayley-Hamilton 定理，$v^T A^k B = 0$ 对所有 $k \geq 0$。从而
    $$
    v^T x(T) = v^T e^{AT} x_0 + \int_0^T v^T e^{A(T-\tau)} Bu(\tau) d\tau = v^T e^{AT} x_0
    $$
    与输入 $u$ 无关。因此从 $x_0 = 0$ 出发无法到达满足 $v^T x_f \neq 0$ 的状态 $x_f$。

### PBH 判据

!!! theorem "定理 66.3 (Popov-Belevitch-Hautus (PBH) 判据)"
    系统 $(A, B)$ 能控，当且仅当对 $A$ 的每个特征值 $\lambda$，
    $$
    \operatorname{rank}\begin{pmatrix} A - \lambda I & B \end{pmatrix} = n
    $$

    等价地，不存在 $A^T$ 的左特征向量 $v$（即 $A^T v = \lambda v$，$v \neq 0$）使得 $B^T v = 0$。

??? proof "证明"
    **$\Rightarrow$**：若 $\operatorname{rank}(A - \lambda I, B) < n$，则存在 $v \neq 0$ 使得 $v^T(A - \lambda I) = 0$ 且 $v^T B = 0$。即 $v^T A = \lambda v^T$，从而 $v^T A^k B = \lambda^k v^T B = 0$ 对所有 $k$。因此 $v^T \mathcal{C} = 0$，$\operatorname{rank}(\mathcal{C}) < n$。

    **$\Leftarrow$**：若 $\operatorname{rank}(\mathcal{C}) < n$，则存在 $v \neq 0$ 使得 $v^T \mathcal{C} = 0$。设 $V = \operatorname{Im}(\mathcal{C})$，则 $V$ 是 $A$-不变的（$A \cdot A^k B = A^{k+1}B \in V$）。$A|_V$ 的某个特征值 $\lambda$ 对应的特征向量 $w$ 满足 $A^T w = \lambda w$ 且 $w \perp V$，故 $B^Tw = 0$。

!!! example "例 66.2"
    考虑系统 $A = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$。

    能控性矩阵：$\mathcal{C} = \begin{pmatrix} B & AB \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & 2 \end{pmatrix}$，$\det(\mathcal{C}) = -1 \neq 0$。系统能控。

    PBH 检验：$\lambda_1 = 1$，$\operatorname{rank}\begin{pmatrix} 0 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix} = 2$。$\lambda_2 = 2$，$\operatorname{rank}\begin{pmatrix} -1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = 2$。均满秩，能控。

### 能控标准型

!!! definition "定义 66.6 (能控标准型)"
    若单输入系统 $(A, b)$（$b \in \mathbb{R}^n$）能控，则存在坐标变换将其化为**能控标准型**（controllable canonical form）：
    $$
    \bar{A} = \begin{pmatrix} 0 & 1 & 0 & \cdots & 0 \\ 0 & 0 & 1 & \cdots & 0 \\ \vdots & & & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & 1 \\ -a_0 & -a_1 & -a_2 & \cdots & -a_{n-1} \end{pmatrix}, \quad
    \bar{b} = \begin{pmatrix} 0 \\ 0 \\ \vdots \\ 0 \\ 1 \end{pmatrix}
    $$
    其中 $\chi_A(\lambda) = \lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_1\lambda + a_0$ 是 $A$ 的特征多项式。

---

## 66.3 能观性

<div class="context-flow" markdown>

**核心问题**：从输出测量能否唯一确定系统的初始状态？

</div>

### 能观性矩阵

!!! definition "定义 66.7 (能观性)"
    系统 $(A, C)$ 称为**能观的**（observable），若初始状态 $x_0$ 可以从输出 $y(t)$（$t \in [0, T]$）和输入 $u(t)$ 中唯一确定。

!!! definition "定义 66.8 (能观性矩阵)"
    **能观性矩阵**（observability matrix）定义为
    $$
    \mathcal{O} = \begin{pmatrix} C \\ CA \\ CA^2 \\ \vdots \\ CA^{n-1} \end{pmatrix} \in \mathbb{R}^{np \times n}
    $$

!!! theorem "定理 66.4 (Kalman 能观性秩条件)"
    系统 $(A, C)$ 能观，当且仅当
    $$
    \operatorname{rank}(\mathcal{O}) = \operatorname{rank}\begin{pmatrix} C \\ CA \\ \vdots \\ CA^{n-1} \end{pmatrix} = n
    $$

??? proof "证明"
    对于零输入响应（$u = 0$），$y(t) = Ce^{At}x_0$。

    **充分性**：若 $\operatorname{rank}(\mathcal{O}) = n$，则由 $y(t) = Ce^{At}x_0$，$\dot{y}(t) = CAe^{At}x_0$，...，$y^{(k)}(t) = CA^k e^{At}x_0$。在 $t = 0$ 处取值：
    $$
    \begin{pmatrix} y(0) \\ \dot{y}(0) \\ \vdots \\ y^{(n-1)}(0) \end{pmatrix} = \mathcal{O} \cdot x_0
    $$
    由 $\operatorname{rank}(\mathcal{O}) = n$，$x_0$ 被唯一确定。

    **必要性**：若 $\operatorname{rank}(\mathcal{O}) < n$，存在 $v \neq 0$ 使得 $\mathcal{O}v = 0$，即 $CA^k v = 0$ 对 $k = 0, \ldots, n-1$。由 Cayley-Hamilton，$CA^k v = 0$ 对所有 $k$。因此 $Ce^{At}v = 0$ 对所有 $t$，无法区分 $x_0 = 0$ 和 $x_0 = v$。

### PBH 能观性判据

!!! theorem "定理 66.5 (PBH 能观性判据)"
    $(A, C)$ 能观，当且仅当对 $A$ 的每个特征值 $\lambda$，
    $$
    \operatorname{rank}\begin{pmatrix} A - \lambda I \\ C \end{pmatrix} = n
    $$

!!! example "例 66.3"
    系统 $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}$，$C = \begin{pmatrix} 1 & 0 \end{pmatrix}$。

    $\mathcal{O} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I_2$，满秩。系统能观。

    直观理解：$y = x_1$（直接测量第一个状态），$\dot{y} = x_2$（通过导数得到第二个状态）。

---

## 66.4 能控能观对偶性

<div class="context-flow" markdown>

**核心问题**：能控性和能观性之间有什么对称关系？

</div>

!!! theorem "定理 66.6 (对偶原理)"
    $(A, B)$ 能控当且仅当 $(A^T, B^T)$ 能观。

    等价地，$(A, C)$ 能观当且仅当 $(A^T, C^T)$ 能控。

??? proof "证明"
    $(A, B)$ 的能控性矩阵为 $\mathcal{C} = (B, AB, \ldots, A^{n-1}B)$。

    $(A^T, B^T)$ 的能观性矩阵为
    $$
    \mathcal{O}' = \begin{pmatrix} B^T \\ B^T A^T \\ \vdots \\ B^T (A^T)^{n-1} \end{pmatrix} = \begin{pmatrix} B^T \\ (AB)^T \\ \vdots \\ (A^{n-1}B)^T \end{pmatrix} = \mathcal{C}^T
    $$

    因此 $\operatorname{rank}(\mathcal{C}) = \operatorname{rank}(\mathcal{C}^T) = \operatorname{rank}(\mathcal{O}')$。

!!! note "注"
    对偶原理将能控性问题和能观性问题统一起来：任何关于能控性的定理，通过转置可以立即得到关于能观性的对应定理。这种对偶性在控制系统设计中也有深刻的工程含义——控制器设计和观测器设计本质上是对偶问题。

---

## 66.5 Kalman 分解

<div class="context-flow" markdown>

**核心问题**：如何将一般系统分解为能控/不能控和能观/不能观的子系统？

</div>

### Kalman 典则分解

!!! theorem "定理 66.7 (Kalman 分解)"
    对任意系统 $(A, B, C, D)$，存在坐标变换 $T$ 使得系统在新坐标下具有分块形式：
    $$
    \bar{A} = T^{-1}AT = \begin{pmatrix} A_{co} & 0 & A_{13} & 0 \\ A_{21} & A_{\bar{c}o} & A_{23} & A_{24} \\ 0 & 0 & A_{c\bar{o}} & 0 \\ 0 & 0 & A_{43} & A_{\bar{c}\bar{o}} \end{pmatrix}
    $$
    $$
    \bar{B} = T^{-1}B = \begin{pmatrix} B_{co} \\ 0 \\ B_{c\bar{o}} \\ 0 \end{pmatrix}, \quad
    \bar{C} = CT = \begin{pmatrix} C_{co} & C_{\bar{c}o} & 0 & 0 \end{pmatrix}
    $$

    其中下标含义为：

    - $co$：能控且能观（controllable and observable）
    - $\bar{c}o$：不能控但能观
    - $c\bar{o}$：能控但不能观
    - $\bar{c}\bar{o}$：不能控且不能观

!!! note "注"
    Kalman 分解揭示了一个深刻的事实：只有**能控且能观**的子系统才能从输入-输出关系中完全辨识和控制。传递函数 $G(s) = C(sI-A)^{-1}B + D$ 实际上只取决于 $A_{co}$ 部分。

### 最小实现

!!! definition "定义 66.9 (最小实现)"
    给定传递函数 $G(s)$，其**最小实现**（minimal realization）$(A, B, C, D)$ 是维数最小的状态空间实现。

!!! theorem "定理 66.8 (最小实现的刻画)"
    实现 $(A, B, C, D)$ 是最小的，当且仅当 $(A, B)$ 能控且 $(A, C)$ 能观。最小实现的状态维度等于 $G(s)$ 的 McMillan 度。

---

## 66.6 极点配置

<div class="context-flow" markdown>

**核心问题**：通过线性状态反馈，能否任意指定闭环系统的特征值？

</div>

### 状态反馈

!!! definition "定义 66.10 (线性状态反馈)"
    **线性状态反馈**（linear state feedback）是控制律
    $$
    u(t) = -Kx(t) + r(t)
    $$
    其中 $K \in \mathbb{R}^{m \times n}$ 是**反馈增益矩阵**，$r(t)$ 是参考输入。

    闭环系统为 $\dot{x} = (A - BK)x + Br$，闭环系统矩阵为 $A_{cl} = A - BK$。

### 极点配置定理

!!! theorem "定理 66.9 (极点配置定理)"
    若 $(A, B)$ 能控，则对任意对称的（关于实轴对称的）$n$ 个复数集合 $\{\lambda_1, \ldots, \lambda_n\}$，存在实矩阵 $K \in \mathbb{R}^{m \times n}$，使得 $A - BK$ 的特征值恰好为 $\{\lambda_1, \ldots, \lambda_n\}$。

??? proof "证明"
    **单输入情形**（$m = 1$，$B = b$）：

    将系统化为能控标准型 $(\bar{A}, \bar{b})$。在标准型中，取 $\bar{K} = (\bar{k}_0, \bar{k}_1, \ldots, \bar{k}_{n-1})$，则
    $$
    \bar{A} - \bar{b}\bar{K} = \begin{pmatrix} 0 & 1 & 0 & \cdots & 0 \\ 0 & 0 & 1 & \cdots & 0 \\ \vdots & & & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & 1 \\ -(a_0+\bar{k}_0) & -(a_1+\bar{k}_1) & \cdots & \cdots & -(a_{n-1}+\bar{k}_{n-1}) \end{pmatrix}
    $$
    其特征多项式为 $\lambda^n + (a_{n-1}+\bar{k}_{n-1})\lambda^{n-1} + \cdots + (a_0 + \bar{k}_0)$。

    要使特征值为 $\lambda_1, \ldots, \lambda_n$，令特征多项式等于 $\prod_{i=1}^n (\lambda - \lambda_i) = \lambda^n + \alpha_{n-1}\lambda^{n-1} + \cdots + \alpha_0$。

    解出 $\bar{k}_i = \alpha_i - a_i$，$i = 0, \ldots, n-1$。然后 $K = \bar{K} T^{-1}$（$T$ 是标准型变换矩阵）。

    **多输入情形**的证明更技术化，需要利用 Heymann 引理将多输入问题归约为单输入问题。

### Ackermann 公式

!!! theorem "定理 66.10 (Ackermann 公式)"
    对于单输入能控系统 $(A, b)$，使闭环特征多项式为 $\alpha(\lambda) = \prod_{i=1}^n (\lambda - \lambda_i)$ 的反馈增益为
    $$
    K = e_n^T \mathcal{C}^{-1} \alpha(A)
    $$
    其中 $e_n = (0, \ldots, 0, 1)^T$，$\mathcal{C} = (b, Ab, \ldots, A^{n-1}b)$ 是能控性矩阵，$\alpha(A) = A^n + \alpha_{n-1}A^{n-1} + \cdots + \alpha_0 I$。

!!! example "例 66.4"
    设 $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$。期望闭环极点为 $\lambda_1 = -5, \lambda_2 = -5$。

    $\mathcal{C} = \begin{pmatrix} 0 & 1 \\ 1 & -3 \end{pmatrix}$，$\det(\mathcal{C}) = -1$，能控。

    期望特征多项式：$\alpha(\lambda) = (\lambda + 5)^2 = \lambda^2 + 10\lambda + 25$。

    原特征多项式：$\chi_A(\lambda) = \lambda^2 + 3\lambda + 2$。

    取 $K = (k_1, k_2)$，使 $(A - BK)$ 的特征多项式为 $\lambda^2 + (3+k_2)\lambda + (2+k_1) = \lambda^2 + 10\lambda + 25$。

    解出 $k_2 = 7$，$k_1 = 23$。因此 $K = (23, 7)$。

---

## 66.7 状态观测器

<div class="context-flow" markdown>

**核心问题**：当状态不能直接测量时，如何从输出重构状态？

</div>

### Luenberger 观测器

!!! definition "定义 66.11 (Luenberger 观测器)"
    **Luenberger 观测器**是估计系统状态的动态系统：
    $$
    \dot{\hat{x}}(t) = A\hat{x}(t) + Bu(t) + L(y(t) - C\hat{x}(t))
    $$
    其中 $\hat{x}(t)$ 是状态估计，$L \in \mathbb{R}^{n \times p}$ 是**观测器增益**。

    定义估计误差 $e(t) = x(t) - \hat{x}(t)$，则
    $$
    \dot{e}(t) = (A - LC)e(t)
    $$

!!! theorem "定理 66.11 (观测器极点配置)"
    若 $(A, C)$ 能观，则对任意期望的极点集合，存在 $L$ 使得 $A - LC$ 的特征值为给定值。

??? proof "证明"
    由对偶原理（定理 66.6），$(A, C)$ 能观等价于 $(A^T, C^T)$ 能控。由极点配置定理（定理 66.9），存在 $\bar{K}$ 使得 $A^T - C^T\bar{K}$ 具有期望特征值。取 $L = \bar{K}^T$，则 $A - LC = (A^T - C^T L^T)^T$ 具有相同的特征值。

### 分离原理

!!! theorem "定理 66.12 (分离原理)"
    若使用状态反馈 $u = -K\hat{x}$（$\hat{x}$ 来自 Luenberger 观测器），则闭环系统的特征值为 $A - BK$ 的特征值（控制器极点）与 $A - LC$ 的特征值（观测器极点）的并集。

    因此，控制器和观测器可以**独立设计**。

??? proof "证明"
    闭环系统的增广状态为 $(x, e)^T$，动力学为
    $$
    \begin{pmatrix} \dot{x} \\ \dot{e} \end{pmatrix} = \begin{pmatrix} A - BK & BK \\ 0 & A - LC \end{pmatrix} \begin{pmatrix} x \\ e \end{pmatrix}
    $$
    系统矩阵是上三角分块的，特征值为 $\sigma(A-BK) \cup \sigma(A-LC)$。

---

## 66.8 线性二次调节器 LQR

<div class="context-flow" markdown>

**核心问题**：如何设计最优的线性反馈控制律？

</div>

### LQR 问题

!!! definition "定义 66.12 (LQR 问题)"
    **线性二次调节器**（Linear Quadratic Regulator, LQR）问题是：对系统 $\dot{x} = Ax + Bu$，求控制 $u(\cdot)$ 最小化性能指标
    $$
    J = \int_0^{\infty} \left(x(t)^T Q x(t) + u(t)^T R u(t)\right) dt
    $$
    其中 $Q \succeq 0$（状态惩罚矩阵），$R \succ 0$（控制惩罚矩阵）。

### 代数 Riccati 方程

!!! theorem "定理 66.13 (LQR 最优解)"
    若 $(A, B)$ 能控且 $(A, Q^{1/2})$ 能观，则 LQR 问题的最优控制为线性状态反馈
    $$
    u^*(t) = -K^* x(t), \quad K^* = R^{-1} B^T P
    $$
    其中 $P \succ 0$ 是**代数 Riccati 方程**（Algebraic Riccati Equation, ARE）的唯一正定解：
    $$
    A^T P + PA - PBR^{-1}B^T P + Q = 0
    $$

    最优代价为 $J^* = x_0^T P x_0$。

??? proof "证明"
    使用动态规划。设最优代价函数为 $V(x) = x^T P x$（猜测二次形式）。Hamilton-Jacobi-Bellman (HJB) 方程为
    $$
    0 = \min_u \left\{x^T Qx + u^T Ru + \nabla V^T (Ax + Bu)\right\}
    $$

    代入 $V(x) = x^T Px$，$\nabla V = 2Px$：
    $$
    0 = \min_u \left\{x^T Qx + u^T Ru + 2x^T P(Ax + Bu)\right\}
    $$

    对 $u$ 求导令其为零：$2Ru + 2B^T Px = 0$，得 $u^* = -R^{-1}B^T Px$。

    代回 HJB 方程：
    $$
    0 = x^T Qx + x^T PBR^{-1}B^T Px + 2x^T PAx - 2x^T PBR^{-1}B^T Px
    $$
    $$
    = x^T(Q + PA + A^T P - PBR^{-1}B^T P)x
    $$
    因为对所有 $x$ 成立，得到 ARE：$A^T P + PA - PBR^{-1}B^T P + Q = 0$。

!!! example "例 66.5"
    一维系统 $\dot{x} = ax + bu$，$Q = q > 0$，$R = r > 0$。

    ARE 变为 $2aP - P^2 b^2/r + q = 0$，即 $\frac{b^2}{r}P^2 - 2aP - q = 0$。

    解出 $P = \frac{ar + \sqrt{a^2r^2 + b^2rq}}{b^2} > 0$（取正根）。

    最优增益 $K^* = \frac{b}{r}P = \frac{a + \sqrt{a^2 + b^2q/r}}{b}$。

    闭环极点 $\lambda_{cl} = a - bK^* = -\sqrt{a^2 + b^2q/r} < 0$（稳定）。

!!! note "注"
    LQR 问题是控制理论中线性代数最优美的应用之一：

    - 最优控制律是**线性**状态反馈，增益由 Riccati 方程的解确定；
    - Riccati 方程是**矩阵方程**（Ch20 的内容）；
    - 闭环系统的稳定性由 Riccati 方程的正定解保证；
    - LQR 控制律具有极好的鲁棒性：增益裕度为无穷大，相位裕度至少为 $60°$。

---

## 本章小结

本章展示了线性代数在控制理论中的核心作用。主要结果包括：

1. **状态空间模型** $(A, B, C, D)$ 用四个矩阵完整描述线性时不变系统，解由矩阵指数 $e^{At}$ 给出。

2. **能控性**由 Kalman 秩条件（$\operatorname{rank}(\mathcal{C}) = n$）或 PBH 判据刻画，保证可以通过输入将状态转移到任意位置。

3. **能观性**由 $\operatorname{rank}(\mathcal{O}) = n$ 刻画，保证可以从输出重构状态。

4. **对偶原理**统一了能控性和能观性：$(A, B)$ 能控等价于 $(A^T, B^T)$ 能观。

5. **Kalman 分解**将系统分为四个子系统，揭示了只有能控且能观的部分才能从输入-输出关系中辨识。

6. **极点配置**定理保证能控系统可通过状态反馈任意指定闭环极点。

7. **Luenberger 观测器**利用能观性设计状态估计器，分离原理允许独立设计控制器和观测器。

8. **LQR** 将最优控制归结为代数 Riccati 方程，最优解是线性反馈。

---

## 习题

!!! question "习题 66.1"
    验证系统 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$ 的能控性。

!!! question "习题 66.2"
    设 $A = \begin{pmatrix} \lambda_1 & 0 \\ 0 & \lambda_2 \end{pmatrix}$（$\lambda_1 \neq \lambda_2$），$B = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}$。证明 $(A, B)$ 能控当且仅当 $b_1 \neq 0$ 且 $b_2 \neq 0$。

!!! question "习题 66.3"
    对系统 $A = \begin{pmatrix} -1 & 0 \\ 0 & -2 \end{pmatrix}$，$C = (1, 1)$，判断能观性。

!!! question "习题 66.4"
    验证对偶原理：对例 66.2 中的系统，构造对偶系统 $(A^T, B^T)$ 并验证其能观性。

!!! question "习题 66.5"
    设 $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$。设计状态反馈 $u = -Kx$ 使闭环极点为 $-1 \pm j$。

!!! question "习题 66.6"
    对一维系统 $\dot{x} = 2x + u$，$Q = 1$，$R = 1$，求解 LQR 问题（求 $P$、$K^*$ 和闭环极点）。

!!! question "习题 66.7"
    证明：若 $(A, B)$ 不能控，则存在特征值不能通过状态反馈改变。

!!! question "习题 66.8"
    设计 Luenberger 观测器，使观测器误差动力学 $\dot{e} = (A - LC)e$ 的特征值为 $-10, -10$，其中 $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}$，$C = (1, 0)$。

!!! question "习题 66.9"
    证明 LQR 闭环系统 $A - BK^*$ 是渐近稳定的（提示：用 $V(x) = x^TPx$ 作为 Lyapunov 函数）。

!!! question "习题 66.10"
    对于离散时间系统 $x_{k+1} = Ax_k + Bu_k$，写出离散时间能控性矩阵和 Kalman 秩条件。
