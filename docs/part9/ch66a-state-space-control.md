# 第 66A 章 状态空间与系统结构

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵指数(Ch13) · 矩阵方程(Ch20) · 矩阵稳定性(Ch36) · 奇异值分解(Ch7)

**本章脉络**：状态空间模型 → 矩阵指数解 → 能控性(Kalman 秩条件 / PBH 判据) → 能观性 → 对偶性 → Gramian → Kalman 分解 → 最小实现 → Ho-Kalman 算法 → 平衡实现与模型降阶

**延伸**：能控性/能观性 Gramian 在 Ch66B 的 LQR/LQG 设计中不可或缺；Kalman 分解揭示传递函数的本质维度，平衡截断是现代大规模系统降阶的主力方法

</div>

状态空间方法是 Kalman 于 1960 年建立的现代控制理论基石。一个线性时不变系统由四元组 $(A, B, C, D)$ 完全描述，系统的动力学解由矩阵指数 $e^{At}$ 给出。系统的两个最基本结构性质——能控性与能观性——分别回答"能否通过输入将状态驱动到任意位置"和"能否通过输出唯一确定系统状态"。Kalman 秩条件和 PBH 判据从不同角度刻画这两个性质，而能控性/能观性 Gramian 则提供了能量层面的定量描述。Kalman 典则分解将系统分为四个子系统，揭示只有能控且能观的部分才出现在传递函数中。Ho-Kalman 算法从传递函数恢复最小实现，平衡实现和平衡截断将奇异值分解引入系统降阶。

本章系统地展示这些经典结果背后的线性代数结构。

---

## 66A.1 状态空间模型

<div class="context-flow" markdown>

**核心问题**：如何用矩阵语言统一描述线性动力系统？

</div>

### 连续时间模型

!!! definition "定义 66A.1 (连续时间线性时不变系统)"
    **连续时间线性时不变**（continuous-time linear time-invariant, CT-LTI）系统由以下方程描述：
    $$
    \begin{aligned}
    \dot{x}(t) &= Ax(t) + Bu(t) \\
    y(t) &= Cx(t) + Du(t)
    \end{aligned}
    $$

    其中：

    - $x(t) \in \mathbb{R}^n$ 是**状态向量**（state vector），$n$ 称为系统的**阶**或**状态维度**；
    - $u(t) \in \mathbb{R}^m$ 是**控制输入**（control input）；
    - $y(t) \in \mathbb{R}^p$ 是**输出**（output）；
    - $A \in \mathbb{R}^{n \times n}$ 是**系统矩阵**（system matrix），决定自由响应；
    - $B \in \mathbb{R}^{n \times m}$ 是**输入矩阵**（input matrix），描述输入如何影响状态；
    - $C \in \mathbb{R}^{p \times n}$ 是**输出矩阵**（output matrix），描述状态如何映射到输出；
    - $D \in \mathbb{R}^{p \times m}$ 是**直接传递矩阵**（feedthrough matrix）。

    四元组 $(A, B, C, D)$ 称为系统的**状态空间实现**（state-space realization）。

### 离散时间模型

!!! definition "定义 66A.2 (离散时间线性时不变系统)"
    **离散时间 LTI 系统**为
    $$
    \begin{aligned}
    x_{k+1} &= A x_k + B u_k \\
    y_k &= C x_k + D u_k
    \end{aligned}
    $$

    解为 $x_k = A^k x_0 + \sum_{j=0}^{k-1} A^{k-1-j} B u_j$。

    连续系统通过零阶保持器采样得到离散系统：$A_d = e^{A T_s}$，$B_d = \left(\int_0^{T_s} e^{A\tau}\,d\tau\right) B$，其中 $T_s$ 是采样周期。

!!! example "例 66A.1"
    **弹簧-质量-阻尼系统**：$m\ddot{q} + c\dot{q} + kq = f(t)$。

    取状态 $x = (q, \dot{q})^T$，输入 $u = f$，输出 $y = q$：
    $$
    A = \begin{pmatrix} 0 & 1 \\ -k/m & -c/m \end{pmatrix}, \;
    B = \begin{pmatrix} 0 \\ 1/m \end{pmatrix}, \;
    C = \begin{pmatrix} 1 & 0 \end{pmatrix}, \;
    D = 0
    $$

    系统矩阵 $A$ 的特征值 $\lambda = \frac{-c \pm \sqrt{c^2 - 4mk}}{2m}$ 决定自由振荡行为。

---

## 66A.2 矩阵指数解

<div class="context-flow" markdown>

**核心问题**：给定初始条件和输入，如何求状态方程的精确解？

</div>

!!! theorem "定理 66A.1 (状态方程的解)"
    给定初始条件 $x(0) = x_0$，系统 $\dot{x} = Ax + Bu$ 的解为
    $$
    x(t) = e^{At} x_0 + \int_0^t e^{A(t-\tau)} B u(\tau) \, d\tau
    $$

    相应的输出为
    $$
    y(t) = C e^{At} x_0 + \int_0^t C e^{A(t-\tau)} B u(\tau) \, d\tau + D u(t)
    $$

    第一项 $e^{At}x_0$ 称为**零输入响应**（自由响应），第二项为**零状态响应**（强迫响应）。

??? proof "证明"
    定义变量代换 $z(t) = e^{-At} x(t)$。对 $z(t)$ 求导，利用矩阵指数的求导法则：
    $$
    \dot{z}(t) = -A e^{-At} x(t) + e^{-At} \dot{x}(t)
    $$

    将 $\dot{x}(t) = Ax(t) + Bu(t)$ 代入：
    $$
    \dot{z}(t) = -A e^{-At} x(t) + e^{-At}(Ax(t) + Bu(t)) = e^{-At} Bu(t)
    $$

    两端从 $0$ 到 $t$ 积分：
    $$
    z(t) - z(0) = \int_0^t e^{-A\tau} Bu(\tau) \, d\tau
    $$

    由 $z(0) = e^{0}x(0) = x_0$，得
    $$
    z(t) = x_0 + \int_0^t e^{-A\tau} Bu(\tau) \, d\tau
    $$

    两端左乘 $e^{At}$：
    $$
    x(t) = e^{At}z(t) = e^{At} x_0 + \int_0^t e^{A(t-\tau)} Bu(\tau) \, d\tau
    $$

    输出公式由 $y(t) = Cx(t) + Du(t)$ 直接代入得到。$\blacksquare$

### 传递函数

!!! definition "定义 66A.3 (传递函数)"
    系统 $(A, B, C, D)$ 的**传递函数**（transfer function）是 $p \times m$ 有理矩阵
    $$
    G(s) = C(sI - A)^{-1}B + D
    $$

    传递函数将输入的 Laplace 变换 $U(s)$ 映射到输出的 Laplace 变换 $Y(s)$：$Y(s) = G(s) U(s)$（假设零初始条件）。

    $G(s)$ 的极点是 $(sI - A)^{-1}$ 的极点的子集，即 $A$ 的特征值的子集。只有能控且能观的模态才出现为传递函数的极点。

!!! example "例 66A.2"
    对例 66A.1 中的弹簧-质量-阻尼系统：
    $$
    G(s) = C(sI - A)^{-1}B + D = \frac{1/m}{s^2 + (c/m)s + k/m} = \frac{1}{ms^2 + cs + k}
    $$

    传递函数的极点与 $A$ 的特征值一致（因为系统能控且能观）。

---

## 66A.3 能控性

<div class="context-flow" markdown>

**核心问题**：什么条件下，系统的状态可以通过控制输入从任意初态转移到任意终态？

</div>

### 能控性的定义

!!! definition "定义 66A.4 (能控性)"
    系统 $(A, B)$ 称为**能控的**（controllable），若对任意初始状态 $x_0 \in \mathbb{R}^n$ 和目标状态 $x_f \in \mathbb{R}^n$，存在有限时间 $T > 0$ 和控制输入 $u(\cdot)$，使得系统从 $x(0) = x_0$ 到达 $x(T) = x_f$。

    等价地，$(A, B)$ 能控当且仅当从原点出发可以在有限时间内到达 $\mathbb{R}^n$ 中的任意状态（由线性性，这等价于从任意初态到任意终态）。

### Kalman 能控性秩条件

!!! definition "定义 66A.5 (能控性矩阵)"
    **能控性矩阵**（controllability matrix）定义为
    $$
    \mathcal{C} = \begin{pmatrix} B & AB & A^2B & \cdots & A^{n-1}B \end{pmatrix} \in \mathbb{R}^{n \times nm}
    $$

!!! theorem "定理 66A.2 (Kalman 能控性秩条件)"
    系统 $(A, B)$ 能控，当且仅当能控性矩阵满秩：
    $$
    \operatorname{rank}(\mathcal{C}) = \operatorname{rank}\begin{pmatrix} B & AB & A^2B & \cdots & A^{n-1}B \end{pmatrix} = n
    $$

??? proof "证明"
    **充分性**：设 $\operatorname{rank}(\mathcal{C}) = n$。对于连续时间系统，从原点出发在时间 $[0, T]$ 内可达的状态集为
    $$
    \mathcal{R}_T = \left\{\int_0^T e^{A(T-\tau)} Bu(\tau) \, d\tau : u \in L^2([0,T]; \mathbb{R}^m)\right\}
    $$

    关键步骤是证明 $\mathcal{R}_T = \operatorname{Im}(\mathcal{C})$。注意到矩阵指数可以展开为
    $$
    e^{A(T-\tau)} = \sum_{k=0}^{\infty} \frac{(T-\tau)^k}{k!} A^k
    $$

    因此 $e^{A(T-\tau)}B$ 的每一列都是 $B, AB, A^2B, \ldots$ 的列的线性组合。但由 Cayley-Hamilton 定理，$A^n$ 可以用 $I, A, \ldots, A^{n-1}$ 线性表示，因此 $A^k B$（$k \geq n$）可以用 $B, AB, \ldots, A^{n-1}B$ 的列线性表示。

    从而 $e^{A(T-\tau)}B$ 的列空间包含在 $\operatorname{Im}(\mathcal{C})$ 中，即 $\mathcal{R}_T \subseteq \operatorname{Im}(\mathcal{C})$。

    反之，可以证明 $\operatorname{Im}(\mathcal{C}) \subseteq \mathcal{R}_T$：对任意 $v \in \operatorname{Im}(\mathcal{C})$，选取输入 $u(\tau) = B^T e^{A^T(T-\tau)} W_c(T)^{-1} v$（其中 $W_c(T)$ 是能控性 Gramian，见后文），可以验证 $x(T) = v$。

    因此 $\operatorname{rank}(\mathcal{C}) = n$ 蕴含 $\mathcal{R}_T = \mathbb{R}^n$，即系统能控。

    **必要性**：若 $\operatorname{rank}(\mathcal{C}) < n$，则存在非零向量 $v \in \mathbb{R}^n$ 使得 $v^T A^k B = 0$ 对所有 $k = 0, 1, \ldots, n-1$。由 Cayley-Hamilton 定理推得 $v^T A^k B = 0$ 对所有 $k \geq 0$。从而
    $$
    v^T e^{A(T-\tau)} B = \sum_{k=0}^{\infty} \frac{(T-\tau)^k}{k!} v^T A^k B = 0
    $$

    因此对任意输入 $u(\cdot)$：
    $$
    v^T x(T) = v^T e^{AT} x_0 + \int_0^T v^T e^{A(T-\tau)} Bu(\tau) \, d\tau = v^T e^{AT} x_0
    $$

    这与输入 $u$ 无关。从 $x_0 = 0$ 出发无法到达满足 $v^T x_f \neq 0$ 的状态 $x_f$。$\blacksquare$

### PBH 能控性判据

!!! theorem "定理 66A.3 (Popov-Belevitch-Hautus (PBH) 判据)"
    系统 $(A, B)$ 能控，当且仅当对 $A$ 的每个特征值 $\lambda$，
    $$
    \operatorname{rank}\begin{pmatrix} A - \lambda I & B \end{pmatrix} = n
    $$

    等价地，不存在 $A^T$ 的左特征向量 $v$（即 $A^T v = \lambda v$，$v \neq 0$）使得 $B^T v = 0$。

??? proof "证明"
    **$\Rightarrow$（能控蕴含 PBH 条件）**：用反证法。若存在 $\lambda \in \mathbb{C}$ 使得 $\operatorname{rank}(A - \lambda I, B) < n$，则存在非零 $v$ 满足 $v^*(A - \lambda I) = 0$ 且 $v^*B = 0$（这里 $v^*$ 表示共轭转置）。即 $v^*A = \lambda v^*$，从而
    $$
    v^* A^k B = \lambda^k v^* B = 0, \quad \forall k \geq 0
    $$

    因此 $v^* \mathcal{C} = 0$，即 $\operatorname{rank}(\mathcal{C}) < n$，与能控性矛盾。

    **$\Leftarrow$（PBH 条件蕴含能控）**：用反证法。若 $\operatorname{rank}(\mathcal{C}) < n$，令 $V = \operatorname{Im}(\mathcal{C})$。注意 $V$ 是 $A$-不变子空间（因为 $A \cdot A^k B = A^{k+1}B$，当 $k \leq n-2$ 时直接在 $V$ 中，$k = n-1$ 时由 Cayley-Hamilton 也在 $V$ 中）。

    考虑 $A$ 限制在 $V^\perp$ 上的作用。由于 $V$ 是 $A$-不变的，$V^\perp$ 是 $A^T$-不变的。令 $\mu$ 为 $A^T|_{V^\perp}$ 的一个特征值，$w \in V^\perp$，$w \neq 0$，$A^T w = \mu w$。则 $w \perp V \supseteq \operatorname{Im}(B)$，即 $B^T w = 0$。这意味着
    $$
    \operatorname{rank}\begin{pmatrix} A - \mu I & B \end{pmatrix} < n
    $$

    与 PBH 条件矛盾。$\blacksquare$

!!! example "例 66A.3"
    考虑系统 $A = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$。

    **Kalman 秩条件**：$\mathcal{C} = \begin{pmatrix} B & AB \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & 2 \end{pmatrix}$，$\det(\mathcal{C}) = -1 \neq 0$。系统能控。

    **PBH 检验**：$\lambda_1 = 1$，$\operatorname{rank}\begin{pmatrix} 0 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix} = 2 = n$。$\lambda_2 = 2$，$\operatorname{rank}\begin{pmatrix} -1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = 2 = n$。均满秩，能控。

!!! example "例 66A.4"
    考虑对角系统 $A = \begin{pmatrix} \lambda_1 & 0 \\ 0 & \lambda_2 \end{pmatrix}$（$\lambda_1 \neq \lambda_2$），$B = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}$。

    $\mathcal{C} = \begin{pmatrix} b_1 & \lambda_1 b_1 \\ b_2 & \lambda_2 b_2 \end{pmatrix}$，$\det(\mathcal{C}) = b_1 b_2(\lambda_2 - \lambda_1)$。

    因为 $\lambda_1 \neq \lambda_2$，系统能控当且仅当 $b_1 \neq 0$ 且 $b_2 \neq 0$。

    物理直觉：每个模态必须受到输入的直接激励，否则输入无法影响该模态的状态。

### 能控标准型

!!! definition "定义 66A.6 (能控标准型)"
    若单输入系统 $(A, b)$（$b \in \mathbb{R}^n$）能控，则存在坐标变换将其化为**能控标准型**（controllable canonical form）：
    $$
    \bar{A} = \begin{pmatrix} 0 & 1 & 0 & \cdots & 0 \\ 0 & 0 & 1 & \cdots & 0 \\ \vdots & & & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & 1 \\ -a_0 & -a_1 & -a_2 & \cdots & -a_{n-1} \end{pmatrix}, \quad
    \bar{b} = \begin{pmatrix} 0 \\ 0 \\ \vdots \\ 0 \\ 1 \end{pmatrix}
    $$
    其中 $\chi_A(\lambda) = \lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_1\lambda + a_0$ 是 $A$ 的特征多项式。

    变换矩阵由能控性矩阵和特征多项式的系数决定。

---

## 66A.4 能观性

<div class="context-flow" markdown>

**核心问题**：从输出测量能否唯一确定系统的初始状态？

</div>

### 能观性的定义

!!! definition "定义 66A.7 (能观性)"
    系统 $(A, C)$ 称为**能观的**（observable），若初始状态 $x_0$ 可以从输出 $y(t)$（$t \in [0, T]$）和输入 $u(t)$ 中唯一确定。

    等价地，$(A, C)$ 能观当且仅当 $Ce^{At}x_0 = 0$（$\forall t \geq 0$）蕴含 $x_0 = 0$。

### Kalman 能观性秩条件

!!! definition "定义 66A.8 (能观性矩阵)"
    **能观性矩阵**（observability matrix）定义为
    $$
    \mathcal{O} = \begin{pmatrix} C \\ CA \\ CA^2 \\ \vdots \\ CA^{n-1} \end{pmatrix} \in \mathbb{R}^{np \times n}
    $$

!!! theorem "定理 66A.4 (Kalman 能观性秩条件)"
    系统 $(A, C)$ 能观，当且仅当
    $$
    \operatorname{rank}(\mathcal{O}) = \operatorname{rank}\begin{pmatrix} C \\ CA \\ \vdots \\ CA^{n-1} \end{pmatrix} = n
    $$

??? proof "证明"
    对于零输入响应（$u = 0$），$y(t) = Ce^{At}x_0$。

    **充分性**：若 $\operatorname{rank}(\mathcal{O}) = n$。对 $y(t) = Ce^{At}x_0$ 逐次求导并在 $t = 0$ 处取值：
    $$
    y(0) = Cx_0, \quad \dot{y}(0) = CAx_0, \quad \ldots, \quad y^{(n-1)}(0) = CA^{n-1}x_0
    $$

    写成矩阵形式：
    $$
    \begin{pmatrix} y(0) \\ \dot{y}(0) \\ \vdots \\ y^{(n-1)}(0) \end{pmatrix} = \mathcal{O} \cdot x_0
    $$

    由 $\operatorname{rank}(\mathcal{O}) = n$，$x_0$ 被唯一确定（$x_0 = \mathcal{O}^\dagger \cdot (y(0); \dot{y}(0); \ldots; y^{(n-1)}(0))$）。

    **必要性**：若 $\operatorname{rank}(\mathcal{O}) < n$，存在 $v \neq 0$ 使得 $\mathcal{O}v = 0$，即 $CA^k v = 0$ 对 $k = 0, 1, \ldots, n-1$。由 Cayley-Hamilton 定理，$CA^k v = 0$ 对所有 $k \geq 0$。因此
    $$
    Ce^{At}v = \sum_{k=0}^{\infty} \frac{t^k}{k!} CA^k v = 0, \quad \forall t \geq 0
    $$

    这意味着 $y(t)$ 对初始状态 $x_0 = 0$ 和 $x_0 = v$ 完全相同，无法区分。$\blacksquare$

### PBH 能观性判据

!!! theorem "定理 66A.5 (PBH 能观性判据)"
    $(A, C)$ 能观，当且仅当对 $A$ 的每个特征值 $\lambda$，
    $$
    \operatorname{rank}\begin{pmatrix} A - \lambda I \\ C \end{pmatrix} = n
    $$

    等价地，不存在 $A$ 的特征向量 $v$（$Av = \lambda v$，$v \neq 0$）使得 $Cv = 0$。

??? proof "证明"
    这是能控性 PBH 判据（定理 66A.3）通过对偶性的直接推论。$(A, C)$ 能观等价于 $(A^T, C^T)$ 能控（见定理 66A.7），而 $(A^T, C^T)$ 能控等价于对 $A^T$ 的每个特征值 $\bar{\lambda}$，$\operatorname{rank}(A^T - \bar{\lambda}I, C^T) = n$。取转置并注意 $A^T$ 和 $A$ 的特征值互为共轭（实矩阵时实际相同），即得所需条件。$\blacksquare$

!!! example "例 66A.5"
    系统 $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}$，$C = \begin{pmatrix} 1 & 0 \end{pmatrix}$。

    $\mathcal{O} = \begin{pmatrix} C \\ CA \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I_2$，满秩。系统能观。

    直观理解：$y = x_1$（直接测量第一个状态），$\dot{y} = \dot{x}_1 = x_2$（通过导数得到第二个状态）。

!!! example "例 66A.6"
    不能观的例子：$A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$，$C = \begin{pmatrix} 1 & 0 \end{pmatrix}$。

    $\mathcal{O} = \begin{pmatrix} 1 & 0 \\ 1 & 0 \end{pmatrix}$，$\operatorname{rank}(\mathcal{O}) = 1 < 2$。不能观。

    PBH 检验：$\lambda = 2$ 时，$\operatorname{rank}\begin{pmatrix} -1 & 0 \\ 0 & 0 \\ 1 & 0 \end{pmatrix} = 1 < 2$。特征向量 $v = (0, 1)^T$ 满足 $Cv = 0$——第二个状态完全不可见。

---

## 66A.5 能控能观对偶性

<div class="context-flow" markdown>

**核心问题**：能控性和能观性之间有什么对称关系？

</div>

!!! theorem "定理 66A.6 (对偶原理)"
    $(A, B)$ 能控当且仅当 $(A^T, B^T)$ 能观。

    等价地，$(A, C)$ 能观当且仅当 $(A^T, C^T)$ 能控。

??? proof "证明"
    $(A, B)$ 的能控性矩阵为 $\mathcal{C} = (B, AB, \ldots, A^{n-1}B)$。

    $(A^T, B^T)$ 的能观性矩阵为
    $$
    \mathcal{O}' = \begin{pmatrix} B^T \\ B^T A^T \\ \vdots \\ B^T (A^T)^{n-1} \end{pmatrix} = \begin{pmatrix} B^T \\ (AB)^T \\ \vdots \\ (A^{n-1}B)^T \end{pmatrix} = \mathcal{C}^T
    $$

    因此 $\operatorname{rank}(\mathcal{C}) = \operatorname{rank}(\mathcal{C}^T) = \operatorname{rank}(\mathcal{O}')$。

    同理可证 $(A, C)$ 能观等价于 $(A^T, C^T)$ 能控。$\blacksquare$

!!! note "注"
    对偶原理将能控性问题和能观性问题统一起来：任何关于能控性的定理，通过转置可以立即得到关于能观性的对应定理。这种对偶性在控制系统设计中也有深刻的工程含义——控制器设计（极点配置 $K$ 使 $A - BK$ 稳定）和观测器设计（选择 $L$ 使 $A - LC$ 稳定）本质上是对偶问题。

---

## 66A.6 能控性与能观性 Gramian

<div class="context-flow" markdown>

**核心问题**：能控性和能观性能否定量地刻画——"多大程度上"能控或能观？

</div>

### 能控性 Gramian

!!! definition "定义 66A.9 (能控性 Gramian)"
    系统 $(A, B)$ 的**有限时间能控性 Gramian** 定义为
    $$
    W_c(T) = \int_0^T e^{A\tau} B B^T e^{A^T\tau} \, d\tau
    $$

    若 $A$ 是稳定矩阵（所有特征值具有负实部），**无限时间能控性 Gramian** 定义为
    $$
    W_c = \int_0^{\infty} e^{A\tau} B B^T e^{A^T\tau} \, d\tau
    $$

!!! theorem "定理 66A.7 (能控性 Gramian 的性质)"
    1. $W_c(T) \succeq 0$（半正定），且 $W_c(T) \succ 0$ 当且仅当 $(A, B)$ 能控。
    2. 若 $(A, B)$ 能控，从原点到状态 $x_f$ 在时间 $[0, T]$ 内所需的**最小能量**为
    $$
    E_{\min} = \min_{u} \int_0^T \|u(\tau)\|^2 d\tau = x_f^T W_c(T)^{-1} x_f
    $$
    最优输入为 $u^*(\tau) = B^T e^{A^T(T-\tau)} W_c(T)^{-1} x_f$。
    3. $W_c(T)$ 的特征值越大，对应方向上的控制越"容易"（需要的能量越小）。

??? proof "证明"
    **性质 1 的证明**：$W_c(T) \succeq 0$ 是显然的（对任意 $v$，$v^T W_c(T) v = \int_0^T \|B^T e^{A^T \tau} v\|^2 d\tau \geq 0$）。

    若 $W_c(T) v = 0$，则 $\int_0^T \|B^T e^{A^T\tau} v\|^2 d\tau = 0$，由连续性得 $B^T e^{A^T\tau} v = 0$ 对所有 $\tau \in [0, T]$。在 $\tau = 0$ 处逐次求导得 $B^T (A^T)^k v = 0$ 对所有 $k \geq 0$，即 $v^T A^k B = 0$，即 $v^T \mathcal{C} = 0$。因此 $W_c(T) \succ 0$ 等价于 $\operatorname{rank}(\mathcal{C}) = n$。

    **性质 2 的证明**：这是一个约束优化问题。令 $x_f = \int_0^T e^{A(T-\tau)} Bu(\tau) d\tau$，利用 Lagrange 乘子法或 Cauchy-Schwarz 不等式可得最小能量解。代入 $u^*(\tau) = B^T e^{A^T(T-\tau)} \lambda$ 并由终态条件 $x_f = W_c(T) \lambda$ 得 $\lambda = W_c(T)^{-1} x_f$。$\blacksquare$

### 能观性 Gramian

!!! definition "定义 66A.10 (能观性 Gramian)"
    系统 $(A, C)$ 的**有限时间能观性 Gramian** 定义为
    $$
    W_o(T) = \int_0^T e^{A^T\tau} C^T C e^{A\tau} \, d\tau
    $$

    若 $A$ 稳定，**无限时间能观性 Gramian** 为
    $$
    W_o = \int_0^{\infty} e^{A^T\tau} C^T C e^{A\tau} \, d\tau
    $$

!!! theorem "定理 66A.8 (能观性 Gramian 的性质)"
    1. $W_o(T) \succ 0$ 当且仅当 $(A, C)$ 能观。
    2. $W_o$ 量化了从输出中提取初始状态信息的"容易程度"：
    $$
    \int_0^T \|y(\tau)\|^2 d\tau = x_0^T W_o(T) x_0 \quad (\text{零输入响应})
    $$
    $W_o$ 特征值大的方向对应"容易观测"的状态方向。

### Lyapunov 方程计算 Gramian

!!! theorem "定理 66A.9 (Gramian 的 Lyapunov 方程)"
    若 $A$ 是稳定矩阵，则无限时间 Gramian 分别满足以下**连续 Lyapunov 方程**：
    $$
    A W_c + W_c A^T + BB^T = 0
    $$
    $$
    A^T W_o + W_o A + C^T C = 0
    $$

??? proof "证明"
    对 $W_c = \int_0^{\infty} e^{A\tau} BB^T e^{A^T\tau} d\tau$，两端左乘 $A$ 并加上右乘 $A^T$：
    $$
    AW_c + W_c A^T = \int_0^{\infty} \left(Ae^{A\tau} BB^T e^{A^T\tau} + e^{A\tau} BB^T e^{A^T\tau} A^T\right) d\tau
    $$
    $$
    = \int_0^{\infty} \frac{d}{d\tau}\left(e^{A\tau} BB^T e^{A^T\tau}\right) d\tau = \left[e^{A\tau} BB^T e^{A^T\tau}\right]_0^{\infty}
    $$

    由于 $A$ 稳定，$e^{A\tau} \to 0$（$\tau \to \infty$），故上式等于 $0 - BB^T = -BB^T$。

    因此 $AW_c + W_c A^T + BB^T = 0$。能观性 Gramian 的 Lyapunov 方程类似推导。$\blacksquare$

!!! note "注"
    Lyapunov 方程是线性矩阵方程（Ch20），可以高效求解（例如 Bartels-Stewart 算法，复杂度 $O(n^3)$）。在实际计算中，通常通过求解 Lyapunov 方程来获得 Gramian，而非直接计算积分。

---

## 66A.7 Kalman 典则分解

<div class="context-flow" markdown>

**核心问题**：如何将一般系统分解为能控/不能控和能观/不能观的子系统？

</div>

!!! theorem "定理 66A.10 (Kalman 分解)"
    对任意系统 $(A, B, C, D)$，存在（非奇异）坐标变换 $T$ 使得系统在新坐标下具有四部分分块形式：
    $$
    \bar{A} = T^{-1}AT = \begin{pmatrix} A_{c\bar{o}} & 0 & A_{13} & 0 \\ A_{21} & A_{co} & A_{23} & A_{24} \\ 0 & 0 & A_{\bar{c}\bar{o}} & 0 \\ 0 & 0 & A_{43} & A_{\bar{c}o} \end{pmatrix}
    $$
    $$
    \bar{B} = T^{-1}B = \begin{pmatrix} B_{c\bar{o}} \\ B_{co} \\ 0 \\ 0 \end{pmatrix}, \quad
    \bar{C} = CT = \begin{pmatrix} 0 & C_{co} & 0 & C_{\bar{c}o} \end{pmatrix}
    $$

    其中下标含义为：

    - $co$：能控且能观（controllable and observable）
    - $\bar{c}o$：不能控但能观
    - $c\bar{o}$：能控但不能观
    - $\bar{c}\bar{o}$：不能控且不能观

??? proof "证明"
    **构造过程**：

    **第一步**：令 $V_c = \operatorname{Im}(\mathcal{C})$ 为能控子空间（$A$-不变），$V_o = \ker(\mathcal{O})^\perp$ 为能观子空间的补。更准确地，不能观子空间为 $V_{\bar{o}} = \ker(\mathcal{O})$，它也是 $A$-不变的。

    **第二步**：状态空间分解为四个子空间：
    $$
    \mathbb{R}^n = (V_c \cap V_{\bar{o}}) \oplus (V_c \cap V_{\bar{o}}^\perp) \oplus (V_c^\perp \cap V_{\bar{o}}) \oplus (V_c^\perp \cap V_{\bar{o}}^\perp)
    $$
    （严格来说，这里的"正交补"需要用 $A$-不变子空间的理论更仔细地处理，但核心思想如此。）

    **第三步**：选取每个子空间的基底构成变换矩阵 $T$ 的列。在新基底下，$A$ 的不变性保证了分块上三角结构，$B$ 和 $C$ 的结构由能控/能观子空间的定义决定。

    关键在于：$B$ 的像在 $V_c$ 中（能控子空间），所以 $\bar{B}$ 在不能控的分量上为零；$C$ 在 $V_{\bar{o}} = \ker(\mathcal{O})$ 上为零，所以 $\bar{C}$ 在不能观的分量上为零。$\blacksquare$

!!! note "注"
    Kalman 分解揭示了一个深刻的事实：传递函数 $G(s) = C(sI-A)^{-1}B + D$ 只取决于**能控且能观**的子系统 $(A_{co}, B_{co}, C_{co}, D)$。不能控或不能观的模态在传递函数中被"对消"了——这就是极零对消（pole-zero cancellation）现象的结构性解释。

---

## 66A.8 最小实现

<div class="context-flow" markdown>

**核心问题**：给定传递函数，维数最小的状态空间实现是什么？

</div>

!!! definition "定义 66A.11 (最小实现)"
    给定传递函数 $G(s)$，其**最小实现**（minimal realization）$(A, B, C, D)$ 是维数最小的状态空间实现，即状态维度 $n$ 在所有满足 $G(s) = C(sI-A)^{-1}B + D$ 的实现中最小。

!!! theorem "定理 66A.11 (最小实现的刻画)"
    实现 $(A, B, C, D)$ 是最小的，当且仅当 $(A, B)$ 能控且 $(A, C)$ 能观。

    最小实现的状态维度等于 $G(s)$ 的 **McMillan 度**——传递函数的 Smith-McMillan 形式中分母多项式的次数之和。

    所有最小实现之间通过相似变换联系：若 $(A_1, B_1, C_1, D)$ 和 $(A_2, B_2, C_2, D)$ 都是 $G(s)$ 的最小实现，则存在非奇异矩阵 $T$ 使得 $A_2 = T^{-1}A_1T$，$B_2 = T^{-1}B_1$，$C_2 = C_1 T$。

### Ho-Kalman 算法

!!! definition "定义 66A.12 (Markov 参数与 Hankel 矩阵)"
    系统 $(A, B, C, D)$ 的 **Markov 参数**为
    $$
    h_k = CA^{k-1}B, \quad k = 1, 2, 3, \ldots
    $$

    它们是系统脉冲响应的矩阵系数（$h_k$ 是 $G(s)$ 的 Laurent 展开系数）。

    **Hankel 矩阵** $H_{r,s}$ 由 Markov 参数构成：
    $$
    H_{r,s} = \begin{pmatrix} h_1 & h_2 & \cdots & h_s \\ h_2 & h_3 & \cdots & h_{s+1} \\ \vdots & \vdots & \ddots & \vdots \\ h_r & h_{r+1} & \cdots & h_{r+s-1} \end{pmatrix} = \mathcal{O}_r \cdot \mathcal{C}_s
    $$
    其中 $\mathcal{O}_r$ 是前 $r$ 块行的能观性矩阵，$\mathcal{C}_s$ 是前 $s$ 块列的能控性矩阵。

!!! theorem "定理 66A.12 (Ho-Kalman 算法)"
    给定 Markov 参数 $\{h_k\}_{k=1}^{2N}$，最小实现的维度为 $n = \operatorname{rank}(H_{N,N})$。

    **算法步骤**：

    1. 构造 Hankel 矩阵 $H = H_{N,N}$；
    2. 对 $H$ 做满秩分解或 SVD：$H = U \Sigma V^T$，保留前 $n$ 个非零奇异值，得 $H \approx U_n \Sigma_n V_n^T$；
    3. 令 $\mathcal{O}_n = U_n \Sigma_n^{1/2}$，$\mathcal{C}_n = \Sigma_n^{1/2} V_n^T$；
    4. 从中提取最小实现：
    $$
    C = \mathcal{O}_n \text{ 的前 } p \text{ 行}, \quad B = \mathcal{C}_n \text{ 的前 } m \text{ 列}
    $$
    $$
    A = \mathcal{O}_n^\dagger H' (\mathcal{C}_n)^\dagger
    $$
    其中 $H'$ 是将 $H$ 中的 $h_k$ 替换为 $h_{k+1}$ 得到的移位 Hankel 矩阵。

!!! example "例 66A.7"
    设 SISO 系统的脉冲响应为 $h_1 = 1, h_2 = -1, h_3 = 1, h_4 = -1, \ldots$

    Hankel 矩阵：$H = \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$，$\operatorname{rank}(H) = 1$。

    系统的最小实现为一阶：$A = -1$，$B = 1$，$C = 1$，$D = 0$。

    验证：$CA^{k-1}B = (-1)^{k-1}$，即 $h_k = (-1)^{k-1}$。传递函数 $G(s) = \frac{1}{s+1}$。

---

## 66A.9 平衡实现

<div class="context-flow" markdown>

**核心问题**：是否存在一种"最佳"的坐标系使得能控性和能观性同时得到最清晰的刻画？

</div>

### Hankel 奇异值

!!! definition "定义 66A.13 (Hankel 奇异值)"
    稳定最小实现 $(A, B, C)$ 的 **Hankel 奇异值** $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_n > 0$ 定义为
    $$
    \sigma_i = \sqrt{\lambda_i(W_c W_o)}
    $$

    其中 $W_c, W_o$ 分别是能控性和能观性 Gramian，$\lambda_i(\cdot)$ 表示特征值（按降序排列）。

    Hankel 奇异值是系统的**不变量**——不依赖于状态空间坐标的选择。它们等于系统 Hankel 算子的奇异值。

!!! theorem "定理 66A.13 (Hankel 奇异值的不变性)"
    Hankel 奇异值 $\sigma_i = \sqrt{\lambda_i(W_c W_o)}$ 在相似变换 $x \mapsto Tx$ 下不变。

??? proof "证明"
    设坐标变换 $\bar{x} = T^{-1}x$，则新坐标下的 Gramian 为
    $$
    \bar{W}_c = T^{-1} W_c T^{-T}, \quad \bar{W}_o = T^T W_o T
    $$

    从而 $\bar{W}_c \bar{W}_o = T^{-1} W_c T^{-T} \cdot T^T W_o T = T^{-1} W_c W_o T$。

    因此 $\bar{W}_c \bar{W}_o$ 与 $W_c W_o$ 相似，特征值相同。$\blacksquare$

### 平衡实现的定义与存在性

!!! definition "定义 66A.14 (平衡实现)"
    稳定最小实现 $(A, B, C)$ 称为**平衡实现**（balanced realization），若在该坐标下能控性 Gramian 和能观性 Gramian 相等且为对角矩阵：
    $$
    W_c = W_o = \Sigma = \operatorname{diag}(\sigma_1, \sigma_2, \ldots, \sigma_n)
    $$

    其中 $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_n > 0$ 恰好是 Hankel 奇异值。

!!! theorem "定理 66A.14 (平衡实现的存在性)"
    任何稳定的最小实现都可以通过坐标变换化为平衡实现。

??? proof "证明"
    **构造步骤**：

    **第一步**：求解 Lyapunov 方程得 $W_c$ 和 $W_o$。

    **第二步**：对 $W_c$ 做 Cholesky 分解 $W_c = L_c L_c^T$。

    **第三步**：对 $L_c^T W_o L_c$ 做特征值分解 $L_c^T W_o L_c = U \Lambda^2 U^T$（其中 $\Lambda = \operatorname{diag}(\sigma_1, \ldots, \sigma_n)$）。

    **第四步**：令 $T = L_c U \Lambda^{-1/2}$。

    **验证**：新坐标下
    $$
    \bar{W}_c = T^{-1} W_c T^{-T} = \Lambda^{1/2} U^T L_c^{-1} \cdot L_c L_c^T \cdot L_c^{-T} U \Lambda^{1/2}
    $$
    经仔细化简可得 $\bar{W}_c = \Lambda$。类似地 $\bar{W}_o = \Lambda$。$\blacksquare$

!!! note "注"
    平衡实现与奇异值分解（SVD）有深刻联系。在平衡实现中，Hankel 奇异值 $\sigma_i$ 同时量化了状态 $x_i$ 的"可控程度"和"可观程度"。$\sigma_i$ 大的状态同时容易控制和观测，$\sigma_i$ 小的状态在输入输出行为中贡献很小。

---

## 66A.10 平衡截断与模型降阶

<div class="context-flow" markdown>

**核心问题**：如何在保持输入-输出行为的前提下降低系统维度？

</div>

!!! definition "定义 66A.15 (平衡截断)"
    给定平衡实现 $(A, B, C)$，按 Hankel 奇异值将状态分为"重要"和"不重要"两组：
    $$
    \Sigma = \begin{pmatrix} \Sigma_1 & 0 \\ 0 & \Sigma_2 \end{pmatrix}, \quad
    A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}, \quad
    B = \begin{pmatrix} B_1 \\ B_2 \end{pmatrix}, \quad
    C = \begin{pmatrix} C_1 & C_2 \end{pmatrix}
    $$

    其中 $\Sigma_1 = \operatorname{diag}(\sigma_1, \ldots, \sigma_r)$ 包含较大的 Hankel 奇异值，$\Sigma_2$ 包含较小的。

    **平衡截断**（balanced truncation）保留前 $r$ 个状态得到降阶模型：
    $$
    \hat{A} = A_{11}, \quad \hat{B} = B_1, \quad \hat{C} = C_1
    $$

!!! theorem "定理 66A.15 (平衡截断的误差界)"
    设 $G(s)$ 是原系统的传递函数，$\hat{G}(s)$ 是 $r$ 阶平衡截断的传递函数。则
    $$
    \|G - \hat{G}\|_{H_\infty} \leq 2(\sigma_{r+1} + \sigma_{r+2} + \cdots + \sigma_n)
    $$

    其中 $\|G\|_{H_\infty} = \sup_\omega \bar{\sigma}(G(j\omega))$ 是 $H_\infty$ 范数（最大奇异值在频率上的上确界）。

    此外，截断后的降阶系统仍然是稳定的。

!!! note "注"
    平衡截断是目前最广泛使用的线性系统降阶方法之一。其优势包括：

    - **先验误差界**：截断前就可以通过 Hankel 奇异值估计降阶误差；
    - **保持稳定性**：截断后的系统自动稳定；
    - **与 SVD 的联系**：整个过程本质上是对系统 Hankel 算子做"矩阵 SVD"再截断。

!!! example "例 66A.8"
    考虑一个 100 阶系统，其 Hankel 奇异值为 $\sigma_1 = 10, \sigma_2 = 5, \sigma_3 = 0.1, \sigma_4 = 0.05, \ldots, \sigma_{100} = 10^{-8}$。

    截断到 2 阶，误差界为 $2(\sigma_3 + \cdots + \sigma_{100}) \approx 2 \times 0.15 = 0.3$。

    由于前两个 Hankel 奇异值远大于其余值，2 阶模型就能很好地近似 100 阶系统的输入-输出行为。这在大规模系统（如有限元模型）的降阶中极其实用。

---

## 本章小结

本章展示了线性代数在控制系统结构分析中的核心作用。主要结果包括：

1. **状态空间模型** $(A, B, C, D)$ 用四个矩阵完整描述线性时不变系统，解由矩阵指数 $e^{At}$ 给出。

2. **能控性**由 Kalman 秩条件（$\operatorname{rank}(\mathcal{C}) = n$）或 PBH 判据（$\operatorname{rank}(A - \lambda I, B) = n$）刻画，保证可以通过输入将状态转移到任意位置。

3. **能观性**由 $\operatorname{rank}(\mathcal{O}) = n$ 或 PBH 判据刻画，保证可以从输出重构状态。

4. **对偶原理**统一了能控性和能观性：$(A, B)$ 能控等价于 $(A^T, B^T)$ 能观。

5. **Gramian** 提供了能控性/能观性的定量刻画，满足 Lyapunov 方程，并给出最小控制能量和输出能量的表达式。

6. **Kalman 分解**将系统分为四个子系统，揭示了只有能控且能观的部分才出现在传递函数中。

7. **最小实现**等价于能控且能观的实现，**Ho-Kalman 算法**从 Markov 参数（传递函数）恢复最小实现。

8. **平衡实现**使 $W_c = W_o = \operatorname{diag}(\sigma_1, \ldots, \sigma_n)$，**平衡截断**通过丢弃小 Hankel 奇异值对应的状态实现模型降阶，误差有先验上界。

---

## 习题

!!! question "习题 66A.1"
    验证系统 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$ 的能控性。写出其能控标准型。

!!! question "习题 66A.2"
    设 $A = \begin{pmatrix} \lambda_1 & 0 \\ 0 & \lambda_2 \end{pmatrix}$（$\lambda_1 \neq \lambda_2$），$B = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}$。证明 $(A, B)$ 能控当且仅当 $b_1 \neq 0$ 且 $b_2 \neq 0$。

!!! question "习题 66A.3"
    对系统 $A = \begin{pmatrix} -1 & 0 \\ 0 & -2 \end{pmatrix}$，$C = (1, 1)$，判断能观性。同时用 PBH 判据验证。

!!! question "习题 66A.4"
    验证对偶原理：对例 66A.3 中的系统，构造对偶系统 $(A^T, B^T)$ 并验证其能观性。

!!! question "习题 66A.5"
    对系统 $A = \begin{pmatrix} -1 & 0 \\ 0 & -2 \end{pmatrix}$，$B = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$，求解连续 Lyapunov 方程 $AW_c + W_cA^T + BB^T = 0$ 得到能控性 Gramian $W_c$。

!!! question "习题 66A.6"
    证明 Hankel 奇异值在相似变换下不变，即若 $\bar{A} = T^{-1}AT$，$\bar{B} = T^{-1}B$，$\bar{C} = CT$，则 $\lambda_i(\bar{W}_c \bar{W}_o) = \lambda_i(W_c W_o)$。

!!! question "习题 66A.7"
    对于离散时间系统 $x_{k+1} = Ax_k + Bu_k$，$y_k = Cx_k$，写出离散时间能控性矩阵和 Kalman 秩条件。证明离散时间的能控性 Gramian $W_c = \sum_{k=0}^{\infty} A^k BB^T (A^T)^k$ 满足离散 Lyapunov 方程 $AW_c A^T - W_c + BB^T = 0$。

!!! question "习题 66A.8"
    给定 SISO 系统的前四个 Markov 参数 $h_1 = 2, h_2 = -1, h_3 = 1/2, h_4 = -1/4$。用 Ho-Kalman 算法求最小实现的阶数，并求出一个最小实现 $(A, B, C)$。

!!! question "习题 66A.9"
    设稳定系统有 Hankel 奇异值 $\sigma_1 = 5, \sigma_2 = 3, \sigma_3 = 0.01, \sigma_4 = 0.005$。若用平衡截断将系统降至 2 阶，估计 $H_\infty$ 范数意义下的逼近误差上界。

!!! question "习题 66A.10"
    考虑系统 $A = \begin{pmatrix} -1 & 1 \\ 0 & -2 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$，$C = \begin{pmatrix} 1 & 0 \end{pmatrix}$。

    (a) 判断系统的能控性和能观性。

    (b) 写出传递函数 $G(s)$。

    (c) 判断该实现是否为最小实现。

!!! question "习题 66A.11"
    证明：若系统 $(A, B, C)$ 既能控又能观，且 $A$ 稳定，则 Hankel 矩阵 $H$ 的秩等于系统的阶数 $n$。

!!! question "习题 66A.12"
    设 $(A, B, C)$ 是平衡实现，证明 $A$ 的对角元素 $A_{ii}$ 满足 $A_{ii} \leq -\frac{\|B_i\|^2 + \|C_i\|^2}{4\sigma_i}$，其中 $B_i$ 是 $B$ 的第 $i$ 行，$C_i$ 是 $C$ 的第 $i$ 列。（提示：利用 Lyapunov 方程的对角元素。）
