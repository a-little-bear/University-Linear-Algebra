# 第 68 章 线性代数在机器人学中的应用

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 矩阵指数(Ch13) · 矩阵群(Ch55) · 辛矩阵(Ch53)

**本章脉络**：刚体运动 → $SE(3)$ 与齐次变换矩阵 → 旋量(twist)与螺旋运动 → 指数坐标 → 正运动学(指数积公式) → Jacobi 矩阵与速度 → 逆运动学 → 动力学中的线性代数

**延伸**：现代机器人学完全建立在矩阵 Lie 群 $SE(3)$ 上；SLAM（同步定位与建图）、视觉惯性里程计和姿态估计都需要在 $SE(3)$ 或 $SO(3)$ 上进行优化和滤波

</div>

机器人学是矩阵 Lie 群理论最自然的工程应用。机器人的每一个关节运动都是一个刚体变换，而刚体变换的数学描述恰好是特殊欧几里得群 $SE(3)$——一个矩阵 Lie 群。机器人手臂的正运动学是 $SE(3)$ 中矩阵指数的乘积，Jacobi 矩阵建立了关节速度与末端速度之间的线性映射，逆运动学则需要在 $SE(3)$ 上求解非线性方程。

本章展示 $SE(3)$ 的矩阵理论如何为机器人学提供统一而优雅的数学框架。从刚体运动的齐次矩阵表示出发，经由旋量理论和指数映射，到达正/逆运动学和动力学的矩阵公式化。

---

## 68.1 刚体运动与 $SE(3)$

<div class="context-flow" markdown>

**核心问题**：如何用矩阵精确描述三维空间中的刚体运动？

</div>

### 旋转群 $SO(3)$

!!! definition "定义 68.1 (旋转群 $SO(3)$)"
    三维**特殊正交群**定义为
    $$
    SO(3) = \{R \in \mathbb{R}^{3 \times 3} : R^T R = I, \; \det(R) = 1\}
    $$

    $SO(3)$ 的元素表示三维空间中保持原点不动的旋转。它是一个 3 维 Lie 群。

!!! theorem "定理 68.1 ($SO(3)$ 的基本性质)"
    1. $SO(3)$ 在矩阵乘法下构成群：$R_1, R_2 \in SO(3) \Rightarrow R_1 R_2 \in SO(3)$；
    2. 逆元为转置：$R^{-1} = R^T$；
    3. $SO(3)$ 是紧致的、连通的；
    4. $SO(3)$ 的 Lie 代数为 $\mathfrak{so}(3) = \{\Omega \in \mathbb{R}^{3 \times 3} : \Omega^T = -\Omega\}$，即 $3 \times 3$ 反对称矩阵集合。

### 特殊欧几里得群 $SE(3)$

!!! definition "定义 68.2 (特殊欧几里得群 $SE(3)$)"
    三维**特殊欧几里得群**定义为
    $$
    SE(3) = \left\{T = \begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} \in \mathbb{R}^{4 \times 4} : R \in SO(3), \; \mathbf{p} \in \mathbb{R}^3\right\}
    $$

    $SE(3)$ 的元素称为**齐次变换矩阵**（homogeneous transformation matrix），表示刚体运动：旋转 $R$ 加平移 $\mathbf{p}$。

!!! theorem "定理 68.2 ($SE(3)$ 的群性质)"
    1. **乘法**：$T_1 T_2 = \begin{pmatrix} R_1 R_2 & R_1 \mathbf{p}_2 + \mathbf{p}_1 \\ \mathbf{0}^T & 1 \end{pmatrix} \in SE(3)$；

    2. **逆元**：$T^{-1} = \begin{pmatrix} R^T & -R^T \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix}$；

    3. **单位元**：$I_4$；

    4. $SE(3)$ 是 6 维 Lie 群（3 个旋转自由度 + 3 个平移自由度）。

??? proof "证明"
    **逆元公式**：设 $T = \begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix}$，验证
    $$
    T \cdot T^{-1} = \begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} \begin{pmatrix} R^T & -R^T\mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} = \begin{pmatrix} RR^T & -RR^T\mathbf{p} + \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} = \begin{pmatrix} I & \mathbf{0} \\ \mathbf{0}^T & 1 \end{pmatrix}
    $$

!!! example "例 68.1"
    机器人手臂的末端执行器相对于基座的位姿可以用一个 $SE(3)$ 元素描述：
    $$
    T_{0e} = \begin{pmatrix} R_{0e} & \mathbf{p}_{0e} \\ \mathbf{0}^T & 1 \end{pmatrix}
    $$
    其中 $R_{0e} \in SO(3)$ 描述末端执行器的**姿态**（orientation），$\mathbf{p}_{0e} \in \mathbb{R}^3$ 描述末端执行器的**位置**（position）。

### 坐标系之间的变换

!!! theorem "定理 68.3 (坐标系变换)"
    设 $T_{ab} \in SE(3)$ 表示坐标系 $\{b\}$ 相对于坐标系 $\{a\}$ 的位姿。若点 $\mathbf{q}$ 在坐标系 $\{b\}$ 中的齐次坐标为 $\tilde{\mathbf{q}}_b = (\mathbf{q}_b; 1)$，则它在坐标系 $\{a\}$ 中的坐标为
    $$
    \tilde{\mathbf{q}}_a = T_{ab} \, \tilde{\mathbf{q}}_b
    $$

    若有三个坐标系 $\{a\}, \{b\}, \{c\}$，则
    $$
    T_{ac} = T_{ab} \cdot T_{bc}
    $$
    （坐标系变换的链式法则）。

!!! example "例 68.2"
    设世界坐标系 $\{0\}$、机器人基座坐标系 $\{b\}$、末端执行器坐标系 $\{e\}$。目标点在末端执行器坐标系中的坐标为 $\mathbf{q}_e$。它在世界坐标系中的坐标为
    $$
    \tilde{\mathbf{q}}_0 = T_{0b} \cdot T_{be} \cdot \tilde{\mathbf{q}}_e
    $$

---

## 68.2 旋量与螺旋运动

<div class="context-flow" markdown>

**核心问题**：如何用 Lie 代数语言描述刚体的瞬时运动？

</div>

### $SE(3)$ 的 Lie 代数

!!! definition "定义 68.3 (Lie 代数 $\mathfrak{se}(3)$)"
    $SE(3)$ 的 **Lie 代数** $\mathfrak{se}(3)$ 由以下 $4 \times 4$ 矩阵组成：
    $$
    \mathfrak{se}(3) = \left\{[\xi] = \begin{pmatrix} [\omega]_\times & \mathbf{v} \\ \mathbf{0}^T & 0 \end{pmatrix} \in \mathbb{R}^{4 \times 4} : \omega \in \mathbb{R}^3, \; \mathbf{v} \in \mathbb{R}^3\right\}
    $$

    其中 $[\omega]_\times = \begin{pmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{pmatrix}$ 是 $\omega$ 的反对称矩阵。

    六维向量 $\xi = (\mathbf{v}; \omega) \in \mathbb{R}^6$ 称为**旋量**（twist）。

!!! note "注"
    旋量 $\xi = (\mathbf{v}; \omega)$ 中：

    - $\omega \in \mathbb{R}^3$ 是**角速度**分量；
    - $\mathbf{v} \in \mathbb{R}^3$ 是**线速度**分量。

    不同文献中 $\mathbf{v}$ 和 $\omega$ 的排列顺序可能不同。本书采用 Murray-Li-Sastry 的约定 $\xi = (\mathbf{v}; \omega)$。

### Chasles 定理

!!! theorem "定理 68.4 (Chasles 定理)"
    三维空间中的任意刚体运动都可以分解为绕某条轴旋转加上沿该轴的平移。这条轴称为**螺旋轴**（screw axis），运动称为**螺旋运动**（screw motion）。

    数学上，每个 $T \in SE(3)$（$T \neq I$）都可以写成绕螺旋轴 $\ell$ 旋转角 $\theta$ 并沿 $\ell$ 平移距离 $d$ 的形式。

!!! definition "定义 68.4 (螺旋轴)"
    **螺旋轴**（screw）$\mathcal{S}$ 由以下参数描述：

    - 轴方向 $\hat{\omega} \in \mathbb{R}^3$（$\|\hat{\omega}\| = 1$）；
    - 轴上一点 $\mathbf{q} \in \mathbb{R}^3$；
    - 螺距 $h$（旋转一弧度时沿轴平移的距离）。

    对应的旋量为
    $$
    \xi = \begin{pmatrix} -\hat{\omega} \times \mathbf{q} + h\hat{\omega} \\ \hat{\omega} \end{pmatrix}
    $$

!!! example "例 68.3"
    **纯旋转**（$h = 0$）：绕过点 $\mathbf{q} = (0, 0, 0)$ 的 $z$ 轴旋转。
    $$
    \xi = \begin{pmatrix} \mathbf{0} \\ \hat{e}_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 1 \end{pmatrix}
    $$

    **纯平移**（$\|\omega\| = 0$）：沿方向 $\hat{\mathbf{v}}$ 平移。
    $$
    \xi = \begin{pmatrix} \hat{\mathbf{v}} \\ \mathbf{0} \end{pmatrix}
    $$

---

## 68.3 指数映射与对数映射

<div class="context-flow" markdown>

**核心问题**：如何在 Lie 代数（线性空间）和 Lie 群（流形）之间建立联系？

</div>

### $SO(3)$ 的指数映射

!!! theorem "定理 68.5 (Rodrigues 公式——$SO(3)$ 的指数映射)"
    对 $\omega \in \mathbb{R}^3$，$\hat{\omega} = \omega/\|\omega\|$（$\omega \neq 0$），$\theta = \|\omega\|$：
    $$
    \exp([\omega]_\times) = I + \sin\theta \, [\hat{\omega}]_\times + (1 - \cos\theta) [\hat{\omega}]_\times^2
    $$

    这是从 $\mathfrak{so}(3)$ 到 $SO(3)$ 的指数映射。每个旋转矩阵 $R \in SO(3)$ 都可以写成某个反对称矩阵的指数。

### $SE(3)$ 的指数映射

!!! theorem "定理 68.6 ($SE(3)$ 的指数映射)"
    对旋量 $\xi = (\mathbf{v}; \omega) \in \mathbb{R}^6$，$\theta = \|\omega\|$（$\omega \neq 0$）：
    $$
    \exp([\xi]\theta) = \exp\left(\begin{pmatrix} [\omega]_\times \theta & \mathbf{v}\theta \\ \mathbf{0}^T & 0 \end{pmatrix}\right) = \begin{pmatrix} e^{[\omega]_\times \theta} & J(\omega, \theta) \mathbf{v}\theta \\ \mathbf{0}^T & 1 \end{pmatrix}
    $$

    其中
    $$
    J(\omega, \theta) = I + \frac{1 - \cos\theta}{\theta} [\hat{\omega}]_\times + \frac{\theta - \sin\theta}{\theta} [\hat{\omega}]_\times^2
    $$
    （有时也写成 $J = \frac{\sin\theta}{\theta}I + \frac{1-\cos\theta}{\theta}[\hat{\omega}]_\times + (1 - \frac{\sin\theta}{\theta})\hat{\omega}\hat{\omega}^T$）

    对纯平移情形（$\omega = 0$）：
    $$
    \exp\left(\begin{pmatrix} 0 & \mathbf{v}\theta \\ \mathbf{0}^T & 0 \end{pmatrix}\right) = \begin{pmatrix} I & \mathbf{v}\theta \\ \mathbf{0}^T & 1 \end{pmatrix}
    $$

??? proof "证明"
    利用 $4 \times 4$ 矩阵的幂级数 $\exp([\xi]) = \sum_{k=0}^{\infty} \frac{[\xi]^k}{k!}$。

    关键观察：$[\xi]^k$ 的左上角 $3 \times 3$ 块为 $([\omega]_\times)^k$，而右上角 $3 \times 1$ 块可以通过递推关系计算。

    利用 $([\hat{\omega}]_\times)^3 = -[\hat{\omega}]_\times$（周期性），将无穷级数归结为有限的三角函数表达式。

### 对数映射

!!! theorem "定理 68.7 ($SO(3)$ 的对数映射)"
    给定 $R \in SO(3)$，其对数 $[\omega]_\times = \log(R)$ 为：

    - 若 $R = I$：$\omega = \mathbf{0}$；
    - 若 $\operatorname{tr}(R) \neq -1$：$\theta = \arccos\frac{\operatorname{tr}(R) - 1}{2}$，$[\hat{\omega}]_\times = \frac{R - R^T}{2\sin\theta}$；
    - 若 $\operatorname{tr}(R) = -1$（$\theta = \pi$）：$\hat{\omega}$ 是 $R + I$ 的任意非零列的归一化。

!!! example "例 68.4"
    设 $R = \begin{pmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}$（绕 $z$ 轴旋转 $90°$）。

    $\operatorname{tr}(R) = 0 + 0 + 1 = 1$，$\theta = \arccos\frac{1-1}{2} = \arccos 0 = \frac{\pi}{2}$。

    $[\hat{\omega}]_\times = \frac{R - R^T}{2\sin(\pi/2)} = \frac{1}{2}\begin{pmatrix} 0 & -2 & 0 \\ 2 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}$

    读出 $\hat{\omega} = (0, 0, 1)^T$，符合预期。

---

## 68.4 正运动学

<div class="context-flow" markdown>

**核心问题**：给定关节角度，如何计算末端执行器的位姿？

</div>

### 指数积公式

!!! definition "定义 68.5 (正运动学)"
    **正运动学**（forward kinematics）问题是：给定 $n$ 个关节变量 $\theta = (\theta_1, \ldots, \theta_n)$，计算末端执行器相对于基座的齐次变换矩阵 $T(\theta) \in SE(3)$。

!!! theorem "定理 68.8 (指数积公式, Product of Exponentials)"
    设机器人有 $n$ 个关节，每个关节的旋量（在初始构型下）为 $\xi_i$。当关节变量为 $(\theta_1, \ldots, \theta_n)$ 时，末端执行器的位姿为
    $$
    T(\theta) = e^{[\xi_1]\theta_1} \cdot e^{[\xi_2]\theta_2} \cdots e^{[\xi_n]\theta_n} \cdot T_{home}
    $$

    其中 $T_{home} \in SE(3)$ 是所有关节角为零时的末端位姿。

!!! note "注"
    指数积公式（PoE）相比传统的 Denavit-Hartenberg（DH）参数方法有以下优势：

    - **几何直观**：每个 $e^{[\xi_i]\theta_i}$ 直接对应一个关节绕其螺旋轴的运动；
    - **无坐标系选择问题**：不需要在每个连杆上放置坐标系；
    - **数学优雅**：直接使用 $SE(3)$ 的 Lie 群结构；
    - **奇异性处理自然**：通过 Jacobi 矩阵的秩分析。

!!! example "例 68.5"
    **平面两连杆机器人**：两个旋转关节在 $z$ 轴方向。连杆长度为 $L_1, L_2$。

    初始构型：手臂完全伸直沿 $x$ 轴。$T_{home} = \begin{pmatrix} I & (L_1+L_2)\hat{e}_1 \\ \mathbf{0}^T & 1 \end{pmatrix}$。

    关节 1 的旋量（绕原点的 $z$ 轴）：$\xi_1 = (0, 0, 0, 0, 0, 1)^T$。
    关节 2 的旋量（绕 $(L_1, 0, 0)$ 的 $z$ 轴）：$\xi_2 = (0, L_1, 0, 0, 0, 1)^T$。

    正运动学：$T(\theta_1, \theta_2) = e^{[\xi_1]\theta_1} \cdot e^{[\xi_2]\theta_2} \cdot T_{home}$。

    展开计算：
    $$
    T = \begin{pmatrix} \cos(\theta_1+\theta_2) & -\sin(\theta_1+\theta_2) & 0 & L_1\cos\theta_1 + L_2\cos(\theta_1+\theta_2) \\ \sin(\theta_1+\theta_2) & \cos(\theta_1+\theta_2) & 0 & L_1\sin\theta_1 + L_2\sin(\theta_1+\theta_2) \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

### 与 DH 参数的比较

!!! definition "定义 68.6 (Denavit-Hartenberg 参数)"
    传统的 **DH 参数**方法用四个参数 $(\theta_i, d_i, a_i, \alpha_i)$ 描述每个连杆变换：
    $$
    T_i = R_z(\theta_i) T_z(d_i) T_x(a_i) R_x(\alpha_i)
    $$

    正运动学为 $T_{0n} = T_1 T_2 \cdots T_n$。

!!! note "注"
    DH 方法需要在每个关节处仔细放置坐标系（遵循 DH 约定），且坐标系的选择不唯一（当连续两个轴平行或相交时有自由度）。PoE 公式避免了这些问题，更适合现代机器人学。

---

## 68.5 Jacobi 矩阵与速度

<div class="context-flow" markdown>

**核心问题**：关节速度 $\dot{\theta}$ 与末端执行器速度之间的线性关系是什么？

</div>

### 空间 Jacobi 矩阵

!!! definition "定义 68.7 (空间 Jacobi 矩阵)"
    设末端执行器的位姿为 $T(\theta) = e^{[\xi_1]\theta_1} \cdots e^{[\xi_n]\theta_n} T_{home}$。**空间 Jacobi 矩阵**（space Jacobian）$J_s(\theta) \in \mathbb{R}^{6 \times n}$ 定义为
    $$
    [\mathcal{V}_s] = \dot{T} T^{-1} = \sum_{i=1}^n J_{s,i}(\theta) \dot{\theta}_i
    $$

    其中 $\mathcal{V}_s = (\mathbf{v}_s; \omega_s) \in \mathbb{R}^6$ 是空间旋量速度。

!!! theorem "定理 68.9 (空间 Jacobi 矩阵的列)"
    $J_s$ 的第 $i$ 列为
    $$
    J_{s,i}(\theta) = \operatorname{Ad}(e^{[\xi_1]\theta_1} \cdots e^{[\xi_{i-1}]\theta_{i-1}}) \xi_i
    $$

    其中 $\operatorname{Ad}(T)$ 是伴随表示（adjoint representation）：
    $$
    \operatorname{Ad}\begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} = \begin{pmatrix} R & [\mathbf{p}]_\times R \\ 0 & R \end{pmatrix} \in \mathbb{R}^{6 \times 6}
    $$

    特别地，$J_{s,1} = \xi_1$（第一列不依赖关节角）。

### 体 Jacobi 矩阵

!!! definition "定义 68.8 (体 Jacobi 矩阵)"
    **体 Jacobi 矩阵**（body Jacobian）$J_b(\theta)$ 将关节速度映射到末端执行器坐标系中的旋量速度：
    $$
    \mathcal{V}_b = J_b(\theta) \dot{\theta}
    $$

    空间 Jacobi 和体 Jacobi 的关系为
    $$
    J_s(\theta) = \operatorname{Ad}(T(\theta)) J_b(\theta)
    $$

### 速度关系

!!! theorem "定理 68.10 (速度映射)"
    末端执行器的旋量速度与关节速度之间的关系为线性映射：
    $$
    \mathcal{V} = J(\theta) \dot{\theta}
    $$

    这里 $J(\theta)$ 是 $6 \times n$ 矩阵（$n$ 为关节数）。当 $n = 6$ 时 $J$ 为方阵。

### 奇异性

!!! definition "定义 68.9 (运动学奇异性)"
    构型 $\theta^*$ 称为**运动学奇异**（kinematic singularity），若 $\operatorname{rank}(J(\theta^*)) < \min(6, n)$。

    在奇异构型处，末端执行器失去了某些方向上的运动自由度。

!!! example "例 68.6"
    平面两连杆机器人的 Jacobi 矩阵（位置部分）：
    $$
    J_p(\theta) = \begin{pmatrix} -L_1\sin\theta_1 - L_2\sin(\theta_1+\theta_2) & -L_2\sin(\theta_1+\theta_2) \\ L_1\cos\theta_1 + L_2\cos(\theta_1+\theta_2) & L_2\cos(\theta_1+\theta_2) \end{pmatrix}
    $$

    $\det(J_p) = L_1 L_2 \sin\theta_2$。当 $\theta_2 = 0$ 或 $\theta_2 = \pi$（手臂完全伸直或完全折叠）时奇异。

### 操作性椭球

!!! definition "定义 68.10 (操作性椭球)"
    给定关节速度约束 $\|\dot{\theta}\| \leq 1$，末端执行器速度的可达集为
    $$
    \{\mathcal{V} = J\dot{\theta} : \|\dot{\theta}\| \leq 1\}
    $$
    这是一个椭球（当 $J$ 列满秩时），称为**操作性椭球**（manipulability ellipsoid）。

    椭球的半轴由 $J$ 的奇异值决定。**操作性指标** $w(\theta) = \sqrt{\det(J J^T)} = \prod_i \sigma_i$ 量化了构型处的操作灵活性。$w = 0$ 对应奇异构型。

---

## 68.6 逆运动学

<div class="context-flow" markdown>

**核心问题**：给定期望的末端执行器位姿，如何求解关节角度？

</div>

### 问题定义

!!! definition "定义 68.11 (逆运动学)"
    **逆运动学**（inverse kinematics, IK）问题是：给定期望的末端位姿 $T_d \in SE(3)$，求关节角 $\theta \in \mathbb{R}^n$ 使得 $T(\theta) = T_d$。

    这是一个关于 $\theta$ 的非线性方程组。与正运动学不同，逆运动学通常**没有**解析解（除非机器人具有特殊结构），且解可能**不唯一**或**不存在**。

### Newton-Raphson 方法

!!! theorem "定理 68.11 ($SE(3)$ 上的 Newton-Raphson)"
    定义误差旋量 $[\mathcal{V}_b] = \log(T(\theta)^{-1} T_d)$（体坐标下的误差）。Newton-Raphson 迭代为：

    1. 计算误差旋量 $\mathcal{V}_b = \text{vee}(\log(T(\theta_k)^{-1} T_d))$；
    2. 若 $\|\mathcal{V}_b\| < \varepsilon$，停止；
    3. 更新关节角 $\theta_{k+1} = \theta_k + J_b(\theta_k)^{-1} \mathcal{V}_b$（若 $n = 6$ 且 $J_b$ 可逆）。

!!! note "注"
    当 $n \neq 6$ 时：

    - **冗余机器人**（$n > 6$）：使用伪逆 $J_b^+ = J_b^T(J_b J_b^T)^{-1}$，给出最小范数关节速度解。冗余自由度可用于避障或优化其他目标；
    - **欠约束机器人**（$n < 6$）：只能在末端速度的子空间中运动，使用最小二乘近似。

### 几何方法

!!! theorem "定理 68.12 (Pieper 定理)"
    若六自由度机器人的后三个旋转轴交于一点（**球形手腕**），则逆运动学有解析解：

    1. 手腕中心的位置由前三个关节确定（位置子问题）；
    2. 末端姿态由后三个关节确定（姿态子问题——Euler 角分解）。

!!! example "例 68.7"
    **平面两连杆的逆运动学**：给定期望的末端位置 $(x_d, y_d)$。

    由正运动学 $x = L_1\cos\theta_1 + L_2\cos(\theta_1+\theta_2)$，$y = L_1\sin\theta_1 + L_2\sin(\theta_1+\theta_2)$。

    利用 $x^2 + y^2 = L_1^2 + L_2^2 + 2L_1L_2\cos\theta_2$，解出
    $$
    \cos\theta_2 = \frac{x_d^2 + y_d^2 - L_1^2 - L_2^2}{2L_1L_2}
    $$

    由此得 $\theta_2 = \pm\arccos(\cdot)$（两个解——"肘上"和"肘下"）。然后 $\theta_1$ 由三角恒等式确定。解存在当且仅当 $|L_1 - L_2| \leq \sqrt{x_d^2+y_d^2} \leq L_1 + L_2$。

---

## 68.7 动力学

<div class="context-flow" markdown>

**核心问题**：机器人的运动方程具有什么矩阵结构？

</div>

### Euler-Lagrange 方程

!!! theorem "定理 68.13 (机器人动力学的矩阵形式)"
    $n$ 自由度机器人的 Euler-Lagrange 运动方程为
    $$
    M(\theta)\ddot{\theta} + C(\theta, \dot{\theta})\dot{\theta} + G(\theta) = \tau
    $$

    其中：

    - $M(\theta) \in \mathbb{R}^{n \times n}$ 是**质量矩阵**（mass matrix）——正定对称矩阵；
    - $C(\theta, \dot{\theta}) \in \mathbb{R}^{n \times n}$ 是 **Coriolis 和离心力矩阵**；
    - $G(\theta) \in \mathbb{R}^n$ 是**重力项**；
    - $\tau \in \mathbb{R}^n$ 是关节力矩（输入）。

!!! theorem "定理 68.14 (质量矩阵的性质)"
    质量矩阵 $M(\theta)$ 具有以下性质：

    1. **对称性**：$M(\theta) = M(\theta)^T$；
    2. **正定性**：$M(\theta) \succ 0$ 对所有 $\theta$；
    3. **有界性**：存在常数 $0 < m_1 \leq m_2 < \infty$ 使得 $m_1 I \preceq M(\theta) \preceq m_2 I$；
    4. $M(\theta)$ 与 Jacobi 矩阵的关系：$M(\theta) = \sum_{i=1}^n \left(m_i J_{v_i}^T J_{v_i} + J_{\omega_i}^T \mathcal{I}_i J_{\omega_i}\right)$。

!!! definition "定义 68.12 (Christoffel 符号与 Coriolis 矩阵)"
    Coriolis 矩阵 $C(\theta, \dot{\theta})$ 可以用 Christoffel 符号 $c_{ijk}$ 定义：
    $$
    C_{ij}(\theta, \dot{\theta}) = \sum_{k=1}^n c_{ijk}(\theta) \dot{\theta}_k
    $$
    其中
    $$
    c_{ijk} = \frac{1}{2}\left(\frac{\partial M_{ij}}{\partial \theta_k} + \frac{\partial M_{ik}}{\partial \theta_j} - \frac{\partial M_{jk}}{\partial \theta_i}\right)
    $$

!!! theorem "定理 68.15 ($\dot{M} - 2C$ 的反对称性)"
    若 $C$ 按 Christoffel 符号定义，则 $\dot{M}(\theta) - 2C(\theta, \dot{\theta})$ 是反对称矩阵：
    $$
    x^T(\dot{M} - 2C)x = 0, \quad \forall x \in \mathbb{R}^n
    $$

    这一性质在自适应控制和被动性分析中至关重要。

??? proof "证明"
    $\dot{M}_{ij} = \sum_k \frac{\partial M_{ij}}{\partial \theta_k}\dot{\theta}_k$。而 $2C_{ij} = \sum_k \left(\frac{\partial M_{ij}}{\partial \theta_k} + \frac{\partial M_{ik}}{\partial \theta_j} - \frac{\partial M_{jk}}{\partial \theta_i}\right)\dot{\theta}_k$。

    因此 $(\dot{M} - 2C)_{ij} = \sum_k \left(\frac{\partial M_{jk}}{\partial \theta_i} - \frac{\partial M_{ik}}{\partial \theta_j}\right)\dot{\theta}_k = -(\dot{M} - 2C)_{ji}$，即反对称。

!!! example "例 68.8"
    **平面两连杆机器人的动力学**：
    $$
    M(\theta) = \begin{pmatrix} m_1 L_1^2 + m_2(L_1^2 + 2L_1L_2\cos\theta_2 + L_2^2) & m_2(L_1L_2\cos\theta_2 + L_2^2) \\ m_2(L_1L_2\cos\theta_2 + L_2^2) & m_2 L_2^2 \end{pmatrix}
    $$

    Coriolis 矩阵：
    $$
    C(\theta, \dot{\theta}) = \begin{pmatrix} -m_2 L_1 L_2 \sin\theta_2 \, \dot{\theta}_2 & -m_2 L_1 L_2 \sin\theta_2(\dot{\theta}_1 + \dot{\theta}_2) \\ m_2 L_1 L_2 \sin\theta_2 \, \dot{\theta}_1 & 0 \end{pmatrix}
    $$

    重力项：
    $$
    G(\theta) = \begin{pmatrix} (m_1 + m_2)gL_1\cos\theta_1 + m_2 gL_2\cos(\theta_1+\theta_2) \\ m_2 g L_2 \cos(\theta_1 + \theta_2) \end{pmatrix}
    $$

---

## 68.8 姿态估计与 SLAM

<div class="context-flow" markdown>

**核心问题**：如何在 $SO(3)$ 和 $SE(3)$ 上进行估计和优化？

</div>

### 旋转平均

!!! definition "定义 68.13 (旋转平均)"
    给定一组旋转 $R_1, \ldots, R_N \in SO(3)$（可能带噪声），其 **Karcher 均值**（或 Frechet 均值）定义为
    $$
    \bar{R} = \arg\min_{R \in SO(3)} \sum_{i=1}^N d(R, R_i)^2
    $$
    其中 $d(R_1, R_2) = \|\log(R_1^T R_2)\|_F$ 是 $SO(3)$ 上的测地距离。

!!! theorem "定理 68.16 (旋转平均的迭代算法)"
    Karcher 均值可以通过以下迭代计算：

    1. 初始化 $\bar{R}_0$（如取 $R_1$）；
    2. 计算 $\delta = \frac{1}{N}\sum_{i=1}^N \log(\bar{R}_k^T R_i)$；
    3. 更新 $\bar{R}_{k+1} = \bar{R}_k \exp(\delta)$；
    4. 重复直到收敛。

    这本质上是 $SO(3)$ 流形上的梯度下降算法。

### 位姿图优化

!!! definition "定义 68.14 (位姿图优化)"
    在 SLAM（同步定位与建图）中，**位姿图优化**（pose graph optimization）问题为：给定一组相对位姿测量 $\tilde{T}_{ij} \in SE(3)$（来自传感器），求绝对位姿 $T_1, \ldots, T_N \in SE(3)$ 使得
    $$
    \min_{T_1, \ldots, T_N \in SE(3)} \sum_{(i,j) \in \mathcal{E}} \|\log(T_i^{-1} T_j \cdot \tilde{T}_{ij}^{-1})\|_\Sigma^2
    $$

    其中 $\|\xi\|_\Sigma^2 = \xi^T \Sigma^{-1} \xi$ 是 Mahalanobis 范数（$\Sigma$ 是协方差矩阵）。

!!! note "注"
    位姿图优化是一个在 $SE(3)^N$ 上的非线性最小二乘问题。通过在当前估计处线性化（使用 $SE(3)$ 的指数映射），可以将每次迭代归结为一个大型稀疏线性方程组的求解。这正是 g2o、GTSAM 等 SLAM 框架的核心算法。

### 协方差传播

!!! theorem "定理 68.17 (SE(3) 上的协方差传播)"
    若位姿 $T \in SE(3)$ 的不确定性由 Lie 代数中的协方差矩阵 $\Sigma_T \in \mathbb{R}^{6 \times 6}$ 描述，则复合位姿 $T_{12} = T_1 T_2$ 的协方差近似为
    $$
    \Sigma_{12} \approx \Sigma_1 + \operatorname{Ad}(T_1) \Sigma_2 \operatorname{Ad}(T_1)^T
    $$

    其中 $\operatorname{Ad}(T_1) \in \mathbb{R}^{6 \times 6}$ 是伴随表示矩阵。

!!! note "注"
    这一公式是 Kalman 滤波在 $SE(3)$ 上的推广的基础。由于 $SE(3)$ 是非线性流形，协方差传播需要使用伴随表示来处理坐标系的变化。

---

## 本章小结

本章展示了矩阵 Lie 群理论在机器人学中的核心应用。主要内容包括：

1. **刚体运动**由 $SE(3)$ 元素描述——$4 \times 4$ 齐次变换矩阵，包含旋转 $R \in SO(3)$ 和平移 $\mathbf{p} \in \mathbb{R}^3$。

2. **旋量**是 $SE(3)$ 的 Lie 代数 $\mathfrak{se}(3)$ 的元素，描述瞬时刚体运动。Chasles 定理保证每个刚体运动都是螺旋运动。

3. **指数映射**建立了 $\mathfrak{se}(3)$（线性空间）与 $SE(3)$（流形）之间的桥梁，其闭合形式由 Rodrigues 公式给出。

4. **正运动学**的指数积公式 $T = e^{[\xi_1]\theta_1} \cdots e^{[\xi_n]\theta_n} T_{home}$ 直接利用 $SE(3)$ 的群结构。

5. **Jacobi 矩阵**建立了关节速度与末端速度之间的线性映射 $\mathcal{V} = J(\theta)\dot{\theta}$。奇异构型对应 $J$ 的秩不足。

6. **逆运动学**通过 $SE(3)$ 上的 Newton-Raphson 方法或特殊结构的几何方法求解。

7. **动力学**方程 $M(\theta)\ddot{\theta} + C(\theta,\dot{\theta})\dot{\theta} + G(\theta) = \tau$ 具有丰富的矩阵结构，如质量矩阵的正定性和 $\dot{M}-2C$ 的反对称性。

8. **SLAM 和姿态估计**需要在 $SO(3)$ 和 $SE(3)$ 上进行优化和统计推断，核心工具是指数/对数映射和伴随表示。

---

## 习题

!!! question "习题 68.1"
    验证 $SE(3)$ 乘法的结合律：$(T_1 T_2)T_3 = T_1(T_2 T_3)$。

!!! question "习题 68.2"
    设 $T = \begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix}$，其中 $R$ 是绕 $z$ 轴旋转 $90°$，$\mathbf{p} = (1, 0, 0)^T$。计算 $T^{-1}$。

!!! question "习题 68.3"
    对旋量 $\xi = (0, 0, 0, 0, 0, 1)^T$，计算 $e^{[\xi]\theta}$ 并解释其几何含义。

!!! question "习题 68.4"
    写出三连杆平面机器人（三个旋转关节，连杆长度 $L_1, L_2, L_3$）的指数积公式。

!!! question "习题 68.5"
    对例 68.6 中的两连杆机器人，求所有奇异构型并解释其物理含义。

!!! question "习题 68.6"
    设 $R_1$ 和 $R_2$ 是两个旋转矩阵。证明 $d(R_1, R_2) = \|\log(R_1^T R_2)\|_F$ 是一个度量。

!!! question "习题 68.7"
    推导平面两连杆机器人的逆运动学解析解（给定末端位置 $(x_d, y_d)$，求 $\theta_1, \theta_2$）。

!!! question "习题 68.8"
    证明质量矩阵 $M(\theta)$ 是正定的（提示：利用动能的物理意义 $T = \frac{1}{2}\dot{\theta}^T M(\theta)\dot{\theta} > 0$）。

!!! question "习题 68.9"
    对三维刚体的角速度 $\omega$，验证 $[\omega]_\times^3 = -\|\omega\|^2 [\omega]_\times$。

!!! question "习题 68.10"
    设两个位姿测量 $T_1, T_2 \in SE(3)$ 有协方差 $\Sigma_1, \Sigma_2 \in \mathbb{R}^{6 \times 6}$。利用协方差传播公式计算 $T_1 T_2$ 的协方差。
