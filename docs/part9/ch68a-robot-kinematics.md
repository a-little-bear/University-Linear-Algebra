# 第 68A 章 机器人运动学

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 矩阵指数(Ch13) · 矩阵群(Ch55) · 特征值(Ch6) · 奇异值分解(Ch7)

**本章脉络**：$SO(3)$ 与 $SE(3)$ 群结构 → 坐标系变换与链式法则 → Chasles 定理(螺旋运动) → Rodrigues 指数映射 → $SE(3)$ 指数映射 → 对数映射 → Lie 括号与伴随表示 → 指数积公式(PoE) → 空间/体 Jacobi 矩阵 → 力旋量(wrench) → 逆运动学(Newton-Raphson / Pieper / Paden-Kahan) → 操作性椭球

**延伸**：运动学是动力学(Ch68B)的基础；Lie 群/Lie 代数结构在 SLAM、姿态估计(Ch68B)和计算机视觉中无处不在

</div>

机器人运动学是矩阵 Lie 群理论最自然的工程应用。机器人的每一个关节运动都是一个刚体变换，而刚体变换的数学描述恰好是特殊欧几里得群 $SE(3)$——一个 $4 \times 4$ 矩阵 Lie 群。机器人手臂的正运动学是 $SE(3)$ 中矩阵指数的乘积（指数积公式），Jacobi 矩阵建立了关节速度与末端速度之间的线性映射，逆运动学则需要在 $SE(3)$ 上求解非线性方程。

本章从刚体运动的群结构出发，经由旋量理论、指数/对数映射和 Lie 代数工具，系统发展正/逆运动学的数学框架。

---

## 68A.1 旋转群 $SO(3)$

<div class="context-flow" markdown>

**核心问题**：如何用矩阵精确描述三维空间中的旋转？

</div>

### 定义与基本性质

!!! definition "定义 68A.1 (旋转群 $SO(3)$)"
    三维**特殊正交群**定义为
    $$
    SO(3) = \{R \in \mathbb{R}^{3 \times 3} : R^T R = I, \; \det(R) = 1\}
    $$

    $SO(3)$ 的元素表示三维空间中保持原点不动的旋转。它是一个 3 维紧致连通 Lie 群。

!!! theorem "定理 68A.1 ($SO(3)$ 的基本性质)"
    1. **群封闭性**：$R_1, R_2 \in SO(3) \Rightarrow R_1 R_2 \in SO(3)$；
    2. **逆元为转置**：$R^{-1} = R^T$；
    3. **行列式为 1**：$\det(R) = 1$（区别于包含反射的 $O(3)$）；
    4. **特征值**：$R$ 的特征值为 $1, e^{i\theta}, e^{-i\theta}$（$\theta$ 为旋转角），因此 $\operatorname{tr}(R) = 1 + 2\cos\theta$；
    5. **Lie 代数**：$\mathfrak{so}(3) = \{\Omega \in \mathbb{R}^{3 \times 3} : \Omega^T = -\Omega\}$，即 $3 \times 3$ 反对称矩阵集合，维度为 3。

??? proof "证明"
    **群封闭性**：$(R_1R_2)^T(R_1R_2) = R_2^TR_1^TR_1R_2 = R_2^TR_2 = I$，$\det(R_1R_2) = \det(R_1)\det(R_2) = 1$。

    **逆元**：$R^TR = I$ 意味着 $R^{-1} = R^T$。验证 $R^T \in SO(3)$：$(R^T)^TR^T = RR^T = I$（由 $R^TR = I$ 两端左乘 $R$、右乘 $R^T$），$\det(R^T) = \det(R) = 1$。

    **特征值**：$R \in SO(3)$ 是实正交矩阵，特征值的模为 1。实矩阵的复特征值成共轭对。$\det(R) = 1$ 意味着特征值之积为 1。三个模为 1 的数之积为 1，其中复共轭对为 $e^{\pm i\theta}$，则第三个必须为 $1$。迹为特征值之和 $= 1 + 2\cos\theta$。

    **Lie 代数**：$\mathfrak{so}(3)$ 是 $SO(3)$ 在单位元处的切空间。设 $R(t)$ 是 $SO(3)$ 中过 $I$ 的曲线（$R(0) = I$），由 $R(t)^TR(t) = I$ 求导：$\dot{R}(0)^T + \dot{R}(0) = 0$，即 $\dot{R}(0)$ 反对称。反之，任何反对称矩阵 $\Omega$ 生成的 $e^{\Omega t} \in SO(3)$。$\blacksquare$

### 反对称矩阵与叉积

!!! definition "定义 68A.2 (hat 映射)"
    向量 $\omega = (\omega_1, \omega_2, \omega_3)^T \in \mathbb{R}^3$ 与反对称矩阵之间的**hat 映射** $(\cdot)^\wedge : \mathbb{R}^3 \to \mathfrak{so}(3)$ 定义为
    $$
    [\omega]_\times \equiv \omega^\wedge = \begin{pmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{pmatrix}
    $$

    其逆映射 $(\cdot)^\vee : \mathfrak{so}(3) \to \mathbb{R}^3$ 称为 **vee 映射**。

    基本性质：$[\omega]_\times v = \omega \times v$（叉积）。

!!! theorem "定理 68A.2 (反对称矩阵的代数性质)"
    对 $\omega, v \in \mathbb{R}^3$，$R \in SO(3)$：

    1. $[\omega]_\times^2 = \omega\omega^T - \|\omega\|^2 I$；
    2. $[\omega]_\times^3 = -\|\omega\|^2 [\omega]_\times$（周期性）；
    3. $R[\omega]_\times R^T = [R\omega]_\times$（旋转的伴随作用）；
    4. $[\omega]_\times [v]_\times - [v]_\times [\omega]_\times = [\omega \times v]_\times$（Lie 括号）。

---

## 68A.2 特殊欧几里得群 $SE(3)$

<div class="context-flow" markdown>

**核心问题**：如何统一描述三维空间中的旋转和平移？

</div>

### 定义与群性质

!!! definition "定义 68A.3 (特殊欧几里得群 $SE(3)$)"
    三维**特殊欧几里得群**定义为
    $$
    SE(3) = \left\{T = \begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} \in \mathbb{R}^{4 \times 4} : R \in SO(3), \; \mathbf{p} \in \mathbb{R}^3\right\}
    $$

    $SE(3)$ 的元素称为**齐次变换矩阵**（homogeneous transformation matrix），表示刚体运动：旋转 $R$ 加平移 $\mathbf{p}$。

!!! theorem "定理 68A.3 ($SE(3)$ 的群性质)"
    1. **乘法**：$T_1 T_2 = \begin{pmatrix} R_1 R_2 & R_1 \mathbf{p}_2 + \mathbf{p}_1 \\ \mathbf{0}^T & 1 \end{pmatrix} \in SE(3)$；

    2. **逆元**：$T^{-1} = \begin{pmatrix} R^T & -R^T \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix}$；

    3. **单位元**：$I_4$；

    4. $SE(3)$ 是 6 维 Lie 群（3 个旋转自由度 + 3 个平移自由度）。

??? proof "证明"
    **乘法封闭性**：$R_1R_2 \in SO(3)$，$R_1\mathbf{p}_2 + \mathbf{p}_1 \in \mathbb{R}^3$，所以乘积具有 $SE(3)$ 的形式。

    **逆元公式**：验证
    $$
    T \cdot T^{-1} = \begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} \begin{pmatrix} R^T & -R^T\mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} = \begin{pmatrix} RR^T & -RR^T\mathbf{p} + \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} = \begin{pmatrix} I & \mathbf{0} \\ \mathbf{0}^T & 1 \end{pmatrix} = I_4
    $$

    **结合律**：矩阵乘法本身满足结合律。$\blacksquare$

### 坐标系变换

!!! theorem "定理 68A.4 (坐标系变换与链式法则)"
    设 $T_{ab} \in SE(3)$ 表示坐标系 $\{b\}$ 相对于坐标系 $\{a\}$ 的位姿。若点 $\mathbf{q}$ 在坐标系 $\{b\}$ 中的齐次坐标为 $\tilde{\mathbf{q}}_b = (\mathbf{q}_b; 1)$，则它在坐标系 $\{a\}$ 中的坐标为
    $$
    \tilde{\mathbf{q}}_a = T_{ab} \, \tilde{\mathbf{q}}_b
    $$

    若有三个坐标系 $\{a\}, \{b\}, \{c\}$，则
    $$
    T_{ac} = T_{ab} \cdot T_{bc}
    $$
    （坐标系变换的链式法则——下标"内消"）。

!!! example "例 68A.1"
    机器人手臂的末端执行器相对于基座的位姿：
    $$
    T_{0e} = \begin{pmatrix} R_{0e} & \mathbf{p}_{0e} \\ \mathbf{0}^T & 1 \end{pmatrix}
    $$

    若世界坐标系 $\{w\}$、基座 $\{0\}$、末端 $\{e\}$、工具 $\{t\}$，则工具在世界坐标系中的位姿为
    $$
    T_{wt} = T_{w0} \cdot T_{0e} \cdot T_{et}
    $$

---

## 68A.3 旋量与螺旋运动

<div class="context-flow" markdown>

**核心问题**：如何用 Lie 代数语言描述刚体的瞬时运动？

</div>

### Lie 代数 $\mathfrak{se}(3)$

!!! definition "定义 68A.4 (Lie 代数 $\mathfrak{se}(3)$)"
    $SE(3)$ 的 **Lie 代数** $\mathfrak{se}(3)$ 由以下 $4 \times 4$ 矩阵组成：
    $$
    \mathfrak{se}(3) = \left\{[\xi]^\wedge = \begin{pmatrix} [\omega]_\times & \mathbf{v} \\ \mathbf{0}^T & 0 \end{pmatrix} \in \mathbb{R}^{4 \times 4} : \omega \in \mathbb{R}^3, \; \mathbf{v} \in \mathbb{R}^3\right\}
    $$

    六维向量 $\xi = (\mathbf{v}; \omega) \in \mathbb{R}^6$ 称为**旋量**（twist）。

    - $\omega \in \mathbb{R}^3$ 是**角速度**分量；
    - $\mathbf{v} \in \mathbb{R}^3$ 是**线速度**分量（参考点在原点时）。

### Chasles 定理

!!! theorem "定理 68A.5 (Chasles 定理)"
    三维空间中的任意刚体运动（$T \in SE(3)$，$T \neq I$）都可以分解为绕某条轴旋转加上沿该轴的平移。这条轴称为**螺旋轴**（screw axis），运动称为**螺旋运动**（screw motion）。

    更精确地，每个 $T \in SE(3)$（$T \neq I$）都可以写成绕螺旋轴 $\ell$ 旋转角 $\theta$ 并沿 $\ell$ 平移距离 $d$ 的形式。

!!! definition "定义 68A.5 (螺旋轴)"
    **螺旋轴**（screw）$\mathcal{S}$ 由以下参数描述：

    - 轴方向 $\hat{\omega} \in \mathbb{R}^3$（$\|\hat{\omega}\| = 1$）；
    - 轴上一点 $\mathbf{q} \in \mathbb{R}^3$；
    - 螺距 $h$（旋转一弧度时沿轴平移的距离）。

    对应的单位旋量为
    $$
    \xi = \begin{pmatrix} -\hat{\omega} \times \mathbf{q} + h\hat{\omega} \\ \hat{\omega} \end{pmatrix}
    $$

    螺旋运动为 $T = e^{[\xi]^\wedge \theta}$，其中 $\theta$ 是旋转角度。

!!! example "例 68A.2"
    **纯旋转**（$h = 0$）：绕过原点的 $z$ 轴旋转。
    $$
    \xi = \begin{pmatrix} \mathbf{0} \\ \hat{e}_3 \end{pmatrix} = (0, 0, 0, 0, 0, 1)^T
    $$

    **纯平移**（$\|\omega\| = 0$）：沿方向 $\hat{\mathbf{v}}$ 平移距离 $d$。
    $$
    \xi = \begin{pmatrix} \hat{\mathbf{v}} \\ \mathbf{0} \end{pmatrix}, \quad T = \begin{pmatrix} I & d\hat{\mathbf{v}} \\ \mathbf{0}^T & 1 \end{pmatrix}
    $$

    **一般螺旋运动**（$h \neq 0$）：绕过点 $\mathbf{q} = (1, 0, 0)^T$ 的 $z$ 轴旋转并沿 $z$ 方向平移（$h = 0.5$）。
    $$
    \xi = \begin{pmatrix} 0 \\ 1 \\ 0.5 \\ 0 \\ 0 \\ 1 \end{pmatrix}
    $$

---

## 68A.4 指数映射

<div class="context-flow" markdown>

**核心问题**：如何在 Lie 代数（线性空间）和 Lie 群（流形）之间建立联系？

</div>

### $SO(3)$ 的指数映射：Rodrigues 公式

!!! theorem "定理 68A.6 (Rodrigues 公式——$SO(3)$ 的指数映射)"
    对 $\omega \in \mathbb{R}^3$，$\hat{\omega} = \omega/\|\omega\|$（$\omega \neq 0$），$\theta = \|\omega\|$：
    $$
    \exp([\omega]_\times) = I + \sin\theta \, [\hat{\omega}]_\times + (1 - \cos\theta) [\hat{\omega}]_\times^2
    $$

    等价地写成：
    $$
    e^{[\hat{\omega}]_\times \theta} = I + \sin\theta \, [\hat{\omega}]_\times + (1 - \cos\theta)(\hat{\omega}\hat{\omega}^T - I) = \cos\theta \cdot I + (1-\cos\theta)\hat{\omega}\hat{\omega}^T + \sin\theta \, [\hat{\omega}]_\times
    $$

    这是从 $\mathfrak{so}(3)$ 到 $SO(3)$ 的指数映射，几何含义为绕轴 $\hat{\omega}$ 旋转角度 $\theta$。

??? proof "证明"
    利用幂级数 $e^{[\hat{\omega}]_\times \theta} = \sum_{k=0}^{\infty} \frac{\theta^k}{k!} [\hat{\omega}]_\times^k$。

    **关键递推关系**（$\|\hat{\omega}\| = 1$）：
    $$
    [\hat{\omega}]_\times^2 = \hat{\omega}\hat{\omega}^T - I, \quad [\hat{\omega}]_\times^3 = -[\hat{\omega}]_\times
    $$

    因此对 $k \geq 1$：
    $$
    [\hat{\omega}]_\times^{2k} = (-1)^{k+1} [\hat{\omega}]_\times^2, \quad [\hat{\omega}]_\times^{2k+1} = (-1)^k [\hat{\omega}]_\times
    $$

    代入幂级数：
    $$
    e^{[\hat{\omega}]_\times \theta} = I + \left(\sum_{k=0}^{\infty} \frac{(-1)^k \theta^{2k+1}}{(2k+1)!}\right) [\hat{\omega}]_\times + \left(\sum_{k=1}^{\infty} \frac{(-1)^{k+1} \theta^{2k}}{(2k)!}\right) [\hat{\omega}]_\times^2
    $$
    $$
    = I + \sin\theta \, [\hat{\omega}]_\times + (1 - \cos\theta) [\hat{\omega}]_\times^2
    $$

    **验证结果在 $SO(3)$ 中**：设 $R = e^{[\hat{\omega}]_\times \theta}$。由 $[\hat{\omega}]_\times^T = -[\hat{\omega}]_\times$：
    $$
    R^T = I - \sin\theta \, [\hat{\omega}]_\times + (1-\cos\theta)[\hat{\omega}]_\times^2
    $$
    $$
    R^TR = (I + (1-\cos\theta)[\hat{\omega}]_\times^2)^2 + \sin^2\theta \, [\hat{\omega}]_\times^2 \cdot (-[\hat{\omega}]_\times^2)
    $$
    经仔细展开利用 $[\hat{\omega}]_\times^4 = -[\hat{\omega}]_\times^2$ 可验证 $R^TR = I$。$\det(R) = 1$ 因为 $\det(e^X) = e^{\operatorname{tr}(X)} = e^0 = 1$。$\blacksquare$

### $SE(3)$ 的指数映射

!!! theorem "定理 68A.7 ($SE(3)$ 的指数映射)"
    对旋量 $\xi = (\mathbf{v}; \omega) \in \mathbb{R}^6$，$\theta = \|\omega\|$（$\omega \neq 0$）：
    $$
    \exp([\xi]^\wedge \theta) = \begin{pmatrix} e^{[\omega]_\times \theta} & J(\hat{\omega}, \theta) \mathbf{v}\theta \\ \mathbf{0}^T & 1 \end{pmatrix}
    $$

    其中 $J(\hat{\omega}, \theta)$ 是 **左 Jacobi 矩阵**：
    $$
    J(\hat{\omega}, \theta) = \frac{\sin\theta}{\theta}I + \left(1 - \frac{\sin\theta}{\theta}\right)\hat{\omega}\hat{\omega}^T + \frac{1-\cos\theta}{\theta}[\hat{\omega}]_\times
    $$

    等价形式：
    $$
    J = I + \frac{1 - \cos\theta}{\theta} [\hat{\omega}]_\times + \frac{\theta - \sin\theta}{\theta} [\hat{\omega}]_\times^2
    $$

    对纯平移情形（$\omega = 0$）：
    $$
    \exp\left(\begin{pmatrix} 0 & \mathbf{v}\theta \\ \mathbf{0}^T & 0 \end{pmatrix}\right) = \begin{pmatrix} I & \mathbf{v}\theta \\ \mathbf{0}^T & 1 \end{pmatrix}
    $$

??? proof "证明"
    设 $[\xi]^\wedge = \begin{pmatrix} [\omega]_\times & \mathbf{v} \\ \mathbf{0}^T & 0 \end{pmatrix}$。利用 $4 \times 4$ 矩阵的幂级数：

    $$
    ([\xi]^\wedge)^k = \begin{pmatrix} [\omega]_\times^k & [\omega]_\times^{k-1}\mathbf{v} \\ \mathbf{0}^T & 0 \end{pmatrix}, \quad k \geq 1
    $$

    因此
    $$
    e^{[\xi]^\wedge \theta} = I_4 + \sum_{k=1}^{\infty} \frac{\theta^k}{k!} ([\xi]^\wedge)^k = \begin{pmatrix} e^{[\omega]_\times\theta} & \left(\sum_{k=1}^{\infty} \frac{\theta^k}{k!} [\omega]_\times^{k-1}\right)\mathbf{v} \\ \mathbf{0}^T & 1 \end{pmatrix}
    $$

    右上角的求和 $\sum_{k=1}^{\infty} \frac{\theta^k}{k!} [\omega]_\times^{k-1} = \theta J(\hat{\omega}, \theta)$（利用 $[\hat{\omega}]_\times$ 的周期性展开并识别三角函数）。$\blacksquare$

---

## 68A.5 对数映射

<div class="context-flow" markdown>

**核心问题**：如何从旋转矩阵或齐次变换矩阵恢复 Lie 代数元素？

</div>

### $SO(3)$ 的对数映射

!!! theorem "定理 68A.8 ($SO(3)$ 的对数映射)"
    给定 $R \in SO(3)$，其对数 $[\omega]_\times = \log(R)$ 为：

    - 若 $R = I$：$\omega = \mathbf{0}$；
    - 若 $\operatorname{tr}(R) \neq -1$（即 $\theta \neq \pi$）：
    $$
    \theta = \arccos\frac{\operatorname{tr}(R) - 1}{2}, \quad [\hat{\omega}]_\times = \frac{R - R^T}{2\sin\theta}
    $$
    - 若 $\operatorname{tr}(R) = -1$（$\theta = \pi$）：$\hat{\omega}$ 是 $R + I$ 的任意非零列的归一化。

!!! example "例 68A.3"
    设 $R = \begin{pmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}$（绕 $z$ 轴旋转 $90°$）。

    $\operatorname{tr}(R) = 0 + 0 + 1 = 1$，$\theta = \arccos\frac{1-1}{2} = \arccos 0 = \frac{\pi}{2}$。

    $[\hat{\omega}]_\times = \frac{R - R^T}{2\sin(\pi/2)} = \frac{1}{2}\begin{pmatrix} 0 & -2 & 0 \\ 2 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}$

    读出 $\hat{\omega} = (0, 0, 1)^T$，确实是 $z$ 轴方向。$\omega = \theta \hat{\omega} = (\pi/2)(0, 0, 1)^T$。

### $SE(3)$ 的对数映射

!!! theorem "定理 68A.9 ($SE(3)$ 的对数映射)"
    给定 $T = \begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} \in SE(3)$，其对数 $[\xi]^\wedge \theta = \log(T)$ 的计算步骤为：

    1. 由 $R$ 计算 $\hat{\omega}$ 和 $\theta$（使用 $SO(3)$ 对数映射）；
    2. 若 $\theta \neq 0$，计算左 Jacobi 矩阵 $J(\hat{\omega}, \theta)$；
    3. 求解 $\mathbf{v} = J^{-1} \mathbf{p}/\theta$，其中
    $$
    J^{-1} = \frac{1}{\theta}\left(I - \frac{\theta}{2}[\hat{\omega}]_\times + \left(1 - \frac{\theta}{2}\cot\frac{\theta}{2}\right)[\hat{\omega}]_\times^2\right)
    $$
    4. 旋量为 $\xi = (\mathbf{v}; \hat{\omega})$。

    若 $R = I$（纯平移），则 $\omega = \mathbf{0}$，$\mathbf{v} = \mathbf{p}/\|\mathbf{p}\|$，$\theta = \|\mathbf{p}\|$。

!!! example "例 68A.4"
    设 $T = \begin{pmatrix} 0 & -1 & 0 & 1 \\ 1 & 0 & 0 & 2 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$。

    $R$ 部分同例 68A.3：$\hat{\omega} = (0,0,1)^T$，$\theta = \pi/2$。

    计算 $J^{-1}$，然后 $\mathbf{v} = J^{-1}(1, 2, 0)^T / (\pi/2)$，得到对应的旋量 $\xi$。

---

## 68A.6 Lie 括号与伴随表示

<div class="context-flow" markdown>

**核心问题**：$\mathfrak{se}(3)$ 上的代数运算如何刻画刚体运动的微分结构？

</div>

### Lie 括号

!!! definition "定义 68A.6 (Lie 括号 / 交换子)"
    $\mathfrak{se}(3)$ 上的 **Lie 括号**（交换子）定义为
    $$
    [[\xi_1]^\wedge, [\xi_2]^\wedge] = [\xi_1]^\wedge [\xi_2]^\wedge - [\xi_2]^\wedge [\xi_1]^\wedge
    $$

    结果仍在 $\mathfrak{se}(3)$ 中。用六维向量表示：若 $\xi_1 = (\mathbf{v}_1; \omega_1)$，$\xi_2 = (\mathbf{v}_2; \omega_2)$，则
    $$
    [\xi_1, \xi_2] = \begin{pmatrix} \omega_1 \times \mathbf{v}_2 - \omega_2 \times \mathbf{v}_1 \\ \omega_1 \times \omega_2 \end{pmatrix}
    $$

!!! theorem "定理 68A.10 (Lie 括号的性质)"
    1. **双线性**：$[\alpha\xi_1 + \beta\xi_2, \xi_3] = \alpha[\xi_1, \xi_3] + \beta[\xi_2, \xi_3]$；
    2. **反对称性**：$[\xi_1, \xi_2] = -[\xi_2, \xi_1]$；
    3. **Jacobi 恒等式**：$[\xi_1, [\xi_2, \xi_3]] + [\xi_2, [\xi_3, \xi_1]] + [\xi_3, [\xi_1, \xi_2]] = 0$；
    4. 对 $\mathfrak{so}(3)$ 子代数：$[[\omega_1]_\times, [\omega_2]_\times] = [\omega_1 \times \omega_2]_\times$（Lie 括号就是叉积）。

??? proof "证明"
    **Jacobi 恒等式**：直接由矩阵交换子的定义验证。设 $A = [\xi_1]^\wedge$，$B = [\xi_2]^\wedge$，$C = [\xi_3]^\wedge$。

    $$
    [A, [B, C]] + [B, [C, A]] + [C, [A, B]]
    $$
    $$
    = A(BC - CB) - (BC - CB)A + B(CA - AC) - (CA - AC)B + C(AB - BA) - (AB - BA)C
    $$

    展开后所有 24 项（$ABC, ACB, BAC, BCA, CAB, CBA$ 各出现 4 次，正负相消）恰好为零。$\blacksquare$

### 伴随表示

!!! definition "定义 68A.7 (大伴随表示 $\operatorname{Ad}$)"
    $SE(3)$ 的**伴随表示**（adjoint representation）将群元素 $T \in SE(3)$ 映射为 $\mathfrak{se}(3)$ 上的线性变换：
    $$
    \operatorname{Ad}_T : \mathfrak{se}(3) \to \mathfrak{se}(3), \quad \operatorname{Ad}_T([\xi]^\wedge) = T [\xi]^\wedge T^{-1}
    $$

    用 $6 \times 6$ 矩阵表示：
    $$
    \operatorname{Ad}(T) = \operatorname{Ad}\begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix} = \begin{pmatrix} R & [\mathbf{p}]_\times R \\ 0 & R \end{pmatrix} \in \mathbb{R}^{6 \times 6}
    $$

    即 $\operatorname{Ad}(T) \xi$ 将旋量从一个坐标系变换到另一个坐标系。

!!! theorem "定理 68A.11 (伴随表示的性质)"
    1. **群同态**：$\operatorname{Ad}(T_1 T_2) = \operatorname{Ad}(T_1) \operatorname{Ad}(T_2)$；
    2. **逆元**：$\operatorname{Ad}(T^{-1}) = \operatorname{Ad}(T)^{-1}$；
    3. **指数兼容**：$\operatorname{Ad}(e^{[\xi]^\wedge}) = e^{\operatorname{ad}_\xi}$，其中 $\operatorname{ad}_\xi$ 是小伴随映射。

!!! definition "定义 68A.8 (小伴随表示 $\operatorname{ad}$)"
    **小伴随映射** $\operatorname{ad}_\xi : \mathfrak{se}(3) \to \mathfrak{se}(3)$ 定义为 Lie 括号：
    $$
    \operatorname{ad}_\xi(\eta) = [\xi, \eta]
    $$

    用 $6 \times 6$ 矩阵表示（$\xi = (\mathbf{v}; \omega)$）：
    $$
    \operatorname{ad}_\xi = [\operatorname{ad}_\xi] = \begin{pmatrix} [\omega]_\times & [\mathbf{v}]_\times \\ 0 & [\omega]_\times \end{pmatrix} \in \mathbb{R}^{6 \times 6}
    $$

    $\operatorname{ad}$ 是 $\operatorname{Ad}$ 在单位元处的微分：$\operatorname{ad}_\xi = \frac{d}{dt}\bigg|_{t=0} \operatorname{Ad}(e^{[\xi]^\wedge t})$。

---

## 68A.7 正运动学：指数积公式

<div class="context-flow" markdown>

**核心问题**：给定关节角度，如何计算末端执行器的位姿？

</div>

### Product of Exponentials (PoE) 公式

!!! definition "定义 68A.9 (正运动学问题)"
    **正运动学**（forward kinematics）是：给定 $n$ 个关节变量 $\theta = (\theta_1, \ldots, \theta_n)$，计算末端执行器相对于基座的齐次变换矩阵 $T(\theta) \in SE(3)$。

!!! theorem "定理 68A.12 (指数积公式)"
    设机器人有 $n$ 个关节，每个关节在初始构型下的旋量为 $\xi_i \in \mathbb{R}^6$。当关节变量为 $(\theta_1, \ldots, \theta_n)$ 时，末端执行器的位姿为
    $$
    T(\theta) = e^{[\xi_1]^\wedge\theta_1} \cdot e^{[\xi_2]^\wedge\theta_2} \cdots e^{[\xi_n]^\wedge\theta_n} \cdot T_{home}
    $$

    其中 $T_{home} \in SE(3)$ 是所有关节角为零时的末端位姿。

!!! note "注"
    指数积公式（PoE）相比传统的 Denavit-Hartenberg（DH）参数方法有以下优势：

    - **几何直观**：每个 $e^{[\xi_i]^\wedge\theta_i}$ 直接对应一个关节绕其螺旋轴的运动；
    - **无坐标系选择问题**：不需要在每个连杆上放置坐标系；
    - **数学优雅**：直接使用 $SE(3)$ 的 Lie 群结构；
    - **奇异性处理自然**：通过 Jacobi 矩阵的秩分析。

!!! example "例 68A.5"
    **平面两连杆机器人**：两个旋转关节在 $z$ 轴方向，连杆长度为 $L_1, L_2$。

    初始构型：手臂完全伸直沿 $x$ 轴。$T_{home} = \begin{pmatrix} I & (L_1+L_2)\hat{e}_1 \\ \mathbf{0}^T & 1 \end{pmatrix}$。

    关节 1 的旋量（绕原点的 $z$ 轴）：$\xi_1 = (0, 0, 0, 0, 0, 1)^T$。

    关节 2 的旋量（绕 $(L_1, 0, 0)$ 的 $z$ 轴）：$\xi_2 = (0, L_1, 0, 0, 0, 1)^T$。

    正运动学：$T(\theta_1, \theta_2) = e^{[\xi_1]^\wedge\theta_1} \cdot e^{[\xi_2]^\wedge\theta_2} \cdot T_{home}$。

    展开计算得：
    $$
    T = \begin{pmatrix} \cos(\theta_1+\theta_2) & -\sin(\theta_1+\theta_2) & 0 & L_1\cos\theta_1 + L_2\cos(\theta_1+\theta_2) \\ \sin(\theta_1+\theta_2) & \cos(\theta_1+\theta_2) & 0 & L_1\sin\theta_1 + L_2\sin(\theta_1+\theta_2) \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

---

## 68A.8 Jacobi 矩阵与速度映射

<div class="context-flow" markdown>

**核心问题**：关节速度 $\dot{\theta}$ 与末端执行器速度之间的线性关系是什么？

</div>

### 空间 Jacobi 矩阵

!!! definition "定义 68A.10 (空间 Jacobi 矩阵)"
    设末端位姿 $T(\theta) = e^{[\xi_1]^\wedge\theta_1} \cdots e^{[\xi_n]^\wedge\theta_n} T_{home}$。**空间 Jacobi 矩阵** $J_s(\theta) \in \mathbb{R}^{6 \times n}$ 定义为
    $$
    \mathcal{V}_s = J_s(\theta) \dot{\theta}
    $$
    其中 $\mathcal{V}_s = (\mathbf{v}_s; \omega_s)$ 是空间坐标系下的旋量速度，由 $[\mathcal{V}_s]^\wedge = \dot{T}T^{-1}$ 定义。

!!! theorem "定理 68A.13 (空间 Jacobi 矩阵的列)"
    $J_s$ 的第 $i$ 列为
    $$
    J_{s,i}(\theta) = \operatorname{Ad}(e^{[\xi_1]^\wedge\theta_1} \cdots e^{[\xi_{i-1}]^\wedge\theta_{i-1}}) \xi_i
    $$

    特别地，$J_{s,1} = \xi_1$（第一列不依赖关节角）。

??? proof "证明"
    对 $T(\theta) = e^{[\xi_1]^\wedge\theta_1} \cdots e^{[\xi_n]^\wedge\theta_n} T_{home}$ 关于时间求导：
    $$
    \dot{T} = \sum_{i=1}^n e^{[\xi_1]^\wedge\theta_1} \cdots e^{[\xi_{i-1}]^\wedge\theta_{i-1}} [\xi_i]^\wedge e^{[\xi_i]^\wedge\theta_i} \cdots e^{[\xi_n]^\wedge\theta_n} T_{home} \cdot \dot{\theta}_i
    $$

    左乘 $T^{-1}$ 不便，改为右乘：
    $$
    \dot{T}T^{-1} = \sum_{i=1}^n e^{[\xi_1]^\wedge\theta_1} \cdots e^{[\xi_{i-1}]^\wedge\theta_{i-1}} [\xi_i]^\wedge \left(e^{[\xi_1]^\wedge\theta_1} \cdots e^{[\xi_{i-1}]^\wedge\theta_{i-1}}\right)^{-1} \dot{\theta}_i
    $$

    利用 $T[\xi]^\wedge T^{-1} = [\operatorname{Ad}(T)\xi]^\wedge$，得
    $$
    [\mathcal{V}_s]^\wedge = \sum_{i=1}^n [\operatorname{Ad}(e^{[\xi_1]^\wedge\theta_1} \cdots e^{[\xi_{i-1}]^\wedge\theta_{i-1}}) \xi_i]^\wedge \dot{\theta}_i
    $$

    读出 $J_{s,i} = \operatorname{Ad}(e^{[\xi_1]^\wedge\theta_1} \cdots e^{[\xi_{i-1}]^\wedge\theta_{i-1}}) \xi_i$。$\blacksquare$

### 体 Jacobi 矩阵

!!! definition "定义 68A.11 (体 Jacobi 矩阵)"
    **体 Jacobi 矩阵** $J_b(\theta)$ 将关节速度映射到末端执行器坐标系中的旋量速度：
    $$
    \mathcal{V}_b = J_b(\theta) \dot{\theta}, \quad [\mathcal{V}_b]^\wedge = T^{-1}\dot{T}
    $$

    空间 Jacobi 和体 Jacobi 的关系为
    $$
    J_s(\theta) = \operatorname{Ad}(T(\theta)) J_b(\theta)
    $$

### 力旋量（Wrench）

!!! definition "定义 68A.12 (力旋量)"
    **力旋量**（wrench）$\mathcal{F} = (\mathbf{f}; \boldsymbol{\tau}_0) \in \mathbb{R}^6$ 是旋量的**对偶**（dual），其中 $\mathbf{f}$ 是力，$\boldsymbol{\tau}_0$ 是相对于参考点的力矩。

    力旋量做的瞬时功为
    $$
    P = \mathcal{F}^T \mathcal{V} = \mathbf{f}^T \mathbf{v} + \boldsymbol{\tau}_0^T \omega
    $$

!!! theorem "定理 68A.14 (力旋量的变换)"
    力旋量在坐标系变换 $T$ 下的变换为
    $$
    \mathcal{F}_a = \operatorname{Ad}(T_{ab})^{-T} \mathcal{F}_b
    $$

    这是因为功率不变：$\mathcal{F}_a^T \mathcal{V}_a = \mathcal{F}_b^T \mathcal{V}_b$，而 $\mathcal{V}_a = \operatorname{Ad}(T_{ab})\mathcal{V}_b$。

    关节力矩 $\tau \in \mathbb{R}^n$ 与末端力旋量 $\mathcal{F}$ 的关系为（虚功原理）：
    $$
    \tau = J(\theta)^T \mathcal{F}
    $$

### 运动学奇异性

!!! definition "定义 68A.13 (运动学奇异性)"
    构型 $\theta^*$ 称为**运动学奇异**（kinematic singularity），若 $\operatorname{rank}(J(\theta^*)) < \min(6, n)$。

    在奇异构型处，末端执行器失去了某些方向上的运动自由度。

!!! example "例 68A.6"
    平面两连杆机器人的 Jacobi 矩阵（位置部分）：
    $$
    J_p(\theta) = \begin{pmatrix} -L_1\sin\theta_1 - L_2\sin(\theta_1+\theta_2) & -L_2\sin(\theta_1+\theta_2) \\ L_1\cos\theta_1 + L_2\cos(\theta_1+\theta_2) & L_2\cos(\theta_1+\theta_2) \end{pmatrix}
    $$

    $\det(J_p) = L_1 L_2 \sin\theta_2$。当 $\theta_2 = 0$（完全伸直）或 $\theta_2 = \pi$（完全折叠）时奇异。

---

## 68A.9 操作性椭球与灵巧性

<div class="context-flow" markdown>

**核心问题**：如何量化机器人在不同构型处的运动能力？

</div>

!!! definition "定义 68A.14 (操作性椭球)"
    给定关节速度约束 $\|\dot{\theta}\| \leq 1$，末端执行器速度的可达集为
    $$
    \mathcal{E} = \{\mathcal{V} = J\dot{\theta} : \|\dot{\theta}\| \leq 1\}
    $$
    当 $J$ 列满秩时，这是一个椭球，称为**操作性椭球**（manipulability ellipsoid）。

    椭球的半轴方向由 $J$ 的左奇异向量决定，半轴长度为奇异值 $\sigma_1 \geq \cdots \geq \sigma_r$。

!!! definition "定义 68A.15 (操作性指标)"
    **Yoshikawa 操作性指标**：
    $$
    w(\theta) = \sqrt{\det(J J^T)} = \prod_{i=1}^r \sigma_i
    $$

    $w = 0$ 对应奇异构型。$w$ 越大，椭球体积越大，操作灵活性越好。

!!! definition "定义 68A.16 (条件数灵巧性)"
    **条件数灵巧指标**：
    $$
    \kappa(\theta) = \frac{\sigma_{\max}(J)}{\sigma_{\min}(J)}
    $$

    $\kappa = 1$ 时操作性椭球为球——所有方向运动能力相同，称为**各向同性**（isotropic）构型。$\kappa$ 越大，椭球越"扁"，某些方向运动困难。

!!! example "例 68A.7"
    对平面两连杆机器人（$L_1 = L_2 = 1$），$J_p^T J_p$ 的特征值决定操作性。当 $\theta_2 = \pi/2$ 时 $\kappa$ 较小（接近各向同性），当 $\theta_2 \to 0$ 或 $\theta_2 \to \pi$ 时 $\kappa \to \infty$（趋向奇异）。

---

## 68A.10 逆运动学

<div class="context-flow" markdown>

**核心问题**：给定期望的末端执行器位姿，如何求解关节角度？

</div>

### Newton-Raphson 方法

!!! definition "定义 68A.17 (逆运动学问题)"
    **逆运动学**（inverse kinematics, IK）：给定期望末端位姿 $T_d \in SE(3)$，求关节角 $\theta \in \mathbb{R}^n$ 使得 $T(\theta) = T_d$。

    这是关于 $\theta$ 的非线性方程组。通常**没有**解析解，且解可能**不唯一**或**不存在**。

!!! theorem "定理 68A.15 ($SE(3)$ 上的 Newton-Raphson 方法)"
    定义体坐标误差旋量 $[\mathcal{V}_b]^\wedge = \log(T(\theta)^{-1} T_d)$。迭代算法为：

    1. 计算误差旋量 $\mathcal{V}_b = \operatorname{vee}(\log(T(\theta_k)^{-1} T_d))$；
    2. 若 $\|\omega_b\| < \varepsilon_\omega$ 且 $\|\mathbf{v}_b\| < \varepsilon_v$，停止；
    3. 更新关节角：
        - $n = 6$ 且 $J_b$ 可逆：$\theta_{k+1} = \theta_k + J_b(\theta_k)^{-1} \mathcal{V}_b$
        - $n > 6$（冗余）：$\theta_{k+1} = \theta_k + J_b^+ \mathcal{V}_b$，$J_b^+ = J_b^T(J_b J_b^T)^{-1}$
        - $n < 6$（欠约束）：$\theta_{k+1} = \theta_k + J_b^T(J_b J_b^T)^{-1} \mathcal{V}_b$（最小二乘）

!!! note "注"
    此方法直接在 $SE(3)$ 流形上工作——误差通过对数映射在 Lie 代数中表示，更新通过 Jacobi 矩阵在关节空间中进行。这比在欧几里得空间中分别处理位置和姿态更自然、收敛更快。

### Pieper 定理

!!! theorem "定理 68A.16 (Pieper 定理)"
    若六自由度机器人的后三个旋转轴交于一点（**球形手腕**），则逆运动学有解析解：

    1. **位置子问题**：手腕中心的位置 $\mathbf{p}_w = \mathbf{p}_d - d_6 R_d \hat{e}_z$ 由前三个关节确定；
    2. **姿态子问题**：末端姿态由后三个关节确定（Euler 角分解）。

### Paden-Kahan 子问题

!!! definition "定义 68A.18 (Paden-Kahan 子问题)"
    Paden-Kahan 子问题是可以解析求解的基本几何 IK 子问题：

    **子问题 1**（绕单轴旋转一个点到目标点）：给定轴 $\hat{\omega}$、点 $\mathbf{p}$ 和目标 $\mathbf{q}$（$\|\mathbf{p}\| = \|\mathbf{q}\|$ 关于轴的投影），求 $\theta$ 使得 $e^{[\hat{\omega}]_\times \theta} \mathbf{p} = \mathbf{q}$。

    **子问题 2**（绕两个相交轴旋转一个点到目标点）：求 $\theta_1, \theta_2$ 使得 $e^{[\hat{\omega}_1]_\times \theta_1} e^{[\hat{\omega}_2]_\times \theta_2} \mathbf{p} = \mathbf{q}$。

    **子问题 3**（绕单轴旋转到与目标点给定距离）：求 $\theta$ 使得 $\|e^{[\hat{\omega}]_\times \theta}\mathbf{p} - \mathbf{q}\| = d$。

!!! theorem "定理 68A.17 (Paden-Kahan 子问题的解)"
    **子问题 1 的解**：将 $\mathbf{p}$ 和 $\mathbf{q}$ 投影到垂直于 $\hat{\omega}$ 的平面，$\theta$ 为投影向量之间的有符号角。
    $$
    \theta = \operatorname{atan2}(\hat{\omega}^T(\mathbf{p}' \times \mathbf{q}'), \; \mathbf{p}'^T \mathbf{q}')
    $$
    其中 $\mathbf{p}' = \mathbf{p} - \hat{\omega}\hat{\omega}^T\mathbf{p}$，$\mathbf{q}' = \mathbf{q} - \hat{\omega}\hat{\omega}^T\mathbf{q}$。

    **子问题 2** 有 0 或 2 个解（通过将问题归约为两个子问题 1）。

    **子问题 3** 有 0、1 或 2 个解（归约为平面几何问题）。

!!! example "例 68A.8"
    **平面两连杆的逆运动学**（即子问题 2 + 子问题 3 的组合）：

    给定期望末端位置 $(x_d, y_d)$。由正运动学：
    $$
    x = L_1\cos\theta_1 + L_2\cos(\theta_1+\theta_2), \quad y = L_1\sin\theta_1 + L_2\sin(\theta_1+\theta_2)
    $$

    利用 $x^2 + y^2 = L_1^2 + L_2^2 + 2L_1L_2\cos\theta_2$：
    $$
    \cos\theta_2 = \frac{x_d^2 + y_d^2 - L_1^2 - L_2^2}{2L_1L_2}
    $$

    $\theta_2 = \pm\arccos(\cdot)$（两个解——"肘上"和"肘下"）。$\theta_1$ 随后由 atan2 确定。

    解存在当且仅当 $|L_1 - L_2| \leq \sqrt{x_d^2+y_d^2} \leq L_1 + L_2$。

---

## 本章小结

本章展示了矩阵 Lie 群理论在机器人运动学中的核心应用。主要内容包括：

1. **旋转群 $SO(3)$** 描述三维旋转，逆元为转置，Lie 代数 $\mathfrak{so}(3)$ 是反对称矩阵空间。

2. **特殊欧几里得群 $SE(3)$** 用 $4 \times 4$ 齐次变换矩阵统一描述旋转和平移，坐标系变换满足链式法则。

3. **旋量**是 $\mathfrak{se}(3)$ 的元素，描述刚体瞬时运动。**Chasles 定理**保证每个刚体运动都是螺旋运动。

4. **Rodrigues 公式**给出 $SO(3)$ 的指数映射闭合形式，$SE(3)$ 的指数映射利用左 Jacobi 矩阵 $J$。对数映射是其逆运算。

5. **Lie 括号**和**伴随表示**（$\operatorname{Ad}$ 和 $\operatorname{ad}$）是 Lie 代数的核心工具，伴随表示用于旋量在不同坐标系间的变换。

6. **指数积公式** $T = e^{[\xi_1]^\wedge\theta_1} \cdots e^{[\xi_n]^\wedge\theta_n} T_{home}$ 直接利用 $SE(3)$ 群结构描述正运动学。

7. **Jacobi 矩阵**建立关节速度与末端旋量速度的线性映射。**力旋量**是旋量的对偶，通过 $\operatorname{Ad}^{-T}$ 变换。

8. **逆运动学**通过 $SE(3)$ 上的 Newton-Raphson 方法（数值）或 **Pieper 定理**、**Paden-Kahan 子问题**（解析/几何）求解。

9. **操作性椭球**和条件数量化了构型处的运动灵巧性，与 Jacobi 矩阵的奇异值直接相关。

---

## 习题

!!! question "习题 68A.1"
    验证 $SE(3)$ 乘法的结合律：$(T_1 T_2)T_3 = T_1(T_2 T_3)$。

!!! question "习题 68A.2"
    设 $T = \begin{pmatrix} R & \mathbf{p} \\ \mathbf{0}^T & 1 \end{pmatrix}$，其中 $R$ 是绕 $z$ 轴旋转 $90°$，$\mathbf{p} = (1, 0, 0)^T$。计算 $T^{-1}$ 并验证 $TT^{-1} = I_4$。

!!! question "习题 68A.3"
    对旋量 $\xi = (0, 0, 0, 0, 0, 1)^T$，计算 $e^{[\xi]^\wedge\theta}$ 并解释其几何含义。

!!! question "习题 68A.4"
    对三维刚体的角速度 $\omega$，验证 $[\omega]_\times^3 = -\|\omega\|^2 [\omega]_\times$。

!!! question "习题 68A.5"
    写出三连杆平面机器人（三个旋转关节，连杆长度 $L_1, L_2, L_3$）的指数积公式，并计算其空间 Jacobi 矩阵。

!!! question "习题 68A.6"
    对例 68A.6 中的两连杆机器人，求所有奇异构型。在奇异构型处，操作性椭球退化为什么形状？

!!! question "习题 68A.7"
    验证伴随表示的群同态性质：对 $T_1, T_2 \in SE(3)$，直接计算验证 $\operatorname{Ad}(T_1T_2) = \operatorname{Ad}(T_1)\operatorname{Ad}(T_2)$。

!!! question "习题 68A.8"
    推导平面两连杆机器人的逆运动学解析解。给定 $L_1 = 1$，$L_2 = 1$，$(x_d, y_d) = (1, 1)$，求所有可能的 $(\theta_1, \theta_2)$。

!!! question "习题 68A.9"
    证明 $d(R_1, R_2) = \|\log(R_1^T R_2)\|_F$ 是 $SO(3)$ 上的度量（满足非负性、对称性、三角不等式）。

!!! question "习题 68A.10"
    设 $\xi_1 = (1, 0, 0, 0, 0, 1)^T$，$\xi_2 = (0, 1, 0, 0, 0, 1)^T$。计算 Lie 括号 $[\xi_1, \xi_2]$ 并解释其几何含义。

!!! question "习题 68A.11"
    对 SCARA 机器人（两个旋转关节绕 $z$ 轴，一个棱柱关节沿 $z$ 轴，一个旋转关节绕 $z$ 轴），写出指数积公式。分析其奇异构型。

!!! question "习题 68A.12"
    给定 $T_d \in SE(3)$，用 Newton-Raphson 方法求解两连杆平面机器人的逆运动学。写出迭代步骤，并讨论初始猜测对收敛的影响。

!!! question "习题 68A.13"
    验证 Jacobi 恒等式：对 $\mathfrak{se}(3)$ 中的三个元素 $\xi_1, \xi_2, \xi_3$，直接计算验证 $[\xi_1, [\xi_2, \xi_3]] + [\xi_2, [\xi_3, \xi_1]] + [\xi_3, [\xi_1, \xi_2]] = 0$。

!!! question "习题 68A.14"
    设六自由度机器人在某构型处的 Jacobi 矩阵 $J$ 的奇异值为 $\sigma_1 = 5, \sigma_2 = 3, \sigma_3 = 2, \sigma_4 = 1, \sigma_5 = 0.5, \sigma_6 = 0.1$。计算操作性指标 $w$ 和条件数 $\kappa$。讨论该构型的运动学特性。
