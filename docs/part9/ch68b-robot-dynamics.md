# 第 68B 章 机器人动力学与状态估计

<div class="context-flow" markdown>

**前置**：$SO(3)$ 与 $SE(3)$ 群结构(Ch68A) · 旋量与指数映射(Ch68A) · Lie 括号与伴随表示(Ch68A) · Jacobi 矩阵(Ch68A) · 矩阵指数(Ch13) · 奇异值分解(Ch7) · 四元数(Ch56)

**本章脉络**：Euler-Lagrange 动力学 → 质量矩阵性质 → $\dot{M}-2C$ 反对称性 → 递归 Newton-Euler 算法 → 计算力矩控制 → BCH 公式 → 对偶四元数 → 旋转平均(Karcher 均值) → $SE(3)$ 协方差传播 → 位姿图优化(SLAM) → 流形上的扩展 Kalman 滤波

**延伸**：动力学是机器人控制（力控、阻抗控制、自适应控制）的基础；状态估计方法推广到视觉惯性里程计(VIO)和自动驾驶

</div>

机器人动力学将牛顿力学与矩阵 Lie 群理论结合。Euler-Lagrange 方程 $M(\theta)\ddot{\theta} + C(\theta, \dot{\theta})\dot{\theta} + G(\theta) = \tau$ 具有丰富的矩阵结构——质量矩阵的正定性和 $\dot{M}-2C$ 的反对称性是被动性分析和控制设计的基础。递归 Newton-Euler 算法是计算动力学的最高效方法。

在状态估计方面，现代机器人学需要在 $SO(3)$ 和 $SE(3)$ 流形上进行优化和滤波。Baker-Campbell-Hausdorff 公式描述了 Lie 代数元素组合的规律，旋转平均和位姿图优化在 SLAM 中至关重要，流形上的扩展 Kalman 滤波是传感器融合的核心工具。

---

## 68B.1 Euler-Lagrange 动力学

<div class="context-flow" markdown>

**核心问题**：机器人的运动方程具有什么矩阵结构？

</div>

### 运动方程

!!! theorem "定理 68B.1 (机器人动力学的矩阵形式)"
    $n$ 自由度机器人的 Euler-Lagrange 运动方程为
    $$
    M(\theta)\ddot{\theta} + C(\theta, \dot{\theta})\dot{\theta} + G(\theta) = \tau
    $$

    其中：

    - $M(\theta) \in \mathbb{R}^{n \times n}$ 是**质量矩阵**（mass matrix / inertia matrix）——正定对称矩阵；
    - $C(\theta, \dot{\theta}) \in \mathbb{R}^{n \times n}$ 是 **Coriolis 和离心力矩阵**；
    - $G(\theta) \in \mathbb{R}^n$ 是**重力项**；
    - $\tau \in \mathbb{R}^n$ 是关节力矩（控制输入）。

    此方程由 Lagrange 函数 $\mathcal{L} = T - V$ 导出，其中动能 $T = \frac{1}{2}\dot{\theta}^T M(\theta)\dot{\theta}$，势能 $V = V(\theta)$。

### 质量矩阵的性质

!!! theorem "定理 68B.2 (质量矩阵的性质)"
    质量矩阵 $M(\theta)$ 具有以下性质：

    1. **对称性**：$M(\theta) = M(\theta)^T$；
    2. **正定性**：$M(\theta) \succ 0$ 对所有 $\theta$；
    3. **有界性**：存在常数 $0 < m_1 \leq m_2 < \infty$ 使得 $m_1 I \preceq M(\theta) \preceq m_2 I$；
    4. **与 Jacobi 矩阵的关系**：
    $$
    M(\theta) = \sum_{i=1}^n \left(m_i J_{v_i}(\theta)^T J_{v_i}(\theta) + J_{\omega_i}(\theta)^T \mathcal{I}_i J_{\omega_i}(\theta)\right)
    $$
    其中 $J_{v_i}$ 和 $J_{\omega_i}$ 分别是第 $i$ 个连杆质心的线速度和角速度 Jacobi 矩阵，$\mathcal{I}_i$ 是连杆惯性张量。

??? proof "证明"
    **正定性**：动能 $T = \frac{1}{2}\dot{\theta}^T M(\theta)\dot{\theta}$ 对非零 $\dot{\theta}$ 必须严格正（物理上，运动的刚体必须有正动能）。

    更严格地，由性质 4 中的公式：
    $$
    \dot{\theta}^T M(\theta) \dot{\theta} = \sum_{i=1}^n \left(m_i \|J_{v_i}\dot{\theta}\|^2 + \dot{\theta}^T J_{\omega_i}^T \mathcal{I}_i J_{\omega_i}\dot{\theta}\right)
    $$

    每一项非负。若 $\dot{\theta}^T M \dot{\theta} = 0$，则所有连杆的线速度和角速度都为零，即 $J_{v_i}\dot{\theta} = 0$ 且 $J_{\omega_i}\dot{\theta} = 0$ 对所有 $i$。对于非退化的机器人构型，这蕴含 $\dot{\theta} = 0$（否则某个连杆必然在运动）。

    **对称性**：由公式直接可见，$M$ 是对称矩阵之和。

    **有界性**：由 $M(\theta)$ 关于 $\theta$ 的连续性和紧集上的极值定理。对于旋转关节，$\theta$ 在 $[0, 2\pi]^n$ 上变化（紧集），$M(\theta)$ 的特征值连续且正，因此有正的下界和有限的上界。$\blacksquare$

### Coriolis 矩阵

!!! definition "定义 68B.1 (Christoffel 符号与 Coriolis 矩阵)"
    Coriolis 矩阵 $C(\theta, \dot{\theta})$ 可以用 **Christoffel 符号** $c_{ijk}$ 定义：
    $$
    C_{ij}(\theta, \dot{\theta}) = \sum_{k=1}^n c_{ijk}(\theta) \dot{\theta}_k
    $$
    其中
    $$
    c_{ijk} = \frac{1}{2}\left(\frac{\partial M_{ij}}{\partial \theta_k} + \frac{\partial M_{ik}}{\partial \theta_j} - \frac{\partial M_{jk}}{\partial \theta_i}\right)
    $$

    此定义确保 Euler-Lagrange 方程的正确性，同时赋予 $C$ 重要的反对称性质。

### $\dot{M} - 2C$ 的反对称性

!!! theorem "定理 68B.3 ($\dot{M} - 2C$ 的反对称性)"
    若 $C$ 按 Christoffel 符号定义，则 $N(\theta, \dot{\theta}) = \dot{M}(\theta) - 2C(\theta, \dot{\theta})$ 是反对称矩阵：
    $$
    x^T(\dot{M} - 2C)x = 0, \quad \forall x \in \mathbb{R}^n
    $$

    即 $N = N^T$... 不，准确地说 $N^T = -N$（反对称）。

??? proof "证明"
    计算 $N$ 的 $(i,j)$ 元素。$\dot{M}$ 的 $(i,j)$ 元素为
    $$
    \dot{M}_{ij} = \sum_k \frac{\partial M_{ij}}{\partial \theta_k}\dot{\theta}_k
    $$

    $2C$ 的 $(i,j)$ 元素为
    $$
    2C_{ij} = \sum_k \left(\frac{\partial M_{ij}}{\partial \theta_k} + \frac{\partial M_{ik}}{\partial \theta_j} - \frac{\partial M_{jk}}{\partial \theta_i}\right)\dot{\theta}_k
    $$

    因此
    $$
    N_{ij} = \dot{M}_{ij} - 2C_{ij} = \sum_k \left(\frac{\partial M_{jk}}{\partial \theta_i} - \frac{\partial M_{ik}}{\partial \theta_j}\right)\dot{\theta}_k
    $$

    交换 $i, j$：
    $$
    N_{ji} = \sum_k \left(\frac{\partial M_{ik}}{\partial \theta_j} - \frac{\partial M_{jk}}{\partial \theta_i}\right)\dot{\theta}_k = -N_{ij}
    $$

    因此 $N^T = -N$，即 $\dot{M} - 2C$ 是反对称的。

    推论：$x^T N x = x^T(-N^T)x = -(x^T N x)$，故 $x^T(\dot{M} - 2C)x = 0$。$\blacksquare$

!!! note "注"
    $\dot{M} - 2C$ 的反对称性意味着系统具有**被动性**（passivity）：沿运动轨线，
    $$
    \frac{d}{dt}\left(\frac{1}{2}\dot{\theta}^T M \dot{\theta}\right) = \dot{\theta}^T M\ddot{\theta} + \frac{1}{2}\dot{\theta}^T\dot{M}\dot{\theta} = \dot{\theta}^T(\tau - G) + \frac{1}{2}\dot{\theta}^T(\dot{M} - 2C)\dot{\theta} = \dot{\theta}^T(\tau - G)
    $$

    动能的变化率等于外力（减去重力）做的功率。这一性质在自适应控制、阻抗控制和力控制设计中至关重要。

!!! example "例 68B.1"
    **平面两连杆机器人的动力学**（质量集中在连杆末端）：
    $$
    M(\theta) = \begin{pmatrix} m_1 L_1^2 + m_2(L_1^2 + 2L_1L_2\cos\theta_2 + L_2^2) & m_2(L_1L_2\cos\theta_2 + L_2^2) \\ m_2(L_1L_2\cos\theta_2 + L_2^2) & m_2 L_2^2 \end{pmatrix}
    $$

    Coriolis 矩阵：
    $$
    C(\theta, \dot{\theta}) = \begin{pmatrix} -m_2 L_1 L_2 \sin\theta_2 \, \dot{\theta}_2 & -m_2 L_1 L_2 \sin\theta_2(\dot{\theta}_1 + \dot{\theta}_2) \\ m_2 L_1 L_2 \sin\theta_2 \, \dot{\theta}_1 & 0 \end{pmatrix}
    $$

    验证：$\dot{M} - 2C$ 的 $(1,2)$ 元素为 $-2m_2L_1L_2\sin\theta_2\,\dot{\theta}_2 - 2(-m_2L_1L_2\sin\theta_2(\dot{\theta}_1+\dot{\theta}_2)) = 2m_2L_1L_2\sin\theta_2\,\dot{\theta}_1$，$(2,1)$ 元素为 $-2m_2L_1L_2\sin\theta_2\,\dot{\theta}_1$。确实反对称。

---

## 68B.2 递归 Newton-Euler 算法

<div class="context-flow" markdown>

**核心问题**：如何高效地计算机器人的逆动力学（给定运动求力矩）？

</div>

!!! definition "定义 68B.2 (逆动力学问题)"
    **逆动力学**（inverse dynamics）：给定关节角 $\theta$、角速度 $\dot{\theta}$ 和角加速度 $\ddot{\theta}$，计算所需的关节力矩 $\tau$：
    $$
    \tau = M(\theta)\ddot{\theta} + C(\theta, \dot{\theta})\dot{\theta} + G(\theta)
    $$

!!! theorem "定理 68B.4 (递归 Newton-Euler 算法)"
    逆动力学可以通过以下 $O(n)$ 复杂度的递归算法计算（$n$ 为关节数）：

    **前向递推**（从基座到末端，$i = 1, \ldots, n$）：

    计算每个连杆的角速度、角加速度和线加速度：
    $$
    \omega_i = R_{i,i-1}^T \omega_{i-1} + \dot{\theta}_i \hat{z}_i
    $$
    $$
    \dot{\omega}_i = R_{i,i-1}^T \dot{\omega}_{i-1} + R_{i,i-1}^T \omega_{i-1} \times \dot{\theta}_i \hat{z}_i + \ddot{\theta}_i \hat{z}_i
    $$
    $$
    \ddot{p}_i = R_{i,i-1}^T(\ddot{p}_{i-1} + \dot{\omega}_{i-1} \times r_{i-1,i} + \omega_{i-1} \times (\omega_{i-1} \times r_{i-1,i}))
    $$

    初始条件：$\omega_0 = 0$，$\dot{\omega}_0 = 0$，$\ddot{p}_0 = -\mathbf{g}$（重力加速度项）。

    **后向递推**（从末端到基座，$i = n, \ldots, 1$）：

    计算每个连杆上的力和力矩：
    $$
    f_i = m_i \ddot{p}_{c_i} + R_{i,i+1} f_{i+1}
    $$
    $$
    n_i = \mathcal{I}_i \dot{\omega}_i + \omega_i \times \mathcal{I}_i \omega_i + R_{i,i+1} n_{i+1} + r_{i,c_i} \times m_i \ddot{p}_{c_i} + r_{i,i+1} \times R_{i,i+1} f_{i+1}
    $$

    关节力矩：$\tau_i = n_i^T \hat{z}_i$（旋转关节）或 $\tau_i = f_i^T \hat{z}_i$（棱柱关节）。

    末端条件：$f_{n+1} = 0$，$n_{n+1} = 0$（或由末端外力给定）。

!!! note "注"
    递归 Newton-Euler 算法的计算复杂度为 $O(n)$，远优于直接构造 $M, C, G$ 矩阵再做矩阵-向量乘法的 $O(n^3)$ 方法。它是工业机器人控制器中实时计算逆动力学的标准方法。

    该算法也可用于高效计算 $M(\theta)$（对 $\ddot{\theta} = e_i$，$\dot{\theta} = 0$ 调用 $n$ 次得到 $M$ 的第 $i$ 列）和 $C\dot{\theta} + G$（令 $\ddot{\theta} = 0$）。

---

## 68B.3 计算力矩控制

<div class="context-flow" markdown>

**核心问题**：如何利用精确的动力学模型实现关节空间的精确轨迹跟踪？

</div>

!!! definition "定义 68B.3 (计算力矩控制 / 逆动力学控制)"
    **计算力矩控制**（computed torque control）选取控制律
    $$
    \tau = M(\theta)\left(\ddot{\theta}_d + K_D(\dot{\theta}_d - \dot{\theta}) + K_P(\theta_d - \theta)\right) + C(\theta, \dot{\theta})\dot{\theta} + G(\theta)
    $$

    其中 $\theta_d(t)$ 是期望轨迹，$K_P, K_D \succ 0$ 是增益矩阵。

!!! theorem "定理 68B.5 (计算力矩控制的稳定性)"
    设 $e = \theta_d - \theta$ 为跟踪误差。代入计算力矩控制律后，误差动力学为
    $$
    \ddot{e} + K_D \dot{e} + K_P e = 0
    $$

    这是线性的常系数 ODE。选择 $K_P, K_D$ 使特征值在左半平面，则 $e(t) \to 0$（指数收敛）。

??? proof "证明"
    将控制律代入动力学方程 $M\ddot{\theta} + C\dot{\theta} + G = \tau$：
    $$
    M\ddot{\theta} + C\dot{\theta} + G = M(\ddot{\theta}_d + K_D\dot{e} + K_Pe) + C\dot{\theta} + G
    $$

    消去 $C\dot{\theta} + G$，得
    $$
    M(\ddot{\theta} - \ddot{\theta}_d - K_D\dot{e} - K_Pe) = 0
    $$

    由 $M \succ 0$，$\ddot{\theta} - \ddot{\theta}_d = -K_D\dot{e} - K_Pe$，即 $\ddot{e} + K_D\dot{e} + K_Pe = 0$。$\blacksquare$

!!! note "注"
    计算力矩控制的关键优势是将非线性动力学精确**反馈线性化**——消去了所有非线性项 $C\dot{\theta} + G$，使误差动力学成为线性的。代价是需要精确的动力学模型和高效的实时逆动力学计算（递归 Newton-Euler 算法）。

---

## 68B.4 Baker-Campbell-Hausdorff 公式

<div class="context-flow" markdown>

**核心问题**：两个 Lie 代数元素的指数之积如何用单个 Lie 代数元素的指数表示？

</div>

!!! theorem "定理 68B.6 (Baker-Campbell-Hausdorff (BCH) 公式)"
    对 Lie 代数元素 $X, Y$（足够小时），$e^X e^Y = e^Z$，其中
    $$
    Z = X + Y + \frac{1}{2}[X, Y] + \frac{1}{12}([X, [X, Y]] + [Y, [Y, X]]) + \cdots
    $$

    前几项展开：
    $$
    Z = X + Y + \frac{1}{2}[X, Y] + \frac{1}{12}[X, [X, Y]] - \frac{1}{12}[Y, [X, Y]] + \cdots
    $$

    **一阶近似**（$X, Y$ 都很小时）：$Z \approx X + Y$（指数映射在零点附近近似为加法）。

    **二阶近似**：$Z \approx X + Y + \frac{1}{2}[X, Y]$。

!!! note "注"
    对 $SO(3)$：$[X, Y] = XY - YX$ 对应角速度的叉积。BCH 公式解释了为什么旋转不满足交换律——$e^X e^Y \neq e^Y e^X$，差异由 Lie 括号 $[X, Y]$ 度量。

    对 $SE(3)$：BCH 公式用于组合两个小运动（如连续两帧之间的位姿增量），在 SLAM 和里程计中频繁使用。

!!! theorem "定理 68B.7 (BCH 公式的收敛性)"
    BCH 级数在 $\|X\| + \|Y\| < \log 2$ 时收敛（对矩阵 Lie 群，使用算子范数）。

    在机器人学应用中，通常只使用前两三项，因为旋量增量 $X, Y$ 都很小（单个时间步长的运动增量）。

!!! example "例 68B.2"
    在 $SO(3)$ 上，设 $X = \theta_1 [\hat{e}_3]_\times$（绕 $z$ 轴旋转 $\theta_1$），$Y = \theta_2 [\hat{e}_1]_\times$（绕 $x$ 轴旋转 $\theta_2$）。

    $[X, Y] = \theta_1\theta_2([\hat{e}_3]_\times[\hat{e}_1]_\times - [\hat{e}_1]_\times[\hat{e}_3]_\times) = \theta_1\theta_2[\hat{e}_3 \times \hat{e}_1]_\times = \theta_1\theta_2[\hat{e}_2]_\times$

    二阶近似：$Z \approx \theta_1[\hat{e}_3]_\times + \theta_2[\hat{e}_1]_\times + \frac{\theta_1\theta_2}{2}[\hat{e}_2]_\times$

    这表明先绕 $z$ 轴再绕 $x$ 轴旋转，与简单地绕两轴之和旋转相比，多了一个绕 $y$ 轴的修正项——这就是旋转不可交换性的定量体现。

---

## 68B.5 对偶四元数与刚体运动

<div class="context-flow" markdown>

**核心问题**：除了 $4 \times 4$ 齐次矩阵，有没有更紧凑的刚体运动表示？

</div>

### 四元数回顾

!!! definition "定义 68B.4 (单位四元数与旋转)"
    **单位四元数** $q = (q_0, \mathbf{q}) \in \mathbb{H}$（$\|q\| = q_0^2 + \|\mathbf{q}\|^2 = 1$）表示旋转：
    $$
    q = \cos\frac{\theta}{2} + \sin\frac{\theta}{2}\hat{\omega}
    $$

    旋转作用：$\mathbf{p}' = q \mathbf{p} q^*$（四元数乘法）。

    单位四元数群 $S^3$ 到 $SO(3)$ 的映射是 2-1 同态：$q$ 和 $-q$ 对应同一旋转。

### 对偶四元数

!!! definition "定义 68B.5 (对偶四元数)"
    **对偶四元数**（dual quaternion）定义为
    $$
    \hat{q} = q_r + \varepsilon q_d
    $$

    其中 $q_r, q_d \in \mathbb{H}$ 是普通四元数，$\varepsilon$ 是**对偶单位**（$\varepsilon^2 = 0$，$\varepsilon \neq 0$）。

    表示刚体运动（旋转 $q_r$ + 平移 $\mathbf{t}$）的对偶四元数为
    $$
    \hat{q} = q_r + \frac{\varepsilon}{2} t \cdot q_r
    $$

    其中 $t = (0, \mathbf{t})$ 是纯四元数。

!!! theorem "定理 68B.8 (对偶四元数的运算)"
    1. **乘法**：$\hat{q}_1 \hat{q}_2 = (q_{r1}q_{r2}) + \varepsilon(q_{r1}q_{d2} + q_{d1}q_{r2})$，对应刚体运动的复合。
    2. **共轭**：$\hat{q}^* = q_r^* + \varepsilon q_d^*$。
    3. **点的变换**：点 $\mathbf{p}$ 的对偶四元数表示为 $\hat{p} = 1 + \varepsilon(0, \mathbf{p})$，变换后 $\hat{p}' = \hat{q}\hat{p}\hat{q}^*$。
    4. **与 $SE(3)$ 的关系**：单位对偶四元数群同构于 $SE(3)$ 的双覆盖。

!!! note "注"
    对偶四元数相比 $4 \times 4$ 矩阵的优势：

    - **存储紧凑**：8 个参数（带 2 个约束）vs. 16 个参数（带 10 个约束）；
    - **插值自然**：可以使用 ScLERP（螺旋线性插值）在两个刚体运动之间平滑插值；
    - **无奇异性**：不存在万向锁问题（Euler 角的弊病）；
    - **乘法高效**：四元数乘法比矩阵乘法更快。

    劣势：不如矩阵直观，文献中不如齐次矩阵普及。

---

## 68B.6 旋转平均（Karcher 均值）

<div class="context-flow" markdown>

**核心问题**：如何对一组带噪声的旋转测量求"平均旋转"？

</div>

!!! definition "定义 68B.6 (旋转的 Karcher 均值)"
    给定一组旋转 $R_1, \ldots, R_N \in SO(3)$（可能带权重 $w_i > 0$），其 **Karcher 均值**（Frechet 均值）定义为
    $$
    \bar{R} = \arg\min_{R \in SO(3)} \sum_{i=1}^N w_i \, d(R, R_i)^2
    $$

    其中 $d(R_1, R_2) = \|\log(R_1^T R_2)\|_F$ 是 $SO(3)$ 上的测地距离（双不变度量）。

!!! theorem "定理 68B.9 (旋转平均的迭代算法)"
    Karcher 均值可以通过以下**Riemannian 梯度下降**迭代计算：

    1. 初始化 $\bar{R}_0$（如取 $R_1$）；
    2. 计算切向量（负梯度方向）：
    $$
    \delta = \sum_{i=1}^N w_i \log(\bar{R}_k^T R_i) \bigg/ \sum_{i=1}^N w_i
    $$
    3. 更新 $\bar{R}_{k+1} = \bar{R}_k \exp(\delta)$；
    4. 重复直到 $\|\delta\|_F < \varepsilon$。

    当所有 $R_i$ 在 $\bar{R}$ 的一个测地球内时，算法线性收敛。

??? proof "证明"
    目标函数 $f(R) = \frac{1}{2}\sum_i w_i \|\log(R^T R_i)\|_F^2$。

    在 $R$ 处参数化：$R(\delta) = R \exp(\delta)$，$\delta \in \mathfrak{so}(3)$。

    $f(R\exp(\delta))$ 关于 $\delta$ 的一阶展开（利用 $\log((R\exp(\delta))^T R_i) = \log(\exp(-\delta) R^TR_i) \approx \log(R^TR_i) - \delta$，对小 $\delta$）：

    梯度为 $\nabla f = -\sum_i w_i \log(R^T R_i)$。令梯度为零的一阶必要条件为 $\sum_i w_i \log(\bar{R}^T R_i) = 0$。

    迭代步 $\delta = -\nabla f / \sum w_i$ 是归一化的负梯度方向。$\blacksquare$

!!! example "例 68B.3"
    设 $N = 3$，$R_1$ 绕 $z$ 轴旋转 $5°$，$R_2$ 绕 $z$ 轴旋转 $-3°$，$R_3$ 绕 $z$ 轴旋转 $10°$（等权重 $w_i = 1$）。

    所有旋转都在 $z$ 轴附近，Karcher 均值近似为绕 $z$ 轴旋转 $(5° - 3° + 10°)/3 = 4°$。

    对于大角度或不同轴的旋转，简单的"角度平均"不适用——必须使用 Riemannian 梯度下降。

---

## 68B.7 $SE(3)$ 上的协方差传播

<div class="context-flow" markdown>

**核心问题**：两个不确定的刚体运动复合后，不确定性如何传播？

</div>

!!! definition "定义 68B.7 ($SE(3)$ 上的不确定性模型)"
    位姿 $T \in SE(3)$ 的不确定性用 Lie 代数中的随机扰动描述：
    $$
    T = \bar{T} \exp([\epsilon]^\wedge), \quad \epsilon \sim \mathcal{N}(0, \Sigma)
    $$

    其中 $\bar{T}$ 是均值位姿，$\epsilon \in \mathbb{R}^6$ 是小扰动，$\Sigma \in \mathbb{R}^{6 \times 6}$ 是**协方差矩阵**。

!!! theorem "定理 68B.10 ($SE(3)$ 上的协方差传播)"
    若位姿 $T_1 = \bar{T}_1 \exp([\epsilon_1]^\wedge)$，$T_2 = \bar{T}_2 \exp([\epsilon_2]^\wedge)$，$\epsilon_1, \epsilon_2$ 独立且均值为零，则复合位姿 $T_{12} = T_1 T_2$ 的协方差近似为
    $$
    \Sigma_{12} \approx \Sigma_1 + \operatorname{Ad}(\bar{T}_1) \Sigma_2 \operatorname{Ad}(\bar{T}_1)^T
    $$

    其中 $\operatorname{Ad}(\bar{T}_1) \in \mathbb{R}^{6 \times 6}$ 是伴随表示矩阵。

??? proof "证明"
    $$
    T_{12} = T_1 T_2 = \bar{T}_1 \exp([\epsilon_1]^\wedge) \bar{T}_2 \exp([\epsilon_2]^\wedge)
    $$

    利用 $\exp([\epsilon_1]^\wedge) \bar{T}_2 = \bar{T}_2 \exp([\operatorname{Ad}(\bar{T}_2^{-1})\epsilon_1]^\wedge)$（一阶近似）：

    实际上更直接地：
    $$
    T_{12} = \bar{T}_1\bar{T}_2 \cdot (\bar{T}_1\bar{T}_2)^{-1} T_1 T_2 = \bar{T}_{12} \cdot \bar{T}_2^{-1}\exp([\epsilon_1]^\wedge)\bar{T}_2 \exp([\epsilon_2]^\wedge)
    $$

    但这不完全正确。正确的一阶分析：

    将 $T_1 T_2$ 写成 $\bar{T}_1\bar{T}_2 \exp([\epsilon_{12}]^\wedge)$，需要找 $\epsilon_{12}$。

    $T_1T_2 = \bar{T}_1 e^{[\epsilon_1]^\wedge}\bar{T}_2 e^{[\epsilon_2]^\wedge}$。设 $\bar{T}_{12} = \bar{T}_1\bar{T}_2$。

    $$
    \bar{T}_{12}^{-1}T_1T_2 = \bar{T}_2^{-1}\bar{T}_1^{-1}\bar{T}_1 e^{[\epsilon_1]^\wedge}\bar{T}_2 e^{[\epsilon_2]^\wedge} = \bar{T}_2^{-1}e^{[\epsilon_1]^\wedge}\bar{T}_2 e^{[\epsilon_2]^\wedge}
    $$

    利用 $\bar{T}_2^{-1} e^{[\epsilon_1]^\wedge}\bar{T}_2 = e^{[\operatorname{Ad}(\bar{T}_2^{-1})\epsilon_1]^\wedge}$。

    一阶近似（BCH 公式的一阶项）：
    $$
    e^{[\operatorname{Ad}(\bar{T}_2^{-1})\epsilon_1]^\wedge}e^{[\epsilon_2]^\wedge} \approx e^{[\operatorname{Ad}(\bar{T}_2^{-1})\epsilon_1 + \epsilon_2]^\wedge}
    $$

    因此 $\epsilon_{12} \approx \operatorname{Ad}(\bar{T}_2^{-1})\epsilon_1 + \epsilon_2$。

    协方差：$\Sigma_{12} = \operatorname{Ad}(\bar{T}_2^{-1})\Sigma_1 \operatorname{Ad}(\bar{T}_2^{-1})^T + \Sigma_2$。

    若采用"左扰动"模型（$T = \exp([\epsilon]^\wedge)\bar{T}$），则公式变为 $\Sigma_{12} \approx \Sigma_1 + \operatorname{Ad}(\bar{T}_1)\Sigma_2\operatorname{Ad}(\bar{T}_1)^T$。$\blacksquare$

!!! note "注"
    协方差传播公式是 $SE(3)$ 上 Kalman 滤波的基础。伴随表示 $\operatorname{Ad}(T)$ 在这里扮演了欧几里得空间中坐标旋转矩阵的角色——将一个坐标系中的不确定性"旋转"到另一个坐标系。

---

## 68B.8 位姿图优化与 SLAM

<div class="context-flow" markdown>

**核心问题**：如何从一组含噪声的相对位姿测量中恢复全局一致的位姿估计？

</div>

### 问题定义

!!! definition "定义 68B.8 (位姿图优化)"
    在 SLAM（同步定位与建图）中，机器人在 $N$ 个位置采集数据，得到一组相对位姿测量 $\{\tilde{T}_{ij}\}_{(i,j) \in \mathcal{E}}$（来自激光、视觉或 IMU）。**位姿图优化**求绝对位姿 $T_1, \ldots, T_N \in SE(3)$ 使得
    $$
    \min_{T_1, \ldots, T_N \in SE(3)} \sum_{(i,j) \in \mathcal{E}} \|\log(T_i^{-1} T_j \cdot \tilde{T}_{ij}^{-1})\|_{\Sigma_{ij}}^2
    $$

    其中 $\|\xi\|_\Sigma^2 = \xi^T \Sigma^{-1} \xi$ 是 Mahalanobis 范数，$\Sigma_{ij}$ 是测量协方差。

### Gauss-Newton 方法

!!! theorem "定理 68B.11 (位姿图的 Gauss-Newton 迭代)"
    在当前估计 $\{T_i^{(k)}\}$ 处线性化。设增量 $\delta_i \in \mathbb{R}^6$（$T_i \leftarrow T_i^{(k)} \exp([\delta_i]^\wedge)$），残差为
    $$
    r_{ij}(\delta) = \log\left((T_i^{(k)} e^{[\delta_i]^\wedge})^{-1} (T_j^{(k)} e^{[\delta_j]^\wedge}) \tilde{T}_{ij}^{-1}\right)^\vee
    $$

    一阶展开：$r_{ij}(\delta) \approx r_{ij}^{(k)} + J_{ij}^i \delta_i + J_{ij}^j \delta_j$，其中 $J_{ij}^i, J_{ij}^j \in \mathbb{R}^{6 \times 6}$ 是 Jacobi 矩阵。

    Gauss-Newton 方程为大型稀疏线性系统：
    $$
    H \delta^* = -b
    $$
    其中 $H = \sum_{(i,j)} J_{ij}^T \Sigma_{ij}^{-1} J_{ij}$ 是近似 Hessian（$6N \times 6N$ 稀疏矩阵），$b = \sum_{(i,j)} J_{ij}^T \Sigma_{ij}^{-1} r_{ij}^{(k)}$。

    更新：$T_i^{(k+1)} = T_i^{(k)} \exp([\delta_i^*]^\wedge)$。

!!! note "注"
    $H$ 是**稀疏**的——只有直接观测连接的位姿对 $(i,j) \in \mathcal{E}$ 贡献非零块。这使得可以使用稀疏 Cholesky 分解或共轭梯度法高效求解。g2o、GTSAM、Ceres Solver 等开源框架都实现了这一算法。

    位姿图优化的核心数学结构是：**在 Lie 群流形上做非线性最小二乘**，每次迭代通过指数/对数映射在 Lie 代数中线性化，然后求解一个大型稀疏线性方程组。

!!! example "例 68B.4"
    一个简单的 SLAM 问题：机器人走一个正方形回路，有 4 个位姿和 5 条边（包括一条回环闭合边）。

    回环闭合边提供了全局约束，使得位姿估计全局一致。没有回环闭合时，误差会累积（里程计漂移）。位姿图优化通过同时优化所有位姿来分配回环闭合的校正量。

---

## 68B.9 流形上的扩展 Kalman 滤波

<div class="context-flow" markdown>

**核心问题**：如何在 $SO(3)$ 或 $SE(3)$ 上进行递推状态估计？

</div>

### 标准 EKF 的问题

!!! note "注"
    标准的扩展 Kalman 滤波（EKF）假设状态空间是欧几里得空间 $\mathbb{R}^n$。直接将旋转矩阵的 9 个元素或 Euler 角作为状态会导致问题：

    - 旋转矩阵有 6 个正交约束，EKF 更新后可能不满足；
    - Euler 角有万向锁奇异性。

    解决方案是**在 Lie 代数中定义误差状态**，利用流形结构进行 EKF。

### $SO(3)$ 上的 EKF

!!! definition "定义 68B.9 (SO(3) 上的误差状态 EKF)"
    状态为旋转 $R \in SO(3)$。定义误差状态 $\delta\phi \in \mathbb{R}^3$：
    $$
    R = \hat{R} \exp([\delta\phi]_\times)
    $$

    其中 $\hat{R}$ 是当前估计，$\delta\phi$ 是 Lie 代数中的小扰动。

!!! theorem "定理 68B.12 ($SO(3)$ 上的 EKF 算法)"
    **系统模型**：$R_{k+1} = R_k \exp([\omega_k \Delta t]_\times) \exp([\mathbf{w}_k]_\times)$（$\mathbf{w}_k$ 是过程噪声）。

    **观测模型**：例如加速度计测量重力方向 $\mathbf{a}_k = R_k^T \mathbf{g} + \mathbf{v}_k$。

    **预测步**：
    $$
    \hat{R}_{k+1|k} = \hat{R}_{k|k} \exp([\omega_k \Delta t]_\times)
    $$
    $$
    P_{k+1|k} = F_k P_{k|k} F_k^T + Q_k
    $$
    其中 $F_k \approx \exp(-[\omega_k\Delta t]_\times)$（误差状态转移矩阵），$Q_k$ 是过程噪声协方差。

    **更新步**：
    $$
    \mathbf{z}_k = \mathbf{a}_k - \hat{R}_{k+1|k}^T \mathbf{g} \quad (\text{观测残差})
    $$
    $$
    H_k = [\hat{R}_{k+1|k}^T \mathbf{g}]_\times \quad (\text{观测 Jacobi 矩阵})
    $$
    $$
    K_k = P_{k+1|k} H_k^T (H_k P_{k+1|k} H_k^T + V_k)^{-1} \quad (\text{Kalman 增益})
    $$
    $$
    \delta\phi = K_k \mathbf{z}_k
    $$
    $$
    \hat{R}_{k+1|k+1} = \hat{R}_{k+1|k} \exp([\delta\phi]_\times) \quad (\text{状态更新——在流形上})
    $$
    $$
    P_{k+1|k+1} = (I - K_k H_k) P_{k+1|k}
    $$

### $SE(3)$ 上的 EKF

!!! definition "定义 68B.10 ($SE(3)$ 上的误差状态 EKF)"
    状态为位姿 $T \in SE(3)$。误差状态 $\delta\xi \in \mathbb{R}^6$：
    $$
    T = \hat{T} \exp([\delta\xi]^\wedge)
    $$

    算法结构与 $SO(3)$ 类似，但误差状态为 6 维，Jacobi 矩阵使用 $\operatorname{Ad}$ 和 $\operatorname{ad}$ 映射。

!!! theorem "定理 68B.13 ($SE(3)$ EKF 的预测步)"
    若运动模型为 $T_{k+1} = T_k \exp([\xi_k \Delta t]^\wedge) \exp([\mathbf{w}_k]^\wedge)$，则：
    $$
    \hat{T}_{k+1|k} = \hat{T}_{k|k} \exp([\xi_k \Delta t]^\wedge)
    $$
    $$
    P_{k+1|k} = F_k P_{k|k} F_k^T + \operatorname{Ad}(\exp(-[\xi_k\Delta t]^\wedge)) Q_k \operatorname{Ad}(\exp(-[\xi_k\Delta t]^\wedge))^T
    $$
    其中 $F_k = \operatorname{Ad}(\exp(-[\xi_k\Delta t]^\wedge))$。

!!! note "注"
    流形 EKF 的关键创新是：

    - **状态本身在流形上**（$SO(3)$ 或 $SE(3)$），保证了旋转矩阵始终是正交的；
    - **误差状态在切空间（Lie 代数）中**，是线性空间，可以应用标准 Kalman 滤波；
    - **状态更新通过指数映射回到流形**，保持了群结构。

    这种方法在惯性导航系统（INS）、视觉惯性里程计（VIO）和卫星姿态确定中是标准做法。

---

## 68B.10 应用实例

### 视觉惯性里程计

!!! example "例 68B.5"
    **视觉惯性里程计**（Visual-Inertial Odometry, VIO）是将 IMU（惯性测量单元）和相机数据融合估计机器人位姿的系统。

    - **IMU 预测**：陀螺仪给出角速度 $\omega$，加速度计给出线加速度 $\mathbf{a}$。预测步在 $SO(3) \times \mathbb{R}^3$（或 $SE(3)$）上进行积分。
    - **视觉更新**：从图像中提取特征点，通过三角化和重投影误差提供位姿和路标点的观测。
    - **状态估计**：使用流形 EKF 或优化方法（滑动窗口优化），核心是在 $SE(3)$ 上的 Gauss-Newton 迭代。

### 机械臂力矩前馈

!!! example "例 68B.6"
    **工业机器人轨迹跟踪**：给定期望关节轨迹 $\theta_d(t)$，使用计算力矩控制：

    1. 离线规划轨迹 $\theta_d(t), \dot{\theta}_d(t), \ddot{\theta}_d(t)$；
    2. 实时计算前馈力矩 $\tau_{ff} = M(\theta)\ddot{\theta}_d + C(\theta, \dot{\theta})\dot{\theta}_d + G(\theta)$（递归 Newton-Euler）；
    3. 叠加反馈修正 $\tau_{fb} = K_P(\theta_d - \theta) + K_D(\dot{\theta}_d - \dot{\theta})$；
    4. 总控制力矩 $\tau = \tau_{ff} + \tau_{fb}$。

    前馈补偿了非线性动力学，反馈修正模型误差和扰动。

---

## 本章小结

本章展示了矩阵理论在机器人动力学和状态估计中的深层应用。主要内容包括：

1. **Euler-Lagrange 动力学** $M\ddot{\theta} + C\dot{\theta} + G = \tau$ 具有丰富的矩阵结构：质量矩阵 $M$ 正定对称，$\dot{M} - 2C$ 反对称（被动性）。

2. **递归 Newton-Euler 算法**以 $O(n)$ 复杂度计算逆动力学，是实时控制的基石。

3. **计算力矩控制**利用精确的动力学模型实现反馈线性化，使跟踪误差动力学成为线性。

4. **BCH 公式**描述了 Lie 代数元素组合的规律，$e^X e^Y = e^{X + Y + \frac{1}{2}[X,Y] + \cdots}$，解释了旋转的不可交换性。

5. **对偶四元数**提供了比 $4 \times 4$ 矩阵更紧凑的刚体运动表示，支持自然的插值。

6. **旋转平均**通过 $SO(3)$ 上的 Riemannian 梯度下降计算 Karcher 均值。

7. **$SE(3)$ 协方差传播**利用伴随表示处理不确定性在坐标系变换下的传播。

8. **位姿图优化**（SLAM）在 $SE(3)^N$ 上做非线性最小二乘，每次迭代归结为大型稀疏线性方程组。

9. **流形 EKF** 在 Lie 代数中定义误差状态，通过指数映射更新流形上的状态，广泛用于 IMU/VIO 系统。

---

## 习题

!!! question "习题 68B.1"
    证明质量矩阵 $M(\theta)$ 是正定的。（提示：利用动能的物理意义 $T = \frac{1}{2}\dot{\theta}^T M(\theta)\dot{\theta} > 0$ 对非零 $\dot{\theta}$。）

!!! question "习题 68B.2"
    对平面两连杆机器人（例 68B.1），验证 $\dot{M} - 2C$ 是反对称矩阵。

!!! question "习题 68B.3"
    设两个位姿测量 $T_1, T_2 \in SE(3)$ 有协方差 $\Sigma_1, \Sigma_2 \in \mathbb{R}^{6 \times 6}$。利用协方差传播公式计算 $T_1 T_2$ 的协方差。

!!! question "习题 68B.4"
    对 $\mathfrak{so}(3)$，设 $X = \theta_1[\hat{e}_1]_\times$，$Y = \theta_2[\hat{e}_2]_\times$。计算 BCH 公式的前三项 $Z = X + Y + \frac{1}{2}[X,Y] + \cdots$。当 $\theta_1 = \theta_2 = 0.1$ rad 时，比较 $e^Z$ 与 $e^X e^Y$ 的差异。

!!! question "习题 68B.5"
    写出三连杆平面机器人的质量矩阵 $M(\theta)$，并验证其正定性。

!!! question "习题 68B.6"
    对计算力矩控制，设 $K_P = \operatorname{diag}(100, 100)$，$K_D = \operatorname{diag}(20, 20)$。分析误差动力学 $\ddot{e} + K_D\dot{e} + K_Pe = 0$ 的特征值和阻尼比。

!!! question "习题 68B.7"
    推导递归 Newton-Euler 算法中前向递推的角加速度公式 $\dot{\omega}_i$。

!!! question "习题 68B.8"
    设 $N = 4$ 个旋转 $R_i$，分别为绕 $z$ 轴旋转 $0°, 30°, 60°, 90°$。用旋转平均迭代算法计算 Karcher 均值。初始值取 $\bar{R}_0 = R_1$，迭代至收敛。

!!! question "习题 68B.9"
    对偶四元数乘法：设 $\hat{q}_1$ 表示绕 $z$ 轴旋转 $90°$ 后平移 $(1, 0, 0)$，$\hat{q}_2$ 表示绕 $x$ 轴旋转 $45°$。计算 $\hat{q}_1\hat{q}_2$ 并解释其几何含义。

!!! question "习题 68B.10"
    考虑一个简单的 SLAM 问题：3 个位姿，3 条边（形成三角形）。每条边的测量为 $\tilde{T}_{12}, \tilde{T}_{23}, \tilde{T}_{13}$，协方差均为 $\sigma^2 I_6$。

    (a) 写出位姿图优化的目标函数。

    (b) 固定 $T_1 = I$（消除全局自由度），对 $T_2, T_3$ 写出 Gauss-Newton 方程的结构。

    (c) 讨论稀疏性：Hessian 矩阵 $H$ 的非零块模式是什么？

!!! question "习题 68B.11"
    设陀螺仪测量角速度 $\omega_k = \bar{\omega}_k + \mathbf{b}_k + \mathbf{w}_k$，其中 $\bar{\omega}_k$ 是真实角速度，$\mathbf{b}_k$ 是缓变偏置（bias），$\mathbf{w}_k \sim \mathcal{N}(0, \sigma_g^2 I)$ 是噪声。设计一个 $SO(3) \times \mathbb{R}^3$ 上的 EKF，状态为 $(R, \mathbf{b})$，估计旋转和陀螺偏置。

!!! question "习题 68B.12"
    证明 $\dot{M} - 2C$ 的反对称性蕴含被动性：沿运动轨线，
    $$
    \frac{d}{dt}\left(\frac{1}{2}\dot{\theta}^T M\dot{\theta}\right) = \dot{\theta}^T(\tau - G)
    $$
    并解释其物理含义（能量守恒的推广）。

!!! question "习题 68B.13"
    比较 $SE(3)$ 上的两种误差定义："左误差" $\epsilon_L = \log(\hat{T}^{-1}T)^\vee$ 和"右误差" $\epsilon_R = \log(T\hat{T}^{-1})^\vee$。它们之间有什么关系？在位姿图优化和 EKF 中分别使用哪种更合适？

!!! question "习题 68B.14"
    （综合题）设一个移动机器人在平面上运动（$SE(2)$），状态为 $(x, y, \theta)$。运动模型为差速驱动（给定左右轮速度），观测为 GPS 位置 $(x_{gps}, y_{gps})$。

    (a) 将运动模型写成 $SE(2)$ 上的递推形式。

    (b) 设计误差状态 EKF，写出预测和更新步的公式。

    (c) 讨论与标准 EKF（直接用 $(x, y, \theta)$ 作为状态）的区别。
