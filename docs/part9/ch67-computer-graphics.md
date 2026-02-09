# 第 67 章 线性代数在计算机图形学中的应用

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 线性变换(Ch5) · 正交矩阵(Ch7) · 特征值(Ch6)

**本章脉络**：齐次坐标 → 仿射变换 → 旋转表示(矩阵/四元数/Euler角/SLERP) → 模型-视图-投影矩阵 → 透视与正交投影 → 视口变换 → 透视校正插值 → 反射矩阵 → 法向量变换 → 颜色空间变换 → 样条曲面

**延伸**：计算机图形学是线性代数最直观的应用领域之一；GPU 硬件本质上是大规模矩阵乘法加速器；现代实时渲染管线中每一帧都执行数百万次 4×4 矩阵变换

</div>

计算机图形学是线性代数最直观、最生动的应用领域。从三维物体的建模到屏幕上像素的最终着色，整个渲染管线中的每一个步骤都依赖矩阵变换。平移、旋转、缩放、投影——这些几何操作都可以统一表示为 $4 \times 4$ 矩阵乘法。现代 GPU（图形处理单元）的硬件架构本质上就是为了高效执行大规模并行矩阵运算而设计的。

本章从齐次坐标的引入开始，系统展示线性代数在计算机图形学中的核心角色。我们将看到，通过引入一个额外的维度，仿射变换（包括平移）都可以用线性变换来表示，从而所有几何变换都统一为矩阵乘法。

---

## 67.1 齐次坐标

<div class="context-flow" markdown>

**核心问题**：如何将平移（非线性操作）纳入线性变换的统一框架？

</div>

### 动机

在标准的三维坐标中，平移操作 $\mathbf{p}' = \mathbf{p} + \mathbf{t}$ 不能表示为矩阵乘法 $\mathbf{p}' = M\mathbf{p}$（因为平移不是线性变换——它不保持原点不动）。然而，通过引入**齐次坐标**，我们可以将平移也表示为矩阵乘法。

!!! definition "定义 67.1 (齐次坐标)"
    三维空间中的**点** $\mathbf{p} = (x, y, z)$ 用四维齐次坐标表示为
    $$
    \tilde{\mathbf{p}} = \begin{pmatrix} x \\ y \\ z \\ 1 \end{pmatrix}
    $$

    更一般地，齐次坐标 $(x, y, z, w)$（$w \neq 0$）表示三维点 $(x/w, y/w, z/w)$。

    三维**方向向量** $\mathbf{v} = (v_x, v_y, v_z)$ 用齐次坐标表示为
    $$
    \tilde{\mathbf{v}} = \begin{pmatrix} v_x \\ v_y \\ v_z \\ 0 \end{pmatrix}
    $$

!!! note "注"
    点和方向的区别在第四个分量：点的 $w = 1$，方向的 $w = 0$。这反映了一个深刻的几何事实：

    - 两个点相减得到方向向量：$(x_1, y_1, z_1, 1) - (x_2, y_2, z_2, 1) = (x_1-x_2, y_1-y_2, z_1-z_2, 0)$；
    - 点加方向向量得到点：$(x, y, z, 1) + (v_x, v_y, v_z, 0) = (x+v_x, y+v_y, z+v_z, 1)$；
    - 两个方向相加得到方向：$w = 0 + 0 = 0$。

### 为什么用四维表示三维？

!!! theorem "定理 67.1 (齐次坐标的统一性)"
    在齐次坐标下，所有仿射变换（平移、旋转、缩放、剪切及其组合）都可以表示为 $4 \times 4$ 矩阵乘法。

    仿射变换 $\mathbf{p}' = M\mathbf{p} + \mathbf{t}$（其中 $M$ 是 $3 \times 3$ 矩阵，$\mathbf{t}$ 是平移向量）在齐次坐标下写为
    $$
    \begin{pmatrix} \mathbf{p}' \\ 1 \end{pmatrix} = \begin{pmatrix} M & \mathbf{t} \\ \mathbf{0}^T & 1 \end{pmatrix} \begin{pmatrix} \mathbf{p} \\ 1 \end{pmatrix}
    $$

!!! example "例 67.1"
    平移 $(x, y, z) \mapsto (x + t_x, y + t_y, z + t_z)$ 对应矩阵
    $$
    T(t_x, t_y, t_z) = \begin{pmatrix} 1 & 0 & 0 & t_x \\ 0 & 1 & 0 & t_y \\ 0 & 0 & 1 & t_z \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

    验证：$T \cdot (x, y, z, 1)^T = (x+t_x, y+t_y, z+t_z, 1)^T$。

    注意 $T \cdot (v_x, v_y, v_z, 0)^T = (v_x, v_y, v_z, 0)^T$——方向向量不受平移影响，符合几何直觉。

---

## 67.2 基本几何变换

<div class="context-flow" markdown>

**核心问题**：如何用 $4 \times 4$ 矩阵表示所有基本几何变换？

</div>

### 缩放

!!! definition "定义 67.2 (缩放矩阵)"
    沿坐标轴的缩放变换 $(x, y, z) \mapsto (s_x x, s_y y, s_z z)$ 对应矩阵
    $$
    S(s_x, s_y, s_z) = \begin{pmatrix} s_x & 0 & 0 & 0 \\ 0 & s_y & 0 & 0 \\ 0 & 0 & s_z & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

    均匀缩放（$s_x = s_y = s_z = s$）为 $S(s, s, s)$。

### 旋转

!!! definition "定义 67.3 (坐标轴旋转矩阵)"
    绕 $z$ 轴旋转角 $\theta$：
    $$
    R_z(\theta) = \begin{pmatrix} \cos\theta & -\sin\theta & 0 & 0 \\ \sin\theta & \cos\theta & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

    绕 $x$ 轴旋转角 $\theta$：
    $$
    R_x(\theta) = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos\theta & -\sin\theta & 0 \\ 0 & \sin\theta & \cos\theta & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

    绕 $y$ 轴旋转角 $\theta$：
    $$
    R_y(\theta) = \begin{pmatrix} \cos\theta & 0 & \sin\theta & 0 \\ 0 & 1 & 0 & 0 \\ -\sin\theta & 0 & \cos\theta & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

### 剪切

!!! definition "定义 67.4 (剪切矩阵)"
    $xy$ 剪切（$x$ 方向随 $y$ 变化）：
    $$
    H_{xy}(s) = \begin{pmatrix} 1 & s & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

### 变换的组合

!!! theorem "定理 67.2 (变换组合的非交换性)"
    几何变换的组合对应矩阵的乘法。但矩阵乘法**不满足交换律**，因此变换的顺序至关重要。

    若先应用变换 $M_1$，再应用 $M_2$，则组合变换为 $M = M_2 M_1$（注意顺序：后应用的矩阵在左边）。

!!! example "例 67.2"
    先绕原点旋转 $90°$（$R_z(\pi/2)$），再平移 $(1, 0, 0)$（$T(1,0,0)$）：
    $$
    M_1 = T \cdot R_z = \begin{pmatrix} 0 & -1 & 0 & 1 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

    先平移 $(1, 0, 0)$，再绕原点旋转 $90°$：
    $$
    M_2 = R_z \cdot T = \begin{pmatrix} 0 & -1 & 0 & 0 \\ 1 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

    $M_1 \neq M_2$：顺序不同，结果不同。

!!! example "例 67.3"
    **绕任意点 $\mathbf{p}$ 旋转**：先平移使 $\mathbf{p}$ 到原点，旋转，再平移回去：
    $$
    M = T(\mathbf{p}) \cdot R \cdot T(-\mathbf{p})
    $$

---

## 67.3 旋转表示

<div class="context-flow" markdown>

**核心问题**：如何高效且无奇异地表示三维旋转？

</div>

### Rodrigues 旋转公式

!!! theorem "定理 67.3 (Rodrigues 公式)"
    绕单位轴 $\hat{\mathbf{k}} = (k_x, k_y, k_z)$ 旋转角 $\theta$ 的旋转矩阵为
    $$
    R = I + \sin\theta \, K + (1 - \cos\theta) K^2
    $$
    其中 $K$ 是轴向量的反对称矩阵（叉积矩阵）：
    $$
    K = [\hat{\mathbf{k}}]_\times = \begin{pmatrix} 0 & -k_z & k_y \\ k_z & 0 & -k_x \\ -k_y & k_x & 0 \end{pmatrix}
    $$

??? proof "证明"
    将任意向量 $\mathbf{v}$ 分解为平行于 $\hat{\mathbf{k}}$ 和垂直于 $\hat{\mathbf{k}}$ 的分量：
    $$
    \mathbf{v}_\parallel = (\hat{\mathbf{k}} \cdot \mathbf{v})\hat{\mathbf{k}} = (\hat{\mathbf{k}}\hat{\mathbf{k}}^T)\mathbf{v}, \quad \mathbf{v}_\perp = \mathbf{v} - \mathbf{v}_\parallel
    $$

    旋转后 $\mathbf{v}_\parallel$ 不变，$\mathbf{v}_\perp$ 在其所在平面中旋转角 $\theta$：
    $$
    R\mathbf{v} = \mathbf{v}_\parallel + \cos\theta \, \mathbf{v}_\perp + \sin\theta \, (\hat{\mathbf{k}} \times \mathbf{v}_\perp)
    $$

    利用 $K\mathbf{v} = \hat{\mathbf{k}} \times \mathbf{v}$，$K^2\mathbf{v} = \hat{\mathbf{k}} \times (\hat{\mathbf{k}} \times \mathbf{v}) = (\hat{\mathbf{k}}\hat{\mathbf{k}}^T - I)\mathbf{v} = \mathbf{v}_\parallel - \mathbf{v}$：
    $$
    R\mathbf{v} = \mathbf{v}_\parallel + \cos\theta(\mathbf{v} - \mathbf{v}_\parallel) + \sin\theta \, K\mathbf{v}
    $$
    $$
    = \mathbf{v} + (1-\cos\theta)(\mathbf{v}_\parallel - \mathbf{v}) + \sin\theta \, K\mathbf{v}
    $$
    $$
    = \mathbf{v} + (1-\cos\theta) K^2\mathbf{v} + \sin\theta \, K\mathbf{v}
    $$
    即 $R = I + \sin\theta \, K + (1 - \cos\theta)K^2$。

### Euler 角与万向锁

!!! definition "定义 67.5 (Euler 角)"
    **Euler 角**表示法将任意旋转分解为三次绕坐标轴的旋转：
    $$
    R = R_z(\psi) R_y(\theta) R_x(\phi)
    $$
    其中 $\phi$ 为横滚（roll），$\theta$ 为俯仰（pitch），$\psi$ 为偏航（yaw）。

    共有 12 种不同的 Euler 角约定（取决于旋转轴的选择和顺序）。

!!! definition "定义 67.6 (万向锁)"
    **万向锁**（gimbal lock）发生在中间轴的旋转角达到 $\pm 90°$ 时。例如当 $\theta = \pm\pi/2$ 时，第一次和第三次旋转的效果退化为绕同一轴，失去一个自由度。

    数学上，当 $\theta = \pi/2$ 时：
    $$
    R_z(\psi)R_y(\pi/2)R_x(\phi) = \begin{pmatrix} 0 & \sin(\phi-\psi) & \cos(\phi-\psi) \\ 0 & \cos(\phi-\psi) & -\sin(\phi-\psi) \\ -1 & 0 & 0 \end{pmatrix}
    $$
    只依赖 $\phi - \psi$，因此 $\phi$ 和 $\psi$ 不再独立。

### 四元数旋转

!!! definition "定义 67.7 (单位四元数与旋转)"
    **单位四元数** $\mathbf{q} = w + xi + yj + zk$（$w^2 + x^2 + y^2 + z^2 = 1$）表示绕轴 $\hat{\mathbf{k}} = (x, y, z)/\|(x,y,z)\|$ 旋转角 $\theta$ 的旋转，其中
    $$
    \mathbf{q} = \cos\frac{\theta}{2} + \sin\frac{\theta}{2}(k_x i + k_y j + k_z k)
    $$

    点 $\mathbf{p} = (p_x, p_y, p_z)$ 经旋转后的结果为
    $$
    \mathbf{p}' = \mathbf{q} \, \mathbf{p} \, \mathbf{q}^{-1}
    $$
    其中 $\mathbf{p}$ 被视为纯四元数 $p_x i + p_y j + p_z k$。

!!! theorem "定理 67.4 (四元数到旋转矩阵)"
    单位四元数 $\mathbf{q} = (w, x, y, z)$ 对应的旋转矩阵为
    $$
    R = \begin{pmatrix} 1 - 2(y^2+z^2) & 2(xy-wz) & 2(xz+wy) \\ 2(xy+wz) & 1-2(x^2+z^2) & 2(yz-wx) \\ 2(xz-wy) & 2(yz+wx) & 1-2(x^2+y^2) \end{pmatrix}
    $$

### 旋转表示的比较

!!! note "注"
    三种旋转表示方法的比较：

    | 特性 | 旋转矩阵 | Euler 角 | 四元数 |
    |------|---------|---------|--------|
    | 参数数目 | 9（6个约束） | 3 | 4（1个约束） |
    | 万向锁 | 无 | 有 | 无 |
    | 插值 | 困难 | 困难 | SLERP 插值 |
    | 组合 | 矩阵乘法 | 非平凡 | 四元数乘法 |
    | 内存 | 9 floats | 3 floats | 4 floats |
    | 应用旋转 | 矩阵-向量乘 | 需转换 | 四元数乘法 |

### SLERP 球面线性插值

在动画和运动控制中，需要在两个旋转之间进行平滑插值。四元数的 **SLERP**（Spherical Linear Interpolation，球面线性插值）是解决这一问题的标准方法。

!!! definition "定义 67.17 (SLERP)"
    给定两个单位四元数 $q_0$ 和 $q_1$，它们之间的**球面线性插值**定义为

    $$\mathrm{Slerp}(q_0, q_1; t) = q_0 \frac{\sin((1-t)\theta)}{\sin\theta} + q_1 \frac{\sin(t\theta)}{\sin\theta}, \quad t \in [0, 1],$$

    其中 $\theta = \arccos(q_0 \cdot q_1)$ 是两个四元数在 $S^3$ 上的夹角（$q_0 \cdot q_1 = w_0 w_1 + x_0 x_1 + y_0 y_1 + z_0 z_1$）。

    当 $t = 0$ 时得到 $q_0$，$t = 1$ 时得到 $q_1$。

??? proof "推导"
    **从 $S^3$ 上的测地线推导。** 单位四元数构成三维球面 $S^3$。$S^3$ 上连接 $q_0$ 和 $q_1$ 的最短路径（测地线）是大圆弧。

    在二维情形下类比：$S^1$ 上从 $p_0$ 到 $p_1$ 的大圆弧参数化为 $p(t) = \alpha(t) p_0 + \beta(t) p_1$，要求 $|p(t)| = 1$ 且 $p(0) = p_0$，$p(1) = p_1$。

    设 $\cos\theta = q_0 \cdot q_1$。寻找系数 $\alpha(t), \beta(t)$ 使得 $q(t) = \alpha(t) q_0 + \beta(t) q_1$ 满足：
    - $q(0) = q_0$：$\alpha(0) = 1, \beta(0) = 0$。
    - $q(1) = q_1$：$\alpha(1) = 0, \beta(1) = 1$。
    - $|q(t)| = 1$ 且参数 $t$ 是角度的线性函数（等角速度）。

    从测地线方程 $q(t) = q_0 \cos(t\theta) + \frac{q_1 - q_0\cos\theta}{\sin\theta}\sin(t\theta)$ 出发，整理得：

    $$q(t) = q_0\frac{\cos(t\theta)\sin\theta - \sin(t\theta)\cos\theta + \sin(t\theta)\cos\theta}{\sin\theta} + q_1\frac{\sin(t\theta)}{\sin\theta}$$

    $$= q_0\frac{\sin((1-t)\theta)}{\sin\theta} + q_1\frac{\sin(t\theta)}{\sin\theta}. \quad \blacksquare$$

!!! theorem "定理 67.11 (SLERP 的性质)"
    **(a)** **等角速度**：$\mathrm{Slerp}(q_0, q_1; t)$ 沿 $S^3$ 大圆弧以恒定角速度运动。具体地，$q(t)$ 与 $q(t + \Delta t)$ 之间的角度为 $|\Delta t| \cdot \theta$，与 $t$ 无关。

    **(b)** **最短路径**：SLERP 给出 $S^3$ 上连接 $q_0$ 和 $q_1$ 的最短测地线（前提是 $\theta \leq \pi$）。

    **(c)** **旋转不变性**：$\mathrm{Slerp}(pq_0, pq_1; t) = p \cdot \mathrm{Slerp}(q_0, q_1; t)$，即 SLERP 与左乘旋转交换。

!!! note "注"
    **与 NLERP 的比较。** 归一化线性插值（Normalized Linear Interpolation）定义为

    $$\mathrm{Nlerp}(q_0, q_1; t) = \frac{(1-t)q_0 + tq_1}{\|(1-t)q_0 + tq_1\|}.$$

    NLERP 计算简单（无三角函数），但不具有等角速度：在 $t = 0$ 和 $t = 1$ 附近速度慢，$t = 0.5$ 附近速度快。当 $\theta$ 较小时（$\theta < 20°$），NLERP 近似于 SLERP；$\theta$ 较大时速度不均匀性明显。游戏引擎中常对小角度用 NLERP、大角度用 SLERP 以平衡性能和质量。

!!! note "注"
    **实现注意事项：双重覆盖与反足问题。** 单位四元数到旋转的映射是 2-1 的：$q$ 和 $-q$ 表示同一旋转（$SU(2) \to SO(3)$ 的核为 $\{\pm 1\}$）。因此在进行 SLERP 之前，必须检查 $q_0 \cdot q_1$ 的符号：

    - 若 $q_0 \cdot q_1 < 0$，将 $q_1$ 替换为 $-q_1$（表示同一旋转，但选择短弧而非长弧）。
    - 这保证 $\theta \leq \pi/2$，SLERP 走最短路径。

    当 $\theta \approx 0$ 时（$\sin\theta \approx 0$），SLERP 公式有数值不稳定性。此时应退化为线性插值 $(1-t)q_0 + tq_1$（无需归一化，因为已经几乎是单位四元数）。

---

## 67.4 模型-视图-投影管线

<div class="context-flow" markdown>

**核心问题**：如何将三维场景中的物体变换到二维屏幕上？

</div>

### MVP 管线

!!! definition "定义 67.8 (MVP 管线)"
    图形渲染管线中，三维顶点 $\mathbf{v}$ 经过三个矩阵变换到达屏幕空间：
    $$
    \mathbf{v}_{clip} = P \cdot V \cdot M \cdot \mathbf{v}_{local}
    $$

    其中：

    - **M**（Model 矩阵）：将物体从**局部坐标系**（object space）变换到**世界坐标系**（world space）；
    - **V**（View 矩阵）：将世界坐标系变换到**相机坐标系**（camera/view space）；
    - **P**（Projection 矩阵）：将三维相机坐标变换到**裁剪坐标**（clip space），实现投影。

### 模型矩阵

!!! definition "定义 67.9 (模型矩阵)"
    模型矩阵 $M$ 通常是旋转、缩放和平移的组合：
    $$
    M = T \cdot R \cdot S
    $$
    将物体放置在世界中的正确位置、方向和大小。

### 视图矩阵

!!! theorem "定理 67.5 (LookAt 视图矩阵)"
    给定相机位置 $\mathbf{eye}$、目标点 $\mathbf{center}$ 和上方向 $\mathbf{up}$，视图矩阵为
    $$
    V = \begin{pmatrix} r_x & r_y & r_z & -\mathbf{r} \cdot \mathbf{eye} \\ u_x & u_y & u_z & -\mathbf{u} \cdot \mathbf{eye} \\ -f_x & -f_y & -f_z & \mathbf{f} \cdot \mathbf{eye} \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

    其中
    $$
    \mathbf{f} = \frac{\mathbf{center} - \mathbf{eye}}{\|\mathbf{center} - \mathbf{eye}\|}, \quad
    \mathbf{r} = \frac{\mathbf{f} \times \mathbf{up}}{\|\mathbf{f} \times \mathbf{up}\|}, \quad
    \mathbf{u} = \mathbf{r} \times \mathbf{f}
    $$

!!! note "注"
    视图矩阵本质上是一个刚体变换的逆：它将世界坐标系"移动"到与相机坐标系对齐。$V = (T_{eye} \cdot R_{cam})^{-1} = R_{cam}^T \cdot T_{-eye}$。

---

## 67.5 透视与正交投影

<div class="context-flow" markdown>

**核心问题**：如何将三维场景投影到二维平面上？

</div>

### 透视投影

!!! definition "定义 67.10 (透视投影矩阵)"
    给定视场角（field of view）$\text{fov}$、宽高比 $\text{aspect} = w/h$、近裁剪面 $n$ 和远裁剪面 $f$，**透视投影矩阵**为
    $$
    P_{persp} = \begin{pmatrix} \frac{1}{\text{aspect} \cdot \tan(\text{fov}/2)} & 0 & 0 & 0 \\ 0 & \frac{1}{\tan(\text{fov}/2)} & 0 & 0 \\ 0 & 0 & \frac{-(f+n)}{f-n} & \frac{-2fn}{f-n} \\ 0 & 0 & -1 & 0 \end{pmatrix}
    $$

!!! note "注"
    透视投影的关键在于最后一行 $(0, 0, -1, 0)$，它使得 $w_{clip} = -z_{eye}$（注意 $z_{eye} < 0$ 在相机前方）。齐次除法 $(x_{clip}/w_{clip}, y_{clip}/w_{clip}, z_{clip}/w_{clip})$ 自动实现了**近大远小**的透视效果：距离相机越远的物体（$|z|$ 越大），除以越大的 $w$，因而在屏幕上显得越小。

!!! example "例 67.4"
    设 $\text{fov} = 90°$，$\text{aspect} = 16/9$，$n = 0.1$，$f = 100$。则 $\tan(45°) = 1$，
    $$
    P = \begin{pmatrix} 9/16 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & -\frac{100.1}{99.9} & -\frac{20}{99.9} \\ 0 & 0 & -1 & 0 \end{pmatrix}
    $$

    点 $(1, 0, -5, 1)$（相机前方 5 个单位处）变换为：
    $$
    P \cdot (1, 0, -5, 1)^T = (9/16, 0, 5.0 \cdot \frac{100.1}{99.9} - \frac{20}{99.9}, 5)^T
    $$
    齐次除法后 $x_{ndc} = 9/(16 \cdot 5) = 9/80 \approx 0.1125$。

### 正交投影

!!! definition "定义 67.11 (正交投影矩阵)"
    正交投影将视景体 $[l, r] \times [b, t] \times [n, f]$ 映射到标准立方体 $[-1, 1]^3$：
    $$
    P_{ortho} = \begin{pmatrix} \frac{2}{r-l} & 0 & 0 & -\frac{r+l}{r-l} \\ 0 & \frac{2}{t-b} & 0 & -\frac{t+b}{t-b} \\ 0 & 0 & \frac{-2}{f-n} & -\frac{f+n}{f-n} \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

!!! note "注"
    正交投影保持 $w = 1$（最后一行是 $(0,0,0,1)$），因此没有透视缩小效果。正交投影用于 CAD、建筑制图和 2D 游戏。

### 深度缓冲

!!! theorem "定理 67.6 (深度值的非线性性)"
    透视投影后的深度值 $z_{ndc} = z_{clip}/w_{clip}$ 关于原始深度 $z_{eye}$ 是**非线性**的：
    $$
    z_{ndc} = \frac{f + n}{f - n} + \frac{2fn}{(f - n) z_{eye}}
    $$

    这导致近处的深度精度远高于远处。深度缓冲的精度不均匀是图形学中"z-fighting"伪影的根本原因。

### 视口变换

MVP 管线将顶点变换到裁剪空间（clip space），经齐次除法后进入标准化设备坐标（NDC，范围 $[-1, 1]^3$）。但最终需要将 NDC 映射到屏幕像素坐标——这就是**视口变换**。

!!! definition "定义 67.18 (视口变换)"
    设屏幕窗口的左下角为 $(x_0, y_0)$，宽度为 $w$，高度为 $h$，深度范围为 $[d_{near}, d_{far}]$（通常 $[0, 1]$）。视口变换将 NDC 坐标 $(x_{ndc}, y_{ndc}, z_{ndc}) \in [-1,1]^3$ 映射到屏幕坐标 $(x_s, y_s, z_s)$：

    $$\begin{pmatrix} x_s \\ y_s \\ z_s \end{pmatrix} = \begin{pmatrix} \frac{w}{2} x_{ndc} + x_0 + \frac{w}{2} \\ \frac{h}{2} y_{ndc} + y_0 + \frac{h}{2} \\ \frac{d_{far} - d_{near}}{2} z_{ndc} + \frac{d_{far} + d_{near}}{2} \end{pmatrix}$$

    用矩阵形式表示（齐次坐标）：

    $$V_{port} = \begin{pmatrix} w/2 & 0 & 0 & x_0 + w/2 \\ 0 & h/2 & 0 & y_0 + h/2 \\ 0 & 0 & \frac{d_{far} - d_{near}}{2} & \frac{d_{far} + d_{near}}{2} \\ 0 & 0 & 0 & 1 \end{pmatrix}$$

!!! note "注"
    完整的顶点变换管线为：

    $$\mathbf{v}_{local} \xrightarrow{M} \mathbf{v}_{world} \xrightarrow{V} \mathbf{v}_{eye} \xrightarrow{P} \mathbf{v}_{clip} \xrightarrow{\div w} \mathbf{v}_{ndc} \xrightarrow{V_{port}} \mathbf{v}_{screen}$$

    其中齐次除法 $\div w$ 是唯一的非线性步骤。所有其他步骤都是矩阵乘法。

### 透视校正插值

!!! definition "定义 67.19 (重心坐标)"
    屏幕空间中三角形的三个顶点 $\mathbf{p}_0, \mathbf{p}_1, \mathbf{p}_2$ 内部的点 $\mathbf{p}$ 可以用**重心坐标** $(\lambda_0, \lambda_1, \lambda_2)$ 表示：

    $$\mathbf{p} = \lambda_0 \mathbf{p}_0 + \lambda_1 \mathbf{p}_1 + \lambda_2 \mathbf{p}_2, \quad \lambda_0 + \lambda_1 + \lambda_2 = 1, \quad \lambda_i \geq 0.$$

!!! theorem "定理 67.12 (透视校正插值)"
    透视投影后，在屏幕空间三角形内对顶点属性（如纹理坐标、颜色）进行插值时，**不能**简单地用屏幕空间重心坐标做线性插值——这会导致失真（因为透视除法是非线性的）。

    正确的**透视校正插值**公式：设三角形顶点的齐次坐标的 $w$ 分量分别为 $w_0, w_1, w_2$（即 clip space 的 $w = -z_{eye}$），顶点属性值为 $a_0, a_1, a_2$，屏幕空间重心坐标为 $(\lambda_0, \lambda_1, \lambda_2)$。则像素处的正确属性值为

    $$a = \frac{\lambda_0 a_0 / w_0 + \lambda_1 a_1 / w_1 + \lambda_2 a_2 / w_2}{\lambda_0 / w_0 + \lambda_1 / w_1 + \lambda_2 / w_2}.$$

??? proof "证明"
    在 clip space 中，属性沿三角形线性变化：$a(\mathbf{v}) = \mu_0 a_0 + \mu_1 a_1 + \mu_2 a_2$，其中 $(\mu_0, \mu_1, \mu_2)$ 是 clip space 的重心坐标。

    齐次除法将 clip space 的点 $\tilde{\mathbf{v}}_i = (x_i, y_i, z_i, w_i)$ 映射到 NDC 的 $\mathbf{v}_i = (x_i/w_i, y_i/w_i, z_i/w_i)$。

    设屏幕空间的重心坐标为 $\lambda_i$，clip space 的重心坐标为 $\mu_i$。由于齐次除法的非线性，$\lambda_i \neq \mu_i$。它们的关系为：

    $$\mu_i = \frac{\lambda_i / w_i}{\sum_j \lambda_j / w_j}.$$

    将此代入 $a = \sum_i \mu_i a_i$ 即得透视校正公式。

    推导 $\mu_i$ 与 $\lambda_i$ 的关系：在 clip space 中点 $\tilde{\mathbf{v}} = \mu_0 \tilde{\mathbf{v}}_0 + \mu_1 \tilde{\mathbf{v}}_1 + \mu_2 \tilde{\mathbf{v}}_2$，其第四分量 $w = \mu_0 w_0 + \mu_1 w_1 + \mu_2 w_2$。齐次除法后 $\mathbf{v} = \tilde{\mathbf{v}} / w$。在屏幕空间中 $\mathbf{v} = \lambda_0 \mathbf{v}_0 + \lambda_1 \mathbf{v}_1 + \lambda_2 \mathbf{v}_2$，故 $\tilde{\mathbf{v}}/w = \lambda_0 \tilde{\mathbf{v}}_0/w_0 + \lambda_1 \tilde{\mathbf{v}}_1/w_1 + \lambda_2 \tilde{\mathbf{v}}_2/w_2$。与 $\tilde{\mathbf{v}} = w(\lambda_0 \tilde{\mathbf{v}}_0/w_0 + \cdots)$ 和 $\tilde{\mathbf{v}} = \sum \mu_i \tilde{\mathbf{v}}_i$ 对比，得 $\mu_i = w \lambda_i / w_i$。由 $\sum \mu_i = 1$ 得 $w = 1/\sum_j(\lambda_j / w_j)$。$\blacksquare$

!!! note "注"
    现代 GPU 硬件在光栅化阶段自动执行透视校正插值。顶点着色器输出的属性值在传递给片段着色器之前，由硬件按上述公式进行透视校正。程序员通常无需手动实现，但理解其数学原理对于调试渲染错误和实现自定义光栅化至关重要。

---

## 67.6 反射矩阵

<div class="context-flow" markdown>

**核心问题**：如何用矩阵表示反射变换？反射与旋转有什么关系？

</div>

反射是与旋转同样基本的几何变换，广泛应用于镜面效果、对称建模和法向量计算。

!!! definition "定义 67.20 (Householder 反射矩阵)"
    给定单位法向量 $\mathbf{n} \in \mathbb{R}^3$（$\|\mathbf{n}\| = 1$），关于通过原点、以 $\mathbf{n}$ 为法线的平面的**反射矩阵**（Householder 矩阵）为

    $$H = I - 2\mathbf{n}\mathbf{n}^T.$$

    对任意向量 $\mathbf{v}$：
    $$H\mathbf{v} = \mathbf{v} - 2(\mathbf{n} \cdot \mathbf{v})\mathbf{n},$$
    即 $\mathbf{v}$ 的平行于 $\mathbf{n}$ 的分量取反，垂直于 $\mathbf{n}$ 的分量不变。

!!! theorem "定理 67.13 (反射矩阵的性质)"
    **(a)** $H$ 是对称的：$H^T = H$。

    **(b)** $H$ 是正交的：$H^T H = I$。

    **(c)** $H$ 是对合的（自逆的）：$H^2 = I$，即 $H^{-1} = H$。

    **(d)** $\det H = -1$（反射是保距变换但不保持定向——它是"反向"正交变换）。

    **(e)** $H$ 的特征值为 $-1$（对应 $\mathbf{n}$，一重）和 $+1$（对应 $\mathbf{n}^\perp$ 平面，$n-1$ 重）。

??? proof "证明"
    **(b)** $H^TH = (I - 2\mathbf{n}\mathbf{n}^T)(I - 2\mathbf{n}\mathbf{n}^T) = I - 4\mathbf{n}\mathbf{n}^T + 4\mathbf{n}(\mathbf{n}^T\mathbf{n})\mathbf{n}^T = I - 4\mathbf{n}\mathbf{n}^T + 4\mathbf{n}\mathbf{n}^T = I$。

    **(d)** $\det H = \det(I - 2\mathbf{n}\mathbf{n}^T)$。由矩阵行列式引理 $\det(I + \mathbf{u}\mathbf{v}^T) = 1 + \mathbf{v}^T\mathbf{u}$，得 $\det H = 1 + (-2)\mathbf{n}^T\mathbf{n} = 1 - 2 = -1$。$\blacksquare$

!!! definition "定义 67.21 (一般平面反射)"
    关于通过点 $\mathbf{p}$、法向量为 $\mathbf{n}$ 的平面的反射，在齐次坐标下为

    $$H_{gen} = T(\mathbf{p}) \begin{pmatrix} I - 2\mathbf{n}\mathbf{n}^T & \mathbf{0} \\ \mathbf{0}^T & 1 \end{pmatrix} T(-\mathbf{p})$$

    即先平移使 $\mathbf{p}$ 到原点，做反射，再平移回去。

!!! theorem "定理 67.14 (Cartan-Dieudonne 定理)"
    $\mathbb{R}^n$ 中的每个正交变换都可以表示为**至多 $n$ 次反射的复合**。

    特别地：
    - 二维旋转 = 两次反射的复合。
    - 三维旋转 = 两次反射的复合（绕两个反射平面的交线旋转，旋转角为两平面夹角的两倍）。

    这意味着旋转可以由反射"生成"，这在数值线性代数中有重要应用（Householder QR 分解）。

!!! example "例 67.5"
    **二维反射生成旋转。** 设 $H_1$ 是关于 $x$ 轴的反射（$\mathbf{n}_1 = (0, 1)^T$），$H_2$ 是关于与 $x$ 轴成角 $\theta/2$ 的直线的反射（$\mathbf{n}_2 = (-\sin(\theta/2), \cos(\theta/2))^T$）。则

    $$H_2 H_1 = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} = R(\theta).$$

    复合两次反射（$\det = (-1)^2 = 1$）得到旋转，旋转角等于两反射轴夹角的两倍。

---

## 67.7 法向量变换

<div class="context-flow" markdown>

**核心问题**：当物体被变换时，表面法向量应该如何变换？

</div>

### 法向量不随模型矩阵变换

!!! theorem "定理 67.7 (法向量变换)"
    若物体顶点按矩阵 $M$（$3 \times 3$ 部分）变换，则表面法向量应按 $(M^{-1})^T$ 变换：
    $$
    \mathbf{n}' = (M^{-1})^T \mathbf{n}
    $$

    特别地，若 $M$ 是正交矩阵（$M^{-1} = M^T$），则 $(M^{-1})^T = M$，法向量和顶点用同一矩阵变换。

??? proof "证明"
    设表面上的切向量为 $\mathbf{t}$，法向量 $\mathbf{n}$ 满足 $\mathbf{n}^T \mathbf{t} = 0$。

    变换后切向量为 $\mathbf{t}' = M\mathbf{t}$。我们需要变换后的法向量 $\mathbf{n}'$ 满足 $(\mathbf{n}')^T \mathbf{t}' = 0$。

    设 $\mathbf{n}' = N\mathbf{n}$，则 $(\mathbf{n}')^T \mathbf{t}' = \mathbf{n}^T N^T M \mathbf{t}$。要求这对所有切向量 $\mathbf{t}$ 为零，需要 $N^T M = I$（或更一般地 $N^T M = cI$），即 $N = (M^{-1})^T = (M^T)^{-1}$。

!!! example "例 67.6"
    **非均匀缩放的例子**：设 $M = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$（$x$ 方向拉伸 2 倍，$y$ 不变）。

    考虑直线 $y = x$，切向量 $\mathbf{t} = (1, 1)^T$，法向量 $\mathbf{n} = (1, -1)^T$（$\mathbf{n}^T\mathbf{t} = 0$）。

    若错误地用 $M$ 变换法向量：$M\mathbf{n} = (2, -1)^T$。变换后的切向量 $M\mathbf{t} = (2, 1)^T$。检验：$(2, -1) \cdot (2, 1) = 3 \neq 0$。法向量不再垂直于切向量！

    正确变换：$(M^{-1})^T = \begin{pmatrix} 1/2 & 0 \\ 0 & 1 \end{pmatrix}$，$\mathbf{n}' = (1/2, -1)^T$。检验：$(1/2, -1) \cdot (2, 1) = 0$。正确。

---

## 67.8 颜色空间变换

<div class="context-flow" markdown>

**核心问题**：不同颜色表示之间的转换如何用线性代数来描述？

</div>

### RGB 与 YCbCr

!!! definition "定义 67.12 (RGB 到 YCbCr 变换)"
    **YCbCr** 颜色空间（用于视频压缩）与 **RGB** 之间的转换是线性变换：
    $$
    \begin{pmatrix} Y \\ C_b \\ C_r \end{pmatrix} = \begin{pmatrix} 0.299 & 0.587 & 0.114 \\ -0.169 & -0.331 & 0.500 \\ 0.500 & -0.419 & -0.081 \end{pmatrix} \begin{pmatrix} R \\ G \\ B \end{pmatrix} + \begin{pmatrix} 0 \\ 128 \\ 128 \end{pmatrix}
    $$

    其中 $Y$ 是亮度分量，$C_b, C_r$ 是色度分量。

!!! note "注"
    $Y$ 分量的系数 $(0.299, 0.587, 0.114)$ 反映了人眼对不同颜色的敏感度：绿色最敏感，红色次之，蓝色最不敏感。JPEG 压缩利用这一事实，对色度分量 $C_b, C_r$ 进行更大幅度的下采样。

### 线性 RGB 与 sRGB

!!! definition "定义 67.13 (Gamma 校正)"
    物理上线性的 RGB 值 $C_{linear}$ 与用于显示的 sRGB 值 $C_{sRGB}$ 之间的关系为非线性的 gamma 函数：
    $$
    C_{sRGB} = \begin{cases} 12.92 \, C_{linear} & \text{若 } C_{linear} \leq 0.0031308 \\ 1.055 \, C_{linear}^{1/2.4} - 0.055 & \text{若 } C_{linear} > 0.0031308 \end{cases}
    $$

    这不是线性变换，但在图形管线中必须正确处理——在线性空间中做光照计算，在 sRGB 空间中做显示。

### 色彩空间之间的矩阵变换

!!! theorem "定理 67.8 (XYZ 到 sRGB 矩阵)"
    CIE XYZ 颜色空间到线性 sRGB 的变换为
    $$
    \begin{pmatrix} R \\ G \\ B \end{pmatrix} = \begin{pmatrix} 3.2406 & -1.5372 & -0.4986 \\ -0.9689 & 1.8758 & 0.0415 \\ 0.0557 & -0.2040 & 1.0570 \end{pmatrix} \begin{pmatrix} X \\ Y \\ Z \end{pmatrix}
    $$

    这个 $3 \times 3$ 矩阵完全由 sRGB 的三原色色度坐标和白点（D65）确定。

---

## 67.9 Bezier 曲线与样条

<div class="context-flow" markdown>

**核心问题**：如何用矩阵形式表示和计算参数曲线？

</div>

### Bezier 曲线的矩阵形式

!!! definition "定义 67.14 (三次 Bezier 曲线)"
    给定四个控制点 $\mathbf{P}_0, \mathbf{P}_1, \mathbf{P}_2, \mathbf{P}_3$，三次 Bezier 曲线为
    $$
    \mathbf{C}(t) = (1-t)^3 \mathbf{P}_0 + 3t(1-t)^2 \mathbf{P}_1 + 3t^2(1-t) \mathbf{P}_2 + t^3 \mathbf{P}_3, \quad t \in [0, 1]
    $$

    矩阵形式为
    $$
    \mathbf{C}(t) = \begin{pmatrix} 1 & t & t^2 & t^3 \end{pmatrix} \begin{pmatrix} 1 & 0 & 0 & 0 \\ -3 & 3 & 0 & 0 \\ 3 & -6 & 3 & 0 \\ -1 & 3 & -3 & 1 \end{pmatrix} \begin{pmatrix} \mathbf{P}_0 \\ \mathbf{P}_1 \\ \mathbf{P}_2 \\ \mathbf{P}_3 \end{pmatrix}
    $$

    即 $\mathbf{C}(t) = \mathbf{t}^T M_{Bez} \mathbf{P}$，其中 $\mathbf{t} = (1, t, t^2, t^3)^T$。

!!! theorem "定理 67.9 (Bezier 矩阵的性质)"
    Bezier 基矩阵 $M_{Bez}$ 满足：

    1. $M_{Bez}$ 的每行之和分别为 $1, 0, 0, 0$（保证 $\mathbf{C}(0) = \mathbf{P}_0$）；
    2. $M_{Bez}$ 是下三角矩阵，行列式 $\det(M_{Bez}) = 1 \cdot 3 \cdot 3 \cdot 1 = 9$（非奇异）。

### B 样条基矩阵

!!! definition "定义 67.15 (均匀三次 B 样条)"
    均匀三次 B 样条曲线段的矩阵形式为
    $$
    \mathbf{C}(t) = \mathbf{t}^T M_{B} \mathbf{P}
    $$
    其中 B 样条基矩阵为
    $$
    M_B = \frac{1}{6}\begin{pmatrix} 1 & 4 & 1 & 0 \\ -3 & 0 & 3 & 0 \\ 3 & -6 & 3 & 0 \\ -1 & 3 & -3 & 1 \end{pmatrix}
    $$

!!! note "注"
    B 样条与 Bezier 曲线的关键区别在于 B 样条具有**局部支撑性**：移动一个控制点只影响附近的曲线段。两者通过基矩阵的变换相互关联：$M_B = M_{Bez} \cdot T_{B \to Bez}$。

### de Casteljau 算法

!!! theorem "定理 67.10 (de Casteljau 算法)"
    Bezier 曲线在参数 $t$ 处的值可以通过递归线性插值计算：
    $$
    \mathbf{P}_i^{(r)} = (1-t)\mathbf{P}_i^{(r-1)} + t \, \mathbf{P}_{i+1}^{(r-1)}, \quad r = 1, \ldots, n
    $$
    从 $\mathbf{P}_i^{(0)} = \mathbf{P}_i$ 开始，最终 $\mathbf{P}_0^{(n)} = \mathbf{C}(t)$。

    这个过程在数值上比直接计算更稳定，且中间结果提供了曲线的细分。

### 曲面片

!!! definition "定义 67.16 (双三次 Bezier 曲面片)"
    双三次 Bezier 曲面片由 $4 \times 4$ 的控制点网格 $\mathbf{P}_{ij}$ 定义：
    $$
    \mathbf{S}(u, v) = \mathbf{u}^T M_{Bez} \, \mathbf{G} \, M_{Bez}^T \mathbf{v}
    $$
    其中 $\mathbf{u} = (1, u, u^2, u^3)^T$，$\mathbf{v} = (1, v, v^2, v^3)^T$，$\mathbf{G}$ 是 $4 \times 4$ 控制点矩阵。

!!! example "例 67.7"
    **茶壶模型**（Utah teapot）——计算机图形学的标志性模型——由 32 个双三次 Bezier 曲面片组成，每个片由 16 个控制点定义。整个茶壶共有约 306 个控制点（因为相邻曲面片共享边界控制点）。

---

## 本章小结

本章展示了线性代数在计算机图形学中无处不在的应用。主要内容包括：

1. **齐次坐标**通过引入第四个分量 $w$，将所有仿射变换（包括平移）统一为 $4 \times 4$ 矩阵乘法。

2. **基本几何变换**（平移、旋转、缩放、剪切）都有简洁的矩阵表示，变换组合对应矩阵乘法。顺序至关重要，因为矩阵乘法不满足交换律。

3. **旋转表示**有三种主要方式：旋转矩阵（Rodrigues 公式）、Euler 角（有万向锁问题）和四元数（无奇异性，适合插值）。**SLERP** 球面线性插值提供等角速度的最短路径旋转插值。

4. **MVP 管线**将顶点从局部空间经由模型、视图、投影三个矩阵变换到裁剪空间，经齐次除法到 NDC，最后通过**视口变换**映射到屏幕像素坐标。

5. **透视投影**通过齐次除法自动实现近大远小的效果，深度值的非线性分布导致远处精度不足。**透视校正插值**通过 $1/w$ 加权修正屏幕空间重心坐标的非线性失真。

6. **反射矩阵** $H = I - 2\mathbf{n}\mathbf{n}^T$ 是正交但行列式为 $-1$ 的变换。Cartan-Dieudonne 定理表明每个正交变换都可以分解为至多 $n$ 次反射。

7. **法向量变换**使用 $(M^{-1})^T$ 而非 $M$，这是许多图形学初学者容易忽略的重要细节。

8. **颜色空间变换**（RGB-YCbCr、XYZ-sRGB）是矩阵乘法（加上 gamma 校正这一非线性步骤）。

9. **Bezier 曲线和 B 样条**都有优雅的矩阵表示，使得曲线计算和连续性分析都可以用矩阵语言进行。

---

## 习题

!!! question "习题 67.1"
    写出将物体绕点 $(1, 2, 3)$ 旋转 $90°$（绕 $z$ 轴）的 $4 \times 4$ 变换矩阵。

!!! question "习题 67.2"
    证明两个平移矩阵 $T(\mathbf{a})$ 和 $T(\mathbf{b})$ 的乘积是 $T(\mathbf{a} + \mathbf{b})$。

!!! question "习题 67.3"
    将旋转轴 $\hat{\mathbf{k}} = (1, 1, 1)/\sqrt{3}$、旋转角 $\theta = 120°$ 代入 Rodrigues 公式，求旋转矩阵。

!!! question "习题 67.4"
    验证四元数 $\mathbf{q} = (\cos(\pi/4), 0, 0, \sin(\pi/4))$ 对应绕 $z$ 轴旋转 $90°$。将其转换为旋转矩阵。

!!! question "习题 67.5"
    推导正交投影矩阵 $P_{ortho}$（将视景体 $[l,r] \times [b,t] \times [-n,-f]$ 映射到 $[-1,1]^3$）。

!!! question "习题 67.6"
    设 $M = \begin{pmatrix} 1 & 1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$（$xy$ 剪切）。计算法向量变换矩阵 $(M^{-1})^T$。

!!! question "习题 67.7"
    证明均匀缩放矩阵 $S(s, s, s)$ 的法向量变换矩阵仍然是 $S(s, s, s)$（方向不变，只是缩放法向量长度）。

!!! question "习题 67.8"
    计算 Bezier 基矩阵 $M_{Bez}$ 的逆，并解释其几何含义。

!!! question "习题 67.9"
    证明三次 Bezier 曲线 $\mathbf{C}(t)$ 满足 $\mathbf{C}(0) = \mathbf{P}_0$ 且 $\mathbf{C}(1) = \mathbf{P}_3$。

!!! question "习题 67.10"
    给定 $2 \times 2$ 的二次 Bezier 曲面片（$3 \times 3$ 控制点网格），写出矩阵形式。
