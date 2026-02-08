# 第 67 章 线性代数在计算机图形学中的应用

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 线性变换(Ch5) · 正交矩阵(Ch7) · 特征值(Ch6)

**本章脉络**：齐次坐标 → 仿射变换 → 旋转表示(矩阵/四元数/Euler角) → 模型-视图-投影矩阵 → 透视与正交投影 → 法向量变换 → 颜色空间变换 → 样条曲面

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

---

## 67.6 法向量变换

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

!!! example "例 67.5"
    **非均匀缩放的例子**：设 $M = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$（$x$ 方向拉伸 2 倍，$y$ 不变）。

    考虑直线 $y = x$，切向量 $\mathbf{t} = (1, 1)^T$，法向量 $\mathbf{n} = (1, -1)^T$（$\mathbf{n}^T\mathbf{t} = 0$）。

    若错误地用 $M$ 变换法向量：$M\mathbf{n} = (2, -1)^T$。变换后的切向量 $M\mathbf{t} = (2, 1)^T$。检验：$(2, -1) \cdot (2, 1) = 3 \neq 0$。法向量不再垂直于切向量！

    正确变换：$(M^{-1})^T = \begin{pmatrix} 1/2 & 0 \\ 0 & 1 \end{pmatrix}$，$\mathbf{n}' = (1/2, -1)^T$。检验：$(1/2, -1) \cdot (2, 1) = 0$。正确。

---

## 67.7 颜色空间变换

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

## 67.8 Bezier 曲线与样条

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

!!! example "例 67.6"
    **茶壶模型**（Utah teapot）——计算机图形学的标志性模型——由 32 个双三次 Bezier 曲面片组成，每个片由 16 个控制点定义。整个茶壶共有约 306 个控制点（因为相邻曲面片共享边界控制点）。

---

## 本章小结

本章展示了线性代数在计算机图形学中无处不在的应用。主要内容包括：

1. **齐次坐标**通过引入第四个分量 $w$，将所有仿射变换（包括平移）统一为 $4 \times 4$ 矩阵乘法。

2. **基本几何变换**（平移、旋转、缩放、剪切）都有简洁的矩阵表示，变换组合对应矩阵乘法。顺序至关重要，因为矩阵乘法不满足交换律。

3. **旋转表示**有三种主要方式：旋转矩阵（Rodrigues 公式）、Euler 角（有万向锁问题）和四元数（无奇异性，适合插值）。

4. **MVP 管线**将顶点从局部空间经由模型、视图、投影三个矩阵变换到裁剪空间，是实时渲染的核心。

5. **透视投影**通过齐次除法自动实现近大远小的效果，深度值的非线性分布导致远处精度不足。

6. **法向量变换**使用 $(M^{-1})^T$ 而非 $M$，这是许多图形学初学者容易忽略的重要细节。

7. **颜色空间变换**（RGB-YCbCr、XYZ-sRGB）是矩阵乘法（加上 gamma 校正这一非线性步骤）。

8. **Bezier 曲线和 B 样条**都有优雅的矩阵表示，使得曲线计算和连续性分析都可以用矩阵语言进行。

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
