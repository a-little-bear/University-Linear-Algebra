# 第 26 章 线性代数在微分方程中的应用

<div class="context-flow" markdown>

**前置**：特征值与对角化(Ch6) · Jordan 标准形(Ch12) · 矩阵函数(Ch13)

**脉络**：$\mathbf{x}' = A\mathbf{x}$(齐次系统) → $e^{At}$(矩阵指数) → 特征值实部(稳定性) → 相平面(几何分类) → 友矩阵(高阶化为系统) → 变参法(非齐次) → Sturm-Liouville(PDE) → Floquet(周期系统)

</div>

线性微分方程是线性代数最经典也最重要的应用领域之一。线性常微分方程组的解可以用矩阵指数 $e^{At}$ 完整描述，而系统的长期行为（稳定性、振荡、衰减）完全由系数矩阵 $A$ 的特征值决定。本章从一阶常系数线性系统出发，依次展开矩阵指数、稳定性分析、相平面分类、高阶方程、非齐次方程、偏微分方程中的特征值问题，以及周期系数系统的 Floquet 理论。

---

## 26.1 线性常微分方程组

<div class="context-flow" markdown>

**核心结构**：$\mathbf{x}' = A\mathbf{x}$ 的解空间是 $n$ 维线性空间 → 基本解矩阵 $\Phi(t)$ 的每一列是一个解 → $\Phi(t) = e^{At}$ 当 $\Phi(0) = I$

**链接**：Ch4 线性空间的维数定理在此直接体现

</div>

线性常微分方程组是动力系统、控制理论、物理学的基本模型。

!!! definition "定义 26.1 (线性常系数 ODE 系统)"
    设 $A \in \mathbb{R}^{n \times n}$ 为常数矩阵，**线性常系数常微分方程组**（linear constant-coefficient ODE system）为

    $$
    \mathbf{x}'(t) = A\mathbf{x}(t), \quad \mathbf{x}(0) = \mathbf{x}_0,
    $$

    其中 $\mathbf{x}(t) \in \mathbb{R}^n$ 为状态向量，$t \ge 0$。

!!! definition "定义 26.2 (基本解矩阵)"
    $n \times n$ 矩阵值函数 $\Phi(t)$ 称为方程 $\mathbf{x}' = A\mathbf{x}$ 的**基本解矩阵**（fundamental matrix），若

    $$
    \Phi'(t) = A\Phi(t), \quad \det \Phi(t) \neq 0, \quad \forall\, t.
    $$

    若进一步 $\Phi(0) = I$，则称 $\Phi(t)$ 为**主基本解矩阵**（principal fundamental matrix）或**状态转移矩阵**。

!!! theorem "定理 26.1 (解空间结构)"
    齐次方程 $\mathbf{x}' = A\mathbf{x}$ 的解集构成 $\mathbb{R}^n$ 上的 $n$ 维线性空间。若 $\Phi(t)$ 为任一基本解矩阵，则一般解为

    $$
    \mathbf{x}(t) = \Phi(t)\mathbf{c}, \quad \mathbf{c} \in \mathbb{R}^n.
    $$

    满足初始条件 $\mathbf{x}(0) = \mathbf{x}_0$ 的唯一解为 $\mathbf{x}(t) = \Phi(t)\Phi(0)^{-1}\mathbf{x}_0$。

??? proof "证明"
    **线性性**：若 $\mathbf{x}_1(t)$ 和 $\mathbf{x}_2(t)$ 是解，则 $\alpha \mathbf{x}_1(t) + \beta \mathbf{x}_2(t)$ 也是解，因为

    $$
    (\alpha \mathbf{x}_1 + \beta \mathbf{x}_2)' = \alpha A\mathbf{x}_1 + \beta A\mathbf{x}_2 = A(\alpha \mathbf{x}_1 + \beta \mathbf{x}_2).
    $$

    **维数为 $n$**：ODE 存在唯一性定理保证对每个初始值 $\mathbf{x}_0 \in \mathbb{R}^n$ 存在唯一解。初始值到解的映射 $\mathbf{x}_0 \mapsto \mathbf{x}(t)$ 是线性双射，因此解空间与 $\mathbb{R}^n$ 同构，维数为 $n$。

    基本解矩阵 $\Phi(t)$ 的 $n$ 列线性无关（因 $\det \Phi(t) \neq 0$），构成解空间的一组基。$\blacksquare$

!!! example "例 26.1"
    **二维线性系统的解。** 考虑

    $$
    \mathbf{x}' = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix} \mathbf{x}.
    $$

    特征方程 $\lambda^2 + 3\lambda + 2 = 0$ 给出 $\lambda_1 = -1$，$\lambda_2 = -2$。对应特征向量为 $\mathbf{v}_1 = (1, -1)^T$，$\mathbf{v}_2 = (1, -2)^T$。一般解为

    $$
    \mathbf{x}(t) = c_1 e^{-t} \begin{pmatrix} 1 \\ -1 \end{pmatrix} + c_2 e^{-2t} \begin{pmatrix} 1 \\ -2 \end{pmatrix}.
    $$

    基本解矩阵为 $\Phi(t) = \begin{pmatrix} e^{-t} & e^{-2t} \\ -e^{-t} & -2e^{-2t} \end{pmatrix}$，满足 $\Phi'(t) = A\Phi(t)$。

!!! theorem "定理 26.2 (Liouville 公式)"
    设 $\Phi(t)$ 为 $\mathbf{x}' = A(t)\mathbf{x}$ 的基本解矩阵（$A(t)$ 可为时变），则

    $$
    \det \Phi(t) = \det \Phi(t_0) \cdot \exp\!\left(\int_{t_0}^{t} \operatorname{tr} A(s)\, ds\right).
    $$

    特别地，对常系数系统 $A(t) \equiv A$，有 $\det \Phi(t) = \det \Phi(0) \cdot e^{t\, \operatorname{tr}(A)}$。

??? proof "证明"
    记 $W(t) = \det \Phi(t)$（Wronskian）。利用行列式对矩阵求导的公式，

    $$
    W'(t) = \sum_{i=1}^{n} \det \Phi_i(t),
    $$

    其中 $\Phi_i(t)$ 是将 $\Phi(t)$ 的第 $i$ 行替换为其导数行。由 $\Phi' = A\Phi$，第 $i$ 行的导数等于 $A$ 的第 $i$ 行与 $\Phi$ 的乘积。展开后只有对角项 $a_{ii}$ 对行列式有贡献，因此 $W'(t) = \operatorname{tr}(A(t)) W(t)$。解此标量 ODE 即得结论。$\blacksquare$

---

## 26.2 矩阵指数与解的表示

<div class="context-flow" markdown>

**核心公式**：$e^{At} = \sum_{k=0}^{\infty} \frac{(At)^k}{k!}$ → 可对角化时 $e^{At} = P\, \text{diag}(e^{\lambda_i t})\, P^{-1}$ → Jordan 块时含 $t^k e^{\lambda t}$ 项

**链接**：Ch13 矩阵函数理论的最重要特例

</div>

矩阵指数是线性 ODE 理论的核心工具，将矩阵代数与微分方程完美统一。

!!! definition "定义 26.3 (矩阵指数)"
    对 $A \in \mathbb{R}^{n \times n}$，**矩阵指数**（matrix exponential）定义为

    $$
    e^{A} = \sum_{k=0}^{\infty} \frac{A^k}{k!} = I + A + \frac{A^2}{2!} + \frac{A^3}{3!} + \cdots
    $$

    此级数对所有 $A$ 绝对收敛。

!!! theorem "定理 26.3 (矩阵指数与 ODE 的关系)"
    $\mathbf{x}' = A\mathbf{x}$，$\mathbf{x}(0) = \mathbf{x}_0$ 的唯一解为

    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}_0.
    $$

    $e^{At}$ 是主基本解矩阵：$\frac{d}{dt}e^{At} = Ae^{At}$，$e^{A \cdot 0} = I$。

??? proof "证明"
    逐项求导：

    $$
    \frac{d}{dt}e^{At} = \frac{d}{dt}\sum_{k=0}^{\infty} \frac{(At)^k}{k!} = \sum_{k=1}^{\infty} \frac{A^k t^{k-1}}{(k-1)!} = A\sum_{k=1}^{\infty} \frac{(At)^{k-1}}{(k-1)!} = Ae^{At}.
    $$

    因此 $\mathbf{x}(t) = e^{At}\mathbf{x}_0$ 满足 $\mathbf{x}' = Ae^{At}\mathbf{x}_0 = A\mathbf{x}(t)$ 以及 $\mathbf{x}(0) = I\mathbf{x}_0 = \mathbf{x}_0$。由唯一性定理，这是唯一解。$\blacksquare$

!!! theorem "定理 26.4 (矩阵指数的计算)"
    设 $A = PJP^{-1}$ 为 Jordan 分解，$J = \operatorname{diag}(J_1, \ldots, J_s)$，其中 $J_i$ 为 Jordan 块。则

    $$
    e^{At} = P\, \operatorname{diag}(e^{J_1 t}, \ldots, e^{J_s t})\, P^{-1}.
    $$

    对 $k \times k$ Jordan 块 $J_i = \lambda_i I + N$（$N$ 为幂零），有

    $$
    e^{J_i t} = e^{\lambda_i t} \begin{pmatrix} 1 & t & \frac{t^2}{2!} & \cdots & \frac{t^{k-1}}{(k-1)!} \\ 0 & 1 & t & \cdots & \frac{t^{k-2}}{(k-2)!} \\ \vdots & & \ddots & & \vdots \\ 0 & 0 & \cdots & 1 & t \\ 0 & 0 & \cdots & 0 & 1 \end{pmatrix}.
    $$

??? proof "证明"
    由 $A = PJP^{-1}$ 得 $A^k = PJ^k P^{-1}$，故

    $$
    e^{At} = \sum_{k=0}^{\infty} \frac{(PJP^{-1})^k t^k}{k!} = P\left(\sum_{k=0}^{\infty} \frac{J^k t^k}{k!}\right)P^{-1} = Pe^{Jt}P^{-1}.
    $$

    由于 $J$ 为分块对角，$e^{Jt} = \operatorname{diag}(e^{J_1 t}, \ldots, e^{J_s t})$。对 Jordan 块 $J_i = \lambda_i I + N$，由于 $\lambda_i I$ 和 $N$ 可交换，$e^{J_i t} = e^{\lambda_i t} e^{Nt}$，其中 $N^k = 0$（$k \ge $ 块大小），故 $e^{Nt}$ 为有限和，给出上三角矩阵形式。$\blacksquare$

!!! example "例 26.2"
    **含重特征值系统的矩阵指数。** 设

    $$
    A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}.
    $$

    $A$ 有唯一特征值 $\lambda = 2$，Jordan 块为 $J = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$（$P = I$）。

    $$
    e^{At} = e^{2t}\begin{pmatrix} 1 & t \\ 0 & 1 \end{pmatrix}.
    $$

    初值 $\mathbf{x}_0 = (1, 3)^T$ 的解为

    $$
    \mathbf{x}(t) = e^{2t}\begin{pmatrix} 1 + 3t \\ 3 \end{pmatrix}.
    $$

    解中出现 $te^{2t}$ 项，这是 Jordan 块（非对角化）的标志性行为。

!!! definition "定义 26.4 (矩阵指数的基本性质)"
    矩阵指数满足以下性质：

    1. $e^{O} = I$（$O$ 为零矩阵）。
    2. $(e^{A})^{-1} = e^{-A}$。
    3. $e^{(s+t)A} = e^{sA}e^{tA}$（半群性质）。
    4. 若 $AB = BA$，则 $e^{A+B} = e^{A}e^{B}$。
    5. $\det(e^{A}) = e^{\operatorname{tr}(A)}$。
    6. 一般地，$e^{A+B} \neq e^{A}e^{B}$。

!!! example "例 26.3"
    **矩阵指数的非交换性。** 令 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}$。

    $e^A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$，$e^B = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}$，$e^A e^B = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$。

    而 $A + B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ 的特征值为 $\pm 1$，

    $$
    e^{A+B} = \begin{pmatrix} \cosh 1 & \sinh 1 \\ \sinh 1 & \cosh 1 \end{pmatrix} \neq e^A e^B.
    $$

    $AB \neq BA$，故 Baker-Campbell-Hausdorff 公式给出非平凡修正项。

---

## 26.3 稳定性分析

<div class="context-flow" markdown>

**判据链**：所有 $\operatorname{Re}(\lambda_i) < 0$ ⇔ 渐近稳定 → Lyapunov 方程 $A^T P + PA = -Q$ 有正定解 ⇔ Hurwitz 矩阵

**应用**：控制理论中系统是否收敛、信号处理中滤波器是否稳定

</div>

线性系统的稳定性完全由系数矩阵的特征值决定，这是线性代数在动力系统中最深刻的应用之一。

!!! definition "定义 26.5 (稳定性分类)"
    对系统 $\mathbf{x}' = A\mathbf{x}$，平衡点 $\mathbf{x} = \mathbf{0}$ 称为：

    - **渐近稳定**（asymptotically stable）：若对所有初值 $\mathbf{x}_0$，$\lim_{t \to \infty} \mathbf{x}(t) = \mathbf{0}$；
    - **稳定**（stable / Lyapunov stable）：若对任意 $\varepsilon > 0$，存在 $\delta > 0$ 使得 $\|\mathbf{x}_0\| < \delta$ 蕴含 $\|\mathbf{x}(t)\| < \varepsilon$，$\forall\, t \ge 0$；
    - **不稳定**（unstable）：若不稳定。

!!! theorem "定理 26.5 (特征值稳定性准则)"
    对常系数系统 $\mathbf{x}' = A\mathbf{x}$：

    1. **渐近稳定** $\iff$ $A$ 的所有特征值满足 $\operatorname{Re}(\lambda_i) < 0$（$A$ 为 **Hurwitz 矩阵**）。
    2. **稳定** $\iff$ 所有特征值满足 $\operatorname{Re}(\lambda_i) \le 0$，且 $\operatorname{Re}(\lambda_i) = 0$ 的特征值的 Jordan 块大小为 $1$。
    3. **不稳定** $\iff$ 存在 $\operatorname{Re}(\lambda_i) > 0$，或存在 $\operatorname{Re}(\lambda_i) = 0$ 且 Jordan 块大小 $> 1$。

??? proof "证明"
    由 $\mathbf{x}(t) = e^{At}\mathbf{x}_0$ 和 Jordan 分解 $e^{At} = Pe^{Jt}P^{-1}$，每个 Jordan 块的贡献为 $t^j e^{\lambda_i t}$（$0 \le j < $ 块大小）。

    - 若 $\operatorname{Re}(\lambda_i) < 0$，则 $|t^j e^{\lambda_i t}| = t^j e^{\operatorname{Re}(\lambda_i) t} \to 0$（$t \to \infty$），无论 $j$ 多大。
    - 若 $\operatorname{Re}(\lambda_i) = 0$ 且 $j = 0$，则 $|e^{\lambda_i t}| = 1$，有界但不趋零。
    - 若 $\operatorname{Re}(\lambda_i) = 0$ 且 $j \ge 1$，则 $|t^j e^{\lambda_i t}| = t^j \to \infty$。
    - 若 $\operatorname{Re}(\lambda_i) > 0$，则 $|e^{\lambda_i t}| \to \infty$。$\blacksquare$

!!! theorem "定理 26.6 (Lyapunov 稳定性定理)"
    $A$ 为 Hurwitz 矩阵当且仅当对任意正定矩阵 $Q \succ 0$，**Lyapunov 方程**

    $$
    A^T P + PA = -Q
    $$

    有唯一正定解 $P \succ 0$。此时 $V(\mathbf{x}) = \mathbf{x}^T P \mathbf{x}$ 是系统的 Lyapunov 函数。

??? proof "证明"
    **充分性**：若 $P \succ 0$ 满足 $A^T P + PA = -Q \prec 0$，定义 $V(\mathbf{x}) = \mathbf{x}^T P \mathbf{x} > 0$。沿轨迹求导：

    $$
    \dot{V} = \mathbf{x}'^T P \mathbf{x} + \mathbf{x}^T P \mathbf{x}' = \mathbf{x}^T(A^T P + PA)\mathbf{x} = -\mathbf{x}^T Q \mathbf{x} < 0.
    $$

    **必要性**：若 $A$ 为 Hurwitz 矩阵，定义 $P = \int_0^{\infty} e^{A^T t} Q e^{At}\, dt$（收敛因 $\operatorname{Re}(\lambda_i) < 0$）。直接验证 $A^T P + PA = -Q$ 且 $P \succ 0$。唯一性由 Lyapunov 方程的线性性和 Hurwitz 条件保证。$\blacksquare$

!!! example "例 26.4"
    **电路系统的稳定性分析。** 一个 RLC 电路的状态方程为

    $$
    A = \begin{pmatrix} 0 & 1 \\ -\frac{1}{LC} & -\frac{R}{L} \end{pmatrix}.
    $$

    特征方程为 $\lambda^2 + \frac{R}{L}\lambda + \frac{1}{LC} = 0$，特征值为

    $$
    \lambda_{1,2} = -\frac{R}{2L} \pm \sqrt{\frac{R^2}{4L^2} - \frac{1}{LC}}.
    $$

    由于 $R, L, C > 0$，有 $\operatorname{Re}(\lambda_{1,2}) < 0$，系统渐近稳定。物理意义：电阻耗散能量，电流和电压最终衰减至零。当 $R^2 > 4L/C$ 时为过阻尼（两个实负特征值），$R^2 < 4L/C$ 时为欠阻尼（共轭复特征值，衰减振荡），$R^2 = 4L/C$ 时为临界阻尼。

---

## 26.4 相平面分析

<div class="context-flow" markdown>

**分类**：$2 \times 2$ 系统的完整分类依赖于 $\operatorname{tr}(A)$ 和 $\det(A)$ → 特征值符号和虚实性决定结点/鞍点/焦点/中心

**几何**：轨迹的定性形状完全由特征值和特征向量决定

</div>

二维线性系统 $\mathbf{x}' = A\mathbf{x}$ 的相平面（phase plane）分析是理解高维系统的起点。

!!! definition "定义 26.6 (临界点分类)"
    对二维系统 $\mathbf{x}' = A\mathbf{x}$（$A \in \mathbb{R}^{2 \times 2}$，$\det A \neq 0$），设特征值为 $\lambda_1, \lambda_2$，原点的临界点类型为：

    | 特征值条件 | 类型 |
    |:---:|:---:|
    | $\lambda_1, \lambda_2 \in \mathbb{R}$，同号 | **结点**（node） |
    | $\lambda_1, \lambda_2 \in \mathbb{R}$，异号 | **鞍点**（saddle） |
    | $\lambda_{1,2} = \alpha \pm \beta i$，$\alpha \neq 0$ | **焦点/螺旋点**（spiral/focus） |
    | $\lambda_{1,2} = \pm \beta i$，$\alpha = 0$ | **中心**（center） |
    | $\lambda_1 = \lambda_2 \in \mathbb{R}$，非对角化 | **退化结点/星形结点** |

!!! theorem "定理 26.7 (迹-行列式分类)"
    设 $\tau = \operatorname{tr}(A)$，$\delta = \det(A)$，$\Delta = \tau^2 - 4\delta$。则：

    - $\delta < 0$：**鞍点**（不稳定）。
    - $\delta > 0$，$\Delta > 0$，$\tau < 0$：**稳定结点**。
    - $\delta > 0$，$\Delta > 0$，$\tau > 0$：**不稳定结点**。
    - $\delta > 0$，$\Delta < 0$，$\tau < 0$：**稳定焦点**。
    - $\delta > 0$，$\Delta < 0$，$\tau > 0$：**不稳定焦点**。
    - $\delta > 0$，$\tau = 0$：**中心**。

??? proof "证明"
    特征值为 $\lambda_{1,2} = \frac{\tau \pm \sqrt{\Delta}}{2}$，其中 $\lambda_1 \lambda_2 = \delta$，$\lambda_1 + \lambda_2 = \tau$。

    - $\delta < 0$ 意味着 $\lambda_1, \lambda_2$ 异号。
    - $\delta > 0$，$\Delta > 0$ 时为两个同号实根，$\tau$ 的正负决定稳定性。
    - $\delta > 0$，$\Delta < 0$ 时为共轭复根 $\alpha \pm \beta i$（$\alpha = \tau/2$），$\alpha$ 的正负决定稳定性。
    - $\tau = 0$ 时 $\alpha = 0$，为纯虚特征值。$\blacksquare$

!!! example "例 26.5"
    **弹簧-阻尼系统的相平面。** 弹簧质量系统 $mx'' + cx' + kx = 0$ 化为

    $$
    \begin{pmatrix} x' \\ v' \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ -k/m & -c/m \end{pmatrix} \begin{pmatrix} x \\ v \end{pmatrix}.
    $$

    此处 $\tau = -c/m < 0$，$\delta = k/m > 0$，$\Delta = c^2/m^2 - 4k/m$。

    - **欠阻尼**（$c^2 < 4mk$）：$\Delta < 0$，稳定焦点。轨迹以螺旋线趋向原点。
    - **过阻尼**（$c^2 > 4mk$）：$\Delta > 0$，稳定结点。轨迹沿特征方向单调趋向原点。
    - **无阻尼**（$c = 0$）：$\tau = 0$，中心。轨迹为以原点为中心的椭圆，对应能量守恒。

---

## 26.5 高阶线性 ODE 与友矩阵

<div class="context-flow" markdown>

**化归**：$n$ 阶线性 ODE → 一阶 $n \times n$ 系统 → 友矩阵 $C$ 的特征多项式 = 原方程的特征多项式

**链接**：Ch6 特征多项式理论在此直接应用

</div>

高阶线性 ODE 可以通过引入状态变量转化为一阶系统，对应的系数矩阵具有特殊的友矩阵结构。

!!! definition "定义 26.7 (高阶线性 ODE)"
    $n$ 阶线性常系数 ODE 为

    $$
    y^{(n)} + a_{n-1}y^{(n-1)} + \cdots + a_1 y' + a_0 y = 0.
    $$

!!! definition "定义 26.8 (友矩阵)"
    上述方程对应的**友矩阵**（companion matrix）为

    $$
    C = \begin{pmatrix} 0 & 1 & 0 & \cdots & 0 \\ 0 & 0 & 1 & \cdots & 0 \\ \vdots & & & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & 1 \\ -a_0 & -a_1 & -a_2 & \cdots & -a_{n-1} \end{pmatrix} \in \mathbb{R}^{n \times n}.
    $$

    $C$ 的特征多项式恰为 $\lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_1\lambda + a_0$。

!!! theorem "定理 26.8 (高阶 ODE 与一阶系统的等价性)"
    令 $x_1 = y$，$x_2 = y'$，$\ldots$，$x_n = y^{(n-1)}$，则高阶 ODE 等价于一阶系统

    $$
    \mathbf{x}' = C\mathbf{x},
    $$

    其中 $\mathbf{x} = (x_1, \ldots, x_n)^T$，$C$ 为友矩阵。$C$ 的特征值即为原方程的特征根，$C$ 的 Jordan 结构决定解的形式。

??? proof "证明"
    由定义，$x_i' = x_{i+1}$（$1 \le i \le n-1$），且

    $$
    x_n' = y^{(n)} = -a_{n-1}y^{(n-1)} - \cdots - a_0 y = -a_{n-1}x_n - \cdots - a_0 x_1.
    $$

    这恰好是 $\mathbf{x}' = C\mathbf{x}$ 的形式。$C$ 的特征多项式：沿最后一行展开行列式 $\det(\lambda I - C)$，利用上方的 $1$ 结构逐步化简，得 $\lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_0$。$\blacksquare$

!!! example "例 26.6"
    **四阶方程化为系统。** 方程 $y^{(4)} + 2y''' + 3y'' + 2y' + y = 0$ 对应友矩阵

    $$
    C = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ -1 & -2 & -3 & -2 \end{pmatrix}.
    $$

    特征多项式 $\lambda^4 + 2\lambda^3 + 3\lambda^2 + 2\lambda + 1 = (\lambda^2 + \lambda + 1)^2$。特征值为 $\lambda = \frac{-1 \pm i\sqrt{3}}{2}$（二重）。解包含 $e^{-t/2}\cos(\frac{\sqrt{3}}{2}t)$，$e^{-t/2}\sin(\frac{\sqrt{3}}{2}t)$ 及其乘以 $t$ 的项。

---

## 26.6 非齐次方程与变参法

<div class="context-flow" markdown>

**非齐次通解** = 齐次通解 + 特解 → 特解由 $e^{At}$ 的卷积积分给出（变参法）→ 实质是 Green 函数/脉冲响应

**链接**：Ch4 商空间和仿射子空间的直接体现

</div>

非齐次线性系统 $\mathbf{x}' = A\mathbf{x} + \mathbf{f}(t)$ 的求解依赖矩阵指数和变参法。

!!! definition "定义 26.9 (非齐次线性 ODE 系统)"
    **非齐次线性 ODE 系统**为

    $$
    \mathbf{x}'(t) = A\mathbf{x}(t) + \mathbf{f}(t), \quad \mathbf{x}(0) = \mathbf{x}_0,
    $$

    其中 $\mathbf{f}(t) \in \mathbb{R}^n$ 为给定的连续驱动项（forcing term）。

!!! theorem "定理 26.9 (变参法/常数变易公式)"
    非齐次系统 $\mathbf{x}' = A\mathbf{x} + \mathbf{f}(t)$，$\mathbf{x}(0) = \mathbf{x}_0$ 的解为

    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)}\mathbf{f}(s)\, ds.
    $$

    第一项为齐次解（自由响应），第二项为特解（受迫响应），即矩阵指数与驱动项的卷积。

??? proof "证明"
    设 $\mathbf{x}(t) = e^{At}\mathbf{c}(t)$（变参法的核心 Ansatz）。代入方程：

    $$
    Ae^{At}\mathbf{c}(t) + e^{At}\mathbf{c}'(t) = Ae^{At}\mathbf{c}(t) + \mathbf{f}(t).
    $$

    消去 $Ae^{At}\mathbf{c}(t)$ 得 $e^{At}\mathbf{c}'(t) = \mathbf{f}(t)$，即 $\mathbf{c}'(t) = e^{-At}\mathbf{f}(t)$。积分：

    $$
    \mathbf{c}(t) = \mathbf{x}_0 + \int_0^t e^{-As}\mathbf{f}(s)\, ds.
    $$

    因此 $\mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)}\mathbf{f}(s)\, ds$。$\blacksquare$

!!! example "例 26.7"
    **受迫振动系统。** 考虑系统

    $$
    \mathbf{x}' = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}\mathbf{x} + \begin{pmatrix} 0 \\ \cos t \end{pmatrix}, \quad \mathbf{x}(0) = \mathbf{0}.
    $$

    齐次部分有特征值 $\pm i$（中心），$e^{At} = \begin{pmatrix} \cos t & \sin t \\ -\sin t & \cos t \end{pmatrix}$。

    受迫响应为

    $$
    \int_0^t \begin{pmatrix} \cos(t-s) & \sin(t-s) \\ -\sin(t-s) & \cos(t-s) \end{pmatrix} \begin{pmatrix} 0 \\ \cos s \end{pmatrix} ds.
    $$

    计算第一分量：$\int_0^t \sin(t-s)\cos s\, ds = \frac{t}{2}\sin t$（共振）。由于驱动频率与固有频率相同，振幅线性增长——这是共振现象的线性代数解释。

---

## 26.7 线性 PDE 中的特征值问题

<div class="context-flow" markdown>

**分离变量** → 空间部分为特征值问题 $\mathcal{L}u = \lambda u$ → Sturm-Liouville 理论保证实特征值、正交特征函数 → 无穷维内积空间上的谱定理

**链接**：Ch8 内积空间谱定理的无穷维推广

</div>

偏微分方程中的分离变量法将 PDE 化为 ODE 特征值问题，Sturm-Liouville 理论是有限维谱定理的无穷维推广。

!!! definition "定义 26.10 (Sturm-Liouville 问题)"
    **Sturm-Liouville 问题**为求 $\lambda$ 和非零函数 $y(x)$（$x \in [a, b]$）使得

    $$
    -(p(x)y')' + q(x)y = \lambda w(x) y,
    $$

    加上边界条件（如 $y(a) = y(b) = 0$），其中 $p(x) > 0$，$w(x) > 0$ 为权函数。

!!! theorem "定理 26.10 (Sturm-Liouville 谱定理)"
    在适当正则性条件下，Sturm-Liouville 问题具有如下性质：

    1. 存在可数无穷个实特征值 $\lambda_1 < \lambda_2 < \cdots$，$\lambda_n \to \infty$。
    2. 对应的特征函数 $\{y_n(x)\}$ 关于权函数 $w(x)$ 正交：

        $$
        \int_a^b y_m(x) y_n(x) w(x)\, dx = 0, \quad m \neq n.
        $$

    3. $\{y_n\}$ 构成 $L^2_w([a, b])$ 的完备正交基（类比有限维对称矩阵的正交特征基）。

??? proof "证明"
    **自伴性**：定义算子 $\mathcal{L}y = \frac{1}{w}[-(py')' + qy]$。通过分部积分和边界条件，可以验证 $\langle \mathcal{L}y, z \rangle_w = \langle y, \mathcal{L}z \rangle_w$，即 $\mathcal{L}$ 是 $L^2_w$ 内积下的自伴算子。

    **正交性**：若 $\mathcal{L}y_m = \lambda_m y_m$ 和 $\mathcal{L}y_n = \lambda_n y_n$（$\lambda_m \neq \lambda_n$），则

    $$
    \lambda_m \langle y_m, y_n \rangle_w = \langle \mathcal{L}y_m, y_n \rangle_w = \langle y_m, \mathcal{L}y_n \rangle_w = \lambda_n \langle y_m, y_n \rangle_w.
    $$

    因此 $(\lambda_m - \lambda_n)\langle y_m, y_n \rangle_w = 0$，即 $\langle y_m, y_n \rangle_w = 0$。这与有限维对称矩阵不同特征值对应特征向量正交的证明完全平行。$\blacksquare$

!!! example "例 26.8"
    **热传导方程的分离变量。** 考虑一维热方程

    $$
    u_t = u_{xx}, \quad u(0, t) = u(\pi, t) = 0, \quad u(x, 0) = f(x).
    $$

    令 $u(x, t) = X(x)T(t)$，分离变量得 $T'/T = X''/X = -\lambda$。

    **空间问题**：$X'' = -\lambda X$，$X(0) = X(\pi) = 0$。这是 $p = w = 1$，$q = 0$ 的 Sturm-Liouville 问题。特征值 $\lambda_n = n^2$，特征函数 $X_n = \sin(nx)$（$n = 1, 2, 3, \ldots$）。

    **时间问题**：$T_n' = -n^2 T_n$，解为 $T_n = e^{-n^2 t}$。

    **完整解**：

    $$
    u(x, t) = \sum_{n=1}^{\infty} b_n e^{-n^2 t} \sin(nx), \quad b_n = \frac{2}{\pi}\int_0^{\pi} f(x)\sin(nx)\, dx.
    $$

    每个模态 $n$ 以速率 $n^2$ 衰减——高频分量衰减更快，这解释了热传导的平滑效应。特征值 $n^2$ 完全决定了扩散的时间尺度。

---

## 26.8 Floquet 理论

<div class="context-flow" markdown>

**问题**：$\mathbf{x}' = A(t)\mathbf{x}$，$A(t+T) = A(t)$ → 解不再是 $e^{At}$ → Floquet 定理：$\Phi(t) = P(t)e^{Bt}$（周期部分 $\times$ 指数部分）→ 稳定性由 $B$ 的特征值决定

**应用**：参数共振（Mathieu 方程）、晶格中的 Bloch 定理

</div>

当系数矩阵具有周期性 $A(t+T) = A(t)$ 时，Floquet 理论将周期系统化为常系数问题。

!!! definition "定义 26.11 (单值矩阵)"
    对周期系统 $\mathbf{x}' = A(t)\mathbf{x}$（$A(t+T) = A(t)$），设 $\Phi(t)$ 为满足 $\Phi(0) = I$ 的基本解矩阵。**单值矩阵**（monodromy matrix）定义为

    $$
    M = \Phi(T).
    $$

    $M$ 的特征值称为 **Floquet 乘子**（Floquet multipliers）。

!!! definition "定义 26.12 (Floquet 指数)"
    若 $\rho$ 为 Floquet 乘子，则 $\mu = \frac{1}{T}\ln \rho$（取适当分支）称为 **Floquet 指数**（Floquet exponent）。它们是等价常系数系统的"有效特征值"。

!!! theorem "定理 26.11 (Floquet 定理)"
    设 $A(t)$ 为连续的 $T$-周期矩阵。则基本解矩阵 $\Phi(t)$（$\Phi(0) = I$）可以写为

    $$
    \Phi(t) = P(t) e^{Bt},
    $$

    其中 $P(t+T) = P(t)$ 为 $T$-周期非奇异矩阵，$B$ 为常数矩阵，满足 $e^{BT} = M = \Phi(T)$。

??? proof "证明"
    **关键观察**：$\Psi(t) = \Phi(t+T)$ 也是基本解矩阵，因为

    $$
    \Psi'(t) = \Phi'(t+T) = A(t+T)\Phi(t+T) = A(t)\Psi(t).
    $$

    由基本解矩阵的唯一性（到右乘常数矩阵），$\Phi(t+T) = \Phi(t)M$。

    选取 $B$ 使得 $e^{BT} = M$（总存在，取矩阵对数）。定义 $P(t) = \Phi(t)e^{-Bt}$。则

    $$
    P(t+T) = \Phi(t+T)e^{-B(t+T)} = \Phi(t)Me^{-BT}e^{-Bt} = \Phi(t)e^{-Bt} = P(t).
    $$

    因此 $P(t)$ 是 $T$-周期的。$\blacksquare$

!!! theorem "定理 26.12 (Floquet 稳定性准则)"
    周期系统 $\mathbf{x}' = A(t)\mathbf{x}$ 的平衡点渐近稳定当且仅当所有 Floquet 乘子满足 $|\rho_i| < 1$（等价地，所有 Floquet 指数满足 $\operatorname{Re}(\mu_i) < 0$）。

??? proof "证明"
    由 $\Phi(nT) = M^n$ 和 $\Phi(t) = P(t)e^{Bt}$，$\|\Phi(nT)\| = \|M^n\|$。$M^n \to 0$ 当且仅当 $M$ 的所有特征值（即 Floquet 乘子）的模小于 $1$。在周期之间，$P(t)$ 有界，因此 $\|\Phi(t)\| \to 0$ 当且仅当 $|\rho_i| < 1$。$\blacksquare$

!!! example "例 26.9"
    **Mathieu 方程与参数共振。** Mathieu 方程

    $$
    y'' + (\delta + \varepsilon \cos 2t) y = 0
    $$

    化为系统 $\mathbf{x}' = A(t)\mathbf{x}$，$A(t) = \begin{pmatrix} 0 & 1 \\ -\delta - \varepsilon\cos 2t & 0 \end{pmatrix}$，周期 $T = \pi$。

    参数平面 $(\delta, \varepsilon)$ 中存在**稳定区**（$|\rho_i| \le 1$）和**不稳定区/共振舌**（$|\rho_i| > 1$）。当 $\varepsilon = 0$ 时，不稳定区退化为点 $\delta = n^2$（$n = 1, 2, \ldots$）。当 $\varepsilon > 0$ 时，不稳定区扩展为"舌"形区域。

    物理应用：倒置摆在适当频率和振幅的垂直振动下可以稳定（Kapitza 摆）——Floquet 分析给出精确的稳定条件。

!!! example "例 26.10"
    **Floquet 理论与 Bloch 定理。** 在固体物理中，晶格中电子的 Schrodinger 方程

    $$
    -\psi''(x) + V(x)\psi(x) = E\psi(x), \quad V(x+a) = V(x)
    $$

    是周期系数 ODE。Floquet 定理保证存在 Bloch 波形式的解

    $$
    \psi(x) = e^{ikx} u(x), \quad u(x+a) = u(x),
    $$

    其中 $k$ 为波数（Floquet 指数的虚部）。不同 $k$ 值对应的能量 $E(k)$ 构成**能带结构**，能带间隙对应 Floquet 乘子模为 $1$ 的边界。

## 练习题

1. **[一阶方程] 对于一阶常系数齐次线性微分方程组 $\mathbf{x}' = A\mathbf{x}$，其通解可以表示为什么形式？**
   ??? success "参考答案"
       通解可以表示为 $\mathbf{x}(t) = e^{At}\mathbf{c}$，其中 $\mathbf{c}$ 是由初始条件 $\mathbf{x}(0)$ 决定的常数向量，而 $e^{At}$ 是系统的基本解矩阵（状态转移矩阵）。

2. **[相图] 若 $2 \times 2$ 矩阵 $A$ 的特征值为一对纯虚数 $\pm i\omega$，此时系统的相图呈现什么形状？称为什么类型的奇点？**
   ??? success "参考答案"
       相图呈现出一系列同心的闭合椭圆轨迹。原点被称为**中心（Center）**，这代表一个无阻尼的简谐振荡系统。

3. **[稳定性] 一个非线性系统 $\mathbf{x}' = \mathbf{f}(\mathbf{x})$ 在其平衡点 $\mathbf{x}^*$ 处渐近稳定的充分条件是什么？**
   ??? success "参考答案"
       在 $\mathbf{x}^*$ 处计算系统的雅可比矩阵 $J = \nabla \mathbf{f}(\mathbf{x}^*)$。如果 $J$ 的所有特征值的实部都严格小于 0，则根据李雅普诺夫第一方法，系统在平衡点局部渐近稳定。

4. **[非齐次] 求解非齐次线性系统 $\mathbf{x}' = A\mathbf{x} + \mathbf{f}(t)$ 所使用的方法叫什么？其积分公式的形式是什么？**
   ??? success "参考答案"
       叫做常数变易法（Variation of Parameters）或杜哈梅原理（Duhamel's Principle）。其解为 $\mathbf{x}(t) = e^{At}\mathbf{x}(0) + \int_0^t e^{A(t-s)}\mathbf{f}(s) ds$。

5. **[高阶降维] 如何将一个高阶的标量微分方程 $y^{(n)} + a_{n-1}y^{(n-1)} + \cdots + a_0 y = 0$ 转化为一阶的线性方程组？**
   ??? success "参考答案"
       通过引入状态变量 $x_1 = y, x_2 = y', \dots, x_n = y^{(n-1)}$，原方程可化为 $\mathbf{x}' = C\mathbf{x}$，其中系数矩阵 $C$ 被称为友矩阵（Companion Matrix）。

6. **[解空间] 为什么 $n$ 阶常系数线性齐次微分方程的解集构成一个 $n$ 维向量空间？**
   ??? success "参考答案"
       因为微分算子是线性的，解的线性组合仍然是解（叠加原理）。且该解空间由 $n$ 个线性无关的基本解（如 $e^{\lambda_i t}$）张成，基础解系的存在性由特征方程保证。

7. **[矩阵指数] 当矩阵 $A$ 包含代数重数大于 1 且不可对角化的特征值时（即存在大小大于 1 的 Jordan 块），其微分方程组的解中会出现什么特殊的项？**
   ??? success "参考答案"
       会出现时间 $t$ 的多项式与指数函数的乘积，例如 $t e^{\lambda t}, t^2 e^{\lambda t}$。这对应于物理上的临界阻尼或多重共振态。

8. **[边值问题] Sturm-Liouville 问题的特征函数在什么内积下是正交的？**
   ??? success "参考答案"
       在带有特定权重函数 $w(x)$ 的函数空间内积 $\langle u, v \rangle = \int_a^b u(x)v(x)w(x) dx$ 下是正交的。

9. **[偏微分方程] 利用分离变量法求解一维热传导方程 $u_t = \alpha u_{xx}$ 时，空间部分和时间部分的解分别是什么形式？**
   ??? success "参考答案"
       空间部分产生一个 Sturm-Liouville 特征值问题，其解通常是由正弦和余弦构成的傅里叶基（空间谐波）；时间部分产生一阶常微分方程，其解是指数衰减函数 $e^{-\alpha \lambda_n t}$。

10. **[爱因斯坦思考题] 为什么牛顿第二定律（引力定律）和薛定谔方程等最基本的物理定律，都是以（偏）微分方程加线性代数的形式写成的？这反映了宇宙的什么哲学？**
    ??? success "参考答案"
        微分方程反映了宇宙运转的“局部性”与“因果性”——未来的状态仅仅平滑地取决于当前状态和当下的变化率，而没有超距的跳跃。而线性代数（特别是算子与本征值）反映了复杂现象背后的“叠加原理”。物理定律以这种形式呈现，表明宇宙看似混沌的演化，在无穷小的切空间底层，其实只是由几个简单的、可线性叠加的“本征模式”在各自独立地随着时间旋转或衰减。

## 本章小结

本章详细探讨了线性代数在常微分方程（ODE）和偏微分方程（PDE）中的核心应用，主要内容包括：

1. **一阶线性系统与矩阵指数**：确立了 $e^{At}$ 作为一阶常系数齐次系统 $\mathbf{x}'=A\mathbf{x}$ 统一解算子的地位，并推导了非齐次系统的杜哈梅积分公式。
2. **相图与局部稳定性**：通过分析 $2 \times 2$ 矩阵特征值的实虚部符号，对动力系统的平衡点（结点、焦点、鞍点、中心）进行了几何分类。对于非线性系统，雅可比矩阵在平衡点处的特征值决定了系统的局部渐近稳定性。
3. **高阶方程与友矩阵**：展示了如何通过状态变量的扩充，将任意高阶常系数标量方程降阶打击为一阶系统，建立了特征多项式与友矩阵特征值的对应关系。
4. **Sturm-Liouville 理论**：将有限维对称矩阵的谱定理推广到无穷维的微分算子，证明了在此框架下的特征函数必然构成一组完备的正交基。
5. **分离变量法**：利用微分算子的本征展开，将偏微分方程（如热传导和波动方程）巧妙地转化为无数个彼此独立的常微分方程，彰显了线性叠加原理在物理学中的无上威力。
6. **Floquet 理论**：处理了周期系数的微分系统，将参数共振和能带理论的机制归结为单值矩阵特征值（Floquet乘子）的分析。
