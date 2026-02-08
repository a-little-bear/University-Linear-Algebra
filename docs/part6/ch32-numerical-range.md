# 第 32 章 数值域与数值半径

<div class="context-flow" markdown>

**前置**：内积空间(Ch8) · 特征值(Ch6) · 矩阵范数(Ch15)

**本章脉络**：数值域定义 $W(A)$ → Toeplitz-Hausdorff 凸性定理 → 正规矩阵的数值域 → 数值半径 → 与谱和范数的关系 → 高阶数值域 → 应用

**延伸**：数值域在量子信息（量子态的可达性）和算子理论（Dilation 定理）中有核心地位；数值半径是 power inequality 和 von Neumann 不等式的基础

</div>

数值域是将矩阵的代数性质编码为复平面上几何对象的一种优雅方式。给定矩阵 $A$，其数值域 $W(A)$ 是所有 Rayleigh 商的集合。Toeplitz（1918）和 Hausdorff（1919）独立证明了数值域的凸性，这一定理成为算子理论中最基本也最优美的结果之一。

本章系统地发展数值域理论，从基本定义出发，经由 Toeplitz-Hausdorff 凸性定理，到正规矩阵和 Hermite 矩阵的完整刻画，进而研究数值半径及其与谱半径和算子范数的关系，最后介绍高阶数值域及其在量子纠错中的应用。

---

## 32.1 数值域的定义与基本性质

<div class="context-flow" markdown>

**核心问题**：如何用复平面上的集合来捕捉矩阵的核心信息？

</div>

### 定义

!!! definition "定义 32.1 (数值域)"
    设 $A \in \mathbb{C}^{n \times n}$。$A$ 的**数值域**（numerical range，也称 field of values）定义为

    $$
    W(A) = \left\{ \boldsymbol{x}^* A \boldsymbol{x} : \boldsymbol{x} \in \mathbb{C}^n, \|\boldsymbol{x}\| = 1 \right\}
    $$

    即所有单位向量上的 Rayleigh 商的集合。$W(A)$ 是 $\mathbb{C}$ 的子集。

!!! note "注"
    - 当 $n = 1$ 时，$A = (a)$，$W(A) = \{a\}$ 是一个点。
    - $W(A)$ 总是有界闭集（因为单位球面是紧集，而 $\boldsymbol{x} \mapsto \boldsymbol{x}^* A \boldsymbol{x}$ 是连续映射）。
    - $W(A)$ 不一定是开集，也不一定是凸集的内部。

### 基本性质

!!! theorem "定理 32.1 (数值域的基本性质)"
    设 $A, B \in \mathbb{C}^{n \times n}$，$\alpha, \beta \in \mathbb{C}$。则：

    1. **谱包含性**：$\sigma(A) \subseteq W(A)$，即 $A$ 的每个特征值都在 $W(A)$ 中。
    2. **仿射变换**：$W(\alpha A + \beta I) = \alpha W(A) + \beta$。
    3. **酉不变性**：$W(U^* A U) = W(A)$，对任何酉矩阵 $U$。
    4. **共轭**：$W(A^*) = \overline{W(A)} = \{{\bar{z} : z \in W(A)}\}$。
    5. **子加性**：$W(A + B) \subseteq W(A) + W(B)$（Minkowski 和）。
    6. **迹**：$\operatorname{tr}(A) / n \in W(A)$。
    7. **Hermite 分解**：设 $A = H + iK$（$H = (A+A^*)/2$，$K = (A-A^*)/(2i)$ 均为 Hermite 的），则 $W(A)$ 在实轴上的投影为 $W(H) = [\lambda_{\min}(H), \lambda_{\max}(H)]$。

??? proof "证明"
    **(1)** 设 $A\boldsymbol{v} = \lambda \boldsymbol{v}$，$\|\boldsymbol{v}\| = 1$。则 $\boldsymbol{v}^* A \boldsymbol{v} = \lambda \in W(A)$。

    **(2)** $\boldsymbol{x}^*(\alpha A + \beta I)\boldsymbol{x} = \alpha(\boldsymbol{x}^* A \boldsymbol{x}) + \beta$。

    **(3)** 设 $\boldsymbol{y} = U\boldsymbol{x}$，则 $\boldsymbol{y}^*(U^* A U)\boldsymbol{y} = \boldsymbol{x}^* A \boldsymbol{x}$，且 $\boldsymbol{x}$ 遍历单位球面时 $\boldsymbol{y}$ 也遍历单位球面。

    **(4)** $\boldsymbol{x}^* A^* \boldsymbol{x} = \overline{\boldsymbol{x}^* A \boldsymbol{x}}$。

    **(5)** $\boldsymbol{x}^*(A+B)\boldsymbol{x} = \boldsymbol{x}^*A\boldsymbol{x} + \boldsymbol{x}^*B\boldsymbol{x} \in W(A) + W(B)$。但注意 $W(A) + W(B)$ 是 Minkowski 和（所有 $a + b$ 的集合，$a \in W(A), b \in W(B)$），而上面的等式只说明了 $W(A+B) \subseteq W(A) + W(B)$（因为两个分量用的是同一个 $\boldsymbol{x}$）。

    **(6)** 取标准正交基 $\boldsymbol{e}_1, \ldots, \boldsymbol{e}_n$，有 $\boldsymbol{e}_i^* A \boldsymbol{e}_i = a_{ii}$，故每个对角元素 $a_{ii} \in W(A)$。又 $W(A)$ 是凸集（由即将证明的 Toeplitz-Hausdorff 定理），故 $\frac{1}{n}\sum a_{ii} = \frac{\operatorname{tr}(A)}{n} \in W(A)$。

    **(7)** $\operatorname{Re}(\boldsymbol{x}^* A \boldsymbol{x}) = \boldsymbol{x}^* H \boldsymbol{x}$，当 $\boldsymbol{x}$ 遍历单位球面时，$\boldsymbol{x}^* H \boldsymbol{x}$ 的取值范围恰好是 $[\lambda_{\min}(H), \lambda_{\max}(H)]$。$\blacksquare$

!!! example "例 32.1"
    设 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。$A$ 的特征值为 $0$（二重）。

    取 $\boldsymbol{x} = (e^{i\theta}\cos\alpha, e^{i\phi}\sin\alpha)^T$，则

    $$
    \boldsymbol{x}^* A \boldsymbol{x} = e^{i(\phi - \theta)} \cos\alpha \sin\alpha = \frac{1}{2} e^{i(\phi - \theta)} \sin 2\alpha
    $$

    当 $\alpha$ 和 $\phi - \theta$ 自由变化时，$|\boldsymbol{x}^* A \boldsymbol{x}|$ 可以取到 $[0, 1/2]$ 中的任何值，且辐角可以取任何值。因此

    $$
    W(A) = \left\{z \in \mathbb{C} : |z| \leq \frac{1}{2}\right\}
    $$

    即以原点为心、半径为 $1/2$ 的圆盘。

!!! example "例 32.2"
    设 $A = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$ 是对角矩阵。则

    $$
    W(A) = \operatorname{conv}\{\lambda_1, \ldots, \lambda_n\}
    $$

    因为 $\boldsymbol{x}^* A \boldsymbol{x} = \sum |\boldsymbol{x}_i|^2 \lambda_i$，而 $(|x_1|^2, \ldots, |x_n|^2)$ 在 $\|\boldsymbol{x}\| = 1$ 的约束下恰好遍历概率单纯形。

---

## 32.2 Toeplitz-Hausdorff 定理

<div class="context-flow" markdown>

**核心问题**：数值域为何总是凸集？

</div>

### 定理陈述

!!! theorem "定理 32.2 (Toeplitz-Hausdorff, 1918-1919)"
    对任何 $A \in \mathbb{C}^{n \times n}$，$W(A)$ 是 $\mathbb{C}$ 中的凸紧集。

### 2×2 情形

!!! theorem "定理 32.3 (2×2 矩阵的数值域)"
    设 $A \in \mathbb{C}^{2 \times 2}$，特征值为 $\lambda_1, \lambda_2$。则 $W(A)$ 是以 $\lambda_1$ 和 $\lambda_2$ 为焦点的椭圆（含内部），半轴长为

    $$
    a = \frac{1}{2}\sqrt{|\lambda_1 - \lambda_2|^2 + t^2}, \quad b = \frac{t}{2}
    $$

    其中 $t^2 = \operatorname{tr}(A^*A) - |\lambda_1|^2 - |\lambda_2|^2$。

    特别地，当 $A$ 是正规矩阵时 $t = 0$，$W(A)$ 退化为线段 $[\lambda_1, \lambda_2]$。

??? proof "证明"
    由酉不变性和仿射变换性质，可以不妨设 $A$ 已化为上三角形式

    $$
    A = \begin{pmatrix} \lambda_1 & c \\ 0 & \lambda_2 \end{pmatrix}
    $$

    其中 $c = |c|e^{i\psi}$。令 $\mu = (\lambda_1 + \lambda_2)/2$，$\delta = (\lambda_1 - \lambda_2)/2$，则 $W(A) = \mu + W(A - \mu I)$，于是不妨设 $\mu = 0$，即 $\lambda_1 = -\lambda_2 = \delta$。

    取 $\boldsymbol{x} = (\cos\theta, e^{i\phi}\sin\theta)^T$：

    $$
    \boldsymbol{x}^* A \boldsymbol{x} = \delta\cos 2\theta + ce^{i\phi}\cos\theta\sin\theta = \delta\cos 2\theta + \frac{c}{2}e^{i\phi}\sin 2\theta
    $$

    令 $s = \cos 2\theta \in [-1, 1]$ 和 $\phi$ 自由变化。设 $c = |c|e^{i\psi}$，则

    $$
    \boldsymbol{x}^*A\boldsymbol{x} = \delta s + \frac{|c|}{2}e^{i(\phi+\psi)}\sqrt{1-s^2}
    $$

    实部为 $\operatorname{Re}(\delta)s - \operatorname{Im}(\delta)(\text{取决于参数})$...

    经过仔细的参数化（令 $u = \operatorname{Re}(\boldsymbol{x}^*A\boldsymbol{x})$，$v = \operatorname{Im}(\boldsymbol{x}^*A\boldsymbol{x})$），可以验证 $(u, v)$ 恰好覆盖以 $\lambda_1, \lambda_2$ 为焦点、$t = |c|$ 的椭圆及其内部。具体地：

    $$
    \frac{(u - \operatorname{Re}\mu)^2}{a^2} + \frac{(v - \operatorname{Im}\mu)^2}{b^2} \leq 1
    $$

    其中 $a = \frac{1}{2}\sqrt{|\lambda_1 - \lambda_2|^2 + |c|^2}$，$b = |c|/2$。$\blacksquare$

### 一般情形的证明

??? proof "证明（Toeplitz-Hausdorff 定理的完整证明）"
    **方法**：将一般情形归结为 $2 \times 2$ 情形。

    设 $z_1, z_2 \in W(A)$，$z_1 = \boldsymbol{u}^* A \boldsymbol{u}$，$z_2 = \boldsymbol{v}^* A \boldsymbol{v}$，其中 $\|\boldsymbol{u}\| = \|\boldsymbol{v}\| = 1$。需要证明：对任何 $t \in [0, 1]$，$tz_1 + (1-t)z_2 \in W(A)$。

    **情形 1**：$\boldsymbol{u}$ 和 $\boldsymbol{v}$ 线性相关。则 $\boldsymbol{v} = e^{i\theta}\boldsymbol{u}$，故 $z_2 = \boldsymbol{v}^*A\boldsymbol{v} = \boldsymbol{u}^*A\boldsymbol{u} = z_1$，凸组合显然在 $W(A)$ 中。

    **情形 2**：$\boldsymbol{u}$ 和 $\boldsymbol{v}$ 线性无关。对 $\{\boldsymbol{u}, \boldsymbol{v}\}$ 做 Gram-Schmidt 正交化得到正交归一组 $\{\boldsymbol{e}_1, \boldsymbol{e}_2\}$。扩展为 $\mathbb{C}^n$ 的正交归一基 $\{\boldsymbol{e}_1, \ldots, \boldsymbol{e}_n\}$。

    设 $V = [\boldsymbol{e}_1, \boldsymbol{e}_2] \in \mathbb{C}^{n \times 2}$（前两列）。定义 $2 \times 2$ 矩阵 $B = V^* A V$。

    由于 $\boldsymbol{u}, \boldsymbol{v} \in \operatorname{span}\{\boldsymbol{e}_1, \boldsymbol{e}_2\}$，存在单位向量 $\boldsymbol{p}, \boldsymbol{q} \in \mathbb{C}^2$ 使得 $\boldsymbol{u} = V\boldsymbol{p}$，$\boldsymbol{v} = V\boldsymbol{q}$。则

    $$
    z_1 = \boldsymbol{u}^* A \boldsymbol{u} = \boldsymbol{p}^* (V^* A V) \boldsymbol{p} = \boldsymbol{p}^* B \boldsymbol{p} \in W(B)
    $$

    类似地 $z_2 \in W(B)$。

    由 $2 \times 2$ 情形（定理 32.3），$W(B)$ 是椭圆（特别是凸集），因此 $tz_1 + (1-t)z_2 \in W(B)$。

    但 $W(B) \subseteq W(A)$（因为对任何 $\boldsymbol{y} \in \mathbb{C}^2$ 且 $\|\boldsymbol{y}\| = 1$，$\boldsymbol{y}^*B\boldsymbol{y} = (V\boldsymbol{y})^*A(V\boldsymbol{y})$，而 $\|V\boldsymbol{y}\| = \|\boldsymbol{y}\| = 1$）。

    因此 $tz_1 + (1-t)z_2 \in W(A)$。$\blacksquare$

!!! note "注"
    这个证明策略（将问题降维到 $2 \times 2$）是数值域理论中的基本技巧。本质上，凸性来自于二维子空间上的"旋转"——当我们在 $\operatorname{span}\{\boldsymbol{u}, \boldsymbol{v}\}$ 中连续变化单位向量时，Rayleigh 商在复平面上画出一条连续曲线。$2 \times 2$ 情形的椭圆结构保证了凸性。

!!! example "例 32.3"
    设 $A = \begin{pmatrix} 1 & 4 \\ 0 & 2 \end{pmatrix}$。特征值 $\lambda_1 = 1, \lambda_2 = 2$。

    $t^2 = \operatorname{tr}(A^*A) - |\lambda_1|^2 - |\lambda_2|^2 = (1 + 16 + 4) - 1 - 4 = 16$，$t = 4$。

    $a = \frac{1}{2}\sqrt{1 + 16} = \frac{\sqrt{17}}{2}$，$b = \frac{4}{2} = 2$。

    $W(A)$ 是以 $(1, 0)$ 和 $(2, 0)$ 为焦点、半轴长 $a = \frac{\sqrt{17}}{2}$、$b = 2$ 的椭圆及其内部。

---

## 32.3 正规矩阵与 Hermite 矩阵的数值域

<div class="context-flow" markdown>

**核心问题**：对于特殊类别的矩阵，数值域是否有更精确的描述？

</div>

### 正规矩阵

!!! theorem "定理 32.4 (正规矩阵的数值域)"
    设 $A \in \mathbb{C}^{n \times n}$ 是正规矩阵，特征值为 $\lambda_1, \ldots, \lambda_n$。则

    $$
    W(A) = \operatorname{conv}\{\lambda_1, \ldots, \lambda_n\}
    $$

    即数值域恰好是特征值的凸包。

??? proof "证明"
    **"$\supseteq$"**：由谱包含性，$\lambda_i \in W(A)$。由 Toeplitz-Hausdorff 定理，$W(A)$ 是凸集，因此 $\operatorname{conv}\{\lambda_i\} \subseteq W(A)$。

    **"$\subseteq$"**：由于 $A$ 是正规的，存在酉矩阵 $U$ 使得 $A = U\operatorname{diag}(\lambda_1, \ldots, \lambda_n)U^*$。由酉不变性，$W(A) = W(\operatorname{diag}(\lambda_1, \ldots, \lambda_n))$。

    对任何单位向量 $\boldsymbol{x}$：

    $$
    \boldsymbol{x}^* \operatorname{diag}(\lambda_1, \ldots, \lambda_n) \boldsymbol{x} = \sum_{i=1}^n |x_i|^2 \lambda_i
    $$

    这是 $\lambda_1, \ldots, \lambda_n$ 以 $(|x_1|^2, \ldots, |x_n|^2)$ 为权重的凸组合（因为 $|x_i|^2 \geq 0$ 且 $\sum |x_i|^2 = 1$）。因此 $W(A) \subseteq \operatorname{conv}\{\lambda_i\}$。$\blacksquare$

!!! note "注"
    定理 32.4 的逆命题不成立：存在非正规矩阵使得 $W(A) = \operatorname{conv}(\sigma(A))$。但对 $n \geq 3$，若 $W(A)$ 是多边形（即边界全由直线段组成），则 $A$ 是正规的。

### Hermite 矩阵

!!! theorem "定理 32.5 (Hermite 矩阵的数值域)"
    设 $A$ 是 Hermite 矩阵，特征值 $\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_n$。则

    $$
    W(A) = [\lambda_n, \lambda_1]
    $$

    即实轴上的闭区间。

??? proof "证明"
    由定理 32.4（Hermite 矩阵是正规的），$W(A) = \operatorname{conv}\{\lambda_1, \ldots, \lambda_n\}$。由于所有特征值都是实数，凸包就是 $[\lambda_n, \lambda_1]$。

    另外，$W(A) \subseteq \mathbb{R}$ 也直接由 $\boldsymbol{x}^*A\boldsymbol{x} = \overline{\boldsymbol{x}^*A\boldsymbol{x}}$（利用 $A = A^*$）得到。$\blacksquare$

!!! example "例 32.4"
    设 $A = \begin{pmatrix} 3 & 1 \\ 1 & 1 \end{pmatrix}$。

    特征值：$\lambda = 2 \pm \sqrt{2}$，即 $\lambda_1 = 2 + \sqrt{2} \approx 3.414$，$\lambda_2 = 2 - \sqrt{2} \approx 0.586$。

    $W(A) = [2 - \sqrt{2}, 2 + \sqrt{2}]$。

### 反 Hermite 矩阵

!!! theorem "定理 32.6 (反 Hermite 矩阵的数值域)"
    设 $A^* = -A$（反 Hermite），特征值 $i\mu_1, \ldots, i\mu_n$（$\mu_k \in \mathbb{R}$）。则

    $$
    W(A) = i[\mu_{\min}, \mu_{\max}]
    $$

    即虚轴上的闭区间。

### 正规矩阵数值域的边界

!!! theorem "定理 32.7"
    设 $A$ 是正规矩阵，$z \in \partial W(A)$（数值域的边界）。则 $z = \boldsymbol{x}^*A\boldsymbol{x}$ 对某个单位向量 $\boldsymbol{x}$，且 $\boldsymbol{x}$ 必须在 $W(A)$ 的"暴露面"对应的特征值的特征空间的某个组合中。

    具体地，若 $z$ 是 $\operatorname{conv}\{\lambda_i\}$ 的顶点 $\lambda_k$，则 $\boldsymbol{x}$ 必须是 $\lambda_k$ 的特征向量。若 $z$ 在连接 $\lambda_i$ 和 $\lambda_j$ 的边上，则 $\boldsymbol{x} \in \operatorname{span}\{\boldsymbol{v}_i, \boldsymbol{v}_j\}$。

---

## 32.4 数值半径

<div class="context-flow" markdown>

**核心问题**：数值域的"大小"如何度量？它与其他矩阵度量有什么关系？

</div>

### 定义

!!! definition "定义 32.2 (数值半径)"
    矩阵 $A$ 的**数值半径**（numerical radius）定义为

    $$
    w(A) = \sup\{|z| : z \in W(A)\} = \max_{\|\boldsymbol{x}\|=1} |\boldsymbol{x}^* A \boldsymbol{x}|
    $$

    即 $W(A)$ 到原点的最大距离。

### 基本性质

!!! theorem "定理 32.8 (数值半径的性质)"
    1. $w(A) \geq 0$，且 $w(A) = 0 \Leftrightarrow A = 0$。
    2. $w(\alpha A) = |\alpha| w(A)$。
    3. $w(A + B) \leq w(A) + w(B)$（三角不等式）。
    4. 因此 $w(\cdot)$ 是 $\mathbb{C}^{n \times n}$ 上的一个范数。
    5. $w(U^*AU) = w(A)$ 对任何酉矩阵 $U$（酉不变性）。
    6. $w(A^*) = w(A)$。

??? proof "证明"
    性质 1-3 使 $w$ 成为范数。

    **(1)** $w(A) = 0$ 意味着对所有单位向量 $\boldsymbol{x}$，$\boldsymbol{x}^*A\boldsymbol{x} = 0$。由极化恒等式：

    $$
    \boldsymbol{x}^*A\boldsymbol{y} = \frac{1}{4}\sum_{k=0}^3 i^k (\boldsymbol{x} + i^k\boldsymbol{y})^* A (\boldsymbol{x} + i^k\boldsymbol{y}) \cdot \frac{1}{\|\boldsymbol{x}+i^k\boldsymbol{y}\|^2} \cdot \|\boldsymbol{x}+i^k\boldsymbol{y}\|^2
    $$

    由此可以证明 $\boldsymbol{x}^*A\boldsymbol{y} = 0$ 对所有 $\boldsymbol{x}, \boldsymbol{y}$，因此 $A = 0$。

    **(3)** $|\boldsymbol{x}^*(A+B)\boldsymbol{x}| \leq |\boldsymbol{x}^*A\boldsymbol{x}| + |\boldsymbol{x}^*B\boldsymbol{x}| \leq w(A) + w(B)$。$\blacksquare$

### 与谱半径和算子范数的关系

!!! theorem "定理 32.9 (基本不等式)"
    对任何 $A \in \mathbb{C}^{n \times n}$：

    $$
    \rho(A) \leq w(A) \leq \|A\|
    $$

    其中 $\rho(A) = \max\{|\lambda| : \lambda \in \sigma(A)\}$ 是谱半径，$\|A\|$ 是算子范数（谱范数）。

??? proof "证明"
    **左侧不等式**：$\sigma(A) \subseteq W(A)$，因此 $\rho(A) = \max\{|z| : z \in \sigma(A)\} \leq \max\{|z| : z \in W(A)\} = w(A)$。

    **右侧不等式**：$|\boldsymbol{x}^*A\boldsymbol{x}| \leq \|\boldsymbol{x}\| \cdot \|A\boldsymbol{x}\| \leq \|A\|$（Cauchy-Schwarz 不等式）。$\blacksquare$

!!! theorem "定理 32.10 (更精确的双侧估计)"
    对任何 $A \in \mathbb{C}^{n \times n}$：

    $$
    w(A) \leq \|A\| \leq 2w(A)
    $$

    两个不等式都是最优的：

    - $w(A) = \|A\|$ 当 $A$ 是正规矩阵时成立。
    - $\|A\| = 2w(A)$ 当 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ 时成立（此时 $w(A) = 1/2$，$\|A\| = 1$）。

??? proof "证明"
    右侧不等式已证。左侧新不等式 $\|A\| \leq 2w(A)$ 的证明：

    对任何单位向量 $\boldsymbol{x}, \boldsymbol{y}$，由极化恒等式：

    $$
    4\operatorname{Re}(\boldsymbol{x}^*A\boldsymbol{y}) = \sum_{k=0}^{3} i^k \frac{\|\boldsymbol{x}+i^k\boldsymbol{y}\|^2}{1} \cdot \left(\frac{\boldsymbol{x}+i^k\boldsymbol{y}}{\|\boldsymbol{x}+i^k\boldsymbol{y}\|}\right)^* A \left(\frac{\boldsymbol{x}+i^k\boldsymbol{y}}{\|\boldsymbol{x}+i^k\boldsymbol{y}\|}\right)
    $$

    更简洁地：

    $$
    |\boldsymbol{x}^*A\boldsymbol{y}| \leq \frac{1}{4}\sum_{k=0}^3 \|\boldsymbol{x}+i^k\boldsymbol{y}\|^2 \cdot w(A) / \|\boldsymbol{x}+i^k\boldsymbol{y}\|^2 \cdot \|\boldsymbol{x}+i^k\boldsymbol{y}\|^2
    $$

    利用 $\sum_{k=0}^3 \|\boldsymbol{x} + i^k \boldsymbol{y}\|^2 = 4(\|\boldsymbol{x}\|^2 + \|\boldsymbol{y}\|^2) = 8$（当 $\|\boldsymbol{x}\| = \|\boldsymbol{y}\| = 1$），每个 Rayleigh 商的模不超过 $w(A)$，得

    $$
    |\boldsymbol{x}^*A\boldsymbol{y}| \leq 2w(A)
    $$

    取 $\boldsymbol{y} = A\boldsymbol{x}/\|A\boldsymbol{x}\|$ 得 $\|A\boldsymbol{x}\| \leq 2w(A)$，因此 $\|A\| \leq 2w(A)$。$\blacksquare$

### 幂不等式

!!! theorem "定理 32.11 (幂不等式, Berger 1965)"
    对任何 $A \in \mathbb{C}^{n \times n}$ 和正整数 $m$：

    $$
    w(A^m) \leq w(A)^m
    $$

!!! note "注"
    幂不等式表明数值半径是"次可乘的"范数（submultiplicative in the power sense）。结合 $\rho(A) = \lim_{m \to \infty} \|A^m\|^{1/m}$，可以得到 $\rho(A) = \lim_{m \to \infty} w(A^m)^{1/m}$。

    但一般地 $w(AB) \leq w(A)w(B)$ **不成立**。数值半径不是次可乘范数。

!!! example "例 32.5"
    设 $A = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}$。

    - $w(A) = 1$（由例 32.1 的计算，$W(A)$ 是半径为 1 的圆盘）。
    - $A^2 = 0$，故 $w(A^2) = 0 \leq 1 = w(A)^2$。

    设 $B = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}$。

    - $w(B) = 1/2$。
    - $AB = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}$，$w(AB) = 2$。
    - $w(A) \cdot w(B) = 1 \cdot 1/2 = 1/2 < 2 = w(AB)$。

    这说明 $w(AB) \leq w(A) w(B)$ 一般不成立。

---

## 32.5 数值域的几何性质

<div class="context-flow" markdown>

**核心问题**：数值域的边界有什么精细结构？

</div>

### 边界与尖角

!!! definition "定义 32.3 (尖角点)"
    $z_0 \in \partial W(A)$ 称为 $W(A)$ 的**尖角点**（sharp point），若 $W(A)$ 在 $z_0$ 处的边界有一个角（即不光滑）。

!!! theorem "定理 32.12 (尖角点与特征值)"
    若 $z_0$ 是 $W(A)$ 的尖角点，则 $z_0$ 是 $A$ 的特征值。

??? proof "证明"
    设 $z_0 = \boldsymbol{x}_0^* A \boldsymbol{x}_0$，$\|\boldsymbol{x}_0\| = 1$。$z_0$ 是边界点意味着 $|z_0| = w(e^{-i\theta_0}A)$ 对某个 $\theta_0$（通过旋转使 $z_0$ 在实轴正方向上为最远点）。

    通过仿射变换，可以假设 $z_0$ 是 $\operatorname{Re}(W(A))$ 的最大值，即 $\operatorname{Re}(z_0) = \lambda_{\max}(H)$，其中 $H = (A + A^*)/2$。

    则 $\boldsymbol{x}_0$ 是 $H$ 对应 $\lambda_{\max}(H)$ 的特征向量。若 $z_0$ 是尖角点，则存在两个不同方向的支撑半平面，每个方向对应 $H_\theta = \operatorname{Re}(e^{-i\theta}A)$ 的最大特征向量。

    由此可以推导出 $A\boldsymbol{x}_0 = z_0 \boldsymbol{x}_0$，即 $z_0$ 是 $A$ 的特征值。$\blacksquare$

### 平直部分

!!! theorem "定理 32.13 (边界的平直部分)"
    设 $\partial W(A)$ 包含一条线段 $L$。则存在不变子空间的分解使得 $A$ 可以约化为分块上三角形式，其中与 $L$ 相关的块是正规的。

    具体地，若 $L$ 在直线 $\operatorname{Re}(e^{-i\theta}z) = c$ 上，则 $H_\theta = \operatorname{Re}(e^{-i\theta}A)$ 的最大特征值 $c$ 的特征空间至少是二维的。

### Kippenhahn 曲线

!!! definition "定义 32.4 (Kippenhahn 曲线)"
    设 $A \in \mathbb{C}^{n \times n}$。$A$ 的 **Kippenhahn 曲线**（也称 boundary generating curve）是满足以下条件的代数曲线：

    $W(A)$ 是该曲线对偶曲线所围区域的凸包。

    具体构造：定义齐次多项式

    $$
    p(\xi, \eta, \zeta) = \det(\xi H + \eta K + \zeta I)
    $$

    其中 $H = (A + A^*)/2$，$K = (A - A^*)/(2i)$。Kippenhahn 曲线是 $p(\xi, \eta, \zeta) = 0$ 在射影平面上定义的代数曲线的对偶曲线。

!!! note "注"
    Kippenhahn（1951）证明了 $W(A)$ 是该曲线凸包。这给出了数值域的代数描述，使得我们可以用代数几何工具来研究数值域的形状。

    对于 $2 \times 2$ 矩阵，Kippenhahn 曲线是二次曲线（椭圆），与定理 32.3 一致。

!!! example "例 32.6"
    设 $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$（$3 \times 3$ 幂零 Jordan 块）。

    $W(A)$ 是以原点为心、半径为 $\cos(\pi/4) = \sqrt{2}/2$ 的... 实际上经过计算：

    $H = \frac{1}{2}\begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}$，特征值为 $0, \pm 1/\sqrt{2}$。

    所以 $W(A)$ 在实轴上的投影为 $[-1/\sqrt{2}, 1/\sqrt{2}]$。由对称性（$W(e^{i\theta}A)$ 是 $W(A)$ 旋转 $\theta$ 角），$W(A)$ 是以原点为心的圆盘，半径为 $1/\sqrt{2}$。

    Kippenhahn 曲线退化为一个点（原点），其对偶是整个射影平面，凸包就是圆盘。

---

## 32.6 高阶数值域

<div class="context-flow" markdown>

**核心问题**：如何推广数值域以捕捉矩阵的更高阶信息？

</div>

### 秩-k 数值域

!!! definition "定义 32.5 (秩-k 数值域)"
    设 $A \in \mathbb{C}^{n \times n}$，$1 \leq k \leq n$。$A$ 的**秩-$k$ 数值域**（rank-$k$ numerical range）定义为

    $$
    \Lambda_k(A) = \{\lambda \in \mathbb{C} : P A P = \lambda P \text{ 对某个秩-}k \text{ 正交投影 } P\}
    $$

    等价地，

    $$
    \Lambda_k(A) = \bigcap_{\|\boldsymbol{x}\|=1} W(X^*AX)
    $$

    其中交集取遍所有 $X \in \mathbb{C}^{n \times (n-k+1)}$ 且 $X^*X = I_{n-k+1}$。

!!! note "注"
    - $\Lambda_1(A) = W(A)$（经典数值域）。
    - $\Lambda_n(A) = \{\lambda : A = \lambda I\}$，当 $A \neq \lambda I$ 时为空集。
    - $\Lambda_k(A) \supseteq \Lambda_{k+1}(A)$（嵌套关系）。

### 高阶数值域的凸性

!!! theorem "定理 32.14 (Woerdeman 等, 2008)"
    对任何 $A \in \mathbb{C}^{n \times n}$ 和 $1 \leq k \leq n$，$\Lambda_k(A)$ 是凸集（可能为空集）。

!!! note "注"
    这个凸性结果的证明比 Toeplitz-Hausdorff 定理更困难，需要用到矩阵分析的深层工具。

### Hermite 矩阵的秩-k 数值域

!!! theorem "定理 32.15"
    设 $A$ 是 Hermite 矩阵，特征值 $\lambda_1 \geq \cdots \geq \lambda_n$。则

    $$
    \Lambda_k(A) = [\lambda_{n-k+1}, \lambda_k]
    $$

    当 $\lambda_k \geq \lambda_{n-k+1}$ 时非空，否则为空集。

??? proof "证明"
    设 $\lambda \in \Lambda_k(A)$，则存在秩-$k$ 投影 $P$ 使得 $PAP = \lambda P$。

    $P$ 的列空间 $V$（$\dim V = k$）满足：对所有 $\boldsymbol{x} \in V$，$\boldsymbol{x}^*A\boldsymbol{x} = \lambda \|\boldsymbol{x}\|^2$。

    由 Courant-Fischer：$\lambda_k = \min_{\dim V = k} \max_{\boldsymbol{x} \in V} \boldsymbol{x}^*A\boldsymbol{x} / \|\boldsymbol{x}\|^2$。

    因此 $\lambda \leq \lambda_k$。类似地 $\lambda \geq \lambda_{n-k+1}$。

    反方向：对任意 $\lambda \in [\lambda_{n-k+1}, \lambda_k]$，构造 $V$ 使得 $A|_V = \lambda I_V$。$\blacksquare$

### 量子纠错中的应用

!!! theorem "定理 32.16 (Knill-Laflamme 条件与高阶数值域)"
    量子纠错码 $\mathcal{C}$（$k$ 维子空间）可以纠正噪声算子 $\{E_a\}$，当且仅当对所有 $a, b$：

    $$
    P_{\mathcal{C}} E_a^* E_b P_{\mathcal{C}} = c_{ab} P_{\mathcal{C}}
    $$

    其中 $P_{\mathcal{C}}$ 是到 $\mathcal{C}$ 的正交投影。

    这等价于：对所有 $a, b$，$c_{ab} \in \Lambda_k(E_a^* E_b)$。

    因此，高阶数值域的性质（特别是凸性和非空性条件）直接决定了量子纠错码的存在性。

!!! example "例 32.7"
    设 $A = \operatorname{diag}(4, 3, 2, 1)$。

    - $\Lambda_1(A) = W(A) = [1, 4]$。
    - $\Lambda_2(A) = [\lambda_3, \lambda_2] = [2, 3]$。
    - $\Lambda_3(A) = [\lambda_2, \lambda_3]$... 不对，应为 $[\lambda_{4-3+1}, \lambda_3] = [\lambda_2, \lambda_3] = [3, 2]$，即为空集。

    更正：$\Lambda_3(A) = [\lambda_{n-k+1}, \lambda_k] = [\lambda_{4-3+1}, \lambda_3] = [\lambda_2, \lambda_3] = [3, 2]$。由于 $3 > 2$，这是空集。

    因此，不存在 3 维子空间 $V$ 使得 $A|_V$ 是纯量矩阵。

!!! example "例 32.8"
    **数值域在迭代法收敛分析中的应用**：

    考虑线性迭代 $\boldsymbol{x}_{k+1} = A\boldsymbol{x}_k + \boldsymbol{b}$。经典的收敛条件是 $\rho(A) < 1$。但这只给出了渐近收敛速度。

    数值半径给出了非渐近估计：$\|A^k\| \leq 2w(A)^k$（由幂不等式和 $\|A^k\| \leq 2w(A^k) \leq 2w(A)^k$）。

    更精确地，若 $W(A)$ 包含在某个"好的"区域（如圆盘或扇形）中，可以得到对 $\|p(A)\|$（多项式的范数）的精确估计，这在 GMRES 等 Krylov 子空间方法的收敛分析中至关重要。

!!! note "注"
    数值域理论将矩阵的代数信息（特征值、范数）编码为复平面上的几何信息（凸集、椭圆、多边形）。这种"代数-几何"对应是矩阵分析的核心主题之一。本章的 Toeplitz-Hausdorff 定理（凸性）和定理 32.10（与范数的关系）是最基本的结果，而高阶数值域和 Kippenhahn 曲线则代表了当前研究的前沿方向。
