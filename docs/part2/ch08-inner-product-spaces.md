# 第 8 章 内积空间

<div class="context-flow" markdown>

**前置**：Ch7 $\mathbb{R}^n$ 的正交性

**本章脉络**：内积公理 → 范数/Cauchy-Schwarz → 正交补/投影 → 伴随算子 → 正规/自伴算子 → **谱定理**

**延伸**：内积空间推广到无穷维得到 Hilbert 空间（$L^2$ 空间、Sobolev 空间）；在量子力学中态空间是 Hilbert 空间；Fourier 分析本质上是 $L^2$ 空间中的正交展开；再生核 Hilbert 空间（RKHS）是核方法与机器学习的理论基础

</div>

内积空间（inner product space）是线性代数中最富有几何直觉的结构之一。通过在向量空间上引入内积运算，我们可以严格地定义长度、角度、正交性等几何概念，从而将代数结构与几何直觉深度融合。本章从内积的公理化定义出发，系统地建立范数与距离理论，研究正交性与正交分解，深入探讨伴随算子、正规算子和自伴算子的性质，最终推导出线性代数中最深刻的结果之一——谱定理。

---

## 8.1 内积的定义

<div class="context-flow" markdown>

将 $\mathbb{R}^n$ 的点积 $\mathbf{x}^T\mathbf{y}$ **公理化**：正定性 + 对称性（或共轭对称性）+ 线性性 → 适用于函数空间等一般空间

</div>

在欧几里得空间 $\mathbb{R}^n$ 中，点积 $\mathbf{x} \cdot \mathbf{y} = \sum_{i=1}^n x_i y_i$ 赋予了向量以长度和角度的概念。内积是点积的推广，它将这种几何结构引入到一般的向量空间中。

### 实内积空间

!!! definition "定义 8.1 (实内积)"
    设 $V$ 是 $\mathbb{R}$ 上的向量空间。$V$ 上的一个**实内积**（real inner product）是一个映射 $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{R}$，满足以下四条公理：对任意 $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ 和 $\alpha \in \mathbb{R}$，

    1. **正定性**（positive definiteness）：$\langle \mathbf{v}, \mathbf{v} \rangle \geq 0$，且 $\langle \mathbf{v}, \mathbf{v} \rangle = 0$ 当且仅当 $\mathbf{v} = \mathbf{0}$；
    2. **对称性**（symmetry）：$\langle \mathbf{u}, \mathbf{v} \rangle = \langle \mathbf{v}, \mathbf{u} \rangle$；
    3. **第一变元的齐次性**（homogeneity）：$\langle \alpha\mathbf{u}, \mathbf{v} \rangle = \alpha \langle \mathbf{u}, \mathbf{v} \rangle$；
    4. **第一变元的可加性**（additivity）：$\langle \mathbf{u} + \mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{u}, \mathbf{w} \rangle + \langle \mathbf{v}, \mathbf{w} \rangle$。

    配备了实内积的实向量空间 $(V, \langle \cdot, \cdot \rangle)$ 称为**实内积空间**。

!!! note "注"
    公理 3 和 4 合在一起等价于内积关于第一变元是线性的：$\langle \alpha\mathbf{u} + \beta\mathbf{v}, \mathbf{w} \rangle = \alpha\langle \mathbf{u}, \mathbf{w} \rangle + \beta\langle \mathbf{v}, \mathbf{w} \rangle$。由对称性可知，实内积关于第二变元也是线性的，即实内积是**双线性**的（bilinear）。

### 复内积空间

!!! definition "定义 8.2 (复内积)"
    设 $V$ 是 $\mathbb{C}$ 上的向量空间。$V$ 上的一个**复内积**（complex inner product，或 Hermitian 内积）是一个映射 $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{C}$，满足以下四条公理：对任意 $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ 和 $\alpha \in \mathbb{C}$，

    1. **正定性**：$\langle \mathbf{v}, \mathbf{v} \rangle \geq 0$（即 $\langle \mathbf{v}, \mathbf{v} \rangle \in \mathbb{R}$ 且非负），且 $\langle \mathbf{v}, \mathbf{v} \rangle = 0$ 当且仅当 $\mathbf{v} = \mathbf{0}$；
    2. **共轭对称性**（conjugate symmetry）：$\langle \mathbf{u}, \mathbf{v} \rangle = \overline{\langle \mathbf{v}, \mathbf{u} \rangle}$；
    3. **第一变元的齐次性**：$\langle \alpha\mathbf{u}, \mathbf{v} \rangle = \alpha \langle \mathbf{u}, \mathbf{v} \rangle$；
    4. **第一变元的可加性**：$\langle \mathbf{u} + \mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{u}, \mathbf{w} \rangle + \langle \mathbf{v}, \mathbf{w} \rangle$。

    配备了复内积的复向量空间 $(V, \langle \cdot, \cdot \rangle)$ 称为**复内积空间**。

!!! note "注"
    在复内积空间中，内积关于第一变元是线性的，但关于第二变元是**共轭线性**的（conjugate linear）：$\langle \mathbf{u}, \alpha\mathbf{v} \rangle = \bar{\alpha}\langle \mathbf{u}, \mathbf{v} \rangle$。这种性质称为**半双线性**（sesquilinear）。注意在某些教材（如物理文献）中采用第二变元线性的约定，这里我们遵循数学中更常见的第一变元线性约定。

!!! definition "定义 8.3 (内积空间)"
    **内积空间**（inner product space）是配备了内积的向量空间 $(V, \langle \cdot, \cdot \rangle)$。当基域为 $\mathbb{R}$ 时称为实内积空间，当基域为 $\mathbb{C}$ 时称为复内积空间。有限维的实内积空间称为**欧几里得空间**（Euclidean space），有限维的复内积空间称为**酉空间**（unitary space）。

!!! example "例 8.1"
    **$\mathbb{R}^n$ 上的标准内积。** 对 $\mathbf{x} = (x_1, \ldots, x_n)^T, \mathbf{y} = (y_1, \ldots, y_n)^T \in \mathbb{R}^n$，定义

    $$\langle \mathbf{x}, \mathbf{y} \rangle = \mathbf{x}^T\mathbf{y} = \sum_{i=1}^n x_i y_i$$

    这是最常见的内积，即标准点积。

!!! example "例 8.2"
    **$\mathbb{C}^n$ 上的标准内积。** 对 $\mathbf{x} = (x_1, \ldots, x_n)^T, \mathbf{y} = (y_1, \ldots, y_n)^T \in \mathbb{C}^n$，定义

    $$\langle \mathbf{x}, \mathbf{y} \rangle = \mathbf{x}^H\mathbf{y} = \sum_{i=1}^n x_i \overline{y_i}$$

    其中 $\mathbf{x}^H = \overline{\mathbf{x}}^T$ 是共轭转置。

!!! example "例 8.3"
    **连续函数空间上的内积。** 在实向量空间 $C[a, b]$（区间 $[a, b]$ 上的连续实值函数）上，定义

    $$\langle f, g \rangle = \int_a^b f(t)g(t) \, dt$$

    可以验证这满足实内积的四条公理。特别地，正定性由以下事实保证：若 $f$ 连续且 $\int_a^b f(t)^2 \, dt = 0$，则 $f \equiv 0$。

!!! example "例 8.4"
    **加权内积。** 设 $w_1, w_2, \ldots, w_n > 0$ 为正实数（称为权重），对 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$ 定义

    $$\langle \mathbf{x}, \mathbf{y} \rangle_w = \sum_{i=1}^n w_i x_i y_i$$

    这是 $\mathbb{R}^n$ 上的一个内积，称为加权内积（weighted inner product）。

!!! proposition "命题 8.1"
    设 $(V, \langle \cdot, \cdot \rangle)$ 是内积空间。对任意 $\mathbf{u}, \mathbf{v} \in V$，有

    $$\langle \mathbf{u}, \mathbf{0} \rangle = \langle \mathbf{0}, \mathbf{v} \rangle = 0$$

??? proof "证明"
    由线性性，$\langle \mathbf{u}, \mathbf{0} \rangle = \langle \mathbf{u}, 0 \cdot \mathbf{v} \rangle = \bar{0} \langle \mathbf{u}, \mathbf{v} \rangle = 0$。类似地，$\langle \mathbf{0}, \mathbf{v} \rangle = \langle 0 \cdot \mathbf{u}, \mathbf{v} \rangle = 0 \cdot \langle \mathbf{u}, \mathbf{v} \rangle = 0$。$\blacksquare$

---

## 8.2 范数与距离

<div class="context-flow" markdown>

内积 → 范数 $\|\mathbf{v}\| = \sqrt{\langle \mathbf{v},\mathbf{v}\rangle}$ → **Cauchy-Schwarz**（内积空间最基本不等式） → 三角不等式 → 距离

关键洞察：平行四边形恒等式刻画了"哪些范数来自内积"

</div>

内积自然地诱导出向量的"长度"概念（范数）和向量之间的"距离"概念。

!!! definition "定义 8.4 (诱导范数)"
    设 $(V, \langle \cdot, \cdot \rangle)$ 是内积空间。由内积**诱导的范数**（induced norm）定义为

    $$\|\mathbf{v}\| = \sqrt{\langle \mathbf{v}, \mathbf{v} \rangle}, \quad \forall \mathbf{v} \in V$$

    向量 $\mathbf{v}$ 与 $\mathbf{w}$ 之间的**距离**（distance）定义为 $d(\mathbf{v}, \mathbf{w}) = \|\mathbf{v} - \mathbf{w}\|$。

!!! definition "定义 8.5 (单位向量)"
    若 $\|\mathbf{v}\| = 1$，则称 $\mathbf{v}$ 为**单位向量**（unit vector）。对任意非零向量 $\mathbf{v}$，向量 $\hat{\mathbf{v}} = \frac{\mathbf{v}}{\|\mathbf{v}\|}$ 是单位向量，将 $\mathbf{v}$ 变换为 $\hat{\mathbf{v}}$ 的过程称为**单位化**（normalization）。

### Cauchy-Schwarz 不等式

<div class="context-flow" markdown>

**洞察**：Cauchy-Schwarz 的证明核心是对 $\|\mathbf{u} - t\mathbf{v}\|^2 \ge 0$ 选取最优 $t$ ——本质是**投影**的最优性

</div>

!!! theorem "定理 8.1 (Cauchy-Schwarz 不等式)"
    设 $(V, \langle \cdot, \cdot \rangle)$ 是内积空间。对任意 $\mathbf{u}, \mathbf{v} \in V$，有

    $$|\langle \mathbf{u}, \mathbf{v} \rangle| \leq \|\mathbf{u}\| \cdot \|\mathbf{v}\|$$

    等号成立当且仅当 $\mathbf{u}$ 与 $\mathbf{v}$ 线性相关。

??? proof "证明"
    若 $\mathbf{v} = \mathbf{0}$，则两边均为 $0$，不等式显然成立。

    设 $\mathbf{v} \neq \mathbf{0}$。对任意标量 $t$，由正定性有

    $$0 \leq \|\mathbf{u} - t\mathbf{v}\|^2 = \langle \mathbf{u} - t\mathbf{v}, \mathbf{u} - t\mathbf{v} \rangle = \|\mathbf{u}\|^2 - t\langle \mathbf{v}, \mathbf{u}\rangle - \bar{t}\langle \mathbf{u}, \mathbf{v}\rangle + |t|^2\|\mathbf{v}\|^2$$

    取 $t = \dfrac{\langle \mathbf{u}, \mathbf{v}\rangle}{\|\mathbf{v}\|^2}$，代入得

    $$0 \leq \|\mathbf{u}\|^2 - \frac{|\langle \mathbf{u}, \mathbf{v}\rangle|^2}{\|\mathbf{v}\|^2}$$

    整理即得 $|\langle \mathbf{u}, \mathbf{v}\rangle|^2 \leq \|\mathbf{u}\|^2 \|\mathbf{v}\|^2$，两边开方即得所求不等式。

    等号成立当且仅当 $\|\mathbf{u} - t\mathbf{v}\|^2 = 0$，即 $\mathbf{u} = t\mathbf{v}$，此时 $\mathbf{u}, \mathbf{v}$ 线性相关。$\blacksquare$

!!! theorem "定理 8.2 (三角不等式)"
    设 $(V, \langle \cdot, \cdot \rangle)$ 是内积空间。对任意 $\mathbf{u}, \mathbf{v} \in V$，有

    $$\|\mathbf{u} + \mathbf{v}\| \leq \|\mathbf{u}\| + \|\mathbf{v}\|$$

??? proof "证明"

    $$\|\mathbf{u} + \mathbf{v}\|^2 = \langle \mathbf{u}+\mathbf{v}, \mathbf{u}+\mathbf{v}\rangle = \|\mathbf{u}\|^2 + \langle \mathbf{u}, \mathbf{v}\rangle + \langle \mathbf{v}, \mathbf{u}\rangle + \|\mathbf{v}\|^2$$

    $$= \|\mathbf{u}\|^2 + 2\operatorname{Re}\langle \mathbf{u}, \mathbf{v}\rangle + \|\mathbf{v}\|^2 \leq \|\mathbf{u}\|^2 + 2|\langle \mathbf{u}, \mathbf{v}\rangle| + \|\mathbf{v}\|^2$$

    由 Cauchy-Schwarz 不等式，$|\langle \mathbf{u}, \mathbf{v}\rangle| \leq \|\mathbf{u}\|\|\mathbf{v}\|$，因此

    $$\|\mathbf{u}+\mathbf{v}\|^2 \leq \|\mathbf{u}\|^2 + 2\|\mathbf{u}\|\|\mathbf{v}\| + \|\mathbf{v}\|^2 = (\|\mathbf{u}\| + \|\mathbf{v}\|)^2$$

    两边开方即得。$\blacksquare$

!!! theorem "定理 8.3 (平行四边形恒等式)"
    设 $(V, \langle \cdot, \cdot \rangle)$ 是内积空间。对任意 $\mathbf{u}, \mathbf{v} \in V$，有

    $$\|\mathbf{u} + \mathbf{v}\|^2 + \|\mathbf{u} - \mathbf{v}\|^2 = 2(\|\mathbf{u}\|^2 + \|\mathbf{v}\|^2)$$

??? proof "证明"
    直接计算：

    $$\|\mathbf{u}+\mathbf{v}\|^2 = \|\mathbf{u}\|^2 + 2\operatorname{Re}\langle \mathbf{u}, \mathbf{v}\rangle + \|\mathbf{v}\|^2$$

    $$\|\mathbf{u}-\mathbf{v}\|^2 = \|\mathbf{u}\|^2 - 2\operatorname{Re}\langle \mathbf{u}, \mathbf{v}\rangle + \|\mathbf{v}\|^2$$

    两式相加，中间项消去，即得 $\|\mathbf{u}+\mathbf{v}\|^2 + \|\mathbf{u}-\mathbf{v}\|^2 = 2\|\mathbf{u}\|^2 + 2\|\mathbf{v}\|^2$。$\blacksquare$

!!! note "注"
    平行四边形恒等式的几何意义是：平行四边形两条对角线的平方和等于四条边的平方和。反过来，一个范数满足平行四边形恒等式是它能够由某个内积诱导的必要充分条件（Jordan-von Neumann 定理）。

!!! theorem "定理 8.4 (极化恒等式)"
    在实内积空间中，内积可由范数恢复：

    $$\langle \mathbf{u}, \mathbf{v} \rangle = \frac{1}{4}\left(\|\mathbf{u}+\mathbf{v}\|^2 - \|\mathbf{u}-\mathbf{v}\|^2\right)$$

    在复内积空间中：

    $$\langle \mathbf{u}, \mathbf{v} \rangle = \frac{1}{4}\sum_{k=0}^{3} i^k \|\mathbf{u} + i^k\mathbf{v}\|^2 = \frac{1}{4}\left(\|\mathbf{u}+\mathbf{v}\|^2 - \|\mathbf{u}-\mathbf{v}\|^2 + i\|\mathbf{u}+i\mathbf{v}\|^2 - i\|\mathbf{u}-i\mathbf{v}\|^2\right)$$

??? proof "证明"
    对实情形：展开 $\|\mathbf{u}+\mathbf{v}\|^2 - \|\mathbf{u}-\mathbf{v}\|^2$：

    $$\|\mathbf{u}+\mathbf{v}\|^2 - \|\mathbf{u}-\mathbf{v}\|^2 = (\|\mathbf{u}\|^2 + 2\langle \mathbf{u},\mathbf{v}\rangle + \|\mathbf{v}\|^2) - (\|\mathbf{u}\|^2 - 2\langle \mathbf{u},\mathbf{v}\rangle + \|\mathbf{v}\|^2) = 4\langle \mathbf{u},\mathbf{v}\rangle$$

    复情形类似，将四个展开式加权求和后，实部与虚部恰好组合成 $4\langle \mathbf{u},\mathbf{v}\rangle$。$\blacksquare$

!!! example "例 8.5"
    在 $\mathbb{R}^3$ 中取标准内积，设 $\mathbf{u} = (1, 2, 3)^T$，$\mathbf{v} = (4, -1, 2)^T$。则

    $$\langle \mathbf{u}, \mathbf{v} \rangle = 1 \cdot 4 + 2 \cdot (-1) + 3 \cdot 2 = 8$$

    $$\|\mathbf{u}\| = \sqrt{1+4+9} = \sqrt{14}, \quad \|\mathbf{v}\| = \sqrt{16+1+4} = \sqrt{21}$$

    验证 Cauchy-Schwarz 不等式：$|8| = 8 \leq \sqrt{14} \cdot \sqrt{21} = \sqrt{294} \approx 17.15$，成立。

    $\mathbf{u}$ 与 $\mathbf{v}$ 的夹角 $\theta$ 满足 $\cos\theta = \dfrac{\langle \mathbf{u}, \mathbf{v}\rangle}{\|\mathbf{u}\|\|\mathbf{v}\|} = \dfrac{8}{\sqrt{294}}$。

---

## 8.3 正交与正交补

<div class="context-flow" markdown>

$\langle \mathbf{u},\mathbf{v}\rangle = 0$ 定义正交 → 正交集自动线性无关 → **正交补** $W^\perp$ 使 $V = W \oplus W^\perp$（直和分解）

</div>

正交性是内积空间中最重要的概念之一，它推广了"垂直"的几何概念。

!!! definition "定义 8.6 (正交)"
    设 $(V, \langle \cdot, \cdot \rangle)$ 是内积空间。

    - 两个向量 $\mathbf{u}, \mathbf{v} \in V$ 称为**正交的**（orthogonal），若 $\langle \mathbf{u}, \mathbf{v} \rangle = 0$，记作 $\mathbf{u} \perp \mathbf{v}$。
    - 向量 $\mathbf{v}$ 与集合 $S \subseteq V$ 正交，若 $\mathbf{v}$ 与 $S$ 中每个向量正交。
    - 集合 $S$ 称为**正交集**（orthogonal set），若 $S$ 中任意两个不同向量都正交。
    - 集合 $S$ 称为**标准正交集**（orthonormal set），若 $S$ 是正交集且每个向量都是单位向量。

!!! definition "定义 8.7 (正交补)"
    设 $W$ 是内积空间 $V$ 的子空间。$W$ 的**正交补**（orthogonal complement）定义为

    $$W^\perp = \{\mathbf{v} \in V : \langle \mathbf{v}, \mathbf{w} \rangle = 0, \, \forall \mathbf{w} \in W\}$$

!!! theorem "定理 8.5 (正交补的性质)"
    设 $W$ 是有限维内积空间 $V$ 的子空间。则

    1. $W^\perp$ 是 $V$ 的子空间；
    2. $W \cap W^\perp = \{\mathbf{0}\}$；
    3. $(W^\perp)^\perp = W$；
    4. $\dim W + \dim W^\perp = \dim V$。

??? proof "证明"
    (1) 若 $\mathbf{u}, \mathbf{v} \in W^\perp$，$\alpha, \beta$ 为标量，则对任意 $\mathbf{w} \in W$：

    $$\langle \alpha\mathbf{u} + \beta\mathbf{v}, \mathbf{w}\rangle = \alpha\langle \mathbf{u}, \mathbf{w}\rangle + \beta\langle \mathbf{v}, \mathbf{w}\rangle = 0$$

    故 $\alpha\mathbf{u}+\beta\mathbf{v} \in W^\perp$，因此 $W^\perp$ 是子空间。

    (2) 若 $\mathbf{v} \in W \cap W^\perp$，则 $\mathbf{v} \in W$ 且 $\langle \mathbf{v}, \mathbf{w}\rangle = 0$ 对所有 $\mathbf{w} \in W$ 成立。特别地取 $\mathbf{w} = \mathbf{v}$，得 $\langle \mathbf{v}, \mathbf{v}\rangle = 0$，由正定性知 $\mathbf{v} = \mathbf{0}$。

    (3) 和 (4) 将在正交投影定理之后证明。$\blacksquare$

!!! theorem "定理 8.6 (正交直和分解)"
    设 $V$ 是有限维内积空间，$W$ 是 $V$ 的子空间。则

    $$V = W \oplus W^\perp$$

    即任意 $\mathbf{v} \in V$ 可以唯一地写成 $\mathbf{v} = \mathbf{w} + \mathbf{w}'$，其中 $\mathbf{w} \in W$，$\mathbf{w}' \in W^\perp$。

??? proof "证明"
    设 $\{\mathbf{e}_1, \ldots, \mathbf{e}_k\}$ 是 $W$ 的一个标准正交基。对任意 $\mathbf{v} \in V$，令

    $$\mathbf{w} = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i \rangle \mathbf{e}_i, \quad \mathbf{w}' = \mathbf{v} - \mathbf{w}$$

    显然 $\mathbf{w} \in W$。验证 $\mathbf{w}' \in W^\perp$：对任意 $j = 1, \ldots, k$，

    $$\langle \mathbf{w}', \mathbf{e}_j \rangle = \langle \mathbf{v} - \mathbf{w}, \mathbf{e}_j\rangle = \langle \mathbf{v}, \mathbf{e}_j\rangle - \sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i\rangle\langle \mathbf{e}_i, \mathbf{e}_j\rangle = \langle \mathbf{v}, \mathbf{e}_j\rangle - \langle \mathbf{v}, \mathbf{e}_j\rangle = 0$$

    因此 $\mathbf{w}' \perp \mathbf{e}_j$ 对所有 $j$ 成立，从而 $\mathbf{w}' \in W^\perp$。故 $\mathbf{v} = \mathbf{w} + \mathbf{w}'$ 是所求的分解。

    唯一性由 $W \cap W^\perp = \{\mathbf{0}\}$ 保证：若 $\mathbf{v} = \mathbf{w}_1 + \mathbf{w}_1' = \mathbf{w}_2 + \mathbf{w}_2'$，则 $\mathbf{w}_1 - \mathbf{w}_2 = \mathbf{w}_2' - \mathbf{w}_1' \in W \cap W^\perp = \{\mathbf{0}\}$。$\blacksquare$

!!! proposition "命题 8.2 (正交集的线性无关性)"
    设 $S = \{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ 是内积空间 $V$ 中由非零向量组成的正交集。则 $S$ 是线性无关的。

??? proof "证明"
    设 $\alpha_1\mathbf{v}_1 + \cdots + \alpha_k\mathbf{v}_k = \mathbf{0}$。对任意 $j$，用 $\mathbf{v}_j$ 取内积：

    $$0 = \left\langle \sum_{i=1}^k \alpha_i\mathbf{v}_i, \mathbf{v}_j \right\rangle = \sum_{i=1}^k \alpha_i \langle \mathbf{v}_i, \mathbf{v}_j\rangle = \alpha_j \|\mathbf{v}_j\|^2$$

    因 $\mathbf{v}_j \neq \mathbf{0}$，故 $\|\mathbf{v}_j\|^2 > 0$，从而 $\alpha_j = 0$。$\blacksquare$

!!! example "例 8.6"
    在 $\mathbb{R}^3$ 中，设 $W = \operatorname{span}\{(1, 1, 0)^T, (0, 1, 1)^T\}$。求 $W^\perp$。

    设 $\mathbf{v} = (x, y, z)^T \in W^\perp$，则

    $$\langle \mathbf{v}, (1,1,0)^T\rangle = x + y = 0, \quad \langle \mathbf{v}, (0,1,1)^T\rangle = y + z = 0$$

    解得 $y = -x$，$z = -y = x$。故 $W^\perp = \operatorname{span}\{(1, -1, 1)^T\}$。

    验证：$\dim W + \dim W^\perp = 2 + 1 = 3 = \dim \mathbb{R}^3$。

---

## 8.4 正交投影定理

<div class="context-flow" markdown>

正交直和 $V = W \oplus W^\perp$ → 投影 $\operatorname{proj}_W$ 是**最佳逼近**（距离最小） → 投影算子兼具幂等性 $P^2 = P$ 与自伴性 $P^* = P$

</div>

!!! definition "定义 8.8 (正交投影)"
    设 $V$ 是有限维内积空间，$W$ 是 $V$ 的子空间。根据正交直和分解 $V = W \oplus W^\perp$，对任意 $\mathbf{v} \in V$，定义**正交投影**（orthogonal projection）$\operatorname{proj}_W: V \to V$ 为

    $$\operatorname{proj}_W(\mathbf{v}) = \mathbf{w}$$

    其中 $\mathbf{v} = \mathbf{w} + \mathbf{w}'$，$\mathbf{w} \in W$，$\mathbf{w}' \in W^\perp$。

    若 $\{\mathbf{e}_1, \ldots, \mathbf{e}_k\}$ 是 $W$ 的标准正交基，则

    $$\operatorname{proj}_W(\mathbf{v}) = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i \rangle \mathbf{e}_i$$

<div class="context-flow" markdown>

**洞察**：最佳逼近 = 正交投影，证明核心是勾股定理 $\|\mathbf{v}-\mathbf{w}\|^2 = \|\mathbf{v}-\hat{\mathbf{v}}\|^2 + \|\hat{\mathbf{v}}-\mathbf{w}\|^2$ → 直通 Ch11 最小二乘

</div>

!!! theorem "定理 8.7 (最佳逼近定理)"
    设 $V$ 是内积空间，$W$ 是 $V$ 的有限维子空间。对任意 $\mathbf{v} \in V$，$\operatorname{proj}_W(\mathbf{v})$ 是 $W$ 中距 $\mathbf{v}$ 最近的向量。即对任意 $\mathbf{w} \in W$，

    $$\|\mathbf{v} - \operatorname{proj}_W(\mathbf{v})\| \leq \|\mathbf{v} - \mathbf{w}\|$$

    等号成立当且仅当 $\mathbf{w} = \operatorname{proj}_W(\mathbf{v})$。

??? proof "证明"
    记 $\hat{\mathbf{v}} = \operatorname{proj}_W(\mathbf{v})$。对任意 $\mathbf{w} \in W$，

    $$\|\mathbf{v} - \mathbf{w}\|^2 = \|(\mathbf{v} - \hat{\mathbf{v}}) + (\hat{\mathbf{v}} - \mathbf{w})\|^2$$

    注意 $\mathbf{v} - \hat{\mathbf{v}} \in W^\perp$ 且 $\hat{\mathbf{v}} - \mathbf{w} \in W$，因此两者正交。由勾股定理：

    $$\|\mathbf{v} - \mathbf{w}\|^2 = \|\mathbf{v} - \hat{\mathbf{v}}\|^2 + \|\hat{\mathbf{v}} - \mathbf{w}\|^2 \geq \|\mathbf{v} - \hat{\mathbf{v}}\|^2$$

    等号成立当且仅当 $\|\hat{\mathbf{v}} - \mathbf{w}\|^2 = 0$，即 $\mathbf{w} = \hat{\mathbf{v}}$。$\blacksquare$

!!! proposition "命题 8.3 (正交投影算子的性质)"
    设 $P = \operatorname{proj}_W$ 是正交投影算子。则

    1. $P$ 是线性算子；
    2. $P^2 = P$（幂等性）；
    3. $\operatorname{Im}(P) = W$，$\operatorname{Ker}(P) = W^\perp$；
    4. $P^* = P$（自伴性）；
    5. 对任意 $\mathbf{v} \in V$，$\|\mathbf{v}\|^2 = \|P\mathbf{v}\|^2 + \|\mathbf{v} - P\mathbf{v}\|^2$。

??? proof "证明"
    (1) 设 $\mathbf{v}_1 = \mathbf{w}_1 + \mathbf{w}_1'$，$\mathbf{v}_2 = \mathbf{w}_2 + \mathbf{w}_2'$，则 $\alpha\mathbf{v}_1 + \beta\mathbf{v}_2 = (\alpha\mathbf{w}_1 + \beta\mathbf{w}_2) + (\alpha\mathbf{w}_1' + \beta\mathbf{w}_2')$，故 $P(\alpha\mathbf{v}_1 + \beta\mathbf{v}_2) = \alpha\mathbf{w}_1 + \beta\mathbf{w}_2 = \alpha P\mathbf{v}_1 + \beta P\mathbf{v}_2$。

    (2) 对 $\mathbf{v} = \mathbf{w} + \mathbf{w}'$，$P\mathbf{v} = \mathbf{w} \in W$，故 $P(P\mathbf{v}) = P\mathbf{w} = \mathbf{w} = P\mathbf{v}$。

    (3) 显然。

    (4) 设 $\mathbf{u} = \mathbf{w}_1 + \mathbf{w}_1'$，$\mathbf{v} = \mathbf{w}_2 + \mathbf{w}_2'$。则 $\langle P\mathbf{u}, \mathbf{v}\rangle = \langle \mathbf{w}_1, \mathbf{w}_2 + \mathbf{w}_2'\rangle = \langle \mathbf{w}_1, \mathbf{w}_2\rangle$ 且 $\langle \mathbf{u}, P\mathbf{v}\rangle = \langle \mathbf{w}_1 + \mathbf{w}_1', \mathbf{w}_2\rangle = \langle \mathbf{w}_1, \mathbf{w}_2\rangle$。

    (5) 由 $\mathbf{v} = P\mathbf{v} + (\mathbf{v} - P\mathbf{v})$ 且两者正交，由勾股定理即得。$\blacksquare$

!!! example "例 8.7"
    设 $W = \operatorname{span}\left\{\mathbf{e}_1 = \frac{1}{\sqrt{2}}(1, 1, 0)^T, \, \mathbf{e}_2 = \frac{1}{\sqrt{6}}(1, -1, 2)^T\right\}$ 是 $\mathbb{R}^3$ 中的子空间（$\mathbf{e}_1, \mathbf{e}_2$ 已是标准正交基）。求 $\mathbf{v} = (1, 2, 3)^T$ 在 $W$ 上的正交投影。

    $$\operatorname{proj}_W(\mathbf{v}) = \langle \mathbf{v}, \mathbf{e}_1\rangle \mathbf{e}_1 + \langle \mathbf{v}, \mathbf{e}_2\rangle \mathbf{e}_2$$

    $$\langle \mathbf{v}, \mathbf{e}_1\rangle = \frac{1}{\sqrt{2}}(1+2) = \frac{3}{\sqrt{2}}$$

    $$\langle \mathbf{v}, \mathbf{e}_2\rangle = \frac{1}{\sqrt{6}}(1-2+6) = \frac{5}{\sqrt{6}}$$

    $$\operatorname{proj}_W(\mathbf{v}) = \frac{3}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\\0\end{pmatrix} + \frac{5}{\sqrt{6}} \cdot \frac{1}{\sqrt{6}}\begin{pmatrix}1\\-1\\2\end{pmatrix} = \frac{3}{2}\begin{pmatrix}1\\1\\0\end{pmatrix} + \frac{5}{6}\begin{pmatrix}1\\-1\\2\end{pmatrix} = \begin{pmatrix}7/3\\2/3\\5/3\end{pmatrix}$$

---

## 8.5 Bessel 不等式与 Parseval 等式

<div class="context-flow" markdown>

投影只捕获"部分能量" → **Bessel**：$\sum |\langle \mathbf{v},\mathbf{e}_i\rangle|^2 \le \|\mathbf{v}\|^2$；补全为基时取等 → **Parseval** → 通向 Fourier 分析

</div>

!!! theorem "定理 8.8 (Bessel 不等式)"
    设 $(V, \langle \cdot, \cdot \rangle)$ 是内积空间，$\{\mathbf{e}_1, \ldots, \mathbf{e}_k\}$ 是 $V$ 中的标准正交集。则对任意 $\mathbf{v} \in V$，

    $$\sum_{i=1}^k |\langle \mathbf{v}, \mathbf{e}_i \rangle|^2 \leq \|\mathbf{v}\|^2$$

??? proof "证明"
    令 $W = \operatorname{span}\{\mathbf{e}_1, \ldots, \mathbf{e}_k\}$，$\hat{\mathbf{v}} = \operatorname{proj}_W(\mathbf{v}) = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i\rangle\mathbf{e}_i$。

    由勾股定理：

    $$\|\mathbf{v}\|^2 = \|\hat{\mathbf{v}}\|^2 + \|\mathbf{v} - \hat{\mathbf{v}}\|^2 \geq \|\hat{\mathbf{v}}\|^2$$

    而 $\|\hat{\mathbf{v}}\|^2 = \left\|\sum_{i=1}^k \langle \mathbf{v}, \mathbf{e}_i\rangle \mathbf{e}_i\right\|^2 = \sum_{i=1}^k |\langle \mathbf{v}, \mathbf{e}_i\rangle|^2$（利用标准正交性）。

    故 $\sum_{i=1}^k |\langle \mathbf{v}, \mathbf{e}_i\rangle|^2 \leq \|\mathbf{v}\|^2$。$\blacksquare$

!!! theorem "定理 8.9 (Parseval 等式)"
    设 $V$ 是有限维内积空间，$\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 是 $V$ 的一个标准正交基。则对任意 $\mathbf{v} \in V$，

    $$\sum_{i=1}^n |\langle \mathbf{v}, \mathbf{e}_i \rangle|^2 = \|\mathbf{v}\|^2$$

    更一般地，对任意 $\mathbf{u}, \mathbf{v} \in V$，

    $$\langle \mathbf{u}, \mathbf{v} \rangle = \sum_{i=1}^n \langle \mathbf{u}, \mathbf{e}_i \rangle \overline{\langle \mathbf{v}, \mathbf{e}_i \rangle}$$

??? proof "证明"
    因 $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 是 $V$ 的标准正交基，$\mathbf{v} = \sum_{i=1}^n \langle \mathbf{v}, \mathbf{e}_i\rangle \mathbf{e}_i$，从而

    $$\|\mathbf{v}\|^2 = \left\langle \sum_{i=1}^n \langle \mathbf{v}, \mathbf{e}_i\rangle \mathbf{e}_i, \sum_{j=1}^n \langle \mathbf{v}, \mathbf{e}_j\rangle \mathbf{e}_j \right\rangle = \sum_{i=1}^n \sum_{j=1}^n \langle \mathbf{v}, \mathbf{e}_i\rangle \overline{\langle \mathbf{v}, \mathbf{e}_j\rangle}\langle \mathbf{e}_i, \mathbf{e}_j\rangle = \sum_{i=1}^n |\langle \mathbf{v}, \mathbf{e}_i\rangle|^2$$

    一般的恒等式（Parseval 恒等式的一般形式）类似证明。$\blacksquare$

!!! note "注"
    Parseval 等式可以理解为：向量在标准正交基下的坐标分量的模的平方和等于向量的范数的平方。Bessel 不等式则说明，如果标准正交集不是基（即不张成整个空间），则"能量"（模的平方）不能被完全捕获。

!!! example "例 8.8"
    在 $\mathbb{R}^3$ 中取标准正交基 $\{\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3\}$，设 $\mathbf{v} = (3, 4, 0)^T$。

    Parseval 等式：$|\langle \mathbf{v}, \mathbf{e}_1\rangle|^2 + |\langle \mathbf{v}, \mathbf{e}_2\rangle|^2 + |\langle \mathbf{v}, \mathbf{e}_3\rangle|^2 = 9 + 16 + 0 = 25 = \|\mathbf{v}\|^2$。

    仅用前两个基向量（Bessel 不等式）：$|\langle \mathbf{v}, \mathbf{e}_1\rangle|^2 + |\langle \mathbf{v}, \mathbf{e}_2\rangle|^2 = 25 \leq 25 = \|\mathbf{v}\|^2$。此例中恰好取等，因为 $\mathbf{v}$ 本身在前两个分量张成的子空间中。

---

## 8.6 伴随算子

<div class="context-flow" markdown>

将内积"转移"到另一侧：$\langle T\mathbf{v},\mathbf{w}\rangle = \langle \mathbf{v},T^*\mathbf{w}\rangle$ → 标准正交基下 $T^*$ 的矩阵 = $A^H$ → 连接到 Ch10 的分解理论

</div>

!!! definition "定义 8.9 (伴随算子)"
    设 $V$ 是有限维内积空间，$T: V \to V$ 是线性算子。$T$ 的**伴随算子**（adjoint operator）$T^*: V \to V$ 是满足以下条件的唯一线性算子：

    $$\langle T\mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{v}, T^*\mathbf{w} \rangle, \quad \forall \mathbf{v}, \mathbf{w} \in V$$

!!! theorem "定理 8.10 (伴随算子的存在唯一性)"
    设 $V$ 是有限维内积空间，$T: V \to V$ 是线性算子。则 $T^*$ 存在且唯一。

??? proof "证明"
    **存在性：** 设 $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 是 $V$ 的标准正交基。$T$ 在此基下的矩阵为 $A = (a_{ij})$，其中 $a_{ij} = \langle T\mathbf{e}_j, \mathbf{e}_i\rangle$。

    定义 $T^*$ 为以矩阵 $B = (b_{ij})$ 表示的算子，其中 $b_{ij} = \overline{a_{ji}}$，即 $B = A^H$（共轭转置）。则

    $$\langle T\mathbf{v}, \mathbf{w}\rangle = (A\mathbf{v})^H \mathbf{w} = \mathbf{v}^H A^H \mathbf{w} = \mathbf{v}^H (B\mathbf{w}) = \langle \mathbf{v}, T^*\mathbf{w}\rangle$$

    **唯一性：** 若 $S$ 也满足 $\langle T\mathbf{v}, \mathbf{w}\rangle = \langle \mathbf{v}, S\mathbf{w}\rangle$，则 $\langle \mathbf{v}, T^*\mathbf{w} - S\mathbf{w}\rangle = 0$ 对所有 $\mathbf{v}$ 成立。取 $\mathbf{v} = T^*\mathbf{w} - S\mathbf{w}$，得 $T^*\mathbf{w} = S\mathbf{w}$ 对所有 $\mathbf{w}$ 成立。$\blacksquare$

!!! theorem "定理 8.11 (伴随算子的性质)"
    设 $V$ 是有限维内积空间，$S, T: V \to V$ 是线性算子，$\alpha$ 是标量。则

    1. $(S + T)^* = S^* + T^*$
    2. $(\alpha T)^* = \bar{\alpha} T^*$
    3. $(ST)^* = T^* S^*$
    4. $(T^*)^* = T$
    5. $\operatorname{Ker}(T^*) = (\operatorname{Im}(T))^\perp$
    6. $\operatorname{Im}(T^*) = (\operatorname{Ker}(T))^\perp$

??? proof "证明"
    我们证明性质 5，其他性质类似。

    $\mathbf{w} \in \operatorname{Ker}(T^*)$ $\Leftrightarrow$ $T^*\mathbf{w} = \mathbf{0}$ $\Leftrightarrow$ $\langle \mathbf{v}, T^*\mathbf{w}\rangle = 0, \, \forall \mathbf{v} \in V$ $\Leftrightarrow$ $\langle T\mathbf{v}, \mathbf{w}\rangle = 0, \, \forall \mathbf{v} \in V$ $\Leftrightarrow$ $\mathbf{w} \in (\operatorname{Im}(T))^\perp$。$\blacksquare$

!!! note "注"
    在标准正交基下，$T$ 的矩阵为 $A$ 时，$T^*$ 的矩阵为 $A^H = \overline{A}^T$（共轭转置）。对实内积空间，$T^*$ 的矩阵就是 $A^T$。

!!! example "例 8.9"
    设 $T: \mathbb{C}^2 \to \mathbb{C}^2$ 在标准基下的矩阵为

    $$A = \begin{pmatrix} 1+i & 2 \\ 3i & 4-i \end{pmatrix}$$

    则 $T^*$ 的矩阵为

    $$A^H = \overline{A}^T = \begin{pmatrix} 1-i & -3i \\ 2 & 4+i \end{pmatrix}$$

!!! example "例 8.10"
    设 $V = C[0,1]$ 带内积 $\langle f, g\rangle = \int_0^1 f(t)\overline{g(t)} \, dt$。定义线性算子 $T: V \to V$，$(Tf)(t) = tf(t)$（乘以自变量）。则 $T^* = T$，因为

    $$\langle Tf, g\rangle = \int_0^1 tf(t)\overline{g(t)} \, dt = \int_0^1 f(t)\overline{tg(t)} \, dt = \langle f, Tg\rangle$$

    其中最后一步利用了 $t$ 是实数。

---

## 8.7 正规算子

<div class="context-flow" markdown>

$TT^* = T^*T$（与伴随可交换）→ 不同特征值的特征向量自动正交 → 正规 = 可**酉对角化**的充要条件（复情形）

</div>

!!! definition "定义 8.10 (正规算子)"
    设 $V$ 是有限维内积空间，$T: V \to V$ 是线性算子。若

    $$TT^* = T^*T$$

    则称 $T$ 是**正规算子**（normal operator）。相应地，若 $n \times n$ 矩阵 $A$ 满足 $AA^H = A^HA$，则称 $A$ 为**正规矩阵**（normal matrix）。

!!! note "注"
    自伴算子（$T^* = T$）、酉算子（$T^*T = I$）、反自伴算子（$T^* = -T$）都是正规算子的特例。

!!! theorem "定理 8.12 (正规算子的等价条件)"
    设 $V$ 是有限维内积空间，$T: V \to V$ 是线性算子。则以下条件等价：

    1. $T$ 是正规的；
    2. $\|T\mathbf{v}\| = \|T^*\mathbf{v}\|$ 对所有 $\mathbf{v} \in V$ 成立；
    3. $T - \lambda I$ 对任意标量 $\lambda$ 都是正规的；
    4. 若 $T\mathbf{v} = \lambda\mathbf{v}$，则 $T^*\mathbf{v} = \bar{\lambda}\mathbf{v}$；
    5. $T$ 的属于不同特征值的特征向量正交。

??? proof "证明"
    **(1)$\Rightarrow$(2)：** $\|T\mathbf{v}\|^2 = \langle T\mathbf{v}, T\mathbf{v}\rangle = \langle T^*T\mathbf{v}, \mathbf{v}\rangle = \langle TT^*\mathbf{v}, \mathbf{v}\rangle = \langle T^*\mathbf{v}, T^*\mathbf{v}\rangle = \|T^*\mathbf{v}\|^2$。

    **(2)$\Rightarrow$(1)：** 条件 (2) 意味着 $\langle (T^*T - TT^*)\mathbf{v}, \mathbf{v}\rangle = 0$ 对所有 $\mathbf{v}$ 成立。由于 $T^*T - TT^*$ 是自伴的，这推出 $T^*T - TT^* = 0$。

    **(1)$\Rightarrow$(3)：** $(T-\lambda I)^* = T^* - \bar{\lambda}I$。直接验证：

    $$(T-\lambda I)(T-\lambda I)^* = TT^* - \lambda T^* - \bar{\lambda}T + |\lambda|^2 I$$

    $$(T-\lambda I)^*(T-\lambda I) = T^*T - \bar{\lambda}T - \lambda T^* + |\lambda|^2 I$$

    由 $TT^* = T^*T$ 知两者相等。

    **(3)$\Rightarrow$(4)：** 若 $T\mathbf{v} = \lambda\mathbf{v}$，即 $(T-\lambda I)\mathbf{v} = \mathbf{0}$。由 (3) 知 $T-\lambda I$ 正规，由 (2) 得 $\|(T-\lambda I)^*\mathbf{v}\| = \|(T-\lambda I)\mathbf{v}\| = 0$，即 $(T^*-\bar{\lambda}I)\mathbf{v} = \mathbf{0}$，故 $T^*\mathbf{v} = \bar{\lambda}\mathbf{v}$。

    **(4)$\Rightarrow$(5)：** 设 $T\mathbf{v}_1 = \lambda_1\mathbf{v}_1$，$T\mathbf{v}_2 = \lambda_2\mathbf{v}_2$，$\lambda_1 \neq \lambda_2$。由 (4) 知 $T^*\mathbf{v}_1 = \bar{\lambda}_1\mathbf{v}_1$。则

    $$\lambda_1\langle \mathbf{v}_1, \mathbf{v}_2\rangle = \langle T\mathbf{v}_1, \mathbf{v}_2\rangle = \langle \mathbf{v}_1, T^*\mathbf{v}_2\rangle = \langle \mathbf{v}_1, \bar{\lambda}_2\mathbf{v}_2\rangle = \lambda_2\langle \mathbf{v}_1, \mathbf{v}_2\rangle$$

    故 $(\lambda_1 - \lambda_2)\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$。因 $\lambda_1 \neq \lambda_2$，得 $\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$。$\blacksquare$

!!! example "例 8.11"
    矩阵 $A = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$ 是正规矩阵，因为

    $$AA^T = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$$

    $$A^TA = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$$

    但 $A$ 不是对称矩阵（$A \neq A^T$），也不是正交矩阵（$A^TA \neq I$）。

---

## 8.8 自伴算子

<div class="context-flow" markdown>

正规算子的最重要特例：$T^* = T$ → 特征值**全为实数** → 不同特征值的特征向量正交 → 直通谱定理

</div>

!!! definition "定义 8.11 (自伴算子)"
    设 $V$ 是有限维内积空间，$T: V \to V$ 是线性算子。若 $T^* = T$，即

    $$\langle T\mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{v}, T\mathbf{w} \rangle, \quad \forall \mathbf{v}, \mathbf{w} \in V$$

    则称 $T$ 是**自伴算子**（self-adjoint operator）。在实内积空间中也称为**对称算子**。相应地，满足 $A^H = A$ 的矩阵称为 **Hermitian 矩阵**（实情形即对称矩阵 $A^T = A$）。

!!! theorem "定理 8.13 (自伴算子的特征值都是实数)"
    设 $V$ 是有限维内积空间（实或复），$T: V \to V$ 是自伴算子。则 $T$ 的所有特征值都是实数。

??? proof "证明"
    设 $T\mathbf{v} = \lambda\mathbf{v}$，$\mathbf{v} \neq \mathbf{0}$。则

    $$\lambda\|\mathbf{v}\|^2 = \lambda\langle \mathbf{v}, \mathbf{v}\rangle = \langle \lambda\mathbf{v}, \mathbf{v}\rangle = \langle T\mathbf{v}, \mathbf{v}\rangle = \langle \mathbf{v}, T\mathbf{v}\rangle = \langle \mathbf{v}, \lambda\mathbf{v}\rangle = \bar{\lambda}\langle \mathbf{v}, \mathbf{v}\rangle = \bar{\lambda}\|\mathbf{v}\|^2$$

    因 $\|\mathbf{v}\|^2 > 0$，故 $\lambda = \bar{\lambda}$，即 $\lambda \in \mathbb{R}$。$\blacksquare$

!!! theorem "定理 8.14 (自伴算子属于不同特征值的特征向量正交)"
    设 $T$ 是自伴算子，$\lambda_1 \neq \lambda_2$ 是 $T$ 的两个不同特征值，$\mathbf{v}_1, \mathbf{v}_2$ 分别是对应的特征向量。则 $\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$。

??? proof "证明"
    由定理 8.12 性质 (5)（自伴算子是正规算子的特例），或直接证明：

    $$\lambda_1\langle \mathbf{v}_1, \mathbf{v}_2\rangle = \langle T\mathbf{v}_1, \mathbf{v}_2\rangle = \langle \mathbf{v}_1, T\mathbf{v}_2\rangle = \lambda_2\langle \mathbf{v}_1, \mathbf{v}_2\rangle$$

    其中最后一步用了 $\lambda_2 \in \mathbb{R}$（故 $\bar{\lambda}_2 = \lambda_2$）。因此 $(\lambda_1 - \lambda_2)\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$，由 $\lambda_1 \neq \lambda_2$ 得 $\langle \mathbf{v}_1, \mathbf{v}_2\rangle = 0$。$\blacksquare$

!!! proposition "命题 8.4"
    设 $V$ 是有限维**复**内积空间，$T: V \to V$ 是线性算子。若 $\langle T\mathbf{v}, \mathbf{v}\rangle = 0$ 对所有 $\mathbf{v} \in V$ 成立，则 $T = 0$。

??? proof "证明"
    对任意 $\mathbf{u}, \mathbf{v} \in V$，利用极化恒等式：

    $$\langle T\mathbf{u}, \mathbf{v}\rangle = \frac{1}{4}\sum_{k=0}^3 i^k \langle T(\mathbf{u}+i^k\mathbf{v}), \mathbf{u}+i^k\mathbf{v}\rangle = 0$$

    取 $\mathbf{v} = T\mathbf{u}$ 得 $\|T\mathbf{u}\|^2 = 0$，故 $T\mathbf{u} = \mathbf{0}$ 对所有 $\mathbf{u}$ 成立。$\blacksquare$

!!! note "注"
    命题 8.4 在实内积空间中不成立。例如 $T: \mathbb{R}^2 \to \mathbb{R}^2$，$T(x,y)^T = (-y, x)^T$（旋转 90 度），满足 $\langle T\mathbf{v}, \mathbf{v}\rangle = 0$ 对所有 $\mathbf{v}$ 成立，但 $T \neq 0$。

!!! example "例 8.12"
    矩阵 $A = \begin{pmatrix} 2 & 1+i \\ 1-i & 3 \end{pmatrix}$ 是 Hermitian 矩阵，因为 $A^H = A$。其特征多项式为

    $$\det(A - \lambda I) = (2-\lambda)(3-\lambda) - |1+i|^2 = \lambda^2 - 5\lambda + 4 = (\lambda - 1)(\lambda - 4)$$

    特征值为 $\lambda_1 = 1, \lambda_2 = 4$，均为实数。

---

## 8.9 谱定理

<div class="context-flow" markdown>

**全章高潮**：内积 + 自伴/正规 → 存在**标准正交特征基** → $A = Q\Lambda Q^T$（实对称）/ $A = U\Lambda U^H$（正规）

→ Ch9 用此对角化二次型，Ch10 用此做谱分解，Ch11 SVD 本质是"两侧的谱定理"

</div>

谱定理是线性代数中最重要、最深刻的定理之一。它指出正规算子（包括自伴算子）可以通过标准正交基对角化。

### 实谱定理

<div class="context-flow" markdown>

**洞察**：谱定理的归纳证明关键——找到一个特征向量后，其正交补是**不变子空间**，从而可递归降维

</div>

!!! theorem "定理 8.15 (实对称矩阵的谱定理)"
    设 $A$ 是 $n \times n$ 实对称矩阵（即 $A^T = A$）。则

    1. $A$ 的所有特征值都是实数；
    2. $A$ 的属于不同特征值的特征向量正交；
    3. 存在正交矩阵 $Q$（即 $Q^TQ = I$）使得 $A = Q\Lambda Q^T$，其中 $\Lambda$ 是对角矩阵。

    等价地，$\mathbb{R}^n$ 存在由 $A$ 的特征向量组成的标准正交基。

??? proof "证明"
    性质 1 和 2 已在定理 8.13 和 8.14 中证明。下面证明性质 3。

    对 $n$ 用数学归纳法。$n=1$ 时显然。设结论对 $n-1$ 阶矩阵成立。

    **关键步骤：** 首先需要证明实对称矩阵 $A$ 有实特征值。$A$ 的特征多项式 $\det(A - \lambda I)$ 是 $n$ 次实系数多项式。在 $\mathbb{C}$ 上它有 $n$ 个根（计重数），而由定理 8.13，所有根都是实数。

    设 $\lambda_1$ 是 $A$ 的一个（实）特征值，$\mathbf{q}_1$ 是对应的单位特征向量。将 $\mathbf{q}_1$ 扩充为 $\mathbb{R}^n$ 的标准正交基 $\{\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n\}$，令 $Q_1 = [\mathbf{q}_1 | \mathbf{q}_2 | \cdots | \mathbf{q}_n]$。则

    $$Q_1^T A Q_1 = \begin{pmatrix} \lambda_1 & \mathbf{b}^T \\ \mathbf{0} & A_1 \end{pmatrix}$$

    由 $A$ 的对称性，$Q_1^TAQ_1$ 也是对称的，故 $\mathbf{b} = \mathbf{0}$ 且 $A_1$ 是 $(n-1)\times(n-1)$ 实对称矩阵。

    由归纳假设，存在 $(n-1)$ 阶正交矩阵 $P_1$ 使 $P_1^TA_1P_1 = \operatorname{diag}(\lambda_2, \ldots, \lambda_n)$。令

    $$P = \begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & P_1 \end{pmatrix}$$

    则 $Q = Q_1P$ 是正交矩阵且 $Q^TAQ = \operatorname{diag}(\lambda_1, \lambda_2, \ldots, \lambda_n)$。$\blacksquare$

### 复谱定理

!!! theorem "定理 8.16 (Hermitian 矩阵的谱定理)"
    设 $A$ 是 $n \times n$ Hermitian 矩阵（即 $A^H = A$）。则存在酉矩阵 $U$（即 $U^HU = I$）使得

    $$A = U\Lambda U^H$$

    其中 $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，$\lambda_i$ 均为实数。

??? proof "证明"
    证明与实对称矩阵的情形完全类似，只需将"正交"替换为"酉"，"转置"替换为"共轭转置"。$\blacksquare$

!!! theorem "定理 8.17 (正规算子的谱定理——复情形)"
    设 $V$ 是有限维**复**内积空间，$T: V \to V$ 是线性算子。则 $T$ 是正规的当且仅当存在 $V$ 的标准正交基使得 $T$ 在此基下的矩阵是对角矩阵。

    等价地，$n \times n$ 复矩阵 $A$ 是正规的当且仅当存在酉矩阵 $U$ 使得

    $$A = U\Lambda U^H$$

    其中 $\Lambda$ 是对角矩阵（对角元素可以是复数）。

??? proof "证明"
    **充分性：** 若 $A = UDU^H$（$D$ 为对角矩阵），则 $A^H = UD^HU^H$，从而

    $$AA^H = UDD^HU^H = UD^HDU^H = A^HA$$

    （因为对角矩阵之间一定可交换。）

    **必要性：** 对 $n$ 用归纳法。$n = 1$ 显然。设对 $n-1$ 成立。

    由代数基本定理，$T$ 在 $\mathbb{C}$ 上有特征值 $\lambda_1$。设 $\mathbf{e}_1$ 是对应的单位特征向量。

    令 $W = \operatorname{span}\{\mathbf{e}_1\}$。由定理 8.12 性质 (4)，$T^*\mathbf{e}_1 = \bar{\lambda}_1\mathbf{e}_1$，因此 $W$ 是 $T^*$-不变的。

    我们证明 $W^\perp$ 是 $T$-不变的：若 $\mathbf{v} \in W^\perp$，则对任意 $\mathbf{w} \in W$，

    $$\langle T\mathbf{v}, \mathbf{w}\rangle = \langle \mathbf{v}, T^*\mathbf{w}\rangle$$

    因 $W$ 是 $T^*$-不变的，$T^*\mathbf{w} \in W$，故 $\langle \mathbf{v}, T^*\mathbf{w}\rangle = 0$。从而 $T\mathbf{v} \in W^\perp$。

    $T|_{W^\perp}$ 是 $W^\perp$ 上的正规算子（继承正规性）。由归纳假设，存在 $W^\perp$ 的标准正交基 $\{\mathbf{e}_2, \ldots, \mathbf{e}_n\}$ 使 $T|_{W^\perp}$ 在此基下为对角矩阵。则 $\{\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\}$ 是 $V$ 的标准正交基，$T$ 在此基下为对角矩阵。$\blacksquare$

!!! note "注"
    **实情形的谱定理**需要额外条件。一般的实正规矩阵不一定可以实正交对角化（例如旋转矩阵有复特征值）。实正规矩阵可以正交相似于如下分块对角形：对角块是 $1\times 1$ 实数或 $2\times 2$ 旋转-缩放块 $\begin{pmatrix} a & -b \\ b & a \end{pmatrix}$（$b \neq 0$），对应于共轭复特征值对 $a \pm bi$。

!!! example "例 8.13"
    对实对称矩阵 $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ 进行正交对角化。

    特征多项式：$\det(A - \lambda I) = (2-\lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda-1)(\lambda-3)$。

    $\lambda_1 = 1$：$(A-I)\mathbf{x} = \mathbf{0}$ 给出 $\mathbf{v}_1 = \frac{1}{\sqrt{2}}(1, -1)^T$。

    $\lambda_2 = 3$：$(A-3I)\mathbf{x} = \mathbf{0}$ 给出 $\mathbf{v}_2 = \frac{1}{\sqrt{2}}(1, 1)^T$。

    验证正交性：$\langle \mathbf{v}_1, \mathbf{v}_2\rangle = \frac{1}{2}(1-1) = 0$。

    令 $Q = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$，则 $A = Q\begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}Q^T$。

!!! example "例 8.14"
    设 $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$，$A$ 是实正规矩阵（$AA^T = A^TA = I$），但其特征值为 $\pm i$，不是实数。因此 $A$ 不能实正交对角化，但可以酉对角化：

    $$A = U\begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}U^H, \quad U = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ i & -i \end{pmatrix}$$

---

## 8.10 酉算子与正交算子

<div class="context-flow" markdown>

$T^*T = I$ ↔ 保内积/保范数/保标准正交基 → 特征值 $|\lambda|=1$ → 正交矩阵 $\det = \pm 1$（旋转或反射） → Ch10 QR 分解/Ch11 SVD 的"旋转部分"

</div>

!!! definition "定义 8.12 (酉算子与正交算子)"
    设 $V$ 是有限维内积空间，$T: V \to V$ 是线性算子。

    - 若 $V$ 是复内积空间且 $T^*T = TT^* = I$（即 $T^* = T^{-1}$），则称 $T$ 为**酉算子**（unitary operator）。
    - 若 $V$ 是实内积空间且 $T^*T = TT^* = I$，则称 $T$ 为**正交算子**（orthogonal operator）。

!!! definition "定义 8.13 (酉矩阵与正交矩阵)"
    - $n \times n$ 复矩阵 $U$ 满足 $U^HU = I$（等价于 $UU^H = I$）称为**酉矩阵**（unitary matrix）。
    - $n \times n$ 实矩阵 $Q$ 满足 $Q^TQ = I$（等价于 $QQ^T = I$）称为**正交矩阵**（orthogonal matrix）。

!!! theorem "定理 8.18 (酉/正交算子的等价条件)"
    设 $V$ 是有限维内积空间，$T: V \to V$ 是线性算子。以下条件等价：

    1. $T$ 是酉（正交）算子；
    2. $T$ 保持内积：$\langle T\mathbf{u}, T\mathbf{v}\rangle = \langle \mathbf{u}, \mathbf{v}\rangle$，$\forall \mathbf{u}, \mathbf{v} \in V$；
    3. $T$ 保持范数（等距映射）：$\|T\mathbf{v}\| = \|\mathbf{v}\|$，$\forall \mathbf{v} \in V$；
    4. $T$ 将标准正交基映射为标准正交基；
    5. $T$ 的列向量构成标准正交集。

??? proof "证明"
    **(1)$\Rightarrow$(2)：** $\langle T\mathbf{u}, T\mathbf{v}\rangle = \langle \mathbf{u}, T^*T\mathbf{v}\rangle = \langle \mathbf{u}, I\mathbf{v}\rangle = \langle \mathbf{u}, \mathbf{v}\rangle$。

    **(2)$\Rightarrow$(3)：** 取 $\mathbf{u} = \mathbf{v}$。

    **(3)$\Rightarrow$(1)：** 由极化恒等式，保持范数蕴含保持内积。即 $\langle T^*T\mathbf{v}, \mathbf{v}\rangle = \langle \mathbf{v}, \mathbf{v}\rangle$ 对所有 $\mathbf{v}$ 成立，故 $\langle (T^*T - I)\mathbf{v}, \mathbf{v}\rangle = 0$。因 $T^*T - I$ 自伴，由命题 8.4 的实情形类似推理（对自伴算子成立）得 $T^*T = I$。

    **(2)$\Leftrightarrow$(4)：** 设 $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 是标准正交基，则 $\langle T\mathbf{e}_i, T\mathbf{e}_j\rangle = \langle \mathbf{e}_i, \mathbf{e}_j\rangle = \delta_{ij}$，即 $\{T\mathbf{e}_1, \ldots, T\mathbf{e}_n\}$ 是标准正交基。反之亦然。

    **(4)$\Leftrightarrow$(5)：** 直接由矩阵列向量与基的像的关系。$\blacksquare$

!!! theorem "定理 8.19 (酉/正交算子的特征值)"
    设 $T$ 是酉（正交）算子。则 $T$ 的所有特征值的模等于 $1$，即 $|\lambda| = 1$。

??? proof "证明"
    设 $T\mathbf{v} = \lambda\mathbf{v}$，$\mathbf{v} \neq \mathbf{0}$。由保范性：

    $$\|\mathbf{v}\| = \|T\mathbf{v}\| = \|\lambda\mathbf{v}\| = |\lambda|\|\mathbf{v}\|$$

    因 $\|\mathbf{v}\| > 0$，故 $|\lambda| = 1$。$\blacksquare$

!!! proposition "命题 8.5 (正交矩阵的行列式)"
    若 $Q$ 是正交矩阵，则 $\det(Q) = \pm 1$。若 $U$ 是酉矩阵，则 $|\det(U)| = 1$。

??? proof "证明"
    $Q^TQ = I$ $\Rightarrow$ $\det(Q^T)\det(Q) = 1$ $\Rightarrow$ $(\det Q)^2 = 1$ $\Rightarrow$ $\det Q = \pm 1$。

    对酉矩阵：$U^HU = I$ $\Rightarrow$ $\overline{\det U} \cdot \det U = 1$ $\Rightarrow$ $|\det U|^2 = 1$ $\Rightarrow$ $|\det U| = 1$。$\blacksquare$

!!! example "例 8.15"
    **$2 \times 2$ 正交矩阵的一般形式。** 每个 $2 \times 2$ 正交矩阵要么是旋转矩阵

    $$R_\theta = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}, \quad \det R_\theta = 1$$

    要么是反射矩阵

    $$S_\theta = \begin{pmatrix} \cos\theta & \sin\theta \\ \sin\theta & -\cos\theta \end{pmatrix}, \quad \det S_\theta = -1$$

    $\det Q = 1$ 的正交矩阵称为**特殊正交矩阵**，它们构成特殊正交群 $SO(n)$。
