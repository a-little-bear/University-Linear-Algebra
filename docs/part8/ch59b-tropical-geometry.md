# 第 59B 章 热带几何初步

<div class="context-flow" markdown>

**前置**：热带半环与热带矩阵论(Ch59A) · 凸几何(Ch64) · 代数几何基础

**本章脉络**：热带凸性 → 热带凸包 → 热带多面体 → 热带半空间与超平面 → 热带多项式与超曲面 → Kapranov 定理 → 热带 Bezout 定理 → Maslov 退量子化 → 阿米巴与 Ronkin 函数 → 热带 Grassmannian → 热带线性空间 → 系统发生学应用 → 热带代数曲线 → 镜像对称联系

**延伸**：热带几何将代数几何中的簇"热带化"为分段线性对象，为枚举几何中的深层计数问题（如 Mikhalkin 的热带曲线计数）提供了组合方法；与镜像对称和弦理论有深层联系

</div>

热带几何是近二十年来发展最为迅速的数学分支之一。它将经典代数几何中的多项式方程和代数簇"退化"为分段线性对象，在保留关键拓扑和组合信息的同时，将困难的代数问题转化为（相对）易处理的组合问题。热带几何的成功案例包括 Mikhalkin 的热带曲线计数定理（给出了一般位置曲线计数的组合公式）、热带 Grassmannian 与系统树空间的联系、以及在 Hodge 猜想等前沿问题上的新进展。

本章介绍热带几何的基本概念：热带凸性、热带多项式与超曲面、阿米巴理论、以及向更高级主题的初步探索。

---

## 59B.1 热带凸性

<div class="context-flow" markdown>

**核心问题**：热带半环中是否存在类似经典凸集和凸包的概念？热带凸性有何独特的几何结构？

</div>

!!! definition "定义 59B.1 (热带凸组合)"
    $\mathbb{T}^n$ 中点 $x_1, \ldots, x_k$ 的**热带凸组合**是

    $$\bigoplus_{i=1}^k \lambda_i \odot x_i = \Bigl(\max_{i}(\lambda_i + x_{i,1}),\, \ldots,\, \max_{i}(\lambda_i + x_{i,n})\Bigr),$$

    其中 $\lambda_1, \ldots, \lambda_k \in \mathbb{T}$ 且 $\bigoplus_i \lambda_i = \max_i \lambda_i = 0$（热带归一化条件）。

    经典凸组合 $\sum \alpha_i x_i$（$\alpha_i \geq 0$，$\sum \alpha_i = 1$）在热带化后变为取逐分量的加权最大值，权重归一化条件变为 $\max \lambda_i = 0$。

!!! definition "定义 59B.2 (热带凸集与热带凸包)"
    集合 $S \subset \mathbb{T}^n$ 称为**热带凸的**，如果 $S$ 在热带凸组合下封闭：对任意 $x, y \in S$ 和 $\lambda, \mu \in \mathbb{T}$（$\max(\lambda, \mu) = 0$），有 $\lambda \odot x \oplus \mu \odot y \in S$。

    点集 $S \subset \mathbb{T}^n$ 的**热带凸包**是 $S$ 中所有有限子集的热带凸组合的集合：

    $$\mathrm{tconv}(S) = \Bigl\{\bigoplus_{i=1}^k \lambda_i \odot x_i : x_i \in S,\, \lambda_i \in \mathbb{T},\, \max_i \lambda_i = 0,\, k \in \mathbb{N}\Bigr\}.$$

    热带凸包是包含 $S$ 的最小热带凸集。

!!! theorem "定理 59B.1 (热带凸集的基本性质)"
    (a) 热带凸集的交仍然是热带凸的。

    (b) $\mathbb{T}^n$ 中的热带凸集在投影空间 $\mathbb{TP}^{n-1} = (\mathbb{T}^n \setminus \{(-\infty)^n\})/\sim$（模去热带标量乘法）中是**分段线性**的。

    (c) 有限点集的热带凸包是热带凸多面体（定义见 59B.2 节）。

!!! example "例 59B.1"
    在 $\mathbb{T}^2/\mathbb{R}\mathbf{1}$ 中（即 $\mathbb{R}^2$ 模去平移 $(t, t)$，用第一坐标减第二坐标来参数化），考虑两点 $p = (0, 0)$ 和 $q = (2, 0)$ 的热带线段。

    热带线段 $\overline{pq}^{\mathrm{trop}}$ = $\{(\max(\lambda, \mu+2), \max(\lambda, \mu)) : \max(\lambda, \mu) = 0\}$。

    当 $\lambda = 0, \mu \leq 0$ 时（$\lambda$ 主导）：$(\max(0, \mu+2), 0)$。
    - 若 $\mu + 2 \leq 0$（即 $\mu \leq -2$）：点为 $(0, 0) = p$。
    - 若 $\mu + 2 > 0$（即 $-2 < \mu \leq 0$）：点为 $(\mu+2, 0)$，从 $(0,0)$ 到 $(2,0)$。

    当 $\mu = 0, \lambda \leq 0$ 时（$\mu$ 主导）：$(\max(\lambda, 2), 0)$。
    - 若 $\lambda < 2$：点为 $(2, 0) = q$。

    在商空间中，热带线段由一段水平线段 $\{(t, 0) : 0 \leq t \leq 2\}$ 构成——与经典线段一致。但对更复杂的配置，热带线段会展现出分段线性的弯折行为。

---

## 59B.2 热带多面体

<div class="context-flow" markdown>

**核心问题**：热带凸多面体有怎样的组合结构？与经典多面体理论有何异同？

</div>

!!! definition "定义 59B.3 (热带多面体)"
    **热带凸多面体**（tropical polytope）是有限点集在 $\mathbb{TP}^{n-1}$ 中的热带凸包。等价地，它是有限多个热带半空间的交。

    $\mathbb{TP}^{d-1}$ 中由 $n$ 个点 $v_1, \ldots, v_n$ 生成的热带多面体是

    $$P = \mathrm{tconv}(v_1, \ldots, v_n) = \Bigl\{\bigoplus_{i=1}^n \lambda_i \odot v_i : \max_i \lambda_i = 0\Bigr\}.$$

!!! theorem "定理 59B.2 (热带多面体的结构，Develin-Sturmfels)"
    设 $P = \mathrm{tconv}(v_1, \ldots, v_n) \subset \mathbb{TP}^{d-1}$。则：

    (a) $P$ 是普通拓扑空间中的**多面体复形**（polyhedral complex），其面由**类型**（type）标记。

    (b) 点 $x \in P$ 的类型 $T(x) = (T_1, \ldots, T_d)$ 定义为：$j \in T_i$ 当且仅当 $\lambda_j + (v_j)_i = x_i$（第 $j$ 个生成元的第 $i$ 坐标是"活跃"的）。

    (c) 类型相同的点构成 $P$ 的一个（开）面。

    (d) 面偏序给出热带多面体的组合结构，它与**正则三角剖分**和**混合细分**有深层联系。

!!! definition "定义 59B.4 (热带顶点)"
    热带多面体 $P$ 的**热带顶点**（pseudovertex 或 tropical vertex）是不能写成 $P$ 中其他两点的严格热带凸组合的点。

    Develin 和 Sturmfels 证明：$\mathrm{tconv}(v_1, \ldots, v_n)$ 的热带顶点集是 $\{v_1, \ldots, v_n\}$ 的某个子集（经可能的热带投影简化后）。

!!! example "例 59B.2"
    在 $\mathbb{TP}^1$（热带投影线，可以参数化为 $\mathbb{R}$）中，三点 $a, b, c \in \mathbb{R}$（$a < b < c$）的热带凸包就是闭区间 $[a, c]$——与经典凸包相同。

    在 $\mathbb{TP}^2$ 中，三点 $p = (0,0,0), q = (a,0,0), r = (0,b,0)$（在商空间中用前两个坐标减第三个来参数化）的热带凸包是一个"热带三角形"，其形状依赖于 $a, b$ 的值。当 $a, b > 0$ 时，热带三角形有 6 条边（3 条有界边和 3 对射线），与经典三角形的 3 条边不同。

---

## 59B.3 热带半空间与超平面

!!! definition "定义 59B.5 (热带超平面)"
    由热带线性形式 $f(x) = \bigoplus_{i=1}^n a_i \odot x_i = \max_i(a_i + x_i)$ 定义的**热带超平面**是 $f$ 不光滑的点集：

    $$H_f = \{x \in \mathbb{R}^n : \text{最大值 } \max_i(a_i + x_i) \text{ 由至少两个 } i \text{ 达到}\}.$$

    在 $\mathbb{R}^n/\mathbb{R}\mathbf{1}$ 中，$H_f$ 是一个维数为 $n-2$ 的分段线性复形。

!!! definition "定义 59B.6 (热带半空间)"
    热带超平面 $H_f$ 将补集分成若干连通分量，每个分量称为一个**扇区**（sector）。扇区 $S_i$ 对应于 $\max_j(a_j + x_j)$ 由 $j = i$ 唯一达到的区域：

    $$S_i = \{x : a_i + x_i > a_j + x_j,\, \forall j \neq i\}.$$

    **热带半空间**是一个或多个扇区的闭包的并，具有形式

    $$\{x : a_i + x_i \geq a_j + x_j\} = \{x : x_i - x_j \geq a_j - a_i\}.$$

!!! theorem "定理 59B.3 (热带凸集的外部描述)"
    热带凸集（在 $\mathbb{TP}^{d-1}$ 中）是**分段线性**的，可以表示为有限多个热带半空间的交：

    $$C = \bigcap_{k} \{x \in \mathbb{TP}^{d-1} : a_k \odot x_{i_k} \oplus b_k \odot x_{j_k} \leq c_k \odot x_{\ell_k}\}.$$

    这是经典 Minkowski-Weyl 定理（凸多面体的内部描述与外部描述的等价性）的热带类比。

!!! example "例 59B.3"
    在 $\mathbb{R}^3/\mathbb{R}\mathbf{1}$ 中，热带线性形式 $f(x) = \max(x_1, x_2, x_3)$ 定义的热带超平面是

    $$H = \{(x_1, x_2, x_3) : \text{某两个坐标相等且不小于第三个}\}.$$

    在商空间 $\mathbb{R}^2$（令 $x_3 = 0$）中，$H$ 由三条射线组成：

    - $\{(t, t) : t \geq 0\}$（$x_1 = x_2 \geq x_3 = 0$）
    - $\{(t, 0) : t \leq 0\}$（$x_1 \leq 0 = x_2 = x_3$，即 $x_2 = x_3 \geq x_1$）
    - $\{(0, t) : t \leq 0\}$（$x_1 = x_3 \geq x_2$）

    三条射线从原点出发，形成"Y"形，将平面分成三个扇区。

---

## 59B.4 热带多项式与超曲面

<div class="context-flow" markdown>

**核心问题**：经典代数几何中的概念如何"热带化"？热带簇有怎样的组合结构？

</div>

!!! definition "定义 59B.7 (热带多项式)"
    $\mathbb{T}$ 上的**热带多项式**是

    $$p(x_1, \ldots, x_n) = \bigoplus_{i \in S} a_i \odot x_1^{\odot i_1} \odot \cdots \odot x_n^{\odot i_n} = \max_{i \in S}\Bigl(a_i + \sum_{k=1}^n i_k x_k\Bigr),$$

    其中 $S$ 是有限支撑集，$i = (i_1, \ldots, i_n)$ 是多重指标，$a_i \in \mathbb{T}$。热带多项式是 $(x_1, \ldots, x_n)$ 的**分段线性凸函数**（是有限多个仿射函数的逐点最大值）。

!!! definition "定义 59B.8 (热带超曲面)"
    热带多项式 $p$ 的**热带超曲面**（或热带簇）$\mathcal{T}(p)$ 是 $p$ **不光滑的点**的集合，即最大值由两个或更多项同时达到的点：

    $$\mathcal{T}(p) = \{x \in \mathbb{R}^n : \text{最大值在 } p(x) \text{中由至少两项达到}\}.$$

    $\mathcal{T}(p)$ 是 $\mathbb{R}^n$ 中的一个**加权平衡多面体复形**（weighted balanced polyhedral complex）。

!!! theorem "定理 59B.4 (热带超曲面的结构定理)"
    设 $p$ 是支撑集为 $S$ 的热带多项式。则：

    (a) $\mathcal{T}(p)$ 是 $\mathbb{R}^n$ 中余维 1 的多面体复形。

    (b) $\mathcal{T}(p)$ 的每个极大面携带正整数**权重**（重数），满足**平衡条件**（balancing condition）：在每个余维 2 的面上，相邻极大面的加权法向量之和为零。

    (c) $\mathcal{T}(p)$ 与 $p$ 的**Newton 多面体** $\mathrm{New}(p) = \mathrm{conv}(S)$ 的**正则细分**对偶。

!!! example "例 59B.4"
    考虑热带二次多项式

    $$p(x) = 3 \odot x^{\odot 2} \oplus 5 \odot x \oplus 2 = \max(3 + 2x,\, 5 + x,\, 2).$$

    三条线 $y_1 = 3+2x$，$y_2 = 5+x$，$y_3 = 2$ 的交点：

    - $y_1 = y_2$：$3+2x = 5+x \Rightarrow x = 2$。
    - $y_2 = y_3$：$5+x = 2 \Rightarrow x = -3$。
    - 在 $x = 2$ 处：$y_1 = y_2 = 7 > y_3 = 2$（活跃项为 $y_1, y_2$）。
    - 在 $x = -3$ 处：$y_2 = y_3 = 2 > y_1 = -3$（活跃项为 $y_2, y_3$）。

    热带根为 $x = 2$ 和 $x = -3$（对应两个"弯折点"）。类比：经典二次多项式有两个根，热带二次多项式也恰好有两个根（计重数）。

!!! example "例 59B.5"
    **热带直线。** $\mathbb{R}^2$ 中的热带直线由热带一次多项式定义：

    $$L(x,y) = a \odot x \oplus b \odot y \oplus c = \max(a+x, b+y, c).$$

    热带直线 $\mathcal{T}(L)$ 由三条射线构成（从拐点 $(c-a, c-b)$ 出发），分别向 $(1,0)$、$(0,1)$、$(-1,-1)$ 方向延伸。这是一个"Y"形图。

    **平衡条件**：三个方向的原始向量之和 $(1,0) + (0,1) + (-1,-1) = (0,0)$，满足平衡。

---

## 59B.5 Kapranov 定理

!!! theorem "定理 59B.5 (Kapranov 定理)"
    设 $\mathbb{K}$ 是具有非 Archimedean 赋值 $\mathrm{val}: \mathbb{K}^* \to \mathbb{R}$ 的代数闭域（例如 Puiseux 级数域 $\mathbb{C}\{\{t\}\}$，赋值为 $t$ 的阶）。

    设 $f(x_1, \ldots, x_n) = \sum_{i \in S} c_i x^i \in \mathbb{K}[x_1, \ldots, x_n]$ 是多项式，$V(f) = \{z \in (\mathbb{K}^*)^n : f(z) = 0\}$ 是超曲面。

    则 $V(f)$ 在赋值映射下的像 $\mathrm{val}(V(f)) = \{(\mathrm{val}(z_1), \ldots, \mathrm{val}(z_n)) : z \in V(f)\}$ 恰好等于热带超曲面 $\mathcal{T}(\mathrm{trop}(f))$，其中

    $$\mathrm{trop}(f)(w) = \max_{i \in S}(\mathrm{val}(c_i) + i \cdot w).$$

??? proof "证明（思路）"
    **"$\subset$" 方向。** 设 $z \in V(f)$，$w = \mathrm{val}(z)$。由 $f(z) = 0$，即 $\sum_i c_i z^i = 0$。

    设 $M = \max_i \mathrm{val}(c_i z^i) = \max_i(\mathrm{val}(c_i) + i \cdot w) = \mathrm{trop}(f)(w)$。

    如果最大值仅由一个项 $i^*$ 达到，则 $\mathrm{val}(c_{i^*} z^{i^*}) > \mathrm{val}(c_i z^i)$ 对所有 $i \neq i^*$。由非 Archimedean 赋值的超度量性质，$\mathrm{val}(\sum_i c_i z^i) = M > -\infty$，即 $f(z) \neq 0$，矛盾。因此 $w \in \mathcal{T}(\mathrm{trop}(f))$。

    **"$\supset$" 方向。** 设 $w \in \mathcal{T}(\mathrm{trop}(f))$。需要构造 $z \in (\mathbb{K}^*)^n$ 使得 $\mathrm{val}(z) = w$ 且 $f(z) = 0$。

    当 $\mathbb{K}$ 代数闭时（特别是 Puiseux 级数域），可以利用代数闭性和 Hensel 引理来提升。直觉上，$w \in \mathcal{T}$ 意味着 $f$ 的主项有消去，因此"接近"零点。代数闭性保证了这种"接近"可以精确到真正的零点。

    精确证明需要利用赋值代数和 Newton 多面体的精细分析。对一元多项式的情形，证明更为直接——热带根恰好是 Newton 多边形的斜率，由 Newton-Puiseux 定理，这些斜率对应经典根的赋值。$\blacksquare$

Kapranov 定理是热带几何的基石。它说明热带化不仅仅是形式上的类比，而是精确地捕捉了代数簇在赋值映射下的"影子"。

---

## 59B.6 热带 Bezout 定理

!!! theorem "定理 59B.6 (热带 Bezout 定理)"
    设 $p, q$ 是 $\mathbb{T}$ 上两个一元热带多项式，次数分别为 $m, n$。则 $\mathcal{T}(p) \cap \mathcal{T}(q)$ 的交点数（计重数）至多为 $mn$。

    更一般地，$\mathbb{R}^2$ 中两条热带曲线 $\mathcal{T}(f)$ 和 $\mathcal{T}(g)$（次数分别为 $m, n$）在"一般位置"时的稳定交点数（计适当的热带交叉重数）恰好为 $mn$。

    热带交叉重数的定义：设 $\mathcal{T}(f)$ 和 $\mathcal{T}(g)$ 横截相交于点 $p$，$p$ 处 $\mathcal{T}(f)$ 的边 $e_1$（权重 $w_1$，原始方向 $v_1$）与 $\mathcal{T}(g)$ 的边 $e_2$（权重 $w_2$，原始方向 $v_2$）交于 $p$。则 $p$ 处的交叉重数为

    $$m(p) = w_1 w_2 \cdot |\det(v_1, v_2)|.$$

!!! example "例 59B.6"
    两条热带直线 $L_1 = \max(x, y, 0)$ 和 $L_2 = \max(x+1, y+2, 3)$ 在 $\mathbb{R}^2$ 中的交。

    $L_1$ 的拐点为 $(0, 0)$，$L_2$ 的拐点为 $(3-1, 3-2) = (2, 1)$。两个"Y"形图。

    由热带 Bezout 定理，$1 \times 1 = 1$，应恰有 1 个交叉点（计重数）。

---

## 59B.7 Maslov 退量子化

<div class="context-flow" markdown>

**核心问题**：热带代数与经典代数之间有什么深层的哲学联系？热带化过程是否有物理直觉？

</div>

!!! definition "定义 59B.9 (Maslov 退量子化)"
    **Maslov 退量子化**（dequantization）提供了从经典代数到热带代数的连续变形。

    对 $\hbar > 0$，在 $\mathbb{R}_{> 0}$ 上定义变形运算：

    $$a \oplus_\hbar b = \hbar \log(e^{a/\hbar} + e^{b/\hbar}), \quad a \odot_\hbar b = a + b.$$

    当 $\hbar > 0$ 时，$(\mathbb{R}, \oplus_\hbar, \odot_\hbar)$ 与 $(\mathbb{R}_{>0}, +, \times)$ 同构（通过映射 $x \mapsto e^{x/\hbar}$）。

    关键观察：当 $\hbar \to 0^+$ 时，

    $$\lim_{\hbar \to 0^+} \hbar \log(e^{a/\hbar} + e^{b/\hbar}) = \max(a, b).$$

    因此 $\oplus_\hbar \to \oplus = \max$，而 $\odot_\hbar = \odot = +$ 不变。

!!! theorem "定理 59B.7 (Maslov 退量子化原理)"
    热带半环 $(\mathbb{R} \cup \{-\infty\}, \max, +)$ 是经典半环 $(\mathbb{R}_{\geq 0}, +, \times)$ 的**退量子化极限**：

    $$(\mathbb{R}_{\geq 0}, +, \times) \xrightarrow[\hbar \to 0]{} (\mathbb{R} \cup \{-\infty\}, \max, +).$$

    这个过程可以理解为：参数 $\hbar$ 类似于物理中的 Planck 常数。经典数学（$+, \times$）对应"量子"情形（$\hbar > 0$），而热带数学（$\max, +$）对应"经典极限"（$\hbar \to 0$）。

    这种视角解释了为什么热带几何中的对象是分段线性的——它们是经典代数对象在 $\hbar \to 0$ 极限下的"骨架"或"经典影子"。

??? proof "证明"
    只需验证极限。设 $a \geq b$（不失一般性），则

    $$\hbar \log(e^{a/\hbar} + e^{b/\hbar}) = \hbar \log\bigl(e^{a/\hbar}(1 + e^{(b-a)/\hbar})\bigr) = a + \hbar \log(1 + e^{(b-a)/\hbar}).$$

    由于 $b - a \leq 0$，当 $\hbar \to 0^+$ 时，$e^{(b-a)/\hbar} \to 0$（当 $b < a$），因此

    $$\hbar \log(1 + e^{(b-a)/\hbar}) \to 0.$$

    极限为 $a = \max(a, b)$。当 $a = b$ 时，$\hbar \log(1 + 1) = \hbar \log 2 \to 0$，极限仍为 $a = \max(a,b)$。$\blacksquare$

!!! example "例 59B.7"
    考虑经典多项式 $f(x) = e^{3/\hbar} x^2 + e^{5/\hbar} x + e^{2/\hbar}$。在对数坐标（$x = e^{w/\hbar}$）下，

    $$f(e^{w/\hbar}) = e^{(3+2w)/\hbar} + e^{(5+w)/\hbar} + e^{2/\hbar}.$$

    $$\hbar \log f(e^{w/\hbar}) = \hbar \log\bigl(e^{(3+2w)/\hbar} + e^{(5+w)/\hbar} + e^{2/\hbar}\bigr) \xrightarrow{\hbar \to 0} \max(3+2w, 5+w, 2).$$

    这就是例 59B.4 中的热带二次多项式。经典多项式的零点在退量子化极限下变为热带超曲面。

---

## 59B.8 阿米巴

<div class="context-flow" markdown>

**核心问题**：代数簇的阿米巴是什么？它如何连接经典代数几何与热带几何？

</div>

!!! definition "定义 59B.10 (阿米巴)"
    设 $V \subset (\mathbb{C}^*)^n$ 是代数簇。$V$ 的**阿米巴**（amoeba）是对数映射的像：

    $$\mathcal{A}(V) = \{(\log|z_1|, \ldots, \log|z_n|) : (z_1, \ldots, z_n) \in V\} \subset \mathbb{R}^n.$$

    阿米巴这一概念由 Gelfand、Kapranov 和 Zelevinsky 在 1994 年引入。

!!! definition "定义 59B.11 (Ronkin 函数)"
    与阿米巴 $\mathcal{A}(V(f))$ 相关的 **Ronkin 函数** $N_f: \mathbb{R}^n \to \mathbb{R}$ 定义为

    $$N_f(x) = \frac{1}{(2\pi i)^n} \int_{|z_i| = e^{x_i}} \log|f(z)| \frac{dz_1}{z_1} \cdots \frac{dz_n}{z_n}.$$

    Ronkin 函数是凸函数，且在阿米巴补集的每个连通分量上是仿射的。

!!! theorem "定理 59B.8 (阿米巴与热带极限)"
    设 $f_t(z) = \sum_{i} c_i t^{a_i} z^i$（$t > 0$ 为参数），$V_t = V(f_t) \subset (\mathbb{C}^*)^n$。定义**缩放阿米巴**

    $$\mathcal{A}_t(V_t) = \{(\log_t|z_1|, \ldots, \log_t|z_n|) : z \in V_t\}.$$

    则当 $t \to \infty$ 时（Hausdorff 度量意义下），

    $$\lim_{t \to \infty} \mathcal{A}_t(V_t) = \mathcal{T}(\mathrm{trop}(f)),$$

    其中 $\mathrm{trop}(f)(w) = \max_i(a_i + i \cdot w)$。

    换言之，热带超曲面是代数簇阿米巴在参数趋于极限时的"骨架"。

!!! example "例 59B.8"
    直线 $f(z_1, z_2) = z_1 + z_2 + 1 = 0$ 在 $(\mathbb{C}^*)^2$ 中的阿米巴。

    令 $z_1 = e^{x_1 + i\theta_1}$，$z_2 = e^{x_2 + i\theta_2}$。条件 $z_1 + z_2 + 1 = 0$ 的模约束给出阿米巴的形状：一个具有三个"触手"的区域，向 $(-\infty, 0)$、$(0, -\infty)$、$(+\infty, +\infty)$ 方向无限延伸。

    阿米巴的补集有三个凸连通分量，分别对应 Newton 多面体 $\mathrm{conv}\{(1,0), (0,1), (0,0)\}$ 的三个顶点。

    热带极限 $\mathcal{T}(\max(x_1, x_2, 0))$ 是三条射线从原点出发的"Y"形图，恰好是阿米巴的"骨架"。

---

## 59B.9 热带 Grassmannian

<div class="context-flow" markdown>

**核心问题**：经典 Grassmannian 的热带化是什么？它有怎样的组合结构？

</div>

!!! definition "定义 59B.12 (热带 Grassmannian)"
    经典 Grassmannian $\mathrm{Gr}(k, n)$ 参数化 $\mathbb{K}^n$ 中的 $k$ 维子空间。通过 Pluecker 嵌入，$\mathrm{Gr}(k, n)$ 可以视为射影空间 $\mathbb{P}^{\binom{n}{k}-1}$ 中由 Pluecker 关系定义的子簇。

    **热带 Grassmannian** $\mathrm{TGr}(k, n)$ 是 Pluecker 关系的热带化所定义的热带簇。具体地，热带 Pluecker 坐标是 $p_I \in \mathbb{T}$（$I \subset [n], |I| = k$），满足热带 Pluecker 关系。

    对 $\mathrm{TGr}(2, n)$（热带直线的空间），热带 Pluecker 关系为：对所有 $i < j < k < l$，

    $$\max(p_{ij} + p_{kl},\, p_{ik} + p_{jl},\, p_{il} + p_{jk})$$

    中最大值由至少两项达到。

!!! theorem "定理 59B.9 (Speyer-Sturmfels, 2004)"
    (a) $\mathrm{TGr}(2, n)$ 与空间 $\mathbb{R}^{\binom{n}{2}}/\mathbb{R}\mathbf{1}$ 中的**系统树空间**（space of phylogenetic trees）同胚：$\mathrm{TGr}(2, n)$ 的点一一对应于 $n$ 叶有权度量树。

    (b) $\mathrm{TGr}(2, n)$ 是一个纯维数 $2n - 3$ 的多面体扇（polyhedral fan），其极大锥与 $n$ 叶树的拓扑类型一一对应。

    (c) 一般的 $\mathrm{TGr}(k, n)$ 是 Dressian $\mathrm{Dr}(k, n)$ 的子集：$\mathrm{TGr}(k, n) \subset \mathrm{Dr}(k, n)$，且一般而言是真子集。

!!! example "例 59B.9"
    **$\mathrm{TGr}(2, 4)$。** 有 $\binom{4}{2} = 6$ 个 Pluecker 坐标 $p_{12}, p_{13}, p_{14}, p_{23}, p_{24}, p_{34}$。

    热带 Pluecker 关系只有一个（对应 $(i,j,k,l) = (1,2,3,4)$）：

    $$\max(p_{12}+p_{34},\, p_{13}+p_{24},\, p_{14}+p_{23})$$

    中最大值至少由两项达到。

    在商空间 $\mathbb{R}^6/\mathbb{R}\mathbf{1}$ 中，$\mathrm{TGr}(2,4)$ 是 3 维的多面体扇，有 3 个极大锥，对应 4 叶树的 3 种拓扑类型（3 种方式将 4 个叶子配对）。

---

## 59B.10 热带线性空间

!!! definition "定义 59B.13 (热带线性空间)"
    $\mathbb{TP}^{n-1}$ 中的**热带线性空间**是由热带线性形式的热带超曲面的适当交来定义的。等价地，$k$ 维热带线性空间可以用热带 Pluecker 坐标来参数化。

    更精确地，给定热带 Pluecker 向量 $p \in \mathrm{TGr}(k, n)$，对应的热带线性空间 $L_p$ 是

    $$L_p = \{x \in \mathbb{TP}^{n-1} : \text{对所有 } |J| = n-k,\, \max_{i \notin J}(p_{J \cup \{i\}} + x_i) \text{ 由至少两个 } i \text{ 达到}\}.$$

!!! theorem "定理 59B.10 (热带线性空间的结构)"
    热带线性空间 $L_p$（对应 $p \in \mathrm{TGr}(k, n)$）是 $\mathbb{TP}^{n-1}$ 中的一个纯维数 $k-1$ 的平衡多面体复形。

    当 $p$ 来自经典 Grassmannian 的赋值（即 $p \in \mathrm{TGr}(k,n)$ 而非仅在 Dressian $\mathrm{Dr}(k,n)$ 中时），$L_p$ 具有更好的组合性质——它是一个热带簇，满足 Kapranov 定理的多维推广。

---

## 59B.11 系统发生学应用

<div class="context-flow" markdown>

**核心问题**：热带几何如何应用于进化生物学中的系统树重构？

</div>

!!! definition "定义 59B.14 (系统发生树空间)"
    $n$ 个物种的**系统发生树空间**（phylogenetic tree space）$\mathcal{T}_n$ 是所有 $n$ 叶有根（或无根）加权树的参数空间。

    无根树的表示：每棵 $n$ 叶无根树 $T$ 由其拓扑结构（树形图）和边权 $w_e > 0$ 确定。树诱导了叶子集 $[n]$ 上的度量 $d_{ij} = \sum_{e \in \mathrm{path}(i,j)} w_e$。

!!! theorem "定理 59B.11 (系统树空间与热带 Grassmannian)"
    (a) $n$ 叶有权度量树的空间（模去全局平移）与 $\mathrm{TGr}(2, n)$ 同胚（Speyer-Sturmfels）。

    (b) 给定叶子间距离矩阵 $D = (d_{ij})$，对应的热带 Pluecker 坐标为 $p_{ij} = -d_{ij}$（取负，因为树度量对应 min-plus 结构）。

    (c) 四点条件：$(d_{ij})$ 是树度量当且仅当对所有四元组 $\{i,j,k,l\}$，

    $$\max(d_{ij}+d_{kl},\, d_{ik}+d_{jl},\, d_{il}+d_{jk})$$

    中最大值由至少两项达到。这恰是热带 Pluecker 关系。

!!! example "例 59B.10"
    考虑 4 个物种的距离矩阵：

    $$D = \begin{pmatrix}0 & 3 & 5 & 5\\3 & 0 & 4 & 4\\5 & 4 & 0 & 2\\5 & 4 & 2 & 0\end{pmatrix}.$$

    热带 Pluecker 坐标：$p_{12} = -3, p_{13} = -5, p_{14} = -5, p_{23} = -4, p_{24} = -4, p_{34} = -2$。

    验证四点条件（$(1,2,3,4)$）：

    $$\max(p_{12}+p_{34}, p_{13}+p_{24}, p_{14}+p_{23}) = \max(-3+(-2), -5+(-4), -5+(-4)) = \max(-5, -9, -9) = -5.$$

    最大值由 $p_{12}+p_{34} = -5$ 唯一达到。等等，这不满足"至少两项达到"的条件。

    重新检查：这意味着 $D$ 不是精确的树度量。实际上这正是四点条件中"至少两项等于最大值"的要求。如果取 $d_{13} = d_{14} = 5$，$d_{23} = d_{24} = 4$，那么 $p_{13}+p_{24} = p_{14}+p_{23} = -9$，两者相等但不是最大值。调整 $d_{34} = 4$：$p_{12}+p_{34} = -7, p_{13}+p_{24} = -9, p_{14}+p_{23} = -9$，仍是 $p_{12}+p_{34}$ 唯一最大。

    取精确树度量的例子：设树为 $((1,2), (3,4))$，内部边权 $a$，叶边权 $w_1, w_2, w_3, w_4$。则 $d_{12} = w_1+w_2, d_{34} = w_3+w_4, d_{13} = w_1+a+w_3, d_{24} = w_2+a+w_4$。此时 $d_{13}+d_{24} = d_{14}+d_{23} = w_1+w_2+w_3+w_4+2a$，而 $d_{12}+d_{34} = w_1+w_2+w_3+w_4$。

    四点条件：$\max(d_{12}+d_{34}, d_{13}+d_{24}, d_{14}+d_{23}) = d_{13}+d_{24} = d_{14}+d_{23}$，后两项相等且为最大值。满足条件。树的拓扑是 $\{1,2\}|\{3,4\}$。

---

## 59B.12 热带代数曲线

<div class="context-flow" markdown>

**核心问题**：热带曲线有怎样的组合结构？如何用热带曲线计算枚举几何中的不变量？

</div>

!!! definition "定义 59B.15 (热带曲线)"
    $\mathbb{R}^2$ 中的**热带曲线**是由一个度为 $d$ 的热带多项式 $p(x, y) = \max_{(i,j) \in S}(a_{ij} + ix + jy)$ 的热带超曲面 $\mathcal{T}(p)$。

    热带曲线 $\Gamma = \mathcal{T}(p)$ 是一个加权平衡图（weighted balanced graph），其特征包括：

    - **有界边**：连接两个内部顶点的线段。
    - **无界边**（射线）：从内部顶点向无穷远延伸。
    - **权重**：每条边携带正整数权重。
    - **平衡条件**：在每个内部顶点，相邻边的加权方向向量之和为零。
    - **亏格**：$g(\Gamma) = $ 图 $\Gamma$ 的第一 Betti 数（环数）。

!!! theorem "定理 59B.12 (热带曲线的亏格-次数关系)"
    度为 $d$ 的光滑热带曲线（Newton 多面体为 $d\Delta_2$，其中 $\Delta_2$ 是标准单纯形）的亏格满足

    $$g \leq \frac{(d-1)(d-2)}{2},$$

    等号在"一般位置"时成立。这与经典代数曲线的亏格公式完全一致。

!!! theorem "定理 59B.13 (Mikhalkin 对应定理, 2005)"
    设 $N_d$ 是通过 $3d - 1$ 个一般位置点的 $\mathbb{P}^2$ 中度为 $d$ 的有理曲线数目（经典枚举几何中的 Gromov-Witten 不变量）。则

    $$N_d = \sum_{\Gamma} \mathrm{mult}(\Gamma),$$

    其中求和遍历通过对应 $3d-1$ 个热带化点的所有度为 $d$、亏格为 $0$ 的热带曲线 $\Gamma$，$\mathrm{mult}(\Gamma)$ 是 $\Gamma$ 的**热带重数**（由各顶点处的局部重数之积定义）。

    这一定理将困难的枚举几何问题（计数代数曲线）归结为组合问题（计数热带曲线），是热带几何最惊人的应用之一。

!!! example "例 59B.11"
    **$d = 1$**：通过 2 个点的直线数 $N_1 = 1$。热带化：恰好一条热带直线通过 2 个一般位置的热带点。

    **$d = 2$**：通过 5 个点的圆锥曲线数 $N_2 = 1$。

    **$d = 3$**：通过 8 个点的有理三次曲线数 $N_3 = 12$。Mikhalkin 的定理将这一经典结果归结为对热带三次曲线的组合计数。

---

## 59B.13 与镜像对称和枚举几何的联系

<div class="context-flow" markdown>

**核心问题**：热带几何在现代数学物理中扮演什么角色？

</div>

!!! definition "定义 59B.16 (SYZ 猜想与热带几何)"
    **Strominger-Yau-Zaslow (SYZ) 猜想**（1996）是镜像对称的几何解释：Calabi-Yau 流形 $X$ 的镜像 $\check{X}$ 可以通过对 $X$ 的 Lagrangian 环面纤维化进行对偶来获得。

    热带几何在 SYZ 猜想中的角色：

    (a) Calabi-Yau 流形在退化极限下（"大复结构极限"）趋于一个热带流形（分段线性空间）。

    (b) 热带流形的组合结构编码了镜像对称的关键数据。

    (c) Gross-Siebert 纲领利用热带几何系统地构造镜像对并证明镜像对称的若干情形。

!!! theorem "定理 59B.14 (热带几何与枚举不变量)"
    热带方法在枚举几何中的主要成果包括：

    (a) **Mikhalkin 对应定理**（定理 59B.13）：将平面曲线计数转化为热带曲线计数。

    (b) **热带 Hurwitz 数**：分支覆盖的计数问题可以热带化，得到组合公式。

    (c) **热带 Welschinger 不变量**：实代数几何中的枚举不变量（实有理曲线计数）也可以通过热带方法计算。

    (d) **高维推广**：Nishinou-Siebert 将 Mikhalkin 的对应定理推广到环面簇上的高维曲线计数。

这些结果表明，热带几何不仅是经典代数几何的"影子"，更是一个独立而强大的计算工具。

---

## 59B.14 习题

**热带凸性**

!!! example "习题 59B.1"
    在 $\mathbb{T}^3/\mathbb{R}\mathbf{1}$ 中，计算点 $p = (0, 0, 0)$、$q = (3, 0, 0)$、$r = (0, 2, 0)$ 的热带凸包的所有类型区域。

!!! example "习题 59B.2"
    证明：热带凸集的交仍是热带凸的。（提示：直接验证定义。）

!!! example "习题 59B.3"
    在 $\mathbb{TP}^1$ 中（即 $\mathbb{R}$），证明任意有限点集的热带凸包就是经典的闭区间 $[\min, \max]$。

**热带多项式与超曲面**

!!! example "习题 59B.4"
    画出以下热带多项式的超曲面：

    (a) $p(x, y) = \max(x, y, 0)$（热带直线）。

    (b) $p(x, y) = \max(2x, x+y, 2y, x, y, 0)$（热带圆锥曲线）。

    标注各边的权重和方向。

!!! example "习题 59B.5"
    证明：一元度为 $d$ 的热带多项式恰好有 $d$ 个根（计重数）。（提示：分析分段线性函数的弯折点。）

!!! example "习题 59B.6"
    对热带多项式 $p(x) = 1 \odot x^{\odot 3} \oplus 2 \odot x^{\odot 2} \oplus 0 \odot x \oplus (-1) = \max(1+3x, 2+2x, x, -1)$，找出所有热带根及其重数。

**Kapranov 定理与阿米巴**

!!! example "习题 59B.7"
    考虑多项式 $f(x) = t^3 x^2 + t^5 x + t^2$（$t$ 为 Puiseux 级数的形式参数）。

    (a) 写出对应的热带多项式。

    (b) 找出热带根。

    (c) 利用 Kapranov 定理，说明 $f$ 的经典根的赋值。

!!! example "习题 59B.8"
    画出曲线 $f(z_1, z_2) = 1 + z_1 + z_2 + z_1 z_2 = 0$ 的阿米巴的大致形状。其热带极限是什么？

**热带 Grassmannian 与系统发生学**

!!! example "习题 59B.9"
    对 5 个物种，给出一个树度量矩阵 $D$，验证所有四点条件成立。画出对应的系统发生树。

!!! example "习题 59B.10"
    证明：$n$ 叶无根二叉树的拓扑类型数为 $(2n-5)!! = 1 \cdot 3 \cdot 5 \cdots (2n-5)$。说明这与 $\mathrm{TGr}(2, n)$ 的极大锥数的关系。

!!! example "习题 59B.11"
    $\mathrm{TGr}(2, 5)$ 中有多少种树拓扑类型？对应多少个极大锥？

**Maslov 退量子化**

!!! example "习题 59B.12"
    验证 Maslov 退量子化中 $\oplus_\hbar$ 运算的结合律和交换律。

!!! example "习题 59B.13"
    对三个正实数 $a, b, c$，计算 $\lim_{\hbar \to 0^+} \hbar \log(e^{a/\hbar} + e^{b/\hbar} + e^{c/\hbar})$，并解释其几何意义。

**热带曲线**

!!! example "习题 59B.14"
    构造一条度为 $3$、亏格为 $1$ 的热带曲线（在 $\mathbb{R}^2$ 中）。验证平衡条件在每个顶点成立。

!!! example "习题 59B.15"
    利用 Mikhalkin 对应定理的思想，通过热带方法验证 $N_2 = 1$（通过 5 个一般位置点恰有 1 条圆锥曲线）。

**理论证明**

!!! example "习题 59B.16"
    证明热带超曲面 $\mathcal{T}(p)$ 是闭集。（提示：$\mathcal{T}(p)$ 是两个凸函数差为零的集合。）

!!! example "习题 59B.17"
    证明平衡条件：设 $p(x,y)$ 是热带多项式，$v$ 是 $\mathcal{T}(p)$ 的一个内部顶点。证明从 $v$ 出发的各边的加权原始方向向量之和为零向量。（提示：利用热带多项式与 Newton 多面体的对偶关系。）

!!! example "习题 59B.18"
    设 $f(x,y) = \max(a + 2x, b + x + y, c + 2y, d + x, e + y, f)$ 是一个度为 $2$ 的热带多项式。当 $a = b = c = d = e = f = 0$ 时，画出 $\mathcal{T}(f)$ 并标注所有顶点、边和权重。验证亏格公式 $g \leq (2-1)(2-2)/2 = 0$。

---

**本章要点总结：**

1. 热带凸性将经典凸几何中的线性组合替换为热带凸组合（加权取最大值），热带凸集是分段线性的。
2. 热带多面体具有丰富的组合结构，由类型（type）系统描述，与正则三角剖分对偶。
3. 热带超平面是分段线性的"Y"形或更高维类比物，热带半空间是它的闭半部分。
4. 热带多项式是分段线性凸函数，热带超曲面是其不光滑点集——加权平衡多面体复形。
5. Kapranov 定理建立了热带簇与代数簇赋值像之间的精确对应，是热带几何的理论基石。
6. 热带 Bezout 定理给出热带曲线的交点计数公式，与经典 Bezout 定理一致。
7. Maslov 退量子化提供了从经典代数到热带代数的连续变形，参数 $\hbar \to 0$ 给出热带极限。
8. 阿米巴是代数簇的对数映射像，热带簇是阿米巴的极限骨架，Ronkin 函数描述其凸结构。
9. 热带 Grassmannian $\mathrm{TGr}(2,n)$ 与系统发生树空间同胚，为进化生物学提供了代数几何工具。
10. 热带代数曲线的组合结构（加权平衡图）使枚举几何中的曲线计数问题变得可处理（Mikhalkin 对应定理）。
11. 热带几何与镜像对称（SYZ 猜想、Gross-Siebert 纲领）有深层联系，是现代数学物理的活跃前沿。
