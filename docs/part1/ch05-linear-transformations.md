# 第 5 章 线性变换

线性变换（linear transformation）是线性代数中最核心的概念之一，它将向量空间的代数结构与几何直觉统一起来。在前几章中，我们已经看到矩阵乘法 $A\mathbf{x}$ 将一个向量映射为另一个向量；本章将从抽象的角度认识到，这种映射的本质——保持加法和标量乘法——才是最关键的性质。我们将系统地研究线性变换的定义与性质、核与像、秩-零化度定理、矩阵表示、基变换与相似矩阵、复合与逆、同构以及不变子空间等内容，为后续的特征值理论和内积空间理论打下坚实基础。

---

## 5.1 线性变换的定义

!!! definition "定义 5.1 (线性变换)"
    设 $V$ 和 $W$ 为数域 $\mathbb{F}$ 上的向量空间。映射 $T: V \to W$ 称为从 $V$ 到 $W$ 的一个**线性变换**（linear transformation），若对所有 $\mathbf{u}, \mathbf{v} \in V$ 和所有标量 $c \in \mathbb{F}$，满足：

    1. **可加性**（additivity）：$T(\mathbf{u} + \mathbf{v}) = T(\mathbf{u}) + T(\mathbf{v})$
    2. **齐次性**（homogeneity）：$T(c\mathbf{v}) = cT(\mathbf{v})$

    当 $W = V$ 时，$T$ 也称为 $V$ 上的**线性算子**（linear operator）。

!!! theorem "定理 5.1 (线性变换的等价条件)"
    设 $T: V \to W$ 为映射，则以下条件等价：

    1. $T$ 是线性变换。
    2. 对所有 $\mathbf{u}, \mathbf{v} \in V$ 和 $c, d \in \mathbb{F}$，有 $T(c\mathbf{u} + d\mathbf{v}) = cT(\mathbf{u}) + dT(\mathbf{v})$。
    3. 对任意有限个向量 $\mathbf{v}_1, \ldots, \mathbf{v}_k \in V$ 和标量 $c_1, \ldots, c_k \in \mathbb{F}$，有

    $$T\left(\sum_{i=1}^{k} c_i \mathbf{v}_i\right) = \sum_{i=1}^{k} c_i T(\mathbf{v}_i)$$

??? proof "证明"
    $(1) \Rightarrow (2)$：由可加性，$T(c\mathbf{u} + d\mathbf{v}) = T(c\mathbf{u}) + T(d\mathbf{v})$；再由齐次性，$= cT(\mathbf{u}) + dT(\mathbf{v})$。

    $(2) \Rightarrow (3)$：对 $k$ 进行归纳。$k = 1$ 时，取 $d = 0$ 即得 $T(c_1 \mathbf{v}_1) = c_1 T(\mathbf{v}_1)$。假设 $k-1$ 时成立，则

    $$T\left(\sum_{i=1}^{k} c_i \mathbf{v}_i\right) = T\left(\sum_{i=1}^{k-1} c_i \mathbf{v}_i + c_k \mathbf{v}_k\right) = T\left(\sum_{i=1}^{k-1} c_i \mathbf{v}_i\right) + c_k T(\mathbf{v}_k) = \sum_{i=1}^{k} c_i T(\mathbf{v}_i)$$

    $(3) \Rightarrow (1)$：取 $k = 2$，$c_1 = c_2 = 1$ 得可加性；取 $k = 1$ 得齐次性。 $\blacksquare$

!!! proposition "命题 5.1 (线性变换的基本性质)"
    设 $T: V \to W$ 为线性变换，则：

    1. $T(\mathbf{0}_V) = \mathbf{0}_W$
    2. $T(-\mathbf{v}) = -T(\mathbf{v})$
    3. $T(\mathbf{u} - \mathbf{v}) = T(\mathbf{u}) - T(\mathbf{v})$

??? proof "证明"
    1. $T(\mathbf{0}) = T(0 \cdot \mathbf{v}) = 0 \cdot T(\mathbf{v}) = \mathbf{0}$。
    2. $T(-\mathbf{v}) = T((-1)\mathbf{v}) = (-1)T(\mathbf{v}) = -T(\mathbf{v})$。
    3. $T(\mathbf{u} - \mathbf{v}) = T(\mathbf{u} + (-\mathbf{v})) = T(\mathbf{u}) + T(-\mathbf{v}) = T(\mathbf{u}) - T(\mathbf{v})$。 $\blacksquare$

!!! example "例 5.1"
    **旋转变换**：$T: \mathbb{R}^2 \to \mathbb{R}^2$ 将平面上的向量逆时针旋转角度 $\theta$：

    $$T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} x\cos\theta - y\sin\theta \\ x\sin\theta + y\cos\theta \end{pmatrix}$$

    可以验证 $T$ 满足可加性和齐次性，因此是线性变换。其几何意义是将整个平面刚性地旋转 $\theta$ 角。

!!! example "例 5.2"
    **正交投影**：$T: \mathbb{R}^3 \to \mathbb{R}^3$ 将空间中的向量投影到 $xy$ 平面：

    $$T\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} x \\ y \\ 0 \end{pmatrix}$$

    验证：$T(c\mathbf{u} + d\mathbf{v}) = cT(\mathbf{u}) + dT(\mathbf{v})$，因此 $T$ 是线性变换。

!!! example "例 5.3"
    **微分算子**：设 $V = C^1(\mathbb{R})$ 为所有连续可微函数的空间，$W = C(\mathbb{R})$ 为所有连续函数的空间。定义 $D: V \to W$，$D(f) = f'$。由导数的线性性：

    $$D(cf + dg) = (cf + dg)' = cf' + dg' = cD(f) + dD(g)$$

    因此微分算子 $D$ 是线性变换。

!!! example "例 5.4"
    **非线性映射的反例**：映射 $T: \mathbb{R}^2 \to \mathbb{R}^2$，$T(\mathbf{v}) = \mathbf{v} + \begin{pmatrix} 1 \\ 0 \end{pmatrix}$（平移）不是线性变换，因为 $T(\mathbf{0}) = \begin{pmatrix} 1 \\ 0 \end{pmatrix} \neq \mathbf{0}$。

!!! theorem "定理 5.2 (线性变换由基上的值唯一确定)"
    设 $V$ 为有限维向量空间，$\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 是 $V$ 的一组基，$W$ 为任意向量空间。对任意给定的 $\mathbf{w}_1, \ldots, \mathbf{w}_n \in W$，存在唯一的线性变换 $T: V \to W$ 使得

    $$T(\mathbf{v}_i) = \mathbf{w}_i, \quad i = 1, 2, \ldots, n$$

??? proof "证明"
    **存在性**：任意 $\mathbf{v} \in V$ 可唯一写成 $\mathbf{v} = c_1 \mathbf{v}_1 + \cdots + c_n \mathbf{v}_n$。定义

    $$T(\mathbf{v}) = c_1 \mathbf{w}_1 + \cdots + c_n \mathbf{w}_n$$

    验证线性性：设 $\mathbf{u} = a_1 \mathbf{v}_1 + \cdots + a_n \mathbf{v}_n$，$\mathbf{v} = b_1 \mathbf{v}_1 + \cdots + b_n \mathbf{v}_n$，则

    $$T(\alpha\mathbf{u} + \beta\mathbf{v}) = T\left(\sum_i (\alpha a_i + \beta b_i)\mathbf{v}_i\right) = \sum_i (\alpha a_i + \beta b_i)\mathbf{w}_i = \alpha \sum_i a_i \mathbf{w}_i + \beta \sum_i b_i \mathbf{w}_i = \alpha T(\mathbf{u}) + \beta T(\mathbf{v})$$

    **唯一性**：若 $S: V \to W$ 也满足 $S(\mathbf{v}_i) = \mathbf{w}_i$，则对任意 $\mathbf{v} = \sum c_i \mathbf{v}_i$，

    $$S(\mathbf{v}) = \sum c_i S(\mathbf{v}_i) = \sum c_i \mathbf{w}_i = T(\mathbf{v})$$

    故 $S = T$。 $\blacksquare$

---

## 5.2 核与像

!!! definition "定义 5.2 (核)"
    设 $T: V \to W$ 为线性变换。$T$ 的**核**（kernel），也称**零空间**（null space），定义为

    $$\ker(T) = \{\mathbf{v} \in V : T(\mathbf{v}) = \mathbf{0}\}$$

!!! definition "定义 5.3 (像)"
    设 $T: V \to W$ 为线性变换。$T$ 的**像**（image），也称**值域**（range），定义为

    $$\operatorname{im}(T) = \{T(\mathbf{v}) : \mathbf{v} \in V\} = \{w \in W : \exists \mathbf{v} \in V,\; T(\mathbf{v}) = \mathbf{w}\}$$

!!! theorem "定理 5.3 (核与像是子空间)"
    设 $T: V \to W$ 为线性变换，则：

    1. $\ker(T)$ 是 $V$ 的子空间。
    2. $\operatorname{im}(T)$ 是 $W$ 的子空间。

??? proof "证明"
    **(1)** 验证 $\ker(T)$ 是子空间：

    - 非空：由命题 5.1，$T(\mathbf{0}) = \mathbf{0}$，故 $\mathbf{0} \in \ker(T)$。
    - 封闭性：设 $\mathbf{u}, \mathbf{v} \in \ker(T)$，$c \in \mathbb{F}$。则 $T(\mathbf{u} + \mathbf{v}) = T(\mathbf{u}) + T(\mathbf{v}) = \mathbf{0} + \mathbf{0} = \mathbf{0}$，故 $\mathbf{u} + \mathbf{v} \in \ker(T)$。$T(c\mathbf{u}) = cT(\mathbf{u}) = c\mathbf{0} = \mathbf{0}$，故 $c\mathbf{u} \in \ker(T)$。

    **(2)** 验证 $\operatorname{im}(T)$ 是子空间：

    - 非空：$T(\mathbf{0}) = \mathbf{0} \in \operatorname{im}(T)$。
    - 封闭性：设 $\mathbf{w}_1 = T(\mathbf{v}_1), \mathbf{w}_2 = T(\mathbf{v}_2) \in \operatorname{im}(T)$，$c \in \mathbb{F}$。则 $\mathbf{w}_1 + \mathbf{w}_2 = T(\mathbf{v}_1) + T(\mathbf{v}_2) = T(\mathbf{v}_1 + \mathbf{v}_2) \in \operatorname{im}(T)$。$c\mathbf{w}_1 = cT(\mathbf{v}_1) = T(c\mathbf{v}_1) \in \operatorname{im}(T)$。 $\blacksquare$

!!! definition "定义 5.4 (零化度与秩)"
    设 $T: V \to W$ 为有限维向量空间之间的线性变换。

    - $T$ 的**零化度**（nullity）定义为 $\operatorname{nullity}(T) = \dim(\ker(T))$
    - $T$ 的**秩**（rank）定义为 $\operatorname{rank}(T) = \dim(\operatorname{im}(T))$

!!! theorem "定理 5.4 (单射与满射的刻画)"
    设 $T: V \to W$ 为线性变换，则：

    1. $T$ 是**单射**（injective）当且仅当 $\ker(T) = \{\mathbf{0}\}$。
    2. $T$ 是**满射**（surjective）当且仅当 $\operatorname{im}(T) = W$。

??? proof "证明"
    **(1)** $(\Rightarrow)$ 若 $T$ 单射，设 $\mathbf{v} \in \ker(T)$，则 $T(\mathbf{v}) = \mathbf{0} = T(\mathbf{0})$，由单射性 $\mathbf{v} = \mathbf{0}$。

    $(\Leftarrow)$ 若 $\ker(T) = \{\mathbf{0}\}$，设 $T(\mathbf{u}) = T(\mathbf{v})$，则 $T(\mathbf{u} - \mathbf{v}) = \mathbf{0}$，即 $\mathbf{u} - \mathbf{v} \in \ker(T) = \{\mathbf{0}\}$，故 $\mathbf{u} = \mathbf{v}$。

    **(2)** 由像的定义直接可得。 $\blacksquare$

!!! example "例 5.5"
    设 $T: \mathbb{R}^3 \to \mathbb{R}^2$ 定义为 $T\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} x - y \\ 2y + z \end{pmatrix}$。

    求 $\ker(T)$：解 $T(\mathbf{v}) = \mathbf{0}$，即 $x - y = 0$，$2y + z = 0$。设 $y = t$，则 $x = t$，$z = -2t$，故

    $$\ker(T) = \left\{ t\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix} : t \in \mathbb{R} \right\}$$

    $\operatorname{nullity}(T) = 1$。

    求 $\operatorname{im}(T)$：$T(\mathbf{v}) = x\begin{pmatrix} 1 \\ 0 \end{pmatrix} + y\begin{pmatrix} -1 \\ 2 \end{pmatrix} + z\begin{pmatrix} 0 \\ 1 \end{pmatrix}$。由于 $\begin{pmatrix} 1 \\ 0 \end{pmatrix}$ 和 $\begin{pmatrix} -1 \\ 2 \end{pmatrix}$ 线性无关，故 $\operatorname{im}(T) = \mathbb{R}^2$，$\operatorname{rank}(T) = 2$。

    验证：$\operatorname{nullity}(T) + \operatorname{rank}(T) = 1 + 2 = 3 = \dim(\mathbb{R}^3)$。

---

## 5.3 秩-零化度定理

!!! theorem "定理 5.5 (秩-零化度定理 / Rank-Nullity Theorem)"
    设 $V$ 为有限维向量空间，$T: V \to W$ 为线性变换，则

    $$\dim(V) = \operatorname{rank}(T) + \operatorname{nullity}(T) = \dim(\operatorname{im}(T)) + \dim(\ker(T))$$

??? proof "证明"
    设 $\dim(V) = n$，$\dim(\ker(T)) = r$。设 $\{\mathbf{u}_1, \ldots, \mathbf{u}_r\}$ 是 $\ker(T)$ 的一组基。将其扩充为 $V$ 的基 $\{\mathbf{u}_1, \ldots, \mathbf{u}_r, \mathbf{v}_1, \ldots, \mathbf{v}_s\}$，其中 $r + s = n$。

    我们证明 $\{T(\mathbf{v}_1), \ldots, T(\mathbf{v}_s)\}$ 是 $\operatorname{im}(T)$ 的基。

    **生成性**：任意 $\mathbf{w} \in \operatorname{im}(T)$，存在 $\mathbf{v} = a_1 \mathbf{u}_1 + \cdots + a_r \mathbf{u}_r + b_1 \mathbf{v}_1 + \cdots + b_s \mathbf{v}_s$ 使得 $T(\mathbf{v}) = \mathbf{w}$。由于 $T(\mathbf{u}_i) = \mathbf{0}$，故

    $$\mathbf{w} = T(\mathbf{v}) = b_1 T(\mathbf{v}_1) + \cdots + b_s T(\mathbf{v}_s)$$

    **线性无关性**：设 $c_1 T(\mathbf{v}_1) + \cdots + c_s T(\mathbf{v}_s) = \mathbf{0}$，则 $T(c_1 \mathbf{v}_1 + \cdots + c_s \mathbf{v}_s) = \mathbf{0}$，即 $c_1 \mathbf{v}_1 + \cdots + c_s \mathbf{v}_s \in \ker(T)$。因此存在标量 $d_1, \ldots, d_r$ 使得

    $$c_1 \mathbf{v}_1 + \cdots + c_s \mathbf{v}_s = d_1 \mathbf{u}_1 + \cdots + d_r \mathbf{u}_r$$

    即 $d_1 \mathbf{u}_1 + \cdots + d_r \mathbf{u}_r - c_1 \mathbf{v}_1 - \cdots - c_s \mathbf{v}_s = \mathbf{0}$。由于 $\{\mathbf{u}_1, \ldots, \mathbf{u}_r, \mathbf{v}_1, \ldots, \mathbf{v}_s\}$ 是 $V$ 的基（线性无关），故所有系数为零，特别地 $c_1 = \cdots = c_s = 0$。

    因此 $\dim(\operatorname{im}(T)) = s = n - r = \dim(V) - \dim(\ker(T))$。 $\blacksquare$

!!! corollary "推论 5.1"
    设 $T: V \to W$ 为线性变换，$\dim(V) = n$，$\dim(W) = m$，则：

    1. $\operatorname{rank}(T) \leq \min(n, m)$
    2. 若 $n = m$，则 $T$ 单射 $\Leftrightarrow$ $T$ 满射 $\Leftrightarrow$ $T$ 双射

??? proof "证明"
    1. $\operatorname{rank}(T) = \dim(\operatorname{im}(T)) \leq \dim(W) = m$，且 $\operatorname{rank}(T) = n - \operatorname{nullity}(T) \leq n$。

    2. 当 $n = m$ 时，$T$ 单射 $\Leftrightarrow$ $\operatorname{nullity}(T) = 0$ $\Leftrightarrow$ $\operatorname{rank}(T) = n = m$ $\Leftrightarrow$ $\operatorname{im}(T) = W$ $\Leftrightarrow$ $T$ 满射。 $\blacksquare$

!!! example "例 5.6"
    设 $T: \mathbb{P}_3 \to \mathbb{P}_2$ 为微分算子 $T(p) = p'$。则 $\ker(T)$ 是常数多项式的集合（即 $\mathbb{P}_0$），$\dim(\ker(T)) = 1$。由秩-零化度定理：

    $$\operatorname{rank}(T) = \dim(\mathbb{P}_3) - \operatorname{nullity}(T) = 4 - 1 = 3 = \dim(\mathbb{P}_2)$$

    故 $T$ 是满射。

---

## 5.4 线性变换的矩阵表示

!!! definition "定义 5.5 (线性变换的矩阵表示)"
    设 $V$ 为 $n$ 维向量空间，$W$ 为 $m$ 维向量空间，$\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 为 $V$ 的基，$\mathcal{C} = \{\mathbf{w}_1, \ldots, \mathbf{w}_m\}$ 为 $W$ 的基。设 $T: V \to W$ 为线性变换。

    对每个 $j = 1, \ldots, n$，将 $T(\mathbf{v}_j)$ 用基 $\mathcal{C}$ 表示：

    $$T(\mathbf{v}_j) = a_{1j}\mathbf{w}_1 + a_{2j}\mathbf{w}_2 + \cdots + a_{mj}\mathbf{w}_m$$

    则 $m \times n$ 矩阵

    $$[T]_{\mathcal{B}}^{\mathcal{C}} = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}$$

    称为 $T$ 关于基 $\mathcal{B}$ 和 $\mathcal{C}$ 的**矩阵表示**（matrix representation）。

!!! theorem "定理 5.6 (矩阵表示与坐标向量)"
    设符号同定义 5.5，若 $\mathbf{v} \in V$ 在基 $\mathcal{B}$ 下的坐标向量为 $[\mathbf{v}]_{\mathcal{B}} \in \mathbb{F}^n$，则

    $$[T(\mathbf{v})]_{\mathcal{C}} = [T]_{\mathcal{B}}^{\mathcal{C}} \cdot [\mathbf{v}]_{\mathcal{B}}$$

    即线性变换在坐标向量上的作用等价于矩阵乘法。

??? proof "证明"
    设 $\mathbf{v} = c_1 \mathbf{v}_1 + \cdots + c_n \mathbf{v}_n$，则 $[\mathbf{v}]_{\mathcal{B}} = (c_1, \ldots, c_n)^T$。

    $$T(\mathbf{v}) = \sum_{j=1}^n c_j T(\mathbf{v}_j) = \sum_{j=1}^n c_j \left(\sum_{i=1}^m a_{ij} \mathbf{w}_i\right) = \sum_{i=1}^m \left(\sum_{j=1}^n a_{ij} c_j\right) \mathbf{w}_i$$

    因此 $[T(\mathbf{v})]_{\mathcal{C}}$ 的第 $i$ 个分量为 $\sum_{j=1}^n a_{ij} c_j$，这正是矩阵 $[T]_{\mathcal{B}}^{\mathcal{C}}$ 与向量 $[\mathbf{v}]_{\mathcal{B}}$ 乘积的第 $i$ 个分量。 $\blacksquare$

!!! example "例 5.7"
    设 $T: \mathbb{P}_2 \to \mathbb{P}_1$ 为微分算子 $T(p) = p'$。取 $\mathbb{P}_2$ 的基 $\mathcal{B} = \{1, t, t^2\}$，$\mathbb{P}_1$ 的基 $\mathcal{C} = \{1, t\}$。

    $$T(1) = 0 = 0 \cdot 1 + 0 \cdot t, \quad T(t) = 1 = 1 \cdot 1 + 0 \cdot t, \quad T(t^2) = 2t = 0 \cdot 1 + 2 \cdot t$$

    故

    $$[T]_{\mathcal{B}}^{\mathcal{C}} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$$

    验证：$p(t) = 3 + 2t - t^2$，$[p]_{\mathcal{B}} = (3, 2, -1)^T$。

    $$[T]_{\mathcal{B}}^{\mathcal{C}} [p]_{\mathcal{B}} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} 3 \\ 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 2 \\ -2 \end{pmatrix}$$

    即 $T(p) = p' = 2 - 2t$，其在基 $\mathcal{C}$ 下的坐标为 $(2, -2)^T$，正确。

---

## 5.5 基变换与相似矩阵

!!! definition "定义 5.6 (过渡矩阵)"
    设 $V$ 为 $n$ 维向量空间，$\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 和 $\mathcal{B}' = \{\mathbf{v}_1', \ldots, \mathbf{v}_n'\}$ 是 $V$ 的两组基。从 $\mathcal{B}$ 到 $\mathcal{B}'$ 的**过渡矩阵**（transition matrix / change of basis matrix）$P$ 是 $n \times n$ 矩阵，其第 $j$ 列为 $\mathbf{v}_j$ 在基 $\mathcal{B}'$ 下的坐标向量：

    $$P = \Big([\mathbf{v}_1]_{\mathcal{B}'} \;\; [\mathbf{v}_2]_{\mathcal{B}'} \;\; \cdots \;\; [\mathbf{v}_n]_{\mathcal{B}'}\Big)$$

    则对任意 $\mathbf{v} \in V$，有 $[\mathbf{v}]_{\mathcal{B}'} = P [\mathbf{v}]_{\mathcal{B}}$。（注：部分教材定义方向相反，此处约定需注意。）

!!! definition "定义 5.7 (相似矩阵)"
    设 $A, B$ 为 $n$ 阶方阵。若存在可逆矩阵 $P$ 使得

    $$B = P^{-1}AP$$

    则称 $A$ 与 $B$ **相似**（similar），记为 $A \sim B$。

!!! theorem "定理 5.7 (基变换公式)"
    设 $T: V \to V$ 为线性算子，$\mathcal{B}$ 和 $\mathcal{B}'$ 是 $V$ 的两组基，$P$ 为从 $\mathcal{B}$ 到 $\mathcal{B}'$ 的过渡矩阵，则

    $$[T]_{\mathcal{B}'} = P^{-1} [T]_{\mathcal{B}} P$$

    即 $T$ 在不同基下的矩阵表示是相似的。

??? proof "证明"
    设 $A = [T]_{\mathcal{B}}$，$B = [T]_{\mathcal{B}'}$。对任意 $\mathbf{v} \in V$，

    $$[T(\mathbf{v})]_{\mathcal{B}} = A [\mathbf{v}]_{\mathcal{B}}$$

    又 $[\mathbf{v}]_{\mathcal{B}'} = P[\mathbf{v}]_{\mathcal{B}}$，即 $[\mathbf{v}]_{\mathcal{B}} = P^{-1}[\mathbf{v}]_{\mathcal{B}'}$，因此

    $$[T(\mathbf{v})]_{\mathcal{B}'} = P [T(\mathbf{v})]_{\mathcal{B}} = P A [\mathbf{v}]_{\mathcal{B}} = P A P^{-1} [\mathbf{v}]_{\mathcal{B}'}$$

    由唯一性，$B = PAP^{-1}$。

    **注意**：此处过渡矩阵的定义决定了公式形式。若采用 $[\mathbf{v}]_{\mathcal{B}} = P[\mathbf{v}]_{\mathcal{B}'}$ 的约定，则 $B = P^{-1}AP$。不同教材的约定不同，本质相同。 $\blacksquare$

!!! theorem "定理 5.8 (相似是等价关系)"
    矩阵的相似关系满足：

    1. **自反性**：$A \sim A$（取 $P = I$）
    2. **对称性**：若 $A \sim B$，则 $B \sim A$
    3. **传递性**：若 $A \sim B$ 且 $B \sim C$，则 $A \sim C$

??? proof "证明"
    1. $A = I^{-1}AI$。
    2. 若 $B = P^{-1}AP$，则 $A = PBP^{-1} = (P^{-1})^{-1}B(P^{-1})$，取 $Q = P^{-1}$ 即可。
    3. 若 $B = P^{-1}AP$，$C = Q^{-1}BQ$，则 $C = Q^{-1}P^{-1}APQ = (PQ)^{-1}A(PQ)$。 $\blacksquare$

!!! example "例 5.8"
    设 $T: \mathbb{R}^2 \to \mathbb{R}^2$，$T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2x + y \\ 3y \end{pmatrix}$。

    取标准基 $\mathcal{B} = \{\mathbf{e}_1, \mathbf{e}_2\}$，则 $[T]_{\mathcal{B}} = A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}$。

    取新基 $\mathcal{B}' = \left\{\begin{pmatrix} 1 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 \\ 1 \end{pmatrix}\right\}$。过渡矩阵 $P = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$，$P^{-1} = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$。

    $$[T]_{\mathcal{B}'} = P^{-1}AP = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$$

    在新基下，$T$ 的矩阵变为对角矩阵，说明 $\mathcal{B}'$ 是由 $T$ 的特征向量构成的基。

---

## 5.6 线性变换的复合与逆

!!! theorem "定理 5.9 (线性变换的复合仍为线性变换)"
    设 $T: U \to V$ 和 $S: V \to W$ 为线性变换，则复合映射 $S \circ T: U \to W$ 也是线性变换，且

    $$[S \circ T]_{\mathcal{B}}^{\mathcal{D}} = [S]_{\mathcal{C}}^{\mathcal{D}} \cdot [T]_{\mathcal{B}}^{\mathcal{C}}$$

    其中 $\mathcal{B}, \mathcal{C}, \mathcal{D}$ 分别为 $U, V, W$ 的基。

??? proof "证明"
    **线性性**：对任意 $\mathbf{u}_1, \mathbf{u}_2 \in U$ 和 $c \in \mathbb{F}$，

    $$(S \circ T)(c\mathbf{u}_1 + \mathbf{u}_2) = S(T(c\mathbf{u}_1 + \mathbf{u}_2)) = S(cT(\mathbf{u}_1) + T(\mathbf{u}_2)) = cS(T(\mathbf{u}_1)) + S(T(\mathbf{u}_2)) = c(S \circ T)(\mathbf{u}_1) + (S \circ T)(\mathbf{u}_2)$$

    **矩阵等式**：对任意 $\mathbf{u} \in U$，

    $$[(S \circ T)(\mathbf{u})]_{\mathcal{D}} = [S(T(\mathbf{u}))]_{\mathcal{D}} = [S]_{\mathcal{C}}^{\mathcal{D}} [T(\mathbf{u})]_{\mathcal{C}} = [S]_{\mathcal{C}}^{\mathcal{D}} [T]_{\mathcal{B}}^{\mathcal{C}} [\mathbf{u}]_{\mathcal{B}}$$

    由唯一性得 $[S \circ T]_{\mathcal{B}}^{\mathcal{D}} = [S]_{\mathcal{C}}^{\mathcal{D}} [T]_{\mathcal{B}}^{\mathcal{C}}$。 $\blacksquare$

!!! definition "定义 5.8 (可逆线性变换)"
    线性变换 $T: V \to W$ 称为**可逆的**（invertible），若存在线性变换 $T^{-1}: W \to V$ 使得

    $$T^{-1} \circ T = I_V, \quad T \circ T^{-1} = I_W$$

    其中 $I_V, I_W$ 分别为 $V, W$ 上的恒等变换。

!!! theorem "定理 5.10 (可逆性的等价条件)"
    设 $T: V \to W$ 为线性变换，$\dim(V) = \dim(W) = n$，则以下条件等价：

    1. $T$ 可逆。
    2. $T$ 是双射。
    3. $\ker(T) = \{\mathbf{0}\}$。
    4. $\operatorname{im}(T) = W$。
    5. $T$ 将 $V$ 的基映射为 $W$ 的基。
    6. $[T]_{\mathcal{B}}^{\mathcal{C}}$ 是可逆矩阵（对任意基 $\mathcal{B}, \mathcal{C}$）。

??? proof "证明"
    $(1) \Leftrightarrow (2)$：可逆当且仅当双射，这是映射论的基本事实。

    $(2) \Leftrightarrow (3) \Leftrightarrow (4)$：由定理 5.4 和推论 5.1，当 $\dim(V) = \dim(W)$ 时，单射、满射和双射等价。

    $(2) \Rightarrow (5)$：设 $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 为 $V$ 的基。由满射，$\{T(\mathbf{v}_1), \ldots, T(\mathbf{v}_n)\}$ 生成 $W$。由单射，若 $\sum c_i T(\mathbf{v}_i) = \mathbf{0}$，则 $T(\sum c_i \mathbf{v}_i) = \mathbf{0}$，故 $\sum c_i \mathbf{v}_i = \mathbf{0}$，因此 $c_i = 0$。

    $(5) \Rightarrow (4)$：若 $T$ 将基映射为基，则 $\operatorname{im}(T) \supseteq \operatorname{span}\{T(\mathbf{v}_1), \ldots, T(\mathbf{v}_n)\} = W$。

    $(3) \Leftrightarrow (6)$：$\ker(T) = \{\mathbf{0}\}$ 当且仅当 $[T]_{\mathcal{B}}^{\mathcal{C}} \mathbf{x} = \mathbf{0}$ 仅有零解，即矩阵可逆。 $\blacksquare$

!!! example "例 5.9"
    设 $T: \mathbb{R}^2 \to \mathbb{R}^2$，$T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2x - y \\ x + 3y \end{pmatrix}$。标准矩阵为 $A = \begin{pmatrix} 2 & -1 \\ 1 & 3 \end{pmatrix}$。

    $\det(A) = 6 - (-1) = 7 \neq 0$，故 $T$ 可逆。

    $$A^{-1} = \frac{1}{7}\begin{pmatrix} 3 & 1 \\ -1 & 2 \end{pmatrix}$$

    因此 $T^{-1}\begin{pmatrix} x \\ y \end{pmatrix} = \frac{1}{7}\begin{pmatrix} 3x + y \\ -x + 2y \end{pmatrix}$。

---

## 5.7 同构

!!! definition "定义 5.9 (同构)"
    若存在可逆线性变换 $T: V \to W$，则称 $V$ 与 $W$ **同构**（isomorphic），记为 $V \cong W$。此时 $T$ 称为从 $V$ 到 $W$ 的一个**同构映射**（isomorphism）。

!!! theorem "定理 5.11 (有限维向量空间的分类定理)"
    设 $V, W$ 为数域 $\mathbb{F}$ 上的有限维向量空间，则

    $$V \cong W \quad \Longleftrightarrow \quad \dim(V) = \dim(W)$$

    特别地，每个 $n$ 维向量空间都同构于 $\mathbb{F}^n$。

??? proof "证明"
    $(\Rightarrow)$ 设 $T: V \to W$ 为同构映射，$\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 为 $V$ 的基。由定理 5.10 的条件 (5)，$\{T(\mathbf{v}_1), \ldots, T(\mathbf{v}_n)\}$ 是 $W$ 的基，故 $\dim(W) = n = \dim(V)$。

    $(\Leftarrow)$ 设 $\dim(V) = \dim(W) = n$，取 $V$ 的基 $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$，$W$ 的基 $\mathcal{C} = \{\mathbf{w}_1, \ldots, \mathbf{w}_n\}$。由定理 5.2，存在唯一的线性变换 $T$ 使得 $T(\mathbf{v}_i) = \mathbf{w}_i$。由于 $T$ 将基映射为基，由定理 5.10，$T$ 可逆。

    特别地，$V$ 的坐标映射 $\phi_{\mathcal{B}}: V \to \mathbb{F}^n$，$\phi_{\mathcal{B}}(\mathbf{v}) = [\mathbf{v}]_{\mathcal{B}}$ 是同构映射。 $\blacksquare$

!!! note "注"
    同构定理说明：在有限维的框架下，向量空间的"本质"完全由维数决定。所有 $n$ 维向量空间——无论是 $\mathbb{R}^n$、$\mathbb{P}_{n-1}$、$\mathbb{R}^{m \times k}$（当 $mk = n$ 时）——从线性代数的角度看都是"相同的"。这是线性代数理论的强大之处，也解释了为何我们总可以用矩阵和向量来处理抽象问题。

!!! example "例 5.10"
    $\mathbb{R}^{2 \times 2}$（$2 \times 2$ 实矩阵的空间）与 $\mathbb{R}^4$ 同构。一个具体的同构映射为

    $$T\begin{pmatrix} a & b \\ c & d \end{pmatrix} = \begin{pmatrix} a \\ b \\ c \\ d \end{pmatrix}$$

    这个映射就是将矩阵"拉直"为列向量（向量化，vectorization）。

!!! theorem "定理 5.12 (同构保持线性结构)"
    设 $T: V \to W$ 为同构映射，则：

    1. $T$ 将线性无关集映射为线性无关集。
    2. $T$ 将生成集映射为生成集。
    3. $T$ 将 $V$ 的基映射为 $W$ 的基。
    4. 对任意子空间 $U \leq V$，$T(U)$ 是 $W$ 的子空间，且 $\dim(T(U)) = \dim(U)$。

??? proof "证明"
    1. 设 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ 线性无关。若 $\sum c_i T(\mathbf{v}_i) = \mathbf{0}$，则 $T(\sum c_i \mathbf{v}_i) = \mathbf{0}$。由 $T$ 单射，$\sum c_i \mathbf{v}_i = \mathbf{0}$，故 $c_i = 0$。

    2. 设 $\operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\} = V$。对任意 $\mathbf{w} \in W$，由 $T$ 满射，存在 $\mathbf{v} = \sum c_i \mathbf{v}_i \in V$ 使得 $T(\mathbf{v}) = \mathbf{w}$，即 $\mathbf{w} = \sum c_i T(\mathbf{v}_i)$。

    3. 由 1 和 2 立得。

    4. $T(U)$ 的子空间性由 $T$ 的线性性保证。$T$ 限制在 $U$ 上仍为单射，故 $\ker(T|_U) = \{\mathbf{0}\}$，由秩-零化度定理，$\dim(T(U)) = \dim(U)$。 $\blacksquare$

---

## 5.8 不变子空间

!!! definition "定义 5.10 (不变子空间)"
    设 $T: V \to V$ 为线性算子，$U$ 为 $V$ 的子空间。若 $T(U) \subseteq U$（即对任意 $\mathbf{u} \in U$，有 $T(\mathbf{u}) \in U$），则称 $U$ 为 $T$ 的**不变子空间**（invariant subspace）。

!!! proposition "命题 5.2 (平凡不变子空间)"
    设 $T: V \to V$ 为线性算子，则：

    1. $\{\mathbf{0}\}$ 和 $V$ 是 $T$ 的不变子空间。
    2. $\ker(T)$ 是 $T$ 的不变子空间。
    3. $\operatorname{im}(T)$ 是 $T$ 的不变子空间（注：$\operatorname{im}(T) \subseteq V$ 因为 $T$ 是算子）。

??? proof "证明"
    1. 显然 $T(\{\mathbf{0}\}) = \{\mathbf{0}\} \subseteq \{\mathbf{0}\}$，$T(V) \subseteq V$。
    2. 若 $\mathbf{u} \in \ker(T)$，则 $T(\mathbf{u}) = \mathbf{0} \in \ker(T)$。
    3. 若 $\mathbf{w} \in \operatorname{im}(T)$，则 $\mathbf{w} = T(\mathbf{v})$ 对某 $\mathbf{v} \in V$，而 $T(\mathbf{w}) = T(T(\mathbf{v})) \in \operatorname{im}(T)$。 $\blacksquare$

!!! theorem "定理 5.13 (不变子空间与分块对角矩阵)"
    设 $T: V \to V$ 为线性算子，$\dim(V) = n$。若 $V = U_1 \oplus U_2$（直和分解），$\dim(U_1) = k$，$\dim(U_2) = n - k$，且 $U_1, U_2$ 都是 $T$ 的不变子空间，则存在 $V$ 的基 $\mathcal{B}$ 使得

    $$[T]_{\mathcal{B}} = \begin{pmatrix} A_1 & 0 \\ 0 & A_2 \end{pmatrix}$$

    其中 $A_1$ 为 $k \times k$ 矩阵（$T|_{U_1}$ 的矩阵表示），$A_2$ 为 $(n-k) \times (n-k)$ 矩阵（$T|_{U_2}$ 的矩阵表示）。

??? proof "证明"
    取 $U_1$ 的基 $\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$ 和 $U_2$ 的基 $\{\mathbf{u}_{k+1}, \ldots, \mathbf{u}_n\}$，合在一起构成 $V$ 的基 $\mathcal{B}$。

    由于 $U_1$ 是不变子空间，对 $j \leq k$，$T(\mathbf{u}_j) \in U_1$，因此 $T(\mathbf{u}_j)$ 可表示为 $\mathbf{u}_1, \ldots, \mathbf{u}_k$ 的线性组合，对 $\mathbf{u}_{k+1}, \ldots, \mathbf{u}_n$ 的系数为零。

    同理，对 $j > k$，$T(\mathbf{u}_j) \in U_2$，对 $\mathbf{u}_1, \ldots, \mathbf{u}_k$ 的系数为零。

    因此矩阵具有分块对角形式。 $\blacksquare$

!!! theorem "定理 5.14 (特征值与不变子空间)"
    设 $T: V \to V$ 为线性算子，$\lambda$ 为标量，则以下等价：

    1. $\lambda$ 是 $T$ 的特征值（即存在非零 $\mathbf{v}$ 使得 $T(\mathbf{v}) = \lambda\mathbf{v}$）。
    2. $T$ 存在一维不变子空间 $U = \operatorname{span}\{\mathbf{v}\}$，其中 $T|_U = \lambda \cdot I_U$。

??? proof "证明"
    $(1) \Rightarrow (2)$：设 $T(\mathbf{v}) = \lambda\mathbf{v}$，$\mathbf{v} \neq \mathbf{0}$。令 $U = \operatorname{span}\{\mathbf{v}\}$，则对任意 $c\mathbf{v} \in U$，$T(c\mathbf{v}) = cT(\mathbf{v}) = c\lambda\mathbf{v} = \lambda(c\mathbf{v}) \in U$。

    $(2) \Rightarrow (1)$：设 $U = \operatorname{span}\{\mathbf{v}\}$ 为一维不变子空间。则 $T(\mathbf{v}) \in U$，即 $T(\mathbf{v}) = \lambda\mathbf{v}$ 对某标量 $\lambda$。由 $\mathbf{v} \neq \mathbf{0}$，$\lambda$ 是 $T$ 的特征值。 $\blacksquare$

!!! example "例 5.11"
    设 $V = \mathbb{R}^2$，$T$ 为关于 $x$ 轴的反射：$T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} x \\ -y \end{pmatrix}$。

    - $U_1 = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 0 \end{pmatrix}\right\}$（$x$ 轴）是 $T$ 的不变子空间，$T|_{U_1} = I$。
    - $U_2 = \operatorname{span}\left\{\begin{pmatrix} 0 \\ 1 \end{pmatrix}\right\}$（$y$ 轴）是 $T$ 的不变子空间，$T|_{U_2} = -I$。
    - $V = U_1 \oplus U_2$，在基 $\{\mathbf{e}_1, \mathbf{e}_2\}$ 下 $[T] = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$，为对角矩阵。

!!! example "例 5.12"
    设 $T: \mathbb{R}^2 \to \mathbb{R}^2$ 为逆时针旋转 $\pi/2$，$T\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} -y \\ x \end{pmatrix}$。

    $T$ 在 $\mathbb{R}^2$ 上没有一维不变子空间（因为旋转 $90°$ 不会把任何直线映射到自身），因此 $T$ 在 $\mathbb{R}$ 上没有特征值。但在 $\mathbb{C}$ 上，$T$ 的特征值为 $\pm i$。

    这说明不变子空间（和特征值）的存在性可能依赖于底层数域的选取。

---

## 本章小结

本章从抽象角度研究了向量空间之间保持线性结构的映射——线性变换。核心内容包括：

1. **线性变换的定义**：保持加法和标量乘法的映射，由基上的值唯一确定。
2. **核与像**：分别刻画了线性变换"损失"和"覆盖"了多少信息。
3. **秩-零化度定理**：$\dim(V) = \operatorname{rank}(T) + \operatorname{nullity}(T)$，是有限维线性代数最基本的等式之一。
4. **矩阵表示**：线性变换在给定基下等价于矩阵乘法，搭建了抽象理论与计算的桥梁。
5. **基变换与相似**：同一线性算子在不同基下的矩阵通过 $B = P^{-1}AP$ 联系。
6. **复合与逆**：复合对应矩阵乘法，可逆当且仅当双射。
7. **同构**：有限维向量空间由维数完全分类。
8. **不变子空间**：为矩阵的对角化和分块结构提供了理论基础，直接通向特征值理论。
