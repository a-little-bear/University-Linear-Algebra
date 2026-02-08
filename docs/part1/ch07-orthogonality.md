# 第 7 章 正交性与最小二乘

正交性（orthogonality）是线性代数中几何直觉与代数理论交汇的核心概念。在 $\mathbb{R}^n$ 中，内积赋予了向量"长度"和"角度"的概念，使得我们能够讨论正交、投影、最佳逼近等问题。本章将从内积与范数出发，系统地研究正交向量、正交补、正交投影、Gram-Schmidt 正交化过程、正交矩阵、最小二乘法以及 QR 分解。这些工具不仅具有深刻的理论意义，在数值线性代数、统计学、信号处理和数据科学中也有着极其广泛的应用。

---

## 7.1 $\mathbb{R}^n$ 中的内积与范数

!!! definition "定义 7.1 (内积)"
    $\mathbb{R}^n$ 中向量 $\mathbf{u} = (u_1, \ldots, u_n)^T$ 和 $\mathbf{v} = (v_1, \ldots, v_n)^T$ 的**内积**（inner product），也称**点积**（dot product），定义为

    $$\langle \mathbf{u}, \mathbf{v} \rangle = \mathbf{u} \cdot \mathbf{v} = \mathbf{u}^T\mathbf{v} = \sum_{i=1}^n u_i v_i$$

!!! definition "定义 7.2 (范数)"
    向量 $\mathbf{v} \in \mathbb{R}^n$ 的**范数**（norm）或**长度**（length）定义为

    $$\|\mathbf{v}\| = \sqrt{\langle \mathbf{v}, \mathbf{v} \rangle} = \sqrt{\mathbf{v}^T\mathbf{v}} = \sqrt{\sum_{i=1}^n v_i^2}$$

    若 $\|\mathbf{v}\| = 1$，则称 $\mathbf{v}$ 为**单位向量**（unit vector）。将非零向量 $\mathbf{v}$ 变为同方向单位向量的过程称为**单位化**（normalization）：$\hat{\mathbf{v}} = \frac{\mathbf{v}}{\|\mathbf{v}\|}$。

!!! definition "定义 7.3 (距离)"
    $\mathbb{R}^n$ 中两点 $\mathbf{u}$ 和 $\mathbf{v}$ 的**距离**（distance）定义为

    $$d(\mathbf{u}, \mathbf{v}) = \|\mathbf{u} - \mathbf{v}\|$$

!!! proposition "命题 7.1 (内积的基本性质)"
    对任意 $\mathbf{u}, \mathbf{v}, \mathbf{w} \in \mathbb{R}^n$ 和 $c \in \mathbb{R}$：

    1. $\langle \mathbf{u}, \mathbf{v} \rangle = \langle \mathbf{v}, \mathbf{u} \rangle$（对称性）
    2. $\langle \mathbf{u} + \mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{u}, \mathbf{w} \rangle + \langle \mathbf{v}, \mathbf{w} \rangle$（可加性）
    3. $\langle c\mathbf{u}, \mathbf{v} \rangle = c\langle \mathbf{u}, \mathbf{v} \rangle$（齐次性）
    4. $\langle \mathbf{v}, \mathbf{v} \rangle \geq 0$，等号成立当且仅当 $\mathbf{v} = \mathbf{0}$（正定性）

!!! theorem "定理 7.1 (Cauchy-Schwarz 不等式)"
    对任意 $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$，

    $$|\langle \mathbf{u}, \mathbf{v} \rangle| \leq \|\mathbf{u}\| \cdot \|\mathbf{v}\|$$

    等号成立当且仅当 $\mathbf{u}$ 与 $\mathbf{v}$ 线性相关（即一个是另一个的标量倍）。

??? proof "证明"
    若 $\mathbf{v} = \mathbf{0}$，两边均为 $0$，不等式成立。设 $\mathbf{v} \neq \mathbf{0}$。

    对任意 $t \in \mathbb{R}$，考虑

    $$0 \leq \|\mathbf{u} - t\mathbf{v}\|^2 = \langle \mathbf{u} - t\mathbf{v}, \mathbf{u} - t\mathbf{v} \rangle = \|\mathbf{u}\|^2 - 2t\langle \mathbf{u}, \mathbf{v} \rangle + t^2\|\mathbf{v}\|^2$$

    这是关于 $t$ 的二次函数，其最小值在 $t = \frac{\langle \mathbf{u}, \mathbf{v} \rangle}{\|\mathbf{v}\|^2}$ 处取到。代入：

    $$0 \leq \|\mathbf{u}\|^2 - \frac{\langle \mathbf{u}, \mathbf{v} \rangle^2}{\|\mathbf{v}\|^2}$$

    即 $\langle \mathbf{u}, \mathbf{v} \rangle^2 \leq \|\mathbf{u}\|^2 \|\mathbf{v}\|^2$，两边取平方根得 $|\langle \mathbf{u}, \mathbf{v} \rangle| \leq \|\mathbf{u}\| \cdot \|\mathbf{v}\|$。

    等号成立 $\Leftrightarrow$ $\|\mathbf{u} - t\mathbf{v}\|^2 = 0$ $\Leftrightarrow$ $\mathbf{u} = t\mathbf{v}$ $\Leftrightarrow$ $\mathbf{u}, \mathbf{v}$ 线性相关。 $\blacksquare$

!!! theorem "定理 7.2 (三角不等式)"
    对任意 $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$，

    $$\|\mathbf{u} + \mathbf{v}\| \leq \|\mathbf{u}\| + \|\mathbf{v}\|$$

??? proof "证明"
    $$\|\mathbf{u} + \mathbf{v}\|^2 = \langle \mathbf{u} + \mathbf{v}, \mathbf{u} + \mathbf{v} \rangle = \|\mathbf{u}\|^2 + 2\langle \mathbf{u}, \mathbf{v} \rangle + \|\mathbf{v}\|^2$$

    由 Cauchy-Schwarz 不等式：

    $$\leq \|\mathbf{u}\|^2 + 2\|\mathbf{u}\|\|\mathbf{v}\| + \|\mathbf{v}\|^2 = (\|\mathbf{u}\| + \|\mathbf{v}\|)^2$$

    两边取平方根（两边均非负）即得。 $\blacksquare$

!!! example "例 7.1"
    设 $\mathbf{u} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$，$\mathbf{v} = \begin{pmatrix} 4 \\ -1 \\ 2 \end{pmatrix}$。

    $\langle \mathbf{u}, \mathbf{v} \rangle = 1 \cdot 4 + 2 \cdot (-1) + 3 \cdot 2 = 4 - 2 + 6 = 8$。

    $\|\mathbf{u}\| = \sqrt{1 + 4 + 9} = \sqrt{14}$，$\|\mathbf{v}\| = \sqrt{16 + 1 + 4} = \sqrt{21}$。

    验证 Cauchy-Schwarz：$|8| = 8 \leq \sqrt{14} \cdot \sqrt{21} = \sqrt{294} \approx 17.15$。成立。

    两向量的夹角：$\cos\theta = \frac{\langle \mathbf{u}, \mathbf{v} \rangle}{\|\mathbf{u}\|\|\mathbf{v}\|} = \frac{8}{\sqrt{294}}$，$\theta = \arccos\frac{8}{\sqrt{294}} \approx 62.2°$。

---

## 7.2 正交向量与正交集

!!! definition "定义 7.4 (正交)"
    两个向量 $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$ 称为**正交的**（orthogonal），若 $\langle \mathbf{u}, \mathbf{v} \rangle = 0$，记为 $\mathbf{u} \perp \mathbf{v}$。

!!! definition "定义 7.5 (正交集与正交基)"
    向量集合 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ 称为**正交集**（orthogonal set），若其中任意两个不同向量正交：

    $$\langle \mathbf{v}_i, \mathbf{v}_j \rangle = 0, \quad \forall\, i \neq j$$

    若进一步每个向量都是单位向量（$\|\mathbf{v}_i\| = 1$），则称之为**标准正交集**（orthonormal set）。

    若正交集（或标准正交集）同时构成某子空间的基，则分别称为**正交基**（orthogonal basis）和**标准正交基**（orthonormal basis）。

!!! theorem "定理 7.3 (正交集线性无关)"
    若 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ 是由非零向量组成的正交集，则该集合线性无关。

??? proof "证明"
    设 $c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_k\mathbf{v}_k = \mathbf{0}$。对任意 $j$，两边与 $\mathbf{v}_j$ 做内积：

    $$0 = \langle \mathbf{0}, \mathbf{v}_j \rangle = \left\langle \sum_{i=1}^k c_i\mathbf{v}_i, \mathbf{v}_j \right\rangle = \sum_{i=1}^k c_i \langle \mathbf{v}_i, \mathbf{v}_j \rangle = c_j \langle \mathbf{v}_j, \mathbf{v}_j \rangle = c_j \|\mathbf{v}_j\|^2$$

    由 $\mathbf{v}_j \neq \mathbf{0}$，$\|\mathbf{v}_j\|^2 > 0$，故 $c_j = 0$。 $\blacksquare$

!!! theorem "定理 7.4 (正交基下的坐标公式)"
    设 $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 为 $V$ 的正交基，则对任意 $\mathbf{w} \in V$，

    $$\mathbf{w} = \sum_{i=1}^n \frac{\langle \mathbf{w}, \mathbf{v}_i \rangle}{\langle \mathbf{v}_i, \mathbf{v}_i \rangle}\mathbf{v}_i = \sum_{i=1}^n \frac{\langle \mathbf{w}, \mathbf{v}_i \rangle}{\|\mathbf{v}_i\|^2}\mathbf{v}_i$$

    若基为标准正交基，则公式简化为 $\mathbf{w} = \sum_{i=1}^n \langle \mathbf{w}, \mathbf{v}_i \rangle \mathbf{v}_i$。

??? proof "证明"
    设 $\mathbf{w} = c_1\mathbf{v}_1 + \cdots + c_n\mathbf{v}_n$。两边与 $\mathbf{v}_j$ 做内积：

    $$\langle \mathbf{w}, \mathbf{v}_j \rangle = c_j \langle \mathbf{v}_j, \mathbf{v}_j \rangle = c_j\|\mathbf{v}_j\|^2$$

    故 $c_j = \frac{\langle \mathbf{w}, \mathbf{v}_j \rangle}{\|\mathbf{v}_j\|^2}$。 $\blacksquare$

!!! example "例 7.2"
    $\mathbb{R}^3$ 中，$\mathbf{v}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$，$\mathbf{v}_2 = \begin{pmatrix} -1 \\ 1 \\ 0 \end{pmatrix}$，$\mathbf{v}_3 = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$。

    验证正交性：$\langle \mathbf{v}_1, \mathbf{v}_2 \rangle = -1 + 1 + 0 = 0$，$\langle \mathbf{v}_1, \mathbf{v}_3 \rangle = 0$，$\langle \mathbf{v}_2, \mathbf{v}_3 \rangle = 0$。

    故 $\{\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3\}$ 是 $\mathbb{R}^3$ 的正交基。

    将 $\mathbf{w} = \begin{pmatrix} 3 \\ 5 \\ 7 \end{pmatrix}$ 用此正交基表示：

    $$c_1 = \frac{\langle \mathbf{w}, \mathbf{v}_1 \rangle}{\|\mathbf{v}_1\|^2} = \frac{3 + 5}{2} = 4, \quad c_2 = \frac{\langle \mathbf{w}, \mathbf{v}_2 \rangle}{\|\mathbf{v}_2\|^2} = \frac{-3 + 5}{2} = 1, \quad c_3 = \frac{\langle \mathbf{w}, \mathbf{v}_3 \rangle}{\|\mathbf{v}_3\|^2} = \frac{7}{1} = 7$$

    验证：$4\mathbf{v}_1 + 1\mathbf{v}_2 + 7\mathbf{v}_3 = \begin{pmatrix} 4 - 1 \\ 4 + 1 \\ 7 \end{pmatrix} = \begin{pmatrix} 3 \\ 5 \\ 7 \end{pmatrix} = \mathbf{w}$。正确。

---

## 7.3 正交补

!!! definition "定义 7.6 (正交补)"
    设 $W$ 为 $\mathbb{R}^n$ 的子空间。$W$ 的**正交补**（orthogonal complement）定义为

    $$W^\perp = \{\mathbf{v} \in \mathbb{R}^n : \langle \mathbf{v}, \mathbf{w} \rangle = 0, \;\forall\, \mathbf{w} \in W\}$$

!!! theorem "定理 7.5 (正交补是子空间)"
    $W^\perp$ 是 $\mathbb{R}^n$ 的子空间。

??? proof "证明"
    $\mathbf{0} \in W^\perp$（因为 $\langle \mathbf{0}, \mathbf{w} \rangle = 0$ 对所有 $\mathbf{w}$）。

    设 $\mathbf{u}, \mathbf{v} \in W^\perp$，$c \in \mathbb{R}$。对任意 $\mathbf{w} \in W$：

    $\langle \mathbf{u} + \mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{u}, \mathbf{w} \rangle + \langle \mathbf{v}, \mathbf{w} \rangle = 0 + 0 = 0$，故 $\mathbf{u} + \mathbf{v} \in W^\perp$。

    $\langle c\mathbf{u}, \mathbf{w} \rangle = c\langle \mathbf{u}, \mathbf{w} \rangle = 0$，故 $c\mathbf{u} \in W^\perp$。 $\blacksquare$

!!! theorem "定理 7.6 (正交分解定理)"
    设 $W$ 为 $\mathbb{R}^n$ 的子空间，则

    $$\mathbb{R}^n = W \oplus W^\perp$$

    即每个向量 $\mathbf{v} \in \mathbb{R}^n$ 可唯一地分解为

    $$\mathbf{v} = \mathbf{w} + \mathbf{w}^\perp, \quad \mathbf{w} \in W, \;\; \mathbf{w}^\perp \in W^\perp$$

??? proof "证明"
    **存在性**：取 $W$ 的一组标准正交基 $\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$。定义

    $$\mathbf{w} = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{u}_i \rangle \mathbf{u}_i \in W, \quad \mathbf{w}^\perp = \mathbf{v} - \mathbf{w}$$

    需验证 $\mathbf{w}^\perp \in W^\perp$：对任意 $j$，

    $$\langle \mathbf{w}^\perp, \mathbf{u}_j \rangle = \langle \mathbf{v} - \mathbf{w}, \mathbf{u}_j \rangle = \langle \mathbf{v}, \mathbf{u}_j \rangle - \sum_{i=1}^k \langle \mathbf{v}, \mathbf{u}_i \rangle \langle \mathbf{u}_i, \mathbf{u}_j \rangle = \langle \mathbf{v}, \mathbf{u}_j \rangle - \langle \mathbf{v}, \mathbf{u}_j \rangle = 0$$

    因此 $\mathbf{w}^\perp$ 与 $W$ 的基正交，即 $\mathbf{w}^\perp \in W^\perp$。

    **唯一性**：若 $\mathbf{v} = \mathbf{w}_1 + \mathbf{w}_1^\perp = \mathbf{w}_2 + \mathbf{w}_2^\perp$，则 $\mathbf{w}_1 - \mathbf{w}_2 = \mathbf{w}_2^\perp - \mathbf{w}_1^\perp$。左边属于 $W$，右边属于 $W^\perp$，因此 $\mathbf{w}_1 - \mathbf{w}_2 \in W \cap W^\perp$。对 $\mathbf{x} \in W \cap W^\perp$，$\langle \mathbf{x}, \mathbf{x} \rangle = 0$，故 $\mathbf{x} = \mathbf{0}$。即 $\mathbf{w}_1 = \mathbf{w}_2$，$\mathbf{w}_1^\perp = \mathbf{w}_2^\perp$。

    **维数关系**：$\dim(W) + \dim(W^\perp) = n$。 $\blacksquare$

!!! corollary "推论 7.1"
    $(W^\perp)^\perp = W$。

??? proof "证明"
    显然 $W \subseteq (W^\perp)^\perp$（$W$ 中的每个向量都与 $W^\perp$ 中的每个向量正交）。又 $\dim((W^\perp)^\perp) = n - \dim(W^\perp) = n - (n - \dim W) = \dim W$，维数相同的子空间且有包含关系，故相等。 $\blacksquare$

!!! theorem "定理 7.7 (行空间与零空间的正交关系)"
    设 $A$ 为 $m \times n$ 矩阵，则：

    1. $\operatorname{Row}(A)^\perp = \operatorname{Null}(A)$，即行空间的正交补是零空间。
    2. $\operatorname{Col}(A)^\perp = \operatorname{Null}(A^T)$，即列空间的正交补是左零空间。

??? proof "证明"
    **1.** $\mathbf{x} \in \operatorname{Null}(A)$ $\Leftrightarrow$ $A\mathbf{x} = \mathbf{0}$ $\Leftrightarrow$ $A$ 的每一行与 $\mathbf{x}$ 的内积为零 $\Leftrightarrow$ $\mathbf{x}$ 与 $A$ 的每一行正交 $\Leftrightarrow$ $\mathbf{x} \in \operatorname{Row}(A)^\perp$。

    **2.** 将 1 应用于 $A^T$：$\operatorname{Row}(A^T)^\perp = \operatorname{Null}(A^T)$，即 $\operatorname{Col}(A)^\perp = \operatorname{Null}(A^T)$。 $\blacksquare$

!!! example "例 7.3"
    设 $W = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}, \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}\right\} \subseteq \mathbb{R}^3$。求 $W^\perp$。

    $W^\perp = \operatorname{Null}(A)$，其中 $A = \begin{pmatrix} 1 & 1 & 1 \\ 1 & -1 & 0 \end{pmatrix}$。

    行化简：$\begin{pmatrix} 1 & 1 & 1 \\ 0 & -2 & -1 \end{pmatrix} \to \begin{pmatrix} 1 & 0 & 1/2 \\ 0 & 1 & 1/2 \end{pmatrix}$。

    $x_3 = t$ 自由，$x_1 = -t/2$，$x_2 = -t/2$。

    $$W^\perp = \operatorname{span}\left\{\begin{pmatrix} -1 \\ -1 \\ 2 \end{pmatrix}\right\} = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}\right\}$$

    验证：$\dim(W) + \dim(W^\perp) = 2 + 1 = 3 = \dim(\mathbb{R}^3)$。

---

## 7.4 正交投影

!!! definition "定义 7.7 (正交投影)"
    设 $W$ 为 $\mathbb{R}^n$ 的子空间。向量 $\mathbf{v}$ 在 $W$ 上的**正交投影**（orthogonal projection）定义为正交分解 $\mathbf{v} = \mathbf{w} + \mathbf{w}^\perp$ 中的 $\mathbf{w}$ 分量，记为

    $$\operatorname{proj}_W(\mathbf{v}) = \mathbf{w}$$

!!! theorem "定理 7.8 (正交投影公式)"
    设 $\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$ 为 $W$ 的一组正交基，则

    $$\operatorname{proj}_W(\mathbf{v}) = \sum_{i=1}^k \frac{\langle \mathbf{v}, \mathbf{u}_i \rangle}{\langle \mathbf{u}_i, \mathbf{u}_i \rangle}\mathbf{u}_i$$

    若 $\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$ 为标准正交基，则

    $$\operatorname{proj}_W(\mathbf{v}) = \sum_{i=1}^k \langle \mathbf{v}, \mathbf{u}_i \rangle \mathbf{u}_i$$

??? proof "证明"
    由定理 7.6 的证明中的构造即可得到。需要验证的核心是 $\mathbf{v} - \operatorname{proj}_W(\mathbf{v}) \in W^\perp$，这在该证明中已经完成。 $\blacksquare$

!!! theorem "定理 7.9 (最佳逼近定理)"
    设 $W$ 为 $\mathbb{R}^n$ 的子空间，$\mathbf{v} \in \mathbb{R}^n$。则 $\operatorname{proj}_W(\mathbf{v})$ 是 $W$ 中距离 $\mathbf{v}$ 最近的点。即对任意 $\mathbf{w} \in W$，

    $$\|\mathbf{v} - \operatorname{proj}_W(\mathbf{v})\| \leq \|\mathbf{v} - \mathbf{w}\|$$

    等号成立当且仅当 $\mathbf{w} = \operatorname{proj}_W(\mathbf{v})$。

??? proof "证明"
    设 $\hat{\mathbf{v}} = \operatorname{proj}_W(\mathbf{v})$。对任意 $\mathbf{w} \in W$，

    $$\|\mathbf{v} - \mathbf{w}\|^2 = \|(\mathbf{v} - \hat{\mathbf{v}}) + (\hat{\mathbf{v}} - \mathbf{w})\|^2$$

    注意 $\mathbf{v} - \hat{\mathbf{v}} \in W^\perp$ 且 $\hat{\mathbf{v}} - \mathbf{w} \in W$，两者正交，由勾股定理：

    $$= \|\mathbf{v} - \hat{\mathbf{v}}\|^2 + \|\hat{\mathbf{v}} - \mathbf{w}\|^2 \geq \|\mathbf{v} - \hat{\mathbf{v}}\|^2$$

    等号成立 $\Leftrightarrow$ $\|\hat{\mathbf{v}} - \mathbf{w}\| = 0$ $\Leftrightarrow$ $\mathbf{w} = \hat{\mathbf{v}}$。 $\blacksquare$

!!! definition "定义 7.8 (投影矩阵)"
    设 $W$ 的列为矩阵 $A$ 的列空间，即 $W = \operatorname{Col}(A)$。若 $A$ 的列线性无关，则从 $\mathbb{R}^n$ 到 $W$ 的正交投影矩阵为

    $$P = A(A^TA)^{-1}A^T$$

    该矩阵满足 $\operatorname{proj}_W(\mathbf{v}) = P\mathbf{v}$。

!!! proposition "命题 7.2 (投影矩阵的性质)"
    投影矩阵 $P = A(A^TA)^{-1}A^T$ 满足：

    1. $P^2 = P$（幂等性）
    2. $P^T = P$（对称性）
    3. $\operatorname{rank}(P) = \operatorname{rank}(A) = \dim(W)$

??? proof "证明"
    1. $P^2 = A(A^TA)^{-1}A^T \cdot A(A^TA)^{-1}A^T = A(A^TA)^{-1}(A^TA)(A^TA)^{-1}A^T = A(A^TA)^{-1}A^T = P$。

    2. $P^T = (A(A^TA)^{-1}A^T)^T = A((A^TA)^{-1})^T A^T = A((A^TA)^T)^{-1}A^T = A(A^TA)^{-1}A^T = P$（因为 $A^TA$ 对称，其逆也对称）。

    3. $P$ 的列空间等于 $\operatorname{Col}(A)$，故 $\operatorname{rank}(P) = \operatorname{rank}(A)$。 $\blacksquare$

!!! example "例 7.4"
    将 $\mathbf{v} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$ 投影到 $W = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}\right\}$ 上。

    注意这两个向量不正交（内积为 $0 + 0 + 1 = 1 \neq 0$），使用投影矩阵公式。

    $A = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{pmatrix}$，$A^TA = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$，$(A^TA)^{-1} = \frac{1}{3}\begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$。

    $$\hat{\mathbf{v}} = A(A^TA)^{-1}A^T\mathbf{v} = A \cdot \frac{1}{3}\begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}\begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$$

    $A^T\mathbf{v} = \begin{pmatrix} 4 \\ 5 \end{pmatrix}$，$(A^TA)^{-1}A^T\mathbf{v} = \frac{1}{3}\begin{pmatrix} 3 \\ 6 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$。

    $$\hat{\mathbf{v}} = A\begin{pmatrix} 1 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$$

    巧合的是 $\mathbf{v}$ 本身就在 $W$ 中！验证：$\mathbf{v} = 1 \cdot \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix} + 2 \cdot \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$。

---

## 7.5 Gram-Schmidt 正交化

!!! theorem "定理 7.10 (Gram-Schmidt 正交化过程)"
    设 $\{\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_k\}$ 为子空间 $W$ 的一组基。定义：

    $$\mathbf{v}_1 = \mathbf{x}_1$$

    $$\mathbf{v}_2 = \mathbf{x}_2 - \frac{\langle \mathbf{x}_2, \mathbf{v}_1 \rangle}{\langle \mathbf{v}_1, \mathbf{v}_1 \rangle}\mathbf{v}_1$$

    $$\mathbf{v}_3 = \mathbf{x}_3 - \frac{\langle \mathbf{x}_3, \mathbf{v}_1 \rangle}{\langle \mathbf{v}_1, \mathbf{v}_1 \rangle}\mathbf{v}_1 - \frac{\langle \mathbf{x}_3, \mathbf{v}_2 \rangle}{\langle \mathbf{v}_2, \mathbf{v}_2 \rangle}\mathbf{v}_2$$

    一般地，

    $$\mathbf{v}_j = \mathbf{x}_j - \sum_{i=1}^{j-1} \frac{\langle \mathbf{x}_j, \mathbf{v}_i \rangle}{\langle \mathbf{v}_i, \mathbf{v}_i \rangle}\mathbf{v}_i, \quad j = 2, 3, \ldots, k$$

    则 $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\}$ 是 $W$ 的正交基。进一步单位化 $\mathbf{q}_i = \frac{\mathbf{v}_i}{\|\mathbf{v}_i\|}$ 得到标准正交基。

    此外，对每个 $j$，$\operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_j\} = \operatorname{span}\{\mathbf{x}_1, \ldots, \mathbf{x}_j\}$。

??? proof "证明"
    对 $j$ 进行归纳。

    **基础**：$j = 1$，$\{\mathbf{v}_1\} = \{\mathbf{x}_1\}$ 是单元素集，自然正交。

    **归纳步**：假设 $\{\mathbf{v}_1, \ldots, \mathbf{v}_{j-1}\}$ 两两正交且 $\operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_{j-1}\} = \operatorname{span}\{\mathbf{x}_1, \ldots, \mathbf{x}_{j-1}\}$。

    由定义，$\mathbf{v}_j = \mathbf{x}_j - \operatorname{proj}_{W_{j-1}}(\mathbf{x}_j)$，其中 $W_{j-1} = \operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_{j-1}\}$。

    验证正交性：对任意 $\ell < j$，

    $$\langle \mathbf{v}_j, \mathbf{v}_\ell \rangle = \langle \mathbf{x}_j, \mathbf{v}_\ell \rangle - \sum_{i=1}^{j-1}\frac{\langle \mathbf{x}_j, \mathbf{v}_i \rangle}{\|\mathbf{v}_i\|^2}\langle \mathbf{v}_i, \mathbf{v}_\ell \rangle = \langle \mathbf{x}_j, \mathbf{v}_\ell \rangle - \frac{\langle \mathbf{x}_j, \mathbf{v}_\ell \rangle}{\|\mathbf{v}_\ell\|^2}\|\mathbf{v}_\ell\|^2 = 0$$

    最后一步用到了 $\langle \mathbf{v}_i, \mathbf{v}_\ell \rangle = 0$（$i \neq \ell$）。

    $\mathbf{v}_j \neq \mathbf{0}$：若 $\mathbf{v}_j = \mathbf{0}$，则 $\mathbf{x}_j = \operatorname{proj}_{W_{j-1}}(\mathbf{x}_j) \in W_{j-1} = \operatorname{span}\{\mathbf{x}_1, \ldots, \mathbf{x}_{j-1}\}$，与 $\{\mathbf{x}_1, \ldots, \mathbf{x}_k\}$ 线性无关矛盾。

    生成空间相等：$\mathbf{v}_j \in \operatorname{span}\{\mathbf{x}_1, \ldots, \mathbf{x}_j\}$（由定义），$\mathbf{x}_j \in \operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_j\}$（由定义可解出），故两个 span 相等。 $\blacksquare$

!!! example "例 7.5"
    对 $\mathbb{R}^3$ 中的基 $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$，$\mathbf{x}_2 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}$，$\mathbf{x}_3 = \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}$ 进行 Gram-Schmidt 正交化。

    **步骤 1**：$\mathbf{v}_1 = \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$。

    **步骤 2**：

    $$\frac{\langle \mathbf{x}_2, \mathbf{v}_1 \rangle}{\|\mathbf{v}_1\|^2} = \frac{1 \cdot 1 + 0 \cdot 1 + 1 \cdot 0}{1 + 1 + 0} = \frac{1}{2}$$

    $$\mathbf{v}_2 = \mathbf{x}_2 - \frac{1}{2}\mathbf{v}_1 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix} - \frac{1}{2}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix}$$

    验证：$\langle \mathbf{v}_1, \mathbf{v}_2 \rangle = 1/2 - 1/2 + 0 = 0$。

    **步骤 3**：

    $$\frac{\langle \mathbf{x}_3, \mathbf{v}_1 \rangle}{\|\mathbf{v}_1\|^2} = \frac{0 + 1 + 0}{2} = \frac{1}{2}, \quad \frac{\langle \mathbf{x}_3, \mathbf{v}_2 \rangle}{\|\mathbf{v}_2\|^2} = \frac{0 - 1/2 + 1}{1/4 + 1/4 + 1} = \frac{1/2}{3/2} = \frac{1}{3}$$

    $$\mathbf{v}_3 = \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix} - \frac{1}{2}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} - \frac{1}{3}\begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix} = \begin{pmatrix} -2/3 \\ 2/3 \\ 2/3 \end{pmatrix}$$

    验证：$\langle \mathbf{v}_1, \mathbf{v}_3 \rangle = -2/3 + 2/3 + 0 = 0$，$\langle \mathbf{v}_2, \mathbf{v}_3 \rangle = -1/3 - 1/3 + 2/3 = 0$。

    单位化：

    $$\mathbf{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad \mathbf{q}_2 = \frac{1}{\sqrt{3/2}}\begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix} = \sqrt{\frac{2}{3}}\begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix}, \quad \mathbf{q}_3 = \frac{1}{\sqrt{4/3}}\begin{pmatrix} -2/3 \\ 2/3 \\ 2/3 \end{pmatrix}$$

---

## 7.6 正交矩阵

!!! definition "定义 7.9 (正交矩阵)"
    实方阵 $Q$ 称为**正交矩阵**（orthogonal matrix），若

    $$Q^TQ = QQ^T = I$$

    即 $Q^{-1} = Q^T$。

!!! theorem "定理 7.11 (正交矩阵的等价刻画)"
    设 $Q$ 为 $n$ 阶实方阵，则以下条件等价：

    1. $Q$ 是正交矩阵。
    2. $Q$ 的列向量构成 $\mathbb{R}^n$ 的标准正交基。
    3. $Q$ 的行向量构成 $\mathbb{R}^n$ 的标准正交基。
    4. 对任意 $\mathbf{x} \in \mathbb{R}^n$，$\|Q\mathbf{x}\| = \|\mathbf{x}\|$（保范性 / 等距性）。
    5. 对任意 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$，$\langle Q\mathbf{x}, Q\mathbf{y} \rangle = \langle \mathbf{x}, \mathbf{y} \rangle$（保内积性）。

??? proof "证明"
    $(1) \Leftrightarrow (2)$：$Q^TQ = I$ 意味着 $(Q^TQ)_{ij} = \mathbf{q}_i^T\mathbf{q}_j = \delta_{ij}$（Kronecker delta），即列向量标准正交。

    $(1) \Leftrightarrow (3)$：$QQ^T = I$ 意味着行向量标准正交。

    $(1) \Rightarrow (5)$：$\langle Q\mathbf{x}, Q\mathbf{y} \rangle = (Q\mathbf{x})^T(Q\mathbf{y}) = \mathbf{x}^TQ^TQ\mathbf{y} = \mathbf{x}^T\mathbf{y} = \langle \mathbf{x}, \mathbf{y} \rangle$。

    $(5) \Rightarrow (4)$：取 $\mathbf{y} = \mathbf{x}$，$\|Q\mathbf{x}\|^2 = \langle Q\mathbf{x}, Q\mathbf{x} \rangle = \langle \mathbf{x}, \mathbf{x} \rangle = \|\mathbf{x}\|^2$。

    $(4) \Rightarrow (1)$：由 $\|Q\mathbf{x}\|^2 = \|\mathbf{x}\|^2$ 对所有 $\mathbf{x}$ 成立，即 $\mathbf{x}^TQ^TQ\mathbf{x} = \mathbf{x}^T\mathbf{x}$ 对所有 $\mathbf{x}$ 成立，故 $\mathbf{x}^T(Q^TQ - I)\mathbf{x} = 0$ 对所有 $\mathbf{x}$。由于 $Q^TQ - I$ 对称，这意味着 $Q^TQ - I = O$，即 $Q^TQ = I$。 $\blacksquare$

!!! theorem "定理 7.12 (正交矩阵的性质)"
    设 $Q$ 为正交矩阵，则：

    1. $\det(Q) = \pm 1$
    2. $Q^{-1} = Q^T$ 也是正交矩阵
    3. 若 $Q_1, Q_2$ 都是正交矩阵，则 $Q_1Q_2$ 也是正交矩阵
    4. $Q$ 的特征值的模为 $1$（即 $|\lambda| = 1$）

??? proof "证明"
    1. $\det(Q^TQ) = \det(I) = 1$。$\det(Q^T)\det(Q) = (\det Q)^2 = 1$，故 $\det Q = \pm 1$。

    2. $(Q^T)^T Q^T = QQ^T = I$。

    3. $(Q_1Q_2)^T(Q_1Q_2) = Q_2^TQ_1^TQ_1Q_2 = Q_2^TQ_2 = I$。

    4. 设 $Q\mathbf{v} = \lambda\mathbf{v}$（$\mathbf{v} \neq \mathbf{0}$，$\lambda \in \mathbb{C}$）。取共轭转置范数：$\|\mathbf{v}\|^2 = \overline{\mathbf{v}}^T\mathbf{v}$，$\|Q\mathbf{v}\|^2 = \overline{(Q\mathbf{v})}^T(Q\mathbf{v}) = |\lambda|^2\overline{\mathbf{v}}^T\mathbf{v} = |\lambda|^2\|\mathbf{v}\|^2$。由等距性 $\|Q\mathbf{v}\| = \|\mathbf{v}\|$，得 $|\lambda| = 1$。 $\blacksquare$

!!! example "例 7.6"
    旋转矩阵 $Q = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$ 是正交矩阵。

    验证：$Q^TQ = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$。

    $\det Q = \cos^2\theta + \sin^2\theta = 1 > 0$，因此这是一个**旋转**（保持方向的正交变换）。

    反射矩阵 $H = \begin{pmatrix} \cos 2\alpha & \sin 2\alpha \\ \sin 2\alpha & -\cos 2\alpha \end{pmatrix}$（关于角度 $\alpha$ 的直线反射）也是正交矩阵，但 $\det H = -1$（改变方向）。

---

## 7.7 最小二乘法

当线性方程组 $A\mathbf{x} = \mathbf{b}$ 无解（即 $\mathbf{b} \notin \operatorname{Col}(A)$）时，我们寻找使误差最小的近似解——这就是最小二乘法（least squares method）的核心思想。

!!! definition "定义 7.10 (最小二乘解)"
    设 $A$ 为 $m \times n$ 矩阵，$\mathbf{b} \in \mathbb{R}^m$。向量 $\hat{\mathbf{x}} \in \mathbb{R}^n$ 称为 $A\mathbf{x} = \mathbf{b}$ 的**最小二乘解**（least squares solution），若

    $$\|A\hat{\mathbf{x}} - \mathbf{b}\| \leq \|A\mathbf{x} - \mathbf{b}\|, \quad \forall\, \mathbf{x} \in \mathbb{R}^n$$

    即 $\hat{\mathbf{x}}$ 使残差 $\|\mathbf{b} - A\mathbf{x}\|$ 达到最小。

!!! theorem "定理 7.13 (法方程 / Normal Equations)"
    $\hat{\mathbf{x}}$ 是 $A\mathbf{x} = \mathbf{b}$ 的最小二乘解当且仅当 $\hat{\mathbf{x}}$ 满足**法方程**：

    $$A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$$

    若 $A$ 的列线性无关，则 $A^TA$ 可逆，最小二乘解唯一：

    $$\hat{\mathbf{x}} = (A^TA)^{-1}A^T\mathbf{b}$$

??? proof "证明"
    **几何论证**：$A\hat{\mathbf{x}}$ 是 $\mathbf{b}$ 在 $\operatorname{Col}(A)$ 上的正交投影。由最佳逼近定理（定理 7.9），$\|A\hat{\mathbf{x}} - \mathbf{b}\|$ 最小当且仅当 $A\hat{\mathbf{x}} = \operatorname{proj}_{\operatorname{Col}(A)}(\mathbf{b})$，即

    $$\mathbf{b} - A\hat{\mathbf{x}} \in \operatorname{Col}(A)^\perp = \operatorname{Null}(A^T)$$

    即 $A^T(\mathbf{b} - A\hat{\mathbf{x}}) = \mathbf{0}$，化简得 $A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$。

    **唯一性**：若 $A$ 的列线性无关，则 $\operatorname{Null}(A) = \{\mathbf{0}\}$。由 $A^TA\mathbf{x} = \mathbf{0}$ $\Rightarrow$ $\mathbf{x}^TA^TA\mathbf{x} = 0$ $\Rightarrow$ $\|A\mathbf{x}\|^2 = 0$ $\Rightarrow$ $A\mathbf{x} = \mathbf{0}$ $\Rightarrow$ $\mathbf{x} = \mathbf{0}$，故 $A^TA$ 可逆。 $\blacksquare$

!!! example "例 7.7"
    用最小二乘法拟合直线 $y = c_0 + c_1 x$ 通过数据点 $(1, 1), (2, 3), (3, 2), (4, 5)$。

    构建矩阵方程 $A\mathbf{c} = \mathbf{b}$：

    $$A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \\ 1 & 4 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 1 \\ 3 \\ 2 \\ 5 \end{pmatrix}, \quad \mathbf{c} = \begin{pmatrix} c_0 \\ c_1 \end{pmatrix}$$

    法方程 $A^TA\mathbf{c} = A^T\mathbf{b}$：

    $$A^TA = \begin{pmatrix} 4 & 10 \\ 10 & 30 \end{pmatrix}, \quad A^T\mathbf{b} = \begin{pmatrix} 11 \\ 33 \end{pmatrix}$$

    解方程：$\begin{pmatrix} 4 & 10 \\ 10 & 30 \end{pmatrix}\begin{pmatrix} c_0 \\ c_1 \end{pmatrix} = \begin{pmatrix} 11 \\ 33 \end{pmatrix}$。

    $\det(A^TA) = 120 - 100 = 20$。

    $$(A^TA)^{-1} = \frac{1}{20}\begin{pmatrix} 30 & -10 \\ -10 & 4 \end{pmatrix}$$

    $$\hat{\mathbf{c}} = \frac{1}{20}\begin{pmatrix} 30 & -10 \\ -10 & 4 \end{pmatrix}\begin{pmatrix} 11 \\ 33 \end{pmatrix} = \frac{1}{20}\begin{pmatrix} 330 - 330 \\ -110 + 132 \end{pmatrix} = \frac{1}{20}\begin{pmatrix} 0 \\ 22 \end{pmatrix} = \begin{pmatrix} 0 \\ 1.1 \end{pmatrix}$$

    最小二乘拟合直线为 $y = 1.1x$。

!!! example "例 7.8"
    **最小二乘的几何解释**：$A\mathbf{x} = \mathbf{b}$ 无解时，$\mathbf{b}$ 不在 $\operatorname{Col}(A)$ 中。最小二乘解使 $A\hat{\mathbf{x}} = \operatorname{proj}_{\operatorname{Col}(A)}(\mathbf{b})$，即将 $\mathbf{b}$ "投影" 到 $A$ 的列空间上。残差向量 $\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}}$ 垂直于 $\operatorname{Col}(A)$，这正是法方程 $A^T\mathbf{r} = \mathbf{0}$ 的含义。

    设 $A = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$，$\mathbf{b} = \begin{pmatrix} 1 \\ 2 \\ 4 \end{pmatrix}$。

    $A^TA = 3$，$A^T\mathbf{b} = 7$，$\hat{x} = 7/3$。

    $A\hat{x} = \begin{pmatrix} 7/3 \\ 7/3 \\ 7/3 \end{pmatrix}$，即 $1, 2, 4$ 的平均值 $7/3$，这就是"常数拟合"的最小二乘解。

---

## 7.8 QR 分解

!!! definition "定义 7.11 (QR 分解)"
    设 $A$ 为 $m \times n$ 矩阵（$m \geq n$）且列线性无关。$A$ 的 **QR 分解**（QR factorization / QR decomposition）是将 $A$ 写成

    $$A = QR$$

    其中 $Q$ 为 $m \times n$ 矩阵，其列构成标准正交集（$Q^TQ = I_n$），$R$ 为 $n \times n$ 上三角矩阵且对角元素为正。

!!! theorem "定理 7.14 (QR 分解的存在性与唯一性)"
    设 $A$ 为 $m \times n$ 矩阵（$m \geq n$），列线性无关，则 $A$ 存在唯一的 QR 分解 $A = QR$，其中 $Q$ 的列是标准正交的，$R$ 为上三角矩阵且对角元素为正。

??? proof "证明"
    **存在性**：对 $A$ 的列 $\{\mathbf{a}_1, \ldots, \mathbf{a}_n\}$ 施加 Gram-Schmidt 正交化，得到标准正交集 $\{\mathbf{q}_1, \ldots, \mathbf{q}_n\}$。

    由 Gram-Schmidt 过程的性质，对每个 $j$：

    $$\operatorname{span}\{\mathbf{q}_1, \ldots, \mathbf{q}_j\} = \operatorname{span}\{\mathbf{a}_1, \ldots, \mathbf{a}_j\}$$

    因此 $\mathbf{a}_j \in \operatorname{span}\{\mathbf{q}_1, \ldots, \mathbf{q}_j\}$，可以写成

    $$\mathbf{a}_j = r_{1j}\mathbf{q}_1 + r_{2j}\mathbf{q}_2 + \cdots + r_{jj}\mathbf{q}_j$$

    其中 $r_{ij} = \langle \mathbf{a}_j, \mathbf{q}_i \rangle$（$i \leq j$），$r_{ij} = 0$（$i > j$）。

    令 $Q = (\mathbf{q}_1 \;\; \cdots \;\; \mathbf{q}_n)$，$R = (r_{ij})$，则 $A = QR$，且 $R$ 为上三角矩阵。

    对角元素 $r_{jj} = \langle \mathbf{a}_j, \mathbf{q}_j \rangle$。由 Gram-Schmidt 过程的构造，$r_{jj} = \|\mathbf{v}_j\| > 0$（其中 $\mathbf{v}_j$ 是正交化后、单位化前的向量）。

    **唯一性**：设 $A = Q_1R_1 = Q_2R_2$ 为两个 QR 分解。则 $Q_2^TQ_1 = R_2R_1^{-1}$。左边的列是标准正交的（$Q_2^TQ_1$ 是正交矩阵），右边是上三角矩阵。正交且上三角的矩阵必为对角矩阵，且对角元素的模为 $1$。由 $R_1, R_2$ 对角元素为正，$R_2R_1^{-1}$ 对角元素也为正，故对角元素为 $1$，即 $Q_1 = Q_2$，$R_1 = R_2$。 $\blacksquare$

!!! theorem "定理 7.15 (利用 QR 分解求最小二乘解)"
    若 $A = QR$（QR 分解），则 $A\mathbf{x} = \mathbf{b}$ 的最小二乘解为

    $$\hat{\mathbf{x}} = R^{-1}Q^T\mathbf{b}$$

    即只需解上三角方程组 $R\hat{\mathbf{x}} = Q^T\mathbf{b}$。

??? proof "证明"
    由法方程 $A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$。代入 $A = QR$：

    $$(QR)^T(QR)\hat{\mathbf{x}} = (QR)^T\mathbf{b}$$
    $$R^TQ^TQR\hat{\mathbf{x}} = R^TQ^T\mathbf{b}$$
    $$R^TR\hat{\mathbf{x}} = R^TQ^T\mathbf{b}$$

    由 $R$ 可逆（对角元素为正），$R^T$ 可逆，左乘 $(R^T)^{-1}$：

    $$R\hat{\mathbf{x}} = Q^T\mathbf{b}$$

    这是上三角方程组，可以用回代法高效求解。 $\blacksquare$

!!! note "注"
    QR 分解在数值计算中比法方程更加稳定。法方程需要计算 $A^TA$，会放大舍入误差（条件数平方），而 QR 分解避免了这一问题。因此在实际数值计算中，QR 分解是求解最小二乘问题的首选方法。

!!! example "例 7.9"
    求 $A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \\ 0 & 1 \end{pmatrix}$ 的 QR 分解。

    $\mathbf{a}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$，$\mathbf{a}_2 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}$。

    **Gram-Schmidt**：

    $\mathbf{v}_1 = \mathbf{a}_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$，$\|\mathbf{v}_1\| = \sqrt{2}$。

    $\mathbf{v}_2 = \mathbf{a}_2 - \frac{\langle \mathbf{a}_2, \mathbf{v}_1 \rangle}{\|\mathbf{v}_1\|^2}\mathbf{v}_1 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix} - \frac{1}{2}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix}$，$\|\mathbf{v}_2\| = \sqrt{1/4 + 1/4 + 1} = \sqrt{3/2}$。

    **单位化**：

    $$\mathbf{q}_1 = \frac{\mathbf{v}_1}{\|\mathbf{v}_1\|} = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad \mathbf{q}_2 = \frac{\mathbf{v}_2}{\|\mathbf{v}_2\|} = \sqrt{\frac{2}{3}}\begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix} = \frac{1}{\sqrt{6}}\begin{pmatrix} 1 \\ -1 \\ 2 \end{pmatrix}$$

    **求 $R$**：

    $$r_{11} = \langle \mathbf{a}_1, \mathbf{q}_1 \rangle = \frac{1}{\sqrt{2}}(1 + 1) = \sqrt{2}$$

    $$r_{12} = \langle \mathbf{a}_2, \mathbf{q}_1 \rangle = \frac{1}{\sqrt{2}}(1 + 0) = \frac{1}{\sqrt{2}}$$

    $$r_{22} = \langle \mathbf{a}_2, \mathbf{q}_2 \rangle = \frac{1}{\sqrt{6}}(1 + 0 + 2) = \frac{3}{\sqrt{6}} = \sqrt{\frac{3}{2}}$$

    因此：

    $$Q = \begin{pmatrix} \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{6}} \\ \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{6}} \\ 0 & \frac{2}{\sqrt{6}} \end{pmatrix}, \quad R = \begin{pmatrix} \sqrt{2} & \frac{1}{\sqrt{2}} \\ 0 & \sqrt{\frac{3}{2}} \end{pmatrix}$$

    验证 $A = QR$：

    $$QR = \begin{pmatrix} \frac{1}{\sqrt{2}}\sqrt{2} + \frac{1}{\sqrt{6}} \cdot 0 & \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}} + \frac{1}{\sqrt{6}}\sqrt{\frac{3}{2}} \\ \frac{1}{\sqrt{2}}\sqrt{2} - \frac{1}{\sqrt{6}} \cdot 0 & \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}} - \frac{1}{\sqrt{6}}\sqrt{\frac{3}{2}} \\ 0 + \frac{2}{\sqrt{6}} \cdot 0 & 0 + \frac{2}{\sqrt{6}}\sqrt{\frac{3}{2}} \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 0 \\ 0 & 1 \end{pmatrix}$$

    验证通过。

!!! example "例 7.10"
    利用 QR 分解求解例 7.7 中的最小二乘问题。

    $A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \\ 1 & 4 \end{pmatrix}$。对列进行 Gram-Schmidt 正交化：

    $\mathbf{a}_1 = \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}$，$\|\mathbf{a}_1\| = 2$，$\mathbf{q}_1 = \frac{1}{2}\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}$。

    $\mathbf{v}_2 = \begin{pmatrix} 1 \\ 2 \\ 3 \\ 4 \end{pmatrix} - \frac{\langle \mathbf{a}_2, \mathbf{q}_1 \rangle}{\|\mathbf{q}_1\|^2}\mathbf{q}_1 \cdot \|\mathbf{q}_1\|$。先计算 $\langle \mathbf{a}_2, \mathbf{a}_1 \rangle / \|\mathbf{a}_1\|^2 = 10/4 = 5/2$。

    $\mathbf{v}_2 = \begin{pmatrix} 1 \\ 2 \\ 3 \\ 4 \end{pmatrix} - \frac{5}{2}\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} -3/2 \\ -1/2 \\ 1/2 \\ 3/2 \end{pmatrix}$，$\|\mathbf{v}_2\| = \sqrt{9/4 + 1/4 + 1/4 + 9/4} = \sqrt{5}$。

    $\mathbf{q}_2 = \frac{1}{\sqrt{5}}\begin{pmatrix} -3/2 \\ -1/2 \\ 1/2 \\ 3/2 \end{pmatrix}$。

    $R = \begin{pmatrix} 2 & 5/2 \cdot 2/1 \\ 0 & \sqrt{5} \end{pmatrix}$。更准确地说：$r_{11} = \|\mathbf{a}_1\| = 2$，$r_{12} = \langle \mathbf{a}_2, \mathbf{q}_1 \rangle = 10/2 = 5$，$r_{22} = \|\mathbf{v}_2\| = \sqrt{5}$。

    $R = \begin{pmatrix} 2 & 5 \\ 0 & \sqrt{5} \end{pmatrix}$。

    $Q^T\mathbf{b} = \begin{pmatrix} \mathbf{q}_1^T\mathbf{b} \\ \mathbf{q}_2^T\mathbf{b} \end{pmatrix} = \begin{pmatrix} (1+3+2+5)/2 \\ (-3/2-1/2\cdot 3+1/2\cdot 2+3/2\cdot 5)/\sqrt{5} \end{pmatrix} = \begin{pmatrix} 11/2 \\ (-3/2-3/2+1+15/2)/\sqrt{5} \end{pmatrix}$

    计算 $\mathbf{q}_2^T\mathbf{b}$：$\frac{1}{\sqrt{5}}(-3/2 \cdot 1 + (-1/2) \cdot 3 + 1/2 \cdot 2 + 3/2 \cdot 5) = \frac{1}{\sqrt{5}}(-3/2 - 3/2 + 1 + 15/2) = \frac{1}{\sqrt{5}} \cdot \frac{-3-3+2+15}{2} = \frac{11}{2\sqrt{5}}$。

    解 $R\hat{\mathbf{c}} = Q^T\mathbf{b}$：

    $$\begin{pmatrix} 2 & 5 \\ 0 & \sqrt{5} \end{pmatrix}\begin{pmatrix} c_0 \\ c_1 \end{pmatrix} = \begin{pmatrix} 11/2 \\ 11/(2\sqrt{5}) \end{pmatrix}$$

    从第二行：$\sqrt{5}\, c_1 = \frac{11}{2\sqrt{5}}$，$c_1 = \frac{11}{10} = 1.1$。

    从第一行：$2c_0 + 5 \cdot 1.1 = 5.5$，$2c_0 = 0$，$c_0 = 0$。

    结果 $\hat{\mathbf{c}} = \begin{pmatrix} 0 \\ 1.1 \end{pmatrix}$，与例 7.7 一致。

---

## 本章小结

本章系统地研究了 $\mathbb{R}^n$ 中的正交性理论及其应用，主要内容包括：

1. **内积与范数**：内积赋予向量空间几何结构，Cauchy-Schwarz 不等式和三角不等式是最基本的不等式。
2. **正交集**：由非零向量组成的正交集自动线性无关，正交基使得坐标计算极为简便。
3. **正交补**：$\mathbb{R}^n = W \oplus W^\perp$ 保证了正交分解的存在性和唯一性。
4. **正交投影**：投影公式 $P = A(A^TA)^{-1}A^T$ 给出了最佳逼近。
5. **Gram-Schmidt 正交化**：将任意基转化为正交基（标准正交基）的算法。
6. **正交矩阵**：$Q^TQ = I$ 的矩阵保持内积和范数，代表刚性变换（旋转或反射）。
7. **最小二乘法**：通过法方程 $A^TA\hat{\mathbf{x}} = A^T\mathbf{b}$ 求解超定方程组的最优近似解。
8. **QR 分解**：$A = QR$ 将矩阵分解为正交部分和上三角部分，是数值稳定的最小二乘算法的基础。
