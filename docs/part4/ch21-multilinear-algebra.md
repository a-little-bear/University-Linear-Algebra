# 第 21 章 多线性代数与张量

多线性代数是线性代数的自然推广，它研究多个向量空间之间的多线性关系。张量（tensor）作为多线性代数的核心对象，不仅为微分几何、广义相对论提供了基本语言，也在近年来的机器学习、量子计算和数据科学中扮演着越来越重要的角色。本章从对偶空间出发，系统建立多线性映射、张量积、外代数等基本理论，并介绍张量分解等现代应用方向。

## 21.1 对偶空间

### 线性泛函与对偶空间

!!! definition "定义 21.1 (线性泛函)"
    设 $V$ 是数域 $\mathbb{F}$（$\mathbb{F} = \mathbb{R}$ 或 $\mathbb{C}$）上的向量空间。一个**线性泛函**（linear functional）是从 $V$ 到 $\mathbb{F}$ 的线性映射 $f: V \to \mathbb{F}$，即满足：

    $$f(\alpha \mathbf{u} + \beta \mathbf{v}) = \alpha f(\mathbf{u}) + \beta f(\mathbf{v}), \quad \forall \mathbf{u}, \mathbf{v} \in V, \ \alpha, \beta \in \mathbb{F}.$$

!!! definition "定义 21.2 (对偶空间)"
    设 $V$ 是 $\mathbb{F}$ 上的向量空间。$V$ 的**对偶空间**（dual space）定义为 $V$ 上所有线性泛函构成的集合，记作 $V^*$：

    $$V^* = \text{Hom}(V, \mathbb{F}) = \{ f: V \to \mathbb{F} \mid f \text{ 是线性映射} \}.$$

    $V^*$ 在逐点加法和标量乘法下构成 $\mathbb{F}$ 上的向量空间。

### 对偶基

!!! definition "定义 21.3 (对偶基)"
    设 $V$ 是 $n$ 维向量空间，$\{\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\}$ 是 $V$ 的一组基。**对偶基**（dual basis）$\{\mathbf{e}^1, \mathbf{e}^2, \ldots, \mathbf{e}^n\}$ 是 $V^*$ 中满足以下条件的一组线性泛函：

    $$\mathbf{e}^i(\mathbf{e}_j) = \delta^i_j = \begin{cases} 1, & i = j, \\ 0, & i \neq j. \end{cases}$$

    其中 $\delta^i_j$ 是 Kronecker delta。

!!! theorem "定理 21.1 (对偶基的存在性与唯一性)"
    设 $V$ 是有限维向量空间，$\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 是 $V$ 的一组基，则对偶基 $\{\mathbf{e}^1, \ldots, \mathbf{e}^n\}$ 存在且唯一，并且构成 $V^*$ 的一组基。特别地，$\dim V^* = \dim V$。

??? proof "证明"
    **存在性**：对每个 $i = 1, \ldots, n$，定义 $\mathbf{e}^i: V \to \mathbb{F}$ 如下。对任意 $\mathbf{v} = \sum_{j=1}^n v_j \mathbf{e}_j \in V$，令

    $$\mathbf{e}^i(\mathbf{v}) = v_i.$$

    显然 $\mathbf{e}^i$ 是线性的，且 $\mathbf{e}^i(\mathbf{e}_j) = \delta^i_j$。

    **唯一性**：若 $f \in V^*$ 满足 $f(\mathbf{e}_j) = \delta^i_j$，则对任意 $\mathbf{v} = \sum_j v_j \mathbf{e}_j$，

    $$f(\mathbf{v}) = \sum_j v_j f(\mathbf{e}_j) = \sum_j v_j \delta^i_j = v_i = \mathbf{e}^i(\mathbf{v}).$$

    故 $f = \mathbf{e}^i$。

    **$\{\mathbf{e}^1, \ldots, \mathbf{e}^n\}$ 是 $V^*$ 的基**：

    - *线性无关*：设 $\sum_i c_i \mathbf{e}^i = 0$。将此等式作用于 $\mathbf{e}_j$，得 $c_j = 0$，$\forall j$。
    - *生成 $V^*$*：设 $f \in V^*$，令 $c_i = f(\mathbf{e}_i)$。则 $f$ 和 $\sum_i c_i \mathbf{e}^i$ 在基向量上的值相同，由线性性知它们相等，即 $f = \sum_i c_i \mathbf{e}^i$。

!!! example "例 21.1"
    设 $V = \mathbb{R}^3$，取标准基 $\{\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3\}$。对偶基为 $\{\mathbf{e}^1, \mathbf{e}^2, \mathbf{e}^3\}$，其中

    $$\mathbf{e}^i(x_1, x_2, x_3) = x_i.$$

    即 $\mathbf{e}^i$ 是"提取第 $i$ 个坐标"的投影函数。

    例如 $\mathbf{e}^2(3, -1, 5) = -1$。

### 双对偶空间与自然同构

!!! definition "定义 21.4 (双对偶空间)"
    $V$ 的**双对偶空间**（double dual space）定义为 $(V^*)^*$，记作 $V^{**}$。$V^{**}$ 的元素是 $V^*$ 上的线性泛函。

!!! theorem "定理 21.2 (自然同构)"
    设 $V$ 是有限维向量空间。定义映射 $\Phi: V \to V^{**}$，

    $$\Phi(\mathbf{v})(f) = f(\mathbf{v}), \quad \forall \mathbf{v} \in V, \ f \in V^*.$$

    则 $\Phi$ 是一个**自然同构**（canonical isomorphism），即 $\Phi$ 是不依赖于基的选取的向量空间同构。

??? proof "证明"
    **$\Phi$ 是线性映射**：对任意 $\alpha, \beta \in \mathbb{F}$，$\mathbf{u}, \mathbf{v} \in V$，$f \in V^*$，

    $$\Phi(\alpha \mathbf{u} + \beta \mathbf{v})(f) = f(\alpha \mathbf{u} + \beta \mathbf{v}) = \alpha f(\mathbf{u}) + \beta f(\mathbf{v}) = (\alpha \Phi(\mathbf{u}) + \beta \Phi(\mathbf{v}))(f).$$

    **$\Phi$ 是单射**：设 $\Phi(\mathbf{v}) = 0$，即对所有 $f \in V^*$，$f(\mathbf{v}) = 0$。取对偶基 $\mathbf{e}^i$，有 $\mathbf{e}^i(\mathbf{v}) = v_i = 0$，故 $\mathbf{v} = \mathbf{0}$。

    **$\Phi$ 是同构**：因为 $\dim V = \dim V^{**}$（$V$ 有限维），而 $\Phi$ 是单射，故 $\Phi$ 是同构。

    **自然性**：$\Phi$ 的定义不涉及任何基的选取，因此是"自然的"。在范畴论的意义下，$\Phi$ 构成从恒等函子到双对偶函子的自然变换。

!!! note "注"
    当 $V$ 是无限维空间时，$V$ 与 $V^*$ 不再同构（$V^*$ 可能"更大"），自然映射 $\Phi: V \to V^{**}$ 仍是单射但不再是满射。这一区别在泛函分析中尤为重要，引出了自反空间（reflexive space）的概念。

!!! example "例 21.2"
    设 $V = \mathbb{R}^2$，基为 $\{\mathbf{e}_1, \mathbf{e}_2\}$。对偶基为 $\{\mathbf{e}^1, \mathbf{e}^2\}$。对 $\mathbf{v} = 3\mathbf{e}_1 + 2\mathbf{e}_2$，

    $$\Phi(\mathbf{v})(\mathbf{e}^1) = \mathbf{e}^1(\mathbf{v}) = 3, \quad \Phi(\mathbf{v})(\mathbf{e}^2) = \mathbf{e}^2(\mathbf{v}) = 2.$$

    因此 $\Phi(\mathbf{v}) = 3(\mathbf{e}^1)^* + 2(\mathbf{e}^2)^*$，其中 $\{(\mathbf{e}^1)^*, (\mathbf{e}^2)^*\}$ 是 $V^{**}$ 中 $\{\mathbf{e}^1, \mathbf{e}^2\}$ 的对偶基。

## 21.2 多线性映射

### 双线性映射

!!! definition "定义 21.5 (双线性映射)"
    设 $V, W, U$ 是 $\mathbb{F}$ 上的向量空间。映射 $B: V \times W \to U$ 称为**双线性映射**（bilinear map），如果它关于每个变量都是线性的：

    - 固定 $\mathbf{w} \in W$，映射 $\mathbf{v} \mapsto B(\mathbf{v}, \mathbf{w})$ 是线性的；
    - 固定 $\mathbf{v} \in V$，映射 $\mathbf{w} \mapsto B(\mathbf{v}, \mathbf{w})$ 是线性的。

    即对任意 $\alpha, \beta \in \mathbb{F}$：

    $$B(\alpha \mathbf{v}_1 + \beta \mathbf{v}_2, \mathbf{w}) = \alpha B(\mathbf{v}_1, \mathbf{w}) + \beta B(\mathbf{v}_2, \mathbf{w}),$$

    $$B(\mathbf{v}, \alpha \mathbf{w}_1 + \beta \mathbf{w}_2) = \alpha B(\mathbf{v}, \mathbf{w}_1) + \beta B(\mathbf{v}, \mathbf{w}_2).$$

!!! example "例 21.3"
    **内积是双线性映射**。设 $V$ 是实内积空间，则内积 $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{R}$ 是双线性映射。

    **矩阵乘法是双线性映射**。映射 $B: \mathbb{R}^{m \times n} \times \mathbb{R}^{n \times p} \to \mathbb{R}^{m \times p}$，$B(A, B) = AB$ 是双线性的。

    **行列式作为行向量的函数**。将 $n \times n$ 矩阵看作 $n$ 个行向量 $(\mathbf{r}_1, \ldots, \mathbf{r}_n)$，则 $\det(\mathbf{r}_1, \ldots, \mathbf{r}_n)$ 关于每个 $\mathbf{r}_i$ 是线性的（即行列式是多线性映射）。

### 多线性映射

!!! definition "定义 21.6 (多线性映射)"
    设 $V_1, V_2, \ldots, V_k, U$ 是 $\mathbb{F}$ 上的向量空间。映射 $f: V_1 \times V_2 \times \cdots \times V_k \to U$ 称为 **$k$-线性映射**（$k$-linear map 或 multilinear map），如果固定其余 $k-1$ 个变量后，$f$ 关于每个变量都是线性的。

    所有从 $V_1 \times \cdots \times V_k$ 到 $U$ 的 $k$-线性映射构成的集合记为 $\mathcal{L}^k(V_1, \ldots, V_k; U)$。当 $U = \mathbb{F}$ 时，称为 **$k$-线性型**（$k$-linear form）。

!!! theorem "定理 21.3 (多线性映射由基上的值决定)"
    设 $V_i$ 是有限维向量空间，$\dim V_i = n_i$，$\{\mathbf{e}^{(i)}_1, \ldots, \mathbf{e}^{(i)}_{n_i}\}$ 是 $V_i$ 的基（$i = 1, \ldots, k$）。则 $k$-线性映射 $f: V_1 \times \cdots \times V_k \to U$ 完全由它在各基向量组上的值

    $$f(\mathbf{e}^{(1)}_{j_1}, \mathbf{e}^{(2)}_{j_2}, \ldots, \mathbf{e}^{(k)}_{j_k}), \quad 1 \leq j_i \leq n_i$$

    唯一确定。这些值共有 $n_1 n_2 \cdots n_k$ 个，可以自由指定。从而

    $$\dim \mathcal{L}^k(V_1, \ldots, V_k; U) = n_1 n_2 \cdots n_k \cdot \dim U.$$

??? proof "证明"
    设 $\mathbf{v}_i = \sum_{j_i=1}^{n_i} v^{j_i}_i \mathbf{e}^{(i)}_{j_i} \in V_i$。由多线性性，

    $$f(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \sum_{j_1=1}^{n_1} \cdots \sum_{j_k=1}^{n_k} v^{j_1}_1 \cdots v^{j_k}_k \, f(\mathbf{e}^{(1)}_{j_1}, \ldots, \mathbf{e}^{(k)}_{j_k}).$$

    这表明 $f$ 完全由 $f(\mathbf{e}^{(1)}_{j_1}, \ldots, \mathbf{e}^{(k)}_{j_k})$ 的取值决定。

    反之，对任意给定的值 $\mathbf{u}_{j_1 \cdots j_k} \in U$，可以用上述公式定义 $f$，直接验证可知 $f$ 是 $k$-线性的。

    当 $U = \mathbb{F}$ 时，$\dim \mathcal{L}^k(V_1, \ldots, V_k; \mathbb{F}) = n_1 n_2 \cdots n_k$。

## 21.3 张量积

### 泛性质定义

!!! definition "定义 21.7 (张量积的泛性质)"
    设 $V, W$ 是 $\mathbb{F}$ 上的向量空间。$V$ 和 $W$ 的**张量积**（tensor product）是一个向量空间 $V \otimes W$ 连同一个双线性映射 $\otimes: V \times W \to V \otimes W$，使得以下泛性质成立：

    对于任意向量空间 $U$ 和任意双线性映射 $B: V \times W \to U$，存在唯一的线性映射 $\tilde{B}: V \otimes W \to U$，使得 $B = \tilde{B} \circ \otimes$，即以下图表交换：

    $$V \times W \xrightarrow{\otimes} V \otimes W$$

    $$\searrow^{B} \quad \downarrow^{\tilde{B}}$$

    $$\qquad \quad U$$

    元素 $\mathbf{v} \otimes \mathbf{w}$（$\mathbf{v} \in V, \mathbf{w} \in W$）称为**简单张量**（simple tensor 或 decomposable tensor）。

!!! theorem "定理 21.4 (张量积的存在性与唯一性)"
    对任意向量空间 $V, W$，张量积 $V \otimes W$ 存在且在同构意义下唯一。

??? proof "证明"
    **存在性（构造法）**：设 $F(V \times W)$ 是以 $V \times W$ 中所有元素 $(\mathbf{v}, \mathbf{w})$ 为基的自由向量空间。设 $R$ 是由以下形式的元素生成的子空间：

    - $(\alpha \mathbf{v}_1 + \beta \mathbf{v}_2, \mathbf{w}) - \alpha(\mathbf{v}_1, \mathbf{w}) - \beta(\mathbf{v}_2, \mathbf{w})$
    - $(\mathbf{v}, \alpha \mathbf{w}_1 + \beta \mathbf{w}_2) - \alpha(\mathbf{v}, \mathbf{w}_1) - \beta(\mathbf{v}, \mathbf{w}_2)$

    定义 $V \otimes W = F(V \times W) / R$，并令 $\mathbf{v} \otimes \mathbf{w}$ 为 $(\mathbf{v}, \mathbf{w})$ 在商空间中的像。则映射 $\otimes: V \times W \to V \otimes W$ 是双线性的，且满足泛性质。

    **唯一性**：设 $(T_1, \otimes_1)$ 和 $(T_2, \otimes_2)$ 均满足泛性质。由 $T_1$ 的泛性质，对双线性映射 $\otimes_2$，存在唯一的线性映射 $\varphi: T_1 \to T_2$，使得 $\otimes_2 = \varphi \circ \otimes_1$。同理存在 $\psi: T_2 \to T_1$，使得 $\otimes_1 = \psi \circ \otimes_2$。则 $\psi \circ \varphi \circ \otimes_1 = \otimes_1$，由泛性质的唯一性，$\psi \circ \varphi = \text{id}_{T_1}$。同理 $\varphi \circ \psi = \text{id}_{T_2}$。故 $T_1 \cong T_2$。

### 张量积的基与维数

!!! theorem "定理 21.5 (张量积的基)"
    设 $\{\mathbf{e}_1, \ldots, \mathbf{e}_m\}$ 是 $V$ 的基，$\{\mathbf{f}_1, \ldots, \mathbf{f}_n\}$ 是 $W$ 的基，则

    $$\{\mathbf{e}_i \otimes \mathbf{f}_j : 1 \leq i \leq m, \ 1 \leq j \leq n\}$$

    是 $V \otimes W$ 的基。特别地，

    $$\dim(V \otimes W) = \dim V \cdot \dim W.$$

??? proof "证明"
    **生成**：$V \otimes W$ 由简单张量 $\mathbf{v} \otimes \mathbf{w}$ 的线性组合生成。设 $\mathbf{v} = \sum_i a_i \mathbf{e}_i$，$\mathbf{w} = \sum_j b_j \mathbf{f}_j$，则

    $$\mathbf{v} \otimes \mathbf{w} = \left(\sum_i a_i \mathbf{e}_i\right) \otimes \left(\sum_j b_j \mathbf{f}_j\right) = \sum_{i,j} a_i b_j \, \mathbf{e}_i \otimes \mathbf{f}_j.$$

    **线性无关**：设 $\sum_{i,j} c_{ij} \, \mathbf{e}_i \otimes \mathbf{f}_j = 0$。对每对 $(p, q)$，取线性泛函 $\mathbf{e}^p \in V^*$，$\mathbf{f}^q \in W^*$，定义双线性映射 $B_{pq}(\mathbf{v}, \mathbf{w}) = \mathbf{e}^p(\mathbf{v}) \mathbf{f}^q(\mathbf{w})$。由泛性质，存在线性映射 $\tilde{B}_{pq}: V \otimes W \to \mathbb{F}$，使得 $\tilde{B}_{pq}(\mathbf{e}_i \otimes \mathbf{f}_j) = \delta^p_i \delta^q_j$。将 $\tilde{B}_{pq}$ 应用于等式两端，得 $c_{pq} = 0$。

!!! example "例 21.4"
    设 $V = \mathbb{R}^2$，$W = \mathbb{R}^3$。则 $V \otimes W$ 是 $6$ 维空间。取 $V$ 的标准基 $\{\mathbf{e}_1, \mathbf{e}_2\}$，$W$ 的标准基 $\{\mathbf{f}_1, \mathbf{f}_2, \mathbf{f}_3\}$，则 $V \otimes W$ 的基为：

    $$\{\mathbf{e}_1 \otimes \mathbf{f}_1, \mathbf{e}_1 \otimes \mathbf{f}_2, \mathbf{e}_1 \otimes \mathbf{f}_3, \mathbf{e}_2 \otimes \mathbf{f}_1, \mathbf{e}_2 \otimes \mathbf{f}_2, \mathbf{e}_2 \otimes \mathbf{f}_3\}.$$

    $V \otimes W$ 与 $\mathbb{R}^{2 \times 3}$（$2 \times 3$ 矩阵空间）同构：$\mathbf{v} \otimes \mathbf{w} \mapsto \mathbf{v}\mathbf{w}^T$。

    例如 $\begin{pmatrix} 1 \\ 2 \end{pmatrix} \otimes \begin{pmatrix} 3 \\ 0 \\ -1 \end{pmatrix} \mapsto \begin{pmatrix} 3 & 0 & -1 \\ 6 & 0 & -2 \end{pmatrix}$。

!!! proposition "命题 21.1 (张量积的基本性质)"
    张量积满足以下性质：

    1. **结合律**：$(V \otimes W) \otimes U \cong V \otimes (W \otimes U)$；
    2. **交换律**：$V \otimes W \cong W \otimes V$；
    3. **分配律**：$V \otimes (W \oplus U) \cong (V \otimes W) \oplus (V \otimes U)$；
    4. **与标量域的张量积**：$V \otimes \mathbb{F} \cong V$；
    5. **与对偶空间的关系**：$V^* \otimes W \cong \text{Hom}(V, W)$。

??? proof "证明"
    我们证明第 5 条。定义映射 $\Phi: V^* \otimes W \to \text{Hom}(V, W)$，在简单张量上定义为

    $$\Phi(f \otimes \mathbf{w})(\mathbf{v}) = f(\mathbf{v}) \mathbf{w}, \quad f \in V^*, \ \mathbf{w} \in W, \ \mathbf{v} \in V.$$

    由张量积的泛性质，映射 $(f, \mathbf{w}) \mapsto [\ \mathbf{v} \mapsto f(\mathbf{v})\mathbf{w}\ ]$ 是双线性的，故 $\Phi$ 存在且是线性的。

    取 $V$ 的基 $\{\mathbf{e}_i\}$ 和 $W$ 的基 $\{\mathbf{f}_j\}$，对偶基 $\{\mathbf{e}^i\}$。则 $\Phi(\mathbf{e}^i \otimes \mathbf{f}_j)$ 是将 $\mathbf{e}_i$ 映为 $\mathbf{f}_j$，其余基向量映为 $\mathbf{0}$ 的线性映射，这恰好是 $\text{Hom}(V, W)$ 的标准基。因此 $\Phi$ 将基映为基，是同构。

## 21.4 张量的分量表示

### 指标记法

在物理学和工程学中，张量通常用**指标记法**（index notation）来表示。设 $V$ 是 $n$ 维向量空间，基为 $\{\mathbf{e}_i\}$，对偶基为 $\{\mathbf{e}^i\}$。

!!! definition "定义 21.8 ($(r, s)$ 型张量)"
    一个 **$(r, s)$ 型张量**（type $(r, s)$ tensor）是一个多线性映射

    $$T: \underbrace{V^* \times \cdots \times V^*}_{r} \times \underbrace{V \times \cdots \times V}_{s} \to \mathbb{F}.$$

    在给定基下，$T$ 的**分量**（components）为

    $$T^{i_1 \cdots i_r}_{j_1 \cdots j_s} = T(\mathbf{e}^{i_1}, \ldots, \mathbf{e}^{i_r}, \mathbf{e}_{j_1}, \ldots, \mathbf{e}_{j_s}).$$

    上标称为**逆变指标**（contravariant indices），下标称为**协变指标**（covariant indices）。

### Einstein 求和约定

!!! definition "定义 21.9 (Einstein 求和约定)"
    **Einstein 求和约定**（Einstein summation convention）规定：当一个指标在表达式中同时作为上标和下标出现时（称为**哑指标**），自动对该指标求和。例如：

    $$a^i b_i \equiv \sum_{i=1}^n a^i b_i, \quad T^i_{\ j} v^j \equiv \sum_{j=1}^n T^i_{\ j} v^j.$$

!!! example "例 21.5"
    使用 Einstein 约定，以下是一些常见的张量运算：

    1. **向量的分量**：$\mathbf{v} = v^i \mathbf{e}_i$（对 $i$ 求和）。
    2. **线性泛函的作用**：$f(\mathbf{v}) = f_i v^i$。
    3. **矩阵乘向量**：$(A\mathbf{v})^i = A^i_{\ j} v^j$。
    4. **矩阵乘法**：$(AB)^i_{\ k} = A^i_{\ j} B^j_{\ k}$。
    5. **迹**：$\text{tr}(A) = A^i_{\ i}$。

!!! theorem "定理 21.6 (张量分量的变换律)"
    设 $\{\mathbf{e}_i\}$ 和 $\{\tilde{\mathbf{e}}_i\}$ 是 $V$ 的两组基，基变换矩阵为 $P$：$\tilde{\mathbf{e}}_j = P^i_{\ j} \mathbf{e}_i$。则 $(r, s)$ 型张量 $T$ 的分量变换律为：

    $$\tilde{T}^{i_1 \cdots i_r}_{j_1 \cdots j_s} = (P^{-1})^{i_1}_{\ k_1} \cdots (P^{-1})^{i_r}_{\ k_r} \, P^{l_1}_{\ j_1} \cdots P^{l_s}_{\ j_s} \, T^{k_1 \cdots k_r}_{l_1 \cdots l_s}.$$

??? proof "证明"
    对偶基的变换为 $\tilde{\mathbf{e}}^i = (P^{-1})^i_{\ j} \mathbf{e}^j$（利用 $\tilde{\mathbf{e}}^i(\tilde{\mathbf{e}}_k) = \delta^i_k$ 可导出）。因此

    $$\tilde{T}^{i_1 \cdots i_r}_{j_1 \cdots j_s} = T(\tilde{\mathbf{e}}^{i_1}, \ldots, \tilde{\mathbf{e}}^{i_r}, \tilde{\mathbf{e}}_{j_1}, \ldots, \tilde{\mathbf{e}}_{j_s}).$$

    将 $\tilde{\mathbf{e}}^{i_\alpha} = (P^{-1})^{i_\alpha}_{\ k_\alpha} \mathbf{e}^{k_\alpha}$ 和 $\tilde{\mathbf{e}}_{j_\beta} = P^{l_\beta}_{\ j_\beta} \mathbf{e}_{l_\beta}$ 代入，由 $T$ 的多线性性即得结果。

## 21.5 对称张量与反对称张量

!!! definition "定义 21.10 (对称张量与反对称张量)"
    设 $T \in V^{\otimes k} = \underbrace{V \otimes \cdots \otimes V}_{k}$ 是一个 $k$ 阶张量。

    - $T$ 是**对称的**（symmetric），如果对任意置换 $\sigma \in S_k$，$T(\mathbf{v}_{\sigma(1)}, \ldots, \mathbf{v}_{\sigma(k)}) = T(\mathbf{v}_1, \ldots, \mathbf{v}_k)$；分量上即 $T_{i_{\sigma(1)} \cdots i_{\sigma(k)}} = T_{i_1 \cdots i_k}$。
    - $T$ 是**反对称的**（antisymmetric / alternating），如果 $T(\mathbf{v}_{\sigma(1)}, \ldots, \mathbf{v}_{\sigma(k)}) = \text{sgn}(\sigma) \, T(\mathbf{v}_1, \ldots, \mathbf{v}_k)$；分量上即 $T_{i_{\sigma(1)} \cdots i_{\sigma(k)}} = \text{sgn}(\sigma) \, T_{i_1 \cdots i_k}$。

!!! definition "定义 21.11 (对称化与反对称化算子)"
    **对称化算子**（symmetrization operator）$\text{Sym}$ 和**反对称化算子**（antisymmetrization operator）$\text{Alt}$ 定义为：

    $$\text{Sym}(T)(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \frac{1}{k!} \sum_{\sigma \in S_k} T(\mathbf{v}_{\sigma(1)}, \ldots, \mathbf{v}_{\sigma(k)}),$$

    $$\text{Alt}(T)(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \frac{1}{k!} \sum_{\sigma \in S_k} \text{sgn}(\sigma) \, T(\mathbf{v}_{\sigma(1)}, \ldots, \mathbf{v}_{\sigma(k)}).$$

!!! theorem "定理 21.7 (对称化与反对称化的性质)"
    1. $\text{Sym}$ 和 $\text{Alt}$ 是从 $V^{\otimes k}$ 到自身的投影算子（即 $\text{Sym}^2 = \text{Sym}$，$\text{Alt}^2 = \text{Alt}$）。
    2. $\text{Sym}(T)$ 是对称张量，$\text{Alt}(T)$ 是反对称张量。
    3. $T$ 是对称的当且仅当 $\text{Sym}(T) = T$；$T$ 是反对称的当且仅当 $\text{Alt}(T) = T$。

??? proof "证明"
    证明 $\text{Alt}$ 是投影算子。设 $S = \text{Alt}(T)$，则

    $$\text{Alt}(S)(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \frac{1}{k!} \sum_{\tau \in S_k} \text{sgn}(\tau) S(\mathbf{v}_{\tau(1)}, \ldots, \mathbf{v}_{\tau(k)}).$$

    由于 $S = \text{Alt}(T)$ 是反对称的，$S(\mathbf{v}_{\tau(1)}, \ldots, \mathbf{v}_{\tau(k)}) = \text{sgn}(\tau) S(\mathbf{v}_1, \ldots, \mathbf{v}_k)$。故

    $$\text{Alt}(S)(\mathbf{v}_1, \ldots, \mathbf{v}_k) = \frac{1}{k!} \sum_{\tau \in S_k} \text{sgn}(\tau)^2 S(\mathbf{v}_1, \ldots, \mathbf{v}_k) = S(\mathbf{v}_1, \ldots, \mathbf{v}_k).$$

    即 $\text{Alt}^2 = \text{Alt}$。其余性质类似可证。

!!! example "例 21.6"
    设 $T: \mathbb{R}^2 \times \mathbb{R}^2 \to \mathbb{R}$ 定义为 $T(\mathbf{u}, \mathbf{v}) = u_1 v_2$。则：

    $$\text{Sym}(T)(\mathbf{u}, \mathbf{v}) = \frac{1}{2}(u_1 v_2 + u_2 v_1),$$

    $$\text{Alt}(T)(\mathbf{u}, \mathbf{v}) = \frac{1}{2}(u_1 v_2 - u_2 v_1).$$

    注意 $\text{Alt}(T)(\mathbf{u}, \mathbf{v})$ 恰好是 $\frac{1}{2} \det \begin{pmatrix} u_1 & u_2 \\ v_1 & v_2 \end{pmatrix}$。

## 21.6 外代数

### 楔积与外幂

!!! definition "定义 21.12 (外积/楔积)"
    设 $V$ 是 $n$ 维向量空间。对 $\omega \in \Lambda^p(V^*)$，$\eta \in \Lambda^q(V^*)$（其中 $\Lambda^k(V^*)$ 是 $V^*$ 上的 $k$ 阶反对称张量空间），定义**楔积**（wedge product）$\omega \wedge \eta \in \Lambda^{p+q}(V^*)$ 为：

    $$\omega \wedge \eta = \frac{(p+q)!}{p! \, q!} \, \text{Alt}(\omega \otimes \eta).$$

    对于向量空间 $V$ 本身，**外幂**（exterior power）$\Lambda^k(V)$ 定义为 $V^{\otimes k}$ 中所有反对称张量构成的子空间。对 $\mathbf{v}_1, \ldots, \mathbf{v}_k \in V$，定义

    $$\mathbf{v}_1 \wedge \mathbf{v}_2 \wedge \cdots \wedge \mathbf{v}_k = \sum_{\sigma \in S_k} \text{sgn}(\sigma) \, \mathbf{v}_{\sigma(1)} \otimes \cdots \otimes \mathbf{v}_{\sigma(k)}.$$

!!! theorem "定理 21.8 (楔积的性质)"
    楔积满足以下性质：

    1. **双线性性**：$(\alpha \omega_1 + \beta \omega_2) \wedge \eta = \alpha (\omega_1 \wedge \eta) + \beta (\omega_2 \wedge \eta)$；
    2. **结合律**：$(\omega \wedge \eta) \wedge \zeta = \omega \wedge (\eta \wedge \zeta)$；
    3. **分次交换律**：$\omega \wedge \eta = (-1)^{pq} \eta \wedge \omega$，其中 $\omega \in \Lambda^p$，$\eta \in \Lambda^q$；
    4. **反对称性**：$\mathbf{v} \wedge \mathbf{v} = 0$，$\forall \mathbf{v} \in V$。

??? proof "证明"
    证明第 3 条。设 $\sigma_0$ 是将 $(1, \ldots, p, p+1, \ldots, p+q)$ 变为 $(p+1, \ldots, p+q, 1, \ldots, p)$ 的置换。这个置换需要 $pq$ 次对换（将 $q$ 个元素依次移过 $p$ 个元素），因此 $\text{sgn}(\sigma_0) = (-1)^{pq}$。

    由反对称张量的性质：

    $$(\omega \wedge \eta)(\mathbf{v}_{\sigma_0(1)}, \ldots, \mathbf{v}_{\sigma_0(p+q)}) = (-1)^{pq} (\omega \wedge \eta)(\mathbf{v}_1, \ldots, \mathbf{v}_{p+q}).$$

    而 $(\mathbf{v}_{\sigma_0(1)}, \ldots, \mathbf{v}_{\sigma_0(p+q)}) = (\mathbf{v}_{p+1}, \ldots, \mathbf{v}_{p+q}, \mathbf{v}_1, \ldots, \mathbf{v}_p)$，由楔积定义，这恰好给出 $(\eta \wedge \omega)(\mathbf{v}_1, \ldots, \mathbf{v}_{p+q})$。

!!! theorem "定理 21.9 (外幂的维数)"
    设 $\dim V = n$，则

    $$\dim \Lambda^k(V) = \binom{n}{k}.$$

    特别地，$\Lambda^0(V) = \mathbb{F}$，$\Lambda^1(V) = V$，$\Lambda^n(V)$ 是 $1$ 维的，$\Lambda^k(V) = \{0\}$（$k > n$）。

??? proof "证明"
    设 $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 是 $V$ 的基。则

    $$\{\mathbf{e}_{i_1} \wedge \mathbf{e}_{i_2} \wedge \cdots \wedge \mathbf{e}_{i_k} : 1 \leq i_1 < i_2 < \cdots < i_k \leq n\}$$

    构成 $\Lambda^k(V)$ 的基（可直接验证线性无关和生成）。这些基向量的个数为 $\binom{n}{k}$。

    当 $k > n$ 时，在 $\mathbf{e}_{i_1} \wedge \cdots \wedge \mathbf{e}_{i_k}$ 中必有重复的基向量，由反对称性得 $\mathbf{e}_{i_1} \wedge \cdots \wedge \mathbf{e}_{i_k} = 0$。

!!! definition "定义 21.13 (外代数)"
    $V$ 上的**外代数**（exterior algebra）定义为

    $$\Lambda(V) = \bigoplus_{k=0}^{n} \Lambda^k(V),$$

    配以楔积运算。$\Lambda(V)$ 是一个分次反交换代数（graded anti-commutative algebra），总维数为 $\sum_{k=0}^n \binom{n}{k} = 2^n$。

### 行列式作为外积

!!! theorem "定理 21.10 (行列式的外积解释)"
    设 $\mathbf{v}_1, \ldots, \mathbf{v}_n \in V$（$\dim V = n$），且 $\mathbf{v}_i = \sum_j a_{ji} \mathbf{e}_j$（即 $a_{ji}$ 是 $\mathbf{v}_i$ 的第 $j$ 个分量）。则

    $$\mathbf{v}_1 \wedge \mathbf{v}_2 \wedge \cdots \wedge \mathbf{v}_n = \det(A) \, \mathbf{e}_1 \wedge \mathbf{e}_2 \wedge \cdots \wedge \mathbf{e}_n,$$

    其中 $A = (a_{ij})$ 是以 $\mathbf{v}_j$ 为列的矩阵。

??? proof "证明"
    由楔积的多线性性和反对称性：

    $$\mathbf{v}_1 \wedge \cdots \wedge \mathbf{v}_n = \left(\sum_{i_1} a_{i_1 1} \mathbf{e}_{i_1}\right) \wedge \cdots \wedge \left(\sum_{i_n} a_{i_n n} \mathbf{e}_{i_n}\right)$$

    $$= \sum_{i_1, \ldots, i_n} a_{i_1 1} \cdots a_{i_n n} \, \mathbf{e}_{i_1} \wedge \cdots \wedge \mathbf{e}_{i_n}.$$

    由于 $\mathbf{e}_{i_1} \wedge \cdots \wedge \mathbf{e}_{i_n} = 0$ 若 $(i_1, \ldots, i_n)$ 中有重复指标，非零项仅当 $(i_1, \ldots, i_n)$ 是 $(1, \ldots, n)$ 的置换 $\sigma$ 时出现。此时 $\mathbf{e}_{\sigma(1)} \wedge \cdots \wedge \mathbf{e}_{\sigma(n)} = \text{sgn}(\sigma) \, \mathbf{e}_1 \wedge \cdots \wedge \mathbf{e}_n$。故

    $$\mathbf{v}_1 \wedge \cdots \wedge \mathbf{v}_n = \left(\sum_{\sigma \in S_n} \text{sgn}(\sigma) \, a_{\sigma(1),1} \cdots a_{\sigma(n),n}\right) \mathbf{e}_1 \wedge \cdots \wedge \mathbf{e}_n = \det(A) \, \mathbf{e}_1 \wedge \cdots \wedge \mathbf{e}_n.$$

!!! example "例 21.7"
    设 $V = \mathbb{R}^3$，$\mathbf{v}_1 = \begin{pmatrix} 1 \\ 2 \\ 0 \end{pmatrix}$，$\mathbf{v}_2 = \begin{pmatrix} 0 \\ 1 \\ 3 \end{pmatrix}$。则

    $$\mathbf{v}_1 \wedge \mathbf{v}_2 = (1 \cdot 1 - 2 \cdot 0)(\mathbf{e}_1 \wedge \mathbf{e}_2) + (1 \cdot 3 - 0 \cdot 0)(\mathbf{e}_1 \wedge \mathbf{e}_3) + (2 \cdot 3 - 0 \cdot 1)(\mathbf{e}_2 \wedge \mathbf{e}_3)$$

    $$= \mathbf{e}_1 \wedge \mathbf{e}_2 + 3 \, \mathbf{e}_1 \wedge \mathbf{e}_3 + 6 \, \mathbf{e}_2 \wedge \mathbf{e}_3.$$

    这对应于 $\mathbf{v}_1 \times \mathbf{v}_2 = (6, -3, 1)^T$（外积与叉积的联系）。

## 21.7 张量分解

张量分解（tensor decomposition）是将高阶张量表示为简单张量之和的技术，在信号处理、机器学习和化学计量学中有广泛应用。

### CP 分解

!!! definition "定义 21.14 (CP 分解)"
    **CP 分解**（Canonical Polyadic decomposition，也称 CANDECOMP/PARAFAC 分解）将一个 $k$ 阶张量 $\mathcal{T} \in \mathbb{R}^{n_1 \times n_2 \times \cdots \times n_k}$ 表示为秩一张量之和：

    $$\mathcal{T} = \sum_{r=1}^{R} \lambda_r \, \mathbf{a}^{(1)}_r \otimes \mathbf{a}^{(2)}_r \otimes \cdots \otimes \mathbf{a}^{(k)}_r,$$

    其中 $\lambda_r \in \mathbb{R}$，$\mathbf{a}^{(i)}_r \in \mathbb{R}^{n_i}$。使上式成立的最小正整数 $R$ 称为 $\mathcal{T}$ 的**张量秩**（tensor rank）。

!!! theorem "定理 21.11 (张量秩与矩阵秩的区别)"
    与矩阵秩不同，张量秩具有以下性质：

    1. 张量秩可能超过其任意模展开（mode unfolding）的矩阵秩；
    2. 在实数域和复数域上，同一张量的秩可能不同；
    3. 计算张量秩是 NP-hard 问题。

??? proof "证明"
    第 1 条的经典反例：考虑 $2 \times 2 \times 2$ 张量 $\mathcal{T}$，其分量为 $\mathcal{T}_{111} = \mathcal{T}_{222} = 1$，其余分量为 $0$。可验证其三个模展开的矩阵秩均为 $2$，但 $\mathcal{T}$ 的张量秩为 $2$（在实数域上）。

    第 2 条的例子：存在张量使得实秩为 $3$ 但复秩为 $2$。例如某些 $2 \times 2 \times 2$ 张量的实秩为 $3$，而在 $\mathbb{C}$ 上可分解为 $2$ 个秩一张量之和。

    第 3 条已由 Hillar 和 Lim (2013) 严格证明。

### Tucker 分解

!!! definition "定义 21.15 (Tucker 分解)"
    **Tucker 分解**（Tucker decomposition）将 $k$ 阶张量表示为一个核心张量（core tensor）与各模因子矩阵的乘积：

    $$\mathcal{T} = \mathcal{G} \times_1 U^{(1)} \times_2 U^{(2)} \times_3 \cdots \times_k U^{(k)},$$

    其中 $\mathcal{G} \in \mathbb{R}^{r_1 \times r_2 \times \cdots \times r_k}$ 是核心张量，$U^{(i)} \in \mathbb{R}^{n_i \times r_i}$ 是因子矩阵，$\times_i$ 表示第 $i$ 模的矩阵乘积。

    在分量形式下：

    $$\mathcal{T}_{i_1 i_2 \cdots i_k} = \sum_{j_1=1}^{r_1} \cdots \sum_{j_k=1}^{r_k} \mathcal{G}_{j_1 j_2 \cdots j_k} \, U^{(1)}_{i_1 j_1} U^{(2)}_{i_2 j_2} \cdots U^{(k)}_{i_k j_k}.$$

!!! example "例 21.8"
    考虑一个 $3 \times 4 \times 2$ 的三阶张量 $\mathcal{T}$，它的 CP 分解为

    $$\mathcal{T} = \mathbf{a}_1 \otimes \mathbf{b}_1 \otimes \mathbf{c}_1 + \mathbf{a}_2 \otimes \mathbf{b}_2 \otimes \mathbf{c}_2,$$

    其中 $\mathbf{a}_r \in \mathbb{R}^3$，$\mathbf{b}_r \in \mathbb{R}^4$，$\mathbf{c}_r \in \mathbb{R}^2$。分量形式为

    $$\mathcal{T}_{ijk} = a_{1i} b_{1j} c_{1k} + a_{2i} b_{2j} c_{2k}, \quad 1 \leq i \leq 3, \ 1 \leq j \leq 4, \ 1 \leq k \leq 2.$$

    这是一个秩为 $2$ 的张量。Tucker 分解则可以取 $U^{(1)} \in \mathbb{R}^{3 \times 2}$，$U^{(2)} \in \mathbb{R}^{4 \times 2}$，$U^{(3)} \in \mathbb{R}^{2 \times 2}$，核心张量 $\mathcal{G} \in \mathbb{R}^{2 \times 2 \times 2}$。

## 21.8 张量的应用

### 多维数据分析

张量方法在数据科学中的应用越来越广泛。多维数据（如用户-物品-时间推荐数据、彩色图像序列、脑电信号等）自然地用张量来表示。

!!! example "例 21.9"
    **推荐系统中的张量分解**：设有 $m$ 个用户、$n$ 个物品和 $p$ 个时间点。用户 $i$ 在时间 $k$ 对物品 $j$ 的评分构成一个三阶张量 $\mathcal{T} \in \mathbb{R}^{m \times n \times p}$。CP 分解将其表示为

    $$\mathcal{T} \approx \sum_{r=1}^{R} \mathbf{u}_r \otimes \mathbf{v}_r \otimes \mathbf{w}_r,$$

    其中 $\mathbf{u}_r \in \mathbb{R}^m$ 编码用户特征，$\mathbf{v}_r \in \mathbb{R}^n$ 编码物品特征，$\mathbf{w}_r \in \mathbb{R}^p$ 编码时间模式。低秩近似自然地进行数据补全和降维。

### 量子信息中的张量

!!! example "例 21.10"
    **量子纠缠态**。在量子计算中，$n$ 个量子比特（qubit）的状态空间是 $(\mathbb{C}^2)^{\otimes n}$，即 $n$ 个 $\mathbb{C}^2$ 的张量积，维数为 $2^n$。一般的量子态 $|\psi\rangle \in (\mathbb{C}^2)^{\otimes n}$ 可以写为

    $$|\psi\rangle = \sum_{i_1, \ldots, i_n \in \{0,1\}} c_{i_1 \cdots i_n} |i_1\rangle \otimes \cdots \otimes |i_n\rangle.$$

    当 $|\psi\rangle$ 不能写成 $|\phi_1\rangle \otimes \cdots \otimes |\phi_n\rangle$ 的形式时，称该态是**纠缠的**（entangled）。例如，Bell 态

    $$|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$$

    是纠缠态。判定纠缠与张量秩密切相关。

!!! proposition "命题 21.2 (纠缠与张量秩)"
    纯态 $|\psi\rangle \in \mathcal{H}_A \otimes \mathcal{H}_B$ 是可分的（separable）当且仅当它作为张量的秩为 $1$，即可以写成 $|\psi\rangle = |\phi_A\rangle \otimes |\phi_B\rangle$。纠缠态对应张量秩大于 $1$ 的情形。

??? proof "证明"
    这是张量积定义的直接推论。$|\psi\rangle$ 可分意味着 $|\psi\rangle = |\phi_A\rangle \otimes |\phi_B\rangle$，这恰好是一个简单（秩一）张量。反之，若 $|\psi\rangle$ 的张量秩为 $1$，则由定义它可以写成一个简单张量，即可分态。
