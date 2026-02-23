# 第 4 章 向量空间

<div class="context-flow" markdown>

**前置**：第 1 章解集结构 · 第 2 章 $\mathbb{R}^{m \times n}$ 的运算律

**本章脉络**：8 条公理 → 子空间 → 线性组合/张成 → 线性无关 → 基与维数 → 坐标/过渡矩阵 → 四个基本子空间 → 秩-零化度定理

**延伸**：向量空间概念推广到无穷维：Hilbert 空间（量子力学的数学框架）、Banach 空间（泛函分析）；拓扑向量空间结合了代数结构与拓扑结构；模论将向量空间从域上推广到环上

</div>

向量空间（vector space）是线性代数的核心抽象概念。在前几章中，我们主要在 $\mathbb{R}^n$ 的具体框架下讨论问题。本章将引入向量空间的公理化定义，使得线性代数的理论适用于远比 $\mathbb{R}^n$ 更广泛的对象——矩阵、多项式、函数等。我们将系统讨论子空间、线性无关、基与维数等基本概念，建立秩-零化度定理，并介绍四个基本子空间。

---

## 4.1 向量空间的定义与公理

<div class="context-flow" markdown>

**公理化跃迁**：第 2 章 $\mathbb{R}^{m \times n}$ 满足的 8 条运算律 → 抽象为任意集合上的公理 → "向量"不再局限于箭头或列向量

</div>

!!! definition "定义 4.1 (向量空间)"
    设 $\mathbb{F}$ 为一个数域（通常取 $\mathbb{R}$ 或 $\mathbb{C}$），$V$ 为一个非空集合，在 $V$ 上定义了**加法**运算 $+: V \times V \to V$ 和**标量乘法**运算 $\cdot: \mathbb{F} \times V \to V$。若以下 8 条公理对所有 $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ 和所有 $c, d \in \mathbb{F}$ 成立，则称 $V$ 为 $\mathbb{F}$ 上的一个**向量空间**（vector space），$V$ 的元素称为**向量**（vector），$\mathbb{F}$ 的元素称为**标量**（scalar）：

    **加法公理**：

    1. $\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}$（交换律）
    2. $(\mathbf{u} + \mathbf{v}) + \mathbf{w} = \mathbf{u} + (\mathbf{v} + \mathbf{w})$（结合律）
    3. 存在**零向量** $\mathbf{0} \in V$，使得 $\mathbf{v} + \mathbf{0} = \mathbf{v}$（零元素）
    4. 对每个 $\mathbf{v} \in V$，存在 $-\mathbf{v} \in V$，使得 $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$（加法逆元）

    **标量乘法公理**：

    5. $c(\mathbf{u} + \mathbf{v}) = c\mathbf{u} + c\mathbf{v}$（对向量加法的分配律）
    6. $(c + d)\mathbf{v} = c\mathbf{v} + d\mathbf{v}$（对标量加法的分配律）
    7. $c(d\mathbf{v}) = (cd)\mathbf{v}$（标量乘法的结合律）
    8. $1\mathbf{v} = \mathbf{v}$（标量乘法的单位元）

!!! proposition "命题 4.1 (向量空间的基本推论)"
    在向量空间 $V$ 中：

    1. 零向量唯一。
    2. 每个向量的加法逆元唯一。
    3. $0\mathbf{v} = \mathbf{0}$（标量 $0$ 乘以任何向量得零向量）。
    4. $c\mathbf{0} = \mathbf{0}$（任何标量乘以零向量得零向量）。
    5. $(-1)\mathbf{v} = -\mathbf{v}$。

??? proof "证明"
    3. $0\mathbf{v} = (0+0)\mathbf{v} = 0\mathbf{v} + 0\mathbf{v}$。两边加上 $0\mathbf{v}$ 的加法逆元，得 $\mathbf{0} = 0\mathbf{v}$。

    5. $\mathbf{v} + (-1)\mathbf{v} = 1\mathbf{v} + (-1)\mathbf{v} = (1 + (-1))\mathbf{v} = 0\mathbf{v} = \mathbf{0}$。因此 $(-1)\mathbf{v}$ 是 $\mathbf{v}$ 的加法逆元，由唯一性得 $(-1)\mathbf{v} = -\mathbf{v}$。$\blacksquare$

---

## 4.2 常见向量空间举例

<div class="context-flow" markdown>

**公理的力量**：$\mathbb{R}^n$、矩阵空间 $\mathbb{R}^{m \times n}$、多项式空间 $\mathbb{P}_n$、函数空间 $C(\mathbb{R})$ → 形态各异但共享同一代数结构

</div>

!!! example "例 4.1"
    **$\mathbb{R}^n$ 空间**：$n$ 维实列向量的集合 $\mathbb{R}^n$，按分量进行加法和标量乘法，是最基本的向量空间。

!!! example "例 4.2"
    **矩阵空间 $\mathbb{R}^{m \times n}$**：全体 $m \times n$ 实矩阵在矩阵加法和标量乘法下构成向量空间。

!!! example "例 4.3"
    **多项式空间 $\mathbb{P}_n$**：次数不超过 $n$ 的所有实系数多项式构成向量空间，其中加法和标量乘法定义为多项式的通常加法和数乘。零向量为零多项式。

    全体实系数多项式的集合记为 $\mathbb{P}$，同样构成向量空间（无限维）。

!!! example "例 4.4"
    **函数空间 $\mathcal{F}(\mathbb{R}, \mathbb{R})$**：从 $\mathbb{R}$ 到 $\mathbb{R}$ 的所有函数的集合，按点态加法 $(f+g)(x) = f(x) + g(x)$ 和标量乘法 $(cf)(x) = cf(x)$ 构成向量空间。连续函数空间 $C(\mathbb{R})$、可微函数空间 $C^1(\mathbb{R})$ 等是其子空间。

!!! example "例 4.5"
    **零空间 $\{0\}$**：仅含零向量的集合构成向量空间，称为**零空间**或**平凡向量空间**，维数为 $0$。

---

## 4.3 子空间

<div class="context-flow" markdown>

**"空间中的空间"**：对加法和数乘封闭即为子空间 → 第 1 章齐次解集是子空间的原型 → 必须过原点（包含 $\mathbf{0}$）

</div>

!!! definition "定义 4.2 (子空间)"
    设 $V$ 为向量空间，$W \subseteq V$ 为非空子集。若 $W$ 在 $V$ 的加法和标量乘法下本身也构成向量空间，则称 $W$ 为 $V$ 的一个**子空间**（subspace）。

!!! theorem "定理 4.1 (子空间判定定理)"
    设 $V$ 为向量空间，$W \subseteq V$ 为非空子集，则 $W$ 是 $V$ 的子空间当且仅当满足以下两个条件：

    1. **对加法封闭**：若 $\mathbf{u}, \mathbf{v} \in W$，则 $\mathbf{u} + \mathbf{v} \in W$。
    2. **对标量乘法封闭**：若 $\mathbf{v} \in W$，$c \in \mathbb{F}$，则 $c\mathbf{v} \in W$。

    等价地，$W$ 是子空间当且仅当：$W$ 非空，且对任意 $\mathbf{u}, \mathbf{v} \in W$，$c, d \in \mathbb{F}$，有 $c\mathbf{u} + d\mathbf{v} \in W$（对线性组合封闭）。

??? proof "证明"
    $(\Rightarrow)$ 显然，向量空间对加法和标量乘法封闭。

    $(\Leftarrow)$ 需要验证 8 条公理。交换律、结合律、分配律等继承自 $V$。由条件 2 取 $c = 0$ 得 $\mathbf{0} = 0\mathbf{v} \in W$（零元素存在）。取 $c = -1$ 得 $-\mathbf{v} \in W$（加法逆元存在）。$\blacksquare$

!!! example "例 4.6"
    以下是 $\mathbb{R}^3$ 的子空间：

    - $\{\mathbf{0}\}$（零子空间）。
    - 过原点的直线 $\{t\mathbf{v} : t \in \mathbb{R}\}$，其中 $\mathbf{v} \neq \mathbf{0}$。
    - 过原点的平面 $\{s\mathbf{u} + t\mathbf{v} : s, t \in \mathbb{R}\}$，其中 $\mathbf{u}, \mathbf{v}$ 线性无关。
    - $\mathbb{R}^3$ 本身。

    不过原点的直线或平面**不是**子空间（不包含零向量）。

---

## 4.4 线性组合与张成空间

<div class="context-flow" markdown>

**生成子空间**：给定向量组 → 所有线性组合构成的集合 = 包含这些向量的**最小子空间** → $A\mathbf{x}=\mathbf{b}$ 有解 $\Leftrightarrow$ $\mathbf{b} \in \operatorname{Span}\{\text{列向量}\}$

</div>

!!! definition "定义 4.3 (线性组合)"
    设 $V$ 为向量空间，$\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k \in V$。形如

    $$
    c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_k\mathbf{v}_k \quad (c_i \in \mathbb{F})
    $$

    的向量称为 $\mathbf{v}_1, \ldots, \mathbf{v}_k$ 的一个**线性组合**（linear combination）。

!!! definition "定义 4.4 (张成空间)"
    向量组 $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\}$ 的所有线性组合构成的集合称为该向量组的**张成空间**（span），记为

    $$
    \operatorname{Span}\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\} = \left\{\sum_{i=1}^k c_i\mathbf{v}_i : c_1, \ldots, c_k \in \mathbb{F}\right\}.
    $$

    若 $\operatorname{Span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\} = V$，则称 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ **张成** $V$（或**生成** $V$）。

!!! theorem "定理 4.2 (张成空间是子空间)"
    $\operatorname{Span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ 是 $V$ 的子空间，而且是包含 $\mathbf{v}_1, \ldots, \mathbf{v}_k$ 的最小子空间。

??? proof "证明"
    设 $W = \operatorname{Span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$。取所有系数为零得零向量，故 $W$ 非空。

    设 $\mathbf{u} = \sum a_i\mathbf{v}_i \in W$，$\mathbf{w} = \sum b_i\mathbf{v}_i \in W$，$c, d \in \mathbb{F}$，则

    $$
    c\mathbf{u} + d\mathbf{w} = \sum(ca_i + db_i)\mathbf{v}_i \in W.
    $$

    故 $W$ 对线性组合封闭，是子空间。任何包含 $\mathbf{v}_1, \ldots, \mathbf{v}_k$ 的子空间必然包含它们的所有线性组合，因此 $W$ 是最小的。$\blacksquare$

---

## 4.5 线性无关与线性相关

<div class="context-flow" markdown>

**无冗余 = 线性无关**：$\sum c_i \mathbf{v}_i = \mathbf{0}$ 只有平凡解 → 没有向量是其余的线性组合 → 在 $\mathbb{R}^n$ 中等价于矩阵列满秩（第 2 章）

</div>

!!! definition "定义 4.5 (线性无关与线性相关)"
    向量组 $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\}$ 称为**线性无关的**（linearly independent），若方程

    $$
    c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_k\mathbf{v}_k = \mathbf{0}
    $$

    只有平凡解 $c_1 = c_2 = \cdots = c_k = 0$。

    否则，称向量组**线性相关**（linearly dependent），即存在不全为零的标量 $c_1, \ldots, c_k$ 使得 $\sum c_i\mathbf{v}_i = \mathbf{0}$。

!!! theorem "定理 4.3 (线性相关的等价条件)"
    向量组 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$（$k \ge 2$）线性相关当且仅当其中至少有一个向量可以表示为其余向量的线性组合。

??? proof "证明"
    $(\Rightarrow)$ 设 $\sum c_i\mathbf{v}_i = \mathbf{0}$，其中某个 $c_j \neq 0$。则 $\mathbf{v}_j = -\frac{1}{c_j}\sum_{i \neq j} c_i\mathbf{v}_i$。

    $(\Leftarrow)$ 设 $\mathbf{v}_j = \sum_{i \neq j} a_i\mathbf{v}_i$，则 $\sum_{i \neq j} a_i\mathbf{v}_i - \mathbf{v}_j = \mathbf{0}$，其中 $\mathbf{v}_j$ 的系数为 $-1 \neq 0$。$\blacksquare$

!!! proposition "命题 4.2 (线性无关的性质)"
    1. 单个非零向量 $\{\mathbf{v}\}$（$\mathbf{v} \neq \mathbf{0}$）是线性无关的。
    2. 含有零向量的向量组必线性相关。
    3. 线性无关向量组的任何子集仍线性无关。
    4. 若 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ 线性无关且 $\mathbf{v}_{k+1} \notin \operatorname{Span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$，则 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{v}_{k+1}\}$ 也线性无关。

!!! example "例 4.7"
    在 $\mathbb{R}^3$ 中，判断向量组 $\{(1,0,1)^T,\; (0,1,1)^T,\; (1,1,0)^T\}$ 是否线性无关。

    构造矩阵并求行阶梯形：

    $$
    \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 1 & 1 & 0 \end{pmatrix}
    \xrightarrow{R_3 - R_1}
    \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 0 & 1 & -1 \end{pmatrix}
    \xrightarrow{R_3 - R_2}
    \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & -2 \end{pmatrix}
    $$

    有 $3$ 个主元，齐次方程只有平凡解，因此三个向量线性无关。

---

## 4.6 基与维数

<div class="context-flow" markdown>

**基 = 线性无关 + 张成**：每个向量唯一表示为基向量的线性组合 → 所有基大小相同（维数的良定义性） → $\dim(V)$ 是向量空间的终极不变量

</div>

!!! definition "定义 4.6 (基)"
    设 $V$ 为向量空间，向量组 $\mathcal{B} = \{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n\}$ 称为 $V$ 的一组**基**（basis），若满足：

    1. $\mathcal{B}$ 线性无关。
    2. $\operatorname{Span}(\mathcal{B}) = V$（$\mathcal{B}$ 张成 $V$）。

!!! example "例 4.8"
    $\mathbb{R}^n$ 的**标准基**（standard basis）为 $\{\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\}$，其中 $\mathbf{e}_i$ 的第 $i$ 个分量为 $1$、其余为 $0$。

    $\mathbb{P}_2$ 的标准基为 $\{1, x, x^2\}$。

    $\mathbb{R}^{2 \times 2}$ 的标准基为 $\left\{\begin{pmatrix}1&0\\0&0\end{pmatrix}, \begin{pmatrix}0&1\\0&0\end{pmatrix}, \begin{pmatrix}0&0\\1&0\end{pmatrix}, \begin{pmatrix}0&0\\0&1\end{pmatrix}\right\}$。

!!! theorem "定理 4.4 (基的唯一表示性)"
    $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 是向量空间 $V$ 的基当且仅当 $V$ 中每个向量都能唯一地表示为 $\mathcal{B}$ 中向量的线性组合。

??? proof "证明"
    $(\Rightarrow)$ **存在性**由张成性保证。**唯一性**：设 $\mathbf{v} = \sum a_i\mathbf{v}_i = \sum b_i\mathbf{v}_i$，则 $\sum(a_i - b_i)\mathbf{v}_i = \mathbf{0}$。由线性无关性，$a_i - b_i = 0$，即 $a_i = b_i$。

    $(\Leftarrow)$ 唯一表示蕴含张成性（每个向量至少有一种表示）和线性无关性（$\mathbf{0}$ 的唯一表示为全零系数）。$\blacksquare$

!!! theorem "定理 4.5 (维数的良定义性)"
    有限维向量空间 $V$ 的任意两组基包含相同数量的向量。

??? proof "证明"
    设 $\mathcal{B}_1 = \{\mathbf{u}_1, \ldots, \mathbf{u}_m\}$ 和 $\mathcal{B}_2 = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 都是 $V$ 的基。

    由于 $\mathcal{B}_1$ 张成 $V$ 且 $\mathcal{B}_2$ 线性无关，由 Steinitz 替换引理可得 $n \le m$。

    对称地，由于 $\mathcal{B}_2$ 张成 $V$ 且 $\mathcal{B}_1$ 线性无关，可得 $m \le n$。

    因此 $m = n$。$\blacksquare$

!!! definition "定义 4.7 (维数)"
    有限维向量空间 $V$ 的基中包含的向量个数称为 $V$ 的**维数**（dimension），记为 $\dim(V)$。零空间 $\{\mathbf{0}\}$ 的维数定义为 $0$。

!!! theorem "定理 4.6 (基的判定)"
    设 $\dim(V) = n$，则：

    1. $V$ 中任何超过 $n$ 个向量的集合必线性相关。
    2. $V$ 中任何少于 $n$ 个向量的集合不能张成 $V$。
    3. $V$ 中 $n$ 个线性无关的向量构成 $V$ 的基。
    4. $V$ 中张成 $V$ 的 $n$ 个向量构成 $V$ 的基。

---

## 4.7 坐标与坐标变换

<div class="context-flow" markdown>

**抽象到具体的桥梁**：选定基 → 每个向量对应唯一坐标向量 $[\mathbf{x}]_\mathcal{B} \in \mathbb{R}^n$ → 换基 = 过渡矩阵 $P$ → 第 5 章基变换 $\to$ 相似矩阵

</div>

!!! definition "定义 4.8 (坐标向量)"
    设 $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 为向量空间 $V$ 的一组有序基，$\mathbf{x} \in V$ 可唯一表示为

    $$
    \mathbf{x} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_n\mathbf{v}_n.
    $$

    则列向量 $[\mathbf{x}]_\mathcal{B} = (c_1, c_2, \ldots, c_n)^T \in \mathbb{R}^n$ 称为 $\mathbf{x}$ 关于基 $\mathcal{B}$ 的**坐标向量**（coordinate vector）。

!!! definition "定义 4.9 (过渡矩阵)"
    设 $\mathcal{B} = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 和 $\mathcal{B}' = \{\mathbf{w}_1, \ldots, \mathbf{w}_n\}$ 为 $V$ 的两组基。将 $\mathcal{B}'$ 中的每个向量用 $\mathcal{B}$ 来表示：

    $$
    \mathbf{w}_j = \sum_{i=1}^n p_{ij}\mathbf{v}_i, \quad j = 1, \ldots, n.
    $$

    矩阵 $P = (p_{ij})$ 称为从基 $\mathcal{B}'$ 到基 $\mathcal{B}$ 的**过渡矩阵**（transition matrix / change-of-basis matrix），满足

    $$
    [\mathbf{x}]_\mathcal{B} = P [\mathbf{x}]_{\mathcal{B}'}.
    $$

!!! theorem "定理 4.7 (过渡矩阵可逆)"
    过渡矩阵 $P$ 是可逆的，且 $P^{-1}$ 是从 $\mathcal{B}$ 到 $\mathcal{B}'$ 的过渡矩阵。

!!! example "例 4.9"
    在 $\mathbb{R}^2$ 中，设 $\mathcal{B} = \{\mathbf{e}_1, \mathbf{e}_2\}$（标准基），$\mathcal{B}' = \{(1,1)^T, (1,-1)^T\}$。

    则 $\mathbf{w}_1 = 1\cdot\mathbf{e}_1 + 1\cdot\mathbf{e}_2$，$\mathbf{w}_2 = 1\cdot\mathbf{e}_1 + (-1)\cdot\mathbf{e}_2$，过渡矩阵为

    $$
    P = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}.
    $$

    对于向量 $\mathbf{x} = (3, 1)^T$，其在 $\mathcal{B}'$ 下的坐标为 $[\mathbf{x}]_{\mathcal{B}'} = P^{-1}[\mathbf{x}]_\mathcal{B} = \frac{1}{-2}\begin{pmatrix} -1 & -1 \\ -1 & 1 \end{pmatrix}\begin{pmatrix} 3 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$。

    验证：$2(1,1)^T + 1(1,-1)^T = (3, 1)^T$。

---

## 4.8 行空间、列空间与零空间

<div class="context-flow" markdown>

**矩阵的四个基本子空间**：$\operatorname{Col}(A)$, $\operatorname{Row}(A)$, $\operatorname{Null}(A)$, $\operatorname{Null}(A^T)$ → 维数由 $\operatorname{rank}(A) = r$ 统一控制 → 正交补关系见第 7 章

</div>

!!! definition "定义 4.10 (四个基本子空间)"
    设 $A$ 为 $m \times n$ 矩阵，定义以下四个子空间：

    1. **列空间**（column space）：$\operatorname{Col}(A) = \{A\mathbf{x} : \mathbf{x} \in \mathbb{R}^n\} \subseteq \mathbb{R}^m$，即 $A$ 的列向量的张成空间。
    2. **行空间**（row space）：$\operatorname{Row}(A) = \operatorname{Col}(A^T) \subseteq \mathbb{R}^n$，即 $A$ 的行向量的张成空间。
    3. **零空间/核**（null space / kernel）：$\operatorname{Null}(A) = \{\mathbf{x} \in \mathbb{R}^n : A\mathbf{x} = \mathbf{0}\} \subseteq \mathbb{R}^n$。
    4. **左零空间**（left null space）：$\operatorname{Null}(A^T) = \{\mathbf{y} \in \mathbb{R}^m : A^T\mathbf{y} = \mathbf{0}\} \subseteq \mathbb{R}^m$。

!!! theorem "定理 4.8 (四个子空间的维数)"
    设 $A$ 为 $m \times n$ 矩阵，$\operatorname{rank}(A) = r$，则：

    1. $\dim(\operatorname{Col}(A)) = r$
    2. $\dim(\operatorname{Row}(A)) = r$
    3. $\dim(\operatorname{Null}(A)) = n - r$
    4. $\dim(\operatorname{Null}(A^T)) = m - r$

!!! theorem "定理 4.9 (正交补关系)"
    1. $\operatorname{Row}(A) \oplus \operatorname{Null}(A) = \mathbb{R}^n$，即行空间与零空间互为正交补。
    2. $\operatorname{Col}(A) \oplus \operatorname{Null}(A^T) = \mathbb{R}^m$，即列空间与左零空间互为正交补。

!!! example "例 4.10"
    设 $A = \begin{pmatrix} 1 & 2 & 1 & 3 \\ 2 & 4 & 3 & 7 \\ 1 & 2 & 2 & 4 \end{pmatrix}$。

    化为简化行阶梯形：

    $$
    R = \begin{pmatrix} 1 & 2 & 0 & 2 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 0 \end{pmatrix}
    $$

    $\operatorname{rank}(A) = 2$。

    - **列空间** $\operatorname{Col}(A)$：主元列为第 1、3 列，基为 $\{(1,2,1)^T,\; (1,3,2)^T\}$，维数 $2$。
    - **行空间** $\operatorname{Row}(A)$：$R$ 的非零行为行空间的基：$\{(1,2,0,2),\; (0,0,1,1)\}$，维数 $2$。
    - **零空间** $\operatorname{Null}(A)$：自由变量 $x_2 = s$，$x_4 = t$，基为 $\{(-2,1,0,0)^T,\; (-2,0,-1,1)^T\}$，维数 $2$。
    - **左零空间** $\operatorname{Null}(A^T)$：维数 $3 - 2 = 1$。

---

## 4.9 秩-零化度定理（维数公式）

<div class="context-flow" markdown>

**维数守恒**：$\operatorname{rank}(A) + \operatorname{nullity}(A) = n$ → 线性映射"保留的维度"与"丢失的维度"之和 = 定义域维度 → 第 5 章秩-零化度定理的矩阵版本

</div>

!!! theorem "定理 4.10 (秩-零化度定理)"
    设 $A$ 为 $m \times n$ 矩阵，则

    $$
    \operatorname{rank}(A) + \operatorname{nullity}(A) = n,
    $$

    其中 $\operatorname{nullity}(A) = \dim(\operatorname{Null}(A))$ 称为 $A$ 的**零化度**（nullity）。

    更一般地，若 $T: V \to W$ 为有限维向量空间之间的线性变换，则

    $$
    \dim(\ker T) + \dim(\operatorname{im} T) = \dim(V).
    $$

??? proof "证明"
    设 $\dim(V) = n$，$\dim(\ker T) = k$。取 $\ker T$ 的基 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$，将其扩展为 $V$ 的基 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{v}_{k+1}, \ldots, \mathbf{v}_n\}$。

    **断言**：$\{T(\mathbf{v}_{k+1}), \ldots, T(\mathbf{v}_n)\}$ 是 $\operatorname{im}(T)$ 的基。

    **张成性**：对任意 $\mathbf{w} \in \operatorname{im}(T)$，存在 $\mathbf{v} = \sum_{i=1}^n c_i\mathbf{v}_i \in V$ 使得 $T(\mathbf{v}) = \mathbf{w}$。则

    $$
    \mathbf{w} = T\left(\sum_{i=1}^n c_i\mathbf{v}_i\right) = \sum_{i=1}^k c_i T(\mathbf{v}_i) + \sum_{i=k+1}^n c_i T(\mathbf{v}_i) = \sum_{i=k+1}^n c_i T(\mathbf{v}_i),
    $$

    因为 $T(\mathbf{v}_i) = \mathbf{0}$（$i = 1, \ldots, k$）。

    **线性无关性**：设 $\sum_{i=k+1}^n c_i T(\mathbf{v}_i) = \mathbf{0}$，则 $T\left(\sum_{i=k+1}^n c_i\mathbf{v}_i\right) = \mathbf{0}$，即 $\sum_{i=k+1}^n c_i\mathbf{v}_i \in \ker T$。因此 $\sum_{i=k+1}^n c_i\mathbf{v}_i = \sum_{i=1}^k d_i\mathbf{v}_i$，即 $\sum_{i=1}^k d_i\mathbf{v}_i - \sum_{i=k+1}^n c_i\mathbf{v}_i = \mathbf{0}$。由 $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 的线性无关性，所有系数为零，特别地 $c_{k+1} = \cdots = c_n = 0$。

    因此 $\dim(\operatorname{im} T) = n - k = \dim(V) - \dim(\ker T)$。$\blacksquare$

!!! example "例 4.11"
    设 $A$ 为 $5 \times 8$ 矩阵，$\operatorname{rank}(A) = 3$。则 $\operatorname{nullity}(A) = 8 - 3 = 5$。方程组 $A\mathbf{x} = \mathbf{0}$ 的解空间是 $\mathbb{R}^8$ 的 $5$ 维子空间，其基础解系包含 $5$ 个向量。

---

## 4.10 子空间的和与直和

<div class="context-flow" markdown>

**子空间的代数运算**：$W_1 + W_2$ 是包含两者的最小子空间 → $W_1 \cap W_2 = \{\mathbf{0}\}$ 时为**直和** → 维数公式：$\dim(W_1+W_2) = \dim W_1 + \dim W_2 - \dim(W_1 \cap W_2)$

</div>

!!! definition "定义 4.11 (子空间的和)"
    设 $W_1, W_2$ 为向量空间 $V$ 的子空间，定义它们的**和**（sum）为

    $$
    W_1 + W_2 = \{\mathbf{w}_1 + \mathbf{w}_2 : \mathbf{w}_1 \in W_1, \mathbf{w}_2 \in W_2\}.
    $$

    $W_1 + W_2$ 也是 $V$ 的子空间，且是包含 $W_1 \cup W_2$ 的最小子空间。

!!! definition "定义 4.12 (直和)"
    若 $W_1 \cap W_2 = \{\mathbf{0}\}$，则称和 $W_1 + W_2$ 为**直和**（direct sum），记为 $W_1 \oplus W_2$。此时 $W_1 + W_2$ 中的每个向量可以**唯一**地分解为 $\mathbf{w}_1 + \mathbf{w}_2$（$\mathbf{w}_i \in W_i$）。

!!! theorem "定理 4.11 (直和的等价条件)"
    设 $W_1, W_2$ 为 $V$ 的子空间，以下条件等价：

    1. $W_1 + W_2$ 是直和。
    2. $W_1 \cap W_2 = \{\mathbf{0}\}$。
    3. $W_1 + W_2$ 中每个向量的分解 $\mathbf{w}_1 + \mathbf{w}_2$（$\mathbf{w}_i \in W_i$）唯一。
    4. 若 $\mathbf{w}_1 + \mathbf{w}_2 = \mathbf{0}$（$\mathbf{w}_i \in W_i$），则 $\mathbf{w}_1 = \mathbf{w}_2 = \mathbf{0}$。

??? proof "证明"
    $(2) \Rightarrow (3)$：设 $\mathbf{v} = \mathbf{w}_1 + \mathbf{w}_2 = \mathbf{w}_1' + \mathbf{w}_2'$。则 $\mathbf{w}_1 - \mathbf{w}_1' = \mathbf{w}_2' - \mathbf{w}_2 \in W_1 \cap W_2 = \{\mathbf{0}\}$，故 $\mathbf{w}_1 = \mathbf{w}_1'$，$\mathbf{w}_2 = \mathbf{w}_2'$。

    $(3) \Rightarrow (4)$：$\mathbf{0} = \mathbf{0} + \mathbf{0}$ 是一种分解，由唯一性即得。

    $(4) \Rightarrow (2)$：设 $\mathbf{v} \in W_1 \cap W_2$，则 $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$（$\mathbf{v} \in W_1$，$-\mathbf{v} \in W_2$），由条件 4 得 $\mathbf{v} = \mathbf{0}$。$\blacksquare$

!!! theorem "定理 4.12 (维数公式)"
    设 $W_1, W_2$ 为有限维向量空间 $V$ 的子空间，则

    $$
    \dim(W_1 + W_2) = \dim(W_1) + \dim(W_2) - \dim(W_1 \cap W_2).
    $$

    特别地，若 $W_1 + W_2$ 是直和，则 $\dim(W_1 \oplus W_2) = \dim(W_1) + \dim(W_2)$。

??? proof "证明"
    设 $\dim(W_1 \cap W_2) = k$，取 $W_1 \cap W_2$ 的基 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$。

    将此基扩展为 $W_1$ 的基 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{u}_1, \ldots, \mathbf{u}_p\}$，其中 $p = \dim(W_1) - k$。

    将此基扩展为 $W_2$ 的基 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{w}_1, \ldots, \mathbf{w}_q\}$，其中 $q = \dim(W_2) - k$。

    可以验证 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{u}_1, \ldots, \mathbf{u}_p, \mathbf{w}_1, \ldots, \mathbf{w}_q\}$ 是 $W_1 + W_2$ 的基。因此

    $$
    \dim(W_1 + W_2) = k + p + q = \dim(W_1) + \dim(W_2) - \dim(W_1 \cap W_2). \quad \blacksquare
    $$

!!! example "例 4.12"
    在 $\mathbb{R}^4$ 中，设 $W_1 = \operatorname{Span}\{(1,0,1,0)^T, (0,1,0,1)^T\}$，$W_2 = \operatorname{Span}\{(1,1,0,0)^T, (0,0,1,1)^T\}$。

    $\dim(W_1) = \dim(W_2) = 2$。为判断是否直和，求 $W_1 \cap W_2$：设 $a(1,0,1,0)^T + b(0,1,0,1)^T = c(1,1,0,0)^T + d(0,0,1,1)^T$，得方程组

    $$
    a = c, \quad b = c, \quad a = d, \quad b = d.
    $$

    因此 $a = b = c = d$，$W_1 \cap W_2 = \operatorname{Span}\{(1,1,1,1)^T\}$，$\dim(W_1 \cap W_2) = 1$。

    $\dim(W_1 + W_2) = 2 + 2 - 1 = 3$，不是直和。

## 练习题

1. **[概念] 为什么零向量必须包含在任何向量子空间中？**
   ??? success "参考答案"
       因为子空间必须对标量乘法封闭。对任意 $\mathbf{v} \in W$，取标量 $0$，有 $0 \cdot \mathbf{v} = \mathbf{0} \in W$。

2. **[基础] 判断所有满足 $x + y = 1$ 的向量 $(x,y)^T$ 是否构成 $\mathbb{R}^2$ 的子空间。**
   ??? success "参考答案"
       否，它不包含零向量 $(0,0)^T$（因为 $0 + 0 \neq 1$）。

3. **[直和] 设 $V = U \oplus W$。证明零向量有且仅有唯一的分解 $\mathbf{0} = \mathbf{u} + \mathbf{w}$，其中 $\mathbf{u} \in U, \mathbf{w} \in W$。**
   ??? success "参考答案"
       显然 $\mathbf{u}=\mathbf{0}, \mathbf{w}=\mathbf{0}$ 是一个分解。若存在其他分解，则 $\mathbf{u} = -\mathbf{w} \in U \cap W$。因为是直和，$U \cap W = \{\mathbf{0}\}$，故唯一。

4. **[维数] 设 $\mathbb{P}_n$ 是所有次数不超过 $n$ 的实系数多项式构成的空间。它的维数是多少？**
   ??? success "参考答案"
       $n+1$。它的标准基为 $\{1, x, x^2, \dots, x^n\}$，共有 $n+1$ 个元素。

5. **[无穷维] 举出一个无穷维向量空间的例子。**
   ??? success "参考答案"
       实数轴上所有连续函数的集合 $C(\mathbb{R})$。或者全体多项式的空间 $\mathbb{P}$。

6. **[域] $\mathbb{C}$ 作为一个复向量空间，其维数是多少？作为一个实向量空间，其维数又是多少？**
   ??? success "参考答案"
       作为复向量空间维数为 1（基为 $\{1\}$）。作为实向量空间维数为 2（基为 $\{1, i\}$）。

7. **[线性相关] 证明：包含零向量的任何向量组必定是线性相关的。**
   ??? success "参考答案"
       设组为 $\{\mathbf{0}, \mathbf{v}_1, \dots\}$，我们可以取系数 $1 \cdot \mathbf{0} + 0 \cdot \mathbf{v}_1 + \dots = \mathbf{0}$，存在不全为零的系数使得组合为零，故线性相关。

8. **[交与和] 在 $\mathbb{R}^3$ 中，两个不同的二维平面（子空间）的交集维数是多少？**
   ??? success "参考答案"
       由维数公式 $\dim(U+V) = \dim(U) + \dim(V) - \dim(U \cap V)$，即 $3 \ge 2 + 2 - \dim(U \cap V)$，故交集维数至少为 1（通常是一条穿过原点的直线）。

9. **[坐标] 在基 $\mathcal{B} = \{\mathbf{e}_1, \mathbf{e}_2\}$ 下，向量 $\mathbf{v}$ 的坐标是 $(2, -1)^T$。这意味着什么？**
   ??? success "参考答案"
       这意味着 $\mathbf{v}$ 可以唯一的表示为基向量的线性组合：$\mathbf{v} = 2\mathbf{e}_1 - 1\mathbf{e}_2$。

10. **[爱因斯坦思考题] 物理学中的“张量”（Tensor）本质上是不依赖于特定坐标系的几何对象。向量空间中基的变换矩阵如何反映这一相对论精神？**
    ??? success "参考答案"
        向量本身 $\mathbf{v}$ 是客观存在的实体（如速度）。当我们改变观察者的参考系（改变基）时，坐标会发生“逆变”（反向变换），以抵消基向量的“协变”（正向变换）。这种 $P^{-1}$ 与 $P$ 的完美抵消（$\mathbf{v} = P \mathbf{x}_{new} = P(P^{-1}\mathbf{x}_{old}) = \mathbf{x}_{old}$），体现了物理规律在所有参考系下形式不变的核心思想。

## 本章小结

本章建立了线性代数的核心抽象结构——向量空间，主要内容包括：

1. **向量空间的公理化定义**：通过 8 条公理将向量的概念从 $\mathbb{R}^n$ 推广到抽象空间。
2. **子空间与生成集**：子空间是对线性运算封闭的子集；生成集是通过线性组合构造子空间的工具。
3. **线性无关与基**：线性无关去除了生成集中的冗余；基是极大的线性无关集，也是极小的生成集，为空间提供了坐标系。
4. **维数与坐标**：维数是基中向量的个数，是有限维空间分类的唯一不变量。坐标映射给出了抽象空间与 $\mathbb{R}^n$ 之间的同构。
5. **四个基本子空间**：列空间、行空间、零空间和左零空间彻底揭示了矩阵方程 $A\mathbf{x} = \mathbf{b}$ 的解的几何结构。
6. **基变换**：过渡矩阵 $P$ 给出了同一向量在不同坐标系下的坐标联系 $\mathbf{x} = P \mathbf{x}'$。
