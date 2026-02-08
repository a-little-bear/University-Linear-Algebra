# 第 6 章 特征值与特征向量

<div class="context-flow" markdown>

**前置**：第 3 章行列式 · 第 5 章相似矩阵 · 不变子空间

**本章脉络**：$A\mathbf{v}=\lambda\mathbf{v}$ → 特征多项式 $\det(A-\lambda I)=0$ → 代数/几何重数 → 特征空间 → 对角化条件 → 实对称矩阵谱定理 → Cayley-Hamilton → 相似不变量

**延伸**：特征值在 Google PageRank（Perron 向量）、量子力学（能量本征态）、振动分析（固有频率）、人口模型（Leslie 矩阵）中有直接应用；无穷维推广为谱理论，连续谱与点谱的区分是量子力学的核心问题

</div>

特征值（eigenvalue）与特征向量（eigenvector）是线性代数中最深刻、应用最广泛的概念之一。直观地说，线性变换的特征向量是那些"方向不变"的向量——变换只改变了它们的长度（伸缩）。这一简单的思想引出了丰富的理论：特征多项式、对角化、谱定理等。特征值理论不仅在纯数学中具有核心地位，在量子力学、振动分析、主成分分析、Google PageRank 算法等领域也有着广泛的应用。本章将系统地研究特征值与特征向量的定义、计算、性质和应用。

---

## 6.1 特征值与特征向量的定义

<div class="context-flow" markdown>

**从不变子空间到特征值**：第 5 章一维不变子空间 $\operatorname{span}\{\mathbf{v}\}$ → $T(\mathbf{v})=\lambda\mathbf{v}$ → $\lambda=0$ 时 $\mathbf{v} \in \ker A$（$A$ 不可逆），$\mathbf{v} \neq \mathbf{0}$ 是本质要求

</div>

!!! definition "定义 6.1 (特征值与特征向量)"
    设 $A$ 为 $n$ 阶方阵（或 $T: V \to V$ 为线性算子）。若存在标量 $\lambda \in \mathbb{F}$ 和**非零**向量 $\mathbf{v} \in \mathbb{F}^n$（或 $\mathbf{v} \in V$，$\mathbf{v} \neq \mathbf{0}$），使得

    $$A\mathbf{v} = \lambda\mathbf{v}$$

    则称 $\lambda$ 为 $A$ 的一个**特征值**（eigenvalue），$\mathbf{v}$ 为 $A$ 的属于特征值 $\lambda$ 的一个**特征向量**（eigenvector）。

!!! note "注"
    在 $A\mathbf{v} = \lambda\mathbf{v}$ 中，要求 $\mathbf{v} \neq \mathbf{0}$ 是本质的——零向量对任何 $\lambda$ 都满足该等式，但不提供任何信息。而 $\lambda = 0$ 是允许的；此时 $A\mathbf{v} = \mathbf{0}$，即 $\mathbf{v} \in \ker(A)$。因此 $0$ 是特征值当且仅当 $A$ 不可逆。

!!! proposition "命题 6.1 (特征值的几何意义)"
    设 $\lambda$ 为 $A$ 的特征值，$\mathbf{v}$ 为对应的特征向量，则：

    - 若 $\lambda > 1$：$\mathbf{v}$ 方向被"拉伸"。
    - 若 $0 < \lambda < 1$：$\mathbf{v}$ 方向被"压缩"。
    - 若 $\lambda = 1$：$\mathbf{v}$ 方向不变（$\mathbf{v}$ 是不动点方向）。
    - 若 $\lambda = 0$：$\mathbf{v}$ 方向被映射到零。
    - 若 $\lambda < 0$：$\mathbf{v}$ 方向被反转并缩放。

!!! example "例 6.1"
    设 $A = \begin{pmatrix} 3 & 1 \\ 0 & 2 \end{pmatrix}$。取 $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$，则

    $$A\mathbf{v}_1 = \begin{pmatrix} 3 \\ 0 \end{pmatrix} = 3\begin{pmatrix} 1 \\ 0 \end{pmatrix} = 3\mathbf{v}_1$$

    故 $\lambda_1 = 3$ 是特征值，$\mathbf{v}_1$ 是对应的特征向量。

    取 $\mathbf{v}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$，则

    $$A\mathbf{v}_2 = \begin{pmatrix} 2 \\ -2 \end{pmatrix} = 2\begin{pmatrix} 1 \\ -1 \end{pmatrix} = 2\mathbf{v}_2$$

    故 $\lambda_2 = 2$ 是特征值，$\mathbf{v}_2$ 是对应的特征向量。

---

## 6.2 特征多项式

<div class="context-flow" markdown>

**将特征值问题转化为求根**：$A\mathbf{v}=\lambda\mathbf{v}$ → $(A-\lambda I)\mathbf{v}=\mathbf{0}$ 有非零解 → $\det(A-\lambda I)=0$（第 3 章） → $n$ 次多项式有 $n$ 个复根

</div>

!!! definition "定义 6.2 (特征多项式)"
    设 $A$ 为 $n$ 阶方阵，则

    $$p(\lambda) = \det(A - \lambda I)$$

    称为 $A$ 的**特征多项式**（characteristic polynomial）。$A$ 的特征值恰好是特征多项式的根。

!!! note "注"
    部分教材定义特征多项式为 $\det(\lambda I - A)$，二者仅差一个 $(-1)^n$ 的符号因子，根相同。本书采用 $\det(A - \lambda I)$ 的约定。

!!! theorem "定理 6.1 (特征值与特征多项式)"
    $\lambda_0$ 是 $A$ 的特征值当且仅当 $\det(A - \lambda_0 I) = 0$。

??? proof "证明"
    $\lambda_0$ 是特征值 $\Leftrightarrow$ 存在非零 $\mathbf{v}$ 使得 $A\mathbf{v} = \lambda_0 \mathbf{v}$ $\Leftrightarrow$ $(A - \lambda_0 I)\mathbf{v} = \mathbf{0}$ 有非零解 $\Leftrightarrow$ $A - \lambda_0 I$ 不可逆 $\Leftrightarrow$ $\det(A - \lambda_0 I) = 0$。 $\blacksquare$

!!! theorem "定理 6.2 (特征多项式的次数与首项)"
    设 $A$ 为 $n$ 阶方阵，则 $p(\lambda) = \det(A - \lambda I)$ 是 $\lambda$ 的 $n$ 次多项式：

    $$p(\lambda) = (-1)^n \lambda^n + (-1)^{n-1}(\operatorname{tr}A)\lambda^{n-1} + \cdots + \det A$$

    其中 $\operatorname{tr}A = a_{11} + a_{22} + \cdots + a_{nn}$ 为 $A$ 的**迹**（trace）。

??? proof "证明"
    行列式 $\det(A - \lambda I)$ 展开后，$\lambda^n$ 项来自主对角线的乘积：

    $$(a_{11} - \lambda)(a_{22} - \lambda)\cdots(a_{nn} - \lambda) = (-\lambda)^n + (-\lambda)^{n-1}(a_{11} + \cdots + a_{nn}) + \cdots$$

    其他项至多包含 $n-2$ 个对角元素因子，因此对 $\lambda$ 的次数贡献至多为 $n-2$。故最高次项为 $(-1)^n \lambda^n$，次高次项系数为 $(-1)^{n-1}\operatorname{tr}A$。

    常数项 $p(0) = \det(A - 0 \cdot I) = \det A$。 $\blacksquare$

!!! example "例 6.2"
    求 $A = \begin{pmatrix} 4 & -2 \\ 1 & 1 \end{pmatrix}$ 的特征值。

    $$p(\lambda) = \det(A - \lambda I) = \det\begin{pmatrix} 4 - \lambda & -2 \\ 1 & 1 - \lambda \end{pmatrix} = (4-\lambda)(1-\lambda) + 2 = \lambda^2 - 5\lambda + 6 = (\lambda - 2)(\lambda - 3)$$

    特征值为 $\lambda_1 = 2$，$\lambda_2 = 3$。

    验证：$\operatorname{tr}A = 4 + 1 = 5 = 2 + 3$，$\det A = 4 \cdot 1 - (-2) \cdot 1 = 6 = 2 \times 3$。

!!! example "例 6.3"
    求 $A = \begin{pmatrix} 2 & 1 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 3 \end{pmatrix}$ 的特征值。

    $$p(\lambda) = \det(A - \lambda I) = (2-\lambda)(2-\lambda)(3-\lambda) = (2-\lambda)^2(3-\lambda)$$

    特征值为 $\lambda_1 = 2$（二重根）和 $\lambda_2 = 3$（一重根）。

---

## 6.3 特征值的性质

<div class="context-flow" markdown>

**代数与几何的张力**：$1 \le \operatorname{GM}(\lambda) \le \operatorname{AM}(\lambda)$ → 二者相等时可对角化 → $\operatorname{tr}A = \sum\lambda_i$，$\det A = \prod\lambda_i$（连接第 3 章行列式与特征值）

</div>

!!! definition "定义 6.3 (代数重数)"
    特征值 $\lambda_0$ 作为特征多项式 $p(\lambda)$ 的根的重数，称为 $\lambda_0$ 的**代数重数**（algebraic multiplicity），记为 $\operatorname{AM}(\lambda_0)$。

!!! definition "定义 6.4 (几何重数)"
    特征值 $\lambda_0$ 对应的特征空间 $E_{\lambda_0} = \ker(A - \lambda_0 I)$ 的维数，称为 $\lambda_0$ 的**几何重数**（geometric multiplicity），记为 $\operatorname{GM}(\lambda_0)$。

!!! theorem "定理 6.3 (迹与行列式)"
    设 $A$ 为 $n$ 阶方阵，特征值（含重数）为 $\lambda_1, \lambda_2, \ldots, \lambda_n$，则

    $$\operatorname{tr}A = \sum_{i=1}^n \lambda_i, \qquad \det A = \prod_{i=1}^n \lambda_i$$

??? proof "证明"
    特征多项式可以写成

    $$p(\lambda) = (-1)^n(\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n)$$

    展开后，$\lambda^{n-1}$ 的系数为 $(-1)^n \cdot (-1)(\lambda_1 + \cdots + \lambda_n) = (-1)^{n-1}(\lambda_1 + \cdots + \lambda_n)$。由定理 6.2，此系数也等于 $(-1)^{n-1}\operatorname{tr}A$，故 $\operatorname{tr}A = \lambda_1 + \cdots + \lambda_n$。

    令 $\lambda = 0$：$p(0) = (-1)^n(-\lambda_1)(-\lambda_2)\cdots(-\lambda_n) = \lambda_1\lambda_2\cdots\lambda_n$。又 $p(0) = \det A$。 $\blacksquare$

!!! theorem "定理 6.4 (几何重数与代数重数的关系)"
    设 $\lambda_0$ 为 $A$ 的特征值，则

    $$1 \leq \operatorname{GM}(\lambda_0) \leq \operatorname{AM}(\lambda_0)$$

??? proof "证明"
    $\operatorname{GM}(\lambda_0) \geq 1$ 是因为特征值必有非零特征向量，故特征空间维数至少为 $1$。

    设 $\operatorname{GM}(\lambda_0) = k$，取特征空间 $E_{\lambda_0}$ 的一组基 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$，扩充为 $\mathbb{F}^n$ 的基 $\{\mathbf{v}_1, \ldots, \mathbf{v}_k, \mathbf{v}_{k+1}, \ldots, \mathbf{v}_n\}$。

    设 $P = (\mathbf{v}_1 \;\; \cdots \;\; \mathbf{v}_n)$，则

    $$P^{-1}AP = \begin{pmatrix} \lambda_0 I_k & B \\ 0 & C \end{pmatrix}$$

    其中 $I_k$ 为 $k$ 阶单位矩阵。相似矩阵有相同的特征多项式：

    $$p(\lambda) = \det(P^{-1}AP - \lambda I) = (\lambda_0 - \lambda)^k \det(C - \lambda I)$$

    因此 $(\lambda_0 - \lambda)^k$ 整除 $p(\lambda)$，即 $\operatorname{AM}(\lambda_0) \geq k = \operatorname{GM}(\lambda_0)$。 $\blacksquare$

!!! theorem "定理 6.5 (不同特征值的特征向量线性无关)"
    设 $\lambda_1, \lambda_2, \ldots, \lambda_k$ 为 $A$ 的互不相同的特征值，$\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k$ 分别为对应的特征向量，则 $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\}$ 线性无关。

??? proof "证明"
    对 $k$ 进行归纳。$k = 1$ 时，$\mathbf{v}_1 \neq \mathbf{0}$ 显然线性无关。

    假设 $k - 1$ 时成立，考虑 $k$ 的情形。设

    $$c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_k\mathbf{v}_k = \mathbf{0} \quad \cdots (*)$$

    对 $(*)$ 两边施加 $A$：

    $$c_1\lambda_1\mathbf{v}_1 + c_2\lambda_2\mathbf{v}_2 + \cdots + c_k\lambda_k\mathbf{v}_k = \mathbf{0} \quad \cdots (**)$$

    将 $(**)$ 减去 $\lambda_k$ 倍的 $(*)$：

    $$c_1(\lambda_1 - \lambda_k)\mathbf{v}_1 + c_2(\lambda_2 - \lambda_k)\mathbf{v}_2 + \cdots + c_{k-1}(\lambda_{k-1} - \lambda_k)\mathbf{v}_{k-1} = \mathbf{0}$$

    由归纳假设，$\{\mathbf{v}_1, \ldots, \mathbf{v}_{k-1}\}$ 线性无关，且 $\lambda_i - \lambda_k \neq 0$（$i < k$），故 $c_i(\lambda_i - \lambda_k) = 0$，即 $c_i = 0$（$i = 1, \ldots, k-1$）。回代入 $(*)$ 得 $c_k\mathbf{v}_k = \mathbf{0}$，由 $\mathbf{v}_k \neq \mathbf{0}$ 得 $c_k = 0$。 $\blacksquare$

!!! example "例 6.4"
    设 $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$。特征多项式为 $p(\lambda) = \lambda^2 + 1$。

    在 $\mathbb{R}$ 上无实根，因此 $A$ 没有实特征值。但在 $\mathbb{C}$ 上，特征值为 $\lambda_1 = i$，$\lambda_2 = -i$。

    验证：$\operatorname{tr}A = 0 = i + (-i)$，$\det A = 1 = i \cdot (-i)$。

    这说明实矩阵可能没有实特征值，但在复数域上，$n$ 阶矩阵必有 $n$ 个特征值（含重数，代数基本定理）。

---

## 6.4 特征空间

<div class="context-flow" markdown>

**特征空间 = $\ker(A-\lambda I)$**：第 1 章齐次方程组的解空间 → 其维数 = 几何重数 → 不同特征值的特征空间**直和**（第 4 章）

</div>

!!! definition "定义 6.5 (特征空间)"
    设 $\lambda_0$ 为 $A$ 的特征值，则

    $$E_{\lambda_0} = \ker(A - \lambda_0 I) = \{\mathbf{v} \in \mathbb{F}^n : A\mathbf{v} = \lambda_0\mathbf{v}\}$$

    称为 $A$ 的属于 $\lambda_0$ 的**特征空间**（eigenspace）。

!!! theorem "定理 6.6 (特征空间是子空间)"
    $E_{\lambda_0}$ 是 $\mathbb{F}^n$ 的子空间。

??? proof "证明"
    $E_{\lambda_0} = \ker(A - \lambda_0 I)$ 是齐次线性方程组的解空间，由第 4 章的结论，它是 $\mathbb{F}^n$ 的子空间。或直接验证：$\mathbf{0} \in E_{\lambda_0}$（虽然 $\mathbf{0}$ 不是特征向量，但在特征空间中）；若 $\mathbf{u}, \mathbf{v} \in E_{\lambda_0}$，则 $A(\mathbf{u} + \mathbf{v}) = A\mathbf{u} + A\mathbf{v} = \lambda_0\mathbf{u} + \lambda_0\mathbf{v} = \lambda_0(\mathbf{u} + \mathbf{v})$；若 $c \in \mathbb{F}$，则 $A(c\mathbf{v}) = cA\mathbf{v} = c\lambda_0\mathbf{v} = \lambda_0(c\mathbf{v})$。 $\blacksquare$

!!! example "例 6.5"
    设 $A = \begin{pmatrix} 5 & -2 & 0 \\ 0 & 3 & 0 \\ 4 & -2 & 3 \end{pmatrix}$。求特征值与特征空间。

    **特征多项式**：

    $$p(\lambda) = \det(A - \lambda I) = \det\begin{pmatrix} 5-\lambda & -2 & 0 \\ 0 & 3-\lambda & 0 \\ 4 & -2 & 3-\lambda \end{pmatrix}$$

    按第三列展开：$p(\lambda) = (3 - \lambda)\det\begin{pmatrix} 5-\lambda & -2 \\ 0 & 3-\lambda \end{pmatrix} = (3-\lambda)^2(5-\lambda)$。

    特征值为 $\lambda_1 = 3$（$\operatorname{AM} = 2$），$\lambda_2 = 5$（$\operatorname{AM} = 1$）。

    **$\lambda_1 = 3$ 的特征空间**：

    $$A - 3I = \begin{pmatrix} 2 & -2 & 0 \\ 0 & 0 & 0 \\ 4 & -2 & 0 \end{pmatrix} \xrightarrow{\text{行化简}} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$$

    等等，让我们仔细计算。$R_3 - 2R_1$：$\begin{pmatrix} 2 & -2 & 0 \\ 0 & 0 & 0 \\ 0 & 2 & 0 \end{pmatrix}$，然后 $R_1/2$：$\begin{pmatrix} 1 & -1 & 0 \\ 0 & 0 & 0 \\ 0 & 2 & 0 \end{pmatrix}$，$R_3/2$：$\begin{pmatrix} 1 & -1 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 0 \end{pmatrix}$，交换 $R_2, R_3$ 后：$\begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$，最终 $\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$。

    故 $x_1 = x_2 = 0$，$x_3 = t$ 自由。$E_3 = \operatorname{span}\left\{\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}\right\}$，$\operatorname{GM}(3) = 1 < 2 = \operatorname{AM}(3)$。

    **$\lambda_2 = 5$ 的特征空间**：

    $$A - 5I = \begin{pmatrix} 0 & -2 & 0 \\ 0 & -2 & 0 \\ 4 & -2 & -2 \end{pmatrix} \xrightarrow{\text{行化简}} \begin{pmatrix} 1 & 0 & -\frac{1}{2} \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$$

    故 $x_1 = \frac{1}{2}t$，$x_2 = 0$，$x_3 = t$。$E_5 = \operatorname{span}\left\{\begin{pmatrix} 1 \\ 0 \\ 2 \end{pmatrix}\right\}$，$\operatorname{GM}(5) = 1 = \operatorname{AM}(5)$。

---

## 6.5 对角化

<div class="context-flow" markdown>

**线性代数的核心目标之一**：$A = PDP^{-1}$ → $P$ 的列是特征向量，$D$ 的对角元是特征值 → 可对角化 $\Leftrightarrow$ $n$ 个线性无关特征向量 $\Leftrightarrow$ $\forall i:\operatorname{GM}=\operatorname{AM}$

</div>

!!! definition "定义 6.6 (可对角化)"
    $n$ 阶方阵 $A$ 称为**可对角化的**（diagonalizable），若存在可逆矩阵 $P$ 和对角矩阵 $D$ 使得

    $$A = PDP^{-1}$$

    等价地，$D = P^{-1}AP$，即 $A$ 相似于对角矩阵。

!!! theorem "定理 6.7 (对角化的充要条件)"
    $n$ 阶方阵 $A$ 可对角化当且仅当以下等价条件之一成立：

    1. $A$ 有 $n$ 个线性无关的特征向量。
    2. 对 $A$ 的每个特征值 $\lambda_i$，$\operatorname{GM}(\lambda_i) = \operatorname{AM}(\lambda_i)$。
    3. 特征空间的直和等于全空间：$\bigoplus_i E_{\lambda_i} = \mathbb{F}^n$。

??? proof "证明"
    **$A$ 可对角化 $\Leftrightarrow$ 条件 1**：

    $(\Rightarrow)$ 若 $A = PDP^{-1}$，设 $P = (\mathbf{p}_1 \;\; \cdots \;\; \mathbf{p}_n)$，$D = \operatorname{diag}(d_1, \ldots, d_n)$。则 $AP = PD$，即 $A\mathbf{p}_j = d_j \mathbf{p}_j$。因此 $\mathbf{p}_1, \ldots, \mathbf{p}_n$ 是 $A$ 的特征向量，且由 $P$ 可逆知它们线性无关。

    $(\Leftarrow)$ 设 $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ 是 $n$ 个线性无关的特征向量，$A\mathbf{v}_j = \lambda_j\mathbf{v}_j$。令 $P = (\mathbf{v}_1 \;\; \cdots \;\; \mathbf{v}_n)$，$D = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，则 $AP = PD$，即 $A = PDP^{-1}$。

    **条件 1 $\Leftrightarrow$ 条件 2**：$A$ 有 $n$ 个线性无关的特征向量 $\Leftrightarrow$ $\sum_i \operatorname{GM}(\lambda_i) = n = \sum_i \operatorname{AM}(\lambda_i)$。由 $\operatorname{GM}(\lambda_i) \leq \operatorname{AM}(\lambda_i)$ 和总和相等，必有 $\operatorname{GM}(\lambda_i) = \operatorname{AM}(\lambda_i)$ 对每个 $i$ 成立。 $\blacksquare$

!!! corollary "推论 6.1"
    若 $n$ 阶方阵 $A$ 有 $n$ 个互不相同的特征值，则 $A$ 可对角化。

??? proof "证明"
    由定理 6.5，不同特征值的特征向量线性无关，$n$ 个不同特征值给出 $n$ 个线性无关的特征向量。由定理 6.7 条件 1，$A$ 可对角化。 $\blacksquare$

!!! example "例 6.6"
    对角化 $A = \begin{pmatrix} 4 & -2 \\ 1 & 1 \end{pmatrix}$。

    由例 6.2，$\lambda_1 = 2$，$\lambda_2 = 3$。

    $\lambda_1 = 2$：$(A - 2I)\mathbf{v} = \mathbf{0}$，$\begin{pmatrix} 2 & -2 \\ 1 & -1 \end{pmatrix}\mathbf{v} = \mathbf{0}$，得 $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。

    $\lambda_2 = 3$：$(A - 3I)\mathbf{v} = \mathbf{0}$，$\begin{pmatrix} 1 & -2 \\ 1 & -2 \end{pmatrix}\mathbf{v} = \mathbf{0}$，得 $\mathbf{v}_2 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$。

    令 $P = \begin{pmatrix} 1 & 2 \\ 1 & 1 \end{pmatrix}$，$D = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$，则 $A = PDP^{-1}$。

    **应用**：$A^k = PD^kP^{-1} = P\begin{pmatrix} 2^k & 0 \\ 0 & 3^k \end{pmatrix}P^{-1}$，可以快速计算矩阵的高次幂。

!!! example "例 6.7"
    $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ 不可对角化。

    $p(\lambda) = (2 - \lambda)^2$，唯一特征值 $\lambda = 2$，$\operatorname{AM}(2) = 2$。

    $A - 2I = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$\operatorname{rank}(A - 2I) = 1$，$\operatorname{GM}(2) = 2 - 1 = 1 < 2 = \operatorname{AM}(2)$。

    $A$ 没有两个线性无关的特征向量，因此不可对角化。

---

## 6.6 实对称矩阵的特征值

<div class="context-flow" markdown>

**对称矩阵的完美性质**：特征值全实 → 不同特征值的特征向量正交 → **谱定理**：$A=QDQ^T$（正交对角化） → 连接第 7 章正交矩阵

</div>

!!! definition "定义 6.7 (实对称矩阵)"
    实方阵 $A$ 称为**对称的**（symmetric），若 $A^T = A$。

!!! theorem "定理 6.8 (实对称矩阵的特征值为实数)"
    设 $A$ 为实对称矩阵，则 $A$ 的所有特征值都是实数。

??? proof "证明"
    设 $\lambda$ 为 $A$ 的特征值（可能为复数），$\mathbf{v} \in \mathbb{C}^n$ 为对应的特征向量。则 $A\mathbf{v} = \lambda\mathbf{v}$。

    取共轭转置：$\overline{\mathbf{v}}^T A^T = \overline{\lambda} \overline{\mathbf{v}}^T$（因 $A$ 为实矩阵，$\overline{A} = A$）。又 $A^T = A$，故 $\overline{\mathbf{v}}^T A = \overline{\lambda} \overline{\mathbf{v}}^T$。

    计算 $\overline{\mathbf{v}}^T A \mathbf{v}$：

    - 从左边：$\overline{\mathbf{v}}^T (A\mathbf{v}) = \overline{\mathbf{v}}^T (\lambda\mathbf{v}) = \lambda (\overline{\mathbf{v}}^T \mathbf{v})$
    - 从右边：$(\overline{\mathbf{v}}^T A)\mathbf{v} = \overline{\lambda} (\overline{\mathbf{v}}^T \mathbf{v})$

    因此 $\lambda (\overline{\mathbf{v}}^T \mathbf{v}) = \overline{\lambda} (\overline{\mathbf{v}}^T \mathbf{v})$。又 $\overline{\mathbf{v}}^T \mathbf{v} = \sum |v_i|^2 > 0$（因 $\mathbf{v} \neq \mathbf{0}$），故 $\lambda = \overline{\lambda}$，即 $\lambda \in \mathbb{R}$。 $\blacksquare$

!!! theorem "定理 6.9 (实对称矩阵不同特征值的特征向量正交)"
    设 $A$ 为实对称矩阵，$\lambda_1 \neq \lambda_2$ 为 $A$ 的两个不同特征值，$\mathbf{v}_1, \mathbf{v}_2$ 分别为对应的特征向量，则 $\mathbf{v}_1 \perp \mathbf{v}_2$（即 $\mathbf{v}_1^T \mathbf{v}_2 = 0$）。

??? proof "证明"
    $\lambda_1(\mathbf{v}_1^T\mathbf{v}_2) = (\lambda_1\mathbf{v}_1)^T\mathbf{v}_2 = (A\mathbf{v}_1)^T\mathbf{v}_2 = \mathbf{v}_1^T A^T \mathbf{v}_2 = \mathbf{v}_1^T A \mathbf{v}_2 = \mathbf{v}_1^T(\lambda_2\mathbf{v}_2) = \lambda_2(\mathbf{v}_1^T\mathbf{v}_2)$

    因此 $(\lambda_1 - \lambda_2)(\mathbf{v}_1^T\mathbf{v}_2) = 0$。由 $\lambda_1 \neq \lambda_2$，得 $\mathbf{v}_1^T\mathbf{v}_2 = 0$。 $\blacksquare$

!!! theorem "定理 6.10 (谱定理 / 正交对角化)"
    实对称矩阵 $A$ 必可正交对角化，即存在正交矩阵 $Q$（$Q^TQ = I$）和对角矩阵 $D$ 使得

    $$A = QDQ^T$$

??? proof "证明"
    （概要）由代数基本定理，$A$ 的特征多项式在 $\mathbb{C}$ 上有 $n$ 个根。由定理 6.8，这些根都是实数。

    对 $n$ 进行归纳。$n = 1$ 时显然。设 $n - 1$ 时成立。

    取 $A$ 的一个特征值 $\lambda_1$ 和对应的单位特征向量 $\mathbf{q}_1$。将 $\mathbf{q}_1$ 扩充为 $\mathbb{R}^n$ 的标准正交基 $\{\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n\}$。令 $Q_1 = (\mathbf{q}_1 \;\; \cdots \;\; \mathbf{q}_n)$，则

    $$Q_1^T A Q_1 = \begin{pmatrix} \lambda_1 & \mathbf{b}^T \\ \mathbf{0} & A_1 \end{pmatrix}$$

    由 $A$ 对称，$Q_1^T A Q_1$ 也对称，故 $\mathbf{b} = \mathbf{0}$ 且 $A_1$ 是 $(n-1)$ 阶实对称矩阵。由归纳假设，$A_1$ 可正交对角化。由此可构造 $A$ 的正交对角化。 $\blacksquare$

!!! example "例 6.8"
    正交对角化 $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。

    $p(\lambda) = (2-\lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda - 1)(\lambda - 3)$。

    $\lambda_1 = 1$：$\mathbf{v}_1 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$，单位化 $\mathbf{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix}$。

    $\lambda_2 = 3$：$\mathbf{v}_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$，单位化 $\mathbf{q}_2 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}$。

    验证 $\mathbf{q}_1^T\mathbf{q}_2 = 0$。

    $$Q = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}, \quad D = \begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}, \quad A = QDQ^T$$

---

## 6.7 Cayley-Hamilton 定理

<div class="context-flow" markdown>

**矩阵满足自己的特征多项式**：$p(A)=O$ → 推论：$A^{-1}$ 可表示为 $A$ 的多项式 → 暗示矩阵代数中 $A$ 的幂次是有限维的

</div>

!!! theorem "定理 6.11 (Cayley-Hamilton 定理)"
    设 $A$ 为 $n$ 阶方阵，$p(\lambda) = \det(A - \lambda I)$ 为 $A$ 的特征多项式，则

    $$p(A) = O$$

    即将 $\lambda$ 替换为 $A$，$p(\lambda)$ 中的常数项乘以 $I$，所得矩阵为零矩阵。

!!! note "注"
    Cayley-Hamilton 定理并**不**是简单地将 $\lambda$ 代入 $\det(A - \lambda I) = 0$ 然后令 $\lambda = A$。行列式是一个标量值函数，而这里需要对**矩阵多项式**求值。定理的正确含义是：若 $p(\lambda) = c_0 + c_1\lambda + \cdots + c_n\lambda^n$，则 $c_0 I + c_1 A + \cdots + c_n A^n = O$。

??? proof "证明"
    考虑矩阵 $B(\lambda) = \operatorname{adj}(A - \lambda I)$（$A - \lambda I$ 的伴随矩阵）。$B(\lambda)$ 的每个元素是 $\lambda$ 的至多 $n-1$ 次多项式，因此可以写成

    $$B(\lambda) = B_0 + B_1\lambda + B_2\lambda^2 + \cdots + B_{n-1}\lambda^{n-1}$$

    其中 $B_0, B_1, \ldots, B_{n-1}$ 为常数矩阵。

    由伴随矩阵的性质：$(A - \lambda I)B(\lambda) = \det(A - \lambda I) \cdot I = p(\lambda)I$。

    设 $p(\lambda) = c_0 + c_1\lambda + \cdots + c_n\lambda^n$。将左边展开并比较 $\lambda$ 的各次幂系数：

    $$AB_0 = c_0 I$$
    $$AB_1 - B_0 = c_1 I$$
    $$AB_2 - B_1 = c_2 I$$
    $$\vdots$$
    $$AB_{n-1} - B_{n-2} = c_{n-1} I$$
    $$-B_{n-1} = c_n I$$

    分别将上述等式左乘 $I, A, A^2, \ldots, A^{n-1}, A^n$ 并求和：

    $$AB_0 + A(AB_1 - B_0) + A^2(AB_2 - B_1) + \cdots + A^n(-B_{n-1})$$
    $$= c_0 I + c_1 A + c_2 A^2 + \cdots + c_n A^n = p(A)$$

    左边展开后，所有含 $B_i$ 的项逐对抵消（伸缩求和），结果为 $O$。因此 $p(A) = O$。 $\blacksquare$

!!! example "例 6.9"
    验证 Cayley-Hamilton 定理：$A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$。

    $p(\lambda) = \det(A - \lambda I) = (1-\lambda)(4-\lambda) - 6 = \lambda^2 - 5\lambda - 2$。

    $$A^2 = \begin{pmatrix} 7 & 10 \\ 15 & 22 \end{pmatrix}$$

    $$p(A) = A^2 - 5A - 2I = \begin{pmatrix} 7 & 10 \\ 15 & 22 \end{pmatrix} - \begin{pmatrix} 5 & 10 \\ 15 & 20 \end{pmatrix} - \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$$

    验证通过。

!!! corollary "推论 6.2"
    设 $A$ 为 $n$ 阶可逆方阵，特征多项式为 $p(\lambda) = c_0 + c_1\lambda + \cdots + c_n\lambda^n$。由 $A$ 可逆知 $c_0 = p(0) = \det A \neq 0$。由 Cayley-Hamilton 定理：

    $$c_0 I + c_1 A + \cdots + c_n A^n = O$$

    两边乘以 $A^{-1}$：

    $$A^{-1} = -\frac{1}{c_0}(c_1 I + c_2 A + \cdots + c_n A^{n-1})$$

    即 $A^{-1}$ 可以表示为 $A$ 的多项式。

!!! example "例 6.10"
    利用推论 6.2 求 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ 的逆。

    $p(\lambda) = \lambda^2 - 5\lambda - 2$，$c_0 = -2$，$c_1 = -5$，$c_2 = 1$。

    $$A^{-1} = -\frac{1}{-2}(-5I + A) = \frac{1}{2}\begin{pmatrix} -4 & 2 \\ 3 & -1 \end{pmatrix} = \begin{pmatrix} -2 & 1 \\ 3/2 & -1/2 \end{pmatrix}$$

    验证：$AA^{-1} = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} -2 & 1 \\ 3/2 & -1/2 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$。正确。

---

## 6.8 相似矩阵的性质

<div class="context-flow" markdown>

**相似 = 同一线性算子的不同表示**（第 5 章） → 共享：特征多项式、迹、行列式、秩 → 但特征值相同不一定相似（需考虑 Jordan 结构）

</div>

!!! definition "定义 6.8 (相似不变量)"
    若性质或量在相似变换 $A \mapsto P^{-1}AP$ 下不变，则称其为**相似不变量**（similarity invariant）。

!!! theorem "定理 6.12 (相似矩阵的基本性质)"
    若 $A \sim B$（$B = P^{-1}AP$），则：

    1. $\det A = \det B$
    2. $\operatorname{tr} A = \operatorname{tr} B$
    3. $\operatorname{rank} A = \operatorname{rank} B$
    4. $A$ 与 $B$ 有相同的特征多项式（从而有相同的特征值，含重数）
    5. $A$ 可逆当且仅当 $B$ 可逆
    6. $A^k \sim B^k$，对所有正整数 $k$

??? proof "证明"
    **1.** $\det B = \det(P^{-1}AP) = \det(P^{-1})\det(A)\det(P) = \frac{1}{\det P} \cdot \det A \cdot \det P = \det A$。

    **2.** 利用迹的循环性质 $\operatorname{tr}(XY) = \operatorname{tr}(YX)$：$\operatorname{tr}B = \operatorname{tr}(P^{-1}AP) = \operatorname{tr}(APP^{-1}) = \operatorname{tr}A$。

    **3.** $P^{-1}$ 和 $P$ 都可逆，相似变换不改变秩。

    **4.** $\det(B - \lambda I) = \det(P^{-1}AP - \lambda P^{-1}IP) = \det(P^{-1}(A - \lambda I)P) = \det(A - \lambda I)$。

    **5.** 由 1，$\det A = \det B$，故 $A$ 可逆 $\Leftrightarrow$ $\det A \neq 0$ $\Leftrightarrow$ $\det B \neq 0$ $\Leftrightarrow$ $B$ 可逆。

    **6.** $B^k = (P^{-1}AP)^k = P^{-1}A^k P$。 $\blacksquare$

!!! note "注"
    相似是比相等更弱的关系——相似矩阵共享许多重要性质，但不一定相等。相似不变量包括但不限于：特征值、特征多项式、迹、行列式、秩、最小多项式、Jordan 标准形等。然而，**特征值相同并不一定意味着矩阵相似**。例如 $\begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$ 和 $\begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ 有相同的特征值但不相似。

!!! example "例 6.11"
    设 $A = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 3 & -4 \\ 1 & -1 \end{pmatrix}$。判断 $A$ 与 $B$ 是否相似。

    $p_A(\lambda) = (1 - \lambda)(3 - \lambda) = \lambda^2 - 4\lambda + 3$。

    $p_B(\lambda) = (3 - \lambda)(-1 - \lambda) + 4 = \lambda^2 - 2\lambda + 1 = (\lambda - 1)^2$。

    $A$ 的特征值为 $1, 3$；$B$ 的特征值为 $1, 1$。特征多项式不同，故 $A \not\sim B$。

!!! example "例 6.12"
    设 $A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$，$B = \begin{pmatrix} 2 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$。

    $A$ 和 $B$ 都是对角矩阵，特征多项式均为 $(1 - \lambda)^2(2 - \lambda)$，相同。

    它们是否相似？由于对角矩阵相似当且仅当对角元素是同一组数的排列（无需考虑顺序），$A$ 和 $B$ 确实相似。实际上取置换矩阵 $P = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix}$（交换第 1 和第 3 行/列），则 $P^{-1}AP = B$。

---

## 本章小结

本章系统研究了特征值与特征向量理论，主要内容包括：

1. **特征值与特征向量的定义**：$A\mathbf{v} = \lambda\mathbf{v}$ 揭示了线性变换的"不变方向"。
2. **特征多项式**：通过 $\det(A - \lambda I) = 0$ 求特征值，将代数问题转化为多项式求根。
3. **特征值的性质**：迹等于特征值之和，行列式等于特征值之积；代数重数不小于几何重数。
4. **特征空间**：属于同一特征值的所有特征向量（加上零向量）构成子空间。
5. **对角化**：矩阵可对角化的充要条件是每个特征值的几何重数等于代数重数。
6. **实对称矩阵**：特征值必为实数，不同特征值的特征向量正交，可正交对角化（谱定理）。
7. **Cayley-Hamilton 定理**：每个方阵满足其自身的特征多项式。
8. **相似不变量**：相似矩阵共享特征多项式、迹、行列式、秩等性质。
