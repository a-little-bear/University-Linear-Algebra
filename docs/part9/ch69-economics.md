# 第 69 章 线性代数在经济学中的应用

<div class="context-flow" markdown>

**前置**：线性方程组(Ch1) · 矩阵运算(Ch2) · 特征值(Ch6) · 非负矩阵(Ch17) · M-矩阵(Ch38)

**本章脉络**：Leontief 投入产出模型 $\to$ Hawkins-Simon 条件 $\to$ Sraffa 模型 $\to$ 投入产出分析的结构 $\to$ 线性交换模型 $\to$ 博弈论中的线性代数 $\to$ Von Neumann 极小极大定理 $\to$ 线性规划对偶

**延伸**：投入产出分析（Leontief 诺贝尔经济学奖）是宏观经济政策制定的标准工具；博弈论的矩阵方法是拍卖设计、机制设计和市场均衡分析的基础

</div>

经济学的核心问题之一是理解经济系统中各部门之间的相互依赖关系。Wassily Leontief 于 1930 年代提出的投入产出模型，将一个经济体抽象为若干生产部门之间的线性关系，由此奠定了计量经济学的一个重要分支。这一模型的数学本质完全建立在矩阵代数之上：消费矩阵描述部门间的技术关系，而 $(I - C)^{-1}$ 的存在性与非负性则决定了经济系统的可行性。

与此同时，博弈论——尤其是 Von Neumann 的极小极大定理——同样深深植根于线性代数。二人零和博弈的解可以通过线性规划来计算，而线性规划的对偶理论又为经济变量赋予了深刻的"影子价格"解释。

本章将系统地展示线性代数如何成为经济分析的语言。

---

## 69.1 Leontief 投入产出模型

<div class="context-flow" markdown>

**核心问题**：一个由 $n$ 个部门组成的经济体中，每个部门既是其他部门的供应者，也是其他部门的消费者。给定外部需求，各部门的总产出应该是多少？

</div>

投入产出分析是 Leontief 在哈佛大学期间发展起来的，他用美国经济数据验证了这个理论框架，并因此获得了 1973 年诺贝尔经济学奖。

!!! definition "定义 69.1 (开放 Leontief 模型)"
    设经济体由 $n$ 个部门组成。令 $x_i$ 为第 $i$ 个部门的总产出，$d_i$ 为第 $i$ 个部门面对的外部（最终）需求，$c_{ij}$ 为生产一个单位的第 $j$ 个部门的产品所需要消耗的第 $i$ 个部门的产品数量。定义：

    - **产出向量**：$\mathbf{x} = (x_1, x_2, \ldots, x_n)^T \in \mathbb{R}^n$
    - **需求向量**：$\mathbf{d} = (d_1, d_2, \ldots, d_n)^T \in \mathbb{R}^n$
    - **消费矩阵**（或技术系数矩阵）：$C = (c_{ij}) \in \mathbb{R}^{n \times n}$

    则经济均衡方程为：

    $$\mathbf{x} = C\mathbf{x} + \mathbf{d}$$

    即每个部门的总产出等于其被其他部门消耗的部分加上最终需求。

消费矩阵 $C$ 的列向量 $\mathbf{c}_j$ 描述了生产一个单位第 $j$ 个部门产品所需的各部门投入。因此 $C\mathbf{x}$ 的第 $i$ 个分量 $\sum_j c_{ij} x_j$ 就是第 $i$ 个部门的产品被所有部门的生产活动消耗掉的总量。

!!! theorem "定理 69.1 (Leontief 逆矩阵)"
    若 $(I - C)$ 可逆，则开放 Leontief 模型的解为：

    $$\mathbf{x} = (I - C)^{-1} \mathbf{d}$$

    矩阵 $(I - C)^{-1}$ 称为 **Leontief 逆矩阵**。

??? proof "证明"
    由均衡方程 $\mathbf{x} = C\mathbf{x} + \mathbf{d}$，移项得：

    $$(I - C)\mathbf{x} = \mathbf{d}$$

    若 $(I - C)$ 可逆，则两端左乘 $(I - C)^{-1}$ 即得：

    $$\mathbf{x} = (I - C)^{-1}\mathbf{d}$$

    $\blacksquare$

!!! theorem "定理 69.2 (Neumann 级数展开)"
    若 $C \geq 0$（元素非负）且 $\rho(C) < 1$（谱半径严格小于1），则 $(I - C)^{-1}$ 存在且：

    $$(I - C)^{-1} = \sum_{k=0}^{\infty} C^k = I + C + C^2 + C^3 + \cdots$$

    且 $(I - C)^{-1} \geq 0$。

??? proof "证明"
    因为 $\rho(C) < 1$，级数 $\sum_{k=0}^{\infty} C^k$ 收敛。设 $S_m = \sum_{k=0}^{m} C^k$，则：

    $$(I - C) S_m = I - C^{m+1}$$

    由于 $\rho(C) < 1$，$\|C^{m+1}\| \leq \|C\|^{m+1} \to 0$（对于任意与谱半径相容的矩阵范数），故 $C^{m+1} \to 0$。因此：

    $$(I - C) \lim_{m \to \infty} S_m = I$$

    即 $(I - C)^{-1} = \sum_{k=0}^{\infty} C^k$。

    由于 $C \geq 0$，每个 $C^k \geq 0$，故 $(I - C)^{-1} \geq 0$。

    $\blacksquare$

这个级数展开有深刻的经济含义。$C^k \mathbf{d}$ 表示为满足最终需求 $\mathbf{d}$ 而引发的第 $k$ 轮间接需求（乘数效应）：

- $\mathbf{d}$：第 0 轮，直接的最终需求
- $C\mathbf{d}$：第 1 轮，为生产 $\mathbf{d}$ 所需的中间投入
- $C^2\mathbf{d}$：第 2 轮，为生产第 1 轮投入所需的投入
- 以此类推...

总产出 $\mathbf{x} = (I + C + C^2 + \cdots)\mathbf{d}$ 是所有这些轮次的累积。

!!! example "例 69.1"
    考虑一个两部门经济体（农业和制造业），消费矩阵为：

    $$C = \begin{pmatrix} 0.2 & 0.3 \\ 0.4 & 0.1 \end{pmatrix}$$

    外部需求为 $\mathbf{d} = \begin{pmatrix} 10 \\ 20 \end{pmatrix}$（单位：亿元）。

    **求解**：

    $$I - C = \begin{pmatrix} 0.8 & -0.3 \\ -0.4 & 0.9 \end{pmatrix}$$

    $$\det(I - C) = 0.8 \times 0.9 - (-0.3)(-0.4) = 0.72 - 0.12 = 0.60$$

    $$(I - C)^{-1} = \frac{1}{0.60} \begin{pmatrix} 0.9 & 0.3 \\ 0.4 & 0.8 \end{pmatrix} = \begin{pmatrix} 1.50 & 0.50 \\ 0.667 & 1.333 \end{pmatrix}$$

    因此：

    $$\mathbf{x} = (I - C)^{-1}\mathbf{d} = \begin{pmatrix} 1.50 & 0.50 \\ 0.667 & 1.333 \end{pmatrix}\begin{pmatrix} 10 \\ 20 \end{pmatrix} = \begin{pmatrix} 25.0 \\ 33.33 \end{pmatrix}$$

    农业部门需要总产出 25 亿元，制造业部门需要 33.33 亿元。注意 $(I - C)^{-1}$ 的元素全部为正，这保证了对任意非负需求，产出也是非负的。

!!! example "例 69.2"
    **Leontief 逆矩阵的经济解释**。$(I-C)^{-1}$ 的第 $(i,j)$ 元素表示：当第 $j$ 个部门的最终需求增加 1 个单位时，第 $i$ 个部门的总产出需要增加多少。

    在上例中，$(I-C)^{-1}$ 的 $(1,2)$ 元素为 $0.50$，意味着制造业的最终需求每增加 1 亿元，农业的总产出需要增加 0.50 亿元（因为制造业生产需要农业原料作为投入）。

---

## 69.2 Hawkins-Simon 条件

<div class="context-flow" markdown>

**核心问题**：什么条件下，投入产出模型对任意非负需求 $\mathbf{d} \geq 0$ 都有非负解 $\mathbf{x} \geq 0$？

</div>

!!! definition "定义 69.2 (生产性矩阵)"
    消费矩阵 $C \geq 0$ 称为**生产性的**（productive），若 $(I - C)^{-1}$ 存在且 $(I - C)^{-1} \geq 0$。

!!! theorem "定理 69.3 (Hawkins-Simon 条件)"
    设 $C \geq 0$ 为 $n \times n$ 消费矩阵，$B = I - C$。以下条件等价：

    1. $C$ 是生产性的（即 $(I-C)^{-1}$ 存在且非负）
    2. $B = I - C$ 的所有顺序主子式为正：

    $$b_{11} > 0, \quad \begin{vmatrix} b_{11} & b_{12} \\ b_{21} & b_{22} \end{vmatrix} > 0, \quad \ldots, \quad \det(B) > 0$$

    3. $\rho(C) < 1$（$C$ 的谱半径严格小于 1）
    4. $B = I - C$ 是 M-矩阵
    5. 存在 $\mathbf{x} > 0$ 使得 $B\mathbf{x} > 0$

??? proof "证明"
    我们证明 $(1) \Leftrightarrow (3) \Leftrightarrow (4) \Leftrightarrow (2)$ 和 $(3) \Leftrightarrow (5)$。

    **(3) $\Rightarrow$ (1)**：若 $\rho(C) < 1$，由定理 69.2，$(I-C)^{-1} = \sum_{k=0}^{\infty} C^k \geq 0$。

    **(1) $\Rightarrow$ (3)**：设 $(I-C)^{-1} \geq 0$。假设 $\rho(C) \geq 1$，则 $\rho(C)$ 是 $C$ 的特征值（由 Perron-Frobenius 定理对非负矩阵成立），设对应特征向量为 $\mathbf{v} \geq 0, \mathbf{v} \neq 0$。则 $(I-C)\mathbf{v} = (1 - \rho(C))\mathbf{v}$，当 $\rho(C) = 1$ 时 $(I-C)\mathbf{v} = 0$，即 $(I-C)$ 不可逆，矛盾。当 $\rho(C) > 1$ 时，$(I-C)\mathbf{v} = (1-\rho(C))\mathbf{v} < 0$，则 $(I-C)^{-1}$ 将负向量映到非负向量，但 $(I-C)^{-1}[(1-\rho(C))\mathbf{v}] = \mathbf{v} \geq 0$ 而 $(1-\rho(C))\mathbf{v} \leq 0$，当 $(I-C)^{-1} \geq 0$ 时这要求 $(1-\rho(C))\mathbf{v} \geq 0$，与 $\rho(C)>1$ 矛盾。故 $\rho(C) < 1$。

    **(3) $\Leftrightarrow$ (4)**：$B = I - C$ 是 Z-矩阵（非对角元素非正，因 $C \geq 0$）。由 M-矩阵理论（Ch38），Z-矩阵 $B$ 是 M-矩阵当且仅当 $B^{-1} \geq 0$，当且仅当 $B$ 的特征值实部均为正，当且仅当 $\rho(C) < 1$。

    **(4) $\Leftrightarrow$ (2)**：由 M-矩阵理论，M-矩阵的所有顺序主子式为正（这是 M-矩阵的一个等价刻画）。

    **(3) $\Rightarrow$ (5)**：取 $\mathbf{x} = (I-C)^{-1}\mathbf{e}$，其中 $\mathbf{e} = (1,\ldots,1)^T$。由 $(I-C)^{-1} \geq 0$ 且列和为正（因为 $(I-C)^{-1}\mathbf{e} = \mathbf{e} + C\mathbf{e} + C^2\mathbf{e} + \cdots > 0$），得 $\mathbf{x} > 0$，且 $B\mathbf{x} = \mathbf{e} > 0$。

    **(5) $\Rightarrow$ (3)**：若存在 $\mathbf{x} > 0$ 使 $(I-C)\mathbf{x} > 0$，即 $\mathbf{x} > C\mathbf{x}$。设 $\rho(C) \geq 1$，由 Perron-Frobenius 定理存在非负特征向量 $\mathbf{v}$ 使 $C\mathbf{v} = \rho(C)\mathbf{v} \geq \mathbf{v}$。但这与存在 $\mathbf{x} > 0$ 使 $\mathbf{x} > C\mathbf{x}$ 的结论矛盾（可通过适当的比较论证得出）。故 $\rho(C) < 1$。

    $\blacksquare$

!!! example "例 69.3"
    验证例 69.1 中消费矩阵的 Hawkins-Simon 条件。

    $$B = I - C = \begin{pmatrix} 0.8 & -0.3 \\ -0.4 & 0.9 \end{pmatrix}$$

    - 第一个顺序主子式：$b_{11} = 0.8 > 0$ ✓
    - 第二个顺序主子式：$\det(B) = 0.72 - 0.12 = 0.60 > 0$ ✓

    Hawkins-Simon 条件满足，经济系统是生产性的。

    **经济解释**：$b_{11} = 1 - c_{11} = 0.8 > 0$ 意味着农业部门不需要消耗超过自身产出的投入来运转；$\det(B) > 0$ 意味着整个两部门经济在考虑相互依赖后仍有正的净产出。这就是"每个子经济体都必须是生产性的"这一直觉的精确数学表述。

---

## 69.3 封闭 Leontief 模型

<div class="context-flow" markdown>

**核心问题**：当没有外部需求（所有部门的产出都被经济体内部消耗）时，均衡产出如何确定？

</div>

!!! definition "定义 69.3 (封闭 Leontief 模型)"
    在封闭模型中，劳动力部门也被纳入消费矩阵，外部需求为零：

    $$\mathbf{x} = C\mathbf{x}$$

    即 $\mathbf{x}$ 是 $C$ 的属于特征值 $\lambda = 1$ 的特征向量。

在封闭模型中，消费矩阵 $C$ 的每一列之和等于 1（所有产出都被消耗），因此 $C$ 是列随机矩阵。

!!! theorem "定理 69.4 (封闭模型的均衡存在性)"
    设 $C \geq 0$ 为不可约的列随机矩阵（每列之和等于 1）。则：

    1. $\lambda = 1$ 是 $C$ 的特征值，且 $\rho(C) = 1$
    2. 存在唯一（在标量倍数意义下）的正特征向量 $\mathbf{x} > 0$ 使得 $C\mathbf{x} = \mathbf{x}$

??? proof "证明"
    由于 $C$ 是列随机矩阵，$\mathbf{e}^T C = \mathbf{e}^T$（其中 $\mathbf{e} = (1,\ldots,1)^T$），即 $\mathbf{e}^T$ 是 $C^T$ 属于特征值 1 的左特征向量。因此 1 是 $C^T$（从而也是 $C$）的特征值。

    又 $C \geq 0$ 且不可约，由 Perron-Frobenius 定理（Ch17），$\rho(C)$ 是 $C$ 的特征值且对应正特征向量。由于 $C$ 列随机，$\rho(C) \geq 1$（因为 1 已经是特征值）。又由 Perron-Frobenius 定理，$\rho(C) \leq \max_j \sum_i c_{ij} = 1$。故 $\rho(C) = 1$。

    由 Perron-Frobenius 定理的唯一性部分，对应于 $\rho(C) = 1$ 的正特征向量在标量倍数意义下唯一。

    $\blacksquare$

!!! example "例 69.4"
    三部门封闭经济的消费矩阵：

    $$C = \begin{pmatrix} 0.2 & 0.3 & 0.3 \\ 0.5 & 0.2 & 0.3 \\ 0.3 & 0.5 & 0.4 \end{pmatrix}$$

    每列之和为 1，是列随机矩阵。求解 $(C - I)\mathbf{x} = \mathbf{0}$：

    $$C - I = \begin{pmatrix} -0.8 & 0.3 & 0.3 \\ 0.5 & -0.8 & 0.3 \\ 0.3 & 0.5 & -0.6 \end{pmatrix}$$

    通过行化简（注意 $\det(C-I)=0$，秩为 2）：

    取 $x_3 = t$ 为自由变量，解得 $\mathbf{x} = t \begin{pmatrix} 0.78 \\ 0.82 \\ 1.00 \end{pmatrix}$（近似值）。

    均衡产出比例约为 $0.78 : 0.82 : 1.00$，即第三部门的产出最大。通过选取适当的标量 $t$，可以确定各部门的绝对产出水平。

---

## 69.4 Sraffa 模型与价格方程

<div class="context-flow" markdown>

**核心问题**：给定技术（投入产出关系）和工资水平，各商品的均衡价格如何确定？利润率与工资之间有什么关系？

</div>

Piero Sraffa 在其 1960 年的经典著作《用商品生产商品》中提出了一种价格决定理论，与 Leontief 模型互为对偶。

!!! definition "定义 69.4 (Sraffa 价格方程)"
    设 $A$ 为投入产出系数矩阵，$\mathbf{l}$ 为劳动系数向量（$l_j$ 为生产一单位第 $j$ 种商品所需的劳动），$r$ 为均一利润率，$w$ 为工资率，$\mathbf{p}$ 为价格向量。则 Sraffa 价格方程为：

    $$\mathbf{p}^T = \mathbf{p}^T(1 + r)A + w\mathbf{l}^T$$

    即每种商品的价格等于其投入成本（按利润率加成）加上劳动成本。

!!! theorem "定理 69.5 (利润率-工资率的权衡)"
    在 Sraffa 模型中，若 $A$ 是生产性的（$\rho(A) < 1$），则：

    1. 最大利润率 $R$（当 $w = 0$ 时）满足 $(1 + R) \rho(A) = 1$，即 $R = \frac{1}{\rho(A)} - 1$
    2. 当 $r = 0$ 时，$\mathbf{p}^T = w\mathbf{l}^T(I - A)^{-1}$（劳动价值论）
    3. 利润率 $r$ 与工资率 $w$ 之间存在反向关系

??? proof "证明"
    **(1)** 当 $w = 0$ 时，$\mathbf{p}^T = (1+r)\mathbf{p}^T A$，即 $\frac{1}{1+r}$ 是 $A$ 的特征值，$\mathbf{p}^T$ 是对应的左特征向量。为使 $\mathbf{p} > 0$，由 Perron-Frobenius 定理，需要 $\frac{1}{1+r} = \rho(A)$，故 $r = R = \frac{1}{\rho(A)} - 1$。

    **(2)** 当 $r = 0$ 时，$\mathbf{p}^T = \mathbf{p}^T A + w\mathbf{l}^T$，即 $\mathbf{p}^T(I - A) = w\mathbf{l}^T$，故 $\mathbf{p}^T = w\mathbf{l}^T(I - A)^{-1}$。

    **(3)** 从价格方程 $\mathbf{p}^T[I - (1+r)A] = w\mathbf{l}^T$，当 $r < R$ 时 $I - (1+r)A$ 仍然是 M-矩阵（可逆且逆非负），故 $\mathbf{p}^T = w\mathbf{l}^T[I - (1+r)A]^{-1}$。固定价格标准化条件（如 $\mathbf{p}^T\mathbf{q} = 1$ 对某标准商品篮子 $\mathbf{q}$），$w$ 随 $r$ 增大而减小。

    $\blacksquare$

!!! definition "定义 69.5 (标准商品)"
    Sraffa 定义**标准商品**为投入产出矩阵 $A$ 的 Perron 特征向量 $\mathbf{q}^*$（满足 $A\mathbf{q}^* = \rho(A)\mathbf{q}^*$）。以标准商品为价格计量单位时，利润率-工资率关系简化为线性关系：

    $$w = 1 - \frac{r}{R}$$

!!! example "例 69.5"
    设两部门经济的投入矩阵为 $A = \begin{pmatrix} 0.1 & 0.2 \\ 0.3 & 0.2 \end{pmatrix}$，劳动系数 $\mathbf{l}^T = (0.5, 0.4)$。

    首先计算 $\rho(A)$：$A$ 的特征多项式为 $\lambda^2 - 0.3\lambda - 0.04 = 0$，解得 $\lambda = \frac{0.3 + \sqrt{0.09 + 0.16}}{2} = \frac{0.3 + 0.5}{2} = 0.4$。

    最大利润率 $R = \frac{1}{0.4} - 1 = 1.5$，即 150%。

---

## 69.5 线性交换模型

<div class="context-flow" markdown>

**核心问题**：在一个纯交换经济中，没有生产活动，各经济主体通过交换各自的初始禀赋来实现均衡。均衡价格如何用线性代数来确定？

</div>

!!! definition "定义 69.6 (线性交换模型)"
    设有 $n$ 种商品和 $n$ 个交易者。交换矩阵 $E = (e_{ij})$ 定义为：$e_{ij}$ 是交易者 $j$ 用其收入购买商品 $i$ 的比例。则均衡价格向量 $\mathbf{p}$ 满足：

    $$\mathbf{p} = E\mathbf{p}$$

    即 $\mathbf{p}$ 是 $E$ 的属于特征值 1 的特征向量。

注意交换矩阵 $E$ 是列随机矩阵（每列之和为 1，因为每个交易者的全部收入都用于购买各种商品）。

!!! theorem "定理 69.6 (交换均衡的存在性)"
    若交换矩阵 $E \geq 0$ 不可约且列随机，则存在唯一的正均衡价格向量 $\mathbf{p} > 0$（在标量倍数意义下）。

??? proof "证明"
    这直接来自 Perron-Frobenius 定理。$E$ 不可约且非负，$\rho(E) = 1$（因为 $E$ 列随机），故存在唯一正特征向量 $\mathbf{p} > 0$ 使得 $E\mathbf{p} = \mathbf{p}$。

    $\blacksquare$

!!! example "例 69.6"
    三种商品的交换矩阵：

    $$E = \begin{pmatrix} 0.5 & 0.2 & 0.1 \\ 0.3 & 0.6 & 0.3 \\ 0.2 & 0.2 & 0.6 \end{pmatrix}$$

    每列之和为 1。解方程 $(E - I)\mathbf{p} = \mathbf{0}$：

    $$E - I = \begin{pmatrix} -0.5 & 0.2 & 0.1 \\ 0.3 & -0.4 & 0.3 \\ 0.2 & 0.2 & -0.4 \end{pmatrix}$$

    化简后得均衡价格比例约为 $\mathbf{p} \propto (0.48, 0.92, 0.68)^T$。标准化为 $p_1 + p_2 + p_3 = 1$ 得 $\mathbf{p} \approx (0.231, 0.442, 0.327)^T$。

---

## 69.6 二人零和博弈

<div class="context-flow" markdown>

**核心问题**：在完全对立的竞争情形（一方的收益恰好是另一方的损失）中，理性的参与者应该如何决策？

</div>

!!! definition "定义 69.7 (二人零和博弈)"
    二人零和博弈由**支付矩阵** $A \in \mathbb{R}^{m \times n}$ 定义。玩家 I 有 $m$ 个纯策略，玩家 II 有 $n$ 个纯策略。当玩家 I 选择策略 $i$、玩家 II 选择策略 $j$ 时，玩家 I 获得支付 $a_{ij}$，玩家 II 获得 $-a_{ij}$。

    **混合策略**：玩家 I 的混合策略是概率向量 $\mathbf{x} \in \Delta_m = \{\mathbf{x} \in \mathbb{R}^m : \mathbf{x} \geq 0, \sum x_i = 1\}$；玩家 II 的混合策略是 $\mathbf{y} \in \Delta_n$。

    混合策略下的期望支付为 $\mathbf{x}^T A \mathbf{y}$。

!!! definition "定义 69.8 (鞍点)"
    若存在 $(\mathbf{x}^*, \mathbf{y}^*)$ 使得对所有 $\mathbf{x} \in \Delta_m$ 和 $\mathbf{y} \in \Delta_n$：

    $$\mathbf{x}^T A \mathbf{y}^* \leq (\mathbf{x}^*)^T A \mathbf{y}^* \leq (\mathbf{x}^*)^T A \mathbf{y}$$

    则 $(\mathbf{x}^*, \mathbf{y}^*)$ 称为鞍点，$v = (\mathbf{x}^*)^T A \mathbf{y}^*$ 称为博弈值。

!!! theorem "定理 69.7 (Von Neumann 极小极大定理)"
    对任意 $m \times n$ 实矩阵 $A$：

    $$\max_{\mathbf{x} \in \Delta_m} \min_{\mathbf{y} \in \Delta_n} \mathbf{x}^T A \mathbf{y} = \min_{\mathbf{y} \in \Delta_n} \max_{\mathbf{x} \in \Delta_m} \mathbf{x}^T A \mathbf{y}$$

    即博弈必存在混合策略鞍点。

??? proof "证明"
    **步骤 1**：不等式 $\max \min \leq \min \max$ 总是成立的。对任意 $\mathbf{x}, \mathbf{y}$：

    $$\min_{\mathbf{y}'} \mathbf{x}^T A \mathbf{y}' \leq \mathbf{x}^T A \mathbf{y} \leq \max_{\mathbf{x}'} (\mathbf{x}')^T A \mathbf{y}$$

    对左端取 $\max_{\mathbf{x}}$，右端取 $\min_{\mathbf{y}}$，即得弱不等式。

    **步骤 2**：反向不等式的证明。定义：

    $$\underline{v} = \max_{\mathbf{x} \in \Delta_m} \min_{\mathbf{y} \in \Delta_n} \mathbf{x}^T A \mathbf{y}, \quad \overline{v} = \min_{\mathbf{y} \in \Delta_n} \max_{\mathbf{x} \in \Delta_m} \mathbf{x}^T A \mathbf{y}$$

    注意 $\min_{\mathbf{y}} \mathbf{x}^T A \mathbf{y} = \min_j (\mathbf{x}^T A)_j$（因为 $\mathbf{y}$ 的最优选择是集中于使 $\mathbf{x}^T A \mathbf{y}$ 最小的纯策略）。类似地，$\max_{\mathbf{x}} \mathbf{x}^T A \mathbf{y} = \max_i (A\mathbf{y})_i$。

    定义集合 $S = \{(\mathbf{x}^T A)^T : \mathbf{x} \in \Delta_m\} \subset \mathbb{R}^n$ 和 $T = \{\mathbf{z} \in \mathbb{R}^n : z_j \leq \underline{v}, \forall j\}$。若 $S \cap T = \emptyset$，则由凸集分离定理，存在超平面 $\mathbf{y}$ 分离 $S$ 和 $T$，但这会导致 $\max_{\mathbf{x}} \mathbf{x}^T A \mathbf{y} \leq \underline{v}$ 对某个 $\mathbf{y} \in \Delta_n$，从而 $\overline{v} \leq \underline{v}$，结合步骤 1 得 $\underline{v} = \overline{v}$。若 $S \cap T \neq \emptyset$，同样得 $\underline{v} = \overline{v}$。

    详细的严格证明可以通过线性规划的强对偶定理得到（见 69.7 节）。

    $\blacksquare$

!!! example "例 69.7"
    考虑支付矩阵：

    $$A = \begin{pmatrix} 3 & -1 \\ -2 & 4 \end{pmatrix}$$

    **求解混合策略均衡**。设玩家 I 的混合策略为 $\mathbf{x} = (p, 1-p)^T$，玩家 II 的为 $\mathbf{y} = (q, 1-q)^T$。

    玩家 I 使玩家 II 无差异：

    $$\mathbf{x}^T A \mathbf{e}_1 = \mathbf{x}^T A \mathbf{e}_2$$

    $$3p - 2(1-p) = -p + 4(1-p)$$

    $$5p - 2 = -5p + 4 \implies 10p = 6 \implies p = 0.6$$

    类似地，$q = 0.6$。

    博弈值 $v = \mathbf{x}^T A \mathbf{y} = 0.6 \times 0.6 \times 3 + 0.6 \times 0.4 \times (-1) + 0.4 \times 0.6 \times (-2) + 0.4 \times 0.4 \times 4 = 1.08 - 0.24 - 0.48 + 0.64 = 1.0$。

---

## 69.7 线性规划与对偶

<div class="context-flow" markdown>

**核心问题**：如何用矩阵语言表达线性规划问题？对偶问题的经济意义是什么？

</div>

!!! definition "定义 69.9 (线性规划的标准形式)"
    **原问题**（Primal）：

    $$\max \quad \mathbf{c}^T\mathbf{x} \quad \text{s.t.} \quad A\mathbf{x} \leq \mathbf{b}, \quad \mathbf{x} \geq \mathbf{0}$$

    其中 $\mathbf{c} \in \mathbb{R}^n$ 为目标函数系数，$A \in \mathbb{R}^{m \times n}$ 为约束矩阵，$\mathbf{b} \in \mathbb{R}^m$ 为资源限制。

    **对偶问题**（Dual）：

    $$\min \quad \mathbf{b}^T\mathbf{y} \quad \text{s.t.} \quad A^T\mathbf{y} \geq \mathbf{c}, \quad \mathbf{y} \geq \mathbf{0}$$

!!! theorem "定理 69.8 (弱对偶定理)"
    若 $\mathbf{x}$ 是原问题的可行解，$\mathbf{y}$ 是对偶问题的可行解，则：

    $$\mathbf{c}^T\mathbf{x} \leq \mathbf{b}^T\mathbf{y}$$

??? proof "证明"
    $$\mathbf{c}^T\mathbf{x} \leq (A^T\mathbf{y})^T\mathbf{x} = \mathbf{y}^T A\mathbf{x} \leq \mathbf{y}^T\mathbf{b} = \mathbf{b}^T\mathbf{y}$$

    第一个不等式由 $A^T\mathbf{y} \geq \mathbf{c}$ 且 $\mathbf{x} \geq 0$ 得出，第二个不等式由 $A\mathbf{x} \leq \mathbf{b}$ 且 $\mathbf{y} \geq 0$ 得出。

    $\blacksquare$

!!! theorem "定理 69.9 (强对偶定理)"
    若原问题和对偶问题都有可行解，则它们都有最优解，且最优值相等：

    $$\max_{\mathbf{x}} \mathbf{c}^T\mathbf{x} = \min_{\mathbf{y}} \mathbf{b}^T\mathbf{y}$$

!!! theorem "定理 69.10 (互补松弛条件)"
    $\mathbf{x}^*$ 和 $\mathbf{y}^*$ 分别是原问题和对偶问题的最优解，当且仅当它们是可行的且满足：

    $$y_i^* \left(b_i - \sum_j a_{ij} x_j^*\right) = 0, \quad \forall i$$

    $$x_j^* \left(\sum_i a_{ij} y_i^* - c_j\right) = 0, \quad \forall j$$

??? proof "证明"
    由强对偶定理，$\mathbf{c}^T\mathbf{x}^* = \mathbf{b}^T\mathbf{y}^*$。在弱对偶定理的证明中：

    $$\mathbf{c}^T\mathbf{x}^* \leq (\mathbf{y}^*)^T A\mathbf{x}^* \leq (\mathbf{y}^*)^T\mathbf{b}$$

    由于左端等于右端，中间的两个不等式必须取等号。

    第一个等号要求 $(A^T\mathbf{y}^* - \mathbf{c})^T\mathbf{x}^* = 0$。由于两个因子都非负，必须 $x_j^*(A^T\mathbf{y}^* - \mathbf{c})_j = 0$ 对每个 $j$。

    第二个等号要求 $(\mathbf{y}^*)^T(\mathbf{b} - A\mathbf{x}^*) = 0$，类似得 $y_i^*(\mathbf{b} - A\mathbf{x}^*)_i = 0$ 对每个 $i$。

    $\blacksquare$

!!! definition "定义 69.10 (影子价格)"
    对偶变量 $y_i^*$ 称为第 $i$ 个约束对应的**影子价格**（shadow price）。$y_i^*$ 表示将第 $i$ 种资源的供给量 $b_i$ 增加一个单位时，目标函数最优值的增加量：

    $$y_i^* = \frac{\partial v^*}{\partial b_i}$$

    其中 $v^* = \max \mathbf{c}^T\mathbf{x}$ 是最优值。

!!! example "例 69.8"
    一个工厂生产两种产品。利润分别为 5 和 4 元/件。需要两种资源，供给量分别为 6 和 8 单位。

    原问题：

    $$\max \quad 5x_1 + 4x_2 \quad \text{s.t.} \quad x_1 + x_2 \leq 6, \quad 2x_1 + x_2 \leq 8, \quad x_1, x_2 \geq 0$$

    对偶问题：

    $$\min \quad 6y_1 + 8y_2 \quad \text{s.t.} \quad y_1 + 2y_2 \geq 5, \quad y_1 + y_2 \geq 4, \quad y_1, y_2 \geq 0$$

    解原问题：可行域的顶点为 $(0,0), (4,0), (2,4), (0,6)$。在各顶点处：$0, 20, 26, 24$。最优解 $\mathbf{x}^* = (2, 4)^T$，最优值 26。

    解对偶问题：$\mathbf{y}^* = (3, 1)^T$，最优值 $6 \times 3 + 8 \times 1 = 26$。

    影子价格：第一种资源的影子价格为 3，即增加 1 单位第一种资源，利润增加 3 元。

!!! theorem "定理 69.11 (博弈论与线性规划的等价)"
    二人零和博弈的求解可以转化为线性规划问题。给定支付矩阵 $A$，玩家 I 的最优混合策略可通过以下线性规划求解：

    $$\max \quad v \quad \text{s.t.} \quad A^T\mathbf{x} \geq v\mathbf{e}, \quad \mathbf{e}^T\mathbf{x} = 1, \quad \mathbf{x} \geq \mathbf{0}$$

    其对偶问题给出玩家 II 的最优策略。Von Neumann 极小极大定理是强对偶定理的推论。

---

## 69.8 Markov 链在经济学中的应用

<div class="context-flow" markdown>

**核心问题**：如何用 Markov 链模型描述经济状态的动态演化？收入分配的长期趋势如何？

</div>

!!! definition "定义 69.11 (收入流动矩阵)"
    **收入流动矩阵**（Income Mobility Matrix）$P$ 是一个行随机矩阵，其中 $P_{ij}$ 表示当前处于收入阶层 $i$ 的个体在下一代处于阶层 $j$ 的概率。

!!! theorem "定理 69.12 (收入分布的收敛)"
    若收入流动矩阵 $P$ 不可约且非周期，则存在唯一的稳态分布 $\boldsymbol{\pi}$ 满足 $\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T$，且对任意初始分布 $\boldsymbol{\pi}_0$：

    $$\lim_{n \to \infty} \boldsymbol{\pi}_0^T P^n = \boldsymbol{\pi}^T$$

    $\boldsymbol{\pi}$ 描述长期均衡下的收入分布。

??? proof "证明"
    这是 Markov 链遍历定理的直接应用（详见 Ch71）。$P$ 不可约且非周期保证了遍历性，收敛性由 $P$ 的谱分解和 $|\lambda_2| < 1$ 保证。

    $\blacksquare$

!!! example "例 69.9"
    一个简化的三阶层收入流动矩阵（低、中、高收入）：

    $$P = \begin{pmatrix} 0.7 & 0.2 & 0.1 \\ 0.1 & 0.7 & 0.2 \\ 0.05 & 0.15 & 0.8 \end{pmatrix}$$

    解稳态方程 $\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T$，$\boldsymbol{\pi}^T \mathbf{e} = 1$：

    设 $\boldsymbol{\pi} = (\pi_1, \pi_2, \pi_3)^T$，由 $\boldsymbol{\pi}^T(P - I) = \mathbf{0}^T$：

    $$-0.3\pi_1 + 0.1\pi_2 + 0.05\pi_3 = 0$$
    $$0.2\pi_1 - 0.3\pi_2 + 0.15\pi_3 = 0$$
    $$\pi_1 + \pi_2 + \pi_3 = 1$$

    解得 $\boldsymbol{\pi} \approx (0.172, 0.345, 0.483)^T$。

    **经济解释**：长期均衡下，约 17.2% 的人口处于低收入阶层，34.5% 处于中等收入，48.3% 处于高收入。注意这与收入流动矩阵中高收入的"粘性"（$P_{33} = 0.8$）一致——高收入阶层的代际延续性最强。

!!! example "例 69.10"
    **收入流动性的度量**。收入流动矩阵 $P$ 的第二大特征值 $\lambda_2$ 度量了收入流动的速度。$\lambda_2$ 越小，收敛到稳态分布越快，即社会流动性越强。

    若 $P$ 为单位矩阵（完全不流动），则 $\lambda_2 = 1$。若 $P$ 的每一行都相同（完全流动），则 $\lambda_2 = 0$。

    对上例，可以计算 $P$ 的特征值为 $1, 0.68, 0.52$，第二大特征值 $\lambda_2 = 0.68$，说明代际收入流动性为中等水平。

!!! theorem "定理 69.13 (收入不平等的动态)"
    设 $\mathbf{w}$ 为各收入阶层的收入水平向量，$\boldsymbol{\pi}_t$ 为第 $t$ 代的收入分布。则第 $t$ 代的平均收入为：

    $$\bar{w}_t = \boldsymbol{\pi}_t^T \mathbf{w} = \boldsymbol{\pi}_0^T P^t \mathbf{w}$$

    长期平均收入为 $\bar{w}_\infty = \boldsymbol{\pi}^T \mathbf{w}$。收入不平等（如基尼系数）的动态变化也可以通过 $P$ 的谱结构来分析。

??? proof "证明"
    由 $\boldsymbol{\pi}_t^T = \boldsymbol{\pi}_0^T P^t$ 直接得到 $\bar{w}_t = \boldsymbol{\pi}_0^T P^t \mathbf{w}$。

    由遍历定理 $P^t \to \mathbf{e}\boldsymbol{\pi}^T$，故：

    $$\bar{w}_\infty = \lim_{t \to \infty} \boldsymbol{\pi}_0^T P^t \mathbf{w} = \boldsymbol{\pi}_0^T \mathbf{e} \boldsymbol{\pi}^T \mathbf{w} = \boldsymbol{\pi}^T \mathbf{w}$$

    其中用到了 $\boldsymbol{\pi}_0^T \mathbf{e} = 1$。

    $\blacksquare$

---

## 本章小结

本章展示了线性代数在经济学中的四个核心应用领域：

1. **投入产出分析**：Leontief 模型将经济部门间的相互依赖关系编码为消费矩阵 $C$，Hawkins-Simon 条件保证了经济系统的可行性。Neumann 级数 $(I-C)^{-1} = I + C + C^2 + \cdots$ 揭示了经济乘数效应的数学本质。

2. **价格理论**：Sraffa 模型从对偶角度研究价格决定，标准商品（Perron 特征向量）将利润率-工资率关系线性化。

3. **博弈论**：Von Neumann 极小极大定理保证了二人零和博弈混合策略均衡的存在性，其证明与线性规划的强对偶定理密切相关。

4. **动态经济分析**：Markov 链的收入流动矩阵描述了代际收入分布的演化，遍历定理保证了长期均衡分布的存在性。

这些应用的共同特征是：经济问题的结构可以自然地用非负矩阵、随机矩阵和 M-矩阵的语言来表达，而 Perron-Frobenius 定理和 M-矩阵理论则为经济均衡的存在性、唯一性和稳定性提供了严格的数学基础。
