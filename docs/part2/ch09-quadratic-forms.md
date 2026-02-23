# 第 9 章 二次型

<div class="context-flow" markdown>

**前置**：Ch8 对称矩阵谱定理

**本章脉络**：$\mathbf{x}^TA\mathbf{x}$ → 配方法/正交法化标准形 → **惯性定理**（签名不变） → 正定性判别 → 几何（椭球/双曲面）

**延伸**：二次型在优化（二次规划）、统计（$\chi^2$ 分布）、微分几何（黎曼度量）中无处不在；双线性型的反对称情形给出辛形式，辛几何是经典力学（Hamilton 方程）和量子场论的数学语言；Hermite 型是量子力学中可观测量的数学结构

</div>

二次型（quadratic form）是二次齐次多项式的代数理论，它与对称矩阵和内积空间有着深刻的联系。二次型的研究不仅是线性代数的重要组成部分，而且在微分几何、优化理论、统计学和物理学中有着广泛的应用。本章将系统地研究二次型的定义、化简方法、惯性定理以及正定性判别等核心理论。

---

## 9.1 二次型的定义

<div class="context-flow" markdown>

**对称矩阵** $A$ ↔ 二次型 $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ 是一一对应；交叉项 $x_ix_j$ 拆分为 $a_{ij} = a_{ji}$

</div>

!!! definition "定义 9.1 (二次型)"
    设 $\mathbb{F} = \mathbb{R}$（或 $\mathbb{C}$）。$n$ 个变量 $x_1, x_2, \ldots, x_n$ 上的**二次型**（quadratic form）是如下形式的二次齐次多项式：

    $$Q(x_1, x_2, \ldots, x_n) = \sum_{i=1}^n \sum_{j=1}^n a_{ij} x_i x_j$$

    其中 $a_{ij} \in \mathbb{F}$。用向量记号，设 $\mathbf{x} = (x_1, \ldots, x_n)^T$，则

    $$Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}$$

    其中 $A = (a_{ij})_{n \times n}$。

!!! definition "定义 9.2 (二次型的矩阵)"
    对于二次型 $Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}$，我们总可以假设矩阵 $A$ 是**对称的**。事实上，对任意矩阵 $A$，$\mathbf{x}^T A \mathbf{x} = \mathbf{x}^T \left(\frac{A + A^T}{2}\right) \mathbf{x}$，而 $\frac{A + A^T}{2}$ 是对称矩阵。

    称对称矩阵 $A$ 为二次型 $Q$ 的**矩阵**，$\operatorname{rank}(A)$ 称为二次型 $Q$ 的**秩**（rank）。

!!! theorem "定理 9.1 (二次型与对称矩阵的一一对应)"
    $n$ 元实二次型 $Q(\mathbf{x})$ 与 $n$ 阶实对称矩阵 $A$ 之间存在一一对应关系：$Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}$。

??? proof "证明"
    给定对称矩阵 $A = (a_{ij})$，$\mathbf{x}^T A \mathbf{x} = \sum_{i,j} a_{ij}x_ix_j$ 是二次型。

    反之，给定二次型 $Q(\mathbf{x}) = \sum_{i \leq j} c_{ij}x_ix_j$（其中 $c_{ii}$ 是 $x_i^2$ 的系数，$c_{ij}$（$i < j$）是 $x_ix_j$ 的系数），令对称矩阵 $A$ 的元素为 $a_{ii} = c_{ii}$，$a_{ij} = a_{ji} = c_{ij}/2$（$i < j$），则 $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$。

    唯一性：若 $\mathbf{x}^TA\mathbf{x} = \mathbf{x}^TB\mathbf{x}$ 对所有 $\mathbf{x}$ 成立，$A, B$ 均对称，则 $\mathbf{x}^T(A-B)\mathbf{x} = 0$ 对所有 $\mathbf{x}$ 成立。取 $\mathbf{x} = \mathbf{e}_i$ 得 $a_{ii} = b_{ii}$；取 $\mathbf{x} = \mathbf{e}_i + \mathbf{e}_j$ 得 $a_{ij} + a_{ji} = b_{ij} + b_{ji}$，由对称性得 $a_{ij} = b_{ij}$。$\blacksquare$

!!! example "例 9.1"
    二次型 $Q(x_1, x_2, x_3) = 2x_1^2 + 3x_2^2 - x_3^2 + 4x_1x_2 - 6x_1x_3 + 2x_2x_3$ 的对称矩阵为

    $$A = \begin{pmatrix} 2 & 2 & -3 \\ 2 & 3 & 1 \\ -3 & 1 & -1 \end{pmatrix}$$

    注意交叉项 $4x_1x_2$ 拆分为 $a_{12} = a_{21} = 2$。

!!! example "例 9.2"
    $\mathbb{R}^2$ 上的二次型 $Q(x_1, x_2) = x_1^2 + x_2^2$ 对应矩阵 $A = I_2$，几何上表示以原点为圆心的圆。二次型 $Q(x_1, x_2) = x_1^2 - x_2^2$ 对应矩阵 $A = \operatorname{diag}(1, -1)$，几何上表示双曲线。

---

## 9.2 二次型的标准形

<div class="context-flow" markdown>

消去交叉项 → 只留 $d_i y_i^2$ → **配方法**（Lagrange）是构造性工具，任意二次型均可化标准形

</div>

!!! definition "定义 9.3 (标准形)"
    如果二次型 $Q(\mathbf{x})$ 只含平方项（没有交叉项），即

    $$Q(\mathbf{x}) = d_1 x_1^2 + d_2 x_2^2 + \cdots + d_n x_n^2$$

    则称 $Q$ 为**标准形**（canonical form，或对角形），其矩阵为对角矩阵 $\operatorname{diag}(d_1, \ldots, d_n)$。

!!! definition "定义 9.4 (非退化线性替换)"
    设 $\mathbf{x} = C\mathbf{y}$，其中 $C$ 是 $n \times n$ 可逆矩阵。这称为**非退化线性替换**（nonsingular linear substitution）。在此替换下，

    $$Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x} = (C\mathbf{y})^TA(C\mathbf{y}) = \mathbf{y}^T(C^TAC)\mathbf{y}$$

    新二次型的矩阵为 $B = C^TAC$。

### 配方法

!!! theorem "定理 9.2 (配方法化标准形)"
    任意实二次型 $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$（$A$ 为对称矩阵）都可以通过非退化线性替换化为标准形。

??? proof "证明"
    用配方法（拉格朗日配方法）。分两种情况：

    **情况 1：** 若某个 $a_{ii} \neq 0$（不妨设 $a_{11} \neq 0$），则将含 $x_1$ 的所有项配方：

    $$Q = a_{11}\left(x_1 + \frac{a_{12}}{a_{11}}x_2 + \cdots + \frac{a_{1n}}{a_{11}}x_n\right)^2 + Q_1(x_2, \ldots, x_n)$$

    令 $y_1 = x_1 + \frac{a_{12}}{a_{11}}x_2 + \cdots + \frac{a_{1n}}{a_{11}}x_n$，$y_i = x_i$（$i \geq 2$），此为非退化线性替换。对 $Q_1$ 递归进行配方。

    **情况 2：** 若所有 $a_{ii} = 0$，但存在 $a_{ij} \neq 0$（$i \neq j$）。令 $x_i = y_i + y_j$，$x_j = y_i - y_j$，其余 $x_k = y_k$，则 $2a_{ij}x_ix_j = 2a_{ij}(y_i^2 - y_j^2)$，产生平方项，回到情况 1。

    由归纳法，经过有限步后，$Q$ 化为标准形。$\blacksquare$

!!! example "例 9.3"
    用配方法化标准形：$Q(x_1, x_2, x_3) = 2x_1x_2 + 2x_1x_3 - 6x_2x_3$。

    所有平方项系数为零（情况 2）。令 $x_1 = y_1 + y_2$，$x_2 = y_1 - y_2$，$x_3 = y_3$：

    $$Q = 2(y_1+y_2)(y_1-y_2) + 2(y_1+y_2)y_3 - 6(y_1-y_2)y_3$$

    $$= 2y_1^2 - 2y_2^2 + 2y_1y_3 + 2y_2y_3 - 6y_1y_3 + 6y_2y_3$$

    $$= 2y_1^2 - 4y_1y_3 - 2y_2^2 + 8y_2y_3$$

    对 $y_1$ 配方：$2(y_1^2 - 2y_1y_3) = 2(y_1 - y_3)^2 - 2y_3^2$。

    对 $y_2$ 配方：$-2(y_2^2 - 4y_2y_3) = -2(y_2 - 2y_3)^2 + 8y_3^2$。

    $$Q = 2(y_1 - y_3)^2 - 2(y_2 - 2y_3)^2 + 6y_3^2$$

    令 $z_1 = y_1 - y_3$，$z_2 = y_2 - 2y_3$，$z_3 = y_3$，得标准形 $Q = 2z_1^2 - 2z_2^2 + 6z_3^2$。

---

## 9.3 正交化法化标准形

<div class="context-flow" markdown>

配方法的替换矩阵不唯一 → 用 Ch8 **谱定理** $A = Q\Lambda Q^T$ 做正交替换 → 标准形系数 = 特征值，变换 = 等距

</div>

配方法得到的标准形依赖于配方的顺序，变换矩阵不唯一。正交化法（利用谱定理）给出了一种"最自然"的标准形化简方法。

!!! theorem "定理 9.3 (正交对角化法)"
    设 $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$，$A$ 为 $n$ 阶实对称矩阵。由谱定理，存在正交矩阵 $P$ 使得

    $$P^TAP = \Lambda = \operatorname{diag}(\lambda_1, \lambda_2, \ldots, \lambda_n)$$

    其中 $\lambda_1, \ldots, \lambda_n$ 是 $A$ 的特征值。令 $\mathbf{x} = P\mathbf{y}$（正交替换），则

    $$Q = \mathbf{y}^T\Lambda\mathbf{y} = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \cdots + \lambda_n y_n^2$$

??? proof "证明"
    由实对称矩阵的谱定理（定理 8.15），存在正交矩阵 $P$（其列为 $A$ 的标准正交特征向量）使得 $P^TAP = \Lambda$。在替换 $\mathbf{x} = P\mathbf{y}$ 下：

    $$Q(\mathbf{x}) = (P\mathbf{y})^T A (P\mathbf{y}) = \mathbf{y}^T(P^TAP)\mathbf{y} = \mathbf{y}^T\Lambda\mathbf{y} = \sum_{i=1}^n \lambda_i y_i^2$$

    $\blacksquare$

!!! note "注"
    正交替换的优点是保持向量的长度和角度不变（它是等距变换），因此在几何应用中特别有意义。正交化法得到的标准形中的系数恰好是特征值，这在理论上非常自然。

!!! example "例 9.4"
    用正交替换化标准形：$Q(x_1, x_2) = 5x_1^2 + 4x_1x_2 + 8x_2^2$。

    对称矩阵 $A = \begin{pmatrix} 5 & 2 \\ 2 & 8 \end{pmatrix}$。

    特征多项式：$\det(A - \lambda I) = (5-\lambda)(8-\lambda) - 4 = \lambda^2 - 13\lambda + 36 = (\lambda - 4)(\lambda - 9)$。

    $\lambda_1 = 4$：$(A - 4I)\mathbf{x} = \mathbf{0}$ $\Rightarrow$ $\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\mathbf{x} = \mathbf{0}$，$\mathbf{v}_1 = \frac{1}{\sqrt{5}}(-2, 1)^T$。

    $\lambda_2 = 9$：$(A - 9I)\mathbf{x} = \mathbf{0}$ $\Rightarrow$ $\begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\mathbf{x} = \mathbf{0}$，$\mathbf{v}_2 = \frac{1}{\sqrt{5}}(1, 2)^T$。

    正交矩阵 $P = \frac{1}{\sqrt{5}}\begin{pmatrix} -2 & 1 \\ 1 & 2 \end{pmatrix}$，令 $\mathbf{x} = P\mathbf{y}$，得 $Q = 4y_1^2 + 9y_2^2$。

---

## 9.4 惯性定理

<div class="context-flow" markdown>

标准形的系数可变，但**正系数个数 $p$ 和负系数个数 $q$ 不变**——Sylvester 惯性定律是二次型理论的核心不变量

</div>

!!! definition "定义 9.5 (惯性指数)"
    设实二次型 $Q(\mathbf{x})$ 经非退化线性替换化为标准形

    $$Q = d_1 y_1^2 + d_2 y_2^2 + \cdots + d_r y_r^2$$

    其中 $d_i \neq 0$（$i = 1, \ldots, r$），$r = \operatorname{rank}(A)$。设其中正系数的个数为 $p$，负系数的个数为 $q = r - p$。则

    - $p$ 称为二次型的**正惯性指数**（positive index of inertia）；
    - $q$ 称为二次型的**负惯性指数**（negative index of inertia）；
    - $(p, q)$ 称为二次型的**符号差**或**签名**（signature）。

<div class="context-flow" markdown>

**洞察**：证明的核心是**维数论证**——$V_1 \cap V_2 \neq \{0\}$（$\dim V_1 + \dim V_2 > n$）导出矛盾，这一技巧在 Ch11 Eckart-Young 中再现

</div>

!!! theorem "定理 9.4 (Sylvester 惯性定律)"
    实二次型的标准形中正系数的个数 $p$ 和负系数的个数 $q$ 是不变的，与化标准形时所用的非退化线性替换无关。即 $p$ 和 $q$ 仅由二次型本身决定。

??? proof "证明"
    设二次型 $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ 经两种不同的非退化线性替换 $\mathbf{x} = C_1\mathbf{y}$ 和 $\mathbf{x} = C_2\mathbf{z}$ 分别化为标准形：

    $$Q = y_1^2 + \cdots + y_p^2 - y_{p+1}^2 - \cdots - y_r^2$$

    $$Q = z_1^2 + \cdots + z_s^2 - z_{s+1}^2 - \cdots - z_r^2$$

    （不失一般性，可将系数化为 $\pm 1$。）

    用反证法，假设 $p \neq s$，不妨设 $p > s$。

    设 $\mathbf{y} = C_1^{-1}\mathbf{x}$，$\mathbf{z} = C_2^{-1}\mathbf{x}$。考虑以下两个子空间：

    - $V_1 = \{\mathbf{x} \in \mathbb{R}^n : y_{p+1} = \cdots = y_n = 0\}$，$\dim V_1 = p$；
    - $V_2 = \{\mathbf{x} \in \mathbb{R}^n : z_1 = \cdots = z_s = 0\}$，$\dim V_2 = n - s$。

    由于 $\dim V_1 + \dim V_2 = p + (n-s) > n$（因 $p > s$），故 $V_1 \cap V_2 \neq \{\mathbf{0}\}$。

    设 $\mathbf{0} \neq \mathbf{x}_0 \in V_1 \cap V_2$。

    - 在 $V_1$ 中：$Q(\mathbf{x}_0) = y_1^2 + \cdots + y_p^2 > 0$（因 $\mathbf{x}_0 \neq \mathbf{0}$ 意味着至少一个 $y_i \neq 0$，$i \leq p$）；
    - 在 $V_2$ 中：$Q(\mathbf{x}_0) = -z_{s+1}^2 - \cdots - z_r^2 \leq 0$。

    矛盾！故 $p = s$。$\blacksquare$

!!! corollary "推论 9.1"
    两个实二次型等价（即可通过非退化线性替换相互转化）当且仅当它们有相同的秩和相同的正惯性指数（或等价地，相同的签名）。

!!! proposition "命题 9.1"
    实对称矩阵 $A$ 的正惯性指数等于 $A$ 的正特征值的个数（计重数），负惯性指数等于 $A$ 的负特征值的个数（计重数）。

??? proof "证明"
    由正交对角化，$Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ 经正交替换可化为 $\lambda_1 y_1^2 + \cdots + \lambda_n y_n^2$，其中 $\lambda_i$ 是 $A$ 的特征值。由惯性定律，正惯性指数等于正特征值个数，负惯性指数等于负特征值个数。$\blacksquare$

!!! example "例 9.5"
    二次型 $Q = x_1^2 + 4x_1x_2 + 4x_2^2 + 2x_3^2$ 的矩阵为

    $$A = \begin{pmatrix} 1 & 2 & 0 \\ 2 & 4 & 0 \\ 0 & 0 & 2 \end{pmatrix}$$

    特征多项式：$\det(A - \lambda I) = (2-\lambda)[(1-\lambda)(4-\lambda) - 4] = (2-\lambda)(\lambda^2 - 5\lambda) = -\lambda(2-\lambda)(\lambda - 5)$。

    特征值为 $\lambda_1 = 0, \lambda_2 = 2, \lambda_3 = 5$。正惯性指数 $p = 2$，负惯性指数 $q = 0$，秩 $r = 2$。

---

## 9.5 合同变换与合同矩阵

<div class="context-flow" markdown>

**相似** $P^{-1}AP$（保特征值）vs **合同** $C^TAC$（保签名）——合同是二次型的自然等价关系，对称性和秩不变但特征值可变

</div>

!!! definition "定义 9.6 (合同)"
    设 $A, B$ 是 $n$ 阶实方阵。若存在可逆矩阵 $C$ 使得

    $$B = C^TAC$$

    则称 $A$ 与 $B$ **合同**（congruent），记作 $A \simeq B$。$\mathbf{x} \mapsto C\mathbf{x}$ 称为**合同变换**（congruence transformation）。

!!! proposition "命题 9.2 (合同是等价关系)"
    矩阵的合同关系是等价关系，即满足：

    1. **自反性**：$A \simeq A$（取 $C = I$）；
    2. **对称性**：若 $A \simeq B$，则 $B \simeq A$；
    3. **传递性**：若 $A \simeq B$，$B \simeq D$，则 $A \simeq D$。

??? proof "证明"
    (2) 若 $B = C^TAC$，则 $A = (C^{-1})^T B C^{-1} = (C^{-1})^T B (C^{-1})$，故 $B \simeq A$。

    (3) 若 $B = C_1^TAC_1$，$D = C_2^TBC_2$，则 $D = C_2^T(C_1^TAC_1)C_2 = (C_1C_2)^T A (C_1C_2)$。$\blacksquare$

!!! theorem "定理 9.5 (合同标准形)"
    任意 $n$ 阶实对称矩阵 $A$（秩为 $r$）合同于

    $$\begin{pmatrix} I_p & & \\ & -I_q & \\ & & O_{n-r} \end{pmatrix}$$

    其中 $p$ 是正惯性指数，$q = r - p$ 是负惯性指数。此标准形由惯性定律保证唯一。

??? proof "证明"
    由配方法（定理 9.2），存在可逆矩阵 $C_1$ 使 $C_1^TAC_1 = \operatorname{diag}(d_1, \ldots, d_r, 0, \ldots, 0)$，其中 $d_i \neq 0$。不妨设 $d_1, \ldots, d_p > 0$，$d_{p+1}, \ldots, d_r < 0$。令

    $$C_2 = \operatorname{diag}\left(\frac{1}{\sqrt{|d_1|}}, \ldots, \frac{1}{\sqrt{|d_r|}}, 1, \ldots, 1\right)$$

    则 $C_2^T(\operatorname{diag}(d_1, \ldots, d_r, 0, \ldots, 0))C_2 = \operatorname{diag}(1, \ldots, 1, -1, \ldots, -1, 0, \ldots, 0)$。

    取 $C = C_1C_2$，即得所求。$\blacksquare$

!!! note "注"
    合同关系保持对称性和秩：若 $A$ 对称且 $B = C^TAC$，则 $B$ 也对称，且 $\operatorname{rank}(B) = \operatorname{rank}(A)$。但合同关系**不**保持特征值。相比之下，相似关系保持特征值但不一定保持对称性。

!!! example "例 9.6"
    矩阵 $A = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$ 的特征值为 $3$ 和 $-1$，故正惯性指数 $p = 1$，负惯性指数 $q = 1$。$A$ 合同于 $\begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$。

    验证：取 $C = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$（对应配方 $(x_1+2x_2)^2 - 3x_2^2$ 中再缩放），可以具体计算出合同变换矩阵。

---

## 9.6 正定二次型与正定矩阵

<div class="context-flow" markdown>

签名 $(n,0)$ ↔ 所有特征值 > 0 ↔ $A = C^TC$ ↔ 顺序主子式全正 → 正定矩阵直通 Ch10 **Cholesky 分解** $A = LL^T$

</div>

正定性是二次型和对称矩阵最重要的性质之一，在优化、统计、微分方程等领域有核心地位。

!!! definition "定义 9.7 (正定性分类)"
    设 $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ 是实二次型，$A$ 为 $n$ 阶实对称矩阵。称 $Q$（或 $A$）为：

    - **正定的**（positive definite）：若对所有 $\mathbf{x} \neq \mathbf{0}$ 有 $Q(\mathbf{x}) > 0$；
    - **半正定的**（positive semidefinite）：若对所有 $\mathbf{x}$ 有 $Q(\mathbf{x}) \geq 0$；
    - **负定的**（negative definite）：若对所有 $\mathbf{x} \neq \mathbf{0}$ 有 $Q(\mathbf{x}) < 0$；
    - **半负定的**（negative semidefinite）：若对所有 $\mathbf{x}$ 有 $Q(\mathbf{x}) \leq 0$；
    - **不定的**（indefinite）：若 $Q$ 既取正值也取负值。

<div class="context-flow" markdown>

**洞察**：五种等价条件统一了代数（特征值）、几何（$A = C^TC$）和组合（顺序主子式）三个视角

</div>

!!! theorem "定理 9.6 (正定的等价条件)"
    设 $A$ 是 $n$ 阶实对称矩阵。以下条件等价：

    1. $A$ 是正定的；
    2. $A$ 的所有特征值 $\lambda_1, \ldots, \lambda_n$ 都大于零；
    3. $A$ 的正惯性指数 $p = n$（即 $Q$ 的标准形为 $y_1^2 + \cdots + y_n^2$）；
    4. 存在可逆矩阵 $C$ 使得 $A = C^TC$；
    5. $A$ 的所有**顺序主子式**（leading principal minors）都大于零。

??? proof "证明"
    **(1)$\Leftrightarrow$(2)：** 若 $A$ 正定，设 $A\mathbf{v} = \lambda\mathbf{v}$，$\mathbf{v} \neq \mathbf{0}$，则 $\lambda\|\mathbf{v}\|^2 = \mathbf{v}^TA\mathbf{v} > 0$，故 $\lambda > 0$。反之，若所有特征值大于零，$Q(\mathbf{x}) = \mathbf{y}^T\Lambda\mathbf{y} = \sum \lambda_i y_i^2 > 0$（$\mathbf{x} \neq \mathbf{0}$）。

    **(2)$\Leftrightarrow$(3)：** 正惯性指数等于正特征值个数（命题 9.1）。

    **(1)$\Rightarrow$(4)：** 由正交对角化 $A = P\Lambda P^T$，$\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，$\lambda_i > 0$。令 $C = \Lambda^{1/2}P^T$（其中 $\Lambda^{1/2} = \operatorname{diag}(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_n})$），则 $C^TC = P\Lambda^{1/2}\Lambda^{1/2}P^T = P\Lambda P^T = A$。

    **(4)$\Rightarrow$(1)：** $\mathbf{x}^TA\mathbf{x} = \mathbf{x}^TC^TC\mathbf{x} = \|C\mathbf{x}\|^2 \geq 0$，且等号成立当且仅当 $C\mathbf{x} = \mathbf{0}$，即 $\mathbf{x} = \mathbf{0}$（因 $C$ 可逆）。

    **(1)$\Leftrightarrow$(5)：** 这是 Sylvester 判据，证明如下。设 $A$ 的 $k$ 阶顺序主子矩阵为 $A_k$，顺序主子式为 $\Delta_k = \det(A_k)$。

    若 $A$ 正定，则 $A_k$ 也正定（对 $\mathbf{x} = (x_1, \ldots, x_k, 0, \ldots, 0)^T$，$\mathbf{x}^TA\mathbf{x} = \mathbf{y}^TA_k\mathbf{y}$ 其中 $\mathbf{y} = (x_1, \ldots, x_k)^T$）。$A_k$ 正定 $\Rightarrow$ 特征值均正 $\Rightarrow$ $\Delta_k = \prod \lambda_i^{(k)} > 0$。

    反之，对 $n$ 用归纳法证明。$n=1$ 时 $\Delta_1 = a_{11} > 0$ 即正定。设 $n-1$ 时成立。由 $\Delta_1, \ldots, \Delta_{n-1} > 0$，$A_{n-1}$ 正定。利用 Schur 补可以证明 $A$ 正定。$\blacksquare$

!!! theorem "定理 9.7 (半正定的等价条件)"
    设 $A$ 是 $n$ 阶实对称矩阵。以下条件等价：

    1. $A$ 是半正定的；
    2. $A$ 的所有特征值 $\lambda_i \geq 0$；
    3. 存在矩阵 $B$（不一定可逆）使得 $A = B^TB$；
    4. $A$ 的所有**主子式**（principal minors，不仅是顺序主子式）都非负。

??? proof "证明"
    (1)$\Leftrightarrow$(2) 和 (1)$\Leftrightarrow$(3) 的证明类似正定情形。

    (4) 的必要性：$A$ 半正定 $\Rightarrow$ 每个主子矩阵也半正定 $\Rightarrow$ 主子式（= 特征值之积）$\geq 0$。$\blacksquare$

!!! note "注"
    判别负定性：$A$ 负定当且仅当 $-A$ 正定，等价于 $(-1)^k\Delta_k > 0$（$k = 1, \ldots, n$），即奇数阶顺序主子式为负，偶数阶顺序主子式为正。

!!! example "例 9.7"
    判断矩阵 $A = \begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}$ 的正定性。

    $\Delta_1 = 2 > 0$，$\Delta_2 = \det\begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix} = 3 > 0$，$\Delta_3 = \det(A) = 2(4-1) - (-1)(-2) = 6 - 2 = 4 > 0$。

    所有顺序主子式为正，故 $A$ 正定。

!!! example "例 9.8"
    判断 $A = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$ 的正定性。

    $\Delta_1 = 1 > 0$，$\Delta_2 = 1 - 4 = -3 < 0$。

    $A$ 不是正定的。事实上，$A$ 的特征值为 $3$ 和 $-1$，故 $A$ 是不定的。

---

## 9.7 二次型的几何意义

<div class="context-flow" markdown>

$\mathbf{x}^TA\mathbf{x} = c$ 定义二次曲面 → 正交替换沿**特征向量方向**（主轴）消去交叉项 → 曲面类型由签名 $(p,q)$ 决定

</div>

二次型在几何上描述了二次曲线和二次曲面。

!!! definition "定义 9.8 (二次曲面)"
    $\mathbb{R}^n$ 中的方程 $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x} = c$（$c$ 为常数）定义的集合称为**二次曲面**（quadric surface）。在 $\mathbb{R}^2$ 中为**二次曲线**（conic section）。

!!! theorem "定理 9.8 (主轴定理)"
    设 $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ 是实二次型，$A$ 为 $n$ 阶实对称矩阵。通过正交替换 $\mathbf{x} = P\mathbf{y}$（$P$ 的列为 $A$ 的标准正交特征向量），二次型化为标准形

    $$Q = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \cdots + \lambda_n y_n^2$$

    新坐标轴 $y_1, y_2, \ldots, y_n$ 的方向（即 $P$ 的列向量）称为二次曲面的**主轴**（principal axes）。

??? proof "证明"
    直接由谱定理和正交替换得到。正交替换相当于旋转坐标系，使得新坐标轴沿特征向量方向。在新坐标下，二次型没有交叉项，从而二次曲面的方程取最简形式。$\blacksquare$

!!! example "例 9.9"
    分类 $\mathbb{R}^2$ 中的二次曲线 $5x_1^2 + 4x_1x_2 + 8x_2^2 = 36$。

    由例 9.4，正交替换后得 $4y_1^2 + 9y_2^2 = 36$，即 $\dfrac{y_1^2}{9} + \dfrac{y_2^2}{4} = 1$。这是以 $y_1$ 轴和 $y_2$ 轴为主轴的**椭圆**，半长轴 $a = 3$，半短轴 $b = 2$。

    主轴方向为 $A$ 的特征向量 $\mathbf{v}_1 = \frac{1}{\sqrt{5}}(-2, 1)^T$（对应 $\lambda = 4$）和 $\mathbf{v}_2 = \frac{1}{\sqrt{5}}(1, 2)^T$（对应 $\lambda = 9$）。

!!! example "例 9.10"
    $\mathbb{R}^3$ 中二次曲面的标准分类（设 $Q = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \lambda_3 y_3^2 = 1$）：

    - **椭球面**：$\lambda_1, \lambda_2, \lambda_3 > 0$，如 $\frac{y_1^2}{a^2} + \frac{y_2^2}{b^2} + \frac{y_3^2}{c^2} = 1$；
    - **单叶双曲面**：两正一负，如 $\frac{y_1^2}{a^2} + \frac{y_2^2}{b^2} - \frac{y_3^2}{c^2} = 1$；
    - **双叶双曲面**：一正两负，如 $\frac{y_1^2}{a^2} - \frac{y_2^2}{b^2} - \frac{y_3^2}{c^2} = 1$；
    - 若 $Q = 0$（齐次情形）则为**二次锥面**。

    二次曲面的类型由二次型的**签名** $(p, q)$ 决定。

---

## 9.8 双线性型

<div class="context-flow" markdown>

**前置**：二次型 $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ 是"一元函数" → 双线性型 $f(\mathbf{x}, \mathbf{y}) = \mathbf{x}^TA\mathbf{y}$ 是"二元函数" → 通过**极化恒等式**二次型可恢复出双线性型 → 对称/反对称双线性型引出正交群与辛群

</div>

二次型是"对角线上"的值 $Q(\mathbf{x}) = f(\mathbf{x}, \mathbf{x})$。要全面理解二次型，必须先理解更一般的双线性型。双线性型是线性代数中最基本的"两个向量之间的标量函数"，内积、行列式、面积形式等都是它的特例。

!!! definition "定义 9.9 (双线性型)"
    设 $V$ 是域 $\mathbb{F}$ 上的 $n$ 维向量空间。映射 $f: V \times V \to \mathbb{F}$ 称为 $V$ 上的**双线性型**（bilinear form），如果 $f$ 对每个变量都是线性的：

    - 对第一个变量：$f(\alpha \mathbf{x}_1 + \beta \mathbf{x}_2, \mathbf{y}) = \alpha f(\mathbf{x}_1, \mathbf{y}) + \beta f(\mathbf{x}_2, \mathbf{y})$
    - 对第二个变量：$f(\mathbf{x}, \alpha \mathbf{y}_1 + \beta \mathbf{y}_2) = \alpha f(\mathbf{x}, \mathbf{y}_1) + \beta f(\mathbf{x}, \mathbf{y}_2)$

    对所有 $\mathbf{x}, \mathbf{x}_1, \mathbf{x}_2, \mathbf{y}, \mathbf{y}_1, \mathbf{y}_2 \in V$ 和 $\alpha, \beta \in \mathbb{F}$ 成立。

!!! definition "定义 9.10 (双线性型的矩阵)"
    设 $\mathcal{B} = \{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 是 $V$ 的一组基。双线性型 $f$ 在基 $\mathcal{B}$ 下的**矩阵**（又称 **Gram 矩阵**）定义为

    $$
    A = (a_{ij})_{n \times n}, \quad a_{ij} = f(\mathbf{e}_i, \mathbf{e}_j)
    $$

    若 $\mathbf{x} = \sum x_i \mathbf{e}_i$，$\mathbf{y} = \sum y_j \mathbf{e}_j$，则

    $$
    f(\mathbf{x}, \mathbf{y}) = \mathbf{x}^T A \mathbf{y} = \sum_{i=1}^n \sum_{j=1}^n a_{ij} x_i y_j
    $$

!!! theorem "定理 9.9 (基变换下的矩阵变换)"
    设双线性型 $f$ 在基 $\mathcal{B}$ 下的矩阵为 $A$，在基 $\mathcal{B}'$ 下的矩阵为 $A'$。若从 $\mathcal{B}$ 到 $\mathcal{B}'$ 的过渡矩阵为 $C$，则

    $$
    A' = C^T A C
    $$

??? proof "证明"
    设 $\mathcal{B}' = \{\mathbf{e}_1', \ldots, \mathbf{e}_n'\}$，$(\mathbf{e}_1', \ldots, \mathbf{e}_n') = (\mathbf{e}_1, \ldots, \mathbf{e}_n)C$。
    则

    $$
    a_{ij}' = f(\mathbf{e}_i', \mathbf{e}_j') = f\!\left(\sum_k c_{ki}\mathbf{e}_k, \sum_l c_{lj}\mathbf{e}_l\right) = \sum_{k,l} c_{ki} a_{kl} c_{lj}
    $$

    写成矩阵形式即 $A' = C^TAC$。$\blacksquare$

!!! definition "定义 9.11 (双线性型的秩与非退化性)"
    双线性型 $f$ 的**秩**定义为其 Gram 矩阵的秩（由定理 9.9，这与基的选取无关）。

    若 $\operatorname{rank}(f) = n$（即 Gram 矩阵可逆），则称 $f$ 是**非退化的**（nondegenerate）。等价地，$f$ 非退化当且仅当：若对所有 $\mathbf{y} \in V$ 有 $f(\mathbf{x}, \mathbf{y}) = 0$，则 $\mathbf{x} = \mathbf{0}$。

!!! theorem "定理 9.10 (极化恒等式)"
    设 $\operatorname{char}(\mathbb{F}) \neq 2$。对称双线性型 $f$ 与其关联的二次型 $Q(\mathbf{x}) = f(\mathbf{x}, \mathbf{x})$ 通过**极化恒等式**相互确定：

    $$
    f(\mathbf{x}, \mathbf{y}) = \frac{1}{2}\bigl[Q(\mathbf{x} + \mathbf{y}) - Q(\mathbf{x}) - Q(\mathbf{y})\bigr]
    $$

??? proof "证明"
    展开 $Q(\mathbf{x}+\mathbf{y}) = f(\mathbf{x}+\mathbf{y}, \mathbf{x}+\mathbf{y})$：

    $$
    Q(\mathbf{x}+\mathbf{y}) = f(\mathbf{x}, \mathbf{x}) + f(\mathbf{x}, \mathbf{y}) + f(\mathbf{y}, \mathbf{x}) + f(\mathbf{y}, \mathbf{y})
    $$

    由对称性 $f(\mathbf{x}, \mathbf{y}) = f(\mathbf{y}, \mathbf{x})$，故 $Q(\mathbf{x}+\mathbf{y}) = Q(\mathbf{x}) + 2f(\mathbf{x}, \mathbf{y}) + Q(\mathbf{y})$。解出 $f(\mathbf{x}, \mathbf{y})$ 即得。$\blacksquare$

!!! definition "定义 9.12 (对称与反对称双线性型)"
    双线性型 $f$ 称为

    - **对称的**（symmetric）：若 $f(\mathbf{x}, \mathbf{y}) = f(\mathbf{y}, \mathbf{x})$，对所有 $\mathbf{x}, \mathbf{y}$；
    - **反对称的**（antisymmetric / skew-symmetric）：若 $f(\mathbf{x}, \mathbf{y}) = -f(\mathbf{y}, \mathbf{x})$，对所有 $\mathbf{x}, \mathbf{y}$。

    对称双线性型的矩阵满足 $A = A^T$；反对称双线性型的矩阵满足 $A = -A^T$（且对角元为零）。

!!! theorem "定理 9.11 (双线性型的分解)"
    设 $\operatorname{char}(\mathbb{F}) \neq 2$。任意双线性型 $f$ 可唯一分解为对称部分与反对称部分之和：

    $$
    f = f_s + f_a, \quad f_s(\mathbf{x},\mathbf{y}) = \frac{f(\mathbf{x},\mathbf{y}) + f(\mathbf{y},\mathbf{x})}{2}, \quad f_a(\mathbf{x},\mathbf{y}) = \frac{f(\mathbf{x},\mathbf{y}) - f(\mathbf{y},\mathbf{x})}{2}
    $$

??? proof "证明"
    直接验证 $f_s$ 对称、$f_a$ 反对称、$f = f_s + f_a$。唯一性：若 $f = g_s + g_a$（$g_s$ 对称，$g_a$ 反对称），则 $f(\mathbf{x},\mathbf{y}) + f(\mathbf{y},\mathbf{x}) = 2g_s(\mathbf{x},\mathbf{y})$，故 $g_s = f_s$，从而 $g_a = f_a$。$\blacksquare$

!!! theorem "定理 9.12 (对称双线性型的标准形)"
    设 $\operatorname{char}(\mathbb{F}) \neq 2$，$f$ 是 $V$ 上的对称双线性型。则存在 $V$ 的一组基 $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 使得

    $$
    f(\mathbf{e}_i, \mathbf{e}_j) = 0 \quad (i \neq j)
    $$

    即 $f$ 在该基下的矩阵为对角矩阵。

??? proof "证明"
    对 $n$ 用归纳法。$n = 1$ 时显然。设对 $n-1$ 维空间成立。

    若 $f \equiv 0$，取任意基即可。否则，存在 $\mathbf{v}$ 使 $f(\mathbf{v}, \mathbf{v}) \neq 0$（若所有 $f(\mathbf{v}, \mathbf{v}) = 0$，由极化恒等式 $f \equiv 0$，矛盾）。令 $\mathbf{e}_1 = \mathbf{v}$，$d = f(\mathbf{e}_1, \mathbf{e}_1) \neq 0$。

    令 $W = \{\mathbf{w} \in V : f(\mathbf{e}_1, \mathbf{w}) = 0\}$。则 $V = \operatorname{span}\{\mathbf{e}_1\} \oplus W$（对任意 $\mathbf{x} \in V$，令 $\mathbf{w} = \mathbf{x} - \frac{f(\mathbf{e}_1, \mathbf{x})}{d}\mathbf{e}_1$，则 $f(\mathbf{e}_1, \mathbf{w}) = 0$）。

    由归纳假设，$f|_W$ 可在某基下对角化。合并 $\mathbf{e}_1$ 和 $W$ 的基即得结论。$\blacksquare$

!!! definition "定义 9.13 (正交补空间)"
    设 $f$ 是 $V$ 上的双线性型，$S \subseteq V$。$S$ 关于 $f$ 的**左正交补**和**右正交补**分别为

    $$
    S^{\perp_L} = \{\mathbf{x} \in V : f(\mathbf{x}, \mathbf{s}) = 0, \forall \mathbf{s} \in S\}, \quad S^{\perp_R} = \{\mathbf{y} \in V : f(\mathbf{s}, \mathbf{y}) = 0, \forall \mathbf{s} \in S\}
    $$

    若 $f$ 对称，则 $S^{\perp_L} = S^{\perp_R}$，简记为 $S^\perp$。$V$ 的**根**（radical）定义为 $\operatorname{rad}(f) = V^\perp$。$f$ 非退化当且仅当 $\operatorname{rad}(f) = \{\mathbf{0}\}$。

!!! example "例 9.11"
    $\mathbb{R}^3$ 上双线性型 $f(\mathbf{x}, \mathbf{y}) = x_1y_1 + x_1y_2 + x_2y_1 + 3x_2y_2 - x_3y_3$ 的 Gram 矩阵为

    $$
    A = \begin{pmatrix} 1 & 1 & 0 \\ 1 & 3 & 0 \\ 0 & 0 & -1 \end{pmatrix}
    $$

    $\det(A) = (3-1)(-1) = -2 \neq 0$，故 $f$ 非退化。$A$ 对称，故 $f$ 是对称双线性型。

!!! example "例 9.12"
    $\mathbb{R}^3$ 上的反对称双线性型 $f(\mathbf{x}, \mathbf{y}) = x_1y_2 - x_2y_1$，矩阵为

    $$
    A = \begin{pmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
    $$

    $\operatorname{rank}(A) = 2$，$\operatorname{rad}(f) = \operatorname{span}\{\mathbf{e}_3\}$，$f$ 退化。

!!! example "例 9.13"
    验证极化恒等式。设 $Q(x_1,x_2) = 2x_1^2 + 3x_1x_2 + x_2^2$，对应对称矩阵 $A = \begin{pmatrix} 2 & 3/2 \\ 3/2 & 1 \end{pmatrix}$。

    取 $\mathbf{x} = (1,0)^T$，$\mathbf{y} = (0,1)^T$：

    $$
    Q(\mathbf{x}+\mathbf{y}) = Q(1,1) = 2 + 3 + 1 = 6, \quad Q(\mathbf{x}) = 2, \quad Q(\mathbf{y}) = 1
    $$

    $$
    f(\mathbf{x},\mathbf{y}) = \frac{6-2-1}{2} = \frac{3}{2} = a_{12} \checkmark
    $$

!!! theorem "定理 9.13 (反对称双线性型的标准形)"
    设 $f$ 是有限维向量空间 $V$ 上的反对称双线性型。则存在基 $\{\mathbf{e}_1, \mathbf{f}_1, \ldots, \mathbf{e}_m, \mathbf{f}_m, \mathbf{g}_1, \ldots, \mathbf{g}_k\}$ 使得

    $$
    f(\mathbf{e}_i, \mathbf{f}_j) = \delta_{ij}, \quad f(\mathbf{e}_i, \mathbf{e}_j) = f(\mathbf{f}_i, \mathbf{f}_j) = 0
    $$

    且 $f(\mathbf{g}_s, \cdot) \equiv 0$。特别地，反对称双线性型的秩必为偶数 $2m$。

??? proof "证明"
    若 $f \equiv 0$，结论显然。否则存在 $\mathbf{x}, \mathbf{y}$ 使 $f(\mathbf{x}, \mathbf{y}) = c \neq 0$。令 $\mathbf{e}_1 = \mathbf{x}$，$\mathbf{f}_1 = \mathbf{y}/c$，则 $f(\mathbf{e}_1, \mathbf{f}_1) = 1$。

    令 $W = \{\mathbf{w} : f(\mathbf{e}_1, \mathbf{w}) = 0 \text{ 且 } f(\mathbf{f}_1, \mathbf{w}) = 0\}$。可验证 $V = \operatorname{span}\{\mathbf{e}_1, \mathbf{f}_1\} \oplus W$（对任意 $\mathbf{v}$，令 $\mathbf{w} = \mathbf{v} - f(\mathbf{v}, \mathbf{f}_1)\mathbf{e}_1 + f(\mathbf{v}, \mathbf{e}_1)\mathbf{f}_1$，验证 $\mathbf{w} \in W$）。

    对 $W$ 上的 $f|_W$ 递归应用即得。$\blacksquare$

---

## 9.9 辛空间

<div class="context-flow" markdown>

反对称非退化双线性型 = **辛形式** → 辛空间必偶维 → **Darboux 定理**：所有同维辛空间等价 → 辛几何是 Hamilton 力学的数学语言

</div>

辛空间（symplectic space）是装备了非退化反对称双线性型的向量空间，它在经典力学（Hamilton 系统）、量子力学和现代微分几何中有着基础性的地位。

!!! definition "定义 9.14 (辛形式与辛空间)"
    设 $V$ 是域 $\mathbb{F}$（$\operatorname{char}(\mathbb{F}) \neq 2$）上的有限维向量空间。$V$ 上的**辛形式**（symplectic form）是一个非退化的反对称双线性型 $\omega: V \times V \to \mathbb{F}$。配备辛形式的向量空间 $(V, \omega)$ 称为**辛空间**（symplectic space）。

!!! theorem "定理 9.14 (辛空间的维数)"
    辛空间的维数必为偶数。

??? proof "证明"
    设 $(V, \omega)$ 是辛空间，$\dim V = n$。$\omega$ 在任意基下的矩阵 $A$ 满足 $A^T = -A$（反对称），故

    $$
    \det(A) = \det(A^T) = \det(-A) = (-1)^n \det(A)
    $$

    由于 $\omega$ 非退化，$\det(A) \neq 0$，故 $(-1)^n = 1$，即 $n$ 为偶数。$\blacksquare$

!!! definition "定义 9.15 (辛基)"
    设 $(V, \omega)$ 是 $2n$ 维辛空间。$V$ 的一组基 $\{\mathbf{e}_1, \ldots, \mathbf{e}_n, \mathbf{f}_1, \ldots, \mathbf{f}_n\}$ 称为**辛基**（symplectic basis / Darboux basis），如果

    $$
    \omega(\mathbf{e}_i, \mathbf{e}_j) = 0, \quad \omega(\mathbf{f}_i, \mathbf{f}_j) = 0, \quad \omega(\mathbf{e}_i, \mathbf{f}_j) = \delta_{ij}
    $$

    在辛基下，$\omega$ 的矩阵为标准辛矩阵

    $$
    J_{2n} = \begin{pmatrix} O_n & I_n \\ -I_n & O_n \end{pmatrix}
    $$

!!! theorem "定理 9.15 (Darboux 定理——辛空间的标准形)"
    每个辛空间 $(V, \omega)$ 都存在辛基。因此所有 $2n$ 维辛空间（在同一域上）彼此同构。

??? proof "证明"
    这是定理 9.13 在非退化情形（$k = 0$）的直接推论。由于 $\omega$ 非退化，秩为 $\dim V = 2m$，定理 9.13 给出基 $\{\mathbf{e}_1, \mathbf{f}_1, \ldots, \mathbf{e}_m, \mathbf{f}_m\}$ 满足 $\omega(\mathbf{e}_i, \mathbf{f}_j) = \delta_{ij}$，$\omega(\mathbf{e}_i, \mathbf{e}_j) = \omega(\mathbf{f}_i, \mathbf{f}_j) = 0$。这恰是辛基。$\blacksquare$

!!! definition "定义 9.16 (辛矩阵与辛群)"
    $2n$ 阶实方阵 $M$ 称为**辛矩阵**（symplectic matrix），如果

    $$
    M^T J_{2n} M = J_{2n}
    $$

    所有 $2n$ 阶辛矩阵在矩阵乘法下构成群，称为**辛群**（symplectic group），记作 $\operatorname{Sp}(2n, \mathbb{F})$。

!!! theorem "定理 9.16 (辛矩阵的行列式)"
    辛矩阵的行列式为 $1$。

??? proof "证明"
    由 $M^T J M = J$，取行列式得 $\det(M)^2 \det(J) = \det(J)$。由于 $\det(J) = 1$（可直接计算或利用 $J^2 = -I$ 得 $\det(J)^2 = 1$），故 $\det(M)^2 = 1$，$\det(M) = \pm 1$。

    为证 $\det(M) = 1$，注意辛群 $\operatorname{Sp}(2n, \mathbb{R})$ 是连通的（可以连续地将 $M$ 变形为 $I_{2n}$），而行列式是连续函数，$\det(I_{2n}) = 1$，故 $\det(M) = 1$。

    另一个代数证明：将 $M$ 写成分块 $M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$，由 $M^TJM = J$ 得 $A^TD - C^TB = I$，因此 $\det(M)$ 的正负性可由 Pfaffian 论证确定为 $+1$。$\blacksquare$

!!! definition "定义 9.17 (Lagrange 子空间)"
    设 $(V, \omega)$ 是 $2n$ 维辛空间。子空间 $L \subseteq V$ 称为 **Lagrange 子空间**（Lagrangian subspace），如果 $L = L^\perp$（关于 $\omega$），即

    $$
    \omega(\mathbf{x}, \mathbf{y}) = 0, \quad \forall \mathbf{x}, \mathbf{y} \in L
    $$

    且 $\dim L = n$（达到最大各向同性子空间的维数）。

!!! theorem "定理 9.17 (Lagrange 子空间的维数)"
    设 $(V, \omega)$ 是 $2n$ 维辛空间，$L$ 是各向同性子空间（即 $\omega|_{L \times L} = 0$）。则 $\dim L \leq n$，且等号成立当且仅当 $L$ 是 Lagrange 子空间。

??? proof "证明"
    令 $L^\perp = \{\mathbf{v} \in V : \omega(\mathbf{v}, \mathbf{l}) = 0, \forall \mathbf{l} \in L\}$。由 $\omega$ 非退化，映射 $V \to V^*$，$\mathbf{v} \mapsto \omega(\mathbf{v}, \cdot)$ 是同构，故 $\dim L^\perp = 2n - \dim L$。

    各向同性意味着 $L \subseteq L^\perp$，故 $\dim L \leq \dim L^\perp = 2n - \dim L$，即 $\dim L \leq n$。

    等号成立当且仅当 $L = L^\perp$，即 $L$ 是 Lagrange 子空间。$\blacksquare$

!!! example "例 9.14"
    标准辛空间 $(\mathbb{R}^4, \omega)$，辛形式由矩阵 $J_4 = \begin{pmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ -1 & 0 & 0 & 0 \\ 0 & -1 & 0 & 0 \end{pmatrix}$ 定义。

    标准辛基为 $\mathbf{e}_1 = (1,0,0,0)^T$，$\mathbf{e}_2 = (0,1,0,0)^T$，$\mathbf{f}_1 = (0,0,1,0)^T$，$\mathbf{f}_2 = (0,0,0,1)^T$。

    $L_1 = \operatorname{span}\{\mathbf{e}_1, \mathbf{e}_2\}$ 是 Lagrange 子空间：$\omega(\mathbf{e}_1, \mathbf{e}_2) = 0$，$\dim L_1 = 2 = n$。

    $L_2 = \operatorname{span}\{\mathbf{f}_1, \mathbf{f}_2\}$ 也是 Lagrange 子空间。

!!! example "例 9.15"
    辛矩阵的例子。矩阵 $M = \begin{pmatrix} a & 0 & 0 & b \\ 0 & d & -c & 0 \\ 0 & c & d & 0 \\ -b & 0 & 0 & a \end{pmatrix}$（其中 $a^2 + b^2 = 1$，$c^2 + d^2 = 1$）是 $4 \times 4$ 辛矩阵。验证 $M^TJ_4M = J_4$ 可直接分块计算。$\det(M) = (a^2+b^2)(c^2+d^2) = 1$。

!!! example "例 9.16"
    Hamilton 力学中的应用。经典力学的相空间 $\mathbb{R}^{2n}$ 以广义坐标 $(q_1, \ldots, q_n)$ 和广义动量 $(p_1, \ldots, p_n)$ 为坐标。辛形式为

    $$
    \omega = \sum_{i=1}^n dp_i \wedge dq_i
    $$

    Hamilton 运动方程 $\dot{q}_i = \frac{\partial H}{\partial p_i}$，$\dot{p}_i = -\frac{\partial H}{\partial q_i}$ 可紧凑地写为

    $$
    \dot{\mathbf{z}} = J_{2n} \nabla H(\mathbf{z})
    $$

    其中 $\mathbf{z} = (q_1, \ldots, q_n, p_1, \ldots, p_n)^T$。正则变换（保持 Hamilton 方程形式的变量替换）恰对应辛矩阵。

---

## 9.10 Hermite 型

<div class="context-flow" markdown>

**实** $\to$ **复**：内积从双线性变为**半双线性** → Hermite 型 $H(\mathbf{x}) = \mathbf{x}^*A\mathbf{x}$（$A = A^*$ Hermite 矩阵）→ **惯性定理**在复数域上的推广 → 签名仍是完全不变量

</div>

当基础域从 $\mathbb{R}$ 扩展到 $\mathbb{C}$ 时，对称双线性型的自然推广是 Hermite 型。Hermite 型在量子力学（可观测量对应 Hermite 算子）和酉几何中有着根本性的作用。

!!! definition "定义 9.18 (半双线性型)"
    设 $V$ 是复向量空间。映射 $f: V \times V \to \mathbb{C}$ 称为**半双线性型**（sesquilinear form），如果

    - 对第一个变量是共轭线性的：$f(\alpha\mathbf{x}_1 + \beta\mathbf{x}_2, \mathbf{y}) = \bar{\alpha}f(\mathbf{x}_1, \mathbf{y}) + \bar{\beta}f(\mathbf{x}_2, \mathbf{y})$
    - 对第二个变量是线性的：$f(\mathbf{x}, \alpha\mathbf{y}_1 + \beta\mathbf{y}_2) = \alpha f(\mathbf{x}, \mathbf{y}_1) + \beta f(\mathbf{x}, \mathbf{y}_2)$

    （此处采用物理学惯例，第一个变量取共轭。数学文献中有时反过来。）

!!! definition "定义 9.19 (Hermite 型)"
    半双线性型 $f$ 称为 **Hermite 型**（Hermitian form），如果满足共轭对称性：

    $$
    f(\mathbf{x}, \mathbf{y}) = \overline{f(\mathbf{y}, \mathbf{x})}, \quad \forall \mathbf{x}, \mathbf{y} \in V
    $$

    特别地，$f(\mathbf{x}, \mathbf{x}) = \overline{f(\mathbf{x}, \mathbf{x})}$，故 $f(\mathbf{x}, \mathbf{x}) \in \mathbb{R}$。

!!! definition "定义 9.20 (Hermite 矩阵表示)"
    设 $\mathcal{B} = \{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 是 $V$ 的基，Hermite 型 $f$ 的矩阵 $A = (a_{ij})$，$a_{ij} = f(\mathbf{e}_i, \mathbf{e}_j)$。则

    $$
    f(\mathbf{x}, \mathbf{y}) = \mathbf{x}^* A \mathbf{y} = \sum_{i,j} \bar{x}_i a_{ij} y_j
    $$

    其中 $\mathbf{x}^* = \bar{\mathbf{x}}^T$ 是共轭转置。$A$ 满足 $A^* = A$（Hermite 矩阵）。

!!! theorem "定理 9.18 (Hermite 型的基变换)"
    设 Hermite 型 $f$ 在基 $\mathcal{B}$ 下矩阵为 $A$，在基 $\mathcal{B}'$ 下矩阵为 $A'$，过渡矩阵为 $C$，则

    $$
    A' = C^* A C
    $$

    即 $A'$ 与 $A$ **共轭合同**（$^*$-congruent）。

??? proof "证明"
    设 $\mathcal{B}' = \{\mathbf{e}_1', \ldots, \mathbf{e}_n'\}$，$(\mathbf{e}_1', \ldots, \mathbf{e}_n') = (\mathbf{e}_1, \ldots, \mathbf{e}_n)C$。则

    $$
    a_{ij}' = f(\mathbf{e}_i', \mathbf{e}_j') = f\!\left(\sum_k c_{ki}\mathbf{e}_k, \sum_l c_{lj}\mathbf{e}_l\right) = \sum_{k,l} \bar{c}_{ki} a_{kl} c_{lj}
    $$

    即 $A' = C^*AC$。$\blacksquare$

!!! theorem "定理 9.19 (Hermite 型的标准形)"
    设 $f$ 是复向量空间 $V$ 上的 Hermite 型。则存在 $V$ 的基使 $f$ 的矩阵为对角矩阵 $\operatorname{diag}(d_1, \ldots, d_n)$，其中 $d_i \in \mathbb{R}$。

    进一步，经过适当缩放，可化为

    $$
    \operatorname{diag}(\underbrace{1, \ldots, 1}_{p}, \underbrace{-1, \ldots, -1}_{q}, \underbrace{0, \ldots, 0}_{n-r})
    $$

??? proof "证明"
    类似实对称双线性型的配方法。若存在 $\mathbf{v}$ 使 $f(\mathbf{v}, \mathbf{v}) \neq 0$（此值为实数），令 $\mathbf{e}_1 = \mathbf{v}$，取正交补 $W = \{\mathbf{w} : f(\mathbf{e}_1, \mathbf{w}) = 0\}$，由归纳法对 $W$ 对角化。

    若对所有 $\mathbf{v}$ 有 $f(\mathbf{v}, \mathbf{v}) = 0$，由极化恒等式的 Hermite 版本

    $$
    f(\mathbf{x}, \mathbf{y}) = \frac{1}{4}\bigl[f(\mathbf{x}+\mathbf{y}, \mathbf{x}+\mathbf{y}) - f(\mathbf{x}-\mathbf{y}, \mathbf{x}-\mathbf{y}) + if(\mathbf{x}+i\mathbf{y}, \mathbf{x}+i\mathbf{y}) - if(\mathbf{x}-i\mathbf{y}, \mathbf{x}-i\mathbf{y})\bigr]
    $$

    可知 $f \equiv 0$。故归纳法可以进行。

    缩放：若 $d_k > 0$，令 $\mathbf{e}_k' = \mathbf{e}_k/\sqrt{d_k}$；若 $d_k < 0$，令 $\mathbf{e}_k' = \mathbf{e}_k/\sqrt{|d_k|}$。$\blacksquare$

!!! theorem "定理 9.20 (Hermite 型的惯性定律)"
    Hermite 型的标准形中正项个数 $p$ 和负项个数 $q$ 是不变量，不依赖于基的选取。$(p, q)$ 称为 Hermite 型的**签名**。

??? proof "证明"
    证明与实二次型的 Sylvester 惯性定律（定理 9.4）完全类似。假设两种标准形有不同的正项个数 $p > s$，构造子空间 $V_1$（维数 $p$，在其上 $f > 0$）和 $V_2$（维数 $n-s$，在其上 $f \leq 0$）。维数论证 $p + (n-s) > n$ 给出 $V_1 \cap V_2 \neq \{\mathbf{0}\}$，导出矛盾。$\blacksquare$

!!! theorem "定理 9.21 (Hermite 矩阵的谱定理)"
    Hermite 矩阵 $A$（$A^* = A$）的所有特征值都是实数，且 $A$ 可被酉矩阵对角化：存在酉矩阵 $U$ 使

    $$
    U^* A U = \operatorname{diag}(\lambda_1, \ldots, \lambda_n), \quad \lambda_i \in \mathbb{R}
    $$

??? proof "证明"
    **特征值为实数**：设 $A\mathbf{v} = \lambda\mathbf{v}$，$\mathbf{v} \neq \mathbf{0}$。则 $\lambda \mathbf{v}^*\mathbf{v} = \mathbf{v}^*A\mathbf{v} = (A\mathbf{v})^*\mathbf{v} = \bar{\lambda}\mathbf{v}^*\mathbf{v}$（利用 $A^* = A$），故 $\lambda = \bar{\lambda}$，$\lambda \in \mathbb{R}$。

    **酉对角化**：不同特征值的特征向量正交（证明同实情形）。对每个特征空间进行 Gram-Schmidt 正交化，合并得到酉矩阵 $U$。$\blacksquare$

!!! definition "定义 9.21 (正定 Hermite 型)"
    Hermite 型 $f$ 称为**正定的**，若 $f(\mathbf{x}, \mathbf{x}) > 0$ 对所有 $\mathbf{x} \neq \mathbf{0}$ 成立。等价条件与实情形类似：

    - Hermite 矩阵 $A$ 的所有特征值大于零；
    - 存在可逆矩阵 $C$ 使 $A = C^*C$；
    - $A$ 的所有顺序主子式大于零。

!!! example "例 9.17"
    Hermite 型 $f(\mathbf{x}, \mathbf{y}) = 2\bar{x}_1y_1 + (1+i)\bar{x}_1y_2 + (1-i)\bar{x}_2y_1 + 3\bar{x}_2y_2$ 的 Hermite 矩阵为

    $$
    A = \begin{pmatrix} 2 & 1+i \\ 1-i & 3 \end{pmatrix}
    $$

    验证 $A^* = A$：$\overline{(1+i)} = 1-i = a_{21}$ $\checkmark$。$\Delta_1 = 2 > 0$，$\Delta_2 = 6 - |1+i|^2 = 6 - 2 = 4 > 0$，故 $A$ 正定。

!!! example "例 9.18"
    找 Hermite 矩阵 $A = \begin{pmatrix} 1 & i \\ -i & 1 \end{pmatrix}$ 的特征值和酉对角化。

    特征多项式：$(1-\lambda)^2 - i(-i) = (1-\lambda)^2 - 1 = \lambda^2 - 2\lambda = \lambda(\lambda - 2)$。

    $\lambda_1 = 0$：$\mathbf{v}_1 = \frac{1}{\sqrt{2}}(i, 1)^T$。

    $\lambda_2 = 2$：$\mathbf{v}_2 = \frac{1}{\sqrt{2}}(-i, 1)^T$。

    酉矩阵 $U = \frac{1}{\sqrt{2}}\begin{pmatrix} i & -i \\ 1 & 1 \end{pmatrix}$，$U^*AU = \begin{pmatrix} 0 & 0 \\ 0 & 2 \end{pmatrix}$。

    签名 $(p, q) = (1, 0)$，秩 $r = 1$，$A$ 半正定但不正定。

!!! example "例 9.19"
    反 Hermite 型（anti-Hermitian / skew-Hermitian）。若半双线性型 $f$ 满足 $f(\mathbf{x}, \mathbf{y}) = -\overline{f(\mathbf{y}, \mathbf{x})}$，则 $f(\mathbf{x}, \mathbf{x})$ 纯虚。此时矩阵 $A$ 满足 $A^* = -A$（反 Hermite 矩阵），特征值全为纯虚数。

    例如 $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ 是反 Hermite 矩阵（也是实反对称矩阵），特征值为 $\pm i$。

## 练习题

1. **[配方] 配方法化二次型为标准形的几何意义是什么？**
   ??? success "参考答案"
       它是非正交的基变换，相当于对原空间进行扭曲、伸缩和平移，使得在新的斜坐标系下，各个坐标轴之间“正交”（即交叉项消失）。

2. **[正定] 一个 $2 \times 2$ 矩阵，对角元都是正数，它一定是正定矩阵吗？**
   ??? success "参考答案"
       不一定。例如 $\begin{pmatrix} 1 & 3 \\ 3 & 1 \end{pmatrix}$，尽管对角元为 1，但其行列式为 $1 - 9 = -8 < 0$，特征值为 4 和 -2，是不定的。

3. **[合同] 矩阵 $A$ 与 $B$ 合同（$B = C^TAC$）和相似（$B = P^{-1}AP$）有什么根本区别？**
   ??? success "参考答案"
       合同保持正负惯性指数（签名）不变，对应的几何是基变换下二次型系数矩阵的变化；相似保持特征值（谱）不变，对应的几何是基变换下线性映射矩阵的变化。正交相似（$C^T=C^{-1}$）则同时满足这两者。

4. **[惯性定理] 二次型的正负惯性指数之和等于什么？**
   ??? success "参考答案"
       等于该二次型的矩阵的秩 $\operatorname{rank}(A)$。

5. **[正定] 证明正定矩阵必可逆。**
   ??? success "参考答案"
       正定矩阵的所有特征值 $\lambda_i > 0$，因此行列式 $\det(A) = \prod \lambda_i > 0$，所以满秩且可逆。

6. **[双线性] 双线性型和内积有什么区别？**
   ??? success "参考答案"
       内积是一种特殊的对称双线性型，它额外要求满足“正定性”。一般的双线性型可以是退化的、不定的甚至是反对称的。

7. **[辛] 什么是辛空间？它的维数有什么限制？**
   ??? success "参考答案"
       装备了非退化、反对称双线性型（即辛形式）的向量空间。因为反对称非退化矩阵的行列式要求维数必须是偶数。

8. **[主子式] 判断 $A$ 负定（Negative Definite）的顺序主子式条件是什么？**
   ??? success "参考答案"
       奇数阶主子式为负，偶数阶主子式为正。即符号呈现 $-, +, -, + \dots$ 交替。

9. **[分解] 半双线性型中的“半”（sesqui-）指的是什么？**
   ??? success "参考答案"
       指其中一个变量（通常是第一个）是“共轭线性”的，提取复数标量时带有共轭符号，而另一个变量是严格线性的，加起来算“一阶半”线性。

10. **[爱因斯坦思考题] 狭义相对论的时空区间 $\Delta s^2 = -c^2\Delta t^2 + \Delta x^2 + \Delta y^2 + \Delta z^2$ 是一个怎样的二次型？它的“合同不变量”保证了什么？**
    ??? success "参考答案"
        它是一个**不定二次型**，签名是 $(3, 1)$ 或 $(1, 3)$。Sylvester 惯性定理保证了，无论使用哪种洛伦兹变换，这个（1个时间负号，3个空间正号）的签名永不改变，这在几何上确保了时间流向和因果律的拓扑结构不可被颠覆。

## 本章小结

本章将对称矩阵的研究扩展到二次型和双线性型，主要内容包括：

1. **二次型及其矩阵表示**：建立了一般实二次型与对称矩阵之间的一一对应。
2. **化标准形**：拉格朗日配方法提供了一种代数途径，而利用谱定理进行正交替换则提供了一种保持几何长度的“最自然”标准形。
3. **惯性定律**：证明了二次型在合同变换下的核心不变量——签名（正、负惯性指数），从而刻画了二次曲面的拓扑性质。
4. **正定性**：通过特征值、合同分解、顺序主子式等五个维度彻底澄清了正定/半正定矩阵的充要条件。
5. **双线性型与辛空间**：将纯代数讨论推广到了一般的对称和反对称形式，并简要介绍了经典力学的底层数学架构——辛空间。
