# 第 9 章 二次型

<div class="context-flow" markdown>

**前置**：Ch8 对称矩阵谱定理 · **本章脉络**：$\mathbf{x}^TA\mathbf{x}$ → 配方法/正交法化标准形 → **惯性定理**（签名不变） → 正定性判别 → 几何（椭球/双曲面）
本质：二次型是对称矩阵的"标量指纹"——签名 $(p,q)$ 完全决定等价类，正定性决定几何形状

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
