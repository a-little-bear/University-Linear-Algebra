# 第 2 章 矩阵与矩阵运算

<div class="context-flow" markdown>

**前置**：第 1 章增广矩阵 · 初等行变换 · **本章脉络**：矩阵定义 → 加法/数乘 → 矩阵乘法（$A\mathbf{x}=\mathbf{b}$ 的抽象） → 转置 → 逆矩阵 → 初等矩阵 → 分块矩阵 → 秩
**一句话本质**：将第 1 章的行变换操作升级为**代数对象**——矩阵本身成为可以运算的"元素"

</div>

矩阵（matrix）是线性代数的核心语言。它不仅是线性方程组的紧凑表示工具，更是描述线性变换的基本载体。本章将系统地介绍矩阵的基本概念和各种运算——加法、标量乘法、矩阵乘法、转置和求逆，并深入讨论初等矩阵、分块矩阵技术以及矩阵的秩等重要概念。

---

## 2.1 矩阵的定义与基本概念

<div class="context-flow" markdown>

**从方程组到矩阵**：第 1 章中系数排列成的矩形阵列 → 现在赋予它独立的数学身份 → 零矩阵、单位矩阵、对角矩阵等特殊角色

</div>

!!! definition "定义 2.1 (矩阵)"
    一个 $m \times n$ **矩阵**（matrix）是由 $m$ 行 $n$ 列元素排列成的矩形阵列，记为

    $$
    A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}
    $$

    简记为 $A = (a_{ij})_{m \times n}$ 或 $A = (a_{ij})$。元素 $a_{ij}$ 位于第 $i$ 行第 $j$ 列。全体 $m \times n$ 实矩阵的集合记为 $\mathbb{R}^{m \times n}$，全体 $m \times n$ 复矩阵的集合记为 $\mathbb{C}^{m \times n}$。

!!! definition "定义 2.2 (特殊矩阵)"
    设 $A = (a_{ij})$ 为矩阵，以下是几类重要的特殊矩阵：

    1. **零矩阵**（zero matrix）：所有元素为零的矩阵，记为 $O$ 或 $O_{m \times n}$。
    2. **方阵**（square matrix）：行数等于列数的矩阵，即 $m = n$，称为 $n$ 阶方阵。
    3. **单位矩阵**（identity matrix）：$n$ 阶方阵，主对角线元素全为 $1$，其余为 $0$，记为 $I_n$ 或 $I$：

        $$
        I_n = \begin{pmatrix} 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & 1 \end{pmatrix}
        $$

    4. **对角矩阵**（diagonal matrix）：$a_{ij} = 0$（$i \ne j$），记为 $\operatorname{diag}(d_1, d_2, \ldots, d_n)$。
    5. **上三角矩阵**（upper triangular matrix）：$a_{ij} = 0$（$i > j$），即主对角线下方的元素全为零。
    6. **下三角矩阵**（lower triangular matrix）：$a_{ij} = 0$（$i < j$），即主对角线上方的元素全为零。
    7. **对称矩阵**（symmetric matrix）：$a_{ij} = a_{ji}$（对所有 $i, j$），即 $A = A^T$。
    8. **反对称矩阵**（skew-symmetric / antisymmetric matrix）：$a_{ij} = -a_{ji}$，即 $A = -A^T$。特别地，反对称矩阵的对角线元素必为零。

---

## 2.2 矩阵的加法与标量乘法

<div class="context-flow" markdown>

**逐元素运算**：对应位置相加/数乘 → 满足 8 条公理 → $\mathbb{R}^{m \times n}$ 本身构成**向量空间**（→ 第 4 章的重要例子）

</div>

!!! definition "定义 2.3 (矩阵加法与标量乘法)"
    设 $A = (a_{ij})$ 和 $B = (b_{ij})$ 为同型矩阵（即都是 $m \times n$ 矩阵），$c$ 为标量，定义：

    1. **矩阵加法**：$A + B = (a_{ij} + b_{ij})$，即对应位置的元素相加。
    2. **标量乘法**：$cA = (ca_{ij})$，即每个元素乘以 $c$。

!!! theorem "定理 2.1 (矩阵加法与标量乘法的运算律)"
    设 $A, B, C$ 为同型矩阵，$c, d$ 为标量，则：

    1. $A + B = B + A$（加法交换律）
    2. $(A + B) + C = A + (B + C)$（加法结合律）
    3. $A + O = A$（零矩阵为加法单位元）
    4. $A + (-A) = O$（加法逆元）
    5. $c(A + B) = cA + cB$（标量乘法对矩阵加法的分配律）
    6. $(c + d)A = cA + dA$（标量加法对矩阵的分配律）
    7. $(cd)A = c(dA)$（标量乘法的结合律）
    8. $1 \cdot A = A$

??? proof "证明"
    以上性质均可由矩阵元素的实数运算律直接推得。以性质 5 为例：

    $(c(A + B))_{ij} = c(a_{ij} + b_{ij}) = ca_{ij} + cb_{ij} = (cA)_{ij} + (cB)_{ij} = (cA + cB)_{ij}$。

    由于对所有 $i, j$ 都成立，故 $c(A + B) = cA + cB$。其余类似。$\blacksquare$

!!! note "注"
    上述 8 条性质表明，全体 $m \times n$ 矩阵在矩阵加法和标量乘法下构成一个向量空间 $\mathbb{R}^{m \times n}$，其维数为 $mn$。

---

## 2.3 矩阵乘法

<div class="context-flow" markdown>

**核心运算**：$A\mathbf{x} = x_1\mathbf{a}_1 + \cdots + x_n\mathbf{a}_n$（列的线性组合） → 矩阵乘法 = 线性变换的复合（→ 第 5 章） → 注意：**不满足交换律**

</div>

!!! definition "定义 2.4 (矩阵乘法)"
    设 $A = (a_{ij})$ 为 $m \times p$ 矩阵，$B = (b_{jk})$ 为 $p \times n$ 矩阵，则 $A$ 与 $B$ 的**乘积**（matrix product）$C = AB$ 为 $m \times n$ 矩阵，其第 $(i, k)$ 元素为

    $$
    c_{ik} = \sum_{j=1}^{p} a_{ij} b_{jk} = a_{i1}b_{1k} + a_{i2}b_{2k} + \cdots + a_{ip}b_{pk}.
    $$

    即 $C$ 的第 $(i, k)$ 元素是 $A$ 的第 $i$ 行与 $B$ 的第 $k$ 列的**内积**。

!!! note "注"
    矩阵乘法 $AB$ 有意义的前提是 $A$ 的列数等于 $B$ 的行数。若 $A$ 为 $m \times p$ 矩阵、$B$ 为 $p \times n$ 矩阵，则 $AB$ 为 $m \times n$ 矩阵。

!!! theorem "定理 2.2 (矩阵乘法的运算律)"
    设矩阵的大小使得以下运算有意义，$c$ 为标量，则：

    1. $A(BC) = (AB)C$（结合律）
    2. $A(B + C) = AB + AC$（左分配律）
    3. $(A + B)C = AC + BC$（右分配律）
    4. $c(AB) = (cA)B = A(cB)$
    5. $I_m A = A = AI_n$（其中 $A$ 为 $m \times n$ 矩阵）

??? proof "证明"
    以结合律为例。设 $A = (a_{ij})_{m \times p}$，$B = (b_{jk})_{p \times q}$，$C = (c_{kl})_{q \times n}$。

    $(A(BC))_{il} = \sum_{j=1}^p a_{ij}(BC)_{jl} = \sum_{j=1}^p a_{ij}\left(\sum_{k=1}^q b_{jk}c_{kl}\right) = \sum_{j=1}^p \sum_{k=1}^q a_{ij}b_{jk}c_{kl}$

    $((AB)C)_{il} = \sum_{k=1}^q (AB)_{ik}c_{kl} = \sum_{k=1}^q \left(\sum_{j=1}^p a_{ij}b_{jk}\right)c_{kl} = \sum_{k=1}^q \sum_{j=1}^p a_{ij}b_{jk}c_{kl}$

    由有限求和可交换求和顺序，两者相等。$\blacksquare$

!!! theorem "定理 2.3 (矩阵乘法不满足交换律)"
    矩阵乘法一般**不满足交换律**，即 $AB \neq BA$（即使两个乘积都有意义）。

!!! example "例 2.1"
    设

    $$
    A = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}, \quad B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
    $$

    则

    $$
    AB = \begin{pmatrix} 2 & 1 \\ 1 & 0 \end{pmatrix}, \quad BA = \begin{pmatrix} 0 & 1 \\ 1 & 2 \end{pmatrix}
    $$

    显然 $AB \neq BA$。

!!! example "例 2.2"
    矩阵乘法的另一个与实数运算不同之处：$AB = O$ 不能推出 $A = O$ 或 $B = O$。例如

    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \quad B = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}, \quad AB = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix} = O.
    $$

!!! example "例 2.3"
    利用矩阵乘法，线性方程组 $A\mathbf{x} = \mathbf{b}$ 可以理解为：向量 $\mathbf{b}$ 是 $A$ 的列向量 $\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n$ 的线性组合。具体地，若 $A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)$，则

    $$
    A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n.
    $$

---

## 2.4 矩阵的转置

<div class="context-flow" markdown>

**行列互换**：$(AB)^T = B^TA^T$（顺序反转！） → 对称矩阵 $A = A^T$ 将在第 6 章谱定理和第 7 章正交性中扮演核心角色

</div>

!!! definition "定义 2.5 (转置矩阵)"
    设 $A = (a_{ij})$ 为 $m \times n$ 矩阵，$A$ 的**转置**（transpose）记为 $A^T$，是将 $A$ 的行变为列（或列变为行）所得的 $n \times m$ 矩阵，即 $(A^T)_{ij} = a_{ji}$。

!!! theorem "定理 2.4 (转置的性质)"
    设 $A, B$ 为适当大小的矩阵，$c$ 为标量，则：

    1. $(A^T)^T = A$
    2. $(A + B)^T = A^T + B^T$
    3. $(cA)^T = cA^T$
    4. $(AB)^T = B^T A^T$（注意顺序反转）

??? proof "证明"
    重点证明性质 4。设 $A$ 为 $m \times p$ 矩阵，$B$ 为 $p \times n$ 矩阵，则 $AB$ 为 $m \times n$ 矩阵，$(AB)^T$ 为 $n \times m$ 矩阵。

    $((AB)^T)_{ji} = (AB)_{ij} = \sum_{k=1}^p a_{ik}b_{kj}$

    $(B^T A^T)_{ji} = \sum_{k=1}^p (B^T)_{jk}(A^T)_{ki} = \sum_{k=1}^p b_{kj}a_{ik} = \sum_{k=1}^p a_{ik}b_{kj}$

    两式相等，故 $(AB)^T = B^T A^T$。$\blacksquare$

!!! proposition "命题 2.1 (对称矩阵与反对称矩阵的性质)"
    1. 若 $A$ 为方阵，则 $A + A^T$ 为对称矩阵，$A - A^T$ 为反对称矩阵。
    2. 任何方阵 $A$ 可以唯一地分解为一个对称矩阵与一个反对称矩阵之和：

    $$
    A = \frac{A + A^T}{2} + \frac{A - A^T}{2}.
    $$

??? proof "证明"
    1. $(A + A^T)^T = A^T + (A^T)^T = A^T + A = A + A^T$，故对称。$(A - A^T)^T = A^T - A = -(A - A^T)$，故反对称。

    2. 上式显然成立。唯一性：设 $A = S + K$（$S$ 对称，$K$ 反对称），则 $A^T = S - K$，故 $S = \frac{A + A^T}{2}$，$K = \frac{A - A^T}{2}$。$\blacksquare$

---

## 2.5 逆矩阵

<div class="context-flow" markdown>

**可逆性 = 方程组有唯一解**：$A$ 可逆 $\Leftrightarrow$ $A\mathbf{x}=\mathbf{b}$ 对任意 $\mathbf{b}$ 有唯一解 → 定理 2.7 给出 9 个等价条件，贯穿前 4 章的核心概念

</div>

!!! definition "定义 2.6 (可逆矩阵)"
    设 $A$ 为 $n$ 阶方阵，若存在 $n$ 阶方阵 $B$ 使得

    $$
    AB = BA = I_n,
    $$

    则称 $A$ 是**可逆的**（invertible）或**非奇异的**（nonsingular），$B$ 称为 $A$ 的**逆矩阵**（inverse matrix），记为 $A^{-1}$。若 $A$ 不可逆，则称 $A$ 为**奇异的**（singular）。

!!! theorem "定理 2.5 (逆矩阵的唯一性)"
    若矩阵 $A$ 可逆，则其逆矩阵唯一。

??? proof "证明"
    设 $B$ 和 $C$ 都是 $A$ 的逆矩阵，则

    $$
    B = BI = B(AC) = (BA)C = IC = C.
    $$

    因此 $B = C$，逆矩阵唯一。$\blacksquare$

!!! theorem "定理 2.6 (逆矩阵的性质)"
    设 $A, B$ 为 $n$ 阶可逆矩阵，$c$ 为非零标量，则：

    1. $(A^{-1})^{-1} = A$
    2. $(AB)^{-1} = B^{-1}A^{-1}$（注意顺序反转）
    3. $(A^T)^{-1} = (A^{-1})^T$
    4. $(cA)^{-1} = \frac{1}{c}A^{-1}$

??? proof "证明"
    以性质 2 为例：

    $(AB)(B^{-1}A^{-1}) = A(BB^{-1})A^{-1} = AIA^{-1} = AA^{-1} = I$

    类似可验证 $(B^{-1}A^{-1})(AB) = I$，因此 $(AB)^{-1} = B^{-1}A^{-1}$。$\blacksquare$

!!! theorem "定理 2.7 (可逆矩阵的等价条件)"
    设 $A$ 为 $n$ 阶方阵，以下条件等价：

    1. $A$ 可逆。
    2. $A$ 行等价于 $I_n$（$A$ 的简化行阶梯形为 $I_n$）。
    3. $A$ 有 $n$ 个主元位置。
    4. 齐次方程组 $A\mathbf{x} = \mathbf{0}$ 只有平凡解。
    5. 对任意 $\mathbf{b} \in \mathbb{R}^n$，方程组 $A\mathbf{x} = \mathbf{b}$ 有唯一解。
    6. $A$ 的列向量线性无关。
    7. $A$ 的列向量张成 $\mathbb{R}^n$。
    8. $\det(A) \neq 0$。
    9. $\operatorname{rank}(A) = n$。

### 求逆方法一：初等变换法

将 $A$ 与 $I_n$ 并列，构成 $n \times 2n$ 增广矩阵 $[A \mid I]$，对其施行行变换将左半部分化为 $I_n$。若成功，则右半部分即为 $A^{-1}$：

$$
[A \mid I] \xrightarrow{\text{行变换}} [I \mid A^{-1}].
$$

!!! example "例 2.4"
    求矩阵的逆：

    $$
    A = \begin{pmatrix} 1 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 3 \end{pmatrix}
    $$

    构造增广矩阵并进行行变换：

    $$
    \left(\begin{array}{ccc|ccc} 1 & 2 & 1 & 1 & 0 & 0 \\ 2 & 5 & 3 & 0 & 1 & 0 \\ 1 & 3 & 3 & 0 & 0 & 1 \end{array}\right)
    \xrightarrow{\text{行变换}}
    \left(\begin{array}{ccc|ccc} 1 & 0 & 0 & 6 & -3 & 1 \\ 0 & 1 & 0 & -3 & 2 & -1 \\ 0 & 0 & 1 & 1 & -1 & 1 \end{array}\right)
    $$

    因此

    $$
    A^{-1} = \begin{pmatrix} 6 & -3 & 1 \\ -3 & 2 & -1 \\ 1 & -1 & 1 \end{pmatrix}.
    $$

### 求逆方法二：伴随矩阵法

对于 $n$ 阶可逆矩阵 $A$，其逆矩阵也可以通过伴随矩阵来计算（将在行列式章节中详细介绍）：

$$
A^{-1} = \frac{1}{\det(A)} \operatorname{adj}(A),
$$

其中 $\operatorname{adj}(A)$ 是 $A$ 的**伴随矩阵**（adjugate matrix），即代数余子式矩阵的转置。

!!! example "例 2.5"
    对于 $2 \times 2$ 矩阵 $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$，若 $ad - bc \neq 0$，则

    $$
    A^{-1} = \frac{1}{ad - bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}.
    $$

---

## 2.6 初等矩阵

<div class="context-flow" markdown>

**行变换的矩阵化**：第 1 章的三种行变换 → 对 $I$ 做一次行变换得**初等矩阵** $E$ → 左乘 $E$ = 做一次行变换 → $A$ 可逆 $\Leftrightarrow$ $A$ 是初等矩阵之积

</div>

!!! definition "定义 2.7 (初等矩阵)"
    对单位矩阵 $I_n$ 施行**一次**初等行变换所得的矩阵称为**初等矩阵**（elementary matrix）。对应三种初等行变换，有三类初等矩阵：

    1. **行交换矩阵** $E(i, j)$：交换 $I_n$ 的第 $i$ 行和第 $j$ 行。
    2. **行倍乘矩阵** $E(i(c))$：将 $I_n$ 的第 $i$ 行乘以非零常数 $c$。
    3. **行倍加矩阵** $E(i, j(c))$：将 $I_n$ 的第 $j$ 行的 $c$ 倍加到第 $i$ 行上。

!!! theorem "定理 2.8 (初等矩阵与行变换)"
    对 $m \times n$ 矩阵 $A$ 施行一次初等行变换，等价于用相应的 $m$ 阶初等矩阵**左乘** $A$。

    类似地，对 $A$ 施行一次初等**列**变换，等价于用相应的 $n$ 阶初等矩阵**右乘** $A$。

??? proof "证明"
    以行倍加为例。设 $E = E(i, j(c))$，则 $E$ 的第 $i$ 行为 $\mathbf{e}_i + c\mathbf{e}_j$（$\mathbf{e}_k$ 表示第 $k$ 个标准基向量），其余行与 $I$ 相同。

    $(EA)$ 的第 $i$ 行 $= (\mathbf{e}_i + c\mathbf{e}_j)A = A$ 的第 $i$ 行 $+ c \cdot A$ 的第 $j$ 行。

    $(EA)$ 的第 $k$ 行（$k \neq i$）$= \mathbf{e}_k A = A$ 的第 $k$ 行。

    这正是行倍加 $R_i + cR_j \to R_i$ 的效果。$\blacksquare$

!!! theorem "定理 2.9 (初等矩阵可逆)"
    每个初等矩阵都是可逆的，且其逆矩阵仍为同类型的初等矩阵：

    1. $E(i, j)^{-1} = E(i, j)$
    2. $E(i(c))^{-1} = E(i(1/c))$
    3. $E(i, j(c))^{-1} = E(i, j(-c))$

!!! theorem "定理 2.10 (可逆矩阵是初等矩阵的乘积)"
    $n$ 阶矩阵 $A$ 可逆当且仅当 $A$ 可以表示为有限个初等矩阵的乘积。

??? proof "证明"
    $(\Rightarrow)$ 若 $A$ 可逆，则 $A$ 行等价于 $I_n$，即存在初等矩阵 $E_1, E_2, \ldots, E_k$ 使得 $E_k \cdots E_2 E_1 A = I$。因此 $A = E_1^{-1} E_2^{-1} \cdots E_k^{-1}$，而每个 $E_i^{-1}$ 也是初等矩阵。

    $(\Leftarrow)$ 初等矩阵都是可逆的，可逆矩阵的乘积仍可逆。$\blacksquare$

!!! example "例 2.6"
    将矩阵 $A = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}$ 表示为初等矩阵的乘积。

    将 $A$ 化为 $I$：

    $$
    \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix} \xrightarrow{R_2 - 3R_1} \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} \xrightarrow{R_1 - 2R_2} \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
    $$

    即 $E(1,2(-2)) \cdot E(2,1(-3)) \cdot A = I$，所以

    $$
    A = E(2,1(-3))^{-1} \cdot E(1,2(-2))^{-1} = E(2,1(3)) \cdot E(1,2(2)) = \begin{pmatrix} 1 & 0 \\ 3 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}.
    $$

---

## 2.7 分块矩阵

<div class="context-flow" markdown>

**分而治之**：大矩阵按块运算 → 分块对角矩阵使行列式/求逆分解为独立子问题 → 分块视角贯穿后续章节（第 6 章对角化、第 5 章不变子空间）

</div>

当矩阵规模较大时，将其分割为若干较小的子矩阵（**块**，block）往往能简化分析和计算。

!!! definition "定义 2.8 (分块矩阵)"
    将矩阵 $A$ 用水平线和垂直线分割成若干子矩阵，所得结果称为 $A$ 的一个**分块**（partitioning）。每个子矩阵称为 $A$ 的一个**块**（block）或**子块**。例如将 $m \times n$ 矩阵 $A$ 分块为

    $$
    A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}
    $$

    其中 $A_{11}$ 为 $m_1 \times n_1$ 矩阵，$A_{12}$ 为 $m_1 \times n_2$ 矩阵，等等，且 $m_1 + m_2 = m$，$n_1 + n_2 = n$。

!!! theorem "定理 2.11 (分块矩阵的乘法)"
    若矩阵 $A$ 和 $B$ 的分块方式使得对应的乘积有意义（即 $A$ 的列分法与 $B$ 的行分法一致），则分块矩阵的乘法可以像标量元素的矩阵乘法一样进行。

    例如，设

    $$
    A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}, \quad
    B = \begin{pmatrix} B_{11} & B_{12} \\ B_{21} & B_{22} \end{pmatrix}
    $$

    分块相容，则

    $$
    AB = \begin{pmatrix} A_{11}B_{11} + A_{12}B_{21} & A_{11}B_{12} + A_{12}B_{22} \\ A_{21}B_{11} + A_{22}B_{21} & A_{21}B_{12} + A_{22}B_{22} \end{pmatrix}.
    $$

!!! definition "定义 2.9 (分块对角矩阵)"
    形如

    $$
    A = \begin{pmatrix} A_1 & & \\ & A_2 & \\ & & \ddots & \\ & & & A_k \end{pmatrix}
    $$

    的矩阵称为**分块对角矩阵**（block diagonal matrix），其中 $A_1, A_2, \ldots, A_k$ 为方阵，对角块之外的元素全为零。记为 $A = \operatorname{diag}(A_1, A_2, \ldots, A_k)$。

!!! proposition "命题 2.2 (分块对角矩阵的性质)"
    设 $A = \operatorname{diag}(A_1, A_2, \ldots, A_k)$，则：

    1. $\det(A) = \det(A_1)\det(A_2)\cdots\det(A_k)$。
    2. $A$ 可逆当且仅当每个 $A_i$ 都可逆，此时 $A^{-1} = \operatorname{diag}(A_1^{-1}, A_2^{-1}, \ldots, A_k^{-1})$。

!!! example "例 2.7"
    设

    $$
    A = \begin{pmatrix} A_{11} & A_{12} \\ O & A_{22} \end{pmatrix}
    $$

    为分块上三角矩阵，其中 $A_{11}$ 和 $A_{22}$ 为方阵。若 $A_{11}$ 和 $A_{22}$ 都可逆，则 $A$ 可逆，且

    $$
    A^{-1} = \begin{pmatrix} A_{11}^{-1} & -A_{11}^{-1}A_{12}A_{22}^{-1} \\ O & A_{22}^{-1} \end{pmatrix}.
    $$

    验证：

    $$
    \begin{pmatrix} A_{11} & A_{12} \\ O & A_{22} \end{pmatrix}\begin{pmatrix} A_{11}^{-1} & -A_{11}^{-1}A_{12}A_{22}^{-1} \\ O & A_{22}^{-1} \end{pmatrix} = \begin{pmatrix} I & A_{12}A_{22}^{-1} - A_{12}A_{22}^{-1} \\ O & I \end{pmatrix} = \begin{pmatrix} I & O \\ O & I \end{pmatrix}.
    $$

---

## 2.8 矩阵的秩

<div class="context-flow" markdown>

**秩 = 矩阵的"有效维度"**：REF 中非零行数 = 列秩 = 行秩 → $\operatorname{rank}(A)$ 统一刻画了方程组的自由度 → Sylvester 不等式控制乘积的秩

</div>

!!! definition "定义 2.10 (矩阵的秩)"
    矩阵 $A$ 的**秩**（rank）记为 $\operatorname{rank}(A)$ 或 $r(A)$，定义为 $A$ 的行阶梯形中非零行的个数，即主元的个数。

    等价地，$\operatorname{rank}(A)$ 等于 $A$ 的列向量组的极大线性无关组的大小（**列秩**），也等于 $A$ 的行向量组的极大线性无关组的大小（**行秩**）。

!!! theorem "定理 2.12 (行秩 = 列秩)"
    对任意矩阵 $A$，其行秩等于列秩。

??? proof "证明"
    设 $A$ 为 $m \times n$ 矩阵，$\operatorname{rank}(A) = r$。将 $A$ 化为行阶梯形 $R$，$R$ 有 $r$ 个非零行，故行秩 $\le r$。由于行变换不改变行空间（行空间由行向量张成），$A$ 的行秩等于 $R$ 的行秩，等于 $r$（$R$ 的 $r$ 个非零行线性无关）。

    另一方面，行变换不改变列向量的线性相关性（因为不改变齐次方程 $A\mathbf{x} = \mathbf{0}$ 的解集），所以 $A$ 和 $R$ 有相同的列秩。$R$ 有 $r$ 个主元列，对应 $r$ 个线性无关的列向量，故列秩也等于 $r$。$\blacksquare$

!!! theorem "定理 2.13 (秩的性质)"
    设 $A$ 为 $m \times n$ 矩阵，$B$ 为 $n \times p$ 矩阵，则：

    1. $0 \le \operatorname{rank}(A) \le \min(m, n)$。
    2. $\operatorname{rank}(A) = \operatorname{rank}(A^T)$。
    3. $\operatorname{rank}(AB) \le \min(\operatorname{rank}(A), \operatorname{rank}(B))$。
    4. 若 $P$ 为 $m$ 阶可逆矩阵，$Q$ 为 $n$ 阶可逆矩阵，则 $\operatorname{rank}(PAQ) = \operatorname{rank}(A)$。
    5. $\operatorname{rank}(A + B) \le \operatorname{rank}(A) + \operatorname{rank}(B)$（当 $A, B$ 同型时）。

!!! theorem "定理 2.14 (Sylvester 秩不等式)"
    设 $A$ 为 $m \times n$ 矩阵，$B$ 为 $n \times p$ 矩阵，则

    $$
    \operatorname{rank}(A) + \operatorname{rank}(B) - n \le \operatorname{rank}(AB).
    $$

??? proof "证明"
    考虑 $B$ 的列空间和 $A$ 的零空间。$B$ 的列空间维数为 $\operatorname{rank}(B)$，$A$ 的零空间维数为 $n - \operatorname{rank}(A)$。

    $AB$ 的列空间由 $A$ 作用在 $B$ 的列空间上得到。$B$ 的列空间中被 $A$ 映射为零的向量恰好是 $B$ 的列空间与 $A$ 的零空间的交集。

    由维数公式：

    $$
    \operatorname{rank}(AB) = \dim(A(\operatorname{Col}(B))) = \operatorname{rank}(B) - \dim(\operatorname{Col}(B) \cap \ker(A))
    $$

    $$
    \ge \operatorname{rank}(B) - \dim(\ker(A)) = \operatorname{rank}(B) - (n - \operatorname{rank}(A))
    $$

    $$
    = \operatorname{rank}(A) + \operatorname{rank}(B) - n. \quad \blacksquare
    $$

!!! example "例 2.8"
    求矩阵的秩：

    $$
    A = \begin{pmatrix} 1 & 2 & 3 & 4 \\ 2 & 4 & 7 & 9 \\ 1 & 2 & 4 & 5 \end{pmatrix}
    $$

    行变换：

    $$
    \begin{pmatrix} 1 & 2 & 3 & 4 \\ 2 & 4 & 7 & 9 \\ 1 & 2 & 4 & 5 \end{pmatrix}
    \xrightarrow{R_2 - 2R_1, \; R_3 - R_1}
    \begin{pmatrix} 1 & 2 & 3 & 4 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 1 & 1 \end{pmatrix}
    \xrightarrow{R_3 - R_2}
    \begin{pmatrix} 1 & 2 & 3 & 4 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 0 \end{pmatrix}
    $$

    有 $2$ 个非零行，故 $\operatorname{rank}(A) = 2$。

!!! example "例 2.9"
    设 $A$ 为 $3 \times 5$ 矩阵，$B$ 为 $5 \times 4$ 矩阵，$\operatorname{rank}(A) = 3$，$\operatorname{rank}(B) = 4$。由 Sylvester 不等式：

    $$
    \operatorname{rank}(AB) \ge 3 + 4 - 5 = 2.
    $$

    又 $\operatorname{rank}(AB) \le \min(3, 4) = 3$，因此 $2 \le \operatorname{rank}(AB) \le 3$。
