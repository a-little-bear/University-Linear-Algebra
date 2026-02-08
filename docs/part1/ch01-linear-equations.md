# 第 1 章 线性方程组

线性方程组（system of linear equations）是线性代数最基本的研究对象之一。从初等数学中解二元一次方程组开始，我们就已经在接触线性方程组的求解问题。在本章中，我们将系统地建立线性方程组的理论框架：引入增广矩阵和初等行变换，发展高斯消元法这一强有力的算法工具，深入分析解的存在性、唯一性及其结构。这些内容构成了后续所有章节的基石。

---

## 1.1 线性方程与线性方程组

### 基本定义

!!! definition "定义 1.1 (线性方程)"
    含有 $n$ 个未知量 $x_1, x_2, \ldots, x_n$ 的**线性方程**（linear equation）是指形如

    $$
    a_1 x_1 + a_2 x_2 + \cdots + a_n x_n = b
    $$

    的方程，其中 $a_1, a_2, \ldots, a_n$ 和 $b$ 是已知的实数（或复数）常量，分别称为方程的**系数**（coefficient）和**常数项**（constant term）。

线性方程的关键特征在于：每个未知量只以一次幂出现，且未知量之间没有乘积项。例如 $3x_1 - 2x_2 + x_3 = 7$ 是线性方程，而 $x_1 x_2 + x_3 = 1$ 或 $\sin(x_1) + x_2 = 0$ 则不是。

!!! definition "定义 1.2 (线性方程组)"
    由 $m$ 个含 $n$ 个未知量的线性方程所构成的集合称为一个 $m \times n$ **线性方程组**（system of linear equations），一般形式为：

    $$
    \begin{cases}
    a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n = b_1 \\
    a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n = b_2 \\
    \quad \vdots \\
    a_{m1}x_1 + a_{m2}x_2 + \cdots + a_{mn}x_n = b_m
    \end{cases}
    $$

    其中 $a_{ij}$（$1 \le i \le m, 1 \le j \le n$）为系数，$b_i$ 为常数项。

!!! definition "定义 1.3 (解与解集)"
    使线性方程组中所有方程同时成立的一组数 $(s_1, s_2, \ldots, s_n)$ 称为该方程组的一个**解**（solution）。方程组所有解的集合称为**解集**（solution set）。若两个线性方程组具有相同的解集，则称它们是**等价的**（equivalent）。

### 几何解释

**二维情形（两个未知量）**：每个线性方程 $a_1 x_1 + a_2 x_2 = b$ 在 $\mathbb{R}^2$ 平面上表示一条直线。两个方程的解集对应两条直线的交集，有三种情形：

1. **唯一解**：两条直线相交于一点（两条直线斜率不同）。
2. **无穷多解**：两条直线重合（两方程成比例）。
3. **无解**：两条直线平行但不重合。

**三维情形（三个未知量）**：每个线性方程 $a_1 x_1 + a_2 x_2 + a_3 x_3 = b$ 在 $\mathbb{R}^3$ 中表示一个平面。三个方程的解集是三个平面的公共交集，可能是一个点、一条直线、一个平面，或者空集。

!!! example "例 1.1"
    考虑方程组

    $$
    \begin{cases} x_1 + x_2 = 3 \\ 2x_1 - x_2 = 0 \end{cases}
    $$

    第一个方程在平面上表示经过 $(3,0)$ 和 $(0,3)$ 的直线，第二个方程表示经过原点斜率为 $2$ 的直线。两条直线交于一点 $(1, 2)$，即唯一解 $x_1 = 1, x_2 = 2$。

!!! example "例 1.2"
    考虑三维方程组

    $$
    \begin{cases} x_1 + x_2 + x_3 = 1 \\ x_1 - x_2 + x_3 = 3 \end{cases}
    $$

    这两个方程各表示 $\mathbb{R}^3$ 中的一个平面。两个不平行的平面交于一条直线，因此该方程组有无穷多解。由两方程相减得 $2x_2 = -2$，即 $x_2 = -1$。代入第一个方程得 $x_1 + x_3 = 2$。设 $x_3 = t$（自由参数），则通解为 $x_1 = 2 - t,\; x_2 = -1,\; x_3 = t$，$t \in \mathbb{R}$。

---

## 1.2 增广矩阵与初等行变换

在求解线性方程组时，未知量的符号 $x_1, x_2, \ldots, x_n$ 仅起占位作用，真正决定解的是系数和常数项。因此，我们可以将方程组的全部信息压缩到一个矩阵中。

!!! definition "定义 1.4 (系数矩阵与增广矩阵)"
    对于线性方程组

    $$
    \begin{cases}
    a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n = b_1 \\
    a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n = b_2 \\
    \quad \vdots \\
    a_{m1}x_1 + a_{m2}x_2 + \cdots + a_{mn}x_n = b_m
    \end{cases}
    $$

    其**系数矩阵**（coefficient matrix）为

    $$
    A = \begin{pmatrix}
    a_{11} & a_{12} & \cdots & a_{1n} \\
    a_{21} & a_{22} & \cdots & a_{2n} \\
    \vdots & \vdots & \ddots & \vdots \\
    a_{m1} & a_{m2} & \cdots & a_{mn}
    \end{pmatrix}
    $$

    其**增广矩阵**（augmented matrix）为

    $$
    [A \mid \mathbf{b}] = \left(\begin{array}{cccc|c}
    a_{11} & a_{12} & \cdots & a_{1n} & b_1 \\
    a_{21} & a_{22} & \cdots & a_{2n} & b_2 \\
    \vdots & \vdots & \ddots & \vdots & \vdots \\
    a_{m1} & a_{m2} & \cdots & a_{mn} & b_m
    \end{array}\right)
    $$

方程组还可以紧凑地写成矩阵方程 $A\mathbf{x} = \mathbf{b}$，其中 $\mathbf{x} = (x_1, x_2, \ldots, x_n)^T$，$\mathbf{b} = (b_1, b_2, \ldots, b_m)^T$。

!!! definition "定义 1.5 (初等行变换)"
    对矩阵施行的以下三种操作称为**初等行变换**（elementary row operations）：

    1. **行交换**（row interchange）：交换矩阵的第 $i$ 行和第 $j$ 行，记作 $R_i \leftrightarrow R_j$。
    2. **行倍乘**（row scaling）：将第 $i$ 行乘以非零常数 $c$，记作 $cR_i \to R_i$。
    3. **行倍加**（row replacement）：将第 $j$ 行的 $c$ 倍加到第 $i$ 行上，记作 $R_i + cR_j \to R_i$。

!!! theorem "定理 1.1 (行变换保持解集不变)"
    对线性方程组的增广矩阵施行有限次初等行变换后得到的新方程组与原方程组等价，即它们具有相同的解集。

??? proof "证明"
    每种初等行变换都是可逆的：

    - 行交换 $R_i \leftrightarrow R_j$ 的逆操作仍是 $R_i \leftrightarrow R_j$。
    - 行倍乘 $cR_i \to R_i$（$c \neq 0$）的逆操作是 $\frac{1}{c}R_i \to R_i$。
    - 行倍加 $R_i + cR_j \to R_i$ 的逆操作是 $R_i - cR_j \to R_i$。

    因此，原方程组的任意解满足新方程组的所有方程，反之亦然。从而两个方程组的解集相同。$\blacksquare$

---

## 1.3 高斯消元法与高斯-Jordan 消元法

### 高斯消元法

**高斯消元法**（Gaussian elimination）是一种系统化的求解线性方程组的算法，其基本思想是通过初等行变换将增广矩阵化为行阶梯形，然后通过回代求解。

**算法步骤**：

1. 写出增广矩阵 $[A \mid \mathbf{b}]$。
2. 从最左列开始，找到该列中第一个非零元素所在的行（如果该列全为零，则移到下一列）。必要时进行行交换，将该非零元素移至当前处理行。该非零元素称为**主元**（pivot）。
3. 利用行倍加操作，将主元下方该列的所有元素消为零。
4. 对下一行和下一列重复步骤 2–3，直到所有行处理完毕。
5. 得到行阶梯形矩阵后，从最后一个含主元的行开始**回代**（back substitution），逐步求出各未知量。

!!! example "例 1.3"
    用高斯消元法解方程组

    $$
    \begin{cases} x_1 + 2x_2 + x_3 = 3 \\ 2x_1 + 5x_2 + 2x_3 = 7 \\ x_1 + 3x_2 + 3x_3 = 5 \end{cases}
    $$

    写出增广矩阵并进行行变换：

    $$
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    2 & 5 & 2 & 7 \\
    1 & 3 & 3 & 5
    \end{array}\right)
    \xrightarrow{R_2 - 2R_1 \to R_2}
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    1 & 3 & 3 & 5
    \end{array}\right)
    $$

    $$
    \xrightarrow{R_3 - R_1 \to R_3}
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    0 & 1 & 2 & 2
    \end{array}\right)
    \xrightarrow{R_3 - R_2 \to R_3}
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 2 & 1
    \end{array}\right)
    $$

    回代：由第三行 $2x_3 = 1$ 得 $x_3 = \frac{1}{2}$；由第二行 $x_2 = 1$；由第一行 $x_1 + 2(1) + \frac{1}{2} = 3$ 得 $x_1 = \frac{1}{2}$。

    解为 $x_1 = \frac{1}{2},\; x_2 = 1,\; x_3 = \frac{1}{2}$。

### 高斯-Jordan 消元法

**高斯-Jordan 消元法**（Gauss-Jordan elimination）在高斯消元的基础上进一步操作：不仅消去主元下方的元素，还消去主元上方的元素，并将每个主元化为 $1$，最终将增广矩阵化为**简化行阶梯形**。

!!! example "例 1.4"
    继续例 1.3 的结果，进行高斯-Jordan 消元：

    $$
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 2 & 1
    \end{array}\right)
    \xrightarrow{\frac{1}{2}R_3 \to R_3}
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 1 & \frac{1}{2}
    \end{array}\right)
    $$

    $$
    \xrightarrow{R_1 - R_3 \to R_1}
    \left(\begin{array}{ccc|c}
    1 & 2 & 0 & \frac{5}{2} \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 1 & \frac{1}{2}
    \end{array}\right)
    \xrightarrow{R_1 - 2R_2 \to R_1}
    \left(\begin{array}{ccc|c}
    1 & 0 & 0 & \frac{1}{2} \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 1 & \frac{1}{2}
    \end{array}\right)
    $$

    直接读出解 $x_1 = \frac{1}{2},\; x_2 = 1,\; x_3 = \frac{1}{2}$，无需回代。

---

## 1.4 行阶梯形与简化行阶梯形

!!! definition "定义 1.6 (行阶梯形矩阵)"
    矩阵处于**行阶梯形**（row echelon form, REF）当且仅当满足以下条件：

    1. 所有全零行位于矩阵的底部。
    2. 每个非零行的首个非零元素（称为**主元**或**领头元素**，leading entry / pivot）严格位于上一行主元的右侧。

!!! definition "定义 1.7 (简化行阶梯形矩阵)"
    矩阵处于**简化行阶梯形**（reduced row echelon form, RREF）当且仅当满足以下条件：

    1. 它是行阶梯形。
    2. 每个主元等于 $1$（称为**主 1**，leading 1）。
    3. 每个主元所在列的其余元素全为 $0$。

!!! example "例 1.5"
    以下矩阵处于行阶梯形但不是简化行阶梯形：

    $$
    \begin{pmatrix} 2 & 1 & 3 & 4 \\ 0 & 0 & 3 & 1 \\ 0 & 0 & 0 & 0 \end{pmatrix}
    $$

    以下矩阵处于简化行阶梯形：

    $$
    \begin{pmatrix} 1 & 0 & 0 & 2 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & 3 \end{pmatrix}, \qquad
    \begin{pmatrix} 1 & 3 & 0 & 5 \\ 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 \end{pmatrix}
    $$

!!! theorem "定理 1.2 (简化行阶梯形的唯一性)"
    每个矩阵都行等价于唯一的一个简化行阶梯形矩阵。

??? proof "证明"
    **存在性**：通过高斯-Jordan 消元法，任意矩阵都可以经过有限次初等行变换化为简化行阶梯形。

    **唯一性**：设矩阵 $A$ 行等价于两个简化行阶梯形矩阵 $R_1$ 和 $R_2$。由于行等价的矩阵对应同解方程组，因此 $R_1$ 和 $R_2$ 对应的齐次方程组具有完全相同的解集。

    我们证明 $R_1 = R_2$。首先，$R_1$ 和 $R_2$ 的主元列位置相同（因为主元列恰好对应基本变量，而基本变量由解集唯一确定）。其次，对于每个主元列，简化行阶梯形中该列只有主元位置为 $1$、其余为 $0$，所以对应列相同。对于自由变量列，其值由方程组的解集唯一确定。综合可得 $R_1 = R_2$。$\blacksquare$

!!! note "注"
    行阶梯形**不是**唯一的。同一个矩阵可以经过不同的初等行变换序列化为不同的行阶梯形。但简化行阶梯形是唯一的，这是一个重要的性质。

---

## 1.5 解的存在性与唯一性

!!! definition "定义 1.8 (相容与不相容)"
    若线性方程组有解（至少一个），则称该方程组是**相容的**（consistent）；若无解，则称为**不相容的**（inconsistent）。

!!! definition "定义 1.9 (主元位置与主元列)"
    矩阵中在行阶梯形中对应主元的位置称为**主元位置**（pivot position）。含有主元位置的列称为**主元列**（pivot column）。对应主元列的未知量称为**基本变量**（basic variable），其余未知量称为**自由变量**（free variable）。

!!! theorem "定理 1.3 (解的存在性)"
    线性方程组 $A\mathbf{x} = \mathbf{b}$ 相容的充要条件是：其增广矩阵 $[A \mid \mathbf{b}]$ 的行阶梯形中，$\mathbf{b}$ 所在列不是主元列。

    等价地，方程组相容当且仅当增广矩阵的行阶梯形中不出现形如 $(0\; 0\; \cdots\; 0 \mid c)$（其中 $c \neq 0$）的行。

??? proof "证明"
    若行阶梯形中出现 $(0\; 0\; \cdots\; 0 \mid c)$（$c \ne 0$），则对应方程为 $0 = c$，矛盾，方程组无解。

    若不出现这样的行，则每个非零行都至少有一个系数列中的主元，可以通过回代（将自由变量赋任意值）求得至少一个解，方程组相容。$\blacksquare$

!!! theorem "定理 1.4 (解的唯一性)"
    若相容的线性方程组没有自由变量（即每一列都是主元列），则方程组有唯一解；若存在自由变量，则方程组有无穷多解。

??? proof "证明"
    若没有自由变量，简化行阶梯形中每个未知量都是基本变量，对应列只有一个 $1$，解被唯一确定。

    若存在自由变量，令该自由变量取不同的实数值就得到不同的解，因此解有无穷多个。由于自由变量可取 $\mathbb{R}$ 中的任意值，解集是无穷的。$\blacksquare$

!!! proposition "命题 1.1"
    线性方程组的解的情况恰好有三种：**无解**、**唯一解**、**无穷多解**。不存在恰好有有限多个（大于1个）解的线性方程组。

??? proof "证明"
    假设方程组有两个不同的解 $\mathbf{x}_1$ 和 $\mathbf{x}_2$。对任意实数 $t$，令 $\mathbf{x}(t) = (1-t)\mathbf{x}_1 + t\mathbf{x}_2$。由于

    $$
    A\mathbf{x}(t) = (1-t)A\mathbf{x}_1 + tA\mathbf{x}_2 = (1-t)\mathbf{b} + t\mathbf{b} = \mathbf{b},
    $$

    所以 $\mathbf{x}(t)$ 也是方程组的解。当 $t$ 取遍所有实数时，得到无穷多个不同的解。$\blacksquare$

!!! example "例 1.6"
    判断下列方程组的解的情况：

    $$
    \begin{cases} x_1 + 2x_2 - x_3 = 1 \\ 2x_1 + 4x_2 - 2x_3 = 3 \end{cases}
    $$

    增广矩阵为

    $$
    \left(\begin{array}{ccc|c} 1 & 2 & -1 & 1 \\ 2 & 4 & -2 & 3 \end{array}\right) \xrightarrow{R_2 - 2R_1 \to R_2} \left(\begin{array}{ccc|c} 1 & 2 & -1 & 1 \\ 0 & 0 & 0 & 1 \end{array}\right)
    $$

    第二行为 $(0\; 0\; 0 \mid 1)$，对应方程 $0 = 1$，矛盾。方程组**无解**。

!!! example "例 1.7"
    解方程组

    $$
    \begin{cases} x_1 - 2x_2 + x_3 + x_4 = 3 \\ 2x_1 - 4x_2 + 3x_4 = 2 \\ -x_1 + 2x_2 - x_3 + 2x_4 = -1 \end{cases}
    $$

    增广矩阵化为简化行阶梯形：

    $$
    \left(\begin{array}{cccc|c} 1 & -2 & 1 & 1 & 3 \\ 2 & -4 & 0 & 3 & 2 \\ -1 & 2 & -1 & 2 & -1 \end{array}\right) \xrightarrow{\text{行变换}} \left(\begin{array}{cccc|c} 1 & -2 & 0 & 0 & 5 \\ 0 & 0 & 1 & 0 & -2 \\ 0 & 0 & 0 & 1 & -2 \end{array}\right)
    $$

    主元列为第 1、3、4 列，对应基本变量 $x_1, x_3, x_4$。第 2 列为非主元列，$x_2$ 是自由变量。设 $x_2 = t$，则

    $$
    x_1 = 5 + 2t, \quad x_2 = t, \quad x_3 = -2, \quad x_4 = -2.
    $$

---

## 1.6 齐次线性方程组

!!! definition "定义 1.10 (齐次线性方程组)"
    若线性方程组中所有常数项 $b_i = 0$，即方程组为

    $$
    A\mathbf{x} = \mathbf{0},
    $$

    则称其为**齐次线性方程组**（homogeneous system of linear equations）。否则称为**非齐次线性方程组**（nonhomogeneous system）。

齐次方程组总有一个解 $\mathbf{x} = \mathbf{0}$，称为**平凡解**（trivial solution）。非零解称为**非平凡解**（nontrivial solution）。

!!! theorem "定理 1.5 (非平凡解的存在条件)"
    齐次线性方程组 $A\mathbf{x} = \mathbf{0}$ 有非平凡解的充要条件是：方程组中存在自由变量。

    特别地，若方程的个数 $m$ 小于未知量的个数 $n$（即 $m < n$），则齐次方程组必有非平凡解。

??? proof "证明"
    齐次方程组总是相容的（$\mathbf{x} = \mathbf{0}$ 是解）。因此由定理 1.4，有非平凡解当且仅当存在自由变量。

    对于 $m < n$ 的情形：增广矩阵为 $m \times (n+1)$ 矩阵，系数矩阵为 $m \times n$ 矩阵。主元最多有 $m$ 个（每行至多一个主元），因此至少有 $n - m > 0$ 个自由变量，从而有非平凡解。$\blacksquare$

!!! theorem "定理 1.6 (齐次方程组解空间的性质)"
    齐次线性方程组 $A\mathbf{x} = \mathbf{0}$ 的解集 $S$ 满足：

    1. $\mathbf{0} \in S$。
    2. 若 $\mathbf{x}_1, \mathbf{x}_2 \in S$，则 $\mathbf{x}_1 + \mathbf{x}_2 \in S$（对加法封闭）。
    3. 若 $\mathbf{x}_1 \in S$，$c \in \mathbb{R}$，则 $c\mathbf{x}_1 \in S$（对标量乘法封闭）。

    因此 $S$ 构成 $\mathbb{R}^n$ 的一个**子空间**（subspace），称为 $A$ 的**零空间**（null space）或**核**。

??? proof "证明"
    1. $A\mathbf{0} = \mathbf{0}$，所以 $\mathbf{0} \in S$。

    2. 若 $A\mathbf{x}_1 = \mathbf{0}$ 且 $A\mathbf{x}_2 = \mathbf{0}$，则 $A(\mathbf{x}_1 + \mathbf{x}_2) = A\mathbf{x}_1 + A\mathbf{x}_2 = \mathbf{0} + \mathbf{0} = \mathbf{0}$。

    3. 若 $A\mathbf{x}_1 = \mathbf{0}$，则 $A(c\mathbf{x}_1) = cA\mathbf{x}_1 = c\mathbf{0} = \mathbf{0}$。$\blacksquare$

!!! example "例 1.8"
    求齐次方程组的通解：

    $$
    \begin{cases} x_1 + 2x_2 - x_3 + x_4 = 0 \\ 2x_1 + 4x_2 + x_3 - 2x_4 = 0 \\ 3x_1 + 6x_2 + 2x_4 = 0 \end{cases}
    $$

    增广矩阵化为简化行阶梯形：

    $$
    \left(\begin{array}{cccc|c} 1 & 2 & -1 & 1 & 0 \\ 2 & 4 & 1 & -2 & 0 \\ 3 & 6 & 0 & 2 & 0 \end{array}\right)
    \xrightarrow{\text{行变换}}
    \left(\begin{array}{cccc|c} 1 & 2 & 0 & -\frac{1}{3} & 0 \\ 0 & 0 & 1 & -\frac{4}{3} & 0 \\ 0 & 0 & 0 & 0 & 0 \end{array}\right)
    $$

    主元列为第 1、3 列。自由变量：$x_2 = s$，$x_4 = t$。通解：

    $$
    \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix}
    = s \begin{pmatrix} -2 \\ 1 \\ 0 \\ 0 \end{pmatrix}
    + t \begin{pmatrix} \frac{1}{3} \\ 0 \\ \frac{4}{3} \\ 1 \end{pmatrix}, \quad s, t \in \mathbb{R}.
    $$

---

## 1.7 解的结构

!!! theorem "定理 1.7 (非齐次方程组解的结构)"
    设 $A\mathbf{x} = \mathbf{b}$ 为一个相容的非齐次线性方程组，$\mathbf{x}_0$ 为其一个**特解**（particular solution），则方程组的通解可以表示为

    $$
    \mathbf{x} = \mathbf{x}_0 + \mathbf{x}_h,
    $$

    其中 $\mathbf{x}_h$ 是对应齐次方程组 $A\mathbf{x} = \mathbf{0}$ 的通解。

    换言之，**非齐次方程组的通解 = 一个特解 + 齐次方程组的通解**。

??? proof "证明"
    **（通解属于上述形式）** 设 $\mathbf{x}_1$ 是 $A\mathbf{x} = \mathbf{b}$ 的任意解，则

    $$
    A(\mathbf{x}_1 - \mathbf{x}_0) = A\mathbf{x}_1 - A\mathbf{x}_0 = \mathbf{b} - \mathbf{b} = \mathbf{0},
    $$

    故 $\mathbf{x}_h = \mathbf{x}_1 - \mathbf{x}_0$ 是齐次方程组的解，因而 $\mathbf{x}_1 = \mathbf{x}_0 + \mathbf{x}_h$。

    **（上述形式都是解）** 反过来，若 $A\mathbf{x}_h = \mathbf{0}$，则

    $$
    A(\mathbf{x}_0 + \mathbf{x}_h) = A\mathbf{x}_0 + A\mathbf{x}_h = \mathbf{b} + \mathbf{0} = \mathbf{b},
    $$

    故 $\mathbf{x}_0 + \mathbf{x}_h$ 确为非齐次方程组的解。$\blacksquare$

!!! note "注"
    此定理揭示了线性方程组解集的几何结构。齐次方程组 $A\mathbf{x} = \mathbf{0}$ 的解集是经过原点的子空间（直线、平面等），而非齐次方程组 $A\mathbf{x} = \mathbf{b}$ 的解集是该子空间的一个**平移**（仿射子空间），即将子空间整体平移到特解 $\mathbf{x}_0$ 的位置。

!!! example "例 1.9"
    解方程组

    $$
    \begin{cases} x_1 + 2x_2 + x_3 = 4 \\ 2x_1 + 4x_2 + 3x_3 = 9 \end{cases}
    $$

    增广矩阵化简：

    $$
    \left(\begin{array}{ccc|c} 1 & 2 & 1 & 4 \\ 2 & 4 & 3 & 9 \end{array}\right)
    \xrightarrow{R_2 - 2R_1}
    \left(\begin{array}{ccc|c} 1 & 2 & 1 & 4 \\ 0 & 0 & 1 & 1 \end{array}\right)
    \xrightarrow{R_1 - R_2}
    \left(\begin{array}{ccc|c} 1 & 2 & 0 & 3 \\ 0 & 0 & 1 & 1 \end{array}\right)
    $$

    **特解**（令自由变量 $x_2 = 0$）：$\mathbf{x}_0 = (3, 0, 1)^T$。

    **齐次通解**（$x_2 = t$）：$\mathbf{x}_h = t(-2, 1, 0)^T$。

    **通解**：

    $$
    \mathbf{x} = \begin{pmatrix} 3 \\ 0 \\ 1 \end{pmatrix} + t \begin{pmatrix} -2 \\ 1 \\ 0 \end{pmatrix}, \quad t \in \mathbb{R}.
    $$

!!! example "例 1.10"
    设 $4 \times 5$ 矩阵 $A$ 的秩为 $3$，$A\mathbf{x} = \mathbf{b}$ 有解。试分析解的结构。

    由于 $A$ 有 $5$ 个未知量，秩为 $3$（即有 $3$ 个主元列），所以有 $5 - 3 = 2$ 个自由变量。齐次方程组 $A\mathbf{x} = \mathbf{0}$ 的解空间是 $\mathbb{R}^5$ 的 $2$ 维子空间，设其基础解系为 $\{\mathbf{x}_1, \mathbf{x}_2\}$。

    非齐次方程组的通解为

    $$
    \mathbf{x} = \mathbf{x}_0 + c_1 \mathbf{x}_1 + c_2 \mathbf{x}_2, \quad c_1, c_2 \in \mathbb{R},
    $$

    其中 $\mathbf{x}_0$ 是某个特解。几何上，这是 $\mathbb{R}^5$ 中一个 $2$ 维仿射子空间（平面的推广）。
