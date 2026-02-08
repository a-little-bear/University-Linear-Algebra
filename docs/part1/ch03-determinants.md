# 第 3 章 行列式

<div class="context-flow" markdown>

**前置**：第 2 章方阵 · 可逆性 · 初等矩阵 · **本章脉络**：排列与逆序数 → 行列式定义 → 性质（行变换效应） → 余子式展开 → 伴随矩阵求逆 → 行列式乘法公式 → Cramer 法则 → 几何意义
**一句话本质**：行列式是方阵到标量的唯一交替多重线性函数——$\det(A) \neq 0 \Leftrightarrow A$ 可逆，$|\det(A)|$ = 体积缩放因子

</div>

行列式（determinant）是方阵特有的一个标量值，蕴含着矩阵的大量信息——可逆性、线性变换对体积的缩放效应、特征值之积等。本章从排列与逆序数出发建立行列式的严格定义，系统推导行列式的基本性质，介绍按行（列）展开定理与各种计算方法，最后证明 Cramer 法则并阐述行列式的几何意义。

---

## 3.1 排列与逆序数

<div class="context-flow" markdown>

**行列式的组合基础**：$n!$ 个排列的符号（奇/偶）决定了行列式定义中每项的正负号

</div>

!!! definition "定义 3.1 (排列)"
    由 $1, 2, \ldots, n$ 组成的一个有序排列 $\sigma = (j_1, j_2, \ldots, j_n)$ 称为 $\{1, 2, \ldots, n\}$ 的一个**排列**（permutation）。$n$ 个元素的全部排列共有 $n!$ 个。

!!! definition "定义 3.2 (逆序与逆序数)"
    在排列 $\sigma = (j_1, j_2, \ldots, j_n)$ 中，若 $s < t$ 但 $j_s > j_t$，则称 $(j_s, j_t)$ 为一个**逆序**（inversion）。排列 $\sigma$ 中所有逆序的总数称为 $\sigma$ 的**逆序数**（inversion number），记为 $\tau(\sigma)$ 或 $\operatorname{inv}(\sigma)$。

    - 若 $\tau(\sigma)$ 为偶数，则称 $\sigma$ 为**偶排列**（even permutation）。
    - 若 $\tau(\sigma)$ 为奇数，则称 $\sigma$ 为**奇排列**（odd permutation）。

    排列 $\sigma$ 的**符号**（sign / signature）定义为 $\operatorname{sgn}(\sigma) = (-1)^{\tau(\sigma)}$。

!!! example "例 3.1"
    排列 $(3, 1, 2)$ 的逆序为 $(3,1)$ 和 $(3,2)$，逆序数为 $2$，是偶排列，符号为 $+1$。

    排列 $(3, 2, 1)$ 的逆序为 $(3,2)$、$(3,1)$、$(2,1)$，逆序数为 $3$，是奇排列，符号为 $-1$。

!!! proposition "命题 3.1"
    在排列中交换（对换）两个元素的位置，逆序数的奇偶性改变。因此排列的奇偶性在一次对换后翻转。

??? proof "证明"
    **相邻对换**：交换相邻的 $j_k$ 和 $j_{k+1}$，仅影响这一对的逆序关系，逆序数变化 $\pm 1$，奇偶性改变。

    **一般对换**：交换位置 $s$ 和 $t$（$s < t$）可以分解为 $2(t-s)-1$ 次相邻对换（先把 $j_t$ 移到位置 $s$ 需要 $t-s$ 次，再把原来的 $j_s$ 从位置 $s+1$ 移到位置 $t$ 需要 $t-s-1$ 次），共 $2(t-s)-1$ 为奇数次相邻对换，因此奇偶性改变。$\blacksquare$

---

## 3.2 行列式的定义

<div class="context-flow" markdown>

**从排列到定义**：$\det(A) = \sum_\sigma \operatorname{sgn}(\sigma) \prod_i a_{i\sigma(i)}$ → 每项从每行每列各取恰好一个元素 → $n!$ 项求和

</div>

!!! definition "定义 3.3 ($n$ 阶行列式)"
    设 $A = (a_{ij})$ 为 $n$ 阶方阵，$A$ 的**行列式**（determinant）定义为

    $$
    \det(A) = \sum_{\sigma \in S_n} \operatorname{sgn}(\sigma) \cdot a_{1\sigma(1)} a_{2\sigma(2)} \cdots a_{n\sigma(n)},
    $$

    其中求和遍历 $\{1, 2, \ldots, n\}$ 的全部 $n!$ 个排列 $\sigma$，$\operatorname{sgn}(\sigma) = (-1)^{\tau(\sigma)}$。

    行列式也记为 $|A|$ 或

    $$
    \begin{vmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{n1} & a_{n2} & \cdots & a_{nn} \end{vmatrix}.
    $$

!!! example "例 3.2"
    $2$ 阶行列式：

    $$
    \begin{vmatrix} a & b \\ c & d \end{vmatrix} = ad - bc.
    $$

    $3$ 阶行列式（Sarrus 法则）：

    $$
    \begin{vmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{vmatrix} = a_{11}a_{22}a_{33} + a_{12}a_{23}a_{31} + a_{13}a_{21}a_{32} - a_{13}a_{22}a_{31} - a_{12}a_{21}a_{33} - a_{11}a_{23}a_{32}.
    $$

!!! example "例 3.3"
    计算

    $$
    \begin{vmatrix} 2 & 1 & 3 \\ 0 & -1 & 2 \\ 1 & 0 & 1 \end{vmatrix} = 2 \cdot(-1)\cdot 1 + 1\cdot 2\cdot 1 + 3\cdot 0\cdot 0 - 3\cdot(-1)\cdot 1 - 1\cdot 0\cdot 1 - 2\cdot 2\cdot 0 = -2 + 2 + 0 + 3 - 0 - 0 = 3.
    $$

---

## 3.3 行列式的性质

<div class="context-flow" markdown>

**计算基础**：行交换变号 → 行倍乘乘 $c$ → 行倍加不变 → 三者对应第 1 章三种初等行变换 → $\det(AB)=\det(A)\det(B)$ 是乘法公式的精髓

</div>

以下性质构成行列式计算的理论基础。

!!! theorem "定理 3.1 (转置不变性)"
    $\det(A^T) = \det(A)$。

??? proof "证明"
    $$
    \det(A^T) = \sum_{\sigma} \operatorname{sgn}(\sigma) a_{\sigma(1)1}a_{\sigma(2)2}\cdots a_{\sigma(n)n}.
    $$

    令 $\tau = \sigma^{-1}$，则 $\sigma(i) = j$ 等价于 $\tau(j) = i$，且 $\operatorname{sgn}(\tau) = \operatorname{sgn}(\sigma)$。

    $$
    \det(A^T) = \sum_{\tau} \operatorname{sgn}(\tau) a_{1\tau(1)}a_{2\tau(2)}\cdots a_{n\tau(n)} = \det(A). \quad \blacksquare
    $$

!!! note "注"
    由转置不变性，行列式关于行成立的性质对列也成立，反之亦然。

!!! theorem "定理 3.2 (行交换变号)"
    交换矩阵的两行（或两列），行列式变号。

??? proof "证明"
    交换第 $i$ 行和第 $j$ 行等价于在行列式定义中将排列 $\sigma$ 复合一个对换，由命题 3.1，排列符号取反，故行列式变号。$\blacksquare$

!!! theorem "定理 3.3 (行的齐次性)"
    若矩阵的某一行（或某一列）的所有元素都乘以常数 $c$，则行列式乘以 $c$。即

    $$
    \det(\ldots, c\mathbf{r}_i, \ldots) = c \cdot \det(\ldots, \mathbf{r}_i, \ldots),
    $$

    其中 $\mathbf{r}_i$ 表示第 $i$ 行。

??? proof "证明"
    在定义式中，每项恰包含第 $i$ 行的一个元素 $a_{i\sigma(i)}$。该元素被乘以 $c$ 后，每项都被乘以 $c$，故行列式乘以 $c$。$\blacksquare$

!!! corollary "推论 3.1"
    $\det(cA) = c^n \det(A)$，其中 $A$ 为 $n$ 阶方阵。

!!! theorem "定理 3.4 (行的加法性)"
    若矩阵的第 $i$ 行可以写成两个行向量之和 $\mathbf{r}_i = \mathbf{r}_i' + \mathbf{r}_i''$，则

    $$
    \det(\ldots, \mathbf{r}_i' + \mathbf{r}_i'', \ldots) = \det(\ldots, \mathbf{r}_i', \ldots) + \det(\ldots, \mathbf{r}_i'', \ldots).
    $$

??? proof "证明"
    在定义式中，第 $i$ 行元素为 $a_{i\sigma(i)}' + a_{i\sigma(i)}''$，展开后得到两部分之和，恰好对应右边两个行列式。$\blacksquare$

!!! theorem "定理 3.5 (两行相同则行列式为零)"
    若矩阵有两行（或两列）相同，则行列式为 $0$。

??? proof "证明"
    设第 $i$ 行和第 $j$ 行相同。交换这两行，矩阵不变但行列式变号（定理 3.2），故 $\det(A) = -\det(A)$，得 $2\det(A) = 0$，即 $\det(A) = 0$。$\blacksquare$

!!! theorem "定理 3.6 (行倍加不改变行列式)"
    将矩阵的第 $j$ 行的 $c$ 倍加到第 $i$ 行上，行列式不变。

??? proof "证明"
    $$
    \det(\ldots, \mathbf{r}_i + c\mathbf{r}_j, \ldots, \mathbf{r}_j, \ldots) = \det(\ldots, \mathbf{r}_i, \ldots, \mathbf{r}_j, \ldots) + c \cdot \det(\ldots, \mathbf{r}_j, \ldots, \mathbf{r}_j, \ldots).
    $$

    第二项中第 $i$ 行和第 $j$ 行都是 $\mathbf{r}_j$，由定理 3.5 为零。$\blacksquare$

!!! theorem "定理 3.7 (三角矩阵的行列式)"
    上三角矩阵（或下三角矩阵）的行列式等于其主对角线元素的乘积：

    $$
    \det(A) = a_{11}a_{22}\cdots a_{nn}.
    $$

??? proof "证明"
    以上三角矩阵为例。在行列式定义式 $\sum_\sigma \operatorname{sgn}(\sigma) a_{1\sigma(1)} \cdots a_{n\sigma(n)}$ 中，由于 $a_{ij} = 0$（$i > j$），要使乘积 $a_{1\sigma(1)}\cdots a_{n\sigma(n)} \neq 0$，需要对每个 $i$ 都有 $\sigma(i) \ge i$。又因为 $\sigma$ 是排列（双射），唯一可能是 $\sigma(i) = i$（恒等排列），其符号为 $+1$。因此 $\det(A) = a_{11}a_{22}\cdots a_{nn}$。$\blacksquare$

<div class="context-flow" markdown>

**关键洞察**：$\det(AB)=\det(A)\det(B)$ 使行列式成为从矩阵乘法群到标量乘法群的**同态** → 第 6 章推论：$\det A = \prod \lambda_i$

</div>

!!! theorem "定理 3.8 (行列式的乘法公式)"
    设 $A, B$ 为 $n$ 阶方阵，则

    $$
    \det(AB) = \det(A)\det(B).
    $$

??? proof "证明"
    **情形一**：$A$ 不可逆。则 $\operatorname{rank}(A) < n$，$AB$ 的列都是 $A$ 的列向量的线性组合，$\operatorname{rank}(AB) \le \operatorname{rank}(A) < n$，所以 $AB$ 不可逆，$\det(AB) = 0$。又 $\det(A) = 0$，故 $\det(AB) = 0 = \det(A)\det(B)$。

    **情形二**：$A$ 可逆。$A$ 可表示为初等矩阵的乘积 $A = E_1 E_2 \cdots E_k$。对初等矩阵，可以直接验证 $\det(EA) = \det(E)\det(A)$（分三种初等矩阵讨论，利用定理 3.2、3.3、3.6）。反复应用得

    $$
    \det(AB) = \det(E_1 E_2 \cdots E_k B) = \det(E_1)\det(E_2)\cdots\det(E_k)\det(B) = \det(A)\det(B). \quad \blacksquare
    $$

---

## 3.4 余子式与代数余子式

<div class="context-flow" markdown>

**降阶工具**：删去第 $i$ 行第 $j$ 列得 $(n-1)$ 阶子式 → 带符号 $(-1)^{i+j}$ 即为代数余子式 → 为 Laplace 展开和伴随矩阵做准备

</div>

!!! definition "定义 3.4 (余子式与代数余子式)"
    设 $A$ 为 $n$ 阶方阵。去掉 $A$ 的第 $i$ 行和第 $j$ 列后所得的 $(n-1)$ 阶子方阵的行列式称为元素 $a_{ij}$ 的**余子式**（minor），记为 $M_{ij}$。

    $a_{ij}$ 的**代数余子式**（cofactor）定义为

    $$
    A_{ij} = (-1)^{i+j} M_{ij}.
    $$

!!! example "例 3.4"
    设 $A = \begin{pmatrix} 2 & 1 & 3 \\ 0 & -1 & 2 \\ 1 & 0 & 1 \end{pmatrix}$。

    元素 $a_{11} = 2$ 的余子式为 $M_{11} = \begin{vmatrix} -1 & 2 \\ 0 & 1 \end{vmatrix} = -1$，代数余子式为 $A_{11} = (-1)^{1+1}(-1) = -1$。

    元素 $a_{12} = 1$ 的余子式为 $M_{12} = \begin{vmatrix} 0 & 2 \\ 1 & 1 \end{vmatrix} = -2$，代数余子式为 $A_{12} = (-1)^{1+2}(-2) = 2$。

---

## 3.5 行列式按行（列）展开

<div class="context-flow" markdown>

**递归计算**：$\det(A) = \sum_j a_{ij}A_{ij}$（按第 $i$ 行展开） → 异行展开为零 → 综合得 $A \cdot \operatorname{adj}(A) = \det(A) \cdot I$ → 由此推出 $A^{-1} = \frac{1}{\det A}\operatorname{adj}(A)$

</div>

!!! theorem "定理 3.9 (Laplace 展开定理)"
    $n$ 阶行列式 $\det(A)$ 可以按第 $i$ 行展开为

    $$
    \det(A) = \sum_{j=1}^n a_{ij} A_{ij} = a_{i1}A_{i1} + a_{i2}A_{i2} + \cdots + a_{in}A_{in},
    $$

    也可以按第 $j$ 列展开为

    $$
    \det(A) = \sum_{i=1}^n a_{ij} A_{ij} = a_{1j}A_{1j} + a_{2j}A_{2j} + \cdots + a_{nj}A_{nj}.
    $$

??? proof "证明"
    按第 $i$ 行展开。在行列式定义式中，将所有项按 $a_{i\sigma(i)}$ 的值分组。对于 $\sigma(i) = j$，将第 $i$ 行移到第一行（需 $i-1$ 次行交换），第 $j$ 列移到第一列（需 $j-1$ 次列交换），行列式变为 $(-1)^{(i-1)+(j-1)} = (-1)^{i+j}$ 倍。移除第一行第一列后，剩余部分恰好是 $M_{ij}$ 对应的 $(n-1)$ 阶行列式。因此

    $$
    \det(A) = \sum_{j=1}^n a_{ij} (-1)^{i+j} M_{ij} = \sum_{j=1}^n a_{ij} A_{ij}. \quad \blacksquare
    $$

!!! theorem "定理 3.10 (异行/列展开为零)"
    一行的元素与另一行对应的代数余子式的乘积之和为零：

    $$
    \sum_{k=1}^n a_{ik} A_{jk} = 0 \quad (i \neq j).
    $$

??? proof "证明"
    构造矩阵 $B$：将 $A$ 的第 $j$ 行替换为第 $i$ 行（其他行不变）。则 $B$ 有两行相同（第 $i$ 行和第 $j$ 行），故 $\det(B) = 0$。将 $\det(B)$ 按第 $j$ 行展开，注意 $B$ 的第 $j$ 行元素为 $a_{ik}$（$k = 1, \ldots, n$），而 $B$ 的第 $j$ 行代数余子式与 $A$ 的第 $j$ 行代数余子式相同（因为删去第 $j$ 行后，$B$ 与 $A$ 在其他行上一致），得

    $$
    0 = \det(B) = \sum_{k=1}^n a_{ik} A_{jk}. \quad \blacksquare
    $$

综合定理 3.9 和 3.10，可以写成

$$
\sum_{k=1}^n a_{ik}A_{jk} = \delta_{ij} \det(A),
$$

其中 $\delta_{ij}$ 为 Kronecker 符号。这等价于矩阵等式 $A \cdot \operatorname{adj}(A) = \det(A) \cdot I$。

!!! definition "定义 3.5 (伴随矩阵)"
    $n$ 阶方阵 $A$ 的**伴随矩阵**（adjugate matrix / classical adjoint）定义为代数余子式矩阵的转置：

    $$
    \operatorname{adj}(A) = (A_{ji})_{n \times n},
    $$

    即 $\operatorname{adj}(A)$ 的第 $(i, j)$ 元素为 $A_{ji}$。

!!! theorem "定理 3.11 (伴随矩阵与逆矩阵)"
    设 $A$ 为 $n$ 阶方阵，则

    $$
    A \cdot \operatorname{adj}(A) = \operatorname{adj}(A) \cdot A = \det(A) \cdot I.
    $$

    因此，若 $\det(A) \neq 0$，则 $A$ 可逆且

    $$
    A^{-1} = \frac{1}{\det(A)} \operatorname{adj}(A).
    $$

!!! example "例 3.5"
    用伴随矩阵法求逆矩阵。设 $A = \begin{pmatrix} 1 & 2 \\ 3 & 5 \end{pmatrix}$。

    $\det(A) = 5 - 6 = -1$。代数余子式：$A_{11} = 5$，$A_{12} = -3$，$A_{21} = -2$，$A_{22} = 1$。

    $$
    \operatorname{adj}(A) = \begin{pmatrix} 5 & -2 \\ -3 & 1 \end{pmatrix}, \quad A^{-1} = \frac{1}{-1}\begin{pmatrix} 5 & -2 \\ -3 & 1 \end{pmatrix} = \begin{pmatrix} -5 & 2 \\ 3 & -1 \end{pmatrix}.
    $$

---

## 3.6 行列式的计算方法

<div class="context-flow" markdown>

**实用方法**：三角化（利用性质 3.2–3.6） · 递推法（利用展开建立递推关系） · Vandermonde 行列式（经典公式） · 分块行列式

</div>

### 三角化法

利用行变换将矩阵化为上三角形（注意行交换变号、行倍乘乘以倍数），再利用定理 3.7 直接求行列式。

!!! example "例 3.6"
    计算

    $$
    D = \begin{vmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 0 \end{vmatrix}
    $$

    $$
    \xrightarrow{R_2 - 4R_1, \; R_3 - 7R_1}
    \begin{vmatrix} 1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & -6 & -21 \end{vmatrix}
    \xrightarrow{R_3 - 2R_2}
    \begin{vmatrix} 1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & 0 & -9 \end{vmatrix} = 1 \cdot (-3) \cdot (-9) = 27.
    $$

### 递推法

对于具有特殊规律的行列式，可以建立递推关系来求解。

!!! example "例 3.7"
    设 $D_n$ 为 $n$ 阶行列式

    $$
    D_n = \begin{vmatrix} 2 & 1 & 0 & \cdots & 0 & 0 \\ 1 & 2 & 1 & \cdots & 0 & 0 \\ 0 & 1 & 2 & \cdots & 0 & 0 \\ \vdots & & & \ddots & & \vdots \\ 0 & 0 & 0 & \cdots & 2 & 1 \\ 0 & 0 & 0 & \cdots & 1 & 2 \end{vmatrix}
    $$

    按第一行展开：$D_n = 2D_{n-1} - D_{n-2}$。初始值 $D_1 = 2$，$D_2 = 3$。

    特征方程 $x^2 - 2x + 1 = 0$ 有重根 $x = 1$，通解 $D_n = (c_1 + c_2 n) \cdot 1^n = c_1 + c_2 n$。

    由 $D_1 = 2, D_2 = 3$ 解得 $c_1 = 1, c_2 = 1$，故 $D_n = n + 1$。

### Vandermonde 行列式

!!! theorem "定理 3.12 (Vandermonde 行列式)"
    **Vandermonde 行列式**（Vandermonde determinant）为

    $$
    V_n = \begin{vmatrix} 1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\ 1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 1 & x_n & x_n^2 & \cdots & x_n^{n-1} \end{vmatrix} = \prod_{1 \le i < j \le n} (x_j - x_i).
    $$

??? proof "证明"
    对 $n$ 进行归纳。$n = 2$ 时，$V_2 = x_2 - x_1$，成立。

    设 $n-1$ 阶时公式成立。对 $n$ 阶 Vandermonde 行列式，从最后一行开始，每行减去上一行的 $x_1$ 倍（即 $R_n - x_1 R_{n-1}$，$R_{n-1} - x_1 R_{n-2}$，...，$R_2 - x_1 R_1$）：

    $$
    V_n = \begin{vmatrix} 1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\ 0 & x_2 - x_1 & x_2(x_2 - x_1) & \cdots & x_2^{n-2}(x_2 - x_1) \\ 0 & x_3 - x_1 & x_3(x_3 - x_1) & \cdots & x_3^{n-2}(x_3 - x_1) \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & x_n - x_1 & x_n(x_n - x_1) & \cdots & x_n^{n-2}(x_n - x_1) \end{vmatrix}
    $$

    按第一列展开，并从第 $i$ 行（$i = 2, \ldots, n$）提出公因子 $(x_i - x_1)$：

    $$
    V_n = \prod_{i=2}^n (x_i - x_1) \cdot V_{n-1}(x_2, x_3, \ldots, x_n)
    $$

    由归纳假设，$V_{n-1}(x_2, \ldots, x_n) = \prod_{2 \le i < j \le n}(x_j - x_i)$，因此

    $$
    V_n = \prod_{i=2}^n(x_i - x_1) \cdot \prod_{2 \le i < j \le n}(x_j - x_i) = \prod_{1 \le i < j \le n}(x_j - x_i). \quad \blacksquare
    $$

### 分块行列式

!!! theorem "定理 3.13 (分块三角矩阵的行列式)"
    设

    $$
    M = \begin{pmatrix} A & B \\ O & D \end{pmatrix} \quad \text{或} \quad M = \begin{pmatrix} A & O \\ C & D \end{pmatrix},
    $$

    其中 $A$、$D$ 为方阵，则 $\det(M) = \det(A)\det(D)$。

!!! example "例 3.8"
    计算分块矩阵的行列式：

    $$
    M = \begin{pmatrix} 1 & 2 & 0 & 0 \\ 3 & 4 & 0 & 0 \\ 5 & 6 & 1 & 1 \\ 7 & 8 & 0 & 2 \end{pmatrix}
    $$

    这是分块下三角矩阵 $\begin{pmatrix} A & O \\ C & D \end{pmatrix}$，其中 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，$D = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}$。

    $\det(M) = \det(A)\det(D) = (4 - 6)(2 - 0) = (-2)(2) = -4$。

---

## 3.7 Cramer 法则

<div class="context-flow" markdown>

**行列式解方程**：$x_j = \det(A_j)/\det(A)$ → 理论意义 > 计算意义（$n+1$ 个行列式太贵） → 实际求解仍用第 1 章高斯消元

</div>

!!! theorem "定理 3.14 (Cramer 法则)"
    设 $A$ 为 $n$ 阶可逆矩阵（即 $\det(A) \neq 0$），则线性方程组 $A\mathbf{x} = \mathbf{b}$ 有唯一解，且

    $$
    x_j = \frac{\det(A_j)}{\det(A)}, \quad j = 1, 2, \ldots, n,
    $$

    其中 $A_j$ 是将 $A$ 的第 $j$ 列替换为 $\mathbf{b}$ 所得的矩阵。

??? proof "证明"
    由 $A$ 可逆，方程有唯一解 $\mathbf{x} = A^{-1}\mathbf{b}$。利用伴随矩阵公式：

    $$
    x_j = (A^{-1}\mathbf{b})_j = \frac{1}{\det(A)}\sum_{i=1}^n A_{ji} b_i.
    $$

    而 $\sum_{i=1}^n A_{ji} b_i$ 恰好是将 $A$ 的第 $j$ 列替换为 $\mathbf{b}$ 后的矩阵 $A_j$ 按第 $j$ 列展开的结果，即 $\det(A_j)$。因此 $x_j = \frac{\det(A_j)}{\det(A)}$。$\blacksquare$

!!! note "注"
    Cramer 法则的理论意义大于计算意义。对于大型方程组，高斯消元法的效率远高于 Cramer 法则（后者需要计算 $n+1$ 个 $n$ 阶行列式）。

!!! example "例 3.9"
    用 Cramer 法则解方程组

    $$
    \begin{cases} 2x_1 + x_2 = 5 \\ 3x_1 + 2x_2 = 8 \end{cases}
    $$

    $$
    \det(A) = \begin{vmatrix} 2 & 1 \\ 3 & 2 \end{vmatrix} = 1, \quad \det(A_1) = \begin{vmatrix} 5 & 1 \\ 8 & 2 \end{vmatrix} = 2, \quad \det(A_2) = \begin{vmatrix} 2 & 5 \\ 3 & 8 \end{vmatrix} = 1.
    $$

    $$
    x_1 = \frac{2}{1} = 2, \quad x_2 = \frac{1}{1} = 1.
    $$

---

## 3.8 行列式的几何意义

<div class="context-flow" markdown>

**代数与几何的桥梁**：$|\det A|$ = 列向量张成的平行体的体积 → $\det A$ 的符号 = 定向（右手/左手系） → 线性变换面积缩放 $|\det A|$ 倍（→ 第 5 章）

</div>

行列式与几何中的面积和体积有着深刻的联系。

!!! theorem "定理 3.15 (行列式与面积/体积)"
    设 $\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n \in \mathbb{R}^n$ 为 $n$ 个列向量，构成矩阵 $A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)$。

    1. **二维**：以 $\mathbf{a}_1, \mathbf{a}_2 \in \mathbb{R}^2$ 为邻边的平行四边形的**有向面积**等于 $\det(A) = \begin{vmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{vmatrix}$。

    2. **三维**：以 $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3 \in \mathbb{R}^3$ 为邻边的平行六面体的**有向体积**等于 $\det(A)$。

    3. **一般情形**：以 $\mathbf{a}_1, \ldots, \mathbf{a}_n$ 为邻边的 $n$ 维平行体的 $n$ 维体积的绝对值为 $|\det(A)|$。

!!! note "注"
    $\det(A) > 0$ 表示向量组 $\{\mathbf{a}_1, \ldots, \mathbf{a}_n\}$ 保持标准基的定向（右手系），$\det(A) < 0$ 表示反转定向。$\det(A) = 0$ 表示向量组线性相关，"平行体"退化为低维对象，体积为零。

!!! example "例 3.10"
    向量 $\mathbf{a}_1 = (3, 0)^T$ 和 $\mathbf{a}_2 = (1, 2)^T$ 张成的平行四边形的面积为

    $$
    |\det(A)| = \left|\begin{vmatrix} 3 & 1 \\ 0 & 2 \end{vmatrix}\right| = |6 - 0| = 6.
    $$

    几何上，这个平行四边形底边长 $3$、高 $2$，面积确实为 $6$。

!!! example "例 3.11"
    线性变换 $T: \mathbb{R}^2 \to \mathbb{R}^2$，$T(\mathbf{x}) = A\mathbf{x}$，将平面上的区域面积缩放为原来的 $|\det(A)|$ 倍。例如，若 $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$，则 $\det(A) = 6$，意味着 $T$ 将任何区域的面积放大到原来的 $6$ 倍。
