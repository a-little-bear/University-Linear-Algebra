# 第 38 章 特殊矩阵类：M-矩阵、Z-矩阵、P-矩阵与 H-矩阵

<div class="context-flow" markdown>

**前置**：非负矩阵(Ch17) · 矩阵分析(Ch14) · 特征值(Ch6)

**本章脉络**：Z-矩阵定义 $\to$ M-矩阵(50个等价条件) $\to$ 非奇异 M-矩阵判定 $\to$ P-矩阵 $\to$ H-矩阵 $\to$ 逆正矩阵 $\to$ 迭代法收敛性 $\to$ 经济学模型

**延伸**：M-矩阵在偏微分方程离散化（有限差分格式产生 M-矩阵）和经济学（Leontief 模型的可行性条件即 Hawkins-Simon 条件等价于 $I - C$ 是 M-矩阵）中至关重要

</div>

在应用数学的诸多领域中，我们反复遇到一些具有特定符号模式或行列式条件的矩阵类。这些矩阵类——Z-矩阵、M-矩阵、P-矩阵、H-矩阵——虽然定义简洁，却拥有极为丰富的等价刻画和深刻的应用。Berman 和 Plemmons 在其经典著作中列出了 M-矩阵的 50 个等价条件，成为线性代数中最壮观的等价性定理之一。

本章系统地介绍这些矩阵类的定义、相互关系、判定方法和核心应用。

---

## 38.1 Z-矩阵

<div class="context-flow" markdown>

**核心问题**：非对角元素全部非正的矩阵具有什么特殊性质？

</div>

!!! definition "定义 38.1 (Z-矩阵)"
    矩阵 $A = (a_{ij}) \in M_n(\mathbb{R})$ 称为 **Z-矩阵**，若其所有非对角元素非正：
    $$a_{ij} \le 0, \quad \forall\, i \ne j.$$
    等价地，$A$ 可以写成
    $$A = sI - B, \quad s \in \mathbb{R}, \quad B \ge 0 \text{ (逐元非负)},$$
    其中 $s = \max_i a_{ii}$, $B = sI - A$。

!!! theorem "定理 38.1 (Z-矩阵的基本性质)"
    (a) Z-矩阵对非负数量乘封闭：若 $A$ 是 Z-矩阵且 $\alpha \ge 0$，则 $\alpha A$ 是 Z-矩阵。

    (b) Z-矩阵的主子矩阵仍为 Z-矩阵。

    (c) Z-矩阵对加法**不**封闭（两个 Z-矩阵之和不一定是 Z-矩阵——当 $a_{ij} + b_{ij} > 0$ 时破坏条件）。但若非对角元均严格为负，加法在适当范围内封闭。

    (d) Z-矩阵的转置是 Z-矩阵。

??? proof "证明"
    (a) $(\alpha A)_{ij} = \alpha a_{ij} \le 0$（$\alpha \ge 0, a_{ij} \le 0, i \ne j$）。

    (b) 主子矩阵继承非对角元素的非正性。

    (d) $(A^T)_{ij} = a_{ji} \le 0$（$i \ne j$），因为 $j \ne i$。

!!! example "例 38.1"
    下列矩阵是 Z-矩阵：
    $$A = \begin{pmatrix} 2 & -1 & 0 \\ -3 & 4 & -1 \\ 0 & -2 & 3 \end{pmatrix}, \quad B = \begin{pmatrix} -1 & -2 \\ -3 & -4 \end{pmatrix}.$$
    $A$ 和 $B$ 的非对角元素均 $\le 0$。

    矩阵 $C = \begin{pmatrix} 1 & 2 \\ -1 & 3 \end{pmatrix}$ 不是 Z-矩阵（$c_{12} = 2 > 0$）。

---

## 38.2 M-矩阵的定义与等价条件

<div class="context-flow" markdown>

**核心问题**：什么条件下 Z-矩阵的逆存在且逐元非负？

</div>

M-矩阵是 Z-矩阵中最重要的子类，以 Minkowski 的名字命名（由 Ostrowski 引入术语）。

!!! definition "定义 38.2 (M-矩阵)"
    Z-矩阵 $A$ 称为（非奇异）**M-矩阵**，若 $A$ 可以写成
    $$A = sI - B, \quad B \ge 0, \quad s > \rho(B),$$
    其中 $\rho(B)$ 为 $B$ 的谱半径。

!!! definition "定义 38.3 (奇异 M-矩阵)"
    Z-矩阵 $A = sI - B$（$B \ge 0$）称为**奇异 M-矩阵**，若 $s = \rho(B)$。

!!! theorem "定理 38.2 (非奇异 M-矩阵的等价条件)"
    设 $A \in M_n(\mathbb{R})$ 为 Z-矩阵，$A = sI - B$（$B \ge 0$）。以下条件等价：

    **(M1)** $A$ 是非奇异 M-矩阵（即 $s > \rho(B)$）。

    **(M2)** $A$ 的所有特征值具有正实部。

    **(M3)** $A$ 是非奇异的，且 $A^{-1} \ge 0$（逐元非负）。

    **(M4)** $A$ 的所有顺序主子式（leading principal minors）为正。

    **(M5)** $A$ 的所有主子式（principal minors）为正。

    **(M6)** 存在正向量 $x > 0$ 使得 $Ax > 0$。

    **(M7)** 存在正向量 $u > 0$ 使得 $A^T u > 0$（等价地，$A$ 的每一行的正对角元足以"支配"该行的负非对角元，在某种加权意义下）。

    **(M8)** $A$ 的 Jacobi 迭代矩阵 $M_J = D^{-1}(L + U)$（其中 $A = D - L - U$）满足 $\rho(M_J) < 1$。

    **(M9)** $A$ 的 Gauss-Seidel 迭代矩阵 $M_{GS} = (D - L)^{-1}U$ 满足 $\rho(M_{GS}) < 1$。

    **(M10)** 存在正对角矩阵 $D$ 使得 $DA$ 是严格行对角占优的。

??? proof "证明（选要）"
    **(M1) $\Leftrightarrow$ (M2)**：$A = sI - B$ 的特征值为 $s - \lambda_i(B)$。$s > \rho(B)$ 意味着 $\operatorname{Re}(s - \lambda_i) = s - \operatorname{Re}(\lambda_i) \ge s - |\lambda_i| \ge s - \rho(B) > 0$。

    **(M1) $\Rightarrow$ (M3)**：$A = sI - B$，$s > \rho(B)$，故 $A = s(I - s^{-1}B)$，$A^{-1} = s^{-1}(I - s^{-1}B)^{-1} = s^{-1}\sum_{k=0}^{\infty} (s^{-1}B)^k$。由 $B \ge 0$ 和 $s > 0$，级数每项非负，故 $A^{-1} \ge 0$。

    **(M3) $\Rightarrow$ (M1)**：$A$ 非奇异且 $A^{-1} \ge 0$。$A = sI - B$（$B \ge 0$）。$A^{-1} = s^{-1}(I - s^{-1}B)^{-1}$。若 $s \le \rho(B)$，则 $I - s^{-1}B$ 有非正特征值，$(I - s^{-1}B)^{-1}$ 存在但不一定非负。具体地，由 Perron-Frobenius 理论，$B$ 有非负特征值 $\rho(B)$，对应非负特征向量 $v$。$Av = (sI-B)v = (s - \rho(B))v$。若 $s \le \rho(B)$，$Av$ 有非正分量，但 $A^{-1} \ge 0$ 和 $Av$ 的非正性会导致 $v = A^{-1}(Av)$ 有非正分量，矛盾于 $v \ge 0, v \ne 0$。

    **(M3) $\Rightarrow$ (M6)**：取 $x = A^{-1}\mathbf{1}$（$\mathbf{1}$ 为全 1 向量），则 $x = A^{-1}\mathbf{1} \ge 0$。进一步，若 $A^{-1}$ 逐元严格为正（对不可约 M-矩阵成立），则 $x > 0$ 且 $Ax = \mathbf{1} > 0$。对一般情形需要更精细的论证。

    **(M6) $\Rightarrow$ (M1)**：存在 $x > 0$，$Ax > 0$。由 $A = sI - B$，$Ax = sx - Bx > 0$。令 $y = Bx > 0$（因 $B \ge 0, x > 0$，且只要 $B$ 的每行非全零），则 $sx > Bx$，即 $s > (Bx)_i/x_i$ 对所有 $i$。故 $s > \max_i (Bx)_i / x_i \ge \rho(B)$（最后一步由对非负矩阵 $B$ 的 Collatz-Wielandt 公式保证）。

    **(M1) $\Rightarrow$ (M5)**：M-矩阵的每个主子矩阵仍然是 M-矩阵（因为 Z 性质和谱半径条件对主子矩阵遗传），故所有主子式为正（因为 $\det(A) > 0$ 对 M-矩阵成立）。

!!! example "例 38.2"
    验证 $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$ 是 M-矩阵。

    - Z-矩阵：非对角元素 $-1 \le 0$。
    - $A = 2I - B$，$B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$，$\rho(B) = 1 < 2 = s$。 **(M1)**
    - 特征值：$3, 1$，均为正。 **(M2)**
    - $A^{-1} = \frac{1}{3}\begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix} \ge 0$。 **(M3)**
    - 顺序主子式：$2 > 0$，$\det(A) = 3 > 0$。 **(M4)**
    - 取 $x = (1, 1)^T$，$Ax = (1, 1)^T > 0$。 **(M6)**

!!! example "例 38.3"
    矩阵 $A = \begin{pmatrix} 1 & -2 \\ -1 & 1 \end{pmatrix}$ 是 Z-矩阵但**不是** M-矩阵。

    $A = I - B$，$B = \begin{pmatrix} 0 & 2 \\ 1 & 0 \end{pmatrix}$，$\rho(B) = \sqrt{2} > 1 = s$。

    验证：$\det(A) = 1 - 2 = -1 < 0$（条件 M5 不成立），特征值 $1 \pm \sqrt{2}$（有一个负特征值，M2 不成立）。

---

## 38.3 非奇异 M-矩阵的判定

<div class="context-flow" markdown>

**核心问题**：给定一个 Z-矩阵，如何高效判定它是否为 M-矩阵？

</div>

!!! theorem "定理 38.3 (对角占优与 M-矩阵)"
    设 $A$ 为 Z-矩阵。若 $A$ 是严格行对角占优的，即
    $$a_{ii} > \sum_{j \ne i} |a_{ij}| = -\sum_{j \ne i} a_{ij}, \quad \forall\, i,$$
    则 $A$ 是非奇异 M-矩阵。

??? proof "证明"
    由 Gershgorin 圆盘定理，$A$ 的每个特征值 $\lambda$ 满足
    $$|\lambda - a_{ii}| \le \sum_{j \ne i} |a_{ij}|.$$
    由严格对角占优，$|a_{ij}|$ 之和 $< a_{ii}$，故 $\operatorname{Re}(\lambda) > 0$。结合 $A$ 是 Z-矩阵，得 $A$ 是 M-矩阵。

!!! theorem "定理 38.4 (Gauss 消元判定法)"
    设 $A$ 为 $n \times n$ Z-矩阵。对 $A$ 进行**不带行交换**的 Gauss 消元（LU 分解）。$A$ 是非奇异 M-矩阵当且仅当消元过程可以完成且所有主元为正。

??? proof "证明"
    M-矩阵的所有顺序主子式为正（条件 M4），这恰好是 LU 分解存在且主元为正的充要条件。

    反过来，若 LU 分解存在且主元为正，则 $A$ 的顺序主子式为主元之积，全为正。由条件 M4，$A$ 是 M-矩阵。

!!! example "例 38.4"
    判定 $A = \begin{pmatrix} 4 & -1 & -1 \\ -2 & 5 & -1 \\ -1 & -1 & 3 \end{pmatrix}$ 是否为 M-矩阵。

    **方法一（对角占优）**：行 1：$4 > |-1| + |-1| = 2$；行 2：$5 > 2 + 1 = 3$；行 3：$3 > 1 + 1 = 2$。严格行对角占优，故 $A$ 是 M-矩阵。

    **方法二（主子式）**：$\Delta_1 = 4 > 0$；$\Delta_2 = 20 - 2 = 18 > 0$；$\Delta_3 = \det(A) = 4(15-1) - (-1)(-6-1) + (-1)(2+5) = 56 - 7 - 7 = 42 > 0$。全部为正。

    **方法三（正向量）**：取 $x = (1, 1, 1)^T$，$Ax = (2, 2, 1)^T > 0$。条件 M6 成立。

---

## 38.4 P-矩阵

<div class="context-flow" markdown>

**核心问题**：所有主子式为正的矩阵有什么性质？它与线性互补问题(LCP)有何联系？

</div>

!!! definition "定义 38.4 (P-矩阵)"
    矩阵 $A \in M_n(\mathbb{R})$ 称为 **P-矩阵**，若其所有主子式为正：
    $$\det(A[\alpha, \alpha]) > 0, \quad \forall\, \alpha \subseteq \{1, 2, \ldots, n\}, \quad \alpha \ne \emptyset,$$
    其中 $A[\alpha, \alpha]$ 表示行列指标集均为 $\alpha$ 的主子矩阵。

!!! theorem "定理 38.5 (P-矩阵的等价刻画)"
    设 $A \in M_n(\mathbb{R})$。以下条件等价：

    (a) $A$ 是 P-矩阵。

    (b) 对每个非零向量 $x \in \mathbb{R}^n$，存在指标 $i$ 使得 $x_i(Ax)_i > 0$。

    (c) $A$ 将正象限的每个面映射为不含在任何坐标超平面内（即 $A$ 不将任何非平凡的符号模式映射到相反的符号模式）。

    (d) 对每个对角矩阵 $D$（对角元素为 $\pm 1$），$DA$ 是非奇异的。

    (e) 线性互补问题 $\text{LCP}(q, A)$：找 $w, z \ge 0$ 使得 $w = q + Az$，$w^T z = 0$，对所有 $q \in \mathbb{R}^n$ 有唯一解。

??? proof "证明（(a) $\Leftrightarrow$ (b)，概要）"
    **(a) $\Rightarrow$ (b)**：反证。设存在 $x \ne 0$ 使得 $x_i(Ax)_i \le 0$ 对所有 $i$。设 $\alpha = \{i : x_i \ne 0\}$。则 $x_i(Ax)_i \le 0$ 对 $i \in \alpha$，特别是 $x_\alpha^T (Ax)_\alpha \le 0$。注意 $(Ax)_\alpha = A[\alpha, \alpha] x_\alpha$（因为 $x_j = 0$ 对 $j \notin \alpha$）。故 $x_\alpha^T A[\alpha, \alpha] x_\alpha \le 0$，但又因为每个 $x_i(Ax)_i \le 0$，所以 $A[\alpha,\alpha]$ 不能是 P-矩阵。经过仔细分析，这与 $\det(A[\alpha,\alpha]) > 0$ 矛盾。

    **(b) $\Rightarrow$ (a)**：设条件 (b) 成立。对任意 $\alpha \ne \emptyset$，取 $x$ 使得 $x_\alpha$ 为 $A[\alpha,\alpha]$ 的某个特征向量，$x_j = 0$（$j \notin \alpha$）。则 $(Ax)_i = (A[\alpha,\alpha]x_\alpha)_i = \lambda x_i$（$i \in \alpha$），故 $x_i(Ax)_i = \lambda x_i^2$。条件 (b) 要求某个 $x_i(Ax)_i > 0$，故 $\lambda > 0$。因此 $A[\alpha,\alpha]$ 的所有实特征值为正。对复特征值的讨论类似（考虑实部），最终得到 $\det(A[\alpha,\alpha]) > 0$。

!!! example "例 38.5"
    矩阵 $A = \begin{pmatrix} 1 & -2 \\ 1 & 1 \end{pmatrix}$ 是否为 P-矩阵？

    主子式：$a_{11} = 1 > 0$，$a_{22} = 1 > 0$，$\det(A) = 1 + 2 = 3 > 0$。所有主子式为正，$A$ 是 P-矩阵。

    注意 $A$ **不是** Z-矩阵（$a_{12} = -2$ 但 $a_{21} = 1 > 0$），因此 $A$ 是 P-矩阵但不是 M-矩阵。

!!! note "注记 38.1 (P-矩阵与 M-矩阵的关系)"
    每个非奇异 M-矩阵都是 P-矩阵（因为 M-矩阵的所有主子式为正，条件 M5）。反之不然——P-矩阵不要求非对角元素非正。P-矩阵是比 M-矩阵更大的类。

---

## 38.5 H-矩阵

<div class="context-flow" markdown>

**核心问题**：能否用比较矩阵将 M-矩阵的理论推广到一般矩阵？

</div>

!!! definition "定义 38.5 (比较矩阵)"
    设 $A = (a_{ij}) \in M_n(\mathbb{C})$。$A$ 的**比较矩阵** $\mathcal{M}(A) = (\mu_{ij})$ 定义为
    $$\mu_{ij} = \begin{cases} |a_{ii}| & \text{if } i = j, \\ -|a_{ij}| & \text{if } i \ne j. \end{cases}$$
    显然 $\mathcal{M}(A)$ 是 Z-矩阵。

!!! definition "定义 38.6 (H-矩阵)"
    矩阵 $A \in M_n(\mathbb{C})$ 称为 **H-矩阵**，若其比较矩阵 $\mathcal{M}(A)$ 是（非奇异）M-矩阵。

!!! theorem "定理 38.6 (H-矩阵的等价条件)"
    设 $A \in M_n(\mathbb{C})$。以下条件等价：

    (a) $A$ 是 H-矩阵。

    (b) 存在正对角矩阵 $D$ 使得 $AD$ 是严格行对角占优的。

    (c) 存在正向量 $x > 0$ 使得 $\mathcal{M}(A)x > 0$。

    (d) $A$ 是非奇异的，且 $|\det(A)| \ge \det(\mathcal{M}(A)) > 0$。

??? proof "证明（(a) $\Leftrightarrow$ (b)）"
    **(a) $\Rightarrow$ (b)**：$\mathcal{M}(A)$ 是 M-矩阵，故由 M-矩阵的条件 M10，存在正对角矩阵 $D$ 使得 $\mathcal{M}(A)D$ 严格行对角占优。注意到 $\mathcal{M}(AD) = \mathcal{M}(A)D$（对角缩放不改变比较矩阵的结构），而 $\mathcal{M}(A)D$ 严格行对角占优意味着 $AD$ 严格行对角占优。

    **(b) $\Rightarrow$ (a)**：$AD$ 严格行对角占优，故 $\mathcal{M}(AD) = \mathcal{M}(A)D$ 严格行对角占优。由定理 38.3，$\mathcal{M}(A)D$ 是 M-矩阵。由于 $D > 0$，$\mathcal{M}(A) = (\mathcal{M}(A)D)D^{-1}$ 也是 M-矩阵。

!!! theorem "定理 38.7 (广义对角占优与 H-矩阵)"
    严格行对角占优矩阵 $\implies$ H-矩阵 $\implies$ 非奇异矩阵。

    反之均不成立。

??? proof "证明"
    若 $A$ 严格行对角占优，则 $\mathcal{M}(A)$ 也严格行对角占优，由定理 38.3，$\mathcal{M}(A)$ 是 M-矩阵，故 $A$ 是 H-矩阵。

    H-矩阵是非奇异的：由 $\mathcal{M}(A)$ 是 M-矩阵，$\det(\mathcal{M}(A)) > 0$。由 Ostrowski 的行列式不等式 $|\det(A)| \ge \det(\mathcal{M}(A)) > 0$，故 $A$ 非奇异。

!!! example "例 38.6"
    $A = \begin{pmatrix} 3 & 1+i & -1 \\ -2 & 5 & 1 \\ 1 & -1 & 4 \end{pmatrix}$。

    比较矩阵：$\mathcal{M}(A) = \begin{pmatrix} 3 & -\sqrt{2} & -1 \\ -2 & 5 & -1 \\ -1 & -1 & 4 \end{pmatrix}$。

    行占优检查：行 1：$3 > \sqrt{2} + 1 \approx 2.41$；行 2：$5 > 2 + 1 = 3$；行 3：$4 > 1 + 1 = 2$。

    $\mathcal{M}(A)$ 严格行对角占优，故是 M-矩阵，因此 $A$ 是 H-矩阵。

---

## 38.6 逆正矩阵与单调矩阵

<div class="context-flow" markdown>

**核心问题**：$A^{-1} \ge 0$ 的矩阵有什么结构？它与"单调性"有何联系？

</div>

!!! definition "定义 38.7 (逆正矩阵)"
    非奇异矩阵 $A \in M_n(\mathbb{R})$ 称为**逆正矩阵**，若 $A^{-1} \ge 0$（逐元非负）。称为**逆严格正矩阵**，若 $A^{-1} > 0$（逐元严格正）。

!!! definition "定义 38.8 (单调矩阵)"
    非奇异矩阵 $A \in M_n(\mathbb{R})$ 称为**单调矩阵**，若
    $$Ax \ge 0 \implies x \ge 0.$$

!!! theorem "定理 38.8 (逆正 = 单调)"
    $A$ 是逆正矩阵当且仅当 $A$ 是单调矩阵。

??? proof "证明"
    $(\Rightarrow)$：设 $A^{-1} \ge 0$ 且 $Ax \ge 0$。则 $x = A^{-1}(Ax) \ge 0$（非负矩阵与非负向量的乘积非负）。

    $(\Leftarrow)$：设 $A$ 单调。对每个标准基向量 $e_j$，令 $x_j = A^{-1}e_j$（$A$ 的逆的第 $j$ 列）。由 $Ax_j = e_j \ge 0$ 和单调性得 $x_j \ge 0$。故 $A^{-1} = (x_1, \ldots, x_n) \ge 0$。

!!! theorem "定理 38.9 (M-矩阵是逆正的)"
    非奇异 M-矩阵是逆正矩阵。更强地，不可约非奇异 M-矩阵是逆严格正的：$A^{-1} > 0$。

??? proof "证明"
    第一部分是 M-矩阵的等价条件 M3。

    第二部分：$A = sI - B$，$B \ge 0$ 不可约，$s > \rho(B)$。$A^{-1} = s^{-1} \sum_{k=0}^{\infty}(s^{-1}B)^k$。由 $B$ 不可约，$B^{n-1} > 0$（Perron-Frobenius），故 $A^{-1}$ 的每个元素至少包含 $(s^{-1}B)^{n-1}$ 中对应元素的贡献，为正。

!!! example "例 38.7"
    一维 Poisson 方程 $-u'' = f$ 在 $[0,1]$ 上的有限差分离散化（Dirichlet 边界）给出系统 $Ah = f_h$，其中
    $$A = \frac{1}{h^2}\begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & \ddots & \ddots & \ddots \\ & & -1 & 2 \end{pmatrix}.$$
    $A$ 是三对角 M-矩阵（不可约），故 $A^{-1} > 0$。这意味着 $f \ge 0$（源项非负）$\implies$ $u_h \ge 0$（数值解非负），保持了连续问题的**极大值原理**。这是 M-矩阵在 PDE 离散化中的核心应用。

---

## 38.7 迭代法的收敛性

<div class="context-flow" markdown>

**核心问题**：Jacobi 和 Gauss-Seidel 迭代法在什么条件下收敛？M-矩阵保证收敛吗？

</div>

!!! definition "定义 38.9 (正则分裂)"
    矩阵 $A$ 的分裂 $A = M - N$ 称为**正则分裂**，若 $M$ 非奇异、$M^{-1} \ge 0$、$N \ge 0$。

!!! theorem "定理 38.10 (正则分裂的收敛性)"
    设 $A = M - N$ 为正则分裂。则迭代 $x^{(k+1)} = M^{-1}N x^{(k)} + M^{-1}b$ 收敛（对任意初值）当且仅当 $A$ 是非奇异的且 $A^{-1} \ge 0$（即 $A$ 是逆正的）。

    当条件满足时，$\rho(M^{-1}N) < 1$。

??? proof "证明"
    设 $T = M^{-1}N = M^{-1}(M - A) = I - M^{-1}A$。迭代收敛 $\iff$ $\rho(T) < 1$。

    **充分性**：$A$ 逆正，$M^{-1} \ge 0, N \ge 0$。$T = M^{-1}N \ge 0$（非负矩阵）。由 Perron-Frobenius 理论，$\rho(T)$ 是 $T$ 的特征值。需要证明 $\rho(T) < 1$。

    取 $x^* = A^{-1}b \ge 0$（假设 $b \ge 0$）。则 $Ax^* = b$，$Mx^* = Nx^* + b$，$x^* = Tx^* + M^{-1}b$。由 $T \ge 0$ 和 Perron-Frobenius，$\rho(T) \le \max_i (Tx^*)_i / x^*_i$（当 $x^* > 0$ 时）。但 $(Tx^*)_i = x^*_i - (M^{-1}b)_i < x^*_i$（当 $M^{-1}b > 0$），故 $\rho(T) < 1$。

    一般情形的严格证明见 Varga 的著作。

!!! theorem "定理 38.11 (M-矩阵的 Jacobi 和 Gauss-Seidel 收敛)"
    设 $A$ 是非奇异 M-矩阵。则：

    (a) Jacobi 迭代收敛：$\rho(D^{-1}(L+U)) < 1$；

    (b) Gauss-Seidel 迭代收敛：$\rho((D-L)^{-1}U) < 1$；

    (c) 若 $A$ 对称（从而正定——因为对称 M-矩阵的特征值为正），Gauss-Seidel 收敛速度至少是 Jacobi 的两倍（对谱半径而言）。

??? proof "证明（(a)的证明）"
    $A = D - (L + U)$，$D$ 为 $A$ 的对角部分（$d_{ii} = a_{ii} > 0$），$L + U = D - A$。由 $A$ 是 M-矩阵，$D > 0$，$L + U \ge 0$（因为 $A$ 的非对角元素 $\le 0$，故 $D - A$ 的非对角元素 $= -a_{ij} \ge 0$）。

    分裂 $A = D - (L+U)$ 是正则分裂：$D^{-1} \ge 0$（对角正矩阵的逆），$L + U \ge 0$。由定理 38.10 和 $A$ 逆正（M-矩阵的性质 M3），$\rho(D^{-1}(L+U)) < 1$。

!!! example "例 38.8"
    对 $A = \begin{pmatrix} 4 & -1 \\ -1 & 4 \end{pmatrix}$（M-矩阵），Jacobi 迭代矩阵
    $$T_J = D^{-1}(L+U) = \begin{pmatrix} 1/4 & 0 \\ 0 & 1/4 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 1/4 \\ 1/4 & 0 \end{pmatrix}.$$
    $\rho(T_J) = 1/4 < 1$，收敛。

    Gauss-Seidel 迭代矩阵
    $$T_{GS} = (D - L)^{-1}U = \begin{pmatrix} 4 & 0 \\ -1 & 4 \end{pmatrix}^{-1}\begin{pmatrix} 0 & -1 \\ 0 & 0 \end{pmatrix}^{(-)} .$$
    注意此处 $L, U$ 是 $A$ 的严格下、上三角部分的**负**。修正后：$D - L = \begin{pmatrix} 4 & 0 \\ -1 & 4 \end{pmatrix}$（这里 $L$ 的元素为 $-(-1) = 1$，则 $D - L = \begin{pmatrix} 4 & 0 \\ -1 & 4 \end{pmatrix}$），$U$ 的非零元素为 $1$（来自 $-a_{12} = 1$），故 $U$ 矩阵为 $\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。

    $T_{GS} = \begin{pmatrix} 4 & 0 \\ -1 & 4 \end{pmatrix}^{-1}\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} = \frac{1}{16}\begin{pmatrix} 4 & 0 \\ 1 & 4 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 1/4 \\ 0 & 1/16 \end{pmatrix}.$

    $\rho(T_{GS}) = 1/16 = \rho(T_J)^2$，验证了 Gauss-Seidel 收敛速度是 Jacobi 的两倍（在谱半径的对数意义下）。

---

## 38.8 应用：Leontief 投入产出模型

<div class="context-flow" markdown>

**核心问题**：经济系统何时能满足所有部门的最终需求？M-矩阵如何刻画经济的"可行性"？

</div>

!!! definition "定义 38.10 (Leontief 投入产出模型)"
    设经济由 $n$ 个产业部门组成。$a_{ij} \ge 0$ 表示第 $j$ 部门生产单位产品所需第 $i$ 部门的投入量。矩阵 $C = (a_{ij}) \ge 0$ 称为**投入系数矩阵**（或技术矩阵）。若 $x = (x_1, \ldots, x_n)^T$ 为各部门的总产出向量，$d = (d_1, \ldots, d_n)^T \ge 0$ 为最终需求向量，则平衡方程为
    $$x = Cx + d,$$
    即 $(I - C)x = d$。

!!! theorem "定理 38.12 (Hawkins-Simon 条件)"
    设 $C \ge 0$ 为投入系数矩阵。以下条件等价：

    (a) 对每个最终需求 $d \ge 0$，平衡方程 $(I - C)x = d$ 有唯一非负解 $x \ge 0$。

    (b) $I - C$ 是非奇异 M-矩阵。

    (c) **Hawkins-Simon 条件**：$I - C$ 的所有顺序主子式为正：
    $$\det((I-C)[1:k, 1:k]) > 0, \quad k = 1, 2, \ldots, n.$$

    (d) $\rho(C) < 1$（投入系数矩阵的谱半径小于 1）。

    (e) $(I - C)^{-1} \ge 0$（**Leontief 逆**非负）。

??? proof "证明"
    $I - C$ 是 Z-矩阵（非对角元素 $-a_{ij} \le 0$），且 $I - C = I - C$（即 $s = 1, B = C$）。

    **(a) $\Leftrightarrow$ (e)**：$x = (I-C)^{-1}d \ge 0$ 对所有 $d \ge 0$ 当且仅当 $(I-C)^{-1} \ge 0$。

    **(b) $\Leftrightarrow$ (d)**：$I - C = I - C$ 是 M-矩阵 $\iff$ $1 > \rho(C)$。

    **(b) $\Leftrightarrow$ (e)**：M-矩阵条件 M3。

    **(b) $\Leftrightarrow$ (c)**：M-矩阵条件 M4。

!!! example "例 38.9"
    三部门经济，投入系数矩阵
    $$C = \begin{pmatrix} 0.2 & 0.3 & 0.1 \\ 0.1 & 0.1 & 0.2 \\ 0.2 & 0.1 & 0.1 \end{pmatrix}.$$

    $I - C = \begin{pmatrix} 0.8 & -0.3 & -0.1 \\ -0.1 & 0.9 & -0.2 \\ -0.2 & -0.1 & 0.9 \end{pmatrix}$。

    Hawkins-Simon 检验：$\Delta_1 = 0.8 > 0$，$\Delta_2 = 0.72 - 0.03 = 0.69 > 0$，$\Delta_3 = \det(I-C)$。

    $\det(I-C) = 0.8(0.81 - 0.02) - (-0.3)(-0.09 - 0.04) + (-0.1)(0.01 + 0.18)$
    $= 0.8 \times 0.79 - 0.3 \times 0.13 - 0.1 \times 0.19 = 0.632 - 0.039 - 0.019 = 0.574 > 0$。

    所有顺序主子式为正，$I - C$ 是 M-矩阵。经济可以满足任意非负最终需求。

    **经济解读**：Hawkins-Simon 条件 $\Delta_k > 0$ 意味着前 $k$ 个部门的"子经济"自身是有生产能力的——它们的总产出能覆盖相互之间的投入需求并有剩余。这是经济可行性的自然条件。

!!! theorem "定理 38.13 (Leontief 乘数)"
    若 $I - C$ 是 M-矩阵，则 Leontief 逆 $(I - C)^{-1} = \sum_{k=0}^{\infty} C^k$ 的 $(i,j)$ 元素 $l_{ij}$ 表示第 $j$ 部门的最终需求增加 1 单位时，第 $i$ 部门的总产出增加量。$l_{ij} \ge \delta_{ij}$（至少等于 Kronecker $\delta$），反映了经济的**乘数效应**。

??? proof "证明"
    $x = (I-C)^{-1}d$，$\partial x_i / \partial d_j = [(I-C)^{-1}]_{ij} = l_{ij}$。$(I-C)^{-1} = I + C + C^2 + \cdots \ge I$，故 $l_{ij} \ge \delta_{ij}$。

!!! note "注记 38.2 (M-矩阵在其他领域的出现)"
    除了经济学的 Leontief 模型外，M-矩阵还出现在以下领域：

    - **偏微分方程**：椭圆型 PDE 的有限差分和有限元离散化（在适当的网格条件下）产生 M-矩阵，保证离散极大值原理。
    - **马尔可夫链**：率矩阵（生成元）$Q$ 的负（$-Q$）是奇异 M-矩阵。
    - **种群动力学**：竞争 Lotka-Volterra 系统的群落矩阵在一定条件下是 M-矩阵。
    - **电网分析**：节点导纳矩阵在无互感时是 M-矩阵。

---

## 本章小结

| 矩阵类 | 定义 | 核心性质 | 关系 |
|--------|------|----------|------|
| Z-矩阵 | 非对角元 $\le 0$ | $A = sI - B$（$B \ge 0$） | 最大类 |
| M-矩阵 | Z-矩阵 + $s > \rho(B)$ | $A^{-1} \ge 0$；所有主子式 $> 0$；50 个等价条件 | $\subset$ Z-矩阵 |
| P-矩阵 | 所有主子式 $> 0$ | LCP 唯一可解 | $\supset$ M-矩阵 |
| H-矩阵 | 比较矩阵是 M-矩阵 | 广义对角占优 | M-矩阵的推广 |
| 逆正矩阵 | $A^{-1} \ge 0$ | 等价于单调矩阵 | M-矩阵 $\subset$ 逆正 |

这些矩阵类之间的包含关系为：
$$\text{严格对角占优 Z-矩阵} \subset \text{M-矩阵} \subset \text{Z-矩阵} \cap \text{P-矩阵} \subset \text{P-矩阵}$$
$$\text{M-矩阵} \subset \text{逆正矩阵}, \quad \text{对角占优矩阵} \subset \text{H-矩阵}.$$
