# 第 61 章 逆特征值问题

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 正定矩阵(Ch16) · 非负矩阵(Ch17) · Jordan形(Ch12)

**本章脉络**：逆特征值问题的一般框架 → 对称逆特征值问题 → 非负逆特征值问题(NIEP) → Jacobi 逆问题 → Toeplitz 逆问题 → 随机矩阵的逆问题 → 数值方法 → 开放问题

**延伸**：逆特征值问题在结构工程（质量-弹簧系统的设计和模型修正）、控制理论（极点配置的矩阵实现）和分子光谱学（从光谱数据重构分子结构）中有直接应用

</div>

逆特征值问题（Inverse Eigenvalue Problem, IEP）是矩阵理论中一类基本问题：给定谱数据（特征值、部分特征向量等），构造满足特定结构约束的矩阵。与正问题（给定矩阵求特征值）相反，逆问题从"果"推"因"，通常更困难，且涉及存在性、唯一性和数值构造三个层面的挑战。

逆特征值问题的研究可以追溯到 19 世纪 Sturm-Liouville 理论对微分算子逆问题的研究。在有限维情形中，问题的组合与代数结构使得它既具有深刻的理论内涵，也有丰富的应用背景——从弹簧-质量系统的设计到分子光谱学中的结构重构。

---

## 61.1 逆特征值问题的一般框架

<div class="context-flow" markdown>

**核心问题**：逆特征值问题的一般形式是什么？有哪些主要变体？

</div>

!!! definition "定义 61.1 (逆特征值问题)"
    **逆特征值问题**的一般形式为：给定谱数据 $\Lambda$ 和矩阵结构类 $\mathcal{S}$，判断是否存在 $A \in \mathcal{S}$ 使得 $\sigma(A) = \Lambda$。若存在，给出构造方法。

    形式化地：给定 $\Lambda = \{\lambda_1, \ldots, \lambda_n\} \subset \mathbb{C}$（含重数）和结构约束集合 $\mathcal{S} \subset M_n(\mathbb{F})$，求

    $$A \in \mathcal{S} \quad \text{s.t.} \quad \sigma(A) = \Lambda.$$

!!! definition "定义 61.2 (IEP 的主要变体)"
    根据结构约束 $\mathcal{S}$ 的不同，逆特征值问题有以下主要变体：

    **(a) 对称 IEP**（SIEP）：$\mathcal{S}$ = 实对称矩阵。

    **(b) 非负 IEP**（NIEP）：$\mathcal{S}$ = 非负矩阵（所有元素 $\geq 0$）。

    **(c) 随机 IEP**：$\mathcal{S}$ = 行随机矩阵（每行元素之和为 1，元素非负）。

    **(d) Jacobi IEP**：$\mathcal{S}$ = Jacobi 矩阵（对称三对角，正的次对角线元素）。

    **(e) Toeplitz IEP**：$\mathcal{S}$ = Toeplitz 矩阵。

    **(f) 带结构 IEP**：$\mathcal{S}$ = 具有规定稀疏模式的矩阵。

    **(g) 部分 IEP**：只给定部分谱数据（部分特征值和/或部分特征向量信息）。

!!! definition "定义 61.3 (IEP 的三个层面)"
    逆特征值问题的研究涉及三个层面：

    - **存在性**：给定 $\Lambda$，是否存在 $A \in \mathcal{S}$ 具有谱 $\Lambda$？
    - **唯一性**：如果存在，是否唯一？需要多少额外数据来确定唯一解？
    - **构造性**：如何高效地计算出满足条件的 $A$？

!!! example "例 61.1"
    **最简单的 IEP。** 给定 $\lambda_1, \ldots, \lambda_n \in \mathbb{C}$，构造一个矩阵 $A$ 使得 $\sigma(A) = \{\lambda_1, \ldots, \lambda_n\}$。

    无结构约束时，这是平凡的：取 $A = \mathrm{diag}(\lambda_1, \ldots, \lambda_n)$。

    问题的困难完全来自结构约束 $\mathcal{S}$。

---

## 61.2 对称逆特征值问题

<div class="context-flow" markdown>

**核心问题**：给定实数 $\lambda_1 \leq \cdots \leq \lambda_n$，何时能构造具有这些特征值的实对称矩阵？如果还有额外的结构约束呢？

</div>

!!! theorem "定理 61.1 (对称 IEP 的可解性)"
    给定任意实数 $\lambda_1, \ldots, \lambda_n \in \mathbb{R}$，总存在 $n \times n$ 实对称矩阵 $A$ 使得 $\sigma(A) = \{\lambda_1, \ldots, \lambda_n\}$。

??? proof "证明"
    取任意正交矩阵 $Q \in O(n)$（例如 $Q = I$），令

    $$A = Q \,\mathrm{diag}(\lambda_1, \ldots, \lambda_n)\, Q^\top.$$

    则 $A$ 是实对称矩阵（$A^\top = Q\Lambda Q^\top = A$），且 $\sigma(A) = \{\lambda_1, \ldots, \lambda_n\}$。$\blacksquare$

这个无约束情形是平凡的。困难在于带有**额外约束**的对称 IEP。

!!! definition "定义 61.4 (带规定对角元的对称 IEP)"
    给定实数 $\lambda_1, \ldots, \lambda_n$ 和 $d_1, \ldots, d_n$，构造实对称矩阵 $A$ 使得

    $$\sigma(A) = \{\lambda_1, \ldots, \lambda_n\}, \quad a_{ii} = d_i, \quad i = 1, \ldots, n.$$

!!! theorem "定理 61.2 (Schur-Horn 定理)"
    给定 $\lambda_1 \geq \cdots \geq \lambda_n$ 和 $d_1 \geq \cdots \geq d_n$，存在实对称矩阵 $A$ 使得 $\sigma(A) = \{\lambda_i\}$ 且 $\mathrm{diag}(A) = (d_1, \ldots, d_n)^\top$，当且仅当 $d$ 被 $\lambda$ **优超**（majorized by）：

    $$\sum_{i=1}^k d_i \leq \sum_{i=1}^k \lambda_i, \quad k = 1, \ldots, n-1,$$

    $$\sum_{i=1}^n d_i = \sum_{i=1}^n \lambda_i.$$

??? proof "证明"
    **必要性。** 设 $A = Q\Lambda Q^\top$，$Q$ 正交。则 $a_{ii} = d_i = \sum_j q_{ij}^2 \lambda_j$。由于 $(q_{i1}^2, \ldots, q_{in}^2)$ 是概率向量，$d_i$ 是 $\lambda_1, \ldots, \lambda_n$ 的凸组合。由 Schur 的定理，对角向量 $d$ 属于 $\lambda$ 的置换多面体（permutahedron），即 $d \prec \lambda$（$d$ 被 $\lambda$ 优超）。

    **充分性。** Horn (1954) 通过构造性证明完成。核心工具是 Givens 旋转：如果 $\lambda = (\lambda_1, \ldots, \lambda_n)$ 且 $d \prec \lambda$，可以通过一系列 $2 \times 2$ 旋转将 $\mathrm{diag}(\lambda)$ 变换为具有规定对角元的对称矩阵。

    具体地，给定两个特征值 $\lambda_i > \lambda_j$ 和两个目标对角元 $d_i, d_j$，满足 $d_i + d_j = \lambda_i + \lambda_j$ 且 $\lambda_j \leq d_i, d_j \leq \lambda_i$，可以找到角度 $\theta$ 使得 Givens 旋转 $G_{ij}(\theta)$ 将对角矩阵的第 $i, j$ 个对角元变为 $d_i, d_j$：

    $$\begin{pmatrix}\cos\theta & -\sin\theta\\\sin\theta & \cos\theta\end{pmatrix}\begin{pmatrix}\lambda_i & 0\\0 & \lambda_j\end{pmatrix}\begin{pmatrix}\cos\theta & \sin\theta\\-\sin\theta & \cos\theta\end{pmatrix} = \begin{pmatrix}d_i & *\\* & d_j\end{pmatrix}.$$

    逐步进行这样的旋转即可达到目标。$\blacksquare$

!!! example "例 61.2"
    $\lambda = (5, 3, 1)$，$d = (4, 3, 2)$。验证优超条件：

    - $4 \leq 5$：成立。
    - $4+3 = 7 \leq 5+3 = 8$：成立。
    - $4+3+2 = 9 = 5+3+1 = 9$：成立。

    因此存在 $3 \times 3$ 实对称矩阵，特征值为 $5,3,1$，对角元为 $4,3,2$。

!!! theorem "定理 61.3 (带结构约束的对称 IEP)"
    **带规定带宽的对称 IEP** 是更困难的问题。设 $\mathcal{S}$ 为 $n \times n$ 带宽 $\beta$ 的实对称矩阵集合（$a_{ij} = 0$ 当 $|i-j| > \beta$）。

    给定 $\lambda_1, \ldots, \lambda_n \in \mathbb{R}$，当 $\beta \geq n-1$ 时（无额外约束），问题总有解。当 $\beta = 1$（三对角）时，问题归结为 Jacobi 逆问题（见 61.4 节）。当 $1 < \beta < n-1$ 时，存在性条件不完全清楚。

---

## 61.3 非负逆特征值问题 (NIEP)

<div class="context-flow" markdown>

**核心问题**：给定复数 $\lambda_1, \ldots, \lambda_n$，何时存在 $n \times n$ 非负矩阵具有这些特征值？这一问题为何如此困难？

</div>

!!! definition "定义 61.5 (非负逆特征值问题)"
    **非负逆特征值问题**（NIEP）：给定 $\Lambda = \{\lambda_1, \ldots, \lambda_n\} \subset \mathbb{C}$（含重数），判断是否存在非负矩阵 $A \in \mathbb{R}_{\geq 0}^{n \times n}$ 使得 $\sigma(A) = \Lambda$。

    如果存在这样的 $A$，则称 $\Lambda$ 是**可实现的**（realizable）。

!!! theorem "定理 61.4 (NIEP 的必要条件)"
    如果 $\Lambda = \{\lambda_1, \ldots, \lambda_n\}$ 是可实现的（$\lambda_1 \geq |\lambda_i|$，$\forall i$），则：

    **(a) Perron 条件**：$\lambda_1 \in \mathbb{R}$，$\lambda_1 \geq 0$，且 $\lambda_1 \geq |\lambda_i|$，$\forall i$。

    **(b) 迹条件**：$s_k = \sum_{i=1}^n \lambda_i^k \geq 0$，$\forall k \geq 1$（即所有幂和为非负实数）。

    **(c) 共轭对称**：非实特征值成共轭对出现。

    **(d) Newton 不等式**：$s_k^2 \leq n \cdot s_{2k}$，$\forall k \geq 1$。

    **(e) Loewy-London 条件**（$n = 3$）：$s_1^3 + 2s_3 \geq 3s_1 s_2$。

??? proof "证明"
    **(a)** 由 Perron-Frobenius 定理，非负矩阵的谱半径是特征值（Perron 根），且为最大模特征值。

    **(b)** $s_k = \mathrm{tr}(A^k)$。$A^k$ 是非负矩阵（非负矩阵的幂仍非负），其迹为对角元之和，非负。

    **(c)** 非负矩阵是实矩阵，实矩阵的特征多项式系数为实数，因此非实特征值成共轭对。

    **(d)** 由 Cauchy-Schwarz 不等式：$s_k^2 = (\sum \lambda_i^k)^2 \leq n \sum |\lambda_i|^{2k} \leq n \sum \lambda_i^{2k} = n \cdot s_{2k}$（利用 $|\lambda_i^k|^2 = |\lambda_i|^{2k} \leq \lambda_i^{2k}$，后者在非负矩阵情形需要更仔细的论证）。

    实际上更精确的论证是：$s_k = \mathrm{tr}(A^k)$，由非负矩阵的性质，$s_{2k} = \mathrm{tr}(A^{2k}) = \mathrm{tr}((A^k)^\top A^k) = \|A^k\|_F^2 \geq 0$，且 $s_k^2 = (\mathrm{tr}(A^k))^2 \leq n \cdot \mathrm{tr}((A^k)^2) = n \cdot s_{2k}$（后一步是 Cauchy-Schwarz 对迹的应用）。$\blacksquare$

!!! theorem "定理 61.5 (NIEP 的充分条件——Suleimanova)"
    设 $\lambda_1 > 0 > \lambda_2 \geq \cdots \geq \lambda_n$（实特征值，一正 $n-1$ 负）。如果 $\sum_{i=1}^n \lambda_i \geq 0$，则 $\Lambda$ 是可实现的。

??? proof "证明"
    构造**伴随矩阵**（companion matrix）形式的非负矩阵。令

    $$A = \begin{pmatrix}0 & 0 & \cdots & 0 & c_n\\1 & 0 & \cdots & 0 & c_{n-1}\\0 & 1 & \cdots & 0 & c_{n-2}\\\vdots & & \ddots & & \vdots\\0 & 0 & \cdots & 1 & c_1\end{pmatrix},$$

    其特征多项式为 $p(\lambda) = \lambda^n - c_1\lambda^{n-1} - \cdots - c_n$。

    当 $\lambda_i$ 满足条件时，Newton 恒等式保证 $c_1 = \sum \lambda_i \geq 0$，且通过精心选择可以使所有 $c_i \geq 0$。详细验证需要用到实特征值的 Vieta 关系和符号条件。

    Suleimanova (1949) 的原始构造使用了秩-1 修正：取 $A = \mathrm{diag}(\lambda_2, \ldots, \lambda_n) + \mathbf{1}v^\top$，选择适当的 $v$ 使得 $A$ 非负且新增特征值为 $\lambda_1$。$\blacksquare$

!!! theorem "定理 61.6 (NIEP 对 $n \leq 4$ 完全解决)"
    **(a)** $n = 1$：$\lambda_1 \geq 0$ 当且仅当可实现。

    **(b)** $n = 2$：$\{\lambda_1, \lambda_2\}$ 可实现当且仅当 $\lambda_1 \geq 0$，$\lambda_1 \geq |\lambda_2|$，$\lambda_1 + \lambda_2 \geq 0$。

    **(c)** $n = 3$：Loewy 和 London (1978) 给出了完整的充要条件。

    **(d)** $n = 4$：Meehan (1998) 和其他人给出了完整刻画。

    **(e)** $n \geq 5$：NIEP **仍是开放问题**。

!!! example "例 61.3"
    **$n = 3$ 的具体验证与构造。** $\Lambda = \{4, -1, -2\}$。

    **第一步：验证必要条件。**

    - Perron 条件：$4 \geq |-1|, |-2|$，成立。
    - $s_1 = 4 + (-1) + (-2) = 1 \geq 0$。
    - $s_2 = 16 + 1 + 4 = 21 \geq 0$。
    - $s_3 = 64 + (-1) + (-8) = 55 \geq 0$。
    - Loewy-London：$s_1^3 + 2s_3 = 1 + 110 = 111 \geq 3 s_1 s_2 = 63$。成立。

    所有必要条件满足。又由 Suleimanova 充分条件（一正两负，$s_1 = 1 > 0$），$\Lambda$ 可实现。

    **第二步：Suleimanova 秩-1 构造。** 取对角矩阵 $D = \mathrm{diag}(-1, -2, 0)$，令

    $$A = D + \mathbf{e} v^\top = \begin{pmatrix}-1 & 0 & 0\\0 & -2 & 0\\0 & 0 & 0\end{pmatrix} + \begin{pmatrix}1\\1\\1\end{pmatrix}\begin{pmatrix}a & b & c\end{pmatrix} = \begin{pmatrix}a-1 & b & c\\a & b-2 & c\\a & b & c\end{pmatrix}.$$

    秩-1 修正的特征值分析：$D$ 的特征值为 $\{-1, -2, 0\}$，秩-1 修正 $\mathbf{e}v^\top$ 保留 $\ker(v^\top)$ 中的特征值不变。具体地，$v^\top$ 对 $D$ 的特征向量 $e_i$ 的作用为 $v^\top e_i = v_i$，因此修正后的特征值由**世俗方程**（secular equation）决定：

    $$1 = \sum_{i=1}^{3} \frac{v_i}{\lambda - d_i} = \frac{a}{\lambda + 1} + \frac{b}{\lambda + 2} + \frac{c}{\lambda}.$$

    我们需要此方程的根为 $\lambda = 4, -1, -2$。但这里 $d_1 = -1, d_2 = -2$ 与目标特征值 $-1, -2$ 重合，这意味着需要 $v_1 = a = 0$ 或 $v_2 = b = 0$ 来保留对应特征值。

    实际上更直接的方法：令 $D$ 的对角元素为我们要保留的特征值。取 $D = \mathrm{diag}(-1, -2, 0)$，要保留 $-1$ 和 $-2$，令 $v = (0, 0, c)^\top$（使得 $v^\top e_1 = v^\top e_2 = 0$，即 $-1, -2$ 对应的特征向量不受影响）。则

    $$A = \begin{pmatrix}-1 & 0 & c\\0 & -2 & c\\0 & 0 & c\end{pmatrix}.$$

    此上三角矩阵的特征值为对角元素 $\{-1, -2, c\}$。令 $c = 4$：

    $$A = \begin{pmatrix}-1 & 0 & 4\\0 & -2 & 4\\0 & 0 & 4\end{pmatrix}.$$

    但 $A$ 含负元素，不是非负矩阵。

    **第三步：使用正确的 Suleimanova 构造。** 令 $D = \mathrm{diag}(\lambda_2, \lambda_3) = \mathrm{diag}(-1, -2)$，构造

    $$A = \begin{pmatrix}-1 & 0 & 0\\0 & -2 & 0\\0 & 0 & 0\end{pmatrix} + \begin{pmatrix}1\\1\\1\end{pmatrix}\begin{pmatrix}a & b & c\end{pmatrix}$$

    并要求 (i) $A$ 非负，(ii) $\sigma(A) = \{4, -1, -2\}$。

    条件 (ii)：$\mathrm{tr}(A) = (a-1) + (b-2) + c = a + b + c - 3 = 4 + (-1) + (-2) = 1$，所以 $a + b + c = 4$。

    秩-1 修正理论：$D' = \mathrm{diag}(-1, -2, 0)$ 经过秩-1 修正 $\mathbf{e}v^\top$（$\mathbf{e} = (1,1,1)^\top$，$v = (a,b,c)^\top$）后，新特征值满足世俗方程。由于 $\mathbf{e}$ 不是 $D'$ 的特征向量，需要更仔细的分析。

    换一种方法：直接选取参数并验证。取 $a = 1, b = 2, c = 1$（满足 $a+b+c = 4$）：

    $$A = \begin{pmatrix}0 & 2 & 1\\1 & 0 & 1\\1 & 2 & 1\end{pmatrix}.$$

    验证非负性：所有元素 $\geq 0$，成立。

    特征多项式：$\det(\lambda I - A) = \det\begin{pmatrix}\lambda & -2 & -1\\-1 & \lambda & -1\\-1 & -2 & \lambda-1\end{pmatrix}$。

    展开：$\lambda[\lambda(\lambda-1) - (-1)(-2)] - (-2)[(-1)(\lambda-1) - (-1)(-1)] + (-1)[(-1)(-2) - \lambda(-1)]$

    $= \lambda[\lambda^2 - \lambda - 2] + 2[-\lambda + 1 - 1] + (-1)[2 + \lambda]$

    $= \lambda^3 - \lambda^2 - 2\lambda + 2(-\lambda) + (-2 - \lambda)$

    $= \lambda^3 - \lambda^2 - 2\lambda - 2\lambda - 2 - \lambda$

    $= \lambda^3 - \lambda^2 - 5\lambda - 2$。

    检验 $\lambda = 4$：$64 - 16 - 20 - 2 = 26 \neq 0$。不正确。

    **第四步：系统性构造。** 使用伴随矩阵（companion matrix）方法。目标特征多项式为

    $$p(\lambda) = (\lambda - 4)(\lambda + 1)(\lambda + 2) = \lambda^3 - \lambda^2 - 6\lambda - 8.$$

    展开验证：$(\lambda - 4)(\lambda + 1) = \lambda^2 - 3\lambda - 4$，$(\lambda^2 - 3\lambda - 4)(\lambda + 2) = \lambda^3 + 2\lambda^2 - 3\lambda^2 - 6\lambda - 4\lambda - 8 = \lambda^3 - \lambda^2 - 10\lambda - 8$。

    修正：$(\lambda^2 - 3\lambda - 4)(\lambda + 2) = \lambda^3 - 3\lambda^2 - 4\lambda + 2\lambda^2 - 6\lambda - 8 = \lambda^3 - \lambda^2 - 10\lambda - 8$。

    因此 $p(\lambda) = \lambda^3 - \lambda^2 - 10\lambda - 8$，$s_1 = 1, s_2 = 21, s_3 = 55$，迹为 $1$。

    非负伴随矩阵：

    $$C = \begin{pmatrix}0 & 0 & 8\\1 & 0 & 10\\0 & 1 & 1\end{pmatrix}.$$

    特征多项式为 $\lambda^3 - 1 \cdot \lambda^2 - 10\lambda - 8 = p(\lambda)$。所有元素非负。特征值 $\{4, -1, -2\}$。

    验证 $\lambda = 4$：$64 - 16 - 40 - 8 = 0$。$\checkmark$

    验证 $\lambda = -1$：$-1 - 1 + 10 - 8 = 0$。$\checkmark$

    验证 $\lambda = -2$：$-8 - 4 + 20 - 8 = 0$。$\checkmark$

    因此 $A = C = \begin{pmatrix}0 & 0 & 8\\1 & 0 & 10\\0 & 1 & 1\end{pmatrix}$ 是一个 $3 \times 3$ 非负矩阵，特征值恰为 $\{4, -1, -2\}$。

---

## 61.4 Jacobi 逆问题

<div class="context-flow" markdown>

**核心问题**：如何从两个交替的谱恢复一个 Jacobi 矩阵？这与质量-弹簧系统有什么关系？

</div>

!!! definition "定义 61.6 (Jacobi 矩阵)"
    **Jacobi 矩阵**是 $n \times n$ 实对称三对角矩阵，次对角线元素为正：

    $$J = \begin{pmatrix}a_1 & b_1 & & \\b_1 & a_2 & b_2 & \\ & \ddots & \ddots & \ddots \\ & & b_{n-1} & a_n\end{pmatrix}, \quad b_i > 0.$$

!!! theorem "定理 61.7 (Jacobi 矩阵的特征值交错)"
    设 $J_n$ 是 $n \times n$ Jacobi 矩阵，$J_{n-1}$ 是删去最后一行和最后一列得到的 $(n-1) \times (n-1)$ 前导子矩阵。若 $\sigma(J_n) = \{\lambda_1 < \cdots < \lambda_n\}$，$\sigma(J_{n-1}) = \{\mu_1 < \cdots < \mu_{n-1}\}$，则特征值**严格交错**：

    $$\lambda_1 < \mu_1 < \lambda_2 < \mu_2 < \cdots < \mu_{n-1} < \lambda_n.$$

??? proof "证明"
    这是 Cauchy 交错定理的特殊情形，对对称三对角矩阵严格成立（因为 $b_i > 0$ 保证了不出现相等的情况）。

    设 $p_k(\lambda) = \det(\lambda I_k - J_k)$ 为 $J_k$ 的特征多项式。三对角结构给出递推关系

    $$p_k(\lambda) = (\lambda - a_k)p_{k-1}(\lambda) - b_{k-1}^2 p_{k-2}(\lambda).$$

    设 $\lambda_i$ 是 $p_n$ 的根。由递推关系，$p_{n-1}(\lambda_i) = -\frac{p_n(\lambda_i)}{(\cdots)} + \cdots$，可以证明 $p_{n-1}$ 在 $J_n$ 的相邻特征值之间恰好变号一次，因此 $J_{n-1}$ 的特征值严格交错于 $J_n$ 的特征值之间。$\blacksquare$

!!! theorem "定理 61.8 (Hochstadt 唯一性定理, 1974)"
    给定两组交错的实数

    $$\lambda_1 < \mu_1 < \lambda_2 < \mu_2 < \cdots < \mu_{n-1} < \lambda_n,$$

    存在**唯一**的 $n \times n$ Jacobi 矩阵 $J$，使得 $\sigma(J) = \{\lambda_i\}$ 且 $\sigma(J_{n-1}) = \{\mu_i\}$（$J_{n-1}$ 为前导子矩阵）。

??? proof "证明"
    **构造算法（de Boor-Golub 算法）。** 利用 Lanczos 三对角化过程的逆过程。

    **第一步：从谱数据构造谱测度。** 定义离散测度

    $$d\alpha(\lambda) = \sum_{i=1}^n w_i \delta(\lambda - \lambda_i),$$

    其中权重 $w_i > 0$ 由条件 $\sigma(J_{n-1}) = \{\mu_j\}$ 唯一确定：

    $$w_i = \frac{\prod_{j=1}^{n-1}(\lambda_i - \mu_j)}{\prod_{j \neq i}(\lambda_i - \lambda_j)}.$$

    交错条件保证 $w_i > 0$。

    **第二步：正交多项式。** 对测度 $d\alpha$ 进行 Gram-Schmidt 正交化，得到正交多项式 $p_0, p_1, \ldots, p_{n-1}$。

    **第三步：三项递推给出 Jacobi 矩阵。** 正交多项式满足三项递推

    $$\lambda p_k(\lambda) = b_k p_{k+1}(\lambda) + a_{k+1} p_k(\lambda) + b_{k-1} p_{k-1}(\lambda),$$

    递推系数 $(a_i, b_i)$ 就是 Jacobi 矩阵的元素。唯一性由正交多项式的唯一性保证。$\blacksquare$

!!! example "例 61.4"
    **质量-弹簧系统。** 考虑 $n$ 个质量通过弹簧串联连接，一端固定：

    $$m_1 \xleftrightarrow{k_1} m_2 \xleftrightarrow{k_2} \cdots \xleftrightarrow{k_{n-1}} m_n \xleftrightarrow{k_n} \text{墙}$$

    无阻尼自由振动方程 $M\ddot{u} + Ku = 0$（$M$ 为质量矩阵，$K$ 为刚度矩阵）。令 $A = M^{-1/2}KM^{-1/2}$，则 $A$ 是 Jacobi 矩阵，特征值 $\omega_i^2$ 是振动频率的平方。

    **逆问题**：从测量到的振动频率 $\omega_1, \ldots, \omega_n$（$J_n$ 的特征值）和移除最后一个质量后的频率 $\omega'_1, \ldots, \omega'_{n-1}$（$J_{n-1}$ 的特征值），唯一确定系统参数 $m_i, k_i$。

---

## 61.5 Toeplitz 逆问题

<div class="context-flow" markdown>

**核心问题**：给定特征值，何时能构造具有这些特征值的 Toeplitz 矩阵？

</div>

!!! definition "定义 61.7 (Toeplitz 矩阵)"
    **Toeplitz 矩阵**是沿对角线元素相同的矩阵：

    $$T = \begin{pmatrix}t_0 & t_{-1} & t_{-2} & \cdots\\t_1 & t_0 & t_{-1} & \cdots\\t_2 & t_1 & t_0 & \cdots\\\vdots & & & \ddots\end{pmatrix}.$$

    实对称 Toeplitz 矩阵形如 $T_{ij} = t_{|i-j|}$，由 $n$ 个参数 $t_0, t_1, \ldots, t_{n-1}$ 确定。

!!! theorem "定理 61.9 (Toeplitz IEP 的必要条件——Landau)"
    给定 $\lambda_1 \leq \cdots \leq \lambda_n$，如果存在 $n \times n$ 实对称 Toeplitz 矩阵 $T$ 使得 $\sigma(T) = \{\lambda_i\}$，则

    **(a)** $\mathrm{tr}(T) = nt_0 = \sum \lambda_i$，确定 $t_0 = \frac{1}{n}\sum \lambda_i$。

    **(b)** $\mathrm{tr}(T^2) = \sum \lambda_i^2$ 给出关于 $t_0, \ldots, t_{n-1}$ 的一个方程。

    **(c)** 更高阶矩 $s_k = \sum \lambda_i^k = \mathrm{tr}(T^k)$ 给出进一步的约束。

    由于 Toeplitz 矩阵有 $n$ 个自由参数（$t_0, \ldots, t_{n-1}$），而规定 $n$ 个特征值也给出 $n$ 个方程（通过 Newton 恒等式），问题在参数计数上是"恰好确定的"。

!!! theorem "定理 61.10 (Toeplitz IEP 的存在性)"
    **(a)** 对 $n \leq 4$，Toeplitz IEP 对所有合法的实谱 $\{\lambda_i\}$ 有解。

    **(b)** 对 $n = 5$，存在反例：某些实谱不能由 $5 \times 5$ 实对称 Toeplitz 矩阵实现。

    **(c)** 一般 $n$ 的完整存在性条件是**开放问题**。

!!! example "例 61.5"
    **$n = 3$ 的 Toeplitz IEP。** 给定 $\lambda_1, \lambda_2, \lambda_3 \in \mathbb{R}$。实对称 Toeplitz 矩阵

    $$T = \begin{pmatrix}a & b & c\\b & a & b\\c & b & a\end{pmatrix}.$$

    特征多项式 $\det(\lambda I - T) = (\lambda - a - b - c)(\lambda - a + c)^2 - \cdots$。实际上 $T$ 的特征值可以用循环矩阵的技巧或直接计算。对于 $3 \times 3$ 对称 Toeplitz，有 3 个参数 $a, b, c$ 和 3 个特征值方程，通常有解。

---

## 61.6 随机矩阵的逆问题

<div class="context-flow" markdown>

**核心问题**：给定特征值，何时能构造行随机矩阵或双随机矩阵？

</div>

!!! definition "定义 61.8 (行随机矩阵的逆特征值问题)"
    **行随机 IEP**：给定 $\Lambda = \{1, \lambda_2, \ldots, \lambda_n\}$（1 是 Perron 根），构造行随机矩阵（非负，每行和为 1）$A$，使得 $\sigma(A) = \Lambda$。

!!! theorem "定理 61.11 (Kellogg 条件)"
    $\Lambda = \{1, \lambda_2, \ldots, \lambda_n\}$ 是 $n \times n$ 行随机矩阵的谱的必要条件包括：

    (a) $|\lambda_i| \leq 1$，$\forall i$。

    (b) 非实特征值成共轭对。

    (c) $s_k = 1 + \sum_{i=2}^n \lambda_i^k \geq 0$，$\forall k \geq 1$。

    (d) 对 $n = 3$：$\Lambda = \{1, \lambda_2, \lambda_3\}$ 可实现当且仅当 $\lambda_2, \lambda_3$ 在 $\mathbb{C}$ 的 Kellogg 区域内。

!!! definition "定义 61.9 (Perfect-Mirsky 区域)"
    对 $n \times n$ 双随机矩阵（非负，每行每列和为 1），可能的谱所构成的集合称为 **Perfect-Mirsky 区域** $\Pi_n$。

    $\Pi_n$ 是 $\mathbb{C}^{n-1}$ 的一个子集（因为一个特征值固定为 1）。对 $n = 2$，$\Pi_2 = [-1, 1]$。对 $n = 3$，$\Pi_3$ 的完整描述由 Perfect 和 Mirsky (1965) 给出。对一般 $n$，$\Pi_n$ 的完整刻画是开放问题。

!!! example "例 61.6"
    **$n = 3$ 的双随机矩阵。** $\Lambda = \{1, -0.5, -0.5\}$。构造

    $$A = \frac{1}{3}\begin{pmatrix}1 & 1 & 1\\1 & 1 & 1\\1 & 1 & 1\end{pmatrix} + \frac{2}{3}\begin{pmatrix}? & ? & ?\\? & ? & ?\\? & ? & ?\end{pmatrix}.$$

    实际上取循环矩阵 $A = \begin{pmatrix}a & b & c\\c & a & b\\b & c & a\end{pmatrix}$，$a + b + c = 1$。特征值为 $1$，$a + b\omega + c\omega^2$，$a + b\omega^2 + c\omega$（$\omega = e^{2\pi i/3}$）。要使后两个特征值都是 $-0.5$（实数），需要 $a + b\omega + c\omega^2 = -0.5$，即 $a - (b+c)/2 + (b-c)\sqrt{3}i/2 = -0.5$。实部条件：$a - (1-a)/2 = -0.5$，$3a/2 - 1/2 = -0.5$，$a = 0$。虚部条件：$b = c$。则 $b = c = 1/2$。

    $$A = \begin{pmatrix}0 & 1/2 & 1/2\\1/2 & 0 & 1/2\\1/2 & 1/2 & 0\end{pmatrix}.$$

    验证：行和 = 1，列和 = 1，非负。特征值 $\{1, -1/2, -1/2\}$。成功。

---

## 61.7 数值方法

<div class="context-flow" markdown>

**核心问题**：如何数值地求解逆特征值问题？有哪些有效的算法？

</div>

!!! definition "定义 61.10 (逆特征值问题的优化表述)"
    许多 IEP 可以表述为优化问题。例如，对称 IEP（带结构约束 $\mathcal{S}$）：

    $$\min_{A \in \mathcal{S}} \sum_{i=1}^n (\lambda_i(A) - \lambda_i^*)^2,$$

    其中 $\lambda_i^*$ 是目标特征值，$\lambda_i(A)$ 是 $A$ 的特征值（升序排列）。

!!! theorem "定理 61.12 (提升-投影方法)"
    **提升-投影**（lift-and-project）方法是求解结构化 IEP 的一类经典方法：

    1. **提升步**：在无结构约束的空间中构造具有目标谱的矩阵 $A = Q\Lambda Q^\top$。
    2. **投影步**：将 $A$ 投影到结构约束集 $\mathcal{S}$ 上：$\hat{A} = \Pi_{\mathcal{S}}(A)$。
    3. **更新步**：更新正交矩阵 $Q$，使得 $Q\Lambda Q^\top$ 更接近 $\mathcal{S}$。
    4. 重复直到收敛。

!!! theorem "定理 61.13 (Newton 方法在流形上)"
    对光滑参数化的 IEP，可以在**矩阵流形**上使用 Newton 方法。

    设结构约束将矩阵参数化为 $A(p)$（$p \in \mathbb{R}^m$ 为参数向量）。定义映射

    $$F(p) = (\lambda_1(A(p)) - \lambda_1^*,\, \ldots,\, \lambda_n(A(p)) - \lambda_n^*) \in \mathbb{R}^n.$$

    Newton 迭代：$p^{(k+1)} = p^{(k)} - J_F(p^{(k)})^{-1} F(p^{(k)})$，其中 $J_F$ 是 Jacobian 矩阵。

    Jacobian 的计算需要特征值对参数的导数，由特征值微扰理论给出：

    $$\frac{\partial \lambda_i}{\partial p_j} = q_i^\top \frac{\partial A}{\partial p_j} q_i,$$

    其中 $q_i$ 是 $\lambda_i$ 对应的单位特征向量。

!!! definition "定义 61.11 (同伦延拓方法)"
    **同伦延拓**（homotopy continuation）方法：构造从已知解到目标问题的连续路径。

    设 $\Lambda_0$ 是一个已知有解的目标谱（$A_0 \in \mathcal{S}$，$\sigma(A_0) = \Lambda_0$），$\Lambda_1$ 是真正的目标谱。定义同伦

    $$H(A, t) = \sigma(A) - [(1-t)\Lambda_0 + t\Lambda_1] = 0, \quad t \in [0, 1].$$

    从 $t = 0$（已知解 $A_0$）出发，沿 $t$ 的增加追踪解曲线 $A(t)$，直到 $t = 1$ 得到目标解。

!!! example "例 61.7"
    **Jacobi IEP 的数值构造。** 给定 $\lambda_1 = 1, \lambda_2 = 3, \lambda_3 = 6$ 和 $\mu_1 = 2, \mu_2 = 4$（交错条件 $1 < 2 < 3 < 4 < 6$ 满足）。

    权重：$w_1 = \frac{(1-2)(1-4)}{(1-3)(1-6)} = \frac{(-1)(-3)}{(-2)(-5)} = \frac{3}{10}$。

    $w_2 = \frac{(3-2)(3-4)}{(3-1)(3-6)} = \frac{(1)(-1)}{(2)(-3)} = \frac{1}{6}$。

    $w_3 = \frac{(6-2)(6-4)}{(6-1)(6-3)} = \frac{(4)(2)}{(5)(3)} = \frac{8}{15}$。

    验证 $w_1 + w_2 + w_3 = 3/10 + 1/6 + 8/15 = 9/30 + 5/30 + 16/30 = 30/30 = 1$。

    然后通过对测度 $\sum w_i \delta(\lambda - \lambda_i)$ 做 Gram-Schmidt 正交化，得到正交多项式的三项递推系数，即为 Jacobi 矩阵的元素。

---

## 61.8 特征值与奇异值的关系

<div class="context-flow" markdown>

**核心问题**：矩阵的特征值和奇异值之间有什么约束关系？给定特征值和奇异值，何时存在相应的矩阵？

</div>

特征值和奇异值是矩阵的两组基本不变量。特征值描述线性变换的"缩放-旋转"行为，奇异值描述其"伸缩"行为。它们之间的关系由一系列深刻的不等式刻画。

!!! theorem "定理 61.14 (Weyl-Horn 定理)"
    设 $A \in M_n(\mathbb{C})$ 的特征值为 $\lambda_1, \ldots, \lambda_n$（按模递减排列：$|\lambda_1| \geq \cdots \geq |\lambda_n|$），奇异值为 $s_1 \geq \cdots \geq s_n \geq 0$。则：

    **(a) Weyl 乘积不等式**：对所有 $k = 1, \ldots, n$，

    $$\prod_{i=1}^k |\lambda_i| \leq \prod_{i=1}^k s_i,$$

    且 $k = n$ 时取等号：$\prod_{i=1}^n |\lambda_i| = \prod_{i=1}^n s_i = |\det(A)|$。

    **(b) Horn 对数优超条件**：设 $\alpha_i = \log|\lambda_i|$，$\beta_i = \log s_i$（当 $\lambda_i \neq 0, s_i > 0$ 时）。则 $\alpha \prec_w \beta$（弱优超），即

    $$\sum_{i=1}^k \log|\lambda_i| \leq \sum_{i=1}^k \log s_i, \quad k = 1, \ldots, n,$$

    且 $k = n$ 时取等号。

    **(c) Horn 的充分性 (1954)**：反过来，给定 $\lambda_1, \ldots, \lambda_n \in \mathbb{C}$ 和 $s_1 \geq \cdots \geq s_n \geq 0$，如果满足上述条件 (a)（所有乘积不等式和最终等式），则存在 $A \in M_n(\mathbb{C})$ 具有特征值 $\{\lambda_i\}$ 和奇异值 $\{s_i\}$。

??? proof "证明"
    **(a) 必要性。** 设 $A = U\Sigma V^*$ 是 $A$ 的奇异值分解，$A = QTQ^*$ 是 Schur 分解（$T$ 上三角，对角元为 $\lambda_i$）。

    考虑 $A$ 的 $k$ 阶复合矩阵（compound matrix）$C_k(A)$。$C_k(A)$ 的特征值为 $\{\lambda_{i_1}\cdots\lambda_{i_k} : i_1 < \cdots < i_k\}$，谱半径为 $|\lambda_1 \cdots \lambda_k|$（因为特征值按模递减排列）。$C_k(A)$ 的算子范数为 $s_1 \cdots s_k$。由谱半径不超过算子范数：

    $$|\lambda_1 \cdots \lambda_k| \leq s_1 \cdots s_k.$$

    等式情形 $k = n$：$|\lambda_1 \cdots \lambda_n| = |\det(A)| = s_1 \cdots s_n$。

    **(c) 充分性（思路）。** Horn 的证明是构造性的。核心步骤：

    1. 先处理 $n = 2$ 的情形：给定 $|\lambda_1| \geq |\lambda_2|$，$s_1 \geq s_2 \geq 0$，$|\lambda_1| \leq s_1$，$|\lambda_1\lambda_2| = s_1 s_2$，构造 $2 \times 2$ 矩阵。这可以通过直接参数化 $A = \begin{pmatrix}a & b\\c & d\end{pmatrix}$ 并求解方程组来实现。

    2. 对一般 $n$，用归纳法。关键引理：给定满足条件的 $(\lambda_i, s_i)$，可以选择 $A$ 的第一行，使得删去第一行第一列后的 $(n-1) \times (n-1)$ 子问题仍满足条件。$\blacksquare$

!!! theorem "定理 61.15 (Mirsky 定理, 1958)"
    设 $\lambda_1, \ldots, \lambda_n \in \mathbb{C}$（给定特征值）和 $s_1 \geq \cdots \geq s_n \geq 0$（给定奇异值）。存在 $n \times n$ 复矩阵 $A$ 使得 $\sigma(A) = \{\lambda_i\}$ 且 $A$ 的奇异值为 $\{s_i\}$，当且仅当以下条件**全部**成立：

    (a) $\prod_{i=1}^k |\lambda_{\pi(i)}| \leq \prod_{i=1}^k s_i$，对所有 $k = 1, \ldots, n-1$ 和特征值的所有排列 $\pi$（使得 $|\lambda_{\pi(1)}| \geq \cdots \geq |\lambda_{\pi(n)}|$，即按模递减排列后的乘积不等式）。

    (b) $\prod_{i=1}^n |\lambda_i| = \prod_{i=1}^n s_i$。

    (c) 非实特征值成共轭对。

    实际上条件 (a)(b) 等价于 Weyl-Horn 条件（定理 61.14），条件 (c) 在考虑实矩阵时额外需要。对复矩阵，条件 (a)(b) 已经充分。

!!! example "例 61.9"
    **特征值-奇异值兼容性检验。** 给定 $\lambda_1 = 3, \lambda_2 = -1$，$s_1 = 3, s_2 = 1$。

    - $|\lambda_1| = 3 \leq s_1 = 3$：成立。
    - $|\lambda_1||\lambda_2| = 3 = s_1 s_2 = 3$：成立（等号）。

    构造：需要 $2 \times 2$ 实矩阵 $A$，特征值 $3, -1$，奇异值 $3, 1$。

    取 $A = \begin{pmatrix}a & b\\c & d\end{pmatrix}$。$a + d = 3 + (-1) = 2$，$ad - bc = 3 \cdot (-1) = -3$。
    $A^\top A$ 的特征值为 $s_1^2 = 9, s_2^2 = 1$，所以 $\mathrm{tr}(A^\top A) = 10$，$\det(A^\top A) = 9$。
    $a^2 + b^2 + c^2 + d^2 = 10$，$(ad-bc)^2 = 9$。

    取 $A = \begin{pmatrix}2 & 1\\1 & 0\end{pmatrix}$：特征值 $1 \pm \sqrt{2}$——不对。

    取对角化方法：$A = PDP^{-1}$，$D = \mathrm{diag}(3, -1)$。取 $P = \begin{pmatrix}\cos\theta & -\sin\theta\\\sin\theta & \cos\theta\end{pmatrix}$，则 $A = P\begin{pmatrix}3&0\\0&-1\end{pmatrix}P^\top = \begin{pmatrix}3\cos^2\theta - \sin^2\theta & 4\sin\theta\cos\theta\\4\sin\theta\cos\theta & 3\sin^2\theta - \cos^2\theta\end{pmatrix}$。

    $A^\top A$ 的特征值（即奇异值平方）为 $9, 1$ 当且仅当 $\mathrm{tr}(A^\top A) = 10$ 且 $\det(A^\top A) = 9$。由于 $A$ 是对称的（$P$ 正交且 $D$ 对称），$A^\top A = A^2 = P D^2 P^\top$，奇异值为 $|3| = 3$ 和 $|-1| = 1$。因此取 $\theta = 0$：$A = \mathrm{diag}(3, -1)$，奇异值 $3, 1$。成功。

---

## 61.9 Friedland 充分条件与 NIEP 的进展

<div class="context-flow" markdown>

**核心问题**：在 Suleimanova 条件之外，还有哪些 NIEP 可实现性的充分条件？

</div>

!!! theorem "定理 61.16 (Friedland 充分条件, 1978)"
    设 $\Lambda = \{\lambda_1, \ldots, \lambda_n\} \subset \mathbb{C}$，$\lambda_1 \geq |\lambda_i|$（Perron 根条件），非实特征值成共轭对。如果

    $$\lambda_1 \geq \sum_{i=2}^n |\lambda_i|,$$

    则 $\Lambda$ 是可实现的（存在 $n \times n$ 非负矩阵具有谱 $\Lambda$）。

??? proof "证明"
    构造伴随矩阵。设特征多项式为 $p(\lambda) = \prod_{i=1}^n (\lambda - \lambda_i) = \lambda^n - c_1\lambda^{n-1} - \cdots - c_n$。

    $c_1 = \sum \lambda_i = s_1$。由 Newton 恒等式和条件 $\lambda_1 \geq \sum_{i \geq 2}|\lambda_i|$，可以证明所有 $c_i \geq 0$。

    关键不等式：$c_k = e_k(\lambda_1, \ldots, \lambda_n)$（初等对称函数的符号调整）。由 $\lambda_1$ 的主导性，$\lambda_1$ 的贡献控制了其余特征值的影响，保证非负性。

    具体地，伴随矩阵

    $$C = \begin{pmatrix}0 & 0 & \cdots & 0 & c_n\\1 & 0 & \cdots & 0 & c_{n-1}\\0 & 1 & \cdots & 0 & c_{n-2}\\\vdots & & \ddots & & \vdots\\0 & 0 & \cdots & 1 & c_1\end{pmatrix}$$

    的所有元素在条件 $\lambda_1 \geq \sum_{i \geq 2}|\lambda_i|$ 下非负。$\blacksquare$

!!! theorem "定理 61.17 (Kellogg-Stephens 充分条件)"
    设 $\Lambda = \{\lambda_1, \ldots, \lambda_n\} \subset \mathbb{C}$，$\lambda_1 \geq 0$ 为 Perron 根。如果

    $$s_k(\Lambda) = \sum_{i=1}^n \lambda_i^k \geq 0, \quad \forall k = 1, 2, \ldots, 2n-2,$$

    且 $\lambda_1 \geq |\lambda_i|$，非实特征值成共轭对，则 $\Lambda$ 是可实现的。

    这一条件比 Suleimanova 条件更一般（允许复特征值），但仍不是充要的。

!!! theorem "定理 61.18 (Laffey-Smigoc 充分条件, 2006)"
    设 $\Lambda = \{\lambda_1, \lambda_2, \ldots, \lambda_n\}$ 满足 $\lambda_1 \geq |\lambda_i|$ 和共轭对称。如果存在**可实现的**分划 $\Lambda = \Lambda_1 \cup \Lambda_2$（$\Lambda_1 \cap \Lambda_2 = \emptyset$），即 $\Lambda_1$ 和 $\Lambda_2$ 分别可被非负矩阵 $A_1, A_2$ 实现，且 Perron 根满足适当的兼容条件，则 $\Lambda$ 也是可实现的。

    此"可实现分划"方法为 NIEP 提供了递归构造工具。

!!! example "例 61.10"
    **Friedland 条件的应用。** $\Lambda = \{6, -1, -2, -3\}$。

    检验 Friedland 条件：$\lambda_1 = 6 \geq |-1| + |-2| + |-3| = 6$。等号成立。

    Friedland 条件要求 $\lambda_1 \geq \sum_{i \geq 2}|\lambda_i|$（含等号），因此适用。

    构造特征多项式：$p(\lambda) = (\lambda-6)(\lambda+1)(\lambda+2)(\lambda+3) = (\lambda-6)(\lambda^3 + 6\lambda^2 + 11\lambda + 6)$
    $= \lambda^4 + 6\lambda^3 + 11\lambda^2 + 6\lambda - 6\lambda^3 - 36\lambda^2 - 66\lambda - 36$
    $= \lambda^4 - 25\lambda^2 - 60\lambda - 36$。

    即 $\lambda^4 - 0\cdot\lambda^3 - 25\lambda^2 - 60\lambda - 36$。$c_1 = 0, c_2 = 25, c_3 = 60, c_4 = 36$，均 $\geq 0$。

    伴随矩阵：

    $$C = \begin{pmatrix}0 & 0 & 0 & 36\\1 & 0 & 0 & 60\\0 & 1 & 0 & 25\\0 & 0 & 1 & 0\end{pmatrix}.$$

    所有元素非负，特征值为 $\{6, -1, -2, -3\}$。

---

## 61.10 逆奇异值问题

<div class="context-flow" markdown>

**核心问题**：给定奇异值和结构约束，何时能构造相应的矩阵？

</div>

!!! definition "定义 61.13 (逆奇异值问题)"
    **逆奇异值问题**（Inverse Singular Value Problem, ISVP）：给定 $s_1 \geq s_2 \geq \cdots \geq s_n \geq 0$ 和矩阵结构类 $\mathcal{S}$，构造 $A \in \mathcal{S}$ 使得 $A$ 的奇异值为 $\{s_i\}$。

!!! theorem "定理 61.19 (无约束 ISVP)"
    给定任意 $s_1 \geq \cdots \geq s_n \geq 0$，总存在 $n \times n$ 实矩阵 $A$ 具有这些奇异值。

??? proof "证明"
    取 $A = \mathrm{diag}(s_1, \ldots, s_n)$。$\blacksquare$

    更一般地，$A = U \,\mathrm{diag}(s_1, \ldots, s_n)\, V^\top$ 对任意正交矩阵 $U, V$ 均可。

!!! theorem "定理 61.20 (带规定对角元的 ISVP)"
    给定 $s_1 \geq \cdots \geq s_n \geq 0$ 和 $d_1, \ldots, d_n \in \mathbb{R}$，存在实矩阵 $A$ 使得 $A$ 的奇异值为 $\{s_i\}$ 且对角元素为 $\{d_i\}$，当且仅当

    $$(d_1, \ldots, d_n) \in \mathrm{conv}\{(\pm s_{\sigma(1)}, \ldots, \pm s_{\sigma(n)}) : \sigma \in S_n\}.$$

    即对角元素向量属于奇异值的所有带符号置换形成的凸包。

!!! definition "定义 61.14 (Jacobi 型 ISVP)"
    **Jacobi 逆奇异值问题**：给定奇异值 $s_1 \geq \cdots \geq s_n > 0$，构造具有这些奇异值的 **Jacobi 矩阵**（对称正定三对角矩阵）。

    此问题等价于：给定 $s_i^2$（即 $A^2 = A^\top A = A^2$ 的特征值），构造对称正定三对角 $A$，使得 $A^2$ 的特征值为 $\{s_i^2\}$。

!!! theorem "定理 61.21 (Jacobi ISVP 的可解性)"
    给定 $s_1 \geq \cdots \geq s_n > 0$，如果 $s_i$ 两两不同，则存在（不一定唯一的）对称正定三对角矩阵 $J$，使得 $J$ 的奇异值为 $\{s_i\}$。

    构造方法类似于 Jacobi IEP：先将 $\{s_i^2\}$ 作为 $J^2$ 的目标特征值，利用 $J^2$ 的五对角结构和 Lanczos 过程来恢复 $J$。

!!! example "例 61.11"
    **简单的 ISVP 构造。** 给定奇异值 $s_1 = 5, s_2 = 3$，且要求 $A$ 为 $2 \times 2$ 上三角矩阵。

    设 $A = \begin{pmatrix}a & b\\0 & d\end{pmatrix}$。$A^\top A = \begin{pmatrix}a^2 & ab\\ab & b^2+d^2\end{pmatrix}$。

    条件：$\mathrm{tr}(A^\top A) = a^2 + b^2 + d^2 = 25 + 9 = 34$，$\det(A^\top A) = a^2 d^2 = 225$。

    由 $ad = 15$（取正值），$d = 15/a$。$a^2 + b^2 + 225/a^2 = 34$。

    取 $a = 5$：$25 + b^2 + 9 = 34$，$b^2 = 0$，$b = 0$。$A = \mathrm{diag}(5, 3)$。

    取 $a = 3$：$9 + b^2 + 25 = 34$，$b^2 = 0$，$b = 0$。$A = \mathrm{diag}(3, 5)$。

    取 $a = 4$：$16 + b^2 + 225/16 = 34$，$b^2 = 34 - 16 - 14.0625 = 3.9375$，$b = \sqrt{63/16} = 3\sqrt{7}/4$。

    $A = \begin{pmatrix}4 & 3\sqrt{7}/4\\0 & 15/4\end{pmatrix}$，奇异值为 $5, 3$。

---

## 61.11 开放问题

<div class="context-flow" markdown>

**核心问题**：逆特征值理论中有哪些重要的未解决问题？

</div>

!!! definition "定义 61.15 (主要开放问题列表)"
    以下是逆特征值理论中最重要的开放问题：

    **(1) 一般 NIEP**：对 $n \geq 5$，给出非负矩阵可实现谱的完整充要条件。这是该领域最著名的开放问题，自 1949 年 Suleimanova 开始研究以来已超过 75 年。

    **(2) 实非负 IEP (RNIEP)**：仅考虑实特征值的 NIEP。$n \geq 5$ 的完整刻画未知。

    **(3) 对称非负 IEP (SNIEP)**：对称非负矩阵的 IEP。SNIEP 与 RNIEP 在 $n \leq 4$ 时等价，但在 $n = 5$ 时不等价（Johnson 等人的例子）。

    **(4) Toeplitz IEP**：一般 $n$ 的实对称 Toeplitz 矩阵 IEP 的完整存在性条件。

    **(5) Perfect-Mirsky 猜想**：完整刻画 $n \times n$ 双随机矩阵可能的谱的区域 $\Pi_n$。

    **(6) 结构化 ISVP**：对一般带宽约束的逆奇异值问题的完整可解性条件。

!!! theorem "定理 61.22 (已知的部分结果)"
    **(a) Boyle-Handelman 定理 (1991)**：$\Lambda = \{\lambda_1, \ldots, \lambda_n\}$ 是某个非负整数矩阵的非零特征值的充要条件是：$\lambda_1$ 是 Perron 根，且对所有 $k \geq 1$，$s_k \geq 0$，且 $\Lambda$ 关于复共轭封闭。但这不直接解决 NIEP，因为矩阵大小不固定。

    **(b) Johnson-Loewy-London 猜想**：RNIEP 和 SNIEP 在所有 $n$ 上等价。这在 $n \leq 4$ 时已证明，$n = 5$ 时已被 Johnson 等人否定。

    **(c) Weyl-Horn 定理的完备性**：Weyl-Horn 条件完全刻画了特征值-奇异值的兼容关系（对无结构约束的复矩阵），但带结构约束的联合逆问题（同时规定特征值、奇异值和矩阵结构）仍大量未解。

!!! example "例 61.12"
    **SNIEP 与 RNIEP 不等价的例子 ($n = 5$)。** Johnson, Loewy 和 London 证明存在谱 $\Lambda = \{5 + 5\epsilon, 5 - 5\epsilon, -2 + \epsilon, -2 + \epsilon, -6 - 2\epsilon\}$（适当选择 $\epsilon$），可以被非负矩阵实现但不能被对称非负矩阵实现。

    这说明对称性约束与非负约束之间的微妙差异——非负性允许的谱范围严格大于对称非负性允许的谱范围（在 $n \geq 5$ 时）。

---

**本章要点总结：**

1. 逆特征值问题从谱数据出发构造具有特定结构的矩阵，涉及存在性、唯一性和构造性三个层面。
2. 无约束对称 IEP 总有解（$A = Q\Lambda Q^\top$）；带对角约束时由 Schur-Horn 优超条件刻画。
3. NIEP 是该领域最著名的开放问题：$n \leq 4$ 已完全解决，$n \geq 5$ 至今未解。
4. Jacobi IEP 由两组交错谱唯一确定，de Boor-Golub 算法通过正交多项式给出构造。
5. Toeplitz IEP 和随机矩阵 IEP 的一般完整存在性条件也是开放问题。
6. 数值方法包括提升-投影、Newton 迭代和同伦延拓。
7. Weyl-Horn 定理完全刻画了特征值与奇异值的兼容关系；Mirsky 定理给出联合存在性的充要条件。
8. Friedland 条件和 Laffey-Smigoc 分划方法为 NIEP 提供了 Suleimanova 之外的充分条件。
9. 逆奇异值问题是逆特征值问题的自然伴侣，带结构约束时同样困难。
10. 逆特征值问题在结构工程（弹簧系统设计）、控制理论（极点配置）和分子光谱学中有直接应用。

## 练习题

1. **[概念] 什么是“逆特征值问题”（IEP）？它与普通的特征值问题在逻辑方向上有何不同？**
   ??? success "参考答案"
       普通特征值问题是“给定矩阵，求特征值”；逆特征值问题是“给定特征值（及其他谱数据），构造满足特定结构约束（如对称、非负、三对角等）的矩阵”。IEP 本质上是矩阵论中的参数重构或系统设计问题。

2. **[基础] 构造一个 $2 \times 2$ 实对称矩阵，使其特征值为 $5$ 和 $1$。**
   ??? success "参考答案"
       最简单的解是对角矩阵 $\operatorname{diag}(5, 1)$。
       另一个非平凡解：取 $\theta = 45^\circ$，进行旋转变换：
       $A = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} \begin{pmatrix} 5 & 0 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix} = \begin{pmatrix} 3 & 2 \\ 2 & 3 \end{pmatrix}$。

3. **[Schur-Horn] 判定是否存在特征值为 $\{4, 2\}$ 且对角元为 $\{3, 3\}$ 的实对称矩阵。**
   ??? success "参考答案"
       检查优超条件：$3 \le 4$ 且 $3+3 = 4+2 = 6$。条件满足，根据 Schur-Horn 定理，这样的矩阵确实存在。例如 $\begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$。

4. **[NIEP] 为什么特征值集合 $\{3, 1, 1\}$ 不能由任何非负矩阵实现？**
   ??? success "参考答案"
       因为它不满足 Perron 条件。非负矩阵的最大模特征值（Perron 根）必须是唯一的（若不可约）或至少在模长并列时具有特定结构。更直接地，根据轨迹条件 $s_k \ge 0$，虽然 $s_1=5 > 0$，但非负矩阵通常需要包含负实部或复特征值来“抵消”对角线上的正值，纯正特征值集合通常很难由非负非对角阵实现（除非是对角阵）。

5. **[Jacobi] 给出两个交错的谱 $\lambda = \{1, 3, 5\}$ 和 $\mu = \{2, 4\}$。它们能唯一确定一个 Jacobi 矩阵吗？**
   ??? success "参考答案"
       可以。根据 Hochstadt 唯一性定理，给定严格交错的 $n$ 阶谱及其 $n-1$ 阶主子阵谱，存在唯一的 Jacobi 矩阵。这对应于物理上通过两组频率重构质量-弹簧链。

6. **[计算] 设 $A$ 是非负矩阵，特征值为 $\{2, -1, -1\}$。写出其对应的伴随矩阵形式。**
   ??? success "参考答案"
       特征多项式为 $p(\lambda) = (\lambda-2)(\lambda+1)^2 = \lambda^3 - 3\lambda - 2$。
       对应的伴随矩阵为 $\begin{pmatrix} 0 & 0 & 2 \\ 1 & 0 & 3 \\ 0 & 1 & 0 \end{pmatrix}$。
       所有元素非负，验证成功。

7. **[Weyl-Horn] 若矩阵的特征值为 $\{2, i, -i\}$，其奇异值最小可能是多少？**
   ??? success "参考答案"
       特征值模长为 $\{2, 1, 1\}$。根据 Weyl 乘积不等式，$\prod s_i = \prod |\lambda_i| = 2$。
       由于 $s_i \ge |\lambda_i|$（按序排列），最小的奇异值集合就是其模长本身：$\{2, 1, 1\}$。

8. **[应用] 在结构动力学中，逆特征值问题通常用来解决什么实际问题？**
   ??? success "参考答案"
       通常用于“模型修正”或“振动设计”。例如，通过已知的实验模态频率，反推系统各部件的刚度系数或质量分布，从而建立准确的数学模型。

9. **[数值] 简述“提升-投影”（Lift-and-Project）方法求解 IEP 的基本步骤。**
   ??? success "参考答案"
       1. 提升：在谱空间中寻找具有目标特征值的矩阵（如 $Q \Lambda Q^T$）。
       2. 投影：将该矩阵投影到满足结构约束的集合 $\mathcal{S}$ 中（如强制非对角线某些元素为 0）。
       3. 迭代：交替进行这两步，直到矩阵既满足谱要求又满足结构要求。

10. **[爱因斯坦思考题] 爱因斯坦的方程 $G_{\mu\nu} = 8\pi T_{\mu\nu}$ 将时空的几何结构（左端）与物质能量分布（右端）联系起来。如果将物质分布看作矩阵的结构约束，而将时空的曲率模态看作特征值，广义相对论是否可以被视为一个关于宇宙的“逆特征值问题”？**
    ??? success "参考答案"
        这是一个深刻的类比。在广义相对论中，我们观测到光线的偏移或引力波的频率（特征值数据），然后去推断宇宙中质量和能量的分布（矩阵的条目）。宇宙就像一个巨大的、具有特定对称性约束的算子，其物理响应（谱）揭示了它内在的构造。逆特征值问题的复杂性——从不唯一的解到存在性的严苛条件——恰恰反映了物理学中从现象推导本质的艰辛。这告诉我们，宇宙的“结构”不仅决定了它的“节奏”，而我们作为观察者，正是在通过捕捉这些“节奏”来重构那不可见的、宏大的“代数结构”。

## 本章小结

本章探讨了从谱数据反向重构矩阵结构的深层理论——逆特征值问题（IEP），主要内容包括：

1. **问题的层次性**：确立了 IEP 在存在性、唯一性和数值构造三个维度的研究范式，区分了对称、非负、带状等多种结构约束。
2. **SIEP 与 NIEP**：详细论述了对称逆问题的 Schur-Horn 判据，并揭示了非负逆问题（NIEP）作为长期未解难题的复杂性与魅力。
3. **Jacobi 与结构化 IEP**：展示了 Jacobi 矩阵如何通过特征值交错属性实现参数的唯一重构，及其在物理振动系统中的直接映射。
4. **特征值与奇异值的纽带**：利用 Weyl-Horn 不等式建立了矩阵两组最核心不变量之间的代数约束关系。
5. **现代算法与应用**：介绍了提升-投影、同伦延拓等数值求解技术，并探讨了 IEP 在结构工程、控制理论及分子光谱分析中的实战价值。

