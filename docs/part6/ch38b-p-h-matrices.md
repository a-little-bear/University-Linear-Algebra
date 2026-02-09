# 第 38B 章 P-矩阵、H-矩阵与相关矩阵类

<div class="context-flow" markdown>

**前置**：非负矩阵与 Perron-Frobenius 理论(Ch17) · Z-矩阵与 M-矩阵(Ch38A) · 矩阵分析(Ch14) · 特征值理论(Ch6) · 线性互补问题基础

**本章脉络**：P-矩阵定义与等价刻画 $\to$ LCP 唯一可解性 $\to$ P$_0$-矩阵 $\to$ 比较矩阵与 H-矩阵 $\to$ H-矩阵等价条件 $\to$ 对角占优层次 $\to$ N-矩阵 $\to$ 半正矩阵 $\to$ Ostrowski-Reich 定理 $\to$ 广义 M-矩阵（复情形） $\to$ 全正矩阵联系 $\to$ 迭代法与数值稳定性应用

**延伸**：P-矩阵在数学规划（LCP 理论、变分不等式）、博弈论（纳什均衡存在性）中有核心地位；H-矩阵将 M-矩阵理论推广到复矩阵和一般符号模式，在迭代法的收敛性分析、预处理技术、区间数学中广泛应用

</div>

上一章系统发展了 Z-矩阵与 M-矩阵的理论——这类矩阵要求非对角元素全部非正。本章转向更广泛的矩阵类。P-矩阵要求所有主子式为正，不限制符号模式；H-矩阵通过比较矩阵将 M-矩阵的思想推广到一般矩阵（包括复矩阵）。这些矩阵类在线性互补问题、迭代法收敛性分析和数值稳定性中有深刻的应用。

本章还介绍 P$_0$-矩阵、N-矩阵、半正矩阵等相关矩阵类，构建一个完整的特殊矩阵类谱系。

---

## 38B.1 P-矩阵

!!! definition "定义 38B.1 (P-矩阵)"
    矩阵 $A \in M_n(\mathbb{R})$ 称为 **P-矩阵**，若其所有主子式为正：
    $$\det(A[\alpha, \alpha]) > 0, \quad \forall\, \alpha \subseteq \{1, 2, \ldots, n\}, \quad \alpha \ne \emptyset,$$
    其中 $A[\alpha, \alpha]$ 表示行列指标集均为 $\alpha$ 的主子矩阵。

    P-矩阵全体记为 $\mathcal{P}_n$。

!!! theorem "定理 38B.1 (P-矩阵的等价刻画)"
    设 $A \in M_n(\mathbb{R})$。以下条件等价：

    **(P1)** $A$ 是 P-矩阵（所有主子式为正）。

    **(P2)** 对每个非零向量 $x \in \mathbb{R}^n$，存在指标 $i$ 使得 $x_i(Ax)_i > 0$。

    **(P3)** 对每个对角矩阵 $D$（对角元素为 $\pm 1$），$DA$ 非奇异。

    **(P4)** $A$ 不将 $\mathbb{R}^n$ 中任何非平凡的正交象限映射到其对立象限：不存在非零 $x$ 使得对所有 $i$，$x_i(Ax)_i \le 0$。

    **(P5)** 线性互补问题 $\operatorname{LCP}(q, A)$：对所有 $q \in \mathbb{R}^n$，找 $w, z \ge 0$ 使得 $w = q + Az$，$w^T z = 0$，有唯一解。

    **(P6)** $A$ 的所有实特征值为正，且 $A$ 的每个主子矩阵也有此性质。

    **(P7)** 存在正对角矩阵 $D_1, D_2$ 使得 $D_1 A D_2 + D_2 A^T D_1$ 正定。

??? proof "证明"
    **(P1) $\Rightarrow$ (P2)**：反证法。设存在 $x \ne 0$ 使得对所有 $i$，$x_i(Ax)_i \le 0$。设 $\alpha = \{i : x_i \ne 0\}$（非空）。对 $i \in \alpha$，$x_i(Ax)_i \le 0$。

    注意 $(Ax)_i = \sum_{j \in \alpha} a_{ij}x_j$（因为 $x_j = 0$ 对 $j \notin \alpha$），即 $(Ax)_\alpha = A[\alpha,\alpha]x_\alpha$。

    定义对角矩阵 $S = \operatorname{diag}(\operatorname{sgn}(x_i))_{i \in \alpha}$（$\operatorname{sgn}(x_i) = \pm 1$），令 $y = |x_\alpha|$（逐元取绝对值），则 $x_\alpha = Sy$，$y > 0$。

    $x_i(Ax)_i \le 0$ 对 $i \in \alpha$ 意味着 $y^T S \cdot A[\alpha,\alpha] \cdot Sy \le 0$，即 $y^T (SA[\alpha,\alpha]S) y \le 0$（$y > 0$）。

    设 $B = SA[\alpha,\alpha]S$。$\det(B) = (\det S)^2 \det(A[\alpha,\alpha]) = \det(A[\alpha,\alpha]) > 0$（由 P1）。但 $y^T B y \le 0$（$y > 0$）意味着 $B$ 不是正定的。

    然而 $B$ 的主子式与 $A[\alpha,\alpha]$ 的主子式在符号上可能不同（$S$ 的缩放改变了主子式），所以需要更精细的分析。

    实际上，$B[\beta,\beta] = S[\beta,\beta] A[\alpha,\alpha][\beta,\beta] S[\beta,\beta]$（对 $\beta \subseteq \alpha$），$\det(B[\beta,\beta]) = \det(A[\alpha,\alpha][\beta,\beta]) > 0$。所以 $B$ 也是 P-矩阵。

    $B$ 是 P-矩阵意味着对任何非零 $y$，存在 $i$ 使 $y_i(By)_i > 0$（由 P2 对 P-矩阵 $B$——这是循环论证）。

    采用更直接的论证：由归纳法对 $|\alpha|$ 归纳。当 $|\alpha| = 1$ 时，$x_i(Ax)_i = x_i a_{ii} x_i = a_{ii}x_i^2$。若 $a_{ii} > 0$（P1 蕴含），则 $x_i(Ax)_i > 0$，矛盾。

    归纳步：设对 $|\alpha| \le k-1$ 成立。取 $|\alpha| = k$。由 $x_i(Ax)_i \le 0$ 对所有 $i \in \alpha$，求和得 $x_\alpha^T A[\alpha,\alpha] x_\alpha \le 0$。由 P1，$A[\alpha,\alpha]$ 的所有主子式为正，特别是 $\det(A[\alpha,\alpha]) > 0$。由 $A[\alpha,\alpha]$ 的特征值分析（所有主子式为正意味着特征值的某些性质），可以得到矛盾。

    完整的证明利用了 Fiedler-Pták 的度论证（degree argument）或 Gale-Nikaido 定理。

    **(P2) $\Rightarrow$ (P1)**：对任意非空 $\alpha$，取 $A[\alpha,\alpha]$。设 $\lambda$ 为 $A[\alpha,\alpha]$ 的特征值，$A[\alpha,\alpha]v = \lambda v$（$v \ne 0$）。定义 $x \in \mathbb{R}^n$ 使 $x_\alpha = v$, $x_j = 0$（$j \notin \alpha$）。则 $(Ax)_i = (A[\alpha,\alpha]v)_i = \lambda v_i$ 对 $i \in \alpha$。

    由 P2，存在 $i \in \alpha$ 使 $x_i(Ax)_i = v_i \cdot \lambda v_i = \lambda|v_i|^2 > 0$（取 $v$ 为实特征向量时）。故实特征值 $\lambda > 0$。

    $A[\alpha,\alpha]$ 的所有实特征值为正。复特征值成共轭对，对行列式贡献 $|\lambda|^2 > 0$。故 $\det(A[\alpha,\alpha]) > 0$。

    **(P1) $\Rightarrow$ (P3)**：设 $D = \operatorname{diag}(d_1,\ldots,d_n)$，$d_i = \pm 1$。$(DA)[\alpha,\alpha] = D[\alpha,\alpha]A[\alpha,\alpha]$，$\det((DA)[\alpha,\alpha]) = \prod_{i \in \alpha}d_i \cdot \det(A[\alpha,\alpha])$。但这不保证正性。

    更正确地：P3 等价于说 $A$ 的每个主子矩阵非奇异。设 $DA$ 奇异，$DAx = 0$，$x \ne 0$。则 $Ax = D^{-1} \cdot 0 = 0$。但 $\det(A) > 0$（P1），故 $A$ 非奇异，$x = 0$，矛盾。

    等等，这只说明了 $A$ 本身非奇异。P3 的精确含义是对**每个**符号矩阵 $D$，$DA$ 非奇异。这需要更深入的分析。

    实际上 P3 的等价性来自于：$DA$ 奇异当且仅当存在非零 $x$ 使 $DAx = 0$，即 $d_i(Ax)_i = 0$ 对所有 $i$。若选择 $d_i = \operatorname{sgn}(x_i(Ax)_i)$（当 $x_i(Ax)_i \ne 0$），则...更标准的证明是通过 (P2) 作为中介。

    **(P1) $\Leftrightarrow$ (P5)** 的证明涉及互补枢轴（complementary pivoting）理论，见 Cottle-Pang-Stone (1992)。关键思想：LCP(q,A) 的解对应于 $2^n$ 个互补锥的交集。P-矩阵保证每个互补锥恰好覆盖 $\mathbb{R}^n$ 的一个部分，不重叠不遗漏。

    **(P1) $\Rightarrow$ (P6)**：P-矩阵的每个主子矩阵也是 P-矩阵（主子式的正性对主子矩阵遗传），每个主子矩阵的实特征值为正（由 P2 的论证）。

!!! example "例 38B.1"
    矩阵 $A = \begin{pmatrix} 1 & -2 \\ 1 & 1 \end{pmatrix}$ 是否为 P-矩阵？

    主子式：$a_{11} = 1 > 0$，$a_{22} = 1 > 0$，$\det(A) = 1 + 2 = 3 > 0$。所有主子式为正，$A$ 是 P-矩阵。

    注意 $A$ **不是** Z-矩阵（$a_{21} = 1 > 0$），更不是 M-矩阵。P-矩阵不限制符号模式。

!!! example "例 38B.2"
    矩阵 $B = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$ 不是 P-矩阵。

    $\det(B) = 1 - 4 = -3 < 0$。

    验证 P2 不成立：取 $x = (1, -1)^T$，$Bx = (-1, 1)^T$，$x_1(Bx)_1 = -1 < 0$，$x_2(Bx)_2 = -1 < 0$。

!!! theorem "定理 38B.2 (P-矩阵与 M-矩阵的关系)"
    (a) 每个非奇异 M-矩阵都是 P-矩阵。

    (b) 反之不然：P-矩阵不要求非对角元素非正。

    (c) P-矩阵 $\cap$ Z-矩阵 $=$ 非奇异 M-矩阵。

??? proof "证明"
    (a) 非奇异 M-矩阵的条件 M5 直接给出所有主子式为正。

    (b) 例 38B.1 给出了一个是 P-矩阵但不是 Z-矩阵（从而不是 M-矩阵）的例子。

    (c) 设 $A$ 是 P-矩阵且是 Z-矩阵。所有主子式为正（P-矩阵），特别是 $\det(A) > 0$，故 $A$ 非奇异。又 $A$ 是 Z-矩阵，$A = sI - B$（$B \ge 0$）。所有主子式为正意味着条件 M5，故 $A$ 是非奇异 M-矩阵。

    反过来，非奇异 M-矩阵满足 M5（所有主子式正）且是 Z-矩阵。$\blacksquare$

---

## 38B.2 线性互补问题（LCP）

!!! definition "定义 38B.2 (线性互补问题)"
    给定 $A \in M_n(\mathbb{R})$ 和 $q \in \mathbb{R}^n$，**线性互补问题** $\operatorname{LCP}(q, A)$ 是求 $w, z \in \mathbb{R}^n$ 使得
    $$w = q + Az, \quad w \ge 0, \quad z \ge 0, \quad w^T z = 0.$$
    互补条件 $w^T z = 0$ 意味着对每个 $i$，$w_i = 0$ 或 $z_i = 0$（或两者均为零）。

!!! theorem "定理 38B.3 (P-矩阵与 LCP 的唯一可解性)"
    $A$ 是 P-矩阵当且仅当对每个 $q \in \mathbb{R}^n$，$\operatorname{LCP}(q, A)$ 有唯一解。

??? proof "证明"
    **充分性** ($A$ 是 P-矩阵 $\Rightarrow$ LCP 唯一可解)：

    **存在性**：由度论证（degree-theoretic argument）。定义映射 $F: \mathbb{R}^n \to \mathbb{R}^n$，$F_i(z) = \min(z_i, (q+Az)_i)$。$F$ 的零点即 LCP 的解。P-矩阵的性质保证 $F$ 是本性映射（proper map），由拓扑度理论，$F$ 有零点。

    **唯一性**：设 $(w^1, z^1)$ 和 $(w^2, z^2)$ 是两个解。令 $\bar{z} = z^1 - z^2$。对每个 $i$，由互补条件分析四种情形：

    - 若 $z^1_i > 0, z^2_i > 0$，则 $w^1_i = w^2_i = 0$，$(q + Az^1)_i = (q + Az^2)_i = 0$，$(A\bar{z})_i = 0$。
    - 若 $z^1_i > 0, z^2_i = 0$，则 $w^1_i = 0, w^2_i \ge 0$，$\bar{z}_i = z^1_i > 0$，$(A\bar{z})_i = -w^2_i \le 0$，故 $\bar{z}_i(A\bar{z})_i \le 0$。
    - 若 $z^1_i = 0, z^2_i > 0$，类似得 $\bar{z}_i(A\bar{z})_i \le 0$。
    - 若 $z^1_i = z^2_i = 0$，则 $\bar{z}_i = 0$，$\bar{z}_i(A\bar{z})_i = 0$。

    综上，对所有 $i$，$\bar{z}_i(A\bar{z})_i \le 0$。由 P-矩阵条件 P2，$\bar{z} = 0$，即 $z^1 = z^2$。

    **必要性**（LCP 唯一可解 $\Rightarrow$ P-矩阵）：反证。若 $A$ 不是 P-矩阵，存在 $\alpha \ne \emptyset$ 使 $\det(A[\alpha,\alpha]) \le 0$。可以构造 $q$ 使得 LCP 有多个解或无解。具体构造见 Cottle-Pang-Stone (1992)。$\blacksquare$

!!! example "例 38B.3 (LCP 的几何解释)"
    对 $n = 2$，$A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$（M-矩阵，从而 P-矩阵），$q = \begin{pmatrix} -1 \\ -1 \end{pmatrix}$。

    $\operatorname{LCP}(q, A)$：$w = q + Az$，$w, z \ge 0$，$w^T z = 0$。

    由互补性，$\{1,2\}$ 分为 $\alpha = \{i: z_i > 0\}$（$w_i = 0$）和 $\beta = \{i: w_i > 0\}$（$z_i = 0$）。

    尝试 $\alpha = \{1,2\}$：$w = 0$，$q + Az = 0$，$z = -A^{-1}q = \frac{1}{3}\begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \ge 0$。解为 $z = (1,1)^T$，$w = (0,0)^T$。

---

## 38B.3 P$_0$-矩阵

!!! definition "定义 38B.3 (P$_0$-矩阵)"
    矩阵 $A \in M_n(\mathbb{R})$ 称为 **P$_0$-矩阵**，若其所有主子式非负：
    $$\det(A[\alpha, \alpha]) \ge 0, \quad \forall\, \alpha \subseteq \{1, 2, \ldots, n\}, \quad \alpha \ne \emptyset.$$

!!! theorem "定理 38B.4 (P$_0$-矩阵的性质)"
    (a) 每个 P-矩阵是 P$_0$-矩阵。

    (b) P$_0$-矩阵是 P-矩阵的闭包：$A$ 是 P$_0$-矩阵当且仅当对每个 $\epsilon > 0$，$A + \epsilon I$ 是 P-矩阵。

    (c) P$_0$-矩阵可以是奇异的。

    (d) $A$ 是 P$_0$-矩阵当且仅当对每个非零 $x$，存在 $i$ 使得 $x_i \ne 0$ 且 $x_i(Ax)_i \ge 0$。

    (e) 奇异 M-矩阵（$K_0$ 类）是 P$_0$-矩阵。

    (f) P$_0$-矩阵的主子矩阵仍为 P$_0$-矩阵。

??? proof "证明"
    (a) 显然（$> 0$ 蕴含 $\ge 0$）。

    (b) **必要性**：设 $A$ 是 P$_0$-矩阵。$A + \epsilon I$ 的主子矩阵为 $A[\alpha,\alpha] + \epsilon I_{|\alpha|}$。设 $\lambda_1, \ldots, \lambda_k$ 为 $A[\alpha,\alpha]$ 的特征值，则 $A[\alpha,\alpha] + \epsilon I$ 的特征值为 $\lambda_i + \epsilon$。$\det(A[\alpha,\alpha] + \epsilon I) = \prod(\lambda_i + \epsilon)$。P$_0$-矩阵的特征值满足 $\operatorname{Re}(\lambda_i) \ge 0$（需要论证），加上 $\epsilon > 0$ 后 $\operatorname{Re}(\lambda_i + \epsilon) > 0$。具体地，由条件 (d) 对 P$_0$-矩阵，每个主子矩阵也满足相同条件，经过分析可得其实特征值非负。故 $\det(A[\alpha,\alpha]+\epsilon I) > 0$。

    **充分性**：若对每个 $\epsilon > 0$，$A + \epsilon I$ 是 P-矩阵，则 $\det(A[\alpha,\alpha]+\epsilon I) > 0$ 对所有 $\epsilon > 0$。令 $\epsilon \to 0$，得 $\det(A[\alpha,\alpha]) \ge 0$。

    (c) $A = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$ 是 P$_0$-矩阵（所有主子式 $= 0 \ge 0$）且奇异。

    (d) 类似 P-矩阵条件 P2 的弱化。

    (e) 奇异 M-矩阵 $A = \rho(B)I - B$（$B \ge 0$）的主子矩阵（不是 $A$ 本身的）是非奇异 M-矩阵（主子式 $> 0$），$\det(A) = 0 \ge 0$。对更一般的主子矩阵，$\det(A[\alpha,\alpha]) \ge 0$（可由 $A + \epsilon I$ 是非奇异 M-矩阵取极限得到）。$\blacksquare$

!!! example "例 38B.4"
    $A = \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$ 是 P$_0$-矩阵但不是 P-矩阵。

    $a_{11} = 1 > 0$，$a_{22} = 1 > 0$，$\det(A) = 0 \ge 0$（但不 $> 0$）。

    这也是奇异 M-矩阵：$A = I - B$，$B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$，$\rho(B) = 1 = s$。

---

## 38B.4 比较矩阵与 H-矩阵

!!! definition "定义 38B.4 (比较矩阵)"
    设 $A = (a_{ij}) \in M_n(\mathbb{C})$。$A$ 的**比较矩阵**（comparison matrix）$\mathcal{M}(A) = (\mu_{ij})$ 定义为
    $$\mu_{ij} = \begin{cases} |a_{ii}| & \text{if } i = j, \\ -|a_{ij}| & \text{if } i \ne j. \end{cases}$$
    显然 $\mathcal{M}(A)$ 是 Z-矩阵。

!!! theorem "定理 38B.5 (比较矩阵的性质)"
    (a) $\mathcal{M}(A)$ 是 Z-矩阵。

    (b) $\mathcal{M}(\alpha A) = |\alpha| \mathcal{M}(A)$（$\alpha \in \mathbb{C}$）。

    (c) 若 $D$ 是正对角矩阵，则 $\mathcal{M}(DA) = D\mathcal{M}(A)$ 且 $\mathcal{M}(AD) = \mathcal{M}(A)D$。

    (d) $|\det(A)| \ge \det(\mathcal{M}(A))$（Ostrowski 行列式不等式）。

    (e) $\mathcal{M}(A)$ 的非对角元素取 $A$ 的对应元素的模的负值，故 $\mathcal{M}(A)$ 的非对角元"最坏化"了 $A$ 的非对角贡献。

??? proof "证明"
    (a) $\mu_{ij} = -|a_{ij}| \le 0$（$i \ne j$）。

    (b) $(\mathcal{M}(\alpha A))_{ii} = |\alpha a_{ii}| = |\alpha||a_{ii}|$，$(\mathcal{M}(\alpha A))_{ij} = -|\alpha a_{ij}| = -|\alpha||a_{ij}|$（$i \ne j$）。

    (c) $(DA)_{ij} = d_i a_{ij}$。$|d_i a_{ij}| = d_i |a_{ij}|$（$d_i > 0$）。$(D\mathcal{M}(A))_{ij} = d_i \mu_{ij}$。当 $i = j$：$d_i|a_{ii}|$；当 $i \ne j$：$-d_i|a_{ij}|$。与 $\mathcal{M}(DA)$ 一致。

    (d) 这是 Ostrowski 的经典结果。证明利用了行列式的多线性和三角不等式的逐步应用。对 $n = 1$，$|\det(A)| = |a_{11}| = \det(\mathcal{M}(A))$。对一般 $n$，通过 Laplace 展开和归纳可以建立不等式。完整证明见 Horn-Johnson (2013)。

!!! definition "定义 38B.5 (H-矩阵)"
    矩阵 $A \in M_n(\mathbb{C})$ 称为 **H-矩阵**，若其比较矩阵 $\mathcal{M}(A)$ 是非奇异 M-矩阵。

    H-矩阵全体记为 $\mathcal{H}_n$。

!!! theorem "定理 38B.6 (H-矩阵的等价条件)"
    设 $A \in M_n(\mathbb{C})$。以下条件等价：

    **(H1)** $A$ 是 H-矩阵（$\mathcal{M}(A)$ 是非奇异 M-矩阵）。

    **(H2)** 存在正对角矩阵 $D$ 使得 $AD$ 是严格行对角占优的。

    **(H3)** 存在正向量 $x > 0$ 使得 $\mathcal{M}(A)x > 0$，即
    $$|a_{ii}|x_i > \sum_{j \ne i} |a_{ij}|x_j, \quad \forall\, i.$$

    **(H4)** $A$ 非奇异，且 $|\det(A)| \ge \det(\mathcal{M}(A)) > 0$。

    **(H5)** 对每个非负对角矩阵 $\Delta \ge 0$，$A + \Delta$ 非奇异。

    **(H6)** $\mathcal{M}(A)^{-1} \ge 0$。

    **(H7)** 写 $\mathcal{M}(A) = sI - B$（$B \ge 0$），则 $s > \rho(B)$。

    **(H8)** 对 $A$ 的每个主子矩阵 $A[\alpha,\alpha]$，$\mathcal{M}(A[\alpha,\alpha])$ 是 M-矩阵。

??? proof "证明"
    H1 $\Leftrightarrow$ H6 $\Leftrightarrow$ H7 直接来自 M-矩阵的定义和等价条件。

    **(H1) $\Leftrightarrow$ (H2)**：

    $(\Rightarrow)$：$\mathcal{M}(A)$ 是 M-矩阵。由 M-矩阵条件 M10，存在正对角矩阵 $D$ 使 $D\mathcal{M}(A)$ 严格行对角占优。但 $D\mathcal{M}(A) = \mathcal{M}(DA)$（定理 38B.5(c)），所以 $DA$ 严格行对角占优。取 $D' = D$，$AD'$ 也可以通过右乘对角矩阵实现（$\mathcal{M}(A)D = \mathcal{M}(AD)$），故存在正对角 $D$ 使 $AD$ 严格行对角占优。

    $(\Leftarrow)$：$AD$ 严格行对角占优意味着 $\mathcal{M}(AD) = \mathcal{M}(A)D$ 严格行对角占优。$\mathcal{M}(A)D$ 是 Z-矩阵且严格行对角占优，由定理 38A.4，$\mathcal{M}(A)D$ 是 M-矩阵。$D > 0$ 可逆，$\mathcal{M}(A) = (\mathcal{M}(A)D)D^{-1}$ 也是 M-矩阵（M-矩阵右乘正对角矩阵的逆仍为 M-矩阵，因为 Z-性质和谱半径条件均保持）。

    **(H1) $\Leftrightarrow$ (H3)**：$\mathcal{M}(A)$ 是 M-矩阵 $\Leftrightarrow$ 条件 M6：存在 $x > 0$ 使 $\mathcal{M}(A)x > 0$。展开：$(\mathcal{M}(A)x)_i = |a_{ii}|x_i - \sum_{j \ne i}|a_{ij}|x_j > 0$。

    **(H1) $\Rightarrow$ (H4)**：$\mathcal{M}(A)$ 是 M-矩阵意味着 $\det(\mathcal{M}(A)) > 0$（条件 M5）。由 Ostrowski 不等式（定理 38B.5(d)），$|\det(A)| \ge \det(\mathcal{M}(A)) > 0$，故 $A$ 非奇异。

    **(H1) $\Rightarrow$ (H5)**：设 $\Delta \ge 0$ 为非负对角矩阵。$\mathcal{M}(A + \Delta) = \mathcal{M}(A) + \Delta$（非对角元素不变，对角元素 $|a_{ii} + \delta_i| \ge |a_{ii}|$——当 $\delta_i \ge 0$ 且 $a_{ii}$ 实正时等号和不等号需分析...实际上比较矩阵对一般复矩阵加非负对角的行为需要注意）。

    更直接地：$\mathcal{M}(A)$ 是 M-矩阵意味着条件 M13：对每个非负对角 $\Delta$，$\mathcal{M}(A) + \Delta$ 非奇异。由 $|\det(A + \Delta)| \ge \det(\mathcal{M}(A+\Delta)) \ge \det(\mathcal{M}(A) + \Delta) > 0$（需要 $\mathcal{M}(A+\Delta) \ge \mathcal{M}(A) + \Delta$ 的逐元不等式...这在一般情况下需要假设对角元素为正实数）。当 $A$ 的对角元素为正实数时（H-矩阵的常见情形），$\mathcal{M}(A + \Delta) = \mathcal{M}(A) + \Delta$，结论成立。

    **(H1) $\Leftrightarrow$ (H8)**：M-矩阵条件 M21 的推广。$\blacksquare$

!!! example "例 38B.5"
    $A = \begin{pmatrix} 3 & 1+i & -1 \\ -2 & 5 & 1 \\ 1 & -1 & 4 \end{pmatrix}$。

    比较矩阵：$\mathcal{M}(A) = \begin{pmatrix} 3 & -\sqrt{2} & -1 \\ -2 & 5 & -1 \\ -1 & -1 & 4 \end{pmatrix}$。

    行占优检查：行 1：$3 > \sqrt{2} + 1 \approx 2.41$；行 2：$5 > 2 + 1 = 3$；行 3：$4 > 1 + 1 = 2$。$\mathcal{M}(A)$ 严格行对角占优，故是 M-矩阵，$A$ 是 H-矩阵。

!!! example "例 38B.6"
    $A = \begin{pmatrix} 2 & -3 \\ 1 & 2 \end{pmatrix}$。$\mathcal{M}(A) = \begin{pmatrix} 2 & -3 \\ -1 & 2 \end{pmatrix}$。$\det(\mathcal{M}(A)) = 4 - 3 = 1 > 0$。

    $\mathcal{M}(A)$ 不是严格行对角占优（行 1：$2 < 3$），但仍可能是 M-矩阵。检验：$\mathcal{M}(A) = 2I - B$，$B = \begin{pmatrix} 0 & 3 \\ 1 & 0 \end{pmatrix}$，$\rho(B) = \sqrt{3} \approx 1.73 < 2$。故 $\mathcal{M}(A)$ 是 M-矩阵，$A$ 是 H-矩阵。

    验证 H3：取 $x = \mathcal{M}(A)^{-1}\mathbf{1} = \begin{pmatrix} 2 & 3 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix}/1 = \begin{pmatrix} 5 \\ 3 \end{pmatrix} > 0$。$\mathcal{M}(A)x = \begin{pmatrix} 2\cdot 5 - 3\cdot 3 \\ -1\cdot 5 + 2\cdot 3 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} > 0$。

---

## 38B.5 对角占优层次

!!! definition "定义 38B.6 (对角占优的各种定义)"
    设 $A = (a_{ij}) \in M_n(\mathbb{C})$。

    (a) $A$ 称为**严格行对角占优**的，若
    $$|a_{ii}| > \sum_{j \ne i} |a_{ij}|, \quad \forall\, i.$$

    (b) $A$ 称为**严格列对角占优**的，若 $A^T$ 严格行对角占优。

    (c) $A$ 称为**不可约对角占优**的，若 $A$ 不可约且
    $$|a_{ii}| \ge \sum_{j \ne i} |a_{ij}|, \quad \forall\, i,$$
    且至少对一个 $i$ 严格不等号成立。

    (d) $A$ 称为**广义严格对角占优**的（或 Ostrowski 意义下的对角占优），若存在正对角矩阵 $D$ 使得 $AD$ 严格行对角占优。

!!! theorem "定理 38B.7 (对角占优层次)"
    以下包含关系成立：
    $$\text{严格对角占优} \subset \text{不可约对角占优} \cup \text{严格对角占优} \subset \text{H-矩阵} \subset \text{非奇异矩阵}.$$

    更精确地：

    (a) 严格行对角占优矩阵是 H-矩阵。

    (b) 不可约对角占优矩阵是 H-矩阵。

    (c) H-矩阵是非奇异矩阵。

    (d) 反之均不成立。

??? proof "证明"
    (a) $A$ 严格行对角占优 $\implies$ $\mathcal{M}(A)$ 严格行对角占优（$|a_{ii}| > \sum|a_{ij}|$ 直接转化为 $\mathcal{M}(A)$ 的行占优条件）$\implies$ $\mathcal{M}(A)$ 是 M-矩阵（定理 38A.4）$\implies$ $A$ 是 H-矩阵。

    (b) $A$ 不可约对角占优。$\mathcal{M}(A)$ 不可约（$A$ 不可约意味着 $\mathcal{M}(A)$ 不可约——不可约性只取决于零/非零模式，$\mathcal{M}(A)$ 和 $A$ 有相同的零/非零模式）且弱行对角占优（$\ge$）且至少一行严格。

    $\mathcal{M}(A) = sI - B$（$B \ge 0$ 不可约）。$\mathcal{M}(A)\mathbf{1} \ge 0$（弱对角占优），至少一个分量 $> 0$。设 $y = \mathcal{M}(A)\mathbf{1} \ge 0, y \ne 0$。$\mathbf{1} > 0$ 且 $\mathcal{M}(A)\mathbf{1} = y$。需要证明 $\rho(B) < s$。

    $B\mathbf{1} = s\mathbf{1} - y$。$s - (B\mathbf{1})_i/1 = y_i \ge 0$（对所有 $i$），至少一个 $> 0$。由 Collatz-Wielandt 公式对不可约非负矩阵：$\rho(B) \le \max_i (B\mathbf{1})_i = \max_i(s - y_i) < s$（因为至少一个 $y_i > 0$，而不可约性保证 $\rho(B) < \max$）。

    更精确地：由不可约 Perron-Frobenius，$\rho(B) < \max_i (Bx)_i / x_i$ 当且仅当 $Bx \ne \rho(B)x$。取 $x = \mathbf{1}$，$(B\mathbf{1})_i = s - y_i$。若 $B\mathbf{1} = \rho(B)\mathbf{1}$，则 $\rho(B) = s - y_i$ 对所有 $i$，故 $y_i$ 对所有 $i$ 相同。但至少一个 $y_i > 0$ 和一个 $y_j \ge 0$...如果所有 $y_i$ 相同且正，则 $\rho(B) = s - y_1 < s$。如果不全相同（某些为零），则 $B\mathbf{1} \ne c\mathbf{1}$，$\rho(B) < \max_i(B\mathbf{1})_i/1_i = s$。

    综合：$\rho(B) < s$，$\mathcal{M}(A)$ 是非奇异 M-矩阵。

    (c) 由 H-矩阵条件 H4。

    (d) 反例见例 38B.6（H-矩阵但不严格对角占优）。非奇异但不是 H-矩阵的例子：$A = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$，$\mathcal{M}(A) = \begin{pmatrix} 1 & -2 \\ 0 & 1 \end{pmatrix}$，$\det(\mathcal{M}(A)) = 1 > 0$。实际上这里 $\mathcal{M}(A) = I - B$，$B = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}$，$\rho(B) = 0 < 1$，所以 $\mathcal{M}(A)$ 是 M-矩阵，$A$ 是 H-矩阵。取更好的反例：$A = \begin{pmatrix} 1 & 3 \\ 3 & 1 \end{pmatrix}$，$\mathcal{M}(A) = \begin{pmatrix} 1 & -3 \\ -3 & 1 \end{pmatrix}$，$\det(\mathcal{M}(A)) = 1-9 = -8 < 0$。$\mathcal{M}(A)$ 不是 M-矩阵，故 $A$ 不是 H-矩阵。但 $\det(A) = 1-9 = -8 \ne 0$，$A$ 非奇异。$\blacksquare$

!!! example "例 38B.7"
    矩阵类包含关系示意：

    $$\text{对称正定} \cap \text{Z-矩阵} \subset \text{对称 M-矩阵} \subset \text{M-矩阵} \subset \text{P-矩阵}$$
    $$\text{M-矩阵} \subset \text{H-矩阵} \cap \text{Z-矩阵}$$
    $$\text{严格对角占优} \subset \text{H-矩阵} \subset \text{非奇异矩阵}$$

---

## 38B.6 N-矩阵

!!! definition "定义 38B.7 (N-矩阵)"
    矩阵 $A \in M_n(\mathbb{R})$ 称为 **N-矩阵**，若其所有主子式为负：
    $$\det(A[\alpha, \alpha]) < 0, \quad \forall\, \alpha \subseteq \{1, 2, \ldots, n\}, \quad \alpha \ne \emptyset.$$
    这要求所有对角元素 $a_{ii} < 0$。

!!! definition "定义 38B.8 (N$_0$-矩阵)"
    $A$ 称为 **N$_0$-矩阵**，若所有主子式非正。

!!! theorem "定理 38B.8 (N-矩阵的性质)"
    (a) $A$ 是 N-矩阵当且仅当 $-A$ 是 P-矩阵且 $n$ 为奇数时 $\det(-A) > 0$...更精确地：$A$ 是 N-矩阵当且仅当对每个非空 $\alpha$，$(-1)^{|\alpha|}\det(A[\alpha,\alpha]) > 0$，这不完全等价于 $-A$ 是 P-矩阵。

    实际上：$\det((-A)[\alpha,\alpha]) = (-1)^{|\alpha|}\det(A[\alpha,\alpha])$。若 $A$ 是 N-矩阵（$\det(A[\alpha,\alpha]) < 0$），则 $\det((-A)[\alpha,\alpha]) = (-1)^{|\alpha|} \cdot (\text{负数})$。当 $|\alpha|$ 为奇数时为正，当 $|\alpha|$ 为偶数时为负。所以 $-A$ **不是** P-矩阵（偶数阶主子式为负）。

    (b) N-矩阵只有在 $n = 1$ 时简单（$a_{11} < 0$ 即可）。对 $n \ge 2$，N-矩阵的存在要求更精细的条件。

    (c) 若 $A$ 是 N-矩阵，则 $A$ 的所有对角元素为负。

    (d) 若 $A$ 是 N-矩阵且 $n$ 为偶数，则 $A$ 非奇异（$\det(A) < 0$）。若 $n$ 为奇数，$\det(A) < 0$，也非奇异。所以 N-矩阵总是非奇异的。

??? proof "证明"
    (c) 取 $\alpha = \{i\}$，$\det(A[\{i\},\{i\}]) = a_{ii} < 0$。

    (d) $\det(A) = \det(A[\{1,\ldots,n\},\{1,\ldots,n\}]) < 0 \ne 0$。

!!! example "例 38B.8"
    $A = \begin{pmatrix} -1 & 2 \\ 3 & -2 \end{pmatrix}$。

    $a_{11} = -1 < 0$，$a_{22} = -2 < 0$，$\det(A) = 2 - 6 = -4 < 0$。所有主子式为负，$A$ 是 N-矩阵。

!!! theorem "定理 38B.9 (N-矩阵与 M-矩阵的关系)"
    若 $A$ 是 Z-矩阵（非对角元 $\le 0$）且是 N-矩阵（所有主子式 $< 0$），则 $-A$ 是 Z-矩阵且所有对角元素为正。$-A$ 的主子式为 $(-1)^{|\alpha|}\det(A[\alpha,\alpha])$。由 $\det(A[\alpha,\alpha]) < 0$，奇数阶主子式 $(-1)^k \cdot (\text{负}) = (-1)^{k+1}$...

    这种情形比较特殊。一般来说，N-矩阵与 M-矩阵没有简单的包含关系——它们的主子式符号恰好相反。

---

## 38B.7 半正矩阵

!!! definition "定义 38B.9 (半正矩阵)"
    矩阵 $A \in M_n(\mathbb{R})$ 称为**半正**的（semipositive），若存在 $x \ge 0$（$x \ne 0$）使得 $Ax > 0$。

    $A$ 称为**最小半正**的（minimally semipositive），若 $A$ 是半正的，但 $A$ 的任何真列子矩阵都不是半正的。

!!! theorem "定理 38B.10 (半正矩阵的性质)"
    (a) 若存在 $x > 0$（严格正）使得 $Ax > 0$，则 $A$ 是半正的。

    (b) M-矩阵是半正的（条件 M6）。

    (c) 非奇异 M-矩阵是半正的且其半正性通过严格正向量 $x > 0$ 实现。

    (d) $A$ 半正意味着 $A$ 有正的列依赖：$A$ 的列的某个非负线性组合是正向量。

    (e) 半正矩阵不一定非奇异。

??? proof "证明"
    (a) $x > 0 \ge 0$，$x \ne 0$，$Ax > 0$。

    (b) M-矩阵条件 M6。

    (c) M-矩阵条件 M6 给出 $x > 0$，$Ax > 0$。

    (d) $Ax > 0$ 等价于 $\sum_j x_j a_j > 0$（$a_j$ 为 $A$ 的第 $j$ 列），$x \ge 0$, $x \ne 0$。

    (e) 反例：$A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$，取 $x = (1,0)^T$，$Ax = (1,1)^T > 0$。但 $\det A = 0$。$\blacksquare$

!!! example "例 38B.9"
    $A = \begin{pmatrix} 1 & -1 & 2 \\ 0 & 1 & 1 \\ -1 & 0 & 1 \end{pmatrix}$。取 $x = (1, 1, 1)^T$，$Ax = (2, 2, 0)^T$。不满足 $Ax > 0$。

    取 $x = (2, 1, 1)^T$，$Ax = (3, 2, -1)^T$。不满足。

    取 $x = (3, 1, 1)^T$，$Ax = (4, 2, -2)^T$。不满足。

    取 $x = (1, 0, 1)^T$，$Ax = (3, 1, 0)^T$。不满足（第三分量为零）。

    实际上需要找 $x \ge 0$ 使 $Ax > 0$。通过线性规划或直接搜索可以判定半正性。

---

## 38B.8 Ostrowski-Reich 定理

!!! theorem "定理 38B.11 (Ostrowski-Reich 定理)"
    设 $A \in M_n(\mathbb{R})$ 为对称矩阵。SOR 迭代（以松弛参数 $\omega \in (0, 2)$）对 $Ax = b$ 收敛（即 $\rho(T_\omega) < 1$）当且仅当 $A$ 正定。

    对 H-矩阵的推广：设 $A$ 是 H-矩阵。则对所有 $\omega \in (0, 1]$，SOR 迭代收敛。

??? proof "证明"
    **对称情形（Ostrowski-Reich 原始定理）**：

    $A = D - L - L^T$（$D$ 正对角，$L$ 严格下三角）。SOR 迭代矩阵
    $$T_\omega = (D - \omega L)^{-1}((1-\omega)D + \omega L^T).$$

    **充分性**：设 $A$ 正定。定义 $G_\omega = (D - \omega L)^T D^{-1}(D - \omega L)$。

    关键恒等式：$G_\omega - T_\omega^T G_\omega T_\omega = \omega(2 - \omega)(D - \omega L)^T D^{-1} A D^{-1}(D - \omega L)$。

    这可以通过直接代入验证（虽然计算繁杂）。

    当 $0 < \omega < 2$ 时，$\omega(2-\omega) > 0$。$A$ 正定意味着 $D^{-1/2}AD^{-1/2}$ 正定，故右端矩阵是正半定的（实际上是正定的，因为 $(D - \omega L)$ 非奇异）。

    $G_\omega$ 是正定的（$(D - \omega L)$ 非奇异，$D > 0$）。

    $G_\omega - T_\omega^T G_\omega T_\omega$ 正定意味着 $T_\omega$ 相对于内积 $\langle x, G_\omega x \rangle$ 是严格压缩的：$\|T_\omega x\|_{G_\omega}^2 = x^T T_\omega^T G_\omega T_\omega x < x^T G_\omega x = \|x\|_{G_\omega}^2$（对 $x \ne 0$）。故 $\rho(T_\omega) < 1$。

    **必要性**：设 $\rho(T_\omega) < 1$。由 Kahan 的结论，$\det(T_\omega) = (1-\omega)^n$，故需要 $|1-\omega| < 1$，即 $0 < \omega < 2$。

    $(I - T_\omega)^{-1}$ 存在。$I - T_\omega = (D - \omega L)^{-1}[(D - \omega L) - (1-\omega)D - \omega L^T] = (D-\omega L)^{-1}[\omega D - \omega L - \omega L^T] = \omega(D-\omega L)^{-1}A$。

    故 $(I - T_\omega)^{-1} = \omega^{-1}A^{-1}(D - \omega L)$，$A^{-1}$ 存在。

    需证 $A$ 正定。取 $\omega$ 使 $\rho(T_\omega) < 1$，由收敛性和对称性的结合，$A$ 的二次型 $x^T A x$ 可以通过迭代的能量递减性质得到正定性。具体地，由恒等式，$G_\omega - T_\omega^T G_\omega T_\omega = \omega(2-\omega)C^T A C$（$C = D^{-1}(D-\omega L)$）正半定。$G_\omega$ 正定且 $\rho(T_\omega) < 1$ 意味着 $G_\omega - T_\omega^T G_\omega T_\omega$ 正定，从而 $A$ 正定。$\blacksquare$

!!! theorem "定理 38B.12 (H-矩阵的 SOR 收敛)"
    设 $A$ 是 H-矩阵。则：

    (a) 对 $\omega \in (0, 1]$，SOR 迭代矩阵 $T_\omega$ 满足 $\rho(T_\omega) < 1$。

    (b) Gauss-Seidel 迭代（$\omega = 1$）收敛。

    (c) Jacobi 迭代收敛。

??? proof "证明"
    $\mathcal{M}(A)$ 是 M-矩阵。由 M-矩阵的迭代收敛定理（定理 38A.12），$\mathcal{M}(A)$ 的 Jacobi 和 Gauss-Seidel 迭代收敛。

    关键的比较引理：对一般矩阵 $A$，$\rho(|M^{-1}||N|) \le \rho(\mathcal{M}(M)^{-1}\mathcal{M}(N)_+)$...

    更直接地：$A$ 是 H-矩阵意味着存在正对角 $D$ 使 $AD$ 严格行对角占优。$AD$ 的 Jacobi 迭代矩阵谱半径 $< 1$（由对角占优），而对角缩放不改变 SOR（$\omega = 1$）的谱半径。

    对 $\omega \in (0, 1]$，利用 $T_\omega$ 的元素被 $\mathcal{M}(A)$ 的对应 SOR 矩阵控制（逐元），由非负矩阵的谱半径单调性，$\rho(T_\omega(A)) \le \rho(T_\omega(\mathcal{M}(A))) < 1$。$\blacksquare$

---

## 38B.9 广义 M-矩阵（复情形）

!!! definition "定义 38B.10 (广义 M-矩阵)"
    矩阵 $A \in M_n(\mathbb{C})$ 称为**广义 M-矩阵**（GM-矩阵），若 $A$ 可以写成
    $$A = sI - B, \quad B \ge 0, \quad s > \rho(B),$$
    其中 $s \in \mathbb{C}$，$|s| > \rho(B)$，且 $B \ge 0$。

    更常见的定义：$A \in M_n(\mathbb{C})$ 称为广义 M-矩阵，若 $A$ 是 H-矩阵且 $a_{ii} > 0$ 对所有 $i$。

!!! theorem "定理 38B.13 (广义 M-矩阵的性质)"
    (a) 实广义 M-矩阵即通常的 M-矩阵。

    (b) 广义 M-矩阵是 H-矩阵。

    (c) 广义 M-矩阵的比较矩阵是 M-矩阵（由 H-矩阵定义）。

    (d) 广义 M-矩阵 $A$ 满足 $|A^{-1}| \le \mathcal{M}(A)^{-1}$（逐元），其中 $|A^{-1}|$ 表示逐元取模。

??? proof "证明"
    (d) 这是 H-矩阵理论中的重要不等式。由 $A = D_A(I - D_A^{-1}(D_A - A))$，其中 $D_A = \operatorname{diag}(a_{11},\ldots,a_{nn})$。令 $J = D_A^{-1}(D_A - A)$，则 $|J_{ij}| \le |D_A^{-1}|_{ii} |a_{ij}|$（$i \ne j$），$J_{ii} = 0$。

    $A^{-1} = (I - J)^{-1}D_A^{-1} = \sum_{k=0}^{\infty}J^k D_A^{-1}$。$|A^{-1}| \le \sum |J|^k |D_A^{-1}| = (I - |J|)^{-1}|D_A|^{-1}$（逐元）。

    $\mathcal{M}(A) = |D_A|(I - |J'|)$（适当定义 $J'$），$\mathcal{M}(A)^{-1} = (I-|J'|)^{-1}|D_A|^{-1}$。

    精确的论证需要仔细处理对角缩放，但核心思想是 $|A^{-1}|$ 的每个元素被 $\mathcal{M}(A)^{-1}$ 的对应元素控制。$\blacksquare$

---

## 38B.10 全正矩阵的联系

!!! definition "定义 38B.11 (全正矩阵，简述)"
    矩阵 $A \in M_n(\mathbb{R})$ 称为**全正矩阵**（totally positive, TP），若 $A$ 的所有子式（不仅是主子式）为非负的。称为**严格全正**（STP），若所有子式严格为正。

!!! theorem "定理 38B.14 (全正矩阵与 P-矩阵的关系)"
    (a) 严格全正矩阵是 P-矩阵（因为所有主子式为正，作为所有子式为正的特例）。

    (b) 反之不然：P-矩阵只要求主子式为正，不要求非主子式为正。

    (c) 全正矩阵的逆（若存在）具有棋盘符号模式：$(A^{-1})_{ij}$ 的符号为 $(-1)^{i+j}$。特别地，若 $A$ 是全正的 Z-矩阵（可能在适当变换后），则 $A^{-1}$ 也有特殊结构。

    (d) $n \times n$ 振荡矩阵（oscillation matrix）是全正矩阵的一个重要子类，其特征值全部为正单实数，特征向量具有 Sturm 振荡性质。

??? proof "证明"
    (a) 全正矩阵的所有子式 $\ge 0$，严格全正矩阵的所有子式 $> 0$。主子式是子式的特例，故严格全正 $\implies$ 所有主子式 $> 0$ $\implies$ P-矩阵。

    (b) 取 $A = \begin{pmatrix} 1 & -2 \\ 1 & 1 \end{pmatrix}$，P-矩阵（例 38B.1），但 $a_{12} = -2 < 0$，$1 \times 1$ 子式 $a_{12} < 0$，不是全正矩阵。

    (c) 这是全正矩阵理论的经典结果，详见第 39 章。$\blacksquare$

---

## 38B.11 迭代法与数值稳定性应用

!!! theorem "定理 38B.15 (H-矩阵的 LU 分解稳定性)"
    设 $A$ 是 H-矩阵。则：

    (a) 不带行交换的 Gauss 消元（LU 分解）可以完成，且不会遇到零主元。

    (b) 消元过程中的增长因子受控：中间元素的模不超过 $\mathcal{M}(A)^{-1}$ 的相应元素乘以初始元素的模。

    (c) H-矩阵的 LU 分解在数值上是稳定的（无需选主元）。

??? proof "证明"
    (a) $\mathcal{M}(A)$ 是 M-矩阵，条件 M4 保证 $\mathcal{M}(A)$ 的所有顺序主子式为正。由 Ostrowski 不等式，$A$ 的顺序主子式也非零。故 LU 分解存在。

    (b) 消元过程中，Schur 补 $A/A[1:k,1:k]$ 的比较矩阵被 $\mathcal{M}(A)$ 的 Schur 补控制（逐元）。由 $\mathcal{M}(A)$ 的 Schur 补是 M-矩阵（定理 38A.7），增长因子受控。

    (c) 由 (b)，增长因子有界，消元数值稳定。$\blacksquare$

!!! theorem "定理 38B.16 (H-矩阵的区间分析)"
    设 $A$ 是 H-矩阵，$A x = b$。若 $A$ 的元素有扰动 $|\Delta A| \le E$（$E \ge 0$），$b$ 的扰动 $|\Delta b| \le f$（$f \ge 0$），则解的扰动满足
    $$|\Delta x| \le \mathcal{M}(A)^{-1}(E|x| + f).$$
    这为 H-矩阵线性系统的解提供了严格的误差界。

??? proof "证明"
    $(A + \Delta A)(x + \Delta x) = b + \Delta b$。$A\Delta x = \Delta b - \Delta A \cdot x - \Delta A \cdot \Delta x$。忽略高阶项：$\Delta x \approx A^{-1}(\Delta b - \Delta A \cdot x)$。

    $|\Delta x| \le |A^{-1}|(|\Delta b| + |\Delta A||x|) \le |A^{-1}|(f + E|x|)$。

    由广义 M-矩阵性质（定理 38B.13(d)），$|A^{-1}| \le \mathcal{M}(A)^{-1}$。

    故 $|\Delta x| \le \mathcal{M}(A)^{-1}(E|x| + f)$。$\blacksquare$

!!! example "例 38B.10 (预处理与 H-矩阵)"
    在迭代法中，若系数矩阵 $A$ 不是 H-矩阵，可以通过预处理将其转化为 H-矩阵。设 $M$ 为预处理矩阵，$M^{-1}A$ 的比较矩阵若是 M-矩阵，则 $M^{-1}A$ 是 H-矩阵，迭代法对预处理后的系统 $M^{-1}Ax = M^{-1}b$ 收敛。

    选择好的预处理 $M$ 使得 $\mathcal{M}(M^{-1}A)$ 的谱半径（去掉对角后除以对角）尽可能小，是数值线性代数中的重要课题。

---

## 本章小结

| 矩阵类 | 定义 | 核心性质 | 与 M-矩阵的关系 |
|--------|------|----------|-----------------|
| P-矩阵 | 所有主子式 $> 0$ | LCP 唯一可解；P2 符号条件 | $\supset$ M-矩阵 |
| P$_0$-矩阵 | 所有主子式 $\ge 0$ | P-矩阵的闭包 | $\supset$ 奇异 M-矩阵 |
| H-矩阵 | $\mathcal{M}(A)$ 是 M-矩阵 | 广义对角占优；数值稳定 | M-矩阵的推广（允许一般符号和复数） |
| N-矩阵 | 所有主子式 $< 0$ | 对角元全负 | 与 M-矩阵互补 |
| 半正矩阵 | $\exists x \ge 0: Ax > 0$ | M-矩阵的必要条件 | M-矩阵 $\subset$ 半正 |
| 广义 M-矩阵 | 复 H-矩阵 + 正对角元 | $|A^{-1}| \le \mathcal{M}(A)^{-1}$ | 复数推广 |

层次关系总结：
$$\text{对称正定} \cap \text{Z-矩阵} \subset \text{M-矩阵} \subset \text{P-矩阵} \cap \text{Z-矩阵}$$
$$\text{M-矩阵} \subset \text{H-矩阵} \cap \text{Z-矩阵}, \quad \text{严格对角占优} \subset \text{H-矩阵}$$
$$\text{P-矩阵} \supset \text{严格全正矩阵}, \quad \text{P}_0\text{-矩阵} \supset \text{P-矩阵}$$

---

## 习题

!!! example "例 38B.E1"
    证明：若 $A$ 是 P-矩阵，则 $A^T$ 也是 P-矩阵。

    **提示**：$\det(A^T[\alpha,\alpha]) = \det((A[\alpha,\alpha])^T) = \det(A[\alpha,\alpha])$。

!!! example "例 38B.E2"
    设 $A$ 是 P-矩阵。证明 $A$ 的所有实特征值为正。P-矩阵的特征值是否一定有正实部？给出例子。

!!! example "例 38B.E3"
    （LCP）对 $A = \begin{pmatrix} 1 & -1 \\ 1 & 2 \end{pmatrix}$，$q = \begin{pmatrix} -3 \\ -2 \end{pmatrix}$。

    (a) 验证 $A$ 是 P-矩阵。

    (b) 求 $\operatorname{LCP}(q, A)$ 的唯一解。

    (c) 给出解的几何解释。

!!! example "例 38B.E4"
    证明 P$_0$-矩阵是 P-矩阵的闭包：$A$ 是 P$_0$-矩阵当且仅当对每个 $\epsilon > 0$，$A + \epsilon I$ 是 P-矩阵。

!!! example "例 38B.E5"
    设 $A \in M_3(\mathbb{C})$，
    $$A = \begin{pmatrix} 4 & 1-i & -1 \\ -1+2i & 5 & i \\ 2 & -1 & 3 \end{pmatrix}.$$

    (a) 计算比较矩阵 $\mathcal{M}(A)$。

    (b) 判断 $\mathcal{M}(A)$ 是否为 M-矩阵。

    (c) 判断 $A$ 是否为 H-矩阵。

    (d) 若是 H-矩阵，求正对角矩阵 $D$ 使 $AD$ 严格行对角占优。

!!! example "例 38B.E6"
    构造一个 $3 \times 3$ 的 H-矩阵，其不是严格行对角占优的，也不是 Z-矩阵。

!!! example "例 38B.E7"
    （Ostrowski-Reich 定理验证）对对称正定矩阵 $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$，

    (a) 对 $\omega = 0.5, 1.0, 1.5$，计算 SOR 迭代矩阵 $T_\omega$ 及其谱半径。

    (b) 验证 $\rho(T_\omega) < 1$ 对 $\omega \in (0, 2)$ 成立。

    (c) 求最优松弛参数 $\omega^*$（使 $\rho(T_\omega)$ 最小）。

!!! example "例 38B.E8"
    设 $A$ 是 H-矩阵，$\Delta A$ 是满足 $|\Delta A| \le 0.01 \cdot |A|$ 的扰动。利用定理 38B.16 估计 $(A + \Delta A)^{-1}$ 存在的条件。

!!! example "例 38B.E9"
    证明以下包含关系，并对每个包含给出反例说明严格性：

    (a) M-矩阵 $\subset$ P-矩阵。

    (b) M-矩阵 $\subset$ H-矩阵 $\cap$ Z-矩阵。

    (c) 严格对角占优 $\subset$ H-矩阵。

    (d) P-矩阵 $\not\subset$ H-矩阵 且 H-矩阵 $\not\subset$ P-矩阵。

!!! example "例 38B.E10"
    （N-矩阵）设 $A = \begin{pmatrix} -2 & 1 \\ 3 & -1 \end{pmatrix}$。

    (a) 验证 $A$ 是 N-矩阵。

    (b) $A$ 是否为 Z-矩阵？

    (c) 计算 $A^{-1}$ 并分析其符号模式。

!!! example "例 38B.E11"
    （广义 M-矩阵）设 $A = \begin{pmatrix} 2 & -1+i \\ -1-i & 3 \end{pmatrix}$。

    (a) 证明 $A$ 是 H-矩阵。

    (b) 计算 $|A^{-1}|$ 和 $\mathcal{M}(A)^{-1}$，验证 $|A^{-1}| \le \mathcal{M}(A)^{-1}$。

!!! example "例 38B.E12"
    设 $A$ 是 $n \times n$ P-矩阵。证明：

    (a) $A$ 的 Schur 补（关于任何非奇异主子矩阵 $A_{11}$）仍为 P-矩阵。

    (b) 与 M-矩阵的 Schur 补定理（定理 38A.7）比较。P-矩阵的 Schur 补定理不需要 Z-矩阵条件。

    **提示**：利用条件 P2 和 Schur 补的二次型表示。
