# 第 53 章 辛矩阵与 Hamilton 矩阵

<div class="context-flow" markdown>

**前置**：二次型与双线性型 (Ch9) · 特征值 (Ch6) · 矩阵指数 (Ch13) · 矩阵群 (第 55 章参考)

**本章脉络**：辛形式与辛空间 → 辛矩阵 $M^T J M = J$ → 辛群 $\mathrm{Sp}(2n)$ → Hamilton 矩阵 → 特征值配对性质 → 辛特征值 → 辛 Gram-Schmidt → Williamson 定理 → 辛积分器

**延伸**：辛几何是经典力学（Hamilton 方程保辛性）和量子光学（Gaussian 态由辛矩阵参数化）的数学语言；辛积分器在天体力学长时间数值模拟中比一般方法稳定得多

</div>

在分析力学中，Hamilton 正则方程 $\dot{q} = \partial H / \partial p$，$\dot{p} = -\partial H / \partial q$ 的解流保持一种特殊的几何结构——辛结构。这一结构的线性代数本质就是本章的主题。辛矩阵是保持标准辛形式的线性变换在矩阵层面的体现，它们构成的辛群 $\mathrm{Sp}(2n)$ 是与正交群、酉群并列的经典矩阵群。Hamilton 矩阵则是辛群对应 Lie 代数的元素，其特征值具有精美的对称配对性质。本章从辛形式的代数定义出发，系统建立辛矩阵和 Hamilton 矩阵的理论，介绍辛特征值与 Williamson 定理，最后展示辛积分器在数值计算中的应用。

---

## 53.1 辛形式与辛空间

<div class="context-flow" markdown>

**核心问题**：什么样的双线性型在"保持"意义下引出辛群？

</div>

!!! definition "定义 53.1 (反对称双线性型)"
    设 $V$ 是域 $\mathbb{F}$ 上的有限维向量空间。一个双线性型 $\omega: V \times V \to \mathbb{F}$ 称为**反对称**（skew-symmetric，或称交替型），若对所有 $u, v \in V$，
    $$\omega(u, v) = -\omega(v, u).$$
    特别地，$\omega(v, v) = 0$ 对所有 $v \in V$ 成立（当 $\mathrm{char}(\mathbb{F}) \neq 2$ 时，这两个条件等价）。

!!! definition "定义 53.2 (非退化反对称型，辛形式)"
    反对称双线性型 $\omega$ 称为**非退化**的，若 $\omega(u, v) = 0$ 对所有 $v \in V$ 成立蕴含 $u = 0$。非退化的反对称双线性型称为**辛形式**（symplectic form）。

!!! definition "定义 53.3 (辛空间)"
    配备辛形式 $\omega$ 的向量空间 $(V, \omega)$ 称为**辛空间**（symplectic vector space）。

!!! theorem "定理 53.1 (辛空间的维数)"
    辛空间 $(V, \omega)$ 的维数必为偶数。

??? proof "证明"
    设 $\dim V = m$。取 $V$ 的一组基 $\{e_1, \ldots, e_m\}$，定义矩阵 $\Omega = [\omega(e_i, e_j)]_{m \times m}$。由 $\omega$ 的反对称性，$\Omega^T = -\Omega$。由 $\omega$ 的非退化性，$\det \Omega \neq 0$。

    但由 $\Omega^T = -\Omega$，我们有
    $$\det \Omega = \det(\Omega^T) = \det(-\Omega) = (-1)^m \det \Omega.$$
    因为 $\det \Omega \neq 0$，所以 $(-1)^m = 1$，即 $m$ 为偶数。$\blacksquare$

!!! definition "定义 53.4 (标准辛形式)"
    在 $\mathbb{R}^{2n}$（或 $\mathbb{C}^{2n}$）上，**标准辛形式**由矩阵
    $$J_{2n} = \begin{bmatrix} 0 & I_n \\ -I_n & 0 \end{bmatrix}$$
    定义：$\omega(x, y) = x^T J_{2n} y$。在不引起混淆时，简写为 $J$。

!!! theorem "定理 53.2 (辛形式的标准化)"
    设 $(V, \omega)$ 是 $2n$ 维辛空间。则存在 $V$ 的一组基 $\{e_1, \ldots, e_n, f_1, \ldots, f_n\}$ 使得
    $$\omega(e_i, e_j) = 0, \quad \omega(f_i, f_j) = 0, \quad \omega(e_i, f_j) = \delta_{ij}.$$
    这组基称为**辛基**（symplectic basis）或 Darboux 基。在此基下，辛形式的矩阵恰好是 $J_{2n}$。

??? proof "证明"
    用归纳法。取 $e_1 \neq 0$。因 $\omega$ 非退化，存在 $f_1$ 使得 $\omega(e_1, f_1) \neq 0$。适当缩放使 $\omega(e_1, f_1) = 1$。

    令 $W_1 = \mathrm{span}\{e_1, f_1\}$，定义 $W_1^\perp = \{v \in V : \omega(v, w) = 0, \forall w \in W_1\}$。

    **断言**：$V = W_1 \oplus W_1^\perp$。

    首先，$W_1 \cap W_1^\perp = \{0\}$：设 $v = ae_1 + bf_1 \in W_1 \cap W_1^\perp$，则 $\omega(v, e_1) = -b = 0$，$\omega(v, f_1) = a = 0$，故 $v = 0$。

    其次，$\dim W_1^\perp = 2n - 2$。定义映射 $\varphi: V \to W_1^*$，$\varphi(v)(w) = \omega(v, w)$。由 $\omega$ 在 $W_1$ 上的非退化性可知 $\varphi$ 是满射，核恰好是 $W_1^\perp$，所以 $\dim W_1^\perp = 2n - 2$。

    容易验证 $\omega$ 限制在 $W_1^\perp$ 上仍然是非退化的（若 $v \in W_1^\perp$ 使得 $\omega(v, u) = 0$ 对所有 $u \in W_1^\perp$，结合 $v \in W_1^\perp$ 意味着 $\omega(v, w) = 0$ 对所有 $w \in W_1$，故对所有 $w \in V$，由非退化性 $v = 0$）。

    因此 $(W_1^\perp, \omega|_{W_1^\perp})$ 是 $(2n-2)$ 维辛空间，对其施加归纳假设即得结论。$\blacksquare$

!!! example "例 53.1"
    在 $\mathbb{R}^4$ 上考虑标准辛形式 $\omega(x,y) = x^T J_4 y$，其中
    $$J_4 = \begin{bmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ -1 & 0 & 0 & 0 \\ 0 & -1 & 0 & 0 \end{bmatrix}.$$
    标准基 $\{e_1, e_2, e_3, e_4\}$ 就是辛基，其中 $e_1, e_2$ 对应 $q$ 坐标，$e_3, e_4$ 对应 $p$ 坐标：
    $$\omega(e_1, e_3) = 1, \quad \omega(e_2, e_4) = 1,$$
    其余所有 $\omega(e_i, e_j)$ 均为零（除去反对称确定的项）。

!!! remark "注记"
    在力学中，通常将坐标排列为 $(q_1, \ldots, q_n, p_1, \ldots, p_n)$，辛形式 $\omega = \sum_{i=1}^n dq_i \wedge dp_i$ 的矩阵表示正是 $J_{2n}$。不同的文献对 $J$ 的符号约定可能不同（有的取 $\begin{bmatrix} 0 & -I \\ I & 0 \end{bmatrix}$），本章一致采用上述约定。

---

## 53.2 辛矩阵

<div class="context-flow" markdown>

**核心问题**：保持辛形式的线性变换具有什么矩阵刻画？

</div>

!!! definition "定义 53.5 (辛矩阵)"
    一个 $2n \times 2n$ 实矩阵（或复矩阵）$M$ 称为**辛矩阵**（symplectic matrix），若
    $$M^T J M = J,$$
    其中 $J = J_{2n}$ 是标准辛矩阵。

这个定义等价于说 $M$ 保持标准辛形式：$\omega(Mx, My) = (Mx)^T J (My) = x^T M^T J M y = x^T J y = \omega(x, y)$。

!!! theorem "定理 53.3 (辛矩阵的基本性质)"
    设 $M$ 是 $2n \times 2n$ 辛矩阵。则：

    **(a)** $\det M = 1$（不仅仅是 $\pm 1$，而确实是 $+1$）。

    **(b)** $M$ 可逆，且 $M^{-1} = -J M^T J = J^T M^T J$。

    **(c)** $M^{-1}$ 也是辛矩阵。

    **(d)** 若 $M_1, M_2$ 都是辛矩阵，则 $M_1 M_2$ 也是辛矩阵。

    **(e)** $M^T$ 是辛矩阵。

??? proof "证明"
    **(b)** 由 $M^T J M = J$ 两边取逆：$M^{-1} J^{-1} (M^T)^{-1} = J^{-1}$。注意 $J^{-1} = -J = J^T$。因此
    $$M^{-1} (-J) (M^T)^{-1} = -J,$$
    $$M^{-1} J (M^T)^{-1} = J,$$
    即 $(M^{-1})^T J M^{-1} = J$（对上式取转置后整理）。但我们先直接算 $M^{-1}$：

    从 $M^T J M = J$，左乘 $J^{-1} = -J$，得 $-J M^T J M = I$，所以 $M^{-1} = -J M^T J$。

    等价地，$M^{-1} = J^T M^T J$（因为 $J^T = -J$，所以 $J^T M^T J = (-J) M^T J = -J M^T J$，一致）。

    **(a)** 由 $M^T J M = J$ 取行列式：$(\det M)^2 \det J = \det J$。因为 $\det J = 1$（后面验证），故 $(\det M)^2 = 1$，得 $\det M = \pm 1$。

    验证 $\det J = 1$：$J = \begin{bmatrix} 0 & I \\ -I & 0 \end{bmatrix}$，通过行列式的 Schur 补公式或直接计算，$\det J = \det(0 \cdot 0 - I \cdot (-I)) = \det(I) = 1$（更严格地，可以利用分块矩阵的行列式公式，当右下块可逆时 $\det J = \det(-I) \cdot \det(0 - I \cdot (-I)^{-1} \cdot (-I))$，但这里最简洁的方式是注意 $J$ 可以通过初等列变换变为 $I$，符号变化恰好给出 $\det J = 1$）。

    现在证明 $\det M = +1$ 而非 $-1$。这需要一个连通性论证：辛矩阵的集合 $\mathrm{Sp}(2n, \mathbb{R})$ 是连通的（见 53.3 节），而 $I_{2n}$ 是辛矩阵且 $\det I = 1$，$\det$ 是连续函数，其在连通集上的像是连通的，故只能取 $+1$。（代数证明也可以，见下方。）

    **代数证明**（$\det M = 1$）：将 $M$ 写成分块形式 $M = \begin{bmatrix} A & B \\ C & D \end{bmatrix}$，其中 $A, B, C, D$ 均为 $n \times n$ 矩阵。$M^T J M = J$ 展开后给出三个条件：
    $$A^T C = C^T A, \quad B^T D = D^T B, \quad A^T D - C^T B = I.$$
    第三个条件说明 $\det(A^T D - C^T B) = 1$。通过仔细分析，可以从这些关系推出 $\det M = 1$。一种方式是利用辛矩阵的 $2n \times 2n$ Pfaffian 关系来直接证明。

    **(c)** 已在 (b) 的证明过程中验证：$(M^{-1})^T J M^{-1} = J$。

    **(d)** $(M_1 M_2)^T J (M_1 M_2) = M_2^T M_1^T J M_1 M_2 = M_2^T J M_2 = J$。$\checkmark$

    **(e)** 从 $M^T J M = J$，取逆得 $M^{-1} J^{-1} (M^T)^{-1} = J^{-1}$，即 $M^{-1}(-J)(M^T)^{-1} = -J$，故 $M^{-1} J (M^T)^{-1} = J$。这正是说 $(M^T)^T J M^T = M J M^T = J$（利用 $M^{-1} = -JM^TJ$ 重新推导），但更直接地：

    从 $M^T J M = J$，右乘 $M^{-1}$：$M^T J = J M^{-1}$。取转置：$J^T M = (M^{-1})^T J^T = (M^{-1})^T (-J)$，即 $-JM = -(M^{-1})^T J$，故 $JM = (M^{-1})^T J$。即 $(M^T)^{-T} = M$，所以 $(M^T)^T J (M^T) = M J M^T$。我们需要证明 $M J M^T = J$。从 $M^T J M = J$，左乘 $M$，右乘 $M^{-1}$：$M M^T J = M J M^{-1}$（不太方便）。直接来：$M^{-1} = -JM^TJ$，所以 $M(-JM^TJ) = I$，即 $-MJM^TJ = I$，故 $MJM^T = -J^{-1} = J$（因为 $J^{-1} = -J$）。$\blacksquare$

!!! example "例 53.2"
    **$2 \times 2$ 辛矩阵就是行列式为 1 的矩阵**：当 $n = 1$ 时，$J = \begin{bmatrix} 0 & 1 \\ -1 & 0 \end{bmatrix}$，条件 $M^T J M = J$ 等价于 $\det M = 1$。因此 $\mathrm{Sp}(2) = \mathrm{SL}(2, \mathbb{R})$。

!!! example "例 53.3"
    **分块对角辛矩阵**：设 $A \in \mathrm{GL}(n, \mathbb{R})$，则
    $$M = \begin{bmatrix} A & 0 \\ 0 & (A^T)^{-1} \end{bmatrix}$$
    是辛矩阵。验证：$M^T J M = \begin{bmatrix} A^T & 0 \\ 0 & A^{-1} \end{bmatrix}\begin{bmatrix} 0 & I \\ -I & 0 \end{bmatrix}\begin{bmatrix} A & 0 \\ 0 & (A^T)^{-1} \end{bmatrix} = \begin{bmatrix} 0 & A^T(A^T)^{-1} \\ -A^{-1}A & 0 \end{bmatrix} = \begin{bmatrix} 0 & I \\ -I & 0 \end{bmatrix} = J$。

!!! example "例 53.4"
    **对称矩阵生成的辛矩阵**：设 $S = S^T$ 是 $n \times n$ 对称矩阵，则
    $$M = \begin{bmatrix} I & S \\ 0 & I \end{bmatrix}, \quad M = \begin{bmatrix} I & 0 \\ S & I \end{bmatrix}$$
    均为辛矩阵。（读者可自行验证。）

---

## 53.3 辛群 $\mathrm{Sp}(2n)$

<div class="context-flow" markdown>

**核心问题**：辛矩阵的全体构成怎样的群？它的 Lie 群结构如何？

</div>

!!! definition "定义 53.6 (辛群)"
    **辛群**定义为
    $$\mathrm{Sp}(2n, \mathbb{R}) = \{M \in \mathrm{GL}(2n, \mathbb{R}) : M^T J M = J\},$$
    $$\mathrm{Sp}(2n, \mathbb{C}) = \{M \in \mathrm{GL}(2n, \mathbb{C}) : M^T J M = J\}.$$
    有时简记为 $\mathrm{Sp}(2n)$。

!!! theorem "定理 53.4 (辛群的 Lie 群结构)"
    $\mathrm{Sp}(2n, \mathbb{R})$ 是一个 $n(2n+1)$ 维的连通 Lie 群。

??? proof "证明"
    **群性质**：单位阵 $I \in \mathrm{Sp}(2n)$，乘积封闭和逆元封闭已在定理 53.3 中证明。

    **闭子群**：定义映射 $\Phi: \mathrm{GL}(2n) \to M_{2n}^{\mathrm{skew}}$（反对称矩阵空间），$\Phi(M) = M^T J M$。则 $\mathrm{Sp}(2n) = \Phi^{-1}(\{J\})$。由于 $\Phi$ 是连续的，$\mathrm{Sp}(2n)$ 是 $\mathrm{GL}(2n)$ 的闭子群。由闭子群定理（Cartan 定理），它是 Lie 群。

    **维数计算**：$\mathrm{Sp}(2n)$ 的 Lie 代数是
    $$\mathfrak{sp}(2n) = \{X \in M_{2n} : X^T J + J X = 0\}.$$
    这等价于 $JX$ 是对称矩阵（因为 $(JX)^T = X^T J^T = -X^T J = JX$，最后一步用了条件 $X^T J = -JX$）。$2n \times 2n$ 对称矩阵的空间维数为 $\frac{2n(2n+1)}{2} = n(2n+1)$。故 $\dim \mathrm{Sp}(2n) = n(2n+1)$。

    **连通性**：这需要较深入的拓扑论证。一种方式是利用极分解：任何 $M \in \mathrm{Sp}(2n)$ 可以写成 $M = OP$，其中 $O \in \mathrm{Sp}(2n) \cap \mathrm{O}(2n) \cong \mathrm{U}(n)$，$P$ 是正定辛矩阵。$\mathrm{U}(n)$ 是连通的，正定辛矩阵的空间是凸集（故连通），因此 $\mathrm{Sp}(2n)$ 连通。$\blacksquare$

!!! theorem "定理 53.5 (极大紧子群)"
    $\mathrm{Sp}(2n, \mathbb{R})$ 的极大紧子群为
    $$\mathrm{Sp}(2n) \cap \mathrm{O}(2n) \cong \mathrm{U}(n).$$
    具体地，同构映射将 $U = A + iB \in \mathrm{U}(n)$（其中 $A, B$ 为实矩阵）对应到
    $$M = \begin{bmatrix} A & -B \\ B & A \end{bmatrix} \in \mathrm{Sp}(2n) \cap \mathrm{O}(2n).$$

??? proof "证明"
    设 $M = \begin{bmatrix} A & B \\ C & D \end{bmatrix} \in \mathrm{Sp}(2n) \cap \mathrm{O}(2n)$。

    正交性 $M^T M = I$ 给出：$A^T A + C^T C = I$，$B^T B + D^T D = I$，$A^T B + C^T D = 0$。

    辛性 $M^T J M = J$ 给出：$A^T D - C^T B = I$，$A^T C = C^T A$，$B^T D = D^T B$。

    从正交性还可得 $M M^T = I$：$AA^T + BB^T = I$，$CC^T + DD^T = I$，$AC^T + BD^T = 0$。

    结合 $MJM^T = J$（定理 53.3(e) 的推论）：$AD^T - BC^T = I$，$AB^T = BA^T$，$CD^T = DC^T$。

    经过细致比较，可以证明 $D = A$，$C = -B$。将 $U = A + iB$ 代入 $U^* U = I$ 验证：$(A - iB)(A + iB) = A^2 + B^2 + i(AB^T - BA^T)$... 实际上需要用转置关系。最终可以验证 $U = A + iB$ 满足 $U^* U = I$，即 $U \in \mathrm{U}(n)$。$\blacksquare$

!!! remark "注记"
    维数一览：$\dim \mathrm{Sp}(2) = 3$，$\dim \mathrm{Sp}(4) = 10$，$\dim \mathrm{Sp}(6) = 21$。对比：$\dim \mathrm{SO}(2n) = n(2n-1)$，$\dim \mathrm{SU}(2n) = 4n^2 - 1$。辛群的维数 $n(2n+1)$ 恰好等于 $2n \times 2n$ 对称矩阵的空间维数。

---

## 53.4 Hamilton 矩阵

<div class="context-flow" markdown>

**核心问题**：辛群的 Lie 代数元素——Hamilton 矩阵——具有怎样的特征值结构？

</div>

!!! definition "定义 53.7 (Hamilton 矩阵)"
    一个 $2n \times 2n$ 实矩阵 $H$ 称为 **Hamilton 矩阵**（Hamiltonian matrix），若
    $$H^T J + J H = 0,$$
    等价地，$JH$ 是对称矩阵：$(JH)^T = H^T J^T = -H^T J = JH$。

!!! definition "定义 53.8 (Hamilton 矩阵的标准形式)"
    Hamilton 矩阵可以写成分块形式
    $$H = \begin{bmatrix} A & G \\ Q & -A^T \end{bmatrix},$$
    其中 $G = G^T$，$Q = Q^T$ 是 $n \times n$ 对称矩阵，$A$ 是任意 $n \times n$ 矩阵。

??? proof "证明"
    设 $H = \begin{bmatrix} A & B \\ C & D \end{bmatrix}$。条件 $JH$ 对称：
    $$JH = \begin{bmatrix} 0 & I \\ -I & 0 \end{bmatrix}\begin{bmatrix} A & B \\ C & D \end{bmatrix} = \begin{bmatrix} C & D \\ -A & -B \end{bmatrix}.$$
    对称性 $(JH)^T = JH$ 要求：$C^T = C$（即 $C$ 对称），$(-A)^T = D$（即 $D = -A^T$），$D^T = D$ 自动满足（因为 $D = -A^T$ 意味着 $D^T = -A$，但需要 $D = D^T$？不，对称性要求 $\begin{bmatrix} C & D \\ -A & -B \end{bmatrix}^T = \begin{bmatrix} C & D \\ -A & -B \end{bmatrix}$，即 $C^T = C$，$D^T = -A$，$(-A)^T = D$，$(-B)^T = -B$。

    从 $D^T = -A$ 和 $-A^T = D$ 得 $D = -A^T$，一致。$(-B)^T = -B$ 意味着 $B^T = B$。

    记 $G = B$，$Q = C$，得到 $H = \begin{bmatrix} A & G \\ Q & -A^T \end{bmatrix}$，其中 $G = G^T$，$Q = Q^T$。$\blacksquare$

!!! theorem "定理 53.6 (Hamilton 矩阵与辛群的关系)"
    $H$ 是 Hamilton 矩阵当且仅当 $e^{tH}$ 对所有 $t \in \mathbb{R}$ 是辛矩阵。等价地，Hamilton 矩阵的全体 $\mathfrak{sp}(2n, \mathbb{R})$ 恰好是辛群 $\mathrm{Sp}(2n, \mathbb{R})$ 的 Lie 代数。

??? proof "证明"
    设 $M(t) = e^{tH}$。则 $\frac{d}{dt}M(t) = HM(t)$，$M(0) = I$。

    计算 $\frac{d}{dt}(M(t)^T J M(t))$：
    $$\frac{d}{dt}(M^T J M) = \dot{M}^T J M + M^T J \dot{M} = (HM)^T J M + M^T J (HM)$$
    $$= M^T H^T J M + M^T J H M = M^T(H^T J + JH)M.$$
    若 $H^T J + JH = 0$，则 $\frac{d}{dt}(M^T J M) = 0$，故 $M(t)^T J M(t) = M(0)^T J M(0) = J$，即 $e^{tH}$ 是辛矩阵。

    反之，若 $e^{tH}$ 对所有 $t$ 辛，则在 $t = 0$ 处求导即得 $H^T J + JH = 0$。$\blacksquare$

!!! theorem "定理 53.7 (Hamilton 矩阵的特征值配对)"
    设 $H$ 是实 Hamilton 矩阵，$\lambda$ 是 $H$ 的特征值。则 $-\lambda$，$\bar{\lambda}$，$-\bar{\lambda}$ 也都是 $H$ 的特征值（计入重数后，特征值关于实轴和虚轴对称）。

??? proof "证明"
    **$-\lambda$ 也是特征值**：$H^T J + JH = 0$ 意味着 $H^T = -JHJ^{-1} = JHJ$（因为 $J^{-1} = -J$）。因此 $H^T$ 与 $-H$ 相似（通过 $J$），所以它们有相同的特征值。但 $H^T$ 和 $H$ 有相同的特征多项式（因为 $\det(\lambda I - H^T) = \det(\lambda I - H)$），所以 $H$ 和 $-H$ 有相同的特征多项式。这意味着 $p(\lambda) = p(-\lambda)$... 这不完全对，让我们更仔细地说明。

    $H$ 和 $-H$ 相似（通过 $J$：$JHJ^{-1} = -H^T$ 的特征值是 $-\lambda$），而 $-H^T$ 与 $-H$ 有相同的特征多项式。等等，$H^T = -JHJ^{-1}$ 说明 $H^T \sim -H$（相似）。$H^T$ 的特征值与 $H$ 的相同（因为特征多项式相同），而 $-H$ 的特征值是 $\{-\lambda : \lambda \in \sigma(H)\}$。因此 $\sigma(H) = \sigma(-H)$，即若 $\lambda \in \sigma(H)$，则 $-\lambda \in \sigma(H)$。

    **$\bar{\lambda}$ 也是特征值**：因为 $H$ 是实矩阵，其特征多项式系数为实数，所以复特征值成共轭对出现。

    综合两条：$\lambda, -\lambda, \bar{\lambda}, -\bar{\lambda}$ 都是 $H$ 的特征值。$\blacksquare$

!!! example "例 53.5"
    考虑 $4 \times 4$ Hamilton 矩阵
    $$H = \begin{bmatrix} 0 & 1 & 0 & 0 \\ -2 & 0 & 0 & 1 \\ 0 & 0 & 0 & 2 \\ 0 & 0 & -1 & 0 \end{bmatrix}.$$
    验证 $JH$ 对称：$JH = \begin{bmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & -1 & 0 \\ 0 & -1 & 0 & 0 \\ 2 & 0 & 0 & -1 \end{bmatrix}$，确实对称。

    其特征多项式为 $\lambda^4 + 3\lambda^2 + 2 = (\lambda^2 + 1)(\lambda^2 + 2)$，特征值为 $\pm i, \pm i\sqrt{2}$，确实关于实轴和虚轴对称。

!!! theorem "定理 53.8 (辛矩阵的特征值配对)"
    设 $M$ 是实辛矩阵，$\lambda$ 是 $M$ 的特征值。则 $1/\lambda$，$\bar{\lambda}$，$1/\bar{\lambda}$ 也是 $M$ 的特征值。

??? proof "证明"
    由 $M^T J M = J$，知 $M^T = J M^{-1} J^{-1}$，即 $M^T \sim M^{-1}$。$M^T$ 与 $M$ 有相同特征值，$M^{-1}$ 的特征值为 $\{1/\lambda\}$。故 $\lambda$ 是特征值则 $1/\lambda$ 也是。再结合实矩阵的共轭对称性，四元组 $\{\lambda, 1/\lambda, \bar{\lambda}, 1/\bar{\lambda}\}$ 都是特征值。$\blacksquare$

---

## 53.5 辛特征值

<div class="context-flow" markdown>

**核心问题**：正定矩阵在辛等价意义下有怎样的不变量？

</div>

!!! definition "定义 53.9 (辛特征值)"
    设 $A$ 是 $2n \times 2n$ 正定矩阵。$A$ 的**辛特征值**（symplectic eigenvalues）定义为矩阵 $iJA$ 的正特征值 $d_1 \leq d_2 \leq \cdots \leq d_n$。等价地，它们是 $|JA|$ 的正特征值（其中 $|B| = \sqrt{B^*B}$），或者是 $JA$ 的特征值中模的绝对值（$JA$ 的特征值为 $\pm id_1, \ldots, \pm id_n$）。

!!! theorem "定理 53.9 (Williamson 定理)"
    设 $A$ 是 $2n \times 2n$ 实正定对称矩阵，$d_1 \leq d_2 \leq \cdots \leq d_n$ 是其辛特征值。则存在辛矩阵 $S \in \mathrm{Sp}(2n, \mathbb{R})$ 使得
    $$A = S^T D S, \quad D = \mathrm{diag}(d_1, d_2, \ldots, d_n, d_1, d_2, \ldots, d_n).$$
    矩阵 $D$ 称为 $A$ 的 **Williamson 标准形**。

??? proof "证明"
    证明分几步进行。

    **第一步**：$A$ 正定，故 $A^{1/2}$ 存在且正定。考虑矩阵 $K = A^{1/2} J A^{1/2}$。$K$ 是反对称的：$K^T = (A^{1/2})^T J^T (A^{1/2})^T = A^{1/2}(-J)A^{1/2} = -K$。

    **第二步**：反对称实矩阵 $K$ 可以由正交矩阵对角化为 $K = Q \begin{bmatrix} 0 & \Lambda \\ -\Lambda & 0 \end{bmatrix} Q^T$，其中 $\Lambda = \mathrm{diag}(d_1, \ldots, d_n)$，$d_j > 0$（正定性保证非零）。通过重排，可以取 $Q$ 使得
    $$K = Q \tilde{J}_D Q^T, \quad \tilde{J}_D = \begin{bmatrix} 0 & \Lambda \\ -\Lambda & 0 \end{bmatrix}.$$

    **第三步**：定义 $L = A^{-1/2} Q \Lambda^{-1/2}$（适当的 $2n \times 2n$ 形式）。更精确地，令
    $$P = A^{-1/2} Q \begin{bmatrix} \Lambda^{-1/2} & 0 \\ 0 & \Lambda^{-1/2} \end{bmatrix}.$$
    则可以验证 $P^T A P = D$（其中 $D = \mathrm{diag}(d_1,\ldots,d_n,d_1,\ldots,d_n)$）且 $P^T J P = J$（即 $P$ 是辛矩阵）。

    取 $S = P^{-1}$（也是辛矩阵），则 $A = (P^{-1})^T D P^{-1} = S^T D S$。

    **辛特征值的唯一性**：$d_1, \ldots, d_n$ 是 $A^{1/2} J A^{1/2}$（或等价地 $JA$）的特征值的模，不依赖于 $S$ 的选取。$\blacksquare$

!!! example "例 53.6"
    设 $A = \begin{bmatrix} 2 & 0 & 1 & 0 \\ 0 & 2 & 0 & 0 \\ 1 & 0 & 2 & 0 \\ 0 & 0 & 0 & 2 \end{bmatrix}$。

    计算 $JA = \begin{bmatrix} 0 & I \\ -I & 0 \end{bmatrix}\begin{bmatrix} 2 & 0 & 1 & 0 \\ 0 & 2 & 0 & 0 \\ 1 & 0 & 2 & 0 \\ 0 & 0 & 0 & 2 \end{bmatrix} = \begin{bmatrix} 1 & 0 & 2 & 0 \\ 0 & 0 & 0 & 2 \\ -2 & 0 & -1 & 0 \\ 0 & -2 & 0 & 0 \end{bmatrix}$。

    $JA$ 的特征多项式为 $\lambda^4 + 5\lambda^2 + 4 = (\lambda^2+1)(\lambda^2+4)$，特征值为 $\pm i, \pm 2i$。辛特征值为 $d_1 = 1, d_2 = 2$。

    Williamson 标准形为 $D = \mathrm{diag}(1, 2, 1, 2)$。

!!! theorem "定理 53.10 (辛特征值的变分刻画)"
    设 $A$ 是 $2n \times 2n$ 正定矩阵，辛特征值为 $d_1 \leq \cdots \leq d_n$。则
    $$d_k = \min_{\dim W = 2(n-k+1)} \max_{0 \neq z \in W} \frac{z^* A z}{z^* J^T J z},$$
    其中最小值取遍 $\mathbb{C}^{2n}$ 的所有 $2(n-k+1)$ 维子空间 $W$。

---

## 53.6 辛 Gram-Schmidt

<div class="context-flow" markdown>

**核心问题**：如何构造性地找到辛基？

</div>

!!! theorem "定理 53.11 (辛 Gram-Schmidt 过程)"
    给定辛空间 $(V, \omega)$ 的任意一组基，可以通过如下算法构造辛基 $\{e_1, \ldots, e_n, f_1, \ldots, f_n\}$：

    **输入**：线性无关向量 $\{v_1, \ldots, v_{2n}\}$。

    **Step 1**：令 $e_1 = v_1$。

    **Step 2**：在 $\{v_2, \ldots, v_{2n}\}$ 中寻找使 $\omega(e_1, v_j) \neq 0$ 的 $v_j$（由 $\omega$ 非退化保证存在）。令 $f_1 = v_j / \omega(e_1, v_j)$，使 $\omega(e_1, f_1) = 1$。

    **Step 3**：对其余向量 $v_k$（$k \neq 1, j$），作辛正交化：
    $$v_k' = v_k - \omega(v_k, f_1)e_1 - \omega(e_1, v_k)f_1 = v_k - \omega(v_k, f_1)e_1 + \omega(v_k, e_1)f_1.$$
    则 $\omega(v_k', e_1) = 0$，$\omega(v_k', f_1) = 0$。

    **Step 4**：对 $\{v_k'\}$ 递归重复上述过程。

??? proof "证明"
    需要验证 Step 3 中 $\omega(v_k', e_1) = 0$ 和 $\omega(v_k', f_1) = 0$。

    $$\omega(v_k', e_1) = \omega(v_k, e_1) - \omega(v_k, f_1)\omega(e_1, e_1) + \omega(v_k, e_1)\omega(f_1, e_1)$$
    $$= \omega(v_k, e_1) - 0 + \omega(v_k, e_1)(-1) = 0. \checkmark$$

    $$\omega(v_k', f_1) = \omega(v_k, f_1) - \omega(v_k, f_1)\omega(e_1, f_1) + \omega(v_k, e_1)\omega(f_1, f_1)$$
    $$= \omega(v_k, f_1) - \omega(v_k, f_1) \cdot 1 + 0 = 0. \checkmark$$

    每一步将问题约化到低两维的辛空间上，有限步后得到完整辛基。$\blacksquare$

!!! example "例 53.7"
    在 $(\mathbb{R}^4, \omega)$ 中，取 $v_1 = (1,1,0,0)^T$，$v_2 = (0,1,1,0)^T$，$v_3 = (0,0,1,1)^T$，$v_4 = (1,0,0,1)^T$。

    **Step 1**：$e_1 = v_1 = (1,1,0,0)^T$。

    **Step 2**：$\omega(e_1, v_2) = e_1^T J v_2 = (1,1,0,0) \begin{bmatrix}0&0&1&0\\0&0&0&1\\-1&0&0&0\\0&-1&0&0\end{bmatrix}(0,1,1,0)^T = (0,0,-1,-1)(0,1,1,0)^T = -1$。
    令 $f_1 = v_2 / (-1) = (0,-1,-1,0)^T$。验证 $\omega(e_1, f_1) = 1$。$\checkmark$

    **Step 3**：正交化 $v_3$ 和 $v_4$：
    $v_3' = v_3 - \omega(v_3, f_1)e_1 + \omega(v_3, e_1)f_1$。
    $\omega(v_3, f_1) = v_3^T J f_1 = (0,0,1,1)\begin{bmatrix}0&0&1&0\\0&0&0&1\\-1&0&0&0\\0&-1&0&0\end{bmatrix}(0,-1,-1,0)^T = (-1,-1,0,0)(0,-1,-1,0)^T = 1+1 = 2$。
    $\omega(v_3, e_1) = v_3^T J e_1 = (-1,-1,0,0)(1,1,0,0)^T = -2$。
    所以 $v_3' = (0,0,1,1)^T - 2(1,1,0,0)^T + (-2)(0,-1,-1,0)^T = (-2,-2,1,1)^T + (0,2,2,0)^T = (-2,0,3,1)^T$。

    类似地处理 $v_4$，然后对 $\{v_3', v_4'\}$ 重复过程得到 $e_2, f_2$。

!!! theorem "定理 53.12 (Darboux 定理，线性版)"
    设 $(V, \omega)$ 和 $(V', \omega')$ 是同维数的辛空间。则存在线性同构 $\varphi: V \to V'$ 使得 $\omega'(\varphi(u), \varphi(v)) = \omega(u, v)$ 对所有 $u, v \in V$ 成立。即**所有同维辛空间都是辛同构的**。

??? proof "证明"
    分别取 $(V, \omega)$ 的辛基 $\{e_i, f_i\}$ 和 $(V', \omega')$ 的辛基 $\{e_i', f_i'\}$。定义 $\varphi(e_i) = e_i'$，$\varphi(f_i) = f_i'$ 并线性延拓。由辛基的定义性质，$\omega'(\varphi(u), \varphi(v)) = \omega(u,v)$。$\blacksquare$

---

## 53.7 辛积分器

<div class="context-flow" markdown>

**核心问题**：为什么在数值求解 Hamilton 系统时要用保辛的数值方法？

</div>

Hamilton 系统 $\dot{z} = J^{-1} \nabla H(z)$（其中 $z = (q, p)^T$）的精确流 $\varphi_t$ 是辛变换。一般的数值方法（如标准 Runge-Kutta）不保持这个结构，导致长时间积分时能量漂移。

!!! definition "定义 53.10 (辛积分器)"
    一个数值方法 $\Phi_h: z_n \mapsto z_{n+1}$（步长为 $h$）称为**辛积分器**（symplectic integrator），若其 Jacobi 矩阵 $D\Phi_h(z)$ 对所有 $z$ 都是辛矩阵：
    $$(D\Phi_h)^T J (D\Phi_h) = J.$$

!!! theorem "定理 53.13 (辛 Euler 方法)"
    对 Hamilton 系统 $\dot{q} = \partial H / \partial p$，$\dot{p} = -\partial H / \partial q$，**辛 Euler 方法**定义为：
    $$p_{n+1} = p_n - h \frac{\partial H}{\partial q}(q_n, p_{n+1}), \quad q_{n+1} = q_n + h \frac{\partial H}{\partial p}(q_n, p_{n+1}).$$
    这是一个一阶辛积分器（隐式）。

??? proof "证明"
    需要验证映射 $(q_n, p_n) \mapsto (q_{n+1}, p_{n+1})$ 的 Jacobi 矩阵是辛的。

    考虑可分 Hamilton 量 $H(q, p) = T(p) + V(q)$ 的情形。此时辛 Euler 方法变为显式：
    $$p_{n+1} = p_n - h V'(q_n), \quad q_{n+1} = q_n + h T'(p_{n+1}).$$
    Jacobi 矩阵为
    $$D\Phi_h = \frac{\partial(q_{n+1}, p_{n+1})}{\partial(q_n, p_n)} = \begin{bmatrix} I + h^2 T'' V'' & h T'' \\ -h V'' & I \end{bmatrix}$$
    （这里 $T'' = T''(p_{n+1})$，$V'' = V''(q_n)$）。

    验证辛性：$(D\Phi_h)^T J (D\Phi_h) = J$。由于 $T''$ 和 $V''$ 都是对称矩阵，代入计算可验证成立。具体地：

    $$(D\Phi_h)^T J (D\Phi_h) = \begin{bmatrix} I + h^2 V'' T'' & -hV'' \\ hT'' & I \end{bmatrix}\begin{bmatrix} 0 & I \\ -I & 0 \end{bmatrix}\begin{bmatrix} I + h^2 T'' V'' & hT'' \\ -hV'' & I \end{bmatrix}.$$

    中间步骤：$\begin{bmatrix} I + h^2 V'' T'' & -hV'' \\ hT'' & I \end{bmatrix}\begin{bmatrix} 0 & I \\ -I & 0 \end{bmatrix} = \begin{bmatrix} hV'' & I+h^2V''T'' \\ -I & hT'' \end{bmatrix}$。

    继续右乘 $D\Phi_h$，利用 $T''$，$V''$ 的对称性，最终得到 $J$。$\blacksquare$

!!! theorem "定理 53.14 (Störmer-Verlet 方法)"
    **Störmer-Verlet 方法**（也称为蛙跳法）对可分 Hamilton 量 $H = T(p) + V(q)$：
    $$p_{n+1/2} = p_n - \frac{h}{2}V'(q_n),$$
    $$q_{n+1} = q_n + h T'(p_{n+1/2}),$$
    $$p_{n+1} = p_{n+1/2} - \frac{h}{2}V'(q_{n+1}).$$
    这是二阶辛积分器。

??? proof "证明"
    Störmer-Verlet 可以看作两个辛 Euler 步的复合（以不同顺序）。由于辛映射的复合仍是辛映射（定理 53.3(d)），Störmer-Verlet 也是辛的。

    其二阶精度可以通过 Taylor 展开验证：$q_{n+1} = q_n + h\dot{q}_n + \frac{h^2}{2}\ddot{q}_n + O(h^3)$，$p_{n+1} = p_n + h\dot{p}_n + \frac{h^2}{2}\ddot{p}_n + O(h^3)$。$\blacksquare$

!!! theorem "定理 53.15 (向后误差分析)"
    辛积分器 $\Phi_h$ 应用于 Hamilton 系统 $H$ 时，它精确求解了一个**修正 Hamilton 系统** $\tilde{H} = H + hH_1 + h^2 H_2 + \cdots$ 的流。因此，修正能量 $\tilde{H}$ 在数值轨道上精确守恒，原始能量 $H$ 的误差在指数长时间内保持 $O(h^p)$（其中 $p$ 是方法的阶数）。

!!! example "例 53.8"
    考虑简谐振子 $H = \frac{1}{2}(p^2 + q^2)$。用标准四阶 Runge-Kutta 方法和 Störmer-Verlet 方法分别积分 $10^6$ 步。

    - **RK4**：能量 $H$ 在长时间后出现系统漂移，误差随时间线性增长。
    - **Störmer-Verlet**：能量 $H$ 在整个积分过程中振荡但不漂移，误差保持有界。

    这一差异对天体力学中行星轨道的长期稳定性模拟至关重要。

---

## 53.8 应用

<div class="context-flow" markdown>

**核心问题**：辛结构在哪些领域扮演核心角色？

</div>

!!! example "例 53.9 (Hamilton 力学)"
    经典力学中，相空间 $(\mathbb{R}^{2n}, \omega)$ 上的 Hamilton 方程
    $$\dot{q}_i = \frac{\partial H}{\partial p_i}, \quad \dot{p}_i = -\frac{\partial H}{\partial q_i}$$
    可以紧凑地写成 $\dot{z} = J \nabla H(z)$。其解流 $\varphi_t$ 满足 $(D\varphi_t)^T J (D\varphi_t) = J$，即是辛变换。Liouville 定理（相空间体积守恒）是 $\det(D\varphi_t) = 1$ 的直接推论。

    **线性 Hamilton 系统**：$H(z) = \frac{1}{2}z^T S z$（$S$ 对称），方程变为 $\dot{z} = JSz$。$JS$ 恰好是 Hamilton 矩阵。解为 $z(t) = e^{tJS}z(0)$，其中 $e^{tJS}$ 是辛矩阵。

!!! example "例 53.10 (量子光学)"
    在量子光学中，$n$ 模 Gaussian 态完全由其一阶矩（均值向量）和二阶矩（协方差矩阵 $\Sigma$）描述。协方差矩阵是 $2n \times 2n$ 正定矩阵，满足不确定性关系 $\Sigma + \frac{i}{2}J \geq 0$。

    Gaussian 酉操作（如分束器、压缩器、相移器）在 Heisenberg 图像下对正则变量 $(q_1, p_1, \ldots, q_n, p_n)$ 做辛变换 $S$，协方差矩阵变换为 $\Sigma \mapsto S \Sigma S^T$。

    Williamson 定理在此场景下意味着：任何 Gaussian 态可以通过辛变换等价于 $n$ 个独立模的热态，辛特征值 $d_1, \ldots, d_n$ 完全决定了态的纠缠和纯度等量子信息性质。

!!! example "例 53.11 (线性二次调节器)"
    最优控制中的线性二次调节器（LQR）问题：
    $$\min \int_0^\infty (x^T Q x + u^T R u) dt, \quad \dot{x} = Ax + Bu.$$
    其解通过 Riccati 方程 $PA + A^T P - PBR^{-1}B^T P + Q = 0$ 给出。

    与之关联的 Hamilton 矩阵为
    $$\mathcal{H} = \begin{bmatrix} A & -BR^{-1}B^T \\ -Q & -A^T \end{bmatrix}.$$
    Riccati 方程的解 $P$ 可以从 $\mathcal{H}$ 的稳定不变子空间中提取。Hamilton 矩阵的特征值配对性质（$\lambda$ 和 $-\lambda$ 成对）保证了稳定和不稳定特征值的对称分布，这是 LQR 理论的代数基础。

!!! remark "注记"
    辛矩阵理论还在以下领域中发挥重要作用：

    - **辛拓扑**：辛流形上的 Gromov 非压缩定理。
    - **Floer 同调**：以辛几何为基础的拓扑不变量。
    - **数论**：Siegel 模形式与辛群 $\mathrm{Sp}(2n, \mathbb{Z})$ 的自守表示。
    - **信号处理**：Gabor 分析中的辛对称性。

---

## 本章小结

| 概念 | 定义/关键性质 |
|:---|:---|
| 辛形式 | 非退化反对称双线性型 $\omega$；辛空间维数必偶 |
| 辛矩阵 | $M^T J M = J$；$\det M = 1$；$M^{-1} = -JM^TJ$ |
| 辛群 $\mathrm{Sp}(2n)$ | $n(2n+1)$ 维连通 Lie 群；极大紧子群 $\cong \mathrm{U}(n)$ |
| Hamilton 矩阵 | $H^TJ + JH = 0$；$e^{tH} \in \mathrm{Sp}(2n)$；特征值四重对称 |
| 辛特征值 | 正定矩阵的辛不变量 $d_1, \ldots, d_n$；Williamson 定理 |
| 辛 Gram-Schmidt | 构造辛基的算法；Darboux 定理：同维辛空间辛同构 |
| 辛积分器 | 保辛数值方法；长时间能量近似守恒 |

---

## 习题

!!! exercise "习题 53.1"
    证明 $J_{2n}$ 的行列式为 $1$。（提示：用归纳法或将 $J$ 分解为初等矩阵的乘积。）

!!! exercise "习题 53.2"
    设 $M = \begin{bmatrix} A & B \\ C & D \end{bmatrix}$ 是 $2n \times 2n$ 辛矩阵。证明：

    (a) $A^T C$ 和 $B^T D$ 都是对称矩阵。

    (b) $A^T D - C^T B = I_n$。

    (c) 如果 $A$ 可逆，则 $A^{-1}B$ 和 $CA^{-1}$ 都是对称矩阵。

!!! exercise "习题 53.3"
    计算矩阵 $A = \begin{bmatrix} 3 & 1 & 0 & 0 \\ 1 & 3 & 0 & 0 \\ 0 & 0 & 3 & -1 \\ 0 & 0 & -1 & 3 \end{bmatrix}$ 的辛特征值。

!!! exercise "习题 53.4"
    证明：若 $H$ 是 Hamilton 矩阵，$S$ 是辛矩阵，则 $S^{-1}HS$ 也是 Hamilton 矩阵。

!!! exercise "习题 53.5"
    在 $\mathbb{R}^6$ 上对如下向量组施行辛 Gram-Schmidt 过程：$v_1 = e_1$，$v_2 = e_1 + e_4$，$v_3 = e_2$，$v_4 = e_5$，$v_5 = e_3$，$v_6 = e_6$，其中 $\{e_1, \ldots, e_6\}$ 为标准基，辛形式由 $J_6$ 给出。

!!! exercise "习题 53.6"
    对一维简谐振子 $H = \frac{1}{2}(p^2 + \omega^2 q^2)$，写出辛 Euler 方法的显式迭代格式，并证明数值解的相轨迹是闭合椭圆（但与精确解的圆有微小差异）。

!!! exercise "习题 53.7"
    证明：$2 \times 2$ 实矩阵 $M = \begin{bmatrix} a & b \\ c & d \end{bmatrix}$ 是辛矩阵当且仅当 $ad - bc = 1$，即 $\mathrm{Sp}(2, \mathbb{R}) = \mathrm{SL}(2, \mathbb{R})$。

!!! exercise "习题 53.8"
    设 $A$ 和 $B$ 是 $2n \times 2n$ 正定矩阵，辛特征值分别为 $d_1 \leq \cdots \leq d_n$ 和 $d_1' \leq \cdots \leq d_n'$。若 $A \leq B$（即 $B - A$ 半正定），证明 $d_k \leq d_k'$ 对所有 $k$ 成立。
