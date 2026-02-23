# 第 13B 章 λ-矩阵与有理标准形

<div class="context-flow" markdown>

**前置**：Ch5 特征值与特征多项式 · Ch7 Jordan 标准形

**本章脉络**：$\lambda$-矩阵 → 初等变换 → Smith 标准形 → 不变因子/初等因子 → 友矩阵与**有理标准形** → Jordan 形的统一推导 → 矩阵相似的完全不变量

**延伸**：有理标准形在控制理论中对应系统的可观测标准形（companion matrix 实现）；Smith 标准形用于整数矩阵的分类（Abel 群的结构定理）和代数 K-理论

</div>

Jordan 标准形是矩阵理论中最精细的相似分类工具，但它本质上依赖于特征多项式在基础域上的完全分解。当基础域不是代数闭域时（例如 $\mathbb{R}$ 或 $\mathbb{Q}$），Jordan 标准形可能不存在。$\lambda$-矩阵理论提供了一套不依赖特征值分解的工具，由此引出的**有理标准形**（rational canonical form）在任意域上都成立，是矩阵相似问题的终极回答。

---

## 13B.1 λ-矩阵

<div class="context-flow" markdown>

**数值矩阵** $A$ → **多项式矩阵** $A(\lambda)$：元素从 $\mathbb{F}$ 扩展到 $\mathbb{F}[\lambda]$ → 特征矩阵 $\lambda I - A$ 是最重要的 $\lambda$-矩阵

</div>

!!! definition "定义 13B.1 (λ-矩阵)"
    设 $\mathbb{F}$ 是域。以 $\mathbb{F}[\lambda]$（$\mathbb{F}$ 上一元多项式环）中元素为元素的 $m \times n$ 矩阵

    $$
    A(\lambda) = (a_{ij}(\lambda))_{m \times n}
    $$

    称为 **$\lambda$-矩阵**（polynomial matrix / $\lambda$-matrix）。若 $A(\lambda)$ 中元素的最高次数为 $s$，则可写为

    $$
    A(\lambda) = A_s\lambda^s + A_{s-1}\lambda^{s-1} + \cdots + A_1\lambda + A_0
    $$

    其中 $A_i \in \mathbb{F}^{m \times n}$ 为数值矩阵系数。

!!! definition "定义 13B.2 (特征矩阵)"
    设 $A \in \mathbb{F}^{n \times n}$。矩阵

    $$
    \lambda I - A
    $$

    称为 $A$ 的**特征矩阵**（characteristic matrix）。它是最重要的 $\lambda$-矩阵，其行列式 $\det(\lambda I - A)$ 即特征多项式。

!!! definition "定义 13B.3 (λ-矩阵的秩)"
    $\lambda$-矩阵 $A(\lambda)$ 的**秩**定义为其非零子式的最高阶数，记作 $\operatorname{rank} A(\lambda)$。若 $A(\lambda)$ 是 $n \times n$ 方阵且 $\operatorname{rank} A(\lambda) = n$（即 $\det A(\lambda) \not\equiv 0$），则称 $A(\lambda)$ 是**非奇异的**。

!!! theorem "定理 13B.1 (λ-矩阵的可逆性)"
    $n$ 阶 $\lambda$-矩阵 $A(\lambda)$ 可逆（即存在 $\lambda$-矩阵 $B(\lambda)$ 使 $A(\lambda)B(\lambda) = I$）当且仅当 $\det A(\lambda)$ 是非零常数。

??? proof "证明"
    充分性：若 $\det A(\lambda) = c \neq 0$（常数），则 $A(\lambda)$ 的伴随矩阵 $\operatorname{adj}(A(\lambda))$ 也是 $\lambda$-矩阵，$A(\lambda)^{-1} = \frac{1}{c}\operatorname{adj}(A(\lambda))$。

    必要性：若 $A(\lambda)B(\lambda) = I$，取行列式得 $\det A(\lambda) \cdot \det B(\lambda) = 1$。由于 $\det A(\lambda)$ 和 $\det B(\lambda)$ 都是 $\mathbb{F}[\lambda]$ 中的多项式，乘积为常数 $1$，故两者都必须是非零常数。$\blacksquare$

!!! example "例 13B.1"
    $\lambda$-矩阵 $A(\lambda) = \begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$ 是非奇异的，$\det A(\lambda) = \lambda^2 \neq 0$。但 $A(\lambda)$ 不可逆，因为 $\lambda^2$ 不是常数。

    $\lambda$-矩阵 $B(\lambda) = \begin{pmatrix} 1 & \lambda \\ 0 & 1 \end{pmatrix}$ 可逆，$\det B(\lambda) = 1$，$B(\lambda)^{-1} = \begin{pmatrix} 1 & -\lambda \\ 0 & 1 \end{pmatrix}$。

!!! example "例 13B.2"
    矩阵 $A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}$ 的特征矩阵为

    $$
    \lambda I - A = \begin{pmatrix} \lambda - 2 & -1 \\ 0 & \lambda - 3 \end{pmatrix}
    $$

    $\det(\lambda I - A) = (\lambda - 2)(\lambda - 3) = \lambda^2 - 5\lambda + 6$。

---

## 13B.2 λ-矩阵的初等变换

<div class="context-flow" markdown>

数值矩阵的初等行/列变换推广到 $\lambda$-矩阵 → 注意倍乘因子必须是**非零常数**（不是多项式）→ 三类初等变换保持**等价性**

</div>

!!! definition "定义 13B.4 (λ-矩阵的初等变换)"
    $\lambda$-矩阵的**初等变换**（elementary operations）包括以下三类：

    1. **交换**两行（列）：$r_i \leftrightarrow r_j$；
    2. **用非零常数** $c \in \mathbb{F} \setminus \{0\}$ **乘**某行（列）：$r_i \to c \cdot r_i$；
    3. **将某行（列）的多项式倍**加到另一行（列）：$r_i \to r_i + p(\lambda) r_j$（$i \neq j$，$p(\lambda) \in \mathbb{F}[\lambda]$）。

    注意第 2 类中的乘数必须是非零**常数**，而非任意多项式。

!!! definition "定义 13B.5 (等价)"
    若 $\lambda$-矩阵 $A(\lambda)$ 经过有限次初等变换可化为 $B(\lambda)$，则称 $A(\lambda)$ 与 $B(\lambda)$ **等价**（equivalent），记作 $A(\lambda) \sim B(\lambda)$。等价地，存在可逆 $\lambda$-矩阵 $P(\lambda), Q(\lambda)$ 使得

    $$
    B(\lambda) = P(\lambda) A(\lambda) Q(\lambda)
    $$

!!! theorem "定理 13B.2 (等价是等价关系)"
    $\lambda$-矩阵的等价关系满足自反性、对称性和传递性。

??? proof "证明"
    自反性：取 $P = Q = I$。对称性：若 $B = PAQ$，则 $A = P^{-1}BQ^{-1}$，$P^{-1}$ 和 $Q^{-1}$ 都是可逆 $\lambda$-矩阵。传递性：若 $B = P_1AQ_1$，$C = P_2BQ_2$，则 $C = (P_2P_1)A(Q_1Q_2)$。$\blacksquare$

!!! example "例 13B.3"
    对 $\lambda$-矩阵 $A(\lambda) = \begin{pmatrix} \lambda - 2 & -1 \\ 0 & \lambda - 3 \end{pmatrix}$ 进行初等变换。

    $r_1 \to r_1 + \frac{1}{\lambda-3} r_2$？不行——$\frac{1}{\lambda-3}$ 不是多项式。

    正确的做法是利用多项式除法。例如 $c_2 \to c_2 + c_1$：

    $$
    \begin{pmatrix} \lambda - 2 & \lambda - 3 \\ 0 & \lambda - 3 \end{pmatrix}
    $$

    再 $c_1 \to c_1 - c_2$：$\begin{pmatrix} 1 & \lambda - 3 \\ -(\lambda-3) & \lambda - 3 \end{pmatrix}$。继续变换可化为 Smith 标准形。

---

## 13B.3 Smith 标准形

<div class="context-flow" markdown>

每个 $\lambda$-矩阵都等价于唯一的**对角形** → 对角元素满足**整除链** $d_1 \mid d_2 \mid \cdots$ → 这是 $\lambda$-矩阵理论的核心定理

</div>

!!! definition "定义 13B.6 (Smith 标准形)"
    $m \times n$（$m \leq n$） $\lambda$-矩阵 $A(\lambda)$（秩为 $r$）的 **Smith 标准形**（Smith normal form）是如下形式的对角矩阵：

    $$
    S(\lambda) = \begin{pmatrix} d_1(\lambda) & & & & \\ & d_2(\lambda) & & & \\ & & \ddots & & \\ & & & d_r(\lambda) & \\ & & & & O \end{pmatrix}
    $$

    其中 $d_i(\lambda)$ 是首一多项式（首项系数为 1），且满足**整除链**：

    $$
    d_1(\lambda) \mid d_2(\lambda) \mid \cdots \mid d_r(\lambda)
    $$

!!! theorem "定理 13B.3 (Smith 标准形的存在性与唯一性)"
    每个非零 $\lambda$-矩阵 $A(\lambda)$ 都等价于唯一的 Smith 标准形。

??? proof "证明"
    **存在性**：通过初等变换构造性地化简。

    **步骤 1**：通过行列交换，将 $A(\lambda)$ 中次数最低的非零元素移到 $(1,1)$ 位置。

    **步骤 2**：若 $(1,1)$ 元素不能整除第一行或第一列的某个元素，用带余除法和初等变换降低 $(1,1)$ 元素的次数。重复此过程直到 $(1,1)$ 元素整除第一行和第一列的所有元素。

    **步骤 3**：用 $(1,1)$ 元素消去第一行和第一列的其他元素。若 $(1,1)$ 元素不能整除某个非第一行第一列的元素，将该元素加到第一行，回到步骤 2。

    **步骤 4**：最终 $(1,1)$ 元素整除所有其他元素。设为 $d_1(\lambda)$（首一化），对右下角的子矩阵递归进行。

    **唯一性**：由行列式因子的唯一性（定理 13B.8）保证。$\blacksquare$

!!! theorem "定理 13B.4 (初等变换保持等价)"
    初等变换不改变 $\lambda$-矩阵的 Smith 标准形。等价地，$A(\lambda) \sim B(\lambda)$ 当且仅当它们有相同的 Smith 标准形。

??? proof "证明"
    每个初等变换都可以表示为左乘或右乘一个可逆 $\lambda$-矩阵（初等矩阵），而可逆 $\lambda$-矩阵的行列式为非零常数，不改变各阶行列式因子（见 13B.5 节），从而 Smith 标准形不变。$\blacksquare$

!!! example "例 13B.4"
    求 $A(\lambda) = \begin{pmatrix} \lambda & \lambda^2 \\ 1 & \lambda \end{pmatrix}$ 的 Smith 标准形。

    $r_1 \leftrightarrow r_2$：$\begin{pmatrix} 1 & \lambda \\ \lambda & \lambda^2 \end{pmatrix}$。

    $r_2 \to r_2 - \lambda r_1$：$\begin{pmatrix} 1 & \lambda \\ 0 & 0 \end{pmatrix}$。

    $c_2 \to c_2 - \lambda c_1$：$\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。

    Smith 标准形为 $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，秩为 $1$。

!!! example "例 13B.5"
    求 $A(\lambda) = \begin{pmatrix} \lambda - 1 & 0 \\ 0 & (\lambda-1)^2 \end{pmatrix}$ 的 Smith 标准形。

    已经是对角形，且 $(\lambda-1) \mid (\lambda-1)^2$，故 Smith 标准形就是

    $$
    S(\lambda) = \begin{pmatrix} \lambda - 1 & 0 \\ 0 & (\lambda-1)^2 \end{pmatrix}
    $$

    $d_1(\lambda) = \lambda - 1$，$d_2(\lambda) = (\lambda-1)^2$。

---

## 13B.4 不变因子与初等因子

<div class="context-flow" markdown>

Smith 标准形的对角元 $d_1, \ldots, d_r$ = **不变因子** → 每个不变因子分解为不可约因式的幂 = **初等因子** → 初等因子完全决定 Jordan 块的结构

</div>

!!! definition "定义 13B.7 (不变因子)"
    $\lambda$-矩阵 $A(\lambda)$（秩 $r$）的 Smith 标准形中的对角元 $d_1(\lambda), d_2(\lambda), \ldots, d_r(\lambda)$ 称为 $A(\lambda)$ 的**不变因子**（invariant factors）。不变因子满足整除链 $d_1 \mid d_2 \mid \cdots \mid d_r$。

!!! definition "定义 13B.8 (初等因子)"
    将每个不变因子 $d_k(\lambda)$ 在 $\mathbb{F}[\lambda]$ 中分解为不可约多项式的幂次之积：

    $$
    d_k(\lambda) = p_1(\lambda)^{e_{k1}} p_2(\lambda)^{e_{k2}} \cdots p_s(\lambda)^{e_{ks}}
    $$

    其中所有 $e_{ki} \geq 0$。所有满足 $e_{ki} \geq 1$ 的因式 $p_i(\lambda)^{e_{ki}}$ 合在一起（按不可约因式分组），称为 $A(\lambda)$ 的**初等因子**（elementary divisors）。

!!! theorem "定理 13B.5 (不变因子与初等因子的关系)"
    设 $A(\lambda)$ 的不变因子为 $d_1(\lambda), \ldots, d_r(\lambda)$，设 $\mathbb{F}[\lambda]$ 中涉及的不可约多项式为 $p_1(\lambda), \ldots, p_s(\lambda)$。则

    $$
    d_k(\lambda) = \prod_{j=1}^s p_j(\lambda)^{e_{kj}}
    $$

    其中 $0 \leq e_{1j} \leq e_{2j} \leq \cdots \leq e_{rj}$（整除链条件）。

    反之，给定初等因子组，可以唯一地恢复出不变因子：对每个不可约多项式 $p_j$，将其幂次 $e_{1j} \leq \cdots \leq e_{rj}$ 从最后一个不变因子向前分配。

??? proof "证明"
    由整除链 $d_1 \mid d_2 \mid \cdots \mid d_r$ 直接得到 $e_{1j} \leq e_{2j} \leq \cdots \leq e_{rj}$。

    反向恢复：最后一个不变因子 $d_r = \prod_j p_j^{e_{rj}}$（取每个不可约因式的最高次幂）。$d_{r-1} = \prod_j p_j^{e_{r-1,j}}$（取次高次幂），以此类推。唯一性由分解的唯一性保证。$\blacksquare$

!!! theorem "定理 13B.6 (特征矩阵的不变因子与特征多项式)"
    设 $A \in \mathbb{F}^{n \times n}$，$\lambda I - A$ 的不变因子为 $d_1(\lambda), \ldots, d_n(\lambda)$。则

    $$
    d_1(\lambda) d_2(\lambda) \cdots d_n(\lambda) = \det(\lambda I - A) = p_A(\lambda)
    $$

    即不变因子之积等于特征多项式。最后一个不变因子 $d_n(\lambda)$ 等于最小多项式 $m_A(\lambda)$。

??? proof "证明"
    Smith 标准形 $S(\lambda)$ 与 $\lambda I - A$ 等价，$\det S(\lambda) = c \cdot \det(\lambda I - A)$（$c$ 为非零常数）。由于 $S(\lambda) = \operatorname{diag}(d_1, \ldots, d_n)$ 且所有 $d_i$ 首一，$c = 1$。

    对于最小多项式，$d_n(\lambda) = m_A(\lambda)$ 的证明需要利用不变因子与零化多项式的关系：$d_n(\lambda)$ 整除所有零化 $A$ 的多项式，且 $d_n(A) = 0$（由 Cayley-Hamilton 定理的推广）。$\blacksquare$

!!! example "例 13B.6"
    设 $A = \begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}$。特征矩阵

    $$
    \lambda I - A = \begin{pmatrix} \lambda & 0 & 0 \\ -1 & \lambda & 0 \\ 0 & -1 & \lambda \end{pmatrix}
    $$

    初等变换化简：$c_1 \to c_1 + \lambda c_2$，$c_2 \to c_2 + \lambda c_3$：

    $$
    \begin{pmatrix} \lambda^2 - 1 & \lambda^2 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & \lambda \end{pmatrix} \to \cdots \to \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & \lambda^3 \end{pmatrix}
    $$

    不变因子：$d_1 = 1$，$d_2 = 1$，$d_3 = \lambda^3$。特征多项式 $p_A(\lambda) = \lambda^3$，最小多项式 $m_A(\lambda) = \lambda^3$。

    在 $\mathbb{C}$ 上的初等因子为 $\lambda^3$（唯一一个），对应 $3 \times 3$ Jordan 块 $J_3(0)$。

!!! example "例 13B.7"
    设矩阵 $A$ 的特征矩阵 $\lambda I - A$ 的 Smith 标准形为

    $$
    \operatorname{diag}(1, 1, \lambda - 1, (\lambda-1)(\lambda-2)^2)
    $$

    不变因子：$d_1 = d_2 = 1$，$d_3 = \lambda - 1$，$d_4 = (\lambda-1)(\lambda-2)^2$。

    特征多项式：$p_A(\lambda) = (\lambda-1)^2(\lambda-2)^2$。

    最小多项式：$m_A(\lambda) = (\lambda-1)(\lambda-2)^2$。

    初等因子（在 $\mathbb{Q}$ 或 $\mathbb{R}$ 或 $\mathbb{C}$ 上）：$(\lambda-1)^1$（来自 $d_3$），$(\lambda-1)^1$（来自 $d_4$），$(\lambda-2)^2$（来自 $d_4$）。

---

## 13B.5 行列式因子

<div class="context-flow" markdown>

行列式因子 $D_k(\lambda)$ = 所有 $k$ 阶子式的最大公因式 → $d_k = D_k/D_{k-1}$ → 行列式因子提供了**计算不变因子的实用方法**，无需做初等变换

</div>

!!! definition "定义 13B.9 (行列式因子)"
    设 $A(\lambda)$ 是秩为 $r$ 的 $\lambda$-矩阵。对 $k = 1, 2, \ldots, r$，$A(\lambda)$ 的所有 $k$ 阶子式的首一最大公因式 $D_k(\lambda)$ 称为 $A(\lambda)$ 的 **$k$ 阶行列式因子**（$k$-th determinantal divisor）。约定 $D_0(\lambda) = 1$。

!!! theorem "定理 13B.7 (行列式因子的整除性)"
    行列式因子满足 $D_k(\lambda) \mid D_{k+1}(\lambda)$（$k = 1, \ldots, r-1$）。更精确地说，$D_k(\lambda)$ 整除 $A(\lambda)$ 的每个 $k$ 阶子式。

??? proof "证明"
    每个 $(k+1)$ 阶子式按某行展开，是若干 $k$ 阶子式的 $\mathbb{F}[\lambda]$-线性组合。因此 $D_k(\lambda)$ 整除每个 $(k+1)$ 阶子式，从而 $D_k(\lambda) \mid D_{k+1}(\lambda)$。$\blacksquare$

!!! theorem "定理 13B.8 (不变因子与行列式因子的关系)"
    $\lambda$-矩阵 $A(\lambda)$ 的不变因子 $d_k(\lambda)$ 与行列式因子 $D_k(\lambda)$ 的关系为

    $$
    d_k(\lambda) = \frac{D_k(\lambda)}{D_{k-1}(\lambda)}, \quad k = 1, 2, \ldots, r
    $$

    其中 $D_0(\lambda) = 1$。

??? proof "证明"
    设 $A(\lambda)$ 的 Smith 标准形为 $S(\lambda) = \operatorname{diag}(d_1, \ldots, d_r, 0, \ldots, 0)$。由于初等变换不改变行列式因子（可逆 $\lambda$-矩阵的行列式为非零常数，Binet-Cauchy 公式保证子式的 GCD 不变），$A(\lambda)$ 与 $S(\lambda)$ 有相同的行列式因子。

    $S(\lambda)$ 的 $k$ 阶子式中非零的仅有 $\prod_{i \in I} d_i(\lambda)$（$|I| = k$，$I \subseteq \{1, \ldots, r\}$）。由整除链 $d_1 \mid d_2 \mid \cdots$，最大公因式为 $d_1 d_2 \cdots d_k$。故

    $$
    D_k(\lambda) = d_1(\lambda)d_2(\lambda)\cdots d_k(\lambda)
    $$

    从而 $d_k(\lambda) = D_k/D_{k-1}$。$\blacksquare$

!!! corollary "推论 13B.1"
    等价的 $\lambda$-矩阵有相同的行列式因子、不变因子和初等因子。反之，相同的不变因子（或等价地，相同的行列式因子）保证 $\lambda$-矩阵等价。

!!! example "例 13B.8"
    用行列式因子求矩阵 $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ 的不变因子。

    特征矩阵 $\lambda I - A = \begin{pmatrix} \lambda - 2 & -1 \\ -1 & \lambda - 2 \end{pmatrix}$。

    $D_1(\lambda)$：所有 $1$ 阶子式为 $\lambda - 2, -1, -1, \lambda - 2$。$\gcd = 1$。

    $D_2(\lambda)$：$\det(\lambda I - A) = (\lambda-2)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda-1)(\lambda-3)$。

    不变因子：$d_1 = D_1/D_0 = 1$，$d_2 = D_2/D_1 = (\lambda-1)(\lambda-3)$。

!!! example "例 13B.9"
    求 $A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$ 的行列式因子和不变因子。

    $\lambda I - A = \operatorname{diag}(\lambda-1, \lambda-1, \lambda-2)$。

    $D_1 = \gcd(\lambda-1, \lambda-1, \lambda-2) = 1$。

    $D_2$：$2$ 阶子式为 $(\lambda-1)^2, (\lambda-1)(\lambda-2), (\lambda-1)(\lambda-2)$。$D_2 = \gcd = \lambda - 1$。

    $D_3 = (\lambda-1)^2(\lambda-2)$。

    不变因子：$d_1 = 1$，$d_2 = \lambda - 1$，$d_3 = (\lambda-1)(\lambda-2)$。

---

## 13B.6 友矩阵与有理标准形

<div class="context-flow" markdown>

$n$ 次多项式 $d(\lambda)$ $\to$ **友矩阵** $C(d)$：其特征多项式和最小多项式都等于 $d(\lambda)$ → 矩阵 $A$ 相似于友矩阵的直和 $\Leftrightarrow$ 有理标准形 → **任意域上成立**

</div>

!!! definition "定义 13B.10 (友矩阵)"
    设 $d(\lambda) = \lambda^m + a_{m-1}\lambda^{m-1} + \cdots + a_1\lambda + a_0$ 是 $\mathbb{F}[\lambda]$ 中的首一多项式。其**友矩阵**（companion matrix）定义为

    $$
    C(d) = \begin{pmatrix} 0 & 0 & \cdots & 0 & -a_0 \\ 1 & 0 & \cdots & 0 & -a_1 \\ 0 & 1 & \cdots & 0 & -a_2 \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & 1 & -a_{m-1} \end{pmatrix}_{m \times m}
    $$

!!! theorem "定理 13B.9 (友矩阵的性质)"
    设 $d(\lambda) = \lambda^m + a_{m-1}\lambda^{m-1} + \cdots + a_0$ 是首一多项式，$C = C(d)$ 是其友矩阵。则

    1. $C$ 的特征多项式为 $p_C(\lambda) = d(\lambda)$；
    2. $C$ 的最小多项式为 $m_C(\lambda) = d(\lambda)$；
    3. $\lambda I - C$ 的 Smith 标准形为 $\operatorname{diag}(1, 1, \ldots, 1, d(\lambda))$。

??? proof "证明"
    **(1)** 沿最后一列展开 $\det(\lambda I - C)$：

    $$
    \det(\lambda I - C) = \det\begin{pmatrix} \lambda & 0 & \cdots & 0 & a_0 \\ -1 & \lambda & \cdots & 0 & a_1 \\ 0 & -1 & \cdots & 0 & a_2 \\ \vdots & & \ddots & & \vdots \\ 0 & 0 & \cdots & -1 & \lambda + a_{m-1} \end{pmatrix}
    $$

    沿第一列展开并递推，得 $\det(\lambda I - C) = \lambda^m + a_{m-1}\lambda^{m-1} + \cdots + a_0 = d(\lambda)$。

    **(2)** 设 $\mathbf{e}_1 = (1, 0, \ldots, 0)^T$。直接计算可得 $C\mathbf{e}_1 = \mathbf{e}_2$，$C^2\mathbf{e}_1 = \mathbf{e}_3$，...，$C^{m-1}\mathbf{e}_1 = \mathbf{e}_m$。因此 $\{\mathbf{e}_1, C\mathbf{e}_1, \ldots, C^{m-1}\mathbf{e}_1\}$ 构成 $\mathbb{F}^m$ 的基。

    若 $f(C) = 0$，$\deg f < m$，则 $f(C)\mathbf{e}_1 = f_0\mathbf{e}_1 + f_1C\mathbf{e}_1 + \cdots + f_{m-1}C^{m-1}\mathbf{e}_1 = \mathbf{0}$，由基的线性无关性得 $f \equiv 0$。故最小多项式的次数至少为 $m$，由 Cayley-Hamilton 定理得 $m_C = p_C = d$。

    **(3)** 由 (1)(2) 直接得出。$\blacksquare$

!!! theorem "定理 13B.10 (有理标准形)"
    设 $A \in \mathbb{F}^{n \times n}$，$\lambda I - A$ 的不变因子为 $d_1(\lambda), \ldots, d_n(\lambda)$，其中 $d_1 = \cdots = d_{n-s} = 1$，$\deg d_{n-s+1} \geq 1, \ldots, \deg d_n \geq 1$。令非平凡的不变因子为 $f_1 = d_{n-s+1}, \ldots, f_s = d_n$（$f_1 \mid f_2 \mid \cdots \mid f_s$）。则 $A$ 相似于

    $$
    R = \begin{pmatrix} C(f_1) & & \\ & C(f_2) & \\ & & \ddots & \\ & & & C(f_s) \end{pmatrix}
    $$

    此即 $A$ 的**有理标准形**（rational canonical form），也称 **Frobenius 标准形**。它在任意域 $\mathbb{F}$ 上都成立。

??? proof "证明"
    $\lambda I - A$ 与 $\lambda I - R$ 有相同的 Smith 标准形。事实上，$\lambda I - R = \operatorname{diag}(\lambda I - C(f_1), \ldots, \lambda I - C(f_s))$，而 $\lambda I - C(f_k)$ 的 Smith 标准形为 $\operatorname{diag}(1, \ldots, 1, f_k(\lambda))$（定理 13B.9）。组合后得到的 Smith 标准形恰好是 $\operatorname{diag}(1, \ldots, 1, f_1, \ldots, f_s)$，与 $\lambda I - A$ 的 Smith 标准形一致。

    由 Smith 标准形唯一确定等价类（推论 13B.1），$\lambda I - A \sim \lambda I - R$，即存在可逆 $\lambda$-矩阵 $P(\lambda), Q(\lambda)$ 使 $P(\lambda)(\lambda I - A)Q(\lambda) = \lambda I - R$。

    由此可以证明（需要更精细的论证，利用 $\lambda$-矩阵等价与数值矩阵相似的桥梁定理）$A$ 与 $R$ 相似。$\blacksquare$

!!! theorem "定理 13B.11 (有理标准形的唯一性)"
    有理标准形在相似意义下是唯一的。即若 $A$ 同时相似于两个有理标准形 $R_1$ 和 $R_2$，则 $R_1 = R_2$（友矩阵块的顺序一致）。

??? proof "证明"
    $A \sim R_1$ 和 $A \sim R_2$ 意味着 $\lambda I - R_1$ 和 $\lambda I - R_2$ 有相同的 Smith 标准形。由友矩阵的 Smith 标准形（定理 13B.9），非平凡不变因子就是构成友矩阵块的多项式。由 Smith 标准形的唯一性，这些多项式相同，从而 $R_1 = R_2$。$\blacksquare$

!!! example "例 13B.10"
    设 $4$ 阶矩阵 $A$ 的不变因子为 $d_1 = 1$，$d_2 = 1$，$d_3 = \lambda - 1$，$d_4 = (\lambda-1)(\lambda-2)^2$。

    非平凡不变因子：$f_1 = \lambda - 1$，$f_2 = (\lambda-1)(\lambda-2)^2 = \lambda^3 - 5\lambda^2 + 8\lambda - 4$。

    $C(f_1) = (1)$（$1 \times 1$ 矩阵）。

    $$
    C(f_2) = \begin{pmatrix} 0 & 0 & 4 \\ 1 & 0 & -8 \\ 0 & 1 & 5 \end{pmatrix}
    $$

    有理标准形：

    $$
    R = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 4 \\ 0 & 1 & 0 & -8 \\ 0 & 0 & 1 & 5 \end{pmatrix}
    $$

!!! example "例 13B.11"
    求矩阵 $A = \begin{pmatrix} 0 & 1 \\ -6 & 5 \end{pmatrix}$ 的有理标准形。

    $p_A(\lambda) = \lambda^2 - 5\lambda + 6 = (\lambda-2)(\lambda-3)$。

    $D_1 = \gcd(\lambda, -1, 1, 5-\lambda) = 1$，$D_2 = (\lambda-2)(\lambda-3)$。

    不变因子：$d_1 = 1$，$d_2 = (\lambda-2)(\lambda-3)$。非平凡不变因子 $f_1 = \lambda^2 - 5\lambda + 6$。

    $$
    C(f_1) = \begin{pmatrix} 0 & -6 \\ 1 & 5 \end{pmatrix}
    $$

    因此 $A$ 相似于 $R = \begin{pmatrix} 0 & -6 \\ 1 & 5 \end{pmatrix}$。事实上取 $P = \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix}$，$P^{-1}AP = R$。

---

## 13B.7 有理标准形与 Jordan 标准形的关系

<div class="context-flow" markdown>

在代数闭域 $\mathbb{C}$ 上：不变因子完全分解 → 初等因子 $(\lambda - \lambda_i)^{e_i}$ → 友矩阵 $C((\lambda-\lambda_i)^{e_i})$ 相似于 Jordan 块 $J_{e_i}(\lambda_i)$ → 有理标准形化为 **Jordan 标准形**

</div>

!!! theorem "定理 13B.12 (友矩阵与 Jordan 块的关系)"
    在代数闭域 $\mathbb{F}$ 上，友矩阵 $C((\lambda - a)^m)$ 相似于 Jordan 块 $J_m(a)$。

??? proof "证明"
    $C((\lambda - a)^m)$ 的特征多项式和最小多项式都是 $(\lambda - a)^m$。而 $J_m(a)$ 的特征多项式和最小多项式也是 $(\lambda - a)^m$。两者的 Smith 标准形相同，都是 $\operatorname{diag}(1, \ldots, 1, (\lambda-a)^m)$。由有理标准形的唯一性，$C((\lambda-a)^m) \sim J_m(a)$。

    具体地，设 $C = C((\lambda-a)^m)$，$\mathbf{e}_1 = (1,0,\ldots,0)^T$。令 $\mathbf{v}_m = \mathbf{e}_1$，$\mathbf{v}_{m-1} = (C - aI)\mathbf{v}_m$，...，$\mathbf{v}_1 = (C-aI)^{m-1}\mathbf{v}_m$。则 $P = (\mathbf{v}_1, \ldots, \mathbf{v}_m)$ 满足 $P^{-1}CP = J_m(a)$。$\blacksquare$

!!! theorem "定理 13B.13 (从有理标准形到 Jordan 标准形)"
    设 $A \in \mathbb{C}^{n \times n}$，不变因子分解为初等因子后，$A$ 的有理标准形中每个友矩阵 $C(f_k)$ 可进一步相似变换为 Jordan 块的直和：

    若 $f_k(\lambda) = \prod_{i} (\lambda - \lambda_i)^{e_{ki}}$，则

    $$
    C(f_k) \sim \bigoplus_{i} J_{e_{ki}}(\lambda_i)
    $$

    因此 $A$ 的 Jordan 标准形为

    $$
    J = \bigoplus_{k,i} J_{e_{ki}}(\lambda_i)
    $$

??? proof "证明"
    对 $C(f_k)$，其特征多项式和最小多项式都是 $f_k(\lambda)$。由于 $f_k$ 可分解为互素因子之积 $f_k = \prod_i (\lambda - \lambda_i)^{e_{ki}}$，由初等因子理论（或直接用 $\mathbb{F}[\lambda]$-模的结构定理），$C(f_k)$ 相似于 $\bigoplus_i C((\lambda-\lambda_i)^{e_{ki}})$。

    由定理 13B.12，$C((\lambda-\lambda_i)^{e_{ki}}) \sim J_{e_{ki}}(\lambda_i)$。$\blacksquare$

!!! example "例 13B.12"
    接例 13B.10。$A$ 的不变因子为 $d_3 = \lambda-1$，$d_4 = (\lambda-1)(\lambda-2)^2$。

    初等因子：$(\lambda-1)^1$（来自 $d_3$），$(\lambda-1)^1$（来自 $d_4$），$(\lambda-2)^2$（来自 $d_4$）。

    有理标准形 → Jordan 标准形：

    $$
    R = \begin{pmatrix} C(\lambda-1) & \\ & C((\lambda-1)(\lambda-2)^2) \end{pmatrix} \sim \begin{pmatrix} 1 & & & \\ & 1 & & \\ & & 2 & 1 \\ & & 0 & 2 \end{pmatrix} = J
    $$

!!! example "例 13B.13"
    考虑 $A \in \mathbb{R}^{4 \times 4}$，不变因子为 $d_1 = d_2 = 1$，$d_3 = \lambda^2+1$，$d_4 = (\lambda^2+1)^2$。

    在 $\mathbb{R}$ 上，$\lambda^2 + 1$ 不可约，有理标准形为

    $$
    R = \begin{pmatrix} C(\lambda^2+1) & \\ & C((\lambda^2+1)^2) \end{pmatrix} = \begin{pmatrix} 0 & -1 & & \\ 1 & 0 & & \\ & & 0 & 0 & 0 & 1 \\ & & 1 & 0 & 0 & 0 \\ & & 0 & 1 & 0 & -2 \\ & & 0 & 0 & 1 & 0 \end{pmatrix}
    $$

    这是 $\mathbb{R}$ 上的"最简形式"——无法进一步简化，因为特征值 $\pm i$ 不在 $\mathbb{R}$ 中。

    在 $\mathbb{C}$ 上，$\lambda^2+1 = (\lambda-i)(\lambda+i)$，初等因子为 $(\lambda-i)^1, (\lambda+i)^1, (\lambda-i)^2, (\lambda+i)^2$，Jordan 标准形为 $\operatorname{diag}(i, -i, J_2(i), J_2(-i))$。

---

## 13B.8 矩阵相似的完全不变量

<div class="context-flow" markdown>

**终极判据**：$A \sim B$ $\Leftrightarrow$ $\lambda I - A$ 与 $\lambda I - B$ 有相同的 Smith 标准形 $\Leftrightarrow$ 相同的不变因子 $\Leftrightarrow$ 相同的初等因子组 → 在任意域上一劳永逸地解决矩阵相似问题

</div>

!!! definition "定义 13B.11 (完全不变量)"
    若矩阵的某组不变量满足：$A \sim B$ 当且仅当 $A$ 和 $B$ 具有相同的该组不变量，则称该组不变量为矩阵相似的**完全不变量**（complete system of invariants）。

!!! theorem "定理 13B.14 (矩阵相似的判定)"
    设 $A, B \in \mathbb{F}^{n \times n}$。以下条件等价：

    1. $A$ 与 $B$ 相似（存在可逆 $P$ 使 $B = P^{-1}AP$）；
    2. $\lambda I - A$ 与 $\lambda I - B$ 等价（作为 $\lambda$-矩阵）；
    3. $\lambda I - A$ 与 $\lambda I - B$ 有相同的 Smith 标准形；
    4. $\lambda I - A$ 与 $\lambda I - B$ 有相同的行列式因子；
    5. $\lambda I - A$ 与 $\lambda I - B$ 有相同的不变因子；
    6. $\lambda I - A$ 与 $\lambda I - B$ 有相同的初等因子组；
    7. $A$ 与 $B$ 有相同的有理标准形。

??? proof "证明"
    **(1)$\Rightarrow$(2)：** 若 $B = P^{-1}AP$（$P$ 为常数可逆矩阵），则

    $$
    \lambda I - B = \lambda I - P^{-1}AP = P^{-1}(\lambda I - A)P
    $$

    $P$ 和 $P^{-1}$ 都是可逆 $\lambda$-矩阵（行列式为非零常数），故 $\lambda I - A \sim \lambda I - B$。

    **(2)$\Leftrightarrow$(3)$\Leftrightarrow$(4)$\Leftrightarrow$(5)$\Leftrightarrow$(6)：** 由 Smith 标准形的唯一性和行列式因子/不变因子/初等因子的一一对应关系。

    **(2)$\Rightarrow$(1)：** 这是最非平凡的部分。需要证明 $\lambda$-矩阵的等价蕴含数值矩阵的相似。设 $P(\lambda)(\lambda I - A)Q(\lambda) = \lambda I - B$。将 $P(\lambda)$ 和 $Q(\lambda)$ 按 $(\lambda I - A)$ 做带余除法进行分析，可以构造出使 $P^{-1}AP = B$ 的常数矩阵 $P$。

    **(5)$\Leftrightarrow$(7)：** 有理标准形由不变因子唯一确定。$\blacksquare$

!!! theorem "定理 13B.15 (各类不变量的层次)"
    对于 $A \in \mathbb{F}^{n \times n}$，以下不变量的信息量递增：

    - **特征多项式** $p_A(\lambda)$：等于不变因子之积 $\prod d_k$。仅知特征多项式**不足以**判定相似；
    - **最小多项式** $m_A(\lambda)$：等于最后一个不变因子 $d_n$。特征多项式 + 最小多项式仍不足以判定相似；
    - **不变因子组** $\{d_1, \ldots, d_n\}$：**完全不变量**，完全判定相似。

??? proof "证明"
    反例说明前两项不充分：

    $A_1 = \operatorname{diag}(1, 1, 2)$ 和 $A_2 = \begin{pmatrix} 1 & 1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$ 有相同的特征多项式 $(\lambda-1)^2(\lambda-2)$，但不变因子不同：$A_1$ 的不变因子为 $1, \lambda-1, (\lambda-1)(\lambda-2)$；$A_2$ 的不变因子为 $1, 1, (\lambda-1)^2(\lambda-2)$。故 $A_1 \not\sim A_2$。

    两者有不同的最小多项式：$m_{A_1} = (\lambda-1)(\lambda-2)$，$m_{A_2} = (\lambda-1)^2(\lambda-2)$。但即使特征多项式和最小多项式都相同，不变因子仍可能不同（对 $n \geq 4$ 可构造反例）。$\blacksquare$

!!! theorem "定理 13B.16 (Cayley-Hamilton 定理的精化)"
    设 $A$ 的不变因子为 $d_1 \mid d_2 \mid \cdots \mid d_n$。则

    1. $d_n(A) = 0$（最小多项式零化 $A$）；
    2. 任意 $f(\lambda) \in \mathbb{F}[\lambda]$ 满足 $f(A) = 0$ 当且仅当 $d_n \mid f$；
    3. $p_A(\lambda) = d_1 \cdots d_n$ 且 $d_n \mid p_A$（Cayley-Hamilton 定理 $p_A(A) = 0$ 是 $d_n(A) = 0$ 的推论）。

??? proof "证明"
    (1) $A$ 相似于有理标准形 $R = \bigoplus C(f_k)$。由定理 13B.9，$f_k(C(f_k)) = 0$。由整除链 $f_k \mid f_s = d_n$，故 $d_n(C(f_k)) = 0$。因此 $d_n(R) = 0$，由相似性 $d_n(A) = 0$。

    (2) $f(A) = 0$ $\Leftrightarrow$ $f(R) = 0$ $\Leftrightarrow$ $f(C(f_k)) = 0$ 对所有 $k$ $\Leftrightarrow$ $f_k \mid f$ 对所有 $k$ $\Leftrightarrow$ $d_n \mid f$（因为 $d_n = f_s$ 是最大的不变因子）。$\blacksquare$

!!! example "例 13B.14"
    判断 $A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 2 \end{pmatrix}$ 和 $B = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 2 \end{pmatrix}$ 是否相似。

    $A$：$D_1 = 1$，$D_2 = \gcd((\lambda-1)(\lambda-2), (\lambda-1)(\lambda-2), (\lambda-2)^2) = \lambda-2$。$D_3 = (\lambda-1)(\lambda-2)^2$。

    不变因子：$d_1 = 1$，$d_2 = \lambda-2$，$d_3 = (\lambda-1)(\lambda-2)$。

    $B$：$D_1 = \gcd(\lambda-1, \lambda-2, -1, \ldots) = 1$。

    $\lambda I - B$ 的 $2$ 阶子式包括 $(\lambda-1)(\lambda-2)$，$-(\lambda-1)$，$(\lambda-2)^2 + 0 = (\lambda-2)^2$ 等。$D_2 = \gcd = 1$（因为 $-(\lambda-1)$ 和 $(\lambda-2)$ 互素）。

    等等，重新计算。$\lambda I - B = \begin{pmatrix} \lambda-1 & 0 & 0 \\ 0 & \lambda-2 & -1 \\ 0 & 0 & \lambda-2 \end{pmatrix}$。

    $2$ 阶子式：$(\lambda-1)(\lambda-2)$（取行12列12），$-(\lambda-1) \cdot 0 = 0$（取行13列12... 不对），逐个列出：行$\{1,2\}$列$\{1,2\}$: $(\lambda-1)(\lambda-2)$；行$\{1,2\}$列$\{1,3\}$: $-(\lambda-1)$；行$\{1,2\}$列$\{2,3\}$: $0$；行$\{1,3\}$列$\{1,2\}$: $0$；行$\{1,3\}$列$\{1,3\}$: $(\lambda-1)(\lambda-2)$；行$\{1,3\}$列$\{2,3\}$: $0$；行$\{2,3\}$列$\{1,2\}$: $0$；行$\{2,3\}$列$\{1,3\}$: $0$；行$\{2,3\}$列$\{2,3\}$: $(\lambda-2)^2$。

    非零子式：$(\lambda-1)(\lambda-2)$，$-(\lambda-1)$，$(\lambda-1)(\lambda-2)$，$(\lambda-2)^2$。$D_2 = 1$（因为 $\gcd$ 包含 $-(\lambda-1)$ 和 $(\lambda-2)^2$，它们互素）。

    $D_3 = (\lambda-1)(\lambda-2)^2$。

    $B$ 的不变因子：$d_1 = 1$，$d_2 = 1$，$d_3 = (\lambda-1)(\lambda-2)^2$。

    $A$ 与 $B$ 的不变因子不同，故 $A \not\sim B$。

    直观理解：$A$ 的特征值 $2$ 有两个线性无关的特征向量（$A$ 在 $\lambda=2$ 处可对角化），而 $B$ 的特征值 $2$ 只有一个线性无关的特征向量（$B$ 在 $\lambda=2$ 处有非平凡的 Jordan 块）。

!!! example "例 13B.15"
    以下两个 $4 \times 4$ 矩阵有相同的特征多项式 $(\lambda-1)^4$ 和相同的最小多项式 $(\lambda-1)^2$，但不相似：

    $$
    A_1 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix}, \quad A_2 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
    $$

    $A_1$ 的不变因子：$1, 1, (\lambda-1)^2, (\lambda-1)^2$。

    $A_2$ 的不变因子：$1, 1, \lambda-1, (\lambda-1)^3$。

    等等——$A_2$ 的最小多项式应该是 $(\lambda-1)^2$ 吗？$(A_2 - I)^2 = \begin{pmatrix} 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix}$？$(A_2 - I) = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix}$，$(A_2-I)^2 = 0$。是的，最小多项式是 $(\lambda-1)^2$。

    正确的不变因子：$A_1$：$d_1 = d_2 = 1$，$d_3 = (\lambda-1)^2$，$d_4 = (\lambda-1)^2$。$A_2$：$d_1 = d_2 = d_3 = 1$，$d_4 = (\lambda-1)^4$。但 $d_4 = (\lambda-1)^4$ 意味着最小多项式为 $(\lambda-1)^4$，矛盾。

    修正：$A_2 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$ 的 Jordan 形是 $J_2(1) \oplus J_1(1) \oplus J_1(1)$。最小多项式是 $(\lambda-1)^2$。不变因子为 $1, \lambda-1, (\lambda-1), (\lambda-1)^2$。验证乘积：$(\lambda-1) \cdot (\lambda-1) \cdot (\lambda-1)^2 = (\lambda-1)^4$ $\checkmark$，整除链 $(\lambda-1) \mid (\lambda-1) \mid (\lambda-1)^2$ $\checkmark$。

    所以 $A_1$ 和 $A_2$ 确实有相同的特征多项式和最小多项式，但不变因子不同：

    - $A_1$（$J_2(1) \oplus J_2(1)$）：不变因子 $1, 1, (\lambda-1)^2, (\lambda-1)^2$
    - $A_2$（$J_2(1) \oplus J_1(1) \oplus J_1(1)$）：不变因子 $1, \lambda-1, \lambda-1, (\lambda-1)^2$

    这说明特征多项式和最小多项式的组合仍然不是完全不变量。

## 练习题

1. **[概念] 什么是 $\lambda$-矩阵？它与数值矩阵的主要区别是什么？**
   ??? success "参考答案"
       $\lambda$-矩阵是以多项式为元素的矩阵。与数值矩阵不同，它的初等变换乘数必须是**非零常数**（多项式环的单位），且其可逆性要求行列式为非零常数而非仅仅非零。

2. **[Smith标准形] 计算 $\lambda$-矩阵 $\begin{pmatrix} \lambda & 1 \\ \lambda^2 & \lambda \end{pmatrix}$ 的 Smith 标准形。**
   ??? success "参考答案"
       通过 $r_2 \to r_2 - \lambda r_1$，矩阵变为 $\begin{pmatrix} \lambda & 1 \\ 0 & 0 \end{pmatrix}$。再交换行列可化为 $\operatorname{diag}(1, 0)$。注意其秩为 1。

3. **[不变因子] 矩阵 $A$ 的不变因子满足整除链 $d_1 \mid d_2 \mid \dots \mid d_n$。最后一个不变因子 $d_n(\lambda)$ 具有什么重要的解析意义？**
   ??? success "参考答案"
       $d_n(\lambda)$ 恰好等于矩阵 $A$ 的**最小多项式** $m_A(\lambda)$。

4. **[初等因子] 设特征矩阵的 Smith 标准形对角元为 $\{1, 1, \lambda-1, (\lambda-1)(\lambda-2)^2\}$。列出其初等因子组。**
   ??? success "参考答案"
       初等因子是每一个不变因子分解出的最高幂次：$(\lambda-1)$ [来自第3个], $(\lambda-1)$ [来自第4个], $(\lambda-2)^2$ [来自第4个]。

5. **[有理标准形] 为什么我们说有理标准形（Rational Canonical Form）比 Jordan 标准形更具有“普适性”？**
   ??? success "参考答案"
       因为 Jordan 形的构造依赖于特征值在基础域上的存在（代数闭性）。如果矩阵在有理数域 $\mathbb{Q}$ 上且特征值是无理数，Jordan 形就无法在 $\mathbb{Q}$ 上定义。而有理标准形只涉及多项式的系数，在任何域上都存在。

6. **[友矩阵] 写出多项式 $p(\lambda) = \lambda^3 - 2\lambda + 5$ 的友矩阵。**
   ??? success "参考答案"
       $\begin{pmatrix} 0 & 0 & -5 \\ 1 & 0 & 2 \\ 0 & 1 & 0 \end{pmatrix}$。注意常数项在最后一行或最后一列（取决于约定），这里采用常见列形式。

7. **[完全不变量] 判定两个矩阵相似的充要条件是什么（利用 $\lambda$-矩阵理论）？**
   ??? success "参考答案"
       它们的特征矩阵 $\lambda I - A$ 和 $\lambda I - B$ 拥有完全相同的 **Smith 标准形**（或者等价地，拥有相同的不变因子组）。

8. **[判定] 两个 $3 \times 3$ 矩阵拥有相同的特征多项式和相同的最小多项式，它们一定相似吗？如果是 $4 \times 4$ 呢？**
   ??? success "参考答案"
       对于 $n \le 3$，特征多项式和最小多项式足以确定不变因子，故一定相似。但对于 $n \ge 4$，存在反例（如两个矩阵都有最小多项式 $(\lambda-1)^2$ 和特征多项式 $(\lambda-1)^4$，但一个包含两个 $2 \times 2$ 块，另一个包含一个 $2 \times 2$ 和两个 $1 \times 1$ 块）。

9. **[Jordan转化] 初等因子 $(\lambda - 5)^3$ 对应 Jordan 标准形中的哪一部分？**
   ??? success "参考答案"
       对应一个 $3 \times 3$ 的 Jordan 块 $J_3(5)$。

10. **[爱因斯坦思考题] 宇宙中的基本粒子由其“量子数”（如自旋、电荷）唯一标识。在 λ-矩阵理论中，什么对象扮演了类似的“矩阵指纹”角色，使得我们无论如何旋转坐标系（基变换）都能识别出它？**
    ??? success "参考答案"
        是**不变因子（Invariant Factors）**。它们是矩阵在相似变换下的“遗传密码”。无论矩阵表现得多么错综复杂，只要我们进入 $\lambda$ 空间（多项式环），通过 Smith 标准形就能剥离掉坐标系的伪装，露出那组唯一确定的多项式链。它们是线性变换最内禀、最本质的代数基因。

## 本章小结

本章通过将数值矩阵提升为多项式矩阵，彻底解决了一般域上的矩阵相似分类问题：

1. **λ-矩阵理论**：引入了以多项式为元素的矩阵及其初等变换，建立了比数值矩阵更广阔的代数框架。
2. **Smith 标准形**：证明了任意 λ-矩阵都等价于唯一的对角形，其对角元（不变因子）构成的整除链是矩阵的深度不变量。
3. **相似的本质判据**：确立了矩阵相似等价于其特征矩阵等价的真理，指出不变因子组是判定相似的完全不变量。
4. **有理标准形**：构造了由友矩阵组成的块对角矩阵，实现了在不依赖特征值分解的前提下，对任意域上矩阵的最简分解。
5. **理论的统一**：展示了初等因子如何作为桥梁，将抽象的有理标准形与直观的 Jordan 标准形统一起来，完美终结了有限维线性算子的结构分析。

