# 第 59 章 热带线性代数

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 特征值(Ch6) · 图论(Ch27)

**本章脉络**：热带半环 $(\mathbb{R}\cup\{-\infty\}, \max, +)$ → 热带矩阵乘法 → 热带行列式 → 热带特征值 → 临界图 → 最短路径 → 调度问题 → 热带凸性 → 热带代数几何初步

**延伸**：热带几何将代数几何中的簇"热带化"为分段线性对象，为计数代数曲线提供了组合工具；离散事件系统（制造、交通网络）中的最优调度问题可自然地表述为热带线性方程组

</div>

热带线性代数是经典线性代数在一个奇异代数结构——**热带半环**——上的对应物。在热带半环中，加法被取最大值运算替代，乘法被普通加法替代。这一看似"古怪"的重新定义绝非无用的抽象游戏。热带代数为最短路径问题、调度优化、离散事件系统提供了统一的代数框架，而热带几何则将代数几何中的多项式簇"退化"为分段线性对象，为枚举几何中的深层计数问题（如 Mikhalkin 的热带曲线计数）提供了惊人的组合方法。

"热带"一词来源于巴西数学家 Imre Simon 的工作，他是该领域的先驱之一。

---

## 59.1 热带半环

<div class="context-flow" markdown>

**核心问题**：什么是热带半环？它与通常的实数域有什么本质区别？为什么热带代数不构成环？

</div>

!!! definition "定义 59.1 (热带半环)"
    **热带半环**（tropical semiring）是集合 $\mathbb{T} = \mathbb{R} \cup \{-\infty\}$，配备两种运算：

    - **热带加法**：$a \oplus b = \max(a, b)$。
    - **热带乘法**：$a \odot b = a + b$（普通加法）。

    热带加法的单位元为 $\mathbb{0} = -\infty$（因为 $\max(a, -\infty) = a$），热带乘法的单位元为 $\mathbb{1} = 0$（因为 $a + 0 = a$）。

!!! theorem "定理 59.1 (热带半环的代数性质)"
    $(\mathbb{T}, \oplus, \odot)$ 满足：

    (a) $(\mathbb{T}, \oplus)$ 是交换幺半群，单位元 $-\infty$。

    (b) $(\mathbb{T}, \odot)$ 是交换幺半群，单位元 $0$。

    (c) 分配律：$a \odot (b \oplus c) = (a \odot b) \oplus (a \odot c)$，即 $a + \max(b,c) = \max(a+b, a+c)$。

    (d) $\mathbb{0} = -\infty$ 是零化子：$a \odot (-\infty) = a + (-\infty) = -\infty$。

    (e) $(\mathbb{T}, \oplus, \odot)$ **不是环**：热带加法没有逆元。$\max(a, b) = -\infty$ 要求 $a = b = -\infty$，因此对 $a \neq -\infty$，不存在 $b$ 使得 $a \oplus b = -\infty$。

??? proof "证明"
    **(a)-(d)** 直接验证。

    **(e)** 假设对某 $a \in \mathbb{R}$ 存在加法逆元 $b$，使得 $\max(a, b) = -\infty$。由 $\max$ 的性质，这要求 $a \leq -\infty$ 且 $b \leq -\infty$，即 $a = b = -\infty$。但我们假设 $a \in \mathbb{R}$，矛盾。

    因此热带加法 $\oplus$ 没有逆元，$\mathbb{T}$ 不构成环，只构成半环。$\blacksquare$

!!! definition "定义 59.2 (min-plus 半环)"
    **min-plus 热带半环**是 $\mathbb{T}_{\min} = \mathbb{R} \cup \{+\infty\}$，加法为 $a \oplus_{\min} b = \min(a,b)$，乘法为 $a \odot b = a + b$。加法单位元为 $+\infty$，乘法单位元为 $0$。

    max-plus 和 min-plus 两种约定在文献中都被广泛使用。本章采用 max-plus 约定。两者通过 $a \mapsto -a$ 相互转换。

!!! example "例 59.1"
    在热带半环中：

    - $3 \oplus 5 = \max(3, 5) = 5$
    - $3 \odot 5 = 3 + 5 = 8$
    - $2 \odot (3 \oplus 7) = 2 + \max(3, 7) = 2 + 7 = 9$
    - $(2 \odot 3) \oplus (2 \odot 7) = \max(5, 9) = 9$（验证分配律）
    - $(-\infty) \oplus 3 = 3$，$(-\infty) \odot 3 = -\infty$（单位元和零化子）

!!! definition "定义 59.3 (热带幂)"
    在热带半环中，$a$ 的 $n$ 次热带幂为

    $$a^{\odot n} = \underbrace{a \odot a \odot \cdots \odot a}_{n} = n \cdot a \quad (\text{普通乘法}).$$

    热带"多项式" $p(x) = \bigoplus_{i=0}^n a_i \odot x^{\odot i} = \max_{0 \leq i \leq n}(a_i + i \cdot x)$ 是 $x$ 的分段线性凸函数。

---

## 59.2 热带矩阵运算

<div class="context-flow" markdown>

**核心问题**：如何定义热带矩阵的加法和乘法？热带矩阵乘法与经典矩阵乘法有何区别？

</div>

!!! definition "定义 59.4 (热带矩阵运算)"
    设 $A = (a_{ij}) \in \mathbb{T}^{m \times p}$，$B = (b_{jk}) \in \mathbb{T}^{p \times n}$。

    **热带矩阵加法**：$(A \oplus B)_{ij} = a_{ij} \oplus b_{ij} = \max(a_{ij}, b_{ij})$（要求 $A, B$ 同型）。

    **热带矩阵乘法**：

    $$(A \odot B)_{ik} = \bigoplus_{j=1}^p a_{ij} \odot b_{jk} = \max_{1 \leq j \leq p}(a_{ij} + b_{jk}).$$

!!! definition "定义 59.5 (热带单位矩阵)"
    $n \times n$ 的**热带单位矩阵**为

    $$I_n^{\oplus} = \begin{pmatrix} 0 & -\infty & \cdots & -\infty \\ -\infty & 0 & \cdots & -\infty \\ \vdots & \ddots & \ddots & \vdots \\ -\infty & \cdots & -\infty & 0\end{pmatrix},$$

    即对角元素为 $\mathbb{1} = 0$，非对角元素为 $\mathbb{0} = -\infty$。满足 $A \odot I^{\oplus} = I^{\oplus} \odot A = A$。

!!! theorem "定理 59.2 (热带矩阵乘法的性质)"
    热带矩阵乘法满足：

    (a) **结合律**：$(A \odot B) \odot C = A \odot (B \odot C)$。

    (b) **分配律**：$A \odot (B \oplus C) = (A \odot B) \oplus (A \odot C)$。

    (c) **不满足消去律**：$A \odot B = A \odot C$ 不能推出 $B = C$。

??? proof "证明"
    **(a)** $((A \odot B) \odot C)_{il} = \max_k(\max_j(a_{ij}+b_{jk}) + c_{kl}) = \max_{j,k}(a_{ij}+b_{jk}+c_{kl})$。
    同理 $(A \odot (B \odot C))_{il} = \max_j(a_{ij}+\max_k(b_{jk}+c_{kl})) = \max_{j,k}(a_{ij}+b_{jk}+c_{kl})$。两者相等。

    **(b)** $(A \odot (B \oplus C))_{ik} = \max_j(a_{ij}+\max(b_{jk},c_{jk})) = \max_j\max(a_{ij}+b_{jk}, a_{ij}+c_{jk})$
    $= \max(\max_j(a_{ij}+b_{jk}), \max_j(a_{ij}+c_{jk})) = ((A\odot B)\oplus(A\odot C))_{ik}$。

    **(c)** 反例：$A = (0)$，$B = (1)$，$C = (0)$。$A \odot B = (1) = A \odot C$ 不成立，但考虑 $A = (5, 3)$，$B = \binom{1}{2}$，$C = \binom{0}{3}$。$A \odot B = \max(5+1, 3+2) = 6 = \max(5+0, 3+3) = A \odot C$，但 $B \neq C$。$\blacksquare$

!!! example "例 59.2"
    设 $A = \begin{pmatrix}3 & 2\\1 & 4\end{pmatrix}$，$B = \begin{pmatrix}1 & 0\\2 & 3\end{pmatrix}$（元素在 $\mathbb{T}$ 中）。

    $$(A \odot B)_{11} = \max(3+1, 2+2) = \max(4, 4) = 4,$$
    $$(A \odot B)_{12} = \max(3+0, 2+3) = \max(3, 5) = 5,$$
    $$(A \odot B)_{21} = \max(1+1, 4+2) = \max(2, 6) = 6,$$
    $$(A \odot B)_{22} = \max(1+0, 4+3) = \max(1, 7) = 7.$$

    因此 $A \odot B = \begin{pmatrix}4 & 5\\6 & 7\end{pmatrix}$。

!!! definition "定义 59.6 (热带矩阵幂)"
    $n \times n$ 热带矩阵 $A$ 的热带 $k$ 次幂定义为

    $$A^{\odot k} = \underbrace{A \odot A \odot \cdots \odot A}_{k},$$

    $(A^{\odot k})_{ij} = \max_{i_1, \ldots, i_{k-1}}(a_{i,i_1} + a_{i_1,i_2} + \cdots + a_{i_{k-1},j})$。

    这就是有向加权图中从 $i$ 到 $j$ 的所有**恰好 $k$ 步**路径中权重之和的最大值。

---

## 59.3 热带行列式与永久式

<div class="context-flow" markdown>

**核心问题**：热带行列式如何定义？它与经典行列式和指派问题有什么关系？

</div>

!!! definition "定义 59.7 (热带行列式)"
    $n \times n$ 热带矩阵 $A$ 的**热带行列式**定义为

    $$\mathrm{tdet}(A) = \bigoplus_{\sigma \in S_n} \bigodot_{i=1}^n a_{i,\sigma(i)} = \max_{\sigma \in S_n} \sum_{i=1}^n a_{i,\sigma(i)},$$

    其中 $S_n$ 是 $n$ 元置换群。

!!! theorem "定理 59.3 (热带行列式 = 热带永久式)"
    在经典代数中，行列式和永久式不同（行列式含符号因子 $\mathrm{sgn}(\sigma)$）。但在热带半环中，

    $$\mathrm{tdet}(A) = \mathrm{tperm}(A) = \max_{\sigma \in S_n}\sum_{i=1}^n a_{i,\sigma(i)}.$$

    原因是热带加法 $\oplus = \max$ **是幂等的**（$a \oplus a = a$），因此符号因子无法在热带半环中定义（$-\infty$ 是加法单位元，不存在"$-a$"），行列式与永久式的区别消失。

??? proof "证明"
    经典行列式 $\det(A) = \sum_{\sigma} \mathrm{sgn}(\sigma)\prod_i a_{i,\sigma(i)}$。热带化后 $\sum \to \max$，$\prod \to \sum$，$\mathrm{sgn}(\sigma) \in \{+1, -1\}$。

    但 $\mathrm{sgn}(\sigma)$ 在热带半环中无意义——热带加法没有逆元，不存在"$-1$"。因此热带行列式自然地定义为不带符号的版本，即永久式的热带化。$\blacksquare$

!!! theorem "定理 59.4 (热带行列式与指派问题)"
    热带行列式的计算等价于**最优指派问题**（assignment problem）：

    $$\mathrm{tdet}(A) = \max_{\sigma \in S_n}\sum_{i=1}^n a_{i,\sigma(i)}$$

    是将 $n$ 个工人分配到 $n$ 个任务的最大权重完美匹配。

    最优指派问题可以用 **Hungarian 算法**（匈牙利算法）在 $O(n^3)$ 时间内求解，因此热带行列式可在多项式时间内计算——这与经典永久式（#P 困难）形成鲜明对比。

!!! example "例 59.3"
    设 $A = \begin{pmatrix}5 & 3 & 2\\4 & 6 & 1\\3 & 2 & 7\end{pmatrix}$。

    $\mathrm{tdet}(A) = \max_\sigma \sum_i a_{i,\sigma(i)}$：

    - $\sigma = (1,2,3)$（恒等）：$5+6+7 = 18$
    - $\sigma = (1,3,2)$：$5+1+2 = 8$
    - $\sigma = (2,1,3)$：$3+4+7 = 14$
    - $\sigma = (2,3,1)$：$3+1+3 = 7$
    - $\sigma = (3,1,2)$：$2+4+2 = 8$
    - $\sigma = (3,2,1)$：$2+6+3 = 11$

    因此 $\mathrm{tdet}(A) = 18$，由恒等置换达到。

!!! definition "定义 59.8 (热带奇异性)"
    热带矩阵 $A$ 称为**热带奇异的**，如果热带行列式的最大值由**两个或更多**置换同时达到。否则 $A$ 称为**热带非奇异的**。

    在例 59.3 中，最大值仅由恒等置换达到，因此 $A$ 是热带非奇异的。

---

## 59.4 热带特征值

<div class="context-flow" markdown>

**核心问题**：如何定义和计算热带矩阵的特征值？热带特征值与有向图的圈结构有何关系？

</div>

!!! definition "定义 59.9 (热带特征值问题)"
    设 $A \in \mathbb{T}^{n \times n}$。标量 $\lambda \in \mathbb{T}$ 和非零向量 $x \in \mathbb{T}^n$（$x \neq (-\infty, \ldots, -\infty)$）称为 $A$ 的**热带特征值**和**热带特征向量**，如果

    $$A \odot x = \lambda \odot x,$$

    即对所有 $i = 1, \ldots, n$：

    $$\max_{j=1,\ldots,n}(a_{ij} + x_j) = \lambda + x_i.$$

!!! definition "定义 59.10 (有向图的最大环均值)"
    设 $G(A)$ 是 $A$ 的**关联有向图**：顶点集 $\{1, \ldots, n\}$，边 $(i,j)$ 存在当且仅当 $a_{ij} \neq -\infty$，权重为 $a_{ij}$。

    有向图中一个环 $\gamma = (i_1, i_2, \ldots, i_k, i_1)$ 的**环均值**为

    $$\mu(\gamma) = \frac{a_{i_1 i_2} + a_{i_2 i_3} + \cdots + a_{i_k i_1}}{k}.$$

    **最大环均值**为 $\lambda^*(A) = \max_\gamma \mu(\gamma)$，最大取遍 $G(A)$ 的所有有向环。

!!! theorem "定理 59.5 (热带 Perron-Frobenius 定理)"
    设 $A \in \mathbb{T}^{n \times n}$，其关联有向图 $G(A)$ 是**强连通的**（从任意顶点到任意其他顶点都存在有向路径）。则：

    (a) $A$ 有唯一的热带特征值 $\lambda = \lambda^*(A)$（最大环均值）。

    (b) 热带特征向量 $x$ 存在（但不一定唯一）。

??? proof "证明"
    **(a) 存在性。** 设 $\lambda = \lambda^*(A)$。构造 $A_\lambda = A - \lambda \odot J$（热带意义下减去 $\lambda$），即 $(A_\lambda)_{ij} = a_{ij} - \lambda$。则 $A_\lambda$ 的所有环均值 $\leq 0$，且至少有一个环均值 $= 0$（达到最大环均值的临界环）。

    定义 $x_i = \max_{j \in C}$ 从临界环上某点 $j$ 到 $i$ 的最大权重路径（在 $A_\lambda$ 中）。可以验证 $x$ 满足热带特征方程。

    **(a) 唯一性。** 假设 $\lambda'$ 也是特征值，$A \odot x' = \lambda' \odot x'$。设 $\lambda' > \lambda$。考虑矩阵 $B = A - \lambda'$（逐元素减），则所有环均值 $< 0$。但 $B \odot x' = 0 \odot x' = x'$，即 $\max_j(b_{ij} + x'_j) = x'_i$，这要求存在从 $i$ 到某 $j$ 的边使得 $b_{ij} + x'_j = x'_i$。沿这样的边走下去最终形成环，但环的总权重为负，矛盾。

    因此 $\lambda^*(A)$ 是唯一的热带特征值。$\blacksquare$

!!! definition "定义 59.11 (临界图)"
    矩阵 $A$ 的**临界图** $G_c(A)$ 是由达到最大环均值的所有环及其上的顶点和边构成的子图。

    临界图在调度理论中有重要意义：临界图上的环对应系统的"瓶颈"环路。

!!! example "例 59.4"
    设 $A = \begin{pmatrix}2 & 5\\3 & 1\end{pmatrix}$。

    $G(A)$ 有四条边：$(1,1)$ 权 2，$(1,2)$ 权 5，$(2,1)$ 权 3，$(2,2)$ 权 1。

    环及其均值：

    - 自环 $(1,1)$：均值 $2/1 = 2$。
    - 自环 $(2,2)$：均值 $1/1 = 1$。
    - 环 $(1,2,1)$：均值 $(5+3)/2 = 4$。

    最大环均值 $\lambda = 4$。验证：取 $x = (0, -1)$（热带意义下），

    $$\max(2+0, 5+(-1)) = \max(2, 4) = 4 = 4 + 0, \quad \checkmark$$
    $$\max(3+0, 1+(-1)) = \max(3, 0) = 3 \neq 4 + (-1) = 3. \quad \checkmark$$

    因此 $\lambda = 4$，$x = (0, -1)$ 是特征对。

---

## 59.5 最短路径与 Kleene 星

<div class="context-flow" markdown>

**核心问题**：热带矩阵幂与最短/最长路径问题有什么关系？Kleene 星算子如何给出全源最短路径？

</div>

!!! definition "定义 59.12 (Kleene 星)"
    设 $A \in \mathbb{T}^{n \times n}$。$A$ 的 **Kleene 星**（或闭包）定义为

    $$A^* = I^{\oplus} \oplus A \oplus A^{\odot 2} \oplus A^{\odot 3} \oplus \cdots = \bigoplus_{k=0}^{\infty} A^{\odot k},$$

    即 $(A^*)_{ij} = \sup_{k \geq 0} (A^{\odot k})_{ij}$。这在热带半环中是从 $i$ 到 $j$ 的所有路径（包括零步路径 $i = j$）中权重之和的最大值。

!!! theorem "定理 59.6 (Kleene 星的收敛性)"
    $A^*$ 在 $\mathbb{T}$ 中有限当且仅当 $A$ 的所有有向环的权重之和 $\leq 0$（在 max-plus 意义下，即 $\lambda^*(A) \leq 0$）。此时

    $$A^* = I^{\oplus} \oplus A \oplus A^{\odot 2} \oplus \cdots \oplus A^{\odot (n-1)}.$$

??? proof "证明"
    如果存在权重之和为正的环 $\gamma$，则沿 $\gamma$ 无限绕行可使路径权重趋于 $+\infty$，$A^*$ 发散。

    反之，如果所有环权重 $\leq 0$，则最优路径不需要经过重复顶点（重复会产生非正环，不改善或恶化权重）。$n$ 个顶点的图中无重复顶点的路径最多 $n-1$ 步，因此 $A^* = \bigoplus_{k=0}^{n-1} A^{\odot k}$。$\blacksquare$

!!! theorem "定理 59.7 (Floyd-Warshall 与热带矩阵幂)"
    设 $D \in \mathbb{T}_{\min}^{n \times n}$ 是有向加权图的邻接矩阵（使用 min-plus 半环）。则：

    (a) $D^{\odot k}_{\min}$ 的 $(i,j)$ 元素是从 $i$ 到 $j$ 的恰好 $k$ 步最短路径长度。

    (b) $D^*_{\min} = \bigoplus_{k=0}^{n-1} D^{\odot k}_{\min}$ 的 $(i,j)$ 元素是从 $i$ 到 $j$ 的最短路径长度（假设无负环）。

    (c) **Floyd-Warshall 算法**可以视为逐步消元计算 $D^*_{\min}$ 的过程。

!!! example "例 59.5"
    **最长路径的热带矩阵计算。** 考虑有向图，邻接矩阵（max-plus）

    $$A = \begin{pmatrix}-\infty & 3 & -\infty\\-\infty & -\infty & 2\\-\infty & -\infty & -\infty\end{pmatrix}.$$

    表示边 $(1,2)$ 权 3，$(2,3)$ 权 2。

    $A^{\odot 2}$：$(A^{\odot 2})_{13} = \max(a_{11}+a_{13}, a_{12}+a_{23}, a_{13}+a_{33}) = \max(-\infty, 3+2, -\infty) = 5$。

    $A^*$：$(A^*)_{13} = \max(0, -\infty, 5) = 5$（从 1 到 3 的最长路径权重为 $3 + 2 = 5$）。

!!! theorem "定理 59.8 (Bellman-Ford 方程的热带表述)"
    设 $d = (d_1, \ldots, d_n)^\top$ 是从源点 $s$ 到各点的最短距离向量（min-plus 约定）。Bellman-Ford 方程

    $$d_i = \min\!\bigl(0 \cdot [i=s],\, \min_{j: (j,i) \in E}(d_j + w_{ji})\bigr)$$

    可以写为热带线性方程

    $$d = D_{\min} \odot d \oplus e_s,$$

    其中 $e_s$ 是第 $s$ 个分量为 $0$、其余为 $+\infty$ 的向量。迭代 $d^{(k+1)} = D_{\min} \odot d^{(k)} \oplus e_s$ 就是 Bellman-Ford 松弛过程。

---

## 59.6 调度与离散事件系统

<div class="context-flow" markdown>

**核心问题**：热带线性代数如何为调度问题和离散事件系统提供统一的数学框架？

</div>

!!! example "例 59.6"
    **火车时刻表的热带模型。** 考虑一个简单的铁路网络，有两个站 $A$、$B$，火车在两站之间往返。设：

    - $x_1(k)$ = 第 $k$ 次从 $A$ 出发的时间。
    - $x_2(k)$ = 第 $k$ 次从 $B$ 出发的时间。
    - 从 $A$ 到 $B$ 需时 3（单位），从 $B$ 到 $A$ 需时 2。
    - 在 $A$ 站最少停留 1，在 $B$ 站最少停留 1。

    约束：$x_1(k+1) \geq x_2(k) + 2$（从 $B$ 返回后），$x_1(k+1) \geq x_1(k) + 1$（在 $A$ 站停留）。
    类似地，$x_2(k+1) \geq x_1(k) + 3$，$x_2(k+1) \geq x_2(k) + 1$。

    用热带记号：$x(k+1) = A \odot x(k)$，其中 $A = \begin{pmatrix}1 & 2\\3 & 1\end{pmatrix}$。

    热带特征值 $\lambda = \max(1/1, 1/1, (2+3)/2) = 5/2 = 2.5$。这意味着系统的渐近周期为 $2.5$：稳定运行后每 $2.5$ 个时间单位发一班车。

!!! theorem "定理 59.9 (热带线性递推的渐近行为)"
    设 $A \in \mathbb{T}^{n \times n}$，$G(A)$ 强连通。考虑热带线性递推

    $$x(k+1) = A \odot x(k), \quad x(0) = x_0.$$

    则对任意初始向量 $x_0$（各分量有限），

    $$\lim_{k \to \infty} \frac{x_i(k)}{k} = \lambda \quad \forall\, i = 1, \ldots, n,$$

    其中 $\lambda = \lambda^*(A)$ 是最大环均值。更精确地，存在 $K$ 和周期 $p$ 使得对所有 $k \geq K$：

    $$x(k+p) = \lambda^{\odot p} \odot x(k) \quad (\text{即 } x_i(k+p) = p\lambda + x_i(k)).$$

??? proof "证明"
    核心思想：$x(k) = A^{\odot k} \odot x_0$。$(A^{\odot k})_{ij}$ 是从 $j$ 到 $i$ 的所有 $k$ 步路径中权重之和的最大值。

    当 $k$ 充分大时，最优 $k$ 步路径包含许多绕临界环的"循环"。每绕一次临界环（设长度为 $\ell$），路径权重增加 $\ell \lambda$。因此 $(A^{\odot k})_{ij} \approx k\lambda +$ 常数项。

    精确证明需要分析 $A$ 的临界图和非临界部分的结构，可得到周期 $p$ 是临界环长度的最小公倍数。$\blacksquare$

!!! example "例 59.7"
    **制造系统调度。** 一条装配线有 3 个工序，工序 $i$ 的处理时间为 $p_i$，工序间有先后约束。设

    $$A = \begin{pmatrix}p_1 & d_{12} & -\infty\\-\infty & p_2 & d_{23}\\d_{31} & -\infty & p_3\end{pmatrix},$$

    其中 $d_{ij}$ 是工序 $i$ 完成后到工序 $j$ 开始的延迟。最大环均值给出了系统的最小循环时间（瓶颈）。

---

## 59.7 热带凸性

<div class="context-flow" markdown>

**核心问题**：热带半环中是否存在类似经典凸集和凸包的概念？

</div>

!!! definition "定义 59.13 (热带凸组合)"
    $\mathbb{T}^n$ 中点 $x_1, \ldots, x_k$ 的**热带凸组合**是

    $$\bigoplus_{i=1}^k \lambda_i \odot x_i = \Bigl(\max_{i}(\lambda_i + x_{i,1}),\, \ldots,\, \max_{i}(\lambda_i + x_{i,n})\Bigr),$$

    其中 $\lambda_1, \ldots, \lambda_k \in \mathbb{T}$ 且 $\bigoplus_i \lambda_i = \max_i \lambda_i = 0$（热带归一化条件）。

!!! definition "定义 59.14 (热带凸包)"
    点集 $S \subset \mathbb{T}^n$ 的**热带凸包**是 $S$ 中所有有限子集的热带凸组合的集合：

    $$\mathrm{tconv}(S) = \Bigl\{\bigoplus_{i=1}^k \lambda_i \odot x_i : x_i \in S,\, \lambda_i \in \mathbb{T},\, \max_i \lambda_i = 0,\, k \in \mathbb{N}\Bigr\}.$$

!!! theorem "定理 59.10 (热带凸集的几何)"
    热带凸集（在 $\mathbb{R}^n/\mathbb{R}\mathbf{1}$ 的商空间中，即模去同时平移）是**分段线性**的。热带凸多面体可以用有限多个热带半空间的交来描述：

    $$\{x \in \mathbb{T}^n : a_i \odot x_{j_i} \oplus b_i \odot x_{k_i} \leq c_i \odot x_{\ell_i}\}.$$

!!! example "例 59.8"
    在 $\mathbb{T}^2/\mathbb{R}\mathbf{1}$ 中（即 $\mathbb{R}^2$ 模去平移 $(t, t)$），三点 $p = (0,0)$、$q = (2,0)$、$r = (0,3)$ 的热带凸包是一个"热带三角形"。

    热带线段 $\overline{pq}^{\mathrm{trop}}$ = $\{(\max(\lambda, \mu+2), \max(\lambda, \mu)) : \max(\lambda, \mu) = 0\}$，在商空间中对应两条射线与一段水平线段的并。

---

## 59.8 热带代数几何初步

<div class="context-flow" markdown>

**核心问题**：经典代数几何中的概念如何"热带化"？热带簇有怎样的组合结构？

</div>

!!! definition "定义 59.15 (热带多项式)"
    $\mathbb{T}$ 上的**热带多项式**是

    $$p(x_1, \ldots, x_n) = \bigoplus_{i \in S} a_i \odot x_1^{\odot i_1} \odot \cdots \odot x_n^{\odot i_n} = \max_{i \in S}\Bigl(a_i + \sum_{k=1}^n i_k x_k\Bigr),$$

    其中 $S$ 是有限支撑集，$i = (i_1, \ldots, i_n)$ 是多重指标。热带多项式是 $(x_1, \ldots, x_n)$ 的**分段线性凸函数**。

!!! definition "定义 59.16 (热带超曲面)"
    热带多项式 $p$ 的**热带超曲面**（或热带簇）$\mathcal{T}(p)$ 是 $p$ **不光滑的点**的集合，即最大值由两个或更多项同时达到的点：

    $$\mathcal{T}(p) = \{x \in \mathbb{R}^n : \text{最大值在 } p(x) \text{中由至少两项达到}\}.$$

!!! theorem "定理 59.11 (Kapranov 定理)"
    设 $f(x) = \sum_{i \in S} c_i x^i$ 是一个一元多项式（在 Puiseux 级数域 $\mathbb{K}$ 上），$\mathrm{val}: \mathbb{K}^* \to \mathbb{R}$ 是赋值映射。则 $f$ 的零点 $\alpha_1, \ldots, \alpha_d$ 的**赋值的多重集** $\{\mathrm{val}(\alpha_1), \ldots, \mathrm{val}(\alpha_d)\}$ 恰好是热带多项式 $\mathrm{trop}(f)(x) = \max_i(\mathrm{val}(c_i) + ix)$ 的热带根（超曲面中的点）的多重集。

!!! example "例 59.9"
    考虑热带二次多项式

    $$p(x) = 3 \odot x^{\odot 2} \oplus 5 \odot x \oplus 2 = \max(3 + 2x,\, 5 + x,\, 2).$$

    三条线 $y_1 = 3+2x$，$y_2 = 5+x$，$y_3 = 2$ 的交点：

    - $y_1 = y_2$：$3+2x = 5+x \Rightarrow x = 2$。
    - $y_2 = y_3$：$5+x = 2 \Rightarrow x = -3$。

    热带根为 $x = 2$ 和 $x = -3$（对应两个"弯折点"）。

!!! theorem "定理 59.12 (热带 Bezout 定理)"
    设 $p, q$ 是 $\mathbb{T}$ 上两个一元热带多项式，次数分别为 $m, n$。则 $\mathcal{T}(p) \cap \mathcal{T}(q)$ 的交点数（计重数）至多为 $mn$。

    更一般地，$\mathbb{R}^2$ 中两条热带曲线（次数为 $m, n$）的交点数（在适当的重数定义下）恰好为 $mn$（当曲线"一般位置"时）。

!!! definition "定义 59.17 (阿米巴)"
    设 $V \subset (\mathbb{C}^*)^n$ 是代数簇。$V$ 的**阿米巴**（amoeba）是

    $$\mathcal{A}(V) = \{(\log|z_1|, \ldots, \log|z_n|) : (z_1, \ldots, z_n) \in V\} \subset \mathbb{R}^n.$$

    Mikhalkin 和 Rullgard 证明：当参数趋于某个极限时，代数曲线的阿米巴"收缩"到对应的热带曲线。这为热带几何与经典代数几何之间建立了深层联系。

!!! example "例 59.10"
    直线 $z_1 + z_2 + 1 = 0$ 在 $(\mathbb{C}^*)^2$ 中的阿米巴是 $\mathbb{R}^2$ 中一个具有三个"触手"的有界区域（向三个方向无限延伸）。其热带极限是三条射线从原点出发：

    $$\mathcal{T}(\max(x, y, 0)) = \{(x, y) : \text{最大值非唯一}\},$$

    由三段构成：$\{x = y \geq 0\}$，$\{x = 0 \geq y\}$，$\{y = 0 \geq x\}$，形成"Y"形图。

---

**本章要点总结：**

1. 热带半环 $(\mathbb{R} \cup \{-\infty\}, \max, +)$ 不是环（加法无逆元），但具有丰富的代数结构。
2. 热带矩阵乘法 $(A \odot B)_{ik} = \max_j(a_{ij}+b_{jk})$ 自然地编码路径优化问题。
3. 热带行列式等于指派问题的最优值，可在多项式时间计算。
4. 热带特征值是有向图的最大环均值，唯一确定（强连通图）。
5. Kleene 星给出全源最长/最短路径，Floyd-Warshall 是其矩阵计算方法。
6. 离散事件系统（调度、交通）的动态方程自然表述为热带线性递推。
7. 热带代数几何将代数簇"退化"为分段线性对象，是现代数学的活跃前沿。
