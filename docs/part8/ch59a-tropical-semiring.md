# 第 59A 章 热带半环与热带矩阵论

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 特征值(Ch6) · 图论(Ch27) · Perron-Frobenius 定理(Ch17)

**本章脉络**：热带半环 $(\mathbb{R}\cup\{-\infty\}, \max, +)$ → min-plus 变体 → 热带矩阵乘法与加权有向图 → 热带行列式 = 最优指派 → 热带奇异性 → 热带特征值 → 热带 Perron-Frobenius 定理 → 临界图 → 周期性定理 → Karp 算法 → Kleene 星与收敛 → Floyd-Warshall · Bellman-Ford → 热带线性方程组 → 热带 Cramer 法则 → 残差理论 → 热带秩 → 离散事件系统

**延伸**：热带半环上的矩阵论为最短/最长路径问题、调度优化、离散事件系统（制造、交通网络）提供了统一的代数框架；热带行列式与最优指派问题的等价性将组合优化与代数联系起来；热带几何初步见第 59B 章

</div>

热带线性代数是经典线性代数在一个奇特的代数结构——**热带半环**——上的对应物。在热带半环中，加法被取最大值运算替代，乘法被普通加法替代。这一看似"古怪"的重新定义绝非抽象游戏。热带代数为最短路径问题、调度优化、离散事件系统提供了统一的代数框架，而热带矩阵论中的行列式、特征值、线性方程组等概念都获得了深刻的组合解释。

"热带"一词来源于巴西数学家 Imre Simon 的工作，他是该领域的先驱之一。法国数学家以此命名致敬他所在的热带地区。

---

## 59A.1 热带半环

<div class="context-flow" markdown>

**核心问题**：什么是热带半环？它与通常的实数域有什么本质区别？为什么热带代数不构成环？

</div>

!!! definition "定义 59A.1 (热带半环)"
    **热带半环**（tropical semiring）是集合 $\mathbb{T} = \mathbb{R} \cup \{-\infty\}$，配备两种运算：

    - **热带加法**：$a \oplus b = \max(a, b)$。
    - **热带乘法**：$a \odot b = a + b$（普通加法）。

    热带加法的单位元为 $\mathbb{0} = -\infty$（因为 $\max(a, -\infty) = a$），热带乘法的单位元为 $\mathbb{1} = 0$（因为 $a + 0 = a$）。

!!! theorem "定理 59A.1 (热带半环的代数性质)"
    $(\mathbb{T}, \oplus, \odot)$ 满足：

    (a) $(\mathbb{T}, \oplus)$ 是交换幺半群，单位元 $-\infty$。

    (b) $(\mathbb{T}, \odot)$ 是交换幺半群，单位元 $0$。

    (c) 分配律：$a \odot (b \oplus c) = (a \odot b) \oplus (a \odot c)$，即 $a + \max(b,c) = \max(a+b, a+c)$。

    (d) $\mathbb{0} = -\infty$ 是零化子：$a \odot (-\infty) = a + (-\infty) = -\infty$。

    (e) **幂等性**：$a \oplus a = \max(a,a) = a$。热带加法是幂等运算。

    (f) $(\mathbb{T}, \oplus, \odot)$ **不是环**：热带加法没有逆元。$\max(a, b) = -\infty$ 要求 $a = b = -\infty$，因此对 $a \neq -\infty$，不存在 $b$ 使得 $a \oplus b = -\infty$。

??? proof "证明"
    **(a)-(d)** 直接验证各公理。

    - 交换律：$\max(a,b) = \max(b,a)$，$a+b = b+a$。
    - 结合律：$\max(a,\max(b,c)) = \max(\max(a,b),c) = \max(a,b,c)$，$(a+b)+c = a+(b+c)$。
    - 单位元：$\max(a,-\infty) = a$，$a+0 = a$。
    - 分配律：$a + \max(b,c) = \max(a+b, a+c)$（因加法是单调的）。
    - 零化子：$a + (-\infty) = -\infty$。

    **(e)** $\max(a,a) = a$ 对所有 $a$ 显然成立。

    **(f)** 假设对某 $a \in \mathbb{R}$ 存在加法逆元 $b$，使得 $\max(a, b) = -\infty$。由 $\max$ 的性质，这要求 $a \leq -\infty$ 且 $b \leq -\infty$，即 $a = b = -\infty$。但我们假设 $a \in \mathbb{R}$，矛盾。

    因此热带加法 $\oplus$ 没有逆元，$\mathbb{T}$ 不构成环，只构成半环。$\blacksquare$

热带半环的幂等性 $a \oplus a = a$ 是其与经典代数最根本的区别。在经典代数中 $a + a = 2a \neq a$（当 $a \neq 0$ 时），而在热带代数中"加两次和加一次一样"。这一特性导致许多经典代数中的概念（如加法逆元、行列式的符号）在热带世界中失效或需要重新诠释。

---

## 59A.2 Min-plus 半环

!!! definition "定义 59A.2 (min-plus 半环)"
    **min-plus 热带半环**是 $\mathbb{T}_{\min} = \mathbb{R} \cup \{+\infty\}$，加法为 $a \oplus_{\min} b = \min(a,b)$，乘法为 $a \odot b = a + b$。加法单位元为 $+\infty$，乘法单位元为 $0$。

    max-plus 和 min-plus 两种约定在文献中都被广泛使用。两者通过取负映射 $a \mapsto -a$ 相互转换：

    $$\min(a,b) = -\max(-a,-b).$$

    在图论中，min-plus 约定更自然地对应**最短路径**问题，而 max-plus 约定更自然地对应**最长路径**和**调度**问题。本章主要采用 max-plus 约定，但在讨论最短路径算法时会切换到 min-plus。

!!! example "例 59A.1"
    在热带半环（max-plus）中：

    - $3 \oplus 5 = \max(3, 5) = 5$
    - $3 \odot 5 = 3 + 5 = 8$
    - $2 \odot (3 \oplus 7) = 2 + \max(3, 7) = 2 + 7 = 9$
    - $(2 \odot 3) \oplus (2 \odot 7) = \max(5, 9) = 9$（验证分配律）
    - $(-\infty) \oplus 3 = 3$，$(-\infty) \odot 3 = -\infty$（单位元和零化子）

    在 min-plus 半环中：

    - $3 \oplus_{\min} 5 = \min(3, 5) = 3$
    - $3 \odot 5 = 3 + 5 = 8$（乘法相同）
    - $(+\infty) \oplus_{\min} 3 = 3$，$(+\infty) \odot 3 = +\infty$

!!! definition "定义 59A.3 (热带幂)"
    在热带半环中，$a$ 的 $n$ 次热带幂为

    $$a^{\odot n} = \underbrace{a \odot a \odot \cdots \odot a}_{n} = n \cdot a \quad (\text{普通乘法}).$$

    热带"多项式" $p(x) = \bigoplus_{i=0}^n a_i \odot x^{\odot i} = \max_{0 \leq i \leq n}(a_i + i \cdot x)$ 是 $x$ 的分段线性凸函数。

---

## 59A.3 热带矩阵运算

<div class="context-flow" markdown>

**核心问题**：如何定义热带矩阵的加法和乘法？热带矩阵乘法与加权有向图有何关系？

</div>

!!! definition "定义 59A.4 (热带矩阵运算)"
    设 $A = (a_{ij}) \in \mathbb{T}^{m \times p}$，$B = (b_{jk}) \in \mathbb{T}^{p \times n}$。

    **热带矩阵加法**：$(A \oplus B)_{ij} = a_{ij} \oplus b_{ij} = \max(a_{ij}, b_{ij})$（要求 $A, B$ 同型）。

    **热带矩阵乘法**：

    $$(A \odot B)_{ik} = \bigoplus_{j=1}^p a_{ij} \odot b_{jk} = \max_{1 \leq j \leq p}(a_{ij} + b_{jk}).$$

!!! definition "定义 59A.5 (热带单位矩阵)"
    $n \times n$ 的**热带单位矩阵**为

    $$I_n^{\oplus} = \begin{pmatrix} 0 & -\infty & \cdots & -\infty \\ -\infty & 0 & \cdots & -\infty \\ \vdots & \ddots & \ddots & \vdots \\ -\infty & \cdots & -\infty & 0\end{pmatrix},$$

    即对角元素为 $\mathbb{1} = 0$，非对角元素为 $\mathbb{0} = -\infty$。满足 $A \odot I^{\oplus} = I^{\oplus} \odot A = A$。

!!! theorem "定理 59A.2 (热带矩阵乘法的性质)"
    热带矩阵乘法满足：

    (a) **结合律**：$(A \odot B) \odot C = A \odot (B \odot C)$。

    (b) **分配律**：$A \odot (B \oplus C) = (A \odot B) \oplus (A \odot C)$。

    (c) **不满足消去律**：$A \odot B = A \odot C$ 不能推出 $B = C$。

??? proof "证明"
    **(a)** $((A \odot B) \odot C)_{il} = \max_k(\max_j(a_{ij}+b_{jk}) + c_{kl}) = \max_{j,k}(a_{ij}+b_{jk}+c_{kl})$。
    同理 $(A \odot (B \odot C))_{il} = \max_j(a_{ij}+\max_k(b_{jk}+c_{kl})) = \max_{j,k}(a_{ij}+b_{jk}+c_{kl})$。两者相等。

    **(b)** $(A \odot (B \oplus C))_{ik} = \max_j(a_{ij}+\max(b_{jk},c_{jk})) = \max_j\max(a_{ij}+b_{jk}, a_{ij}+c_{jk})$
    $= \max(\max_j(a_{ij}+b_{jk}), \max_j(a_{ij}+c_{jk})) = ((A\odot B)\oplus(A\odot C))_{ik}$。

    **(c)** 反例：取 $A = (5, 3)$，$B = \binom{1}{2}$，$C = \binom{0}{3}$。$A \odot B = \max(5+1, 3+2) = 6 = \max(5+0, 3+3) = A \odot C$，但 $B \neq C$。$\blacksquare$

!!! example "例 59A.2"
    设 $A = \begin{pmatrix}3 & 2\\1 & 4\end{pmatrix}$，$B = \begin{pmatrix}1 & 0\\2 & 3\end{pmatrix}$（元素在 $\mathbb{T}$ 中）。

    $$(A \odot B)_{11} = \max(3+1, 2+2) = \max(4, 4) = 4,$$
    $$(A \odot B)_{12} = \max(3+0, 2+3) = \max(3, 5) = 5,$$
    $$(A \odot B)_{21} = \max(1+1, 4+2) = \max(2, 6) = 6,$$
    $$(A \odot B)_{22} = \max(1+0, 4+3) = \max(1, 7) = 7.$$

    因此 $A \odot B = \begin{pmatrix}4 & 5\\6 & 7\end{pmatrix}$。

!!! definition "定义 59A.6 (热带矩阵幂与加权有向图)"
    $n \times n$ 热带矩阵 $A$ 的热带 $k$ 次幂定义为

    $$A^{\odot k} = \underbrace{A \odot A \odot \cdots \odot A}_{k},$$

    $(A^{\odot k})_{ij} = \max_{i_1, \ldots, i_{k-1}}(a_{i,i_1} + a_{i_1,i_2} + \cdots + a_{i_{k-1},j})$。

    **图论解释**：将 $A$ 视为有向加权图 $G(A)$ 的邻接矩阵（顶点集 $\{1,\ldots,n\}$，边 $(i,j)$ 权重 $a_{ij}$，$a_{ij} = -\infty$ 表示无边）。则 $(A^{\odot k})_{ij}$ 是从 $i$ 到 $j$ 的所有**恰好 $k$ 步**路径中权重之和的最大值。

    这一对应关系是热带矩阵论与组合优化之间桥梁的基石。

---

## 59A.4 热带行列式与最优指派

<div class="context-flow" markdown>

**核心问题**：热带行列式如何定义？它与经典行列式和指派问题有什么关系？

</div>

!!! definition "定义 59A.7 (热带行列式)"
    $n \times n$ 热带矩阵 $A$ 的**热带行列式**定义为

    $$\mathrm{tdet}(A) = \bigoplus_{\sigma \in S_n} \bigodot_{i=1}^n a_{i,\sigma(i)} = \max_{\sigma \in S_n} \sum_{i=1}^n a_{i,\sigma(i)},$$

    其中 $S_n$ 是 $n$ 元置换群。

!!! theorem "定理 59A.3 (热带行列式 = 热带永久式)"
    在经典代数中，行列式和永久式不同（行列式含符号因子 $\mathrm{sgn}(\sigma)$）。但在热带半环中，

    $$\mathrm{tdet}(A) = \mathrm{tperm}(A) = \max_{\sigma \in S_n}\sum_{i=1}^n a_{i,\sigma(i)}.$$

    原因是热带加法 $\oplus = \max$ 是幂等的（$a \oplus a = a$），因此符号因子无法在热带半环中定义（$-\infty$ 是加法单位元，不存在"$-a$"），行列式与永久式的区别消失。

??? proof "证明"
    经典行列式 $\det(A) = \sum_{\sigma} \mathrm{sgn}(\sigma)\prod_i a_{i,\sigma(i)}$。热带化后 $\sum \to \max$，$\prod \to \sum$，$\mathrm{sgn}(\sigma) \in \{+1, -1\}$。

    但 $\mathrm{sgn}(\sigma)$ 在热带半环中无意义——热带加法没有逆元，不存在"$-1$"这个加法逆元素。更精确地说，热带化 $+1$ 和 $-1$ 时，$+1$ 对应热带乘法单位元 $0$，而 $-1$ 在热带半环中没有对应物。

    因此热带行列式自然地定义为不带符号的版本，即永久式的热带化：

    $$\mathrm{tdet}(A) = \max_{\sigma \in S_n} \sum_{i=1}^n a_{i,\sigma(i)}.$$

    $\blacksquare$

!!! theorem "定理 59A.4 (热带行列式与最优指派问题)"
    热带行列式的计算等价于**最优指派问题**（assignment problem）：

    $$\mathrm{tdet}(A) = \max_{\sigma \in S_n}\sum_{i=1}^n a_{i,\sigma(i)}$$

    是将 $n$ 个工人分配到 $n$ 个任务的最大权重完美匹配。

    最优指派问题可以用 **Hungarian 算法**（匈牙利算法）在 $O(n^3)$ 时间内求解，因此热带行列式可在多项式时间内计算——这与经典永久式（#P 困难）形成鲜明对比。

??? proof "证明"
    热带行列式的定义

    $$\mathrm{tdet}(A) = \max_{\sigma \in S_n} \sum_{i=1}^n a_{i,\sigma(i)}$$

    恰好是二部图上最大权重完美匹配问题的形式化：一部分顶点为行 $\{1,\ldots,n\}$，另一部分为列 $\{1,\ldots,n\}$，行 $i$ 到列 $j$ 的边权为 $a_{ij}$，置换 $\sigma$ 对应完美匹配 $i \mapsto \sigma(i)$。

    Hungarian 算法通过交替维护对偶可行解和原始匹配来在 $O(n^3)$ 时间内求解。关键在于该问题可以表述为线性规划：

    $$\max \sum_{i,j} a_{ij} x_{ij} \quad \text{s.t.} \quad \sum_j x_{ij} = 1,\; \sum_i x_{ij} = 1,\; x_{ij} \geq 0,$$

    其约束矩阵是全幺模的（totally unimodular），因此线性松弛的最优解自动为整数解（即置换矩阵），且可用组合方法高效求解。

    相比之下，经典永久式 $\mathrm{perm}(A) = \sum_{\sigma} \prod_i a_{i,\sigma(i)}$ 是 #P 完全问题（Valiant, 1979），不存在已知的多项式时间精确算法。热带化将"求和"变为"取最大"，将困难的计数问题转化为可处理的优化问题。$\blacksquare$

!!! example "例 59A.3"
    设 $A = \begin{pmatrix}5 & 3 & 2\\4 & 6 & 1\\3 & 2 & 7\end{pmatrix}$。

    $\mathrm{tdet}(A) = \max_\sigma \sum_i a_{i,\sigma(i)}$：

    - $\sigma = \mathrm{id} = (1,2,3)$：$5+6+7 = 18$
    - $\sigma = (1,3,2)$：$5+1+2 = 8$
    - $\sigma = (2,1,3)$：$3+4+7 = 14$
    - $\sigma = (2,3,1)$：$3+1+3 = 7$
    - $\sigma = (3,1,2)$：$2+4+2 = 8$
    - $\sigma = (3,2,1)$：$2+6+3 = 11$

    因此 $\mathrm{tdet}(A) = 18$，由恒等置换达到。

---

## 59A.5 热带奇异性

!!! definition "定义 59A.8 (热带奇异性)"
    热带矩阵 $A$ 称为**热带奇异的**（tropically singular），如果热带行列式的最大值由**两个或更多**置换同时达到。否则 $A$ 称为**热带非奇异的**（tropically nonsingular）。

    形式化地：设 $S^* = \{\sigma \in S_n : \sum_i a_{i,\sigma(i)} = \mathrm{tdet}(A)\}$。

    - $|S^*| = 1$：$A$ 热带非奇异。
    - $|S^*| \geq 2$：$A$ 热带奇异。

    在例 59A.3 中，最大值仅由恒等置换达到，因此 $A$ 是热带非奇异的。

!!! example "例 59A.4"
    考虑 $B = \begin{pmatrix}3 & 2\\1 & 4\end{pmatrix}$。

    - $\sigma = \mathrm{id}$：$3 + 4 = 7$。
    - $\sigma = (12)$：$2 + 1 = 3$。

    最大值仅由 $\mathrm{id}$ 达到，$B$ 是热带非奇异的。

    现在考虑 $C = \begin{pmatrix}3 & 2\\2 & 1\end{pmatrix}$。

    - $\sigma = \mathrm{id}$：$3 + 1 = 4$。
    - $\sigma = (12)$：$2 + 2 = 4$。

    两个置换同时达到最大值 $4$，因此 $C$ 是热带奇异的。

    热带奇异性的图论解释：在二部图的最大权重匹配问题中，热带奇异意味着最优匹配不唯一。

---

## 59A.6 热带特征值

<div class="context-flow" markdown>

**核心问题**：如何定义和计算热带矩阵的特征值？热带特征值与有向图的圈结构有何关系？

</div>

!!! definition "定义 59A.9 (热带特征值问题)"
    设 $A \in \mathbb{T}^{n \times n}$。标量 $\lambda \in \mathbb{T}$ 和非零向量 $x \in \mathbb{T}^n$（$x \neq (-\infty, \ldots, -\infty)$）称为 $A$ 的**热带特征值**和**热带特征向量**，如果

    $$A \odot x = \lambda \odot x,$$

    即对所有 $i = 1, \ldots, n$：

    $$\max_{j=1,\ldots,n}(a_{ij} + x_j) = \lambda + x_i.$$

!!! definition "定义 59A.10 (有向图的最大环均值)"
    设 $G(A)$ 是 $A$ 的**关联有向图**：顶点集 $\{1, \ldots, n\}$，边 $(i,j)$ 存在当且仅当 $a_{ij} \neq -\infty$，权重为 $a_{ij}$。

    有向图中一个环 $\gamma = (i_1, i_2, \ldots, i_k, i_1)$ 的**环均值**（cycle mean）为

    $$\mu(\gamma) = \frac{a_{i_1 i_2} + a_{i_2 i_3} + \cdots + a_{i_k i_1}}{k}.$$

    **最大环均值**（maximum cycle mean）为 $\lambda^*(A) = \max_\gamma \mu(\gamma)$，最大取遍 $G(A)$ 的所有有向环。

---

## 59A.7 热带 Perron-Frobenius 定理

!!! theorem "定理 59A.5 (热带 Perron-Frobenius 定理)"
    设 $A \in \mathbb{T}^{n \times n}$，其关联有向图 $G(A)$ 是**强连通的**（从任意顶点到任意其他顶点都存在有向路径）。则：

    (a) $A$ 有唯一的热带特征值 $\lambda = \lambda^*(A)$（最大环均值）。

    (b) 热带特征向量 $x$ 存在（但不一定唯一，在热带投影空间中是凸多面体）。

    (c) 对任意有限初始向量 $x_0$，序列 $A^{\odot k} \odot x_0$ 的增长率渐近为 $\lambda$。

??? proof "证明"
    我们详细展开这一基本定理的证明。

    **第一步：标准化。** 不失一般性，令 $B = A \ominus \lambda$，其中 $(B)_{ij} = a_{ij} - \lambda$（逐元素减去 $\lambda$，普通减法）。则 $B$ 的最大环均值为 $0$，且 $\lambda$ 是 $A$ 的特征值当且仅当 $0$ 是 $B$ 的特征值。因此只需证明：当 $\lambda^*(A) = 0$ 时，$0$ 是唯一特征值。

    **第二步：存在性。** 假设 $\lambda^*(A) = 0$。定义矩阵

    $$A^+ = A \oplus A^{\odot 2} \oplus \cdots \oplus A^{\odot n} = \bigoplus_{k=1}^n A^{\odot k}.$$

    由于所有环均值 $\leq 0$，最大权重路径不需要绕正环（不存在），因此 $A^+$ 有限。$(A^+)_{ij}$ 是从 $i$ 到 $j$ 的所有长度 $\leq n$ 的路径中权重之和的最大值。

    取临界环上的某个顶点 $c$（即存在经过 $c$ 的环，其环均值恰为 $0$）。定义 $x_i = (A^+)_{ic}$，即从 $i$ 到 $c$ 的最大权重路径。我们证明 $x$ 满足 $A \odot x = x$（即 $\lambda = 0$ 的特征方程）。

    对任意 $i$，$(A \odot x)_i = \max_j(a_{ij} + x_j) = \max_j(a_{ij} + (A^+)_{jc})$。由于 $A^+ = A \oplus A \odot A^+$（截断到有限步时的递推关系），有

    $$(A \odot A^+)_{ic} = \max_j(a_{ij} + (A^+)_{jc}) \leq (A^+)_{ic} = x_i,$$

    其中不等式来自于最大环均值 $\leq 0$（多走一步不会改善）。另一方面，$x_i = (A^+)_{ic} \leq \max_j(a_{ij} + (A^+)_{jc})$，因为 $A^+$ 中从 $i$ 到 $c$ 的最优路径的第一步必须走到某个 $j$。因此等式成立：$(A \odot x)_i = x_i$。

    **第三步：唯一性。** 假设 $\lambda' \neq 0$ 也是特征值，$A \odot y = \lambda' \odot y$ 对某非零向量 $y$。

    **情形 1：$\lambda' > 0$。** 令 $B = A - \lambda'$（逐元素减去 $\lambda'$），则 $B$ 的所有环均值 $\leq -\lambda' < 0$（因为 $A$ 的最大环均值为 $0$）。但 $B \odot y = y$，即 $\max_j(b_{ij} + y_j) = y_i$。

    选取使 $y_i$ 最大的分量 $i^* = \arg\max_i y_i$（存在，因为 $y$ 非零有限）。由 $B \odot y = y$，存在 $j$ 使得 $b_{i^*j} + y_j = y_{i^*}$。由 $y_{j} \leq y_{i^*}$，得 $b_{i^*j} \geq 0$。

    从 $j$ 出发继续此过程：存在 $k$ 使得 $b_{jk} + y_k = y_j$，因此 $b_{jk} \geq y_j - y_k \geq 0 - (y_{i^*} - y_j) \geq -(y_{i^*} - y_j)$。

    由于图是强连通的且顶点有限，这条链最终形成环。设环为 $i_1, i_2, \ldots, i_\ell, i_1$，则沿环有

    $$b_{i_1 i_2} + b_{i_2 i_3} + \cdots + b_{i_\ell i_1} = (y_{i_1} - y_{i_2}) + (y_{i_2} - y_{i_3}) + \cdots + (y_{i_\ell} - y_{i_1}) = 0.$$

    这意味着 $B$ 有一个环的权重之和为 $0$，即环均值为 $0/\ell = 0$。但我们已知 $B$ 的所有环均值 $\leq -\lambda' < 0$，矛盾。

    **情形 2：$\lambda' < 0$。** 令 $B = A - \lambda'$（逐元素减 $\lambda'$），此时 $B$ 的最大环均值为 $-\lambda' > 0$。$B \odot y = y$。

    取 $y_{i^*} = \min_i y_i$（最小分量）。由 $\max_j(b_{i^*j} + y_j) = y_{i^*}$，对所有 $j$ 有 $b_{i^*j} + y_j \leq y_{i^*}$，即 $b_{i^*j} \leq y_{i^*} - y_j \leq 0$。

    由强连通性，从 $i^*$ 出发可达所有顶点。沿类似的论证可以证明所有 $b_{ij} \leq 0$，这意味着 $B$ 的所有环权重 $\leq 0$，即 $B$ 的最大环均值 $\leq 0$，与 $-\lambda' > 0$ 矛盾。

    因此 $\lambda^*(A) = 0$ 是唯一的热带特征值。对一般的 $A$（$\lambda^*(A) = \lambda$），通过标准化回去，$\lambda$ 是唯一的热带特征值。

    **第四步 (b)：特征向量空间的结构。** 特征向量集合

    $$V(\lambda) = \{x \in \mathbb{T}^n \setminus \{(-\infty)^n\} : A \odot x = \lambda \odot x\}$$

    在热带投影空间 $\mathbb{TP}^{n-1} = (\mathbb{T}^n \setminus \{(-\infty)^n\})/\sim$（其中 $x \sim y \iff x = \alpha \odot y$ 对某 $\alpha \in \mathbb{R}$）中是一个热带凸多面体。其维度和结构由临界图决定（见定义 59A.11）。

    **(c)** 渐近增长率的证明：$x(k) = A^{\odot k} \odot x_0$。$(A^{\odot k})_{ij}$ 是从 $j$ 到 $i$ 的 $k$ 步最大权重路径。当 $k$ 充分大时，最优 $k$ 步路径的主要贡献来自反复绕临界环。每绕一次长度为 $\ell$ 的临界环，权重增加 $\ell \cdot \lambda$。因此 $(A^{\odot k})_{ij} \sim k\lambda + O(1)$，从而 $x_i(k)/k \to \lambda$。$\blacksquare$

---

## 59A.8 临界图

!!! definition "定义 59A.11 (临界图)"
    矩阵 $A$ 的**临界图** $G_c(A)$ 是由达到最大环均值的所有环及其上的顶点和边构成的子图。形式化地：

    $$G_c(A) = \bigcup_{\gamma : \mu(\gamma) = \lambda^*(A)} \gamma.$$

    临界图的顶点称为**临界顶点**，边称为**临界边**。

    临界图在调度理论中有重要意义：临界图上的环对应系统的"瓶颈"环路，改善系统性能的关键在于改善临界环。

!!! example "例 59A.5"
    设 $A = \begin{pmatrix}2 & 5\\3 & 1\end{pmatrix}$。

    $G(A)$ 有四条边：$(1,1)$ 权 2，$(1,2)$ 权 5，$(2,1)$ 权 3，$(2,2)$ 权 1。

    环及其均值：

    - 自环 $(1,1)$：均值 $2/1 = 2$。
    - 自环 $(2,2)$：均值 $1/1 = 1$。
    - 环 $(1,2,1)$：均值 $(5+3)/2 = 4$。

    最大环均值 $\lambda = 4$。临界图 $G_c(A)$ 由环 $(1,2,1)$ 构成，包含顶点 $\{1,2\}$ 和边 $\{(1,2), (2,1)\}$。

    验证特征方程：取 $x = (0, -1)$，

    $$\max(2+0, 5+(-1)) = \max(2, 4) = 4 = 4 + 0, \quad \checkmark$$
    $$\max(3+0, 1+(-1)) = \max(3, 0) = 3 = 4 + (-1) = 3. \quad \checkmark$$

    因此 $\lambda = 4$，$x = (0, -1)$ 是特征对。

---

## 59A.9 周期性定理

<div class="context-flow" markdown>

**核心问题**：热带矩阵幂序列 $A, A^{\odot 2}, A^{\odot 3}, \ldots$ 最终呈现怎样的周期行为？

</div>

!!! theorem "定理 59A.6 (周期性定理 / Cyclicity Theorem)"
    设 $A \in \mathbb{T}^{n \times n}$，$G(A)$ 强连通，$\lambda = \lambda^*(A)$ 为最大环均值。则存在正整数 $K$（暂态长度）和 $\gamma$（周期），使得对所有 $k \geq K$：

    $$A^{\odot(k+\gamma)} = \lambda^{\odot \gamma} \odot A^{\odot k},$$

    即 $(A^{\odot(k+\gamma)})_{ij} = \gamma\lambda + (A^{\odot k})_{ij}$ 对所有 $i, j$。

    周期 $\gamma$ 等于临界图 $G_c(A)$ 中所有环长度的**最大公约数**：

    $$\gamma = \gcd\{|\gamma_c| : \gamma_c \text{ 是 } G_c(A) \text{ 中的环}\}.$$

    特别地，如果临界图包含长度互素的环，则 $\gamma = 1$（即最终完全周期行为退化为逐步线性增长）。

??? proof "证明"
    **第一步：归约。** 不妨设 $\lambda = 0$（否则用 $B_{ij} = a_{ij} - \lambda$ 替换，$B$ 的最大环均值为 $0$，且 $A^{\odot k} = \lambda^{\odot k} \odot B^{\odot k}$，即 $(A^{\odot k})_{ij} = k\lambda + (B^{\odot k})_{ij}$）。

    在 $\lambda = 0$ 的情形下，需证 $B^{\odot(k+\gamma)} = B^{\odot k}$，即矩阵幂最终周期。

    **第二步：路径分析。** $(B^{\odot k})_{ij}$ 是从 $i$ 到 $j$ 的 $k$ 步最大权重路径。由于所有环均值 $\leq 0$，最优路径中每绕一个非临界环都会损失权重，而绕临界环则权重不变。因此，对充分大的 $k$，最优 $k$ 步路径的形式是：一段从 $i$ 到临界图的路径，在临界图上绕若干圈，再一段从临界图到 $j$ 的路径。

    **第三步：周期性。** 临界图上的绕圈行为具有周期性。设临界图中的环长度集合为 $\{l_1, l_2, \ldots\}$，$\gamma = \gcd(l_1, l_2, \ldots)$。由数论中的基本结果，对充分大的整数 $m$，$m$ 能表示为 $\{l_i\}$ 的非负整数线性组合当且仅当 $\gamma | m$。

    当 $k$ 充分大时，从 $i$ 到 $j$ 的 $k$ 步最优路径主要由临界图上的 $k - O(1)$ 步构成。增加 $\gamma$ 步等于在临界图上多绕一些环（总长度恰为 $\gamma$），而临界环的权重之和为 $0$（因为环均值为 $0$），因此权重不变：$(B^{\odot(k+\gamma)})_{ij} = (B^{\odot k})_{ij}$。

    **第四步：暂态长度。** 暂态长度 $K$ 的上界可以用图的直径和临界图的结构来估计。一个经典上界是 $K \leq n^2$（Schwarz, 1970）。更精确的界与 Wielandt 数有关。$\blacksquare$

!!! example "例 59A.6"
    设 $A = \begin{pmatrix}0 & 1 & -\infty\\-\infty & 0 & 2\\3 & -\infty & 0\end{pmatrix}$。

    环 $(1,2,3,1)$：权重 $1 + 2 + 3 = 6$，均值 $6/3 = 2$。
    自环 $(1,1), (2,2), (3,3)$：均值 $0, 0, 0$。

    最大环均值 $\lambda = 2$。临界图仅含环 $(1,2,3,1)$，长度为 $3$。

    周期 $\gamma = \gcd(3) = 3$。这意味着序列 $A^{\odot k}$（标准化后）以周期 $3$ 最终重复。

    实际上，如果存在自环且自环也在临界图中（均值 $= \lambda$），则 $\gamma = \gcd(3, 1) = 1$。但本例中自环均值 $0 < 2 = \lambda$，不在临界图中。

---

## 59A.10 Karp 算法

<div class="context-flow" markdown>

**核心问题**：如何高效计算最大环均值（即热带特征值）？

</div>

!!! theorem "定理 59A.7 (Karp 算法)"
    设 $A \in \mathbb{T}^{n \times n}$，$G(A)$ 强连通。最大环均值可以通过以下公式计算：

    $$\lambda^*(A) = \max_{1 \leq i \leq n} \min_{0 \leq k \leq n-1} \frac{(A^{\odot n})_{ii} - (A^{\odot k})_{ii}}{n - k}.$$

    更一般地，对任意固定顶点 $j$：

    $$\lambda^*(A) = \max_{1 \leq i \leq n} \min_{0 \leq k \leq n-1} \frac{(A^{\odot n})_{ij} - (A^{\odot k})_{ij}}{n - k}.$$

    该算法的时间复杂度为 $O(n^3)$（计算所有 $A^{\odot k}$ 需 $O(n^3)$，然后取极值需 $O(n^2)$）。

??? proof "证明"
    **核心思想**：对任意 $i, j$，$(A^{\odot m})_{ij}$ 是从 $j$ 到 $i$ 的 $m$ 步最大权重路径。对充分大的 $m$，最优路径的增长率趋于 $\lambda$。

    设 $F_m(i) = (A^{\odot m})_{ij}$（固定 $j$）。考虑 $n$ 步路径。$n$ 步路径至少经过 $n+1$ 个顶点（含重复），因此包含至少一个环。

    **上界**：对任意环 $\gamma$，$\mu(\gamma) \leq \lambda$。$n$ 步路径可以分解为一条无环路径（$\leq n-1$ 步）和若干环。若无环部分为 $k$ 步（$0 \leq k \leq n-1$），权重 $\leq F_k(i)$，环的总步数为 $n-k$，总权重 $\leq (n-k)\lambda$。因此

    $$F_n(i) \leq F_k(i) + (n-k)\lambda, \quad \forall k \in \{0, \ldots, n-1\}.$$

    取对 $k$ 的最大：$\lambda \geq \min_k \frac{F_n(i) - F_k(i)}{n-k}$。再取对 $i$ 的最大：$\lambda \geq \max_i \min_k \frac{F_n(i) - F_k(i)}{n-k}$。

    **下界**：取临界环上的顶点 $i^*$，环长 $\ell$，环均值 $\lambda$。$n$ 步从 $j$ 到 $i^*$ 的最优路径，先走 $k^*$ 步到 $i^*$，然后绕环 $\lfloor(n-k^*)/\ell\rfloor$ 次，剩余步数走附近。精确分析可以证明

    $$\min_k \frac{F_n(i^*) - F_k(i^*)}{n-k} \geq \lambda.$$

    综合得

    $$\lambda = \max_i \min_k \frac{F_n(i) - F_k(i)}{n-k}.$$

    $\blacksquare$

!!! example "例 59A.7"
    对例 59A.5 中的 $A = \begin{pmatrix}2 & 5\\3 & 1\end{pmatrix}$。

    $A^{\odot 0} = I^{\oplus} = \begin{pmatrix}0 & -\infty\\-\infty & 0\end{pmatrix}$。

    $A^{\odot 1} = A = \begin{pmatrix}2 & 5\\3 & 1\end{pmatrix}$。

    $A^{\odot 2}$：$(A^{\odot 2})_{11} = \max(2+2, 5+3) = 8$，$(A^{\odot 2})_{12} = \max(2+5, 5+1) = 7$，$(A^{\odot 2})_{21} = \max(3+2, 1+3) = 5$，$(A^{\odot 2})_{22} = \max(3+5, 1+1) = 8$。

    对 $j = 1$（取第一列）：$F_0(1) = 0, F_0(2) = -\infty; F_1(1) = 2, F_1(2) = 3; F_2(1) = 8, F_2(2) = 5$。

    对 $i = 1$：$\min_k \frac{F_2(1)-F_k(1)}{2-k} = \min\left(\frac{8-0}{2}, \frac{8-2}{1}\right) = \min(4, 6) = 4$。

    对 $i = 2$：$\min_k \frac{F_2(2)-F_k(2)}{2-k} = \min\left(\frac{5-(-\infty)}{2}, \frac{5-3}{1}\right) = \min(+\infty, 2) = 2$。

    $\lambda = \max(4, 2) = 4$。与直接计算一致。

---

## 59A.11 Kleene 星与收敛

<div class="context-flow" markdown>

**核心问题**：热带矩阵幂与最短/最长路径问题有什么关系？Kleene 星算子如何给出全源最短路径？

</div>

!!! definition "定义 59A.12 (Kleene 星)"
    设 $A \in \mathbb{T}^{n \times n}$。$A$ 的 **Kleene 星**（或闭包）定义为

    $$A^* = I^{\oplus} \oplus A \oplus A^{\odot 2} \oplus A^{\odot 3} \oplus \cdots = \bigoplus_{k=0}^{\infty} A^{\odot k},$$

    即 $(A^*)_{ij} = \sup_{k \geq 0} (A^{\odot k})_{ij}$。这在热带半环中是从 $i$ 到 $j$ 的所有路径（包括零步路径 $i = j$）中权重之和的最大值。

!!! theorem "定理 59A.8 (Kleene 星的收敛性)"
    $A^*$ 在 $\mathbb{T}$ 中有限当且仅当 $A$ 的所有有向环的权重之和 $\leq 0$（在 max-plus 意义下，即 $\lambda^*(A) \leq 0$）。此时

    $$A^* = I^{\oplus} \oplus A \oplus A^{\odot 2} \oplus \cdots \oplus A^{\odot (n-1)}.$$

??? proof "证明"
    如果存在权重之和为正的环 $\gamma$，则沿 $\gamma$ 无限绕行可使路径权重趋于 $+\infty$，$A^*$ 发散。

    反之，如果所有环权重 $\leq 0$，则最优路径不需要经过重复顶点（重复会产生非正环，不改善或恶化权重）。$n$ 个顶点的图中无重复顶点的路径最多 $n-1$ 步，因此 $A^* = \bigoplus_{k=0}^{n-1} A^{\odot k}$。$\blacksquare$

---

## 59A.12 Floyd-Warshall 与 Bellman-Ford 的热带解释

!!! theorem "定理 59A.9 (Floyd-Warshall 与热带矩阵幂)"
    设 $D \in \mathbb{T}_{\min}^{n \times n}$ 是有向加权图的邻接矩阵（使用 min-plus 半环）。则：

    (a) $D^{\odot k}_{\min}$ 的 $(i,j)$ 元素是从 $i$ 到 $j$ 的恰好 $k$ 步最短路径长度。

    (b) $D^*_{\min} = \bigoplus_{k=0}^{n-1} D^{\odot k}_{\min}$ 的 $(i,j)$ 元素是从 $i$ 到 $j$ 的最短路径长度（假设无负环）。

    (c) **Floyd-Warshall 算法**可以视为逐步消元计算 $D^*_{\min}$ 的过程。

    具体地，Floyd-Warshall 的更新

    $$d_{ij}^{(k)} = \min(d_{ij}^{(k-1)},\, d_{ik}^{(k-1)} + d_{kj}^{(k-1)})$$

    恰好对应热带高斯消元：用中间顶点 $k$ "消去"间接路径，逐步构建 $D^*_{\min}$。

!!! example "例 59A.8"
    **最长路径的热带矩阵计算。** 考虑有向图，邻接矩阵（max-plus）

    $$A = \begin{pmatrix}-\infty & 3 & -\infty\\-\infty & -\infty & 2\\-\infty & -\infty & -\infty\end{pmatrix}.$$

    表示边 $(1,2)$ 权 3，$(2,3)$ 权 2。

    $A^{\odot 2}$：$(A^{\odot 2})_{13} = \max(a_{11}+a_{13}, a_{12}+a_{23}, a_{13}+a_{33}) = \max(-\infty, 3+2, -\infty) = 5$。

    $A^*$：$(A^*)_{13} = \max(0, -\infty, 5) = 5$（从 1 到 3 的最长路径权重为 $3 + 2 = 5$）。

!!! theorem "定理 59A.10 (Bellman-Ford 方程的热带表述)"
    设 $d = (d_1, \ldots, d_n)^\top$ 是从源点 $s$ 到各点的最短距离向量（min-plus 约定）。Bellman-Ford 方程

    $$d_i = \min\!\bigl(0 \cdot [i=s],\, \min_{j: (j,i) \in E}(d_j + w_{ji})\bigr)$$

    可以写为热带线性方程

    $$d = D_{\min} \odot d \oplus_{\min} e_s,$$

    其中 $e_s$ 是第 $s$ 个分量为 $0$、其余为 $+\infty$ 的向量。迭代 $d^{(k+1)} = D_{\min} \odot d^{(k)} \oplus_{\min} e_s$ 就是 Bellman-Ford 松弛过程。

    从热带代数视角看，Bellman-Ford 是在 min-plus 半环上解线性不动点方程的 Jacobi 迭代。

---

## 59A.13 热带线性方程组

<div class="context-flow" markdown>

**核心问题**：如何求解热带线性方程组 $A \odot x = b$？何时有解？

</div>

!!! definition "定义 59A.13 (热带线性方程组)"
    热带线性方程组 $A \odot x = b$，其中 $A \in \mathbb{T}^{m \times n}$，$b \in \mathbb{T}^m$，是指

    $$\max_{1 \leq j \leq n}(a_{ij} + x_j) = b_i, \quad i = 1, \ldots, m.$$

    这是一个**分段线性**方程组。注意，由于热带加法的幂等性，这与经典线性方程组有根本不同。

!!! theorem "定理 59A.11 (热带线性方程组的可解性)"
    热带方程组 $A \odot x = b$ 的可解性问题与经典情形有根本差异：

    (a) **不等式形式**：$A \odot x \leq b$（即 $a_{ij} + x_j \leq b_i$，$\forall i,j$）总有解，最大解为

    $$x_j^* = \min_{1 \leq i \leq m}(b_i - a_{ij}).$$

    (b) **等式形式**：$A \odot x = b$ 有解当且仅当 $x^* = (-A^\top \backslash b)$（定义见残差理论）满足 $A \odot x^* = b$。

    (c) 如果有解，$x^*$ 是**最大解**（分量最大的解）。

??? proof "证明"
    **(a)** $A \odot x \leq b$ 等价于 $a_{ij} + x_j \leq b_i$ 对所有 $i, j$，即 $x_j \leq b_i - a_{ij}$ 对所有 $i$。因此 $x_j \leq \min_i(b_i - a_{ij}) = x_j^*$。$x^*$ 自身满足不等式是直接验证。

    **(b)** 考虑候选解 $x_j^* = \min_i(b_i - a_{ij})$。由 (a)，$A \odot x^* \leq b$。$A \odot x = b$ 有解当且仅当 $A \odot x^* \geq b$，即 $\max_j(a_{ij} + x_j^*) \geq b_i$ 对所有 $i$。由于 $A \odot x^* \leq b$ 已成立，条件简化为 $A \odot x^* = b$。

    **(c)** 如果 $x$ 是 $A \odot x = b$ 的任何解，则 $A \odot x \leq b$（因为等式意味着不等式），由 (a) 的论证得 $x_j \leq x_j^*$。$\blacksquare$

!!! example "例 59A.9"
    热带方程 $A \odot x = b$，其中

    $$A = \begin{pmatrix}3 & 1\\2 & 4\end{pmatrix}, \quad b = \begin{pmatrix}5\\7\end{pmatrix}.$$

    候选最大解：$x_1^* = \min(5-3, 7-2) = \min(2, 5) = 2$，$x_2^* = \min(5-1, 7-4) = \min(4, 3) = 3$。

    验证：$A \odot x^* = \begin{pmatrix}\max(3+2, 1+3)\\\max(2+2, 4+3)\end{pmatrix} = \begin{pmatrix}\max(5,4)\\\max(4,7)\end{pmatrix} = \begin{pmatrix}5\\7\end{pmatrix} = b$。

    成功。$x^* = (2, 3)^\top$ 是最大解。

---

## 59A.14 热带 Cramer 法则

!!! theorem "定理 59A.12 (热带 Cramer 法则)"
    设 $A \in \mathbb{T}^{n \times n}$ 是热带非奇异的，$b \in \mathbb{T}^n$。考虑热带线性方程组 $A \odot x = b$。设 $A_j(b)$ 是将 $A$ 的第 $j$ 列替换为 $b$ 得到的矩阵。则：

    (a) 如果方程组有解，则

    $$x_j = \mathrm{tdet}(A_j(b)) \odot (\mathrm{tdet}(A))^{\odot(-1)} = \mathrm{tdet}(A_j(b)) - \mathrm{tdet}(A).$$

    注意 $(\mathrm{tdet}(A))^{\odot(-1)} = -\mathrm{tdet}(A)$（热带乘法逆就是取负）。

    (b) 更一般地，即使方程组无精确解，上式给出的 $x$ 是**最佳近似解**（在某种热带范数意义下）。

??? proof "证明"
    热带 Cramer 法则是经典 Cramer 法则的热带类比。经典 Cramer 法则 $x_j = \det(A_j(b))/\det(A)$，热带化后

    $$x_j = \mathrm{tdet}(A_j(b)) \oslash \mathrm{tdet}(A) = \mathrm{tdet}(A_j(b)) - \mathrm{tdet}(A).$$

    其中 $\oslash$ 是热带除法（即普通减法）。

    证明需要展开 $\mathrm{tdet}(A_j(b))$。设 $\sigma^*$ 是达到 $\mathrm{tdet}(A)$ 的最优置换（唯一，因为 $A$ 非奇异）。在 $A_j(b)$ 中，最优置换将第 $j$ 列的 $b$ 贡献包含进来。当 $b = A \odot x$ 时，$b_i = \max_k(a_{ik} + x_k)$。

    对热带非奇异矩阵，最优置换的唯一性保证了展开式的精确对应。详细的组合论证可以验证公式在有解情况下成立。$\blacksquare$

!!! example "例 59A.10"
    设 $A = \begin{pmatrix}5 & 3\\2 & 6\end{pmatrix}$，$b = \begin{pmatrix}8\\10\end{pmatrix}$。

    $\mathrm{tdet}(A) = \max(5+6, 3+2) = \max(11, 5) = 11$（由恒等置换达到，唯一，故非奇异）。

    $A_1(b) = \begin{pmatrix}8 & 3\\10 & 6\end{pmatrix}$，$\mathrm{tdet}(A_1(b)) = \max(8+6, 3+10) = \max(14, 13) = 14$。

    $A_2(b) = \begin{pmatrix}5 & 8\\2 & 10\end{pmatrix}$，$\mathrm{tdet}(A_2(b)) = \max(5+10, 8+2) = \max(15, 10) = 15$。

    $x_1 = 14 - 11 = 3$，$x_2 = 15 - 11 = 4$。

    验证：$A \odot x = \begin{pmatrix}\max(5+3, 3+4)\\\max(2+3, 6+4)\end{pmatrix} = \begin{pmatrix}\max(8,7)\\\max(5,10)\end{pmatrix} = \begin{pmatrix}8\\10\end{pmatrix} = b$。正确。

---

## 59A.15 残差理论

<div class="context-flow" markdown>

**核心问题**：热带半环中没有减法和除法（加法无逆元），如何进行"除法"运算？

</div>

!!! definition "定义 59A.14 (残差 / Residuation)"
    在热带半环中，虽然不存在加法逆元，但存在**残差**（residuation）运算，作为"除法"的替代。

    对 $a, b \in \mathbb{T}$，定义**左残差**

    $$a \backslash b = \sup\{x \in \mathbb{T} : a \odot x \leq b\} = b - a \quad (\text{当 } a \neq -\infty).$$

    对矩阵 $A \in \mathbb{T}^{m \times n}$ 和 $b \in \mathbb{T}^m$，定义

    $$(A \backslash b)_j = \min_{1 \leq i \leq m}(b_i - a_{ij}) = \bigwedge_{i=1}^m (a_{ij} \backslash b_i).$$

    这正是定理 59A.11 中热带不等式 $A \odot x \leq b$ 的最大解。

!!! theorem "定理 59A.13 (残差的性质)"
    设 $A \in \mathbb{T}^{m \times n}$，$b \in \mathbb{T}^m$，$x^* = A \backslash b$。则：

    (a) $A \odot x^* \leq b$。

    (b) 对任意 $x$，$A \odot x \leq b \implies x \leq x^*$（分量逐个比较）。

    (c) $A \odot (A \backslash b) \leq b \leq A \backslash (A \odot b)$（对向量 $b$）。

    (d) $A \odot (A \backslash (A \odot x)) = A \odot x$（投影性质）。

    这些性质表明残差给出了热带线性不等式的"最佳可能"解，类似于经典情形中的伪逆。

??? proof "证明"
    **(a)** $(A \odot x^*)_i = \max_j(a_{ij} + x_j^*) = \max_j(a_{ij} + \min_k(b_k - a_{kj})) \leq \max_j(a_{ij} + (b_i - a_{ij})) = b_i$。

    **(b)** 若 $A \odot x \leq b$，则 $a_{ij} + x_j \leq b_i$ 对所有 $i$，故 $x_j \leq \min_i(b_i - a_{ij}) = x_j^*$。

    **(c)** 第一个不等式就是 (a)。第二个：$(A \backslash (A \odot b))_j = \min_i((A \odot b)_i - a_{ij}) = \min_i(\max_k(a_{ik}+b_k) - a_{ij}) \geq \min_i(a_{ij}+b_j - a_{ij}) = b_j$。

    **(d)** 由 (c) 的对偶形式和等式条件的验证。$\blacksquare$

---

## 59A.16 热带秩

<div class="context-flow" markdown>

**核心问题**：热带矩阵的"秩"如何定义？为什么有多种不同的定义？

</div>

在经典线性代数中，矩阵的秩有多种等价定义（列秩 = 行秩 = 最大非零子式的阶数）。但在热带代数中，这些定义**不再等价**，产生了多种不同的"秩"概念。

!!! definition "定义 59A.15 (热带秩的多种定义)"
    设 $A \in \mathbb{T}^{m \times n}$。

    **(a) Barvinok 秩**（因子秩）：$\mathrm{rk}_B(A) = \min\{r : A = U \odot V,\, U \in \mathbb{T}^{m \times r},\, V \in \mathbb{T}^{r \times n}\}$，即 $A$ 的最小热带分解秩。

    **(b) Kapranov 秩**：$\mathrm{rk}_K(A)$ 是使得 $A$ 能被赋值映射 $\mathrm{val}: \mathbb{K}^{m \times n} \to \mathbb{T}^{m \times n}$ 从经典矩阵 $\tilde{A}$（秩 $\leq r$）提升的最小 $r$。即 $\mathrm{rk}_K(A) = \min\{\mathrm{rk}(\tilde{A}) : \mathrm{val}(\tilde{A}) = A\}$。

    **(c) 热带秩**（tropical rank）：$\mathrm{rk}_T(A) = $ 最大的 $r$ 使得 $A$ 有一个 $r \times r$ 热带非奇异子矩阵。

!!! theorem "定理 59A.14 (热带秩的不等式)"
    对任意 $A \in \mathbb{T}^{m \times n}$：

    $$\mathrm{rk}_T(A) \leq \mathrm{rk}_K(A) \leq \mathrm{rk}_B(A).$$

    这些不等式一般是**严格的**。特别地：

    (a) Barvinok 秩与 Kapranov 秩可以相差任意大。

    (b) 热带秩与 Kapranov 秩也可以不同。

    (c) Barvinok 秩的计算是 NP 困难的。

!!! example "例 59A.11"
    考虑 $3 \times 3$ 矩阵

    $$A = \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}.$$

    所有元素为 $0$（热带乘法单位元）。

    **热带秩**：$\mathrm{tdet}(A) = \max_\sigma(0+0+0) = 0$。所有 $6$ 个置换都给出相同的值 $0$，因此 $A$ 是热带奇异的。但任意 $2 \times 2$ 子矩阵的热带行列式 $= \max(0+0, 0+0) = 0$，也是奇异的（两个置换同时达到）。而 $1 \times 1$ 子矩阵（单元素 $0$）是非奇异的。所以 $\mathrm{rk}_T(A) = 1$。

    **Barvinok 秩**：$A = \begin{pmatrix}0\\0\\0\end{pmatrix} \odot (0, 0, 0) = \mathbf{0} \odot \mathbf{0}^\top$。所以 $\mathrm{rk}_B(A) = 1$。

    **Kapranov 秩**：需要找经典矩阵 $\tilde{A}$，使得每个元素的赋值为 $0$（即每个元素是单位量级的）。例如取 $\tilde{A} = \begin{pmatrix}1&1&1\\1&1&1\\1&1&1\end{pmatrix}$，经典秩为 $1$，$\mathrm{val}(\tilde{A}) = A$。所以 $\mathrm{rk}_K(A) = 1$。

    本例中三种秩相同。但对更复杂的矩阵，三者可以不同。

!!! example "例 59A.12"
    Develin-Santos-Sturmfels 的经典反例。考虑

    $$A = \begin{pmatrix}0 & 0 & 0\\0 & 0 & -1\\0 & -1 & 0\end{pmatrix}.$$

    可以验证 $\mathrm{rk}_T(A) = 2$（某 $2\times2$ 子矩阵非奇异，但 $3\times3$ 奇异）。然而 $\mathrm{rk}_K(A) = 3$（无法从经典秩 $\leq 2$ 的矩阵提升），因此 $\mathrm{rk}_T(A) < \mathrm{rk}_K(A)$。

---

## 59A.17 离散事件系统

<div class="context-flow" markdown>

**核心问题**：热带线性代数如何为调度问题和离散事件系统提供统一的数学框架？

</div>

!!! example "例 59A.13"
    **火车时刻表的热带模型。** 考虑一个简单的铁路网络，有两个站 $A$、$B$，火车在两站之间往返。设：

    - $x_1(k)$ = 第 $k$ 次从 $A$ 出发的时间。
    - $x_2(k)$ = 第 $k$ 次从 $B$ 出发的时间。
    - 从 $A$ 到 $B$ 需时 3（单位），从 $B$ 到 $A$ 需时 2。
    - 在 $A$ 站最少停留 1，在 $B$ 站最少停留 1。

    约束：$x_1(k+1) \geq x_2(k) + 2$（从 $B$ 返回后），$x_1(k+1) \geq x_1(k) + 1$（在 $A$ 站停留）。
    类似地，$x_2(k+1) \geq x_1(k) + 3$，$x_2(k+1) \geq x_2(k) + 1$。

    用热带记号：$x(k+1) = A \odot x(k)$，其中 $A = \begin{pmatrix}1 & 2\\3 & 1\end{pmatrix}$。

    热带特征值 $\lambda = \max(1/1, 1/1, (2+3)/2) = 5/2 = 2.5$。这意味着系统的渐近周期为 $2.5$：稳定运行后每 $2.5$ 个时间单位发一班车。

!!! theorem "定理 59A.15 (热带线性递推的渐近行为)"
    设 $A \in \mathbb{T}^{n \times n}$，$G(A)$ 强连通。考虑热带线性递推

    $$x(k+1) = A \odot x(k), \quad x(0) = x_0.$$

    则对任意初始向量 $x_0$（各分量有限），

    $$\lim_{k \to \infty} \frac{x_i(k)}{k} = \lambda \quad \forall\, i = 1, \ldots, n,$$

    其中 $\lambda = \lambda^*(A)$ 是最大环均值。更精确地，由周期性定理（定理 59A.6），存在 $K$ 和周期 $\gamma$ 使得对所有 $k \geq K$：

    $$x(k+\gamma) = \lambda^{\odot \gamma} \odot x(k) \quad (\text{即 } x_i(k+\gamma) = \gamma\lambda + x_i(k)).$$

??? proof "证明"
    核心思想：$x(k) = A^{\odot k} \odot x_0$。$(A^{\odot k})_{ij}$ 是从 $j$ 到 $i$ 的所有 $k$ 步路径中权重之和的最大值。

    当 $k$ 充分大时，最优 $k$ 步路径包含许多绕临界环的"循环"。每绕一次临界环（设长度为 $\ell$），路径权重增加 $\ell \lambda$。因此 $(A^{\odot k})_{ij} \approx k\lambda +$ 常数项。

    精确证明依赖周期性定理：对 $k \geq K$，$A^{\odot(k+\gamma)} = \lambda^{\odot\gamma} \odot A^{\odot k}$，因此

    $$x_i(k+\gamma) = (A^{\odot(k+\gamma)} \odot x_0)_i = (\lambda^{\odot\gamma} \odot A^{\odot k} \odot x_0)_i = \gamma\lambda + x_i(k).$$

    从而 $x_i(k)/k \to \lambda$。$\blacksquare$

!!! example "例 59A.14"
    **制造系统调度。** 一条装配线有 3 个工序，工序 $i$ 的处理时间为 $p_i$，工序间有先后约束。设

    $$A = \begin{pmatrix}p_1 & d_{12} & -\infty\\-\infty & p_2 & d_{23}\\d_{31} & -\infty & p_3\end{pmatrix},$$

    其中 $d_{ij}$ 是工序 $i$ 完成后到工序 $j$ 开始的延迟。最大环均值给出了系统的最小循环时间（瓶颈）。

    例如取 $p_1 = 3, p_2 = 2, p_3 = 4, d_{12} = 1, d_{23} = 1, d_{31} = 2$。

    环分析：
    - 自环 $(1,1)$：均值 $3$。
    - 自环 $(2,2)$：均值 $2$。
    - 自环 $(3,3)$：均值 $4$。
    - 环 $(1,2,3,1)$：权重 $1+1+2 = 4$，均值 $4/3 \approx 1.33$。

    最大环均值 $\lambda = 4$（工序 3 的处理时间是瓶颈）。要改善系统性能，必须缩短工序 3 的时间。

---

## 59A.18 习题

**基础计算**

!!! example "习题 59A.1"
    在热带半环中计算：

    (a) $(3 \oplus 7) \odot (2 \oplus 5)$。

    (b) $\begin{pmatrix}1 & 3\\4 & 2\end{pmatrix} \odot \begin{pmatrix}2 & 0\\1 & 3\end{pmatrix}$。

    (c) $A^{\odot 3}$，其中 $A = \begin{pmatrix}0 & 1\\2 & 0\end{pmatrix}$。

!!! example "习题 59A.2"
    计算下列矩阵的热带行列式，并判断是否热带奇异：

    (a) $\begin{pmatrix}4 & 1\\2 & 3\end{pmatrix}$。

    (b) $\begin{pmatrix}1 & 2 & 0\\3 & 0 & 1\\0 & 1 & 2\end{pmatrix}$。

!!! example "习题 59A.3"
    证明：在热带半环中，$a \oplus b = b$ 当且仅当 $a \leq b$（普通序）。由此推出热带加法给 $\mathbb{T}$ 赋予了自然的偏序结构。

**热带特征值**

!!! example "习题 59A.4"
    计算下列矩阵的热带特征值和一组热带特征向量：

    (a) $A = \begin{pmatrix}3 & 1\\4 & 2\end{pmatrix}$。

    (b) $A = \begin{pmatrix}0 & 2 & -\infty\\3 & 0 & 1\\-\infty & 4 & 0\end{pmatrix}$。

!!! example "习题 59A.5"
    证明：如果 $A$ 的关联有向图不是强连通的，$A$ 可能有多个热带特征值。构造一个 $3 \times 3$ 的例子。

!!! example "习题 59A.6"
    对例 59A.13 中的火车调度问题，如果从 $A$ 到 $B$ 的行车时间从 $3$ 改为 $4$，新的系统周期是什么？临界图如何变化？

**周期性与算法**

!!! example "习题 59A.7"
    设 $A = \begin{pmatrix}0 & 3 & -\infty\\-\infty & 0 & 2\\1 & -\infty & 0\end{pmatrix}$。

    (a) 计算 $A^{\odot k}$（$k = 1, 2, 3, 4, 5$）。

    (b) 找到暂态长度 $K$ 和周期 $\gamma$，验证周期性定理。

    (c) 用 Karp 算法计算最大环均值。

!!! example "习题 59A.8"
    给定有向图（min-plus 意义）的邻接矩阵

    $$D = \begin{pmatrix}+\infty & 2 & 5\\+\infty & +\infty & 1\\3 & +\infty & +\infty\end{pmatrix}.$$

    (a) 用热带矩阵幂计算所有点对之间的最短距离。

    (b) 验证结果与 Floyd-Warshall 算法一致。

**热带线性方程组**

!!! example "习题 59A.9"
    求解热带线性方程组 $A \odot x = b$：

    $$A = \begin{pmatrix}2 & 1 & 3\\0 & 4 & 1\\3 & 2 & 0\end{pmatrix}, \quad b = \begin{pmatrix}6\\8\\7\end{pmatrix}.$$

    先计算候选最大解 $x^*$，然后验证 $A \odot x^* = b$ 是否成立。

!!! example "习题 59A.10"
    设 $A$ 是 $n \times n$ 热带非奇异矩阵。证明热带 Cramer 法则给出的解确实满足方程 $A \odot x = b$。

**热带秩**

!!! example "习题 59A.11"
    计算下列矩阵的热带秩和 Barvinok 秩：

    $$A = \begin{pmatrix}0 & 1 & 2\\1 & 2 & 3\\2 & 3 & 4\end{pmatrix}.$$

!!! example "习题 59A.12"
    证明：$n \times n$ 热带矩阵 $A$ 的热带秩 $\mathrm{rk}_T(A) = n$ 当且仅当 $A$ 是热带非奇异的。

**理论证明**

!!! example "习题 59A.13"
    证明热带矩阵乘法与 min-plus 矩阵乘法通过取负映射 $A \mapsto -A$ 互相转换：

    $$(-A) \odot_{\min} (-B) = -(A \odot_{\max} B).$$

!!! example "习题 59A.14"
    设 $A \in \mathbb{T}^{n \times n}$ 的最大环均值为 $\lambda$。证明 $A \oplus \lambda \odot I^{\oplus}$ 的最大环均值也是 $\lambda$。

!!! example "习题 59A.15"
    证明残差的投影性质：对任意 $A \in \mathbb{T}^{m \times n}$ 和 $x \in \mathbb{T}^n$，$A \odot (A \backslash (A \odot x)) = A \odot x$。

---

**本章要点总结：**

1. 热带半环 $(\mathbb{R} \cup \{-\infty\}, \max, +)$ 不是环（加法无逆元），但具有丰富的代数结构，幂等性 $a \oplus a = a$ 是其本质特征。
2. 热带矩阵乘法 $(A \odot B)_{ik} = \max_j(a_{ij}+b_{jk})$ 自然地编码路径优化问题。
3. 热带行列式等于最优指派问题的最优值，可在多项式时间 $O(n^3)$ 计算。
4. 热带特征值是有向图的最大环均值，在强连通图上唯一确定（热带 Perron-Frobenius 定理）。
5. 周期性定理揭示热带矩阵幂的最终周期行为，周期由临界环长度的最大公约数决定。
6. Karp 算法在 $O(n^3)$ 时间内计算最大环均值。
7. Kleene 星给出全源最长/最短路径，Floyd-Warshall 和 Bellman-Ford 都有自然的热带解释。
8. 热带线性方程组 $A \odot x = b$ 的可解性通过残差理论分析，最大解有显式公式。
9. 热带 Cramer 法则给出非奇异情形的解公式。
10. 热带矩阵有多种不等价的秩概念（Barvinok 秩、Kapranov 秩、热带秩），反映了热带代数的深层组合结构。
11. 离散事件系统（调度、交通）的动态方程自然表述为热带线性递推，热带特征值给出系统的渐近周期。

## 练习题

1. **[基础] 在热带半环中计算 $(2 \oplus 5) \odot (3 \oplus 4)$。**
   ??? success "参考答案"
       $(2 \oplus 5) = \max(2, 5) = 5$。
       $(3 \oplus 4) = \max(3, 4) = 4$。
       $5 \odot 4 = 5 + 4 = 9$。

2. **[矩阵乘法] 设 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$。计算热带矩阵积 $A \odot B$。**
   ??? success "参考答案"
       $(A \odot B)_1 = \max(1+0, 2+1) = \max(1, 3) = 3$。
       $(A \odot B)_2 = \max(3+0, 4+1) = \max(3, 5) = 5$。
       故 $A \odot B = \begin{pmatrix} 3 \\ 5 \end{pmatrix}$。

3. **[行列式] 计算 $2 \times 2$ 热带矩阵 $A = \begin{pmatrix} 3 & 1 \\ 2 & 4 \end{pmatrix}$ 的热带行列式。它是热带奇异的吗？**
   ??? success "参考答案"
       $\operatorname{tdet}(A) = \max(3+4, 1+2) = \max(7, 3) = 7$。
       由于最大值 7 由唯一的置换（恒等置换）达到，因此该矩阵是热带非奇异的。

4. **[特征值] 什么是热带矩阵的特征值？它与关联图的环路有什么关系？**
   ??? success "参考答案"
       热带特征值 $\lambda$ 满足 $A \odot \mathbf{x} = \lambda \odot \mathbf{x}$。在强连通图中，它等于图的所有有向环的**最大平均权重**。

5. **[计算] 求矩阵 $A = \begin{pmatrix} 2 & 5 \\ 3 & 1 \end{pmatrix}$ 的热带特征值。**
   ??? success "参考答案"
       环路包括：自环 (1,1) 权 2；自环 (2,2) 权 1；双向环 (1,2,1) 权 $(5+3)/2 = 4$。
       最大环均值为 4，故热带特征值为 4。

6. **[路径] Kleene 星算子 $A^* = I \oplus A \oplus A^{\odot 2} \oplus \dots$ 在图论中代表什么？**
   ??? success "参考答案"
       代表全源最长路径（max-plus）或全源最短路径（min-plus）。$(A^*)_{ij}$ 是从顶点 $i$ 到 $j$ 的所有路径权重的最大值。

7. **[线性方程] 热带方程 $A \odot \mathbf{x} = \mathbf{b}$ 与经典线性方程组的主要区别是什么？**
   ??? success "参考答案"
       热带方程通常没有精确解。我们通常寻找“最大次解” $\mathbf{x}^* = A \backslash \mathbf{b}$（即满足 $A \odot \mathbf{x} \le \mathbf{b}$ 的最大向量）。只有当 $A \odot \mathbf{x}^* = \mathbf{b}$ 成立时，原方程才有精确解。

8. **[调度] 在火车时刻表模型中，热带特征值代表什么物理量？**
   ??? success "参考答案"
       代表系统的**最小稳定循环周期**。它决定了在不发生连锁延误的前提下，系统最快能以多高的频率运行。

9. **[秩] 为什么热带矩阵有多种不同的“秩”定义？**
   ??? success "参考答案"
       因为在热带半环中，行列式定义、线性无关定义和因子分解定义不再相互等价。热带秩（最大非奇异子阵）通常小于或等于因子秩（Barvinok 秩）。

10. **[爱因斯坦思考题] 爱因斯坦在考虑经典力学的统计极限时，路径积分往往被最速路径（极值路径）主导。热带代数将“求和”变为“取极值”，将“乘法”变为“加法”。这种“对数变换”后的线性化逻辑，如何反映了从波动（叠加）到粒子（轨道优化）的物理演化？**
    ??? success "参考答案"
        热带代数本质上是**经典极限下的量子力学**。在量子力学中，我们对所有路径的相位进行求和；但在宏观极限下，只有作用量（权重之和）最小/最大的那条路径由于相位干涉而存留。热带代数舍弃了精细的干涉项，只保留了最强的路径。这种从“求和”到“取最”的跨越，正是从概率云坍缩到确定性物理轨道的代数过程。它体现了一种极端的效率原则：宇宙在宏观尺度上不再玩弄概率的叠加，而只执行最优化的路径。

## 本章小结

本章将线性代数移植到了非传统的幂等代数结构——热带半环之上，展示了其在组合优化中的巨大威力：

1. **热带算术**：确立了 $(\max, +)$ 运算体系，将复杂的非线性优化问题转化为具有线性形式的代数问题。
2. **组合行列式**：证明了热带行列式的计算等价于最优指派问题，从而将行列式论与图匹配理论完美挂钩。
3. **图论特征值**：揭示了热带特征值作为图最大环均值的本质，提供了描述调度系统瓶颈的精确代数指标。
4. **路径闭包与 Kleene 星**：利用矩阵幂建立了全源路径优化的代数表达，统一了 Floyd-Warshall 等经典算法。
5. **残差与系统分析**：通过残差理论解决了热带线性方程组的求解与近似问题，为离散事件动力系统的稳定性与周期性分析提供了坚实的数学框架。

