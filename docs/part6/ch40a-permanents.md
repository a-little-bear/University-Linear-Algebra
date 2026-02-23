# 第 40A 章 积和式

<div class="context-flow" markdown>

**前置**：行列式(Ch3) · 矩阵不等式(Ch18) · Majorization 与双随机矩阵(Ch31) · 计算复杂性基础

**本章脉络**：积和式定义与行列式对比 $\to$ 基本代数性质 $\to$ 半正定矩阵的积和式 $\to$ 混合判别式 $\to$ Van der Waerden 猜想与 Egorychev-Falikman 定理 $\to$ Hadamard-积和式不等式 $\to$ Bregman-Minc 不等式 $\to$ Marcus-Merris 不等式 $\to$ Schrijver 不等式 $\to$ London 不等式 $\to$ Valiant 的 $\#P$-完全性 $\to$ Ryser 公式 $\to$ Gurvits 容量方法 $\to$ Hafnian $\to$ 完美匹配计数与 BosonSampling

**延伸**：积和式（permanent）是组合数学、量子光学和计算复杂性的交汇点；$\#P$-完全性揭示了行列式（多项式时间）与积和式之间的根本计算鸿沟；Van der Waerden 猜想的证明是凸分析与组合学交叉的里程碑；BosonSampling 将积和式的难度转化为量子优越性的论证

</div>

行列式 $\det(A) = \sum_{\sigma \in S_n} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}$ 是矩阵论中最基本的量。如果去掉求和中的符号因子 $\operatorname{sgn}(\sigma)$，得到的量就是**积和式**（permanent）。这一看似微小的改变却带来了本质性的差异：行列式可以在 $O(n^3)$ 时间内计算（通过 Gauss 消元），而积和式的计算是 $\#P$-完全的——被认为比 $NP$-完全问题还要困难。

积和式最早由 Cauchy 和 Binet 在 19 世纪初引入，其名称"permanent"由 Cauchy (1812) 命名——因为与行列式不同，交换两行不改变其值，它是"永久不变的"。此后，积和式在组合学（完美匹配计数）、概率论（随机矩阵）和量子计算（玻色子采样）中发展出深刻的应用。

本章系统发展积和式的代数性质、不等式理论和计算复杂性，并介绍与之密切相关的 Hafnian 和混合判别式。

---

## 40A.1 积和式的定义

<div class="context-flow" markdown>

**核心问题**：积和式是什么？它与行列式有什么相似和不同之处？

</div>

!!! definition "定义 40A.1 (积和式)"
    设 $A = (a_{ij}) \in M_n(\mathbb{C})$。$A$ 的**积和式**（permanent）定义为
    $$\operatorname{perm}(A) = \sum_{\sigma \in S_n} \prod_{i=1}^n a_{i,\sigma(i)},$$
    其中求和遍历 $n$ 元对称群 $S_n$ 的所有 $n!$ 个置换。

    对 $m \times n$ 矩阵（$m \le n$），积和式定义为
    $$\operatorname{perm}(A) = \sum_{\substack{S \subseteq \{1,\ldots,n\} \\ |S| = m}} \sum_{\sigma: \{1,\ldots,m\} \to S \text{ 双射}} \prod_{i=1}^m a_{i,\sigma(i)}.$$

积和式与行列式的公式极为相似，但一个符号的差异导致了代数性质和计算复杂性的根本不同。下面的对比表总结了这些差异：

| 性质 | 行列式 $\det(A)$ | 积和式 $\operatorname{perm}(A)$ |
|------|-----------------|-------------------------------|
| 定义 | $\sum_\sigma \operatorname{sgn}(\sigma) \prod a_{i,\sigma(i)}$ | $\sum_\sigma \prod a_{i,\sigma(i)}$ |
| 多重线性 | 是（每行/列） | 是（每行/列） |
| 交换两行 | 变号 | 不变 |
| 两行相同 | $= 0$ | $\ne 0$（一般地） |
| 矩阵乘积 | $\det(AB) = \det(A)\det(B)$ | $\operatorname{perm}(AB) \ne \operatorname{perm}(A)\operatorname{perm}(B)$ |
| 相似不变 | $\det(P^{-1}AP) = \det(A)$ | 一般不成立 |
| 行运算化简 | 行倍加不变 | 改变值 |
| 计算复杂度 | $O(n^3)$（Gauss 消元） | $\#P$-完全（最快精确 $O(2^n n)$） |
| 组合意义 | 有向图权重行列式 | 二部图完美匹配计数 |
| 对非负矩阵 | 可正可负 | $\ge 0$ |

!!! example "例 40A.1"
    $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$。

    $$\det(A) = ad - bc, \quad \operatorname{perm}(A) = ad + bc.$$

    对 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$：$\det(A) = -2$，$\operatorname{perm}(A) = 4 + 6 = 10$。

!!! example "例 40A.2"
    全 1 矩阵 $J_n$（所有元素等于 1）：

    $$\operatorname{perm}(J_n) = \sum_{\sigma \in S_n} 1 = n!, \quad \det(J_n) = 0 \text{ (当 $n \ge 2$)}.$$

    这个例子清晰地展示了两者的本质差异：行列式因行相同而为零，积和式则对每个置换贡献 1。

!!! example "例 40A.3"
    $$\operatorname{perm}\begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix} = 1\cdot5\cdot9 + 1\cdot6\cdot8 + 2\cdot4\cdot9 + 2\cdot6\cdot7 + 3\cdot4\cdot8 + 3\cdot5\cdot7 = 45 + 48 + 72 + 84 + 96 + 105 = 450.$$

    作为对比，$\det = 1(45-48) - 2(36-42) + 3(32-35) = -3+12-9 = 0$。

!!! example "例 40A.4"
    单位矩阵 $I_n$ 的积和式：$\operatorname{perm}(I_n) = 1$（只有恒等置换贡献非零项）。

    更一般地，对任何置换矩阵 $P_\sigma$，$\operatorname{perm}(P_\sigma) = 1$。

---

## 40A.2 基本性质

<div class="context-flow" markdown>

**核心问题**：积和式有哪些代数性质？它与行列式共享哪些性质，哪些不同？

</div>

!!! theorem "定理 40A.1 (积和式的基本性质)"
    设 $A = (a_{ij}) \in M_n(\mathbb{C})$。

    (a) **多重线性**：$\operatorname{perm}(A)$ 关于 $A$ 的每一行（及每一列）是线性函数。即若第 $i$ 行 $\mathbf{r}_i = \alpha \mathbf{u} + \beta \mathbf{v}$，则
    $$\operatorname{perm}(A) = \alpha \operatorname{perm}(A_\mathbf{u}) + \beta \operatorname{perm}(A_\mathbf{v}),$$
    其中 $A_\mathbf{u}$、$A_\mathbf{v}$ 分别是将第 $i$ 行替换为 $\mathbf{u}$、$\mathbf{v}$ 后的矩阵。

    (b) **行置换不变性**：交换 $A$ 的任意两行，$\operatorname{perm}(A)$ 不变。更一般地，对任何 $n$ 阶置换矩阵 $P$，$\operatorname{perm}(PA) = \operatorname{perm}(A)$。

    (c) **列置换不变性**：交换 $A$ 的任意两列，$\operatorname{perm}(A)$ 不变。即 $\operatorname{perm}(AQ) = \operatorname{perm}(A)$ 对任何置换矩阵 $Q$。

    (d) **Laplace 展开**：按第 $i$ 行展开，
    $$\operatorname{perm}(A) = \sum_{j=1}^n a_{ij} \operatorname{perm}(A(i|j)),$$
    其中 $A(i|j)$ 是去掉第 $i$ 行第 $j$ 列的 $(n-1) \times (n-1)$ 子矩阵。注意：**没有** $(-1)^{i+j}$ 符号因子。

    类似地，按第 $j$ 列展开：
    $$\operatorname{perm}(A) = \sum_{i=1}^n a_{ij} \operatorname{perm}(A(i|j)).$$

    (e) **转置不变**：$\operatorname{perm}(A^T) = \operatorname{perm}(A)$。

    (f) **非负性**：若 $A \ge 0$（所有元素非负），则 $\operatorname{perm}(A) \ge 0$。若 $A > 0$（所有元素严格正），则 $\operatorname{perm}(A) > 0$。

    (g) **对角缩放**：设 $D_1 = \operatorname{diag}(d_1, \ldots, d_n)$，$D_2 = \operatorname{diag}(e_1, \ldots, e_n)$，则
    $$\operatorname{perm}(D_1 A D_2) = \left(\prod_{i=1}^n d_i\right) \left(\prod_{j=1}^n e_j\right) \operatorname{perm}(A).$$

??? proof "证明"
    (a) 每个求和项 $\prod_{k=1}^n a_{k,\sigma(k)}$ 关于第 $i$ 行是线性的，因为只含一个因子 $a_{i,\sigma(i)}$。将求和项中的 $a_{i,\sigma(i)} = \alpha u_{\sigma(i)} + \beta v_{\sigma(i)}$ 代入并利用有限和的线性性即得。

    (b) 交换第 $i, j$ 行等价于对置换进行预复合：用对换 $\tau = (ij)$ 替换 $\sigma$ 为 $\tau\sigma$。由于 $\sigma \mapsto \tau\sigma$ 是 $S_n$ 到自身的双射，求和不变。注意在行列式的情形中，$\operatorname{sgn}(\tau\sigma) = \operatorname{sgn}(\tau)\operatorname{sgn}(\sigma) = -\operatorname{sgn}(\sigma)$，所以行列式变号。

    (c) 类似地，交换列对应于后复合对换。

    (d) 将 $\prod_{k=1}^n a_{k,\sigma(k)}$ 按 $\sigma(i) = j$ 分组：固定 $j$，满足 $\sigma(i) = j$ 的置换 $\sigma$ 与 $S_{n-1}$ 的置换（作用在 $\{1,\ldots,n\} \setminus \{i\} \to \{1,\ldots,n\} \setminus \{j\}$ 上）一一对应。因此
    $$\sum_{\sigma: \sigma(i)=j} \prod_{k=1}^n a_{k,\sigma(k)} = a_{ij} \operatorname{perm}(A(i|j)).$$
    对 $j$ 求和即得。

    (e) $\operatorname{perm}(A^T) = \sum_\sigma \prod_i a_{\sigma(i),i} = \sum_\sigma \prod_i a_{i,\sigma^{-1}(i)}$。令 $\tau = \sigma^{-1}$，当 $\sigma$ 遍历 $S_n$ 时 $\tau$ 也遍历 $S_n$，故 $= \sum_{\tau} \prod_i a_{i,\tau(i)} = \operatorname{perm}(A)$。

    (f) 若 $A \ge 0$，则每个求和项 $\prod a_{i,\sigma(i)} \ge 0$，故和非负。若 $A > 0$，则恒等置换的贡献 $\prod a_{ii} > 0$。

    (g) $(D_1 A D_2)_{ij} = d_i a_{ij} e_j$，故 $\prod_i (D_1 A D_2)_{i,\sigma(i)} = \prod_i d_i a_{i,\sigma(i)} e_{\sigma(i)} = (\prod d_i)(\prod e_j) \prod a_{i,\sigma(i)}$（利用 $\sigma$ 是 $\{1,\ldots,n\}$ 的置换）。

!!! example "例 40A.5"
    验证 Laplace 展开。对 $A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix}$，按第 1 行展开：

    $$\operatorname{perm}(A) = 1 \cdot \operatorname{perm}\begin{pmatrix} 5 & 6 \\ 8 & 9 \end{pmatrix} + 2 \cdot \operatorname{perm}\begin{pmatrix} 4 & 6 \\ 7 & 9 \end{pmatrix} + 3 \cdot \operatorname{perm}\begin{pmatrix} 4 & 5 \\ 7 & 8 \end{pmatrix}$$

    $$= 1(45+48) + 2(36+42) + 3(32+35) = 93 + 156 + 201 = 450.$$

---

## 40A.3 半正定 Hermite 矩阵的积和式

<div class="context-flow" markdown>

**核心问题**：当矩阵具有半正定结构时，积和式有何特殊行为？

</div>

!!! theorem "定理 40A.2 (半正定矩阵的积和式非负性)"
    设 $A \in M_n(\mathbb{C})$ 为半正定 Hermite 矩阵。则 $\operatorname{perm}(A) \ge 0$。

    若 $A$ 为正定 Hermite 矩阵，则 $\operatorname{perm}(A) > 0$。

??? proof "证明"
    设 $A$ 为半正定 Hermite 矩阵，则 $A = B^*B$ 对某个矩阵 $B$。由 Cauchy-Binet 型公式的积和式版本：

    $$\operatorname{perm}(B^*B) = \sum_{1 \le j_1, \ldots, j_n \le m} \left| \prod_{i=1}^n b_{j_i, i} \right|_{\text{perm}}$$

    更精确地，利用 Schur 幂积（power product）和特征标理论，可以证明
    $$\operatorname{perm}(A) = \operatorname{perm}(B^*B) = \sum_{\mathbf{j}} |\operatorname{perm}(B[\mathbf{j}, :])|^2 \ge 0,$$
    其中求和遍历所有从 $\{1, \ldots, m\}$ 中取 $n$ 个指标（允许重复）的多重指标 $\mathbf{j}$。

    对正定情形，$A$ 的对角元素 $a_{ii} > 0$，恒等置换的贡献 $\prod a_{ii} > 0$，且其余偶数置换的贡献也非负，因此 $\operatorname{perm}(A) > 0$。

    **另证**（直接方法）：对半正定 Hermite 矩阵 $A$，利用 Fischer 不等式的推广可以证明对任何置换 $\sigma$，如果将 $\sigma$ 分解为不相交的循环 $(i_1, \ldots, i_k)$，则相应的"循环积" $a_{i_1 i_2} a_{i_2 i_3} \cdots a_{i_k i_1}$ 的实部非负。从而 $\operatorname{perm}(A)$ 的所有贡献项求和后为非负实数。

!!! example "例 40A.6"
    $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ 是正定矩阵（特征值为 $3, 1$）。

    $\operatorname{perm}(A) = 2 \cdot 2 + 1 \cdot 1 = 5 > 0$。

!!! example "例 40A.7"
    $A = \begin{pmatrix} 1 & i \\ -i & 1 \end{pmatrix}$ 是半正定 Hermite 矩阵（特征值为 $2, 0$）。

    $\operatorname{perm}(A) = 1 \cdot 1 + i \cdot (-i) = 1 + 1 = 2 \ge 0$。

---

## 40A.4 混合判别式

<div class="context-flow" markdown>

**核心问题**：混合判别式如何推广积和式？它在 Van der Waerden 猜想的证明中扮演什么角色？

</div>

混合判别式（mixed discriminant）是积和式在多个矩阵上的自然推广，它在 Van der Waerden 猜想的证明中起到了关键作用。

!!! definition "定义 40A.2 (混合判别式)"
    设 $A_1, A_2, \ldots, A_n \in M_n(\mathbb{C})$ 为 $n$ 个 $n \times n$ 矩阵。**混合判别式**（mixed discriminant）定义为
    $$D(A_1, A_2, \ldots, A_n) = \frac{1}{n!} \sum_{\sigma \in S_n} \det\bigl(A_{\sigma(1)}[\cdot, 1], A_{\sigma(2)}[\cdot, 2], \ldots, A_{\sigma(n)}[\cdot, n]\bigr),$$
    其中 $A_k[\cdot, j]$ 表示 $A_k$ 的第 $j$ 列，右端的行列式是将第 $j$ 列取自矩阵 $A_{\sigma(j)}$ 所构成的矩阵的行列式。

    等价地，混合判别式可以表示为
    $$D(A_1, \ldots, A_n) = \frac{\partial^n}{\partial t_1 \cdots \partial t_n} \det(t_1 A_1 + t_2 A_2 + \cdots + t_n A_n) \bigg|_{t_1 = \cdots = t_n = 0} \cdot \frac{1}{n!}.$$

!!! theorem "定理 40A.3 (混合判别式与积和式的关系)"
    设 $A = (a_{ij}) \in M_n(\mathbb{C})$，令 $A_i = \operatorname{diag}(a_{i1}, a_{i2}, \ldots, a_{in})$（$A$ 的第 $i$ 行构成的对角矩阵）。则
    $$\operatorname{perm}(A) = n! \cdot D(A_1, A_2, \ldots, A_n).$$

??? proof "证明"
    由混合判别式的定义，
    $$D(A_1, \ldots, A_n) = \frac{1}{n!} \sum_{\sigma \in S_n} \det\bigl(\operatorname{diag}(a_{\sigma(1),1}, \ldots, a_{\sigma(1),n})[\cdot,1], \ldots\bigr).$$
    由于每个 $A_i$ 是对角矩阵，将第 $j$ 列从 $A_{\sigma(j)}$ 中取出得到 $a_{\sigma(j),j} \mathbf{e}_j$。因此构成的矩阵是对角矩阵 $\operatorname{diag}(a_{\sigma(1),1}, a_{\sigma(2),2}, \ldots, a_{\sigma(n),n})$，其行列式为 $\prod_{j=1}^n a_{\sigma(j),j}$。

    因此 $D(A_1, \ldots, A_n) = \frac{1}{n!} \sum_\sigma \prod_j a_{\sigma(j),j} = \frac{1}{n!} \operatorname{perm}(A^T) = \frac{1}{n!} \operatorname{perm}(A)$。

!!! theorem "定理 40A.4 (混合判别式的基本性质)"
    设 $A_1, \ldots, A_n, B \in M_n(\mathbb{C})$ 为 Hermite 矩阵。

    (a) **多重线性**：$D$ 关于每个参数是线性函数。

    (b) **对称性**：$D$ 关于参数的任何置换不变。

    (c) **非负性**：若 $A_1, \ldots, A_n$ 均为半正定 Hermite 矩阵，则 $D(A_1, \ldots, A_n) \ge 0$。

    (d) **Alexandrov 不等式**：若 $A_1, \ldots, A_n$ 为半正定 Hermite 矩阵，则
    $$D(A_1, A_2, A_3, \ldots, A_n)^2 \ge D(A_1, A_1, A_3, \ldots, A_n) \cdot D(A_2, A_2, A_3, \ldots, A_n).$$

??? proof "证明"
    (a) 由行列式关于每一列的线性性和求和的线性性立即得到。

    (b) 对参数的置换 $\tau$ 等价于在求和中用 $\sigma \circ \tau$ 替换 $\sigma$，由 $S_n$ 的群结构不改变求和。

    (c) 这是 Alexandrov (1938) 的结果。证明的核心思想：设 $A_i$ 半正定，则 $\det(t_1 A_1 + \cdots + t_n A_n)$ 在 $t_i \ge 0$ 时是关于 $(t_1, \ldots, t_n)$ 的非负多项式。混合判别式作为该多项式的（归一化）混合偏导数，其非负性可以通过超曲面理论中的 Alexandrov-Fenchel 型论证得到。

    (d) 这是 Alexandrov 不等式的核心形式，类似于 Brunn-Minkowski 不等式的矩阵版本。证明使用了关于半正定矩阵线性组合的行列式的对数凹性。

---

## 40A.5 Van der Waerden 猜想

<div class="context-flow" markdown>

**核心问题**：双随机矩阵的积和式的最小值是什么？何时达到？

</div>

!!! definition "定义 40A.3 (双随机矩阵与 Birkhoff 多面体)"
    非负矩阵 $A \ge 0$ 称为**双随机矩阵**（doubly stochastic matrix），若每行每列之和均为 1：
    $$\sum_{j=1}^n a_{ij} = 1 \quad \forall\, i, \qquad \sum_{i=1}^n a_{ij} = 1 \quad \forall\, j.$$

    全体 $n$ 阶双随机矩阵构成凸集 $\Omega_n$，称为 **Birkhoff 多面体**。由 Birkhoff-von Neumann 定理（定理 31.X），$\Omega_n$ 的顶点恰好是 $n!$ 个置换矩阵。

!!! theorem "定理 40A.5 (Van der Waerden 猜想 / Egorychev-Falikman 定理, 1981)"
    对所有 $n \times n$ 双随机矩阵 $A \in \Omega_n$，
    $$\operatorname{perm}(A) \ge \frac{n!}{n^n}.$$
    等号成立当且仅当 $A = J_n/n$（即 $A$ 为所有元素等于 $1/n$ 的矩阵）。

这是组合数学中最著名的猜想之一，由 Van der Waerden 于 1926 年提出，经过半个多世纪的努力，最终在 1981 年由 G.P. Egorychev 和 D.I. Falikman 独立证明。

??? proof "证明"
    **Egorychev 的证明**基于混合判别式和 Alexandrov 不等式，分为以下几个步骤。

    **步骤一（最小值的存在性）**：$\operatorname{perm}$ 是 $\Omega_n$ 上的连续函数，$\Omega_n$ 是 $\mathbb{R}^{n^2}$ 中的紧凸集。由 Weierstrass 极值定理，$\operatorname{perm}$ 在 $\Omega_n$ 上取到最小值。设最小值在 $A^* \in \Omega_n$ 处取到。

    **步骤二（最小值点的内部性）**：首先证明 $A^*$ 的所有元素严格为正（$A^* > 0$），即 $A^*$ 位于 $\Omega_n$ 的内部。

    假设 $A^*$ 有某个元素 $a^*_{ij} = 0$。利用 Birkhoff 多面体的结构和积和式的多重线性性，可以构造 $\Omega_n$ 中的路径使积和式严格减小，这与 $A^*$ 是最小值点矛盾。

    **步骤三（余积和式的均匀性）**：在 $A^* > 0$ 的条件下，利用 Lagrange 乘子法。$\operatorname{perm}$ 关于 $a_{ij}$ 的偏导数为**余积和式**（permanent cofactor）$\operatorname{perm}(A^*(i|j))$。双随机矩阵的约束条件为行和与列和均等于 1。

    在最小值点，KKT 条件给出：存在常数 $\lambda_1, \ldots, \lambda_n$（行约束乘子）和 $\mu_1, \ldots, \mu_n$（列约束乘子），使得
    $$\operatorname{perm}(A^*(i|j)) = \lambda_i + \mu_j \quad \forall\, i, j.$$
    利用 $A^*$ 的双随机性和 Laplace 展开 $\operatorname{perm}(A^*) = \sum_j a^*_{ij} \operatorname{perm}(A^*(i|j))$，可以推导出
    $$\operatorname{perm}(A^*(i|j)) = c \quad \text{对所有 } i, j$$
    为某个常数 $c$。

    **步骤四（利用 Alexandrov 不等式）**：这是证明的核心。将积和式表示为混合判别式（定理 40A.3），对 $A^*$ 的行向量构成的对角矩阵 $H_i = \operatorname{diag}(a^*_{i1}, \ldots, a^*_{in})$，
    $$\operatorname{perm}(A^*) = n! \cdot D(H_1, \ldots, H_n).$$

    Alexandrov 不等式（定理 40A.4(d)）给出
    $$D(H_i, H_j, H_3, \ldots, H_n)^2 \ge D(H_i, H_i, H_3, \ldots, H_n) \cdot D(H_j, H_j, H_3, \ldots, H_n).$$

    由步骤三中余积和式的均匀性条件，结合混合判别式的对称性和多重线性性，通过精细的归纳论证，可以证明在最小值点处
    $$H_1 = H_2 = \cdots = H_n = \frac{1}{n} I_n,$$
    即 $A^*$ 的每一行都是 $(1/n, \ldots, 1/n)$，因此 $A^* = J_n/n$。

    **步骤五（计算最小值）**：
    $$\operatorname{perm}(J_n/n) = \frac{1}{n^n} \operatorname{perm}(J_n) = \frac{n!}{n^n}.$$

!!! example "例 40A.8"
    $n = 2$：$\Omega_2$ 中的双随机矩阵为 $A_t = \begin{pmatrix} t & 1-t \\ 1-t & t \end{pmatrix}$（$0 \le t \le 1$）。

    $$\operatorname{perm}(A_t) = t^2 + (1-t)^2 = 2t^2 - 2t + 1.$$

    最小值在 $t = 1/2$ 时取到，$\operatorname{perm}(A_{1/2}) = 1/2 = 2!/2^2$。

!!! example "例 40A.9"
    $n = 3$：$\operatorname{perm}(J_3/3) = 3!/3^3 = 6/27 = 2/9 \approx 0.222$。

    对比：置换矩阵 $P$ 的积和式为 1，远大于 $2/9$。均匀矩阵 $J_n/n$ 确实给出了最小值。

!!! note "注记 40A.1 (Van der Waerden 猜想的历史)"
    Van der Waerden 猜想提出于 1926 年，是组合数学中最著名的猜想之一。在被证明之前，许多数学家建立了部分结果：

    - Marcus-Newman (1959)：证明了 $n = 3$ 的情形。
    - London (1971)：证明了下界 $n!/((n-1)^n + n - 1)$。
    - Friedland (1979)：证明了 $\operatorname{perm}(A) \ge e^{-n}$ 对双随机矩阵成立。

    1981 年，Egorychev 和 Falikman 独立地给出了完整证明。Egorychev 的证明使用了混合判别式和 Alexandrov 不等式；Falikman 的证明使用了类似但略有不同的凸分析工具。两人因此获得了 1982 年的 Fulkerson 奖。

---

## 40A.6 积和式的不等式

<div class="context-flow" markdown>

**核心问题**：积和式满足哪些重要的上界和下界不等式？

</div>

### Hadamard-积和式不等式

!!! theorem "定理 40A.6 (Hadamard-积和式不等式)"
    设 $A \ge 0$ 为 $n \times n$ 非负矩阵，行和为 $r_i = \sum_{j=1}^n a_{ij}$。则
    $$\operatorname{perm}(A) \le \prod_{i=1}^n r_i.$$
    等号成立当且仅当 $A$ 的每一行至多有一个非零元素，或 $A$ 的所有行成比例。

??? proof "证明"
    将 $A$ 的每一行归一化：设 $r_i > 0$（若 $r_i = 0$ 则该行全零，$\operatorname{perm}(A) = 0$，不等式平凡成立）。令 $\tilde{a}_{ij} = a_{ij}/r_i$，则 $\tilde{A}$ 是行随机矩阵（每行之和为 1）。

    由积和式的多重线性性（定理 40A.1(a)和(g)），
    $$\operatorname{perm}(A) = \prod_{i=1}^n r_i \cdot \operatorname{perm}(\tilde{A}).$$

    因此只需证明对行随机非负矩阵 $\tilde{A}$，$\operatorname{perm}(\tilde{A}) \le 1$。

    这由 AM-GM 不等式得到。对每个置换 $\sigma$，
    $$\prod_{i=1}^n \tilde{a}_{i,\sigma(i)} \le \frac{1}{n!} \sum_{\sigma \in S_n} \prod_{i=1}^n \tilde{a}_{i,\sigma(i)}$$
    这不直接有用。换一种方法：对每一行，$\sum_j \tilde{a}_{ij} = 1$。注意
    $$\operatorname{perm}(\tilde{A}) = \sum_\sigma \prod_i \tilde{a}_{i,\sigma(i)} \le \prod_i \left(\sum_j \tilde{a}_{ij}\right) = 1,$$
    其中最后一步的不等式可以通过归纳法证明：对 $n = 1$ 显然成立；对 $n \ge 2$，利用 Laplace 展开和归纳假设。

### Bregman-Minc 不等式

!!! theorem "定理 40A.7 (Bregman-Minc 不等式, 1973)"
    设 $A$ 为 $n \times n$ $(0,1)$-矩阵（所有元素为 0 或 1），行和为 $r_i$。则
    $$\operatorname{perm}(A) \le \prod_{i=1}^n (r_i!)^{1/r_i}.$$
    等号成立当且仅当 $A$（经行列重排后）是若干全 1 方阵的直和。

??? proof "证明"
    Minc 于 1963 年提出此猜想，Bregman 于 1973 年利用双随机矩阵与熵的联系给出了证明。我们概述 Radhakrishnan (1997) 的信息论证明，这是目前最简洁的证明。

    将 $A$ 视为二部图 $G$ 的邻接矩阵。$\operatorname{perm}(A)$ 等于 $G$ 的完美匹配数 $M$。设 $\sigma$ 为均匀随机选取的完美匹配。

    定义随机变量 $X_j = \sigma(j)$（第 $j$ 个左顶点的匹配对象）。则

    $$\log M = H(X_1, X_2, \ldots, X_n) = \sum_{j=1}^n H(X_j \mid X_1, \ldots, X_{j-1}),$$

    其中 $H$ 表示 Shannon 熵。关键观察：给定 $X_1, \ldots, X_{j-1}$，$X_j$ 的取值范围至多是第 $j$ 行中值为 1 的列中，去除已被使用的列。

    通过仔细估计条件熵，利用 $\log$ 的凹性和对称性论证：选取行的随机排列 $\pi$ 并取期望，可以证明
    $$H(X_{\pi(j)} \mid X_{\pi(1)}, \ldots, X_{\pi(j-1)}) \le \frac{1}{r_j} \log(r_j!),$$
    对 $j$ 求和即得 $\log M \le \sum_j \frac{1}{r_j} \log(r_j!)$，即 $M \le \prod (r_j!)^{1/r_j}$。

!!! example "例 40A.10"
    对 $J_3$（$3 \times 3$ 全 1 矩阵），$r_i = 3$。Bregman-Minc 上界为 $(3!)^{3/3} = 6$。实际 $\operatorname{perm}(J_3) = 6$，等号成立。

!!! example "例 40A.11"
    对 $A = \begin{pmatrix} 1 & 1 & 0 \\ 1 & 0 & 1 \\ 0 & 1 & 1 \end{pmatrix}$，$r_i = 2$。Bregman-Minc 上界为 $(2!)^{3/2} = 2^{3/2} = 2\sqrt{2} \approx 2.83$。实际 $\operatorname{perm}(A) = 2 \le 2.83$。

### Marcus-Merris 不等式

!!! theorem "定理 40A.8 (Marcus-Merris 不等式)"
    设 $A \in M_n(\mathbb{C})$ 为半正定 Hermite 矩阵。则
    $$\operatorname{perm}(A) \ge \det(A).$$
    等号 $\operatorname{perm}(A) = \det(A)$ 成立当且仅当 $A$ 为对角矩阵或 $A$ 有零特征值。

??? proof "证明"
    $$\operatorname{perm}(A) - \det(A) = \sum_{\sigma \in S_n} (1 - \operatorname{sgn}(\sigma)) \prod_{i=1}^n a_{i,\sigma(i)}.$$

    对偶数置换（$\operatorname{sgn}(\sigma) = 1$），贡献为 0；对奇数置换（$\operatorname{sgn}(\sigma) = -1$），贡献为 $2\prod_i a_{i,\sigma(i)}$。因此

    $$\operatorname{perm}(A) - \det(A) = 2 \sum_{\sigma: \operatorname{sgn}(\sigma) = -1} \prod_{i=1}^n a_{i,\sigma(i)}.$$

    对半正定 Hermite 矩阵 $A$，我们需要证明每个奇数置换 $\sigma$ 对应的积 $\prod_i a_{i,\sigma(i)}$ 非负（当 $A$ 为实矩阵时）。

    对实半正定矩阵，$a_{ij} = a_{ji}$，且主子式非负。将置换 $\sigma$ 分解为不相交的循环 $(i_1, i_2, \ldots, i_k)$。对应的积为
    $$a_{i_1, i_2} a_{i_2, i_3} \cdots a_{i_k, i_1}.$$

    对长度为 2 的循环（对换 $(i,j)$），积为 $a_{ij}^2 \ge 0$。对更长的循环，利用半正定矩阵的元素满足 $|a_{ij}| \le \sqrt{a_{ii} a_{jj}}$（Cauchy-Schwarz），通过 AM-GM 不等式可以证明循环积的乘积非负。

    综合所有循环的贡献，$\operatorname{perm}(A) - \det(A) \ge 0$。

    等号条件：若 $A$ 为对角矩阵，则只有恒等置换贡献非零项，$\operatorname{perm}(A) = \det(A) = \prod a_{ii}$。

!!! example "例 40A.12"
    $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$（正定）。$\det(A) = 3$，$\operatorname{perm}(A) = 5 \ge 3$。

    $A = \begin{pmatrix} 3 & 0 \\ 0 & 5 \end{pmatrix}$（对角）。$\det(A) = \operatorname{perm}(A) = 15$，等号成立。

### Schrijver 不等式

!!! theorem "定理 40A.9 (Schrijver 不等式, 1998)"
    设 $A$ 为 $n \times n$ $(0,1)$-矩阵，每行每列之和均为 $d$（即 $A$ 为 $d$-正则二部图的邻接矩阵）。则
    $$\operatorname{perm}(A) \ge \prod_{i=1}^n \frac{(d - n + i)^+}{i} = \frac{((d-n+1)^+)!}{(d!)^0} \cdots$$

    更精确地，Schrijver 的下界为
    $$\operatorname{perm}(A) \ge \left(\frac{d-1}{d}\right)^{(d-1)n} \cdot d^n.$$

    即 $d$-正则 $(0,1)$-矩阵的积和式满足
    $$\operatorname{perm}(A) \ge \left(\frac{(d-1)^{d-1}}{d^{d-2}}\right)^n.$$

??? proof "证明"
    Schrijver 的证明使用了 Van der Waerden 猜想（定理 40A.5）作为核心工具。

    设 $A$ 为 $d$-正则 $(0,1)$-矩阵，则 $A/d$ 为双随机矩阵。由 Van der Waerden 猜想的推广（对具有特定支撑结构的双随机矩阵），Schrijver 利用了以下关键引理：

    **引理**：设 $B$ 为双随机矩阵且 $B$ 的支撑（非零元素位置）包含在 $d$-正则 $(0,1)$-矩阵 $A$ 的支撑中。则 $\operatorname{perm}(B) \ge \operatorname{perm}(A/d)$，在 $B = A/d$ 时取到。

    对 $B = A/d$，$\operatorname{perm}(A/d) = d^{-n} \operatorname{perm}(A)$。将此与 Egorychev-Falikman 定理的精细化版本结合，通过对 $d$-正则二部图的结构分析（利用 Hall 定理和交替路径），得到
    $$\operatorname{perm}(A/d) \ge \left(\frac{d-1}{d}\right)^{(d-1)n} = \frac{(d-1)^{(d-1)n}}{d^{(d-1)n}},$$
    从而 $\operatorname{perm}(A) \ge d^n \cdot \frac{(d-1)^{(d-1)n}}{d^{(d-1)n}} = \frac{(d-1)^{(d-1)n}}{d^{(d-2)n}}$。

!!! note "注记 40A.2 (Schrijver 不等式的意义)"
    Schrijver 不等式给出了 $(0,1)$-矩阵积和式的**下界**，与 Bregman-Minc 的上界互补。对 $d$-正则情形，两个界分别为：

    - 下界（Schrijver）：$\operatorname{perm}(A) \ge \left(\frac{(d-1)^{d-1}}{d^{d-2}}\right)^n$
    - 上界（Bregman-Minc）：$\operatorname{perm}(A) \le (d!)^{n/d}$

    当 $d$ 固定而 $n \to \infty$ 时，这两个界的比值是有界的。

### London 不等式与 Birkhoff 多面体上的单调性

!!! theorem "定理 40A.10 (London 不等式, 1971)"
    设 $A, B \in \Omega_n$（双随机矩阵），且 $A = \alpha B + (1-\alpha) J_n/n$，其中 $0 \le \alpha \le 1$。则
    $$\operatorname{perm}(A) \ge \alpha^n \operatorname{perm}(B) + (1-\alpha^n) \frac{n!}{n^n}.$$

    特别地，沿 Birkhoff 多面体中从 $J_n/n$ 到任何双随机矩阵 $B$ 的线段，积和式是单调递增的。

??? proof "证明"
    利用积和式的多重线性性。$A$ 的第 $i$ 行为 $\alpha \mathbf{b}_i + (1-\alpha) \frac{1}{n}\mathbf{1}^T$。将 $\operatorname{perm}(A)$ 按多重线性性展开：

    $$\operatorname{perm}(A) = \sum_{S \subseteq \{1,\ldots,n\}} \alpha^{|S|}(1-\alpha)^{n-|S|} \operatorname{perm}(M_S),$$

    其中 $M_S$ 是取第 $i$ 行为 $\mathbf{b}_i$（当 $i \in S$）或 $\frac{1}{n}\mathbf{1}^T$（当 $i \notin S$）的矩阵。

    当 $S = \{1,\ldots,n\}$ 时，$M_S = B$，贡献 $\alpha^n \operatorname{perm}(B)$。

    对 $S \ne \{1,\ldots,n\}$，利用 Van der Waerden 猜想的推论，每个 $\operatorname{perm}(M_S)$ 都有下界。通过仔细估计每一项，可以得到所需的不等式。

    单调性的结论：设 $A_t = (1-t) J_n/n + t B$（$0 \le t \le 1$）。则 $\frac{d}{dt} \operatorname{perm}(A_t) \ge 0$，因为 $\operatorname{perm}(A_t)$ 是 $t$ 的多项式，在 $t = 0$ 处取最小值 $n!/n^n$（由 Van der Waerden 猜想），且 $\operatorname{perm}(A_1) = \operatorname{perm}(B) \ge n!/n^n$。

!!! example "例 40A.13"
    $n = 2$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$，$A_t = \begin{pmatrix} (1+t)/2 & (1-t)/2 \\ (1-t)/2 & (1+t)/2 \end{pmatrix}$。

    $\operatorname{perm}(A_t) = ((1+t)/2)^2 + ((1-t)/2)^2 = (1+t^2)/2$。

    这在 $t = 0$ 时取最小值 $1/2 = 2!/2^2$，在 $t = 1$ 时取值 $1 = \operatorname{perm}(B)$，确实单调递增。

---

## 40A.7 计算复杂性

<div class="context-flow" markdown>

**核心问题**：为什么积和式比行列式难算得多？最快的精确算法是什么？

</div>

### Valiant 的 $\#P$-完全性

!!! theorem "定理 40A.11 (Valiant, 1979: $\#P$-完全性)"
    计算 $(0,1)$-矩阵的积和式是 $\#P$-完全问题。即使限制为 $(0,1)$-矩阵（每个元素为 0 或 1），计算积和式仍然是 $\#P$-完全的。

??? proof "证明（概要）"
    **$\#P$ 的定义**：$\#P$ 是包含所有如下形式问题的复杂类——给定 $NP$ 关系 $R$（即对 $(x, y)$ 的验证可在多项式时间内完成），问满足 $R(x, y)$ 的见证 $y$ 有多少个。

    **积和式属于 $\#P$**：$\operatorname{perm}(A) = \sum_\sigma \prod a_{i,\sigma(i)}$。对 $(0,1)$-矩阵，这等于满足 $\prod a_{i,\sigma(i)} = 1$ 的置换 $\sigma$ 的数量，即二部图的完美匹配数。验证一个给定的 $\sigma$ 是否为完美匹配可以在多项式时间内完成，故积和式属于 $\#P$。

    **$\#P$-硬度**：Valiant 的证明通过归约证明计算积和式至少与计数任何 $NP$ 问题的解一样困难。具体地，他证明了计算 $\#3\text{-SAT}$（3-SAT 公式的满足赋值数）可以多项式时间归约到计算某个 $(0,1)$-矩阵的积和式。

    归约的关键步骤：

    1. 将 $\#3\text{-SAT}$ 实例编码为一个带权二部图。
    2. 将带权图中可能出现的权重 $-1$ 通过"gadget"构造消除，转化为 $(0,1)$-矩阵的积和式。
    3. 这一步需要精巧的构造：利用 $a - b = a + 3b - 4b$ 的恒等式和模 4 运算来处理负权重。

    因此，除非 $P = \#P$，否则不存在计算积和式的多项式时间算法。

!!! note "注记 40A.3 ($\#P$ 与 $NP$ 的关系)"
    $\#P$ 是**计数**复杂类——不仅问"是否存在"某个对象，而是问"有多少个"。$\#P$-完全问题至少与 $NP$-完全问题一样难，且在大多数复杂性理论假设下严格更难。

    Toda 定理 (1991) 指出 $PH \subseteq P^{\#P}$，即多项式层次可以用一次 $\#P$ 预言机访问来模拟。这意味着如果积和式可以高效计算，则整个多项式层次都会坍缩。

    行列式与积和式之间的计算鸿沟是理论计算机科学中最深刻的现象之一：

    - $\det$：多项式时间（$O(n^3)$ 或 $O(n^\omega)$，$\omega \approx 2.37$）；
    - $\operatorname{perm}$：$\#P$-完全（除非 $P = \#P$，没有多项式时间算法）。

    这种差异的根源在于行列式中的符号交替允许大量消除（cancellation），使得 Gauss 消元成为可能；而积和式中所有项同号，没有消除发生。

### Ryser 公式

!!! theorem "定理 40A.12 (Ryser 公式, 1963)"
    对 $A = (a_{ij}) \in M_n(\mathbb{C})$，
    $$\operatorname{perm}(A) = (-1)^n \sum_{S \subseteq \{1,\ldots,n\}} (-1)^{|S|} \prod_{i=1}^n \sum_{j \in S} a_{ij}.$$
    此公式的计算复杂度为 $O(2^n n)$，远优于直接定义的 $O(n! \cdot n)$。

??? proof "证明"
    证明基于**容斥原理**（inclusion-exclusion principle）。

    **第一步：建立函数与双射的关系。** 对 $S \subseteq \{1, \ldots, n\}$，定义
    $$P(S) = \prod_{i=1}^n \left(\sum_{j \in S} a_{ij}\right).$$
    展开乘积，$P(S) = \sum_{f: \{1,\ldots,n\} \to S} \prod_{i=1}^n a_{i,f(i)}$，其中求和遍历所有从 $\{1,\ldots,n\}$ 到 $S$ 的函数 $f$（允许非单射）。

    **第二步：提取双射的贡献。** 积和式 $\operatorname{perm}(A) = \sum_{\sigma \in S_n} \prod_i a_{i,\sigma(i)}$，其中求和仅遍历**双射**（置换）。

    对函数 $f: \{1,\ldots,n\} \to \{1,\ldots,n\}$，设其像集为 $\operatorname{Im}(f)$。$f$ 是双射当且仅当 $|\operatorname{Im}(f)| = n$。

    由容斥原理，"$f$ 的像集恰好为 $\{1,\ldots,n\}$" 可以表示为：

    $$[\operatorname{Im}(f) = \{1,\ldots,n\}] = \sum_{T \subseteq \{1,\ldots,n\}} (-1)^{n - |T|} [\operatorname{Im}(f) \subseteq T].$$

    **第三步：求和。** 对所有函数 $f$ 的贡献求和：

    $$\operatorname{perm}(A) = \sum_f [\operatorname{Im}(f) = \{1,\ldots,n\}] \prod_i a_{i,f(i)}$$

    $$= \sum_{T \subseteq \{1,\ldots,n\}} (-1)^{n-|T|} \sum_{f: \{1,\ldots,n\} \to T} \prod_i a_{i,f(i)}$$

    $$= \sum_{T \subseteq \{1,\ldots,n\}} (-1)^{n-|T|} P(T)$$

    $$= (-1)^n \sum_{T \subseteq \{1,\ldots,n\}} (-1)^{|T|} P(T)$$

    $$= (-1)^n \sum_{S \subseteq \{1,\ldots,n\}} (-1)^{|S|} \prod_{i=1}^n \sum_{j \in S} a_{ij}.$$

    **第四步：复杂度分析。** $\{1,\ldots,n\}$ 有 $2^n$ 个子集。对每个子集 $S$，计算 $P(S) = \prod_i (\sum_{j \in S} a_{ij})$ 需要 $O(n)$ 时间（若利用 Gray 码遍历子集，每次仅改变一个元素，则每个行和只需 $O(1)$ 时间更新）。因此总复杂度为 $O(2^n n)$。

!!! example "例 40A.14"
    用 Ryser 公式计算 $\operatorname{perm}\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$。

    $n = 2$。$S$ 遍历 $\{1,2\}$ 的子集：

    - $S = \emptyset$：$\prod r_i(\emptyset) = 0 \cdot 0 = 0$。
    - $S = \{1\}$：$r_1 = a_{11} = 1, r_2 = a_{21} = 3$，积 $= 3$。
    - $S = \{2\}$：$r_1 = a_{12} = 2, r_2 = a_{22} = 4$，积 $= 8$。
    - $S = \{1,2\}$：$r_1 = 3, r_2 = 7$，积 $= 21$。

    $\operatorname{perm} = (-1)^2[(-1)^0 \cdot 0 + (-1)^1(3 + 8) + (-1)^2 \cdot 21] = 0 - 11 + 21 = 10$。

    验证：$1 \cdot 4 + 2 \cdot 3 = 10$。正确。

!!! example "例 40A.15"
    用 Ryser 公式计算 $\operatorname{perm}\begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 1 \\ 1 & 1 & 1 \end{pmatrix}$。

    $n = 3$，对每个子集 $S$，$r_i(S) = |S|$（因为所有元素为 1），故 $P(S) = |S|^3$。

    $$\operatorname{perm} = (-1)^3 \sum_{k=0}^{3} (-1)^k \binom{3}{k} k^3 = -[0 - 3 \cdot 1 + 3 \cdot 8 - 27] = -[0 - 3 + 24 - 27] = -(-6) = 6.$$

    正确：$\operatorname{perm}(J_3) = 3! = 6$。

### Gurvits 容量方法

!!! note "注记 40A.4 (Gurvits 的容量方法)"
    Gurvits (2006) 提出了一种基于**容量**（capacity）的积和式下界方法，给出了 Van der Waerden 猜想的另一种优雅证明。

!!! definition "定义 40A.4 (多项式的容量)"
    对具有非负系数的齐次多项式 $p(x_1, \ldots, x_n) = \sum_\alpha c_\alpha x^\alpha$（$\alpha$ 为多重指标），定义其**容量**（capacity）为
    $$\operatorname{cap}(p) = \inf_{x_1, \ldots, x_n > 0} \frac{p(x_1, \ldots, x_n)}{\prod_{i=1}^n x_i^{d_i}},$$
    其中 $d_i$ 为 $p$ 关于 $x_i$ 的次数（对齐次多项式，$\sum d_i = \deg p$）。

    对 $n$ 次齐次 $n$ 元多项式，$d_i$ 可以视为 $p$ 的"期望次数"分布。

!!! theorem "定理 40A.13 (Gurvits 容量下界)"
    设 $A \ge 0$ 为 $n \times n$ 非负矩阵，定义生成多项式
    $$p_A(x_1, \ldots, x_n) = \prod_{i=1}^n \left(\sum_{j=1}^n a_{ij} x_j\right).$$
    则 $\operatorname{perm}(A)$ 是 $p_A$ 中 $x_1 x_2 \cdots x_n$ 项的系数乘以 $n!$，且
    $$\operatorname{perm}(A) \ge \operatorname{cap}(p_A).$$

    对双随机矩阵 $A$，$\operatorname{cap}(p_A) \ge n!/n^n$，从而给出 Van der Waerden 猜想的另一种证明。

    Gurvits 的方法还推广到了**稳定多项式**（stable polynomials）理论，在组合学和优化中有广泛应用。

!!! note "注记 40A.5 (近似算法)"
    虽然精确计算积和式是 $\#P$-完全的，Jerrum, Sinclair 和 Vigoda (2004) 给出了一个**完全多项式随机近似方案**（FPRAS）：对非负矩阵 $A \ge 0$，可以在多项式时间内以高概率输出 $\operatorname{perm}(A)$ 的 $(1 \pm \epsilon)$ 倍近似值。该算法基于马尔可夫链蒙特卡罗（MCMC）方法。

    Gurvits 的容量方法还给出了一个**确定性**的 $e^n$ 近似因子：$\operatorname{cap}(p_A) \le \operatorname{perm}(A) \le e^n \operatorname{cap}(p_A)$。

---

## 40A.8 Hafnian

<div class="context-flow" markdown>

**核心问题**：如何将积和式的完美匹配计数推广到非二部图？

</div>

!!! definition "定义 40A.5 (Hafnian)"
    对 $2n \times 2n$ 对称矩阵 $A = (a_{ij})$（$a_{ii} = 0$），**Hafnian** 定义为
    $$\operatorname{hf}(A) = \sum_{M \in \mathcal{M}_{2n}} \prod_{\{i,j\} \in M} a_{ij},$$
    其中 $\mathcal{M}_{2n}$ 是 $\{1, \ldots, 2n\}$ 的所有**完美匹配**（即将 $2n$ 个元素配成 $n$ 对的方式）的集合，$|\mathcal{M}_{2n}| = (2n-1)!! = (2n)!/(2^n n!)$。

    等价地，
    $$\operatorname{hf}(A) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \prod_{i=1}^n a_{\sigma(2i-1), \sigma(2i)}.$$

!!! theorem "定理 40A.14 (Hafnian 的基本性质)"
    (a) 若 $G$ 为简单图，$A$ 为其邻接矩阵（$a_{ij} = 1$ 当 $\{i,j\} \in E(G)$，$a_{ii} = 0$），则 $\operatorname{hf}(A)$ 等于 $G$ 的**完美匹配数**。

    (b) 对二部图 $G = (U \cup V, E)$（$|U| = |V| = n$），若将 $G$ 表示为 $2n$ 个顶点的图，其邻接矩阵为 $\begin{pmatrix} 0 & B \\ B^T & 0 \end{pmatrix}$，则
    $$\operatorname{hf}\begin{pmatrix} 0 & B \\ B^T & 0 \end{pmatrix} = \operatorname{perm}(B).$$

    (c) 计算 Hafnian 也是 $\#P$-完全问题。最快的精确算法复杂度为 $O(2^n n^2)$（利用类似 Ryser 公式的容斥方法）。

    (d) **Loop Hafnian**：对不要求 $a_{ii} = 0$ 的对称矩阵，定义带环 Hafnian（loop Hafnian），将完美匹配推广为允许"自环"的匹配。

??? proof "证明"
    (a) 与积和式计数二部图完美匹配的证明完全类似。$\operatorname{hf}(A)$ 中的每一项 $\prod_{\{i,j\} \in M} a_{ij}$ 在 $M$ 是 $G$ 的完美匹配时等于 1（所有边都存在），否则为 0。

    (b) 设 $U = \{u_1, \ldots, u_n\}$，$V = \{v_1, \ldots, v_n\}$。$2n$ 个顶点的完美匹配将每个 $u_i$ 匹配到某个 $v_j$（由二部结构，$U$ 内部和 $V$ 内部无边，所以匹配只能跨越 $U$ 和 $V$）。这样的匹配一一对应于置换 $\sigma \in S_n$（$u_i$ 匹配 $v_{\sigma(i)}$），从而 $\operatorname{hf} = \sum_\sigma \prod b_{i,\sigma(i)} = \operatorname{perm}(B)$。

!!! example "例 40A.16"
    完全图 $K_4$ 的邻接矩阵 $A = J_4 - I_4$。$K_4$ 的完美匹配数：

    将 4 个顶点配成 2 对，有 $3!! = 3$ 种方式：$\{1,2\}\{3,4\}$，$\{1,3\}\{2,4\}$，$\{1,4\}\{2,3\}$。

    $\operatorname{hf}(A) = 3$。

!!! example "例 40A.17"
    六边形（6 个顶点的环图 $C_6$）的完美匹配数。

    $C_6$ 的邻接矩阵 $A$ 中，$a_{ij} = 1$ 当 $|i-j| = 1 \pmod{6}$。

    $C_6$ 的完美匹配：将 6 个相邻顶点配对。通过枚举可知有 $\operatorname{hf}(A) = 3$ 个完美匹配。

---

## 40A.9 完美匹配计数与 BosonSampling

<div class="context-flow" markdown>

**核心问题**：积和式在图论和量子计算中有何应用？

</div>

### 完美匹配计数

!!! theorem "定理 40A.15 (积和式 = 二部图完美匹配数)"
    设 $G = (U \cup V, E)$ 为二部图，$|U| = |V| = n$，邻接矩阵 $B = (b_{ij})$（$b_{ij} = 1$ 当 $(u_i, v_j) \in E$，否则为 0）。则 $G$ 的**完美匹配**（perfect matching）的数量等于 $\operatorname{perm}(B)$。

??? proof "证明"
    完美匹配是边集 $M \subseteq E$，使得 $U$ 和 $V$ 中每个顶点恰被 $M$ 中一条边覆盖。这等价于选择一个双射 $\sigma: \{1,\ldots,n\} \to \{1,\ldots,n\}$（$u_i$ 匹配 $v_{\sigma(i)}$），使得 $b_{i,\sigma(i)} = 1$ 对所有 $i$ 成立。

    因此完美匹配数为
    $$|\{\sigma \in S_n : b_{i,\sigma(i)} = 1, \forall i\}| = \sum_{\sigma \in S_n} \prod_{i=1}^n b_{i,\sigma(i)} = \operatorname{perm}(B).$$

!!! example "例 40A.18"
    完全二部图 $K_{3,3}$ 的邻接矩阵为 $J_3$。$\operatorname{perm}(J_3) = 3! = 6$，即 $K_{3,3}$ 有 6 个完美匹配，一一对应于 $S_3$ 的 6 个置换。

!!! example "例 40A.19"
    二部图 $G$（$n = 3$），邻接矩阵
    $$B = \begin{pmatrix} 1 & 1 & 0 \\ 1 & 0 & 1 \\ 0 & 1 & 1 \end{pmatrix}.$$
    逐一验证 $S_3$ 中 6 个置换：

    - $\sigma = \text{id} = (1,2,3)$: $b_{11}b_{22}b_{33} = 1 \cdot 0 \cdot 1 = 0$
    - $\sigma = (1,3,2)$: $b_{11}b_{23}b_{32} = 1 \cdot 1 \cdot 1 = 1$
    - $\sigma = (2,1,3)$: $b_{12}b_{21}b_{33} = 1 \cdot 1 \cdot 1 = 1$
    - $\sigma = (2,3,1)$: $b_{12}b_{23}b_{31} = 1 \cdot 1 \cdot 0 = 0$
    - $\sigma = (3,1,2)$: $b_{13}b_{21}b_{32} = 0 \cdot 1 \cdot 1 = 0$
    - $\sigma = (3,2,1)$: $b_{13}b_{22}b_{31} = 0 \cdot 0 \cdot 0 = 0$

    $\operatorname{perm}(B) = 2$。$G$ 恰有 2 个完美匹配。

### BosonSampling

!!! note "注记 40A.6 (玻色子采样与量子计算)"
    2011 年，Aaronson 和 Arkhipov 提出了**玻色子采样**（BosonSampling）问题：模拟 $n$ 个全同玻色子（如光子）通过 $m$ 模线性光学网络的输出分布。

    **物理设置**：$n$ 个光子输入一个由 $m \times m$ 酉矩阵 $U$ 描述的线性光学干涉仪。输出状态的概率幅与 $U$ 的某个 $n \times n$ 子矩阵的积和式成正比：
    $$\Pr(\text{output } \mathbf{s}) \propto \frac{|\operatorname{perm}(U_S)|^2}{\prod_i s_i!},$$
    其中 $U_S$ 是从 $U$ 中选取特定行和列构成的子矩阵，$\mathbf{s}$ 为输出光子数分布。

    **计算意义**：Aaronson 和 Arkhipov 证明，若经典计算机能够高效模拟 BosonSampling 的输出分布（即使是近似采样），则多项式层次将坍缩到第三层——这在复杂性理论中被认为是极不可能的。

    因此，BosonSampling 为**量子优越性**（quantum advantage）提供了一条可能的途径：即使在缺乏通用量子计算机的情况下，线性光学系统也可能执行经典计算机无法高效完成的采样任务。

    2020 年，中国科学技术大学的"九章"量子计算原型机首次在实验上演示了具有量子计算优势的高斯玻色采样，使用了 76 个光子。

!!! note "注记 40A.7 (高斯 BosonSampling 与 Hafnian)"
    高斯玻色采样（Gaussian BosonSampling）是 BosonSampling 的变体，使用压缩态光源而非单光子源。在这种情形中，输出概率与 Hafnian（而非积和式）相关：
    $$\Pr(\text{output}) \propto |\operatorname{hf}(A_S)|^2,$$
    其中 $A_S$ 是与高斯态协方差矩阵相关的某个子矩阵。

    由于 Hafnian 的计算同样是 $\#P$-完全的，高斯 BosonSampling 也提供了量子优越性的论证。

---

## 40A.10 习题

!!! exercise "习题 40A.1"
    计算以下矩阵的积和式：

    (a) $A = \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 1 & 1 & 0 \end{pmatrix}$

    (b) $A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 2 & 1 \\ 0 & 1 & 2 \end{pmatrix}$

    (c) $A = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & 1 & 1 & 1 \\ 1 & 1 & 1 & 1 \\ 1 & 1 & 1 & 1 \end{pmatrix}$

!!! exercise "习题 40A.2"
    证明：若 $A$ 为 $n \times n$ 上三角矩阵，则 $\operatorname{perm}(A) = \prod_{i=1}^n a_{ii}$（与行列式相同）。

!!! exercise "习题 40A.3"
    设 $A$ 为 $n \times n$ 非负矩阵，$D = \operatorname{diag}(a_{11}, \ldots, a_{nn})$。证明 $\operatorname{perm}(A) \ge \operatorname{perm}(D) = \prod a_{ii}$。

!!! exercise "习题 40A.4"
    用 Ryser 公式计算 $\operatorname{perm}\begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix}$，并验证结果。

!!! exercise "习题 40A.5"
    设 $A \in \Omega_3$（$3 \times 3$ 双随机矩阵），$A = \begin{pmatrix} 1/3 & 1/3 & 1/3 \\ 1/3 & 1/3 & 1/3 \\ 1/3 & 1/3 & 1/3 \end{pmatrix}$。验证 $\operatorname{perm}(A) = 2/9 = 3!/3^3$。

!!! exercise "习题 40A.6"
    证明：对 $n \times n$ 置换矩阵 $P$ 和任意矩阵 $A$，
    $$\operatorname{perm}(PA) = \operatorname{perm}(A), \quad \operatorname{perm}(AP) = \operatorname{perm}(A).$$

!!! exercise "习题 40A.7"
    设 $A$ 为 $n \times n$ $(0,1)$-矩阵，每行之和为 $r$，每列之和为 $r$（$r$-正则）。

    (a) 证明 $\operatorname{perm}(A) \ge 1$（提示：Hall 定理）。

    (b) 利用 Van der Waerden 猜想证明 $\operatorname{perm}(A) \ge r^n \cdot n!/n^n$。

    (c) 利用 Bregman-Minc 不等式证明 $\operatorname{perm}(A) \le (r!)^{n/r}$。

!!! exercise "习题 40A.8"
    (混合判别式) 设 $A_1 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$A_2 = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$。

    (a) 计算 $D(A_1, A_2)$ 和 $D(A_1, A_1)$。

    (b) 验证 Alexandrov 不等式 $D(A_1, A_2)^2 \ge D(A_1, A_1) \cdot D(A_2, A_2)$。

!!! exercise "习题 40A.9"
    证明 Hadamard-积和式不等式的等号条件：若 $A \ge 0$ 且 $\operatorname{perm}(A) = \prod_i r_i$（$r_i$ 为行和），则 $A$ 的每一行至多有一个非零元素。

!!! exercise "习题 40A.10"
    设 $A$ 为 $2n \times 2n$ 对称矩阵，$a_{ii} = 0$。

    (a) 计算 $K_4$（完全图）的 Hafnian。

    (b) 证明 Petersen 图有恰好 2000 个完美匹配。（提示：Petersen 图有 10 个顶点。）

!!! exercise "习题 40A.11"
    (Gurvits 容量) 对 $2 \times 2$ 双随机矩阵 $A_t = \begin{pmatrix} t & 1-t \\ 1-t & t \end{pmatrix}$（$0 \le t \le 1$），计算生成多项式 $p_{A_t}(x_1, x_2)$ 的容量，并验证 $\operatorname{cap}(p_{A_t}) \le \operatorname{perm}(A_t)$。

!!! exercise "习题 40A.12"
    证明：计算 $n \times n$ 矩阵（元素为 $\{0, 1\}$）的积和式不能用少于 $2^{n-1}$ 次加减法完成。（提示：考虑积和式作为多重线性函数的单项式数量。）

!!! exercise "习题 40A.13"
    (London 不等式的验证) 对 $n = 3$，取 $B$ 为 $3 \times 3$ 置换矩阵 $\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$。

    (a) 计算 $A_\alpha = \alpha B + (1-\alpha) J_3/3$ 的积和式（作为 $\alpha$ 的函数）。

    (b) 验证 $\operatorname{perm}(A_\alpha)$ 在 $[0, 1]$ 上关于 $\alpha$ 单调递增。

!!! exercise "习题 40A.14"
    设 $A$ 为 $n \times n$ 半正定 Hermite 矩阵，特征值为 $\lambda_1 \ge \cdots \ge \lambda_n \ge 0$。

    (a) 证明 $\operatorname{perm}(A) \ge \prod_{i=1}^n \lambda_i = \det(A)$（Marcus-Merris 不等式的特征值形式）。

    (b) 证明 $\operatorname{perm}(A) \le \prod_{i=1}^n (\operatorname{tr}(A) - \lambda_i + \lambda_i) = (\operatorname{tr} A)^n / n^n \cdot n!$ 不成立。找到正确的上界。

!!! exercise "习题 40A.15"
    (BosonSampling) 设 $U = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$（Hadamard 矩阵）。

    两个光子从输入端口 1 和 2 进入。计算两个光子分别从端口 1 和 2 输出的概率（与 $|\operatorname{perm}(U)|^2$ 成正比）。解释 Hong-Ou-Mandel 效应。

## 练习题

1. **[对比] 计算 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ 的积和式 $\operatorname{perm}(A)$。它与行列式 $\det(A)$ 有什么区别？**
   ??? success "参考答案"
       $\operatorname{perm}(A) = 1 \times 4 + 2 \times 3 = 10$。
       $\det(A) = 1 \times 4 - 2 \times 3 = -2$。
       积和式通过取消求和项中的符号因子 $\operatorname{sgn}(\sigma)$，使得所有项均为正（对非负阵），但失去了 Gauss 消元的代数简化性质。

2. **[组合] 证明：对 $(0,1)$-矩阵，积和式等于对应二部图的完美匹配数量。**
   ??? success "参考答案"
       积和式展开式中的每一项 $\prod a_{i, \sigma(i)}$ 只有在所有选中的 $a_{i, \sigma(i)}$ 均为 1 时才为 1，否则为 0。在二部图中，这意味着每个左部顶点 $i$ 必须连接到一个唯一的右部顶点 $\sigma(i)$，且该边必须存在。这正是完美匹配的定义。

3. **[不变性] 证明：交换矩阵的任意两行（或两列），积和式的值保持不变。**
   ??? success "参考答案"
       交换两行相当于对所有置换 $\sigma$ 进行一次对换复合。由于置换群 $S_n$ 对复合运算封闭，且积和式定义中没有 $\operatorname{sgn}(\sigma)$ 因子，求和项只是重新排列，总和不变。

4. **[Van der Waerden] 对于 $3 \times 3$ 的双随机矩阵，其积和式的最小值是多少？在什么矩阵处取得？**
   ??? success "参考答案"
       最小值为 $3!/3^3 = 6/27 = 2/9 \approx 0.222$。在所有元素均为 $1/3$ 的均匀矩阵 $J_3/3$ 处取得。

5. **[Marcus-Merris] 验证矩阵 $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ 满足 $\operatorname{perm}(A) \ge \det(A)$。**
   ??? success "参考答案"
       $\operatorname{perm}(A) = 4 + 1 = 5$。
       $\det(A) = 4 - 1 = 3$。
       显然 $5 \ge 3$。对于半正定矩阵，积和式总是大于或等于行列式。

6. **[计算] 利用 Ryser 公式计算 $\operatorname{perm}\begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 1 \\ 1 & 1 & 1 \end{pmatrix}$。**
   ??? success "参考答案"
       根据 Ryser 公式：$\operatorname{perm}(A) = (-1)^n \sum_{S \subseteq \{1,\dots,n\}} (-1)^{|S|} (\text{行和之积})_S$。
       对于 $J_3$，子集大小为 $k$ 的项对应的行和之积为 $k^3$。
       $\operatorname{perm}(J_3) = -[0 - 3(1)^3 + 3(2)^3 - 1(3)^3] = -[-3 + 24 - 27] = -[-6] = 6$。
       验证：$3! = 6$。正确。

7. **[复杂性] 为什么我们不能用 Gauss 消元法来计算积和式？**
   ??? success "参考答案"
       因为 Gauss 消元的核心是“行倍加变换保持行列式不变”，这依赖于交换行变号以及两行相同行列式为零的性质。积和式不具备这些性质（交换行值不变，两行相同值通常不为零），因此任何试图消元的尝试都会彻底改变积和式的值。

8. **[Hafnian] Hafnian 与积和式有什么联系？**
   ??? success "参考答案"
       Hafnian 是将完美匹配计数推广到非二部图的工具。如果将一个 $n \times n$ 的积和式问题嵌入到一个特定的 $2n \times 2n$ 对称块矩阵中（见定理 40A.14），该矩阵的 Hafnian 就精确等于原矩阵的积和式。

9. **[BosonSampling] 简述为什么计算积和式的难度与量子计算的优越性有关。**
   ??? success "参考答案"
       在玻色子采样实验中，全同玻色子的输出概率幅由干涉仪矩阵的子矩阵积和式决定。由于积和式是 $\#P$-完全的（极其难算），而量子系统可以通过自然的演化“算出”这个结果，这种实验如果在大规模下成功演示，就证明了量子系统在处理特定计算任务上远超经典计算机。

10. **[爱因斯坦思考题] 玻色-爱因斯坦统计描述的是“不可分辨”的玻色子。为什么不可分辨性会导致求和式中 $\operatorname{sgn}(\sigma)$ 消失，从而将计算从简单的行列式（费米子）变成了困难的积和式（玻色子）？这反映了自然界怎样的统计对称性差异？**
    ??? success "参考答案"
        费米子满足泡利不相容原理，其波函数在交换粒子时是反对称的（由行列式描述），这种相消相长使得系统呈现出某种“排他性”秩序，代数上易于处理。而玻色子是全同且社交的，其波函数在交换粒子时完全对称（由积和式描述），这意味着所有可能的交换路径都要以“正向叠加”的方式累加。这种累加消除了抵消机制，使得系统的微观态空间呈指数级复杂化。从爱因斯坦的视角看，宇宙通过对称性的选择（费米 vs 玻色），在最底层决定了信息的计算复杂性：排他性带来高效（行列式），而纠缠与共存带来复杂（积和式）。

## 本章小结

本章深入研究了行列式的“对称版”——积和式（Permanent），揭示了其在代数、组合与计算领域的独特地位：

1. **基本定义与性质**：确立了积和式作为不带符号因子的置换和，强调了其多重线性与置换不变性。
2. **计算鸿沟**：通过 Valiant 的 $\#P$-完全性定理，阐明了积和式在计算上远比行列式困难的本质原因，并介绍了 Ryser 公式等精确算法。
3. **极值与不等式**：详细讨论了 Van der Waerden 猜想的证明及其几何意义，以及 Marcus-Merris、Bregman-Minc 等限制性边界。
4. **组合意义**：将积和式与二部图完美匹配计数建立了精确的一一对应关系。
5. **现代物理关联**：展示了积和式如何作为玻色子统计的数学语言，成为当代量子优越性论证（如玻色采样）的核心工具。

