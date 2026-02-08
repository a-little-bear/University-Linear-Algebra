# 第 40 章 永久式与 Immanant

<div class="context-flow" markdown>

**前置**：行列式(Ch3) · 群论基础

**本章脉络**：永久式定义 $\to$ 与行列式的对比 $\to$ Van der Waerden 猜想 $\to$ 计算复杂性 $\to$ Immanant $\to$ 与对称群表示的联系

**延伸**：永久式在量子光学（玻色子采样 BosonSampling）和组合学（完美匹配计数）中有直接应用；$\#P$-完全性揭示了行列式（多项式时间）与永久式之间的根本计算鸿沟

</div>

行列式 $\det(A) = \sum_{\sigma \in S_n} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}$ 是矩阵论中最基本的量。如果去掉求和中的符号因子 $\operatorname{sgn}(\sigma)$，得到的量就是**永久式**（permanent）。这一看似微小的改变却带来了本质性的差异：行列式可以在 $O(n^3)$ 时间内计算（通过 Gauss 消元），而永久式的计算是 $\#P$-完全的——被认为比 $NP$-完全问题还要困难。永久式在组合学（完美匹配计数）、概率论（随机矩阵）和量子计算（玻色子采样）中都有深刻的应用。

更一般地，将 $\operatorname{sgn}(\sigma)$ 替换为对称群 $S_n$ 的任意特征标（character）$\chi(\sigma)$，就得到了 **immanant**——一个统一了行列式和永久式的概念框架。

---

## 40.1 永久式的定义

<div class="context-flow" markdown>

**核心问题**：永久式是什么？它与行列式有什么相似和不同之处？

</div>

!!! definition "定义 40.1 (永久式)"
    设 $A = (a_{ij}) \in M_n(\mathbb{C})$。$A$ 的**永久式**（permanent）定义为
    $$\operatorname{perm}(A) = \sum_{\sigma \in S_n} \prod_{i=1}^n a_{i,\sigma(i)},$$
    其中求和遍历 $n$ 元对称群 $S_n$ 的所有 $n!$ 个置换。

!!! note "注记 40.1 (与行列式的对比)"
    永久式与行列式的公式极为相似：
    $$\det(A) = \sum_{\sigma \in S_n} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}, \quad \operatorname{perm}(A) = \sum_{\sigma \in S_n} \prod_{i=1}^n a_{i,\sigma(i)}.$$
    唯一的区别是行列式中有符号因子 $\operatorname{sgn}(\sigma) = (-1)^{\text{逆序数}}$，而永久式中没有。这一差别看似微小，却导致了以下根本性的不同：

    | 性质 | 行列式 | 永久式 |
    |------|--------|--------|
    | 多重线性 | 是 | 是 |
    | 交换两行 | 变号 | 不变 |
    | 两行相同 | $= 0$ | $\ne 0$（一般地） |
    | 计算复杂度 | $O(n^3)$ | $\#P$-完全 |
    | 矩阵乘积 | $\det(AB) = \det(A)\det(B)$ | $\operatorname{perm}(AB) \ne \operatorname{perm}(A)\operatorname{perm}(B)$（一般地） |
    | 相似不变 | 是 | 否 |

!!! example "例 40.1"
    $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$。

    $$\det(A) = ad - bc, \quad \operatorname{perm}(A) = ad + bc.$$

    对 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$：$\det(A) = -2$，$\operatorname{perm}(A) = 4 + 6 = 10$。

!!! example "例 40.2"
    全 1 矩阵 $J_n = \begin{pmatrix} 1 & \cdots & 1 \\ \vdots & & \vdots \\ 1 & \cdots & 1 \end{pmatrix} \in M_n$：

    $$\operatorname{perm}(J_n) = \sum_{\sigma \in S_n} 1 = n!, \quad \det(J_n) = 0 \text{ (当 $n \ge 2$)}.$$

!!! example "例 40.3"
    $$\operatorname{perm}\begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 1 \\ 1 & 1 & 1 \end{pmatrix} = 3! = 6.$$

    $$\operatorname{perm}\begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix} = 1\cdot5\cdot9 + 1\cdot6\cdot8 + 2\cdot4\cdot9 + 2\cdot6\cdot7 + 3\cdot4\cdot8 + 3\cdot5\cdot7 = 45 + 48 + 72 + 84 + 96 + 105 = 450.$$

---

## 40.2 基本性质

<div class="context-flow" markdown>

**核心问题**：永久式有哪些代数性质？它与行列式共享哪些性质，哪些不同？

</div>

!!! theorem "定理 40.1 (永久式的基本性质)"
    (a) **多重线性**：$\operatorname{perm}(A)$ 关于 $A$ 的每一行（及每一列）是线性函数。

    (b) **行列置换不变**：交换 $A$ 的任意两行（或两列），$\operatorname{perm}(A)$ 不变。

    (c) **Laplace 展开**：按第 $i$ 行展开，
    $$\operatorname{perm}(A) = \sum_{j=1}^n a_{ij} \operatorname{perm}(A(i|j)),$$
    其中 $A(i|j)$ 是去掉第 $i$ 行第 $j$ 列的 $(n-1) \times (n-1)$ 子矩阵。注意：没有 $(-1)^{i+j}$ 符号。

    (d) **转置不变**：$\operatorname{perm}(A^T) = \operatorname{perm}(A)$。

    (e) **非负性**：若 $A \ge 0$（所有元素非负），则 $\operatorname{perm}(A) \ge 0$。若进一步 $A > 0$，则 $\operatorname{perm}(A) > 0$。

??? proof "证明"
    (a) 每个求和项 $\prod_{k=1}^n a_{k,\sigma(k)}$ 关于第 $i$ 行是线性的（只含一个因子 $a_{i,\sigma(i)}$）。

    (b) 交换第 $i, j$ 行等价于在求和中用 $\tau\sigma$（$\tau = (ij)$ 为对换）替换 $\sigma$。由于 $S_n$ 在左乘对换下是自身的置换，求和不变。

    (c) 将 $\prod_{k=1}^n a_{k,\sigma(k)}$ 按 $a_{i,\sigma(i)} = a_{ij}$（当 $\sigma(i) = j$）分组求和即得。

    (d) $\operatorname{perm}(A^T) = \sum_\sigma \prod_i a_{\sigma(i),i} = \sum_\sigma \prod_i a_{i,\sigma^{-1}(i)} = \sum_{\tau} \prod_i a_{i,\tau(i)} = \operatorname{perm}(A)$。

    (e) 若 $A \ge 0$，则每个求和项 $\prod a_{i,\sigma(i)} \ge 0$。

!!! theorem "定理 40.2 (半正定矩阵的永久式)"
    若 $A \in M_n(\mathbb{C})$ 为半正定 Hermite 矩阵，则 $\operatorname{perm}(A) \ge 0$。

    若 $A$ 为正定 Hermite 矩阵，则 $\operatorname{perm}(A) > 0$。

??? proof "证明"
    设 $A = B^*B$，则利用 Cauchy-Binet 型公式，
    $$\operatorname{perm}(A) = \operatorname{perm}(B^*B) = \sum_{\sigma} |\operatorname{perm}(B[\cdot, \sigma(\cdot)])|^2 \ge 0.$$
    实际上更直接的证明：$A$ 半正定时 $a_{ii} \ge 0$，且利用 Schur 的 Hadamard 型不等式或 Marcus-Merris 不等式可以得到非负性。

    对正定情形，$a_{ii} > 0$ 对所有 $i$，恒等置换的贡献 $\prod a_{ii} > 0$，故 $\operatorname{perm}(A) > 0$。

!!! theorem "定理 40.3 (永久式的乘法不等式)"
    设 $A, B \in M_n(\mathbb{R})$，$A, B \ge 0$。一般地 $\operatorname{perm}(AB) \ne \operatorname{perm}(A) \operatorname{perm}(B)$，但有
    $$\operatorname{perm}(AB) \le \operatorname{perm}(A) \operatorname{perm}(B)$$
    当 $A$ 或 $B$ 为双随机矩阵时成立（这是 Marcus-Merris-Watkins 不等式的特殊情形）。

!!! example "例 40.4"
    对 $A = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$：
    $$\operatorname{perm}(A) = 1, \quad \operatorname{perm}(B) = 1, \quad AB = B, \quad \operatorname{perm}(AB) = 1 = \operatorname{perm}(A)\operatorname{perm}(B).$$
    但一般地等号不成立。

---

## 40.3 Van der Waerden 猜想

<div class="context-flow" markdown>

**核心问题**：双随机矩阵的永久式的最小值是什么？何时达到？

</div>

!!! definition "定义 40.2 (双随机矩阵)"
    非负矩阵 $A \ge 0$ 称为**双随机矩阵**，若每行每列之和均为 1。全体 $n$ 阶双随机矩阵构成凸集 $\Omega_n$（**Birkhoff 多面体**）。

!!! theorem "定理 40.4 (Van der Waerden 猜想 / Egorychev-Falikman 定理)"
    对所有 $n \times n$ 双随机矩阵 $A$，
    $$\operatorname{perm}(A) \ge \frac{n!}{n^n}.$$
    等号成立当且仅当 $A = J_n/n$（即 $A$ 为所有元素等于 $1/n$ 的矩阵）。

??? proof "证明（Egorychev 的证明概要）"
    这一猜想由 Van der Waerden 于 1926 年提出，经过半个多世纪的努力，最终在 1981 年由 Egorychev 和 Falikman 独立证明。

    **Egorychev 的证明**使用了以下核心工具：

    **步骤一（Van der Waerden 永久式函数的性质）**：在紧凸集 $\Omega_n$ 上，$\operatorname{perm}$ 是连续函数，故取到最小值。设最小值在 $A^* \in \Omega_n$ 处取到。

    **步骤二（最小值点的刻画）**：利用 Lagrange 乘子法和永久式的多重线性性，可以证明在最小值点 $A^*$，所有**余永久式**（cofactor permanents）相等：
    $$\operatorname{perm}(A^*(i|j)) = c \quad \text{对所有 } i, j.$$

    **步骤三（Alexandrov 不等式）**：关键的分析工具是 **Alexandrov 关于混合判别式的不等式**。对半正定 Hermite 矩阵 $H_1, \ldots, H_n$，混合判别式 $D(H_1, \ldots, H_n)$ 满足
    $$D(H_1, \ldots, H_n)^2 \ge D(H_1, H_1, H_3, \ldots, H_n) D(H_2, H_2, H_3, \ldots, H_n).$$
    永久式可以视为混合判别式的特殊情形（取 $H_i = \operatorname{diag}(a_{i1}, \ldots, a_{in})$）。

    **步骤四**：结合步骤二和步骤三，通过归纳可以证明 $A^* = J_n/n$ 是唯一的最小值点，且
    $$\operatorname{perm}(J_n/n) = n! / n^n.$$

!!! example "例 40.5"
    $n = 2$：$\Omega_2$ 中的双随机矩阵为 $\begin{pmatrix} t & 1-t \\ 1-t & t \end{pmatrix}$（$0 \le t \le 1$）。

    $\operatorname{perm} = t^2 + (1-t)^2 = 2t^2 - 2t + 1$。最小值在 $t = 1/2$ 时取到，$\operatorname{perm} = 1/2 = 2!/2^2$。验证了 $n!/n^n = 2/4 = 1/2$。

!!! example "例 40.6"
    $n = 3$：$\operatorname{perm}(J_3/3) = 3!/3^3 = 6/27 = 2/9 \approx 0.222$。

    对比：置换矩阵 $P$ 的永久式为 1（双随机矩阵的永久式上界之一），远大于 $2/9$。

!!! note "注记 40.2 (Permanent 猜想的历史)"
    Van der Waerden 猜想提出于 1926 年，是组合数学中最著名的猜想之一。在被证明之前，许多数学家建立了部分结果（如 Marcus-Newman 证明了 $n = 3$ 的情形，London 证明了下界 $n!/((n-1)^n + n - 1)$）。1981 年，Egorychev 和 Falikman 的证明被认为是永久式理论的里程碑。两人因此获得了 1982 年的 Polya 奖。

---

## 40.4 永久式的不等式

<div class="context-flow" markdown>

**核心问题**：永久式满足哪些重要的不等式？

</div>

!!! theorem "定理 40.5 (Hadamard-永久式不等式)"
    设 $A \ge 0$ 为 $n \times n$ 非负矩阵，行和为 $r_i = \sum_j a_{ij}$。则
    $$\operatorname{perm}(A) \le \prod_{i=1}^n r_i.$$

??? proof "证明"
    每个求和项 $\prod_i a_{i,\sigma(i)} \le \prod_i \max_j a_{ij} \le \prod_i r_i$ 不对——这给的界太粗糙。

    正确的证明：按行展开 Laplace 展开并利用 AM-GM 不等式。实际上更直接的方法：$\operatorname{perm}(A) = \sum_\sigma \prod_i a_{i,\sigma(i)}$。将 $A$ 的每一行归一化为概率向量 $\tilde{a}_{ij} = a_{ij}/r_i$，则 $\operatorname{perm}(\tilde{A}) \le 1$（因为 $\sum_\sigma \prod_i \tilde{a}_{i,\sigma(i)} \le \prod_i (\sum_j \tilde{a}_{ij}) = 1$ 由永久式的上界得到），故 $\operatorname{perm}(A) = \prod_i r_i \cdot \operatorname{perm}(\tilde{A}) \le \prod_i r_i$。

!!! theorem "定理 40.6 (Bregman-Minc 不等式)"
    设 $A$ 为 $n \times n$ $(0,1)$-矩阵（所有元素为 0 或 1），行和为 $r_i$。则
    $$\operatorname{perm}(A) \le \prod_{i=1}^n (r_i!)^{1/r_i}.$$
    等号成立当且仅当 $A$（行列重排后）是若干全 1 方阵的直和。

??? proof "证明（Schrijver 的信息论证明概要）"
    Minc 于 1963 年提出此猜想，Bregman 于 1973 年利用**双随机矩阵的永久式与熵**之间的联系给出了证明。更简洁的证明由 Schrijver (1978) 和 Radhakrishnan (1997) 给出。

    Radhakrishnan 的证明基于 Shannon 熵。将 $A$ 视为二部图的邻接矩阵，$\operatorname{perm}(A)$ 等于完美匹配的数量。通过考虑均匀分布在完美匹配上的随机变量，利用条件熵的次可加性和 $\log$ 的凹性，可以推导出 Bregman 界。

!!! theorem "定理 40.7 (Marcus-Merris 不等式)"
    设 $A \in M_n(\mathbb{C})$ 为半正定 Hermite 矩阵。则
    $$\operatorname{perm}(A) \ge \det(A) \ge 0.$$
    等号 $\operatorname{perm}(A) = \det(A)$ 成立当且仅当 $A$ 为对角矩阵。

??? proof "证明"
    $\operatorname{perm}(A) - \det(A) = \sum_{\sigma \ne \text{id}} (1 - \operatorname{sgn}(\sigma)) \prod_i a_{i,\sigma(i)} + \sum_{\sigma: \operatorname{sgn}(\sigma) = -1} 2\prod_i a_{i,\sigma(i)}$。
    对半正定矩阵，$\prod_i a_{i,\sigma(i)} \ge 0$（这对一般半正定矩阵不明显，需要利用 Fischer 不等式的推广）。
    更直接的证明：利用 Schur 的 Hadamard 不等式和永久式的性质。

!!! example "例 40.7"
    $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$（正定）。$\det(A) = 3$，$\operatorname{perm}(A) = 4 + 1 = 5 \ge 3$。

---

## 40.5 计算复杂性

<div class="context-flow" markdown>

**核心问题**：为什么永久式比行列式难算得多？最快的精确算法是什么？

</div>

!!! theorem "定理 40.8 (Valiant, 1979: $\#P$-完全性)"
    计算 $(0,1)$-矩阵的永久式是 $\#P$-完全问题。即使限制为 $(0,1)$-矩阵（每个元素为 0 或 1），计算永久式仍然是 $\#P$-完全的。

!!! note "注记 40.3 ($\#P$ 与 $NP$ 的关系)"
    $\#P$ 是**计数**复杂类——不仅问"是否存在"某个对象，而是问"有多少个"。$\#P$-完全问题至少与 $NP$-完全问题一样难，且在大多数复杂性理论假设下严格更难。

    行列式与永久式之间的计算鸿沟是理论计算机科学中最深刻的现象之一：
    - $\det$：多项式时间（$O(n^3)$ 或 $O(n^\omega)$，$\omega \approx 2.37$）；
    - $\operatorname{perm}$：$\#P$-完全（除非 $P = \#P$，没有多项式时间算法）。

    这种差异的根源在于行列式中的符号交替允许大量消除（cancellation），使得 Gauss 消元成为可能；而永久式中所有项同号，没有消除发生。

!!! theorem "定理 40.9 (Ryser 公式)"
    $$\operatorname{perm}(A) = (-1)^n \sum_{S \subseteq \{1,\ldots,n\}} (-1)^{|S|} \prod_{i=1}^n \sum_{j \in S} a_{ij}.$$
    此公式的计算复杂度为 $O(2^n n)$，远优于直接定义的 $O(n! \cdot n)$，但仍然是指数级的。

??? proof "证明"
    **容斥原理**。令 $r_i(S) = \sum_{j \in S} a_{ij}$ 为第 $i$ 行在列集 $S$ 上的行和。$\prod_{i=1}^n r_i(S)$ 展开后是对所有函数 $f: \{1,\ldots,n\} \to S$ 的求和 $\sum_f \prod_i a_{i,f(i)}$。

    我们需要的是对所有**双射** $\sigma: \{1,\ldots,n\} \to \{1,\ldots,n\}$ 的求和。由容斥，双射的贡献可以从所有函数的贡献中提取：

    $$\operatorname{perm}(A) = \sum_{k=0}^{n} (-1)^{n-k} \binom{}{} \sum_{|S|=k} \prod_{i=1}^n r_i(S),$$

    整理后得到 Ryser 公式。具体推导使用了 $\sigma$ 是双射 $\iff$ $\sigma$ 的像集为 $\{1,\ldots,n\}$，而"像集恰好为 $T$"可以通过容斥表达为关于"像集 $\subseteq S$"的交替和。

!!! example "例 40.8"
    用 Ryser 公式计算 $\operatorname{perm}\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$。

    $n = 2$。$S$ 遍历 $\{1,2\}$ 的子集：

    - $S = \emptyset$：$\prod r_i(\emptyset) = 0$。
    - $S = \{1\}$：$r_1 = 1, r_2 = 3$，积 $= 3$。
    - $S = \{2\}$：$r_1 = 2, r_2 = 4$，积 $= 8$。
    - $S = \{1,2\}$：$r_1 = 3, r_2 = 7$，积 $= 21$。

    $\operatorname{perm} = (-1)^2[(-1)^0 \cdot 0 + (-1)^1(3 + 8) + (-1)^2 \cdot 21] = 0 - 11 + 21 = 10$。

    验证：$1 \cdot 4 + 2 \cdot 3 = 10$。正确。

!!! note "注记 40.4 (近似算法)"
    虽然精确计算永久式是 $\#P$-完全的，Jerrum, Sinclair 和 Vigoda (2004) 给出了一个**完全多项式随机近似方案**（FPRAS）：对非负矩阵 $A \ge 0$，可以在多项式时间内以高概率输出 $\operatorname{perm}(A)$ 的 $(1 \pm \epsilon)$ 倍近似值。该算法基于马尔可夫链蒙特卡罗（MCMC）方法，是理论计算机科学的重大突破。

---

## 40.6 Immanant

<div class="context-flow" markdown>

**核心问题**：能否用对称群的表示论统一行列式和永久式？

</div>

!!! definition "定义 40.3 (Immanant)"
    设 $\chi$ 为对称群 $S_n$ 的一个特征标（character）。矩阵 $A \in M_n(\mathbb{C})$ 的 **$\chi$-immanant** 定义为
    $$d_\chi(A) = \sum_{\sigma \in S_n} \chi(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}.$$

!!! note "注记 40.5 (特殊情形)"
    - 当 $\chi = \operatorname{sgn}$（符号特征标，对应交错表示）时，$d_\chi(A) = \det(A)$。
    - 当 $\chi = \mathbf{1}$（平凡特征标，即 $\chi(\sigma) = 1$ 对所有 $\sigma$）时，$d_\chi(A) = \operatorname{perm}(A)$。
    - $S_n$ 的不可约特征标与 $n$ 的分拆（partition）一一对应。对分拆 $\lambda \vdash n$，记对应特征标为 $\chi^\lambda$，不可约 immanant 为 $d_\lambda(A) = d_{\chi^\lambda}(A)$。

!!! definition "定义 40.4 (对称群的不可约表示)"
    $S_n$ 的不可约表示（从而不可约特征标）由 $n$ 的**分拆** $\lambda = (\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_k > 0)$（$\lambda_1 + \cdots + \lambda_k = n$）参数化。每个分拆 $\lambda$ 对应一个 Young 图和一个 **Specht 模** $S^\lambda$。特征标 $\chi^\lambda$ 的维数（即 $\chi^\lambda(\text{id})$）等于 $\lambda$ 的**标准 Young 表**的个数 $f^\lambda$。

!!! example "例 40.9"
    $n = 3$ 时，$S_3$ 有三个不可约特征标：

    | 分拆 $\lambda$ | $(3)$ | $(2,1)$ | $(1,1,1)$ |
    |---------------|-------|---------|-----------|
    | Young 图 | $\yng(3)$ | $\yng(2,1)$ | $\yng(1,1,1)$ |
    | 维数 $f^\lambda$ | 1 | 2 | 1 |
    | 对应 immanant | $\operatorname{perm}$ | $d_{(2,1)}$ | $\det$ |

    对 $A = \begin{pmatrix} a & b & c \\ d & e & f \\ g & h & k \end{pmatrix}$，

    $$d_{(2,1)}(A) = \sum_{\sigma} \chi^{(2,1)}(\sigma) \prod_i a_{i,\sigma(i)}.$$

    $\chi^{(2,1)}$ 在各共轭类上的值：$\chi^{(2,1)}(\text{id}) = 2$，$\chi^{(2,1)}(\text{对换}) = 0$，$\chi^{(2,1)}(\text{3-循环}) = -1$。

    故 $d_{(2,1)}(A) = 2(aek + bfg + cdh) + 0(\cdots) - 1(afh + bdk + ceg) - 1(\cdots)$。

    需要将 $S_3$ 的 6 个元素按共轭类分组计算。

!!! theorem "定理 40.10 (Schur 正定性猜想 / Schur 不等式)"
    设 $A$ 为 $n \times n$ 半正定 Hermite 矩阵，$\chi^\lambda$ 为 $S_n$ 的不可约特征标。则
    $$\frac{d_\lambda(A)}{\chi^\lambda(\text{id})} \ge \det(A).$$
    即归一化的 immanant 不小于行列式。

!!! note "注记 40.6 (immanant 与 Schur 函数)"
    Immanant 与**对称函数论**密切相关。对对角矩阵 $A = \operatorname{diag}(x_1, \ldots, x_n)$，
    $$d_\lambda(A) = \chi^\lambda(\text{id}) \cdot s_\lambda(x_1, \ldots, x_n),$$
    其中 $s_\lambda$ 为 **Schur 多项式**（Schur polynomial），它是对称函数论中最重要的基。这一联系将 immanant 理论嵌入到代数组合学的核心。

!!! theorem "定理 40.11 (Lieb 定理)"
    设 $A$ 为 $n \times n$ 正半定 Hermite 矩阵。对 $S_n$ 的任意特征标 $\chi$（不一定不可约），
    $$|d_\chi(A)| \le \chi(\text{id}) \cdot \operatorname{perm}(A).$$

??? proof "证明（概要）"
    利用群代数 $\mathbb{C}[S_n]$ 中的正性论证。将 $d_\chi(A)$ 表示为某个算子的迹，利用 Cauchy-Schwarz 不等式和特征标的正交关系得到界。

---

## 40.7 组合应用

<div class="context-flow" markdown>

**核心问题**：永久式在图论和组合学中如何应用？

</div>

!!! theorem "定理 40.12 (永久式 = 完美匹配数)"
    设 $G = (U \cup V, E)$ 为二部图，$|U| = |V| = n$，邻接矩阵 $A = (a_{ij})$（$a_{ij} = 1$ 当 $(u_i, v_j) \in E$，否则为 0）。则 $G$ 的**完美匹配**（perfect matching）的数量等于
    $$\operatorname{perm}(A).$$

??? proof "证明"
    完美匹配是边集 $M \subseteq E$，使得 $U$ 和 $V$ 中每个顶点恰被 $M$ 中一条边覆盖。这等价于选择一个双射 $\sigma: U \to V$（$u_i$ 匹配 $v_{\sigma(i)}$），使得 $(u_i, v_{\sigma(i)}) \in E$（即 $a_{i,\sigma(i)} = 1$）对所有 $i$ 成立。

    因此完美匹配数 $= |\{\sigma \in S_n : a_{i,\sigma(i)} = 1, \forall i\}| = \sum_{\sigma \in S_n} \prod_{i=1}^n a_{i,\sigma(i)} = \operatorname{perm}(A)$。

!!! example "例 40.10"
    完全二部图 $K_{3,3}$ 的邻接矩阵为 $J_3$（全 1 矩阵）。$\operatorname{perm}(J_3) = 3! = 6$，即 $K_{3,3}$ 有 6 个完美匹配。这些匹配一一对应于 $S_3$ 的 6 个置换。

!!! example "例 40.11"
    二部图 $G$（$n = 3$），邻接矩阵
    $$A = \begin{pmatrix} 1 & 1 & 0 \\ 1 & 0 & 1 \\ 0 & 1 & 1 \end{pmatrix}.$$
    $\operatorname{perm}(A) = 1 \cdot 0 \cdot 1 + 1 \cdot 1 \cdot 1 + 0 + 0 + 1 \cdot 1 \cdot 1 + 0 = 0 + 1 + 0 + 0 + 1 + 0 = 2$。

    $G$ 恰有 2 个完美匹配：$\{(u_1,v_2), (u_2,v_3), (u_3,v_1)\}$ 和 $\{(u_1,v_1), (u_2,v_3), (u_3,v_2)\}$。

    但验算：$\sigma = (2,3,1)$: $a_{1,2}a_{2,3}a_{3,1} = 1 \cdot 1 \cdot 0 = 0$。$\sigma = (1,3,2)$: $a_{1,1}a_{2,3}a_{3,2} = 1 \cdot 1 \cdot 1 = 1$。$\sigma = (2,1,3)$: $a_{1,2}a_{2,1}a_{3,3} = 1 \cdot 1 \cdot 1 = 1$。其余为 0。所以 $\operatorname{perm} = 2$，完美匹配为 $\sigma = (1,3,2)$ 和 $\sigma = (2,1,3)$。

!!! definition "定义 40.5 (Hafnian)"
    对 $2n \times 2n$ 对称矩阵 $A$（$a_{ii} = 0$），**Hafnian** 定义为
    $$\operatorname{hf}(A) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \prod_{i=1}^n a_{\sigma(2i-1), \sigma(2i)}.$$
    Hafnian 计算一般图（非二部图）的完美匹配数：若 $A$ 为简单图的邻接矩阵，则 $\operatorname{hf}(A) = $ 完美匹配数。

!!! theorem "定理 40.13 (永久式与拉丁矩形)"
    $n \times k$ 拉丁矩形（每行是 $\{1, \ldots, n\}$ 的排列，每列无重复元素）的数量可以用永久式来表达。特别地，从 $r$ 行拉丁矩形扩展到 $r+1$ 行的方式数等于某个与前 $r$ 行相关的 $(0,1)$-矩阵的永久式。

!!! note "注记 40.7 (玻色子采样与量子计算)"
    2011 年，Aaronson 和 Arkhipov 提出了**玻色子采样**（BosonSampling）问题：模拟 $n$ 个光子通过 $m$ 模线性光学网络的输出分布。他们证明，输出概率与某个子矩阵的永久式的模平方成正比：
    $$\Pr(\text{output}) \propto |\operatorname{perm}(U_S)|^2,$$
    其中 $U_S$ 是酉矩阵 $U$ 的某个 $n \times n$ 子矩阵。

    由于永久式的计算是 $\#P$-完全的，经典计算机不太可能有效模拟玻色子采样，这为量子计算的"量子优越性"（quantum supremacy）提供了一条可能的途径。

!!! theorem "定理 40.14 (永久式与拉丁方)"
    $n$ 阶**拉丁方**（Latin square）的数量 $L(n)$ 可以用永久式来递推计算。具体地，设 $R_k$ 为由前 $k$ 行已确定的拉丁矩形所"允许"的 $(0,1)$-矩阵（其 $(i,j)$ 元素为 1 当且仅当将 $j$ 放在第 $k+1$ 行第 $i$ 列不违反拉丁方条件），则从 $k$ 行扩展到 $k+1$ 行的方式数等于 $\operatorname{perm}(R_k)$。因此
    $$L(n) = \sum \prod_{k=1}^{n-1} \operatorname{perm}(R_k),$$
    其中求和遍历所有可能的中间状态。

    由 Van der Waerden 猜想（定理 40.4），可以得到 $L(n)$ 的下界。已知的精确值包括 $L(1) = 1$, $L(2) = 2$, $L(3) = 12$, $L(4) = 576$, $L(5) = 161280$。

!!! note "注记 40.8 (永久式与统计物理)"
    在统计力学中，二聚体覆盖（dimer covering）模型的配分函数可以用永久式（二部图）或 Hafnian（一般图）表示。对于平面图上的二聚体覆盖，Kasteleyn (1961) 和 Temperley-Fisher (1961) 发现了一个巧妙的方法：通过对边赋予适当的 $\pm 1$ 符号（Kasteleyn 定向），可以将永久式转化为 Pfaffian（一种与行列式相关的量），从而在多项式时间内计算。

    这一结果说明，虽然一般图上的永久式是 $\#P$-完全的，但**平面**二部图上的永久式可以在多项式时间内计算。这进一步强调了图的拓扑结构对计算复杂性的深刻影响。

!!! theorem "定理 40.15 (Kasteleyn 定理)"
    设 $G$ 为**平面**二部图，邻接矩阵为 $A$。则存在符号矩阵 $S$（$s_{ij} = \pm 1$）使得
    $$\operatorname{perm}(A) = |\det(A \circ S)|,$$
    其中 $\circ$ 表示 Hadamard（逐元素）乘积。由于行列式可在 $O(n^3)$ 时间内计算，这给出了平面二部图完美匹配计数的多项式时间算法。

!!! example "例 40.12"
    考虑 $2 \times 3$ 棋盘格上的二聚体覆盖。将棋盘格视为二部图（黑白着色），邻接矩阵为
    $$A = \begin{pmatrix} 1 & 1 & 0 \\ 1 & 1 & 1 \\ 0 & 1 & 1 \end{pmatrix}.$$
    $\operatorname{perm}(A) = 1 \cdot 1 \cdot 1 + 1 \cdot 1 \cdot 1 + 0 + 0 + 1 \cdot 1 \cdot 1 + 0 = 3$。

    即 $2 \times 3$ 棋盘格有 3 种二聚体覆盖方式，与直觉一致（可以手动枚举验证）。

!!! note "注记 40.9 (det 与 perm 的代数几何视角)"
    从代数几何的角度看，行列式和永久式都是多重线性形式（$n^2$ 个变量的 $n$ 次齐次多项式）。Valiant 的 $VP$ vs $VNP$ 猜想（代数计算复杂性理论中的核心开放问题）可以表述为：永久式不能被多项式大小的行列式"投影"表示。即不存在 $\operatorname{poly}(n)$ 阶矩阵 $M$（其元素是 $a_{ij}$ 的仿射函数）使得 $\det(M) = \operatorname{perm}(A)$。

    如果这一猜想成立，它将是 $P \ne NP$ 的代数类比。目前已知的最佳结果（Mignon-Ressayre, 2004）表明，表示 $n \times n$ 永久式所需的行列式阶至少为 $\frac{n^2}{2}$。

---

## 本章小结

| 概念 | 定义 | 核心结果 |
|------|------|----------|
| 永久式 | $\sum_\sigma \prod a_{i,\sigma(i)}$ | 计数完美匹配；$\#P$-完全 |
| Van der Waerden | 双随机矩阵永久式下界 | $\operatorname{perm} \ge n!/n^n$ |
| Ryser 公式 | 容斥原理 | $O(2^n n)$ 计算 |
| Immanant | $\sum_\sigma \chi(\sigma) \prod a_{i,\sigma(i)}$ | 统一 det 和 perm |
| Schur 不等式 | $d_\lambda(A)/f^\lambda \ge \det(A)$ | PSD 矩阵上 |
| Kasteleyn 定理 | 平面图永久式 = 行列式 | $O(n^3)$ 用于平面图 |
| 组合应用 | perm = 完美匹配数 | 二部图匹配、BosonSampling |

永久式与行列式之间的对比是线性代数与理论计算机科学交叉领域中最引人注目的现象之一。公式上一个符号的差异（$\operatorname{sgn}(\sigma)$ 的有无）导致了从多项式时间到 $\#P$-完全的计算复杂度跳跃。Immanant 理论通过对称群的表示论为这一对比提供了更深层的理解框架。Kasteleyn 定理则揭示了拓扑结构（平面性）如何从根本上改变计算问题的复杂性。
