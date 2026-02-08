# 第 52 章 有限域上的线性代数

<div class="context-flow" markdown>

**前置**：向量空间 (Ch4) · 行列式 (Ch3) · 多项式代数 (Ch0)

**本章脉络**：有限域 $\mathrm{GF}(q)$ → $\mathrm{GF}(q)$ 上的向量空间 → 子空间计数（Gauss 二项式系数）→ $\mathrm{GL}(n,q)$ 的阶 → 有限域上的标准形 → 线性码 → 生成矩阵与校验矩阵

**延伸**：有限域线性代数是编码理论（线性码、Reed-Solomon 码）和密码学（AES 在 $\mathrm{GF}(2^8)$ 上运算、椭圆曲线密码）的数学基础；$q$-模拟将组合恒等式推广到 $\mathrm{GF}(q)$ 上的计数

</div>

实数域 $\mathbb{R}$ 和复数域 $\mathbb{C}$ 是线性代数最常用的标量域，但绝非唯一的选择。当我们将标量域取为**有限域**——只有有限个元素的域——线性代数的面貌发生了深刻的变化。虽然基本定理（维数公式、秩-零度定理等）依然成立，但新的**组合**现象涌现出来：子空间的个数是有限的（可以精确计数），一般线性群的阶可以用显式公式计算，而向量空间本身成为了信息传输中**纠错编码**的载体。

有限域 $\mathrm{GF}(q)$（其中 $q = p^n$，$p$ 为素数）上的线性代数是现代编码理论和密码学的数学核心。从 QR 码中的 Reed-Solomon 编码到 AES 加密算法中的 $\mathrm{GF}(2^8)$ 运算，有限域无处不在。

---

## 52.1 有限域回顾

<div class="context-flow" markdown>

**核心问题**：有限域的基本结构是什么？它们如何构造？

</div>

!!! definition "定义 52.1 (有限域)"
    **有限域**（finite field）是只有有限个元素的域。有限域也称为 **Galois 域**，记为 $\mathrm{GF}(q)$ 或 $\mathbb{F}_q$。

!!! theorem "定理 52.1 (有限域的存在与唯一性)"
    1. 有限域的元素个数必为素数幂：$q = p^n$，其中 $p$ 为素数，$n \geq 1$。
    2. 对每个素数幂 $q = p^n$，恰好存在一个（在同构意义下）$q$ 元域 $\mathrm{GF}(q)$。
    3. $\mathrm{GF}(q)$ 的特征为 $p$：$\underbrace{1 + 1 + \cdots + 1}_{p} = 0$。
    4. $\mathrm{GF}(q)$ 的乘法群 $\mathrm{GF}(q)^{\times} = \mathrm{GF}(q) \setminus \{0\}$ 是**循环群**，阶为 $q - 1$。

??? proof "证明要点"
    **(1)** 设 $F$ 是有限域，$\operatorname{char}(F) = p$（必为素数）。$F$ 包含素域 $\mathbb{F}_p = \mathbb{Z}/p\mathbb{Z}$ 作为子域。$F$ 是 $\mathbb{F}_p$ 上的有限维向量空间，设维数为 $n$，则 $|F| = p^n$。

    **(2)** $\mathrm{GF}(q)$ 可以构造为 $\mathbb{F}_p[x]/(f(x))$ 的分裂域，其中 $f(x) = x^q - x$。也可以取 $f(x)$ 为 $\mathbb{F}_p$ 上任何 $n$ 次不可约多项式。

    **(4)** 有限域的乘法群是有限 Abel 群。设其指数为 $e$（即最大阶元素的阶），则每个非零元素都是 $x^e - 1 = 0$ 的根。但 $x^e - 1$ 在域上至多有 $e$ 个根，故 $q - 1 \leq e$。另一方面 $e \mid q-1$（Lagrange 定理）。因此 $e = q-1$，乘法群是循环的。

!!! definition "定义 52.2 (有限域的构造)"
    **方法 1（素域）：** $\mathrm{GF}(p) = \mathbb{Z}/p\mathbb{Z} = \{0, 1, 2, \ldots, p-1\}$，加法和乘法模 $p$。

    **方法 2（扩域）：** 选取 $\mathbb{F}_p$ 上 $n$ 次不可约多项式 $f(x)$。则

    $$\mathrm{GF}(p^n) = \mathbb{F}_p[x]/(f(x)) = \{a_0 + a_1\alpha + \cdots + a_{n-1}\alpha^{n-1} : a_i \in \mathbb{F}_p\},$$

    其中 $\alpha = \bar{x}$ 是 $f$ 的根。加法逐系数进行，乘法后对 $f(\alpha) = 0$ 化简。

!!! example "例 52.1 (具体有限域)"
    **(a) $\mathrm{GF}(2) = \{0, 1\}$：** 加法表和乘法表：

    | $+$ | 0 | 1 |  | $\times$ | 0 | 1 |
    |:---:|:---:|:---:| |:---:|:---:|:---:|
    | 0 | 0 | 1 |  | 0 | 0 | 0 |
    | 1 | 1 | 0 |  | 1 | 0 | 1 |

    注意 $1 + 1 = 0$（特征 $2$）。

    **(b) $\mathrm{GF}(4)$：** $\mathbb{F}_2$ 上 $x^2 + x + 1$ 不可约。$\mathrm{GF}(4) = \{0, 1, \alpha, \alpha+1\}$，其中 $\alpha^2 = \alpha + 1$。

    乘法表（非零元素）：

    | $\times$ | 1 | $\alpha$ | $\alpha+1$ |
    |:---:|:---:|:---:|:---:|
    | 1 | 1 | $\alpha$ | $\alpha+1$ |
    | $\alpha$ | $\alpha$ | $\alpha+1$ | $1$ |
    | $\alpha+1$ | $\alpha+1$ | $1$ | $\alpha$ |

    $\mathrm{GF}(4)^{\times} = \{1, \alpha, \alpha^2 = \alpha+1\}$ 是 $3$ 阶循环群，$\alpha$ 是生成元。

    **(c) $\mathrm{GF}(2^8)$：** AES 加密使用的域。取不可约多项式 $f(x) = x^8 + x^4 + x^3 + x + 1$。每个元素可以表示为 $8$ 位二进制串（一个字节）。加法是逐位异或（XOR），乘法是多项式乘法模 $f(x)$。

!!! definition "定义 52.3 (Frobenius 自同构)"
    $\mathrm{GF}(p^n)$ 上的 **Frobenius 自同构**是映射

    $$\varphi: \mathrm{GF}(p^n) \to \mathrm{GF}(p^n), \quad \varphi(a) = a^p.$$

    这确实是域自同构：$\varphi(a+b) = (a+b)^p = a^p + b^p$（在特征 $p$ 中！），$\varphi(ab) = (ab)^p = a^p b^p$。

    $\mathrm{Gal}(\mathrm{GF}(p^n)/\mathrm{GF}(p)) = \langle \varphi \rangle \cong \mathbb{Z}/n\mathbb{Z}$（由 $\varphi$ 生成的 $n$ 阶循环群）。

---

## 52.2 $\mathrm{GF}(q)$ 上的向量空间

<div class="context-flow" markdown>

**核心问题**：有限域上的向量空间有多少元素？与实向量空间有什么相似与不同？

</div>

!!! definition "定义 52.4 ($\mathrm{GF}(q)^n$)"
    $\mathrm{GF}(q)$ 上的 $n$ 维向量空间

    $$\mathrm{GF}(q)^n = \{(x_1, x_2, \ldots, x_n) : x_i \in \mathrm{GF}(q)\}$$

    有 $q^n$ 个元素。标准基 $\{e_1, \ldots, e_n\}$ 与实向量空间相同。

!!! theorem "定理 52.2 (有限向量空间的基本性质)"
    $\mathrm{GF}(q)^n$ 的以下性质与一般向量空间相同：

    1. 维数定理：$\dim(U + W) + \dim(U \cap W) = \dim U + \dim W$；
    2. 秩-零度定理：$\dim \ker T + \dim \operatorname{im} T = \dim V$；
    3. 每个子空间都有补空间；
    4. 矩阵的行秩等于列秩。

    但有限性带来的新特征：

    5. $k$ 维子空间 $U \subseteq \mathrm{GF}(q)^n$ 恰有 $q^k$ 个元素。
    6. 子空间的个数是有限的（下节精确计算）。
    7. $\mathrm{GF}(q)^n$ 上不存在非退化的正定内积（当 $\operatorname{char} = 2$ 时，$\langle v, v \rangle = 0$ 对所有 $v$）。

!!! example "例 52.2 ($\mathrm{GF}(2)^3$ 的子空间)"
    $\mathrm{GF}(2)^3$ 有 $2^3 = 8$ 个元素。

    **$0$ 维子空间**：$\{(0,0,0)\}$，共 $1$ 个。

    **$1$ 维子空间**（过原点的"直线"）：每条直线包含 $\{0, v\}$（$v \neq 0$）。$7$ 个非零向量，每条直线包含 $2-1=1$ 个非零向量（$\mathrm{GF}(2)$ 中标量只有 $0, 1$），故共 $7$ 条直线。

    **$2$ 维子空间**（过原点的"平面"）：每个 $2$ 维子空间有 $2^2 = 4$ 个元素。共 $7$ 个（对偶性：$\binom{3}{1}_2 = \binom{3}{2}_2 = 7$）。

    **$3$ 维子空间**：$\mathrm{GF}(2)^3$ 本身，共 $1$ 个。

    这 $1 + 7 + 7 + 1 = 16$ 个子空间构成 $\mathrm{GF}(2)^3$ 的**子空间格**（subspace lattice）。

---

## 52.3 子空间计数与 Gauss 二项式系数

<div class="context-flow" markdown>

**核心问题**：$\mathrm{GF}(q)^n$ 中 $k$ 维子空间有多少个？这个计数公式与普通二项式系数有什么关系？

</div>

!!! definition "定义 52.5 (Gauss 二项式系数)"
    **Gauss 二项式系数**（或 $q$-二项式系数）定义为

    $$\binom{n}{k}_q = \frac{(q^n - 1)(q^n - q)(q^n - q^2) \cdots (q^n - q^{k-1})}{(q^k - 1)(q^k - q)(q^k - q^2) \cdots (q^k - q^{k-1})} = \prod_{i=0}^{k-1} \frac{q^{n-i} - 1}{q^{k-i} - 1},$$

    其中 $0 \leq k \leq n$，$q$ 为素数幂。约定 $\binom{n}{0}_q = 1$。

!!! theorem "定理 52.3 (Gauss 二项式系数 = 子空间个数)"
    $\mathrm{GF}(q)^n$ 中 $k$ 维子空间的个数恰为 $\binom{n}{k}_q$。

??? proof "证明"
    **计数方法**：$k$ 维子空间由 $k$ 个线性无关的向量张成。

    **分子**：选取 $k$ 个有序线性无关向量的方式数：
    - 第 $1$ 个向量：$q^n - 1$ 种选择（非零）。
    - 第 $2$ 个向量：$q^n - q$ 种选择（不在第 $1$ 个向量张成的 $1$ 维空间中，该空间有 $q$ 个元素）。
    - 第 $i$ 个向量：$q^n - q^{i-1}$ 种选择（不在前 $i-1$ 个向量张成的 $q^{i-1}$ 元空间中）。
    - 共 $(q^n-1)(q^n-q)\cdots(q^n-q^{k-1})$。

    **分母**：不同的有序 $k$ 元组可能张成同一个子空间。同一 $k$ 维子空间 $W$ 中选取 $k$ 个有序线性无关向量的方式数为 $(q^k-1)(q^k-q)\cdots(q^k-q^{k-1})$（即 $W \cong \mathrm{GF}(q)^k$ 中有序基的个数）。

    因此子空间个数 $= \frac{(q^n-1)(q^n-q)\cdots(q^n-q^{k-1})}{(q^k-1)(q^k-q)\cdots(q^k-q^{k-1})} = \binom{n}{k}_q$。

!!! theorem "定理 52.4 (Gauss 二项式系数的性质)"
    1. **$q \to 1$ 极限**：$\lim_{q \to 1} \binom{n}{k}_q = \binom{n}{k}$（普通二项式系数）。
    2. **对称性**：$\binom{n}{k}_q = \binom{n}{n-k}_q$。
    3. **$q$-Pascal 恒等式**：$\binom{n}{k}_q = \binom{n-1}{k-1}_q + q^k \binom{n-1}{k}_q$。
    4. **$q$-Vandermonde 恒等式**：$\binom{m+n}{k}_q = \sum_{j=0}^{k} q^{j(m-k+j)} \binom{m}{k-j}_q \binom{n}{j}_q$。
    5. $\binom{n}{k}_q$ 是 $q$ 的多项式（$q$-多项式），次数为 $k(n-k)$。

??? proof "证明（部分）"
    **(1)** $\frac{q^m - 1}{q - 1} = 1 + q + q^2 + \cdots + q^{m-1} \to m$（当 $q \to 1$）。故

    $$\binom{n}{k}_q = \prod_{i=0}^{k-1} \frac{q^{n-i}-1}{q^{k-i}-1} = \prod_{i=0}^{k-1} \frac{[n-i]_q}{[k-i]_q} \to \prod_{i=0}^{k-1} \frac{n-i}{k-i} = \frac{n!}{k!(n-k)!} = \binom{n}{k}.$$

    **(2)** 由 $k$ 维子空间和 $(n-k)$ 维子空间的自然双射（取正交补，或等价地取零化子）给出。

    **(3)** 组合证明：固定超平面 $H = \mathrm{GF}(q)^{n-1} \times \{0\}$。$\mathrm{GF}(q)^n$ 的 $k$ 维子空间 $W$ 分两类：

    - $W \subseteq H$：$H$ 中 $k$ 维子空间有 $\binom{n-1}{k}_q$ 个；
    - $W \not\subseteq H$：$W \cap H$ 是 $(k-1)$ 维子空间，$H$ 中 $(k-1)$ 维子空间有 $\binom{n-1}{k-1}_q$ 个。对每个这样的子空间，$W$ 的选取还需要在 $H$ 的余维（由一个 $\mathrm{GF}(q)^{n-1}$ 的补向量确定）中选择一个分量，给出 $q^{k-1}$ 个额外选择... 但精确分析给出系数 $q^k$。

!!! example "例 52.3 (Gauss 二项式系数的计算)"
    **(a)** $\mathrm{GF}(2)^4$ 中 $2$ 维子空间的个数：

    $$\binom{4}{2}_2 = \frac{(2^4-1)(2^4-2)}{(2^2-1)(2^2-2)} = \frac{15 \cdot 14}{3 \cdot 2} = \frac{210}{6} = 35.$$

    **(b)** $\mathrm{GF}(3)^3$ 中 $1$ 维子空间的个数：

    $$\binom{3}{1}_3 = \frac{3^3-1}{3-1} = \frac{26}{2} = 13.$$

    这 $13$ 条"直线"加上原点构成射影平面 $\mathrm{PG}(2,3)$，恰有 $13$ 个点（$= \frac{3^3-1}{3-1}$）和 $13$ 条线。

    **(c)** 验证 $q \to 1$：$\binom{4}{2}_q \to \binom{4}{2} = 6$。$\binom{4}{2}_q = \frac{(q^4-1)(q^3-1)}{(q^2-1)(q-1)} = \frac{[4]_q [3]_q}{[2]_q [1]_q}$。当 $q = 1$：$\frac{4 \cdot 3}{2 \cdot 1} = 6$。

---

## 52.4 一般线性群 $\mathrm{GL}(n,q)$

<div class="context-flow" markdown>

**核心问题**：$\mathrm{GF}(q)^n$ 的自同构群有多少元素？它的结构是什么？

</div>

!!! theorem "定理 52.5 ($\mathrm{GL}(n,q)$ 的阶)"
    一般线性群 $\mathrm{GL}(n,q) = \{A \in M_n(\mathrm{GF}(q)) : \det A \neq 0\}$ 的阶为

    $$|\mathrm{GL}(n,q)| = \prod_{i=0}^{n-1} (q^n - q^i) = (q^n-1)(q^n-q)(q^n-q^2)\cdots(q^n-q^{n-1}).$$

??? proof "证明"
    可逆矩阵的列构成 $\mathrm{GF}(q)^n$ 的基（有序基）。选取有序基的方式数即为答案：

    - 第 $1$ 列：$q^n - 1$（任意非零向量）。
    - 第 $2$ 列：$q^n - q$（不在第 $1$ 列张成的 $1$ 维空间中，该空间有 $q$ 个元素）。
    - 第 $i+1$ 列：$q^n - q^i$（不在前 $i$ 列张成的 $q^i$ 元子空间中）。

    共 $\prod_{i=0}^{n-1}(q^n - q^i)$。

!!! example "例 52.4 (低阶一般线性群)"
    **(a)** $|\mathrm{GL}(2,2)| = (4-1)(4-2) = 3 \cdot 2 = 6$。$\mathrm{GL}(2,2) \cong S_3$（对称群），因为它置换 $\mathrm{GF}(2)^2$ 的 $3$ 个非零向量。

    **(b)** $|\mathrm{GL}(2,3)| = (9-1)(9-3) = 8 \cdot 6 = 48$。

    **(c)** $|\mathrm{GL}(3,2)| = (8-1)(8-2)(8-4) = 7 \cdot 6 \cdot 4 = 168$。$\mathrm{GL}(3,2) = \mathrm{SL}(3,2)$（因为 $\mathrm{GF}(2)^{\times} = \{1\}$），这是最小的单群之一（同构于 $\mathrm{PSL}(2,7)$）。

!!! theorem "定理 52.6 (相关群的阶)"
    1. **特殊线性群**：$|\mathrm{SL}(n,q)| = \frac{|\mathrm{GL}(n,q)|}{q-1} = \frac{1}{q-1}\prod_{i=0}^{n-1}(q^n-q^i)$。
    2. **上三角可逆矩阵群（Borel 子群）**：$|B(n,q)| = (q-1)^n \cdot q^{n(n-1)/2}$。
    3. **上三角幂幺矩阵群**（对角线全为 $1$）：$|U(n,q)| = q^{n(n-1)/2}$。
    4. **射影一般线性群**：$|\mathrm{PGL}(n,q)| = \frac{|\mathrm{GL}(n,q)|}{q-1}$。

!!! definition "定义 52.6 ($\mathrm{GL}(n,q)$ 的 Bruhat 分解)"
    $\mathrm{GL}(n,q)$ 有 **Bruhat 分解**：

    $$\mathrm{GL}(n,q) = \bigsqcup_{w \in S_n} BwB,$$

    其中 $B$ 是上三角可逆矩阵群（Borel 子群），$w$ 视为置换矩阵，$S_n$ 是 $n$ 次对称群。这个分解在表示论中极为重要。

---

## 52.5 有限域上的标准形

<div class="context-flow" markdown>

**核心问题**：有限域上矩阵的标准形与实/复数域有何不同？

</div>

!!! theorem "定理 52.7 (有限域上的有理标准形)"
    **有理标准形**（rational canonical form）对任意域都成立，包括有限域。设 $A \in M_n(\mathrm{GF}(q))$。则 $A$ 相似于

    $$\begin{pmatrix} C(d_1) & & \\ & \ddots & \\ & & C(d_k) \end{pmatrix},$$

    其中 $d_1 \mid d_2 \mid \cdots \mid d_k$ 是 $\mathrm{GF}(q)[x]$ 中的首一多项式，$C(d_i)$ 是友矩阵。

!!! theorem "定理 52.8 (有限域上的 Jordan 形)"
    **Jordan 标准形在一般有限域上不存在**。Jordan 形要求特征多项式在基域上完全分裂（即分解为一次因子之积），这对 $\mathrm{GF}(q)$（$q$ 非无穷）一般不成立。

    **解决方案**：

    1. 在代数闭包 $\overline{\mathrm{GF}(q)}$ 上取 Jordan 形；
    2. 使用**有理标准形**（对任意域成立）；
    3. 对特征多项式在 $\mathrm{GF}(q)$ 上完全分裂的矩阵，Jordan 形在 $\mathrm{GF}(q)$ 上存在。

!!! example "例 52.5 (Jordan 形不存在的例子)"
    设 $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \in M_2(\mathrm{GF}(3))$。

    特征多项式 $\chi_A(x) = x^2 - 1 = (x-1)(x+1) = (x-1)(x-2)$（在 $\mathrm{GF}(3)$ 中 $-1 = 2$）。

    特征多项式完全分裂，故 Jordan 形在 $\mathrm{GF}(3)$ 上存在：$A$ 相似于 $\operatorname{diag}(1, 2)$。

    **另一个例子**：$B = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \in M_2(\mathrm{GF}(3))$。

    $\chi_B(x) = x^2 + 1$。在 $\mathrm{GF}(3)$ 中，$x^2 + 1$ 无根（$0^2+1=1$，$1^2+1=2$，$2^2+1=2$），故不可约。Jordan 形在 $\mathrm{GF}(3)$ 上不存在。有理标准形为 $B$ 本身（友矩阵）。

    在 $\mathrm{GF}(9) = \mathrm{GF}(3)(\alpha)$（$\alpha^2 = -1$）上，$x^2+1 = (x-\alpha)(x+\alpha)$，$B$ 可对角化为 $\operatorname{diag}(\alpha, -\alpha)$。

!!! definition "定义 52.7 (Frobenius 自同构对特征值的作用)"
    设 $A \in M_n(\mathrm{GF}(q))$，$\lambda \in \overline{\mathrm{GF}(q)}$ 是 $A$ 的特征值。则 $\lambda^q$ 也是 $A$ 的特征值（因为 Frobenius 自同构 $x \mapsto x^q$ 保持特征多项式的系数不变）。

    因此 $A$ 的特征值在 $\overline{\mathrm{GF}(q)}$ 中以 **Frobenius 轨道**出现：$\{\lambda, \lambda^q, \lambda^{q^2}, \ldots\}$。每个轨道对应 $\mathrm{GF}(q)$ 上的一个不可约因子。

---

## 52.6 线性码

<div class="context-flow" markdown>

**核心问题**：如何用线性代数构建纠错码？

</div>

!!! definition "定义 52.8 (线性码)"
    **$[n, k, d]_q$ 线性码** $C$ 是 $\mathrm{GF}(q)^n$ 的一个 $k$ 维子空间，其中 $d$ 是**最小距离**：

    $$d = d(C) = \min\{d_H(c_1, c_2) : c_1, c_2 \in C, \, c_1 \neq c_2\} = \min\{w_H(c) : c \in C, \, c \neq 0\},$$

    这里 $d_H(x, y) = |\{i : x_i \neq y_i\}|$ 是**Hamming 距离**，$w_H(c) = |\{i : c_i \neq 0\}|$ 是**Hamming 重量**。

    参数含义：$n$ = 码长，$k$ = 信息位数（编码 $q^k$ 个消息），$d$ = 最小距离（可纠正 $\lfloor(d-1)/2\rfloor$ 个错误）。

    **码率** $R = k/n$（信息传输效率），**相对距离** $\delta = d/n$（纠错能力）。

!!! definition "定义 52.9 (生成矩阵与校验矩阵)"
    设 $C$ 是 $[n, k]_q$ 线性码。

    - **生成矩阵** $G$ 是 $k \times n$ 矩阵，其行构成 $C$ 的一组基。$C = \{xG : x \in \mathrm{GF}(q)^k\}$。
    - **校验矩阵**（奇偶校验矩阵）$H$ 是 $(n-k) \times n$ 矩阵，满足 $C = \ker H = \{c \in \mathrm{GF}(q)^n : Hc^T = 0\}$。

    等价地，$G$ 的行空间 $= C = H$ 的零空间，$GH^T = 0$。

!!! theorem "定理 52.9 (最小距离与校验矩阵)"
    设 $H$ 是 $[n,k,d]_q$ 码 $C$ 的校验矩阵。则

    $$d = \min\{j : H \text{ 的某 } j \text{ 列线性相关}\}.$$

    等价地，$d \geq t + 1$ 当且仅当 $H$ 的任意 $t$ 列线性无关。

??? proof "证明"
    $c \in C$ 当且仅当 $Hc^T = 0$，即 $c$ 的各分量是 $H$ 的列的一个线性组合等于零。$w_H(c) = j$ 意味着这个线性组合涉及 $j$ 列。因此 $d = \min\{w_H(c) : c \neq 0, Hc^T = 0\} = \min\{j : H$ 的某 $j$ 列线性相关$\}$。

!!! definition "定义 52.10 (系统形式)"
    通过行变换，$G$ 可化为**系统形式**（标准形）：

    $$G = [I_k \mid P],$$

    其中 $I_k$ 是 $k \times k$ 单位矩阵，$P$ 是 $k \times (n-k)$ 矩阵。此时校验矩阵为

    $$H = [-P^T \mid I_{n-k}].$$

    验证：$GH^T = I_k(-P) + P \cdot I_{n-k} = -P + P = 0$（在 $\mathrm{GF}(2)$ 上 $-P = P$，简化为 $GH^T = 0$）。

!!! definition "定义 52.11 (伴随式译码)"
    接收到 $y = c + e$（$c$ 是发送的码字，$e$ 是错误向量）。**伴随式**（syndrome）为

    $$s = Hy^T = H(c + e)^T = Hc^T + He^T = He^T.$$

    伴随式只取决于错误 $e$，不取决于发送的码字。若 $s = 0$，则 $y \in C$，假设无错误（或检测不到的错误）。若 $s \neq 0$，在 $s$ 对应的陪集中找最小重量的向量（**陪集首**），作为 $\hat{e}$，译码为 $\hat{c} = y - \hat{e}$。

!!! example "例 52.6 ($[7, 4, 3]_2$ Hamming 码)"
    校验矩阵（列是 $1$ 到 $7$ 的二进制表示）：

    $$H = \begin{pmatrix} 1 & 0 & 1 & 0 & 1 & 0 & 1 \\ 0 & 1 & 1 & 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 & 1 & 1 & 1 \end{pmatrix}.$$

    $H$ 是 $3 \times 7$ 矩阵，$n - k = 3$，故 $k = 4$。任意两列线性无关（任意两个不同的 $3$ 位二进制数不成比例），但有三列线性相关（例如第 $1, 2, 3$ 列：$(1,0,0)^T + (0,1,0)^T + (1,1,0)^T = (0,0,0)^T$），故 $d = 3$。

    可纠正 $\lfloor(3-1)/2\rfloor = 1$ 个错误。伴随式 $s = He^T$ 等于出错位置的二进制表示。

---

## 52.7 重要线性码族

<div class="context-flow" markdown>

**核心问题**：有哪些重要的线性码构造？它们如何用有限域线性代数描述？

</div>

!!! definition "定义 52.12 (Hamming 码)"
    **$[2^r-1, 2^r-1-r, 3]_2$ Hamming 码**：校验矩阵 $H$ 是 $r \times (2^r-1)$ 矩阵，其列是 $\mathrm{GF}(2)^r$ 的所有非零向量。

    - 码长 $n = 2^r - 1$；
    - 维数 $k = 2^r - 1 - r$；
    - 最小距离 $d = 3$（任意两列不同，故线性无关；但存在三列之和为零）；
    - 可纠正 $1$ 位错误（完美码）。

    推广到 $q$-ary：$[(q^r-1)/(q-1), (q^r-1)/(q-1)-r, 3]_q$ Hamming 码，校验矩阵的列是 $\mathrm{GF}(q)^r$ 中所有 $1$ 维子空间的代表。

!!! definition "定义 52.13 (对偶码)"
    $[n, k]_q$ 线性码 $C$ 的**对偶码**为

    $$C^{\perp} = \{x \in \mathrm{GF}(q)^n : \langle x, c \rangle = 0, \, \forall c \in C\},$$

    其中 $\langle x, c \rangle = \sum_i x_i c_i$ 是标准内积。$C^{\perp}$ 是 $[n, n-k]_q$ 线性码。

    若 $G$ 是 $C$ 的生成矩阵，则 $G$ 也是 $C^{\perp}$ 的校验矩阵；$H$ 也是 $C^{\perp}$ 的生成矩阵。

    **自对偶码**：$C = C^{\perp}$（需要 $k = n/2$）。

!!! theorem "定理 52.10 (MacWilliams 恒等式)"
    设 $C$ 是 $[n, k]_q$ 线性码。定义 $C$ 的**重量枚举多项式**为

    $$W_C(x, y) = \sum_{c \in C} x^{n - w_H(c)} y^{w_H(c)} = \sum_{i=0}^{n} A_i x^{n-i} y^i,$$

    其中 $A_i = |\{c \in C : w_H(c) = i\}|$。则 $C^{\perp}$ 的重量枚举多项式为

    $$W_{C^{\perp}}(x, y) = \frac{1}{|C|} W_C(x + (q-1)y, \, x - y).$$

    MacWilliams 恒等式将码 $C$ 的重量分布与对偶码 $C^{\perp}$ 的重量分布联系起来。

!!! definition "定义 52.14 (Reed-Muller 码)"
    **Reed-Muller 码** $\mathrm{RM}(r, m)$ 定义在 $\mathrm{GF}(2)$ 上：

    - 码长 $n = 2^m$；
    - 维数 $k = \sum_{i=0}^{r} \binom{m}{i}$；
    - 最小距离 $d = 2^{m-r}$。

    $\mathrm{RM}(r,m)$ 由 $m$ 元 $\mathrm{GF}(2)$ 多项式（次数 $\leq r$）在 $\mathrm{GF}(2)^m$ 的所有点上的求值构成。

    特别地：
    - $\mathrm{RM}(0, m) = \{0\ldots0, 1\ldots1\}$（重复码）；
    - $\mathrm{RM}(1, m)$ 包含所有仿射函数的求值（一阶 RM 码，用于 Mariner 任务的图像传输）；
    - $\mathrm{RM}(m, m) = \mathrm{GF}(2)^{2^m}$（全空间）；
    - $\mathrm{RM}(m-1, m) = \mathrm{RM}(0,m)^{\perp}$（偶重量码）。

!!! example "例 52.7 (一阶 Reed-Muller 码 $\mathrm{RM}(1, 3)$)"
    $m = 3$，$r = 1$。码长 $n = 8$，$k = 1 + 3 = 4$，$d = 2^2 = 4$。这是 $[8, 4, 4]_2$ 码。

    $\mathrm{GF}(2)^3$ 的 $8$ 个点（按字典序）：$(0,0,0), (0,0,1), \ldots, (1,1,1)$。

    基码字（对应函数 $1, x_1, x_2, x_3$）：

    - $f = 1$：$(1,1,1,1,1,1,1,1)$
    - $f = x_1$：$(0,0,0,0,1,1,1,1)$
    - $f = x_2$：$(0,0,1,1,0,0,1,1)$
    - $f = x_3$：$(0,1,0,1,0,1,0,1)$

---

## 52.8 应用

<div class="context-flow" markdown>

**核心问题**：有限域线性代数在密码学和通信中有哪些具体应用？

</div>

!!! example "例 52.8 (AES 中的 $\mathrm{GF}(2^8)$ 运算)"
    **高级加密标准**（AES）的核心运算在 $\mathrm{GF}(2^8)$ 上进行，使用不可约多项式 $f(x) = x^8 + x^4 + x^3 + x + 1$。

    **S-box（替代盒）**：将输入字节 $a \in \mathrm{GF}(2^8)$ 映射为 $a^{-1}$（在 $\mathrm{GF}(2^8)$ 中求逆，$0 \mapsto 0$），然后进行仿射变换。

    **MixColumns**：将 $4$ 个字节视为 $\mathrm{GF}(2^8)[x]/(x^4+1)$ 中的元素，左乘固定矩阵

    $$\begin{pmatrix} 02 & 03 & 01 & 01 \\ 01 & 02 & 03 & 01 \\ 01 & 01 & 02 & 03 \\ 03 & 01 & 01 & 02 \end{pmatrix}$$

    （元素为 $\mathrm{GF}(2^8)$ 中的元素，十六进制表示）。这个矩阵在 $\mathrm{GF}(2^8)$ 上可逆，保证了解密的可行性。

!!! example "例 52.9 (Reed-Solomon 码在 QR 码中的应用)"
    **QR 码**使用 $\mathrm{GF}(2^8)$ 上的 Reed-Solomon 码。设 $\alpha$ 是 $\mathrm{GF}(2^8)$ 的本原元（$\alpha^{255} = 1$）。

    **$[n, k, d]_{256}$ Reed-Solomon 码**的生成多项式为

    $$g(x) = (x - \alpha)(x - \alpha^2) \cdots (x - \alpha^{d-1}).$$

    码字是 $\mathrm{GF}(2^8)^n$ 中满足 $g(x) \mid c(x)$ 的向量。最小距离恰为 $d$，达到了 **Singleton 界** $d \leq n - k + 1$，因此是 **MDS 码**（最大距离可分码）。

    **RS 码的应用**：
    - CD/DVD 使用双层交错 RS 码（CIRC）；
    - QR 码使用 RS 码实现最高 30% 的容错率；
    - 深空通信（旅行者号探测器）；
    - 数字电视（DVB 标准）。

!!! example "例 52.10 (LDPC 码)"
    **低密度奇偶校验码**（LDPC, Low-Density Parity-Check Code）由 Gallager 于 1960 年提出。其校验矩阵 $H$ 是**稀疏矩阵**——大部分元素为零。

    $\mathrm{GF}(2)$ 上的 LDPC 码 $C$ 定义为 $C = \ker H$，其中 $H$ 的每行恰有 $d_c$ 个 $1$（行重均匀），每列恰有 $d_v$ 个 $1$（列重均匀）。这样的码称为 **$(d_v, d_c)$-正则 LDPC 码**。

    LDPC 码通过迭代译码（信念传播算法）可以逼近 **Shannon 限**，在现代通信标准（5G NR、Wi-Fi 6）中广泛使用。

!!! theorem "定理 52.11 (Singleton 界)"
    任何 $[n, k, d]_q$ 线性码满足

    $$d \leq n - k + 1.$$

    达到此界的码称为 **MDS 码**（maximum distance separable code）。Reed-Solomon 码是 MDS 码的典型例子。

??? proof "证明"
    删除码字的任意 $d-1$ 个坐标，得到的缩短码仍然是单射的（因为任意两个不同码字在至少 $d$ 个位置不同，删 $d-1$ 个位置后仍然不同）。缩短后的码长为 $n - (d-1)$，维数仍为 $k$。由于码字数 $q^k \leq q^{n-d+1}$，得 $k \leq n - d + 1$，即 $d \leq n - k + 1$。

!!! theorem "定理 52.12 (Hamming 界)"
    任何 $[n, k, d]_q$ 码（$d = 2t + 1$）满足

    $$q^k \sum_{i=0}^{t} \binom{n}{i}(q-1)^i \leq q^n.$$

    达到此界的码称为**完美码**。Hamming 码和 Golay 码是完美码。

!!! example "例 52.11 (有限域线性代数的更多应用)"
    | 应用领域 | 有限域 | 线性代数工具 |
    |:---|:---|:---|
    | AES 加密 | $\mathrm{GF}(2^8)$ | 矩阵乘法、逆矩阵 |
    | RS 纠错码 | $\mathrm{GF}(2^8)$ | 多项式、Vandermonde 矩阵 |
    | QR 码 | $\mathrm{GF}(2^8)$ | RS 码 |
    | 椭圆曲线密码 | $\mathrm{GF}(p)$, $\mathrm{GF}(2^n)$ | 椭圆曲线上的线性代数 |
    | LDPC 码 | $\mathrm{GF}(2)$ | 稀疏矩阵、零空间 |
    | 秘密共享 | $\mathrm{GF}(p)$ | Lagrange 插值 |
    | 网络编码 | $\mathrm{GF}(q)$ | 矩阵秩、子空间 |

---

### 本章总结

有限域上的线性代数保留了向量空间理论的代数框架，但引入了丰富的组合结构：

- **子空间的个数**由 Gauss 二项式系数 $\binom{n}{k}_q$ 给出，它是普通二项式系数的 $q$-模拟；
- **$\mathrm{GL}(n,q)$ 的阶**有精确的闭合公式，反映了有限域上"非退化"条件的组合意义；
- **标准形理论**中，Jordan 形需要代数闭域，但有理标准形对任意域成立；
- **线性码**将子空间的代数结构转化为纠错能力，通过生成矩阵和校验矩阵的对偶关系实现编码和译码。

这些工具在信息时代的应用——从手机通信到量子计算——使得有限域线性代数成为最具实际影响力的数学分支之一。

---

### 习题

!!! exercise "习题 52.1"
    计算 $\mathrm{GF}(2)^5$ 中 $2$ 维子空间的个数 $\binom{5}{2}_2$。

!!! exercise "习题 52.2"
    构造 $\mathrm{GF}(2^3)$，使用不可约多项式 $x^3 + x + 1$。列出所有 $8$ 个元素及其乘法表。

!!! exercise "习题 52.3"
    计算 $|\mathrm{GL}(4, 2)|$ 并分解为素因子。

!!! exercise "习题 52.4"
    设 $C$ 是 $[7, 4, 3]_2$ Hamming 码。若接收到 $y = (1,1,0,0,1,0,1)$，计算伴随式 $s = Hy^T$ 并确定错误位置。

!!! exercise "习题 52.5"
    证明 Hamming 码 $[2^r-1, 2^r-r-1, 3]_2$ 是完美码（即达到 Hamming 界）。

!!! exercise "习题 52.6"
    构造一个 $[6, 3, 3]_2$ 线性码，给出其生成矩阵和校验矩阵。这个码是否是 MDS 码？
