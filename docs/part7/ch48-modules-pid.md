# 第 48 章 主理想整环上的模

<div class="context-flow" markdown>

**前置**：向量空间 (Ch4) · 多项式代数 (Ch0) · Jordan 形 (Ch12) · 有理标准形 (Ch13B)

**本章脉络**：模的定义 → 自由模 → 子模 → 商模 → 挠模 → PID 上有限生成模的结构定理 → 不变因子与初等因子 → 应用到 $\mathbb{F}[x]$-模 → 重新推导 Jordan 形与有理标准形

**延伸**：模论是交换代数和代数几何的基础语言；有限生成 Abel 群基本定理是 $\mathbb{Z}$-模的特例；层论中的 sheaf of modules 将模论推广到几何对象上

</div>

向量空间理论的核心在于域 $\mathbb{F}$ 的每个非零元素都可逆——正是这一性质保证了维数的良定义性和基的存在性。当我们将标量环从域放宽到一般的交换环 $R$ 时，便进入了**模论**（module theory）的范畴。模论不仅是向量空间理论的自然推广，更为我们提供了一个统一的框架：**有限维向量空间的线性变换的分类问题**（Jordan 标准形、有理标准形）与**有限生成 Abel 群的结构定理**，在模论的视角下不过是同一个定理的不同特例。

本章将聚焦于**主理想整环**（principal ideal domain, PID）上的有限生成模。PID 是一类特别良好的环——包括整数环 $\mathbb{Z}$ 和域上的一元多项式环 $\mathbb{F}[x]$——在其上可以建立类似于向量空间的结构理论。最终我们将看到，著名的 Jordan 标准形和有理标准形，本质上是 $\mathbb{F}[x]$-模的结构定理的直接推论。

---

## 48.1 模的基本概念

<div class="context-flow" markdown>

**核心问题**：如何将向量空间的定义从域推广到一般环？推广后会失去哪些性质？

</div>

!!! definition "定义 48.1 (左 $R$-模)"
    设 $R$ 是一个含幺环（不一定交换）。一个**左 $R$-模**是一个 Abel 群 $(M, +)$ 连同一个标量乘法映射 $R \times M \to M$，记为 $(r, m) \mapsto rm$，满足以下公理：对所有 $r, s \in R$ 和 $m, n \in M$，

    1. **分配律 I**：$r(m + n) = rm + rn$；
    2. **分配律 II**：$(r + s)m = rm + sm$；
    3. **结合律**：$(rs)m = r(sm)$；
    4. **幺元作用**：$1_R \cdot m = m$。

!!! definition "定义 48.2 (右 $R$-模)"
    类似地，**右 $R$-模**是一个 Abel 群 $(M, +)$ 连同映射 $M \times R \to M$，记为 $(m, r) \mapsto mr$，满足相应的公理。当 $R$ 是交换环时，左 $R$-模和右 $R$-模没有本质区别。本章主要讨论交换环上的模，因此简称为 **$R$-模**。

!!! example "例 48.1 (模的基本例子)"
    **(a) 向量空间作为 $\mathbb{F}$-模。** 设 $\mathbb{F}$ 是一个域，$V$ 是 $\mathbb{F}$ 上的向量空间。那么 $V$ 自然就是一个 $\mathbb{F}$-模。事实上，$\mathbb{F}$-模与 $\mathbb{F}$-向量空间完全等价。

    **(b) Abel 群作为 $\mathbb{Z}$-模。** 任何 Abel 群 $(G, +)$ 都自然地成为 $\mathbb{Z}$-模：对 $n \in \mathbb{Z}$，$g \in G$，定义

    $$ng = \underbrace{g + g + \cdots + g}_{n \text{ 个}}, \quad n > 0; \qquad 0 \cdot g = 0; \qquad ng = -((-n)g), \quad n < 0.$$

    反之，每个 $\mathbb{Z}$-模的底层 Abel 群就是它本身。因此 **$\mathbb{Z}$-模与 Abel 群等价**。

    **(c) 线性变换诱导的 $\mathbb{F}[x]$-模。** 设 $V$ 是有限维 $\mathbb{F}$-向量空间，$T: V \to V$ 是线性变换。定义 $\mathbb{F}[x]$ 在 $V$ 上的作用：对 $f(x) = a_0 + a_1 x + \cdots + a_n x^n \in \mathbb{F}[x]$ 和 $v \in V$，

    $$f(x) \cdot v = a_0 v + a_1 T(v) + a_2 T^2(v) + \cdots + a_n T^n(v) = f(T)(v).$$

    容易验证这使得 $V$ 成为 $\mathbb{F}[x]$-模。这个构造是本章的关键——**线性变换 $T$ 的全部信息都编码在 $\mathbb{F}[x]$-模结构中**。

!!! example "例 48.2 (理想作为模)"
    设 $R$ 是交换环，$I \subseteq R$ 是一个理想。那么 $I$ 自然是一个 $R$-模（作为 $R$ 的子模）。特别地，$R$ 本身是一个 $R$-模。

!!! definition "定义 48.3 ($R$-模同态)"
    设 $M, N$ 是 $R$-模。映射 $\varphi: M \to N$ 称为 **$R$-模同态**（或 $R$-线性映射），若对所有 $r \in R$，$m, m' \in M$：

    1. $\varphi(m + m') = \varphi(m) + \varphi(m')$；
    2. $\varphi(rm) = r\varphi(m)$。

    所有从 $M$ 到 $N$ 的 $R$-模同态构成的集合记为 $\operatorname{Hom}_R(M, N)$，它自然是一个 $R$-模。若 $\varphi$ 是双射，则称为**同构**，记 $M \cong N$。

!!! theorem "定理 48.1 (向量空间 vs. 模的关键区别)"
    设 $R$ 是交换环（不一定是域）。$R$-模与 $R$-向量空间（当 $R$ 是域时）的关键区别在于：

    1. **不是每个模都有基**（即不是每个模都是自由模）；
    2. **$rm = 0$ 不一定意味着 $r = 0$ 或 $m = 0$**（挠元素的存在）；
    3. **子模不一定有补**（即不一定存在直和分解）；
    4. **模同态的像不一定是直和项**。

??? proof "说明"
    **(1)** $\mathbb{Z}$-模 $\mathbb{Z}/2\mathbb{Z}$ 没有基：任何元素 $\bar{a}$ 满足 $2\bar{a} = \bar{0}$，故 $\{\bar{a}\}$ 不构成基（系数不唯一：$0 \cdot \bar{a} = 2 \cdot \bar{a} = \bar{0}$）。

    **(2)** 在 $\mathbb{Z}/6\mathbb{Z}$ 中，$2 \cdot \bar{3} = \bar{0}$，但 $2 \neq 0$，$\bar{3} \neq \bar{0}$。

    **(3)** $\mathbb{Z}$-模 $\mathbb{Z}$ 中，子模 $2\mathbb{Z}$ 没有补模 $N$ 使得 $\mathbb{Z} = 2\mathbb{Z} \oplus N$（因为 $\mathbb{Z}/2\mathbb{Z}$ 不可嵌入 $\mathbb{Z}$）。

    **(4)** 包含映射 $2\mathbb{Z} \hookrightarrow \mathbb{Z}$ 的像 $2\mathbb{Z}$ 不是 $\mathbb{Z}$ 的直和项。

---

## 48.2 子模与商模

<div class="context-flow" markdown>

**核心问题**：如何定义模的子对象和商对象？同构定理如何推广？

</div>

!!! definition "定义 48.4 (子模)"
    设 $M$ 是 $R$-模。子集 $N \subseteq M$ 称为**子模**，若 $N$ 在加法和标量乘法下封闭：

    1. $0 \in N$；
    2. 若 $m, n \in N$，则 $m + n \in N$；
    3. 若 $r \in R$，$n \in N$，则 $rn \in N$。

!!! definition "定义 48.5 (生成集与有限生成模)"
    设 $M$ 是 $R$-模，$S \subseteq M$。$S$ **生成**的子模定义为

    $$\langle S \rangle = \left\{ \sum_{i=1}^{k} r_i s_i : k \in \mathbb{N},\, r_i \in R,\, s_i \in S \right\}.$$

    若存在有限集 $S = \{m_1, \ldots, m_n\}$ 使得 $M = \langle S \rangle$，则称 $M$ 是**有限生成**的，记 $M = Rm_1 + Rm_2 + \cdots + Rm_n$。

!!! definition "定义 48.6 (商模)"
    设 $N$ 是 $R$-模 $M$ 的子模。**商模** $M/N$ 定义为 Abel 群 $M/N$（按 $N$ 的陪集分类），其上的 $R$-作用为

    $$r(m + N) = rm + N, \quad r \in R,\, m \in M.$$

    自然投射 $\pi: M \to M/N$，$\pi(m) = m + N$ 是满的 $R$-模同态，其核为 $N$。

!!! theorem "定理 48.2 (模的同构定理)"
    **(第一同构定理)** 设 $\varphi: M \to N$ 是 $R$-模同态。则

    $$M / \ker \varphi \cong \operatorname{im} \varphi.$$

    **(第二同构定理)** 设 $A, B$ 是 $M$ 的子模。则

    $$(A + B) / B \cong A / (A \cap B).$$

    **(第三同构定理)** 设 $N \subseteq L$ 是 $M$ 的子模。则

    $$(M/N) / (L/N) \cong M/L.$$

??? proof "证明"
    **第一同构定理：** 定义 $\bar{\varphi}: M/\ker\varphi \to \operatorname{im}\varphi$，$\bar{\varphi}(m + \ker\varphi) = \varphi(m)$。

    **良定义性：** 若 $m + \ker\varphi = m' + \ker\varphi$，则 $m - m' \in \ker\varphi$，故 $\varphi(m) = \varphi(m')$。

    **$R$-线性：** $\bar{\varphi}(r(m+\ker\varphi) + (m'+\ker\varphi)) = \bar{\varphi}((rm+m')+\ker\varphi) = \varphi(rm+m') = r\varphi(m)+\varphi(m') = r\bar{\varphi}(m+\ker\varphi) + \bar{\varphi}(m'+\ker\varphi)$。

    **单射：** 若 $\bar{\varphi}(m+\ker\varphi) = 0$，则 $\varphi(m) = 0$，即 $m \in \ker\varphi$，故 $m + \ker\varphi = \ker\varphi$。

    **满射：** 由 $\operatorname{im}\varphi$ 的定义显然。

    第二和第三同构定理的证明与群的情形完全类似，此处从略。

!!! definition "定义 48.7 (直和)"
    设 $M_1, \ldots, M_k$ 是 $R$-模。它们的**外直和**定义为

    $$M_1 \oplus M_2 \oplus \cdots \oplus M_k = \{(m_1, m_2, \ldots, m_k) : m_i \in M_i\},$$

    其上的 $R$-作用分量逐一进行。

    若 $M$ 的子模 $N_1, \ldots, N_k$ 满足：(i) $M = N_1 + N_2 + \cdots + N_k$，且 (ii) 对每个 $i$，$N_i \cap (N_1 + \cdots + N_{i-1} + N_{i+1} + \cdots + N_k) = \{0\}$，则称 $M$ 是 $N_1, \ldots, N_k$ 的**内直和**，记 $M = N_1 \oplus N_2 \oplus \cdots \oplus N_k$。

!!! example "例 48.3 (商模的计算)"
    **(a)** $\mathbb{Z}/n\mathbb{Z}$ 是 $\mathbb{Z}$-模 $\mathbb{Z}$ 关于子模 $n\mathbb{Z}$ 的商模。

    **(b)** 设 $R = \mathbb{F}[x]$，$I = (f(x))$ 是 $f(x)$ 生成的理想（$f$ 次数为 $d$）。则 $R/I = \mathbb{F}[x]/(f(x))$ 作为 $\mathbb{F}$-向量空间有基 $\{1, \bar{x}, \bar{x}^2, \ldots, \bar{x}^{d-1}\}$，维数为 $d = \deg f$。

    **(c)** 设 $V$ 是 $n$ 维 $\mathbb{F}$-向量空间，$T: V \to V$ 是线性变换，$V$ 按例 48.1(c) 成为 $\mathbb{F}[x]$-模。$V$ 的 $T$-不变子空间恰好是 $V$ 的 $\mathbb{F}[x]$-子模。

---

## 48.3 自由模与秩

<div class="context-flow" markdown>

**核心问题**：哪些模最像向量空间？PID 上的自由模有什么特殊性质？

</div>

!!! definition "定义 48.8 (自由模)"
    $R$-模 $F$ 称为**自由模**，若 $F$ 有一个**基**（basis），即存在子集 $B \subseteq F$ 使得：

    1. $B$ 生成 $F$：$F = \langle B \rangle$；
    2. $B$ 线性无关：若 $\sum_{i=1}^{k} r_i b_i = 0$（$r_i \in R$，$b_i \in B$ 两两不同），则 $r_1 = r_2 = \cdots = r_k = 0$。

    基为 $n$ 元集的自由模同构于 $R^n$，其中 $n$ 称为自由模的**秩**（rank）。

!!! theorem "定理 48.3 (秩的良定义性)"
    设 $R$ 是交换环（具有不变基数性质，IBN）。若 $R^m \cong R^n$，则 $m = n$。特别地，交换环上自由模的秩是良定义的。

??? proof "证明"
    设 $\mathfrak{m}$ 是 $R$ 的极大理想（由 Zorn 引理保证存在）。取商得到域 $k = R/\mathfrak{m}$。

    若 $R^m \cong R^n$，对两边模掉 $\mathfrak{m}$ 得

    $$(R/\mathfrak{m})^m \cong (R/\mathfrak{m})^n, \quad \text{即} \quad k^m \cong k^n.$$

    由向量空间维数的唯一性，$m = n$。

!!! theorem "定理 48.4 (PID 上自由模的子模是自由的)"
    设 $R$ 是主理想整环，$F$ 是秩为 $n$ 的自由 $R$-模。则 $F$ 的每个子模 $N$ 也是自由的，且 $\operatorname{rank}(N) \leq n$。

??? proof "证明"
    对 $n$ 进行归纳。

    **$n = 1$：** $F = R$，子模 $N$ 是 $R$ 的理想。由于 $R$ 是 PID，$N = (d)$ 对某个 $d \in R$。若 $d = 0$，则 $N = \{0\}$ 是秩 $0$ 的自由模。若 $d \neq 0$，则 $N \cong R$（通过映射 $r \mapsto rd$，注意 $R$ 是整环故此映射为单射），是秩 $1$ 的自由模。

    **归纳步骤：** 设结论对秩 $< n$ 的自由模成立。令 $F = R^n$，$\pi: F \to R$ 为投射到最后一个分量的映射。$\pi(N)$ 是 $R$ 的理想，由 PID 性质，$\pi(N) = (d)$。

    若 $d = 0$，则 $N \subseteq \ker \pi = R^{n-1}$，由归纳假设 $N$ 自由且秩 $\leq n - 1$。

    若 $d \neq 0$，选取 $e \in N$ 使得 $\pi(e) = d$。对任意 $m \in N$，$\pi(m) = rd$ 对某 $r \in R$，故 $m - re \in N \cap \ker\pi$。令 $N' = N \cap \ker\pi \subseteq R^{n-1}$，由归纳假设 $N'$ 自由。

    断言 $N = N' \oplus Re$：$N = N' + Re$ 由上述论证可得；若 $n' + re = 0$（$n' \in N'$），则 $\pi(n' + re) = rd = 0$，由 $R$ 是整环和 $d \neq 0$ 得 $r = 0$，进而 $n' = 0$。故 $N$ 自由，$\operatorname{rank}(N) = \operatorname{rank}(N') + 1 \leq (n-1) + 1 = n$。

!!! example "例 48.4 (自由模与非自由模)"
    **(a)** $\mathbb{Z}^n$ 是秩为 $n$ 的自由 $\mathbb{Z}$-模。$\mathbb{Z}^n$ 的每个子群（子模）都是自由 Abel 群，秩不超过 $n$。

    **(b)** $\mathbb{Z}/n\mathbb{Z}$（$n \geq 2$）不是自由 $\mathbb{Z}$-模——它没有基。

    **(c)** $\mathbb{F}[x]^n$ 是秩为 $n$ 的自由 $\mathbb{F}[x]$-模。

    **(d)** 设 $R = \mathbb{Z}[x, y]$（不是 PID），理想 $I = (x, y)$ 是 $R$ 的子模，但 $I$ 不是自由 $R$-模。这表明定理 48.4 中 PID 的条件不可省略。

---

## 48.4 挠模与挠元素

<div class="context-flow" markdown>

**核心问题**：模中什么元素的行为与向量空间不同？如何刻画这些"病态"元素？

</div>

!!! definition "定义 48.9 (挠元素与挠子模)"
    设 $R$ 是整环，$M$ 是 $R$-模。元素 $m \in M$ 称为**挠元素**（torsion element），若存在非零 $r \in R$ 使得 $rm = 0$。$m$ 的**零化子**定义为

    $$\operatorname{Ann}(m) = \{r \in R : rm = 0\},$$

    这是 $R$ 的一个理想。$M$ 的**挠子模**定义为

    $$\operatorname{Tor}(M) = \{m \in M : \exists\, r \in R \setminus \{0\},\, rm = 0\}.$$

    若 $\operatorname{Tor}(M) = M$，则称 $M$ 是**挠模**（torsion module）。若 $\operatorname{Tor}(M) = \{0\}$，则称 $M$ 是**无挠模**（torsion-free module）。

!!! theorem "定理 48.5 (挠子模是子模)"
    设 $R$ 是整环，$M$ 是 $R$-模。则 $\operatorname{Tor}(M)$ 是 $M$ 的子模，且 $M / \operatorname{Tor}(M)$ 是无挠的。

??? proof "证明"
    设 $m, n \in \operatorname{Tor}(M)$，$r, s \in R \setminus \{0\}$ 使得 $rm = 0$，$sn = 0$。则 $rs \neq 0$（$R$ 是整环），且

    $$rs(m + n) = s(rm) + r(sn) = 0,$$

    故 $m + n \in \operatorname{Tor}(M)$。对任意 $a \in R$，$r(am) = a(rm) = 0$，故 $am \in \operatorname{Tor}(M)$。因此 $\operatorname{Tor}(M)$ 是子模。

    设 $\bar{m} = m + \operatorname{Tor}(M) \in M/\operatorname{Tor}(M)$ 是挠元素，即存在 $r \neq 0$ 使得 $r\bar{m} = \bar{0}$，即 $rm \in \operatorname{Tor}(M)$。则存在 $s \neq 0$ 使得 $s(rm) = (sr)m = 0$。由于 $sr \neq 0$，$m \in \operatorname{Tor}(M)$，即 $\bar{m} = \bar{0}$。

!!! example "例 48.5 (挠模的例子)"
    **(a)** $\mathbb{Z}/n\mathbb{Z}$ 是挠 $\mathbb{Z}$-模：每个元素 $\bar{a}$ 满足 $n\bar{a} = \bar{0}$。

    **(b)** $\mathbb{Z}$ 本身是无挠 $\mathbb{Z}$-模。

    **(c)** $\mathbb{Z} \oplus \mathbb{Z}/6\mathbb{Z}$ 既不是挠模也不是无挠模。其挠子模为 $\{0\} \oplus \mathbb{Z}/6\mathbb{Z}$。

    **(d)** 设 $V$ 是有限维 $\mathbb{F}$-向量空间，$T \in \operatorname{End}(V)$，视 $V$ 为 $\mathbb{F}[x]$-模。由于 $V$ 有限维，对每个 $v \in V$，向量 $v, Tv, T^2v, \ldots$ 必线性相关，故存在非零多项式 $f(x) \in \mathbb{F}[x]$ 使得 $f(T)v = 0$。因此 **$V$ 作为 $\mathbb{F}[x]$-模是挠模**。

!!! definition "定义 48.10 (零化子理想)"
    设 $M$ 是 $R$-模。$M$ 的**零化子**定义为

    $$\operatorname{Ann}(M) = \{r \in R : rm = 0, \forall m \in M\}.$$

    这是 $R$ 的一个理想。若 $R$ 是 PID 且 $M$ 是有限生成挠模，则 $\operatorname{Ann}(M) = (d)$ 对某个 $d \in R$。

---

## 48.5 PID 上有限生成模的结构定理

<div class="context-flow" markdown>

**核心问题**：PID 上的有限生成模有怎样的分类？这是本章的中心定理。

</div>

!!! theorem "定理 48.6 (PID 上有限生成模的结构定理——不变因子形式)"
    设 $R$ 是主理想整环，$M$ 是有限生成 $R$-模。则存在唯一的非负整数 $r \geq 0$ 和非零非单位元素 $d_1, d_2, \ldots, d_k \in R$（在相伴意义下唯一），满足

    $$d_1 \mid d_2 \mid \cdots \mid d_k,$$

    使得

    $$M \cong R^r \oplus R/(d_1) \oplus R/(d_2) \oplus \cdots \oplus R/(d_k).$$

    其中 $r$ 称为 $M$ 的**自由秩**（free rank），$d_1, \ldots, d_k$ 称为 $M$ 的**不变因子**（invariant factors）。

该定理的证明需要几个步骤。我们通过 **Smith 标准形**给出一个构造性的证明。

!!! definition "定义 48.11 (Smith 标准形)"
    设 $R$ 是 PID，$A$ 是 $m \times n$ 的 $R$-矩阵。$A$ 的 **Smith 标准形**是一个对角矩阵

    $$S = \begin{pmatrix} d_1 & & & \\ & d_2 & & \\ & & \ddots & \\ & & & d_s \\ & & & & 0 \\ & & & & & \ddots \end{pmatrix}, \quad d_1 \mid d_2 \mid \cdots \mid d_s,$$

    使得 $S = PAQ$，其中 $P \in \operatorname{GL}_m(R)$，$Q \in \operatorname{GL}_n(R)$ 是可逆矩阵。

!!! theorem "定理 48.7 (Smith 标准形的存在性)"
    设 $R$ 是 PID。任意 $m \times n$ 的 $R$-矩阵 $A$ 都有 Smith 标准形，且不变因子 $d_1, \ldots, d_s$（在相伴意义下）唯一确定。

??? proof "证明（存在性，构造法）"
    通过以下三类**初等行列变换**（对应于左乘或右乘初等矩阵，这些矩阵在 $R$ 上可逆）：

    1. 交换两行（或两列）；
    2. 将一行（或一列）加上另一行（或列）的 $R$-倍数；
    3. 将一行（或一列）乘以 $R$ 的可逆元（单位）。

    **步骤 1：** 若 $A = 0$，已是 Smith 形。否则，通过行列交换，使左上角元素 $a_{11}$ 非零，且在所有非零元素中的某种度量（例如在 $R = \mathbb{Z}$ 中取绝对值，在 $R = \mathbb{F}[x]$ 中取次数）下最小。

    **步骤 2：** 若 $a_{11} \nmid a_{1j}$（某 $j > 1$），则 $\gcd(a_{11}, a_{1j}) = d$（PID 中 gcd 存在），通过初等变换可将 $a_{11}$ 替换为 $d$（Bezout 等式 $d = ua_{11} + va_{1j}$）。类似地处理 $a_{i1}$。反复进行直到 $a_{11}$ 整除第一行和第一列的所有元素。

    **步骤 3：** 用 $a_{11}$ 消去第一行和第一列的其余元素，得到

    $$A' = \begin{pmatrix} a_{11} & 0 \\ 0 & A'' \end{pmatrix}.$$

    **步骤 4：** 若 $a_{11}$ 不整除 $A''$ 的某个元素 $a''_{ij}$，将第 $i+1$ 行加到第一行，回到步骤 2。最终 $a_{11}$ 整除 $A''$ 的所有元素。

    **步骤 5：** 对 $A''$ 递归进行上述操作。由于 $R$ 中的整除链是升链，且 PID 满足升链条件（PID 是 Noether 环），此过程必终止。

??? proof "结构定理的证明"
    设 $M$ 是有限生成 $R$-模，生成元为 $m_1, \ldots, m_n$。定义满同态

    $$\varphi: R^n \to M, \quad \varphi(r_1, \ldots, r_n) = r_1 m_1 + \cdots + r_n m_n.$$

    令 $K = \ker \varphi$。由定理 48.4，$K$ 是自由 $R$-模，秩为某个 $s \leq n$。选取 $K$ 的基 $k_1, \ldots, k_s$，将它们作为列写成 $n \times s$ 矩阵 $A$（即关系矩阵）。

    对 $A$ 做 Smith 标准形：$PAQ = S$，其中 $P \in \operatorname{GL}_n(R)$，$Q \in \operatorname{GL}_s(R)$。$P$ 对应 $R^n$ 基的变换，$Q$ 对应 $K$ 基的变换。在新基下，关系矩阵变为 $S = \operatorname{diag}(d_1, \ldots, d_s, 0, \ldots, 0)$，$d_1 \mid d_2 \mid \cdots \mid d_s$。

    这意味着

    $$M \cong R^n / K \cong R/(d_1) \oplus \cdots \oplus R/(d_s) \oplus R^{n-s}.$$

    去掉 $d_i$ 为单位的因子（此时 $R/(d_i) = 0$），令 $r = n - s$（自由秩），得到结构定理的结论。

    **唯一性**可通过分析 $p$-部分（对每个素元 $p$，考虑 $M$ 的 $p$-挠子模）来证明。

!!! example "例 48.6 (Smith 标准形的计算)"
    对 $\mathbb{Z}$-矩阵

    $$A = \begin{pmatrix} 2 & 4 & 4 \\ -6 & 6 & 12 \\ 10 & -4 & -16 \end{pmatrix},$$

    通过初等行列变换化为 Smith 标准形：

    $$S = \begin{pmatrix} 2 & 0 & 0 \\ 0 & 6 & 0 \\ 0 & 0 & 12 \end{pmatrix}.$$

    不变因子为 $d_1 = 2$，$d_2 = 6$，$d_3 = 12$，满足 $2 \mid 6 \mid 12$。

    对应的商模为

    $$\mathbb{Z}^3 / \operatorname{im}(A) \cong \mathbb{Z}/2\mathbb{Z} \oplus \mathbb{Z}/6\mathbb{Z} \oplus \mathbb{Z}/12\mathbb{Z}.$$

---

## 48.6 不变因子与初等因子

<div class="context-flow" markdown>

**核心问题**：不变因子分解和初等因子分解有什么关系？各有何优势？

</div>

!!! theorem "定理 48.8 (PID 上有限生成模的结构定理——初等因子形式)"
    设 $R$ 是主理想整环，$M$ 是有限生成 $R$-模。则

    $$M \cong R^r \oplus R/(p_1^{a_1}) \oplus R/(p_2^{a_2}) \oplus \cdots \oplus R/(p_t^{a_t}),$$

    其中 $p_1, \ldots, p_t$ 是 $R$ 中的素元（不一定两两不同），$a_1, \ldots, a_t \geq 1$。素幂 $p_1^{a_1}, \ldots, p_t^{a_t}$（在相伴和排列意义下）唯一确定，称为 $M$ 的**初等因子**（elementary divisors）。

??? proof "证明"
    关键步骤是**中国剩余定理**：若 $d = p_1^{a_1} p_2^{a_2} \cdots p_l^{a_l}$（$p_i$ 两两不相伴的素元），则

    $$R/(d) \cong R/(p_1^{a_1}) \oplus R/(p_2^{a_2}) \oplus \cdots \oplus R/(p_l^{a_l}).$$

    将不变因子形式中的每个 $R/(d_i)$ 按此方式分解，即得初等因子形式。

    反之，从初等因子可恢复不变因子：将初等因子按素元分组，$d_k$（最大的不变因子）是所有不同素元的最高幂次之积，$d_{k-1}$ 是次高幂次之积，依此类推。

!!! example "例 48.7 (不变因子与初等因子的转换)"
    考虑 $\mathbb{Z}$-模 $M$ 的不变因子为 $d_1 = 12$，$d_2 = 180$。则

    $$M \cong \mathbb{Z}/12\mathbb{Z} \oplus \mathbb{Z}/180\mathbb{Z}.$$

    将不变因子分解为素幂：

    - $12 = 2^2 \cdot 3$
    - $180 = 2^2 \cdot 3^2 \cdot 5$

    初等因子为 $2^2, 3, 2^2, 3^2, 5$，故

    $$M \cong \mathbb{Z}/4\mathbb{Z} \oplus \mathbb{Z}/3\mathbb{Z} \oplus \mathbb{Z}/4\mathbb{Z} \oplus \mathbb{Z}/9\mathbb{Z} \oplus \mathbb{Z}/5\mathbb{Z}.$$

    反之，从初等因子恢复不变因子：按素数 $2, 3, 5$ 分组，幂次分别为 $(2,2)$，$(1,2)$，$(1)$。最大不变因子 $d_2 = 2^2 \cdot 3^2 \cdot 5 = 180$，次大 $d_1 = 2^2 \cdot 3 = 12$。

!!! theorem "定理 48.9 (循环模的分解)"
    设 $R$ 是 PID，$d \in R$ 非零非单位，$d = u \cdot p_1^{a_1} \cdots p_l^{a_l}$（$u$ 为单位，$p_i$ 两两不相伴的素元）。则

    $$R/(d) \cong R/(p_1^{a_1}) \oplus \cdots \oplus R/(p_l^{a_l}).$$

    进一步，$R/(p^a)$ 是不可分解模（不能写成两个非零模的直和）。

??? proof "证明"
    **中国剩余定理的模版本：** 令 $d = d_1 d_2$，$\gcd(d_1, d_2) = 1$（即 $(d_1) + (d_2) = R$）。定义

    $$\varphi: R/(d) \to R/(d_1) \oplus R/(d_2), \quad \bar{r} \mapsto (\bar{r}_1, \bar{r}_2).$$

    由 Bezout 等式 $s d_1 + t d_2 = 1$，可构造逆映射 $(\bar{a}, \bar{b}) \mapsto \overline{a t d_2 + b s d_1}$，故 $\varphi$ 是同构。对多个互素因子递归应用即可。

    **不可分解性：** 若 $R/(p^a) \cong A \oplus B$，$A, B \neq 0$，则 $|A| \cdot |B| = |R/(p^a)|$（元素个数意义下），且 $A, B$ 的零化子理想分别为 $(p^s)$，$(p^t)$，$s, t < a$。但 $p^{\max(s,t)}$ 零化整个 $R/(p^a)$，需 $\max(s,t) \geq a$，矛盾。（严格论证需用素元的幂次分析。）

!!! example "例 48.8 (小阶Abel群的分类)"
    阶为 $n$ 的有限 Abel 群的分类等价于 $n$ 的分拆对应的初等因子分解：

    | 阶 $n$ | 不变因子 | 初等因子 | 群 |
    |:---:|:---|:---|:---|
    | 4 | (4) | $(2^2)$ | $\mathbb{Z}/4\mathbb{Z}$ |
    | 4 | (2, 2) | $(2, 2)$ | $\mathbb{Z}/2\mathbb{Z} \oplus \mathbb{Z}/2\mathbb{Z}$ |
    | 8 | (8) | $(2^3)$ | $\mathbb{Z}/8\mathbb{Z}$ |
    | 8 | (2, 4) | $(2, 2^2)$ | $\mathbb{Z}/2\mathbb{Z} \oplus \mathbb{Z}/4\mathbb{Z}$ |
    | 8 | (2, 2, 2) | $(2, 2, 2)$ | $(\mathbb{Z}/2\mathbb{Z})^3$ |
    | 12 | (12) | $(2^2, 3)$ | $\mathbb{Z}/12\mathbb{Z}$ |
    | 12 | (2, 6) | $(2, 2, 3)$ | $\mathbb{Z}/2\mathbb{Z} \oplus \mathbb{Z}/6\mathbb{Z}$ |

---

## 48.7 应用：重新推导 Jordan 形

<div class="context-flow" markdown>

**核心问题**：如何通过 $\mathbb{F}[x]$-模的结构定理统一推导有理标准形和 Jordan 标准形？

</div>

设 $V$ 是 $n$ 维 $\mathbb{F}$-向量空间，$T: V \to V$ 是线性变换。回顾例 48.1(c)，$V$ 通过 $f(x) \cdot v = f(T)(v)$ 成为 $\mathbb{F}[x]$-模。

!!! theorem "定理 48.10 ($\mathbb{F}[x]$-模的结构定理与有理标准形)"
    设 $V$ 是 $n$ 维 $\mathbb{F}$-向量空间，$T \in \operatorname{End}(V)$。视 $V$ 为 $\mathbb{F}[x]$-模。则自由秩 $r = 0$（因为 $V$ 是挠 $\mathbb{F}[x]$-模），且

    $$V \cong \mathbb{F}[x]/(d_1(x)) \oplus \mathbb{F}[x]/(d_2(x)) \oplus \cdots \oplus \mathbb{F}[x]/(d_k(x)),$$

    其中 $d_1(x) \mid d_2(x) \mid \cdots \mid d_k(x)$ 是首一多项式，$\deg d_i \geq 1$，且 $\sum_{i=1}^k \deg d_i = n$。

    不变因子 $d_1, \ldots, d_k$ 就是 $T$ 的不变因子。$d_k(x)$ 是 $T$ 的极小多项式 $m_T(x)$，$d_1(x) d_2(x) \cdots d_k(x)$ 是 $T$ 的特征多项式 $\chi_T(x)$。

??? proof "推导有理标准形"
    考虑循环模 $C_i = \mathbb{F}[x]/(d_i(x))$，其中 $d_i(x) = x^{n_i} + a_{n_i-1}x^{n_i-1} + \cdots + a_1 x + a_0$（$n_i = \deg d_i$）。$C_i$ 作为 $\mathbb{F}$-向量空间有基 $\{1, \bar{x}, \bar{x}^2, \ldots, \bar{x}^{n_i-1}\}$。

    $T$ 在 $C_i$ 上的作用就是"乘以 $x$"：

    $$T(\bar{x}^j) = \bar{x}^{j+1}, \quad j = 0, 1, \ldots, n_i - 2,$$
    $$T(\bar{x}^{n_i-1}) = -a_0 - a_1 \bar{x} - \cdots - a_{n_i-1}\bar{x}^{n_i-1}.$$

    因此 $T$ 在基 $\{1, \bar{x}, \ldots, \bar{x}^{n_i-1}\}$ 下的矩阵是**友矩阵**（companion matrix）：

    $$C(d_i) = \begin{pmatrix} 0 & 0 & \cdots & 0 & -a_0 \\ 1 & 0 & \cdots & 0 & -a_1 \\ 0 & 1 & \cdots & 0 & -a_2 \\ \vdots & & \ddots & & \vdots \\ 0 & 0 & \cdots & 1 & -a_{n_i-1} \end{pmatrix}.$$

    因此 $T$ 的矩阵在适当基下为分块对角矩阵

    $$\begin{pmatrix} C(d_1) & & \\ & \ddots & \\ & & C(d_k) \end{pmatrix},$$

    这就是 $T$ 的**有理标准形**（rational canonical form）。

!!! theorem "定理 48.11 (Jordan 标准形的模论推导)"
    设 $\mathbb{F}$ 是代数闭域（例如 $\mathbb{C}$）。将每个不变因子 $d_i(x)$ 分解为一次因子的幂：

    $$d_i(x) = (x - \lambda_{i,1})^{e_{i,1}} (x - \lambda_{i,2})^{e_{i,2}} \cdots (x - \lambda_{i,l_i})^{e_{i,l_i}}.$$

    由中国剩余定理，

    $$\mathbb{F}[x]/(d_i(x)) \cong \bigoplus_{j=1}^{l_i} \mathbb{F}[x]/((x - \lambda_{i,j})^{e_{i,j}}).$$

    因此

    $$V \cong \bigoplus_{i,j} \mathbb{F}[x]/((x - \lambda)^{e}),$$

    其中 $(x - \lambda)^e$ 遍历所有初等因子。

??? proof "从 $\mathbb{F}[x]/((x-\lambda)^e)$ 到 Jordan 块"
    考虑 $W = \mathbb{F}[x]/((x - \lambda)^e)$，它是 $e$ 维 $\mathbb{F}$-向量空间。取基

    $$v_1 = \overline{(x-\lambda)^{e-1}}, \quad v_2 = \overline{(x-\lambda)^{e-2}}, \quad \ldots, \quad v_e = \bar{1}.$$

    $T$（即乘以 $x = (x - \lambda) + \lambda$）在此基下的作用：

    $$T(v_k) = x \cdot v_k = ((x-\lambda) + \lambda) \cdot v_k = \lambda v_k + (x-\lambda) \cdot v_k.$$

    注意 $(x-\lambda) \cdot v_k = \overline{(x-\lambda)^{e-k+1}}$：

    - 若 $k > 1$：$(x-\lambda) \cdot v_k = \overline{(x-\lambda)^{e-k+1}} = v_{k-1}$；
    - 若 $k = 1$：$(x-\lambda) \cdot v_1 = \overline{(x-\lambda)^e} = 0$。

    因此

    $$T(v_1) = \lambda v_1, \quad T(v_k) = v_{k-1} + \lambda v_k, \quad k = 2, \ldots, e.$$

    $T$ 在基 $\{v_1, v_2, \ldots, v_e\}$ 下的矩阵为

    $$J_e(\lambda) = \begin{pmatrix} \lambda & 1 & & \\ & \lambda & 1 & \\ & & \ddots & 1 \\ & & & \lambda \end{pmatrix},$$

    这正是 $e \times e$ 的 **Jordan 块**。因此 $T$ 的矩阵在适当基下为 Jordan 标准形

    $$J = \begin{pmatrix} J_{e_1}(\lambda_1) & & \\ & J_{e_2}(\lambda_2) & \\ & & \ddots \end{pmatrix}.$$

!!! example "例 48.9 (从不变因子到 Jordan 形)"
    设 $T$ 的不变因子为 $d_1(x) = (x-1)(x-2)$，$d_2(x) = (x-1)^2(x-2)$。则：

    **有理标准形：**

    $$\begin{pmatrix} C(d_1) & \\ & C(d_2) \end{pmatrix} = \begin{pmatrix} 0 & -2 & & \\ 1 & 3 & & \\ & & 0 & 0 & 2 \\ & & 1 & 0 & -5 \\ & & 0 & 1 & 4 \end{pmatrix}.$$

    **初等因子：** $(x-1), (x-2), (x-1)^2, (x-2)$。

    **Jordan 标准形：**

    $$J = \begin{pmatrix} 1 & & & & \\ & 2 & & & \\ & & 1 & 1 & \\ & & & 1 & \\ & & & & 2 \end{pmatrix}.$$

    特征多项式 $\chi_T(x) = d_1(x) d_2(x) = (x-1)^3(x-2)^2$，极小多项式 $m_T(x) = d_2(x) = (x-1)^2(x-2)$。

!!! theorem "定理 48.12 (模论视角下的 Cayley-Hamilton 定理)"
    设 $V$ 是 $n$ 维 $\mathbb{F}$-向量空间，$T \in \operatorname{End}(V)$。则 $\chi_T(T) = 0$。

    **模论证明：** 视 $V$ 为 $\mathbb{F}[x]$-模。$\operatorname{Ann}(V) = (m_T(x))$（极小多项式）。由结构定理，$\chi_T(x) = \prod_i d_i(x)$，而 $d_i \mid d_k = m_T$ 对所有 $i$。故 $\chi_T(x) \in \operatorname{Ann}(V)$（因为 $m_T \mid \chi_T$），即 $\chi_T(T) = 0$。

    事实上，$m_T(x) \mid \chi_T(x)$ 可以更精确地看到：$\chi_T = d_1 d_2 \cdots d_k$，$m_T = d_k$，而 $d_1 \mid d_2 \mid \cdots \mid d_k$，故 $d_1 d_2 \cdots d_{k-1} \mid d_k^{k-1}$，因此 $m_T \mid \chi_T$。

---

## 48.8 应用：有限生成 Abel 群

<div class="context-flow" markdown>

**核心问题**：有限生成 Abel 群基本定理如何作为 $\mathbb{Z}$-模结构定理的特例得出？

</div>

!!! theorem "定理 48.13 (有限生成 Abel 群基本定理)"
    设 $G$ 是有限生成 Abel 群。则

    **(不变因子形式)** 存在唯一的 $r \geq 0$ 和整数 $d_1, \ldots, d_k \geq 2$，$d_1 \mid d_2 \mid \cdots \mid d_k$，使得

    $$G \cong \mathbb{Z}^r \oplus \mathbb{Z}/d_1\mathbb{Z} \oplus \mathbb{Z}/d_2\mathbb{Z} \oplus \cdots \oplus \mathbb{Z}/d_k\mathbb{Z}.$$

    **(初等因子形式)** 存在唯一的 $r \geq 0$ 和素幂 $p_1^{a_1}, \ldots, p_t^{a_t}$（在排列意义下唯一），使得

    $$G \cong \mathbb{Z}^r \oplus \mathbb{Z}/p_1^{a_1}\mathbb{Z} \oplus \cdots \oplus \mathbb{Z}/p_t^{a_t}\mathbb{Z}.$$

    $r$ 称为 $G$ 的**秩**（或 Betti 数），挠部分 $G_{\mathrm{tor}} = \mathbb{Z}/d_1\mathbb{Z} \oplus \cdots \oplus \mathbb{Z}/d_k\mathbb{Z}$ 是有限群。

??? proof "证明"
    这正是定理 48.6 和定理 48.8 在 $R = \mathbb{Z}$ 的情形。$\mathbb{Z}$ 是 PID，有限生成 Abel 群就是有限生成 $\mathbb{Z}$-模，直接应用结构定理即可。

!!! example "例 48.10 (有限 Abel 群的分类)"
    **(a)** 阶为 $36 = 2^2 \cdot 3^2$ 的 Abel 群的分类：

    初等因子的可能组合（对 $2$ 的幂次：$\{4\}$ 或 $\{2,2\}$；对 $3$ 的幂次：$\{9\}$ 或 $\{3,3\}$）：

    | 初等因子 | 不变因子 | 群 |
    |:---|:---|:---|
    | $4, 9$ | $(36)$ | $\mathbb{Z}/36\mathbb{Z}$ |
    | $4, 3, 3$ | $(3, 12)$ | $\mathbb{Z}/3\mathbb{Z} \oplus \mathbb{Z}/12\mathbb{Z}$ |
    | $2, 2, 9$ | $(2, 18)$ | $\mathbb{Z}/2\mathbb{Z} \oplus \mathbb{Z}/18\mathbb{Z}$ |
    | $2, 2, 3, 3$ | $(6, 6)$ | $\mathbb{Z}/6\mathbb{Z} \oplus \mathbb{Z}/6\mathbb{Z}$ |

    共 4 种互不同构的阶为 36 的 Abel 群。

    **(b)** 秩为 2、挠部分同构于 $\mathbb{Z}/12\mathbb{Z}$ 的有限生成 Abel 群为

    $$G \cong \mathbb{Z}^2 \oplus \mathbb{Z}/12\mathbb{Z} \cong \mathbb{Z}^2 \oplus \mathbb{Z}/4\mathbb{Z} \oplus \mathbb{Z}/3\mathbb{Z}.$$

!!! example "例 48.11 (应用：确定群结构)"
    设 $G$ 是阶为 $200 = 2^3 \cdot 5^2$ 的有限 Abel 群，且 $G$ 恰好有 $7$ 个阶为 $2$ 的元素。确定 $G$ 的结构。

    **解：** 阶为 2 的元素（加上单位元）构成 $G$ 的 $2$-挠子群 $G[2] = \{g \in G : 2g = 0\}$。$|G[2]| = 8$（包含 $7$ 个阶为 $2$ 的元素和单位元）。

    $G$ 的 $2$-Sylow 部分 $G_2$ 同构于 $\mathbb{Z}/2^{a_1}\mathbb{Z} \oplus \cdots \oplus \mathbb{Z}/2^{a_l}\mathbb{Z}$，$\sum a_i = 3$。$G_2[2] \cong (\mathbb{Z}/2\mathbb{Z})^l$，故 $|G_2[2]| = 2^l$。要求 $2^l = 8$，即 $l = 3$，故 $a_1 = a_2 = a_3 = 1$，$G_2 \cong (\mathbb{Z}/2\mathbb{Z})^3$。

    $G$ 的 $5$-Sylow 部分 $G_5$ 同构于 $\mathbb{Z}/5^{b_1}\mathbb{Z} \oplus \cdots$，$\sum b_j = 2$。两种可能：$G_5 \cong \mathbb{Z}/25\mathbb{Z}$ 或 $G_5 \cong (\mathbb{Z}/5\mathbb{Z})^2$。条件不足以区分，故 $G$ 有两种可能：

    $$G \cong (\mathbb{Z}/2\mathbb{Z})^3 \oplus \mathbb{Z}/25\mathbb{Z} \quad \text{或} \quad G \cong (\mathbb{Z}/2\mathbb{Z})^3 \oplus (\mathbb{Z}/5\mathbb{Z})^2.$$

---

### 本章总结

模论为线性代数提供了一个更高的视角。PID 上有限生成模的结构定理是核心结果，它统一了：

| 环 $R$ | 有限生成 $R$-模 | 结构定理的含义 |
|:---:|:---:|:---|
| $\mathbb{Z}$ | 有限生成 Abel 群 | 有限生成 Abel 群基本定理 |
| $\mathbb{F}[x]$（通过 $T$）| 线性变换的表示空间 | Jordan 标准形 / 有理标准形 |
| $\mathbb{F}[x]$（通用）| $\mathbb{F}[x]$-模 | Smith 标准形 |

这种统一性正是抽象代数的力量所在：看似不同的具体问题，在适当的抽象层次上不过是同一个定理的不同面貌。

---

### 习题

!!! exercise "习题 48.1"
    设 $R = \mathbb{Z}$，$M = \mathbb{Z}^3$，$N$ 是由 $(2, 0, 0)$，$(0, 6, 0)$，$(0, 0, 8)$ 生成的子模。求 $M/N$ 的不变因子分解和初等因子分解。

!!! exercise "习题 48.2"
    设 $T: \mathbb{R}^4 \to \mathbb{R}^4$ 的特征多项式为 $(x-1)^2(x+1)^2$，极小多项式为 $(x-1)^2(x+1)$。求 $T$ 的不变因子、有理标准形和 Jordan 标准形。

!!! exercise "习题 48.3"
    证明：设 $R$ 是 PID，$M$ 是有限生成无挠 $R$-模，则 $M$ 是自由的。

!!! exercise "习题 48.4"
    设 $G$ 是阶为 $72 = 2^3 \cdot 3^2$ 的有限 Abel 群。列出所有可能的同构类（不变因子形式和初等因子形式）。

!!! exercise "习题 48.5"
    对矩阵 $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \in M_2(\mathbb{R})$，将 $\mathbb{R}^2$ 视为 $\mathbb{R}[x]$-模（通过 $A$）。求其不变因子分解。将域扩大到 $\mathbb{C}$ 后，求初等因子和 Jordan 形。

!!! exercise "习题 48.6"
    证明：PID 上有限生成模 $M$ 的零化子 $\operatorname{Ann}(M)$ 由最大的不变因子生成。
