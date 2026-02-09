# 第 40B 章 不变量与广义矩阵函数

<div class="context-flow" markdown>

**前置**：行列式(Ch3) · 积和式(Ch40A) · 群论基础 · 对称群的表示论

**本章脉络**：对称群 $S_n$ 的表示论入门（Young 图、不可约特征标） $\to$ Immanant 定义 $\to$ Merris 广义矩阵函数 $\to$ Schur 不等式 $\to$ Lieb 定理与积和式支配猜想 $\to$ Stembridge 猜想 $\to$ Immanant 与 Schur 函数的联系 $\to$ 拉丁方与积和式 $\to$ Kasteleyn 定理与 FKT 算法 $\to$ VP 与 VNP 猜想

**延伸**：Immanant 通过对称群的表示论统一了行列式与积和式，是代数组合学的核心对象；Schur 正定性和积和式支配猜想深刻联系了矩阵论与对称函数论；Kasteleyn 定理揭示了平面性如何从根本上降低计算复杂性；VP vs VNP 猜想是代数计算复杂性的中心问题

</div>

行列式和积和式分别对应于对称群 $S_n$ 的两个"极端"特征标——符号特征标 $\operatorname{sgn}$ 和平凡特征标 $\mathbf{1}$。一个自然的问题是：如果使用 $S_n$ 的**任意**特征标 $\chi$ 来替换 $\operatorname{sgn}(\sigma)$，会得到什么？答案就是 **immanant**——一个统一了行列式和积和式的概念框架，深植于对称群的表示论之中。

本章首先建立对称群表示论的必要背景（Young 图、Specht 模、不可约特征标），然后系统发展 immanant 理论及其推广——Merris 广义矩阵函数。主要结果包括 Schur 正定性不等式、Lieb 的积和式支配定理和 Stembridge 猜想。最后，我们讨论拉丁方与积和式的组合联系、Kasteleyn 定理及代数计算复杂性中的 VP vs VNP 猜想。

---

## 40B.1 对称群 $S_n$ 的表示论入门

<div class="context-flow" markdown>

**核心问题**：理解 immanant 需要对称群表示论的哪些基本知识？

</div>

本节提供理解 immanant 理论所需的对称群表示论的自包含入门。

### 分拆与 Young 图

!!! definition "定义 40B.1 (整数分拆)"
    正整数 $n$ 的一个**分拆**（partition）是一个非递增正整数序列 $\lambda = (\lambda_1, \lambda_2, \ldots, \lambda_k)$，满足 $\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_k > 0$ 且 $\lambda_1 + \lambda_2 + \cdots + \lambda_k = n$。记为 $\lambda \vdash n$。

    $k$ 称为 $\lambda$ 的**长度**（length），记为 $\ell(\lambda) = k$。

!!! example "例 40B.1"
    $n = 4$ 的分拆有 5 个：$(4)$，$(3,1)$，$(2,2)$，$(2,1,1)$，$(1,1,1,1)$。

    $n = 5$ 的分拆有 7 个：$(5)$，$(4,1)$，$(3,2)$，$(3,1,1)$，$(2,2,1)$，$(2,1,1,1)$，$(1,1,1,1,1)$。

!!! definition "定义 40B.2 (Young 图)"
    分拆 $\lambda = (\lambda_1, \ldots, \lambda_k) \vdash n$ 的 **Young 图**（Young diagram，也称 Ferrers 图）是一个由 $n$ 个方格组成的左对齐图形，第 $i$ 行有 $\lambda_i$ 个方格。

    例如 $\lambda = (3, 2, 1) \vdash 6$ 的 Young 图为：
    $$\yng(3,2,1) \quad = \quad \begin{array}{|c|c|c|}\hline \phantom{x} & \phantom{x} & \phantom{x} \\\hline \phantom{x} & \phantom{x} & \multicolumn{1}{c}{} \\\hline \phantom{x} & \multicolumn{2}{c}{} \\\hline \end{array}$$

!!! definition "定义 40B.3 (共轭分拆)"
    分拆 $\lambda$ 的**共轭分拆** $\lambda'$ 是将 $\lambda$ 的 Young 图沿主对角线转置后得到的分拆。即 $\lambda'_j = |\{i : \lambda_i \ge j\}|$。

    例如 $(3,2,1)$ 的共轭为 $(3,2,1)$（自共轭），$(4,1)$ 的共轭为 $(2,1,1,1)$。

### 标准 Young 表与 Hook 长度公式

!!! definition "定义 40B.4 (标准 Young 表)"
    分拆 $\lambda \vdash n$ 的**标准 Young 表**（standard Young tableau，SYT）是将 $1, 2, \ldots, n$ 分别填入 $\lambda$ 的 Young 图的各格中，使得每行从左到右递增、每列从上到下递增的填法。

    分拆 $\lambda$ 的标准 Young 表的数量记为 $f^\lambda$。

!!! example "例 40B.2"
    $\lambda = (2,1) \vdash 3$ 有 $f^{(2,1)} = 2$ 个标准 Young 表：

    $$\begin{array}{|c|c|}\hline 1 & 2 \\\hline 3 & \multicolumn{1}{c}{} \\\hline \end{array} \qquad \begin{array}{|c|c|}\hline 1 & 3 \\\hline 2 & \multicolumn{1}{c}{} \\\hline \end{array}$$

!!! theorem "定理 40B.1 (Hook 长度公式)"
    对 $\lambda \vdash n$，标准 Young 表的数量为
    $$f^\lambda = \frac{n!}{\prod_{(i,j) \in \lambda} h(i,j)},$$
    其中 $h(i,j) = \lambda_i - j + \lambda'_j - i + 1$ 是位置 $(i,j)$ 的 **hook 长度**（即从 $(i,j)$ 出发，向右和向下延伸的方格总数加 1）。

!!! example "例 40B.3"
    $\lambda = (3,2) \vdash 5$。各位置的 hook 长度：

    $$\begin{array}{|c|c|c|}\hline 4 & 3 & 1 \\\hline 2 & 1 & \multicolumn{1}{c}{} \\\hline \end{array}$$

    $f^{(3,2)} = \frac{5!}{4 \cdot 3 \cdot 1 \cdot 2 \cdot 1} = \frac{120}{24} = 5$。

### 不可约表示与特征标

!!! definition "定义 40B.5 (对称群的不可约表示)"
    $S_n$ 的不可约复表示（从而不可约特征标）由 $n$ 的**分拆**参数化。对每个分拆 $\lambda \vdash n$，存在唯一（同构意义下）的不可约表示，称为 **Specht 模** $S^\lambda$。

    - $S^\lambda$ 的维数为 $f^\lambda$（标准 Young 表的数量）。
    - 对应的不可约特征标记为 $\chi^\lambda$，满足 $\chi^\lambda(\text{id}) = f^\lambda$。

!!! theorem "定理 40B.2 (对称群特征标的基本性质)"
    (a) $S_n$ 的不可约特征标 $\{\chi^\lambda : \lambda \vdash n\}$ 构成 $S_n$ 上类函数空间的正交基。

    (b) **维数公式**：$\sum_{\lambda \vdash n} (f^\lambda)^2 = n!$。

    (c) **特征标表**：$\chi^\lambda(\sigma)$ 仅取决于 $\sigma$ 的共轭类（即循环类型），可以排列成**特征标表**。特征标值均为整数。

    (d) **两个极端分拆**：

    - $\lambda = (n)$（一行）：$S^{(n)}$ 是**平凡表示**，$\chi^{(n)}(\sigma) = 1$ 对所有 $\sigma$。
    - $\lambda = (1^n) = (1,1,\ldots,1)$（一列）：$S^{(1^n)}$ 是**符号表示**，$\chi^{(1^n)}(\sigma) = \operatorname{sgn}(\sigma)$。

    (e) **Frobenius 互反律**：$\chi^\lambda$ 可以通过 Murnaghan-Nakayama 规则或 Frobenius 特征标公式从 Young 图的组合数据计算。

!!! example "例 40B.4"
    $S_3$ 的特征标表。$S_3$ 有 3 个共轭类和 3 个不可约表示：

    | 分拆 $\lambda$ | $\chi^\lambda(\text{id})$ | $\chi^\lambda((12))$ | $\chi^\lambda((123))$ |
    |---------------|-------------------------|--------------------|-----------------------|
    | $(3)$ | $1$ | $1$ | $1$ |
    | $(2,1)$ | $2$ | $0$ | $-1$ |
    | $(1,1,1)$ | $1$ | $-1$ | $1$ |

    其中共轭类的代表元分别为：恒等置换（循环类型 $(1,1,1)$）、对换（循环类型 $(2,1)$）、3-循环（循环类型 $(3)$）。

!!! example "例 40B.5"
    $S_4$ 的特征标表。$S_4$ 有 5 个共轭类（对应 $4$ 的 5 个分拆 $(4), (3,1), (2,2), (2,1,1), (1,1,1,1)$）：

    | $\lambda$ | id | $(12)$ | $(12)(34)$ | $(123)$ | $(1234)$ |
    |-----------|-----|--------|-----------|---------|----------|
    | $(4)$ | $1$ | $1$ | $1$ | $1$ | $1$ |
    | $(3,1)$ | $3$ | $1$ | $-1$ | $0$ | $-1$ |
    | $(2,2)$ | $2$ | $0$ | $2$ | $-1$ | $0$ |
    | $(2,1,1)$ | $3$ | $-1$ | $-1$ | $0$ | $1$ |
    | $(1,1,1,1)$ | $1$ | $-1$ | $1$ | $1$ | $-1$ |

### Murnaghan-Nakayama 规则

!!! theorem "定理 40B.3 (Murnaghan-Nakayama 规则)"
    设 $\sigma \in S_n$ 的循环类型为 $\mu = (\mu_1, \mu_2, \ldots, \mu_\ell)$。则
    $$\chi^\lambda(\sigma) = \sum_T (-1)^{\operatorname{ht}(T)},$$
    其中求和遍历所有**边缘带表** $T$（border-strip tableau）：从 $\lambda$ 的 Young 图出发，依次移除长度为 $\mu_1, \mu_2, \ldots, \mu_\ell$ 的**边缘带**（border strip，也称 rim hook），每步移除后剩余的仍是有效的 Young 图。$\operatorname{ht}(T)$ 是所有被移除边缘带的**高度**（= 行数 $-1$）之和。

---

## 40B.2 Immanant 的定义

<div class="context-flow" markdown>

**核心问题**：如何用对称群的表示论统一行列式和积和式？

</div>

!!! definition "定义 40B.6 (Immanant)"
    设 $\chi$ 为对称群 $S_n$ 的一个特征标（character）。矩阵 $A \in M_n(\mathbb{C})$ 的 **$\chi$-immanant** 定义为
    $$d_\chi(A) = \sum_{\sigma \in S_n} \chi(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}.$$

    当 $\chi = \chi^\lambda$ 为分拆 $\lambda \vdash n$ 对应的不可约特征标时，记 $d_\lambda(A) = d_{\chi^\lambda}(A)$。

!!! note "注记 40B.1 (Immanant 的特殊情形)"
    - 当 $\chi = \operatorname{sgn} = \chi^{(1^n)}$（符号特征标，对应交错表示）时，$d_\chi(A) = \det(A)$。
    - 当 $\chi = \mathbf{1} = \chi^{(n)}$（平凡特征标）时，$d_\chi(A) = \operatorname{perm}(A)$。
    - $S_n$ 的不可约特征标与 $n$ 的分拆一一对应。因此 $n$ 的每个分拆 $\lambda$ 给出一个不可约 immanant $d_\lambda(A)$。

    Immanant 的概念由 Littlewood 和 Richardson 在 1934 年引入。

!!! theorem "定理 40B.4 (Immanant 的基本性质)"
    设 $\chi$ 为 $S_n$ 的特征标，$A \in M_n(\mathbb{C})$。

    (a) **多重线性**：$d_\chi(A)$ 关于 $A$ 的每一行（及每一列）是线性函数。

    (b) **转置**：$d_\chi(A^T) = d_\chi(A)$。

    (c) **对角缩放**：设 $D_1 = \operatorname{diag}(d_1, \ldots, d_n)$，$D_2 = \operatorname{diag}(e_1, \ldots, e_n)$。则
    $$d_\chi(D_1 A D_2) = \left(\prod_{i=1}^n d_i e_i\right) d_\chi(A).$$

    (d) **置换行列**：对置换矩阵 $P_\tau$、$P_\rho$（$\tau, \rho \in S_n$），
    $$d_\chi(P_\tau A P_\rho) = d_\chi(A).$$
    这是因为 $\chi$ 是类函数（仅取决于共轭类），而行列置换对应于 $\sigma \mapsto \tau^{-1} \sigma \rho$ 的重参数化。

    (e) **对角矩阵**：$d_\chi(\operatorname{diag}(x_1, \ldots, x_n)) = \chi(\text{id}) \cdot x_1 x_2 \cdots x_n$。

??? proof "证明"
    (a) 与积和式的多重线性性证明完全一致——每个求和项 $\chi(\sigma) \prod a_{i,\sigma(i)}$ 关于第 $i$ 行是线性的。

    (b) $d_\chi(A^T) = \sum_\sigma \chi(\sigma) \prod a_{\sigma(i),i} = \sum_\sigma \chi(\sigma) \prod a_{i,\sigma^{-1}(i)} = \sum_\tau \chi(\tau^{-1}) \prod a_{i,\tau(i)}$。由于特征标满足 $\chi(\tau^{-1}) = \overline{\chi(\tau)}$，对实特征标（$S_n$ 的特征标值均为整数）有 $\chi(\tau^{-1}) = \chi(\tau)$，故 $d_\chi(A^T) = d_\chi(A)$。

    (c) $(D_1 A D_2)_{ij} = d_i a_{ij} e_j$，因此 $\prod_i (D_1 A D_2)_{i,\sigma(i)} = (\prod d_i)(\prod e_{\sigma(i)}) \prod a_{i,\sigma(i)} = (\prod d_i)(\prod e_j) \prod a_{i,\sigma(i)}$。

    (d) $(P_\tau A P_\rho)_{ij} = a_{\tau^{-1}(i), \rho(j)}$。因此 $\prod_i (P_\tau A P_\rho)_{i,\sigma(i)} = \prod_i a_{\tau^{-1}(i), \rho(\sigma(i))} = \prod_k a_{k, \rho\sigma\tau(k)}$（令 $k = \tau^{-1}(i)$）。令 $\phi = \rho \sigma \tau$，则求和变为 $\sum_\phi \chi(\tau^{-1} \rho^{-1} \phi) \prod a_{k,\phi(k)}$。利用 $\chi$ 是类函数和共轭不变性即可完成。

    (e) $d_\chi(\operatorname{diag}(x_1, \ldots, x_n)) = \sum_\sigma \chi(\sigma) \prod x_{\sigma(i)} \delta_{i,\sigma(i)}$。只有 $\sigma = \text{id}$ 时 $\prod \delta_{i,\sigma(i)} = 1$，故结果为 $\chi(\text{id}) \prod x_i$。

    等一下——这不对。$\operatorname{diag}(x_1,\ldots,x_n)$ 的 $(i,j)$ 元素为 $x_i \delta_{ij}$，所以 $\prod_i a_{i,\sigma(i)} = \prod_i x_i \delta_{i,\sigma(i)}$，仅当 $\sigma = \text{id}$ 时非零。故 $d_\chi(\operatorname{diag}(x_1,\ldots,x_n)) = \chi(\text{id}) \prod x_i$。

!!! example "例 40B.6"
    $n = 3$，$A = \begin{pmatrix} a & b & c \\ d & e & f \\ g & h & k \end{pmatrix}$。

    利用 $S_3$ 的特征标表（例 40B.4）计算 $d_{(2,1)}(A)$：

    $S_3$ 的元素、特征标值 $\chi^{(2,1)}$ 和对应乘积：

    | $\sigma$ | 循环类型 | $\chi^{(2,1)}(\sigma)$ | $\prod a_{i,\sigma(i)}$ |
    |----------|---------|----------------------|------------------------|
    | $()$ | $(1,1,1)$ | $2$ | $aek$ |
    | $(12)$ | $(2,1)$ | $0$ | $bdk$ |
    | $(13)$ | $(2,1)$ | $0$ | $ceh$ |
    | $(23)$ | $(2,1)$ | $0$ | $afh$ |
    | $(123)$ | $(3)$ | $-1$ | $bfg$ |
    | $(132)$ | $(3)$ | $-1$ | $cdh$ |

    $$d_{(2,1)}(A) = 2aek + 0 \cdot bdk + 0 \cdot ceh + 0 \cdot afh - bfg - cdh = 2aek - bfg - cdh.$$

    对比：$\det(A) = aek - afh - bdk + bfg + cdh - ceg$，$\operatorname{perm}(A) = aek + afh + bdk + bfg + cdh + ceg$。

!!! example "例 40B.7"
    对 $A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix}$：

    - $\operatorname{perm}(A) = d_{(3)}(A) = 450$
    - $d_{(2,1)}(A) = 2 \cdot 1 \cdot 5 \cdot 9 - 2 \cdot 6 \cdot 7 - 3 \cdot 4 \cdot 8 = 90 - 84 - 96 = -90$
    - $\det(A) = d_{(1,1,1)}(A) = 0$

---

## 40B.3 Merris 广义矩阵函数

<div class="context-flow" markdown>

**核心问题**：Immanant 可以推广到什么程度？

</div>

!!! definition "定义 40B.7 (广义矩阵函数)"
    设 $G$ 为 $S_n$ 的子群，$\chi$ 为 $G$ 的特征标。矩阵 $A \in M_n(\mathbb{C})$ 的 **$(G, \chi)$-广义矩阵函数**（generalized matrix function，Merris 意义下的）定义为
    $$d_\chi^G(A) = \sum_{\sigma \in G} \chi(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}.$$

    当 $G = S_n$ 时，退化为普通的 immanant。

!!! example "例 40B.8"
    几个特殊情形：

    (a) $G = S_n$，$\chi = \operatorname{sgn}$：$d_\chi^G(A) = \det(A)$。

    (b) $G = S_n$，$\chi = \mathbf{1}$：$d_\chi^G(A) = \operatorname{perm}(A)$。

    (c) $G = A_n$（交替群），$\chi = \mathbf{1}$：
    $$d_\chi^{A_n}(A) = \sum_{\sigma \in A_n} \prod_{i=1}^n a_{i,\sigma(i)},$$
    即只对偶数置换求和。注意 $\det(A) + \operatorname{perm}(A) = 2 d_{\mathbf{1}}^{A_n}(A)$（当取 $G = S_n$ 时的偶数置换贡献的两倍）。

    (d) $G = \{e\}$（平凡群）：$d_\chi^G(A) = \prod_{i=1}^n a_{ii}$（对角线乘积）。

!!! theorem "定理 40B.5 (广义矩阵函数的性质)"
    设 $G \le S_n$，$\chi$ 为 $G$ 的特征标。

    (a) $d_\chi^G(A)$ 关于 $A$ 的每一行是线性函数。

    (b) 若 $D = \operatorname{diag}(d_1, \ldots, d_n)$，则 $d_\chi^G(DA) = (\prod d_i) d_\chi^G(A)$。

    (c) 若 $A$ 为半正定 Hermite 矩阵且 $\chi$ 为不可约特征标，则 $d_\chi^G(A) \ge 0$（Schur 的推广）。

---

## 40B.4 Schur 正定性不等式

<div class="context-flow" markdown>

**核心问题**：半正定矩阵的 immanant 与行列式有什么定量关系？

</div>

!!! theorem "定理 40B.6 (Schur 不等式, 1918)"
    设 $A \in M_n(\mathbb{C})$ 为半正定 Hermite 矩阵，$\chi^\lambda$ 为 $S_n$ 的不可约特征标（对应分拆 $\lambda \vdash n$）。则
    $$\frac{d_\lambda(A)}{\chi^\lambda(\text{id})} \ge \det(A).$$
    即归一化的不可约 immanant 不小于行列式。等价地，
    $$d_\lambda(A) \ge f^\lambda \cdot \det(A),$$
    其中 $f^\lambda = \chi^\lambda(\text{id})$ 为标准 Young 表的数量。

??? proof "证明"
    **Schur 的原始证明（1918）** 使用了群代数 $\mathbb{C}[S_n]$ 中的正性论证。

    **步骤一：群代数元素。** 在群代数 $\mathbb{C}[S_n]$ 中定义
    $$e_\lambda = \frac{f^\lambda}{n!} \sum_{\sigma \in S_n} \chi^\lambda(\sigma) \sigma.$$
    由 Schur 正交关系，$e_\lambda$ 是群代数的**中心幂等元**（central idempotent），满足 $e_\lambda^2 = e_\lambda$ 和 $e_\lambda e_\mu = 0$（$\lambda \ne \mu$）。

    **步骤二：正定性。** 对半正定 Hermite 矩阵 $A = B^*B$，定义
    $$T = \sum_{\sigma \in S_n} \left(\prod_{i=1}^n a_{i,\sigma(i)}\right) \sigma^{-1} \in \mathbb{C}[S_n].$$
    可以验证 $T$ 是群代数中的正元素（即 $T = \sum_g c_g g$ 其中 $(c_g)$ 对应的矩阵在正则表示下半正定）。

    **步骤三：迹计算。** 在 $S^\lambda$ 上取迹：
    $$\operatorname{tr}_{S^\lambda}(T \cdot e_\lambda) = \frac{f^\lambda}{n!} \sum_{\sigma, \tau} \chi^\lambda(\tau) \left(\prod_i a_{i,\sigma(i)}\right) \operatorname{tr}_{S^\lambda}(\sigma^{-1} \tau).$$
    利用 Schur 正交关系 $\sum_\sigma \chi^\lambda(\sigma^{-1} \tau) = \frac{n!}{f^\lambda} \delta_\lambda$ 化简，最终得到
    $$d_\lambda(A) = \operatorname{tr}_{S^\lambda}(\text{某正算子}) \ge 0.$$

    **步骤四：与行列式的比较。** 利用 $\det(A) = d_{(1^n)}(A)$ 和正元素 $T$ 在不同不可约表示上的特征值关系，可以证明 $d_\lambda(A)/f^\lambda$ 是 $T$ 在表示 $S^\lambda$ 上的"平均特征值"，而 $\det(A)$ 是 $T$ 在一维符号表示上的值。由正性和表示论的结构，前者不小于后者。

!!! example "例 40B.9"
    $n = 3$，$A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 2 & 1 \\ 0 & 1 & 2 \end{pmatrix}$（正定，特征值 $2 \pm \sqrt{2}$ 和 $2$）。

    $\det(A) = 4$。

    不可约 immanant：

    - $d_{(3)}(A) = \operatorname{perm}(A)$。$\operatorname{perm}(A) = 2 \cdot 2 \cdot 2 + 2 \cdot 1 \cdot 1 + 1 \cdot 1 \cdot 2 + 1 \cdot 1 \cdot 0 + 0 \cdot 1 \cdot 2 + 0 \cdot 2 \cdot 1 = 8 + 2 + 2 + 0 + 0 + 0 = 12$。$f^{(3)} = 1$。$d_{(3)}/f^{(3)} = 12 \ge 4 = \det(A)$。

    - $d_{(2,1)}(A) = 2 \cdot (2)(2)(2) - (1)(1)(0) - (0)(1)(1) = 16 - 0 - 0 = 16$。

      让我们更仔细地计算：$d_{(2,1)}(A) = 2 \cdot a_{11}a_{22}a_{33} + 0 \cdot (对换项) - 1 \cdot (3\text{-循环项})$。

      3-循环 $(123)$: $a_{12}a_{23}a_{31} = 1 \cdot 1 \cdot 0 = 0$。
      3-循环 $(132)$: $a_{13}a_{21}a_{32} = 0 \cdot 1 \cdot 1 = 0$。

      $d_{(2,1)}(A) = 2 \cdot 8 - 0 - 0 = 16$。$f^{(2,1)} = 2$。$d_{(2,1)}/f^{(2,1)} = 8 \ge 4 = \det(A)$。

    - $d_{(1,1,1)}(A) = \det(A) = 4$。$f^{(1,1,1)} = 1$。$d_{(1,1,1)}/f^{(1,1,1)} = 4 = \det(A)$。等号成立。

    验证了 Schur 不等式。

---

## 40B.5 Lieb 定理与积和式支配猜想

<div class="context-flow" markdown>

**核心问题**：在所有 immanant 中，积和式是否总是最大的？

</div>

!!! theorem "定理 40B.7 (Lieb 定理, 1966)"
    设 $A \in M_n(\mathbb{C})$ 为半正定 Hermite 矩阵。对 $S_n$ 的任意特征标 $\chi$（不一定不可约），
    $$|d_\chi(A)| \le \chi(\text{id}) \cdot \operatorname{perm}(A).$$

??? proof "证明"
    **步骤一：分解为不可约特征标。** 设 $\chi = \sum_\lambda m_\lambda \chi^\lambda$（$m_\lambda \ge 0$ 为非负整数）。则 $\chi(\text{id}) = \sum_\lambda m_\lambda f^\lambda$，且
    $$d_\chi(A) = \sum_\lambda m_\lambda d_\lambda(A).$$

    **步骤二：不可约情形。** 对不可约特征标 $\chi^\lambda$，需要证明
    $$d_\lambda(A) \le f^\lambda \cdot \operatorname{perm}(A).$$

    这等价于证明
    $$\frac{d_\lambda(A)}{f^\lambda} \le \operatorname{perm}(A).$$

    **步骤三：利用群代数的正性。** 在群代数 $\mathbb{C}[S_n]$ 中，积和式对应于平凡表示 $\chi^{(n)}$，行列式对应于符号表示 $\chi^{(1^n)}$。

    定义正元素 $T = \sum_\sigma (\prod a_{i,\sigma(i)}) \sigma^{-1}$（如 Schur 不等式的证明中）。则 $d_\lambda(A)/f^\lambda$ 等于 $T$ 在表示 $S^\lambda$ 上的平均特征值（Rayleigh 商），而 $\operatorname{perm}(A)$ 等于 $T$ 在平凡表示上的值（也是 $T$ 的最大特征值之一）。

    关键引理：对半正定 $A$，元素 $T$ 在群代数 $\mathbb{C}[S_n]$ 的正则表示下是正半定算子。利用正则表示的分解 $\mathbb{C}[S_n] \cong \bigoplus_\lambda (S^\lambda)^{\oplus f^\lambda}$，$T$ 在每个不可约分块 $S^\lambda$ 上的作用的迹为 $d_\lambda(A)$。

    由 Schur 的正性论证（利用 $A$ 半正定时 $T$ 的正性），$T$ 的最大特征值在平凡表示上实现，给出 $\operatorname{perm}(A)$。从而每个 $d_\lambda(A)/f^\lambda$ 不超过最大特征值 $\operatorname{perm}(A)$。

    **步骤四：回到一般特征标。**
    $$|d_\chi(A)| \le \sum_\lambda m_\lambda |d_\lambda(A)| = \sum_\lambda m_\lambda d_\lambda(A) \le \sum_\lambda m_\lambda f^\lambda \operatorname{perm}(A) = \chi(\text{id}) \operatorname{perm}(A).$$
    （第一个等号利用了半正定矩阵的不可约 immanant 非负。）

!!! definition "定义 40B.8 (积和式支配猜想)"
    **积和式支配猜想**（Permanent Dominance Conjecture, Lieb）：对半正定 Hermite 矩阵 $A$ 和 $S_n$ 的任意不可约特征标 $\chi^\lambda$，
    $$\frac{d_\lambda(A)}{f^\lambda} \le \operatorname{perm}(A).$$

    Lieb 定理（定理 40B.7）证明了这个猜想。因此"猜想"在这里是历史名称——它现在是一个定理。

    积和式支配猜想与 Schur 不等式结合，给出了 immanant 的完整层次结构：对半正定 $A$，
    $$\det(A) \le \frac{d_\lambda(A)}{f^\lambda} \le \operatorname{perm}(A).$$

!!! note "注记 40B.2 (积和式支配的精细化)"
    积和式支配猜想在更精细的层面上可以表述为：对分拆偏序（dominance order）$\lambda \unrhd \mu$ 中 $\lambda$ 支配 $\mu$ 时，半正定矩阵 $A$ 是否满足 $d_\lambda(A)/f^\lambda \ge d_\mu(A)/f^\mu$？这一更强的猜想至今尚未完全解决。

---

## 40B.6 Stembridge 猜想

<div class="context-flow" markdown>

**核心问题**：对特殊类型的矩阵（Jacobi-Trudi 矩阵），immanant 有什么更精确的结论？

</div>

!!! definition "定义 40B.9 (Jacobi-Trudi 矩阵)"
    设 $\lambda = (\lambda_1, \ldots, \lambda_n) \vdash N$ 为分拆（允许 $\lambda_i = 0$），$h_k$ 为第 $k$ 个**完全齐次对称函数**（complete homogeneous symmetric function）。**Jacobi-Trudi 矩阵**定义为
    $$H_\lambda = (h_{\lambda_i - i + j})_{1 \le i, j \le n},$$
    其中约定 $h_k = 0$（$k < 0$），$h_0 = 1$。

    由 Jacobi-Trudi 恒等式，$s_\lambda = \det(H_\lambda)$，其中 $s_\lambda$ 为 Schur 函数。

!!! theorem "定理 40B.8 (Stembridge 猜想/定理, 1992)"
    设 $H_\mu$ 为分拆 $\mu$ 对应的 Jacobi-Trudi 矩阵。则对 $S_n$ 的任意不可约特征标 $\chi^\lambda$，
    $$d_\lambda(H_\mu)$$
    作为对称函数，是 Schur 函数的**非负线性组合**（即 Schur 正的）。

    等价地，$d_\lambda(H_\mu) = \sum_\nu c^\nu_{\lambda,\mu} s_\nu$，其中 $c^\nu_{\lambda,\mu} \ge 0$。

    这一猜想由 Stembridge 于 1992 年提出。对多种特殊情形已被证明：

    - $\lambda = (n)$：$d_{(n)}(H_\mu) = \operatorname{perm}(H_\mu)$，其 Schur 正性由 Stembridge 本人证明。
    - $\lambda = (n-1, 1)$：即"标准"immanant，由 Greene (1992) 和 Haiman (1993) 证明。
    - 一般的 hook 形分拆 $\lambda = (n-k, 1^k)$：由 Haiman 证明。

!!! note "注记 40B.3 (Stembridge 猜想的意义)"
    Stembridge 猜想将 immanant 理论与**对称函数论**中最精细的正性问题联系起来。Schur 正性是比普通非负性更强的条件——它不仅要求系数非负，而且要求在 Schur 基下展开的系数非负。

    这一猜想的完整证明将统一 immanant 理论中的多个已知结果，并为对称群表示论与矩阵不等式的交互提供更深层的理解。

---

## 40B.7 Immanant 与 Schur 函数

<div class="context-flow" markdown>

**核心问题**：Immanant 与对称函数论中的 Schur 函数有什么深层联系？

</div>

### Schur 函数回顾

!!! definition "定义 40B.10 (Schur 函数)"
    对分拆 $\lambda = (\lambda_1 \ge \cdots \ge \lambda_n \ge 0)$（允许零部分使长度为 $n$），**Schur 多项式**（Schur polynomial）在 $n$ 个变量 $x_1, \ldots, x_n$ 上定义为
    $$s_\lambda(x_1, \ldots, x_n) = \frac{\det(x_j^{\lambda_i + n - i})_{1 \le i,j \le n}}{\det(x_j^{n-i})_{1 \le i,j \le n}} = \frac{\det(x_j^{\lambda_i + n - i})}{\prod_{i < j}(x_i - x_j)}.$$
    分母是 Vandermonde 行列式 $\Delta(x) = \prod_{i<j}(x_i - x_j)$。

    等价地，$s_\lambda$ 可以定义为半标准 Young 表的生成函数：
    $$s_\lambda(x_1, \ldots, x_n) = \sum_{T \in \text{SSYT}(\lambda, n)} x^T,$$
    其中 SSYT$(\lambda, n)$ 是形状为 $\lambda$ 且元素来自 $\{1, \ldots, n\}$ 的半标准 Young 表的集合，$x^T = \prod_i x_i^{c_i(T)}$（$c_i(T)$ 为 $i$ 在 $T$ 中出现的次数）。

### Immanant 与 Schur 函数的联系

!!! theorem "定理 40B.9 (Immanant 对对角矩阵的值)"
    设 $A = \operatorname{diag}(x_1, \ldots, x_n)$。则对分拆 $\lambda \vdash n$，
    $$d_\lambda(A) = f^\lambda \cdot e_1^n(x_1, \ldots, x_n) = \chi^\lambda(\text{id}) \cdot x_1 x_2 \cdots x_n,$$
    这是 $x_1 \cdots x_n$ 的 $f^\lambda$ 倍。

    更一般地，对**幂和矩阵** $A = (x_j^{k_i})$（即 $a_{ij} = x_j^{k_i}$），immanant 与 Schur 函数之间有深刻的联系。

!!! theorem "定理 40B.10 (Immanant 的 Schur 函数展开)"
    设 $A = (a_{ij})$ 为矩阵，其元素可以表示为对称函数。对分拆 $\lambda \vdash n$，归一化 immanant
    $$\bar{d}_\lambda(A) = \frac{d_\lambda(A)}{f^\lambda}$$
    与 Schur 函数之间的关系由以下恒等式给出（Littlewood-Richardson 理论框架）：

    在群代数 $\mathbb{C}[S_n]$ 中，Immanant 对应的投影算子为
    $$\hat{e}_\lambda = \frac{f^\lambda}{n!} \sum_{\sigma \in S_n} \chi^\lambda(\sigma) \sigma.$$

    当 $A$ 的元素具有组合意义时（如完全齐次对称函数），$d_\lambda(A)$ 可以用 Schur 函数展开，展开系数由 Littlewood-Richardson 系数和特征标值决定。

!!! theorem "定理 40B.11 (Schur 函数与 Immanant 的联系：Jacobi-Trudi 形式)"
    对 $n$ 个变量的 Schur 多项式，由 Jacobi-Trudi 恒等式：
    $$s_\lambda(x_1, \ldots, x_n) = \det(h_{\lambda_i - i + j}(x_1, \ldots, x_n))_{1 \le i,j \le n},$$
    其中 $h_k$ 为完全齐次对称多项式。

    将行列式替换为 immanant，得到**广义 Schur 函数**：
    $$s_\lambda^{(\mu)}(x) = \frac{1}{f^\mu} d_\mu\left((h_{\lambda_i - i + j})_{i,j}\right).$$

    这是 Stembridge 猜想中研究的核心对象。当 $\mu = (1^n)$ 时，$s_\lambda^{(1^n)} = s_\lambda$（普通 Schur 函数）。当 $\mu = (n)$ 时，$s_\lambda^{(n)} = \operatorname{perm}(H_\lambda)/1 = \operatorname{perm}(H_\lambda)$。

!!! note "注记 40B.4 (表示论视角)"
    从表示论的角度看，immanant 可以理解为 $\operatorname{GL}_n(\mathbb{C})$ 的张量表示 $V^{\otimes n}$ 中，在 $S_n$ 的不同不可约表示上的投影。

    具体地，$n$ 阶张量空间 $V^{\otimes n}$（$V = \mathbb{C}^n$）在 $\operatorname{GL}_n$ 和 $S_n$ 的联合作用下分解为
    $$V^{\otimes n} \cong \bigoplus_{\lambda \vdash n} V_\lambda \otimes S^\lambda,$$
    其中 $V_\lambda$ 是 $\operatorname{GL}_n$ 的不可约表示（对应权 $\lambda$），$S^\lambda$ 是 $S_n$ 的 Specht 模（**Schur-Weyl 对偶**）。

    行列式对应于反对称张量 $\bigwedge^n V = V_{(1^n)}$，积和式对应于对称张量 $\operatorname{Sym}^n V = V_{(n)}$，一般的 immanant 对应于 $V_\lambda$。

    Schur 函数 $s_\lambda(x_1, \ldots, x_n)$ 恰好是 $V_\lambda$ 的特征标——即 $\operatorname{GL}_n$ 在 $V_\lambda$ 上作用的迹关于对角元素的多项式。

---

## 40B.8 拉丁方与积和式

<div class="context-flow" markdown>

**核心问题**：积和式如何与拉丁方的计数联系？

</div>

!!! definition "定义 40B.11 (拉丁方与拉丁矩形)"
    (a) $n$ 阶**拉丁方**（Latin square）是一个 $n \times n$ 的数组，其中每一行和每一列都恰好是 $\{1, 2, \ldots, n\}$ 的一个排列。

    (b) $r \times n$ **拉丁矩形**（$r \le n$）是一个 $r \times n$ 的数组，每一行是 $\{1, \ldots, n\}$ 的排列，每一列元素互不相同。

!!! theorem "定理 40B.12 (拉丁矩形扩展与积和式)"
    设 $L$ 为 $r \times n$ 拉丁矩形（$r < n$）。定义 $(0,1)$-矩阵 $A_L \in M_n$ 如下：
    $$a_{ij} = \begin{cases} 1, & \text{若 } j \text{ 不在 } L \text{ 的第 } i \text{ 列中出现} \\ 0, & \text{否则} \end{cases} \quad (i = 1, \ldots, n).$$

    实际上，$A_L$ 是 $(n-r)$-正则的 $(0,1)$-矩阵（每行每列之和为 $n-r$）。从 $r$ 行拉丁矩形扩展到 $r+1$ 行的方式数等于 $\operatorname{perm}(A_L)$。

??? proof "证明"
    扩展 $L$ 的第 $r+1$ 行意味着选择一个排列 $\sigma \in S_n$，使得 $\sigma(i)$ 不在 $L$ 的第 $i$ 列中出现。这等价于选择 $\sigma$ 使得 $a_{i,\sigma(i)} = 1$ 对所有 $i$ 成立，即选择 $A_L$ 的一个完美匹配。由积和式与完美匹配的对应关系，方式数等于 $\operatorname{perm}(A_L)$。

!!! theorem "定理 40B.13 (拉丁方计数的递推)"
    $n$ 阶拉丁方的数量 $L(n)$ 满足
    $$L(n) = n! \cdot \prod_{r=1}^{n-1} \operatorname{perm}(A_{L_r}),$$
    其中（非精确表述）求和/乘积遍历所有可能的中间拉丁矩形 $L_r$。

    更精确地，
    $$L(n) = n! \cdot (n-1)! \cdot R(n),$$
    其中 $R(n)$ 为**约化拉丁方**（第一行和第一列固定为 $1, 2, \ldots, n$）的数量。$R(n)$ 可以用积和式递推：
    $$R(n) = \sum_{\text{configurations}} \prod_{r=2}^{n-1} \operatorname{perm}(B_r),$$
    其中 $B_r$ 是由已填入的前 $r$ 行确定的 $(0,1)$-矩阵。

    已知精确值：$L(1) = 1$, $L(2) = 2$, $L(3) = 12$, $L(4) = 576$, $L(5) = 161280$, $L(6) = 812851200$。

!!! example "例 40B.10"
    $n = 3$。一个 $1 \times 3$ 拉丁矩形 $L = (1, 2, 3)$。

    $A_L = \begin{pmatrix} 0 & 1 & 1 \\ 1 & 0 & 1 \\ 1 & 1 & 0 \end{pmatrix}$（第 $i$ 列缺少数字 $L_{1,i}$）。

    $\operatorname{perm}(A_L) = 0 + 0 + 0 + 1 + 1 + 0 = 2$。所以有 2 种方式扩展到 2 行拉丁矩形。

    继续扩展到第 3 行：对每种 2 行拉丁矩形，对应的 $(0,1)$-矩阵是置换矩阵（每行每列恰有 1 个 1），积和式为 1。

    因此从 $(1,2,3)$ 出发可以构造 $2 \times 1 = 2$ 个约化拉丁方。总共 $L(3) = 3! \cdot 2 = 12$。

!!! note "注记 40B.5 (Van der Waerden 猜想与拉丁方)"
    由 Van der Waerden 猜想（定理 40A.5），$(n-r)$-正则 $(0,1)$-矩阵的积和式有下界。这给出了 $L(n)$ 的下界估计。利用 Stirling 近似，可以得到
    $$L(n) \ge \left(\frac{n}{e^2}\right)^{n^2} \cdot (\text{低阶项}).$$

---

## 40B.9 Kasteleyn 定理与 FKT 算法

<div class="context-flow" markdown>

**核心问题**：平面图上的完美匹配计数是否可以高效完成？

</div>

虽然一般图上的完美匹配计数是 $\#P$-完全的，Kasteleyn (1961) 和 Temperley-Fisher (1961) 发现了一个深刻的结果：**平面图**上的完美匹配数可以在多项式时间内计算。

### Pfaffian 回顾

!!! definition "定义 40B.12 (Pfaffian)"
    对 $2n \times 2n$ 反对称矩阵 $B$（$B^T = -B$），**Pfaffian** 定义为
    $$\operatorname{pf}(B) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \operatorname{sgn}(\sigma) \prod_{i=1}^n b_{\sigma(2i-1), \sigma(2i)}.$$

    Pfaffian 与行列式的关系为 $\operatorname{pf}(B)^2 = \det(B)$。因此 $|\operatorname{pf}(B)| = \sqrt{|\det(B)|}$，可以在 $O(n^3)$ 时间内计算。

!!! note "注记 40B.6 (Pfaffian 与 Hafnian 的对比)"
    Pfaffian 和 Hafnian 的定义非常相似：

    - **Pfaffian**：$\operatorname{pf}(B) = \sum_{M \in \mathcal{M}_{2n}} \operatorname{sgn}(M) \prod_{\{i,j\} \in M} b_{ij}$（带符号）
    - **Hafnian**：$\operatorname{hf}(A) = \sum_{M \in \mathcal{M}_{2n}} \prod_{\{i,j\} \in M} a_{ij}$（无符号）

    其中 $\operatorname{sgn}(M)$ 是匹配 $M$ 对应的置换的符号。

    这完全类似于行列式与积和式的关系——Pfaffian 是"带符号的 Hafnian"。由于符号的存在，Pfaffian 可以利用消除法（通过 $\det(B) = \operatorname{pf}(B)^2$）高效计算，而 Hafnian 是 $\#P$-完全的。

### Kasteleyn 定理

!!! theorem "定理 40B.14 (Kasteleyn 定理, 1961)"
    设 $G$ 为**平面**图（可以画在平面上无交叉边）。则存在 $G$ 的一个**定向** $\vec{G}$（对每条边指定方向），使得
    $$\text{完美匹配数} = |\operatorname{pf}(B_{\vec{G}})|,$$
    其中 $B_{\vec{G}}$ 是定向图 $\vec{G}$ 的反对称邻接矩阵：
    $$(B_{\vec{G}})_{ij} = \begin{cases} +1, & \text{若 } (i \to j) \in \vec{G} \\ -1, & \text{若 } (j \to i) \in \vec{G} \\ 0, & \text{若 } \{i,j\} \notin E(G) \end{cases}$$

    这样的定向称为 **Kasteleyn 定向**（Pfaffian orientation）。

??? proof "证明"
    **关键思想**：Hafnian 与 Pfaffian 的差异在于符号。如果我们能选择边的方向，使得对每个完美匹配 $M$，对应的符号 $\operatorname{sgn}(M)$ 都相同（都为 $+1$ 或都为 $-1$），则 $|\operatorname{pf}(B_{\vec{G}})| = \operatorname{hf}(A) = $ 完美匹配数。

    **Kasteleyn 条件**：定向 $\vec{G}$ 是 Kasteleyn 定向，当且仅当 $G$ 的每个（长度为偶数的）面的边界上，顺时针方向的边数为奇数。

    **为什么平面图可以满足此条件**：对平面图，利用其面结构（由 Euler 公式 $V - E + F = 2$），可以贪心地逐面调整边的方向。从任意定向开始，对每个面检查条件；若不满足，翻转该面边界上的一条边。对平面图，这些条件是"几乎独立的"（因为面之间只共享边界边），可以通过系统的方法（例如使用生成树互补的方式）构造满足条件的定向。

    **非平面图**：对非平面图（例如完全二部图 $K_{3,3}$ 或 Petersen 图），一般不存在 Kasteleyn 定向。实际上，Kasteleyn 定向的存在性与图的 genus（亏格）有关：亏格为 $g$ 的图需要 $4^g$ 个 Pfaffian 的组合来计算完美匹配数。

!!! example "例 40B.11"
    $2 \times 3$ 棋盘格（6 个顶点的网格图），是平面二部图。

    构造 Kasteleyn 定向：将水平边向右定向，垂直边交替向上向下定向，使得每个面（$1 \times 1$ 正方形）的边界上有奇数条顺时针边。

    反对称邻接矩阵 $B \in M_6$，$|\operatorname{pf}(B)| = \sqrt{|\det(B)|}$。

    $\det(B)$ 可以直接计算（$O(6^3)$ 时间），得到完美匹配数为 3。

### FKT 算法

!!! theorem "定理 40B.15 (FKT 算法)"
    **Fisher-Kasteleyn-Temperley (FKT) 算法**：给定 $n$ 个顶点的平面图 $G$，计算其完美匹配数的算法如下：

    1. **平面嵌入**：计算 $G$ 的一个平面嵌入，确定所有面。时间 $O(n)$。
    2. **Kasteleyn 定向**：构造满足 Kasteleyn 条件的定向 $\vec{G}$。时间 $O(n)$。
    3. **计算 Pfaffian**：构造反对称邻接矩阵 $B_{\vec{G}}$，计算 $\det(B_{\vec{G}})$。时间 $O(n^3)$。
    4. **输出**：$\sqrt{|\det(B_{\vec{G}})|}$ 即为完美匹配数。

    总时间复杂度：$O(n^3)$。

!!! example "例 40B.12"
    $m \times n$ 棋盘格上的二聚体覆盖（domino tiling）数。

    $2 \times n$ 棋盘格：覆盖数为 Fibonacci 数 $F_{n+1}$。例如 $2 \times 3$ 格有 $F_4 = 3$ 种覆盖。

    $m \times n$ 棋盘格（$mn$ 为偶数）的二聚体覆盖数由 Kasteleyn-Temperley-Fisher 公式给出：
    $$\prod_{j=1}^{\lceil m/2 \rceil} \prod_{k=1}^{\lceil n/2 \rceil} \left(4\cos^2\frac{\pi j}{m+1} + 4\cos^2\frac{\pi k}{n+1}\right).$$

    例如 $8 \times 8$ 棋盘格有 $12988816$ 种二聚体覆盖。

---

## 40B.10 VP 与 VNP 猜想

<div class="context-flow" markdown>

**核心问题**：行列式与积和式之间的计算差异有什么代数几何层面的深层含义？

</div>

!!! definition "定义 40B.13 (代数计算复杂类 VP 和 VNP)"
    **Valiant 的代数计算复杂类**（1979）：

    (a) **VP**（Valiant's $P$-class）：多项式族 $(f_n)_{n \ge 1}$（$f_n$ 是 $\operatorname{poly}(n)$ 个变量的 $\operatorname{poly}(n)$ 次多项式）属于 VP，若 $f_n$ 可以由 $\operatorname{poly}(n)$ 大小的**代数电路**（arithmetic circuit，仅使用 $+, \times$ 和常数）计算。

    (b) **VNP**（Valiant's $NP$-class）：多项式族 $(f_n)$ 属于 VNP，若
    $$f_n(x_1, \ldots, x_m) = \sum_{e \in \{0,1\}^t} g_n(x_1, \ldots, x_m, e_1, \ldots, e_t),$$
    其中 $g_n \in$ VP 且 $t = \operatorname{poly}(n)$。即 $f_n$ 是某个 VP 多项式对指数多个布尔变量的求和。

!!! theorem "定理 40B.16 (Valiant 的 VP 与 VNP 结果)"
    (a) $n \times n$ 矩阵的行列式 $\det_n$ 属于 VP。（由 Gauss 消元的代数电路实现给出。）

    (b) $n \times n$ 矩阵的积和式 $\operatorname{perm}_n$ 属于 VNP。

    (c) $\operatorname{perm}_n$ 是 **VNP-完全**的：VNP 中的任何多项式族都可以通过 VP 归约到积和式。

!!! definition "定义 40B.14 (VP vs VNP 猜想)"
    **VP $\ne$ VNP 猜想**（Valiant, 1979）：积和式不属于 VP。即不存在 $\operatorname{poly}(n)$ 大小的代数电路计算 $n \times n$ 矩阵的积和式。

    等价表述：**不存在** $\operatorname{poly}(n)$ 阶矩阵 $M$（其元素是 $a_{ij}$ 和额外辅助变量的仿射函数）使得 $\det(M) = \operatorname{perm}(A)$。

    这是 $P \ne NP$ 猜想的**代数类比**，被认为是代数计算复杂性理论中最核心的开放问题。

!!! note "注记 40B.7 (VP vs VNP 的已知进展)"
    - **Mignon-Ressayre (2004)**：证明了表示 $n \times n$ 积和式所需的行列式阶至少为 $n^2/2$。即若 $\det(M) = \operatorname{perm}_n$，则 $M$ 的阶至少为 $\Omega(n^2)$。
    - **Landsberg-Manivel-Ressayre (2013)**：利用表示论方法改进了下界。
    - **Grochow-Moore (2014)**：利用对称性论证给出了新的限制。
    - 目前已知的最好下界离 $\operatorname{poly}(n)$ 还有很大差距。完全证明 VP $\ne$ VNP 需要本质性的新思想。

!!! note "注记 40B.8 (几何复杂性理论)"
    Mulmuley 和 Sohoni (2001) 提出了**几何复杂性理论**（Geometric Complexity Theory, GCT）纲领，试图利用代数几何和表示论来证明 VP $\ne$ VNP。核心思想是：

    1. 行列式 $\det_m$ 在 $\operatorname{GL}_{m^2}$ 的作用下生成一个射影簇 $\overline{\operatorname{GL} \cdot \det_m}$（行列式的轨道闭包）。
    2. 若 $\operatorname{perm}_n$ 不在这个闭包中（对所有 $m = \operatorname{poly}(n)$），则 VP $\ne$ VNP。
    3. 两个轨道闭包的区分可以通过它们的**坐标环**（coordinate ring）中的不可约表示的差异来实现。

    这一纲领将 VP vs VNP 问题转化为表示论中的重数问题，使得 immanant 理论、Schur 函数和 plethysm 等工具成为潜在的攻击手段。

---

## 40B.11 习题

!!! exercise "习题 40B.1"
    列出 $n = 5$ 的所有分拆，画出对应的 Young 图，并计算每个分拆的标准 Young 表数量 $f^\lambda$（利用 Hook 长度公式）。验证 $\sum (f^\lambda)^2 = 5! = 120$。

!!! exercise "习题 40B.2"
    写出 $S_4$ 的完整特征标表（5 个共轭类，5 个不可约特征标）。验证正交关系：

    (a) 行正交：$\frac{1}{24} \sum_\sigma \chi^\lambda(\sigma) \overline{\chi^\mu(\sigma)} = \delta_{\lambda\mu}$。

    (b) 列正交：$\sum_\lambda \chi^\lambda(C_1) \overline{\chi^\lambda(C_2)} = \frac{|S_4|}{|C_1|} \delta_{C_1, C_2}$。

!!! exercise "习题 40B.3"
    对 $A = \begin{pmatrix} 1 & 2 & 3 & 4 \\ 5 & 6 & 7 & 8 \\ 9 & 10 & 11 & 12 \\ 13 & 14 & 15 & 16 \end{pmatrix}$，计算所有 5 个不可约 immanant $d_\lambda(A)$（$\lambda \vdash 4$）。

!!! exercise "习题 40B.4"
    (Schur 不等式验证) 设 $A = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$。

    (a) 计算 $\det(A)$ 和 $\operatorname{perm}(A)$。
    (b) 验证 $\operatorname{perm}(A) \ge \det(A)$。
    (c) $S_2$ 只有两个不可约表示。Schur 不等式在此情形变成什么？

!!! exercise "习题 40B.5"
    利用 Murnaghan-Nakayama 规则计算 $\chi^{(3,2)}$ 在 $S_5$ 的各共轭类上的值。

!!! exercise "习题 40B.6"
    证明：对任意矩阵 $A \in M_n(\mathbb{C})$，
    $$\sum_{\lambda \vdash n} f^\lambda \cdot d_\lambda(A) = n! \cdot \prod_{i=1}^n a_{ii}.$$
    （提示：利用 $\sum_\lambda f^\lambda \chi^\lambda(\sigma) = 0$（$\sigma \ne \text{id}$）和 $\sum_\lambda (f^\lambda)^2 = n!$。）

!!! exercise "习题 40B.7"
    (Lieb 定理应用) 设 $A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 2 & 1 \\ 0 & 1 & 2 \end{pmatrix}$（正定）。

    (a) 计算 $\operatorname{perm}(A)$，$\det(A)$，$d_{(2,1)}(A)$。
    (b) 验证 $\det(A) \le d_{(2,1)}(A)/f^{(2,1)} \le \operatorname{perm}(A)$。

!!! exercise "习题 40B.8"
    (拉丁方) 构造所有 $3$ 阶约化拉丁方（第一行和第一列固定为 $1, 2, 3$），验证 $R(3) = 1$，从而 $L(3) = 3! \cdot 2! \cdot 1 = 12$。

!!! exercise "习题 40B.9"
    (Kasteleyn 定理) 对 $2 \times 4$ 棋盘格（8 个顶点）：

    (a) 画出对应的平面二部图。
    (b) 构造一个 Kasteleyn 定向。
    (c) 写出反对称邻接矩阵 $B$，计算 $\det(B)$，从而得到二聚体覆盖数。
    (d) 手动枚举所有覆盖方式，验证结果。

!!! exercise "习题 40B.10"
    (FKT 算法) 证明：$2 \times n$ 棋盘格的二聚体覆盖数等于 Fibonacci 数 $F_{n+1}$。（提示：建立递推关系。）

!!! exercise "习题 40B.11"
    (VP vs VNP) 证明：$2 \times 2$ 积和式 $\operatorname{perm}_2 = a_{11}a_{22} + a_{12}a_{21}$ 可以表示为 $3 \times 3$ 矩阵的行列式（即找到 $M \in M_3$ 使得 $\det(M) = \operatorname{perm}_2$，其中 $M$ 的元素是 $a_{ij}$ 的仿射函数）。

!!! exercise "习题 40B.12"
    (Schur 函数与 Immanant) 对 $n = 3$ 和对角矩阵 $A = \operatorname{diag}(x, y, z)$：

    (a) 验证 $d_{(3)}(A) = xyz$，$d_{(1,1,1)}(A) = xyz$，$d_{(2,1)}(A) = 2xyz$。
    (b) 解释为什么所有不可约 immanant 在对角矩阵上都是 $f^\lambda \cdot x_1 \cdots x_n$。

!!! exercise "习题 40B.13"
    (Stembridge 猜想特殊情形) 对 $\lambda = (2,1)$（$n = 3$），计算 Jacobi-Trudi 矩阵 $H_{(2,1)} = \begin{pmatrix} h_2 & h_3 \\ 1 & h_1 \end{pmatrix}$ 的积和式 $\operatorname{perm}(H_{(2,1)}) = h_1 h_2 + h_3$，并将其展开为 Schur 函数的线性组合，验证系数非负。

!!! exercise "习题 40B.14"
    证明：对 $n \times n$ 矩阵 $A$ 和 $S_n$ 的任意特征标 $\chi$，
    $$d_\chi(A \otimes B) = d_\chi(A) \cdot d_\chi(B)$$
    **不**一般成立。构造反例。

!!! exercise "习题 40B.15"
    (综合) 设 $G$ 为 Petersen 图（10 个顶点，15 条边，3-正则）。

    (a) 解释为什么 FKT 算法不能直接用于 $G$（提示：$G$ 是非平面图）。
    (b) 利用 $G$ 的对称性和分支定界方法，计算 $G$ 的完美匹配数。
    (c) Petersen 图的完美匹配数与其邻接矩阵（视为一般图）的 Hafnian 有什么关系？
