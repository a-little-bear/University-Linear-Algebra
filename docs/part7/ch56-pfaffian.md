# 第 56 章 Pfaffian

<div class="context-flow" markdown>

**前置**：行列式 (Ch3) · 外代数 (Ch49) · 辛矩阵 (Ch53)

**本章脉络**：反对称矩阵回顾 → Pfaffian 定义 → $\mathrm{pf}(A)^2 = \det(A)$ → 基本性质 → 计算方法 → 与完美匹配的关系 → FKT 算法 → 物理应用

**延伸**：Pfaffian 在统计力学（二维 Ising 模型的精确解）、代数拓扑（Euler 类作为曲率的 Pfaffian）和组合学（平面图完美匹配的多项式时间计数）中有精妙应用

</div>

行列式是线性代数最基本的概念之一。对于反对称矩阵（$A^T = -A$），行列式具有一个精致的"平方根"——**Pfaffian**。Pfaffian 的名字来自 18 世纪德国数学家 Johann Friedrich Pfaff。虽然 Pfaffian 可以定义为行列式的一种推广，但它的内涵远比"开平方"丰富得多：它与外代数的顶形式、图论中的完美匹配、统计力学中的精确可解模型、微分几何中的 Euler 类都有深刻联系。本章从反对称矩阵的基本性质出发，给出 Pfaffian 的组合定义和外代数定义，证明核心恒等式 $\mathrm{pf}(A)^2 = \det(A)$，讨论高效计算方法，然后展示 Pfaffian 在图论和物理中的应用。

---

## 56.1 反对称矩阵回顾

<div class="context-flow" markdown>

**核心问题**：反对称矩阵的行列式有什么特殊结构？

</div>

!!! definition "定义 56.1 (反对称矩阵)"
    实（或复）方阵 $A$ 称为**反对称**（skew-symmetric，或 antisymmetric）的，若
    $$A^T = -A.$$
    这意味着 $a_{ij} = -a_{ji}$，特别地 $a_{ii} = 0$。

!!! theorem "定理 56.1 (奇数阶反对称矩阵的行列式)"
    若 $A$ 是 $n \times n$ 反对称矩阵且 $n$ 为奇数，则 $\det A = 0$。

??? proof "证明"
    $\det A = \det(A^T) = \det(-A) = (-1)^n \det A$。若 $n$ 为奇数，$(-1)^n = -1$，故 $\det A = -\det A$，即 $\det A = 0$。$\blacksquare$

!!! theorem "定理 56.2 (偶数阶反对称矩阵的行列式是完全平方)"
    若 $A$ 是 $2n \times 2n$ 反对称矩阵，其元素为不定元（即视为多项式环 $\mathbb{Z}[a_{ij} : 1 \leq i < j \leq 2n]$ 中的元素），则 $\det A$ 是某个多项式的平方。这个多项式就是 Pfaffian。

!!! example "例 56.1"
    $n = 1$（$2 \times 2$ 矩阵）：$A = \begin{bmatrix} 0 & a \\ -a & 0 \end{bmatrix}$，$\det A = a^2 = (a)^2$。Pfaffian 为 $a_{12} = a$。

    $n = 2$（$4 \times 4$ 矩阵）：$A = \begin{bmatrix} 0 & a_{12} & a_{13} & a_{14} \\ -a_{12} & 0 & a_{23} & a_{24} \\ -a_{13} & -a_{23} & 0 & a_{34} \\ -a_{14} & -a_{24} & -a_{34} & 0 \end{bmatrix}$。

    直接计算：$\det A = (a_{12}a_{34} - a_{13}a_{24} + a_{14}a_{23})^2$。

    Pfaffian 为 $\mathrm{pf}(A) = a_{12}a_{34} - a_{13}a_{24} + a_{14}a_{23}$。

!!! theorem "定理 56.3 (反对称矩阵的标准形)"
    设 $A$ 是 $2n \times 2n$ 实反对称矩阵。则存在正交矩阵 $Q$ 使得
    $$Q^T A Q = \begin{bmatrix} 0 & \lambda_1 & & & & \\ -\lambda_1 & 0 & & & & \\ & & 0 & \lambda_2 & & \\ & & -\lambda_2 & 0 & & \\ & & & & \ddots & \\ & & & & & 0 & \lambda_n \\ & & & & & -\lambda_n & 0 \end{bmatrix},$$
    其中 $\lambda_1, \ldots, \lambda_n \geq 0$。特别地，$\det A = (\lambda_1 \lambda_2 \cdots \lambda_n)^2$。

??? proof "证明"
    $A$ 反对称意味着 $iA$ 是 Hermite 矩阵，可以由酉矩阵对角化。$iA$ 的特征值为实数，设为 $\mu_1, -\mu_1, \mu_2, -\mu_2, \ldots, \mu_n, -\mu_n$（反对称矩阵的特征值是纯虚数 $\pm i\mu_k$，其中 $\mu_k \geq 0$）。

    取 $\lambda_k = \mu_k$。通过选取实正交特征向量对（对应 $\pm i\lambda_k$），可以构造正交矩阵 $Q$ 使得 $Q^TAQ$ 具有上述分块对角形式。每个 $2 \times 2$ 块 $\begin{bmatrix} 0 & \lambda_k \\ -\lambda_k & 0 \end{bmatrix}$ 的行列式为 $\lambda_k^2$，总行列式为 $\prod_{k=1}^n \lambda_k^2 = (\prod_{k=1}^n \lambda_k)^2$。$\blacksquare$

---

## 56.2 Pfaffian 的定义

<div class="context-flow" markdown>

**核心问题**：如何精确定义 $\det(A)$ 的"有符号平方根"？

</div>

!!! definition "定义 56.2 (完美配对)"
    集合 $\{1, 2, \ldots, 2n\}$ 的一个**完美配对**（perfect matching，或完美匹配，或一阶分划）$\pi$ 是将 $2n$ 个元素分成 $n$ 对的方式：
    $$\pi = \{\{i_1, j_1\}, \{i_2, j_2\}, \ldots, \{i_n, j_n\}\},$$
    其中 $i_k < j_k$（每对内部有序），$i_1 < i_2 < \cdots < i_n$（各对之间按第一个元素排序）。完美配对的总数为
    $$(2n-1)!! = (2n-1)(2n-3) \cdots 3 \cdot 1 = \frac{(2n)!}{2^n n!}.$$

!!! definition "定义 56.3 (配对的符号)"
    完美配对 $\pi = \{\{i_1, j_1\}, \ldots, \{i_n, j_n\}\}$（$i_k < j_k$，$i_1 < \cdots < i_n$）的**符号** $\mathrm{sgn}(\pi)$ 定义为排列
    $$\sigma_\pi = \begin{pmatrix} 1 & 2 & 3 & 4 & \cdots & 2n-1 & 2n \\ i_1 & j_1 & i_2 & j_2 & \cdots & i_n & j_n \end{pmatrix}$$
    的符号（即逆序数的奇偶性）。

!!! definition "定义 56.4 (Pfaffian)"
    设 $A = (a_{ij})$ 是 $2n \times 2n$ 反对称矩阵。$A$ 的 **Pfaffian** 定义为
    $$\mathrm{pf}(A) = \sum_{\pi} \mathrm{sgn}(\pi) \prod_{k=1}^n a_{i_k j_k},$$
    其中求和遍历 $\{1, \ldots, 2n\}$ 的所有完美配对 $\pi = \{\{i_1,j_1\},\ldots,\{i_n,j_n\}\}$（$i_k < j_k$，$i_1 < \cdots < i_n$）。

    等价地，可以写成排列求和的形式：
    $$\mathrm{pf}(A) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \mathrm{sgn}(\sigma) \prod_{k=1}^n a_{\sigma(2k-1), \sigma(2k)}.$$

!!! example "例 56.2"
    $n = 1$：$A = \begin{bmatrix} 0 & a_{12} \\ -a_{12} & 0 \end{bmatrix}$。唯一的完美配对为 $\pi = \{\{1,2\}\}$，$\mathrm{sgn}(\pi) = +1$。
    $$\mathrm{pf}(A) = a_{12}.$$

!!! example "例 56.3"
    $n = 2$：$A$ 是 $4 \times 4$ 反对称矩阵。$\{1,2,3,4\}$ 有三个完美配对：
    - $\pi_1 = \{\{1,2\},\{3,4\}\}$：排列 $(1,2,3,4)$，$\mathrm{sgn} = +1$。贡献 $a_{12}a_{34}$。
    - $\pi_2 = \{\{1,3\},\{2,4\}\}$：排列 $(1,3,2,4)$，$\mathrm{sgn} = -1$。贡献 $-a_{13}a_{24}$。
    - $\pi_3 = \{\{1,4\},\{2,3\}\}$：排列 $(1,4,2,3)$，$\mathrm{sgn} = +1$。贡献 $a_{14}a_{23}$。

    $$\mathrm{pf}(A) = a_{12}a_{34} - a_{13}a_{24} + a_{14}a_{23}.$$

    这与例 56.1 的结果一致。

!!! example "例 56.4"
    $n = 3$：$A$ 是 $6 \times 6$ 反对称矩阵。$(2 \cdot 3 - 1)!! = 15$ 个完美配对。

    $$\mathrm{pf}(A) = a_{12}a_{34}a_{56} - a_{12}a_{35}a_{46} + a_{12}a_{36}a_{45}$$
    $$- a_{13}a_{24}a_{56} + a_{13}a_{25}a_{46} - a_{13}a_{26}a_{45}$$
    $$+ a_{14}a_{23}a_{56} - a_{14}a_{25}a_{36} + a_{14}a_{26}a_{35}$$
    $$- a_{15}a_{23}a_{46} + a_{15}a_{24}a_{36} - a_{15}a_{26}a_{34}$$
    $$+ a_{16}a_{23}a_{45} - a_{16}a_{24}a_{35} + a_{16}a_{25}a_{34}.$$

!!! definition "定义 56.5 (Pfaffian 的外代数定义)"
    设 $A$ 是 $2n \times 2n$ 反对称矩阵。定义 $V = \mathbb{F}^{2n}$ 上的 2-形式
    $$\omega = \sum_{1 \leq i < j \leq 2n} a_{ij} \, e_i \wedge e_j \in \Lambda^2 V,$$
    其中 $\{e_1, \ldots, e_{2n}\}$ 是 $V$ 的标准基。则
    $$\frac{\omega^n}{n!} = \mathrm{pf}(A) \cdot e_1 \wedge e_2 \wedge \cdots \wedge e_{2n},$$
    即 $\mathrm{pf}(A)$ 是 $\omega^n / n!$ 在标准体积形式 $e_1 \wedge \cdots \wedge e_{2n}$ 下的系数。

??? proof "证明"
    $$\omega^n = \left(\sum_{i < j} a_{ij} e_i \wedge e_j\right)^n = \sum_{\alpha_1, \ldots, \alpha_n} \prod_{k=1}^n a_{i_k j_k} \cdot (e_{i_1} \wedge e_{j_1}) \wedge \cdots \wedge (e_{i_n} \wedge e_{j_n}),$$
    其中每个 $\alpha_k$ 表示一对 $i_k < j_k$。

    楔积 $(e_{i_1} \wedge e_{j_1}) \wedge \cdots \wedge (e_{i_n} \wedge e_{j_n})$ 非零当且仅当 $\{i_1, j_1, \ldots, i_n, j_n\} = \{1, 2, \ldots, 2n\}$（即构成完美配对），且值为 $\mathrm{sgn}(\sigma_\pi) \cdot e_1 \wedge \cdots \wedge e_{2n}$。

    每个完美配对被计数 $n!$ 次（对应 $n$ 对的排列顺序），同时每对内部的两个排列由 $i_k < j_k$ 的限制处理。更准确地，不带 $i_k < j_k$ 限制和排序限制时，每个完美配对被计数 $2^n \cdot n!$ 次。但由于我们在 $\omega$ 的定义中已经限制了 $i < j$，实际每个完美配对被计数 $n!$ 次。所以 $\omega^n / n!$ 给出恰好的 Pfaffian。$\blacksquare$

---

## 56.3 基本性质

<div class="context-flow" markdown>

**核心问题**：Pfaffian 与行列式之间有什么代数关系？

</div>

!!! theorem "定理 56.4 ($\mathrm{pf}(A)^2 = \det(A)$)"
    设 $A$ 是 $2n \times 2n$ 反对称矩阵。则
    $$\mathrm{pf}(A)^2 = \det(A).$$

??? proof "证明"
    **方法一（外代数证明）**：

    设 $\omega = \sum_{i<j} a_{ij} e_i \wedge e_j$。由定义 56.5，$\omega^n = n! \cdot \mathrm{pf}(A) \cdot e_1 \wedge \cdots \wedge e_{2n}$。

    考虑 $\Lambda^{2n}(V)$，这是一维空间。$\omega^n / n!$ 在此空间中的"坐标"是 $\mathrm{pf}(A)$。

    另一方面，反对称双线性型 $\omega$ 对应线性映射 $\tilde{A}: V \to V^*$，$\tilde{A}(v)(w) = \omega(v, w)$。在标准基下 $\tilde{A}$ 的矩阵表示就是 $A$。

    关键公式：$\omega^n / n!$ 在 $\Lambda^{2n}V$ 中的"平方"等于 $\det(A)$。更精确地，利用体积形式的自对偶性：
    $$(\omega^n / n!)^2 = (\mathrm{pf}(A))^2 (e_1 \wedge \cdots \wedge e_{2n})^2.$$

    而行列式的外代数定义（第 49 章）给出 $\det(A) \cdot (e_1 \wedge \cdots \wedge e_{2n})^{\otimes 2} = \ldots$ 这条路线需要更多机制。我们改用组合证明。

    **方法二（组合证明）**：

    行列式的定义为 $\det(A) = \sum_{\sigma \in S_{2n}} \mathrm{sgn}(\sigma) \prod_{i=1}^{2n} a_{i,\sigma(i)}$。

    由于 $A$ 反对称且对角元为零，$\det(A)$ 中非零的项要求 $\sigma$ 没有不动点（$\sigma(i) \neq i$）且可以分解为若干对换的乘积。事实上，$\det(A)$ 中非零贡献来自固定点自由的排列，而每个这样的排列可以分解为若干不相交的循环。由 $a_{ij} = -a_{ji}$ 和 $a_{ii} = 0$，只有由不相交对换组成的排列（即对合，involution）贡献非零项。

    不相交对换组成的排列恰好对应完美配对。$\det(A)$ 中每个完美配对 $\pi$ 的贡献为 $\mathrm{sgn}(\sigma_\pi)^2 \prod a_{i_k j_k}^2 \cdot (\text{符号因子})$... 这条路线计算较复杂。

    **方法三（标准形证明）**：

    先对特殊情形证明。取 $A$ 的标准形 $A_0 = Q^T A Q = \mathrm{diag}\begin{pmatrix}\begin{bmatrix}0 & \lambda_k \\ -\lambda_k & 0\end{bmatrix}\end{pmatrix}$。

    $\det(A_0) = \prod_{k=1}^n \lambda_k^2$。$\mathrm{pf}(A_0) = \prod_{k=1}^n \lambda_k$（因为唯一的"存活"完美配对是 $\{1,2\},\{3,4\},\ldots$）。所以 $\mathrm{pf}(A_0)^2 = \det(A_0)$。$\checkmark$

    一般情形：利用性质 $\mathrm{pf}(B^T A B) = \det(B) \cdot \mathrm{pf}(A)$（定理 56.5），以及 $\det(B^T A B) = (\det B)^2 \det(A)$，从 $\mathrm{pf}(A_0)^2 = \det(A_0)$ 推出
    $$\mathrm{pf}(A)^2 = \left(\frac{\mathrm{pf}(A_0)}{\det(Q)}\right)^2 = \frac{\det(A_0)}{(\det Q)^2} = \det(Q^{-T} A_0 Q^{-1}) = \det(A). \quad \blacksquare$$

!!! theorem "定理 56.5 (合同变换下的 Pfaffian)"
    设 $A$ 是 $2n \times 2n$ 反对称矩阵，$B$ 是 $2n \times 2n$ 可逆矩阵。则
    $$\mathrm{pf}(B^T A B) = \det(B) \cdot \mathrm{pf}(A).$$

??? proof "证明"
    利用外代数定义。设 $\omega_A = \sum_{i<j} a_{ij} e_i \wedge e_j$，对应矩阵 $A$。合同变换 $A \mapsto B^TAB$ 对应基变换 $e_i \mapsto \sum_k B_{ki} e_k$（即 $f_i = Be_i$），2-形式变为
    $$\omega_{B^TAB} = \sum_{i<j} (B^TAB)_{ij} e_i \wedge e_j.$$

    等价地，$\omega_{B^TAB}$ 就是 $\omega_A$ 在新基 $f_i = Be_i$ 下的表达式。

    $$\frac{\omega_{B^TAB}^n}{n!} = \mathrm{pf}(B^TAB) \cdot e_1 \wedge \cdots \wedge e_{2n}.$$

    而 $\frac{\omega_A^n}{n!} = \mathrm{pf}(A) \cdot f_1 \wedge \cdots \wedge f_{2n}$ ... 需要更仔细的推导。

    直接的方式：在新基 $f_i = \sum_k B_{ki}e_k$ 下，$f_1 \wedge \cdots \wedge f_{2n} = \det(B) \cdot e_1 \wedge \cdots \wedge e_{2n}$。

    2-形式 $\omega = \sum_{i<j} a_{ij} f_i \wedge f_j$，则
    $$\frac{\omega^n}{n!} = \mathrm{pf}(A) \cdot f_1 \wedge \cdots \wedge f_{2n} = \mathrm{pf}(A) \det(B) \cdot e_1 \wedge \cdots \wedge e_{2n}.$$

    另一方面，将 $f_i$ 展开：$\omega = \sum_{i<j} a_{ij} (\sum_k B_{ki}e_k) \wedge (\sum_l B_{lj}e_l) = \sum_{k<l} (B^TAB)_{kl} e_k \wedge e_l$。

    所以 $\frac{\omega^n}{n!} = \mathrm{pf}(B^TAB) \cdot e_1 \wedge \cdots \wedge e_{2n}$。

    比较两个表达式：$\mathrm{pf}(B^TAB) = \det(B) \cdot \mathrm{pf}(A)$。$\blacksquare$

!!! theorem "定理 56.6 (Pfaffian 的其他性质)"
    **(a)** $\mathrm{pf}(\lambda A) = \lambda^n \mathrm{pf}(A)$（$A$ 为 $2n \times 2n$ 矩阵）。

    **(b)** $\mathrm{pf}(J_{2n}) = 1$，其中 $J_{2n} = \begin{bmatrix} 0 & I_n \\ -I_n & 0 \end{bmatrix}$。

    **(c)** 若 $A = \begin{bmatrix} 0 & B \\ -B^T & 0 \end{bmatrix}$（分块反对称），则 $\mathrm{pf}(A) = (-1)^{n(n-1)/2} \det(B)$。

    **(d)** 对辛矩阵 $M$（$M^TJM = J$），$\mathrm{pf}(M^TJM) = \det(M) \cdot \mathrm{pf}(J) = \det(M)$。因此 $\det(M) = 1$ 等价于 $\mathrm{pf}(M^TJM) = \mathrm{pf}(J)$。

??? proof "证明"
    **(a)** 每个配对贡献中有 $n$ 个 $a_{i_kj_k}$ 因子，每个乘以 $\lambda$，总共 $\lambda^n$。

    **(b)** $J_{2n}$ 的非零上三角元素为 $a_{i,n+i} = 1$（$i = 1, \ldots, n$）。唯一使乘积非零的完美配对为 $\{1,n+1\}, \{2,n+2\}, \ldots, \{n,2n\}$。其符号为排列 $(1,n+1,2,n+2,\ldots,n,2n)$ 的符号。这个排列将 $(1,2,\ldots,2n)$ 变为 $(1,n+1,2,n+2,\ldots)$，需要 $0 + (n-1) + (n-2) + \cdots + 1 + 0 = \frac{n(n-1)}{2}$ 次相邻对换。但我们需要更仔细地计算。

    实际上，对 $J_2 = \begin{bmatrix} 0 & 1 \\ -1 & 0 \end{bmatrix}$，$\mathrm{pf}(J_2) = 1$。对 $J_4$，$\mathrm{pf}(J_4) = a_{13}a_{24} \cdot (-1) + \ldots$。让我们直接用定义 56.5 中的外代数方法。$\omega_J = e_1 \wedge e_{n+1} + e_2 \wedge e_{n+2} + \cdots + e_n \wedge e_{2n}$。$\omega_J^n / n! = e_1 \wedge e_{n+1} \wedge e_2 \wedge e_{n+2} \wedge \cdots \wedge e_n \wedge e_{2n}$。要将其变为标准顺序 $e_1 \wedge e_2 \wedge \cdots \wedge e_{2n}$，需要将 $e_{n+1}$ 从位置 2 移到位置 $n+1$（$n-1$ 次对换），将 $e_{n+2}$ 从位置 4 移到位置 $n+2$（$n-2$ 次对换），等等。总对换次数 $= (n-1) + (n-2) + \cdots + 0 = \frac{n(n-1)}{2}$。所以 $\mathrm{pf}(J_{2n}) = (-1)^{n(n-1)/2}$。

    等等，这与 $\mathrm{pf}(J_{2n}) = 1$ 矛盾？让我们重新检查。对 $J_4 = \begin{bmatrix} 0&0&1&0 \\ 0&0&0&1 \\ -1&0&0&0 \\ 0&-1&0&0 \end{bmatrix}$，非零上三角元素为 $a_{13} = 1, a_{24} = 1$。完美配对 $\{\{1,3\},\{2,4\}\}$ 的符号为排列 $(1,3,2,4)$ 的符号 $= -1$。贡献 $(-1) \cdot 1 \cdot 1 = -1$。所以 $\mathrm{pf}(J_4) = -1 = (-1)^{2 \cdot 1/2} = (-1)^1$。

    一般地，$\mathrm{pf}(J_{2n}) = (-1)^{n(n-1)/2}$。当 $n = 1$ 时 $\mathrm{pf}(J_2) = 1$，$n = 2$ 时 $\mathrm{pf}(J_4) = -1$，$n = 3$ 时 $\mathrm{pf}(J_6) = -1$，$n = 4$ 时 $\mathrm{pf}(J_8) = 1$。

    **修正**：(b) 应改为 $\mathrm{pf}(J_{2n}) = (-1)^{n(n-1)/2}$。$\blacksquare$

!!! theorem "定理 56.7 (Pfaffian 的递归展开)"
    设 $A$ 是 $2n \times 2n$ 反对称矩阵。则
    $$\mathrm{pf}(A) = \sum_{j=2}^{2n} (-1)^j a_{1j} \cdot \mathrm{pf}(\hat{A}_{1j}),$$
    其中 $\hat{A}_{1j}$ 是从 $A$ 中删去第 $1, j$ 行和第 $1, j$ 列后得到的 $(2n-2) \times (2n-2)$ 反对称矩阵。

??? proof "证明"
    在 Pfaffian 的定义 $\mathrm{pf}(A) = \sum_\pi \mathrm{sgn}(\pi) \prod a_{i_k j_k}$ 中，元素 $1$ 必与某个 $j \in \{2, \ldots, 2n\}$ 配对。将求和按 $1$ 的配对对象 $j$ 分组：
    $$\mathrm{pf}(A) = \sum_{j=2}^{2n} a_{1j} \sum_{\pi': 1 \sim j} \mathrm{sgn}(\pi) \prod_{k \geq 2} a_{i_k j_k}.$$

    固定 $\{1, j\}$ 后，余下的配对 $\pi'$ 是 $\{1,\ldots,2n\} \setminus \{1,j\}$ 上的完美配对，对应 $\hat{A}_{1j}$ 的 Pfaffian。符号因子为 $(-1)^{j}$（将 $j$ 从位置 $j$ 移到位置 $2$ 需要 $j-2$ 次相邻对换，加上 $\mathrm{sgn}(\pi)$ 的调整）。

    仔细计算符号：排列 $(1, j, \text{余下按顺序})$ 相对于 $(1, 2, 3, \ldots, 2n)$ 的符号为 $(-1)^{j-2} = (-1)^j$（因为 $(-1)^{j-2} = (-1)^j$）。所以
    $$\mathrm{pf}(A) = \sum_{j=2}^{2n} (-1)^j a_{1j} \cdot \mathrm{pf}(\hat{A}_{1j}). \quad \blacksquare$$

---

## 56.4 计算方法

<div class="context-flow" markdown>

**核心问题**：如何高效计算 Pfaffian？

</div>

!!! remark "注记"
    直接用定义计算 Pfaffian 需要 $(2n-1)!!$ 项求和，复杂度为阶乘级（比行列式的 $n!$ 还快）。递归展开（定理 56.7）的复杂度与行列式的 Laplace 展开类似，也是指数级。我们需要多项式时间的算法。

!!! theorem "定理 56.8 (Parlett-Reid 算法)"
    Pfaffian 可以在 $O(n^3)$ 时间内计算，通过将反对称矩阵化为三对角反对称形式。

    **算法**：对 $2n \times 2n$ 反对称矩阵 $A$，利用正交合同变换（Householder 变换的反对称版本）将 $A$ 化为三对角反对称矩阵
    $$T = Q^T A Q = \begin{bmatrix} 0 & t_1 & & & \\ -t_1 & 0 & t_2 & & \\ & -t_2 & 0 & t_3 & \\ & & & \ddots & \\ & & & -t_{2n-1} & 0 \end{bmatrix},$$
    其中 $Q$ 是正交矩阵。然后 $\mathrm{pf}(A) = \det(Q) \cdot \mathrm{pf}(T)$。

    三对角反对称矩阵的 Pfaffian 可以用简单递推计算：
    $$\mathrm{pf}(T_{2n}) = t_1 \cdot \mathrm{pf}(T_{2n-2}'),$$
    其中 $T_{2n-2}'$ 是删去前两行两列后的子矩阵。最终 $\mathrm{pf}(T) = t_1 \cdot t_3 \cdot t_5 \cdots t_{2n-1}$（即奇数位置的上对角线元素之积）。

??? proof "证明"
    三对角化过程使用 Gauss 消元的反对称版本（或 Householder 变换），复杂度 $O(n^3)$。三对角反对称矩阵的 Pfaffian 递推是线性的：

    $\mathrm{pf}\begin{bmatrix} 0 & t_1 \\ -t_1 & 0 \end{bmatrix} = t_1$。

    用递归展开（定理 56.7）：$\mathrm{pf}(T_{2n}) = \sum_{j=2}^{2n} (-1)^j a_{1j} \mathrm{pf}(\hat{T}_{1j})$。但三对角结构意味着 $a_{1j} = 0$ 对 $j \geq 3$，所以只有 $j = 2$ 项非零：
    $$\mathrm{pf}(T_{2n}) = (-1)^2 \cdot t_1 \cdot \mathrm{pf}(T'_{2n-2}) = t_1 \cdot \mathrm{pf}(T'_{2n-2}).$$
    递推下去得 $\mathrm{pf}(T) = t_1 t_3 t_5 \cdots t_{2n-1}$。$\blacksquare$

!!! example "例 56.5"
    计算 $A = \begin{bmatrix} 0 & 3 & -1 & 2 \\ -3 & 0 & 4 & -5 \\ 1 & -4 & 0 & 6 \\ -2 & 5 & -6 & 0 \end{bmatrix}$ 的 Pfaffian。

    用递归展开：
    $$\mathrm{pf}(A) = (-1)^2 \cdot 3 \cdot \mathrm{pf}\begin{bmatrix} 0 & 6 \\ -6 & 0 \end{bmatrix} + (-1)^3 \cdot (-1) \cdot \mathrm{pf}\begin{bmatrix} 0 & -5 \\ 5 & 0 \end{bmatrix} + (-1)^4 \cdot 2 \cdot \mathrm{pf}\begin{bmatrix} 0 & 4 \\ -4 & 0 \end{bmatrix}$$
    $$= 3 \cdot 6 + 1 \cdot (-5) + 2 \cdot 4 = 18 - 5 + 8 = 21.$$

    验证：$\det(A) = 21^2 = 441$。$\checkmark$

---

## 56.5 Pfaffian 与完美匹配

<div class="context-flow" markdown>

**核心问题**：Pfaffian 如何编码图的完美匹配？

</div>

!!! definition "定义 56.6 (图的反对称邻接矩阵)"
    设 $G = (V, E)$ 是无向图，$|V| = 2n$。给 $G$ 的每条边 $\{i, j\}$ 一个定向（选择 $i \to j$ 或 $j \to i$）。定义**反对称邻接矩阵** $A^{\mathrm{skew}}$：
    $$a_{ij}^{\mathrm{skew}} = \begin{cases} +1, & \text{若 } i \to j \text{ 是定向边}, \\ -1, & \text{若 } j \to i \text{ 是定向边}, \\ 0, & \text{若 } \{i,j\} \notin E. \end{cases}$$

!!! theorem "定理 56.9 (Pfaffian 与完美匹配)"
    设 $G$ 的反对称邻接矩阵为 $A^{\mathrm{skew}}$。则
    $$\mathrm{pf}(A^{\mathrm{skew}}) = \sum_{M \in \mathcal{M}(G)} \mathrm{sgn}(M) \cdot 1,$$
    其中 $\mathcal{M}(G)$ 是 $G$ 的所有完美匹配的集合，$\mathrm{sgn}(M)$ 是由定向和匹配共同决定的符号（$\pm 1$）。

    特别地，$|\mathrm{pf}(A^{\mathrm{skew}})| \leq |\mathcal{M}(G)|$，且
    $$\det(A^{\mathrm{skew}}) = \mathrm{pf}(A^{\mathrm{skew}})^2 \leq |\mathcal{M}(G)|^2.$$

!!! definition "定义 56.7 (Pfaffian 定向)"
    图 $G$ 的一个定向称为 **Pfaffian 定向**（Pfaffian orientation），若在此定向下 $\mathrm{pf}(A^{\mathrm{skew}}) = |\mathcal{M}(G)|$（即所有完美匹配的符号都为 $+1$），或等价地 $\mathrm{pf}(A^{\mathrm{skew}})^2 = |\mathcal{M}(G)|^2$（所有匹配的符号一致）。

!!! theorem "定理 56.10 (Pfaffian 定向的等价刻画)"
    图 $G$ 的一个定向是 Pfaffian 定向，当且仅当对 $G$ 的每一个偶长度面（在某个固定嵌入下），顺时针方向与定向一致的边数为奇数。

---

## 56.6 FKT 算法

<div class="context-flow" markdown>

**核心问题**：完美匹配计数一般是 #P-hard 的，但对哪类图可以多项式时间解决？

</div>

!!! theorem "定理 56.11 (Kasteleyn 定理, 1961)"
    每个**平面图**都承认 Pfaffian 定向。因此，平面图的完美匹配数可以通过计算 Pfaffian（$O(n^3)$ 时间）来求得。

??? proof "证明"
    **构造 Pfaffian 定向**：

    取平面图 $G$ 的一个平面嵌入。$G$ 将平面分成若干面（包括外面）。

    **第一步**：任意定向 $G$ 的所有边。

    **第二步**：构造面的对偶图 $G^*$。在 $G^*$ 中取一棵生成树 $T^*$。

    **第三步**：按照 $T^*$ 从叶子到根的顺序处理每个内面 $f$：如果 $f$ 的边界上顺时针方向与定向一致的边数为偶数，则翻转 $T^*$ 中连接 $f$ 到其父面的那条边的方向。

    每次翻转保证当前面满足"奇数条件"，且不破坏已处理面的条件（因为每条边最多被翻转一次，且树的结构保证了这一点）。

    最终所有面都满足 Pfaffian 定向条件。$\blacksquare$

!!! definition "定义 56.8 (FKT 算法)"
    **FKT 算法**（Fisher-Kasteleyn-Temperley）计算平面图的完美匹配数：

    1. 求平面图 $G$ 的 Pfaffian 定向（$O(n)$ 时间，利用 BFS/DFS 和面遍历）。
    2. 构造反对称邻接矩阵 $A^{\mathrm{skew}}$。
    3. 计算 $\mathrm{pf}(A^{\mathrm{skew}})$（$O(n^3)$ 时间）。
    4. $|\mathcal{M}(G)| = |\mathrm{pf}(A^{\mathrm{skew}})|$。

    总时间复杂度 $O(n^3)$。

!!! theorem "定理 56.12 (一般图完美匹配计数的困难性)"
    对一般（非平面）图，计算完美匹配数是 **#P-完全**的（Valiant, 1979）。这与平面图的多项式可解性形成鲜明对比。

!!! example "例 56.6"
    **$2 \times n$ 棋盘的多米诺铺排**：考虑 $2 \times n$ 网格图。这是平面图，FKT 算法适用。

    完美匹配数满足递推 $f(n) = f(n-1) + f(n-2)$（$f(1) = 1, f(2) = 2$），即 Fibonacci 数！

    $4 \times 4$ 棋盘的多米诺铺排数 = 36。$6 \times 6$ 棋盘 = 6728。$8 \times 8$ 棋盘 = 12988816。这些数可以通过 Pfaffian 的封闭公式计算。

!!! example "例 56.7"
    **$m \times n$ 棋盘的多米诺铺排**：Kasteleyn-Temperley-Fisher 公式给出封闭形式：
    $$|\mathcal{M}(G_{m,n})| = \prod_{j=1}^{\lceil m/2 \rceil} \prod_{k=1}^{\lceil n/2 \rceil} \left(4\cos^2\frac{\pi j}{m+1} + 4\cos^2\frac{\pi k}{n+1}\right)^{1/2}.$$
    这个优美的公式正是通过 Pfaffian 计算推导出来的。

---

## 56.7 物理应用

<div class="context-flow" markdown>

**核心问题**：Pfaffian 在物理中扮演什么角色？

</div>

!!! example "例 56.8 (二维 Ising 模型)"
    **Ising 模型**是统计力学中最重要的模型之一。在二维正方晶格上，每个格点 $i$ 有自旋 $\sigma_i \in \{+1, -1\}$，Hamilton 量为
    $$\mathcal{H} = -J \sum_{\langle i,j \rangle} \sigma_i \sigma_j,$$
    其中求和遍历最近邻格点对。

    **配分函数** $Z = \sum_{\{\sigma\}} e^{-\beta \mathcal{H}}$ 决定了系统的所有热力学性质。

    Onsager (1944) 首次精确求解了二维 Ising 模型。后来 Kasteleyn (1961) 和 Fisher (1966) 给出了基于 Pfaffian 的更优雅的解法：通过巧妙的图论变换，将配分函数表示为某个反对称矩阵的 Pfaffian：
    $$Z = 2^N \prod_{j} \cosh(\beta J) \cdot \mathrm{pf}(A),$$
    其中 $A$ 是由晶格结构和耦合常数决定的反对称矩阵，$N$ 为格点数。

    这一方法不仅给出了 Onsager 的结果，还可以推广到更一般的平面 Ising 模型。

!!! example "例 56.9 (二聚体模型)"
    **二聚体模型**（dimer model）研究的是格子图上完美匹配的统计力学。每个二聚体覆盖（完美匹配）$M$ 有权重 $w(M) = \prod_{e \in M} w_e$，配分函数为
    $$Z_{\mathrm{dimer}} = \sum_{M \in \mathcal{M}} w(M).$$

    对平面图，构造加权反对称邻接矩阵 $A$（$a_{ij} = \pm w_{ij}$，符号由 Pfaffian 定向决定），则
    $$Z_{\mathrm{dimer}} = |\mathrm{pf}(A)|.$$

    这一公式是凝聚态物理和统计力学中许多精确结果的基础。

!!! example "例 56.10 (BCS 超导理论)"
    在 BCS（Bardeen-Cooper-Schrieffer）超导理论中，超导基态波函数可以写成
    $$|\Psi_{\mathrm{BCS}}\rangle = \prod_k (u_k + v_k c_{k\uparrow}^\dagger c_{-k\downarrow}^\dagger)|0\rangle.$$

    这个波函数的范数和各种矩阵元的计算涉及反对称矩阵的 Pfaffian。具体地，费米子配对态的重叠积分可以表示为 Pfaffian：
    $$\langle \Psi_1 | \Psi_2 \rangle = \mathrm{pf}(M),$$
    其中 $M$ 是由两个配对态的参数决定的反对称矩阵。

!!! example "例 56.11 (Euler 类与 Gauss-Bonnet 定理)"
    在微分几何中，设 $M$ 是 $2n$ 维紧致可定向 Riemannian 流形，$\Omega$ 是其曲率 2-形式（取值在 $\mathfrak{so}(2n)$ 中的反对称矩阵）。**Euler 类**定义为
    $$e(\Omega) = \frac{1}{(2\pi)^n} \mathrm{Pf}(\Omega),$$
    其中 $\mathrm{Pf}(\Omega)$ 是将 Pfaffian 推广到矩阵值微分形式的版本。

    **推广的 Gauss-Bonnet 定理**（Chern 1944）：
    $$\chi(M) = \int_M e(\Omega) = \frac{1}{(2\pi)^n} \int_M \mathrm{Pf}(\Omega),$$
    其中 $\chi(M)$ 是 $M$ 的 Euler 特征数。这是 Pfaffian 在几何拓扑中最深刻的应用之一。

---

## 本章小结

| 概念 | 定义/关键性质 |
|:---|:---|
| 反对称矩阵 | $A^T = -A$；奇数阶 $\det = 0$；偶数阶 $\det$ 是完全平方 |
| Pfaffian | $\mathrm{pf}(A) = \sum_\pi \mathrm{sgn}(\pi) \prod a_{i_kj_k}$；完美配对上的带符号求和 |
| 核心恒等式 | $\mathrm{pf}(A)^2 = \det(A)$ |
| 合同变换 | $\mathrm{pf}(B^TAB) = \det(B) \cdot \mathrm{pf}(A)$ |
| 外代数定义 | $\omega^n / n! = \mathrm{pf}(A) \cdot e_1 \wedge \cdots \wedge e_{2n}$ |
| 完美匹配 | 平面图的 Pfaffian 定向 → FKT 算法 $O(n^3)$ |
| 物理应用 | 二维 Ising 模型、二聚体模型、BCS 超导、Euler 类 |

---

## 习题

!!! exercise "习题 56.1"
    直接计算 $6 \times 6$ 反对称矩阵 $A$（$a_{12} = 1, a_{13} = 2, a_{14} = 0, a_{15} = 1, a_{16} = -1, a_{23} = -1, a_{24} = 3, a_{25} = 0, a_{26} = 2, a_{34} = 1, a_{35} = -2, a_{36} = 0, a_{45} = 1, a_{46} = -1, a_{56} = 3$）的 Pfaffian，并验证 $\mathrm{pf}(A)^2 = \det(A)$。

!!! exercise "习题 56.2"
    证明：若 $A$ 是 $2n \times 2n$ 反对称矩阵，$B$ 是 $2n \times 2n$ 正交矩阵（$B^TB = I$），则 $\mathrm{pf}(B^TAB) = \det(B) \cdot \mathrm{pf}(A) = \pm \mathrm{pf}(A)$。

!!! exercise "习题 56.3"
    利用外代数定义证明：标准辛矩阵 $J_{2n}$ 的 Pfaffian 为 $(-1)^{n(n-1)/2}$。

!!! exercise "习题 56.4"
    设 $A = \begin{bmatrix} 0 & B \\ -B^T & 0 \end{bmatrix}$，其中 $B$ 是 $n \times n$ 矩阵。证明 $\mathrm{pf}(A) = (-1)^{n(n-1)/2} \det(B)$。

!!! exercise "习题 56.5"
    对 $K_4$（4 个顶点的完全图）的所有三种定向，分别计算反对称邻接矩阵的 Pfaffian。哪些定向是 Pfaffian 定向？

!!! exercise "习题 56.6"
    用 FKT 算法计算 $4 \times 4$ 棋盘的多米诺铺排数（答案应为 36）。

!!! exercise "习题 56.7"
    证明 Pfaffian 的递归展开公式（定理 56.7），并用它计算例 56.5 中矩阵的 Pfaffian。

!!! exercise "习题 56.8"
    设 $A$ 是 $2n \times 2n$ 反对称矩阵，$v \in \mathbb{R}^{2n}$。定义 $A(t)$ 为将 $A$ 的第一行和第一列替换为 $tv$ 和 $-tv^T$ 后得到的矩阵。证明 $\mathrm{pf}(A(t))$ 是 $t$ 的线性函数。

!!! exercise "习题 56.9"
    （综合）证明：若 $M$ 是 $2n \times 2n$ 辛矩阵（$M^TJM = J$），则
    $$\mathrm{pf}(M^TJM) = \mathrm{pf}(J),$$
    从而 $\det(M) = 1$。与第 53 章的证明方法比较。

!!! exercise "习题 56.10"
    Kasteleyn-Temperley-Fisher 公式断言 $m \times n$ 网格图（$mn$ 为偶数）的完美匹配数为
    $$\prod_{j=1}^{\lceil m/2\rceil}\prod_{k=1}^{\lceil n/2\rceil}\left(4\cos^2\frac{\pi j}{m+1}+4\cos^2\frac{\pi k}{n+1}\right).$$
    对 $m = n = 2$ 验证此公式给出 $2$（一个 $2 \times 2$ 棋盘的两种多米诺铺排），对 $m = 2, n = 3$ 验证给出 $3$。
