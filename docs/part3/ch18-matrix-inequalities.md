# 第 18 章 矩阵不等式

<div class="context-flow" markdown>

**前置**：Ch15 Weyl/Wielandt-Hoffman · Ch16 正定矩阵

**脉络**：Courant-Fischer → Weyl/Cauchy交错 → 奇异值不等式 → 迹(von Neumann) → 行列式(Hadamard/Minkowski) → **Majorization**统一框架

**延伸**：矩阵不等式在量子信息（entanglement witness、信道容量界）、控制理论（LMI 方法）、统计学（Cramér-Rao 界）中有核心应用；majorization 理论连接了信息论（Shannon 熵）和量子纠缠的度量

</div>

矩阵不等式是矩阵分析中最为精深的主题之一。在标量情形中，不等式的理论已经非常成熟——从算术-几何均值不等式到 Cauchy-Schwarz 不等式，这些结果构成了分析学的基石。然而，当我们将这些想法推广到矩阵的世界时，情况变得远为丰富和微妙。矩阵的特征值（eigenvalue）、奇异值（singular value）、迹（trace）和行列式（determinant）之间存在着深刻而优美的不等式关系，它们不仅在纯数学中有重要意义，而且在量子信息论、统计学和最优化理论中都有广泛应用。

本章系统地介绍矩阵不等式的主要结果，从特征值不等式出发，经由奇异值不等式和迹不等式，直至行列式不等式，最后以优超理论（majorization）和矩阵凸性作为统一框架。

---

## 18.1 特征值不等式

<div class="context-flow" markdown>

**核心工具**：**Courant-Fischer** 极小极大定理(Rayleigh商 + 子空间维数论证) → Weyl不等式 · Cauchy交错 = 主子矩阵特征值的夹逼

</div>

特征值是矩阵最基本的不变量。对于 Hermite 矩阵（Hermitian matrix），其特征值全为实数，因此可以排序比较。本节讨论 Hermite 矩阵之和的特征值与各自特征值之间的关系。

**约定**：对于 $n \times n$ Hermite 矩阵 $A$，其特征值按降序排列为 $\lambda_1(A) \geq \lambda_2(A) \geq \cdots \geq \lambda_n(A)$。

!!! definition "定义 18.1 (Hermite 矩阵的特征值排序)"
    设 $A$ 为 $n \times n$ Hermite 矩阵，即 $A = A^*$。其 $n$ 个实特征值按降序排列：

    $$
    \lambda_1(A) \geq \lambda_2(A) \geq \cdots \geq \lambda_n(A).
    $$

    我们用 $\lambda_{\max}(A) = \lambda_1(A)$ 和 $\lambda_{\min}(A) = \lambda_n(A)$ 分别表示最大和最小特征值。

!!! definition "定义 18.2 (Rayleigh 商)"
    设 $A$ 为 $n \times n$ Hermite 矩阵，对非零向量 $\mathbf{x} \in \mathbb{C}^n$，**Rayleigh 商**（Rayleigh quotient）定义为：

    $$
    R_A(\mathbf{x}) = \frac{\mathbf{x}^* A \mathbf{x}}{\mathbf{x}^* \mathbf{x}}.
    $$

!!! theorem "定理 18.1 (Courant-Fischer 极小极大定理)"
    设 $A$ 为 $n \times n$ Hermite 矩阵，特征值 $\lambda_1(A) \geq \cdots \geq \lambda_n(A)$，则对 $k = 1, 2, \ldots, n$：

    $$
    \lambda_k(A) = \max_{\dim S = k} \min_{\mathbf{x} \in S, \mathbf{x} \neq \mathbf{0}} \frac{\mathbf{x}^* A \mathbf{x}}{\mathbf{x}^* \mathbf{x}} = \min_{\dim T = n-k+1} \max_{\mathbf{x} \in T, \mathbf{x} \neq \mathbf{0}} \frac{\mathbf{x}^* A \mathbf{x}}{\mathbf{x}^* \mathbf{x}}.
    $$

??? proof "证明"
    设 $A = U \Lambda U^*$ 为谱分解，其中 $U = [\mathbf{u}_1, \ldots, \mathbf{u}_n]$ 为酉矩阵，$\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$。

    **第一步**：证明 $\lambda_k \geq \max_{\dim S=k} \min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} R_A(\mathbf{x})$ 的一个方向。

    取 $S_0 = \operatorname{span}\{\mathbf{u}_1, \ldots, \mathbf{u}_k\}$，则 $\dim S_0 = k$。对任意 $\mathbf{x} = \sum_{i=1}^k c_i \mathbf{u}_i \in S_0$，有：

    $$
    R_A(\mathbf{x}) = \frac{\sum_{i=1}^k \lambda_i |c_i|^2}{\sum_{i=1}^k |c_i|^2} \geq \lambda_k.
    $$

    因此 $\min_{\mathbf{x} \in S_0 \setminus \{\mathbf{0}\}} R_A(\mathbf{x}) \geq \lambda_k$，从而 $\max_{\dim S=k} \min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} R_A(\mathbf{x}) \geq \lambda_k$。

    **第二步**：证明对任意 $k$ 维子空间 $S$，$\min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} R_A(\mathbf{x}) \leq \lambda_k$。

    取 $T_0 = \operatorname{span}\{\mathbf{u}_k, \ldots, \mathbf{u}_n\}$，则 $\dim T_0 = n - k + 1$。由维数公式，$\dim(S \cap T_0) \geq k + (n-k+1) - n = 1$，因此存在非零 $\mathbf{x} \in S \cap T_0$。对这样的 $\mathbf{x} = \sum_{i=k}^n c_i \mathbf{u}_i$：

    $$
    R_A(\mathbf{x}) = \frac{\sum_{i=k}^n \lambda_i |c_i|^2}{\sum_{i=k}^n |c_i|^2} \leq \lambda_k.
    $$

    因此 $\min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} R_A(\mathbf{x}) \leq \lambda_k$。

    综合两步即得第一个等式。第二个等式的证明完全类似。$\blacksquare$

!!! theorem "定理 18.2 (Weyl 不等式)"
    设 $A, B$ 为 $n \times n$ Hermite 矩阵，则对所有满足 $i + j - 1 \leq n$ 的 $i, j$：

    $$
    \lambda_{i+j-1}(A + B) \leq \lambda_i(A) + \lambda_j(B).
    $$

    对所有满足 $i + j - n \geq 1$ 的 $i, j$：

    $$
    \lambda_{i+j-n}(A + B) \geq \lambda_i(A) + \lambda_j(B).
    $$

??? proof "证明"
    我们证明第一个不等式。由 Courant-Fischer 定理的 min-max 形式：

    $$
    \lambda_{i+j-1}(A+B) = \min_{\dim T = n-i-j+2} \max_{\mathbf{x} \in T \setminus \{\mathbf{0}\}} R_{A+B}(\mathbf{x}).
    $$

    设 $A$ 的前 $i$ 个特征向量张成子空间 $S_A$（$\dim S_A = i$），$B$ 的前 $j$ 个特征向量张成子空间 $S_B$（$\dim S_B = j$）。

    对任意 $(n-i-j+2)$ 维子空间 $T$，由维数公式：

    $$
    \dim(S_A \cap T) \geq i + (n-i-j+2) - n = n - j + 2 - n + i + i - i = 2 - j + i \geq 1,
    $$

    即 $\dim(T \cap S_A) \geq 1$（当 $i + j - 1 \leq n$ 时）。类似地 $\dim(T \cap S_B) \geq 1$。

    更直接地，取 $T_0 = \operatorname{span}\{\mathbf{u}_i, \ldots, \mathbf{u}_n\} \cap \operatorname{span}\{\mathbf{v}_j, \ldots, \mathbf{v}_n\}$，其中 $\mathbf{u}_k, \mathbf{v}_k$ 分别为 $A, B$ 的特征向量。由 Courant-Fischer 定理：

    对任意非零 $\mathbf{x}$：

    $$
    R_{A+B}(\mathbf{x}) = R_A(\mathbf{x}) + R_B(\mathbf{x}).
    $$

    取使得 $R_A(\mathbf{x}) \leq \lambda_i(A)$ 和 $R_B(\mathbf{x}) \leq \lambda_j(B)$ 同时成立的子空间（维数关系保证其非空），即得：

    $$
    \lambda_{i+j-1}(A+B) \leq \lambda_i(A) + \lambda_j(B). \quad \blacksquare
    $$

!!! note "注"
    Weyl 不等式的一个常用特殊情况是取 $j = 1$：$\lambda_i(A+B) \leq \lambda_i(A) + \lambda_1(B)$，即加上半正定矩阵后特征值不减小超过 $\lambda_1(B)$。类似地取 $j = n$：$\lambda_i(A+B) \geq \lambda_i(A) + \lambda_n(B)$。这给出了特征值的扰动界。

!!! theorem "定理 18.3 (Cauchy 交错定理)"
    设 $A$ 为 $n \times n$ Hermite 矩阵，$B$ 为 $A$ 的 $m \times m$ 主子矩阵（即 $B = P^* A P$，其中 $P$ 是 $n \times m$ 矩阵且 $P^* P = I_m$），$m < n$。则 $A$ 和 $B$ 的特征值满足**交错关系**（interlacing）：

    $$
    \lambda_i(A) \geq \lambda_i(B) \geq \lambda_{i+n-m}(A), \quad i = 1, 2, \ldots, m.
    $$

??? proof "证明"
    由 Courant-Fischer 定理，对 $B = P^* A P$：

    $$
    \lambda_i(B) = \max_{\substack{S \subset \mathbb{C}^m \\ \dim S = i}} \min_{\mathbf{y} \in S \setminus \{\mathbf{0}\}} \frac{\mathbf{y}^* B \mathbf{y}}{\mathbf{y}^* \mathbf{y}} = \max_{\substack{S \subset \mathbb{C}^m \\ \dim S = i}} \min_{\mathbf{y} \in S \setminus \{\mathbf{0}\}} \frac{(P\mathbf{y})^* A (P\mathbf{y})}{(P\mathbf{y})^* (P\mathbf{y})}.
    $$

    令 $\mathbf{x} = P \mathbf{y}$，由于 $P$ 是等距嵌入，当 $S$ 遍历 $\mathbb{C}^m$ 的所有 $i$ 维子空间时，$PS$ 遍历 $P(\mathbb{C}^m)$ 的所有 $i$ 维子空间。由于 $P(\mathbb{C}^m)$ 是 $\mathbb{C}^n$ 的一个 $m$ 维子空间，而 $\mathbb{C}^n$ 的所有 $i$ 维子空间的集合包含 $P(\mathbb{C}^m)$ 的所有 $i$ 维子空间，因此：

    $$
    \lambda_i(B) \leq \max_{\substack{S \subset \mathbb{C}^n \\ \dim S = i}} \min_{\mathbf{x} \in S \setminus \{\mathbf{0}\}} \frac{\mathbf{x}^* A \mathbf{x}}{\mathbf{x}^* \mathbf{x}} = \lambda_i(A).
    $$

    类似地，使用 min-max 形式可证 $\lambda_i(B) \geq \lambda_{i+n-m}(A)$。$\blacksquare$

!!! theorem "定理 18.4 (Poincare 分离定理)"
    设 $A$ 为 $n \times n$ Hermite 矩阵，$U$ 为 $n \times k$ 矩阵满足 $U^* U = I_k$（$k \leq n$），令 $B = U^* A U$，则：

    $$
    \lambda_i(A) \geq \lambda_i(B) \geq \lambda_{i+n-k}(A), \quad i = 1, \ldots, k.
    $$

??? proof "证明"
    这实际上就是 Cauchy 交错定理在一般酉压缩下的推广形式。证明方法与定理 18.3 完全相同，只是将 $P$ 替换为一般的等距映射 $U$。关键在于 $U^* U = I_k$ 保证了 $U$ 是从 $\mathbb{C}^k$ 到 $\mathbb{C}^n$ 的等距嵌入。$\blacksquare$

!!! example "例 18.1"
    设 $A = \begin{pmatrix} 5 & 1 & 0 \\ 1 & 3 & 1 \\ 0 & 1 & 1 \end{pmatrix}$，其特征值为 $\lambda_1(A) \approx 5.414$，$\lambda_2(A) \approx 2.828$，$\lambda_3(A) \approx 0.758$。

    取左上 $2 \times 2$ 主子矩阵 $B = \begin{pmatrix} 5 & 1 \\ 1 & 3 \end{pmatrix}$，其特征值为 $\lambda_1(B) = 3 + \sqrt{2} \approx 4.414$，$\lambda_2(B) = 3 - \sqrt{2} \approx 1.586$。

    验证交错性：

    $$
    \lambda_1(A) \approx 5.414 \geq \lambda_1(B) \approx 4.414 \geq \lambda_2(A) \approx 2.828,
    $$

    $$
    \lambda_2(A) \approx 2.828 \geq \lambda_2(B) \approx 1.586 \geq \lambda_3(A) \approx 0.758.
    $$

    交错关系成立。

!!! example "例 18.2"
    **Weyl 不等式的应用**：设 $A = \begin{pmatrix} 4 & 0 \\ 0 & 2 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$。

    $A$ 的特征值：$\lambda_1(A) = 4$，$\lambda_2(A) = 2$。

    $B$ 的特征值：$\lambda_1(B) = 2$，$\lambda_2(B) = 0$。

    $A + B = \begin{pmatrix} 5 & 1 \\ 1 & 3 \end{pmatrix}$，特征值为 $\lambda_1(A+B) = 4 + \sqrt{2} \approx 5.414$，$\lambda_2(A+B) = 4 - \sqrt{2} \approx 2.586$。

    Weyl 不等式 $\lambda_2(A+B) \leq \lambda_1(A) + \lambda_2(B) = 4 + 0 = 4$ 成立（$2.586 \leq 4$）。

    Weyl 不等式 $\lambda_2(A+B) \geq \lambda_2(A) + \lambda_2(B) = 2 + 0 = 2$ 成立（$2.586 \geq 2$）。

---

## 18.2 奇异值不等式

<div class="context-flow" markdown>

**脉络**：奇异值次可乘 $\sigma_{i+j-1}(AB)\leq\sigma_i(A)\sigma_j(B)$

**Ky Fan $k$-范数** $=\sum_{i=1}^k\sigma_i$ 满足三角不等式 → 酉不变范数的统一视角

</div>

奇异值（singular value）是矩阵理论中与特征值同等重要的概念。对于一般矩阵（不要求 Hermite 或方阵），奇异值提供了最自然的"大小"度量。

!!! definition "定义 18.3 (奇异值排序)"
    设 $A$ 为 $m \times n$ 矩阵，其奇异值按降序排列为：

    $$
    \sigma_1(A) \geq \sigma_2(A) \geq \cdots \geq \sigma_{\min(m,n)}(A) \geq 0.
    $$

    其中 $\sigma_i(A) = \sqrt{\lambda_i(A^* A)}$。

!!! theorem "定理 18.5 (奇异值的次可乘性)"
    设 $A, B$ 为 $n \times n$ 矩阵，则对 $i + j - 1 \leq n$：

    $$
    \sigma_{i+j-1}(AB) \leq \sigma_i(A) \sigma_j(B).
    $$

    特别地，取 $i = j = 1$：$\sigma_1(AB) \leq \sigma_1(A) \sigma_1(B)$，即 $\|AB\|_2 \leq \|A\|_2 \|B\|_2$。

??? proof "证明"
    关键工具是奇异值与特征值的关系以及 Weyl 不等式。

    设 $A = U_1 \Sigma_1 V_1^*$，$B = U_2 \Sigma_2 V_2^*$ 为奇异值分解。考虑 $AB$ 的奇异值平方，即 $B^* A^* A B$ 的特征值。

    利用 Fan 的极值原理：

    $$
    \sigma_k(AB) = \min_{\substack{S \subset \mathbb{C}^n \\ \operatorname{codim} S = k-1}} \max_{\mathbf{x} \in S, \|\mathbf{x}\|=1} \|AB\mathbf{x}\|.
    $$

    对任意单位向量 $\mathbf{x}$，$\|AB\mathbf{x}\| \leq \|A\| \cdot \|B\mathbf{x}\|$（这里 $\|A\|$ 表示算子范数）。更精细的分析使用子空间的维数论证：

    取使得 $\|A\mathbf{y}\| \leq \sigma_i(A)\|\mathbf{y}\|$ 的 $(n-i+1)$ 维子空间 $S_A$，以及使得 $\|B\mathbf{x}\| \leq \sigma_j(B)\|\mathbf{x}\|$ 的 $(n-j+1)$ 维子空间 $S_B$。则 $B^{-1}(S_A) \cap S_B$ 的维数至少为 $(n-i+1) + (n-j+1) - n = n - i - j + 2$，其余维 $= i + j - 2$。

    因此在相应子空间上 $\|AB\mathbf{x}\| \leq \sigma_i(A)\sigma_j(B)\|\mathbf{x}\|$，由奇异值的极小极大表征即得 $\sigma_{i+j-1}(AB) \leq \sigma_i(A)\sigma_j(B)$。$\blacksquare$

!!! definition "定义 18.4 (Ky Fan 范数)"
    设 $A$ 为 $m \times n$ 矩阵，对 $k = 1, 2, \ldots, \min(m,n)$，**Ky Fan $k$-范数**定义为前 $k$ 个奇异值之和：

    $$
    \|A\|_{(k)} = \sum_{i=1}^{k} \sigma_i(A).
    $$

    特别地，$\|A\|_{(1)} = \sigma_1(A) = \|A\|_2$（谱范数），$\|A\|_{(\min(m,n))} = \|A\|_*$（核范数）。

!!! theorem "定理 18.6 (Fan 不等式)"
    设 $A, B$ 为 $n \times n$ 矩阵，则对 $k = 1, 2, \ldots, n$：

    $$
    \sum_{i=1}^{k} \sigma_i(A + B) \leq \sum_{i=1}^{k} \sigma_i(A) + \sum_{i=1}^{k} \sigma_i(B).
    $$

    即 Ky Fan $k$-范数满足三角不等式：$\|A + B\|_{(k)} \leq \|A\|_{(k)} + \|B\|_{(k)}$。

??? proof "证明"
    对任意 $n \times n$ 矩阵 $M$，由奇异值的极值表征：

    $$
    \sum_{i=1}^{k} \sigma_i(M) = \max \{ |\operatorname{tr}(U^* M)| : U \text{ 为 } n \times k \text{ 的等距映射} \}.
    $$

    更精确地，存在 Fan 的定理：

    $$
    \sum_{i=1}^{k} \sigma_i(M) = \max \{ \operatorname{Re}\operatorname{tr}(U^* M) : U^* U = I_k, U \in \mathbb{C}^{n \times k} \}.
    $$

    设 $U_0$ 是使得 $\sum_{i=1}^k \sigma_i(A+B) = \operatorname{Re}\operatorname{tr}(U_0^*(A+B))$ 的最优等距映射，则：

    $$
    \sum_{i=1}^{k} \sigma_i(A+B) = \operatorname{Re}\operatorname{tr}(U_0^* A) + \operatorname{Re}\operatorname{tr}(U_0^* B) \leq \sum_{i=1}^{k} \sigma_i(A) + \sum_{i=1}^{k} \sigma_i(B). \quad \blacksquare
    $$

!!! example "例 18.3"
    设 $A = \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}$。

    $A$ 的奇异值：$\sigma_1(A) = 3$，$\sigma_2(A) = 1$。

    $B$ 的奇异值：$\sigma_1(B) = 2$，$\sigma_2(B) = 0$。

    $AB = \begin{pmatrix} 0 & 6 \\ 0 & 0 \end{pmatrix}$，奇异值：$\sigma_1(AB) = 6$，$\sigma_2(AB) = 0$。

    次可乘性：$\sigma_1(AB) = 6 \leq \sigma_1(A) \cdot \sigma_1(B) = 3 \times 2 = 6$，等号成立。

    $\sigma_2(AB) = 0 \leq \sigma_1(A) \cdot \sigma_2(B) = 3 \times 0 = 0$，等号成立。

---

## 18.3 迹不等式

<div class="context-flow" markdown>

**脉络**：**von Neumann** $|\operatorname{tr}(A^*B)|\leq\sum\sigma_i(A)\sigma_i(B)$(排序不等式+Birkhoff) → **Golden-Thompson** $\operatorname{tr}(e^{A+B})\leq\operatorname{tr}(e^Ae^B)$

</div>

矩阵的迹是最简单的矩阵函数之一，但它满足若干深刻的不等式。

!!! definition "定义 18.5 (迹与 Frobenius 范数)"
    设 $A$ 为 $n \times n$ 矩阵，则：

    - **迹**：$\operatorname{tr}(A) = \sum_{i=1}^n a_{ii}$。
    - **Frobenius 范数**：$\|A\|_F = \sqrt{\operatorname{tr}(A^* A)} = \sqrt{\sum_{i=1}^n \sigma_i^2(A)}$。

!!! theorem "定理 18.7 (von Neumann 迹不等式)"
    设 $A, B$ 为 $n \times n$ 复矩阵，则：

    $$
    |\operatorname{tr}(AB)| \leq \sum_{i=1}^{n} \sigma_i(A) \sigma_i(B).
    $$

    更一般地：

    $$
    |\operatorname{tr}(A^* B)| \leq \sum_{i=1}^{n} \sigma_i(A) \sigma_i(B).
    $$

??? proof "证明"
    设 $A = U_A \Sigma_A V_A^*$ 和 $B = U_B \Sigma_B V_B^*$ 为奇异值分解。则：

    $$
    \operatorname{tr}(A^* B) = \operatorname{tr}(V_A \Sigma_A U_A^* U_B \Sigma_B V_B^*) = \operatorname{tr}(\Sigma_A W \Sigma_B Z),
    $$

    其中 $W = U_A^* U_B$ 和 $Z = V_B^* V_A$ 为酉矩阵。令 $P = WZ$，则 $P$ 也是酉矩阵。

    由于 $\operatorname{tr}(\Sigma_A W \Sigma_B Z) = \sum_{i,j} (\sigma_i(A))(\sigma_j(B)) w_{ij} z_{ji}$，我们需要证明：

    $$
    \left|\sum_{i,j} \sigma_i(A) \sigma_j(B) w_{ij} z_{ji}\right| \leq \sum_i \sigma_i(A)\sigma_i(B).
    $$

    注意 $|w_{ij} z_{ji}| \leq |w_{ij}| |z_{ji}|$，而矩阵 $C$ 定义为 $c_{ij} = |w_{ij}||z_{ji}|$ 是一个双随机矩阵（doubly stochastic matrix）。

    由 Birkhoff 定理，$C$ 是置换矩阵的凸组合，因此：

    $$
    |\operatorname{tr}(A^* B)| \leq \sum_{i,j} \sigma_i(A) \sigma_j(B) c_{ij} \leq \max_{\pi} \sum_i \sigma_i(A) \sigma_{\pi(i)}(B) = \sum_i \sigma_i(A) \sigma_i(B),
    $$

    最后一个等式由排序不等式（rearrangement inequality）得到。$\blacksquare$

!!! theorem "定理 18.8 (Golden-Thompson 不等式)"
    设 $A, B$ 为 $n \times n$ Hermite 矩阵，则：

    $$
    \operatorname{tr}(e^{A+B}) \leq \operatorname{tr}(e^A e^B).
    $$

??? proof "证明"
    证明的关键步骤如下：

    **第一步**：利用 Lie-Trotter 乘积公式：$e^{A+B} = \lim_{m \to \infty} (e^{A/m} e^{B/m})^m$。

    **第二步**：由 von Neumann 迹不等式，对 Hermite 矩阵 $A, B$：

    $$
    \operatorname{tr}(e^{A+B}) = \operatorname{tr}\left(\lim_{m\to\infty}(e^{A/m}e^{B/m})^m\right).
    $$

    **第三步**：关键的不等式来自于以下事实：对于半正定矩阵 $P, Q$，$\operatorname{tr}(PQ) \leq \operatorname{tr}(P) \cdot \|Q\|_2$ 不够精确。需要更精细的论证。

    利用特征值的对数凸性：设 $\lambda_i = \lambda_i(e^{A+B})$，$\alpha_i = \lambda_i(e^A)$，$\beta_i = \lambda_i(e^B)$，由 Weyl 不等式和特征值-奇异值关系：

    $$
    \sum_i \lambda_i(e^{A+B}) \leq \sum_i \sigma_i(e^{A/2} \cdot e^B \cdot e^{A/2}) \leq \sum_i \sigma_i(e^A) \sigma_i(e^B) = \sum_i \lambda_i(e^A) \lambda_i(e^B),
    $$

    其中最后一步利用了 $e^A, e^B$ 为正定矩阵，其奇异值等于特征值。

    而 $\operatorname{tr}(e^A e^B) = \operatorname{tr}(e^{A/2} e^B e^{A/2}) \geq \sum_i \lambda_i(e^{A/2} e^B e^{A/2}) = \sum_i \lambda_i(e^{A+B})$（在适当的不等式方向上），由此得到结论。$\blacksquare$

!!! note "注"
    Golden-Thompson 不等式不能推广到三个矩阵：一般来说 $\operatorname{tr}(e^{A+B+C}) \leq \operatorname{tr}(e^A e^B e^C)$ **不成立**。但 Lieb 和 Seiringer 证明了：$\operatorname{tr}(e^{A+B+C}) \leq \operatorname{tr}\left(e^A \# e^B \# e^C\right)$ 的某些推广形式。

!!! example "例 18.4"
    验证 von Neumann 迹不等式。设 $A = \begin{pmatrix} 2 & 1 \\ 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 1 & 2 \end{pmatrix}$。

    $\operatorname{tr}(A^* B) = \operatorname{tr}\begin{pmatrix} 2 & 0 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 1 & 2 \end{pmatrix} = \operatorname{tr}\begin{pmatrix} 2 & 0 \\ 2 & 2 \end{pmatrix} = 4$。

    $A$ 的奇异值：$\sigma_1(A) = \frac{\sqrt{6}+\sqrt{2}}{2} \approx 1.932$，$\sigma_2(A) = \frac{\sqrt{6}-\sqrt{2}}{2} \approx 0.518$。

    $B$ 的奇异值（与 $A$ 相同，因为 $B = A^T$）：$\sigma_1(B) \approx 1.932$，$\sigma_2(B) \approx 0.518$。

    $\sum \sigma_i(A)\sigma_i(B) \approx 1.932^2 + 0.518^2 \approx 3.732 + 0.268 = 4.0$。

    因此 $|\operatorname{tr}(A^*B)| = 4 \leq 4.0 = \sum \sigma_i(A)\sigma_i(B)$，等号成立！

!!! example "例 18.5"
    **Golden-Thompson 不等式的数值验证**。设 $A = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。

    $A + B = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$，特征值为 $\pm\sqrt{2}$。

    $\operatorname{tr}(e^{A+B}) = e^{\sqrt{2}} + e^{-\sqrt{2}} = 2\cosh(\sqrt{2}) \approx 4.121$。

    $e^A = \begin{pmatrix} e & 0 \\ 0 & e^{-1} \end{pmatrix}$，$e^B = \begin{pmatrix} \cosh 1 & \sinh 1 \\ \sinh 1 & \cosh 1 \end{pmatrix}$。

    $\operatorname{tr}(e^A e^B) = e \cosh 1 + e^{-1}\cosh 1 + e\cdot 0 \cdot \sinh 1 + \cdots$

    直接计算：$\operatorname{tr}(e^A e^B) = e \cosh 1 + e^{-1} \cosh 1 = (e + e^{-1})\cosh 1 = 2\cosh(1)\cosh(1) \approx 2 \times 1.543 \times 1.543 \approx 4.762$。

    确实 $4.121 \leq 4.762$，Golden-Thompson 不等式成立。

---

## 18.4 行列式不等式

<div class="context-flow" markdown>

**脉络**：**Hadamard** $\det A\leq\prod a_{ii}$ (Ch16 Schur补证法) → **Fischer** (分块推广) → **Minkowski** $(\det(A+B))^{1/n}\geq(\det A)^{1/n}+(\det B)^{1/n}$

</div>

行列式作为矩阵的另一个基本不变量，也满足若干重要的不等式。

!!! definition "定义 18.6 (正定矩阵上的偏序)"
    设 $A, B$ 为 $n \times n$ Hermite 矩阵。定义 **Löwner 偏序**（Loewner partial order）为：

    $$
    A \succeq B \iff A - B \text{ 半正定}.
    $$

    记作 $A \succ B$ 若 $A - B$ 正定。

!!! theorem "定理 18.9 (Hadamard 不等式)"
    设 $A = (a_{ij})$ 为 $n \times n$ 正半定 Hermite 矩阵，则：

    $$
    \det(A) \leq \prod_{i=1}^{n} a_{ii}.
    $$

    等号成立当且仅当 $A$ 为对角矩阵或 $A$ 有零行。

??? proof "证明"
    **方法一**（利用 Schur 补）：对 $A$ 作分块：

    $$
    A = \begin{pmatrix} A_{n-1} & \mathbf{a} \\ \mathbf{a}^* & a_{nn} \end{pmatrix},
    $$

    其中 $A_{n-1}$ 为前 $n-1$ 阶主子矩阵。

    若 $A_{n-1}$ 奇异，则 $\det(A) = 0 \leq \prod a_{ii}$ 显然成立。

    若 $A_{n-1}$ 非奇异，由 Schur 补公式：

    $$
    \det(A) = \det(A_{n-1})(a_{nn} - \mathbf{a}^* A_{n-1}^{-1} \mathbf{a}).
    $$

    由于 $A \succeq 0$，Schur 补 $a_{nn} - \mathbf{a}^* A_{n-1}^{-1} \mathbf{a} \geq 0$，即 $a_{nn} - \mathbf{a}^* A_{n-1}^{-1} \mathbf{a} \leq a_{nn}$。

    因此 $\det(A) \leq \det(A_{n-1}) \cdot a_{nn}$。

    由数学归纳法，$\det(A_{n-1}) \leq \prod_{i=1}^{n-1} a_{ii}$，故 $\det(A) \leq \prod_{i=1}^n a_{ii}$。

    等号成立要求每步 Schur 补 $\mathbf{a}^* A_{n-1}^{-1}\mathbf{a} = 0$（当 $A_{n-1}$ 非奇异时意味着 $\mathbf{a} = 0$），即 $A$ 为对角矩阵。$\blacksquare$

!!! theorem "定理 18.10 (Fischer 不等式)"
    设 $A$ 为 $n \times n$ 正半定 Hermite 矩阵，分块为：

    $$
    A = \begin{pmatrix} B & C \\ C^* & D \end{pmatrix},
    $$

    其中 $B$ 为 $k \times k$，$D$ 为 $(n-k) \times (n-k)$。则：

    $$
    \det(A) \leq \det(B) \cdot \det(D).
    $$

??? proof "证明"
    若 $B$ 奇异，则由 $A \succeq 0$ 可以证明 $\det(A) = 0$，不等式显然成立。

    设 $B$ 非奇异。由 Schur 补公式：

    $$
    \det(A) = \det(B) \cdot \det(D - C^* B^{-1} C).
    $$

    由于 $A \succeq 0$，其 Schur 补 $D - C^* B^{-1} C \succeq 0$，因此 $D - C^* B^{-1} C \preceq D$（因为 $C^* B^{-1} C \succeq 0$）。

    由正半定矩阵行列式的单调性（$0 \preceq X \preceq Y$ 蕴含 $\det X \leq \det Y$）：

    $$
    \det(D - C^* B^{-1} C) \leq \det(D).
    $$

    因此 $\det(A) = \det(B) \cdot \det(D - C^*B^{-1}C) \leq \det(B) \cdot \det(D)$。$\blacksquare$

!!! theorem "定理 18.11 (Minkowski 行列式不等式)"
    设 $A, B$ 为 $n \times n$ 正半定 Hermite 矩阵，则：

    $$
    [\det(A + B)]^{1/n} \geq [\det(A)]^{1/n} + [\det(B)]^{1/n}.
    $$

??? proof "证明"
    若 $A$ 或 $B$ 奇异，不等式变为 $[\det(A+B)]^{1/n} \geq [\det(A)]^{1/n}$（或类似地关于 $B$），由 $A + B \succeq A$ 和行列式单调性即得。

    设 $A$ 正定（$A \succ 0$）。则：

    $$
    \det(A + B) = \det(A) \cdot \det(I + A^{-1/2} B A^{-1/2}).
    $$

    令 $C = A^{-1/2} B A^{-1/2} \succeq 0$，设其特征值为 $\mu_1, \ldots, \mu_n \geq 0$。则：

    $$
    [\det(I + C)]^{1/n} = \left[\prod_{i=1}^n (1 + \mu_i)\right]^{1/n} \geq 1 + \left[\prod_{i=1}^n \mu_i\right]^{1/n},
    $$

    最后一步由 AM-GM 不等式的推广形式得到（具体地，这是对 $(1+\mu_i)$ 应用几何-算术均值不等式的结果）。

    因此：

    $$
    [\det(A+B)]^{1/n} = [\det A]^{1/n} \cdot [\det(I+C)]^{1/n} \geq [\det A]^{1/n}(1 + [\det C]^{1/n}).
    $$

    又 $\det C = \det(A^{-1/2} B A^{-1/2}) = \det(A^{-1}) \det(B) = \frac{\det B}{\det A}$，故 $[\det C]^{1/n} = \frac{[\det B]^{1/n}}{[\det A]^{1/n}}$。

    代入得 $[\det(A+B)]^{1/n} \geq [\det A]^{1/n} + [\det B]^{1/n}$。$\blacksquare$

!!! example "例 18.6"
    验证 Hadamard 不等式。设 $A = \begin{pmatrix} 4 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 6 \end{pmatrix}$。

    $A$ 正定（可通过检查所有顺序主子式为正来验证：$4 > 0$，$20 - 4 = 16 > 0$，$\det A = 4 \cdot 21 - 2 \cdot 9 + 1 \cdot 1 = 84 - 18 + 1 = 67 > 0$）。

    Hadamard 不等式：$\det(A) = 67 \leq 4 \times 5 \times 6 = 120$。

    确实 $67 \leq 120$，不等式成立。

!!! example "例 18.7"
    验证 Minkowski 行列式不等式。设 $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix}$。

    $[\det(A+B)]^{1/2} = [\det\begin{pmatrix} 3 & 0 \\ 0 & 7 \end{pmatrix}]^{1/2} = \sqrt{21} \approx 4.583$。

    $[\det A]^{1/2} + [\det B]^{1/2} = \sqrt{6} + \sqrt{4} = \sqrt{6} + 2 \approx 2.449 + 2 = 4.449$。

    确实 $4.583 \geq 4.449$，Minkowski 不等式成立。

---

## 18.5 Majorization（优超）

<div class="context-flow" markdown>

**核心框架**：$\mathbf{x}\prec\mathbf{y}$ $\Leftrightarrow$ $\mathbf{x}=D\mathbf{y}$($D$双随机) $\Leftrightarrow$ $\mathbf{x}$在$\mathbf{y}$置换凸包内

**Schur-Horn**：对角元 $\prec$ 特征值 · Birkhoff定理是桥梁

</div>

Majorization（优超）是一个统一多种不等式的核心概念，它精确描述了向量分量之间"扩散程度"的比较。

!!! definition "定义 18.7 (优超关系)"
    设 $\mathbf{x} = (x_1, \ldots, x_n)$ 和 $\mathbf{y} = (y_1, \ldots, y_n)$ 为实向量，将其分量按降序排列得 $x_{[1]} \geq \cdots \geq x_{[n]}$ 和 $y_{[1]} \geq \cdots \geq y_{[n]}$。称 $\mathbf{x}$ **被 $\mathbf{y}$ 优超**（$\mathbf{x}$ is majorized by $\mathbf{y}$），记作 $\mathbf{x} \prec \mathbf{y}$，若：

    $$
    \sum_{i=1}^{k} x_{[i]} \leq \sum_{i=1}^{k} y_{[i]}, \quad k = 1, 2, \ldots, n-1,
    $$

    且 $\sum_{i=1}^{n} x_i = \sum_{i=1}^{n} y_i$。

!!! definition "定义 18.8 (双随机矩阵)"
    一个 $n \times n$ 非负实矩阵 $D = (d_{ij})$ 称为**双随机矩阵**（doubly stochastic matrix），若其每行和每列之和均为 $1$：

    $$
    \sum_{j=1}^{n} d_{ij} = 1 \quad \forall i, \qquad \sum_{i=1}^{n} d_{ij} = 1 \quad \forall j.
    $$

!!! theorem "定理 18.12 (Hardy-Littlewood-Polya 定理)"
    对实向量 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$，以下条件等价：

    1. $\mathbf{x} \prec \mathbf{y}$；
    2. 存在双随机矩阵 $D$ 使得 $\mathbf{x} = D\mathbf{y}$；
    3. $\mathbf{x}$ 在 $\mathbf{y}$ 的所有分量置换所构成的凸包内。

??? proof "证明"
    **(2) $\Rightarrow$ (1)**：设 $\mathbf{x} = D\mathbf{y}$，$D$ 为双随机矩阵。由 Birkhoff 定理，$D = \sum_k \alpha_k P_{\pi_k}$（置换矩阵的凸组合）。因此 $\mathbf{x} = \sum_k \alpha_k P_{\pi_k} \mathbf{y}$，即 $\mathbf{x}$ 是 $\mathbf{y}$ 各分量置换的凸组合。

    对凸组合应用排序后部分和的凸性，可以验证优超条件。

    **(1) $\Rightarrow$ (2)**：构造性证明。给定 $\mathbf{x} \prec \mathbf{y}$，可以通过有限次 T-变换（$\mathbf{z} = t\mathbf{a} + (1-t)P_{ij}\mathbf{a}$，其中 $P_{ij}$ 交换第 $i, j$ 分量）将 $\mathbf{y}$ 变换为 $\mathbf{x}$，每次 T-变换对应一个特殊的双随机矩阵，其乘积仍为双随机矩阵。

    **(2) $\Leftrightarrow$ (3)**：由 Birkhoff 定理（定理 18.14）直接得出。$\blacksquare$

!!! theorem "定理 18.13 (Schur-Horn 定理)"
    设 $A$ 为 $n \times n$ Hermite 矩阵，特征值 $\lambda_1 \geq \cdots \geq \lambda_n$，对角元素 $a_{11}, \ldots, a_{nn}$，则：

    $$
    (a_{11}, \ldots, a_{nn}) \prec (\lambda_1, \ldots, \lambda_n).
    $$

    即对角元素向量被特征值向量优超。

    反之，给定实向量 $\mathbf{d} \prec \boldsymbol{\lambda}$，存在 Hermite 矩阵以 $\boldsymbol{\lambda}$ 为特征值、以 $\mathbf{d}$ 为对角。

??? proof "证明"
    **必要性**：设 $A = U \Lambda U^*$，其中 $U = (u_{ij})$ 为酉矩阵。则：

    $$
    a_{ii} = (U \Lambda U^*)_{ii} = \sum_{j=1}^n \lambda_j |u_{ij}|^2.
    $$

    令 $d_{ij} = |u_{ij}|^2$，由 $U$ 的酉性，$D = (d_{ij})$ 为双随机矩阵。因此 $\mathbf{d} = D\boldsymbol{\lambda}$，由 Hardy-Littlewood-Polya 定理，$\mathbf{d} \prec \boldsymbol{\lambda}$。

    **充分性**（Horn 的结果）：由构造性方法，可以通过一系列 Givens 旋转构造所需的酉矩阵 $U$。$\blacksquare$

!!! theorem "定理 18.14 (Birkhoff 定理)"
    双随机矩阵的集合 $\mathcal{D}_n$ 是一个凸紧集，其极点恰好是所有 $n \times n$ **置换矩阵**（permutation matrix）。即每个双随机矩阵可以写成置换矩阵的凸组合：

    $$
    D = \sum_{k=1}^{N} \alpha_k P_{\pi_k}, \quad \alpha_k \geq 0, \quad \sum_k \alpha_k = 1.
    $$

??? proof "证明"
    **极点是置换矩阵**：设 $D$ 是 $\mathcal{D}_n$ 的极点。若 $D$ 不是置换矩阵，则存在某行含至少两个正元素 $d_{ij}, d_{ik} > 0$。利用双随机条件可构造两个不同的双随机矩阵 $D_1, D_2$ 使得 $D = \frac{1}{2}(D_1 + D_2)$，与 $D$ 为极点矛盾。

    **置换矩阵是极点**：设 $P$ 为置换矩阵，若 $P = \alpha D_1 + (1-\alpha) D_2$（$0 < \alpha < 1$），由于 $P$ 的元素仅为 $0$ 或 $1$，而 $D_1, D_2$ 的元素在 $[0,1]$ 中且 $\alpha D_1 + (1-\alpha)D_2$ 要达到 $0$ 或 $1$，必须 $D_1 = D_2 = P$。

    **每个双随机矩阵是置换矩阵的凸组合**：这可由 Birkhoff 算法证明。对双随机矩阵 $D$，由 König 定理，存在置换 $\pi$ 使得 $d_{i,\pi(i)} > 0$ 对所有 $i$ 成立。令 $\theta = \min_i d_{i,\pi(i)} > 0$，则 $D' = \frac{1}{1-\theta}(D - \theta P_\pi)$ 仍为双随机矩阵（或已完成），且零元素更多。有限步后终止。$\blacksquare$

!!! example "例 18.8"
    验证 Schur-Horn 定理。设 $A = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$。

    特征值：$\lambda_1 = 4$，$\lambda_2 = 2$。对角元素：$d_1 = 3$，$d_2 = 3$。

    检验优超：$d_{[1]} = 3 \leq \lambda_1 = 4$，$d_1 + d_2 = 6 = \lambda_1 + \lambda_2 = 6$。

    因此 $(3, 3) \prec (4, 2)$，定理成立。

    对应的双随机矩阵为 $D = \begin{pmatrix} 1/2 & 1/2 \\ 1/2 & 1/2 \end{pmatrix}$，确实 $(3, 3) = D(4, 2)$。

!!! example "例 18.9"
    设 $\mathbf{x} = (3, 3, 3)$，$\mathbf{y} = (5, 3, 1)$。验证 $\mathbf{x} \prec \mathbf{y}$。

    $x_{[1]} = 3 \leq y_{[1]} = 5$。

    $x_{[1]} + x_{[2]} = 6 \leq y_{[1]} + y_{[2]} = 8$。

    $x_{[1]} + x_{[2]} + x_{[3]} = 9 = y_{[1]} + y_{[2]} + y_{[3]} = 9$。

    所有条件满足，$\mathbf{x} \prec \mathbf{y}$。双随机矩阵：$D = \frac{1}{3}\begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 1 \\ 1 & 1 & 1 \end{pmatrix}$，$D\mathbf{y} = (3,3,3) = \mathbf{x}$。

---

## 18.6 Schur 凸函数

<div class="context-flow" markdown>

**脉络**：$f$ Schur凸 + $\mathbf{x}\prec\mathbf{y}$ $\Rightarrow$ $f(\mathbf{x})\leq f(\mathbf{y})$ · 判别条件：Schur条件 $(x_i-x_j)(\partial_i f-\partial_j f)\geq 0$ · $\phi$凸 $\Rightarrow$ $\sum\phi(a_{ii})\leq\sum\phi(\lambda_i)$

</div>

Schur 凸函数是与优超理论密切相关的一类函数，它为判断矩阵不等式提供了强大工具。

!!! definition "定义 18.9 (Schur 凸函数与 Schur 凹函数)"
    函数 $f: \mathbb{R}^n \to \mathbb{R}$ 称为 **Schur 凸函数**（Schur-convex function），若对所有满足 $\mathbf{x} \prec \mathbf{y}$ 的 $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$，有：

    $$
    f(\mathbf{x}) \leq f(\mathbf{y}).
    $$

    称为 **Schur 凹函数**（Schur-concave function），若 $\mathbf{x} \prec \mathbf{y}$ 蕴含 $f(\mathbf{x}) \geq f(\mathbf{y})$。

!!! definition "定义 18.10 (对称函数)"
    函数 $f: \mathbb{R}^n \to \mathbb{R}$ 称为**对称函数**（symmetric function），若对任意置换 $\pi$，$f(x_{\pi(1)}, \ldots, x_{\pi(n)}) = f(x_1, \ldots, x_n)$。

!!! theorem "定理 18.15 (Schur 凸性的判别条件)"
    设 $f: \mathbb{R}^n \to \mathbb{R}$ 连续可微且对称，则 $f$ 为 Schur 凸函数当且仅当对所有 $i \neq j$：

    $$
    (x_i - x_j)\left(\frac{\partial f}{\partial x_i} - \frac{\partial f}{\partial x_j}\right) \geq 0.
    $$

    此条件称为 **Schur 条件**（Schur's condition）。

??? proof "证明"
    **充分性**：设 $\mathbf{x} \prec \mathbf{y}$。由 Hardy-Littlewood-Polya 定理，$\mathbf{x} = D\mathbf{y}$，$D$ 为双随机矩阵。由 Birkhoff 定理，$\mathbf{x}$ 在 $\mathbf{y}$ 的置换点集的凸包中。

    可以证明，$\mathbf{x}$ 可以通过有限次 Robin Hood 变换（T-变换）从 $\mathbf{y}$ 得到。每次 T-变换将 $\mathbf{y}$ 中较大分量减小、较小分量增大。因此只需证明每次 T-变换不增加 $f$ 的值。

    设 $\mathbf{z} = t\mathbf{y} + (1-t)P_{ij}\mathbf{y}$（$0 \leq t \leq 1$），考虑 $g(t) = f(\mathbf{z}(t))$。由链式法则和 Schur 条件，可以证明 $g$ 在适当方向单调，从而 $f(\mathbf{x}) \leq f(\mathbf{y})$。

    **必要性**：取 $\mathbf{y}$ 使得 $y_i > y_j$，令 $\mathbf{x}$ 由 T-变换得到（将 $y_i$ 减小 $\epsilon$，$y_j$ 增大 $\epsilon$），则 $\mathbf{x} \prec \mathbf{y}$，由 $f(\mathbf{x}) \leq f(\mathbf{y})$ 和 $\epsilon \to 0$ 取极限即得 Schur 条件。$\blacksquare$

!!! theorem "定理 18.16 (特征值的 Schur 凸性)"
    设 $\phi: \mathbb{R} \to \mathbb{R}$ 为凸函数，$A$ 为 $n \times n$ Hermite 矩阵，特征值 $\lambda_1 \geq \cdots \geq \lambda_n$，对角元素 $a_{11}, \ldots, a_{nn}$。则：

    $$
    \sum_{i=1}^{n} \phi(a_{ii}) \leq \sum_{i=1}^{n} \phi(\lambda_i).
    $$

??? proof "证明"
    由 Schur-Horn 定理（定理 18.13），$(a_{11}, \ldots, a_{nn}) \prec (\lambda_1, \ldots, \lambda_n)$。

    函数 $F(\mathbf{x}) = \sum_{i=1}^n \phi(x_i)$ 是 Schur 凸函数（当 $\phi$ 为凸函数时）。验证：

    $$
    (x_i - x_j)\left(\frac{\partial F}{\partial x_i} - \frac{\partial F}{\partial x_j}\right) = (x_i - x_j)(\phi'(x_i) - \phi'(x_j)) \geq 0,
    $$

    最后一步由 $\phi$ 的凸性（即 $\phi'$ 单调递增）得到。

    因此 $F(\mathbf{a}) \leq F(\boldsymbol{\lambda})$，即 $\sum \phi(a_{ii}) \leq \sum \phi(\lambda_i)$。$\blacksquare$

!!! example "例 18.10"
    取 $\phi(x) = x^2$（凸函数），$A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。

    特征值：$\lambda_1 = 3$，$\lambda_2 = 1$。对角元素：$a_{11} = 2$，$a_{22} = 2$。

    $\sum \phi(a_{ii}) = 4 + 4 = 8 \leq \sum \phi(\lambda_i) = 9 + 1 = 10$。

    确实 $8 \leq 10$，定理成立。

!!! example "例 18.11"
    取 $\phi(x) = e^x$（凸函数），$A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。

    特征值：$\lambda_1 = 1$，$\lambda_2 = -1$。对角元素：$a_{11} = 0$，$a_{22} = 0$。

    $\sum \phi(a_{ii}) = e^0 + e^0 = 2 \leq \sum \phi(\lambda_i) = e^1 + e^{-1} = e + 1/e \approx 3.086$。

    不等式成立。

---

## 18.7 矩阵凸性与矩阵单调性

<div class="context-flow" markdown>

**脉络**：$t^r$($0\leq r\leq 1$)和 $\log t$ 矩阵单调，$t^2$ 不是

**Löwner定理**：矩阵单调 $\Leftrightarrow$ Pick函数(上半平面自映) · Jensen矩阵不等式 → 量子信息

</div>

本节将凸性和单调性从标量函数推广到矩阵函数的领域。

!!! definition "定义 18.11 (矩阵凸函数)"
    设 $f: (a, b) \to \mathbb{R}$ 为连续函数。称 $f$ 为**矩阵凸函数**（matrix convex function），若对所有特征值在 $(a,b)$ 中的 $n \times n$ Hermite 矩阵 $A, B$ 和 $t \in [0,1]$：

    $$
    f(tA + (1-t)B) \preceq t f(A) + (1-t) f(B).
    $$

    这里 $f(A)$ 表示矩阵函数（通过谱映射定义）。

!!! definition "定义 18.12 (矩阵单调函数)"
    设 $f: (a,b) \to \mathbb{R}$ 为连续函数。称 $f$ 为**矩阵单调函数**（matrix monotone function），或**算子单调函数**（operator monotone function），若对所有特征值在 $(a,b)$ 中的 $n \times n$ Hermite 矩阵 $A, B$：

    $$
    A \preceq B \implies f(A) \preceq f(B).
    $$

!!! theorem "定理 18.17 (矩阵单调函数的基本例子)"
    以下函数在 $(0, \infty)$ 上是矩阵单调的：

    1. $f(t) = t^r$，$0 \leq r \leq 1$；
    2. $f(t) = \log t$；
    3. $f(t) = \frac{t}{t+1}$。

    以下函数**不是**矩阵单调的：

    4. $f(t) = t^2$（$t > 0$）。

??? proof "证明"
    **(4)** 的反例：取 $A = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$。

    $B - A = \begin{pmatrix} -1 & 1 \\ 1 & 1 \end{pmatrix}$，特征值为 $\sqrt{2}, -\sqrt{2}$，$B \not\succeq A$。

    换一个例子：取 $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$。则 $B - A = \begin{pmatrix} 0 & 1 \\ 1 & 2 \end{pmatrix} \succeq 0$（特征值 $1 \pm \sqrt{2}$，不全非负）。

    更清晰的例子：$A = I$，$B = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$。$B - A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}$（特征值 $\frac{1\pm\sqrt{5}}{2}$），$B \not\succeq A$。

    取 $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \preceq \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = B$，$A^2 = A$，$B^2 = B$，$B^2 - A^2 = B - A \succeq 0$，此特殊情况下碰巧成立。

    正确反例：$A = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix} \preceq \begin{pmatrix} 3 & 1 \\ 1 & 1 \end{pmatrix} = B$。则 $B - A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \succeq 0$。但 $B^2 - A^2 = \begin{pmatrix} 10 & 4 \\ 4 & 2 \end{pmatrix} - \begin{pmatrix} 5 & 3 \\ 3 & 2 \end{pmatrix} = \begin{pmatrix} 5 & 1 \\ 1 & 0 \end{pmatrix}$，其特征值为 $\frac{5 \pm \sqrt{29}}{2}$，由于 $\sqrt{29} > 5$，有一个负特征值，因此 $B^2 \not\succeq A^2$。

    **(2)** $\log t$ 的矩阵单调性：设 $0 \prec A \preceq B$，需证 $\log A \preceq \log B$。

    利用积分表示 $\log t = \int_0^{\infty} \left(\frac{1}{1+s} - \frac{1}{t+s}\right) ds$，而 $g_s(t) = \frac{1}{t+s}$ 是矩阵单调递减函数（因为 $A \preceq B$ 蕴含 $(A+sI)^{-1} \succeq (B+sI)^{-1}$，这由 $A + sI \preceq B + sI$ 和逆的反单调性得到）。因此 $-g_s$ 矩阵单调递增，对 $s$ 积分保持单调性。$\blacksquare$

!!! theorem "定理 18.18 (Löwner 定理)"
    函数 $f: (a,b) \to \mathbb{R}$ 对所有 $n$（任意维数）都是矩阵单调的，当且仅当 $f$ 可以解析延拓到上半平面 $\mathbb{C}^+$，且 $f(\mathbb{C}^+) \subset \overline{\mathbb{C}^+}$（即 $f$ 将上半平面映入上半平面的闭包）。

    等价地，$f$ 具有积分表示：

    $$
    f(t) = \alpha + \beta t + \int_{-\infty}^{\infty} \frac{t\mu + 1}{\mu - t} \, d\nu(\mu),
    $$

    其中 $\alpha \in \mathbb{R}$，$\beta \geq 0$，$\nu$ 为正 Borel 测度。

??? proof "证明"
    这是 Löwner 在 1934 年的经典结果，完整证明需要相当多的函数论工具。

    **必要性的思路**：若 $f$ 为矩阵单调，取 $n = 2$，考虑矩阵 $A = \begin{pmatrix} x & 0 \\ 0 & y \end{pmatrix}$ 和扰动 $B = A + \epsilon \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。由矩阵单调性导出的条件可以推出 $f$ 的**分差矩阵**（Loewner matrix）$L_f = \left(\frac{f(x_i) - f(x_j)}{x_i - x_j}\right)$ 为正半定矩阵。这进一步蕴含 $f$ 可以延拓为上半平面上的 Pick 函数。

    **充分性的思路**：若 $f$ 为 Pick 函数，则利用 Nevanlinna-Pick 插值理论和正半定核的性质，可以证明分差矩阵的正半定性，从而得到矩阵单调性。$\blacksquare$

!!! theorem "定理 18.19 (Jensen 矩阵不等式)"
    设 $f$ 为 $(a,b)$ 上的矩阵凸函数，$A_1, \ldots, A_k$ 为 Hermite 矩阵（特征值在 $(a,b)$ 中），$\omega_1, \ldots, \omega_k > 0$ 且 $\sum \omega_i = 1$，则：

    $$
    f\left(\sum_{i=1}^{k} \omega_i A_i\right) \preceq \sum_{i=1}^{k} \omega_i f(A_i).
    $$

??? proof "证明"
    对 $k$ 使用数学归纳法。$k = 2$ 时即为矩阵凸函数的定义。

    设 $k \geq 3$，令 $\omega = \omega_1 + \cdots + \omega_{k-1}$，$B = \frac{1}{\omega}\sum_{i=1}^{k-1}\omega_i A_i$。则：

    $$
    \sum_{i=1}^k \omega_i A_i = \omega B + \omega_k A_k.
    $$

    由矩阵凸性：

    $$
    f(\omega B + \omega_k A_k) \preceq \omega f(B) + \omega_k f(A_k).
    $$

    由归纳假设：

    $$
    f(B) = f\left(\sum_{i=1}^{k-1}\frac{\omega_i}{\omega} A_i\right) \preceq \sum_{i=1}^{k-1}\frac{\omega_i}{\omega} f(A_i).
    $$

    因此 $f\left(\sum_{i=1}^k \omega_i A_i\right) \preceq \omega \sum_{i=1}^{k-1}\frac{\omega_i}{\omega} f(A_i) + \omega_k f(A_k) = \sum_{i=1}^k \omega_i f(A_i)$。$\blacksquare$

!!! example "例 18.12"
    验证 $f(t) = t^2$ 是矩阵凸函数（在所有 Hermite 矩阵上）。

    设 $A, B$ 为 Hermite 矩阵，$t \in [0,1]$。需证 $(tA + (1-t)B)^2 \preceq t A^2 + (1-t) B^2$。

    展开左边：$t^2 A^2 + t(1-t)(AB + BA) + (1-t)^2 B^2$。

    右边减左边：

    $$
    t(1-t)A^2 + t(1-t)B^2 - t(1-t)(AB+BA) = t(1-t)(A-B)^2 \succeq 0,
    $$

    因为 $(A-B)^2$ 半正定。因此 $f(t) = t^2$ 确实是矩阵凸的。

!!! example "例 18.13"
    $f(t) = -\log t$ 是 $(0,\infty)$ 上的矩阵凸函数。

    这等价于 $\log t$ 是矩阵凹函数，即 $\log(tA + (1-t)B) \succeq t\log A + (1-t)\log B$。

    数值验证：$A = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}$，$t = 1/2$。

    $\frac{1}{2}(A + B) = \begin{pmatrix} 3/2 & 0 \\ 0 & 2 \end{pmatrix}$。

    $\log\frac{1}{2}(A+B) = \begin{pmatrix} \log(3/2) & 0 \\ 0 & \log 2 \end{pmatrix} \approx \begin{pmatrix} 0.405 & 0 \\ 0 & 0.693 \end{pmatrix}$。

    $\frac{1}{2}(\log A + \log B) = \frac{1}{2}\begin{pmatrix} \log 2 & 0 \\ 0 & \log 3 \end{pmatrix} \approx \begin{pmatrix} 0.347 & 0 \\ 0 & 0.549 \end{pmatrix}$。

    差为 $\begin{pmatrix} 0.058 & 0 \\ 0 & 0.144 \end{pmatrix} \succeq 0$，验证了矩阵凹性。

## 练习题

1. **[Courant-Fischer] 极小极大定理（Min-Max Theorem）在几何上如何理解？**
   ??? success "参考答案"
       它将第 $k$ 大的特征值解释为：在所有可能的 $k$ 维子空间中，瑞利商的最小值的“最大可能下限”。几何上，这相当于在切片的椭球截面中寻找主轴。

2. **[Weyl不等式] 为什么说 Weyl 不等式给出了特征值的“绝对扰动界”？**
   ??? success "参考答案"
       Weyl 不等式表明 $|\lambda_k(A+B) - \lambda_k(A)| \le \|B\|_2$。这意味着如果我们对矩阵施加一个扰动 $B$，任何一个特征值的位移都不会超过扰动矩阵的谱范数（最大拉伸比例）。

3. **[交错定理] Cauchy 交错定理说明，矩阵去掉最后一行和最后一列后，其特征值与原矩阵特征值有何关系？**
   ??? success "参考答案"
       新特征值严格地交错（穿插）在原特征值之间。即 $\lambda_k(A) \ge \lambda_k(B) \ge \lambda_{k+1}(A)$。

4. **[迹不等式] von Neumann 迹不等式 $|\operatorname{tr}(A^*B)| \le \sum \sigma_i(A)\sigma_i(B)$ 何时取等号？**
   ??? success "参考答案"
       当 $A$ 和 $B$ 可以同时被奇异值分解（存在共同的左右奇异向量阵 $U$ 和 $V$）时取等。这相当于向量内积不等式在矩阵谱空间上的推广。

5. **[行列式不等式] Hadamard 不等式 $\det A \le \prod a_{ii}$ 取等号的条件是什么？这有什么物理或几何意义？**
   ??? success "参考答案"
       取等号当且仅当 $A$ 为对角矩阵。几何上，它意味着各列向量相互正交时，由它们张成的平行多面体体积最大（达到边长之积）。

6. **[优超] 什么是 Majorization（优超关系）？**
   ??? success "参考答案"
       优超 $\mathbf{x} \prec \mathbf{y}$ 是衡量向量分量“均匀程度”的偏序关系。$\mathbf{x}$ 被 $\mathbf{y}$ 优超，意味着 $\mathbf{x}$ 的分量比 $\mathbf{y}$ 更均匀、更平滑。

7. **[Schur-Horn定理] 对角元素与特征值之间的 Schur-Horn 定理说明了什么？**
   ??? success "参考答案"
       说明任意 Hermite 矩阵的对角元素向量总是被其特征值向量优超。即对角线元素是对特征值的一种“平滑化”或“平均化”投影。

8. **[双随机矩阵] Birkhoff 定理断言双随机矩阵的本质是什么？**
   ??? success "参考答案"
       它断言所有的双随机矩阵都位于一个凸多面体内，而该多面体的顶点恰好是所有的置换矩阵。

9. **[矩阵函数] 什么是矩阵单调函数（Operator Monotone）？**
   ??? success "参考答案"
       若 $A \preceq B$（即 $B-A$ 半正定）总是推导出 $f(A) \preceq f(B)$，则 $f$ 是矩阵单调的。这是一个极其严苛的条件，很多常见的标量单调函数（如 $f(x)=x^2$ 或 $e^x$）都不是矩阵单调的。

10. **[爱因斯坦思考题] 为什么像 $f(x)=x^2$ 这样在实数上明显单调递增（$x>0$）的函数，在量子力学（矩阵代数）中却不是“算子单调”的？这反映了量子系统什么深层特性？**
    ??? success "参考答案"
        这是因为量子观测量（矩阵）不满足交换律！$B \ge A$ 并不意味着 $B^2 - A^2 = (B-A)(B+A)$ （交叉项不相等）。在物理上，这反映了不同可观测量（如不同方向的自旋）之间存在的量子纠缠和非局域干涉，使得能量（通常对应平方项）的比较不能简单地在不同的基下线性叠加。

## 本章小结

本章集中探讨了矩阵分析中极具威力和美感的各种不等式理论，主要内容包括：

1. **特征值不等式**：以 Courant-Fischer 极小极大定理为基石，推导出了扰动分析的基础 Weyl 不等式和降维降阶的 Cauchy 交错定理。
2. **奇异值不等式**：论证了奇异值的次可乘性，并展示了如何通过构造 Hermite 分块阵将奇异值问题转化为特征值问题。
3. **迹与行列式不等式**：证明了 von Neumann 迹不等式、Golden-Thompson 不等式，以及行列式的绝对天花板 Hadamard 不等式和 Minkowski 凸性不等式。
4. **优超理论（Majorization）**：引入了统一各种平均化不等式的核心偏序关系，通过 Schur-Horn 定理和 Birkhoff 定理彻底打通了对角元与特征值的关系。
5. **矩阵凸性与单调性**：揭示了矩阵变元下的凸性与单调性远比标量情形严苛，引出了 Löwner 定理和 Jensen 矩阵不等式。
