# 第 10 章 矩阵分解

<div class="context-flow" markdown>

**前置**：Ch8 正交性/谱定理 · Ch9 正定性 · **本章脉络**：LU（消元） → PLU（选主元） → Cholesky（正定） → QR（正交化） → Schur（三角化） → 谱分解 → 极分解
本质：将矩阵拆为"结构简单的因子之积"——每种分解回答不同的问题

</div>

矩阵分解（matrix decomposition / matrix factorization）是将一个矩阵表示为若干特殊矩阵之积的技术，是数值线性代数和应用数学的核心工具。不同的矩阵分解适用于不同的应用场景：LU 分解用于高效求解线性方程组，QR 分解用于最小二乘问题，Schur 分解揭示矩阵的深层结构，谱分解刻画了对称矩阵和正规矩阵的本质。本章系统地介绍各种重要的矩阵分解方法。

---

## 10.1 LU 分解

<div class="context-flow" markdown>

高斯消元的矩阵化表述：$A = LU$（下三角 × 上三角） → 将 $O(n^3)$ 的分解做一次，之后每个右端 $\mathbf{b}$ 只需 $O(n^2)$

</div>

!!! definition "定义 10.1 (LU 分解)"
    设 $A$ 是 $n \times n$ 矩阵。若存在 $n \times n$ 下三角矩阵 $L$（对角元素全为 1）和上三角矩阵 $U$ 使得

    $$A = LU$$

    则称 $A = LU$ 为 $A$ 的 **LU 分解**（LU decomposition / LU factorization）。$L$ 称为**单位下三角矩阵**（unit lower triangular matrix），$U$ 称为**上三角矩阵**（upper triangular matrix）。

!!! theorem "定理 10.1 (LU 分解的存在唯一性)"
    设 $A$ 是 $n \times n$ 矩阵。$A$ 具有 LU 分解当且仅当 $A$ 的所有前 $k$ 阶顺序主子矩阵 $A_k$（$k = 1, \ldots, n-1$）都是非奇异的，即 $\det(A_k) \neq 0$。若 $A$ 可逆，则 LU 分解唯一。

??? proof "证明"
    **必要性：** 设 $A = LU$，$L, U$ 分别为下三角和上三角。$A_k = L_k U_k$，其中 $L_k, U_k$ 分别是 $L, U$ 的前 $k$ 阶主子矩阵。$L_k$ 是单位下三角，$\det(L_k) = 1 \neq 0$。若 $\det(A_k) = 0$，则 $\det(U_k) = 0$，这意味着 $U$ 的前 $k$ 个对角元中有零。但高斯消元过程要求主元非零，矛盾。

    **充分性：** 用高斯消元法的语言。由 $\det(A_k) \neq 0$ 保证消元过程中所有主元非零，不需要行交换。将消元过程记录为 $A = L_1^{-1}L_2^{-1}\cdots L_{n-1}^{-1}U$，其中 $L_i$ 是初等下三角矩阵。$L = L_1^{-1}L_2^{-1}\cdots L_{n-1}^{-1}$ 仍是单位下三角矩阵。

    **唯一性（$A$ 可逆时）：** 设 $A = L_1U_1 = L_2U_2$，则 $L_2^{-1}L_1 = U_2U_1^{-1}$。左边是单位下三角，右边是上三角，因此两边都等于单位矩阵 $I$。$\blacksquare$

!!! example "例 10.1"
    对矩阵 $A = \begin{pmatrix} 2 & 1 & 1 \\ 4 & 3 & 3 \\ 8 & 7 & 9 \end{pmatrix}$ 进行 LU 分解。

    高斯消元：$R_2 \leftarrow R_2 - 2R_1$，$R_3 \leftarrow R_3 - 4R_1$：

    $$\begin{pmatrix} 2 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 3 & 5 \end{pmatrix}$$

    $R_3 \leftarrow R_3 - 3R_2$：

    $$U = \begin{pmatrix} 2 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix}$$

    乘数矩阵 $L = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 4 & 3 & 1 \end{pmatrix}$。

    验证：$LU = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 4 & 3 & 1 \end{pmatrix}\begin{pmatrix} 2 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix} = \begin{pmatrix} 2 & 1 & 1 \\ 4 & 3 & 3 \\ 8 & 7 & 9 \end{pmatrix} = A$。

!!! note "注"
    LU 分解的主要应用是高效求解线性方程组 $A\mathbf{x} = \mathbf{b}$：先解 $L\mathbf{y} = \mathbf{b}$（前代法），再解 $U\mathbf{x} = \mathbf{y}$（回代法）。计算量为 $O(n^2)$（假设 $L, U$ 已知），而分解本身的计算量为 $O(\frac{2}{3}n^3)$。

---

## 10.2 PLU 分解

<div class="context-flow" markdown>

LU 要求顺序主子式非零 → 引入**置换矩阵** $P$（行交换）→ $PA = LU$ 对**任意可逆矩阵**成立

</div>

!!! definition "定义 10.2 (置换矩阵)"
    **置换矩阵**（permutation matrix）是恰好在每行每列有且仅有一个 $1$、其余元素为 $0$ 的方阵。它对应于单位矩阵的行（列）的一个排列。

!!! proposition "命题 10.1 (置换矩阵的性质)"
    设 $P$ 是置换矩阵。则

    1. $P$ 是正交矩阵：$P^TP = PP^T = I$，即 $P^{-1} = P^T$；
    2. $\det(P) = \pm 1$；
    3. 置换矩阵的乘积仍是置换矩阵。

??? proof "证明"
    (1) $P$ 的行（列）向量是标准基向量的排列，因此是标准正交集。

    (2) 由 $\det(P^TP) = (\det P)^2 = 1$ 得 $\det P = \pm 1$。

    (3) 两个排列的复合仍是排列。$\blacksquare$

!!! theorem "定理 10.2 (PLU 分解)"
    任意 $n \times n$ 可逆矩阵 $A$ 都可以分解为

    $$PA = LU$$

    其中 $P$ 是置换矩阵，$L$ 是单位下三角矩阵，$U$ 是上三角矩阵。等价地，

    $$A = P^{-1}LU = P^TLU$$

??? proof "证明"
    在高斯消元过程中，若某步的主元为零，可以通过行交换（部分选主元）使主元非零。所有行交换操作对应一个置换矩阵 $P$。对 $PA$ 进行不需要行交换的高斯消元，得到 $PA = LU$。

    更精确地说：对可逆矩阵 $A$，每步消元时选择当前列中绝对值最大的元素作为主元（部分选主元策略），通过行交换将其移至对角位置。全部行交换的效果由 $P$ 记录。由于 $A$ 可逆，消元过程一定能完成。$\blacksquare$

!!! example "例 10.2"
    对 $A = \begin{pmatrix} 0 & 2 & 1 \\ 1 & 1 & 0 \\ 2 & 0 & 3 \end{pmatrix}$ 进行 PLU 分解。

    第一列主元为 $0$，需要行交换。交换第 1 行和第 3 行：

    $$P_1A = \begin{pmatrix} 2 & 0 & 3 \\ 1 & 1 & 0 \\ 0 & 2 & 1 \end{pmatrix}$$

    $R_2 \leftarrow R_2 - \frac{1}{2}R_1$：

    $$\begin{pmatrix} 2 & 0 & 3 \\ 0 & 1 & -3/2 \\ 0 & 2 & 1 \end{pmatrix}$$

    $R_3 \leftarrow R_3 - 2R_2$：

    $$U = \begin{pmatrix} 2 & 0 & 3 \\ 0 & 1 & -3/2 \\ 0 & 0 & 4 \end{pmatrix}$$

    因此 $P = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}$，$L = \begin{pmatrix} 1 & 0 & 0 \\ 1/2 & 1 & 0 \\ 0 & 2 & 1 \end{pmatrix}$，$PA = LU$。

---

## 10.3 Cholesky 分解

<div class="context-flow" markdown>

Ch9 正定 $A = C^TC$ → 取 $C = L^T$ 得 $A = LL^T$ → 计算量仅 $\frac{1}{3}n^3$（LU 的一半），无需选主元，数值稳定

</div>

!!! definition "定义 10.3 (Cholesky 分解)"
    设 $A$ 是 $n \times n$ 实对称正定矩阵。$A$ 的 **Cholesky 分解**（Cholesky decomposition）是

    $$A = LL^T$$

    其中 $L$ 是对角元素全为正的下三角矩阵。对复 Hermitian 正定矩阵，$A = LL^H$。

!!! theorem "定理 10.3 (Cholesky 分解的存在唯一性)"
    设 $A$ 是 $n \times n$ 实对称正定矩阵。则存在唯一的对角元素全为正的下三角矩阵 $L$ 使得 $A = LL^T$。

??? proof "证明"
    **存在性：** 由 LU 分解（正定矩阵的顺序主子式全为正，保证 LU 分解存在）：$A = L_0 U_0$。令 $D = \operatorname{diag}(u_{11}, \ldots, u_{nn})$（$U_0$ 的对角线），则 $U_0 = DU_1$，其中 $U_1$ 是单位上三角。

    由 $A = A^T$，有 $L_0DU_1 = U_1^TDL_0^T$。由唯一性可推出 $U_1 = L_0^T$，故 $A = L_0DL_0^T$。

    正定性保证 $u_{ii} = d_i > 0$。令 $L = L_0 D^{1/2}$（其中 $D^{1/2} = \operatorname{diag}(\sqrt{d_1}, \ldots, \sqrt{d_n})$），则 $A = L_0 D^{1/2}(D^{1/2})^TL_0^T = LL^T$。

    **唯一性：** 设 $A = L_1L_1^T = L_2L_2^T$。则 $L_2^{-1}L_1 = L_2^T L_1^{-T}$。左边是下三角，右边是上三角，故两边都是对角矩阵 $D$。由 $L_2^{-1}L_1 = D$ 和 $L_2^TL_1^{-T} = D$，得 $D = D^T = D$，且 $DD^T = I$。由对角元素为正，$D = I$，故 $L_1 = L_2$。$\blacksquare$

!!! example "例 10.3"
    对 $A = \begin{pmatrix} 4 & 2 & -2 \\ 2 & 10 & 2 \\ -2 & 2 & 5 \end{pmatrix}$ 进行 Cholesky 分解。

    设 $L = \begin{pmatrix} l_{11} & 0 & 0 \\ l_{21} & l_{22} & 0 \\ l_{31} & l_{32} & l_{33} \end{pmatrix}$，由 $A = LL^T$：

    - $l_{11}^2 = 4 \Rightarrow l_{11} = 2$
    - $l_{21}l_{11} = 2 \Rightarrow l_{21} = 1$
    - $l_{31}l_{11} = -2 \Rightarrow l_{31} = -1$
    - $l_{21}^2 + l_{22}^2 = 10 \Rightarrow l_{22}^2 = 9 \Rightarrow l_{22} = 3$
    - $l_{31}l_{21} + l_{32}l_{22} = 2 \Rightarrow l_{32} = 1$
    - $l_{31}^2 + l_{32}^2 + l_{33}^2 = 5 \Rightarrow l_{33}^2 = 3 \Rightarrow l_{33} = \sqrt{3}$

    $$L = \begin{pmatrix} 2 & 0 & 0 \\ 1 & 3 & 0 \\ -1 & 1 & \sqrt{3} \end{pmatrix}$$

!!! note "注"
    Cholesky 分解的计算量约为 $\frac{1}{3}n^3$，是 LU 分解的一半。它在数值计算中非常稳定，不需要选主元。Cholesky 分解也常用于判定矩阵是否正定：若算法顺利完成（所有对角元素为正），则矩阵正定。

---

## 10.4 QR 分解

<div class="context-flow" markdown>

正交 × 上三角：$A = QR$ → 三种构造法：**Gram-Schmidt**（Ch8 正交化）、**Householder** 反射、**Givens** 旋转 → 数值最小二乘的核心工具

</div>

!!! definition "定义 10.4 (QR 分解)"
    设 $A$ 是 $m \times n$ 矩阵（$m \geq n$）。$A$ 的 **QR 分解**（QR decomposition）是

    $$A = QR$$

    其中 $Q$ 是 $m \times m$ 正交（或酉）矩阵，$R$ 是 $m \times n$ 上三角矩阵。

    **缩减形式**（thin QR）：$A = Q_1R_1$，其中 $Q_1$ 是 $m \times n$（列正交，$Q_1^TQ_1 = I_n$），$R_1$ 是 $n \times n$ 上三角。

!!! theorem "定理 10.4 (QR 分解的存在性)"
    任意 $m \times n$ 矩阵 $A$（$m \geq n$）都具有 QR 分解。若 $A$ 列满秩且要求 $R$ 的对角元素为正，则缩减 QR 分解唯一。

??? proof "证明"
    **存在性（Gram-Schmidt 方法）：** 设 $A = [\mathbf{a}_1 | \cdots | \mathbf{a}_n]$。对列向量 $\mathbf{a}_1, \ldots, \mathbf{a}_n$ 施以 Gram-Schmidt 正交化：

    $$\mathbf{u}_1 = \mathbf{a}_1, \quad \mathbf{e}_1 = \frac{\mathbf{u}_1}{\|\mathbf{u}_1\|}$$

    $$\mathbf{u}_k = \mathbf{a}_k - \sum_{j=1}^{k-1}\langle \mathbf{a}_k, \mathbf{e}_j\rangle \mathbf{e}_j, \quad \mathbf{e}_k = \frac{\mathbf{u}_k}{\|\mathbf{u}_k\|}$$

    则 $\mathbf{a}_k = \sum_{j=1}^{k}\langle \mathbf{a}_k, \mathbf{e}_j\rangle \mathbf{e}_j$，写成矩阵形式即 $A = Q_1R_1$。

    **唯一性：** 若 $A = Q_1R_1 = Q_2R_2$（$R_i$ 对角元素为正），则 $Q_2^TQ_1 = R_2R_1^{-1}$。左边列正交，右边上三角。$Q_2^TQ_1$ 同时是正交矩阵和上三角矩阵，因此是对角矩阵，且对角元素为 $\pm 1$。由 $R_1, R_2$ 对角元素为正，该对角矩阵为 $I$。$\blacksquare$

### Gram-Schmidt 正交化

!!! definition "定义 10.5 (Gram-Schmidt 正交化)"
    设 $\{\mathbf{a}_1, \ldots, \mathbf{a}_n\}$ 是线性无关向量组。**Gram-Schmidt 正交化**（Gram-Schmidt orthogonalization）过程为：

    $$\mathbf{u}_1 = \mathbf{a}_1, \quad \mathbf{e}_1 = \frac{\mathbf{u}_1}{\|\mathbf{u}_1\|}$$

    对 $k = 2, \ldots, n$：

    $$\mathbf{u}_k = \mathbf{a}_k - \sum_{j=1}^{k-1} \frac{\langle \mathbf{a}_k, \mathbf{u}_j \rangle}{\langle \mathbf{u}_j, \mathbf{u}_j \rangle} \mathbf{u}_j, \quad \mathbf{e}_k = \frac{\mathbf{u}_k}{\|\mathbf{u}_k\|}$$

    则 $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ 是标准正交组，且对每个 $k$，$\operatorname{span}\{\mathbf{e}_1, \ldots, \mathbf{e}_k\} = \operatorname{span}\{\mathbf{a}_1, \ldots, \mathbf{a}_k\}$。

### Householder 变换

!!! definition "定义 10.6 (Householder 反射)"
    设 $\mathbf{v} \in \mathbb{R}^n$ 是非零向量。**Householder 反射**（Householder reflection）定义为

    $$H = I - 2\frac{\mathbf{v}\mathbf{v}^T}{\mathbf{v}^T\mathbf{v}}$$

    $H$ 是对称正交矩阵（$H^T = H$，$H^2 = I$），几何上表示关于 $\mathbf{v}$ 的正交补超平面的反射。

!!! proposition "命题 10.2"
    Householder 反射矩阵 $H$ 满足：

    1. $H = H^T$（对称）；
    2. $H^2 = I$（对合）；
    3. $H$ 是正交矩阵，$\det(H) = -1$。

??? proof "证明"
    (1) $(I - 2\frac{\mathbf{v}\mathbf{v}^T}{\mathbf{v}^T\mathbf{v}})^T = I - 2\frac{\mathbf{v}\mathbf{v}^T}{\mathbf{v}^T\mathbf{v}} = H$。

    (2) $H^2 = (I - 2\frac{\mathbf{v}\mathbf{v}^T}{\|\mathbf{v}\|^2})^2 = I - 4\frac{\mathbf{v}\mathbf{v}^T}{\|\mathbf{v}\|^2} + 4\frac{\mathbf{v}(\mathbf{v}^T\mathbf{v})\mathbf{v}^T}{\|\mathbf{v}\|^4} = I$。

    (3) 由 (1) 和 (2) 得 $HH^T = H^2 = I$。$\det(H) = -1$ 因为 $H$ 有一个特征值为 $-1$，其余 $n-1$ 个特征值为 $1$。$\blacksquare$

!!! theorem "定理 10.5 (Householder QR)"
    对任意 $m \times n$ 矩阵 $A$（$m \geq n$），存在 Householder 反射 $H_1, H_2, \ldots, H_n$ 使得

    $$H_n \cdots H_2 H_1 A = R$$

    为上三角矩阵。因此 $A = QR$，其中 $Q = H_1 H_2 \cdots H_n$。

??? proof "证明"
    **构造 $H_1$：** 设 $\mathbf{a}_1$ 是 $A$ 的第一列。选择 $\mathbf{v}_1 = \mathbf{a}_1 + \operatorname{sign}(a_{11})\|\mathbf{a}_1\|\mathbf{e}_1$，令 $H_1 = I - 2\frac{\mathbf{v}_1\mathbf{v}_1^T}{\mathbf{v}_1^T\mathbf{v}_1}$，则 $H_1\mathbf{a}_1 = -\operatorname{sign}(a_{11})\|\mathbf{a}_1\|\mathbf{e}_1$，即第一列被消为只有第一个分量非零。

    对 $H_1A$ 的右下角 $(m-1)\times(n-1)$ 子矩阵重复此过程。将 $(m-1)$ 维的 Householder 矩阵嵌入 $m$ 维（左上角补 $1$）得到 $H_2$。如此迭代 $n$ 步。$\blacksquare$

### Givens 旋转

!!! definition "定义 10.7 (Givens 旋转)"
    **Givens 旋转**（Givens rotation）$G(i, j, \theta)$ 是单位矩阵在第 $(i,i), (i,j), (j,i), (j,j)$ 位置修改为

    $$g_{ii} = g_{jj} = \cos\theta, \quad g_{ij} = -\sin\theta, \quad g_{ji} = \sin\theta$$

    其余位置不变。$G(i,j,\theta)$ 是正交矩阵，表示在 $(i,j)$ 平面内旋转 $\theta$ 角。

!!! example "例 10.4"
    用 Gram-Schmidt 方法对 $A = \begin{pmatrix} 1 & 1 & 0 \\ 1 & 0 & 1 \\ 0 & 1 & 1 \end{pmatrix}$ 进行 QR 分解。

    $\mathbf{a}_1 = (1,1,0)^T$，$\|\mathbf{a}_1\| = \sqrt{2}$，$\mathbf{e}_1 = \frac{1}{\sqrt{2}}(1,1,0)^T$。

    $\mathbf{u}_2 = \mathbf{a}_2 - \langle \mathbf{a}_2, \mathbf{e}_1\rangle \mathbf{e}_1 = (1,0,1)^T - \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}}(1,1,0)^T = (1,0,1)^T - \frac{1}{2}(1,1,0)^T = \frac{1}{2}(1,-1,2)^T$。

    $\|\mathbf{u}_2\| = \frac{1}{2}\sqrt{6}$，$\mathbf{e}_2 = \frac{1}{\sqrt{6}}(1,-1,2)^T$。

    $\mathbf{u}_3 = \mathbf{a}_3 - \langle \mathbf{a}_3, \mathbf{e}_1\rangle \mathbf{e}_1 - \langle \mathbf{a}_3, \mathbf{e}_2\rangle \mathbf{e}_2$。

    $\langle \mathbf{a}_3, \mathbf{e}_1\rangle = \frac{1}{\sqrt{2}}$，$\langle \mathbf{a}_3, \mathbf{e}_2\rangle = \frac{1}{\sqrt{6}}(0-1+2) = \frac{1}{\sqrt{6}}$。

    $\mathbf{u}_3 = (0,1,1)^T - \frac{1}{2}(1,1,0)^T - \frac{1}{6}(1,-1,2)^T = \frac{1}{3}(-1,1,1)^T \cdot \frac{3}{2}$。

    经计算，$\mathbf{u}_3 = (-\frac{2}{3}, \frac{2}{3}, \frac{2}{3})^T$，$\|\mathbf{u}_3\| = \frac{2}{\sqrt{3}}$，$\mathbf{e}_3 = \frac{1}{\sqrt{3}}(-1, 1, 1)^T$。

    $$Q = \begin{pmatrix} 1/\sqrt{2} & 1/\sqrt{6} & -1/\sqrt{3} \\ 1/\sqrt{2} & -1/\sqrt{6} & 1/\sqrt{3} \\ 0 & 2/\sqrt{6} & 1/\sqrt{3} \end{pmatrix}$$

    $$R = Q^TA = \begin{pmatrix} \sqrt{2} & 1/\sqrt{2} & 1/\sqrt{2} \\ 0 & \sqrt{6}/2 & 1/\sqrt{6} \\ 0 & 0 & 2/\sqrt{3} \end{pmatrix}$$

---

## 10.5 Schur 分解

<div class="context-flow" markdown>

任意方阵可**酉三角化** $A = UTU^H$ → 对角元 = 特征值 → 比对角化更一般（不要求可对角化） → 正规矩阵时 $T$ 退化为对角

</div>

<div class="context-flow" markdown>

**洞察**：Schur 分解的存在性证明与谱定理同构——找特征向量、正交补不变、归纳降维——但不需正规性假设

</div>

!!! theorem "定理 10.6 (Schur 分解)"
    设 $A$ 是 $n \times n$ 复矩阵。则存在酉矩阵 $U$ 使得

    $$U^H A U = T$$

    其中 $T$ 是上三角矩阵，其对角元素是 $A$ 的特征值（按任意指定顺序排列）。等价地，$A = UTU^H$。

??? proof "证明"
    对 $n$ 用数学归纳法。$n = 1$ 时显然。设结论对 $n-1$ 阶矩阵成立。

    由代数基本定理，$A$ 有特征值 $\lambda_1$。设 $\mathbf{u}_1$ 是对应的单位特征向量。将 $\mathbf{u}_1$ 扩充为 $\mathbb{C}^n$ 的标准正交基 $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$，令 $U_1 = [\mathbf{u}_1 | \cdots | \mathbf{u}_n]$。则

    $$U_1^H A U_1 = \begin{pmatrix} \lambda_1 & \mathbf{b}^H \\ \mathbf{0} & A_1 \end{pmatrix}$$

    其中 $A_1$ 是 $(n-1) \times (n-1)$ 矩阵（第一列为零向量是因为 $\mathbf{u}_j^HA\mathbf{u}_1 = \lambda_1 \mathbf{u}_j^H\mathbf{u}_1 = 0$，$j \geq 2$）。

    由归纳假设，存在 $(n-1)$ 阶酉矩阵 $V_1$ 使 $V_1^HA_1V_1 = T_1$ 为上三角。令

    $$V = \begin{pmatrix} 1 & \mathbf{0}^H \\ \mathbf{0} & V_1 \end{pmatrix}$$

    则 $U = U_1V$ 是酉矩阵，且 $U^HAU$ 是上三角矩阵，对角元素为 $A$ 的特征值。$\blacksquare$

!!! theorem "定理 10.7 (实 Schur 分解)"
    设 $A$ 是 $n \times n$ 实矩阵。则存在正交矩阵 $Q$ 使得

    $$Q^TAQ = T$$

    其中 $T$ 是准上三角矩阵（quasi-upper triangular），即对角块是 $1\times 1$ 块（对应实特征值）或 $2\times 2$ 块 $\begin{pmatrix} a & b \\ c & d \end{pmatrix}$（对应共轭复特征值对）。

??? proof "证明"
    思路类似复 Schur 分解。对实特征值，可像复情形一样处理。对共轭复特征值对 $\lambda = a \pm bi$，对应两个共轭复特征向量 $\mathbf{v} \pm i\mathbf{w}$，其实部和虚部 $\mathbf{v}, \mathbf{w}$ 张成一个 $2$ 维实不变子空间，在此子空间上 $A$ 的限制给出 $2\times 2$ 块。$\blacksquare$

!!! example "例 10.5"
    对 $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$（旋转 90 度矩阵），其特征值为 $\pm i$。

    复 Schur 分解：$U = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ i & -i \end{pmatrix}$，$U^HAU = \begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}$。

    实 Schur 分解：$A$ 本身已经是准上三角形式（一个 $2\times 2$ 块）：$Q = I$，$T = A$。

---

## 10.6 谱分解

<div class="context-flow" markdown>

$A = \sum \lambda_i P_i$（特征值 × 正交投影之和） → 将矩阵拆为**独立模式的叠加** → Ch11 SVD 是非对称版本，Ch13 $f(A) = \sum f(\lambda_i)P_i$

</div>

!!! definition "定义 10.8 (谱分解)"
    设 $A$ 是可对角化的矩阵。若 $A$ 有特征值 $\lambda_1, \ldots, \lambda_k$（互不相同），对应的特征空间投影为 $P_1, \ldots, P_k$，则

    $$A = \lambda_1 P_1 + \lambda_2 P_2 + \cdots + \lambda_k P_k$$

    称为 $A$ 的**谱分解**（spectral decomposition）。

<div class="context-flow" markdown>

**洞察**：$A = \sum \lambda_i \mathbf{q}_i\mathbf{q}_i^T$ 中每个 $\mathbf{q}_i\mathbf{q}_i^T$ 是秩一正交投影——矩阵被分解为"频率分量"，这是 Fourier 分析的有限维版本

</div>

!!! theorem "定理 10.8 (对称矩阵的谱分解)"
    设 $A$ 是 $n \times n$ 实对称矩阵，特征值为 $\lambda_1, \ldots, \lambda_n$（允许重复），对应的标准正交特征向量为 $\mathbf{q}_1, \ldots, \mathbf{q}_n$。则

    $$A = Q\Lambda Q^T = \sum_{i=1}^n \lambda_i \mathbf{q}_i \mathbf{q}_i^T$$

    若将相同特征值的项合并，设不同特征值为 $\mu_1, \ldots, \mu_k$，对应的特征空间正交投影为 $P_j = \sum_{\lambda_i = \mu_j} \mathbf{q}_i\mathbf{q}_i^T$，则

    $$A = \sum_{j=1}^k \mu_j P_j, \quad \text{其中} \quad \sum_{j=1}^k P_j = I, \quad P_j P_l = \delta_{jl} P_j$$

??? proof "证明"
    由谱定理，$A = Q\Lambda Q^T$。展开矩阵乘法：

    $$A = \sum_{i=1}^n \lambda_i \mathbf{q}_i \mathbf{q}_i^T$$

    验证 $P_j$ 的性质：$P_j$ 是到特征空间 $E_{\mu_j}$ 上的正交投影，因此 $P_j^2 = P_j$，$P_j^T = P_j$，$P_jP_l = 0$（$j \neq l$，因不同特征空间正交），$\sum P_j = I$（因特征向量组成完整的标准正交基）。$\blacksquare$

!!! theorem "定理 10.9 (正规矩阵的谱分解)"
    设 $A$ 是 $n \times n$ 正规矩阵，特征值为 $\lambda_1, \ldots, \lambda_n$，对应的标准正交特征向量为 $\mathbf{u}_1, \ldots, \mathbf{u}_n$。则

    $$A = U\Lambda U^H = \sum_{i=1}^n \lambda_i \mathbf{u}_i \mathbf{u}_i^H$$

??? proof "证明"
    由正规矩阵的谱定理（定理 8.17），$A$ 可酉对角化为 $A = U\Lambda U^H$。展开即得。$\blacksquare$

!!! example "例 10.6"
    对 $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ 进行谱分解。

    特征值 $\lambda_1 = 1$，$\mathbf{q}_1 = \frac{1}{\sqrt{2}}(1,-1)^T$；$\lambda_2 = 3$，$\mathbf{q}_2 = \frac{1}{\sqrt{2}}(1,1)^T$。

    $$P_1 = \mathbf{q}_1\mathbf{q}_1^T = \frac{1}{2}\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}, \quad P_2 = \mathbf{q}_2\mathbf{q}_2^T = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$$

    $$A = 1 \cdot P_1 + 3 \cdot P_2 = \frac{1}{2}\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix} + \frac{3}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$$

    验证：$P_1 + P_2 = I$，$P_1P_2 = 0$。

---

## 10.7 极分解

<div class="context-flow" markdown>

$A = UP$：酉 × 正定 ↔ 复数极坐标 $z = e^{i\theta}|z|$ → 通过 SVD 计算：$U = VW^H$, $P = W\Sigma W^H$ → 分离"旋转"与"伸缩"

</div>

!!! definition "定义 10.9 (极分解)"
    设 $A$ 是 $n \times n$ 可逆复矩阵。$A$ 的**极分解**（polar decomposition）是

    $$A = UP$$

    其中 $U$ 是酉矩阵（$U^HU = I$），$P$ 是正定 Hermitian 矩阵（$P^H = P$，$P > 0$）。这称为**右极分解**。类似地，$A = P'U$（$P'$ 正定）称为**左极分解**。

!!! theorem "定理 10.10 (极分解的存在唯一性)"
    (1) 任意 $n \times n$ 可逆矩阵 $A$ 都有唯一的右极分解 $A = UP$（$U$ 酉，$P$ 正定 Hermitian）和唯一的左极分解 $A = P'U$。

    (2) 更一般地，任意 $m \times n$ 矩阵 $A$ 都可以分解为 $A = UP$，其中 $U$ 的列是正交的（$U^HU = I_n$），$P$ 是半正定 Hermitian 矩阵。$P$ 唯一确定（$P = \sqrt{A^HA}$），但 $U$ 在 $A$ 不满秩时不唯一。

??? proof "证明"
    **存在性：** $A^HA$ 是半正定 Hermitian 矩阵。令 $P = (A^HA)^{1/2}$（正半定 Hermitian 平方根，由谱定理可构造）。

    若 $A$ 可逆，则 $A^HA$ 正定，$P$ 正定可逆。令 $U = AP^{-1}$。验证 $U$ 是酉矩阵：

    $$U^HU = (AP^{-1})^H(AP^{-1}) = P^{-H}A^HAP^{-1} = P^{-1}P^2P^{-1} = I$$

    **唯一性（可逆情形）：** 设 $A = U_1P_1 = U_2P_2$，则 $A^HA = P_1^2 = P_2^2$。由正定矩阵正平方根的唯一性（见定理 13.7），$P_1 = P_2$。从而 $U_1 = U_2$。

    **左极分解：** $P' = (AA^H)^{1/2}$，$U = (P')^{-1}A$。$\blacksquare$

!!! note "注"
    极分解的名称来自与复数极坐标表示 $z = e^{i\theta}|z|$ 的类比：$U$ 对应"旋转"（$e^{i\theta}$），$P$ 对应"伸缩"（$|z|$）。极分解可以通过 SVD 方便地计算：若 $A = V\Sigma W^H$（SVD），则 $U = VW^H$，$P = W\Sigma W^H$。

!!! example "例 10.7"
    对 $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ 进行极分解。

    $A^TA = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$，特征值为 $\frac{3 \pm \sqrt{5}}{2}$。

    设 $\lambda_1 = \frac{3+\sqrt{5}}{2}$，$\lambda_2 = \frac{3-\sqrt{5}}{2}$。$P = (A^TA)^{1/2}$，即 $P$ 的特征值为 $\sqrt{\lambda_1}, \sqrt{\lambda_2}$。

    $P = Q\operatorname{diag}(\sqrt{\lambda_1}, \sqrt{\lambda_2})Q^T$（$Q$ 由 $A^TA$ 的标准正交特征向量组成）。

    $U = AP^{-1}$ 是正交矩阵。

    具体数值：$P \approx \begin{pmatrix} 1.0515 & 0.4851 \\ 0.4851 & 1.2965 \end{pmatrix}$，$U \approx \begin{pmatrix} 0.8507 & 0.5257 \\ -0.5257 & 0.8507 \end{pmatrix}$。
