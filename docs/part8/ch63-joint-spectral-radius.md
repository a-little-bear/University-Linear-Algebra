# 第 63 章 联合谱半径与同时三角化

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵范数(Ch15) · 矩阵分析(Ch14)

**本章脉络**：同时三角化条件 → McCoy 定理 → 联合谱半径定义 → Berger-Wang 定理 → 广义谱半径 → 可计算性 → 切换系统稳定性 → 小波正则性

**延伸**：联合谱半径在自动控制（切换线性系统的稳定性判定）、小波理论（细分格式的收敛性和光滑性判定）和编码理论（联合谱半径与码的容量）中有核心应用

</div>

当我们从单一矩阵的谱理论走向**矩阵集合**的谱理论时，许多经典概念都获得了深刻的推广。单一矩阵的谱半径 $\rho(A)$ 刻画了 $A^k$ 的增长行为；而一组矩阵 $\Sigma = \{A_1, A_2, \ldots, A_m\}$ 的**联合谱半径** $\rho(\Sigma)$ 则刻画了所有可能的矩阵乘积 $A_{i_1} A_{i_2} \cdots A_{i_k}$ 的最大增长率。这一概念由 Rota 和 Strang 于 1960 年引入，此后在动力系统、小波理论和自动控制中获得了广泛而深刻的应用。

在研究联合谱半径之前，我们首先需要理解矩阵集合的基本结构性质——同时三角化和同时对角化。这些经典结果不仅是联合谱分析的基础，也在表示论和 Lie 代数理论中占据核心地位。

本章从同时三角化的 McCoy 定理出发，经由联合谱半径的定义和 Berger-Wang 定理，到达计算方法与应用。

---

## 63.1 同时三角化

<div class="context-flow" markdown>

**核心问题**：什么条件下，一组矩阵可以被同一个可逆矩阵同时化为上三角形？

</div>

### 基本概念

对于单个矩阵，Schur 分解告诉我们：在代数闭域上，任意矩阵都可以酉三角化。但如果我们有**一组**矩阵，能否找到**一个**公共的可逆矩阵，同时将它们化为上三角形？

!!! definition "定义 63.1 (同时三角化)"
    设 $\mathcal{F} = \{A_1, A_2, \ldots, A_m\} \subseteq M_n(\mathbb{F})$ 是一组 $n \times n$ 矩阵。称 $\mathcal{F}$ 可以**同时三角化**（simultaneously triangularizable），若存在可逆矩阵 $P \in M_n(\mathbb{F})$，使得所有 $P^{-1} A_i P$（$i = 1, \ldots, m$）均为上三角矩阵。

    等价地，存在 $\mathbb{F}^n$ 的一个完全旗（complete flag）
    $$
    \{0\} = V_0 \subset V_1 \subset V_2 \subset \cdots \subset V_n = \mathbb{F}^n, \quad \dim V_k = k
    $$
    使得每个 $V_k$ 都是 $\mathcal{F}$ 中所有矩阵的不变子空间。

!!! note "注"
    同时三角化要求所有矩阵共享**同一个**完全旗的不变子空间链，这比每个矩阵单独可三角化强得多。

### McCoy 定理

McCoy 定理给出了代数闭域上同时三角化的完整刻画。

!!! theorem "定理 63.1 (McCoy 定理)"
    设 $\mathbb{F}$ 是代数闭域，$A_1, A_2, \ldots, A_m \in M_n(\mathbb{F})$。则 $\{A_1, \ldots, A_m\}$ 可同时三角化，当且仅当对任意非交换多项式 $p(x_1, \ldots, x_m)$，矩阵 $p(A_1, \ldots, A_m)$ 是幂零的蕴含 $p$ 在所有上三角矩阵上取值也为幂零的。

    更实用的等价条件为：$\{A_1, \ldots, A_m\}$ 可同时三角化，当且仅当对任意非交换多项式 $p(x_1, \ldots, x_m)$，矩阵
    $$
    [A_i, p(A_1, \ldots, A_m)]
    $$
    对所有 $i$ 都是幂零的（其中 $[X, Y] = XY - YX$ 是换位子）。

??? proof "证明"
    **必要性**：若所有 $A_i$ 可同时三角化为上三角矩阵 $T_i = P^{-1}A_iP$，则 $P^{-1}p(A_1, \ldots, A_m)P = p(T_1, \ldots, T_m)$。上三角矩阵的任意多项式仍为上三角矩阵，而两个上三角矩阵的换位子 $[T_i, p(T_1, \ldots, T_m)]$ 是严格上三角的（对角线为零），因此是幂零的。

    **充分性**：使用归纳法对 $n$ 进行归纳。$n = 1$ 时显然成立。

    设 $n \geq 2$，假设命题对所有维度小于 $n$ 的情形成立。关键步骤是证明存在一个公共的一维不变子空间（即公共特征向量）。

    取 $A_1$ 的某个特征值 $\lambda$，设 $V_\lambda = \ker(A_1 - \lambda I)$ 是对应的特征空间。我们需要证明 $V_\lambda$ 中存在其他 $A_i$ 的公共特征向量。

    条件保证了 $[A_i, A_1 - \lambda I]$ 限制在 $V_\lambda$ 上的行为是"良好"的。利用换位子的幂零性和归纳假设，可以在 $V_\lambda$ 中找到所有 $A_i$ 的公共不变一维子空间。

    找到公共特征向量 $v_1$ 后，选取将 $v_1$ 作为第一个基向量的基，所有 $A_i$ 在此基下具有分块形式
    $$
    A_i = \begin{pmatrix} \lambda_i & * \\ 0 & A_i' \end{pmatrix}
    $$
    其中 $A_i' \in M_{n-1}(\mathbb{F})$ 继承了原始条件。由归纳假设完成证明。

### 可交换矩阵的同时三角化

最重要的特殊情形是可交换矩阵族。

!!! theorem "定理 63.2 (可交换矩阵的同时三角化)"
    设 $\mathbb{F}$ 是代数闭域。若 $A_1, A_2, \ldots, A_m \in M_n(\mathbb{F})$ 两两可交换（即 $A_i A_j = A_j A_i$ 对所有 $i, j$ 成立），则它们可以同时三角化。

??? proof "证明"
    由 McCoy 定理，只需验证：对所有 $i$ 和所有非交换多项式 $p$，$[A_i, p(A_1, \ldots, A_m)]$ 是幂零的。

    但由于 $A_i$ 两两可交换，任意多项式 $p(A_1, \ldots, A_m)$ 的计算与变量顺序无关，且 $A_i$ 与 $p(A_1, \ldots, A_m)$ 可交换，从而 $[A_i, p(A_1, \ldots, A_m)] = 0$，当然是幂零的。

    也可以直接证明：对 $n$ 归纳。因为 $A_1, \ldots, A_m$ 两两可交换且 $\mathbb{F}$ 代数闭，$A_1$ 的每个特征空间都是 $A_2, \ldots, A_m$ 的不变子空间。在 $A_1$ 的某一特征空间中，限制 $A_2, \ldots, A_m$ 仍然两两可交换，由归纳法即得公共特征向量。

!!! example "例 63.1"
    设
    $$
    A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}, \quad
    B = \begin{pmatrix} 1 & -1 \\ 0 & 4 \end{pmatrix}
    $$

    验证 $AB = BA$：
    $$
    AB = \begin{pmatrix} 2 & 2 \\ 0 & 12 \end{pmatrix}, \quad
    BA = \begin{pmatrix} 2 & -2 \\ 0 & 12 \end{pmatrix}
    $$
    $AB \neq BA$，所以 $A, B$ 不可交换。但它们已经是上三角的，因此确实可以同时三角化（以 $P = I$ 为变换矩阵）。

    这说明同时三角化的条件比可交换性弱。

### Lie 定理

同时三角化的一个深刻推广来自 Lie 代数理论。

!!! theorem "定理 63.3 (Lie 定理)"
    设 $\mathbb{F}$ 是特征零的代数闭域，$\mathfrak{g} \subseteq \mathfrak{gl}_n(\mathbb{F})$ 是一个可解（solvable）Lie 代数。则 $\mathfrak{g}$ 中的所有矩阵可以同时三角化。

    特别地，若 $\{A_1, \ldots, A_m\}$ 生成的 Lie 代数是可解的，则 $\{A_1, \ldots, A_m\}$ 可同时三角化。

!!! note "注"
    Lie 代数 $\mathfrak{g}$ 是**可解的**，若其导列（derived series）$\mathfrak{g}^{(0)} = \mathfrak{g}$，$\mathfrak{g}^{(k+1)} = [\mathfrak{g}^{(k)}, \mathfrak{g}^{(k)}]$ 最终到达零。两两可交换的矩阵族生成的 Lie 代数是 Abel 的，当然可解。因此 Lie 定理推广了定理 63.2。

---

## 63.2 同时对角化

<div class="context-flow" markdown>

**核心问题**：什么条件下，一组矩阵可以被同一个可逆矩阵同时对角化？

</div>

### 同时对角化的充要条件

!!! definition "定义 63.2 (同时对角化)"
    设 $A_1, A_2, \ldots, A_m \in M_n(\mathbb{F})$。称它们可**同时对角化**（simultaneously diagonalizable），若存在可逆矩阵 $P$ 使得 $P^{-1}A_iP$ 对所有 $i$ 都是对角矩阵。

!!! theorem "定理 63.4 (同时对角化的充要条件)"
    矩阵 $A_1, A_2, \ldots, A_m \in M_n(\mathbb{C})$ 可同时对角化，当且仅当以下两个条件同时成立：

    1. 每个 $A_i$ 都是可对角化的（即每个 $A_i$ 的最小多项式无重根）；
    2. 它们两两可交换：$A_i A_j = A_j A_i$ 对所有 $i, j$ 成立。

??? proof "证明"
    **必要性**：若存在 $P$ 使得 $D_i = P^{-1}A_iP$ 均为对角矩阵，则每个 $A_i$ 可对角化（因为 $A_i = P D_i P^{-1}$），且
    $$
    A_i A_j = P D_i P^{-1} \cdot P D_j P^{-1} = P D_i D_j P^{-1} = P D_j D_i P^{-1} = A_j A_i
    $$
    因为对角矩阵之间总是可交换的。

    **充分性**：对 $m$ 归纳。$m = 1$ 时就是单个矩阵的对角化。

    设 $m \geq 2$，假设 $A_1, \ldots, A_{m-1}$ 已经可以同时对角化（由归纳假设）。设 $P$ 使得 $P^{-1}A_iP = D_i$ 为对角矩阵（$i = 1, \ldots, m-1$）。

    考虑 $B = P^{-1}A_m P$。由于 $A_m$ 与每个 $A_i$ 可交换，$B$ 与每个 $D_i$ 可交换。

    设 $D_1$ 的不同对角元素为 $\lambda_1, \ldots, \lambda_r$。由 $BD_1 = D_1B$，$B$ 必须保持 $D_1$ 的每个特征空间不变。因此 $B$ 具有分块对角形式
    $$
    B = \operatorname{diag}(B_1, B_2, \ldots, B_r)
    $$
    其中 $B_k$ 是 $B$ 在 $D_1$ 的特征值 $\lambda_k$ 对应的特征空间上的限制。

    由于 $A_m$ 可对角化，每个 $B_k$ 也可对角化。设 $Q_k$ 对角化 $B_k$，令 $Q = \operatorname{diag}(Q_1, \ldots, Q_r)$。则 $Q^{-1}BQ$ 是对角矩阵，且 $Q^{-1}D_iQ = D_i$（因为 $Q$ 保持 $D_1$ 的特征空间不变，从而保持 $D_i$ 在这些子空间上的分块结构）。

    因此 $(PQ)^{-1}A_i(PQ)$ 对所有 $i = 1, \ldots, m$ 都是对角矩阵。

!!! example "例 63.2"
    设
    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}, \quad
    B = \begin{pmatrix} 3 & 0 \\ 0 & 5 \end{pmatrix}
    $$
    $A, B$ 都已经是对角矩阵且可交换。它们由 $P = I$ 同时对角化。

!!! example "例 63.3"
    设
    $$
    A = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}, \quad
    B = \begin{pmatrix} 3 & -2 \\ -2 & 3 \end{pmatrix}
    $$
    验证 $AB = BA$（均为实对称矩阵，可交换性由谱定理保证）。共同的特征向量为 $v_1 = (1,1)^T/\sqrt{2}$（$Av_1 = 9v_1, Bv_1 = v_1$）和 $v_2 = (1,-1)^T/\sqrt{2}$（$Av_2 = v_2, Bv_2 = 5v_2$）。

    取 $P = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$，则
    $$
    P^T A P = \begin{pmatrix} 9 & 0 \\ 0 & 1 \end{pmatrix}, \quad
    P^T B P = \begin{pmatrix} 1 & 0 \\ 0 & 5 \end{pmatrix}
    $$

### 正规矩阵的同时对角化

!!! theorem "定理 63.5 (正规矩阵的同时酉对角化)"
    可交换的正规矩阵族可以同时**酉**对角化。即若 $A_1, \ldots, A_m$ 都是正规的（$A_i^* A_i = A_i A_i^*$）且两两可交换，则存在酉矩阵 $U$ 使得 $U^* A_i U$ 对所有 $i$ 都是对角矩阵。

!!! note "注"
    实对称矩阵当然是正规的，因此可交换的实对称矩阵族可以被同一个正交矩阵同时正交对角化。这一结论在量子力学中尤为重要：可交换的可观测量（Hermite 算子）具有公共的本征态基底。

---

## 63.3 联合谱半径

<div class="context-flow" markdown>

**核心问题**：如何量化一组矩阵的所有可能乘积的最大增长率？

</div>

### 定义与动机

对于单个矩阵 $A$，谱半径 $\rho(A) = \max |\lambda_i(A)|$ 通过 Gelfand 公式等于 $\lim_{k\to\infty} \|A^k\|^{1/k}$。当我们有一组矩阵时，需要考虑所有可能的乘积序列。

!!! definition "定义 63.3 (联合谱半径)"
    设 $\Sigma = \{A_1, A_2, \ldots, A_m\} \subset M_n(\mathbb{C})$ 是一个有限矩阵集合，$\|\cdot\|$ 是 $M_n(\mathbb{C})$ 上的某个矩阵范数。$\Sigma$ 的**联合谱半径**（joint spectral radius, JSR）定义为
    $$
    \rho(\Sigma) = \lim_{k \to \infty} \max_{(i_1, \ldots, i_k) \in \{1,\ldots,m\}^k} \|A_{i_1} A_{i_2} \cdots A_{i_k}\|^{1/k}
    $$

    等价地，定义 $\Sigma^k = \{A_{i_1} \cdots A_{i_k} : i_j \in \{1,\ldots,m\}\}$ 为所有长度 $k$ 的乘积的集合，则
    $$
    \rho(\Sigma) = \lim_{k \to \infty} \left(\max_{A \in \Sigma^k} \|A\|\right)^{1/k}
    $$

!!! theorem "定理 63.6 (联合谱半径的良定义性)"
    上述定义中的极限存在，且与所选矩阵范数无关。

??? proof "证明"
    **极限存在**：设 $\hat{\rho}_k = \max_{A \in \Sigma^k} \|A\|^{1/k}$。对于 $k = p + q$，我们有
    $$
    \max_{A \in \Sigma^{p+q}} \|A\| \leq \max_{B \in \Sigma^p} \|B\| \cdot \max_{C \in \Sigma^q} \|C\|
    $$
    因此 $a_k = \log \max_{A \in \Sigma^k} \|A\|$ 是次可加序列（$a_{p+q} \leq a_p + a_q$）。由 Fekete 引理，$\lim_{k\to\infty} a_k/k = \inf_{k} a_k/k$ 存在。从而 $\lim_{k\to\infty} \hat{\rho}_k$ 存在。

    **范数无关性**：设 $\|\cdot\|_a$ 和 $\|\cdot\|_b$ 是两个矩阵范数。由范数等价性，存在常数 $c, C > 0$ 使得 $c\|A\|_a \leq \|A\|_b \leq C\|A\|_a$ 对所有 $A$ 成立。则
    $$
    c^{1/k} \hat{\rho}_k^{(a)} \leq \hat{\rho}_k^{(b)} \leq C^{1/k} \hat{\rho}_k^{(a)}
    $$
    令 $k \to \infty$，因 $c^{1/k}, C^{1/k} \to 1$，两个范数给出相同的极限。

!!! example "例 63.4"
    设 $\Sigma = \{A_1, A_2\}$，其中
    $$
    A_1 = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}, \quad
    A_2 = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}
    $$

    $\rho(A_1) = \rho(A_2) = 1$（均为 Jordan 块），但乘积
    $$
    A_1 A_2 = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}
    $$
    的特征值为 $\frac{3 \pm \sqrt{5}}{2}$，其谱半径 $\rho(A_1 A_2) = \frac{3+\sqrt{5}}{2} \approx 2.618 > 1$。

    因此 $\rho(\Sigma) \geq \rho(A_1 A_2)^{1/2} \approx 1.618$，即黄金比例 $\phi$。实际上 $\rho(\Sigma) = \phi$。

### 基本性质

!!! theorem "定理 63.7 (联合谱半径的基本性质)"
    设 $\Sigma, \Sigma' \subset M_n(\mathbb{C})$ 为有限矩阵集合。则：

    1. **非负性**：$\rho(\Sigma) \geq 0$，且 $\rho(\Sigma) = 0$ 当且仅当 $\Sigma$ 中所有矩阵都是幂零的且满足有界长度乘积性质；
    2. **齐次性**：$\rho(\alpha \Sigma) = |\alpha| \rho(\Sigma)$，其中 $\alpha\Sigma = \{\alpha A : A \in \Sigma\}$；
    3. **次可乘性**：$\rho(\Sigma \cdot \Sigma') \leq \rho(\Sigma) \cdot \rho(\Sigma')$；
    4. **单调性**：若 $\Sigma \subseteq \Sigma'$，则 $\rho(\Sigma) \leq \rho(\Sigma')$；
    5. **连续性**：$\rho(\cdot)$ 作为 $M_n(\mathbb{C})^m$ 到 $\mathbb{R}$ 的函数是连续的；
    6. **退化**：当 $\Sigma = \{A\}$ 时，$\rho(\Sigma) = \rho(A)$（退化为经典谱半径）。

---

## 63.4 Berger-Wang 定理

<div class="context-flow" markdown>

**核心问题**：联合谱半径能否仅用特征值来刻画，而不依赖范数？

</div>

### 广义谱半径

!!! definition "定义 63.4 (广义谱半径)"
    设 $\Sigma = \{A_1, \ldots, A_m\} \subset M_n(\mathbb{C})$。$\Sigma$ 的**广义谱半径**（generalized spectral radius）定义为
    $$
    \hat{\rho}(\Sigma) = \limsup_{k \to \infty} \max_{A \in \Sigma^k} \rho(A)^{1/k}
    $$
    其中 $\rho(A)$ 是矩阵 $A$ 的（经典）谱半径。

!!! note "注"
    广义谱半径 $\hat{\rho}(\Sigma)$ 完全由矩阵乘积的特征值决定，不涉及任何范数的选取。然而，其定义中使用了 $\limsup$ 而非 $\lim$，因为极限不一定存在。

### Berger-Wang 定理

!!! theorem "定理 63.8 (Berger-Wang 定理, 1992)"
    设 $\Sigma \subset M_n(\mathbb{C})$ 是有界矩阵集合。则联合谱半径等于广义谱半径：
    $$
    \rho(\Sigma) = \hat{\rho}(\Sigma) = \limsup_{k \to \infty} \max_{A \in \Sigma^k} \rho(A)^{1/k}
    $$

??? proof "证明"
    **不等式 $\hat{\rho}(\Sigma) \leq \rho(\Sigma)$**：对任意 $A \in \Sigma^k$，
    $$
    \rho(A) \leq \|A\| \leq \max_{B \in \Sigma^k} \|B\|
    $$
    因此 $\max_{A \in \Sigma^k} \rho(A)^{1/k} \leq \left(\max_{B \in \Sigma^k} \|B\|\right)^{1/k}$。取 $\limsup$ 即得 $\hat{\rho}(\Sigma) \leq \rho(\Sigma)$。

    **不等式 $\rho(\Sigma) \leq \hat{\rho}(\Sigma)$**：这是定理的深刻之处。证明需要用到 Gelfand 公式的一个推广。

    设 $\hat{\rho} = \hat{\rho}(\Sigma)$。对任意 $\varepsilon > 0$，需要证明 $\rho(\Sigma) \leq \hat{\rho} + \varepsilon$。

    由 $\hat{\rho}$ 的定义，存在 $k_0$ 使得对所有 $k \geq k_0$ 和所有 $A \in \Sigma^k$，$\rho(A) \leq (\hat{\rho} + \varepsilon/2)^k$。

    利用矩阵范数与谱半径的关系（Gelfand 公式的推广），可以构造一个适当的范数 $\|\cdot\|_\varepsilon$，使得 $\|A\|_\varepsilon \leq (\hat{\rho} + \varepsilon)$ 对所有 $A \in \Sigma$ 成立。由此 $\rho(\Sigma) \leq \hat{\rho} + \varepsilon$。由 $\varepsilon$ 的任意性得到 $\rho(\Sigma) \leq \hat{\rho}(\Sigma)$。

    完整证明需要 Berger 和 Wang 原始论文中更精细的分析，涉及矩阵半群的结构理论。

!!! example "例 63.5"
    回到例 63.4 中的 $\Sigma = \{A_1, A_2\}$。Berger-Wang 定理保证
    $$
    \rho(\Sigma) = \limsup_{k\to\infty} \max_{A \in \Sigma^k} \rho(A)^{1/k}
    $$

    由于 $\rho(A_1 A_2) = \phi^2$，我们有 $\rho(\Sigma) \geq \phi^{2 \cdot 1/2} = \phi$。进一步的分析（利用所有长度 $k$ 的乘积的特征值）可以证明 $\rho(\Sigma) = \phi$。

### 有限性猜想

!!! definition "定义 63.5 (有限性猜想)"
    **有限性猜想**（Finiteness Conjecture）断言：对任意有限矩阵集 $\Sigma$，存在有限长度的乘积 $A_{i_1} \cdots A_{i_k}$ 使得
    $$
    \rho(\Sigma) = \rho(A_{i_1} \cdots A_{i_k})^{1/k}
    $$

    即联合谱半径可由某个有限长度的乘积精确实现。

!!! note "注"
    有限性猜想曾被认为是正确的，但 Bousch-Mairesse（2002）和 Blondel-Theys-Vladimirov（2003）构造了反例。然而，对于许多重要的特殊情形（如非负矩阵对），有限性猜想是成立的。

---

## 63.5 计算方法与复杂性

<div class="context-flow" markdown>

**核心问题**：给定一组矩阵，如何有效地计算或近似其联合谱半径？

</div>

### 不可判定性

!!! theorem "定理 63.9 (JSR 的不可判定性)"
    以下问题是算法不可判定的：给定一组整数矩阵 $\Sigma = \{A_1, \ldots, A_m\} \subset M_n(\mathbb{Z})$ 和一个有理数 $r$，判断 $\rho(\Sigma) \leq r$ 还是 $\rho(\Sigma) > r$。

    更具体地，判断 $\rho(\Sigma) \leq 1$（即切换系统的稳定性）对于一般整数矩阵集是不可判定的。

!!! note "注"
    不可判定性是通过将有界性问题归约到 Post 对应问题（PCP）或 Turing 机的停机问题来证明的。这意味着不存在通用算法在有限步内精确计算联合谱半径。

### 近似算法

尽管精确计算不可判定，但联合谱半径可以被**近似**到任意精度。

!!! definition "定义 63.6 (多面体范数方法)"
    **多面体范数方法**通过寻找一个多面体范数 $\|\cdot\|_P$（其单位球是多面体 $P$）来逼近联合谱半径。若存在 $\gamma > 0$ 使得对所有 $A \in \Sigma$，$\|Ax\|_P \leq \gamma \|x\|_P$，则 $\rho(\Sigma) \leq \gamma$。

    核心思想是寻找一个**不变多面体**（invariant polytope）$P$ 使得 $A_i P \subseteq \gamma P$ 对所有 $i$ 成立。

!!! theorem "定理 63.10 (SOS 方法)"
    **平方和**（Sum of Squares, SOS）方法利用半定规划（SDP）来构造 Lyapunov 函数。对于给定的次数 $2d$，寻找正定齐次多项式 $V(x) = x^T Q x$（或更高次的）使得
    $$
    V(A_i x) \leq \gamma^2 V(x), \quad \forall x, \; \forall i
    $$
    这可以写成 LMI：$A_i^T Q A_i \preceq \gamma^2 Q$，$Q \succ 0$。

    若此 LMI 可行，则 $\rho(\Sigma) \leq \gamma$。

!!! example "例 63.6"
    对于 $\Sigma = \{A_1, A_2\}$，其中 $A_1 = \begin{pmatrix} 0.5 & 0.6 \\ 0 & 0.5 \end{pmatrix}$，$A_2 = \begin{pmatrix} 0.5 & 0 \\ 0.6 & 0.5 \end{pmatrix}$。

    使用 SOS 方法（$d=1$，即二次 Lyapunov 函数），求解 SDP 可得 $Q = I$ 验证 $A_i^T A_i \preceq \gamma^2 I$，给出上界 $\gamma = \|A_i\|_2 \leq 0.5 + 0.6 = 1.1$（粗略估计）。更精确的 SDP 计算可给出更紧的上界。

    下界可通过计算有限长度乘积的谱半径得到：
    $\rho(A_1)^{1/1} = 0.5$，$\rho(A_1 A_2)^{1/2} \approx 0.781$，$\rho(A_1 A_2 A_1)^{1/3} \approx 0.773$。

### 特殊情形的精确算法

!!! theorem "定理 63.11 (非负矩阵对的 JSR)"
    对于一对非负矩阵 $\Sigma = \{A_1, A_2\}$，$A_1, A_2 \geq 0$，联合谱半径可以在有限步内精确计算。此外，有限性猜想对非负矩阵集成立。

---

## 63.6 切换线性系统的稳定性

<div class="context-flow" markdown>

**核心问题**：切换线性系统 $x_{k+1} = A_{\sigma(k)} x_k$ 在所有切换策略下是否稳定？

</div>

### 切换系统模型

!!! definition "定义 63.7 (离散时间切换线性系统)"
    给定矩阵集 $\Sigma = \{A_1, \ldots, A_m\} \subset M_n(\mathbb{R})$ 和切换信号 $\sigma: \mathbb{N} \to \{1, \ldots, m\}$，**离散时间切换线性系统**定义为
    $$
    x_{k+1} = A_{\sigma(k)} x_k, \quad k = 0, 1, 2, \ldots
    $$
    其中 $x_k \in \mathbb{R}^n$ 是状态向量，$\sigma(k)$ 在每个时刻 $k$ 选择一个矩阵。

!!! definition "定义 63.8 (绝对渐近稳定性)"
    切换系统称为**绝对渐近稳定**（absolutely asymptotically stable），若对**所有**初始条件 $x_0$ 和**所有**切换信号 $\sigma$，状态 $x_k \to 0$（$k \to \infty$）。

### 稳定性判据

!!! theorem "定理 63.12 (切换系统稳定性与 JSR)"
    离散时间切换线性系统 $x_{k+1} = A_{\sigma(k)} x_k$ 绝对渐近稳定，当且仅当
    $$
    \rho(\Sigma) < 1
    $$

??? proof "证明"
    **充分性**：若 $\rho(\Sigma) < 1$，则存在 $\gamma < 1$ 和矩阵范数 $\|\cdot\|$（可取极端范数）使得 $\|A\| \leq \gamma$ 对所有 $A \in \Sigma$ 成立。从而
    $$
    \|x_k\| = \|A_{\sigma(k-1)} \cdots A_{\sigma(0)} x_0\| \leq \gamma^k \|x_0\| \to 0
    $$

    **必要性**：若 $\rho(\Sigma) \geq 1$，由定义存在序列 $(i_1^{(k)}, \ldots, i_k^{(k)})$ 使得 $\|A_{i_1^{(k)}} \cdots A_{i_k^{(k)}}\|^{1/k} \to \rho(\Sigma) \geq 1$。选取使乘积范数最大的切换序列，可以构造不收敛到零的轨道。

### 公共 Lyapunov 函数

!!! definition "定义 63.9 (公共二次 Lyapunov 函数)"
    矩阵集 $\Sigma = \{A_1, \ldots, A_m\}$ 具有**公共二次 Lyapunov 函数**（Common Quadratic Lyapunov Function, CQLF），若存在正定矩阵 $P \succ 0$ 使得
    $$
    A_i^T P A_i \prec P, \quad \forall i = 1, \ldots, m
    $$

    此时 $V(x) = x^T P x$ 是公共 Lyapunov 函数。

!!! theorem "定理 63.13 (CQLF 与稳定性)"
    若 $\Sigma$ 具有 CQLF，则 $\rho(\Sigma) < 1$。反之不然：存在 $\rho(\Sigma) < 1$ 但不存在 CQLF 的矩阵集。

!!! example "例 63.7"
    设 $A_1 = \begin{pmatrix} 0.5 & 0.3 \\ 0 & 0.4 \end{pmatrix}$，$A_2 = \begin{pmatrix} 0.4 & 0 \\ 0.3 & 0.5 \end{pmatrix}$。

    由于 $\rho(A_1) = 0.5 < 1$，$\rho(A_2) = 0.5 < 1$，但单独的谱半径不足以保证切换系统稳定。

    取 $P = I$，检验：
    $$
    A_1^T A_1 = \begin{pmatrix} 0.25 & 0.15 \\ 0.15 & 0.25 \end{pmatrix} \prec I, \quad
    A_2^T A_2 = \begin{pmatrix} 0.25 & 0.15 \\ 0.15 & 0.34 \end{pmatrix} \prec I
    $$
    因此 $P = I$ 是 CQLF，系统在所有切换下稳定。

!!! example "例 63.8"
    设 $A_1 = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$A_2 = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}$。

    $\rho(A_1) = \rho(A_2) = 0$（均为幂零矩阵），但 $A_1 A_2 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$\rho(A_1 A_2) = 1$。

    切换序列 $\sigma = (1, 2, 1, 2, \ldots)$ 给出 $x_{2k} = (A_1 A_2)^k x_0$，不收敛到零。因此虽然每个矩阵的谱半径为零，$\rho(\Sigma) = 1$，系统不稳定。

---

## 63.7 小波中的联合谱半径

<div class="context-flow" markdown>

**核心问题**：如何用联合谱半径来判断小波细分格式的收敛性和尺度函数的光滑性？

</div>

### 细分格式

!!! definition "定义 63.10 (二进细分格式)"
    给定一组**细分系数** $\{c_0, c_1, \ldots, c_N\} \subset \mathbb{R}$，二进细分格式（binary subdivision scheme）定义迭代过程：给定控制点序列 $\{f_i^{(j)}\}$，下一层的控制点由
    $$
    f_i^{(j+1)} = \sum_{k=0}^{N} c_k \, f_{i-k}^{(j)}
    $$
    定义，其中偶数和奇数位置分别用不同的规则。

    等价地，定义两个**转移矩阵**（transition matrices）$T_0$ 和 $T_1$：
    $$
    (T_\epsilon)_{ij} = c_{2i - j + \epsilon}, \quad \epsilon = 0, 1
    $$

### 收敛性判据

!!! theorem "定理 63.14 (细分格式的收敛性)"
    二进细分格式收敛（即存在连续的极限函数），当且仅当
    $$
    \rho(\{T_0|_V, T_1|_V\}) < 1
    $$
    其中 $V$ 是差分空间（由向量 $e_i - e_{i+1}$ 张成的子空间），$T_0|_V, T_1|_V$ 是转移矩阵在 $V$ 上的限制。

!!! example "例 63.9"
    **Daubechies-Lagarias 定理**：考虑系数 $c = (1/2, 1, 1/2)$（线性 B 样条的细分系数）。转移矩阵在差分空间上的限制为
    $$
    T_0|_V = \begin{pmatrix} 1/2 \end{pmatrix}, \quad T_1|_V = \begin{pmatrix} 1/2 \end{pmatrix}
    $$
    $\rho(\{T_0|_V, T_1|_V\}) = 1/2 < 1$，因此格式收敛。极限函数是线性 B 样条（帽子函数）。

### 光滑性判据

!!! theorem "定理 63.15 (尺度函数的 Sobolev 正则性)"
    设细分格式的转移矩阵在 $s$ 阶差分空间上的限制为 $\{T_0^{(s)}, T_1^{(s)}\}$。尺度函数 $\phi$ 属于 Sobolev 空间 $W^{s,p}$，当且仅当
    $$
    \rho_p(\{T_0^{(s)}, T_1^{(s)}\}) < 2^{s - 1 + 1/p}
    $$
    其中 $\rho_p$ 是 $L^p$ 联合谱半径。

    特别地，$\phi$ 的 Holder 正则性指数为
    $$
    \alpha = s - \log_2 \rho(\{T_0^{(s)}, T_1^{(s)}\})
    $$
    其中 $s$ 是格式的逼近阶数。

!!! example "例 63.10"
    对于 Daubechies 小波 $D_4$（4 系数），细分系数为
    $$
    c_0 = \frac{1+\sqrt{3}}{4\sqrt{2}}, \; c_1 = \frac{3+\sqrt{3}}{4\sqrt{2}}, \; c_2 = \frac{3-\sqrt{3}}{4\sqrt{2}}, \; c_3 = \frac{1-\sqrt{3}}{4\sqrt{2}}
    $$

    计算限制转移矩阵在一阶差分空间上的联合谱半径，可以确定尺度函数 $\phi$ 的 Holder 正则性约为 $\alpha \approx 1.09$，即 $\phi$ 略优于 Lipschitz 连续。

### 多尺度分析中的 JSR

!!! note "注"
    联合谱半径在小波理论中的应用由 Daubechies 和 Lagarias（1991, 1992）系统建立。他们证明了**二分量细分**的收敛性和光滑性完全由转移矩阵的 JSR 控制。这一工作开创了联合谱半径理论的应用研究。

    后续工作将这一框架推广到：

    - **多元小波**（multivariate wavelets）：膨胀矩阵确定多个转移矩阵；
    - **向量值细分**（vector subdivision）：每个控制点是向量，转移矩阵更大；
    - **非平稳细分**（non-stationary subdivision）：每层使用不同的细分规则，需要计算无限矩阵集的 JSR。

---

## 本章小结

本章研究了矩阵集合的谱理论。主要结果包括：

1. **同时三角化**：McCoy 定理给出了代数闭域上同时三角化的充要条件；可交换矩阵族和可解 Lie 代数都可同时三角化（Lie 定理）。

2. **同时对角化**：可交换且各自可对角化的矩阵族可同时对角化。对正规矩阵可以做到同时酉对角化。

3. **联合谱半径** $\rho(\Sigma)$ 量化了矩阵集合所有可能乘积的最大增长率。它由 Rota-Strang 引入，是谱半径概念的自然推广。

4. **Berger-Wang 定理**建立了联合谱半径与广义谱半径的相等，将"范数量"与"特征值量"联系起来。

5. **计算复杂性**：一般矩阵集的 JSR 是不可判定的，但可以通过多面体范数和 SOS/SDP 方法任意精度地逼近。

6. **切换系统**：切换线性系统的绝对稳定性等价于 $\rho(\Sigma) < 1$。

7. **小波理论**：细分格式的收敛性和尺度函数的光滑性由限制转移矩阵的 JSR 控制。

---

## 习题

!!! question "习题 63.1"
    证明两个上三角矩阵 $A, B \in M_n(\mathbb{C})$ 的换位子 $[A, B] = AB - BA$ 是严格上三角的。

!!! question "习题 63.2"
    设 $A, B \in M_2(\mathbb{C})$ 且 $AB = BA$。直接证明（不使用 McCoy 定理）$A, B$ 可以同时三角化。

!!! question "习题 63.3"
    设 $A_1, \ldots, A_m$ 都是实对称矩阵。证明它们可同时（正交）对角化当且仅当两两可交换。

!!! question "习题 63.4"
    计算 $\Sigma = \left\{\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 1 & 1 \end{pmatrix}\right\}$ 的联合谱半径。

!!! question "习题 63.5"
    设 $\Sigma = \{\alpha A : \alpha \in [0, 1]\}$，其中 $A$ 为固定矩阵。证明 $\rho(\Sigma) = \rho(A)$。

!!! question "习题 63.6"
    设 $A_1, A_2$ 是两个 $2 \times 2$ 幂零矩阵。证明 $\rho(\{A_1, A_2\}) = \sqrt{\rho(A_1 A_2)}$。

!!! question "习题 63.7"
    证明：若 $\Sigma$ 中所有矩阵两两可交换，则 $\rho(\Sigma) = \max_{A \in \Sigma} \rho(A)$。

!!! question "习题 63.8"
    给出一个切换系统 $x_{k+1} = A_{\sigma(k)} x_k$ 的例子，其中每个 $A_i$ 都是 Schur 稳定的（$\rho(A_i) < 1$），但切换系统不稳定（$\rho(\Sigma) \geq 1$）。

!!! question "习题 63.9"
    证明 Berger-Wang 定理中不等式 $\hat{\rho}(\Sigma) \leq \rho(\Sigma)$ 的方向。

!!! question "习题 63.10"
    对于系数 $c = (1/4, 1/2, 1/4)$ 的细分格式，计算转移矩阵并判断收敛性。
