# 第 39 章 全正矩阵

<div class="context-flow" markdown>

**前置**：行列式(Ch3) · 特征值(Ch6) · 非负矩阵(Ch17)

**本章脉络**：全正(TP)与全非负(TN)定义 $\to$ 基本性质 $\to$ 振荡矩阵 $\to$ 特征值性质 $\to$ 双对角分解 $\to$ 平面网络 $\to$ 与 cluster algebra 的联系

**延伸**：全正性在逼近论（B-样条的全正性保证形状保持）、组合学（cluster algebra 与正 Grassmannian）、概率论（Polya 频率序列）中有深刻应用

</div>

全正矩阵是矩阵论中一个优美而深刻的专题。一个矩阵如果它的**所有**子式（不仅仅是主子式）都为正，就称为全正的。这个看似极端的条件蕴含着令人惊讶的丰富结构：特征值全为正且互不相同，特征向量具有严格的振荡性质，矩阵可以分解为非负初等双对角矩阵的乘积。全正矩阵理论的根源可以追溯到 Schoenberg 和 Gantmacher-Krein 在 20 世纪 30-50 年代的工作，而近年来通过 Fomin 和 Zelevinsky 的 cluster algebra 理论又获得了全新的活力。

---

## 39.1 全正与全非负矩阵的定义

<div class="context-flow" markdown>

**核心问题**：什么样的矩阵的所有子式都为正（或非负）？这类矩阵是否常见？

</div>

!!! definition "定义 39.1 (全正矩阵, TP)"
    矩阵 $A \in M_{m \times n}(\mathbb{R})$ 称为**全正的**（Totally Positive, TP），若 $A$ 的所有子式（minor）为正：对所有 $1 \le i_1 < \cdots < i_k \le m$，$1 \le j_1 < \cdots < j_k \le n$，$1 \le k \le \min(m, n)$，
    $$\det(A[\{i_1, \ldots, i_k\}, \{j_1, \ldots, j_k\}]) > 0.$$

!!! definition "定义 39.2 (全非负矩阵, TN)"
    矩阵 $A \in M_{m \times n}(\mathbb{R})$ 称为**全非负的**（Totally Nonnegative, TN），若 $A$ 的所有子式非负：
    $$\det(A[\{i_1, \ldots, i_k\}, \{j_1, \ldots, j_k\}]) \ge 0.$$

!!! note "注记 39.1 (术语说明)"
    文献中有时使用不同的术语：Karlin 使用"totally positive"表示全非负（所有子式 $\ge 0$），"strictly totally positive"表示全正（所有子式 $> 0$）。本书采用 Ando, Pinkus 等人的较现代约定：TP = 全正（严格），TN = 全非负。读者在查阅文献时需注意术语差异。

!!! example "例 39.1"
    (a) 任何正实数 $(a)$（$1 \times 1$ 矩阵）是 TP 的。

    (b) $A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$：子式为 $a_{11} = 1 > 0$, $a_{12} = 1 > 0$, $a_{21} = 1 > 0$, $a_{22} = 2 > 0$, $\det(A) = 1 > 0$。$A$ 是 TP 的。

    (c) $A = \begin{pmatrix} 1 & 2 & 3 \\ 2 & 5 & 8 \\ 3 & 8 & 14 \end{pmatrix}$：需要检验所有 $\binom{3}{1}^2 + \binom{3}{2}^2 + 1 = 9 + 9 + 1 = 19$ 个子式。所有元素为正。$2 \times 2$ 子式包括 $\det\begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix} = 1$, $\det\begin{pmatrix} 1 & 3 \\ 2 & 8 \end{pmatrix} = 2$, $\det\begin{pmatrix} 2 & 3 \\ 5 & 8 \end{pmatrix} = 1$, $\det\begin{pmatrix} 1 & 2 \\ 3 & 8 \end{pmatrix} = 2$, $\det\begin{pmatrix} 1 & 3 \\ 3 & 14 \end{pmatrix} = 5$, $\det\begin{pmatrix} 2 & 3 \\ 8 & 14 \end{pmatrix} = 4$, 以及类似地 $\det\begin{pmatrix} 2 & 5 \\ 3 & 8 \end{pmatrix} = 1$, $\det\begin{pmatrix} 2 & 8 \\ 3 & 14 \end{pmatrix} = 4$, $\det\begin{pmatrix} 5 & 8 \\ 8 & 14 \end{pmatrix} = 6$。$\det(A) = 1$。全部为正。$A$ 是 TP 的。

    (d) 指数 Vandermonde 矩阵 $A = (e^{x_i t_j})$（$x_1 < \cdots < x_m$，$t_1 < \cdots < t_n$）是 TP 的。这是全正矩阵理论的一个经典事实。

!!! example "例 39.2 (Pascal 矩阵)"
    Pascal 矩阵 $P = \left(\binom{i+j-2}{i-1}\right)_{i,j=1}^n$ 是全非负的（实际上是全正的）。例如
    $$P_4 = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & 2 & 3 & 4 \\ 1 & 3 & 6 & 10 \\ 1 & 4 & 10 & 20 \end{pmatrix}.$$
    这可以通过组合恒等式 $\binom{m}{k} = \sum \binom{m_1}{k_1}\binom{m_2}{k_2}$ 和 Lindstrom-Gessel-Viennot 引理来证明。

---

## 39.2 基本性质

<div class="context-flow" markdown>

**核心问题**：全正性在矩阵运算下如何保持？

</div>

!!! theorem "定理 39.1 (TP 矩阵的乘积)"
    若 $A \in M_{m \times p}(\mathbb{R})$ 和 $B \in M_{p \times n}(\mathbb{R})$ 都是 TP 的，则 $AB$ 是 TP 的。

??? proof "证明"
    这是 **Cauchy-Binet 公式**的直接推论。对 $AB$ 的任意 $k$ 阶子式：
    $$\det((AB)[\alpha, \beta]) = \sum_{\gamma} \det(A[\alpha, \gamma]) \det(B[\gamma, \beta]),$$
    其中求和遍历所有 $k$ 元子集 $\gamma \subseteq \{1, \ldots, p\}$。

    由 $A$ 和 $B$ 的 TP 性，每个求和项 $\det(A[\alpha, \gamma]) \det(B[\gamma, \beta]) > 0$，故总和为正。

!!! theorem "定理 39.2 (TP/TN 矩阵的基本性质)"
    (a) 若 $A$ 是 TP (TN)，则 $A^T$ 是 TP (TN)。

    (b) 若 $A$ 是 $n \times n$ TP 矩阵，则 $A$ 非奇异且 $A^{-1}$ 具有棋盘符号模式：$(A^{-1})_{ij}$ 的符号为 $(-1)^{i+j}$（对角及相邻位置为正或负交替）。

    (c) TP 矩阵的行列**不可以**自由排列：行或列的置换一般会破坏全正性。

    (d) 若 $D_1, D_2$ 为正对角矩阵，$A$ 为 TP，则 $D_1 A D_2$ 为 TP。

??? proof "证明"
    (a) $\det(A^T[\alpha, \beta]) = \det(A[\beta, \alpha]) > 0$。

    (b) 由 Cramer 法则，$(A^{-1})_{ij} = (-1)^{i+j} \det(A[\hat{j}, \hat{i}]) / \det(A)$，其中 $A[\hat{j}, \hat{i}]$ 是去掉第 $j$ 行第 $i$ 列的子矩阵。$A$ 的所有子式为正，故 $\det(A[\hat{j}, \hat{i}]) > 0$，$(A^{-1})_{ij}$ 的符号恰为 $(-1)^{i+j}$。

    (d) 对角缩放不改变子式的符号：$\det((D_1AD_2)[\alpha,\beta]) = \prod_{i \in \alpha}(d_1)_i \cdot \det(A[\alpha,\beta]) \cdot \prod_{j \in \beta}(d_2)_j > 0$。

!!! theorem "定理 39.3 (TN 矩阵的子矩阵)"
    TN 矩阵的任何子矩阵（不一定是主子矩阵）仍然是 TN 的。TP 矩阵的连续行列子矩阵是 TP 的。

??? proof "证明"
    子矩阵的子式是原矩阵的子式，继承非负性。对于 TP，连续行列子矩阵的子式仍是原矩阵的子式，继承正性。对非连续子矩阵，TP 性不一定保持（因为原矩阵的某些子式可能不出现在子矩阵中）。实际上 TP 矩阵的任何子矩阵也是 TP 的，这是因为 TP 矩阵的任意 $k$ 行 $k$ 列的子式都为正。

---

## 39.3 振荡矩阵

<div class="context-flow" markdown>

**核心问题**：什么条件下全非负矩阵的某个幂次变为全正的？

</div>

!!! definition "定义 39.3 (振荡矩阵)"
    方阵 $A \in M_n(\mathbb{R})$ 称为**振荡矩阵**（oscillatory matrix），若 $A$ 是全非负的且存在正整数 $p$ 使得 $A^p$ 是全正的。

!!! theorem "定理 39.4 (Gantmacher-Krein 振荡矩阵判据)"
    全非负矩阵 $A \in M_n(\mathbb{R})$ 是振荡矩阵当且仅当以下两个条件同时成立：

    (a) $A$ 非奇异（$\det(A) > 0$）；

    (b) $A$ 是**全不可约的**：对 $i = 1, \ldots, n-1$，有 $a_{i,i+1} > 0$ 且 $a_{i+1,i} > 0$。

    当这两个条件成立时，$A^{n-1}$ 就已经是全正的。

??? proof "证明（概要）"
    **必要性**：若 $A^p$ TP，则 $\det(A^p) = (\det A)^p > 0$，故 $\det(A) > 0$。全不可约性来自对相邻行列的 $2 \times 2$ 子式的分析。

    **充分性**（核心思路）：条件 (b) 保证 $A$ 的图（将 $A$ 视为有向图的邻接矩阵）是强连通的。类似于 Perron-Frobenius 理论中不可约非负矩阵的论证，全不可约性加上全非负性可以通过反复相乘"传播正性"：$A^k$ 的越来越多的子式变为严格正的，最终在 $k = n - 1$ 步后所有子式都为正。

    严格的证明需要利用 Cauchy-Binet 公式和一个关于全非负矩阵中零子式的结构定理（"消灭"零子式的归纳论证）。

!!! example "例 39.3"
    三对角矩阵 $A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 2 & 1 \\ 0 & 1 & 2 \end{pmatrix}$。

    检验 TN 性：所有元素 $\ge 0$。$2 \times 2$ 子式：$\det\begin{pmatrix} 2 & 1 \\ 1 & 2\end{pmatrix} = 3$，$\det\begin{pmatrix} 2 & 0 \\ 1 & 1\end{pmatrix} = 2$，$\det\begin{pmatrix} 1 & 0 \\ 2 & 1\end{pmatrix} = 1$，$\det\begin{pmatrix} 1 & 1 \\ 1 & 2\end{pmatrix} = 1$，以及对称情形。$\det(A) = 4$。所有子式 $\ge 0$，且 $\det(A) > 0$。

    全不可约性：$a_{12} = 1 > 0$, $a_{21} = 1 > 0$, $a_{23} = 1 > 0$, $a_{32} = 1 > 0$。

    由 Gantmacher-Krein 定理，$A$ 是振荡矩阵，$A^2$ 应为 TP。

    $A^2 = \begin{pmatrix} 5 & 4 & 1 \\ 4 & 6 & 4 \\ 1 & 4 & 5 \end{pmatrix}$。所有元素为正。$2 \times 2$ 子式：如 $\det\begin{pmatrix}5&4\\4&6\end{pmatrix} = 14 > 0$，$\det\begin{pmatrix}5&1\\4&4\end{pmatrix} = 16 > 0$，等等。$\det(A^2) = 16 > 0$。$A^2$ 确实是 TP 的。

!!! note "注记 39.1b (振荡矩阵的物理意义)"
    振荡矩阵的名称来源于力学中的**小振动理论**。Gantmacher 和 Krein 研究了弹性体（如弦、杆、梁）的离散化模型，发现刚度矩阵和质量矩阵的组合往往是振荡矩阵。振荡矩阵的特征值性质（正、简单、特征向量有确定的振荡模式）恰好对应于物理系统的振动模态：第 $k$ 个模态恰有 $k-1$ 个节点（零点），这与实验观测完全一致。

    例如，弦振动的离散化产生三对角振荡矩阵，其特征向量是离散正弦函数，恰好有定理 39.6 所预测的符号变化数。

!!! example "例 39.3b"
    考虑 $4 \times 4$ Jacobi 矩阵（对称三对角矩阵，次对角元素为正）：
    $$A = \begin{pmatrix} 2 & 1 & 0 & 0 \\ 1 & 2 & 1 & 0 \\ 0 & 1 & 2 & 1 \\ 0 & 0 & 1 & 2 \end{pmatrix}.$$
    这是振荡矩阵（TN + 非奇异 + 全不可约）。由定理 39.4，$A^3$ 为 TP。其特征值 $\lambda_k = 2 + 2\cos(k\pi/5)$（$k = 1,2,3,4$）全部为正且互不相同。

---

## 39.4 特征值与特征向量性质

<div class="context-flow" markdown>

**核心问题**：全正矩阵的特征值和特征向量有什么特殊的定性性质？

</div>

这是 Gantmacher 和 Krein 的经典理论的核心，揭示了全正矩阵具有"最佳可能"的谱性质。

!!! theorem "定理 39.5 (TP 矩阵的特征值)"
    设 $A \in M_n(\mathbb{R})$ 为全正矩阵。则 $A$ 的特征值 $\lambda_1 > \lambda_2 > \cdots > \lambda_n > 0$ 全为正实数且互不相同。

??? proof "证明（概要）"
    **步骤一（特征值为正实数）**：$A$ 是 TP 的，故所有主子式为正，特别是 $A$ 的迹 $> 0$、行列式 $> 0$。但这不足以证明所有特征值为正。

    关键工具是**复合矩阵**（compound matrix）。$A$ 的 $k$ 阶复合矩阵 $C_k(A)$ 的元素是 $A$ 的所有 $k$ 阶子式。若 $A$ 是 TP 的，则 $C_k(A) > 0$（逐元为正）。由 Perron-Frobenius 定理，$C_k(A)$ 有正的 Perron 特征值，等于 $A$ 的前 $k$ 个最大特征值之积 $\lambda_1 \lambda_2 \cdots \lambda_k$。由此可递推得出每个 $\lambda_k > 0$。

    **步骤二（特征值互不相同）**：若 $\lambda_k = \lambda_{k+1}$，则 $C_k(A)$ 的 Perron 特征值 $\lambda_1 \cdots \lambda_k$ 将是 $C_{k+1}(A)$ 对应特征值的因子，但详细的 Perron-Frobenius 分析表明这要求 $C_k(A)$ 具有重复的主特征值，与 $C_k(A) > 0$（不可约且基元）矛盾。

!!! theorem "定理 39.6 (特征向量的振荡性质)"
    设 $A \in M_n(\mathbb{R})$ 为 TP 矩阵，特征值 $\lambda_1 > \lambda_2 > \cdots > \lambda_n > 0$，对应特征向量 $v^{(1)}, v^{(2)}, \ldots, v^{(n)}$。则 $v^{(k)}$ 恰有 $k - 1$ 个**符号变化**（sign changes）。

    更精确地，设 $v = (v_1, \ldots, v_n)^T \ne 0$，定义 $v$ 的**符号变化数** $S^-(v)$ 为将零分量去除后，相邻非零分量符号不同的次数。则
    $$S^-(v^{(k)}) = k - 1, \quad k = 1, 2, \ldots, n.$$

??? proof "证明（概要）"
    对 $k = 1$：$v^{(1)}$ 是 Perron 向量，由于 $A > 0$（TP 蕴含所有元素为正），Perron-Frobenius 定理保证 $v^{(1)} > 0$（所有分量正），故 $S^-(v^{(1)}) = 0$。

    对一般 $k$，证明利用 TP 矩阵的**变差减少性**（variation diminishing property）：若 $x$ 有 $s$ 个符号变化，则 $Ax$ 至多有 $s$ 个符号变化。结合特征值的交错性质和线性无关性，可以证明 $v^{(k)}$ 恰有 $k-1$ 个符号变化。

!!! definition "定义 39.4 (变差减少性)"
    矩阵 $A \in M_{m \times n}(\mathbb{R})$ 称为**变差减少的**（variation diminishing），若对所有 $x \in \mathbb{R}^n$，
    $$S^-(Ax) \le S^-(x),$$
    其中 $S^-(v)$ 计算向量 $v$ 的符号变化数（零分量用使符号变化数最大的方式处理）。

!!! theorem "定理 39.7 (Schoenberg 变差减少定理)"
    矩阵 $A$ 是全非负的当且仅当 $A$ 是变差减少的。

??? proof "证明（必要性方向概要）"
    设 $A$ TN，$x \in \mathbb{R}^n$，$y = Ax$。需证 $S^-(y) \le S^-(x)$。

    不妨设 $x$ 有 $s$ 个符号变化。将 $\{1, \ldots, n\}$ 分为 $s+1$ 个连续块 $B_1, \ldots, B_{s+1}$，使得 $x$ 在每个块内同号。

    若 $y$ 有 $\ge s + 1$ 个符号变化，则存在 $s + 2$ 个指标 $i_1 < \cdots < i_{s+2}$ 使得 $y$ 在这些位置交替变号。构造特定的线性组合并利用 Cauchy-Binet 和 TN 性质导出矛盾。

!!! example "例 39.4"
    验证定理 39.5 对 $A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$（TP）的适用性。

    特征多项式：$\lambda^2 - 3\lambda + 1 = 0$，$\lambda_{1,2} = \frac{3 \pm \sqrt{5}}{2}$。$\lambda_1 \approx 2.618 > \lambda_2 \approx 0.382 > 0$。两个特征值均为正且不同，验证了定理 39.5。

    $v^{(1)} = (1, \frac{1+\sqrt{5}}{2})^T$（分量同号，$S^- = 0$），$v^{(2)} = (1, \frac{1-\sqrt{5}}{2})^T$（分量异号，$S^- = 1$）。验证了定理 39.6。

---

## 39.5 双对角分解

<div class="context-flow" markdown>

**核心问题**：全非负矩阵能否分解为更简单的全非负矩阵的乘积？

</div>

!!! definition "定义 39.5 (非负初等双对角矩阵)"
    **非负初等下双对角矩阵**是形如
    $$E_k(m) = I + m \cdot e_{k+1} e_k^T, \quad m \ge 0$$
    的矩阵（在 $(k+1, k)$ 位置有一个非负元素 $m$，其余为单位矩阵）。类似地定义上双对角矩阵 $\hat{E}_k(m) = I + m \cdot e_k e_{k+1}^T$。

!!! theorem "定理 39.8 (Loewner-Whitney 定理; 双对角分解)"
    矩阵 $A \in M_n(\mathbb{R})$ 是全非负的（TN）当且仅当 $A$ 可以写成
    $$A = \prod_{\text{有序}} E_{k_i}(m_i) \cdot D \cdot \prod_{\text{有序}} \hat{E}_{l_j}(\hat{m}_j),$$
    即非负初等下双对角矩阵的乘积、正对角矩阵、非负初等上双对角矩阵的乘积。

    更精确地，$A$ 可以分解为
    $$A = L_1 L_2 \cdots L_{n-1} \cdot D \cdot U_{n-1} \cdots U_2 U_1,$$
    其中每个 $L_k$ 是非负下双对角矩阵（仅第 $k$ 条次对角线非零），$U_k$ 是非负上双对角矩阵，$D$ 为非负对角矩阵。

??? proof "证明（概要）"
    **必要性**（TN $\Rightarrow$ 有此分解）：这可以通过对 $A$ 进行 **Neville 消元**（也称全正性保持的 Gauss 消元）来证明。与通常的 Gauss 消元不同，Neville 消元使用相邻行的操作：$\text{row}_i \leftarrow \text{row}_i - m \cdot \text{row}_{i-1}$（$m \ge 0$），这等价于左乘 $E_{i-1}(-m)^{-1} = E_{i-1}(m)^{-1}$。TN 性保证所有乘子 $m \ge 0$。经过 $\binom{n}{2}$ 步消元，$A$ 被化为上三角形式，反转得到所需分解。

    **充分性**：每个 $E_k(m)$（$m \ge 0$）都是 TN 的（因为它是单位矩阵加一个非负秩 1 矩阵，所有子式 $\ge 0$）。TN 矩阵的乘积是 TN 的（定理 39.1 的 TN 版本）。

!!! example "例 39.5"
    对 TN 矩阵 $A = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}$ 进行双对角分解。

    $\det(A) = 7 - 6 = 1 > 0$。$A$ 实际上是 TP 的。

    Neville 消元：$\text{row}_2 \leftarrow \text{row}_2 - 3 \cdot \text{row}_1$：得 $\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$。

    即 $E_1(-3) A = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$，$A = E_1(3) \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = E_1(3) \cdot D \cdot \hat{E}_1(2)$，

    其中 $D = I$，$E_1(3) = \begin{pmatrix} 1 & 0 \\ 3 & 1 \end{pmatrix}$，$\hat{E}_1(2) = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$。

    验证：$E_1(3) \cdot I \cdot \hat{E}_1(2) = \begin{pmatrix} 1 & 0 \\ 3 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix} = A$。正确。

!!! note "注记 39.2 (平面网络解释)"
    双对角分解有一个优美的**平面网络**（planar network）解释。每个初等双对角矩阵对应网络中的一条"斜边"，权重为 $m$。矩阵 $A$ 的 $(i,j)$ 元素等于从源点 $j$ 到汇点 $i$ 的所有路径的权重之和（权重为路径上所有边权的乘积）。TN 条件等价于所有路径权重非负（因为 $m \ge 0$），而任意子式等于 Lindstrom-Gessel-Viennot (LGV) 引理给出的不相交路径族的权重之和。

!!! theorem "定理 39.8b (双对角分解的参数化)"
    对 $n \times n$ 非奇异 TN 矩阵，双对角分解中的参数个数为 $\binom{n}{2}$（下三角部分）+ $n$（对角线）+ $\binom{n}{2}$（上三角部分）= $n^2$ 个非负参数。这 $n^2$ 个参数提供了 TN 矩阵空间的一个**全局参数化**——它将非奇异 TN 矩阵的集合与 $\mathbb{R}_{>0}^n \times \mathbb{R}_{\ge 0}^{n(n-1)}$ 建立了一一对应（对角参数为正，其余为非负）。

    对 TP 矩阵，所有 $n^2$ 个参数均严格为正：$m_i > 0, d_j > 0, \hat{m}_k > 0$。

!!! example "例 39.5b"
    对 $3 \times 3$ TN 矩阵，双对角分解为
    $$A = E_1(m_1) E_2(m_2) E_1(m_3) \cdot D \cdot \hat{E}_1(\hat{m}_3) \hat{E}_2(\hat{m}_2) \hat{E}_1(\hat{m}_1),$$
    其中 $D = \operatorname{diag}(d_1, d_2, d_3)$。9 个参数：$m_1, m_2, m_3 \ge 0$（下三角），$d_1, d_2, d_3 > 0$（对角），$\hat{m}_1, \hat{m}_2, \hat{m}_3 \ge 0$（上三角）。

    展开：
    $$E_1(m_1) = \begin{pmatrix} 1 & 0 & 0 \\ m_1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad E_2(m_2) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & m_2 & 1 \end{pmatrix}, \quad E_1(m_3) = \begin{pmatrix} 1 & 0 & 0 \\ m_3 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}.$$
    类似地定义上双对角因子。

---

## 39.6 判定算法

<div class="context-flow" markdown>

**核心问题**：$n \times n$ 矩阵有 $\sum_{k=1}^n \binom{n}{k}^2$ 个子式，如何高效判定全正性？

</div>

一个 $n \times n$ 矩阵共有 $\sum_{k=1}^n \binom{n}{k}^2 = \binom{2n}{n} - 1$ 个子式，数量随 $n$ 指数增长。幸运的是，不需要检验所有子式。

!!! theorem "定理 39.9 (初始子式判据)"
    矩阵 $A \in M_n(\mathbb{R})$（$A$ 的所有元素为正）是 TP 的当且仅当 $A$ 的所有**初始子式**（initial minors）为正。初始子式是指行指标集和列指标集均为**连续的**（即形如 $\{i, i+1, \ldots, i+k-1\}$）子式。

    初始子式的数量为 $\binom{n+1}{2}^2$ 级别——一个关于 $n$ 的多项式，远少于指数级的总子式数量。更精确地，只需检验 $n^2$ 个**实心子式**（contiguous minors，行列指标均连续）。

!!! definition "定义 39.6 (Neville 消元)"
    **Neville 消元**是一种专门用于全非负矩阵的消元算法，仅使用相邻行的行操作。算法过程如下：

    对 $k = 1, 2, \ldots, n-1$：
    对 $i = n, n-1, \ldots, k+1$：
    $$a_{ij}^{(\text{new})} = a_{ij} - \frac{a_{i,k}}{a_{i-1,k}} a_{i-1,j}, \quad j = k, k+1, \ldots, n.$$

    若 $A$ 为 TN，则所有乘子 $a_{i,k}/a_{i-1,k} \ge 0$，且所有中间量 $a_{ij}^{(\text{new})} \ge 0$。

!!! theorem "定理 39.10 (TN 判定)"
    矩阵 $A$（所有元素 $\ge 0$）是全非负的当且仅当 Neville 消元可以完成（不出现除以零或负数的情况），且所有中间量非负。

!!! example "例 39.6"
    判定 $A = \begin{pmatrix} 1 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 2 \end{pmatrix}$ 是否 TN。

    所有元素 $\ge 0$。进行 Neville 消元：

    第 1 步（$k = 1$）：$\text{row}_3 \leftarrow \text{row}_3 - (1/2)\text{row}_2$：$(1, 3, 2) - (1/2)(2, 5, 3) = (0, 1/2, 1/2)$。
    $\text{row}_2 \leftarrow \text{row}_2 - (2/1)\text{row}_1$：$(2, 5, 3) - 2(1, 2, 1) = (0, 1, 1)$。

    当前矩阵：$\begin{pmatrix} 1 & 2 & 1 \\ 0 & 1 & 1 \\ 0 & 1/2 & 1/2 \end{pmatrix}$。

    第 2 步（$k = 2$）：$\text{row}_3 \leftarrow \text{row}_3 - (1/2)/1 \cdot \text{row}_2$：$(0, 1/2, 1/2) - (1/2)(0, 1, 1) = (0, 0, 0)$。

    所有中间量 $\ge 0$，消元完成。$A$ 是 TN。$\det(A) = 0$，所以 $A$ 是**奇异的** TN 矩阵（不是 TP）。

---

## 39.7 应用

<div class="context-flow" markdown>

**核心问题**：全正矩阵在逼近论、概率论和代数组合学中有哪些重要应用？

</div>

!!! theorem "定理 39.11 (B-样条的全正性)"
    设 $\{N_{i,k}(t)\}$ 为以节点序列 $\{t_i\}$ 定义的 $k$ 阶 B-样条基函数。则**配置矩阵**
    $$A_{ij} = N_{j,k}(\tau_i)$$
    （在适当条件下）是全正的。这一全正性是 B-样条曲线具有**保形性**（shape-preserving property）——即变差减少性——的数学根源。

!!! note "注记 39.3 (B-样条的变差减少性质)"
    设 $c = (c_1, \ldots, c_n)^T$ 为 B-样条控制点的纵坐标，曲线为 $f(t) = \sum_j c_j N_{j,k}(t)$。则 $f$ 的符号变化数不超过 $c$ 的符号变化数。直观地说，B-样条曲线不会比控制多边形振荡更多。这正是变差减少性质（定理 39.7），来源于配置矩阵的全正性。

!!! definition "定义 39.7 (Polya 频率序列)"
    实数序列 $(a_0, a_1, a_2, \ldots)$ 称为 **Polya 频率序列**（PF 序列），若无穷 Toeplitz 矩阵
    $$T = \begin{pmatrix} a_0 & 0 & 0 & \cdots \\ a_1 & a_0 & 0 & \cdots \\ a_2 & a_1 & a_0 & \cdots \\ \vdots & & & \ddots \end{pmatrix}$$
    是全非负的。

!!! theorem "定理 39.12 (Polya 频率序列的生成函数)"
    序列 $(a_k)_{k \ge 0}$ 是 PF 序列当且仅当其生成函数具有形式
    $$\sum_{k=0}^{\infty} a_k z^k = C e^{\gamma z} \prod_{i=1}^{\infty} \frac{1 + \alpha_i z}{1 - \beta_i z},$$
    其中 $C > 0$，$\gamma \ge 0$，$\alpha_i, \beta_i \ge 0$，$\sum(\alpha_i + \beta_i) < \infty$。

!!! example "例 39.7"
    序列 $a_k = 1/k!$ 是 PF 序列，其生成函数为 $e^z = \sum z^k/k!$。这对应于 Polya 频率函数 $f(x) = e^{-x}$（$x \ge 0$），在统计学中与指数分布密切相关。

!!! note "注记 39.4 (Cluster algebra 与全正性)"
    Fomin 和 Zelevinsky 在 2000 年代引入 cluster algebra 的动机之一就是理解"全正性"的代数-组合结构。**正 Grassmannian** $\operatorname{Gr}_{\ge 0}(k, n)$——即所有 Plucker 坐标非负的 Grassmannian 点——的单胞体分解（cell decomposition）与 cluster algebra 的种子（seed）一一对应。这一发现将全正矩阵理论与热带几何、散射振幅（amplituhedron）等前沿课题联系起来。

    在物理学中，Arkani-Hamed 等人发现，$\mathcal{N} = 4$ 超对称 Yang-Mills 理论的散射振幅可以用正 Grassmannian 的几何来计算，而全正性在其中扮演核心角色。

!!! theorem "定理 39.13 (Lindstrom-Gessel-Viennot 引理与全正性)"
    设 $G$ 为带权有向无环图（DAG），$A_1, \ldots, A_n$ 为源点集，$B_1, \ldots, B_n$ 为汇点集。定义**路径矩阵** $M = (m_{ij})$，其中 $m_{ij}$ 等于从 $A_i$ 到 $B_j$ 的所有路径的权重之和。则 $M$ 的 $k$ 阶子式
    $$\det(M[\alpha, \beta])$$
    等于从 $\{A_{\alpha_1}, \ldots, A_{\alpha_k}\}$ 到 $\{B_{\beta_1}, \ldots, B_{\beta_k}\}$ 的**不相交路径族**（non-intersecting path families）的带符号权重之和。

    特别地，若图 $G$ 是**平面的**，则所有不相交路径族的符号均为正，从而 $M$ 的所有子式非负——即 $M$ 为 TN 矩阵。若进一步每对源-汇之间都有正权路径，则 $M$ 为 TP。

??? proof "证明（概要）"
    LGV 引理的证明基于一个优美的对合（involution）论证。考虑从 $\{A_{\alpha_i}\}$ 到 $\{B_{\beta_{\sigma(i)}}\}$ 的所有路径组（对所有置换 $\sigma$）。对任意相交的路径对，构造一个"交换操作"：在第一个交叉点处交换两条路径的尾部。这一操作是一个符号翻转的对合（$\operatorname{sgn}(\sigma)$ 变号），因此相交路径的贡献两两消除，只留下不相交路径族的贡献，且这些贡献的符号恰好对应于 $\operatorname{sgn}(\sigma)$。

    对平面图，不相交路径族只能对应恒等置换（平面性排除了路径的"交叉"），故所有贡献均为正。

!!! example "例 39.8 (Pascal 矩阵的 LGV 解释)"
    考虑整数格点上从 $(0, i-1)$ 到 $(j-1, n-1)$ 的格路径（只允许向右或向上移动），权重均为 1。从 $(0, i-1)$ 到 $(j-1, n-1)$ 的路径数为 $\binom{(j-1)+(n-1)-(i-1)}{j-1} = \binom{n+j-i-1}{j-1}$。

    对 Pascal 矩阵的特殊参数选择，路径矩阵恰好是 Pascal 矩阵。由 LGV 引理和平面性，Pascal 矩阵是全非负的。事实上，由于每对源汇之间都有路径，Pascal 矩阵是全正的。

!!! note "注记 39.5 (全正性的检验复杂度)"
    虽然 $n \times n$ 矩阵的子式总数为 $\binom{2n}{n} - 1 = O(4^n/\sqrt{n})$（指数增长），但全正性的检验可以在**多项式时间**内完成。具体方法有：

    - **Neville 消元**：$O(n^3)$ 次运算。在消元过程中检验所有中间量的符号。
    - **初始子式法**：只需检验 $O(n^2)$ 个连续子式。
    - **双对角分解法**：计算双对角分解并检验所有参数非负，$O(n^2)$。

    这些多项式时间算法的存在性是全正矩阵理论中的一个重要结果，它使得全正性在实际应用中是可计算的。

!!! theorem "定理 39.14 (全正核与全正性的无穷维推广)"
    函数 $K(x, y)$（$x \in X, y \in Y$）称为**全正核**（TP kernel），若对所有 $x_1 < \cdots < x_n \in X$ 和 $y_1 < \cdots < y_n \in Y$，矩阵 $(K(x_i, y_j))$ 都是全正的。

    经典例子包括：

    - **指数核**：$K(x, y) = e^{xy}$（$x, y \in \mathbb{R}$）。
    - **Gaussian 核**：$K(x, y) = e^{-(x-y)^2}$（严格地说是全非负核）。
    - **Cauchy 核**：$K(x, y) = 1/(x + y)$（$x, y > 0$）。

!!! example "例 39.9 (指数核的全正性验证)"
    取 $x_1 < x_2$，$y_1 < y_2$，矩阵 $M = \begin{pmatrix} e^{x_1 y_1} & e^{x_1 y_2} \\ e^{x_2 y_1} & e^{x_2 y_2} \end{pmatrix}$。

    $$\det(M) = e^{x_1 y_1 + x_2 y_2} - e^{x_1 y_2 + x_2 y_1} = e^{x_1 y_1 + x_2 y_2}\left(1 - e^{(x_1 - x_2)(y_2 - y_1)}\right).$$

    由 $x_1 < x_2$，$y_1 < y_2$，$(x_1 - x_2)(y_2 - y_1) < 0$，故 $e^{(x_1-x_2)(y_2-y_1)} < 1$，$\det(M) > 0$。

    对一般 $n$，这是 Vandermonde 行列式的指数版本，可以通过类似的方法或由一般理论（指数函数构成 Chebyshev 系统）证明全正性。

!!! note "注记 39.6 (全正矩阵与概率论)"
    全正性在概率论中有重要应用。Karlin 证明了，如果随机过程的转移核是全正的，则该过程具有良好的"单调性"性质——状态之间的序关系在转移下保持。这类过程称为**随机单调的**（stochastically monotone），在排队论和可靠性理论中有广泛应用。

    特别地，生灭过程（birth-death process）的转移矩阵是全非负的（实际上是振荡矩阵），这使得生灭过程的谱分析具有定理 39.5-39.6 所描述的优美性质：特征值全部实、正、简单，且特征向量具有严格的振荡模式。

---

## 本章小结

| 概念 | 定义 | 关键性质 |
|------|------|----------|
| 全正 (TP) | 所有子式 $> 0$ | 特征值正、简单；变差减少 |
| 全非负 (TN) | 所有子式 $\ge 0$ | 双对角分解；平面网络 |
| 振荡矩阵 | TN + 某幂次为 TP | 非奇异 TN + 全不可约 |
| 变差减少 | $S^-(Ax) \le S^-(x)$ | 等价于 TN |
| 双对角分解 | TN = 初等双对角之积 | Loewner-Whitney 定理 |
| LGV 引理 | 路径矩阵的子式 = 不相交路径族 | 平面图 $\Rightarrow$ TN |
| TP 核 | $K(x,y)$，所有有限抽样矩阵 TP | 指数核、Cauchy 核 |
| Polya 频率 | 无穷 Toeplitz TN 矩阵 | 生成函数的解析刻画 |

全正矩阵理论将矩阵的代数性质（行列式的符号）、分析性质（特征值和特征向量的定性行为）、组合性质（平面网络、路径计数）编织成一个和谐的整体。

从逼近论中的 B-样条到物理学中的散射振幅，从概率论中的 Polya 频率序列到代数组合学中的 cluster algebra，全正性是一个不断被重新发现和深化的数学主题。它联结了分析、代数、组合和几何等多个数学分支，是线性代数最富有生命力的研究领域之一。

## 练习题

1. **[概念] 什么是全正矩阵（TP）？它与正矩阵（Positive Matrix）的主要区别是什么？**
   ??? success "参考答案"
       正矩阵仅要求矩阵元素全部为正（即一阶子式为正）；而全正矩阵要求**所有阶数**的子式（Minor）都必须为正。TP 是一个远强于正矩阵的约束。

2. **[基础] 验证 $A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$ 是否为全正矩阵。**
   ??? success "参考答案"
       1. 一阶子式：$1, 1, 1, 2$，全正。
       2. 二阶子式：$\det(A) = 2 - 1 = 1 > 0$。
       所有子式均为正，故 $A$ 是全正矩阵。

3. **[特征值] 根据 Gantmacher-Krein 定理，全正矩阵的特征值具有什么极端优美的性质？**
   ??? success "参考答案"
       全正矩阵的所有特征值都是**互不相同**的**严格正实数**：$\lambda_1 > \lambda_2 > \dots > \lambda_n > 0$。这在实矩阵中是非常罕见的强稳定性。

4. **[特征向量] 设 $A$ 是 $n$ 阶全正矩阵，$v^{(k)}$ 是对应于第 $k$ 大特征值的特征向量。$v^{(k)}$ 的分量符号变化次数是多少？**
   ??? success "参考答案"
       $v^{(k)}$ 的分量恰好有 $k-1$ 次符号变化。例如，主特征向量 $v^{(1)}$ 没有符号变化（全正或全负），而最小特征值的特征向量 $v^{(n)}$ 则会剧烈振荡 $n-1$ 次。

5. **[乘积] 证明：两个全正矩阵的乘积仍然是全正矩阵。**
   ??? success "参考答案"
       利用 **Cauchy-Binet 公式**：$\det((AB)[\alpha, \beta]) = \sum_{\gamma} \det(A[\alpha, \gamma]) \det(B[\gamma, \beta])$。由于右侧求和中的每一项都是两个正子式的乘积，结果必然为正。

6. **[变差减少] 什么是全非负矩阵（TN）的“变差减少”（Variation Diminishing）性质？它在物理上意味着什么？**
   ??? success "参考答案"
       性质：若 $y = Ax$，则 $y$ 的分量符号变化次数不会超过 $x$ 的分量符号变化次数（$S^-(y) \le S^-(x)$）。物理意义：全非负变换起到了“平滑滤波”的作用，它不会向系统中引入额外的振荡，反映了某种扩散或能量耗散过程。

7. **[Pascal] 证明 Pascal 矩阵 $P_{ij} = \binom{i+j-2}{i-1}$ 是全正的。**
   ??? success "参考答案"
       可以利用 Lindström-Gessel-Viennot 引理，将 Pascal 矩阵看作格点上不相交路径的计数矩阵。由于路径始终向右或向上，不可能交叉（平面性），其所有子式必然非负且可通过路径构造证明严格为正。

8. **[判定] 高效判定一个矩阵是否全正，最少需要检查多少个子式？**
   ??? success "参考答案"
       不需要检查所有指数级数量的子式。利用 Neville 消元法或特定的初始子式判据，只需要检查 $O(n^2)$ 个子式（如所有的实心子式/Contiguous Minors）即可。

9. **[振荡矩阵] 什么是振荡矩阵（Oscillatory Matrix）？它与 TP 有什么关系？**
   ??? success "参考答案"
       一个 TN 矩阵 $A$ 称为振荡矩阵，若其非奇异且其次对角线元素为正。根据 Gantmacher-Krein 定理，这种矩阵的某个幂次 $A^p$ 必定是全正矩阵。

10. **[爱因斯坦思考题] 爱因斯坦曾在《几何学与经验》中探讨物理世界的连续性。在全正核（TP Kernel）理论中，函数 $K(x,y)$ 被称为全正的，若对其任意有限采样的矩阵都是全正的。为什么这种性质被称为“离散与连续的完美桥梁”？它如何体现了自然界中“有序性”的跨尺度保持？**
    ??? success "参考答案"
        全正性要求在任何尺度（$n$ 阶子式）和任何采样点下都保持行列式的单调性（正号）。这意味着系统在微观采样和宏观整体上都遵循相同的“非交叠”规律（如不相交路径）。这种跨尺度的有序性使得我们可以通过有限维全正矩阵的性质，去推断连续物理介质（如热传导、概率分布）的解析性质。全正性是数学中的一种“全息”属性：局部的小子式决定了全局的谱振荡模式，完美契合了爱因斯坦追求的宇宙规律的统一性与自恰性。

## 本章小结

本章系统探讨了矩阵论中各阶子式皆具有确定符号模式的特殊矩阵类——全正矩阵，主要内容包括：

1. **基本定义与层级**：区分了全正（TP）与全非负（TN）的概念，确立了子式正性作为核心评判准则。
2. **谱理论的巅峰**：介绍了全正矩阵特征值互异且严格正、特征向量具有严格振荡模式（节点定理）的深刻性质，解释了其在物理振动模态分析中的基石作用。
3. **代数分解与参数化**：推导了双对角分解（Loewner-Whitney 定理），将复杂的 TP 矩阵简化为初等非负双对角矩阵的乘积，并揭示了其与平面网络路径计数的几何关联。
4. **分析与组合应用**：展示了全正性在 B-样条保形逼近、Polya 频率序列、以及现代代数组合学（如 Cluster 代数与正 Grassmannian）中的广泛渗透。
5. **计算稳定性**：论证了变差减少性质如何确保物理演化中的平滑性，并提供了高效判定全正性的 $O(n^3)$ 算法。

