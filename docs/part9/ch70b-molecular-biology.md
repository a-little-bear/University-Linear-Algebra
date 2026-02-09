# 第 70B 章 分子生物学与系统生物学中的线性代数

<div class="context-flow" markdown>

**前置**：特征值与特征向量(Ch6) · 矩阵指数(Ch13) · 奇异值分解(Ch7) · 正定矩阵(Ch16) · 线性规划(Ch65) · 概率与统计(Ch27)

**本章脉络**：序列比对评分矩阵(PAM, BLOSUM) $\to$ 系统发育树(UPGMA) $\to$ DNA 进化的 Markov 模型(Jukes-Cantor, Kimura, GTR) $\to$ 隐 Markov 模型(基因预测) $\to$ 基因组学中的 PCA $\to$ 化学计量矩阵与流量平衡分析 $\to$ 基因调控网络(邻接矩阵、Jacobi 稳定性) $\to$ Boolean 网络 $\to$ 网络推断(协方差/精度矩阵)

**延伸**：评分矩阵是 BLAST 等序列比对工具的理论基础；DNA 进化模型是分子系统发育学的核心；HMM 是基因预测和序列标注的标准方法；PCA 在人类遗传学中揭示了种群结构；流量平衡分析是代谢工程和合成生物学的重要工具；基因调控网络分析推动了系统生物学的发展

</div>

分子生物学和系统生物学是线性代数方法最富成果的应用领域之一。从蛋白质序列比对到基因组进化，从代谢网络到基因调控，矩阵方法渗透在每个层面。评分矩阵本质上是突变概率矩阵的对数；DNA 进化模型由 $4 \times 4$ 速率矩阵和矩阵指数刻画；隐 Markov 模型将基因预测转化为矩阵乘法；化学计量矩阵的零空间确定了代谢稳态；基因调控网络的稳定性由 Jacobi 矩阵的特征值决定。

本章系统展示这些联系，从分子层面到系统层面。

---

## 70B.1 序列比对评分矩阵

<div class="context-flow" markdown>

**核心问题**：比较两条蛋白质序列时，氨基酸之间的替换评分应该如何确定？评分矩阵的数学本质是什么？

</div>

!!! definition "定义 70B.1 (PAM 矩阵)"
    **PAM**（Point Accepted Mutation）矩阵由 Dayhoff 等人于 1978 年提出。设 $M$ 为 $20 \times 20$ 的氨基酸突变概率矩阵（$M_{ij}$ 为氨基酸 $i$ 在一个进化时间单位内突变为氨基酸 $j$ 的条件概率），$\mathbf{q}$ 为背景氨基酸频率。则 PAM-$k$ 评分矩阵定义为：

    $$S_{ij}^{(k)} = \log_2 \frac{(M^k)_{ij}}{q_j}$$

    其中 $M^k$ 是 $M$ 的 $k$ 次方（$k$ 个进化时间单位后的突变概率），$q_j$ 是随机配对下观察到氨基酸 $j$ 的概率。

    PAM-1 对应于 1% 的氨基酸被替换的进化距离；PAM-250 对应于约 80% 的位点已发生替换（多次替换叠加），适用于检测远缘同源序列。

!!! theorem "定理 70B.1 (评分矩阵的矩阵对数形式)"
    PAM-1 评分矩阵可以写为：

    $$S^{(1)} = \log_2 M - \mathbf{e}\log_2(\mathbf{q})^T$$

    PAM-$k$ 评分矩阵满足：

    $$S^{(k)} = \log_2 M^k - \mathbf{e}\log_2(\mathbf{q})^T$$

    其中对数逐元素取，$\mathbf{e}$ 为全 1 向量。

??? proof "证明"
    由定义，$S_{ij}^{(k)} = \log_2 (M^k)_{ij} - \log_2 q_j$。

    写成矩阵形式：$S^{(k)}$ 的 $(i,j)$ 元素为 $[\log_2 M^k]_{ij} - \log_2 q_j$。

    $\log_2 q_j$ 只依赖于列指标 $j$，因此可以写成 $\mathbf{e} \cdot [\log_2(\mathbf{q})^T]$，其中 $\mathbf{e} = (1, \ldots, 1)^T$。

    注意 $M^k$ 的计算可以通过 $M$ 的谱分解来实现：若 $M = P \Lambda P^{-1}$，则 $M^k = P \Lambda^k P^{-1}$。这对大 $k$ 的计算非常高效。

    $\blacksquare$

!!! definition "定义 70B.2 (BLOSUM 矩阵)"
    **BLOSUM**（BLOcks SUbstitution Matrix）由 Henikoff 和 Henikoff 于 1992 年直接从序列比对数据库中估计，不依赖于进化模型的外推。BLOSUM-$n$ 矩阵基于序列相似度阈值为 $n\%$ 的聚类：

    $$S_{ij} = \log_2 \frac{p_{ij}}{q_i q_j}$$

    其中 $p_{ij}$ 是在比对中观察到氨基酸对 $(i,j)$ 的频率，$q_i q_j$ 是随机情况下的期望频率。

    这是一个**对数优势比**（log-odds ratio）矩阵：$S_{ij} > 0$ 表示 $(i,j)$ 替换比随机预期更频繁（保守替换），$S_{ij} < 0$ 表示更少见。

!!! theorem "定理 70B.2 (评分矩阵的谱结构与比对能力)"
    有效的氨基酸评分矩阵 $S$ 应具有以下谱性质：

    1. $S$ 是对称矩阵（在适当加权下）
    2. $S$ 有少数几个大的正特征值，对应于保守的物理化学性质维度（如疏水性、电荷、分子大小）
    3. $S$ 的正特征值对应的特征向量揭示了氨基酸的物理化学聚类

??? proof "证明"
    定义加权评分矩阵 $\tilde{S}_{ij} = q_i^{1/2} S_{ij} q_j^{1/2} = q_i^{1/2} q_j^{1/2} \log_2(p_{ij}/(q_i q_j))$。

    由于 $p_{ij} = p_{ji}$（替换频率对称），$\tilde{S}$ 是对称矩阵，因此可以正交对角化。

    正特征值对应于"比随机更偏好的替换方向"。在 BLOSUM62 中，最大正特征值对应的特征向量与疏水性高度相关，第二大正特征值与分子大小相关。

    比对算法的区分能力取决于 $S$ 的**相对熵**（期望得分）：$H = \sum_{ij} p_{ij} S_{ij}$。有效矩阵要求 $H > 0$。$H$ 与 $S$ 的正特征值之和有关。

    $\blacksquare$

!!! example "例 70B.1"
    BLOSUM62 矩阵的部分条目（以 1/2 bit 为单位）：

    |   | A  | R  | D  | C  |
    |---|----|----|----|----|
    | A |  4 | -1 | -2 |  0 |
    | R | -1 |  5 | -2 | -3 |
    | D | -2 | -2 |  6 | -3 |
    | C |  0 | -3 | -3 |  9 |

    对角线元素（自身替换）总是最高的。半胱氨酸（C）的自替换得分最高（9），因为半胱氨酸高度保守（通过二硫键固定三维结构）。

---

## 70B.2 系统发育树

<div class="context-flow" markdown>

**核心问题**：如何从物种之间的遗传距离矩阵推断进化关系（系统发育树）？

</div>

!!! definition "定义 70B.3 (距离矩阵)"
    **距离矩阵** $D \in \mathbb{R}^{n \times n}$ 是对称非负矩阵，$D_{ij}$ 表示物种 $i$ 和 $j$ 之间的进化距离，$D_{ii} = 0$。若 $D$ 满足三角不等式 $D_{ij} \leq D_{ik} + D_{kj}$，则称 $D$ 为度量矩阵。

!!! definition "定义 70B.4 (超度量矩阵)"
    若对所有 $i, j, k$：

    $$D_{ij} \leq \max(D_{ik}, D_{kj})$$

    则 $D$ 称为**超度量矩阵**。超度量条件等价于：对任意三个物种，三个两两距离中最大的两个相等。

    超度量矩阵对应于**分子钟假设**下的系统发育树（所有叶节点到根的距离相等，即所有谱系以相同速率进化）。

!!! theorem "定理 70B.3 (UPGMA 算法的正确性)"
    若真实距离矩阵 $D$ 是超度量的，则 UPGMA（Unweighted Pair Group Method with Arithmetic Mean）算法能正确重建系统发育树。

??? proof "证明"
    UPGMA 算法的步骤：

    1. 找到 $D$ 中最小的非对角元素 $D_{ij}$
    2. 将物种 $i$ 和 $j$ 合并为一个簇，分支长度为 $D_{ij}/2$
    3. 更新距离矩阵：新簇 $\{i,j\}$ 与其他物种 $k$ 的距离为 $D_{\{i,j\},k} = (D_{ik} + D_{jk})/2$
    4. 重复直到所有物种合并

    **正确性证明**：需要证明在超度量条件下，$D$ 中距离最小的两个物种在真实树中是最近邻（姊妹类群）。

    反证法：设 $D_{ij}$ 是最小距离，但 $i$ 的真正最近邻是 $k \neq j$。在超度量树中，$i$ 和 $k$ 的最近公共祖先比 $i$ 和 $j$ 的更近，因此 $D_{ik} < D_{ij}$，这与 $D_{ij}$ 是最小距离矛盾。

    合并后距离的更新保持超度量性质：设 $D$ 为超度量，$D_{ij}$ 为最小距离。对任意 $k, l$（$k, l \notin \{i,j\}$）：

    - $D_{\{i,j\},k} = (D_{ik} + D_{jk})/2$。由超度量性，$D_{ik} = D_{jk}$（因为 $D_{ij} \leq \max(D_{ik}, D_{jk})$，且 $D_{ij}$ 最小，所以 $D_{ik} = D_{jk}$）。故 $D_{\{i,j\},k} = D_{ik} = D_{jk}$。
    - 新距离矩阵仍然是超度量的（可用归纳法验证三角不等式的超度量版本）。

    因此算法的每一步都正确选择了最近邻，递归正确。

    $\blacksquare$

!!! definition "定义 70B.5 (邻接法)"
    **邻接法**（Neighbor-Joining，NJ）不要求超度量条件，适用于进化速率不均匀的情形。对每个物种 $i$，定义：

    $$u_i = \frac{1}{n-2} \sum_{k=1}^{n} D_{ik}$$

    选择使 $D_{ij} - u_i - u_j$ 最小的配对 $(i,j)$。分支长度为：

    $$v_i = \frac{1}{2}(D_{ij} + u_i - u_j), \quad v_j = D_{ij} - v_i$$

!!! example "例 70B.2"
    四个物种的距离矩阵：

    $$D = \begin{pmatrix} 0 & 2 & 4 & 4 \\ 2 & 0 & 4 & 4 \\ 4 & 4 & 0 & 2 \\ 4 & 4 & 2 & 0 \end{pmatrix}$$

    验证超度量性：对任意三元组，最大的两个距离相等。例如 $(1,2,3)$：$D_{12}=2, D_{13}=4, D_{23}=4$，最大两个相等。

    UPGMA 首先合并最近的一对。$D_{12} = D_{34} = 2$ 是最小距离，合并 $(1,2)$ 和 $(3,4)$，分支长度为 1。更新距离后，两个簇之间的距离为 4，最终树为 $((1,2),(3,4))$，内部分支长度为 $4/2 - 1 = 1$。

---

## 70B.3 DNA 进化的 Markov 模型

<div class="context-flow" markdown>

**核心问题**：DNA 序列在进化过程中的碱基替换如何用矩阵指数来建模？不同的进化模型对应什么样的速率矩阵？

</div>

DNA 由四种碱基组成：A（腺嘌呤）、C（胞嘧啶）、G（鸟嘌呤）、T（胸腺嘧啶）。在进化过程中，碱基会发生替换。连续时间 Markov 链提供了描述这一过程的自然框架。

!!! definition "定义 70B.6 (核苷酸替换的速率矩阵)"
    碱基替换过程由 $4 \times 4$ **速率矩阵**（rate matrix）$Q$ 描述：

    - $Q_{ij} \geq 0$（$i \neq j$）：从碱基 $i$ 到碱基 $j$ 的瞬时替换速率
    - $Q_{ii} = -\sum_{j \neq i} Q_{ij}$：行和为 0

    经过时间 $t$ 后的**替换概率矩阵**为矩阵指数：

    $$P(t) = e^{Qt}$$

    $P_{ij}(t)$ 是经过 $t$ 个进化时间单位后碱基 $i$ 变为碱基 $j$ 的概率。

!!! definition "定义 70B.7 (Jukes-Cantor 模型)"
    **Jukes-Cantor 模型**（JC69）是最简单的核苷酸替换模型，假设所有替换等概率：

    $$Q_{JC} = \alpha \begin{pmatrix} -3 & 1 & 1 & 1 \\ 1 & -3 & 1 & 1 \\ 1 & 1 & -3 & 1 \\ 1 & 1 & 1 & -3 \end{pmatrix} = \alpha(J - 4I)$$

    其中 $J$ 为全 1 矩阵，$\alpha > 0$ 为替换速率参数。

!!! theorem "定理 70B.4 (Jukes-Cantor 转移概率)"
    Jukes-Cantor 模型的转移概率矩阵为：

    $$P(t) = e^{Q_{JC}t} = \begin{pmatrix} p_0(t) & p_1(t) & p_1(t) & p_1(t) \\ p_1(t) & p_0(t) & p_1(t) & p_1(t) \\ p_1(t) & p_1(t) & p_0(t) & p_1(t) \\ p_1(t) & p_1(t) & p_1(t) & p_0(t) \end{pmatrix}$$

    其中：

    $$p_0(t) = \frac{1}{4} + \frac{3}{4}e^{-4\alpha t}, \quad p_1(t) = \frac{1}{4} - \frac{1}{4}e^{-4\alpha t}$$

    Jukes-Cantor 进化距离为 $d = -\frac{3}{4}\ln(1 - \frac{4}{3}p)$，其中 $p$ 为观察到的替换比例。

??? proof "证明"
    $Q_{JC} = \alpha(J - 4I)$。$J$ 的特征值为 4（对应特征向量 $\mathbf{e} = (1,1,1,1)^T$）和 0（三重，对应正交补空间）。

    因此 $Q_{JC}$ 的特征值为 $\alpha(4 - 4) = 0$（单重）和 $\alpha(0 - 4) = -4\alpha$（三重）。

    矩阵指数：

    $$e^{Q_{JC}t} = \frac{1}{4}\mathbf{e}\mathbf{e}^T \cdot e^{0 \cdot t} + \left(I - \frac{1}{4}\mathbf{e}\mathbf{e}^T\right) \cdot e^{-4\alpha t}$$

    $$= \frac{1}{4}J + \left(I - \frac{1}{4}J\right)e^{-4\alpha t}$$

    对角元素：$P_{ii}(t) = 1/4 + 3/4 \cdot e^{-4\alpha t} = p_0(t)$。

    非对角元素：$P_{ij}(t) = 1/4 - 1/4 \cdot e^{-4\alpha t} = p_1(t)$。

    进化距离 $d = 3\alpha t$（每个位点期望的总替换数），观察到的替换比例 $p = 3p_1(t) = 3/4(1 - e^{-4\alpha t})$。解出 $e^{-4\alpha t} = 1 - 4p/3$，故 $d = 3\alpha t = -3/4 \cdot \ln(1 - 4p/3)$。

    $\blacksquare$

!!! definition "定义 70B.8 (Kimura 双参数模型)"
    **Kimura 模型**（K80）区分转换（transition，嘌呤-嘌呤或嘧啶-嘧啶替换 A$\leftrightarrow$G, C$\leftrightarrow$T）和颠换（transversion，嘌呤-嘧啶替换），因为转换在自然界中更频繁：

    $$Q_K = \begin{pmatrix} * & \beta & \alpha & \beta \\ \beta & * & \beta & \alpha \\ \alpha & \beta & * & \beta \\ \beta & \alpha & \beta & * \end{pmatrix}$$

    其中 $\alpha$ 为转换速率，$\beta$ 为颠换速率（$\alpha > \beta$），$*$ 表示使行和为 0 的值（$= -\alpha - 2\beta$）。碱基顺序为 A, C, G, T。

!!! definition "定义 70B.9 (GTR 模型)"
    **GTR**（General Time Reversible）模型是最一般的时间可逆核苷酸替换模型：

    $$Q_{GTR} = \begin{pmatrix} * & a\pi_C & b\pi_G & c\pi_T \\ a\pi_A & * & d\pi_G & e\pi_T \\ b\pi_A & d\pi_C & * & f\pi_T \\ c\pi_A & e\pi_C & f\pi_G & * \end{pmatrix}$$

    其中 $\pi_A, \pi_C, \pi_G, \pi_T$ 为平衡碱基频率（$\sum \pi_i = 1$），$a, b, c, d, e, f$ 为 6 个交换速率参数。

    **时间可逆条件**（detailed balance）：$\pi_i Q_{ij} = \pi_j Q_{ji}$，即 $Q$ 关于 $\text{diag}(\boldsymbol{\pi})$ 可逆。

!!! theorem "定理 70B.5 (GTR 速率矩阵的对称化与特征值)"
    定义 $\Pi = \text{diag}(\pi_A, \pi_C, \pi_G, \pi_T)$。GTR 速率矩阵 $Q$ 满足 $\Pi Q = Q^T \Pi$（时间可逆性）。因此：

    1. $\Pi^{1/2} Q \Pi^{-1/2}$ 是对称矩阵
    2. $Q$ 的特征值全为实数
    3. $Q$ 有一个零特征值（对应平衡分布 $\boldsymbol{\pi}$），其余特征值为负

??? proof "证明"
    **(1)** 时间可逆性 $\pi_i Q_{ij} = \pi_j Q_{ji}$ 意味着 $\Pi Q = (\Pi Q)^T$，即 $\Pi Q$ 是对称矩阵。令 $\tilde{Q} = \Pi^{1/2} Q \Pi^{-1/2}$，则 $\tilde{Q} = \Pi^{-1/2}(\Pi Q)\Pi^{-1/2}$，这是对称矩阵 $\Pi Q$ 的合同变换，仍然对称。

    **(2)** 对称矩阵的特征值全为实数。$\tilde{Q}$ 与 $Q$ 相似（$\tilde{Q} = \Pi^{1/2} Q \Pi^{-1/2}$），故特征值相同。

    **(3)** $Q\boldsymbol{\pi} = 0$（$\boldsymbol{\pi}$ 是平衡分布，因此 $Q$ 的行和为 0 意味着 $Q\mathbf{e} = 0$...更准确地说，$\boldsymbol{\pi}^T Q = 0$ 因为 $\boldsymbol{\pi}$ 是左特征向量）。实际上，$Q$ 的行和为 0 意味着 $\mathbf{e}$ 是右零空间的元素：$Q\mathbf{e} = 0$（不是 $\boldsymbol{\pi}$）。而 $\boldsymbol{\pi}^T$ 是左零空间：$\boldsymbol{\pi}^T Q = 0^T$（这由 detailed balance 得出）。

    其余特征值为负：$\tilde{Q}$ 是行和非正的对称矩阵（对角占优），由 Gershgorin 定理，特征值 $\leq 0$。零特征值简单（因为 $Q$ 不可约），故其余特征值严格为负。

    $\blacksquare$

!!! example "例 70B.3"
    Kimura 模型中，设转换/颠换比 $\kappa = \alpha/\beta = 2$，$\beta = 1$，$\alpha = 2$。

    $$Q_K = \begin{pmatrix} -4 & 1 & 2 & 1 \\ 1 & -4 & 1 & 2 \\ 2 & 1 & -4 & 1 \\ 1 & 2 & 1 & -4 \end{pmatrix}$$

    特征值：$0, -4, -4, -8$（注意 $-4$ 是二重的）。

    转移概率：$P(t) = e^{Qt}$ 的对角元素为 $p_0(t) = 1/4 + 1/4 e^{-4t} + 1/2 e^{-8t}$。

---

## 70B.4 隐 Markov 模型在生物信息学中的应用

<div class="context-flow" markdown>

**核心问题**：如何利用隐 Markov 模型（HMM）进行基因预测、序列标注和蛋白质家族建模？其核心算法如何归结为矩阵运算？

</div>

隐 Markov 模型是生物信息学中最重要的统计工具之一。基因预测（区分编码区和非编码区）、蛋白质家族建模（profile HMM）和多序列比对都依赖 HMM。

!!! definition "定义 70B.10 (隐 Markov 模型)"
    一个**隐 Markov 模型** $\lambda = (\mathbf{A}, \mathbf{B}, \boldsymbol{\pi})$ 由以下要素构成：

    - **隐状态集** $\{1, 2, \ldots, N\}$：如"外显子"、"内含子"、"基因间区"
    - **观测符号集** $\{o_1, \ldots, o_M\}$：如 DNA 碱基 $\{A, C, G, T\}$
    - **转移矩阵** $\mathbf{A} \in \mathbb{R}^{N \times N}$：$A_{ij} = P(\text{时刻 } t+1 \text{ 处于状态 } j \mid \text{时刻 } t \text{ 处于状态 } i)$
    - **发射矩阵** $\mathbf{B} \in \mathbb{R}^{N \times M}$：$B_{ik} = P(\text{观测到 } o_k \mid \text{隐状态为 } i)$
    - **初始分布** $\boldsymbol{\pi} \in \mathbb{R}^N$：$\pi_i = P(\text{初始状态为 } i)$

!!! theorem "定理 70B.6 (前向算法与矩阵乘法)"
    给定观测序列 $O = (o_{t_1}, o_{t_2}, \ldots, o_{t_T})$，观测似然 $P(O \mid \lambda)$ 可以通过**前向算法**计算，其核心是矩阵乘法：

    定义对角发射矩阵 $D_t = \text{diag}(B_{1,t_t}, B_{2,t_t}, \ldots, B_{N,t_t})$（观测到符号 $o_{t_t}$ 时各状态的发射概率）。

    则：

    $$P(O \mid \lambda) = \boldsymbol{\pi}^T D_1 \mathbf{A} D_2 \mathbf{A} D_3 \cdots \mathbf{A} D_T \mathbf{e}$$

    其中 $\mathbf{e} = (1, \ldots, 1)^T$。

    计算复杂度为 $O(N^2 T)$（而暴力枚举所有状态路径的复杂度为 $O(N^T)$）。

??? proof "证明"
    定义前向变量 $\alpha_t(i) = P(o_{t_1}, \ldots, o_{t_t}, \text{时刻 } t \text{ 状态} = i \mid \lambda)$。

    **初始化**：$\alpha_1(i) = \pi_i B_{i,t_1}$，向量形式 $\boldsymbol{\alpha}_1 = D_1 \boldsymbol{\pi}$（逐元素乘）。

    **递推**：$\alpha_{t+1}(j) = \left[\sum_{i=1}^N \alpha_t(i) A_{ij}\right] B_{j,t_{t+1}}$

    向量形式：$\boldsymbol{\alpha}_{t+1} = D_{t+1} \mathbf{A}^T \boldsymbol{\alpha}_t$。

    **终止**：$P(O \mid \lambda) = \sum_i \alpha_T(i) = \mathbf{e}^T \boldsymbol{\alpha}_T$。

    展开递推：$\boldsymbol{\alpha}_T = D_T \mathbf{A}^T D_{T-1} \mathbf{A}^T \cdots D_2 \mathbf{A}^T D_1 \boldsymbol{\pi}$。

    因此 $P(O \mid \lambda) = \mathbf{e}^T D_T \mathbf{A}^T D_{T-1} \cdots \mathbf{A}^T D_1 \boldsymbol{\pi} = \boldsymbol{\pi}^T D_1 \mathbf{A} D_2 \cdots \mathbf{A} D_T \mathbf{e}$。

    每步矩阵-向量乘法 $O(N^2)$，共 $T$ 步，总复杂度 $O(N^2 T)$。

    $\blacksquare$

!!! definition "定义 70B.11 (Viterbi 算法——最优路径解码)"
    **Viterbi 算法**寻找使 $P(O, \text{path} \mid \lambda)$ 最大化的隐状态路径。其结构与前向算法类似，但将求和替换为取最大值：

    $$\delta_{t+1}(j) = \max_i [\delta_t(i) \cdot A_{ij}] \cdot B_{j,t_{t+1}}$$

    在对数空间中，这变为**最大值加法**（max-plus 代数）下的矩阵乘法：

    $$\log \delta_{t+1}(j) = \max_i [\log \delta_t(i) + \log A_{ij}] + \log B_{j,t_{t+1}}$$

    这是 $(\max, +)$ 半环上的矩阵运算。

!!! example "例 70B.4"
    **简化的基因预测 HMM**。两个隐状态：$E$（外显子）和 $I$（内含子），观测符号为 $\{A, C, G, T\}$。

    转移矩阵：$\mathbf{A} = \begin{pmatrix} 0.95 & 0.05 \\ 0.1 & 0.9 \end{pmatrix}$（外显子倾向于保持，内含子也倾向于保持）

    发射矩阵：

    $$\mathbf{B} = \begin{pmatrix} 0.25 & 0.25 & 0.25 & 0.25 \\ 0.40 & 0.10 & 0.10 & 0.40 \end{pmatrix}$$

    外显子中碱基均匀分布（GC 含量约 50%），内含子中 AT 富集。

    对给定的 DNA 序列 ATGCAT...，前向算法计算总似然，Viterbi 算法给出最可能的外显子/内含子标注。

---

## 70B.5 基因组学中的主成分分析

<div class="context-flow" markdown>

**核心问题**：如何利用 PCA 从基因型数据中揭示人群的遗传结构？

</div>

2006 年，Patterson, Price 和 Reich 的 EIGENSOFT 方法展示了 PCA 在人类遗传学中的强大应用。

!!! definition "定义 70B.12 (基因型矩阵与 PCA)"
    设 $X \in \mathbb{R}^{n \times p}$ 为**基因型矩阵**，其中 $n$ 为个体数，$p$ 为 SNP（单核苷酸多态性）标记数。$X_{ij} \in \{0, 1, 2\}$ 表示第 $i$ 个个体在第 $j$ 个 SNP 位点的等位基因计数。

    **标准化**：$\tilde{X}_{ij} = (X_{ij} - 2p_j) / \sqrt{2p_j(1-p_j)}$，其中 $p_j$ 为第 $j$ 个 SNP 的等位基因频率。

    **PCA**：对 $\tilde{X}$ 进行奇异值分解 $\tilde{X} = U\Sigma V^T$，或等价地对 $n \times n$ 遗传相关矩阵

    $$G = \frac{1}{p}\tilde{X}\tilde{X}^T$$

    进行特征值分解。

!!! theorem "定理 70B.7 (PCA 揭示种群结构)"
    在 $F_{ST}$ 模型（Wright 的群体遗传结构模型）下，若 $n$ 个个体来自 $K$ 个不同的亚群，则遗传相关矩阵 $G$ 的前 $K-1$ 个主成分近似地将个体按亚群分离。

    具体地，设 $\boldsymbol{\pi}_i$ 为第 $i$ 个个体的祖先成分向量（$\boldsymbol{\pi}_i \in \mathbb{R}^K$，$\sum_k \pi_{ik} = 1$），则 $G \approx \Pi F \Pi^T$，其中 $\Pi$ 为 $n \times K$ 的祖先成分矩阵，$F$ 为 $K \times K$ 的 $F_{ST}$ 矩阵。

??? proof "证明"
    在 Balding-Nichols 模型下，SNP 位点 $j$ 在亚群 $k$ 中的等位基因频率为 $p_{jk}$，满足 $E[p_{jk}] = p_j$ 和 $\text{Var}(p_{jk}) = F_{ST,k} p_j(1-p_j)$。

    个体 $i$ 的标准化基因型向量为 $\tilde{\mathbf{x}}_i \in \mathbb{R}^p$。其期望（在给定祖先成分 $\boldsymbol{\pi}_i$ 下）取决于祖先成分。

    遗传相关矩阵的 $(i,i')$ 元素为：

    $$G_{ii'} = \frac{1}{p} \tilde{\mathbf{x}}_i^T \tilde{\mathbf{x}}_{i'} \approx \sum_{k=1}^{K} \pi_{ik} \pi_{i'k} F_{ST,k}$$

    这正是 $(\Pi F \Pi^T)_{ii'}$（当 $F$ 为对角 $F_{ST}$ 矩阵时）。

    $\Pi F \Pi^T$ 的秩至多为 $K$，因此其前 $K$ 个特征值非零（实际上由于约束 $\sum_k \pi_{ik} = 1$，有效秩为 $K-1$）。

    因此 $G$ 的前 $K-1$ 个主成分捕获了种群结构的主要信息。

    $\blacksquare$

!!! example "例 70B.5"
    **人类遗传学中的经典应用**。对来自欧洲不同国家的 3,000 个个体进行约 500,000 个 SNP 的基因分型，PCA 的前两个主成分（PC1 和 PC2）能够几乎完美地重建欧洲地理地图——遗传变异的主要方向与地理距离高度相关。

    这是因为人类迁徙的隔离-距离效应导致了连续的种群结构，PCA 的前两个主成分分别捕获了南北和东西方向的遗传梯度。

    在实际应用中，PCA 主成分常用作 GWAS（全基因组关联研究）中的混杂因素校正协变量，以避免种群分层导致的假阳性。

---

## 70B.6 化学计量矩阵与流量平衡分析

<div class="context-flow" markdown>

**核心问题**：如何用线性代数建模细胞内的代谢网络？代谢稳态的数学结构是什么？

</div>

代谢网络是细胞内化学反应的集合。流量平衡分析（Flux Balance Analysis, FBA）是系统生物学和代谢工程中最重要的建模方法之一。

!!! definition "定义 70B.13 (化学计量矩阵)"
    一个代谢网络有 $m$ 种代谢物和 $n$ 个反应。**化学计量矩阵** $S \in \mathbb{R}^{m \times n}$ 定义为：

    - $S_{ij} > 0$：反应 $j$ 产生代谢物 $i$（产物）
    - $S_{ij} < 0$：反应 $j$ 消耗代谢物 $i$（底物）
    - $S_{ij} = 0$：反应 $j$ 不涉及代谢物 $i$

    代谢物浓度向量 $\mathbf{x} \in \mathbb{R}^m$ 的动力学方程为：

    $$\frac{d\mathbf{x}}{dt} = S \mathbf{v}$$

    其中 $\mathbf{v} \in \mathbb{R}^n$ 为**通量向量**（flux vector），$v_j$ 为第 $j$ 个反应的速率。

!!! theorem "定理 70B.8 (稳态通量空间)"
    在代谢稳态（$d\mathbf{x}/dt = 0$）下，通量向量 $\mathbf{v}$ 满足：

    $$S\mathbf{v} = 0$$

    即稳态通量空间是化学计量矩阵 $S$ 的**零空间**（null space）：

    $$\mathcal{V}_{\text{steady}} = \text{Null}(S) = \{\mathbf{v} \in \mathbb{R}^n : S\mathbf{v} = 0\}$$

    该空间的维数为 $n - \text{rank}(S)$，即反应数减去独立代谢平衡方程数。

??? proof "证明"
    稳态条件直接给出 $S\mathbf{v} = 0$。

    $\text{Null}(S)$ 是 $\mathbb{R}^n$ 的子空间，维数由秩-零度定理给出：

    $$\dim(\text{Null}(S)) = n - \text{rank}(S)$$

    在典型的代谢网络中，$n > m$（反应数多于代谢物数），且 $S$ 的秩通常接近 $m$（大多数代谢平衡约束是独立的），因此零空间维数约为 $n - m > 0$——即存在一个多维的可行通量空间。

    $\blacksquare$

!!! definition "定义 70B.14 (流量平衡分析)"
    **流量平衡分析**（FBA）在稳态通量空间中，通过线性规划找到最优通量分布：

    $$\max_{\mathbf{v}} \quad \mathbf{c}^T \mathbf{v}$$

    $$\text{s.t.} \quad S\mathbf{v} = 0, \quad \mathbf{v}_{\min} \leq \mathbf{v} \leq \mathbf{v}_{\max}$$

    其中：

    - $\mathbf{c}^T \mathbf{v}$ 为目标函数（通常为生物质合成通量，即细胞生长速率）
    - $S\mathbf{v} = 0$：稳态约束
    - $\mathbf{v}_{\min}, \mathbf{v}_{\max}$：通量上下界（如不可逆反应 $v_j \geq 0$，底物摄取速率有上限）

!!! example "例 70B.6"
    **简化的中心碳代谢网络**。3 种代谢物（A, B, C），5 个反应：

    - $R_1$：$\to$ A（底物摄取）
    - $R_2$：A $\to$ B
    - $R_3$：A $\to$ C
    - $R_4$：B $\to$（产物输出）
    - $R_5$：C $\to$（生物质合成）

    化学计量矩阵：

    $$S = \begin{pmatrix} 1 & -1 & -1 & 0 & 0 \\ 0 & 1 & 0 & -1 & 0 \\ 0 & 0 & 1 & 0 & -1 \end{pmatrix}$$

    $\text{rank}(S) = 3$，$\dim(\text{Null}(S)) = 5 - 3 = 2$。

    稳态约束：$v_1 = v_2 + v_3$，$v_2 = v_4$，$v_3 = v_5$。

    以 $v_2, v_3$ 为自由变量，其余通量由稳态条件唯一确定。

    FBA：最大化 $v_5$（生物质），约束 $v_1 \leq 10$（底物摄取限制），$v_j \geq 0$。最优解为 $v_3 = v_5 = 10$，$v_2 = v_4 = 0$——将所有底物导向生物质合成。

---

## 70B.7 基因调控网络

<div class="context-flow" markdown>

**核心问题**：基因之间的调控关系（激活/抑制）如何用矩阵来表示和分析？网络的稳态有什么生物学意义？

</div>

!!! definition "定义 70B.15 (线性基因调控网络)"
    在简化模型中，$n$ 个基因的表达水平 $\mathbf{x}(t) = (x_1(t), \ldots, x_n(t))^T$ 满足线性微分方程：

    $$\dot{\mathbf{x}}(t) = W\mathbf{x}(t) + \mathbf{b}$$

    其中 $W \in \mathbb{R}^{n \times n}$ 为**调控矩阵**（$W_{ij} > 0$ 表示基因 $j$ 激活基因 $i$，$W_{ij} < 0$ 表示抑制），$\mathbf{b}$ 为基础表达率。

!!! theorem "定理 70B.9 (基因网络的稳定性)"
    若调控矩阵 $W$ 的所有特征值具有负实部（$W$ 是 Hurwitz 稳定的），则系统有唯一稳态：

    $$\mathbf{x}^* = -W^{-1}\mathbf{b}$$

    且 $\mathbf{x}^*$ 是全局渐近稳定的。

??? proof "证明"
    稳态条件 $\dot{\mathbf{x}} = 0$ 给出 $W\mathbf{x}^* + \mathbf{b} = 0$，即 $\mathbf{x}^* = -W^{-1}\mathbf{b}$（$W$ 可逆，因为 0 不是特征值）。

    令 $\mathbf{y} = \mathbf{x} - \mathbf{x}^*$，则 $\dot{\mathbf{y}} = W\mathbf{y}$，解为 $\mathbf{y}(t) = e^{Wt}\mathbf{y}(0)$。

    由于 $W$ 的所有特征值 $\lambda_k$ 实部为负，$e^{Wt}$ 的每一项 $e^{\lambda_k t}$ 都趋于 0：

    $$\|e^{Wt}\| \leq C e^{-\alpha t} \to 0 \quad (t \to \infty)$$

    其中 $\alpha = \min_k |\text{Re}(\lambda_k)| > 0$。故 $\mathbf{y}(t) \to 0$，即 $\mathbf{x}(t) \to \mathbf{x}^*$。

    这对任意初始条件 $\mathbf{x}(0)$ 成立，因此稳态是全局渐近稳定的。

    $\blacksquare$

!!! theorem "定理 70B.10 (非线性网络的局部稳定性——Jacobi 矩阵分析)"
    对更一般的非线性基因调控网络 $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x})$，稳态 $\mathbf{x}^*$（$\mathbf{f}(\mathbf{x}^*) = 0$）的局部稳定性由 Jacobi 矩阵决定：

    $$J = \left(\frac{\partial f_i}{\partial x_j}\right)\bigg|_{\mathbf{x}^*}$$

    若 $J$ 的所有特征值实部为负，则 $\mathbf{x}^*$ 局部渐近稳定（Hartman-Grobman 定理）。若存在正实部特征值，则不稳定。

??? proof "证明"
    在 $\mathbf{x}^*$ 附近，$\mathbf{f}(\mathbf{x}) = \mathbf{f}(\mathbf{x}^*) + J(\mathbf{x} - \mathbf{x}^*) + O(\|\mathbf{x} - \mathbf{x}^*\|^2)$。

    令 $\mathbf{y} = \mathbf{x} - \mathbf{x}^*$，线性化系统为 $\dot{\mathbf{y}} = J\mathbf{y}$。

    Hartman-Grobman 定理保证：若 $J$ 没有纯虚特征值（双曲稳态），则非线性系统在 $\mathbf{x}^*$ 附近的定性行为与线性化系统相同。

    因此，$J$ 的所有特征值实部为负 $\implies$ $\mathbf{x}^*$ 局部渐近稳定。

    $\blacksquare$

!!! example "例 70B.7"
    **Repressilator（抑制振荡器）**。三个基因形成环状抑制网络：

    $$x_1(t+1) = \overline{x_3(t)}, \quad x_2(t+1) = \overline{x_1(t)}, \quad x_3(t+1) = \overline{x_2(t)}$$

    用调控矩阵表示线性化版本：

    $$W = \begin{pmatrix} -1 & 0 & -\alpha \\ -\alpha & -1 & 0 \\ 0 & -\alpha & -1 \end{pmatrix}$$

    其中 $\alpha > 0$ 为抑制强度，对角线 $-1$ 为自降解。

    $W$ 的特征方程：$(\lambda + 1)^3 + \alpha^3 = 0$，故 $\lambda + 1 = -\alpha \omega^k$（$k = 0, 1, 2$，$\omega = e^{2\pi i/3}$）。

    特征值 $\lambda_k = -1 - \alpha \omega^k$：

    - $\lambda_0 = -1 - \alpha$：总是实负
    - $\lambda_1 = -1 - \alpha e^{2\pi i/3} = -1 + \alpha/2 + i\alpha\sqrt{3}/2$
    - $\lambda_2 = \overline{\lambda_1}$

    $\text{Re}(\lambda_1) = -1 + \alpha/2$。当 $\alpha > 2$ 时，$\text{Re}(\lambda_1) > 0$，稳态不稳定，对应于振荡行为——这正是 Repressilator 的设计原理。

---

## 70B.8 Boolean 网络

<div class="context-flow" markdown>

**核心问题**：基因调控网络的离散模型——Boolean 网络——与线性代数有什么联系？

</div>

!!! definition "定义 70B.16 (Boolean 网络模型)"
    在 **Boolean 网络**中，每个基因的状态为 $x_i \in \{0, 1\}$（关闭/开启），更新规则为布尔函数：

    $$x_i(t+1) = f_i(x_1(t), x_2(t), \ldots, x_n(t))$$

    网络的状态空间有 $2^n$ 个状态。网络最终会进入**吸引子**——固定点（稳态）或极限环（振荡）。

!!! theorem "定理 70B.11 (Boolean 网络的吸引子存在性)"
    一个 $n$ 基因的 Boolean 网络有有限的状态空间 $\{0,1\}^n$（$2^n$ 个状态），因此其动力学最终必进入长度有限的吸引子。

    - **固定点**（周期 1）对应于网络的稳定表达模式（如特定的细胞类型）
    - **极限环**（周期 $> 1$）对应于振荡表达模式（如细胞周期）

??? proof "证明"
    Boolean 网络的状态转移函数 $F: \{0,1\}^n \to \{0,1\}^n$ 定义了有限集合上的确定性映射。

    考虑轨道 $\mathbf{x}_0, \mathbf{x}_1 = F(\mathbf{x}_0), \mathbf{x}_2 = F(\mathbf{x}_1), \ldots$。由鸽巢原理，在 $2^n + 1$ 步内必有重复：$\mathbf{x}_i = \mathbf{x}_j$（$i < j$）。此后轨道以周期 $j - i$ 循环。

    **与线性代数的联系**：在 $\mathbb{F}_2$（二元域）上，若布尔函数是线性的（$f_i = \bigoplus_{j \in S_i} x_j$），则网络动态可以写成 $\mathbf{x}(t+1) = W \mathbf{x}(t) \pmod{2}$，其中 $W \in \mathbb{F}_2^{n \times n}$。此时吸引子的周期整除 $W$ 在 $\mathbb{F}_2$ 上的阶（即满足 $W^k = I$ 的最小正整数 $k$），固定点恰好是 $(W - I)\mathbf{x} = 0$ 在 $\mathbb{F}_2$ 上的解空间，其大小为 $2^{n - \text{rank}_{\mathbb{F}_2}(W - I)}$。

    $\blacksquare$

!!! example "例 70B.8"
    **线性 Boolean 网络的例子**。$n = 3$，更新规则 $\mathbf{x}(t+1) = W\mathbf{x}(t) \pmod{2}$，其中

    $$W = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix} \pmod{2}$$

    $W$ 是循环置换矩阵。$W^3 = I \pmod{2}$，因此所有轨道的周期整除 3。

    固定点：$(W - I)\mathbf{x} = 0 \pmod{2}$，即 $\begin{pmatrix} 1 & 1 & 0 \\ 0 & 1 & 1 \\ 1 & 0 & 1 \end{pmatrix}\mathbf{x} = 0 \pmod{2}$。

    该矩阵的 $\mathbb{F}_2$-秩为 2，零空间维数为 1，给出 2 个固定点：$(0,0,0)$ 和 $(1,1,1)$。其余 6 个状态组成两个长度为 3 的极限环。

---

## 70B.9 网络推断：协方差与精度矩阵

<div class="context-flow" markdown>

**核心问题**：如何从基因表达数据推断基因之间的调控关系？精度矩阵（逆协方差矩阵）有什么特殊意义？

</div>

!!! definition "定义 70B.17 (从表达数据推断网络)"
    给定 $m$ 个条件下 $n$ 个基因的表达数据 $X \in \mathbb{R}^{n \times m}$，推断调控矩阵 $W$。

    **回归方法**：假设 $\dot{\mathbf{x}} \approx (\mathbf{x}_{t+1} - \mathbf{x}_t)/\Delta t$，则：

    $$\frac{X_{\text{next}} - X_{\text{curr}}}{\Delta t} \approx W X_{\text{curr}} + \mathbf{b}\mathbf{e}^T$$

    令 $\dot{X} = (X_{\text{next}} - X_{\text{curr}})/\Delta t$，最小二乘解为：

    $$\hat{W} = \dot{X} X_{\text{curr}}^T (X_{\text{curr}} X_{\text{curr}}^T)^{-1}$$

!!! definition "定义 70B.18 (精度矩阵与条件独立性)"
    对稳态基因表达数据，设 $\mathbf{x} \sim N(\boldsymbol{\mu}, \Sigma)$，**精度矩阵**定义为：

    $$\Theta = \Sigma^{-1}$$

    精度矩阵的关键性质：**$\Theta_{ij} = 0$ 当且仅当 $x_i$ 和 $x_j$ 在给定所有其他基因的表达水平后条件独立**。

    因此，精度矩阵的非零模式直接给出了基因调控网络的**图结构**——边 $(i,j)$ 存在当且仅当 $\Theta_{ij} \neq 0$。

!!! theorem "定理 70B.12 (精度矩阵与条件独立性)"
    设 $\mathbf{x} \sim N(\boldsymbol{\mu}, \Sigma)$，$\Theta = \Sigma^{-1}$。则：

    $$x_i \perp x_j \mid \mathbf{x}_{\setminus\{i,j\}} \iff \Theta_{ij} = 0$$

    其中 $\mathbf{x}_{\setminus\{i,j\}}$ 表示除 $x_i, x_j$ 外的所有分量。

??? proof "证明"
    在多元正态分布中，$x_i$ 和 $x_j$ 给定 $\mathbf{x}_{\setminus\{i,j\}}$ 后的条件分布仍为二元正态。条件相关系数为：

    $$\rho_{ij \mid \text{rest}} = -\frac{\Theta_{ij}}{\sqrt{\Theta_{ii}\Theta_{jj}}}$$

    这可以通过 Schur 补公式推导：$\Sigma$ 的 $(i,j)$ 子矩阵的 Schur 补与 $\Theta$ 的对应元素有直接关系。

    具体地，将 $\Sigma$ 的逆分块：$\Theta_{ij}$ 是 $x_i$ 对 $\mathbf{x}_{\setminus\{i\}}$ 回归后残差与 $x_j$ 对 $\mathbf{x}_{\setminus\{j\}}$ 回归后残差的（带符号的）条件协方差。

    因此 $\Theta_{ij} = 0 \iff \rho_{ij \mid \text{rest}} = 0 \iff$ 条件独立（多元正态下）。

    $\blacksquare$

!!! definition "定义 70B.19 (Graphical LASSO)"
    在高维情形（$n \gg m$）中，样本协方差矩阵 $\hat{\Sigma}$ 不可逆或估计不稳定。**Graphical LASSO** 通过 $L_1$ 惩罚估计稀疏的精度矩阵：

    $$\hat{\Theta} = \arg\min_{\Theta \succ 0} \left[\text{tr}(\hat{\Sigma}\Theta) - \log\det\Theta + \rho \sum_{i \neq j} |\Theta_{ij}|\right]$$

    其中 $\rho > 0$ 为正则化参数。$L_1$ 惩罚使得 $\hat{\Theta}$ 的许多非对角元素恰好为 0，从而给出稀疏网络。

    该优化问题是凸的（$-\log\det$ 是凸函数在正定锥上的限制），可以高效求解。

!!! example "例 70B.9"
    **酵母基因调控网络推断**。从 $m = 100$ 个微阵列实验中测量 $n = 500$ 个基因的表达水平。

    直接估计 $500 \times 500$ 的精度矩阵需要至少 500 个样本（否则 $\hat{\Sigma}$ 奇异）。Graphical LASSO 通过稀疏性假设克服了这个限制。

    选取适当的 $\rho$ 后，$\hat{\Theta}$ 中约 2% 的非对角元素非零，对应于约 2,500 条调控边。这些边可以与已知的转录因子-靶基因关系数据库进行比较验证。

---

## 习题

!!! exercise "习题 70B.1"
    设氨基酸突变概率矩阵 $M$ 的特征值为 $1, \lambda_2, \ldots, \lambda_{20}$（$|\lambda_k| < 1$ 对 $k \geq 2$）。证明当 $k \to \infty$ 时，PAM-$k$ 评分矩阵趋于常数矩阵 $S_{ij}^{(\infty)} = \log_2(\pi_j / q_j)$，其中 $\boldsymbol{\pi}$ 是 $M$ 的平稳分布。

!!! exercise "习题 70B.2"
    对 Jukes-Cantor 模型，证明：

    (a) 当 $t \to \infty$ 时，$P(t) \to \frac{1}{4}J$（等概率矩阵），即所有碱基频率趋于 1/4。

    (b) 进化距离 $d$ 与观察到的替换比例 $p$ 的关系 $d = -3/4 \ln(1 - 4p/3)$ 要求 $p < 3/4$。解释这个上限的含义。

!!! exercise "习题 70B.3"
    对 Kimura 双参数模型，利用速率矩阵 $Q_K$ 的特征值分解，导出转移概率矩阵 $P(t)$ 的显式表达式。验证当 $\alpha = \beta$ 时退化为 Jukes-Cantor 模型。

!!! exercise "习题 70B.4"
    **GTR 模型的参数计数**。GTR 模型有多少个自由参数？（提示：考虑 $\pi_A + \pi_C + \pi_G + \pi_T = 1$ 的约束以及速率矩阵的归一化。）

    列出以下模型的参数数量层级：JC69 $\subset$ K80 $\subset$ HKY85 $\subset$ GTR。

!!! exercise "习题 70B.5"
    一个简化的 HMM 有两个隐状态和观测符号 $\{A, B\}$。

    $$\mathbf{A} = \begin{pmatrix} 0.7 & 0.3 \\ 0.4 & 0.6 \end{pmatrix}, \quad \mathbf{B} = \begin{pmatrix} 0.9 & 0.1 \\ 0.2 & 0.8 \end{pmatrix}, \quad \boldsymbol{\pi} = \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}$$

    对观测序列 $O = (A, B, A)$，用前向算法计算 $P(O \mid \lambda)$。将计算过程写成矩阵乘法形式。

!!! exercise "习题 70B.6"
    证明化学计量矩阵的零空间中每个向量对应于一种可能的稳态通量分布。如果 $\mathbf{v}_1$ 和 $\mathbf{v}_2$ 都是可行通量（满足 $S\mathbf{v} = 0$ 且 $\mathbf{v} \geq 0$），$\alpha \mathbf{v}_1 + (1-\alpha)\mathbf{v}_2$（$\alpha \in [0,1]$）是否仍然可行？可行通量的集合是什么几何形状？

!!! exercise "习题 70B.7"
    对 Repressilator 模型（例 70B.7），

    (a) 当 $\alpha = 1$ 时，所有特征值实部为负，验证稳态是稳定的。

    (b) 精确计算使 Hopf 分岔发生的临界 $\alpha$ 值（即 $\text{Re}(\lambda) = 0$ 的条件）。

    (c) 讨论 $n$ 基因环状抑制网络（$n$ 为奇数）的类似分析。

!!! exercise "习题 70B.8"
    对 $3 \times 3$ 精度矩阵 $\Theta = \begin{pmatrix} 2 & -1 & 0 \\ -1 & 3 & -0.5 \\ 0 & -0.5 & 1.5 \end{pmatrix}$：

    (a) 画出对应的条件独立图（基因 1 和基因 3 条件独立吗？）

    (b) 计算协方差矩阵 $\Sigma = \Theta^{-1}$，并验证 $\Sigma_{13} \neq 0$（即边际相关但条件独立）

    (c) 计算所有条件相关系数 $\rho_{ij|\text{rest}}$

!!! exercise "习题 70B.9"
    **PCA 与 $F_{ST}$**。考虑两个等大小的亚群，$F_{ST} = 0.01$，$p = 100,000$ 个独立 SNP。

    (a) 估计遗传相关矩阵 $G$ 的最大特征值大约是多少？（提示：信号强度与 $p \cdot F_{ST}$ 有关。）

    (b) 当 $n$ 个个体中 $n/2$ 来自每个亚群时，PCA 的第一主成分应该能将两个亚群分开。需要多少个 SNP 才能保证统计显著性？

!!! exercise "习题 70B.10"
    在 Graphical LASSO 中，证明目标函数

    $$L(\Theta) = \text{tr}(\hat{\Sigma}\Theta) - \log\det\Theta + \rho \sum_{i \neq j} |\Theta_{ij}|$$

    在 $\Theta \succ 0$ 上是凸函数。（提示：$\text{tr}(\hat{\Sigma}\Theta)$ 是线性的，$-\log\det\Theta$ 在正定锥上是凸的，$|\Theta_{ij}|$ 是凸的。）
