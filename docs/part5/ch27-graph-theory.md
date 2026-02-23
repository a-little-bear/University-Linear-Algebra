# 第 27 章 线性代数在图论与网络中的应用

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 对称矩阵谱定理(Ch7-8) · 非负矩阵(Ch17)

**脉络**：邻接/关联/Laplacian 矩阵(图的线性代数表示) → 图谱(谱图论) → Laplacian 与连通性(Kirchhoff 矩阵树定理) → Cheeger 不等式(谱聚类) → 随机游走/PageRank(Markov 链) → 扩展图(Ramanujan 界) → 图着色(特征值界) → 网络流(LP 对偶)

</div>

图论与线性代数的交汇催生了谱图论（spectral graph theory），其核心思想是将图的组合结构编码为矩阵，再利用矩阵的谱性质揭示图的全局结构。本章从图的矩阵表示出发，依次展开图谱理论、Laplacian 矩阵与连通性、Cheeger 不等式与图分割、随机游走与 PageRank、扩展图、图着色、网络流与线性规划。

---

## 27.1 图的矩阵表示

<div class="context-flow" markdown>

**三种矩阵**：邻接矩阵 $A$（对称 ↔ 无向图） → 关联矩阵 $B$（顶点-边关系） → Laplacian $L = D - A$（最重要）

**链接**：Ch7 对称矩阵的谱性质直接应用于图的邻接矩阵和 Laplacian

</div>

图的组合信息可以完整地编码为矩阵，不同矩阵表示各有优势。

!!! definition "定义 27.1 (邻接矩阵)"
    设 $G = (V, E)$ 为有 $n$ 个顶点的简单无向图。$G$ 的**邻接矩阵**（adjacency matrix）$A \in \mathbb{R}^{n \times n}$ 定义为

    $$
    A_{ij} = \begin{cases} 1 & \text{若 } \{i, j\} \in E, \\ 0 & \text{否则}. \end{cases}
    $$

    $A$ 为对称矩阵（$A = A^T$），对角线为零（无自环）。对加权图，$A_{ij} = w_{ij} \ge 0$。

!!! definition "定义 27.2 (度矩阵与 Laplacian 矩阵)"
    **度矩阵**（degree matrix）$D = \operatorname{diag}(d_1, \ldots, d_n)$，$d_i = \sum_j A_{ij}$。**Laplacian 矩阵**（Laplacian matrix）定义为

    $$
    L = D - A.
    $$

    $L$ 对称半正定，且 $L\mathbf{1} = \mathbf{0}$（$\mathbf{1}$ 为全一向量），故 $0$ 总是 $L$ 的特征值。**规范化 Laplacian** 为

    $$
    \mathcal{L} = D^{-1/2} L D^{-1/2} = I - D^{-1/2} A D^{-1/2}.
    $$

!!! definition "定义 27.3 (关联矩阵)"
    给定有向图 $G$ 的一个定向，**关联矩阵**（incidence matrix）$B \in \mathbb{R}^{n \times m}$（$m = |E|$）定义为

    $$
    B_{ve} = \begin{cases} +1 & \text{若 } v \text{ 为边 } e \text{ 的尾}, \\ -1 & \text{若 } v \text{ 为边 } e \text{ 的头}, \\ 0 & \text{否则}. \end{cases}
    $$

    关键性质：$L = BB^T$（与定向无关）。

!!! theorem "定理 27.1 (Laplacian 的二次形式)"
    对任意 $\mathbf{x} \in \mathbb{R}^n$，

    $$
    \mathbf{x}^T L \mathbf{x} = \sum_{\{i,j\} \in E} w_{ij}(x_i - x_j)^2.
    $$

    这直接证明 $L \succeq 0$（半正定），且 $\mathbf{x}^T L \mathbf{x} = 0$ 当且仅当 $\mathbf{x}$ 在每个连通分量上为常数。

??? proof "证明"
    由 $L = D - A$，

    $$
    \mathbf{x}^T L \mathbf{x} = \mathbf{x}^T D \mathbf{x} - \mathbf{x}^T A \mathbf{x} = \sum_i d_i x_i^2 - \sum_{\{i,j\} \in E} 2w_{ij} x_i x_j.
    $$

    由 $d_i = \sum_j w_{ij}$，

    $$
    \sum_i d_i x_i^2 = \sum_{\{i,j\} \in E} w_{ij}(x_i^2 + x_j^2).
    $$

    因此 $\mathbf{x}^T L \mathbf{x} = \sum_{\{i,j\} \in E} w_{ij}(x_i^2 - 2x_i x_j + x_j^2) = \sum_{\{i,j\} \in E} w_{ij}(x_i - x_j)^2 \ge 0$。

    等号成立当且仅当对每条边 $\{i, j\}$ 有 $x_i = x_j$，即 $\mathbf{x}$ 在每个连通分量上为常数。$\blacksquare$

!!! example "例 27.1"
    **路径图和环图的矩阵表示。** 路径图 $P_4$（4 个顶点的路径 $1 - 2 - 3 - 4$）：

    $$
    A = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 1 & 0 & 1 & 0 \\ 0 & 1 & 0 & 1 \\ 0 & 0 & 1 & 0 \end{pmatrix}, \quad L = \begin{pmatrix} 1 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 1 \end{pmatrix}.
    $$

    $L$ 的特征值为 $0, 2 - \sqrt{2}, 2, 2 + \sqrt{2}$。零特征值的重数为 $1$，确认 $P_4$ 是连通图。

---

## 27.2 图的谱

<div class="context-flow" markdown>

**核心概念**：图的谱 = 邻接矩阵（或 Laplacian）的特征值集合 → 谱携带丰富的图结构信息（正则性、二部性、直径等）

**链接**：Ch7 对称矩阵特征值的极值性质（Courant-Fischer）在谱图论中反复出现

</div>

图的谱（spectrum）是其矩阵表示的特征值集合，蕴含了图的丰富结构信息。

!!! definition "定义 27.4 (图的谱)"
    设 $G$ 为 $n$ 阶简单图。$G$ 的**谱**（spectrum）为邻接矩阵 $A$ 的特征值（含重数）：

    $$
    \lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_n.
    $$

    **Laplacian 谱**为 $L$ 的特征值：$0 = \mu_1 \le \mu_2 \le \cdots \le \mu_n$。

!!! theorem "定理 27.2 (谱的基本性质)"
    设 $G$ 为 $n$ 阶简单图，邻接谱为 $\lambda_1 \ge \cdots \ge \lambda_n$：

    1. $\sum_{i=1}^n \lambda_i = \operatorname{tr}(A) = 0$。
    2. $\sum_{i=1}^n \lambda_i^2 = \operatorname{tr}(A^2) = 2|E|$。
    3. $\sum_{i=1}^n \lambda_i^3 = \operatorname{tr}(A^3) = 6 \times (\text{三角形数量})$。
    4. $\lambda_1 \ge \bar{d}$（平均度），等号在正则图时成立。
    5. $G$ 是二部图 $\iff$ 谱关于 $0$ 对称（$\lambda_i = -\lambda_{n+1-i}$）。

??? proof "证明"
    (1)-(3)：$\operatorname{tr}(A^k) = \sum_i \lambda_i^k$，而 $(A^k)_{ii}$ 计数从 $i$ 出发长度为 $k$ 的闭途径数。$k = 1$：无自环故 $(A)_{ii} = 0$。$k = 2$：$(A^2)_{ii} = d_i$，求和得 $2|E|$。$k = 3$：$(A^3)_{ii}$ 计数经过 $i$ 的三角形数的 $2$ 倍，总和为 $6 \times$ 三角形数。

    (4)：$\lambda_1 = \max_{\|\mathbf{x}\| = 1} \mathbf{x}^T A \mathbf{x} \ge \frac{1}{n}\mathbf{1}^T A \mathbf{1} = \frac{1}{n}\sum_i d_i = \bar{d}$。

    (5)：若 $G$ 为二部图（$V = V_1 \cup V_2$），令 $S$ 为对角矩阵，$S_{ii} = +1$（$i \in V_1$），$S_{ii} = -1$（$i \in V_2$）。则 $SAS = -A$，故 $\lambda$ 是特征值蕴含 $-\lambda$ 也是。$\blacksquare$

!!! example "例 27.2"
    **完全图与 Petersen 图的谱。** 完全图 $K_n$：$A = J - I$（$J$ 为全一矩阵）。$J$ 的特征值为 $n$（重数 $1$）和 $0$（重数 $n-1$），故 $A$ 的谱为 $n - 1$（重数 $1$）和 $-1$（重数 $n-1$）。

    Petersen 图（$10$ 顶点 $3$-正则图）的邻接谱为 $3, 1, 1, 1, 1, 1, -2, -2, -2, -2$。$\lambda_1 = 3 = d$（正则），$\lambda_2 = 1$ 较小，意味着 Petersen 图具有良好的扩展性。

---

## 27.3 Laplacian 矩阵与连通性

<div class="context-flow" markdown>

**核心结论**：$\mu_2 > 0$ ⇔ 图连通 → $\mu_2$（Fiedler 值/代数连通度）度量连通的"强度" → Kirchhoff 定理：生成树数 = $\frac{1}{n}\lambda_1(L)\cdots\lambda_{n-1}(L)$

**链接**：Ch3 行列式理论在 Kirchhoff 定理中的核心应用

</div>

Laplacian 矩阵的谱完整地刻画了图的连通性。

!!! definition "定义 27.5 (代数连通度)"
    $L$ 的第二小特征值 $\mu_2$ 称为图的**代数连通度**（algebraic connectivity）或 **Fiedler 值**。对应的特征向量称为 **Fiedler 向量**。

!!! theorem "定理 27.3 (Laplacian 与连通性)"
    设 $L$ 的特征值为 $0 = \mu_1 \le \mu_2 \le \cdots \le \mu_n$：

    1. $\mu_1 = 0$ 的重数等于图的连通分量数 $k$。
    2. $G$ 连通 $\iff$ $\mu_2 > 0$。
    3. 对连通图，$\mu_2 \le \frac{n}{n-1}\min_i d_i$。

??? proof "证明"
    $L\mathbf{x} = \mathbf{0}$ 的解空间由各连通分量的指示向量张成。若 $G$ 有 $k$ 个连通分量，定义 $\mathbf{x}^{(j)}$ 为第 $j$ 个分量的指示向量，则 $L\mathbf{x}^{(j)} = \mathbf{0}$（每条边的两个端点在同一分量中，故 $(x_i - x_j)^2 = 0$）。这 $k$ 个向量线性无关，因此 $\ker(L)$ 的维数为 $k$，即零特征值重数为 $k$。

    特别地，$k = 1$（连通）当且仅当 $\mu_2 > 0$。$\blacksquare$

!!! theorem "定理 27.4 (Kirchhoff 矩阵树定理)"
    设 $G$ 为连通图，$L$ 为其 Laplacian。$G$ 的生成树数量 $\tau(G)$ 为

    $$
    \tau(G) = \frac{1}{n}\mu_2 \mu_3 \cdots \mu_n = \frac{1}{n}\prod_{i=2}^{n} \mu_i.
    $$

    等价地，$\tau(G)$ 等于 $L$ 的任意 $(n-1) \times (n-1)$ 余子式。

??? proof "证明"
    定义 $L$ 删去第 $i$ 行和第 $i$ 列后的矩阵为 $L_i$。Kirchhoff 定理断言 $\tau(G) = \det(L_i)$（对任意 $i$）。

    由 $L$ 的特征值，$\det(L_i)$ 可以通过矩阵-树定理的组合证明得到。另一方面，注意到

    $$
    \det(\lambda I - L) = \lambda \prod_{i=2}^n (\lambda - \mu_i),
    $$

    而 $\det(\lambda I - L)$ 的 $\lambda$ 项系数为 $(-1)^{n-1} \sum_i \det(L_i) / n \cdot n$。利用特征多项式展开和 Cauchy-Binet 公式可以证明 $\det(L_i) = \frac{1}{n}\prod_{i=2}^n \mu_i$。$\blacksquare$

!!! example "例 27.3"
    **完全图的生成树计数。** $K_n$ 的 Laplacian $L = nI - J$，特征值为 $0$（重数 $1$）和 $n$（重数 $n-1$）。由 Kirchhoff 定理，

    $$
    \tau(K_n) = \frac{1}{n} \cdot n^{n-1} = n^{n-2}.
    $$

    这就是著名的 **Cayley 公式**。例如 $\tau(K_4) = 4^2 = 16$。

---

## 27.4 Cheeger 不等式与图分割

<div class="context-flow" markdown>

**核心不等式**：$\frac{h^2}{2d_{\max}} \le \mu_2 \le 2h$ → Fiedler 值 $\mu_2$ 与 Cheeger 常数 $h$ 等价 → 谱聚类算法：用 Fiedler 向量将顶点分为两组

**应用**：图像分割、社区发现、数据聚类

</div>

Cheeger 不等式将代数量（$\mu_2$）与组合量（图的最优割）联系起来，是谱聚类的理论基础。

!!! definition "定义 27.6 (等周常数/Cheeger 常数)"
    图 $G$ 的 **Cheeger 常数**（或等周常数/conductance）定义为

    $$
    h(G) = \min_{S \subset V,\, 0 < |S| \le n/2} \frac{|\partial(S)|}{\min(|S|, |V \setminus S|)},
    $$

    其中 $\partial(S) = \{\{u, v\} \in E : u \in S, v \notin S\}$ 为 $S$ 的边界。$h(G)$ 度量将图切割为两部分的最小代价。

!!! theorem "定理 27.5 (离散 Cheeger 不等式)"
    对 $d$-正则图 $G$，设 $\mu_2$ 为规范化 Laplacian $\mathcal{L}$ 的第二小特征值，$h$ 为 Cheeger 常数，则

    $$
    \frac{h^2}{2} \le \mu_2 \le 2h.
    $$

    对一般图（使用 $L$ 和适当定义的 $h$），类似不等式成立。

??? proof "证明"
    **上界** $\mu_2 \le 2h$：取 $S$ 为实现 $h$ 的集合，构造测试向量 $\mathbf{x}$ 使得 $x_i = |V \setminus S|$（$i \in S$），$x_i = -|S|$（$i \notin S$），$\mathbf{x} \perp \mathbf{1}$。由 Rayleigh 商

    $$
    \mu_2 \le \frac{\mathbf{x}^T L \mathbf{x}}{\mathbf{x}^T D \mathbf{x}} = \frac{n^2 |\partial(S)|}{d \cdot |S| \cdot |V \setminus S|} \le \frac{2|\partial(S)|}{d \cdot |S|} = \frac{2h}{d} \cdot d = 2h.
    $$

    **下界** $h^2/2 \le \mu_2$：设 $\mathbf{f}$ 为 $\mu_2$ 对应的特征向量。对 $\mathbf{f}$ 的分量排序后，利用水平集分析和 Cauchy-Schwarz 不等式建立下界。$\blacksquare$

!!! example "例 27.4"
    **谱聚类算法。** 给定图 $G$（或数据点的相似度图），谱聚类的步骤为：

    1. 计算 Laplacian $L = D - A$（或规范化 Laplacian $\mathcal{L}$）。
    2. 求 $\mathcal{L}$ 的最小 $k$ 个特征值对应的特征向量 $\mathbf{u}_1, \ldots, \mathbf{u}_k$。
    3. 构造矩阵 $U \in \mathbb{R}^{n \times k}$，每行归一化。
    4. 对 $U$ 的行向量进行 $k$-means 聚类。

    Fiedler 向量（$\mu_2$ 对应的特征向量）的分量正负自然将顶点分为两组。对社交网络数据，谱聚类能有效发现社区结构。

---

## 27.5 随机游走与 PageRank

<div class="context-flow" markdown>

**链条**：图上随机游走 → 转移矩阵 $P = D^{-1}A$ → Perron-Frobenius 定理(Ch17) → 平稳分布 → Google PageRank = 阻尼随机游走的平稳分布

**链接**：Ch17 非负矩阵理论直接应用

</div>

图上的随机游走将图论问题转化为随机矩阵分析，PageRank 是其最著名的应用。

!!! definition "定义 27.7 (图上的随机游走)"
    在无向图 $G$ 上，**简单随机游走**（simple random walk）在每一步从当前顶点 $i$ 等概率移动到其邻居 $j$。转移概率矩阵为

    $$
    P = D^{-1}A, \quad P_{ij} = \frac{A_{ij}}{d_i}.
    $$

    $P$ 为行随机矩阵（每行和为 $1$）。

!!! definition "定义 27.8 (PageRank)"
    对有向图 $G$，**PageRank** 向量 $\boldsymbol{\pi}$ 是修改后的随机游走的平稳分布。Google 矩阵定义为

    $$
    M = \alpha P + (1 - \alpha) \frac{1}{n} \mathbf{1}\mathbf{1}^T,
    $$

    其中 $\alpha \in (0, 1)$ 为阻尼因子（通常 $\alpha = 0.85$），$P$ 为列随机的链接矩阵。$\boldsymbol{\pi}$ 满足 $M\boldsymbol{\pi} = \boldsymbol{\pi}$，$\boldsymbol{\pi}^T \mathbf{1} = 1$。

!!! theorem "定理 27.6 (PageRank 的存在唯一性)"
    对 $0 < \alpha < 1$，Google 矩阵 $M$ 是严格正矩阵（所有元素 $> 0$），因此：

    1. $M$ 有唯一的特征值 $1$（Perron-Frobenius 定理），对应的特征向量 $\boldsymbol{\pi} > \mathbf{0}$（所有分量正）。
    2. 幂迭代 $\boldsymbol{\pi}^{(k+1)} = M\boldsymbol{\pi}^{(k)}$ 以速率 $\alpha$ 几何收敛到 $\boldsymbol{\pi}$。
    3. $M$ 的第二大特征值模 $|\lambda_2| \le \alpha < 1$。

??? proof "证明"
    $M$ 的每个元素至少为 $(1-\alpha)/n > 0$，故 $M$ 为严格正矩阵。由 Perron-Frobenius 定理（Ch17），$M$ 有唯一的最大特征值 $1$（因 $M$ 为列随机或行随机），对应正特征向量。

    收敛速率由谱间隙 $1 - |\lambda_2|$ 控制。由于 $M = \alpha P + (1-\alpha)J/n$，其中 $P$ 的谱半径为 $1$，$J/n$ 的非零特征值为 $1$（重数 $1$），可以证明 $|\lambda_2(M)| \le \alpha$。$\blacksquare$

!!! example "例 27.5"
    **简单网络的 PageRank 计算。** 考虑三页面网络：页面 1 链接到 2 和 3，页面 2 链接到 3，页面 3 链接到 1。列随机矩阵为

    $$
    P = \begin{pmatrix} 0 & 0 & 1 \\ 1/2 & 0 & 0 \\ 1/2 & 1 & 0 \end{pmatrix}.
    $$

    取 $\alpha = 0.85$，$M = 0.85P + 0.05 \cdot \mathbf{1}\mathbf{1}^T$。幂迭代从 $\boldsymbol{\pi}^{(0)} = (1/3, 1/3, 1/3)^T$ 出发，收敛到 $\boldsymbol{\pi} \approx (0.387, 0.214, 0.399)^T$。页面 3 的 PageRank 最高（被页面 1 和 2 链接），页面 1 次之（被页面 3 链接，而 3 的权重高）。

---

## 27.6 扩展图

<div class="context-flow" markdown>

**定义**：$d$-正则图的 $\lambda_2(A) \le d - \varepsilon$ → 好的扩展图 = 谱间隙大 → Alon-Boppana 界：$\lambda_2 \ge 2\sqrt{d-1} - o(1)$ → Ramanujan 图达到此界

**应用**：纠错码、去随机化、网络设计

</div>

扩展图（expander graphs）是具有良好连通性和伪随机性质的稀疏图，其定义和分析本质上依赖于谱理论。

!!! definition "定义 27.9 (谱扩展图)"
    $d$-正则图 $G$ 称为 $(n, d, \lambda)$**-图**，若 $A$ 的第二大特征值（绝对值）$\lambda = \max(|\lambda_2|, |\lambda_n|) \le \lambda$。$\lambda$ 越小，扩展性越好。**谱间隙**（spectral gap）定义为 $d - \lambda$。

!!! theorem "定理 27.7 (Expander Mixing Lemma)"
    对 $(n, d, \lambda)$-图 $G$，任意两个顶点子集 $S, T \subseteq V$，边数满足

    $$
    \left| e(S, T) - \frac{d \cdot |S| \cdot |T|}{n} \right| \le \lambda \sqrt{|S| \cdot |T|},
    $$

    其中 $e(S, T)$ 为 $S$ 到 $T$ 的边数。

??? proof "证明"
    设 $\mathbf{x} = \mathbf{1}_S$，$\mathbf{y} = \mathbf{1}_T$。则 $e(S, T) = \mathbf{x}^T A \mathbf{y}$。将 $\mathbf{x}$ 和 $\mathbf{y}$ 分解为 $A$ 的特征向量：$\mathbf{x} = \sum \hat{x}_i \mathbf{v}_i$，$\mathbf{y} = \sum \hat{y}_i \mathbf{v}_i$。

    $$
    \mathbf{x}^T A \mathbf{y} = \sum_i \lambda_i \hat{x}_i \hat{y}_i = d \hat{x}_1 \hat{y}_1 + \sum_{i \ge 2} \lambda_i \hat{x}_i \hat{y}_i.
    $$

    $\hat{x}_1 = \langle \mathbf{x}, \frac{\mathbf{1}}{\sqrt{n}} \rangle = |S|/\sqrt{n}$，类似 $\hat{y}_1 = |T|/\sqrt{n}$。第一项为 $d|S||T|/n$。第二项由 Cauchy-Schwarz 控制：

    $$
    \left|\sum_{i \ge 2} \lambda_i \hat{x}_i \hat{y}_i\right| \le \lambda \sqrt{\sum_{i \ge 2} \hat{x}_i^2} \sqrt{\sum_{i \ge 2} \hat{y}_i^2} \le \lambda \|\mathbf{x}\| \|\mathbf{y}\| = \lambda\sqrt{|S| \cdot |T|}.
    $$

    $\blacksquare$

!!! theorem "定理 27.8 (Alon-Boppana 界)"
    对无穷族 $d$-正则图 $\{G_n\}$（$|V(G_n)| \to \infty$），

    $$
    \liminf_{n \to \infty} \lambda_2(G_n) \ge 2\sqrt{d - 1}.
    $$

    达到此界的图（$\lambda \le 2\sqrt{d-1}$）称为 **Ramanujan 图**。

??? proof "证明"
    （证明梗概）利用 $d$-正则无穷树 $T_d$ 的谱。$T_d$ 的邻接算子的谱为 $[-2\sqrt{d-1}, 2\sqrt{d-1}]$。有限 $d$-正则图可以"局部"近似 $T_d$（高 girth 时），因此其非平凡特征值不能全部远离此区间。精确证明使用 trace 方法和组合计数。$\blacksquare$

!!! example "例 27.6"
    **Ramanujan 图的构造。** Lubotzky-Phillips-Sarnak (LPS) 图是 Ramanujan 图的经典构造。对素数 $p \equiv 1 \pmod{4}$ 和 $q$（$q \ne p$），LPS 图是 $(p+1)$-正则的 Cayley 图，顶点集为 $\text{PSL}(2, \mathbb{F}_q)$。

    $\lambda_2 \le 2\sqrt{p}$，满足 Ramanujan 界。这些图有 $O(q^3)$ 个顶点，$(p+1)$-正则，且谱间隙为 $p + 1 - 2\sqrt{p} = (\sqrt{p} - 1)^2$，几乎是最优的稀疏扩展图。

---

## 27.7 图着色与特征值

<div class="context-flow" markdown>

**界**：$\chi(G) \ge 1 + \lambda_1 / |\lambda_n|$（Hoffman 界）→ 特征值给出色数下界 → 独立数上界：$\alpha(G) \le n \cdot |\lambda_n| / (\lambda_1 + |\lambda_n|)$

**应用**：组合优化中的松弛界

</div>

图的特征值可以为色数和独立数提供有效的界。

!!! definition "定义 27.10 (色数)"
    图 $G$ 的**色数**（chromatic number）$\chi(G)$ 是使得 $G$ 可正常着色（相邻顶点不同色）的最少颜色数。

!!! theorem "定理 27.9 (Hoffman 色数界)"
    对非空图 $G$，

    $$
    \chi(G) \ge 1 + \frac{\lambda_1}{|\lambda_n|} = 1 - \frac{\lambda_1}{\lambda_n},
    $$

    其中 $\lambda_1$ 和 $\lambda_n$ 分别为 $A$ 的最大和最小特征值。

??? proof "证明"
    设 $G$ 可正常 $k$-着色，颜色类为 $C_1, \ldots, C_k$（独立集）。对每个颜色类 $C_j$，其诱导子图无边，因此

    $$
    \sum_{i, j \in C_s} A_{ij} = 0, \quad \forall s.
    $$

    令 $\mathbf{x}^{(s)}$ 为 $C_s$ 的（适当中心化的）指示向量。利用 $A$ 的二次形式和 Rayleigh 商，可以证明 $\lambda_n \le -\lambda_1 / (k-1)$，即 $k \ge 1 + \lambda_1/|\lambda_n|$。$\blacksquare$

!!! theorem "定理 27.10 (Hoffman 独立数界)"
    对 $d$-正则图 $G$，最大独立集大小满足

    $$
    \alpha(G) \le \frac{n \cdot |\lambda_n|}{\lambda_1 + |\lambda_n|} = \frac{n \cdot |\lambda_n|}{d + |\lambda_n|}.
    $$

??? proof "证明"
    设 $S$ 为独立集，$|S| = \alpha$。令 $\mathbf{x} = \mathbf{1}_S - \frac{\alpha}{n}\mathbf{1}$（$\mathbf{x} \perp \mathbf{1}$）。由 $S$ 为独立集，$\mathbf{1}_S^T A \mathbf{1}_S = 0$。

    $$
    \mathbf{x}^T A \mathbf{x} = -\frac{2\alpha}{n}\mathbf{1}_S^T A \mathbf{1} + \frac{\alpha^2}{n^2}\mathbf{1}^T A \mathbf{1} = -\frac{2\alpha d\alpha}{n} + \frac{\alpha^2 dn}{n^2} = -\frac{\alpha^2 d}{n}.
    $$

    由 $\mathbf{x} \perp \mathbf{1}$，$\mathbf{x}^T A \mathbf{x} \ge \lambda_n \|\mathbf{x}\|^2 = \lambda_n (\alpha - \alpha^2/n)$。因此

    $$
    -\frac{\alpha^2 d}{n} \ge \lambda_n \alpha (1 - \alpha/n),
    $$

    化简得 $\alpha \le n|\lambda_n| / (d + |\lambda_n|)$。$\blacksquare$

!!! example "例 27.7"
    **Kneser 图的色数。** Kneser 图 $K(n, k)$ 的顶点为 $\{1, \ldots, n\}$ 的所有 $k$-元子集，两个顶点相邻当且仅当对应子集不相交。

    $K(5, 2)$ 即 Petersen 图，谱为 $3, 1^5, (-2)^4$。Hoffman 界给出 $\chi \ge 1 + 3/2 = 2.5$，即 $\chi \ge 3$。实际上 $\chi(K(5, 2)) = 3$（Lovasz 证明的 Kneser 猜想），谱界在此恰好紧密。

---

## 27.8 网络流与线性规划

<div class="context-flow" markdown>

**模型**：网络流 = 图上的线性约束优化 → 最大流/最小割 = LP 对偶的特例 → 关联矩阵 $B$ 在流守恒约束中的核心角色

**链接**：Ch25 线性规划理论在图上的具体化

</div>

网络流问题是图论与线性规划的经典交汇点。

!!! definition "定义 27.11 (网络流)"
    **网络流问题**：给定有向图 $G = (V, E)$，源 $s$，汇 $t$，容量函数 $c : E \to \mathbb{R}_{\ge 0}$。**流** $\mathbf{f} \in \mathbb{R}^m$ 满足：

    1. **容量约束**：$0 \le f_e \le c_e$，$\forall e \in E$。
    2. **流守恒**：$\sum_{e \in \delta^+(v)} f_e = \sum_{e \in \delta^-(v)} f_e$，$\forall v \neq s, t$。

    用关联矩阵 $B$，流守恒写为 $B\mathbf{f} = \mathbf{b}$（$b_s = -F$，$b_t = F$，其余为 $0$）。

!!! theorem "定理 27.11 (最大流-最小割定理)"
    最大流值等于最小割容量：

    $$
    \max_{\mathbf{f}} F = \min_{S:\, s \in S,\, t \notin S} \sum_{e \in \delta^+(S)} c_e.
    $$

    这是线性规划对偶定理在网络上的体现：最大流 LP 的对偶恰好是最小割 LP。

??? proof "证明"
    **弱对偶**（$\le$）：任何流 $\mathbf{f}$ 和割 $(S, \bar{S})$，$F = \sum_{e \in \delta^+(S)} f_e - \sum_{e \in \delta^-(S)} f_e \le \sum_{e \in \delta^+(S)} c_e$。

    **强对偶**（$=$）：由线性规划强对偶定理，或由 Ford-Fulkerson 增广路算法的终止性证明。当不存在从 $s$ 到 $t$ 的增广路时，令 $S$ 为从 $s$ 在残余图中可达的顶点集，则 $(S, \bar{S})$ 是最小割，且当前流为最大流。$\blacksquare$

!!! theorem "定理 27.12 (全幺模性与整数流)"
    有向图的关联矩阵 $B$ 是**全幺模**（totally unimodular）的：$B$ 的每个方阵子矩阵的行列式为 $0$、$+1$ 或 $-1$。因此，当容量 $c_e$ 为整数时，最大流 LP 的最优解自动为整数。

??? proof "证明"
    对 $B$ 的 $k \times k$ 子矩阵 $B'$ 进行归纳。$k = 1$：$B'$ 的元素为 $0, \pm 1$。对 $k > 1$：若 $B'$ 某列全零，行列式为 $0$。若某列有恰好一个非零元素，按该列展开，归结为 $(k-1)$ 阶子矩阵。若某列有两个非零元素（$+1$ 和 $-1$），则该列求和为零，行和线性相关可用于化简，再归纳。

    全幺模性保证 $B\mathbf{f} = \mathbf{b}$ 在整数 $\mathbf{b}$ 下有整数基本可行解，因此 LP 最优解为整数。$\blacksquare$

!!! example "例 27.8"
    **最大流的线性规划表述。** 网络有 $4$ 个顶点 $\{s, a, b, t\}$，边和容量为 $s \to a$ (3), $s \to b$ (2), $a \to b$ (1), $a \to t$ (2), $b \to t$ (3)。LP 表述：

    $$
    \max F, \quad \text{s.t.} \quad B\mathbf{f} = F(\mathbf{e}_t - \mathbf{e}_s), \quad \mathbf{0} \le \mathbf{f} \le \mathbf{c}.
    $$

    最大流 $F = 4$（$s \to a$: 3, $s \to b$: 1, $a \to b$: 1, $a \to t$: 2, $b \to t$: 2）。最小割为 $\{s, a\}$，割容量 $= 1 + 2 + 1 = 4$。对偶间隙为零，印证强对偶定理。关联矩阵的全幺模性保证了整数流的存在。

## 练习题

1. **[邻接矩阵] 对于一个无向图，如果它的邻接矩阵 $A$ 的平方 $A^2$ 的主对角线元素全为 0，这说明该图有什么特征？**

   ??? success "参考答案"
       这意味着图没有任何边。因为 $(A^2)_{ii}$ 表示顶点 $i$ 到自身的长度为 2 的路径数，即顶点 $i$ 的度数。度数全为 0 说明是孤立点图。

2. **[邻接矩阵] 如何通过邻接矩阵 $A$ 判断一个图是否包含三角形？**

   ??? success "参考答案"
       计算 $\operatorname{tr}(A^3)$。如果该迹大于 0，说明存在三角形（且三角形的总数等于 $\operatorname{tr}(A^3) / 6$）。

3. **[关联矩阵] 无向图的关联矩阵 $B$（顶点-边矩阵）的列向量之和有什么特点？**

   ??? success "参考答案"
       每一列恰好有两个 1，其余全为 0（因为每条边连接两个顶点）。因此其行向量之和为 $(2, 2, \ldots, 2)^T$。

4. **[拉普拉斯矩阵] 证明任何无向图的拉普拉斯矩阵 $L$ 总是半正定的。**

   ??? success "参考答案"
       对任意向量 $\mathbf{x}$，二次型 $\mathbf{x}^T L \mathbf{x} = \sum_{(i,j) \in E} (x_i - x_j)^2$ 是一组平方和，必然大于等于 0。

5. **[拉普拉斯矩阵] 为什么全 1 向量 $\mathbf{1}$ 总是拉普拉斯矩阵 $L$ 的特征向量？对应的特征值是多少？**

   ??? success "参考答案"
       因为 $L$ 的每一行的和都等于顶点的度数减去该顶点在非对角线上的 1 的个数（即度数），所以行和总为 0。这意味着 $L \mathbf{1} = \mathbf{0}$，对应的特征值是 0。

6. **[图的连通性] 图的连通分支数 $k$ 与拉普拉斯矩阵 $L$ 的哪个代数特征直接相关？**

   ??? success "参考答案"
       连通分支数 $k$ 恰好等于 $L$ 的特征值 0 的代数重数（也是几何重数），即 $L$ 的零空间的维数。

7. **[代数连通度] Fiedler 值（第二小特征值 $\lambda_2$）为 0 意味着什么？**

   ??? success "参考答案"
       意味着图不是连通的（它至少有两个独立的连通分支，从而至少有两个特征值为 0）。

8. **[二分图] 如果一个图是二分图（Bipartite graph），它的邻接矩阵谱有什么极强的对称性？**

   ??? success "参考答案"
       它的谱关于原点完全对称。即如果 $\lambda$ 是一个特征值，那么 $-\lambda$ 也必然是特征值，并且具有相同的重数。

9. **[全幺模矩阵] 在网络流问题中，为什么使用线性规划求解最大流时，即便不强制加上整数约束，求出的最优解（如果容量为整数）也自然全是整数？**

   ??? success "参考答案"
       因为有向图的关联矩阵是全幺模矩阵（Totally Unimodular Matrix），这意味着线性规划单纯形表中的所有基矩阵的行列式都是 $\pm 1$。由克莱姆法则可知，常数项为整数时，求出的解必定为整数。

10. **[爱因斯坦思考题] 为什么我们在描述像互联网（PageRank）、社交网络或者神经网络这样极度复杂的系统时，依然可以仅仅使用图的拉普拉斯矩阵的几个特征向量就能完成聚类（Spectral Clustering）？这揭示了复杂网络的什么本质？**

   ??? success "参考答案"
        这揭示了复杂网络在宏观层面的“低频主导”本质。就像琴弦振动一样，拉普拉斯矩阵的特征向量就是网络上的“驻波”。最小的非零特征值（Fiedler值）对应的特征向量，描述了网络在拉扯中最容易断裂的“低频瓶颈”。自然界和人类社会形成的复杂图结构，其信息传递和结构切分往往受制于极少数的宏观能量极小模式，这正是代数图论赋予我们的“宏观透视眼”。

## 本章小结

本章将线性代数作为解码工具引入了图论领域，主要内容包括：

1. **图的矩阵表示**：系统介绍了邻接矩阵、关联矩阵及其几何意义，将图上的路径计数和度序列等问题转化为纯粹的矩阵运算。
2. **拉普拉斯矩阵**：作为图论最重要的矩阵，推导了 $L = D - A = B B^T$ 的性质，确立了其二次型的能量物理意义（Dirichlet 能量）。
3. **连通性与 Fiedler 向量**：证明了特征值 $0$ 的重数等于连通分支数，并详述了代数连通度（$\lambda_2$）在谱聚类和图分割中的决定性作用。
4. **二分图与循环图的谱**：揭示了图的宏观拓扑结构（如二部性）如何严格限制了邻接矩阵谱的对称性。
5. **极值界与独立集**：利用特征值给出了关于图的最大独立集大小和色数的 Hoffmann 界等经典结论。
6. **网络流与全幺模性**：将图上的物质传输问题抽象为线性规划，并通过关联矩阵的全幺模性证明了离散优化与连续优化在图论中的完美重合。
