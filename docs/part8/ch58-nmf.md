# 第 58 章 非负矩阵分解

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · SVD(Ch11) · 非负矩阵(Ch17) · 优化(Ch25)

**本章脉络**：NMF 问题定义 → 与 SVD/PCA 的区别 → 乘性更新算法（Lee-Seung） → 交替最小二乘 → 唯一性条件 → 非负秩 → 稀疏/正则化 NMF → 应用

**延伸**：NMF 在文本挖掘（主题模型的矩阵视角）、音乐信息检索（音源分离）、高光谱成像（端元提取）和生物信息学（基因表达谱分析）中被广泛使用

</div>

非负矩阵分解（Nonnegative Matrix Factorization, NMF）是线性代数与优化理论交汇处的一个重要课题。给定一个非负矩阵 $V \in \mathbb{R}_{\geq 0}^{m \times n}$，NMF 旨在找到两个非负矩阵 $W \in \mathbb{R}_{\geq 0}^{m \times r}$、$H \in \mathbb{R}_{\geq 0}^{r \times n}$ 使得 $V \approx WH$。这一看似简单的约束条件——所有矩阵元素非负——产生了与 SVD 截然不同的分解特性：NMF 自然地给出"部分-整体"的表示，使得数据的每个分量都可以解释为基元素的加法叠加，而不涉及减法抵消。

Lee 和 Seung 在 1999 年 Nature 论文中将 NMF 引入机器学习社区，展示了它在人脸识别中学到的"部件"表示（眼睛、鼻子、嘴巴）远比 PCA 的"特征脸"（全脸的叠加与相消）更具可解释性。此后 NMF 成为数据科学中的标准工具之一。

---

## 58.1 NMF 问题定义

<div class="context-flow" markdown>

**核心问题**：如何形式化地定义非负矩阵分解问题？有哪些常用的目标函数？

</div>

!!! definition "定义 58.1 (非负矩阵分解)"
    给定非负矩阵 $V \in \mathbb{R}_{\geq 0}^{m \times n}$ 和正整数 $r$（$r \leq \min(m,n)$），**非负矩阵分解**（NMF）问题是求解

    $$\min_{W \geq 0,\, H \geq 0} \, D(V,\, WH),$$

    其中 $W \in \mathbb{R}_{\geq 0}^{m \times r}$，$H \in \mathbb{R}_{\geq 0}^{r \times n}$，$D$ 是某种距离或散度度量。

!!! definition "定义 58.2 (Frobenius 范数目标)"
    最常用的目标函数是 Frobenius 范数的平方：

    $$D_F(V, WH) = \|V - WH\|_F^2 = \sum_{i,j} (V_{ij} - (WH)_{ij})^2.$$

    此时 NMF 问题等价于

    $$\min_{W \geq 0,\, H \geq 0}\, \|V - WH\|_F^2.$$

!!! definition "定义 58.3 (KL 散度目标)"
    另一种常用目标是广义 Kullback-Leibler 散度：

    $$D_{KL}(V \| WH) = \sum_{i,j}\Bigl(V_{ij}\ln\frac{V_{ij}}{(WH)_{ij}} - V_{ij} + (WH)_{ij}\Bigr).$$

    当 $V$ 的列是（非归一化的）概率分布时，KL 散度目标更为自然。

!!! definition "定义 58.4 (Itakura-Saito 散度)"
    在音频处理中常用的 IS 散度：

    $$D_{IS}(V \| WH) = \sum_{i,j}\Bigl(\frac{V_{ij}}{(WH)_{ij}} - \ln\frac{V_{ij}}{(WH)_{ij}} - 1\Bigr).$$

    这三种目标函数都是 $\beta$-散度 $D_\beta$ 在 $\beta = 2, 1, 0$ 时的特殊情形。

!!! theorem "定理 58.1 (NMF 的 NP 困难性)"
    判定给定非负矩阵 $V$ 是否存在精确的秩-$r$ 非负分解 $V = WH$（$W \geq 0$，$H \geq 0$）是 **NP 困难**的（Vavasis, 2009）。

??? proof "证明"
    Vavasis 通过从 3-SAT 问题的一个变体进行多项式时间归约来证明。核心思想是将布尔满足问题编码为非负矩阵是否具有特定非负秩的判定问题。详细的归约构造较为技术性，这里略去完整证明。关键观察是：非负性约束使得NMF的组合结构远比普通矩阵分解更为复杂。$\blacksquare$

!!! example "例 58.1"
    考虑简单的 $3 \times 3$ 非负矩阵

    $$V = \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 1 & 1 & 0 \end{pmatrix}.$$

    $\mathrm{rank}(V) = 3$，但如果允许近似，NMF 可以用 $r = 2$ 来近似 $V$：

    $$W = \begin{pmatrix} 1 & 0.5 \\ 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad H = \begin{pmatrix} 1 & 0 & 0.5 \\ 0 & 1 & 1 \end{pmatrix},$$

    则 $WH = \begin{pmatrix} 1 & 0.5 & 1 \\ 0 & 1 & 1 \\ 1 & 0 & 0.5\end{pmatrix}$，近似误差 $\|V-WH\|_F^2 = 0.5$。

---

## 58.2 与 SVD/PCA 的本质区别

<div class="context-flow" markdown>

**核心问题**：NMF 与经典的 SVD/PCA 分解有什么根本性的区别？非负性约束带来了怎样的表示特性？

</div>

!!! theorem "定理 58.2 (SVD 的最优性)"
    截断 SVD 给出 Frobenius 范数下的最优低秩近似（Eckart-Young 定理）：对任意 $m \times n$ 矩阵 $A$，秩-$r$ 最优近似为

    $$A_r = \sum_{i=1}^r \sigma_i u_i v_i^\top, \quad \min_{\mathrm{rank}(B)\leq r}\|A - B\|_F^2 = \sum_{i=r+1}^{\min(m,n)} \sigma_i^2.$$

    但 SVD 的因子 $u_i, v_i$ 通常包含正负元素。

NMF 与 SVD 的核心区别在于表示的"语义"：

**1. 全局表示 vs. 部分表示。** SVD 的每个基向量 $u_i$ 是全局模式——它描述了数据在某个方向上的变化，包含正负分量。要重构一个数据点，需要加上某些基并减去另一些基。NMF 的基向量 $w_j$（$W$ 的列）是非负的"部件"，重构过程只涉及加法：$v_i \approx \sum_j h_{ji} w_j$，$h_{ji} \geq 0$。

**2. 正交性 vs. 可解释性。** SVD 的基向量是正交的，这在数学上优美但缺乏直观解释。NMF 的基向量通常不正交，但每个基向量对应一个可识别的"部件"。

!!! example "例 58.2"
    **人脸分解的对比。** 设 $V \in \mathbb{R}_{\geq 0}^{m \times n}$，每列是一张 $\sqrt{m} \times \sqrt{m}$ 人脸图像的像素向量。

    - **PCA/SVD**：特征脸 $u_1, u_2, \ldots$ 是全脸的模式，包含正负像素值，每张脸 $v_j \approx \bar{v} + \sum_i c_i u_i$。
    - **NMF**：基图像 $w_1, w_2, \ldots$ 是局部部件（眼睛区域、鼻子区域、嘴巴区域等），每张脸 $v_j \approx \sum_i h_{ij} w_i$，$h_{ij} \geq 0$。

    NMF 的这种部件表示使得我们可以说"这张脸 = 0.8 × 眼睛模板 + 1.2 × 鼻子模板 + 0.6 × 嘴巴模板"。

!!! theorem "定理 58.3 (NMF 近似的次优性)"
    对于非负矩阵 $V \geq 0$，NMF 给出的秩-$r$ 近似一般**劣于** SVD 的秩-$r$ 近似（在 Frobenius 范数意义下），即

    $$\min_{W \geq 0, H \geq 0}\|V - WH\|_F^2 \geq \|V - V_r\|_F^2 = \sum_{i>r}\sigma_i^2,$$

    其中 $V_r$ 是 SVD 截断。等号成立当且仅当 $V_r$ 本身可以表示为非负分解。

??? proof "证明"
    这直接来自 Eckart-Young 定理：$V_r$ 是所有秩不超过 $r$ 的矩阵中 Frobenius 范数意义下最接近 $V$ 的。$WH$（$W \geq 0, H \geq 0$）的秩不超过 $r$，所以 $\|V - WH\|_F \geq \|V - V_r\|_F$。$\blacksquare$

---

## 58.3 乘性更新算法

<div class="context-flow" markdown>

**核心问题**：如何设计保持非负性的迭代算法来求解 NMF？Lee 和 Seung 的乘性更新规则是怎样从梯度下降推导出来的？

</div>

!!! theorem "定理 58.4 (Lee-Seung 乘性更新规则——Frobenius 范数)"
    对于目标函数 $\|V - WH\|_F^2$，以下乘性更新规则使目标函数单调不增：

    $$H_{aj} \leftarrow H_{aj} \cdot \frac{(W^\top V)_{aj}}{(W^\top WH)_{aj}}, \qquad W_{ia} \leftarrow W_{ia} \cdot \frac{(VH^\top)_{ia}}{(WHH^\top)_{ia}}.$$

??? proof "证明"
    **推导思路：修正的梯度下降。**

    考虑固定 $W$，对 $H$ 优化 $f(H) = \|V - WH\|_F^2$。梯度为

    $$\nabla_H f = -2W^\top V + 2W^\top WH = 2(W^\top WH - W^\top V).$$

    将其分解为正部和负部：

    $$[\nabla_H f]^+ = 2W^\top WH, \quad [\nabla_H f]^- = 2W^\top V.$$

    标准梯度下降 $H \leftarrow H - \eta \nabla_H f$ 不能保证非负性。Lee 和 Seung 的关键思想是选取逐元素的自适应步长：

    $$\eta_{aj} = \frac{H_{aj}}{[\nabla_H f]^+_{aj}} = \frac{H_{aj}}{2(W^\top WH)_{aj}}.$$

    代入梯度下降公式：

    $$H_{aj} \leftarrow H_{aj} - \frac{H_{aj}}{2(W^\top WH)_{aj}} \cdot 2\bigl[(W^\top WH)_{aj} - (W^\top V)_{aj}\bigr]$$

    $$= H_{aj} \cdot \frac{(W^\top V)_{aj}}{(W^\top WH)_{aj}}.$$

    **单调性证明。** Lee 和 Seung 使用辅助函数方法（auxiliary function method）来严格证明单调性。定义辅助函数 $G(H, \tilde{H})$ 满足：

    (i) $G(H, \tilde{H}) \geq f(H)$，$\forall H, \tilde{H}$；

    (ii) $G(H, H) = f(H)$。

    则更新 $H^{(t+1)} = \arg\min_H G(H, H^{(t)})$ 保证 $f(H^{(t+1)}) \leq G(H^{(t+1)}, H^{(t)}) \leq G(H^{(t)}, H^{(t)}) = f(H^{(t)})$。

    对于 Frobenius 范数目标，合适的辅助函数为

    $$G(H, \tilde{H}) = f(\tilde{H}) + \mathrm{tr}\bigl[\nabla f(\tilde{H})^\top (H - \tilde{H})\bigr] + \sum_{a,j}\frac{(W^\top W \tilde{H})_{aj}}{\tilde{H}_{aj}}(H_{aj} - \tilde{H}_{aj})^2.$$

    其中最后一项用对角矩阵替换了 Hessian $W^\top W$（这是因为 $(W^\top W)_{ab} \leq \delta_{ab}(W^\top W\tilde{H})_{aj}/\tilde{H}_{aj}$ 在适当意义下成立）。

    对 $H$ 求导令其为零，恰好得到乘性更新规则。$\blacksquare$

!!! theorem "定理 58.5 (Lee-Seung 乘性更新规则——KL 散度)"
    对于 KL 散度目标 $D_{KL}(V \| WH)$，乘性更新规则为

    $$H_{aj} \leftarrow H_{aj} \cdot \frac{\sum_i W_{ia} V_{ij}/(WH)_{ij}}{\sum_i W_{ia}}, \qquad W_{ia} \leftarrow W_{ia} \cdot \frac{\sum_j H_{aj} V_{ij}/(WH)_{ij}}{\sum_j H_{aj}}.$$

!!! theorem "定理 58.6 (乘性更新的收敛性)"
    Lee-Seung 乘性更新序列 $\{(W^{(t)}, H^{(t)})\}$ 满足：

    (a) 目标函数序列 $\{f(W^{(t)}, H^{(t)})\}$ 单调不增且有下界（$\geq 0$），因此收敛；

    (b) 序列的每个聚点都是目标函数的**稳定点**（满足 KKT 条件）。

    但乘性更新不保证收敛到全局最优或局部最小值（可能收敛到鞍点）。

!!! example "例 58.3"
    **简单的乘性更新迭代。** 设 $V = \begin{pmatrix}5 & 3\\4 & 2\end{pmatrix}$，$r = 1$。初始化 $W^{(0)} = \begin{pmatrix}1\\1\end{pmatrix}$，$H^{(0)} = (1, 1)$。

    第一步更新 $H$：$W^\top V = (9, 5)$，$W^\top WH = (2, 2)$。

    $$H^{(1)} = H^{(0)} \odot \frac{W^\top V}{W^\top WH^{(0)}} = (1, 1) \odot \frac{(9, 5)}{(2, 2)} = (4.5, 2.5).$$

    第一步更新 $W$：$VH^{(1)\top} = \begin{pmatrix}5\cdot4.5+3\cdot2.5\\4\cdot4.5+2\cdot2.5\end{pmatrix} = \begin{pmatrix}30\\23\end{pmatrix}$。

    $WH^{(1)}H^{(1)\top} = \begin{pmatrix}1\\1\end{pmatrix}(4.5^2+2.5^2) = \begin{pmatrix}26.5\\26.5\end{pmatrix}$。

    $$W^{(1)} = W^{(0)} \odot \frac{VH^{(1)\top}}{W^{(0)}H^{(1)}H^{(1)\top}} = \begin{pmatrix}30/26.5\\23/26.5\end{pmatrix} \approx \begin{pmatrix}1.13\\0.87\end{pmatrix}.$$

    反复迭代将逐渐收敛。

---

## 58.4 交替非负最小二乘

<div class="context-flow" markdown>

**核心问题**：除了乘性更新，是否有更高效或收敛更快的 NMF 算法？交替最小二乘方法如何保持非负性？

</div>

!!! definition "定义 58.5 (交替非负最小二乘, ANLS)"
    ANLS 方法在每次迭代中交替求解两个非负最小二乘子问题：

    $$H^{(t+1)} = \arg\min_{H \geq 0} \|V - W^{(t)} H\|_F^2,$$

    $$W^{(t+1)} = \arg\min_{W \geq 0} \|V - W H^{(t+1)}\|_F^2.$$

!!! theorem "定理 58.7 (非负最小二乘子问题)"
    固定 $W$，子问题 $\min_{H \geq 0} \|V - WH\|_F^2$ 可以**按列分解**为 $n$ 个独立的非负最小二乘问题：

    $$\min_{h_j \geq 0} \|v_j - W h_j\|_2^2, \quad j = 1, \ldots, n,$$

    其中 $v_j$ 和 $h_j$ 分别是 $V$ 和 $H$ 的第 $j$ 列。每个子问题是凸的二次规划。

??? proof "证明"
    $\|V - WH\|_F^2 = \sum_{j=1}^n \|v_j - Wh_j\|_2^2$，各列独立，所以可以分开优化。每个 $\min_{h_j \geq 0}\|v_j - Wh_j\|_2^2$ 是凸二次函数在非负象限上的最小化，即标准的非负最小二乘问题。$\blacksquare$

!!! theorem "定理 58.8 (非负最小二乘的 KKT 条件)"
    $h^* \geq 0$ 是 $\min_{h \geq 0}\|v - Wh\|_2^2$ 的解当且仅当满足 KKT 条件：

    $$W^\top(v - Wh^*) \geq 0, \quad h^* \geq 0, \quad h^*_a \cdot [W^\top(v - Wh^*)]_a = 0, \quad \forall a.$$

    即：梯度的每个分量要么为零（当 $h^*_a > 0$ 时），要么非负（当 $h^*_a = 0$ 时）。

!!! definition "定义 58.6 (活跃集方法)"
    **Lawson-Hanson 活跃集方法**是求解非负最小二乘的经典算法：

    1. 维护活跃集 $\mathcal{A}$（$h_a = 0$ 的指标集）和自由集 $\mathcal{F}$（$h_a > 0$ 的指标集）。
    2. 在自由集上求解无约束最小二乘。
    3. 如果有自由变量变负，将其移入活跃集；如果梯度指示某活跃变量应变正，将其移入自由集。
    4. 反复直到 KKT 条件满足。

!!! theorem "定理 58.9 (ANLS 的收敛性)"
    ANLS 方法生成的序列 $\{(W^{(t)}, H^{(t)})\}$ 满足：

    (a) 目标函数 $\|V - W^{(t)}H^{(t)}\|_F^2$ 单调不增；

    (b) 在温和的正则性条件下（例如 $W^\top W$ 正定），序列的每个聚点都满足 KKT 条件。

!!! example "例 58.4"
    **ANLS vs. 乘性更新的收敛速度比较。** 在典型的文本数据实验中（$m = 5000$ 词汇，$n = 1000$ 文档，$r = 50$），ANLS 通常在 20-50 次迭代后收敛到高精度解，而乘性更新可能需要数百次迭代。但 ANLS 每次迭代的计算量更大（需要求解非负最小二乘），所以总体时间不一定更快。

---

## 58.5 唯一性条件

<div class="context-flow" markdown>

**核心问题**：NMF 分解是否唯一？什么条件下唯一？非唯一性对应用有何影响？

</div>

NMF 分解一般**不唯一**。最明显的非唯一性来自缩放和置换。

!!! theorem "定理 58.10 (缩放和置换不确定性)"
    如果 $V = WH$，则对任意正对角矩阵 $D = \mathrm{diag}(d_1, \ldots, d_r)$（$d_i > 0$）和置换矩阵 $P$，

    $$V = (WDP)(P^\top D^{-1}H) = \tilde{W}\tilde{H},$$

    其中 $\tilde{W} = WDP \geq 0$，$\tilde{H} = P^\top D^{-1}H \geq 0$。

??? proof "证明"
    $\tilde{W}\tilde{H} = WDP \cdot P^\top D^{-1}H = W \cdot DD^{-1} \cdot H = WH = V$。非负性由 $D > 0$、$P \geq 0$ 保证。$\blacksquare$

除了这种平凡的不确定性外，NMF 可能有本质上不同的分解。

!!! definition "定义 58.7 (可分矩阵)"
    非负矩阵 $V \in \mathbb{R}_{\geq 0}^{m \times n}$ 称为**可分的**（separable），如果存在秩-$r$ 非负分解 $V = WH$，使得 $W$ 的列集合包含 $V$ 的 $r$ 个列。等价地，存在指标集 $\mathcal{K} \subset \{1,\ldots,n\}$，$|\mathcal{K}| = r$，使得

    $$V = V(:, \mathcal{K})\, H, \quad H \geq 0.$$

!!! theorem "定理 58.11 (可分 NMF 的唯一性——Donoho-Stodden)"
    设 $V = WH \geq 0$，$\mathrm{rank}(V) = r$。如果 $V$ 是可分的，且 $W$ 的列在**简单体条件**下是唯一可确定的：$W$ 的每一列与 $V$ 的某一列成正比，且 $H$ 的每一行都有一个分量足够大（使得对应的 $V$ 列"接近纯"），则 $W$ 在缩放和置换意义下唯一。

    具体地，Donoho 和 Stodden (2003) 证明：如果数据点在非负象限中的分布满足"简单体几何"条件——数据的凸包是一个单纯形，且 $W$ 的列对应单纯形的顶点——则 NMF 在缩放和置换意义下唯一。

!!! theorem "定理 58.12 (充分稀疏条件下的唯一性)"
    设 $V = WH$，$W \in \mathbb{R}_{\geq 0}^{m \times r}$，$H \in \mathbb{R}_{\geq 0}^{r \times n}$。如果 $W$ 和 $H$ 满足以下**充分稀疏条件**：

    (a) $W$ 的每一列都有足够多的零元素（$W$ 足够稀疏）；

    (b) $H$ 是满秩的（$\mathrm{rank}(H) = r$）；

    则在缩放和置换意义下，NMF 分解是唯一的。

!!! example "例 58.5"
    考虑 $V = \begin{pmatrix}1 & 0 & 1\\0 & 1 & 1\end{pmatrix}$。这是可分矩阵，$r = 2$：

    $$V = \begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}\begin{pmatrix}1 & 0 & 1\\0 & 1 & 1\end{pmatrix} = I \cdot V.$$

    由于 $W = I$ 的列对应 $V$ 的前两列，且 $H = V$ 每行都有非零元素，此分解在缩放和置换意义下唯一。

---

## 58.6 非负秩

<div class="context-flow" markdown>

**核心问题**：非负矩阵分解中的最小分解秩是多少？非负秩与普通秩有何关系？

</div>

!!! definition "定义 58.8 (非负秩)"
    非负矩阵 $V \in \mathbb{R}_{\geq 0}^{m \times n}$ 的**非负秩**定义为

    $$\mathrm{rank}_+(V) = \min\{r \in \mathbb{N} : \exists\, W \in \mathbb{R}_{\geq 0}^{m \times r},\, H \in \mathbb{R}_{\geq 0}^{r \times n},\, V = WH\}.$$

!!! theorem "定理 58.13 (非负秩与普通秩的关系)"
    对任意非负矩阵 $V \in \mathbb{R}_{\geq 0}^{m \times n}$：

    (a) $\mathrm{rank}_+(V) \geq \mathrm{rank}(V)$。

    (b) $\mathrm{rank}_+(V) \leq \min(m, n)$。

    (c) 不等式 (a) 可以严格成立，且差距可以任意大。

??? proof "证明"
    **(a)** 若 $V = WH$，$W \in \mathbb{R}^{m \times r}$，$H \in \mathbb{R}^{r \times n}$，则 $\mathrm{rank}(V) \leq \min(\mathrm{rank}(W), \mathrm{rank}(H)) \leq r$。对所有可行的 $r$ 取最小值，得 $\mathrm{rank}(V) \leq \mathrm{rank}_+(V)$。

    **(b)** 将 $V$ 的每列表示为标准基向量的非负组合：$V = I_m \cdot V$，给出 $\mathrm{rank}_+(V) \leq m$。类似地 $\mathrm{rank}_+(V) \leq n$。

    **(c)** 经典的例子是以下矩阵（见例 58.6）。$\blacksquare$

!!! example "例 58.6"
    **非负秩严格大于普通秩的经典例子。** 考虑 $4 \times 4$ 的"正斜杠"矩阵

    $$V = \begin{pmatrix}1 & 1 & 0 & 0\\1 & 0 & 1 & 0\\0 & 1 & 0 & 1\\0 & 0 & 1 & 1\end{pmatrix}.$$

    可以验证 $\mathrm{rank}(V) = 3$。但 $\mathrm{rank}_+(V) = 4$：不存在 $W \in \mathbb{R}_{\geq 0}^{4 \times 3}$ 和 $H \in \mathbb{R}_{\geq 0}^{3 \times 4}$ 使得 $V = WH$。

    直觉上，$V$ 的四个列在非负象限中构成一个"不可简化"的结构——无法找到三个非负向量的非负组合来同时表示全部四列。

!!! theorem "定理 58.14 (计算非负秩的困难性)"
    计算非负矩阵的非负秩是 **NP 困难**的。更精确地：

    (a) 判定 $\mathrm{rank}_+(V) \leq r$ 是 NP 完全的（当 $r$ 是输入的一部分时）。

    (b) 存在 $n \times n$ 非负矩阵使得 $\mathrm{rank}_+(V) = 2^{\Omega(n^{1/3})}$ 而 $\mathrm{rank}(V) = O(n)$。

!!! definition "定义 58.9 (布尔秩)"
    非负矩阵理论中另一个相关概念是**布尔秩**：$V$ 的布尔秩是最小的 $r$，使得存在 0-1 矩阵 $W \in \{0,1\}^{m \times r}$，$H \in \{0,1\}^{r \times n}$，满足 $V = W \odot_{\mathrm{Bool}} H$（布尔矩阵乘法，$\lor$ 代替 $+$，$\land$ 代替 $\times$）。布尔秩与非负秩有类似的复杂性。

---

## 58.7 稀疏与正则化 NMF

<div class="context-flow" markdown>

**核心问题**：如何通过正则化来改善 NMF 的唯一性和可解释性？常见的正则化策略有哪些？

</div>

!!! definition "定义 58.10 (稀疏 NMF)"
    **稀疏 NMF** 在标准 NMF 目标上加入 $L_1$ 正则化项以促进稀疏性：

    $$\min_{W \geq 0, H \geq 0}\, \|V - WH\|_F^2 + \alpha \sum_{i,a}|W_{ia}| + \beta \sum_{a,j}|H_{aj}|.$$

    由于 $W, H \geq 0$，$L_1$ 范数简化为元素之和：$\sum |W_{ia}| = \sum W_{ia} = \mathbf{1}^\top W \mathbf{1}$。

!!! theorem "定理 58.15 (稀疏 NMF 的乘性更新)"
    对稀疏 NMF 目标，修改后的乘性更新规则为

    $$H_{aj} \leftarrow H_{aj} \cdot \frac{(W^\top V)_{aj}}{(W^\top WH + \beta)_{aj}}, \qquad W_{ia} \leftarrow W_{ia} \cdot \frac{(VH^\top)_{ia}}{(WHH^\top + \alpha)_{ia}}.$$

??? proof "证明"
    对 $H$ 的梯度为 $\nabla_H f = 2(W^\top WH - W^\top V) + \beta \mathbf{1}$。正部为 $2W^\top WH + \beta$，负部为 $2W^\top V$。使用与标准乘性更新相同的自适应步长策略即得。$\blacksquare$

!!! definition "定义 58.11 (正交 NMF)"
    **正交 NMF**（ONMF）在 NMF 的基础上添加正交约束：

    $$\min_{W \geq 0, H \geq 0}\, \|V - WH\|_F^2 \quad \text{s.t.} \quad H H^\top = I_r.$$

    正交 NMF 与 $k$-均值聚类有密切关系：当 $H$ 的每行恰好有一个非零元素时，NMF 退化为将 $V$ 的列聚类为 $r$ 组。

!!! theorem "定理 58.16 (ONMF 与聚类的等价性)"
    当 $H \geq 0$ 且 $HH^\top = I$ 时，$H$ 的每一行至多有一个非零元素。因此，$V \approx WH$ 意味着 $V$ 的每一列被"指派"到 $W$ 的某一列，这等价于 $k$-均值聚类。

??? proof "证明"
    设 $H$ 的第 $a$ 行为 $h_a^\top \in \mathbb{R}^{1 \times n}$。正交条件 $HH^\top = I$ 意味着 $h_a^\top h_b = \delta_{ab}$。非负性 $H \geq 0$ 意味着 $h_a$ 的分量非负。

    $h_a^\top h_b = 0$（$a \neq b$）加上非负性意味着 $h_a$ 和 $h_b$ 的支撑集不相交。因此 $H$ 的列（对应数据点）在行方向上最多在一个位置非零，即每个数据点恰好被指派到一个聚类。$\blacksquare$

!!! definition "定义 58.12 (半非负矩阵分解)"
    **半 NMF**（Semi-NMF）放松了对 $W$ 的非负约束：

    $$\min_{W, H \geq 0}\, \|V - WH\|_F^2, \quad W \in \mathbb{R}^{m \times r}, \quad H \in \mathbb{R}_{\geq 0}^{r \times n}.$$

    当数据矩阵 $V$ 本身可能包含负元素时（如中心化后的数据），半 NMF 更为合适。

!!! definition "定义 58.13 (对称 NMF)"
    当 $V$ 是对称非负矩阵时，**对称 NMF** 寻找

    $$\min_{W \geq 0}\, \|V - WW^\top\|_F^2.$$

    对称 NMF 与图聚类和社区发现密切相关。

!!! example "例 58.7"
    **稀疏正则化的效果。** 对文本数据矩阵 $V$（词频矩阵），增加 $H$ 上的 $L_1$ 正则化 $\beta = 0.1$ 后：

    - 未正则化：$H$ 平均每行 80% 元素非零，主题表示模糊。
    - $L_1$ 正则化：$H$ 平均每行 20% 元素非零，每个文档只与少数主题关联，语义更清晰。

---

## 58.8 应用

<div class="context-flow" markdown>

**核心问题**：NMF 在实际应用中如何发挥其"部分-整体"表示的优势？

</div>

### 58.8.1 主题模型

!!! example "例 58.8"
    **文本挖掘中的主题建模。** 设 $V \in \mathbb{R}_{\geq 0}^{m \times n}$ 是**词-文档矩阵**：$V_{ij}$ 是第 $i$ 个词在第 $j$ 个文档中的出现频率（TF-IDF 权重）。

    NMF 分解 $V \approx WH$ 的解释：

    - $W$ 的第 $a$ 列 $w_a \in \mathbb{R}_{\geq 0}^m$ 是第 $a$ 个**主题**的词分布：$w_{ia}$ 越大，词 $i$ 在主题 $a$ 中越重要。
    - $H$ 的第 $j$ 列 $h_j \in \mathbb{R}_{\geq 0}^r$ 是文档 $j$ 的**主题混合系数**：$h_{aj}$ 越大，文档 $j$ 越多地涉及主题 $a$。
    - 重构 $v_j \approx Wh_j = \sum_{a=1}^r h_{aj} w_a$：文档 $j$ 的词频近似为各主题词分布的加权和。

    NMF 的非负性保证了所有权重都是非负的，使得"文档 = 主题的加权混合"这一解释在语义上合理。

### 58.8.2 音源分离

!!! example "例 58.9"
    **音乐信号的音源分离。** 给定混合音频信号的短时 Fourier 变换幅度谱 $V \in \mathbb{R}_{\geq 0}^{m \times n}$（$m$ 个频率 bin，$n$ 个时间帧），NMF 分解

    $$V \approx WH$$

    - $W$ 的每一列是一个**频谱模板**（spectral template），代表一种音源（如钢琴、人声、鼓声）的典型频率特征。
    - $H$ 的每一行是对应音源的**时间激活模式**（temporal activation）。
    - 通过 Wiener 滤波：分离后第 $a$ 个音源的幅度谱为 $V \odot \frac{w_a h_a^\top}{WH}$（逐元素运算）。

    选择 $r$ 等于音源数（如 $r = 3$），KL 散度目标通常比 Frobenius 范数更适合音频数据。

### 58.8.3 高光谱解混

!!! example "例 58.10"
    **高光谱图像解混。** 高光谱图像 $V \in \mathbb{R}_{\geq 0}^{m \times n}$（$m$ 个光谱波段，$n$ 个像素）的每列 $v_j$ 是像素 $j$ 的光谱向量。由于每个像素可能包含多种地物类型：

    $$v_j \approx \sum_{a=1}^r h_{aj}\, w_a, \quad h_{aj} \geq 0,$$

    其中 $w_a$ 是第 $a$ 种**端元**（endmember）的纯光谱，$h_{aj}$ 是其在像素 $j$ 中的**丰度**（abundance）。通常还要求 $\sum_a h_{aj} = 1$（丰度约束）。

    NMF 的非负性约束完美匹配了物理约束（反射率和丰度都非负）。

### 58.8.4 人脸部件分解

!!! example "例 58.11"
    **Lee-Seung 的经典人脸实验。** 对 ORL 人脸数据集（$m = 92 \times 112 = 10304$ 像素，$n = 400$ 张脸，$r = 49$ 个基）：

    - PCA 得到的 49 个特征脸是全脸的"鬼影"模式，每个基包含正负像素。
    - NMF 得到的 49 个基图像是局部化的面部部件：眉毛、眼睛轮廓、鼻尖、嘴角等。

    这种部件表示对于遮挡鲁棒性更强——如果眼睛被遮挡，只影响与眼睛相关的基的系数，而不影响其他部件。

---

**本章要点总结：**

1. NMF 通过非负性约束实现"部分-整体"的加法表示，与 SVD/PCA 的"全局-对消"表示形成对比。
2. Lee-Seung 乘性更新是最经典的 NMF 算法，保证非负性和目标函数单调下降。
3. ANLS 在每步求解凸子问题，收敛性更好但每步计算量更大。
4. NMF 一般不唯一；可分性和稀疏性是保证唯一性的重要条件。
5. 非负秩可以远大于普通秩，且计算是 NP 困难的。
6. 稀疏、正交、对称等正则化变体在特定应用中各有优势。
7. NMF 在主题建模、音源分离、高光谱解混和人脸分析中有广泛应用。
