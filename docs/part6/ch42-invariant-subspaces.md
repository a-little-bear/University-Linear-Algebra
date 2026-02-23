# 第 42 章 不变子空间与扰动

<div class="context-flow" markdown>

**前置**：向量空间(Ch4) · 线性变换(Ch5) · 特征值(Ch6) · Jordan形(Ch12) · 范数与扰动(Ch15)

**本章脉络**：不变子空间定义 → 不变子空间格 → 超不变子空间 → 约化子空间 → 互补不变子空间 → 子空间之间的角度 → Davis-Kahan sin Θ 定理 → 谱投影扰动 → Wedin sin Θ 定理 → Stewart tan Θ 定理 → Rosenblum 定理 → 无穷维不变子空间问题

**延伸**：Davis-Kahan 定理是现代统计学和机器学习中 PCA 扰动分析的理论基石；不变子空间理论在算子代数（von Neumann 代数中的投影格）中有无穷维推广

</div>

不变子空间是线性代数中最核心的结构概念之一。当我们研究一个线性变换 $T: V \to V$ 时，寻找使 $T$ 的作用"封闭"的子空间，本质上是在寻找将 $T$ 分解为更简单部分的途径。特征空间是最基本的不变子空间，而 Jordan 标准形理论则展示了如何通过广义特征空间实现矩阵的精细分解。

本章从不变子空间的基本定义出发，建立不变子空间格的代数结构，深入讨论超不变子空间与约化子空间这两类特殊的不变子空间，然后转向一个截然不同但极其重要的方向——子空间之间的"距离"和"角度"。Davis-Kahan sin Θ 定理将不变子空间的扰动与矩阵的扰动精确联系起来，它不仅是矩阵扰动理论的瑰宝，更是现代数据科学中主成分分析（PCA）理论保证的数学基础。

---

## 42.1 不变子空间的定义与基本性质

<div class="context-flow" markdown>

**核心问题**：什么样的子空间在线性变换的作用下保持"封闭"？不变子空间的全体构成怎样的代数结构？

</div>

!!! definition "定义 42.1 ($T$-不变子空间)"
    设 $V$ 是域 $\mathbb{F}$ 上的向量空间，$T: V \to V$ 是线性变换。子空间 $\mathcal{M} \subseteq V$ 称为 **$T$-不变子空间**（$T$-invariant subspace），如果
    $$T(\mathcal{M}) \subseteq \mathcal{M},$$
    即对任意 $v \in \mathcal{M}$，有 $Tv \in \mathcal{M}$。

    等价地，$T$ 在 $\mathcal{M}$ 上的限制 $T|_{\mathcal{M}}: \mathcal{M} \to \mathcal{M}$ 是良定义的线性变换。

在矩阵语言中，如果 $A \in \mathbb{F}^{n \times n}$ 且 $\mathcal{M}$ 是 $\mathbb{F}^n$ 的一个 $k$ 维 $A$-不变子空间，取 $\mathcal{M}$ 的一组基排成矩阵 $X \in \mathbb{F}^{n \times k}$，则 $A$-不变性等价于存在矩阵 $B \in \mathbb{F}^{k \times k}$ 使得

$$AX = XB.$$

矩阵 $B$ 就是 $A|_{\mathcal{M}}$ 在所选基下的矩阵表示。

!!! example "例 42.1 (基本不变子空间)"
    设 $A \in \mathbb{C}^{n \times n}$。以下子空间都是 $A$-不变的：

    1. **平凡子空间**：$\{0\}$ 和 $\mathbb{C}^n$。
    2. **特征空间**：若 $\lambda$ 是 $A$ 的特征值，则 $\ker(A - \lambda I)$ 是 $A$-不变的。
    3. **广义特征空间**：$\ker(A - \lambda I)^k$ 对每个 $k \geq 1$ 都是 $A$-不变的。
    4. **值域与核**：$\operatorname{Im}(A)$ 和 $\ker(A)$ 都是 $A$-不变的。
    5. **$A$-循环子空间**：对任意 $v \in \mathbb{C}^n$，$\mathcal{K}(A, v) = \operatorname{span}\{v, Av, A^2v, \ldots\}$ 是 $A$-不变的。
    6. **多项式子空间**：若 $p(t)$ 是多项式，则 $\ker(p(A))$ 和 $\operatorname{Im}(p(A))$ 都是 $A$-不变的。

!!! theorem "定理 42.1 (不变子空间的等价刻画)"
    设 $A \in \mathbb{F}^{n \times n}$，$\mathcal{M}$ 是 $\mathbb{F}^n$ 的 $k$ 维子空间。以下条件等价：

    1. $\mathcal{M}$ 是 $A$-不变的。
    2. 存在 $B \in \mathbb{F}^{k \times k}$，使得 $AX = XB$，其中 $X$ 的列是 $\mathcal{M}$ 的一组基。
    3. 存在可逆矩阵 $P$ 使得 $P^{-1}AP$ 具有分块上三角形式
    $$P^{-1}AP = \begin{pmatrix} B & C \\ 0 & D \end{pmatrix},$$
    其中 $B \in \mathbb{F}^{k \times k}$，且 $P$ 的前 $k$ 列张成 $\mathcal{M}$。

??? proof "证明"
    **(1) $\Rightarrow$ (2)**：设 $\{x_1, \ldots, x_k\}$ 是 $\mathcal{M}$ 的一组基。由 $A$-不变性，每个 $Ax_j \in \mathcal{M}$，因此
    $$Ax_j = \sum_{i=1}^{k} b_{ij} x_i, \quad j = 1, \ldots, k.$$
    令 $X = (x_1, \ldots, x_k)$，$B = (b_{ij})$，则 $AX = XB$。

    **(2) $\Rightarrow$ (3)**：将 $X$ 的列扩充为 $\mathbb{F}^n$ 的一组基，设扩充部分构成矩阵 $Y$。令 $P = (X \mid Y)$，则 $P$ 可逆，且
    $$P^{-1}AP = \begin{pmatrix} B & C \\ 0 & D \end{pmatrix},$$
    这是因为 $AP = P \cdot P^{-1}AP$ 的前 $k$ 列给出 $AX = XB + Y \cdot 0 = XB$。

    **(3) $\Rightarrow$ (1)**：设 $P = (X \mid Y)$，其中 $X$ 的列张成 $\mathcal{M}$。由 $P^{-1}AP$ 的形式，$AX = XB$，因此对任意 $v = Xc \in \mathcal{M}$，$Av = AXc = XBc \in \mathcal{M}$。$\blacksquare$

!!! definition "定义 42.2 (不变子空间格 $\operatorname{Lat}(T)$)"
    设 $T: V \to V$ 是线性变换。$T$ 的所有不变子空间构成的集合，按子空间包含关系 $\subseteq$ 构成偏序集，记为 $\operatorname{Lat}(T)$，称为 $T$ 的**不变子空间格**（invariant subspace lattice）。

    在 $\operatorname{Lat}(T)$ 中，格运算定义为：

    - **交**（meet）：$\mathcal{M}_1 \wedge \mathcal{M}_2 = \mathcal{M}_1 \cap \mathcal{M}_2$，
    - **并**（join）：$\mathcal{M}_1 \vee \mathcal{M}_2 = \mathcal{M}_1 + \mathcal{M}_2$。

!!! theorem "定理 42.2 ($\operatorname{Lat}(T)$ 构成完备格)"
    对任意线性变换 $T: V \to V$，$\operatorname{Lat}(T)$ 在上述运算下构成完备格。具体地：

    1. 若 $\{\mathcal{M}_\alpha\}_{\alpha \in \Lambda}$ 是 $T$-不变子空间的任意族，则 $\bigcap_\alpha \mathcal{M}_\alpha$ 和 $\sum_\alpha \mathcal{M}_\alpha$ 都是 $T$-不变的。
    2. 最小元素是 $\{0\}$，最大元素是 $V$。

??? proof "证明"
    (1) 设 $v \in \bigcap_\alpha \mathcal{M}_\alpha$。则对每个 $\alpha$，$v \in \mathcal{M}_\alpha$，因此 $Tv \in \mathcal{M}_\alpha$（因为 $\mathcal{M}_\alpha$ 是 $T$-不变的）。从而 $Tv \in \bigcap_\alpha \mathcal{M}_\alpha$。

    对于和空间，设 $v = \sum_{i=1}^{m} v_{\alpha_i}$，其中 $v_{\alpha_i} \in \mathcal{M}_{\alpha_i}$。则 $Tv = \sum_{i=1}^{m} Tv_{\alpha_i}$，其中每个 $Tv_{\alpha_i} \in \mathcal{M}_{\alpha_i}$。因此 $Tv \in \sum_\alpha \mathcal{M}_\alpha$。

    (2) 显然 $T(\{0\}) = \{0\} \subseteq \{0\}$，且 $T(V) \subseteq V$。$\blacksquare$

!!! example "例 42.2 (不变子空间格的例子)"
    考虑 $A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}$ 在 $\mathbb{C}^2$ 上的作用。

    特征值为 $\lambda_1 = 2$，$\lambda_2 = 3$。对应特征向量分别为 $v_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$，$v_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。

    $\operatorname{Lat}(A) = \bigl\{ \{0\},\; \operatorname{span}\{v_1\},\; \operatorname{span}\{v_2\},\; \mathbb{C}^2 \bigr\}$。

    这是因为 $\mathbb{C}^2$ 中任何一维子空间若是 $A$-不变的，其生成元必须是特征向量。由于 $A$ 有两个不同的特征值，恰好有两个一维不变子空间。

    而对于 $B = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$（单一特征值 $\lambda = 2$，Jordan 块大小 2），唯一的一维不变子空间是 $\operatorname{span}\left\{\begin{pmatrix} 1 \\ 0 \end{pmatrix}\right\}$，因此
    $$\operatorname{Lat}(B) = \left\{ \{0\},\; \operatorname{span}\left\{\begin{pmatrix} 1 \\ 0 \end{pmatrix}\right\},\; \mathbb{C}^2 \right\}.$$

---

## 42.2 超不变子空间

<div class="context-flow" markdown>

**核心问题**：哪些不变子空间在所有与 $T$ 交换的算子作用下都保持不变？这类子空间有何结构？

</div>

!!! definition "定义 42.3 (超不变子空间)"
    设 $T: V \to V$ 是线性变换。子空间 $\mathcal{M} \subseteq V$ 称为 $T$ 的**超不变子空间**（hyperinvariant subspace），如果对任意与 $T$ 交换的线性变换 $S$（即 $ST = TS$），$\mathcal{M}$ 都是 $S$-不变的。

    等价地，$\mathcal{M}$ 对 $T$ 的**中心化子**（centralizer）$\mathcal{C}(T) = \{S \in \operatorname{End}(V) : ST = TS\}$ 中的每个元素都不变。

超不变子空间显然是不变子空间（因为 $T \in \mathcal{C}(T)$），但反之不一定成立。超不变子空间代表了 $T$ 的"内禀"结构，不依赖于特定的交换算子选取。

!!! theorem "定理 42.3 (超不变子空间的 Jordan 形刻画)"
    设 $A \in \mathbb{C}^{n \times n}$。子空间 $\mathcal{M}$ 是 $A$ 的超不变子空间，当且仅当 $\mathcal{M}$ 可以表示为以下形式：

    对 $A$ 的每个不同特征值 $\lambda_i$（$i = 1, \ldots, s$），设 $A$ 在 $\lambda_i$ 处的 Jordan 块大小为 $n_{i,1} \geq n_{i,2} \geq \cdots \geq n_{i,r_i}$。则超不变子空间 $\mathcal{M}$ 必然是某些**根子空间**（root subspaces）的直和：
    $$\mathcal{M} = \bigoplus_{i=1}^{s} \mathcal{M}_i,$$
    其中每个 $\mathcal{M}_i$ 是 $\ker(A - \lambda_i I)^{k_i}$ 对某个 $0 \leq k_i \leq n_{i,1}$ 的子空间，并且 $\mathcal{M}_i$ 必须由完整的 Jordan 链（或其截断）生成。

??? proof "证明"
    关键思路是利用 Jordan 标准形将问题简化。

    不妨设 $A$ 已化为 Jordan 标准形 $J$。首先，根子空间分解 $V = \bigoplus_{i=1}^{s} V_i$（其中 $V_i = \ker(A - \lambda_i I)^{m_i}$，$m_i$ 是 $\lambda_i$ 的代数重数）给出的每个 $V_i$ 都是超不变的——这是因为任何与 $A$ 交换的 $S$ 必须保持 $V_i$ 不变（若 $(A - \lambda_i I)^{m_i} v = 0$，则 $(A - \lambda_i I)^{m_i} Sv = S(A - \lambda_i I)^{m_i} v = 0$）。

    因此 $\mathcal{M} = \bigoplus_{i=1}^{s} (\mathcal{M} \cap V_i)$，每个分量 $\mathcal{M} \cap V_i$ 必须是 $A|_{V_i}$ 的超不变子空间。这将问题归结为单一特征值的情形。

    对于单一特征值 $\lambda$（不妨设 $\lambda = 0$），$A$ 是幂零矩阵。此时需要证明超不变子空间恰好是 $\ker(A^k)$ 的形式。这可以通过分析幂零矩阵中心化子的结构来完成：与幂零 Jordan 块交换的矩阵的精确描述表明，只有 $\ker(A^k)$ 这样的子空间才能同时在所有交换矩阵作用下保持不变。$\blacksquare$

!!! example "例 42.3"
    设 $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$（$3 \times 3$ 幂零 Jordan 块）。

    **不变子空间**：$\{0\}$，$\operatorname{span}\{e_1\}$，$\operatorname{span}\{e_1, e_2\}$，$\mathbb{C}^3$。

    **超不变子空间**：$\{0\} = \ker(A^0)$，$\operatorname{span}\{e_1\} = \ker(A)$，$\operatorname{span}\{e_1, e_2\} = \ker(A^2)$，$\mathbb{C}^3 = \ker(A^3)$。

    在这个例子中，所有不变子空间恰好都是超不变的。但如果 Jordan 形有多个相同特征值的 Jordan 块，情况会更复杂。例如 $A = \operatorname{diag}(J_2(0), J_1(0))$，其中 $J_2(0) = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，则 $\operatorname{span}\{e_1, e_3\} = \ker(A)$ 是超不变的，但 $\operatorname{span}\{e_1\}$ 和 $\operatorname{span}\{e_3\}$ 虽然是不变的，却不是超不变的。

---

## 42.3 约化子空间与互补不变子空间

<div class="context-flow" markdown>

**核心问题**：何时一个不变子空间能够将线性变换真正"分解"为两个独立部分？

</div>

!!! definition "定义 42.4 (约化子空间)"
    设 $V$ 是有限维内积空间，$T: V \to V$ 是线性变换。子空间 $\mathcal{M}$ 称为 $T$ 的**约化子空间**（reducing subspace），如果 $\mathcal{M}$ 和 $\mathcal{M}^\perp$ 都是 $T$-不变的。

    等价地，$\mathcal{M}$ 约化 $T$ 当且仅当 $T$ 与正交投影 $P_{\mathcal{M}}$ 交换：$TP_{\mathcal{M}} = P_{\mathcal{M}}T$。

!!! theorem "定理 42.4 (约化子空间的等价条件)"
    设 $A \in \mathbb{C}^{n \times n}$，$\mathcal{M}$ 是 $\mathbb{C}^n$ 的子空间。以下条件等价：

    1. $\mathcal{M}$ 同时是 $A$-不变和 $A^*$-不变的。
    2. $\mathcal{M}$ 和 $\mathcal{M}^\perp$ 都是 $A$-不变的。
    3. $AP_{\mathcal{M}} = P_{\mathcal{M}} A$，其中 $P_{\mathcal{M}}$ 是到 $\mathcal{M}$ 的正交投影。
    4. 存在酉矩阵 $U$ 使得 $U^*AU = \begin{pmatrix} B & 0 \\ 0 & D \end{pmatrix}$，其中 $U$ 的前 $k$ 列（$k = \dim \mathcal{M}$）张成 $\mathcal{M}$。

??? proof "证明"
    **(1) $\Leftrightarrow$ (2)**：若 $\mathcal{M}$ 是 $A$-不变的，则对 $v \in \mathcal{M}$，$w \in \mathcal{M}^\perp$，有 $\langle w, Av \rangle = 0$。因此 $\langle A^*w, v \rangle = 0$，即 $A^*w \in \mathcal{M}^\perp$。这说明 $\mathcal{M}^\perp$ 是 $A^*$-不变的。类似地，$\mathcal{M}^\perp$ 的 $A$-不变性等价于 $\mathcal{M}$ 的 $A^*$-不变性。

    **(2) $\Leftrightarrow$ (3)**：设 $P = P_{\mathcal{M}}$。对任意 $v \in \mathbb{C}^n$，$v = Pv + (I-P)v$，其中 $Pv \in \mathcal{M}$，$(I-P)v \in \mathcal{M}^\perp$。若 $\mathcal{M}$ 和 $\mathcal{M}^\perp$ 都是 $A$-不变的，则 $APv \in \mathcal{M}$，$A(I-P)v \in \mathcal{M}^\perp$。因此 $PAv = PAPv + PA(I-P)v = APv + 0 = APv$，即 $PA = AP$。反之亦然。

    **(3) $\Leftrightarrow$ (4)**：这是直接推论，取 $\mathcal{M}$ 和 $\mathcal{M}^\perp$ 的标准正交基排成酉矩阵 $U$ 即可。$\blacksquare$

!!! definition "定义 42.5 (互补不变子空间)"
    设 $T: V \to V$ 是线性变换。若 $V = \mathcal{M} \oplus \mathcal{N}$，且 $\mathcal{M}$ 和 $\mathcal{N}$ 都是 $T$-不变的，则称 $\mathcal{N}$ 是 $\mathcal{M}$ 关于 $T$ 的**互补不变子空间**（complementary invariant subspace）。

!!! theorem "定理 42.5 (互补不变子空间的存在性)"
    设 $A \in \mathbb{C}^{n \times n}$，$\mathcal{M}$ 是 $A$-不变子空间。

    1. **一般情形**：互补不变子空间不一定存在。
    2. **半单情形**：若 $A$ 是可对角化的（半单的），则互补不变子空间总是存在的。
    3. **谱分离情形**：若 $\sigma(A|_{\mathcal{M}}) \cap \sigma(A|_{V/\mathcal{M}}) = \emptyset$（$\mathcal{M}$ 上和商空间上的谱不相交），则互补不变子空间存在且唯一。

??? proof "证明"
    **(1)** 反例：设 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$\mathcal{M} = \operatorname{span}\{e_1\} = \ker(A)$。若存在 $A$-不变的 $\mathcal{N}$ 使得 $\mathbb{C}^2 = \mathcal{M} \oplus \mathcal{N}$，则 $\mathcal{N}$ 是一维的，设 $\mathcal{N} = \operatorname{span}\{v\}$，$v = \begin{pmatrix} a \\ b \end{pmatrix}$，$b \neq 0$。$A$-不变性要求 $Av = \begin{pmatrix} b \\ 0 \end{pmatrix} \in \mathcal{N}$，但 $\begin{pmatrix} b \\ 0 \end{pmatrix} \in \mathcal{M}$，$b \neq 0$，矛盾于 $\mathcal{M} \cap \mathcal{N} = \{0\}$。

    **(3)** 由谱分离条件，可以用 Riesz 投影（环路积分）构造：设 $\Gamma$ 是围绕 $\sigma(A|_{\mathcal{M}})$ 但不包含 $\sigma(A|_{V/\mathcal{M}})$ 的简单闭曲线。定义谱投影
    $$P = \frac{1}{2\pi i} \oint_\Gamma (zI - A)^{-1} dz.$$
    则 $P^2 = P$，$AP = PA$，$\operatorname{Im}(P) = \mathcal{M}$，且 $\mathcal{N} = \ker(P)$ 是唯一的互补不变子空间。$\blacksquare$

!!! example "例 42.4"
    设 $A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 2 \end{pmatrix}$。

    - $\mathcal{M}_1 = \operatorname{span}\{e_1\}$（特征值 $1$）和 $\mathcal{M}_2 = \operatorname{span}\{e_2, e_3\}$（特征值 $2$）的谱不相交，因此 $\mathbb{C}^3 = \mathcal{M}_1 \oplus \mathcal{M}_2$ 是互补不变分解。

    - 但 $\mathcal{M}_3 = \operatorname{span}\{e_2\}$（$\ker(A - 2I)$ 中的子空间）是 $A$-不变的，其补空间中的谱与 $\mathcal{M}_3$ 上的谱重叠（都包含特征值 $2$）。可以验证不存在一维的互补不变子空间。

---

## 42.4 不变子空间与 Jordan 形

<div class="context-flow" markdown>

**核心问题**：Jordan 标准形如何完整描述所有不变子空间？

</div>

!!! theorem "定理 42.6 (Jordan 链与不变子空间)"
    设 $A \in \mathbb{C}^{n \times n}$ 的 Jordan 标准形中，特征值 $\lambda$ 对应的 Jordan 块为 $J_{n_1}(\lambda), J_{n_2}(\lambda), \ldots, J_{n_r}(\lambda)$，$n_1 \geq n_2 \geq \cdots \geq n_r$。设第 $j$ 个 Jordan 块对应的 Jordan 链为 $\{v_1^{(j)}, v_2^{(j)}, \ldots, v_{n_j}^{(j)}\}$，满足
    $$(A - \lambda I)v_k^{(j)} = v_{k-1}^{(j)}, \quad v_0^{(j)} = 0.$$

    则：

    1. 每条 Jordan 链的前 $k$ 个向量 $\operatorname{span}\{v_1^{(j)}, \ldots, v_k^{(j)}\}$（$1 \leq k \leq n_j$）张成 $A$-不变子空间。
    2. **谱不变子空间**：$\mathcal{V}_\lambda = \ker(A - \lambda I)^{n_1}$ 是 $A$-不变的，且 $\sigma(A|_{\mathcal{V}_\lambda}) = \{\lambda\}$。
    3. $\mathbb{C}^n = \bigoplus_{\lambda \in \sigma(A)} \mathcal{V}_\lambda$ 是谱不变子空间的直和分解。

??? proof "证明"
    (1) 设 $\mathcal{M}_k^{(j)} = \operatorname{span}\{v_1^{(j)}, \ldots, v_k^{(j)}\}$。对任意 $v = \sum_{i=1}^{k} c_i v_i^{(j)}$，
    $$Av = \sum_{i=1}^{k} c_i Av_i^{(j)} = \sum_{i=1}^{k} c_i (\lambda v_i^{(j)} + v_{i-1}^{(j)}) = \lambda v + \sum_{i=1}^{k} c_i v_{i-1}^{(j)}.$$
    由于 $v_{i-1}^{(j)} \in \mathcal{M}_k^{(j)}$（对 $i \leq k$），可得 $Av \in \mathcal{M}_k^{(j)}$。

    (2) 和 (3) 直接由 Jordan 标准形的结构和第 12 章的结论得出。$\blacksquare$

!!! definition "定义 42.6 (根子空间)"
    设 $A \in \mathbb{C}^{n \times n}$，$\lambda \in \sigma(A)$。定义升链
    $$\{0\} \subseteq \ker(A - \lambda I) \subseteq \ker(A - \lambda I)^2 \subseteq \cdots$$
    该链在某个 $m \leq n$ 处稳定，即 $\ker(A - \lambda I)^m = \ker(A - \lambda I)^{m+1} = \cdots$。稳定后的空间 $\ker(A - \lambda I)^m$ 就是特征值 $\lambda$ 的**根子空间**（root subspace），也称为**广义特征空间**。

!!! example "例 42.5"
    考虑 $A = \operatorname{diag}(J_3(0), J_2(0), J_1(1))$，其中
    $$J_3(0) = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}, \quad J_2(0) = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \quad J_1(1) = (1).$$

    特征值 $\lambda = 0$ 的根子空间升链：

    - $\ker(A) = \operatorname{span}\{e_1, e_4, e_6\}$（注意 $e_6$ 不在 $\ker(A)$ 中，因为 $Ae_6 = e_6 \neq 0$）。

      更准确地：$\ker(A) = \operatorname{span}\{e_1, e_4\}$（只看 $\lambda = 0$ 的部分），$\dim = 2$。

    - $\ker(A^2)|_{\mathcal{V}_0} = \operatorname{span}\{e_1, e_2, e_4, e_5\}$，$\dim = 4$。

    - $\ker(A^3)|_{\mathcal{V}_0} = \operatorname{span}\{e_1, e_2, e_3, e_4, e_5\}$，$\dim = 5$。

    每层核空间的维数差 $2, 2, 1$ 给出 Weyr 特征（参见第 44 章）。

---

## 42.5 子空间之间的角度

<div class="context-flow" markdown>

**核心问题**：如何量化两个子空间之间的"距离"或"接近程度"？

</div>

!!! definition "定义 42.7 (典则角度)"
    设 $\mathcal{F}$ 和 $\mathcal{G}$ 是 $\mathbb{R}^n$（或 $\mathbb{C}^n$）的子空间，$\dim \mathcal{F} = p$，$\dim \mathcal{G} = q$，$p \leq q$。$\mathcal{F}$ 与 $\mathcal{G}$ 之间的**典则角度**（canonical angles，也称主角）$\theta_1, \theta_2, \ldots, \theta_p$（$0 \leq \theta_1 \leq \theta_2 \leq \cdots \leq \theta_p \leq \pi/2$）递归定义如下：

    $$\cos \theta_k = \max_{\substack{u \in \mathcal{F}, \|u\| = 1 \\ u \perp u_1, \ldots, u_{k-1}}} \max_{\substack{v \in \mathcal{G}, \|v\| = 1 \\ v \perp v_1, \ldots, v_{k-1}}} |u^* v|,$$

    其中 $u_k$ 和 $v_k$ 是达到上述最大值的向量。

!!! theorem "定理 42.7 (典则角度的 SVD 刻画)"
    设 $P_{\mathcal{F}}$ 和 $P_{\mathcal{G}}$ 分别是到 $\mathcal{F}$ 和 $\mathcal{G}$ 的正交投影。则：

    1. $P_{\mathcal{F}} P_{\mathcal{G}}$ 的非零奇异值恰好是 $\cos \theta_1 \geq \cos \theta_2 \geq \cdots \geq \cos \theta_p > 0$（若某些 $\theta_k = \pi/2$ 则对应奇异值为 $0$）。

    2. 等价地，设 $F$ 和 $G$ 分别是 $\mathcal{F}$ 和 $\mathcal{G}$ 的标准正交基矩阵（列向量为基），则 $F^* G$ 的奇异值给出 $\cos \theta_k$。

??? proof "证明"
    设 $F \in \mathbb{C}^{n \times p}$ 和 $G \in \mathbb{C}^{n \times q}$ 是 $\mathcal{F}$ 和 $\mathcal{G}$ 的标准正交基矩阵。则 $P_{\mathcal{F}} = FF^*$，$P_{\mathcal{G}} = GG^*$。

    矩阵 $F^*G \in \mathbb{C}^{p \times q}$ 的 SVD 为 $F^*G = U \Sigma V^*$。设 $\hat{F} = FU$，$\hat{G} = GV$，则 $\hat{F}^* \hat{G} = U^* F^* G V = \Sigma$。

    这意味着 $\hat{f}_k^* \hat{g}_j = \sigma_k \delta_{kj}$，其中 $\hat{f}_k$ 和 $\hat{g}_k$ 是 $\hat{F}$ 和 $\hat{G}$ 的列。由于 $\hat{F}$ 和 $\hat{G}$ 仍是 $\mathcal{F}$ 和 $\mathcal{G}$ 的标准正交基（酉变换保持正交性），比较递归定义即得 $\sigma_k = \cos \theta_k$。$\blacksquare$

!!! definition "定义 42.8 ($\sin \Theta$ 矩阵和 $\tan \Theta$ 矩阵)"
    基于典则角度 $\theta_1, \ldots, \theta_p$，定义对角矩阵：

    $$\sin \Theta(\mathcal{F}, \mathcal{G}) = \operatorname{diag}(\sin \theta_1, \sin \theta_2, \ldots, \sin \theta_p),$$
    $$\cos \Theta(\mathcal{F}, \mathcal{G}) = \operatorname{diag}(\cos \theta_1, \cos \theta_2, \ldots, \cos \theta_p),$$
    $$\tan \Theta(\mathcal{F}, \mathcal{G}) = \operatorname{diag}(\tan \theta_1, \tan \theta_2, \ldots, \tan \theta_p).$$

    当两个子空间维数相同（$p = q$）时，有：
    $$\|\sin \Theta(\mathcal{F}, \mathcal{G})\|_2 = \|P_{\mathcal{F}} - P_{\mathcal{G}}\|_2.$$

!!! theorem "定理 42.8 ($\sin \Theta$ 的几何意义)"
    设 $\dim \mathcal{F} = \dim \mathcal{G} = p$。则：

    $$\|P_{\mathcal{F}} - P_{\mathcal{G}}\|_2 = \sin \theta_p,$$
    $$\|P_{\mathcal{F}} - P_{\mathcal{G}}\|_F = \left(\sum_{k=1}^{p} \sin^2 \theta_k\right)^{1/2} = \|\sin \Theta\|_F.$$

??? proof "证明"
    利用 $P_{\mathcal{F}} - P_{\mathcal{G}}$ 的奇异值等于 $\sin \theta_k$（每个出现两次）的经典结论。

    具体地，取 $\mathcal{F}$ 和 $\mathcal{G}$ 的典则基 $\hat{F}$ 和 $\hat{G}$，构造酉矩阵 $Q = (\hat{F} \mid \hat{F}_\perp \mid \hat{G}_\perp)$（适当扩充），可以将 $P_{\mathcal{F}} - P_{\mathcal{G}}$ 化为显式分块形式，从中读出奇异值。$\blacksquare$

!!! example "例 42.6"
    设 $\mathcal{F} = \operatorname{span}\{e_1\}$ 和 $\mathcal{G} = \operatorname{span}\left\{\frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}\right\}$ 是 $\mathbb{R}^2$ 中的一维子空间。

    $F^*G = (1, 0) \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \frac{1}{\sqrt{2}}$。

    因此 $\cos \theta_1 = \frac{1}{\sqrt{2}}$，$\theta_1 = \frac{\pi}{4}$。这正是两条直线之间的夹角。

    $\|P_{\mathcal{F}} - P_{\mathcal{G}}\|_2 = \sin \frac{\pi}{4} = \frac{1}{\sqrt{2}} \approx 0.707$。

---

## 42.6 Davis-Kahan sin Θ 定理

<div class="context-flow" markdown>

**核心问题**：当 Hermite 矩阵被扰动时，其不变子空间会偏转多少角度？

</div>

这是矩阵扰动理论中最重要的定理之一。它将矩阵的扰动大小与不变子空间的偏转角度定量联系起来。

!!! theorem "定理 42.9 (Davis-Kahan sin Θ 定理)"
    设 $A, \tilde{A} \in \mathbb{C}^{n \times n}$ 是 Hermite 矩阵，$E = \tilde{A} - A$。设 $\mathcal{V}$ 是 $A$ 的某个不变子空间（对应谱的一部分 $\Sigma_1 \subseteq \sigma(A)$），$\tilde{\mathcal{V}}$ 是 $\tilde{A}$ 的对应不变子空间（对应 $\tilde{\Sigma}_1 \subseteq \sigma(\tilde{A})$）。设

    $$\delta = \min_{\lambda \in \Sigma_1, \mu \in \sigma(A) \setminus \Sigma_1} |\lambda - \mu|$$

    是谱间隙（spectral gap）。若 $\delta > 0$，则

    $$\|\sin \Theta(\mathcal{V}, \tilde{\mathcal{V}})\|_F \leq \frac{\|E\|_F}{\delta}, \qquad \|\sin \Theta(\mathcal{V}, \tilde{\mathcal{V}})\|_2 \leq \frac{\|E\|_2}{\delta}.$$

    更精确地，对于任意酉不变范数 $\|\cdot\|$：

    $$\|\sin \Theta(\mathcal{V}, \tilde{\mathcal{V}})\| \leq \frac{\|R\|}{\delta},$$

    其中 $R = \tilde{A}\hat{V} - \hat{V}(\hat{V}^* \tilde{A} \hat{V})$ 是**残差**，$\hat{V}$ 是 $\mathcal{V}$ 的标准正交基矩阵。

??? proof "证明"
    **关键思路**：利用 Sylvester 方程建立 $\sin \Theta$ 与扰动的关系。

    设 $U = (\hat{V} \mid \hat{V}_\perp)$ 是酉矩阵，其中 $\hat{V}$ 的列张成 $\mathcal{V}$，$\hat{V}_\perp$ 的列张成 $\mathcal{V}^\perp$。类似地，$\tilde{U} = (\hat{\tilde{V}} \mid \hat{\tilde{V}}_\perp)$。

    由于 $\mathcal{V}$ 是 $A$ 的不变子空间，
    $$U^* A U = \begin{pmatrix} A_1 & 0 \\ 0 & A_2 \end{pmatrix},$$
    其中 $\sigma(A_1) = \Sigma_1$，$\sigma(A_2) = \sigma(A) \setminus \Sigma_1$。

    设 $\tilde{V}_\perp^* \hat{V} = S$（这个矩阵与 $\sin \Theta$ 密切相关——实际上 $S$ 的奇异值就是 $\sin \theta_k$）。

    由 $\tilde{A} \hat{\tilde{V}} = \hat{\tilde{V}} \tilde{A}_1$ 可得 $\hat{\tilde{V}}_\perp^* \tilde{A} \hat{\tilde{V}} = 0$。利用 $\tilde{A} = A + E$，

    $$\hat{\tilde{V}}_\perp^* (A + E) \hat{\tilde{V}} = 0.$$

    将 $\hat{\tilde{V}}$ 分解到 $\mathcal{V}$ 和 $\mathcal{V}^\perp$ 上，设 $\hat{V}^* \hat{\tilde{V}} = C$（奇异值为 $\cos \theta_k$），$\hat{V}_\perp^* \hat{\tilde{V}} = S'$（奇异值为 $\sin \theta_k$）。

    经过代数运算，可以得到 Sylvester 方程
    $$A_2 S' - S' A_1 = -\hat{V}_\perp^* E \hat{\tilde{V}} + O(\|E\|^2).$$

    由 Sylvester 方程的解估计（利用谱间隙 $\delta$）：
    $$\|S'\| \leq \frac{\|\hat{V}_\perp^* E \hat{\tilde{V}}\|}{\delta} \leq \frac{\|E\|}{\delta}.$$

    由于 $S'$ 的奇异值是 $\sin \theta_k$，这就给出了所需的估计。

    严格的证明需要更仔细地处理高阶项，使用如下更精确的 Sylvester 方程：

    $$A_2 (\hat{V}_\perp^* \hat{\tilde{V}}) - (\hat{V}_\perp^* \hat{\tilde{V}}) \tilde{A}_1 = -\hat{V}_\perp^* E \hat{\tilde{V}}.$$

    此方程精确成立（无高阶项），因为 $A \hat{V}_\perp = \hat{V}_\perp A_2$ 且 $\tilde{A} \hat{\tilde{V}} = \hat{\tilde{V}} \tilde{A}_1$。由 $\sigma(A_2) \cap \sigma(\tilde{A}_1) = \emptyset$（当 $\|E\| < \delta/2$ 时成立），Sylvester 方程有唯一解，且
    $$\|\hat{V}_\perp^* \hat{\tilde{V}}\| \leq \frac{\|\hat{V}_\perp^* E \hat{\tilde{V}}\|}{\min_{\lambda \in \sigma(A_2), \mu \in \sigma(\tilde{A}_1)} |\lambda - \mu|} \leq \frac{\|E\|}{\delta - \|E\|}.$$

    当 $\|E\|$ 相对于 $\delta$ 足够小时，这给出 $\|\sin \Theta\| \lesssim \|E\|/\delta$。$\blacksquare$

!!! example "例 42.7 (PCA 扰动分析)"
    在主成分分析（PCA）中，总体协方差矩阵 $\Sigma$ 的前 $k$ 个特征向量张成的子空间 $\mathcal{V}_k$ 是我们感兴趣的对象。样本协方差矩阵 $\hat{\Sigma}$ 是 $\Sigma$ 的扰动估计。

    设 $\lambda_1 \geq \cdots \geq \lambda_n$ 是 $\Sigma$ 的特征值。前 $k$ 个特征向量的谱间隙为
    $$\delta = \lambda_k - \lambda_{k+1}.$$

    Davis-Kahan 定理给出：
    $$\|\sin \Theta(\mathcal{V}_k, \hat{\mathcal{V}}_k)\|_F \leq \frac{\|\hat{\Sigma} - \Sigma\|_F}{\lambda_k - \lambda_{k+1}}.$$

    **数值例子**：设 $n = 100$，$k = 3$，$\lambda_3 = 5$，$\lambda_4 = 1$（谱间隙 $\delta = 4$），$\|\hat{\Sigma} - \Sigma\|_F = 0.8$。则

    $$\|\sin \Theta\|_F \leq \frac{0.8}{4} = 0.2,$$

    即前 3 个主成分方向的偏转（以 Frobenius 范数度量的 $\sin \Theta$）不超过 $0.2$。

!!! example "例 42.8 (具体矩阵的 Davis-Kahan 估计)"
    设
    $$A = \begin{pmatrix} 5 & 0 \\ 0 & 1 \end{pmatrix}, \quad \tilde{A} = \begin{pmatrix} 5 & 0.3 \\ 0.3 & 1 \end{pmatrix}.$$

    $E = \tilde{A} - A = \begin{pmatrix} 0 & 0.3 \\ 0.3 & 0 \end{pmatrix}$，$\|E\|_2 = 0.3$。

    谱间隙 $\delta = 5 - 1 = 4$。$A$ 的特征向量 $v_1 = e_1$ 张成 $\mathcal{V} = \operatorname{span}\{e_1\}$。

    $\tilde{A}$ 的特征值为 $\frac{6 \pm \sqrt{16.36}}{2} \approx 5.0225, 0.9775$，对应特征向量 $\tilde{v}_1 \approx \begin{pmatrix} 0.9978 \\ 0.0665 \end{pmatrix}$。

    实际 $\sin \theta_1 = |e_1^* \tilde{v}_{1,\perp}| \approx 0.0665$。

    Davis-Kahan 上界为 $\frac{0.3}{4} = 0.075$。

    估计相当紧：$0.0665 \leq 0.075$。

---

## 42.7 谱投影的扰动

<div class="context-flow" markdown>

**核心问题**：当矩阵被扰动时，谱投影（由环路积分定义的投影算子）如何变化？

</div>

!!! definition "定义 42.9 (谱投影)"
    设 $A \in \mathbb{C}^{n \times n}$，$\Gamma$ 是复平面中的简单闭曲线（正向），使得 $\sigma(A)$ 的一部分 $\Sigma_1$ 在 $\Gamma$ 内部，其余部分 $\Sigma_2 = \sigma(A) \setminus \Sigma_1$ 在 $\Gamma$ 外部。则

    $$P = \frac{1}{2\pi i} \oint_\Gamma (zI - A)^{-1} \, dz$$

    称为 $A$ 关于 $\Sigma_1$ 的**谱投影**（spectral projection）。

!!! theorem "定理 42.10 (谱投影的基本性质)"
    上述谱投影 $P$ 满足：

    1. $P^2 = P$（$P$ 是幂等的，即投影算子）。
    2. $AP = PA$（$P$ 与 $A$ 交换）。
    3. $\operatorname{Im}(P)$ 是 $A$ 的不变子空间，$\sigma(A|_{\operatorname{Im}(P)}) = \Sigma_1$。
    4. $\ker(P)$ 也是 $A$ 的不变子空间，$\sigma(A|_{\ker(P)}) = \Sigma_2$。
    5. 若 $A$ 是 Hermite 矩阵，则 $P$ 是正交投影（$P = P^*$）。

??? proof "证明"
    **(1)** 设 $\Gamma_1$ 和 $\Gamma_2$ 是两条嵌套的等价曲线（$\Gamma_1$ 在 $\Gamma_2$ 内部，但都包围 $\Sigma_1$）。利用预解式恒等式
    $$(zI - A)^{-1} - (wI - A)^{-1} = (w - z)(zI - A)^{-1}(wI - A)^{-1},$$
    计算
    \begin{align}
    P^2 &= \frac{1}{(2\pi i)^2} \oint_{\Gamma_1} \oint_{\Gamma_2} (zI - A)^{-1}(wI - A)^{-1} \, dw \, dz \\
    &= \frac{1}{(2\pi i)^2} \oint_{\Gamma_1} \oint_{\Gamma_2} \frac{(zI - A)^{-1} - (wI - A)^{-1}}{w - z} \, dw \, dz.
    \end{align}

    对 $w$ 积分时，$(wI - A)^{-1}$ 在 $\Gamma_2$ 内部关于 $w$ 全纯（因为 $z$ 在 $\Gamma_2$ 内部，是唯一的极点），利用留数定理得到 $P^2 = P$。

    **(2)** $A$ 与 $(zI - A)^{-1}$ 交换，因此与积分交换。

    **(5)** 当 $A = A^*$ 时，$(zI - A)^{-1}$ 的共轭满足 $\overline{(zI - A)^{-1}} = (\bar{z}I - A)^{-1}$，利用谱在实轴上的对称性可得 $P^* = P$。$\blacksquare$

!!! theorem "定理 42.11 (Kato 扰动定理)"
    设 $A(t) = A + tE$ 是矩阵的解析扰动（$t \in \mathbb{C}$，$|t|$ 足够小）。设 $\Gamma$ 是围绕 $\sigma(A)$ 的一部分 $\Sigma_1$ 的简单闭曲线，且 $\operatorname{dist}(\Gamma, \sigma(A)) = d > 0$。

    当 $|t| \cdot \|E\| < d$ 时，$A(t)$ 的谱在 $\Gamma$ 内部的部分 $\Sigma_1(t)$ 与 $\Sigma_1$ 具有相同数量的特征值（计重数），且谱投影

    $$P(t) = \frac{1}{2\pi i} \oint_\Gamma (zI - A(t))^{-1} \, dz$$

    是 $t$ 的解析函数。特别地：

    $$P(t) = P + t P^{(1)} + t^2 P^{(2)} + \cdots,$$

    其中

    $$P^{(1)} = -\frac{1}{2\pi i} \oint_\Gamma (zI - A)^{-1} E (zI - A)^{-1} \, dz.$$

??? proof "证明"
    关键是预解式的 Neumann 级数展开。当 $|t| \cdot \|E\| < d$ 时，对 $z \in \Gamma$：
    $$(zI - A(t))^{-1} = (zI - A - tE)^{-1} = \left[(zI - A)(I - t(zI - A)^{-1}E)\right]^{-1}.$$

    由于 $\|t(zI - A)^{-1}E\| \leq |t| \cdot \|E\| / d < 1$，Neumann 级数
    $$(zI - A(t))^{-1} = \sum_{k=0}^{\infty} t^k \left[(zI - A)^{-1} E\right]^k (zI - A)^{-1}$$
    一致收敛。逐项积分即得 $P(t)$ 的级数展开。$\blacksquare$

!!! example "例 42.9 (一阶谱投影扰动)"
    设 $A = \operatorname{diag}(\lambda_1, \lambda_2)$（$\lambda_1 \neq \lambda_2$），$E = \begin{pmatrix} e_{11} & e_{12} \\ e_{21} & e_{22} \end{pmatrix}$。

    $P = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$ 是 $A$ 关于 $\{\lambda_1\}$ 的谱投影。

    一阶扰动：取 $\Gamma$ 围绕 $\lambda_1$，
    $$P^{(1)} = -\frac{1}{2\pi i} \oint_\Gamma \begin{pmatrix} (z - \lambda_1)^{-1} & 0 \\ 0 & (z - \lambda_2)^{-1} \end{pmatrix} E \begin{pmatrix} (z - \lambda_1)^{-1} & 0 \\ 0 & (z - \lambda_2)^{-1} \end{pmatrix} dz.$$

    计算留数（$z = \lambda_1$ 处）：
    $$P^{(1)} = \begin{pmatrix} 0 & \frac{e_{12}}{\lambda_1 - \lambda_2} \\[4pt] \frac{e_{21}}{\lambda_2 - \lambda_1} & 0 \end{pmatrix} = \begin{pmatrix} 0 & \frac{e_{12}}{\lambda_1 - \lambda_2} \\[4pt] -\frac{e_{21}}{\lambda_1 - \lambda_2} & 0 \end{pmatrix}.$$

    可以验证：$P^{(1)}P + PP^{(1)} = P^{(1)}$（由 $P(t)^2 = P(t)$ 的一阶条件），且

    $$\|P^{(1)}\|_2 = \frac{\max(|e_{12}|, |e_{21}|)}{|\lambda_1 - \lambda_2|},$$

    谱间隙 $|\lambda_1 - \lambda_2|$ 出现在分母中，与 Davis-Kahan 定理一致。

!!! theorem "定理 42.12 (谱投影距离与 sin Θ 的关系)"
    设 $A$ 和 $\tilde{A}$ 是 Hermite 矩阵，$P$ 和 $\tilde{P}$ 分别是关于谱的对应部分的谱投影（正交投影），$\dim \operatorname{Im}(P) = \dim \operatorname{Im}(\tilde{P}) = k$。则：

    $$\|P - \tilde{P}\|_2 = \|\sin \Theta(\operatorname{Im}(P), \operatorname{Im}(\tilde{P}))\|_2 = \sin \theta_{\max},$$

    $$\|P - \tilde{P}\|_F = \sqrt{2} \|\sin \Theta(\operatorname{Im}(P), \operatorname{Im}(\tilde{P}))\|_F.$$

    因此 Davis-Kahan 定理等价于谱投影的扰动界。

??? proof "证明"
    $P - \tilde{P}$ 是 Hermite 矩阵（因为 $P$ 和 $\tilde{P}$ 都是 Hermite 的）。通过在适当的典则基下将 $P$ 和 $\tilde{P}$ 同时分块对角化，可以证明 $P - \tilde{P}$ 的特征值恰好是 $\pm \sin \theta_k$（每个出现一次）和若干个 $0$。

    具体地，在典则基下：
    $$P = \begin{pmatrix} I & 0 \\ 0 & 0 \end{pmatrix}, \quad \tilde{P} = \begin{pmatrix} C^2 & CS \\ CS & S^2 \end{pmatrix},$$
    其中 $C = \cos \Theta$，$S = \sin \Theta$ 是对角矩阵。因此
    $$P - \tilde{P} = \begin{pmatrix} S^2 & -CS \\ -CS & -S^2 \end{pmatrix},$$
    其特征值为 $\pm \sin \theta_k$。从而 $\|P - \tilde{P}\|_2 = \max_k |\sin \theta_k| = \sin \theta_{\max}$，$\|P - \tilde{P}\|_F = \sqrt{2 \sum_k \sin^2 \theta_k} = \sqrt{2} \|\sin \Theta\|_F$。$\blacksquare$

!!! example "例 42.10 (综合应用)"
    考虑对称矩阵
    $$A = \begin{pmatrix} 10 & 0 & 0 \\ 0 & 5 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$
    及其扰动 $\tilde{A} = A + \epsilon E$，其中 $E$ 是随机对称矩阵，$\|E\|_2 = 1$，$\epsilon = 0.1$。

    - 关于特征值 $\{10\}$（$\mathcal{V} = \operatorname{span}\{e_1\}$），谱间隙 $\delta_1 = 10 - 5 = 5$。
      Davis-Kahan 估计：$\sin \theta_1 \leq \epsilon / \delta_1 = 0.1 / 5 = 0.02$。

    - 关于特征值 $\{10, 5\}$（$\mathcal{V} = \operatorname{span}\{e_1, e_2\}$），谱间隙 $\delta_2 = 5 - 1 = 4$。
      Davis-Kahan 估计：$\|\sin \Theta\|_2 \leq 0.1 / 4 = 0.025$。

    这些估计表明，谱间隙越大，不变子空间对扰动越稳定。

---

## 42.8 Wedin sin Θ 定理

<div class="context-flow" markdown>

**核心问题**：Davis-Kahan 定理适用于 Hermite 矩阵的特征子空间扰动。对于一般矩阵的**奇异子空间**（左奇异向量和右奇异向量张成的子空间），扰动界如何建立？

</div>

Wedin sin Θ 定理是 Davis-Kahan 定理在 SVD 框架下的推广，由 Per-Åke Wedin 于 1972 年建立。它在 PCA 扰动分析、低秩矩阵逼近、信号处理等领域具有核心地位。

!!! theorem "定理 42.13 (Wedin sin Θ 定理)"
    设 $A, \tilde{A} \in \mathbb{C}^{m \times n}$，$E = \tilde{A} - A$。设 $A$ 的 SVD 中，奇异值分为两组：$\Sigma_1 = \operatorname{diag}(\sigma_1, \ldots, \sigma_k)$（对应左奇异向量矩阵 $U_1$ 和右奇异向量矩阵 $V_1$）及其余 $\Sigma_2$。类似地，$\tilde{A}$ 的对应分组为 $\tilde{\Sigma}_1, \tilde{U}_1, \tilde{V}_1$。

    设奇异值间隙
    $$\delta = \min\bigl(\min_{i \leq k, j > k} |\sigma_i - \tilde{\sigma}_j|,\; \min_{i \leq k} \sigma_i\bigr) > 0.$$

    更标准的表述：设 $\alpha = \min_{i \leq k} \sigma_i$ 且 $\tilde{\beta} = \max_{j > k} \tilde{\sigma}_j$，若 $\alpha > \tilde{\beta}$，则

    $$\max\!\bigl(\|\sin \Theta(U_1, \tilde{U}_1)\|_F,\; \|\sin \Theta(V_1, \tilde{V}_1)\|_F\bigr) \leq \frac{\max(\|EU_1\|_F, \|E^* V_1\|_F)}{\alpha - \tilde{\beta}}.$$

    对谱范数同样成立：

    $$\max\!\bigl(\|\sin \Theta(U_1, \tilde{U}_1)\|_2,\; \|\sin \Theta(V_1, \tilde{V}_1)\|_2\bigr) \leq \frac{\max(\|E\|_2, \|E\|_2)}{\alpha - \tilde{\beta}} = \frac{\|E\|_2}{\alpha - \tilde{\beta}}.$$

??? proof "证明"
    **CS 分解方法。** 设 $A$ 的紧 SVD 分组为

    $$A = (U_1 \mid U_2)\begin{pmatrix} \Sigma_1 & 0 \\ 0 & \Sigma_2 \end{pmatrix}(V_1 \mid V_2)^*, \quad \tilde{A} = (\tilde{U}_1 \mid \tilde{U}_2)\begin{pmatrix} \tilde{\Sigma}_1 & 0 \\ 0 & \tilde{\Sigma}_2 \end{pmatrix}(\tilde{V}_1 \mid \tilde{V}_2)^*.$$

    由 $\tilde{A} \tilde{V}_1 = \tilde{U}_1 \tilde{\Sigma}_1$，左乘 $U_2^*$ 得

    $$U_2^* \tilde{A} \tilde{V}_1 = U_2^* \tilde{U}_1 \tilde{\Sigma}_1.$$

    另一方面，$A^* \tilde{U}_1 = V \Sigma^* U^* \tilde{U}_1$，因此

    $$\tilde{A}^* \tilde{U}_1 = \tilde{V}_1 \tilde{\Sigma}_1, \quad V_2^* \tilde{A}^* \tilde{U}_1 = V_2^* \tilde{V}_1 \tilde{\Sigma}_1.$$

    由 $\tilde{A} = A + E$，展开 $U_2^* \tilde{A} \tilde{V}_1$：

    $$U_2^* (A + E) \tilde{V}_1 = U_2^* A \tilde{V}_1 + U_2^* E \tilde{V}_1.$$

    注意 $U_2^* A = \Sigma_2 V_2^*$（由 SVD），因此

    $$\Sigma_2 V_2^* \tilde{V}_1 + U_2^* E \tilde{V}_1 = U_2^* \tilde{U}_1 \tilde{\Sigma}_1.$$

    这给出 Sylvester 型方程：

    $$(U_2^* \tilde{U}_1) \tilde{\Sigma}_1 - \Sigma_2 (V_2^* \tilde{V}_1) = U_2^* E \tilde{V}_1.$$

    类似地，从 $\tilde{A}^* = A^* + E^*$ 出发可得

    $$(V_2^* \tilde{V}_1) \tilde{\Sigma}_1 - \Sigma_2^* (U_2^* \tilde{U}_1) = V_2^* E^* \tilde{U}_1.$$

    设 $P = U_2^* \tilde{U}_1$，$Q = V_2^* \tilde{V}_1$。$P$ 的奇异值是 $\sin \theta_j(U_1, \tilde{U}_1)$，$Q$ 的奇异值是 $\sin \theta_j(V_1, \tilde{V}_1)$。

    将两个方程合并写成分块 Sylvester 方程后，利用 $\alpha > \tilde{\beta}$ 条件（确保 $\tilde{\Sigma}_1$ 的奇异值与 $\Sigma_2$ 的奇异值分离），由 Sylvester 方程解的范数估计得

    $$\max(\|P\|, \|Q\|) \leq \frac{\max(\|U_2^* E \tilde{V}_1\|, \|V_2^* E^* \tilde{U}_1\|)}{\alpha - \tilde{\beta}} \leq \frac{\max(\|E\|, \|E\|)}{\alpha - \tilde{\beta}}.$$

    由此得到 $\sin \Theta$ 界。$\blacksquare$

!!! example "例 42.11 (PCA 扰动分析——Wedin 定理的应用)"
    在主成分分析中，数据矩阵 $X \in \mathbb{R}^{n \times p}$（$n$ 个样本，$p$ 个特征）的前 $k$ 个右奇异向量张成**主成分子空间** $\mathcal{V}_k$。在统计应用中，$X = X_0 + E$，其中 $X_0$ 是信号部分，$E$ 是噪声。

    设 $X_0$ 的前 $k$ 个奇异值为 $\sigma_1 \geq \cdots \geq \sigma_k$，第 $k+1$ 个奇异值为 $\sigma_{k+1}$（若 $X_0$ 是秩 $k$ 的则 $\sigma_{k+1} = 0$）。由 Weyl 不等式，$\tilde{\sigma}_{k+1} \leq \sigma_{k+1} + \|E\|_2$。

    **信噪比条件**：若 $\sigma_k \gg \|E\|_2 + \sigma_{k+1}$，则奇异值间隙 $\alpha - \tilde{\beta} \approx \sigma_k - \sigma_{k+1} - \|E\|_2$，Wedin 定理给出

    $$\|\sin \Theta(\mathcal{V}_k, \tilde{\mathcal{V}}_k)\|_F \leq \frac{\|E\|_F}{\sigma_k - \sigma_{k+1} - \|E\|_2}.$$

    **数值例子**：设 $n = 500$，$p = 200$，$k = 5$，$\sigma_5 = 20$，$\sigma_6 = 0$（秩 5 信号），$\|E\|_2 = 5$，$\|E\|_F = 30$。则

    $$\|\sin \Theta\|_F \leq \frac{30}{20 - 0 - 5} = 2.0.$$

    由于 $\|\sin \Theta\|_F \leq \sqrt{k} = \sqrt{5} \approx 2.24$（自然上界），此界是有信息量的。若信号更强（$\sigma_5 = 50$），则界改善为 $30/45 \approx 0.67$。

    此分析在高维统计中极为重要：它告诉我们 PCA 所估计的子空间与真实信号子空间之间的偏差，直接受信噪比（$\sigma_k / \|E\|$）控制。

---

## 42.9 Stewart 的 tan Θ 定理

<div class="context-flow" markdown>

**核心问题**：能否得到比 sin Θ 更紧的子空间扰动界？

</div>

G. W. Stewart 的 tan Θ 定理在某些情形下给出比 sin Θ 定理更强的结果，因为 $\tan \theta \geq \sin \theta$（当 $\theta \in [0, \pi/2)$ 时）意味着 tan Θ 界控制了 sin Θ 界。

!!! theorem "定理 42.14 (Stewart tan Θ 定理)"
    设 $A \in \mathbb{C}^{n \times n}$ 是 Hermite 矩阵，其不变子空间分解为

    $$U^* A U = \begin{pmatrix} A_1 & 0 \\ 0 & A_2 \end{pmatrix},$$

    其中 $U = (\hat{V} \mid \hat{V}_\perp)$ 是酉矩阵。设 $R \in \mathbb{C}^{(n-k) \times k}$ 满足 $\hat{V}_\perp^* \tilde{A} \hat{V} = R$（残差矩阵），且谱分离条件 $\delta = \inf\{|\lambda - \mu| : \lambda \in \sigma(A_1), \mu \in \sigma(A_2)\} > 0$ 成立。

    若 $\tilde{\mathcal{V}}$ 是 $\tilde{A} = A + E$ 的对应不变子空间，则

    $$\|\tan \Theta(\mathcal{V}, \tilde{\mathcal{V}})\| \leq \frac{\|R\|}{\delta}$$

    对任意酉不变范数成立，其中 $R = \hat{V}_\perp^* \tilde{A} \hat{V}$ 是残差。

??? proof "证明"
    设 $\tilde{\mathcal{V}}$ 的标准正交基为 $\hat{\tilde{V}}$。将 $\hat{\tilde{V}}$ 分解在 $\mathcal{V}$ 和 $\mathcal{V}^\perp$ 上：

    $$\hat{\tilde{V}} = \hat{V} C + \hat{V}_\perp S',$$

    其中 $C = \hat{V}^* \hat{\tilde{V}}$，$S' = \hat{V}_\perp^* \hat{\tilde{V}}$。$C$ 的奇异值是 $\cos \theta_j$，$S'$ 的奇异值是 $\sin \theta_j$。

    由 $\tilde{A} \hat{\tilde{V}} = \hat{\tilde{V}} \tilde{A}_1$，左乘 $\hat{V}_\perp^*$：

    $$\hat{V}_\perp^* \tilde{A} \hat{\tilde{V}} = S' \tilde{A}_1,$$

    即 $\hat{V}_\perp^* \tilde{A} (\hat{V} C + \hat{V}_\perp S') = S' \tilde{A}_1$。

    展开并利用 $A_2 S' - S' \tilde{A}_1 = -\hat{V}_\perp^* E \hat{\tilde{V}}$ 的 Sylvester 方程结构，结合 $R = \hat{V}_\perp^* \tilde{A} \hat{V}$，可得

    $$S' = RC^{-1} \cdot (\text{bounded factor involving Sylvester equation}).$$

    更直接地，$\tan \Theta$ 矩阵的奇异值为 $\sin \theta_j / \cos \theta_j$，即 $S' C^{-1}$ 的奇异值（当 $\cos \theta_j > 0$）。由精确的 Sylvester 方程

    $$A_2 (S' C^{-1}) - (S' C^{-1}) \tilde{A}_1 = -R,$$

    其解范数受控于 $\|R\| / \delta$，从而 $\|\tan \Theta\| \leq \|R\| / \delta$。$\blacksquare$

tan Θ 定理的优势在于：$\|R\|$ 是**残差**的大小，而非整个扰动 $\|E\|$。当 $\tilde{\mathcal{V}}$ 是 $\tilde{A}$ 的精确不变子空间时，$R = 0$，界为零——这是自然的。在迭代算法（如子空间迭代法）中，残差通常比扰动更小，因此 tan Θ 定理给出更紧的收敛保证。

---

## 42.10 Rosenblum 定理与 Sylvester 方程

<div class="context-flow" markdown>

**核心问题**：Davis-Kahan 和 Wedin 定理的证明中反复出现 Sylvester 方程 $AX - XB = C$。该方程的可解性条件是什么？

</div>

!!! theorem "定理 42.15 (Rosenblum 定理)"
    设 $A \in \mathbb{C}^{m \times m}$，$B \in \mathbb{C}^{n \times n}$。Sylvester 方程

    $$AX - XB = C$$

    对任意 $C \in \mathbb{C}^{m \times n}$ 有唯一解 $X$，当且仅当 $\sigma(A) \cap \sigma(B) = \emptyset$（即 $A$ 与 $B$ 没有公共特征值）。

??? proof "证明"
    **必要性**：若 $\sigma(A) \cap \sigma(B) \neq \emptyset$，设 $\lambda$ 是公共特征值。存在 $u \neq 0$ 使得 $Au = \lambda u$，$v \neq 0$ 使得 $Bv = \lambda v$（或 $B^* \bar{v} = \bar{\lambda} \bar{v}$）。

    考虑齐次方程 $AX - XB = 0$。取 $X_0 = u v^*$，则

    $$AX_0 - X_0 B = Au v^* - u v^* B = \lambda u v^* - u (\lambda v^*) = 0.$$

    因此 $X_0 \neq 0$ 是齐次方程的非平凡解，方程 $AX - XB = C$ 的解不唯一。

    **充分性**：定义线性算子 $\mathcal{T}: \mathbb{C}^{m \times n} \to \mathbb{C}^{m \times n}$，$\mathcal{T}(X) = AX - XB$。$\mathcal{T}$ 的特征值恰好是 $\{\lambda_i - \mu_j : \lambda_i \in \sigma(A), \mu_j \in \sigma(B)\}$。

    当 $\sigma(A) \cap \sigma(B) = \emptyset$ 时，$0 \notin \sigma(\mathcal{T})$，因此 $\mathcal{T}$ 可逆，方程对任意 $C$ 有唯一解。

    解的显式表示可通过环路积分给出：

    $$X = \frac{1}{2\pi i} \oint_\Gamma (zI - A)^{-1} C (zI - B)^{-1} dz,$$

    其中 $\Gamma$ 是围绕 $\sigma(A)$ 但不包含 $\sigma(B)$ 的简单闭曲线。$\blacksquare$

Rosenblum 定理给出了不变子空间扰动理论中的核心代数工具。Davis-Kahan 定理证明中出现的方程 $A_2 S - S \tilde{A}_1 = -\hat{V}_\perp^* E \hat{\tilde{V}}$ 正是 Sylvester 方程，其可解性由谱间隙 $\delta > 0$ 保证。解的范数估计 $\|X\| \leq \|C\| / \delta$ 直接给出了 sin Θ 界。

---

## 42.11 无穷维不变子空间问题

**不变子空间问题**（invariant subspace problem）是泛函分析中最著名的未解决问题之一：

> 设 $\mathcal{H}$ 是可分的无穷维复 Hilbert 空间，$T: \mathcal{H} \to \mathcal{H}$ 是有界线性算子。$T$ 是否一定存在非平凡的闭不变子空间？

这里"非平凡"指既不是 $\{0\}$ 也不是 $\mathcal{H}$ 的闭子空间。

在有限维情形中，答案是肯定的——每个特征值对应的特征空间就是非平凡不变子空间（因为代数闭域上的线性变换总有特征值）。但在无穷维情形中，算子可能没有特征值，问题变得深刻得多。

已知的部分结果包括：

- **紧算子**：Aronszajn-Smith 定理（1954）证明了紧算子在可分 Hilbert 空间上总有非平凡闭不变子空间。
- **正规算子**：由谱定理，正规算子的谱投影提供了丰富的不变子空间。
- **反例在 Banach 空间中存在**：Enflo（1976，发表于 1987）和 Read（1984）分别构造了在某些 Banach 空间上没有非平凡闭不变子空间的有界线性算子。Read 甚至在 $\ell^1$ 上给出了反例。
- 然而，对于 **Hilbert 空间**上的一般有界算子，问题仍然开放。

此问题与算子代数、复分析（解析函数空间上的算子）、遍历理论等领域有深刻的联系。

---

## 练习题

1. **[基础] 证明：线性变换 $T$ 的核 $\ker(T)$ 和像 $\operatorname{Im}(T)$ 总是 $T$-不变子空间。**
   ??? success "参考答案"
       - 若 $v \in \ker(T)$，则 $Tv = 0$。由于 $0$ 属于任何子空间，故 $0 \in \ker(T)$，即 $Tv \in \ker(T)$。
       - 若 $v \in \operatorname{Im}(T)$，则 $Tv$ 显然属于 $T$ 的像集，即 $Tv \in \operatorname{Im}(T)$。

2. **[特征向量] 证明：任何 1 维不变子空间必由一个特征向量张成。**
   ??? success "参考答案"
       设 $\mathcal{M} = \operatorname{span}\{v\}$（$v \ne 0$）。若 $\mathcal{M}$ 是不变的，则 $Tv \in \mathcal{M}$，即 $Tv = \lambda v$ 对某个标量 $\lambda$ 成立。这正是特征向量的定义。

3. **[分块形式] 若 $A$ 有一个 $k$ 维不变子空间，证明 $A$ 相似于一个分块上三角矩阵。**
   ??? success "参考答案"
       取该子空间的一组基并扩充为全空间的基。在此基下，前 $k$ 列在第 $k$ 行以下全为零（因为变换结果仍落在前 $k$ 个基向量的张成空间内），故呈现 $\begin{pmatrix} B & C \\ 0 & D \end{pmatrix}$ 的形式。

4. **[超不变] 定义超不变子空间。是否每个不变子空间都是超不变的？**
   ??? success "参考答案"
       超不变子空间是指在所有与 $T$ 交换的算子作用下都保持不变的子空间。并非所有不变子空间都是超不变的。例如当 $T=I$ 时，所有子空间都是不变的，但只有 $\{0\}$ 和 $V$ 是超不变的。

5. **[约化子空间] 证明：$\mathcal{M}$ 约化 $T$ 当且仅当 $\mathcal{M}$ 和 $\mathcal{M}^\perp$ 都是 $T$-不变的。**
   ??? success "参考答案"
       这是约化子空间的定义。在投影算子语言中，这等价于 $T$ 与到 $\mathcal{M}$ 的正交投影算子交换。

6. **[典则角度] 求 $\mathbb{R}^2$ 中 $\mathcal{F} = \operatorname{span}\{e_1\}$ 与 $\mathcal{G} = \operatorname{span}\{(1, 1)^T\}$ 之间的典则角度。**
   ??? success "参考答案"
       角度的余弦值等于正交投影乘积的奇异值。$\cos \theta = \frac{|(1,0) \cdot (1,1)|}{1 \cdot \sqrt{2}} = \frac{1}{\sqrt{2}}$。故 $\theta = \pi/4$（即 45 度）。

7. **[sin Θ 界] 利用 Davis-Kahan 定理，估计 $A = \operatorname{diag}(10, 1)$ 的主特征向量在受到扰动 $E = \begin{pmatrix} 0 & 0.1 \\ 0.1 & 0 \end{pmatrix}$ 后的偏转。**
   ??? success "参考答案"
       谱间隙 $\delta = 10 - 1 = 9$。扰动范数 $\|E\|_2 = 0.1$。根据定理，$\sin \theta \le 0.1 / 9 \approx 0.011$。

8. **[谱投影] 描述与孤立特征值 $\lambda$ 相关的谱投影。**
   ??? success "参考答案"
       它由围道积分 $P = \frac{1}{2\pi i} \oint_\Gamma (zI - A)^{-1} dz$ 给出，其中 $\Gamma$ 仅包围 $\lambda$。对于可对角化矩阵，它等于到该特征值对应特征空间的投影。

9. **[Wedin定理] 何时使用 Wedin sin Θ 定理而不是 Davis-Kahan 定理？**
   ??? success "参考答案"
       Wedin 定理用于研究**奇异值分解（SVD）**中的奇异子空间扰动，而 Davis-Kahan 定理专门用于 **Hermite 矩阵**的特征子空间扰动。

10. **[Rosenblum定理] Sylvester 方程 $AX - XB = C$ 何时有唯一解？**
    ??? success "参考答案"
        当且仅当 $A$ 和 $B$ 没有公共特征值时，即 $\sigma(A) \cap \sigma(B) = \emptyset$。

## 本章小结

本章探讨了线性算子的稳定性与结构：

1. **子空间动力学**：定义了不变、超不变和约化子空间，作为算子分解的代数单元。
2. **几何度量**：引入了典则角度，定量化了线性子空间之间的“距离”。
3. **扰动稳健性**：解读了 Davis-Kahan 和 Wedin 定理，将谱间隙与特征向量及奇异向量的稳定性相联系。
4. **分析工具**：利用谱投影和 Rosenblum 定理，为子空间分析提供了严谨的数学框架。

