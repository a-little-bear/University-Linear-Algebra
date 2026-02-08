# 第 54 章 Quiver 表示

<div class="context-flow" markdown>

**前置**：向量空间 (Ch4) · 线性变换 (Ch5) · Jordan 形 (Ch12)

**本章脉络**：Quiver 定义 → 表示 → 态射与直和 → 不可分解表示 → Gabriel 定理（Dynkin 图 $\leftrightarrow$ 有限型） → Auslander-Reiten 理论初步 → 持久同调中的应用

**延伸**：Quiver 表示论统一了大量矩阵分类问题；Gabriel 定理揭示了 ADE 分类（与单 Lie 代数、有限子群、奇点分类的惊人对应）；持久模是拓扑数据分析（TDA）的代数基础

</div>

线性代数中的许多问题——单个矩阵的相似分类、矩阵对的同时分类、矩阵束的分析——表面上各不相同，但都可以统一到一个优雅的框架中：**Quiver 表示**。一个 quiver（箭图）就是一个有向图；它的表示是在每个顶点上放一个向量空间，在每条箭头上放一个线性映射。这个看似简单的定义蕴含了极其丰富的结构：Gabriel 于 1972 年证明了一个深刻的定理——quiver 具有有限多个不可分解表示，当且仅当其底图是 Dynkin 图（$A_n, D_n, E_6, E_7, E_8$）。这一结果将 quiver 表示论与 Lie 理论、代数几何、奇点理论等领域深刻地联系了起来。近年来，quiver 表示论还在拓扑数据分析（TDA）中找到了重要应用：持久同调的代数结构正是 $A_n$ 型 quiver 的表示。

---

## 54.1 Quiver 的定义

<div class="context-flow" markdown>

**核心问题**：如何用图论语言为矩阵分类问题提供统一框架？

</div>

!!! definition "定义 54.1 (Quiver)"
    一个 **quiver**（箭图）$Q = (Q_0, Q_1, s, t)$ 由以下数据组成：

    - $Q_0$：**顶点集**（vertex set），一个有限集。
    - $Q_1$：**箭头集**（arrow set），一个有限集。
    - $s: Q_1 \to Q_0$：**源映射**（source map），将每条箭头映到它的起点。
    - $t: Q_1 \to Q_0$：**靶映射**（target map），将每条箭头映到它的终点。

    对于箭头 $\alpha \in Q_1$，记 $\alpha: s(\alpha) \to t(\alpha)$。

!!! remark "注记"
    Quiver 允许自环（$s(\alpha) = t(\alpha)$）和重边（不同的 $\alpha, \beta$ 可以有 $s(\alpha) = s(\beta)$，$t(\alpha) = t(\beta)$）。Quiver 可以看作有限有向图的同义词，但在表示论中习惯使用 "quiver" 这一术语（源自 Gabriel 的法语 "carquois"，即箭袋）。

!!! example "例 54.1 (基本 quiver)"
    **(a) Jordan quiver** $L_1$：一个顶点、一条自环。
    $$\bullet \circlearrowleft$$

    **(b) Kronecker quiver** $K_2$：两个顶点、两条同向箭头。
    $$\bullet \rightrightarrows \bullet$$

    **(c) $A_n$ 型 quiver**（线性定向）：$n$ 个顶点排成一行，每对相邻顶点之间一条箭头。
    $$\bullet \to \bullet \to \cdots \to \bullet$$

    **(d) $D_4$ 型 quiver**（一种定向）：
    $$\begin{array}{ccc} & \bullet & \\ & \downarrow & \\ \bullet \to & \bullet & \leftarrow \bullet \end{array}$$

    **(e) 循环 quiver** $\tilde{A}_n$：$n+1$ 个顶点排成环。

!!! definition "定义 54.2 (Quiver 的底图)"
    Quiver $Q$ 的**底图**（underlying graph）$\bar{Q}$ 是忘掉箭头方向后得到的无向图。若 $\bar{Q}$ 连通，则称 $Q$ 是**连通**的。

!!! definition "定义 54.3 (路)"
    Quiver $Q$ 中的一条**路**（path）是箭头的有限序列 $\alpha_k \alpha_{k-1} \cdots \alpha_1$，其中 $t(\alpha_i) = s(\alpha_{i+1})$（$i = 1, \ldots, k-1$）。路的**长度**为 $k$。每个顶点 $i$ 对应一条长度为 $0$ 的**平凡路** $e_i$。

!!! definition "定义 54.4 (路代数)"
    Quiver $Q$ 的**路代数** $\mathbb{F}Q$ 是以 $Q$ 的所有路为基，乘法为路的复合（若首尾不相接则乘积为零）的结合代数。平凡路 $e_i$ 是（局部）单位元：$e_i \cdot \alpha = \alpha$ 若 $s(\alpha) = i$，否则为零。

---

## 54.2 Quiver 的表示

<div class="context-flow" markdown>

**核心问题**：如何在 quiver 上"放置"线性代数数据？

</div>

!!! definition "定义 54.5 (Quiver 的表示)"
    域 $\mathbb{F}$ 上 quiver $Q = (Q_0, Q_1, s, t)$ 的一个**表示** $V = (V_i, f_\alpha)$ 由以下数据组成：

    - 对每个顶点 $i \in Q_0$，一个有限维 $\mathbb{F}$-向量空间 $V_i$。
    - 对每条箭头 $\alpha: i \to j$，一个线性映射 $f_\alpha: V_i \to V_j$。

    表示 $V$ 的**维数向量**（dimension vector）为 $\mathbf{d} = (\dim V_i)_{i \in Q_0} \in \mathbb{Z}_{\geq 0}^{Q_0}$。

!!! example "例 54.2 (Jordan quiver 的表示)"
    Jordan quiver $L_1$ 的一个表示 = 一个有限维向量空间 $V$ 加上一个自同态 $f: V \to V$。即单个矩阵（在选定基后）。$L_1$ 的表示分类 = 矩阵的相似分类 = Jordan 标准形。

!!! example "例 54.3 (Kronecker quiver 的表示)"
    Kronecker quiver 的一个表示 = 两个向量空间 $V_1, V_2$ 加上两个线性映射 $f, g: V_1 \to V_2$。在选定基后，这等价于一对矩阵 $(A, B)$ 在左乘可逆矩阵（基变换 $V_2$）和右乘可逆矩阵（基变换 $V_1$）下的等价分类。这就是**矩阵束**（matrix pencil）$A + \lambda B$ 的 Kronecker 标准形问题。

!!! example "例 54.4 ($A_3$ 型 quiver 的表示)"
    考虑 $Q: 1 \to 2 \to 3$。一个表示是 $(V_1, V_2, V_3, f: V_1 \to V_2, g: V_2 \to V_3)$。维数向量 $\mathbf{d} = (d_1, d_2, d_3)$。选定各空间的基后，$f$ 和 $g$ 分别由 $d_2 \times d_1$ 和 $d_3 \times d_2$ 矩阵表示，等价关系为同时基变换。

!!! definition "定义 54.6 (零表示)"
    **零表示** $0$ 是每个顶点放零空间 $V_i = 0$ 的表示。

!!! definition "定义 54.7 (简单表示)"
    对每个顶点 $i \in Q_0$，**简单表示** $S_i$ 定义为：$(S_i)_j = \delta_{ij}\mathbb{F}$（即仅在顶点 $i$ 上放一维空间，其余顶点放零空间），所有箭头对应的线性映射为零映射。

---

## 54.3 态射与直和

<div class="context-flow" markdown>

**核心问题**：如何定义 quiver 表示之间的"结构保持映射"？如何分解表示？

</div>

!!! definition "定义 54.8 (表示的态射)"
    设 $V = (V_i, f_\alpha)$ 和 $W = (W_i, g_\alpha)$ 是 quiver $Q$ 的两个表示。从 $V$ 到 $W$ 的一个**态射**（morphism）$\varphi: V \to W$ 是一族线性映射 $\{\varphi_i: V_i \to W_i\}_{i \in Q_0}$，使得对每条箭头 $\alpha: i \to j$，下图交换：
    $$\begin{CD}
    V_i @>{f_\alpha}>> V_j \\
    @V{\varphi_i}VV @VV{\varphi_j}V \\
    W_i @>{g_\alpha}>> W_j
    \end{CD}$$
    即 $\varphi_j \circ f_\alpha = g_\alpha \circ \varphi_i$ 对所有 $\alpha: i \to j$ 成立。

!!! definition "定义 54.9 (同构)"
    态射 $\varphi: V \to W$ 是**同构**（isomorphism），若每个 $\varphi_i$ 都是线性同构。此时记 $V \cong W$。

!!! definition "定义 54.10 (表示的直和)"
    表示 $V = (V_i, f_\alpha)$ 和 $W = (W_i, g_\alpha)$ 的**直和** $V \oplus W$ 定义为：

    - 顶点 $i$ 上的空间为 $V_i \oplus W_i$。
    - 箭头 $\alpha: i \to j$ 上的映射为 $f_\alpha \oplus g_\alpha: V_i \oplus W_i \to V_j \oplus W_j$。

!!! definition "定义 54.11 (不可分解表示)"
    非零表示 $V$ 称为**不可分解**（indecomposable）的，若 $V \cong V' \oplus V''$ 蕴含 $V' = 0$ 或 $V'' = 0$。

!!! theorem "定理 54.1 (Krull-Schmidt 定理)"
    quiver $Q$ 的任何有限维表示 $V$ 可以分解为不可分解表示的直和：
    $$V \cong V_1^{m_1} \oplus V_2^{m_2} \oplus \cdots \oplus V_k^{m_k},$$
    其中 $V_1, \ldots, V_k$ 是两两不同构的不可分解表示。分解是唯一的（在同构和重排序意义下）。

??? proof "证明"
    这是更一般的 Krull-Schmidt 定理在有限长度模范畴中的特例。

    **存在性**：对表示的维数向量的总维数 $\sum_i \dim V_i$ 进行归纳。若 $V$ 不可分解，则已经完成。否则 $V \cong V' \oplus V''$，其中 $V', V''$ 非零，总维数均严格小于 $V$，由归纳假设可分解。

    **唯一性**：关键是 Fitting 引理：对不可分解表示 $V$，$\mathrm{End}(V)$ 是局部环（即每个非同构的自同态都是幂零的）。有了这一点，唯一性可由标准的 Krull-Schmidt 论证得出：设
    $$V_1 \oplus V_2 \oplus \cdots \oplus V_m \cong W_1 \oplus W_2 \oplus \cdots \oplus W_l$$
    都是不可分解分解。设 $\varphi$ 是此同构，考虑 $\pi_{W_1} \circ \varphi \circ \iota_{V_1}: V_1 \to W_1$（其中 $\pi, \iota$ 是投影和嵌入）。利用 $\mathrm{End}(V_1)$ 的局部性，可以证明存在某个 $W_j$ 使得 $V_1 \cong W_j$，然后消去并归纳。$\blacksquare$

!!! definition "定义 54.12 (子表示)"
    $V$ 的**子表示** $W$ 是一族子空间 $\{W_i \subseteq V_i\}_{i \in Q_0}$，使得 $f_\alpha(W_{s(\alpha)}) \subseteq W_{t(\alpha)}$ 对所有 $\alpha \in Q_1$。

!!! definition "定义 54.13 (商表示)"
    给定子表示 $W \subseteq V$，**商表示** $V/W$ 定义为：顶点 $i$ 上的空间为 $V_i / W_i$，箭头 $\alpha$ 上的映射为 $\bar{f}_\alpha: V_i/W_i \to V_j/W_j$（良定义性由子表示条件保证）。

!!! theorem "定理 54.2 (Schur 引理)"
    设 $V, W$ 是不可分解表示。则：

    (a) 若 $V \not\cong W$，则 $\mathrm{Hom}(V, W)$ 中的任何态射都不是同构（更强地，若基域代数闭，则 $\mathrm{Hom}(V, W) = 0$ 当 $V \not\cong W$）。

    (b) 若基域 $\mathbb{F}$ 代数闭，则 $\mathrm{End}(V) / \mathrm{rad}(\mathrm{End}(V)) \cong \mathbb{F}$。

??? proof "证明"
    (a) 设 $\varphi: V \to W$ 是态射，$V$ 不可分解。则 $\ker \varphi$ 是 $V$ 的子表示，$\mathrm{Im}\,\varphi$ 是 $W$ 的子表示。若 $\varphi$ 是同构，则 $V \cong W$，矛盾。

    (b) 由 Fitting 引理，$\mathrm{End}(V)$ 中的非同构元素构成唯一的极大理想（即根 $\mathrm{rad}$）。商是除法代数。若 $\mathbb{F}$ 代数闭，有限维除法代数只有 $\mathbb{F}$ 自身。$\blacksquare$

---

## 54.4 经典矩阵问题作为 Quiver 表示

<div class="context-flow" markdown>

**核心问题**：Quiver 表示如何统一经典的矩阵分类问题？

</div>

!!! example "例 54.5 (单矩阵相似 = Jordan quiver)"
    **问题**：给定 $n \times n$ 矩阵 $A$，在相似变换 $A \mapsto PAP^{-1}$ 下分类。

    **Quiver 翻译**：Jordan quiver $L_1$（一个顶点一条自环）的表示 = $(\mathbb{F}^n, A)$。态射 = 可逆矩阵 $P$ 使得 $PA = AP$... 不，态射 $\varphi: (\mathbb{F}^n, A) \to (\mathbb{F}^n, B)$ 要求 $\varphi A = B \varphi$，即 $B = \varphi A \varphi^{-1}$。同构类 = 相似类。

    **结论**：Jordan 标准形 = Jordan quiver 的不可分解表示分类。维数为 $n$ 的不可分解表示恰好是 $n \times n$ Jordan 块 $J_n(\lambda)$（$\lambda \in \mathbb{F}$）。

!!! example "例 54.6 (矩阵对的同时等价 = Kronecker quiver)"
    **问题**：给定矩阵对 $(A, B)$，$A, B \in M_{m \times n}(\mathbb{F})$，在 $(A, B) \mapsto (PAQ^{-1}, PBQ^{-1})$ 下分类。

    **Quiver 翻译**：Kronecker quiver $1 \rightrightarrows 2$ 的维数向量 $(n, m)$ 的表示 = $(\mathbb{F}^n, \mathbb{F}^m, A, B)$。

!!! example "例 54.7 (子空间问题 = $A_n$ 型 quiver)"
    **问题**：给定向量空间 $V$ 和子空间链 $V_1 \subseteq V_2 \subseteq \cdots \subseteq V_n = V$，在 $V$ 的基变换下分类。

    **Quiver 翻译**：$A_n$ 型 quiver $1 \to 2 \to \cdots \to n$ 的表示，其中所有箭头对应的线性映射为嵌入。

!!! example "例 54.8 (四子空间问题 = $D_4$ 型 quiver)"
    **问题**：给定向量空间 $V$ 和四个子空间 $U_1, U_2, U_3, U_4 \subseteq V$，在 $V$ 的基变换下分类。

    **Quiver 翻译**：$D_4$ 型 quiver（四个叶子顶点各有一条箭头指向中心顶点）的表示。这个问题有有限多个不可分解表示——恰好 12 个（正根的个数）。

!!! theorem "定理 54.3 (表示空间)"
    固定 quiver $Q$ 和维数向量 $\mathbf{d} = (d_i)_{i \in Q_0}$。$Q$ 的所有维数向量为 $\mathbf{d}$ 的表示构成仿射空间
    $$\mathrm{Rep}(Q, \mathbf{d}) = \prod_{\alpha: i \to j} \mathrm{Hom}(\mathbb{F}^{d_i}, \mathbb{F}^{d_j}) \cong \prod_{\alpha: i \to j} M_{d_j \times d_i}(\mathbb{F}),$$
    维数为 $\sum_{\alpha: i \to j} d_i d_j$。

    群 $\mathrm{GL}(\mathbf{d}) = \prod_{i \in Q_0} \mathrm{GL}(d_i, \mathbb{F})$ 通过同时基变换作用在 $\mathrm{Rep}(Q, \mathbf{d})$ 上：
    $$(g_i)_{i \in Q_0} \cdot (f_\alpha)_{\alpha \in Q_1} = (g_{t(\alpha)} f_\alpha g_{s(\alpha)}^{-1})_{\alpha \in Q_1}.$$
    **同构类 = $\mathrm{GL}(\mathbf{d})$-轨道**。

---

## 54.5 Gabriel 定理

<div class="context-flow" markdown>

**核心问题**：哪些 quiver 具有有限多个不可分解表示（有限表示型）？

</div>

!!! definition "定义 54.14 (有限表示型)"
    Quiver $Q$ 称为**有限表示型**（finite representation type），若 $Q$ 的不可分解表示（在同构意义下）只有有限多个。

!!! definition "定义 54.15 (Dynkin 图)"
    **Dynkin 图**是以下无向图：

    - $A_n$（$n \geq 1$）：$n$ 个顶点的链。$\bullet - \bullet - \cdots - \bullet$
    - $D_n$（$n \geq 4$）：$n-1$ 个顶点的链，加上从第 $n-2$ 个顶点分出的一个额外顶点。
    - $E_6$：6 个顶点，形如 $\bullet - \bullet - \overset{\bullet}{|} - \bullet - \bullet$（从中间分出一个支链）。
    - $E_7$：7 个顶点。
    - $E_8$：8 个顶点。

!!! definition "定义 54.16 (根系与正根)"
    对 Dynkin 图 $\Gamma$，令 $C = 2I - A$（其中 $A$ 是邻接矩阵）为 **Cartan 矩阵**。$\Gamma$ 的**根系** $\Phi$ 是 $\mathbb{Z}^n$ 中满足一定条件的向量集（由 Weyl 群的轨道生成）。**正根** $\Phi^+$ 是那些所有分量非负的根。

    $|\Phi^+|$ 的数量：$|A_n^+| = \binom{n+1}{2}$，$|D_n^+| = n(n-1)$，$|E_6^+| = 36$，$|E_7^+| = 63$，$|E_8^+| = 120$。

!!! theorem "定理 54.4 (Gabriel 定理, 1972)"
    设 $Q$ 是连通 quiver（无定向环），$\mathbb{F}$ 是代数闭域。则以下等价：

    **(a)** $Q$ 是有限表示型。

    **(b)** $Q$ 的底图 $\bar{Q}$ 是 Dynkin 图（$A_n, D_n, E_6, E_7, E_8$ 之一）。

    而且，当 $Q$ 是有限表示型时，不可分解表示的同构类与 $\bar{Q}$ 对应的根系的正根之间存在**一一对应**：不可分解表示 $V$ 的维数向量 $\mathbf{d}(V)$ 恰好遍历所有正根。

??? proof "证明"
    完整证明需要较多的代数背景，这里给出关键思路。

    **(b) $\Rightarrow$ (a)** 的证明思路：

    **第一步（Tits 二次型）**：定义 **Tits 二次型** $q: \mathbb{Z}^{Q_0} \to \mathbb{Z}$：
    $$q(\mathbf{d}) = \sum_{i \in Q_0} d_i^2 - \sum_{\alpha: i \to j} d_i d_j.$$
    这个二次型度量了"参数空间维数"和"对称群维数"的差：
    $$q(\mathbf{d}) = \dim \mathrm{GL}(\mathbf{d}) - \dim \mathrm{Rep}(Q, \mathbf{d}) + 1$$
    （减去了对角标量矩阵的一维）。

    **关键事实**：$\bar{Q}$ 是 Dynkin 图当且仅当 $q$ 正定。

    **第二步**：若 $q$ 正定，则 $q(\mathbf{d}) \geq 1$ 对所有非零 $\mathbf{d} \geq 0$。这意味着 $\dim \mathrm{GL}(\mathbf{d}) > \dim \mathrm{Rep}(Q, \mathbf{d})$（"对称性太多"），群轨道数有限。

    更精确地，可以证明：对每个正根 $\mathbf{d}$，$\mathrm{Rep}(Q, \mathbf{d})$ 中恰好有一个 $\mathrm{GL}(\mathbf{d})$-轨道是开稠密的（唯一不可分解表示）。对非正根的维数向量 $\mathbf{d}$，不存在不可分解表示。

    **第三步（反射函子）**：Bernstein-Gelfand-Ponomarev 引入的**反射函子**（reflection functors）提供了沿 Weyl 群作用在正根之间移动的工具。源点反射 $\sigma_i^+$ 和汇点反射 $\sigma_i^-$ 实现了不可分解表示之间的双射（除了简单表示 $S_i$）。

    **(a) $\Rightarrow$ (b)**：若 $\bar{Q}$ 不是 Dynkin 图，则要么包含扩展 Dynkin 图 $\tilde{A}_n, \tilde{D}_n, \tilde{E}_6, \tilde{E}_7, \tilde{E}_8$ 作为子图，要么包含更大的图。对扩展 Dynkin 图，$q$ 半正定且存在非零向量 $\mathbf{d}$ 使 $q(\mathbf{d}) = 0$，由此可以构造出参数族的不可分解表示（无穷多个）。$\blacksquare$

!!! example "例 54.9 ($A_3$ 型 quiver 的不可分解表示)"
    $Q: 1 \to 2 \to 3$。正根为：$(1,0,0), (0,1,0), (0,0,1), (1,1,0), (0,1,1), (1,1,1)$，共 $\binom{4}{2} = 6$ 个。

    对应的不可分解表示（在同构意义下唯一）：

    | 正根 | 表示 |
    |:-----|:-----|
    | $(1,0,0)$ | $\mathbb{F} \xrightarrow{0} 0 \xrightarrow{0} 0$ |
    | $(0,1,0)$ | $0 \xrightarrow{0} \mathbb{F} \xrightarrow{0} 0$ |
    | $(0,0,1)$ | $0 \xrightarrow{0} 0 \xrightarrow{0} \mathbb{F}$ |
    | $(1,1,0)$ | $\mathbb{F} \xrightarrow{1} \mathbb{F} \xrightarrow{0} 0$ |
    | $(0,1,1)$ | $0 \xrightarrow{0} \mathbb{F} \xrightarrow{1} \mathbb{F}$ |
    | $(1,1,1)$ | $\mathbb{F} \xrightarrow{1} \mathbb{F} \xrightarrow{1} \mathbb{F}$ |

!!! example "例 54.10 ($D_4$ 型 quiver 的不可分解表示)"
    $D_4$ 有 $|D_4^+| = 4 \times 3 = 12$ 个正根，因此有 12 个不可分解表示。维数向量包括：

    - 4 个简单表示：$(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)$。
    - 维数向量为 $(0,1,1,0), (0,1,0,1), (0,0,1,1), (1,1,0,0), (1,0,1,0), (1,0,0,1)$ 的 6 个表示。
    - 维数向量 $(1,1,1,1)$ 的 1 个表示。
    - 维数向量 $(1,2,1,1)$（或其排列）的 1 个表示。

    （具体取决于箭头方向。）

!!! theorem "定理 54.5 (扩展 Dynkin 图与驯型)"
    若 $\bar{Q}$ 是**扩展 Dynkin 图**（$\tilde{A}_n, \tilde{D}_n, \tilde{E}_6, \tilde{E}_7, \tilde{E}_8$），则 $Q$ 是**驯型**（tame representation type）：不可分解表示在每个维数向量下至多构成有限多个一参数族。

    若 $\bar{Q}$ 既非 Dynkin 图也非扩展 Dynkin 图，则 $Q$ 是**野型**（wild representation type）：其表示分类与任意有限维代数的表示分类一样困难。

---

## 54.6 Auslander-Reiten 箭图

<div class="context-flow" markdown>

**核心问题**：如何系统地描述不可分解表示之间的"不可约"态射？

</div>

!!! definition "定义 54.17 (不可约态射)"
    态射 $f: X \to Y$（$X, Y$ 不可分解）称为**不可约**（irreducible）的，若 $f$ 不是同构，且对任何分解 $f = gh$，要么 $g$ 是分裂满射，要么 $h$ 是分裂单射。

!!! definition "定义 54.18 (几乎分裂序列)"
    一个正合序列
    $$0 \to X \xrightarrow{f} Y \xrightarrow{g} Z \to 0$$
    称为**几乎分裂序列**（almost split sequence，或 Auslander-Reiten 序列），若：

    (i) $X, Z$ 不可分解。
    (ii) $f$ 不是分裂单射（即序列不分裂）。
    (iii) 对任何不是分裂满射的 $h: W \to Z$，存在 $h': W \to Y$ 使得 $gh' = h$。

!!! theorem "定理 54.6 (Auslander-Reiten 定理)"
    设 $Q$ 是没有定向环的有限 quiver，$Z$ 是不可分解表示且不是投射表示。则存在唯一（在同构意义下）的几乎分裂序列 $0 \to \tau Z \to Y \to Z \to 0$，其中 $\tau Z$ 称为 $Z$ 的 **Auslander-Reiten 平移**（AR translate）。

!!! definition "定义 54.19 (Auslander-Reiten 箭图)"
    Quiver $Q$ 的 **Auslander-Reiten 箭图**（AR quiver）$\Gamma_Q$ 定义为：

    - 顶点：$Q$ 的不可分解表示的同构类。
    - 箭头：$[X] \to [Y]$（当存在不可约态射 $X \to Y$ 时）。

!!! example "例 54.11 ($A_3$ 型 quiver $1 \to 2 \to 3$ 的 AR 箭图)"
    记 6 个不可分解表示为 $P_1 = (1,1,1)$，$P_2 = (0,1,1)$，$P_3 = (0,0,1)$（投射不可分解），$I_1 = (1,0,0)$，$I_2 = (1,1,0)$，$I_3 = (1,1,1)$。注意 $P_1 = I_3$。

    用维数向量标记，AR 箭图为：
    ```
    (1,0,0) → (1,1,0) → (1,1,1)
         ↗         ↗
    (0,1,0) → (0,1,1)
         ↗
    (0,0,1)
    ```

    AR 平移为：$\tau(1,1,1) = (0,0,1)$，$\tau(1,1,0) = (0,1,0) $，$\tau(0,1,1) = (1,0,0)$。投射表示 $(0,0,1), (0,1,0), (1,0,0)$... 让我们重新整理。

    对 quiver $1 \to 2 \to 3$，投射不可分解表示是 $P(1) = (1,1,1)$（维数向量），$P(2) = (0,1,1)$，$P(3) = (0,0,1)$。入射不可分解表示是 $I(1) = (1,0,0)$，$I(2) = (1,1,0)$，$I(3) = (1,1,1)$。

    几乎分裂序列：
    $$0 \to P(2) \to P(1) \oplus P(3) \to I(2) \to 0$$
    即 $0 \to (0,1,1) \to (1,1,1) \oplus (0,0,1) \to (1,1,0) \to 0$。

    这给出 $\tau(I(2)) = P(2)$，即 $\tau(1,1,0) = (0,1,1)$。

!!! remark "注记"
    AR 箭图在有限表示型情形完全决定了表示范畴的结构。对 $A_n$ 型 quiver，AR 箭图的形状与 $n$ 阶三角网格相关。

---

## 54.7 持久模与拓扑数据分析

<div class="context-flow" markdown>

**核心问题**：Quiver 表示论如何与拓扑数据分析中的持久同调联系？

</div>

!!! definition "定义 54.20 (持久模)"
    一个**持久模**（persistence module）是 $A_n$ 型 quiver $1 \to 2 \to \cdots \to n$ 的一个表示 $(V_1, V_2, \ldots, V_n; f_1, f_2, \ldots, f_{n-1})$，其中 $f_i: V_i \to V_{i+1}$ 是线性映射。

    在拓扑数据分析中，$V_i$ 通常是某个随参数 $\varepsilon_1 < \varepsilon_2 < \cdots < \varepsilon_n$ 变化的简单复形 $K_{\varepsilon_i}$ 的同调群 $H_k(K_{\varepsilon_i}; \mathbb{F})$，$f_i$ 是包含映射诱导的同调映射。

!!! theorem "定理 54.7 (区间分解定理)"
    设 $\mathbb{F}$ 是域。$A_n$ 型 quiver 的任何有限维表示可以唯一分解为**区间表示**的直和：
    $$V \cong \bigoplus_{[b, d]} \mathbb{I}_{[b,d]}^{m_{b,d}},$$
    其中 $\mathbb{I}_{[b,d]}$（$1 \leq b \leq d \leq n$）是区间 $[b, d]$ 上的不可分解表示：

    $$(\mathbb{I}_{[b,d]})_i = \begin{cases} \mathbb{F}, & b \leq i \leq d, \\ 0, & \text{otherwise}, \end{cases}$$

    箭头 $i \to i+1$ 的映射在 $b \leq i < d$ 时为恒等映射 $\mathrm{id}_\mathbb{F}$，否则为零映射。

??? proof "证明"
    这是 Gabriel 定理对 $A_n$ 型 quiver 的特例。$A_n$ 的正根恰好是 $\{e_b + e_{b+1} + \cdots + e_d : 1 \leq b \leq d \leq n\}$，共 $\binom{n+1}{2}$ 个，与区间 $[b, d]$ 一一对应。

    也可以直接证明：对 $A_n$ 的表示 $(V_i, f_i)$，通过对维数向量的归纳和选取适当的基来实现区间分解。关键观察是可以选取每个 $V_i$ 的一组基，使得每个基向量要么映射到下一空间的基向量（$f_i(v) = w$），要么映射到零。按基向量的"生存区间"分组，即得区间分解。$\blacksquare$

!!! definition "定义 54.21 (条形码和持久图)"
    持久模 $V$ 的区间分解
    $$V \cong \bigoplus_{k=1}^N \mathbb{I}_{[b_k, d_k]}$$
    可以用两种等价方式可视化：

    - **条形码**（barcode）：每个区间 $[b_k, d_k]$ 画一条水平线段。
    - **持久图**（persistence diagram）：每个区间 $[b_k, d_k]$ 对应平面上的点 $(b_k, d_k)$。

    长条（$d_k - b_k$ 大的区间）代表**持久**的拓扑特征（真实的信号），短条代表**短暂**的特征（噪声）。

!!! example "例 54.12 (点云的持久同调)"
    给定平面上的点云 $\{x_1, \ldots, x_N\} \subset \mathbb{R}^2$，取参数 $\varepsilon_1 < \varepsilon_2 < \cdots < \varepsilon_n$。

    对每个 $\varepsilon_i$，构造 **Vietoris-Rips 复形** $\mathrm{VR}_{\varepsilon_i}$：以点为顶点，若 $\|x_j - x_k\| \leq \varepsilon_i$ 则连边，高维单形类似。

    包含关系 $\mathrm{VR}_{\varepsilon_1} \subseteq \mathrm{VR}_{\varepsilon_2} \subseteq \cdots \subseteq \mathrm{VR}_{\varepsilon_n}$ 给出持久模
    $$H_k(\mathrm{VR}_{\varepsilon_1}) \to H_k(\mathrm{VR}_{\varepsilon_2}) \to \cdots \to H_k(\mathrm{VR}_{\varepsilon_n}).$$

    若点云采样自一个圆，则 $H_1$ 的条形码中会有一条明显的长条（对应圆的 1 维洞），其余为短条（噪声）。

!!! theorem "定理 54.8 (持久图的稳定性)"
    设 $f, g: X \to \mathbb{R}$ 是拓扑空间 $X$ 上的驯函数，$\mathrm{dgm}(f)$ 和 $\mathrm{dgm}(g)$ 分别是它们的持久图。则
    $$d_B(\mathrm{dgm}(f), \mathrm{dgm}(g)) \leq \|f - g\|_\infty,$$
    其中 $d_B$ 是**瓶颈距离**（bottleneck distance）。即：小的函数扰动只产生小的持久图变化。

??? proof "证明"
    这个定理的完整证明需要使用交错引理（interpolation lemma）和代数稳定性定理。核心思想是：

    函数值的 $\delta$-扰动最多使每个区间的端点移动 $\delta$，同时最多产生长度 $\leq 2\delta$ 的新区间或消灭长度 $\leq 2\delta$ 的旧区间。

    代数层面，这通过**交错模**（interleaving）来形式化：若 $\|f - g\|_\infty \leq \delta$，则 $f$ 和 $g$ 的持久模是 $\delta$-交错的，而 $\delta$-交错意味着瓶颈距离 $\leq \delta$。$\blacksquare$

!!! remark "注记"
    区间分解定理仅对 $A_n$ 型 quiver（即线性序的持久模）成立。对**多参数持久模**（multiparameter persistence，对应格形 quiver 的表示），一般没有这样简洁的分解，这是当前 TDA 研究的活跃方向。从 quiver 表示论的视角看，这是因为 $A_n$ 是 Dynkin 型（有限表示型），而二维网格 quiver 是野型。

---

## 本章小结

| 概念 | 定义/关键性质 |
|:---|:---|
| Quiver | 有向图 $Q = (Q_0, Q_1, s, t)$ |
| 表示 | 顶点 $\mapsto$ 向量空间，箭头 $\mapsto$ 线性映射 |
| 态射 | 交换图的线性映射族 |
| 不可分解表示 | 不能非平凡直和分解 |
| Gabriel 定理 | 有限型 $\Leftrightarrow$ Dynkin 图；不可分解 $\leftrightarrow$ 正根 |
| AR 箭图 | 不可分解表示间不可约态射的结构图 |
| 持久模 | $A_n$ 型表示；区间分解定理；条形码/持久图 |

---

## 习题

!!! exercise "习题 54.1"
    列出 $A_4$ 型 quiver $1 \to 2 \to 3 \to 4$ 的所有不可分解表示（给出维数向量和具体的线性映射），并验证共有 $\binom{5}{2} = 10$ 个。

!!! exercise "习题 54.2"
    对 Kronecker quiver $K_2$（两个顶点，两条同向箭头），描述维数向量 $(1,1)$ 的所有不可分解表示。说明为什么有无穷多个（参数化为 $\mathbb{P}^1$），从而 $K_2$ 不是有限表示型。验证 $K_2$ 的底图是扩展 Dynkin 图 $\tilde{A}_1$。

!!! exercise "习题 54.3"
    证明：$A_n$ 型 quiver 的 Tits 二次型 $q(d_1, \ldots, d_n) = \sum_{i=1}^n d_i^2 - \sum_{i=1}^{n-1} d_i d_{i+1}$ 是正定的。（提示：写成 $q = \frac{1}{2}[(d_1-d_2)^2 + (d_2-d_3)^2 + \cdots + d_1^2 + d_n^2]$ 或直接计算对应的矩阵。）

!!! exercise "习题 54.4"
    设 $Q$ 是 $D_4$ 型 quiver（三个叶子顶点各有一条箭头指向中心），中心为顶点 $0$，叶子为 $1, 2, 3$。

    (a) 计算 Tits 二次型并验证它正定。

    (b) 列出 $D_4$ 的所有 12 个正根。

    (c) 对维数向量 $(1,1,1,2)$（$d_0 = 2$，$d_1 = d_2 = d_3 = 1$），写出不可分解表示的具体形式。

!!! exercise "习题 54.5"
    构造一个 $A_5$ 型持久模的具体例子，使其条形码恰好由区间 $[1,3], [2,5], [4,4]$ 组成。写出每个顶点上的向量空间维数和每条箭头的线性映射矩阵。

!!! exercise "习题 54.6"
    证明 Schur 引理（定理 54.2(a)）：若 $V, W$ 是不可分解表示且 $V \not\cong W$，则从 $V$ 到 $W$ 不存在同构的态射。

!!! exercise "习题 54.7"
    验证 $\tilde{A}_2$（循环 quiver $1 \to 2 \to 3 \to 1$）不是有限表示型：对维数向量 $(1,1,1)$，描述所有不可分解表示，说明它们构成一个一参数族（参数化为 $\mathbb{F} \setminus \{0\}$ 或 $\mathbb{P}^1$ 减去三个点）。

!!! exercise "习题 54.8"
    （开放性）多参数持久模——即 $\mathbb{Z}^2$-指标的持久模——为什么没有类似区间分解定理的结果？从 quiver 表示论的角度解释原因。
