# 第 49 章 外代数与 Grassmannian

<div class="context-flow" markdown>

**前置**：多线性代数 (Ch21) · 行列式 (Ch03) · 子空间与基 (Ch04)

**本章脉络**：从多重线性到反对称性 $\to$ 外积 (Exterior Product) $\wedge$ 的定义与公理 $\to$ 外代数 (Exterior Algebra) / Grassmann 代数 $\Lambda(V)$ $\to$ 基的构造与维数 $\binom{n}{k}$ $\to$ 行列式的外代数起源 $\to$ 可解张量 (Decomposable Tensors) 与子空间的对应 $\to$ Grassmannian 流形 $Gr(k, V)$ 的代数描述 $\to$ Plücker 坐标与 Plücker 关系式 $\to$ 应用：微分形式、电磁场张量、子空间聚类、射影几何

**延伸**：外代数是描述“有向体积”的自然语言；它通过引入反对称性，将子空间本身看作代数对象，是现代微分几何、广义相对论以及理论计算机科学中网络流分析的核心

</div>

在经典代数中，我们研究向量的组合。在**外代数**（Exterior Algebra）中，我们将视角提升到“向量张成的平行多面体”。通过引入**外积**（楔积） $\wedge$，我们能够以纯代数的方式处理面积、体积及其推广。这种理论最终引出了 **Grassmannian**——即所有 $k$ 维子空间构成的流形，为我们在算子空间中进行“子空间运算”提供了终极框架。

---

## 49.1 外积与外代数

!!! definition "定义 49.1 (外积)"
    向量 $v, w \in V$ 的**外积** $v \wedge w$ 满足：
    1.  **双线性**：$(av_1+bv_2) \wedge w = a(v_1 \wedge w) + b(v_2 \wedge w)$。
    2.  **反对称性**：$v \wedge w = -w \wedge v$。
    3.  **幂零性**：$v \wedge v = 0$。

!!! note "外代数 $\Lambda(V)$"
    所有 $k$ 阶外幂 $\Lambda^k(V)$ 的直和构成了外代数。其总维数为 $2^n$。

---

## 49.2 几何意义与子空间

!!! theorem "定理 49.1 (线性相关判定)"
    向量组 $\{v_1, \ldots, v_k\}$ 线性无关，当且仅当 $v_1 \wedge v_2 \wedge \cdots \wedge v_k \neq 0$。
    **几何直观**：非零的 $k$ 阶外积代表了一个具有确定朝向和 $k$ 维体积的子空间片段。

---

## 49.3 Grassmannian 与 Plücker 坐标

!!! definition "定义 49.2 (Grassmannian)"
    $Gr(k, V)$ 是 $V$ 的所有 $k$ 维子空间构成的集合。
    每一个 $k$ 维子空间都可以由其基向量的外积 $v_1 \wedge \cdots \wedge v_k$（一个纯张量）唯一表示（至多差一个常数倍）。

!!! technique "Plücker 坐标"
    将 $v_1 \wedge \cdots \wedge v_k$ 在 $\Lambda^k(V)$ 的基下展开，所得的系数称为 **Plücker 坐标**。它们满足一组复杂的二次等式，称为 **Plücker 关系式**。

---

## 练习题

**1. [基础] 计算 $(e_1 + e_2) \wedge (e_1 - e_2)$。**

??? success "参考答案"
    **计算步骤：**
    1. 利用分配律：$e_1 \wedge e_1 - e_1 \wedge e_2 + e_2 \wedge e_1 - e_2 \wedge e_2$。
    2. 利用幂零性：$e_1 \wedge e_1 = 0, e_2 \wedge e_2 = 0$。
    3. 利用反对称性：$e_2 \wedge e_1 = -e_1 \wedge e_2$。
    4. 代入：$0 - e_1 \wedge e_2 - e_1 \wedge e_2 - 0 = -2 e_1 \wedge e_2$。
    **结论**：结果为 $-2 e_1 \wedge e_2$。

**2. [维数] 若 $\dim V = 4$，计算 $\Lambda^2(V)$ 的维数并列出基。**

??? success "参考答案"
    **计算：**
    1. 维数为 $\binom{4}{2} = \frac{4 \times 3}{2 \times 1} = 6$。
    2. 设基为 $\{e_1, e_2, e_3, e_4\}$。
    **基向量**：$\{e_1 \wedge e_2, e_1 \wedge e_3, e_1 \wedge e_4, e_2 \wedge e_3, e_2 \wedge e_4, e_3 \wedge e_4\}$。

**3. [行列式] 证明 $n$ 阶方阵的行列式是 $n$ 个列向量的外积。**

??? success "参考答案"
    **代数映射：**
    1. 考虑 $v_1 \wedge v_2 \wedge \cdots \wedge v_n$。
    2. 将 $v_j = \sum a_{ij} e_i$ 代入并利用外积的反对称性展开。
    3. 只有当指标 $(i_1, \ldots, i_n)$ 是 $\{1, \ldots, n\}$ 的排列时，项才非零。
    4. 每一项的系数正好带有一个排列的符号 $\operatorname{sgn}(\sigma)$。
    **结论**：$v_1 \wedge \cdots \wedge v_n = \det(A) (e_1 \wedge \cdots \wedge e_n)$。这揭示了行列式的最本质定义：它是顶级外幂的缩放因子。

**4. [判定] 判定 $v = e_1 \wedge e_2 + e_3 \wedge e_4$ 是否是可分解的（即是否能写成 $u \wedge w$）。**

??? success "参考答案"
    **Plücker 判定：**
    1. 在 4 维空间中，2 阶外幂是可分解的当且仅当 $v \wedge v = 0$。
    2. 计算 $v \wedge v = (e_1 \wedge e_2 + e_3 \wedge e_4) \wedge (e_1 \wedge e_2 + e_3 \wedge e_4)$。
    3. 展开得：$e_1 \wedge e_2 \wedge e_1 \wedge e_2 + e_1 \wedge e_2 \wedge e_3 \wedge e_4 + e_3 \wedge e_4 \wedge e_1 \wedge e_2 + e_3 \wedge e_4 \wedge e_3 \wedge e_4$。
    4. 重复项归零，余下 $2 e_1 \wedge e_2 \wedge e_3 \wedge e_4$。
    **结论**：由于 $v \wedge v \neq 0$，该元素不可分解。这意味着它不对应于单一的 2 维子空间，而是子空间的叠加。

**5. [性质] 证明：若 $v_1, \ldots, v_k$ 线性相关，则 $v_1 \wedge \cdots \wedge v_k = 0$。**

??? success "参考答案"
    **证明：**
    1. 若线性相关，则存在某个向量可由其余向量线性表示，设 $v_1 = \sum_{i=2}^k c_i v_i$。
    2. 代入外积：$(\sum c_i v_i) \wedge v_2 \wedge \cdots \wedge v_k$。
    3. 利用分配律，每一项都包含重复的向量（如 $c_2 v_2 \wedge v_2 \wedge \cdots$）。
    4. 由幂零性，所有项均为 0。

**6. [Grassmannian] $Gr(1, V)$ 对应于什么几何对象？**

??? success "参考答案"
    **结论：射影空间 $P(V)$。**
    $Gr(1, V)$ 代表所有过原点的直线。在几何上，这正是射影空间（Projective Space）的定义。

**7. [Plücker关系] 写出 $Gr(2, 4)$ 的唯一一个 Plücker 关系式。**

??? success "参考答案"
    **公式：**
    设坐标为 $p_{ij}$。
    $p_{12}p_{34} - p_{13}p_{24} + p_{14}p_{23} = 0$。
    这就是判定一个 6 维向量是否代表一个 2 维子空间的代数方程。

**8. [对偶] $\Lambda^k(V)$ 与 $\Lambda^{n-k}(V)$ 之间有什么关系？**

??? success "参考答案"
    **结论：它们是同构的。**
    **理由**：维数相等 $\binom{n}{k} = \binom{n}{n-k}$。
    在配备了内积的空间中，这种对应关系由 **Hodge 对偶**（Hodge Star Operator）算子 $\star$ 给出。

**9. [计算] 在 $\mathbb{R}^3$ 中， $v \wedge w$ 与叉积 $v \times w$ 有何联系？**

??? success "参考答案"
    **联系：**
    外积 $v \wedge w$ 是一个 2 阶张量（代表有向面积元），属于 $\Lambda^2(\mathbb{R}^3)$。
    叉积 $v \times w$ 是一个向量，属于 $\mathbb{R}^3$。
    在 3 维空间中，通过 Hodge 对偶，可以将面积元唯一映射到其法向量上。因此，叉积本质上是外积在 3 维空间的对偶表现。

**10. [应用] 简述外代数在电磁学中的应用。**

??? success "参考答案"
    在相对论电动力学中，电场和磁场被统一为一个 2 阶反对称张量（法拉第张量 $F$），它其实是 4 维时空外代数中的一个元素。麦克斯韦方程组可以极其简洁地写为 $dF = 0$ 和 $d{\star F} = J$，这完美体现了外代数在描述通量与环流方面的天然优势。

## 本章小结

外代数是几何直观的代数化顶峰：

1.  **子空间的代数化**：通过外积，我们将抽象的“子空间”转化为具体的“代数元素”，实现了从研究点到研究空间碎片的跨越。
2.  **体积的算子化**：行列式的起源和性质在外代数框架下得到了最彻底的解释，证明了反对称性是多维测度论的代数本质。
3.  **流形的基石**：Grassmannian 流形及其 Plücker 坐标为现代几何与拓扑提供了精细的局部刻画工具，是连接纯代数与现代物理（如弦论）的重要纽带。
