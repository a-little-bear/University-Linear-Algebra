# 第 45A 章 完全正矩阵 (Completely Positive)

<div class="context-flow" markdown>

**前置**：正定矩阵 (Ch16) · 非负矩阵 (Ch17) · 矩阵分解 (Ch10)

**本章脉络**：从正定性到非负性分解 $\to$ 完全正矩阵 (CP Matrices) 的定义 $\to$ 与一般正定非负矩阵 (Copositive) 的区别 $\to$ 核心判据：存在非负分解 $A = BB^T$ $\to$ CP-秩 (CP-rank) 的定义与界限 $\to$ 与非负矩阵分解 (NMF) 的关系 $\to$ 应用：图论中的团数估计、组合优化中的半正定规划 (SDP) 松弛、概率分布的混合模型

**延伸**：完全正矩阵是正定矩阵空间中具有“非负灵魂”的子集；它要求算子不仅在能量上为正，且必须能由纯正的非负基向量构造而成，是连接连续凸优化与离散图论的绝佳桥梁

</div>

在 Ch16 中我们学到，如果 $A = BB^T$，则 $A$ 是正定的。如果在此基础上，我们进一步要求 $B$ 的每一个元素都必须是非负的（即 $B \ge 0$），我们便得到了一类特殊的矩阵——**完全正矩阵**（Completely Positive Matrices, CP）。尽管这类矩阵看起来只是在分解上加了一个小小的限制，但其背后的复杂性极高，且与图论中的许多难题（如寻找最大团）有着深刻的代数联系。

---

## 45A.1 定义与核心判据

!!! definition "定义 45A.1 (完全正矩阵)"
    方阵 $A$ 称为 **完全正矩阵**（Completely Positive），如果存在一个 $n \times k$ 的**非负矩阵** $B$ 使得：
    $$A = B B^T, \quad B_{ij} \ge 0$$

!!! note "辨析：CP vs. PSD+Nonnegative"
    - **正定非负阵**：既正定又每一项都非负（即 $A \succeq 0$ 且 $A \ge 0$）。
    - **完全正矩阵**：不仅 $A \succeq 0$ 且 $A \ge 0$，还要求其分解因子也非负。
    **结论**：所有的 CP 矩阵都是 PSD 且非负的，但反之不一定成立（对于 $n \ge 5$）。

---

## 45A.2 CP-秩 (CP-rank)

!!! definition "定义 45A.2 (CP-秩)"
    使得 $A = BB^T$ 成立的非负矩阵 $B$ 的最小列数 $k$ 称为 $A$ 的 **CP-秩**。
    **界限**：一般地， $r \le \text{cp-rank}(A) \le n(n+1)/2$。计算精确的 CP-秩是一个极具挑战性的 NP-难问题。

---

## 45A.3 图论应用

!!! technique "应用：图的完全正性"
    每一个图 $G$ 都有一个关联的矩阵类。若一个图对应的所有 PSD 非负矩阵都是完全正的，则该图称为 **CP-图**。
    **重要结论**：一个图是 CP-图当且仅当它不包含长度大于 4 的奇环（即它是二部图或某些特定的弦图）。

---

## 练习题

**1. [基础] 判定 $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$ 是否为完全正矩阵。**

??? success "参考答案"
    **构造分解：**
    1. 观察到 $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \end{pmatrix}$。
    2. 分解因子 $B = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。
    3. 由于 $B$ 的所有元素均为正数。
    **结论**：是的，它是完全正矩阵。其 CP-秩为 1。

**2. [对比] 举出一个 PSD 且非负，但不是 CP 的矩阵（针对 $n=5$）。**

??? success "参考答案"
    **说明：**
    这是全正矩阵理论中最著名的发现之一。对于 $n \le 4$， $A \succeq 0$ 且 $A \ge 0$ 蕴含 $A$ 是 CP 的。
    但在 $n=5$ 时，存在如 **Horn 矩阵** 的变体，虽然矩阵每一项都大于 0 且特征值非负，但由于其关联图包含一个 5 阶奇环（五边形），无法找到非负的 $BB^T$ 分解。

**3. [性质] 证明：完全正矩阵的对角线元素必须非负。**

??? success "参考答案"
    **证明：**
    1. 由定义 $a_{ii} = \sum_{j=1}^k B_{ij}^2$。
    2. 由于实数的平方非负，故 $a_{ii} \ge 0$。
    事实上，由于 CP 阵是半正定的，这一性质天然满足。

**4. [CP-秩] 计算 $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ 的 CP-秩。**

??? success "参考答案"
    **分析：**
    1. 尝试 $k=2$ 的分解：$A = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}^T$。这个分解含有负数，不满足。
    2. 重新构造：$A = \begin{pmatrix} 1 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 0 & 1 \end{pmatrix}$ 也不对。
    3. 注意：对于 $2 \times 2$ 阵，$A \ge 0$ 且 $A \succeq 0$ 即为 CP。
    4. 事实上，该阵可写为 $\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix} + \operatorname{diag}(1, 1) \ldots$
    5. 一个简单的 CP 分解是 $B = \begin{pmatrix} \sqrt{1.5} & \sqrt{0.5} \\ \sqrt{0.5} & \sqrt{1.5} \end{pmatrix}$。
    **结论**：其 CP-秩等于普通秩 2。

**5. [判定] 证明：若 $A$ 是 CP 矩阵，则 $A$ 必须是半正定的。**

??? success "参考答案"
    **证明：**
    1. 对于任何向量 $x$， $x^T A x = x^T (BB^T) x = (x^T B)(B^T x) = \|B^T x\|^2$。
    2. 范数的平方总是 $\ge 0$。
    **结论**：CP 矩阵是 PSD 矩阵的一个真子集。

**6. [阿达马积] 两个 CP 矩阵的 Hadamard 积还是 CP 矩阵吗？**

??? success "参考答案"
    **是的。**
    **理由**：若 $A = BB^T$ 且 $C = DD^T$，则 $A \circ C$ 可以通过将 $B$ 的行与 $D$ 的行进行 Kronecker 积构造出新的非负分解因子。这一性质是 CP 锥在 Hadamard 运算下封闭的代数保证。

**7. [应用] 为什么 CP 矩阵与“最大团问题”有关？**

??? success "参考答案"
    **解释：**
    图 $G$ 的团数（Clique Number） $\omega(G)$ 可以通过求解一个关于 CP 矩阵的优化问题精确得到：
    $\omega(G) = \max \{ \mathbf{1}^T X \mathbf{1} : \operatorname{tr}(X)=1, X \in CP, X_{ij}=0 \text{ if } (i,j) \notin E \}$。
    这一等价性将组合爆炸难题转化为了连续空间的矩阵判定问题。

**8. [NMF] CP 分解与非负矩阵分解（NMF）有什么联系？**

??? success "参考答案"
    **联系：**
    CP 分解是 NMF 的一种特殊对称形式。在 NMF 中我们分解 $V \approx WH$，而在 CP 中我们要求分解因子的对称性 $V = BB^T$。CP 理论为 NMF 的存在性、唯一性和秩的估计提供了深层的代数支撑。

**9. [凸性] 证明所有 $n \times n$ CP 矩阵构成的集合是一个凸锥。**

??? success "参考答案"
    **证明思路：**
    1. 若 $A, C$ 为 CP，则 $A = BB^T, C = DD^T$。
    2. 对于 $\alpha, \beta > 0$， $\alpha A + \beta C = \begin{pmatrix} \sqrt{\alpha}B & \sqrt{\beta}D \end{pmatrix} \begin{pmatrix} \sqrt{\alpha}B & \sqrt{\beta}D \end{pmatrix}^T$。
    3. 新的分解因子显然是非负的。
    **结论**：CP 矩阵集对加法和正数乘法封闭，故为凸锥。

**10. [极限] 随着 $n$ 增加，判定一个矩阵是否为 CP 的计算复杂度如何变化？**

??? success "参考答案"
    **结论**：变得**极其困难 (NP-hard)**。
    虽然对于小维度有简单的解析准则，但对于一般维度，目前尚无已知的多项式时间算法能准确判定一个 PSD 非负阵是否具备非负分解。

## 本章小结

完全正矩阵是线性代数中最具“结构约束”的算子类别之一：

1.  **分解的纯净性**：它要求算子的正性必须能溯源到非负的原子成分，这种“源头的清白”导致了其在描述概率混合和组合计数中的特殊地位。
2.  **图论的代数化**：CP 理论证明了图的拓扑结构（如奇环的存在）能够跨越领域，直接决定其关联矩阵的代数分解性质。
3.  **优化的终极挑战**：作为半正定规划中最具吸引力但也最难处理的锥之一，CP 矩阵的研究标志着我们从“计算矩阵”走向“理解矩阵构成规律”的高级阶段。
