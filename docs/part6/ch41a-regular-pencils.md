# 第 41A 章 正则矩阵束

<div class="context-flow" markdown>

**前置**：特征值 (Ch06) · λ-矩阵 (Ch13B) · 广义特征值问题

**本章脉络**：从单一矩阵到矩阵对 $\to$ 矩阵束 (Matrix Pencil) $A - \lambda B$ 的定义 $\to$ 正则矩阵束 (Regular Pencils) 与奇异矩阵束 $\to$ 广义特征值与广义特征向量 $\to$ 广义特征方程 $\det(A - \lambda B) = 0$ $\to$  Weierstrass 标准形（正则情形） $\to$ 无穷大特征值的处理 $\to$ 应用：结构动力学、广义线性系统状态空间分析

**延伸**：矩阵束理论是研究广义线性动力系统（如包含代数约束的微分方程组 DAEs）的数学基础；它将算子的谱理论从 $A$ 推广到了两个算子 $A, B$ 的相对演化规律，是复杂系统稳定性分析的关键

</div>

在经典特征值问题中，我们研究 $A\mathbf{x} = \lambda \mathbf{x}$。但在许多工程问题（如有限元分析）中，方程具有形式 $A\mathbf{x} = \lambda B\mathbf{x}$。**矩阵束**（Matrix Pencil）$A - \lambda B$ 正是描述这类相对特征关系的工具。当特征方程不恒等于零时，我们称其为**正则矩阵束**。本章将介绍正则矩阵束的分解理论及其在动力学系统中的支配作用。

---

## 41A.1 矩阵束的基本概念

!!! definition "定义 41A.1 (矩阵束)"
    两个 $m \times n$ 矩阵 $A, B$ 定义的集合 $\{ A - \lambda B : \lambda \in \mathbb{C} \}$ 称为 **矩阵束**。

!!! definition "定义 41A.2 (正则矩阵束)"
    若 $A, B$ 是同阶方阵，且特征多项式 $p(\lambda) = \det(A - \lambda B)$ 不恒为零，则称该矩阵束是**正则的**（Regular）。否则称为**奇异的**（Singular）。

---

## 41A.2 广义特征值与标准形

!!! definition "定义 41A.3 (广义特征值)"
    满足 $\det(A - \lambda B) = 0$ 的标量 $\lambda$ 称为**有限广义特征值**。
    若 $\det(B) = 0$，则该矩阵束还可能具有**无穷大特征值** $\lambda = \infty$。

!!! theorem "定理 41A.1 (Weierstrass 标准形)"
    对于正则矩阵束 $A - \lambda B$，存在非奇异阵 $P, Q$ 使得：
    $$P(A - \lambda B)Q = \operatorname{diag}(J - \lambda I, I - \lambda N)$$
    - $J$ 是 Jordan 标准形，对应有限特征值。
    - $N$ 是幂零 Jordan 块，对应无穷大特征值。

---

## 41A.3 动力学应用

!!! technique "应用：振动分析"
    在机械振动中，方程为 $M \ddot{x} + K x = 0$。假设解为 $x = e^{i\omega t} v$，则得到广义特征值问题 $Kv = \omega^2 Mv$。这里 $K$ 是刚度阵，$M$ 是质量阵。

---

## 练习题

**1. [计算] 计算矩阵束 $A - \lambda B$ 的特征多项式，其中 $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}, B = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。**

??? success "参考答案"
    **计算步骤：**
    1. 写出差分阵：$A - \lambda B = \begin{pmatrix} 1-\lambda & 0 \\ 0 & 2 \end{pmatrix}$。
    2. 计算行列式：$\det(A - \lambda B) = (1-\lambda) \cdot 2 = 2 - 2\lambda$。
    **结论**：特征多项式为 $2 - 2\lambda$。

**2. [特征值] 求上题中的广义特征值（包含无穷大特征值）。**

??? success "参考答案"
    **分析：**
    1. 有限特征值：令 $2 - 2\lambda = 0 \implies \lambda = 1$。
    2. 无穷大特征值检查：由于 $\det(B) = \det \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = 0$，存在无穷大特征值。
    3. 次数判定：由于 $n=2$ 但特征多项式次数为 1，缺失的 1 个根即为 $\infty$。
    **结论**：特征值为 $\{1, \infty\}$。

**3. [正则判定] 判定 $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, B = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ 是否为正则矩阵束。**

??? success "参考答案"
    **计算：**
    $A - \lambda B = \begin{pmatrix} 1 & -\lambda \\ 0 & 0 \end{pmatrix}$。
    $\det(A - \lambda B) = 1 \cdot 0 - (-\lambda) \cdot 0 = 0$。
    **结论**：由于行列式恒为 0，该矩阵束是**奇异的**（Singular）。

**4. [性质] 若 $B$ 可逆，广义特征值问题与普通特征值问题有何联系？**

??? success "参考答案"
    **结论：**
    若 $B$ 可逆， $A\mathbf{x} = \lambda B\mathbf{x}$ 等价于 **$(B^{-1}A)\mathbf{x} = \lambda \mathbf{x}$**。
    此时所有广义特征值都是有限的，且正是矩阵 $B^{-1}A$ 的普通特征值。

**5. [无穷大] 矩阵束具有无穷大特征值的物理意义是什么？**

??? success "参考答案"
    **解释：**
    在微分代数系统 (DAEs) 中，无穷大特征值通常对应于**代数约束**。
    它们代表了响应速度无限快的变量（即瞬时响应），或者是系统阶数发生退化的表现。在 Weierstrass 标准形中，它们由幂零部分 $N$ 描述。

**6. [计算] 求 $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ 相对于 $B = I$ 的广义特征值。**

??? success "参考答案"
    **计算：**
    由于 $B=I$，这退化为 $A$ 的普通特征值问题。
    $\det(A - \lambda I) = \lambda^2 - 1 = 0 \implies \lambda = \pm 1$。

**7. [Weierstrass] 正则矩阵束的标准形中，$I - \lambda N$ 块的 $N$ 满足什么性质？**

??? success "参考答案"
    **性质：**
    $N$ 是一个**幂零矩阵**（即存在 $k$ 使得 $N^k = O$）。
    其非零元素仅分布在对角线的上方（次对角线），对应于 $\infty$ 特征值的 Jordan 链。

**8. [稳定性] 若有限广义特征值的实部均小于 0，且无 $\infty$ 特征值，该系统稳定吗？**

??? success "参考答案"
    **是的。**
    这保证了微分部分的演化是衰减的，且没有代数约束导致的脉冲项（impulse）或阶数不匹配。

**9. [对角化] 什么条件下两个矩阵 $A, B$ 可以同时对角化？**

??? success "参考答案"
    **结论：**
    若 $A, B$ 均为厄米阵且其中一个（通常是 $B$）是正定的，则它们可以同时对角化。
    **意义**：这正是机械振动中“振型分解”的代数基础。

**10. [应用] 在控制理论中，广义特征值如何用于确定系统的传输零点？**

??? success "参考答案"
    **联系：**
    系统的零点定义为使系统矩阵束（Rosenbrock 矩阵）亏秩的复数值 $s$。
    这本质上是在特定的广义矩阵束中寻找广义特征值的过程。

## 本章小结

正则矩阵束理论是广义线性系统的终极图谱：

1.  **谱的相对性**：它将特征值从算子的固有属性提升为算子间相互干扰的度量，建立了描述物理系统相对演化的通用框架。
2.  **无穷大的解析**：通过引入 $\infty$ 特征值和幂零结构 $N$，矩阵束理论完美刻画了连续系统中的突变、约束与奇异行为。
3.  **结构的标准化**：Weierstrass 标准形为处理复杂的微分代数方程组提供了坐标变换的终极参考，实现了动力学系统各模态的彻底解耦。
