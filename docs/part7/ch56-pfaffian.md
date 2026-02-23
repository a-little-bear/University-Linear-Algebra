# 第 56 章 Pfaffian

<div class="context-flow" markdown>

**前置**：行列式 (Ch03) · 积和式 (Ch40A) · 反对称矩阵 (Ch02)

**本章脉络**：反对称矩阵的特殊行列式 $\to$ Pfaffian 的定义与多项式构造 $\to$ 基本性质：$\operatorname{Pf}(A)^2 = \det(A)$ $\to$ Pfaffian 的递归展开与指标匹配 $\to$ 代数恒等式：$\operatorname{Pf}(M^T A M) = \det(M)\operatorname{Pf}(A)$ $\to$ Pfaffian 与图中完美匹配的关系 $\to$ 应用：统计力学中的 Ising 模型、平面图的完美匹配计数（FKT 算法）、量子力学中的超对称

**延伸**：Pfaffian 是反对称矩阵行列式的“代数平方根”；它不仅填补了反对称系统在行列式开方后的符号模糊性，还通过其独特的指标匹配结构，成为了连接矩阵代数与统计物理中二聚体问题的核心纽带

</div>

在讨论对称矩阵时，我们关注特征值。但在讨论**反对称矩阵**（Skew-symmetric Matrices，$A^T = -A$）时，一个令人惊讶的现象出现了：其行列式总是一个完全平方式。**Pfaffian** 正是这个平方根的代数表达。它在形式上比行列式更紧凑，且在处理图论中的匹配问题以及描述物理系统中的费米子对时，具有无可替代的作用。

---

## 56.1 Pfaffian 的定义

!!! definition "定义 56.1 (Pfaffian)"
    对于 $2n$ 阶反对称矩阵 $A$，其 **Pfaffian** 定义为：
    $$\operatorname{Pf}(A) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{\sigma(2i-1), \sigma(2i)}$$
    它是一个关于矩阵元素的 $n$ 次齐次多项式。

!!! theorem "定理 56.1 (核心恒等式)"
    对于任何反对称矩阵 $A$：
    $$\operatorname{det}(A) = [\operatorname{Pf}(A)]^2$$
    这意味着 Pfaffian 确定了反对称矩阵行列式的符号根。

---

## 56.2 基本性质

!!! note "代数算律"
    1.  **缩放**：$\operatorname{Pf}(\lambda A) = \lambda^n \operatorname{Pf}(A)$。
    2.  **基变换**：$\operatorname{Pf}(M A M^T) = \det(M) \operatorname{Pf}(A)$。
    3.  **分块对角阵**：$\operatorname{Pf}(\operatorname{diag}(A, B)) = \operatorname{Pf}(A)\operatorname{Pf}(B)$。

---

## 56.3 图论与匹配

!!! technique "应用：平面图匹配计数"
    FKT 算法证明了，对于平面图，可以通过给边分配特定的方向（Pfaffian 取向），使得其关联反对称矩阵的 Pfaffian 恰好等于图中**完美匹配的数量**。
    这实现了从指数级计数难题到多项式级行列式运算的惊人转化。

---

## 练习题

**1. [基础] 计算 $2 \times 2$ 反对称阵 $A = \begin{pmatrix} 0 & a \\ -a & 0 \end{pmatrix}$ 的 Pfaffian。**

??? success "参考答案"
    **计算步骤：**
    1. 计算行列式：$\det(A) = 0 \cdot 0 - a(-a) = a^2$。
    2. 根据定义 $\operatorname{Pf}(A)^2 = \det(A)$，故 $\operatorname{Pf}(A) = \pm a$。
    3. 根据 2 阶显式定义：$\operatorname{Pf}(A) = a_{12}$。
    **结论**：$\operatorname{Pf}(A) = a$。

**2. [维度] 证明：奇数阶反对称矩阵的 Pfaffian 没有定义（或者说等于 0）。**

??? success "参考答案"
    **证明：**
    1. 若 $A$ 是 $n$ 阶反对称阵且 $n$ 为奇数。
    2. $\det(A) = \det(A^T) = \det(-A) = (-1)^n \det(A) = -\det(A)$。
    3. 故 $\det(A) = 0$。
    4. 由于 $\operatorname{Pf}(A)^2 = \det(A)$，若强行定义，其值必为 0。
    **结论**：Pfaffian 理论主要关注偶数维空间。

**3. [计算] 计算标准辛阵 $J = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} \oplus \cdots \oplus \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ 的 Pfaffian。**

??? success "参考答案"
    **利用性质：**
    1. 分块对角阵的 Pfaffian 等于各块之积。
    2. 每一块 $J_2 = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ 的 Pfaffian 为 1。
    **结论**：$\operatorname{Pf}(J) = 1 \cdot 1 \cdots 1 = 1$。这也是辛几何中体积形式的正向定义。

**4. [性质] 若将反对称矩阵的第 1 行与第 2 行交换（同时第 1 列与第 2 列也对应交换），Pfaffian 如何变化？**

??? success "参考答案"
    **结论：变号（乘以 -1）。**
    **理由**：这种操作相当于 $M A M^T$ 变换，其中 $M$ 是交换矩阵且 $\det(M) = -1$。
    根据公式 $\operatorname{Pf}(M A M^T) = \det(M) \operatorname{Pf}(A)$，结果变号。

**5. [递归] 写出 $4 \times 4$ Pfaffian 沿第一行的展开式。**

??? success "参考答案"
    **公式：**
    $\operatorname{Pf}(A) = a_{12} a_{34} - a_{13} a_{24} + a_{14} a_{23}$。
    这展示了 Pfaffian 如何通过指标的完美匹配（二聚体）来覆盖所有元素。

**6. [关系] 比较 $\operatorname{perm}(A), \det(A), \operatorname{Pf}(A)$ 在计算复杂度上的区别。**

??? success "参考答案"
    **对比：**
    - $\det(A)$：$O(n^3)$，多项式级。
    - $\operatorname{perm}(A)$：$O(n 2^n)$，#P-完全（极难）。
    - $\operatorname{Pf}(A)$：$O(n^3)$（通过类似于高斯消元的方法），多项式级。
    **意义**：Pfaffian 奇迹般地在保持组合计数功能的同时，拥有与行列式同等的计算效率。

**7. [应用] 什么是 Ising 模型的“Pfaffian 方法”？**

??? success "参考答案"
    在统计力学中，二维 Ising 模型的配分函数可以表示为一个庞大反对称矩阵的 Pfaffian。
    这一发现使得物理学家能精确求解该模型，揭示了相变的代数本质。

**8. [性质] 证明：$\operatorname{Pf}(A \otimes J_2) = (\det A)$ 并不成立，正确的联系是什么？**

??? success "参考答案"
    通常涉及的是 $A \otimes J_2$ 这种形式在描述辛结构时的行列式根。
    更常用的结论是 $\operatorname{Pf}(A \oplus B) = \operatorname{Pf}(A)\operatorname{Pf}(B)$。

**9. [判定] 判定 $\begin{pmatrix} 0 & 1 & 2 & 3 \\ -1 & 0 & 1 & 2 \\ -2 & -1 & 0 & 1 \\ -3 & -2 & -1 & 0 \end{pmatrix}$ 的 Pfaffian 值。**

??? success "参考答案"
    **计算：**
    利用 4 阶公式：$a_{12}a_{34} - a_{13}a_{24} + a_{14}a_{23}$。
    $1 \cdot 1 - 2 \cdot 2 + 3 \cdot 1 = 1 - 4 + 3 = 0$。
    **验证**：该矩阵秩为 2，行列式为 0，故 $\operatorname{Pf}=0$ 正确。

**10. [应用] 为什么 Pfaffian 在描述超对称（Supersymmetry）中很重要？**

??? success "参考答案"
    超对称涉及费米子（反对称性）与玻色子（对称性）的转换。
    在计算费米子路径积分时，结果往往表现为算子反对称部分的 Pfaffian。这保证了物理概率幅在各种对称性操作下的符号一致性。

## 本章小结

Pfaffian 是反对称代数的核心算子：

1.  **平方根的优雅**：它证明了反对称系统的整体度量（行列式）具有一种更基本的、半阶的代数根，消除了开方后的符号歧义。
2.  **计数的捷径**：通过将复杂的组合匹配问题转化为多项式级别的 Pfaffian 运算，它为解决图论难题提供了最强有力的代数支点。
3.  **物理的对映**：作为费米子统计和 Ising 模型的代数载体，Pfaffian 揭示了自然界中微观对称性与宏观相变规律之间的深刻联系。
