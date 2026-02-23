# 第 02 章 矩阵与矩阵运算

<div class="context-flow" markdown>

**前置**：线性方程组 (Ch01)

**本章脉络**：矩阵定义与表示 $\to$ 基本运算（加法、数乘、乘法） $\to$ 矩阵乘法的非交换性 $\to$ 转置运算与性质 $\to$ 特殊矩阵（单位阵、对角阵、三角阵、对称阵） $\to$ 初等矩阵与行变换的代数化 $\to$ 逆矩阵定义与性质 $\to$ 高斯-约当求逆法 $\to$ 分块矩阵运算 $\to$ 矩阵的迹

**延伸**：矩阵不仅是数据的容器，更是线性空间的映射算子；矩阵乘法的定义反映了线性变换的复合 (Ch05)

</div>

如果说线性方程组是线性代数的语言，那么矩阵就是它的符号系统。矩阵将复杂的线性关系浓缩为简洁的矩形阵列，并赋予了其一套精妙的代数规则。本章将确立矩阵运算的标准公理，并深入探讨逆矩阵这一核心工具。

---

## 02.1 矩阵的基本定义与运算

!!! definition "定义 02.1 (矩阵)"
    一个 $m \times n$ 的**矩阵**（Matrix）是由 $m \cdot n$ 个元素排成的 $m$ 行 $n$ 列的矩形阵列。常用大写字母 $A, B$ 表示。

!!! definition "定义 02.2 (矩阵乘法)"
    若 $A$ 是 $m \times n$ 矩阵，$B$ 是 $n \times p$ 矩阵，则积 $C = AB$ 是 $m \times p$ 矩阵，其条目为：
    $$c_{ij} = \sum_{k=1}^n a_{ik} b_{kj}$$
    **警告**：矩阵乘法一般不满足交换律，即 $AB \neq BA$。

---

## 02.2 特殊矩阵类

!!! definition "定义 02.3 (特殊矩阵)"
    1.  **单位矩阵 $I$**：对角线全为 1，其余为 0。满足 $AI = IA = A$。
    2.  **对称矩阵**：满足 $A^T = A$。
    3.  **反对称矩阵**：满足 $A^T = -A$。
    4.  **三角矩阵**：上三角矩阵（主对角线下方全为 0）或下三角矩阵。

---

## 02.3 初等矩阵与逆矩阵

!!! definition "定义 02.4 (逆矩阵)"
    对于 $n$ 阶方阵 $A$，若存在方阵 $B$ 使得 $AB = BA = I$，则称 $A$ 是**可逆的**（或非奇异的），$B$ 称为 $A$ 的**逆矩阵**，记作 $A^{-1}$。

!!! theorem "定理 02.1 (逆矩阵的性质)"
    1.  $(A^{-1})^{-1} = A$
    2.  $(AB)^{-1} = B^{-1} A^{-1}$（穿脱法则）
    3.  $(A^T)^{-1} = (A^{-1})^T$

!!! algorithm "算法 02.1 (高斯-约当求逆法)"
    构造分块矩阵 $[A | I]$，通过初等行变换将其左侧化为 $I$，则右侧即为 $A^{-1}$：
    $$[A | I] \xrightarrow{\text{row operations}} [I | A^{-1}]$$

---

## 02.4 分块矩阵

!!! technique "技术：分块运算"
    对于大规模矩阵，常将其划分为较小的子块。若分块尺寸匹配，分块矩阵的加法和乘法规则与普通矩阵完全一致。这在分布式计算和稀疏矩阵处理中至关重要。

---

## 练习题

1. **[基础] 已知 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}, B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。计算 $AB$ 和 $BA$。**
   ??? success "参考答案"
       $AB = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}, BA = \begin{pmatrix} 3 & 4 \\ 1 & 2 \end{pmatrix}$。可见 $AB \neq BA$。

2. **[单位矩阵] 证明对于任何 $n \times n$ 矩阵 $A$，都有 $AI = IA = A$。**
   ??? success "参考答案"
       利用乘法定义：$(AI)_{ij} = \sum a_{ik} \delta_{kj}$。由于 $\delta_{kj}$ 仅在 $k=j$ 时为 1，故结果为 $a_{ij}$。

3. **[转置] 已知 $(AB)^T = B^T A^T$。利用此性质求 $(A^T B)^T$。**
   ??? success "参考答案"
       $(A^T B)^T = B^T (A^T)^T = B^T A$。

4. **[对称性] 若 $A$ 是对称矩阵，证明 $A^2$ 也是对称矩阵。**
   ??? success "参考答案"
       $(A^2)^T = (AA)^T = A^T A^T = AA = A^2$。

5. **[逆矩阵] 求 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ 的逆矩阵。**
   ??? success "参考答案"
       $\det(A) = 4-6 = -2$。$A^{-1} = \frac{1}{-2} \begin{pmatrix} 4 & -2 \\ -3 & 1 \end{pmatrix} = \begin{pmatrix} -2 & 1 \\ 1.5 & -0.5 \end{pmatrix}$。

6. **[幂运算] 计算 $A^k$ 其中 $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$。**
   ??? success "参考答案"
       $A^2 = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}, A^3 = \begin{pmatrix} 1 & 3 \\ 0 & 1 \end{pmatrix}$。归纳得 $A^k = \begin{pmatrix} 1 & k \\ 0 & 1 \end{pmatrix}$。

7. **[迹] 证明 $\operatorname{tr}(AB) = \operatorname{tr}(BA)$。**
   ??? success "参考答案"
       $\operatorname{tr}(AB) = \sum_i \sum_k a_{ik}b_{ki} = \sum_k \sum_i b_{ki}a_{ik} = \operatorname{tr}(BA)$。

8. **[初等矩阵] 左乘一个初等矩阵 $E$ 相当于对 $A$ 执行什么操作？**
   ??? success "参考答案"
       相当于对 $A$ 执行对应的初等行变换。

9. **[分块] 计算 $\begin{pmatrix} I & A \\ 0 & I \end{pmatrix} \begin{pmatrix} I & -A \\ 0 & I \end{pmatrix}$。**
   ??? success "参考答案"
       $\begin{pmatrix} I & -A+A \\ 0 & I \end{pmatrix} = \begin{pmatrix} I & 0 \\ 0 & I \end{pmatrix}$。这说明该矩阵的逆是将其右上角变号。

10. **[秩初步] 矩阵 $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$ 的秩是多少？**
    ??? success "参考答案"
        秩为 1。它只有一个非零行。

## 本章小结

矩阵运算构建了线性代数的计算语法：

1.  **非交换性**：这是矩阵乘法与标量乘法最本质的区别，决定了算子作用的顺序不可随意调换。
2.  **结构保全**：转置和求逆操作保持了矩阵内部的逻辑一致性，是解决算子方程的基础。
3.  **计算简化**：特殊矩阵（单位、对角、分块）极大简化了复杂系统的分析，是数值算法优化的核心切入点。
