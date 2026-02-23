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

**1. [基础] 已知 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}, B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。计算 $AB$ 和 $BA$。**

??? success "参考答案"
    **计算 $AB$：**
    $AB = \begin{pmatrix} 1\cdot 0 + 2\cdot 1 & 1\cdot 1 + 2\cdot 0 \\ 3\cdot 0 + 4\cdot 1 & 3\cdot 1 + 4\cdot 0 \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}$。
    
    **计算 $BA$：**
    $BA = \begin{pmatrix} 0\cdot 1 + 1\cdot 3 & 0\cdot 2 + 1\cdot 4 \\ 1\cdot 1 + 0\cdot 3 & 1\cdot 2 + 0\cdot 4 \end{pmatrix} = \begin{pmatrix} 3 & 4 \\ 1 & 2 \end{pmatrix}$。
    
    **结论：**
    $AB \neq BA$，这验证了矩阵乘法的非交换性。

**2. [单位矩阵] 证明对于任何 $n \times n$ 矩阵 $A$，都有 $AI = IA = A$。**

??? success "参考答案"
    **证明：**
    1. 考虑 $AI$ 的分量 $(AI)_{ij} = \sum_{k=1}^n a_{ik} \delta_{kj}$。
    2. 其中 $\delta_{kj}$ 是 Kronecker 符号，仅当 $k=j$ 时为 1，其余为 0。
    3. 因此求和式中只有 $k=j$ 的项保留，即 $(AI)_{ij} = a_{ij} \cdot 1 = a_{ij}$。
    4. 同理可证 $(IA)_{ij} = a_{ij}$。
    5. 故 $AI = IA = A$。

**3. [转置] 已知 $(AB)^T = B^T A^T$。利用此性质求 $(A^T B)^T$。**

??? success "参考答案"
    **推导：**
    1. 根据转置的分配律，积的转置等于转置之积的倒序。
    2. 令 $M = A^T, N = B$。则 $(MN)^T = N^T M^T$。
    3. 代入：$(A^T B)^T = B^T (A^T)^T$。
    4. 由于转置的转置是其自身：$(A^T)^T = A$。
    5. **结论：** $(A^T B)^T = B^T A$。

**4. [对称性] 若 $A$ 是对称矩阵，证明 $A^2$ 也是对称矩阵。**

??? success "参考答案"
    **证明：**
    1. 对称矩阵满足 $A^T = A$。
    2. 我们需要证明 $(A^2)^T = A^2$。
    3. 利用性质 $(AB)^T = B^T A^T$：$(AA)^T = A^T A^T$。
    4. 代入 $A^T = A$ 得：$(AA)^T = AA = A^2$。
    5. 因此 $A^2$ 是对称矩阵。

**5. [逆矩阵] 求 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ 的逆矩阵。**

??? success "参考答案"
    **步骤 1：计算行列式。**
    $\det(A) = 1\cdot 4 - 2\cdot 3 = 4 - 6 = -2$。
    由于 $\det(A) \neq 0$，逆矩阵存在。
    
    **步骤 2：利用 $2 \times 2$ 求逆公式。**
    对于 $\begin{pmatrix} a & b \\ c & d \end{pmatrix}$，逆为 $\frac{1}{ad-bc} \begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$。
    $A^{-1} = \frac{1}{-2} \begin{pmatrix} 4 & -2 \\ -3 & 1 \end{pmatrix} = \begin{pmatrix} -2 & 1 \\ 1.5 & -0.5 \end{pmatrix}$。

**6. [幂运算] 计算 $A^k$ 其中 $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$。**

??? success "参考答案"
    **观察前几次幂：**
    $A^2 = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$。
    $A^3 = A^2 A = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 0 & 1 \end{pmatrix}$。
    
    **归纳规律：**
    可以推测 $A^k = \begin{pmatrix} 1 & k \\ 0 & 1 \end{pmatrix}$。
    该矩阵常用于描述线性系统中的平移累积。

**7. [迹] 证明 $\operatorname{tr}(AB) = \operatorname{tr}(BA)$。**

??? success "参考答案"
    **证明：**
    1. 根据迹的定义：$\operatorname{tr}(M) = \sum M_{ii}$。
    2. $(AB)_{ii} = \sum_j a_{ij} b_{ji}$。
    3. 故 $\operatorname{tr}(AB) = \sum_i \sum_j a_{ij} b_{ji}$。
    4. $(BA)_{jj} = \sum_i b_{ji} a_{ij}$。
    5. 故 $\operatorname{tr}(BA) = \sum_j \sum_i b_{ji} a_{ij}$。
    6. 由于标量乘法满足交换律且有限项求和顺序可交换，二者相等。

**8. [初等矩阵] 左乘一个初等矩阵 $E$ 相当于对 $A$ 执行什么操作？**

??? success "参考答案"
    **结论：**
    左乘 $E$ 相当于对 $A$ 执行**一次对应的初等行变换**。
    - 若 $E$ 是交换两行的单位阵，则 $EA$ 交换 $A$ 的对应行。
    - 若 $E$ 是某行乘 $k$ 的单位阵，则 $EA$ 使 $A$ 的对应行乘 $k$。
    这一性质是高斯消元法能够被“矩阵化”的核心原因。

**9. [分块] 计算 $\begin{pmatrix} I & A \\ 0 & I \end{pmatrix} \begin{pmatrix} I & -A \\ 0 & I \end{pmatrix}$。**

??? success "参考答案"
    **分块乘法步骤：**
    结果矩阵的分块元素为：
    - 左上：$I\cdot I + A\cdot 0 = I$
    - 右上：$I\cdot(-A) + A\cdot I = -A + A = 0$
    - 左下：$0\cdot I + I\cdot 0 = 0$
    - 右下：$0\cdot(-A) + I\cdot I = I$
    
    **结论：**
    结果为单位阵 $\begin{pmatrix} I & 0 \\ 0 & I \end{pmatrix}$。这说明该矩阵的逆就是将其右上角变号。

**10. [秩初步] 矩阵 $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$ 的秩是多少？**

??? success "参考答案"
    **分析：**
    1. 矩阵的秩定义为最大线性无关行（或列）的个数。
    2. 第一行为 $(1, 0)$，非零。
    3. 第二行为 $(0, 0)$，零行。
    4. 显然只有第一行是独立贡献的。
    **结论：** 秩为 1。该矩阵将整个平面投影到了 $x$ 轴上。
