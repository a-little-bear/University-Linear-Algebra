# 第 02 章 矩阵与矩阵运算

<div class="context-flow" markdown>

**前置**：代数运算

**本章脉络**：矩阵定义与表示 → 加法与数乘 → 矩阵乘法（非交换性） → 转置运算 → 特殊矩阵（单位阵、对角阵、对称阵） → 逆矩阵 $A^{-1}$ → 矩阵幂与多项式

**延伸**：矩阵不仅是数据的表格，更是线性空间的映射算子

</div>

矩阵是线性代数的计算核心。它将向量的线性变换浓缩为一个矩形阵列。矩阵乘法的引入不仅仅是规则的堆叠，它反映了线性映射的复合。

---

## 02.1 矩阵运算规则

!!! definition "定义 02.1 (矩阵乘法)"
    若 $A$ 是 $m \times n$ 矩阵，$B$ 是 $n \times p$ 矩阵，则积 $C = AB$ 是 $m \times p$ 矩阵，其条目为：
    $$c_{ij} = \sum_{k=1}^n a_{ik} b_{kj}$$
    注意：一般情况下 $AB \neq BA$。

---

## 练习题

1. **[基础运算] 已知 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}, B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。计算 $AB$ 和 $BA$。**
   ??? success "参考答案"
       $AB = \begin{pmatrix} 1(0)+2(1) & 1(1)+2(0) \\ 3(0)+4(1) & 3(1)+4(0) \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}$。
       $BA = \begin{pmatrix} 0(1)+1(3) & 0(2)+1(4) \\ 1(1)+0(3) & 1(2)+0(4) \end{pmatrix} = \begin{pmatrix} 3 & 4 \\ 1 & 2 \end{pmatrix}$。
       可见 $AB \neq BA$。

2. **[单位矩阵] 证明对于任何 $n \times n$ 矩阵 $A$，都有 $AI = IA = A$。**
   ??? success "参考答案"
       单位矩阵 $I$ 的条目 $\delta_{ij}$ 满足：当 $i=j$ 时为 1，否则为 0。代入乘法公式得 $(AI)_{ij} = \sum a_{ik} \delta_{kj} = a_{ij}$。

3. **[转置] 已知 $(AB)^T = B^T A^T$。利用此性质求 $(A^T B)^T$。**
   ??? success "参考答案"
       $(A^T B)^T = B^T (A^T)^T = B^T A$。

4. **[对称矩阵] 若 $A$ 是对称矩阵，证明 $A^2$ 也是对称矩阵。**
   ??? success "参考答案"
       对称阵满足 $A^T = A$。则 $(A^2)^T = (AA)^T = A^T A^T = AA = A^2$。故 $A^2$ 是对称的。

5. **[逆矩阵初步] 求 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ 的逆矩阵。**
   ??? success "参考答案"
       利用 $2 \times 2$ 公式：$A^{-1} = \frac{1}{ad-bc} \begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$。
       $\det A = 1(4)-2(3) = -2$。
       $A^{-1} = \frac{1}{-2} \begin{pmatrix} 4 & -2 \\ -3 & 1 \end{pmatrix} = \begin{pmatrix} -2 & 1 \\ 1.5 & -0.5 \end{pmatrix}$。

6. **[幂运算] 计算 $A^k$ 其中 $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$。**
   ??? success "参考答案"
       $A^2 = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}, A^3 = \begin{pmatrix} 1 & 3 \\ 0 & 1 \end{pmatrix}$。
       数学归纳法可得 $A^k = \begin{pmatrix} 1 & k \\ 0 & 1 \end{pmatrix}$。

7. **[矩阵方程] 若 $AX = B$ 且 $A$ 可逆，求 $X$。**
   ??? success "参考答案"
       两边左乘 $A^{-1}$，得 $A^{-1}AX = A^{-1}B \implies X = A^{-1}B$。注意不能右乘。

8. **[秩初步] 矩阵 $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$ 的秩是多少？**
   ??? success "参考答案"
       秩为 1。它只有一个非零行（或一个线性无关的列）。

9. **[正交性初步] 验证 $\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$ 是正交矩阵（即 $AA^T = I$）。**
   ??? success "参考答案"
       $A^T = \begin{pmatrix} \cos & \sin \\ -\sin & \cos \end{pmatrix}$。乘法得对角元为 $\cos^2+\sin^2=1$，非对角元为 $\cos\sin-\sin\cos=0$。故 $AA^T = I$。

10. **[对角阵] 两个对角矩阵相乘的结果仍是对角阵吗？其条目有何规律？**
    ??? success "参考答案"
        是的。结果是对角元对应相乘的对角矩阵：$\operatorname{diag}(a_i) \operatorname{diag}(b_i) = \operatorname{diag}(a_i b_i)$。

## 本章小结

矩阵运算构建了线性代数的计算语法：

1. **非交换性**：这是矩阵乘法与标量乘法最本质的区别。
2. **结构保全**：转置和求逆操作保持了矩阵内部的代数逻辑。
3. **运算简化**：特殊矩阵（单位、对角）极大简化了复杂系统的分析。
