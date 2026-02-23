# 第 16 章 正定矩阵

<div class="context-flow" markdown>

**前置**：二次型 (Ch09) · 矩阵分解 (Ch10) · 特征值 (Ch06)

**本章脉络**：正定 (PD) 与半正定 (PSD) 矩阵的定义 $\to$ 五大等价判定准则（特征值、主子式、二次型、Cholesky、Gram 矩阵） $\to$ 正定矩阵的代数性质 $\to$ Schur 补与正定性判定 $\to$ Löwner 偏序 ($\succeq$) $\to$ 矩阵平方根 $\to$ 统计学中的协方差矩阵 $\to$ 最优化中的凸性判定

**延伸**：正定矩阵是凸优化的几何核心；它不仅是标量“正数”在矩阵维度的延伸，更是现代控制理论、金融风险建模以及机器学习损失函数的基础

</div>

正定矩阵（Positive Definite Matrix）是线性代数中最受青睐的一类矩阵。它们在结构上极其对称，在谱上极其纯净（全为正），且具有完美的几何稳定性。在物理学中，它们代表了系统的能量基态；在统计学中，它们刻画了变量间的协方差；在数学中，它们是定义测度与距离的天然载体。

---

## 16.1 定义与等价判定

!!! definition "定义 16.1 (正定与半正定)"
    对于实对称矩阵 $A$：
    1.  **正定 (PD)**：对所有非零向量 $\mathbf{x}$，均有 $\mathbf{x}^T A \mathbf{x} > 0$。记作 $A \succ 0$。
    2.  **半正定 (PSD)**：对所有 $\mathbf{x}$，均有 $\mathbf{x}^T A \mathbf{x} \ge 0$。记作 $A \succeq 0$。

!!! theorem "定理 16.1 (五大判定准则)"
    对于实对称矩阵 $A$，以下条件等价于 $A \succ 0$：
    1.  **特征值**：$A$ 的所有特征值 $\lambda_i > 0$。
    2.  **主子式 (Sylvester)**：$A$ 的所有顺序主子式均大于 0。
    3.  **Cholesky**：存在唯一的对角元全正的下三角阵 $L$ 使得 $A = LL^T$。
    4.  **Gram 矩阵**：存在列满秩矩阵 $B$ 使得 $A = B^T B$。
    5.  **能量**：二次型 $Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}$ 的图像是开口向上的超抛物面。

---

## 16.2 Schur 补与正定性

!!! theorem "定理 16.2 (分块矩阵正定性)"
    设分块矩阵 $M = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$。则：
    $M \succ 0 \iff A \succ 0$ 且 **Schur 补** $S = C - B^T A^{-1} B \succ 0$。
    **应用**：这是求解带约束最优化问题和分析大型电力系统的关键技术。

---

## 16.3 Löwner 偏序

!!! definition "定义 16.2 (Löwner 偏序)"
    对于对称矩阵 $A, B$，定义 $A \succeq B \iff A - B \succeq 0$。
    **性质**：
    1.  若 $A \succeq B$ 且 $C \succeq 0$，则 $A + C \succeq B$。
    2.  若 $A \succeq B \succ 0$，则 $B^{-1} \succeq A^{-1} \succ 0$（逆算子反号）。

---

## 练习题

1. **[判定] 判定 $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$ 是否正定。**

   ??? success "参考答案"
       顺序主子式：$D_1 = 2 > 0$，$D_2 = 4 - 1 = 3 > 0$。故正定。

2. **[特征值] 若 $A$ 是正定矩阵，其行列式 $\det(A)$ 一定大于 0 吗？**

   ??? success "参考答案"
       是的。$\det(A) = \prod \lambda_i$，因为所有 $\lambda_i > 0$，故乘积必为正。

3. **[性质] 证明：若 $A \succ 0$ 且 $B \succ 0$，则 $A + B \succ 0$。**

   ??? success "参考答案"
       $\mathbf{x}^T (A+B) \mathbf{x} = \mathbf{x}^T A \mathbf{x} + \mathbf{x}^T B \mathbf{x} > 0 + 0 = 0$。

4. **[逆矩阵] 若 $A \succ 0$，证明 $A^{-1} \succ 0$。**

   ??? success "参考答案"
       $A$ 的特征值为 $\lambda_i > 0$，则 $A^{-1}$ 的特征值为 $1/\lambda_i > 0$。

5. **[Gram矩阵] 已知 $A = B^T B$，若 $B$ 不是列满秩，则 $A$ 是什么矩阵？**

   ??? success "参考答案"
       半正定矩阵（PSD）。存在非零 $x$ 使得 $Bx = 0 \implies x^T A x = \|Bx\|^2 = 0$。

6. **[对角元] 证明：正定矩阵的主对角线元素必为正。**

   ??? success "参考答案"
       取 $\mathbf{x} = \mathbf{e}_i$，则 $\mathbf{e}_i^T A \mathbf{e}_i = a_{ii} > 0$。

7. **[Schur补] 判定 $\begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}$ 是否正定（使用 Schur 补）。**

   ??? success "参考答案"
       $A = (1) \succ 0$。Schur 补 $S = 5 - 2(1)^{-1}2 = 1 > 0$。故矩阵正定。

8. **[平方根] 若 $A \succeq 0$，证明存在唯一的 $B \succeq 0$ 使得 $B^2 = A$。**

   ??? success "参考答案"
       利用对角化 $A = Q \Lambda Q^T$，取 $B = Q \Lambda^{1/2} Q^T$ 即可。

9. **[偏序] 若 $A \succeq B \succ 0$，证明 $\operatorname{tr}(A) \ge \operatorname{tr}(B)$。**

   ??? success "参考答案"
       $A - B \succeq 0 \implies \operatorname{tr}(A - B) = \sum \lambda_i(A-B) \ge 0 \implies \operatorname{tr}(A) \ge \operatorname{tr}(B)$。

10. **[应用] 为什么协方差矩阵总是半正定的？**

   ??? success "参考答案"
        协方差矩阵可以表示为 $E[(\mathbf{x}-\mu)(\mathbf{x}-\mu)^T]$，对于任何 $v$，其二次型 $v^T \Sigma v = E[(v^T(\mathbf{x}-\mu))^2] \ge 0$，代表方差是非负的。

## 本章小结

正定矩阵构建了高维空间的凸性几何：

1.  **多维的正值性**：正定性是实数“大于零”概念在算子层面的完美延伸，它确立了能量、概率和距离的合法性。
2.  **判定的多样性**：从微观的元素规律（主子式）到宏观的能量表现（二次型），再到内部的结构分解（Cholesky），多维度的等价性为不同领域的应用提供了便利工具。
3.  **计算的优越性**：正定矩阵在数值计算中具有天然的“单调稳定性”，Löwner 偏序的引入使得我们能像处理实数不等式一样处理矩阵函数，构成了算子分析的核心。
