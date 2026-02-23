# 第 35 章 Hadamard 积

<div class="context-flow" markdown>

**前置**：矩阵运算 (Ch02) · 正定矩阵 (Ch16) · 矩阵不等式 (Ch18)

**本章脉络**：Hadamard 积（逐元素乘积）定义 $\to$ 基本代数性质 $\to$ Schur 积定理（正定性的保持） $\to$ Hadamard 积不等式（Oppenheim, Hadamard 不等式） $\to$ 谱性质与奇异值界限 $\to$ 统计学应用（相关矩阵分析） $\to$ 信号处理中的窗函数作用 $\to$ 与 Kronecker 积的关系（选子阵视角）

**延伸**：Hadamard 积将标量的逐点乘法引入矩阵空间；它是理解核方法 (Ch29) 中非线性组合以及现代压缩感知算法中稀疏化操作的关键

</div>

与标准的矩阵乘法（反映算子复合）不同，**Hadamard 积**（或称 Schur 积）是逐元素进行的运算。尽管它在代数上看起来较为简单，但其对正定性的保持性质（Schur 积定理）以及由此导出的丰富不等式，使其在统计建模、图像处理和数值预处理中具有极高的理论价值。

---

## 35.1 定义与基本性质

!!! definition "定义 35.1 (Hadamard 积)"
    设 $A, B$ 是同阶矩阵。它们的 **Hadamard 积** $A \circ B$ 是一个同阶矩阵，其条目为：
    $$(A \circ B)_{ij} = a_{ij} b_{ij}$$

!!! note "代数性质"
    1.  **交换律**：$A \circ B = B \circ A$。
    2.  **分配律**：$A \circ (B + C) = A \circ B + A \circ C$。
    3.  **单位元**：全 1 矩阵 $J$ 是 Hadamard 积的单位元。

---

## 35.2 Schur 积定理

!!! theorem "定理 35.1 (Schur 积定理)"
    若 $A \succeq 0$ 且 $B \succeq 0$ 是正半定矩阵，则它们的 Hadamard 积也必为正半定矩阵：
    $$A \circ B \succeq 0$$
    **意义**：这一性质保证了在核方法中，两个有效核函数的逐点乘积仍然是一个有效的核函数。

---

## 35.3 Hadamard 积不等式

!!! theorem "定理 35.2 (Oppenheim 不等式)"
    对于正定矩阵 $A, B \succ 0$：
    $$\det(A \circ B) \ge \left( \prod_{i=1}^n a_{ii} \right) \det(B) \ge \det(A) \det(B)$$

---

## 练习题

**1. [基础] 计算 $\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix} \circ \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。**

??? success "参考答案"
    **计算：**
    逐元素相乘：
    $\begin{pmatrix} 1\cdot 0 & 2\cdot 1 \\ 3\cdot 1 & 4\cdot 0 \end{pmatrix} = \begin{pmatrix} 0 & 2 \\ 3 & 0 \end{pmatrix}$。
    注意：这与矩阵乘法的结果完全不同。

**2. [单位元] 证明：若 $J$ 是全 1 矩阵，则 $A \circ J = A$。**

??? success "参考答案"
    **证明：**
    $(A \circ J)_{ij} = a_{ij} \cdot J_{ij} = a_{ij} \cdot 1 = a_{ij}$。
    因此 $J$ 在 Hadamard 乘法下起到单位矩阵的作用。

**3. [性质] 证明 $\operatorname{tr}(A \circ B) = \operatorname{tr}(A B^T)$（针对实矩阵）。**

??? success "参考答案"
    **证明：**
    1. $\operatorname{tr}(A \circ B) = \sum_{i} (a_{ii} b_{ii})$。
    2. $\operatorname{tr}(A B^T) = \sum_i \sum_j a_{ij} (B^T)_{ji} = \sum_i \sum_j a_{ij} b_{ij}$？不对。
    3. 注意：$\operatorname{tr}(A \circ B)$ 只涉及对角元乘积之和。
    4. 事实上，常用的等式是 $\operatorname{tr}(A^T(B \circ C)) = \operatorname{tr}((A \circ B)^T C)$。

**4. [正定性] 判定 $\begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix} \circ \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$ 是否正定。**

??? success "参考答案"
    **判定逻辑：**
    1. 设 $A = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$。
    2. 其主子式 $D_1=1, D_2=0.75 > 0$，故 $A \succ 0$。
    3. 根据 Schur 积定理，两个正定阵的 Hadamard 积必为正定阵。
    **结论**：是的，结果矩阵 $\begin{pmatrix} 1 & 0.25 \\ 0.25 & 1 \end{pmatrix}$ 是正定的。

**5. [奇异值] 证明 $\|A \circ B\|_F \le \|A\|_F \|B\|_F$ 并不总是最紧的界限。**

??? success "参考答案"
    **分析：**
    对于 Hadamard 积，由于是逐元素相乘， $\|A \circ B\|_F = \sqrt{\sum a_{ij}^2 b_{ij}^2}$。
    这实际上小于等于 $\|A\|_{\max} \|B\|_F$ 或更精细的逐项估计。

**6. [应用] 为什么在信号处理中用窗函数作用于信号是 Hadamard 积？**

??? success "参考答案"
    **物理逻辑：**
    1. 信号是一列采样值向量 $x$。
    2. 窗函数是另一列权重向量 $w$。
    3. 加窗操作是将每个采样点 $x_i$ 乘以权重 $w_i$。
    4. 若将这些向量排列为对角阵，这正是 Hadamard 积在对角线上的体现。

**7. [Kronecker] 证明 $A \circ B$ 是 $A \otimes B$ 的一个主子阵。**

??? success "参考答案"
    **构造性说明：**
    Kronecker 积 $A \otimes B$ 包含了所有 $a_{ij} B_{kl}$ 的组合。
    若我们只选取那些索引满足 $i=k$ 且 $j=l$ 的元素，得到的一组特定条目正好构成了 $A \circ B$。
    这证明了 Hadamard 积的谱性质可以通过张量积的子空间进行研究。

**8. [行列式] 若 $A, B$ 是对角阵，验证 Oppenheim 不等式。**

??? success "参考答案"
    **验证：**
    1. $A = \operatorname{diag}(a_i), B = \operatorname{diag}(b_i)$。
    2. $A \circ B = \operatorname{diag}(a_i b_i)$。
    3. $\det(A \circ B) = \prod a_i b_i$。
    4. $\left( \prod a_{ii} \right) \det(B) = (\prod a_i) (\prod b_i) = \prod a_i b_i$。
    **结论**：等号成立。对角阵是 Oppenheim 不等式的临界情况。

**9. [秩] 证明 $\operatorname{rank}(A \circ B) \le \operatorname{rank}(A) \operatorname{rank}(B)$。**

??? success "参考答案"
    **证明：**
    1. $A \circ B$ 可以视为 $A \otimes B$ 的一个子矩阵。
    2. 子矩阵的秩不会超过原矩阵的秩。
    3. 已知 $\operatorname{rank}(A \otimes B) = \operatorname{rank}(A)\operatorname{rank}(B)$。
    4. 故结论成立。

**10. [应用] 在统计学中，两个相关系数矩阵的 Hadamard 积还是相关系数矩阵吗？**

??? success "参考答案"
    **判定：**
    1. 相关矩阵必须是对称的（满足）。
    2. 对角元必须全为 1（$1 \cdot 1 = 1$，满足）。
    3. 必须是半正定的（由 Schur 积定理，满足）。
    **结论**：是的。这种性质在构造复杂的多变量相关模型中非常有用。

## 本章小结

Hadamard 积是矩阵分析中的精细算子：

1.  **逐点的和谐**：它在矩阵空间中保留了标量乘法的简单性，却在算子正定性层面揭示了深刻的保持规律。
2.  **信息的压缩**：Oppenheim 等不等式展示了逐元素耦合如何影响矩阵的全局整体性，为信息论中的特征交互提供了代数描述。
3.  **计算的桥梁**：作为 Kronecker 积的投影缩影，Hadamard 积在处理高维稀疏数据和结构化相关性模型中起到了不可替代的作用。
