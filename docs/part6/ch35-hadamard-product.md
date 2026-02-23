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

!!! theorem "定理 35.3 (Hadamard 不等式)"
    作为 Schur 积的一个特例，对于 $A \succ 0$：
    $$\det(A) \le \prod_{i=1}^n a_{ii}$$
    这可以看作 $A$ 与单位阵 $I$ 的某种广义 Hadamard 交互。

---

## 练习题

1. **[基础] 计算 $\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix} \circ \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。**

   ??? success "参考答案"
       $\begin{pmatrix} 0 & 2 \\ 3 & 0 \end{pmatrix}$。

2. **[Schur积] 若 $A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$，求 $A \circ A$。**

   ??? success "参考答案"
       仍然是 $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$。

3. **[性质] 证明 $\operatorname{tr}(A \circ B) = \operatorname{tr}(A^T B)$。**

   ??? success "参考答案"
       $\sum a_{ii} b_{ii} = \sum a_{ji} b_{ij}$ (对于对称阵) = $\sum (A^T)_{ij} b_{ij}$。

4. **[正定性] 判定 $\begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix} \circ \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$ 是否正定。**

   ??? success "参考答案"
       由 Schur 积定理，两个正定阵的 Hadamard 积必正定。

5. **[奇异值] 证明 $\|A \circ B\|_2 \le \|A\|_2 \|B\|_2$ 是否成立？**

   ??? success "参考答案"
       不一定。通常成立的是关于 Frobenius 范数的不等式或特定的谱界限。

6. **[应用] 为什么在信号处理中用窗函数作用于信号是 Hadamard 积？**

   ??? success "参考答案"
       窗函数是对信号在时域进行逐点的加权削减，这完全符合逐元素相乘的定义。

7. **[Kronecker] Hadamard 积 $A \circ B$ 是 $A \otimes B$ 的子矩阵吗？**

   ??? success "参考答案"
       是的。它恰好是 $A \otimes B$ 中对应于特定行和列索引的选取（主子阵）。

8. **[行列式] 若 $A, B$ 是对角阵，验证 Oppenheim 不等式。**

   ??? success "参考答案"
       此时 $A \circ B = AB$，$\det(AB) = \det(A)\det(B)$，等号成立。

9. **[rank] 证明 $\operatorname{rank}(A \circ B) \le \operatorname{rank}(A) \operatorname{rank}(B)$。**

   ??? success "参考答案"
       利用 Hadamard 积作为 Kronecker 积子阵的性质，其秩不超过母矩阵的秩。

10. **[相关矩阵] 两个相关系数矩阵的 Hadamard 积还是相关系数矩阵吗？**

   ??? success "参考答案"
        是的。它保持了对称性、主对角线全为 1 以及正半定性（由 Schur 积定理）。

## 本章小结

Hadamard 积是矩阵分析中的精细算子：

1.  **逐点的和谐**：它在矩阵空间中保留了标量乘法的简单性，却在算子正定性层面揭示了深刻的非平凡保持规律（Schur 积定理）。
2.  **信息的压缩**：Oppenheim 等不等式展示了逐元素耦合如何影响矩阵的全局整体性（如行列式），为信息论中的特征交互提供了代数描述。
3.  **计算的桥梁**：作为 Kronecker 积的投影缩影，Hadamard 积在处理高维稀疏数据和结构化相关性模型中起到了不可替代的降维作用。
