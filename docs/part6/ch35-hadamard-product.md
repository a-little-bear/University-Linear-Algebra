# 第 35 章 Hadamard 积与 Schur 积定理

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 正定性(Ch16) · 内积(Ch8)

**本章脉络**：Hadamard (Schur) 积定义 → 性质 → Schur 积定理 (保持正定性) → Oppenheim 不等式 → Hadamard 积与秩 → 谱范数界限 → 在核方法与统计学中的应用

**延伸**：Hadamard 积是神经网络中逐元素运算的基础，也是机器学习中构造正定核函数的核心

</div>

**Hadamard 积**（也称为 Schur 积）是同型矩阵之间的逐元素乘法。虽然标准矩阵乘法代表线性映射的复合，但 Hadamard 积产生于需要特定条目缩放的场景，如图像处理中的掩模操作或统计学中的核函数。该领域最深刻的结果是 **Schur 积定理**，它指出两个半正定矩阵的 Hadamard 积仍然是半正定的。

---

## 35.1 定义与 Schur 积定理

!!! definition "定义 35.1 (Hadamard 积)"
    $A = (a_{ij})$ 与 $B = (b_{ij})$ 的 Hadamard 积记作 $A \circ B$，其定义为：
    $$(A \circ B)_{ij} = a_{ij} b_{ij}$$

!!! theorem "定理 35.1 (Schur 积定理)"
    若 $A \succeq 0$ 且 $B \succeq 0$ 是 $n \times n$ 矩阵，则 $A \circ B \succeq 0$。

---

## 练习题

1. **[基础] 计算 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ 与 $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ 的 Hadamard 积。**
   ??? success "参考答案"
       $A \circ B = \begin{pmatrix} 1 \cdot 0 & 2 \cdot 1 \\ 3 \cdot 1 & 4 \cdot 0 \end{pmatrix} = \begin{pmatrix} 0 & 2 \\ 3 & 0 \end{pmatrix}$。

2. **[Schur定理] 利用 $A \succeq 0$ 可写为秩-1 矩阵 $v_k v_k^*$ 之和的性质证明 Schur 积定理。**
   ??? success "参考答案"
       设 $A = \sum v_k v_k^*$ 且 $B = \sum w_l w_l^*$。则 $A \circ B = \sum_{k,l} (v_k v_k^*) \circ (w_l w_l^*) = \sum_{k,l} (v_k \circ w_l) (v_k \circ w_l)^*$。由于每一项都是秩-1 的半正定阵，其和必然半正定。

3. **[Oppenheim] 叙述关于 $A, B \succeq 0$ 的 Oppenheim 不等式。**
   ??? success "参考答案"
       $\det(A \circ B) \ge \det A \cdot \prod b_{ii}$。这为 Hadamard 积的行列式提供了下界估计。

4. **[特征值] $A \circ B$ 的特征值与 $A, B$ 的特征值有何关系？**
   ??? success "参考答案"
       虽然 $\rho(A \circ B) \le \rho(A) \rho(B)$ 通常不成立，但对于半正定阵，满足 $\lambda_{\max}(A \circ B) \le \lambda_{\max}(A) \cdot \max b_{ii}$。

5. **[谱范数] 证明 $\|A \circ B\| \le \|A\| \cdot \|B\|$。**
   ??? success "参考答案"
       利用 $A \circ B$ 是 Kronecker 积 $A \otimes B$ 的主子阵这一事实。由于 $A \otimes B$ 的谱范数等于 $\|A\| \cdot \|B\|$，且取子阵不会增加谱范数，故不等式成立。

6. **[秩] $A \circ B$ 的秩最大可能是多少？**
   ??? success "参考答案"
       $\operatorname{rank}(A \circ B) \le \operatorname{rank}(A) \cdot \operatorname{rank}(B)$。

7. **[核方法] 为什么 Hadamard 积对高斯核函数很重要？**
   ??? success "参考答案"
       高斯核 $K(x, y) = \exp(- \gamma \|x-y\|^2)$ 可以视为多个简单核函数的 Hadamard 积。Schur 积定理保证了有效核函数的乘积仍然是有效且正定的核函数。

8. **[对角线] 利用 Hadamard 积表示乘积 $AB$ 的对角线。**
   ??? success "参考答案"
       $\operatorname{diag}(AB) = (A \circ B^T) \mathbf{1}$。

9. **[单位元] Hadamard 积的单位元矩阵是什么？**
   ??? success "参考答案"
       全 1 矩阵 $J$（即所有元素均为 1 的矩阵）。

10. **[迹] 证明 Schur-Hadamard 迹恒等式 $\operatorname{tr}((A \circ B)C) = \operatorname{tr}(A(B \circ C))$。**
    ??? success "参考答案"
        展开左侧：$\sum_{i,j} a_{ij} b_{ij} c_{ji}$。展开右侧：$\sum_{i,j} a_{ij} (b_{ij} c_{ij})$（假设 $B, C$ 对称）。两侧均等于元素对应乘积的总和。

## 本章小结

本章探讨了矩阵的逐元素微积分：

1. **分析特性**：确立了 Schur 积定理作为 Hadamard 积定义性的正定保持属性。
2. **行列式界限**：利用 Oppenheim 不等式约束了 Hadamard 变换后算子的空间体积。
3. **子空间动力学**：分析了秩与范数关系，将 Hadamard 积定位为 Kronecker 积的一种收缩形式。
4. **统计关联**：强调了其在协方差建模和正定核函数构造中的基础地位。
