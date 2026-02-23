# 第 38B 章 P-矩阵、H-矩阵与相关矩阵类

<div class="context-flow" markdown>

**前置**：非负矩阵与 Perron-Frobenius 理论(Ch17) · Z-矩阵与 M-矩阵(Ch38A) · 矩阵分析(Ch14) · 特征值理论(Ch6) · 线性互补问题基础

**本章脉络**：P-矩阵定义与等价刻画 $\to$ LCP 唯一可解性 $\to$ P$_0$-矩阵 $\to$ 比较矩阵与 H-矩阵 $\to$ 对角占优层次 $\to$ N-矩阵 $\to$ 半正矩阵 $\to$ Ostrowski-Reich 定理 $\to$ 迭代法与数值稳定性应用

**延伸**：P-矩阵在数学规划（LCP 理论、变分不等式）中有核心地位；H-矩阵将 M-矩阵理论推广到复矩阵和一般符号模式，广泛应用于迭代法的收敛性分析

</div>

上一章系统发展了 Z-矩阵与 M-矩阵的理论。本章转向更广泛的矩阵类。**P-矩阵**要求所有主子式为正，不限制符号模式；**H-矩阵**通过比较矩阵将 M-矩阵的思想推广到一般复矩阵。这些矩阵类在线性互补问题（LCP）和数值稳定性中有深刻的应用。

---

## 38B.1 P-矩阵

!!! definition "定义 38B.1 (P-矩阵)"
    若方阵 $A$ 的所有主子式均为正：
    $$\det(A[\alpha, \alpha]) > 0, \quad \forall\, \alpha \subseteq \{1, \dots, n\}, \alpha \ne \emptyset,$$
    则称 $A$ 为 **P-矩阵**。

!!! theorem "定理 38B.3 (LCP 的唯一可解性)"
    $A$ 是 P-矩阵，当且仅当对每个 $q \in \mathbb{R}^n$，线性互补问题 $\operatorname{LCP}(q, A)$ 有唯一解。

---

## 练习题

1. **[基础] 判定 $A = \begin{pmatrix} 1 & -2 \\ 1 & 1 \end{pmatrix}$ 是否为 P-矩阵。**
   ??? success "参考答案"
       主子式：$a_{11}=1 > 0$，$a_{22}=1 > 0$，且 $\det A = 1 - (-2) = 3 > 0$。由于所有主子式均为正，故 $A$ 是 P-矩阵。注意它不是 Z-矩阵。

2. **[转置] 证明：若 $A$ 是 P-矩阵，则 $A^T$ 也是 P-矩阵。**
   ??? success "参考答案"
       $A^T$ 的主子式是 $(A[\alpha, \alpha])^T$ 的行列式，根据行列式转置不变性，它等于 $A[\alpha, \alpha]$ 的行列式。由于后者为正，故 $A^T$ 的所有主子式也为正。

3. **[实特征值] 证明 P-矩阵的所有实特征值均为正。**
   ??? success "参考答案"
       设 $\lambda$ 为实特征值，$v$ 为对应特征向量。根据 P-矩阵性质 P2，存在指标 $i$ 使得 $v_i(Av)_i = \lambda v_i^2 > 0$。由于 $v_i^2 \ge 0$ 且 $v$ 非零，必有 $\lambda > 0$。

4. **[M-矩阵] P-矩阵类与 Z-矩阵类的交集是什么？**
   ??? success "参考答案"
       交集恰好是非奇异 M-矩阵类。

5. **[比较矩阵] 计算 $A = \begin{pmatrix} 4 & 1-i \\ -1+2i & 5 \end{pmatrix}$ 的比较矩阵 $\mathcal{M}(A)$。**
   ??? success "参考答案"
       $\mathcal{M}(A)_{ii} = |a_{ii}|$，$\mathcal{M}(A)_{ij} = -|a_{ij}|$。
       $\mathcal{M}(A) = \begin{pmatrix} 4 & -\sqrt{2} \\ -\sqrt{5} & 5 \end{pmatrix}$。

6. **[H-矩阵判据] 判定上题中的 $A$ 是否为 H-矩阵。**
   ??? success "参考答案"
       $A$ 是 H-矩阵当且仅当 $\mathcal{M}(A)$ 是 M-矩阵。计算 $\det \mathcal{M}(A) = 20 - \sqrt{10} > 0$。对角元为正且行列式为正，故其为非奇异 M-矩阵。因此 $A$ 是 H-矩阵。

7. **[稳定性] 叙述关于 SOR 迭代的 Ostrowski-Reich 定理。**
   ??? success "参考答案"
       对于对称矩阵，当且仅当 $A$ 正定时，SOR 迭代对 $\omega \in (0, 2)$ 收敛。对于 H-矩阵，SOR 迭代对 $\omega \in (0, 1]$ 保证收敛。

8. **[N-矩阵] 定义 N-矩阵。**
   ??? success "参考答案"
       所有主子式均为负的矩阵。例如 $A = \begin{pmatrix} -1 & 2 \\ 3 & -2 \end{pmatrix}$。

9. **[半正矩阵] 什么是半正矩阵（Semipositive Matrix）？**
   ??? success "参考答案"
       存在向量 $x \ge 0$ 且 $x \ne 0$，使得 $Ax > 0$。每一个非奇异 M-矩阵都是半正矩阵。

10. **[LCP示例] 求解 $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$ 且 $q = \begin{pmatrix} -1 \\ -1 \end{pmatrix}$ 的 $\operatorname{LCP}(q, A)$。**
    ??? success "参考答案"
        设 $z > 0$，解 $Az = -q \implies z = A^{-1} \begin{pmatrix} 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 & 1/2 \\ 1/2 & 1 \end{pmatrix} \begin{pmatrix} 1 \\ 1 \end{pmatrix} \frac{2}{3} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。解为 $z = (1, 1)^T, w = (0, 0)^T$。

## 本章小结

本章根据子结构的性质对矩阵进行了分类：

1. **符号无关的正性**：通过主子式定义了 P-矩阵，确立了其在 LCP 可解性中的地位。
2. **比较理论**：利用数值占优关系将 H-矩阵作为 M-矩阵的复数域推广。
3. **收敛层次**：建立了特殊矩阵类与迭代求解器数值稳定性及误差界之间的联系。
