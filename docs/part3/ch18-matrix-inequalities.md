# 第 18 章 矩阵不等式

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 特征值(Ch6) · 范数(Ch15) · 矩阵分析(Ch14)

**本章脉络**：Löwner 偏序 $(\succeq)$ → 正定矩阵的不等式 (Weyl, Ky Fan) → 奇异值不等式 → 范数不等式 (Hadamard, Cauchy-Schwarz) → 迹不等式 (Golden-Thompson) → 特征值之和的变分刻画 → 算子单调与算子凸函数

**延伸**：矩阵不等式是凸优化、信息论以及量子计算的分析基础，它将标量序推广到了复杂的算子流形上

</div>

矩阵不等式研究方阵之间基于半正定性质的偏序关系。与标量不等式相比，由于矩阵的非对角项存在，不等式的建立往往需要通过谱分解或变分原理来实现。

---

## 18.1 Löwner 序与谱不等式

!!! definition "定义 18.1 (Löwner 偏序)"
    设 $A, B$ 是 Hermitian 矩阵。称 $A$ 在 Löwner 序下大于等于 $B$（记作 $A \succeq B$），若 $A - B$ 是半正定矩阵。

!!! theorem "定理 18.3 (Weyl 不等式)"
    设 $\lambda_i(A)$ 为按降序排列的特征值。对于 Hermitian 矩阵 $A, B$：
    $$\lambda_{i+j-1}(A+B) \le \lambda_i(A) + \lambda_j(B)$$

---

## 练习题

1. **[基础性质] 证明：若 $A \succeq B \succeq 0$，则对任何矩阵 $C$，都有 $C^* A C \succeq C^* B C$。**
   ??? success "参考答案"
       考虑二次型：$x^* (C^* A C - C^* B C) x = (Cx)^* (A-B) (Cx)$。
       令 $y = Cx$。由于 $A-B \succeq 0$，故 $y^* (A-B) y \ge 0$。
       因为对所有 $x$ 成立，故 $C^* A C \succeq C^* B C$。

2. **[特征值单调性] 若 $A \succeq B$，证明 $\lambda_i(A) \ge \lambda_i(B)$ 对所有 $i$ 成立。**
   ??? success "参考答案"
       这是由特征值的 Courant-Fischer 变分刻画（极小极大原理）得出的。因为对任意子空间，$x^T Ax \ge x^T Bx$，其最大值/最小值也必然满足相应的序关系。

3. **[求逆] 设 $A \succeq B \succ 0$，证明 $B^{-1} \succeq A^{-1}$。**
   ??? success "参考答案"
       利用合同变换。$A \succeq B \implies B^{-1/2} A B^{-1/2} \succeq I$。
       设 $B^{-1/2} A B^{-1/2} = X$，则 $X \succeq I \implies X^{-1} \preceq I$。
       即 $(B^{-1/2} A B^{-1/2})^{-1} \preceq I \implies B^{1/2} A^{-1} B^{1/2} \preceq I$。
       两侧乘以 $B^{-1/2}$ 得 $A^{-1} \preceq B^{-1}$。

4. **[Hadamard] 叙述 Hadamard 不等式并给出几何解释。**
   ??? success "参考答案"
       $\det(A) \le \prod a_{ii}$（对正定阵）。几何上，这意味着由列向量构成的平行多面体的体积，在列向量互相正交时（对角阵）达到最大。

5. **[迹不等式] 证明 $\operatorname{tr}(AB) \le \operatorname{tr}(A) \operatorname{tr}(B)$ 对 $A, B \succeq 0$ 成立吗？**
   ??? success "参考答案"
       不成立。反例：$A = B = I_{2 \times 2}$。$\operatorname{tr}(I^2) = 2$，而 $\operatorname{tr}(I)\operatorname{tr}(I) = 4$。虽然 $2 \le 4$ 成立，但这不是通用结论。正确的结论通常是 $\operatorname{tr}(AB) \le \lambda_{\max}(A) \operatorname{tr}(B)$。

6. **[Ky Fan] 叙述 Ky Fan $k$-范数之和不等式。**
   ??? success "参考答案"
       $\sum_{i=1}^k \lambda_i(A+B) \le \sum_{i=1}^k \lambda_i(A) + \sum_{i=1}^k \lambda_i(B)$。这反映了特征值和在加法下的次可加性。

7. **[Golden-Thompson] 叙述 Golden-Thompson 不等式。**
   ??? success "参考答案"
       $\operatorname{tr}(e^{A+B}) \le \operatorname{tr}(e^A e^B)$ 对所有 Hermitian 矩阵 $A, B$ 成立。这在统计力学中衡量了自由能的界限。

8. **[算子单调] 举出一个在标量域递增但在算子序下不单调的函数。**
   ??? success "参考答案"
       $f(t) = t^2$。尽管当 $a > b > 0$ 时 $a^2 > b^2$，但存在 $A \succeq B \succeq 0$ 使得 $A^2 \nsucceq B^2$。只有如 $\sqrt{t}, \log t, 1/t$ 等特定的“算子单调函数”能保持 Löwner 序。

9. **[Fiedler] 叙述关于对称阵乘积特征值的 Fiedler 不等式。**
   ??? success "参考答案"
       $\lambda(A+B)$ 被 $\lambda(A) + \lambda(B)$ 在 Majorization（优超）意义下约束。

10. **[应用] 矩阵不等式在压缩感知 (Compressed Sensing) 中有何作用？**
    ??? success "参考答案"
        用于证明受限等距性质 (RIP)。通过约束观察矩阵的奇异值波动范围，保证了在降维映射后，稀疏信号的能量能够被基本保留，从而实现精确重构。

## 本章小结

矩阵不等式是线性代数的“软约束”理论：

1. **序的推广**：Löwner 序将线性比较从直线推向了半正定锥。
2. **变分本质**：所有的谱不等式背后都隐藏着极值优化问题的影子。
3. **分析刚性**：不等式确立了矩阵在加法、乘法和函数作用下的结构演化边界。
