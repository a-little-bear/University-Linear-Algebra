# 第 18 章 矩阵不等式

<div class="context-flow" markdown>

**前置**：正定矩阵 (Ch16) · 范数 (Ch15) · 奇异值分解 (Ch11)

**本章脉络**：从标量不等式到算子不等式 $\to$ 特征值不等式（Weyl 不等式、交错定理） $\to$ 行列式不等式（Hadamard, Fischer 不等式） $\to$ 迹不等式（von Neumann, Golden-Thompson） $\to$ 奇异值不等式（Ky Fan 范数） $\to$ 优序理论 (Majorization) $\to$ 算子单调与算子凸函数初步 (Ch46)

**延伸**：矩阵不等式是信息论（熵的凹性）、压缩感知以及量子力学中不确定性原理的数学支柱；它将精确的“等于”转化为受限的“包含”或“界限”

</div>

矩阵不等式是矩阵分析中最精妙的分支。它研究的不是矩阵的精确值，而是矩阵性质（如特征值、奇异值、迹）之间的制约关系。正如实数轴上的不等式刻画了量的相对大小，矩阵不等式刻画了算子能量和信息量的相对分布。

---

## 18.1 特征值不等式

!!! theorem "定理 18.1 (Weyl 不等式)"
    设 $A, B$ 是 Hermite 矩阵，$C = A + B$。设特征值按降序排列，则对所有 $j+k-1 \le n$：
    $$\lambda_{j+k-1}(A+B) \le \lambda_j(A) + \lambda_k(B)$$
    **物理意义**：扰动对系统特征值（能级）的影响受到扰动矩阵谱规模的严格限制。

!!! theorem "定理 18.2 (Cauchy 交错定理)"
    设 $B$ 是 $n$ 阶 Hermite 矩阵 $A$ 的 $n-1$ 阶主子阵。则：
    $$\lambda_1(A) \ge \lambda_1(B) \ge \lambda_2(A) \ge \lambda_2(B) \ge \cdots \ge \lambda_{n-1}(B) \ge \lambda_n(A)$$

---

## 18.2 行列式与迹不等式

!!! theorem "定理 18.3 (Hadamard 不等式)"
    对于任意正定矩阵 $A \succ 0$：
    $$\det(A) \le \prod_{i=1}^n a_{ii}$$
    等号成立当且仅当 $A$ 为对角阵。
    **几何解释**：平行多面体的体积小于等于其棱长乘积（仅在正交时取等）。

!!! theorem "定理 18.4 (Golden-Thompson 不等式)"
    对于 Hermite 矩阵 $A, B$：
    $$\operatorname{tr}(e^{A+B}) \le \operatorname{tr}(e^A e^B)$$
    这是量子统计力学中的核心不等式。

---

## 18.3 优序理论 (Majorization)

!!! definition "定义 18.1 (优序 $\prec$)"
    设 $x, y \in \mathbb{R}^n$。若 $\sum_{i=1}^k x_{(i)} \le \sum_{i=1}^k y_{(i)}$ 对 $k=1,\ldots,n-1$ 成立且总和相等，则称 $y$ **优于** $x$，记作 $x \prec y$。
    **Schur-Horn 定理**：Hermite 矩阵的对角元向量被其特征值向量优序：$diag(A) \prec \lambda(A)$。

---

## 练习题

1. **[Weyl] 已知 $\|E\|_2 = 0.1$，若 $A$ 的特征值为 5，则 $A+E$ 的特征值在什么区间？**
   ??? success "参考答案"
       由 Weyl 不等式 $|\lambda_i(A+E) - \lambda_i(A)| \le \|E\|_2$，故特征值位于 $[4.9, 5.1]$。

2. **[Hadamard] 计算 $\begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ 的行列式并验证 Hadamard 不等式。**
   ??? success "参考答案"
       $\det = 3$。对角元乘积 $2 \cdot 2 = 4$。验证 $3 \le 4$。

3. **[交错] 若 $3 \times 3$ 阵特征值为 10, 5, 1，其 $2 \times 2$ 主子阵的最大特征值可能为 12 吗？**
   ??? success "参考答案"
       不可能。由交错定理 $\lambda_1(B) \le \lambda_1(A) = 10$。

4. **[迹] 证明对于正定阵 $A, B$，$\operatorname{tr}(AB) \le \operatorname{tr}(A)\operatorname{tr}(B)$。**
   ??? success "参考答案"
       $\operatorname{tr}(AB) = \sum \lambda_i(AB) \le \sum \sigma_i(A)\sigma_i(B) \le (\sum \sigma_i(A))(\sum \sigma_i(B)) = \operatorname{tr}(A)\operatorname{tr}(B)$。

5. **[优序] 判定向量 $(1, 1)$ 与 $(2, 0)$ 的优序关系。**
   ??? success "参考答案"
       $(1, 1) \prec (2, 0)$。因为 $1 < 2$ 且 $1+1 = 2+0$。

6. **[Fischer] 叙述 Fischer 不等式。**
   ??? success "参考答案"
       对于分块正定阵 $\begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$，$\det \begin{pmatrix} A & B \\ B^T & C \end{pmatrix} \le \det(A)\det(C)$。

7. **[Ky Fan] 什么是 Ky Fan $k$-范数？**
   ??? success "参考答案"
       前 $k$ 个最大奇异值之和：$\|A\|_{(k)} = \sum_{i=1}^k \sigma_i(A)$。

8. **[算术几何] 证明 $\det(A)^{1/n} \le \frac{1}{n} \operatorname{tr}(A)$ 对 $A \succ 0$ 成立。**
   ??? success "参考答案"
       这是特征值的算术-几何平均不等式：$(\prod \lambda_i)^{1/n} \le \frac{1}{n} \sum \lambda_i$。

9. **[凸性] 证明映射 $A \mapsto \log \det A$ 在正定锥上是凹的。**
   ??? success "参考答案"
       这等价于对 $\det(\lambda A + (1-\lambda)B) \ge (\det A)^\lambda (\det B)^{1-\lambda}$ 的验证，是 Brunn-Minkowski 不等式的矩阵版。

10. **[应用] 矩阵不等式在量子信息中有什么用？**
    ??? success "参考答案"
        用于证明量子熵的强次可加性，确立量子通信容量的上限。

## 本章小结

矩阵不等式确立了线性系统的“能量边界”：

1.  **谱的稳定性**：Weyl 不等式和交错定理证明了矩阵特征值具有极强的几何惯性，微小的结构变动只能引发可控的谱漂移。
2.  **信息的极值**：Hadamard 和迹不等式揭示了矩阵在非对角化（耦合）状态下信息的损失规律，为信息论中的熵估计提供了代数上限。
3.  **分布的量化**：优序理论提供了一种比较向量“分散程度”的有力工具，揭示了矩阵对角元素与其谱之间深刻的包容关系。
