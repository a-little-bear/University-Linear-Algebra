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

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

## 本章小结

矩阵不等式确立了线性系统的“能量边界”：

1.  **谱的稳定性**：Weyl 不等式和交错定理证明了矩阵特征值具有极强的几何惯性，微小的结构变动只能引发可控的谱漂移。
2.  **信息的极值**：Hadamard 和迹不等式揭示了矩阵在非对角化（耦合）状态下信息的损失规律，为信息论中的熵估计提供了代数上限。
3.  **分布的量化**：优序理论提供了一种比较向量“分散程度”的有力工具，揭示了矩阵对角元素与其谱之间深刻的包容关系。
