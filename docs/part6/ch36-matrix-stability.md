# 第 36 章 矩阵稳定性与惯性

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵指数(Ch13) · 微分方程(Ch26) · Lyapunov 方程(Ch20)

**本章脉络**：连续系统稳定性 (Hurwitz) → 离散系统稳定性 (Schur) → Lyapunov 稳定性定理 → Lyapunov 方程 $A^T P + PA = -Q$ → Routh-Hurwitz 准则 → Jury 稳定性判据 → 稳定性半径 → 鲁棒稳定性

**延伸**：矩阵稳定性是控制理论(Ch66)的代数前提，也是分析经济、生物系统中平衡点动态行为的基础

</div>

矩阵稳定性关注线性动力系统轨迹的长期行为。如果一个矩阵的所有特征值均位于复平面的左半开平面，则称其为 **Hurwitz 稳定的**，这保证了微分方程 $\dot{x} = Ax$ 的解会衰减到零。如果所有特征值均位于单位圆内部，则称其为 **Schur 稳定的**，这保证了差分方程 $x_{k+1} = Ax_k$ 的稳定性。**Lyapunov 方程**提供了一种无需显式计算特征值，而是通过正定矩阵来验证稳定性的代数方法。

---

## 36.1 Hurwitz 与 Schur 稳定性

!!! definition "定义 36.1 (Hurwitz 矩阵)"
    若 $n \times n$ 复矩阵 $A$ 的所有特征值 $\lambda_i$ 均满足 $\operatorname{Re}(\lambda_i) < 0$，则称 $A$ 是 **Hurwitz 稳定的**。

!!! definition "定义 36.2 (Schur 矩阵)"
    若 $n \times n$ 复矩阵 $A$ 的所有特征值 $\lambda_i$ 均满足 $|\lambda_i| < 1$，则称 $A$ 是 **Schur 稳定的**。

!!! theorem "定理 36.1 (Lyapunov 稳定性定理)"
    $A$ 是 Hurwitz 稳定的，当且仅当对于任何给定的 $Q \succ 0$，都存在唯一的 $P \succ 0$ 满足 Lyapunov 方程：
    $$A^T P + PA = -Q$$

---

## 练习题

1. **[基础] 判定 $A = \begin{pmatrix} -1 & 10 \\ 0 & -2 \end{pmatrix}$ 是否为 Hurwitz 稳定的。**
   ??? success "参考答案"
       是的。上三角矩阵的特征值即为其对角元：$\{-1, -2\}$。由于实部均小于 0，故 $A$ 是 Hurwitz 稳定的。

2. **[Schur稳定性] 判定 $A = \begin{pmatrix} 0.5 & 0.5 \\ 0.5 & 0.5 \end{pmatrix}$ 是否为 Schur 稳定的。**
   ??? success "参考答案"
       特征方程为 $\lambda^2 - \lambda = 0$，解得特征值为 $\{1, 0\}$。因为存在特征值在单位圆上（$|1|=1$），故 $A$ 不是 Schur 稳定的（属于临界稳定）。

3. **[Lyapunov方程] 为什么 Lyapunov 方程中的 $P$ 必须是正定的？**
   ??? success "参考答案"
       二次型 $V(x) = x^T P x$ 充当了系统的“能量函数”。$P \succ 0$ 确保了当 $x \neq 0$ 时能量为正，且随 $\|x\|$ 趋于无穷而趋于无穷。方程 $\dot{V} = -x^T Q x < 0$ 保证了能量沿系统轨迹严格递减，直到系统回到原点。

4. **[迹与行列式] $\operatorname{tr}(A)$ 和 $\det(A)$ 对 Hurwitz 稳定性有何暗示？**
   ??? success "参考答案"
       对于 Hurwitz 矩阵，必须满足 $\operatorname{tr}(A) = \sum \operatorname{Re}(\lambda_i) < 0$ 且 $(-1)^n \det A > 0$。这些是稳定的必要条件，但对于 $n \ge 3$ 的情况并非充分条件。

5. **[Routh-Hurwitz] 对于 $2 \times 2$ 实矩阵，写出其 Hurwitz 稳定的充要条件。**
   ??? success "参考答案"
       设特征多项式为 $\lambda^2 + a_1 \lambda + a_0$。稳定要求 $a_1 > 0$ 且 $a_0 > 0$。对应矩阵参数即 $\operatorname{tr}(A) < 0$ 且 $\det(A) > 0$。

6. **[双线性变换] 如何将 Hurwitz 稳定性问题转化为 Schur 稳定性问题？**
   ??? success "参考答案"
       利用 Cayley 变换（或称双线性变换） $s = \frac{z-1}{z+1}$。该映射将 $s$ 平面的左半部分映射到 $z$ 平面的单位圆内部。

7. **[正矩阵] 证明 Metzler 矩阵（非对角元 $\ge 0$）是 Hurwitz 的，当且仅当存在正向量 $d > 0$ 满足 $Ad < 0$。**
   ??? success "参考答案"
       这是 M-矩阵的基本性质。对于正系统，稳定性可以通过线性 Lyapunov 函数 $V(x) = d^T x$ 来验证，这比通用的二次型判据更简单。

8. **[稳定性半径] 定义复稳定性半径 $r_{\mathbb{C}}(A)$。**
   ??? success "参考答案"
       $r_{\mathbb{C}}(A) = \min \{ \|\Delta\| : A+\Delta \text{ 不稳定} \}$。根据小增益定理，它等于 $1 / \sup_{\omega} \|(i\omega I - A)^{-1}\|_2$。

9. **[对易族] 若 $A$ 和 $B$ 均 Hurwitz 稳定且 $AB=BA$，证明 $A+B$ 也是稳定的。**
   ??? success "参考答案"
       由于 $A, B$ 对易，它们可以同时三角化。$A+B$ 的特征值集合为 $\{\lambda_i(A) + \lambda_i(B)\}$。两个实部为负的复数之和，其实部必然仍为负。

10. **[离散Lyapunov] 写出用于判定 Schur 稳定性的 Lyapunov 方程形式。**
    ??? success "参考答案"
        $A^T P A - P = -Q$。若 $Q \succ 0$，则 $P \succ 0$ 存在的充要条件是 $A$ 为 Schur 稳定的。这捕捉了离散时间步之间的能量差：$V(x_{k+1}) - V(x_k) = -x_k^T Q x_k$。

## 本章小结

本章确立了线性系统渐近收敛的代数判据：

1. **谱域界定**：根据特征值相对于虚轴和单位圆的位置，定义了 Hurwitz 和 Schur 稳定性。
2. **Lyapunov 方法**：将稳定性转化为线性矩阵方程的正定解存在性问题，避免了直接计算特征值。
3. **多项式准则**：探讨了从特征多项式系数直接判定稳定性的 Routh-Hurwitz 和 Jury 检验。
4. **鲁棒性度量**：引入稳定性半径量化了系统在遭受模型摄动时维持稳定性的边界距离。
