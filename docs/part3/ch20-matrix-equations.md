# 第 20 章 矩阵方程

<div class="context-flow" markdown>

**前置**：Kronecker 积 (Ch19) · 矩阵分析 (Ch14) · 特征值 (Ch06)

**本章脉络**：线性矩阵方程概览 $\to$ 基础方程 $AX=B$ 与 $AXB=C$ $\to$ Sylvester 方程 ($AX+XB=C$) $\to$ Lyapunov 方程 ($AX+XA^T=Q$) $\to$ 解的存在性与唯一性准则（谱分离条件） $\to$ 代数 Riccati 方程 (ARE) $\to$ 连续与离散情形对比 $\to$ 数值算法初步（Bartels-Stewart 算法） $\to$ 控制理论应用（LQR, 稳定性判定）

**延伸**：矩阵方程是连接现代控制理论与数值线性代数的纽带；它是分析线性系统可观性、可控性以及求解最优控制策略的核心数学工具

</div>

当未知数本身是一个矩阵时，我们称其为矩阵方程。矩阵方程不仅是线性算子理论的延伸，更是动力系统平衡点分析、控制增益计算以及数值模拟的直接产物。本章将从最简单的线性耦合方程出发，逐步深入到复杂的非线性 Riccati 方程。

---

## 20.1 线性矩阵方程

!!! definition "定义 20.1 (基础方程)"
    1.  **左乘方程 $AX = B$**：有解当且仅当 $B$ 的列属于 $A$ 的列空间。
    2.  **双边方程 $AXB = C$**：可利用 Kronecker 积化为 $(B^T \otimes A) \operatorname{vec}(X) = \operatorname{vec}(C)$。

---

## 20.2 Sylvester 与 Lyapunov 方程

!!! definition "定义 20.2 (Sylvester 方程)"
    $$AX + XB = C$$
    其中 $A, B, C$ 为给定方阵。

!!! theorem "定理 20.1 (唯一可解性准则)"
    Sylvester 方程 $AX + XB = C$ 对任意 $C$ 有唯一解 $\iff$ $\sigma(A) \cap \sigma(-B) = \emptyset$（即 $A$ 与 $-B$ 无公共特征值）。

!!! definition "定义 20.3 (Lyapunov 方程)"
    $$AX + XA^T = Q$$
    这是 Sylvester 方程的特殊对称形式，用于判定动力系统的稳定性。若 $A$ 是稳定阵（谱在左半平面），则对任何 $Q \prec 0$，方程有唯一正定解 $X \succ 0$。

---

## 20.3 代数 Riccati 方程 (ARE)

!!! definition "定义 20.4 (连续时间 ARE)"
    $$A^T X + XA - X B R^{-1} B^T X + Q = 0$$
    这是一个二次非线性矩阵方程。
    **应用**：求解线性二次调节器 (LQR) 问题的最优反馈增益 $K = R^{-1} B^T X$。

---

## 20.4 数值算法

!!! algorithm "算法 20.1 (Bartels-Stewart 算法)"
    用于高效求解线性矩阵方程：
    1.  对 $A$ 和 $B$ 进行 Schur 分解（三角化）。
    2.  对变换后的方程进行前代/回代求解。
    3.  反变换得到原方程的解。
    **复杂度**：$O(n^3)$，比直接向量化求解（$O(n^6)$）快得多。

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

矩阵方程是高级线性系统的代数语言：

1.  **解的耦合性**：Sylvester 方程刻画了两个不同算子间的线性干涉，其谱分离条件揭示了系统共振与否的本质。
2.  **能量与稳定**：Lyapunov 方程通过二次型矩阵建立了代数与解析稳定性的桥梁，是控制工程中 Lyapunov 第二方法的算子表达。
3.  **最优性寻找**：Riccati 方程展示了非线性结构如何自然地出现在变分和最优控制问题中，其求解标志着从线性系统分析到线性系统综合的跨越。
