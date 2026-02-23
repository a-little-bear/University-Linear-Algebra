# 第 36 章 矩阵稳定性与惯性

<div class="context-flow" markdown>

**前置**：特征值 (Ch06) · 矩阵分析 (Ch14) · 矩阵方程 (Ch20) · 正定矩阵 (Ch16)

**本章脉络**：稳定性的物理动机 $\to$ Hurwitz 稳定性（连续系统） $\to$ Schur 稳定性（离散系统） $\to$ Lyapunov 稳定性定理（正定矩阵判据） $\to$ 矩阵的惯性 (Inertia) 定义 $\to$ Sylvester 惯性定理推广 $\to$ Routh-Hurwitz 准则 $\to$ D-稳定性与 P-稳定性 $\to$ 应用：控制系统稳态分析、人口模型演化

**延伸**：矩阵稳定性是判定动力系统（从气象模型到经济循环）平衡点是否具备“自我恢复能力”的代数准则；它是控制理论 (Ch66) 的灵魂

</div>

如果一个物理系统受到微小扰动后能回到平衡状态，我们称该系统是稳定的。在线性模型中，这一物理属性被完全编码为系数矩阵的谱分布。本章将确立判定矩阵稳定性的代数准则，并引入“惯性”这一深刻的拓扑不变量来描述谱在复平面上的分布。

---

## 36.1 Hurwitz 与 Schur 稳定性

!!! definition "定义 36.1 (Hurwitz 稳定性)"
    方阵 $A \in M_n(\mathbb{C})$ 称为 **Hurwitz 稳定的**（或渐近稳定的），如果其所有特征值的实部均为负：
    $$\operatorname{Re}(\lambda_i) < 0, \quad \forall i$$
    **物理意义**：连续系统 $\dot{\mathbf{x}} = A\mathbf{x}$ 的解在 $t \to \infty$ 时趋于零。

!!! definition "定义 36.2 (Schur 稳定性)"
    方阵 $A$ 称为 **Schur 稳定的**，如果其所有特征值的模均小于 1：
    $$|\lambda_i| < 1, \quad \forall i$$
    **物理意义**：离散系统 $\mathbf{x}_{k+1} = A\mathbf{x}_k$ 的解在 $k \to \infty$ 时趋于零。

---

## 36.2 Lyapunov 稳定性定理

!!! theorem "定理 36.1 (Lyapunov 判据)"
    方阵 $A$ 是 Hurwitz 稳定的，当且仅当对于任何正定矩阵 $Q \succ 0$，以下 **Lyapunov 方程** 有唯一正定解 $P \succ 0$：
    $$A^T P + PA = -Q$$
    **意义**：这一结论将“谱分布”这一全局信息转化为“矩阵方程”这一代数计算，避免了显式求解特征值。

---

## 36.3 矩阵的惯性 (Inertia)

!!! definition "定义 36.3 (矩阵的惯性)"
    矩阵 $A$ 的**惯性**是一个三元组 $\operatorname{In}(A) = (\pi, \nu, \delta)$：
    - $\pi$：具有正实部的特征值个数。
    - $\nu$：具有负实部的特征值个数。
    - $\delta$：具有零实部的特征值个数。
    **性质**：$A$ 是 Hurwitz 稳定的 $\iff \operatorname{In}(A) = (0, n, 0)$。

---

## 36.4 Routh-Hurwitz 准则

!!! technique "技术：Routh 表"
    对于给定的特征多项式 $p(\lambda) = \sum a_k \lambda^k$，无需解方程，通过构造 Routh 表并检查其第一列符号变化次数，即可确定具有正实部特征值的个数。

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

矩阵稳定性是分析动态系统的核心判据：

1.  **谱的拓扑分类**：Hurwitz 和 Schur 稳定性分别定义了连续与离散演化下，“秩序”战胜“混沌”的代数边界。
2.  **能量的单调性**：Lyapunov 定理证明了稳定性本质上是某种广义能量（二次型）随时间的单调递减，将分析学问题转化为正定矩阵的代数问题。
3.  **结构的鲁棒性**：通过惯性和 Routh 表等工具，我们能在不求解精确特征值的情况下，预判系统在参数波动下的生存空间，确立了控制设计的稳健性基准。
