# 第 49 章 外代数与 Grassmannian

<div class="context-flow" markdown>

**前置**：多线性代数与张量 (Ch21) · 行列式 (Ch03) · 向量空间 (Ch04)

**本章脉络**：反对称多线性型 $\to$ 楔积 (Wedge Product) 定义与性质 $\to$ 外幂空间 $\Lambda^k(V)$ 的维数与基 $\to$ 复合矩阵 (Compound Matrices) $\to$ 柯西-比内公式 (Binet-Cauchy) 的本质 $\to$ Grassmann 流行 $Gr(k, V)$ 定义 $\to$ Plücker 嵌入与 Plücker 坐标 $\to$ 关系：外代数作为行列式的理论家园 $\to$ 应用：微分形式、组合几何

**延伸**：外代数是现代微分几何（微分形式）和拓扑学的基础；Plücker 坐标将子空间的研究转化为了射影空间中代数簇的研究，是代数几何的经典范例

</div>

行列式的代数本质是什么？为什么它是反对称的？**外代数**（Exterior Algebra，又称 Grassmann 代数）给出了终极答案。通过引入**楔积**（$\wedge$），我们不仅能形式化描述“有向体积”，还能将 $k$ 维子空间视为一个单一的代数对象。本章将揭示隐藏在行列式背后的几何代数结构。

---

## 49.1 楔积与外幂

!!! definition "定义 49.1 (楔积 $\wedge$)"
    设 $V$ 是向量空间。对于 $u, v \in V$，其**楔积** $u \wedge v$ 满足：
    1.  **双线性**：关于 $u$ 和 $v$ 均满足分配律。
    2.  **反对称性**：$v \wedge u = -(u \wedge v)$。
    3.  **零性质**：$v \wedge v = 0$。

!!! theorem "定理 49.1 (外幂空间的维数)"
    由 $V$ 中 $k$ 个向量的楔积张成的空间称为 **$k$ 阶外幂**，记作 $\Lambda^k(V)$。
    若 $\dim V = n$，则 $\dim \Lambda^k(V) = \binom{n}{k}$。
    当 $k=n$ 时，$\Lambda^n(V)$ 是 1 维的，这正是行列式的本质。

---

## 49.2 复合矩阵与 Binet-Cauchy

!!! definition "定义 49.2 (复合矩阵 $C_k(A)$)"
    对于 $m \times n$ 矩阵 $A$，其 **$k$ 阶复合矩阵** $C_k(A)$ 是一个 $\binom{m}{k} \times \binom{n}{k}$ 矩阵，其元素是 $A$ 的所有 $k$ 阶子式。

!!! theorem "定理 49.2 (算子视角)"
    复合矩阵 $C_k(A)$ 是线性变换 $A$ 在外幂空间 $\Lambda^k(V)$ 上的诱导变换。
    由此可得 **Binet-Cauchy 公式**：$C_k(AB) = C_k(A)C_k(B)$。

---

## 49.3 Grassmannian 与 Plücker 嵌入

!!! definition "定义 49.3 (Grassmannian)"
    $Gr(k, V)$ 是 $V$ 中所有 $k$ 维子空间构成的集合。它是一个紧致光滑流形。

!!! technique "Plücker 嵌入"
    将子空间 $W = \operatorname{span}\{w_1, \ldots, w_k\}$ 映为外幂空间中的一个点：
    $$\Phi(W) = [w_1 \wedge w_2 \wedge \cdots \wedge w_k] \in \mathbb{P}(\Lambda^k(V))$$
    对应的坐标称为 **Plücker 坐标**。它们必须满足一组二次关系式（Plücker 关系式）。

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

外代数是行列式与几何子空间的终极代数语言：


****：通过引入 $\wedge$ 运算，外代数将几何上的“方向”与“体积”精确地转化为代数上的正负号切换，确立了描述有向几何量的标准框架。

****：复合矩阵理论证明了线性变换在更高阶张量空间上的作用仍然具有完美的结构（Binet-Cauchy），这是研究复杂行列式恒等式的通用利器。

****：Grassmannian 与 Plücker 嵌入将动态的子空间选取转化为了静态的射影点，标志着从初等线性代数向代数几何与流形理论的跨越。
