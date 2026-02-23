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

1. **[基础] 在 $\mathbb{R}^3$ 中，计算 $(e_1 + e_2) \wedge (e_2 + e_3)$。**

   ??? success "参考答案"
       $= e_1 \wedge e_2 + e_1 \wedge e_3 + e_2 \wedge e_2 + e_2 \wedge e_3 = e_1 \wedge e_2 - e_3 \wedge e_1 + e_2 \wedge e_3$。

2. **[维数] 若 $\dim V = 4$，求 $\Lambda^2(V)$ 的维数。**

   ??? success "参考答案"
       $\binom{4}{2} = 6$。

3. **[独立性] 证明 $v_1, \ldots, v_k$ 线性无关当且仅当 $v_1 \wedge \cdots \wedge v_k \neq 0$。**

   ??? success "参考答案"
       若相关，其中一个向量可由其余表示，由反对称性楔积必为 0。反之，若无关，可扩充为基，其楔积为基向量之一，非零。

4. **[行列式] 证明 $n$ 阶方阵 $A$ 的行列式满足 $Av_1 \wedge \cdots \wedge Av_n = \det(A)(v_1 \wedge \cdots \wedge v_n)$。**

   ??? success "参考答案"
       这是 $\Lambda^n(V)$ 作为一维空间的性质，任何线性变换在上面的作用都是标量乘法，该标量即为行列式。

5. **[复合矩阵] 写出 $\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ 的 2 阶复合矩阵。**

   ??? success "参考答案"
       只有一个元素，即其 2 阶子式（行列式）：$(-2)$。

6. **[Plücker] 判定向量 $(1, 0, 0, 0) \wedge (0, 1, 0, 0)$ 在 $\mathbb{P}(\Lambda^2(\mathbb{R}^4))$ 中的 Plücker 坐标。**

   ??? success "参考答案"
       在基 $\{e_i \wedge e_j\}$ 下，只有 $e_1 \wedge e_2$ 的分量为 1，其余为 0。坐标为 $(1, 0, 0, 0, 0, 0)$。

7. **[性质] 证明 $(u \wedge v) \wedge w = u \wedge (v \wedge w)$。**

   ??? success "参考答案"
       楔积满足结合律，这是由张量积的结合律及反对称化算子的性质保证的。

8. **[迹] 证明 $\operatorname{tr}(C_k(A))$ 是 $A$ 的 $k$ 阶特征值乘积之和。**

   ??? success "参考答案"
       $C_k(A)$ 的特征值是 $A$ 特征值的所有 $k$ 个组合的乘积。其迹即为这些乘积之和。

9. **[关系] 为什么说外代数是微分形式的基础？**

   ??? success "参考答案"
       微分形式是流形上每一点切空间的对偶空间的外幂空间的截面。楔积对应于形式的乘法。

10. **[应用] 什么是 Plücker 关系式？**

   ??? success "参考答案"
        是一组二次齐次方程，定义了哪些张量可以表示为简单张量（即对应于一个子空间）。对于 $Gr(2, 4)$，它是 $p_{12}p_{34} - p_{13}p_{24} + p_{14}p_{23} = 0$。

## 本章小结

外代数是行列式与几何子空间的终极代数语言：

1.  **反对称的威力**：通过引入 $\wedge$ 运算，外代数将几何上的“方向”与“体积”精确地转化为代数上的正负号切换，确立了描述有向几何量的标准框架。
2.  **算子的升维**：复合矩阵理论证明了线性变换在更高阶张量空间上的作用仍然具有完美的结构（Binet-Cauchy），这是研究复杂行列式恒等式的通用利器。
3.  **空间的对象化**：Grassmannian 与 Plücker 嵌入将动态的子空间选取转化为了静态的射影点，标志着从初等线性代数向代数几何与流形理论的跨越。
