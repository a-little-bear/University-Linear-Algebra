# 第 49 章 外代数与 Grassmannian

<div class="context-flow" markdown>

**前置**：向量空间 (Ch4) · 线性变换 (Ch5) · 多线性代数 (Ch21) · 行列式 (Ch3)

**本章脉络**：外积（楔积）→ 外幂空间 $\Lambda^k(V)$ → 交替多线性映射 → 行列式作为顶外幂 → 复合矩阵 → Grassmannian → Plücker 坐标 → 外代数的万有性质

**延伸**：外代数是微分形式（现代微分几何和物理学的语言）的代数基础；Grassmannian 是代数几何和优化（子空间跟踪）中的核心对象；复合矩阵给出了特征值乘积的不等式

</div>

向量空间中的标量积 $\langle u, v \rangle$ 衡量两个向量的"对齐程度"。但很多几何对象——面积、体积、通量——本质上需要一种**反对称**的乘法运算。两个向量 $u, v$ 张成的平行四边形的（有向）面积，不能用点积或矩阵乘法自然表达，却恰好由**楔积**（wedge product）$u \wedge v$ 完美捕捉。

---

## 49.1 外积与外幂空间

!!! definition "定义 49.1 (外积 / 楔积)"
    楔积 $u \wedge v$ 满足：
    1. **双线性**；
    2. **反对称性**：$u \wedge u = 0$，对所有 $u \in V$。
    推出：$u \wedge v = -(v \wedge u)$。

!!! theorem "定理 49.1 (维数)"
    $n$ 维向量空间 $V$ 的 $k$ 次外幂空间 $\Lambda^k(V)$ 的维数为 $\binom{n}{k}$。其基元素为 $e_{i_1} \wedge \dots \wedge e_{i_k}$ ($i_1 < \dots < i_k$)。

---

## 49.2 行列式的外代数定义

!!! theorem "定理 49.3"
    线性变换 $T: V \to V$ 在唯一的顶外幂 $\Lambda^n(V)$（维数为 1）上的作用标量即为 $\det(T)$：
    $$T(v_1) \wedge \dots \wedge T(v_n) = \det(T) (v_1 \wedge \dots \wedge v_n)$$

---

## 49.3 Grassmannian 与 Plücker 嵌入

!!! definition "定义 49.8 (Grassmannian)"
    $\operatorname{Gr}(k, V)$ 是 $V$ 中所有 $k$ 维子空间的集合。通过 Plücker 嵌入：
    $$\iota: \operatorname{Gr}(k, V) \to \mathbb{P}(\Lambda^k(V)), \quad \langle v_1, \dots, v_k \rangle \mapsto [v_1 \wedge \dots \wedge v_k]$$
    Grassmannian 成为射影空间中的代数簇。其维数为 $k(n-k)$。

---

## 练习题

1. **[几何] 两个向量 $u, v \in \mathbb{R}^3$ 的外积 $u \wedge v$ 的维数是多少？其分量与叉积 $u \times v$ 有什么联系？**
   ??? success "参考答案"
       维数为 $\binom{3}{2} = 3$。在 $\mathbb{R}^3$ 中，外积的分量 $(u_1 v_2 - u_2 v_1, u_1 v_3 - u_3 v_1, u_2 v_3 - u_3 v_2)$ 恰好对应于叉积的分量（符号可能因约定而异）。这种对应是通过 Hodge 对偶 $\star$ 实现的。

2. **[代数] 证明：$u \wedge v = 0$ 当且仅当向量 $u$ 和 $v$ 线性相关。**
   ??? success "参考答案"
       若 $u, v$ 相关，则 $v = cu$，故 $u \wedge v = c(u \wedge u) = 0$。若无关，则 $\{u, v\}$ 可扩充为基，由外代数基的性质，$u \wedge v$ 是基元素之一，非零。

3. **[计算] 计算 $3 \times 3$ 矩阵 $A = \operatorname{diag}(1, 2, 3)$ 的二次复合矩阵 $C_2(A)$。**
   ??? success "参考答案"
       $C_2(A)$ 是 $3 \times 3$ 矩阵，特征值是 $A$ 特征值的两两乘积：$\{1 \cdot 2, 1 \cdot 3, 2 \cdot 3\} = \{2, 3, 6\}$。故 $C_2(A) = \operatorname{diag}(2, 3, 6)$。

4. **[性质] 为什么在 $n$ 维空间中，任何 $n+1$ 个向量的楔积必为零？**
   ??? success "参考答案"
       因为外代数是交替的。$n+1$ 个向量在 $n$ 维空间中必然线性相关。其中一个向量可表示为其余向量的组合，代入楔积后会出现重复项，根据 $v \wedge v = 0$ 规律，整体必为 0。

5. **[可分解性] 判定 $\omega = e_1 \wedge e_2 + e_3 \wedge e_4 \in \Lambda^2(\mathbb{R}^4)$ 是否可分解为两个向量的积？**
   ??? success "参考答案"
       不可分解。计算 $\omega \wedge \omega = 2 e_1 \wedge e_2 \wedge e_3 \wedge e_4 \neq 0$。而对于任何可分解形式 $\alpha = u \wedge v$，必有 $\alpha \wedge \alpha = 0$。

6. **[Plücker] 什么是 Plücker 嵌入？它将 $k$ 维子空间映到了哪里？**
   ??? success "参考答案"
       Plücker 嵌入将 Grassmannian $\operatorname{Gr}(k, n)$（所有 $k$ 维子空间的集合）映入射影空间 $\mathbb{P}(\Lambda^k V)$。其映像是所有“可分解”的 $k$-向量构成的射影子簇。

7. **[行列式] 证明：$\det(AB) = \det(A)\det(B)$（利用外代数视角）。**
   ??? success "参考答案"
       线性变换 $T$ 的行列式被定义为它在顶外幂 $\Lambda^n V$ 上的作用标量。由于 $\Lambda^n(AB) = \Lambda^n(A) \circ \Lambda^n(B)$，作用在基元素上的标量自然满足乘法性质。

8. **[Grassmannian] 计算 $\operatorname{Gr}(2, 4)$ 的流形维数。**
   ??? success "参考答案"
       $\dim = k(n-k) = 2(4-2) = 4$。这对应于 Plücker 空间 $\mathbb{P}^5$ 中的一个 4 维二次超曲面。

9. **[缩并] 向量对 $k$-形式的缩并（Contraction）$\iota_v$ 具有什么代数性质？**
   ??? success "参考答案"
       它是一个次数为 -1 的反导子（Antiderivation），满足 $\iota_v^2 = 0$。它在微分几何中用于将向量场与微分形式结合产生低阶形式。

10. **[爱因斯坦思考题] 爱因斯坦的微分几何中，曲率张量可以通过微分形式的楔积来表达。为什么说“外代数是体积与定向的自然语言”，而普通的矩阵乘法不是？**
    ??? success "参考答案"
        矩阵乘法通常不区分“长度”与“面积”。而外代数通过反对称性（$u \wedge v = -v \wedge u$）内生了“定向”的概念。楔积的系数正好对应于子空间在各坐标平面上的投影面积。这种结构使得物理定律（如麦克斯韦方程组的微分形式表达）在不依赖具体度规的情况下依然保持形式简洁。它直接操作几何对象本身——线、面、体——而非它们的数值投影，完美契合了广义相对论中“背景独立性”的哲学。

## 本章小结

本章系统构建了处理反对称性质的代数框架——外代数（Exterior Algebra），其核心价值包括：

1. **楔积与反对称性**：引入了向量间的反对称乘法 $u \wedge v$，为描述有向面积和体积提供了最基本的代数算子。
2. **外幂空间 $\Lambda^k(V)$**：确立了由 $k$-向量构成的空间及其组合维数 $\binom{n}{k}$，并区分了可分解与不可分解形式。
3. **行列式的几何起源**：通过顶外幂 $\Lambda^n(V)$ 的标量映射，重新定义了行列式，从而使行列式的各种复合性质（如乘法法则）获得了范畴论层次的简洁证明。
4. **复合矩阵 $C_k(A)$**：展示了线性变换如何诱导外幂空间之间的映射，并揭示了子式与特征值乘积之间的深层联系。
5. **Grassmannian 几何**：利用 Plücker 嵌入将子空间集合转化为射影空间中的几何对象，通过 Plücker 关系刻画了子空间的代数约束。
6. **物理与几何算子**：介绍了 Hodge 对偶和缩并运算，为微分流形上的 Hodge 分解、Lie 导数及场论提供了代数地基。
