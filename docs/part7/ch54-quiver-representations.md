# 第 54 章 Quiver 表示论

<div class="context-flow" markdown>

**前置**：向量空间(Ch4) · 线性变换(Ch5) · 同时三角化(Ch63A) · 图论基础(Ch27)

**本章脉络**：Quiver（有向图）定义 → 路径代数 → 作为向量空间与映射的 Quiver 表示 → 表示间的态射 → 不可分解表示 → Gabriel 定理 → 与 Dynkin 图的联系 → 根系 → 在表示论中的应用

**延伸**：Quiver 表示论提供了一种几何化和图示化的方法来研究线性映射系统，将线性代数中的多样问题统一在单一的视觉框架下

</div>

一个 **Quiver** 简单来说就是一个有向图。Quiver 的一个**表示**为每个顶点分配一个向量空间，并为每条边分配一个线性映射。这种设置将对单一线性算子的研究（对应于只有一个顶点和一条自环的 Quiver）推广到了对相互依赖的空间和映射的整个配置的研究。Gabriel 定理惊人地将这些表示的复杂度与 Lie 理论中的经典 Dynkin 图联系了起来。

---

## 54.1 Quiver 及其表示

!!! definition "定义 54.1 (Quiver 表示)"
    域 $K$ 上 Quiver $Q = (V_0, E)$ 的一个表示 $V$ 是一组集合：
    - 向量空间集 $\{V_i\}_{i \in V_0}$；
    - 线性映射集 $\{\phi_\alpha: V_i \to V_j\}_{\alpha: i \to j \in E}$。

!!! theorem "定理 54.1 (Gabriel 定理, 1972)"
    一个 Quiver $Q$ 只有有限个不可分解表示（在同构意义下），当且仅当其底图是 $A_n, D_n, E_6, E_7, E_8$ 型的 **Dynkin 图**。

---

## 练习题

1. **[单环] 描述只有一个顶点和一条自环的 Quiver 的表示。**
   ??? success "参考答案"
       这对应于一个向量空间 $V$ 和一个自同态 $T: V \to V$。分类这些表示等价于 Jordan 标准形问题。不可分解表示恰好对应于各个 Jordan 块。

2. **[路径代数] 什么是 Quiver 的路径代数 $KQ$？**
   ??? success "参考答案"
       $KQ$ 是以 $Q$ 中所有路径（包括长度为 0 的路径）为基底的 $K$-向量空间。乘法定义为路径的级联。一个 Quiver 表示等价于其路径代数上的一个模。

3. **[维数向量] 定义表示的维数向量。**
   ??? success "参考答案"
       维数向量 $\mathbf{d} = (\dim V_1, \dots, \dim V_n)$ 记录了每个顶点处向量空间的大小。对于有限型 Quiver，不可分解表示的维数向量恰好是关联根系的正根。

4. **[态射] 两个表示 $V$ 和 $W$ 何时同构？**
   ??? success "参考答案"
       若存在一组可逆线性映射 $f_i: V_i \to W_i$，使得对于每条边 $\alpha: i \to j$，图表均对易：$f_j \circ \phi_\alpha = \psi_\alpha \circ f_i$，则称它们同构。

5. **[A2 Quiver] 找出 Quiver $1 \to 2$ ($A_2$ 型) 的所有不可分解表示。**
   ??? success "参考答案"
       共有三个：$(K \to 0)$，$(0 \to K)$，以及 $(K \xrightarrow{\text{id}} K)$。它们分别对应根 $(1,0), (0,1), (1,1)$。

6. **[Dynkin图] 画出 Dynkin 图 $A_3$ 并写出其对应的 Quiver。**
   ??? success "参考答案"
       $A_3$ 是 $1-2-3$。对应的 Quiver 可以是 $1 \to 2 \to 3$ 或 $1 \leftarrow 2 \to 3$。Gabriel 定理指出不可分解表示的数量与边的定向无关。

7. **[Kronecker] 分析 Kronecker Quiver（两个顶点，两条平行边 $1 \rightrightarrows 2$）。**
   ??? success "参考答案"
       这不是 Dynkin 图（属于 $\tilde{A}_1$ 型或欧几里得型）。它有无穷多个不可分解表示，对应于矩阵束 $(A, B)$ 理论和 Kronecker 标准形。

8. **[直和] 定义两个 Quiver 表示的直和。**
   ??? success "参考答案"
       在每个顶点处取 $(V \oplus W)_i = V_i \oplus W_i$，在每条边处，映射为分块对角矩阵 $\begin{pmatrix} \phi_\alpha & 0 \\ 0 & \psi_\alpha \end{pmatrix}$。这允许将复杂系统分解为更简单的不可分解单元。

9. **[子表示] 什么是 $V$ 的子表示？**
   ??? success "参考答案"
       子表示 $U \subseteq V$ 是一组子空间 $U_i \subseteq V_i$，使得对于每条边 $\alpha: i \to j$，映射 $\phi_\alpha$ 将 $U_i$ 映入 $U_j$。

10. **[反射函子] 简述 BGP 反射函子的作用。**
    ??? success "参考答案"
        反射函子提供了一种将一个 Quiver 的表示映射到在某个顶点处反转了所有边定向的另一个 Quiver 的表示的方法。这允许通过遍历根系来归纳证明 Gabriel 定理。

## 本章小结

本章通过基于图的表示统一了线性代数：

1. **图示综合**：将空间与映射的配置形式化为 Quiver 表示，将算子理论扩展到了网络结构。
2. **Gabriel 对称性**：揭示了有限表示类型与 Lie 理论中 Dynkin 图之间的深层联系。
3. **路径微积分**：将表示与路径代数挂钩，为基于图的映射提供了结合代数框架。
4. **根系对应**：确立了不可分解表示与根系几何之间的双射关系。
