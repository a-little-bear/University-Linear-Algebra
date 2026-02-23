# 第 49 章 外代数与 Grassmann 几何

<div class="context-flow" markdown>

**前置**：向量空间(Ch4) · 行列式(Ch3) · 对偶空间(Ch13a) · 多线性代数(Ch21)

**本章脉络**：外积 $\wedge$ → 反对称性与交替张量 → 分级代数 → 外幂 $\Lambda^k(V)$ → 与行列式的关系 → Plücker 嵌入 → Grassmann 流形 $\operatorname{Gr}(k, n)$ → Hodge 星算子 → 内积/缩并

**延伸**：外代数是流形上的微积分（微分形式）和电磁学（麦克斯韦方程组）的数学语言

</div>

外代数，也称为 **Grassmann 代数**，将叉积和行列式的概念推广到了更高维度。通过引入**外积**（楔积 $\wedge$），它提供了一种与坐标无关的方法来描述有向面积、体积和线性子空间。这一代数结构是现代微分几何和研究 Grassmannian 流形的基础。

---

## 49.1 外积与交替张量

!!! definition "定义 49.1 (外积)"
    两个向量的外积 $u \wedge v$ 满足**反对称性**：
    $$u \wedge v = -v \wedge u$$
    由此推论，对任意向量 $v$，有 $v \wedge v = 0$。

!!! theorem "定理 49.1 (线性无关性与 $\wedge$)"
    向量集 $\{v_1, \dots, v_k\}$ 线性无关，当且仅当它们的外积非零：
    $$v_1 \wedge v_2 \wedge \dots \wedge v_k \neq 0$$

---

## 练习题

1. **[反对称性] 在 $\mathbb{R}^2$ 中计算 $(e_1 + e_2) \wedge (e_1 - e_2)$。**
   ??? success "参考答案"
       利用线性性和反对称性：$(e_1 + e_2) \wedge (e_1 - e_2) = e_1 \wedge e_1 - e_1 \wedge e_2 + e_2 \wedge e_1 - e_2 \wedge e_2 = 0 - e_1 \wedge e_2 - e_1 \wedge e_2 - 0 = -2(e_1 \wedge e_2)$。

2. **[维数] 计算 $n$ 维向量空间的 $k$ 阶外幂 $\Lambda^k(V)$ 的维数。**
   ??? success "参考答案"
       维数为组合数 $\binom{n}{k}$。这代表从 $n$ 个基向量中选取 $k$ 个来构造楔积基向量的方法数。

3. **[行列式] 证明对于线性算子 $T$，其在最高阶外幂 $\Lambda^n(V)$ 上的作用等价于乘以 $\det T$。**
   ??? success "参考答案"
       设 $\{e_1, \dots, e_n\}$ 为一组基。$T(e_1) \wedge \dots \wedge T(e_n) = (\det T) e_1 \wedge \dots \wedge e_n$。这一与坐标无关的定义是现代行列式理论的基础。

4. **[线性相关] 证明 $v \wedge w = 0$ 蕴含 $v$ 和 $w$ 线性相关。**
   ??? success "参考答案"
       若 $v, w$ 线性无关，则它们可以扩充为一组基。此时 $v \wedge w$ 是 $\Lambda^2(V)$ 的一个基元素，必非零。由逆否命题得证。

5. **[Plücker坐标] 定义由 $\{v_1, \dots, v_k\}$ 张成的 $k$ 维子空间的 Plücker 坐标。**
   ??? success "参考答案"
       Plücker 坐标是楔积 $v_1 \wedge \dots \wedge v_k$ 在 $\Lambda^k(V)$ 的基底下的系数。它们在相差一个标量倍数的意义下唯一确定了该子空间。

6. **[Grassmannian] 解释 Grassmannian $\operatorname{Gr}(k, n)$ 如何嵌入到射影空间 $\mathbb{P}(\Lambda^k(V))$ 中。**
   ??? success "参考答案"
       这是 Plücker 嵌入。它将每个 $k$ 维子空间映射到由其基向量楔积所确定的 $\Lambda^k(V)$ 中的直线。其像是满足 Plücker 关系的射影代数簇。

7. **[内积] 定义内积（收缩） $i_v \omega$ 及其几何意义。**
   ??? success "参考答案"
       $i_v$ 是一个 -1 阶的反导数，它将一个向量“代入”一个交替形式。几何上，它表示通过固定一个方向将体积形式降维为低维的面积形式。

8. **[Hodge星] 描述有向内积空间中的 Hodge 星算子 $*: \Lambda^k(V) \to \Lambda^{n-k}(V)$。**
   ??? success "参考答案"
       Hodge 星将一个 $k$-向量映射到其正交的 $(n-k)$-补。它满足 $\alpha \wedge *\beta = \langle \alpha, \beta \rangle \omega$，其中 $\omega$ 是体积形式。

9. **[复合映射] 若 $T: V \to V$ 的特征值为 $\lambda_1, \dots, \lambda_n$，诱导算子 $\Lambda^k(T)$ 的特征值是什么？**
   ??? success "参考答案"
       特征值是 $T$ 的 $n$ 个特征值中取 $k$ 个的所有可能乘积：$\{\lambda_{i_1} \dots \lambda_{i_k} : 1 \le i_1 < \dots < i_k \le n\}$。

10. **[可分解性] 定义“可分解” $k$-向量，并在 $\Lambda^2(\mathbb{R}^4)$ 中举出一个不可分解的例子。**
    ??? success "参考答案"
        若 $k$-向量能写成 $v_1 \wedge \dots \wedge v_k$ 的形式，则称其为可分解的。在 $\mathbb{R}^4$ 中，$\omega = e_1 \wedge e_2 + e_3 \wedge e_4$ 是不可分解的，因为它不满足 Plücker 关系 $\omega \wedge \omega = 0$（计算得 $\omega \wedge \omega = 2 e_1 \wedge e_2 \wedge e_3 \wedge e_4 \neq 0$）。

## 本章小结

本章通过外代数形式化了有向体积和子空间的几何：

1. **楔积运算**：定义了外积作为构造交替多线性形式的基础操作。
2. **结构对偶性**：探讨了子空间与可分解 $k$-向量之间的关系，引入了 Plücker 嵌入。
3. **流形地基**：确立了微分形式、Hodge 对偶以及广义 Stokes 定理的代数前提。
4. **谱映射**：展示了线性算子如何提升到外幂空间，为行列式和特征多项式提供了深层视角。
