# 第 08 章 内积空间

<div class="context-flow" markdown>

**前置**：向量空间(Ch4) · 正交性(Ch7)

**本章脉络**：抽象内积定义 → 范数与度量 → Cauchy-Schwarz 不等式 → 夹角与正交性 → 完备性与 Hilbert 空间初步 → 伴随算子 $T^*$ → 正规算子与酉算子 → 谱定理

**延伸**：内积空间为线性空间引入了连续性与几何形状，是泛函分析的有限维缩影

</div>

内积空间是向量空间的一种强化。通过定义内积，我们不仅能谈论线性组合，还能谈论长度、角度和逼近。它是现代信号处理、量子力学和最优化理论的数学地基。

---

## 08.1 定义与核心性质

!!! definition "定义 08.1 (内积)"
    域 $F$ 上的向量空间 $V$ 配备内积 $\langle \cdot, \cdot \rangle$，若满足：
    1. **正定性**：$\langle v, v \rangle \ge 0$，且为 0 当且仅当 $v=0$。
    2. **共轭对称性**：$\langle u, v \rangle = \overline{\langle v, u \rangle}$。
    3. **左线性**：$\langle au+bv, w \rangle = a\langle u, w \rangle + b\langle v, w \rangle$。

---

## 练习题

1. **[公理判定] 验证 $\mathbb{R}^2$ 上的运算 $\langle u, v \rangle = u_1 v_1 + 2u_2 v_2$ 是否构成一个内积。**
   ??? success "参考答案"
       是的。它是标准的点积的加权版本。由于系数 1 和 2 均为正，它满足正定性。线性性和对称性显然成立。这被称为加权内积。

2. **[Cauchy-Schwarz] 证明：$|\langle u, v \rangle| \le \|u\| \|v\|$。**
   ??? success "参考答案"
       考虑函数 $f(t) = \|u - tv\|^2 = \langle u-tv, u-tv \rangle \ge 0$。展开得 $t^2 \|v\|^2 - 2t \operatorname{Re}\langle u, v \rangle + \|u\|^2 \ge 0$。该二次式的判别式必须小于等于 0，整理即得 Cauchy-Schwarz 不等式。

3. **[伴随算子] 设 $T$ 在标准基下的矩阵为 $A$。证明其伴随算子 $T^*$ 的矩阵为 $A^*$。**
   ??? success "参考答案"
       由伴随算子定义：$\langle Tv, w \rangle = \langle v, T^* w \rangle$。
       左侧 $= (Av)^* w = v^* A^* w$。
       右侧 $= v^* (T^* w)$。
       由于对所有 $v, w$ 成立，故 $T^*$ 对应的矩阵必为 $A^*$（共轭转置）。

4. **[酉算子] 证明酉算子保持内积不变，即 $\langle Uu, Uv \rangle = \langle u, v \rangle$。**
   ??? success "参考答案"
       $\langle Uu, Uv \rangle = \langle u, U^* U v \rangle$。
       由于 $U$ 是酉算子，$U^* U = I$。
       故 $\langle u, Iv \rangle = \langle u, v \rangle$。这解释了为什么旋转和反射不改变角度和距离。

5. **[正规算子] 判定 $A = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$ 是否为正规矩阵。**
   ??? success "参考答案"
       $A^* = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix}$。
       $AA^* = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$。
       $A^*A = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$。
       $AA^* = A^*A$，故 $A$ 是正规矩阵。

6. **[计算] 在多项式内积 $\langle f, g \rangle = \int_0^1 f(x)g(x) dx$ 下，计算 $f(x)=1$ 与 $g(x)=x$ 的内积。**
   ??? success "参考答案"
       $\langle 1, x \rangle = \int_0^1 x dx = [\frac{1}{2}x^2]_0^1 = 1/2$。

7. **[极化恒等式] 写出实内积空间中利用范数表达内积的公式。**
   ??? success "参考答案"
       $\langle u, v \rangle = \frac{1}{4}(\|u+v\|^2 - \|u-v\|^2)$。这说明在实空间中，长度信息足以重建角度信息。

8. **[投影定理] 设 $W$ 是 $V$ 的子空间。证明 $V = W \oplus W^\perp$。**
   ??? success "参考答案"
       对任何 $v \in V$，取 $v$ 在 $W$ 上的正交投影 $p = \operatorname{proj}_W v$。令 $e = v - p$。显然 $p \in W$ 且 $\langle e, w \rangle = 0$ 对所有 $w \in W$ 成立，故 $e \in W^\perp$。这种分解是唯一的。

9. **[Schur定理应用] 每一个复方阵是否都酉相似于一个上三角矩阵？**
   ??? success "参考答案"
       是的。这是 Schur 分解定理。它保证了通过酉变换（数值稳定性好），我们可以揭示任何矩阵的特征值（位于对角线上）。

10. **[谱定理] 证明：若 $T$ 是自伴算子，则其特征值必为实数。**
    ??? success "参考答案"
        设 $Tv = \lambda v$，则 $\langle Tv, v \rangle = \langle \lambda v, v \rangle = \lambda \|v\|^2$。
        同时 $\langle Tv, v \rangle = \langle v, T^* v \rangle = \langle v, Tv \rangle = \langle v, \lambda v \rangle = \bar{\lambda} \|v\|^2$。
        由于 $v \neq 0$，必有 $\lambda = \bar{\lambda}$。

## 本章小结

内积空间确立了线性代数的度量衡：

1. **几何化**：引入了范数、角度和投影，使抽象空间具有了物理直观。
2. **算子对称性**：自伴、正规和酉算子的分类构成了谱理论的核心。
3. **逼近本质**：正交投影是解决一切线性逼近和误差最小化问题的数学终点。
