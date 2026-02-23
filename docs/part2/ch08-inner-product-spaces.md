# 第 08 章 内积空间

<div class="context-flow" markdown>

**前置**：向量空间 (Ch04) · 正交性基础 (Ch07)

**本章脉络**：抽象内积定义 $\to$ 内积空间公理 $\to$ Cauchy-Schwarz 不等式 $\to$ 范数与距离的推广 $\to$ 伴随算子 $T^*$ $\to$ 正规算子、自伴算子（Hermite）与酉算子 $\to$ 谱定理（谱分解） $\to$ 正算子 $\to$ 完备性与 Hilbert 空间初步

**延伸**：内积空间为线性空间引入了连续性与几何形状，是泛函分析的有限维缩影；谱定理是量子力学中“可观测量”理论的数学核心

</div>

内积空间是向量空间的一种强化。通过定义抽象内积，我们不仅能谈论线性组合，还能谈论长度、角度和逼近。它是现代信号处理、量子力学和最优化理论的数学地基。

---

## 08.1 定义与核心性质

!!! definition "定义 08.1 (内积)"
    域 $F$ 上的向量空间 $V$ 配备内积 $\langle \cdot, \cdot \rangle$，若满足：
    1.  **正定性**：$\langle v, v \rangle \ge 0$，且为 0 当且仅当 $v=0$。
    2.  **共轭对称性**：$\langle u, v \rangle = \overline{\langle v, u \rangle}$。
    3.  **左线性**：$\langle au+bv, w \rangle = a\langle u, w \rangle + b\langle v, w \rangle$。

!!! theorem "定理 08.1 (Cauchy-Schwarz 不等式)"
    对于内积空间中的任意向量 $u, v$：
    $$|\langle u, v \rangle| \le \|u\| \|v\|$$
    等号成立当且仅当 $u, v$ 线性相关。

---

## 08.2 伴随算子

!!! definition "定义 08.2 (伴随算子 $T^*$)"
    设 $T$ 是 $V$ 上的线性算子。若存在算子 $T^*$ 使得对所有 $u, v \in V$ 均有：
    $$\langle Tu, v \rangle = \langle u, T^* v \rangle$$
    则称 $T^*$ 为 $T$ 的**伴随算子**。在标准基下，其矩阵表示为原矩阵的共轭转置 $A^*$。

---

## 08.3 算子分类与谱定理

!!! definition "定义 08.3 (特殊算子类)"
    1.  **自伴算子 (Self-adjoint)**：$T^* = T$。其特征值必为实数。
    2.  **酉算子 (Unitary)**：$T^* T = TT^* = I$。保持内积和范数不变。
    3.  **正规算子 (Normal)**：$T^* T = TT^*$。

!!! theorem "定理 08.2 (谱定理)"
    算子 $T$ 可以在标准正交基下对角化 $\iff$ $T$ 是正规算子。
    特别地，若 $T$ 是自伴算子，则其特征值全为实数，且存在由特征向量构成的标准正交基。

---

## 练习题

1. **[公理判定] 验证 $\mathbb{R}^2$ 上的运算 $\langle u, v \rangle = u_1 v_1 + 2u_2 v_2$ 是否构成一个内积。**
   ??? success "参考答案"
       是的。它是加权点积。由于权重均为正，它满足正定性；线性性和对称性显然成立。

2. **[Cauchy-Schwarz] 证明：$|\langle u, v \rangle| \le \|u\| \|v\|$。**
   ??? success "参考答案"
       考虑 $\|u - tv\|^2 \ge 0$，展开并取 $t = \langle u, v \rangle / \|v\|^2$ 即可。

3. **[伴随算子] 设 $T$ 在标准基下的矩阵为 $A$。证明 $T^*$ 的矩阵为 $A^*$。**
   ??? success "参考答案"
       $\langle Au, v \rangle = (Au)^* v = u^* A^* v = \langle u, A^* v \rangle$。由伴随算子定义即得。

4. **[酉算子] 证明酉算子保持内积不变。**
   ??? success "参考答案"
       $\langle Uu, Uv \rangle = \langle u, U^* U v \rangle = \langle u, Iv \rangle = \langle u, v \rangle$。

5. **[正规算子] 判定 $A = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$ 是否为正规矩阵。**
   ??? success "参考答案"
       $AA^* = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$，$A^*A = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$。相等，故是正规矩阵。

6. **[计算] 在多项式内积 $\langle f, g \rangle = \int_0^1 f(x)g(x) dx$ 下，计算 $1$ 与 $x$ 的内积。**
   ??? success "参考答案"
       $\int_0^1 x dx = 1/2$。

7. **[极化恒等式] 写出实内积空间中利用范数表达内积的公式。**
   ??? success "参考答案"
       $\langle u, v \rangle = \frac{1}{4}(\|u+v\|^2 - \|u-v\|^2)$。

8. **[投影定理] 证明 $V = W \oplus W^\perp$。**
   ??? success "参考答案"
       对任何 $v$，取 $v$ 在 $W$ 上的正交投影 $p$，则 $v = p + (v-p)$。容易验证 $v-p \in W^\perp$。

9. **[Schur定理] 每个复方阵都酉相似于什么矩阵？**
   ??? success "参考答案"
       上三角矩阵。

10. **[自伴算子] 证明自伴算子的特征值必为实数。**
    ??? success "参考答案"
        $\lambda \|v\|^2 = \langle Tv, v \rangle = \langle v, Tv \rangle = \bar{\lambda} \|v\|^2 \implies \lambda = \bar{\lambda}$。

## 本章小结

内积空间确立了线性代数的度量衡：

1.  **几何化**：引入了范数、角度和投影，使抽象空间具有了物理直观和逼近能力。
2.  **算子对称性**：自伴、正规和酉算子的分类构成了谱理论的核心，揭示了线性算子最稳定的结构。
3.  **逼近本质**：正交投影是解决一切线性逼近和误差最小化问题的数学终点。
