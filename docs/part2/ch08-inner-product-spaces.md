# 第 08 章 内积空间

<div class="context-flow" markdown>

**前置**：向量空间 (Ch04) · 正交性基础 (Ch07)

**本章脉络**：抽象内积定义 $\to$ 内积空间公理 $\to$ Cauchy-Schwarz 不等式及其几何证明 $\to$ 范数与距离的推广 $\to$ 正交补空间 $W^\perp$ $\to$ 伴随算子 $T^*$ 的代数本质 $\to$ 正规算子、自伴算子（Hermite）与酉算子 $\to$ 谱定理（Spectral Theorem）的深入剖析 $\to$ 正算子与矩阵平方根 $\to$ 完备性初步与有限维 Hilbert 空间

**延伸**：内积空间为线性空间引入了连续性与几何形状，是泛函分析的有限维缩影；谱定理是量子力学中“可观测量”理论的数学核心，揭示了算子最纯净的结构

</div>

内积空间是向量空间的一种强化。通过定义抽象内积，我们不仅能谈论线性组合，还能谈论长度、角度和逼近。它是现代信号处理、量子力学和最优化理论的数学地基。本章将超越 $\mathbb{R}^n$ 中的点积，建立起一套通用的度量几何框架。

---

## 08.1 抽象内积的定义

!!! definition "定义 08.1 (内积)"
    域 $F$（$\mathbb{R}$ 或 $\mathbb{C}$）上的向量空间 $V$ 配备内积 $\langle \cdot, \cdot \rangle$，若满足以下公理：
    1.  **正定性**：$\langle v, v \rangle \ge 0$，且为 0 当且仅当 $v=0$。
    2.  **共轭对称性**：$\langle u, v \rangle = \overline{\langle v, u \rangle}$。
    3.  **左线性**：$\langle au+bv, w \rangle = a\langle u, w \rangle + b\langle v, w \rangle$。

!!! note "导出范数"
    任何内积都自然导出一个范数 $\|v\| = \sqrt{\langle v, v \rangle}$。

---

## 08.2 核心不等式与正交补

!!! theorem "定理 08.1 (Cauchy-Schwarz 不等式)"
    对于内积空间中的任意向量 $u, v$：
    $$|\langle u, v \rangle| \le \|u\| \|v\|$$
    等号成立当且仅当 $u, v$ 线性相关。该不等式是定义向量间“夹角”的基础。

!!! definition "定义 08.2 (正交补 $W^\perp$)"
    设 $W$ 是 $V$ 的子空间。$W$ 的**正交补**定义为与 $W$ 中所有向量都正交的向量集合：
    $$W^\perp = \{ v \in V : \langle v, w \rangle = 0, \forall w \in W \}$$
    **性质**：$V = W \oplus W^\perp$（投影定理）。

---

## 08.3 伴随算子与特殊算子类

!!! definition "定义 08.3 (伴随算子 $T^*$)"
    设 $T$ 是 $V$ 上的线性算子。若存在算子 $T^*$ 使得对所有 $u, v \in V$ 均有：
    $$\langle Tu, v \rangle = \langle u, T^* v \rangle$$
    则称 $T^*$ 为 $T$ 的**伴随算子**。在标准正交基下，其矩阵表示为 $A^*$（共轭转置）。

!!! theorem "定理 08.2 (算子分类)"
    1.  **自伴算子 (Self-adjoint)**：$T^* = T$。其特征值全为实数。
    2.  **酉算子 (Unitary)**：$T^* T = TT^* = I$。保持长度和角度不变。
    3.  **正规算子 (Normal)**：$T^* T = TT^*$。

---

## 08.4 谱定理 (The Spectral Theorem)

!!! theorem "定理 08.3 (复谱定理)"
    算子 $T \in \mathcal{L}(V)$ 是正规算子，当且仅当存在 $V$ 的一组标准正交基由 $T$ 的特征向量构成。
    **物理意义**：这意味着正规算子可以通过坐标轴的旋转被完全解耦为独立的方向缩放。

---

## 练习题

**1. [公理判定] 验证 $\mathbb{R}^2$ 上的运算 $\langle u, v \rangle = u_1 v_1 + 2u_2 v_2$ 是否构成一个内积。**

??? success "参考答案"
    **验证步骤：**
    1. **正定性**：$\langle u, u \rangle = u_1^2 + 2u_2^2$。由于平方项非负且系数为正，结果必 $\ge 0$。且仅当 $u_1=u_2=0$ 时为 0。
    2. **对称性**：$\langle u, v \rangle = u_1 v_1 + 2u_2 v_2 = v_1 u_1 + 2v_2 u_2 = \langle v, u \rangle$。
    3. **线性性**：显然对于 $u$ 的各分量是线性的。
    **结论**：这是一个加权内积，满足所有公理。它定义的几何形状是一个以 $\sqrt{2}$ 为轴缩放的椭圆球面。

**2. [Cauchy-Schwarz] 证明：$|\langle u, v \rangle| \le \|u\| \|v\|$。**

??? success "参考答案"
    **证明过程：**
    1. 若 $v=0$，等式两边均为 0，结论成立。
    2. 若 $v \neq 0$，考虑向量 $u - \lambda v$ 的范数平方（必非负）：
       $\|u - \lambda v\|^2 = \langle u - \lambda v, u - \lambda v \rangle \ge 0$。
    3. 展开得：$\langle u, u \rangle - \bar{\lambda}\langle u, v \rangle - \lambda\langle v, u \rangle + |\lambda|^2 \langle v, v \rangle \ge 0$。
    4. 选取 $\lambda = \frac{\langle u, v \rangle}{\langle v, v \rangle}$（即 $u$ 在 $v$ 上的投影系数）。
    5. 代入计算：$\|u\|^2 - \frac{|\langle u, v \rangle|^2}{\|v\|^2} \ge 0$。
    6. 两边同乘 $\|v\|^2$ 并开方，即得 $|\langle u, v \rangle| \le \|u\| \|v\|$。

**3. [伴随算子] 设 $T$ 在标准基下的矩阵为 $A$。证明 $T^*$ 的矩阵为 $A^*$。**

??? success "参考答案"
    **证明：**
    1. 利用内积的矩阵表达：$\langle Tu, v \rangle = (Au)^* v = u^* A^* v$。
    2. 根据伴随定义：$\langle u, T^* v \rangle = u^* [T^*]_B v$。
    3. 对所有 $u, v$ 均成立，故对应的矩阵必须相等。
    **结论**：伴随算子的矩阵表示是原矩阵的共轭转置。

**4. [酉算子] 证明酉算子保持内积不变，即 $\langle Uu, Uv \rangle = \langle u, v \rangle$。**

??? success "参考答案"
    **推导：**
    1. $\langle Uu, Uv \rangle = \langle u, U^* U v \rangle$（利用伴随算子定义）。
    2. 根据酉算子定义：$U^* U = I$。
    3. $= \langle u, Iv \rangle = \langle u, v \rangle$。
    **几何意义**：酉变换（或实数域的正交变换）是刚体运动，不改变空间的长度和角度。

**5. [正规算子] 判定 $A = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$ 是否为正规矩阵。**

??? success "参考答案"
    **计算验证：**
    1. $A^* = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix}$。
    2. $AA^* = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = 2I$。
    3. $A^*A = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = 2I$。
    **结论**：由于 $AA^* = A^*A$，该矩阵是正规的。根据谱定理，它可以在复数域下酉对角化。

**6. [计算] 在多项式空间 $P_1$ 内积 $\langle f, g \rangle = \int_0^1 f(x)g(x) dx$ 下，计算 $1$ 与 $x$ 的内积及其夹角。**

??? success "参考答案"
    **计算步骤：**
    1. 内积：$\langle 1, x \rangle = \int_0^1 1 \cdot x dx = [x^2/2]_0^1 = 0.5$。
    2. 范数：$\|1\|^2 = \int_0^1 1 dx = 1$；$\|x\|^2 = \int_0^1 x^2 dx = 1/3$。
    3. 余弦值：$\cos\theta = \frac{0.5}{1 \cdot \sqrt{1/3}} = \frac{\sqrt{3}}{2}$。
    **结论**：夹角 $\theta = 30^\circ$ 或 $\pi/6$。

**7. [极化恒等式] 写出实内积空间中利用范数表达内积的公式。**

??? success "参考答案"
    **公式：**
    $$\langle u, v \rangle = \frac{1}{4}(\|u+v\|^2 - \|u-v\|^2)$$
    **证明：**
    展开右侧：$\frac{1}{4}(\langle u+v, u+v \rangle - \langle u-v, u-v \rangle) = \frac{1}{4}( (\|u\|^2 + 2\langle u, v \rangle + \|v\|^2) - (\|u\|^2 - 2\langle u, v \rangle + \|v\|^2) ) = \langle u, v \rangle$。
    这说明内积信息完全蕴含在范数（长度）信息中。

**8. [投影定理] 证明 $V = W \oplus W^\perp$。**

??? success "参考答案"
    **证明思路：**
    1. **存在性**：对任何 $v \in V$，设其在 $W$ 上的正交投影为 $p = \sum \langle v, e_i \rangle e_i$。显然 $p \in W$。
    2. 令 $z = v - p$。可以验证对所有 $e_i$，$\langle z, e_i \rangle = \langle v, e_i \rangle - \langle v, e_i \rangle = 0$。故 $z \in W^\perp$。
    3. 因此 $v = p + z$，即 $V = W + W^\perp$。
    4. **唯一性**：若 $x \in W \cap W^\perp$，则 $\langle x, x \rangle = 0 \implies x = 0$。
    **结论**：空间可以唯一分解为子空间及其正交补的直和。

**9. [Schur定理] 每个复方阵都酉相似于什么矩阵？其对角元是什么？**

??? success "参考答案"
    **结论：**
    每个复方阵都**酉相似于一个上三角矩阵**。
    其对角线上的元素正是该矩阵的**所有特征值**（计重数）。
    这是谱定理的基础，证明了通过旋转基，我们可以让任何算子表现为一种“分层”的结构。

**10. [自伴算子] 证明自伴算子的特征值必为实数。**

??? success "参考答案"
    **证明：**
    1. 设 $Tv = \lambda v$ 且 $v \neq 0$。
    2. 考虑 $\langle Tv, v \rangle = \langle \lambda v, v \rangle = \lambda \|v\|^2$。
    3. 同时也等于 $\langle v, T^* v \rangle$。由于 $T^*=T$，故等于 $\langle v, Tv \rangle = \langle v, \lambda v \rangle = \bar{\lambda} \|v\|^2$。
    4. 因此 $\lambda \|v\|^2 = \bar{\lambda} \|v\|^2$。由于 $\|v\|^2 > 0$，必有 $\lambda = \bar{\lambda}$。
    **结论**：自伴算子的特征值只能是实数。
