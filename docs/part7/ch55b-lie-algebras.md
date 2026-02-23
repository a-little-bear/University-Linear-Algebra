# 第 55B 章 李代数与无限小生成元

<div class="context-flow" markdown>

**前置**：Ch55 矩阵群 · Ch13 矩阵指数 $e^A$ · Ch21 张量

**脉络**：群在单位元的切空间 $	o$ 李括号 $[A, B] = AB - BA$ $	o$ 指数映射 $	o$ 结构常数 $	o$ 伴随表示

**延伸**：爱因斯坦的广义相对论依赖于时空的微分对称性，而量子力学的角动量理论（自旋）本质上就是 $SU(2)$ 李代数的表示。掌握了李代数，你就掌握了物理学中“守恒量”（如能量、动量）的数学起源。

</div>

如果说矩阵群（Lie Group）是某种平滑的几何形状（流形），那么**李代数（Lie Algebra）**就是这个形状在起始点（单位元）处的“切平面”。爱因斯坦的直觉告诉我们：与其研究复杂的旋转运动，不如研究旋转的“趋势”——即角速度。这种将全局对称性转化为线性空间操作的思路，是现代数学物理的巅峰。

---

## 55B.1 从指数映射到李代数

<div class="context-flow" markdown>

**核心思想**：若 $M(t)$ 是群中的一条曲线且 $M(0)=I$，则其导数 $X = M'(0)$ 捕捉了群的局部结构。

</div>

!!! definition "定义 55B.1 (矩阵群的李代数)"
    设 $G$ 是一个矩阵李群（如 $SO(n), SU(n)$）。其对应的**李代数** $\mathfrak{g}$ 是所有满足以下条件的矩阵 $X$ 的集合：
    对于所有实数 $t$，矩阵指数 $e^{tX}$ 始终属于群 $G$。
    
    $$ \mathfrak{g} = \{ X \in \mathbb{C}^{n 	imes n} : e^{tX} \in G, \forall t \in \mathbb{R} \} $$

!!! theorem "定理 55B.1 (经典群的李代数刻画)"
    1. **$SO(n)$ 的李代数 $\mathfrak{so}(n)$**：由所有**反对称矩阵**（$X^T = -X$）组成。
    2. **$SU(n)$ 的李代数 $\mathfrak{su}(n)$**：由所有**无迹厄米特矩阵**（$X^\dagger = -X, \operatorname{tr}(X)=0$）组成。
    3. **$SL(n, \mathbb{R})$ 的李代数 $\mathfrak{sl}(n)$**：由所有**无迹矩阵**（$\operatorname{tr}(X)=0$）组成。

??? proof "证明 (以 $SO(n)$ 为例)"
    若 $e^{tX} \in SO(n)$，则 $(e^{tX})^T (e^{tX}) = I$。
    对 $t$ 求导并在 $t=0$ 处求值：
    $\frac{d}{dt} (e^{tX^T} e^{tX})|_{t=0} = X^T e^0 e^0 + e^0 e^0 X = X^T + X = 0$。
    故 $X^T = -X$。反之亦然。$\blacksquare$

---

## 55B.2 李括号与代数结构

李代数不仅仅是一个向量空间，它还拥有一种特殊的“乘法”——**李括号**。

!!! definition "定义 55B.2 (李括号 / 交换子)"
    对于 $\mathfrak{g}$ 中的任意两个元素 $X, Y$，其**李括号**（Lie Bracket）定义为：
    
    $$ [X, Y] = XY - YX $$
    
    李代数 $\mathfrak{g}$ 在此运算下满足：
    1. **双线性性**：$[\alpha X + \beta Y, Z] = \alpha [X, Z] + \beta [Y, Z]$。
    2. **反对称性**：$[X, Y] = -[Y, X]$。
    3. **雅可比恒等式 (Jacobi Identity)**：$[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$。

!!! note "爱因斯坦的思考"
    为什么物理学家总是说“算符不交换”？因为 $[X, Y]$ 衡量了群结构非交换性的第一阶近似。在量子力学中，这直接对应于海森堡不确定性原理。

---

## 55B.3 结构常数与基

作为一个向量空间，李代数可以选取一组基 $\{ T_a \}$。此时李括号的结果必然可以由基向量线性表示：

$$ [T_a, T_b] = f_{ab}^{\ \ c} T_c $$

这里的 $f_{ab}^{\ \ c}$ 称为该李代数的**结构常数**。它们包含了群的全部遗传信息。

!!! example "例 55B.1 (三维旋转群 $\mathfrak{so}(3)$)"
    $\mathfrak{so}(3)$ 的一组标准基（生成元）是：
    $L_x = \begin{pmatrix} 0 & 0 & 0 \ 0 & 0 & -1 \ 0 & 1 & 0 \end{pmatrix}, L_y = \begin{pmatrix} 0 & 0 & 1 \ 0 & 0 & 0 \ -1 & 0 & 0 \end{pmatrix}, L_z = \begin{pmatrix} 0 & -1 & 0 \ 1 & 0 & 0 \ 0 & 0 & 0 \end{pmatrix}$
    它们满足：$[L_x, L_y] = L_z, [L_y, L_z] = L_x, [L_z, L_x] = L_y$。
    这就是我们在物理中学到的**角动量对易关系**。

---

## 10 道深度练习题 (含爱因斯坦直觉启发)

1. **[基础] 证明反对称矩阵的迹必为零。**
   ??? success "参考答案"
       设 $X^T = -X$。则 $\operatorname{tr}(X) = \operatorname{tr}(X^T) = \operatorname{tr}(-X) = -\operatorname{tr}(X)$。故 $2\operatorname{tr}(X) = 0 \Rightarrow \operatorname{tr}(X) = 0$。物理意义：旋转变换（由反对称阵生成）不改变空间体积。

2. **[理解] 给定 $X \in \mathfrak{g}$，证明对任何标量 $s, t$ 有 $[sX, tX] = 0$。**
   ??? success "参考答案"
       $[sX, tX] = (sX)(tX) - (tX)(sX) = st(X^2 - X^2) = 0$。这说明沿着同一个生成元前进的变换总是可交换的（它们属于群的同一个一参数子群）。

3. **[几何] 证明如果矩阵 $X$ 是反对称的，那么 $e^X$ 必然是正交矩阵。**
   ??? success "参考答案"
       $(e^X)^T = e^{X^T} = e^{-X}$。因此 $(e^X)^T e^X = e^{-X} e^X = e^0 = I$。这证明了 $\mathfrak{so}(n)$ 确实生成了旋转群。

4. **[代数] 验证雅可比恒等式对任意三个矩阵 $X, Y, Z$ 成立。**
   ??? success "参考答案"
       直接展开 $[X, [Y, Z]] = X(YZ-ZY) - (YZ-ZY)X = XYZ - XZY - YZX + ZYX$。
       将三个项全部展开，你会发现 12 个项两两抵消。

5. **[物理] 设 $\sigma_x, \sigma_y, \sigma_z$ 为泡利矩阵。证明 $i\sigma_x, i\sigma_y, i\sigma_z$ 构成了 $\mathfrak{su}(2)$ 的一组基。**
   ??? success "参考答案"
       泡利矩阵是厄米特的（$A^\dagger = A$），故 $X = i\sigma$ 满足 $X^\dagger = -i\sigma = -X$，且 $\operatorname{tr}(\sigma)=0$。计算对易关系可知它们满足 $\mathfrak{su}(2)$ 的代数结构。

6. **[转换] 已知旋转矩阵 $R(	heta) = \begin{pmatrix} \cos	heta & -\sin	heta \ \sin	heta & \cos	heta \end{pmatrix}$，求其生成元 $X$。**
   ??? success "参考答案"
       $X = \frac{dR}{d	heta}|_{	heta=0} = \begin{pmatrix} -\sin 0 & -\cos 0 \ \cos 0 & -\sin 0 \end{pmatrix} = \begin{pmatrix} 0 & -1 \ 1 & 0 \end{pmatrix}$。
       验证：$e^{	heta X} = I + 	heta X + \frac{	heta^2 X^2}{2!} + \dots$ 展开后即得原旋转矩阵。

7. **[深度] 证明若 $[X, Y] = 0$，则 $e^{X+Y} = e^X e^Y$。**
   ??? success "参考答案"
       利用幂级数展开。当 $X, Y$ 交换时，二项式定理 $(X+Y)^n = \sum \binom{n}{k} X^k Y^{n-k}$ 成立。
       这在物理上意味着：如果两个守恒量的生成元可交换，它们生成的对称性变换可以独立进行。

8. **[伴随] 定义 $\operatorname{ad}_X(Y) = [X, Y]$。证明 $\operatorname{ad}_X$ 是李代数上的一个线性变换。**
   ??? success "参考答案"
       线性性来自李括号对第二变元的线性性。$\operatorname{ad}_X$ 称为伴随表示，它描述了李代数如何“在自己内部旋转”。

9. **[结构常数] 证明结构常数 $f_{ab}^{\ \ c}$ 必然关于前两个指标反对称。**
   ??? success "参考答案"
       由 $[T_a, T_b] = -[T_b, T_a]$ 立即得出。

10. **[爱因斯坦挑战] 闵可夫斯基时空的洛伦兹群 $SO(1,3)$ 的李代数由哪些矩阵组成？**
    ??? success "参考答案"
        满足 $(\Lambda e^{tX})^T \eta (\Lambda e^{tX}) = \eta$。对 $t$ 求导得 $X^T \eta + \eta X = 0$。
        即 $\eta X$ 必须是一个反对称矩阵。这类矩阵包含 3 个空间旋转生成元和 3 个洛伦兹推升（Boost）生成元。
