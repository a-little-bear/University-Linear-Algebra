# 第 53B 章 闵可夫斯基空间与洛伦兹群

<div class="context-flow" markdown>

**前置**：Ch8 内积空间 · Ch9 二次型的惯性定理 · Ch55 矩阵群

**脉络**：正定内积放宽 $	o$ 不定内积空间 $	o$ 时空的"光锥"几何 $	o$ 保度量矩阵：洛伦兹群 $O(1,3)$ $	o$ 深层代数结构：$SL(2,\mathbb{C})$ 与旋量映射

**延伸**：从爱因斯坦的视角看，线性代数不仅仅是解方程的工具，它是**时空本身的几何语言**。正定内积定义了平坦空间的欧氏几何，而不定内积定义了狭义相对论的闵可夫斯基几何。本章的数学直接导出了光速不变原理、时间膨胀，并为狄拉克方程（反物质的发现）奠定了代数基础。

</div>

如果让爱因斯坦来审视线性代数，他会敏锐地指出：传统的内积空间（第 8 章）要求 $\langle \mathbf{v}, \mathbf{v} angle > 0$，这只描述了没有时间维度的绝对静止空间。真实的宇宙时空需要一种允许 $\langle \mathbf{v}, \mathbf{v} angle \leq 0$ 的"不定内积"。这种轻微的公理放宽，催生了极其绚丽的几何结构——光锥、洛伦兹变换以及隐藏在时空背后的旋量（Spinor）代数。

本章将从代数和几何的双重直觉出发，探索不定内积空间，并展示线性代数如何完美地编码了狭义相对论的物理定律。

---

## 53B.1 不定内积空间与度量符号

<div class="context-flow" markdown>

**突破点**：放弃正定性公理 $\langle \mathbf{v},\mathbf{v}angle > 0$ → 引入**度量张量** $\eta$ → 向量的"长度平方"可以为负或零！

</div>

!!! definition "定义 53B.1 (不定内积与伪欧氏空间)"
    设 $V$ 是 $n$ 维实向量空间。一个**非退化对称双线性型** $\langle \cdot, \cdot angle : V 	imes V 	o \mathbb{R}$ 称为 $V$ 上的**不定内积**。
    
    与标准内积（第 8 章）的唯一区别在于，我们**放弃了正定性公理**。存在非零向量 $\mathbf{v}$ 使得 $\langle \mathbf{v}, \mathbf{v} angle \leq 0$。配备了不定内积的空间称为**伪欧氏空间**（Pseudo-Euclidean space）。

根据第 9 章的 Sylvester 惯性定理，任何对称双线性型在适当的基（正交基）下，其矩阵可以对角化为仅包含 $+1, -1$ 和 $0$ 的对角阵。因为我们要求非退化，所以没有 $0$。

!!! definition "定义 53B.2 (度量符号 Signature)"
    如果在某组正交基下，内积的矩阵 $\eta$ 对角线上有 $p$ 个 $+1$ 和 $q$ 个 $-1$（$p+q=n$），则称该不定内积的**符号差（Signature）**为 $(p, q)$。
    
    $$ \eta = \operatorname{diag}(\underbrace{1, \ldots, 1}_{p}, \underbrace{-1, \ldots, -1}_{q}) $$

!!! note "爱因斯坦的洞察"
    物理学中，三维空间是 $(3,0)$ 的欧氏空间。当你加入一维时间，时空就变成了一个符号差为 $(3,1)$ 或 $(1,3)$ 的伪欧氏空间。这就是大名鼎鼎的**闵可夫斯基空间（Minkowski space）**。

---

## 53B.2 闵可夫斯基时空的几何：光锥

<div class="context-flow" markdown>

**因果性的几何化**：$\langle \mathbf{v},\mathbf{v}angle$ 的符号决定了向量是"空间"的、"时间"的，还是"光"的 → **非零向量可以与自己正交！**

</div>

设 $M^{1,3}$ 为四维时空，坐标记为 $\mathbf{x} = (t, x, y, z)^T$（采用自然单位制 $c=1$）。我们采用符号差 $(-, +, +, +)$ 的度量矩阵：

$$ \eta = \begin{pmatrix} -1 & 0 & 0 & 0 \ 0 & 1 & 0 & 0 \ 0 & 0 & 1 & 0 \ 0 & 0 & 0 & 1 \end{pmatrix} $$

两个时空向量 $\mathbf{v}_1 = (t_1, \mathbf{r}_1)^T$ 和 $\mathbf{v}_2 = (t_2, \mathbf{r}_2)^T$ 的内积为：
$$ \langle \mathbf{v}_1, \mathbf{v}_2 angle = \mathbf{v}_1^T \eta \mathbf{v}_2 = -t_1 t_2 + x_1 x_2 + y_1 y_2 + z_1 z_2 $$

向量 $\mathbf{v}$ 的"长度平方"（时空区间） $\Delta s^2 = \langle \mathbf{v}, \mathbf{v} angle = -t^2 + x^2 + y^2 + z^2$ 根据其符号，将整个时空划分为三个区域：

!!! definition "定义 53B.3 (时空向量的三分法)"
    1. **类时向量 (Timelike)**：$\langle \mathbf{v}, \mathbf{v} angle < 0$。代表有质量物体的运动轨迹。它们总是在光速以内。
    2. **类空向量 (Spacelike)**：$\langle \mathbf{v}, \mathbf{v} angle > 0$。代表同一时刻不同空间点的连线。因果关系无法跨越类空区间传递。
    3. **类光向量 (Lightlike / Null)**：$\langle \mathbf{v}, \mathbf{v} angle = 0$ 且 $\mathbf{v} 
eq \mathbf{0}$。代表光子的轨迹。所有类光向量构成的集合称为**光锥（Light Cone）**。

!!! theorem "定理 53B.1 (不可思议的光锥正交性)"
    在闵可夫斯基空间中，**类光向量与自身正交**（$\langle \mathbf{v}, \mathbf{v} angle = 0$）。
    此外，如果两个类光向量 $\mathbf{u}, \mathbf{v}$ 相互正交（$\langle \mathbf{u}, \mathbf{v} angle = 0$），那么它们必然在同一直线上（即线性相关）。

??? proof "证明"
    第一句话是定义的直接重述。这在欧氏空间中是不可想象的，但在不定内积中，一个向量在时间轴的投影平方抵消了在空间轴的投影平方。
    第二句话：设 $\mathbf{u} = (u_0, \vec{u})$, $\mathbf{v} = (v_0, \vec{v})$ 且均为类光向量，即 $u_0^2 = |\vec{u}|^2$，$v_0^2 = |\vec{v}|^2$。
    若 $\langle \mathbf{u}, \mathbf{v} angle = -u_0 v_0 + \vec{u} \cdot \vec{v} = 0$，则 $\vec{u} \cdot \vec{v} = u_0 v_0$。
    由欧氏空间的 Cauchy-Schwarz 不等式：$|\vec{u} \cdot \vec{v}| \leq |\vec{u}| |\vec{v}| = |u_0| |v_0|$。
    上式取到了等号，说明 $\vec{u}$ 和 $\vec{v}$ 必须平行，进而整个四维向量 $\mathbf{u}$ 与 $\mathbf{v}$ 线性相关。$\blacksquare$

---

## 53B.3 洛伦兹群 $O(1,3)$：时空的旋转

<div class="context-flow" markdown>

**寻找不变量**：正如正交矩阵 $Q^T I Q = I$ 保持空间长度不变，**洛伦兹矩阵** $\Lambda^T \eta \Lambda = \eta$ 保持时空因果结构不变

</div>

当坐标系变换时，爱因斯坦假定物理定律的形式不变，即"光速在任何惯性系下相同"。代数上，这意味着变换矩阵必须保持闵可夫斯基内积不变。

!!! definition "定义 53B.4 (洛伦兹群)"
    一个 $4 	imes 4$ 实矩阵 $\Lambda$ 称为**洛伦兹变换**（Lorentz transformation），如果它保持度量 $\eta$ 不变：
    
    $$ \Lambda^T \eta \Lambda = \eta $$
    
    所有这种矩阵构成的群称为**洛伦兹群**，记为 $O(1,3)$。

!!! proposition "命题 53B.1 (洛伦兹矩阵的代数性质)"
    设 $\Lambda \in O(1,3)$，则：
    1. $\det(\Lambda) = \pm 1$。
    2. 时间分量的放缩因子 $|\Lambda^0_{\ 0}| \geq 1$。

??? proof "证明"
    (1) 对 $\Lambda^T \eta \Lambda = \eta$ 两边取行列式：$\det(\Lambda)^2 \det(\eta) = \det(\eta)$。因为 $\det(\eta) = -1 
eq 0$，故 $\det(\Lambda)^2 = 1$。
    
    (2) 考察矩阵相乘的 $(0,0)$ 分量：
    $$ (\Lambda^T \eta \Lambda)^0_{\ 0} = -(\Lambda^0_{\ 0})^2 + (\Lambda^1_{\ 0})^2 + (\Lambda^2_{\ 0})^2 + (\Lambda^3_{\ 0})^2 = \eta_{00} = -1 $$
    因此 $(\Lambda^0_{\ 0})^2 = 1 + (\Lambda^1_{\ 0})^2 + (\Lambda^2_{\ 0})^2 + (\Lambda^3_{\ 0})^2 \geq 1$。$\blacksquare$

这说明洛伦兹群在拓扑上分为四个不连通的"岛屿"（由 $\det \Lambda = \pm 1$ 和 $\Lambda^0_{\ 0} \geq 1$ 还是 $\leq -1$ 划分）。满足 $\det \Lambda = 1$ 且 $\Lambda^0_{\ 0} \geq 1$ 的分支称为**固有正时洛伦兹群 $SO^+(1,3)$**，它代表了我们在物理上能实际实现的连续旋转和加速。

!!! example "例 53B.1 (双曲旋转 / 洛伦兹推升 Boost)"
    考虑仅在 $t-x$ 平面内发生的变换：
    $$ \Lambda = \begin{pmatrix} \cosh \phi & \sinh \phi & 0 & 0 \ \sinh \phi & \cosh \phi & 0 & 0 \ 0 & 0 & 1 & 0 \ 0 & 0 & 0 & 1 \end{pmatrix} $$
    容易验证 $\Lambda^T \eta \Lambda = \eta$（利用双曲恒等式 $\cosh^2\phi - \sinh^2\phi = 1$）。
    这在物理上代表一个沿着 $x$ 轴以速度 $v = 	anh \phi$ 运动的参考系变换。**洛伦兹推升，本质上就是时空平面内的一个伪旋转（双曲旋转）！**

---

## 53B.4 终极洞察：旋量与 $SL(2,\mathbb{C})$

<div class="context-flow" markdown>

**最深层的代数奇迹**：四维时空向量 $\cong$ $2	imes 2$ 厄米特矩阵 → 时空内积 $\cong$ 负的矩阵行列式 → $SO^+(1,3)$ 完美同构于 $SL(2,\mathbb{C})$ 的作用！

</div>

爱因斯坦的理论虽然极美，但直到保罗·狄拉克将一种名为**旋量（Spinor）**的纯代数对象引入时空，物理学才迎来了真正的完整。这里的关键，是一个将四维向量映射为 $2 	imes 2$ 矩阵的绝妙同构。

!!! theorem "定理 53B.2 (时空向量的厄米特矩阵表示)"
    将四维闵可夫斯基向量 $\mathbf{x} = (t, x, y, z)^T$ 映射为一个 $2 	imes 2$ 的复矩阵 $X$：
    
    $$ X = \begin{pmatrix} t + z & x - iy \ x + i y & t - z \end{pmatrix} $$
    
    则：
    1. $X$ 是一个**厄米特矩阵**（Hermitian matrix，$X^\dagger = X$）。任何 $2 	imes 2$ 厄米特矩阵都唯一对应一个四维实向量。
    2. $X$ 的**行列式等于时空区间的负值**：
       $$ \det(X) = (t+z)(t-z) - (x-iy)(x+iy) = t^2 - z^2 - x^2 - y^2 = -\langle \mathbf{x}, \mathbf{x} angle $$

*(注：这利用了泡利矩阵 $\sigma_1, \sigma_2, \sigma_3$ 和单位阵 $I$，即 $X = t I + x \sigma_1 + y \sigma_2 + z \sigma_3$。)*

!!! theorem "定理 53B.3 ($SL(2,\mathbb{C})$ 到洛伦兹群的同态)"
    设 $S \in SL(2,\mathbb{C})$（即行列式为 $1$ 的 $2 	imes 2$ 复可逆矩阵）。
    定义 $S$ 对上述厄米特矩阵 $X$ 的作用为：
    
    $$ X' = S X S^\dagger $$
    
    则：
    1. $X'$ 仍然是厄米特矩阵（即它对应一个新的四维实向量 $\mathbf{x}'$）。
    2. $\det(X') = \det(S)\det(X)\det(S^\dagger) = 1 \cdot \det(X) \cdot 1 = \det(X)$。
    
    **结论：** 这个映射 $X \mapsto X'$ 保持了行列式不变，从而**完全保持了闵可夫斯基内积（$\Delta s^2$）不变**！因此，每一个 $S \in SL(2,\mathbb{C})$ 都定义了一个洛伦兹变换 $\Lambda_S \in SO^+(1,3)$。

!!! note "物理学的基石：旋量"
    仔细观察上述变换：时空向量（$X$）被 $S$ 在左边乘一次，在右边被 $S^\dagger$ 乘一次。
    但是，如果存在某种"半个"时空向量，它只被 $S$ 从左边乘一次： $\psi \mapsto S \psi$，这种对象叫什么？
    **这就是旋量（Spinor）**。
    因为 $S$ 和 $-S$ 给出的洛伦兹变换 $\Lambda_S = \Lambda_{-S}$ 完全相同（二次方消去了符号），但它们对旋量的作用差一个负号。这就是为什么**自旋为 $1/2$ 的电子在空间中旋转 $360^\circ$ 后状态会反号，而旋转 $720^\circ$ 才恢复原状**！

## 练习题

1. **[基础] 在度量 $\eta = \operatorname{diag}(-1, 1, 1, 1)$ 下，计算向量 $\mathbf{v} = (5, 3, 4, 0)^T$ 的范数平方。它是哪类向量？**
   ??? success "参考答案"
       $\langle \mathbf{v}, \mathbf{v} \rangle = -(5^2) + 3^2 + 4^2 + 0^2 = -25 + 9 + 16 = 0$。
       它是一个**类光向量**。物理上，这可以代表一个沿特定方向飞行的光子。

2. **[概念] 证明两个类时向量之和未必是类时向量（举反例）。**
   ??? success "参考答案"
       设 $\mathbf{u} = (1, 0, 0, 0)^T$（类时，$-1 < 0$）和 $\mathbf{v} = (-1, 0, 0, 0)^T$（类时，$-1 < 0$）。
       其和 $\mathbf{u} + \mathbf{v} = \mathbf{0}$，零向量不是类时向量（其长度为 0）。
       更具启发性的例子：取 $\mathbf{u} = (1.1, 1, 0, 0)^T$ 和 $\mathbf{v} = (1.1, -1, 0, 0)^T$。

3. **[证明] 证明在闵可夫斯基空间中，不存在三个相互正交的类时向量。**
   ??? success "参考答案"
       类时向量指向光锥内部。由于时间维只有一维，任何两个正交的向量中，最多只能有一个是类时的（它占据了时间轴的“配额”）。

4. **[矩阵] 证明洛伦兹变换 $\Lambda$ 保持时空间距不变，即 $\langle \Lambda \mathbf{x}, \Lambda \mathbf{y} \rangle = \langle \mathbf{x}, \mathbf{y} \rangle$。**
   ??? success "参考答案"
       $\langle \Lambda \mathbf{x}, \Lambda \mathbf{y} \rangle = (\Lambda \mathbf{x})^T \eta (\Lambda \mathbf{y}) = \mathbf{x}^T (\Lambda^T \eta \Lambda) \mathbf{y}$。
       根据洛伦兹矩阵的定义 $\Lambda^T \eta \Lambda = \eta$，上式等于 $\mathbf{x}^T \eta \mathbf{y} = \langle \mathbf{x}, \mathbf{y} \rangle$。

5. **[双曲旋转] 计算例 53B.1 中洛伦兹推升矩阵的行列式。**
   ??? success "参考答案"
       $\det \Lambda = \cosh^2 \phi - \sinh^2 \phi = 1$。
       这符合固有洛伦兹变换的要求。

6. **[旋量] 已知时空点 $(t, x, y, z) = (1, 0, 0, 1)$，求其对应的 $2 \times 2$ 厄米特矩阵 $X$。**
   ??? success "参考答案"
       根据公式 $X = \begin{pmatrix} t+z & x-iy \\ x+iy & t-z \end{pmatrix} = \begin{pmatrix} 1+1 & 0 \\ 0 & 1-1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}$。
       验证：$\det(X) = 0$。对应时空区间 $-(1^2) + 0 + 0 + 1^2 = 0$（类光）。

7. **[群论] 证明若 $\Lambda_1, \Lambda_2$ 是洛伦兹矩阵，则其乘积 $\Lambda_1 \Lambda_2$ 也是洛伦兹矩阵。**
   ??? success "参考答案"
       $(\Lambda_1 \Lambda_2)^T \eta (\Lambda_1 \Lambda_2) = \Lambda_2^T (\Lambda_1^T \eta \Lambda_1) \Lambda_2 = \Lambda_2^T \eta \Lambda_2 = \eta$。
       这证明了洛伦兹变换确实构成一个群。

8. **[因果性] 证明若 $\mathbf{v}$ 是类时向量，$\Lambda$ 是正时洛伦兹变换（$\Lambda^0_{\ 0} \geq 1$），则 $\Lambda \mathbf{v}$ 仍然是类时向量。**
   ??? success "参考答案"
       洛伦兹变换保持内积符号，故 $\Lambda \mathbf{v}$ 长度平方仍为负。正时性质保证了变换不会将未来倒转为过去。

9. **[光速] 证明类光向量在任何洛伦兹变换下仍然是类光的。**
   ??? success "参考答案"
       $\langle \Lambda \mathbf{v}, \Lambda \mathbf{v} \rangle = \langle \mathbf{v}, \mathbf{v} \rangle = 0$。
       这在物理上意味着：光速在所有参考系中保持不变。

10. **[极限挑战] 设 $X = tI + \vec{x} \cdot \vec{\sigma}$。证明对任意 $S \in SL(2, \mathbb{C})$，变换 $X' = S X S^\dagger$ 总能保持 $\det(X)$ 不变。**
    ??? success "参考答案"
        $\det(X') = \det(S) \det(X) \det(S^\dagger) = 1 \cdot \det(X) \cdot \overline{\det(S)} = 1 \cdot \det(X) \cdot 1 = \det(X)$。
        这是 $SL(2,\mathbb{C})$ 作为洛伦兹群的双覆盖（Double Cover）的代数核心。

## 本章小结

从爱因斯坦和狄拉克的眼光重新审视线性代数，我们发现：

1. **放弃正定性**：将内积拓展为不定内积，得到了描述时空的伪欧氏空间（闵可夫斯基空间）。
2. **零范数不为零向量**：类光向量（光子）在这个度量下"长度"为零。
3. **等距群的拓展**：正交群 $O(n)$ 推广为洛伦兹群 $O(1,3)$，普通旋转推广为包含"双曲旋转"此时空推升。
4. **最深邃的同构**：四维时空的洛伦兹变换，在代数上竟然是复平面上 $2 	imes 2$ 行列式为 1 的矩阵群 $SL(2,\mathbb{C})$ 的伴随作用。这不仅是数学上的巧合，更是宇宙中一切物质粒子（费米子）存在的基础。
