# 第 47B 章 Fréchet 导数与灵敏度分析

<div class="context-flow" markdown>

**前置**：矩阵函数(Ch13) · 矩阵范数(Ch15) · 内积空间(Ch8) · 矩阵微积分基础(Ch47A) · 特征值(Ch6) · Sylvester 方程(Ch21)

**本章脉络**：Fréchet 导数定义与唯一性 $\to$ 基本性质（线性性、乘积规则、链式法则） $\to$ 经典矩阵函数的 Fréchet 导数（求逆、指数、对数、平方根） $\to$ Daleckii-Krein 定理与除商差分 $\to$ 特征值与奇异值的导数 $\to$ 矩阵分解的导数 $\to$ 条件数理论 $\to$ 高阶 Fréchet 导数 $\to$ 链式法则与矩阵方程的灵敏度分析

**延伸**：Fréchet 导数是矩阵流形上优化（Ch24）和自动微分（AD）的理论基础；Daleckii-Krein 定理将矩阵函数的导数与除商差分统一起来，是谱理论中最深刻的结果之一；矩阵分解的导数在机器学习中的可微编程和灵敏度分析中至关重要

</div>

在第 47A 章中，我们建立了标量函数对向量和矩阵的求导规则，并利用 Vec 算子与 Kronecker 积处理了矩阵对矩阵的导数。然而，对于矩阵函数——如矩阵指数 $e^A$、矩阵对数 $\log A$、矩阵平方根 $A^{1/2}$——我们需要一种更强大的导数理论。

Fréchet 导数将微积分推广到 Banach 空间之间的映射，在矩阵分析的语境中，它将矩阵函数 $f: \mathbb{C}^{n \times n} \to \mathbb{C}^{n \times n}$ 的导数定义为一个线性算子。这一框架不仅统一了前章的所有求导规则，更为矩阵函数的灵敏度分析、条件数理论和数值算法提供了坚实的理论基础。

---

## 47B.1 Fréchet 导数定义

<div class="context-flow" markdown>

**核心问题**：矩阵函数 $f(A)$（如 $e^A$、$\log A$、$A^{1/2}$）的导数如何严格定义？Fréchet 导数与 Gâteaux 导数有何区别？

</div>

!!! definition "定义 47B.1 (Fréchet 导数)"
    设 $\mathcal{U} \subset \mathbb{C}^{n \times n}$ 为开集，$f: \mathcal{U} \to \mathbb{C}^{n \times n}$ 为矩阵函数。$f$ 在 $A \in \mathcal{U}$ 处**Fréchet 可微**，若存在有界线性映射 $L_f(A): \mathbb{C}^{n \times n} \to \mathbb{C}^{n \times n}$，使得

    $$f(A + E) = f(A) + L_f(A)[E] + o(\|E\|),$$

    即 $\lim_{\|E\| \to 0} \dfrac{\|f(A + E) - f(A) - L_f(A)[E]\|}{\|E\|} = 0$。

    $L_f(A)$ 称为 $f$ 在 $A$ 处的 **Fréchet 导数**。它是一个从 $\mathbb{C}^{n \times n}$ 到 $\mathbb{C}^{n \times n}$ 的**线性算子**（不是一个矩阵，尽管可以用 $n^2 \times n^2$ 矩阵表示）。

!!! definition "定义 47B.2 (Gâteaux 导数)"
    $f$ 在 $A$ 处沿方向 $E$ 的 **Gâteaux 导数**（方向导数）定义为

    $$D_f(A)[E] = \lim_{t \to 0} \frac{f(A + tE) - f(A)}{t} = \frac{d}{dt} f(A + tE) \bigg|_{t=0},$$

    前提是该极限存在。

!!! theorem "定理 47B.1 (Fréchet 导数的唯一性)"
    若 $f$ 在 $A$ 处 Fréchet 可微，则 Fréchet 导数 $L_f(A)$ 是唯一的。

??? proof "证明"
    设 $L_1$ 和 $L_2$ 都是 $f$ 在 $A$ 处的 Fréchet 导数。则对任意 $E$：

    $$\|L_1[E] - L_2[E]\| = \|L_1[E] - (f(A+E) - f(A)) + (f(A+E) - f(A)) - L_2[E]\|$$
    $$\leq \|f(A+E) - f(A) - L_1[E]\| + \|f(A+E) - f(A) - L_2[E]\| = o(\|E\|).$$

    取 $E = t\hat{E}$（$\|\hat{E}\| = 1$），令 $t \to 0$：$|t| \cdot \|L_1[\hat{E}] - L_2[\hat{E}]\| = o(|t|)$，即 $\|L_1[\hat{E}] - L_2[\hat{E}]\| = 0$。

    由 $\hat{E}$ 的任意性，$L_1 = L_2$。 $\blacksquare$

!!! theorem "定理 47B.2 (Fréchet 导数与 Gâteaux 导数的关系)"
    1. 若 $f$ 在 $A$ 处 Fréchet 可微，则 $f$ 在 $A$ 处沿所有方向 Gâteaux 可微，且 $D_f(A)[E] = L_f(A)[E]$。
    2. 反之不一定成立：Gâteaux 可微不蕴含 Fréchet 可微。
    3. **充分条件**：若 $D_f(A)[E]$ 对所有 $E$ 存在，关于 $E$ 线性，且关于 $A$ **连续**（即 $A \mapsto D_f(A)$ 在算子范数下连续），则 $f$ 在 $A$ 处 Fréchet 可微。

??? proof "证明"
    **(1)** 取 $E = t\hat{E}$：

    $$\frac{f(A + t\hat{E}) - f(A)}{t} = L_f(A)[\hat{E}] + \frac{o(|t|)}{t}.$$

    令 $t \to 0$，右侧趋于 $L_f(A)[\hat{E}]$，即 Gâteaux 导数存在且等于 Fréchet 导数。

    **(2)** 反例：考虑 $\mathbb{R}^2 \to \mathbb{R}$ 的函数 $f(x, y) = \dfrac{x^2 y}{x^4 + y^2}$（$(x,y) \ne 0$），$f(0,0) = 0$。在原点沿任何方向的方向导数都存在且等于 $0$，但 $f$ 在原点不连续（沿抛物线 $y = x^2$ 趋近时极限为 $1/2$），因此不 Fréchet 可微。矩阵函数中也可构造类似的例子。

    **(3)** 这是有限维空间上的标准结果。关键是利用 Gâteaux 导数的连续性来控制非线性余项，证明它是 $o(\|E\|)$ 而非仅仅沿每条直线是 $o(|t|)$。 $\blacksquare$

!!! example "例 47B.1 (线性映射的 Fréchet 导数)"
    若 $f(A) = PAQ$（$P, Q$ 为常矩阵），则 $f(A+E) - f(A) = PEQ$，因此

    $$L_f(A)[E] = PEQ,$$

    且余项为零。线性映射的 Fréchet 导数就是映射本身。

---

## 47B.2 基本性质

<div class="context-flow" markdown>

**核心问题**：Fréchet 导数满足哪些运算规则？如何用这些规则计算复合函数的导数？

</div>

!!! theorem "定理 47B.3 (Fréchet 导数的基本性质)"
    设 $f, g: \mathcal{U} \to \mathbb{C}^{n \times n}$ 在 $A$ 处 Fréchet 可微。

    1. **线性性**：$L_{\alpha f + \beta g}(A)[E] = \alpha L_f(A)[E] + \beta L_g(A)[E]$，$\alpha, \beta \in \mathbb{C}$。
    2. **乘积规则**：$L_{fg}(A)[E] = L_f(A)[E] \cdot g(A) + f(A) \cdot L_g(A)[E]$。
    3. **链式法则**：若 $h = g \circ f$，$g$ 在 $f(A)$ 处 Fréchet 可微，则
       $$L_h(A)[E] = L_g(f(A))\big[L_f(A)[E]\big].$$
    4. **转置规则**：$L_{f^T}(A)[E] = (L_f(A)[E])^T$（对实矩阵）。
       对复矩阵的共轭转置：$L_{f^*}(A)[E] = (L_f(A)[E])^*$。

??? proof "证明"
    **(1)** 直接由线性映射之和的定义。

    **(2)** 令 $h(A) = f(A)g(A)$。

    $$h(A+E) = f(A+E)g(A+E) = (f(A) + L_f[E] + o(\|E\|))(g(A) + L_g[E] + o(\|E\|))$$
    $$= f(A)g(A) + L_f[E] \cdot g(A) + f(A) \cdot L_g[E] + O(\|E\|^2) + o(\|E\|)$$
    $$= h(A) + \big(L_f[E] \cdot g(A) + f(A) \cdot L_g[E]\big) + o(\|E\|).$$

    因此 $L_h(A)[E] = L_f(A)[E] \cdot g(A) + f(A) \cdot L_g(A)[E]$。

    **(3)** 链式法则：

    $$g(f(A+E)) = g(f(A) + L_f(A)[E] + o(\|E\|))$$
    $$= g(f(A)) + L_g(f(A))\big[L_f(A)[E] + o(\|E\|)\big] + o(\|L_f(A)[E] + o(\|E\|)\|)$$
    $$= g(f(A)) + L_g(f(A))[L_f(A)[E]] + o(\|E\|).$$

    最后一步利用了 $L_g$ 的有界性：$\|L_g(f(A))[o(\|E\|)]\| \leq \|L_g\| \cdot o(\|E\|) = o(\|E\|)$，以及 $\|L_f(A)[E]\| \leq \|L_f\| \cdot \|E\|$，所以 $o(\|L_f[E]\|) = o(\|E\|)$。

    **(4)** 由 $f(A+E)^T = f(A)^T + (L_f(A)[E])^T + o(\|E\|)$，直接读出 Fréchet 导数。 $\blacksquare$

!!! theorem "定理 47B.4 (矩阵幂的 Fréchet 导数)"
    对 $f(A) = A^k$（$k \geq 1$ 为正整数）：

    $$L_{A^k}(A)[E] = \sum_{j=0}^{k-1} A^j E A^{k-1-j}.$$

??? proof "证明"
    对 $k$ 进行归纳。$k = 1$：$L_{A}(A)[E] = E$，显然。

    设 $k-1$ 时成立。对 $A^k = A \cdot A^{k-1}$，由乘积规则：

    $$L_{A^k}(A)[E] = E \cdot A^{k-1} + A \cdot L_{A^{k-1}}(A)[E] = EA^{k-1} + A\sum_{j=0}^{k-2} A^j E A^{k-2-j}$$
    $$= EA^{k-1} + \sum_{j=0}^{k-2} A^{j+1} E A^{k-2-j} = EA^{k-1} + \sum_{j=1}^{k-1} A^j E A^{k-1-j} = \sum_{j=0}^{k-1} A^j E A^{k-1-j}.$$

    归纳完成。 $\blacksquare$

!!! example "例 47B.2 (矩阵多项式的 Fréchet 导数)"
    设 $p(A) = \sum_{k=0}^{m} c_k A^k$ 是矩阵多项式。由线性性和定理 47B.4：

    $$L_p(A)[E] = \sum_{k=1}^{m} c_k \sum_{j=0}^{k-1} A^j E A^{k-1-j}.$$

    当 $A$ 可对角化为 $A = V \operatorname{diag}(\lambda_i) V^{-1}$ 时：

    $$V^{-1} L_p(A)[E] V = \left[\frac{p(\lambda_i) - p(\lambda_j)}{\lambda_i - \lambda_j}\right] \odot (V^{-1}EV),$$

    其中 $\odot$ 是 Hadamard 积，$\dfrac{p(\lambda_i) - p(\lambda_j)}{\lambda_i - \lambda_j}$ 在 $\lambda_i = \lambda_j$ 时理解为 $p'(\lambda_i)$。这就是 Daleckii-Krein 定理的多项式情形。

---

## 47B.3 经典矩阵函数的 Fréchet 导数

<div class="context-flow" markdown>

**核心问题**：$A^{-1}$、$e^A$、$\log A$、$A^{1/2}$ 的 Fréchet 导数是什么？

</div>

!!! theorem "定理 47B.5 (矩阵求逆的 Fréchet 导数)"
    设 $f(A) = A^{-1}$（$A$ 可逆）。则

    $$L_{A^{-1}}(A)[E] = -A^{-1}EA^{-1}.$$

??? proof "证明"
    $(A + E)^{-1} = (A(I + A^{-1}E))^{-1} = (I + A^{-1}E)^{-1}A^{-1}$。

    当 $\|A^{-1}E\| < 1$ 时，$(I + A^{-1}E)^{-1} = I - A^{-1}E + O(\|E\|^2)$。

    因此 $(A+E)^{-1} = A^{-1} - A^{-1}EA^{-1} + O(\|E\|^2)$。

    余项估计：$\|(A+E)^{-1} - A^{-1} + A^{-1}EA^{-1}\| = \|A^{-1}E A^{-1}E(I+A^{-1}E)^{-1}A^{-1}\| = O(\|E\|^2)$。

    因此 $L_{A^{-1}}(A)[E] = -A^{-1}EA^{-1}$。 $\blacksquare$

!!! theorem "定理 47B.6 (矩阵指数的 Fréchet 导数)"
    设 $f(A) = e^A$。则

    $$L_{e^A}(A)[E] = \int_0^1 e^{sA} E \, e^{(1-s)A} \, ds.$$

??? proof "证明"
    **方法一（Duhamel 公式）**：设 $F(t) = e^{(A+tE)}$，我们要求 $F'(0)$。

    定义 $G(t) = e^{-A} F(t) = e^{-A} e^{A+tE}$。则 $G(0) = I$。

    $G(t)$ 满足微分方程：

    $$G'(t) = e^{-A} \frac{d}{dt} e^{(A+tE)} = e^{-A} \cdot E \cdot e^{(A+tE)} \cdot (\text{需要更仔细})$$

    这里不能直接用 $\dfrac{d}{dt}e^{Bt} = Be^{Bt}$ 的公式，因为 $A$ 和 $E$ 一般不对易。

    改用另一种方法。考虑辅助函数 $\Phi(t, s) = e^{s(A+tE)}$，它满足

    $$\frac{\partial \Phi}{\partial s} = (A+tE)\Phi, \quad \Phi(t, 0) = I.$$

    对 $t$ 求导：设 $\Psi(t, s) = \dfrac{\partial \Phi}{\partial t}$，则

    $$\frac{\partial \Psi}{\partial s} = (A+tE)\Psi + E\Phi, \quad \Psi(t, 0) = 0.$$

    在 $t = 0$ 时，$\Phi(0, s) = e^{sA}$，$\Psi$ 满足

    $$\frac{\partial \Psi}{\partial s} = A\Psi + Ee^{sA}, \quad \Psi(0, 0) = 0.$$

    这是关于 $\Psi$ 的线性常微分方程，利用常数变易法（variation of constants）求解：

    $$\Psi(0, s) = e^{sA} \int_0^s e^{-\tau A} E e^{\tau A} d\tau.$$

    取 $s = 1$：

    $$F'(0) = \Psi(0, 1) = e^A \int_0^1 e^{-\tau A} E e^{\tau A} d\tau = \int_0^1 e^{(1-\tau)A} E e^{\tau A} d\tau.$$

    换元 $u = 1 - \tau$... 不，直接整理：令 $s = 1 - \tau$，则

    $$F'(0) = \int_0^1 e^{sA} E e^{(1-s)A} ds.$$

    **余项估计**：$e^{A+tE} = e^A + t\int_0^1 e^{sA}Ee^{(1-s)A}ds + O(t^2)$，其中 $O(t^2)$ 项可通过迭代 Duhamel 公式获得精确界。具体地：

    $$\left\|e^{A+tE} - e^A - t\int_0^1 e^{sA}Ee^{(1-s)A}ds\right\| \leq \frac{t^2 \|E\|^2}{2} e^{\|A\| + |t|\|E\|}.$$

    这表明余项确实是 $O(t^2)$，即 $o(t)$，从而 Fréchet 导数（而非仅 Gâteaux 导数）存在。 $\blacksquare$

!!! theorem "定理 47B.7 (矩阵对数的 Fréchet 导数)"
    设 $A$ 无非正实特征值，$f(A) = \log A$（主对数）。则

    $$L_{\log}(A)[E] = \int_0^1 (sA + (1-s)I)^{-1} E \, (sA + (1-s)I)^{-1} \cdot s \, ds \cdot A... $$

    更标准的形式：利用积分表示 $\log A = \int_0^{\infty} \left[(tI + I)^{-1} - (tI + A)^{-1}\right] dt$，

    $$L_{\log}(A)[E] = \int_0^{\infty} (tI + A)^{-1} E (tI + A)^{-1} dt.$$

    或者等价地，利用参数化 $\log A = \int_0^1 (I - s(I - A))^{-1}(I - A)^{-1}... $

    最简洁的形式：当 $A$ 可对角化为 $A = V\operatorname{diag}(\lambda_i)V^{-1}$ 时，

    $$V^{-1}L_{\log}(A)[E]V = \left[\frac{\log \lambda_i - \log \lambda_j}{\lambda_i - \lambda_j}\right] \odot (V^{-1}EV),$$

    其中 $\dfrac{\log \lambda_i - \log \lambda_j}{\lambda_i - \lambda_j}$ 在 $\lambda_i = \lambda_j$ 时理解为 $1/\lambda_i$。

??? proof "证明"
    由 $\log(e^A) = A$，两侧取 Fréchet 导数，利用链式法则：

    $$L_{\log}(e^A)\left[L_{e^A}(A)[E]\right] = E.$$

    因此 $L_{\log}(e^A) = (L_{e^A}(A))^{-1}$，即对数的 Fréchet 导数是指数的 Fréchet 导数的逆算子。

    积分表示：由 $\log A = \int_0^{\infty}[(1+t)^{-1}I - (A+tI)^{-1}]dt$（Cauchy 积分），对 $A$ 取 Fréchet 导数：

    $$L_{\log}(A)[E] = \int_0^{\infty} L_{(A+tI)^{-1}}(A)[E] dt = \int_0^{\infty} (A+tI)^{-1} E (A+tI)^{-1} dt.$$

    这里利用了 $L_{(A+tI)^{-1}}(A)[E] = -(A+tI)^{-1}E(A+tI)^{-1}$（由定理 47B.5），以及被积函数前面的负号与 $-$ 号抵消。积分收敛性由 $A$ 的特征值有正实部保证。 $\blacksquare$

!!! theorem "定理 47B.8 (矩阵平方根的 Fréchet 导数)"
    设 $A$ 无非正实特征值，$S = A^{1/2}$（主平方根）。则 $L_{A^{1/2}}(A)[E] = X$，其中 $X$ 是 **Sylvester 方程**

    $$SX + XS = E$$

    的唯一解。

??? proof "证明"
    由 $f(A)^2 = A$，两侧取 Fréchet 导数。设 $S = A^{1/2}$，$X = L_{A^{1/2}}(A)[E]$。

    对 $g(S) = S^2$ 利用乘积规则：$L_{S^2}(S)[X] = XS + SX$。

    由链式法则：$L_{S^2}(S)[L_{A^{1/2}}(A)[E]] = L_A(A)[E] = E$。

    即 $XS + SX = E$。

    由于 $S = A^{1/2}$ 的特征值 $\sqrt{\lambda_i}$ 都有正实部（因为 $A$ 无非正实特征值），$\sigma(S) \cap \sigma(-S) = \emptyset$，因此 Sylvester 方程有唯一解。 $\blacksquare$

!!! example "例 47B.3 (矩阵指数的 Fréchet 导数计算)"
    设 $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$，$E = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。

    $$L_{e^A}(A)[E] = \int_0^1 e^{sA} E \, e^{(1-s)A} ds = \int_0^1 \begin{pmatrix} e^s & 0 \\ 0 & e^{2s} \end{pmatrix} \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} \begin{pmatrix} e^{1-s} & 0 \\ 0 & e^{2(1-s)} \end{pmatrix} ds$$

    $$= \int_0^1 \begin{pmatrix} 0 & e^s \\ 0 & 0 \end{pmatrix} \begin{pmatrix} e^{1-s} & 0 \\ 0 & e^{2-2s} \end{pmatrix} ds = \int_0^1 \begin{pmatrix} 0 & e^{2-s} \\ 0 & 0 \end{pmatrix} ds$$

    $$= \begin{pmatrix} 0 & [-e^{2-s}]_0^1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & e^2 - e \\ 0 & 0 \end{pmatrix}.$$

    这符合 Daleckii-Krein 公式：$\dfrac{e^{\lambda_1} - e^{\lambda_2}}{\lambda_1 - \lambda_2} = \dfrac{e^1 - e^2}{1 - 2} = e^2 - e$。

---

## 47B.4 Daleckii-Krein 定理

<div class="context-flow" markdown>

**核心问题**：是否有一个统一的公式，用矩阵的特征值和特征向量来表达任意矩阵函数的 Fréchet 导数？

</div>

Daleckii-Krein 定理是矩阵函数 Fréchet 导数理论的核心结果。它将 Fréchet 导数表示为**除商差分**（divided differences）与特征投影的组合，提供了统一的理论框架。

!!! definition "定义 47B.3 (除商差分)"
    设 $f$ 是定义在包含 $\lambda_i, \lambda_j$ 的区域上的函数。$f$ 的**一阶除商差分**（first divided difference）定义为

    $$f[\lambda_i, \lambda_j] = \begin{cases} \dfrac{f(\lambda_i) - f(\lambda_j)}{\lambda_i - \lambda_j}, & \lambda_i \ne \lambda_j, \\ f'(\lambda_i), & \lambda_i = \lambda_j. \end{cases}$$

    $f[\lambda_i, \lambda_j]$ 关于两个变量是对称的，且当 $f$ 解析时关于 $(\lambda_i, \lambda_j)$ 连续。

    **高阶除商差分**递归定义：

    $$f[\lambda_0, \lambda_1, \ldots, \lambda_k] = \frac{f[\lambda_1, \ldots, \lambda_k] - f[\lambda_0, \ldots, \lambda_{k-1}]}{\lambda_k - \lambda_0}.$$

!!! theorem "定理 47B.9 (Daleckii-Krein 定理——可对角化情形)"
    设 $A \in \mathbb{C}^{n \times n}$ 可对角化：$A = V \operatorname{diag}(\lambda_1, \ldots, \lambda_n) V^{-1}$，其中 $V = (v_1, \ldots, v_n)$ 的列是右特征向量。设 $w_i^* = e_i^T V^{-1}$ 是对应的左特征向量（$w_i^* v_j = \delta_{ij}$）。

    设 $f$ 在 $A$ 的特征值上可微。则

    $$L_f(A)[E] = \sum_{i=1}^{n}\sum_{j=1}^{n} f[\lambda_i, \lambda_j] \, (w_i^* E v_j) \, v_i w_j^*,$$

    或等价地，

    $$V^{-1} L_f(A)[E] V = f^{[1]}(\Lambda) \odot (V^{-1}EV),$$

    其中 $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，$f^{[1]}(\Lambda)$ 是**除商差分矩阵**（Loewner matrix）：

    $$(f^{[1]}(\Lambda))_{ij} = f[\lambda_i, \lambda_j] = \begin{cases} \dfrac{f(\lambda_i) - f(\lambda_j)}{\lambda_i - \lambda_j}, & \lambda_i \ne \lambda_j, \\ f'(\lambda_i), & \lambda_i = \lambda_j, \end{cases}$$

    $\odot$ 表示 Hadamard 积。

??? proof "证明"
    **第一步**：多项式情形。

    对 $f(A) = A^k$，由定理 47B.4：

    $$L_{A^k}(A)[E] = \sum_{j=0}^{k-1} A^j E A^{k-1-j}.$$

    设 $A = V\Lambda V^{-1}$，$\tilde{E} = V^{-1}EV$。则

    $$V^{-1} L_{A^k}(A)[E] V = \sum_{j=0}^{k-1} \Lambda^j \tilde{E} \Lambda^{k-1-j}.$$

    计算第 $(p,q)$ 元素：

    $$\left(\sum_{j=0}^{k-1} \Lambda^j \tilde{E} \Lambda^{k-1-j}\right)_{pq} = \tilde{E}_{pq} \sum_{j=0}^{k-1} \lambda_p^j \lambda_q^{k-1-j}.$$

    当 $\lambda_p \ne \lambda_q$ 时，$\sum_{j=0}^{k-1} \lambda_p^j \lambda_q^{k-1-j} = \dfrac{\lambda_p^k - \lambda_q^k}{\lambda_p - \lambda_q}$（几何级数求和）。

    当 $\lambda_p = \lambda_q$ 时，$\sum_{j=0}^{k-1} \lambda_p^{k-1} = k\lambda_p^{k-1} = f'(\lambda_p)$（$f(z) = z^k$）。

    因此多项式情形成立：$V^{-1}L_{A^k}(A)[E]V = f^{[1]}(\Lambda) \odot \tilde{E}$。

    **第二步**：推广到一般解析函数。

    设 $f$ 在包含 $\sigma(A)$ 的开集上解析。则 $f$ 可以在该区域上一致逼近为多项式序列（Runge 定理）。由 Fréchet 导数的连续性（关于 $f$ 的一致逼近），极限交换给出一般情形的结论。

    更直接地：$f(A)$ 可用 Cauchy 积分表示为 $f(A) = \frac{1}{2\pi i}\oint_\Gamma f(z)(zI - A)^{-1}dz$，取 Fréchet 导数：

    $$L_f(A)[E] = \frac{1}{2\pi i}\oint_\Gamma f(z)(zI - A)^{-1}E(zI-A)^{-1}dz.$$

    代入 $A = V\Lambda V^{-1}$：

    $$(zI - A)^{-1} = V(zI - \Lambda)^{-1}V^{-1},$$

    $$V^{-1}L_f(A)[E]V = \frac{1}{2\pi i}\oint_\Gamma f(z)(zI-\Lambda)^{-1}\tilde{E}(zI-\Lambda)^{-1}dz.$$

    第 $(p,q)$ 元素：

    $$\frac{\tilde{E}_{pq}}{2\pi i}\oint_\Gamma \frac{f(z)}{(z-\lambda_p)(z-\lambda_q)}dz.$$

    当 $\lambda_p \ne \lambda_q$ 时，用部分分式分解和留数定理：

    $$\frac{1}{2\pi i}\oint_\Gamma \frac{f(z)}{(z-\lambda_p)(z-\lambda_q)}dz = \frac{f(\lambda_p) - f(\lambda_q)}{\lambda_p - \lambda_q} = f[\lambda_p, \lambda_q].$$

    当 $\lambda_p = \lambda_q$ 时：

    $$\frac{1}{2\pi i}\oint_\Gamma \frac{f(z)}{(z-\lambda_p)^2}dz = f'(\lambda_p) = f[\lambda_p, \lambda_p].$$

    因此 $V^{-1}L_f(A)[E]V = f^{[1]}(\Lambda) \odot \tilde{E}$。 $\blacksquare$

!!! theorem "定理 47B.10 (Daleckii-Krein 定理——一般情形)"
    设 $A \in \mathbb{C}^{n \times n}$（不一定可对角化），Jordan 标准形为 $A = VJV^{-1}$，其中 $J = \operatorname{diag}(J_1, \ldots, J_s)$ 是 Jordan 块。设 $f$ 在 $\sigma(A)$ 的邻域上足够光滑。则

    $$V^{-1}L_f(A)[E]V = f^{[1]}(J) \star (V^{-1}EV),$$

    其中 $\star$ 是一种推广的 Hadamard 型运算，$f^{[1]}(J)$ 的分块 $(p,q)$ 定义为：

    若 $J_p$ 是 $m_p \times m_p$ Jordan 块（特征值 $\lambda_p$），$J_q$ 是 $m_q \times m_q$ Jordan 块（特征值 $\lambda_q$），则 $f^{[1]}(J)$ 在第 $(p,q)$ 分块位置是 $m_p \times m_q$ 矩阵，其 $(r,s)$ 元素为高阶除商差分

    $$f[\underbrace{\lambda_p, \ldots, \lambda_p}_{r \text{ 个}}, \underbrace{\lambda_q, \ldots, \lambda_q}_{s \text{ 个}}].$$

    此定理的完整证明需要 Jordan 块上 Cauchy 积分的留数计算。

!!! example "例 47B.4 (不同函数的除商差分矩阵)"
    设 $A$ 有特征值 $\lambda_1, \lambda_2$（$\lambda_1 \ne \lambda_2$）。

    (a) $f(z) = e^z$：$f[\lambda_1, \lambda_2] = \dfrac{e^{\lambda_1} - e^{\lambda_2}}{\lambda_1 - \lambda_2}$，$f[\lambda_i, \lambda_i] = e^{\lambda_i}$。

    (b) $f(z) = \log z$：$f[\lambda_1, \lambda_2] = \dfrac{\log\lambda_1 - \log\lambda_2}{\lambda_1 - \lambda_2}$，$f[\lambda_i, \lambda_i] = 1/\lambda_i$。

    (c) $f(z) = z^{1/2}$：$f[\lambda_1, \lambda_2] = \dfrac{\sqrt{\lambda_1} - \sqrt{\lambda_2}}{\lambda_1 - \lambda_2} = \dfrac{1}{\sqrt{\lambda_1} + \sqrt{\lambda_2}}$，$f[\lambda_i, \lambda_i] = \dfrac{1}{2\sqrt{\lambda_i}}$。

    (d) $f(z) = z^{-1}$：$f[\lambda_1, \lambda_2] = \dfrac{1/\lambda_1 - 1/\lambda_2}{\lambda_1 - \lambda_2} = -\dfrac{1}{\lambda_1\lambda_2}$，$f[\lambda_i, \lambda_i] = -1/\lambda_i^2$。

---

## 47B.5 特征值和奇异值的导数

<div class="context-flow" markdown>

**核心问题**：矩阵的特征值 $\lambda_i(A)$ 和奇异值 $\sigma_i(A)$ 如何对 $A$ 的扰动做出响应？

</div>

特征值和奇异值作为 $A$ 的函数，其微分理论是矩阵扰动分析的核心。然而，特征值一般不是矩阵元素的可微函数（例如重特征值分裂时），需要格外小心。

!!! theorem "定理 47B.11 (单特征值的导数)"
    设 $A(t)$ 是依赖于实参数 $t$ 的矩阵族，$A(0) = A$。设 $\lambda_0$ 是 $A$ 的**单特征值**（代数重数为 $1$），对应右特征向量 $v$（$Av = \lambda_0 v$）和左特征向量 $w$（$w^*A = \lambda_0 w^*$），归一化使得 $w^*v = 1$。

    则存在 $t$ 的邻域内的解析函数 $\lambda(t)$ 满足 $\lambda(0) = \lambda_0$，且

    $$\lambda'(0) = w^* A'(0) v.$$

    更一般地，若 $A$ 受矩阵 $E$ 的扰动（$A \to A + \varepsilon E$），则

    $$\frac{d\lambda}{d\varepsilon}\bigg|_{\varepsilon=0} = w^* E v.$$

??? proof "证明"
    由 $(A + \varepsilon E)v(\varepsilon) = \lambda(\varepsilon)v(\varepsilon)$，对 $\varepsilon$ 求导并在 $\varepsilon = 0$ 处取值：

    $$Ev + Av'(0) = \lambda'(0)v + \lambda_0 v'(0).$$

    两侧左乘 $w^*$：

    $$w^*Ev + w^*Av'(0) = \lambda'(0)w^*v + \lambda_0 w^*v'(0).$$

    由 $w^*A = \lambda_0 w^*$，得 $w^*Av'(0) = \lambda_0 w^*v'(0)$，两项抵消：

    $$w^*Ev = \lambda'(0) \cdot 1.$$

    因此 $\lambda'(0) = w^*Ev$。 $\blacksquare$

!!! theorem "定理 47B.12 (Hermite 矩阵特征值的导数)"
    设 $A(t)$ 是 Hermite 矩阵族（$A(t) = A(t)^*$），特征值 $\lambda_1(t) \geq \lambda_2(t) \geq \cdots \geq \lambda_n(t)$。

    1. 每个 $\lambda_i(t)$ 是 $t$ 的**连续**函数。
    2. 若 $\lambda_i(0)$ 是单特征值，对应单位特征向量 $v_i$，则
       $$\lambda_i'(0) = v_i^* A'(0) v_i.$$
    3. 即使 $\lambda_i(0)$ 是重特征值（重数 $r$），$\lambda_i(t)$ 仍然是 Lipschitz 连续的，且
       $$\sum_{k=i}^{i+r-1} \lambda_k'(0) = \operatorname{tr}(P_i A'(0)),$$
       其中 $P_i$ 是对应于 $\lambda_i(0)$ 的特征投影。

??? proof "证明"
    **(2)** 对 Hermite 矩阵，左特征向量 $w = v$（$w^*A = \lambda w^* \Leftrightarrow Aw = \lambda w$）。由定理 47B.11，$\lambda'(0) = v^*A'(0)v$。

    **(3)** 设 $\lambda_0$ 是 $A(0)$ 的 $r$ 重特征值，对应特征空间的正交基 $\{v_1, \ldots, v_r\}$，$P_i = \sum_{k=1}^r v_k v_k^*$。

    约化到特征空间：$r$ 个扰动后的特征值 $\lambda_{i_1}(t), \ldots, \lambda_{i_r}(t)$ 的导数是 $r \times r$ 矩阵 $P_i A'(0) P_i$ 限制在特征空间上的特征值。总和为 $\operatorname{tr}(P_i A'(0) P_i) = \operatorname{tr}(P_i A'(0))$。 $\blacksquare$

!!! theorem "定理 47B.13 (奇异值的导数)"
    设 $A(t)$ 是 $m \times n$ 矩阵族，$\sigma_i(0)$ 是 $A(0) = A$ 的**单奇异值**（$\sigma_i > 0$），对应左、右奇异向量 $u_i, v_i$（$Av_i = \sigma_i u_i$，$A^*u_i = \sigma_i v_i$）。则

    $$\sigma_i'(0) = \operatorname{Re}(u_i^* A'(0) v_i).$$

    对实矩阵：$\sigma_i'(0) = u_i^T A'(0) v_i$。

??? proof "证明"
    方法一：利用特征值与奇异值的关系。$\sigma_i^2$ 是 $A^*A$（或 $AA^*$）的特征值。

    $\dfrac{d(\sigma_i^2)}{dt} = 2\sigma_i \sigma_i'$。

    另一方面，$\sigma_i^2$ 是 Hermite 矩阵 $B(t) = A(t)^*A(t)$ 的特征值，对应特征向量 $v_i$。由定理 47B.12：

    $$\frac{d(\sigma_i^2)}{dt}\bigg|_0 = v_i^* B'(0) v_i = v_i^*(A'^* A + A^*A')v_i = v_i^*A'^*Av_i + v_i^*A^*A'v_i.$$

    $v_i^*A^*A'v_i = (Av_i)^*A'v_i = \sigma_i u_i^*A'v_i$。

    $v_i^*A'^*Av_i = (A'v_i)^*(Av_i) = \sigma_i(A'v_i)^*u_i = \sigma_i \overline{u_i^*A'v_i}$。

    因此 $2\sigma_i\sigma_i' = \sigma_i(u_i^*A'v_i + \overline{u_i^*A'v_i}) = 2\sigma_i\operatorname{Re}(u_i^*A'v_i)$。

    消去 $\sigma_i$（$\sigma_i > 0$）：$\sigma_i' = \operatorname{Re}(u_i^*A'v_i)$。 $\blacksquare$

!!! example "例 47B.5 (特征值灵敏度)"
    设 $A = \begin{pmatrix} 1 & 1000 \\ 0 & 2 \end{pmatrix}$，$E = \begin{pmatrix} 0 & 0 \\ \varepsilon & 0 \end{pmatrix}$。

    $A$ 的特征值 $\lambda_1 = 1, \lambda_2 = 2$。右特征向量 $v_1 = \binom{1}{0}$，$v_2 = \binom{1000}{1}$。左特征向量 $w_1^T = (1, -1000)/1$，$w_2^T = (0, 1)/1$。（归一化使 $w_i^*v_i = 1$。）

    $\lambda_1$ 的导数：$w_1^*Ev_1 = (1, -1000)\binom{0}{\varepsilon} = -1000\varepsilon$。

    虽然扰动 $E$ 很小，$\lambda_1$ 的变化被放大了 $1000$ 倍——这是因为 $A$ 高度非正规（特征向量 $v_1, v_2$ 几乎平行）。

---

## 47B.6 矩阵分解的导数

<div class="context-flow" markdown>

**核心问题**：经典矩阵分解（Cholesky、QR、LU）的因子如何随矩阵的扰动而变化？

</div>

矩阵分解的导数在自动微分（automatic differentiation）和灵敏度分析中至关重要。如果我们的计算流程中包含矩阵分解步骤，反向传播需要通过分解"求导"。

!!! theorem "定理 47B.14 (Cholesky 分解的导数)"
    设 $A = LL^T$ 是对称正定矩阵的 Cholesky 分解（$L$ 下三角，对角线元素为正）。则 $L$ 对 $A$ 的 Fréchet 导数 $dL$ 满足

    $$dA = dL \cdot L^T + L \cdot dL^T.$$

    设 $S = L^{-1} dA \, L^{-T}$，$\Phi = L^{-1} dL$。则

    $$S = \Phi + \Phi^T,$$

    其中 $\Phi$ 是下三角矩阵。因此

    $$\Phi = \operatorname{tril}(S) - \frac{1}{2}\operatorname{diag}(S),$$

    即 $\Phi$ 的严格下三角部分等于 $S$ 的严格下三角部分，对角线部分等于 $S$ 对角线的一半。

    因此 $dL = L\Phi = L\left(\operatorname{tril}(S) - \frac{1}{2}\operatorname{diag}(S)\right)$，其中 $S = L^{-1}dA\,L^{-T}$。

??? proof "证明"
    由 $A = LL^T$，取微分：$dA = dL \cdot L^T + L \cdot dL^T$。

    两侧左乘 $L^{-1}$，右乘 $L^{-T}$：

    $$S := L^{-1} dA \, L^{-T} = L^{-1}dL + (L^{-1}dL)^T = \Phi + \Phi^T.$$

    $\Phi = L^{-1}dL$ 是下三角矩阵（因为 $L$ 下三角，$dL$ 也下三角，$L^{-1}$ 也下三角）。

    将 $S = \Phi + \Phi^T$ 分解为下三角与上三角部分：$S$ 的严格下三角 $= \Phi$ 的严格下三角，$S$ 的对角线 $= 2\Phi$ 的对角线，$S$ 的严格上三角 $= \Phi^T$ 的严格上三角。

    因此 $\Phi_{ij} = S_{ij}$（$i > j$），$\Phi_{ii} = S_{ii}/2$。 $\blacksquare$

!!! theorem "定理 47B.15 (QR 分解的导数)"
    设 $A = QR$ 是 $m \times n$ 矩阵（$m \geq n$，$A$ 列满秩）的薄 QR 分解（$Q \in \mathbb{R}^{m \times n}$ 列正交，$R \in \mathbb{R}^{n \times n}$ 上三角，对角线为正）。则

    $$dA = dQ \cdot R + Q \cdot dR.$$

    设 $\Omega = Q^T dQ$（反对称矩阵，因为 $d(Q^TQ) = 0 \Rightarrow \Omega + \Omega^T = 0$），$M = Q^T dA \, R^{-1} = \Omega + dR \cdot R^{-1}$。

    则 $\Omega$ 是 $M$ 的严格下三角减去严格上三角：$\Omega = \operatorname{tril}(M, -1) - \operatorname{triu}(M, 1)$，$dR \cdot R^{-1}$ 是 $M$ 的上三角部分（含对角线）加上对角线以下的负转置：

    $$dR = (\operatorname{triu}(M) + \operatorname{diag}(M)/2 ...)$$

    更精确地：$dR = \operatorname{triu}(Q^TdA) \cdot$... 由于公式较复杂，直接给出结论：

    $$dR = \operatorname{triu}(R^{-T}(Q^TdA)^T + Q^TdA R^{-1})... $$

    简化表述：设 $W = Q^TdA$，则 $W = \Omega R + dR$。$\Omega$ 反对称，$dR$ 上三角。从 $WR^{-1} = \Omega + dRR^{-1}$ 中，$\Omega$ 的元素可由 $WR^{-1}$ 的反对称部分确定，$dR$ 由上三角部分确定。

!!! theorem "定理 47B.16 (LU 分解的导数)"
    设 $A = LU$ 是 $n \times n$ 可逆矩阵的 LU 分解（$L$ 单位下三角，$U$ 上三角）。假设所有顺序主子式非零。则

    $$dA = dL \cdot U + L \cdot dU.$$

    设 $M = L^{-1}dA\,U^{-1}$。则

    $$M = L^{-1}dL + dU \cdot U^{-1},$$

    其中 $L^{-1}dL$ 严格下三角（$L$ 的对角线为 $1$，故 $dL$ 对角线为 $0$），$dU \cdot U^{-1}$ 上三角。

    因此 $L^{-1}dL = \operatorname{tril}(M, -1)$，$dU \cdot U^{-1} = \operatorname{triu}(M)$。

    $dL = L \cdot \operatorname{tril}(M, -1)$，$dU = \operatorname{triu}(M) \cdot U$。

??? proof "证明"
    $dA = dL \cdot U + L \cdot dU$。左乘 $L^{-1}$，右乘 $U^{-1}$：

    $$M := L^{-1}dA\,U^{-1} = L^{-1}dL + dU\,U^{-1}.$$

    $L^{-1}dL$：$L$ 单位下三角 $\Rightarrow$ $L^{-1}$ 单位下三角 $\Rightarrow$ $L^{-1}dL$ 严格下三角。

    $dU \cdot U^{-1}$：$U$ 上三角 $\Rightarrow$ $U^{-1}$ 上三角 $\Rightarrow$ $dU \cdot U^{-1}$ 上三角。

    严格下三角 + 上三角 = 全矩阵。由于分解唯一（下三角部分和上三角部分不重叠，除对角线归属于上三角），可以直接读出。 $\blacksquare$

!!! example "例 47B.6 (Cholesky 分解的反向传播)"
    在机器学习中，设损失函数 $\ell$ 依赖于 Cholesky 因子 $L$，即 $\ell = \ell(L(A))$。

    已知 $\bar{L} = \dfrac{\partial \ell}{\partial L}$（由后续层反向传播得到），求 $\bar{A} = \dfrac{\partial \ell}{\partial A}$。

    由链式法则 $d\ell = \operatorname{tr}(\bar{L}^T dL)$，以及 $dL = L\Phi$（$\Phi = \operatorname{tril}(L^{-1}dA\,L^{-T}) - \frac{1}{2}\operatorname{diag}(L^{-1}dA\,L^{-T})$），

    $d\ell = \operatorname{tr}(\bar{L}^T L \Phi)$。设 $P = L^T\bar{L}$，则 $d\ell = \operatorname{tr}(P\Phi)$。

    由对称性可推导出 $\bar{A} = L^{-T}(\bar{P})L^{-1}$，其中 $\bar{P} = \operatorname{tril}(P) + \operatorname{tril}(P)^T - \operatorname{diag}(P)$... 具体公式取决于约定，但核心思路是通过定理 47B.14 的微分关系来"反转"梯度传播。

---

## 47B.7 条件数

<div class="context-flow" markdown>

**核心问题**：矩阵函数对输入扰动有多敏感？

</div>

!!! definition "定义 47B.4 (矩阵函数的条件数)"
    矩阵函数 $f$ 在 $A$ 处的**（相对）条件数**定义为

    $$\operatorname{cond}(f, A) = \lim_{\varepsilon \to 0} \sup_{\|E\| \leq \varepsilon} \frac{\|f(A + E) - f(A)\| / \|f(A)\|}{\|E\| / \|A\|} = \frac{\|L_f(A)\| \cdot \|A\|}{\|f(A)\|},$$

    其中 $\|L_f(A)\| = \sup_{\|E\| = 1} \|L_f(A)[E]\|$ 是 Fréchet 导数作为线性映射的**算子范数**。

    条件数衡量了 $A$ 的相对扰动导致 $f(A)$ 的相对扰动的最大放大倍数。

!!! theorem "定理 47B.17 (矩阵求逆的条件数)"
    $$\operatorname{cond}(A^{-1}, A) = \kappa(A),$$

    其中 $\kappa(A) = \|A\| \cdot \|A^{-1}\|$ 是矩阵的条件数。

??? proof "证明"
    $L_{A^{-1}}(A)[E] = -A^{-1}EA^{-1}$。

    $\|L_{A^{-1}}(A)\| = \sup_{\|E\|=1}\|A^{-1}EA^{-1}\| \leq \|A^{-1}\|^2$。

    取 $E = \|A^{-1}\|^{-1} u v^*$，其中 $u, v$ 使得 $\|A^{-1}E A^{-1}\| = \|A^{-1}\|^2$（可以选取使 $A^{-1}$ 达到范数的方向），得 $\|L\| = \|A^{-1}\|^2$。

    因此

    $$\operatorname{cond}(A^{-1}, A) = \frac{\|A^{-1}\|^2 \cdot \|A\|}{\|A^{-1}\|} = \|A^{-1}\| \cdot \|A\| = \kappa(A).$$

    即**矩阵求逆的条件数等于矩阵本身的条件数**。 $\blacksquare$

!!! theorem "定理 47B.18 (矩阵指数的条件数)"
    $$\operatorname{cond}(e^A, A) = \frac{\|L_{e^A}(A)\| \cdot \|A\|}{\|e^A\|}.$$

    对**正规矩阵** $A$（特征值 $\lambda_1, \ldots, \lambda_n$）：

    $$\|L_{e^A}(A)\| = \max_{i,j} |e[\lambda_i, \lambda_j]| = \max_{i,j} \left|\frac{e^{\lambda_i} - e^{\lambda_j}}{\lambda_i - \lambda_j}\right|,$$

    其中约定 $e[\lambda, \lambda] = e^{\lambda}$。

    对实对称矩阵 $A$（实特征值 $\alpha_1 \geq \cdots \geq \alpha_n$）：$\|e^A\| = e^{\alpha_1}$，$\|L_{e^A}\| = e^{\alpha_1}$（当所有特征值非负时），因此 $\operatorname{cond}(e^A, A) = \|A\|$。

!!! theorem "定理 47B.19 (矩阵平方根的条件数)"
    设 $A$ 为正定矩阵，特征值 $\lambda_1 \geq \cdots \geq \lambda_n > 0$。

    $$\operatorname{cond}(A^{1/2}, A) = \frac{1}{2}\frac{\sqrt{\lambda_1}}{\sqrt{\lambda_n}} \cdot \frac{\lambda_1}{\sqrt{\lambda_1}} \cdot ... $$

    更精确地：$\|L_{A^{1/2}}\| = \dfrac{1}{2\sqrt{\lambda_n}}$（最坏情况），因此

    $$\operatorname{cond}(A^{1/2}, A) = \frac{1}{2\sqrt{\lambda_n}} \cdot \frac{\lambda_1}{\sqrt{\lambda_1}} = \frac{\sqrt{\lambda_1}}{2\sqrt{\lambda_n}} = \frac{1}{2}\sqrt{\kappa(A)}.$$

??? proof "证明"
    $L_{A^{1/2}}(A)[E] = X$，其中 $SX + XS = E$（$S = A^{1/2}$）。

    $S$ 正定，特征值 $\sqrt{\lambda_i}$。Sylvester 方程的解可以表示为

    $$X = \int_0^{\infty} e^{-tS} E \, e^{-tS} dt.$$

    $\|X\| \leq \int_0^{\infty} \|e^{-tS}\|^2 dt \cdot \|E\| = \int_0^{\infty} e^{-2t\sqrt{\lambda_n}} dt \cdot \|E\| = \frac{1}{2\sqrt{\lambda_n}} \|E\|$。

    等号可以达到（取 $E = v_n v_n^T$，$v_n$ 是 $\lambda_n$ 的特征向量）。

    因此 $\|L_{A^{1/2}}\| = \dfrac{1}{2\sqrt{\lambda_n}}$。

    $\operatorname{cond}(A^{1/2}, A) = \dfrac{1}{2\sqrt{\lambda_n}} \cdot \dfrac{\|A\|}{\|A^{1/2}\|} = \dfrac{1}{2\sqrt{\lambda_n}} \cdot \dfrac{\lambda_1}{\sqrt{\lambda_1}} = \dfrac{\sqrt{\lambda_1}}{2\sqrt{\lambda_n}} = \dfrac{1}{2}\sqrt{\kappa(A)}$。 $\blacksquare$

!!! example "例 47B.7 (条件数的实际意义)"
    设 $A = \operatorname{diag}(1, 10^{-8})$，计算 $A^{1/2} = \operatorname{diag}(1, 10^{-4})$。

    $\kappa(A) = 10^8$，$\operatorname{cond}(A^{1/2}, A) = \frac{1}{2}\sqrt{10^8} = 5000$。

    $A$ 的 $10^{-16}$（双精度机器精度）的相对扰动导致 $A^{1/2}$ 约 $5000 \times 10^{-16} = 5 \times 10^{-13}$ 的相对误差——仍然可以接受。

    但若 $A = \operatorname{diag}(1, 10^{-16})$，则 $\operatorname{cond}(A^{1/2}, A) = \frac{1}{2} \times 10^8$，$A^{1/2}$ 的相对误差可达 $5 \times 10^{-9}$，失去约 $8$ 位有效数字。

---

## 47B.8 高阶 Fréchet 导数

<div class="context-flow" markdown>

**核心问题**：Fréchet 导数的"导数"是什么？高阶 Fréchet 导数在矩阵 Newton 法中扮演什么角色？

</div>

!!! definition "定义 47B.5 (二阶 Fréchet 导数)"
    设 $f$ 在 $A$ 处 Fréchet 可微，$L_f(A)$ 也关于 $A$ Fréchet 可微。$f$ 在 $A$ 处的**二阶 Fréchet 导数**是双线性映射 $L_f^{(2)}(A): \mathbb{C}^{n \times n} \times \mathbb{C}^{n \times n} \to \mathbb{C}^{n \times n}$，满足

    $$L_f(A + F)[E] = L_f(A)[E] + L_f^{(2)}(A)[F, E] + o(\|F\|),$$

    即 $L_f^{(2)}(A)[F, E] = \lim_{t \to 0} \dfrac{L_f(A + tF)[E] - L_f(A)[E]}{t}$。

    等价地，$f$ 的二阶 Taylor 展开：

    $$f(A + E) = f(A) + L_f(A)[E] + \frac{1}{2}L_f^{(2)}(A)[E, E] + O(\|E\|^3).$$

!!! theorem "定理 47B.20 (经典矩阵函数的二阶 Fréchet 导数)"
    1. **矩阵求逆**：$L^{(2)}_{A^{-1}}(A)[E, F] = A^{-1}EA^{-1}FA^{-1} + A^{-1}FA^{-1}EA^{-1}$。
    2. **矩阵指数**：$L^{(2)}_{e^A}(A)[E, F] = \int_0^1 \int_0^s e^{uA} E e^{(s-u)A} F e^{(1-s)A} du\,ds + (E \leftrightarrow F)$。

??? proof "证明"
    **(1)** $L_{A^{-1}}(A)[E] = -A^{-1}EA^{-1}$。对 $A$ 取 Fréchet 导数：

    $$\frac{d}{dt}\big|_0 L_{(A+tF)^{-1}}(A+tF)[E]$$
    $$= \frac{d}{dt}\big|_0 \left(-(A+tF)^{-1}E(A+tF)^{-1}\right)$$
    $$= A^{-1}FA^{-1}EA^{-1} + A^{-1}EA^{-1}FA^{-1}.$$

    这就是 $L^{(2)}_{A^{-1}}(A)[F, E]$。注意它关于 $E, F$ 是对称的（双线性映射的对称性）。 $\blacksquare$

!!! example "例 47B.8 (矩阵 Newton 法)"
    求解矩阵方程 $F(X) = 0$（如 $X^2 - A = 0$ 求矩阵平方根）。

    Newton 迭代：$X_{k+1} = X_k - (L_F(X_k))^{-1}[F(X_k)]$。

    对 $X^2 = A$：$F(X) = X^2 - A$，$L_F(X)[H] = XH + HX$。

    Newton 步：求 $H$ 使 $X_k H + H X_k = -(X_k^2 - A)$，然后 $X_{k+1} = X_k + H$。

    这就是 Sylvester 方程！每步 Newton 迭代需要求解一个 Sylvester 方程。

    **收敛性**：由于二阶 Fréchet 导数是有界的，Newton 法在解的邻域内二次收敛：

    $$\|X_{k+1} - A^{1/2}\| \leq C \|X_k - A^{1/2}\|^2.$$

---

## 47B.9 链式法则与应用

<div class="context-flow" markdown>

**核心问题**：如何将 Fréchet 导数的链式法则应用于矩阵方程的灵敏度分析和深度学习的反向传播？

</div>

!!! theorem "定理 47B.21 (Fréchet 导数的链式法则)"
    设 $h = g \circ f$，其中 $f$ 和 $g$ 是矩阵函数。则

    $$L_h(A)[E] = L_g(f(A))\big[L_f(A)[E]\big].$$

    即复合函数的 Fréchet 导数是各步 Fréchet 导数的**复合**（不是矩阵乘法，因为 Fréchet 导数是线性映射/算子）。

!!! theorem "定理 47B.22 (矩阵隐函数定理)"
    设隐式矩阵方程 $F(X, P) = 0$（$X$ 是未知矩阵，$P$ 是参数矩阵），$F$ 关于 $X$ 的 Fréchet 导数 $L_F^X$ 在解 $(X_0, P_0)$ 处作为线性算子可逆。则在 $(X_0, P_0)$ 的邻域内，$X$ 可表示为 $P$ 的可微函数 $X = X(P)$，且

    $$L_X(P_0)[dP] = -\left(L_F^X(X_0, P_0)\right)^{-1}\left[L_F^P(X_0, P_0)[dP]\right],$$

    其中 $L_F^X$ 和 $L_F^P$ 分别是 $F$ 对 $X$ 和 $P$ 的偏 Fréchet 导数。

??? proof "证明"
    由 $F(X(P), P) = 0$，两侧对 $P$ 取 Fréchet 导数（利用链式法则）：

    $$L_F^X[L_X(P)[dP]] + L_F^P[dP] = 0.$$

    由 $L_F^X$ 可逆：

    $$L_X(P)[dP] = -(L_F^X)^{-1}[L_F^P[dP]]. \quad \blacksquare$$

!!! example "例 47B.9 (Sylvester 方程的灵敏度分析)"
    考虑 Sylvester 方程 $AX + XB = C$，其解 $X = X(A, B, C)$。

    **对 $A$ 的灵敏度**：$F(X; A) = AX + XB - C$。

    $L_F^X[dX] = A \, dX + dX \cdot B$（这是 Sylvester 算子 $\mathcal{T}(dX) = A\,dX + dX B$）。

    $L_F^A[dA] = dA \cdot X$。

    由隐函数定理：$\mathcal{T}(dX) = -dA \cdot X$，即 $A \, dX + dX \cdot B = -dA \cdot X$。

    Vec 化：$(I \otimes A + B^T \otimes I)\operatorname{vec}(dX) = -(X^T \otimes I)\operatorname{vec}(dA)$。

    $$\frac{\partial \operatorname{vec}(X)}{\partial \operatorname{vec}(A)^T} = -(I \otimes A + B^T \otimes I)^{-1}(X^T \otimes I).$$

    条件数：$\operatorname{cond}(X, A) = \dfrac{\|(I \otimes A + B^T \otimes I)^{-1}\| \cdot \|X\| \cdot \|A\|}{\|X\|} = \dfrac{\|\mathcal{T}^{-1}\| \cdot \|X\| \cdot \|A\|}{\|X\|}$。

    当 $\sigma(A)$ 和 $\sigma(-B)$ 的距离 $\operatorname{sep}(A, -B) = \|\mathcal{T}^{-1}\|^{-1}$ 很小时，条件数很大。

!!! example "例 47B.10 (Lyapunov 方程的灵敏度)"
    离散 Lyapunov 方程 $AXA^T - X = -Q$（$Q$ 正定）。

    $F(X; A, Q) = AXA^T - X + Q = 0$。

    对 $A$：$L_F^A[dA] = dA \cdot X A^T + AX(dA)^T$。

    $\mathcal{L}(dX) := A\,dX\,A^T - dX$（Stein 算子）。

    $\mathcal{L}(dX) = -dA \cdot XA^T - AX\,dA^T$。

    条件数与 $\rho(A)$（$A$ 的谱半径）密切相关：当 $\rho(A) \to 1$ 时（系统接近不稳定），条件数趋于无穷。

!!! example "例 47B.11 (反向传播的矩阵视角)"
    考虑深度学习中的单层前向传播：

    $$Z = XW + B, \quad H = \sigma(Z), \quad L = \ell(H, Y),$$

    其中 $X$ 是输入，$W$ 是权重，$B$ 是偏置，$\sigma$ 是逐元素激活函数，$\ell$ 是损失函数。

    **前向传播**（正序计算）：$X \to Z \to H \to L$。

    **反向传播**（逆序计算 Fréchet 导数）：

    1. $\bar{H} = \dfrac{\partial L}{\partial H}$（损失对激活的导数）。
    2. $\bar{Z} = \bar{H} \odot \sigma'(Z)$（逐元素函数的 Fréchet 导数）。
    3. $\bar{W} = X^T \bar{Z}$（线性函数的 Fréchet 导数的伴随）。

    在 Fréchet 导数框架中，每一步都是 Fréchet 导数的伴随算子：

    - $L_\ell(H)[dH] = \operatorname{tr}(\bar{H}^T dH)$（标量对矩阵的 Fréchet 导数）。
    - $L_\sigma(Z)[dZ] = \sigma'(Z) \odot dZ$（逐元素函数）。
    - $L_{Z}(W)[dW] = X \, dW$（线性函数）。

    **伴随（adjoint）**：$L_{Z}(W)^*[\bar{Z}] = X^T \bar{Z}$，这正是反向传播中 $\bar{W}$ 的计算公式。

    链式法则 $\bar{W} = L_Z^*[L_\sigma^*[\bar{H}]] = X^T(\bar{H} \odot \sigma'(Z))$ 就是反向传播公式。

!!! example "例 47B.12 (总结：矩阵微积分体系)"
    ```
    第 47A 章（基础工具）：
      标量 → 向量：梯度 ∇f
      向量 → 向量：Jacobi 矩阵 J_f
      标量 → 矩阵：∂f/∂X
      矩阵 → 矩阵：Vec 算子 + Kronecker 积

    第 47B 章（Fréchet 导数理论）：
      矩阵函数的导数：Fréchet 导数 L_f(A)[E]
      核心定理：Daleckii-Krein（除商差分表示）
      灵敏度分析：条件数 cond(f, A) = ‖L_f‖·‖A‖/‖f(A)‖
      应用：矩阵方程、自动微分、Newton 法

    统一工具：
      - 微分形式 df = tr((∂f/∂X)ᵀ dX)
      - 链式法则：L_{g∘f}(A)[E] = L_g(f(A))[L_f(A)[E]]
      - 隐函数定理：矩阵方程灵敏度的统一框架

    核心应用：
      - 最优化（梯度下降、Newton 法）
      - 机器学习（反向传播、可微编程）
      - 灵敏度分析（条件数、扰动界）
      - 矩阵方程的扰动理论
    ```

---

## 练习

!!! exercise "习题 47B.1"
    设 $f(A) = A^3$。

    (a) 利用乘积规则求 $L_f(A)[E]$。

    (b) 利用 Daleckii-Krein 定理验证你的结果（假设 $A$ 可对角化）。

!!! exercise "习题 47B.2"
    设 $f(A) = A^{-2}$。利用链式法则（$A^{-2} = (A^{-1})^2$ 或 $A^{-2} = (A^2)^{-1}$）求 $L_f(A)[E]$，并验证两种方法给出相同结果。

!!! exercise "习题 47B.3"
    设 $A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}$，$E = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}$。

    (a) 计算 $L_{e^A}(A)[E]$（利用积分公式或 Daleckii-Krein 定理）。

    (b) 计算 $L_{A^{-1}}(A)[E]$。

    (c) 用数值方法验证：$\dfrac{e^{A+0.001E} - e^A}{0.001} \approx L_{e^A}(A)[E]$。

!!! exercise "习题 47B.4"
    证明矩阵对数的条件数公式：对正定矩阵 $A$（特征值 $\lambda_1 \geq \lambda_n > 0$），

    $$\operatorname{cond}(\log A, A) = \frac{\max(|\log\lambda_1|, |\log\lambda_n|)}{\max_{i,j}\frac{|\log\lambda_i - \log\lambda_j|}{|\lambda_i - \lambda_j|} \cdot \lambda_1}...$$

    提示：利用 Daleckii-Krein 定理计算 $\|L_{\log}\|$。

!!! exercise "习题 47B.5"
    设 $A(t) = \begin{pmatrix} 1+t & t \\ t & 1-t \end{pmatrix}$。

    (a) 求 $A(0)$ 的特征值和特征向量。

    (b) 利用定理 47B.12 求 $\lambda_i'(0)$。

    (c) 直接计算 $A(t)$ 的特征值 $\lambda_i(t)$（解特征方程），并验证 (b) 的结果。

!!! exercise "习题 47B.6"
    设 $A = \begin{pmatrix} 4 & 2 \\ 2 & 3 \end{pmatrix}$（正定），$A = LL^T$ 的 Cholesky 分解为 $L = \begin{pmatrix} 2 & 0 \\ 1 & \sqrt{2} \end{pmatrix}$。

    设 $dA = \begin{pmatrix} 0.01 & 0 \\ 0 & 0.01 \end{pmatrix}$。

    利用定理 47B.14 计算 $dL$。

!!! exercise "习题 47B.7"
    证明：对 Hermite 矩阵族 $A(t)$，$\lambda_{\max}(A(t))$ 是凸函数（假设 $A(t)$ 关于 $t$ 仿射）。

    提示：利用 $\lambda_{\max}(A) = \max_{\|x\|=1} x^*Ax$ 和凸函数的逐点上确界仍是凸函数。

!!! exercise "习题 47B.8"
    设 Sylvester 方程 $AX + XB = C$ 的条件数为 $\kappa_{\text{Syl}} = \dfrac{\|A\| + \|B\|}{\operatorname{sep}(A, -B)}$，其中 $\operatorname{sep}(A, -B) = \min_{\|X\|_F=1}\|AX+XB\|_F$。

    (a) 证明 $\operatorname{sep}(A, -B) = \sigma_{\min}(I \otimes A + B^T \otimes I)$。

    (b) 对 $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$，$B = \begin{pmatrix} -1.01 & 0 \\ 0 & -3 \end{pmatrix}$，计算 $\operatorname{sep}(A, -B)$ 并评估方程的条件。

!!! exercise "习题 47B.9"
    利用二阶 Fréchet 导数，证明矩阵 Newton 迭代

    $$X_{k+1} = X_k - (L_F(X_k))^{-1}[F(X_k)]$$

    在解 $X^*$ 的邻域内二次收敛：$\|X_{k+1} - X^*\| \leq C\|X_k - X^*\|^2$。

    提示：将 $F(X^*)$ 在 $X_k$ 处 Taylor 展开到二阶。

!!! exercise "习题 47B.10"
    设 $f(A) = e^A$，$g(A) = \log A$。验证链式法则：$L_{g \circ f}(A)[E] = E$（因为 $\log(e^A) = A$），即

    $$L_{\log}(e^A)\left[\int_0^1 e^{sA}Ee^{(1-s)A}ds\right] = E.$$

    这说明 $L_{\log}(e^A)$ 是 $L_{e^A}(A)$ 的逆算子。

!!! exercise "习题 47B.11"
    **(综合题)** 考虑矩阵 Riccati 方程 $A^TX + XA - XBR^{-1}B^TX + Q = 0$，其中 $A, B, Q, R$ 为适当维度的矩阵，$Q, R$ 对称正定。

    (a) 写出 $F(X; A) = A^TX + XA - XBR^{-1}B^TX + Q$ 对 $X$ 的 Fréchet 导数 $L_F^X[dX]$。

    (b) 写出对 $A$ 的 Fréchet 导数 $L_F^A[dA]$。

    (c) 由隐函数定理，$dX$ 满足什么方程？这个方程的结构是什么？

    (d) 讨论该方程条件数与闭环矩阵 $A - BR^{-1}B^TX$ 的稳定性之间的关系。

## 练习题

1. **[概念] 什么是 Fréchet 导数？它与 Gâteaux 导数的主要区别是什么？**
   ??? success "参考答案"
       Fréchet 导数是一个线性算子 $L_f(A)$，满足 $f(A+E) = f(A) + L_f(A)[E] + o(\|E\|)$。Gâteaux 导数仅要求沿特定方向的极限存在，而 Fréchet 导数要求在全方向上的一致逼近。Fréchet 可微必导致 Gâteaux 可微，且导数相等。

2. **[求逆] 利用定义证明矩阵求逆 $f(A) = A^{-1}$ 的 Fréchet 导数为 $L_{A^{-1}}(A)[E] = -A^{-1}EA^{-1}$。**
   ??? success "参考答案"
       $(A+E)^{-1} - A^{-1} = (A(I+A^{-1}E))^{-1} - A^{-1} = (I+A^{-1}E)^{-1}A^{-1} - A^{-1}$。
       利用级数展开 $(I+X)^{-1} \approx I - X$，得：
       $(I - A^{-1}E + O(\|E\|^2))A^{-1} - A^{-1} = -A^{-1}EA^{-1} + O(\|E\|^2)$。
       线性主部即为所求。

3. **[Daleckii-Krein] 若 $A = \operatorname{diag}(\lambda_1, \lambda_2)$，写出矩阵指数 $e^A$ 的 Fréchet 导数在方向 $E = \begin{pmatrix} e_{11} & e_{12} \\ e_{21} & e_{22} \end{pmatrix}$ 上的作用。**
   ??? success "参考答案"
       根据公式 $L_f(A)[E] = f^{[1]}(\Lambda) \odot E$：
       $L_{e^A}(A)[E] = \begin{pmatrix} e^{\lambda_1} e_{11} & \frac{e^{\lambda_1}-e^{\lambda_2}}{\lambda_1-\lambda_2} e_{12} \\ \frac{e^{\lambda_2}-e^{\lambda_1}}{\lambda_2-\lambda_1} e_{21} & e^{\lambda_2} e_{22} \end{pmatrix}$。

4. **[特征值] 证明单特征值 $\lambda$ 的导数满足 $d\lambda = \mathbf{w}^H (dA) \mathbf{v}$，其中 $\mathbf{v}, \mathbf{w}$ 分别是右和左特征向量。**
   ??? success "参考答案"
       对 $A\mathbf{v} = \lambda \mathbf{v}$ 求微分：$(dA)\mathbf{v} + A(d\mathbf{v}) = (d\lambda)\mathbf{v} + \lambda(d\mathbf{v})$。
       两侧左乘 $\mathbf{w}^H$，利用 $\mathbf{w}^H A = \lambda \mathbf{w}^H$：
       $\mathbf{w}^H (dA) \mathbf{v} + \lambda \mathbf{w}^H d\mathbf{v} = d\lambda (\mathbf{w}^H \mathbf{v}) + \lambda \mathbf{w}^H d\mathbf{v}$。
       消去相同项并设 $\mathbf{w}^H \mathbf{v} = 1$ 即得结果。

5. **[平方根] 矩阵平方根 $X = A^{1/2}$ 的导数 $dX$ 满足什么方程？**
   ??? success "参考答案"
       满足 Sylvester 方程：$X(dX) + (dX)X = dA$。这是由于对 $X^2 = A$ 两侧求微分得到的。

6. **[条件数] 为什么矩阵求逆的条件数 $\operatorname{cond}(A^{-1}, A)$ 正好等于矩阵本身的条件数 $\kappa(A)$？**
   ??? success "参考答案"
       因为 $\|L_{A^{-1}}\| = \|A^{-1}\|^2$（算子范数）。根据条件数定义：
       $\frac{\|L\| \cdot \|A\|}{\|f(A)\|} = \frac{\|A^{-1}\|^2 \cdot \|A\|}{\|A^{-1}\|} = \|A^{-1}\| \cdot \|A\| = \kappa(A)$。

7. **[链式法则] 设 $h(A) = \exp(A^2)$。利用链式法则写出其 Fréchet 导数。**
   ??? success "参考答案"
       设 $f(A) = A^2$, $g(X) = e^X$。则 $L_h(A)[E] = L_g(A^2)[L_f(A)[E]]$。
       其中 $L_f(A)[E] = AE + EA$，将其作为整体代入指数函数的积分或除商公式中。

8. **[分解] 在 Cholesky 分解 $A = LL^T$ 中，已知 $dA$，求 $dL$ 的基本思路是什么？**
   ??? success "参考答案"
       思路是利用 $dA = (dL)L^T + L(dL^T)$。通过左乘 $L^{-1}$ 和右乘 $L^{-T}$ 将其对称化，然后利用 $L^{-1}dL$ 是下三角矩阵的性质，提取出对称矩阵的下三角部分（对角线减半）。

9. **[SVD导数] 简述奇异值 $\sigma_i$ 的导数公式。**
   ??? success "参考答案"
       $d\sigma_i = \operatorname{Re}(\mathbf{u}_i^H (dA) \mathbf{v}_i)$，其中 $\mathbf{u}_i, \mathbf{v}_i$ 是对应的左右奇异向量。这说明奇异值的变化只取决于扰动在奇异向量方向上的投影。

10. **[爱因斯坦思考题] 爱因斯坦在发展广义相对论时，使用了黎曼几何中的联络（Connection）来描述度规张量的变化。Fréchet 导数在算子空间中扮演了类似的角色：它描述了当算子“坐标”（矩阵元素）微动时，算子整体（如指数映射）如何“弯曲”。为什么在研究如矩阵 Riccati 方程等非线性算子方程时，我们必须使用 Fréchet 导数而非简单的偏导数？**
    ??? success "参考答案"
        因为非线性算子方程（如 $A^T X + XA - XBR^{-1}B^T X + Q = 0$）中的变量 $X$ 处于一个具有特定结构的流形或 Banach 空间中。简单的偏导数只能处理“分量”的变化，而 Fréchet 导数捕捉的是算子作为**整体线性映射**的内在灵敏度。它不依赖于具体的矩阵表示，能够直接揭示系统稳定性（如闭环极点）与解的变化率之间的几何联系。这正是爱因斯坦“几何即物理”思想在算子理论中的体现：系统的微分结构决定了它的物理响应。

## 本章小结

本章将微积分的触角延伸至算子空间，建立了严谨的矩阵函数导数理论——Fréchet 导数，主要内容包括：

1. **Fréchet 导数体系**：定义了作为线性算子的 Fréchet 导数，并对比了其与方向导数（Gâteaux 导数）的层级关系。
2. **运算律与经典函数**：推导了求逆、指数、对数及平方根等核心矩阵函数的导数公式，揭示了指数函数的积分表示与平方根的 Sylvester 方程表示。
3. **Daleckii-Krein 定理**：给出了基于特征值除商差分的统一导数表达形式，这是矩阵分析中连接解析性质与谱性质的枢纽。
4. **扰动与灵敏度**：详细论述了特征值、奇异值以及各种矩阵分解（LU, QR, Cholesky）的微分性质，确立了矩阵函数条件数作为衡量不稳定性量度的权威地位。
5. **现代计算连接**：展示了 Fréchet 导数如何在隐函数定理下统一处理矩阵方程的灵敏度，并为深度学习中的矩阵反向传播提供了坚实的算子代数基础。

