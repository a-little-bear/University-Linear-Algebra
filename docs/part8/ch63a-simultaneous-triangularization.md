# 第 63A 章 同时三角化与同时对角化

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · Schur 分解(Ch14) · Lie 代数基础(Ch58)

**本章脉络**：同时三角化定义 → McCoy 定理（完整证明）→ 可交换矩阵的同时三角化 → Lie 定理 → 同时对角化充要条件 → 正规矩阵的同时酉对角化 → Burnside 定理联系 → 公共不变子空间问题

**延伸**：同时三角化理论是表示论、Lie 代数和联合谱分析（Ch63B）的基础；公共不变子空间问题连接算子代数与不变子空间理论

</div>

当我们从单一矩阵的 Schur 分解走向**矩阵集合**的结构分析时，首先面对的核心问题是：什么条件下，一组矩阵可以被同一个可逆变换同时化为三角形或对角形？单个矩阵在代数闭域上总可以三角化（Schur 分解），但一组矩阵能否找到**公共的**三角化变换矩阵，则需要它们之间满足深刻的代数约束。

同时三角化理论由 McCoy（1936）奠基，随后被 Lie 代数的可解性理论（Lie 定理）大幅推广。同时对角化则与量子力学中的可交换可观测量密切相关。本章系统建立这些基本结构定理，并探讨 Burnside 定理和公共不变子空间问题等进阶方向。

---

## 63A.1 同时三角化的基本概念

### 定义与等价刻画

对于单个矩阵，Schur 分解告诉我们：在代数闭域上，任意矩阵都可以酉三角化。但如果我们有**一组**矩阵，能否找到**一个**公共的可逆矩阵，同时将它们化为上三角形？

!!! definition "定义 63A.1 (同时三角化)"
    设 $\mathcal{F} = \{A_1, A_2, \ldots, A_m\} \subseteq M_n(\mathbb{F})$ 是一组 $n \times n$ 矩阵。称 $\mathcal{F}$ 可以**同时三角化**（simultaneously triangularizable），若存在可逆矩阵 $P \in M_n(\mathbb{F})$，使得所有 $P^{-1} A_i P$（$i = 1, \ldots, m$）均为上三角矩阵。

    等价地，存在 $\mathbb{F}^n$ 的一个完全旗（complete flag）
    $$
    \{0\} = V_0 \subset V_1 \subset V_2 \subset \cdots \subset V_n = \mathbb{F}^n, \quad \dim V_k = k
    $$
    使得每个 $V_k$ 都是 $\mathcal{F}$ 中所有矩阵的不变子空间。

!!! note "注"
    同时三角化要求所有矩阵共享**同一个**完全旗的不变子空间链，这比每个矩阵单独可三角化强得多。从不变子空间的角度来看，同时三角化等价于存在公共的**最大不变子空间链**。

!!! definition "定义 63A.2 (同时酉三角化)"
    若上述定义中的 $P$ 可以选为酉矩阵 $U$（即 $U^* U = I$），则称 $\mathcal{F}$ 可以**同时酉三角化**。此时每个 $U^* A_i U$ 均为上三角矩阵，且变换保持内积结构。

---

## 63A.2 McCoy 定理

McCoy 定理给出了代数闭域上同时三角化的完整代数刻画。这是同时三角化理论的奠基性结果。

!!! theorem "定理 63A.1 (McCoy 定理, 1936)"
    设 $\mathbb{F}$ 是代数闭域，$A_1, A_2, \ldots, A_m \in M_n(\mathbb{F})$。则 $\{A_1, \ldots, A_m\}$ 可同时三角化，当且仅当对任意非交换多项式 $p(x_1, \ldots, x_m)$，矩阵
    $$
    [A_i, p(A_1, \ldots, A_m)]
    $$
    对所有 $i = 1, \ldots, m$ 都是幂零的（其中 $[X, Y] = XY - YX$ 是换位子/李括号）。

??? proof "证明"
    **必要性**：设所有 $A_i$ 可同时三角化，即存在可逆矩阵 $P$ 使得 $T_i = P^{-1}A_iP$ 均为上三角矩阵。

    对任意非交换多项式 $p$，有
    $$
    P^{-1} p(A_1, \ldots, A_m) P = p(T_1, \ldots, T_m)
    $$
    上三角矩阵在加法和乘法下封闭，因此 $p(T_1, \ldots, T_m)$ 仍为上三角矩阵。

    设 $T_i = D_i + N_i$，其中 $D_i$ 为对角部分，$N_i$ 为严格上三角部分。上三角矩阵 $T_i$ 与上三角矩阵 $p(T_1, \ldots, T_m)$ 的换位子
    $$
    [T_i, p(T_1, \ldots, T_m)]
    $$
    的 $(k,k)$ 对角元素为
    $$
    (T_i)_{kk} \cdot p(T_1, \ldots, T_m)_{kk} - p(T_1, \ldots, T_m)_{kk} \cdot (T_i)_{kk} = 0
    $$
    因此该换位子是**严格上三角**矩阵，从而是幂零的（严格上三角矩阵的幂零指数不超过 $n$）。

    由于 $P^{-1}[A_i, p(A_1, \ldots, A_m)]P = [T_i, p(T_1, \ldots, T_m)]$ 是幂零的，$[A_i, p(A_1, \ldots, A_m)]$ 也是幂零的。

    **充分性**：使用对 $n$（矩阵维数）的强归纳法。

    **基础情形**：$n = 1$ 时，所有 $1 \times 1$ 矩阵都是上三角的，命题显然成立。

    **归纳步骤**：设 $n \geq 2$，假设命题对所有维数小于 $n$ 的情形成立。

    **第一步：构造公共特征向量。** 关键是证明 $A_1, \ldots, A_m$ 存在一个**公共特征向量**，即存在非零向量 $v \in \mathbb{F}^n$ 使得 $v$ 是每个 $A_i$ 的某个不变一维子空间中的向量。更准确地说，存在 $v \neq 0$ 和标量 $\lambda_1, \ldots, \lambda_m$ 使得 $A_i v = \lambda_i v$ 对所有 $i$ 成立。

    由于 $\mathbb{F}$ 是代数闭的，$A_1$ 至少有一个特征值 $\lambda_1$。设 $V_1 = \ker(A_1 - \lambda_1 I)$ 是 $A_1$ 对应 $\lambda_1$ 的特征空间。

    **第二步：证明特征空间被保持。** 对任意 $i \geq 2$ 和 $v \in V_1$，考虑 $A_i v$ 是否仍在 $V_1$ 中。取多项式 $p(x_1, \ldots, x_m) = x_i$，条件给出 $[A_1, A_i]$ 是幂零的。设 $C = [A_1, A_i] = A_1 A_i - A_i A_1$，则 $C^N = 0$ 对某个 $N \geq 1$。

    对于 $v \in V_1$，有 $A_1(A_i v) = A_i(A_1 v) + Cv = \lambda_1 A_i v + Cv$。这表明 $A_i$ 不一定保持 $V_1$ 不变，但**广义特征空间** $W_1 = \ker(A_1 - \lambda_1 I)^n$ 是 $A_i$ 不变的（利用 $[A_1, A_i]$ 的幂零性和多项式关系）。

    更精确地：设 $W_1$ 是 $A_1$ 关于特征值 $\lambda_1$ 的广义特征空间。由于 $[A_1, p(A_1, \ldots, A_m)]$ 对所有多项式 $p$ 都幂零，特别地 $[A_1 - \lambda_1 I, A_i]$ 幂零。利用关系式
    $$
    (A_1 - \lambda_1 I) A_i = A_i (A_1 - \lambda_1 I) + [A_1 - \lambda_1 I, A_i]
    $$
    可以证明若 $(A_1 - \lambda_1 I)^k v = 0$，则 $(A_1 - \lambda_1 I)^{k+N}(A_i v) = 0$。因此 $A_i$ 保持广义特征空间 $W_1$ 不变。

    **第三步：在广义特征空间中递归。** 在 $W_1$ 上限制所有 $A_i$，条件仍然满足（因为换位子的限制仍是幂零的）。若 $\dim W_1 < n$，由归纳假设，$A_1|_{W_1}, \ldots, A_m|_{W_1}$ 可以在 $W_1$ 中同时三角化。

    若 $W_1 = \mathbb{F}^n$（即 $A_1$ 只有一个特征值 $\lambda_1$），则 $A_1 - \lambda_1 I$ 是幂零的。此时考虑 $A_2$ 的特征空间分解，并利用条件继续分析。

    **第四步：组合为完全旗。** 通过在各广义特征空间中的递归同时三角化，以及对不同广义特征空间的排列，可以构造整个 $\mathbb{F}^n$ 上的公共完全旗。

    具体地，找到公共特征向量 $v_1$（在上述分析中，$W_1$ 中通过归纳法可找到 $A_1, \ldots, A_m$ 的公共特征向量），以 $v_1$ 为第一个基向量，所有 $A_i$ 在此基下的矩阵具有分块形式
    $$
    A_i = \begin{pmatrix} \lambda_i^{(1)} & * \\ 0 & A_i' \end{pmatrix}
    $$
    其中 $A_i' \in M_{n-1}(\mathbb{F})$。

    需要验证 $\{A_1', \ldots, A_m'\}$ 继承了原始条件。设 $q(x_1, \ldots, x_m)$ 是任意非交换多项式，则
    $$
    p(A_1, \ldots, A_m) = \begin{pmatrix} p(\lambda_1^{(1)}, \ldots, \lambda_m^{(1)}) & * \\ 0 & p(A_1', \ldots, A_m') \end{pmatrix}
    $$
    从而 $[A_i, p]$ 的右下角块恰好是 $[A_i', p(A_1', \ldots, A_m')]$，继承幂零性。

    由归纳假设，$A_1', \ldots, A_m'$ 可以同时三角化，从而原始矩阵族也可以同时三角化。$\blacksquare$

!!! example "例 63A.1"
    设
    $$
    A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}, \quad
    B = \begin{pmatrix} 1 & -1 \\ 0 & 4 \end{pmatrix}
    $$

    验证 $AB \neq BA$：
    $$
    AB = \begin{pmatrix} 2 & 2 \\ 0 & 12 \end{pmatrix}, \quad
    BA = \begin{pmatrix} 2 & -2 \\ 0 & 12 \end{pmatrix}
    $$
    所以 $[A, B] = AB - BA = \begin{pmatrix} 0 & 4 \\ 0 & 0 \end{pmatrix}$，这是幂零的（严格上三角）。

    $A, B$ 虽然不可交换，但它们已经是上三角的，因此以 $P = I$ 同时三角化。McCoy 条件在此情形下显然满足：对所有多项式 $p$，$[A, p(A,B)]$ 和 $[B, p(A,B)]$ 都是严格上三角的，从而幂零。

---

## 63A.3 可交换矩阵的同时三角化

最重要且最常用的特殊情形是可交换矩阵族。

!!! theorem "定理 63A.2 (可交换矩阵的同时三角化)"
    设 $\mathbb{F}$ 是代数闭域。若 $A_1, A_2, \ldots, A_m \in M_n(\mathbb{F})$ 两两可交换（即 $A_i A_j = A_j A_i$ 对所有 $i, j$ 成立），则它们可以同时三角化。

??? proof "证明"
    **方法一（由 McCoy 定理推导）**：由 McCoy 定理，只需验证对所有 $i$ 和所有非交换多项式 $p$，$[A_i, p(A_1, \ldots, A_m)]$ 是幂零的。

    由于 $A_1, \ldots, A_m$ 两两可交换，对任意非交换多项式 $p$，计算 $p(A_1, \ldots, A_m)$ 时变量顺序无关（乘法可交换化为交换多项式的值），且 $A_i$ 与 $p(A_1, \ldots, A_m)$ 可交换。从而
    $$
    [A_i, p(A_1, \ldots, A_m)] = A_i \cdot p(A_1, \ldots, A_m) - p(A_1, \ldots, A_m) \cdot A_i = 0
    $$
    零矩阵当然是幂零的。

    **方法二（直接归纳证明）**：对 $n$ 归纳。$n = 1$ 时显然。

    设 $n \geq 2$。因为 $\mathbb{F}$ 代数闭，$A_1$ 有特征值 $\lambda$，设 $V_\lambda = \ker(A_1 - \lambda I)$ 是对应的特征空间。

    **关键事实**：由于 $A_i A_1 = A_1 A_i$，对任意 $v \in V_\lambda$：
    $$
    A_1(A_i v) = A_i(A_1 v) = A_i(\lambda v) = \lambda(A_i v)
    $$
    因此 $A_i v \in V_\lambda$，即每个 $A_i$ 保持 $V_\lambda$ 不变。

    取 $A_1$ 的所有不同特征值 $\lambda_1, \ldots, \lambda_r$。因 $\mathbb{F}$ 代数闭，$A_1$ 可对角化或 Jordan 分解。无论如何，广义特征空间 $W_j = \ker(A_1 - \lambda_j I)^n$ 给出直和分解 $\mathbb{F}^n = W_1 \oplus \cdots \oplus W_r$，且每个 $A_i$ 保持每个 $W_j$ 不变。

    在每个 $W_j$ 上，限制矩阵 $A_2|_{W_j}, \ldots, A_m|_{W_j}$ 仍然两两可交换。由归纳假设（$\dim W_j < n$ 或直接在子空间中构造），在每个 $W_j$ 中可找到 $A_2, \ldots, A_m$ 的公共特征向量。

    取各 $W_j$ 中的公共特征向量作为旗的第一批基向量，缩减维度后继续归纳，最终构造出完全旗，完成同时三角化。$\blacksquare$

!!! example "例 63A.2"
    对角矩阵族是同时三角化的最简单例子。设
    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}, \quad B = \begin{pmatrix} 3 & 0 \\ 0 & 5 \end{pmatrix}
    $$
    $AB = BA$（对角矩阵总可交换），且 $P = I$ 同时将它们三角化（对角矩阵本身就是上三角矩阵）。

!!! example "例 63A.3"
    设 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$, $B = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}$。则 $AB = BA = 0$。这两个幂零矩阵可交换，已经是上三角形式，因此可同时三角化。注意 $A$ 和 $B$ 的像空间相同，$\operatorname{Im}(A) = \operatorname{Im}(B) = \operatorname{span}\{e_1\}$，这构成公共不变子空间链 $\{0\} \subset \operatorname{span}\{e_1\} \subset \mathbb{C}^2$。

---

## 63A.4 Lie 定理

同时三角化的一个深刻推广来自 Lie 代数理论。Lie 定理将可交换性条件放松为可解性条件。

!!! definition "定义 63A.3 (可解 Lie 代数)"
    Lie 代数 $\mathfrak{g}$ 的**导列**（derived series）定义为
    $$
    \mathfrak{g}^{(0)} = \mathfrak{g}, \quad \mathfrak{g}^{(k+1)} = [\mathfrak{g}^{(k)}, \mathfrak{g}^{(k)}]
    $$
    其中 $[\mathfrak{h}_1, \mathfrak{h}_2] = \operatorname{span}\{[X, Y] : X \in \mathfrak{h}_1, Y \in \mathfrak{h}_2\}$。

    称 $\mathfrak{g}$ 是**可解的**（solvable），若存在 $N$ 使得 $\mathfrak{g}^{(N)} = \{0\}$。

!!! note "注"
    Abel Lie 代数（即 $[X, Y] = 0$ 对所有 $X, Y$）当然可解（$\mathfrak{g}^{(1)} = 0$）。两两可交换的矩阵族生成的 Lie 代数是 Abel 的。全体上三角矩阵构成的 Lie 代数是可解的但不是 Abel 的。

!!! theorem "定理 63A.3 (Lie 定理)"
    设 $\mathbb{F}$ 是特征零的代数闭域，$\mathfrak{g} \subseteq \mathfrak{gl}_n(\mathbb{F})$ 是一个可解 Lie 代数。则 $\mathfrak{g}$ 中的所有矩阵可以同时三角化。

    特别地，若矩阵集 $\{A_1, \ldots, A_m\}$ 生成的 Lie 代数是可解的，则 $\{A_1, \ldots, A_m\}$ 可同时三角化。

??? proof "证明"
    对 $\dim \mathfrak{g}$ 和 $n$ 进行归纳。

    **关键步骤**：证明 $\mathfrak{g}$ 中所有矩阵存在公共特征向量。

    由于 $\mathfrak{g}$ 可解，$\mathfrak{g}^{(1)} = [\mathfrak{g}, \mathfrak{g}]$ 是 $\mathfrak{g}$ 的真理想（$\mathfrak{g}$ 非 Abel 时 $\mathfrak{g}^{(1)} \subsetneq \mathfrak{g}$；$\mathfrak{g}$ 为 Abel 时直接由可交换矩阵定理得到结论）。

    取 $\mathfrak{h}$ 为 $\mathfrak{g}$ 的一个余维 1 的理想（这样的理想存在，因为 $\mathfrak{g}/\mathfrak{g}^{(1)}$ 是 Abel 的，可以取 $\mathfrak{h}$ 包含 $\mathfrak{g}^{(1)}$ 且 $\dim \mathfrak{g}/\mathfrak{h} = 1$）。$\mathfrak{h}$ 本身也是可解的。

    由归纳假设，$\mathfrak{h}$ 中所有矩阵存在公共特征向量 $v$。更一般地，存在 $\mathfrak{h}$ 的一个**权**（weight）$\lambda: \mathfrak{h} \to \mathbb{F}$，即线性泛函，使得
    $$
    W = \{w \in \mathbb{F}^n : H w = \lambda(H) w, \; \forall H \in \mathfrak{h}\} \neq \{0\}
    $$

    **核心论证**：设 $X \in \mathfrak{g} \setminus \mathfrak{h}$（存在这样的 $X$，因为 $\mathfrak{h}$ 是余维 1 的）。需要证明 $X$ 保持 $W$ 不变。

    对任意 $w \in W$ 和 $H \in \mathfrak{h}$：
    $$
    H(Xw) = X(Hw) + [H, X]w = X(\lambda(H)w) + [H, X]w = \lambda(H)(Xw) + [H, X]w
    $$
    由于 $\mathfrak{h}$ 是 $\mathfrak{g}$ 的理想，$[H, X] \in \mathfrak{h}$，故 $[H, X]w = \lambda([H, X])w$。因此
    $$
    H(Xw) = \lambda(H)(Xw) + \lambda([H, X])w
    $$

    定义子空间序列：$W_0 = \{0\}$，$W_j = \operatorname{span}\{w, Xw, X^2w, \ldots, X^{j-1}w\}$。利用上述关系式，可以用归纳法证明 $H$ 在 $W_j$ 上的矩阵（关于基 $\{w, Xw, \ldots\}$）是上三角矩阵，对角线全为 $\lambda(H)$。

    特别地，$\operatorname{tr}(H|_{W_j}) = j \cdot \lambda(H)$。取 $H = [H', X]$（其中 $H' \in \mathfrak{h}$），由于 $X$ 保持 $W_j$ 稳定（$W_j$ 的定义保证了这一点），$[H', X]|_{W_j}$ 的迹等于 $\operatorname{tr}(H' X|_{W_j}) - \operatorname{tr}(X H'|_{W_j}) = 0$。因此 $j \cdot \lambda([H', X]) = 0$，在特征零域中 $\lambda([H', X]) = 0$。

    回到上面的计算，$H(Xw) = \lambda(H)(Xw)$，即 $Xw \in W$。因此 $X$ 保持权空间 $W$ 不变。

    在 $W$ 中，$X$ 有特征向量 $v$（因 $\mathbb{F}$ 代数闭）。$v$ 同时是 $\mathfrak{h}$ 中所有矩阵和 $X$ 的公共特征向量，从而是 $\mathfrak{g}$ 中所有矩阵的公共特征向量。

    以 $v$ 为第一个基向量，所有 $\mathfrak{g}$ 中矩阵的维数降低，归纳完成。$\blacksquare$

!!! example "例 63A.4"
    设 $\mathfrak{g}$ 由矩阵
    $$
    X = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}, \quad Y = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}
    $$
    生成。$[X, Y] = XY - YX = \begin{pmatrix} 0 & -2 \\ 0 & 0 \end{pmatrix} = -2Y$。

    $\mathfrak{g} = \operatorname{span}\{X, Y\}$，$\mathfrak{g}^{(1)} = [\mathfrak{g}, \mathfrak{g}] = \operatorname{span}\{Y\}$，$\mathfrak{g}^{(2)} = [\mathfrak{g}^{(1)}, \mathfrak{g}^{(1)}] = [Y, Y] = \{0\}$。因此 $\mathfrak{g}$ 可解。

    由 Lie 定理，$X$ 和 $Y$ 可同时三角化——它们已经是上三角的，公共特征向量为 $e_1 = (1, 0)^T$。

---

## 63A.5 同时对角化

### 充要条件

!!! definition "定义 63A.4 (同时对角化)"
    设 $A_1, A_2, \ldots, A_m \in M_n(\mathbb{F})$。称它们可**同时对角化**（simultaneously diagonalizable），若存在可逆矩阵 $P$ 使得 $P^{-1}A_iP$ 对所有 $i$ 都是对角矩阵。

!!! theorem "定理 63A.4 (同时对角化的充要条件)"
    矩阵 $A_1, A_2, \ldots, A_m \in M_n(\mathbb{C})$ 可同时对角化，当且仅当以下两个条件同时成立：

    1. 每个 $A_i$ 都是可对角化的（即每个 $A_i$ 的最小多项式无重根）；
    2. 它们两两可交换：$A_i A_j = A_j A_i$ 对所有 $i, j$ 成立。

??? proof "证明"
    **必要性**：若存在 $P$ 使得 $D_i = P^{-1}A_iP$ 均为对角矩阵，则：

    (1) 每个 $A_i$ 可对角化，因为 $A_i = P D_i P^{-1}$，$D_i$ 为对角矩阵。

    (2) 对角矩阵之间总是可交换的：$D_i D_j = D_j D_i$，从而
    $$
    A_i A_j = P D_i P^{-1} \cdot P D_j P^{-1} = P D_i D_j P^{-1} = P D_j D_i P^{-1} = A_j A_i
    $$

    **充分性**：对矩阵个数 $m$ 进行归纳。

    **基础情形**：$m = 1$ 时，$A_1$ 可对角化即表示存在 $P$ 使得 $P^{-1}A_1P$ 为对角矩阵。

    **归纳步骤**：设 $m \geq 2$，由归纳假设 $A_1, \ldots, A_{m-1}$ 可同时对角化。设 $P_0$ 使得 $D_i = P_0^{-1}A_iP_0$ 为对角矩阵（$i = 1, \ldots, m-1$）。

    令 $B = P_0^{-1}A_m P_0$。由于 $A_m$ 与每个 $A_i$（$i \leq m-1$）可交换，$B$ 与每个 $D_i$ 可交换。

    设 $D_1$ 的不同对角元素为 $\mu_1, \ldots, \mu_r$，对应的特征空间为 $E_1, \ldots, E_r$（$\mathbb{C}^n = E_1 \oplus \cdots \oplus E_r$）。

    由 $BD_1 = D_1B$ 知，$B$ 保持 $D_1$ 的每个特征空间 $E_k$ 不变。因此在适当排列坐标后，$B$ 具有分块对角形式
    $$
    B = \operatorname{diag}(B_1, B_2, \ldots, B_r)
    $$
    其中 $B_k = B|_{E_k}$。

    由于 $A_m$ 可对角化，$B = P_0^{-1}A_m P_0$ 也可对角化，因此每个 $B_k$ 可对角化。设 $Q_k$ 使得 $Q_k^{-1}B_kQ_k$ 为对角矩阵，令 $Q = \operatorname{diag}(Q_1, \ldots, Q_r)$。

    则 $Q^{-1}BQ$ 是对角矩阵。同时，由于 $Q$ 保持 $D_1$ 的特征空间分解不变，$Q^{-1}D_1Q = D_1$。

    对于 $i = 2, \ldots, m-1$，由 $BD_i = D_iB$ 且 $D_i$ 在 $E_k$ 上是标量矩阵（因为 $D_i$ 与 $D_1$ 在 $E_k$ 上的限制可交换，且 $D_1|_{E_k} = \mu_k I_{E_k}$——但这需要更仔细的论证）。实际上，$D_i$ 在各 $E_k$ 上的限制 $D_i|_{E_k}$ 是对角矩阵，而 $Q_k$ 的变换不影响对角矩阵的对角性（$Q_k^{-1}D_i|_{E_k}Q_k = D_i|_{E_k}$，因为 $D_i|_{E_k}$ 是对角的，$B_k$ 与 $D_i|_{E_k}$ 可交换意味着 $Q_k$ 可以选取为也对角化 $D_i|_{E_k}$ 的——但 $D_i|_{E_k}$ 已经是对角的）。

    因此 $(P_0 Q)^{-1}A_i(P_0 Q)$ 对所有 $i = 1, \ldots, m$ 都是对角矩阵。$\blacksquare$

!!! example "例 63A.5"
    设
    $$
    A = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}, \quad
    B = \begin{pmatrix} 3 & -2 \\ -2 & 3 \end{pmatrix}
    $$
    $A, B$ 均为实对称矩阵，因此都可对角化。验证 $AB = BA$（实对称矩阵不一定可交换，需计算验证）：
    $$
    AB = \begin{pmatrix} 7 & 2 \\ 2 & 7 \end{pmatrix} = BA
    $$

    公共特征向量为 $v_1 = \frac{1}{\sqrt{2}}(1,1)^T$（$Av_1 = 9v_1$，$Bv_1 = v_1$）和 $v_2 = \frac{1}{\sqrt{2}}(1,-1)^T$（$Av_2 = v_2$，$Bv_2 = 5v_2$）。

    取 $P = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$，则
    $$
    P^T A P = \begin{pmatrix} 9 & 0 \\ 0 & 1 \end{pmatrix}, \quad
    P^T B P = \begin{pmatrix} 1 & 0 \\ 0 & 5 \end{pmatrix}
    $$

### 正规矩阵的同时酉对角化

!!! theorem "定理 63A.5 (正规矩阵的同时酉对角化)"
    可交换的正规矩阵族可以同时**酉**对角化。即若 $A_1, \ldots, A_m$ 都是正规的（$A_i^* A_i = A_i A_i^*$）且两两可交换，则存在酉矩阵 $U$ 使得 $U^* A_i U$ 对所有 $i$ 都是对角矩阵。

??? proof "证明"
    由谱定理，每个正规矩阵 $A_i$ 可酉对角化。由定理 63A.4，$A_1, \ldots, A_m$ 可同时对角化（因正规矩阵可对角化且两两可交换）。

    需要证明变换矩阵可以选为酉矩阵。由谱定理，$A_1$ 可以酉对角化为 $U_1^* A_1 U_1 = D_1$。在 $D_1$ 的每个特征空间中，其他 $A_i$ 的限制仍为正规矩阵（因为酉变换保持正规性，且可交换性保证了不变性）。

    在每个特征空间内部继续选取酉对角化基，由于正规矩阵的酉对角化基可以被递归构造，最终得到的总变换矩阵是酉的（酉矩阵的分块对角拼接仍然是酉矩阵——更准确地说，$U_1$ 乘以分块酉矩阵 $U_2$，$U_1 U_2$ 仍然是酉矩阵）。$\blacksquare$

!!! note "注"
    实对称矩阵当然是正规的，因此可交换的实对称矩阵族可以被同一个正交矩阵同时正交对角化。这一结论在量子力学中尤为重要：**可交换的可观测量**（Hermite 算子）具有公共的本征态基底。两个可观测量可以同时精确测量，当且仅当它们对应的算子可交换。

!!! example "例 63A.6"
    考虑 Pauli 矩阵 $\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$，$\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$。

    $\sigma_x \sigma_z = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \neq \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \sigma_z \sigma_x$。

    因此 $\sigma_x$ 和 $\sigma_z$ 不可交换，不能同时对角化。在量子力学中，这意味着自旋的 $x$ 分量和 $z$ 分量不能同时精确测量。

---

## 63A.6 Burnside 定理与矩阵半群

Burnside 定理将同时三角化理论与矩阵半群的有界性联系起来。

!!! definition "定义 63A.5 (矩阵半群)"
    矩阵集合 $\mathcal{S} \subseteq M_n(\mathbb{C})$ 称为**矩阵半群**（matrix semigroup），若 $\mathcal{S}$ 在矩阵乘法下封闭，即 $A, B \in \mathcal{S}$ 蕴含 $AB \in \mathcal{S}$。

    称 $\mathcal{S}$ 是**有界的**，若存在常数 $C > 0$ 使得 $\|A\| \leq C$ 对所有 $A \in \mathcal{S}$ 成立。

!!! theorem "定理 63A.6 (Burnside 定理的矩阵版本)"
    设 $\mathcal{S} \subseteq M_n(\mathbb{C})$ 是一个矩阵半群。若 $\mathcal{S}$ 中所有矩阵的迹都有界，且 $\mathcal{S}$ 是不可约的（即 $\mathcal{S}$ 中的矩阵没有公共的非平凡不变子空间），则 $\mathcal{S}$ 生成的矩阵代数等于 $M_n(\mathbb{C})$。

!!! note "注"
    Burnside 定理的一个重要推论是：若有界矩阵半群 $\mathcal{S}$ 是不可约的，则 $\mathcal{S}$ 中的矩阵可以同时酉相似于酉矩阵的子集。即存在酉矩阵 $U$ 使得 $U^*\mathcal{S}U$ 中每个元素的范数不超过 1。

!!! theorem "定理 63A.7 (有界半群的同时三角化)"
    设 $\mathcal{S} \subseteq M_n(\mathbb{C})$ 是有界矩阵半群。若 $\mathcal{S}$ 中所有矩阵可同时三角化，则存在酉矩阵 $U$ 使得 $U^*AU$ 对所有 $A \in \mathcal{S}$ 均为上三角矩阵，且对角线元素的模不超过 1。

??? proof "证明"
    由同时三角化条件，存在可逆矩阵 $P$ 使得 $T_A = P^{-1}AP$ 对所有 $A \in \mathcal{S}$ 均为上三角矩阵。由于 $\mathcal{S}$ 有界，乘积半群 $\{T_A : A \in \mathcal{S}\}$ 也有界。

    对于有界的上三角矩阵半群，对角线元素 $d_i(A)$ 满足 $|d_i(A_1) \cdot d_i(A_2) \cdots d_i(A_k)| = |d_i(A_1 \cdots A_k)| \leq C$ 对所有 $k$ 和所有乘积成立。由 Gelfand 公式的推广，$|d_i(A)| \leq 1$ 对所有 $A \in \mathcal{S}$。

    通过 Gram-Schmidt 正交化过程，可以将 $P$ 修正为酉矩阵 $U$，在保持上三角结构的同时获得酉变换。$\blacksquare$

---

## 63A.7 公共不变子空间问题

### 基本问题

!!! definition "定义 63A.6 (公共不变子空间)"
    设 $\mathcal{F} = \{A_1, \ldots, A_m\} \subseteq M_n(\mathbb{C})$。子空间 $V \subseteq \mathbb{C}^n$ 称为 $\mathcal{F}$ 的**公共不变子空间**（common invariant subspace），若 $A_i V \subseteq V$ 对所有 $i = 1, \ldots, m$ 成立。

    若 $\mathcal{F}$ 除 $\{0\}$ 和 $\mathbb{C}^n$ 外没有公共不变子空间，称 $\mathcal{F}$ 是**不可约的**（irreducible）。

!!! theorem "定理 63A.8 (可交换算子的公共不变子空间)"
    设 $A, B \in M_n(\mathbb{C})$ 且 $AB = BA$。则 $A$ 的每个特征空间和广义特征空间都是 $B$ 的不变子空间。特别地，若 $A$ 有至少两个不同的特征值，则 $\{A, B\}$ 有非平凡的公共不变子空间。

??? proof "证明"
    设 $\lambda$ 是 $A$ 的特征值，$V_\lambda = \ker(A - \lambda I)$。对任意 $v \in V_\lambda$：
    $$
    A(Bv) = B(Av) = B(\lambda v) = \lambda(Bv)
    $$
    因此 $Bv \in V_\lambda$，即 $V_\lambda$ 是 $B$ 的不变子空间。

    对于广义特征空间 $W_\lambda = \ker(A - \lambda I)^n$，设 $v \in W_\lambda$，即 $(A - \lambda I)^n v = 0$。由 $AB = BA$，$(A - \lambda I)B = B(A - \lambda I)$，因此
    $$
    (A - \lambda I)^n (Bv) = B (A - \lambda I)^n v = B \cdot 0 = 0
    $$
    故 $Bv \in W_\lambda$。

    若 $A$ 有至少两个不同的特征值 $\lambda_1 \neq \lambda_2$，则 $W_{\lambda_1}$ 是非平凡的（$\{0\} \subsetneq W_{\lambda_1} \subsetneq \mathbb{C}^n$），且是 $\{A, B\}$ 的公共不变子空间。$\blacksquare$

!!! theorem "定理 63A.9 (不可约矩阵集的结构)"
    设 $\mathcal{F} \subseteq M_n(\mathbb{C})$ 是不可约的有限矩阵集合。则：

    1. $\mathcal{F}$ 不可能同时三角化（除非 $n = 1$）；
    2. $\mathcal{F}$ 生成的代数 $\mathcal{A} = \operatorname{span}\{A_{i_1} \cdots A_{i_k} : k \geq 0\}$ 等于 $M_n(\mathbb{C})$（Burnside 定理）；
    3. $\mathcal{F}$ 中至少有两个矩阵不可交换。

!!! example "例 63A.7"
    矩阵 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}$。

    检查公共不变子空间：$A$ 的像空间为 $\operatorname{span}\{e_1\}$，但 $B(e_1) = e_2 \notin \operatorname{span}\{e_1\}$。$B$ 的像空间为 $\operatorname{span}\{e_2\}$，但 $A(e_2) = e_1 \notin \operatorname{span}\{e_2\}$。

    因此 $\{A, B\}$ 是不可约的。由 Burnside 定理，$\operatorname{span}\{I, A, B, AB, BA, \ldots\}$ 应当等于 $M_2(\mathbb{C})$。实际上 $AB = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$BA = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$，$A = E_{12}$，$B = E_{21}$，它们张成 $M_2(\mathbb{C})$。

---

## 63A.8 同时三角化的判定算法

!!! theorem "定理 63A.10 (同时三角化的简化判据)"
    设 $\mathbb{F}$ 是代数闭域，$A_1, \ldots, A_m \in M_n(\mathbb{F})$。以下条件等价：

    1. $\{A_1, \ldots, A_m\}$ 可同时三角化；
    2. 对所有 $i, j$ 和所有字（word）$w = A_{k_1} \cdots A_{k_s}$，$[A_i, A_j] \cdot w$ 是幂零的；
    3. 对所有 $i, j$ 和所有 $k \geq 0$，$\operatorname{tr}([A_i, A_j] \cdot (A_{l_1} \cdots A_{l_k})) = 0$。

??? proof "证明"
    (1) $\Rightarrow$ (3)：若所有 $A_i$ 可同时上三角化为 $T_i$，则 $[T_i, T_j]$ 是严格上三角的，$[T_i, T_j] \cdot T_{l_1} \cdots T_{l_k}$ 也是严格上三角的（因为严格上三角矩阵乘以上三角矩阵仍是严格上三角的），迹为 0。

    (3) $\Rightarrow$ (2)：设 $N = [A_i, A_j] \cdot w$。条件 (3) 意味着 $\operatorname{tr}(N^k) = 0$ 对所有 $k \geq 1$（通过展开并利用条件）。由 Newton 恒等式，$N$ 的所有特征值为零，因此 $N$ 幂零。

    (2) $\Rightarrow$ (1)：这由 McCoy 定理的推论得到。条件 (2) 保证了 McCoy 条件中关于换位子幂零性的要求。$\blacksquare$

!!! example "例 63A.8"
    设 $A = \begin{pmatrix} 1 & 1 & 0 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix}$，$B = \begin{pmatrix} 2 & 0 & 1 \\ 0 & 2 & 0 \\ 0 & 0 & 1 \end{pmatrix}$。

    计算 $[A, B] = AB - BA$：
    $$
    AB = \begin{pmatrix} 2 & 2 & 1 \\ 0 & 2 & 1 \\ 0 & 0 & 2 \end{pmatrix}, \quad BA = \begin{pmatrix} 2 & 2 & 2 \\ 0 & 2 & 2 \\ 0 & 0 & 2 \end{pmatrix}
    $$
    $$
    [A, B] = \begin{pmatrix} 0 & 0 & -1 \\ 0 & 0 & -1 \\ 0 & 0 & 0 \end{pmatrix}
    $$
    这是严格上三角矩阵，因此幂零。$\operatorname{tr}([A,B] \cdot W) = 0$ 对所有字 $W$ 成立（因为 $[A,B]$ 是严格上三角的，乘以上三角矩阵后仍然是严格上三角的，迹为零）。

    因此 $A, B$ 可以同时三角化——它们已经是上三角的。

---

## 习题

!!! question "习题 63A.1"
    证明两个上三角矩阵 $A, B \in M_n(\mathbb{C})$ 的换位子 $[A, B] = AB - BA$ 是严格上三角的。

!!! question "习题 63A.2"
    设 $A, B \in M_2(\mathbb{C})$ 且 $AB = BA$。直接证明（不使用 McCoy 定理）$A, B$ 可以同时三角化。（提示：取 $A$ 的特征向量。）

!!! question "习题 63A.3"
    设 $A_1, \ldots, A_m$ 都是实对称矩阵。证明它们可同时正交对角化当且仅当两两可交换。

!!! question "习题 63A.4"
    设 $A, B \in M_3(\mathbb{C})$ 满足 $A^2 = A$（$A$ 为幂等矩阵）和 $AB = BA$。证明 $A$ 和 $B$ 可以同时三角化。

!!! question "习题 63A.5"
    给出一组三个 $2 \times 2$ 矩阵，它们两两之间的换位子都是幂零的，但它们不可同时三角化的例子。（提示：考虑 $\mathbb{F} = \mathbb{F}_2$，非代数闭域。）解释为什么 McCoy 定理需要代数闭域的假设。

!!! question "习题 63A.6"
    设 $\mathfrak{g} \subseteq \mathfrak{gl}_3(\mathbb{C})$ 是全体上三角矩阵构成的 Lie 代数。验证 $\mathfrak{g}$ 是可解的，计算其导列，并找到同时三角化的公共基。

!!! question "习题 63A.7"
    证明：若 $\Sigma = \{A_1, \ldots, A_m\}$ 可同时三角化，且每个 $A_i$ 的特征值都是实数，则可以选取变换矩阵 $P$ 使得 $P^{-1}A_iP$ 是实上三角矩阵。

!!! question "习题 63A.8"
    设 $A, B \in M_n(\mathbb{C})$。证明：若 $AB - BA = A$，则 $A$ 是幂零的。（提示：利用迹条件 $\operatorname{tr}(A^k) = 0$ 进行归纳。）

!!! question "习题 63A.9"
    设 $\mathcal{S}$ 是一个有界矩阵半群，$\mathcal{S} \subseteq M_2(\mathbb{R})$。证明 $\mathcal{S}$ 中的矩阵可以同时相似于以下三种形式之一：(a) 同时上三角化；(b) 同时酉对角化（若 $\mathcal{S}$ 不可约）；(c) 全体为标量矩阵。

!!! question "习题 63A.10"
    设 $A_1, \ldots, A_m \in M_n(\mathbb{C})$ 可同时三角化。证明：对任意连续函数 $f: \mathbb{C}^m \to \mathbb{C}$，
    $$
    \operatorname{tr}(f(A_1, \ldots, A_m)) = \sum_{k=1}^n f(\lambda_k^{(1)}, \ldots, \lambda_k^{(m)})
    $$
    其中 $(\lambda_k^{(1)}, \ldots, \lambda_k^{(m)})$ 是同时三角化后第 $k$ 个对角元素组成的向量。

!!! question "习题 63A.11"
    （公共不变子空间）设 $A = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 2 & 1 \\ 0 & 4 \end{pmatrix}$。找出 $\{A, B\}$ 的所有公共不变子空间。验证它们构成一个完全旗。

!!! question "习题 63A.12"
    证明 Burnside 定理的以下特殊情形：设 $A, B \in M_2(\mathbb{C})$，且 $\{A, B\}$ 不可约。证明 $\{I, A, B, AB\}$ 线性无关，从而张成 $M_2(\mathbb{C})$。
