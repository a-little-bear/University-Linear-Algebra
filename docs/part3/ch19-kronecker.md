# 第 19 章 Kronecker 积与 Vec 算子

<div class="context-flow" markdown>

**前置**：矩阵乘法/特征值(Ch6) · **脉络**：$A\otimes B$ 组合两个空间 → **Vec算子**将矩阵方程 $AXB=C$ 变为 $(B^T\otimes A)\operatorname{vec}(X)=\operatorname{vec}(C)$ → Ch20 Sylvester/Lyapunov求解

</div>

在线性代数的实际应用中，我们经常遇到需要将矩阵方程转化为向量方程的情形，或者需要构造具有特殊结构的大矩阵。**Kronecker 积**（Kronecker product）和 **Vec 算子**（vectorization operator）正是处理这类问题的核心工具。Kronecker 积提供了一种系统的方式来"组合"两个矩阵空间上的线性映射，而 Vec 算子则将矩阵"拉直"为向量，使得矩阵方程可以借助 Kronecker 积转化为标准的线性方程组。

本章从 Kronecker 积的定义和基本性质出发，介绍 Vec 算子及其核心公式，讨论置换矩阵（commutation matrix）的作用，然后展示如何利用这些工具求解矩阵方程，最后介绍 Kronecker 和及其与矩阵指数的关系。

---

## 19.1 Kronecker 积定义

<div class="context-flow" markdown>

**本质**：$A\otimes B$ 是 $mp\times nq$ 分块矩阵，第$(i,j)$块为 $a_{ij}B$ · 不满足交换律，但置换相似($B\otimes A = P(A\otimes B)P^T$)

</div>

!!! definition "定义 19.1 (Kronecker 积)"
    设 $A = (a_{ij})$ 为 $m \times n$ 矩阵，$B$ 为 $p \times q$ 矩阵。$A$ 与 $B$ 的 **Kronecker 积**（又称**张量积**，tensor product），记作 $A \otimes B$，定义为 $mp \times nq$ **分块矩阵**：
    $$
    A \otimes B = \begin{pmatrix}
    a_{11}B & a_{12}B & \cdots & a_{1n}B \\
    a_{21}B & a_{22}B & \cdots & a_{2n}B \\
    \vdots & \vdots & \ddots & \vdots \\
    a_{m1}B & a_{m2}B & \cdots & a_{mn}B
    \end{pmatrix}.
    $$

!!! note "注"
    Kronecker 积不满足交换律：一般 $A \otimes B \neq B \otimes A$，但它们是**置换相似**的（permutation similar），即存在置换矩阵 $P$ 使得 $B \otimes A = P(A \otimes B)P^T$。这个置换矩阵就是后面要讨论的置换矩阵（commutation matrix）。

!!! example "例 19.1"
    设 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix}$，则：
    $$
    A \otimes B = \begin{pmatrix}
    1 \cdot \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix} & 2 \cdot \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix} \\[6pt]
    3 \cdot \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix} & 4 \cdot \begin{pmatrix} 0 & 5 \\ 6 & 7 \end{pmatrix}
    \end{pmatrix}
    = \begin{pmatrix}
    0 & 5 & 0 & 10 \\
    6 & 7 & 12 & 14 \\
    0 & 15 & 0 & 20 \\
    18 & 21 & 24 & 28
    \end{pmatrix}.
    $$

---

## 19.2 Kronecker 积性质

<div class="context-flow" markdown>

**核心性质**：**混合积** $(A\otimes B)(C\otimes D)=(AC)\otimes(BD)$ → 推出逆、转置、行列式公式 · $\det(A\otimes B)=(\det A)^n(\det B)^m$

</div>

Kronecker 积具有丰富而优美的代数性质，使得它成为矩阵理论中不可或缺的工具。

!!! theorem "定理 19.1 (混合积性质)"
    设 $A, C$ 为可相乘的矩阵对，$B, D$ 为可相乘的矩阵对，则：
    $$
    (A \otimes B)(C \otimes D) = (AC) \otimes (BD).
    $$
    此性质称为**混合积性质**（mixed-product property）。

??? proof "证明"
    设 $A$ 为 $m \times n$，$B$ 为 $p \times q$，$C$ 为 $n \times r$，$D$ 为 $q \times s$。则 $A \otimes B$ 为 $mp \times nq$，$C \otimes D$ 为 $nq \times rs$，乘积 $(A \otimes B)(C \otimes D)$ 为 $mp \times rs$。

    $(A \otimes B)$ 的第 $(i,j)$ 块（$p \times q$ 大小）为 $a_{ij}B$，$(C \otimes D)$ 的第 $(j,k)$ 块（$q \times s$ 大小）为 $c_{jk}D$。

    乘积的第 $(i,k)$ 块为：
    $$
    \sum_{j=1}^{n} (a_{ij}B)(c_{jk}D) = \sum_{j=1}^{n} a_{ij}c_{jk}(BD) = \left(\sum_{j=1}^n a_{ij}c_{jk}\right)(BD) = (AC)_{ik}(BD).
    $$

    这正是 $(AC) \otimes (BD)$ 的第 $(i,k)$ 块。因此 $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$。$\blacksquare$

!!! theorem "定理 19.2 (Kronecker 积的基本代数性质)"
    设 $A, B, C$ 为适当大小的矩阵，$\alpha$ 为标量。则：

    1. **结合律**：$(A \otimes B) \otimes C = A \otimes (B \otimes C)$。
    2. **分配律**：$A \otimes (B + C) = A \otimes B + A \otimes C$，$(A + B) \otimes C = A \otimes C + B \otimes C$。
    3. **标量乘法**：$(\alpha A) \otimes B = A \otimes (\alpha B) = \alpha(A \otimes B)$。
    4. **转置**：$(A \otimes B)^T = A^T \otimes B^T$。
    5. **共轭转置**：$(A \otimes B)^* = A^* \otimes B^*$。
    6. **逆**：若 $A, B$ 均可逆，则 $(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$。

??? proof "证明"
    我们证明第 6 条。由混合积性质：
    $$
    (A \otimes B)(A^{-1} \otimes B^{-1}) = (AA^{-1}) \otimes (BB^{-1}) = I_m \otimes I_p = I_{mp}.
    $$

    类似地 $(A^{-1} \otimes B^{-1})(A \otimes B) = I_{mp}$。因此 $(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$。

    其他性质可由定义直接验证。$\blacksquare$

!!! theorem "定理 19.3 (Kronecker 积的迹、行列式与秩)"
    设 $A$ 为 $m \times m$ 矩阵，$B$ 为 $n \times n$ 矩阵。则：

    1. **迹**：$\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \cdot \operatorname{tr}(B)$。
    2. **行列式**：$\det(A \otimes B) = (\det A)^n (\det B)^m$。
    3. **秩**：$\operatorname{rank}(A \otimes B) = \operatorname{rank}(A) \cdot \operatorname{rank}(B)$。

??? proof "证明"
    **(1) 迹**：$A \otimes B$ 的对角块为 $a_{ii}B$（$i = 1,\ldots,m$），因此：
    $$
    \operatorname{tr}(A \otimes B) = \sum_{i=1}^m \operatorname{tr}(a_{ii}B) = \sum_{i=1}^m a_{ii} \operatorname{tr}(B) = \operatorname{tr}(A) \cdot \operatorname{tr}(B).
    $$

    **(2) 行列式**：利用混合积性质和分块对角化。设 $A$ 有特征值 $\lambda_1, \ldots, \lambda_m$，$B$ 有特征值 $\mu_1, \ldots, \mu_n$（计入重数）。由定理 19.5（后面将证明），$A \otimes B$ 的特征值为 $\{\lambda_i \mu_j\}$。因此：
    $$
    \det(A \otimes B) = \prod_{i=1}^m \prod_{j=1}^n \lambda_i \mu_j = \left(\prod_{i=1}^m \lambda_i\right)^n \left(\prod_{j=1}^n \mu_j\right)^m = (\det A)^n (\det B)^m.
    $$

    **(3) 秩**：设 $\operatorname{rank}(A) = r$，$\operatorname{rank}(B) = s$。存在可逆矩阵 $P, Q, U, V$ 使得 $A = P \begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} Q$，$B = U \begin{pmatrix} I_s & 0 \\ 0 & 0 \end{pmatrix} V$。

    则 $A \otimes B = (P \otimes U) \left(\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} \otimes \begin{pmatrix} I_s & 0 \\ 0 & 0 \end{pmatrix}\right) (Q \otimes V)$。

    中间矩阵通过直接计算可知其秩为 $rs$，而 $P \otimes U$ 和 $Q \otimes V$ 可逆，故 $\operatorname{rank}(A \otimes B) = rs$。$\blacksquare$

!!! example "例 19.2"
    设 $A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$。

    **迹**：$\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \cdot \operatorname{tr}(B) = 5 \times 3 = 15$。

    直接验证：$A \otimes B$ 的对角元素为 $2 \cdot 1, 2 \cdot 2, 3 \cdot 1, 3 \cdot 2 = 2, 4, 3, 6$，迹为 $2+4+3+6=15$。

    **行列式**：$\det(A \otimes B) = (\det A)^2 (\det B)^2 = 6^2 \cdot 2^2 = 36 \cdot 4 = 144$。

    **秩**：$\operatorname{rank}(A \otimes B) = 2 \times 2 = 4$（满秩）。

!!! definition "定义 19.2 (Kronecker 幂)"
    对方阵 $A$，定义 $k$ 次 **Kronecker 幂**为：
    $$
    A^{\otimes k} = \underbrace{A \otimes A \otimes \cdots \otimes A}_{k \text{ 个}}.
    $$
    若 $A$ 为 $n \times n$，则 $A^{\otimes k}$ 为 $n^k \times n^k$。

---

## 19.3 Vec 算子

<div class="context-flow" markdown>

**桥梁公式**：$\operatorname{vec}(AXB)=(B^T\otimes A)\operatorname{vec}(X)$ · 特例：$\operatorname{vec}(AX)=(I\otimes A)\operatorname{vec}(X)$，$\operatorname{vec}(XB)=(B^T\otimes I)\operatorname{vec}(X)$

</div>

Vec 算子将矩阵按列堆叠为向量，是连接矩阵方程与向量方程的桥梁。

!!! definition "定义 19.3 (Vec 算子)"
    设 $A = (\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n)$ 为 $m \times n$ 矩阵，其中 $\mathbf{a}_j$ 为第 $j$ 列。**Vec 算子**（vectorization）将 $A$ 映为 $mn \times 1$ 列向量：
    $$
    \operatorname{vec}(A) = \begin{pmatrix} \mathbf{a}_1 \\ \mathbf{a}_2 \\ \vdots \\ \mathbf{a}_n \end{pmatrix}.
    $$
    即将 $A$ 的各列从左到右依次堆叠。

!!! theorem "定理 19.4 (Vec 算子的核心公式)"
    设 $A$ 为 $m \times n$ 矩阵，$X$ 为 $n \times p$ 矩阵，$B$ 为 $p \times q$ 矩阵。则：
    $$
    \operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X).
    $$

??? proof "证明"
    **方法一**（利用列向量）：令 $B = (\mathbf{b}_1, \ldots, \mathbf{b}_q)$，$Y = AXB$，则 $Y$ 的第 $j$ 列为：
    $$
    \mathbf{y}_j = AX\mathbf{b}_j = A \sum_{k=1}^p b_{kj} \mathbf{x}_k = \sum_{k=1}^p b_{kj} A \mathbf{x}_k,
    $$
    其中 $\mathbf{x}_k$ 为 $X$ 的第 $k$ 列。

    因此：
    $$
    \operatorname{vec}(Y) = \begin{pmatrix} \mathbf{y}_1 \\ \vdots \\ \mathbf{y}_q \end{pmatrix} = \begin{pmatrix} \sum_k b_{k1} A\mathbf{x}_k \\ \vdots \\ \sum_k b_{kq} A\mathbf{x}_k \end{pmatrix}.
    $$

    另一方面：
    $$
    (B^T \otimes A)\operatorname{vec}(X) = \begin{pmatrix} b_{11}A & b_{21}A & \cdots & b_{p1}A \\ b_{12}A & b_{22}A & \cdots & b_{p2}A \\ \vdots & & \ddots & \vdots \\ b_{1q}A & b_{2q}A & \cdots & b_{pq}A \end{pmatrix} \begin{pmatrix} \mathbf{x}_1 \\ \mathbf{x}_2 \\ \vdots \\ \mathbf{x}_p \end{pmatrix}.
    $$

    第 $j$ 块为 $\sum_{k=1}^p b_{kj} A \mathbf{x}_k = \mathbf{y}_j$。两者一致。$\blacksquare$

    **方法二**（利用基矩阵）：设 $E_{ij}$ 为第 $(i,j)$ 位置为 1 其余为 0 的矩阵。注意 $\operatorname{vec}(AE_{ij}B)$ 的线性性，只需对 $X = E_{ij}$ 验证即可。$AE_{ij}B$ 的第 $(r,s)$ 元素为 $a_{ri}b_{js}$，这与 $(B^T \otimes A)$ 的相应元素一致。

!!! theorem "定理 19.5 (Vec 算子的特殊情形)"
    以下是定理 19.4 的几个常用特殊情形：

    1. $\operatorname{vec}(AX) = (I \otimes A)\operatorname{vec}(X)$（取 $B = I$）；
    2. $\operatorname{vec}(XB) = (B^T \otimes I)\operatorname{vec}(X)$（取 $A = I$）；
    3. $\operatorname{vec}(A\mathbf{x}\mathbf{b}^T) = (\mathbf{b} \otimes A)\mathbf{x}$（$X$ 为列向量时）；
    4. $\operatorname{vec}(\alpha A) = \alpha \operatorname{vec}(A)$。

??? proof "证明"
    均为定理 19.4 的直接推论，将相应的 $A$、$B$ 或 $X$ 取为特殊矩阵（单位阵、向量等）即可。例如：

    (1) 取 $B = I_p$，则 $\operatorname{vec}(AXI) = (I^T \otimes A)\operatorname{vec}(X) = (I \otimes A)\operatorname{vec}(X)$。

    (2) 取 $A = I_m$，则 $\operatorname{vec}(IXB) = (B^T \otimes I)\operatorname{vec}(X)$。$\blacksquare$

!!! definition "定义 19.4 (half-vectorization 算子)"
    对 $n \times n$ 对称矩阵 $A$，**半向量化算子** $\operatorname{vech}(A)$ 将 $A$ 的下三角部分（含对角线）按列堆叠为 $\frac{n(n+1)}{2} \times 1$ 向量：
    $$
    \operatorname{vech}(A) = (a_{11}, a_{21}, \ldots, a_{n1}, a_{22}, a_{32}, \ldots, a_{nn})^T.
    $$

!!! example "例 19.3"
    设 $X = \begin{pmatrix} x_{11} & x_{12} \\ x_{21} & x_{22} \end{pmatrix}$，$A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$。

    $\operatorname{vec}(X) = (x_{11}, x_{21}, x_{12}, x_{22})^T$。

    $AX = \begin{pmatrix} ax_{11}+bx_{21} & ax_{12}+bx_{22} \\ cx_{11}+dx_{21} & cx_{12}+dx_{22} \end{pmatrix}$。

    $\operatorname{vec}(AX) = (ax_{11}+bx_{21},\; cx_{11}+dx_{21},\; ax_{12}+bx_{22},\; cx_{12}+dx_{22})^T$。

    $(I_2 \otimes A)\operatorname{vec}(X) = \begin{pmatrix} A & 0 \\ 0 & A \end{pmatrix}\begin{pmatrix} x_{11} \\ x_{21} \\ x_{12} \\ x_{22} \end{pmatrix} = \begin{pmatrix} ax_{11}+bx_{21} \\ cx_{11}+dx_{21} \\ ax_{12}+bx_{22} \\ cx_{12}+dx_{22} \end{pmatrix}$。

    两者相等，验证了 $\operatorname{vec}(AX) = (I \otimes A)\operatorname{vec}(X)$。

!!! example "例 19.4"
    利用 Vec 公式求解 $AXB = C$。

    设 $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$，$B = \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix}$，$C = \begin{pmatrix} 6 & 2 \\ 0 & 8 \end{pmatrix}$。

    向量化：$(B^T \otimes A)\operatorname{vec}(X) = \operatorname{vec}(C)$。

    $B^T \otimes A = \begin{pmatrix} 3 & 0 \\ 0 & 6 \end{pmatrix} \otimes$... 更准确地：

    $B^T = \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix}$，故 $B^T \otimes A = \begin{pmatrix} 3A & 0 \\ 0 & A \end{pmatrix} = \begin{pmatrix} 3 & 0 & 0 & 0 \\ 0 & 6 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 2 \end{pmatrix}$。

    $\operatorname{vec}(C) = (6, 0, 2, 8)^T$。

    解：$\operatorname{vec}(X) = (B^T \otimes A)^{-1}\operatorname{vec}(C) = (2, 0, 2, 4)^T$。

    因此 $X = \begin{pmatrix} 2 & 2 \\ 0 & 4 \end{pmatrix}$。

    验证：$AXB = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}\begin{pmatrix} 2 & 2 \\ 0 & 4 \end{pmatrix}\begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 2 \\ 0 & 8 \end{pmatrix}\begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 6 & 2 \\ 0 & 8 \end{pmatrix} = C$。正确。

---

## 19.4 置换矩阵（Commutation Matrix）

<div class="context-flow" markdown>

**作用**：$K_{m,n}\operatorname{vec}(A)=\operatorname{vec}(A^T)$ · 实现 $A\otimes B \leftrightarrow B\otimes A$ 的指标重排 · $K_{n,n}^2=I$(对合)

</div>

Kronecker 积不满足交换律，但两种顺序的 Kronecker 积通过一个特殊的置换矩阵相联系。

!!! definition "定义 19.5 (置换矩阵 / 交换矩阵)"
    **置换矩阵**（commutation matrix）$K_{m,n}$ 是 $mn \times mn$ 的置换矩阵，满足对任意 $m \times n$ 矩阵 $A$：
    $$
    K_{m,n} \operatorname{vec}(A) = \operatorname{vec}(A^T).
    $$
    等价地，$K_{m,n}$ 可以用基矩阵表示为：
    $$
    K_{m,n} = \sum_{i=1}^{m} \sum_{j=1}^{n} E_{ij} \otimes E_{ji},
    $$
    其中 $E_{ij}$ 为 $m \times n$ 的基矩阵（第 $(i,j)$ 元素为 1，其余为 0），$E_{ji}$ 为 $n \times m$ 的基矩阵。

!!! theorem "定理 19.6 (置换矩阵的性质)"
    置换矩阵 $K_{m,n}$ 具有以下性质：

    1. $K_{m,n}^T = K_{m,n}^{-1} = K_{n,m}$。
    2. $K_{m,n}(A \otimes B)K_{p,q} = B \otimes A$（适当大小下）。特别地，$K_{m,n}(A \otimes B) = (B \otimes A)K_{p,q}$。
    3. $K_{n,n}$ 对称且正交，$K_{n,n}^2 = I_{n^2}$（对合矩阵）。
    4. $K_{1,n} = K_{n,1} = I_n$。
    5. $(A \otimes B) = K_{p,m}(B \otimes A)K_{n,q}$，其中 $A$ 为 $m \times n$，$B$ 为 $p \times q$。

??? proof "证明"
    **(1)**：对任意 $m \times n$ 矩阵 $A$，$K_{m,n}\operatorname{vec}(A) = \operatorname{vec}(A^T)$。将 $A$ 替换为 $A^T$（$n \times m$）：$K_{n,m}\operatorname{vec}(A^T) = \operatorname{vec}(A)$。因此 $K_{n,m}K_{m,n}\operatorname{vec}(A) = \operatorname{vec}(A)$ 对所有 $A$ 成立，故 $K_{n,m}K_{m,n} = I_{mn}$，即 $K_{m,n}^{-1} = K_{n,m}$。

    又 $K_{m,n}$ 为置换矩阵，故 $K_{m,n}^T = K_{m,n}^{-1} = K_{n,m}$。

    **(2)**：对任意 $n \times q$ 矩阵 $X$，取 $p \times m$ 矩阵的情形：
    $$
    K_{m,n}(A \otimes B)\operatorname{vec}(X) = K_{m,n}\operatorname{vec}(BXA^T) \quad \text{(需要调整大小)}.
    $$

    更严格地，对 $A$ 为 $m \times n$，$B$ 为 $p \times q$：考虑 $(A \otimes B)$ 作用于 $\operatorname{vec}(X)$（$X$ 为 $q \times 1$ 向量不合适）。

    采用直接验证：由定义 $K_{m,n}$ 将 $\operatorname{vec}$ 中按列排列的元素重排为按行排列，$(A \otimes B)$ 的第 $((i-1)p+k, (j-1)q+l)$ 元素为 $a_{ij}b_{kl}$，而 $(B \otimes A)$ 的第 $((k-1)m+i, (l-1)n+j)$ 元素也是 $b_{kl}a_{ij}$。置换矩阵恰好实现了指标 $((i-1)p+k) \leftrightarrow ((k-1)m+i)$ 的重排。$\blacksquare$

!!! theorem "定理 19.7 (Vec 与转置的关系)"
    对任意 $m \times n$ 矩阵 $A$：
    $$
    \operatorname{vec}(A^T) = K_{m,n}\operatorname{vec}(A).
    $$
    进而，对矩阵乘积的转置：
    $$
    \operatorname{vec}((AXB)^T) = \operatorname{vec}(B^T X^T A^T) = (A \otimes B^T)\operatorname{vec}(X^T) = (A \otimes B^T)K_{n,p}\operatorname{vec}(X).
    $$

??? proof "证明"
    第一个等式即为 $K_{m,n}$ 的定义。第二个等式利用了 $\operatorname{vec}(B^T X^T A^T) = (A \otimes B^T)\operatorname{vec}(X^T)$（由定理 19.4），再用 $\operatorname{vec}(X^T) = K_{n,p}\operatorname{vec}(X)$。$\blacksquare$

!!! example "例 19.5"
    构造 $K_{2,3}$。需要找到 $6 \times 6$ 置换矩阵，使得对任意 $2 \times 3$ 矩阵 $A = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \end{pmatrix}$：

    $\operatorname{vec}(A) = (a_{11}, a_{21}, a_{12}, a_{22}, a_{13}, a_{23})^T$。

    $\operatorname{vec}(A^T) = (a_{11}, a_{12}, a_{13}, a_{21}, a_{22}, a_{23})^T$。

    因此 $K_{2,3}$ 将位置 $(1,2,3,4,5,6)$ 映射到 $(1,3,5,2,4,6)$：
    $$
    K_{2,3} = \begin{pmatrix}
    1 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 1 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 1 & 0 \\
    0 & 1 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 1 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 1
    \end{pmatrix}.
    $$

---

## 19.5 Kronecker 积与矩阵方程

<div class="context-flow" markdown>

**应用**：$\sum A_kXB_k=C$ → $(\sum B_k^T\otimes A_k)\operatorname{vec}(X)=\operatorname{vec}(C)$ · **Sylvester** $AX+XB=C$：可解 $\Leftrightarrow$ $\sigma(A)\cap\sigma(-B)=\emptyset$ → Ch20

</div>

Kronecker 积和 Vec 算子的一个最重要的应用是将矩阵方程转化为标准的线性方程组。

!!! definition "定义 19.6 (线性矩阵方程)"
    形如
    $$
    \sum_{k=1}^{K} A_k X B_k = C
    $$
    的方程称为**线性矩阵方程**（linear matrix equation），其中 $A_k, B_k, C$ 已知，$X$ 为未知矩阵。

!!! theorem "定理 19.8 (矩阵方程的向量化)"
    线性矩阵方程 $\sum_{k=1}^K A_k X B_k = C$ 等价于向量方程：
    $$
    \left(\sum_{k=1}^{K} B_k^T \otimes A_k\right) \operatorname{vec}(X) = \operatorname{vec}(C).
    $$

??? proof "证明"
    对每一项 $A_k X B_k$ 应用定理 19.4：
    $$
    \operatorname{vec}(A_k X B_k) = (B_k^T \otimes A_k)\operatorname{vec}(X).
    $$

    对方程两边取 Vec：
    $$
    \operatorname{vec}\left(\sum_{k=1}^K A_k X B_k\right) = \sum_{k=1}^K \operatorname{vec}(A_k X B_k) = \sum_{k=1}^K (B_k^T \otimes A_k)\operatorname{vec}(X) = \left(\sum_{k=1}^K B_k^T \otimes A_k\right)\operatorname{vec}(X).
    $$

    右边 $\operatorname{vec}(C) = \operatorname{vec}(C)$。$\blacksquare$

!!! theorem "定理 19.9 (Sylvester 方程的 Kronecker 积形式)"
    Sylvester 方程 $AX + XB = C$（其中 $A$ 为 $m \times m$，$B$ 为 $n \times n$，$X, C$ 为 $m \times n$）等价于：
    $$
    (I_n \otimes A + B^T \otimes I_m)\operatorname{vec}(X) = \operatorname{vec}(C).
    $$
    该方程有唯一解当且仅当 $I_n \otimes A + B^T \otimes I_m$ 非奇异，即 $A$ 与 $-B$ 无公共特征值。

??? proof "证明"
    由定理 19.8，取 $K = 2$，$A_1 = A$，$B_1 = I$，$A_2 = I$，$B_2 = B$：
    $$
    (I^T \otimes A + B^T \otimes I)\operatorname{vec}(X) = (I_n \otimes A + B^T \otimes I_m)\operatorname{vec}(X) = \operatorname{vec}(C).
    $$

    $I_n \otimes A + B^T \otimes I_m$ 的特征值为 $\{\lambda_i(A) + \lambda_j(B) : i = 1,\ldots,m;\; j=1,\ldots,n\}$（见定理 19.12），因此它非奇异当且仅当 $\lambda_i(A) + \lambda_j(B) \neq 0$ 对所有 $i, j$ 成立，即 $A$ 与 $-B$ 无公共特征值。$\blacksquare$

!!! example "例 19.6"
    求解 Sylvester 方程 $AX + XB = C$，其中：
    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}, \quad B = (3), \quad C = \begin{pmatrix} 4 \\ 10 \end{pmatrix}.
    $$
    这里 $A$ 为 $2 \times 2$，$B$ 为 $1 \times 1$（标量 3），$X$ 为 $2 \times 1$。

    向量化：$(I_1 \otimes A + B^T \otimes I_2)\operatorname{vec}(X) = \operatorname{vec}(C)$。

    $I_1 \otimes A = A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$，$B^T \otimes I_2 = 3I_2 = \begin{pmatrix} 3 & 0 \\ 0 & 3 \end{pmatrix}$。

    系数矩阵：$\begin{pmatrix} 4 & 0 \\ 0 & 5 \end{pmatrix}$。

    $\operatorname{vec}(X) = \begin{pmatrix} 4 & 0 \\ 0 & 5 \end{pmatrix}^{-1}\begin{pmatrix} 4 \\ 10 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$。

    因此 $X = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$。

    验证：$AX + XB = \begin{pmatrix} 1 \\ 4 \end{pmatrix} + \begin{pmatrix} 3 \\ 6 \end{pmatrix} = \begin{pmatrix} 4 \\ 10 \end{pmatrix} = C$。正确。

!!! example "例 19.7"
    利用 Kronecker 积判断矩阵方程的可解性。

    考虑 $AX - XA = C$（其中 $A$ 为 $n \times n$）。向量化为 $(I \otimes A - A^T \otimes I)\operatorname{vec}(X) = \operatorname{vec}(C)$。

    若 $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$，则系数矩阵的特征值为 $\{\lambda_i - \lambda_j\} = \{0, -1, 1, 0\}$。

    由于存在零特征值（$\lambda_1 - \lambda_1 = 0$ 和 $\lambda_2 - \lambda_2 = 0$），系数矩阵奇异，因此方程 $AX - XA = C$ 不对所有 $C$ 有解。

    具体地，$C$ 必须满足 $\operatorname{tr}(C) = 0$（因为 $\operatorname{tr}(AX - XA) = 0$），且对角元素为 0 也是必要条件（对于对角 $A$）。

---

## 19.6 Kronecker 积的特征值分解

<div class="context-flow" markdown>

**洞察**：$\sigma(A\otimes B)=\{\lambda_i\mu_j\}$，特征向量 $\mathbf{u}_i\otimes\mathbf{v}_j$ · SVD也乘积化：$\sigma_k(A\otimes B)=\{\sigma_i(A)\sigma_j(B)\}$

</div>

!!! definition "定义 19.7 (Kronecker 积的谱)"
    设 $A$ 为 $m \times m$ 矩阵，特征值 $\lambda_1, \ldots, \lambda_m$；$B$ 为 $n \times n$ 矩阵，特征值 $\mu_1, \ldots, \mu_n$。则 $A \otimes B$ 的 $mn$ 个特征值为：
    $$
    \sigma(A \otimes B) = \{\lambda_i \mu_j : i = 1,\ldots,m;\; j = 1,\ldots,n\}.
    $$

!!! theorem "定理 19.10 (Kronecker 积的特征值与特征向量)"
    设 $A\mathbf{u} = \lambda\mathbf{u}$，$B\mathbf{v} = \mu\mathbf{v}$，则：
    $$
    (A \otimes B)(\mathbf{u} \otimes \mathbf{v}) = \lambda\mu(\mathbf{u} \otimes \mathbf{v}).
    $$
    即 $\mathbf{u} \otimes \mathbf{v}$ 是 $A \otimes B$ 对应特征值 $\lambda\mu$ 的特征向量。

    若 $A$ 和 $B$ 均可对角化，$A = P \operatorname{diag}(\lambda_1,\ldots,\lambda_m) P^{-1}$，$B = Q \operatorname{diag}(\mu_1,\ldots,\mu_n) Q^{-1}$，则：
    $$
    A \otimes B = (P \otimes Q) \operatorname{diag}(\lambda_1\mu_1, \lambda_1\mu_2, \ldots, \lambda_m\mu_n) (P \otimes Q)^{-1}.
    $$

??? proof "证明"
    由混合积性质：
    $$
    (A \otimes B)(\mathbf{u} \otimes \mathbf{v}) = (A\mathbf{u}) \otimes (B\mathbf{v}) = (\lambda\mathbf{u}) \otimes (\mu\mathbf{v}) = \lambda\mu(\mathbf{u} \otimes \mathbf{v}).
    $$

    对于对角化的情形：
    $$
    A \otimes B = (P \Lambda_A P^{-1}) \otimes (Q \Lambda_B Q^{-1}) = (P \otimes Q)(\Lambda_A \otimes \Lambda_B)(P^{-1} \otimes Q^{-1}).
    $$

    由 $(P \otimes Q)^{-1} = P^{-1} \otimes Q^{-1}$，以及 $\Lambda_A \otimes \Lambda_B$ 为对角矩阵（对角元素为 $\lambda_i \mu_j$），即得结论。$\blacksquare$

!!! theorem "定理 19.11 (Kronecker 积的奇异值分解)"
    设 $A = U_A \Sigma_A V_A^*$，$B = U_B \Sigma_B V_B^*$ 为 $A, B$ 的奇异值分解，则：
    $$
    A \otimes B = (U_A \otimes U_B)(\Sigma_A \otimes \Sigma_B)(V_A \otimes V_B)^*.
    $$
    因此 $A \otimes B$ 的奇异值为 $\{\sigma_i(A)\sigma_j(B)\}$。

??? proof "证明"
    由混合积性质：
    $$
    A \otimes B = (U_A \Sigma_A V_A^*) \otimes (U_B \Sigma_B V_B^*) = (U_A \otimes U_B)(\Sigma_A \otimes \Sigma_B)(V_A^* \otimes V_B^*).
    $$

    由 $(V_A^* \otimes V_B^*) = (V_A \otimes V_B)^*$，且 $U_A \otimes U_B$ 和 $V_A \otimes V_B$ 均为酉矩阵（因为酉矩阵的 Kronecker 积仍为酉矩阵），$\Sigma_A \otimes \Sigma_B$ 为非负对角矩阵，这正是奇异值分解。$\blacksquare$

!!! example "例 19.8"
    设 $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix}$。

    $A$ 的特征值：$\lambda_1 = 2$，$\lambda_2 = 3$。$B$ 的特征值：$\mu_1 = 1$，$\mu_2 = 4$。

    $A \otimes B$ 的特征值：$\{2 \times 1, 2 \times 4, 3 \times 1, 3 \times 4\} = \{2, 8, 3, 12\}$。

    直接计算 $A \otimes B = \begin{pmatrix} 2 & 0 & 0 & 0 \\ 0 & 8 & 0 & 0 \\ 0 & 0 & 3 & 0 \\ 0 & 0 & 0 & 12 \end{pmatrix}$，对角矩阵，特征值为 $2, 8, 3, 12$。验证成功。

---

## 19.7 Kronecker 和

<div class="context-flow" markdown>

**脉络**：$A\oplus B = A\otimes I+I\otimes B$，特征值 $\lambda_i+\mu_j$ · **关键等式**：$e^{A\oplus B}=e^A\otimes e^B$(因 $A\otimes I$ 与 $I\otimes B$ 可交换) → Ch20 Lyapunov方程

</div>

Kronecker 和是与 Kronecker 积密切相关的另一种运算，它与矩阵指数和 Lyapunov 方程有着深刻联系。

!!! definition "定义 19.8 (Kronecker 和)"
    设 $A$ 为 $m \times m$ 矩阵，$B$ 为 $n \times n$ 矩阵。$A$ 与 $B$ 的 **Kronecker 和**（Kronecker sum）定义为：
    $$
    A \oplus B = A \otimes I_n + I_m \otimes B.
    $$
    它是 $mn \times mn$ 矩阵。

!!! theorem "定理 19.12 (Kronecker 和的特征值)"
    设 $A$ 的特征值为 $\lambda_1, \ldots, \lambda_m$，$B$ 的特征值为 $\mu_1, \ldots, \mu_n$。则 $A \oplus B$ 的 $mn$ 个特征值为：
    $$
    \sigma(A \oplus B) = \{\lambda_i + \mu_j : i = 1,\ldots,m;\; j = 1,\ldots,n\}.
    $$
    对应的特征向量为 $\mathbf{u}_i \otimes \mathbf{v}_j$。

??? proof "证明"
    设 $A\mathbf{u}_i = \lambda_i \mathbf{u}_i$，$B\mathbf{v}_j = \mu_j \mathbf{v}_j$。则：
    $$
    (A \oplus B)(\mathbf{u}_i \otimes \mathbf{v}_j) = (A \otimes I_n + I_m \otimes B)(\mathbf{u}_i \otimes \mathbf{v}_j)
    $$
    $$
    = (A\mathbf{u}_i) \otimes (I_n\mathbf{v}_j) + (I_m\mathbf{u}_i) \otimes (B\mathbf{v}_j) = \lambda_i(\mathbf{u}_i \otimes \mathbf{v}_j) + \mu_j(\mathbf{u}_i \otimes \mathbf{v}_j) = (\lambda_i + \mu_j)(\mathbf{u}_i \otimes \mathbf{v}_j).
    $$

    当 $A, B$ 可对角化时，这给出了所有 $mn$ 个特征值。一般情况需要利用 Jordan 标准形，但结论相同（计入代数重数）。$\blacksquare$

!!! theorem "定理 19.13 (Kronecker 和与矩阵指数)"
    设 $A$ 为 $m \times m$ 矩阵，$B$ 为 $n \times n$ 矩阵。则：
    $$
    e^{A \oplus B} = e^A \otimes e^B.
    $$

??? proof "证明"
    关键观察是 $A \otimes I$ 和 $I \otimes B$ **可交换**：
    $$
    (A \otimes I)(I \otimes B) = A \otimes B = (I \otimes B)(A \otimes I).
    $$

    由于 $A \oplus B = A \otimes I + I \otimes B$，且这两个矩阵可交换，因此矩阵指数满足：
    $$
    e^{A \oplus B} = e^{A \otimes I + I \otimes B} = e^{A \otimes I} \cdot e^{I \otimes B}.
    $$

    又 $e^{A \otimes I} = \sum_{k=0}^{\infty} \frac{(A \otimes I)^k}{k!} = \sum_{k=0}^{\infty} \frac{A^k \otimes I}{k!} = \left(\sum_{k=0}^{\infty}\frac{A^k}{k!}\right) \otimes I = e^A \otimes I$。

    类似地 $e^{I \otimes B} = I \otimes e^B$。

    因此 $e^{A \oplus B} = (e^A \otimes I)(I \otimes e^B) = e^A \otimes e^B$。$\blacksquare$

!!! theorem "定理 19.14 (Kronecker 和与 Lyapunov 方程)"
    Lyapunov 方程 $AX + XA^T = C$ 等价于：
    $$
    (A \oplus A^T)\operatorname{vec}(X) = (I \otimes A + A^* \otimes I)\operatorname{vec}(X) = \operatorname{vec}(C),
    $$
    其中我们注意到 $(A^T)^T = A$，故 $I \otimes A + (A^T)^T \otimes I = I \otimes A + A \otimes I$。

    更准确地写：向量化后为 $(I_n \otimes A + \bar{A} \otimes I_n)\operatorname{vec}(X) = \operatorname{vec}(C)$（实数情形下 $\bar{A} = A$）。

    $A \oplus A^T$ 的特征值为 $\lambda_i(A) + \lambda_j(A^T) = \lambda_i(A) + \lambda_j(A)$，方程有唯一解当且仅当 $\lambda_i(A) + \lambda_j(A) \neq 0$ 对所有 $i,j$ 成立。

??? proof "证明"
    对 $AX + XA^T = C$ 两边取 Vec：
    $$
    \operatorname{vec}(AX) + \operatorname{vec}(XA^T) = \operatorname{vec}(C).
    $$

    由定理 19.5：$\operatorname{vec}(AX) = (I \otimes A)\operatorname{vec}(X)$，$\operatorname{vec}(XA^T) = ((A^T)^T \otimes I)\operatorname{vec}(X) = (A \otimes I)\operatorname{vec}(X)$。

    因此 $(I \otimes A + A \otimes I)\operatorname{vec}(X) = (A \oplus A)\operatorname{vec}(X) = \operatorname{vec}(C)$。

    注意这里 $A \oplus A = A \otimes I + I \otimes A$（在实数情形下，$A^T$ 的转置回来就是 $A$）。特征值的结论由定理 19.12 直接得出。$\blacksquare$

!!! example "例 19.9"
    设 $A = \begin{pmatrix} -1 & 0 \\ 0 & -2 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$。

    $A \oplus B = A \otimes I_2 + I_2 \otimes B = \begin{pmatrix} -1 & 0 & 0 & 0 \\ 0 & -1 & 0 & 0 \\ 0 & 0 & -2 & 0 \\ 0 & 0 & 0 & -2 \end{pmatrix} + \begin{pmatrix} 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & -1 & 0 \end{pmatrix} = \begin{pmatrix} -1 & 1 & 0 & 0 \\ -1 & -1 & 0 & 0 \\ 0 & 0 & -2 & 1 \\ 0 & 0 & -1 & -2 \end{pmatrix}$。

    $A$ 的特征值：$-1, -2$。$B$ 的特征值：$i, -i$。

    $A \oplus B$ 的特征值：$\{-1+i, -1-i, -2+i, -2-i\}$。

    验证：$\begin{pmatrix} -1 & 1 \\ -1 & -1 \end{pmatrix}$ 的特征值为 $-1 \pm i$，$\begin{pmatrix} -2 & 1 \\ -1 & -2 \end{pmatrix}$ 的特征值为 $-2 \pm i$。正确。

!!! example "例 19.10"
    验证 $e^{A \oplus B} = e^A \otimes e^B$。

    取 $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$。

    $e^A = \begin{pmatrix} e & 0 \\ 0 & 1 \end{pmatrix}$，$e^B = \begin{pmatrix} 1 & 0 \\ 0 & e \end{pmatrix}$。

    $e^A \otimes e^B = \begin{pmatrix} e & 0 & 0 & 0 \\ 0 & e^2 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & e \end{pmatrix}$。

    $A \oplus B = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix} + \begin{pmatrix} 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 2 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$。

    $e^{A \oplus B} = \begin{pmatrix} e & 0 & 0 & 0 \\ 0 & e^2 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & e \end{pmatrix} = e^A \otimes e^B$。验证成功。
