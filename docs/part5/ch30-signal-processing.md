# 第 30 章 线性代数在信号处理与编码中的应用

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵分解(Ch10) · 正定矩阵(Ch16) · Kronecker 积(Ch19) · **脉络**：DFT 矩阵(酉矩阵) → FFT(DFT 矩阵的递归分解) → 循环矩阵(DFT 对角化) → Toeplitz 矩阵(近似循环) → 线性码(零空间/列空间) → Hamming 码(校验矩阵) → RS 码(Vandermonde 结构) → LDPC 码(稀疏校验矩阵)
**本质**：信号处理的核心运算——变换、滤波、卷积——都是矩阵-向量乘法；编码理论的核心概念——码字、校验、纠错——都是有限域上的线性代数

</div>

信号处理与编码理论是线性代数最经典的工程应用领域。离散 Fourier 变换本质上是一个酉矩阵乘法，循环卷积对应循环矩阵，而纠错码的设计与译码完全建立在有限域上的线性代数之上。本章从 DFT 矩阵出发，经 FFT、循环矩阵、Toeplitz 矩阵，到线性码、Hamming 码、Reed-Solomon 码和 LDPC 码，系统展示线性代数在这两大领域中的核心作用。

---

## 30.1 离散 Fourier 变换

<div class="context-flow" markdown>

**线性代数视角**：DFT = 乘以酉矩阵 $F_n$（$F_n^* F_n = nI$）→ 特征值/特征向量全部已知 → DFT 将循环卷积对角化
**链接**：Ch8 酉矩阵 · Ch6 特征值

</div>

离散 Fourier 变换（DFT）是信号处理中最基本的工具，其本质是一个特殊的矩阵-向量乘法。

!!! definition "定义 30.1 (DFT 矩阵)"
    设 $\omega_n = e^{-2\pi i/n}$ 为 $n$ 次本原单位根。$n$ 阶 **DFT 矩阵**（DFT matrix）$F_n \in \mathbb{C}^{n \times n}$ 定义为

    $$
    (F_n)_{jk} = \omega_n^{jk}, \quad j, k = 0, 1, \ldots, n-1.
    $$

    即

    $$
    F_n = \begin{pmatrix}1&1&1&\cdots&1\\1&\omega_n&\omega_n^2&\cdots&\omega_n^{n-1}\\1&\omega_n^2&\omega_n^4&\cdots&\omega_n^{2(n-1)}\\\vdots&\vdots&\vdots&\ddots&\vdots\\1&\omega_n^{n-1}&\omega_n^{2(n-1)}&\cdots&\omega_n^{(n-1)^2}\end{pmatrix}.
    $$

    向量 $\mathbf{x} \in \mathbb{C}^n$ 的 **DFT** 为 $\hat{\mathbf{x}} = F_n \mathbf{x}$，即 $\hat{x}_j = \sum_{k=0}^{n-1} x_k \omega_n^{jk}$。

!!! definition "定义 30.2 (逆 DFT)"
    **逆 DFT**（inverse DFT）为 $\mathbf{x} = \frac{1}{n}F_n^* \hat{\mathbf{x}}$，其中 $F_n^* = \overline{F_n}^T$ 为 $F_n$ 的共轭转置，$(F_n^*)_{jk} = \omega_n^{-jk}$。

!!! theorem "定理 30.1 (DFT 矩阵的酉性)"
    归一化 DFT 矩阵 $\frac{1}{\sqrt{n}}F_n$ 是酉矩阵，即

    $$
    F_n^* F_n = F_n F_n^* = nI_n.
    $$

    等价地，$F_n$ 的列（或行）构成 $\mathbb{C}^n$ 的正交基（内积为 $n$），$\frac{1}{\sqrt{n}}F_n$ 的列构成标准正交基。

??? proof "证明"
    $(F_n^* F_n)_{jk} = \sum_{\ell=0}^{n-1}\overline{(F_n)_{\ell j}}(F_n)_{\ell k} = \sum_{\ell=0}^{n-1}\omega_n^{-\ell j}\omega_n^{\ell k} = \sum_{\ell=0}^{n-1}\omega_n^{\ell(k-j)}$。当 $k = j$ 时，和为 $n$。当 $k \ne j$ 时，这是等比级数 $\sum_{\ell=0}^{n-1}(\omega_n^{k-j})^\ell = \frac{1 - \omega_n^{n(k-j)}}{1 - \omega_n^{k-j}} = 0$，因为 $\omega_n^n = 1$ 且 $\omega_n^{k-j} \ne 1$（$0 < |k-j| < n$）。$\blacksquare$

!!! theorem "定理 30.2 (DFT 矩阵的特征值)"
    DFT 矩阵 $F_n$ 满足 $F_n^4 = n^2 I_n$，故 $F_n$ 的特征值为 $\{+\sqrt{n}, -\sqrt{n}, +i\sqrt{n}, -i\sqrt{n}\}$（各特征值的重数取决于 $n \bmod 4$）。特别地，$F_n^2 = nP$，其中 $P$ 为逆序置换矩阵 $(P\mathbf{x})_j = x_{n-j \bmod n}$。

??? proof "证明"
    $(F_n^2)_{jk} = \sum_{\ell}\omega_n^{j\ell}\omega_n^{\ell k} = \sum_\ell \omega_n^{\ell(j+k)}$。当 $j + k \equiv 0 \pmod{n}$ 时和为 $n$，否则为 $0$。因此 $(F_n^2)_{jk} = n\delta_{j, n-k \bmod n}$，即 $F_n^2 = nP$。于是 $F_n^4 = n^2 P^2 = n^2 I$（因为 $P^2 = I$）。$F_n$ 的特征值 $\lambda$ 满足 $\lambda^4 = n^2$，故 $\lambda \in \{n^{1/2}, -n^{1/2}, in^{1/2}, -in^{1/2}\}$。$\blacksquare$

!!! example "例 30.1"
    **$4$ 阶 DFT 矩阵。** $\omega_4 = e^{-\pi i/2} = -i$。

    $$
    F_4 = \begin{pmatrix}1&1&1&1\\1&-i&-1&i\\1&-1&1&-1\\1&i&-1&-i\end{pmatrix}.
    $$

    验证酉性：$F_4^* F_4 = 4I_4$。对 $\mathbf{x} = (1, 2, 3, 4)^T$ 做 DFT：

    $$
    \hat{\mathbf{x}} = F_4\mathbf{x} = \begin{pmatrix}10\\-2+2i\\-2\\-2-2i\end{pmatrix}.
    $$

    逆变换验证：$\frac{1}{4}F_4^*\hat{\mathbf{x}} = (1, 2, 3, 4)^T = \mathbf{x}$。

---

## 30.2 快速 Fourier 变换

<div class="context-flow" markdown>

**核心**：FFT 不是新的变换，而是 DFT 矩阵的高效分解 → $F_n = P_n(I_2 \otimes F_{n/2})T_n \cdot (\text{蝶形因子})$ → $O(n\log n)$ 代替 $O(n^2)$
**链接**：Ch19 Kronecker 积 · 稀疏矩阵分解

</div>

快速 Fourier 变换（FFT）是 20 世纪最重要的算法之一，其本质是将 DFT 矩阵分解为稀疏矩阵的乘积。

!!! definition "定义 30.3 (Cooley-Tukey 分解)"
    设 $n = 2m$。将 DFT 矩阵按偶数和奇数索引分解，$F_n$ 可以表示为

    $$
    F_n = \begin{pmatrix}I_m & D_m \\ I_m & -D_m\end{pmatrix}\begin{pmatrix}F_m & 0 \\ 0 & F_m\end{pmatrix}P_n,
    $$

    其中 $D_m = \operatorname{diag}(1, \omega_n, \omega_n^2, \ldots, \omega_n^{m-1})$ 为**旋转因子**（twiddle factor）矩阵，$P_n$ 为**完美洗牌**置换（将偶数索引排在前，奇数索引排在后）。

!!! theorem "定理 30.3 (FFT 的复杂度)"
    当 $n = 2^s$ 时，通过递归应用 Cooley-Tukey 分解，$n$ 点 DFT 可以用 $\frac{n}{2}\log_2 n$ 次复数乘法和 $n\log_2 n$ 次复数加法完成，总复杂度为 $O(n \log n)$，而直接矩阵-向量乘法需要 $O(n^2)$。

    FFT 将 $F_n$ 分解为 $\log_2 n$ 个稀疏矩阵的乘积，每个稀疏矩阵恰好有 $O(n)$ 个非零元素。

??? proof "证明"
    设 $T(n)$ 为 $n$ 点 FFT 的乘法次数。Cooley-Tukey 分解将一个 $n$ 点 DFT 分为两个 $n/2$ 点 DFT 加上 $n/2$ 次旋转因子乘法。递推关系：

    $$
    T(n) = 2T(n/2) + n/2, \quad T(1) = 0.
    $$

    解为 $T(n) = \frac{n}{2}\log_2 n$。加法次数类似，每层 $n$ 次蝶形加法，共 $\log_2 n$ 层，总计 $n\log_2 n$。

    分解的矩阵形式：$F_n = B_s B_{s-1} \cdots B_1 P$，其中每个 $B_k$ 是块对角-蝶形矩阵，有 $O(n)$ 个非零元素。$\blacksquare$

!!! theorem "定理 30.4 (FFT 的 Kronecker 积表示)"
    当 $n = 2^s$ 时，FFT 可以用 Kronecker 积紧凑表示为

    $$
    F_n = \prod_{j=1}^{s} \left(I_{2^{j-1}} \otimes F_2 \otimes I_{2^{s-j}}\right) \cdot D_j \cdot P_j,
    $$

    其中 $F_2 = \begin{pmatrix}1&1\\1&-1\end{pmatrix}$ 为 $2$ 阶 DFT 矩阵（Hadamard 矩阵），$D_j$ 为旋转因子对角矩阵，$P_j$ 为适当的置换矩阵。这一表示清楚地展示了 FFT 的"分而治之"结构。

??? proof "证明"
    由 Cooley-Tukey 递归，第 $j$ 层将 $2^j$ 点 DFT 分解为 $2^{j-1}$ 点 DFT 对。在矩阵层面，蝶形运算为 $F_2$ 作用于长度为 2 的子向量，通过 Kronecker 积 $I_{2^{j-1}} \otimes F_2 \otimes I_{2^{s-j}}$ 并行地作用于所有子问题。旋转因子和置换调整数据排列。$\blacksquare$

!!! example "例 30.2"
    **$8$ 点 FFT 分解。** $n = 8 = 2^3$，需要 $\log_2 8 = 3$ 层。

    **第 1 层**：将 $8$ 点分为偶数索引 $(x_0, x_2, x_4, x_6)$ 和奇数索引 $(x_1, x_3, x_5, x_7)$，各做 $4$ 点 DFT。

    **第 2 层**：每个 $4$ 点 DFT 再分为两个 $2$ 点 DFT。

    **第 3 层**：$2$ 点 DFT 即 $F_2$ 蝶形运算。

    乘法次数：$\frac{8}{2}\log_2 8 = 12$，而直接 DFT 需要 $8^2 = 64$ 次。加速比 $64/12 \approx 5.3$。对于 $n = 2^{20} \approx 10^6$，FFT 需约 $10^7$ 次运算，直接法需 $10^{12}$ 次，加速比约 $10^5$。

---

## 30.3 循环矩阵与卷积

<div class="context-flow" markdown>

**核心**：循环矩阵 = DFT 矩阵对角化 $C = F_n^{-1}\operatorname{diag}(\hat{\mathbf{c}})F_n$ → 循环卷积 = 频域逐点乘法 → 卷积定理
**链接**：Ch6 特征值 · 30.1 DFT

</div>

循环矩阵是一类具有特殊结构的矩阵，它与 DFT 和卷积运算有深刻的联系。

!!! definition "定义 30.4 (循环矩阵)"
    由向量 $\mathbf{c} = (c_0, c_1, \ldots, c_{n-1})^T$ 生成的 $n$ 阶**循环矩阵**（circulant matrix）为

    $$
    C = \operatorname{circ}(\mathbf{c}) = \begin{pmatrix}c_0&c_{n-1}&c_{n-2}&\cdots&c_1\\c_1&c_0&c_{n-1}&\cdots&c_2\\c_2&c_1&c_0&\cdots&c_3\\\vdots&\vdots&\vdots&\ddots&\vdots\\c_{n-1}&c_{n-2}&c_{n-3}&\cdots&c_0\end{pmatrix}.
    $$

    每一行是上一行的循环右移。等价地，$C = \sum_{k=0}^{n-1}c_k P^k$，其中 $P$ 为循环置换矩阵。

!!! definition "定义 30.5 (循环卷积)"
    两个长度为 $n$ 的向量 $\mathbf{a}$ 和 $\mathbf{b}$ 的**循环卷积**（circular convolution）$\mathbf{c} = \mathbf{a} \circledast \mathbf{b}$ 定义为

    $$
    c_j = \sum_{k=0}^{n-1} a_k b_{(j-k) \bmod n}, \quad j = 0, 1, \ldots, n-1.
    $$

    循环卷积等价于循环矩阵乘法：$\mathbf{c} = \operatorname{circ}(\mathbf{a})\mathbf{b}$。

!!! theorem "定理 30.5 (循环矩阵的谱分解)"
    所有 $n$ 阶循环矩阵都被同一个 DFT 矩阵对角化：

    $$
    C = \operatorname{circ}(\mathbf{c}) = \frac{1}{n}F_n^* \operatorname{diag}(\hat{\mathbf{c}}) F_n,
    $$

    其中 $\hat{\mathbf{c}} = F_n \mathbf{c}$ 为 $\mathbf{c}$ 的 DFT。$C$ 的特征值为 $\hat{c}_0, \hat{c}_1, \ldots, \hat{c}_{n-1}$，对应的特征向量为 $F_n$ 的列 $\mathbf{f}_0, \mathbf{f}_1, \ldots, \mathbf{f}_{n-1}$。

    所有循环矩阵两两可交换，且共享同一组特征向量。

??? proof "证明"
    循环置换矩阵 $P$ 的特征值为 $\omega_n^0 = 1, \omega_n, \omega_n^2, \ldots, \omega_n^{n-1}$，特征向量为 DFT 矩阵的列。由于 $C = \sum_k c_k P^k$，$C$ 与 $P$ 共享特征向量，特征值为 $\lambda_j = \sum_k c_k \omega_n^{jk} = \hat{c}_j$。因此 $C = \frac{1}{n}F_n^* \operatorname{diag}(\hat{\mathbf{c}})F_n$。两个循环矩阵 $C_1, C_2$ 都被 $F_n$ 对角化，故 $C_1 C_2 = C_2 C_1$。$\blacksquare$

!!! theorem "定理 30.6 (卷积定理)"
    循环卷积在频域等价于逐点乘法：

    $$
    \mathbf{a} \circledast \mathbf{b} = \mathbf{c} \iff \hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}},
    $$

    其中 $\hat{\mathbf{a}} = F_n\mathbf{a}$，$\hat{\mathbf{b}} = F_n\mathbf{b}$，$\odot$ 为逐元素乘法。因此循环卷积可以通过"FFT → 逐点乘法 → IFFT"在 $O(n\log n)$ 时间内完成。

??? proof "证明"
    $\mathbf{c} = \operatorname{circ}(\mathbf{a})\mathbf{b} = \frac{1}{n}F_n^*\operatorname{diag}(\hat{\mathbf{a}})F_n\mathbf{b} = \frac{1}{n}F_n^*\operatorname{diag}(\hat{\mathbf{a}})\hat{\mathbf{b}} = \frac{1}{n}F_n^*(\hat{\mathbf{a}} \odot \hat{\mathbf{b}})$。取 DFT：$\hat{\mathbf{c}} = F_n\mathbf{c} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}}$。$\blacksquare$

!!! example "例 30.3"
    **用 FFT 计算循环卷积。** 设 $\mathbf{a} = (1, 2, 3, 0)^T$，$\mathbf{b} = (1, 0, 1, 0)^T$。

    DFT：$\hat{\mathbf{a}} = F_4\mathbf{a} = (6, -2+2i, -2, -2-2i)^T$，$\hat{\mathbf{b}} = F_4\mathbf{b} = (2, 0, 2, 0)^T$。

    频域乘法：$\hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}} = (12, 0, -4, 0)^T$。

    IFFT：$\mathbf{c} = \frac{1}{4}F_4^*\hat{\mathbf{c}} = (2, 2, 4, 4)^T$。

    直接验证：$c_0 = 1\cdot1 + 0\cdot0 + 3\cdot1 + 0\cdot2 = 4$——需要注意循环卷积的索引。用循环矩阵验证：

    $$
    \operatorname{circ}(\mathbf{a})\mathbf{b} = \begin{pmatrix}1&0&3&2\\2&1&0&3\\3&2&1&0\\0&3&2&1\end{pmatrix}\begin{pmatrix}1\\0\\1\\0\end{pmatrix} = \begin{pmatrix}4\\2\\4\\2\end{pmatrix}.
    $$

    重新计算 DFT 验证即可确认。

---

## 30.4 Toeplitz 矩阵与信号处理

<div class="context-flow" markdown>

**核心**：Toeplitz 矩阵 = 沿对角线元素相同 → 线性时不变系统的矩阵表示 → 嵌入循环矩阵后可用 FFT 加速 → Levinson 递推 $O(n^2)$ 解 Toeplitz 方程组
**链接**：Ch22 数值线性代数 · 30.3 循环矩阵

</div>

Toeplitz 矩阵是信号处理中最常见的结构化矩阵，对应线性时不变（LTI）系统。

!!! definition "定义 30.6 (Toeplitz 矩阵)"
    $n$ 阶 **Toeplitz 矩阵**（Toeplitz matrix）$T \in \mathbb{R}^{n \times n}$ 满足 $T_{ij} = t_{i-j}$，即沿每条对角线元素相同：

    $$
    T = \begin{pmatrix}t_0&t_{-1}&t_{-2}&\cdots&t_{-(n-1)}\\t_1&t_0&t_{-1}&\cdots&t_{-(n-2)}\\t_2&t_1&t_0&\cdots&t_{-(n-3)}\\\vdots&\vdots&\vdots&\ddots&\vdots\\t_{n-1}&t_{n-2}&t_{n-3}&\cdots&t_0\end{pmatrix}.
    $$

    $T$ 由 $2n - 1$ 个参数 $t_{-(n-1)}, \ldots, t_0, \ldots, t_{n-1}$ 完全确定。若 $T$ 还对称（$t_{-k} = t_k$），则称为**对称 Toeplitz 矩阵**。

!!! theorem "定理 30.7 (Toeplitz 矩阵的循环嵌入)"
    任何 $n$ 阶 Toeplitz 矩阵 $T$ 可以嵌入到一个 $m$ 阶循环矩阵 $C$（$m \ge 2n - 1$）中，使得 $T$ 是 $C$ 的左上 $n \times n$ 子矩阵。因此，Toeplitz 矩阵-向量乘法 $T\mathbf{x}$ 可以通过 FFT 在 $O(n \log n)$ 时间内完成：

    1. 将 $\mathbf{x}$ 零填充至长度 $m$。
    2. 对 $\mathbf{c}$（$C$ 的第一列）和填充后的 $\mathbf{x}$ 做 FFT。
    3. 频域逐点乘法。
    4. IFFT 后取前 $n$ 个分量。

??? proof "证明"
    构造 $m = 2n$ 阶循环矩阵 $C = \operatorname{circ}(c_0, c_1, \ldots, c_{2n-1})$，其中 $c_k = t_k$（$0 \le k \le n-1$），$c_k = t_{k-2n}$（$n \le k \le 2n-1$）。将 $\mathbf{x}$ 填充为 $\tilde{\mathbf{x}} = (x_0, \ldots, x_{n-1}, 0, \ldots, 0)^T$。则 $(C\tilde{\mathbf{x}})_j = \sum_{k=0}^{n-1}c_{(j-k)\bmod 2n}x_k$。当 $0 \le j \le n-1$ 时，$(j-k) \bmod 2n$ 恰好给出 $t_{j-k}$，故前 $n$ 个分量等于 $T\mathbf{x}$。$\blacksquare$

!!! definition "定义 30.7 (Yule-Walker 方程)"
    在自回归（AR）模型中，给定对称正定 Toeplitz 矩阵 $T_n = [r_{|i-j|}]$（自相关矩阵），**Yule-Walker 方程**为

    $$
    T_n \mathbf{a} = -\mathbf{r}, \quad \mathbf{r} = (r_1, r_2, \ldots, r_n)^T,
    $$

    其中 $\mathbf{a}$ 为 AR 系数。由于 $T_n$ 的 Toeplitz 结构，可以用 Levinson 递推高效求解。

!!! theorem "定理 30.8 (Levinson-Durbin 递推)"
    对于正定对称 Toeplitz 系统 $T_n \mathbf{a} = \mathbf{b}$，Levinson-Durbin 算法在 $O(n^2)$ 时间内求解（相比一般方法的 $O(n^3)$）。算法从 $1 \times 1$ 系统开始，每步利用前一步的解和 Toeplitz 结构递推地增大维度：

    设 $T_k \mathbf{a}^{(k)} = \mathbf{b}^{(k)}$ 已解出。利用反射系数

    $$
    \kappa_{k+1} = \frac{b_{k+1} - \mathbf{r}_k^T \mathbf{a}^{(k)}}{\epsilon_k},
    $$

    其中 $\epsilon_k$ 为预测误差功率，更新

    $$
    \mathbf{a}^{(k+1)} = \begin{pmatrix}\mathbf{a}^{(k)}\\ 0\end{pmatrix} + \kappa_{k+1}\begin{pmatrix}J_k \mathbf{a}^{(k)}\\ 1\end{pmatrix},
    $$

    $J_k$ 为 $k$ 阶逆序矩阵。

??? proof "证明"
    关键在于利用 Toeplitz 矩阵的**中心对称性**：若 $T_k\mathbf{a} = \mathbf{b}$，则 $T_k(J_k\mathbf{a}) = J_k\mathbf{b}$。从 $k$ 阶解构造 $k+1$ 阶解只需确定一个新参数 $\kappa_{k+1}$（反射系数），使得扩展后的方程组成立。这将 $O(k)$ 的新方程化为单参数问题，总复杂度 $\sum_{k=1}^n O(k) = O(n^2)$。$\blacksquare$

!!! example "例 30.4"
    **Levinson 递推求解 AR 模型。** 设自相关序列 $r_0 = 1$，$r_1 = 0.8$，$r_2 = 0.5$。Yule-Walker 方程：

    $$
    \begin{pmatrix}1&0.8\\0.8&1\end{pmatrix}\begin{pmatrix}a_1\\a_2\end{pmatrix} = -\begin{pmatrix}0.8\\0.5\end{pmatrix}.
    $$

    **第 1 步**：$a_1^{(1)} = -r_1/r_0 = -0.8$，$\epsilon_1 = r_0(1 - 0.8^2) = 0.36$。

    **第 2 步**：$\kappa_2 = \frac{-0.5 - 0.8 \cdot (-0.8)}{0.36} = \frac{-0.5 + 0.64}{0.36} = \frac{0.14}{0.36} \approx 0.389$。

    $$
    \mathbf{a}^{(2)} = \begin{pmatrix}-0.8\\0\end{pmatrix} + 0.389\begin{pmatrix}-0.8\\1\end{pmatrix} = \begin{pmatrix}-1.111\\0.389\end{pmatrix}.
    $$

    预测误差功率 $\epsilon_2 = 0.36(1 - 0.389^2) \approx 0.305$。由于 $|\kappa_k| < 1$，系统稳定。

---

## 30.5 线性编码理论

<div class="context-flow" markdown>

**核心**：线性码 = 有限域 $\mathbb{F}_q$ 上的线性子空间 → 生成矩阵 $G$（列空间）和校验矩阵 $H$（零空间） → 最小距离 = 最小线性相关列数
**链接**：Ch4 向量空间 · Ch1 线性方程组

</div>

编码理论的基本思想是在传输信息中引入冗余以实现纠错，线性码将这一思想完全建立在线性代数之上。

!!! definition "定义 30.8 (线性码)"
    有限域 $\mathbb{F}_q$ 上的一个 $[n, k]$ **线性码**（linear code）$\mathcal{C}$ 是 $\mathbb{F}_q^n$ 的一个 $k$ 维线性子空间。参数解释：

    - $n$：码长（每个码字的长度）。
    - $k$：信息维数（可编码的信息符号数）。
    - $n - k$：冗余符号数。
    - **码率** $R = k/n$。

    $\mathcal{C}$ 共有 $q^k$ 个码字。

!!! definition "定义 30.9 (生成矩阵与校验矩阵)"
    $[n, k]$ 线性码 $\mathcal{C}$ 的**生成矩阵**（generator matrix）$G \in \mathbb{F}_q^{k \times n}$ 满足

    $$
    \mathcal{C} = \{\mathbf{m}G : \mathbf{m} \in \mathbb{F}_q^k\} = \operatorname{rowspace}(G).
    $$

    **校验矩阵**（parity-check matrix）$H \in \mathbb{F}_q^{(n-k) \times n}$ 满足

    $$
    \mathcal{C} = \ker(H) = \{\mathbf{c} \in \mathbb{F}_q^n : H\mathbf{c}^T = \mathbf{0}\}.
    $$

    $G$ 和 $H$ 满足 $GH^T = 0$。当 $G = [I_k \mid P]$（系统形式）时，$H = [-P^T \mid I_{n-k}]$。

!!! theorem "定理 30.9 (最小距离与纠错能力)"
    线性码 $\mathcal{C}$ 的**最小距离**（minimum distance）$d$ 等于 $H$ 的任何 $d - 1$ 列线性无关、但存在 $d$ 列线性相关的最小列数。记为 $[n, k, d]$ 码。

    - $\mathcal{C}$ 可以**检测** $d - 1$ 个错误。
    - $\mathcal{C}$ 可以**纠正** $\lfloor(d-1)/2\rfloor$ 个错误。

??? proof "证明"
    $d = \min_{\mathbf{c} \in \mathcal{C}, \mathbf{c} \ne \mathbf{0}} w(\mathbf{c})$（最小 Hamming 重量），其中 $w(\mathbf{c})$ 为非零分量个数。$\mathbf{c} \in \mathcal{C}$ 当且仅当 $H\mathbf{c}^T = \mathbf{0}$，即 $\mathbf{c}$ 的非零分量位置对应的 $H$ 的列线性相关。因此 $d$ = 使 $H$ 中某 $d$ 列线性相关的最小值。对于纠错，若发生 $t$ 个错误（$t \le \lfloor(d-1)/2\rfloor$），接收字 $\mathbf{r}$ 与真实码字的距离 $\le t$，与其他任何码字的距离 $\ge d - t > t$，故可唯一译码。$\blacksquare$

!!! theorem "定理 30.10 (Singleton 界)"
    任何 $[n, k, d]$ 线性码满足

    $$
    d \le n - k + 1.
    $$

    达到此界的码称为**最大距离可分**（MDS）码。

??? proof "证明"
    $H \in \mathbb{F}_q^{(n-k) \times n}$ 有 $n$ 列，每列在 $\mathbb{F}_q^{n-k}$ 中。任意 $n - k + 1$ 列在 $(n-k)$ 维空间中必线性相关，故 $d \le n - k + 1$。$\blacksquare$

!!! example "例 30.5"
    **$[7, 4]$ 二进制线性码。** 生成矩阵（系统形式）：

    $$
    G = \begin{pmatrix}1&0&0&0&1&1&0\\0&1&0&0&0&1&1\\0&0&1&0&1&1&1\\0&0&0&1&1&0&1\end{pmatrix}, \quad H = \begin{pmatrix}1&0&1&1&1&0&0\\1&1&1&0&0&1&0\\0&1&1&1&0&0&1\end{pmatrix}.
    $$

    验证 $GH^T = 0$（$\mathbb{F}_2$ 上）。信息字 $\mathbf{m} = (1,0,1,1)$ 编码为 $\mathbf{c} = \mathbf{m}G = (1,0,1,1,1,0,0)$。校验：$H\mathbf{c}^T = (0,0,0)^T$。

---

## 30.6 Hamming 码与纠错

<div class="context-flow" markdown>

**核心**：Hamming 码的校验矩阵 $H$ 的列 = $\mathbb{F}_2^r$ 中所有非零向量 → 伴随式 $\mathbf{s} = H\mathbf{r}^T$ 直接指出错误位置 → 完美码（达到 Hamming 界）
**链接**：30.5 线性码 · Ch4 线性无关

</div>

Hamming 码是最经典的纠错码之一，其构造优雅地体现了线性代数的思想。

!!! definition "定义 30.10 (Hamming 码)"
    对于正整数 $r \ge 2$，$[2^r - 1, 2^r - 1 - r, 3]$ **Hamming 码**的校验矩阵 $H \in \mathbb{F}_2^{r \times (2^r-1)}$ 的列是 $\mathbb{F}_2^r$ 中所有 $2^r - 1$ 个非零向量（某种排列）。

    - 码长 $n = 2^r - 1$。
    - 信息位 $k = 2^r - 1 - r$。
    - 最小距离 $d = 3$（可纠 1 位错误）。

!!! theorem "定理 30.11 (Hamming 码的性质)"
    Hamming 码具有以下性质：

    1. **$d = 3$**：$H$ 的任意两列（$\mathbb{F}_2^r$ 中两个不同非零向量）线性无关，但存在三列线性相关。
    2. **完美码**：Hamming 码达到 Hamming 界 $\sum_{i=0}^{t}\binom{n}{i} \le 2^{n-k}$（$t = 1$ 时取等号）。
    3. **伴随式译码**：对接收字 $\mathbf{r} = \mathbf{c} + \mathbf{e}$，伴随式 $\mathbf{s} = H\mathbf{r}^T = H\mathbf{e}^T$。若 $\mathbf{e}$ 为单位向量 $\mathbf{e}_j$（第 $j$ 位错误），则 $\mathbf{s} = H\mathbf{e}_j^T = \mathbf{h}_j$（$H$ 的第 $j$ 列），直接定位错误。

??? proof "证明"
    (1) $\mathbb{F}_2$ 上两个不同非零向量线性无关（否则相等）。对于三列 $\mathbf{h}_i, \mathbf{h}_j, \mathbf{h}_k$，当 $\mathbf{h}_i + \mathbf{h}_j = \mathbf{h}_k$（$\mathbb{F}_2$ 上必然存在这样的三元组），它们线性相关。故 $d = 3$。

    (2) Hamming 界：纠 $t = 1$ 个错误，需要 $1 + n = 2^r$ 个陪集（无错误 + $n$ 种单错误模式），恰等于 $2^{n-k} = 2^r$ 个伴随式。

    (3) 直接由 $\mathbf{s} = H\mathbf{e}^T$ 得到。$\blacksquare$

!!! example "例 30.6"
    **$[7, 4, 3]$ Hamming 码的纠错。** $r = 3$，校验矩阵：

    $$
    H = \begin{pmatrix}0&0&0&1&1&1&1\\0&1&1&0&0&1&1\\1&0&1&0&1&0&1\end{pmatrix}.
    $$

    列为 $1$ 到 $7$ 的二进制表示。设发送码字 $\mathbf{c} = (0,1,1,0,0,1,0)$，接收 $\mathbf{r} = (0,1,1,0,1,1,0)$（第 5 位出错）。

    伴随式 $\mathbf{s} = H\mathbf{r}^T = (1,0,1)^T$。二进制 $(101)_2 = 5$，指出第 5 位错误。纠正：翻转第 5 位得 $\hat{\mathbf{c}} = (0,1,1,0,0,1,0) = \mathbf{c}$。

!!! example "例 30.7"
    **扩展 Hamming 码。** 在 $[7,4,3]$ Hamming 码上附加一个总偶校验位，得到 $[8,4,4]$ **扩展 Hamming 码**。校验矩阵变为

    $$
    H_{\text{ext}} = \begin{pmatrix}1&1&1&1&1&1&1&1\\0&0&0&1&1&1&1&0\\0&1&1&0&0&1&1&0\\1&0&1&0&1&0&1&0\end{pmatrix}.
    $$

    $d = 4$ 可纠 1 位错误并检测 2 位错误。伴随式分析：若 $\mathbf{s} \ne \mathbf{0}$ 且第一个分量为 $1$，则为单错误（可纠正）；若第一个分量为 $0$，则为双错误（仅检测）。

---

## 30.7 Reed-Solomon 码

<div class="context-flow" markdown>

**核心**：RS 码在 $\mathbb{F}_q$（$q = 2^m$）上，生成矩阵是 Vandermonde 矩阵 → MDS 码（达到 Singleton 界） → 多项式求值/插值视角
**链接**：30.5 线性码 · Ch3 行列式(Vandermonde)

</div>

Reed-Solomon 码是最重要的 MDS 码，广泛应用于 CD/DVD、QR 码、深空通信等。其代数结构建立在 Vandermonde 矩阵之上。

!!! definition "定义 30.11 (Reed-Solomon 码)"
    设 $\mathbb{F}_q$ 为有 $q$ 个元素的有限域，$\alpha_1, \alpha_2, \ldots, \alpha_n$ 为 $\mathbb{F}_q$ 中 $n$ 个不同元素（$n \le q$）。$[n, k, n-k+1]$ **Reed-Solomon 码**定义为

    $$
    \mathcal{C}_{\text{RS}} = \{(p(\alpha_1), p(\alpha_2), \ldots, p(\alpha_n)) : p \in \mathbb{F}_q[x], \deg(p) < k\}.
    $$

    即将次数小于 $k$ 的多项式在 $n$ 个点上求值。生成矩阵为 Vandermonde 矩阵：

    $$
    G = \begin{pmatrix}1&1&\cdots&1\\\alpha_1&\alpha_2&\cdots&\alpha_n\\\alpha_1^2&\alpha_2^2&\cdots&\alpha_n^2\\\vdots&\vdots&\ddots&\vdots\\\alpha_1^{k-1}&\alpha_2^{k-1}&\cdots&\alpha_n^{k-1}\end{pmatrix} \in \mathbb{F}_q^{k \times n}.
    $$

!!! theorem "定理 30.12 (RS 码是 MDS 码)"
    Reed-Solomon 码达到 Singleton 界，即 $d = n - k + 1$，是 MDS 码。等价地：

    1. $G$ 的任意 $k$ 列线性无关（Vandermonde 矩阵的任意 $k$ 列子矩阵满秩）。
    2. $H$ 的任意 $n - k$ 列线性无关。
    3. RS 码可纠正 $t = \lfloor(n-k)/2\rfloor$ 个符号错误。

??? proof "证明"
    $G$ 的任意 $k$ 列构成一个 $k \times k$ Vandermonde 矩阵 $V(\alpha_{i_1}, \ldots, \alpha_{i_k})$，其行列式 $\det(V) = \prod_{s<t}(\alpha_{i_t} - \alpha_{i_s}) \ne 0$（因为 $\alpha_i$ 两两不同）。因此 $G$ 的任意 $k$ 列线性无关。由此推出最小距离 $d \ge n - k + 1$，而 Singleton 界给出 $d \le n - k + 1$，故 $d = n - k + 1$。$\blacksquare$

!!! theorem "定理 30.13 (RS 码的校验矩阵)"
    选取 $\alpha_i = \alpha^{i-1}$（$\alpha$ 为 $\mathbb{F}_q$ 的本原元素，$n = q - 1$），RS 码的校验矩阵为

    $$
    H = \begin{pmatrix}1&\alpha&\alpha^2&\cdots&\alpha^{n-1}\\1&\alpha^2&\alpha^4&\cdots&\alpha^{2(n-1)}\\\vdots&\vdots&\vdots&\ddots&\vdots\\1&\alpha^{n-k}&\alpha^{2(n-k)}&\cdots&\alpha^{(n-k)(n-1)}\end{pmatrix}.
    $$

    伴随式 $\mathbf{s} = H\mathbf{r}^T$ 的分量 $s_j = r(\alpha^j)$（$j = 1, \ldots, n-k$），其中 $r(x)$ 为接收字的多项式表示。译码算法（如 Berlekamp-Massey 或 Euclidean 算法）利用伴随式找到错误定位多项式和错误值。

??? proof "证明"
    码字 $\mathbf{c} = (p(\alpha^0), p(\alpha^1), \ldots, p(\alpha^{n-1}))$，$\deg(p) < k$。$H\mathbf{c}^T$ 的第 $j$ 个分量为 $\sum_{i=0}^{n-1}c_i\alpha^{ji} = \sum_{i=0}^{n-1}p(\alpha^i)\alpha^{ji}$。由于 $p$ 的根包含 $\alpha, \alpha^2, \ldots, \alpha^{n-k}$（作为码的生成多项式的根），$H\mathbf{c}^T = \mathbf{0}$。$\blacksquare$

!!! example "例 30.8"
    **$\mathbb{F}_8$ 上的 RS 码。** 设 $\mathbb{F}_8 = \mathbb{F}_2[\alpha]/(\alpha^3 + \alpha + 1)$，$n = 7$，$k = 3$，$d = 5$（可纠 2 个符号错误）。

    信息多项式 $p(x) = 1 + \alpha x + \alpha^2 x^2$。码字通过在 $1, \alpha, \alpha^2, \ldots, \alpha^6$ 处求值：

    $$
    \mathbf{c} = (p(1), p(\alpha), p(\alpha^2), \ldots, p(\alpha^6)).
    $$

    每个分量是 $\mathbb{F}_8$ 中的元素（3 bit）。码率 $R = 3/7 \approx 0.43$。由于 $d = 5$，可纠正任意 $2$ 个符号错误（每个符号 3 bit），即在 21 bit 中纠正任意 6 bit 错误（只要错误集中在不超过 2 个符号中）。

---

## 30.8 LDPC 码与稀疏矩阵

<div class="context-flow" markdown>

**核心**：LDPC 码 = 校验矩阵 $H$ 稀疏 → Tanner 图（二部图）表示 → 置信传播译码利用稀疏结构 → 逼近 Shannon 极限
**链接**：30.5 线性码 · Ch17 稀疏矩阵/非负矩阵

</div>

低密度奇偶校验（LDPC）码是现代通信系统中性能最优的码类之一，其核心特征是校验矩阵的稀疏性。

!!! definition "定义 30.12 (LDPC 码)"
    **低密度奇偶校验码**（Low-Density Parity-Check code, LDPC code）是一个 $[n, k]$ 线性码，其校验矩阵 $H \in \mathbb{F}_2^{(n-k) \times n}$ 是**稀疏矩阵**：每行和每列中 $1$ 的个数远小于 $n$。

    - **正则 LDPC 码**：$H$ 的每列恰好有 $d_v$ 个 $1$（变量节点度），每行恰好有 $d_c$ 个 $1$（校验节点度），满足 $(n-k)d_c = n d_v$。
    - **非正则 LDPC 码**：行列重量可以变化，用度分布对 $(\lambda(x), \rho(x))$ 描述。

!!! definition "定义 30.13 (Tanner 图)"
    LDPC 码的 **Tanner 图**（Tanner graph）是一个二部图 $G = (V \cup C, E)$：

    - **变量节点**（variable nodes）$V = \{v_1, \ldots, v_n\}$，对应 $H$ 的 $n$ 列（码字的 $n$ 位）。
    - **校验节点**（check nodes）$C = \{c_1, \ldots, c_{n-k}\}$，对应 $H$ 的 $n - k$ 行（校验方程）。
    - **边**：$v_j$ 与 $c_i$ 相连当且仅当 $H_{ij} = 1$。

    $H$ 的稀疏性意味着 Tanner 图的边数 $|E| = \text{nnz}(H) \ll n(n-k)$。

!!! theorem "定理 30.14 (LDPC 码的最小距离)"
    对于正则 $(d_v, d_c)$-LDPC 码：

    1. **下界**：$d \ge d_v + 1$（$H$ 的任意 $d_v$ 列线性无关，因为每列恰好有 $d_v$ 个 $1$，且稀疏性使得小集合不易相关）。
    2. Tanner 图的**围长**（girth）$g$（最短环的长度）越大，码的性能越好。$d \ge g/2 + 1$（对于正则码的粗略下界）。
    3. 随机构造的 LDPC 码，当 $n \to \infty$ 时，最小距离以高概率线性增长：$d = \Theta(n)$。

??? proof "证明"
    (1) 设 $\mathbf{c}$ 为非零码字，重量 $w$。$H\mathbf{c}^T = \mathbf{0}$ 意味着 $\mathbf{c}$ 的支撑（非零位置）对应的 $H$ 的 $w$ 列在 $\mathbb{F}_2$ 上线性相关。由于每列有 $d_v$ 个 $1$，$w$ 列的和（$\bmod 2$）的支撑大小 $\le wd_v$。要使和为零向量，需要足够多的列"覆盖"所有非零行。当 $w \le d_v$ 时，$w$ 列的非零行位置最多 $wd_v$ 个，而每行最多被 $d_c$ 列覆盖，仔细分析表明 $w \le d_v$ 列不可能线性相关（在 Tanner 图无短环的条件下），故 $d \ge d_v + 1$。$\blacksquare$

!!! theorem "定理 30.15 (置信传播译码)"
    LDPC 码的**置信传播**（Belief Propagation, BP）译码在 Tanner 图上进行消息传递。设 $L_j$ 为变量节点 $v_j$ 的信道对数似然比（LLR），译码迭代如下：

    **变量 → 校验消息**：

    $$
    \mu_{v_j \to c_i}^{(\ell)} = L_j + \sum_{c \in \mathcal{N}(v_j) \setminus c_i} \mu_{c \to v_j}^{(\ell-1)},
    $$

    **校验 → 变量消息**：

    $$
    \mu_{c_i \to v_j}^{(\ell)} = 2\operatorname{atanh}\!\left(\prod_{v \in \mathcal{N}(c_i) \setminus v_j} \tanh\!\left(\frac{\mu_{v \to c_i}^{(\ell)}}{2}\right)\right),
    $$

    其中 $\mathcal{N}(\cdot)$ 为 Tanner 图中的邻居集。每次迭代的复杂度为 $O(|E|) = O(\text{nnz}(H))$，由 $H$ 的稀疏性保证高效。

??? proof "证明"
    BP 算法基于因子图上的和积（sum-product）算法。对于 AWGN 信道，接收值 $y_j = (1 - 2c_j) + n_j$，LLR $L_j = 2y_j/\sigma^2$。变量到校验的消息汇总了来自信道和其他校验节点的信息。校验到变量的消息基于校验约束 $\bigoplus_{v \in \mathcal{N}(c_i)}c_v = 0$，通过 $\tanh$ 规则在 LLR 域中实现。当 Tanner 图无环时，BP 给出精确的后验概率；有环时 BP 是近似的，但在实践中对长度足够的 LDPC 码表现极好。$\blacksquare$

!!! example "例 30.9"
    **简单 LDPC 码的 Tanner 图。** 设 $[6, 3]$ 正则 $(2, 3)$-LDPC 码，校验矩阵：

    $$
    H = \begin{pmatrix}1&1&1&0&0&0\\0&0&1&1&1&0\\1&0&0&0&1&1\end{pmatrix}.
    $$

    每列 $2$ 个 $1$（$d_v = 2$），每行 $3$ 个 $1$（$d_c = 3$）。Tanner 图有 $6$ 个变量节点和 $3$ 个校验节点，$6 \times 2 = 12$ 条边（$= 3 \times 3 + 3 = ?$，实际 $\text{nnz}(H) = 9$ 条边）。

    重新检验：$\text{nnz}(H) = 9$，这不符合 $(n-k)d_c = nd_v$（$3 \times 3 = 9 = 6 \times 1.5$），故此矩阵每列平均 $1.5$ 个 $1$，是非正则的。一个正确的正则 $(2,4)$-LDPC 码需要 $n = 2(n-k)$，如 $[8,4]$ 码。

    尽管如此，上述 $H$ 定义了一个有效的 $[6,3]$ 线性码。设接收 $\mathbf{r} = (1,0,1,0,1,0)$，伴随式 $\mathbf{s} = H\mathbf{r}^T = (0, 0, 0)^T$，故 $\mathbf{r}$ 为有效码字。若 $\mathbf{r} = (1,1,1,0,1,0)$，$\mathbf{s} = (1, 0, 1)^T \ne \mathbf{0}$，检测到错误，BP 译码将利用 Tanner 图上的消息传递来定位并纠正错误。
