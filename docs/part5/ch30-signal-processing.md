# 第 30 章 线性代数在信号处理与编码中的应用

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵分解(Ch10-11) · 正定矩阵(Ch16) · Kronecker 积(Ch19)

**脉络**：DFT 矩阵(酉变换) → FFT(稀疏矩阵分解) → 循环矩阵与卷积(DFT 对角化) → 滤波与频域分析(频率响应) → 采样定理与插值(带限信号) → 压缩感知(稀疏恢复/RIP) → 线性纠错码(有限域上的线性代数) → LDPC 码(稀疏校验矩阵) → 小波变换(多分辨率分析)

</div>

信号处理与编码理论是线性代数最经典的工程应用领域。离散 Fourier 变换本质上是酉矩阵乘法，循环卷积对应循环矩阵，滤波是频域中的对角矩阵运算，采样与插值由 Nyquist-Shannon 定理精确刻画。在编码侧，纠错码的设计与译码完全建立在有限域上的线性代数之上；LDPC 码利用稀疏矩阵结构实现逼近 Shannon 极限的性能。本章最后介绍压缩感知和小波变换，展示线性代数在现代信号处理中的前沿应用。

---

## 30.1 离散 Fourier 变换

<div class="context-flow" markdown>

**线性代数视角**：DFT = 乘以酉矩阵 $F_n$（$F_n^* F_n = nI$）→ 特征值/特征向量全部已知 → DFT 将循环卷积对角化 → FFT 将 $F_n$ 分解为 $\log_2 n$ 个稀疏矩阵的乘积

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

!!! definition "定义 30.2 (逆 DFT 与 FFT)"
    **逆 DFT**（inverse DFT）为 $\mathbf{x} = \frac{1}{n}F_n^* \hat{\mathbf{x}}$，其中 $F_n^* = \overline{F_n}^T$。

    当 $n = 2^s$ 时，**快速 Fourier 变换**（FFT）通过 Cooley-Tukey 分解

    $$
    F_n = \begin{pmatrix}I_m & D_m \\ I_m & -D_m\end{pmatrix}\begin{pmatrix}F_m & 0 \\ 0 & F_m\end{pmatrix}P_n, \quad m = n/2,
    $$

    将 $n$ 点 DFT 分解为两个 $n/2$ 点 DFT 加旋转因子乘法。其中 $D_m = \operatorname{diag}(1, \omega_n, \ldots, \omega_n^{m-1})$，$P_n$ 为完美洗牌置换。递推后总复杂度 $O(n \log n)$。

!!! theorem "定理 30.1 (DFT 矩阵的酉性)"
    归一化 DFT 矩阵 $\frac{1}{\sqrt{n}}F_n$ 是酉矩阵，即

    $$
    F_n^* F_n = F_n F_n^* = nI_n.
    $$

    等价地，$F_n$ 的列（或行）构成 $\mathbb{C}^n$ 的正交基（内积为 $n$），$\frac{1}{\sqrt{n}}F_n$ 的列构成标准正交基。

??? proof "证明"
    $(F_n^* F_n)_{jk} = \sum_{\ell=0}^{n-1}\overline{(F_n)_{\ell j}}(F_n)_{\ell k} = \sum_{\ell=0}^{n-1}\omega_n^{-\ell j}\omega_n^{\ell k} = \sum_{\ell=0}^{n-1}\omega_n^{\ell(k-j)}$。当 $k = j$ 时，和为 $n$。当 $k \ne j$ 时，这是等比级数 $\sum_{\ell=0}^{n-1}(\omega_n^{k-j})^\ell = \frac{1 - \omega_n^{n(k-j)}}{1 - \omega_n^{k-j}} = 0$，因为 $\omega_n^n = 1$ 且 $\omega_n^{k-j} \ne 1$（$0 < |k-j| < n$）。$\blacksquare$

!!! theorem "定理 30.2 (DFT 矩阵的特征值与 FFT 复杂度)"
    1. $F_n$ 满足 $F_n^4 = n^2 I_n$，故特征值为 $\{+\sqrt{n}, -\sqrt{n}, +i\sqrt{n}, -i\sqrt{n}\}$。
    2. 当 $n = 2^s$ 时，FFT 的乘法次数为 $\frac{n}{2}\log_2 n$，加法次数为 $n\log_2 n$，总复杂度 $O(n\log n)$。FFT 将 $F_n$ 分解为 $\log_2 n$ 个稀疏矩阵的乘积，每个稀疏矩阵有 $O(n)$ 个非零元素。

??? proof "证明"
    (1) $(F_n^2)_{jk} = \sum_\ell \omega_n^{\ell(j+k)}$。当 $j + k \equiv 0 \pmod{n}$ 时和为 $n$，否则为 $0$，故 $F_n^2 = nP$（$P$ 为逆序置换矩阵）。$F_n^4 = n^2 P^2 = n^2 I$，特征值 $\lambda$ 满足 $\lambda^4 = n^2$。

    (2) 设 $T(n)$ 为乘法次数。Cooley-Tukey 分解给出 $T(n) = 2T(n/2) + n/2$，$T(1) = 0$，解为 $T(n) = \frac{n}{2}\log_2 n$。矩阵分解形式 $F_n = B_s B_{s-1} \cdots B_1 P$，每个 $B_k$ 有 $O(n)$ 个非零元素。$\blacksquare$

!!! example "例 30.1"
    **$4$ 阶 DFT 矩阵与 FFT。** $\omega_4 = e^{-\pi i/2} = -i$。

    $$
    F_4 = \begin{pmatrix}1&1&1&1\\1&-i&-1&i\\1&-1&1&-1\\1&i&-1&-i\end{pmatrix}.
    $$

    验证：$F_4^* F_4 = 4I_4$。对 $\mathbf{x} = (1, 2, 3, 4)^T$ 做 DFT：

    $$
    \hat{\mathbf{x}} = F_4\mathbf{x} = \begin{pmatrix}10\\-2+2i\\-2\\-2-2i\end{pmatrix}.
    $$

    FFT 分解为 2 层蝶形运算，仅需 $\frac{4}{2}\log_2 4 = 4$ 次乘法（直接法需 $16$ 次）。对于 $n = 2^{20}$，FFT 加速比约 $10^5$。

---

## 30.2 循环矩阵与卷积

<div class="context-flow" markdown>

**核心**：循环矩阵 = DFT 矩阵对角化 $C = \frac{1}{n}F_n^*\operatorname{diag}(\hat{\mathbf{c}})F_n$ → 循环卷积 = 频域逐点乘法 → 卷积定理 → "FFT → 逐点乘 → IFFT" $O(n\log n)$ 加速

**链接**：Ch6 特征值 · 30.1 DFT

</div>

循环矩阵是一类具有特殊结构的矩阵，它与 DFT 和卷积运算有深刻的联系。

!!! definition "定义 30.3 (循环矩阵)"
    由向量 $\mathbf{c} = (c_0, c_1, \ldots, c_{n-1})^T$ 生成的 $n$ 阶**循环矩阵**（circulant matrix）为

    $$
    C = \operatorname{circ}(\mathbf{c}) = \begin{pmatrix}c_0&c_{n-1}&c_{n-2}&\cdots&c_1\\c_1&c_0&c_{n-1}&\cdots&c_2\\c_2&c_1&c_0&\cdots&c_3\\\vdots&\vdots&\vdots&\ddots&\vdots\\c_{n-1}&c_{n-2}&c_{n-3}&\cdots&c_0\end{pmatrix}.
    $$

    每一行是上一行的循环右移。等价地，$C = \sum_{k=0}^{n-1}c_k P^k$，其中 $P$ 为循环置换矩阵。

!!! definition "定义 30.4 (循环卷积)"
    两个长度为 $n$ 的向量 $\mathbf{a}$ 和 $\mathbf{b}$ 的**循环卷积**（circular convolution）$\mathbf{c} = \mathbf{a} \circledast \mathbf{b}$ 定义为

    $$
    c_j = \sum_{k=0}^{n-1} a_k b_{(j-k) \bmod n}, \quad j = 0, 1, \ldots, n-1.
    $$

    循环卷积等价于循环矩阵乘法：$\mathbf{c} = \operatorname{circ}(\mathbf{a})\mathbf{b}$。

!!! theorem "定理 30.3 (循环矩阵的谱分解与卷积定理)"
    1. **谱分解**：所有 $n$ 阶循环矩阵都被同一个 DFT 矩阵对角化：

        $$
        C = \operatorname{circ}(\mathbf{c}) = \frac{1}{n}F_n^* \operatorname{diag}(\hat{\mathbf{c}}) F_n,
        $$

        其中 $\hat{\mathbf{c}} = F_n \mathbf{c}$。$C$ 的特征值为 $\hat{c}_0, \hat{c}_1, \ldots, \hat{c}_{n-1}$，特征向量为 $F_n$ 的列。

    2. **卷积定理**：循环卷积在频域等价于逐点乘法：$\mathbf{a} \circledast \mathbf{b} = \mathbf{c} \iff \hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}}$，其中 $\odot$ 为逐元素乘法。

    因此循环卷积可通过 "FFT → 逐点乘法 → IFFT" 在 $O(n\log n)$ 时间内完成。所有循环矩阵两两可交换。

??? proof "证明"
    (1) 循环置换矩阵 $P$ 的特征值为 $1, \omega_n, \omega_n^2, \ldots, \omega_n^{n-1}$，特征向量为 DFT 矩阵的列。由于 $C = \sum_k c_k P^k$，$C$ 与 $P$ 共享特征向量，特征值为 $\lambda_j = \sum_k c_k \omega_n^{jk} = \hat{c}_j$。两个循环矩阵都被 $F_n$ 对角化，故两两可交换。

    (2) $\mathbf{c} = \operatorname{circ}(\mathbf{a})\mathbf{b} = \frac{1}{n}F_n^*\operatorname{diag}(\hat{\mathbf{a}})F_n\mathbf{b} = \frac{1}{n}F_n^*(\hat{\mathbf{a}} \odot \hat{\mathbf{b}})$。取 DFT：$\hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}}$。$\blacksquare$

!!! example "例 30.2"
    **用 FFT 计算循环卷积。** 设 $\mathbf{a} = (1, 2, 3, 0)^T$，$\mathbf{b} = (1, 0, 1, 0)^T$。

    DFT：$\hat{\mathbf{a}} = F_4\mathbf{a} = (6, -2+2i, -2, -2-2i)^T$，$\hat{\mathbf{b}} = F_4\mathbf{b} = (2, 0, 2, 0)^T$。

    频域乘法：$\hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}} = (12, 0, -4, 0)^T$。

    IFFT：$\mathbf{c} = \frac{1}{4}F_4^*\hat{\mathbf{c}} = (2, 2, 4, 4)^T$。

    用循环矩阵验证：$\operatorname{circ}(\mathbf{a})\mathbf{b} = \begin{pmatrix}1&0&3&2\\2&1&0&3\\3&2&1&0\\0&3&2&1\end{pmatrix}\begin{pmatrix}1\\0\\1\\0\end{pmatrix} = \begin{pmatrix}4\\2\\4\\2\end{pmatrix}$。

---

## 30.3 滤波与频域分析

<div class="context-flow" markdown>

**核心**：线性时不变（LTI）系统 = Toeplitz 矩阵乘法 → 频率响应 $H(\omega)$ = 传递函数 → FIR 滤波器 = 有限卷积 → IIR 滤波器 = 递归差分方程 → Toeplitz 矩阵的循环嵌入与 FFT 加速

**链接**：30.2 循环矩阵 · Ch22 数值线性代数

</div>

滤波是信号处理的核心操作，其数学本质是矩阵-向量乘法。

!!! definition "定义 30.5 (Toeplitz 矩阵与 LTI 系统)"
    $n$ 阶 **Toeplitz 矩阵**（Toeplitz matrix）$T \in \mathbb{R}^{n \times n}$ 满足 $T_{ij} = t_{i-j}$，即沿每条对角线元素相同：

    $$
    T = \begin{pmatrix}t_0&t_{-1}&t_{-2}&\cdots&t_{-(n-1)}\\t_1&t_0&t_{-1}&\cdots&t_{-(n-2)}\\t_2&t_1&t_0&\cdots&t_{-(n-3)}\\\vdots&\vdots&\vdots&\ddots&\vdots\\t_{n-1}&t_{n-2}&t_{n-3}&\cdots&t_0\end{pmatrix}.
    $$

    Toeplitz 矩阵是**线性时不变**（LTI）系统的矩阵表示：输出 $y_j = \sum_k t_{j-k} x_k$ 是输入 $\mathbf{x}$ 与脉冲响应 $\{t_k\}$ 的卷积。

!!! definition "定义 30.6 (频率响应与滤波器)"
    LTI 系统的**频率响应**（frequency response）$H(\omega) = \sum_{k} h_k e^{-i\omega k}$ 是脉冲响应的 DTFT。

    - **FIR 滤波器**（有限脉冲响应）：$h_k$ 仅对有限个 $k$ 非零，$y = H * x$ 为有限卷积。矩阵表示为带状 Toeplitz 矩阵。
    - **IIR 滤波器**（无限脉冲响应）：由递归差分方程 $\sum_{k=0}^{N} a_k y_{n-k} = \sum_{k=0}^{M} b_k x_{n-k}$ 定义。传递函数为有理函数 $H(z) = B(z)/A(z)$，稳定性条件为 $A(z)$ 的根均在单位圆内。

!!! theorem "定理 30.4 (Toeplitz 矩阵的循环嵌入)"
    任何 $n$ 阶 Toeplitz 矩阵 $T$ 可以嵌入到 $m$ 阶循环矩阵 $C$（$m \ge 2n - 1$）中，使得 $T$ 是 $C$ 的左上 $n \times n$ 子矩阵。因此，Toeplitz 矩阵-向量乘法 $T\mathbf{x}$ 可通过 FFT 在 $O(n \log n)$ 时间内完成：

    1. 将 $\mathbf{x}$ 零填充至长度 $m$；
    2. 对 $\mathbf{c}$（$C$ 的第一列）和填充后的 $\mathbf{x}$ 做 FFT；
    3. 频域逐点乘法；
    4. IFFT 后取前 $n$ 个分量。

??? proof "证明"
    构造 $m = 2n$ 阶循环矩阵 $C = \operatorname{circ}(c_0, c_1, \ldots, c_{2n-1})$，其中 $c_k = t_k$（$0 \le k \le n-1$），$c_k = t_{k-2n}$（$n \le k \le 2n-1$）。将 $\mathbf{x}$ 填充为 $\tilde{\mathbf{x}} = (x_0, \ldots, x_{n-1}, 0, \ldots, 0)^T$。则 $(C\tilde{\mathbf{x}})_j = \sum_{k=0}^{n-1}c_{(j-k)\bmod 2n}x_k$。当 $0 \le j \le n-1$ 时，$(j-k) \bmod 2n$ 恰好给出 $t_{j-k}$，故前 $n$ 个分量等于 $T\mathbf{x}$。$\blacksquare$

!!! theorem "定理 30.5 (Levinson-Durbin 递推)"
    对于正定对称 Toeplitz 系统 $T_n \mathbf{a} = \mathbf{b}$，Levinson-Durbin 算法利用 Toeplitz 结构的中心对称性，从 $1 \times 1$ 系统递推地扩展维度，在 $O(n^2)$ 时间内求解（相比一般方法的 $O(n^3)$）。关键递推式利用反射系数

    $$
    \kappa_{k+1} = \frac{b_{k+1} - \mathbf{r}_k^T \mathbf{a}^{(k)}}{\epsilon_k},
    $$

    更新 $\mathbf{a}^{(k+1)} = \begin{pmatrix}\mathbf{a}^{(k)}\\ 0\end{pmatrix} + \kappa_{k+1}\begin{pmatrix}J_k \mathbf{a}^{(k)}\\ 1\end{pmatrix}$，其中 $J_k$ 为逆序矩阵，$\epsilon_k$ 为预测误差功率。

??? proof "证明"
    关键在于 Toeplitz 矩阵的中心对称性：若 $T_k\mathbf{a} = \mathbf{b}$，则 $T_k(J_k\mathbf{a}) = J_k\mathbf{b}$。从 $k$ 阶解构造 $k+1$ 阶解只需确定一个新参数 $\kappa_{k+1}$。将 $O(k)$ 的新方程化为单参数问题，总复杂度 $\sum_{k=1}^n O(k) = O(n^2)$。稳定性由反射系数满足 $|\kappa_k| < 1$ 保证。$\blacksquare$

!!! example "例 30.3"
    **FIR 低通滤波器的矩阵表示。** 长度 $3$ 的移动平均滤波器 $h = \frac{1}{3}(1, 1, 1)$，作用于长度 $6$ 的信号 $\mathbf{x}$：

    $$
    T = \frac{1}{3}\begin{pmatrix}1&0&0&0&0&0\\1&1&0&0&0&0\\1&1&1&0&0&0\\0&1&1&1&0&0\\0&0&1&1&1&0\\0&0&0&1&1&1\end{pmatrix}.
    $$

    这是下三角带状 Toeplitz 矩阵。频率响应 $H(\omega) = \frac{1}{3}(1 + e^{-i\omega} + e^{-2i\omega}) = \frac{1}{3}e^{-i\omega}\frac{\sin(3\omega/2)}{\sin(\omega/2)}$，在 $\omega = 0$ 时 $|H| = 1$（通过直流），在高频时 $|H|$ 衰减（低通特性）。

    利用循环嵌入，可将 $T\mathbf{x}$ 的计算转化为 FFT 加速的循环卷积。

---

## 30.4 采样定理与插值

<div class="context-flow" markdown>

**核心**：Shannon-Nyquist 定理 = 带限函数由其采样点唯一确定 → 采样矩阵与 sinc 插值 → 混叠 = 欠采样导致的频谱折叠 → 插值公式的线性代数表示

**链接**：30.1 DFT · Ch11 SVD/伪逆

</div>

采样定理是连接连续信号与离散信号的桥梁，其数学表述本质上是线性代数中的基展开问题。

!!! definition "定义 30.7 (带限信号与采样)"
    **带限信号**（bandlimited signal）$f(t)$ 的 Fourier 变换 $\hat{f}(\omega)$ 满足 $\hat{f}(\omega) = 0$（$|\omega| > W$），即带宽为 $W$。以采样间隔 $T_s = 1/(2W)$ 采样得到离散序列 $f_n = f(nT_s)$。

    **采样算子** $S : f \mapsto (f(nT_s))_{n \in \mathbb{Z}}$ 是从连续函数空间到离散序列空间的线性映射。

!!! theorem "定理 30.6 (Shannon-Nyquist 采样定理)"
    若 $f(t)$ 为带宽 $W$ 的带限信号（$\hat{f}(\omega) = 0$，$|\omega| > W$），则 $f$ 由其采样序列 $f_n = f(n/(2W))$ 唯一确定，且可通过 **sinc 插值**精确恢复：

    $$
    f(t) = \sum_{n=-\infty}^{\infty} f_n \operatorname{sinc}(2Wt - n),
    $$

    其中 $\operatorname{sinc}(x) = \frac{\sin(\pi x)}{\pi x}$。采样频率 $f_s = 2W$ 称为 **Nyquist 率**。

    当采样频率 $f_s < 2W$ 时，发生**混叠**（aliasing）：高频分量被折叠到低频，信号无法精确恢复。

??? proof "证明"
    在频域中，采样将 $\hat{f}(\omega)$ 在频率轴上以 $f_s = 2W$ 为周期进行周期化：$\hat{f}_s(\omega) = \frac{1}{T_s}\sum_k \hat{f}(\omega - kf_s)$。当 $f_s \ge 2W$ 时，各周期化副本不重叠，可以通过理想低通滤波器（$|\omega| \le W$ 时通过，否则截止）恢复 $\hat{f}(\omega)$。理想低通滤波器的时域形式即为 $\operatorname{sinc}$ 函数，故恢复公式为 sinc 插值。当 $f_s < 2W$ 时副本重叠，产生混叠。$\blacksquare$

!!! definition "定义 30.8 (离散 sinc 插值矩阵)"
    对于有限长度 $N$ 的信号，从 $M$ 个采样点插值到 $N$ 个点的插值矩阵 $S \in \mathbb{R}^{N \times M}$ 为

    $$
    S_{jk} = \operatorname{sinc}\!\left(\frac{j - k \cdot N/M}{1}\right), \quad j = 0, \ldots, N-1, \quad k = 0, \ldots, M-1.
    $$

    当 $M \ge N$（过采样）时 $S$ 的列多于行，系统过定，可用最小二乘（伪逆 $S^+$）求解。当 $M < N$（欠采样）时 $S$ 列少于行，系统欠定。

!!! example "例 30.4"
    **从 4 个采样点恢复 8 点信号。** 设带限信号的 DFT 仅在低频分量 $\hat{x}_0, \hat{x}_1, \hat{x}_6, \hat{x}_7$ 非零（带宽 $W = 2$，Nyquist 率 $= 4$）。采样 $4$ 个点得到 $\mathbf{y} = (y_0, y_1, y_2, y_3)^T$。

    恢复过程：对 $\mathbf{y}$ 做 $4$ 点 DFT 得到 $\hat{\mathbf{y}}$；将 $\hat{\mathbf{y}}$ 零填充为 $8$ 点频谱 $\hat{\mathbf{x}} = (\hat{y}_0, \hat{y}_1, 0, 0, 0, 0, \hat{y}_2, \hat{y}_3)^T$；对 $\hat{\mathbf{x}}$ 做 $8$ 点 IFFT 得到恢复信号。

    线性代数视角：这是一个降采样矩阵 $D \in \mathbb{R}^{4 \times 8}$（每隔一行取一行的单位矩阵）与 DFT 矩阵的组合：$\mathbf{y} = D\mathbf{x}$，恢复需要利用 $\mathbf{x}$ 的带限结构（即 DFT 域的稀疏性）。

---

## 30.5 压缩感知

<div class="context-flow" markdown>

**核心**：稀疏信号可以从远少于 Nyquist 率的测量中精确恢复 → 测量矩阵 $A$（$m \ll n$）→ 受限等距性质（RIP）→ $\ell_1$ 最小化（basis pursuit）等价于 $\ell_0$ 最小化 → 恢复条件与矩阵的几何性质

**链接**：Ch11 SVD · Ch10 矩阵范数 · Ch25 线性规划（$\ell_1$ 最小化）

</div>

压缩感知（compressed sensing）是 21 世纪初最具影响力的数学发现之一，它突破了 Shannon-Nyquist 采样定理的限制，利用信号的稀疏性实现亚 Nyquist 采样。

!!! definition "定义 30.9 (稀疏信号与压缩感知问题)"
    向量 $\mathbf{x} \in \mathbb{R}^n$ 称为 $s$**-稀疏**的，若 $\|\mathbf{x}\|_0 := |\{i : x_i \ne 0\}| \le s$。**压缩感知问题**：给定测量矩阵 $A \in \mathbb{R}^{m \times n}$（$m \ll n$）和测量向量 $\mathbf{y} = A\mathbf{x}$，在已知 $\mathbf{x}$ 是 $s$-稀疏的条件下恢复 $\mathbf{x}$。

    直接求解 $\ell_0$ 最小化 $\min \|\mathbf{x}\|_0 \text{ s.t. } A\mathbf{x} = \mathbf{y}$ 是 NP 困难问题。压缩感知的核心发现是：在适当条件下，可以用凸松弛（$\ell_1$ 最小化）精确恢复。

!!! definition "定义 30.10 (受限等距性质)"
    矩阵 $A \in \mathbb{R}^{m \times n}$ 满足 $s$ 阶**受限等距性质**（Restricted Isometry Property, RIP），常数为 $\delta_s \in [0, 1)$，若对所有 $s$-稀疏向量 $\mathbf{x}$，

    $$
    (1 - \delta_s)\|\mathbf{x}\|_2^2 \le \|A\mathbf{x}\|_2^2 \le (1 + \delta_s)\|\mathbf{x}\|_2^2.
    $$

    RIP 要求 $A$ 的任意 $s$ 列构成的子矩阵都近似保持范数，即 $A$ 在所有 $s$-稀疏方向上都近似等距。

!!! theorem "定理 30.7 (Candes-Tao RIP 恢复定理)"
    设 $A$ 满足 $2s$ 阶 RIP，常数 $\delta_{2s} < \sqrt{2} - 1$。则对任何 $s$-稀疏信号 $\mathbf{x}$，从 $\mathbf{y} = A\mathbf{x}$ 通过 **basis pursuit**

    $$
    \min_{\mathbf{z}} \|\mathbf{z}\|_1 \quad \text{s.t.} \quad A\mathbf{z} = \mathbf{y}
    $$

    可以精确恢复 $\mathbf{x}$。即 $\ell_1$ 最小化的解等于 $\ell_0$ 最小化的解。

    对于噪声观测 $\mathbf{y} = A\mathbf{x} + \mathbf{e}$（$\|\mathbf{e}\|_2 \le \epsilon$），使用 **LASSO** 或 **basis pursuit denoising** $\min \|\mathbf{z}\|_1 \text{ s.t. } \|A\mathbf{z} - \mathbf{y}\|_2 \le \epsilon$，恢复误差与噪声水平成比例。

??? proof "证明"
    设 $\mathbf{h} = \hat{\mathbf{x}} - \mathbf{x}$（$\hat{\mathbf{x}}$ 为 $\ell_1$ 最小化的解）。由 $A\mathbf{h} = \mathbf{0}$ 和 $\|\hat{\mathbf{x}}\|_1 \le \|\mathbf{x}\|_1$（最优性），将 $\mathbf{h}$ 按分量大小分解为 $\mathbf{h}_{T_0}$（$\mathbf{x}$ 的支撑上的分量）和 $\mathbf{h}_{T_0^c}$。$\ell_1$ 最优性给出 $\|\mathbf{h}_{T_0^c}\|_1 \le \|\mathbf{h}_{T_0}\|_1$。

    将 $\mathbf{h}_{T_0^c}$ 按大小递减分成大小为 $s$ 的块 $T_1, T_2, \ldots$。RIP 给出 $\|A\mathbf{h}_{T_0 \cup T_1}\|_2^2$ 的上下界，结合 $A\mathbf{h} = \mathbf{0}$，逐步推导出 $\|\mathbf{h}\|_2 = 0$，即 $\hat{\mathbf{x}} = \mathbf{x}$。$\delta_{2s} < \sqrt{2} - 1$ 是使推导闭合的关键条件。$\blacksquare$

!!! theorem "定理 30.8 (随机测量矩阵的 RIP)"
    以下随机矩阵以高概率满足 $s$ 阶 RIP（$\delta_s \le \delta$），只要行数 $m \ge C \cdot s \ln(n/s) / \delta^2$：

    1. **高斯随机矩阵**：$A_{ij} \sim \mathcal{N}(0, 1/m)$ i.i.d.。
    2. **Bernoulli 随机矩阵**：$A_{ij} = \pm 1/\sqrt{m}$ 等概率。
    3. **随机部分 Fourier 矩阵**：从 $F_n$ 中均匀随机选取 $m$ 行并归一化。

    关键结论：$m = O(s \log(n/s))$ 个测量即可恢复 $s$-稀疏的 $n$ 维信号，远少于 Nyquist 要求的 $n$ 个测量。

??? proof "证明"
    对高斯矩阵，固定一个 $s$-稀疏方向 $\mathbf{x}$，$\|A\mathbf{x}\|_2^2$ 是 $m$ 个独立高斯随机变量的和除以 $m$，集中在 $\|\mathbf{x}\|_2^2$ 附近（Chernoff 界）。对所有 $s$-稀疏方向取并界：$\binom{n}{s}$ 种支撑选择，每种上的单位球用 $\epsilon$-网覆盖，总共需 $(9/\delta)^s\binom{n}{s}$ 个点。当 $m \ge C s\ln(n/s)/\delta^2$ 时，并界保证所有稀疏方向上的偏差 $\le \delta$。$\blacksquare$

!!! example "例 30.5"
    **压缩感知恢复稀疏信号。** 设 $n = 100$，$s = 5$（$5$-稀疏信号），测量数 $m = 30$。

    生成 $30 \times 100$ 高斯测量矩阵 $A$，$\mathbf{y} = A\mathbf{x}$。求解 $\ell_1$ 最小化（线性规划）：

    $$
    \min_{\mathbf{z}} \|\mathbf{z}\|_1 \quad \text{s.t.} \quad A\mathbf{z} = \mathbf{y}.
    $$

    理论保证：$m = 30 \ge C \cdot 5 \cdot \ln(100/5) \approx 5C \cdot 3 = 15C$。取 $C = 2$ 时 $m = 30$ 足够。实践中，$30$ 个测量即可精确恢复 $100$ 维的 $5$-稀疏信号，压缩比 $30/100 = 30\%$。

---

## 30.6 纠错码的线性代数基础

<div class="context-flow" markdown>

**核心**：线性码 = 有限域 $\mathbb{F}_q$ 上的线性子空间 → 生成矩阵 $G$（列空间）和校验矩阵 $H$（零空间） → 最小距离 = $H$ 的最小线性相关列数 → Hamming 码与 Reed-Solomon 码

**链接**：Ch4 向量空间 · Ch1 线性方程组

</div>

编码理论的基本思想是在传输信息中引入冗余以实现纠错，线性码将这一思想完全建立在线性代数之上。

!!! definition "定义 30.11 (线性码与生成矩阵/校验矩阵)"
    有限域 $\mathbb{F}_q$ 上的 $[n, k]$ **线性码** $\mathcal{C}$ 是 $\mathbb{F}_q^n$ 的 $k$ 维子空间，有 $q^k$ 个码字。**生成矩阵** $G \in \mathbb{F}_q^{k \times n}$ 满足 $\mathcal{C} = \operatorname{rowspace}(G)$，**校验矩阵** $H \in \mathbb{F}_q^{(n-k) \times n}$ 满足 $\mathcal{C} = \ker(H)$。

    $G$ 和 $H$ 满足 $GH^T = 0$。系统形式：$G = [I_k \mid P]$ 时 $H = [-P^T \mid I_{n-k}]$。

!!! theorem "定理 30.9 (最小距离、Singleton 界与 MDS 码)"
    1. **最小距离**：$[n, k, d]$ 线性码的最小距离 $d$ 等于 $H$ 的任何 $d-1$ 列线性无关但存在 $d$ 列线性相关的最小列数。$\mathcal{C}$ 可检测 $d-1$ 个错误，纠正 $\lfloor(d-1)/2\rfloor$ 个错误。

    2. **Singleton 界**：$d \le n - k + 1$。达到此界的码称为 **MDS 码**。

    3. **Hamming 界**：$\sum_{i=0}^{t}\binom{n}{i}(q-1)^i \le q^{n-k}$，其中 $t = \lfloor(d-1)/2\rfloor$。达到此界的码称为**完美码**。

??? proof "证明"
    (1) $d = \min_{\mathbf{c} \ne \mathbf{0}, \mathbf{c} \in \mathcal{C}} w(\mathbf{c})$（最小 Hamming 重量）。$\mathbf{c} \in \mathcal{C}$ 当且仅当 $H\mathbf{c}^T = \mathbf{0}$，即非零分量对应的 $H$ 列线性相关。

    (2) $H \in \mathbb{F}_q^{(n-k) \times n}$ 的任意 $n-k+1$ 列在 $(n-k)$ 维空间中必线性相关，故 $d \le n-k+1$。

    (3) 纠正 $t$ 个错误的码的球 $B(\mathbf{c}, t)$ 互不相交，每个球包含 $\sum_{i=0}^t \binom{n}{i}(q-1)^i$ 个向量，总数不超过 $q^n$，而有 $q^k$ 个码字，故 $q^k \cdot \sum_{i=0}^t \binom{n}{i}(q-1)^i \le q^n$。$\blacksquare$

!!! definition "定义 30.12 (Hamming 码)"
    对正整数 $r \ge 2$，$[2^r - 1, 2^r - 1 - r, 3]$ **Hamming 码** 的校验矩阵 $H \in \mathbb{F}_2^{r \times (2^r-1)}$ 的列是 $\mathbb{F}_2^r$ 中所有非零向量。Hamming 码可纠 $1$ 位错误，是完美码。

    **伴随式译码**：接收 $\mathbf{r} = \mathbf{c} + \mathbf{e}$，伴随式 $\mathbf{s} = H\mathbf{r}^T = H\mathbf{e}^T$。若 $\mathbf{e} = \mathbf{e}_j$，则 $\mathbf{s} = \mathbf{h}_j$（$H$ 的第 $j$ 列），直接定位错误。

!!! definition "定义 30.13 (Reed-Solomon 码)"
    设 $\alpha_1, \ldots, \alpha_n$ 为 $\mathbb{F}_q$ 中 $n$ 个不同元素。$[n, k, n-k+1]$ **Reed-Solomon 码** 定义为

    $$
    \mathcal{C}_{\text{RS}} = \{(p(\alpha_1), \ldots, p(\alpha_n)) : p \in \mathbb{F}_q[x], \deg(p) < k\}.
    $$

    生成矩阵为 Vandermonde 矩阵。$G$ 的任意 $k$ 列构成 $k \times k$ Vandermonde 矩阵，行列式 $\ne 0$（各点不同），故 RS 码达到 Singleton 界，是 MDS 码。

!!! example "例 30.6"
    **$[7, 4, 3]$ Hamming 码的纠错。** $r = 3$，校验矩阵列为 $1$ 到 $7$ 的二进制表示：

    $$
    H = \begin{pmatrix}0&0&0&1&1&1&1\\0&1&1&0&0&1&1\\1&0&1&0&1&0&1\end{pmatrix}.
    $$

    发送 $\mathbf{c} = (0,1,1,0,0,1,0)$，接收 $\mathbf{r} = (0,1,1,0,1,1,0)$（第 5 位错）。伴随式 $\mathbf{s} = H\mathbf{r}^T = (1,0,1)^T$，$(101)_2 = 5$，定位第 5 位。纠正后得 $\hat{\mathbf{c}} = \mathbf{c}$。

    码率 $R = 4/7 \approx 0.57$。Hamming 界验证：$1 + 7 = 8 = 2^3 = 2^{n-k}$，完美码。

---

## 30.7 LDPC 码与稀疏矩阵

<div class="context-flow" markdown>

**核心**：LDPC 码 = 校验矩阵 $H$ 稀疏 → Tanner 图（二部图）表示 → 置信传播译码利用稀疏结构 → 逼近 Shannon 极限 → 度分布优化

**链接**：30.6 线性码 · Ch17 稀疏矩阵/非负矩阵

</div>

低密度奇偶校验（LDPC）码是现代通信系统中性能最优的码类之一，其核心特征是校验矩阵的稀疏性。

!!! definition "定义 30.14 (LDPC 码)"
    **低密度奇偶校验码**（LDPC code）是 $[n, k]$ 线性码，其校验矩阵 $H \in \mathbb{F}_2^{(n-k) \times n}$ 是**稀疏矩阵**：每行和每列中 $1$ 的个数远小于 $n$。

    - **正则 LDPC 码**：$H$ 的每列恰好有 $d_v$ 个 $1$（变量节点度），每行恰好有 $d_c$ 个 $1$（校验节点度），满足 $(n-k)d_c = n d_v$。
    - **非正则 LDPC 码**：度分布用多项式对 $(\lambda(x), \rho(x))$ 描述，可通过密度进化优化设计。

!!! definition "定义 30.15 (Tanner 图)"
    LDPC 码的 **Tanner 图** 是二部图 $G = (V \cup C, E)$：

    - **变量节点** $V = \{v_1, \ldots, v_n\}$ 对应 $H$ 的 $n$ 列。
    - **校验节点** $C = \{c_1, \ldots, c_{n-k}\}$ 对应 $H$ 的 $n-k$ 行。
    - $v_j$ 与 $c_i$ 相连当且仅当 $H_{ij} = 1$。

    $H$ 的稀疏性意味着 Tanner 图的边数 $|E| = \operatorname{nnz}(H) \ll n(n-k)$。Tanner 图的**围长**（girth，最短环长度）越大，BP 译码性能越好。

!!! theorem "定理 30.10 (置信传播译码)"
    LDPC 码的**置信传播**（BP）译码在 Tanner 图上进行消息传递。设 $L_j$ 为变量节点 $v_j$ 的信道对数似然比（LLR），迭代公式为

    **变量 → 校验消息**：

    $$
    \mu_{v_j \to c_i}^{(\ell)} = L_j + \sum_{c \in \mathcal{N}(v_j) \setminus c_i} \mu_{c \to v_j}^{(\ell-1)},
    $$

    **校验 → 变量消息**：

    $$
    \mu_{c_i \to v_j}^{(\ell)} = 2\operatorname{atanh}\!\left(\prod_{v \in \mathcal{N}(c_i) \setminus v_j} \tanh\!\left(\frac{\mu_{v \to c_i}^{(\ell)}}{2}\right)\right).
    $$

    每次迭代复杂度 $O(|E|) = O(\operatorname{nnz}(H))$，由稀疏性保证高效。当 Tanner 图无环时 BP 给出精确后验概率；有环时为近似算法，但实践中性能优异。

??? proof "证明"
    BP 基于因子图上的和积算法。对 AWGN 信道，LLR $L_j = 2y_j/\sigma^2$。变量到校验的消息汇总来自信道和其他校验的信息。校验到变量的消息基于校验约束 $\bigoplus_{v \in \mathcal{N}(c_i)}c_v = 0$，通过 $\tanh$ 规则在 LLR 域中实现。

    收敛性分析使用**密度进化**（density evolution）：追踪 LLR 消息的概率密度函数在迭代中的演化，可以精确计算给定度分布的 LDPC 码在二进制输入 AWGN 信道上的阈值（接近 Shannon 极限的信噪比）。$\blacksquare$

!!! theorem "定理 30.11 (LDPC 码的最小距离)"
    对正则 $(d_v, d_c)$-LDPC 码：

    1. **下界**：$d \ge d_v + 1$（稀疏性使小集合不易线性相关）。
    2. Tanner 图的围长 $g$ 满足 $d \ge g/2 + 1$（粗略下界）。
    3. 随机构造的 LDPC 码，当 $n \to \infty$ 时，$d = \Theta(n)$ 以高概率成立。

??? proof "证明"
    (1) 设非零码字 $\mathbf{c}$ 的重量为 $w$。$H\mathbf{c}^T = \mathbf{0}$ 要求 $\mathbf{c}$ 支撑对应的 $H$ 的 $w$ 列线性相关。每列 $d_v$ 个 $1$，当 $w \le d_v$ 时，在 Tanner 图无短环条件下不可能线性相关。

    (2) 围长 $g$ 意味着 Tanner 图中最短环涉及至少 $g/2$ 个变量节点，对应最短线性相关列集大小 $\ge g/2$，故 $d \ge g/2 + 1$。

    (3) 随机 LDPC 码的 Tanner 图局部类似树，Gilbert-Varshamov 型论证给出线性增长的最小距离。$\blacksquare$

!!! example "例 30.7"
    **$(3, 6)$-正则 LDPC 码。** 设 $n = 12$，$d_v = 3$，$d_c = 6$，则 $n - k = n \cdot d_v / d_c = 6$，$k = 6$，码率 $R = 1/2$。校验矩阵 $H \in \mathbb{F}_2^{6 \times 12}$，每列 $3$ 个 $1$，每行 $6$ 个 $1$，共 $\operatorname{nnz}(H) = 36$ 个非零元素。

    Tanner 图有 $12$ 个变量节点和 $6$ 个校验节点，$36$ 条边。与稠密 $6 \times 12$ 矩阵的 $72$ 个元素相比，稀疏度 $= 36/72 = 50\%$。对于实际的 LDPC 码（如 5G NR 中的 $n \sim 10^4$），稀疏度通常 $< 1\%$，BP 译码每次迭代极为高效。

    在 AWGN 信道上，经过 $50$--$100$ 次 BP 迭代，性能可达 Shannon 极限的 $0.1$ dB 以内。

---

## 30.8 小波变换

<div class="context-flow" markdown>

**核心**：小波变换 = 多分辨率分析（MRA）→ 尺度函数 $\phi$ 和小波函数 $\psi$ 的正交性 → 离散小波变换 = 滤波器组（低通 $h$ + 高通 $g$）→ 矩阵分解为稀疏正交矩阵的乘积 → Haar 小波 = 最简单的例子

**链接**：Ch8 正交性 · 30.1 DFT(对比) · Ch11 矩阵分解

</div>

小波变换提供了信号的时频联合表示，弥补了 Fourier 变换仅有频率分辨率而无时间分辨率的不足。

!!! definition "定义 30.16 (多分辨率分析)"
    $L^2(\mathbb{R})$ 上的**多分辨率分析**（MRA）是一族嵌套闭子空间 $\cdots \subset V_{-1} \subset V_0 \subset V_1 \subset \cdots$，满足：

    1. $\overline{\bigcup_j V_j} = L^2(\mathbb{R})$，$\bigcap_j V_j = \{0\}$。
    2. $f(t) \in V_j \iff f(2t) \in V_{j+1}$（尺度关系）。
    3. 存在**尺度函数** $\phi \in V_0$，使得 $\{\phi(t - k)\}_{k \in \mathbb{Z}}$ 构成 $V_0$ 的正交基。

    **小波空间** $W_j$ 为 $V_j$ 在 $V_{j+1}$ 中的正交补：$V_{j+1} = V_j \oplus W_j$。**小波函数** $\psi$ 满足 $\{\psi(t - k)\}_{k \in \mathbb{Z}}$ 是 $W_0$ 的正交基。

!!! definition "定义 30.17 (离散小波变换)"
    **离散小波变换**（DWT）通过滤波器组实现。设尺度函数满足**两尺度关系**

    $$
    \phi(t) = \sqrt{2}\sum_k h_k \phi(2t - k),
    $$

    小波函数 $\psi(t) = \sqrt{2}\sum_k g_k \phi(2t - k)$，其中 $g_k = (-1)^k h_{1-k}$。DWT 的一层分解将信号 $\mathbf{x}$ 分为**近似系数**（低频）$\mathbf{a} = H_{\text{low}}\mathbf{x}$ 和**细节系数**（高频）$\mathbf{d} = H_{\text{high}}\mathbf{x}$，其中 $H_{\text{low}}$ 和 $H_{\text{high}}$ 分别对应低通和高通滤波加下采样。

    $J$ 层 DWT 对应矩阵分解 $W = W_J W_{J-1} \cdots W_1$，每个 $W_j$ 是稀疏正交矩阵，总复杂度 $O(n)$（远优于 FFT 的 $O(n\log n)$）。

!!! theorem "定理 30.12 (DWT 的正交性与完美重构)"
    设 $\{h_k\}$ 和 $\{g_k\}$ 满足正交条件

    $$
    \sum_k h_k h_{k-2m} = \delta_{m0}, \quad \sum_k g_k g_{k-2m} = \delta_{m0}, \quad \sum_k h_k g_{k-2m} = 0,
    $$

    则 DWT 矩阵 $W$ 为正交矩阵（$W^T W = I$），信号可通过 $\mathbf{x} = W^T \mathbf{c}$ 完美重构（$\mathbf{c}$ 为小波系数）。

    等价地，正交 MRA 的分解和重构形成完美重构滤波器组，无信息损失。

??? proof "证明"
    $W$ 的每一层 $W_j$ 由低通和高通滤波器加下采样构成。正交条件保证 $W_j^T W_j = I$（每层正交），因此 $W = W_J \cdots W_1$ 也正交。$W^T W = W_1^T \cdots W_J^T W_J \cdots W_1 = I$。

    完美重构：$\mathbf{x} = W^T(W\mathbf{x}) = W^T\mathbf{c}$，对应逐层上采样加滤波的重构过程。$\blacksquare$

!!! theorem "定理 30.13 (Haar 小波)"
    **Haar 小波** 是最简单的正交小波，滤波器系数为 $h_0 = h_1 = 1/\sqrt{2}$，$g_0 = 1/\sqrt{2}$，$g_1 = -1/\sqrt{2}$。对 $n = 2^J$ 长度的信号，Haar DWT 矩阵的一层分解为

    $$
    W_1 = \frac{1}{\sqrt{2}}\begin{pmatrix}1&1&0&0&\cdots\\0&0&1&1&\cdots\\\hline 1&-1&0&0&\cdots\\0&0&1&-1&\cdots\end{pmatrix},
    $$

    上半部分为低通（求和/平均），下半部分为高通（差分）。Haar 变换等价于递归地对相邻元素求平均和差值。

??? proof "证明"
    Haar 尺度函数 $\phi = \mathbf{1}_{[0,1)}$，小波函数 $\psi = \mathbf{1}_{[0,1/2)} - \mathbf{1}_{[1/2,1)}$。两尺度关系 $\phi(t) = \phi(2t) + \phi(2t - 1)$（$h_0 = h_1 = 1/\sqrt{2}$）。

    正交性验证：$\langle \phi(\cdot - k), \phi(\cdot - m) \rangle = \delta_{km}$，$\langle \psi(\cdot - k), \psi(\cdot - m) \rangle = \delta_{km}$，$\langle \phi(\cdot - k), \psi(\cdot - m) \rangle = 0$。

    $W_1$ 的行两两正交且范数为 $1$，故 $W_1$ 正交。$J$ 层分解递归地对近似系数应用 $W_1$，得到完整的 Haar DWT。$\blacksquare$

!!! example "例 30.8"
    **Haar 小波分解。** 对信号 $\mathbf{x} = (4, 6, 10, 2, 3, 1, 7, 5)^T$（$n = 8 = 2^3$）进行 $3$ 层 Haar DWT。

    **第 1 层**：低通 $\mathbf{a}^{(1)} = \frac{1}{\sqrt{2}}(10, 12, 4, 12)^T$（相邻求和），高通 $\mathbf{d}^{(1)} = \frac{1}{\sqrt{2}}(-2, 8, 2, 2)^T$（相邻求差）。

    **第 2 层**：对 $\mathbf{a}^{(1)}$ 继续分解。$\mathbf{a}^{(2)} = \frac{1}{2}(22, 16)^T$，$\mathbf{d}^{(2)} = \frac{1}{2}(-2, -8)^T$。

    **第 3 层**：$\mathbf{a}^{(3)} = \frac{1}{2\sqrt{2}}(38)$，$\mathbf{d}^{(3)} = \frac{1}{2\sqrt{2}}(6)$。

    小波系数向量 $\mathbf{c} = (a^{(3)}, d^{(3)}, d_1^{(2)}, d_2^{(2)}, d_1^{(1)}, d_2^{(1)}, d_3^{(1)}, d_4^{(1)})$。

    **稀疏性**：若信号在某些尺度上是平滑的，则对应的细节系数 $d$ 接近零，可以压缩。JPEG 2000 图像压缩即基于此原理。DWT 的计算量仅为 $O(n)$（每层 $O(n_j)$，$\sum n_j = 2n$），远快于 FFT 的 $O(n \log n)$。

!!! note "DFT 与 DWT 的对比"
    | 特征 | DFT/FFT | DWT |
    |------|---------|-----|
    | 基函数 | 正弦/余弦（全局） | 小波（局部化） |
    | 频率分辨率 | 均匀 | 对数（低频细，高频粗） |
    | 时间分辨率 | 无 | 有（多尺度） |
    | 计算复杂度 | $O(n \log n)$ | $O(n)$ |
    | 矩阵结构 | 酉矩阵 $F_n$ | 正交稀疏矩阵乘积 |
    | 典型应用 | 频谱分析、滤波 | 图像压缩、去噪 |

    两者互补：DFT 适合稳态频率分析，DWT 适合非稳态信号的多尺度分析。在压缩感知中，信号可能在 DFT 域或 DWT 域稀疏，选择合适的变换域是关键。

## 练习题

1. **[DFT矩阵] 离散傅里叶变换（DFT）的矩阵 $F_n$ 是什么类型的矩阵？**

   ??? success "参考答案"
       它是对称的复酉矩阵（在乘以归一化常数 $\frac{1}{\sqrt{n}}$ 后），即 $\frac{1}{n} F_n F_n^H = I$。

2. **[基底] DFT 的行向量（或列向量）在 $\mathbb{C}^n$ 空间中构成了什么？**

   ??? success "参考答案"
       构成了一组极其特殊的标准正交基，即离散时间复简谐波基底 $e^{-i \frac{2\pi}{n} k j}$。

3. **[循环矩阵] 循环矩阵具有什么极其深刻的谱性质？**

   ??? success "参考答案"
       所有相同维度的循环矩阵都可以被同一个矩阵——离散傅里叶变换矩阵 $F_n$ 所对角化。这意味着复简谐波是所有平移不变系统的共同特征向量。

4. **[卷积定理] 离散时间下的循环卷积，在变换到 DFT 频域后变成了什么运算？**

   ??? success "参考答案"
       变成了两个向量频谱的逐元素相乘（Hadamard积）。这是“时域卷积等于频域乘积”在离散有限长序列上的严格代数版本。

5. **[FFT] 为什么快速傅里叶变换（FFT）被誉为 20 世纪最伟大的算法之一？**

   ??? success "参考答案"
       它巧妙利用了旋转因子 $W_n$ 的对称性和周期性，采用分治法（Divide and Conquer），将 DFT 原本 $O(n^2)$ 的稠密矩阵乘法暴力计算，降维打击到了 $O(n \log n)$ 级别。

6. **[滤波器] 线性时不变（LTI）系统的离散冲激响应矩阵是什么形状的？**

   ??? success "参考答案"
       是 Toeplitz 矩阵。如果信号受到周期边界条件的约束，则退化为循环矩阵。

7. **[小波] 小波变换与傅里叶变换在基函数的分布支撑上有什么根本不同？**

   ??? success "参考答案"
       傅里叶的基函数是无限延伸的全局正弦波，在时间上完全无法定位（如一个突刺脉冲）；而小波的基函数在时域上是紧支撑的（局部化的），能够同时且精确地捕捉频率成分和时间跳变信息。

8. **[DWT] 离散小波变换（DWT）的矩阵有什么显著的结构特征？**

   ??? success "参考答案"
       它通常是高度稀疏的正交矩阵（对于正交小波），其变换过程表现为在矩阵乘法中不断迭代的低通与高通滤波下采样操作。

9. **[计算复杂度] 与 FFT 的 $O(n \log n)$ 相比，多尺度 DWT 的时间复杂度可以达到多少？**

   ??? success "参考答案"
       由于在每一尺度上不断按因子 2 进行下采样，总计算量是一个收敛的几何级数，其复杂度可以达到令人惊叹的绝对线性 $O(n)$。

10. **[爱因斯坦思考题] 无论是傅里叶变换（将信号投影到纯频域），还是小波变换（投影到时频多尺度空间），它们在线性代数的最底层其实只是做了一件什么微不足道的小事？这揭示了物理观察的什么原理？**

   ??? success "参考答案"
        它们都只是在做**一次极其平凡的“基变换（坐标系旋转）”**：$\mathbf{y} = Q\mathbf{x}$。它深刻地揭示了：那些看似极其复杂、杂乱无章的物理或电子现象，往往只是因为我们固执地站在了错误的“时域基底”上去观察。只要我们通过正交矩阵旋转到一个与该物理系统内在规律完美共振的“本征基底”（如平移不变的频域基底），宇宙隐藏的和谐法则就会如水晶般清澈地展现在对角线上。

## 本章小结

本章证明了信号处理中那些著名的变换与滤波器，本质上都是应用线性代数中的矩阵分解与基变换理论，主要内容包括：

1. **傅里叶变换的矩阵视角**：从复指数向量的正交性出发，将 DFT 严格定义为一个酉矩阵 $F_n$ 乘法，它是联系时域序列和频域序列的等距变换（保持信号总能量不变，Parseval定理）。
2. **快速傅里叶变换 (FFT)**：通过对 DFT 矩阵的分块分解（基于旋转因子的代数对称性），将算法复杂度从平方级降低到了对数线性级。
3. **循环矩阵与卷积定理**：证明了任意循环矩阵都能被 DFT 矩阵对角化，从而在代数上极其优美地解释了为什么“时域的循环卷积等价于频域的逐元素相乘”。
4. **线性时不变系统 (LTI)**：将信号的滤波操作抽象为与 Toeplitz 矩阵（或循环矩阵）的左乘，为现代通信系统打下了代数基础。
5. **小波变换 (DWT)**：引入了多分辨率分析框架（MRA），利用高通/低通滤波器的正交矩阵树，实现了比傅里叶更加强大的时-频局部化信号压缩，且算法复杂度极低至 $O(n)$。
