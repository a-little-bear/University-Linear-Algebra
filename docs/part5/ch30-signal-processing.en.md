# Chapter 30  Linear Algebra in Signal Processing and Coding

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · matrix decompositions (Ch10-11) · positive definite matrices (Ch16) · Kronecker products (Ch19) · **Arc**: DFT matrix (unitary transform) → FFT (sparse matrix factorization) → circulant matrices and convolution (DFT diagonalization) → filtering and frequency domain analysis (frequency response) → sampling theorem and interpolation (bandlimited signals) → compressed sensing (sparse recovery / RIP) → linear error-correcting codes (linear algebra over finite fields) → LDPC codes (sparse parity-check matrices) → wavelet transform (multiresolution analysis)
**Essence**: The core operations of signal processing — transforms, filtering, convolution — are all matrix-vector multiplications; the core concepts of coding theory — codewords, parity checks, error correction — are all linear algebra over finite fields; compressed sensing and wavelet transforms open new applications of linear algebra in sparsity and multiscale analysis

</div>

Signal processing and coding theory are among the most classical engineering application domains of linear algebra. The discrete Fourier transform is essentially a unitary matrix multiplication, circular convolution corresponds to a circulant matrix, filtering is a diagonal matrix operation in the frequency domain, and sampling and interpolation are precisely characterized by the Nyquist-Shannon theorem. On the coding side, the design and decoding of error-correcting codes are entirely built on linear algebra over finite fields; LDPC codes exploit sparse matrix structure to achieve performance approaching the Shannon limit. This chapter concludes with compressed sensing and wavelet transforms, showcasing the frontiers of linear algebra in modern signal processing.

---

## 30.1 The Discrete Fourier Transform

<div class="context-flow" markdown>

**Linear algebra perspective**: DFT = multiplication by the unitary matrix $F_n$ ($F_n^* F_n = nI$) → eigenvalues/eigenvectors completely known → DFT diagonalizes circular convolution → FFT decomposes $F_n$ into a product of $\log_2 n$ sparse matrices
**Links**: Ch8 unitary matrices · Ch6 eigenvalues

</div>

The discrete Fourier transform (DFT) is the most fundamental tool in signal processing, and its essence is a special matrix-vector multiplication.

!!! definition "Definition 30.1 (DFT Matrix)"
    Let $\omega_n = e^{-2\pi i/n}$ be a primitive $n$-th root of unity. The $n \times n$ **DFT matrix** $F_n \in \mathbb{C}^{n \times n}$ is defined as

    $$
    (F_n)_{jk} = \omega_n^{jk}, \quad j, k = 0, 1, \ldots, n-1.
    $$

    That is,

    $$
    F_n = \begin{pmatrix}1&1&1&\cdots&1\\1&\omega_n&\omega_n^2&\cdots&\omega_n^{n-1}\\1&\omega_n^2&\omega_n^4&\cdots&\omega_n^{2(n-1)}\\\vdots&\vdots&\vdots&\ddots&\vdots\\1&\omega_n^{n-1}&\omega_n^{2(n-1)}&\cdots&\omega_n^{(n-1)^2}\end{pmatrix}.
    $$

    The **DFT** of a vector $\mathbf{x} \in \mathbb{C}^n$ is $\hat{\mathbf{x}} = F_n \mathbf{x}$, i.e., $\hat{x}_j = \sum_{k=0}^{n-1} x_k \omega_n^{jk}$.

!!! definition "Definition 30.2 (Inverse DFT and FFT)"
    The **inverse DFT** is $\mathbf{x} = \frac{1}{n}F_n^* \hat{\mathbf{x}}$, where $F_n^* = \overline{F_n}^T$.

    When $n = 2^s$, the **fast Fourier transform** (FFT) uses the Cooley-Tukey decomposition

    $$
    F_n = \begin{pmatrix}I_m & D_m \\ I_m & -D_m\end{pmatrix}\begin{pmatrix}F_m & 0 \\ 0 & F_m\end{pmatrix}P_n, \quad m = n/2,
    $$

    to decompose the $n$-point DFT into two $n/2$-point DFTs plus twiddle factor multiplications. Here $D_m = \operatorname{diag}(1, \omega_n, \ldots, \omega_n^{m-1})$ and $P_n$ is the perfect shuffle permutation. Recursion yields total complexity $O(n \log n)$.

!!! theorem "Theorem 30.1 (Unitarity of the DFT Matrix)"
    The normalized DFT matrix $\frac{1}{\sqrt{n}}F_n$ is unitary, i.e.,

    $$
    F_n^* F_n = F_n F_n^* = nI_n.
    $$

    Equivalently, the columns (or rows) of $F_n$ form an orthogonal basis for $\mathbb{C}^n$ (with inner product $n$), and the columns of $\frac{1}{\sqrt{n}}F_n$ form an orthonormal basis.

??? proof "Proof"
    $(F_n^* F_n)_{jk} = \sum_{\ell=0}^{n-1}\overline{(F_n)_{\ell j}}(F_n)_{\ell k} = \sum_{\ell=0}^{n-1}\omega_n^{-\ell j}\omega_n^{\ell k} = \sum_{\ell=0}^{n-1}\omega_n^{\ell(k-j)}$. When $k = j$, the sum is $n$. When $k \ne j$, this is a geometric series $\sum_{\ell=0}^{n-1}(\omega_n^{k-j})^\ell = \frac{1 - \omega_n^{n(k-j)}}{1 - \omega_n^{k-j}} = 0$, since $\omega_n^n = 1$ and $\omega_n^{k-j} \ne 1$ ($0 < |k-j| < n$). $\blacksquare$

!!! theorem "Theorem 30.2 (Eigenvalues of the DFT Matrix and FFT Complexity)"
    1. $F_n$ satisfies $F_n^4 = n^2 I_n$, so its eigenvalues are $\{+\sqrt{n}, -\sqrt{n}, +i\sqrt{n}, -i\sqrt{n}\}$.
    2. When $n = 2^s$, the FFT requires $\frac{n}{2}\log_2 n$ complex multiplications and $n\log_2 n$ complex additions, giving total complexity $O(n\log n)$. The FFT decomposes $F_n$ into a product of $\log_2 n$ sparse matrices, each with $O(n)$ nonzero entries.

??? proof "Proof"
    (1) $(F_n^2)_{jk} = \sum_\ell \omega_n^{\ell(j+k)}$. The sum is $n$ when $j + k \equiv 0 \pmod{n}$ and $0$ otherwise, so $F_n^2 = nP$ ($P$ is the reversal permutation matrix). $F_n^4 = n^2 P^2 = n^2 I$, and eigenvalues $\lambda$ satisfy $\lambda^4 = n^2$.

    (2) Let $T(n)$ be the multiplication count. The Cooley-Tukey decomposition gives $T(n) = 2T(n/2) + n/2$, $T(1) = 0$, with solution $T(n) = \frac{n}{2}\log_2 n$. The matrix factorization is $F_n = B_s B_{s-1} \cdots B_1 P$, where each $B_k$ has $O(n)$ nonzero entries. $\blacksquare$

!!! example "Example 30.1"
    **The $4 \times 4$ DFT matrix and FFT.** $\omega_4 = e^{-\pi i/2} = -i$.

    $$
    F_4 = \begin{pmatrix}1&1&1&1\\1&-i&-1&i\\1&-1&1&-1\\1&i&-1&-i\end{pmatrix}.
    $$

    Verification: $F_4^* F_4 = 4I_4$. Applying the DFT to $\mathbf{x} = (1, 2, 3, 4)^T$:

    $$
    \hat{\mathbf{x}} = F_4\mathbf{x} = \begin{pmatrix}10\\-2+2i\\-2\\-2-2i\end{pmatrix}.
    $$

    The FFT decomposes this into 2 stages of butterfly operations, requiring only $\frac{4}{2}\log_2 4 = 4$ multiplications (versus $16$ for direct computation). For $n = 2^{20}$, the FFT speedup ratio is approximately $10^5$.

---

## 30.2 Circulant Matrices and Convolution

<div class="context-flow" markdown>

**Core**: Circulant matrix = diagonalized by the DFT matrix $C = \frac{1}{n}F_n^*\operatorname{diag}(\hat{\mathbf{c}})F_n$ → circular convolution = pointwise multiplication in the frequency domain → convolution theorem → "FFT → pointwise multiply → IFFT" $O(n\log n)$ speedup
**Links**: Ch6 eigenvalues · 30.1 DFT

</div>

Circulant matrices are a class of matrices with special structure, deeply connected to the DFT and convolution.

!!! definition "Definition 30.3 (Circulant Matrix)"
    The $n \times n$ **circulant matrix** generated by the vector $\mathbf{c} = (c_0, c_1, \ldots, c_{n-1})^T$ is

    $$
    C = \operatorname{circ}(\mathbf{c}) = \begin{pmatrix}c_0&c_{n-1}&c_{n-2}&\cdots&c_1\\c_1&c_0&c_{n-1}&\cdots&c_2\\c_2&c_1&c_0&\cdots&c_3\\\vdots&\vdots&\vdots&\ddots&\vdots\\c_{n-1}&c_{n-2}&c_{n-3}&\cdots&c_0\end{pmatrix}.
    $$

    Each row is a cyclic right shift of the previous row. Equivalently, $C = \sum_{k=0}^{n-1}c_k P^k$, where $P$ is the cyclic permutation matrix.

!!! definition "Definition 30.4 (Circular Convolution)"
    The **circular convolution** of two length-$n$ vectors $\mathbf{a}$ and $\mathbf{b}$ is $\mathbf{c} = \mathbf{a} \circledast \mathbf{b}$, defined as

    $$
    c_j = \sum_{k=0}^{n-1} a_k b_{(j-k) \bmod n}, \quad j = 0, 1, \ldots, n-1.
    $$

    Circular convolution is equivalent to circulant matrix multiplication: $\mathbf{c} = \operatorname{circ}(\mathbf{a})\mathbf{b}$.

!!! theorem "Theorem 30.3 (Spectral Decomposition of Circulant Matrices and Convolution Theorem)"
    1. **Spectral decomposition**: All $n \times n$ circulant matrices are diagonalized by the same DFT matrix:

        $$
        C = \operatorname{circ}(\mathbf{c}) = \frac{1}{n}F_n^* \operatorname{diag}(\hat{\mathbf{c}}) F_n,
        $$

        where $\hat{\mathbf{c}} = F_n \mathbf{c}$. The eigenvalues of $C$ are $\hat{c}_0, \hat{c}_1, \ldots, \hat{c}_{n-1}$, with eigenvectors being the columns of $F_n$.

    2. **Convolution theorem**: Circular convolution is equivalent to pointwise multiplication in the frequency domain: $\mathbf{a} \circledast \mathbf{b} = \mathbf{c} \iff \hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}}$, where $\odot$ denotes elementwise multiplication.

    Therefore circular convolution can be computed via "FFT → pointwise multiplication → IFFT" in $O(n\log n)$ time. All circulant matrices commute with each other.

??? proof "Proof"
    (1) The cyclic permutation matrix $P$ has eigenvalues $1, \omega_n, \omega_n^2, \ldots, \omega_n^{n-1}$ with eigenvectors being the columns of the DFT matrix. Since $C = \sum_k c_k P^k$, $C$ shares eigenvectors with $P$, and its eigenvalues are $\lambda_j = \sum_k c_k \omega_n^{jk} = \hat{c}_j$. Two circulant matrices are both diagonalized by $F_n$, hence they commute.

    (2) $\mathbf{c} = \operatorname{circ}(\mathbf{a})\mathbf{b} = \frac{1}{n}F_n^*\operatorname{diag}(\hat{\mathbf{a}})F_n\mathbf{b} = \frac{1}{n}F_n^*(\hat{\mathbf{a}} \odot \hat{\mathbf{b}})$. Taking the DFT: $\hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}}$. $\blacksquare$

!!! example "Example 30.2"
    **Computing circular convolution via FFT.** Let $\mathbf{a} = (1, 2, 3, 0)^T$, $\mathbf{b} = (1, 0, 1, 0)^T$.

    DFT: $\hat{\mathbf{a}} = F_4\mathbf{a} = (6, -2+2i, -2, -2-2i)^T$, $\hat{\mathbf{b}} = F_4\mathbf{b} = (2, 0, 2, 0)^T$.

    Frequency domain multiplication: $\hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}} = (12, 0, -4, 0)^T$.

    IFFT: $\mathbf{c} = \frac{1}{4}F_4^*\hat{\mathbf{c}} = (2, 2, 4, 4)^T$.

    Verification via circulant matrix: $\operatorname{circ}(\mathbf{a})\mathbf{b} = \begin{pmatrix}1&0&3&2\\2&1&0&3\\3&2&1&0\\0&3&2&1\end{pmatrix}\begin{pmatrix}1\\0\\1\\0\end{pmatrix} = \begin{pmatrix}4\\2\\4\\2\end{pmatrix}$.

---

## 30.3 Filtering and Frequency Domain Analysis

<div class="context-flow" markdown>

**Core**: Linear time-invariant (LTI) system = Toeplitz matrix multiplication → frequency response $H(\omega)$ = transfer function → FIR filter = finite convolution → IIR filter = recursive difference equation → circulant embedding of Toeplitz matrices and FFT acceleration
**Links**: 30.2 circulant matrices · Ch22 numerical linear algebra

</div>

Filtering is the core operation of signal processing, and its mathematical essence is matrix-vector multiplication.

!!! definition "Definition 30.5 (Toeplitz Matrix and LTI Systems)"
    An $n \times n$ **Toeplitz matrix** $T \in \mathbb{R}^{n \times n}$ satisfies $T_{ij} = t_{i-j}$, meaning the entries are constant along each diagonal:

    $$
    T = \begin{pmatrix}t_0&t_{-1}&t_{-2}&\cdots&t_{-(n-1)}\\t_1&t_0&t_{-1}&\cdots&t_{-(n-2)}\\t_2&t_1&t_0&\cdots&t_{-(n-3)}\\\vdots&\vdots&\vdots&\ddots&\vdots\\t_{n-1}&t_{n-2}&t_{n-3}&\cdots&t_0\end{pmatrix}.
    $$

    A Toeplitz matrix is the matrix representation of a **linear time-invariant** (LTI) system: the output $y_j = \sum_k t_{j-k} x_k$ is the convolution of the input $\mathbf{x}$ with the impulse response $\{t_k\}$.

!!! definition "Definition 30.6 (Frequency Response and Filters)"
    The **frequency response** of an LTI system is $H(\omega) = \sum_{k} h_k e^{-i\omega k}$, the DTFT of the impulse response.

    - **FIR filter** (finite impulse response): $h_k$ is nonzero for only finitely many $k$; $y = H * x$ is a finite convolution. The matrix representation is a banded Toeplitz matrix.
    - **IIR filter** (infinite impulse response): defined by the recursive difference equation $\sum_{k=0}^{N} a_k y_{n-k} = \sum_{k=0}^{M} b_k x_{n-k}$. The transfer function is a rational function $H(z) = B(z)/A(z)$; the stability condition requires all roots of $A(z)$ to lie inside the unit circle.

!!! theorem "Theorem 30.4 (Circulant Embedding of Toeplitz Matrices)"
    Any $n \times n$ Toeplitz matrix $T$ can be embedded into an $m \times m$ circulant matrix $C$ ($m \ge 2n - 1$), so that $T$ is the upper-left $n \times n$ submatrix of $C$. Therefore, the Toeplitz matrix-vector product $T\mathbf{x}$ can be computed via FFT in $O(n \log n)$ time:

    1. Zero-pad $\mathbf{x}$ to length $m$.
    2. Compute the FFTs of $\mathbf{c}$ (the first column of $C$) and the padded $\mathbf{x}$.
    3. Pointwise multiplication in the frequency domain.
    4. Compute the IFFT and extract the first $n$ components.

??? proof "Proof"
    Construct an $m = 2n$ circulant matrix $C = \operatorname{circ}(c_0, c_1, \ldots, c_{2n-1})$, where $c_k = t_k$ ($0 \le k \le n-1$) and $c_k = t_{k-2n}$ ($n \le k \le 2n-1$). Pad $\mathbf{x}$ to $\tilde{\mathbf{x}} = (x_0, \ldots, x_{n-1}, 0, \ldots, 0)^T$. Then $(C\tilde{\mathbf{x}})_j = \sum_{k=0}^{n-1}c_{(j-k)\bmod 2n}x_k$. For $0 \le j \le n-1$, $(j-k) \bmod 2n$ yields exactly $t_{j-k}$, so the first $n$ components equal $T\mathbf{x}$. $\blacksquare$

!!! theorem "Theorem 30.5 (Levinson-Durbin Recursion)"
    For a positive definite symmetric Toeplitz system $T_n \mathbf{a} = \mathbf{b}$, the Levinson-Durbin algorithm exploits the centrosymmetry of the Toeplitz structure, recursively extending the solution from a $1 \times 1$ system, and solves in $O(n^2)$ time (versus $O(n^3)$ for general methods). The key recursion uses the reflection coefficient

    $$
    \kappa_{k+1} = \frac{b_{k+1} - \mathbf{r}_k^T \mathbf{a}^{(k)}}{\epsilon_k},
    $$

    with the update $\mathbf{a}^{(k+1)} = \begin{pmatrix}\mathbf{a}^{(k)}\\ 0\end{pmatrix} + \kappa_{k+1}\begin{pmatrix}J_k \mathbf{a}^{(k)}\\ 1\end{pmatrix}$, where $J_k$ is the reversal matrix and $\epsilon_k$ is the prediction error power.

??? proof "Proof"
    The key is the centrosymmetry of Toeplitz matrices: if $T_k\mathbf{a} = \mathbf{b}$, then $T_k(J_k\mathbf{a}) = J_k\mathbf{b}$. Extending from a $k$-th order solution to a $(k+1)$-th order solution requires determining only one new parameter $\kappa_{k+1}$. This reduces the $O(k)$ new equations to a single-parameter problem, giving total complexity $\sum_{k=1}^n O(k) = O(n^2)$. Stability is ensured by $|\kappa_k| < 1$. $\blacksquare$

!!! example "Example 30.3"
    **Matrix representation of an FIR lowpass filter.** A length-3 moving average filter $h = \frac{1}{3}(1, 1, 1)$ applied to a signal $\mathbf{x}$ of length $6$:

    $$
    T = \frac{1}{3}\begin{pmatrix}1&0&0&0&0&0\\1&1&0&0&0&0\\1&1&1&0&0&0\\0&1&1&1&0&0\\0&0&1&1&1&0\\0&0&0&1&1&1\end{pmatrix}.
    $$

    This is a lower triangular banded Toeplitz matrix. The frequency response is $H(\omega) = \frac{1}{3}(1 + e^{-i\omega} + e^{-2i\omega}) = \frac{1}{3}e^{-i\omega}\frac{\sin(3\omega/2)}{\sin(\omega/2)}$: $|H| = 1$ at $\omega = 0$ (DC passthrough) and $|H|$ decays at high frequencies (lowpass characteristic).

    Using circulant embedding, the computation of $T\mathbf{x}$ can be converted to an FFT-accelerated circular convolution.

---

## 30.4 Sampling Theorem and Interpolation

<div class="context-flow" markdown>

**Core**: Shannon-Nyquist theorem = a bandlimited function is uniquely determined by its samples → sampling matrix and sinc interpolation → aliasing = spectral folding caused by undersampling → linear algebra formulation of interpolation
**Links**: 30.1 DFT · Ch11 SVD/pseudoinverse

</div>

The sampling theorem bridges continuous and discrete signals, and its mathematical formulation is essentially a basis expansion problem in linear algebra.

!!! definition "Definition 30.7 (Bandlimited Signals and Sampling)"
    A **bandlimited signal** $f(t)$ has Fourier transform $\hat{f}(\omega)$ satisfying $\hat{f}(\omega) = 0$ for $|\omega| > W$, with bandwidth $W$. Sampling at interval $T_s = 1/(2W)$ yields the discrete sequence $f_n = f(nT_s)$.

    The **sampling operator** $S : f \mapsto (f(nT_s))_{n \in \mathbb{Z}}$ is a linear map from the space of continuous functions to the space of discrete sequences.

!!! theorem "Theorem 30.6 (Shannon-Nyquist Sampling Theorem)"
    If $f(t)$ is a bandlimited signal with bandwidth $W$ ($\hat{f}(\omega) = 0$ for $|\omega| > W$), then $f$ is uniquely determined by its samples $f_n = f(n/(2W))$ and can be exactly recovered via **sinc interpolation**:

    $$
    f(t) = \sum_{n=-\infty}^{\infty} f_n \operatorname{sinc}(2Wt - n),
    $$

    where $\operatorname{sinc}(x) = \frac{\sin(\pi x)}{\pi x}$. The sampling frequency $f_s = 2W$ is called the **Nyquist rate**.

    When the sampling frequency $f_s < 2W$, **aliasing** occurs: high-frequency components fold into low frequencies, and the signal cannot be exactly recovered.

??? proof "Proof"
    In the frequency domain, sampling periodizes $\hat{f}(\omega)$ along the frequency axis with period $f_s = 2W$: $\hat{f}_s(\omega) = \frac{1}{T_s}\sum_k \hat{f}(\omega - kf_s)$. When $f_s \ge 2W$, the periodized copies do not overlap, and $\hat{f}(\omega)$ can be recovered through an ideal lowpass filter (passing $|\omega| \le W$, blocking otherwise). The time-domain form of the ideal lowpass filter is the $\operatorname{sinc}$ function, yielding the sinc interpolation formula. When $f_s < 2W$, the copies overlap, producing aliasing. $\blacksquare$

!!! definition "Definition 30.8 (Discrete Sinc Interpolation Matrix)"
    For a signal of finite length $N$, the interpolation matrix $S \in \mathbb{R}^{N \times M}$ from $M$ sample points to $N$ points is

    $$
    S_{jk} = \operatorname{sinc}\!\left(\frac{j - k \cdot N/M}{1}\right), \quad j = 0, \ldots, N-1, \quad k = 0, \ldots, M-1.
    $$

    When $M \ge N$ (oversampling), $S$ has more columns than rows and the system is overdetermined; least squares (pseudoinverse $S^+$) can be used. When $M < N$ (undersampling), $S$ has fewer columns than rows and the system is underdetermined.

!!! example "Example 30.4"
    **Recovering an 8-point signal from 4 samples.** Suppose a bandlimited signal has nonzero DFT only in the low-frequency components $\hat{x}_0, \hat{x}_1, \hat{x}_6, \hat{x}_7$ (bandwidth $W = 2$, Nyquist rate $= 4$). Sampling 4 points yields $\mathbf{y} = (y_0, y_1, y_2, y_3)^T$.

    Recovery: compute the 4-point DFT of $\mathbf{y}$ to get $\hat{\mathbf{y}}$; zero-pad to an 8-point spectrum $\hat{\mathbf{x}} = (\hat{y}_0, \hat{y}_1, 0, 0, 0, 0, \hat{y}_2, \hat{y}_3)^T$; apply the 8-point IFFT to obtain the recovered signal.

    Linear algebra perspective: this involves a downsampling matrix $D \in \mathbb{R}^{4 \times 8}$ (selecting every other row of the identity) combined with the DFT matrix: $\mathbf{y} = D\mathbf{x}$. Recovery exploits the bandlimited structure of $\mathbf{x}$ (i.e., sparsity in the DFT domain).

---

## 30.5 Compressed Sensing

<div class="context-flow" markdown>

**Core**: Sparse signals can be exactly recovered from far fewer measurements than the Nyquist rate → measurement matrix $A$ ($m \ll n$) → Restricted Isometry Property (RIP) → $\ell_1$ minimization (basis pursuit) equivalent to $\ell_0$ minimization → recovery conditions linked to geometric properties of matrices
**Links**: Ch11 SVD · Ch10 matrix norms · Ch25 linear programming ($\ell_1$ minimization)

</div>

Compressed sensing is one of the most influential mathematical discoveries of the early 21st century, breaking through the Shannon-Nyquist sampling barrier by exploiting signal sparsity to achieve sub-Nyquist sampling.

!!! definition "Definition 30.9 (Sparse Signals and the Compressed Sensing Problem)"
    A vector $\mathbf{x} \in \mathbb{R}^n$ is called $s$**-sparse** if $\|\mathbf{x}\|_0 := |\{i : x_i \ne 0\}| \le s$. The **compressed sensing problem**: given a measurement matrix $A \in \mathbb{R}^{m \times n}$ ($m \ll n$) and a measurement vector $\mathbf{y} = A\mathbf{x}$, recover $\mathbf{x}$ knowing that it is $s$-sparse.

    Directly solving the $\ell_0$ minimization $\min \|\mathbf{x}\|_0 \text{ s.t. } A\mathbf{x} = \mathbf{y}$ is NP-hard. The core discovery of compressed sensing is: under appropriate conditions, exact recovery is possible via convex relaxation ($\ell_1$ minimization).

!!! definition "Definition 30.10 (Restricted Isometry Property)"
    A matrix $A \in \mathbb{R}^{m \times n}$ satisfies the $s$-th order **Restricted Isometry Property** (RIP) with constant $\delta_s \in [0, 1)$ if, for all $s$-sparse vectors $\mathbf{x}$,

    $$
    (1 - \delta_s)\|\mathbf{x}\|_2^2 \le \|A\mathbf{x}\|_2^2 \le (1 + \delta_s)\|\mathbf{x}\|_2^2.
    $$

    RIP requires that every submatrix of $A$ formed by $s$ columns approximately preserves norms, i.e., $A$ is approximately isometric in all $s$-sparse directions.

!!! theorem "Theorem 30.7 (Candes-Tao RIP Recovery Theorem)"
    Suppose $A$ satisfies the $2s$-th order RIP with constant $\delta_{2s} < \sqrt{2} - 1$. Then for any $s$-sparse signal $\mathbf{x}$, **basis pursuit**

    $$
    \min_{\mathbf{z}} \|\mathbf{z}\|_1 \quad \text{s.t.} \quad A\mathbf{z} = \mathbf{y}
    $$

    exactly recovers $\mathbf{x}$ from $\mathbf{y} = A\mathbf{x}$. That is, the $\ell_1$ minimizer equals the $\ell_0$ minimizer.

    For noisy observations $\mathbf{y} = A\mathbf{x} + \mathbf{e}$ ($\|\mathbf{e}\|_2 \le \epsilon$), using **LASSO** or **basis pursuit denoising** $\min \|\mathbf{z}\|_1 \text{ s.t. } \|A\mathbf{z} - \mathbf{y}\|_2 \le \epsilon$, the recovery error is proportional to the noise level.

??? proof "Proof"
    Let $\mathbf{h} = \hat{\mathbf{x}} - \mathbf{x}$ ($\hat{\mathbf{x}}$ is the $\ell_1$ minimizer). From $A\mathbf{h} = \mathbf{0}$ and $\|\hat{\mathbf{x}}\|_1 \le \|\mathbf{x}\|_1$ (optimality), decompose $\mathbf{h}$ into $\mathbf{h}_{T_0}$ (components on the support of $\mathbf{x}$) and $\mathbf{h}_{T_0^c}$. The $\ell_1$ optimality yields $\|\mathbf{h}_{T_0^c}\|_1 \le \|\mathbf{h}_{T_0}\|_1$.

    Partition $\mathbf{h}_{T_0^c}$ into blocks $T_1, T_2, \ldots$ of size $s$ in decreasing magnitude order. The RIP provides upper and lower bounds on $\|A\mathbf{h}_{T_0 \cup T_1}\|_2^2$. Combined with $A\mathbf{h} = \mathbf{0}$, this progressively shows $\|\mathbf{h}\|_2 = 0$, i.e., $\hat{\mathbf{x}} = \mathbf{x}$. The condition $\delta_{2s} < \sqrt{2} - 1$ is the critical threshold that closes the argument. $\blacksquare$

!!! theorem "Theorem 30.8 (RIP for Random Measurement Matrices)"
    The following random matrices satisfy the $s$-th order RIP ($\delta_s \le \delta$) with high probability, provided the number of rows $m \ge C \cdot s \ln(n/s) / \delta^2$:

    1. **Gaussian random matrices**: $A_{ij} \sim \mathcal{N}(0, 1/m)$ i.i.d.
    2. **Bernoulli random matrices**: $A_{ij} = \pm 1/\sqrt{m}$ with equal probability.
    3. **Random partial Fourier matrices**: $m$ rows selected uniformly at random from $F_n$, then normalized.

    Key conclusion: $m = O(s \log(n/s))$ measurements suffice to recover an $s$-sparse signal in $\mathbb{R}^n$, far fewer than the $n$ measurements required by Nyquist.

??? proof "Proof"
    For Gaussian matrices, fixing an $s$-sparse direction $\mathbf{x}$, $\|A\mathbf{x}\|_2^2$ is a sum of $m$ independent Gaussian random variables divided by $m$, concentrating around $\|\mathbf{x}\|_2^2$ (Chernoff bound). Taking a union bound over all $s$-sparse directions: $\binom{n}{s}$ choices of support, each covered by an $\epsilon$-net of the unit sphere requiring $(9/\delta)^s\binom{n}{s}$ points total. When $m \ge C s\ln(n/s)/\delta^2$, the union bound ensures deviation $\le \delta$ in all sparse directions. $\blacksquare$

!!! example "Example 30.5"
    **Compressed sensing recovery of a sparse signal.** Let $n = 100$, $s = 5$ ($5$-sparse signal), number of measurements $m = 30$.

    Generate a $30 \times 100$ Gaussian measurement matrix $A$, with $\mathbf{y} = A\mathbf{x}$. Solve the $\ell_1$ minimization (linear program):

    $$
    \min_{\mathbf{z}} \|\mathbf{z}\|_1 \quad \text{s.t.} \quad A\mathbf{z} = \mathbf{y}.
    $$

    Theoretical guarantee: $m = 30 \ge C \cdot 5 \cdot \ln(100/5) \approx 5C \cdot 3 = 15C$. With $C = 2$, $m = 30$ is sufficient. In practice, 30 measurements exactly recover a 5-sparse signal in $\mathbb{R}^{100}$, achieving a compression ratio of $30/100 = 30\%$.

---

## 30.6 Linear Algebra Foundations of Error-Correcting Codes

<div class="context-flow" markdown>

**Core**: Linear code = linear subspace over $\mathbb{F}_q$ → generator matrix $G$ (column space) and parity-check matrix $H$ (null space) → minimum distance = minimum number of linearly dependent columns of $H$ → Hamming codes and Reed-Solomon codes
**Links**: Ch4 vector spaces · Ch1 systems of linear equations

</div>

The fundamental idea of coding theory is to introduce redundancy into transmitted information for error correction. Linear codes build this idea entirely on linear algebra.

!!! definition "Definition 30.11 (Linear Codes, Generator and Parity-Check Matrices)"
    An $[n, k]$ **linear code** $\mathcal{C}$ over $\mathbb{F}_q$ is a $k$-dimensional subspace of $\mathbb{F}_q^n$, containing $q^k$ codewords. The **generator matrix** $G \in \mathbb{F}_q^{k \times n}$ satisfies $\mathcal{C} = \operatorname{rowspace}(G)$; the **parity-check matrix** $H \in \mathbb{F}_q^{(n-k) \times n}$ satisfies $\mathcal{C} = \ker(H)$.

    $G$ and $H$ satisfy $GH^T = 0$. In systematic form: when $G = [I_k \mid P]$, then $H = [-P^T \mid I_{n-k}]$.

!!! theorem "Theorem 30.9 (Minimum Distance, Singleton Bound, and MDS Codes)"
    1. **Minimum distance**: The minimum distance $d$ of an $[n, k, d]$ linear code equals the smallest number of columns of $H$ that are linearly dependent, minus one: any $d-1$ columns of $H$ are linearly independent, but there exist $d$ linearly dependent columns. $\mathcal{C}$ can detect $d-1$ errors and correct $\lfloor(d-1)/2\rfloor$ errors.

    2. **Singleton bound**: $d \le n - k + 1$. Codes achieving this bound are called **MDS codes**.

    3. **Hamming bound**: $\sum_{i=0}^{t}\binom{n}{i}(q-1)^i \le q^{n-k}$, where $t = \lfloor(d-1)/2\rfloor$. Codes achieving this bound are called **perfect codes**.

??? proof "Proof"
    (1) $d = \min_{\mathbf{c} \ne \mathbf{0}, \mathbf{c} \in \mathcal{C}} w(\mathbf{c})$ (minimum Hamming weight). $\mathbf{c} \in \mathcal{C}$ if and only if $H\mathbf{c}^T = \mathbf{0}$, meaning the columns of $H$ corresponding to nonzero positions of $\mathbf{c}$ are linearly dependent.

    (2) $H \in \mathbb{F}_q^{(n-k) \times n}$ has $n$ columns in $(n-k)$-dimensional space, so any $n-k+1$ columns must be linearly dependent, giving $d \le n-k+1$.

    (3) The balls $B(\mathbf{c}, t)$ for different codewords are disjoint, each containing $\sum_{i=0}^t \binom{n}{i}(q-1)^i$ vectors. Since there are $q^k$ codewords and the total space has $q^n$ vectors, $q^k \cdot \sum_{i=0}^t \binom{n}{i}(q-1)^i \le q^n$. $\blacksquare$

!!! definition "Definition 30.12 (Hamming Codes)"
    For a positive integer $r \ge 2$, the $[2^r - 1, 2^r - 1 - r, 3]$ **Hamming code** has a parity-check matrix $H \in \mathbb{F}_2^{r \times (2^r-1)}$ whose columns are all $2^r - 1$ nonzero vectors in $\mathbb{F}_2^r$. Hamming codes can correct $1$-bit errors and are perfect codes.

    **Syndrome decoding**: upon receiving $\mathbf{r} = \mathbf{c} + \mathbf{e}$, the syndrome $\mathbf{s} = H\mathbf{r}^T = H\mathbf{e}^T$. If $\mathbf{e} = \mathbf{e}_j$, then $\mathbf{s} = \mathbf{h}_j$ (the $j$-th column of $H$), directly locating the error.

!!! definition "Definition 30.13 (Reed-Solomon Codes)"
    Let $\alpha_1, \ldots, \alpha_n$ be $n$ distinct elements of $\mathbb{F}_q$. The $[n, k, n-k+1]$ **Reed-Solomon code** is defined as

    $$
    \mathcal{C}_{\text{RS}} = \{(p(\alpha_1), \ldots, p(\alpha_n)) : p \in \mathbb{F}_q[x], \deg(p) < k\}.
    $$

    The generator matrix is a Vandermonde matrix. Any $k$ columns of $G$ form a $k \times k$ Vandermonde matrix with nonzero determinant (since the evaluation points are distinct), so RS codes achieve the Singleton bound and are MDS codes.

!!! example "Example 30.6"
    **Error correction with the $[7, 4, 3]$ Hamming code.** $r = 3$; the columns of the parity-check matrix are the binary representations of $1$ through $7$:

    $$
    H = \begin{pmatrix}0&0&0&1&1&1&1\\0&1&1&0&0&1&1\\1&0&1&0&1&0&1\end{pmatrix}.
    $$

    Transmit $\mathbf{c} = (0,1,1,0,0,1,0)$, receive $\mathbf{r} = (0,1,1,0,1,1,0)$ (error in bit 5). Syndrome $\mathbf{s} = H\mathbf{r}^T = (1,0,1)^T$; $(101)_2 = 5$, locating the error at bit 5. Correction yields $\hat{\mathbf{c}} = \mathbf{c}$.

    Code rate $R = 4/7 \approx 0.57$. Hamming bound verification: $1 + 7 = 8 = 2^3 = 2^{n-k}$, confirming a perfect code.

---

## 30.7 LDPC Codes and Sparse Matrices

<div class="context-flow" markdown>

**Core**: LDPC code = sparse parity-check matrix $H$ → Tanner graph (bipartite graph) representation → belief propagation decoding exploits sparsity → approaches Shannon limit → degree distribution optimization
**Links**: 30.6 linear codes · Ch17 sparse matrices / nonnegative matrices

</div>

Low-Density Parity-Check (LDPC) codes are among the highest-performing code families in modern communications, with sparse parity-check matrices as their defining feature.

!!! definition "Definition 30.14 (LDPC Codes)"
    An **LDPC code** is an $[n, k]$ linear code whose parity-check matrix $H \in \mathbb{F}_2^{(n-k) \times n}$ is **sparse**: the number of $1$s in each row and column is much smaller than $n$.

    - **Regular LDPC code**: each column of $H$ has exactly $d_v$ ones (variable node degree), each row has exactly $d_c$ ones (check node degree), satisfying $(n-k)d_c = n d_v$.
    - **Irregular LDPC code**: the degree distribution is described by polynomial pairs $(\lambda(x), \rho(x))$, which can be optimized via density evolution.

!!! definition "Definition 30.15 (Tanner Graph)"
    The **Tanner graph** of an LDPC code is a bipartite graph $G = (V \cup C, E)$:

    - **Variable nodes** $V = \{v_1, \ldots, v_n\}$ correspond to the $n$ columns of $H$.
    - **Check nodes** $C = \{c_1, \ldots, c_{n-k}\}$ correspond to the $n-k$ rows of $H$.
    - $v_j$ is connected to $c_i$ if and only if $H_{ij} = 1$.

    Sparsity of $H$ means $|E| = \operatorname{nnz}(H) \ll n(n-k)$. The **girth** of the Tanner graph (length of the shortest cycle) -- the larger the girth, the better the BP decoding performance.

!!! theorem "Theorem 30.10 (Belief Propagation Decoding)"
    **Belief propagation** (BP) decoding of LDPC codes performs message passing on the Tanner graph. Let $L_j$ be the channel log-likelihood ratio (LLR) for variable node $v_j$. The iterative formulas are

    **Variable → check messages**:

    $$
    \mu_{v_j \to c_i}^{(\ell)} = L_j + \sum_{c \in \mathcal{N}(v_j) \setminus c_i} \mu_{c \to v_j}^{(\ell-1)},
    $$

    **Check → variable messages**:

    $$
    \mu_{c_i \to v_j}^{(\ell)} = 2\operatorname{atanh}\!\left(\prod_{v \in \mathcal{N}(c_i) \setminus v_j} \tanh\!\left(\frac{\mu_{v \to c_i}^{(\ell)}}{2}\right)\right).
    $$

    Each iteration has complexity $O(|E|) = O(\operatorname{nnz}(H))$, guaranteed efficient by sparsity. When the Tanner graph is cycle-free, BP yields exact posterior probabilities; with cycles it is approximate, but performs excellently in practice.

??? proof "Proof"
    BP is based on the sum-product algorithm on factor graphs. For an AWGN channel, $L_j = 2y_j/\sigma^2$. Variable-to-check messages aggregate information from the channel and other check nodes. Check-to-variable messages are based on the parity constraint $\bigoplus_{v \in \mathcal{N}(c_i)}c_v = 0$, implemented via the $\tanh$ rule in the LLR domain.

    Convergence analysis uses **density evolution**: tracking the probability density function of LLR messages through iterations, one can precisely compute the threshold (SNR approaching the Shannon limit) of an LDPC code with a given degree distribution on a binary-input AWGN channel. $\blacksquare$

!!! theorem "Theorem 30.11 (Minimum Distance of LDPC Codes)"
    For regular $(d_v, d_c)$-LDPC codes:

    1. **Lower bound**: $d \ge d_v + 1$ (sparsity makes small sets unlikely to be linearly dependent).
    2. The girth $g$ of the Tanner graph yields $d \ge g/2 + 1$ (rough lower bound).
    3. Randomly constructed LDPC codes satisfy $d = \Theta(n)$ with high probability as $n \to \infty$.

??? proof "Proof"
    (1) Let $\mathbf{c}$ be a nonzero codeword of weight $w$. $H\mathbf{c}^T = \mathbf{0}$ requires that the $w$ columns of $H$ corresponding to the support of $\mathbf{c}$ are linearly dependent. Each column has $d_v$ ones; when $w \le d_v$, under the assumption that the Tanner graph has no short cycles, linear dependence is impossible.

    (2) Girth $g$ means the shortest cycle in the Tanner graph involves at least $g/2$ variable nodes, corresponding to a minimum linearly dependent column set of size $\ge g/2$, so $d \ge g/2 + 1$.

    (3) The Tanner graph of random LDPC codes is locally tree-like; a Gilbert-Varshamov type argument yields linearly growing minimum distance. $\blacksquare$

!!! example "Example 30.7"
    **A $(3, 6)$-regular LDPC code.** Let $n = 12$, $d_v = 3$, $d_c = 6$. Then $n - k = n \cdot d_v / d_c = 6$, $k = 6$, code rate $R = 1/2$. The parity-check matrix $H \in \mathbb{F}_2^{6 \times 12}$ has $3$ ones per column and $6$ ones per row, totaling $\operatorname{nnz}(H) = 36$ nonzero entries.

    The Tanner graph has $12$ variable nodes and $6$ check nodes, with $36$ edges. Compared to a dense $6 \times 12$ matrix with $72$ entries, the sparsity is $36/72 = 50\%$. For practical LDPC codes (e.g., in 5G NR with $n \sim 10^4$), the sparsity is typically $< 1\%$, making each BP iteration extremely efficient.

    On an AWGN channel, after $50$--$100$ BP iterations, performance can come within $0.1$ dB of the Shannon limit.

---

## 30.8 Wavelet Transform

<div class="context-flow" markdown>

**Core**: Wavelet transform = multiresolution analysis (MRA) → orthogonality of scaling function $\phi$ and wavelet function $\psi$ → discrete wavelet transform = filter bank (lowpass $h$ + highpass $g$) → matrix factorization into products of sparse orthogonal matrices → Haar wavelet = simplest example
**Links**: Ch8 orthogonality · 30.1 DFT (comparison) · Ch11 matrix decompositions

</div>

The wavelet transform provides a joint time-frequency representation of signals, compensating for the Fourier transform's lack of time resolution.

!!! definition "Definition 30.16 (Multiresolution Analysis)"
    A **multiresolution analysis** (MRA) on $L^2(\mathbb{R})$ is a family of nested closed subspaces $\cdots \subset V_{-1} \subset V_0 \subset V_1 \subset \cdots$ satisfying:

    1. $\overline{\bigcup_j V_j} = L^2(\mathbb{R})$ and $\bigcap_j V_j = \{0\}$.
    2. $f(t) \in V_j \iff f(2t) \in V_{j+1}$ (scaling relation).
    3. There exists a **scaling function** $\phi \in V_0$ such that $\{\phi(t - k)\}_{k \in \mathbb{Z}}$ forms an orthogonal basis for $V_0$.

    The **wavelet space** $W_j$ is the orthogonal complement of $V_j$ in $V_{j+1}$: $V_{j+1} = V_j \oplus W_j$. The **wavelet function** $\psi$ satisfies that $\{\psi(t - k)\}_{k \in \mathbb{Z}}$ is an orthogonal basis for $W_0$.

!!! definition "Definition 30.17 (Discrete Wavelet Transform)"
    The **discrete wavelet transform** (DWT) is implemented via filter banks. Given the **two-scale relation** for the scaling function

    $$
    \phi(t) = \sqrt{2}\sum_k h_k \phi(2t - k),
    $$

    the wavelet function is $\psi(t) = \sqrt{2}\sum_k g_k \phi(2t - k)$, where $g_k = (-1)^k h_{1-k}$. One level of DWT decomposition splits the signal $\mathbf{x}$ into **approximation coefficients** (low frequency) $\mathbf{a} = H_{\text{low}}\mathbf{x}$ and **detail coefficients** (high frequency) $\mathbf{d} = H_{\text{high}}\mathbf{x}$, where $H_{\text{low}}$ and $H_{\text{high}}$ correspond to lowpass and highpass filtering followed by downsampling.

    A $J$-level DWT corresponds to the matrix factorization $W = W_J W_{J-1} \cdots W_1$, where each $W_j$ is a sparse orthogonal matrix, with total complexity $O(n)$ (much better than FFT's $O(n\log n)$).

!!! theorem "Theorem 30.12 (Orthogonality and Perfect Reconstruction of DWT)"
    If $\{h_k\}$ and $\{g_k\}$ satisfy the orthogonality conditions

    $$
    \sum_k h_k h_{k-2m} = \delta_{m0}, \quad \sum_k g_k g_{k-2m} = \delta_{m0}, \quad \sum_k h_k g_{k-2m} = 0,
    $$

    then the DWT matrix $W$ is orthogonal ($W^T W = I$), and the signal can be perfectly reconstructed via $\mathbf{x} = W^T \mathbf{c}$ ($\mathbf{c}$ is the wavelet coefficient vector).

    Equivalently, the decomposition and reconstruction of an orthogonal MRA form a perfect reconstruction filter bank with no information loss.

??? proof "Proof"
    Each level $W_j$ of $W$ consists of lowpass and highpass filters followed by downsampling. The orthogonality conditions ensure $W_j^T W_j = I$ (each level is orthogonal), so $W = W_J \cdots W_1$ is also orthogonal: $W^T W = W_1^T \cdots W_J^T W_J \cdots W_1 = I$.

    Perfect reconstruction: $\mathbf{x} = W^T(W\mathbf{x}) = W^T\mathbf{c}$, corresponding to the level-by-level upsampling and filtering reconstruction process. $\blacksquare$

!!! theorem "Theorem 30.13 (Haar Wavelet)"
    The **Haar wavelet** is the simplest orthogonal wavelet, with filter coefficients $h_0 = h_1 = 1/\sqrt{2}$, $g_0 = 1/\sqrt{2}$, $g_1 = -1/\sqrt{2}$. For a signal of length $n = 2^J$, one level of the Haar DWT matrix is

    $$
    W_1 = \frac{1}{\sqrt{2}}\begin{pmatrix}1&1&0&0&\cdots\\0&0&1&1&\cdots\\\hline 1&-1&0&0&\cdots\\0&0&1&-1&\cdots\end{pmatrix},
    $$

    with the upper half performing lowpass (summation/averaging) and the lower half performing highpass (differencing). The Haar transform is equivalent to recursively computing averages and differences of adjacent elements.

??? proof "Proof"
    The Haar scaling function is $\phi = \mathbf{1}_{[0,1)}$; the wavelet function is $\psi = \mathbf{1}_{[0,1/2)} - \mathbf{1}_{[1/2,1)}$. The two-scale relation is $\phi(t) = \phi(2t) + \phi(2t - 1)$ ($h_0 = h_1 = 1/\sqrt{2}$).

    Orthogonality verification: $\langle \phi(\cdot - k), \phi(\cdot - m) \rangle = \delta_{km}$, $\langle \psi(\cdot - k), \psi(\cdot - m) \rangle = \delta_{km}$, $\langle \phi(\cdot - k), \psi(\cdot - m) \rangle = 0$.

    The rows of $W_1$ are pairwise orthogonal with unit norm, so $W_1$ is orthogonal. The $J$-level decomposition recursively applies $W_1$ to the approximation coefficients, yielding the complete Haar DWT. $\blacksquare$

!!! example "Example 30.8"
    **Haar wavelet decomposition.** Apply a 3-level Haar DWT to the signal $\mathbf{x} = (4, 6, 10, 2, 3, 1, 7, 5)^T$ ($n = 8 = 2^3$).

    **Level 1**: Lowpass $\mathbf{a}^{(1)} = \frac{1}{\sqrt{2}}(10, 12, 4, 12)^T$ (pairwise sums), highpass $\mathbf{d}^{(1)} = \frac{1}{\sqrt{2}}(-2, 8, 2, 2)^T$ (pairwise differences).

    **Level 2**: Continue decomposing $\mathbf{a}^{(1)}$. $\mathbf{a}^{(2)} = \frac{1}{2}(22, 16)^T$, $\mathbf{d}^{(2)} = \frac{1}{2}(-2, -8)^T$.

    **Level 3**: $\mathbf{a}^{(3)} = \frac{1}{2\sqrt{2}}(38)$, $\mathbf{d}^{(3)} = \frac{1}{2\sqrt{2}}(6)$.

    Wavelet coefficient vector: $\mathbf{c} = (a^{(3)}, d^{(3)}, d_1^{(2)}, d_2^{(2)}, d_1^{(1)}, d_2^{(1)}, d_3^{(1)}, d_4^{(1)})$.

    **Sparsity**: If the signal is smooth at certain scales, the corresponding detail coefficients $d$ are near zero and can be compressed. JPEG 2000 image compression is based on this principle. The DWT computation requires only $O(n)$ operations (each level costs $O(n_j)$, $\sum n_j = 2n$), much faster than FFT's $O(n \log n)$.

!!! note "Comparison of DFT and DWT"
    | Feature | DFT/FFT | DWT |
    |---------|---------|-----|
    | Basis functions | Sinusoids (global) | Wavelets (localized) |
    | Frequency resolution | Uniform | Logarithmic (fine at low, coarse at high freq.) |
    | Time resolution | None | Yes (multiscale) |
    | Computational complexity | $O(n \log n)$ | $O(n)$ |
    | Matrix structure | Unitary matrix $F_n$ | Product of sparse orthogonal matrices |
    | Typical applications | Spectral analysis, filtering | Image compression, denoising |

    The two are complementary: DFT is suited for stationary frequency analysis, while DWT is suited for multiscale analysis of non-stationary signals. In compressed sensing, the signal may be sparse in either the DFT domain or the DWT domain; choosing the appropriate transform domain is key.
