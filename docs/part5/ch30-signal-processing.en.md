# Chapter 30  Linear Algebra in Signal Processing and Coding

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · matrix decompositions (Ch10) · positive definite matrices (Ch16) · Kronecker product (Ch19) · **Chapter arc**: DFT matrix (unitary) → FFT (recursive factorization of DFT matrix) → circulant matrices (diagonalized by DFT) → Toeplitz matrices (approximately circulant) → linear codes (null space/column space) → Hamming codes (parity-check matrix) → RS codes (Vandermonde structure) → LDPC codes (sparse parity-check matrix)
**Essence**: The core operations of signal processing — transforms, filtering, convolution — are all matrix-vector multiplications; the core concepts of coding theory — codewords, parity checks, error correction — are all linear algebra over finite fields

</div>

Signal processing and coding theory are among the most classical engineering applications of linear algebra. The discrete Fourier transform is essentially a unitary matrix multiplication, circular convolution corresponds to circulant matrices, and the design and decoding of error-correcting codes are built entirely on linear algebra over finite fields. This chapter proceeds from the DFT matrix through FFT, circulant matrices, and Toeplitz matrices to linear codes, Hamming codes, Reed-Solomon codes, and LDPC codes, systematically demonstrating the central role of linear algebra in these two major fields.

---

## 30.1 The Discrete Fourier Transform

<div class="context-flow" markdown>

**Linear algebra perspective**: DFT = multiplication by unitary matrix $F_n$ ($F_n^* F_n = nI$) → eigenvalues/eigenvectors fully known → DFT diagonalizes circular convolution
**Link**: Ch8 unitary matrices · Ch6 eigenvalues

</div>

The Discrete Fourier Transform (DFT) is the most fundamental tool in signal processing, with its essence being a special matrix-vector multiplication.

!!! definition "Definition 30.1 (DFT matrix)"
    Let $\omega_n = e^{-2\pi i/n}$ be a primitive $n$-th root of unity. The $n \times n$ **DFT matrix** $F_n \in \mathbb{C}^{n \times n}$ is defined by

    $$
    (F_n)_{jk} = \omega_n^{jk}, \quad j, k = 0, 1, \ldots, n-1.
    $$

    Explicitly,

    $$
    F_n = \begin{pmatrix}1&1&1&\cdots&1\\1&\omega_n&\omega_n^2&\cdots&\omega_n^{n-1}\\1&\omega_n^2&\omega_n^4&\cdots&\omega_n^{2(n-1)}\\\vdots&\vdots&\vdots&\ddots&\vdots\\1&\omega_n^{n-1}&\omega_n^{2(n-1)}&\cdots&\omega_n^{(n-1)^2}\end{pmatrix}.
    $$

    The **DFT** of a vector $\mathbf{x} \in \mathbb{C}^n$ is $\hat{\mathbf{x}} = F_n \mathbf{x}$, i.e., $\hat{x}_j = \sum_{k=0}^{n-1} x_k \omega_n^{jk}$.

!!! definition "Definition 30.2 (Inverse DFT)"
    The **inverse DFT** is $\mathbf{x} = \frac{1}{n}F_n^* \hat{\mathbf{x}}$, where $F_n^* = \overline{F_n}^T$ is the conjugate transpose of $F_n$, with $(F_n^*)_{jk} = \omega_n^{-jk}$.

!!! theorem "Theorem 30.1 (Unitarity of the DFT matrix)"
    The normalized DFT matrix $\frac{1}{\sqrt{n}}F_n$ is unitary:

    $$
    F_n^* F_n = F_n F_n^* = nI_n.
    $$

    Equivalently, the columns (or rows) of $F_n$ form an orthogonal basis for $\mathbb{C}^n$ (with inner product $n$), and the columns of $\frac{1}{\sqrt{n}}F_n$ form an orthonormal basis.

??? proof "Proof"
    $(F_n^* F_n)_{jk} = \sum_{\ell=0}^{n-1}\overline{(F_n)_{\ell j}}(F_n)_{\ell k} = \sum_{\ell=0}^{n-1}\omega_n^{-\ell j}\omega_n^{\ell k} = \sum_{\ell=0}^{n-1}\omega_n^{\ell(k-j)}$. When $k = j$, the sum is $n$. When $k \ne j$, this is a geometric series $\sum_{\ell=0}^{n-1}(\omega_n^{k-j})^\ell = \frac{1 - \omega_n^{n(k-j)}}{1 - \omega_n^{k-j}} = 0$, since $\omega_n^n = 1$ and $\omega_n^{k-j} \ne 1$ ($0 < |k-j| < n$). $\blacksquare$

!!! theorem "Theorem 30.2 (Eigenvalues of the DFT matrix)"
    The DFT matrix satisfies $F_n^4 = n^2 I_n$, so the eigenvalues of $F_n$ are $\{+\sqrt{n}, -\sqrt{n}, +i\sqrt{n}, -i\sqrt{n}\}$ (with multiplicities depending on $n \bmod 4$). In particular, $F_n^2 = nP$, where $P$ is the reversal permutation matrix $(P\mathbf{x})_j = x_{n-j \bmod n}$.

??? proof "Proof"
    $(F_n^2)_{jk} = \sum_{\ell}\omega_n^{j\ell}\omega_n^{\ell k} = \sum_\ell \omega_n^{\ell(j+k)}$. When $j + k \equiv 0 \pmod{n}$ the sum is $n$, otherwise it is $0$. Thus $(F_n^2)_{jk} = n\delta_{j, n-k \bmod n}$, i.e., $F_n^2 = nP$. Therefore $F_n^4 = n^2 P^2 = n^2 I$ (since $P^2 = I$). The eigenvalues $\lambda$ of $F_n$ satisfy $\lambda^4 = n^2$, giving $\lambda \in \{n^{1/2}, -n^{1/2}, in^{1/2}, -in^{1/2}\}$. $\blacksquare$

!!! example "Example 30.1"
    **The $4 \times 4$ DFT matrix.** $\omega_4 = e^{-\pi i/2} = -i$.

    $$
    F_4 = \begin{pmatrix}1&1&1&1\\1&-i&-1&i\\1&-1&1&-1\\1&i&-1&-i\end{pmatrix}.
    $$

    Verify unitarity: $F_4^* F_4 = 4I_4$. Computing the DFT of $\mathbf{x} = (1, 2, 3, 4)^T$:

    $$
    \hat{\mathbf{x}} = F_4\mathbf{x} = \begin{pmatrix}10\\-2+2i\\-2\\-2-2i\end{pmatrix}.
    $$

    Inverse transform verification: $\frac{1}{4}F_4^*\hat{\mathbf{x}} = (1, 2, 3, 4)^T = \mathbf{x}$.

---

## 30.2 The Fast Fourier Transform

<div class="context-flow" markdown>

**Core idea**: FFT is not a new transform but an efficient factorization of the DFT matrix → $F_n = P_n(I_2 \otimes F_{n/2})T_n \cdot (\text{butterfly factors})$ → $O(n\log n)$ instead of $O(n^2)$
**Link**: Ch19 Kronecker product · sparse matrix factorization

</div>

The Fast Fourier Transform (FFT) is one of the most important algorithms of the 20th century, with its essence being the factorization of the DFT matrix into a product of sparse matrices.

!!! definition "Definition 30.3 (Cooley-Tukey factorization)"
    Let $n = 2m$. Splitting the DFT matrix by even and odd indices, $F_n$ can be expressed as

    $$
    F_n = \begin{pmatrix}I_m & D_m \\ I_m & -D_m\end{pmatrix}\begin{pmatrix}F_m & 0 \\ 0 & F_m\end{pmatrix}P_n,
    $$

    where $D_m = \operatorname{diag}(1, \omega_n, \omega_n^2, \ldots, \omega_n^{m-1})$ is the **twiddle factor** matrix and $P_n$ is the **perfect shuffle** permutation (placing even indices first, odd indices second).

!!! theorem "Theorem 30.3 (Complexity of the FFT)"
    When $n = 2^s$, recursively applying the Cooley-Tukey factorization computes the $n$-point DFT using $\frac{n}{2}\log_2 n$ complex multiplications and $n\log_2 n$ complex additions, for a total complexity of $O(n \log n)$, compared to $O(n^2)$ for direct matrix-vector multiplication.

    The FFT factorizes $F_n$ into a product of $\log_2 n$ sparse matrices, each with exactly $O(n)$ nonzero entries.

??? proof "Proof"
    Let $T(n)$ be the number of multiplications for an $n$-point FFT. The Cooley-Tukey factorization splits one $n$-point DFT into two $n/2$-point DFTs plus $n/2$ twiddle factor multiplications. The recurrence is:

    $$
    T(n) = 2T(n/2) + n/2, \quad T(1) = 0.
    $$

    The solution is $T(n) = \frac{n}{2}\log_2 n$. The addition count is similar: $n$ butterfly additions per stage, $\log_2 n$ stages, totaling $n\log_2 n$.

    In matrix form: $F_n = B_s B_{s-1} \cdots B_1 P$, where each $B_k$ is a block-diagonal butterfly matrix with $O(n)$ nonzero entries. $\blacksquare$

!!! theorem "Theorem 30.4 (Kronecker product representation of FFT)"
    When $n = 2^s$, the FFT can be compactly represented using Kronecker products as

    $$
    F_n = \prod_{j=1}^{s} \left(I_{2^{j-1}} \otimes F_2 \otimes I_{2^{s-j}}\right) \cdot D_j \cdot P_j,
    $$

    where $F_2 = \begin{pmatrix}1&1\\1&-1\end{pmatrix}$ is the $2 \times 2$ DFT matrix (Hadamard matrix), $D_j$ are twiddle factor diagonal matrices, and $P_j$ are appropriate permutation matrices. This representation clearly displays the "divide and conquer" structure of the FFT.

??? proof "Proof"
    By the Cooley-Tukey recursion, stage $j$ decomposes $2^j$-point DFTs into pairs of $2^{j-1}$-point DFTs. At the matrix level, the butterfly operation is $F_2$ acting on length-2 sub-vectors, applied in parallel to all subproblems via the Kronecker product $I_{2^{j-1}} \otimes F_2 \otimes I_{2^{s-j}}$. Twiddle factors and permutations adjust the data arrangement. $\blacksquare$

!!! example "Example 30.2"
    **8-point FFT factorization.** $n = 8 = 2^3$, requiring $\log_2 8 = 3$ stages.

    **Stage 1**: Split 8 points into even indices $(x_0, x_2, x_4, x_6)$ and odd indices $(x_1, x_3, x_5, x_7)$, compute a 4-point DFT on each.

    **Stage 2**: Each 4-point DFT is split into two 2-point DFTs.

    **Stage 3**: The 2-point DFT is the $F_2$ butterfly operation.

    Multiplication count: $\frac{8}{2}\log_2 8 = 12$, versus $8^2 = 64$ for direct DFT. Speedup: $64/12 \approx 5.3\times$. For $n = 2^{20} \approx 10^6$, FFT requires about $10^7$ operations versus $10^{12}$ for the direct method — a speedup of roughly $10^5$.

---

## 30.3 Circulant Matrices and Convolution

<div class="context-flow" markdown>

**Core idea**: Circulant matrix = diagonalized by DFT $C = F_n^{-1}\operatorname{diag}(\hat{\mathbf{c}})F_n$ → circular convolution = pointwise multiplication in frequency domain → convolution theorem
**Link**: Ch6 eigenvalues · 30.1 DFT

</div>

Circulant matrices are a class of specially structured matrices with deep connections to the DFT and convolution.

!!! definition "Definition 30.4 (Circulant matrix)"
    The $n \times n$ **circulant matrix** generated by the vector $\mathbf{c} = (c_0, c_1, \ldots, c_{n-1})^T$ is

    $$
    C = \operatorname{circ}(\mathbf{c}) = \begin{pmatrix}c_0&c_{n-1}&c_{n-2}&\cdots&c_1\\c_1&c_0&c_{n-1}&\cdots&c_2\\c_2&c_1&c_0&\cdots&c_3\\\vdots&\vdots&\vdots&\ddots&\vdots\\c_{n-1}&c_{n-2}&c_{n-3}&\cdots&c_0\end{pmatrix}.
    $$

    Each row is a cyclic right shift of the row above. Equivalently, $C = \sum_{k=0}^{n-1}c_k P^k$, where $P$ is the cyclic permutation matrix.

!!! definition "Definition 30.5 (Circular convolution)"
    The **circular convolution** $\mathbf{c} = \mathbf{a} \circledast \mathbf{b}$ of two length-$n$ vectors $\mathbf{a}$ and $\mathbf{b}$ is defined by

    $$
    c_j = \sum_{k=0}^{n-1} a_k b_{(j-k) \bmod n}, \quad j = 0, 1, \ldots, n-1.
    $$

    Circular convolution is equivalent to circulant matrix multiplication: $\mathbf{c} = \operatorname{circ}(\mathbf{a})\mathbf{b}$.

!!! theorem "Theorem 30.5 (Spectral decomposition of circulant matrices)"
    All $n \times n$ circulant matrices are diagonalized by the same DFT matrix:

    $$
    C = \operatorname{circ}(\mathbf{c}) = \frac{1}{n}F_n^* \operatorname{diag}(\hat{\mathbf{c}}) F_n,
    $$

    where $\hat{\mathbf{c}} = F_n \mathbf{c}$ is the DFT of $\mathbf{c}$. The eigenvalues of $C$ are $\hat{c}_0, \hat{c}_1, \ldots, \hat{c}_{n-1}$, with corresponding eigenvectors being the columns $\mathbf{f}_0, \mathbf{f}_1, \ldots, \mathbf{f}_{n-1}$ of $F_n$.

    All circulant matrices commute pairwise and share the same eigenvectors.

??? proof "Proof"
    The cyclic permutation matrix $P$ has eigenvalues $\omega_n^0 = 1, \omega_n, \omega_n^2, \ldots, \omega_n^{n-1}$ with eigenvectors being the columns of the DFT matrix. Since $C = \sum_k c_k P^k$, $C$ shares eigenvectors with $P$, and its eigenvalues are $\lambda_j = \sum_k c_k \omega_n^{jk} = \hat{c}_j$. Therefore $C = \frac{1}{n}F_n^* \operatorname{diag}(\hat{\mathbf{c}})F_n$. Two circulant matrices $C_1, C_2$ are both diagonalized by $F_n$, so $C_1 C_2 = C_2 C_1$. $\blacksquare$

!!! theorem "Theorem 30.6 (Convolution theorem)"
    Circular convolution is equivalent to pointwise multiplication in the frequency domain:

    $$
    \mathbf{a} \circledast \mathbf{b} = \mathbf{c} \iff \hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}},
    $$

    where $\hat{\mathbf{a}} = F_n\mathbf{a}$, $\hat{\mathbf{b}} = F_n\mathbf{b}$, and $\odot$ denotes elementwise multiplication. Therefore circular convolution can be computed via "FFT → pointwise multiplication → IFFT" in $O(n\log n)$ time.

??? proof "Proof"
    $\mathbf{c} = \operatorname{circ}(\mathbf{a})\mathbf{b} = \frac{1}{n}F_n^*\operatorname{diag}(\hat{\mathbf{a}})F_n\mathbf{b} = \frac{1}{n}F_n^*\operatorname{diag}(\hat{\mathbf{a}})\hat{\mathbf{b}} = \frac{1}{n}F_n^*(\hat{\mathbf{a}} \odot \hat{\mathbf{b}})$. Taking the DFT: $\hat{\mathbf{c}} = F_n\mathbf{c} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}}$. $\blacksquare$

!!! example "Example 30.3"
    **Computing circular convolution via FFT.** Let $\mathbf{a} = (1, 2, 3, 0)^T$, $\mathbf{b} = (1, 0, 1, 0)^T$.

    DFT: $\hat{\mathbf{a}} = F_4\mathbf{a} = (6, -2+2i, -2, -2-2i)^T$, $\hat{\mathbf{b}} = F_4\mathbf{b} = (2, 0, 2, 0)^T$.

    Frequency-domain multiplication: $\hat{\mathbf{c}} = \hat{\mathbf{a}} \odot \hat{\mathbf{b}} = (12, 0, -4, 0)^T$.

    IFFT: $\mathbf{c} = \frac{1}{4}F_4^*\hat{\mathbf{c}} = (2, 2, 4, 4)^T$.

    Direct verification using the circulant matrix:

    $$
    \operatorname{circ}(\mathbf{a})\mathbf{b} = \begin{pmatrix}1&0&3&2\\2&1&0&3\\3&2&1&0\\0&3&2&1\end{pmatrix}\begin{pmatrix}1\\0\\1\\0\end{pmatrix} = \begin{pmatrix}4\\2\\4\\2\end{pmatrix}.
    $$

    Recomputing the DFT confirms the result.

---

## 30.4 Toeplitz Matrices and Signal Processing

<div class="context-flow" markdown>

**Core idea**: Toeplitz matrix = constant along each diagonal → matrix representation of LTI systems → embed in circulant matrix for FFT acceleration → Levinson recursion solves Toeplitz systems in $O(n^2)$
**Link**: Ch22 numerical linear algebra · 30.3 circulant matrices

</div>

Toeplitz matrices are the most common structured matrices in signal processing, corresponding to linear time-invariant (LTI) systems.

!!! definition "Definition 30.6 (Toeplitz matrix)"
    An $n \times n$ **Toeplitz matrix** $T \in \mathbb{R}^{n \times n}$ satisfies $T_{ij} = t_{i-j}$, meaning each diagonal has constant entries:

    $$
    T = \begin{pmatrix}t_0&t_{-1}&t_{-2}&\cdots&t_{-(n-1)}\\t_1&t_0&t_{-1}&\cdots&t_{-(n-2)}\\t_2&t_1&t_0&\cdots&t_{-(n-3)}\\\vdots&\vdots&\vdots&\ddots&\vdots\\t_{n-1}&t_{n-2}&t_{n-3}&\cdots&t_0\end{pmatrix}.
    $$

    $T$ is completely determined by $2n - 1$ parameters $t_{-(n-1)}, \ldots, t_0, \ldots, t_{n-1}$. If $T$ is also symmetric ($t_{-k} = t_k$), it is called a **symmetric Toeplitz matrix**.

!!! theorem "Theorem 30.7 (Circulant embedding of Toeplitz matrices)"
    Any $n \times n$ Toeplitz matrix $T$ can be embedded in an $m \times m$ circulant matrix $C$ ($m \ge 2n - 1$) such that $T$ is the upper-left $n \times n$ submatrix of $C$. Therefore, the Toeplitz matrix-vector product $T\mathbf{x}$ can be computed via FFT in $O(n \log n)$ time:

    1. Zero-pad $\mathbf{x}$ to length $m$.
    2. FFT $\mathbf{c}$ (first column of $C$) and the padded $\mathbf{x}$.
    3. Pointwise multiplication in the frequency domain.
    4. IFFT and take the first $n$ components.

??? proof "Proof"
    Construct the $m = 2n$ circulant matrix $C = \operatorname{circ}(c_0, c_1, \ldots, c_{2n-1})$ with $c_k = t_k$ ($0 \le k \le n-1$) and $c_k = t_{k-2n}$ ($n \le k \le 2n-1$). Pad $\mathbf{x}$ to $\tilde{\mathbf{x}} = (x_0, \ldots, x_{n-1}, 0, \ldots, 0)^T$. Then $(C\tilde{\mathbf{x}})_j = \sum_{k=0}^{n-1}c_{(j-k)\bmod 2n}x_k$. When $0 \le j \le n-1$, $(j-k) \bmod 2n$ yields exactly $t_{j-k}$, so the first $n$ components equal $T\mathbf{x}$. $\blacksquare$

!!! definition "Definition 30.7 (Yule-Walker equations)"
    In autoregressive (AR) modeling, given a symmetric positive definite Toeplitz matrix $T_n = [r_{|i-j|}]$ (the autocorrelation matrix), the **Yule-Walker equations** are

    $$
    T_n \mathbf{a} = -\mathbf{r}, \quad \mathbf{r} = (r_1, r_2, \ldots, r_n)^T,
    $$

    where $\mathbf{a}$ contains the AR coefficients. Due to the Toeplitz structure of $T_n$, this can be solved efficiently using Levinson recursion.

!!! theorem "Theorem 30.8 (Levinson-Durbin recursion)"
    For a positive definite symmetric Toeplitz system $T_n \mathbf{a} = \mathbf{b}$, the Levinson-Durbin algorithm solves it in $O(n^2)$ time (compared to $O(n^3)$ for general methods). The algorithm starts from a $1 \times 1$ system and recursively increases the dimension using the previous solution and the Toeplitz structure:

    Given that $T_k \mathbf{a}^{(k)} = \mathbf{b}^{(k)}$ has been solved, the reflection coefficient is

    $$
    \kappa_{k+1} = \frac{b_{k+1} - \mathbf{r}_k^T \mathbf{a}^{(k)}}{\epsilon_k},
    $$

    where $\epsilon_k$ is the prediction error power. The update is

    $$
    \mathbf{a}^{(k+1)} = \begin{pmatrix}\mathbf{a}^{(k)}\\ 0\end{pmatrix} + \kappa_{k+1}\begin{pmatrix}J_k \mathbf{a}^{(k)}\\ 1\end{pmatrix},
    $$

    where $J_k$ is the $k \times k$ reversal matrix.

??? proof "Proof"
    The key is exploiting the **centrosymmetry** of Toeplitz matrices: if $T_k\mathbf{a} = \mathbf{b}$, then $T_k(J_k\mathbf{a}) = J_k\mathbf{b}$. Constructing the $(k+1)$-dimensional solution from the $k$-dimensional one requires determining only one new parameter $\kappa_{k+1}$ (the reflection coefficient) to satisfy the extended system. This reduces the $O(k)$ new equations to a single-parameter problem, giving total complexity $\sum_{k=1}^n O(k) = O(n^2)$. $\blacksquare$

!!! example "Example 30.4"
    **Levinson recursion for AR modeling.** Let the autocorrelation sequence be $r_0 = 1$, $r_1 = 0.8$, $r_2 = 0.5$. The Yule-Walker equations:

    $$
    \begin{pmatrix}1&0.8\\0.8&1\end{pmatrix}\begin{pmatrix}a_1\\a_2\end{pmatrix} = -\begin{pmatrix}0.8\\0.5\end{pmatrix}.
    $$

    **Step 1**: $a_1^{(1)} = -r_1/r_0 = -0.8$, $\epsilon_1 = r_0(1 - 0.8^2) = 0.36$.

    **Step 2**: $\kappa_2 = \frac{-0.5 - 0.8 \cdot (-0.8)}{0.36} = \frac{-0.5 + 0.64}{0.36} = \frac{0.14}{0.36} \approx 0.389$.

    $$
    \mathbf{a}^{(2)} = \begin{pmatrix}-0.8\\0\end{pmatrix} + 0.389\begin{pmatrix}-0.8\\1\end{pmatrix} = \begin{pmatrix}-1.111\\0.389\end{pmatrix}.
    $$

    Prediction error power: $\epsilon_2 = 0.36(1 - 0.389^2) \approx 0.305$. Since $|\kappa_k| < 1$, the system is stable.

---

## 30.5 Linear Coding Theory

<div class="context-flow" markdown>

**Core idea**: Linear code = linear subspace over finite field $\mathbb{F}_q$ → generator matrix $G$ (column space) and parity-check matrix $H$ (null space) → minimum distance = minimum number of linearly dependent columns
**Link**: Ch4 vector spaces · Ch1 linear systems

</div>

The basic idea of coding theory is to introduce redundancy into transmitted information for error correction, and linear codes place this idea entirely on the foundation of linear algebra.

!!! definition "Definition 30.8 (Linear code)"
    An $[n, k]$ **linear code** $\mathcal{C}$ over the finite field $\mathbb{F}_q$ is a $k$-dimensional linear subspace of $\mathbb{F}_q^n$. Parameter interpretation:

    - $n$: code length (length of each codeword).
    - $k$: information dimension (number of information symbols encoded).
    - $n - k$: number of redundancy symbols.
    - **Code rate** $R = k/n$.

    $\mathcal{C}$ contains $q^k$ codewords.

!!! definition "Definition 30.9 (Generator and parity-check matrices)"
    The **generator matrix** $G \in \mathbb{F}_q^{k \times n}$ of an $[n, k]$ linear code $\mathcal{C}$ satisfies

    $$
    \mathcal{C} = \{\mathbf{m}G : \mathbf{m} \in \mathbb{F}_q^k\} = \operatorname{rowspace}(G).
    $$

    The **parity-check matrix** $H \in \mathbb{F}_q^{(n-k) \times n}$ satisfies

    $$
    \mathcal{C} = \ker(H) = \{\mathbf{c} \in \mathbb{F}_q^n : H\mathbf{c}^T = \mathbf{0}\}.
    $$

    $G$ and $H$ satisfy $GH^T = 0$. In systematic form $G = [I_k \mid P]$, we have $H = [-P^T \mid I_{n-k}]$.

!!! theorem "Theorem 30.9 (Minimum distance and error-correcting capability)"
    The **minimum distance** $d$ of a linear code $\mathcal{C}$ equals the smallest number of columns of $H$ that are linearly dependent: any $d - 1$ columns of $H$ are linearly independent, but some set of $d$ columns is linearly dependent. This is denoted as an $[n, k, d]$ code.

    - $\mathcal{C}$ can **detect** $d - 1$ errors.
    - $\mathcal{C}$ can **correct** $\lfloor(d-1)/2\rfloor$ errors.

??? proof "Proof"
    $d = \min_{\mathbf{c} \in \mathcal{C}, \mathbf{c} \ne \mathbf{0}} w(\mathbf{c})$ (minimum Hamming weight), where $w(\mathbf{c})$ is the number of nonzero components. $\mathbf{c} \in \mathcal{C}$ if and only if $H\mathbf{c}^T = \mathbf{0}$, meaning the columns of $H$ at the positions of the nonzero components of $\mathbf{c}$ are linearly dependent. Therefore $d$ is the minimum number of linearly dependent columns of $H$. For error correction, if $t$ errors occur ($t \le \lfloor(d-1)/2\rfloor$), the received word $\mathbf{r}$ is at distance $\le t$ from the true codeword and at distance $\ge d - t > t$ from any other codeword, enabling unique decoding. $\blacksquare$

!!! theorem "Theorem 30.10 (Singleton bound)"
    Any $[n, k, d]$ linear code satisfies

    $$
    d \le n - k + 1.
    $$

    Codes achieving this bound are called **maximum distance separable** (MDS) codes.

??? proof "Proof"
    $H \in \mathbb{F}_q^{(n-k) \times n}$ has $n$ columns, each in $\mathbb{F}_q^{n-k}$. Any $n - k + 1$ columns in an $(n-k)$-dimensional space must be linearly dependent, so $d \le n - k + 1$. $\blacksquare$

!!! example "Example 30.5"
    **A $[7, 4]$ binary linear code.** Generator matrix (systematic form):

    $$
    G = \begin{pmatrix}1&0&0&0&1&1&0\\0&1&0&0&0&1&1\\0&0&1&0&1&1&1\\0&0&0&1&1&0&1\end{pmatrix}, \quad H = \begin{pmatrix}1&0&1&1&1&0&0\\1&1&1&0&0&1&0\\0&1&1&1&0&0&1\end{pmatrix}.
    $$

    Verify $GH^T = 0$ (over $\mathbb{F}_2$). The information word $\mathbf{m} = (1,0,1,1)$ encodes to $\mathbf{c} = \mathbf{m}G = (1,0,1,1,1,0,0)$. Check: $H\mathbf{c}^T = (0,0,0)^T$.

---

## 30.6 Hamming Codes and Error Correction

<div class="context-flow" markdown>

**Core idea**: Hamming code parity-check matrix $H$ has columns = all nonzero vectors in $\mathbb{F}_2^r$ → syndrome $\mathbf{s} = H\mathbf{r}^T$ directly locates the error position → perfect code (achieves Hamming bound)
**Link**: 30.5 linear codes · Ch4 linear independence

</div>

Hamming codes are among the most classical error-correcting codes, with their construction elegantly embodying ideas from linear algebra.

!!! definition "Definition 30.10 (Hamming code)"
    For a positive integer $r \ge 2$, the $[2^r - 1, 2^r - 1 - r, 3]$ **Hamming code** has a parity-check matrix $H \in \mathbb{F}_2^{r \times (2^r-1)}$ whose columns are all $2^r - 1$ nonzero vectors in $\mathbb{F}_2^r$ (in some ordering).

    - Code length $n = 2^r - 1$.
    - Information bits $k = 2^r - 1 - r$.
    - Minimum distance $d = 3$ (corrects 1-bit errors).

!!! theorem "Theorem 30.11 (Properties of Hamming codes)"
    Hamming codes have the following properties:

    1. **$d = 3$**: Any two columns of $H$ (two distinct nonzero vectors in $\mathbb{F}_2^r$) are linearly independent, but there exist three linearly dependent columns.
    2. **Perfect code**: Hamming codes achieve the Hamming bound $\sum_{i=0}^{t}\binom{n}{i} \le 2^{n-k}$ (with equality when $t = 1$).
    3. **Syndrome decoding**: For received word $\mathbf{r} = \mathbf{c} + \mathbf{e}$, the syndrome is $\mathbf{s} = H\mathbf{r}^T = H\mathbf{e}^T$. If $\mathbf{e}$ is the unit vector $\mathbf{e}_j$ (error in position $j$), then $\mathbf{s} = H\mathbf{e}_j^T = \mathbf{h}_j$ (the $j$-th column of $H$), directly locating the error.

??? proof "Proof"
    (1) Two distinct nonzero vectors over $\mathbb{F}_2$ are linearly independent (otherwise they would be equal). For three columns $\mathbf{h}_i, \mathbf{h}_j, \mathbf{h}_k$ with $\mathbf{h}_i + \mathbf{h}_j = \mathbf{h}_k$ (such triples necessarily exist over $\mathbb{F}_2$), they are linearly dependent. Hence $d = 3$.

    (2) Hamming bound: to correct $t = 1$ error, we need $1 + n = 2^r$ cosets (no error + $n$ single-error patterns), which exactly equals $2^{n-k} = 2^r$ syndromes.

    (3) Follows directly from $\mathbf{s} = H\mathbf{e}^T$. $\blacksquare$

!!! example "Example 30.6"
    **Error correction with the $[7, 4, 3]$ Hamming code.** $r = 3$, parity-check matrix:

    $$
    H = \begin{pmatrix}0&0&0&1&1&1&1\\0&1&1&0&0&1&1\\1&0&1&0&1&0&1\end{pmatrix}.
    $$

    The columns are the binary representations of 1 through 7. Suppose codeword $\mathbf{c} = (0,1,1,0,0,1,0)$ is sent, and $\mathbf{r} = (0,1,1,0,1,1,0)$ is received (error in position 5).

    Syndrome: $\mathbf{s} = H\mathbf{r}^T = (1,0,1)^T$. Binary $(101)_2 = 5$, identifying position 5 as erroneous. Correction: flip bit 5 to get $\hat{\mathbf{c}} = (0,1,1,0,0,1,0) = \mathbf{c}$.

!!! example "Example 30.7"
    **Extended Hamming code.** Appending an overall parity bit to the $[7,4,3]$ Hamming code yields the $[8,4,4]$ **extended Hamming code**. The parity-check matrix becomes

    $$
    H_{\text{ext}} = \begin{pmatrix}1&1&1&1&1&1&1&1\\0&0&0&1&1&1&1&0\\0&1&1&0&0&1&1&0\\1&0&1&0&1&0&1&0\end{pmatrix}.
    $$

    With $d = 4$, it can correct 1-bit errors and detect 2-bit errors. Syndrome analysis: if $\mathbf{s} \ne \mathbf{0}$ and the first component is $1$, it is a single error (correctable); if the first component is $0$, it is a double error (detection only).

---

## 30.7 Reed-Solomon Codes

<div class="context-flow" markdown>

**Core idea**: RS codes over $\mathbb{F}_q$ ($q = 2^m$), generator matrix is a Vandermonde matrix → MDS code (achieves Singleton bound) → polynomial evaluation/interpolation perspective
**Link**: 30.5 linear codes · Ch3 determinants (Vandermonde)

</div>

Reed-Solomon codes are the most important MDS codes, widely used in CDs/DVDs, QR codes, and deep-space communications. Their algebraic structure is built on Vandermonde matrices.

!!! definition "Definition 30.11 (Reed-Solomon code)"
    Let $\mathbb{F}_q$ be a finite field with $q$ elements, and $\alpha_1, \alpha_2, \ldots, \alpha_n$ be $n$ distinct elements of $\mathbb{F}_q$ ($n \le q$). The $[n, k, n-k+1]$ **Reed-Solomon code** is defined as

    $$
    \mathcal{C}_{\text{RS}} = \{(p(\alpha_1), p(\alpha_2), \ldots, p(\alpha_n)) : p \in \mathbb{F}_q[x], \deg(p) < k\}.
    $$

    That is, polynomials of degree less than $k$ evaluated at $n$ points. The generator matrix is a Vandermonde matrix:

    $$
    G = \begin{pmatrix}1&1&\cdots&1\\\alpha_1&\alpha_2&\cdots&\alpha_n\\\alpha_1^2&\alpha_2^2&\cdots&\alpha_n^2\\\vdots&\vdots&\ddots&\vdots\\\alpha_1^{k-1}&\alpha_2^{k-1}&\cdots&\alpha_n^{k-1}\end{pmatrix} \in \mathbb{F}_q^{k \times n}.
    $$

!!! theorem "Theorem 30.12 (RS codes are MDS codes)"
    Reed-Solomon codes achieve the Singleton bound, i.e., $d = n - k + 1$, and are MDS codes. Equivalently:

    1. Any $k$ columns of $G$ are linearly independent (any $k$-column submatrix of a Vandermonde matrix has full rank).
    2. Any $n - k$ columns of $H$ are linearly independent.
    3. RS codes can correct $t = \lfloor(n-k)/2\rfloor$ symbol errors.

??? proof "Proof"
    Any $k$ columns of $G$ form a $k \times k$ Vandermonde matrix $V(\alpha_{i_1}, \ldots, \alpha_{i_k})$ with determinant $\det(V) = \prod_{s<t}(\alpha_{i_t} - \alpha_{i_s}) \ne 0$ (since the $\alpha_i$ are pairwise distinct). Therefore any $k$ columns of $G$ are linearly independent. This implies $d \ge n - k + 1$, and the Singleton bound gives $d \le n - k + 1$, so $d = n - k + 1$. $\blacksquare$

!!! theorem "Theorem 30.13 (Parity-check matrix of RS codes)"
    Choosing $\alpha_i = \alpha^{i-1}$ ($\alpha$ a primitive element of $\mathbb{F}_q$, $n = q - 1$), the parity-check matrix of the RS code is

    $$
    H = \begin{pmatrix}1&\alpha&\alpha^2&\cdots&\alpha^{n-1}\\1&\alpha^2&\alpha^4&\cdots&\alpha^{2(n-1)}\\\vdots&\vdots&\vdots&\ddots&\vdots\\1&\alpha^{n-k}&\alpha^{2(n-k)}&\cdots&\alpha^{(n-k)(n-1)}\end{pmatrix}.
    $$

    The syndrome components $s_j = r(\alpha^j)$ ($j = 1, \ldots, n-k$), where $r(x)$ is the polynomial representation of the received word. Decoding algorithms (such as Berlekamp-Massey or the Euclidean algorithm) use the syndrome to find the error-locator polynomial and error values.

??? proof "Proof"
    A codeword $\mathbf{c} = (p(\alpha^0), p(\alpha^1), \ldots, p(\alpha^{n-1}))$ with $\deg(p) < k$. The $j$-th component of $H\mathbf{c}^T$ is $\sum_{i=0}^{n-1}c_i\alpha^{ji} = \sum_{i=0}^{n-1}p(\alpha^i)\alpha^{ji}$. Since the roots of the code's generator polynomial include $\alpha, \alpha^2, \ldots, \alpha^{n-k}$, we have $H\mathbf{c}^T = \mathbf{0}$. $\blacksquare$

!!! example "Example 30.8"
    **RS code over $\mathbb{F}_8$.** Let $\mathbb{F}_8 = \mathbb{F}_2[\alpha]/(\alpha^3 + \alpha + 1)$, $n = 7$, $k = 3$, $d = 5$ (corrects 2 symbol errors).

    Information polynomial $p(x) = 1 + \alpha x + \alpha^2 x^2$. The codeword is obtained by evaluating at $1, \alpha, \alpha^2, \ldots, \alpha^6$:

    $$
    \mathbf{c} = (p(1), p(\alpha), p(\alpha^2), \ldots, p(\alpha^6)).
    $$

    Each component is an element of $\mathbb{F}_8$ (3 bits). Code rate $R = 3/7 \approx 0.43$. Since $d = 5$, the code corrects any 2 symbol errors (each symbol being 3 bits), meaning it can correct any 6-bit errors within 21 bits as long as the errors are confined to at most 2 symbols.

---

## 30.8 LDPC Codes and Sparse Matrices

<div class="context-flow" markdown>

**Core idea**: LDPC code = parity-check matrix $H$ is sparse → Tanner graph (bipartite graph) representation → belief propagation decoding exploits sparse structure → approaches Shannon limit
**Link**: 30.5 linear codes · Ch17 sparse/nonnegative matrices

</div>

Low-Density Parity-Check (LDPC) codes are among the best-performing code families in modern communication systems, with their defining characteristic being the sparsity of the parity-check matrix.

!!! definition "Definition 30.12 (LDPC code)"
    A **Low-Density Parity-Check code** (LDPC code) is an $[n, k]$ linear code whose parity-check matrix $H \in \mathbb{F}_2^{(n-k) \times n}$ is a **sparse matrix**: the number of 1s in each row and column is much smaller than $n$.

    - **Regular LDPC code**: each column of $H$ has exactly $d_v$ ones (variable node degree) and each row has exactly $d_c$ ones (check node degree), satisfying $(n-k)d_c = n d_v$.
    - **Irregular LDPC code**: row and column weights may vary, described by a degree distribution pair $(\lambda(x), \rho(x))$.

!!! definition "Definition 30.13 (Tanner graph)"
    The **Tanner graph** of an LDPC code is a bipartite graph $G = (V \cup C, E)$:

    - **Variable nodes** $V = \{v_1, \ldots, v_n\}$, corresponding to the $n$ columns of $H$ (the $n$ bits of the codeword).
    - **Check nodes** $C = \{c_1, \ldots, c_{n-k}\}$, corresponding to the $n - k$ rows of $H$ (the parity-check equations).
    - **Edges**: $v_j$ is connected to $c_i$ if and only if $H_{ij} = 1$.

    The sparsity of $H$ means the Tanner graph has $|E| = \text{nnz}(H) \ll n(n-k)$ edges.

!!! theorem "Theorem 30.14 (Minimum distance of LDPC codes)"
    For a regular $(d_v, d_c)$-LDPC code:

    1. **Lower bound**: $d \ge d_v + 1$ (any $d_v$ columns of $H$ are linearly independent since each column has exactly $d_v$ ones and sparsity prevents small sets from being dependent).
    2. The **girth** $g$ of the Tanner graph (length of the shortest cycle) affects performance: larger girth leads to better codes. A rough lower bound is $d \ge g/2 + 1$ (for regular codes).
    3. For randomly constructed LDPC codes, as $n \to \infty$, the minimum distance grows linearly with high probability: $d = \Theta(n)$.

??? proof "Proof"
    (1) Let $\mathbf{c}$ be a nonzero codeword of weight $w$. $H\mathbf{c}^T = \mathbf{0}$ means the $w$ columns of $H$ corresponding to the support of $\mathbf{c}$ are linearly dependent over $\mathbb{F}_2$. Since each column has $d_v$ ones, the sum ($\bmod 2$) of $w$ columns has support size $\le wd_v$. For the sum to be the zero vector requires enough columns to "cover" all nonzero rows. When $w \le d_v$, the nonzero row positions of $w$ columns number at most $wd_v$, and each row is covered by at most $d_c$ columns. Careful analysis shows that $w \le d_v$ columns cannot be linearly dependent (under the condition that the Tanner graph has no short cycles), so $d \ge d_v + 1$. $\blacksquare$

!!! theorem "Theorem 30.15 (Belief propagation decoding)"
    LDPC codes are decoded by **Belief Propagation** (BP), a message-passing algorithm on the Tanner graph. Let $L_j$ be the channel log-likelihood ratio (LLR) for variable node $v_j$. The decoding iterations proceed as:

    **Variable → check messages**:

    $$
    \mu_{v_j \to c_i}^{(\ell)} = L_j + \sum_{c \in \mathcal{N}(v_j) \setminus c_i} \mu_{c \to v_j}^{(\ell-1)},
    $$

    **Check → variable messages**:

    $$
    \mu_{c_i \to v_j}^{(\ell)} = 2\operatorname{atanh}\!\left(\prod_{v \in \mathcal{N}(c_i) \setminus v_j} \tanh\!\left(\frac{\mu_{v \to c_i}^{(\ell)}}{2}\right)\right),
    $$

    where $\mathcal{N}(\cdot)$ denotes the neighbor set in the Tanner graph. Each iteration has complexity $O(|E|) = O(\text{nnz}(H))$, guaranteed to be efficient by the sparsity of $H$.

??? proof "Proof"
    The BP algorithm is based on the sum-product algorithm on factor graphs. For the AWGN channel, the received value is $y_j = (1 - 2c_j) + n_j$, giving LLR $L_j = 2y_j/\sigma^2$. Variable-to-check messages aggregate information from the channel and other check nodes. Check-to-variable messages implement the parity constraint $\bigoplus_{v \in \mathcal{N}(c_i)}c_v = 0$ via the $\tanh$ rule in the LLR domain. When the Tanner graph is cycle-free, BP yields exact posterior probabilities; with cycles, BP is approximate, but in practice it performs excellently for sufficiently long LDPC codes. $\blacksquare$

!!! example "Example 30.9"
    **Tanner graph of a simple LDPC code.** Consider a $[6, 3]$ LDPC code with parity-check matrix:

    $$
    H = \begin{pmatrix}1&1&1&0&0&0\\0&0&1&1&1&0\\1&0&0&0&1&1\end{pmatrix}.
    $$

    Each column has at most $2$ ones, and each row has $3$ ones. The Tanner graph has $6$ variable nodes and $3$ check nodes with $\text{nnz}(H) = 9$ edges.

    This $H$ defines a valid $[6,3]$ linear code. Consider received word $\mathbf{r} = (1,0,1,0,1,0)$: syndrome $\mathbf{s} = H\mathbf{r}^T = (0, 0, 0)^T$, so $\mathbf{r}$ is a valid codeword. If instead $\mathbf{r} = (1,1,1,0,1,0)$, the syndrome is $\mathbf{s} = (1, 0, 1)^T \ne \mathbf{0}$, detecting an error. BP decoding then uses message passing on the Tanner graph to locate and correct the error.
