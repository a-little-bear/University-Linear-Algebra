# Chapter 30: Linear Algebra in Signal Processing and Coding

<div class="context-flow" markdown>

**Prerequisites**: Unitary Matrices (Ch55A) · Circulant Matrices (Ch37A) · Finite Fields (Ch52) · Orthogonality (Ch07)

**Chapter Outline**: Signals as Vectors → Algebraic Description of Linear Time-Invariant (LTI) Systems → Matrix Form of Convolution: Circulant and Toeplitz Matrices → Core Transform: Matrix Representation of the Discrete Fourier Transform (DFT) via Vandermonde Matrices → Spectral Analysis: Frequency as Eigenvalues → Wavelet Transforms and Orthonormal Bases → Linear Algebraic Interpretation of the Sampling Theorem → Coding Theory: Generator and Parity-Check Matrices → Subspace View of Error-Correction → Applications: Audio Compression (MP3), Digital Watermarking, and Channel Coding in 5G

**Extension**: Signal processing is the "dynamic implementation" of linear algebra; it transforms frequency analysis into operator diagonalization and information redundancy into subspace encoding. It is the bridge between physical waves and digital logic.

</div>

In the modern digital world, whether listening to music or sending a text, high-frequency linear algebra operations are occurring behind the scenes. **Signal Processing** treats time series as vectors and filters as matrices. **Coding Theory** protects information by constructing specific subspaces over finite fields. This chapter introduces the matrix essence of the Discrete Fourier Transform and how linear algebra makes communication both fast and accurate.

---

## 30.1 The Matrix Essence of DFT

!!! definition "Definition 30.1 (DFT Matrix $F_n$)"
    The $n$-point Discrete Fourier Transform is a linear transformation with entries given by powers of the root of unity:
    $$(F_n)_{jk} = \frac{1}{\sqrt{n}} \omega^{jk}, \quad \omega = e^{-2\pi i/n}$$
    **Property**: $F_n$ is a **unitary matrix**. This ensures the frequency domain representation fully preserves the energy of the time-domain signal (Parseval's Theorem).

---

## 30.2 Convolution and Circulant Matrices

!!! theorem "Theorem 30.1 (Matrix Convolution Theorem)"
    Discrete convolution is equivalent to multiplying a vector by a **Circulant Matrix**. Since circulant matrices are always diagonalized by the DFT matrix, this explains why "convolution in time equals pointwise multiplication in frequency."

---

## 30.3 Algebraic Structure of Error-Correcting Codes

!!! technique "Subspace Coding"
    An $(n, k)$ linear error-correcting code is a $k$-dimensional subspace of $GF(q)^n$.
    - **Error-Correction Capability**: Determined by the minimum Hamming distance between vectors in the subspace.
    - **Syndrome**: Utilizes the null space property of the parity-check matrix $H$ to locate errors via $H\mathbf{r}^T$.

---

## Exercises

**1. [Basics] Construct the $2 \times 2$ DFT matrix $F_2$.**

??? success "Solution"
    **Steps:**
    1. $n=2, \omega = e^{-i\pi} = -1$.
    2. Indices $j, k \in \{0, 1\}$.
    3. Entries: $F_{00}=1, F_{01}=1, F_{10}=1, F_{11}=-1$.
    4. Normalize by $1/\sqrt{2}$.
    **Conclusion**: $F_2 = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$. Note: This is exactly the Hadamard matrix.

**2. [Spectral Analysis] Prove that the eigenvalues of a circulant matrix are obtained by the DFT of its first row.**

??? success "Solution"
    **Proof Strategy:**
    1. Let $C$ be a circulant matrix. We know $C = F^* \Lambda F$.
    2. The diagonals of $\Lambda$ are the eigenvalues.
    3. The first row $c$ of $C$ is $\mathbf{e}_1^T C$.
    4. Substituting the diagonalization form and using the structure of $F$, one finds that the eigenvalues are precisely the basis change of $c$ under $F$.

**3. [Calculation] For signal vector $x = (1, 0, 1, 0)^T$, find its DC component (0-th frequency) after $F_4$ transform.**

??? success "Solution"
    **Calculation:**
    1. The 0-th frequency corresponds to the first row of the DFT matrix, which is all $1/2$ (after normalization).
    2. $X_0 = \frac{1}{2}(1\cdot 1 + 0\cdot 1 + 1\cdot 1 + 0\cdot 1) = 1$.
    **Conclusion**: The DC component is 1.

**4. [Sampling] What does the Sampling Theorem mean in terms of vector spaces?**

??? success "Solution"
    **Algebraic Reason:**
    It means that if we know a continuous function lies in a low-dimensional subspace (band-limited), then the vector obtained by sampling at specific points is sufficient to uniquely reconstruct the original function. It essentially guarantees that the sampling operator is **injective** on the band-limited subspace.

**5. [Wavelets] Briefly describe the matrix characteristic of the Haar wavelet transform relative to Fourier.**

??? success "Solution"
    The Haar wavelet matrix is sparse and has a recursive block structure.
    **Contrast**: Fourier bases are globally supported (each basis vector affects all time points), while wavelet bases are **locally supported**. This makes wavelets superior for representing abrupt changes and for image compression (JPEG 2000).

**6. [Coding] Given parity-check matrix $H = \begin{pmatrix} 1 & 1 & 1 \end{pmatrix}$ over $GF(2)$, is received vector $\mathbf{r} = (1, 0, 1)$ valid?**

??? success "Solution"
    **Calculate Syndrome:**
    $H \mathbf{r}^T = 1\cdot 1 + 1\cdot 0 + 1\cdot 1 = 1+0+1 = 0 \pmod 2$.
    **Conclusion**: Since the syndrome is 0, the received vector lies in the code subspace and is deemed error-free.

**7. [Property] Prove the inverse of the DFT matrix satisfies $F^{-1} = F^*$.**

??? success "Solution"
    **Proof:**
    1. The $(j, k)$ entry of $F F^*$ is $\frac{1}{n} \sum_m \omega^{jm} \overline{\omega^{km}} = \frac{1}{n} \sum_m \omega^{m(j-k)}$.
    2. Using the root of unity summation property: if $j=k$, the sum is $n$; if $j \neq k$, the sum is 0.
    3. Thus $F F^* = I$.
    **Conclusion**: The DFT matrix is a quintessential unitary matrix.

**8. [Application] What is "Zero Forcing" in channel equalization?**

??? success "Solution"
    **Algebraic Essence:**
    Channel interference is modeled as $y = Hx + n$. Zero forcing cancels interference by multiplying with the **inverse or pseudoinverse** $H^+$ of the channel matrix: $\hat{x} = H^+ y$. This directly applies linear system theory to signal recovery.

**9. [Filters] What shape does an ideal low-pass filter take in the frequency domain?**

??? success "Solution"
    **Conclusion: A Diagonal Matrix.**
    In the frequency basis, filtering (frequency weighting) manifests as $y_f = \Lambda X_f$. For an ideal low-pass, $\Lambda$ has 1s for the first $k$ diagonals and 0s elsewhere. This is essentially an **orthogonal projection operator**.

**10. [Efficiency] Why is the Fast Fourier Transform (FFT) so important?**

??? success "Solution"
    **Efficiency Leap:**
    Standard matrix-vector multiplication is $O(n^2)$. FFT exploits the recursive symmetry of the DFT matrix (Danielson-Lanczos Lemma) to factor it into a product of sparse matrices, reducing complexity to **$O(n \log n)$**. This algebraic acceleration is the cornerstone of modern real-time communication (WiFi, 4G/5G).

## Chapter Summary

Signal processing and coding are the "engineered realizations" of linear algebra:

1.  **Algebraic Frequency**: Through the DFT matrix, we prove that frequency analysis is essentially the eigen-decomposition of a shift operator, providing mathematical tools for understanding wave phenomena.
2.  **Spatial Protection of Information**: Coding theory demonstrates how to design the "distance" between subspaces in finite field vector spaces to achieve algebraic immunity to physical noise.
3.  **Efficiency of Structure**: From the diagonalization of convolution to the divide-and-conquer logic of FFT, linear algebra provides not just the language for description, but also the high-speed algorithms that power the Information Age.
