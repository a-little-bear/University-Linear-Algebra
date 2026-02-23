# Chapter 30: Applications of Linear Algebra in Signal Processing and Coding

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces (Ch8) · Unitary Matrices (Ch8) · Matrix Analysis (Ch14) · Finite Fields (Ch52)

**Chapter Outline**: Basis Expansion of Signals → Discrete Fourier Transform (DFT) and Circulant Matrices → Discrete Cosine Transform (DCT) and Compression → Convolution and Toeplitz Matrices → Algebraic Description of Sampling Theorem → Filter Design and System Functions → Linear Error-Correcting Codes → Generator and Parity-check Matrices → Compressed Sensing and Sparse Bases

**Extension**: Signal processing is a frequency-domain variant of linear algebra; Fourier Transform is essentially a unitary transformation from time-domain basis to frequency-domain (eigenbasis).

</div>

Signal processing and coding theory are the mathematical foundations of the information age. Through basis transformations in linear space, we decompose entangled time-domain signals into clear spectral components. Linear algebra provides the language for signal analysis and ensures info reliability through error-correcting codes.

---

## 30.1 Core Transforms and Algebraic Structures

!!! definition "Definition 30.1 (DFT Matrix)"
    The $n$-th order DFT matrix $F$ has entries $F_{jk} = \frac{1}{\sqrt{n}} \omega^{jk}$, where $\omega = e^{-2\pi i / n}$. It is a unitary matrix.

!!! theorem "Theorem 30.3 (Convolution Theorem and Matrices)"
    Convolution of two sequences corresponds to multiplication of their associated circulant matrices. Since circulant matrices are diagonalized by the DFT matrix, time-domain convolution is equivalent to element-wise multiplication in the frequency domain.

---

## Exercises

1. **[DFT Matrix] Prove that the 2nd-order DFT matrix $F_2 = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \ 1 & -1 \end{pmatrix}$ is unitary.**
   ??? success "Solution"
       $F_2 F_2^* = \frac{1}{2} \begin{pmatrix} 1 & 1 \ 1 & -1 \end{pmatrix} \begin{pmatrix} 1 & 1 \ 1 & -1 \end{pmatrix} = \frac{1}{2} \begin{pmatrix} 2 & 0 \ 0 & 2 \end{pmatrix} = I$.
       Since $F_2^* F_2 = I$, it is unitary.

2. **[Toeplitz Matrix] Why can Linear Time-Invariant (LTI) systems be represented by Toeplitz matrices?**
   ??? success "Solution"
       The core of LTI systems is convolution: $y[n] = \sum h[k]x[n-k]$. This means the response depends only on the relative displacement between input and impulse response. A Toeplitz matrix, where each row is a shift of the previous one, perfectly captures this time-invariance.

3. **[Spectral Energy] Prove the matrix form of Parseval's Theorem: $\|Fx\|_2 = \|x\|_2$.**
   ??? success "Solution"
       Since $F$ is unitary, it preserves vector length (inner product). This means the total energy of a signal in the time domain equals its energy in the frequency domain.

4. **[Data Compression] Why use DCT instead of DFT in JPEG image compression?**
   ??? success "Solution"
       DCT has superior "energy compaction" properties. For typical image blocks, DCT concentrates most energy into a few low-frequency eigenvalues, allowing efficient compression by discarding high-frequency coefficients (components corresponding to small singular values).

5. **[Coding over Finite Fields] Let $G = \begin{pmatrix} 1 & 0 & 1 & 1 \ 0 & 1 & 0 & 1 \end{pmatrix}$ be a generator matrix. Find the codeword for info bits $[1, 1]$.**
   ??? success "Solution"
       Codeword $c = [1, 1] G = [1, 1, 1, 0] \pmod 2$. In linear coding, codewords are linear combinations of generator matrix rows.

6. **[Parity-check Matrix] Given a parity-check matrix $H$, what linear equation must a codeword $c$ satisfy?**
   ??? success "Solution"
       $H c^T = 0$. That is, all valid codewords must belong to the nullspace of the parity-check matrix.

7. **[Matrix Rank and Sampling] In Compressed Sensing, what property must measurement matrix $A$ satisfy to recover signals from few samples?**
   ??? success "Solution"
       It must satisfy the Restricted Isometry Property (RIP). This means any $k$ columns extracted from $A$ form a submatrix that is nearly orthogonal (singular values near 1), ensuring sparse signals don't collapse after dimensionality reduction.

8. **[Fast Transforms] Why is FFT (Fast Fourier Transform) extremely important in engineering?**
   ??? success "Solution"
       Direct matrix-vector multiplication for order $n$ takes $O(n^2)$. FFT exploits the block symmetry and recursive structure of the DFT matrix (Danielson-Lanczos Lemma) to reduce complexity to $O(n \log n)$, making real-time processing possible.

9. **[Sampling Theorem] Describe the sampling process as a projection.**
   ??? success "Solution"
       Sampling can be viewed as projecting a continuous signal onto a subspace spanned by $\operatorname{sinc}$ functions (or impulse sequences). Perfect reconstruction is possible if the original signal lies within the span of this subspace.

10. **[Matrix Spectrum] If the first row of a circulant matrix $C$ is $[c_0, c_1, \dots, c_{n-1}]$, what are its eigenvalues?**
    ??? success "Solution"
        The eigenvalues are precisely the DFT of the sequence $[c_0, \dots, c_{n-1}]$. This establishes a direct algebraic link between circulant structure and frequency response.

## Chapter Summary

Signal processing is the frequency-domain art of linear algebra:

1. **Change of Basis**: All signal analysis is essentially finding the "eigenbasis" that best reveals physical laws.
2. **Structural Symmetry**: Circulant and Toeplitz structures transform time evolution into matrix algebra.
3. **Redundancy and Recovery**: Coding theory utilizes nullspaces and low-rank approximation to establish algebraic boundaries for info survival in noise.
