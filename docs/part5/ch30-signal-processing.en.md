# Chapter 30: Linear Algebra in Signal Processing and Coding

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Decompositions (Ch10-11) · Kronecker Product (Ch19) · Numerical Stability (Ch22)

**Chapter Outline**: The DFT Matrix (Unitary Transforms) → Fast Fourier Transform (FFT) as Matrix Factorization → Circulant Matrices & Convolution → Filtering in the Frequency Domain → Toeplitz Matrices & LTI Systems → Sampling Theorem & Sinc Interpolation → Wavelet Transform & Sparse Representations → Compressed Sensing (Basis Pursuit) → Linear Error-Correcting Codes (Finite Field Algebra)

**Extension**: Signal processing is essentially the art of changing bases; the Fourier transform rotates a signal into the "spectral basis" where convolution becomes simple multiplication.

</div>

Almost everything in modern telecommunications, audio processing, and image compression is built on the matrix representation of signals. The Discrete Fourier Transform (DFT) is a unitary matrix multiplication, and digital filters are linear operators defined by Toeplitz or circulant structures. This chapter explores how linear algebra provides the engine for extracting information from waves and protecting data during transmission.

---

## 30.1 The Discrete Fourier Transform (DFT)

<div class="context-flow" markdown>

**Algebraic Insight**: The DFT is not just a formula; it is a change of basis to a set of orthogonal complex sinusoids.

</div>

!!! definition "Definition 30.1 (The DFT Matrix)"
    The $n \times n$ **DFT matrix** $F_n$ has entries:
    $$(F_n)_{jk} = \omega^{jk}, \quad \omega = e^{-2\pi i / n}$$
    The DFT of a signal $\mathbf{x}$ is $\hat{\mathbf{x}} = F_n \mathbf{x}$.

!!! theorem "Theorem 30.1 (Unitary Property)"
    The normalized DFT matrix $\frac{1}{\sqrt{n}} F_n$ is a **unitary matrix**. This implies that the DFT preserves the total energy of the signal (Parseval's Theorem).

!!! technique "FFT as Matrix Factorization"
    The **Fast Fourier Transform (FFT)** is a method to factor $F_n$ into $O(\log n)$ sparse matrices. This reduces the complexity of signal transformation from $O(n^2)$ to $O(n \log n)$.

---

## 30.2 Circulant Matrices and Convolution

!!! definition "Definition 30.2 (Circulant Matrix)"
    A matrix $C$ is **circulant** if each row is a cyclic shift of the previous one. 
    **Key Fact**: Every circulant matrix is diagonalized by the DFT matrix $F_n$:
    $$C = F_n^{-1} \Lambda F_n$$
    where $\Lambda$ is a diagonal matrix containing the eigenvalues of $C$ (which are exactly the DFT of the first row of $C$).

!!! theorem "Theorem 30.2 (The Convolution Theorem)"
    The cyclic convolution of two signals $\mathbf{x} \circledast \mathbf{h}$ is equivalent to the element-wise multiplication of their DFTs:
    $$\mathcal{F}(\mathbf{x} \circledast \mathbf{h}) = \mathcal{F}(\mathbf{x}) \cdot \mathcal{F}(\mathbf{h})$$
    This is why we filter signals in the frequency domain.

---

## 30.3 Sampling and Interpolation

!!! theorem "Theorem 30.3 (Nyquist-Shannon Sampling Theorem)"
    A band-limited signal can be perfectly reconstructed from its samples if it is sampled at least twice its highest frequency. Reconstruction is a linear process involving **Sinc Interpolation**, which can be expressed as a matrix-vector product with a Sinc matrix.

---

## 30.4 Compressed Sensing

<div class="context-flow" markdown>

**Paradigm Shift**: Can we recover a signal from fewer samples than Nyquist requires? Yes, if the signal is **sparse**.

</div>

!!! technique "Basis Pursuit"
    Given $y = Ax$ where $A$ is a wide measurement matrix ($m \ll n$), we recover $x$ by solving:
    $$\min \|x\|_1 \quad \text{subject to } Ax = y$$
    Linear algebra ensures that if $A$ satisfies the **Restricted Isometry Property** (RIP), the sparse signal is exactly recovered.

---

## 30.5 Linear Error-Correcting Codes

!!! definition "Definition 30.3 (Linear Code)"
    A linear code $\mathcal{C}$ is a subspace of a vector space over a finite field (usually $\mathbb{F}_2$). It is defined by a **Generator Matrix** $G$ (columns span the code) and a **Parity Check Matrix** $H$ ($Hx = 0$ for all codewords).

---

## Exercises


****
??? success "Solution"
     $\omega = e^{-i\pi} = -1$. $F_2 = \begin{pmatrix} \omega^0 & \omega^0 \\ \omega^0 & \omega^1 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$.


****
??? success "Solution"
     The eigenvalues are the DFT coefficients of its first row.


****
??? success "Solution"
     Two FFTs ($O(n \log n)$), one element-wise multiplication ($O(n)$), and one inverse FFT ($O(n \log n)$). Total: $O(n \log n)$.


****
??? success "Solution"
     A Linear Time-Invariant (LTI) system with finite-length signals.


****
??? success "Solution"
     Its sparsity level is $s = 5$.


****
??? success "Solution"
     The number of non-zero entries. For linear codes, the minimum Hamming distance of the code equals the minimum Hamming weight of any non-zero codeword.


****
??? success "Solution"
     Fourier uses global sinusoids (frequency info only); Wavelets use localized pulses (time and frequency info). DWT matrix is typically sparse and orthogonal.


****
??? success "Solution"
     It guarantees that the measurement matrix $A$ preserves the distances between sparse vectors, making reconstruction robust.


****
??? success "Solution"
     $k$.

****
??? success "Solution"
    ## Chapter Summary

Linear algebra is the "mother tongue" of signal processing:


****: Recognized the DFT as a rotation into a spectral basis, where signals are decomposed into fundamental frequencies.

****: Showed how circulant and Toeplitz structures allow for FFT-based algorithms that are exponentially faster than naive methods.

****: Developed the theory of compressed sensing, proving that linear algebra can reconstruct information from sub-Nyquist data streams.

****: Framed error correction as a subspace problem over finite fields, where parity checks are linear constraints.
