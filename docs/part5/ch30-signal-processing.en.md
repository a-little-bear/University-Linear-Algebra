# Chapter 30: Linear Algebra in Signal Processing

<div class="context-flow" markdown>

**Prerequisites**: Matrix Multiplication (Ch2) · Complex Numbers · Toeplitz Matrices (Ch37A) · Orthogonality (Ch7)

**Chapter Outline**: Signals as Vectors → Linear Time-Invariant (LTI) Systems as Toeplitz Matrices → Convolution and Matrix-Vector Multiplication → Discrete Fourier Transform (DFT) as a Basis Change → The Vandermonde structure of the DFT Matrix → Circulant Matrices and Spectral Decomposition → Fast Fourier Transform (FFT) Complexity → Filter Banks and Sparse Sampling

**Extension**: Digital signal processing is the study of linear operators acting on discrete-time sequences; the DFT is the "Eigenbasis" of every shift-invariant system.

</div>

Signal processing is the application of linear algebra to time-varying data. A signal is viewed as a vector, and a **linear filter** is an operator that transforms this vector. The core discovery of this field is that for systems that are **time-invariant**, the complex exponentials (sine and cosine waves) are the universal eigenvectors. This leads to the **Discrete Fourier Transform (DFT)**, a change of basis into the "frequency domain" that diagonalizes all convolution-type operators. This chapter explores the matrix structure of the DFT and the algebraic reasons why the FFT is so efficient.

---

## 30.1 Signals, Systems, and the DFT

!!! definition "Definition 30.1 (LTI System)"
    A linear time-invariant system is an operator $H$ such that if $y = Hx$, then $y$ is the **convolution** of $x$ with the impulse response $h$. In matrix form, $H$ is a **Toeplitz matrix**.

!!! theorem "Theorem 30.1 (The DFT Matrix)"
    The DFT of a signal $x$ is $X = Fx$, where $F$ is the unitary Vandermonde matrix:
    $$F_{nk} = \frac{1}{\sqrt{N}} e^{-j \frac{2\pi}{N} nk}$$
    $F$ diagonalizes all circulant matrices.

---

## Exercises

1. **[Fundamentals] Write the $2 \times 2$ DFT matrix.**
   ??? success "Solution"
       $F_2 = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$. This is the same as the normalized Hadamard matrix.

2. **[Convolution] Why is convolution equivalent to matrix multiplication by a circulant matrix?**
   ??? success "Solution"
       For periodic signals, the output at time $n$ is $y_n = \sum_k h_{n-k} x_k$. The indices $n-k \pmod N$ create a circulant pattern in the matrix $H$.

3. **[DFT Inverse] Show that the inverse of the DFT matrix $F$ is its conjugate transpose $F^*$.**
   ??? success "Solution"
       The columns of $F$ are orthonormal. Since $F$ is a unitary matrix ($F^* F = I$), its inverse is its adjoint.

4. **[Fast Fourier Transform] Compare the complexity of $Fx$ versus the FFT.**
   ??? success "Solution"
       Direct matrix-vector multiplication is $O(N^2)$. The FFT exploits the recursive structure of $F$ to achieve $O(N \log N)$. For $N=10^6$, this is the difference between minutes and milliseconds.

5. **[Filtering] How is a signal filtered in the frequency domain?**
   ??? success "Solution"
       1. Transform signal to frequency domain: $X = Fx$. 2. Multiply by the filter's transfer function (eigenvalues): $Y = \Lambda X$. 3. Transform back: $y = F^* Y$. This is $y = F^* \Lambda F x$.

6. **[Parseval] State Parseval's Theorem in matrix terms.**
   ??? success "Solution"
       $\|x\|^2 = \|Fx\|^2$. Unitary transformations (like the DFT) preserve the total energy of the signal.

7. **[Vandermonde] Why is $F$ a Vandermonde matrix?**
   ??? success "Solution"
       The entries are powers of $w = e^{-j 2\pi/N}$: $F_{nk} = (w^n)^k$. It is the evaluation of a polynomial at the $N$-th roots of unity.

8. **[Toeplitz] Relate LTI systems to Toeplitz matrices.**
   ??? success "Solution"
       Non-periodic convolution is represented by a Toeplitz matrix. While not perfectly diagonalized by the DFT, large Toeplitz matrices are asymptotically diagonalized by it (Szegő theory).

9. **[Sparsity] What does it mean for a signal to be "sparse" in the frequency domain?**
   ??? success "Solution"
       It means $X = Fx$ has only a few non-zero entries. This allows for **Compressed Sensing**—recovering the signal from far fewer samples than the Nyquist limit suggests.

10. **[Sampling] How does sampling relate to the null space of a transformation?**
    ??? success "Solution"
        Under-sampling a signal is a projection that collapses the space. Aliasing occurs when the original signal has a component in the null space of the sampling operator.

## Chapter Summary

This chapter explores the spectral decomposition of signals:

1. **Vectorized Waveforms**: Represented discrete-time signals as vectors and linear filters as shift-invariant operators.
2. **Frequency Duality**: Established the DFT as the universal change of basis that diagonalizes convolution.
3. **Recursive Symmetry**: Analyzed the FFT as the algorithmic exploitation of matrix redundancy.
4. **Energy Conservation**: Utilized the unitary property of the DFT matrix to guarantee the preservation of signal power across domains.
