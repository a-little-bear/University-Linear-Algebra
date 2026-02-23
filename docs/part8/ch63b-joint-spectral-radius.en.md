# Chapter 63B: Joint Spectral Radius

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues and Spectral Radius (Ch6) · Matrix Norms (Ch15) · Simultaneous Triangularization (Ch63A)

**Chapter Outline**: Definition (Rota-Strang) → Norm-independence → Generalized Spectral Radius → Berger-Wang Theorem → Finiteness Conjecture → Undecidability → Extremal Norms (Barabanov) → Switched Linear Systems Stability → Wavelet Regularity

**Extension**: The joint spectral radius (JSR) is the central tool for analyzing the stability of switched systems, the convergence of subdivision schemes in wavelets, and the capacity of codes.

</div>

For a set of matrices $\Sigma = \{A_1, \ldots, A_m\}$, the **joint spectral radius** $\rho(\Sigma)$ characterizes the maximum growth rate of all possible product sequences $A_{i_1} A_{i_2} \cdots A_{i_k}$ as $k \to \infty$. This concept generalizes the spectral radius of a single matrix to matrix families, quantifying the "worst-case" growth rate.

---

## 63B.1 Definitions and the Berger-Wang Theorem

!!! definition "Definition 63B.1 (Joint Spectral Radius)"
    The joint spectral radius of a finite set of matrices $\Sigma$ is defined as:
    $$\rho(\Sigma) = \lim_{k \to \infty} \max_{A \in \Sigma^k} \|A\|^{1/k}$$
    where $\Sigma^k$ is the set of all products of length $k$. This limit exists and is independent of the choice of matrix norm $\|\cdot\|$.

!!! theorem "Theorem 63B.3 (Berger-Wang Theorem, 1992)"
    For finite sets of matrices, the joint spectral radius equals the supremum of the spectral radii of all finite products:
    $$\rho(\Sigma) = \sup_{k \ge 1} \max_{A \in \Sigma^k} \rho(A)^{1/k}$$

---

## Exercises

1. **[Concept] Contrast the joint spectral radius (JSR) with the standard spectral radius.**
   ??? success "Solution"
       The standard spectral radius $\rho(A)$ describes the growth of a single matrix power $A^k$. JSR describes the maximum possible growth rate across all sequences in a family $\Sigma$. JSR captures the aggregate instability of a collection of possible evolution rules.

2. **[Norm Bound] Prove: For any matrix norm $\|\cdot\|$, $\rho(\Sigma) \le \max_{A \in \Sigma} \|A\|$.**
   ??? success "Solution"
       By submultiplicativity, $\|A_{i_1} \dots A_{i_k}\| \le \prod \|A_{i_j}\| \le (\max \|A\|)^k$. Taking the $k$-th root and the limit as $k \to \infty$ yields $\rho(\Sigma) \le \max \|A\|$.

3. **[Calculation] Compute the JSR of $\Sigma = \left\{ \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}, \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix} \right\}$.**
   ??? success "Solution"
       Both matrices have $\rho(A_i) = 1$. However, their product $A_1 A_2 = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$ has eigenvalues $\frac{3\pm\sqrt{5}}{2}$, so $\rho(A_1 A_2) = \frac{3+\sqrt{5}}{2} \approx 2.618$. By Berger-Wang, $\rho(\Sigma) \ge \rho(A_1 A_2)^{1/2} = \frac{1+\sqrt{5}}{2} \approx 1.618$.

4. **[Berger-Wang Significance] Explain why the Berger-Wang Theorem is fundamental for JSR estimation.**
   ??? success "Solution"
       It proves that the external estimate (based on norms) and the internal estimate (based on eigenvalues of products) coincide in the limit. This allows for approximating the JSR by inspecting the spectra of a finite number of matrix products.

5. **[Stability] Why is the absolute stability of a switched linear system $x_{k+1} = A_{\sigma(k)} x_k$ equivalent to $\rho(\Sigma) < 1$?**
   ??? success "Solution"
       If $\rho(\Sigma) \ge 1$, there exists at least one switching sequence $\sigma$ such that the state does not converge to zero. Only when the worst-case growth rate is strictly less than 1 can the system be guaranteed to be safe under arbitrary switching.

6. **[CQLF] Describe the implication of a Common Quadratic Lyapunov Function (CQLF) for the JSR.**
   ??? success "Solution"
       A CQLF exists if there is a $P \succ 0$ such that $A_i^T P A_i \prec P$ for all $i$. This implies $\rho(\Sigma) < 1$. CQLF provides a sufficient (but not necessary) condition for stability that is computationally verifiable via SDP.

7. **[Undecidability] Explain why calculating the JSR is an "undecidable" problem.**
   ??? success "Solution"
       Determining whether $\rho(\Sigma) \le 1$ for integer matrices is equivalent to the Halting Problem. No universal algorithm can decide this property for all inputs in finite time, meaning we must rely on approximation schemes.

8. **[Finiteness] What is the "Finiteness Conjecture"?**
   ??? success "Solution"
       The conjecture proposed that $\rho(\Sigma)$ is always realized by the spectral radius of some finite product in $\Sigma^k$. While true for many classes (e.g., nonnegative matrices), it has been disproven for general matrix sets.

9. **[Wavelet Smoothness] How does the JSR relate to the differentiability of a wavelet?**
   ??? success "Solution"
       The regularity (Hölder exponent) of a scaling function $\phi$ is determined by the JSR of the matrices governing its subdivision scheme. A smaller JSR implies a faster decay of coefficients and a smoother limit function.

10. **[Spectral Gap] Relate the JSR to information spreading in random walks.**
    ??? success "Solution"
        For stochastic matrices, the JSR of the matrices restricted to the mean-zero subspace determines the mixing time. The spectral gap $1 - \rho(\Sigma)$ quantifies how quickly the probability distribution converges to the stationary state.

## Chapter Summary

This chapter examines the aggregate growth rates of matrix families:

1. **Analytical Grounding**: Defined the JSR via Rota-Strang limits and established its norm-independence.
2. **Spectral Equivalence**: Decoded the Berger-Wang theorem, linking norm-based growth to the product eigenvalues.
3. **Stability Power**: Demonstrated the central role of JSR in the absolute stability of switched dynamical systems.
4. **Computational Bounds**: Discussed the undecidability and finiteness properties, highlighting the intrinsic complexity of joint spectral analysis.
