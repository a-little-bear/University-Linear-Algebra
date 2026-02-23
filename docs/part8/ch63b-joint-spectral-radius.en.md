# Chapter 63B: Joint Spectral Radius

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues and Spectral Radius (Ch6) · Matrix Norms (Ch15) · Simultaneous Triangularization (Ch63A)

**Chapter Outline**: Definition (Rota-Strang) → Norm-independence → Generalized Spectral Radius → Berger-Wang Theorem → Finiteness Conjecture → Undecidability → Extremal Norms (Barabanov) → Switched Linear Systems Stability → Wavelet Regularity

**Extension**: The joint spectral radius (JSR) is the central tool for analyzing the stability of switched systems, the convergence of subdivision schemes in wavelets, and the capacity of codes.

</div>

For a set of matrices $\Sigma = \{A_1, \ldots, A_m\}$, the **joint spectral radius** $ho(\Sigma)$ characterizes the maximum growth rate of all possible product sequences $A_{i_1} A_{i_2} \cdots A_{i_k}$ as $k 	o \infty$. This concept generalizes the spectral radius of a single matrix to matrix families, quantifying the "worst-case" instability.

---

## 63B.1 Definitions and the Berger-Wang Theorem

!!! definition "Definition 63B.1 (Joint Spectral Radius)"
    The joint spectral radius of a finite set of matrices $\Sigma$ is defined as:
    $$ho(\Sigma) = \lim_{k 	o \infty} \max_{A \in \Sigma^k} \|A\|^{1/k}$$
    where $\Sigma^k$ is the set of all products of length $k$. This limit exists and is independent of the choice of norm $\|\cdot\|$.

!!! theorem "Theorem 63B.3 (Berger-Wang Theorem, 1992)"
    For finite sets of matrices, the joint spectral radius equals the supremum of the spectral radii of all finite products:
    $$ho(\Sigma) = \sup_{k \ge 1} \max_{A \in \Sigma^k} ho(A)^{1/k}$$

---

## Exercises

1. **[Concept] Contrast the joint spectral radius (JSR) with the standard spectral radius.**
   ??? success "Solution"
       The standard radius $ho(A)$ describes the growth of a single matrix power $A^k$. JSR describes the growth of the *largest* product sequence in a family. JSR is a measure of the aggregate instability of a collection of possible evolution rules.

2. **[Norm Bound] Prove: For any matrix norm $\|\cdot\|$, $ho(\Sigma) \le \max_{A \in \Sigma} \|A\|$.**
   ??? success "Solution"
       By definition, $\|A_{i_1} \dots A_{i_k}\| \le \prod \|A_{i_j}\| \le (\max \|A\|)^k$. Taking the $k$-th root and the limit as $k 	o \infty$ yields the result.

3. **[Calculation] Compute the JSR of $\Sigma = \left\{ \begin{pmatrix} 1 & 1 \ 0 & 1 \end{pmatrix}, \begin{pmatrix} 1 & 0 \ 1 & 1 \end{pmatrix} ight\}$.**
   ??? success "Solution"
       While both matrices have $ho(A_i) = 1$, their product $A_1 A_2 = \begin{pmatrix} 2 & 1 \ 1 & 1 \end{pmatrix}$ has $ho(A_1 A_2) = \frac{3+\sqrt{5}}{2} \approx 2.618$. By Berger-Wang, $ho(\Sigma) \ge ho(A_1 A_2)^{1/2} = \phi \approx 1.618$. One can prove that for this set, $ho(\Sigma) = \phi$ (the golden ratio).

4. **[Berger-Wang Significance] Explain the importance of the Berger-Wang Theorem.**
   ??? success "Solution"
       It proves that the external estimate (based on norms) and the internal estimate (based on eigenvalues of products) converge to the same value. This allows for approximating the JSR by inspecting eigenvalues of finite-length products.

5. **[Stability] Why is the stability of a switched linear system $x_{k+1} = A_{\sigma(k)} x_k$ equivalent to $ho(\Sigma) < 1$?**
   ??? success "Solution"
       If $ho(\Sigma) \ge 1$, there exists a switching sequence $\sigma$ such that the state norm does not decay to zero (and may explode). $ho(\Sigma) < 1$ ensures that for *any* switching rule, the energy of the system strictly decreases in the limit.

6. **[CQLF] What is the implication of having a Common Quadratic Lyapunov Function (CQLF) for a matrix set?**
   ??? success "Solution"
       A CQLF implies $ho(\Sigma) < 1$. It defines an ellipsoidal energy surface that all matrices in the set "contract," providing a geometric certificate of absolute stability.

7. **[Undecidability] Explain why calculating the JSR is considered "undecidable."**
   ??? success "Solution"
       Determining whether $ho(\Sigma) \le 1$ for a set of integer matrices is equivalent to the Halting Problem. This means no general algorithm can produce an exact "yes/no" answer for every input in finite time.

8. **[Finiteness] Define the "Finiteness Conjecture" and its current status.**
   ??? success "Solution"
       The conjecture stated that the JSR is always achieved by the spectral radius of some finite-length product. While true for many practical cases (e.g., nonnegative matrices), it has been proven false in general by counterexamples in $M_2(\mathbb{R})$.

9. **[Wavelets] How does the JSR relate to wavelet regularity?**
   ??? success "Solution"
       The Hölder regularity of a scaling function $\phi$ is determined by the JSR of the transition matrices associated with the subdivision scheme. A smaller JSR corresponds to a smoother (more differentiable) wavelet.

10. **[Spectral Gap] Describe the relationship between the JSR and the "Spectral Gap" in Markov chains.**
    ??? success "Solution"
        In the context of stochastic matrices, the JSR of the matrices restricted to the zero-sum subspace relates to the rate of convergence to equilibrium (mixing time). The gap $1 - ho(\Sigma)$ determines the speed of "information spreading."

## Chapter Summary

This chapter examines the aggregate growth rates of matrix families:

1. **Analytical Grounding**: Defined the JSR via Rota-Strang limits and established its norm-independence.
2. **Spectral Equivalence**: Decoded the Berger-Wang theorem, linking norm-based growth to product eigenvalues.
3. **Stability Applications**: Demonstrated the central role of JSR in the absolute stability of switched dynamical systems.
4. **Computational Limits**: Addressed the undecidability and the finiteness properties, highlighting the complexity inherent in joint spectral analysis.
