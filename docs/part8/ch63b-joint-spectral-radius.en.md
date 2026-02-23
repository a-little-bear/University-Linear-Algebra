# Chapter 63B: Joint Spectral Radius (JSR)

<div class="context-flow" markdown>

**Prerequisites**: Spectral Radius (Ch14) · Matrix Norms (Ch15) · Matrix Analysis (Ch14) · Switched Systems Basics

**Chapter Outline**: From Single Matrices to Matrix Sets → Definition of the Joint Spectral Radius (JSR) → Generalization of Gelfand’s Formula → Relationship with Individual Spectral Radii → Computational Challenges: NP-hardness and Undecidability → Extremal Norms and Barabanov Norms → The Rise and Fall of the Finiteness Conjecture → Numerical Approximation (Polyhedral Norms, SOS Programming) → Applications: Stability of Switched Systems, Fractal Geometry, and Regularity of Wavelets

**Extension**: The Joint Spectral Radius characterizes the long-term growth rate of a system under "worst-case" switching; it reveals that even if every matrix in a set is individually stable, their combinations can still lead to divergence, marking the "red line" for complex dynamical systems.

</div>

When we study a time-varying linear system $\mathbf{x}_{k+1} = A_k \mathbf{x}_k$, where $A_k$ is chosen arbitrarily or randomly from a set $\mathcal{A}$, the traditional concept of spectral radius is no longer sufficient. The **Joint Spectral Radius** (JSR) was introduced to characterize the asymptotic growth rate of such a family of matrices under the worst-possible switching sequence. It is not only central to stability analysis but also bridges algebra, combinatorics, and measure theory.

---

## 63B.1 Definition and Core Formulas

!!! definition "Definition 63B.1 (Joint Spectral Radius)"
    Let $\mathcal{A}$ be a finite set of square matrices. Its **Joint Spectral Radius** $\rho(\mathcal{A})$ is defined as:
    $$\rho(\mathcal{A}) = \lim_{k \to \infty} \max_{A_{i_j} \in \mathcal{A}} \|A_{i_k} \cdots A_{i_1}\|^{1/k}$$
    This limit exists and is independent of the choice of norm $\|\cdot\|$.

!!! theorem "Theorem 63B.1 (Generalized Gelfand Formula)"
    For any matrix norm, we have:
    $$\rho(\mathcal{A}) = \inf_{k \ge 1} \max_{A \in \mathcal{A}^k} \|A\|^{1/k}$$
    This establishes JSR as the intrinsic growth rate of the norm of infinite products.

---

## 63B.2 Stability and Counter-Intuitive Phenomena

!!! note "Stability Mismatch"
    It is possible for $\rho(\mathcal{A}) > 1$ even if every matrix $A \in \mathcal{A}$ satisfies $\rho(A) < 1$.
    **Reason**: The invariant subspaces of different matrices may not align. This leads to transient growth during switching which, over infinite sequences, can accumulate into a global explosion of the state.

---

## 63B.3 Computational Complexity

!!! challenge "NP-hardness and Undecidability"
    Computing $\rho(\mathcal{A})$ is proven to be **NP-hard**. Furthermore, the problem of determining whether $\rho(\mathcal{A}) \le 1$ is undecidable. This implies there is no general polynomial-time algorithm for exact JSR computation. Current research relies on numerical approximations using convex optimization.

---

## 63B.4 Applications: Wavelets and Fractals

!!! technique "Wavelet Regularity"
    In the construction of compactly supported wavelets, the regularity (smoothness) of the wavelet function is directly determined by the joint spectral radius of a set of refinement matrices. A smaller JSR implies a smoother wavelet.

---

## Exercises


****
??? success "Solution"
     It equals the standard spectral radius $\rho(A)$.


****
??? success "Solution"
     Consider a sequence that repeatedly uses the same matrix $A$. Its growth rate is $\rho(A)$. Since the JSR is the supremum over all possible sequences, it must be at least as large as the maximum individual growth rate.


****
??? success "Solution"
     Since $A_1 A_2 = 0$ and $A_2 A_1 = 0$, any product of length $\ge 2$ is the zero matrix. The maximum norm is provided by the individual matrices, so $\rho = 0.8$.


****
??? success "Solution"
     $A_1 = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, A_2 = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}$. $\rho(A_1)=\rho(A_2)=0$. However, $A_1 A_2 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$, which has spectral radius 1. Thus $\rho(\{A_1, A_2\}) = 1^{1/2} = 1$.


****
??? success "Solution"
     $\rho(\mathcal{A}) \le \gamma$. This follows from the sub-multiplicativity of norms: $\|A_k \cdots A_1\| \le \prod \|A_i\| \le \gamma^k$.


****
??? success "Solution"
     It conjectured that for any finite set of matrices, the JSR is always achieved by the spectral radius of some finite product. This has been disproven by counter-examples.


****
??? success "Solution"
     By searching for a common homogeneous polynomial Lyapunov function (which corresponds to a higher-order norm), the estimation of the upper bound of JSR is transformed into a semidefinite program.


****
??? success "Solution"
     It is the necessary and sufficient condition for the switched linear system to be asymptotically stable for **arbitrary** switching sequences.


****
??? success "Solution"
     In the commutative case, $\rho(\mathcal{A}) = \max \rho(A_i)$. Commutativity eliminates the transient growth caused by basis mismatch.

****
??? success "Solution"
    ## Chapter Summary

The Joint Spectral Radius characterizes extreme instability under multi-operator interaction:


****: Proved that local (individual) stability is insufficient for global balance, revealing the geometric cost of non-commutativity in dynamical evolution.

****: The difficulty of computing JSR reflects the essential leap from single to mixed evolution, representing a core barrier in computational algebra for infinite products.

****: Through extremal norms and numerical approximations, JSR theory provides a metric for risk in modern control, signal processing, and fractal research, establishing the algebraic boundaries for the robust operation of complex systems.
